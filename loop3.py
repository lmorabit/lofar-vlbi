#!/usr/bin/env python
# J. Callingham, L. Morabito, N. Jackson 

# Loop 3 for longbaseline survey pipeline following details here: 
# https://docs.google.com/document/d/1qHICQF1IevEISzKQS4z06TzLQrteuKNW_EvuFyIPBHE/edit

import numpy as np
import os
import sys
import glob
import astropy
import lofar.parmdb
import argparse
import time
import pyrap.tables as pt
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt


print("#----------------------------------------------------------#" )

parser = argparse.ArgumentParser(description="Perfroming Loop 3 of the longbaseline pipeline. This loop first selfcals with a \
    few minute solution interval and checks if all stations are coherent. \
    If any stations are not coherent, it attempts to triple the solution interval to detemine whether to use the station in the selfcal solutions. \
    It then passes the selfcal solution to the facet solution to give a new sub-facet solution.")

parser.add_argument('vis', type=str, help='The filename of the visibility you want to read in.')

args = parser.parse_args()

vis = args.vis

def coherence_metric(phasesol_xx,phasesol_yy):

    diff = phasesol_xx - phasesol_yy
    # unwrap phases: it will correct an array of phases modulo 2pi such that all jumps are less than or equal to pi
    diff = np.unwrap(diff)
    # taking mean of absolute difference 
    coh_met = np.nanmean(np.gradient(abs(diff))**2)
    return(coh_met)

def writePhaseOnlyParset( parsetname="phaseonly.parset", incol="DATA", outcol="CORRECTED_DATA", solint=180, outms="." ):
    # open write-only parset; this will overwrite previous files of the same name
    with open(parsetname,"w") as f:
        f.write("msin.datacolumn=%s\n"%incol)
        f.write("msout=%s\n"%outms)
        f.write("msout.datacolumn=%s\n"%outcol)
        f.write("steps=[gaincal]\n")
        f.write("gaincal.usemodelcolumn=true\n")
        f.write("gaincal.caltype=phaseonly\n")
        f.write("gaincal.solint=%i\n"%(solint))
        # f.write("gaincal.maxiter=200\n")
        f.write("gaincal.usebeammodel=false")
        # f.write("gaincal.tolerance=0.0001\n")
        f.write("gaincal.applysolution=True\n")
    f.close()

def imaging(vis,niters,threshold):
    imsize = 4096
    rms_bkg = '-local-rms'
    save_cc = '-save-source-list' # will save in coloumn MODEL_DATA that will be used in the self-cal loop
    imagename = 'test'
     
    # first image is a shallow clean to get a model
    ss = "wsclean -data-column CORRECTED_DATA -scale 0.05asec -size %i %i %s -maxuv-l 800000 -gain 0.1 -mgain 0.85 -stop-negative -niter %i -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold %i -multiscale -auto-mask 7 %s -reorder -name %s %s"%(imsize,imsize,rms_bkg,niters,threshold,save_cc,imagename,vis)
    os.system( ss )

# First shallow cleaning with solutions and 
# using wsclean to generate a sky model from clean components to then use in selfcal loop

imaging(vis,1000,10)

# Running selfcal with a 3 minute solution interval, then tripling solution interval if no good.

for sol_int_loop in [180.,540.,1620.,4860.]:
    #reading in phases 
    instrument_name = vis + '/' + 'instrument'
    pdb = lofar.parmdb.parmdb(instrument_name)
    parms = pdb.getValuesGrid("*")
    key_names = parms.keys()
    antenna_list = np.copy(key_names)
    solint = sol_int_loop # staring solution interval is 3 mins, assuming a interval of 1 s, as per data taken

    for ii in range(len(key_names)):
        string_a = str(key_names[ii])
        split_a  = string_a.split( ":" )
        antenna_list[ii] = split_a[4]

    antenna_list = np.unique(antenna_list)
    no_antenna = len(antenna_list)

    pol_list = ['0:0','1:1'] # 0:0 is XX, 1:1 is YY

    coherence_metric_store = np.zeros(len(antenna_list))
    ind = 0

    for antenna in antenna_list:

        phasesol_xx = parms['Gain:'+pol_list[0]+':Phase:'+ antenna]['values'][::].flatten() # without flattening these are 2D arrays with sol per subband (avg.)
        phasesol_yy = parms['Gain:'+pol_list[1]+':Phase:'+ antenna]['values'][::].flatten()
        coherence_metric_store[ind] = coherence_metric(phasesol_xx,phasesol_yy)
        ind = ind + 1

    # Testing if all stations are less than 1.6, a value set empircally by N. Jackson based on closure phases
    if len(coherence_metric_store[coherence_metric_store >= 1.6]) == 0:
        print 'Success. All stations are coherent. Adopting solution and adding to facet solution to give new sub-facet solution.'
        break
    elif len(coherence_metric_store[coherence_metric_store >= 1.6]) == 14 and solint > 3600.:
        print 'All stations are incoherent. Solution interval for self-cal is now > 1 hr. Do not use. Exiting...'
        break
    elif len(coherence_metric_store[coherence_metric_store >= 1.6]) == 14 and solint < 3600.:
        print 'All stations are incoherent. Tripling solution interval for self-cal and testing for coherency again.'
        writePhaseOnlyParset(solint = sol_int_loop*3)
        print "Performing phase-only self-calibration again." 
        ss = 'DPPP DPPP_calibrate_phaseonly.parset msin=%s'%(vis)
        os.system( ss )
        print "Phase selfcal done. Testing coherence again."
        