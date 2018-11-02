#!/usr/bin/env python
# -* coding: utf-8 -*-

"""
Initialise international stations with amp 1 and phase 0 

Created on Tue Aug 28 2018

@author: Alexander Drabent (parts by Maaijke Mevius) (modified by Leah Morabito)
"""

import argparse
import numpy as np
import os
from losoto.h5parm import h5parm
from losoto.lib_operations import *
import pyrap.tables as pt
import logging

def makesolset(MS, data, solset_name):
    solset = data.makeSolset(solset_name)

    antennaFile = MS + "/ANTENNA"
    logging.info('Collecting information from the ANTENNA table.')
    antennaTable = pt.table(antennaFile, ack=False)
    antennaNames = antennaTable.getcol('NAME')
    antennaPositions = antennaTable.getcol('POSITION')
    antennaTable.close()
    antennaTable = solset.obj._f_get_child('antenna')
    antennaTable.append(zip(*(antennaNames,antennaPositions)))

    fieldFile = MS + "/FIELD"
    logging.info('Collecting information from the FIELD table.')
    fieldTable = pt.table(fieldFile, ack=False)
    phaseDir = fieldTable.getcol('PHASE_DIR')
    pointing = phaseDir[0, 0, :]
    fieldTable.close()

    sourceTable = solset.obj._f_get_child('source')
    # add the field centre, that is also the direction for Gain and Common*
    sourceTable.append([('pointing',pointing)])

    return antennaNames

def main(h5parmfile, MSfiles, solset_in='target', solset_out='targetIS'):

    mslist = MSfiles.lstrip('[').rstrip(']').replace(' ','').replace("'","").split(',')

    if len(mslist) == 0:
        logging.error("Did not find any existing directory in input MS list!")
        return(1)
        pass
    else:
        MSfile = mslist[0]
        pass

    if not os.path.exists(h5parmfile):
        logging.error("H5parm file %s doesn't exist!" % h5parmfile)
        return(1)

    if not os.path.exists(MSfile):
        logging.error("MS file %s doesn't exist!" % MSfile)
        return(1)

    if solset_in == solset_out:
        logging.error("Output solset has to be different from input solset!")
        return(1)

    # Open up the h5parm, get an example value
    data = h5parm(h5parmfile, readonly = False)

    # Create a new solset for the data
    if not solset_out in data.getSolsetNames():
        new_station_names = makesolset(MSfile, data, solset_out)
        pass

    # loading solset
    solset = data.getSolset(solset_in)
    OutSolset = data.getSolset(solset_out)
    station_names = solset.getAnt().keys()

    # copy existing information for the tables
    in_soltabNames = solset.getSoltabNames()
    for soltab_name in in_soltabNames:
	logging.info('Processing solution table %s'%(soltab_name) )
        soltab = solset.getSoltab(soltab_name)
        soltab_axes = soltab.getAxesNames()
        soltab_type = ''.join([i for i in soltab_name if not i.isdigit()])
	if soltab_type == 'RMextract':
	    soltab_type = 'rotationmeasure'
        if len(soltab_axes) == 2:
            out_axes = ['ant','time']
            out_axes_vals = [new_station_names, soltab.time]
        else:
            out_axes = ['ant','time','freq','pol']
            out_axes_vals = [new_station_names, soltab.time, soltab.freq, soltab.pol]
        ## check number of antennas
        if len(new_station_names) == soltab.getAxisLen('ant'):
            ## all stations are present in the solution, just copy
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=soltab_axes, weight=True):
                vals = reorderAxes( vals, soltab_axes, out_axes )
                weights = reorderAxes( weights, soltab_axes, out_axes )
                pass
            new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=soltab_name,
                                        axesNames=out_axes, axesVals=out_axes_vals, vals=vals, weights=weights)
        else:
            ## there are new stations in the measurement set, default values have to be added
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=soltab_axes, weight=True):
                vals = reorderAxes( vals, soltab_axes, out_axes )
                weights = reorderAxes( weights, soltab_axes, out_axes )
                pass
            dimension = np.shape(vals)
            old_stations = soltab.getAxisValues('ant').tolist()

            if len(soltab_axes) == 2:
                ## things with a single value per antenna, time
                new_vals = np.ndarray( shape=(len(new_station_names), dimension[1]) )
                new_weights = np.ndarray( shape=(len(new_station_names), dimension[1]) )
                for i, new_station in enumerate(new_station_names):
                    if new_station in old_stations:
                        ## antenna exists, copy solutions
                        ant_index = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(new_station)
                        new_vals[i,:] = vals[ant_axis,:]
                        new_weights[i,:] = weights[ant_axis,:]
                    else:
                        ## antenna does not exist, create default solutions
                        if soltab_type == 'amplitude':
                            new_vals[i,:] = np.ones(shape=(len(new_station_names), dimension[1]) )
                            new_weights[i,:] = np.ones(shape=(len(new_station_names), dimension[1]) )
                        else:
                            new_vals[i,:] = np.zeros(shape=(len(new_station_names), dimension[1]) )
                            new_weights[i,:] = np.ones(shape=(len(new_station_names), dimension[1]) )
                new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=soltab_name,
                                            axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)

            else:
                ## things with polarisation and frequency information
                new_vals = np.ndarray(shape=(len(new_station_names), dimension[1], dimension[2], dimension[3]) )
                new_weights = np.ndarray(shape=(len(new_station_names), dimension[1], dimension[2], dimension[3]) )
                for i, new_station in enumerate(new_station_names):
                    if new_station in old_stations:
                        ## antenna exists, copy solutions
                        ant_index = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(new_station)
			stuff = vals[ant_index,:,:,:]
                        new_vals[i,:,:,:] = vals[ant_index,:,:,:]
                        new_weights[i,:,:,:] = weights[ant_index,:,:,:]
                    else:
                        ## antenna does not exist, create default solutions
                        if soltab_type == 'amplitude':
                            new_vals[i,:,:,:] = np.ones(shape=(dimension[1], dimension[2], dimension[3]) )
                            new_weights[i,:,:,:] = np.ones(shape=(dimension[1], dimension[2], dimension[3]))
                        else:
                            new_vals[i,:,:,:] = np.zeros(shape=(dimension[1], dimension[2], dimension[3]) )
                            new_weights[i,:,:,:] = np.ones(shape=(dimension[1], dimension[2], dimension[3]) )
                new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=soltab_name,
                                            axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)


        soltab = 0
        pass
    data.close()
    return(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Adds phases and amplitudes to international stations if they do not exist in an h5parm.')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which this action should be performed .')
    parser.add_argument('MSfiles', type=str,
                        help='MS for which the new solset shall be created.')
    parser.add_argument('--solset_in', type=str, default='target',
                        help='Input solution set')
    parser.add_argument('--solset_out', type=str, default='targetIS',
                        help='Output solution set (has to be different from input solution set)')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    MSfiles = args.MSfiles
    h5parmfile = args.h5parm
    main(h5parmfile, MSfiles, solset_in=args.solset_in, solset_out=args.solset_out)

