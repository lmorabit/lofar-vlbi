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

def main(h5parmfile, MSfiles, cal_solset=None, solset_in='target', solset_out='targetIS', do_int_stations=False ):

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

        # loading solset
        solset = data.getSolset(solset_in)
        OutSolset = data.getSolset(solset_out)
        station_names = solset.getAnt().keys()
        in_soltabNames = solset.getSoltabNames()
	for testname in in_soltabNames:
	    if 'phase' in testname:
		phase_soltab = testname
        # get new time axis
        if do_int_stations:
	    if 'phase_soltab' in locals():
                tmp = solset.getSoltab(phase_soltab)
	    else:
		logging.error('Phase solutions do not exist in target solset')
		return(1)
        else:
            tmp = solset.getSoltab('RMextract')
        new_times = tmp.time
        tmp = 0
	
	if cal_solset is not None
		calsols = data.getSolset(cal_solset)
		cal_soltabNames = ['polalign', 'clock', 'bandpass']

		logging.info('Finding median values for calibrator and interpolating to target times')
		for cal_soltab_name in cal_soltabNames:
		    logging.info('Interpolating solution table %s'%(cal_soltab_name))
		    soltab = calsols.getSoltab(cal_soltab_name)
		    soltab_axes = soltab.getAxesNames()
		    soltab_type = ''.join([i for i in cal_soltab_name if not i.isdigit()])
		    if soltab_type == 'RMextract':
			soltab_type = 'rotationmeasure'
		    if soltab_type == 'bandpass':
			soltab_type = 'amplitude'
		    if soltab_type == 'polalign':
			soltab_type = 'phase'
		    # format output axes
		    if len(soltab_axes) == 2:
			out_axes = ['time','ant']
			out_axes_vals = [new_times, soltab.ant]
			out_lens = (len(new_times),len(soltab.ant))
		    elif len(soltab_axes) == 4:
			out_axes = ['time','ant','freq','pol']
			out_axes_vals = [new_times, soltab.ant, soltab.freq, soltab.pol]
			out_lens = (len(new_times),len(soltab.ant),len(soltab.freq),len(soltab.pol))
		    elif len(soltab_axes) == 5:
			out_axes = ['time','freq','ant','dir','pol']
			out_axes_vals = [new_times, soltab.freq, soltab.ant, soltab.dir, soltab.pol]
			out_lens = (len(new_times),len(soltab.freq),len(soltab.ant),len(soltab.dir),len(soltab.pol))

		    # get the values, weights, etc.
		    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=soltab_axes, weight=True ):
			vals = reorderAxes( vals, soltab_axes, out_axes )
			weights = reorderAxes( weights, soltab_axes, out_axes )
			pass
		    # find the time axis
		    time_axis_len = soltab.getAxisLen('time')
		    # get the median values for the antennas
		    med_vals = np.nanmedian(vals, axis=0)
		    # replace vals, weights with the medians
		    new_vals = np.ndarray(shape=out_lens)
		    new_weights = np.ones(shape=out_lens)
		    # make an index array
		    for x in np.arange(len(new_times)):
			if len(soltab_axes) == 2:
			    new_vals[x,:] = med_vals
			elif len(soltab_axes) == 4:
			    new_vals[x,:,:,:] = med_vals
			elif len(soltab_axes) == 5:
			    new_vals[x,:,:,:,:] = med_vals
		    new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=cal_soltab_name,
					    axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)


		    soltab = 0
		    pass

        if do_int_stations:
            # copy existing information for the tables
            for soltab_name in in_soltabNames:
                logging.info('Processing solution table %s'%(soltab_name) )
                soltab = solset.getSoltab(soltab_name)
                soltab_axes = soltab.getAxesNames()
                soltab_type = ''.join([i for i in soltab_name if not i.isdigit()])
                if soltab_type == 'RMextract':
                    soltab_type = 'rotationmeasure'
                if soltab_type == 'bandpass':
                    soltab_type = 'amplitude'
                if soltab_type == 'polalign':
                    soltab_type = 'phase'
	        if 'phase' in soltab_type:
		    soltab_type = 'phase'
                if len(soltab_axes) == 2:
                    out_axes = ['time','ant']
                    out_axes_vals = [soltab.time, new_station_names]
                elif len(soltab_axes) == 4:
                    out_axes = ['time','ant','freq','pol']
                    out_axes_vals = [soltab.time, new_station_names, soltab.freq, soltab.pol]
		elif len(soltab_axes) == 5:
		    out_axes = ['time','freq','ant','dir','pol']
		    out_axes_vals = [soltab.time, soltab.freq, new_station_names, soltab.dir, soltab.pol]
	
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
                        new_vals = np.ndarray( shape=(dimension[0],len(new_station_names)) )
                        new_weights = np.ndarray( shape=new_vals.shape )
                        for i, new_station in enumerate(new_station_names):
                            if new_station in old_stations:
                                ## antenna exists, copy solutions
                                ant_index = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(new_station)
                                new_vals[:,i] = vals[:,ant_index]
                                new_weights[:,i] = weights[:,ant_index]
                            else:
                                ## antenna does not exist, create default solutions
                                if soltab_type == 'amplitude':
                                    new_vals[:,i] = np.ones(shape=(dimension[0],len(new_station_names)) )
                                    new_weights[:,i] = np.ones(shape=(dimension[0],len(new_station_names)) )
                                else:
                                    new_vals[:,i] = np.zeros(dimension[0])
                                    new_weights[:,i] = np.ones(dimension[0])
                        new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=soltab_name,
                                                axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)

                    elif len(soltab_axes) == 4:
			len_time = dimension[0]
			len_freq = dimension[2]
			len_pol = dimension[3]

                        ## things with polarisation and frequency information
                        new_vals = np.ndarray(shape=(len_time,len(new_station_names),len_freq,len_pol))
                        new_weights = np.ndarray(shape=new_vals.shape )
                        for i, new_station in enumerate(new_station_names):
                            if new_station in old_stations:
                                ## antenna exists, copy solutions
                                ant_index = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(new_station)
                                new_vals[:,i,:,:] = vals[:,ant_index,:,:]
                                new_weights[:,i,:,:] = weights[:,ant_index,:,:]
                            else:
                                ## antenna does not exist, create default solutions
                                if soltab_type == 'amplitude':
                                    new_vals[:,i,:,:] = np.ones(shape=(len_time,len_freq,len_pol) )
                                    new_weights[:,i,:,:] = np.ones(shape=(len_time,len_freq,len_pol))
                                else:
                                    new_vals[:,i,:,:] = np.zeros(shape=(len_time,len_freq,len_pol))
                                    new_weights[:,i,:,:] = np.ones(shape=(len_time,len_freq,len_pol) )
                        new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=soltab_name,
                                                axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)
		    elif len(soltab_axes) == 5:
			len_time = dimension[0]
			len_freq = dimension[1]
			len_dir = dimension[3]
			len_pol = dimension[4]

                        ## things with polarisation and frequency information
                        new_vals = np.ndarray(shape=(len_time,len_freq,len(new_station_names),len_dir,len_pol))
                        new_weights = np.ndarray(shape=new_vals.shape )
                        for i, new_station in enumerate(new_station_names):
                            if new_station in old_stations:
                                ## antenna exists, copy solutions
                                ant_index = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(new_station)
                                new_vals[:,:,i,:,:] = vals[:,:,ant_index,:,:]
                                new_weights[:,:,i,:,:] = weights[:,:,ant_index,:,:]
                            else:
                                ## antenna does not exist, create default solutions
                                if soltab_type == 'amplitude':
                                    new_vals[:,:,i,:,:] = np.ones(shape=(len_time,len_freq,len_dir,len_pol) )
                                    new_weights[:,:,i,:,:] = np.ones(shape=(len_time,len_freq,len_dir,len_pol))
                                else:
                                    new_vals[:,:,i,:,:] = np.zeros(shape=(len_time,len_freq,len_dir,len_pol))
                                    new_weights[:,:,i,:,:] = np.ones(shape=(len_time,len_freq,len_dir,len_pol) )
                        new_soltab = OutSolset.makeSoltab(soltype=soltab_type, soltabName=soltab_name,
                                                axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)


                soltab = 0
                pass
            data.close()
            return(0)
    else:
        data.close()
        logging.info('Solutions %s already exist, not adding to h5parm.'%solset_out)
        return(0)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Adds phases and amplitudes to international stations if they do not exist in an h5parm.')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which this action should be performed .')
    parser.add_argument('MSfiles', type=str,
                        help='MS for which the new solset shall be created.')
    parser.add_argument('--cal_solset', type=str, default='calibrator',
			help='Input calibrator solutions')
    parser.add_argument('--solset_in', type=str, default='target',
                        help='Input solution set')
    parser.add_argument('--solset_out', type=str, default='targetIS',
                        help='Output solution set (has to be different from input solution set)')
    parser.add_argument('--do_int_stations', action='store_true', dest='do_int_stations' )

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    MSfiles = args.MSfiles
    h5parmfile = args.h5parm
    main(h5parmfile, MSfiles, cal_solset=args.cal_solset, solset_in=args.solset_in, solset_out=args.solset_out, do_int_stations=args.do_int_stations)

