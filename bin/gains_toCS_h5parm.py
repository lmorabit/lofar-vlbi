#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copy H5parm values from a reference station to all new stations, e.g., from ST001 to all core stations

Created on Tue Aug 28 2018

@author: Alexander Drabent (parts by Maaijke Mevius)
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

def main(h5parmfile, MSfiles, solset_in = 'sol000', solset_out = 'sol001', soltab_list = ['clock000', 'tec000', 'amplitude000'], superstation = 'ST001', restrictToCS = True):

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

    if superstation not in station_names:
        logging.error("Couldn't find station " + str(superstation))
        return(1)
    
    # start copying
    for soltab_name in soltab_list:
        logging.info("Running copySTgains_toCS on: " + soltab_name)
        soltab = solset.getSoltab(soltab_name)
        STindex = soltab.getAxisValues('ant', ignoreSelection = True).tolist().index(superstation)
        
        if 'clock' in soltab_name or 'tec' in soltab_name:
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes = ['time', 'ant'], weight = True):
                vals = reorderAxes( vals, soltab.getAxesNames(), ['time', 'ant'] )
                weights = reorderAxes( weights, soltab.getAxesNames(), ['time', 'ant'] )
                pass
            pass
        elif 'amplitude' in soltab_name:
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['pol','ant','freq', 'time'], weight=True):
                vals = reorderAxes( vals, soltab.getAxesNames(), ['time', 'ant', 'freq', 'pol'] )
                weights = reorderAxes( weights, soltab.getAxesNames(), ['time', 'ant', 'freq', 'pol'] )
                pass
            pass
        else:
            loggin.error('No clock, tec or amplitude soltab has been found or specified.')
            return(1)
            pass
        
        dimension = np.shape(vals)
        
        if 'clock' in soltab_name or 'tec' in soltab_name:
            new_vals    = np.ndarray(shape = (dimension[0], len(new_station_names)))
            new_weights = np.ndarray(shape = (dimension[0], len(new_station_names)))
            pass
        if 'amplitude' in soltab_name:
            new_vals    = np.ndarray(shape = (dimension[0], len(new_station_names), dimension[2], dimension[3]))
            new_weights = np.ndarray(shape = (dimension[0], len(new_station_names), dimension[2], dimension[3]))
            pass
        
        for i, new_station in enumerate(new_station_names):
            if new_station in station_names:
                ant_index        = soltab.getAxisValues('ant', ignoreSelection = True).tolist().index(new_station)
                if 'clock' in soltab_name or 'tec' in soltab_name:
                    new_vals[:,i]    = vals[:,ant_index]
                    new_weights[:,i] = weights[:,ant_index]
                    pass
                if 'amplitude' in soltab_name:
                    new_vals[:,i,:,:]    = vals[:,ant_index,:,:]
                    new_weights[:,i,:,:] = weights[:,ant_index,:,:]
                    pass
                pass
            else:
                if restrictToCS and 'CS' not in new_station:
                    logging.info('RestrictToCS: Omitting station ' + new_station)
                    continue
                    pass
                logging.info( 'Adding ' + str(soltab_name) + ' to ' + new_station)
                if 'clock' in soltab_name or 'tec' in soltab_name:
                    new_vals[:,i]    = vals[:,STindex]
                    new_weights[:,i] = weights[:,STindex]
                    pass
                if 'amplitude' in soltab_name:
                    new_vals[:,i,:,:]    = vals[:,STindex,:,:]
                    new_weights[:,i,:,:] = weights[:,STindex,:,:]
                    pass                 
                pass
            pass
        
        if 'clock' in soltab_name:
            new_soltab = OutSolset.makeSoltab(soltype='clock', soltabName=soltab_name,
                                        axesNames=['time', 'ant'],
                                        axesVals=[soltab.time, new_station_names],
                                        vals=new_vals, weights=new_weights)
            pass
        elif 'tec' in soltab_name:
            new_soltab = OutSolset.makeSoltab(soltype='tec', soltabName=soltab_name,
                                        axesNames=['time', 'ant'],
                                        axesVals=[soltab.time, new_station_names],
                                        vals=new_vals, weights=new_weights)
            pass
        elif 'amplitude' in soltab_name:
            new_soltab = OutSolset.makeSoltab(soltype='amplitude', soltabName=soltab_name,
                                        axesNames=['time', 'ant', 'freq', 'pol'],
                                        axesVals=[soltab.time, new_station_names, soltab.freq, soltab.pol],
                                        vals=new_vals, weights=new_weights)
            pass
        
        soltab = 0
        pass

    return(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Copies calibration values from a reference station to a set of new stations, e.g., from the superterp (ST001) to all core stations.')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which this action should be performed .')
    parser.add_argument('MSfiles', type=str,
                        help='MS for which the new solset shall be created.')
    parser.add_argument('--solset_in', type=str, default='sol000',
                        help='Input solution set')
    parser.add_argument('--solset_out', type=str, default='sol001',
                        help='Output solution set (has to be different from input solution set)')
    parser.add_argument('--soltab_list', '--soltab_list', type=str, default='clock000,tec000,amplitude000',
                        help='Comma-separated list of soltabs to be copied')
    parser.add_argument('--superstation', type=str, default='ST001',
                        help='Reference station from which data should be copied')
    parser.add_argument('--restrictToCS',
                        help='Restrict the copy action to core stations only',
                        action='store_true',dest="restrictToCS")

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    MSfiles = args.MSfiles
    h5parmfile = args.h5parm
    soltablist = args.soltab_list.split(',')

    main(h5parmfile, MSfiles, solset_in=args.solset_in, solset_out=args.solset_out, 
                 soltab_list = soltablist, superstation=args.superstation, restrictToCS=args.restrictToCS)
