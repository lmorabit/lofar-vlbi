
from __future__ import print_function
from functools import partial
from multiprocessing import Pool
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from losoto.lib_operations import reorderAxes
import losoto.h5parm as lh5  # on CEP3, "module load losoto"
import astropy.units as u
import pyrap.tables as pt
import numpy as np
import argparse
import csv
import datetime
import os
import subprocess
import uuid
import glob
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct



def interpolate_nan(x_):
    '''Interpolate NaN values using this answer from Stack Overflow:
    https://stackoverflow.com/a/6520696/6386612. This works even if the first
    or last value is nan or if there are multiple nans. It raises an error if
    all values are nan.
    Args:
    x_ (list or NumPy array): Values to interpolate.
    Returns:
    The interpolated values. (NumPy array)'''

    x_ = np.array(x_)
    if np.isnan(x_).all():  # if all values are nan
        raise ValueError('All values in the array are nan, so interpolation is'
                         ' not possible.')
    nans, x = np.isnan(x_), lambda z: z.nonzero()[0]
    x_[nans] = np.interp(x(nans), x(~nans), x_[~nans])

    return x_


def coherence_metric_TEC(xx):
    '''Calculates the coherence metric by comparing the XX and YY phases.
    Args:
    xx (list or NumPy array): TEC solutions.
    Returns:
    The coherence metric. (float)'''
    
    try:
        xx = interpolate_nan(xx)
    except:
        return np.nan  # if the values are all nan, they cannot be interpolated
                       # so return a coherence value of nan also

    return np.nanmean(np.gradient(abs(np.unwrap(xx))) ** 2.)

def plugin_main( args, **kwargs ):
    """
        given a list of h5parms, finds the one with the most coherent TEC
        
        Parameters
        ----------
        h5pattern : str
            pattern for glob search for h5parms
        ----------

        Returns
        -------
            name of h5parm with most coherent TEC
        -------
    """
    h5pattern = kwargs['h5pattern']
    jobdir = kwargs['jobdir']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    if type(h5pattern) == str:
	tmp = os.path.join(jobdir,h5pattern)
        h5files = glob.glob(tmp)
    else:
        h5files = h5pattern

    median_vals = {}
    for h5parm in h5files:

        h = lh5.h5parm(h5parm)
        solname = h.getSolsetNames()[0]  # set to -1 to use only the last solset
        solset = h.getSolset(solname)
        soltabnames = solset.getSoltabNames()
        soltab = solset.getSoltab('tec000')
        stations = soltab.ant
        source = solset.getSou()  # dictionary
        direction = np.degrees(np.array(source[list(source.keys())[0]]))  # degrees
        generator = soltab.getValuesIter(returnAxes=['freq', 'time'])
        evaluations, temporary = {}, {} # evaluations holds the coherence metrics
        
        for g in generator:
            temporary[g[1]['ant']] = np.squeeze(g[0])

        for station in stations:
            xx = temporary[station]
            try:
                cohs = []
                num_solns, num_freqs = xx.shape
                for i in range(num_freqs):
                    cohs.append(coherence_metric(xx[:, i], xx[:,i]))
                coh = np.mean(cohs)
                print('{} {} coherence: {} ({} frequency axes)'.format(h5parm, station, coh, num_freqs)) 
                evaluations[station] = coh
            except:
                coh = coherence_metric_TEC(xx)
                print('{} {} coherence: {}'.format(h5parm, station, coh))
                evaluations[station] = coh

        h.close()


        median_vals[h5parm] = np.median(evaluations.values())

    best_h5parm = min( median_vals, key=median_vals.get ) 
    print( 'Using %s as the best TEC solutions.'%best_h5parm)

    map_out = DataMap([])
    map_out.data.append(DataProduct( 'localhost', best_h5parm, False ))

    fileid = os.path.join(mapfile_dir, filename) 
    map_out.save(fileid)


    result = {'h5tecfile': best_h5parm, 'mapfile':fileid}

    return result
