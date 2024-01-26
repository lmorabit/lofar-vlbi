#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Search VLASS for a given RA and Dec
Downloads best subtile
Original code by Anna Ho
Edited by Roland Timmerman
"""


#Imports
import os
import sys
import errno
import subprocess
import glob
import argparse
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import casacore.tables as pt
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

#Settings to potentially tweak
summary_file_location = "VLASS_dyn_summary.php"
crop = True
crop_scale = 256
consider_QA_rejected = True


def sign(x):
    """
    Returns the mathematical sign of a number
    Input: x (int/float)
    Returns: +1 or -1
    """
    
    if x<0:
        return -1
    else:
        return 1

def get_tiles():
    """
    Read tiles from tile catalog. If file missing, try wget https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php
    Input: None
    Returns: tile names (numpy str), dec_min (numpy float), dec_max (numpy float), ra_min (numpy float), ra_max (numpy float), observing epoch (numpy str), observation date (numpy str)
    """

    summary_file_location = "VLASS_dyn_summary.php"
    if not os.path.exists(summary_file_location):
        os.system('wget https://raw.githubusercontent.com/lmorabit/lofar-vlbi/master/VLASS_dyn_summary.php')
    
    names, dec_min, dec_max, ra_min, ra_max, epoch, obsdate = np.loadtxt(summary_file_location, skiprows=3, unpack=True, dtype='str', usecols=(0,1,2,3,4,5,6))

    dec_min = dec_min.astype(float)
    dec_max = dec_max.astype(float)
    ra_min = ra_min.astype(float)
    ra_max = ra_max.astype(float)

    return (names, dec_min, dec_max, ra_min, ra_max, epoch, obsdate)


def search_tiles(tiles, c):
    """
    Search the tile catalog for tiles containing the input coordinate
    Input: tile properties (see get_tiles() output), c (SkyCoord object)
    Returns: tile name (str), observing epoch (str), observing date (str)
    """
    ra_h = c.ra.hour
    dec_d = c.dec.deg
    names, dec_min, dec_max, ra_min, ra_max, epochs, obsdate = tiles
    has_dec = np.logical_and(dec_d >= dec_min, dec_d < dec_max)
    has_ra = np.logical_and(ra_h >= ra_min, ra_h < ra_max)
    in_tile = np.logical_and(has_ra, has_dec)
    name = names[in_tile]
    epoch = epochs[in_tile]
    date = obsdate[in_tile]    
    if len(name) == 0:
        raise IndexError("Zero VLASS tiles available for the given coordinate")
    c_grid = SkyCoord(7.5*(ra_min[in_tile]+ra_max[in_tile]), 0.5*(dec_min[in_tile]+dec_max[in_tile]), unit='deg', frame='icrs')
    dist = c_grid.separation(c)
    best_idx = np.argmin(dist)
    return name[best_idx], epoch[best_idx], date[best_idx]


def get_subtiles(tilename, epoch, consider_QA_rejected):
    """
    For a given tile name, get the subtile filenames in the VLASS directory
    Parse those filenames and return a list of subtile RA and Dec.
    RA and Dec returned as a SkyCoord object
    """
    
    #Obtain the HTML for the given tile
    urlpath = urlopen("https://archive-new.nrao.edu/vlass/quicklook/{}v2/{}".format(epoch, tilename))
    string = urlpath.read().decode('utf-8').split("\n")

    if consider_QA_rejected:
        #Obtain the HTML for the QA Rejected
        urlpath_rejected = urlopen("https://archive-new.nrao.edu/vlass/quicklook/{}v2/QA_REJECTED".format(epoch))
        string += urlpath_rejected.read().decode('utf-8').split("\n")

    #Select only the subtile parts
    vals = np.array([val.strip() for val in string if ("href" in val.strip()) and (tilename in val.strip())])
    
    #Select the coordinate part. You want the 'VLASS1.1.ql.T25t12.J150000+603000.10.2048.v1/' bit
    fname = np.array([val.split("\"")[7] for val in vals])

    #Split out the actual coordinate string
    pos_raw = np.array([val.split(".")[4] for val in fname])
    if '-' in pos_raw[0]:
        # dec < 0
        ra_raw = np.array([val.split("-")[0] for val in pos_raw])
        dec_raw = np.array([val.split("-")[1] for val in pos_raw])
    else:
        # dec > 0
        ra_raw = np.array([val.split("+")[0] for val in pos_raw])
        dec_raw = np.array([val.split("+")[1] for val in pos_raw])
    ra = []
    dec = []
    for ii,val in enumerate(ra_raw):
        if val[1:3] == '24':
            rah = '00'
        else:
            rah = val[1:3]
        ra.append("{}h{}m{}s".format(rah, val[3:5], val[5:]))
        dec.append("{}d{}m{}s".format(dec_raw[ii][:2], dec_raw[ii][2:4], dec_raw[ii][4:]))
    ra = np.array(ra)
    dec = np.array(dec)
    c = SkyCoord(ra, dec, frame='icrs')#.directional_offset_by(45*u.deg, 0.75*u.deg)
    return fname, c


def get_cutout(imname, c, crop_scale):
    #Define output name
    output_fits = "vlass_poststamp.fits"
    
    #Get header info
    hdu_list = fits.open(imname)
    header = hdu_list[0].header
    data = hdu_list[0].data[0,0,:,:]
    
    #Obtain header and drop useless axes
    wcs = WCS(header)
    wcs = wcs.dropaxis(2).dropaxis(2)
    
    pixel_coords = skycoord_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg', frame='icrs'), wcs)
    
    if pixel_coords[0] < 0  or pixel_coords[1] < 0 or pixel_coords[0] > data.shape[0]  or pixel_coords[1] > data.shape[1]:
        subprocess.call("rm -f {}".format(imname), shell=True)
        raise Exception("Requested coordinate not within the available subtiles. Consider running with consider_QA_rejected=True to also search additional subtiles which failed initial QA checks")
        
    #Produce a cutout
    cutout = Cutout2D(data, c, (crop_scale, crop_scale), wcs=wcs)
        
    #Update the HDU
    hdu_list[0].data = cutout.data
    new_header = cutout.wcs.to_header()
    hdu_list[0].header.update(new_header)
    hdu_list[0].header.set('NAXIS', 4)
    hdu_list[0].header.insert('NAXIS2', ('NAXIS3', 1), after=True)
    hdu_list[0].header.insert('NAXIS3', ('NAXIS4', 1), after=True)
    hdu_list[0].header.remove('WCSAXES', ignore_missing=True)
    hdu_list[0].header.remove('MJDREF', ignore_missing=True)
    hdu_list[0].header.remove('MJD-OBS', ignore_missing=True)
    
    #Write the new fits
    hdu_list.writeto(output_fits, overwrite=True)
    
    #Cleanup
    subprocess.call("rm -f {}".format(imname), shell=True)

    return output_fits


def search_vlass(c, crop=False, crop_scale=256, consider_QA_rejected=False):
    """ 
    Searches the VLASS catalog for a source

    Parameters
    ----------
    c: coordinates as SkyCoord object
    """
    # Find the VLASS tile
    tiles = get_tiles()
    tilename, epoch, obsdate = search_tiles(tiles, c)

    subtiles, c_tiles = get_subtiles(tilename, epoch, consider_QA_rejected)
    dist = c.separation(c_tiles)
    subtile = subtiles[np.argmin(dist)]

    imname = "{}.I.iter1.image.pbcor.tt0.subim.fits".format(subtile[:-1])
    if len(glob.glob(imname)) == 0:
        url_get = "https://archive-new.nrao.edu/vlass/quicklook/{}v2/{}/{}".format(epoch, tilename, subtile)
        fname = "{}{}".format(url_get, imname)
        subprocess.call("wget {}".format(fname), shell=True)
        if not os.path.exists(subtile) and consider_QA_rejected:
            url_get = "https://archive-new.nrao.edu/vlass/quicklook/{}v2/QA_REJECTED/{}".format(epoch, subtile)
            fname = "{}{}".format(url_get, imname)
            subprocess.call("wget {}".format(fname), shell=True)
    if crop:    
        out = get_cutout(imname, c, crop_scale=crop_scale)
        return out
    else:
        return imname


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        '''
        Searches VLASS for a source.
        User needs to supply RA (in decimal degrees), Dec (in decimal degrees)
        
        Usage: vlass_search.py <RA [deg]> <Dec [deg]> or vlass_search.py <MS>
        ''', formatter_class=argparse.RawTextHelpFormatter)
    
    #Test coords: 232.517541667 60.4054444
    try:
        ra = 15*(float(sys.argv[1]) + float(sys.argv[2])/60 + float(sys.argv[3])/3600)
        dec = float(sys.argv[4]) + sign(float(sys.argv[4]))*float(sys.argv[5])/60 + sign(float(sys.argv[4]))*float(sys.argv[6])/3600
    except IndexError:
        ra = float(sys.argv[1])
        dec = float(sys.argv[2])
    except ValueError:
        pi = 3.14159265358979
        table = pt.table(sys.argv[1]+"::FIELD")
        direction = table.getcol('PHASE_DIR').squeeze()
        ra=(direction[0]%(2*pi))/pi*180
        dec=direction[1]/pi*180
        table.close()
    except Exception:
        raise TypeError("Incorrect inputs given. Usage: vlass_search.py <RA [deg]> <Dec [deg]> or vlass_search.py <MS>")
        
    c = SkyCoord(ra, dec, unit='deg')

    if not glob.glob(summary_file_location):
        os.system('wget https://raw.githubusercontent.com/lmorabit/lofar-vlbi/master/VLASS_dyn_summary.php')
        #raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), summary_file_location)
    
    search_vlass(c, crop=crop, crop_scale=crop_scale, consider_QA_rejected=consider_QA_rejected) 
