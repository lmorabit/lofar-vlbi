#!/usr/bin/env python
#from scipy import stats  
import os, sys, logging
import numpy as np
#import csv
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from lmfit import SkewedGaussianModel
import pyvo as vo
import pyrap.tables as pt
from astropy.table import Table, Column, vstack, unique
import argparse
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

def grab_coo_MS(MS):
    """
    Read the coordinates of a field from one MS corresponding to the selection given in the parameters
    Parameters
    ----------
    MS : str
        Full name (with path) to one MS of the field
    Returns
    -------
    RA, Dec : float,float
        coordinates of the field (RA, Dec in deg , J2000)
    """

    # reading the coordinates ("position") from the MS
    # NB: they are given in rad,rad (J2000) 
    [[[ra,dec]]] = pt.table(MS+'/FIELD', readonly=True, ack=False).getcol('PHASE_DIR')
    # [[[ra,dec]]] = pt.table(MS, readonly=True, ack=False).getcol('PHASE_DIR')

    # RA is stocked in the MS in [-pi;pi]
    # => shift for the negative angles before the conversion to deg (so that RA in [0;2pi])
    if ra<0:
        ra=ra+2*np.pi

    # convert radians to degrees
    ra_deg =  ra/np.pi*180.
    dec_deg = dec/np.pi*180.

    # and sending the coordinates in deg
    return ra_deg,dec_deg

def input2strlist_nomapfile(invar):
    """ from bin/download_IONEX.py
    give the list of MSs from the list provided as a string
    """

    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            str_list = [invar.strip(' \'\"')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list

def my_lotss_catalogue( ms_input, Radius=1.5, bright_limit_Jy=5. ):

    """
    Download the LoTSS skymodel for the target field
    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    Radius : float (default = 1.5)
        Radius for the LOTSS cone search in degrees
    
    """
    print "DOWNLOADING LOTSS Skymodel for the target field"

    # Reading a MS to find the coordinate (pyrap)
    RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    ## this is the tier 1 database to query
    #url = 'http://vo.astron.nl/lofartier1/q/cone/scs.xml'
    # HETDEX database.
    url = 'https://vo.astron.nl/hetdex/lotss-dr1/cone/scs.xml'

    ## query the database
    query = vo.dal.scs.SCSQuery( url )
    query['RA'] = float( RATar )
    query['DEC'] = float( DECTar )
    query.radius = float( Radius )
    t = query.execute()

    ## convert to VO table
    tb = t.votable.to_table()
    flux_sort = tb.argsort('Total_flux')
    tb_sorted = tb[flux_sort[::-1]]
    ## and keep only some of the columns
    tb_final = tb_sorted['Source_Name', 'RA', 'DEC','Total_flux','Peak_flux', 'Major', 'Minor', 'PA', 'DC_Maj', 'DC_Min', 'DC_PA', 'LGZ_Size', 'LGZ_Width', 'LGZ_PA', 'Isl_rms']
    resolved = np.where(is_resolved(tb_final['Total_flux'], tb_final['Peak_flux'], tb_final['Isl_rms']), 'R', 'U')
    tb_final['Resolved'] = resolved
    tb_final.rename_column('Source_Name', 'Source_id')

    return tb_final

def my_lbcs_catalogue( ms_input, Radius=1.5 ):

    """
    Download the LBCS skymodel for the target field
    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    Radius : float (default = 1.5)
        Radius for the LOTSS cone search in degrees
    
    """
    print "DOWNLOADING LBCS Skymodel for the target field"

    # Reading a MS to find the coordinate (pyrap)
    RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    ## this is the tier 1 database to query
    url = 'http://vo.astron.nl/lbcs/lobos/cone/scs.xml'

    ## query the database
    query = vo.dal.scs.SCSQuery( url )
    query['RA'] = float( RATar )
    query['DEC'] = float( DECTar )
    query.radius = float( Radius )
    t = query.execute()

    ## convert to VO table
    tb = t.votable.to_table()

    #### Workaround to sort and pick good calibrator info from tb array ###########
    counts=[]
    P_count  = 0
    for i in tb:
        b = i[5].count('P')      #### Count of 'P' - good baselines       
        counts.append(b)
        if b >=2:
            P_count = P_count + 1    #### To determine how many sources to take 
    print 'Good sources - ' + str(P_count)
    if P_count == 0:
        logging.critical('There are no good LBCS sources within the given radius. Check your source is within the LBCS footprint and increase the search radius. Exiting...')
        return

    inds = np.argsort(counts)
    tb_sorted =tb[inds[::-1]]
    len_array = len(tb_sorted)
    for i in range ((len_array-P_count)):
        len_array-=1
        tb_sorted.remove_row(len_array)

    ## remove duplicates
    tb_tmp = np.array( tb_sorted['raj2000','decj2000'] )
    result = [ idx for idx, item in enumerate( tb_tmp ) if item in tb_tmp[:idx] ]
    tb_sorted.remove_rows(result)

    ## keep only some columns
    tb_out = tb_sorted['raj2000','decj2000','ObsID']

    return tb_out

def find_close_objs(lo, lb, tolerance=5.):

    ## get the RA and DEC for both catalogues
    lotssRA= np.array(lo['RA'])
    lbcsRA= np.array(lb['raj2000'])
    lotssDEC= np.array(lo['DEC'])
    lbcsDEC= np.array(lb['decj2000'])
    ## and the fluxes from LoTSS
    lotssflux=np.array(lo['Total_flux'])
    lotsspeakflux=np.array(lo['Peak_flux'])
    ## source ID, rms, and resolved information
    src_id=np.array(lo['Source_id'])
    rms=np.array(lo['Isl_rms'])
    resolved=np.array(lo['Resolved'])
    lbcsid=np.array(lb['ObsID'])

    ## astropy.tables.Table is fast at adding columns, less fast at adding rows -- since it has to make a copy each time
    Source_id = []
    LBCS_ID = []
    LBCS_RA = []
    LBCS_DEC = []
    LOTSS_RA = []
    LOTSS_DEC = []
    Delta_Position = []
    Total_flux = []
    Peak_flux = []
    peak_to_total = []
    Isl_rms = []
    Resolved = []

    if len(lbcsRA.shape) > 0:
        ## there's more than one LBCS sources, iterate through the list
        for x, val in enumerate(lbcsRA):
            ## find the difference in the positions
            RAres=abs(lotssRA-lbcsRA[x])
            DECres=abs(lotssDEC-lbcsDEC[x])
            delta_position= (((RAres**2)+ (DECres**2))**(0.5))*3600

            ## find the closest match
            c_idx = np.where( delta_position == np.min( delta_position ) )[0]

            ## check if it's less than tolerance (in arcsec)
            if delta_position[c_idx] < tolerance:
	        Source_id.append(src_id[c_idx][0])
		LBCS_ID.append(lbcsid[x])
		LBCS_RA.append(lbcsRA[x])
		LBCS_DEC.append(lbcsDEC[x])
	        LOTSS_RA.append(lotssRA[c_idx][0])
	        LOTSS_DEC.append(lotssDEC[c_idx][0])
	        Delta_Position.append(np.min(delta_position))
	        Total_flux.append(lotssflux[c_idx][0])
	        Peak_flux.append(lotsspeakflux[c_idx][0])
	        peak_to_total.append((lotsspeakflux[c_idx][0]/lotssflux[c_idx][0]))
		Isl_rms.append(rms[c_idx][0])
	        Resolved.append(resolved[c_idx][0])
            else:
                print "No match found in LoTSS for LBCS calibrator %s"%(lbcsid[x])

        result = Table()
        result['Source_id'] = Source_id
	result['LBCS_ID'] = LBCS_ID
        result['LBCS_RA'] = LBCS_RA
        result['LBCS_DEC'] = LBCS_DEC
        result['LOTSS_RA'] = LOTSS_RA
        result['LOTSS_DEC'] = LOTSS_DEC
        result['Delta_Position'] = Delta_Position
        result['Total_flux'] = Total_flux
        result['Peak_flux'] = Peak_flux
        result['peak_to_total'] = peak_to_total
        result['Isl_rms'] = Isl_rms
        result['Resolved'] = Resolved

        pass

    else:
        ## there's only one LBCS source, no iteration needed
        ## find the difference in the positions
        RAres=abs(lotssRA-lbcsRA)
        DECres=abs(lotssDEC-lbcsDEC)
        delta_position= (((RAres**2)+ (DECres**2))**(0.5))*3600

        ## find the closest match
        c_idx = np.where( delta_position == np.min( delta_position ) )[0]

        ## check if it's less than 5 arcsec away
        if delta_position[c_idx] < tolerance:
            result = Table()
            result['Source_id'] = [src_id[c_idx][0]]
            result['LBCS_ID'] = [lbcsid]
            result['LBCS_RA'] = [lbcsRA]
            result['LBCS_DEC'] = [lbcsDEC]
            result['LOTSS_RA'] = [lotssRA[c_idx][0]]
            result['LOTSS_DEC'] = [lotssDEC[c_idx][0]]
            result['Delta_Position'] = [np.min(delta_position)]
            result['Total_flux'] = [lotssflux[c_idx][0]]
            result['Peak_flux'] = [lotsspeakflux[c_idx][0]]
            result['peak_to_total'] = [(lotsspeakflux[c_idx][0]/lotssflux[c_idx][0])]
            result['Isl_rms'] = [rms[c_idx][0]]
            result['Resolved'] = [resolved[c_idx][0]]

        else:
            print "No match found in LoTSS for LBCS calibrator %s"%(lbcsid)

    ## convert result to numpy array
    return result

def plugin_main( args, **kwargs ):

    mapfile_in = kwargs['mapfile_in']
    lotss_radius = kwargs['lotss_radius']
    lbcs_radius  = kwargs['lbcs_radius']
    bright_limit_Jy = float(kwargs['bright_limit_Jy'])
    lotss_result_file = kwargs['lotss_result_file']
    delay_cals_file = kwargs['delay_cals_file']
    subtract_file = kwargs['subtract_file']
    match_tolerance = float(kwargs['match_tolerance'])
    subtract_limit = float(kwargs['subtract_limit'])
    image_limit_Jy = float(kwargs['image_limit_Jy'])
    doDownload = kwargs['doDownload']

    mslist = DataMap.load(mapfile_in)
    MSname = mslist[0].file
    

    if doDownload:
        lotss_catalogue = my_lotss_catalogue( MSname, Radius=lotss_radius, bright_limit_Jy=bright_limit_Jy )
        lbcs_catalogue = my_lbcs_catalogue( MSname, Radius=lbcs_radius )
        result = find_close_objs( lotss_catalogue, lbcs_catalogue, tolerance=match_tolerance )

        ## Need to write the following catalogues:
        ## 1 - delay calibrators -- from lbcs_catalogue
        result.write( delay_cals_file, format='csv' )
	print 'wrote'+delay_cals_file

        ####### sources to subtract
        ## convert Jy to milliJy
        subtract_index = np.where( result['Total_flux'] > subtract_limit*1e3 )[0]
        subtract_cals = result[['Source_id','LOTSS_RA','LOTSS_DEC']][subtract_index]
	#subtract_cals.rename_column('Source_ID','Source_id')
	subtract_cals.rename_column('LOTSS_RA','RA')
	subtract_cals.rename_column('LOTSS_DEC','DEC')
	
        ## also include bright sources from the lotss catalogue
        ## convert Jy to milliJy
        bright_index = np.where( lotss_catalogue['Total_flux'] >= bright_limit_Jy*1e3 )[0]
        subtract_bright = lotss_catalogue[['Source_id','RA','DEC']][bright_index]
	subtract_bright['Source_id'] = subtract_bright['Source_id'].astype(str)
	## change the data type of the source_id column
	subtract_bright['tmp'] = subtract_bright['Source_id'].astype(str)
	subtract_bright.remove_column('Source_id')
	subtract_bright.rename_column('tmp','Source_id')
	
	subtract_sources = vstack( [subtract_cals, subtract_bright])
	subtract_sources = unique( subtract_sources )
	subtract_sources.write( subtract_file, format='csv' )

        ## sources to image -- first remove things that are already in the delay_cals_file and subtract_file
	good_index = [ x for x, src_id in enumerate( lotss_catalogue['Source_id'] ) if src_id not in result['Source_id'] and src_id not in subtract_cals['Source_id'] ]
	
	tmp_cat = lotss_catalogue[good_index]

	## make a flux cut
	image_index = np.where( tmp_cat['Peak_flux'] >= image_limit_Jy*1e3 )[0]
	sources_to_image = tmp_cat[image_index]

        ## find unresolved
        nsrcs = float( len( sources_to_image ) )
	print "There are "+str(nsrcs)+" sources above "+str(image_limit_Jy)+" mJy."
        unresolved_index = np.where( sources_to_image['Resolved'] == 'U' )[0]
        if nsrcs==0:
            print "Did not find any unresolved objects."
        else:
            perc_unres = len( unresolved_index ) / nsrcs * 100.
            print 'Percentage of sources which are unresolved: '+str( perc_unres )

	sources_to_image.write( lotss_result_file, format='csv' )

    else:
	print "doDownload == False, not downloading catalogues!"


    return

 def is_resolved(Sint, Speak, rms):
    """ Determines if a source is resolved or unresolved.
    The calculation is presented in Shimwell et al. 2018 of the LOFAR DR1 paper splash.
    
    Args:
        Sint (float or ndarray): integrated flux density.
        Speak (float or ndarray): peak flux density.
        rms (float or ndarray): local rms around the source.
    Returns:
        resolved (bool or ndarray): True if the source is resolved, False if not.
    """
    resolved = ((Sint / Speak) > 1.25 + 3.1 * (Speak / rms) ** (-0.53))
    return resolved

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( '--lotss_radius', dest='lotss_radius', type=float, help='Radius to search LoTSS', default=5. )
    parser.add_argument( '--lbcs_radius', dest='lbcs_radius', type=float, help='Radius to search LBCS', default=5. )
    parser.add_argument( '--lotss_result_file', dest='lotss_result_file', type=str, help='output file of sources to image', default='lotss_catalogue.csv' )
    parser.add_argument( '--delay_cals_file', dest='delay_cals_file', type=str, help='output file of delay calibrators', default='delay_calibrators.csv' )
    parser.add_argument( '--subtract_file', dest='subtract_file', type=str, help='output file of sources to subtract', default='subtract_sources.csv' )
    parser.add_argument( '--match_tolerance', dest='match_tolerance', type=float, help='radius for matching LBCS to LoTSS [arcsec]', default=5. )
    parser.add_argument( '--subtract_limit', dest='subtract_limit', type=float, help='Flux limit for sources to subtract [Jy]', default=0.5 )
    parser.add_argument( '--bright_limit_Jy', dest='bright_limit_Jy', type=float, help='Flux limit for bright sources [Jy]', default=5.0 )
    parser.add_argument( '--image_limit_Jy', dest='image_limit_Jy', type=float, help='Flux limit for which sources to image [Jy]', default = 0.05 )
    parser.add_argument( 'MSname', type=str, nargs='+', help='Measurement set name (to get pointing center)' )

    args = parser.parse_args()

    plugin_main( args.MSname, lotss_radius=args.lotss_radius, lbcs_radius=args.lbcs_radius, bright_limit_Jy=args.bright_limit_Jy, lotss_result_file=args.lotss_result_file, delay_cals_file=args.delay_cals_file, subtract_file=args.subtract_file, match_tolerance=args.match_tolerance, subtract_limit=args.subtract_limit, image_limit_Jy=args.image_limit_Jy )

