#!/usr/bin/env python
#from scipy import stats  
import os, sys, logging, io
import numpy as np
#import csv
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from lmfit import SkewedGaussianModel
import pyvo as vo
import pyrap.tables as pt
from astropy.table import Table, Column, vstack, unique, hstack
import argparse
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import requests
from astropy.coordinates import SkyCoord
from astropy import units as u


def sum_digits(n):
    s = 0
    c = 0
    for ii in n:
        if ii != '-':
            s += int(ii)
            c += 1
    norm = float(s)/float(c)
    return( norm )

def count_p(n):
    s = 0
    c = 0
    for ii in n:
	if ii != '-':
	    c += 1
	    if ii == 'P':
		s += 1
    norm = float(s)/float(c)
    return(norm)

def count_s(n):
    s = 0
    c = 0
    for ii in n:
        if ii != '-':
            c += 1
            if ii == 'S':
                s += 1
    norm = float(s)/float(c)
    return(norm)

def count_x(n):
    s = 0
    c = 0
    for ii in n:
        if ii != '-':
            c += 1
            if ii == 'X':
                s += 1
    norm = float(s)/float(c)
    return(norm)

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

def my_lotss_catalogue( ms_input, Radius=1.5, bright_limit_Jy=5., outfile='' ):

    """
    Download the LoTSS skymodel for the target field
    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    Radius : float (default = 1.5)
        Radius for the LOTSS cone search in degrees
    
    """
    ## first check if the file already exists
    if os.path.isfile( outfile ):
	print("LOTSS Skymodel for the target field exists on disk, reading in.")
	tb_final = Table.read( outfile, format='csv' )
        if 'Source_Name' in tb_final.colnames:
            tb_final.rename_column('Source_Name','Source_id')
            tb_final.write( outfile, format='csv' )
    else:
        print("DOWNLOADING LOTSS Skymodel for the target field")

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
        try:
            tb = t.votable.to_table()
        except AttributeError:
            # Above statement didn't work, try the alternative.
            tb = t.to_table()
        flux_sort = tb.argsort('Total_flux')
        tb_sorted = tb[flux_sort[::-1]]
        ## and keep only some of the columns
        colnames = tb_sorted.colnames
        ## first check for a resolved column
        if 'Resolved' not in colnames:
            resolved = np.where(is_resolved(tb_sorted['Total_flux'], tb_sorted['Peak_flux'], tb_sorted['Isl_rms']), 'R', 'U')
            tb_sorted['Resolved'] = resolved
        ## rename source_id column if necessary   
        if 'Source_Name' in tb_sorted.colnames:
            tb_sorted.rename_column('Source_Name', 'Source_id')
        keep_cols = ['Source_id', 'RA', 'DEC','Total_flux','Peak_flux', 'Major', 'Minor', 'PA', 'DC_Maj', 'DC_Min', 'DC_PA', 'Isl_rms', 'Resolved']
        if 'LGZ_Size' in colnames:
            keep_cols = keep_cols + ['LGZ_Size', 'LGZ_Width', 'LGZ_PA']

        tb_final = tb_sorted[keep_cols]

        tb_final.write( outfile, format='csv' )

    return tb_final

def my_lbcs_catalogue( ms_input, Radius=1.5, outfile='' ):

    """
    Download the LBCS skymodel for the target field
    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    Radius : float (default = 1.5)
        Radius for the LOTSS cone search in degrees
    
    """
    ## first check if the file already exists
    print( outfile )
    if os.path.isfile( outfile ):
        print("LBCS Skymodel for the target field exists on disk, reading in.")
        tb = Table.read( outfile, format='csv' )
    else:
        print("DOWNLOADING LBCS Skymodel for the target field")

        # Reading a MS to find the coordinate (pyrap)
        RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

	## construct an html query and try to connect
	url = 'https://lofar-surveys.org/lbcs-search.fits?ra=%f&dec=%f&radius=%f' % (float(RATar), float(DECTar), float(Radius))
        connected=False
        while not connected:
            try:
                response = requests.get(url, stream=True,verify=True,timeout=60)
                if response.status_code!=200:
                    print response.headers
                    raise RuntimeError('Code was %i' % response.status_code)
            except requests.exceptions.ConnectionError:
                print 'Connection error! sleeping 30 seconds before retry...'
                sleep(30)
            except (requests.exceptions.Timeout,requests.exceptions.ReadTimeout):
                print 'Timeout! sleeping 30 seconds before retry...'
                sleep(30)
            else:
                connected=True

        mem = io.BytesIO()
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                mem.write(chunk)
        mem.seek(0)
        tb = Table.read(mem)
        mem.close()
        del(response)

        print( len( tb ) )

        #### Workaround to sort and pick good calibrator info from tb array ###########
    
        if not len(tb) > 1:
            logging.critical('There are no LBCS sources within the given radius. Check your source is within the LBCS footprint and increase the search radius. Exiting...')
            return
	else:
            ## calculate the total FT goodness
            ft_total = []
	    for xx in range(len(tb)):
                ft_total.append( sum_digits( tb[xx]['FT_Goodness'] ) )
            ft_col = Column( ft_total, name='FT_total' )
            tb.add_column( ft_col )

        tb.write( outfile, format='csv' )

    return tb

def find_close_objs(lo, lbcs, tolerance=5.):

    ## first filter the LBCS data on Flags
    lbcs_idx = np.where( np.logical_or(lbcs['Flags'] == 'O', lbcs['Flags'] == 'A' ) )
    lbcs = lbcs[lbcs_idx]
    ## get rid of anything with only X's
    nump = []
    nums = []
    for xx in range(len(lbcs)):
        nump.append(count_p(lbcs['Goodness'][xx]))
        nums.append(count_s(lbcs['Goodness'][xx]))
    print( np.array(nump) + np.array(nums) )
    lbcs_idx = np.where( np.array(nump)+np.array(nums) > 0 )[0]
    lbcs = lbcs[lbcs_idx]

    ## get RA and DEC for the catalogues
    lotss_coords = SkyCoord( lo['RA'], lo['DEC'], frame='icrs', unit='deg' )
    lbcs_coords = SkyCoord( lbcs['RA'], lbcs['DEC'], frame='icrs', unit='deg' )

    ## search radius 
    search_rad = 5. / 60. / 60. * u.deg

    ## loop through the lbcs coordinates -- this will be much faster than looping through lotss
    lotss_idx = []
    lbcs_idx = []    
    for xx in range(len(lbcs)):
        seps = lbcs_coords[xx].separation(lotss_coords)
        match_idx = np.where( seps < search_rad )[0]
	if len( match_idx ) == 0:
            # there's no match, move on to the next source
            m_idx = [-1]
            pass
        else:
	    if len( match_idx ) == 1:
                ## there's only one match
                m_idx = match_idx[0]
                lbcs_idx.append(xx)
                lotss_idx.append(m_idx)
            if len( match_idx ) > 1:
                ## there's more than one match, pick the brightest
                tmp = lo[match_idx]
                m_idx = np.where( tmp['Total_flux'] == np.max( tmp['Total_flux'] ) )[0]
                if not isinstance(m_idx,int):
		    m_idx = m_idx[0]
                lbcs_idx.append(xx)
                lotss_idx.append(m_idx) 
    lbcs_matches = lbcs[lbcs_idx]
    lotss_matches = lo[lotss_idx]

    combined = hstack( [lbcs_matches, lotss_matches], join_type='exact' )
    ## check if there are duplicate lbcs observations for a lotss source
    if len( np.unique( combined['Source_id'] ) ) != len( combined ):
        # there are duplicates
        print( 'There are duplicate LBCS sources, selecting the best candidate(s).' )
        src_ids = np.unique( combined['Source_id'] )
        good_idx = []
        for src_id in src_ids:
            idx = np.where( combined['Source_id'] == src_id )[0]
	    if len(idx) > 1:
		## multiple matches found.  Count P's first and then break ties with Goodness_FT 
		num_P = []
	        total_ft = []
		for yy in range( len( idx ) ):
		    tmp = combined[idx[yy]]['Goodness']
		    num_P.append( count_p( tmp ) )

		    tmp = combined[idx[yy]]['FT_Goodness']
                    total_ft.append( sum_digits( tmp ) )
                ## check that the total_ft values are non-zero before looking for a best cal
                if np.max( total_ft ) > 0:
	            ## pick the one with the highest number of P's -- if tie, use total_ft
		    best_idx = np.where( num_P == np.max( num_P ) )[0]  ## this is an array
	            if len( best_idx ) == 1:
		        good_idx.append(idx[best_idx][0])  ## idx[best_idx][0] is a number
		    if len( best_idx ) > 1:
		        currentmax = 0.0 
                        for i in range(0,len(best_idx)):
		            if total_ft[best_idx[i]] > currentmax:
			        currentmax = total_ft[best_idx[i]]
                                ft_idx = i
		        good_idx.append( idx[best_idx[ft_idx]] )
                else:
                    print( 'Duplicate sources have total_ft = 0, removing from results.' )
            else:
		good_idx.append(idx[0])

	result = combined[good_idx]
    else:
	print( 'No duplicate sources found' )
        result = combined
    ## rename RA columns
    result.rename_column('RA_1','RA_LBCS')
    result.rename_column('DEC_1','DEC_LBCS')
    result.rename_column('RA_2','RA_LOTSS')
    result.rename_column('DEC_2','DEC_LOTSS')

    return result

def plugin_main( args, **kwargs ):

    mapfile_in = kwargs['mapfile_in']
    lotss_radius = kwargs['lotss_radius']
    lbcs_radius  = kwargs['lbcs_radius']
    bright_limit_Jy = float(kwargs['bright_limit_Jy'])
    lotss_result_file = kwargs['lotss_result_file']
    lotss_catalogue = kwargs['lotss_catalogue']
    lbcs_catalogue = kwargs['lbcs_catalogue']
    delay_cals_file = kwargs['delay_cals_file']
    subtract_file = kwargs['subtract_file']
    match_tolerance = float(kwargs['match_tolerance'])
    subtract_limit = float(kwargs['subtract_limit'])
    image_limit_Jy = float(kwargs['image_limit_Jy'])
    fail_lotss_ok = kwargs['continue_no_lotss'].lower().capitalize()

    mslist = DataMap.load(mapfile_in)
    MSname = mslist[0].file
    # For testing
    #MSname = kwargs['MSname']
 
    ## first check for a valid delay_calibrator file
    if os.path.isfile(delay_cals_file):
	print( 'Delay calibrators file {:s} exists! returning.'.format(delay_cals_file) )
        return

    ## look for or download LBCS
    print("Attempting to find or download LBCS catalogue.")
    lbcs_catalogue = my_lbcs_catalogue( MSname, Radius=lbcs_radius, outfile=lbcs_catalogue )
    ## look for or download LoTSS
    print("Attempting to find or download LoTSS catalogue.")
    lotss_catalogue = my_lotss_catalogue( MSname, Radius=lotss_radius, bright_limit_Jy=bright_limit_Jy, outfile=lotss_catalogue )

    ## if lbcs exists, and either lotss exists or continue_without_lotss = True, process the catalogue(s).
    ## else provide an error message and stop
    if len(lbcs_catalogue) == 0:
	logging.error('LBCS coverage does not exist, and catalogue not found on disk.')
	return
    if len(lotss_catalogue) == 0 and not fail_lotss_ok:
	logging.error('LoTSS coverage does not exist, and contine_without_lotss is set to False.')
	return 

    ## if the LoTSS catalogue is empty, write out the delay cals only and stop
    if len(lotss_catalogue) == 0:
        print('Target field not in LoTSS coverage yet! Only writing {:s} and best_{:s} based on LBCS'.format(delay_cals_file,delay_cals_file))
        lbcs_catalogue.write(delay_cals_file, format='csv')
        ## pick the best calibrator based on LBCS information only
        RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(MSname)[0])
        ptg_coords = SkyCoord( RATar, DECTar, frame='icrs', unit='deg' )

        src_coords = SkyCoord( lbcs_catalogue['RA'], lbcs_catalogue['DEC'], frame='icrs', unit='deg' )
        separations = src_coords.separation(ptg_coords )
        seps = Column( separations.deg, name='Radius' )
        lbcs_catalogue.add_column( seps )

        ## highest FT_total; if tie use closest to phase centre
        best_idx = np.where( lbcs_catalogue['FT_total'] == np.max( lbcs_catalogue['FT_total'] ) )[0]
        if len( best_idx ) > 1:
            tmp = seps[best_idx]
            new_best_idx = np.where( tmp == np.min(tmp) )[0]
            best_idx = best_idx[new_best_idx]

        ## rename some columns 
        lbcs_catalogue.rename_column('RA','RA_LOTSS')
        lbcs_catalogue.rename_column('DEC','DEC_LOTSS')
        lbcs_catalogue.rename_column('Observation','Source_id')
        ## add in some dummy data
        Total_flux = Column( np.ones(len(lbcs_catalogue)), name='Total_flux' )
        lbcs_catalogue.add_column( Total_flux )
        LGZ_Size = Column( np.ones( len(lbcs_catalogue) )*20., name='LGZ_Size' ) ## set to a default of 20 arcsec
        lbcs_catalogue.add_column( LGZ_Size )

        best_result = lbcs_catalogue[best_idx]
        best_file = delay_cals_file.replace('delay_','best_delay_')
        best_result.write( best_file, format='csv' )
        print( 'Writing best delay calibrator information to file {:s}'.format(best_file) )
        return

    ## else continue 
    result = find_close_objs( lotss_catalogue, lbcs_catalogue, tolerance=match_tolerance )

    ## check if there are any matches
    if len(result) == 0:
        logging.error('LoTSS and LBCS coverage exists, but no matches found. This indicates something went wrong, please check your catalogues.')
        return
    else:
	# find the best delay calibrator
        RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(MSname)[0])
        ptg_coords = SkyCoord( RATar, DECTar, frame='icrs', unit='deg' )

	src_coords = SkyCoord( result['RA_LBCS'], result['DEC_LBCS'], frame='icrs', unit='deg' )
	separations = src_coords.separation(ptg_coords )
        seps = Column( separations.deg, name='Radius' )
        result.add_column( seps )

	radsq = np.array( result['Radius'] )**2.
        fluxjy = np.array( result['Total_flux'] )*1e-3
        ftsq = np.array( result['FT_total'] )**2.

	mystat = [ radsq[xx]/fluxjy[xx]/ftsq[xx] if ftsq[xx] > 0 else 100 for xx in range(len(radsq)) ]
        mycol = Column( mystat, name='Select_stat' )
	result.add_column( mycol )

	## pick the best one
	best_idx = np.where( mystat == np.min( mystat ) )[0]
	tmp = result['Source_id'][best_idx].data
	tmp2 = result['Select_stat'][best_idx].data
	print( 'Best delay calibrator candidate is {:s} with a stastic of {:f}'.format(str(tmp[0]),tmp2[0]) )

        ## Write catalogues
        ## 1 - delay calibrators -- from lbcs_catalogue
        result.write( delay_cals_file, format='csv' )
	print('Writing delay calibrator candidate file {:s}'.format(delay_cals_file))

	## best delay calibrator
	best_result = result[best_idx]
	best_file = delay_cals_file.replace('delay_','best_delay_')
        best_result.write( best_file, format='csv' )
        print( 'Writing best delay calibrator information to file {:s}'.format(best_file) )

        ## convert Jy to milliJy
        subtract_index = np.where( result['Total_flux'] > subtract_limit*1e3 )[0]
        subtract_cals = result[subtract_index]
	
        ## also include bright sources from the lotss catalogue
        ## convert Jy to milliJy
        bright_index = np.where( lotss_catalogue['Total_flux'] >= bright_limit_Jy*1e3 )[0]
        lotss_bright = lotss_catalogue[bright_index]
	## lotss catalogue has units, redefine a table that doesn't
	
	subtract_sources = vstack( [subtract_cals, lotss_bright], join_type='outer' )

	unique_srcs = np.unique( subtract_sources['Source_id'] )
	if len( unique_srcs ) != len( subtract_sources ):
	    ## remove duplicates, keep LBCS information
            ## this is untested and will probably fail ...
	    good_idx = []
	    for src_id in unique_srcs:
		idx = np.where( subtract_sources['Source_id'] == src_id )[0]
		if len( idx ) > 1:
		    tmp = subtract_sources[idx]
		    lbcs_idx = np.where( tmp['RA_LBCS'] > 0 )[0]
		    good_idx.append( idx[lbcs_idx][0] )
		else:
		    good_idx.append( idx[0] )
	    subtract_sources = subtract_sources[good_idx]

	subtract_sources.write( subtract_file, format='csv' )

        ## sources to image -- first remove things that are already in the delay_cals_file and subtract_file
	good_index = [ x for x, src_id in enumerate( lotss_catalogue['Source_id'] ) if src_id not in result['Source_id'] and src_id not in subtract_sources['Source_id'] ]
	
	tmp_cat = lotss_catalogue[good_index]

	## make a flux cut
	image_index = np.where( tmp_cat['Peak_flux'] >= image_limit_Jy*1e3 )[0]
	sources_to_image = tmp_cat[image_index]

        ## find unresolved
        nsrcs = float( len( sources_to_image ) )
	print "There are "+str(nsrcs)+" sources above "+str(image_limit_Jy)+" mJy."
        try:
            unresolved_index = np.where( sources_to_image['Resolved'] == 'U' )[0]
	except:
	    unresolved_index = np.where( is_resolved(sources_to_image['Total_flux'],  sources_to_image['Peak_flux'], sources_to_image['Isl_rms'] ) )[0]
        if nsrcs==0:
            print "Did not find any unresolved objects."
        else:
            perc_unres = len( unresolved_index ) / nsrcs * 100.
            print 'Percentage of sources which are unresolved: '+str( perc_unres )

	sources_to_image.write( lotss_result_file, format='csv' )

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
    parser.add_argument( '--lotss_catalogue', dest='lotss_catalogue', type=str, help='input file for LoTSS catalogue [will be downloaded if does not exist]', default='lotss_catalogue.csv' )
    parser.add_argument( '--lbcs_catalogue', dest='lbcs_catalogue', type=str, help='input file for LBCS catalogue [will be downloaded if does not exist]', default='lbcs_catalogue.csv' )
    parser.add_argument( '--lotss_result_file', dest='lotss_result_file', type=str, help='output file of sources to image', default='image_catalogue.csv' )
    parser.add_argument( '--delay_cals_file', dest='delay_cals_file', type=str, help='output file of delay calibrators', default='delay_calibrators.csv' )
    parser.add_argument( '--subtract_file', dest='subtract_file', type=str, help='output file of sources to subtract', default='subtract_sources.csv' )
    parser.add_argument( '--match_tolerance', dest='match_tolerance', type=float, help='radius for matching LBCS to LoTSS [arcsec]', default=5. )
    parser.add_argument( '--subtract_limit', dest='subtract_limit', type=float, help='Flux limit for sources to subtract [Jy]', default=0.5 )
    parser.add_argument( '--bright_limit_Jy', dest='bright_limit_Jy', type=float, help='Flux limit for bright sources [Jy]', default=5.0 )
    parser.add_argument( '--image_limit_Jy', dest='image_limit_Jy', type=float, help='Flux limit for which sources to image [Jy]', default = 0.05 )
    parser.add_argument( '--continue_no_lotss', dest='continue_no_lotss', action='store_true', help='Continue with the pipeline if no lotss cross-matches can be found?', default = False )
    parser.add_argument( 'MSname', type=str, nargs='+', help='Measurement set name (to get pointing center)' )

    args = parser.parse_args()

    plugin_main( args.MSname, MSname=args.MSname, lotss_radius=args.lotss_radius, lbcs_radius=args.lbcs_radius, bright_limit_Jy=args.bright_limit_Jy, lotss_catalogue=args.lotss_catalogue, lbcs_catalogue=args.lbcs_catalogue, lotss_result_file=args.lotss_result_file, delay_cals_file=args.delay_cals_file, subtract_file=args.subtract_file, match_tolerance=args.match_tolerance, subtract_limit=args.subtract_limit, image_limit_Jy=args.image_limit_Jy, continue_no_lotss = str(args.continue_no_lotss) )

