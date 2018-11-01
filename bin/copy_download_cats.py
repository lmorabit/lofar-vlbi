from scipy import stats  
import os, sys, logging
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from lmfit import SkewedGaussianModel
import pyvo as vo
import pyrap.tables as pt
import pandas
import argparse

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
    url = 'http://vo.astron.nl/lofartier1/q/cone/scs.xml'

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
    tb_final = tb_sorted['Source_id', 'RA', 'DEC','Total_flux','Peak_flux','Isl_rms','Resolved']

    return tb_final

def my_lbcs_catalogue( ms_input, Radius=1.5 ):

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
    Source_ID=np.array(lo['Source_id'])
    rms=np.array(lo['Isl_rms'])
    resolved=np.array(lo['Resolved'])
    lbcsid=np.array(lb['ObsID'])

    result = [] 

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
                resultarray=[Source_ID[c_idx][0],lbcsRA[x],lbcsDEC[x],lotssRA[c_idx][0],lotssDEC[c_idx][0],np.min(delta_position),lotssflux[c_idx][0],lotsspeakflux[c_idx][0],(lotsspeakflux[c_idx][0]/lotssflux[c_idx][0]), rms[c_idx][0], resolved[c_idx][0] ]
                result.append(resultarray)
            else:
                print "No match found in LoTSS for LBCS calibrator %s"%(lbcsid[x])
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
            resultarray=[Source_ID[c_idx][0],lbcsRA,lbcsDEC,lotssRA[c_idx][0],lotssDEC[c_idx][0],np.min(delta_position),lotssflux[c_idx][0],lotsspeakflux[c_idx][0],(lotsspeakflux[c_idx][0]/lotssflux[c_idx][0]), rms[c_idx][0], resolved[c_idx][0] ]
            result.append(resultarray)
        else:
            print "No match found in LoTSS for LBCS calibrator %s"%(lbcsid)

    ## convert result to numpy array
    #result = np.array(result)
    return result

def main( MSname, lotss_radius=5., lbcs_radius=5., bright_limit_Jy=5., lotss_result_file='lotss_catalogue.csv', delay_cals_file='delay_calibrators.csv', match_tolerance=5. ):

    lotss_catalogue = my_lotss_catalogue( MSname, Radius=lotss_radius, bright_limit_Jy=bright_limit_Jy )
    lbcs_catalogue = my_lbcs_catalogue( MSname, Radius=lbcs_radius )
    result = find_close_objs( lotss_catalogue, lbcs_catalogue, tolerance=match_tolerance )
    ## result columns are ['Source_ID','LBCS_RA','LBCS_DEC','LOTSS_RA','LOTSS_DEC','Delta_Position','total_flux','Peak_flux','peak_to_total','Isl_rms','Resolved']

    ## Need to write the following catalogues:
    ## 1 - delay calibrators -- from lbcs_catalogue
    colnames = 'Source_ID,LBCS_RA,LBCS_DEC,LOTSS_RA,LOTSS_DEC,Delta_Position,total_flux,Peak_flux,peak_to_total,Isl_rms,Resolved'
    np.savetxt( delay_cals_file, result, fmt='%s', delimiter=',', header=colnames )
 
    ## find unresolved
    nsrcs = float( len( lotss_catalogue.columns['Resolved'] ) )
    unresolved_index = np.where( tb.columns['Resolved'] == 'U' )[0]
    if nsrcs==0:
        print "No-founded objs in this radius!!"
    else:
        perc_unres = len( unresolved_index ) / nsrcs * 100.
        print 'Percentage of sources which are unresolved: '+str( perc_unres )

    ## also find very bright things
    ## have to convert Jy to milliJy
    bright_index = np.where( tb.columns['Total_flux'] >= bright_limit_Jy*1e3 )[0]
    if len(bright_index) == 0:
	print "No sources are found above %f Jy"%(bright_limit_Jy)
    else:
	btb = tb[bright_index]
	unresolved_index = np.append( unresolved_index, bright_index )
	unresolved_index = np.unique( unresolved_index )

    utb = tb[unresolved_index]



    ## 2 - sources to subtract -- initially delay calibrators + extremely bright sources
    ## 3 - sources to image
   
    return

