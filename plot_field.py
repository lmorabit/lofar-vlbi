#!/usr/bin/env python
import os, sys, logging, io
import numpy as np
import pyvo as vo
import casacore.tables as pt
from astropy.table import Table, Column, vstack, unique, hstack
import argparse
#from lofarpipe.support.data_map import DataMap
#from lofarpipe.support.data_map import DataProduct
import requests
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from time import sleep
from glob import glob


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

def my_lotss_catalogue( RATar, DECTar,  Radius=1.5, bright_limit_Jy=5., faint_limit_Jy = 0.0, outfile='' ):

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
        #RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    ## this is the tier 1 database to query
        #url = 'http://vo.astron.nl/lofartier1/q/cone/scs.xml'
        # HETDEX database.
        #url = 'https://vo.astron.nl/hetdex/lotss-dr1/cone/scs.xml'
        url = 'https://vo.astron.nl/lotss_dr2/q/src_cone/scs.xml'

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
        keep_cols = ['Source_id', 'RA', 'DEC','Total_flux','Peak_flux', 'Majax', 'Minax', 'PA', 'DC_Maj', 'DC_Min', 'DC_PA', 'Isl_rms', 'Resolved']
        if 'LGZ_Size' in colnames:
            keep_cols = keep_cols + ['LGZ_Size', 'LGZ_Width', 'LGZ_PA']

        tb_final = tb_sorted[keep_cols]

        tb_final = tb_final[tb_final['Total_flux'] >= faint_limit_Jy*1e3]

        tb_final = tb_final[tb_final['Total_flux'] <= bright_limit_Jy*1e3]

        tb_final.write( outfile, format='csv' )

    return tb_final

def my_lbcs_catalogue( RATar, DECTar, Radius=1.5, outfile='' ):

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
    #print( outfile )
    if os.path.isfile( outfile ):
        print("LBCS Skymodel for the target field exists on disk, reading in.")
        tb = Table.read( outfile, format='csv' )
        return(tb)
    else:
        print("DOWNLOADING LBCS Skymodel for the target field")

        # Reading a MS to find the coordinate (pyrap)
        #RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    ## construct an html query and try to connect
    #print(RATar, DECTar)
    url = 'https://lofar-surveys.org/lbcs-search.fits?ra=%f&dec=%f&radius=%f' % (float(RATar), float(DECTar), float(Radius))
    connected=False
    while not connected:
            try:
                response = requests.get(url, stream=True,verify=True,timeout=60)
                if response.status_code!=200:
                    print (response.headers)
                    raise RuntimeError('Code was %i' % response.status_code)
            except requests.exceptions.ConnectionError:
                print ('Connection error! sleeping 30 seconds before retry...')
                sleep(30)
            except (requests.exceptions.Timeout,requests.exceptions.ReadTimeout):
                print( 'Timeout! sleeping 30 seconds before retry...')
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

    #print( len( tb ) )

    ### Workaround to sort and pick good calibrator info from tb array ###########
    
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
    #print( np.array(nump) + np.array(nums) )
    lbcs_idx = np.where( np.array(nump)+np.array(nums) > 0 )[0]
    lbcs = lbcs[lbcs_idx]

    ## get RA and DEC for the catalogues
    lotss_coords = SkyCoord( lo['RA'], lo['DEC'], frame='icrs', unit='deg' )
    lbcs_coords = SkyCoord( lbcs['RA'], lbcs['DEC'], frame='icrs', unit='deg' )

    ## search radius 
    search_rad = tolerance / 60. / 60. * u.deg

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
    result.rename_column('RA_2','RA')
    result.rename_column('DEC_2','DEC')

    return result

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

def remove_multiples_position( mycat, racol='RA', decol='DEC' ):
    radecstrings = []
    for i in np.arange(0,len(mycat)):
        radecstrings.append(str(mycat[i][racol]) + str(mycat[i][decol]) )
    radecstrings = np.asarray(radecstrings)
    if len( np.unique( radecstrings ) ) != len( mycat ):
        radecs = np.unique( radecstrings )
        good_idx = []
        for radec in radecs:
            idx = np.where( radecstrings == radec )[0]
            if len(idx) > 1:
                ## multiple matches found.  Count P's first and then break ties with Goodness_FT 
                num_P = []
                total_ft = []
                for yy in range( len( idx ) ):
                    tmp = mycat[idx[yy]]['Goodness']
                    num_P.append( count_p( tmp ) )
                    tmp = mycat[idx[yy]]['FT_Goodness']
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
        mycat = mycat[good_idx]
    else:
        print( 'All LBCS sources are unique' )

    return( mycat )

def mkfits (rasiz,decsiz,imsiz,pixsiz):
    hdu=fits.PrimaryHDU(np.zeros((int(imsiz),int(imsiz))))
    hdu.header.update({'CTYPE1':'RA---SIN'})
    hdu.header.update({'CRVAL1':rasiz})
    hdu.header.update({'CRPIX1':imsiz/2.})
    hdu.header.update({'CDELT1':-pixsiz})
    hdu.header.update({'CTYPE2':'DEC--SIN'})
    hdu.header.update({'CRVAL2':decsiz})
    hdu.header.update({'CRPIX2':imsiz/2.})
    hdu.header.update({'CDELT2':pixsiz})
    hdu.header.update({'EQUINOX':2000.0})
    hdu.header.update({'EPOCH':2000.0})
#    hdu.data = np.random.random(imsiz*imsiz).reshape(imsiz,imsiz)
    if os.path.exists('temp.fits'):
        os.system('rm temp.fits')
    hdu.writeto('temp.fits')
    
def smearing_calculation(nchan = 16, obs_freq=144000000, radii=np.arange(0.001,4,0.00001), resolution=0.3/206265 * 180 / np.pi,
              av_time=1, thresholds = [0.2, 0.4, 0.6, 0.8] ):
    from scipy.special import erf

    bandwidth=1.95e3/nchan

    beta_fac = bandwidth/obs_freq * radii / resolution
    gamma_fac = 2. * np.sqrt( np.log(2.) ) 
    reduction_bandwidth = np.sqrt(np.pi)/(gamma_fac*beta_fac) * erf( gamma_fac*beta_fac/2)
    reduction_time = 1 - 1.22e-9 * ( radii / resolution )**2 * av_time**2
    reduction_total = reduction_bandwidth*reduction_time
    #print(reduction_bandwidth, reduction_total)
    thres_radii = []
    for thresh in thresholds:
        idx = np.argwhere(reduction_total < thresh)[0]
        #print(radii[idx])
        thres_radii.append(radii[idx])
    return thres_radii


def angular_distance(RA1, DEC1, RA2, DEC2):
    return np.sqrt((RA1-RA2)**2 + (DEC1-DEC2)**2)

def smallest_distance(RA, DEC, lbcs_catalogue):
    distances = [angular_distance(RA, DEC, source['RA'], source['DEC']) for source in lbcs_catalogue]

    ids = np.argmin(distances)

    return ids, distances[ids]
    
    
def make_plot(RA, DEC,  lotss_catalogue, extreme_catalogue, lbcs_catalogue, targRA=None, 
          targDEC=None, nchan = 16, av_time = 1, outdir='.'):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    mkfits(RA, DEC, 2048, 2.6367)

    f = fits.open("temp.fits")
    w = WCS(f[0].header)

    lbcs = lbcs_catalogue

    if os.path.exists(lotss_catalogue):
        lotss = Table.read(lotss_catalogue)
        avg_flux = np.median(lotss['Total_flux'])
    scaling = 0.1

    thres_radii = smearing_calculation(nchan = nchan, av_time=av_time)

    plt.figure(figsize = (10,10))

    ax = plt.subplot(projection = w,)

    color = 0.5
    fraction = ["80%" , "60%", "40%", "20%"]
    for i,thresh in enumerate(thres_radii):
        color = color + 0.1
        c = Circle((RA, DEC), thresh, edgecolor=None, facecolor=str(color),
           transform=ax.get_transform('fk5'), zorder = -1)
        ax.text(RA + 0.98*thresh , DEC, fraction[i], 
              transform=ax.get_transform('fk5'), fontsize = "large")
        ax.add_patch(c)

    ax.scatter(RA, DEC, marker= 'x', s = 30, transform=ax.get_transform('fk5'), label = "Ptg Centre")

    if os.path.exists(lotss_catalogue):
        ax.scatter(lotss['RA'], lotss['DEC'], transform=ax.get_transform('fk5'), 
            s = lotss['Total_flux'] * scaling, label = "Potential Targets")
        
    if len(extreme_catalogue) > 0:
        for source in extreme_catalogue:
            vector_orig = [source['RA'] - RA, source['DEC'] - DEC]
            norm = np.sqrt(vector_orig[0]**2 + vector_orig[1]**2)
            vector = [vector_orig[0]/norm, vector_orig[1]/norm]
            alt = [vector[0] -RA, vector[1] - DEC]
            print(vector, alt)
            ax.arrow(RA+2*vector[0], DEC+2*vector[1], 0.5*vector[0], 0.5*vector[1], transform=ax.get_transform('fk5'), 
                width = 0.01, head_width = 0.05, head_length = 0.05, length_includes_head = True) 
            ax.text(RA+2.5*vector[0], DEC+2.5*vector[1], 
                    s = "%.2f"%(source['Total_flux']/1000) + " Jy - %.2f degrees"%norm, 
                    transform=ax.get_transform('fk5'))

    ax.scatter(lbcs['RA'], lbcs['DEC'], transform=ax.get_transform('fk5'), s = 60, label = "LBCS Sources")
    

    c = Circle((RA, DEC), 1.5, edgecolor='yellow', facecolor='none',
           transform=ax.get_transform('fk5'))
    ax.add_patch(c)


    ### Calculate and plot smallest distance to target 
    dist_ids, dist = smallest_distance(RA, DEC, lbcs)

    closest_calib = lbcs[dist_ids]
    ax.plot([RA, closest_calib['RA']], [DEC, closest_calib['DEC']], 
             linestyle = '--', linewidth = 2, transform=ax.get_transform('fk5'), label = "Distance = %.2f degrees"%dist)    

    if targRA != None:
        ax.scatter(targRA, targDEC, marker= 'x', s = 80, 
         transform=ax.get_transform('fk5'), label = "Target")
        dist = angular_distance(targRA, targDEC, closest_calib['RA'], closest_calib['DEC'])

        ax.plot([targRA, closest_calib['RA']], [targDEC, closest_calib['DEC']], 
                 linestyle = '--', color = 'black',  linewidth = 2, transform=ax.get_transform('fk5'), label = "Distance = %.2f degrees"%dist)
    
    for i in range(len(lbcs)):
        ax.text(lbcs[i]['RA']-0.05, lbcs[i]['DEC'], "%.2f Jy"%(lbcs[i]['Total_flux']/1000),  transform=ax.get_transform('fk5'))
    plt.legend(fontsize = 'x-large')

    plt.savefig(os.path.join(outdir,"output.png"))
    os.system('rm temp.fits')

def convert_vlass_fits(fitsfile):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from astropy.visualization import PercentileInterval, imshow_norm
    from astropy.wcs import WCS
    header = fits.open(fitsfile)[0].header
    wcs = WCS(header)

    image_data = fits.getdata(fitsfile)

    # Shape is (1,1, 3722, 3722). Plot the first image
    interval = PercentileInterval(99.9)
    process_data = interval(image_data)
    plt.subplot(projection = wcs, slices=('x','y',0,0))
    imshow_norm(process_data, cmap='gray')

    
    plt.savefig(fitsfile[:-5] + ".png")

def make_html(RATar, DECTar,  lotss_result_file, 
              extreme_catalogue, result, targRA, 
              targDEC,nchan, av_time):
    
    from dash import Dash, dcc, html, Input, Output, no_update
    import plotly.graph_objects as go
#    import pandas as pd
    from dash.exceptions import PreventUpdate

    # Small molcule drugbank dataset
    # Source: https://raw.githubusercontent.com/plotly/dash-sample-apps/main/apps/dash-drug-discovery/data/small_molecule_drugbank.csv'
    data_path = 'delay_calibrators.csv'

    df = Table.read(data_path,format='csv')
#    df = pd.read_csv(data_path, header=0,)

    print(df)
    observation = df["Observation"]

    fig = go.Figure(data=[
        go.Scatter(
            x=df["RA"],
            y=df["DEC"],
            mode="markers",
            customdata = df['Observation'],
            marker=dict(
                colorscale='viridis',
                color=df["Total_flux"],
                size=df["Total_flux"],
                colorbar={"title": "Total<br>Flux[mJy]"},
                line={"color": "#444"},
                reversescale=True,
                sizeref=45,
                sizemode="diameter",
                opacity=0.8,
                
            ),
        )
    ])

    # turn off native plotly.js hover effects - make sure to use
    # hoverinfo="none" rather than "skip" which also halts events.
    #fig.update_traces(hoverinfo="none", hovertemplate=None)

    fig.update_layout(
        xaxis=dict(title='RA'),
        yaxis=dict(title='DEC'),
        plot_bgcolor='rgba(255,255,255,0.1)'
    )

    # Update layout and update traces
    fig.update_layout(clickmode='event+select')
    fig.update_traces(marker_size=20)

    app = Dash(__name__)


    # Create app layout to show dash graph
    app.layout = html.Div(
    [
        dcc.Graph(
            id="graph_interaction",
            figure=fig,
        ),
        html.Img(id='image', src='')
    ]
    )

    # html callback function to hover the data on specific coordinates
    @app.callback(
    Output('image', 'src'),
    Input('graph_interaction', 'hoverData'))
    def open_url(hoverData):
        if hoverData:
            print(hoverData["points"][0])
            obs_id = hoverData["points"][0]["customdata"]
            img_str = "data:./" + obs_id + "_vlass.png;base64"
            return img_str
        else:
            raise PreventUpdate
        
    return app
    

def generate_catalogues( RATar, DECTar, targRA = 0.0, targDEC = 0.0, lotss_radius=1.5, lbcs_radius=1.5, im_radius=1.24,
               bright_limit_Jy=5., lotss_catalogue='lotss_catalogue.csv', lbcs_catalogue='lbcs_catalogue.csv', lotss_result_file='image_catalogue.csv',
               delay_cals_file='delay_calibrators.csv', match_tolerance=5., image_limit_Jy=0.01, continue_no_lotss=False,
                nchan = 16,  av_time = 1., vlass=False, html=False, outdir='.'):

#def plugin_main( RA, DEC, **kwargs ):
#    im_radius = float(kwargs['im_radius'])
#    image_limit_Jy = float(kwargs['image_limit_Jy'])
#    bright_limit_Jy = float(kwargs['bright_limit_Jy'])
#    match_tolerance = float(kwargs['match_tolerance'])
    #mslist = DataMap.load(mapfile_in)
    #MSname = mslist[0].file
    # For testing
    #MSname = kwargs['MSname']

    ## prepend everything with outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    lotss_catalogue = os.path.join(outdir,lotss_catalogue)
    lbcs_catalogue = os.path.join(outdir,lbcs_catalogue)
    lotss_result_file = os.path.join(outdir,lotss_result_file)
    delay_cals_file=os.path.join(outdir,delay_cals_file)

 
    ## first check for a valid delay_calibrator file
    if os.path.isfile(delay_cals_file):
        print( 'Delay calibrators file {:s} exists! returning.'.format(delay_cals_file) )
        return

    ## look for or download LBCS
    print("Attempting to find or download LBCS catalogue.")
    lbcs_catalogue = my_lbcs_catalogue( RATar, DECTar, Radius=lbcs_radius, outfile=lbcs_catalogue )
    ## look for or download LoTSS
    print("Attempting to find or download LoTSS catalogue.")
    lotss_catalogue = my_lotss_catalogue( RATar, DECTar,Radius=lotss_radius, bright_limit_Jy=bright_limit_Jy, faint_limit_Jy = 0.0,
                                         outfile=lotss_catalogue )
    
    print("Finding bright sources outside field")
    extreme_catalogue = my_lotss_catalogue( RATar, DECTar,Radius=10.0, bright_limit_Jy=1000., faint_limit_Jy = 10.0, outfile = os.path.join(outdir,"extreme_catalogue.csv") )
    extreme_catalogue = remove_multiples_position(extreme_catalogue)
    ## if lbcs exists, and either lotss exists or continue_without_lotss = True, continue
    ## else provide an error message and stop
    if len(lbcs_catalogue) == 0:
        logging.error('LBCS coverage does not exist, and catalogue not found on disk.')
        return
    if len(lotss_catalogue) == 0 and not fail_lotss_ok:
        logging.error('LoTSS coverage does not exist, and contine_without_lotss is set to False.')
        return 
    
    


        

    ## if the LoTSS catalogue is empty, write out the delay cals only
    if len(lotss_catalogue) == 0:
        print('Target field not in LoTSS coverage yet! Only writing {:s} based on LBCS'.format(delay_cals_file))

        ## Add the radius from phase centre to the catalogue
        #RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(MSname)[0])
        ptg_coords = SkyCoord( RATar, DECTar, frame='icrs', unit='deg' )

        src_coords = SkyCoord( lbcs_catalogue['RA'], lbcs_catalogue['DEC'], frame='icrs', unit='deg' )
        separations = src_coords.separation(ptg_coords )
        seps = Column( separations.deg, name='Radius' )
        lbcs_catalogue.add_column( seps )

        ## rename the source_id column
        #lbcs_catalogue.rename_column('Observation','Source_id')

        ## add in some dummy data
        Total_flux = Column( np.ones(len(lbcs_catalogue))*1e3, name='Total_flux', unit='mJy' )
        lbcs_catalogue.add_column( Total_flux )
        LGZ_Size = Column( np.ones( len(lbcs_catalogue) )*20., name='LGZ_Size', unit='arcsec' ) ## set to a default of 20 arcsec
        lbcs_catalogue.add_column( LGZ_Size )

        ## remove duplicate sources if necessary 
        lbcs_catalogue = remove_multiples_position( lbcs_catalogue )

        ## order based on radius from the phase centre
        lbcs_catalogue.sort('Radius')

        ## write the catalogue
        lbcs_catalogue.write(delay_cals_file, format='csv')

        result = lbcs_catalogue
    else:
        ## else continue 
        result = find_close_objs( lotss_catalogue, lbcs_catalogue, tolerance=match_tolerance )
        ## check if there are any matches
        if len(result) == 0:
            logging.error('LoTSS and LBCS coverage exists, but no matches found. This indicates something went wrong, please check your catalogues.')
        else:
        # add radius to the catalogue
            #RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(MSname)[0])
            #result = lbcs_catalogue
            ptg_coords = SkyCoord( RATar, DECTar, frame='icrs', unit='deg' )

            src_coords = SkyCoord( result['RA'], result['DEC'], frame='icrs', unit='deg' )
            separations = src_coords.separation(ptg_coords )
            seps = Column( separations.deg, name='Radius' )
            result.add_column( seps )

            ## order by radius from the phase centre
            result.sort( 'Radius' )

            #result.rename_column('Observation','Source_id')

            ## Write catalogues
            ## 1 - delay calibrators -- from lbcs_catalogue
            result.write( delay_cals_file, format='csv' )
            print('Writing delay calibrator candidate file {:s}'.format(delay_cals_file))

            ## sources to image -- first remove things that are already in the delay_cals_file
            good_index = [ x for x, src_id in enumerate( lotss_catalogue['Source_id'] ) if src_id not in result['Source_id'] ]
        
            tmp_cat = lotss_catalogue[good_index]

            ## make a flux cut
            image_index = np.where( tmp_cat['Peak_flux'] >= image_limit_Jy*1e3 )[0]
            flux_cut_sources = tmp_cat[image_index]

            ## make a radius cut
            src_coords = SkyCoord( flux_cut_sources['RA'], flux_cut_sources['DEC'], frame='icrs', unit='deg' )
            separations = src_coords.separation( ptg_coords )
            seps = Column( separations.deg, name='Radius' )
            flux_cut_sources.add_column( seps )
            good_idx = np.where( flux_cut_sources['Radius'] <= im_radius )[0]
            sources_to_image = flux_cut_sources[good_idx]

            nsrcs = float( len( sources_to_image ) )
            print( "There are "+str(len(lbcs_catalogue))+" delay calibrators within " + str(im_radius) + " degrees of the pointing centre")
            print( "There are "+str(nsrcs)+" sources above "+str(image_limit_Jy*1000)+" mJy within "+str(im_radius)+" degrees of the phase centre.")
            sources_to_image.write( lotss_result_file, format='csv' )

    print("Assumed averaging - nchannels: %s; time averaging: %s"%(nchan, av_time))
    make_plot(RATar, DECTar,  lotss_result_file, extreme_catalogue, result, targRA, targDEC,nchan = nchan, av_time = av_time, outdir=outdir)

    
    
    if vlass:
        from vlass_search import search_vlass
        ## Get cutouts of all LBCS sources
        print("Getting cutouts of LBCS sources")
        for i, source in enumerate(result):
            ra, dec = source['RA'], source['DEC']
            c = SkyCoord(ra, dec, unit = (u.deg, u.deg))
            outfile = os.path.join(outdir,"%s_vlass.fits"%source['Observation'])
            print('OUTFILE IS:')
            print(outfile)
            try:
                search_vlass(c, crop = True, crop_scale = 256)
                os.system("mv vlass_post**.fits  %s"%outfile)
                convert_vlass_fits(outfile)
            except:
                pass

    if html:
        app = make_html(RATar, DECTar,  lotss_result_file, extreme_catalogue, result, targRA, targDEC,nchan = nchan, av_time = av_time)
        app.run_server(debug=True, use_reloader=False)
    #return




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( '--output_dir', dest='outdir', type=str, help='directory to save results in [default cwd]', default='.' )
    parser.add_argument( '--lotss_radius', dest='lotss_radius', type=float, help='Radius to search LoTSS', default=1.5 )
    parser.add_argument( '--lbcs_radius', dest='lbcs_radius', type=float, help='Radius to search LBCS', default=1.5 )
    parser.add_argument( '--im_radius', dest='im_radius', type=float, help='Radius in which to image', default=1.24 )
    parser.add_argument( '--lotss_catalogue', dest='lotss_catalogue', type=str, help='input file for LoTSS catalogue [will be downloaded if does not exist]', default='lotss_catalogue.csv' )
    parser.add_argument( '--lbcs_catalogue', dest='lbcs_catalogue', type=str, help='input file for LBCS catalogue [will be downloaded if does not exist]', default='lbcs_catalogue.csv' )
    parser.add_argument( '--lotss_result_file', dest='lotss_result_file', type=str, help='output file of sources to image', default='image_catalogue.csv' )
    parser.add_argument( '--delay_cals_file', dest='delay_cals_file', type=str, help='output file of delay calibrators', default='delay_calibrators.csv' )
    parser.add_argument( '--match_tolerance', dest='match_tolerance', type=float, help='radius for matching LBCS to LoTSS [arcsec]', default=5. )
    parser.add_argument( '--bright_limit_Jy', dest='bright_limit_Jy', type=float, help='Flux limit for bright sources [Jy]', default=5.0 )
    parser.add_argument( '--image_limit_Jy', dest='image_limit_Jy', type=float, help='Flux limit for which sources to image [Jy]', default = 0.01 )
    parser.add_argument( '--continue_no_lotss', dest='continue_no_lotss', action='store_true', help='Continue with the pipeline if no lotss cross-matches can be found?', default = False )
    parser.add_argument( '--targRA', type=float, dest='targRA', help='Target RA in deg' )
    parser.add_argument( '--targDEC', type=float, dest='targDEC', help='Target DEC in deg' )
    parser.add_argument( '--nchan', dest='nchan', type=int, help='Number of frequency channels', default=16 )
    parser.add_argument( '--av_time', dest='av_time', type=float, help='Time averaging', default = 1. )

    parser.add_argument( '--MS', type=str, dest='MS', help='Measurement Set' )
    parser.add_argument( '--vlass', dest='vlass', action='store_true', help='Get VLASS cutouts of delay_calibraors.csv sources', default = False)
    parser.add_argument( '--html', dest='html', action='store_true', help='Create html page of plots', default = False)

    parser.add_argument( '--RA', type=float, dest = 'RA', help='Ptg RA in deg' )
    parser.add_argument( '--DEC', type=float, dest = 'DEC', help='Ptg DEC in deg' )


    args = parser.parse_args()

    if args.MS is not None:
        print("Using MS to get RA and DEC") 
        ptgRA, ptgDEC = grab_coo_MS(args.MS)
    else:
        ptgRA, ptgDEC = args.RA, args.DEC
        
    generate_catalogues( float(ptgRA), float(ptgDEC), targRA = args.targRA, targDEC = args.targDEC, 
         lotss_radius=args.lotss_radius, lbcs_radius=args.lbcs_radius, im_radius=args.im_radius,
           bright_limit_Jy=args.bright_limit_Jy, lotss_catalogue=args.lotss_catalogue, 
           lbcs_catalogue=args.lbcs_catalogue, lotss_result_file=args.lotss_result_file, 
           delay_cals_file=args.delay_cals_file, match_tolerance=args.match_tolerance, 
           image_limit_Jy=args.image_limit_Jy, continue_no_lotss = args.continue_no_lotss,
            nchan = args.nchan,  av_time = args.av_time, vlass=args.vlass, html = args.html, outdir=args.outdir)






### TO DO LIST

# Change size of targets based on flux - Add argument for setting minimum flux - DONE
# Look for bright sources within 10 degrees above a few Jy - use the bright_limit_Jy argument - DONE
# Add colour for quality of LBCS sources - Use fit or the compactness codes - Maybe not 
# Number the calibrators based on brightness - DONE


# Label shading of the smearing - DONE
# Print some statistics about the sources - Number etc. - DONE
# How many sources above the flux limit - DONE
# Target to phase center distance - DONE
# Make all the grey's lighter - DONE
# Print what averaging was assumed - Also add as an optional argument - Both channels and time - DONE
