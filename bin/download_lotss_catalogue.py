#!/usr/bin/env python
import pyvo as vo
import os
import glob
import numpy as np
import pyrap.tables as pt

## written by Leah Morabito May 2017
## helper functions taken from download_tgss_skymodel_target.py in prefactor pipeline


########################################################################
def grab_coo_MS(MS):
    """
    Read the coordinates of a field from one MS corresponding to the selection given in the parameters

    Parameters
    ----------
    MS : str
        Full name (with path) to one MS of the field

    Returns
    -------
    RA, Dec : "tuple"
        coordinates of the field (RA, Dec in deg , J2000)
    """

    # reading the coordinates ("position") from the MS
    # NB: they are given in rad,rad (J2000) 
    [[[ra,dec]]] = pt.table(MS+'/FIELD', readonly=True, ack=False).getcol('PHASE_DIR')

    # RA is stocked in the MS in [-pi;pi]
    # => shift for the negative angles before the conversion to deg (so that RA in [0;2pi])
    if ra<0:
        ra=ra+2*np.pi

    # convert radians to degrees
    ra_deg =  ra/np.pi*180.
    dec_deg = dec/np.pi*180.

    # and sending the coordinates in deg
    return ra_deg,dec_deg

########################################################################
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

def main(ms_input, ResultsFile, Radius=1.5, DoDownload="True", AllFile=None):

    """
    Download the LoTSS skymodel for the target field

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    ResultsFile : str
        Full name (with path) to the skymodel; if YES is true, the LOTSS skymodel will be downloaded here
    Radius : float (default = 1.5)
        Radius for the LOTSS cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Download or not the LOTSS skymodel.
        "Force": download skymodel from LOTSS, delete existing skymodel if needed.
        "True" or "Yes": use existing skymodel file if it exists, download skymodel from 
                         LOTSS if it does not.
        "False" or "No": Do not download skymodel, raise an exception if skymodel
                         file does not exist.
    
    """

    FileExists = os.path.isfile(ResultsFile)
    if (not FileExists and os.path.exists(ResultsFile)):
        raise ValueError("download_lotss_skymodel: WTF! Path: \"%s\" exists but is not a file!"%(ResultsFile))
    download_flag = False
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(ResultsFile)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print "USING the exising skymodel in "+ ResultsFile
            return
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print "USING the exising skymodel in "+ ResultsFile
            return
         else:
            raise ValueError("download_lotss_skymodel: Path: \"%s\" does not exist and LOTSS download is disabled!"%(ResultsFile))

    # If we got here, then we are supposed to download the skymodel.
    assert download_flag == True # Jaja, belts and suspenders...
    print "DOWNLOADING LOTSS Skymodel for the target into "+ ResultsFile

    # Reading a MS to find the coordinate (pyrap)
    #[RATar,DECTar]=grab_coo_MS(input2strlist_nomapfile(ms_input)[0])
    #mypos = ( RATar, DECTar )
    RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    ## this is the tier 1 database to query
    url = 'http://vo.astron.nl/lofartier1/q/cone/scs.xml'

    ## this works
    query = vo.dal.scs.SCSQuery( url )
    query['RA'] = float( RATar )
    query['DEC'] = float( DECTar )
    query.radius = float( Radius )
    t = query.execute()   
    ## this does not
    #t = vo.conesearch( url, pos=mypos, radius=float(Radius) )

    ## convert to VO table
    tb = t.votable.to_table()
    ## find unresolved
    nsrcs = float( len( tb.columns['Resolved'] ) )
    unresolved_index = np.where( tb.columns['Resolved'] == 'U' )[0]
    perc_unres = len( unresolved_index ) / nsrcs * 100.
    print 'Percentage of sources which are unresolved: '+str( perc_unres )
    utb = tb[unresolved_index]
 
    ## sort by flux
    flux_sort = utb.argsort('Total_flux')
    utb_sorted = utb[ flux_sort[::-1] ]
    if AllFile is not None:
        utb_sorted.write( AllFile, format='ascii.csv' )
    
    ## keep only RA, DEC, Source_id and reorder the columns
    utb_sorted.keep_columns(['RA','DEC','Source_id'])

    utb_final = utb_sorted[ 'RA', 'DEC', 'Source_id' ]

    utb_final.write( ResultsFile, format='ascii.csv' )

    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=' Download the TGSS skymodel for the target field')

    parser.add_argument('MSfile', type=str, nargs='+', help='One (or more MSs) for which a LoTSS skymodel will be download.')
    parser.add_argument('--Radius', type=float, help='Radius for the LoTSS cone search in degrees')
    parser.add_argument('--Outfile', type=str, help='Filename to save the results to')

    args = parser.parse_args()
    radius=1.5
    if args.Radius:
        radius=args.Radius

    main(args.MSfile,args.Outfile,Radius=radius)

