#!/usr/bin/env python
import pyvo as vo
import os
import glob
import numpy as np
import argparse
import pyrap.tables as pt

## based on script by Leah Morabito
## updated by Kaspars and Atvars to work for lbcs
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

def main(ms_input, ResultsFile, Radius=1.5, DoDownload="True"):

    """
    Download the LBCS calibrator list for the target field and find things with good P values

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    ResultsFile : str
        Full name (with path) to the skymodel; if YES is true, the LOTSS skymodel will be downloaded here
    Radius : float (default = 1.5)
        Radius for the cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Download or not the LBCS catalogue.
        "Force": download catalogue, delete existing if needed.
        "True" or "Yes": use existing file if it exists, download if it does not.
        "False" or "No": Do not download catalogue, raise an exception if file does not exist.
    
    """

    FileExists = os.path.isfile(ResultsFile)
    if (not FileExists and os.path.exists(ResultsFile)):
        raise ValueError("download_lbcs_catalogue: WTF! Path: \"%s\" exists but is not a file!"%(ResultsFile))
    download_flag = False
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(ResultsFile)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print "USING the exising catalogue in "+ ResultsFile
            return
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print "USING the exising catalogue in "+ ResultsFile
            return
         else:
            raise ValueError("download_lbcs_catalogue: Path: \"%s\" does not exist and LBCS download is disabled!"%(ResultsFile))


    # If we got here, then we are supposed to download the skymodel.
    assert download_flag == True # Jaja, belts and suspenders...
    print "DOWNLOADING LBCS catalogue to "+ ResultsFile

    # Reading a MS to find the coordinate (pyrap)
    RATar, DECTar=grab_coo_MS(input2strlist_nomapfile(ms_input)[0])
    mypos = [( float(RATar), float(DECTar) )]
    #mypos = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    ## this is the tier 1 database to query
    url = 'http://vo.astron.nl/lbcs/lobos/cone/scs.xml'

    ## this works
    query = vo.dal.scs.SCSQuery( url )
    query['RA'] = float( RATar )
    query['DEC'] = float( DECTar )
    query.radius = float( Radius )
    t = query.execute()   
    ## this does not
    #t = vo.conesearch( url, pos=mypos, radius=float(Radius) )

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
    
    ############ Print the new table with 3 columns in it  ########################
    
    tb_out = tb_sorted['raj2000','decj2000','ObsID'] 
    tb_out.write( ResultsFile, format='ascii.csv')    
   
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=' Download LBCS catalogue for a target field and find good calibrators (more than 2 P values)')

    parser.add_argument('MSfile', type=str, nargs='+', help='One (or more MSs) for which an LBCS catalogue will be download.')
    parser.add_argument('--Radius', type=float, help='Radius for the LBCS cone search in degrees')
    parser.add_argument('--Outfile', type=str, help='Filename to save the results to')

    args = parser.parse_args()
    radius=1.5
    if args.Radius:
        radius=args.Radius

    main(args.MSfile,args.Outfile,Radius=radius)
