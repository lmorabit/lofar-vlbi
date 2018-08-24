import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import numpy as np
import glob


# Leah Morabito, May 2017

def plugin_main(args, **kwargs):
    """
    Reads in closure phase file and returns the best delay calibrator mapfile

    Parameters
    ----------
    mapfile_dir : str
        Directory for output mapfile
    closurePhaseMap: str
        Name of output mapfile
    closurePhase_file: str
	Name of file with closure phase scatter

    Returns
    -------
    result : dict
        Output datamap closurePhaseFile

    """
    delaycal_list	= kwargs['delaycals']

    # read the file
    with open( 'closure_phases.txt', 'r' ) as f:
	lines = f.readlines()
    f.close()

    ## get lists of directions and scatter
    direction = []
    scatter = []
    for l in lines:
        direction.append(l.split()[4])
        scatter.append(np.float(l.split()[6]))

    ## convert to numpy arrays
    direction = np.asarray( direction )
    scatter = np.asarray( scatter )

    ## find the minimum scatter
    if len(scatter) > 1:
        min_scatter_index = np.where( scatter == np.min( scatter ) )[0]
        best_calibrator = direction[min_scatter_index[0]]
    else:
        best_calibrator = direction[0][0]

    with open( delaycal_list, 'r' ) as f:
	lines = f.readlines()
    f.close()

    for l in lines:
	tmp = l.split(',')
	if tmp[0] == best_calibrator:
	    cal_ra = tmp[4]
	    cal_dec = tmp[5]
	    cal_total_flux = tmp[7]

    ss = ','.join([cal_ra, cal_dec, best_calibrator, cal_total_flux])

    outfile = 'primary_delay_calibrator.csv'
    with open( outfile, 'w' ) as f:
	f.write(ss)
    f.close()

    result = {'calfile': outfile}  ## add coordinates here

    return result
    

