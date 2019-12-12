import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import numpy as np
import pyrap.tables, math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table


# Colm Coughlan, March 2016
# used in locapi generic pipeline implementation
# Alexander Drabent, January 2017
# added more functionality 
# Leah Morabito May 2017 updated to skip header line in manual target file if it exists


def plugin_main(args, **kwargs):
    """
    Takes in list of targets and returns the appropriate one in a mapfile
    Knows which is the current target by storing target ID in a mapfile
    Outputs an expanded list of the current target

    Parameters
    ----------
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    target_list: str
		List of all targets
	target_id: str
		Current 

    Returns
    -------
    result : dict
        Output datamap filename

    """
    infile_map  = kwargs['mapfile_in']   ## input measurement sets
    mapfile_dir = kwargs['mapfile_dir']  ## mapfile directory 
    filename    = kwargs['filename']   ## output measurement sets
    target_file = kwargs['target_file']  ## file with info on the calibrator
    data        = DataMap.load(infile_map)	# these are actual MS files
    datalist    = [data[i].file for i in xrange(len(data))] ## a list of the measurement sets 

    fileid    = os.path.join(mapfile_dir, filename)	           ## the output measurement sets
    coordfileid   = os.path.join(mapfile_dir, 'coords_'+filename )

    print( fileid )
    print( coordfileid )
    
    map_out = DataMap([])
    map_out_coords = DataMap([])

    ## load an astropy table
    t = Table.read( target_file, format='csv' )
    RA_val = t['RA_LOTSS'].data[0]
    DEC_val = t['DEC_LOTSS'].data[0]
    Source_id = t['Source_id'].data[0]

    ss = '[\"' + str(RA_val) + 'deg\",\"' + str(DEC_val) + 'deg\"]'

    map_out_coords.data.append(DataProduct( ss, Source_id, False ))

    map_out_coords.save(coordfileid)	        # save to file
    current_coords = map_out_coords[0].host	# save direction
    current_name = map_out_coords[0].file      # save filename

    map_out = DataMap([])
    
    for msID, ms_file in enumerate(datalist): 
	map_out.data.append(DataProduct( data[msID].host, '/'.join(data[msID].file.split('/')[:-1]) + '/' + current_name + '_' + data[msID].file.split('/')[-1], data[msID].skip)) 
	pass
      
    map_out.save(fileid)			# save all output measurement sets
    
    result = {'coordfile':coordfileid,'coords':current_coords,'mapfile':fileid}
    return result
