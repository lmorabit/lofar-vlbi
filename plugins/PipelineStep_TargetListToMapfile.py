import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import numpy as np
import pyrap.tables, math
from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord


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
    infile_map  = kwargs['infile']
    mapfile_dir = kwargs['mapfile_dir']
    filename    = kwargs['filename']
    outdir      = kwargs['wd']
    tick        = int(kwargs['counter'])
    data        = DataMap.load(infile_map)	# these are actual MS files
    datalist    = [data[i].file for i in xrange(len(data))] ## use only 20 subbands
    target_file = kwargs['target_file']

    fileid    = os.path.join(mapfile_dir, filename)	           # this file holds all the output measurement sets
    infileid    = os.path.join(mapfile_dir, 'input_' + filename )  # this file holds all the input measurement sets
    bigfileid = os.path.join(mapfile_dir, filename + '_bigfield')  # this big file holds all the directions
    
    ## if tick = 0, need to do the work to make directions files etc., otherwise just update the ticker
    if tick == 0:
        map_out_big = DataMap([])

	## use astropy.io.ascii to load the file
	with open( target_file, 'r' ) as f:
	    lines = f.readlines()
	f.close()
	RA_vals = []
	DEC_vals = []
	source_id = []
	for l in lines:
	    tmp = l.split(',')
	    RA_vals.append(tmp[0])
	    DEC_vals.append(tmp[1])
	    source_id.append(tmp[2])
	## get the coordinates
	for x, src_id in enumerate(source_id):
	    ss = '[\"'+str(RA_vals[x])+'deg\",\"'+str(DEC_vals[x])+'deg\"]'
	    map_out_big.data.append(DataProduct( ss, src_id, False ))

        map_out_big.save(bigfileid)	        # save all directions
        current_coords = map_out_big[0].host	# save current direction
        current_name = map_out_big[0].file      # save current filename
        n = len(map_out_big)
    ## not the first iteration, load the directions
    else:
        data_big = DataMap.load(bigfileid)	# load all directions
        current_coords = data_big[tick].host	# save current direction
        current_name = data_big[tick].file
	n = len(data_big) 			# current progress

    map_out = DataMap([])
    map_out2 = DataMap([])
    
    for msID, ms_file in enumerate(datalist): 
	map_out.data.append(DataProduct( data[msID].host, '/'.join(data[msID].file.split('/')[:-1]) + '/' + current_name + '_' + data[msID].file.split('/')[-1], data[msID].skip)) 
        map_out2.data.append(DataProduct( data[msID].host, data[msID].file, data[msID].skip))
	pass
      
    map_out.save(fileid)			# save all output measurement sets
    map_out2.save(infileid)			# save all input measurement sets
    
    if (tick + 1) == n:     			# check how far the progress is
        do_break = True
    else:
        do_break = False
    result = {'targetlist':bigfileid,'coords':current_coords,'ndir':int(n),'break':do_break,'mapfile':fileid,'inmapfile':infileid}
    return result
