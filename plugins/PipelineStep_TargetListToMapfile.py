import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import numpy as np
import pyrap.tables, math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

def plugin_main(args, **kwargs):
    """
    Takes in a catalogue with a target and returns an appropriate mapfile
    
    Parameters
    ----------
    mapfile_in: str
        
    mapfile_dir: str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    target_file: str
        file containing target info

    Returns
    -------
    result : dict
        Output datamap filename
    
    """
    # parse the inputs
    infile_map = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    target_file = kwargs['target_file']
    all_to_one = kwargs['all_to_one'].lower().capitalize()
    
    # the input data
    data = DataMap.load(infile_map)
    datalist = [ data[i].file for i in xrange(len(data)) ]
    
    # outfile information
    fileid = os.path.join(mapfile_dir, filename)
    coordfileid = os.path.join(mapfile_dir, 'coords_' + filename)

    # initialise the output data map for the coordinates
    map_out_coords = DataMap([])
    # read in the catalogue to get source_id, RA, and DEC
    t = Table.read(target_file, format='csv')
    RA_val = t['RA_LOTSS'].data[0]
    DEC_val = t['DEC_LOTSS'].data[0]
    Source_id = t['Source_id'].data[0]
    if str(Source_id)[0:1] == 'I':
	pass
    elif str(Source_id)[0:1] == 'S':
	pass
    else:
	Source_id = 'S' + str(Source_id)
    # make a string of coordinates for the NDPPP command
    ss = '["' + str(RA_val) + 'deg","' + str(DEC_val) + 'deg"]'
    # save the coordinate information
    map_out_coords.data.append(DataProduct(ss, Source_id, False))
    map_out_coords.save(coordfileid)
    # save the coords to a variable to return
    current_coords = map_out_coords[0].host

    # get the name (source_id)
    current_name = map_out_coords[0].file
    # initialise an output data map
    map_out = DataMap([])
    if all_to_one == 'True':
        msID = 0 
        ms_file = datalist[0]
        map_out.data.append(DataProduct(data[msID].host, '/'.join(data[msID].file.split('/')[:-1]) + '/' + current_name + '_' + data[msID].file.split('/')[-1], data[msID].skip))
    else:
        print( 'HELLO HELLO HELLO' )
        for msID, ms_file in enumerate(datalist):
	    map_out.data.append(DataProduct(data[msID].host, '/'.join(data[msID].file.split('/')[:-1]) + '/' + current_name + '_' + data[msID].file.split('/')[-1], data[msID].skip))
    # save the file
    map_out.save(fileid)
    result = {'coordfile': coordfileid,
     'coords': current_coords,
     'name': current_name,
     'mapfile': fileid}
    return result
