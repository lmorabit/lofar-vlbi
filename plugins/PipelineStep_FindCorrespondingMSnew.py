import os
import numpy as np
import pyrap.tables
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

# Colm Coughlan, March 2016
# Alexander Drabent, January 2017
# used in locapi generic pipeline implementation

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx]):
        return array[idx-1]
    else:
        return array[idx]


def plugin_main(args, **kwargs):
    """
    Find the measurement set closest to a given solution table, suitable for reading station names.

    Parameters
    ----------
    mapfile_ms : str
        Mapfile of the measurement sets
    mapfile_grpd : str
        Mapfile of the (grouped) calibration tables
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """

    mapfile_dir = kwargs['mapfile_dir']
    mapfile_ms = kwargs['mapfile_ms']
    mapfile_grpd = kwargs['mapfile_grpd']
    filename = kwargs['filename']

    result = {}
    data = DataMap.load(mapfile_ms)	# these are actual MS files
    groups = DataMap.load(mapfile_grpd)	# these are probably parmdbs

    datalist = [data[i].file for i in xrange(len(data))]
    grp_list = [groups[i].file for i in xrange(len(groups))]

    frequency_groups = []
    map_out = DataMap([]) 		## list of MS for which there are phase solutions (length valid MS)
    map_out_MStables = DataMap([])  	## list of phase tables to apply (length valid MS)
    map_out_addIS = DataMap([])  	## list of phase tables to apply (length phase tables)
	
    for grp_file in grp_list:
	table = pyrap.tables.table(grp_file, readonly = True)
	frequency_range = list(np.zeros(2))
	frequency_range[0] = float(table.getcol('STARTX')[0])
	frequency_range[1] = float(table.getcol('ENDX')[0])
	frequency_groups.append(frequency_range)
	table.close()
	pass

    tomap_addIS = 0
    for msID, ms_file in enumerate(datalist):
	table = pyrap.tables.table(ms_file + '/SPECTRAL_WINDOW', readonly = True)
	ref_frequency = float(table.getcol('REF_FREQUENCY')[0])
	table.close()
	for groupID,freq_group in enumerate(frequency_groups):
		if freq_group[0] <= ref_frequency < freq_group[1]:
			map_out.data.append( DataProduct( data[msID].host, data[msID].file, data[msID].skip ) )
			map_out_MStables.data.append( DataProduct(groups[groupID].host, groups[groupID].file, (groups[groupID].skip or data[msID].skip)))
			if tomap_addIS <= groupID:
				map_out_addIS.data.append( DataProduct(groups[groupID].host, groups[groupID].file, (groups[groupID].skip or data[msID].skip)))
				tomap_addIS += 1
				pass
			break
			pass
		else:
			continue
			pass
		pass
    	pass


    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result['mapfile'] = fileid

    fileid = os.path.join(mapfile_dir, filename + '_parmdbs')
    map_out_Mstables.save(fileid)
    result['parmdbs'] = fileid

    fileid = os.path.join(mapfile_dir, filename + '_ISparmdbs')
    map_out_addIS.save(fileid) 
    result['ISparmdbs'] = fileid

    return result
    pass
