import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import numpy as np
import pyrap.tables, math
from correlate import *
from astropy import units as u
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
    manual      = kwargs['manual']
    manual      = manual.lower() in ['true','t','1']
    data        = DataMap.load(infile_map)	# these are actual MS files
    datalist    = [data[i].file for i in xrange(len(data))]
    
    if manual:
        target_file = kwargs['target_file']
        pass
    try:
        nP = int(kwargs['nP'])	# if the user has defined a different value, use it
        pass
    except:
        nP = 3
        pass
    try:
        radius = float(kwargs['radius'])	# if the user has defined a different value, use it
        pass
    except:
        radius = 2.5
        pass

    ## if tick = 0, need to do the work to make directions files etc., otherwise just update the ticker
   
    fileid    = os.path.join(mapfile_dir, filename)	           # this file holds all the output measurement sets
    bigfileid = os.path.join(mapfile_dir, filename + '_bigfield')  # this big file holds all the directions
    
    if tick == 0:
        map_out_big = DataMap([])

        if manual:
            with open(target_file, 'r') as f:			# if the user has provided a list of targets, use it: otherwise use Lobos to find good targets
                for line in f:
		   if 'RA' not in line:
                      coords = (line.rstrip('\n')).split(',')
                      map_out_big.data.append(DataProduct( '[\"'+coords[0]+'\",\"'+coords[1]+'\"]' , coords[2], False ))

        else:
            infile = ((DataMap.load(infile_map))[0]).file	# get the actual filename from the map provided
            table = pyrap.tables.table(infile + '/FIELD', readonly = True)
            ra    = math.degrees(float(table.getcol('PHASE_DIR')[0][0][0] ) % (2 * math.pi))
            dec   = math.degrees(float(table.getcol('PHASE_DIR')[0][0][-1]))
            table.close()
            hexcoords = SkyCoord(ra, dec, unit = 'degree', frame='fk5')
            hexcoords = hexcoords.to_string('hmsdms', sep=':')
            hexra = [hexcoords.split(' ')[0] for i in hexcoords]
            hexdec = [hexcoords.split(' ')[1] for i in hexcoords]

            if not os.path.isfile (outdir+'/lobos_stats.sum'):
                os.system ('wget http://www.jb.man.ac.uk/~njj/lobos_stats.sum -P '+outdir)
            lobos = np.loadtxt(outdir+'/lobos_stats.sum',dtype='S')
            for l in lobos:
                newcoords = SkyCoord(l[1],l[2], unit=(u.hourangle, u.deg),frame='fk5')
                new = np.array([newcoords.ra.degree,newcoords.dec.degree])
                try:
                    lobos_coord = np.vstack((lobos_coord,new))
                except:
                    lobos_coord = np.copy(new)
            a = correlate(np.array([[ra,dec]]),0,1,lobos_coord,0,1,radius)
            for i in np.asarray(a[:,1],dtype='int'):
                if lobos[i][5].count('P')>=nP:
                    namera = lobos[i,1].replace(':','').split('.')[0]
                    namedec = lobos[i,2].replace(':','').split('.')[0]
                    dpppra = lobos[i,1].replace(':','h',1).replace(':','m',1)+'s' ## phase-shift RA
                    dpppdec = lobos[i,2].replace(':','d',1).replace(':','m',1)+'s'## phase-shift DEC
                    hemisphere = '-' if '-' in hexdec else '+'
                    outfile = namera+hemisphere+namedec             ## outfilename
                    map_out_big.data.append(DataProduct('[\"'+dpppra+'\",\"'+dpppdec+'\"]', outfile, False ))
                    
        map_out_big.save(bigfileid)	        # save all directions
        current_coords = map_out_big[0].host	# save current direction
        current_name = map_out_big[0].file      # save current filename
        n = len(map_out_big)
    else:
        data_big = DataMap.load(bigfileid)	# load all directions
        current_coords = data_big[tick].host	# save current direction
        current_name = data_big[tick].file
	n = len(data_big) 			# current progress

    map_out = DataMap([])
    
    for msID, ms_file in enumerate(datalist): 
	map_out.data.append(DataProduct( data[msID].host, '/'.join(data[msID].file.split('/')[:-1]) + '/' + current_name + '_' + data[msID].file.split('/')[-1], data[msID].skip)) 
	pass
      
    map_out.save(fileid)			# save all output measurement sets
    
    if (tick + 1) == n:     			# check how far the progress is
        do_break = True
    else:
        do_break = False
    #result = {'targetlist':bigfileid,'cords':current_coords,'cdir':current_name,'cdir_pattern':'*'+current_name+'*','ndir':int(n),'break':do_break}
    result = {'targetlist':bigfileid,'cords':current_coords,'ndir':int(n),'break':do_break,'mapfile':fileid}
    return result
