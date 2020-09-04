#!/usr/bin/python
import numpy as np
import os,glob,sys
import losoto; from losoto.h5parm import h5parm
import argparse
import casacore.tables as ct
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import angle_utilities
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy.table import Table
import bdsf

#  The difmap installation must be built with the modified corplt.c 
#     from  https://github.com/nealjackson/loop3_difmap
#  You will also need PGPLOT_FONT=/soft/pgplot/grfont.dat
#     and PGPLOT_DIR=/soft/pgplot  (or wherever pgplot lives)
#  these are both in the singularity image!!

def make_ant_table( station_list ):

    tied = {'ST001': np.array([3826557.5, 461029.06, 5064908],
                              dtype='float64')}

    core = {'CS001HBA0': np.array([3826896.235, 460979.455, 5064658.203],
                                  dtype='float64'),
            'CS001HBA1': np.array([3826979.384, 460897.597, 5064603.189],
                                  dtype='float64'),
            'CS002HBA0': np.array([3826600.961, 460953.402, 5064881.136],
                                  dtype='float64'),
            'CS002HBA1': np.array([3826565.594, 460958.110, 5064907.258],
                                  dtype='float64'),
            'CS003HBA0': np.array([3826471.348, 461000.138, 5064974.201],
                                  dtype='float64'),
            'CS003HBA1': np.array([3826517.812, 461035.258, 5064936.15],
                                  dtype='float64'),
            'CS004HBA0': np.array([3826585.626, 460865.844, 5064900.561],
                                  dtype='float64'),
            'CS004HBA1': np.array([3826579.486, 460917.48, 5064900.502],
                                  dtype='float64'),
            'CS005HBA0': np.array([3826701.16, 460989.25, 5064802.685],
                                  dtype='float64'),
            'CS005HBA1': np.array([3826631.194, 461021.815, 5064852.259],
                                  dtype='float64'),
            'CS006HBA0': np.array([3826653.783, 461136.440, 5064824.943],
                                  dtype='float64'),
            'CS006HBA1': np.array([3826612.499, 461080.298, 5064861.006],
                                  dtype='float64'),
            'CS007HBA0': np.array([3826478.715, 461083.720, 5064961.117],
                                  dtype='float64'),
            'CS007HBA1': np.array([3826538.021, 461169.731, 5064908.827],
                                  dtype='float64'),
            'CS011HBA0': np.array([3826637.421, 461227.345, 5064829.134],
                                  dtype='float64'),
            'CS011HBA1': np.array([3826648.961, 461354.241, 5064809.003],
                                  dtype='float64'),
            'CS013HBA0': np.array([3826318.954, 460856.125, 5065101.85],
                                  dtype='float64'),
            'CS013HBA1': np.array([3826402.103, 460774.267, 5065046.836],
                                  dtype='float64'),
            'CS017HBA0': np.array([3826405.095, 461507.460, 5064978.083],
                                  dtype='float64'),
            'CS017HBA1': np.array([3826499.783, 461552.498, 5064902.938],
                                  dtype='float64'),
            'CS021HBA0': np.array([3826463.502, 460533.094, 5065022.614],
                                  dtype='float64'),
            'CS021HBA1': np.array([3826368.813, 460488.057, 5065097.759],
                                  dtype='float64'),
            'CS024HBA0': np.array([3827218.193, 461403.898, 5064378.79],
                                  dtype='float64'),
            'CS024HBA1': np.array([3827123.504, 461358.861, 5064453.935],
                                  dtype='float64'),
            'CS026HBA0': np.array([3826418.227, 461805.837, 5064941.199],
                                  dtype='float64'),
            'CS026HBA1': np.array([3826335.078, 461887.696, 5064996.213],
                                  dtype='float64'),
            'CS028HBA0': np.array([3825573.134, 461324.607, 5065619.039],
                                  dtype='float64'),
            'CS028HBA1': np.array([3825656.283, 461242.749, 5065564.025],
                                  dtype='float64'),
            'CS030HBA0': np.array([3826041.577, 460323.374, 5065357.614],
                                  dtype='float64'),
            'CS030HBA1': np.array([3825958.428, 460405.233, 5065412.628],
                                  dtype='float64'),
            'CS031HBA0': np.array([3826383.037, 460279.343, 5065105.85],
                                  dtype='float64'),
            'CS031HBA1': np.array([3826477.725, 460324.381, 5065030.705],
                                  dtype='float64'),
            'CS032HBA0': np.array([3826864.262, 460451.924, 5064730.006],
                                  dtype='float64'),
            'CS032HBA1': np.array([3826947.411, 460370.066, 5064674.992],
                                  dtype='float64'),
            'CS101HBA0': np.array([3825899.977, 461698.906, 5065339.205],
                                  dtype='float64'),
            'CS101HBA1': np.array([3825805.288, 461653.869, 5065414.35],
                                  dtype='float64'),
            'CS103HBA0': np.array([3826331.59, 462759.074, 5064919.62],
                                  dtype='float64'),
            'CS103HBA1': np.array([3826248.441, 462840.933, 5064974.634],
                                  dtype='float64'),
            'CS201HBA0': np.array([3826679.281, 461855.243, 5064741.38],
                                  dtype='float64'),
            'CS201HBA1': np.array([3826690.821, 461982.139, 5064721.249],
                                  dtype='float64'),
            'CS301HBA0': np.array([3827442.564, 461050.814, 5064242.391],
                                  dtype='float64'),
            'CS301HBA1': np.array([3827431.025, 460923.919, 5064262.521],
                                  dtype='float64'),
            'CS302HBA0': np.array([3827973.226, 459728.624, 5063975.3],
                                  dtype='float64'),
            'CS302HBA1': np.array([3827890.077, 459810.483, 5064030.313],
                                  dtype='float64'),
            'CS401HBA0': np.array([3826795.752, 460158.894, 5064808.929],
                                  dtype='float64'),
            'CS401HBA1': np.array([3826784.211, 460031.993, 5064829.062],
                                  dtype='float64'),
            'CS501HBA0': np.array([3825568.82, 460647.62, 5065683.028],
                                  dtype='float64'),
            'CS501HBA1': np.array([3825663.508, 460692.658, 5065607.883],
                                  dtype='float64')}

    antenna_soltab = {'RS106HBA': np.array([3829205.598, 469142.533000,
                                            5062181.002], dtype='float64'),
                      'RS205HBA': np.array([3831479.67, 463487.529000,
                                            5060989.903], dtype='float64'),
                      'RS208HBA': np.array([3847753.31, 466962.809000,
                                            5048397.244], dtype='float64'),
                      'RS210HBA': np.array([3877827.56186, 467536.604956,
                                            5025445.584], dtype='float64'),
                      'RS305HBA': np.array([3828732.721, 454692.421000,
                                            5063850.334], dtype='float64'),
                      'RS306HBA': np.array([3829771.249, 452761.702000,
                                            5063243.181], dtype='float64'),
                      'RS307HBA': np.array([3837964.52, 449627.261000,
                                            5057357.585], dtype='float64'),
                      'RS310HBA': np.array([3845376.29, 413616.564000,
                                            5054796.341], dtype='float64'),
                      'RS404HBA': np.array([0.0, 0.0, 0.0],
                                           dtype='float64'),  # not operational
                      'RS406HBA': np.array([3818424.939, 452020.269000,
                                            5071817.644], dtype='float64'),
                      'RS407HBA': np.array([3811649.455, 453459.894000,
                                            5076728.952], dtype='float64'),
                      'RS409HBA': np.array([3824812.621, 426130.330000,
                                            5069251.754], dtype='float64'),
                      'RS410HBA': np.array([0.0, 0.0, 0.0],
                                           dtype='float64'),  # not operational
                      'RS503HBA': np.array([3824138.566, 459476.972,
                                            5066858.578], dtype='float64'),
                      'RS508HBA': np.array([3797136.484, 463114.447,
                                            5086651.286], dtype='float64'),
                      'RS509HBA': np.array([3783537.525, 450130.064,
                                            5097866.146], dtype='float64'),
                      'DE601HBA': np.array([4034101.522, 487012.757,
                                            4900230.499], dtype='float64'),
                      'DE602HBA': np.array([4152568.006, 828789.153,
                                            4754362.203], dtype='float64'),
                      'DE603HBA': np.array([3940295.706, 816722.865,
                                            4932394.416], dtype='float64'),
                      'DE604HBA': np.array([3796379.823, 877614.13,
                                            5032712.528], dtype='float64'),
                      'DE605HBA': np.array([4005681.02, 450968.643,
                                            4926458.211], dtype='float64'),
                      'FR606HBA': np.array([4324016.708, 165545.525,
                                            4670271.363], dtype='float64'),
                      'SE607HBA': np.array([3370271.657, 712125.881,
                                            5349991.165], dtype='float64'),
                      'UK608HBA': np.array([4008461.941, -100376.609,
                                            4943716.874], dtype='float64'),
                      'DE609HBA': np.array([3727217.673, 655109.175,
                                            5117003.123], dtype='float64'),
                      'PL610HBA': np.array([3738462.416, 1148244.316,
                                            5021710.658], dtype='float64'),
                      'PL611HBA': np.array([3850980.881, 1438994.879,
                                            4860498.993], dtype='float64'),
                      'PL612HBA': np.array([3551481.817, 1334203.573,
                                            5110157.41], dtype='float64'),
                      'IE613HBA': np.array([3801692.0, -528983.94,
                                            5076958.0], dtype='float64')}


    for a in station_list:
        if a[:2] == 'ST':
            antenna_soltab.update(tied)  # there will only be the tied station
        if a[:2] == 'CS':
            antenna_soltab.update(core)
            break  # only add the core stations to the antenna table once

    keys_to_remove = []
    for key in antenna_soltab:
        if key not in station_list:
            keys_to_remove.append(key)

    for k in keys_to_remove:
        antenna_soltab.pop(k, None)

    return( antenna_soltab )


def insert_into_filestem( infile, insert ):
    tmp = infile.split('/')
    mypath = '/'.join(tmp[0:-1])
    myfile = tmp[-1]
    new_file = os.path.join( mypath, insert+myfile )
    return new_file

# given an array of OK channels, write difmap select line
def chan2write (output, a):
    nums=sorted(set(a))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    if len(gaps) > 19:
        ## standard difmap can't handle it
        ## print('There are more than 19 chunks of channels, exiting.')
        ## return(1)
        print('WARNING: There are more than 19 chunks of channels, difmap will only read in the first 19!!!')
        tmp = [ x for i,x in enumerate(gaps) if i < 19 ]
        gaps = tmp
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    seq = np.ravel(list(zip(edges,edges)))+1 # difmap first chan=1 not 0
    chan_file = insert_into_filestem( output, 'dchan_' )
    fo=open(chan_file,'w')
    fo.write('select I,')
    for i in range(len(seq)):
        fo.write('%d'%seq[i])
        if i<len(seq)-1:
            fo.write(',')
    fo.write('\n')
    fo.close()

def fake_telescope( msfile, revert=False ):
    obs_table = ct.table( msfile + '::OBSERVATION', readonly=False )
    if revert:
        obs_table.putcell('TELESCOPE_NAME',0, 'LOFAR')
    else:
        obs_table.putcell('TELESCOPE_NAME',0, 'EVN')
    obs_table.close()

def find_good_chans( fitsfile ):
    f = fits.open( fitsfile, readonly=False )
    ## get data
    mydata = f[0].data
    myvis = np.squeeze( mydata['DATA'] )
    ## structure is: [time,bl,pol,X] where X = [real, imag, weight]
    real = myvis[:,:,0,0] * myvis[:,:,3,0]
    ## ms2uvfits converts nans to zeros, and flagged data to zeros (?) so search for zeros
    goodchans = [ x for x in np.arange(real.shape[1]) if np.sum(real[:,x]) != 0. ]
    return goodchans

# given MS file, return list of good channels
def find_chan_ms (filename,datacolumn="CORRECTED_DATA",flagname="FLAG"):
    # open ms
    ms=ct.table(filename,ack=False,readonly=False)
    # read data, flags
    d=ms.getcol(datacolumn)
    f=ms.getcol(flagname)
    d[f==True]=0.
    # see presence of nonzero vis per chan
    goodchans=[]
    for i in range(d.shape[1]):
        datasum=np.sum(np.abs(d[:,i],dtype=np.float64))
        if datasum!=0.:
            goodchans.append(i)
    ms.putcol(datacolumn,d)
    ms.flush()
    ms.close()
    return goodchans

# convert corplt output to an amp/phase array 
def corplt2array(corpltfile):
    f = open(corpltfile)
    fc = f.readlines()
    f.close()
    ut=[]
    ## get a list of stations and indices
    ss = "grep 'Stn' {:s} > ant_list.txt".format( corpltfile )
    os.system( ss )
    with open( 'ant_list.txt', 'r' ) as f:
        lines = []
        for l in f:
            lines.append(l.rstrip('\n'))
    ant_idx = []
    ant_name = []
    tmp_name = 'fred'
    for l in lines:
        ant_idx.append(l.split()[1])
        tmp = l.split()[3]
        if tmp[0:2] == 'CS':
            if tmp == tmp_name:
                ant_name.append(tmp+'1')
            else:
                ant_name.append(tmp+'0')
                tmp_name = tmp
        else:
            ant_name.append(l.split()[3])
    stn = np.array([],dtype='S')
    for l in fc:
        ls = l.split()
        if not len(ls):
            continue
        if ls[0]=='Stn':
            ## reference antenna list to get the right name
            tmp_idx = ls[1]
            tmp_stn = [ y for x,y in zip(ant_idx,ant_name) if x == tmp_idx ][0]
            stn = np.append(stn,tmp_stn)
        if len(ls)==4 and ls[0] in ['Amp','Phs']:
            ut.append(float(ls[1]))  # faster than numpy append    
    ut = np.unique(np.asarray(ut,dtype=np.float64))
    amp = np.ones((len(stn),len(ut)),dtype=np.float64)*np.nan
    phs = np.ones((len(stn),len(ut)),dtype=np.float64)*np.nan
    amperr = np.ones((len(stn),len(ut)),dtype=np.float64)*np.nan
    phserr = np.ones((len(stn),len(ut)),dtype=np.float64)*np.nan
    stn_idx = 0
    for l in fc:
        if l[:3] in ['Amp','Phs']:
            t_ut,t_val,t_err = np.asarray(l.split()[1:],dtype=np.float64)
            ut_idx = np.argwhere(ut==t_ut)[0][0]
            if l[:3] == 'Amp':
                amp[stn_idx,ut_idx] = t_val
                amperr[stn_idx,ut_idx] = t_err
            else:
                phs[stn_idx,ut_idx] = t_val
                phserr[stn_idx,ut_idx] = t_err
        if l[:3] == 'Stn':
            stn_idx += 1
    return amp,amperr,phs,phserr,ut,stn

def clean_selfcal_loop (fs,weight):
    fs.write('uvw %s\n'%weight)
    fs.write('if(peak(flux,max)/imstat(rms) > clean_sigma)\n')
    fs.write('	repeat;\\\n')
    fs.write('		peakwin 1.5; clean; selfcal\n')
    fs.write('	until(peak(flux,max)/imstat(rms) < clean_sigma)\n')
    fs.write('end if\n')

# Process the input file supplied, converting to FITS if needed, and write
# a dchan_XXX file containing the bad channels.
def file2fits(infile,datacolumn):
    isms = os.path.isdir( infile )
    # Process file itself. If FITS, add '.fits' to the filename and any dchan_
    # file if not present already. Otherwise convert to FITS, replacing any '.ms'
    # or '.MS' extension with '.fits'. Do the same for any associated dchan_ file
    if not isms:
        print( 'File is already in fits format.' )
        if infile[-5:]!='.fits':
            os.system ('mv %s %s.fits'%(infile,infile))
            fitsfile=infile+'.fits'
        else:
            fitsfile=infile
    else:
        print( 'File is in measurement set format, checking to see if it already exists in fits format.' )	
        if infile[-3:] in ['.MS','.ms']:
            fitsfile = infile[:-3]+'.fits' 
        else:
            fitsfile = infile+'.fits'
        if not os.path.isfile( fitsfile ):
            ## check length of fitsfile name
            nchar = len( fitsfile.split('/')[-1] )
	    if nchar > 80:
		tmp = fitsfile.split('/')
		modname = tmp[-1].split('_')[0] + '.fits'
		tmp[-1] = modname
		fitsfile = '/'.join(tmp)
	    ## fake the telescope name so ASTRON-specific format (which converts some real data to zeros) isn't written
            #fake_telescope( infile, revert=False )
            os.system('ms2uvfits in=%s out=%s writesyscal=F'%(infile,fitsfile))
            ## change telescope name back
            #fake_telescope( infile, revert=True )
            ## change telescope name in uvfits file
            

    ## find the good channels (i.e., not empty) from the fits file
    chan_file = insert_into_filestem( fitsfile, 'dchan_' )
    if not os.path.isfile( chan_file ):
        a = find_good_chans( fitsfile )
        chan2write(fitsfile,a)
    return fitsfile

# Write and execute the difmap script
def dif_script (infile,pol='I',aipsno=340,clean_sigma=6,map_size=512,\
                pixel_size=100,obs_length=900,datacolumn='CORRECTED_DATA',\
                startmod=True,do_natural_weighting=True):

    fitsfile = file2fits(infile,datacolumn)
    # Open script and declare variables
    fs = open('dif_script','w')
    fs.write('float clean_sigma; clean_sigma = %f\n'%float(clean_sigma))
    fs.write('float map_size; map_size = %f\n'%float(map_size))
    fs.write('float pixel_size; pixel_size = %f\n'%float(pixel_size))
    fs.write('float obs_length; obs_length = %f\n'%float(obs_length))
    fs.write('float sn_p;\nfloat sn_a\n')
    # Do the imaging with channels selected by find_difmap_chan
    fs.write('observe %s\n'%fitsfile)
    chan_file = insert_into_filestem( fitsfile, 'dchan_' )
    fsel = open (chan_file)
    fs.write('%s'%fsel.readline().replace('I',pol))
    fsel.close()
    fs.write('selftaper 0.5,0.025\n')
    fs.write('mapsize map_size,pixel_size\n')
    fs.write('startmod "",1\n' if startmod else 'clean\n')
    fs.write('peakwin 1.5\nselfcal false,false,0\n')
    clean_selfcal_loop (fs,'2,0')
    clean_selfcal_loop (fs,'2,-1')
    if do_natural_weighting:
        clean_selfcal_loop (fs,'0,-1')
    else:
        clean_selfcal_loop (fs,'2,-1')
    fs.write('gscale\n')
    fs.write('sn_p = peak(flux,max)/imstat(rms)\n')
    fs.write('repeat;\\\n')
    fs.write('    selfcal true,true,obs_length\n')
    fs.write('    sn_a = peak(flux,max)/imstat(rms)\n')
    fs.write('    if(peak(flux,max)/imstat(rms) > clean_sigma)\n')
    fs.write('        if(sn_a <= sn_p)\n')
    fs.write('            peakwin 1.5; clean; selfcal\n')
    fs.write('        obs_length=obs_length/2\n')
    fs.write('            sn_a = sn_p\n')
    fs.write('        end if\n')
    fs.write('        if(sn_a > sn_p)\n')
    fs.write('            clrmod false,true\n')
    fs.write('        obs_length=obs_length/2\n')
    fs.write('            sn_a = sn_p\n')
    fs.write('        end if\n')
    fs.write('    else\n')
    fs.write('        obs_length=obs_length/2\n')
    fs.write('    end if\n')
    fs.write('until(obs_length < 2)\n')
    fs.write('selfcal true,true,0\n')
    fs.write('delwin\n')
    fs.write('clean 1000,0.01\n')
    local_fits = fitsfile.split('/')[-1]
    fs.write('device %s.ps/vps\ncorplot\n'%local_fits.replace('.fits','_auto%s.p'%pol))
    fs.write('device /NULL\nmapl\ncmul=3*imstat(rms)\n')
    fs.write('device %s.ps/vps\nmapl clean,false\n' % local_fits.replace('.fits','_auto%s'%pol))
    fs.write('save %s\n' % local_fits.replace('.fits','_auto%s'%pol) )
    ## make a 6 arcsec ish image
    fs.write('selftaper 0\n')
    fs.write('uvrange 0,0.05\n')
    fs.write('uvtaper 0.1,0.03\n')
    fs.write('device %s_uvcut.ps/vps\n'%local_fits.replace('.fits','_auto%s'%pol))
    fs.write('mapsize 256,1000\n')
    fs.write('restore 6000\n')
    fs.write('mapl clean,false\n')
    fs.write('device /NULL\n' )
    fs.write('wmap %s_uvcut.fits\n'%local_fits.replace('.fits','_auto%s'%pol) )
    fs.write('quit\n')
    fs.close()
    os.system('difmap <dif_script')
    return fitsfile

def offset_by(lon, lat, posang, distance):
    """
    Point with the given offset from the given point.
    Parameters
    ----------
    lon, lat, posang, distance : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the starting point,
        position angle and distance to the final point.
        Quantities should be in angular units; floats in radians.
        Polar points at lat= +/-90 are treated as limit of +/-(90-epsilon) and same lon.
    Returns
    -------
    lon, lat : `~astropy.coordinates.Angle`
        The position of the final point.  If any of the angles are arrays,
        these will contain arrays following the appropriate `numpy` broadcasting rules.
        0 <= lon < 2pi.
    Notes
    -----
    """
    from astropy.coordinates.angles import Angle

    # Calculations are done using the spherical trigonometry sine and cosine rules
    # of the triangle A at North Pole,   B at starting point,   C at final point
    # with angles     A (change in lon), B (posang),            C (not used, but negative reciprocal posang)
    # with sides      a (distance),      b (final co-latitude), c (starting colatitude)
    # B, a, c are knowns; A and b are unknowns
    # https://en.wikipedia.org/wiki/Spherical_trigonometry

    cos_a = np.cos(distance)
    sin_a = np.sin(distance)
    cos_c = np.sin(lat)
    sin_c = np.cos(lat)
    cos_B = np.cos(posang)
    sin_B = np.sin(posang)

    # cosine rule: Know two sides: a,c and included angle: B; get unknown side b
    cos_b = cos_c * cos_a + sin_c * sin_a * cos_B
    # sin_b = np.sqrt(1 - cos_b**2)
    # sine rule and cosine rule for A (using both lets arctan2 pick quadrant).
    # multiplying both sin_A and cos_A by x=sin_b * sin_c prevents /0 errors
    # at poles.  Correct for the x=0 multiplication a few lines down.
    # sin_A/sin_a == sin_B/sin_b    # Sine rule
    xsin_A = sin_a * sin_B * sin_c
    # cos_a == cos_b * cos_c + sin_b * sin_c * cos_A  # cosine rule
    xcos_A = cos_a - cos_b * cos_c

    A = Angle(np.arctan2(xsin_A, xcos_A), u.radian)
    # Treat the poles as if they are infinitesimally far from pole but at given lon
    small_sin_c = sin_c < 1e-12
    if small_sin_c.any():
        # For south pole (cos_c = -1), A = posang; for North pole, A=180 deg - posang
        A_pole = (90*u.deg + cos_c*(90*u.deg-Angle(posang, u.radian))).to(u.rad)
        if A.shape:
            # broadcast to ensure the shape is like that of A, which is also
            # affected by the (possible) shapes of lat, posang, and distance.
            small_sin_c = np.broadcast_to(small_sin_c, A.shape)
            A[small_sin_c] = A_pole[small_sin_c]
        else:
            A = A_pole

    outlon = (Angle(lon, u.radian) + A).wrap_at(360.0*u.deg).to(u.deg)
    outlat = Angle(np.arcsin(cos_b), u.radian).to(u.deg)

    return outlon, outlat

def new_coords( ref_coords, pos_angle, sep ):

        slat = ref_coords.represent_as(UnitSphericalRepresentation).lat
        slon = ref_coords.represent_as(UnitSphericalRepresentation).lon

        newlon, newlat = offset_by(lon=slon, lat=slat, posang=pos_angle, distance=sep)
        return SkyCoord(newlon, newlat, frame=ref_coords.frame)

def read_mod( difmap_mod ):

    ## read in the file
    with open( difmap_mod, 'r' ) as f:
        lines = f.readlines()
    f.close()

    ## get the central coordinates
    coords = lines[0].split()
    ra_cen = ':'.join(coords[3:6]).rstrip(',')
    dec_cen = ':'.join(coords[7:10])
    cen_coords = SkyCoord( ra_cen, dec_cen, frame='icrs', unit=(u.hourangle, u.degree) )

    ## get the model values
    idx = [ xx for xx,val in enumerate(lines) if 'Tentative' in val ][0]
    vals = [ line.rstrip('\n') for line in lines[idx+3:len(lines)] ]
    flux = [ np.float(vv.split()[0]) for vv in vals ]
    radius = [ np.float(vv.split()[1]) for vv in vals ]
    theta = [ np.float(vv.split()[2]) for vv in vals ]

    ## convert radius and theta to RA, DEC
    dirs = [ new_coords( cen_coords, thet*u.deg, rad*1e-3*u.arcsec ) for thet,rad in zip( theta, radius ) ]
    new_dirs = [ str(dir.to_string('hmsdms')) for dir in dirs ]

    cc_dict = {}
    for xx in range(len(dirs)):
        mykey = 'cc' + str(xx)
        cc_dict[mykey] = [flux[xx], dirs[xx] ]

    return( cc_dict )

def find_centre( val_list ):

    ra_vals = []
    dec_vals = []
    for mykey in val_list.keys():
        ra_vals.append(val_list[mykey][1].ra.value)
        dec_vals.append(val_list[mykey][1].dec.value)
    cen_ra = np.mean( ra_vals )
    cen_dec = np.mean( dec_vals )
    return ( cen_ra, cen_dec )

def lof_coords( myskycoord ):
    tmp = str( myskycoord.to_string('hmsdms') )
    tmp = tmp.replace('s','')
    tmp_ra = tmp.split()[0]
    tmp_dec = tmp.split()[1]
    tmp_ra = tmp_ra.replace('h',':')
    tmp_ra = tmp_ra.replace('m',':')
    tmp_dec = tmp_dec.replace('d','.')
    tmp_dec = tmp_dec.replace('m','.')
    tmp_dec = tmp_dec.replace('+','')
    return( tmp_ra + ', ' + tmp_dec )

#  **** Top-level routine to process an input visibility file ***
# Specified file may be fits or MS
# If MS, will be inspected for good channels and converted to FITS for mapping
#    (any .MS or .ms extension converted to .fits, otherwise .fits added)
#    Corrected_data column assumed unless otherwise specified
# If FITS and we have AIPS/Parseltongue, then inspected for good channels
#    in AIPS and mapped, *If no AIPS and no good-channel file, exit with error*.
# We assume that XX and YY have the same number of telescopes and UTs, and also
#    the same bad/missing channels, and image separately for separate XX/YY corrs
def main( infile, clean_sig=6, map_size=512, pix_size=100, obs_length=900, datacolumn='CORRECTED_DATA', startmod=True, verbose=False, pols='I', catalogue=None, naturalwt=True ):

    # current working directory
    current_dir = os.getcwd()
    ## make a working directory, move the data, and chdir
    # get filestem to make unique name
    infile = infile.rstrip('/')
    tmp = infile.split('/')
    filestem = tmp[-1].split('_')[0] 
    # working directory
    work_dir = os.path.join( current_dir, 'difmap_{:s}'.format( filestem ) )
    os.mkdir( work_dir )
    os.system( 'cp -r {:s} {:s}'.format( infile, work_dir ) )
    os.chdir( work_dir )
    ## redefine infile
    infile = glob.glob( os.path.join( work_dir, tmp[-1] ) )[0]

    ## if a catalogue is specified, override the imaging parameters
    if catalogue is not None:
        print( 'Catalogue is specified, reading information to set imaging parameters.' )
        ## get ra and dec
        [[[ra,dec]]] = ct.table( infile + '::FIELD', readonly=True ).getcol('PHASE_DIR' )
        # => shift for the negative angles before the conversion to deg (so that RA in [0;2pi])
        if ra<0:
            ra=ra+2*np.pi
        # convert radians to degrees
        ra_deg =  ra/np.pi*180.
        dec_deg = dec/np.pi*180.
        tgt_coords = SkyCoord( ra_deg, dec_deg, unit='deg' )       

	t = Table.read( catalogue, format='csv' )
	## more flexible RA/DEC column naming
    	mycols = t.colnames
	ra_col = [ val for val in mycols if val == 'RA' ]
	de_col = [ val for val in mycols if val == 'DEC' ]
    	if len(ra_col) == 1:
            ra_col = ra_col[0]
            de_col = de_col[0]
        else:
            ## pick LoTSS position
            ra_col = [ val for val in mycols if val == 'RA_LOTSS' ][0]
            de_col = [ val for val in mycols if val == 'DEC_LOTSS' ][0]
        coords = SkyCoord( t[ra_col], t[de_col], unit='deg' )
        seps = coords.separation(tgt_coords).value
	src_idx = np.where( seps == np.min(seps) )[0]
	src_tmp = t[src_idx]
	## use LGZ_Size if available, otherwise use DC_Maj
	if 'LGZ_Size' in src_tmp.colnames:
	    ## in arcsec
	    size_asec = src_tmp['LGZ_Size'][0]
        else:
            size_asec = src_tmp['DC_Maj'][0] * 60. * 60.

	padding = 1.5
	possible_map_sizes = np.array([512,1024,2048,4096,8192])
        possible_map_asec = possible_map_sizes * float(pix_size) * 1e-3  ## convert to arcsec
	possible_idx = np.where( size_asec*padding <= possible_map_asec )[0] 
        if len( possible_idx ) >= 1:
            map_size = possible_map_sizes[np.min( possible_idx )]
            print( 'Estimated source size {:s}, making image with {:s}x{:s} pixels ({:s}x{:s} arcsec)'.format( str(size_asec), str(map_size), str(map_size), str(map_size*float(pix_size)*1e-3), str(map_size*float(pix_size)*1e-3) ) )
	else:
	    print( 'Image size exceeds {:s} arcseconds! Are you sure you want to make this image?'.format( str(np.max(possible_map_asec)) ) )
            print( 'Image size too large, aborting.' )
            return

        total_flux = src_tmp['Total_flux'][0]

    if pols == 'I':
        ## use stokes I only for self-cal and imaging (best option)
        fitsfile = dif_script( infile, pol=pols, clean_sigma=clean_sig, map_size=map_size, pixel_size=pix_size, obs_length=obs_length, datacolumn=datacolumn, startmod=startmod, do_natural_weighting=naturalwt)
        corpltfile = glob.glob( os.path.join( work_dir, 'CORPLT' ) )[0]
        ampI,amperrI,phsI,phserrI,utI,stnI = corplt2array(corpltfile)
        corpltout = insert_into_filestem( corpltfile.replace('CORPLT','_CORPLT_I'), filestem )
        os.system('mv {:s} {:s}'.format( corpltfile, corpltout ) )
    else:
        ## self-cal the polarisations separately
        ## create a difmap script
        fitsfile = dif_script(infile, pol='XX', clean_sigma=clean_sig, map_size=map_size, pixel_size=pix_size, obs_length=obs_length, datacolumn=datacolumn, startmod=startmod, do_natural_weighting=naturalwt)
        ## plot solutions and get values
        corpltfile = glob.glob( os.path.join( work_dir, 'CORPLT' ) )[0]
        ampXX,amperrXX,phsXX,phserrXX,utXX,stnXX = corplt2array(corpltfile)
        corpltout = insert_into_filestem( corpltfile.replace('CORPLT','_CORPLT_XX'), filestem )
        os.system('mv {:s} {:s}'.format( corpltfile, corpltout ) )
        ## write a difmap script for the YY polarisation and run
        fitsfile = dif_script(infile, pol='YY', clean_sigma=clean_sig, map_size=map_size, pixel_size=pix_size, obs_length=obs_length, datacolumn=datacolumn, startmod=startmod)
        ## plot solutions and get values
        corpltfile = glob.glob( os.path.join( work_dir, 'CORPLT' ) )[0]
        ampYY,amperrYY,phsYY,phserrYY,utYY,stnYY = corplt2array(corpltfile)
        corpltout = insert_into_filestem( corpltfile.replace('CORPLT','_CORPLT_YY'), filestem )
        os.system('mv {:s} {:s}'.format( corpltfile, corpltout ) )

    if pols == 'I':
        ut = utI
        stn = stnI
    else:
        ut = utXX
        stn = stnXX

    ## convert time axis to lofar times
    myms = ct.table( infile )
    lof_times = myms.getcol('TIME')
    myms.close()
    utvec = ut - np.min(ut)
    time_ax = utvec + np.min(lof_times)

    ## get frequency axis
    myspw = ct.table( infile + '::SPECTRAL_WINDOW' )
    freq_ax = np.squeeze( myspw.getcol('CHAN_FREQ') )
    myspw.close()

    if pols == 'I':
	## split the stokes I correction across XX and YY
        tmp_amp = np.rollaxis(np.dstack((np.sqrt(ampI/2.),np.sqrt(ampI/2.))),1,0)
        tmp_phs = np.rollaxis(np.dstack((phsI,phsI)),1,0)
    else:
        ## combine XX and YY information and reformat axes
        tmp_amp = np.rollaxis(np.dstack((ampXX,ampYY)),1,0)
        tmp_phs = np.rollaxis(np.dstack((phsXX,phsYY)),1,0)
    ## expand to fill frequency axis
    tmp_amp2 = np.expand_dims( tmp_amp, axis=1 )
    tmp_phs2 = np.expand_dims( tmp_phs, axis=1 )
    amp = np.repeat( tmp_amp2, len(freq_ax), axis=1 )
    phs = np.repeat( tmp_phs2, len(freq_ax), axis=1 )
    phs = phs * -1.0 ## difmap has a different convention

    ## get antenna information
    new_ants  = make_ant_table( stn )

    ## get pointing information
    ptg = ct.table( infile + '::FIELD' )
    ptg_dir = ptg.getcol('PHASE_DIR')[0]
    ptg.close()
    new_dir = {}
    new_dir[filestem] = ptg_dir[0]

    ## write solutions to an h5parm
    h5parmfile = os.path.join( work_dir, filestem + '_sols.h5' )
    out_h5 = h5parm(h5parmfile, readonly = False)
    out_solset = out_h5.makeSolset(solsetName='sol000')
    antenna_table = out_solset.obj._f_get_child('antenna')
    antenna_table.append(new_ants.items())
    out_solset.obj.source.append(new_dir.items())
    out_solset.makeSoltab('amplitude',axesNames=['time','freq','ant','pol'], axesVals=[time_ax,freq_ax,stn,['XX','YY']], vals=amp, weights=np.ones_like(amp))   
    out_solset.makeSoltab('phase',axesNames=['time','freq','ant','pol'], axesVals=[time_ax,freq_ax,stn,['XX','YY']], vals=phs, weights=np.ones_like(phs))   
    out_h5.close()

    ## apply the solutions
    with open( 'applysols.parset', 'w' ) as f:
	f.write('msin={:s}\n'.format(infile))
        f.write('msin.datacolumn=DATA\n')
        f.write('msout=.\n')
        f.write('msout.datacolumn=CORRECTED_DATA\n')
        f.write('numthreads=6\n')
        f.write('steps=[applycal]\n')
        f.write('applycal.type=applycal\n')
        f.write('applycal.parmdb={:s}\n'.format(h5parmfile))
        f.write('applycal.steps=[applyphs,applyamp]\n')
        f.write('applycal.applyphs.correction=phase000\n')
        f.write('applycal.applyamp.correction=amplitude000\n')
    f.close()

    ss = 'NDPPP applysols.parset > applysols.log 2>&1'
    os.system( ss )

    ## rename files so they have unique names
    diflogs = glob.glob( os.path.join( work_dir, 'dif*' ) )
    for diflog in diflogs:
        diflogout = insert_into_filestem( diflog, filestem+'_' )
        os.system( 'mv {:s} {:s}'.format( diflog, diflogout ) )

    ## Make a BBS format skymodel from the difmap XX model
    ## get reference frequency
    rf = ct.taql('select REF_FREQUENCY from {:s}::SPECTRAL_WINDOW'.format(infile))
    ref_freq = rf.getcol('REF_FREQUENCY')[0]
    # find the model file
    modfiles = glob.glob( filestem + '*.mod' )
    if pols == 'I':
        xx_file = modfiles[0]
    else:
        xx_file = [mf for mf in modfiles if 'XX' in mf ][0]
    ## read the model
    xx_mod = read_mod( xx_file )
    ## find central point
    xx_cen = find_centre( xx_mod )
    ## and get the flux
    xx_flux = []
    for mykey in xx_mod.keys():
        xx_flux.append( xx_mod[mykey][0] )
    ## write the model file
    outfile = os.path.join( work_dir, filestem + '_StokesI.mod' )
    with open( outfile, 'w' ) as f:
        f.write( "# (Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='{:s}', SpectralIndex='[]') = format\n\n".format(str(ref_freq)))
        ## patch name
        patch_cen = lof_coords( SkyCoord( xx_cen[0], xx_cen[1], frame='icrs', unit='deg') )
        f.write( ", , difmap_cc, {:s}\n".format(patch_cen) )
        for mykey in xx_mod.keys():
            flux = xx_mod[mykey][0]
            coords = lof_coords( xx_mod[mykey][1] )
            f.write("{:s}, POINT, difmap_cc, {:s}, {:s}, 0.0, 0.0, 0.0, 0.00000e+00, 0.00000e+00, 0.00000e+00, {:s}, [-0.8]\n".format(mykey, coords, str(flux), str(ref_freq) ) )
    f.close()

    ## convert to a sourcedb
    ss = 'makesourcedb in={:s} out={:s} format="<"'.format( outfile, outfile.replace('mod','skymodel') )
    os.system( ss )

    

    ## run wsclean 
    #wsclean_name = filestem + '_wsclean' 
    ## convert the im params to wsclean format
    #ss = 'wsclean -j 16 -mem 20 -v -reorder -update-model-required -weight uniform -mf-weighting -weighting-rank-filter 3 -name {:s} -size {:s} {:s} -padding 1.4 -scale {:s}asec -channels-out 6 -data-column CORRECTED_DATA -niter 10000 -auto-threshold 3 -auto-mask 5 -mgain 0.8 -join-channels -fit-spectral-pol 3 -fit-beam {:s}'.format( wsclean_name, str(map_size), str(map_size), str(float(pix_size)*0.001), infile )
    #os.system( ss )

    ## TO DO: move final files
    ## h5parm, images, log files, and skymodel
    image_files = glob.glob( os.path.join( work_dir, '*.ps') )
    log_files = glob.glob( os.path.join( work_dir, '*log') )
    #wsclean_ims = glob.glob( os.path.join( work_dir, '*wsclean*MFS*fits' ) )
    myh5parm = glob.glob( os.path.join( work_dir, '*h5' ) )
    skymodel = glob.glob( os.path.join( work_dir, '*skymodel' ) )
    #file_list = image_files + log_files + wsclean_ims + myh5parm + skymodel
    file_list = image_files + log_files + myh5parm + skymodel
    for myfile in file_list:
        ff = myfile.split('/')[-1]
        ss = 'mv {:s} {:s}'.format( myfile, os.path.join( current_dir, ff ) )
        os.system( ss )
    ## move the infile back and rename it
    tmp = infile.split('/')[-1]
    new_file = tmp + '.selfcal'
    selfcal_file = os.path.join( current_dir, new_file )
    os.system( 'cp -r {:s} {:s}'.format( infile, selfcal_file ) )

    print( 'done' )
    
    

if __name__ == "__main__":

    parser=argparse.ArgumentParser("Self-calibration using difmap.")
    ## arguments / flags
    parser.add_argument("-v","--verbose",help="Be verbose, say everything program does. Default is False",required=False,action="store_true")
    parser.add_argument("--clean_sig",help="Clean sigma, default is 6.",required=False,default=6)
    parser.add_argument("--map_size",help="Number of pixels for map size, default is 512.",required=False,default=512)
    parser.add_argument("--pix_size",help="pixel size in units of mas, default 100",required=False,default=100)
    parser.add_argument("--obs_length",help="Observation length, default 900",required=False,default=900)
    parser.add_argument("--colname",type=str,help="Name of the data column. Default is CORRECTED_DATA.",required=False,default="CORRECTED_DATA")
    parser.add_argument("--startmod",help="Generate a starting model. Default is True",required=False,action="store_false")
    parser.add_argument("--pols",type=str,help="polarisations to self-calibrate. Default is Stokes I.",required=False,default="I")
    parser.add_argument("--catalogue",type=str,help="catalogue to help determine imaging parameters.",required=False,default=None)
    parser.add_argument("--naturalwt",type=bool,default=True)

    ## positionals
    parser.add_argument("filename",type=str,help="Name of the measurement set")

    args=parser.parse_args()

    main( args.filename, clean_sig=args.clean_sig, map_size=args.map_size, pix_size=args.pix_size, 
	obs_length=args.obs_length, datacolumn=args.colname, startmod=args.startmod, verbose=args.verbose, pols=args.pols, catalogue=args.catalogue, naturalwt=args.naturalw


