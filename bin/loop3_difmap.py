#!/usr/bin/python
import numpy as np,os,glob,sys
import losoto; from losoto.h5parm import h5parm
from pyrap.tables import table
import argparse
import casacore.tables as ct

#  The difmap installation must be built with the modified corplt.c 
#     from  https://github.com/nealjackson/loop3_difmap
#  You will also need PGPLOT_FONT=/soft/pgplot/grfont.dat
#     and PGPLOT_DIR=/soft/pgplot  (or wherever pgplot lives)

try:
    from AIPSData import AIPSUVData
    from AIPSTask import AIPSTask,AIPSList,AIPSMessageLog
    import Wizardry
    from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
    HAVE_AIPS = True
except:
    HAVE_AIPS = False

def find_first_unflagged_ant( msfile ):
    ant_t = ct.taql( 'select ANTENNA1 from {:s} where all(FLAG) limit 1'.format(msfile) )
    ant_idx = ant_t.getcol('ANTENNA1')[0]
    ant_table = ct.table( msfile + '::ANTENNA' )
    ant_names = ant_table.getcol('NAME')
    first_ant = ant_names[ant_idx]
    ant_table.close()
    return( first_ant )

def make_ant_table( station_list ):

    tied = {'ST001': np.array([3826557.5, 461029.06, 5064908],
                              dtype='float32')}

    core = {'CS001HBA0': np.array([3826896.235, 460979.455, 5064658.203],
                                  dtype='float32'),
            'CS001HBA1': np.array([3826979.384, 460897.597, 5064603.189],
                                  dtype='float32'),
            'CS002HBA0': np.array([3826600.961, 460953.402, 5064881.136],
                                  dtype='float32'),
            'CS002HBA1': np.array([3826565.594, 460958.110, 5064907.258],
                                  dtype='float32'),
            'CS003HBA0': np.array([3826471.348, 461000.138, 5064974.201],
                                  dtype='float32'),
            'CS003HBA1': np.array([3826517.812, 461035.258, 5064936.15],
                                  dtype='float32'),
            'CS004HBA0': np.array([3826585.626, 460865.844, 5064900.561],
                                  dtype='float32'),
            'CS004HBA1': np.array([3826579.486, 460917.48, 5064900.502],
                                  dtype='float32'),
            'CS005HBA0': np.array([3826701.16, 460989.25, 5064802.685],
                                  dtype='float32'),
            'CS005HBA1': np.array([3826631.194, 461021.815, 5064852.259],
                                  dtype='float32'),
            'CS006HBA0': np.array([3826653.783, 461136.440, 5064824.943],
                                  dtype='float32'),
            'CS006HBA1': np.array([3826612.499, 461080.298, 5064861.006],
                                  dtype='float32'),
            'CS007HBA0': np.array([3826478.715, 461083.720, 5064961.117],
                                  dtype='float32'),
            'CS007HBA1': np.array([3826538.021, 461169.731, 5064908.827],
                                  dtype='float32'),
            'CS011HBA0': np.array([3826637.421, 461227.345, 5064829.134],
                                  dtype='float32'),
            'CS011HBA1': np.array([3826648.961, 461354.241, 5064809.003],
                                  dtype='float32'),
            'CS013HBA0': np.array([3826318.954, 460856.125, 5065101.85],
                                  dtype='float32'),
            'CS013HBA1': np.array([3826402.103, 460774.267, 5065046.836],
                                  dtype='float32'),
            'CS017HBA0': np.array([3826405.095, 461507.460, 5064978.083],
                                  dtype='float32'),
            'CS017HBA1': np.array([3826499.783, 461552.498, 5064902.938],
                                  dtype='float32'),
            'CS021HBA0': np.array([3826463.502, 460533.094, 5065022.614],
                                  dtype='float32'),
            'CS021HBA1': np.array([3826368.813, 460488.057, 5065097.759],
                                  dtype='float32'),
            'CS024HBA0': np.array([3827218.193, 461403.898, 5064378.79],
                                  dtype='float32'),
            'CS024HBA1': np.array([3827123.504, 461358.861, 5064453.935],
                                  dtype='float32'),
            'CS026HBA0': np.array([3826418.227, 461805.837, 5064941.199],
                                  dtype='float32'),
            'CS026HBA1': np.array([3826335.078, 461887.696, 5064996.213],
                                  dtype='float32'),
            'CS028HBA0': np.array([3825573.134, 461324.607, 5065619.039],
                                  dtype='float32'),
            'CS028HBA1': np.array([3825656.283, 461242.749, 5065564.025],
                                  dtype='float32'),
            'CS030HBA0': np.array([3826041.577, 460323.374, 5065357.614],
                                  dtype='float32'),
            'CS030HBA1': np.array([3825958.428, 460405.233, 5065412.628],
                                  dtype='float32'),
            'CS031HBA0': np.array([3826383.037, 460279.343, 5065105.85],
                                  dtype='float32'),
            'CS031HBA1': np.array([3826477.725, 460324.381, 5065030.705],
                                  dtype='float32'),
            'CS032HBA0': np.array([3826864.262, 460451.924, 5064730.006],
                                  dtype='float32'),
            'CS032HBA1': np.array([3826947.411, 460370.066, 5064674.992],
                                  dtype='float32'),
            'CS101HBA0': np.array([3825899.977, 461698.906, 5065339.205],
                                  dtype='float32'),
            'CS101HBA1': np.array([3825805.288, 461653.869, 5065414.35],
                                  dtype='float32'),
            'CS103HBA0': np.array([3826331.59, 462759.074, 5064919.62],
                                  dtype='float32'),
            'CS103HBA1': np.array([3826248.441, 462840.933, 5064974.634],
                                  dtype='float32'),
            'CS201HBA0': np.array([3826679.281, 461855.243, 5064741.38],
                                  dtype='float32'),
            'CS201HBA1': np.array([3826690.821, 461982.139, 5064721.249],
                                  dtype='float32'),
            'CS301HBA0': np.array([3827442.564, 461050.814, 5064242.391],
                                  dtype='float32'),
            'CS301HBA1': np.array([3827431.025, 460923.919, 5064262.521],
                                  dtype='float32'),
            'CS302HBA0': np.array([3827973.226, 459728.624, 5063975.3],
                                  dtype='float32'),
            'CS302HBA1': np.array([3827890.077, 459810.483, 5064030.313],
                                  dtype='float32'),
            'CS401HBA0': np.array([3826795.752, 460158.894, 5064808.929],
                                  dtype='float32'),
            'CS401HBA1': np.array([3826784.211, 460031.993, 5064829.062],
                                  dtype='float32'),
            'CS501HBA0': np.array([3825568.82, 460647.62, 5065683.028],
                                  dtype='float32'),
            'CS501HBA1': np.array([3825663.508, 460692.658, 5065607.883],
                                  dtype='float32')}

    antenna_soltab = {'RS106HBA': np.array([3829205.598, 469142.533000,
                                            5062181.002], dtype='float32'),
                      'RS205HBA': np.array([3831479.67, 463487.529000,
                                            5060989.903], dtype='float32'),
                      'RS208HBA': np.array([3847753.31, 466962.809000,
                                            5048397.244], dtype='float32'),
                      'RS210HBA': np.array([3877827.56186, 467536.604956,
                                            5025445.584], dtype='float32'),
                      'RS305HBA': np.array([3828732.721, 454692.421000,
                                            5063850.334], dtype='float32'),
                      'RS306HBA': np.array([3829771.249, 452761.702000,
                                            5063243.181], dtype='float32'),
                      'RS307HBA': np.array([3837964.52, 449627.261000,
                                            5057357.585], dtype='float32'),
                      'RS310HBA': np.array([3845376.29, 413616.564000,
                                            5054796.341], dtype='float32'),
                      'RS404HBA': np.array([0.0, 0.0, 0.0],
                                           dtype='float32'),  # not operational
                      'RS406HBA': np.array([3818424.939, 452020.269000,
                                            5071817.644], dtype='float32'),
                      'RS407HBA': np.array([3811649.455, 453459.894000,
                                            5076728.952], dtype='float32'),
                      'RS409HBA': np.array([3824812.621, 426130.330000,
                                            5069251.754], dtype='float32'),
                      'RS410HBA': np.array([0.0, 0.0, 0.0],
                                           dtype='float32'),  # not operational
                      'RS503HBA': np.array([3824138.566, 459476.972,
                                            5066858.578], dtype='float32'),
                      'RS508HBA': np.array([3797136.484, 463114.447,
                                            5086651.286], dtype='float32'),
                      'RS509HBA': np.array([3783537.525, 450130.064,
                                            5097866.146], dtype='float32'),
                      'DE601HBA': np.array([4034101.522, 487012.757,
                                            4900230.499], dtype='float32'),
                      'DE602HBA': np.array([4152568.006, 828789.153,
                                            4754362.203], dtype='float32'),
                      'DE603HBA': np.array([3940295.706, 816722.865,
                                            4932394.416], dtype='float32'),
                      'DE604HBA': np.array([3796379.823, 877614.13,
                                            5032712.528], dtype='float32'),
                      'DE605HBA': np.array([4005681.02, 450968.643,
                                            4926458.211], dtype='float32'),
                      'FR606HBA': np.array([4324016.708, 165545.525,
                                            4670271.363], dtype='float32'),
                      'SE607HBA': np.array([3370271.657, 712125.881,
                                            5349991.165], dtype='float32'),
                      'UK608HBA': np.array([4008461.941, -100376.609,
                                            4943716.874], dtype='float32'),
                      'DE609HBA': np.array([3727217.673, 655109.175,
                                            5117003.123], dtype='float32'),
                      'PL610HBA': np.array([3738462.416, 1148244.316,
                                            5021710.658], dtype='float32'),
                      'PL611HBA': np.array([3850980.881, 1438994.879,
                                            4860498.993], dtype='float32'),
                      'PL612HBA': np.array([3551481.817, 1334203.573,
                                            5110157.41], dtype='float32'),
                      'IE613HBA': np.array([3801692.0, -528983.94,
                                            5076958.0], dtype='float32')}


    keys_to_remove = []
    for key in antenna_soltab:
        if key not in station_list:
            keys_to_remove.append(key)

    for k in keys_to_remove:
        antenna_soltab.pop(k, None)

    for a in station_list:
        if a[:2] == 'ST':
            antenna_soltab.update(tied)  # there will only be the tied station
        if a[:2] == 'CS':
            antenna_soltab.update(core)
            break  # only add the core stations to the antenna table once

    return( antenna_soltab )

def insert_into_filestem( infile, insert ):
    tmp = infile.split('/')
    mypath = '/'.join(tmp[0:-1])
    myfile = tmp[-1]
    new_file = os.path.join( mypath, insert+myfile )
    return new_file

# given FITS file, return list of good channels
def find_chan_fits (infile,indisk=1,inseq=1):
    try:
        AIPSUVData('DIFMAP','UVDATA',indisk,inseq).zap()
        print ('Removed DIFMAP.UVDATA.%d disk %d'%(inseq,indisk))
    except:
        pass
    fitld=AIPSTask('fitld')
    fitld.datain=infile if '/' in infile else './'+infile
    fitld.outname = 'DIFMAP'
    fitld.outclass = 'UVDATA'
    fitld.outdisk = indisk
    fitld.outseq = inseq
    fitld.go()
    data = WizAIPSUVData('DIFMAP','UVDATA',indisk,inseq)
    for i in data:
        v = i.visibility
        na = np.asarray(np.ravel(v[:,:,0,0]),dtype='bool')
        try:
            a = na | a
        except:
            a = np.copy(na)
    return a

# given an array of OK channels, write difmap select line
def chan2write (output, a):
    nums=sorted(set(a))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
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

# given MS file, return list of good channels
def find_chan_ms (filename,datacolumn="CORRECTED_DATA",flagname="FLAG"):
    # open ms
    ms=table(filename,ack=False)
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
    return goodchans

# convert corplt output to an amp/phase array - note
# that the UTs are arbitrary **** needs fixing ****
def corplt2array():
    f = open('./CORPLT')
    fc = f.readlines()
    f.close()
    ut=[]
    stn = np.array([],dtype='S')
    for l in fc:
        ls = l.split()
        if not len(ls):
            continue
        if ls[0]=='Stn':
            stn = np.append(stn,ls[3])
        if len(ls)==4 and ls[0] in ['Amp','Phs']:
            ut.append(float(ls[1]))  # faster than numpy append    
    ut = np.unique(np.asarray(ut,dtype='f'))
    amp = np.ones((len(stn),len(ut)))*np.nan
    phs = np.ones((len(stn),len(ut)))*np.nan
    amperr = np.ones((len(stn),len(ut)))*np.nan
    phserr = np.ones((len(stn),len(ut)))*np.nan
    stn_idx = 0
    for l in fc:
        if l[:3] in ['Amp','Phs']:
            t_ut,t_val,t_err = np.asarray(l.split()[1:],dtype='float')
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
def file2fits(infile,datacolumn,flagcol):
    isms = os.path.isdir( infile )
    # Process file itself. If FITS, add '.fits' to the filename and any dchan_
    # file if not present already. Otherwise convert to FITS, replacing any '.ms'
    # or '.MS' extension with '.fits'. Do the same for any associated dchan_ file
    if not isms:
        print( 'File is already in fits format.' )
        if infile[-5:]!='.fits':
            os.system ('mv %s %s.fits'%(infile,infile))
            if os.path.isfile('dchan_%s'%infile):
                os.system ('mv dchan_%s dchan_%s.fits'%(infile,infile))
	    else:
		if HAVE_AIPS:
		    a = find_chan_fits( fitsfile )
	        else:
                    print ('FITS file provided but no AIPS available, so cannot')
                    print ('find the flagged channels, aborting. Either use AIPS')
                    print ('or provide a file called dchan_[infile] with a select')
                    print ('command e.g. select I,1,49,51,100 if channel 50 is bad.')
                    exit(0)
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
            os.system('ms2uvfits in=%s out=%s writesyscal=F'%(infile,fitsfile))
    # Check for presence of a dchan_file. If not, use find_chan_fits or
    # find_chan_ms as necessary

    chan_file = insert_into_filestem( fitsfile, 'dchan_' )
    if not os.path.isfile( chan_file ):
        a = find_chan_ms(infile,datacolumn=datacolumn,flagname=flagcol)
        chan2write(fitsfile,a)
    return fitsfile

# Write and execute the difmap script
def dif_script (infile,pol='XX',aipsno=340,clean_sigma=6,map_size=512,\
                pixel_size=100,obs_length=900,datacolumn='CORRECTED_DATA',\
                startmod=True,flagcol='FLAG'):
    fitsfile = file2fits(infile,datacolumn,flagcol)
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
    fs.write('mapsize map_size,pixel_size\n')
    fs.write('startmod "",1\n' if startmod else 'clean\n')
    fs.write('peakwin 1.5\nselfcal false,false,0\n')
    clean_selfcal_loop (fs,'2,0')
    clean_selfcal_loop (fs,'2,-1')
    clean_selfcal_loop (fs,'0,-1')
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
    fs.write('device %s.ps/vps\nmapl clean,false\n' % \
             local_fits.replace('.fits','_auto%s'%pol))
    fs.write('save %s\nquit\n' % local_fits.replace('.fits','_auto%s'%pol))
    fs.close()
    os.system('difmap <dif_script')
    return fitsfile


#  **** Top-level routine to process an input visibility file ***
# Specified file may be fits or MS
# If MS, will be inspected for good channels and converted to FITS for mapping
#    (any .MS or .ms extension converted to .fits, otherwise .fits added)
#    Corrected_data column assumed unless otherwise specified
# If FITS and we have AIPS/Parseltongue, then inspected for good channels
#    in AIPS and mapped, *If no AIPS and no good-channel file, exit with error*.
# We assume that XX and YY have the same number of telescopes and UTs, and also
#    the same bad/missing channels, and image separately for separate XX/YY corrs
def main( infile, insolfile, clean_sig=6, map_size=512, pix_size=100, obs_length=900, datacolumn='CORRECTED_DATA', startmod=True, flagcolumn='FLAG', verbose=False ):

    ## make a working directory, move the data, and chdir
    tmp = infile.split('/')
    filestem = tmp[-1].split('_')[0] 
    mypath = '/'.join(tmp[0:-1])
    work_dir = os.path.join( mypath, 'difmap_{:s}'.format( filestem ) )
    os.mkdir( work_dir )
    os.system( 'mv {:s} {:s}'.format( infile, work_dir ) )
    infile = os.path.join( work_dir, tmp[-1] )
    os.chdir( work_dir )
    os.system('rm CORPLT')
    ## write a difmap script for the XX polarisation and run
    fitsfile = dif_script(infile, pol='XX', clean_sigma=clean_sig, map_size=map_size, pixel_size=pix_size, obs_length=obs_length, datacolumn=datacolumn, startmod=startmod, flagcol=flagcolumn)
    ## plot solutions and get valies
    ampXX,amperrXX,phsXX,phserrXX,utXX,stnXX = corplt2array()
    os.system('mv CORPLT CORPLT_XX')
    ## write a difmap script for the YY polarisation and run
    fitsfile = dif_script(infile, pol='YY', clean_sigma=clean_sig, map_size=map_size, pixel_size=pix_size, obs_length=obs_length, datacolumn=datacolumn, startmod=startmod, flagcol=flagcolumn)
    ## plot solutions and get values
    ampYY,amperrYY,phsYY,phserrYY,utYY,stnYY = corplt2array()
    os.system('mv CORPLT CORPLT_YY')

    ## COPY FREQUENCY AND TIME FROM LB-Delay-Calibration/solutions.h5 target:TGSSphase
    insols = h5parm( insolfile )
    inss = insols.getSolset('target')
    phasename = [ xx for xx in inss.getSoltabNames() if 'extract' not in xx ][0]
    inst = inss.getSoltab(phasename)
    freq_ax = inst.getAxisValues('freq')
    time_ax = inst.getAxisValues('time')
    insols.close()
    if len(utXX) < len(time_ax):
	## interpolate along axis
        ut = np.float64(utXX)
        utmin = np.float64(ut)
        times = ct.taql('select TIME from {:s} limit 1'.format(infile))
        first_time = times.getcol('TIME')[0]
        time_ax = ut - utmin + first_time

    ## combine XX and YY information and reformat axes
    tmp_amp = np.rollaxis(np.dstack((ampXX,ampYY)),1,0)
    tmp_phs = np.rollaxis(np.dstack((phsXX,phsYY)),1,0)
    ## expand to fill frequency axis
    tmp_amp2 = np.expand_dims( tmp_amp, axis=1 )
    tmp_phs2 = np.expand_dims( tmp_phs, axis=1 )
    amp = np.repeat( tmp_amp2, len(freq_ax), axis=1 )
    phs = np.repeat( tmp_phs2, len(freq_ax), axis=1 )
    phs = phs * -1.0 ## difmap has a different convention

    ## re-reference phases to first antenna in measurement set that isn't flagged
    #ref_ant = find_first_unflagged_ant( infile )
    #ref_idx = [ xx for xx, myant in enumerate( stnXX ) if myant == ref_ant ][0]
    #ref_phs = phs[:,:,ref_idx,:]
    #for xx in range(len(stnXX)):
    #    phs[:,:,xx,:] = phs[:,:,xx,:] - ref_phs

    ## get antenna information
    new_ants = make_ant_table( stnXX )

    ## get pointing information
    ptg = ct.table( infile + '::FIELD' )
    ptg_dir = ptg.getcol('PHASE_DIR')[0]
    ptg.close()
    new_dir = {}
    new_dir[filestem] = ptg_dir[0]

    ## write solutions to an h5parm
    h5parmfile = fitsfile.replace('.fits','_auto.h5')
    out_h5 = h5parm(h5parmfile, readonly = False)
    out_solset = out_h5.makeSolset(solsetName='sol000')
    antenna_table = out_solset.obj._f_get_child('antenna')
    antenna_table.append(new_ants.items())
    out_solset.obj.source.append(new_dir.items())
    out_solset.makeSoltab('amplitude',axesNames=['time','freq','ant','pol'], axesVals=[time_ax,freq_ax,stnXX,['XX','YY']], vals=amp, weights=np.ones_like(amp))   
    out_solset.makeSoltab('phase',axesNames=['time','freq','ant','pol'], axesVals=[time_ax,freq_ax,stnYY,['XX','YY']], vals=phs, weights=np.ones_like(phs))   
    out_h5.close()

    ## rename files so they have unique names
    os.system( 'mv CORPLT_XX {:s}_CORPLT_XX'.format( filestem ) )
    os.system( 'mv CORPLT_YY {:s}_CORPLT_YY'.format( filestem ) )
    os.system( 'mv dif_script {:s}_dif_script'.format( filestem ) )
    os.system( 'mv difmap_log {:s}_difmap_log'.format( filestem ) )
    os.system( 'mv difmap_log1 {:s}_difmap_log1'.format( filestem ) )
    

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
    parser.add_argument("--flagcol",type=str,help="Name of the flag column you want to apply. Default is FLAG.",required=False,default="FLAG")

    ## positionals
    parser.add_argument("filename",type=str,help="Name of the measurement set")

    args=parser.parse_args()

    main( args.filename, clean_sig=args.clean_sig, map_size=args.map_size, pix_size=args.pix_size, 
	obs_length=args.obs_length, datacolumn=args.colname, startmod=args.startmod, flagcolumn=args.flagcol, verbose=args.verbose )
