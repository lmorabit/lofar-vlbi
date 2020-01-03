#!/usr/bin/python
import numpy as np,os,glob,sys
import losoto; from losoto.h5parm import h5parm
from pyrap.tables import table
import argparse

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
def main( infile, clean_sig=6, map_size=512, pix_size=100, obs_length=900, datacolumn='CORRECTED_DATA', startmod=True, flagcolumn='FLAG', verbose=False ):

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
    ## combine XX and YY information
    amp = np.rollaxis(np.dstack((ampXX,ampYY)),2,0)
    amperr = np.rollaxis(np.dstack((amperrXX,amperrYY)),2,0)
    phs = np.rollaxis(np.dstack((phsXX,phsYY)),2,0)
    phserr = np.rollaxis(np.dstack((phserrXX,phserrYY)),2,0)
    ut,stn = utXX,stnXX
    ## write solutions to an h5parm
    h5parmfile = fitsfile.replace('.fits','_auto.h5')
    data = h5parm(h5parmfile, readonly = False)
    outSolset = data.makeSolset('sol000')
    outSolset.makeSoltab(soltype='amplitude',soltabName='amplitude000',axesNames=['pol','ant','time'],axesVals=[['XX','YY'],stn,ut],vals=amp,weights=np.ones_like(amp))
    outSolset.makeSoltab(soltype='phase',soltabName='phase000',axesNames=['pol','ant','time'],axesVals=[['XX','YY'],stn,ut],vals=phs,weights=np.ones_like(phs))
    data.close()
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
