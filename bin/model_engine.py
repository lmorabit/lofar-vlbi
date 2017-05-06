#--- model_engine.py: makes a (currently) two-component model
#    from u-v data on a particular closure triangle
#           v.1 Neal Jackson, 2015.09.29
#           v.2 NJ, 2017.01.16 converted to CASA, many changes
#           v.3 NJ, 2017.03.07 parallelised
import astropy,numpy as np,scipy,sys,os,glob,warnings,multiprocessing,matplotlib,idx_tels
from matplotlib import pyplot as plt; from scipy import ndimage,optimize
from correlate import *; from mkgauss import *
try:
    import pyrap; from pyrap import tables as pt
    have_tables = True
except:
    have_tables = False

plt.rcParams['image.origin'],plt.rcParams['image.interpolation']='lower','nearest'
warnings.simplefilter('ignore')
RAD2ARC,LIGHT,FIRSTNPY = 3600.*180./np.pi, 2.99792458E+8, './first_2008.simple.npy'
MX,MY,MF,MW,MR,MP = range(6)
ISPARALLEL = int(os.popen('nproc').read())>4
CASAPY = '/pkg/casa-release-4.7.0-1-el6/bin/casa'

# Make a movie from a set of png files
def movie (fps=4):
    a=np.sort(glob.glob('model_engine_*.png'))
    for i in range(len(a)):
        os.system('convert %s model_engine%03d.jpeg'%(a[i],i))
    command = "mencoder \"mf://model_engine*.jpeg\" -mf fps=%d -o model_engine.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=4800" %fps
    os.system(command)
    os.system('rm model_engine*.jpeg')

# Work out real and imaginary parts of a list of visibilities given model and uvw's
def uvw2reim (uvw, model):
    u,v,re,im = uvw[:,0],uvw[:,1],0.,0.
    for m in model:     # note - sign in m[MX] as this is HA not RA
        cmpamp,cmpphs = m[MF],2.*np.pi*(-u*m[MX]+v*m[MY])/RAD2ARC
        if m[MW]!=0.0:
            sphi,cphi = np.sin(np.deg2rad(m[MP])),np.cos(np.deg2rad(m[MP]))
            tc = np.pi*(m[MW]/RAD2ARC)*np.hypot(v*cphi+u*sphi,m[MR]*(u*cphi-v*sphi))
            cmpamp = m[MF]*np.exp(-0.3696737602*tc*tc)
        re += cmpamp * np.cos(cmpphs)
        im += cmpamp * np.sin(cmpphs)
    return re,im

def taql_calc (vis, vistable, qtab, qtype):
    os.system('taql \'CALC '+qtype+' ([select '+qtab+' from '+vis+\
              '/'+vistable+'])\' >taql_out')
    f=open('taql_out')
    for v in f:
        try:
            val = float(v.rstrip('\n'))
        except:
            pass
    f.close()#; os.system('rm taql_out')
    return val

def taql_num (vis,vistable,qtab):
    os.system('taql \'select '+qtab+' from '+vis+'/'+vistable+'\' >taql_out')
    f=open('taql_out')
    for v in f:
        if 'select result of' in v:
            n = int(v.split('of')[1].split('row')[0])
            break
    return n

def taql_from (vis,vistable,qtab):
    os.system('taql \'select '+qtab+' from '+vis+'/'+vistable+'\' >taql_out')
    f = open('taql_out')
    for v in f:
        pass
    f.close()
    return v.rstrip('\n').rstrip(']').lstrip('[').split(',')

def dget_c (vis, tel1, tel2):
    f,fo = open('model_dget.py'), open('temp.py','w')
    for line in f:
        fo.write(line)
    f.close()
    os.system('rm d.npy; rm ut.npy; rm uvw.npy')
    fo.write('d,ut,uvw=dget (\''+vis+'\','+str(tel1)+','+str(tel2)+')\n')
    fo.write('np.save(\'d\',d)\n')
    fo.write('np.save(\'ut\',ut)\n')
    fo.write('np.save(\'uvw\',uvw)\n')
    fo.close()
    os.system(CASAPY+' --nologger -c temp.py')
    d,ut,uvw = np.load('d.npy'),np.load('ut.npy'),np.load('uvw.npy')
    return d,ut,uvw

def dget_t (vis, tel1, tel2):
    os.system('taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' % (vis, tel1, tel2, 'cl_temp.ms'))
    t = pt.table('cl_temp.ms')
    ut = np.ravel(np.asarray([tuple(each.values()) for each in t.select('TIME')]))
    spw = np.ravel(np.asarray([tuple(each.values()) for each in t.select('DATA_DESC_ID')]))
    dc = t.select('DATA')
    d = np.asarray([tuple(each.values()) for each in dc])[:,0,:,:]
    d = np.swapaxes(d,0,2)     # makes the order pol - chan - time as in casa
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw,u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    if spw.sum():
        for i in np.unique(spw):
            new = np.take(d,np.argwhere(spw==i),axis=2)[:,:,:,0]
            try:
                d_out = np.concatenate((d_out,new),axis=1)
            except:
                d_out = np.copy(new)
        d = d_out
    return d,ut,uvw

def norm(a,isred=True):
    nlim = np.pi
    a = a%(2*nlim) if isred else a
    np.putmask(a,a>nlim,a-2.*nlim)
    np.putmask(a,a<-nlim,a+2.*nlim)
    return a

def getap (d,pol=0):
    ph = np.sum(d[pol],axis=0)
    return np.sum(abs(d[pol]),axis=0)/d.shape[1],np.arctan2(ph.imag,ph.real)


def get_uvw_table (t):
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw,u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    return uvw

# Get data and u-v arrays from a measurement set on a given triangle
def data_extract (vis):
    global uvw01,uvw02,uvw12,cp012,a01,a02,a12
    chw = taql_calc(vis,'SPECTRAL_WINDOW','CHAN_WIDTH','mean')
    ch0 = taql_calc(vis,'SPECTRAL_WINDOW','REF_FREQUENCY','mean')
    nchan = taql_calc(vis,'SPECTRAL_WINDOW','NUM_CHAN','mean')
    sra,sdec = taql_from(vis,'FIELD','PHASE_DIR')
    nspw = taql_num(vis,'SPECTRAL_WINDOW','NUM_CHAN')
    sra = np.asarray(sra.replace('h',' ').replace('m',' ').split(),dtype='f')
    sdec = np.asarray(sdec.replace('d',' ').replace('m',' ').split(),dtype='f')
    ra = 15.0*(sra[0]+sra[1]/60.0+sra[2]/3600.0)
    dec = sdec[0]+sdec[1]/60.0+sdec[2]/3600.0
    wlength = np.mean(LIGHT/(ch0 + chw*np.arange(nspw*nchan)))
    itel = np.array(idx_tels.get_idx_tels (vis,trname))
    # order of telescopes on baselines 0-1, 0-2, 1-2; -1 if in "wrong" order
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    btel = [min(itel[0],itel[1]),max(itel[0],itel[1]),min(itel[0],itel[2]),\
            max(itel[0],itel[2]),min(itel[1],itel[2]),max(itel[1],itel[2])]
    if have_tables:
        d01,ut01,uvw01 = dget_t (vis, btel[0],btel[1])
        d02,ut02,uvw02 = dget_t (vis, btel[2],btel[3])
        d12,ut12,uvw12 = dget_t (vis, btel[4],btel[5])
    else:
        d01,ut01,uvw01 = dget_c (vis, btel[0],btel[1])
        d02,ut02,uvw02 = dget_c (vis, btel[2],btel[3])
        d12,ut12,uvw12 = dget_c (vis, btel[4],btel[5])
    a01,p01 = getap(d01)
    a02,p02 = getap(d02)
    a12,p12 = getap(d12)
    cp012 = otel[0]*p01-otel[1]*p02+otel[2]*p12
    np.putmask(cp012,cp012>np.pi,cp012-2*np.pi)
    np.putmask(cp012,cp012<-np.pi,cp012+2*np.pi)
    uvw01 /= wlength
    uvw02 /= wlength
    uvw12 /= wlength
    print trname,'-> antenna numbers:',itel
    print 'Baseline lengths: %s-%s: %dkm %s-%s: %dkm %s-%s: %dkm' % \
       (trname[0],trname[1],int(np.sqrt((uvw01[0]**2).sum())*wlength/1000),\
        trname[0],trname[2],int(np.sqrt((uvw02[0]**2).sum())*wlength/1000),\
        trname[1],trname[2],int(np.sqrt((uvw12[0]**2).sum())*wlength/1000))
    os.system('rm -fr cl_tmp*.ms')
    return itel,np.mean(wlength),ra,dec

# --------------------------------------------------

def model_extract (model,itel):
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    re01,im01 = uvw2reim (uvw01,model)
    re12,im12 = uvw2reim (uvw12,model)
    re02,im02 = uvw2reim (uvw02,model)
    ph01 = norm(np.arctan2(im01,re01))
    ph02 = norm(np.arctan2(im02,re02))
    ph12 = norm(np.arctan2(im12,re12))
    clph = norm(ph01*otel[0] - ph02*otel[1] + ph12*otel[2])
    return np.hypot(re01,im01),np.hypot(re02,im02),clph

def plotimg (A01,A02,CP012,model,goodness,itel,aplot,gcou):
    ells = []
    for i in model:
        if i[MW]!=0.0:
            ells.append(matplotlib.patches.Ellipse(xy=[i[MX],i[MY]],width=2.*i[MW],\
                height=2.*i[MW]*i[MR],lw=0.5,angle=i[MP]+90.,fill=0.0))
        else:
            ells.append(matplotlib.patches.Ellipse(xy=[i[MX],i[MY]],width=0.1,\
                height=0.1,lw=0.5,angle=0.,fill=1.0,color='red'))
    if plottype in [1,10]:
        plt.subplot2grid((6,5),(0,0),rowspan=2)
        scalarray = np.sort(np.ravel(a01)); ls = len(scalarray)
        amp_min, amp_max = scalarray[int(0.1*ls)], scalarray[int(0.9*ls)]
        plt.imshow(a01,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(0,1),rowspan=2)
        plt.imshow(A01,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(0,2),rowspan=2)
        plt.imshow(a01-A01,aspect='auto',\
                   vmin=2*amp_min-amp_max,vmax=2*amp_max-amp_min)
        plt.subplot2grid((6,5),(2,0),rowspan=2)
        scalarray = np.sort(np.ravel(a02)); ls = len(scalarray)
        amp_min, amp_max = scalarray[int(0.1*ls)], scalarray[int(0.9*ls)]
        plt.imshow(a02,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(2,1),rowspan=2)
        plt.imshow(A02,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(2,2),rowspan=2)
        plt.imshow(a02-A02,aspect='auto',\
                   vmin=2*amp_min-amp_max,vmax=2*amp_max-amp_min)
        plt.subplot2grid((6,5),(4,0),rowspan=2,yticks=[])
        plt.imshow(cp012,aspect='auto',vmin=-180,vmax=180)
        plt.subplot2grid((6,5),(4,1),rowspan=2,yticks=[])
        plt.imshow(CP012,aspect='auto',vmin=-180,vmax=180)
        plt.subplot2grid((6,5),(4,2),rowspan=2,yticks=[])
        plt.imshow(CP012-cp012,aspect='auto')
        plt.subplot2grid((6,5),(3,3),rowspan=3,colspan=3,xticks=[],yticks=[])
        plt.imshow(aplot,cmap=matplotlib.cm.gray_r)
        plt.title('X2=%.3f'%goodness)
    elif plottype in [2,20]:
        plt.subplot2grid((6,5),(0,0),rowspan=2,colspan=3)
        plt.plot(a01,'b-'); plt.plot(A01,'r-')
        plt.legend([trname[0]+'-'+trname[1]],fontsize=6)
        plt.subplot2grid((6,5),(2,0),rowspan=2,colspan=3)
        plt.plot(a02,'b-'); plt.plot(A02,'r-')
        plt.legend([trname[0]+'-'+trname[2]],fontsize=6)
        plt.subplot2grid((6,5),(4,0),rowspan=2,colspan=3)
        plt.plot(cp012,'b-'); plt.plot(CP012,'r-')
        plt.subplot2grid((6,5),(4,3),rowspan=2,colspan=3,xticks=[],yticks=[])
        plt.imshow(aplot,cmap=matplotlib.cm.gray_r)
        plt.title('X2=%.3f'%goodness)
    ax = plt.subplot2grid((6,5),(0,3),rowspan=3,colspan=2)
    for i in range(len(ells)):
        e = ells[i]
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
    ax.set_xlim(-glim,glim);ax.set_ylim(-glim,glim)
    plt.grid()
    if gcou==-1:   # find how many png files, make the next one
        gcou = len(glob.glob('model_engine_*.png'))
    plt.title('model_engine_%03d.png'%gcou)
    if plottype in [10,20]:
        plt.savefig('model_engine_%03d.png'%gcou)

#===================================================
def ndiff (a,b):
    sqd = np.array([])
    for i in range(-len(a)/2,len(a)/2):
        a1 = np.roll(a,i)
        idx1,idx2 = max(0,i),min(len(a),len(a)+i)
        sqd = np.append(sqd, np.mean(((b-a1)[idx1:idx2]**2)))
    return sqd

# this is the bit that may need replacing - maybe fit splines and look at the spline points?

def get_goodness(A01,A02,CP012):
    beta = 0.00001    #   this is a pretty vital parameter
    ascat = np.median(abs(np.gradient(np.ravel(a02))))    
    cscat = np.median(abs(np.gradient(np.ravel(cp012))))
    sq = ndiff(a01*np.nanmean(A01)/np.nanmean(a01),A01)/ascat**2 + \
         ndiff(a02*np.nanmean(A02)/np.nanmean(a02),A02)/ascat**2 + \
         ndiff(cp012,CP012)/cscat**2
    difmin = 0.5*len(a01)-np.argwhere(sq==sq.min())[0][0]
    return sq.min() + beta*difmin**2

def mod_func (x0, *x):
    model,opt,itel,aplot,gcou,iy,ix = x
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = x0
    model = model.reshape(len(model)/6,6)
    A01,A02,CP012 = model_extract (model,itel)
    if ampfiddle:
        A01 *= np.median(a01)/np.median(A01)
        A02 *= np.median(a02)/np.median(A02)
    goodness = get_goodness(A01,A02,CP012)
    if plottype:
        plotimg (A01,A02,CP012,model,goodness,itel,aplot,gcou)
        if plottype in [1,2]:
            plt.draw()
            plt.pause(0.001)
            plt.clf()
    return goodness

def getmodel(coord,beam):
    gsiz = int(gridsize/(bsub*beam))   # original grid very big
    ginc = bsub*beam                   # arcsec/pix grid size
    grid = np.ones((gsiz,gsiz))*np.nan # original grid full of NaN
    first = np.load(FIRSTNPY)
    pflux,pcoord = [],[]
    a = correlate (np.array([coord]),0,1,first,0,1,0.5/60)
    cosdec = np.cos(np.deg2rad(coord[1]))
    print 'Found: %d FIRST sources'%len(a)
    for i in range(len(a)):
        fi = first[int(a[i,1])]
        scoord = astropy.coordinates.SkyCoord(fi[0],fi[1],unit='degree')
        scoord = scoord.to_string(style='hmsdms')
        scoord = scoord.replace('d',':').replace('m',':')
        scoord = scoord.replace('h',':').replace('s','')
        print '%s flux=%.1f shape=(%.1fx%.1f PA%.1f)'%(scoord,fi[2],fi[4],fi[5],fi[6])
        pix = 0.5*gsiz+3600.*(fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        try:
            pcoord.append(fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        except:
            pcoord = (fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        pflux.append(fi[2])
        new = mkgauss([gsiz,gsiz],pix,fi[2],fi[4]/ginc,fi[5]/fi[4],fi[6])
        np.putmask(grid,new>0.004*new.max(),i)  # decrease 0.01 if missing cpts
    # grid now consists of original very big grid, with numbers instead of NaN where
    # the secondary might be
    if pflux:    # shrink the grid around the Gaussian near FIRST sources
        while all(np.isnan(grid[0,:])) and all(np.isnan(grid[-1,:])) and \
              all(np.isnan(grid[:,0])) and all(np.isnan(grid[:,-1])):
            grid = grid[1:-1,1:-1]
    return grid,ginc,grid.shape[0],pflux,pcoord

def parallel_function(f):
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        ncores = max(1,int(os.popen('nproc').read())-1)  # use all cores - 1
        print 'Using',ncores,'cores'
        pool = Pool(processes=ncores) # depends on available cores
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = np.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    from functools import partial
    return partial(easy_parallize, f)

def grid_search_thread (k):
    a =  mod_func(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7])
    f=open('grid_search_thread.log','a')
    f.write('%d %d %f\n'%(k[-2],k[-1],a))
    f.close()
    return a

def grid_search (model,cpt,gridcpt,itel,aplot,gcou,grid,gsiz,ginc,isparallel=ISPARALLEL):
    opt = np.zeros_like(model,dtype='bool')
    args = []
    import copy
    for ix in range(gsiz):
        x = ginc*(ix-gsiz/2.0)   # x,y in arcsec; a in ginc-size pixels
        for iy in range(gsiz):
            if grid[iy,ix] == gridcpt:
                y = ginc*(iy-gsiz/2.0)
                model[cpt][0],model[cpt][1] = x,y
                if isparallel:
                    arg = ([],model,opt,itel,aplot,gcou,iy,ix)
                    args.append(copy.deepcopy(arg))   # otherwise overwrites all elements
                    gcou += 1
                else:
                    aplot[iy][ix] = mod_func ([],model,opt,itel,aplot,gcou,iy,ix)
                    gcou += 1
    if isparallel:
        print 'Starting grid search with',len(args),'points'
        os.system('rm grid_search_thread.log')
        grid_search_thread.parallel = parallel_function(grid_search_thread)
        parallel_result = grid_search_thread.parallel (args)
        for i in range(len(parallel_result)):
            aplot[args[i][-2],args[i][-1]] = parallel_result[i]
    np.putmask(aplot,np.isnan(aplot),np.nanmax(aplot))
    model[cpt,:2] = ginc*(np.asarray(ndimage.measurements.minimum_position \
            (aplot)[::-1])-0.5*np.asarray(grid.shape))
    return model,aplot

def recentroid (model,startcpt,endcpt,ginc,gsiz):
    pos = model[startcpt:endcpt+1,:2]*model[startcpt:endcpt+1,2]
    centroid = np.sum(pos,axis=0)/len(pos)
    model[startcpt:endcpt+1,:2] -= centroid
    return model

def adjust_all (model,itel,aplot,gcou):
    opt = np.ones_like(model,dtype='bool')
    x0 = np.ravel(model)[np.ravel(opt)]
    args = (model,opt,itel,aplot,-1,0,0)
    xopt = optimize.fmin(mod_func, x0, args=args, maxiter=100)
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = xopt
    model = model.reshape(len(model)/6,6)
    return model

def refine_points (model,cpts,itel,aplot,gcou):
    opt = np.zeros_like(model,dtype='bool')
    for i in cpts:
        opt[i,2:] = True
    x0 = np.ravel(model)[np.ravel(opt)]
    args = (model,opt,itel,aplot,-1,0,0)
    xopt = optimize.fmin(mod_func, x0, args=args, maxiter=20)
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = xopt
    model = model.reshape(len(model)/6,6)
    return model

def write_skymodel (ra,dec,model,outname):
    if outname!='':
        f = open(outname,'w')
        f.write ('# (Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, Orientation) = format\n')
    for i in range(len(model)):
        cosd = 3600.*np.cos(np.deg2rad(dec))
        s = astropy.coordinates.SkyCoord(ra-model[i,0]/cosd,dec+model[i,1]/3600,unit='degree')
        s = s.to_string(style='hmsdms')
        sra = s.split()[0].replace('h',':').replace('m',':').replace('s','')
        sdec = s.split()[1].replace('d','.').replace('m','.').replace('s','')
        print 'ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f'%(i,sra,sdec,model[i,MF],\
              model[i,MW],model[i,MW]*model[i,MR],np.rad2deg(model[i,MP]))
        if outname!='':
            f.write('ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f\n'%(i,sra,sdec,model[i,MF],\
                  model[i,MW],model[i,MW]*model[i,MR],np.rad2deg(model[i,MP])))
    if outname!='':
        f.close()

#============= main script ========================
def mainscript(vis,TRNAME,BSUB=0.3,GRIDSIZE=12.0,PLOTTYPE=20,AMPFIDDLE=True,outname='model_engine.sky'):
    global bsub,gridsize,plottype,ampfiddle,glim,trname
    bsub,gridsize,plottype,ampfiddle,trname,gcou = BSUB,GRIDSIZE,PLOTTYPE,AMPFIDDLE,TRNAME,0
    os.system('rm model_engine*.png')
    itel,wv,ra,dec = data_extract (vis)
    s_amp = np.sort(np.ravel(a01)); ls = len(s_amp)
    flux1,flux2 = np.median(s_amp), 0.5*(s_amp[int(0.99*ls)]-s_amp[int(0.01*ls)])
    flux = np.median(a01)
    beam = RAD2ARC/np.nanmax(abs(uvw01)); print 'Beam:',beam,'arcsec'
    grid,ginc,gsiz,pflux,pcoord = getmodel (np.array([ra,dec]),beam)
    if not len(pflux):   # no FIRST source, search the whole grid
        grid = np.zeros_like (grid)
    aplot = np.ones_like(grid)*np.nan
    glim = 0.5*ginc*gsiz
    if plottype in [1,2]:
        plt.ion()
        fig = plt.figure()
    if len(pflux) < 2:    # zero, or one FIRST source
        model = np.array([[0.0,0.0,flux1,0.5,1.0,0.0],[-1.0,-2.0,flux2,0.5,1.0,0.0]])
        model,aplot = grid_search (model,1,0,itel,aplot,gcou,grid,gsiz,ginc)   
# fix one cpt, grid-search posn of 2nd.
        print 'Model after grid search:'
        write_skymodel (ra,dec,model,'')
        model[0,MW] = model[1,MW] = 1.0
        model = refine_points (model,[0,1],itel,aplot,gcou)  
# fit for size/orientation only
        print 'Model after refining points:'
        write_skymodel (ra,dec,model,'')
        model = recentroid (model,0,1,ginc,gsiz)  # put centroid in centre of image
        print 'Model after recentroiding:'
    else:   # >1 FIRST source - not tested yet
        model = np.zeros((len(pflux),6))
        for i in range(len(pflux)):
            model[i,:2] = pcoord[i]
            model[i,3] = pflux[i]*flux/pflux.sum()
            model[i,4:] = [0.5,1.0,0.0]
        model = adjust_all(model,itel,aplot,gcou)
    if plottype in [10,20]:
        movie()
    write_skymodel (ra,dec,model,'')
    write_skymodel (ra,dec,model,outname)
#mainscript('./SIM5.ms', ['ST001','DE601','DE605HBA'])
#mainscript('./PP1_av_corr.ms', ['ST001','DE601HBA','DE605HBA'])
