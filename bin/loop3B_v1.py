#!/usr/bin/env python

import numpy as np,os,sys,glob,astropy,argparse,time
import casacore.tables as casatb
import astropy.io.fits as pyfits
import loop3_serviceB as loop3_service
CCRIT=1.6
TSAMP=4.0    # time per sample. All times are in seconds. TSAMP is passed to NDPPP for writing
             # the parset, since solints in NDPPP parsets are in samples.


# Sean: This does not run in screen because matplotlib crashes as it can't connect to local host.

def cleanup(vis):
    os.system('rm -fr %s_processing'%vis)
    os.system('mkdir %s_processing'%vis)
    h5all = np.sort(glob.glob(vis+'*_*c0.h5'))
    h5A = np.sort(glob.glob(vis+'*A_*c0.h5'))
    h5vis = h5all[~np.in1d(h5all,h5A)]
    os.system('mv %s*.fits %s_processing'%(vis,vis))
    os.system('mv %s*.log %s_processing'%(vis,vis))
    os.system('mv %s*.png %s_processing'%(vis,vis))
    os.system('mv %s*.h5 %s_processing'%(vis,vis))
    calfiles = np.array([])
    if len(h5vis):
        h5vis_head, h5vis_tail = os.path.split(h5vis[-1])  # Sean added
        os.system('cp %s_processing/%s .'%(vis,h5vis_tail))  # Sean edited
        calfiles = np.append(calfiles,h5vis[-1])
    if len(h5A):
        h5A_head, h5A_tail = os.path.split(h5A[-1])  # Sean added
        os.system('cp %s_processing/%s .'%(vis,h5A_tail))  # Sean edited
        calfiles = np.append(calfiles,h5A[-1])
    vis_head, vis_tail = os.path.split(vis)
    os.system('cp %s_processing/%s_output.png .'%(vis,vis_tail))

    return '%s_output.png'%vis, calfiles

def imaging(vis,niters,threshold):
    loop3_service.imagr (vis,cellsize='0.05asec',imsize=1024,maxuvl=1000000,\
          gain=0.1,mgain=0.85,dostopnegative=True,niter=niters,\
          autothreshold=threshold,weightingrankfiltersize=256,\
          weightingrankfilter=3,domultiscale=True,automask=7.0,\
          outname='test',dolocalrms=True)

# Needs an initial model. May be provided as:
#    model=None        Make an image and use that. (Unlikely to be a good idea)
#    model='MODEL'     Look in the MODEL_DATA column
#    model=[filename]  LOFAR sourcedb format e.g. converted from FIRST, LoTSS, EHT imager....
#                      [Not working yet, maybe require calling routine to do this?]
#  NOTE: uses bdsf - version 1.8.13 which loads by default has a conflict with
#  other libraries - may need to unload and use 1.8.10 instead 

def selfcal(vis,model='MODEL',outcal_root='',max_sol=600.0,init_sol=30.0,\
            incol='DATA',outcol='DATA',caltype='P'):
    if not model:
        imaging(vis,1000,10)
    # need a predict step to deal with sourcedb here if necessary
    ptant = casatb.table(vis+'/ANTENNA')
    antenna_list = np.array([],dtype='S')
    for i in ptant.select('NAME'):
        antenna_list = np.append(antenna_list,i.values())
    nant = len(antenna_list)
    if caltype=='P':
        sol_int_range = np.arange(np.ceil(np.log(max_sol/init_sol)/np.log(3.)))
        sol_int_range = np.ceil(init_sol*3.**sol_int_range)
        nsol = len(sol_int_range)
        coh = CCRIT*np.ones((nsol,nant))
        for i in range(nsol):
            solint = sol_int_range[i]
            outcal_root = outcal_root if len(outcal_root) else vis
            outcal = outcal_root+'_c%d.h5'%i
            loop3_service.loop3log (vis,'Beginning pass with solint %.1f sec\n' % (solint))
            loop3_service.calib (vis, solint=solint, outcal=outcal, incol=incol, \
                                 outcol=outcol,solmode='P',tsamp=TSAMP)
            loop3_service.snplt (vis,htab=outcal,outpng=outcal)
            coh[i] = loop3_service.coherence_metric (outcal)
            loop3_service.loop3log(vis,'Coherences: \n')
            for j in range(nant):
                loop3_service.loop3log(vis,'%s:%f '%(antenna_list[j],coh[i,j]),cret=False)
            if len(coh[i][coh[i]>=CCRIT])==0:
                break
    # For each antenna in the antenna list, find the selfcal table with 
    # the shortest solution interval that contains coherent solutions. If 
    # there is no such table, report -1 in order to signal that they should 
    # all be set to zero.
        ncoh = np.ones(nant,dtype=int)*-1
        allcoh = np.ones(nant,dtype=float)*CCRIT
        for i in range(nant):
            try:
                ncoh[i] = np.min(np.ravel(np.argwhere(coh[:,i]<CCRIT)))
                allcoh[i] = coh[:,i][ncoh[i]]
            except:
                pass
    # For each selfcal table containing the shortest solution interval with 
    # coherence on some antennas, replace the entries in the first selfcal 
    # table with the interpolated values from that antenna
        for i in range(1,coh.shape[0]):
            iant = antenna_list[ncoh==i]
            if len(iant):
                loop3_service.clcal (outcal_root+'_c0.h5',outcal_root+'_c%d.h5'%i,\
                                     ant_interp=iant)
    # For each antenna without any coherence at all, zero the phase 
    # solutions for that antenna
        iant = antenna_list[ncoh==-1]
        if len(iant):
            loop3_service.zerosol (outcal_root+'_c0.h5',iant)
    else:    # amplitude selfcal: only one interval
        outcal_root = outcal_root if len(outcal_root) else vis
        outcal = outcal_root+'_c0.h5'%i
        loop3_service.calib (vis, solint=init_sol, outcal=outcal, incol=incol, \
                                 outcol=outcol,solmode='A')
        loop3_service.snplt (vis,htab=outcal,outpng=outcal,soltab='amplitude000')
        allcoh = loop3_service.coherence_metric (outcal)
        loop3_service.loop3log(vis,'Coherences: \n')
        for i in range(nant):
            loop3_service.loop3log(vis,'%s:%f '%(antenna_list[i],allcoh[i]),cret=False)
    # For each antenna without any coherence at all, zero the amp/phase 
    # solutions for that antenna
        iant = antenna_list[allcoh<CCRIT]
        if len(iant):
            loop3_service.zerosol (outcal_root+'_c0.h5',ant=iant)
    # find the maximum baseline length with coherent cal signal
    cohlength = loop3_service.getcoh_baseline (antenna_list,allcoh,CCRIT)
    return allcoh,cohlength

# following is based on Frits's algorithm with measure_statistic

def measure_statistic (filename):
    img = pyfits.open(filename)[0].data.squeeze()
    return abs (img.max()/img.min())

def applycal_split (vis, visA, solset, parmdb, soltab='phase000',\
                    correction='phase000'):
    fo = open('applycal.parset','w')
    fo.write ('msin=%s\n'%vis)
    fo.write ('msout=%s\n'%vis)
    fo.write ('msin.datacolumn=DATA\n')
    fo.write ('msout.datacolumn=CORRECTED_DATA\n')
    fo.write ('steps=[applycal]\n')
    fo.write ('applycal.type=applycal\n')
    fo.write ('applycal.parmdb=%s\n'%parmdb)
    fo.write ('applycal.solset=%s\n'%solset)
    fo.write ('applycal.soltab=%s\n'%soltab)
    fo.write ('applycal.correction=%s\n'%correction)
    fo.close()
    os.system ('NDPPP applycal.parset')
    fo = open('split.parset', 'w')
    fo.write ('msin=%s\n'%vis)
    fo.write ('msin.datacolumn=CORRECTED_DATA\n')
    fo.write ('msout=%s\n'%visA)
    fo.write ('msout.datacolumn=DATA\n')
    fo.write ('steps=[]\n')
    fo.close()
    os.system ('NDPPP split.parset')

def main (vis,strategy='P30,P30,P30,A500,A450,A400',startmod='',ith=5.0,bandwidth='8MHz'):

    ## process arguments
    vis = vis.rstrip('/')  # remove trailing slash so the files are not created inside the MS
    strategy = strategy.split(',')
    bw_val = ''
    bw_unit = ''
    for c in bandwidth:
	try:
	    float(c)
	    bw_val = bw_val + c
	except ValueError:
	    bw_unit = bw_unit + c
    if bw_unit == 'MHz':
	bw_val = float(bw_val)*1e6
    ## get bandwidth of vis
    spec_info = casatb.table( vis + '/SPECTRAL_WINDOW')
    total_bw = spec_info.getcol('TOTAL_BANDWIDTH')[0]
    spec_info.close()
    if total_bw < bw_val:
	wsclean_chans = 0
	mfs = ''
    else:
	wsclean_chans = int( np.ceil(total_bw/bw_val) )
	mfs='-MFS'

    ## make a temporary directory to save things to
    tmp_dir = 'loop3_'+vis.rstrip('.ms').rstrip('.MS')
    os.system('mkdir %s'%tmp_dir)
    os.chdir(tmp_dir)
    os.system('mv ../%s .'%vis)


    import bdsf
    prevstat = 0.0
    cohlength = 2.0E6
    strategy_type = []
    for i in strategy: 
        strategy_type.append(i[0])
    ploop, nloop, snver = strategy_type.count('P'), len(strategy), 0
    #
    # PHASE CALIBRATION - run through ploop iterations, exiting if we have convergence
    #
    for iloop in range(ploop):
        fitsmask = vis+'_%02d-mask.fits'%(iloop-1) if iloop else ''
        if startmod=='' or iloop:
            pstr = '******* PHASE LOOP %d running wsclean ************'%iloop
            loop3_service.loop3log (vis, pstr+'\n')
            loop3_service.imagr(vis,cellsize='0.05asec',domultiscale=True,\
                  outname=vis+'_%02d'%iloop,channelsout=wsclean_chans,robust=-1,\
                  fitsmask=fitsmask,dolocalrms=True,maxuvwm=cohlength)
        else:
            # Need something here to produce an image from startmod
            pass
        pstr='******* PHASE LOOP %d making mask %s_%02d%s-image.fits ********'%(iloop,vis,iloop,mfs)
        loop3_service.loop3log (vis, pstr+'\n')
        stdout = sys.stdout; sys.stdout = open('bdsf_chunterings','a')
        img=bdsf.process_image('%s_%02d%s-image.fits'%(vis,iloop,mfs),atrous_do=True,thresh_isl=ith)
        sys.stdout.close(); sys.stdout = stdout
        img.export_image(img_type='island_mask',outfile='%s_%02d-mask.fits'%(vis,iloop))
        thisstat = measure_statistic(vis+'_%02d%s-image.fits'%(iloop,mfs))
        ####### need a line here to bomb out if no coherence
        # exit loop if clean finishing
        pstr='******* PHASE LOOP %d goodness stat %f ************' % (iloop,thisstat)
        loop3_service.loop3log (vis, pstr+'\n')
        if thisstat-prevstat<0.01:
            pstr='****** EXITING PHASE CAL with diff %f *********'%(thisstat-prevstat)
            loop3_service.loop3log (vis, pstr+'\n')
            break
        else:   
            prevstat = thisstat
            loop3_service.imagr(vis,dopredict=True,fitsmask=fitsmask,autothreshold=3,dolocalrms=True,\
                                robust=-1,outname=vis+'_%02d%s'%(iloop,mfs))
        pstr='******* PHASE LOOP %d making new cal file %s ************' % (iloop,vis+'_%02d'%iloop)
        loop3_service.loop3log (vis, pstr+'\n')
        caltype, sol0 = strategy[iloop][0], float(strategy[iloop][1:])
        coh, cohlength = selfcal(vis,model='MODEL',incol='DATA',outcol='CORRECTED_DATA',\
                      outcal_root=vis+'_%02d'%iloop,caltype=caltype,init_sol=sol0)
        snver = iloop
        pstr='******** END PHASE LOOP %d - coherence on %.1f km **********' % \
              (iloop,cohlength/1000.)
        loop3_service.loop3log (vis, pstr+'\n')
    # Exit at this point if we are not doing amplitude cal
    if ploop == nloop:
        exit()
    #
    # If we are doing amplitude calibration, we now need to apply the 
    # calibration and write a new MS with a DATA column
    visA = vis+'_A'
    # delete all existing files beginning with vis+'_A'
    os.system('rm -fr %s*'%visA)
    pstr='****** APPLYING CALIBRATION TABLE %d\n'%snver
    loop3_service.loop3log (vis, pstr+'\n')
    applycal_split (vis, visA, 'sol000', '%s_%02d_c0.h5' % (vis,snver))
    init_fitsmask = vis+'_%02d-mask.fits'%iloop
    init_img = vis+'_%02d%s-image.fits'%(iloop,mfs)
    pred_img = vis+'_%02d%s'%(iloop,mfs)
    for iloop in range(ploop,nloop):
        fitsmask = init_fitsmask if iloop==ploop else visA+'_%02d-mask.fits'%(iloop-1)
        pstr='******* AMPLITUDE LOOP %d running wsclean ************'%iloop
        loop3_service.loop3log (vis, pstr+'\n')
        loop3_service.imagr(visA,cellsize='0.05asec',domultiscale=True,\
                  outname=visA+'_%02d'%iloop,channelsout=wsclean_chans,robust=-1,\
                  fitsmask=fitsmask,dolocalrms=True,maxuvwm=cohlength)
        image_bdsf = '%s_%02d%s-image.fits'%(visA,iloop,mfs)
        pstr='******* AMPLITUDE LOOP %d making mask %s_%02d%s-image.fits ************'%(iloop,visA,iloop,mfs)
        loop3_service.loop3log (vis, pstr+'\n')
        img=bdsf.process_image(image_bdsf,atrous_do=True,thresh_isl=5.0)
        img.export_image(img_type='island_mask',outfile='%s_%02d-mask.fits'%(visA,iloop))
        thisstat = measure_statistic(visA+'_%02d%s-image.fits'%(iloop,mfs))
        pstr='******* AMPLITUDE LOOP %d goodness stat %f ************' % (iloop,thisstat)
        loop3_service.loop3log (vis, pstr+'\n')
        if iloop!=ploop and thisstat-prevstat<0.01:
            pstr='****** EXITING AMPLITUDE CAL with diff %f *********'%(thisstat-prevstat)
            loop3_service.loop3log (vis, pstr+'\n')
            break
        else:   
            prevstat = thisstat
            loop3_service.imagr(visA,dopredict=True,fitsmask=fitsmask,autothreshold=3,dolocalrms=True,\
                                robust=-1,outname=visA+'_%02d%s'%(iloop,mfs))
        pstr='******* AMPLITUDE LOOP %d making new cal file %s ************' % (iloop,visA+'_%02d'%iloop)
        loop3_service.loop3log (vis, pstr+'\n')
        caltype, sol0 = strategy[iloop][0], float(strategy[iloop][1:])
        coh,cohlength = selfcal(visA,model='MODEL',incol='DATA',outcol='CORRECTED_DATA',\
                      outcal_root=visA+'_%02d'%iloop,caltype=caltype,init_sol=sol0)
        pstr='******** END AMPLITUDE LOOP %d - coherence on %.1f km **********' % \
                      (iloop,cohlength/1000.)
        loop3_service.loop3log (vis, pstr+'\n')


    fitsmask = init_fitsmask if iloop==ploop else visA+'_%02d-mask.fits'%(iloop-1)
    loop3_service.imagr(visA,cellsize='0.05asec',domultiscale=True,\
          outname=visA+'_final',channelsout=wsclean_chans,robust=-1,\
          fitsmask=fitsmask,dolocalrms=True)
    
    loop3_service.montage_plot (vis)
    pngfile, h5files = cleanup (vis)
    # loop3_service.loop3log (vis,'Output png file %s'%pngfile)  # commented out as the log file is already moved so this creates a new one with just one line
    print 'Output calibration tables',h5files
    return pngfile,h5files

#hybridloops(vis='/data020/scratch/sean/fafa/gaga/initsols_combined.ms/')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('vis', type=str, help='measurement set')
    parser.add_argument('--strategy', default='P30,P30,P30,A500,A450,A400', type=str, help='strategy for loops (default:P30,P30,P30,A500,A450,A400)' )
    parser.add_argument('--startmod', default='', type=str, help='starting model')
    parser.add_argument('--ith', default=5.0, type=float, help='threshold for pybdsf island detection, default 5.0')
    parser.add_argument('--bandwidth', default='8MHz', type=str, help='max bandwidth before breaking imaging into channels, default 8MHz')


    args = parser.parse_args()

    main( vis=args.vis, strategy=args.strategy, ith=args.ith, bandwidth=args.bandwidth )

