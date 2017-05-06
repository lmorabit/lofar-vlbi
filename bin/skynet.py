import numpy as np,os,sys,glob,astropy
import closure,model_engine,taql_funcs

# Requires:
# 1) Measurement set
# 2) list of three closure telescopes (e.g. ['TS001','DE601','DE605'])

# For the measuement set, evaluates the closure-scatter. If this is <CTHR, makes a model

# modes: 
# 1) make a short-baseline image (phases already calibrated) and selfcal to that
# 2) start with a point source
# 3) start with something from the model-engine

# Returns: closure threshold value
#          Images will be left in directory depending on CTh and mode

# 4ch/8s 20Gb per source (FOV) 1ch/16s 3 GB/source to go to 5' fields
# eor scripts gives out cal table as a lofar parmdb

def skynet_NDPPP (vis,model,solint=1.0):
    os.system('rm -fr %s/sky\n'%vis)
    os.system ('makesourcedb in=%s_mod out=%s/sky format=\'<\''%(vis,vis))
    with open('NDPPP.parset','w') as f:
        f.write('msin=%s\n'%vis)
        f.write('msin.datacolumn=DATA\n')
        f.write('msout=%s\n'%vis)
        f.write('msout.datacolumn=CORRECTED_DATA\n')
        f.write('steps=[gaincal]\n')
        f.write('gaincal.applysolution=True\n')
        f.write('gaincal.sourcedb = %s/sky\n'%vis)
        f.write('gaincal.maxiter=200\n')
        f.write('gaincal.caltype=phaseonly\n')
        f.write('gaincal.solint=%d\n'%int(solint))
    f.close()
    os.system('NDPPP NDPPP.parset')
        

def skynet (vis, mode=3, closure_tels=['ST001','DE601','DE605'],cthr=1.3):
    ra,dec = taql_funcs.taql_from (vis, 'FIELD', 'PHASE_DIR')
    closure_scatter = closure.closure(vis, closure_tels, plotfile='')
    if closure_scatter > cthr:
        return closure_scatter
    if mode == 1:     # this is the default of the self_calibration_pipeline_V2.
        os.system('python self_calibration_pipeline_V2.py -m '+vis+' -p') # amplitudes?
    elif mode == 2:   # use NDPPP to selfcal the long baselines to a small (0.1") Gaussian
        model_engine.write_skymodel (ra,dec,np.array([0.0,0.0,1.0,0.1,0.0,0.0]),vis+'_mod')
        skynet_NDPPP (vis,vis+'_mod',solint=5)  # timesteps
        os.system('python self_calibration_pipeline_V2.py -d CORRECTED_DATA -m '+vis+' -p')
    elif mode == 3:   # make an engine model and selfcal against this
        model_engine.mainscript (vis,closure_tels,PLOTTYPE=0,outname=vis+'_mod')
        skynet_NDPPP (vis,vis+'_mod',solint=5)
        os.system('python self_calibration_pipeline_V2.py -d CORRECTED_DATA -m '+vis+' -p')
    else:
        return np.nan
#    return closure_scatter
