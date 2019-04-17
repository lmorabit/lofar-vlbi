import numpy as np
import os
from losoto.h5parm import h5parm
from losoto.lib_operations import *
import pyrap.tables as pt
import logging
import argparse

h5parmfile="L1327+5504_phaseonly_ct.h5"
solsetname="sol000"

def readargs():
    parser=argparse.ArgumentParser("Apply the flags of an h5parm solset to its solutions, and set the respective weights to default.")
    parser.add_argument("--h5parm",type=str,help="Name of the h5parm(s) you want to apply weights for.",required=True,nargs="+")
    parser.add_argument("--soltabs",type=str,help="List of parameters you want to apply the flags on. "+\
                        "Default is [\'phase000\', \'tec\', \'phase_offset\', \'clock\'];"+\
                        " phase_offset will actually be skipped.",required=False,nargs="+",default=['phase000', 'tec', 'phase_offset', 'clock'])
    parser.add_argument("--solset",type=str,help="Solset in which your soltabs live. Default is sol000.",default="sol000",required=False)
    parser.add_argument("--suffix",type=str,help="Suffix you want to add for your flagged solsets. Default is _nice.",required=False,default="_nice")
    args=parser.parse_args()
    return vars(args)


def FixEverythingForever(h5parmf,soltabs=['phase000', 'tec', 'phase_offset', 'clock'],solset="sol000",suffix="_nice"):
    data=h5parm(h5parmf,readonly=False)
    if solset not in data.getSolsetNames():
        solset1=solset
        solset=data.getSolSetNames()[0]
        print "Requested solset %s not in h5parm; reading %s instead"%(solset1,solset)
    calsols=data.getSolset(solset)

    if "phase000" not in calsols.getSoltabNames():
        print "need phases, they're not here. Exiting."
        exit        
    else:
        for tab in soltabs:
            soltab=calsols.getSoltab(tab)
            if soltab.name=="phase_offset":
                print "Skipping phase_offset soltab. Curse h5parms"
                continue
            
            soltab=calsols.getSoltab(tab)
            # do the reshape
            soltab_axes = soltab.getAxesNames()
            
            new_times=soltab.time
            if len(soltab_axes) == 2:
                out_axes = ['time','ant']
                out_axes_vals = [new_times, soltab.ant]
                out_lens = (len(new_times),len(soltab.ant))
            elif len(soltab_axes) == 4:
                out_axes = ['time','ant','freq','pol']
                out_axes_vals = [new_times, soltab.ant, soltab.freq, soltab.pol]
                out_lens = (len(new_times),len(soltab.ant),len(soltab.freq),len(soltab.pol))
            elif len(soltab_axes) == 5:
                out_axes = ['time','freq','ant','dir','pol']
                out_axes_vals = [new_times, soltab.freq, soltab.ant, soltab.dir, soltab.pol]
                out_lens = (len(new_times),len(soltab.freq),len(soltab.ant),len(soltab.dir),len(soltab.pol))
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=soltab_axes, weight=True ):
                vals = reorderAxes ( vals, soltab_axes, out_axes )
                weights = reorderAxes( weights, soltab_axes, out_axes )

            # finally, can filter and reweight. Where weights are 0 set to nonzero average, and set associated values to 0.
            mask=(weights==False)
            vals[mask]=0
            weights[mask]=np.mean(weights[mask==False])
            # now write the thing with the suffix required
            newtabname=tab+suffix
            if newtabname in calsols.getSoltabNames():
                arse=calsols.getSoltab(newtabname)
                arse.delete()
            calsols.makeSoltab(soltype=soltab.getType(), soltabName=tab+suffix,
                               axesNames=out_axes, axesVals=out_axes_vals, vals=vals, weights=weights)
    
        print calsols.getSoltabNames()
            
                
if __name__=="__main__":
    args=readargs()
    filenames=args["h5parm"]
    soltabs=args["soltabs"]
    solset=args["solset"]
    suffix=args["suffix"]
    for h5 in filenames: 
        FixEverythingForever(h5,soltabs,solset,suffix)
