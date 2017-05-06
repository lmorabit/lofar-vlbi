#!/usr/bin/env python
import os
import numpy as np

## written by Leah Morabito May 2017


########################################################################

def main(closure_phase_file):

    """
    Find the calibrator with the best closure phase to use for the delay calibration

    Parameters
    ----------
    closure_phase_file : str
	String containg the name of the closure phase file as output by the pipeline

    Returns
    -------
    {'DelayCal':delayCal} : "dict"
	name of the best calibrator
   
    """

    ## read the file with the scatter
    with open( closure_phase_file ) as f:
	lines = f.readlines()
    f.close()

    ## get lists of directions and scatter
    direction = []
    scatter = []
    for l in lines:
	direction.append(l.split()[4])
	scatter.append(np.float(l.split()[6]))

    ## convert to numpy arrays
    direction = np.asarray( direction )
    scatter = np.asarray( scatter )

    ## find the minimum scatter
    if len(scatter) > 1:
        min_scatter_index = np.where( scatter == np.min( scatter ) )[0]
        best_calibrator = direction[min_scatter_index][0]	
    else:
        best_calibrator = direction[0][0]

    delayCal = os.getcwd() + '/' + best_calibrator

    return { 'DelayCal' : delayCal }

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Find best calibrator source')
    parser.add_argument('closurePhaseFile', type=str, help='Name of closure phase file.')
    args = parser.parse_args()
    main( args.closurePhaseFile )

