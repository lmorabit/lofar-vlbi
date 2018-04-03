#!/bin/python
from astropy.io import fits
import numpy as np
import sys

def main(fitsim):

    ## read the data and header
    data, header = fits.getdata(fitsim, header=True)

    ## get the keys
    header_keys = header.keys()

    ## check if there is already the BMAJ keyword
    bmaj_index = [ i for i, val in enumerate(header_keys) if val == 'BMAJ' ]
    ## if not, get it from the history and write it as a key
    if len( bmaj_index == 0 ):
        history_vals = header['HISTORY']
        beam_line = [ hval for hval in history_vals if 'BMAJ' in hval ] 
        tmp = beam_line[-1]
        tmp1 = tmp.split('BMAJ')[-1].lstrip('=').lstrip(' ').split(' ')
        bmaj = np.float(tmp1[0])
        tmp2 = tmp.split('BMIN')[-1].lstrip('=').lstrip(' ').split(' ')
        bmin = np.float(tmp2[0])
        tmp3 = tmp.split('BPA')[-1].lstrip('=').lstrip(' ').split(' ')
        bpa = np.float(tmp3[0])

        header['BMAJ'] = bmaj
        header['BMIN'] = bmin
        header['BPA'] = bpa
        fits.writeto(fitsim, data, header, overwrite=True)


if __name__ == "__main__":
    fitsim = sys.argv[1]
    main(fitsim)
