import argparse
from astropy.table import Table, Column
import numpy as np

def is_resolved(Sint, Speak, rms):
    """ Determines if a source is resolved or unresolved.
    The calculation is presented in Shimwell et al. 2018 of the LOFAR DR1 paper splash.
    
    Args:
        Sint (float or ndarray): integrated flux density.
        Speak (float or ndarray): peak flux density.
        rms (float or ndarray): local rms around the source.
    Returns:
        resolved (bool or ndarray): True if the source is resolved, False if not.
    """
    resolved = ((Sint / Speak) > 1.25 + 3.1 * (Speak / rms) ** (-0.53))
    return resolved

def main( catalogue ):

    ## read the catalogue
    mycat = Table.read( catalogue, format='csv' )

    resolved = np.where( is_resolved( mycat['Total_flux'], mycat['Peak_flux'], mycat['Isl_rms'] ), 'R', 'U' )
    res_column = Column( name='Resolved', data=resolved )
    mycat.add_column( res_column )
    mycat.write( catalogue.replace('.csv', '_rescol.csv' ), format='csv' )

   
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('csvfile', type=str, help='csv catalogue' )
    args = parser.parse_args()

    main( args.csvfile )
