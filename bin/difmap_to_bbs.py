import glob
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import angle_utilities
from astropy.coordinates.representation import UnitSphericalRepresentation
import casacore.tables as ct
import argparse
import os

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

def main( msfile ):

    tmp = msfile.split('/')[-1]
    filestem = tmp.split('_')[0]

    ## get reference frequency
    rf = ct.taql('select REF_FREQUENCY from {:s}::SPECTRAL_WINDOW'.format(msfile))
    ref_freq = rf.getcol('REF_FREQUENCY')[0]

    modfiles = glob.glob( filestem + '*.mod' )
    xx_file = [mf for mf in modfiles if 'XX' in mf ][0]
    yy_file = [mf for mf in modfiles if 'YY' in mf ][0]

    xx_mod = read_mod( xx_file )
    yy_mod = read_mod( yy_file )

    ## find central point
    xx_cen = find_centre( xx_mod )
    yy_cen = find_centre( yy_mod )

    xx_flux = []
    for mykey in xx_mod.keys():
        xx_flux.append( xx_mod[mykey][0] )

    yy_flux = []
    for mykey in yy_mod.keys():
        yy_flux.append( yy_mod[mykey][0] )


    outfile = filestem + '_StokesI.mod'
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Convert a difmap model into bbs format")
    parser.add_argument( "msfile", help="measurement set (to get reference frequency" )
   
    args = parser.parse_args()

    main( args.msfile )

