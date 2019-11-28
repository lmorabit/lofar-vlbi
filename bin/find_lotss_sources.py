# Try to work out what LoTSS sources to use in the calibration of a
# long-baseline field

# We use the LBCS catalogue to find all the sources within 1.5 degrees
# of the pointing centre. Then we need to look for potential LoTSS
# calibrator sources away from those objects.

import matplotlib.pyplot as plt
from download_lbcs_catalogue import main as download_lbcs
from download_lbcs_catalogue import grab_coo_MS as getmspos
from download_lotss_catalogue import main as download_lotss
import numpy as np
from astropy.table import Table,vstack

def separation(c_ra,c_dec,ra,dec):
    # all values in degrees
    return np.sqrt((np.cos(c_dec*np.pi/180.0)*(ra-c_ra))**2.0+(dec-c_dec)**2.0)

def plotcircle(ra,dec,radius,color='black'):
    theta=np.linspace(0,2*np.pi,1000)
    ra_c=ra+(radius/np.cos(dec*np.pi/180))*np.sin(theta)
    dec_c=dec+radius*np.cos(theta)
    plt.plot(ra_c,dec_c,color=color)
    

def locate_lotss_calibrators(ms,lotss_cat=None,outfile=None):
    radius=1.5
    exclude_radius=0.5
    figsize=8
    ra,dec=getmspos(ms)
    download_lbcs(ms,'lbcs.cat',Radius=radius)
    lines=open('lbcs.cat').readlines()
    lbcs_ra=[]
    lbcs_dec=[]
    for l in lines[1:]:
        bits=l.split(',')
        lbcs_ra.append(float(bits[0]))
        lbcs_dec.append(float(bits[1]))
    print ra,dec
    plt.figure(figsize=(figsize,figsize))
    plt.scatter(ra,dec,color='red',marker='+',s=200)
    plotcircle(ra,dec,radius)
    plt.scatter(lbcs_ra,lbcs_dec)
    for r,d in zip(lbcs_ra,lbcs_dec):
        plotcircle(r,d,exclude_radius,color='blue')

    # now look at LoTSS sources
    if lotss_cat is not None:
        print 'Reading user-supplied LOTSS catalogue',lotss_cat
        # use a user-supplied catalogue
        t=Table.read(lotss_cat)
        # cut to the radius we use
        s=separation(ra,dec,t['RA'],t['DEC'])
        t=t[s<radius]
        t['Resolved']=np.where(t['DC_Maj']<(10/3600.0),'U','R')
        t['Total_flux']*=1000
        t['Peak_flux']*=1000
    else:
        # we don't have user-supplied catalogue, query VO
        download_lotss(ms,'lotss.cat',Radius=radius,AllFile='lotss.csv')
        t=Table.read('lotss.csv')

    # exclude all sources already covered by LBCS
    filter=np.array([True]*len(t))
    for r,d in zip(lbcs_ra,lbcs_dec):
        s=separation(r,d,t['RA'],t['DEC'])
        filter&=(s>0.9*exclude_radius)
    t=t[filter]
    
    filter=t['Peak_flux']>100
    #filter&=t['Resolved']=='U'
    t_bright=t[filter]
    t_bright['Status']='Primary'
    print 'Filtered to',len(t_bright),'sources'
    plt.scatter(t_bright['RA'],t_bright['DEC'],color='orange')
    filter=t['Peak_flux']<100
    filter&=t['Peak_flux']>50
    filter&=t['Resolved']=='U'
    t_faint=t[filter]
    t_faint['Status']='Secondary'
    print 'Plus',len(t_faint),'fainter compact sources'
    plt.scatter(t_faint['RA'],t_faint['DEC'],color='yellow')
    
    plt.gca().invert_xaxis()
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    if outfile is not None:
        plt.savefig(outfile)
    else:
        plt.show()
    return vstack((t_bright,t_faint))

if __name__=='__main__':
    t=locate_lotss_calibrators('ILTJ133409.4+550149.0_L401323_SB244_uv.dppp_12723A7FEt_132MHz.msdpppconcat',lotss_cat='P205+55.mosaic.cat.fits',outfile='lotss.pdf')
    t.write('caltable.fits',overwrite=True)
    
    
