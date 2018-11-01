import numpy as np,sys,os,multiprocessing,glob
from AIPS import AIPS,AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from lofipi_aips import *

def zaptabs (inna,incl,indisk,inseq,tabtype,top,bottom):
    uv = AIPSUVData(inna,incl,indisk,inseq)
    for i in range (top,bottom-1,-1):
        try:
            uv.table(tabtype,i).zap()
        except:
            pass

def zapuv (inna, incl, indisk, top, bottom):
    for i in range(top,bottom-1,-1):
        try:
            uv = AIPSUVData(inna,incl,indisk,i).zap()
        except:
            pass

def zapim (inna, incl, indisk, top, bottom):
    for i in range(top,bottom-1,-1):
        try:
            im = AIPSImage(inna,incl,indisk,i).zap()
        except:
            pass
           

def fullcal (fringe_inna, fcl, indisk, refant):
    zaptabs (fringe_inna, fcl, indisk, 1, 'SN', 10, 1)     
    zaptabs (fringe_inna, fcl, indisk, 1, 'CL', 10, 2)     
    zapim (fringe_inna, 'ICL001', indisk, 10, 1)
    zapim (fringe_inna, 'IBM001', indisk, 10, 1)
#    source = AIPSUVData(fringe_inna,fcl,indisk,1).sources[0]
    source = 'BEAM_1'
    pfring (fringe_inna,refant,[0],source,solint=5,\
            weightit=1,zero=5,aipsclass=fcl,\
            suppress_rate=1)                        # delay cal -> SN1
    psnsmo (fringe_inna,fcl,indisk,'VLBI')          # edited delay -> SN2
    pclcal (fringe_inna,indisk,2,aipsclass=fcl)     # SN2 -> CL2
    pcalib (fringe_inna,fcl,indisk,0.5,refant,\
            docalib=1)                              # phase cal w/CL2 -> SN3
    pclcal (fringe_inna,indisk,3,aipsclass=fcl)     # SN3+CL2 ->CL3
    pimagr (fringe_inna,fcl,1)                      # image with CL3 -> ICL001.1
    pcalib (fringe_inna,fcl,indisk,0.5,refant,\
            solmode='A&P',docalib=1,in2name=fringe_inna,\
            in2class='ICL001',in2seq=1,calsour=source)  # CL3/ICL001.1 -> SN4
    pclcal (fringe_inna,indisk,4,aipsclass=fcl)     # SN4 + CL3 -> CL4
    pimagr (fringe_inna,fcl,1)                      # image with CL4 -> ICL001.2
    for i in [4,3,2,1]:
        AIPSUVData(fringe_inna,fcl,indisk,1).table('SN',i).zap()
    for i in [4,3,2]:
        AIPSUVData(fringe_inna,fcl,indisk,1).table('CL',i).zap()
    pfring (fringe_inna,refant,[0],source,solint=5,\
            weightit=1,zero=5,aipsclass=fcl,\
            suppress_rate=1,in2name=fringe_inna,\
            in2class='ICL001',in2seq=2,docalib=-1)  # delay cal w/ICL001.2 -> SN1
    psnsmo (fringe_inna,fcl,indisk,'VLBI')          # edited delay ->SN2
    pclcal (fringe_inna,indisk,2,aipsclass=fcl)     # SN2 -> CL2
    pcalib (fringe_inna,fcl,indisk,0.5,refant,\
            docalib=1,in2name=fringe_inna,\
            in2class='ICL001',in2seq=2)             # phase cal w/CL2,ICL001.2 -> SN3
    pclcal (fringe_inna,indisk,3,aipsclass=fcl)     # SN3+CL2 -> CL3
    pimagr (fringe_inna,fcl,1)                      # image w/CL3 -> ICL001.3
    source = AIPSUVData(fringe_inna,fcl,indisk,1).sources[0]
    pcalib (fringe_inna,fcl,indisk,0.5,refant,\
            solmode='A&P',docalib=1,\
            in2name=fringe_inna, in2class='ICL001',
            in2seq=3,calsour=source)                # CL3/ICL001.3 -> SN4
    pclcal (fringe_inna,indisk,4,aipsclass=fcl)     # SN4 + CL3 -> CL4
    pimagr (fringe_inna,fcl,1)                      # image with CL4 -> ICL001.4
    
