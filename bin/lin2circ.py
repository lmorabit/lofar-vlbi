#!/usr/bin/env python

"""
lin2circ.py
Converts linear XX,XY,YX,YY correlations to circular RR,RL,LR,LL
assuming that the definitions follow the IAU ones, which are
described by Hamaker & Bregman (1996, A&AS, 117, 161).

In particular, I use the coordinate transformation given in Section 3,

          1     /  1   +i  \ 
C_A =  -------  |          |
       sqrt(2)  \  1   -i  /


The Hermitian conjugate of this is,

           1     /   1   1  \ 
C+_A =  -------  |          |
        sqrt(2)  \  -i  +i  /

So V_RL = C_A * V_XY * C+_A

where V_XY is the visibilities in linear coordinates,
      V_RL is the visibilities in circular coordinates.

This reduces to:

RR = XX - iXY + iYX + YY
RL = XX + iXY + iYX - YY
LR = XX - iXY - iYX - YY
LL = XX + iXY - iYX + YY

Version 1.0 written 15 April 2010 by George Heald

Version 1.1 (18 March 2012) adds update of POLARIZATION table

Version 1.2 (July 6, 2012) - added inverse capability, convert from circular
to linear polarization.

Version 1.3 (2014, Reinout van Weeren) - corrected the flux error, factor of 0.5 (/2) was missing. Added option to convert back from circular to linear

Since the transformation matrices are hermitian, C_A ^-1 = C+_A and C+_A^-1 = C_A.
So, we have: 

	V_XY = C+_A*V_RL*C_A
	
which is:

XX =   RR +  RL +  LR +  LL
XY =  iRR - iRL + iLR - iLL
YX = -iRR - iRL + iLR + iLL
YY =   RR -  RL -  LR +  LL
"""

import optparse
import pyrap.tables as pt
import numpy

def main(options):

	cI = numpy.complex(0.,1.)
	
	inms = options.inms
	if inms == '':
			print 'Error: you have to specify an input MS, use -h for help'
			return
	column = options.column
	outcol = options.outcol
	
	t = pt.table(inms, readonly=False, ack=True)
	if options.back:
		lincol = options.lincol
		if lincol not in t.colnames():
				print 'Adding the output linear polarization column',lincol,'to',inms
				coldmi = t.getdminfo(column)
				coldmi['NAME'] = lincol
				t.addcols(pt.maketabdesc(pt.makearrcoldesc(lincol, 0., valuetype='complex', shape=numpy.array(t.getcell(column,0)).shape)), coldmi)

                ### RVW EDIT 2012   
		print 'Reading the input column (circular)', column
		if column not in t.colnames():
			print 'Error: Input column does not exist'
			return
		
		### RVW EDIT 2012 Input column with the -c switch 
		cirdata = t.getcol(column)
		#cirdata = t.getcol(outcol)
		#print 'SHAPE ARRAY', numpy.shape(cirdata)
		
		print 'Computing the linear polarization terms...'
		lindata = numpy.transpose(numpy.array([
				0.5*(cirdata[:,:,0]+cirdata[:,:,1]+cirdata[:,:,2]+cirdata[:,:,3]),
				0.5*(cI*cirdata[:,:,0]-cI*cirdata[:,:,1]+cI*cirdata[:,:,2]-cI*cirdata[:,:,3]),
				0.5*(-cI*cirdata[:,:,0]-cI*cirdata[:,:,1]+cI*cirdata[:,:,2]+cI*cirdata[:,:,3]),
				0.5*(cirdata[:,:,0]-cirdata[:,:,1]-cirdata[:,:,2]+cirdata[:,:,3])]),
				(1,2,0))
		print 'Finishing up...'
		t.putcol(lincol, lindata)
	else:
		if outcol not in t.colnames():
			print 'Adding the output column',outcol,'to',inms
			coldmi = t.getdminfo(column)
			coldmi['NAME'] = outcol
			t.addcols(pt.maketabdesc(pt.makearrcoldesc(outcol, 0., valuetype='complex', shape=numpy.array(t.getcell(column,0)).shape)), coldmi)
		print 'Reading the input column (linear)', column
		data = t.getcol(column)
		print 'Computing the output column'
		outdata = numpy.transpose(numpy.array([
				0.5*(data[:,:,0]-cI*data[:,:,1]+cI*data[:,:,2]+data[:,:,3]),
				0.5*(data[:,:,0]+cI*data[:,:,1]+cI*data[:,:,2]-data[:,:,3]),
				0.5*(data[:,:,0]-cI*data[:,:,1]-cI*data[:,:,2]-data[:,:,3]),
				0.5*(data[:,:,0]+cI*data[:,:,1]-cI*data[:,:,2]+data[:,:,3])]),
				(1,2,0))
		print 'Finishing up...'
		t.putcol(outcol, outdata)
	if options.poltable:
		print 'Updating the POLARIZATION table...'
		tp = pt.table(inms+'/POLARIZATION',readonly=False,ack=True)
		
		### RVW EDIT 2012
		if options.back:
		   tp.putcol('CORR_TYPE',numpy.array([[9,10,11,12]],dtype=numpy.int32)) # FROM CIRC-->LIN
                else:
		   tp.putcol('CORR_TYPE',numpy.array([[5,6,7,8]],dtype=numpy.int32)) # FROM LIN-->CIRC

opt = optparse.OptionParser()
opt.add_option('-i','--inms',help='Input MS [no default]',default='')
opt.add_option('-c','--column',help='Input column [default DATA]',default='DATA')
opt.add_option('-o','--outcol',help='Output column [default DATA_CIRC]',default='DATA_CIRC')
opt.add_option('-p','--poltable',help='Update POLARIZATION table? [default False]',default=False,action='store_true')
opt.add_option('-b','--back',help='Go back to linear polarization [default False]',default=False,action='store_true')
opt.add_option('-l','--lincol',help='Output linear polarization column, if the -b switch is used [default DATA_LIN]; we want to keep the original DATA column',default='DATA_LIN')
options, arguments = opt.parse_args()
main(options)

