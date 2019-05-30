#!/usr/bin/env python
import os
import glob

def main( ms_input ):

    tec_h5parms = glob.glob('S*_tec.h5')

    for tec_h5parm in tec_h5parms:

        filestem = tec_h5parm.replace('.h5','')
        losoto_parset = filestem + '_losoto.parset' 
        imstem = filestem.split('_')[0] + '_' + filestem.split('_')[1] + '_' + filestem.split('_')[2] + '_'
        with open( losoto_parset, 'w' ) as f:
            f.write('[plotTEC]\n')
            f.write('plotFlag = False\n')
            f.write('axesInPlot = [time]\n')
            f.write('prefix = %s\n'%imstem)
            f.write('soltab = sol000/tec000\n')
            f.write('axisInTable = ant\n')
            f.write('operation = PLOT\n')
            f.write('refAnt = ST001')
        f.close()

        os.system( 'losoto -V %s %s > %s 2>&1'%(tec_h5parm, losoto_parset, losoto_parset.replace('parset','log') ) )



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='plot tec solutions for all existing h5parms')
    parser.add_argument('MS_pattern',type=str, help='pattern to search for MS')
    args = parser.parse_args()

    MS_input = glob.glob( args.MS_pattern )

    main( MS_input )
