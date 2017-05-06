# From a MS, return the index numbers of telescopes. Neal Jackson 03/17
# Input telescopes given as an array of strings
# Note that the index numbers are found in one of two places, the NAME and
# STATION columns (latter if it has been written by AIPS, for instance), so
# this routine deals with that too. NB returns the OFFSET (starts at 0) not
# the AIPS telescope number (starts at 1) in this case.
# Output is a list of indices of telescopes. Note that you need to use the
# lower one first if reading data out of a MS.
import os,numpy as np
def get_idx_tels (data, tel):
    os.system('taql \'select NAME from %s/ANTENNA\' >closure_which'%data)
    f = open('closure_which')
    for line in f:
        if 'select' in line or not len(line.rstrip('\n')):
            continue
        try:
            a = int(line) # it's just a number, use STATION column instead
            f.close()
            os.system('taql \'select STATION from %s/ANTENNA\' >closure_which'%data)
            break
        except:
            f.close()
            break
    idx_tels, iline = [-1]*len(tel), 0
    f = open('closure_which')
    for line in f:
        if not 'select' in line:
            for i in range(len(tel)):
                if tel[i]==line[:len(tel[i])]:
                    idx_tels[i] = iline
            iline += 1
    f.close()
    if -1 in idx_tels:
        os.system('cat closure_which')
        print 'Did not find one or more of the telescopes'
        print 'Telescopes present are those in list above'
        return []
    os.system('rm closure_which')
    return idx_tels

