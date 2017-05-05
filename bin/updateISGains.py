#!/usr/bin/env python

import sys
import numpy
import subprocess
import sys
import textwrap
import cmath
from lofar.parmdb import parmdb

class ComplexArray(object):
    def __init__(self):
        raise NotImplementedError

    def get_amp(self):
        return numpy.absolute(numpy.nan_to_num(self.data))
    def set_amp(self, new_amps):
        self.data = numpy.array(new_amps) * self.data / numpy.absolute(self.data)
    amp = property(get_amp, set_amp)
    Ampl = amp

    def get_phase(self):
        return numpy.angle(numpy.nan_to_num(self.data))
    def set_phase(self, new_phase):
        self.data = numpy.vectorize(cmath.rect)(
            numpy.absolute(self.data), numpy.array(new_phase)
        )
    phase = property(get_phase, set_phase)
    Phase = phase

    def get_real(self):
        return numpy.real(numpy.nan_to_num(self.data))
    def set_real(self, new_real):
        self.data = numpy.array(new_real) + 1j * numpy.imag(self.data)
    real = property(get_real, set_real)
    Real = real

    def get_imag(self):
        return numpy.imag(numpy.nan_to_num(self.data))
    def set_imag(self, new_imag):
        self.data = numpy.real(self.data) + 1j *  numpy.array(new_imag)
    imag = property(get_imag, set_imag)
    Imag = imag

    @property
    def writeable(self):
        return dict((key, getattr(self, key)) for key in self.keys)

class RealImagArray(ComplexArray):
    keys = ("Real", "Imag")
    def __init__(self, real, imag):
        self.data = numpy.array(real) + 1j * numpy.array(imag)

class AmplPhaseArray(ComplexArray):
    keys = ("Ampl", "Phase")
    def __init__(self, ampl, phase):
        self.data = numpy.vectorize(cmath.rect)(
            numpy.array(ampl), numpy.array(phase)
        )

class PhaseOnlyArray(ComplexArray):
    keys = ("Phase", )
    def __init__(self, phase):
        self._phase = phase

    def get_zeros(self):
        return numpy.zeros(self.phase.shape)
    amp = property(get_zeros)
    real = property(get_zeros)
    imag = property(get_zeros)

    Ampl, Real, Imag = amp, real, imag

    def get_phase(self):
        return self._phase
    def set_phase(self, new_phase):
        self._phase = new_phase
    phase = property(get_phase, set_phase)
    Phase = phase

ARRAY_TYPES = [RealImagArray, AmplPhaseArray, PhaseOnlyArray]

class WriteableParmDB(parmdb):
    def __init__(self, name):
        super(WriteableParmDB, self).__init__(name)
        self.pdbname = name

    def setValues(self, name, values, start_freq, freqstep, start_time, timestep):
        """
        Write values to the ParmDB.

        Note that values should be a two dimenstional array with the first
        index corresponding to time and the second to time (this is the same
        as returned by ParmDB.getValues()).

        Arguments:

        name       -- Parameter name to write.
        values     -- NumPy array of values to write.
        start_freq -- Frequency at centre of first bin (Hz).
        freqstep   -- Bin-to-bin frequency increment (Hz).
        start_time -- Time at centre of first bin (MJD in seconds).
        timestep   -- Bin-to-bin time increment (s).
        """
        # This is the sequence of commands passed to parmdbm.
        template = '''
            open table="%(pdbname)s"
            remove %(name)s
            add %(name)s type="scalar", nx=%(freq_steps)s, ny=%(time_steps)s, values=%(values)s, domain=[%(start_freq)s, %(end_freq)s, %(start_time)s, %(end_time)s]
            quit
        '''
        template = textwrap.dedent(template).strip()
        time_steps, freq_steps = values.shape

        # For simplicity of user interface, we take the start values of the
        # axes at the centre of the bin, but we actually need to pass the
        # edges of the bin to parmdbm.
        start_time = start_time - timestep/2
        start_freq = start_freq - freqstep/2

        # Substitute appropriate values into the parmdbm command.
        command = template % {
            "pdbname": self.pdbname,
            "name": name,
            "time_steps": time_steps,
            "freq_steps": freq_steps,
            "values": str(list(values.ravel())),
            "start_freq": start_freq,
            "end_freq": start_freq+freqstep*freq_steps,
            "start_time": start_time,
            "end_time": start_time+timestep*time_steps
        }

        # Execute parmdbm and return the result.
        p = subprocess.Popen(['parmdbm'], stdin=subprocess.PIPE)
        p.communicate(command)
        p.wait()
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, cmd)
        else:
            return p.returncode

def list_stations(pdbfile):
    """
    Returns a list of all stations in the parmdb.
    """
    try:
        pdb = WriteableParmDB(pdbfile)
        return sorted(set(name.split(":")[-1] for name in pdb.getNames()))
    finally:
        pdb = None

class StationGain(dict):
    def __init__(self, parmdbfile, station):
        self.station = station
        self.parmdbfile = parmdbfile

        pdb = WriteableParmDB(self.parmdbfile)

        names = pdb.getNames("Gain:*:*:*:%s" % (station,))
        pols = set(":".join(x[1:3]) for x in  (x.split(":") for x in names))
        types = set(x[3] for x in  (x.split(":") for x in names))

        for array_type in ARRAY_TYPES:
            if sorted(types) == sorted(array_type.keys):
                self.array_type = array_type
                break
        assert(hasattr(self, "array_type"))

        for polarization in pols:
            data = []
            for key in self.array_type.keys:
                query = "Gain:%s:%s:%s" % (polarization, key, station)
                data.append(pdb.getValuesGrid(query)[query])

            self.timescale = data[0]['times']
            self.timestep = data[0]['timewidths'][0]

            self.freqscale = data[0]['freqs']
            self.freqstep = data[0]['freqwidths'][0]

            if len(self.array_type.keys) == 2:
                self[polarization] = self.array_type(data[0]["values"], data[1]["values"])
            elif len(self.array_type.keys) == 1:
                self[polarization] = self.array_type(data[0]["values"])
            else:
                raise Exception("Unknown number of keys")
        pdb = None

    def writeout(self):
        pdb = WriteableParmDB(self.parmdbfile)
        for pol, data in self.iteritems():
            for component, value in data.writeable.iteritems():
                pdb.setValues(
                    "Gain:%s:%s:%s" % (pol, component, self.station),
                    value,
                    self.freqscale[0],
                    self.freqstep,
                    self.timescale[0],
                    self.timestep
                )
        pdb = None
        
def main(parmdbfile):

    stations = list_stations(parmdbfile)
    
    amps00 = numpy.zeros(0)
    amps11 = numpy.zeros(0)
    arrayshape = (1)
    for station in stations:
        if 'CS' in station:
            a00 = StationGain(parmdbfile, station)['0:0'].amp
            a11 = StationGain(parmdbfile, station)['1:1'].amp
            arrayshape = a00.shape
            if numpy.count_nonzero(a00) > 0:
                amps00 = numpy.append(amps00,a00)
            if numpy.count_nonzero(a11) > 0:
                amps11 = numpy.append(amps11,a11)
    
    intl00 = numpy.median(amps00)*3.75
    intl11 = numpy.median(amps11)*3.75
    print intl00, intl11
    
    for station in stations:
        if not 'CS' in station and not 'RS' in station:
            print "Updating: " + station
            sg = StationGain(parmdbfile, station)
            if numpy.count_nonzero(sg['0:0'].amp) == 0:
                sg['0:0'].real = numpy.ones(arrayshape)
            sg['0:0'].amp = sg['0:0'].amp*intl00
            if numpy.count_nonzero(sg['1:1'].amp) == 0:
                sg['1:1'].real = numpy.ones(arrayshape)
            sg['1:1'].amp = sg['1:1'].amp*intl11
            sg.writeout()


if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print "Usage: %s <parmdbfile>" % sys.argv[0]
        sys.exit()
        
    main(sys.argv[1])
    
