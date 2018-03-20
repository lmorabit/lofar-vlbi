Code descriptions

find_lotss_sources.py: provides a function

    locate_lotss_calibrators(ms,lotss_cat,outfile)

  -- ms is the name of a reference ms (must be at the field center).
  -- lotss_cat if supplied is the name of a LoTSS PyBDSF
     catalogue. Otherwise the VO will be queried
  -- outfile if supplied is the name of an output file for the
     diagnostic plot, otherwise plt.show() is called

  The function returns an astropy table of potential calibrator
  sources and makes a diagnostic plot. Note that the
  download_lbcs_catalogue code needs to be on PYTHONPATH.

losoto_add_phases.py: provides a function

    combine_phases(reference, differential, output)

  -- reference is the name of the reference HDF5 file
  -- differential is the name of the differential file
  -- output is the name of the outfile HDF5 file (will be overwritten)

  The function adds the phases in differential to the ones in
  reference and writes the result to output.


