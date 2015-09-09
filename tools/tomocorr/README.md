# Tomographic correction

This is now callable as a Python function, but must be built using f2py first.
See 'building the Python interface' below for details.


## Calling the code in Python

	>>> import tomo_predict
	>>> dt = tomo_predict.tomo_predict.predict(1d_model_filename, 3d_model_filename,
		                                       lat, lon, depth)


## TODO

### TauPy doesn't give you actual coordinates
The current implementation of raypath calculation in ObsPy does not give you
ray paths which are in true geographic space; they are only given in distance
and radius.  However, a small routine to 'do the maths' would not be difficult
(for a spherical Earth; elliptical calculation is a bit tricker, but will
presumably be in some Python module out there somewhere).

One would need to be able to do:

	>>> (lon, lat, dep, distance) = path(TauPy_model, evtlon, evtlat, evtdep,
		                                 stalon, stalat)

### `tomocorr.py`

Rewrite `tomocorr.py` after implementing the above.

## Updated Fortran module

`tomo_predict.f90` is an updated version of the original
`tomo_predict_vdh.f` code, which is designed to be callable directly from
Python.  Before first using the code, call

	>>> tomo_predict.tomo_predict.setup(1d_model_file, 3d_model_file)

and then you can do

	>>> dt = tomo_predict.tomo_predict.tomo_delay(lat, lon, dep)

to get the relative predicted travel time anomaly for the arrays describing lon,
lat and depth in degrees and km.

Alternatively, you can just call

	>>> dt = tomo_predict.tomo_predict.predict(1d_model_file, 3d_model_file, lat, lon, dep)

each time, and setup will be done automatically on the first invocation.


## Building Fortran program `tomo_predict_vdh`

This functions in the same way as Ed Garnero's `tomo_predict.f` program,
except it ignores the `top` and `bot` input; though it does read it.

To build, type

	$ make progs

This has been tested and gives the same results as the original, though
it is arguably more accurate, as it does not remove the top 1 km of the
path at the receiver end (and source, if at the surface).


## Building the Python interface

Type

	$ make

or

	$ make libs

You may need to set your f2py, Fortran compiler and favourite options, e.g.:

	$ make libs F2PY=f2py FC="ifort" FFLAGS="-Qweird-intel-options"
