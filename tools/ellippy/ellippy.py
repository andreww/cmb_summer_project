#!/usr/bin/env python

import inspect
import os

import numpy as np

__DATA_DIR = os.path.dirname(os.path.abspath(inspect.getfile(
              inspect.currentframe())))

import _ellip_fort
import _direct_fort

def ellip_setup():

    # Nasty nasty nasty Fortran, with no concept of a path
    startdir = os.getcwd()
    os.chdir(__DATA_DIR)
    _direct_fort.direct()
    os.chdir(startdir)
    

def ellip_correct(src_lat, src_depth, bazim, delta, phase):
   """

      src_lat: source lattitude (in degrees)
      src_depth: source depth (in km)
      bazim: azimuth from source (in degrees)
      delta: epicentral distance (in degrees)
      
   """

   # Convert src_lat to co-lat
   co_lat = (90.0 - src_lat)

   # convert all angular values from degrees to rad
   co_lat = np.radians(co_lat)
   bazim = np.radians(bazim)
   #delta = np.radians(delta)

   # Phase must be 8 chars (at least)
   if len(phase) < 8:
      phase = phase + (' '*(8-len(phase)))

   # call fortran
   # Nasty nasty nasty Fortran, with no concept of a path
   startdir = os.getcwd()
   os.chdir(__DATA_DIR)
   _ellip_fort.ellref(co_lat)
   tcor, abrt = _ellip_fort.ellcor(phase, delta, src_depth,
        co_lat, bazim)
   os.chdir(startdir)

   # Handle "phase not found" errors
   if abrt:
       raise ValueError("Phase " + phase.strip() + " is not in the phase list")

   return tcor

def ellip_src_sta(src_lat, src_lon, src_depth, sta_lat, sta_lon, phase):

    # FIXME: implement this...

    raise NotImplementedError

if __name__ == "__main__":
    # TODO: add doc string and use command line parser
    import sys
    src_lat = float(sys.argv[1]) 
    src_depth = float(sys.argv[2]) 
    bazim = float(sys.argv[3]) 
    delta = float(sys.argv[4]) 
    phase = sys.argv[5]
    tcor = ellip_correct(src_lat, src_depth, bazim, delta, phase)
    print tcor
