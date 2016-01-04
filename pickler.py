#!/usr/bin/env python

# Charles McEachern

# Fall 2015

# Note: This document wraps at column 80. 

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This script reads in Fortran-style data (which is pretty slow), then pickles
# it so that the plotter can have fast access. 

# #############################################################################
# #################################################### Importing Python Modules
# #############################################################################

# Explicit compression doesn't seem to do anything. The pickling protocol
# already delivers as much compression as is available. 
#import gzip

import gc
from sys import argv
import numpy as np
import os
# The cPickle module is faster, but not always available. 
try:
  import cPickle as pickle
except ImportError:
  import pickle
from time import time

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Start here. 
  if len(argv)>1:
    path = argv[1] 
  else:
    path = '/export/scratch/users/mceachern/2015nov23'

  print 'Pickling starting at ', path

  # Starting at the given path, walk through all (dir, subdir, files). Keep
  # directories containing a file 'params.in'.  
  dirs = [ d for d, s, f in os.walk(path) if 'params.in' in f ]

  # It's only worth pickling big files -- time-resolved bulk fields. 
  filenames = ('Bx.out', 'By.out', 'Bz.out', 'Ex.out', 'Ey.out', 'Ez.out', 
               'ExBy.out', 'EyBx.out', 'jz.out')

  # Pickle... everything! 
  for d in dirs:

    # We're going to be reading in TONS of data. Python's automatic garbage
    # collection might not be good enough to handle that, so we periodically
    # force collection. 
    print ''
    print 'Collecting garbage... ', gc.collect()
    print ''

    for f in filenames:

      # -----------------------------------------------------------------------
      # ------------------------------------------------------------ Read Array
      # -----------------------------------------------------------------------

      # Always use absolute paths (except colloquially). 
      filepath = d + '/' + f
      print filepath.replace(path, '').lstrip('/')
      # Make sure the file exists before trying to access it. 
      if not os.path.isfile(filepath):
        print '\tFile not found.'
        continue
      # Read in the array. 
      t0 = time()
      print '\tReading... ',
      arr = readArray(filepath)
      print format(time() - t0, '5.1f') + 's'

      # -----------------------------------------------------------------------
      # ---------------------------------------------------------- Pickle Array
      # -----------------------------------------------------------------------

      picklepath = d + '/' + f[:f.rfind('.')] + '.pkl'
      # Pickle the array. 
      t0 = time()
      with open(picklepath, 'wb') as handle:
        pickle.dump(arr, handle, protocol=-1)
      print '\tPickling... ' + format(time() - t0, '5.1f') + 's'
      # Report the before and after file sizes. 
      print '\t' + str( int( 1e-6*os.path.getsize(filepath) ) ) + ' MB - > ',
      if not os.path.isfile(picklepath):
        print '??? MB'
        continue
      else:
        print str( int( 1e-6*os.path.getsize(picklepath) ) ) + ' MB'
      # Check load time. 
      t0 = time()
      with open(picklepath, 'rb') as handle:
        pklArr = pickle.load(handle)
      print '\tUnpickling... ' + format(time() - t0, '5.1f') + 's'
      # Check load accuracy. 
      if np.all(arr==pklArr):
        print '\tChecking agreement... OK'
      else:
        print '\tChecking agreement... X'

#      # -----------------------------------------------------------------------
#      # ------------------------------------------------------ Zip-Pickle Array
#      # -----------------------------------------------------------------------
#
#      zippath = d + '/' + f[:f.rfind('.')] + '.pkz'
#      # Also try pickling with gzip for added compression. 
#      t0 = time()
#      handle = gzip.open(zippath, 'wb')
#      pickle.dump(arr, handle, protocol=-1)
#      handle.close()
#      print '\tZip-pickling... ' + format(time() - t0, '5.1f') + 's'
#      # Report the before and after file sizes. 
#      print '\t' + str( int( 1e-6*os.path.getsize(filepath) ) ) + ' MB - > ',
#      if not os.path.isfile(zippath):
#        print '??? MB'
#        continue
#      else:
#        print str( int( 1e-6*os.path.getsize(zippath) ) ) + ' MB'
#      # Check load time. 
#      t0 = time()
#      handle = gzip.open(zippath, 'rb')
#      pkzArr = pickle.load(handle)
#      handle.close()
#      print '\tUn-zip-pickling... ' + format(time() - t0, '5.1f') + 's'
#      # Check load accuracy. 
#      if np.all(arr==pkzArr):
#        print '\tChecking agreement... OK'
#      else:
#        print '\tChecking agreement... X'

  return

# #############################################################################
# ############################################################### Read in Array
# #############################################################################

# Turns a list of numbers (1, 2, 3) into the string '1x2x3'. 
def by(x):
  return str( x[0] ) + 'x' + by( x[1:] ) if len(x)>1 else str( x[0] )

# Convert a string to a complex or float. 
def com(x):
  if ',' not in x:
    return float(x)
  else:
    # Shave off the parentheses then split into real and imaginary parts. 
    re, im = x[1:-1].split(',')
    return (float(re) + float(im)*1j)

# Note that we have already checked to make sure that the file exists. 
def readArray(filename):
  # Read in the file as a list of strings. 
  arrayLines = open(filename, 'r').readlines()
  # The first line is the array dimensions. 
  dims = [ int(x) for x in arrayLines.pop(0).split() ]
  # Assemble a one-dimensional array large enough to hold all of the values.
  # (This is much faster than appending as we go.) This means figuring out if
  # we want reals or complexes. 
  firstValue = com( arrayLines[0].split()[0] )
  dtype = np.complex if isinstance(firstValue, np.complex) else np.float
  nVals = np.prod(dims)
  vals = np.empty(nVals, dtype=dtype)
  # Now fill the array with values one at a time. Stop when it's full, or when
  # we run out of values. 
  i = 0
  for line in arrayLines:
    for val in line.split():
      # Check if the array is full. 
      if i==nVals:
        break
      # If it's not, grab the next value. 
      else:
        vals[i] = com(val)
        i = i + 1
  # Before returning the array, reshape and transpose it. Fortran and Python
  # have opposite conventions for which index should change fastest. 
  arr = np.transpose( np.reshape( vals, dims[::-1] ) )
  # Check if the array is full before returning it. If it's not, the run may
  # have crashed or run out of time. Return as many time steps as possible. 
  if i==nVals:
    return arr
  else:
    actualDims = dims[:-1] + [ np.int( i/np.prod(dims[:-1]) ) ]
    print 'WARNING: Found ' + by(actualDims) + ' not ' + by(dims) + '... ',
    return arr[..., :actualDims[-1] ]

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()

