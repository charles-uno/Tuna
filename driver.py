#!/usr/bin/env python

# Charles McEachern

# Fall 2015

# Note: This document wraps at column 80. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

from math import ceil
import os
from shutil import copy, copytree, rmtree
from socket import gethostname as host
from subprocess import Popen, PIPE
from sys import argv
from time import localtime as lt, time

# #############################################################################
# ############################################################# Parameter Entry
# #############################################################################

# If we want to give this run a custom directory name, put that here. 
runDirName = 'tuna'

# Tuna includes default values, which it uses for any parameter not specified. 
parameters = {
              'jdrive':[4e-4], 
#              'bdrive':[10], 
              'tmax':[500],
              'azm':[1, 4, 16, 64],
              'model':[1, 2, 3, 4],
              'fdrive':[0.010, 0.016, 0.022],
              'inertia':[-1]
         }

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():
  # If on Itasca, this grabs a fresh version of the current directory. Then it
  # cleans up any temporary files. 
  print refresh() + '\n'
  # Prepare run directory. 
  global runDirName
  rd = setRunDir(runDirName)
  print rd + '\n'
  # Permute the parameters to get a list of runs to conduct. 
  runs = getRuns(rd)
  # Prepare programming environment. 
  print setEnvironment()
  # Compile binary executable. If the compilation fails, bail. 
  if not setBin():
    print '\tCompilation failed. '
    print '\t' + rd + 'stdoe.txt'
    return
  # Iterate over the runs. 
  for r in runs:
    runRun(r)

  # We have a plotter to look at plots. Get rid of plot stuff from the driver, 
  # so we can call it over an SSH connection without graphics forwarding. 
#  # Let's just look at the Alfven speed profile. 
#  return vPlot()
  # Plot T000 for debugging. 
#  ax = initPlot()
#  showPlot(ax)

  return

# #############################################################################
# ############################################################## Setup for Runs
# #############################################################################

# =============================================================================
# ========================================================= Refresh Source Code
# =============================================================================

# If we're on Itasca, grab a fresh copy of the source code from Eelpout before
# the run. 
def refresh():
  # No refreshing is necessary if we're already on the machine with the most
  # recent version of the code. 
  if local():
    return 'No refresh needed. '
  # Nuke the current directory. 
  for x in os.listdir('.'):
    rm(x)
  # Grab a fresh version of everything from Eelpout.       
  src = 'mceachern@eelpout.spa.umn.edu:~/Desktop/tuna/'
  dest = '/home/lysakrl/mceache1/tuna/'
  # The status bar from scp doesn't get captured nicely. We can worry about
  # that later... it's probably involve the -v flag. 
  bash('scp -r ' + src + '/*.* ' + src + '/models .', out=None, err=None)
  # Clean up any temporary files in the current directory. 
  [ os.remove(f) for f in files() if f.endswith('~') ]
  # Check if this is all we're supposed to do. 
  if 'refresh' in argv:
    print 'Refreshing current directory, then exiting. '
    exit()
  else:
    return 'Refreshing current directory'

# =============================================================================
# ======================================================== Create Run Directory
# =============================================================================

def setRunDir(rdn='runs'):
  sd, rd = srcDir(), runDir(rdn)
  # Create a timestamped run directory. 
  os.mkdir(rd)
  # Copy the source files and ionospheric models to the run directory. 
  copy(sd + 'source.f90', rd)
  copytree(sd + 'models', rd + 'models')
  # Move to the run directory. 
  os.chdir(rd)
  return rd

# =============================================================================
# ========================================================= Create List of Runs
# =============================================================================

# We scramble the given parameters into all possible permutations. Each 
# permutation corresponds to a run. 
def getRuns(rd):
  global parameters
  # For nonzero field driving we need 0 current driving and vice versa. 
  if 'bdrive' in parameters and 'jdrive' in parameters:
    parameters['bdrive'] = ulist( parameters['bdrive'] + [0] )
    parameters['jdrive'] = ulist( parameters['jdrive'] + [0] )
  # Initialize runs object to be bifurcated. 
  runs = [ {} ]
  # for each parameter we're looking at...
  for key in parameters:
    # Bifurcate the existing runs for each value of that parameter.
    runs = [ dsum(r, {key:val}) for r in runs for val in parameters[key] ]
  # Cut runs that have no driving or double driving. 
  runs = [ r for r in runs if nonzero(r, 'bdrive') != nonzero(r, 'jdrive') ]
  # Cut runs with inertial effects but no Boris factor. They will take months. 
  for r in runs[:]:
    if 'inertia' in r and r['inertia']>0:
      if 'epsfac' in r and 0<r['epsfac']<10:
        runs.remove(r)
  # Name each run based on its position in the list. 
  for i, r in enumerate(runs):
    r['name'] = 'T' + znt(i, 3)
    r['dir'] = rd + r['name']
  # Return the list of run dictionaries. 
  return runs

# =============================================================================
# ================================================ Load Programming Environment
# =============================================================================

def setEnvironment():
  # Headers make the error file more legible. 
  append('### Loading Modules', 'stdoe.txt')
  # Load all modules we'll need for compilation. 
  for m in mods():
    # Unload before loading to avoid collisions. 
    module('unload ' + m)
    module('load ' + m)
  # Also set up OpenMP. 
  os.environ['OMP_NUM_THREADS'] = str( nproc() )
  return module('list')

# =============================================================================
# ============================================================= Create Makefile
# =============================================================================

def setMakefile():
  # Compilation flags. Mostly OpenMP stuff. The -heap-arrays flags tells Intel
  # to put arrays on the heap, instead of the stack. The stack may be
  # marginally faster, but it can also cause crashes when arrays get too large. 
  flags = ('-132 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential ' + 
           '-lmkl_core -openmp -heap-arrays 1600')
  append('tuna.bin: source.f90', 'Makefile')
  # The library and include paths are different between Itasca and Eelpout. 
  append('\tifort -o tuna.bin source.f90 ' + incl() + libs() + flags, 
         'Makefile')
  return

# =============================================================================
# =================================================== Compile Binary Executable
# =============================================================================

def setBin():
  # Headers make the error file more legible. 
  append('### Compiling', 'stdoe.txt')
  # Create the Makefile. 
  setMakefile()
  print 'Compiling... ',
  # Compile. Keep track of how long it takes. 
  timer = time()
  out, err = bash('make tuna.bin')
  print format(time() - timer, '5.1f') + 's' + '\n'
  # If no binary executable was created, the compilation failed. 
  return 'tuna.bin' in files()

# #############################################################################
# ############################################################# Conduct One Run
# #############################################################################

# =============================================================================
# ====================================================== Create Parameters File
# =============================================================================

def setParams(run):
  # Write parameter keywords and values out, one per line. 
  for key, val in run.items():
    if key not in ('name', 'dir'):
      append(str(key) + ' = ' + str(val), 'params.in')
  return

# =============================================================================
# =========================================================== Create PBS Script
# =============================================================================

# If running on Itasca, the job is submitted to the queue using a PBS script. 
def setPBS(run):
  # In most cases, a 100s run will complete in an hour. But some parameter
  # combinations run slower than others. We include a factor of 5 to be safe. 
  hours = znt(1 if 'tmax' not in run else ceil(run['tmax']/20.) )
  # Write out the PBS script. 
  append('#!/bin/bash -l', 'tuna.pbs')
  # Indicate the size and length of the job. 
  append('#PBS -l nodes=1:ppn=16,walltime=' + hours + ':00:00', 'tuna.pbs')
  # Get email when the job begins, ends, or aborts. 
  append('#PBS -m abe', 'tuna.pbs')
  append('#PBS -M mceachern@physics.umn.edu', 'tuna.pbs')
  # Use the SandyBridge queue. That's 16-core nodes. 
  append('#PBS -q sb', 'tuna.pbs')
  # Specify run name. 
  append('#PBS -N '+run['name'], 'tuna.pbs')
  # Execute in the directory we've created for this run. 
  append('cd $PBS_O_WORKDIR', 'tuna.pbs')
  # Load the programming environment. 
  append('module purge', 'tuna.pbs')
  append('module load intel', 'tuna.pbs')
  append('module load mkl/10.2.5.035', 'tuna.pbs')
  # Specify the number of threads to use. 
  append('export OMP_NUM_THREADS=16', 'tuna.pbs')
  # Run and time the code. 
  append('time ./tuna.bin 1>>tuna.out 2>>stdoe.txt', 'tuna.pbs')
  return

# =============================================================================
# ================================================================== Run Driver
# =============================================================================

def runRun(run):
  # If the executable doesn't exist, bail. 
  if 'tuna.bin' not in files():
    print 'Skipping ' + run['name'] + ': Executable not found. ' + '\n'
    return False
  # Create a subdirectory for this run. 
  os.mkdir(run['dir'])
  os.chdir(run['dir'])
  # Grab the executable. 
  copy('../tuna.bin', '.')
  # Create input parameter file. 
  setParams(run)
  # If we're on Eelpout, we call the executable directly. 
  if local():
    print 'Running ' + run['name'] + '... ', 
    timer = time()
    out, err = bash('./tuna.bin', out='tuna.out')
    # Print the time it took to execute. 
    print format(time() - timer, '5.1f') + 's' + '\n'
  # On Itasca, we run using the queue. 
  else:
    # Create the queue submission script. 
    setPBS(run)
    print 'Submitting ' + run['name'] + ' to the queue... ', 
    out, err = bash('qsub tuna.pbs')
    # Print the queue id number. 
    print out
  # If anything was sent to stderr, something is wrong. 
  if len(err)>0:
    print 'Something is wrong. '
    print err
    return False
  # Otherwise, let's take a look at the output. 
  else:
    print out
    return True

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# =============================================================================
# ====================================================== Lists and Dictionaries
# =============================================================================

# Combine two dictionaries. Later values win in the case of collisions. 
def dsum(d0, d1):
  return dict( d0.items() + d1.items() )

# Get indeces for a given list. 
def index(x):
  return range( len(x) )

# Get all items present in both lists.
def intersect(l0, l1):
  return [ x for x in l0 if x in l1 ]

# Check if the given key is present and nonzero. 
def nonzero(dic, key):
  return key in dic and dic[key]!=0

# Return a given list with duplicate entries removed. Order may get scrambled. 
def ulist(x):
  return list( set(x) )

# =============================================================================
# ============================================================ Input and Output
# =============================================================================

# Append text to a file. 
def append(text, filename=None, eol='\n'):
  if filename is not None:
    with open(filename, 'a') as f:
      f.write(text + eol)
  return text

# Turns a list of numbers (1, 2, 3) into the string '1x2x3'. 
def by(x):
  return str( x[0] ) + 'x' + by( x[1:] ) if len(x)>1 else str( x[0] )

# Grab the names of all files in a given directory. 
def files(loc='./'):
  if not os.path.isdir(loc):
    return []
  return [ f for f in os.listdir(loc) if os.path.isfile(loc + f) ]

# Timestamp, to prevent runs from colliding. 
def now():
  return ( znt(lt().tm_year) + znt(lt().tm_mon, 2) + znt(lt().tm_mday, 2) + 
           '_' + znt(lt().tm_hour, 2) + znt(lt().tm_min, 2) +
           znt(lt().tm_sec, 2) )

# =============================================================================
# ================================================================= OS Handling
# =============================================================================

# Make a call as if in the terminal. Return stdout and stderr after maybe
# writing them to files. 
def bash(command, out='stdoe.txt', err='stdoe.txt'):
  o, e = Popen(command.split(), stdout=PIPE, stderr=PIPE).communicate()
  return append(o, out), append(e, err)

# For module load, unload, list, avail. 
def module(command, out='stdoe.txt', err='stdoe.txt'):
  o, e = bash(modcmd() + command, out=out, err=err)
  exec o
  return e

# Remove a file or dictionary. 
def rm(path):
  if os.path.isdir(path):
    rmtree(path)
  else:
    os.remove(path)
  return

# Truncates a number and turns it into a zero-padded string. 
def znt(x, width=0):
  return str( int(x) ).zfill(width)

# =============================================================================
# ============================================================ Itasca v Eelpout
# =============================================================================

# Include path for compilation. 
def incl():
  if local():
    return '-I/opt/intel/composer_xe_2011_sp1.9.293/mkl/include/intel64/lp64 '
  else:
    return '-I/soft/intel/mkl/10.2.5.035/include/em64t/lp64 '

# Libs path for compilation.
def libs():
  if local():
    return '-L/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64 '
  else:
    return '-L/soft/intel/mkl/10.2.5.035/lib/em64t '

# Check if we are running locally or on a supercomputer. 
def local():
  return host().endswith('spa.umn.edu')

# Path to the executable that tells Python how to interact with modules. 
def modcmd():
  if host().endswith('spa.umn.edu'):
    return '/usr/bin/modulecmd python '
  else:
    return os.environ['MODULESHOME'] + '/bin/modulecmd python '

# What modules are needed for compilation and execution?
def mods():
  if local():
    return ('intelcce',)
  else:
    return ('intel', 'mkl/10.2.5.035')

# On Eelpout, don't bother parallelizing. On Itasca, use as many cores as possible. 
def nproc():
  return 1 if local() else 16

# Where should the run be carried out? 
def runDir(runDirName):
  if local():
    return '/export/scratch/users/mceachern/' + runDirName + '_' + now() + '/'
  else:
    return '/scratch.global/mceache1/' + runDirName + '_' + now() + '/'

# Where is the source code stored? 
def srcDir():
  if local():
    return '/home/user1/mceachern/Desktop/tuna/'
  else:
    return '/home/lysakrl/mceache1/tuna/'
























'''







# #############################################################################
# ############################################################## Driver Plotter
# #############################################################################

# Let's look at the results quickly, for debugging. This isn't nearly as pretty
# as the real plotter. Itasca hasn't got the plotting libraries, so this whole
# chunk also has to be commented out when not running locally...

# =============================================================================
# ========================================================== Plotting Libraries
# =============================================================================

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.colors import LogNorm

# =============================================================================
# ================================================================== Data Input
# =============================================================================

# Convert a string to a complex or float. 
def com(x):
  if ',' not in x:
    return float(x)
  else:
    # Shave off the parentheses then split into real and imaginary parts. 
    re, im = x[1:-1].split(',')
    return (float(re) + float(im)*1j)

# Read in a real or complex array from file. 
def getArray(path, filename):
  # Filenames are caps insensitive. Find the file we want. 
  matches = [ x for x in os.listdir(path) if x.lower()==filename.lower() ]
  # If the file doesn't exist, bail. 
  if len(matches)==0:
    print 'WARNING: No file found matching \"' + path + '/' + filename + '\"'
    return False
  else:
    filename = matches[0]
    print 'Reading ' + path + '/' + filename
  # Read in the file as a list of strings. 
  arrayLines = open(path + '/' + filename, 'r').readlines()
  # The first line is the array dimensions. 
  dims = [ int(x) for x in arrayLines.pop(0).split() ]
  # Assemble a one-dimensional array large enough to hold all of the values. 
  # (This is much faster than appending as we go.) This means figuring out if
  # we want reals or complexes. 
  firstValue = com( arrayLines[0].split()[0] )
  nVals = np.prod(dims)
  dtype = np.complex if isinstance(firstValue, np.complex) else np.float
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
  # have opposite conventions for which index should change fastest. Note that
  # complex arrays are returned as complexes! 
  arr = np.transpose( np.reshape( vals, dims[::-1] ) )
  # If the array is full, return it. 
  if i==nVals:
    return arr
  # If we're missing time steps, the run may have crashed or run out of time.
  # We still want to look at the values. Return the time steps that happened. 
  elif filename=='t.out':
    print ('WARNING: Expected ' + str(nVals) + ' values for ' + path + '/' +
           filename + ' but found ' + str(i) )
    return arr[:i]
  elif len(dims)==3:
    stepsFound = i/( dims[0]*dims[1] )
    print ('WARNING: Expected ' + by(dims) + ' values for ' + path + '/' + 
           filename + ' but found ' + by( dims[0:2] + [stepsFound] ) )
    return arr[:, :, :stepsFound]
  # If a non-time-resolved array is incomplete, something is very wrong. 
  else: 
    print ('ERROR: Expected ' + str(nVals) + ' values for ' + path + '/' + 
           filename + ' but found ' + str(i) )
    exit()

# =============================================================================
# ====================================================== Title Number Formatter
# =============================================================================

def fmt(x):
  if x==0:
    return '0'
  elif 1e-3<abs(x)<1e3:
    return str( float(x) if float(x)!=int( float(x) ) else int( float(x) ) )
  else:
    snx = format(x, '.0e').replace('e', '\\cdot 10^{') + '}'
    return snx.replace('+0', '+').replace('-0', '-')

# =============================================================================
# ====================================================== Initialize Plot Window
# =============================================================================

# We want to transpose our subplots -- go down the columns instead of across the rows. 
def subplotPos(nRows, nCols, i):
  return nRows, nCols, 1 + (i/nRows) + (i%nRows)*nCols

def initPlot():
  # Set up LaTeX fonts. 
  rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica'], 'size':'18'})
  rc('text', usetex=True)
  rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}')
  # Set up the plot window. 
  fig = plt.figure(figsize=(20, 10), facecolor='white')
  fig.canvas.set_window_title('Tuna Driver')
  # Partition out subplots. 
  nRows = 3
  nCols = 4
  # Return a list of axes. Note that we flip the usual order of indeces. 
  return [ plt.subplot( *subplotPos(nRows, nCols, i) ) for i in range(nRows*nCols) ]

# =============================================================================
# =============================================================== Assemble Data
# =============================================================================

def colorParams(p99):
  return {'levels':np.linspace(-p99, p99, 50, endpoint=True), 
          'cmap':plt.get_cmap('seismic'), 'extend':'both'}


def vPlot():
  # Set up LaTeX fonts. 
  rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica'], 'size':'18'})
  rc('text', usetex=True)
  rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}')
  # Set up the plot window. 
  fig = plt.figure(figsize=(20, 10), facecolor='white')
  fig.canvas.set_window_title('Tuna Driver')
  # Grab the grid and Alfven speed. 
  path = 'T000'
  r, q = getArray(path, 'r.out'), getArray(path, 'q.out')
  X, Z = r*np.sin(q), r*np.cos(q)
  v = getArray(path, 'vA.out')
  # Figure out appropriate levels for the contours. 
  vmin, vmax = np.min(v), np.max(v)
  xmin = int( np.floor( np.log10(vmin) ) )
  xmax = int( np.ceil( np.log10(vmax) ) )
  levels = np.logspace( xmin, xmax, 7 )

  ticks = [ 10**i for i in range(xmin, xmax+1) ]
  labels = [ '$' + format(t, '.0f') + '$' for t in ticks ]

  con = plt.contourf( X, Z, v, levels=levels, norm=LogNorm() )
  cbar = plt.colorbar(con)
  cbar.set_ticks(ticks)
  cbar.set_ticklabels(labels)

  plt.show()





def showPlot(ax):
  # We'll just look at T000. 
  path = 'T000'
  # Grab coordinate arrays. 
  r, q = getArray(path, 'r.out'), getArray(path, 'q.out')
  X, Z = r*np.sin(q), r*np.cos(q)
  # Set the supertitle based on the last time step. 
  t = getArray(path, 't.out')
  plt.suptitle('$\mathrm{Snapshot \;\; at \;\; ' + fmt( float(t[-1]) ) + 's}$', fontsize=24)
  # Grab the data to plot. 
  names = ('Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez')
  data = [ getArray(path, n + '.out') for n in names ]
  # Split into real and imaginary components. 
  realNames = [ '$\mathbb{R} \; ' + n + '$' for n in names ]
  imagNames = [ '$\mathbb{I} \; ' + n + '$' for n in names ]
  realData = [ np.real(d) for d in data ]
  imagData = [ np.imag(d) for d in data ]
  # Assemble the subplots. We don't show a color bar, but we do set the color levels to be
  # symmetric around zero, and indicate the max color value in the title. 
  for i in range( len(names) ):
    # Real values. 
    p99 = np.percentile(np.abs( np.real(data[i]) ), 99)
    ax[i].contourf( X, Z, np.real(data[i][:, :, -1]), **colorParams(p99) )
    ax[i].set_title('$\mathbb{R} \;\;\; ' + names[i] + ' \;\;\; P_{99} = ' + fmt(p99) + '$')
    # Imaginary values. 
    p99 = np.percentile(np.abs( np.imag(data[i]) ), 99)
    ax[6 + i].contourf( X, Z, np.imag(data[i][:, :, -1]), **colorParams(p99) )
    ax[6 + i].set_title('$\mathbb{I} \;\;\; ' + names[i] + ' \;\;\; P_{99} = ' + fmt(p99) + '$')
  # Also, let's draw dipole outlines. 
  [ a.plot( X[i, :], Z[i, :], color=(0, 0, 0) ) for i in (0, -1) for a in ax ]
  [ a.plot( X[:, k], Z[:, k], color=(0, 0, 0) ) for k in (0, -1) for a in ax ]
  # And remove the axis ticks. This plot will be packed enough already. 
  [ ( a.set_xticks( [] ), a.set_yticks( [] ) ) for a in ax ]
  # Show the plot window now that it's all assembled. 
  plt.show()
  return





'''








# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


