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

# If running anything at crazy high resolution, use this flag to just request
# the maximum allowed runtime. 
USE_MAX_TIME = False

# Tuna includes default values, which it uses for any parameter not specified. 
parameters = {
              'jdrive':[4e-4], 
#              'bdrive':[10], 
#              'inertia':[1],
#              'lpp':[4],
              'ldrive':[6],
              'tmax':[300],
              'azm':[1, 2, 4, 8, 16, 32, 64],
              'fdrive':[0.010, 0.013, 0.016, 0.019, 0.022],
              'model':[4]
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
  global USE_MAX_TIME
  # A 100s run should run in about an hour if there are no inertial effects. We
  # ask for four just to be safe. Resolving inertial length scales slows this
  # significantly. We just request the maximum allowed runtime.  
  if USE_MAX_TIME is True:
    hours = '96'
  else:
    hours = znt(1 if 'tmax' not in run else ceil(run['tmax']/25.) )
  # Write out the PBS script. 
  append('#!/bin/bash -l', 'tuna.pbs')
  # Indicate the size and length of the job. 
  append('#PBS -l nodes=1:ppn=16,walltime=' + hours + ':00:00', 'tuna.pbs')
  # Get email when the job begins, ends, or aborts. 
  append('#PBS -m abe', 'tuna.pbs')
  append('#PBS -M mceachern@physics.umn.edu', 'tuna.pbs')
  # Use the SandyBridge queue. That's 16-core nodes. The sb128 queue allows a
  # higher walltime limit, but has fewer nodes. Use it only for inertial runs. 
  if USE_MAX_TIME is True:
    append('#PBS -q sb128', 'tuna.pbs')
  else:
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
#    return os.environ['MODULESHOME'] + '/bin/modulecmd python '
    return 'modulecmd python '

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

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


