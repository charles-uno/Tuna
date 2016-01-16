#!/usr/bin/env python

# Charles McEachern

# Fall 2015

# Note: This document wraps at column 80. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

import gc
import matplotlib
import os
from os.path import basename
# Change matplotlib settings to allow use over SSH without X forwarding. 
if 'DISPLAY' not in os.environ or os.environ['DISPLAY'] is '':
  matplotlib.use('Agg')
from matplotlib.colorbar import ColorbarBase
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import gridspec, rc
from matplotlib.colors import Normalize, LogNorm
from matplotlib.colors import LinearSegmentedColormap as LSC
# To easily draw a half-black-half-white semicircle. 
from matplotlib.patches import Wedge
import numpy as np
# The cPickle module is faster, but not always available. 
try:
  import cPickle as pickle
except ImportError:
  import pickle
from sys import argv, stdout
from time import localtime as lt, sleep, time

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():
  # Parse input from the terminal. 
  flags, paths, names = getArgs()

  # Initialize the Tuna Plotter. 
  TP = tunaPlotter(flags, paths)


  allParams = TP.getAllParams()

  for p in allParams:
    print '\t', p, allParams[p]

  print ''

  # Use the Tuna Plotter to create plots per the names and corresponding 
  # arguments from the terminal. 

  if 'a' in names:
    for model in (1, 2, 3, 4):
      for azm in (1, 8, 64):
        TP.plotAtmosphere(model=model, azm=azm)

  if 'c' in names:
    TP.runCheckup()

  if 'd' in names:
    for model in (1, 2, 3, 4):
      for azm in (1, 8, 64):
        TP.plotDownward(model=model, azm=azm)


  if 'e' in names:
    TP.plotExample()



  if 'f' in names:
    TP.plotFrequencies( *names['f'] )

  if 'g' in names:
    TP.plotGrid( *names['g'] )




  if 'l' in names:
    TP.plotLazy( *names['l'] )

  if 'p' in names:
    TP.paramSummary()

  if 'q' in names:
    for model in (1, 2, 3, 4):
      for azm in (1, 8, 64):
        TP.plotParallel(model=model, azm=azm)

  if 'x' in names:
    TP.pickle()





  '''

  if 'b' in names:
    for fdrive in (12, 14, 17, 20, 25):
      TP.plotB(fdrive)

  if 'c' in names:
    TP.plotConductivity()

  if 'd' in names:
    TP.plotDst()

  if 'e' in names:
    for fdrive in (12, 14, 17, 20, 25):
      TP.plotE(fdrive)

  if 'g' in names:
    TP.plotGrid()


  if 'm' in names:
    TP.plotM(4)

  if 'o' in names:
    TP.plotOne()

  if 'q' in names:
    TP.plotQuick()

  if 's' in names:
    for fdrive in (12, 14, 17, 20, 25):
      TP.plotS(fdrive)

#  for path in paths:
#    TP.plotV(path)
  '''

  return

# #############################################################################
# ######################################################### Plot Window Wrapper
# #############################################################################

# The window wrapper class includes utilities which are specific to the Tuna
# plotter, such as interfacing with the output data. 
class tunaPlotter:

  # ===========================================================================
  # ==================================================== Tuna Plotter Utilities
  # ===========================================================================

  # To avoid spending time reading in the same file more than once, data input
  # and access is centralized. 
  data = {}

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Tuna Plotter
  # ---------------------------------------------------------------------------

  def __init__(self, flags, paths):
    # Keep track of flags from the terminal. 
    self.flags = flags
    # If we'll be saving output, figure out where to put it. 
    self.outDir = '/home/user1/mceachern/Desktop/plots/plots_' + now() + '/'
    # Sort the paths before storing them. 
    self.paths = sorted(paths)
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------ Data Selection
  # ---------------------------------------------------------------------------

  # Look at the parameter input file for a given path and find a value. 
  def getParam(self, path, key):
    # Grab the input parameters. 
    params = '\n'.join( open(path + 'params.in', 'r').readlines() )
    # Make sure the parameter was actually provided in this run. 
    if key not in params:
      return None
    # Split out the line with the key we want. 
    line = params.split(key)[1].split('\n')[0]
    return num( line.strip(' =') )

  # Get a listing of all the parameters used in these runs. Essentially, we're
  # reconstructing the parameters dictionary at the top of the driver script. 
  def getAllParams(self):
    allParams = {}
    # Scroll through all the output directories. 
    for path in self.paths:
      # Grab the parameter input file. 
      paramLines = open(path + 'params.in', 'r').readlines()
      # Split each line into key and value. 
      for line in paramLines:
        # Watch out for input compatible with Bob's code. 
        if ',' in line:
          continue
        key = line.split('=')[0].strip()
        val = num( line.split('=')[-1] )
        # Add the key to our parameter dictionary, if it's not already there. 
        if key not in allParams:
          allParams[key] = []
        # Add the value, if it's not already there. 
        if val not in allParams[key]:
          allParams[key].append(val)
    # Return the parameter dictionary. 
    return allParams

  # Given a set of run parameters, find the desired run path. 
  def getPath(self, **kargs):
    # Check all the paths we have. 
    paths = self.paths
    # Weed out any that don't match. 
    for key in kargs:
      paths = [ p for p in paths if self.getParam(p, key)==kargs[key] ]
    # If there's anything other than exactly one match, something is wrong. 
    if len(paths)<1:
      print 'ERROR: No matching path found for ', kargs
      exit()
    elif len(paths)>1:
      print 'ERROR: Multiple matching paths found for ', kargs
      exit()
    else:
      return paths[0]

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------- Data Access
  # ---------------------------------------------------------------------------

  # If we're only plotting one time slice for each field, we had better make
  # sure it's a good one. 
  def getBestSlice(self, filename, ReIm='real'):
    # Grab the real or imaginary component of the array. 
    if ReIm=='real':
      arr = np.real( self.getArray(filename) )
    else:
      arr = np.imag( self.getArray(filename) )
    # Make sure it's safe to slice the array. 
    if len(arr.shape)>2:
      # Look at each time step's standard deviation. That's probably a good
      # proxy for wave activity. Return the one with the most. 
      nSteps = arr.shape[2]
      stdev = [ np.std( np.abs( arr[:, :, step] ) ) for step in range(nSteps) ]
      bestStep = np.argmax(stdev)
      return arr[:, :, bestStep], bestStep
    # Arrays without steps don't get sliced, but do still get returned. 
    else:
      print 'WARNING: Can\'t slice ' + filename
      return arr, None

  def getArray(self, filename):
    # If we haven't yet read in this array, we need to do that. 
    if filename not in self.data:
      self.data[filename] = readArray(filename)
    # Return the array. 
    return self.data[filename]

  def refresh(self):
    # Release all of the data we read in for this plot. 
    self.data = {}
    # Python's automatic garbage collection is not fast enough to handle the
    # amount of data running through this routine. Every time we make a plot,
    # let's also manually run garbage collection to release memory that's no
    # longer being used. 
    print 'Collecting garbage... ', gc.collect()
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Coordinate Handling
  # ---------------------------------------------------------------------------

  # For setting up the axes on a dipole or unwrapped plot. 
  def getCoords(self, path):
    # We'll be returning a coordinate dictionary to be passed as a list of
    # keyword arguments to the plot window's setParam method. 
    params = {}
    # Grab the coordinates. 
    r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
    # If this run is on the dayside, flip theta (which in effect flips X). 
    if path in self.paths and self.getParam(path, 'model')<3:
      q = -q
    # Dipole plots want GSE X and Z. 
    if '-u' not in self.flags:
      params['X'] = r*np.sin(q)
      params['Y'] = r*np.cos(q)
      params['xLabel'] = 'X_{GSE} \;\; (R_E)'
      params['yLabel'] = 'Z_{GSE} \;\; (R_E)'
      # Dipole plots are outlined. 
      params['outline'] = True
      # Snap the grid to nice round numbers. 
      if 0 < np.min( params['X'] ) <= 2 and 8 < np.max( params['X'] ) <= 12:
        params['xTicks'] = (0, 2, 4, 6, 8, 10)
        params['xTickLabels'] = ('$0$', '$2$', '$4$', '$6$', '$8$', '$10$')
        params['xLimits'] = (0, 10)
      if 3 < -np.min( params['Y'] ) <= 4 and 3 < np.max( params['Y'] ) <= 4:
        params['yTicks'] = (-4, -2, 0, 2, 4)
        params['yTickLabels'] = ('$-4$', '$-2$', '$0$', '$+2$', '$+4$')
        params['yLimits'] = (-4, 4)
    # Unwrapped plots use normalized cos theta and the McIlwain parameter. 
    else:
      params['Y'] = r/np.sin(q)**2
      params['X'] = np.cos(q)/np.sqrt(1 - np.min(r)/params['Y'])
      params['xLabel'] = '\\cos \\theta / \\cos \\theta_0'
#      params['yLabel'] = 'L = \\frac{r}{\\sin^2 \\theta} \;\; (R_E)'
      params['yLabel'] = 'L \;\; (R_E)'
      # The x axis needs to go from -1 to 1, and just needs to be labeled N/S. 
      params['xTicks'] = (-1, 0, 1)
      params['xTickLabels'] = ('$\\mathrm{S}$', '', '$\\mathrm{N}$')
      # In principle Lmin and Lmax could change, so let's only overwrite the
      # ticks and tick labels for the usual values, 1.5 and 10 (respectively). 
      if 1 < np.min( params['Y'] ) <= 2 and 10 <= np.max( params['Y'] ) < 12:
        params['yTicks'] = (2, 4, 6, 8, 10)
        params['yTickLabels'] = ('$2$', '$4$', '$6$', '$8$', '$10$')
    # Sometimes Python gets a little overenthusiastic about significant digits.
    # Round the domain boundaries to round numbers. 
    if 'xLimits' not in params:
      params['xLimits'] = ( float( format(np.min(params['X']), '.2e') ), 
                            float( format(np.max(params['X']), '.2e') ) )
    if 'yLimits' not in params:
      params['yLimits'] = ( float( format(np.min(params['Y']), '.2e') ), 
                            float( format(np.max(params['Y']), '.2e') ) )
    # Return these as a dictionary to be unpacked. 
    return params

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Display Handling
  # ---------------------------------------------------------------------------

  # Either display the plot or save it as an image. 
  def render(self, PW, name='plot.png'):
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # If we're supposed to save the plot as an image, do so. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + name
      return PW.render(filename)
    # Otherwise, show the image window. 
    else:
      return PW.render()

  # ===========================================================================
  # ======================================================== Tuna Plotter Plots
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Parameter Summary
  # ---------------------------------------------------------------------------

  # This method doesn't actually make a plot. It just prints out the params
  # used for each plot. 
  def paramSummary(self):
    # We'll keep a dictionary of dictionaries. 
    paramDict = {}
    # Scroll through the paths and grab the parameters used by each. 
    for path in self.paths:
      # Each path gets a dictionary of parameters. 
      paramDict[path] = {}
      # Grab the parameters input file. 
      paramLines = open(path + 'params.in', 'r').readlines()
      # Lines are always formatted as 'key = val'
      for line in paramLines:
        key, val = [ x.strip() for x in line.split('=') ]
        # Store this in the dictionary. 
        paramDict[path][key] = val
    # Don't bother printing out a parameter shared by all runs. 
    for key, val in paramDict[ self.paths[0] ].items():
      # Check if any path has a different value for this key. 
      for path in self.paths:
        if key not in paramDict[path] or paramDict[path][key]!=val:
          break
      # If all paths have the same value for this key, cut it. 
      else:
        for path in self.paths:
          if key in paramDict[path]:
            del paramDict[path][key]
    # Now print off the parameters used for each path. 
    for path in self.paths:
      # Just keep the immediate directory name as a label. 
      shortPath = path.split('/')[-2].split('_')[0]
      print col(shortPath + ': '),
      # Print out the params as a comma-delimited set of columns. 
      kparams = sorted( paramDict[path].items() )
      print ', '.join( col(key) + col(val) for key, val in kparams )
    return

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Last Run Checkup
  # ---------------------------------------------------------------------------

  # This should make it easy to SSH in and check on the status of the most
  # recent (or ongoing) run.  
  def runCheckup(self):
    # Path to the directory that holds the output directories. 
    root = '/export/scratch/users/mceachern/'
    outDirs = dirs(root)
    # Find the one that's been modified most recently. 
    newDir = max( (os.path.getmtime(d), d) for d in outDirs )[1]
    print ''
    print 'Checking on: ', newDir.replace(root, '')
    print ''
    # Get the paths for each run within the most recent output directory. 
    runDirs = [ d for d in dirs(newDir) if 'params.in' in os.listdir(d) ]
    # Use a dictionary of dictionaries to keep track of parameters. 
    paramDict = {}
    for d in runDirs:
      # Each path gets a dictionary of parameters. 
      paramDict[d] = {}
      # Grab the parameters input file. 
      paramLines = read(d + 'params.in')
      # Lines are always formatted as 'key = val'
      for line in paramLines:
        key, val = [ x.strip() for x in line.split('=') ]
        # Store this in the dictionary. 
        paramDict[d][key] = val
    # Don't bother printing out a parameter shared by all runs. 
    for key, val in paramDict[ runDirs[0] ].items():
      # Check if any path has a different value for this key. 
      for d in runDirs:
        if key not in paramDict[d] or paramDict[d][key]!=val:
          break
      # If all paths have the same value for this key, cut it. 
      else:
        for d in runDirs:
          if key in paramDict[d]:
            del paramDict[d][key]

    # Assume that all runs have the same parameters. 
    header = col('name') + ''.join( col( x[0] ) for x in sorted( paramDict[ runDirs[0] ].items() ) ) + col('time')

    print header

    # Scroll through each run and check on its progress. 
    for d in runDirs:

      name = d.replace(newDir, '').rstrip('/')
      tMax = self.getParam(d, 'tmax')
      tLast = num( read(d + 't.dat')[-1] )

      # Run name. 
      print col(name),

      # Run params. 
      print ''.join( col( x[1] ) for x in sorted( paramDict[d].items() ) ),

      # Status. 
      print col(tLast, 4) + '/ ' + col(tMax, 4),

      out, err = read(d + 'tuna.out'), read(d + 'stdoe.txt')

      if len(err)>0:
        print 'ERROR: ' + err[-1]
      elif 'Finished' in ''.join(out):
        print 'Done. '
      else:
        print 'Waiting...'

    return

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------------- Lazy Plot
  # ---------------------------------------------------------------------------

  def plotLazy(self, step=-1):
    # We will plot the real and imaginary components for each field. 
    fields = ('Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Jz')
    # Each run gets two rows: one real and one imaginary. 
    PW = plotWindow(len(fields), 2*len(self.paths), colorbar='sym', yPad=1)
    # Keep track of the runs we're comparing -- particularly their time stamps. 
    titles = []
    # Scroll through the runs. 
    for runNum, path in enumerate(self.paths):
      # Match each run's name up with its time stamp. 
      timeStamp = format(self.getArray(path + 't.out')[step], '.2f') + 's'
      shortPath = path.split('/')[-2].split('_')[0]
      titles.append('\\mathrm{' + shortPath + '\; at \;' + timeStamp + '}' ) 
      # Handle the coordinate grid and axis labels. 
      PW.setParam( nColors=12, **self.getCoords(path) )
      # Each field gets a column. 
      for col, name in enumerate(fields):
        # Label the column. 
        units = ( ' \\quad \mathrm{(' + 
                  {'B':'nT', 
                   'E':'\\frac{mV}{m}', 
                   'J':'\\frac{\\mu A}{m^2}'}[ name[0] ] +
                  ')}' )
        PW[col].setParam(colLabel = name[0] + '_' + name[1] + units)
        # Grab the data. 
        zComp = self.getArray(path + name + '.out')[:, :, step]
        # Plot real and imaginary components. 
        for ReIm in ('R', 'I'):
          # Figure out which row this plot actually belongs on. 
          row = runNum if ReIm=='R' else len(self.paths) + runNum
          # Label the row. 
          rowLabel = '\\mathbb{' + ReIm + '} \;\; \\mathrm{' + shortPath + '}'
          PW[row].setParam(rowLabel=rowLabel)
          # Grab the appropriate data component. 
          Z = np.real(zComp) if ReIm=='R' else np.imag(zComp)


          zmax = np.max(Z)
          imax = np.unravel_index( np.argmax(Z), Z.shape)
          coords = self.getCoords(path)
          xmax = coords['X'][imax]
          ymax = coords['Y'][imax]


          print 'max ' + name + ' = ', zmax, ' at ', imax, ' -> ', xmax, ', ', ymax


          # Plot the contour. 
          PW[col, row].setContour(Z)
    # Add a title to the top so we know the time stamp(s). 
    PW.setParam( title= ' \\quad '.join(titles) )
    # Render the plot window, either as a window or as an image. 
    return self.render(PW, 'lazy.png')

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------- Parallel Plot
  # ---------------------------------------------------------------------------

  def plotParallel(self, model=1, azm=1, step=-1):

    # Compare no E3, no J3, yes J3. 

#    fields = ('Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Jz')
#    PW = plotWindow(len(fields), 3, colorbar='sym', yPad=1)
    PW = plotWindow(4, 3, colorbar='sym', yPad=1)

    # Magnetic constant in nH/m. 
    mu0 = 1256.63706

    # Iterate through the rows. 
    for row in range(3):
      # -1 means no inertia, and automatic Boris factor. With a manual Boris
      # factor of 1, the parallel electric field doesn't show up at all. 
      inertia = 1 if row==2 else -1
      epsfac = 1 if row==0 else -1
      # Find the appropriate data path. 
      path = self.getPath(inertia=inertia, epsfac=epsfac, model=model, azm=azm)

      # Label the column. 
      rowLabels = ('\\mathrm{No} \\;\\; E_\\parallel', 
                   '\\mathrm{No} \\;\\; J_\\parallel', 
                   '\\mathrm{Yes} \\;\\; J_\\parallel')
      PW[row].setParam( rowLabel=rowLabels[row] )
      # Handle the coordinate grid and axis labels. 
      PW.setParam( nColors=10, **self.getCoords(path) )

      re, im = ' \\mathbb{R}\\mathrm{e} \\; ', ' \\mathbb{I}\\mathrm{m} \\; '

      units = {'B':' \\;\\; \\mathrm{(\\frac{mV}{m})} ', 
               'E':' \\;\\; \\mathrm{(nT)} ',
               'J':' \\;\\; \\mathrm{(\\frac{\\mu A}{m^2})} ',
               'S':' \\;\\; \\mathrm{(\\frac{mW}{m^2})} '}

      # Toroidal Poynting flux. Imaginary components. Conjugate of By. 
      Ex = np.imag( self.getArray(path + 'Ex.dat')[..., step] )
      By = -np.imag( self.getArray(path + 'By.dat')[..., step] )
      PW[0, row].setContour( Ex*By/mu0 )
      PW[0].setParam(colLabel='\\frac{1}{\\mu_0} E_x B_y^* \\quad \\mathrm{(\\frac{mW}{m^2})}')

      # Poloidal Poynting flux. 
      Ey = np.real( self.getArray(path + 'Ey.dat')[..., step] )
      Bx = np.real( self.getArray(path + 'Bx.dat')[..., step] )
      PW[1, row].setContour( -Ey*Bx/mu0 )
      PW[1].setParam(colLabel='-\\frac{1}{\\mu_0} E_y B_x^* \\quad \\mathrm{(\\frac{mW}{m^2})}')

      # Parallel electric field. 
      Ez = np.imag( self.getArray(path + 'Ez.dat')[..., step] )
      PW[2, row].setContour(Ez)
      PW[2].setParam(colLabel='\\mathbb{I}\\mathrm{m} \\; E_z \\quad \\mathrm{(\\frac{mV}{m})}')

      # Parallel current. 
      Jz = np.imag( self.getArray(path + 'Jz.dat')[..., step] )
      PW[3, row].setContour(Jz)
      PW[3].setParam(colLabel='\\mathbb{I}\\mathrm{m} \\; J_z \\quad \\mathrm{(\\frac{\\mu A}{m^2})}')

#      # Go through the columns. 
#      for col, name in enumerate(fields):
#
#        arr = self.getArray(path + name + '.dat')[..., step]
#
#        realMax = np.median( np.abs( np.real(arr) ) )
#        imagMax = np.median( np.abs( np.imag(arr) ) )
#
#        if realMax>imagMax:
#          PW[col, row].setContour( np.real(arr) )
#          colLabel = re + name[0] + '_' + name[1] + units[ name[0] ]
#        else:
#          PW[col, row].setContour( np.imag(arr) )
#          colLabel = im + name[0] + '_' + name[1] + units[ name[0] ]
#        PW[col].setParam(colLabel=colLabel)

    # Read in the data. For the moment, let's look at the last time step. 
    tLabel = str( int( self.getArray(path + 't.dat')[step] ) ) + 's'

    # Plot supertitle. 
    PW.setParam(title='\\mathrm{Model \\;\\; ' + str(model) +
                      ' \\;\\; with \\;\\; } m = \\mathrm{' + str(azm) +
                      ' \\;\\; at \\;\\; } t = \\mathrm{' + tLabel + '}')

    pngName = 'parallel_' + str(model) + '_' + str(azm).zfill(3) + '_' + tLabel + '.png'


    return self.render(PW, pngName)

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Atmospheric Plot
  # ---------------------------------------------------------------------------

  def plotAtmosphere(self, azm=1, model=1):

    fields = ('Br', 'BfE', 'BfI', 'BqE', 'BqI')

    PW = plotWindow(nCols=len(fields), nRows=2, colorbar='sym', yPad=1)

    re, im = ' \\mathbb{R}\\mathrm{e} \\; ', ' \\mathbb{I}\\mathrm{m} \\; '

    nT = ' \\; \\mathrm{(nT)} '

    def texName(name):
      sub = name[1].replace('f', '\\phi').replace('q', '\\theta')
      if not len(name)>2 or name[2]=='I':
        return name[0] + '_' + sub + ' \\mathrm{\\; at \\;} R_I '
      else:
        return name[0] + '_' + sub + ' \\mathrm{\\; at \\;} R_E '

    # Iterate through the rows. 
    for row in range(2):
      # -1 means no inertia, and automatic Boris factor. With a manual Boris
      # factor of 1, the parallel electric field doesn't show up at all. 
      inertia = (-1)**(row+1)
      epsfac = (-1)**(row)

      # Find the appropriate data path. 
      path = self.getPath(inertia=inertia, epsfac=epsfac, model=model, azm=azm)

      # Label the column. 
      rowLabels = ('\\mathrm{Inertia \\; Off}', 
                   '\\mathrm{Inertia \\; On}')
      PW[row].setParam( rowLabel=rowLabels[row] )

      q = self.getArray(path + 'q.dat')
      lat = q[:, 0]*180/np.pi
      t = self.getArray(path + 't.dat')

      for col, name in enumerate(fields):

        PW.setParam(x=t, xLabel='\\mathrm{Time \\;\\; (s)}', y=lat, 
                    yLabel='\\mathrm{Latitude \\;\\; (^\\circ)}')

        # Grab the field data. 
        arr = self.getArray(path + name + '.dat')[:, 0, :]
        # Figure out if we want the real or imaginary component. 
        realMed = np.median( np.abs( np.real(arr) ) )
        imagMed = np.median( np.abs( np.imag(arr) ) )
        if realMed>imagMed:
          PW[col, row].setContour( np.real(arr) )
          colLabel = re + texName(name) + nT
        else:
          PW[col, row].setContour( np.imag(arr) )
          colLabel = im + texName(name) + nT
        PW[col].setParam(colLabel=colLabel)

    # Plot supertitle. 
    PW.setParam(title='\\mathrm{Model \\;\\; ' + str(model) +
                      ' \\;\\; with \\;\\; } m = \\mathrm{' + str(azm) + '}')

    pngName = 'atm_' + str(model) + '_' + str(azm).zfill(3) + '.png'

    return self.render(PW, pngName)



  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Downward Flux Plot
  # ---------------------------------------------------------------------------

  def plotDownward(self, azm=1, model=1):

    # Magnetic constant in nH/m. 
    mu0 = 1256.63706

    PW = plotWindow(nCols=4, nRows=2, colorbar='sym', yPad=1)

    north, south = '\\mathrm{N} \\; ', '\\mathrm{S} \\; '

    torName = '\\frac{1}{\\mu_0} E_x B_y^*'
    polName = '-\\frac{1}{\\mu_0} E_y B_x^*'

    units = ' \\; \\mathrm{(\\frac{mW}{m^2})} '

    # Iterate through the rows. 
    for row in range(2):
      # -1 means no inertia, and automatic Boris factor. With a manual Boris
      # factor of 1, the parallel electric field doesn't show up at all. 
      inertia = (-1)**(row+1)
      epsfac = (-1)**(row)
      # Find the appropriate data path. 
      path = self.getPath(inertia=inertia, epsfac=epsfac, model=model, azm=azm)
      # Label the column. 
      rowLabels = ('\\mathrm{Inertia \\; Off}', 
                   '\\mathrm{Inertia \\; On}')
      PW[row].setParam( rowLabel=rowLabels[row] )
      # Figure out coordinates. 
      q = self.getArray(path + 'q.dat')
      lat = q[:, 0]*180/np.pi
      t = self.getArray(path + 't.dat')

      # Toroidal fields. Imaginary. Complex conjugate for By. 
      Ex = np.imag( self.getArray(path + 'Ex.dat') )
      By = -np.imag( self.getArray(path + 'By.dat') )

      # Poloidal fields. 
      Ey = np.real( self.getArray(path + 'Ey.dat') )
      Bx = np.real( self.getArray(path + 'Bx.dat') )

      # We want the downward Poynting flux at each ionosphere. That means we
      # need to take the edge slices in k, and flip the southern ones. 
      SNtor = Ex[:, 0, :]*By[:, 0, :]/mu0
      SStor = -Ex[:, -1, :]*By[:, -1, :]/mu0
      SNpol = Ey[:, 0, :]*Bx[:, 0, :]/mu0
      SSpol = -Ey[:, -1, :]*Bx[:, -1, :]/mu0

      PW[0, row].setContour(SNtor)
      PW[0].setParam(colLabel=torName + units + north)

      PW[1, row].setContour(SStor)
      PW[1].setParam(colLabel=torName + units + south)

      PW[2, row].setContour(SNpol)
      PW[2].setParam(colLabel=polName + units + north)

      PW[3, row].setContour(SSpol)
      PW[3].setParam(colLabel=polName + units + south)

#      # Toroidal Poynting flux. Imaginary components. Conjugate of By. 
#      Ex = np.imag( self.getArray(path + 'Ex.dat')[..., step] )
#      By = -np.imag( self.getArray(path + 'By.dat')[..., step] )
#      PW[0, row].setContour( Ex*By/mu0 )
#      PW[0].setParam(colLabel='\\frac{1}{\\mu_0} E_x B_y^* \\quad \\mathrm{(\\frac{mW}{m^2})}')

#      # Poloidal Poynting flux. 
#      Ey = np.real( self.getArray(path + 'Ey.dat')[..., step] )
#      Bx = np.real( self.getArray(path + 'Bx.dat')[..., step] )
#      PW[1, row].setContour( -Ey*Bx/mu0 )
#      PW[1].setParam(colLabel='-\\frac{1}{\\mu_0} E_y B_x^* \\quad \\mathrm{(\\frac{mW}{m^2})}')

    PW.setParam(x=t, xLabel='\\mathrm{Time \\;\\; (s)}', y=lat, 
                    yLabel='\\mathrm{|Latitude| \\;\\; (^\\circ)}')

    # Plot supertitle. 
    PW.setParam(title='\\mathrm{Model \\;\\; ' + str(model) +
                      ' \\;\\; with \\;\\; } m = \\mathrm{' + str(azm) + '}')

    pngName = 's_' + str(model) + '_' + str(azm).zfill(3) + '.png'

    return self.render(PW, pngName)



  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------ Frequency Plot
  # ---------------------------------------------------------------------------

  def plotFrequencies(self, boris=1):
    # Plotting dimensionless ratios for each ionospheric profile. 
    models = (1, 2, 3, 4)

#    PW = plotWindow(7, len(models), colorbar='log')
    PW = plotWindow(7, len(models), colorbar='log', yPad=1)

    # Magnetic constant in nH/m. 
    mu0 = 1256.63706
    # Electric constant, in mF/m. 
    eps0 = 8.854e-9
    # The Boris factor is applied to the parallel electric constant. 
    epsPara = eps0*boris
    # Electron charge in MC. 
    qe = 1.60218e-25
    # Electron mass in g. 
    me = 9.10938e-28
    # Earth radius in Mm. 
    RE = 6.378388

    # Wavelength should be order of an RE. 
    k = 1./RE

#    # The fastest driving we have to worry about is a 40s period.  
#    w = 1./40
#    wLabel = str( int(1000*w) ) + ' mHz'

#    # Label the plot and its columns. 
#    if boris==1:
#      borisLabel = ''
#    else:
#      sn = format(boris, '.0e').replace('e', '\\cdot 10^{').replace('+0', '').replace('-0', '-') + '}'
#      if sn[0]=='1':
#        borisLabel = sn.split('cdot')[-1]
#      else:
#        borisLabel = sn

#    PW.setParam(title='\\mathrm{Frequency \;\; Ratios \;\; at \;\; ' + wLabel +
#                      ' \;\; with \;\; \\epsilon_\\parallel = ' + borisLabel +
#                      ' \\epsilon_0}')

#    PW[0].setParam(colLabel='\\omega^2 / \\omega_p^2')
#    PW[1].setParam(colLabel='\\omega \\nu_{\\parallel} / \\omega_p^2')
#    PW[2].setParam(colLabel='\\omega / \\nu_{\\parallel}')
#    PW[3].setParam(colLabel='\\frac{\\sigma_0}{\\epsilon_\\parallel} \\delta \\! t')
#    PW[4].setParam(colLabel='\\frac{\\sigma_H}{\\epsilon_\\bot} \\delta \\! t')
#    PW[5].setParam(colLabel='\\frac{\\sigma_P}{\\epsilon_\\bot} \\delta \\! t')
#    PW[6].setParam(colLabel='\\nu \; \\delta \\! t')

    PW.setParam(title='\\mathrm{Compare \\; to \\; \\sim 10^{-2} Hz}')

    PW[0].setParam(colLabel='\\nu')
    PW[1].setParam(colLabel='\\omega_\\parallel')
    PW[2].setParam(colLabel='c_\\parallel k')
    PW[3].setParam(colLabel='v_A k')
    PW[4].setParam(colLabel='\\sigma_P / \\epsilon_\\bot')
    PW[5].setParam(colLabel='\\sigma_H / \\epsilon_\\bot')
    PW[6].setParam(colLabel='\\sigma_0 / \\epsilon_0')

    # Each row is from a different model. 
    for row, model in enumerate(models):

      # Find the path that goes with this model. 
      path = self.getPath(model=model)

      # Label the row. 
      PW[row].setParam(rowLabel='\\mathrm{Model \; ' + str(model) +  '}')

      # Handle the coordinates and axes. 
      PW.setParam( outline=True, nColors=14, **self.getCoords(path) )

#      # Time step in seconds, based on the grid and on the plasma frequency,
#      # before applying the Boris factor. 
#      dtGrid = num( open(path + 'dt.out').readlines()[0] )
#      dtPlasma = num( open(path + 'dt.out').readlines()[1] )

      # Number density is printed in cm^-3 but we want it in Mm^-3. 
      n = self.getArray(path + 'n.out')*1e24

      # Condictivity was printed in S/m and we want it in mS/m. 
      sig0 = self.getArray(path + 'sig0.out')/1e3
      sigP = self.getArray(path + 'sigP.out')/1e3
      sigH = self.getArray(path + 'sigH.out')/1e3

      # Perpendicular dielectric constant in units of eps0, convert to mF/m. 
      epsPerp = self.getArray(path + 'epsPerp.out')*eps0

      # Compute the plasma frequency. 
      wp = np.sqrt( n*qe**2 / (me*epsPara) )

      # Compute the collision frequency. 
      nu = n*qe**2 / (me*sig0)

      # Compute the Boris-adjusted speed of light. 
      cc = np.ones(n.shape)/(mu0*epsPara)

      # Compute the Alfven speed. 
      vv = 1/(mu0*epsPerp)


      clip = 1e8

#      PW[0, row].setContour(nu)
#      PW[1, row].setContour(wp)
#      PW[2, row].setContour(np.sqrt(cc)*k)
#      PW[3, row].setContour(np.sqrt(vv)*k)
      PW[4, row].setContour(sigP/epsPerp)
      PW[5, row].setContour(sigH/epsPerp)
      PW[6, row].setContour( np.minimum(clip, sig0/epsPara) )

      print 'Max sP/ep = ', np.max(sigP/epsPerp)
      print 'Max sH/ep = ', np.max(sigH/epsPerp)
      print 'Min s0/e0 = ', np.min(sig0/epsPara)

      print 'Max e0/s0 = ', np.max(epsPara/sig0)
      print 'Max nu/wp^2 = ', np.max(nu/wp**2)

      print ''

#      # Integrating factors... 
#      s0 = sig0*dt/epsPara
#      sH = sigH*dt/epsPerp
#      sP = sigP*dt/epsPerp

      # Throw these things on the plot. Let's clip the values that are
      # basically infinity. 
#      clip = 1e2
#      PW[0, row].setContour( np.minimum(clip, (w/wp)**2) )
#      PW[1, row].setContour( np.minimum(clip, w*nu/wp**2) )
#      PW[2, row].setContour( np.minimum(clip, w/nu) )

    # Render the plot window, either as a window or as an image. 
    return self.render(PW, 'cutoff.png')

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------------- Data Pickler
  # ---------------------------------------------------------------------------

  # This routine takes the time to read in all of the arrays from a run, then
  # writes them back out as pickles for easy access later. 
  def pickle(self):
    # Scroll through the output directories. 
    for path in self.paths:
      # Print whenever we go to a new directory. 
      print ''
      print path
      # Find all of the data arrays to be read in, and scroll through them. 
      datFiles = files(path, end='.dat')
      for f in datFiles:
        print '\t' + basename(f)
        # We call readArray directly. The Tuna Plotter wants to store a copy of
        # everything we read in but we don't actually care about this data. 
        readArray(f, indent='\t\t')
    return




  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------------- Grid Plot
  # ---------------------------------------------------------------------------

  def plotGrid(self, stride=1):
    # One plot cell, no color bar. 
    PW = plotWindow(colorbar=None)

#    plt.axis('equal')

    PW.setParam(title='\\mathrm{Dipole \\; Grid}')

    # Just grab the first path. 
    path = self.paths[0]
    # Grab the coordinates. 
    coords = self.getCoords(path)
    PW.setParam(**coords)
    # Draw some lines. 
    x, y = coords['X'], coords['Y']
    n1, n3 = x.shape
    [ PW.setLine( x[i, :], y[i, :] ) for i in range(0, n1, stride) ]
    [ PW.setLine( x[:, k], y[:, k] ) for k in range(0, n3, stride) ]
    # Show the plot. 
    return self.render(PW, 'grid.png')








  # ---------------------------------------------------------------------------
  # -------------------------------------------------------------- Example Plot
  # ---------------------------------------------------------------------------

  def plotExample(self):
    # One plot cell, no color bar. 
    PW = plotWindow(4, 4, colorbar='sym')

    PW.setParam(title='\\mathrm{Dipole \\; Grid}')

    for col, model in enumerate( (1, 2, 3, 4) ):

      for row, azm in enumerate( (1, 8, 64) ):

        path = self.getPath(inertia=1, model=model, azm=azm)

        PW[col, row].setParam( **self.getCoords(path) )







#    # Just grab the first path. 
#    path = self.paths[0]
#    # Grab the coordinates. 
#    coords = self.getCoords(path)
#    PW.setParam(**coords)
#    # Draw some lines. 
#    x, y = coords['X'], coords['Y']
#    n1, n3 = x.shape
#    [ PW.setLine( x[i, :], y[i, :] ) for i in range(0, n1, stride) ]
#    [ PW.setLine( x[:, k], y[:, k] ) for k in range(0, n3, stride) ]


    # Show the plot. 
    return self.render(PW, 'grid.png')



















  '''

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------------ Dst Plot
  # ---------------------------------------------------------------------------

  def plotDst(self):
    # Shorthand for legibility. 
    symh = 'Sym\\text{-}H'
    # Set up the window. 
    PW = plotWindow(2, 1, colorbar=False, xPad=2)
    PW.setTitle('\mathrm{' + symh + ' \;\; Frequency \;\; Breakdown \;\; ' +
                ' for \;\; June \;\; 2013 \;\; Storm}')
    # Grab the data that we downloaded from NASA CDAWeb. 
    filename = '/home/user1/mceachern/Desktop/tuna-old/symh/symh_20130601.txt'
    data = [ line for line in open(filename, 'r').readlines() if line.strip() ]
    t, SYMH = [], []
    for line in data:
      year, day, hour, minute, value = [ int(col) for col in line.split() ]
      t.append( 24*60*day + 60*hour + minute )
      SYMH.append(value) 
    t = np.array(t) - t[0]
    dt, tRange = t[1], t[-1]
    SYMH = np.array(SYMH)
    # Break SYMH down into Fourier modes. 
    nModes = 4000
    amplitudes = np.zeros(nModes, dtype=np.float)
    for m in range(nModes):
      harmonic = np.cos( m*np.pi*t / tRange )
      amplitudes[m] = np.sum( SYMH*harmonic*dt ) / np.sum( dt*harmonic**2 )
    # Plot SYMH. 
    PW.setLine(t, SYMH, pos=(0, 0), color='b')
    # Set title and labels. 
    PW.setTitle( '\mathrm{' + symh + '\;\; Data \;\; from \;\; NASA \;\; ' + 
                 'CDAWeb}', pos=(0, 0) )
    PW.setYlabel( '\mathrm{' + symh + ' \;\; (nT)}', pos=(0, 0) )
    PW.setXlabel( '\mathrm{Time \;\; (Minutes)}', pos=(0, 0) )
    # Optionally, put a few Fourier components over the data. 
    if True:
      nShow = 20
      PW.setTitle( '\mathrm{' + symh + ' \;\; from \;\; NASA \;\; CDAWeb ' +
                   ' \;\; with \;\; ' + str(nShow) + ' \;\; Fourier \;\; ' +  
                   'Modes}', pos=(0, 0) )
      # Reconstruct the data from cosines...
      reconstruction = np.zeros(len(t), dtype=np.float)
      for m in range(nShow):
        harmonic = np.cos( m*np.pi*t / tRange )
        amplitude = np.sum( SYMH*harmonic*dt ) / np.sum( dt*harmonic**2 )
        reconstruction = reconstruction + amplitude*harmonic
      # Add it to the plot. 
      PW.setLine(t, reconstruction, pos=(0, 0), color='r')
    # Plot Fourier amplitudes against their frequencies in mHz. 
    frequencies = np.array( range(nModes) )*1000./(60*2*tRange)
    PW.setLine( frequencies[1:], amplitudes[1:], pos=(1, 0), color='b' )
    PW.setXlog( pos=(1, 0) )
    PW.setYlog( pos=(1, 0) )
    # Set title and labels. 
    PW.setTitle( '\mathrm{' + symh + ' \;\; Fourier \;\; Amplitudes}', 
                 pos=(1, 0) )
    PW.setYlabel( '\mathrm{Amplitude \;\; (nT)}', pos=(1, 0) )
    PW.setXlabel( '\mathrm{Frequency \;\; (mHz)}', pos=(1, 0) )
    # Set limits. They don't quite line up with the data, but that's OK. 
    fMin, fMax = 1e-2, 10
    PW.setXlimits( (fMin, fMax), pos=(1, 0) )
    f = np.linspace(fMin, fMax, 10000)
    # Let's not do a fit, per se -- let's eyeball the top of the distribution. 
    # Work in units of 20mHz. 
    intercept, slope = np.log(1.e-2), -0.9
    scale = 20. 
    # Plot the fit. 
    fit = np.exp(intercept) * np.power(f/scale, slope)
    label = (format( np.exp(intercept), '.2f' ) + '\; \mathrm{nT} ' + 
            '\cdot \left( \\frac{f}{' + format(scale, '.0f') +
            ' \; \mathrm{mHz} } \\right) ^{' + 
            format(slope, '.1f') + '}')
    PW.setLine( f, fit, pos=(1, 0), color='r', label=label)
    # Display the legend. 
    PW.setLegend( pos=(1, 0) )
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      filename = self.outDir + 'SYMH.png'
      PW.save(filename)
      return
    else:
      return PW.show()

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------------ One Plot
  # ---------------------------------------------------------------------------

  # This lets us zoom in on a single contour plot. 
  def plotOne(self, step=12, ReIm='R', name='Ex'):
    # Create the window. 
    PW = plotWindow(colorbar='sym', nColors=12, nTicks=11)
    # We only expect to get one path. If more are given, we take whichever one
    # is listed last, I guess. 
    for p in self.runs:
      path = p
    # Grab the coordinates and label the axes. 
    X, Y = self.getCoords(path)
    xLabel, yLabel = self.getCoordNames()
    PW.setXlabel(xLabel)
    PW.setYlabel(yLabel)
    PW.setTitle( name[0] + '_' + name[1] )
    # Grab the real or imaginary component of a slice of the data. 
    comp = np.real if ReIm=='R' else np.imag
    Z = comp( self.getArray(path + name + '.out')[:, :, step] )
    # Plot the contour. 
    PW.setContour(X, Y, Z)
    # Draw the dipole outline. 
    [ PW.setLine( X[i, :], Y[i, :] ) for i in (0, -1) ]
    [ PW.setLine( X[:, k], Y[:, k] ) for k in (0, -1) ]
    # Matplotlib wants to cut the unwrapped x axis a little short. 
    if '-u' in self.flags:
      PW.setXlimits( (-1, 1) )
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'lazy.png'
      return PW.save(filename)
    else:
      return PW.show()





  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Magnetic Field Plot
  # ---------------------------------------------------------------------------

  def plotB(self, fdrive):
    # We keep track of run parameters as strings. 
    if not isinstance(fdrive, str):
      fdrive = format(fdrive, '.0f')
    # For the Poynting flux, we want to see a 4x4 grid. Drive frequency is
    # constant across the window (so we will make a plot for each frequency). 
    # Modenumber is constant across each row. The columns are dayside toroidal,
    # nightside toroidal, dayside poloidal, nightside poloidal. Each run just
    # shows just the best frame from a run. All plot cells share a color bar. 

    PW = plotWindow(6, 4, colorbar='sym', yPad=1)

    # Figure out what each column will be. 
    columns = [ (s, p) for p in ('Bx', 'By', 'Bz') for s in ('Day', 'Night') ]
    # For the moment, we're only looking at storm time, not calm. 
    conditions = 'Stormtime '
    # Set plot supertitle. 
    PW.setTitle('\mathrm{' + conditions + ' \;\; Magnetic \;\; Field \;\; ' +
                ' Components \;\; (nT) \;\; at \;\; }' + fdrive +
                ' \mathrm{mHz}')
    # Squish the subplots a little bit to make the title stand out. 
    plt.subplots_adjust(top=0.88)
    # Each row is a different modenumber. 
    for row, azm in enumerate( ('1', '4', '16', '64') ):
      # Label each row my modenumber. 
      PW.setRowLabel('m = ' + azm, row)
      # Columns have day/night toroidal/poloidal polarizations. 
      for col, sideComp in enumerate(columns):
        # Find the directory that holds our data. 
        path = self.getRunPath(azm=azm, side=sideComp[0], fdrive=fdrive)
        # Coordinate arrays. Note that the sign of X is flipped. 
        r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
        t = self.getArray(path + 't.out')
        X = -r*np.sin(q) if sideComp[0]=='Day' else r*np.sin(q)
        Z = r*np.cos(q)
        # Note that we want the imaginary component for By. 
        ReIm = 'imag' if sideComp[1]=='By' else 'real'
        # Grab the field values and plot them as a contour. 
        B, step = self.getBestSlice(path + sideComp[1] + '.out', ReIm=ReIm)
        PW.setContour( X, Z, B, (col, row) )


        # Put together a title for this plot from its parameters. 
        title = (sideComp[0] + ' \;\; \mathbb{' + ReIm[0].upper() +
                 '} \;\; B_' + sideComp[1][1])
        PW.setTitle(title, (col, 0) )


        # Draw dipole outline. 
        [ PW.setLine( X[i, :], Z[i, :], (col, row) ) for i in (0, -1) ]
        [ PW.setLine( X[:, k], Z[:, k], (col, row) ) for k in (0, -1) ]
    # Label the axes. 
    PW.setXlabel( 'X_{GSE} \;\; (R_E)')
    PW.setYlabel( 'Z_{GSE} \;\; (R_E)')
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'B_' + fdrive + 'mHz.png'
      return PW.save(filename)
    else:
      return PW.show()

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Electric Field Plot
  # ---------------------------------------------------------------------------

  def plotE(self, fdrive):
    # We keep track of run parameters as strings. 
    if not isinstance(fdrive, str):
      fdrive = format(fdrive, '.0f')
    # For the Poynting flux, we want to see a 4x4 grid. Drive frequency is
    # constant across the window (so we will make a plot for each frequency). 
    # Modenumber is constant across each row. The columns are dayside toroidal,
    # nightside toroidal, dayside poloidal, nightside poloidal. Each run just
    # shows just the best frame from a run. All plot cells share a color bar. 

    PW = plotWindow(4, 4, colorbar='sym', yPad=1)

    # Figure out what each column will be. 
    columns = [ (s, p) for p in ('Ex', 'Ey') for s in ('Day', 'Night') ]
    # For the moment, we're only looking at storm time, not calm. 
    conditions = 'Stormtime '
    # Set plot supertitle. 
    PW.setTitle('\mathrm{' + conditions + ' \;\; Electric \;\; Field \;\; ' +
                ' Components \;\; (\\frac{mV}{m}) \;\; at \;\; }' + fdrive +
                ' \mathrm{mHz}')
    # Squish the subplots a little bit to make the title stand out. 
    plt.subplots_adjust(top=0.88)
    # Each row is a different modenumber. 
    for row, azm in enumerate( ('1', '4', '16', '64') ):
      # Label each row my modenumber. 
      PW.setRowLabel('m = ' + azm, row)
      # Columns have day/night toroidal/poloidal polarizations. 
      for col, sideComp in enumerate(columns):
        # Find the directory that holds our data. 
        path = self.getRunPath(azm=azm, side=sideComp[0], fdrive=fdrive)
        # Coordinate arrays. Note that the sign of X is flipped. 
        r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
        t = self.getArray(path + 't.out')
        X = -r*np.sin(q) if sideComp[0]=='Day' else r*np.sin(q)
        Z = r*np.cos(q)
        # Note that we want the imaginary component for By. 
        ReIm = 'imag' if sideComp[1]=='Ex' else 'real'
        # Grab the field values and plot them as a contour. 
        E, step = self.getBestSlice(path + sideComp[1] + '.out', ReIm=ReIm)
        PW.setContour( X, Z, E, (col, row) )


        # Just column headers. Not cell titles. 

        # Put together a title for this plot from its parameters. 
        title = (sideComp[0] + ' \;\; \mathbb{' + ReIm[0].upper() +
                 '} \;\; E_' + sideComp[1][1])
        PW.setTitle(title, (col, 0) )


        # Draw dipole outline. 
        [ PW.setLine( X[i, :], Z[i, :], (col, row) ) for i in (0, -1) ]
        [ PW.setLine( X[:, k], Z[:, k], (col, row) ) for k in (0, -1) ]
    # Label the axes. 
    PW.setXlabel( 'X_{GSE} \;\; (R_E)')
    PW.setYlabel( 'Z_{GSE} \;\; (R_E)')
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'E_' + fdrive + 'mHz.png'
      return PW.save(filename)
    else:
      return PW.show()





  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Polarization Plot
  # ---------------------------------------------------------------------------

  def plotPolarizations(self, fdrive):
    # We keep track of run parameters as strings. 
    if not isinstance(fdrive, str):
      fdrive = format(fdrive, '.0f')
    # Electric constant, in mF/m. 
    eps0 = 8.854e-9
    # Magnetic constant, in nH/m. 
    mu0 = 1256.63706
    # Geocentric radii to earth's surface and ionospheric boundary, in Mm. 
    RE, RI = 6.378388, 6.478388
    # Let's look at the energy in the toroidal and poloidal modes as a function
    # of time. Plot on a grid -- m number by row, etc. 
    PW = plotWindow(5, 4, colorbar=False, xPad=2, yPad=3)
    conditions = 'Stormtime'


#    # Set plot supertitle. 
#    PW.setTitle('\mathrm{' + conditions + ' \;\; Toroidal \;\; (Red) \;\; ' +
#                ' and \;\; Poloidal \;\; (Blue), \;\; Dayside \;\; (Solid) ' +
#                ' \;\; and \;\; Nightside \;\;  (Dotted) \;\; Energy}')


    # We'll have the plots all share a Y axis. 
    Umin, Umax = 0, 0
    # Set outer axis labels.  
    PW.setYlabel('\mathrm{Energy \;\; (GJ)}')
    PW.setXlabel('\mathrm{Time \;\; (s)}')
    # Label the top of each column. 
    for col, fdrive in enumerate( ('12', '14', '17', '20', '25') ):
      PW.setTitle( fdrive + '\mathrm{mHz}', pos=(col, 0) )
    # Iterate over the plot cells. 
    for row, azm in enumerate( ('1', '4', '16', '64') ):
      # Label each row by modenumber. 
      PW.setRowLabel('m = ' + azm, row)
      for col, fdrive in enumerate( ('12', '14', '17', '20', '25') ):
        # This plot takes a TON of data. Flush the memory regularly. 
        self.refresh()


#        # Each cell gets toroidal and poloidal, day and night. 
#        for side in ('Day', 'Night'):

        # Each cell gets toroidal and poloidal, day and night. 
        for side in ('Night', ):



          # Set plot supertitle. 
          PW.setTitle('\mathrm{' + side + 'side \;\; ' + conditions + ' \;\; Toroidal \;\; (Red) \;\; ' +
                      ' and \;\; Poloidal \;\; (Blue) \;\; Energy}')



          path = self.getRunPath(azm=azm, side=side, fdrive=fdrive)
          # There's a lot to read in for these plots...
          # Perpendicular electric constant, in mF/m. 
          epsp = eps0*( 1e-10 + self.getArray(path + 'epsp.out') ) # F/m
          # Geocentric radius, in Mm. 
          r = RE*self.getArray(path + 'r.out')
          # Colatitude, in radians. 
          q = self.getArray(path + 'q.out')
          # GSE coordinates, in Mm. 
          X, Z = r*np.sin(q), r*np.cos(q)
          # Time, in seconds. 
          t = self.getArray(path + 't.out')
          # Cosine of the invariant latitude. 
          cosq0 = np.sqrt( 1 - RI*np.sin(q)**2/r )
          # Nonorthogonal coordinates. 
          u1 = -RI/r * np.sin(q)**2
          u3 = RI**2/r**2 * np.cos(q)/cosq0
          # Crunch out the Jacobian, in Mm^3. 
          Jacobian = r**6/RI**3 * cosq0/( 1 + 3*np.cos(q)**2 )
          # Take absolute values so we don't have to worry about real vs imaginary.
          # Each field should have one component that is by far dominant. 
          # Electric fields, in mV/m. 
          Ex = np.abs( self.getArray(path + 'Ex.out') )
          Ey = np.abs( self.getArray(path + 'Ey.out') )
          # Magnetic fields, in nT. 
          Bx = np.abs( self.getArray(path + 'Bx.out') )
          By = np.abs( self.getArray(path + 'By.out') )
          # Note that we don't include the parallel magnetic field. Then we would
          # have to worry about the fact that we're looking at perturbations, not
          # at the zeroth-order field. 
          # Compute the energy density. 
          uP, uT = np.zeros(Ex.shape), np.zeros(Ex.shape)
          for step in range(uP.shape[2]):
            uP[:, :, step] = epsp*Ey[:, :, step]**2 + Bx[:, :, step]**2/mu0
            uT[:, :, step] = epsp*Ex[:, :, step]**2 + By[:, :, step]**2/mu0
          # Now integrate the total energy over the spatial grid. 
          UP = np.zeros( len(t) )
          UT = np.zeros( len(t) )
          for i in range(1, q.shape[0]-1):
            for k in range(1, q.shape[1]-1):
              # The Jacobian maps between a volume in nonorthogonal coordinate 
              # space and a volume in physical space. 
              du1 = u1[i+1, k] - u1[i-1, k]/2
              du3 = u3[i, k+1] - u3[i, k-1]/2
              # Add this area's contribution to all time steps. 
              UP = UP + du1*du3*Jacobian[i, k]*2*np.pi*uP[i, k, :]
              UT = UT + du1*du3*Jacobian[i, k]*2*np.pi*uT[i, k, :]
          # Get the extrema, for the purpose of adjusting the Y limits. 
          Umin = min( Umin, np.min(UP), np.min(UT) )
          Umax = max( Umax, np.max(UP), np.max(UT) )
          # Add these lines to the plot. 
          PW.setLine(t, UP, pos=(col, row), color='b', label='\mathrm{Poloidal}')
          PW.setLine(t, UT, pos=(col, row), color='r', label='\mathrm{Toroidal}')
    # Round up to the nearest power of ten? 
    Umax = 10**( np.ceil( np.log10(Umax) ) )
    PW.setYlimits( [10, Umax] )
    PW.setYlog()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'polarizations.png'
      return PW.save(filename)
    else:
      return PW.show()







  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------- Debug Plot
  # ---------------------------------------------------------------------------

  # For a while, runs were using a sloppy approximation for the current
  # driving -- basically neglecting that the rotation due to the Hall
  # conductivity. The conductivity should be zero-ish where we're driving, so
  # we don't expect this to have been a problem. But let's check. 
  def plotDebug(self):
    # This path uses the sloppy, kludgey current driving. 
    oldPath = '/export/scratch/users/mceachern/2015nov23/T022_004_1_016mHz/'
    # This run uses the proper formulation of the current driving. 
    newPath = '/export/scratch/users/mceachern/runs_20151208_192100/T000/'
    PW = plotWindow(2, 2, colorbar='sym')
    for col, path in enumerate( (oldPath, newPath) ):
      PW.setTitle('\mathrm{Driving \;\; Comparison}')
      # Coordinate arrays. Note that the sign of X is flipped. 
      r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
      t = self.getArray(path + 't.out')
      X, Z = r*np.sin(q), r*np.cos(q)
      which = 'Old \;\; ' if path==oldPath else 'New \;\; '
      PW.setTitle( which + 'B_x', (col, 0) )
      Bx = np.real( self.getArray(path + 'Bx.out') )[:, :, 100]
      PW.setContour( X, Z, Bx, (col, 0) )
      PW.setTitle( which + 'E_y', (col, 1) )
      Ey = np.real( self.getArray(path + 'Ey.out') )[:, :, 100]
      PW.setContour( X, Z, Ey, (col, 1) )
      # Draw dipole outline. 
      [ PW.setLine( X[i, :], Z[i, :] ) for i in (0, -1) ]
      [ PW.setLine( X[:, k], Z[:, k] ) for k in (0, -1) ]
      # Label the axes. 
      PW.setXlabel( 'X_{GSE} \;\; (R_E)')
      PW.setYlabel( 'Z_{GSE} \;\; (R_E)')
    # We don't need to output this. Just show it. 
    return PW.show()

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Poynting Flux Plot
  # ---------------------------------------------------------------------------

  def plotS(self, fdrive):
    # We keep track of run parameters as strings. 
    if not isinstance(fdrive, str):
      fdrive = format(fdrive, '.0f')
    # For the Poynting flux, we want to see a 4x4 grid. Drive frequency is
    # constant across the window (so we will make a plot for each frequency). 
    # Modenumber is constant across each row. The columns are dayside toroidal,
    # nightside toroidal, dayside poloidal, nightside poloidal. Each run just
    # shows just the best frame from a run. All plot cells share a color bar. 


    PW = plotWindow(4, 4, colorbar='sym', yPad=1)


    # Figure out what each column will be. 
    columns = [ (s, p) for p in ('ExBy', 'EyBx') for s in ('Day', 'Night') ]
    # For the moment, we're only looking at storm time, not calm. 
    conditions = 'Stormtime '
    # Set plot supertitle. 
    PW.setTitle('\mathrm{' + conditions + ' \;\; Parallel \;\; Poynting \;\; Flux ' +
                ' \;\; (\\frac{mW}{m^2}) \;\; at \;\; }' + fdrive +
                ' \mathrm{mHz}')
    # Squish the subplots a little bit to make the title stand out. 
    plt.subplots_adjust(top=0.88)
    # Each row is a different modenumber. 
    for row, azm in enumerate( ('1', '4', '16', '64') ):
      # Label each row my modenumber. 
      PW.setRowLabel('m = ' + azm, row)
      # Columns have day/night toroidal/poloidal polarizations. 
      for col, sidePol in enumerate(columns):
        # Find the directory that holds our data. 
        path = self.getRunPath(azm=azm, side=sidePol[0], fdrive=fdrive)
        # Coordinate arrays. Note that the sign of X is flipped. 
        r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
        t = self.getArray(path + 't.out')
        X = -r*np.sin(q) if sidePol[0]=='Day' else r*np.sin(q)
        Z = r*np.cos(q)
        # Grab the Poynting flux and put it on a contour plot. 
        S, step = self.getBestSlice(path + sidePol[1] + '.out')

        # Flip the sign. This is a cross product, after all. 
        if sidePol[1]=='EyBx':
          S = -S

        PW.setContour( X, Z, S, (col, row) )


        # Put together a title for this plot from its parameters. 
        if sidePol[1]=='EyBx':
          title = ('\mathrm{' + sidePol[0] + 'side \;\; Poloidal}')
        else:
          title = ('\mathrm{' + sidePol[0] + 'side \;\; Toroidal}')
        PW.setTitle(title, (col, 0) )


        # Draw dipole outline. 
        [ PW.setLine( X[i, :], Z[i, :], (col, row) ) for i in (0, -1) ]
        [ PW.setLine( X[:, k], Z[:, k], (col, row) ) for k in (0, -1) ]
    # Label the axes. 
    PW.setXlabel( 'X_{GSE} \;\; (R_E)')
    PW.setYlabel( 'Z_{GSE} \;\; (R_E)')
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'S_' + fdrive + 'mHz.png'
      return PW.save(filename)
    else:
      return PW.show()

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------------- Grid Plot
  # ---------------------------------------------------------------------------

  def plotGrid(self):
    # All of the grids should be the same, really. 
    path = self.runs.items[0][0]
    PW = plotWindow(colorbar=False)
    # Coordinate arrays. Note that the sign of X is flipped. 
    r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
    X, Z = r*np.sin(q), r*np.cos(q)
    # Title the window. 
    PW.setTitle('\mathrm{Nonorthogonal \;\; Grid \;\; Spacing}')
    # Label the axes. 
    PW.setXlabel( 'X_{GSE} \;\; (R_E)')
    PW.setYlabel( 'Z_{GSE} \;\; (R_E)')
    # Draw the lines. 
    stride = 2
    [ PW.setLine( X[i, :], Z[i, :] ) for i in range(0, X.shape[0], stride) ]
    [ PW.setLine( X[:, k], Z[:, k] ) for k in range(0, X.shape[1], stride) ]
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'grid.png'
      return PW.save(filename)
    else:
      return PW.show()

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Modenumber Plot
  # ---------------------------------------------------------------------------

  # This plot shows off the difference between modenumber as periodicity and
  # modenumber as localization. 
  def plotM(self, azm):
    # The plot window object isn't actually set up to do polar plots...
    plt.close('all')
    # Set fonts for math formatting. Let's also bump up the font size. 
    rc('font', **{'family':'sans-serif', 'sans-serif':['Helvetica'], 
                  'size':'18'})
    rc('text', usetex=True)
    rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}')
    # Set window. 
    self.fig = plt.figure(figsize=(20., 10.), facecolor='white')
    self.fig.canvas.set_window_title('Tuna Plotter')
    # Create axes. 
    axes = [ plt.subplot(121 + i, projection='polar') for i in range(2) ]
    # Make axes pretty. 
    for ax in axes:
      # Create polar subplot. 
      ax.set_rmax(3.)
      # Remove grid and tick labels. 
      ax.grid(False)
      ax.set_xticklabels( [] )
      ax.set_yticklabels( [] )
      # Draw Earth. 
      day = Wedge( (0, 0), 1, 0, 180, fc='w', transform=ax.transData._b )
      night = Wedge( (0, 0), 1, 180, 0, fc='k', transform=ax.transData._b )
      ax.add_artist(day)
      ax.add_artist(night)
      # Draw dotted lines indicating m. 
      r = np.linspace(1, 3, 1000)
      for section in range(azm):
        q = np.ones( (1000,) )*2*(section+0.5)*np.pi/azm
        ax.plot(q, r, 'k--')
    # Draw periodic wave. 
    q = np.linspace(0, 2*np.pi, 1000)
    E = 2 + 0.5*np.sin(azm*q)
    B = 2 + 0.5*np.cos(azm*q)
    axes[0].plot(q, E, color='r', linewidth=3)
    axes[0].plot(q, B, color='b', linewidth=3)
    # Draw localized wave. 
    q = np.linspace(np.pi/2 - np.pi/azm, np.pi/2 + np.pi/azm, 1000)
    E = 2 + 0.5*np.sin(azm*q)
    B = 2 + 0.5*np.cos(azm*q)
    axes[1].plot(q, E, color='r', linewidth=3)
    axes[1].plot(q, B, color='b', linewidth=3)
    # Set titles. 
    plt.suptitle('$\mathrm{Periodicity \;\; vs \;\; Locality}$', fontsize=24)
    axes[0].set_title('$\mathrm{Periodic \;\; Wave \;\; with \;\;} m = ' +
                      str(azm) + '$')
    axes[1].set_title('$\mathrm{Localized \;\; Wave \;\; with \;\;} m = ' +
                      str(azm) + '$')
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      # Come up with a filename for the output. 
      filename = self.outDir + 'azm.png'
      # Make the directory to hold it, if necessary. 
      savePath = os.path.dirname(filename)
      if not os.path.exists(savePath):
        os.makedirs(savePath)
      # Create the image. 
      plt.savefig(filename)
      print 'Saved plot to ' + filename.replace(os.environ['HOME'], '~')
      return
    else:
      return plt.show()

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Conductivity Plot
  # ---------------------------------------------------------------------------

  def plotConductivity(self):
    PW = plotWindow(colorbar=False)
    for side in ('Day', 'Night'):
      path = self.getRunPath(azm='1', side=side, fdrive='25')
      # Earth radius in km. 
      RE = 6378.388
      alt = self.getArray(path + 'r.out')*RE - RE
      for cond in ('sigP', 'sigH', 'sig0'):
        sig = self.getArray(path + cond + '.out')
        color = 'r' if side=='Day' else 'b'
        if cond=='sigH':
          color = color + '--'
        elif cond=='sig0':
          color = color + ':'
        label = '$\mathrm{' + side + 'side} \;\; \sigma_' + cond[-1] + '$'
        PW.setLine(sig[0, :], alt[0, :], color=color, label=label)
    PW.setXlog()
#    PW.setXlimits( (1e-7, 1e3) )
    PW.setXlimits( (1e-7, 1e8) )
    PW.setYlimits( (0, 2000) )
    PW.setTitle('\mathrm{Conductivity \;\; Profiles}')
    PW.setXlabel( '\mathrm{Conductivity \;\; (\\frac{S}{m})}')
    PW.setYlabel( '\mathrm{Altitude \;\; (km)}')
    PW.setLegend()
    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      filename = self.outDir + 'sig.png'
      PW.save(filename)
      return
    else:
      return PW.show()




  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Alfven Speed Plot
  # ---------------------------------------------------------------------------

  def plotV(self, path):

    # Initialize plot window object. 
    PW = plotWindow(2, 2, colorbar=True)

    # Figure out run parameters. 
    model = '???'
    for line in open(path + 'params.in', 'r').readlines():
      if 'model' in line.lower():
        model = ('Dayside' if '1' in line else 'Nightside') + ' \;\; Stormtime'
    # Use those to set the window's supertitle. 
    PW.setTitle(model + '\;\; Dimensionless \;\; Ionosphere \;\; Parameters')
    # Label the axes. 
    PW.setXlabel( '\mathrm{Invariant \;\; Latitude} \;\; (^\circ)')
    PW.setYlabel( '\mathrm{Altitude} \;\; (\mathrm{km})')
    # Grab coordinates. 
    r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
    # Invariant latitude folds the hemispheres together. To avoid overlapping
    # data, let's drop the southern hemisphere values completely. 
    nHalf = r.shape[1]/2 + 1
    # Compute altitude and invariant latitude for the plot axes. 
    RE = 6378.388 # km
    altitude = (r[:, :nHalf] - 1)*RE # km
    invLat = np.empty(altitude.shape)
    for k in range(nHalf):
      invLat[:, k] = 90 - (180/np.pi)*q[:, 0] # degrees. 
    # Put this on a log scale to make better use of space. 
    PW.setYlog()
    # Grab data to plot. Let's make sure everything is in SI units. Also we're
    # careful about nonpositive values since we'll be dividing a lot. 
    vA = 1e3*self.getArray(path + 'vA.out')[:, :nHalf] # m/s
    sig0 = 1e-10 + np.abs( self.getArray(path + 'sig0.out')[:, :nHalf] ) # S/m
    sigH = 1e-10 + np.abs( self.getArray(path + 'sigH.out')[:, :nHalf] ) # S/m
    sigP = 1e-10 + np.abs( self.getArray(path + 'sigP.out')[:, :nHalf] ) # S/m
    eps0 = 8.854e-12 # F/m
    epsPerp = eps0*self.getArray(path + 'epsp.out')[:, :nHalf] # F/m
    c = 3e8 # m/s

    PW.setTitle( 'v_A / c', (0, 0) )
    PW.setContour( invLat, altitude, vA/c, (0, 0) )

    PW.setTitle( '\sigma_H / \sigma_P', (1, 0) )
    PW.setContour( invLat, altitude, sigH/sigP, (1, 0) )

    PW.setTitle( 'c \epsilon_\parallel / \sigma_0 R_E', (0, 1) )
    PW.setContour( invLat, altitude, c*epsPerp/(sig0*RE), (0, 1) )

    PW.setTitle( '\epsilon_\\bot v_A / \sigma_H R_E', (1, 1) )
    PW.setContour( invLat, altitude, (epsPerp*vA)/(sigH*RE), (1, 1) )

    # Show the plot. 
    PW.show()
    return'''






# #############################################################################
# ########################################################## Plot Window Object
# #############################################################################

# The plot window class is meant to be a general-purpose wrapper around 
# Pyplot. It creates a plot window with a grid of subplots, and makes it 
# convenient to add data to them. 
class plotWindow:

  axes = None
  colorAxis = None
  cells = None
  nColors = 8
  pos = None

  # ---------------------------------------------------------------------------
  # ----------------------------- Initialize Plot Window and Space Out Subplots
  # ---------------------------------------------------------------------------

  # At this point we need to know the number of rows, the number of columns, 
  # how they should be spaced, and if there will be a color bar. 
  def __init__(self, nCols=1, nRows=1, colorbar=None, xSize=10, xPad=1, 
               ySize=10, yPad=5, cSize=2, cPad=1, **kargs):
    # Options are 'linear', 'log', and 'loglog'. Boolean True gives a linear
    # scale, while False and None give no color bar. For now, we just care if
    # there will be a color bar at all, so we can make room for it. 
    self.colorbar = colorbar
    # Make sure there isn't anything lingering from a previous plot. 
    plt.close('all')
    # Set fonts for math formatting. Let's also bump up the font size. 
    rc('font', **{'family':'sans-serif', 'sans-serif':['Helvetica'], 
                  'size':'18'})
    rc('text', usetex=True)
    rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}')
    # Figure proportions are determined by the number of subplots. 
#    self.fig = plt.figure(figsize=(15., nRows*10./nCols), facecolor='white')
    self.fig = plt.figure(figsize=(nCols*15./nRows, 10.), facecolor='white')
    self.fig.canvas.set_window_title('Tuna Plotter')
    # The default subplot spacing sometimes overlaps labels, etc. We instead
    # use GridSpec. The subplots (and color bar, if any) are spaced out in
    # terms of a grid of equally-sized tiles. If we want a color bar, we need
    # extra columns to hold it. 
    if self.colorbar:
      cSize, cPad = 2, 1
      tiles = gridspec.GridSpec( (ySize + yPad)*nRows - yPad, 
                                 (xSize + xPad)*nCols - xPad + cPad + cSize )
      # Grab the tiles for the color bar. 
      self.colorAxis = plt.subplot( tiles[:, -cSize:] )
    # If we don't need a color bar, partition the window into a different
    # number of tiles. 
    else:
      tiles = gridspec.GridSpec( (ySize + yPad)*nRows - yPad, 
                                 (xSize + xPad)*nCols - xPad )
    # Create an array of axes, one for each subplot. 
    self.axes = np.empty( (nCols, nRows), dtype=object)
    for c in range(nCols):
      for r in range(nRows):
        xPos, yPos = (xSize + xPad)*c, (ySize + yPad)*r
        xSlice, ySlice = slice(xPos, xPos + xSize), slice(yPos, yPos + ySize)
        self.axes[c, r] = plt.subplot( tiles[ySlice, xSlice] )
    # Create a corresponding array of plot cells. 
    self.cells = np.empty( (nCols, nRows), dtype=object)
    for c in range(nCols):
      for r in range(nRows):
        self.cells[c, r] = plotCell( self.axes[c, r] )
    # If there were any other keyword args in the initialization call, they're
    # presumably parameters to be set for the whole window. 
    for key, val in kargs.items():
      self.setParam( **{key:val} )
    # The plot window is now ready for us to add some data to it! 
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Set Plot Parameters
  # ---------------------------------------------------------------------------

  # Both the plot window and the plot cell objects have setParam methods. The
  # plot window method is used to set whole-window parameters, while the cell
  # method is for affecting a single cell. The __getitem__ method allows users
  # to easily specify which they want. 
  def __getitem__(self, index):
    # To access a row or column of cells, we still want to talk to the window. 
    if isinstance(index, int):
      self.pos = index
      return self
    # We can also access a single cell. 
    return self.cells[index]

  # Plot labels, coordinates, and other parameters are all handled through this
  # function using keyword arguments. Parameters can be applied to the whole
  # window or to just a single cell. 
  def setParam(self, pos=None, **kargs):
    # Use caps-insensitive keys. 
    params = dict( (key.lower(), kargs[key]) for key in kargs )
    # Check if we're using plotWindow[pos].setParam for a row or column. 
    if self.pos is not None:
      # There shouldn't be many things to scroll through here...
      for key, val in params.items():
        # Row label, to the right of the row of subplots. 
        if key=='rowlabel':
          # Plots without very many columns may need to be rearranged. 
          if self.cells.shape[0]==2:
            plt.subplots_adjust(left=0.25)
          # Figure out where the leftmost axes are, and write to the left. 
          loc = self.axes[0, self.pos].get_position()
          plt.figtext(0.3*loc.x0, loc.y0 + 0.5*loc.height, '$' + val + '$',
                      horizontalalignment='center', 
                      verticalalignment='center')
        # A column label is just a title on the top row. This may seem like
        # cheating, but why would we ever want both? 
        elif key=='collabel':
          self.cells[self.pos, 0].setParam(title=val)
        # Notify for anything that needs updating. 
        else:
          print 'UNKNOWN KEY WITH INTEGER POS: ', key
      # The effect only lasts for one call, obviously. 
      self.pos = None
    # If not calling through __getitem__, apply the parameters to the whole
    # plot window. 
    else:
      # Usually, applying something to the whole plot means applying it to each
      # subplot individually. Here, we check for special cases. 
      for key, val in params.items():
        # Number of colors. 
        if key=='ncolors':
          self.nColors = val
        # A positionless title is the supertitle. 
        elif key=='title':
          plt.suptitle('$' + val + '$', fontsize=24)
          # Squish the subplots a little bit to make the title stand out. 
          plt.subplots_adjust(top=0.88)
        # Only put shared x labels on the bottom edge. 
        elif key=='xlabel':
          for cell in self.cells[:, -1]:
            cell.setParam( **{key:val} )
        # Only put shared y labels on the left edge. 
        elif key=='ylabel':
          for cell in self.cells[0, :]:
            cell.setParam( **{key:val} )
        # Everything else just goes to each cell. 
        else:
          for column in self.cells:
            for cell in column:
              cell.setParam( **{key:val} )
    return

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Add Data to Plot
  # ---------------------------------------------------------------------------

  # Add a line to the plot. This can't be (easily) handled with setParam
  # because in addition to the line itself, the line has several other options. 
  def setLine(self, X, Y, pos=None, color='k', label=None):
    if pos is None:
      for column in self.cells:
        [ cell.setLine(X, Y, color=color, label=label) for cell in column ]
    else:
      self.cells[pos].setLine(X, Y, color=color, label=label)
    return

  # If no position is given, all plots will draw the contour. This could be set
  # with the other parameters above, but it's possible that we'll want to add
  # other options to the contour call. 
  def setContour(self, Z, pos=None):
    if pos is None:
      [ cell.setContour(Z) for column in self.cells for cell in column ]
    else:
      self.cells[pos].setContour(Z)
    return

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Render Plot Window
  # ---------------------------------------------------------------------------

  # Once all of the contours are loaded, we can figure out the color levels. 
  def render(self, filename=None):
    # Get the maximum of all contour maxima. 
    vmax = nax( cell.getMax() for column in self.cells for cell in column )
    # The plot colors object constructor creates the color bar, then returns
    # a dictionary of keyword parameters. 
    PC = plotColors(vmax, self.colorAxis, self.colorbar, nColors=self.nColors)
    # We send those parameters to each cell for use in their contourf calls. 
    [ cell.render(**PC) for column in self.cells for cell in column ]
    # If given a filename, save the plot window as an image. 
    if filename is not None:
      # If the output directory doesn't exist, make it. This ensures that we
      # don't create an output directory unless we're actually making output. 
      savePath = os.path.dirname(filename)
      if not os.path.exists(savePath):
        os.makedirs(savePath)
      # Create the image. 
      plt.savefig(filename)
      print 'Saved plot to ' + filename.replace(os.environ['HOME'], '~')
    # Otherwise, display it. 
    else:
      plt.show()
    return

# #############################################################################
# ############################################################ Plot Cell Object
# #############################################################################

class plotCell:

  # A place to hold the coordinates, and the data, in case of a contour. 
  X, Y, Z = None, None, None
  # If any lines are to be drawn on the plot, they are stored here. 
  lines = ()
  # By default, do not use a legend or an outline.
  legend = False
  outline = False

  # ---------------------------------------------------------------------------
  # ------------------------------------- Initialize Plot Cell from Axis Object
  # ---------------------------------------------------------------------------

  def __init__(self, ax):
    self.ax = ax
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Set Plot Parameters
  # ---------------------------------------------------------------------------

  # All parameter handling -- coordinates, labels, etc -- is handled here using
  # keyword arguments. 
  def setParam(self, **kargs):
    for key, val in kargs.items():
      # Caps insensitivity. 
      key = key.lower()
      # Outline of grid (in case of dipole geometry or whatever).  
      if key=='outline':
        self.outline = val
      # Title. 
      elif key=='title':
        self.ax.set_title('$' + val + '$')
      # Array of horizontal coordinate values. 
      elif key=='x':
        self.X = val
      # Horizontal axis label. 
      elif key=='xlabel':
        self.ax.set_xlabel('$' + val + '$')
      # Horizontal axis bounds. 
      elif key.startswith('xlim'):
        self.ax.set_xlim(val)
      # Horizontal axis log scale. 
      elif key=='xlog' and val==True:
        self.ax.set_xscale('log')
      # Horizontal axis tick locations. 
      elif key=='xticks':
        self.ax.set_xticks(val)
      # Horizontal axis tick labels. 
      elif key=='xticklabels':
        self.ax.set_xticklabels(val)
      # Array of vertical coordinate values. 
      elif key=='y':
        self.Y = val
      # Vertical axis label. 
      elif key=='ylabel':
        self.ax.set_ylabel('$' + val + '$')
      # Vertical axis bounds. 
      elif key.startswith('ylim'):
        self.ax.set_ylim(val)
      # Vertical axis log scale. 
      elif key=='ylog' and val==True:
        self.ax.set_yscale('log')
      # Vertical axis tick locations. 
      elif key=='yticks':
        self.ax.set_yticks(val)
      # Vertical axis tick labels. 
      elif key=='yticklabels':
        self.ax.set_yticklabels(val)

      else:
        print 'UNKNOWN KEY: ', key

    return

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Add Data to Plot
  # ---------------------------------------------------------------------------

  def setContour(self, Z):
    self.Z = Z
    return

  def setLine(self, X, Y, color='k', label=None):
    self.lines = self.lines + ( (X, Y, color, label), )
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------ Normalize Data
  # ---------------------------------------------------------------------------

  # All cells share a color bar. We need a minimum if we're using a log scale. 
  def getMin(self):
    return None if self.Z is None else np.min(self.Z)

  # We always need the maximum. 
  def getMax(self):
    return None if self.Z is None else np.max( np.abs(self.Z) )

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # For line plots, all cells should share a vertical scale. But we want to
  # base this on the yticks -- which don't exist until after we make the plot
  # commands -- rather than the actual line values. We should probably just
  # plot lines as they come in. 
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  def getYmin(self):
    return nin( np.min( line[1] ) for line in self.lines )

  def getYmax(self):
    return nax( np.max( line[1] ) for line in self.lines )

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------- Render Cell
  # ---------------------------------------------------------------------------

  def render(self, **colorParams):
    # If this cell has a contour, lay that down first. 
    if self.Z is not None:
      self.ax.contourf(self.X, self.Y, self.Z, **colorParams)
    # On top of that, draw any lines. 
    for line in self.lines:
      label = None if line[3] is None else '$' + line[3] + '$'
      self.ax.plot( line[0], line[1], line[2], label=label )
    # Draw the outline of the grid, usually for dipole plots. 
    if self.outline:
      [ self.ax.plot(self.X[i, :], self.Y[i, :], color='k')for i in (0, -1) ]
      [ self.ax.plot(self.X[:, k], self.Y[:, k], color='k')for k in (0, -1) ]
    # Draw the legend. 
    if self.legend:
      self.ax.legend(loc='best')

    # Remove ticks without removing tick labels. 
#    self.ax.tick_params( width=0 )
#    # Remove plot frames. 
#    for pos in ('top', 'bottom', 'left', 'right'):
#      self.ax.spines[pos].set_visible(False)
    # Only put ticks on the sides with numbers. 
    self.ax.get_xaxis().tick_bottom()
    self.ax.get_yaxis().tick_left()
    # Only put numbers on sides with labels. 
    if not self.ax.get_xlabel():
      self.ax.set_xticklabels( [] )
#      self.ax.set_xticks( [] )
#      self.ax.get_xaxis().set_visible(False)
    if not self.ax.get_ylabel():
      self.ax.set_yticklabels( [] )
#      self.ax.set_yticks( [] )
#      self.ax.get_yaxis().set_visible(False)
    return

# #############################################################################
# ######################################################### Plot Colors Handler
# #############################################################################

# This class figures out the color levels and ticks to be used by all of the
# contour plots. It draws the color bar, then serves as a keyword dictionary
# for the contourf calls. 

class plotColors(dict):

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Initialize Colors
  # ---------------------------------------------------------------------------

  def __init__(self, vmax, ax, colorbar, nColors):
    # Some plots don't have contours. 
    if not vmax or not ax or not colorbar:
      return dict.__init__(self, {})
    # Store the data scale so that it can be hard-wired into our normalization
    # functions. We don't want to pass vmax all over the place. 
    self.vmax = vmax
    self.colorbar = colorbar
    self.nColors = nColors
    self.nTicks = nColors - 1
    # Assemble the keyword parameters in a temporary dictionary. We'll then use
    # the dictionary constructor to build this object based on it. 
    temp = {}
    # Determine location of contour color levels and color bar ticks. 
    if self.colorbar=='log':
      temp['ticks'], temp['levels'] = self.logTicksLevels()
      temp['norm'] = LogNorm()
    elif self.colorbar=='sym':
      temp['ticks'], temp['levels'] = self.symTicksLevels()
      temp['norm'] = Normalize()
    else:
      temp['ticks'], temp['levels'] = self.linTicksLevels()
      temp['norm'] = Normalize()
    # Rework the color map to match the normalization of our ticks and levels. 
    temp['cmap'] = self.getCmap()
    # Draw the color bar. 
    self.setColorbar(ax, **temp)
    # Become a dictionary of color parameters to be used by the contour plots. 
    return dict.__init__(self, temp)

  # ---------------------------------------------------------------------------
  # ----------------------------------- Tick Locations and Contour Color Levels
  # ---------------------------------------------------------------------------

  def linTicksLevels(self):
    ticks = np.linspace( -self.vmax, self.vmax, self.nTicks)
    levels = np.linspace(-self.vmax, self.vmax, self.nColors)
    return ticks, levels

  def logTicksLevels(self):
    # One tick at each order of magnitude. 
    power = int( np.floor( np.log10(self.vmax) ) )
    self.vmin = self.vmax/10**self.nTicks
    ticks = [ 10**(power - i) for i in range(self.nTicks) ]
    logMin, logMax = np.log10(self.vmin), np.log10(self.vmax)
    levels = np.logspace(logMin, logMax, self.nColors)
    return ticks, levels

  def symTicksLevels(self):
    # A tick at zero, then one per order of magnitude. 
    nOrders = (self.nTicks - 1)/2
    power = int( np.floor( np.log10(self.vmax) ) )
    posTicks = [ 10**(power - i) for i in range(nOrders) ]
    # For uniform tick spacing, the log cutoff needs to be a factor of ten
    # smaller than the lowest positive tick. 
    self.vmin = min(posTicks)/10.
    ticks = sorted( posTicks + [0] + [ -t for t in posTicks ] )
    # We figure out color levels by spacing them evenly on the unit interval,
    # then mapping the unit interval to the symlog scale. 
    levels = [ self.symNorm(x) for x in np.linspace(0, 1, self.nColors) ]
    return ticks, levels

  # ---------------------------------------------------------------------------
  # ----------------------------------------------- Data Interval Normalization
  # ---------------------------------------------------------------------------

  # Map from the unit interval to the data scale via linear scale. 
  def linNorm(self, x):
    return self.vmax*(2*x - 1)

  # Map from the data scale to the unit interval via linear scale. 
  def linMron(self, x):
    return 0.5 + 0.5*x/self.vmax

  # Map from the unit interval to the data scale via log scale. 
  def logNorm(self, x):
    return self.vmin*(self.vmax/self.vmin)**x

  # Map from the log scaled data scale to the unit interval. 
  def logMron(self, x):
    return np.log10(x/self.vmin)/np.log10(self.vmax/self.vmin)

  # Map from the unit interval to the data scale via symmetric log scale. 
  def symNorm(self, x):
    if x>0.5:
      return self.vmin*(self.vmax/self.vmin)**(2*x - 1)
    elif x<0.5:
      return -self.vmin*(self.vmax/self.vmin)**(1 - 2*x)
    else:
      return 0

  # Map from the symmetric log scaled data scale to the unit interval. 
  def symMron(self, x):
    if x>self.vmin:
      return 0.5 + 0.5*np.log10(x/self.vmin)/np.log10(self.vmax/self.vmin)
    elif x<-self.vmin:
      return 0.5 - 0.5*np.log10(-x/self.vmin)/np.log10(self.vmax/self.vmin)
    else:
      return 0.5

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------------- Color Map
  # ---------------------------------------------------------------------------

  # A color map is really just a handful of RGB codes defined on the unit
  # interval. To show up nicely, we need to renormalize that unit interval to
  # match the normalization of our ticks and color levels. 
  def getCmap(self):
    # Figure out the unit interval renormalization to use. 
    if self.colorbar=='log':

      # Kinda kludgey. See setColorbar for explanation. 
#      norm = self.logNorm
#      return plt.get_cmap('seismic')
      return None

    elif self.colorbar=='sym':
      norm = self.symNorm
    else:
      norm = self.linNorm
    # Get a fine sampling of the color map on the unit interval. 
    N = 1000
    unitInterval = [ i/(N - 1.) for i in range(N) ]
    cmap = plt.get_cmap('seismic')
    rgb = [ cmap( unitInterval[i] ) for i in range(N) ]
    # Renormalize the unit interval. Make sure we get the end points right. 
    newInterval = [ self.linMron( norm(u) ) for u in unitInterval ]
    newInterval[0], newInterval[-1] = 0., 1.
    # The color dict contains three channels. Each of them is a list of tuples,
    # (u, rgbLeft, rgbRight). Between u values, colors are interpolated
    # linearly. Approaching u from the right, the color should approach
    # rgbRight, and vice versa. Since our color map is smooth and
    # finely-resolved, we don't bother to distinguish.  
    red = [ (newInterval[i], rgb[i][0], rgb[i][0]) for i in range(N) ]
    grn = [ (newInterval[i], rgb[i][1], rgb[i][1]) for i in range(N) ]
    blu = [ (newInterval[i], rgb[i][2], rgb[i][2]) for i in range(N) ]
    # Return a LinearSegmentedColormap built from the color dictionary. We use
    # TONS of samples because the symmetric log normalization devotes very
    # little of the unit interval to zero... but that it the most important bin
    # to have sharply defined. 
    return LSC('myMap', {'red':red, 'green':grn, 'blue':blu}, 1000000)

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------- Set Color Bar
  # ---------------------------------------------------------------------------

  # Without SymLogNorm, we can't very well use the built-in color bar
  # functionality. Instead, we make our own: a tall, narrow contour plot. 
  def setColorbar(self, ax, **colorParams):
    # Unit interval axes. 
    X, Y = np.empty( (2, 1000) ), np.empty( (2, 1000) )
    for i in range(2):
      X[i, :] = i
      Y[i, :] = np.linspace(0, 1, 1000)
    # The contour values are just the Y axis, mapped to the data scale via
    # linear, log, or symmetric log normalizer. We'll also need the inverse
    # normalizers, since we have the tick values in the data scale, and we map
    # them to the Y axis (which is on the unit interval). And the formatters,
    # to make our ticks look pretty in LaTeX. 
    if self.colorbar=='log':
      # This is kludgey right now. Sorry. We can't use a real color bar for the
      # symmetric norm plot, since since SymLogNorm isn't defined. But we can't
      # use a renormalized color map for the log plot due to sampling
      # constraints. This is the odd case out right now. 
      norm, mron, fmt = self.logNorm, self.logMron, self.logFormatter
      ColorbarBase(ax, boundaries=colorParams['levels'],
                   ticks=colorParams['ticks'], norm=colorParams['norm'],
                   cmap=colorParams['cmap'])
      ax.set_yticklabels( [ fmt(t) for t in colorParams['ticks'] ] )
      return
    elif self.colorbar=='sym':
      norm, mron, fmt = self.symNorm, self.symMron, self.symFormatter
    else:
      norm, mron, fmt = self.linNorm, self.linMron, self.linFormatter
    # Draw the contour. 
    Z = np.vectorize(norm)(Y)
    ax.contourf( X, Y, Z, **colorParams)
    # Place the ticks appropriately on the unit interval (Y axis). 
    ax.set_yticks( [ mron(t) for t in colorParams['ticks'] ] )
    # Format tick names nicely. 
    ax.set_yticklabels( [ fmt(t) for t in colorParams['ticks'] ] )
    # Put the color bar ticks on the right, get rid of the ticks on the bottom,
    # and hide the little notches in the color bar. 
    ax.yaxis.tick_right()
    ax.set_xticks( [] )
    ax.tick_params( width=0 )
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------ Tick Name Formatting
  # ---------------------------------------------------------------------------

  # We have to format the ticks for two reasons. First, because the color bar Y
  # axis is on the unit interval, not the data scale (and may not be normalized
  # properly). Second, because that's how we make sure to get dollar signs in
  # there so LaTeX handles the font rendering. 

  def linFormatter(self, x):
    # Zero is always zero. 
    if x==0:
      return '$0$'
    # If our numbers are around order unity, show floats to two significant
    # digits. Otherwise, use scientific notation. 
    elif 1e-3<self.vmax<1e3:
      digs = int( format(abs(x), '.1e').split('e')[0].replace('.', '') )
      power = int( format(abs(x), '.1e').split('e')[1] ) - 1
      return '$' + ( '-' if x<0 else '+' ) + str(digs*10**power) + '$'
    else:
      # Cast the number in scientific notation. 
      s = format(x, '.1e').replace('e', ' \\cdot 10^{') + '}'
      # If the number is positive, throw a plus sign on there. 
      s = '+ ' + s if x>0 else s
      # Before returning, get rid of any extra digits in the exponent. 
      return '$' + s.replace('+0', '').replace('-0', '-') + '$'

  def logFormatter(self, x):
    # Zero is always zero. 
    if x==0:
      return '$0$'
    # Otherwise, just keep the power of ten. 
    return '$ 10^{' + format(np.log10(x), '.0f') + '}$'

  def symFormatter(self, x):
    # Zero is always zero. 
    if x==0:
      return '$0$'
    # Otherwise, just keep the sign and the power of ten. 
    power = format(np.log10( np.abs(x) ), '.0f')
    return '$ ' + ( '-' if x<0 else '+' ) + ' 10^{' + power + '}$'

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Turns a list of numbers (1, 2, 3) into the string '1x2x3'. 
def by(x):
  return str( x[0] ) + 'x' + by( x[1:] ) if len(x)>1 else str( x[0] )

# Left-justify and/or trim a value to fit into a fixed-width column. 
def col(x, width=10):
  return str(x).ljust(width)[:width-1] + ' '

# Convert a string to a complex or float. 
def com(x):
  if ',' not in x:
    return float(x)
  else:
    # Shave off the parentheses then split into real and imaginary parts. 
    re, im = x[1:-1].split(',')
    return (float(re) + float(im)*1j)

# Get a listing of directories at a given path. Remember that we always use
# absolute paths, and that directory paths end in a slash. 
def dirs(path=None):
  # By default, use the current location. 
  p = path if path is not None else os.getcwd() + '/'
  return sorted( p + x + '/' for x in os.listdir(p) if os.path.isdir(p + x) )

# Grab all files in the given directory. Remember that we always use absolute
# paths, and that directory names always end with a slash. 
def files(path=None, end=''):
  p = path if path is not None else os.getcwd()
  # Grab all files in the directory. 
  f = [ p + x for x in os.listdir(p) if os.path.isfile(p + x) ]
  # Only return the ones that end with our desired extension. 
  return sorted( x for x in f if x.endswith(end) )

# Safely ask for the first element of a possibly-empty list or generator. 
def first(x):
  lx = list(x)
  return lx[0] if len(lx)>0 else None

# Safely ask for the minimum of a possibly-empty list or generator. 
def nin(x):
  lx = list(x)
  return min(lx) if len(lx)>0 else None

# Safely ask for the maximum of a possibly-empty list or generator. 
def nax(x):
  lx = list(x)
  return max(lx) if len(lx)>0 else None

# Timestamp for labeling output. 
def now():
  return ( znt(lt().tm_year, 4) + znt(lt().tm_mon, 2) + znt(lt().tm_mday, 2) +
           '_' + znt(lt().tm_hour, 2) + znt(lt().tm_min, 2) +
           znt(lt().tm_sec, 2) )

# Turn a string into a float or integer. 
def num(x):
  return int( float(x) ) if float(x)==int( float(x) ) else float(x)

# Returns a list of non-blank lines from a file. 
def read(filename):
  if not os.path.isfile(filename):
    return []
  else:
    return [ x.strip() for x in open(filename, 'r').readlines() if x.strip() ]

# Turns the number 3 into '003'. If given a float, truncate it. 
def znt(x, width=0):
  return str( int(x) ).zfill(width)

# =============================================================================
# ========================================================= Shell Input Parsing
# =============================================================================

# Grab the arguments from the shell and partition them into the names of plots 
# to make, the paths where those plots should be created, and the flags that
# will affect how those plots are displayed. 
def getArgs():
  flags, paths, names = [], [], []
  # The first argument is the script call. Skip it. 
  for arg in argv[1:]:
    # Anything that starts with a dash is a flag. Use lower case. 
    if arg.startswith('-'):
      flags.append( arg.lower() )
    # Arguments that are directories probably hold data.  
    elif os.path.isdir(arg):
      # Look recursively through the contents of each directory. 
      for root, dirs, files in os.walk(arg):
        # Keep any directory containing a Tuna run. Always use absolute paths,
        # and always end directory names with a slash. 
        if 'params.in' in files:
          paths.append( os.path.abspath(root) + '/' )
    # Plot names are a single character long, and can accept a comma-separated
    # list of arguments (no whitespace) after a colon. 
    else:
      # Plot name with no arguments. 
      if len(arg)==1:
        names.append( ( arg.lower(), () ) )
      # Plot name with arguments. 
      elif arg[1]==':':
        name = arg[0].lower()
        argList = [ num(x) for x in arg[2:].split(',') ]
        names.append( (name, argList) )
  return flags, paths, dict(names)

# =============================================================================
# ============================================================ Data File Access
# =============================================================================

def readArray(filename, indent='\t'):
  # The out prefix is just an older convention for dat files, Fortran output.
  # We ultimately want to use the Python data format, pickles. 
  name = filename[ :filename.rfind('.') ]
  datname, outname, pklname = name + '.dat', name + '.out', name + '.pkl'

  # Rename jz.out to Jz.dat (change in capitalization). 
  outname = outname.replace('Jz', 'jz')

  # If we're looking at an out postfix, move it to dat. 
  if os.path.isfile(outname) and not os.path.exists(datname):
    os.rename(outname, datname)
    print ( indent + 'WARNING: Renamed ' + basename(outname) + ' to ' +
            basename(datname) )

#  # Get the expected array dimensions from the dat file. Note that this seems
#  # to slow things down significantly. Rather than perform this check, let's
#  # just report the dimensions. 
#  if os.path.exists(datname):
#    # Grab the first line of the dat file, while contains the dimensions. 
#    with open(datname, 'r') as f:
#      dimLine = f.readline()
#    # The first line is the array dimensions. 
#    dims = [ int(x) for x in dimLine.split() ]
#  else:
#    print indent + 'WARNING: ' + basename(datname) + ' not found. '
#    dims = (0,)

  # If a pickle is available, read that instead of parsing the Fortran output. 
  if os.path.isfile(pklname):
    print indent + 'Reading ' + col( basename(pklname), 10) + ' ... ',
    stdout.flush()
    # Note how long it takes to read the file. 
    start = time()
    with open(pklname, 'rb') as handle:
      arr = pickle.load(handle)
    print format(time() - start, '5.1f') + 's' + ' ... ' + by(arr.shape)
    return arr
  # If the pickle doesn't exist yet, parse the Fortran output. 
  elif os.path.isfile(datname):
    print indent + 'Reading ' + col( basename(datname), 10) + ' ... ',
    stdout.flush()
    # Note how long it takes to read the file. 
    start = time()
    # Grab the data as a list of strings. 
    with open(datname, 'r') as f:
      arrayLines = f.readlines()
    # The first line is the dimensions of the array. 
    dims = [ int(x) for x in arrayLines.pop(0).split() ]
    # Assemble a one-dimensional array large enough to hold all of the values.
    # (This is much faster than appending as we go.) This means figuring out if
    # we want reals or complexes. We check by looking for a comma, since
    # complex values are listed as ordered pairs. 
    if len(arrayLines)>0 and ',' in arrayLines[0]:
      dtype = np.complex
    else:
      dtype = np.float
    # Create the empty array. 
    nVals = np.prod(dims)
    vals = np.empty(nVals, dtype=dtype)
    # Now fill the array with values one at a time. Stop when it's full, or
    # when we run out of values. 
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
    # Reshape and transpose the array. Fortran and Python have opposite
    # indexing conventions. 
    arr = np.transpose( np.reshape( vals, dims[::-1] ) )
    # Check the array dimensions. If it's not full, the run may have crashed. 
    actualDims = dims[:-1] + [ np.int( i/np.prod(dims[:-1]) ) ]
    # Keep only the time steps that actually got printed. 
    if dims!=actualDims:
      arr = arr[ ..., :actualDims[-1] ]
    # Report how long this all took. 
    print format(time() - start, '5.1f') + 's'
    # Dump a pickle for next time. 
    print indent + 'Creating ' + basename(pklname) + ' ... ',
    stdout.flush()
    start = time()
    with open(pklname, 'wb') as handle:
      pickle.dump(arr, handle, protocol=-1)
    print format(time() - start, '5.1f') + 's'
    # Check if the data dimensions are consistent with the dimensions in the
    # header of the dat file. 
    if by(dims)!=by(arr.shape):
      print ( indent + 'WARNING: Expected ' + by(dims) + ' but found ' +
              by(arr.shape) )
    # When we read in lots of arrays in quick succession, Python's automatic
    # garbage collection gets overwhelmed. We solve that by manually forcing
    # garbage collection. 
    gc.collect()
    # Finally, return the array. 
    return arr
  # If the pickle doesn't exist and there's no Fortran output, return nothing. 
  else:
    print indent + 'WARNING: ' + datname + ' not found. '

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


