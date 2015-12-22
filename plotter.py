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
from sys import argv
from time import localtime as lt, sleep

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():
  # Parse input from the terminal. 
  flags, paths, names = getArgs()

  print 'got flags: ', flags

  print 'got names: ', names

  TP = tunaPlotter(flags, paths)

#  TP.plotDebug()

#  TP.plotGrid()

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

  if 'f' in names:
    TP.plotFactor()

  if 'g' in names:
    TP.plotGrid()

  if 'm' in names:
    TP.plotM(4)

  if 'p' in names:
    TP.plotPolarizations(20)

  if 's' in names:
    for fdrive in (12, 14, 17, 20, 25):
      TP.plotS(fdrive)

#  for path in paths:

#    TP.plotV(path)

  return

# #############################################################################
# ######################################################### Plot Window Wrapper
# #############################################################################

# The window wrapper class includes utilities which are specific to the Tuna
# plotter, such as interfacing with the output data. 
class tunaPlotter:

  # To avoid spending time reading in the same file more than once, data input
  # and access is centralized. Data is deleted when we change directories. 
  data = {}
  path = None

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Tuna Plotter
  # ---------------------------------------------------------------------------

  def __init__(self, flags, paths):
    # Keep track of flags from the terminal. 
    self.flags = flags
    # If we'll be saving output, figure out where to put it. 
    self.outDir = '/home/user1/mceachern/Desktop/plots/plots_' + now() + '/'
    # The plotter needs to know where to find all the data we'll be plotting. 
    self.paths = paths
    # Let's also match up those paths with run parameters. We do this by
    # opening up the params file for each path and taking a peek inside. 
    self.runs = dict( (p, {'azm':'?', 'tau':'?', 'side':'?'} ) for p in paths )
    for p in paths:
      # Grab the entire contents of the parameters input file. 
      params = '\n'.join( open(p + 'params.in', 'r').readlines() )
      # Split out the azimuthal modenumber. 
      self.runs[p]['azm'] = params.split('azm')[1].split('\n')[0].strip(' =')
      # Get the drive frequency. 
      fdrive = float( params.split('fdrive')[1].split('\n')[0].strip(' =') )
      self.runs[p]['fdrive'] = format(1000*fdrive, '.0f')
      # Get the model number and spell out what it means. 
      model = params.split('model')[1].split('\n')[0].strip(' =')
      self.runs[p]['side'] = 'Day' if model=='1' else 'Night'
    return

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------- Data Access
  # ---------------------------------------------------------------------------

  # Given a set of run parameters, find the desired run path. Note that all run
  # parameters are stored as strings! 
  def getRunPath(self, **kargs):
    # Start out with a list of all paths. 
    candidates = [ p for p in self.paths ]
    # For each keyword parameter, remove candidates that don't match. 
    for key in kargs:
      candidates = [ p for p in candidates if self.runs[p][key]==kargs[key] ]
    # If there's anything other than exactly one match, something is wrong. 
    if len(candidates)<1:
      print 'ERROR: No matching path found for ', kargs
      exit()
    elif len(candidates)>1:
      print 'ERROR: Multiple matching paths found for ', kargs
      exit()
    else:
      return candidates[0]

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
    path = self.paths[0]
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
  # -------------------------------------------------
  # ---------------------------------------------------------------------------

  def plotFactor(self):
    PW = plotWindow(2, 1, colorbar='log')

    PW.setTitle( '\mathrm{Minimum \;\; Alfv\\\'en \;\; Frequency \;\; (mHz) }')

    # Cutoff when m=1, 16

    # Dayside and nightside look pretty similar. 

    path = self.getRunPath(azm='1', side='Night', fdrive='25')
    r, q = self.getArray(path + 'r.out'), self.getArray(path + 'q.out')
    X, Z = r*np.sin(q), r*np.cos(q)

    vA = self.getArray(path + 'vA.out') # km/s
    # Earth radius in km. 
    RE = 6378.388

    PW.setTitle( 'm = 1', pos=(0, 0) )
    PW.setContour( X, Z, 1000*vA/(RE*r), pos=(0, 0) )

    PW.setTitle( 'm = 16', pos=(1, 0) )
    PW.setContour( X, Z, 16*1000*vA/(RE*r), pos=(1, 0) )

    [ PW.setLine( X[i, :], Z[i, :] ) for i in (0, -1) ]
    [ PW.setLine( X[:, k], Z[:, k] ) for k in (0, -1) ]

    PW.setXlabel( 'X_{GSE} \;\; (R_E)')
    PW.setYlabel( 'Z_{GSE} \;\; (R_E)')

    # All of the data that the plot needs has been deposited in the plot window
    # object. This object no longer needs to keep track of anything. 
    self.refresh()
    # Either save the plot or show the plot. 
    if '-i' in self.flags:
      filename = self.outDir + 'cutoff.png'
      PW.save(filename)
      return
    else:
      return PW.show()














  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------------ Dst Plot
  # ---------------------------------------------------------------------------



  def plotDst(self):
    PW = plotWindow(2, 1, colorbar=False, xPad=2)

    filename = '/home/user1/mceachern/Desktop/tuna/symh/symh_20130601.txt'

    data = [ line for line in open(filename, 'r').readlines() if line.strip() ]

    t, SYMH = [], []

    for line in data:
      year, day, hour, minute, value = [ int(col) for col in line.split() ]
      t.append( 24*60*day + 60*hour + minute )
      SYMH.append(value) 

    t = np.array(t)
    t = t - t[0]
    SYMH = np.array(SYMH)

    # Plot SYMH. 
    PW.setLine(t, SYMH, pos=(0, 0), color='b')

    # Get a few Fourier components. 
    nModes = 20
    dt = t[1] - t[0]
    tRange = t[-1] - t[0]
    reconstruction = np.zeros(len(t), dtype=np.float)
    for m in range(nModes):
      harmonic = np.cos( m*np.pi*t / tRange )
      amplitude = np.sum( SYMH*harmonic*dt ) / np.sum( dt*harmonic**2 )
      reconstruction = reconstruction + amplitude*harmonic

    PW.setTitle( '\mathrm{SYMH \;\; Frequency \;\; Breakdown \;\; for \;\; June \;\; 2013 \;\; Storm}')

#    # Plot the fourier components. 
#    PW.setLine(t, reconstruction, pos=(0, 0), color='r')

#    PW.setTitle( '\mathrm{SYMH \;\; with \;\; ' + str(nModes) + ' \;\; Fourier \;\; Modes}', pos=(0, 0) )
    PW.setTitle( '\mathrm{SYMH \;\; Data \;\; from \;\; NASA \;\; CDAWeb}', pos=(0, 0) )
    PW.setYlabel( '\mathrm{SYMH \;\; (nT)}', pos=(0, 0) )
    PW.setXlabel( '\mathrm{Time \;\; (Minutes)}', pos=(0, 0) )

    # Now let's get a ton of Fourier components. 
    nModes = 4000

    amplitudes = np.zeros(nModes, dtype=np.float)

    for m in range(nModes):
      harmonic = np.cos( m*np.pi*t / tRange )
      amplitudes[m] = np.sum( SYMH*harmonic*dt ) / np.sum( dt*harmonic**2 )
    # Periods are in minutes. 
    periods = tRange*2./np.array( range(nModes) )
    # To get frequencies in seconds, multiply period by 60 then invert. For
    # mHz, put in a factor of 1000. 
    frequencies = 1000/(60*periods)
    PW.setLine( frequencies[1:], amplitudes[1:], pos=(1, 0), color='b' )
    PW.setXlog( pos=(1, 0) )
    PW.setYlog( pos=(1, 0) )
    # Set title and labels. 
    PW.setTitle( '\mathrm{SYMH \;\; Fourier \;\; Amplitudes}', pos=(1, 0) )
    PW.setYlabel( '\mathrm{Amplitude \;\; (nT)}', pos=(1, 0) )
    PW.setXlabel( '\mathrm{Frequency \;\; (mHz)}', pos=(1, 0) )
    # Set limits. They don't quite line up with the data, but that's OK. 
    fMin, fMax = 1e-2, 10
    PW.setXlimits( (fMin, fMax), pos=(1, 0) )
    f = np.linspace(fMin, fMax, 10000)
    # Let's not do a fit, per se -- let's eyeball the top of the distribution. 
    intercept = np.log(1.e-2)
    slope = -0.9
    # Scale to 20 mHz. 
    scale = 20. 

#    # Scale to a frequency of 16.7 mHz... 1000/(60s). 
#    scale = 1000/60.

    fit = np.exp(intercept) * np.power(f/scale, slope)

    label = ('$ ' + format( np.exp(intercept), '.2f' ) + '\; \mathrm{nT} ' + 
            '\cdot \left( \\frac{f}{' + format(scale, '.0f') +
            ' \; \mathrm{mHz} } \\right) ^{' + 
            format(slope, '.1f') + '}$')

    PW.setLine( f, fit, pos=(1, 0), color='r', label=label)

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
    return






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

  # ---------------------------------------------------------------------------
  # ----------------------------- Initialize Plot Window and Space Out Subplots
  # ---------------------------------------------------------------------------

  def __init__(self, nCols=1, nRows=1, colorbar=None, xPad=1, yPad=5):
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
    # terms of a grid of equally-sized tiles. 
    xSize = 10
    ySize = 10
    # Options are 'linear', 'log', and 'loglog'. Boolean True gives a linear
    # scale, while False and None give no color bar. For now, we just care if
    # there will be a color bar at all, so we can make room for it. 
    self.colorbar = colorbar
    # Set up the tiles. If we want a color bar, we need more columns. 
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
    # The plot window is now ready for us to add some data to it! 
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------ Set Title and Labels
  # ---------------------------------------------------------------------------

  # If we're not given a position, this is meant to be the super title. 
  def setTitle(self, text, pos=None):
    if pos is None:
      return plt.suptitle('$' + text + '$', fontsize=24)
      # Squish the subplots a little bit to make the title stand out. 
      plt.subplots_adjust(top=0.88)
    else:
      return self.cells[pos].setTitle(text)

  # If we're not given a position, apply this label to all edge cells. 
  def setXlabel(self, text, pos=None):
    if pos is None:
      return [ cell.setXlabel(text) for cell in self.cells[:, -1] ]
    else:
      return self.cells[pos].setXlabel(text)

  # If we're not given a position, apply this label to all edge cells. 
  def setYlabel(self, text, pos=None):
    if pos is None:
      return [ cell.setYlabel(text) for cell in self.cells[0, :] ]
    else:
      return self.cells[pos].setYlabel(text)

  # Set a label off to the right of a row. 
  def setRowLabel(self, text, row):
    pos = self.axes[0, row].get_position()
    return plt.figtext(0.3*pos.x0, pos.y0 + 0.5*pos.height, '$' + text + '$',
                       horizontalalignment='center', 
                       verticalalignment='center')

  # ---------------------------------------------------------------------------
  # ----------------------------------------------- Set Contour and Line Arrays
  # ---------------------------------------------------------------------------

  # If no position is given, all plots will draw the contour. 
  def setContour(self, X, Y, Z, pos=None):
    if pos is None:
      for column in self.cells:
        [ cell.setContour(X, Y, Z) for cell in column ]
    else:
      self.cells[pos].setContour(X, Y, Z)
    return

  # If no position is given, all plots will draw the line. 
  def setLine(self, X, Y, pos=None, color='k', label=''):
    if pos is None:
      for column in self.cells:
        [ cell.setLine(X, Y, color=color, label=label) for cell in column ]
    else:
      self.cells[pos].setLine(X, Y, color=color, label=label)
    return

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Set Plot Domain
  # ---------------------------------------------------------------------------

  # If no position is given, adjust limits on all plots. 
  def setXlimits(self, limits, pos=None):
    if pos is None:
      for column in self.cells:
        [ cell.setXlimits(limits) for cell in column ]
    else:
      self.cells[pos].setXlimits(limits)
    return

  # If no position is given, adjust limits on all plots. 
  def setYlimits(self, limits, pos=None):
    if pos is None:
      for column in self.cells:
        [ cell.setYlimits(limits) for cell in column ]
    else:
      self.cells[pos].setYlimits(limits)
    return

  # If no position is given, use a log scale for all plots. 
  def setXlog(self, pos=None):
    if pos is None:
      for column in self.cells:
        [ cell.setXlog() for cell in column ]
    else:
      self.cells[pos].setXlog()
    return

  # If no position is given, use a log scale for all plots. 
  def setYlog(self, pos=None):
    if pos is None:
      for column in self.cells:
        [ cell.setYlog() for cell in column ]
    else:
      self.cells[pos].setYlog()
    return

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Render Plot Window
  # ---------------------------------------------------------------------------

  # If no position is given, add a legend to all plots. 
  def setLegend(self, pos=None):
    if pos is None:
      for column in self.cells:
        [ cell.setLegend() for cell in column ]
    else:
      self.cells[pos].setLegend()
    return

  # Once all of the contours are loaded, we can figure out the color levels. 
  def render(self):
    # Get the maximum of all contour maxima. 
    vmax = nax( cell.getMax() for column in self.cells for cell in column )
    # The plot colors object constructor creates the color bar, then returns
    # a dictionary of keyword parameters. 
    PC = plotColors(vmax, self.colorAxis, self.colorbar)
    # We send those parameters to each cell for use in their contourf calls. 
    return [ cell.render(**PC) for column in self.cells for cell in column ]

  # ---------------------------------------------------------------------------
  # ------------------------------------------------ Show Plot Window on Screen
  # ---------------------------------------------------------------------------

  def show(self):
    # Before showing the window, assemble the contour plots, color bar, etc. 
    self.render()
    return plt.show()

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Save Plot Window as PNG
  # ---------------------------------------------------------------------------

  def save(self, filename):
    # Before saving the image, assemble the contour plots, color bar, etc. 
    self.render()
    # If the output directory doesn't exist, make it. This ensures that we
    # don't create an output directory unless we're actually making output. 
    savePath = os.path.dirname(filename)
    if not os.path.exists(savePath):
      os.makedirs(savePath)
    # Create the image. 
    plt.savefig(filename)
    print 'Saved plot to ' + filename.replace(os.environ['HOME'], '~')
    return

# #############################################################################
# ############################################################ Plot Cell Object
# #############################################################################

class plotCell:

  # A place to hold the coordinates, and the data, in case of a contour. 
  X, Y, Z = None, None, None
  # If any lines are to be drawn on the plot, they are stored here. 
  lines = ()
  # By default, do not use a legend.
  legend = False

  # ---------------------------------------------------------------------------
  # ------------------------------------- Initialize Plot Cell from Axis Object
  # ---------------------------------------------------------------------------

  def __init__(self, ax):
    self.ax = ax
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------ Set Title and Labels
  # ---------------------------------------------------------------------------

  def setTitle(self, text):
    return self.ax.set_title('$' + text + '$')

  def setXlabel(self, text):
    return self.ax.set_xlabel('$' + text + '$')

  def setYlabel(self, text):
    return self.ax.set_ylabel('$' + text + '$')

  # ---------------------------------------------------------------------------
  # ----------------------------------------------- Set Contour and Line Arrays
  # ---------------------------------------------------------------------------

  # Field values are scaled to the unit interval. 
  def setContour(self, X, Y, Z):
    self.X, self.Y, self.Z = X, Y, Z
    return

  def setLine(self, X, Y, color='k', label='$ \mathrm{label} $'):
    self.lines = self.lines + ( (X, Y, color, label), )
    return

  def setLegend(self):
    self.legend = True
    return

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Set Plot Domain
  # ---------------------------------------------------------------------------

  def setXlimits(self, limits):
    return self.ax.set_xlim(limits)

  def setYlimits(self, limits):
    return self.ax.set_ylim(limits)

  def setXlog(self):
    return self.ax.set_xscale('log')

  def setYlog(self):
    return self.ax.set_yscale('log')

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
      self.ax.plot( line[0], line[1], line[2], label=line[3] )
    # Draw the legend. 
    if self.legend:
      self.ax.legend(loc='best')
    # Stylistic stuff... none of this is actually important. 
    # Remove ticks without removing tick labels. 
#    self.ax.tick_params( width=0 )
    # Remove plot frames. 
    for pos in ('top', 'bottom', 'left', 'right'):
      self.ax.spines[pos].set_visible(False)
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

  # Number of color levels for contour plots. Should be even. 
  nColors = 8
  # Number of color bar ticks. Zero is shown if it's odd.  
  nTicks = 7

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Initialize Colors
  # ---------------------------------------------------------------------------

  def __init__(self, vmax=None, ax=None, colorbar=None):
    # Some plots don't have contours. 
    if not vmax or not ax or not colorbar:
      return dict.__init__(self, {})
    # Store the data scale so that it can be hard-wired into our normalization
    # functions. We don't want to pass vmax all over the place. 
    self.vmax = vmax
    self.colorbar = colorbar
    # Assemble the keyword parameters in a temporary dictionary. We'll then use
    # the dictionary constructor to build this object based on it. 
    temp = {}
    # Determine location of contour color levels and color bar ticks. 
    if self.colorbar=='log':

      self.nColors = 12
      self.nTicks = 5

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
    self.vmin = min(posTicks)/10
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
      return plt.get_cmap('seismic')


      N = 1000

      # At first, we just grab from the top half. 
      unitInterval = [ 0.5 + 0.5*i/(N - 1.) for i in range(N) ]
      cmap = plt.get_cmap('seismic')
      rgb = [ cmap( unitInterval[i] ) for i in range(N) ]

      newInterval = [ i/(N - 1.) for i in range(N) ]
      newInterval[0], newInterval[-1] = 0., 1.
      red = [ (newInterval[i], rgb[i][0], rgb[i][0]) for i in range(N) ]
      grn = [ (newInterval[i], rgb[i][1], rgb[i][1]) for i in range(N) ]
      blu = [ (newInterval[i], rgb[i][2], rgb[i][2]) for i in range(N) ]
      return LSC('myMap', {'red':red, 'green':grn, 'blue':blu}, 1000000)







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
    power = format(np.log10( abs(x) ), '.0f')
    return '$ ' + ( '-' if x<0 else '+' ) + ' 10^{' + power + '}$'

# #############################################################################
# ############################################################ Helper Functions
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

# Turns the number 3 into '003'. If given a float, truncate it. 
def znt(x, width=0):
  return str( int(x) ).zfill(width)

# =============================================================================
# ========================================================= Shell Input Parsing
# =============================================================================

# Check if the given path is a data directory. 
def isDataPath(path):
  return os.path.isdir(path) and 'params.in' in os.listdir(path)

# Format a path nicely. Always use absolute paths, and always end directory
# names with a slash. 
def getPath(p):
  if p.startswith('.'):
    return (os.getcwd() + '/' + p[1:] + '/').replace('//', '/')
  else:
    return (p + '/').replace('//', '/')

# Grab the arguments from the shell and partition them into the names of plots 
# to make, the paths where those plots should be created, and the flags that
# will affect how those plots are displayed. 
def getArgs():
  flags, paths, names = [], [], []
  # The first argument is the script call. 
  for arg in argv[1:]:
    # Grab flags. Anything that starts with a dash. 
    if arg.startswith('-'):
      flags.append( arg.lower() )
    # Grab directories. Any directory names that contain data. 
    elif isDataPath(arg):
      paths.append( getPath(arg) )
    elif len(arg)==1:
      names.append( arg.lower() )
    else:
      # Other arguments are probably paths that are not data directories. This
      # happens a lot if called with a wildcard. 
      pass
  return flags, paths, names

# =============================================================================
# ============================================================ Data File Access
# =============================================================================

def readArray(filename):
  # If a pickle is available, use that. 
  pklname = filename[:filename.rfind('.')] + '.pkl'
  if os.path.isfile( pklname ):
    print 'Reading ' + '/'.join(pklname.split('/')[-2:])
    with open(pklname, 'rb') as handle:
      return pickle.load(handle)
  # Make sure the file exists before trying to access it. 
  elif os.path.isfile(filename):
    print 'Reading ' + '/'.join(filename.split('/')[-2:])
  else:
    print 'WARNING: ' + '/'.join(filename.split('/')[-2:]) + ' not found. '
    return None
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
    print '\tWARNING: Expected ' + by(dims) + ' but found ' + by(actualDims)
    return arr[..., :actualDims[-1] ]

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


