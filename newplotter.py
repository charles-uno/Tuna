#!/usr/bin/env python

# Charles McEachern

# Fall 2015

# Note: This document wraps at column 80. 

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This plotter is specifically designed to make plots that look good on paper.
# The width is fixed to remove the need for rescaling, so that fonts show up
# true to size. Additionally, the same amount of space is reserved for the
# title, labels, and color bar in each plot, regardless of how many subplots
# there are (up to 4 columns). 

# #############################################################################
# ##################################################### Import Python Libraries
# #############################################################################

# The cPickle module is faster, but not always available. 
try:
  import cPickle as pickle
except ImportError:
  import pickle
from matplotlib import gridspec, rc
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap as LSC
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import basename
from sys import argv, stdout
from time import localtime as lt, time

from random import choice

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Create the Tuna Plotter object, a Tuna-specific wrapper around the Plot
  # Window object. If the -i flag was given from the terminal, tell it that we
  # want the output to be saved as in image instead of displayed. 
  TP = tunaPlotter('-i' in argv)


#  for kargs in loopover( fdrive=TP.getValues('fdrive'), mode=('BE', 'PT'), pick=False ):
#    TP.plotUlines(**kargs)

  TP.plotUlines(fdrive=0.022, mode='BE')

#  for kargs in loopover( fdrive=TP.getValues('fdrive'), mode=('B', 'E', 'pol', 'tor'), pick=True ):
#  for kargs in loopover( fdrive=TP.getValues('fdrive'), mode=('pol', 'tor'), pick=False ):
#    TP.plotucolor(**kargs)





#  for kargs in loopover( model=(1, 2), fdrive=TP.getValues('fdrive'), compare=('S', 'u'), pick=False ):
#      TP.plotJetc(**kargs)

#  for kargs in loopover( model=(1,2,3,4), mode=('', 'pol', 'tor'), pick=False ):
#    TP.plotLayers(driving='J', **kargs)

#  TP.plotGrid()

#  TP.plotBounce()

#  TP.plotSymh()

#  for model in (1, 2, 3, 4):
#    for fdrive in TP.getValues('fdrive'):
#      for alt in ('RE', 'RI', 'RX'):
#        for azm in TP.getValues('azm'):
#          TP.plotEdge(model=model, fdrive=fdrive, alt=alt, driving='J', azm=azm)

#  for kargs in loopover( model=(1, 2), fdrive=TP.getValues('fdrive'), pick=False ):
#      TP.plotJdotE(inertia=1, **kargs)

#  for model in (1, 2):
#    for fdrive in TP.getValues('fdrive'):
#      for azm in TP.getValues('azm'):
#        TP.plotInertia(model=model, fdrive=fdrive, azm=azm)

#  for path in TP.paths:
#    TP.plotJ(path)

#  TP.plotAlfvenSpeed()

#  for model in (1, 2, 3, 4):
#    TP.plotToroidal(model=model)

#  # Integrating the energy over the whole domain only tells us so much. Let's
#  # take some snapshots of each run, so that we actually know something about
#  # what they look like. 
#  for path in TP.paths:
#    for phase in (False, True):
#      TP.plotS(path, phase=phase)

#  for model in (1, 2, 3, 4):
#    for driving in ('B', 'J'):
#      TP.plotUPUT(model=model, driving=driving)

#  for model in (1, 2, 3, 4):
#    for driving in ('B', 'J'):
#      TP.plotUBUE(model=model, driving=driving)

#  for model in (1,):
#    TP.plotubar(model)

#  TP.plotSigma()



  return








  if 'energy' in argv:
    TP.plotA(filename='UP_UT.pdf')

  if 'snapshot' in argv or 'lazy' in argv:
    TP.plotB()

  if 'ground' in argv:
    TP.plotC()

  if 'rms' in argv:
    TP.plotD()

  if 'u' in argv or 'shells' in argv:
    TP.plotE()

  if 'sigma' in argv or 'conductivity' in argv:
    TP.plotF(filename='sigma.pdf')

  if 'power' in argv:
    TP.plotG()

  if 'drive' in argv or 'driving' in argv:
    TP.plotH()

  if 'compare' in argv:
    TP.plotI(filename='UB_UE.pdf')

  if 'poloidal' in argv:
    TP.plotJ(filename='SP_unwrapped.pdf')


  if 'RI' in argv:
    TP.plotK(filename='SP_RI.pdf')



  return


# Timestamp for labeling output. 
def now():
  return ( znt(lt().tm_year, 4) + znt(lt().tm_mon, 2) + znt(lt().tm_mday, 2) +
           '_' + znt(lt().tm_hour, 2) + znt(lt().tm_min, 2) +
           znt(lt().tm_sec, 2) )

# Turns the number 3 into '003'. If given a float, truncate it. 
def znt(x, width=0):
  return str( int(x) ).zfill(width)


# #############################################################################
# ######################################################### Tuna Plotter Object
# #############################################################################

# The Tuna Plotter is a wrapper around the Plot Window class which includes
# utilities specific to Tuna. (The Plot Window itself is fairly general.)
class tunaPlotter:

  # Keep track of the paths where our data is coming from. 
  paths = None
  # Physical constants, for crunching out the Poynting flux, etc. 
  mu0 = 1256.63706 # nH/m
  eps0 = 8.854e-9 # mF/m
  qe = 1.60218e-25 # MC
  me = 9.10938e-28 # g
  RE = 6.378388 # Mm
  RI = 6.478388 # Mm
  # Remember the last few arrays we've read in to avoid duplicating our work. 
  arrays = None
  # Some runs have been longer than 300 seconds. All runs are at least 300
  # seconds, and that turns out to be a pretty good window for looking at Pc4s. 
  tmax = 300

  # ===========================================================================
  # ======================================================== Initialize Plotter
  # ===========================================================================

  def __init__(self, save=False):
    # Check if we're supposed to be saving this image or showing it. 
    if save:
      self.savepath = '/home/user1/mceachern/Desktop/plots/' + now() + '/'
      os.mkdir(self.savepath)
    else:
      self.savepath = None
    # Check for any path(s) from the terminal. 
    self.setPaths(*argv)
    # If no paths were given, use the default. 
    if not self.paths:
      self.setPaths('/media/My Passport/RUNS/')
    return

  # ===========================================================================
  # ================================================ Finding and Filtering Data
  # ===========================================================================

  # Given one or more paths, find all data directories, as well as the
  # parameters used by each run. 
  def setPaths(self, *args):
    # Note that self.paths isn't a list of paths to run directories -- it's a
    # dictionary, keyed by those paths, which lists the parameters used by each
    # run. We can stoll loop over path in self.paths, but this makes it easier
    # to track down the run we're looking for based on its parameters. 
    self.paths = {}
    # Increment over everything we're passed. 
    for arg in args:
      # Ignore non-paths. 
      if os.path.isdir(arg):
        for root, dirs, files in os.walk(arg):
          if 'params.in' in files:
            # Always use absolute paths. End directory paths with slashes. 
            path = os.path.abspath(root) + '/'
            self.paths[path] = self.getParams(path)
    return

  # Grab all parameters from a parameter input file. 
  def getParams(self, path):
    # Parameters are returned as a dictionary. Make sure that every dictionary
    # contains both bdrive and jdrive (one of which will be overwritten). Also
    # note that the default value is to have inertial effects turned off. 
    params = {'bdrive':0, 'jdrive':0, 'inertia':-1}
    # The input file has one parameter per line, 'key = val'. 
    with open(path + 'params.in', 'r') as paramsfile:
      paramlines = paramsfile.readlines()
    for line in paramlines:
      key, val = [ x.strip() for x in line.split('=') ]
      params[key] = num(val)
    return params

  # Find all values for a given keyword parameter that are present in the runs
  # we're looking at. 
  def getValues(self, key):
    values = []
    # We want to be able to ask about keys that are not present. 
    for params in self.paths.values():
      values.append( params[key] if key in params else None )
    # Make sure we can safely assume the values are in ascending order. 
    return sorted( set(values) )

  # Given some keyword arguments, filter self.paths. 
  def getPath(self, **kargs):
    # Start with all the paths we have, then weed out any that don't match. 
    paths = [ p for p in self.paths ]
    for key, val in kargs.items():
      paths = [ p for p in paths if self.paths[p][key]==val ]
    # If there's anything other than exactly one match, something is wrong. 
    if len(paths)<1:
      print 'WARNING: No matching path found for ', kargs
      return None
    elif len(paths)>1:
      print 'ERROR: Multiple matching paths found for ', kargs
      exit()
    else:
      return paths[0]

  # Find the most recent run, including runs still in progress. 
  def newPath(self):
    root = '/export/scratch/users/mceachern/'
    runs = [ d for d in os.listdir(root) if os.path.isdir(root + d) ]
    # Sort by time of last modification, then return the last path. 
    return max( (os.path.getmtime(d), d) for d in runs )[1]

  # ===========================================================================
  # ==================================================== LaTeX Helper Functions
  # ===========================================================================

  def texText(self, x):
    return '\\operatorname{' + x.replace(' ', '\\;') + '}'

  def texName(self, x):
    # Dictionary of strings we might need. 
    names = {
             # Spell out what each model means. 
             1:self.texText('Active Day'),
             2:self.texText('Quiet Day'),
             3:self.texText('Active Night'),
             4:self.texText('Quiet Night'),
             # Names for fields, axes, etc. 
             'alt':self.texText('Altitude'), 
             'B':self.texText('Compression'), 
             'Bf':'B_\\phi', 
             'BfE':'B_\\phi' + self.texText(' at R_E'), 
             'BfI':'B_\\phi' + self.texText(' at R_I'), 
             'Bq':'B_\\theta', 
             'BqE':'B_\\theta' + self.texText(' at R_E'), 
             'BqI':'B_\\theta' + self.texText(' at R_I'), 
             'Bx':'B_x', 
             'By':'B_y', 
             'Bz':'B_z', 
             'C':'\\cos\\theta / \\cos\\theta_0', 
             'Ex':'E_x', 
             'Ey':'E_y', 
             'Ez':'E_z', 
             'f':self.texText('Frequency'), 
             'imag':self.texText('\\mathbb{I}m'), 
             'J':self.texText('Current'), 
             'Jz':'J_\\parallel', 
             'JzI':'J_\\parallel', 
             'L':self.texText('L'), 
             'L0':self.texText('L'), 
             'lat':self.texText('Latitude'), 
             'lat0':self.texText('Latitude'), 
             'logU':self.texText('Log ') + 'U', 
             'logsigma':self.texText('Log Conductivity'), 
             'logsymh':self.texText('Log Amplitude'), 
             'logalt':self.texText('Log Altitude'), 
             'real':self.texText('\\mathbb{R}e'), 
             'RE':self.texText('Earth\'s Surface'), 
             'RI':self.texText('the Top of the Atmosphere'), 
             'RX':self.texText('the Bottom of the Ionosphere'), 
             'sigma':self.texText('Conductivity'),
             'S':'S_P + S_T',
             'Spol':'S_P',
             'Stor':'S_T',
             'Sym-H':self.texText('Amplitude'),
             'symh':self.texText('Amplitude'),
             't':self.texText('Time'), 
             'u':'u_P + u_T', 
             'upol':'u_P', 
             'utor':'u_T', 
             'U':self.texText('Energy'), 
             'v':self.texText('Alfv\\acute{e}n Speed'), 
             'va':self.texText('Alfv\\acute{e}n Speed'), 
             'X':'X', 
             'Z':'Z'
            }
    return '?' if x not in names else names[x]

  # If there are only one or two columns, we can fit words in the title instead
  # of symbols. 
  def texTitle(self, x):
    # Dictionary of strings we might need. 
    titles = {
              # Spell out what each model means. 
              1:self.texText('Active Day'),
              2:self.texText('Quiet Day'),
              3:self.texText('Active Night'),
              4:self.texText('Quiet Night'),
              # Names for fields, axes, etc. 
              'alt':self.texText('Altitude'), 
              'B':self.texText('Compression'), 
              'f':self.texText('Frequency'), 
              'imag':self.texText('\\mathbb{I}m'), 
              'J':self.texText('Current'), 
              'lat':self.texText('Latitude'), 
              'lat0':self.texText('Latitude'), 
              'logU':self.texText('Log Energy'), 
              'logsigma':self.texText('Log Conductivity'), 
              'logsymh':self.texText('Log Amplitude'), 
              'logalt':self.texText('Log Altitude'), 
              'real':self.texText('\\mathbb{R}e'), 
              'RE':self.texText('Earth\'s Surface'), 
              'RI':self.texText('the Top of the Atmosphere'), 
              'RX':self.texText('the Bottom of the Ionosphere'), 
              'sigma':self.texText('Conductivity'),
              'Sym-H':self.texText('Amplitude'),
              'symh':self.texText('Amplitude'),
              't':self.texText('Time'), 
              'U':self.texText('Energy'), 
              'v':self.texText('Alfv\\acute{e}n Speed'), 
              'va':self.texText('Alfv\\acute{e}n Speed'), 
             }
    return '?' if x not in titles else titles[x]

  def texReIm(self, x):
    if x in ('Ex', 'By', 'Ez', 'Jz', 'BqE', 'BqI'):
      return ' \\mathbb{I}\\mathrm{m}\\;\\; '
    elif x in ('Bx', 'Ey', 'Bz', 'BfE', 'BfI'):
      return ' \\mathbb{R}\\mathrm{e}\\;\\; '
    else:
      return ''

  def texTime(self, x):
    # Time is a float, but we don't need a '.0' appended to everything. 
    strx = str( int(x) ) if int(x)==x else str(x)
    return self.texText(' at ' + strx + 's')

  def texUnit(self, x):
    # Dictionary of units we might care about. 
    units = {
             'alt':'km',
             'B':'nT',
             'Bf':'nT', 
             'Bq':'nT', 
             'Bx':'nT', 
             'By':'nT', 
             'Bz':'nT', 
             'C':None,
             'deg':'^\\circ',
             'E':'\\frac{mV}{m}',
             'Ex':'\\frac{mV}{m}',
             'Ey':'\\frac{mV}{m}',
             'Ez':'\\frac{mV}{m}',
             'f':'mHz',
             'JE':'\\frac{nW}{m^3}',
             'J':'\\frac{\\mu\\!A}{m^2}',
             'Jz':'\\frac{\\mu\\!A}{m^2}',
             'km':'km',
             'L':'R_E',
             'L0':'R_E',
             'lat':'^\\circ',
             'lat0':'^\\circ',
             'logU':'\\frac{GJ}{rad}',
             'logsigma':'\\frac{mS}{m}',
             'logsymh':'nT',
             'logalt':'km',
             'mHz':'mHz',
             'nT':'nT',
             's':'s',
             'sigma':'\\frac{mS}{m}',
             'S':'\\frac{mW}{m^2}',
             'Stor':'\\frac{mW}{m^2}',
             'Spol':'\\frac{mW}{m^2}',
             'symh':'nT',
             't':'s',
             'u':'\\frac{nJ}{m^3}',
             'upol':'\\frac{nJ}{m^3}',
             'utor':'\\frac{nJ}{m^3}',
             'U':'\\frac{GJ}{rad}',
             'v':'\\frac{Mm}{s}',
             'va':'\\frac{Mm}{s}',
             'X':'R_E', 
             'Z':'R_E'
            }
    if x not in units:
      return self.texText(' (?)')
    elif units[x] is None:
      return ''
    else:
      return self.texText(' (' + units[x] + ')')

  def texLabel(self, x, units=True):
    return '{\\displaystyle ' + self.texName(x) + '}' + self.texUnit(x)

  def texFreq(self, x):
    return self.texText(format(1e3*x, '.0f') + 'mHz ')

  # ===========================================================================
  # =============================================================== Data Access
  # ===========================================================================

  # This wrapper on readArray remembers the last few files we've read in, in
  # case we need them again. 
  def readArray(self, filename):
    if self.arrays is None or len(self.arrays)>20:
      self.arrays = {}
    if filename not in self.arrays:
      self.arrays[filename] = readArray(filename)
    return self.arrays[filename]

  # All arrays come in through here. It knows which fields are predominantly
  # imaginary (well, it knows how to look it up), and it knows how to combine
  # fields to get energy density, Poynting flux, etc. This is how we keep the
  # individual plot methods clean. 
  def getArray(self, path, name, makeReal=True):
    # Check if we're looking for something that corresponds to a data file...
    if name in ('BfE', 'BfI', 'BqE', 'BqI', 'Bx', 'By', 'Bz', 'Ex', 'Ey', 
                'Ez', 'JyDrive', 'Jz', 'q', 't'):
      # Usually we just want whichever complex component is larger. 
      if makeReal is True:
        phase = np.imag if '\\mathbb{I}' in self.texReIm(name) else np.real
        arr = phase( self.readArray(path + name + '.dat') )
      else:
        arr = self.readArray(path + name + '.dat')
      # Chop off anything after tmax time steps. 
      return arr[..., :self.tmax] if name=='t' or len(arr.shape)>2 else arr
    # A few quantities get printed out with scale factors. Un-scale them. 
    # Radius in Mm (from RE). 
    elif name=='r':
      return self.RE*self.readArray(path + 'r.dat')
    # Perpendicular electric constant, from units of eps0 to mF/m. 
    elif name=='epsPerp' or name=='epsp':
      return self.eps0*self.readArray(path + 'epsPerp.dat')
    # Conductivities. Printed in S/m but we want mS/m. 
    elif name in ('sigH', 'sigP', 'sig0'):
      return 1e-3*self.readArray(path + name + '.dat')
    # Altitude in km. Note that we read r in RE and we store RE in Mm. 
    elif name=='alt':
      return 1000*(self.getArray(path, 'r') - self.RE)
    # Normalied latitude. 
    elif name=='C':
      cosq = np.cos( self.getArray(path, 'q') )
      cosq0 = cosq[:, 0]
      return cosq/cosq0[:, None]
    # Driving electric field in mV/m. 
    elif name=='EyDrive':
      return self.getArray(path, 'JyDrive')/self.getArray(path, 'epsPerp')
    # Perpendicular currents, computed from electric fields and conductivity. 
    elif name=='Jx':
      Ex, Ey = self.getArray(path, 'Ex', makeReal=makeReal), self.getArray(path, 'Ey', makeReal=makeReal)
      sigP, sigH = self.getArray(path, 'sigP'), self.getArray(path, 'sigH')
      return sigP[:, :, None]*Ex - sigH[:, :, None]*Ey
    elif name=='Jy':
      Ex, Ey = self.getArray(path, 'Ex', makeReal=makeReal), self.getArray(path, 'Ey', makeReal=makeReal)
      sigP, sigH = self.getArray(path, 'sigP'), self.getArray(path, 'sigH')
      return sigH[:, :, None]*Ex + sigP[:, :, None]*Ey
    # Parallel current at the ionosphere. 
    elif name=='JzI':
      return self.getArray(path, 'Jz', makeReal=makeReal)[:, 0, :]
    # McIlwain parameter. 
    elif name=='L':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      return r / ( self.RE * np.sin(q)**2 )
    # McIlwain parameter as a 1D array. 
    elif name=='L0':
      return self.getArray(path, 'L')[:, 0]
    # Latitude in degrees, from colatitude in radians. 
    elif name=='lat':
      return 90 - self.getArray(path, 'q')*180/np.pi
    # Latitude only along the ionospheric boundary. 
    elif name=='lat0':
      return self.getArray(path, 'lat')[:, 0]
    # Differential volume for the grid, based on dipole coordinates and the
    # Jacobian determinant. 
    elif name=='dV':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      # Invariant colatitude. 
      cosq0 = np.sqrt( 1 - self.RI*np.sin(q)**2/r )
      # Dipole coordinates. 
      u1 = -self.RI/(r*self.RE) * np.sin(q)**2
      u3 = self.RI**2/(r*self.RE)**2 * np.cos(q)/cosq0
      # Dipole differentials. Rolling messes up the edges, so we zero them. 
      du1 = ( np.roll(u1, shift=1, axis=0) - np.roll(u1, shift=-1, axis=0) )/2
      du1[0, :], du1[-1, :] = 0, 0
      du3 = ( np.roll(u3, shift=1, axis=1) - np.roll(u3, shift=-1, axis=1) )/2
      du3[:, 0], du3[:, -1] = 0, 0
      # Jacobian determinant. 
      jac = (r*self.RE)**6/self.RI**3 * cosq0/( 1 + 3*np.cos(q)**2 )
      # The Jacobian may be negative. Make sure we return a positive volume. 
      return np.abs( du1*du3*jac )
    # Perpendicular grid spacing. This is cheating a little bit, since the
    # coordinates aren't orthogonal, but should give a decent estimate.  
    elif name=='dx0':
      return self.RI*self.d( self.getArray(path, 'q')[:, 0] )
    # Azimuthal effective grid spacing for taking derivatives. This assumes an
    # azimuthal mode number of 1, and scales linearly. 
    elif name=='dy0':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      return r[:, 0]*np.sin( q[:, 0] )
    # Differential field line length. 
    elif name=='dz':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      dr = self.d(r, axis=1)
      rdq = r*self.d(q, axis=1)
      return np.sqrt( dr**2 + rdq**2 )
    # Field line length right along the ionospheric boundary. 
    elif name=='dz0':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      dr = r[:, 1] - r[:, 0]
      rdq = 0.5*( r[:, 1] + r[:, 0] )*( q[:, 1] - q[:, 0] )
      return np.sqrt( dr**2 + rdq**2 )
    # Alfven bounce frequency. The axis will be set by the lines we draw. 
    elif name=='f':
      return None
    # Toroidal Poynting flux. 
    elif name=='Stor':
      return self.getArray(path, 'Ex', makeReal=makeReal)*np.conj( self.getArray(path, 'By', makeReal=makeReal) )/self.mu0
    # Poloidal Poynting flux. 
    elif name=='Spol':
      return -self.getArray(path, 'Ey', makeReal=makeReal)*np.conj( self.getArray(path, 'Bx', makeReal=makeReal) )/self.mu0
    # Toroidal Poynting flux. 
    elif name=='Sx':
      return self.getArray(path, 'Ey', makeReal=makeReal)*np.conj( self.getArray(path, 'Bz', makeReal=makeReal) )/self.mu0
    # Poloidal Poynting flux. 
    elif name=='Sy':
      return -self.getArray(path, 'Ex', makeReal=makeReal)*np.conj( self.getArray(path, 'Bz', makeReal=makeReal) )/self.mu0
    # Parallel Poynting flux. 
    elif name=='S':
      return self.getArray(path, 'Spol') + self.getArray(path, 'Stor')
    # Sometimes we don't actually need the array, such as when we're grabbing
    # the y axis for a line plot... so the axis will be set by line values. 
    elif name in ('logU', 'U'):
      return None
    # Poloidal magnetic field contribution to the energy density. 
    elif name=='uBx':
      return 0.5*self.getArray(path, 'Bx')**2 / self.mu0
    # Toroidal magnetic field contribution to the energy density. 
    elif name=='uBy':
      return 0.5*self.getArray(path, 'By')**2 / self.mu0
    # Magnetic field contribution to the energy density. 
    elif name=='uB':
      return self.getArray(path, 'uBx') + self.getArray(path, 'uBy')
    # Toroidal electric field contribution to the energy density. 
    elif name=='uEx':
      E, epsPerp = self.getArray(path, 'Ex'), self.getArray(path, 'epsPerp')
      return 0.5*epsPerp[:, :, None]*E[:, :, :]**2
    # Poloidal electric field contribution to the energy density. 
    elif name=='uEy':
      E, epsPerp = self.getArray(path, 'Ey'), self.getArray(path, 'epsPerp')
      return 0.5*epsPerp[:, :, None]*E[:, :, :]**2
    # Magnetic field contribution to the energy density. 
    elif name=='uE':
      return self.getArray(path, 'uEx') + self.getArray(path, 'uEy')
    # Poloidal energy density. 
    elif name=='upol':
      return self.getArray(path, 'uEy') + self.getArray(path, 'uBx')
    # Toroidal energy density. 
    elif name=='utor':
      return self.getArray(path, 'uEx') + self.getArray(path, 'uBy')
    # Total energy density. 
    elif name=='u':
      return self.getArray(path, 'upol') + self.getArray(path, 'utor')
    # Integrated magnetic energy. 
    elif name=='UB':
      ux, uy = self.getArray(path, 'uBx'), self.getArray(path, 'uBy')
      dV = self.getArray(path, 'dV')
      return np.sum( np.sum( (ux + uy)*dV[:, :, None], 1), 0)
    # Integrated electric energy. 
    elif name=='UE':
      ux, uy = self.getArray(path, 'uEx'), self.getArray(path, 'uEy')
      dV = self.getArray(path, 'dV')
      return np.sum( np.sum( (ux + uy)*dV[:, :, None], 1), 0)
    # Integrated poloidal energy. 
    elif name=='Upol':
      ux, uy = self.getArray(path, 'uBx'), self.getArray(path, 'uEy')
      dV = self.getArray(path, 'dV')
      return np.sum( np.sum( (ux + uy)*dV[:, :, None], 1), 0)
    # Integrated toroidal energy. 
    elif name=='Utor':
      ux, uy = self.getArray(path, 'uEx'), self.getArray(path, 'uBy')
      dV = self.getArray(path, 'dV')
      return np.sum( np.sum( (ux + uy)*dV[:, :, None], 1), 0)
    # Alfven speed. 
    elif name=='vA' or name=='va':
      return 1/np.sqrt( self.getArray(path, 'epsp')*self.mu0 )
    # GSE X in RE. 
    elif name=='X':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      return r*np.sin(q)/self.RE
    # GSE Z in RE. 
    elif name=='Z':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      return r*np.cos(q)/self.RE
    # Keep an eye out for typos. 
    else:
      print 'ERROR: Not sure how to get ' + name
      exit()

  # ===========================================================================
  # ============================================================ Data Finagling
  # ===========================================================================

  # Helper for when we need to take a derivative. This gets the difference
  # between adjacent values (which we then have to scale by dx, etc). 
  def d(self, arr, axis=0):
    darr = ( np.roll(arr, shift=-1, axis=axis) - 
             np.roll(arr, shift=1, axis=axis) )/2
    darr[0], darr[-1] = darr[1], darr[-2]
    return darr

  # ===========================================================================
  # ====================================================== Coordinate Shorthand
  # ===========================================================================

  # This function returns a dictionary of keyword arguments meant to be plugged
  # straight into plotWindow.setParams(). 
  def getCoords(self, path, xaxis='X', yaxis='Z', lim=None):
    # The "unwrap" flag overwrites the dipole default coordinates. 
    if '-u' in argv and xaxis=='X' and yaxis=='Z':
      xaxis, yaxis = 'C', 'L'
    coords = { 'x':self.getArray(path, xaxis), 'xlabel':self.texLabel(xaxis),
             'y':self.getArray(path, yaxis), 'ylabel':self.texLabel(yaxis) }
    # Latitude vs altitude isn't much good for plotting the whole dipole. Zoom
    # in on the ionosphere. 
    if xaxis=='lat' and yaxis=='alt':
      # Set the window range based on the latitudes seen at the northern
      # hemisphere boundary, and the altitudes we resolve at those latitudes.
      lat, alt = coords['x'], coords['y']
      xmin, xmax = np.min( lat[:, 0] ), np.max( lat[:, 0] )
      ymin = np.min(alt)
      # Allow the altitude maximum to be overwritten to zoom in. 
      ymax = np.max( np.where(lat>xmin, alt, 0) ) if lim is None else lim
      coords['xlims'], coords['ylims'] = (xmin, xmax), (ymin, ymax)
    # The first time output is at 1s, but we want to start the axis at zero. 
    if xaxis=='t':
      coords['x'] = coords['x'][:self.tmax]
      coords['xlims'] = (0, lim)
    # Dipole plots need outlines drawn on them. 
    if xaxis=='X' and yaxis=='Z':
      coords['outline'] = True
    # If we're looking at the log of the energy, we need a bottom limit. 
    if yaxis=='logU':
#      coords['ylims'] = (1 if lim is None else lim, None)
      # Let's force the plots to all have the same y axis limits. Note that
      # compressional driving imparts a lot more energy than current. 
      if self.paths[path]['bdrive']==0:
        coords['ylims'] = (2, 6)
      else:
        coords['ylims'] = (4, 8)
    # If we're looking at electromagnetic energy on the y axis, we want a log
    # scale, and we also need a minimum. 
    if yaxis=='U':
      coords['ylog'] = True
      coords['ylims'] = (10 if lim is None else lim, None)
    # For line plots, the y axis is specified by the data. 
    if coords['y'] is None:
      del coords['y']
    return coords

  # ===========================================================================
  # ============================== Line Plot of Energy vs Time, Fixed Frequency
  # ===========================================================================

  def plotUlines(self, fdrive, mode):
    azms = self.getValues('azm')
    models = self.getValues('model')
    PW = plotWindow(nrows=len(azms), ncols=len(models), colorbar=False)
    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
    colLabels = [ self.texName(model) for model in models ]

    if mode=='BE':
      modetitle = 'Electric (Blue) and Magnetic (Red)'
    else:
      modetitle = 'Poloidal (Blue) and Toroidal (Red)'

    title = self.texText(modetitle + ' Energy: ' + self.texFreq(fdrive) + 'Current')
    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)
    for row, azm in enumerate(azms):
      for col, model in enumerate(models):
        path = self.getPath(azm=azm, model=model, fdrive=fdrive, bdrive=0, inertia=-1)
        if path is None:
          continue

        PW[row, col].setParams( **self.getCoords(path, 't', 'logU') )

        if mode=='BE':
          U1, U2 = self.getArray(path, 'UE'), self.getArray(path, 'UB')
        else:
          U1, U2 = self.getArray(path, 'Upol'), self.getArray(path, 'Utor')

        t = self.getArray(path, 't')

        PW[row, col].setLine(t, np.log10(U1), 'b')
        PW[row, col].setLine(t, np.log10(U2), 'r')

        if mode=='BE':
          mean1, mean2 = np.mean(U1), np.mean(U2)
          PW[row, col].setLine(t, np.log10(mean1)*np.ones(t.shape), 'b:')
          PW[row, col].setLine(t, np.log10(mean2)*np.ones(t.shape), 'r:')

    if self.savepath is not None:
      return PW.render(self.savepath + 'U_' + mode + '_' + znt(1e3*fdrive, 3) + 'mHz.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ============= Contour Plot of Shell Energy Density vs Time, Fixed Frequency
  # ===========================================================================

  def plotucolor(self, fdrive, mode):
    azms = self.getValues('azm')
    models = self.getValues('model')
    PW = plotWindow(nrows=len(azms), ncols=len(models), colorbar='log', zmax=0.1)
    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
    colLabels = [ self.texName(model) for model in models ]

    modetitle = {'':'', 'B':'Magnetic ', 'E':'Electric ', 'pol':'Poloidal ', 'tor':'Toroidal '}[mode]

    title = self.texText(modetitle + 'Energy Density by L-Shell: ' + self.texFreq(fdrive) + 'Current')
    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)
    for row, azm in enumerate(azms):
      for col, model in enumerate(models):
        path = self.getPath(azm=azm, model=model, fdrive=fdrive, bdrive=0, inertia=-1)
        if path is None:
          continue
        PW[row, col].setParams( **self.getCoords(path, 't', 'L0') )
        u, dV = self.getArray(path, 'u' + mode), self.getArray(path, 'dV')

        dU = u*dV[:, :, None]
        UofL = np.sum(dU, axis=1)

        VofL = np.sum(dV, axis=1)
        # Careful... dV is 0 at the edges. 
        VofL[0], VofL[-1] = VofL[1], VofL[-2]

        uofL = UofL/VofL[:, None]

        PW[row, col].setContour(uofL)

    if self.savepath is not None:
      return PW.render(self.savepath + 'ucolor_' + mode + '_' + znt(1e3*fdrive, 3) + 'mHz.pdf')
    else:
      return PW.render()














  # ===========================================================================
  # ============================================== Parallel Current, etc, at RI
  # ===========================================================================

  def plotJetc(self, model, fdrive, compare):

    azms = self.getValues('azm')

    PW = plotWindow(ncols=4, nrows=len(azms), colorbar='sym', zmax=1)

    colLabels = [ self.texName('real') + self.texLabel('Jz'), 
                  self.texName('imag') + self.texLabel('Jz'), 
                  self.texLabel(compare + 'pol'), 
                  self.texLabel(compare + 'tor')             ]

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

    comparename = 'Energy Density' if compare=='u' else 'Poynting Flux'

    title = self.texText('Field-Aligned Current and ' + comparename + ' at R_I: ' + self.texName(model) + ', ' + self.texFreq(fdrive) + 'Current')

    name = 'J' + compare + '_' + str(model) + '_' + znt(1e3*fdrive, 3) + 'mHz'

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title)

    for row, azm in enumerate(azms):
      path = self.getPath(azm=azm, model=model, fdrive=fdrive, inertia=1)
      [ PW[row, col].setParams( **self.getCoords(path, 't', 'lat0') ) for col in range(4) ]

      Jz = self.getArray(path, 'Jz', makeReal=False)[:, 0, :]
      JzRe, JzIm = np.real(Jz), np.imag(Jz)

      pol = self.getArray(path, compare + 'pol')[:, 0, :]
      tor = self.getArray(path, compare + 'tor')[:, 0, :]

      PW[row, 0].setContour(JzRe)
      PW[row, 1].setContour(JzIm)
      PW[row, 2].setContour(pol)
      PW[row, 3].setContour(tor)

    if self.savepath is not None:
      return PW.render(self.savepath + name + '.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ================================================================== The Grid
  # ===========================================================================

  def plotGrid(self):
    PW = plotWindow(ncols=1, nrows=1, colorbar=False)
    PW.setParams( title=self.texText('Nonorthogonal Dipole Grid') )
    path = self.paths.keys()[0]
    X, Z = self.getArray(path, 'X'), self.getArray(path, 'Z')
    PW.setParams(  **self.getCoords(path) )
    nx, nz = X.shape
    stride = 5
    [ PW.setLine( X[i, :], Z[i, :], 'k' ) for i in range(0, nx, stride) ]
    [ PW.setLine( X[:, k], Z[:, k], 'k' ) for k in range(0, nz, stride) ]
    if self.savepath is not None:
      return PW.render(self.savepath + 'grid.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ================================================= Alfven Bounce Frequencies
  # ===========================================================================

  def plotBounce(self):
    PW = plotWindow(nrows=2, ncols=2, colorbar=None)

    # Set the labels and title. 
    colLabels = [ self.texText('Active'), self.texText('Quiet') ]
    rowLabels = [ self.texText('Day'), self.texText('Night') ]

#    title = self.texText('Alfv\\acute{e}n Bounce Frequencies')
    title = self.texText('ALFVEN BOUNCE STUFF')

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title)

    for i in range(4):

      path = self.getPath(jdrive=0, azm=1, model=i+1, fdrive=0.010)

      PW[i].setParams( **self.getCoords(path, 'L0', 'f') )

      va = self.getArray(path, 'va')

      dz = self.getArray(path, 'dz')

      dt = dz/va

      tbounce = 2*np.sum(dt, axis=1)
      fbounce = 1e3/tbounce

      L0 = self.getArray(path, 'L0')

#      PW[i].setLine(L0, fbounce)

#      PW[i].setParams( ylims=(0, 60) )

      dL = self.d(L0)

      fprime = self.d(fbounce)/dL
#      PW[i].setLine(L0[1:-1], fprime[1:-1])

      omega = 2*np.pi*fbounce
      omegaprime = self.d(omega) / self.d(L0)
      lam = 1/(2*np.pi*L0)
      tau = np.abs( self.d(lam) / self.d(omegaprime) )
      PW[i].setLine(L0[2:-2], 10000*tau[2:-2])

#      Lmin, Lmax = np.min(L0), np.max(L0)
#      Lrange = Lmax - Lmin
#      def harmonic(n):
#        return np.exp( 1j*np.pi*n*(L0 - Lmin)/(Lmax - Lmin) ) / np.sqrt(2*Lrange)
#      def harmonicprime(n):
#        return  ( 1j*np.pi*n/(Lmax - Lmin) )*np.exp( 1j*np.pi*n*(L0 - Lmin)/(Lmax - Lmin) ) / np.sqrt(2*Lrange)
#      nmodes = 1000
#      modes = np.linspace(-30, 30, nmodes)
#      dmode = modes[1] - modes[0]
#      weights = [ np.sum(harmonic(-modes[n])*fbounce*dL) for n in range(nmodes) ]
#      fseries = np.sum( weights[n]*harmonic(modes[n])*dmode for n in range(nmodes) )
#      PW[i].setLine( L0, np.real(fseries) )
#      fprimeseries = np.sum( weights[n]*harmonicprime(modes[n])*dmode for n in range(nmodes) )
#      PW[i].setLine( L0, np.real(fprimeseries) )

#      pc4min = 7*np.ones(L0.shape)
#      pc4max = 25*np.ones(L0.shape)
#      PW[i].setLine(L0, pc4min, 'r:')
#      PW[i].setLine(L0, pc4max, 'r:')

    if self.savepath is not None:
#      return PW.render(self.savepath + 'fa.pdf')
      return PW.render(self.savepath + 'ta.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ====================== Comparison of Runs With and Without Inertial Effects
  # ===========================================================================

  def plotInertia(self, model, fdrive, azm):
    tmax = 300
    fields = ('Bx', 'By', 'Bz')

    inertias = (-1, 1)

    PW = plotWindow(nrows=len(fields), ncols=len(inertias), colorbar='sym', zmax=100)

    rowLabels = [ self.texName(field) for field in fields ]

    colLabels = [ self.texText('No Inertial Effects'), 
                  self.texText('Inertial Effects') ]

    title = self.texText('Magnetic Fields from 300s of ' + self.texFreq(fdrive) + 'Current: ' + self.texName(model) + ', ') + 'm = ' + str(azm)

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title)

    for col, inertia in enumerate(inertias):
      for row, field in enumerate(fields):
        path = self.getPath(azm=azm, model=model, fdrive=fdrive, inertia=inertia, bdrive=0)
        f = self.getArray(path, field)[:, :, tmax - 1]
        PW[row, col].setParams( **self.getCoords(path) )
        PW[row, col].setContour(f)

    name = 'B_' + str(model) + '_' + znt(azm, 3) + '_' + znt(1e3*fdrive, 3) + 'mHz'
    if self.savepath is not None:
      return PW.render(self.savepath + name + '.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ====================================== What Lines Up with Parallel Current?
  # ===========================================================================

  def plotJ(self, path):
    params = self.getParams(path)

    if params['inertia']<0:
      print 'SKIPPING ' + path
      return

    modelname = str(params['model'])
    modeltitle = self.texName(params['model'])

    azmname = znt(params['azm'], 3)
    azmtitle = 'm = ' + str(params['azm'])

    drivename = znt(1e3*params['fdrive'], 3) + 'mHz'
    drivetitle = self.texFreq(params['fdrive']) + self.texText('Current')

    name = 'J_' + modelname + '_' + azmname + '_' + drivename
    title = self.texText('Snapshots at 300s: ' + drivetitle + ', ' + modeltitle + ', ') + azmtitle

    fields = ('Stor', 'utor', 'Jz')

    PW = plotWindow(ncols=1, nrows=len(fields), colorbar='sym', zmax=1)

    PW.setParams(rowlabels=fields, title=title)

    for row, field in enumerate(fields):

      PW[row].setParams( **self.getCoords(path) )

      PW[row].setContour( self.getArray(path, field) )

    # It's easier to flip through a bunch of PNGs than it is to flip through a bunch of PDFs. These images aren't going in the thesis. 
    if self.savepath is not None:
      return PW.render(self.savepath + name + '.png')
    else:
      return PW.render()

  # ===========================================================================
  # ============================== Power Density from J dot E and Poynting Flux
  # ===========================================================================

  def plotJdotE(self, model, fdrive, inertia=1):

    azms = self.getValues('azm')

    PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=10)

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

    colLabels = ( '{\\displaystyle \\partial_\\parallel S_P}', 
                  '{\\displaystyle \\partial_\\parallel S_T}', 
                  '{\\displaystyle J_\\parallel E_\parallel}',
                  '{\\displaystyle \\underline{J}_\\bot \\cdot \\underline{E}_\\bot}'
                 )

    title = self.texText('Ionospheric Power Density: ') + self.texName(model) + self.texFreq(fdrive) + self.texText(' Compression') + self.texUnit('JE')

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title)

    for row, azm in enumerate(azms):
      path = self.getPath(azm=azm, model=model, fdrive=fdrive, inertia=inertia)

      Spol = self.getArray(path, 'Spol')
      Stor = self.getArray(path, 'Stor')

      JEx = self.getArray(path, 'Jx')*self.getArray(path, 'Ex')
      JEy = self.getArray(path, 'Jy')*self.getArray(path, 'Ey')
      JEz = self.getArray(path, 'Jz')*self.getArray(path, 'Ez')

      dx = self.getArray(path, 'dx0')[:, None]
      dy = azm*self.getArray(path, 'dy0')[:, None]
      dz = self.getArray(path, 'dz0')[:, None]

      dxStor = self.d( Stor[:, 0, :] )/dx
      dxSpol = self.d( Spol[:, 0, :] )/dx

      dyStor = Stor[:, 0, :]/dy
      dySpol = Spol[:, 0, :]/dy

      dzStor = ( Stor[:, 1, :] - Stor[:, 0, :] )/dz
      dzSpol = ( Spol[:, 1, :] - Spol[:, 0, :] )/dz

      [ PW[row, col].setParams( **self.getCoords(path, 't', 'lat0') ) for col in range(4) ]

      PW[row, 0].setContour(dxSpol + dySpol + dzSpol)
      PW[row, 1].setContour(dxStor + dyStor + dzStor)
      PW[row, 2].setContour(JEz[:, 0, :])
      PW[row, 3].setContour(JEx[:, 0, :] + JEy[:, 0, :])

    name = 'JE_' + str(model) + '_' + znt(1e3*fdrive, 3) + 'mHz'
    if self.savepath is not None:
      return PW.render(self.savepath + name + '.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # =================================================== Contours of Edge Fields
  # ===========================================================================

  # We can look at the fields themselves or at the phase of those fields. 
  def plotEdge(self, model, azm, fdrive, driving='J', alt='RE'):

    # 500 seconds is too long compared to the decays we're looking at. 
    tmax = 300

    PW = plotWindow(nrows=2, ncols=2, colorbar='sym', zmax=100)

    if alt=='RE':
      fields = ('BqE', 'BfE')
    elif alt=='RI':
      fields = ('BqI', 'BfI')
    else:
      fields = ('Bx', 'By')

    colLabels = [ self.texName('real'), self.texName('imag') ]

    rowLabels = [ self.texName( f[:2] ) for f in fields ]

    title = self.texText( 'Fields at ' + self.texName(alt) + ': ' + self.texName(model) + ', ' + self.texFreq(fdrive) + self.texName(driving) + ', ' ) + 'm = ' + str(azm)

    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)

    drive = { ('bdrive' if driving=='J' else 'jdrive') : 0 }
    path = self.getPath(azm=azm, model=model, fdrive=fdrive, **drive)

    for row, field in enumerate(fields):
      for col, phase in enumerate( (np.real, np.imag) ):

        B = phase( self.getArray(path, field, makeReal=False)[:, 0, :] )

        PW[row, col].setParams( **self.getCoords(path, 't', 'lat0', lim=tmax) )
        PW[row, col].setContour(B)

    if self.savepath is not None:
      name = alt + '_' + str(model) + '_' + znt(azm, 3) + '_' + znt(fdrive*1e3) + 'mHz_' + driving
      return PW.render(self.savepath + name + '.png')
    else:
      return PW.render()

  # ===========================================================================
  # =================================================== Poynting Flux Snapshots
  # ===========================================================================

  def plotS(self, path, phase=False):

    params = self.getParams(path)

    if params['bdrive']==0:
      drivename = znt(1000*params['fdrive'], 2) + 'mHz_j'
      drivetitle = self.texText( format(1000*params['fdrive'], '.0f') + 'mHz Current' )
    else:
      drivename = znt(1000*params['fdrive'], 2) + 'mHz_b'
      drivetitle = self.texText( format(1000*params['fdrive'], '.0f') + 'mHz Compression' )

    modelname = 'model' + str( params['model'] )
    modeltitle = self.texName( params['model'] )

    azmname = 'm' + znt(params['azm'], 2)
    azmtitle = 'm = ' + str( params['azm'] )

    name = modelname + '_' + drivename + ('_phase' if phase is True else '') + '_' + azmname

    ph = 'Phase ' if phase is True else ''

    title = self.texText('Poynting Flux '  + ph + 'Snapshots: ' + drivetitle + ', ' + modeltitle + ', ') + azmtitle

    steps = (239, 249, 259, 269, 279, 289, 299)
    fields = ('Spol', 'Stor', 'Sx', 'Sy')

    colorbar = 'phase' if phase is True else 'sym'
    ncolors = 9 if phase is True else None
    PW = plotWindow(nrows=len(steps), ncols=len(fields), colorbar=colorbar, ncolors=ncolors)

    t = self.getArray(path, 't')

    colLabels = ( self.texText('Poloidal'), 
                  self.texText('Toroidal'), 
                  self.texText('Perpendicular'), 
                  self.texText('Azimuthal') )

    rowLabels = [ self.texText( format(t[s], '.0f') + ' s') for s in steps ]

    zmax = 1 if self.paths[path]['bdrive']==0 else 100

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title, zmax=zmax)

    for col, field in enumerate(fields):
      for row, step in enumerate(steps):

        f = self.getArray(path, field)

        PW[row, col].setParams(  **self.getCoords(path) )

        if phase is True:
          fI, fR = np.imag(f), np.real(f)
          PW[row, col].setContour( f[:, :, step] )
          fphase = np.arctan( np.abs(fI)/(np.abs(fR) + 1e-10) )
          PW[row, col].setContour( fphase[:, :, step] )
        else:
          PW[row, col].setContour( f[:, :, step] )

    # It's easier to flip through a bunch of PNGs than it is to flip through a bunch of PDFs. These images aren't going in the thesis. 
    if self.savepath is not None:
      return PW.render(self.savepath + name + '.png')
    else:
      return PW.render()

  # ===========================================================================
  # ========== Line Plot of Poloidal and Toroidal Energy vs Time, Fixed Profile
  # ===========================================================================

  def plotUPUT(self, model=1, driving='J'):
    # 500 seconds is too long compared to the decays we're looking at. 
    tmax = 300

    azms = self.getValues('azm')

    fdrives = self.getValues('fdrive')

    PW = plotWindow(nrows=len(azms), ncols=len(fdrives), colorbar=False)

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

    driveLabel = 'Current' if driving=='J' else 'Compression'

    colLabels = [ self.texText(format(1e3*f, '.0f') + 'mHz ' + driveLabel) for f in fdrives ]

    title = self.texText( 'Poloidal (Blue) and Toroidal (Red) Energy: ' +
                          self.texName(model) )
    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)

    for row, azm in enumerate(azms):
      for col, fdrive in enumerate(fdrives):

        if driving=='J':
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, bdrive=0)
        else:
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, jdrive=0)

        Utor = np.log10( self.getArray(path, 'Utor') )[:tmax]
        Upol = np.log10( self.getArray(path, 'Upol') )[:tmax]

        t = self.getArray(path, 't')[:tmax]
        coords = self.getCoords(path, 't', 'logU')
        coords['X'] = t

        PW[row, col].setParams( **coords )
        PW[row, col].setLine( t, Utor, 'r')
        PW[row, col].setLine(t, Upol, 'b')

    if self.savepath is not None:
      return PW.render(self.savepath + 'UP_UT_' + driving + '_' + str(model) + '.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ========================= Line Plot of Electric and Magnetic Energy vs Time
  # ===========================================================================

  def plotUBUE(self, model=1, driving='J'):
    # 500 seconds is too long compared to the decays we're looking at. 
    tmax = 300

    azms = self.getValues('azm')

    fdrives = self.getValues('fdrive')

    PW = plotWindow(nrows=len(azms), ncols=len(fdrives), colorbar=False)

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

    driveLabel = 'Current' if driving=='J' else 'Compression'

    colLabels = [ self.texText(format(1e3*f, '.0f') + 'mHz ' + driveLabel) for f in fdrives ]

    title = self.texText( 'Electric (Blue) and Magnetic (Red) Energy: ' +
                          self.texName(model) )
    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)

    for row, azm in enumerate(azms):
      for col, fdrive in enumerate(fdrives):

        if driving=='J':
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, bdrive=0)
        else:
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, jdrive=0)

        UB = np.log10( self.getArray(path, 'UB') )[:tmax]
        UE = np.log10( self.getArray(path, 'UE') )[:tmax]

        t = self.getArray(path, 't')[:tmax]
        coords = self.getCoords(path, 't', 'logU')
        coords['X'] = t

        PW[row, col].setParams( **coords )
        PW[row, col].setLine( t, UB, 'r')
        PW[row, col].setLine(t, UE, 'b')

        meanUB = np.mean(UB)
        meanUE = np.mean(UE)

        PW[row, col].setLine(t, meanUB*np.ones(t.shape), 'r:')
        PW[row, col].setLine(t, meanUE*np.ones(t.shape), 'b:')

    if self.savepath is not None:
      return PW.render(self.savepath + 'UB_UE_' + driving + '_' + str(model) + '.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ========================================= Ionospheric Conductivity Profiles
  # ===========================================================================

  def plotSigma(self):
    # Create the window. 
    PW = plotWindow(nrows=2, ncols=2, colorbar=None)
    # For this plot, we don't actually need the 2D arrays that Tuna spat out.
    # We can just read in the profiles directly. 
    for i in range(4):
      datname = './models/ionpar' + str(i+1) + '.dat'
      with open(datname, 'r') as fileobj:
        lines = fileobj.readlines()
      ionos = np.array( [ [ float(x) for x in l.split() ] for l in lines ] )
      # Chop off the altitudes at 100km and 10000km. 
      bottom = np.argmin(ionos[:, 0]<100)
      if np.max( ionos[:, 0] )>1e4:
        top = np.argmax(ionos[:, 0]>1e4)
      else:
        top = len( ionos[:, 0] )
      # Log altitude, rather than altitude on a log scale. 
      logalt = np.log10( ionos[bottom:top, 0] )
      # We have to worry about zero and negative values for the perpendicular
      # conductivity. Clip at the minimum positive value. 
      sigP = np.abs( 4*np.pi*self.eps0*ionos[bottom:top, 2] )
      sigH = np.abs( 4*np.pi*self.eps0*ionos[bottom:top, 3] )
      minP = np.min( sigP[ np.nonzero(sigP) ] )
      minH = np.min( sigP[ np.nonzero(sigH) ] )
      # Plot log conductivity instead of conductivity on a log scale. 
      logsig0 = np.log10( 1e6/( self.mu0*ionos[bottom:top, 4] ) )
      logsigH = np.log10( np.clip(sigH, minH, np.inf) )
      logsigP = np.log10( np.clip(sigP, minP, np.inf) )
      # Add the lines to the plot. 
      PW[i].setLine(logsigP, logalt, 'r')
      PW[i].setLine(logsigH, logalt, 'b')
      PW[i].setLine(logsig0, logalt, 'g')
    # Set the labels and title. 
    colLabels = [ self.texText('Active'), self.texText('Quiet') ]
    rowLabels = [ self.texText('Day'), self.texText('Night') ]
    title = self.texText('Pedersen (Blue), Hall (Red), and Parallel (Green) ' +
                         'Conductivities')
    PW.setParams(collabels=colLabels, rowlabels=rowLabels, nxticks=5, 
                 xlabel=self.texLabel('logsigma'), 
                 ylabel=self.texLabel('logalt'), title=title)

    if self.savepath is not None:
      return PW.render(self.savepath + 'sigma.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ====================================================== Alfven Speed Profile
  # ===========================================================================

  def plotAlfvenSpeed(self):
    # Create the window. 
    PW = plotWindow(nrows=2, ncols=2, colorbar='log')

    # We don't actually care about fdrive, azm, or driving style. Just model. 

    azm = self.getValues('azm')[0]
    fdrive = self.getValues('fdrive')[0]

    for i in range(4):

      path = self.getPath(azm=azm, fdrive=fdrive, bdrive=0, model=i+1)

      PW[i].setParams( **self.getCoords(path) )

      PW[i].setContour( self.getArray(path, 'va') )

    # Set the labels and title. 
    colLabels = [ self.texText('Active'), self.texText('Quiet') ]
    rowLabels = [ self.texText('Day'), self.texText('Night') ]
    title = self.texText('Alfv\\\'en Speed ') + self.texUnit('v')

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title)

    if self.savepath is not None:
      return PW.render(self.savepath + 'va.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # ============================== Contour Plot of Time-Averaged Energy Density
  # ===========================================================================

  def plotubar(self, model=1, driving='J'):
    # 500 seconds is too long compared to the decays we're looking at. 
    tmax = 300

    azms = self.getValues('azm')

    fdrives = self.getValues('fdrive')

    PW = plotWindow(nrows=len(azms), ncols=len(fdrives), colorbar='sym')

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

    driveLabel = 'Current' if driving=='J' else 'Compression'

    colLabels = [ self.texText(format(1e3*f, '.0f') + 'mHz ' + driveLabel) for f in fdrives ]

    title = self.texText( 'Time-Averaged Energy Density: ' + self.texName(model) + self.texUnit('u') )

    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)

    for row, azm in enumerate(azms):
      for col, fdrive in enumerate(fdrives):

        if driving=='J':
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, bdrive=0)
        else:
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, jdrive=0)

        path = self.getPath(azm=azm, model=model, fdrive=fdrive)

        u = self.getArray(path, 'upol') + self.getArray(path, 'utor')

        PW[row, col].setParams( **self.getCoords(path) )

        PW[row, col].setContour( umean )

    if self.savepath is not None:
      return PW.render(self.savepath + 'ubar_' + driving + '_' + str(model) + '.pdf')
    else:
      return PW.render()

  # ===========================================================================
  # =============== Toroidal Poynting Flux Snapshots -- Failing to Propagate In
  # ===========================================================================

  def plotToroidal(self, model=2, driving='B'):
    # 500 seconds is too long compared to the decays we're looking at. 
    tmax = 300

    azms = self.getValues('azm')

    fdrives = self.getValues('fdrive')

    PW = plotWindow(nrows=len(azms), ncols=len(fdrives), colorbar='sym')

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

    driveLabel = 'Current' if driving=='J' else 'Compression'

    colLabels = [ self.texText(format(1e3*f, '.0f') + 'mHz ' + driveLabel) for f in fdrives ]

    title = self.texText( 'Clockwise Toroidal Poynting Flux at 300s: ' + self.texName(model) + self.texUnit('S') )

    PW.setParams(colLabels=colLabels, rowLabels=rowLabels, title=title)

    for row, azm in enumerate(azms):
      for col, fdrive in enumerate(fdrives):

        if driving=='J':
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, bdrive=0)
        else:
          path = self.getPath(azm=azm, model=model, fdrive=fdrive, jdrive=0)

        Stor = self.getArray(path, 'Stor')

        PW[row, col].setParams( **self.getCoords(path) )

        PW[row, col].setContour( Stor[:, :, 299] )

    if self.savepath is not None:
      return PW.render(self.savepath + 'Stor_' + driving + '_' + str(model) + '.pdf')
    else:
      return PW.render()















  # ===========================================================================
  # ================================================================ Sym-H Plot
  # ===========================================================================

  def plotSymh(self):

    PW = plotWindow(nrows=1, ncols=2, colorbar=False, sharelimits=False)

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

    amplitude = np.zeros(nModes, dtype=np.complex)

    def harmonic(m):
      return np.cos( m*np.pi*t / tRange )
#      return np.exp(1j*np.pi*m*t/tRange)

    for m in range(nModes):
      amplitude[m] = np.sum( SYMH*harmonic(m)*dt ) / np.sum( dt*harmonic(m)*harmonic(-m) )

    title = self.texText('Sym-H Frequency Breakdown for June 2013 Storm')
    colLabels = [ self.texText('Sym-H'), self.texText('Sym-H Fourier Modes') ]

    PW[0].setLine(t, SYMH, 'b')

    # Optionally, add a Fourier reconstitution. 
    if False:
      nShow = 10
      colLabels[0] = colLabels[0] + self.texText(' with ' + str(nShow) + ' Fourier Modes')
#      reconstruction = sum( np.real( amplitude[m]*harmonic(-m) ) for m in range(nShow) )
      reconstruction = np.zeros(len(t), dtype=np.float)
      for m in range(nShow):
        reconstruction + np.real( amplitude[m]*harmonic(-m) )

      PW[0].setLine(t, reconstruction, 'r')

    PW.setParams(title=title, collabels=colLabels, ylabel=self.texLabel('symh') )

    PW[0].setParams( xlabel=self.texText('Time (minutes)') )

    PW[1].setParams(xlabel=self.texText('Frequency (mHz)'), ylabel=self.texLabel('logsymh'), ylimits=(-4, None), xlog=True, legend=True)
#    PW[1].setParams(xlabel=self.texText('Frequency (mHz)'), ylabel=self.texLabel('logsymh'), ylimits=(-4, None), xlog=True)

    frequencies = np.array( range(nModes) )*1000./(60*2*tRange)

    PW[1].setLine(frequencies[1:], np.log10( np.abs( amplitude[1:] ) ), 'b')

    # Set limits. They don't quite line up with the data, but that's OK. 
    fMin, fMax = 1e-2, 10
    PW[1].setParams( xlims=(fMin, fMax) )
    f = np.linspace(fMin, fMax, 10000)
    # Let's not do a fit, per se -- let's eyeball the top of the distribution. 
    # Work in units of 20mHz. 
    intercept, slope = np.log(1.e-2), -1
    scale = 20. 
    # Plot the fit. 
    fit = np.exp(intercept) * np.power(f/scale, slope)

    label = '$' + format(np.exp(intercept), '.2f') + self.texText(' nT') + '\cdot \left( \\frac{f}{' + format(scale, '.0f') + self.texText('mHz') + '} \\right)^{' + format(slope, '.1f') + '}' + '$'
    PW[1].setLine( f, np.log10(fit), color='r', label=label)

    if self.savepath is not None:
      return PW.render(self.savepath + 'symh.pdf')
    else:
      return PW.render()




















'''





  # ===========================================================================
  # ==================================================== Snapshot of All Fields
  # ===========================================================================

  def plotB(self, filename=None):
    # Parameters to be held constant. 
    fdrive = 0.007
    azm = 16
    model = 1
    step = -1
    # Rows and columns to be populated. 
    rows = ('x', 'y', 'z')
    cols = ('B', 'E')
    # Create window. Find data. 
    PW = plotWindow(nrows=len(rows), ncols=len(cols), colorbar='sym')
    path = self.getPath(azm=azm, model=model, fdrive=fdrive)
    # Iterate through rows and columns. 
    for row, xyz in enumerate(rows):
      for col, BE in enumerate(cols):
        PW[row, col].setParams( **self.getCoords(path, 'X', 'Z') )
        PW[row, col].setContour( self.getArray(path, BE + xyz)[:, :, step] )
    # Assemble labels and title. 
    rowLabels = [ r for r in rows ]
    colLabels = ( self.texText('Magnetic Field') + self.texUnit('B'),
                  self.texText('Electric Field') + self.texUnit('E') )
    t = self.getArray(path, 't')
    drive = 'Current' if self.getParam(path, 'jdrive')>0 else 'Compression'
    freq = format(1e3*fdrive, '.0f') + 'mHz'
    title = self.texText( self.texName(model) + ' After ' +
                          format(t[step], '.0f') + 's of ' + freq + ' ' +
                          drive )
    PW.setParams(rowlabels=rowLabels, collabels=colLabels, title=title)
    # Show the plot or save it as an image. 
    return PW.render(filename)

  # ===========================================================================
  # ============================== Power Density from J dot E and Poynting Flux
  # ===========================================================================

  def plotG(self, filename=None):
    # Parameters to be held constant. 
    model = 1
    fdrive = 0.007
    azms = (1, 4, 16, 64)

    # Leave out JzEz. It's tiny. 

    PW = plotWindow(nrows=len(azms), ncols=3, colorbar='sym')

    for row, azm in enumerate(azms):
      path = self.getPath(azm=azm, model=model, fdrive=fdrive)

      JEx = self.getArray(path, 'Jx')*self.getArray(path, 'Ex')
      JEy = self.getArray(path, 'Jy')*self.getArray(path, 'Ey')
      JEz = self.getArray(path, 'Jz')*self.getArray(path, 'Ez')

      Stor = self.getArray(path, 'Stor')
      Spol = self.getArray(path, 'Spol')

      dx = self.getArray(path, 'dx0')[:, None]
      dy = azm*self.getArray(path, 'dy0')[:, None]
      dz = self.getArray(path, 'dz0')[:, None]

      dxStor = self.d( Stor[:, 0, :] )/dx
      dxSpol = self.d( Spol[:, 0, :] )/dx

      dyStor = Stor[:, 0, :]/dy
      dySpol = Spol[:, 0, :]/dy

      dzStor = ( Stor[:, 1, :] - Stor[:, 0, :] )/dz
      dzSpol = ( Spol[:, 1, :] - Spol[:, 0, :] )/dz

      PW[row, 0].setParams( **self.getCoords(path, 't', 'lat0') )
      PW[row, 1].setParams( **self.getCoords(path, 't', 'lat0') )
      PW[row, 2].setParams( **self.getCoords(path, 't', 'lat0') )
#      PW[row, 3].setParams( **self.getCoords(path, 't', 'lat0') )

      PW[row, 0].setContour(dxStor + dyStor + dzStor)
      PW[row, 1].setContour(dxSpol + dySpol + dzSpol)
      PW[row, 2].setContour(JEx[:, 0, :] + JEy[:, 0, :])
#      PW[row, 3].setContour(JEz[:, 0, :])

    rowLabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
    colLabels = ( '\\nabla \\cdot' + self.texName('Stor'), 
                  '\\nabla \\cdot' + self.texName('Spol'), 
                  '\\underline{J}_\\bot \\cdot \\underline{E}_\\bot',
#                  'J_\\parallel E_\parallel' 
                 )

    drive = 'Current' if self.getParam(path, 'jdrive')>0 else 'Compression'
    freq = format(1e3*fdrive, '.0f') + 'mHz'

    title = self.texText( self.texName(model) + ' Power Density at R_I with ' + freq + ' ' + drive ) + self.texUnit('JE')

    PW.setParams(collabels=colLabels, rowlabels=rowLabels, title=title)

    return PW.render(filename)

  # ===========================================================================
  # ==================================================== Driving Electric Field
  # ===========================================================================

  def plotH(self, filename=None):
    path = self.getPath()
    PW = plotWindow(nrows=1, ncols=1, colorbar='sym')
    PW.setParams( title=self.texText('Driving Electric Field') +
                        self.texUnit('E') )
    PW.setParams( **self.getCoords(path) )
    PW.setContour( self.getArray(path, 'EyDrive') )
    return PW.render(filename)

'''





# #############################################################################
# ########################################################## Plot Window Object
# #############################################################################

# The Plot Window is a wrapper around PyPlot which handles the high-level
# commands -- subplot spacing, titles, labels, etc. It contains an array of
# Plot Cell objects, which keep track of the actual data. 

class plotWindow:

  # We keep a color axis for the color bar, a title axis, an array of header
  # axes for column labels, an array of side axes for row labels, a side header
  # axis to label the row labels, and a footer axis just in case. Data axes are
  # stored by the individual Plot Cells. 
  cax = None
  fax = None
  hax = None
  sax = None
  shax = None
  tax = None
  # The Plot Window also holds an array of Plot Cells, one for each data axis. 
  cells = None
  # Keep track of the style of color bar, if any, we'll be using. Use 'log' for
  # a log scale, 'sym' for a symmetric log scale, and 'lin' for a linear scale.
  # For no color bar, use False or None. 
  colorbar = None
  # By default, all axes share limits. 
  sharelimits = True
  # Overwrite the default number of colors. 
  ncolors = 8
  # Overwrite the automatically-determined color bar range. 
  zmaxManual = None

  # ---------------------------------------------------------------------------
  # --------------------------------- Initialize Plot Window and Space Out Axes
  # ---------------------------------------------------------------------------

  def __init__(self, ncols=1, nrows=1, colorbar=None, **kargs):
    # Make sure there's nothing lingering from a previous plot. 
    plt.close('all')
    # Set the font to match LaTeX. 
    rc('font', **{'family':'sans-serif', 'sans-serif':['Helvetica'], 
                  'size':'9'})
    rc('text', usetex=True)
    rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}, ' + 
                              '\usepackage{color}')
    # The window width in inches is fixed to match the size of the page. 
    windowWidth = 5.75
    # The window will be broken up into some number of equally-sized tiles.
    # That's the unit we use to specify the relative sizes of plot elements. 
    sideMargin = 40
    cellPadding = 5
    titleMargin = 15
    headMargin = 10
    footMargin = 15
    # The size of each subplot depends on how many columns there are. The total
    # width of the subplot area (including padding) will always be the same.
    # No more than four columns are allowed. 
    cellWidth = {1:175, 2:80, 3:55, 4:40}[ncols]
    # Cells are proportioned to show a dipole plot, which is 10RE wide and 8RE
    # tall, in proper proportion. 
    cellHeight = 4*cellWidth/5
    # Tally up how many tiles we need. 
    tileWidth = ncols*cellWidth + (ncols-1)*cellPadding + 2*sideMargin
    tileHeight = ( nrows*cellHeight + (nrows-1)*cellPadding + titleMargin +
                   headMargin + footMargin )
    # Set the window size in proportion with the number of tiles we need. This
    # ensures that that tiles are square. 
    windowHeight = tileHeight*windowWidth/tileWidth
    # Create the window. Tell it that we want the subplot area to go all the
    # way to the edges, then break that area up into tiles. 
    fig = plt.figure(figsize=(windowWidth, windowHeight), facecolor='white')
    fig.canvas.set_window_title('Tuna Plotter')
    tiles = gridspec.GridSpec(tileHeight, tileWidth)
    plt.subplots_adjust(bottom=0., left=0., right=1., top=1.)
    # Create a lattice of axes and use it to initialize an array of Plot Cells.
    self.cells = np.empty( (nrows, ncols), dtype=object)
    for row in range(nrows):
      for col in range(ncols):
        xpos = sideMargin + col*(cellWidth + cellPadding)
        ypos = titleMargin + headMargin + row*(cellHeight + cellPadding)
        ax = plt.subplot( tiles[ypos:ypos + cellHeight, 
                                xpos:xpos + cellWidth] )
        self.cells[row, col] = plotCell(ax)
    # Space out the title axis. 
    self.tax = plt.subplot( tiles[:titleMargin, sideMargin:-sideMargin] )
    # Space out an array of side axes to hold row labels. 
    self.sax = np.empty( (nrows,), dtype=object)
    for row in range(nrows):
      ypos = titleMargin + headMargin + row*(cellHeight + cellPadding)
      self.sax[row] = plt.subplot( tiles[ypos:ypos + cellHeight, 
                                         :sideMargin - 3*cellPadding] )
    # Space out an array of header axes on the top to hold column labels. 
    self.hax = np.empty( (ncols,), dtype=object)
    for col in range(ncols):
      xpos = sideMargin + col*(cellWidth + cellPadding)
      self.hax[col] = plt.subplot( tiles[titleMargin:titleMargin + headMargin, 
                                         xpos:xpos + cellWidth] )
    # Side header axis, for the row label header. 
    self.shax = plt.subplot( tiles[titleMargin:titleMargin + headMargin, 
                                   :sideMargin - 3*cellPadding] )
    # The title, header, and side axes are for spacing text, not showing data.
    # The axes themselves need to be hidden. 
    self.tax.axis('off')
    self.shax.axis('off')
    [ x.axis('off') for x in self.sax ]
    [ x.axis('off') for x in self.hax ]
    # If we're supposed to have a color bar, space out a narrow axis for it
    # in the right margin. 
    self.colorbar = colorbar
    if self.colorbar:
      self.cax = plt.subplot( tiles[titleMargin + headMargin:-footMargin, 
                                    -sideMargin + cellPadding:-sideMargin +
                                                              3*cellPadding] )
#    # Otherwise, space out room for a legend box. 
#    else:
#      self.cax = plt.subplot( tiles[titleMargin + headMargin:-footMargin, 
#                                    -sideMargin + cellPadding:-sideMargin +
#                                                              3*cellPadding] )
    # We're done setting up the axes. If we were given any other arguments, 
    # send them to the parameter handler. 
    return self.setParams(**kargs)

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------- Adjust Plot Parameters
  # ---------------------------------------------------------------------------

  def setParams(self, **kargs):
    # Keyword parameters used for centering text in axes. 
    targs = {'x':0.5, 'y':0.5, 'horizontalalignment':'center', 
             'verticalalignment':'center'}
    # Address the parameters one at a time. 
    for key, val in kargs.items():
      # Keys are caps insensitive. 
      key = key.lower()
      # Accept a list of strings as column labels. 
      if key=='collabels':
        for col, label in enumerate(val):
          self.hax[col].text(s='$' + label + '$', **targs)
      # Overwrite the default number of colors. 
      elif key=='ncolors':
        self.ncolors = val
      # Accept a list of strings as row labels. 
      elif key=='rowlabels':
        for row, label in enumerate(val):
          self.sax[row].text(s='$' + label + '$', **targs)
      # Sometimes, we may want to label the row labels. 
      elif key=='rowlabellabel':
        self.shax.text(s='$' + val + '$', **targs)
      # By default, axis limits are common to all plots. 
      elif key=='sharelimits':
        self.sharelimits = bool(val)
      # Accept a string as the window supertitle. 
      elif key=='title':
        self.tax.text(s='$' + val + '$', fontsize=12, **targs)
      # Only the bottom x axes get labels. 
      elif key=='xlabel':
        [ cell.setParams(xlabel=val) for cell in self.cells[-1, :] ]
      # Only the leftmost y axes get labels. 
      elif key=='ylabel':
        [ cell.setParams(ylabel=val) for cell in self.cells[:, 0] ]
      # Overwrite the automatically-determined color bar range. 
      elif key=='zmax':
        self.zmaxManual = val
      # Any other parameters get sent to the cells. 
      else:
        [ cell.setParams( **{key:val} ) for cell in self.cells.flatten() ]
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------------ Add Data
  # ---------------------------------------------------------------------------

  # The Plot Window doesn't actually handle any data. Individual cells should
  # instead be accessed as array entries. 
  def __getitem__(self, index):
    if isinstance(index, int):
      return self.cells.flatten()[index]
    else:
      return self.cells[index]

  # If the window gets passed a contour, just send it along to each cell. 
  def setContour(self, *args, **kargs):
    return [ cell.setContour(*args, **kargs) for cell in self.cells.flatten() ]

  # If the window gets passed a line, just send it along to each cell. 
  def setLine(self, *args, **kargs):
    return [ cell.setLine(*args, **kargs) for cell in self.cells.flatten() ]

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Get Cell Extrema
  # ---------------------------------------------------------------------------

  # We standardize color levels and axis ranges across all cells.

  def xmax(self):
    return max( cell.xmax() for cell in self.cells.flatten() )

  def ymax(self):
    return max( cell.ymax() for cell in self.cells.flatten() )

  def zmax(self):
    return max( cell.zmax() for cell in self.cells.flatten() )

  # Looking at minima is tricky, since some of them may be None, which is
  # smaller than any number. 

  def xmin(self):
    xmn = [ cell.xmin() for cell in self.cells.flatten() ]
    return None if max(xmn) is None else min( x for x in xmn if x is not None )

  def ymin(self):
    ymn = [ cell.ymin() for cell in self.cells.flatten() ]
    return None if max(ymn) is None else min( y for y in ymn if y is not None )

  def zmin(self):
    zmn = [ cell.zmin() for cell in self.cells.flatten() ]
    return None if max(zmn) is None else min( z for z in zmn if z is not None )

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------- Render Window
  # ---------------------------------------------------------------------------

  # Once all of the contours are loaded, we can standardize the plot domain and
  # color levels. 
  def render(self, filename=None):
    # Use the most extreme x and y values to set the plot domain. Snap to
    # integers, but don't round based on tiny bits of numerical noise. 
    if self.sharelimits:
      xmin = np.floor( float( format(self.xmin(), '.4e') ) )
      xmax = np.ceil( float( format(self.xmax(), '.4e') ) )
      ymin = np.floor( float( format(self.ymin(), '.4e') ) )
      ymax = np.ceil( float( format(self.ymax(), '.4e') ) )
      self.setParams( xlims=(xmin, xmax), ylims=(ymin, ymax) )

      print 'ymin, ymax = ', ymin, ymax

      # Try to de-cramp cramped plots a bit. 
      if ymin==2 and ymax==6:
        for cell in self.cells.flatten():
          cell.setParams( yticks=(2, 3, 4, 5, 6), yticklabels=('$2$', '', '$4$', '', '$6$') )

      if ymin==2 and ymax==10:
        for cell in self.cells.flatten():
          cell.setParams( yticks=(2, 4, 6, 8, 10), yticklabels=('$2$', '', '$6$', '', '$10$') )

      # Only the leftmost cells get y axis labels and tick labels. 
      for cell in self.cells[:, 1:].flatten():
        cell.setParams( ylabel='', yticklabels=() )
      # Try to de-cramp cramped plots a bit. 
      if xmin==1 and xmax==300:
        for cell in self.cells.flatten():
          cell.setParams( xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
      # Only the bottom cells get x axis labela and tick labels. 
      for cell in self.cells[:-1, :].flatten():
        cell.setParams( xlabel='', xticklabels=() )
    # Sometimes, oddly, different things might go on different axes. 
    else:
      [ cell.setParams(rightlabel=True) for cell in self.cells[:, -1] ]
    # Use the most extreme contour value among all plots to set the color bar. 
    if self.zmaxManual is not None:
      colors = plotColors(zmax=self.zmaxManual, cax=self.cax, colorbar=self.colorbar, ncolors=self.ncolors)
    else:
      colors = plotColors(zmax=self.zmax(), cax=self.cax, colorbar=self.colorbar, ncolors=self.ncolors)
    # Send the color params to each cell. 
    [ cell.render(**colors) for cellRow in self.cells for cell in cellRow ]
    # If given a filename, save the image. 
    if filename is not None:
      print 'Saving ' + filename
      return plt.savefig(filename)
    # Otherwise, display it. 
    else:
      return plt.show()

# #############################################################################
# ############################################################ Plot Cell Object
# #############################################################################

class plotCell:

  # If this cell contains a contour, we'll need to hold the spatial coordinates
  # and the data values. We also keep room for any arguments for contourf. 
  x, y, z, kargs = None, None, None, None
  # A plot can have any number of lines drawn on it. Those will be stored here. 
  lines = None
  # If we manually set the axis limits, we want to ignore the automatically-set
  # limits that come down the line later. 
  xlims, ylims = (None, None), (None, None)
  # If we're on a log scale, we need different rules for placing ticks. 
  xlog, ylog = False, False
  # Keep track if we're supposed to be tracing the outline of our domain, such
  # as if the data is dipole-shaped. 
  outline = False
  # Cells can be small. Let's try to keep the number of ticks under control.
  nxticks, nyticks = 3, 4
  # The plots are squished pretty close together. Sometimes we want to squeeze
  # an axis in on the right. 
  rightlabel = False
  # By default, have no legend. 
  legend = False

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Initialize Cell
  # ---------------------------------------------------------------------------

  def __init__(self, ax):
    self.ax = ax
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Set Cell Parameters
  # ---------------------------------------------------------------------------

  def setParams(self, **kargs):
    # Scroll through the parameters one at a time. 
    for key, val in kargs.items():
      # Keys are caps insensitive. 
      key = key.lower()
      # Turn on the legend. 
      if key=='legend' and bool(val) is True:
        self.legend = True
      # Sometimes we have to finagle with the number of ticks. 
      elif key=='nxticks':
        self.nxticks = val
      elif key=='nyticks':
        self.nyticks = val
      # Draw an outline around the plot contents. 
      elif key=='outline':
        self.outline = val
      # Squeeze in a label on the right axis. 
      elif key=='rightlabel' and bool(val) is True:
        self.ax.yaxis.tick_right()
        self.ax.yaxis.set_label_position('right')
      # Horizontal axis coordinate. 
      elif key=='x':
        self.x = val
      # Label the horizontal axis. 
      elif key=='xlabel':
        self.ax.set_xlabel('' if not val else '$' + val + '$')
      # Set horizontal axis domain. 
      elif key.startswith('xlim'):
        # If the limits are set manually, we want to ignore the automatic
        # limits that come down the line later. 
        self.xlims = ( val[0] if self.xlims[0] is None else self.xlims[0], 
                       val[1] if self.xlims[1] is None else self.xlims[1] )
      # Set log horizontal scale. 
      elif key=='xlog' and val is True:
        self.ax.set_xscale('log')
        self.xlog = True
      # Set the horizontal axis tick labels. 
      elif key=='xticklabels':
        self.ax.set_xticklabels(val)
      # Set the horizontal axis ticks. 
      elif key=='xticks':
        self.ax.set_xticks(val)
      # Vertical axis coordinate. 
      elif key=='y':
        self.y = val
      # Label the vertical axis. 
      elif key=='ylabel':
        self.ax.set_ylabel('' if not val else '$' + val + '$')
      # Set the vertical axis domain. 
      elif key.startswith('ylim'):
        # If the limits are set manually, we want to ignore the automatic
        # limits that come down the line later. 
        self.ylims = ( val[0] if self.ylims[0] is None else self.ylims[0], 
                       val[1] if self.ylims[1] is None else self.ylims[1] )
      # Set log vertical scale. 
      elif key=='ylog' and val is True:
        self.ax.set_yscale('log')
        self.ylog = True
      # Set the vertical axis tick labels. 
      elif key=='yticklabels':
        self.ax.set_yticklabels(val)
      # Set the vertical axis ticks. 
      elif key=='yticks':
        self.ax.set_yticks(val)
      # Report any unfamiliar parameters. 
      else:
        print 'WARNING: Unknown param ', key, ' = ', val
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------- Set Cell Data
  # ---------------------------------------------------------------------------

  def setContour(self, *args, **kargs):
    # Store any keyword parameters meant for the contourf call. 
    self.kargs = kargs
    # Accept the contour with or without its spatial coordinates. 
    if len(args)==1:
      self.z = args[0]
    elif len(args)==3:
      self.x, self.y, self.z = args
    # If we're passed a weird number of arguments, bail. 
    else:
      print 'ERROR: Illegal number of arguments to plotCell.setContour '
      exit()
    return

  def setLine(self, *args, **kargs):
    # Initialize line list. 
    if self.lines is None:
      self.lines = []
    # Store this line in the list. Worry about the extra arguments later. 
    self.lines.append( (args, kargs) )
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------- Report Cell Extrema
  # ---------------------------------------------------------------------------

  # Cells all share a color bar (for contour plots) and axis limits (for line
  # plots). To manage that, the Plot Window asks each cell for its extrema. 

  def xmax(self):
    amax = None if self.x is None else np.max(self.x)
    if self.lines is None:
      return amax
    else:
      return max( amax, max( max( args[0] ) for args, kargs in self.lines ) )

  def ymax(self):
    amax = None if self.y is None else np.max(self.y)
    if self.lines is None:
      return amax
    else:
      return max( amax, max( max( args[1] ) for args, kargs in self.lines ) )

  def zmax(self):
    return None if self.z is None else np.max(self.z)

  # Minima are tricky, since None counts as smaller than any number. 

  def xmin(self):
    amin = None if self.x is None else np.min(self.x)
    if self.lines is None:
      return amin
    else:
      lmin = min( min( args[0] ) for args, kargs in self.lines )
      return lmin if amin is None else min(amin, lmin)

  def ymin(self):
    amin = None if self.y is None else np.min(self.y)
    if self.lines is None:
      return amin
    else:
      lmin = min( min( args[1] ) for args, kargs in self.lines )
      return lmin if amin is None else min(amin, lmin)

  def zmin(self):
    return None if self.z is None else np.min(self.z)

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Render Plot Cell
  # ---------------------------------------------------------------------------

  def render(self, **colors):
    # If this cell has a contour, lay that down first. 
    if self.z is not None:
      # Use the color params we were passed, but allow keyword arguments from
      # the contour call to overwrite them. 
      kargs = dict( colors.items() + self.kargs.items() )
      self.ax.contourf(self.x, self.y, self.z, **kargs)
    # Optionally, draw the outline of the data. 
    if self.outline:
      [ self.setLine(self.x[i, :], self.y[i, :], 'k') for i in (0, -1) ]
      [ self.setLine(self.x[:, k], self.y[:, k], 'k') for k in (0, -1) ]
    # Draw any lines. 
    if self.lines is not None:
      [ self.ax.plot(*args, **kargs) for args, kargs in self.lines ]
    # Set axis limits. 
    self.ax.set_xlim(self.xlims)
    self.ax.set_ylim(self.ylims)

    # Optionally, turn on the legend. 
    if self.legend:
      self.ax.legend()
#      self.ax.legend(loc="upper left", bbox_to_anchor=[0, 1],
#           ncol=2, shadow=True, title="Legend", fancybox=True)




#    # Try to keep the axis ticks under control. 
#    ymin, ymax = self.ylims
#    if 2 <= (ymax - ymin) < 5:
#      ticks = range(int(ymin), int(ymax)+1)
#      print self.ax.get_yticklabels()
#      if self.ax.get_yticklabels():
#        labels = [ '$' + str(tick) + '$' for tick in ticks ]
#      else:
#        labels = []
#      self.ax.set_yticks(ticks)
#      self.ax.set_yticklabels(labels)

    # These subplots can get cramped. Let's reduce the number of ticks. 
    if not self.xlog:
      if self.xlims==(0, 300):
        pass
#        print 'xlims are 0 to 300'
#        self.ax.set_xticks( (0, 100, 200, 300) )
#        print self.ax.get_xticklabels()
#        if all( label=='' for label in self.ax.get_xticklabels() ):
#          self.ax.set_xticklabels( ('$0$', '', '', '$300$') )
      else:
        self.ax.xaxis.set_major_locator( plt.MaxNLocator(self.nxticks,
                                                       integer=True) )

    if self.ylog:
      pass
    else:
      self.ax.yaxis.set_major_locator( plt.MaxNLocator(self.nyticks, 
                                                       integer=True) )

    self.ax.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
    self.ax.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')

    self.ax.yaxis.get_majorticklabels()[0].set_verticalalignment('bottom')
    self.ax.yaxis.get_majorticklabels()[-1].set_verticalalignment('top')

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

  def __init__(self, zmax, cax, colorbar=None, ncolors=None):
    # Some plots don't have contours. 
    if not zmax or not cax or not colorbar:
      return dict.__init__(self, {})
    # Store the data scale so that it can be hard-wired into our normalization
    # functions. 
    self.colorbar = colorbar
    self.ncolors = ncolors if ncolors is not None else 8
    self.nticks = self.ncolors - 1
    # Assemble the keyword parameters in a temporary dictionary. We'll then use
    # the dictionary constructor to build this object based on it. 
    temp = {}
    # Determine location of contour color levels and color bar ticks. 
    if self.colorbar=='log':
      temp['ticks'], temp['levels'] = self.logTicksLevels(zmax)
      temp['norm'] = LogNorm()
    elif self.colorbar=='sym':
      temp['ticks'], temp['levels'] = self.symTicksLevels(zmax)
      temp['norm'] = Normalize()

    elif self.colorbar=='phase':
      temp['ticks'], temp['levels'] = self.phaseTicksLevels(zmax)
      temp['norm'] = Normalize()

    else:
      temp['ticks'], temp['levels'] = self.linTicksLevels(zmax)
      temp['norm'] = Normalize()
    # Rework the color map to match the normalization of our ticks and levels. 
    temp['cmap'] = self.getCmap()
    # Draw the color bar. 
    self.setColorbar(cax, **temp)
    # Become a dictionary of color parameters to be used by the contour plots. 
    return dict.__init__(self, temp)

  # ---------------------------------------------------------------------------
  # ----------------------------------- Tick Locations and Contour Color Levels
  # ---------------------------------------------------------------------------

  def linTicksLevels(self, zmax):
    # We put zmax at the top of the top color level. Ticks go at the middle of
    # color levels, as they do with symlog plots. 
    self.zmax = zmax
    levels = np.linspace(-self.zmax, self.zmax, self.ncolors)
    ticks = 0.5*( levels[1:] + levels[:-1] )
    # Make sure that the middle tick is exactly zero. 
    ticks[ len(ticks)/2 ] = 0.
    return ticks, levels

  def phaseTicksLevels(self, zmax):
    self.zmax = np.pi/2
    self.zmin = 0
    self.nticks = self.ncolors
    ticks = np.linspace(self.zmin, self.zmax, self.nticks)
    levels = np.linspace(self.zmin, self.zmax, self.ncolors)
    return ticks, levels

  def logTicksLevels(self, zmax):
    # Ticks are located at powers of ten. Color levels are centered on ticks. 
    power = np.ceil(np.log10(zmax) - 0.5)
    # Each color spans a factor of root ten. This is in contrast to the
    # symmetric log scale, where each color was a whole order of magnitude. The
    # goal is to, for each, have the same number of colors and the same number
    # of orders of magnitude. 
    # Symetric log scale with 7 colors will have three positive powers of ten,
    # three negative powers, and zero. The log scale will just have three
    # positive powers. Anything below there will automatically show 0, though
    # it won't be marked explicitly on the color bar. 
    self.zmax = 10.**(power + 0.25)
    self.zmin = self.zmax/10**(self.nticks/2 + 0.5)
    ticks = [ 10**(power - 0.5*i) for i in range(self.nticks) ]
    logMin, logMax = np.log10(self.zmin), np.log10(self.zmax)
    levels = np.logspace(logMin, logMax, self.ncolors)
    return ticks, levels

  def symTicksLevels(self, zmax):
    # Ticks are located at powers of ten. Color levels are centered on ticks. 
    power = np.ceil( np.log10( zmax/np.sqrt(10.) ) )
    self.zmax = np.sqrt(10.)*10**power
    # A tick at zero, then one per order of magnitude. 
    norders = (self.nticks - 1)/2
    pticks = [ 10**(power - i) for i in range(norders) ]
    plevels = [ np.sqrt(10.)*10**(power - i) for i in range(norders + 1) ]
    # For uniform tick spacing, the log cutoff is a factor of 10 below the
    # lowest positive tick. That is, it's where the next tick would be, if we
    # had one more tick. 
    self.zmin = min(pticks)/10.
    ticks = sorted( pticks + [0] + [ -tick for tick in pticks ] )
    levels = sorted( plevels + [ -level for level in plevels ] )
    return ticks, levels

  # ---------------------------------------------------------------------------
  # ----------------------------------------------- Data Interval Normalization
  # ---------------------------------------------------------------------------

  # Map from the unit interval to the data scale via linear scale. 
  def linNorm(self, x):
    return self.zmax*(2*x - 1)

  # Map from the data scale to the unit interval via linear scale. 
  def linMron(self, x):
    return 0.5 + 0.5*x/self.zmax

  # Map from the unit interval to the data scale via linear scale. 
  def phaseNorm(self, x):
    return x*self.zmax

  # Map from the data scale to the unit interval via linear scale. 
  def phaseMron(self, x):
    return x/self.zmax

  # Map from the unit interval to the data scale via log scale. 
  def logNorm(self, x):
    return self.zmin*(self.zmax/self.zmin)**x

  # Map from the log scaled data scale to the unit interval. 
  def logMron(self, x):
    return np.log10(x/self.zmin)/np.log10(self.zmax/self.zmin)

  # Map from the unit interval to the data scale via symmetric log scale. 
  def symNorm(self, x):
    if x>0.5:
      return self.zmin*(self.zmax/self.zmin)**(2*x - 1)
    elif x<0.5:
      return -self.zmin*(self.zmax/self.zmin)**(1 - 2*x)
    else:
      return 0

  # Map from the symmetric log scaled data scale to the unit interval. 
  def symMron(self, x):
    if x>self.zmin:
      return 0.5 + 0.5*np.log10(x/self.zmin)/np.log10(self.zmax/self.zmin)
    elif x<-self.zmin:
      return 0.5 - 0.5*np.log10(-x/self.zmin)/np.log10(self.zmax/self.zmin)
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
#      cmap = plt.get_cmap('seismic')
#      colors = cmap( np.linspace(0.5, 1, cmap.N // 2) )
#      return LSC.from_list('Upper Half', colors)

#      norm = self.logNorm
#      return plt.get_cmap('YlOrRd')
      return None

    elif self.colorbar=='sym':
      norm = self.symNorm
    elif self.colorbar=='phase':
      # The physics machines at the U use an old version of Matplotlib... 1.0.1. Cubehelix was added in 1.5. It can also be obtained here: 
      # https://github.com/jradavenport/cubehelix/blob/master/cubehelix.py
      return None


    else:
      norm = self.linNorm
    # Get a fine sampling of the color map on the unit interval. 
    N = 1000
    unitInterval = [ i/(N - 1.) for i in range(N) ]
    cmap = plt.get_cmap('seismic')
    rgb = [ cmap(u) for u in unitInterval ]

#    if self.colorbar=='log':
#      rgb = [ cmap( 0.5 + 0.5*u ) for u in unitInterval ]



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
  def setColorbar(self, cax, **colorParams):
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

      ColorbarBase(cax, boundaries=colorParams['levels'],
                   ticks=colorParams['ticks'], norm=colorParams['norm'],
                   cmap=colorParams['cmap'])
      cax.set_yticklabels( [ fmt(t) for t in colorParams['ticks'] ] )
      return

    elif self.colorbar=='sym':
      norm, mron, fmt = self.symNorm, self.symMron, self.symFormatter

    elif self.colorbar=='phase':
      norm, mron, fmt = self.phaseNorm, self.phaseMron, self.phaseFormatter


    else:
      norm, mron, fmt = self.linNorm, self.linMron, self.linFormatter
    # Draw the contour. 
    Z = np.vectorize(norm)(Y)
    cax.contourf( X, Y, Z, **colorParams)
    # Place the ticks appropriately on the unit interval (Y axis). 
    cax.set_yticks( [ mron(t) for t in colorParams['ticks'] ] )
    # Format tick names nicely. 
    cax.set_yticklabels( [ fmt(t) for t in colorParams['ticks'] ] )
    # Put the color bar ticks on the right, get rid of the ticks on the bottom,
    # and hide the little notches in the color bar. 
    cax.yaxis.tick_right()
    cax.set_xticks( [] )
    cax.tick_params( width=0 )
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
    # If our numbers are around order unity, the top tick should show two
    # significant figures, and the rest should match that decimal place. 
    elif 1e-3<self.zmax<1e3:
      power = int( format(self.zmax, '.1e').split('e')[-1] )
      digs = max(0, 1 - power)
      sign = '' if x<0 else '+'
      return '$' + sign + format(x, '.' + str(digs) + 'f') + '$'
    else:
      # Cast the number in scientific notation. 
      s = format(x, '.1e').replace('e', ' \\cdot 10^{') + '}'
      # If the number is positive, throw a plus sign on there. 
      s = '+ ' + s if x>0 else s
      # Before returning, get rid of any extra digits in the exponent. 
      return '$' + s.replace('+0', '').replace('-0', '-') + '$'

  def phaseFormatter(self, x):
    # Zero is always zero. 
    if x==0:
      return '$0$'
    else:

      numer, denom = (x/np.pi).as_integer_ratio()
      if numer==1:
        return '${\\displaystyle \\frac{\\pi}{' + str(denom) + '}}$'
      else:
        return '${\\displaystyle \\frac{' + str(numer) + ' \\pi}{' + str(denom) + '}}$'


#      print x, '\t->', (x/np.pi).as_integer_ratio()
#      return '$' + format(x/np.pi, '.2f') + '\\pi$'



  def logFormatter(self, x):
    # Zero is always zero. 
    if x==0:
      return '$0$'
    # Half-power ticks don't get labels. 
    elif format(x, '.1e').startswith('3'):
      return ''
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

# Convert a string to a complex or float. 
def com(x):
  if ',' not in x:
    return float(x)
  else:
    # Shave off the parentheses then split into real and imaginary parts. 
    re, im = x[1:-1].split(',')
    return (float(re) + float(im)*1j)

# Given kargs full of lists, return a list of kargs. Or just one of them. 
def loopover(pick=False, **kargs):
  lo = [ [] ]
  for key, vals in kargs.items():
    lo = [ l + [ (key, v) ] for l in lo for v in vals ]
  if pick is True:
    return [ dict( choice(lo) )  ]
  else:
    return [ dict(l) for l in lo ]

# Turn a string into a float or integer. 
def num(x):
  return int( float(x) ) if float(x)==int( float(x) ) else float(x)

# #############################################################################
# ################################################################ Array Reader
# #############################################################################

# Read a file of values into an array. We expect the first line to give the
# array dimensions. The rest of the file should just be real or complex
# numbers, and is agnostic to whitespace and newlines. 
def readArray(filename):
  # The out prefix is just an older convention for dat files, Fortran output.
  # We ultimately want to use the Python data format, pickles. 
  name = filename[ :filename.rfind('.') ]
  datname, outname, pklname = name + '.dat', name + '.out', name + '.pkl'

  # Allow epsp.dat, but rename to epsPerp.dat. 
  if os.path.exists( datname.replace('epsPerp', 'epsp') ):
    outname = outname.replace('epsPerp.out', 'epsp.dat')

  # Allow jz.out, but rename to Jz.dat (change in capitalization). 
  outname = outname.replace('Jz', 'jz')

  # If we see a old file (.out) move it to the new convention (.dat). 
  if os.path.isfile(outname) and not os.path.exists(datname):
    os.rename(outname, datname)
    print 'Renamed ' + basename(outname) + ' to ' + basename(datname)

  # If a pickle is available, read that instead of parsing the Fortran output. 
  if os.path.isfile(pklname):
    print 'Reading ' + basename(pklname) + ' ... ',
    stdout.flush()
    # Note how long it takes to read the file. 
    start = time()
    with open(pklname, 'rb') as handle:
      inputArray = pickle.load(handle)
    print format(time() - start, '5.1f') + 's' + ' ... ' + by(inputArray.shape)
    return inputArray
  # If the pickle doesn't exist yet, parse the Fortran output. 
  elif os.path.isfile(datname):
    print 'Reading ' + basename(datname) + ' ... ',
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
    inputArray = np.transpose( np.reshape( vals, dims[::-1] ) )
    # Check the array dimensions. If it's not full, the run may have crashed. 
    # In that case we return only the time steps that exist. 
    actualDims = dims[:-1] + [ np.int( i/np.prod(dims[:-1]) ) ]
    if dims!=actualDims:
      inputArray = inputArray[ ..., :actualDims[-1] ]
    # Report how long this all took. 
    print format(time() - start, '5.1f') + 's'
    # Dump a pickle for next time. 
    print '\tCreating ' + basename(pklname) + ' ... ',
    stdout.flush()
    start = time()
    with open(pklname, 'wb') as handle:
      pickle.dump(inputArray, handle, protocol=-1)
    print format(time() - start, '5.1f') + 's'
    return inputArray
  # If the pickle doesn't exist and there's no Fortran output, return nothing. 
  else:
    print 'WARNING: ' + datname + ' not found. '

# #####################################################################
# ################################################### For Importability
# #####################################################################

if __name__=='__main__':
  main()

