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

#import gc
#import matplotlib
## Change matplotlib settings to allow use over SSH without X forwarding. 
#if 'DISPLAY' not in os.environ or os.environ['DISPLAY'] is '':
#  matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from matplotlib.ticker import FuncFormatter
#from matplotlib import gridspec, rc
## To easily draw a half-black-half-white semicircle. 
#from matplotlib.patches import Wedge
#import numpy as np

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

#from matplotlib.ticker import FuncFormatter
#import matplotlib.ticker as ticker

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  TP = tunaPlotter()

  return TP.plot()

# #############################################################################
# ######################################################### Tuna Plotter Object
# #############################################################################

# The Tuna Plotter is a wrapper around the Plot Window class which includes
# utilities specific to Tuna. (The Plot Window itself is fairly general.)
class tunaPlotter:

  paths = None

  # ===========================================================================
  # ============================================================= Path Handling
  # ===========================================================================

  # Given a path, find all data directories. 
  def setPaths(self, path):
    self.paths = []
    for root, dirs, files in os.walk(path):
      if 'params.in' in files:
        self.paths.append( os.path.abspath(root) + '/' )
    return

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

  # Given some keyword arguments, filter self.paths. 
  def getPath(self, **kargs):
    # Start with all the paths we have, then weed out any that don't match. 
    paths = self.paths
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


  # ===========================================================================
  # ==================================================== LaTeX Helper Functions
  # ===========================================================================

  def texText(self, x):
    return '\\mathrm{' + x.replace(' ', ' \\; ') + '}'

  def texName(self, x):
    names = {
             # Spell out what each model means. 
             1:self.texText('Active Dayside '),
             2:self.texText('Quiet Dayside '),
             3:self.texText('Active Nightside '),
             4:self.texText('Quiet Nightside '),
             # Names for fields, axes, etc. 
             'Bx':'B_x', 
             'By':'B_y', 
             'Bz':'B_z', 
             'X':'X', 
             'Z':'Z'
            }
    return '?' if x not in names else names[x]

  def texReIm(self, x):
    if x in ('Ex', 'By', 'Ez', 'Jz'):
      return ' \\mathbb{I}\\mathrm{m}\\;\\; '
    elif x in ('Bx', 'Ey', 'Bz'):
      return ' \\mathbb{R}\\mathrm{e}\\;\\; '
    else:
      return ''

  def texTime(self, x):
    # Time is a float, but we don't need a '.0' appended to everything. 
    strx = str( int(x) ) if int(x)==x else str(x)
    return self.texText(' at ' + strx + 's')

  def texUnit(self, x):
    units = {
             'B':'nT',
             'Bx':'nT', 
             'By':'nT', 
             'Bz':'nT', 
             'nT':'nT',
             's':'s',
             't':'s',
             'X':'R_E', 
             'Z':'R_E'
            }
    return self.texText( ' (' + ( '?' if x not in units else units[x] ) + ')' )

  def texLabel(self, x):
    return self.texReIm(x) + self.texName(x) + self.texUnit(x)

  # ===========================================================================
  # =============================================================== Data Access
  # ===========================================================================


  def getArray(self, path, name):

    if name in ('Bx', 'By', 'Bz', 'r', 'q', 't'):
      return readArray(path + name + '.dat')
    elif name=='X':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      return r*np.sin(q)
    elif name=='Z':
      r, q = self.getArray(path, 'r'), self.getArray(path, 'q')
      return r*np.cos(q)

    else:
      print 'ERROR: Not sure how to get ' + name
      exit()






  # ===========================================================================
  # ====================================================== Coordinate Shorthand
  # ===========================================================================

  def getCoords(self, path, x='X', y='Z'):

    # This function returns a dictionary of keyword arguments meant to be
    # plugged straight into plotWindow.setParams(). 

    params = {'x':self.getArray(path, x), 
              'xlabel':self.texLabel(x),
              'y':self.getArray(path, y), 
              'ylabel':self.texLabel(y)}

    return params


  # ===========================================================================
  # ============================================================= Plot Assembly
  # ===========================================================================

  def plot(self, path='/export/scratch/users/mceachern/parallel_test'):
    # Starting at the given path, find all directories with data. 
    self.setPaths(path)
    # Let's plot magnetic field components at different m values. 
    azms = (1, 8, 64)
    names = ('Bx', 'By', 'Bz')
    # Initialize the plot window. 
    PW = plotWindow(nrows=len(azms), ncols=len(names), colorbar='sym')
    # Choose a model and a time step. 
    model = 1
    step = 99

    # Create row and column labels. 
    rowLabels = [ 'm = ' + str(azm) for azm in azms ]
    colLabels = [ self.texLabel(name) for name in names ]
    PW.setParams(rowLabels=rowLabels, colLabels=colLabels)

    # Loop through the rows and columns. 
    for row, azm in enumerate(azms):
      for col, name in enumerate(names[:-1]):

        # Find the path that matches the parameters we want for this cell. 
        path = self.getPath(inertia=1, azm=azm, model=model)

#        x = self.getArray(path, 'X')
#        y = self.getArray(path, 'Z')
        z = self.getArray(path, name)

        # Fields should be overwhelmingly real or complex. It's safe to label
        # the column based on the first row. 
        if row==0:
          colLabels[col] = z.phase + colLabels[col]

        PW[row, col].setParams( **self.getCoords(path, 'X', 'Z') )
        PW[row, col].setContour( z[:, :, step] )

        t = self.getArray(path, 't')

    title = self.texName(model) + self.texText('Magnetic Field') + self.texTime( t[step] ) + self.texUnit('nT')

    PW.setParams(title=title)


#    xlabel = texLabel('X')
#    ylabel = texLabel('Z')

#    title = conditions + description + time + (units)
#    title = text('Plot Title') + unit('S')
#    rowLabels = [ text( 'Row ' + str(row) ) for row in range(nrows) ]
#    colLabels = [ text( 'Column ' + str(col) ) for col in range(ncols) ]

#    for row in range(nrows):
#      for col in range(ncols):
#        z = np.random.randn(20, 20)
#        PW[row, col].setContour(x, y, z)

    return PW.render()


  '''





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
'''











# #############################################################################
# ########################################################## Plot Window Object
# #############################################################################

# The Plot Window is a wrapper around PyPlot which handles the high-level
# commands -- subplot spacing, titles, labels, etc. It contains an array of
# Plot Cell objects, which keep track of the actual data. 

class plotWindow:

  # We keep a color axis for the color bar, a title axis, an array of header
  # axes for column labels, an array of side axes for row labels, and a footer
  # axis just in case. Data axes are stored by the individual Plot Cells. 
  cax = None
  fax = None
  hax = None
  sax = None
  tax = None
  # The Plot Window also holds an array of Plot Cells, one for each data axis. 
  cells = None
  # Keep track of the style of color bar, if any, we'll be using. Use 'log' for
  # a log scale, 'sym' for a symmetric log scale, and 'lin' for a linear scale.
  # For no color bar, use False or None. 
  colorbar = None

  # ---------------------------------------------------------------------------
  # --------------------------------- Initialize Plot Window and Space Out Axes
  # ---------------------------------------------------------------------------

  def __init__(self, ncols=1, nrows=1, colorbar=None, **kargs):
    # Make sure there's nothing lingering from a previous plot. 
    plt.close('all')
    # Set the font to match LaTeX. 
    rc('font', **{'family':'sans-serif', 'sans-serif':['Helvetica'], 
                  'size':'11'})
    rc('text', usetex=True)
    rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}')
    # The window width in inches is fixed to match the size of the page. 
    windowWidth = 5.75
    # The window will be broken up into some number of equally-sized tiles.
    # That's the unit we use to specify the relative sizes of plot elements. 
    sideMargin = 40
    cellPadding = 5
    titleMargin = 25
    headMargin = 10
    footMargin = 20
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
    # The title, header, and side axes are for spacing text, not showing data.
    # The axes themselves need to be hidden. 
    self.tax.axis('off')
    [ x.axis('off') for x in self.sax ]
    [ x.axis('off') for x in self.hax ]
    # If we're supposed to have a color bar, space out a narrow axis for it
    # in the right margin. 
    self.colorbar = colorbar
    if colorbar:
      self.cax = plt.subplot( tiles[titleMargin + headMargin:-footMargin, 
                                    -sideMargin + cellPadding:-sideMargin +
                                                              3*cellPadding] )
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
      # Accept a list of strings as row labels. 
      elif key=='rowlabels':
        for row, label in enumerate(val):
          self.sax[row].text(s='$' + label + '$', **targs)
      # Accept a string as the window supertitle. 
      elif key=='title':
        self.tax.text(s='$' + val + '$', fontsize=20, **targs)
      # Only the bottom x axes get labels. 
      elif key=='xlabel':
        [ cell.setParams(xlabel=val) for cell in self.cells[-1, :] ]
      # Only the leftmost y axes get labels. 
      elif key=='ylabel':
        [ cell.setParams(ylabel=val) for cell in self.cells[:, 0] ]
      # Report if we see any parameter we're not prepared for. 
      else:
        print 'WARNING: Unknown param ', key, ' = ', val
    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------------------ Add Data
  # ---------------------------------------------------------------------------

  # The Plot Window doesn't actually handle any data. Individual cells should
  # instead be accessed as array entries. 
  def __getitem__(self, index):
    return self.cells[index]

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
    # Remove the text from all axes except the edge ones. 
    for cell in self.cells[:-1, :].flatten():
      cell.setParams( xlabel='', xticklabels=() )
    for cell in self.cells[:, 1:].flatten():
      cell.setParams( ylabel='', yticklabels=() )
    # Use the most extreme x and y values to set the plot domain. This should
    # work for both line and contour plots. 
    xlm = [ np.floor( self.xmin() ), np.ceil( self.xmax() ) ]
    ylm = [ np.floor( self.ymin() ), np.ceil( self.ymax() ) ]
    [ cell.setParams(xlims=xlm, ylims=ylm) for cell in self.cells.flatten() ]
    # Use the most extreme contour value among all plots to set the color bar. 
    colors = plotColors(zmax=self.zmax(), cax=self.cax, colorbar=self.colorbar)
    [ cell.render(**colors) for cellRow in self.cells for cell in cellRow ]



    return plt.show()
#    return plt.savefig('/home/user1/mceachern/Desktop/plots/test.pdf')


'''
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
'''


# #############################################################################
# ############################################################ Plot Cell Object
# #############################################################################

class plotCell:

  # If this cell contains a contour, we'll need to hold the spatial coordinates
  # and the data values. We also keep room for any arguments for contourf. 
  x, y, z, kargs = None, None, None, None
  # A plot can have any number of lines drawn on it. Those will be stored here. 
  lines = None

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
      # Horizontal axis coordinate. 
      if key=='x':
        self.x = val
      # Label the horizontal axis. 
      elif key=='xlabel':
        self.ax.set_xlabel('' if not val else '$' + val + '$')
      # Set horizontal axis domain. 
      elif key.startswith('xlim'):
        self.ax.set_xlim(val)
      # Set the horizontal axis tick labels. 
      elif key=='xticklabels':
        self.ax.set_xticklabels(val)
      # Vertical axis coordinate. 
      elif key=='y':
        self.y = val
      # Label the vertical axis. 
      elif key=='ylabel':
        self.ax.set_ylabel('' if not val else '$' + val + '$')
      # Set vertical axis domain. 
      elif key.startswith('ylim'):
        self.ax.set_ylim(val)
      # Set the horizontal axis tick labels. 
      elif key=='yticklabels':
        self.ax.set_yticklabels(val)
      # Report any unfamiliar parameters. 
      else:
        print 'WARNING: Unknown param ', key, ' = ', val
    return


  '''

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
'''

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
    lmax = None if self.lines is None else max( l[0][0] for l in self.lines )
    amax = None if self.x is None else np.max(self.x)
    return max(lmax, amax)

  def ymax(self):
    lmax = None if self.lines is None else max( l[0][1] for l in self.lines )
    amax = None if self.y is None else np.max(self.y)
    return max(lmax, amax)

  def zmax(self):
    return None if self.z is None else np.max(self.z)

  # Minima are tricky, since None counts as smaller than any number. 

  def xmin(self):
    lmin = None if self.lines is None else min( l[0][0] for l in self.lines )
    amin = None if self.x is None else np.min(self.x)
    return min(amin, lmin) if None not in (amin, lmin) else max(amin, lmin)

  def ymin(self):
    lmin = None if self.lines is None else min( l[0][1] for l in self.lines )
    amin = None if self.y is None else np.min(self.y)
    return min(amin, lmin) if None not in (amin, lmin) else max(amin, lmin)

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
    # Draw any lines. 
    if self.lines is not None:
      [ self.ax.plot(x, y, *args, **kargs) for x, y, args, kargs in self.lines ]

#    # These subplots can get cramped. Let's reduce the number of ticks. 
#    self.ax.xaxis.set_major_locator( plt.MaxNLocator(3) )
#    self.ax.yaxis.set_major_locator( plt.MaxNLocator(3) )

#    # Only put numbers on axes with labels (usually, just the edge axes). 
#    if not self.ax.get_xlabel():
#      self.ax.set_xticklabels( [] )
#    if not self.ax.get_ylabel():
#      self.ax.set_yticklabels( [] )

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

  def __init__(self, zmax, cax, colorbar=None, ncolors=8):
    # Some plots don't have contours. 
    if not zmax or not cax or not colorbar:
      return dict.__init__(self, {})
    # Store the data scale so that it can be hard-wired into our normalization
    # functions. We don't want to pass vmax all over the place. 
    self.zmax = zmax
    self.colorbar = colorbar
    self.ncolors = ncolors
    self.nticks = ncolors - 1
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
    self.setColorbar(cax, **temp)
    # Become a dictionary of color parameters to be used by the contour plots. 
    return dict.__init__(self, temp)

  # ---------------------------------------------------------------------------
  # ----------------------------------- Tick Locations and Contour Color Levels
  # ---------------------------------------------------------------------------

  def linTicksLevels(self):
    ticks = np.linspace( -self.zmax, self.zmax, self.nticks)
    levels = np.linspace(-self.zmax, self.zmax, self.ncolors)
    # Make sure that the middle tick is exactly zero. 
    ticks[ len(ticks)/2 ] = 0.
    return ticks, levels

  def logTicksLevels(self):
    # One tick at each order of magnitude. 
    power = int( np.floor( np.log10(self.zmax) ) )
    self.zmin = self.zmax/10**self.nticks
    ticks = [ 10**(power - i) for i in range(self.nticks) ]
    logMin, logMax = np.log10(self.zmin), np.log10(self.zmax)
    levels = np.logspace(logMin, logMax, self.ncolors)
    return ticks, levels

  def symTicksLevels(self):
    # A tick at zero, then one per order of magnitude. 
    norders = (self.nticks - 1)/2
    power = int( np.floor( np.log10(self.zmax) ) )
    posTicks = [ 10**(power - i) for i in range(norders) ]
    # For uniform tick spacing, the log cutoff needs to be a factor of ten
    # smaller than the lowest positive tick. 
    self.zmin = min(posTicks)/10.
    ticks = sorted( posTicks + [0] + [ -t for t in posTicks ] )
    # We figure out color levels by spacing them evenly on the unit interval,
    # then mapping the unit interval to the symlog scale. 
    levels = [ self.symNorm(x) for x in np.linspace(0, 1, self.ncolors) ]
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
      ax.set_yticklabels( [ fmt(t) for t in colorParams['ticks'] ] )
      return
    elif self.colorbar=='sym':
      norm, mron, fmt = self.symNorm, self.symMron, self.symFormatter
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
    # If our numbers are around order unity, show floats to two significant
    # digits. Otherwise, use scientific notation. 
    elif 1e-3<self.zmax<1e3:
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

# Convert a string to a complex or float. 
def com(x):
  if ',' not in x:
    return float(x)
  else:
    # Shave off the parentheses then split into real and imaginary parts. 
    re, im = x[1:-1].split(',')
    return (float(re) + float(im)*1j)

# Turn a string into a float or integer. 
def num(x):
  return int( float(x) ) if float(x)==int( float(x) ) else float(x)

# #############################################################################
# #################################################### Array Access and Storage
# #############################################################################

# Data is stored in files full of complex numbers. We want to plot arrays of
# real numbers. 

# =============================================================================
# =============================================================== Array Storage
# =============================================================================

# Given a complex array, this class decides whether the real or imaginary
# component is more interesting (based on its median). It keeps only that
# slice, which now acts like a real array, and remembers which slice it kept. 
class arr(np.ndarray):
  # Proper use of __new__, __init__, and __array_finalize__ is tricky. Luckily,
  # there are plenty of examples online. 
  def __new__(cls, inputArray):
    # If this object is complex...
    if np.iscomplexobj(inputArray):
      re = np.median( np.abs( np.real(inputArray) ) )
      im = np.median( np.abs( np.imag(inputArray) ) )
      # Figure out whether the real or imaginary values are more interesting.
      # Add a label that we can just slap right into a plot title. 
      if im>re:
        obj = np.asarray( np.imag(inputArray) ).view(cls)
        obj.phase = ' \\mathbb{I}\\mathrm{m}\\;\\; '
      else:
        obj = np.asarray( np.real(inputArray) ).view(cls)
        obj.phase = ' \\mathbb{R}\\mathrm{e}\\;\\; '
    # If this isn't a complex array, we don't really do anything. 
    else:
      obj = np.asarray(inputArray).view(cls)
      obj.phase = ''
    return obj
  # Finalize the array...
  def __array_finalize__(self, obj):
    self.phase = getattr(obj, 'phase', None)
    return

# =============================================================================
# =========================================================== Array File Parser
# =============================================================================

# Read a file of values into an array. We expect the first line to give the
# array dimensions. The rest of the file should just be real or complex
# numbers, and is agnostic to whitespace and newlines. 
def readArray(filename):
  # The out prefix is just an older convention for dat files, Fortran output.
  # We ultimately want to use the Python data format, pickles. 
  name = filename[ :filename.rfind('.') ]
  datname, outname, pklname = name + '.dat', name + '.out', name + '.pkl'

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
    # Return the array as an arr (see above class definition). 
    return arr(inputArray)
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

#    # Check if the data dimensions are consistent with the dimensions in the
#    # header of the dat file. 
#    if by(dims)!=by(inputArray.shape):
#      print '\tWARNING: Expected ' + by(dims) + ' but found ' + by(inputArray.shape)

    # Return the array as an arr (see above class definition). 
    return arr(inputArray)
  # If the pickle doesn't exist and there's no Fortran output, return nothing. 
  else:
    print 'WARNING: ' + datname + ' not found. '

# #####################################################################
# ################################################### For Importability
# #####################################################################

if __name__=='__main__':
  main()

