#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# Note: This document wraps at column 80. 

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This plotter makes plots. The under-the-hood working are defined in plotmod. 

# #############################################################################
# ##################################################### Import Python Libraries
# #############################################################################

from plotmod import *

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

#  TP = tunaPlotter()

#  for kargs in loopover( fdrive=TP.getValues('fdrive'), mode=('BE', 'PT') ):
#    TP.plotUlines(**kargs)

#  plotSigma()

  for kargs in loopover( fdrive=(0.010, 0.013, 0.016, 0.019, 0.022) ):
    plotGround(**kargs)

  return

# #############################################################################
# ########################################################### Plotting Routines
# #############################################################################

# =============================================================================
# =========================================== Ionospheric Conductivity Profiles
# =============================================================================

def plotSigma():
  # Create the window. 
  PW = plotWindow(nrows=2, ncols=2, colorbar=None)
  # For this plot, we don't actually need the 2D arrays that Tuna spits out. We
  # can just read in the profiles directly. 
  for i in range(4):
    datname = './models/ionpar' + str(i+1) + '.dat'
    with open(datname, 'r') as fileobj:
      lines = fileobj.readlines()
    ionos = np.array( [ [ float(x) for x in l.split() ] for l in lines ] )
    # Log altitude, rather than altitude on a log scale. 
    logalt = np.log10( ionos[:, 0] )
    # Log conductivity instead of conductivity on a log scale. Values are given
    # in mS/m; plot in S/m. 
    logsig0 = np.log10( 1e6/( phys.mu0*ionos[:, 4] ) ) - 3
    logsigH = np.log10( 4*np.pi*phys.eps0*ionos[:, 3] ) - 3
    logsigP = np.log10( 4*np.pi*phys.eps0*ionos[:, 2] ) - 3
    # Add the lines to the plot. 
    PW[i].setLine(logsigP, logalt, 'r')
    PW[i].setLine(logsigH, logalt, 'b')
    PW[i].setLine(logsig0, logalt, 'g')
  # Set the labels and title. 
  colLabels = [ notex('Active'), notex('Quiet') ]
  rowLabels = [ notex('Day'), notex('Night') ]
  title = notex('Pedersen (Blue), Hall (Red), and Parallel (Green) ' + 
                'Conductivities')
  xlabel = notex('Log Conductivity (\\frac{S}{m})')
  xlims = (-15, 5)
  ylabel = notex('Log Altitude (km)')
  ylims = (2, 4)
  PW.setParams(collabels=colLabels, rowlabels=rowLabels, xlabel=xlabel, 
               xlims=xlims, ylabel=ylabel, ylims=ylims, title=title)
  return PW.render()

# =============================================================================
# =================================== Line Plot of Poloidal and Toroidal Energy
# =============================================================================

# Columns are frequency: 13mHz, 16mHz, 19mHz, 22mHz. All current. 
# Columns are modenumber: 1, 2, 4, 8, 16, 32, 64. 
# Fix the model. We'll want a figure for day and one for night. 

# =============================================================================
# ======================== Contour Plot of Poloidal and Toroidal Energy Density
# =============================================================================

# Columns are poloidal and toroidal energy density, averaged over an L-shell. 
# Rows are modenumber: probably 1, 8, 64. Maybe 1, 4, 16, 64. 
# Fix the model. Only model 1 or only model 2. 

# =============================================================================
# =========== Contour Plot of Active/Quiet North/East Dayside Ground Signatures
# =============================================================================

def plotGround(fdrive):
  # The Tuna Plotter is in charge of data access. 
  TP = tunaPlotter()
  azms = TP.getValues('azm')
  # The Plot Window does the actual plotting. 
  PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=10)
  # Set title and labels. 
  title = notex('Magnetic Ground Signatures: ') + tex(fdrive) + notex(' Current')
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = ( tex(1) + tex('Bq'), tex(2) + tex('Bq'),
                tex(1) + tex('Bf'), tex(2) + tex('Bf') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Iterate through the rows. Each row is its own modenumber. 
  for row, azm in enumerate(azms):
    # Grab the data for model 1 and plot it. 
    path = TP.getPath(model=1, fdrive=fdrive, azm=azm)
    PW[row, 0].setContour( TP.getArray(path, 'BqE')[:, 0, :] )
    PW[row, 2].setContour( TP.getArray(path, 'BfE')[:, 0, :] )
    PW[row, 0].setParams( **TP.getCoords(path, 't', 'lat0') )
    PW[row, 2].setParams( **TP.getCoords(path, 't', 'lat0') )
    # Grab the data for model 2 and plot it. 
    path = TP.getPath(model=2, fdrive=fdrive, azm=azm)
    PW[row, 1].setContour( TP.getArray(path, 'BqE')[:, 0, :] )
    PW[row, 3].setContour( TP.getArray(path, 'BfE')[:, 0, :] )
    PW[row, 1].setParams( **TP.getCoords(path, 't', 'lat0') )
    PW[row, 3].setParams( **TP.getCoords(path, 't', 'lat0') )
  # Manually clean up the x axis. 
  PW.setParams( xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )

  return PW.render()







'''

      # Try to de-cramp cramped plots a bit. 
      if ymin==2 and ymax==6:
        for cell in self.cells.flatten():
          cell.setParams( yticks=(2, 3, 4, 5, 6), yticklabels=('$2$', '', '$4$', '', '$6$') )

      if ymin==2 and ymax==10:
        for cell in self.cells.flatten():
          cell.setParams( yticks=(2, 4, 6, 8, 10), yticklabels=('$2$', '', '$6$', '', '$10$') )

      # Try to de-cramp cramped plots a bit. 
      if xmin==1 and xmax==300:
        for cell in self.cells.flatten():
          cell.setParams( xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )

      if all( cell.xlims==(-15, 5) for cell in self.cells.flatten() ):
        for cell in self.cells.flatten():
          cell.setParams( xticks=(-15, -10, -5, 0, 5), xticklabels=('$-15$', '', '$-5$', '', '$5$') )

    if self.xlims==(-3, 1):
      self.setParams( xticks=(-3, -2, -1, 0, 1), xticklabels=('$-3$', '$-2$', '$-1$', '$0$', '$1$') )

'''

# #####################################################################
# ################################################### For Importability
# #####################################################################

if __name__=='__main__':
  main()

