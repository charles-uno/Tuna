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

  # The Tuna Plotter is in charge of data access. 
  TP = tunaPlotter()

#  plotSigma()

#  for kargs in loopover( fdrive=TP.getValues('fdrive') ):
#    plotGround(TP, **kargs)

  for kargs in loopover( model=(1, 2), fdrive=TP.getValues('fdrive') ):
#  for kargs in loopover( model=(1,), fdrive=(0.019,) ):
    plotLayers(TP, **kargs)

#  for kargs in loopover( model=TP.getValues('model') ):
#    plotEnergy(TP, **kargs)

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

def plotEnergy(TP, model):
  # One modenumber per row. One frequency per column, to a maximum of 4. 
  azms = TP.getValues('azm')
  fdrives = TP.getValues('fdrive')[-4:]
  # The Plot Window does the actual plotting. 
  PW = plotWindow( nrows=len(azms), ncols=len(fdrives) )
  # Set title and labels. 
  title = notex('Poloidal (Blue) and Toroidal (Red) Energy: ') + tex(model)
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = [ tex(fdrive) + notex('Current') for fdrive in fdrives ]
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Iterate over the cells. 
  for row, azm in enumerate(azms):
    for col, fdrive in enumerate(fdrives):
      # Find the data for this cell, and plot it. 
      path = TP.getPath(model=model, fdrive=fdrive, azm=azm)
      PW[row, col].setParams( **TP.getCoords(path, 't', 'logU') )
      PW[row, col].setLine( np.log10( TP.getArray(path, 'Upol') ), 'b')
      PW[row, col].setLine( np.log10( TP.getArray(path, 'Utor') ), 'r')
  # Manually clean up the axes. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  PW.setParams( ylims=(2, 6), yticks=(2, 3, 4, 5, 6), yticklabels=('$2$', '', '$4$', '', '$6$') )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'U_' + str(model)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# ======================== Contour Plot of Poloidal and Toroidal Energy Density
# =============================================================================

def plotLayers(TP, model, fdrive):
  # With only two columns, we can't fit 7 rows without deforming them like
  # crazy. Let's keep the frames legible but have fewer of them. 
  azms = (1, 8, 64)
  # The Plot Window does the actual plotting. 
  PW = plotWindow(nrows=len(azms), ncols=2, colorbar='log', zmax=0.1)
  # Set title and labels. 
  title = notex('Energy Density by L-Shell (\\frac{nJ}{m^3}): ') + tex(model) + notex(', ') + tex(fdrive) + notex('Current')
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = ( notex('Poloidal Energy Density'),
                notex('Toroidal Energy Density') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Find the data and set up the axes. 
  for row, azm in enumerate(azms):
    path = TP.getPath(model=model, fdrive=fdrive, azm=azm)
    PW[row, :].setParams( **TP.getCoords(path, 't', 'L0') )
    # Loop over the poloidal and toroidal components. 
    for col, name in enumerate( ('upol', 'utor') ):
      # Grab the energy density and the differential volume. 
      u, dV = TP.getArray(path, name), TP.getArray(path, 'dV')
      dU = u*dV[:, :, None]
      # Compute the mean energy density at each L shell. 
      UofL = np.sum(dU, axis=1)
      VofL = np.sum(dV, axis=1)
      # Careful... dV is 0 at the edges. 
      VofL[0], VofL[-1] = VofL[1], VofL[-2]
      uofL = UofL/VofL[:, None]
      PW[row, col].setContour(uofL)
  # Manually clean up the x axis. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'layers_' + znt(1e3*fdrive) + 'mHz_' + str(model)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# =========== Contour Plot of Active/Quiet North/East Dayside Ground Signatures
# =============================================================================

def plotGround(TP, fdrive):
  azms = TP.getValues('azm')
  # The Plot Window does the actual plotting. 
  PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=10)
  # Set title and labels. 
  title = notex('Magnetic Ground Signatures (nT): ') + tex(fdrive) + notex('Current')
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = ( tex(1) + tex('Bq'), tex(1) + tex('Bf'),
                tex(2) + tex('Bq'), tex(2) + tex('Bf') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Iterate through the rows. Each row is its own modenumber. 
  for row, azm in enumerate(azms):
    # Grab the data for model 1 and plot it. 
    path = TP.getPath(model=1, fdrive=fdrive, azm=azm)
    PW[row, 0].setContour( TP.getArray(path, 'BqE')[:, 0, :] )
    PW[row, 1].setContour( TP.getArray(path, 'BfE')[:, 0, :] )
    PW[row, 0:2].setParams( **TP.getCoords(path, 't', 'lat0') )
    # Grab the data for model 2 and plot it. 
    path = TP.getPath(model=2, fdrive=fdrive, azm=azm)
    PW[row, 2].setContour( TP.getArray(path, 'BqE')[:, 0, :] )
    PW[row, 3].setContour( TP.getArray(path, 'BfE')[:, 0, :] )
    PW[row, 2:4].setParams( **TP.getCoords(path, 't', 'lat0') )
  # Manually clean up the x axis. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'ground_' + znt(1e3*fdrive) + 'mHz'
    return PW.render( TP.savepath + name + '.pdf' )
  else:
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

