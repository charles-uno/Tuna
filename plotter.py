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

  # Plots for the defense...

#  waveProperties(TP)

#  for kargs in loopover( fdrive=(0.010, 0.016, 0.022) ):
#    driveComparison(TP, **kargs)

#  va(TP)

#  fa(TP)

#  path = TP.getPath(azm=32, fdrive=0.022, model=2)
#  plotFields(TP, path=path)

#  energy(TP, model=2, fdrive=0.022)
#  energy(TP, model=4, fdrive=0.013)

  ground(TP, side='night', fdrive=0.016)


#  # Let's just tally up Bz/Bx as a function of m. 
##  for kargs in loopover( fdrive=(0.007, 0.010, 0.013, 0.016, 0.019, 0.022, 0.025) ):
#  for kargs in loopover( model=(1, 2, 3, 4), ldrive=(5,) ):
#    plotCompression(TP, **kargs)

#  # Look at snapshots of a run. 
#  for kargs in loopover( path=TP.paths.keys() ):
#    plotFields(TP, **kargs)

#  # Plot the radial distribution in energy. 
#  for kargs in loopover( mode=('p', 't'), model=(3, 4), ldrive=(5,) ):
#    plotLayers(TP, **kargs)

#  for kargs in loopover( model=(3, 4), ldrive=(5,) ):
#    plotEnergy(TP, **kargs)

#  # Plot magnetic field signatures at the ground. 
#  for kargs in loopover( fdrive=(0.007, 0.010, 0.013, 0.016, 0.019, 0.022, 0.025), side=('day', 'night',), ldrive=(5,) ):
#    plotGround(TP, **kargs)

#  # Do poloidal and toroidal waves have opposite phase offsets?
#  for kargs in loopover( path=TP.paths.keys() ):
#    plotPhase(TP, **kargs)

#  # Schematic illustrating poloidal and toroidal waves. 
#  plotToroidal(TP)
#  plotPoloidal(TP)
#  plotPoloidalToroidal(TP)

#  # Schematic illustrating first and second harmonics. 
#  plotOddEven(TP)

#  # Plasma frequency profile. 
#  plotPlasmaFrequency(TP)

#  # Alfven speed profiles. 
#  plotAlfvenSpeed(TP)

#  # Compressional Alfven frequency cutoff. 
#  plotAlfvenCutoff(TP, model=2)

#  # Snapshots of the electric field components. 
#  plotSnapshots(TP, model=2, fdrive=0.016, azm=16)

#  for kargs in loopover( model=(1, 2), fdrive=(0.016, ), alt=(100, 1000) ):
#    plotSlice(TP, **kargs)

#  # Comparison of J dot E to the divergence of the Poynting flux. 
#  plotDivFlux(TP, model=2, fdrive=0.016)

#  # Zoom way in on runs with and without inertial length scales. 
#  plotZoom(TP, azm=16)

#  # Conductivity profiles. 
#  plotSigma(TP)

#  # The grid. 
#  plotGrid(TP)

#  # Plot the Alfven bounce profiles.
#  plotBounceFrequency(TP)

#  # Sym-H index and its Fourier transform. 
#  plotSymh(TP)

#  # Show waves failing to propagate in when driven from the outer boundary. 
#  plotBdrive(TP, fdrive=0.022)

  return

# #############################################################################
# ########################################################### Plotting Routines
# #############################################################################

# =============================================================================
# ====================================== Harmonic, Polarization, and Modenumber
# =============================================================================

def waveProperties(TP):
  PW = plotWindow(nrows=2, ncols=3, colorbar=None, landscape=True, square=True)
  clabs = ( notex('Harmonic'), notex('Polarization (Side View)'), notex('Modenumber (Top View)') )
  PW.setParams(collabels=clabs)
  lims = (-8, 8)
  PW.setParams(xlims=lims, xticks=lims, xticklabels=(), 
               ylims=lims, yticks=lims, yticklabels=() )


  PW[:, 1].setParams(earth='left')
  PW[:, 2].setParams(earth='top')

  PW[0, 0].setParams( toptext=notex('First / Odd') )
  PW[1, 0].setParams( toptext=notex('Second / Even') )

  PW[0, 1].setParams( toptext=notex('Poloidal') )
  PW[1, 1].setParams( toptext=notex('Toroidal') )

  PW[0, 2].setParams( toptext='m = 4' )
  PW[1, 2].setParams( toptext='m = 16' )

  # Draw first and second harmonics. 
  x = np.linspace(lims[0], lims[1], 1000)
  dx = lims[1] - lims[0]

  e1, e2 = np.cos(np.pi*x/dx), np.sin(2*np.pi*x/dx)
  [ PW[0, 0].setLine(x,  pm*e1, 'r' + pat) for pm, pat in ( (2, ''), (-2, ':') ) ]
  [ PW[1, 0].setLine(x,  pm*e2, 'b' + pat) for pm, pat in ( (2, ':'), (-2, '') ) ]

  # Draw poloidal polarization. 
  L, dL = 6, 0.5
  q0 = np.arcsin( np.sqrt(1./L) )
  q = np.linspace(q0, np.pi - q0, 100)
  # Normalized coordinate for convenience. 
  qp = np.pi*( q - q[0] )/( q[-1] - q[0] )
  r = L*np.sin(q)**2
  dr = dL*np.sin(qp)
  [ PW[0, 1].setLine( -(r + pm*dr)*np.sin(q), (r + pm*dr)*np.cos(q), 'r' + pat) for pm, pat in ( (1, ''), (-1, ':') ) ]
  dr = dL*np.sin(2*qp)
  [ PW[0, 1].setLine( (r + pm*dr)*np.sin(q), (r + pm*dr)*np.cos(q), 'b' + pat) for pm, pat in ( (1, ''), (-1, ':') ) ]

  # Toroidal polarization. 
  ax = PW.cells[1, 1].ax
  x, z = r*np.sin(q), r*np.cos(q)
  for n, color in ( (1, 'r'), (2, 'b') ):
    w = 5*np.sin(n*qp)
    w = np.where( np.abs(w) < 1, np.sign(w), w )
    for i in range( 1, qp.size, 2 ):
      pat = '' if w[i] > 0 else ':'
      ax.plot( (-1)**n * x[i-1:i+2], z[i-1:i+2], color + pat, linewidth=np.abs( w[i] ) )

  # Azimuthal modenumber. 
  q = np.linspace(0, 2*np.pi, 1000)
  r = 4

  for row, azm in enumerate( (4, 16) ):
    dr = np.sin(azm*q)
    PW[row, 2].setLine( (r + dr)*np.sin(q), (r + dr)*np.cos(q), 'k' )
    PW[row, 2].setLine( (r + dr)*np.cos(q), (r + dr)*np.sin(q), 'k:' )


  if TP.savepath is not None:
    return PW.render(TP.savepath + 'properties.pdf')
  else:
    return PW.render()

# =============================================================================
# ======================================== Driving with Compression and Current
# =============================================================================

def driveComparison(TP, fdrive):
  azms = (1, 4, 16, 64)
  model = 2
  PW = plotWindow(nrows=2, ncols=len(azms), colorbar='log', landscape=True)
  # Set title and labels. 
  title = notex('Mean Energy Density: ') + tex(fdrive) + notex('Driving')
  ulab = notex('\\frac{nJ}{m^3}')
  clabs = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  rlabs = [ notex('Compression\nat L=10'), notex('Current\nat L\\sim5') ]
  PW.setParams(title=title, collabels=clabs, rowlabels=rlabs, unitlabel=ulab, zmax=1)

  # Look at the runs driven with magnetic compression.  
  TP.setPaths('/media/My Passport/RUNS/BDRIVE/')
  for col, azm in enumerate(azms):
    path = TP.getPath(model=model, fdrive=fdrive, azm=azm)

    if path is None:
      PW[0, col].setParams( text=notex('Not Found') )
      continue

    PW[0, col].setParams( **TP.getCoords(path) )
    u = np.mean(TP.getArray(path, 'u'), axis=-1)
    PW[0, col].setContour(u)

  # Look at the runs driven with ring current modulation.  
  TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_4/')
  for col, azm in enumerate(azms):
    path = TP.getPath(model=model, fdrive=fdrive, azm=azm)

    if path is None:
      PW[1, col].setParams( text=notex('Not Found') )
      continue

    PW[1, col].setParams( **TP.getCoords(path) )
    u = np.mean(TP.getArray(path, 'u'), axis=-1)
    PW[1, col].setContour(u)

  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'drivers_' + znt(1e3*fdrive) + 'mHz'
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# =================================================== Energy Lines and Contours
# =============================================================================

def energy(TP, model, fdrive):

  azms = (1, 2, 4, 8, 16, 32, 64)

  clabs = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  rlabs = ( 'U_P' + notex(' (Blue)\nvs\n') + 'U_T' + notex(' (Red)'), 
            notex('Poloidal\nEnergy\nDensity'), 
            notex('Toroidal\nEnergy\nDensity') )
  title = notex('Poloidal and Toroidal Energy: ') + tex(model) + notex(', ') + tex(fdrive) + notex('Current')


  PW = plotWindow(nrows=3, ncols=len(azms), colorbar='log', landscape=True, zmax=0.1)

  PW.setParams(collabels=clabs, unitlabel=tex('mW/m^2'), rowlabels=rlabs, title=title)

  # Iterate across seven runs, one at each modenumber. 
  for col, azm in enumerate(azms):
    path = TP.getPath(azm=azm, fdrive=fdrive, model=model)

    # Top row: poloidal and toroidal energy curves. 
    PW[0, col].setParams( **TP.getCoords(path, 't', 'logU') )
    PW[0, col].setLine( np.log10( TP.getArray(path, 'Upol') ), 'b')
    PW[0, col].setLine( np.log10( TP.getArray(path, 'Utor') ), 'r')

    # Second and third rows: poloidal and toroidal energy contours. 
    PW[1:3, col].setParams( **TP.getCoords(path, 't', 'L0') )

    # Grab the energy density. 
    upt = [ TP.getArray(path, name) for name in ('upol', 'utor') ]
    dV = TP.getArray(path, 'dV')
    # Compute differential energy. 
    dUpt = [ u*dV[:, :, None] for u in upt ]
    # Mean energy density in each L-shell. 
    UofLpt = [ np.sum(dU, axis=1) for dU in dUpt ]
    VofL = np.sum(dV, axis=1)
    # Careful... dV is 0 at the edges. 
    VofL[0], VofL[-1] = VofL[1], VofL[-2]

    uofLpt = [ UofL/VofL[:, None] for UofL in UofLpt ]
    [ PW[i+1, col].setContour(uofL) for i, uofL in enumerate(uofLpt) ]

    # Manually clean up the y axes. 
#    PW[0, col].setParams( ylims=(2, 6), yticks=(2, 3, 4, 5, 6), yticklabels=('$2$', '', '$4$', '', '$6$') )
#    PW[1:3, col].setParams( ylims=(2, 10), yticks=(2, 4, 6, 8, 10), yticklabels=('$2$', '', '$6$', '', '$10$') )
    PW[0, col].setParams( ylims=(2, 6), yticks=(2, 3, 4, 5, 6), yticklabels=(), ylabelpad=-2, ylabel=notex('Log') + 'U' )
    PW[1:3, col].setParams( ylims=(2, 10), yticks=(2, 4, 6, 8, 10), yticklabels=(), ylabelpad=-2 )

  # Manually clean up the x axis. 
#  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=() )

  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'energy_' + str(model) + '_' + znt(1e3*fdrive)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()




# =============================================================================
# =================================== Line Plot of Poloidal and Toroidal Energy
# =============================================================================

def plotEnergy(TP, model, lpp=4, ldrive=5):

  if lpp==4 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_4/')
  elif lpp==5 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_5/')
  elif lpp==4 and ldrive==6:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LDRIVE_6/')
  elif lpp==4 and ldrive==7:
    TP.setPaths('/media/My Passport/RUNS/LDRIVE_7/')
  else:
    return

  # One modenumber per row. One frequency per column, to a maximum of 4. 
  azms = (1, 2, 4, 8, 16, 32, 64)
  fdrives=(0.013, 0.016, 0.019, 0.022)

  # The Plot Window does the actual plotting. 
  PW = plotWindow( nrows=len(azms), ncols=len(fdrives) )
  # Set title and labels. 


#  title = notex('Poloidal (Blue) and Toroidal (Red) Energy: ') + tex(model) + notex(', ') + 'L_{PP} = ' + str(lpp) + notex(', ') + 'L_{drive} = ' + znt(ldrive)
  title = notex('Poloidal (Blue) and Toroidal (Red) Energy: ') + tex(model)

  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = [ tex(fdrive) + notex('Current') for fdrive in fdrives ]
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Iterate over the cells. 
  for row, azm in enumerate(azms):
    for col, fdrive in enumerate(fdrives):
      # Find the data for this cell, and plot it. 
      path = TP.getPath(model=model, fdrive=fdrive, azm=azm, lpp=lpp)
      # Make sure this run exists. 
      if path is None:
        PW[row, col].setParams( text=notex('Not Found') )
        continue
      PW[row, col].setParams( **TP.getCoords(path, 't', 'logU') )
      # If there was a crash, don't show the data. Just say "CRASH."
      if TP.getArray(path, 't').size < 300:
        PW[row, col].setParams(text=notex('Unstable'))
        continue
      PW[row, col].setLine( np.log10( TP.getArray(path, 'Upol') ), 'b')
      PW[row, col].setLine( np.log10( TP.getArray(path, 'Utor') ), 'r')
  # Manually clean up the axes. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  PW.setParams( ylims=(2, 6), yticks=(2, 3, 4, 5, 6), yticklabels=('$2$', '', '$4$', '', '$6$') )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'U_' + str(model) + '_' + znt(lpp) + '_' + znt(ldrive)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# ======================== Contour Plot of Poloidal and Toroidal Energy Density
# =============================================================================

def plotLayers(TP, mode, model, lpp=4, ldrive=5):

  if lpp==4 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_4/')
  elif lpp==5 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_5/')
  elif lpp==4 and ldrive==6:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LDRIVE_6/')
  elif lpp==4 and ldrive==7:
    TP.setPaths('/media/My Passport/RUNS/LDRIVE_7/')
  else:
    return

  # With only two columns, we can't fit 7 rows without deforming them like
  # crazy. Let's keep the frames legible but have fewer of them. 
  azms = (1, 2, 4, 8, 16, 32, 64)

  fdrives=(0.013, 0.016, 0.019, 0.022)

  # The Plot Window does the actual plotting. 
  PW = plotWindow(nrows=len(azms), ncols=len(fdrives), colorbar='log', zmax=0.1)
  # Set title and labels. 
  modename = {'p':'Poloidal ', 't':'Toroidal '}[mode]
#  title = notex(modename + 'Energy Density by L-Shell: ') + tex(model) + notex(', ') + 'L_{PP} = ' + str(lpp) + notex(', ') + 'L_{drive} = ' + znt(ldrive)
  title = notex(modename + 'Energy Density by L-Shell: ') + tex(model)
  unitlabel = notex('\\frac{nJ}{m^3}')
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = [ tex(fdrive) + notex('Current') for fdrive in fdrives ]
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels, unitlabel=unitlabel)
  # Find the data and set up the axes. 
  for row, azm in enumerate(azms):
    for col, fdrive in enumerate(fdrives):
      path = TP.getPath(model=model, fdrive=fdrive, azm=azm, lpp=lpp)
      # Make sure this run exists. 
      if path is None:
        PW[row, col].setParams( text=notex('Not Found') )
        continue
      PW[row, col].setParams( **TP.getCoords(path, 't', 'L0') )
      # If there was a crash, don't show the data. Just say "CRASH."
      if TP.getArray(path, 't').size < 300:
        PW[row, col].setParams(text=notex('Unstable'))
        continue
      # Grab the energy density and the differential volume. 
      uname = {'p':'upol', 't':'utor'}[mode]
      u, dV = TP.getArray(path, uname), TP.getArray(path, 'dV')
      dU = u*dV[:, :, None]
      # Compute the mean energy density at each L shell. 
      UofL = np.sum(dU, axis=1)
      VofL = np.sum(dV, axis=1)
      # Careful... dV is 0 at the edges. 
      VofL[0], VofL[-1] = VofL[1], VofL[-2]
      uofL = UofL/VofL[:, None]
      PW[row, col].setContour(uofL)

  # Manually clean up the x axis. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$'),
                ylims=(2, 10), yticks=(2, 4, 6, 8, 10), yticklabels=('$2$', '', '$6$', '', '$10$') )

  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'layers_' + mode + '_' + str(model) + '_' + znt(lpp) + '_' + znt(ldrive)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# =========== Contour Plot of Active/Quiet North/East Dayside Ground Signatures
# =============================================================================

def plotGround(TP, fdrive, side='day', lpp=4, ldrive=5):

  if lpp==4 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_4/')
  elif lpp==4 and ldrive==6 and side=='night':
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LDRIVE_6/')
  elif lpp==4 and ldrive==7 and side=='night':
    TP.setPaths('/media/My Passport/RUNS/LDRIVE_7/')
  else:
    return

  azms = (1, 2, 4, 8, 16, 32, 64)
  # The Plot Window does the actual plotting. 
  PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=100)

  # Active model, quiet model. 
  am, qm = (1, 2) if side=='day' else (3, 4)

  # Set title and labels. 
#  title = notex('Magnetic Ground Signatures: ') + tex(fdrive) + notex('Current') + notex(', ') + 'L_{PP} = ' + str(lpp) + notex(', ') + 'L_{drive} = ' + znt(ldrive)

  sidename = {'day':'Dayside', 'night':'Nightside'}[side]
  title = notex(sidename + ' Magnetic Ground Signatures: ') + tex(fdrive) + notex('Current')

  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = ( tex(am) + tex('Bq'), tex(am) + tex('Bf'),
                tex(qm) + tex('Bq'), tex(qm) + tex('Bf') )
  unitlabel = notex('nT')
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels, unitlabel=unitlabel)
  # Iterate through the rows. Each row is its own modenumber. 
  for row, azm in enumerate(azms):
    for c0, model in enumerate( (am, qm) ):
      path = TP.getPath(model=model, fdrive=fdrive, azm=azm, lpp=lpp)
      for dc, name in enumerate( ('BqE', 'BfE') ):
        col = 2*c0 + dc
        # Make sure this run exists. 
        if path is None:
          PW[row, col].setParams( text=notex('Not Found') )
          continue
        PW[row, col].setParams( **TP.getCoords(path, 't', 'lat0') )
        # If there was a crash, don't show the data. Just say "CRASH."
        if TP.getArray(path, 't').size < 300:
          PW[row, col].setParams(text=notex('Unstable'))
          continue
        B = TP.getArray(path, name)[:, 0, :]
        PW[row, col].setContour(B)
        Bmax = np.max( np.abs( B[5:-5, 5:-5] ) )
        if Bmax > np.sqrt(10.):
          PW[row, col].setParams( bottext=notex('Max: ') + format(Bmax, '.0f') + notex('nT') )
  # Manually clean up the x axis. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'ground_' + znt(1e3*fdrive) + 'mHz_' + side + '_' + znt(lpp) + '_' + znt(ldrive)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()







def ground(TP, fdrive, side='day'):

  azms = (1, 2, 4, 8, 16, 32, 64)
  PW = plotWindow(nrows=4, ncols=len(azms), colorbar='sym', zmax=10, landscape=True)

  # Active model, quiet model. 
  am, qm = (1, 2) if side=='day' else (3, 4)

  sidename = {'day':'Dayside', 'night':'Nightside'}[side]
  title = notex(sidename + ' Magnetic Ground Signatures: ') + tex(fdrive) + notex('Current')

  clabs = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
#  rlabs = ( tex(am) + tex('Bq'),
#            tex(am) + tex('Bf'),
#            tex(qm) + tex('Bq'),
#            tex(qm) + tex('Bf') )

  rlabs = ( notex('Active') + tex('Bq'),
            notex('Active') + tex('Bf'),
            notex('Quiet') + tex('Bq'),
            notex('Quiet') + tex('Bf') )

  ulab = notex('nT')
  PW.setParams(title=title, collabels=clabs, rowlabels=rlabs, unitlabel=ulab)
  # Iterate through the rows. Each row is its own modenumber. 
  for col, azm in enumerate(azms):
    for r0, model in enumerate( (am, qm) ):
      path = TP.getPath(model=model, fdrive=fdrive, azm=azm)
      for dr, name in enumerate( ('BqE', 'BfE') ):
        row = 2*r0 + dr
        # Make sure this run exists. 
        if path is None:
          PW[row, col].setParams( text=notex('Not Found') )
          continue
        PW[row, col].setParams( **TP.getCoords(path, 't', 'lat0') )
        # If there was a crash, don't show the data. Just say "CRASH."
        if TP.getArray(path, 't').size < 300:
          PW[row, col].setParams(text=notex('Unstable'))
          continue
        B = TP.getArray(path, name)[:, 0, :]
        PW[row, col].setContour(B)
        Bmax = np.max( np.abs( B[5:-5, 5:-5] ) )
        if Bmax > np.sqrt(10.):
          PW[row, col].setParams( bottext=notex('Max: ') + format(Bmax, '.0f') + notex('nT') )
  # Manually clean up the x axis. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=(), ylabel='' )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'ground_' + znt(1e3*fdrive) + 'mHz_' + side
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# ================================= Relation Between Modenumber and Compression
# =============================================================================

def plotCompression(TP, model=None, fdrive=None, lpp=4, ldrive=5):
  # One modenumber per row. Columns may be frequencies or models. 
  azms = (1, 2, 4, 8, 16, 32, 64)


  if lpp==4 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_4/')
  elif lpp==5 and ldrive==5:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LPP_5/')
  elif lpp==4 and ldrive==6:
    TP.setPaths('/media/My Passport/RUNS/JDRIVE_LDRIVE_6/')
  elif lpp==4 and ldrive==7:
    TP.setPaths('/media/My Passport/RUNS/LDRIVE_7/')
  else:
    return

  if model is None and fdrive is None:
    print 'What should the columns be? '
    return
  elif model is not None:
    items0 = [ ('model', model) ]
    key = 'fdrive'
    vals = (0.013, 0.016, 0.019, 0.022)
    clabs = [ tex(val) + notex('Current') for val in vals ]
    title = notex('Compressional Coupling to the Poloidal Mode: ') + tex(model)
  elif fdrive is not None:
    items0 = [ ('fdrive', fdrive) ]
    key = 'model'
    vals = (1, 2, 3, 4)
    clabs = [ tex(val) for val in vals ]
    title = notex('Compressional Coupling to the Poloidal Mode: ') + tex(fdrive) + notex('Current')
  PW = plotWindow( nrows=len(azms), ncols=len(vals) )
  # Set title and labels. 
  rlabs = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  PW.setParams(title=title, collabels=clabs, rowlabels=rlabs)

  # Iterate over the cells. 
  for row, azm in enumerate(azms):
    for col, val in enumerate(vals):

      # Find the data for this cell, and plot it. 
      path = TP.getPath(azm=azm, lpp=lpp, ldrive=ldrive, **dict( items0 + [ (key, val) ] ) )
      # Make sure this run exists. 
      if path is None:
        PW[row, col].setParams( text=notex('Not Found') )
        continue

      PW[row, col].setParams( **TP.getCoords(path, 't', 'comp') )
      # If there was a crash, don't show the data. Just say "CRASH."
      if TP.getArray(path, 't').size < 300:
        PW[row, col].setParams(text=notex('Unstable'))
        continue

      RMS = np.sqrt( TP.getArray(path, 'UBz') / TP.getArray(path, 'UBx') )

      PW[row, col].setLine(RMS, 'r')

      dots = np.mean(RMS)*np.ones( RMS.shape )

      PW[row, col].setLine(dots, 'k:')

      PW[row, col].setParams(toptext=tdp(np.mean(RMS)*100.) + '\\%')


  print 'LPP, LDRIVE = ', lpp, ldrive


  # Manually clean up the axes. 
  PW.setParams( xlims=(0, 300), xticks=(0, 100, 200, 300), xticklabels=('$0$', '', '', '$300$') )
  PW.setParams( ylims=(0, 1.5), yticks=(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5), yticklabels=('$0\\%$', '', '', '$75\\%$', '', '', '$150\\%$'), ylabelpad=-3 )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'comp_' + ( str(model) if model is not None else znt(1000*fdrive) + 'mHz' )
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()







# =============================================================================
# ====================================================== Snapshot of All Fields
# =============================================================================

def plotFields(TP, path):
  # What are we looking at? 
  steps = (269, 279, 289, 299)
  names = ('Bx', 'By', 'Bz')
  # Set up the plot window. 
  PW = plotWindow(nrows=len(steps), ncols=len(names), colorbar='sym', zmax=10)
  # Titles and labels. 
  fdrive, model, azm = [ TP.getParams(path)[key] for key in ('fdrive', 'model', 'azm') ]

  title = notex('Magnetic Field Snapshots: ') + tex(model) + notex(', ') + tex(fdrive) + notex('Current, ') +  'm = ' + znt(azm)

  rlabs = [ notex(znt(step + 1) + 's') for step in steps ]


#  clabs = [ tex(name) for name in names ]
  clabs = ( 'B_x' + notex(' (Poloidal)'), 
            'B_y' + notex(' (Toroidal)'), 
            'B_z' + notex(' (Compressional)') )

  PW.setParams(title=title, rowlabels=rlabs, collabels=clabs, unitlabel=notex('nT'))
  # Set up the axes and labels. 
  PW.setParams( **TP.getCoords(path) )
  # Add the data to the plot. 
  for col, name in enumerate(names):
    z = TP.getArray(path, name)
    for row, step in enumerate(steps):
      PW[row, col].setContour( z[:, :, step] )
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'snapshot_' + str(model) + '_' + znt(1000*fdrive) + 'mHz_' + znt(azm, 3)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()


  return PW.render()

# =============================================================================
# ============================================ Finagling with Phase... Not Done
# =============================================================================

def dft(t, x):
  n = t.size
  def omega(m):
    return 2*m*np.pi/n
  def harm(m):
    return np.exp(1j*omega(m)*t)
  weights = np.zeros(n/2, dtype=np.complex)
  for m in range(n/2):
    weights[m] = np.sum( x*harm(-m) ) / np.sum( harm(m)*harm(-m) )
  return np.array( [ omega(m) for m in range(n/2) ] ), weights

def tfd(omega, weights):
  t = np.arange(2*omega.size)
  x = np.zeros(2*omega.size, dtype=np.complex)

  m = np.argmax( np.abs(weights) )

  return t, weights[m]*np.exp(1j*omega[m]*t)

#  for m in range(omega.size)[:10]:
#    x = x + weights[m]*np.exp(1j*omega[m]*t)
#  return t, x


def plotPhase(TP, path):
  PW = plotWindow(nrows=2, ncols=1, colorbar=None)
  fdrive, model, azm = [ TP.getParams(path)[key] for key in ('fdrive', 'model', 'azm') ]

  Ex, Ey, Bx, By = [ TP.getArray(path, name) for name in ('Ey', 'Bx', 'Ex', 'By') ]

  t = TP.getArray(path, 't')
  SP = Ey*Bx
  ST = Ex*By
  n1, n3, nt = SP.shape
  i = np.random.randint(0, n1)
  k = np.random.randint(0, n3/2)

  E = Ey[i, k, :]
  B = Bx[i, k, :]
  PW[0].setLine(t, E, 'b')
  PW[0].setLine(t, B, 'r')

  omega, Ep = dft(t, E)
  omega, Bp = dft(t, B) 

  Sp = Ep*np.conj(Bp)

  m = np.argmax(Sp)

  tp, Ep = tfd( *dft(t, E) )
  PW[0].setLine(tp, Ep, 'b+')
  tp, Bp = tfd( *dft(t, B) )
  PW[0].setLine(tp, Bp, 'r+')
  E = Ex[i, k, :]
  B = By[i, k, :]
  PW[1].setLine(t, E, 'b')
  PW[1].setLine(t, B, 'r')
  tp, Ep = tfd( *dft(t, E) )
  PW[1].setLine(tp, Ep, 'b+')
  tp, Bp = tfd( *dft(t, B) )
  PW[1].setLine(tp, Bp, 'r+')
  PW.setParams( ylims=(-10, 10) )

  X, Z = [ TP.getArray(path, name) for name in ('X', 'Z') ]

  title = tex(fdrive) + '\\qquad m = ' + znt(azm) + '\\qquad ' + tex(model) + '\\qquad{}X = ' + format(X[i, k], '.1f') + '\\qquad{}Z = ' + format(Z[i, k], '.1f')
  rlabs = [ notex(modename) for modename in ('Poloidal', 'Toroidal') ]
  PW.setParams(title=title, rowlabels=rlabs)
  return PW.render()

# =============================================================================
# ================================ Parallel Current and Poynting Flux Snapshots
# =============================================================================

def altslice(arr, r, alt):
  # Altitude in Mm. 
  a = 1e3*(r - phys.RE)
  # Make a new array, skipping the coordinate that indexes along field lines. 
  n1, n3, nt = arr.shape
  ret = np.zeros( (n1, nt) )
  # For each field line, find the k value that is closest to alt. 
  for i in range(n1):
    k = np.argmin( np.abs(a[i, :n3/2] - alt) )
    ret[i, :] = arr[i, k, :]
  return ret

def plotSlice(TP, model, fdrive, alt=100):

  TP.setPaths('/media/My Passport/RUNS/INERTIA/')

  # Set up the window. 
  azms = (1, 4, 16, 64)
  PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=1)
  # Set title and labels. 
  title = ( notex('Current and Poynting Flux at ' + znt(alt) + 'km: ') +
            tex(model) + notex(', ') + tex(fdrive) + notex(' Current') )

  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]

  collabels = ( tex('real') + notex(' ') + 'J_z' + notex(' (\\frac{\\mu{}A}{m^2})'),
                'S_P' + notex(' (\\frac{mW}{m^2})'),
                tex('imag') + notex(' ') + 'J_z' + notex(' (\\frac{\\mu{}A}{m^2})'),
                'S_T' + notex(' (\\frac{mW}{m^2})') )

#  collabels = ( tex('real') + 'J_z' + notex(' (\\frac{\\mu{}A}{m^2})'),
#                tex('real') + 'S_z' + notex(' (\\frac{mW}{m^2})'),
#                tex('imag') + 'J_z' + notex(' (\\frac{\\mu{}A}{m^2})'),
#                tex('imag') + 'S_z' + notex(' (\\frac{mW}{m^2})') )

  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Each row comes from a different run. 
  for row, azm in enumerate(azms):
    path = TP.getPath(model=model, fdrive=fdrive, azm=azm)

    if path is None:
      print 'No such run. '
      return False

    # Take slices of the field at the boundary and a ways up. 
    PW[row, :].setParams( **TP.getCoords(path, 't', 'lat0') )
    r = TP.getArray(path, 'r')

#    RJ = altslice(TP.getArray(path, 'Jz', real=True), r, alt)
#    RS = altslice(TP.getArray(path, 'Sz', real=True), r, alt)
#    IJ = altslice(TP.getArray(path, 'Jz', real=False), r, alt)
#    IS = altslice(TP.getArray(path, 'Sz', real=False), r, alt)

#    # Add them to the plot. 
#    PW[row, 0].setContour(RJ)
#    PW[row, 1].setContour(RS)
#    PW[row, 2].setContour(IJ)
#    PW[row, 3].setContour(IS)

    RJ = altslice(TP.getArray(path, 'Jz', real=True), r, alt)
    SP = altslice(TP.getArray(path, 'Spol'), r, alt)
    IJ = altslice(TP.getArray(path, 'Jz', real=False), r, alt)
    ST = altslice(TP.getArray(path, 'Stor'), r, alt)
    # Add them to the plot. 
    PW[row, 0].setContour(RJ)
    PW[row, 1].setContour(SP)
    PW[row, 2].setContour(IJ)
    PW[row, 3].setContour(ST)

  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'slice_' + znt(1e3*fdrive) + 'mHz_' + str(model) + '_' + znt(alt) + 'km'
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

def plotCurrentFlux(TP, model, fdrive, alt=100):
  TP.setPaths('/media/My Passport/RUNS/INERTIA/')
  # Set up the window. 
  azms = (1, 4, 16, 64)
  PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=1)
  # Set title and labels. 
  title = ( notex('Current and Poynting Flux at ') + 'R_I' + notex(': ') +
            tex(model) + notex(', ') + tex(fdrive) + notex(' Current') )
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = ( tex('real') + notex(' ') + 'J_z' + notex(' (\\frac{\\mu{}A}{m^2})'),
                'S_P' + notex(' (\\frac{mW}{m^2})'),
                tex('imag') + notex(' ') + 'J_z' + notex(' (\\frac{\\mu{}A}{m^2})'),
                'S_T' + notex(' (\\frac{mW}{m^2})') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
  # Each row comes from a different run. 
  for row, azm in enumerate(azms):
    path = TP.getPath(model=model, fdrive=fdrive, azm=azm)
    PW[row, :].setParams( **TP.getCoords(path, 't', 'lat0') )
    # Load up the Poynting flux and the real and imaginary current. 
    RJ = TP.getArray(path, 'Jz', real=True)[:, 0, :]
    SP = TP.getArray(path, 'Spol')[:, 0, :]
    IJ = TP.getArray(path, 'Jz', real=False)[:, 0, :]
    ST = TP.getArray(path, 'Stor')[:, 0, :]
    # Add them to the plot. 
    PW[row, 0].setContour(RJ)
    PW[row, 1].setContour(SP)
    PW[row, 2].setContour(IJ)
    PW[row, 3].setContour(ST)
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'JS_' + znt(1e3*fdrive) + 'mHz_' + str(model)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# ==================================================== Electric Field Snapshots
# =============================================================================

def plotSnapshots(TP, model, fdrive, azm):
  TP.setPaths('/media/My Passport/RUNS/INERTIA/')
  # Each row is a snapshot at a different time. Three rows is probably enough. 
  steps = (259, 279, 299)
  # Plot is handled using a Plot Window object. 
  zmax = 10
  PW = plotWindow(ncols=3, nrows=len(steps), colorbar='sym', zmax=zmax)
  # Grab the path for the run we're looking at. Then grab the field data. 
  path = TP.getPath(model=model, fdrive=fdrive, azm=azm)
  Ex, Ey, Ez = [ TP.getArray(path, name) for name in ('Ex', 'Ey', 'Ez') ]
  # Figure out an appropriate factor by which to scale up Ez. 
  scalefac = np.floor( np.log10( zmax*np.sqrt(10)/np.max(Ez) ) )

  scalefac = 6

  # Set up the plot axes. 
  PW.setParams( **TP.getCoords(path) )
  # Add each contour to the plot. 
  for row, step in enumerate(steps):
    PW[row, 0].setContour( Ex[:, :, step] )
    PW[row, 1].setContour( Ey[:, :, step] )
    PW[row, 2].setContour( Ez[:, :, step]*10**scalefac )
  # Set the plot title and labels. 
  collabels = ( tex('Ex') + notex(' (Toroidal)'), tex('Ey') + notex(' (Poloidal)'), 
               '10^{' + znt(scalefac) + '} \\times ' + tex('Ez') + notex(' (Parallel)') )
  rowlabels = [ notex(str(s+1) + 's') for s in steps ]
  title = ( notex('Electric Field Snapshots: ' + tex(model) + ', ' +
                  tex(fdrive) + 'Current, ') + 'm = ' + znt(azm) )
  unitlabel = notex('\\frac{mV}{m}')
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels, 
               unitlabel=unitlabel)
  # Name this figure, in case we want to save it. 
  name = ('snapshot_' + str(model) + '_' + znt(azm, 3) + '_' +
          znt(1e3*fdrive, 3) + 'mHz')
  if TP.savepath is not None:
    return PW.render(TP.savepath + name + '.pdf')
  else:
    return PW.render()

# =============================================================================
# =================================================== Bounce Frequency Profiles
# =============================================================================

def plotBounceFrequency(TP, allsix=False):

  # Figure out if we want to do all six models, or just the four we really use. 
  if allsix:
    TP.setPaths('/media/My Passport/RUNS/profiles/')
    PW = plotWindow(nrows=3, ncols=2, colorbar=None)
    models = (1, 2, 6, 5, 3, 4)
    rowlabels = ( notex('Day'), notex('Flank'), notex('Night') )
    kargs = {}
  else:
    PW = plotWindow(nrows=2, ncols=2, colorbar=None)
    models = (1, 2, 3, 4)
    rowlabels = ( notex('Day'), notex('Night') )
    kargs = {'fdrive':0.010, 'azm':1}

  # Set the labels and title. 
  collabels = ( notex('Active'), notex('Quiet') )
  title = notex('Alfv\\acute{e}n Bounce Frequency')
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)

  # Set ticks and labels. 
  xlims = (2, 10)
  xticks = (2, 4, 6, 8, 10)
  xticklabels = ('$2$', '', '$6$', '', '$10$')
  ylims = (0, 60)
  yticks = (0, 10, 20, 30, 40, 50, 60)
  yticklabels = ('$0$', '', '$20$', '', '$40$', '', '$60$')
  PW.setParams(xlims=xlims, xticks=xticks, xticklabels=xticklabels, 
               ylims=ylims, yticks=yticks, yticklabels=yticklabels)

  # Draw the lines. 
  for i, model in enumerate(models):
    path = TP.getPath(model=model, **kargs)
    PW[i].setParams( **TP.getCoords(path, 'L0', 'f') )
    va, dz = TP.getArray(path, 'va'), TP.getArray(path, 'dz')
    fbounce = 1e3/( 2*np.sum(dz/va, axis=1) )
    L0 = TP.getArray(path, 'L0')
    PW[i].setLine(L0, fbounce, 'b')
    # Draw Pc4 frequency range. 
    PW[i].setLine(L0, 7*np.ones(L0.shape), 'r:')
    PW[i].setLine(L0, 25*np.ones(L0.shape), 'r:')

  # Show or save the plot. 
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'fa.pdf')
  else:
    return PW.render()



def fa(TP, model=2):
  PW = plotWindow(colorbar=None)
  # Set the labels and title. 
  title = notex('Alfv\\acute{e}n Bounce Frequency: ') + tex(model)
  PW.setParams(title=title)

  # Set ticks and labels. 
  xlims = (2, 10)
  xticks = (2, 4, 6, 8, 10)
  xticklabels = ('$2$', '', '$6$', '', '$10$')
  ylims = (0, 60)
  yticks = (0, 10, 20, 30, 40, 50, 60)
  yticklabels = ('$0$', '', '$20$', '', '$40$', '', '$60$')
  PW.setParams(xlims=xlims, xticks=xticks, xticklabels=xticklabels, 
               ylims=ylims, yticks=yticks, yticklabels=yticklabels)

  # Draw the lines. 
  path = TP.getPath(model=model, fdrive=0.01, azm=1)
  PW.setParams( **TP.getCoords(path, 'L0', 'f') )
  va, dz = TP.getArray(path, 'va'), TP.getArray(path, 'dz')
  fbounce = 1e3/( 2*np.sum(dz/va, axis=1) )
  L0 = TP.getArray(path, 'L0')
  PW.setLine(L0, fbounce, 'b')
  # Draw Pc4 frequency range. 
  PW.setLine(L0, 7*np.ones(L0.shape), 'r:')
  PW.setLine(L0, 25*np.ones(L0.shape), 'r:')

  # Show or save the plot. 
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'fa.pdf')
  else:
    return PW.render()

# =============================================================================
# =============================================================== Toroidal Mode
# =============================================================================

def plotToroidal(TP):
  PW = plotWindow(nrows=1, ncols=-2, colorbar=None, joinlabel=False)
  # Set title and labels. 
  title = notex('Toroidal Resonance Structure')
  collabels = ( notex('First Harmonic'), notex('Second Harmonic') )
  PW.setParams(title=title, collabels=collabels)
  # Set ticks and limits. No tick labels. 
  xtks = np.linspace(-10, 10, 9)
  ytks = np.linspace(-4, 4, 5)
  xtls, ytls = (), ()
  xlms, ylms = (-10, 10), (-4, 4)
  xlbl, ylbl = notex('X'), notex('Z')
  PW.setParams(xlabel=xlbl, xlims=xlms, xticks=xtks, xticklabels=xtls, 
               ylabel=ylbl, ylims=ylms, yticks=ytks, yticklabels=ytls)
  # Draw Earth. This is a bit kludgey. 
  ax = ax = PW.cells.flatten()[0].ax
  ax.add_artist( Wedge( (0, 0), 1, 0, 360, fc='w' ) )
  # Draw a field line, leaving out the part inside Earth. 
  L, dL = 8, 0.2
  q0 = np.arcsin( np.sqrt(1./L) )
  q = np.linspace(q0, np.pi - q0, 150)
  r = L*np.sin(q)**2
  x, z = r*np.sin(q), r*np.cos(q)
  # Normalized coordinate for convenience. 
  u = np.pi*( q - q[0] )/( q[-1] - q[0] )
  for n in (1, 2):
    # Thickness of the line. Bottoms out at plus or minus 1. 
    w = 5*np.sin(n*u)
    w = np.where(np.abs(w)<1, np.sign(w), w)
    for i in range(1, u.size, 2):
      color = 'k' if w[i]>0 else 'k:'
      ax.plot( (-1)**n * x[i-1:i+2], z[i-1:i+2], color, linewidth=np.abs( w[i] ) )
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'toroidal.pdf')
  else:
    return PW.render()

# =============================================================================
# =============================================================== Poloidal Mode
# =============================================================================

def plotPoloidal(TP):
  PW = plotWindow(nrows=1, ncols=-2, colorbar=None, joinlabel=False)
  # Set title and labels. 
  title = notex('Poloidal Resonance Structure')
  collabels = ( notex('First Harmonic'), notex('Second Harmonic') )
  PW.setParams(title=title, collabels=collabels)
  # Set ticks and limits. No tick labels. 
  xtks = np.linspace(-10, 10, 9)
  ytks = np.linspace(-4, 4, 5)
  xtls, ytls = (), ()
  xlms, ylms = (-10, 10), (-4, 4)
  xlbl, ylbl = notex('X'), notex('Z')
  PW.setParams(xlabel=xlbl, xlims=xlms, xticks=xtks, xticklabels=xtls, 
               ylabel=ylbl, ylims=ylms, yticks=ytks, yticklabels=ytls)
  # Draw a field line, leaving out the part inside Earth. 
  L, dL = 8, 0.7
  q0 = np.arcsin( np.sqrt(1./L) )
  q = np.linspace(q0, np.pi - q0, 100)
  r = L*np.sin(q)**2
  x, z = r*np.sin(q), r*np.cos(q)
  [ PW.setLine(sign*x, z, 'k') for sign in (-1, 1) ]
  # Draw Earth. This is a bit kludgey. 
  ax = PW.cells.flatten()[0].ax
  ax.add_artist( Wedge( (0, 0), 1, 0, 360, fc='w' ) )
  # Normalized coordinate for convenience. 
  u = np.pi*( q - q[0] )/( q[-1] - q[0] )
  for i in range(2):
    dr = dL*np.sin( (i+1)*u )
    xm, zm = (-1)**(i+1)*(r - dr)*np.sin(q), (r - dr)*np.cos(q)
    xp, zp = (-1)**(i+1)*(r + dr)*np.sin(q), (r + dr)*np.cos(q)
    PW.setLine(xm, zm, 'b')
    PW.setLine(xp, zp, 'r')
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'poloidal.pdf')
  else:
    return PW.render()

# =============================================================================
# ================================================= Poloidal and Toroidal Modes
# =============================================================================

def plotPoloidalToroidal(TP):
  PW = plotWindow(nrows=2, ncols=-2, colorbar=None, joinlabel=False)
  # Set title and labels. 
  title = notex('Poloidal and Toroidal Polarizations')
  clabs = ( notex('First Harmonic'), notex('Second Harmonic') )
  rlabs = ( notex('Poloidal'), notex('Toroidal') )

  PW.setParams(title=title, collabels=clabs, rowlabels=rlabs)
  # Set ticks and limits. No tick labels. 
#  xtks = np.linspace(-10, 10, 9)
#  ytks = np.linspace(-4, 4, 5)
  xtls, ytls = (), ()
  xlms, ylms = (-10, 10), (-4, 4)
  xlbl, ylbl = notex('X'), notex('Z')
  PW.setParams(xlabel=xlbl, xlims=xlms, xticks=xlms, xticklabels=xtls, 
               ylabel=ylbl, ylims=ylms, yticks=ylms, yticklabels=ytls)
  # Draw a field line, leaving out the part inside Earth. 
  L, dL = 8, 0.7
  q0 = np.arcsin( np.sqrt(1./L) )
  q = np.linspace(q0, np.pi - q0, 100)
  r = L*np.sin(q)**2
  x, z = r*np.sin(q), r*np.cos(q)
  [ PW[0].setLine(sign*x, z, 'k') for sign in (-1, 1) ]
  # Draw Earth. This is a bit kludgey. 
  axes = [ cell.ax for cell in PW.cells.flatten() ]
  [ ax.add_artist( Wedge( (0, 0), 1, 0, 360, fc='w' ) ) for ax in axes ]
  # Normalized coordinate for convenience. 
  u = np.pi*( q - q[0] )/( q[-1] - q[0] )
  for i in range(2):
    dr = dL*np.sin( (i+1)*u )
    xm, zm = (-1)**(i+1)*(r - dr)*np.sin(q), (r - dr)*np.cos(q)
    xp, zp = (-1)**(i+1)*(r + dr)*np.sin(q), (r + dr)*np.cos(q)
    PW[0].setLine(xm, zm, 'k')
    PW[0].setLine(xp, zp, 'k')

  L, dL = 8, 0.2
  q0 = np.arcsin( np.sqrt(1./L) )
  q = np.linspace(q0, np.pi - q0, 150)
  r = L*np.sin(q)**2
  x, z = r*np.sin(q), r*np.cos(q)
  # Normalized coordinate for convenience. 
  u = np.pi*( q - q[0] )/( q[-1] - q[0] )
  for n in (1, 2):
    # Thickness of the line. Bottoms out at plus or minus 1. 
    w = 5*np.sin(n*u)
    w = np.where(np.abs(w)<1, np.sign(w), w)
    for i in range(1, u.size, 2):
      color = 'k' if w[i]>0 else 'k:'
      axes[1].plot( (-1)**n * x[i-1:i+2], z[i-1:i+2], color, linewidth=np.abs( w[i] ) )

  if TP.savepath is not None:
    return PW.render(TP.savepath + 'polarizations.pdf')
  else:
    return PW.render()

# =============================================================================
# ============================================ Odd and Even Harmonic Structures
# =============================================================================

def plotOddEven(TP):
  # Make a triple-wide window. 
  PW = plotWindow(nrows=2, ncols=-3, colorbar=None)
  # Set the ticks, etc, manually. 
  ytks = (0, 0.5, 1)
  ytls = [ '$' + notex(name) + '$' for name in ('S', '0', 'N') ]
  ylms = (0, 1)
  xlms = (-1, 1)
  xtks = -1 + 2*np.arange(1, 18, 2)/18.
  xtls = ( '$-\\pi$', '', '$-\\frac{\\pi}{2}$', '', '$0$', '',
           '$+\\frac{\\pi}{2}$', '', '$+\\pi$' )
  xlabel = notex('Time (\\frac{1}{\\omega})')
  PW.setParams(xlabel=xlabel, xlims=xlms, xticks=xtks, xticklabels=xtls, 
               xtickrelax=True, ylims=ylms, yticks=ytks, yticklabels=ytls)
  # Normalize the waves to the grid spacing. 
  mag = ( xtks[1] - xtks[0] )/2.5
  y = np.linspace(0, 1, 1000)
  # Iterate over the time steps. 
  for x in xtks:
    # Iterate over first and second harmonics. 
    for i in range(2):
      # Plot electric and magnetic perturbations. 
      B = mag*np.cos(np.pi*(i+1)*y)*np.sin( 2*np.pi*x/( xtks[-1] - xtks[0] ) )
      E = mag*np.sin(np.pi*(i+1)*y)*np.cos( 2*np.pi*x/( xtks[-1] - xtks[0] ) )
      PW[i].setLine(x + B, y, 'r')
      PW[i].setLine(x + E, y, 'b')
  # Put in a guide for the eye juse above the equator. 
  gx = np.linspace(xlms[0], xlms[1], 1000)
  gy = np.ones(gx.shape)*0.6
  PW.setLine(gx, gy, 'g:')
  # Plot a drift-bounce path. 
  dbx = xtks[4] - xtks[0]
  for x in xtks[0::8]:
    xup = np.linspace(x - dbx, x, 10)
    xdn = np.linspace(x, x + dbx, 10)
    yup = np.linspace(0.2, 0.8, 10)
    ydn = np.linspace(0.8, 0.2, 10)
    PW[1].setLine(xup, yup, 'm:')
    PW[1].setLine(xdn, ydn, 'm:')
  # Set title. Label rows. 
  title = notex('Electric (Blue) and Magnetic (Red) Harmonic Perturbations')
  rowlabels = ( notex('First') + '$\n$' + notex('Harmonic'), notex('Second') + '$\n$' + notex('Harmonic') )
  PW.setParams(title=title, rowlabels=rowlabels)
  # Create the plot. 
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'harmonics.pdf')
  else:
    return PW.render()

# =============================================================================
# ======================================= Compressional Alfven Cutoff Frequency
# =============================================================================

def plotAlfvenCutoff(TP, model=2):
  PW = plotWindow(nrows=1, ncols=3, colorbar='log', zmax=1e3)
  azms = (1, 8, 64)
  # Frequency doesn't matter. Just pick one. 
  fdrive = TP.getValues('fdrive')[0]
  collabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  title = notex('Compressional Alfv\\acute{e}n Cutoff Frequency')
  unitlabel = notex('mHz')
  PW.setParams(collabels=collabels, title=title, unitlabel=unitlabel)
  for i, azm in enumerate(azms):
    path = TP.getPath(azm=azm, fdrive=fdrive, model=model)
    PW[i].setParams( **TP.getCoords(path) )
    r = TP.getArray(path, 'r')
    va = TP.getArray(path, 'va')
    kperp = azm/(2*np.pi*r)
    PW[i].setContour( 1e3*kperp*va )
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'alfven_cutoff.pdf')
  else:
    return PW.render()

# =============================================================================
# ======================================================= Alfven Speed Profiles
# =============================================================================

def plotAlfvenSpeed(TP):
  PW = plotWindow(nrows=3, ncols=2, colorbar='log')
  # We don't actually care about fdrive, azm, or driving style. Just model. 
  azm = TP.getValues('azm')[0]
  fdrive = TP.getValues('fdrive')[0]
  for i, model in enumerate( (1, 2, 6, 5, 3, 4) ):
    path = TP.getPath(azm=azm, fdrive=fdrive, bdrive=0, model=model)
    PW[i].setParams( **TP.getCoords(path) )
    PW[i].setContour( 1000*TP.getArray(path, 'va') )
  # Set the labels and title. 
  collabels = [ notex('Active'), notex('Quiet') ]
  rowlabels = [ notex('Day'), notex('Flank'), notex('Night') ]
  title = notex('Alfv\\acute{e}n Speed')
  unitlabel = notex('\\frac{km}{s}')
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, 
               unitlabel=unitlabel)
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'va.pdf')
  else:
    return PW.render()



def va(TP, model=2):
  PW = plotWindow(colorbar='log')
  # We don't actually care about fdrive, azm, or driving style. Just model. 
  azm, fdrive = TP.getValues('azm')[0], TP.getValues('fdrive')[0]
  # Grab the profile. 
  path = TP.getPath(azm=azm, fdrive=fdrive, model=model)
  PW.setParams( **TP.getCoords(path) )
  PW.setContour( 1000*TP.getArray(path, 'va') )
  # Set the labels and title. 
  title = notex('Alfv\\acute{e}n Speed Profile: ') + tex(model)
  unitlabel = notex('\\frac{km}{s}')
  PW.setParams(title=title, unitlabel=unitlabel)
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'va.pdf')
  else:
    return PW.render()










# =============================================================================
# =========================================== Sym-H Index in Time and Frequency
# =============================================================================

from time import gmtime
from calendar import timegm

# Returns the time, in seconds, from 1970-01-01. 
def timeint(date=None, time=None):
  # Parse a string of the form hh:mm:ss. 
  if time is None:
    hh, mm, ss = 0, 0, 0
  else:
    # Account for missing colons and/or missing seconds. 
    hms = ( time.replace(':', '') + '00' )[:6]
    hh, mm, ss = int( hms[0:2] ), int( hms[2:4] ), int( hms[4:6] )
  # Parse a string of the form yyyy-mm-dd. If no date is given, use 
  # 1970-01-01 so that the returned value is just in seconds from midnight. 
  if date is None:
    year, mon, day = 1970, 1, 1
  else:
    # Allow missing dashes in the date. 
    ymd = date.replace('-', '')
    year, mon, day = int( ymd[0:4] ), int( ymd[4:6] ), int( ymd[6:8] )
  return timegm( (year, mon, day, hh, mm, ss) )

# Returns strings indicating the date and time. 
def timestr(ti):
  year, mon, day = gmtime(ti).tm_year, gmtime(ti).tm_mon, gmtime(ti).tm_mday
  hh, mm, ss = gmtime(ti).tm_hour, gmtime(ti).tm_min, gmtime(ti).tm_sec
  date = znt(year, 4) + '-' + znt(mon, 2) + '-' + znt(day, 2)
  time = znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
  return date, time

def plotSymh(TP, showline=False):
  # Two cells, side by side, sharing neither axis. 
  PW = plotWindow(nrows=1, ncols=2, colorbar=False)
  # Grab the data that we downloaded from NASA CDAWeb. 
  filename = '/home/user1/mceachern/Desktop/tuna-old/symh/symh_20130601.txt'
  data = [ line for line in open(filename, 'r').readlines() if line.strip() ]
  t, symh = np.zeros( len(data) ), np.zeros( len(data) )
  for i, line in enumerate(data):
    year, day, hour, minute, value = [ int(col) for col in line.split() ]
    t[i] = 24*60*day + 60*hour + minute
    symh[i] = value

  # Plot the Sym-H waveform in units of days. 
  date = ( t - t[0] )/1440.
  # Set the ticks manually. 
  ylabel = notex('Amplitude (nT)')
  ylims = (-150, 50)
  yticks = (-150, -100, -50, 0, 50)
  yticklabels = ('$-150$', '', '$-50$', '', '$+50$')
  xlabel = notex('Time (Days)')
  xlims = (0, 3)
  xticks = (0, 0.5, 1, 1.5, 2, 2.5, 3)
  xticklabels = ('$0$', '', '$1$', '', '$2$', '', '$3$')
  # Plot. 
  PW[0].setParams(xlabel=xlabel, xlims=xlims, xticks=xticks, 
                  xticklabels=xticklabels, ylabel=ylabel, ylims=ylims, 
                  yticks=yticks, yticklabels=yticklabels)
  PW[0].setLine(date, symh, 'k')

  # Get the Fourier transform of the Sym-H waveform. Throw a factor of 1/N in
  # there, since Numpy normalizes backwards. 
  symhfft = np.fft.rfft(symh)[:-1] / date.size
  # Scale days to ks to get frequency in mHz. 
  dt = ( date[1] - date[0] )*86400/1e3
  freq = np.fft.fftfreq(symh.size, d=dt)[:symhfft.size]
  # Plot log-log, skipping the DC offset. 
  logf, logs,  = np.log10( freq[1:] ), np.log10( np.abs(symhfft)[1:] )
  # Set ticks manually. 
  lims = (-3, 1)
  ticks = (-3, -2, -1, 0, 1)
  ticklabels = ('$-3$', '', '$-1$', '', '$+1$')
  ylabel = notex('Log Amplitude (nT)')
  xlabel = notex('Log Frequency (mHz)')
  PW[1].setParams(xlabel=xlabel, xlims=lims, xticks=ticks, 
                  xticklabels=ticklabels, ylabel=ylabel, ylims=lims, 
                  yticks=ticks, yticklabels=ticklabels, axright=True)
  PW[1].setLine(logf, logs)

  # Also plot a fit. Really just a red line along the top of the distribution. 
  power = -1
  intercept = 0.1
  fit = intercept * freq[1:]**power
  PW[1].setLine(logf, np.log10(fit), color='r')

  # Set title and labels. 
  title = notex('Sym-H for June 2013 Storm')
  colLabels = ( notex('Sym-H'), notex('Sym-H Fourier Transform') )
  PW.setParams(title=title, collabels=colLabels)

  if TP.savepath is not None:
    return PW.render(TP.savepath + 'symh.pdf')
  else:
    return PW.render()

# =============================================================================
# ============================================================ Plasma Frequency
# =============================================================================

def plotPlasmaFrequency(TP):
  PW = plotWindow(nrows=1, ncols=1, colorbar='log')
  # Set title and labels. 
  title = notex('Plasma Frequency')
  # Grab whichever path. 
  path = TP.paths.keys()[0]
  X, Z = TP.getArray(path, 'X'), TP.getArray(path, 'Z')
  coords = TP.getCoords(path)
  # This isn't a cramped diagram, so let's touch up the ticks. This should
  # probably go in plotmod, to ensure that this plot always matches the grid. 
  coords['xticks'] = np.arange(11)
  coords['xticklabels'] = [ '$' + znt(x) + '$' for x in coords['xticks'] ]
  coords['xticklabels'][1::2] = ['']*len( coords['xticklabels'][1::2] )
  coords['yticks'] = (-4, -3, -2, -1, 0, 1, 2, 3, 4)
  coords['yticklabels'] = ['']*len( coords['yticks'] )
  coords['yticklabels'][::2] = ('$-4$', '$-2$', '$0$', '$+2$', '$+4$')
  PW.setParams( **coords )
  # Grab the number density, and use it to compute the plasma frequency. 
  n = TP.getArray(path, 'n')
  # Scaling to Hz lines up better with the color cutoffs. 
  PW.setContour( np.sqrt( n*phys.qe**2 / (phys.me*phys.eps0) )/(2*np.pi) )
  unitlabel = notex('Hz')
#  PW.setContour( np.sqrt( n*phys.qe**2 / (phys.me*phys.eps0) ) )
#  unitlabel = notex('\\frac{rad}{s}')
  PW.setParams(title=title, unitlabel=unitlabel)
  # Show or save the plot. 
  if TP.savepath is not None:
    return PW.render( TP.savepath + 'op.pdf' )
  else:
    return PW.render()

# =============================================================================
# ================================= Compressional Waves Failing to Propagate In
# =============================================================================

def plotBdrive(TP, fdrive):
  # Let's look at the runs driven with magnetic compression.  
  TP.setPaths('/media/My Passport/RUNS/BDRIVE/')

  # Set up the window. 
  azms = (1, 4, 16, 64)
  models = (1, 2, 3, 4)

#  PW = plotWindow(nrows=len(azms), ncols=len(models), colorbar='sym')
  PW = plotWindow(nrows=len(azms), ncols=len(models), colorbar='log')

  # Set title and labels. 
  title = notex('Mean Energy Density: ') + tex(fdrive) + notex('Compression')
  unitlabel = notex('\\frac{nJ}{m^3}')
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = [ tex(model) for model in models ]
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels, 
               unitlabel=unitlabel, zmax=1)

  # Each cell is a different run. 
  for row, azm in enumerate(azms):

    for col, model in enumerate(models):

      path = TP.getPath(model=model, fdrive=fdrive, azm=azm)

      PW[:, :].setParams( **TP.getCoords(path) )

#      # Toroidal Poynting flux snapshot? 
#      Stor = TP.getArray(path, 'Stor')[:, :, 299]
#      PW[row, col].setContour(Stor)

      # Mean energy density? 
      u = np.mean(TP.getArray(path, 'u'), axis=-1)
      PW[row, col].setContour(u)

  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'bdrive_' + znt(1e3*fdrive) + 'mHz'
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# ==================================================================== The Grid
# =============================================================================

def plotGrid(TP):
  PW = plotWindow(ncols=1, nrows=1, colorbar=False)
  PW.setParams( title=notex('Nonorthogonal Dipole Grid') )
  path = TP.paths.keys()[0]
  X, Z = TP.getArray(path, 'X'), TP.getArray(path, 'Z')
  coords = TP.getCoords(path)
  # This isn't a cramped diagram, so let's touch up the ticks. 
  coords['xticks'] = np.arange(11)
  coords['xticklabels'] = [ '$' + znt(x) + '$' for x in coords['xticks'] ]
  coords['xticklabels'][1::2] = ['']*len( coords['xticklabels'][1::2] )
  coords['yticks'] = (-4, -3, -2, -1, 0, 1, 2, 3, 4)
  coords['yticklabels'] = ['']*len( coords['yticks'] )
  coords['yticklabels'][::2] = ('$-4$', '$-2$', '$0$', '$+2$', '$+4$')
  PW.setParams( **coords )
  nx, nz = X.shape
  stride = 5
  [ PW.setLine( X[i, :], Z[i, :], 'k' ) for i in range(0, nx, stride) ]
  [ PW.setLine( X[:, k], Z[:, k], 'k' ) for k in range(0, nz, stride) ]
  if TP.savepath is not None:
    return PW.render(TP.savepath + 'grid.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================= Electron Inertial Length Comparison
# =============================================================================

def plotZoom(TP, azm, step=99):
  fdrive = 0.016
  model = 2
  lmin = 5

  PW = plotWindow(ncols=2, nrows=1, colorbar='sym')

  title = notex('Parallel Electric Fields: ' + tex(model) + ', ' +
                str(step+1) + 's \\,of 16mHz Current, ') + 'm = ' + znt(azm)
  unitlabel = notex('\\frac{mV}{m}')
  collabels = ('\\delta x > \\frac{c}{\\omega_P}', 
               '\\delta x < \\frac{c}{\\omega_P}')

  xlims = (62, 68)
  xticks = (62, 63, 64, 65, 66, 67, 68)
  xticklabels = ('$62^\\circ$', '', '', '$65^\\circ$', '', '', '$68^\\circ$')

  ylims = (100, 1000)
  yticks = (100, 250, 400, 550, 700, 850, 1000)
  yticklabels = ('$100$', '', '$400$', '', '$700$', '', '$1000$')
  PW.setParams(xlims=xlims, xticks=xticks, xticklabels=xticklabels,ylims=ylims, 
               yticks=yticks,yticklabels=yticklabels, title=title, 
               collabels=collabels, unitlabel=unitlabel)

  TP.setPaths('/media/My Passport/RUNS/LMIN_LMAX/inertia_on')
  path = TP.getPath(azm=azm)
  Ez = TP.getArray(path, 'Ez')[:, :, step]
  PW[1].setContour(Ez)
  PW.setParams( outline=True, **TP.getCoords(path, 'lat', 'alt') )

  TP.setPaths('/media/My Passport/RUNS/INERTIA')
  path = TP.getPath(azm=azm, model=model, fdrive=fdrive)
  Ez = TP.getArray(path, 'Ez')[:, :, step]
  PW[0].setContour(Ez)
  PW[0].setParams( **TP.getCoords(path, 'lat', 'alt') )

  if TP.savepath is not None:
    return PW.render(TP.savepath + 'inertial_length.pdf')
  else:
    return PW.render()

# =============================================================================
# ============================== Compare Divergence of Poynting Flux to J dot E
# =============================================================================

def plotDivFlux(TP, model, fdrive):
  TP.setPaths('/media/My Passport/RUNS/INERTIA/')
  # Set up the window. 
  azms = (1, 4, 16, 64)
  PW = plotWindow(nrows=len(azms), ncols=4, colorbar='sym', zmax=10)
  # Set title and labels. 
  title = ( notex('Power Density at ') + 'R_I' + notex(': ') + tex(model) +
            notex(', ') + tex(fdrive) + notex(' Current') )
  unitlabel = notex('\\frac{nW}{m^3}')
  rowlabels = [ 'm \\! = \\! ' + str(azm) for azm in azms ]
  collabels = ( '\\nabla_\\bot\\cdot{}\\underline{S}_\\bot', 
                '\\frac{\\partial}{\\partial{}z}S_z', 
                '\\underline{J}_\\bot\\cdot\\underline{E}_\\bot', 
                'J_z\\,E_z' )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels, 
               unitlabel=unitlabel)
  # Each row comes from a different run. 
  for row, azm in enumerate(azms):
    path = TP.getPath(model=model, fdrive=fdrive, azm=azm)
    PW[row, :].setParams( **TP.getCoords(path, 't', 'lat0') )
    # J dot E. 
    JEx = ( TP.getArray(path, 'Jx')*TP.getArray(path, 'Ex') )[:, 0, :]
    JEy = ( TP.getArray(path, 'Jy')*TP.getArray(path, 'Ey') )[:, 0, :]
    JEz = ( TP.getArray(path, 'Jz')*TP.getArray(path, 'Ez') )[:, 0, :]
    # Differential lengths. 
    dx = TP.getArray(path, 'dx0')
    dy = azm*TP.getArray(path, 'dy0')
    dz = TP.getArray(path, 'dz0')
    # Poynting flux.
    Sx = TP.getArray(path, 'Sx')
    Sy = TP.getArray(path, 'Sy')
    Sz = TP.getArray(path, 'Sz')
    # Div S. 
    dSx = ( Sx[2:, 0, :] - Sx[:-2, 0, :] ) / ( 2*dx[1:-1, None] )
    dSy = Sy[:, 0, :]/dy[:, None]
    dSz = ( Sz[:, 1, :] - Sz[:, 0, :] ) / dz[:, None]
    dSxy = dSy
    dSxy[1:-1, :] = dSxy[1:-1, :] + dSx
    # Add them to the plot. 
    PW[row, 0].setContour(dSxy)
    PW[row, 1].setContour(dSz)
    PW[row, 2].setContour(JEx + JEy)
    PW[row, 3].setContour(JEz)
  # Show or save the plot. 
  if TP.savepath is not None:
    name = 'JE_' + znt(1e3*fdrive) + 'mHz_' + str(model)
    return PW.render( TP.savepath + name + '.pdf' )
  else:
    return PW.render()

# =============================================================================
# =========================================== Ionospheric Conductivity Profiles
# =============================================================================

def plotSigma(TP):
  # Create the window. 
  PW = plotWindow(nrows=3, ncols=2, colorbar=None)
  # For this plot, we don't actually need the 2D arrays that Tuna spits out. We
  # can just read in the profiles directly. 
  for i, model in enumerate( (1, 2, 6, 5, 3, 4) ):
    datname = './models/ionpar' + str(model) + '.dat'
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
  rowLabels = [ notex('Day'), notex('Flank'), notex('Night') ]
  title = notex('Pedersen (Blue), Hall (Red), and Parallel (Green) ' + 
                'Conductivities')
  xlabel = notex('Log Conductivity (\\frac{S}{m})')
  ylabel = notex('Log Altitude  (km)')
  # Set the ticks manually. 
  xlims = (-15, 5)
  xticks = (-15, -10, -5, 0, 5)
  xticklabels = ('$-15$', '', '$-5$', '', '$+5$')
  ylims = (2, 4)
  yticks = (2, 2.5, 3, 3.5, 4)
  yticklabels = ('$2$', '', '$3$', '', '$4$')
  PW.setParams(collabels=colLabels, rowlabels=rowLabels, xlabel=xlabel, 
               xlims=xlims, xticks=xticks, xticklabels=xticklabels, 
               ylabel=ylabel, ylims=ylims, yticks=yticks, 
               yticklabels=yticklabels, title=title)
  # Show or save the plot. 
  if TP.savepath is not None:
    return PW.render( TP.savepath + 'sigma.pdf' )
  else:
    return PW.render()

# #####################################################################
# ################################################### For Importability
# #####################################################################

if __name__=='__main__':
  main()

