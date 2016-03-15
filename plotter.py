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

#  # Plasma frequency profile. 
#  plotPlasmaFrequency(TP)

#  # Alfven speed profiles. 
#  plotAlfvenSpeed(TP)

  # Compressional Alfven frequency cutoff. 
  plotAlfvenCutoff(TP, model=2)

#  # Snapshots of the electric field components. 
#  plotSnapshots(TP, model=2, fdrive=0.010, azm=16)

#  # Comparison of parallel current and Poynting flux. 
#  plotCurrentFlux(TP, model=1, fdrive=0.016)

#  # Comparison of J dot E to the divergence of the Poynting flux. 
#  plotDivFlux(TP, model=1, fdrive=0.016)

#  # Zoom way in on runs with and without inertial length scales. 
#  plotZoom(TP, azm=16)

#  # Conductivity profiles. 
#  plotSigma(TP)

#  # The grid. 
#  plotGrid(TP)

#  # Sym-H index and its Fourier transform. 
#  plotSymh(TP)

#  # Show waves failing to propagate in when driven from the outer boundary. 
#  plotBdrive(TP, fdrive=0.022)

#  # Plot magnetic field signatures at the ground. 
#  for kargs in loopover( fdrive=TP.getValues('fdrive') ):
#    plotGround(TP, **kargs)

#  # Plot the radial distribution in energy. 
#  for kargs in loopover( model=(1, 2), fdrive=TP.getValues('fdrive') ):
#  for kargs in loopover( model=(1,), fdrive=(0.019,) ):
#    plotLayers(TP, **kargs)

#  for kargs in loopover( model=TP.getValues('model') ):
#    plotEnergy(TP, **kargs)

  return

# #############################################################################
# ########################################################### Plotting Routines
# #############################################################################

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
  PW = plotWindow(nrows=2, ncols=2, colorbar='log')
  # We don't actually care about fdrive, azm, or driving style. Just model. 
  azm = TP.getValues('azm')[0]
  fdrive = TP.getValues('fdrive')[0]
  for i in range(4):
    path = TP.getPath(azm=azm, fdrive=fdrive, bdrive=0, model=i+1)
    PW[i].setParams( **TP.getCoords(path) )
    PW[i].setContour( 1000*TP.getArray(path, 'va') )
  # Set the labels and title. 
  collabels = [ notex('Active'), notex('Quiet') ]
  rowlabels = [ notex('Day'), notex('Night') ]
  title = notex('Alfv\\acute{e}n Speed')
  unitlabel = notex('\\frac{km}{s}')
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, 
               unitlabel=unitlabel)
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
# ================================ Parallel Current and Poynting Flux Snapshots
# =============================================================================

def plotCurrentFlux(TP, model, fdrive):
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
  collabels = ( tex('Ex'), tex('Ey'), 
               '10^{' + znt(scalefac) + '} \\times ' + tex('Ez') )
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
# =========================================== Ionospheric Conductivity Profiles
# =============================================================================

def plotSigma(TP):
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

# #####################################################################
# ################################################### For Importability
# #####################################################################

if __name__=='__main__':
  main()

