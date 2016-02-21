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

  plotSigma()

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

# Columns are field: active Bq, quiet Bq, active Bf, quiet Bf. 
# Rows are modenumber: 1, 2, 4, 8, 16, 32, 64. 



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


# #############################################################################
# ###################################################### LaTeX Helper Functions
# #############################################################################

def notex(x):
  return '\\operatorname{' + x.replace(' ', '\\;') + '}'

def tex(x):
  # Format frequencies nicely. 
  if isinstance(x, float):
    return notex(format(1e3*x, '.0f') + 'mHz ')
  # Check the input against our dictionary of LaTeX strings. 
  texdict = {
             # Spell out what each model means. 
             1:notex('Active Day'),
             2:notex('Quiet Day'),
             3:notex('Active Night'),
             4:notex('Quiet Night'),
             # Electric and magnetic field components. 
             'Bf':'B_\\phi', 
             'Bq':'B_\\theta', 
             'Bx':'B_x', 
             'By':'B_y', 
             'Bz':'B_z', 
             'Ex':'E_x', 
             'Ey':'E_y', 
             'Ez':'E_z', 
             # Axis labels. 
             'alt':notex('Altitude (km)'), 
             'lat':notex('Latitude (^\\circ)'), 
             't':notex('Time (s)'), 
             'logU':notex('Log_{10}') + 'U', 
             'X':notex('X (R_E)'), 
             'Z':notex('Z (R_E)')
            }
  return '?' if x not in texdict else texdict[x]

# #####################################################################
# ################################################### For Importability
# #####################################################################

if __name__=='__main__':
  main()

