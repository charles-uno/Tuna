#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# Note: This document wraps at column 80. 

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# #############################################################################
# ##################################################### Import Python Libraries
# #############################################################################

from plotmod import *

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

#  TP = tunaPlotter()

#  for kargs in loopover( fdrive=TP.getValues('fdrive'), mode=('BE', 'PT'), pick=('-i' not in argv) ):
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
    # Chop off the altitudes at 100km and 10000km. 
    bottom = np.argmin(ionos[:, 0]<100)
    if np.max( ionos[:, 0] )>1e4:
      top = np.argmax(ionos[:, 0]>1e4)
    else:
      top = len( ionos[:, 0] )
    # Log altitude, rather than altitude on a log scale. 
    logalt = np.log10( ionos[bottom:top, 0] )
    # Log conductivity instead of conductivity on a log scale. Values are given
    # in mS/m; plot in S/m. 
    logsig0 = np.log10( 1e6/( phys.mu0*ionos[bottom:top, 4] ) ) - 3
    logsigH = np.log10( 4*np.pi*phys.eps0*ionos[bottom:top, 3] ) - 3
    logsigP = np.log10( 4*np.pi*phys.eps0*ionos[bottom:top, 2] ) - 3
    # Add the lines to the plot. 
    PW[i].setLine(logsigP, logalt, 'r')
    PW[i].setLine(logsigH, logalt, 'b')
    PW[i].setLine(logsig0, logalt, 'g')
  # Set the labels and title. 

  colLabels = [ self.texText('Active'), self.texText('Quiet') ]

  rowLabels = [ self.texText('Day'), self.texText('Night') ]

  title = self.texText('Pedersen (Blue), Hall (Red), and Parallel (Green) Conductivities')

  xlabel = tex('Log_{10} Conductivity (\\frac{S}{m})')
  xlims = (-15, 5)

  ylabel = tex('Log_{10} Altitude (km)')

  PW.setParams(collabels=colLabels, rowlabels=rowLabels, nxticks=5, xlabel=xlabel, xlims=xlims, ylabel=self.texLabel('logalt'), title=title)

  if self.savepath is not None:
    return PW.render(self.savepath + 'sigma.pdf')
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

