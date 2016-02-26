#!/usr/bin/env python

# Charles McEachern

# Sprint 2016

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

try:
  import cPickle as pickle
except ImportError:
  import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import io
from subprocess import Popen, PIPE

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  for sat in ('a', 'b'):
    raw = '/home/user1/mceachern/Desktop/RBSP/raw_events_' + sat + '.txt'
    out = '/home/user1/mceachern/Desktop/RBSP/events_' + sat + '.txt'
    if os.path.exists(out):
      print 'list of events already exists for rbsp ' + sat
    else:
      print 'creating list of events for rbsp ' + sat + '...',
      events = sum( [ line.lstrip('>').split() for line in read(raw) ], [] )
      print len(events), 'found'
      # Don't worry about end time. Events are 10 minutes long by construction. 
      [ append(x, '/home/user1/mceachern/Desktop/RBSP/events_' + sat + '.txt') for x in events ]


  module('load idl')

  for sat in ('a', 'b')[:1]:
    events = read('/home/user1/mceachern/Desktop/RBSP/events_' + sat + '.txt')[:1]
    for event in events:
      print 'event: ', event
      name = event.replace('/', '_').replace('-', '').replace(':', '')

      '''
      if os.path.exists('temp.pro'):
        os.remove('temp.pro')
#      print '\tcreating pro script to download the data as a sav file'
      append( pro(sat=sat, event=event) )
      print '\tDownloading the data using SPEDAS... '
      out, err = bash('idl -e @temp -IDL_PATH +~/Desktop/RBSP/packages/:<IDL_DEFAULT>')
      print '\tReading the IDL data into Python...'
      sav = io.readsav('/home/user1/mceachern/Desktop/RBSP/temp.sav')
      '''

      pklpath = '/media/My Passport/RBSP/pickles/' + name + '.pkl'

      '''
      # Tuples are safer to pickle than dictionaries. 
      print '\tCreating a pickle: ', pklpath
      with open(pklpath, 'wb') as handle:
        pickle.dump(sav.items(), handle, protocol=-1)
#      print '\tsanity check... '
#      with open(pklpath, 'rb') as handle:
#        x = dict( pickle.load(handle) )
#      print all( sav['time'] == x['time'] )
#      print all( sav['bgse'].flatten() == x['bgse'].flatten() )
      '''

  # Now let's look at an event. 
  with open(pklpath, 'rb') as handle:
    x = dict( pickle.load(handle) )

  t = x['time']
  # Ten minutes compared to a day. 
  N = (600*t.size)/(24*3600)
  t = t[:N] - t[0]

  BX = x['bgse'][0][:N]
  BY = x['bgse'][1][:N]
  BZ = x['bgse'][2][:N]

  BX = BX - np.mean(BX)
  BY = BY - np.mean(BY)
  BZ = BZ - np.mean(BZ)

  plt.plot(t, BX)
  plt.plot(t, BY)
  plt.plot(t, BZ)

  plt.show()


  return

# #############################################################################
# ########################################################### 
# #############################################################################

def pro(sat, event):
  return ( 'rbsp_efw_init\n' +
           'timespan,\'' + event + '\'\n' +
           'rbsp_load_emfisis,probe=' + sat + ',coord=\'gse\',cadence=\'hires\',level=\'l3\'\n' + 
           'get_data,1,time,bgse\n' + 
           'save,time,bgse,filename=\'~/Desktop/RBSP/temp.sav\'' )

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

def append(text, filename=None):
  if filename is not None:
    with open(filename, 'a') as fileobj:
      fileobj.write(text + '\n')
  return text

def read(filename):
  with open(filename, 'r') as fileobj:
    return [ x.strip() for x in fileobj.readlines() ]

def bash(command, save='stdoe.txt'):
  out, err = Popen(command.split(), stdout=PIPE, stderr=PIPE).communicate()
  return append(out, save), append(err, save)

def module(command, save='stdoe.txt'):
  out, err = bash('/usr/bin/modulecmd python ' + command, save=save)
  exec out
  return err

# ########################################################### For Importability
# #############################################################################
# #############################################################################

if __name__=='__main__':
  main()


