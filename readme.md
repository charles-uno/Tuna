# Tuna

Tuna is a 2.5D simulation of electromagnetic waves in Earth's magnetosphere. It's written by Charles McEachern, based on work by Bob Lysak. 

## `source.f90`

The code is written in Fortran. All modules are packed into a single file. 

## `driver.py`

The driver is a Python script which conducts one or more runs with Tuna, based on the input parameters at the top of the script. 

For each run, the driver creates a timestamped output directory. It then compiles the source code, copies over the necessary ionospheric profiles, and creates the input file `params.in`. If running on a local machine, the driver carries out the run and checks for a crash; if run on a supercomputer, it instead creates a PBS script and submits the run to the queue. 

## `pickler.py`

The pickler is a Python script which sits between the driver and the plotter. It reads in Tuna's Fortran-style output, converts it into Numpy arrays, and outputs those arrays as pickle files. This significantly decreases the size of the data, while also dramatically increasing the speed at which data can be read in and plotted. 

## `plotter.py`

The plotter is a Python script which creates visualizations from Tuna's output. It's a bit kludgey right now. 




