# outer __init__.py
'''
Kernel module
@summary: Contains all the routines necessary to the wavelet analysis<br />
          of along-track altimetry data, and to compute diagnostics from it.
@author: Renaud DUSSURGET, LER/PAC IFREMER.
@change: Create in November 2012 by RD.
'''
from runAnalysis import runAnalysis
from external import wavelet
from detectEddies import _2D as _2Ddetection
from getScales import decorrelation_scale, solid_body_scale, eddy_amplitude, cyclone
from bins import grid_space, grid_var, grid_time
#from detectEddies import oneD as eddies_1D