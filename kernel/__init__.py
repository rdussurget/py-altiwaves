# outer __init__.py
from runAnalysis import runAnalysis
import wavelet 
from detectEddies import _2D as _2Ddetection
from getScales import decorrelation_scale, solid_body_scale, eddy_amplitude, cyclone
from bins import grid, grid_var, grid_time
#from detectEddies import oneD as eddies_1D