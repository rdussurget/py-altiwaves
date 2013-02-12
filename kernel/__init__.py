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
from detectEddies import detection
from getScales import get_characteristics, eddy_amplitude, cyclone
from spectrum import spectral_analysis, periodogram_analysis
from bins import bin_space, bin_time#grid_space, grid_var, grid_time
from io import save_analysis, save_detection,save_binning
#from detectEddies import oneD as eddies_1D