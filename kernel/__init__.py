# outer __init__.py
'''
Kernel module
@summary: Contains all the routines necessary to the wavelet analysis<br />
          of along-track altimetry data, and to compute diagnostics from it.
@author: Renaud DUSSURGET, LER/PAC IFREMER.
@change: Create in November 2012 by RD.
@copyright: Renaud Dussurget 2012.
@license: GNU Lesser General Public License
    
    This file is part of PyAltiWAVES.
    
    PyAltiWAVES is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.
    PyAltiWAVES is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
    for more details.
    
    You should have received a copy of the GNU Lesser General Public License along
    with PyAltiWAVES.  If not, see <http://www.gnu.org/licenses/>.
'''
from runAnalysis import runAnalysis
from external import wavelet
from detectEddies import detection
from getScales import get_characteristics, eddy_amplitude, cyclone
from spectrum import spectral_analysis, periodogram_analysis
from bins import bin_space, bin_time#grid_space, grid_var, grid_time
from io import save_analysis, save_detection,save_binning
#from detectEddies import oneD as eddies_1D