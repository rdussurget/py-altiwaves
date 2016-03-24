# -*- coding: latin-1 -*-
'''
kernel.io module
@summary : I/O tools for wavelet analysis
@requires: altimetry.tools.nctools
@since: Created on 6 déc. 2012
@author: rdussurg
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
import os
import numpy as np
from collections import OrderedDict
from altimetry.tools.nctools import nc
if __debug__ : import matplotlib.pyplot as plt

def save_analysis(filename,
                  lon,
                  lat,
                  time,
                  cycle,
                  track,
                  sla,
                  sa_spectrum,
                  sa_lscales,
                  wvsla,
                  daughter,
                  dimensions=None,
                  attributes=None,
                  verbose=1,
                  clobber=True):
    
    #Init netcdf object
    obj=nc(verbose=verbose) #Init NetCDF object
    
    #Setup dimensions if required
    if dimensions is None :
        nx=len(lon)
        ncyc=len(cycle)
        nt=len(track)
        dimensions=OrderedDict({'_ndims':3,'record':0,'cycle':ncyc,'track':nt})
    
    #Setup default attributes
    attrStr=OrderedDict()
    attrStr['description']='Along-track wavelet analysis'
    attrStr['input']=filename
    if attributes is not None : #Append attributes if asked
        for a in [k for k in attributes.keys() if (k is not 'description') and (k is not 'input')] :
            attrStr.update({a:attributes[a]})
    
    #Setup output structure
    dataStr=OrderedDict()
    dataStr['_dimensions']=dimensions
    dataStr['_attributes']=attrStr
    dataStr['record']={'data':np.arange(nx),'long_name':'record_number','units':'1','_dimensions':('record',)}
    dataStr['lon']={'data':lon,'long_name':'longitude','units':'degrees_north','_dimensions':('record',)}
    dataStr['lat']={'data':lat,'long_name':'latitude','units':'degrees_east','_dimensions':('record',)}
    dataStr['time']={'data':time,'long_name':'time of measurement','units':'days since 1950-01-01 00:00:00UTC','_dimensions':('cycle',)}
    dataStr['cycle']={'data':cycle,'long_name':'cycle number','units':'1','_dimensions':('cycle',)}
    dataStr['track']={'data':track,'long_name':'track number','units':'1','_dimensions':('track',)}
    dataStr['sla']={'data':sla,'long_name':'sea_level_anomaly','units':'m','_dimensions':('cycle','record')}
    dataStr['sa_spectrum']={'data':sa_spectrum,'long_name':'scale_averaged_spectrum','units':'m^^2','_dimensions':('cycle','record')}
    dataStr['sa_lscales']={'data':sa_lscales,'long_name':'spatial_lengthscale','units':'km','_dimensions':('cycle','record')}
    dataStr['wvsla']={'data':wvsla,'long_name':'filtered_sea_level_anomaly','units':'m','_dimensions':('cycle','record')}
    dataStr['daughter']={'data':daughter,'long_name':'most_energetic_wavelet','units':'m','_dimensions':('cycle','record')}
    
    #Write data
    res=obj.write(dataStr,filename,clobber=clobber)
    
    return res

def save_detection(filename,
                   eind,
                   lon,
                   lat,
                   amplitude, diameter, relvort, ugdiameter, ugamplitude, rk_relvort, rk_center, self_advect,
                   dimensions=None,
                   attributes=None,
                   verbose=1,
                   clobber=True):
    
    #Init netcdf object
    obj=nc(verbose=verbose) #Init NetCDF object
    
    #Setup dimensions if required
    if dimensions is None :
        nx=len(amplitude)
        dimensions=OrderedDict({'_ndims':1,'record':0})
    
    #Setup default attributes
    attrStr=OrderedDict()
    attrStr['description']='Eddy-like features detected from wavelet analysis'
    attrStr['input']=filename
    if attributes is not None : #Append attributes if asked
        for a in [k for k in attributes.keys() if (k is not 'description') and (k is not 'input')] :
            attrStr.update({a:attributes[a]})
    
    #Setup output structure
    dataStr=OrderedDict()
    dataStr['_dimensions']=dimensions
    dataStr['_attributes']=attrStr
    dataStr['record']={'data':np.arange(nx),'long_name':'record_number','units':'1','_dimensions':('record',)}
    dataStr['xind']={'data':eind[0,:],'long_name':'along_track_index','units':'1','_dimensions':('record',)}
    dataStr['yind']={'data':eind[1,:],'long_name':'time_index','units':'1','_dimensions':('record',)}
    dataStr['lon']={'data':lon[eind[0,:]],'long_name':'longitude','units':'degrees_north','_dimensions':('record',)}
    dataStr['lat']={'data':lat[eind[0,:]],'long_name':'latitude','units':'degrees_east','_dimensions':('record',)}
    dataStr['amplitude']={'data':amplitude,'long_name':'amplitude','units':'cm','_dimensions':('record',)}
    dataStr['diameter']={'data':diameter,'long_name':'diameter','units':'km','_dimensions':('record',)}
    dataStr['relvort']={'data':relvort,'long_name':'relative_vorticity','units':'s-1','_dimensions':('record',)}
    dataStr['ugdiameter']={'data':ugdiameter,'long_name':'eddy_core_diameter','units':'km','_dimensions':('record',)}
    dataStr['ugamplitude']={'data':ugamplitude,'long_name':'eddy_core_amplitude','units':'cm','_dimensions':('record',)}
    dataStr['rk_relvort']={'data':rk_relvort,'long_name':'rankine_eddy_vorticity','units':'s-1','_dimensions':('record',)}
    dataStr['rkxind']={'data':rk_center,'long_name':'rankine_eddy_index','units':'1','_dimensions':('record',)}
    dataStr['rk_lon']={'data':lon[rk_center],'long_name':'rankine_eddy_longitude','units':'degrees_north','_dimensions':('record',)}
    dataStr['rk_lat']={'data':lat[rk_center],'long_name':'rankine_eddy_latitude','units':'degrees_north','_dimensions':('record',)}
    dataStr['advection']={'data':self_advect,'long_name':'eddy_self_advection','units':'m.s-1','_dimensions':('record',)}
    
    #Write data
    res=obj.write(dataStr,filename,clobber=clobber)
    
    return res

def save_binning(filename,
                 blon, blat, hist,  ampmn, lenmn, rvmn, amprms, lenrms, rvrms,  \
                 btime, thist, tampmn, tlenmn, trvmn, tamprms, tlenrms, trvrms, \
#                 blon, blat, hist,  ampmn, lenmn, rvmn, uglmn, wvlenmn, amprms, lenrms, rvrms, uglrms, wvlenrms,
#                 btime, thist, tampmn, tlenmn, trvmn, tuglmn, twvlmn, tamprms, tlenrms, trvrms, tuglrms, twvlrms,
                 dimensions=None,
                 attributes=None,
                 description='Climatology of eddy-like processes variability',
                 verbose=1,
                 clobber=True):
    
    
    #Init netcdf object
    obj=nc(verbose=verbose) #Init NetCDF object
    
    #Setup dimensions if required
    if dimensions is None :
        nx=len(blat)
        nt=len(btime)
        sdimensions=OrderedDict({'_ndims':1,'lat':0})
        tdimensions=OrderedDict({'_ndims':1,'time':0})
    
    #Setup default attributes
    attrStr=OrderedDict()
    attrStr['description']=description+' - space'
    if attributes is not None : #Append attributes if asked
        for a in [k for k in attributes.keys() if (k is not 'description') and (k is not 'input')] :
            attrStr.update({a:attributes[a]})
    
    #Setup output structure
    sStr=OrderedDict()
    sStr['_dimensions']=sdimensions
    sStr['_attributes']=attrStr
    sStr['lat']={'data':blat,'long_name':'latitude','units':'degrees_east','_dimensions':('lat',)}
    sStr['lon']={'data':blon,'long_name':'longitude','units':'degrees_north','_dimensions':('lat',)}
    
    sStr['hist']={'data':hist,'long_name':'spatial_occurence_frequency','units':'%','_dimensions':('lat',)}
    sStr['amplitude']={'data':ampmn,'long_name':'amplitude','units':'cm','_dimensions':('lat',)}
    sStr['diameter']={'data':lenmn,'long_name':'diameter','units':'km','_dimensions':('lat',)}
    sStr['relvort']={'data':rvmn,'long_name':'relative_vorticity','units':'s-1','_dimensions':('lat',)}
    sStr['amplitude_rms']={'data':amprms,'long_name':'RMS_of_amplitude','units':'cm','_dimensions':('lat',)}
    sStr['diameter_rms']={'data':lenrms,'long_name':'RMS_of_diameter','units':'km','_dimensions':('lat',)}
    sStr['relvort_rms']={'data':rvrms,'long_name':'RMS_of_relative_vorticity','units':'s-1','_dimensions':('lat',)}
    
    tStr=OrderedDict()
    attrStr['description']=description+' - time'
    
    tStr['_dimensions']=tdimensions
    tStr['_attributes']=attrStr
    tStr['time']={'data':btime,'long_name':'time','units':'days since 1950-01-01 00:00:00UTC','_dimensions':('time',)}
    
    tStr['hist']={'data':thist,'long_name':'time_occurence_frequency','units':'%','_dimensions':('time',)}
    tStr['amplitude']={'data':tampmn,'long_name':'amplitude','units':'cm','_dimensions':('time',)}
    tStr['diameter']={'data':tlenmn,'long_name':'diameter','units':'km','_dimensions':('time',)}
    tStr['relvort']={'data':trvmn,'long_name':'relative_vorticity','units':'s-1','_dimensions':('time',)}
    tStr['amplitude_rms']={'data':tamprms,'long_name':'RMS_of_amplitude','units':'cm','_dimensions':('time',)}
    tStr['diameter_rms']={'data':tlenrms,'long_name':'RMS_of_diameter','units':'km','_dimensions':('time',)}
    tStr['relvort_rms']={'data':trvrms,'long_name':'RMS_of_relative_vorticity','units':'s-1','_dimensions':('time',)}
    
    #Setup file name
    sfname=os.path.dirname(filename)+os.path.sep+'space_clim.'+os.path.splitext(os.path.basename(filename))[0]+os.path.splitext(os.path.basename(filename))[1]
    tfname=os.path.dirname(filename)+os.path.sep+'time_clim.'+os.path.splitext(os.path.basename(filename))[0]+os.path.splitext(os.path.basename(filename))[1]
    
    #Write data
    res=obj.write(sStr,sfname,clobber=clobber)
    
    return res
