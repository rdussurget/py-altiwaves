# -*- coding: latin-1 -*-
'''
IO
@summary : I/O tools for wavelet analysis
@requires: altimetry.tools.nctools
Created on 6 déc. 2012
@author: rdussurg
'''
import os
import numpy as np
from collections import OrderedDict
from altimetry.tools.nctools import nc

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
                  verbose=1):
    
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
            attrStr.update({a:attrStr[a]})
    
    #Setup output structure
    dataStr=OrderedDict()
    dataStr['_dimensions']=dimensions
    dataStr['_attributes']=attrStr
    dataStr['lon']={'data':lon,'long_name':'longitude','units':'degrees_north','_dimensions':('record',)}
    dataStr['lat']={'data':lat,'long_name':'latitude','units':'degrees_east','_dimensions':('record',)}
    dataStr['time']={'data':time,'long_name':'time of measurement','units':'days since 1950-01-01 00:00:00UTC','_dimensions':('cycle',)}
    dataStr['cycle']={'data':cycle,'long_name':'cycle number','units':'1','_dimensions':('cycle',)}
    dataStr['track']={'data':track,'long_name':'track number','units':'1','_dimensions':('track',)}
    dataStr['sla']={'data':sla,'long_name':'sea_level_anomaly','units':'m','_dimensions':('record','cycle')}
    dataStr['sa_spectrum']={'data':sa_spectrum.data,'long_name':'scale_averaged_spectrum','units':'m^^2','_dimensions':('record','cycle')}
    dataStr['sa_lscales']={'data':sa_lscales.data,'long_name':'spatial_lengthscale','units':'km','_dimensions':('record','cycle')}
    dataStr['wvsla']={'data':wvsla.data,'long_name':'filtered_sea_level_anomaly','units':'m','_dimensions':('record','cycle')}
    dataStr['daughter']={'data':daughter.data,'long_name':'most_energetic_wavelet','units':'m','_dimensions':('record','cycle')}
    
    #Write data
    res=obj.write(dataStr,filename)
    
    return res

def save_detection(filename,
                   eind,
                   lon,
                   lat,
                   amplitude,
                   diameter,
                   relvort,
                   ugdiameter,
                   wvdiameter,
                   dimensions=None,
                   attributes=None,
                   verbose=1):
    
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
            attrStr.update({a:attrStr[a]})
    
    #Setup output structure
    dataStr=OrderedDict()
    dataStr['_dimensions']=dimensions
    dataStr['_attributes']=attrStr
    dataStr['xind']={'data':eind[0,:],'long_name':'along_track_index','units':'1','_dimensions':('record',)}
    dataStr['yind']={'data':eind[1,:],'long_name':'time_index','units':'1','_dimensions':('record',)}
    dataStr['lon']={'data':lon,'long_name':'longitude','units':'degrees_north','_dimensions':('record',)}
    dataStr['lat']={'data':lat,'long_name':'latitude','units':'degrees_east','_dimensions':('record',)}
    dataStr['amplitude']={'data':amplitude,'long_name':'amplitude','units':'m','_dimensions':('record',)}
    dataStr['diameter']={'data':diameter,'long_name':'diameter','units':'km','_dimensions':('record',)}
    dataStr['relvort']={'data':relvort,'long_name':'relative_vorticity','units':'s-1','_dimensions':('record',)}
    dataStr['ugdiameter']={'data':ugdiameter,'long_name':'eddy_core_diameter','units':'km','_dimensions':('record',)}
    dataStr['wvdiameter']={'data':wvdiameter,'long_name':'wavelet_diameter','units':'km','_dimensions':('record',)}
    
    #Write data
    res=obj.write(dataStr,filename)
    
    return res

def save_binning(filename,
                 blon, blat, hist,  ampmn, lenmn, rvmn, uglmn, wvlenmn, amprms, lenrms, rvrms, uglrms, wvlenrms,
                 btime, thist, tampmn, tlenmn, trvmn, tuglmn, twvlmn, tamprms, tlenrms, trvrms, tuglrms, twvlrms,
                 dimensions=None,
                 attributes=None,
                 verbose=1):
    
    
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
    attrStr['description']='Climatology of eddy-like processes spatial variability'
    if attributes is not None : #Append attributes if asked
        for a in [k for k in attributes.keys() if (k is not 'description') and (k is not 'input')] :
            attrStr.update({a:attrStr[a]})
    
    #Setup output structure
    sStr=OrderedDict()
    sStr['_dimensions']=sdimensions
    sStr['_attributes']=attrStr
    sStr['lat']={'data':blat,'long_name':'latitude','units':'degrees_east','_dimensions':('lat',)}
    sStr['lon']={'data':blon,'long_name':'longitude','units':'degrees_north','_dimensions':('lat',)}
    
    sStr['hist']={'data':hist,'long_name':'spatial_occurence_frequency','units':'%','_dimensions':('lat',)}
    sStr['amplitude']={'data':ampmn,'long_name':'amplitude','units':'m','_dimensions':('lat',)}
    sStr['diameter']={'data':lenmn,'long_name':'diameter','units':'km','_dimensions':('lat',)}
    sStr['relvort']={'data':rvmn,'long_name':'relative_vorticity','units':'s-1','_dimensions':('lat',)}
    sStr['ugdiameter']={'data':uglmn,'long_name':'eddy_core_diameter','units':'km','_dimensions':('lat',)}
    sStr['wvdiameter']={'data':wvlenmn,'long_name':'wavelet_diameter','units':'km','_dimensions':('lat',)}
    sStr['amplitude_rms']={'data':amprms,'long_name':'RMS_of_amplitude','units':'m','_dimensions':('lat',)}
    sStr['diameter_rms']={'data':lenrms,'long_name':'RMS_of_diameter','units':'km','_dimensions':('lat',)}
    sStr['relvort_rms']={'data':rvrms,'long_name':'RMS_of_relative_vorticity','units':'s-1','_dimensions':('lat',)}
    sStr['ugdiameter_rms']={'data':uglrms,'long_name':'RMS_of_eddy_core_diameter','units':'km','_dimensions':('lat',)}
    sStr['wvdiameter_rms']={'data':wvlenrms,'long_name':'RMS_of_wavelet_diameter','units':'km','_dimensions':('lat',)}
    
    tStr=OrderedDict()
    attrStr['description']='Climatology of eddy-like processes time variability'
    
    tStr['_dimensions']=tdimensions
    tStr['_attributes']=attrStr
    tStr['time']={'data':btime,'long_name':'time','units':'days since 1950-01-01 00:00:00UTC','_dimensions':('time',)}
    
    tStr['hist']={'data':thist,'long_name':'time_occurence_frequency','units':'%','_dimensions':('time',)}
    tStr['amplitude']={'data':tampmn,'long_name':'amplitude','units':'m','_dimensions':('time',)}
    tStr['diameter']={'data':tlenmn,'long_name':'diameter','units':'km','_dimensions':('time',)}
    tStr['relvort']={'data':trvmn,'long_name':'relative_vorticity','units':'s-1','_dimensions':('time',)}
    tStr['ugdiameter']={'data':tuglmn,'long_name':'eddy_core_diameter','units':'km','_dimensions':('time',)}
    tStr['wvdiameter']={'data':twvlmn,'long_name':'wavelet_diameter','units':'km','_dimensions':('time',)}
    tStr['amplitude_rms']={'data':tamprms,'long_name':'RMS_of_amplitude','units':'m','_dimensions':('time',)}
    tStr['diameter_rms']={'data':tlenrms,'long_name':'RMS_of_diameter','units':'km','_dimensions':('time',)}
    tStr['relvort_rms']={'data':trvrms,'long_name':'RMS_of_relative_vorticity','units':'s-1','_dimensions':('time',)}
    tStr['ugdiameter_rms']={'data':tuglrms,'long_name':'RMS_of_eddy_core_diameter','units':'km','_dimensions':('time',)}
    tStr['wvdiameter_rms']={'data':twvlrms,'long_name':'RMS_of_wavelet_diameter','units':'km','_dimensions':('time',)}
    
    #Setup file name
    sfname=os.path.dirname(filename)+'\\space_clim.'+os.path.splitext(os.path.basename(filename))[0]+os.path.splitext(os.path.basename(filename))[1]
    tfname=os.path.dirname(filename)+'\\time_clim.'+os.path.splitext(os.path.basename(filename))[0]+os.path.splitext(os.path.basename(filename))[1]
    
    #Write data
    res=obj.write(sStr,sfname)
    res=obj.write(tStr,tfname)
    
    return res

