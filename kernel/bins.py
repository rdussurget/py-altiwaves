'''
kernel.bins module
@summary: binning module
@author: Renaud DUSSURGET, LER/PAC IFREMER.
@since: Create in November 2012 by RD.
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
import numpy as np
from altimetry.externals.esutils_stat import histogram
import altimetry.tools as AT
if __debug__ : import matplotlib.pyplot as plt

def grid_time(time,tid,binsize=1):
    '''
    GRID_TIME
    @summary: Regrid the time index regularly.
    @param time {float}: time array (in julian days)
    @param tid {int}: time index
    @return hist, ind, btime : hist is the resulting histogram <br />
                               ind is the index within the histogram <br />
                               btime is the regridded time series.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @change: Create in November 2012 by RD.
    '''
    
    nt=len(time)
    hist,R = histogram(tid,binsize=binsize, min=0, max=nt-1, rev=True, use_weave=False, verbose=0)
    ind = AT.histogram_indices(hist, R)
    
    hist=np.ma.masked_array(hist,mask=hist==0,dtype=np.float32)
    hist.data[hist.mask]=hist.fill_value
    
    btime = np.ma.masked_array(np.zeros(len(hist)),mask=hist.mask,dtype=time.dtype)
#    btime=time.copy()
#    btime.data[:]=btime.fill_value
    for i in np.arange(len(hist)) : btime[i]=np.mean(time[tid[ind[i]]])
#    for i in np.arange(len(hist)) : blon[i]=np.mean(lon[eid[ind[i]]])
    
    
    
#    #normalise histogram with binsize
#    hist/=np.float(binsize)
    
    return hist, ind, btime
    

def grid_space(lon,lat,eid,binsize=7):
    '''
    GRID_SPACE
    @summary: Regrid the space index regularly.
    @param lon {float}: Longitude array
    @param lat {float}: Latitude array 
    @param eid {int}: eddy detection index along reference track
    @return hist, ind, blon, blat : hist is the resulting histogram <br />
                               ind is the index within the histogram <br />
                               blon, blat are the regridded lon/lat series.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @change: Create in November 2012 by RD.
    '''
    
    #Get lon/lat grid
    nx=len(lon)
#    id=np.arange(nx)
    hist,R = histogram(eid,binsize=binsize, rev=True, use_weave=False, verbose=0)
    ind = AT.histogram_indices(hist, R)
    
    hist=np.ma.masked_array(hist,mask=hist==0,dtype=np.float32)
    hist.data[hist.mask]=hist.fill_value
    
    blon = np.ma.masked_array(np.zeros(len(hist)),mask=hist.mask,dtype=lon.dtype)
    blat = np.ma.masked_array(np.zeros(len(hist)),mask=hist.mask,dtype=lon.dtype)
    for i,j in enumerate(hist.mask) :
        if ~j : blon[i]=np.mean([np.max(lon[eid[ind[i]]]),np.min(lon[eid[ind[i]]])])
    for i,j in enumerate(hist.mask) :
        if ~j: blat[i]=np.mean([np.max(lat[eid[ind[i]]]),np.min(lat[eid[ind[i]]])])
    
#    #normalise histogram with binsize
#    hist/=np.float(binsize)
    
    return hist, ind, blon, blat

def grid_var(var,hist,ind,method='mean'):
    '''
    GRID_VAR
    @summary: Compute statistics of a given variable in the gridded (space/time) frame.
    @param var : Variable to analyse
    @param hist {float}: histogram, as return by grid_time and/or grid_space 
    @param ind {int}: index within the histogram
    @return mn, rms : mn is the mean value of the binned variable, and rms its RMS.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @change: Create in November 2012 by RD.
    '''
    
    nx=len(hist)
    mn = np.ma.masked_array(np.zeros(nx),mask=np.ones(nx,dtype=bool))
    md = np.ma.masked_array(np.zeros(nx),mask=np.ones(nx,dtype=bool))
    rms = np.ma.masked_array(np.zeros(nx),mask=np.ones(nx,dtype=bool))
    mn.data[mn.mask]=mn.fill_value
    md.data[md.mask]=md.fill_value
    rms.data[rms.mask]=rms.fill_value
    for i,j in enumerate(np.array(ind)) :
        if ~hist.mask[i] :
            mn[i]=np.mean(var[j])
            md[i]=np.median(var[j])
#            if method == 'mean' : mn[i]=np.mean(var[j])
#            elif method == 'median' : mn[i]=np.median(var[j])
#            else : raise Exception('Unknown averaging method')
            rms[i]=AT.rms(var[j])
    
    if method == 'mean' : return mn, rms 
    elif method == 'median' : return md, rms
    else : raise Exception('Unknown averaging method')

def bin_space(lon,lat,eind,amplitude,diameter,relvort,binsize=7,method='mean',verbose=1):
    
    if verbose >= 1:
        str_header = '\t===Time binning parameters===\n\tbinsize:{0}pts, method:{1} '.format(np.int(binsize),method)    
        print(str_header)
    
    hist, ind, blon, blat = grid_space(lon,lat,eind[0],binsize=binsize)
    ampmn, amprms= grid_var(amplitude,hist,ind,method=method)
    lenmn, lenrms= grid_var(diameter,hist,ind,method=method)
#    wvlenmn, wvlenrms= grid_var(wvdiameter,hist,ind,method=method)
#    uglmn, uglrms= grid_var(ugdiameter,hist,ind,method=method)
    rvmn, rvrms= grid_var(relvort,hist,ind,method=method)
    
    return blon, blat, hist,  ampmn, lenmn, rvmn, amprms, lenrms, rvrms

def bin_time(time,eind,amplitude,diameter,relvort,binsize=1,method='mean',verbose=1):
    
    if verbose >= 1:
        str_header = '\t===Spatial binning parameters===\n\tbinsize:{0}pts, method:{1} '.format(np.int(binsize),method)    
        print(str_header)
    
    date,datetime=AT.cnes_convert(time)
    thist,tind,btime=grid_time(time,eind[1],binsize=binsize)
    trvmn,trvrms=grid_var(relvort,thist,tind,method=method)
    tampmn,tamprms=grid_var(amplitude,thist,tind,method=method)
    tlenmn,tlenrms=grid_var(diameter,thist,tind,method=method)
#    twvlmn,twvlrms=grid_var(wvdiameter,thist,tind,method=method)
#    tuglmn,tuglrms=grid_var(ugdiameter,thist,tind,method=method)
    return datetime, btime, thist, tampmn, tlenmn, trvmn, tamprms, tlenrms, trvrms