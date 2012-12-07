import numpy as np
from altimetry.externals.esutils_stat import histogram
import altimetry.tools as AT

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
    
    btime=time.copy()
    btime[:]=0.
    for i in np.arange(nt) : btime[i]=np.mean(time[tid[ind[i]]])
    
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
    
    blon = np.repeat(np.NaN,len(hist))
    blat = np.repeat(np.NaN,len(hist))
    for i in np.arange(len(hist)) : blon[i]=np.mean(lon[eid[ind[i]]])
    for i in np.arange(len(hist)) : blat[i]=np.mean(lat[eid[ind[i]]])
    
    return hist, ind, blon, blat

def grid_var(var,hist,ind):
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
    mn = np.repeat(np.NaN,nx)
    rms = np.repeat(np.NaN,nx)
    for i in np.arange(nx) : mn[i]=np.mean(var[ind[i]])
    for i in np.arange(nx) : rms[i]=AT.rms(var[ind[i]])
    
    return mn, rms

def bin_space(lon,lat,eind,amplitude,diameter,relvort,ugdiameter,wvdiameter):
    hist, ind, blon, blat = grid_space(lon,lat,eind[0])
    ampmn, amprms= grid_var(amplitude,hist,ind)
    lenmn, lenrms= grid_var(diameter,hist,ind)
    wvlenmn, wvlenrms= grid_var(wvdiameter,hist,ind)
    uglmn, uglrms= grid_var(ugdiameter,hist,ind)
    rvmn, rvrms= grid_var(relvort,hist,ind)
    
    return blon, blat, hist,  ampmn, lenmn, rvmn, uglmn, wvlenmn, amprms, lenrms, rvrms, uglrms, wvlenrms

def bin_time(time,eind,amplitude,diameter,relvort,ugdiameter,wvdiameter):
    date,datetime=AT.cnes_convert(time)
    thist,tind,btime=grid_time(time,eind[1])
    trvmn,trvrms=grid_var(relvort,thist,tind)
    tampmn,tamprms=grid_var(amplitude,thist,tind)
    tlenmn,tlenrms=grid_var(diameter,thist,tind)
    twvlmn,twvlrms=grid_var(wvdiameter,thist,tind)
    tuglmn,tuglrms=grid_var(ugdiameter,thist,tind)
    return datetime, btime, thist, tampmn, tlenmn, trvmn, tuglmn, twvlmn, tamprms, tlenrms, trvrms, tuglrms, twvlrms