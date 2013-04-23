# -*- coding: utf-8 -*-
'''
    TEST_ANALYSIS
    @summary: This is a testing script, which applies the along-track wavelet transform <br />
              to a simulated red noise data set and shows the results.
    @note: The output graph should show two hovmöllers of the scale-averaged spectrum and<br />
           simulated sea level, with detected features marked as black circles.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since: Created in November 2012 by RD.
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
import kernel as ke
import altimetry.tools as AT
import matplotlib.pyplot as plt


if __name__ == "__main__" :
    
    
    #Simulate sea level data from red noise
    #######################################
    
    #Position data set in space and time
    lat = np.arange(43.0,44.0,0.01)
    lon = np.arange(6,6.5,0.005)
    dst=AT.calcul_distance(lat,lon)
    dx=np.median(AT.deriv(dst))
    N=len(lon)
    nt=25
    dt=9.9
    time=22705.0 + np.arange(0,nt)*dt
    dj_factor=20 #Scale factor wrt. number of elements
    
    #Red noise generation (lengthscale > 10 km)
    #sla=np.cumsum(np.random.randn(N*nt)).reshape((nt,N))
    sla=np.ma.array(np.cumsum(np.cumsum(np.random.randn(nt,N),axis=1),axis=0),mask=np.zeros((nt,N),dtype=bool))
    
    #Filter small scales
    for i in np.arange(nt):
        sla[i,:]=AT.loess(sla[i,:], dst, 10.)
    
    #Run wavelet analysis
    #####################
    
    #WV analysis
    sa_spectrum, sa_lscales, wvsla, daughter, params = ke.runAnalysis(lon,lat,time,sla,len_range=[10.,150.],m=0,dj_factor=20)
    
    #Detection of energy peaks on scale-averaged spectrum
    res = ke.detection(sa_spectrum,sa_lscales,params,amplim=1.0,clean=True)#np.ones((5,5),dtype=bool))
    
    #Plot results
    plt.subplot(2,1,1); plt.pcolormesh(dst,time,sa_spectrum); plt.colorbar(); plt.plot(dst[res[0]]+dx/2,time[res[1]]+dt/2,'ok'); plt.subplot(2,1,2); plt.pcolormesh(dst,time,sla - np.repeat(sla.mean(axis=1),N).reshape((nt,N))); plt.colorbar(); plt.plot(dst[res[0]]+dx/2,time[res[1]]+dt/2,'ok'); plt.show()
    
    print 'done'