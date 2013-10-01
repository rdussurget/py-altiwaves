# -*- coding: utf-8 -*-
'''
    TEST_SPECTRAL_ANALYSIS
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

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


if __name__ == "__main__" :
    
    
    #Simulate sea level data from red noise
    #######################################
    
    #Position data set in space and time
    lonlatres=0.001 #Increase this factor if you want faster computation time (lesser resolution)
    nt=50   #Time steps : decrease it to have faster computations
    dj_factor=20 #Scale factor wrt. number of elements :  increase it to get faster computation
    
    lat = np.arange(43.0,44.0,lonlatres*2)
    lon = np.arange(6,6.5,lonlatres)
    dst=AT.calcul_distance(lat,lon)
    dx=np.median(AT.deriv(dst))
    N=len(lon)
    dt=9.9
    time=22705.0 + np.arange(0,nt)*dt
    
    
    #Red noise generation (lengthscale > 10 km)
    #sla=np.cumsum(np.random.randn(N*nt)).reshape((nt,N))
    sla=np.ma.array(np.cumsum(np.cumsum(np.random.randn(nt,N),axis=1),axis=0),mask=np.zeros((nt,N),dtype=bool))
    
    #Filter small scales
    for i in np.arange(nt):
        sla[i,:]=AT.loess(sla[i,:], dst, 10.)
    
    #Run wavelet analysis
    #####################
    
    #WV periodogram analysis   
    perWV = ke.periodogram_analysis(dx, sla.transpose(),res_factor=dj_factor,average=True)
    perpsdWV = perWV['psd']
    pWV = perWV['p']
    DWV = perWV['D']
    
    gx,gy=np.meshgrid(DWV,pWV)

    #3D periodogram
    dum=np.ma.array(np.log10(perpsdWV),mask=(np.log10(perpsdWV) < -6))
    dum.data[dum.mask]=-2
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(68, -30)
    surf = ax.plot_surface(gx, np.log(gy), dum,vmin=-6,vmax=3,cmap=cm.jet,linewidth=0, antialiased=True,shade=True)#, rstride=1, cstride=1, cmap=cm.jet, extend3d=True)
    ax.set_zlim([-6,3])
    ax.set_xlabel('Along-Track Distance(km)')
    ax.set_ylabel('log10(Spatial scale - m)')
    ax.set_zlabel('log10(Power Spectral Density - cm2)')
    plt.show()
    
    print 'done'