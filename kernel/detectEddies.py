# -*- coding: utf-8 -*-
import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import scipy.interpolate
import matplotlib.pyplot as plt
from kernel.getScales import cyclone, decorrelation_scale, solid_body_scale,\
    eddy_amplitude

'''
Created on 9 nov. 2012

@author: rdussurg
'''
def _2D(sa_spectrum, amplim=3., kernel=None): #len_range=[60.,450.], 
    '''
    _2D
    @summary: Eddy detection on both axes of the scale-averaged spectrum.
    @note: This technique was first applied in Le Hénaff et al., 2012. Cyclonic<br />
           activity in the eastern Gulf of Mexico: characterization. Submitted to <br />
           Progress in Oceanography.
    @param sa_spectrum: scale-avergaed spectrum returned by runAnalysis().
    @param sa_lscales: Lengthscale (in km) of the most energetic wavelet returned<br />
           by the wavelet analysis.
    @keyword amplim: Amplitude threshold for eddy detection.
    @keyword kernel: Kernel to pass to maximum_filter1d(). Default is a cross-shaped<br / >
                     kernel.
    @return x, y: Detection locations along X and Y axis of the SA spectrum.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since : November 2012.
    @change: Create in November 2012 by RD.
    '''
    
#  IF (~exist(lclxtrem)) THEN lclxtrem=1B
    sas=sa_spectrum.copy()
    shape=sas.shape
    nt=shape[0]
    npts=shape[1]
    
    #Interpolate on missing points (this should be to avoid peaking anomalies on edges of data gaps).
    if isinstance(sas,np.ma.masked_array) :
        xx=np.arange(npts)
        yy=np.arange(nt)
        xout,yout=np.meshgrid(xx, yy)
        points=zip(*(xout[~sa_spectrum.mask].flatten(), yout[~sa_spectrum.mask].flatten()))
        values=sas.data[~sa_spectrum.mask].flatten()
        xi=zip(*(xout[sa_spectrum.mask].flatten(), yout[sa_spectrum.mask].flatten()))
        sas.data[sa_spectrum.mask]=scipy.interpolate.griddata(points, values, xi, method='linear') #Do not use nearest neighbour with maximum_filter
    
    #define maximum filter kernel
    if kernel is None :
        xs=1 #size wrt center
        ys=2
        kernel=np.zeros((2*ys+1,2*xs+1),dtype=bool)
        kernel[:,xs]=True
        kernel[ys,:]=True #This is a cross-shaped kernel
    
    #    anisotropy=(3,1)
    #    kx,ky= np.mgrid[-xs:xs+1, -ys:ys+1]
    #    kx*=anisotropy[0]
    #    ky*=anisotropy[1]
    #    kernel = np.exp(-(kx**2+ky**2)).transpose() > 0.05 #Kernel valid for distances < 95% of normal distribution
    
    data_max = filters.maximum_filter(sas, footprint=kernel)
    maxima = (sas == data_max) & (data_max > amplim**2)
  
#  fg = np.isfinite(sa_spectrum) & (np.sqrt(sa_spectrum >= amplim)) & (sa_lscales < np.max(len_range)) & (sa_lscales >= np.min(len_range)) 
  
    #Label & enumerate objects
    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
  
    #Get X,Y positions
    x, y = [], []
    for dy,dx in slices:
        nx=(dx.stop - dx.start)
        x_center = (dx.start + dx.stop - 1)/2
        if nx == 1 : x.append(x_center)
        ny=(dy.stop - dy.start)
        y_center = (dy.start + dy.stop - 1)/2    
        if ny == 1 : y.append(y_center)
  
    #Check if peaks are found in unmasked data
    if isinstance(sas,np.ma.masked_array):
        inter=np.array(list(set(zip(*(x,y))).difference(set(zip(*(xout[sa_spectrum.mask].flatten(),yout[sa_spectrum.mask].flatten()))))))
        x=inter[:,0]
        y=inter[:,1]

    return x, y
  
def _1D(sa_spectrum, sa_lscales, win_width=5., amplim=3., len_range=[60.,450.]):
    '''
    _1D
    @summary: Detection of the most energetic eddy along the time axis of the <br />
              scale-averaged spectrum.
    @note: This is the original technique applied in :
           Dussurget, R, F Birol, R.A. Morrow, et P. De Mey. 2011. « Fine Resolution<br />
           Altimetry Data for a Regional Application in the Bay of Biscay ». Marine<br />
           Geodesy 2 (34): 1‑30. doi:10.1080/01490419.2011.584835.
    @warning: This function is currently deprecated. Use _2D instead.
    @param sa_spectrum: scale-avergaed spectrum returned by runAnalysis().
    @param sa_lscales: Lengthscale (in km) of the most energetic wavelet returned<br />
           by the wavelet analysis.
    @keyword amplim: Amplitude threshold for eddy detection.
    @keyword win_width: Window size of the maximum filter.
    @keyword len_range: Range of admitted lengthscales (km). 
    @return x, y: Detection locations along X and Y axis of the SA spectrum.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since : November 2012.
    @change: Create in November 2012 by RD.
    '''
    raise Exception("[ERROR] This function is not available yet and/or deprecated.")
    return

def detection(sa_spectrum,amplim=0.03,twoD=True):
    eind = _2D(sa_spectrum, amplim=amplim) if twoD else _1D(sa_spectrum, amplim=amplim)
    eind = np.squeeze(eind)
    return eind
    