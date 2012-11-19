import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import alti_tools as AT
import scipy.interpolate

import matplotlib.pyplot as plt

'''
Created on 9 nov. 2012

@author: rdussurg
'''
def _2D(sa_spectrum, sa_lscales, win_width=5., amplim=3., len_range=[60.,450.]):
    
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
    xs=1 #size wrt center
    ys=2
    kd=np.zeros((2*ys+1,2*xs+1),dtype=bool)
    kd[:,xs]=True
    kd[ys,:]=True #This is a cross-shaped kernel

#    anisotropy=(3,1)
#    kx,ky= np.mgrid[-xs:xs+1, -ys:ys+1]
#    kx*=anisotropy[0]
#    ky*=anisotropy[1]
#    kd = np.exp(-(kx**2+ky**2)).transpose() > 0.05 #Kernel valid for distances < 95% of normal distribution
    
    data_max = filters.maximum_filter(sas, footprint=kd)
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
  
  
  
#  #Detect peaks on vectorised form of spectrum
#  detx=lclxtrem(REFORM(spec_nonan,npts*nt),win_width,/MAXIMA)
#  dety=lclxtrem(REFORM(TRANSPOSE(spec_nonan),npts*nt),win_width,/MAXIMA)
#  
#  ;Convert from vector to matrix indices
#  detxa=ARRAY_INDICES(spec_nonan,detx)
#  detya=SHIFT(ARRAY_INDICES(TRANSPOSE(spec_nonan),dety),1,0) ;Due to the transposition, x will be in col 1 and y in col 0 -> shift columns
#  
#  ;A is the shortest list & B the longest
#  mndet=MIN([N_ELEMENTS(detx),N_ELEMENTS(dety)],mnid)
#  a=(mnid EQ 0)? TEMPORARY(detxa) : TEMPORARY(detya)
#  b=(mnid EQ 0)? TEMPORARY(detya) : TEMPORARY(detxa)
#  
#  na=N_ELEMENTS(a)/2
#  nb=N_ELEMENTS(b)/2
#  
#  tot=INTARR(nb)
#  nz=0
#  FOR i=0, na -1 DO BEGIN
#    tot[*]=0
#    at=cmreplicate(a[*,i],nb)
#    tot=TOTAL(at EQ b,1)
#    dumw=WHERE(tot EQ 2,ntot)
#    IF (ntot GT 0) THEN BEGIN
#      z=(nz EQ 0)? dumw : [z,dumw]
#      nz+=ntot
#    ENDIF
#  ENDFOR
#
