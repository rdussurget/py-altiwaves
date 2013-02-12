# -*- coding: utf-8 -*-*
'''
detectEddies
@summary: functions for detection of eddy-like features after wavelet analysis  
@author: Renaud DUSSURGET, LER/PAC IFREMER.
@since: Created on 9 nov. 2012
'''
import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import scipy.interpolate
import matplotlib.pyplot as plt
from kernel.getScales import cyclone, decorrelation_scale, solid_body_scale,\
    eddy_amplitude

def _2D(sa_spectrum, amplim=0.04, kernel=None, verbose=1): #len_range=[60.,450.], 
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
    
    #Print detection parameters
    ###########################   
    #Interpolate on missing points (this should be to avoid peaking anomalies on edges of data gaps).
    if isinstance(sas,np.ma.masked_array) :
        xx=np.arange(npts)
        yy=np.arange(nt)
        xout,yout=np.meshgrid(xx, yy)
        points=zip(*(xout[~sa_spectrum.mask].flatten(), yout[~sa_spectrum.mask].flatten()))
        values=sas.data[~sa_spectrum.mask].flatten()
        xi=zip(*(xout[sa_spectrum.mask].flatten(), yout[sa_spectrum.mask].flatten()))
        sas[sa_spectrum.mask]=scipy.interpolate.griddata(points, values, xi, method='linear',fill_value=sas.fill_value) #Do not use nearest neighbour with maximum_filter
        sas.mask[sas.data == sas.fill_value]=True
    
    #define maximum filter kernel
    if kernel is None :
        xs=2 #size wrt center
        ys=2
        kernel=np.zeros((2*ys+1,2*xs+1),dtype=bool)
        kernel[:,xs]=True
        kernel[:,xs-1]=True
        kernel[:,xs+1]=True
        kernel[ys,:]=True
##        kernel[:,xs]=True
#        kernel[ys,:]=True #This is a cross-shaped kernel
#        kernel=~kernel
#        kernel[1:4,:]=True
    
    #    anisotropy=(3,1)
    #    kx,ky= np.mgrid[-xs:xs+1, -ys:ys+1]
    #    kx*=anisotropy[0]
    #    ky*=anisotropy[1]
    #    kernel = np.exp(-(kx**2+ky**2)).transpose() > 0.05 #Kernel valid for distances < 95% of normal distribution
    
    if verbose >= 1 : print('\tkernel shape :\t{0}'.format(kernel.shape))
    if verbose > 1 :
        for i in np.arange(kernel.shape[0]) : print ('\tkernel:\t\t'+str(kernel[i]) if i == 0 else '\t\t\t'+str(kernel[i]))
    
    data_max = filters.maximum_filter(sas, footprint=kernel)
    maxima = (sas == data_max) & (data_max > amplim**2) & (~sa_spectrum.mask)
  
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

def clean_indices(sa_spectrum,sa_lscales,eind,params):
    """
    @summary: Removes not valid indices (masked point within 4 points or other point within lenghtscale)
    @author: Renaud DUSSURGET
    """
    xind=eind[0,:]
    yind=eind[1,:]
    toRm=np.array([])
    
    timescale_width = np.sqrt(4 - np.log(2) - np.log(np.pi)) #THIS IS ONLY VALID FOR wavelet DOG, m=0
    spectral_width = (np.sqrt(4-0.5*np.log(np.pi)))
    N=params['N']
    dj=params['dj']
    s0=params['s0']
    l2s=params['len2scale']
    dt=params['dt']*1e-5
    
    #remove indices with masked data closer than 4 points
    for i,x in enumerate(xind):
        
        #Check structure lengthscale
        s=sa_lscales[yind[i],x]*l2s #spatial lengthscale
        j=np.log((s*1e5)/s0) / (dj*np.log(2)) #Get the wavenumber
        ets=(timescale_width*s)/dt  #This is the e-folding time of the wavelet function
                                    #as a function of scale
        e=np.ceil(1.25*ets) #Add a 25% confidence margin
#        esw=(spectral_width*N*dt)/(s) #This is the spectral uncertainty
#        s1=s0*2**((j-esw)*dj)*1e-5
#        s2=s0*2**((j+esw)*dj)*1e-5
#        e = 0.75*np.floor((s)**2/(np.sqrt(2)*params['N'])) #This is the e-folding time of the wavelet function
#                                                          #as a function of scale, divided by a margin of 75%

#        exp=np.exp
#        sqrt=np.sqrt
#        pi=np.pi
#        s=sa_lscales[6,23]*params['len2scale']
#        d=daughter[6,:]
#        f=np.sqrt(2*np.log(np.sqrt(2*np.pi)))
#        f=2.0*np.sqrt(np.log(np.sqrt(np.pi)))
#        e = (s/dx) * f
#
#        n=np.arange(N)
#        nu=(n*dx)/s
#        k=n
#        omega= k/(N*dx)
#        
#        def psi(nu):return (1.0/sqrt(2.0*pi))*exp(-0.5*(nu**2)) 
#        timescale_width = np.sqrt(4 - np.log(2) - np.log(np.pi))
#        spectral_width = (np.sqrt(4-0.5*np.log(np.pi)))
#        (timescale_width*s)/dx #temporal e-folding time
#        (spectral_width*N*dx)/(s) #spectral e-folding time
#        
#        psi_hat = (1.0/sqrt(sqrt(pi)))*exp(-0.5*((s*omega)**2))
#        def psi_hat(s_omega): return (1.0/sqrt(sqrt(pi)))*exp(-0.5*((s_omega)**2))
#        
#        (sqrt(2)*(s**2))/dx
#        
#        psi=(1.0/sqrt(2.0*pi))*exp(-0.5*(nu**2))

#        e = (s/dt) * f

        #Get surrounding energy peaks for same cycle
        pks=xind[(yind == yind[i]) & (xind != x)]
        pks=np.array(list(set(pks).difference(set(toRm))))
        
        dst=np.abs(pks - x)
        
        #For each more energetic peaks within the e-folding scale, remove current point, otherwise remove secondary peak
        for s in pks[dst < e]:
#            print yind[i], s, x
            if sa_spectrum[yind[i],s] > sa_spectrum[yind[i],x] : toRm=np.append(toRm,i)
            else : toRm=np.append(toRm,s)

        if np.min(np.abs(np.arange(sa_spectrum.shape[1])[sa_spectrum.mask[yind[i],:]] - x)) < 3 :
            toRm = np.append(toRm,i) #Append to index list if true
        
        
        
    eindin=eind.copy()
    eind=np.squeeze([np.delete(xind,toRm),np.delete(yind,toRm)])
#    plt.pcolormesh(sa_spectrum);plt.plot(eindin[0,:],eindin[1,:],'ok');plt.plot(eind[0,:],eind[1,:],'.r');plt.show()
    
    return eind

def detection(sa_spectrum,sa_lscales,params,amplim=0.03,twoD=True,clean=True, verbose=1, **kwargs):
    if verbose >= 1:
        str_header = '\t===Eddy detection parameters===\n\tthreshold:{0}cm, clean:{1}, 2D:{2} '.format(np.int(amplim*100), clean, twoD)    
        print(str_header)
    eind = _2D(sa_spectrum, amplim=amplim, verbose=verbose, **kwargs) if twoD else _1D(sa_spectrum, amplim=amplim, verbose=verbose, **kwargs)
    eind = np.squeeze(eind)
    n_noclean=eind.shape[1]
    if clean : eind = clean_indices(sa_spectrum,sa_lscales, eind,params)
    n_clean=eind.shape[1]
    if verbose >= 1: print '\tDone : {0} peaks found ({1} of {3} ({2}%) rejected)'.format(n_clean,n_noclean - n_clean,np.round(100*np.float(n_noclean - n_clean)/n_noclean).astype(int), n_noclean)
    return eind
    