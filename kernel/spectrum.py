# -*- coding: utf-8 -*-
'''
kernel.spectrum module
@summary: spectral analysis functions
@since: Created on 10 janv. 2013
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

import numpy as np
from altimetry.tools import deriv, interp1d, detrend as detrend_fun
import kernel

if __debug__ : import matplotlib.pyplot as plt

#define constants
km2cm = 1e5
cm2km = km2cm ** (-1.0)

def get_spec(dx,Vin,verbose=False,m=6,gain=1.0,res_factor=10.,integration=False,periodogram=False,mother='morlet'):
    
    
    #Setup Wavelet Transform parameters
    ###################################
    
    if mother == 'morlet':
        scale2len = (4*np.pi) / (m + np.sqrt(2+m*m)) #"fourrier factor" (for Morlet)
        len2scale = 1/scale2len
    elif mother=='dog':
        len2scale = np.sqrt(m+0.5) / (2*np.pi)
        scale2len = (2*np.pi) / np.sqrt(m+0.5) #"fourrier factor"
    else :
        raise Exception('Wavelet function "{0}" is not defined - choose between "morlet" & "dog"')
    
    N=len(Vin)
    
    T = dx*N
    s0 = 2.0 * dx if mother == 'morlet' else (2*dx) * len2scale #smallest wavescale (devided by fourier wavelength)
    dj = res_factor * (dx/T) * len2scale #set scale interval (remember this is scaled by s0^(-1/2) ; no units!
    dj0 =  (dx/T) * len2scale
    J = np.fix((np.log((T*len2scale)/s0) / np.log(2.)) / dj).astype(int) #number of scales from s0 to T
    
    #Setup wavelet object
    exec('mother_obj=kernel.wavelet.{0}({1})'.format(mother,m))
        
    #Run transform
    wave, s, k, coi, daughter, fft, fftfreqs = kernel.wavelet.cwt(Vin, dx, dj, s0, J,mother_obj)
    
    Cd = mother_obj.cdelta
    p = 1. / k
    
    #Mask data
    ##########
    # 1) Out of confidence interval
    mask = np.repeat(coi,J+1).reshape((N,J+1)).transpose() <= np.repeat(p,N).reshape((J+1,N))
#    mask = np.ones(N,J+1)
    
    a = np.ma.array(np.real(wave),mask=mask)             #Real part of the transform
    b = np.ma.array(np.imag(wave),mask=mask)
    csquared = a**2.0 + b**2.0
    
    c = np.sqrt(csquared)
    
    if not periodogram : csquared = csquared.sum(axis=1) / (~mask).sum(axis=1)
    
    # integration of the spectrum
    # shift to have values centered on the intervals
    dk = deriv(k)
    k_ = interp1d(np.arange(J+1),k,np.arange(J)+0.5)  #Shift wavelengths by half the unit
    
    #Spectral integration
    if integration and not periodogram:
        esd = k_*0.0
        psd = k_*0.0
        dk = deriv(k)
        for i in np.arange(len(k_)):
            esd[i] = np.sum((csquared * (N/2.0))[(k > (k_[i]-dk[i])) & (k < (k_[i]+dk[i]))]) / 2.0
            psd[i] = np.sum(csquared[(k > (k_[i]-dk[i])) & (k < (k_[i]+dk[i]))]) / (2.0**2) #This is variance units integration
        fq=k_
    else :
        esd = csquared
        psd = esd.copy()  /2.
        fq=k.copy()
        
    psd = (psd / (dj*res_factor))
    
    #Get frequencies and period  
    p=1/fq
    
    return {'psd':psd,'esd':esd,'fq':fq,'p':p,'gain':gain}

def spectral_analysis(dx,Ain,res_factor=10.0,tapering=None,overlap=None,wsize=None,alpha=3.0,detrend=False,normalise=False,integration=False):
    """
     Spectral_Analysis
     @summary: This function performs a spatial spectral analysis with different options on a time series of SLA profiles.
     @param dx {type:numeric} : sampling distance
     @param Ain {type:numeric} : 2D table of sla data with time along 2nd axis (NXxNT with NX the spatial length and NT the time length)
     @keyword tapering {type:string|bool|nd.array} : apply tapering to the data. <br \>
                    If this keyword is of type bool : apply hamming window. <br \>
                    If this keyword is a string : apply a hamming ('hamm'), hann ('hann'), kaiser-bessel ('kaiser'), kaiser-bessel ('blackman') or no ('none') tapering function. <br \>
                    If this keyword is an nd.array aobject : apply this array as taper.
     @keyword overlap {type:float} : overlap coefficient of the windows (0.75 means 75% overlap).
     @keyword wsize {type:numeric} : size of the sub-segments.
     @keyword normalise {type:bool,default:False} : If True, normalise the spectrum by its overall energy content.
     @keyword detrend {type:bool,default:False} : If True, removes a linear trend to the segmented signal (if tapered) or to the whole signal (if not tapered).
     @keyword integration {type:bool,default:False} : If True, integrate the spectrum between 2 frequencies. 
     @param sla {type:numeric} : data
     @return: a spectrum structrue with Energy Spectral Density ('esd'), Power Spectral Density ('PSD'), frequency ('fq'), wavelength ('p') and tapering parameters.
    
     @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
     @change: Created by RD, December 2012
    """
    
    A=Ain.copy()

    #Check dimensions
    sh = A.shape
    ndims = len(sh)
    N = sh[0] #Time series are found along the last dimension
    
    #If vector, add one dimension
    if ndims == 1 :
        A = A.reshape((N,1))
        sh = A.shape
        ndims = len(sh)
    
    nr = sh[1] #Number of repeats  
    nt = nr
    
#    gain=1.0 #Scaling gain... (used for tapering issues)
    
#    #Get the overall energy
#    spec=get_spec(dx, A[:,0])
#    F=spec['fq']   
#    Eref = ((A[:,0]-A[:,0].mean())**2).sum() #Get the reference energy level
#    ESDref=spec['esd']
#    SFactor=Eref/spec['esd'].sum()
#    ESDref*=SFactor
#    PSDref=spec['psd']*SFactor
#    print 'Check parseval theorem : SUM|Y(f)|�={0}, SUM|y(t)|�={1}'.format(spec['esd'].sum(),((A[:,0]-A[:,0].mean())**2).sum())
    
    #Apply tapering if asked
    ########################
    if tapering is not None:
        
        #Set tapering defaults
        overlap=0.50 if overlap is None else overlap
        wsize=0.5*N if wsize is None else wsize

        #Get time splitting (tapering) parameters
        #########################################
        a = np.float32(wsize)
        b = np.float32(overlap) 
        c = np.float32(N) 
        nn=np.floor((c - (a * b))/(a - (a * b))) #This is the number of segments
        print 'Number of windows :{0}\nTotal windowed points : {1} ({2} missing)\nTotal points : {3}'.format(nn,nn*wsize,N - nn*wsize,N)
        
        ix = np.arange(nn) * ((1.0 - b) * a) #These are the starting points of each segments

        #Moving window
        ##############
        dum = np.zeros((wsize, nn, nr),dtype=np.float64)
        for j in np.arange(nr):
            for i in np.arange(nn): #looping through time to get splitted time series 
                dum[:,i,j] = detrend_fun(np.arange(wsize),A[ix[i] : ix[i] + wsize,j]) if detrend else A[ix[i] : ix[i] + wsize,j]
        
        #Set up tapering window
        #######################
        beta=np.pi*alpha
        hamm = np.hamming(wsize)
        hann = np.hanning(wsize)
        kbess = np.kaiser(wsize,beta)
        blackman = np.blackman(wsize)
        notaper = np.ones(wsize) #overpass tapering option
        gain=1.0
        
        if isinstance(tapering,bool) : which='hamm'
        elif isinstance(tapering,str) :
            if tapering.upper() == 'HAMMING' :
                which='hamm'
                gain=np.sum(hamm)/wsize #0.530416666667
            elif tapering.upper() == 'HANNING' :
                which='hann'
                gain=np.sum(hann)/wsize #0.489583333333
            elif tapering.upper() == 'KAISER' :
                which='kbess'
                gain=np.sum(kbess)/wsize #0.394170357504
            elif tapering.upper() == 'NONE' :
                which='notaper'
                gain=1.0
            elif tapering.upper() == 'BLACKMAN' :
                which='blackman'
                gain=np.sum(blackman)/wsize
            else : raise Exception('Unknown taper {0}'.format(tapering))
        elif isinstance(tapering,np.ndarray) : pass
        else :
            raise Exception('Bad value for tapering keyword')
        if not isinstance(tapering,np.ndarray) : exec('window='+which)
        else : window=tapering
        window = np.repeat(window,nn*nr).reshape((wsize,nn,nr))
    
        #Apply tapering on segmented data
        A=dum.copy()*window
        A=A.reshape(wsize,nr*nn) #Reshapa matrix
        nr=nn*nr
    else :
        if detrend :
            for i in np.arange(nr): A[:,i] = detrend_fun(np.arange(N),A[:,i]) if detrend else A[:,i]
        gain=1.0
    
    #Run transform
    ###############
    for i in np.arange(nr):
        spec=get_spec(dx, A[:,i],integration=integration,gain=gain,res_factor=res_factor)
        if i == 0:
            esd = spec['esd']
            psd = spec['psd']
            fq = spec['fq']
        else : 
            esd = np.append(esd,spec['esd'])
            psd = np.append(psd,spec['psd'])
    
#    factor=((A[:,0]-A[:,0].mean())**2).sum()/spec['esd'].sum()
    
    #Average spectrum
    #################
    nf=len(fq)
    p=1./fq
    esd=esd.reshape(nr,nf)
    psd=psd.reshape(nr,nf)
    esd=(np.sum(esd,axis=0)/nr)#/gain
    psd=(np.sum(psd,axis=0)/nr)#/gain
    
    psd = psd * (gain**0.5)
#    print gain, np.sqrt(gain), gain **2, gain*0.5, gain/2.
#    esd=(np.sum(esd,axis=0))#/gain
#    psd=(np.sum(psd,axis=0))#/gain


    #Normalise by energy content    
    Scaling_Factor=len(fq)/esd.sum()
    if normalise :
        esd*=Scaling_Factor
        psd*=Scaling_Factor
    
    if tapering is not None : return {'params':{'tapering':tapering is not None,'which':which,'wsize':int(wsize),'nwind':int(nn),'overlap':int(100.*overlap),'gain':gain},'psd':psd,'esd':esd,'fq':fq,'p':p}
    else : return {'params':{'tapering':tapering is not None},'psd':psd,'esd':esd,'fq':fq,'p':p}
    
def periodogram_analysis(dx,Ain,res_factor=10.0,detrend=False,normalise=False,average=True):
    """
    Periogram_Analysis
     @summary: This function computes a time averaged wavelet periodogram of SLA profiles.
     
     @param dx {type:numeric} : sampling distance
     @param Ain {type:numeric} : 2D table of sla data with time along 2nd axis (NXxNT with NX the spatial length and NT the time length)
     @keyword normalise {type:bool,default:False} : If True, normalise the spectrum by its overall energy content.
     @keyword detrend {type:bool,default:False} : If True, removes a linear trend to the segmented signal (if tapered) or to the whole signal (if not tapered).
     @return: a 2D spectrum structure ((NX)x(Number of scales))with Energy Spectral Density ('esd'), Power Spectral Density ('PSD'), frequency ('fq'), wavelength ('p') and distance ('D').
    
     @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
     @change: Created by RD, December 2012
    """
    
    A=Ain.copy()

    #Check dimensions
    sh = A.shape
    ndims = len(sh)
    N = sh[0] #Time series are found along the last dimension
    
    D = np.arange(N) * dx
    
    #If vector, add one dimension
    if ndims == 1 :
        A = A.reshape((N,1))
        sh = A.shape
        ndims = len(sh)
    
    nr = sh[1] #Number of repeats  
    nt = nr
    
    
    if detrend :
        for i in np.arange(nr): A[:,i] = detrend_fun(np.arange(N),A[:,i]) if detrend else A[:,i]
    gain=1.0
    
    #Run transform
    ###############
    for i in np.arange(nr):
        spec=get_spec(dx, A[:,i],gain=gain,periodogram=True,res_factor=res_factor)
        if i == 0:
            esd = np.ma.expand_dims(spec['esd'],0)
            psd = np.ma.expand_dims(spec['psd'],0)
            fq = spec['fq']
        else : 
            esd = np.ma.concatenate((esd,np.expand_dims(spec['esd'],0)),axis=0)
            psd = np.ma.concatenate((psd,np.expand_dims(spec['psd'],0)),axis=0)
    
#    factor=((A[:,0]-A[:,0].mean())**2).sum()/spec['esd'].sum()
    
    #Average spectrum
    #################
    nf=len(fq)
    p=1./fq
    
    if average :
        mask=psd.mask[0,:,:]
        esd=(np.sum(esd,axis=0)/nr)#/gain
        psd=(np.sum(psd,axis=0)/nr)#/gain
        #Update mask
        esd.mask=mask
        psd.mask=mask
    
        
    psd = psd * (gain**0.5)
#    print gain, np.sqrt(gain), gain **2, gain*0.5, gain/2.
#    esd=(np.sum(esd,axis=0))#/gain
#    psd=(np.sum(psd,axis=0))#/gain


    #Normalise by energy content    
    Scaling_Factor=len(fq)/esd.sum()
    if normalise :
        esd*=Scaling_Factor
        psd*=Scaling_Factor
    
    return {'psd':psd,'esd':esd,'fq':fq,'p':p,'D':D}