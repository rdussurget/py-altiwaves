from __future__ import print_function
import sys
import numpy as np

import alti_tools as AT
import kernel

import matplotlib.pyplot as plt #Optional - only for debugging


'''
Created on 9 nov. 2012

@author: rdussurg
'''

curper=0

def process(i,N,step=10):
    '''
    process : counter bar
    '''
    global curper
    per = np.int16((100.0 * i)/N)
    if (per >= (curper + 10)) or i == 0 :
        print('o',end='')
        
        curper=per
    
def runAnalysis(lon, lat, time, sla, \
             
                 #Analysis parameters
                 w0=None, \
                 len_range=[60.,450.], \
                 detrend=True, \
                 demean=True, \
                 mother='dog', \
                 high=None, \
                 
                 #Verbose option
                 verbose=True):
    
    '''
    run_analysis :  
    
    @summary: Run wavelet analysis on along-track altimetry sea level data as in Dussurget et al., 2011
    @author : Renaud DUSSURGET, LER/PAC IFREMER (RD)
    @change: Created on 8 Nov. 2011 from IDL version of IDL version of AltiWAVEs/kernel/aw_run_analysis.pro.
    '''
    #define constants
    km2cm = 1e5
    cm2km = km2cm ** (-1.0)
    
    #Size parameters
    nx=len(lon)
    nt=len(time)
        
    #Get distance array
    dst = AT.calcul_distance(lat,lon) #Along-track distance in km
    dstcm = dst * km2cm #convert in centimeters
    
    #Setup Wavelet Transform parameters
    ###################################
    dt = np.median(AT.deriv(dstcm))
    w0 =  6 if mother == 'morlet' else (2 if w0 is None else w0) #Wavenumber

    len2scale = np.sqrt(w0+0.5) / (2*np.pi)
    scale2len = (2*np.pi) / np.sqrt(w0+0.5) #"fourrier factor"
    
    T = max(dstcm)
    s0 = 2.0 * dt if mother == 'morlet' else (2*dt) * len2scale #smallest wavescale (devided by fourier wavelength)
    dj = (dt/T) * len2scale #set scale interval (remember this is scaled by s0^(-1/2) ; no units!
    J = np.fix((np.log((T*len2scale)/s0) / np.log(2.)) / dj).astype(int) #number of scales from s0 to T
#    NJ = np.fix(J+1)

    avg1, avg2 = len_range # Range of periods to average
    slevel = 0.95 # Significance level    
    alpha = 0.0 # Lag-1 autocorrelation for white noise
    
    #Setup wavelet object
    exec('mother_obj=kernel.wavelet.{0}({1})'.format(mother,w0))
    

    #Print wavelet analysis parameters
    ##################################
    if verbose :
        str_header = '\t===Wavelet analysis parameters===\n\twavelet:{0}, degree:{1}, scale range {2} {3} km\n\toptions : '.format(mother, w0, len_range[0], len_range[1])
        str_opt = ','.join([x for x in np.array(['demean','detrend','filter'])[np.array([demean,detrend,filter])==True]])
        print(str_header + str_opt)
        print('\tstatus : ',end='')


    #Process data
    #############
    
    ####Not activated :
#    #Filter data if necessary
#    IF (filter) THEN BEGIN
#        TwoNplusOne = nt ;Number of points of the weighting function (2n + 1)
#        dum = TRANSPOSE(sla)
#        nans=WHERE(~FINITE(dum),nancnt)
#        dum[nans] = 0.
#        dum0_filt = lanczos_filter(dum[0,*],0,high, fit[1],wg_len=TwoNplusOne,filter=filter)
#        resp = TRANSPOSE(filter.response # TRANSPOSE(MAKE_ARRAY(nx,VALUE=1D))) ;Replicate filter
#        filt2 = TRANSPOSE(DOUBLE( FFT(FFT(TEMPORARY(dum),1,DIMENSION=2)*TEMPORARY(resp), -1, DIMENSION=2)))
#        filt2[nans] = !VALUES.D_NAN
#        sla = TEMPORARY(filt2)
    
    #Compute time coverage
    acov=np.nansum(np.isfinite(sla),axis=1)    
    per=100.0*acov.astype(float)/nx
    fg = (per.data >= 0.15)
    count = fg.sum()
    if count == nt :
        enough = np.arange(nt)
    else : enough = np.where(fg)[0]
    
    #Remove mean if necessary
    if (demean) : sla -= np.repeat(np.nansum(sla, axis=1)/nx,nx).reshape((nt,nx))
    if detrend :
        sla[enough,:]= AT.detrend(dst,sla[enough,:])
    
    
    #Setup output variables
    #######################
    sa_spectrum = np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool))
    wvsla =  np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool))
    sa_lscales =  np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool))
    daughtout =  np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool))
    avg_sig = np.arange(nt,dtype=np.float64)
    Cpsi = np.arange(nt,dtype=np.float64) #FFT energy matrix
    lenscale = np.ma.array(np.zeros(nt),mask=np.ones(nt,dtype=bool)) #return corresponding lengthscales
    
    outsla = np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool))
    
    #Setup intermediate variables
    WPower=np.ma.array(np.zeros((J+1,nx),dtype=np.float64),mask=np.ones((J+1,nx),dtype=bool))   # Normalized wavelet power spectrum
    W = WPower.copy()
    daughter = np.ma.array(np.ones((J+1,nx))*np.complex128(0),mask=np.ones((J+1,nx),dtype=bool))
    
    
    #Run transform
    ##############
    
    #Loop on valid cycles
    for i,valid in enumerate(enough):
#        print(i,nx)
        process(i,count)
        
        fg = np.isfinite(sla[valid,:]) if not isinstance(sla, np.ma.masked_array) else ~sla[valid,:].mask
        fgcnt = (~fg).sum()
        
        #Fill gaps when possible
        #Trucate dataset as well to remove edges
        if fgcnt > 0 :
            dum = sla[valid,:]
            dum, dumlon, dumlat, dumind, ngaps, gapedges, gaplen, interpolated = AT.fill_gaps(lat, lon, sla[valid,:], ~fg,remove_edges=True) #Truncate serie if not full
        else : dum = sla[valid,:]
        
        ndum=len(dum)
        
#        if ngaps > 0:
#            if gaplen.max() >= 3 :
#                print('long_gap')
        
        std = dum.std() # Standard deviation
        std2 = std ** 2 # Variance
        
#        dum = (dum - dum.mean())  / std
        
        #Reset intermediate arrays
        scale=0
        per=0
        
        #run transform        
        wave, scale, wavenb, coi_temp, daughter_temp, fft, fftfreqs = kernel.wavelet.cwt(dum, dt, dj, s0, J,
                                                      mother_obj)
        
        #Compute significance
        signif, fft_theor = kernel.wavelet.significance(1.0, dt, scale, 0, alpha,
                        significance_level=slevel, wavelet=mother_obj)
        
        Cd = mother_obj.cdelta
        length = 1. / wavenb
        
        #Compute power (WPower) and magnitude (W) of the transform
        WPower[:]=np.ma.array(0,mask=True)
        W[:]=np.ma.array(0,mask=True)
        daughter[:]=np.ma.array(np.complex128(0),mask=True)
        
        WPower[:,dumind] = (abs(wave)) ** 2     #Wavelet power spectrum
        W[:,dumind] = np.real(wave)             #Real part of the transform
        daughter[:,dumind] = daughter_temp
        
#        WPower = np.ma.array((abs(wave)) ** 2,mask=np.zeros((J+1,nx),dtype=bool))   # Normalized wavelet power spectrum
#        W = np.ma.array(np.real(wave),mask=np.zeros((J+1,nx),dtype=bool))
        sig95 = WPower / (signif * np.ones((nx, 1))).transpose() # Where ratio > 1, power is significant
        
        #Update Cone of Influence
        coi=np.repeat(s0*scale2len,nx)
        coi[dumind]=coi_temp
        coi[~fg[dumind]]=(s0*scale2len) #Put smallest possible value in masked areas (equivalent to 1/wavenb[0])
        
        #convert in km
        lengthkm = length * cm2km
        scalekm = scale * cm2km
        coikm = coi * cm2km
        
        #Mask data
        ##########
        # 1) Out of confidence interval
        # 2) Not significant at 95% (white noise)
        # 3) If confidence interval is too low wrt. smallest integration scale (avg1)
        coimask = np.repeat(coi,J+1).reshape((nx,J+1)).transpose() <= np.repeat(length,nx).reshape((J+1,nx))
        sig95mask = ~(sig95 > 1)
        int_fg = coikm < avg1 #RQ THIS MASK IS SIMILAR AND LESS EFFICIENT THAN PREVIOUS ONE
        int_fg_tab = np.repeat(int_fg,J+1).reshape((nx,J+1)).transpose() #Points where SA spectrum is valid (COI > avg1)
        data_mask = np.ones((J+1,nx),dtype=bool)
        data_mask[:,fg]=False
#        mask=(coimask | sig95mask |  int_fg_tab) #Unsure about sig95 mask
        mask = (coimask | int_fg_tab | data_mask)
#        mask = coimask
        WPower.mask[:]= mask
        W.mask[:]=mask
        
        #Compute the scaled-average spectrum (eddy band (avg1,avg2)- TC98, eq. 24)
        scale_fg = (lengthkm >= avg1) & (lengthkm <= avg2) #Integration scales
        
        
#        avg_sig[valid] = signif[scale_fg].mean()     
        
        #Scale-averaged spectrum
        dum_sa_spec = WPower / (scale * np.ones((nx, 1))).transpose()
        sa_spectrum[valid,:] = ((dj * dt) / Cd) * dum_sa_spec.sum(axis=0) #Integrate Power/scales over scales and normalize
        
        #Wavelet filtering (TC98 Eq. 29)
        dum_wvsla = W / (np.sqrt(scale) * np.ones((nx, 1))).transpose()
        wvsla[valid,:] = ((dj * np.sqrt(dt)) / (Cd * 1.) ) * dum_wvsla.sum(axis=0) #Integrate Power/scales over scales and normalize

#        #Compute Finite Energy constraint (Gu & Philander, 1995)
#        fq = np.fft.fftfreq(nx,dt)
#        psd = (abs(np.fft.fft(dum))**2.) /  fq   
        
#        # Scale average between avg1 and avg2 periods and significance level
#        sel = (lengthkm >= avg1) & (length < avg2)
#        Cdelta = mother.cdelta
#        scale_avg = (scale * np.ones((nx, 1))).transpose()
#        # As in Torrence and Compo (1998) equation 24
#        scale_avg = WPower / scale_avg
#        scale_avg = std2 * dj * dt / Cd * scale_avg[scale_fg, :].sum(axis=0)
        
        #
        scale_avg_signif, tmp = kernel.wavelet.significance(std2, dt, scale, 2, alpha,
                                    significance_level=slevel, dof=[scale[scale_fg][0],
                                    scale[scale_fg][-1]], wavelet=mother_obj)
        
        
        #Get points for which a valid spectrum value exists
        masked_pts=~(np.fix(~mask).sum(axis=0) > 0)
        
        scid = np.ma.array(WPower.argmax(axis=0),mask=masked_pts)
        mxid = sa_spectrum[valid,:].argmax()
        sa_lscales[valid,:]=np.ma.array(lengthkm[scid],mask=scid.mask)
        lenscale[valid]=sa_lscales[valid,:][mxid]
        #Check everyhing (original signal, most energetic lengthscale at each point & scale averaged spectrum)
#        plt.figure(0);plt.plot(np.ma.array(lengthkm[scid],mask=scid.mask)); plt.figure(1); plt.plot(dum);plt.figure(2);plt.plot(sa_spectrum[valid,:]);plt.show()
        
        #Get scaled daughter
        #This formulation normalise the amplitude of the daughter wavelet
        #It may be biased (not detrended)
        ddum = np.real(daughter.data[scid[mxid],:]) / np.real(daughter.data[scid[mxid],:]).max() * np.abs(np.sqrt(sa_spectrum[valid,mxid])) * (np.abs(sla[valid,mxid])/sla[valid,mxid])
        #Get scaled daughter
        daughtout[valid,:]=ddum
    
    
    #Get 1D output results (1 eddy/pass)
    #####
    #TODO : Get this function out of analysis code
    mx = np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool),dtype=np.float)
    truemx = np.ma.array(np.zeros((nt,nx)),mask=np.ones((nt,nx),dtype=bool),dtype=np.float)
    idmx = np.ma.array(np.zeros(nt),mask=np.ones(nt,dtype=bool),dtype=np.float)
    for i in np.arange(nt) :
        ind = sa_spectrum[i,:].argmax()
        mx[i,ind] = sa_spectrum[i,ind]
        truemx[i,ind] = sla[i,ind]
        idmx[i]=ind
    
    
    #find maximums
    mxval = mx.max(axis=1)
    tmxval = truemx.max(axis=1)

    if (verbose) : print('o done\n\t=================================')
    
    return sa_spectrum, sa_lscales, wvsla, daughtout 
#           mxval, idmx, truemx, tmxval, lenscale
#                 outsla=outsla
    
        
        
    