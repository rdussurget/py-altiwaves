# -*- coding: utf-8 -*-
import numpy as np
from  scipy.ndimage.filters import maximum_filter1d
from altimetry.tools import grid_track, geost_1d, deriv

'''
Created on 12 nov. 2012

@author: rdussurg
'''

def cyclone(sla,ind):
    '''
    cyclone
    @summary: Test the rotation sense of detected eddies.
    @param sla: Along-track sea level anomalies for detected events.
    @param ind: Indices of the detected events.
    @return {array}: True if detected events are cyclones.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since : November 2012.
    @change: Create in November 2012 by RD.
    '''
    return sla[ind[1],ind[0]] < 0

def eddy_amplitude(sla,ind):
    '''
    eddy_amplitude
    @summary: Get the eddy amplitude of detected eddies.
    @param sla: Along-track sea level anomalies for detected events.
    @param ind: Indices of the detected events.
    @return {array}: Amplitudes.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since : November 2012.
    @change: Create in November 2012 by RD.
    '''
    return np.abs(sla[ind[1],ind[0]])
    
def solid_body_scale(var,lat,lon,ind):
    '''
    solid_body_scale
    @summary: Compute the diameter of eddy core using maxima of geostrophic velocities<br />
              computed on both sides of the eddy, and computes the equivalent Relative<br />
              Vorticity for a solid body rotating eddy.
    @note: This technique was first applied in Le Hénaff et al., 2012. Cyclonic<br />
           activity in the eastern Gulf of Mexico: characterization. Submitted to <br />
           Progress in Oceanography.
    @note: A 2nd order polynomial over 3-4 points around the velocity maximum is<br />
           computed to better detect its position.
    @note: Geostrophic velocities are computed using the Powell and Leben (2004)<br />
           methodology - powell_leben_filter_km() function. Filtering parameters are<br/>
           p=q=12km on each sides of the point.
           Powell, B.S., et R.R. Leben. 2004. « An Optimal Filter for Geostrophic Mesoscale<br/>
           Currents from Along-Track Satellite Altimetry ». Journal of Atmospheric and<br/>
           Oceanic Technology 21 (10) (octobre 1): 1633‑1642.
    @param var: variable on which to apply the analysis : SLA, wavelet-filtered SLA,<br />
                daughter wavelets, etc...
    @param lon, lat: Longitude/latitude arrays.
    @ind: Indices of detected eddies.
    @return diameter, relvort : Diameter (km) and Relative Vorticity (s-1) of detected eddies.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since : November 2012.
    @change: Create in November 2012 by RD.
    '''
    xid=ind[1]
    yid=ind[0]
    ne=np.size(xid)
    diameter=np.zeros(ne,dtype=np.float64)
    relvort=np.zeros(ne,dtype=np.float64)
    
    for j in np.arange(ne) :
#        print j
        #Extract current SLA profile
        cursla=var[xid[j],:]
        fg=~cursla.mask
        dumy = np.where(np.arange(len(cursla))[fg] == yid[j])[0][0] #update yid with compressed SLA array
        cursla=cursla[fg]
        cursla -= np.median(cursla)
        
        #Fill gaps
        dst, dumlon, dumlat, dumsla, gaplen, ngaps, gapedges, interpolated = grid_track(lat[fg],lon[fg],cursla)
        ugeo=geost_1d(dumlon,dumlat,dumsla,pl04=True)
        
        #Update yid with resampled SLA array
        dumy = np.arange(len(dst))[~interpolated][dumy]
        
#        #Shrink data if there are any gaps greater than 3 points (ie. no valid interpolation)
#        if (gaplen > 3).any() :
#            gapedges=gapedges[:,gaplen > 3] #Remove gaps smaller than 3 points
#            left=gapedges[1,np.where(gapedges[1,:] >= dumy)[0][0]] if (gapedges[1,:] < dumy).any() else 0
#            right=gapedges[0,np.where(gapedges[0,:] >= dumy)[0][0]] if (gapedges[0,:] >= dumy).any() else len(dumsla)-1
#            dumsla=dumsla[left:right+1]
            
        #Get sla profile on both sides of the eddy 
        dumsla_r = dumsla[dumy:]      #np.roll(cursla,-dumy)
        ugeo_r = ugeo[dumy:]
        dumsla_l = dumsla[dumy:0:-1]  #np.roll(cursla[::-1],dumy+1)
        ugeo_l = ugeo[dumy:0:-1]
        nr = len(dumsla_r)
        nl = len(dumsla_l)
        
        #If not enought data on one side, take the opposite
        if nr < 3 :
            dumsla_r = dumsla_l
            ugeo_r = ugeo_l
            nr=nl
        if nl < 3 :
            dumsla_l = dumsla_r
            ugeo_l = ugeo_r
            nl=nr
        
        #Detect local velocity maxima on both sides of eddy
        mx_l = np.where(maximum_filter1d(np.abs(ugeo_l),3) == np.abs(ugeo_l))[0] #Do not take into account eddy center at position 0
        mx_r = np.where(maximum_filter1d(np.abs(ugeo_r),3) == np.abs(ugeo_r))[0]
        mx_l = mx_l[mx_l != 0] #Rq: This happens when peak is found at eddy center... This could possibly avoided?
        mx_r = mx_r[mx_r != 0]
        
        #Replace with data from the opposite side if no maxima are found
        if len(mx_l) == 0 :
            dumsla_l = dumsla_r.copy()
            ugeo_l = ugeo_r.copy()
            nl=nr
            mx_l = mx_r[0]
        else : mx_l = mx_l[0]
            
        if len(mx_r) == 0 :
            dumsla_r = dumsla_l.copy()
            ugeo_r = ugeo_l.copy()
            nr=nl
            mx_r = mx_l[0]
        else : mx_r = mx_r[0]
        
        
        #Fit a 2nd order polynomial on the 4 points surrounding the extrema
        if mx_l != nl - 1 :
            if np.abs(ugeo_l[mx_l-1]) > np.abs(ugeo_l[mx_l+1]) : fit= np.polyfit(dst[mx_l-1:mx_l+3 if mx_l+3 <= nl else nl],ugeo_l[mx_l-1:mx_l+3 if mx_l+3 <= nl else nl], 2)
            else : fit= np.polyfit(dst[mx_l-2 if mx_l >= 2 else 0:mx_l+2],ugeo_l[mx_l-2 if mx_l >= 2 else 0:mx_l+2], 2)
        else : fit=np.polyfit(dst[mx_l-2 if mx_l >= 2 else 0:mx_l+1],ugeo_l[mx_l-2 if mx_l >= 2 else 0:mx_l+1], 2)
        diameter[j] += (-fit[1]) / (2 * fit[0])
        
        if mx_r != nr - 1 :
            if np.abs(ugeo_r[mx_r-1]) > np.abs(ugeo_r[mx_r+1]) : fit= np.polyfit(dst[mx_r-1:mx_r+3 if mx_r+3 <= nr else nr],ugeo_r[mx_r-1:mx_r+3 if mx_r+3 <= nr else nr], 2)
            else : fit= np.polyfit(dst[mx_r-2 if mx_r >= 2 else 0:mx_r+2],ugeo_r[mx_r-2 if mx_r >= 2 else 0:mx_r+2], 2)
        else : fit=np.polyfit(dst[mx_r-2 if mx_r >= 2 else 0:mx_r+1],ugeo_r[mx_r-2 if mx_r >= 2 else 0:mx_r+1], 2)
        diameter[j] += np.abs((-fit[1]) / (2 * fit[0]))
        
        #Compute relative vorticity
        relvort[j] = np.median(np.append(np.abs(ugeo_r[1:mx_r+1])/(dst[1:mx_r+1]*1e3),np.abs(ugeo_l[1:mx_l+1])/(dst[1:mx_l+1]*1e3)))
#        if (dumsla[dumy] > dumsla[dumy-1]) | (dumsla[dumy] > dumsla[dumy+1]) : relvort[j] *= -1 #Inver sign if anticyclonic  
        
    return diameter, relvort



def decorrelation_scale(var,lat,lon,ind):
    '''
    solid_body_scale
    @summary: Compute the decorrelation length-scale of detected eddies.
    @note: This is the original technique applied in :
           Dussurget, R, F Birol, R.A. Morrow, et P. De Mey. 2011. « Fine Resolution<br />
           Altimetry Data for a Regional Application in the Bay of Biscay ». Marine<br />
           Geodesy 2 (34): 1‑30. doi:10.1080/01490419.2011.584835.
    @note: A linear regression is applied between the two points around the decorrelation<br />
           scale to better detect its position.
    @note: If no sufficient data is found on one of both sides, eddy is considered as<br />
           symmetric and scales are thus only computed from one side.
    @param var: variable on which to apply the analysis : SLA, wavelet-filtered SLA,<br />
                daughter wavelets, etc...
    @param lon, lat: Longitude/latitude arrays.
    @ind: Indices of detected eddies.
    @return diameter, symmetric : Diameter (km) of detected eddies, and symmetric flag to<br />
            check whether symmetry assumption was used or not.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @since : November 2012.
    @change: Create in November 2012 by RD.
    '''
    xid=ind[1]
    yid=ind[0]
    ne=np.size(xid)
    diameter=np.zeros(ne,dtype=np.float64)
    symmetric=np.zeros(ne,dtype=bool)
    
    for j in np.arange(ne) :
#        print j
        #Extract current SLA profile
        cursla=var[xid[j],:]
        fg=~cursla.mask
        dumy = np.where(np.arange(len(cursla))[fg] == yid[j])[0][0] #update yid with compressed SLA array
        cursla=cursla[fg]
        cursla -= np.median(cursla)
        
        #Fill gaps
        dst, dumlon, dumlat, dumsla, gaplen, ngaps, gapedges, interpolated = grid_track(lat[fg],lon[fg],cursla)
        
        #Update yid with resampled SLA array
        dumy = np.arange(len(dst))[~interpolated][dumy]
        
        #Shrink data if there are any gaps greater than 3 points (ie. no valid interpolation)
#        if (gaplen > 3).any() :
#            gapedges=gapedges[:,gaplen > 3] #Remove gaps smaller than 3 points
#            left=gapedges[1,np.where(gapedges[1,:] >= dumy)[0][0]] if (gapedges[1,:] < dumy).any() else 0
#            right=gapedges[0,np.where(gapedges[0,:] >= dumy)[0][0]] if (gapedges[0,:] >= dumy).any() else len(dumsla)-1
#            dumsla=dumsla[left:right+1]
            
        #Get sla profile on both sides of the eddy 
        dumsla_r = dumsla[dumy:]      #np.roll(cursla,-dumy)
        dumsla_l = dumsla[dumy:0:-1]  #np.roll(cursla[::-1],dumy+1)
        
        #Compute the auto-correlation on each sides.
        nr = len(dumsla_r)
        nl = len(dumsla_l)
        acorr_l = np.zeros(nl)
        acorr_r = np.zeros(nr)
        lag_r=np.arange(nr)
        lag_l=np.arange(nl)
        for i,l in enumerate(lag_l) : acorr_l[i] = np.corrcoef(dumsla_l, np.roll(dumsla_l,l))[0][1]
        for i,l in enumerate(lag_r) : acorr_r[i] = np.corrcoef(dumsla_r, np.roll(dumsla_r,l))[0][1]
    
        #detect first zero crossing of auto-corr function with the derivative of its absolute
        zc_l = (np.where(deriv(np.abs(acorr_l)) > acorr_l))[0] if nl >= 3 else []
        zc_r = (np.where(deriv(np.abs(acorr_r)) > acorr_r))[0] if nr >= 3 else []
        zer_cross = []
        if len(zc_l) != 0 : zer_cross = np.append(zer_cross,zc_l[0])
        if len(zc_r) != 0 : zer_cross = np.append(zer_cross,zc_r[0])
        
        #Linearly interpolate the auto-correlation function to get the zero-crossing distance
        fit = np.ma.array(np.zeros((2,2)),mask=np.ones((2,2),dtype=bool))
        if len(zc_l) != 0 : fit[0,:]=np.ma.array(np.polyfit(dst[zc_l[0]-1:zc_l[0]+1],acorr_l[zc_l[0]-1:zc_l[0]+1], 1),mask=np.zeros((2),dtype=bool))
        if len(zc_r) != 0 : fit[1,:]=np.ma.array(np.polyfit(dst[zc_r[0]-1:zc_r[0]+1],acorr_r[zc_r[0]-1:zc_r[0]+1], 1),mask=np.zeros((2),dtype=bool))
        
        #If no decorrelation is found on one side, use twice the decorrelation on valid side
        if fit.mask.sum() == 2 :
            fit[fit.mask] = fit[~fit.mask]
            symmetric[j] = True
        
        #Get diameter
        diameter[j]=-((fit[0][1]/fit[0][0]) +(fit[1][1]/fit[1][0]))
    
    return diameter, symmetric

#def getScales(sa_spectrum,sa_lscales, lon, lat, time, sla, wvsla, daughter, id, binsize=7) :
#
#    shape=sa_spectrum.shape
#    nt=shape[0]
#    nx=shape[1]
#    
#    xid=id[1]
#    yid=id[0]
#    ne=np.size(xid)
#    
#    #Regrid results if neceassary
#    #CHECK CONSISTENCY WITH LENGTH COMPUTATION
#    dst=AT.calcul_distance(lat,lon)
#    mndst=np.median(AT.deriv(dst))
#    dst_grid=np.arange(np.ceil(dst.max()/mndst)+1.)*mndst
#  
#    #Check length
#    #Wavelet based method
#    #####################
##    wvdiameter=decorrelation_scale(daughter, dst_grid, id)
#        
#    #Check length using SLA signal directly
#    #This algorithm has a problem as it does not include both faces of eddy
#    #######################################################################
#    diameter=decorrelation_scale(sla, lon, lat, dst_grid, id)
#    
#    #Rebin results
#    ##############
#    dhist,R=histogram(yid, binsize=binsize, rev=True, use_weave=False, verbose=0)
#    ind = AT.histogram_indices(dhist, R)
#    dist_shist = np.repeat(np.NaN,len(dhist))


def get_characteristics(eind,lon,lat,time,sla,wvsla,sa_spectrum):
    #Sort indexes against time
    isort=np.argsort(time[eind[1]])
    eind[0][:]=eind[0][isort]
    eind[1][:]=eind[1][isort]
    
    #Detect eddies and select cyclones
    diameter, symmetric= decorrelation_scale(sla, lat, lon, eind)
    wvdiameter, wvsymmetric = decorrelation_scale(wvsla, lat, lon, eind)
    ugdiameter, relvort = solid_body_scale(sla, lat, lon, eind)
    amplitude = eddy_amplitude(np.sqrt(sa_spectrum), eind)*100.
    return amplitude, diameter, relvort, ugdiameter, wvdiameter