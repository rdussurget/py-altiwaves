import numpy as np
import alti_tools as AT
from esutils_stat import histogram
from  scipy.ndimage.filters import maximum_filter1d

import matplotlib.pyplot as plt

'''
Created on 12 nov. 2012

@author: rdussurg
'''

def cyclone(sla,ind):
    return sla[ind[1],ind[0]] < 0

def eddy_amplitude(var,ind):
    return var[ind[1],ind[0]]
    
def solid_body_scale(var,lat,lon,ind):
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
        dst, dumlon, dumlat, dumsla, gaplen, ngaps, gapedges, interpolated = AT.grid_track(lat[fg],lon[fg],cursla)
        ugeo=AT.geost_1d(dumlon,dumlat,dumsla,pl04=True)
        
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
        dst, dumlon, dumlat, dumsla, gaplen, ngaps, gapedges, interpolated = AT.grid_track(lat[fg],lon[fg],cursla)
        
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
        zc_l = (np.where(AT.deriv(np.abs(acorr_l)) > acorr_l))[0] if nl >= 3 else []
        zc_r = (np.where(AT.deriv(np.abs(acorr_r)) > acorr_r))[0] if nr >= 3 else []
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

def getScales(sa_spectrum,sa_lscales, lon, lat, time, sla, wvsla, daughter, id, binsize=7) :
    
#PRO AW_getScales, xval, lon, lat, xtime, len, idmx, sla, daughtout, $
#  longi=longi, $
#  lati=lati, $
#  lenval=lenval, $
#  ampval=ampval, $
#  cycval=cycval, $
#  ancval=ancval, $
#  hist=hist, $
#  loc=loc, $
#  nx=nx, $
#  diameter=diameter, $
#  wvdiameter=wvdiameter, $
#  bin=bin, $
#  true_anom=true_anom, $
#  cyclone=cyclone

#  IF (~exist(true_anom)) THEN true_anom=0B

#;  id=AW_detectEddies(xval, len, idmx, idmxid=idmxid, rglen=rglen)
    shape=sa_spectrum.shape
    nt=shape[0]
    nx=shape[1]
    
    xid=id[1]
    yid=id[0]
    ne=np.size(xid)
    
    #Regrid results if neceassary
    #CHECK CONSISTENCY WITH LENGTH COMPUTATION
    dst=AT.calcul_distance(lat,lon)
    mndst=np.median(AT.deriv(dst))
    dst_grid=np.arange(np.ceil(dst.max()/mndst)+1.)*mndst
  
    #Check length
    #Wavelet based method
    #####################
#    wvdiameter=decorrelation_scale(daughter, dst_grid, id)
        
    #Check length using SLA signal directly
    #This algorithm has a problem as it does not include both faces of eddy
    #######################################################################
    diameter=decorrelation_scale(sla, lon, lat, dst_grid, id)
#  diameter=DBLARR(nt)
#  cyclone=DBLARR(nt)
#  slaval=DBLARR(nt)
    
    #Rebin results
    ##############
    dhist,R=histogram(yid, binsize=binsize, rev=True, use_weave=False, verbose=0)
    ind = AT.histogram_indices(dhist, R)
    dist_shist = np.repeat(np.NaN,len(dhist))
#    for i in np.arange(len(ind)) : dist_shist[i]=np.mean(np.abs(spd_mat[ind[i]]))
#    dist_centers = np.arange(0,max_dist+max_dist/20.,max_dist/20.) + (max_dist/(2*20.))
#  
#  hist=HISTOGRAM(idmx,BINSIZE=bin,LOCATIONS=loc,REVERSE_INDICES=R)
#  
#  lenval = DBLARR(N_ELEMENTS(loc),/NOZERO)
#  ampval = DBLARR(N_ELEMENTS(loc),/NOZERO)
#  cycval = DBLARR(N_ELEMENTS(loc),/NOZERO)
#  ancval = DBLARR(N_ELEMENTS(loc),/NOZERO)
#  notempty=WHERE(hist GT 0, cnt,COMPLEMENT=empty)
#  IF (cnt NE N_ELEMENTS(hist)) THEN lenval[empty]=!VALUES.D_NAN
#  IF (cnt NE N_ELEMENTS(hist)) THEN  ampval[empty]=!VALUES.D_NAN
#  
#  ;Get current position
#  lati=lat[loc]
#  longi=lon[loc]
#  
#  FOR j = 0, cnt - 1 DO BEGIN
#    k=notempty[j]
#    lenval[k]=MEDIAN(wvdiameter[R[R[k] : R[k+1]-1]],/DOUBLE)
#;    ampval[k]=(~true_anom)? $
#;      MEDIAN(SQRT(xval[R[R[k] : R[k+1]-1]]),/DOUBLE) : $
#;      MEDIAN(ABS(slaval[R[R[k] : R[k+1]-1]]),/DOUBLE)
#    ampval[k]=(~true_anom)? $
#      MEAN(SQRT(xval[R[R[k] : R[k+1]-1]]),/DOUBLE) : $
#      MEAN(ABS(slaval[R[R[k] : R[k+1]-1]]),/DOUBLE)
#    cycval[j]=TOTAL(cyclone[R[R[k] : R[k+1]-1]])
#    ancval[j]=TOTAL(~(cyclone[R[R[k] : R[k+1]-1]]))
#;    
#;    lenval[k]=MEAN(wvdiameter[R[R[k] : R[k+1]-1]],/DOUBLE)
#;    ampval[k]=MEAN(xval[R[R[k] : R[k+1]-1]],/DOUBLE)
#  ENDFOR
#  
#  
#  RETURN
#  
#;  SAVE, longi, lati, lenval, ampval, hist, loc, nx, bin, FILENAME='/home/ctoh/dussurge/IDLWorkspace/Biscay_xtrack/wavelets/analysis_data/MGdata9209/scales/scales_'+label+'_'+STRTRIM(tracks[i],2)+'.DOG'+STRTRIM(m0,2)+'.sav'  
# 
#  
#  
#
#END
