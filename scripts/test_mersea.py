# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mpl_dates
#import enthought.chaco.shell as plt

import alti_tools as AT
import kernel as ke
        
if __name__ == "__main__" :
    
    limit=[40.0,4,44.5,11.0]
    trange_str = ['01/01/2008','08/05/2012']
    verbose = 0
    
    #Data filtering
#    p=14.
    p=27.
#    p=54.
#    p =81.
    q=p
    
    loess = True
    pl04 = True
    
    sat='j2'
 
    trange=[AT.cnes_convert(trange_str[0])[0],AT.cnes_convert(trange_str[1])[0]]

    #Setup directories
    ##################
#    alti_pattern = "C:\\VMShared/data/alti/regional/mersea/{0}_cf/nrt_mersea*.nc".format(sat)
    alti_pattern = "C:\\VMShared/data/alti/regional/mersea-dt/{0}_cf/dt_mersea*.nc".format(sat)

    
    #Load data
    ##########

    alti=AT.alti_data(alti_pattern,limit=limit,verbose=verbose,time_range=trange) #Load data 
    alti.reorder() #2D reordering of the data
    
    track_list = [9]
#    track_list = alti.track_list()
    
    #Loop over tracks
    #################
    for t in track_list:
        
        #Load variables and split by track number
        fg=alti.slice('track',t)
        lon=alti.lon[fg]
        lat=alti.lat[fg]
        dst=AT.calcul_distance(lat,lon)
        N=len(lon)
        dx=np.median(AT.deriv(dst))
        time = AT.grid_time(alti.date[:,fg])
        sla = alti.sla[:,fg]
#        time = np.mean(time,axis=1)
        nt=len(time)
        dt=np.median(AT.deriv(time))

        #Loess filtering of SLA
#        for i in np.arange(nt):
#            sla[i,:]=AT.loess(sla[i,:], dst, 40.)
        
        #Run wavelet analysis
        #####################
        
        #WV analysis
        sa_spectrum, sa_lscales, wvsla, daughter = ke.runAnalysis(lon,lat,time,sla,len_range=[40,150],w0=0)
        
        #Detection of energy peaks on scale-averaged spectrum
        eind = ke._2Ddetection(sa_spectrum, sa_lscales, amplim=0.04, win_width=3)
        
        #Sort indexes against time
        isort=np.argsort(time[eind[0]])
        eind[0][:]=eind[0][isort]
        eind[1][:]=eind[1][isort]
        
        #Get eddy properties
        diameter, symmetric= ke.decorrelation_scale(sla, lat, lon, eind)
        wvdiameter, wvsymmetric = ke.decorrelation_scale(wvsla, lat, lon, eind)
        ugdiameter, relvort = ke.solid_body_scale(sla, lat, lon, eind)
        amplitude = ke.eddy_amplitude(np.sqrt(sa_spectrum), eind)*100.
        cyclone = ke.cyclone(wvsla,eind)
        
#        #Sort data against time
#        isort=np.lexsort((diameter,wvdiameter,ugdiameter,relvort,amplitude,time[eind[0]]))
#        eind[0][:]=eind[0][isort]
#        eind[1][:]=eind[1][isort]
#        diameter[:]=diameter[isort]
#        wvdiameter[:]=wvdiameter[isort]
#        ugdiameter[:]=ugdiameter[isort]
#        relvort[:]=relvort[isort]
#        amplitude[:]=amplitude[isort]
        
        #Rebin results in space
        hist, ind, blon, blat = ke.grid(lon,lat,eind[0])
        ampmn, amprms= ke.grid_var(amplitude,hist,ind)
        lenmn, lenrms= ke.grid_var(diameter,hist,ind)
        wvlenmn, wvlenrms= ke.grid_var(wvdiameter,hist,ind)
        uglmn, uglrms= ke.grid_var(ugdiameter,hist,ind)
        rvmn, rvrms= ke.grid_var(relvort,hist,ind)
        
        #Rebin results in time
        date,datetime=AT.cnes_convert(time)
        thist,tind,btime=ke.grid_time(time,eind[1])
        trvmn,trvrms=ke.grid_var(relvort,thist,tind)
        tampmn,tamprms=ke.grid_var(amplitude,thist,tind)
        tlenmn,tlenrms=ke.grid_var(diameter,thist,tind)
        twvlmn,twvlrms=ke.grid_var(wvdiameter,thist,tind)
        tuglmn,tuglrms=ke.grid_var(ugdiameter,thist,tind)
        
        #Plot results
        #############
        
        #Show SLA and SA spectrum hovmollers with eddy positions overlaid
#        plt.subplot(1,2,1); plt.pcolor(lat,time/365.25+1950,np.sqrt(sa_spectrum)); plt.title('Wavelet scale-averaged spectrum');plt.xlabel('Latitude (�N)');plt.ylabel('Date'); plt.colorbar(); plt.plot(lat[eind[0]]+np.median(AT.deriv(lat))/2,(time/365.25+1950)[eind[1]]+(dt/365.25)/2,'ok',markersize=5); plt.subplot(1,2,2); plt.pcolor(lat,time/365.25+1950,wvsla - np.repeat(wvsla.mean(axis=1),N).reshape((nt,N))); plt.colorbar(); plt.title('wavelet-filtered SLA (m)');plt.xlabel('Latitude (�N)');plt.ylabel('Date');plt.plot(lat[eind[0]]+np.median(AT.deriv(lat))/2,(time/365.25+1950)[eind[1]]+(dt/365.25)/2,'ok',markersize=5); plt.show()
        
        #Show maps of binned eddy properties
        pmap=AT.plot_map(0,0,0,limit=alti.limit,resolution='i')
        pmap.title('Spatial evolution of eddy amplitude (cm)') ;pmap.scatter(blon,blat,ampmn,vmin=5,vmax=10,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Spatial evolution of eddy lengthscale (km)') ;pmap.scatter(blon,blat,lenmn,vmin=50,vmax=80,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Spatial evolution of eddy core scale (km)') ;pmap.scatter(blon,blat,uglmn,vmin=0,vmax=40,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Spatial evolution of relative vorticity (% of f)') ;pmap.scatter(blon,blat,rvmn/AT.coriolis(42.),vmin=0.0,vmax=0.8,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS of eddy amplitude (cm)') ;pmap.scatter(blon,blat,amprms,vmin=0,vmax=5,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS of eddy lengthscale (km)') ;pmap.scatter(blon,blat,lenrms,vmin=0,vmax=50,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS  of eddy core scale (km)') ;pmap.scatter(blon,blat,uglrms,vmin=0,vmax=40,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS of relative vorticity (% of f)') ;pmap.scatter(blon,blat,rvrms/AT.coriolis(42.),vmin=0,vmax=0.5,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Eddy observation frequency (%)') ;pmap.scatter(blon,blat,(100.0*hist)/nt,vmin=0,vmax=50,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        
        #Plot the time evolution of amplitude, scale and relative vorticity
        plt.subplot(3,1,1);plt.plot(datetime,tlenmn);plt.plot(datetime,AT.loess(tlenmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.title('Time evolution of eddy lengthscale');plt.ylim((0,150));
        plt.subplot(3,1,2);plt.plot(datetime,tampmn);plt.plot(datetime,AT.loess(tampmn, time, 90.));plt.ylabel('Amplitude (cm)');plt.title('Time evolution of eddy amplitude');plt.ylim((0,12));
        plt.subplot(3,1,3);plt.plot(datetime,trvmn/AT.coriolis(42.));plt.plot(datetime,AT.loess(trvmn, time, 90.)/AT.coriolis(42.));plt.xlabel('Date');plt.ylabel('Relative vorticity (% of f)');plt.title('Time evolution of the relative vorticity');plt.ylim((0.,1.0));plt.show()
                
        plt.subplot(3,1,1);plt.plot(datetime,tuglmn);plt.plot(datetime,AT.loess(tuglmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.title('Time evolution of eddy lengthscale');plt.ylim((0,50));
        plt.subplot(3,1,2);plt.plot(datetime,tampmn);plt.plot(datetime,AT.loess(tampmn, time, 90.));plt.ylabel('Amplitude (cm)');plt.title('Time evolution of eddy amplitude');plt.ylim((0,12));
        plt.subplot(3,1,3);plt.plot(datetime,trvmn/AT.coriolis(42.));plt.plot(datetime,AT.loess(trvmn, time, 90.)/AT.coriolis(42.));plt.xlabel('Date');plt.ylabel('Relative vorticity (% of f)');plt.title('Time evolution of the relative vorticity');plt.ylim((0.,1.0));plt.show()

        #Compare decorrelations scales from signal and wavelet plus core scale
        plt.subplot(3,1,1);plt.plot(datetime,tlenmn);plt.plot(datetime,AT.loess(tlenmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.title('Time evolution of eddy lengthscale');plt.ylim((0,150));
        plt.subplot(3,1,2);plt.plot(datetime,twvlmn);plt.plot(datetime,AT.loess(twvlmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.ylim((0,150));
        plt.subplot(3,1,3);plt.plot(datetime,tuglmn);plt.plot(datetime,AT.loess(tuglmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.ylim((0,50));plt.show()
        
    print 'done'