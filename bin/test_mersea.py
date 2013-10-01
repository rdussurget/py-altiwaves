# -*- coding: utf-8 -*-
'''
    TEST_MERSEA
    @summary: This is a testing script. It applies the wavelet transform to a DUACS residuals <br />
              data set (as found on AVISO's website), and compute some space/time-averages<br />
              for it.
    @author: Renaud DUSSURGET, LER/PAC IFREMER.
    @change: Create in November 2012 by RD.
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
import matplotlib.pyplot as plt

import altimetry.tools as AT
from altimetry.data import alti_data
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

    alti=alti_data(alti_pattern,limit=limit,verbose=verbose,time_range=trange) #Load data 
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
        sa_spectrum, sa_lscales, wvsla, daughter = ke.runAnalysis(lon,lat,time,sla,len_range=[40,150],m=0)
        
        #Detection of energy peaks on scale-averaged spectrum
        eind = ke.detection(sa_spectrum, sa_lscales, amplim=0.04)
        
        #Sort indexes against time
        isort=np.argsort(time[eind[0]])
        eind[0][:]=eind[0][isort]
        eind[1][:]=eind[1][isort]
        
        #Get eddy properties
        print '\n## Eddy characteristics\n##\t1) Whole dataset'
        amplitude, diameter, relvort, ugdiameter, ugamplitude, rk_relvort, rk_center, rk_diameter, self_advect = \
            ke.get_characteristics(eind,lon,lat,time,sla,wvsla,sa_spectrum,filter=filter,p=40.0)
        blon, blat, hist,  ampmn, lenmn, rvmn, amprms, lenrms, rvrms = \
            ke.bin_space(lon,lat,eind,amplitude,diameter,relvort,method='mean',verbose=verbose,binsize=7)
        datetime, btime, thist, tampmn, tlenmn, trvmn, tamprms, tlenrms, trvrms = \
            ke.bin_time(time,eind,amplitude,diameter,relvort,method='mean',verbose=verbose,binsize=3)

        
        #Plot results
        #############
        
        #Show SLA and SA spectrum hovmollers with eddy positions overlaid
        plt.subplot(1,2,1); plt.pcolor(lat,time/365.25+1950,np.sqrt(sa_spectrum)); plt.title('Wavelet scale-averaged spectrum');plt.xlabel('Latitude (�N)');plt.ylabel('Date'); plt.colorbar(); plt.plot(lat[eind[0]]+np.median(AT.deriv(lat))/2,(time/365.25+1950)[eind[1]]+(dt/365.25)/2,'ok',markersize=2); plt.subplot(1,2,2); plt.pcolor(lat,time/365.25+1950,wvsla - np.repeat(wvsla.mean(axis=1),N).reshape((nt,N))); plt.colorbar(); plt.title('wavelet-filtered SLA (m)');plt.xlabel('Latitude (�N)');plt.ylabel('Date');plt.plot(lat[eind[0]]+np.median(AT.deriv(lat))/2,(time/365.25+1950)[eind[1]]+(dt/365.25)/2,'ok',markersize=2); plt.show()
        
        #Show maps of binned eddy properties
        pmap=AT.plot_map(0,0,0,limit=alti.limit,resolution='i')
        pmap.title('Spatial evolution of eddy amplitude (cm)') ;pmap.scatter(blon,blat,ampmn,vmin=5,vmax=10,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Spatial evolution of eddy lengthscale (km)') ;pmap.scatter(blon,blat,lenmn,vmin=50,vmax=80,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
#         pmap.title('Spatial evolution of eddy core scale (km)') ;pmap.scatter(blon,blat,uglmn,vmin=0,vmax=40,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Spatial evolution of relative vorticity (% of f)') ;pmap.scatter(blon,blat,rvmn/AT.coriolis(42.),vmin=0.0,vmax=0.8,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS of eddy amplitude (cm)') ;pmap.scatter(blon,blat,amprms,vmin=0,vmax=5,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS of eddy lengthscale (km)') ;pmap.scatter(blon,blat,lenrms,vmin=0,vmax=50,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
#         pmap.title('RMS  of eddy core scale (km)') ;pmap.scatter(blon,blat,uglrms,vmin=0,vmax=40,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('RMS of relative vorticity (% of f)') ;pmap.scatter(blon,blat,rvrms/AT.coriolis(42.),vmin=0,vmax=0.5,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        pmap.title('Eddy observation frequency (%)') ;pmap.scatter(blon,blat,(100.0*hist)/nt,vmin=0,vmax=50,s=50); plt.colorbar(); pmap.setup_map(); pmap.show()
        
        #Plot the time evolution of amplitude, scale and relative vorticity
        plt.subplot(3,1,1);plt.plot(datetime,tlenmn);plt.plot(datetime,AT.loess(tlenmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.title('Time evolution of eddy lengthscale');plt.ylim((0,150));
        plt.subplot(3,1,2);plt.plot(datetime,tampmn);plt.plot(datetime,AT.loess(tampmn, time, 90.));plt.ylabel('Amplitude (cm)');plt.title('Time evolution of eddy amplitude');plt.ylim((0,12));
        plt.subplot(3,1,3);plt.plot(datetime,trvmn/AT.coriolis(42.));plt.plot(datetime,AT.loess(trvmn, time, 90.)/AT.coriolis(42.));plt.xlabel('Date');plt.ylabel('Relative vorticity (% of f)');plt.title('Time evolution of the relative vorticity');plt.ylim((0.,1.0));plt.show()
                
#         plt.subplot(3,1,1);plt.plot(datetime,tuglmn);plt.plot(datetime,AT.loess(tuglmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.title('Time evolution of eddy lengthscale');plt.ylim((0,50));
        plt.subplot(3,1,2);plt.plot(datetime,tampmn);plt.plot(datetime,AT.loess(tampmn, time, 90.));plt.ylabel('Amplitude (cm)');plt.title('Time evolution of eddy amplitude');plt.ylim((0,12));
        plt.subplot(3,1,3);plt.plot(datetime,trvmn/AT.coriolis(42.));plt.plot(datetime,AT.loess(trvmn, time, 90.)/AT.coriolis(42.));plt.xlabel('Date');plt.ylabel('Relative vorticity (% of f)');plt.title('Time evolution of the relative vorticity');plt.ylim((0.,1.0));plt.show()

        #Compare decorrelations scales from signal and wavelet plus core scale
        plt.subplot(3,1,1);plt.plot(datetime,tlenmn);plt.plot(datetime,AT.loess(tlenmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.title('Time evolution of eddy lengthscale');plt.ylim((0,150));
#         plt.subplot(3,1,2);plt.plot(datetime,twvlmn);plt.plot(datetime,AT.loess(twvlmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.ylim((0,150));
#         plt.subplot(3,1,3);plt.plot(datetime,tuglmn);plt.plot(datetime,AT.loess(tuglmn, time, 90.));plt.ylabel('Lengthscale (km)');plt.ylim((0,50));plt.show()
        
    print 'done'