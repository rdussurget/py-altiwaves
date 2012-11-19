import numpy as np
import kernel as ke
import alti_tools as AT
import matplotlib.pyplot as plt


#Simulate sea level data from red noise
#######################################

#Position data set in space and time
lat = np.arange(43.0,44.0,0.01)
lon = np.arange(6,6.5,0.005)
dst=AT.calcul_distance(lat,lon)
dx=np.median(AT.deriv(dst))
N=len(lon)
nt=25
dt=9.9
time=22705.0 + np.arange(0,nt)*dt

#Red noise generation (lengthscale > 10 km)
#sla=np.cumsum(np.random.randn(N*nt)).reshape((nt,N))
sla=np.cumsum(np.cumsum(np.random.randn(nt,N),axis=1),axis=0)
for i in np.arange(nt):
    sla[i,:]=AT.loess(sla[i,:], dst, 10.)

#Run wavelet analysis
#####################

#WV analysis
sa_spectrum, sa_lscales, wvsla, daughtout = ke.runAnalysis(lon,lat,time,sla,len_range=[10,150],w0=0)

#Detection of energy peaks on scale-averaged spectrum
res = ke._2Ddetection(sa_spectrum, sa_lscales, amplim=1.0, win_width=5)

#Plot results
plt.subplot(2,1,1); plt.pcolormesh(dst,time,sa_spectrum); plt.colorbar(); plt.plot(dst[res[0]]+dx/2,time[res[1]]+dt/2,'ok'); plt.subplot(2,1,2); plt.pcolormesh(dst,time,sla - np.repeat(sla.mean(axis=1),N).reshape((nt,N))); plt.colorbar(); plt.plot(dst[res[0]]+dx/2,time[res[1]]+dt/2,'ok'); plt.show()

print 'done'