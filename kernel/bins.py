import numpy as np
from esutils_stat import histogram
import alti_tools as AT

def grid_time(time,tid):
    
    nt=len(time)
    hist,R = histogram(tid,binsize=1, min=0, max=nt-1, rev=True, use_weave=False, verbose=0)
    ind = AT.histogram_indices(hist, R)
    
    btime=time.copy()
    btime[:]=0.
    for i in np.arange(nt) : btime[i]=np.mean(time[tid[ind[i]]])
    
    return hist, ind, btime
    

def grid(lon,lat,eid,binsize=7):
    
    #Get lon/lat grid
    nx=len(lon)
#    id=np.arange(nx)
    hist,R = histogram(eid,binsize=binsize, rev=True, use_weave=False, verbose=0)
    ind = AT.histogram_indices(hist, R)
    
    blon = np.repeat(np.NaN,len(hist))
    blat = np.repeat(np.NaN,len(hist))
    for i in np.arange(len(hist)) : blon[i]=np.mean(lon[eid[ind[i]]])
    for i in np.arange(len(hist)) : blat[i]=np.mean(lat[eid[ind[i]]])
    
    return hist, ind, blon, blat

def grid_var(var,hist,ind):
    
    nx=len(hist)
    mn = np.repeat(np.NaN,nx)
    rms = np.repeat(np.NaN,nx)
    for i in np.arange(nx) : mn[i]=np.mean(var[ind[i]])
    for i in np.arange(nx) : rms[i]=AT.rms(var[ind[i]])
    
    return mn, rms
    
    
#    
#    
#    hist,R=histogram(eid,binsize=binsize,LOCATIONS=loc,REVERSE_INDICES=R)
#    
#    #Get corrrspondance between histogram indices and detected eddies
#    nx = len(hist)
#    id = []
#    for i in np.arange(nx) : id.append(tuple(set(ind[i]).intersection(set(eid))))
#    id2=[]
#    for i in np.arange(len(id)) : id2.append([i for j in id[i]])
#    id3=[]
#    id4=[]
#    for i,ii in enumerate(id2) :
#        for k,kk in enumerate(ii) :
#            id3.append(kk)
#            id4.append(list(id[i])[k])
##         id2.append(tuple(set(eid[i]).intersection(set(id))))
#    
#    sid=np.zeros(len(eid),dtype=int)
#    for i,ii in enumerate(id4) : sid[i]=np.where(eid==ii)[0]
#    id5=sid.tolist()
#    for i,ii in enumerate(id2) :
#        dum=[]
#        for k,kk in enumerate(ii):
#            if len(id5) > 0 : dum.append(id5.pop(0))
#            else : break
#        id5.append(dum)
#        
#    mn = np.repeat(np.NaN,nx)
#    md = np.repeat(np.NaN,nx)
#    rms = np.repeat(np.NaN,nx)
#    for i in np.arange(nx) : mn[i]=np.mean(var[id5[i]])
#    for i in np.arange(nx) : md[i]=np.median(var[id5[i]])
#    for i in np.arange(nx) : rms[i]=AT.rms(var[id5[i]])
#    
#    #Get frequency
#    f=[len(i) for i in id5]
#    
#    return mn, md, rms, f