#-------------------------------------------------------------------------------
# Name:        Reading polar volume data
# Purpose:
#
# Author:      heistermann
#
# Created:     14.01.2013
# Copyright:   (c) heistermann 2013
# Licence:     MIT
#-------------------------------------------------------------------------------
#!/usr/bin/env python
 
import wradlib
import numpy as np
import pylab as pl
# just making sure that the plots immediately pop up
pl.interactive(True)
import datetime as dt
import os
import matplotlib.pyplot as plt
 
filename='/home/jruiz/EZE_PPIVol_20151014_232005.hdf'
 
raw = wradlib.io.read_OPERA_hdf5(filename)
 
sitecoords = (raw["where"]["lon"], raw["where"]["lat"],raw["where"]["height"])
 
where = raw["dataset%d/where"%1]
what  = raw["dataset%d/data2/what"%1]
ref = what["offset"] + what["gain"] * raw["dataset%d/data2/data"%1]
 
where = raw["dataset%d/where"%1]
what  = raw["dataset%d/data3/what"%1]
vr  = what["offset"] + what["gain"] * raw["dataset%d/data3/data"%1]
 
nrays = where["nrays"]
nbins = where["nbins"]
rscale = where["rscale"]
rstart = where["rstart"]
r  = np.arange(rstart, rstart+nbins*rscale,rscale)
 
r[0]=1
r=r/1000
 
#plt.show()
 
tmp=np.zeros((nrays,nbins))
tmp[:,:]= ref[:,:] 
tmp[ np.logical_and(vr > -6 , vr < 6) ]=-32
 
cref=np.zeros((nrays,nbins))
cref[:,:]=-32.0
 
offset=100
att=0.01
 
tmp[:,:offset]=-32
 
 
dbm=np.zeros((nbins))
w=np.zeros((nbins))
dbzrayo=np.zeros((nrays,nbins))
zrayo=np.zeros((nrays,nbins))
z    =np.zeros((nrays,nbins))
a=np.zeros(nrays)
b=np.zeros(nrays)
corrcoef=np.zeros(nrays)
raux=np.zeros(nbins)
#correction=ref
#correction[:,:]=0.0
 
refthreshold=7.0
min_sample_size=50
 
#Traducir a fortran? 
 
for i in range(nrays):
    dbm[:]=0.0
    raux[:]=0.0
    contador=0
    for j in range(nbins):
        if tmp[i,j] != -32.0:
           dbm[contador]= tmp[i,j]  - 20*np.log10( r[j] ) - 2*att*r[j]
           raux[contador]=r[j]
           #w[j]=1.0
           contador=contador+1
        #else:
           #w[j]=0.0
           #power[j]= np.power(10,tmp[i,j]/10.0)/np.square(r[j])
    if contador > min_sample_size:
       p=np.polyfit(r[:contador-1],dbm[:contador-1],1)
       tmpcorr=np.corrcoef(r[:contador-1],dbm[:contador-1])
       corrcoef[i]=tmpcorr[0,1]
       a[i]=p[0]
       b[i]=p[1]
    else:
       a[i]=0.0
       b[i]=0.0
            
    for j in range(nbins):
        #if ref[i,j] != -32.0:
        dbzrayo[i,j]=a[i]*r[j]+b[i] + 20*np.log10( r[j] ) + 2*att*r[j]
        if np.abs(dbzrayo[i,j] - ref[i,j]) < refthreshold:
           cref[i,j]=-32.0  
        else:
           zrayo = np.power(10,dbzrayo[i,j]/10.0)
           z     = np.power(10,ref[i,j]/10.0)
     
           if z - zrayo > 0.0:
              #We correct the reflectivity substracting WI-FI power
              cref[i,j]=10*np.log10(z-zrayo)
           else:
              cref[i,j]=-32.0
         
#Additional filter for the remaining echoes
#TODO: consider cyclic boundary conditions in azimuth.
npassfilter=3
for ifilter in range(npassfilter):
    for i in range(nrays):
        if np.logical_and( i> 2 , i< nrays-1):
            for j in range(nbins):
                if np.logical_and(cref[i-1,j]<0,cref[i+1,j]<0): #,cref[i,j]>0):
                    if cref[i,j] > 0:
                        cref[i,j]=-32.0
                if np.logical_and(cref[i-2,j]<0,cref[i+1,j]<0):
                    if cref[i,j] > 0:
                        cref[i,j]=-32.0
                if np.logical_and(cref[i-1,j]>0,cref[i+1,j]>0): #,cref[i,j]==-32.0):
                    if cref[i,j] == -32.0:
                        cref[i,j]=0.5*(cref[i-1,j]+cref[i+1,j])
                if np.logical_and(cref[i-2,j]>0,cref[i+1,j]>0):
                    if cref[i,j] == -32.0:
                        cref[i,j]=0.25*cref[i-2,j]+0.75*cref[i+1,j]
                if np.logical_and(cref[i-3,j]>0,cref[i+1,j]>0):
                    if cref[i,j] == -32.0:
                        cref[i,j]=0.10*cref[i-3,j]+0.90*cref[i+1,j]
               
 
cref[ cref < -32.0]=-32.0
                                 
    
