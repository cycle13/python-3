#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Tue Feb 26 23:30:13 2019

@author: jruiz

#Posibles analisis que se pueden desprender de esto.

Sistematizacion de los resultados (barriendo rangos amplios de media del ensamble de 
qr y su dispersion. Realizando diferentes repeticiones de cada experimento para colectar
estadisticas respecto de los "updates catastroficos" y de las diferencias entre los diferentes metodos)

Implementar un ful EnKF (LETKF EnSRF) para ver no solo que pasa con la media sino tambien
con los diferentes miembros.

Realizar experimentos mas sofisticados en donde existan 2 puntos a los que se hace el update
uno uno con una observacion con un O-B postivio y grande y otra con un O-B pequenio y en donde 
el state tiene una covarianza importante entre las dos (para ver si una obs catastrofica
puede tambien afectar puntos vecinos a pesar de que en esos puntos la obs no conduzca es ese comp
ortamiento).

Realizar experimentos en donde aparte de q se agregue alguna otra variable como W
(se pueden tomar ensambles reales sacados de los experimentos para hacer esto).
Ver el impacto sobre W y sobre otras variables de los "updates catastroficos".


Para sistematizar este script barrer el spacio de media y dispersion y ahi crear
histogramas del valor de q final, dbz final, de la diferencia update pf vs update enkf
de la diferencia update truncado / update sin truncar , calcular la frecuencia de
updates "catastroficos". Para cada combinacion de valores hacer 1e6 repeticiones del 
experimento. 

"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')
sys.path.append('../../common_python/common_letkf/')

import numpy as np
import matplotlib.pyplot as plt
from cletkf           import common_letkf        as cletkf
from cletkf           import common_da           as cda

#Define qr ensemble
    
ens_size =60    

qr_mean = 0.0004
qr_std  = 0.0004

dbz_tr = 10.0
obs=15.0
R_in=5.0 

infl = 0.5e0

dynamic_r = False 
assume_graupel = False
    
qr_ens = qr_mean+qr_std*np.random.randn( ens_size )

qr_ens[ qr_ens < 0.0 ] = np.mean(qr_ens)  #Desviamos el ensamble hacia los valores positivos.
#qr_ens[0]=0.0025   #Agregamos 3 outliers.
#qr_ens[1]=0.0012
#qr_ens[2]=0.0012  
t=287.0
p=100000.0

#Define dbz ensemble    
dbz_ens = np.zeros( np.shape( qr_ens ) )   
for ii in range( ens_size ):
    if assume_graupel :
       dbz_ens[ii]=cda.calc_ref(0.0,0.0,qr_ens[ii],t,p)
    else              :
       dbz_ens[ii]=cda.calc_ref(qr_ens[ii],0.0,0.0,t,p)
    
plt.figure()

dbz_ens_trunc=np.copy( dbz_ens )
dbz_ens_trunc[ dbz_ens_trunc < dbz_tr ] = dbz_tr

plt.plot( qr_ens , dbz_ens_trunc , 'or' )

if dynamic_r :
   #R = np.max([ R_in , 0.5*abs( obs-np.mean(dbz_ens) ) ])
   if np.mean(dbz_ens_trunc) > obs   :
      R=0.35*(np.max(dbz_ens_trunc)-np.min(dbz_ens_trunc))
   else                              :
      R=R_in
      
   R= np.power( R , 2 )
   print( np.sqrt(R) )
else         :
   R=np.power( R_in , 2)   


#Fit linear regression
fit=np.polyfit(qr_ens,dbz_ens, 1)   
fit_trunc=np.polyfit(qr_ens,dbz_ens_trunc, 1)   

x=np.array([-1e-4,np.max(qr_ens)])
#plt.plot(x,x*fit[0]+fit[1],'--b')
plt.plot(x,x*fit_trunc[0]+fit_trunc[1],'--r')
plt.grid()


#Compute updates.

hdxb=np.zeros((1,ens_size))
hdxb[0,:]=dbz_ens_trunc - np.mean(dbz_ens_trunc)
dep =obs - np.mean(dbz_ens_trunc)

minfl=infl


[wa,wamean,pa]=cletkf.letkf_core(ne=ens_size,nobsl=1,hdxb=hdxb,rdiag=R,rloc=1.0,dep=dep,parm_infl=infl,minfl=minfl)

qr_a = np.mean( qr_ens ) * np.ones( ens_size )

for ii in range( ens_size ) :
    for ii2 in range( ens_size ) :

        qr_a[ii] = qr_a[ii] + (qr_ens[ii2]-np.mean(qr_ens))*( wa[ii2,ii] + wamean[ii2] )

        
qr_a[ qr_a < 0 ] = 0.0        

dbz_a = np.zeros( np.shape( qr_a ) )
for ii in range( ens_size ):
   if assume_graupel :
      dbz_a[ii]=cda.calc_ref(0.0,0.0,qr_a[ii],t,p)
   else              :
      dbz_a[ii]=cda.calc_ref(qr_a[ii],0.0,0.0,t,p)

if assume_graupel  :
    dbz_a_mean = cda.calc_ref(0.0,0.0,np.mean(qr_a),t,p)
else               :
    dbz_a_mean = cda.calc_ref(np.mean(qr_a),0.0,0.0,t,p)


plt.plot(np.mean(qr_a),np.mean(dbz_a),'or',markersize=20.0)

plt.plot(np.mean(qr_a),dbz_a_mean,'ob',markersize=20.0)

plt.plot(np.mean(qr_ens),np.mean(dbz_ens_trunc),'ok',markersize=20.0)

#Posterior ensemble.
plt.plot( qr_a , dbz_a , 'ob' )


#Update the ensemble using a particle filter.

w= np.exp( -np.power(dbz_ens - obs , 2) / R_in )
w=w/np.sum(w) #Normalizo los pesos.

qr_a_pf = np.sum( qr_ens * w )
if assume_graupel  :
    dbz_a_pf = cda.calc_ref(0.0,0.0,qr_a_pf,t,p)
else               :
    dbz_a_pf = cda.calc_ref(qr_a_pf,0.0,0.0,t,p)
plt.plot(qr_a_pf,dbz_a_pf,'om',markersize=20.0)


plt.figure()

plt.plot( dbz_ens )
plt.plot( dbz_a )


plt.figure()

plt.plot( qr_ens )
plt.plot( qr_a )

print( 1e4*(qr_mean - np.mean(qr_a) ) , np.mean(qr_a) , dbz_a_mean )
print( 1e4*(qr_mean - np.mean(qr_a_pf) ) , np.mean(qr_a_pf) , dbz_a_pf )
print( 1e4*np.std(qr_ens),1e4*np.std(qr_a))
print( np.std(dbz_ens),np.std(dbz_a))


