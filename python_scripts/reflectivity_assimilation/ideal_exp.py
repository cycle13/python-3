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
import numpy as np
import matplotlib.pyplot as plt

def calc_ref( q , assume_graupel = False )  :
    
    nor=8.0e6        #[m^-4]
    ror=1000.0e0     #[Kg/m3]
    pip=np.power( np.pi , 1.75 )   #factor
    cf =1.0e18 * 720 #factor 
    ro=1.0e0
    
    nog=4.0e4      #[m^-4]
    rog=913.0e0    #[Kg/m3] 

    if q < 1.0e-10   :
        q=1.0e-10
        
    z= cf * ( np.power( ro * q , 1.75 ) )
    z= z / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
  
    if assume_graupel  :
       z= ( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
       z= z * ( ( ro * q ) ** 1.6625 )

    z = 10.0*np.log10( z )
   
    
    return z


#Define qr ensemble
    
ens_size = 60    

qr_mean = 0.0003
qr_std  = 0.0003

dbz_tr = 10.0
obs=10.0
R_in=5.0 
dynamic_r = False

assume_graupel = False
    
qr_ens = qr_mean+qr_std*np.random.randn( ens_size )

qr_ens[ qr_ens < 0.0 ] = np.mean(qr_ens)  #Desviamos el ensamble hacia los valores positivos.
#qr_ens[0]=0.0025   #Agregamos 3 outliers.
#qr_ens[1]=0.0012
#qr_ens[2]=0.0012  


#Define dbz ensemble    
dbz_ens = np.zeros( np.shape( qr_ens ) )   
for ii in range( ens_size ):
    dbz_ens[ii] = calc_ref( qr_ens[ii] , assume_graupel = assume_graupel )
    
plt.figure()
#plt.plot( qr_ens , dbz_ens , 'ob' )




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

cov=np.cov( np.stack((qr_ens,dbz_ens), axis=0) )[0,1]
cov_trunc=np.cov( np.stack((qr_ens,dbz_ens_trunc), axis=0) )[0,1]

qr_var = np.var( qr_ens )
dbz_var = np.var( dbz_ens )
dbz_var_trunc = np.var( dbz_ens_trunc )
dbz_mean = np.mean( dbz_ens )
dbz_mean_trunc = np.mean( dbz_ens_trunc )


qr_update = ( cov / ( dbz_var + R ) ) * ( obs - dbz_mean )

qr_update_trunc = ( cov_trunc / ( dbz_var_trunc + R ) ) * ( obs - dbz_mean_trunc )

qr_a = qr_mean + qr_update

qr_a_trunc = qr_mean + qr_update_trunc

dbz_a = calc_ref( qr_a , assume_graupel = assume_graupel)
dbz_a_trunc = calc_ref( qr_a_trunc , assume_graupel = assume_graupel)

dbz_a_lin= dbz_mean + ( dbz_var / ( dbz_var + R ) ) * ( obs - dbz_mean )
dbz_a_trunc_lin= dbz_mean_trunc + ( dbz_var_trunc / ( dbz_var_trunc + R ) ) * ( obs - dbz_mean_trunc )

plt.plot(qr_a,dbz_a,'ob',markersize=20.0)

plt.plot(qr_a_trunc,dbz_a_trunc,'or',markersize=20.0)


plt.plot(qr_a,dbz_a_lin,'ob',markersize=15.0)
plt.plot(qr_a_trunc,dbz_a_trunc_lin,'or',markersize=15.0)

plt.plot(qr_mean,dbz_mean_trunc,'ok',markersize=10.0)    
plt.plot(qr_mean,dbz_mean,'ok',markersize=10.0)   
plt.plot(np.array([0,np.max(qr_ens)]),np.array([dbz_mean_trunc,dbz_mean_trunc]),'-k',markersize=10.0) 
plt.plot(np.array([0,np.max(qr_ens)]),np.array([dbz_a_trunc,dbz_a_trunc]),'--k',markersize=10.0) 


#Update the ensemble using a particle filter.

w= np.exp( -np.power(dbz_ens - obs , 2) / R_in )
w=w/np.sum(w) #Normalizo los pesos.

qr_a_pf = np.sum( qr_ens * w )
dbz_a_pf = calc_ref( qr_a_pf , assume_graupel = assume_graupel )
plt.plot(qr_a_pf,dbz_a_pf,'om',markersize=20.0)



print( 1e4*(qr_mean - qr_a) , qr_a , dbz_a )
#print( 1e4*(qr_mean - qr_a_trunc) )
print( 1e4*(qr_mean - qr_a_pf) , qr_a_pf , dbz_a_pf )
