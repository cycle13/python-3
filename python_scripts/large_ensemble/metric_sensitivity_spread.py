# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import matplotlib.pyplot as plt
import os

from common_functions import common_functions as comm

nbv=1000
#Estimate number of bins in histogram. 
nbins=int(np.round(np.sqrt(nbv))) #Number of bins to be used in histogram computation.

my_ensemble=np.zeros( (1,1,1,nbv) )
my_undefmask=np.zeros( (1,1,1) , dtype=bool )

npvalues=10
pvalues=np.arange(1/npvalues,3*(1+1/npvalues),1/npvalues)

nsamples = 100

metrics=np.zeros( (nsamples,7,np.size(pvalues)) ) 

sample_mean=np.zeros( (nsamples ,np.size(pvalues)))
sample_std =np.zeros( (nsamples ,np.size(pvalues)))

hist=np.zeros( ( nsamples , nbins , np.size(pvalues) ) )
hist_gaussfit=np.zeros( ( nsamples , nbins , np.size(pvalues) ) )

for ip in range(0,np.size(pvalues))   :
   print('computing for p=' + str(pvalues[ip]) )
   for isample in range( 0 , nsamples )   :
      sample   =np.random.randn( nbv ) * pvalues[ip]
      sample   = ( sample - np.mean(sample) ) / np.std(sample) 
      sample_sq=np.power( np.random.randn( nbv ) , 2)
      sample_sq= ( sample_sq - np.mean( sample_sq ) ) / np.std( sample_sq )
      #Linearly blend the Gaussian and the chi-square distributions 
      my_ensemble[0,0,0,:]= (1-0.6)*sample + 0.6*sample_sq

      sample_mean[isample,ip]=np.mean( np.squeeze( my_ensemble[0,0,0,:] ) )
      sample_std[isample,ip] =np.std ( np.squeeze( my_ensemble[0,0,0,:] ) )

      [ hist[isample,:,ip] , x ]= np.histogram( np.squeeze( my_ensemble )  , bins=nbins , range=(-8.0,8.0) )
      hist[isample,:,ip] = hist[isample,:,ip]/np.sum( hist[isample,:,ip] )

      x_center=0.5*(x[0:-1]+x[1:])
      hist_delta=x[1]-x[0]

      hist_gaussfit[isample,:,ip]=hist_delta*(1/(np.sqrt(2*np.pi*np.power(sample_std[isample,ip],2))))*np.exp( -0.5*np.power( (x_center-sample_mean[isample,ip])/(sample_std[isample,ip]) ,2) )

      metrics[isample,:,ip]=comm.compute_pdf_metric(my_ensemble=my_ensemble,my_undefmask=my_undefmask
                            ,nx=1,ny=1,nz=1,nbv=nbv,undef=-999.0,nbins=nbins)

metrics_mean=np.squeeze( np.mean( metrics , axis=0 ) )
metrics_std =np.squeeze( np.std ( metrics , axis=0 ) )
metrics_max =np.squeeze( np.max ( metrics , axis=0 ) )
metrics_min =np.squeeze( np.min ( metrics , axis=0 ) )

sample_mean_mean=np.mean( sample_mean ) 
sample_std_mean =np.mean( sample_std  )

hist_mean=np.mean( hist , axis = 0 )
hist_gaussfit_mean=np.mean( hist_gaussfit , axis = 0 )

plt.figure()
m0=plt.errorbar( pvalues , metrics_mean[0,:] , yerr=np.squeeze(metrics_std[0,:]), fmt='r-' , label='KLD' )
m1=plt.errorbar( pvalues , metrics_mean[1,:] , yerr=np.squeeze(metrics_std[1,:]), fmt='b-' , label='Hellinger')
m2=plt.errorbar( pvalues , metrics_mean[2,:] , yerr=np.squeeze(metrics_std[2,:]), fmt='m-' , label='Total Variation')
m3=plt.errorbar( pvalues , metrics_mean[3,:] , yerr=np.squeeze(metrics_std[3,:]), fmt='k-' , label='Renyi Divergence')
m4=plt.errorbar( pvalues , metrics_mean[4,:] , yerr=np.squeeze(metrics_std[4,:]), fmt='y-' , label='J-Sh-Divergence')
m5=plt.errorbar( pvalues , metrics_mean[5,:] , yerr=np.squeeze(metrics_std[5,:]), fmt='c-' , label='Bhattacharyy')
m6=plt.errorbar( pvalues , metrics_mean[6,:] , yerr=np.squeeze(metrics_std[6,:]), fmt='v-' , label='Energy Distance')

plt.legend(handles=[m0,m1,m2,m3,m4,m5,m6])

plt.ylim(0,1.0)
plt.xlim(0,1.0)

plt.savefig( './Figura_metricas_sensibilidad_spread.png' )


plt.figure()
m0,=plt.plot( pvalues , metrics_std[0,:]/metrics_mean[0,:] , 'r-' , label='KLD' )
m1,=plt.plot( pvalues , metrics_std[1,:]/metrics_mean[1,:] , 'b-' , label='Hellinger')
m2,=plt.plot( pvalues , metrics_std[2,:]/metrics_mean[2,:] , 'm-' , label='Total Variation')
m3,=plt.plot( pvalues , metrics_std[3,:]/metrics_mean[3,:] , 'k-' , label='Renyi Divergence')
m4,=plt.plot( pvalues , metrics_std[4,:]/metrics_mean[4,:] , 'y-' , label='J-Sh-Divergence')
m5,=plt.plot( pvalues , metrics_std[5,:]/metrics_mean[5,:] , 'c-' , label='Bhattacharyy')
m6,=plt.plot( pvalues , metrics_std[6,:]/metrics_mean[6,:] , 'v-' , label='Energy Distance')

plt.legend(handles=[m0,m1,m2,m3,m4,m5,m6])

plt.ylim(0,0.3)
plt.xlim(0,1.0)

plt.savefig( './Figura_signal2noise_sensibilidad_spread.png' )
plt.close()

for ip in range(0,np.size(pvalues))   :
  plt.figure() 
  plt.plot( x_center , hist_mean[:,ip] , 'r-' )
  plt.plot( x_center , hist_gaussfit_mean[:,ip] , 'k--' )
  plt.savefig( './Figura_histogram_' + str(pvalues[ip]) + '_sensibilidad_spread.png' )  
  plt.close()

# !1-KL divergence
# !2-Hellinger distance
# !3-Total variation distance
# !4-Renyi Divergence
# !5-Jensen Shannon Divergence
# !6-Bhattacharyy distance
# !7-Energy distance.




