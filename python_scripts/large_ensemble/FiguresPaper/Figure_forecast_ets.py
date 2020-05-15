# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pickle as pkl

import common_plot_functions as cpf

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/'

figname='./Figure_forecast_ets'

filename='fcst_scores.pkl'

exps=['LE_D1_1km_5min_OFP_V2','LE_D1_1km_2min_OFP_V2','LE_D1_1km_1min','LE_D1_1km_30sec_OFP_V2','LE_D1_1km_5min_4D_OFP_V2','LE_D1_1km_1min_4D']

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

parameter = dict()   

ensemble_mean = dict()

for iexp , my_exp in enumerate(exps)   :

  #=========================================================
  #  READ THE DATA
  #=========================================================

  my_file=basedir + '/' + my_exp + '/time_mean/' + filename

  filehandler=open( my_file , 'rb' )
  parameter[my_exp] = pkl.load( filehandler ) 

print(' Finish the loop over experiments ' )    

print( np.shape( parameter['LE_D1_1km_5min_OFP_V2']['ref_ets'] ) )

plt.figure()
plt.plot(  parameter['LE_D1_1km_5min_OFP_V2']['mref_csi'][1,1,:]  ,'r')
plt.plot(  parameter['LE_D1_1km_1min']['mref_csi'][1,1,:] ,'b')


#ncols=3
#nrows=2
#fig, axs = plt.subplots( nrows,ncols , subplot_kw=dict(projection=ccrs.PlateCarree() ), figsize=[10,6.5] , sharex = 'col' , sharey = 'row')
#fig.subplots_adjust(wspace=0.075,hspace=0.28,bottom=0.03,left=0.05,right=0.97,top=0.97)
#icol = 0
#irow = 0

#xtick=[134.5,135,135.5,136,136.5,137]
#ytick=[34,34.5,35,35.5] 
#titles = ['(a) - 5MIN ','(b) - 2MIN ','(c) - 1MIN ','(d) - 30SEC ','(e) - 5MIN-4D ','(f) - 1MIN-4D ']

#ax=axs[irow,icol]
plt.show()
plt.savefig( figname + '.png' , format='png' , dpi=300  )
plt.close()


