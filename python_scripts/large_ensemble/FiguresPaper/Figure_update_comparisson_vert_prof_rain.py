#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:54:08 2019

@author: jruiz
"""
import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os
import matplotlib.pyplot as plt
from common_functions import common_functions as cf
from scipy.ndimage.filters import gaussian_filter


import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

figname='./Figure_update_comparisson_vert_prof_rain.eps'

exps=['LE_D1_1km_5min','LE_D1_1km_30sec']


filetype='guesgp'   #analgp , analgz , guesgp , guesgz


variable_combination=[['dbz','dbz','dbz','dbz','v' ,'v' ,'v','v' ],
                      ['tk' ,'v'  ,'qv' ,'w'  ,'tk','v','qv','w' ]]

obs_increment = [5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]
obs_error     = [5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]

nbv=1000

sigma_smooth=2.0

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

profile_mean_rmsd=dict()
profile_mean_rmsu=dict()
profile_kld      =dict()


ctl_dict = ctlr.read_ctl( basedir + '/LE_D1_1km_5min/ctl/moment0001_for.ctl' )

for iexp , my_exp in enumerate(exps)   :

  profile_mean_rmsd[ my_exp ] = []
  profile_mean_rmsu[ my_exp ] = []
  profile_kld[ my_exp ] = []
  
  
  
  #=========================================================
  #  READ THE DATA
  #=========================================================

  for iv , my_obs_inc  in enumerate(obs_increment)  :
      
      var_obs = variable_combination[0][iv]
      var_upd = variable_combination[1][iv]
     
      my_file=basedir + '/' + my_exp + '/time_mean/guesgp/update_mean_rmsd_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      profile_mean_rmsd[my_exp].append( np.nanmean( np.load(my_file)['parameter'] , 1 ) )
      
      my_file=basedir + '/' + my_exp + '/time_mean/guesgp/update_mean_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      profile_mean_rmsu[my_exp].append( np.nanmean( np.load(my_file)['parameter'] , 1 ) )
      
      my_file=basedir + '/' + my_exp + '/time_mean/guesgp/update_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      profile_kld[my_exp].append( np.nanmean( np.load(my_file)['parameter'] , 1 ) )
      
 
print(' Finish the loop over experiments ' )    

#=========================================================================================
#Plot the mean KLD and its standard deviation.
#=========================================================================================

nexp = len( exps )

levels=ctl_dict['vlevels']
levels_str=list()

#To solve the postprocessing error.
levels=np.delete(levels,4,axis=0)
levels[3]=850.0


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


ncols=2
nrows=2

fig, axs = plt.subplots( nrows,ncols , figsize=[10,10] , )

fig.subplots_adjust(wspace=0.01,hspace=0.01,bottom=0.05,left=0.05,right=0.93,top=0.95)

titles = ['(a)','(b)','(c)','(d)']

#for key in kld_time_mean  :

ax = axs[0,0]

p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][0][1:] / profile_mean_rmsu['LE_D1_1km_5min'][0][1:] ) , -np.log( levels[1:] ) , 'b-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][0][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][0][1:] ) , -np.log( levels[1:] ) , 'b--' )
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][4][1:] / profile_mean_rmsu['LE_D1_1km_5min'][4][1:] ) , -np.log( levels[1:] ) , 'r-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][4][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][4][1:] ) , -np.log( levels[1:] ) , 'r--' )


ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
my_levs = levels[[0,3,4,6,8]]
ytick=-np.log(my_levs)
ax.set_yticks(ytick)
#Get the level string list.
levels_str=[]
for ilev in my_levs  :
   levels_str.append( str(int(ilev)) )
ax.set_yticklabels(levels_str,fontsize=10,color='k')
plt.setp(ax.get_xticklabels(), visible=False)
ax.set_xlim(left=0, right=50)
ax.text(1,-5.38,'(a) - $T$',fontsize=20,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})

ax = axs[0,1]

p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][1][1:] / profile_mean_rmsu['LE_D1_1km_5min'][1][1:] ) , -np.log( levels[1:] ) , 'b-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][1][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][1][1:] ) , -np.log( levels[1:] ) , 'b--' )
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][5][1:] / profile_mean_rmsu['LE_D1_1km_5min'][5][1:] ) , -np.log( levels[1:] ) , 'r-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][5][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][5][1:] ) , -np.log( levels[1:] ) , 'r--' )


ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
my_levs = levels[[0,3,4,6,8]]
ytick=-np.log(my_levs)
ax.set_yticks(ytick)
#Get the level string list.
levels_str=[]
for ilev in my_levs  :
   levels_str.append( str(int(ilev)) )
ax.set_yticklabels(levels_str,fontsize=10,color='k')
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_xlim(left=0, right=50)
ax.text(1,-5.38,'(b) - $V$',fontsize=20,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})


ax = axs[1,0]

p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][2][1:] / profile_mean_rmsu['LE_D1_1km_5min'][2][1:] ) , -np.log( levels[1:] ) , 'b-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][2][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][2][1:] ) , -np.log( levels[1:] ) , 'b--' )
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][6][1:] / profile_mean_rmsu['LE_D1_1km_5min'][6][1:] ) , -np.log( levels[1:] ) , 'r-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][6][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][6][1:] ) , -np.log( levels[1:] ) , 'r--' )


ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
my_levs = levels[[0,3,4,6,8]]
ytick=-np.log(my_levs)
ax.set_yticks(ytick)
#Get the level string list.
levels_str=[]
for ilev in my_levs  :
   levels_str.append( str(int(ilev)) )
ax.set_yticklabels(levels_str,fontsize=10,color='k')
ax.set_xlim(left=0, right=50)
ax.text(1,-5.38,'(c) - ${q}_{v}$',fontsize=20,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})


ax = axs[1,1]

p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][3][1:] / profile_mean_rmsu['LE_D1_1km_5min'][3][1:] ) , -np.log( levels[1:] ) , 'b-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][3][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][3][1:] ) , -np.log( levels[1:] ) , 'b--' )
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_5min'][7][1:] / profile_mean_rmsu['LE_D1_1km_5min'][7][1:] ) , -np.log( levels[1:] ) , 'r-')
p=ax.plot( 100*(profile_mean_rmsd['LE_D1_1km_30sec'][7][1:] / profile_mean_rmsu['LE_D1_1km_30sec'][7][1:] ) , -np.log( levels[1:] ) , 'r--' )


ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
my_levs = levels[[0,3,4,6,8]]
ytick=-np.log(my_levs)
ax.set_yticks(ytick)
#Get the level string list.
levels_str=[]
for ilev in my_levs  :
   levels_str.append( str(int(ilev)) )
ax.set_yticklabels(levels_str,fontsize=10,color='k')
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_xlim(left=0, right=50)
ax.text(1,-5.38,'(d) - $W$',fontsize=20,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})



#plt.show()
plt.savefig( figname , format='eps' , dpi=300)
plt.close()
