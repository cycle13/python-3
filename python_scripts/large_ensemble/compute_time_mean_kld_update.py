# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

#Compute the temporal mean of a quantity in this case the mean.

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os

import common_plot_functions as cpf
import common_mask_functions as cmf
from common_functions import common_functions as comm 

file_name='kldistance'

basedir='/home/ra001011/a03471/data/output_data/'

expnames=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

deltat=[300,240,300,300,300,300]

init_date = ['20130713050500','20130713050400','20130713050500','20130713050500','20130713050500','20130713050500']

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

#Define initial and end times using datetime module.
etime = dt.datetime(2013,7,13,5,59,30)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iexp , my_exp in enumerate( expnames ) :

   delta = dt.timedelta(seconds=deltat[iexp])

   itime = dt.datetime.strptime( init_date[iexp] , '%Y%m%d%H%M%S' )

   #Compute the total number of times
   ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

   for my_file_type in filetypes     :

      ctl_file = basedir + '/' + my_exp + '/ctl/' + my_file_type + '.ctl'

      outputdir=basedir + '/' + my_exp + '/time_mean/' + my_file_type + '/'

      if not os.path.exists( outputdir)  :

         os.makedirs( outputdir )

      #=========================================================
      #  READ CTL FILE
      #=========================================================

      ctl_dict = ctlr.read_ctl( ctl_file )

      nx=ctl_dict['nx']
      ny=ctl_dict['nx']
      nlev=len( ctl_dict['full_lev_list'] )
      nt=int(1)             #Force the number of times to be one.
      ctl_dict['nt']=int(1) #Force the number of times to be one.

      undef=np.float32( ctl_dict['undef'] )

      if  ctl_dict['big_endian']   :
         dtypein = '>f4'
         endian='big_endian'
      else                         :
         dtypein = 'f4'
         endian='little_endian'

      if  ctl_dict['sequential']   :
         access='sequential'
      else                         :
         access='direct'

      sequential=ctl_dict['sequential']


      #=========================================================
      #  START TIME LOOP
      #=========================================================

      it=0
 
      #Initialize a temporal array
      tmp_parameter = np.zeros([nx,ny,nlev,ntimes])

      ctime = itime 

      while ( ctime <= etime )  :

         print( ctime )

         print ( 'Reading data')

         #=========================================================
         #  READ THE DATA
         #=========================================================

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/' + file_name +'.grd'

         #Read all the variables and levels at once
         tmp_parameter[:,:,:,it] = np.squeeze( ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) )

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

         gues = np.squeeze( ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) )

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/moment0001.grd'

         anal = np.squeeze( ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) )

         #Remove the grid points in which the prior has not been updated by the observations.
         tmp_parameter[:,:,:,it][ np.abs( gues - anal ) < 1.0e-7  ] = np.nan 

         print( np.nanmin(tmp_parameter[:,:,:,it]) , np.nanmax( tmp_parameter[:,:,:,it] ) )

         ctime = ctime + delta
 
         it = it + 1

      print ( "Finish time loop" )

      parameter = np.squeeze( np.nanmean( tmp_parameter , 3 ) )

      print( np.nanmin(parameter) , np.nanmax(parameter) )

      my_file=outputdir + '/' + file_name + 'update_mean' + '.grd'
      comm.write_data(outfile=my_file,mydata=parameter,nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)
     




