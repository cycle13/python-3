# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
#import matplotlib.pyplot as plt
import datetime as dt
import ctl_reader as ctlr
import binary_io  as bio
import os

from common_functions import common_functions as comm


basedir='/work/hp150019/f22003/output_data/'

#expnames  = ['LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min','LE_D1_1km_1min_4D','LE_D1_1km_5min']
expnames  = ['LE_D1_1km_5min_4D_OFP_V2']
#basedir='/dalab/jruiz/EXPERIMENTS_LARGE_ENSEMBLE/'

#expnames  = ['LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min','LE_D1_1km_1min_4D','LE_D1_1km_5min']
#expnames  = ['LE_D1_1km_30sec']


#expdeltas = [30,30,60,120,60,300]
expdeltas = [300]

#delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds

filetypes=['guesgp','analgp']   #analgp , analgz , guesgp , guesgz

smooth=False             #False- no smooth , True apply smooth
#smooth_lambda=10         #Smooth length scale (in number of grid points)

nbv=1000                 #Total number of ensemble members.

nmoments=4               #Cantidad total de momentos que vamos a calcular.

get_kldistance=True      #Wether we compute or not the Kullback-Leiber distance.

get_histogram=True       #Wether if histogram will be explicitelly calculated and stored.

get_moments=True         #Wether we will be computeing ensemble moments.

get_correlation=True      #Compute correlation between reflectivity, wind and other variables.

itime = dt.datetime(2013,7,13,5,0,0)  #Initial time.
etime = dt.datetime(2013,7,13,6,0,0)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iexp , my_exp_name in enumerate( expnames ) :

   delta=dt.timedelta( seconds = expdeltas[iexp] )

   for my_file_type in filetypes  :

      #=========================================================
      #  READ CTL FILE
      #=========================================================

      ctl_file = basedir + my_exp_name + '/ctl/' + my_file_type + '.ctl' 

      ctl_dict = ctlr.read_ctl( ctl_file )


      #Estimate number of bins in histogram. 
      nbins=int(np.round(np.sqrt(nbv))) #Number of bins to be used in histogram computation.


      nx=np.int(ctl_dict['nx'])
      ny=np.int(ctl_dict['ny'])
      nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
      nlev=ctl_dict['nz']                 #Number of vertical levels for 3D variables.

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

      #All the fields will be read so.
      n_selected_fields = nz
      selected_fields = np.arange(0,nz) + 1

      #=========================================================
      #  START COMPUTATION LOOP
      #=========================================================

      my_ensemble=np.zeros([nx,ny,nz,nbv]).astype('float32')
      #comm.allocate_ensemble(nx=nx,ny=ny,nz=nz,nbv=nbv) #Allocate the ensemble and the undefmask.

      ctime = itime

      while ( ctime <= etime ):

        print(ctime)

        print(' Reading the ensemble ')
 
        my_path=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/' 

        if not os.path.exists(my_path)                       :
           print('[Warning]:Path ' + my_path + ' does not exist: SKIP IT')

        else                                                 :
  
          #Read the ensemble and the undef mask
          [my_ensemble , my_undefmask] = comm.read_ensemble(path=my_path,nx=nx,ny=ny,nbv=nbv
                                                        ,selected_fields=selected_fields
                                                        ,n_selected_fields=n_selected_fields
                                                        ,undef=undef,ie=endian,acc=access)

          if get_moments :
             print("Computing pdf moments" + ctime.strftime("%Y%m%d%H%M%S") )
             #Get the moments of the PDF for different variables and levels.
             my_moments=comm.compute_moments(my_ensemble=my_ensemble,my_undefmask=my_undefmask 
                                        ,nx=nx,ny=ny,nz=nz,nbv=nbv,nmoments=nmoments,undef=undef) 

             for imoment in range (0,nmoments)  :
                momentstr="%04d" % ( imoment + 1 )
                my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/moment' + momentstr + '.grd'
                comm.write_data(outfile=my_file,mydata=my_moments[:,:,:,imoment],nx=nx,ny=ny,nz=nz,ie=endian,acc=access)

          if get_kldistance    :
             print("Computing kl divergence" + ctime.strftime("%Y%m%d%H%M%S") ) 
             kldist=comm.compute_kld(my_ensemble=my_ensemble,my_undefmask=my_undefmask
                                  ,nx=nx,ny=ny,nz=nz,nbv=nbv,undef=undef,nbins=nbins)
             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/kldistance.grd'
             comm.write_data(outfile=my_file,mydata=kldist,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)

          if get_correlation   :
             print("Computing correlation" + ctime.strftime("%Y%m%d%H%M%S") )
             #Each variable in var_list1 will be correlated with all variables in var_list2
             var_list1=['dbz','u','v']
             var_list2=['dbz','u','v','w','tk','qv']
             for my_var1 in var_list1 :
                 [var_si1,var_ei1] = ctlr.get_var_start_end( ctl_dict , my_var1 )
                 nzl=var_ei1 - var_si1 + 1
                 for my_var2 in  var_list2 :
                     if my_var2 != my_var1 :  #Skip correlation of a variable with itself
                        [var_si2,var_ei2] = ctlr.get_var_start_end( ctl_dict , my_var2 )
                        correlation=comm.compute_correlation(my_ensemble1=my_ensemble[:,:,var_si1:var_ei1+1],my_ensemble2=my_ensemble[:,:,var_si2:var_ei2+1]
                                  ,my_undefmask1=my_undefmask[:,:,var_si1:var_ei1+1],my_undefmask2=my_undefmask[:,:,var_si2:var_ei2+1]
                                  ,nx=nx,ny=ny,nz=nzl,nbv=nbv,undef=undef)
                        my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/correlation_'+my_var1+'_'+my_var2+'.grd'
                        comm.write_data(outfile=my_file,mydata=correlation,nx=nx,ny=ny,nz=nzl,ie=endian,acc=access)


          if get_histogram     :
             print("Computing histogram" + ctime.strftime("%Y%m%d%H%M%S") )
             #Compute histogram explicitelly
             [ensmin,ensmax,histogram]=comm.compute_histogram(my_ensemble=my_ensemble,my_undefmask=my_undefmask
                                                         ,nx=nx,ny=ny,nz=nz,nbv=nbv,nbins=nbins,undef=undef)
             #Write ensemble minimum
             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/ensmin.grd' 
             comm.write_data(outfile=my_file,mydata=ensmin,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)
             #Write ensemble maximum
             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/ensmax.grd'
             comm.write_data(outfile=my_file,mydata=ensmax,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)
             #Write histogram 
             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/histogram.grd'
             hist_out = np.zeros( (nx,ny,nz*nbins) )
             ii = 0
             for iz in range( 0 , nz )  :
                for ib in range( 0 , nbins )  :
                   hist_out[:,:,ii]=histogram[:,:,iz,ib]
                   ii = ii + 1

             bio.write_data_direct_woundef(my_file,hist_out,'i2',order='f') #Using low precission integer to store histogram.

             histogram_bimodality = comm.compute_histogram_bimodality_improved(histogram=histogram,my_undefmask=my_undefmask
                                                                    ,nx=nx,ny=ny,nz=nz,nbins=nbins,smoothw=2)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/bimodality_index_improved.grd'
             comm.write_data(outfile=my_file,mydata=histogram_bimodality,nx=nx,ny=ny,nz=nz,ie=endian,acc=access) 

        ctime = ctime + delta


print ( "Finish time loop" )



