# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""


import numpy as np
import matplotlib as plt
import datetime as dt
import binary_io as bio
import bred_vector_functions as bvf
import os 

basedir='/home/jruiz/share/scale_breeding/tmp/'

expname = '/breeding_osaka_pawr_1km/out/'

plotbasedir=basedir + expname + '/plots/'

inibv=1          #Initial bred vector to plot.
endbv=1          #Final bred vector to plot.
iniit=1          #Initial iter to plot.
endit=1          #End iter to plot.

plotlevels=np.array([0,6,13,17])   #Which levels will be plotted.
plotvars=['U','V','W','T','QV']    #Which variables will be plotted.
#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,10,30)  #Initial time.
etime = dt.datetime(2013,7,13,5,31,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

nx=180
ny=180
nz=20

data_pp_o=dict()
data_pn_o=dict()

data_pp_r=dict()
data_pn_r=dict()

bv_o=dict()
bv_r=dict()

ctime=itime + delta

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[0,:,:]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[0,:,:]


nit=endit-iniit+1
nbv=endbv-inibv+1

ntimes= np.round( (itime - etime).seconds / delta.seconds )

#Allocate the bv_norm array.
bv_norm=np.zeros( nbv , ntimes , nit )

for ibv in range (inibv , endbv + 1):

   bvstr="%04d" % ibv

   print( ' Plotting bred vector number ' + bvstr )

   while ( ctime <= etime ):

    for iter in range ( iniit , endit + 1 ):

      iterstr="%04d" % iter
 
      ptime=ctime - delta #Data correspinding to the previous step (to compute bv growth)
 
      print ( 'The date is :', ctime )

      print ( 'Reading the norm')

      myfile=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + bvstr + '/' + 'norm' + iterstr + '/.grd'

      norm=bio.read_norm_scale(myfile,nx,ny,nz,'>f4')


      ctime = ctime + delta

print ( "Finish time loop" )





