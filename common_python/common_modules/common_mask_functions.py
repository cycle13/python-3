import numpy as np
import warnings        #To supress certain warnings.

import sys
sys.path.append('../common_functions/')
sys.path.append('../common_modules/')

from common_functions import common_functions as comm



def lat_lon_to_i_j(lonfield,latfield,lonlist,latlist) :
#Gets the i,j which is closer to a particular lat and lon give a latfield, lonfield.
   npoints=latlist.size

   i=np.zeros(latlist.shape)
   j=np.zeros(latlist.shape)
 
   for ipoint in range(0,npoints) :
 
     dist=np.power( latfield - latlist[ipoint] , 2.0 ) + np.power( lonfield - lonlist[ipoint] , 2.0 )

     #Get the indices of the minimum
     i[ipoint] , j[ipoint] = np.unravel_index(dist.argmin(), dist.shape)

   return i , j 


def get_regional_average_grid( var , xi , xe , yi , ye, zi, ze  , undef ) :
   #Given a 3D variable and a set of regions, get the mean, max and min of the variable over the grid points
   #contained in these regions.
   #Region is defined as a 3D rectangular section of the grid (in grid space). 

   #If input is scalar, then convert it to a list 
   if np.size(xi) == 1 :
      xi = np.array([xi]) ; xe = np.array([xe]) ; yi = np.array([yi]) ; ye = np.array([ye]) ; zi = np.array([zi]) ; ze = np.array([ze])


   nx = np.shape( var )[0]
   ny = np.shape( var )[1]
   if ( np.size ( np.shape ( var ) ) >= 3 ) :
      nz = np.shape( var )[2]
   else :
      nz = 1

   nregs=xi.size  #Get number of regions
  
   var_mean=np.zeros(nregs)
   var_max=np.zeros(nregs)
   var_min=np.zeros(nregs)

   for ireg in range(0,nregs) :
 
      if ( nz > 1 ) :
       
        tmp=var[int(xi[ireg]):int(xe[ireg])+1,
                int(yi[ireg]):int(ye[ireg])+1,
                int(zi[ireg]):int(ze[ireg])+1] 

        tmp[ tmp == undef ] = np.nan

      else : 
 
        tmp=var[int(xi[ireg]):int(xe[ireg])+1,
                int(yi[ireg]):int(ye[ireg])+1]

        tmp[ tmp == undef ] = np.nan

      #Np.nan.. functions produce warnings with the input consist only of nan  values. This warning is ignored. 
      with warnings.catch_warnings()    :
        warnings.simplefilter("ignore", category=RuntimeWarning)

        var_mean[ireg]=np.nanmean( tmp )

        var_max[ireg]=np.nanmax( tmp )
  
        var_min[ireg]=np.nanmin( tmp )

   return var_mean , var_max , var_min 


def distance_range_mask( lon_o , lat_o , max_range , lon , lat )   :
  #Given a lat,lon 2d field, and an origin lon_o , lat_o compute a mask
  #which is true if the grid point distance to the origin is less than max_range (meters)
  #and false otherwise. 

  my_mask=np.zeros( lon.shape , dtype=bool )

  
  for ii in range( 0 , my_mask.shape[0] ) :
     for jj in range( 0 , my_mask.shape[1] ) :

         dist=comm.com_distll_sngl( lon_o,lat_o,lon[ii,jj],lat[ii,jj] ) 

         if ( dist < max_range )  :
            my_mask[ii,jj]=True


  return my_mask

def box_mask( lon_i , lon_e , lat_i , lat_e , lon , lat )   :
   #Give a box defined by 2 lons and 2 lats (rectangular box in the lat lon domain)
   #Return a mask indicating which grid points are within the max.
   #The max is true is the point is within the box and False otherwise.

   return  np.logical_and( np.logical_and( lat >= lat_i , lat <= lat_e ) , np.logical_and( lon >= lon_i , lon <= lon_e ) ) 











    

