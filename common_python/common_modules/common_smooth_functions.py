import numpy as np
from scipy import signal
from scipy import misc



def gaussian_smooth_2d( field , sigma , cutoff )   :
    #A gaussian filter in 2D that is applied independently to each
    #level if the field is 3D (levels are assumed to be the last dimension)
    field=np.squeeze(field)

    filter_size=np.around( 2*cutoff*sigma , 0 ).astype(int)

    filter_field=np.zeros( np.shape( field ) )

    if np.size( np.shape( field ) )  == 2   :
       nz = 1
    else                                    :
       nz = np.shape( field )[2] 

    gaussian_conv_func=np.zeros([ 2*filter_size+1 , 2*filter_size+1 ])
    for ii in range(0,2*filter_size+1) :
       for jj in range(0,2*filter_size+1) :
               gaussian_conv_func[ii,jj] = np.exp( -0.5*(np.power( ii-filter_size,2 ) + np.power( jj-filter_size, 2) ) /np.power(sigma,2)   )

    #Normalize the convolving function.
    gaussian_conv_func=gaussian_conv_func/np.sum( gaussian_conv_func )

    if nz > 1                  :
       for iz in range(0,nz) :
           filter_field[:,:,iz]=signal.convolve2d(np.squeeze(field[:,:,iz]),gaussian_conv_func, boundary='symm', mode='same')
    else                       :
           filter_field=signal.convolve2d(np.squeeze(field),gaussian_conv_func, boundary='symm', mode='same')[:,:]


    return  filter_field
