"""
import numpy as np
import matplotlib.pyplot as plt

UObs = np.zeros((240,240)).flatten()
VObs = np.zeros((240,240)).flatten()
sigma = 5
nx = 240
ny = 240
nens = 100
desv_gauss = 15
"""

def gen_pertfield_ini(nx, ny, UObs, VObs, sigma, desv_gauss, nens):
    import numpy as np 
    import scipy as sp
    from scipy import ndimage
    pert_field = np.zeros(((UObs.size,nens,2))) 

    for i in range(nens):
        pert_field[:,i,0] = sp.ndimage.filters.gaussian_filter(np.random.normal(loc=0.0, scale=sigma, size=(nx,ny)),desv_gauss, mode='constant').flatten()
        pert_field[:,i,1] = sp.ndimage.filters.gaussian_filter(np.random.normal(loc=0.0, scale=sigma, size=(nx,ny)),desv_gauss, mode='constant').flatten()

    #Forzamos a que la media de las perturbaciones sea 0 y a que el sigma sea igual al prescripto.
    pert_field_mean=np.mean( pert_field , 1 )
    pert_field_sigma = np.std( pert_field , 1 )

    for i in range(nens) :
        pert_field[:,i,:] = ( pert_field[:,i,:] - pert_field_mean ) * sigma / pert_field_sigma
        pert_field[:,i,0] = pert_field[:,i,0] + UObs
        pert_field[:,i,1] = pert_field[:,i,1] + VObs

    return pert_field



def gen_pertfield(nx,ny,XF,sigma,desv_gauss,nens):
    import numpy as np
    import scipy as sp
    from scipy import ndimage
    pert_field = np.zeros(((XF[:,0,0].size,nens,2)))

    for i in range(nens):
        pert_field[:,i,0] = sp.ndimage.filters.gaussian_filter(np.random.normal(loc=0.0, scale=sigma, size=(nx,ny)),desv_gauss, mode='constant').flatten()
        pert_field[:,i,1] = sp.ndimage.filters.gaussian_filter(np.random.normal(loc=0.0, scale=sigma, size=(nx,ny)),desv_gauss, mode='constant').flatten()

    #Forzamos a que la media de las perturbaciones sea 0 y a que el sigma sea igual al prescripto.
    pert_field_mean=np.mean( pert_field , 1 )
    pert_field_sigma = np.std( pert_field , 1 )

    for i in range(nens) :
        pert_field[:,i,:] = ( pert_field[:,i,:] - pert_field_mean ) * sigma / pert_field_sigma

        pert_field[:,i,0] = pert_field[:,i,0] + XF[:,i,0]
        pert_field[:,i,1] = pert_field[:,i,1] + XF[:,i,1]

    return pert_field
    
"""
pert_field, pert_smooth_field = gen_pertfield_ini(nx,ny,UObs,VObs,sigma,desv_gauss,nens)
plt.pcolormesh(np.reshape(pert_field[:,0,0],(nx,ny)), cmap='seismic')#, vmin=-0.5, vmax=0.5)
plt.colorbar()
plt.show()
plt.pcolormesh(np.reshape(pert_smooth_field[:,0,1],(nx,ny)), cmap='seismic')#, vmin=-0.5, vmax=0.5)
plt.colorbar()
plt.show()
plt.pcolormesh(np.reshape(pert_smooth_field[:,99,1],(nx,ny)), cmap='seismic')#, vmin=-0.5, vmax=0.5)
plt.colorbar()
plt.show()
"""

