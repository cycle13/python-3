#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

#This model contains routines for ploting diagnostics related to the different
#filters applied in the QC routine.

#Function names are formed as follows: plot_[filter_name] 

def plot_Dealiasing( qc_output , options )  :

    import numpy as np
    import matplotlib.pyplot as plt
  

    filter_name='Dealiasing'

    for ilev in options['plot']['Elevs']  :

       plt.figure(figsize=(8, 8))
       plt.subplot(2,2,1)

       tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
       tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

       plt.pcolor(qc_output['x'][:,:,ielev]/1e3,qc_output['y'][:,:,ielev]/1e3,tmp_cv[:,:,ielev],vmin=ptions['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
       plt.title('Corrected Doppler Velocity')
       plt.colorbar()

       plt.subplot(2,2,2)
       plt.pcolor(qc_output['x'][:,:,ielev]/1e3,qc_output['y'][:,:,ielev]/1e3,tmp_v[:,:,ielev],vmin=ptions['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
       plt.title('Original Doppler Velocity')
       plt.colorbar()

       #plt.subplot(2,2,3)
       #plt.pcolor(qc_output['x'][:,:,ielev]/1e3,qc_output['y'][:,:,ielev]/1e3, ( qc_output['qcv'][:,:,ielev]==options[filter_name]['texture_code'] ).astype(float) )
       #plt.title('Pixels eliminated by texture texture filter')
       #plt.colorbar()

       plt.subplot(2,2,4)
       plt.pcolor(qc_output['x'][:,:,ielev]/1e3,qc_output['y'][:,:,ielev]/1e3, ( qc_output['qcv'][:,:,ielev]==options[filter_name]['code'] ).astype(float) )
       plt.title('Pixels where aliasing was corrected')
       plt.colorbar()

       if show  :

           plt.show()

       figname=options['plot']['FigNamePrefix'] + filter_name + options['plot']['FigNameSufix']
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return



def plot_RhoFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='RhoFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)
                            
    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by RHO Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['rho_smooth'][:,:,elev] ),vmin=0,vmax=1.1 )
    plt.title('Smooth RHO HV')
    plt.colorbar()
                                                                                                    
    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return



def plot_EchoTopFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='EchoTopFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Echo Top Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['echo_top'][:,:,elev] ),vmin=0,vmax=15000 )
    plt.title('Echo Top')
    plt.colorbar()
													
    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_EchoDepthFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='EchoDepthFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Echo Depth Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['echo_top'][:,:,elev] ),vmin=0,vmax=15000 )
    plt.title('Echo Depth')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return

def plot_RefSpeckleFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='RefSpeckleFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Ref. Speckle Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['speckle_ref'][:,:,elev] ),vmin=0,vmax=1 )
    plt.title('Speckle Index')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_DopplerSpeckleFilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='DopplerSpeckleFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, qc_output['speckle_v'][:,:,elev] , vmin=0,vmax=1)
    plt.title('Speckle Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels eliminated by Doppler Speckle Filter')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return



def plot_DopplerTextureFilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='DopplerTextureFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, qc_output['texture_v'][:,:,elev] , vmin=0 , vmax=50 )
    plt.title('Texture Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels eliminated by Texture Filter')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return


def plot_DopplerLocalStdFilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='DopplerLocalStdFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )
    tmp_std=np.ma.masked_array( qc_output['v_std'] , qc_output['v_std'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_std[:,:,elev]  )
    plt.title('Local Std')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels eliminated by Texture Filter')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return

def plot_DopplerSpatialCoherenceFilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='DopplerSpatialCoherenceFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, qc_output['coherence_index'][:,:,elev] )
    plt.title('Coherence Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by DSCF')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return



def plot_DopplerNoiseFilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt

    filter_name='DopplerNoiseFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    tmp_tmp=np.ma.masked_array( qc_output['distance_1'] , qc_output['distance_1'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_tmp[:,:,elev] )
    plt.title('Noise index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels eliminated by Noise Filter')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return



def plot_ReflectivityTextureFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='ReflectivityTextureFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['texture_ref'][:,:,elev] ) )
    plt.title('Texture Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Ref. Texture Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return

def plot_AttenuationFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='AttenuationFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['attenuation'][:,:,elev] ) )
    plt.title('Attenuation')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Attenuation Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return

def plot_BlockingFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='BlockingFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['blocking'][:,:,elev] ) )
    plt.title('Blocking')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Blocking Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_LowElevFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='LowElevFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['smooth_ref'][:,:,elev] ) )
    plt.title('Smooth Ref')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Low Elev. Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_LowDopplerFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='LowDopplerFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['v'][:,:,elev] ) , vmin=-30.0 , vmax=30.0, cmap='pyart_NWSVel' )
    plt.title('Radial Velocity')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Low Doppler Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_InterferenceFilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=False)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='InterferenceFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,  tmp_cref[:,:,elev] - tmp_ref[:,:,elev] )
    plt.title('Reflectivity difference')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Interference Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return





