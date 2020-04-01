import numpy as np
import matplotlib.pyplot as plt

def Data_File_Name( Conf , date )  :

   #Conf : experiment configuration
   #Date : a datetime object containing the date corresponding to the file we want to read.
   date_str = date.strftime(Conf['format_date'])
   
   file_name = '../datos/cappi_2km_resh2km_'+ date_str +'.npy'

   return file_name 


def Plot_Data_Fast( Conf , field , nx , ny , figname=None ) :

    #Reshape and plot any field during the data assimilation process.

    field_2d= np.reshape( field , (nx,ny) ) 

    plt.figure()
    plt.pcolor( field_2d ) 
    plt.colorbar()
    if figname == None :
       plt.show() 
    else               :
       plt.savefig(figname)
    plt.close()
    

def Plot_Data_Fast_2d( Conf , field , figname=None )  :

    plt.figure()
    plt.pcolor( field )
    plt.colorbar()
    if figname == None :
       plt.show()
    else               :
       plt.savefig(figname)
    plt.close()

def Plot_Ref_Wind_Fast( Conf , Ref , Umv , Vmv , skip=5 , figname=None )  :

    Umv = np.reshape( Umv , np.shape( Ref ) )
    Vmv = np.reshape( Vmv , np.shape( Ref ) ) 

    Umv_plot = np.zeros( np.shape( Umv ) ) * np.nan
    Vmv_plot = np.zeros( np.shape( Vmv ) ) * np.nan
 
    Umv_plot[0::5,0::5]=Umv[0::5,0::5]
    Vmv_plot[0::5,0::5]=Vmv[0::5,0::5]

    plt.figure()
    plt.pcolor(Ref)
    plt.colorbar()
    plt.quiver(Umv_plot,Vmv_plot)

    mag= np.sqrt( np.power( Umv , 2 ) + np.power( Vmv , 2 ) )
    max_mag= np.max( mag )

    plt.title('Max mag = ' + str(max_mag) )
    if figname == None  :
       plt.show()
    else                :
       plt.savefig(figname)
    plt.close()




