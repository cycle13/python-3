import matplotlib.pyplot as plt 
import numpy as np
import wrf_module as wrf
import metpy.calc as mpcalc 
from metpy.cbook import get_test_data
from metpy.plots import SkewT , Hodograph
from metpy.units import units

dpi=300

def plot_momentum_equation_v( my_data , plot_path , show=False , force = False ) :
    #Ploteo un cross section vertical de la ecuacion de movimiento (una figura para w y otra para u y v)

    for it in range( my_data['nt'] ) :

       fig_name = plot_path + '/momentum_equation_w_' + my_data['slice_type'] + '_i_' + str(my_data['slice_index']) + '_t' + str(it) + '.png'
       if my_data['slice_type'] == 'vx' :
          dist=np.arange(0,my_data['dy']*my_data['ny'],my_data['dy'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,:,sw,it]
          t=my_data['t'][:,:,sw,it]
          qv=my_data['qv'][:,:,sw,it]
          w=my_data['w'][:,:,sw,it]
          hwind=my_data['v'][:,:,sw,it]
          ref=my_data['ref'][:,:,sw,it]
          bt=my_data['dwdt_bt'][:,:,sw,it]
          bqv=my_data['dwdt_bqv'][:,:,sw,it]
          bqc=my_data['dwdt_bqc'][:,:,sw,it]
          bp =my_data['dwdt_bp'][:,:,sw,it]
          dwdt =my_data['dwdt_loc'][:,:,sw,it]
          pz =my_data['dwdt_pz'][:,:,sw,it]
          adv =my_data['dwdt_adv'][:,:,sw,it]
          ppert = my_data['p'][:,:,sw,it] - my_data['p0'][:,:,sw,it]
       if my_data['slice_type'] == 'vy' :
          dist=np.arange(0,my_data['dx']*my_data['nx'],my_data['dx'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,sw,:,it]
          t=my_data['t'][:,sw,:,it]
          qv=my_data['qv'][:,sw,:,it]
          w=my_data['w'][:,sw,:,it]
          hwind=my_data['u'][:,sw,:,it]
          ref=my_data['ref'][:,sw,:,it]
          bt=my_data['dwdt_bt'][:,sw,:,it]
          bqv=my_data['dwdt_bqv'][:,sw,:,it]
          bqc=my_data['dwdt_bqc'][:,sw,:,it]
          bp =my_data['dwdt_bp'][:,sw,:,it]
          dwdt =my_data['dwdt_loc'][:,sw,:,it]
          pz =my_data['dwdt_pz'][:,sw,:,it]
          adv =my_data['dwdt_adv'][:,sw,:,it]
          ppert = my_data['p'][:,sw,:,it] - my_data['p0'][:,sw,:,it]
      
       ncols=4
       nrows=2
       if my_data['slice_type'] == 'vx' or my_data['slice_type'] == 'vy'  :
          ybound=[0,20000]
       fig, axs = plt.subplots( nrows,ncols , figsize=[15,9] , )
       fig.subplots_adjust(wspace=0.15,hspace=0.1,bottom=0.095,left=0.045,right=0.98,top=0.96)

       #Ploteo la reflectividad y la velocidad vertical.
       ax = axs[0,0]
       clevs1=np.arange(0,70,1)
       my_map = cmap_discretize('gist_ncar',clevs1.size)
       p1=ax.contourf( dist , z , ref , clevs1 , cmap=my_map)
       clevs2=np.array([1.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0])
       p2=ax.contour( dist , z , w , clevs2 , colors='k',linestyles='solid' )
       clevs2=np.array([-30.0,-20.0,-15.0,-10.0,-5.0,-1.0])
       p2=ax.contour( dist , z , w , clevs2 , colors='k',linestyles='dashed' )
       ax.set_ybound( ybound )
       ax.set_title('Reflectividad (sh, dBZ) y W (cont., $ms^{-1}$)')
       cbar_ax = fig.add_axes([0.06, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(ref)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)
       ax.set_xticks([])

       #Fijo los niveles para las componentes del empuje.
       clevs1=np.arange(-0.8,0.8+0.05,0.05)
       my_map = cmap_discretize('bwr',clevs1.size)

       #Ploteo el aporte de T al empuje y la temperatura
       ax = axs[0,1]
       p1=ax.contourf( dist , z , bt , clevs1 , cmap=my_map)
       clevs2=np.arange(100,320,5)
       p2=ax.contour( dist , z , t , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title('$B_T$ (sh., $ms^{-2}$) y T (cont., K)')
       ax.set_xticks([])
       ax.set_yticks([])

       cbar_ax = fig.add_axes([0.55, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(bt)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Ploteo el aporte de qv al empuje y qv.
       ax = axs[0,2]
       p1=ax.contourf( dist , z , bqv , clevs1 , cmap=my_map)
       clevs2=np.arange(0,18,1)
       p2=ax.contour( dist , z , qv * 1000 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title('$B_{q_v}$ (sh., $ms^{-2}$) y $q_v$ (cont., K)')
       ax.set_xticks([])
       ax.set_yticks([])
       
       #Ploteo el aporte de qc al empuje
       ax = axs[0,3]
       p1=ax.contourf( dist , z , bqc , clevs1 , cmap=my_map)
       clevs2=np.arange(0.1,1.0,0.1)
       p2=ax.contour( dist , z , bt + bp + bqv + bqc , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.arange(-1.0,0.0,0.1)
       p2=ax.contour( dist , z , bt + bp + bqv + bqc , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title('$B_{q_c}$ (sh., $ms^{-2}$) y Empuje total (cont., $ms^{-2}$) ')
       ax.set_xticks([])
       ax.set_yticks([])
       
       #Ploteo el aporte de P al empuje
       ax = axs[1,0]
       p1=ax.contourf( dist , z , bp , clevs1 , cmap=my_map) 
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       skipx=3
       skipz=3
       ax.quiver(dist[0::skipz,0::skipx],z[0::skipz,0::skipx],hwind[0::skipz,0::skipx],w[0::skipz,0::skipx])
       ax.set_ybound( ybound )
       ax.set_title('$B_P$ (sh., $ms^{-2}$) y viento(U,W) (vectores)')

       #Ploteo la componente vertical de la fuerza de presion y las perturbaciones de presion.
       ax = axs[1,1]
       p1=ax.contourf( dist , z , pz , clevs1 , cmap=my_map)
       clevs2=np.arange(25,500,50)
       p2=ax.contour( dist , z , ppert , clevs2 , colors='k',linestyles='solid' )
       clevs3=np.arange(-500,0,50)
       p3=ax.contour( dist , z , ppert , clevs3 , colors='k',linestyles='dashed' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{\frac{1}{{\rho}_0}}{{\nabla}_z}{P}^{\prime}$ (sh., $ms^{-2}$) y ${P}^{\prime}$ (cont., Pa)')
       ax.set_yticks([])
       
       #Ploteo el termino advectivo.
       ax = axs[1,2]
       p1=ax.contourf( dist , z , adv , clevs1 , cmap=my_map)
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{V}{\nabla}{w}$ (sh., ms{-2})')
       ax.set_yticks([])
       #Ploteo la aceleracion local y el residuo.
       ax = axs[1,3]
       p1=ax.contourf( dist , z , dwdt , clevs1 , cmap=my_map)
       clevs2=np.arange(0.05,0.5,0.2)
       p2=ax.contour( dist , z , dwdt - adv - pz - bt -bqv -bqc -bp , clevs2 , colors='k' , linestyles='solid' )
       clevs2=np.arange(-0.5,-0.05,0.2)
       p3=ax.contour( dist , z , dwdt - adv - pz - bt -bqv -bqc -bp , clevs2 , colors='k' , linestyles='dashed' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{{\partial}w}{{\partial}t}}}$ (sh., $ms^{-2}$), Residuo (cont. $ms^-2$)')
       ax.set_yticks([])

       if show :
          plt.show()

       plt.savefig(fig_name,dpi=dpi)
       plt.close()


    return

def plot_termo_equation_v( my_data , plot_path , show=False , force = False ) :
    #Ploteo un cross section vertical de la ecuacion de movimiento (una figura para w y otra para u y v)

    for it in range( my_data['nt'] ) :

       fig_name = plot_path + '/termo_equation_' + my_data['slice_type'] + '_i_' + str(my_data['slice_index']) + '_t' + str(it) + '.png'
       if my_data['slice_type'] == 'vx' :
          dist=np.arange(0,my_data['dy']*my_data['ny'],my_data['dy'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,:,sw,it]
          t=my_data['t'][:,:,sw,it]
          theta=my_data['theta'][:,:,sw,it]
          theta_pert=my_data['theta'][:,:,sw,it]-my_data['theta0'][:,:,sw,it]
          qv=my_data['qv'][:,:,sw,it]
          w=my_data['w'][:,:,sw,it]
          u=my_data['u'][:,:,sw,it]
          ref=my_data['ref'][:,:,sw,it]
          dthetadt_loc =my_data['dthetadt_loc'][:,:,sw,it]
          dthetadt_adv =my_data['dthetadt_adv'][:,:,sw,it]
          h_diabatic =my_data['h_diabatic'][:,:,sw,it]                   
          ppert = my_data['p'][:,:,sw,it] - my_data['p0'][:,:,sw,it]
          
       if my_data['slice_type'] == 'vy' :
          dist=np.arange(0,my_data['dx']*my_data['nx'],my_data['dx'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,sw,:,it]
          t=my_data['t'][:,sw,:,it]
          theta=my_data['theta'][:,sw,:,it]
          theta_pert=my_data['theta'][:,sw,:,it]-my_data['theta0'][:,sw,:,it]
          qv=my_data['qv'][:,sw,:,it]
          w=my_data['w'][:,sw,:,it]
          u=my_data['u'][:,sw,:,it]
          ref=my_data['ref'][:,sw,:,it]
          dthetadt_loc =my_data['dthetadt_loc'][:,sw,:,it]
          dthetadt_adv =my_data['dthetadt_adv'][:,sw,:,it]
          h_diabatic =my_data['h_diabatic'][:,sw,:,it]
          ppert = my_data['p'][:,sw,:,it] - my_data['p0'][:,sw,:,it]
          
       ncols=2
       nrows=2
       fig, axs = plt.subplots( nrows,ncols , figsize=[10,8] , )
       fig.subplots_adjust(wspace=0.04,hspace=0.09,bottom=0.1,left=0.06,right=0.98,top=0.96)
       if my_data['slice_type'] == 'vx' or my_data['slice_type'] == 'vy'  :
          ybound=[0,20000]
 
       #Ploteo theta y su perturbacion
       ax = axs[0,0]
       clevs1=np.arange(-15,15.5,0.5)
       theta_pert[theta_pert > np.max(clevs1)]=np.max(clevs1)
       theta_pert[theta_pert < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('RdBu_r',clevs1.size)
       p1=ax.contourf( dist , z , theta_pert , clevs1 , cmap=my_map)
       clevs2=np.arange(270,500,5)
       p2=ax.contour( dist , z , theta , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\theta}^{\prime}$ (sh., $K$) y $\theta$ (cont., $K$)')   
       ax.set_xticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cbar_ax = fig.add_axes([0.06, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(theta_pert)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Defino la escala de colores para la tasa de cambio de theta con el tiempo.
       clevs1=np.arange(-0.15,0.16,0.01) 

       #Ploteo la tasa de cambio local de theta con el tiempo
       ax = axs[0,1]
       dthetadt_loc[dthetadt_loc > np.max(clevs1)]=np.max(clevs1)
       dthetadt_loc[dthetadt_loc < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('bwr',clevs1.size)
       p1=ax.contourf( dist , z , dthetadt_loc , clevs1 , cmap=my_map)
       clevs2=np.arange(270,500,5)
       p2=ax.contour( dist , z , theta , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{\partial{\theta}}{{\partial}t}}$ (sh., $Ks^{-1}$) y $\theta$ (cont., $K$)')
       ax.set_xticks([])
       ax.set_yticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cbar_ax = fig.add_axes([0.55, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(dthetadt_loc)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Ploteo el aporte de qv al empuje y qv.
       ax = axs[1,0]
       my_map = cmap_discretize('bwr',clevs1.size)
       dthetadt_adv[dthetadt_adv > np.max(clevs1)]=np.max(clevs1)
       dthetadt_adv[dthetadt_adv < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dthetadt_adv , clevs1 , cmap=my_map)
       clevs2=np.arange(270,500,5)
       p2=ax.contour( dist , z , theta , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       skipx=3
       skipz=3
       ax.quiver(dist[0::skipz,0::skipx],z[0::skipz,0::skipx],u[0::skipz,0::skipx],w[0::skipz,0::skipx])
       ax.set_ybound( ybound )
       ax.set_title(r'$-V{\nabla}{\theta}$ (sh., $Ks^{-1}$) y $\theta$ (cont., $K$) y viento (U,W) (vectores)')
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)

       #Ploteo el aporte del calor latente. 
       #TODO reemplazar la estimacion del calor latente por el calor latente posta.
       ax = axs[1,1]
       my_map = cmap_discretize('bwr',clevs1.size)
       h_diabatic[h_diabatic > np.max(clevs1)]=np.max(clevs1)
       h_diabatic[h_diabatic < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , h_diabatic , clevs1 , cmap=my_map)
       clevs2=np.array([0.01,0.05,0.1,0.2])
       p2=ax.contour( dist , z , dthetadt_loc - dthetadt_adv - h_diabatic , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.array([-0.2,-0.1,-0.05,-0.01])
       p2=ax.contour( dist , z , dthetadt_loc - dthetadt_adv - h_diabatic , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$\dot{Q}_{lat}$ (sh., $Ks^{-1}$) y ${\frac{d{\theta}}{dt}}_{res}$ (cont., $K$)')
       ax.set_yticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       
       if show :
          plt.show()

       plt.savefig(fig_name,dpi=dpi)
       plt.close()


    return

def plot_vapor_equation_v( my_data , plot_path , show=False , force = False ) :
    #Ploteo un cross section vertical de la ecuacion de movimiento (una figura para w y otra para u y v)

    for it in range( my_data['nt'] ) :

       fig_name = plot_path + '/vapor_equation_' + my_data['slice_type'] + '_i_' + str(my_data['slice_index']) + '_t' + str(it) + '.png'
       if my_data['slice_type'] == 'vx' :
          dist=np.arange(0,my_data['dy']*my_data['ny'],my_data['dy'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,:,sw,it]
          qv=my_data['qv'][:,:,sw,it]
          qv_pert=my_data['qv'][:,:,sw,it]-my_data['qv0'][:,:,sw,it]
          qv=my_data['qv'][:,:,sw,it] 
          w=my_data['w'][:,:,sw,it]
          u=my_data['u'][:,:,sw,it]
          ref=my_data['ref'][:,:,sw,it]
          dqvdt_loc =my_data['dqvdt_loc'][:,:,sw,it]
          dqvdt_adv =my_data['dqvdt_adv'][:,:,sw,it]
          qv_diabatic =my_data['qv_diabatic'][:,:,sw,it]                   
          
       if my_data['slice_type'] == 'vy' :
          dist=np.arange(0,my_data['dx']*my_data['nx'],my_data['dx'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,sw,:,it]
          qv=my_data['qv'][:,sw,:,it]
          qv_pert=my_data['qv'][:,sw,:,it]-my_data['qv0'][:,sw,:,it]
          w=my_data['w'][:,sw,:,it]
          u=my_data['u'][:,sw,:,it]
          ref=my_data['ref'][:,sw,:,it]
          dqvdt_loc =my_data['dqvdt_loc'][:,sw,:,it]
          dqvdt_adv =my_data['dqvdt_adv'][:,sw,:,it]
          qv_diabatic =my_data['qv_diabatic'][:,sw,:,it]
          
       ncols=2
       nrows=2
       fig, axs = plt.subplots( nrows,ncols , figsize=[10,8] , )
       fig.subplots_adjust(wspace=0.04,hspace=0.09,bottom=0.1,left=0.06,right=0.98,top=0.96)
       if my_data['slice_type'] == 'vx' or my_data['slice_type'] == 'vy'  :
          ybound=[0,20000]
 
       #Ploteo theta y su perturbacion
       ax = axs[0,0]
       clevs1=np.arange(-5.0,5.05,0.05)
       qv_pert[qv_pert > np.max(clevs1)]=np.max(clevs1)
       qv_pert[qv_pert < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('RdBu_r',clevs1.size)
       p1=ax.contourf( dist , z , qv_pert * 1000 , clevs1 , cmap=my_map)
       clevs2=np.arange(0,18,1)
       p2=ax.contour( dist , z , qv * 1.0e3 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${q_v}^{\prime}$ (sh., $gKg^{-1}$) y $q_v$ (cont., $gKg^{-1}$)')   
       ax.set_xticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cbar_ax = fig.add_axes([0.06, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(qv_pert)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Defino la escala de colores para la tasa de cambio de theta con el tiempo.
       clevs1=np.arange(-0.05,0.055,0.005) 

       #Ploteo la tasa de cambio local de theta con el tiempo
       ax = axs[0,1]
       dqvdt_loc[dqvdt_loc > np.max(clevs1)]=np.max(clevs1)
       dqvdt_loc[dqvdt_loc < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('bwr',clevs1.size)
       p1=ax.contourf( dist , z , dqvdt_loc * 1.0e3 , clevs1 , cmap=my_map)
       clevs2=np.arange(0,18,1)
       p2=ax.contour( dist , z , qv * 1.0e3 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{\partial{q_v}}{{\partial}t}}$ (sh., $gKg^{-1}s^{-1}$) y $q_v$ (cont., $gKg^{-1}$)')
       ax.set_xticks([])
       ax.set_yticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cbar_ax = fig.add_axes([0.55, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(dqvdt_loc)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Ploteo el aporte de qv al empuje y qv.
       ax = axs[1,0]
       my_map = cmap_discretize('bwr',clevs1.size)
       dqvdt_adv[dqvdt_adv > np.max(clevs1)]=np.max(clevs1)
       dqvdt_adv[dqvdt_adv < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqvdt_adv * 1.0e3 , clevs1 , cmap=my_map)
       clevs2=np.arange(0,18,1)
       p2=ax.contour( dist , z , qv * 1.0e3 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       skipx=3
       skipz=3
       ax.quiver(dist[0::skipz,0::skipx],z[0::skipz,0::skipx],u[0::skipz,0::skipx],w[0::skipz,0::skipx])
       ax.set_ybound( ybound )
       ax.set_title(r'$-V{\nabla}{q_v}$ (sh., $gKg^{-1}$) y $q_v$ (cont., $gKg^{-1}s^{-1}$) y viento (U,W) (vectores)')
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)

       #Ploteo el aporte del calor latente. 
       #TODO reemplazar la estimacion del calor latente por el calor latente posta.
       ax = axs[1,1]
       my_map = cmap_discretize('bwr',clevs1.size)
       qv_diabatic[qv_diabatic > np.max(clevs1)]=np.max(clevs1)
       qv_diabatic[qv_diabatic < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , qv_diabatic * 1.0e3 , clevs1 , cmap=my_map)
       clevs2=np.array([0.005,0.01,0.05,0.1])
       p2=ax.contour( dist , z , ( dqvdt_loc - dqvdt_adv - qv_diabatic ) * 1.0e3 , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.array([-0.1,-0.05,-0.01,-0.005])
       p2=ax.contour( dist , z , ( dqvdt_loc - dqvdt_adv - qv_diabatic ) * 1.0e3 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$cond.$ (sh., $gKg^{-1}s^{-1}$) y ${\frac{d{q_v}}{dt}}_{res}$ (cont., $gKg^{-1}$)')
       ax.set_yticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       
       if show :
          plt.show()

       plt.savefig(fig_name,dpi=150)
       plt.close()

    return


def plot_water_equation_v( my_data , plot_path , show=False , force = False ) :
    #Ploteo un cross section vertical de la ecuacion de movimiento (una figura para w y otra para u y v)

    for it in range( my_data['nt'] ) :

       fig_name = plot_path + '/water_equation_' + my_data['slice_type'] + '_i_' + str(my_data['slice_index']) + '_t' + str(it) + '.png'
       if my_data['slice_type'] == 'vx' :
          dist=np.arange(0,my_data['dy']*my_data['ny'],my_data['dy'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,:,sw,it]
          ref=my_data['ref'][:,:,sw,it]
          qv=my_data['qv'][:,:,sw,it] * 1000.0
          qc=my_data['qc'][:,:,sw,it] * 1000.0
          qr=my_data['qr'][:,:,sw,it] * 1000.0
          qi=my_data['qi'][:,:,sw,it] * 1000.0
          qs=my_data['qs'][:,:,sw,it] * 1000.0
          qg=my_data['qg'][:,:,sw,it] * 1000.0
          dqvdt=( my_data['dqvdt_loc'][:,:,sw,it] - my_data['dqvdt_adv'][:,:,sw,it] ) * 1000.0
          dqcdt=( my_data['dqcdt_loc'][:,:,sw,it] - my_data['dqcdt_adv'][:,:,sw,it] ) * 1000.0
          dqrdt=( my_data['dqrdt_loc'][:,:,sw,it] - my_data['dqrdt_adv'][:,:,sw,it] ) * 1000.0
          dqidt=( my_data['dqidt_loc'][:,:,sw,it] - my_data['dqidt_adv'][:,:,sw,it] ) * 1000.0
          dqsdt=( my_data['dqsdt_loc'][:,:,sw,it] - my_data['dqsdt_adv'][:,:,sw,it] ) * 1000.0
          dqgdt=( my_data['dqgdt_loc'][:,:,sw,it] - my_data['dqgdt_adv'][:,:,sw,it] ) * 1000.0

       if my_data['slice_type'] == 'vy' :
          dist=np.arange(0,my_data['dx']*my_data['nx'],my_data['dx'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,sw,:,it]
          ref=my_data['ref'][:,sw,:,it]
          qv=my_data['qv'][:,sw,:,it] * 1000.0
          qc=my_data['qc'][:,sw,:,it] * 1000.0
          qr=my_data['qr'][:,sw,:,it] * 1000.0
          qi=my_data['qi'][:,sw,:,it] * 1000.0
          qs=my_data['qs'][:,sw,:,it] * 1000.0
          qg=my_data['qg'][:,sw,:,it] * 1000.0
          dqvdt=( my_data['dqvdt_loc'][:,sw,:,it] - my_data['dqvdt_adv'][:,sw,:,it] ) * 1000.0
          dqcdt=( my_data['dqcdt_loc'][:,sw,:,it] - my_data['dqcdt_adv'][:,sw,:,it] ) * 1000.0
          dqrdt=( my_data['dqrdt_loc'][:,sw,:,it] - my_data['dqrdt_adv'][:,sw,:,it] ) * 1000.0
          dqidt=( my_data['dqidt_loc'][:,sw,:,it] - my_data['dqidt_adv'][:,sw,:,it] ) * 1000.0
          dqsdt=( my_data['dqsdt_loc'][:,sw,:,it] - my_data['dqsdt_adv'][:,sw,:,it] ) * 1000.0
          dqgdt=( my_data['dqgdt_loc'][:,sw,:,it] - my_data['dqgdt_adv'][:,sw,:,it] ) * 1000.0

      
       ncols=3
       nrows=2
       
       if my_data['slice_type'] == 'vx' or my_data['slice_type'] == 'vy'  :
          ybound=[0,20000]
          
       #Fijo los niveles para las tasas de cambio local
       clevs1=np.arange(-0.1,0.1+0.01,0.01)
       my_map = cmap_discretize('bwr',clevs1.size)   
       clevs2=np.array([0.5,1.0,2.5,5.0,10.0,15.0,20.0,25.0])
       #clevs2=np.arange(0,18,0.5)
          
       fig, axs = plt.subplots( nrows,ncols , figsize=[15,9] , )
       fig.subplots_adjust(wspace=0.15,hspace=0.1,bottom=0.095,left=0.045,right=0.98,top=0.96)

       #grafico los cambios en q_v y q_v
       ax = axs[0,0]
       my_map = cmap_discretize('bwr',clevs1.size)
       dqvdt[dqvdt > np.max(clevs1)]=np.max(clevs1)
       dqvdt[dqvdt < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqvdt , clevs1 , cmap=my_map )
       p2=ax.contour( dist , z , qv , clevs2 , colors='k',linestyles='solid' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{dq_v}{dt}}$ (sh, $gKg^{-1}s{-1}$) y $q_v$ (cont., $gKg^{-1}$)')
       cbar_ax = fig.add_axes([0.06, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(dqvdt)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)
       ax.set_xticks([])

       #Grafico los cambios en q_cloud y q_cloud
       ax = axs[0,1]
       dqcdt[dqcdt > np.max(clevs1)]=np.max(clevs1)
       dqcdt[dqcdt < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqcdt , clevs1 , cmap=my_map )
       p2=ax.contour( dist , z , qc , clevs2 , colors='k',linestyles='solid' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{dq_c}{dt}}$ (sh, $gKg^{-1}s{-1}$) y $q_c$ (cont., $gKg^{-1}$)')
       ax.set_xticks([])
       ax.set_yticks([])

       #Grafico los cambios en q_rain y q_rain
       ax = axs[0,2]
       dqrdt[dqrdt > np.max(clevs1)]=np.max(clevs1)
       dqrdt[dqrdt < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqrdt , clevs1 , cmap=my_map )
       p2=ax.contour( dist , z , qr , clevs2 , colors='k',linestyles='solid' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{dq_r}{dt}}$ (sh, $gKg^{-1}s{-1}$) y $q_r$ (cont., $gKg^{-1}$)')
       ax.set_xticks([])
       ax.set_yticks([])
            
       #Grafico los cambios en q_ice y q_ice
       ax = axs[1,0]
       dqidt[dqidt > np.max(clevs1)]=np.max(clevs1)
       dqidt[dqidt < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqidt , clevs1 , cmap=my_map )
       p2=ax.contour( dist , z , qi , clevs2 , colors='k',linestyles='solid' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{dq_i}{dt}}$ (sh, $gKg^{-1}s{-1}$) y $q_i$ (cont., $gKg^{-1}$)')
              
       #Grafico los cambios en q_snow y q_snow
       ax = axs[1,1]
       dqsdt[dqsdt > np.max(clevs1)]=np.max(clevs1)
       dqsdt[dqsdt < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqsdt , clevs1 , cmap=my_map )
       p2=ax.contour( dist , z , qs , clevs2 , colors='k',linestyles='solid' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{dq_s}{dt}}$ (sh, $gKg^{-1}s{-1}$) y $q_s$ (cont., $gKg^{-1}$)')
       ax.set_yticks([])
       
       #Grafico los cambios en q_graup y q_graup
       ax = axs[1,2]
       dqgdt[dqgdt > np.max(clevs1)]=np.max(clevs1)
       dqgdt[dqgdt < np.min(clevs1)]=np.min(clevs1)
       p1=ax.contourf( dist , z , dqgdt , clevs1 , cmap=my_map )
       p2=ax.contour( dist , z , qg , clevs2 , colors='k',linestyles='solid' )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${\frac{dq_g}{dt}}$ (sh, $gKg^{-1}s{-1}$) y $q_g$ (cont., $gKg^{-1}$)')
       ax.set_yticks([])
       
       if show :
          plt.show()

       plt.savefig(fig_name,dpi=dpi)
       plt.close()


    return

def plot_ppert_equation_v( my_data , plot_path , show=False , force = False ) :
    #Ploteo un cross section vertical de la ecuacion de movimiento (una figura para w y otra para u y v)

    for it in range( my_data['nt'] ) :

       fig_name = plot_path + '/ppert_equation_' + my_data['slice_type'] + '_i_' + str(my_data['slice_index']) + '_t' + str(it) + '.png'
       if my_data['slice_type'] == 'vx' :
          dist=np.arange(0,my_data['dy']*my_data['ny'],my_data['dy'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,:,sw,it]
          bouy=my_data['bouy'][:,:,sw,it]
          w=my_data['w'][:,:,sw,it]
          hwind=my_data['v'][:,:,sw,it]
          hfppert=my_data['fpperty'][:,:,sw,it]
          vfppert=my_data['fppertz'][:,:,sw,it]
          ref=my_data['ref'][:,:,sw,it]
          lp_ppert=my_data['lp_ppert'][:,:,sw,it]
          lp_spin_tot=my_data['lp_spin_tot'][:,:,sw,it]
          lp_splat_tot=my_data['lp_splat_tot'][:,:,sw,it]
          lp_spin_pert =my_data['lp_spin_pert'][:,:,sw,it]
          lp_splat_pert =my_data['lp_splat_pert'][:,:,sw,it]
          lp_lineal =my_data['lp_lineal'][:,:,sw,it]
          lp_b =my_data['lp_b'][:,:,sw,it]
          ppert = my_data['p'][:,:,sw,it] - my_data['p0'][:,:,sw,it]
          ppert_hy = my_data['ppert_hy'][:,:,sw,it]
          ppert_nhy= my_data['ppert_nhy'][:,:,sw,it]
       if my_data['slice_type'] == 'vy' :
          dist=np.arange(0,my_data['dx']*my_data['nx'],my_data['dx'])
          dist=np.tile(dist,(my_data['nz'],1))
          sw=my_data['slice_width']
          z=my_data['z'][:,sw,:,it]
          bouy=my_data['bouy'][:,sw,:,it]
          w=my_data['w'][:,sw,:,it]
          hwind=my_data['u'][:,sw,:,it]
          hfppert=my_data['fppertx'][:,sw,:,it]
          vfppert=my_data['fppertz'][:,sw,:,it]
          ref=my_data['ref'][:,sw,:,it]
          lp_ppert=my_data['lp_ppert'][:,sw,:,it]
          lp_spin_tot=my_data['lp_spin_tot'][:,sw,:,it]
          lp_splat_tot=my_data['lp_splat_tot'][:,sw,:,it]
          lp_spin_pert =my_data['lp_spin_pert'][:,sw,:,it]
          lp_splat_pert =my_data['lp_splat_pert'][:,sw,:,it]
          lp_lineal =my_data['lp_lineal'][:,sw,:,it]
          lp_b =my_data['lp_b'][:,sw,:,it]
          ppert = my_data['p'][:,sw,:,it] - my_data['p0'][:,sw,:,it]
          ppert_hy = my_data['ppert_hy'][:,sw,:,it]
          ppert_nhy= my_data['ppert_nhy'][:,sw,:,it]
      
       ncols=4
       nrows=2
       if my_data['slice_type'] == 'vx' or my_data['slice_type'] == 'vy'  :
          ybound=[0,20000]
       fig, axs = plt.subplots( nrows,ncols , figsize=[15,9] , )
       fig.subplots_adjust(wspace=0.15,hspace=0.1,bottom=0.095,left=0.045,right=0.98,top=0.96)

       #Ploteo la perturbacion de presion y la fuerza de presion
       ax = axs[0,0]
       clevs1=np.arange(-600,650,50)
       ppert[ppert > np.max(clevs1)]=np.max(clevs1)
       ppert[ppert < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('RdBu_r',clevs1.size)
       p1=ax.contourf( dist , z , ppert , clevs1 , cmap=my_map)
       skipx=3
       skipz=3
       ax.quiver(dist[0::skipz,0::skipx],z[0::skipz,0::skipx],hfppert[0::skipz,0::skipx],vfppert[0::skipz,0::skipx],scale=10)

       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'${P}^{\prime}$ (sh., $Pa$) y $-{\frac{1}{{\rho}_0}}{\nabla}{P}^{\prime}$ (vectores)',fontsize=10)
       ax.set_xticks([])
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cbar_ax = fig.add_axes([0.06, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(ppert)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Ploteo la perturbacion de presion no hidrostatica
       ax = axs[0,1]
       ppert_hy[ppert_hy > np.max(clevs1)]=np.max(clevs1)
       ppert_hy[ppert_hy < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('RdBu_r',clevs1.size)
       p1=ax.contourf( dist , z , ppert_hy , clevs1 , cmap=my_map)
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       skipx=3
       skipz=3
       ax.quiver(dist[0::skipz,0::skipx],z[0::skipz,0::skipx],hwind[0::skipz,0::skipx],w[0::skipz,0::skipx])
       ax.set_ybound( ybound )
       ax.set_title(r'${P}^{\prime}_{hydro}$ (sh., $Pa$) y viento (vectores)',fontsize=10)
       ax.set_xticks([])
       ax.set_yticks([])

       #Ploteo la perturbacion de presion hidrostatica
       ax = axs[0,2]
       ppert_nhy[ppert_nhy > np.max(clevs1)]=np.max(clevs1)
       ppert_nhy[ppert_nhy < np.min(clevs1)]=np.min(clevs1)
       my_map = cmap_discretize('RdBu_r',clevs1.size)
       p1=ax.contourf( dist , z , ppert_nhy , clevs1 , cmap=my_map)
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       skipx=3
       skipz=3
       ax.quiver(dist[0::skipz,0::skipx],z[0::skipz,0::skipx],hwind[0::skipz,0::skipx],w[0::skipz,0::skipx])

       ax.set_ybound( ybound )
       ax.set_title(r'${P}^{\prime}_{nohydro}$ (sh., $Pa$) y viento (vectores) ',fontsize=10)
       ax.set_xticks([])
       ax.set_yticks([])

       #Fijo los niveles para las componentes del empuje.
       clevs1=np.arange(-8.0,8.0+0.05,0.05) 
       my_map = cmap_discretize('bwr',clevs1.size)

       #Ploteo laplaciano de la perturbacion de presion por empuje y empuje.
       ax = axs[0,3]
       lp_b[lp_b > np.max(clevs1)]=np.max(clevs1) 
       lp_b[lp_b < np.min(clevs1)]=np.min(clevs1) 
       p1=ax.contourf( dist , z , -lp_b * 1.0e4 , clevs1 , cmap=my_map)
       clevs2=np.arange(0.1,1.0,0.1)*10.0e-1
       p2=ax.contour( dist , z , bouy , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.arange(-1.0,0.0,0.1)*10.0e-1
       p2=ax.contour( dist , z , bouy , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )

       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{\nabla}^{2}{P}^{\prime}_b$ (sh., $10^{-4} ms^{-2}$) y B (cont., $10^{-1}ms^{-2}$)',fontsize=10)
       ax.set_xticks([])
       ax.set_yticks([])

       cbar_ax = fig.add_axes([0.55, 0.03, 0.4, 0.02])
       m = plt.cm.ScalarMappable(cmap=my_map )
       m.set_array(-lp_b)
       m.set_clim(np.min(clevs1),np.max(clevs1))
       delta= ( np.max(clevs1)-np.min(clevs1) )/ (clevs1.size-1)
       cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(np.min(clevs1),np.max(clevs1)+delta,delta),orientation='horizontal')
       cb.ax.tick_params(labelsize=14)

       #Ploteo laplaciano de la perturbacion de presion dinamica 
       ax = axs[1,0]
       p1=ax.contourf( dist , z , (-lp_spin_tot -lp_splat_tot) * 1.0e4 , clevs1 , cmap=my_map)
       clevs2=np.arange(1.0,11,1.0)
       p2=ax.contour( dist , z , (-lp_spin_pert -lp_splat_pert) * 1.0e4 , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.arange(-10.0,0.0,1.0)
       p2=ax.contour( dist , z , (-lp_spin_pert -lp_splat_pert) * 1.0e4 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{\nabla}^{2}{P}^{\prime}_{dyn}$ pert (sh. $10^{-4}ms^{-2}$) y tot (cont) ($ms^{-2}$)',fontsize=10)
       ax.set_ybound( ybound )
       
       #Ploteo laplaciano de la perturbacion por spin total y asociada al viento perturbado
       ax = axs[1,1]
       p1=ax.contourf( dist , z , -lp_spin_pert * 1.0e4 , clevs1 , cmap=my_map)
       clevs2=np.arange(1.0,11,1.0)
       p2=ax.contour( dist , z , -lp_spin_tot * 1.0e4 , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.arange(-10.0,0.0,1.0)
       p2=ax.contour( dist , z , -lp_spin_tot * 1.0e4 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{\nabla}^{2}{P}^{\prime}_{spin}$ pert (sh.) y tot (cont) ($10^{-4}ms^{-2}$)',fontsize=10)
       ax.set_yticks([])
       
       #Ploteo laplaciano de la perturbacion por splat total y asociada al viento perturbado
       ax = axs[1,2]
       p1=ax.contourf( dist , z , -lp_splat_pert * 1.0e4 , clevs1 , cmap=my_map) 
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       clevs2=np.arange(1.0,11,1.0)
       p2=ax.contour( dist , z , -lp_splat_tot * 1.0e4 , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )
       clevs2=np.arange(-10.0,0.0,1.0)
       p2=ax.contour( dist , z , -lp_splat_tot * 1.0e4 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{\nabla}^{2}{P}^{\prime}_{splat}$ pert (sh.) y tot (cont) ($10^{-4}ms^{-2}$)',fontsize=10)
       ax.set_yticks([])

       #Ploteo el termino lineal.
       ax = axs[1,3]
       p1=ax.contourf( dist , z , -lp_lineal  * 1.0e4, clevs1 , cmap=my_map)
       clevs2=np.arange(-7.0,0.0,1.0)
       p2=ax.contour( dist , z , -lp_ppert * 1.0e4 , clevs2 , colors='k',linestyles='dashed' , linewidths=0.5 )
       clevs2=np.arange(1.0,8.0,1.0)
       p2=ax.contour( dist , z , -lp_ppert * 1.0e4 , clevs2 , colors='k',linestyles='solid' , linewidths=0.5 )

       p3=ax.contour( dist , z , ref , [30.0] , colors='c',linestyles='solid' , linewidths=2.5 )
       ax.set_ybound( ybound )
       ax.set_title(r'$-{\nabla}^{2}{P}^{\prime}_{lineal}$ (sh. $10^{-4}ms^{-2}$) y $-{\nabla}^{2}{P}^{\prime}$ (cont. $10^{-4}ms^{-2}$)',fontsize=10)
       ax.set_yticks([])

       if show :
          plt.show()

       plt.savefig(fig_name,dpi=dpi)
       plt.close()


    return

#Modified from http://wiki.scipy.org/Cookbook/Matplotlib/ColormapTransformations

def plot_sounding( my_data , plot_path , show=False )  :

    fig_name = plot_path + '/profile' + '_i_' + str(my_data['xp']) + '_j_' + str(my_data['yp']) + '_t' + str(my_data['tp']) + '.png'
    #Ejemplo tomado de https://unidata.github.io/MetPy/latest/examples/plots/Skew-T_Layout.html#sphx-glr-examples-plots-skew-t-layout-py
    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(9, 7))
    #add_metpy_logo(fig, 630, 80, size='large')

    # Grid for plots
    #rotation controla la inclinacion de las isotermas. Si es 0 el grafico se asemeja a un emagrama pero 
    #si es 45 entonces se obtiene un skew-T.
    skew = SkewT(fig, rotation=0)

    p=my_data['p'].data * units.hPa
    t=my_data['t'].data * units.degC
    td=my_data['td'] * units.degC
    u=my_data['u'] * units.meter / units.second
    v=my_data['v'] * units.meter / units.second

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.plot(p , t , 'r')
    skew.plot(p , td , 'g')
    skew.plot_barbs(p,u,v)
    skew.ax.set_ylim(1000, 100)

    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    # Calculate LCL height and plot as black dot. Because `p`'s first value is
    # ~1000 mb and its last value is ~250 mb, the `0` index is selected for
    # `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
    # i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
    # should be selected.
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0] , t[0] , td[0]  )
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

    # Calculate full parcel profile and add to plot as black line
    prof = mpcalc.parcel_profile(p , t[0] , td[0] ).to('degC')
    skew.plot(p, prof, 'k', linewidth=2)

    # Shade areas of CAPE and CIN
    skew.shade_cin(p,t,prof)
    skew.shade_cape(p,t,prof)
    # Good bounds for aspect ratio
    skew.ax.set_xlim(-70, 40)

    # Create a hodograph
    ax = fig.add_axes([0.72,0.66,0.2,0.2])
    h = Hodograph(ax, component_range=60.)
    h.add_grid(increment=20)
    h.plot(my_data['u'],my_data['v'])

    if show :
       plt.show()

    plt.savefig(fig_name,dpi=dpi)
    plt.close()


def cmap_map(function,cmap):

    import matplotlib
    import matplotlib.cm     as cm
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker


    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = [x[0] for x in cdict[key]]
    step_list = sum(list(step_dict.values()), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map( reduced_cmap, step_list)))
    new_LUT = np.array(list(map( function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  [x + (x[1], ) for x in list(this_cdict.items())]
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def cmap_xmap(function,cmap):
    """ Applies function, on the indices of colormap cmap. Beware, function
    should map the [0, 1] segment to itself, or you are in for surprises.

    See also cmap_xmap.
    """

    import matplotlib
    import matplotlib.cm     as cm
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    cdict = cmap._segmentdata
    function_to_map = lambda x : (function(x[0]), x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = list(map(function_to_map, cdict[key]))
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), "Resulting indices extend out of the [0, 1] segment."


    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    import matplotlib
    import matplotlib.cm     as cm
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


