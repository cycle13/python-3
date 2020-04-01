import numpy as np
import os 

figconf=dict()


def set_default()   :

  #Reset configuration values to default values.

  figconf['figure']=True
  figconf['show']=False
  figconf['close']=True
  figconf['save']=True
  figconf['dpi']=300
  figconf['figtype']='png'
  figconf['figpath']='./'
  figconf['figname']='my_figure'

  figconf['figsize']=(12,10)

  #Title
  figconf['addtitle']=True
  figconf['title']=''
  figconf['titlefontsize']=20

  #Labels
  figconf['labelfontsize']=20
  figconf['labelfontcolor']='black'
  figconf['xlabel']=[]
  figconf['ylabel']=[]

  #Coastline
  figconf['coastline']=True
  figconf['coastlinewidth']=1
  figconf['coastlineres']='10m'

  #Cartopy
  figconf['projection']='PlateCarree'

  #Plots
  #Filled
  figconf['pcolor']=False
  figconf['contourf']=False
  figconf['shadedmin']=[]
  figconf['shadedmax']=[]
  figconf['shadedcolormap']='Blues'
  figconf['colorbar']=True
  figconf['colorbarfontsize']=15
  figconf['axessize']=[0.1,0.1,0.8,0.8]

  figconf['contourflevels']=[]
  figconf['cflabel']=False
  

  #Cont
  figconf['contour']=False
  figconf['contourlevels']=[]
  figconf['ncontourlevels']=10
  figconf['contourcolormap']=None
  figconf['contourcolor']='k'
  figconf['clabel']=True
  figconf['clabelfont']=12

  #Vectors
  figconf['vector']=False

  figconf['axislabels']=True
  figconf['gridline']=True
  figconf['gridlinewidth']=1
  figconf['gridlinecolor']='Gray'
  figconf['gridlinestyle']='--'

  #Axis ticks
  figconf['xtick']=[]
  figconf['ytick']=[]
  figconf['yticklabel']=[]
  figconf['xticklabel']=[]

  #Axis range
  figconf['axesrange']=[]

  #Lineplots
  figconf['linestyle']=['-']
  figconf['linemarker']=[]
  figconf['linecolor']=['k']

  #Subplot
  figconf['nsubplotrows']=0
  figconf['nsubplotcolumns']=0
  figconf['subplotvmargin']=0.045
  figconf['subplothmargin']=0.045
  figconf['subplothoffset']=0.0
  figconf['subplotvoffset']=0.0

  #Text annotation
  figconf['text']=False
  figconf['textstring']=[]
  figconf['textlocx']=0.5
  figconf['textlocy']=0.5
  figconf['textfontsize']=[]
  figconf['textha']='center'
  figconf['textfontsize']=20
  figconf['textcolor']='k'

 

def  compute_default_limits( x , y , varsh , varc , varv )  :

    if figconf['contour']    :
       if np.size( figconf['contourlevels'] ) == 0  :
           vmax=np.nanmax(varc)
           vmin=np.nanmin(varc)
           figconf['contourlevles']=np.arange(vmin,vmax,(vmax-vmin)/figconf['ncontourlevels'])

    if figconf['pcolor'] or figconf['contourf']  :
       if np.size( figconf['shadedmax'] ) == 0      :
           figconf['shadedmax'] = np.nanmax( varsh )
       if np.size( figconf['shadedmin'] ) == 0      :
           figconf['shadedmin'] = np.nanmax( varsh )


def  plot_x_y_subplot( x , y , varsh , varc , varv , cartopy=False )   :
     import matplotlib
     import matplotlib.cm     as cm
     import matplotlib.pyplot as plt
     import matplotlib.ticker as mticker

     from mpl_toolkits.basemap import Basemap
     import cartopy.crs as ccrs

     if not os.path.exists(figconf['figpath']):
         os.makedirs(figconf['figpath'])


     figconf['figure'] = False

     showfig = figconf['show'] 
     figconf['show'] = False
     figclose = figconf['close']
     figconf['close']  = False

     fig=plt.figure(1,figsize=figconf['figsize'])

     irow=figconf['nsubplotrows']
     icol=1
     icoldelta=1.0/figconf['nsubplotcolumns']
     irowdelta=1.0/figconf['nsubplotrows']

     vmargin=figconf['subplotvmargin']
     hmargin=figconf['subplothmargin']
   
     hoffset=figconf['subplothoffset']
     voffset=figconf['subplotvoffset']


     axessize=figconf['axessize']

     savefig=figconf['save']
     figconf['save']=False

     for isubplot in range( 0 , figconf['nsubplotrows'] * figconf['nsubplotcolumns'] )  :

         sub_varsh=np.squeeze(varsh[:,:,isubplot])
         sub_varc =np.squeeze(varc[:,:,isubplot])
         sub_varv =np.squeeze(varv[:,:,isubplot])
         sub_varc=np.ma.masked_where(np.isnan(sub_varc),sub_varc)
         sub_varsh=np.ma.masked_where(np.isnan(sub_varsh),sub_varsh)
         compute_default_limits( x , y , sub_varsh , sub_varc , sub_varv )

         figconf['axessize']=[icoldelta*(icol-1)+hmargin+hoffset,irowdelta*(irow-1)+vmargin+voffset,irowdelta-2*hmargin,icoldelta-2*vmargin]
 
  
         if figconf['addtitle']   :
            figconf['title'] = figconf['subplottitles'][isubplot]

         tmp_figconf = figconf.copy()

         if type( figconf['colorbar'] ) == list  :
            tmp_figconf['colorbar'] = figconf['colorbar'][isubplot]
         if type( figconf['shadedcolormap'] ) == list :
            tmp_figconf['shadedcolormap'] = figconf['shadedcolormap'][isubplot]
         if type( figconf['shadedmin'] ) == list :
            tmp_figconf['shadedmin'] = figconf['shadedmin'][isubplot]
         if type( figconf['shadedmax'] ) == list :
            tmp_figconf['shadedmax'] = figconf['shadedmax'][isubplot]


         if  cartopy   :
             plot_x_y_cartopy( x , y , sub_varsh , sub_varc , sub_varv , my_figconf=tmp_figconf) 
         else          :
             plot_x_y( x , y , sub_varsh , sub_varc , sub_varv , my_figconf=tmp_figconf )

         icol = icol + 1
         if icol > figconf['nsubplotcolumns']   :
            irow = irow - 1
            icol = 1

     if figconf['text']    :
        plt.figtext( figconf['textlocx'],figconf['textlocy'],figconf['textstring']
            , fontsize=figconf['textfontsize'] , ha=figconf['textha'] 
            , color = figconf['textcolor'] ,bbox={'facecolor':'red', 'alpha':0.5, 'pad':5})

     if savefig   :
         plt.savefig( figconf['figpath'] + figconf['figname'] + '.' + figconf['figtype'] , dpi=figconf['dpi'])
         figconf['save']=True

     if showfig    :
         plt.show()
         figconf['show']=True

     if figclose   :
         plt.close()
         figconf['close']=True

     figconf['axessize']=axessize
  
def  plot_x_y_cartopy( lon , lat , varsh , varc , varv , my_figconf )   :

     import matplotlib
     import matplotlib.cm     as cm
     import matplotlib.pyplot as plt
     import matplotlib.ticker as mticker

     from mpl_toolkits.basemap import Basemap
     import cartopy.crs as ccrs

     if my_figconf == None :
        my_figconf = figconf 
       
     #varsh -> var for shading plot
     #varc  -> var for contour plot
     #varv  -> var for vector plot 
 
     #Estimate some quantities.
     varsh=np.squeeze(varsh)
     varc =np.squeeze(varc)
     varv =np.squeeze(varv)

     varc=np.ma.masked_where(np.isnan(varc),varc)
     varsh=np.ma.masked_where(np.isnan(varsh),varsh)

     if not os.path.exists(my_figconf['figpath']):
         os.makedirs(my_figconf['figpath'])


     compute_default_limits( lon , lat , varsh , varc , varv )


     if my_figconf['figure']  :

          fig=plt.figure(1,figsize=my_figconf['figsize'])


     ax = plt.axes(my_figconf['axessize'],projection=ccrs.PlateCarree())


     if my_figconf['pcolor']   :
          p=ax.pcolor(lon, lat, np.squeeze( varsh ) ,
              transform=ccrs.PlateCarree(),vmin=my_figconf['shadedmin'] , 
              vmax=my_figconf['shadedmax'],cmap=my_figconf['shadedcolormap'] )

     if my_figconf['contourf']  :
          p=ax.contourf(lon, lat, np.squeeze( varsh ) ,
              transform=ccrs.PlateCarree(),vmin=my_figconf['shadedmin'] ,
              vmax=my_figconf['shadedmax'],cmap=my_figconf['shadedcolormap'] )

     if my_figconf['colorbar'] and ( my_figconf['pcolor'] or my_figconf['contourf']  )  :
          cb=plt.colorbar(p)
          cb.ax.tick_params(labelsize=my_figconf['colorbarfontsize'])

     if my_figconf['contour']  :
         if my_figconf['contourcolormap'] == None       :
            c=ax.contour( lon , lat , np.squeeze( varc ) ,
              levels = my_figconf['contourlevels'] , colors=my_figconf['contourcolor'])
         else                                       :
            c=ax.contour( lon , lat , np.squeeze( varc ) , 
              levels = my_figconf['contourlevels'] , cmap=my_figconf['contourcolormap'])

     if my_figconf['coastline']  :
         ax.coastlines(my_figconf['coastlineres'],linewidth=my_figconf['coastlinewidth'])

     if my_figconf['gridline']   :
         gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=my_figconf['axislabels'],
                linewidth=my_figconf['gridlinewidth'], color=my_figconf['gridlinecolor'], 
                alpha=0.5, linestyle=my_figconf['gridlinestyle'])

         if np.size(my_figconf['xtick']) > 0  :
             gl.xlocator=mticker.FixedLocator(my_figconf['xtick'])
         if np.size(my_figconf['ytick']) > 0  :
             gl.ylocator=mticker.FixedLocator(my_figconf['ytick'])

         gl.xlabel_style = {'size': my_figconf['labelfontsize'], 'color': my_figconf['labelfontcolor']}
         gl.ylabel_style = {'size': my_figconf['labelfontsize'], 'color': my_figconf['labelfontcolor']}
         gl.xlabels_top = False
         gl.ylabels_right = False

     if np.size( my_figconf['axesrange'] ) > 0   :

         ax.set_extent(my_figconf['axesrange'], ccrs.PlateCarree())


     if my_figconf['addtitle']   :
         plt.title(my_figconf['title'] , fontsize = my_figconf['titlefontsize'])
     if np.size(my_figconf['xlabel']) > 0  :
         plt.xlabel(my_figconf['xlabel'])
     if np.size(my_figconf['ylabel']) > 0  :
         plt.ylabel(my_figconf['ylabel'])



     if my_figconf['save']   :
         plt.savefig( my_figconf['figpath'] + my_figconf['figname'] + '.' + my_figconf['figtype'] , dpi=my_figconf['dpi'])


     if my_figconf['show']    :
         plt.show()

     if my_figconf['close']   :
         plt.close()


def  plot_x_y( x , y , varsh , varc , varv , my_figconf = None )   :
     #varsh -> var for shading plot
     #varc  -> var for contour plot
     #varv  -> var for vector plot 
     import matplotlib.pyplot as plt

     if my_figconf == None :
        my_figconf = figconf 
 
     #Estimate some quantities.
     varsh=np.squeeze(varsh)
     varc =np.squeeze(varc)
     varv =np.squeeze(varv)

     varc=np.ma.masked_where(np.isnan(varc),varc)
     varsh=np.ma.masked_where(np.isnan(varsh),varsh)

     print(my_figconf['shadedcolormap'])

     compute_default_limits( x , y , varsh , varc , varv )

     if not os.path.exists(my_figconf['figpath']):
         os.makedirs(my_figconf['figpath'])

     if my_figconf['figure']  :

          fig=plt.figure(1,figsize=my_figconf['figsize'])


     ax = plt.axes(my_figconf['axessize'])

     if my_figconf['pcolor']   :
          p=plt.pcolor(x, y, np.squeeze( varsh ) 
             ,vmin=my_figconf['shadedmin'] , 
              vmax=my_figconf['shadedmax'],cmap=my_figconf['shadedcolormap'] )

          p.set_clim( my_figconf['shadedmin'],my_figconf['shadedmax'] )
     if my_figconf['contourf']  :
          if np.size( my_figconf['contourflevels'] == 0 )  :
             p=plt.contourf(x, y, np.squeeze( varsh ) , 
             levels=np.linspace(my_figconf['shadedmin'],my_figconf['shadedmax'],my_figconf['ncontourlevels'])
             #,vmin=my_figconf['shadedmin'] ,
             # vmax=my_figconf['shadedmax']
             ,cmap=my_figconf['shadedcolormap'] )
          
          else                                          :
             p=plt.contourf(x, y, np.squeeze( varsh )
             ,levels=my_figconf['contourflevels'] )
          if my_figconf['cflabel']    :
             plt.clabel(c, inline=1, fontsize=my_figconf['clabelfont'])


          p.set_clim( my_figconf['shadedmin'],my_figconf['shadedmax'] )

     if my_figconf['colorbar'] and ( my_figconf['pcolor'] or my_figconf['contourf']  )  :
          cb=plt.colorbar(p)
          cb.ax.tick_params(labelsize=my_figconf['colorbarfontsize'])

     if my_figconf['contour']  :
         c=ax.contour( x , y , np.squeeze( varc ) , 
              levels = my_figconf['contourlevels'] , cmap=my_figconf['contourcolormap'])
         if my_figconf['clabel']    :
            plt.clabel(c, inline=1, fontsize=my_figconf['clabelfont'])

     if my_figconf['addtitle']   :
         plt.title(my_figconf['title'] , fontsize = my_figconf['titlefontsize'])
     if np.size(my_figconf['xlabel']) > 0  :
         plt.xlabel(my_figconf['xlabel'])
     if np.size(my_figconf['ylabel']) > 0  :
         plt.ylabel(my_figconf['ylabel'])


     if my_figconf['gridline']   :
         ax.grid(color=my_figconf['gridlinecolor'], linestyle=my_figconf['gridlinestyle']
               , linewidth=my_figconf['gridlinewidth'] )

     if np.size(my_figconf['xtick']) > 0  :
         ax.set_xticks(my_figconf['xtick'])
         plt.xticks(fontsize=my_figconf['labelfontsize'])

     if np.size(my_figconf['ytick']) > 0  :
         ax.set_yticks(my_figconf['ytick'])
         plt.yticks(fontsize=my_figconf['labelfontsize'])


     if np.size(my_figconf['xticklabel']) > 0 :
         ax.set_xticklabels(my_figconf['xticklabel'],fontsize=my_figconf['labelfontsize'])
     if np.size(my_figconf['yticklabel']) > 0 :
         ax.set_yticklabels(my_figconf['yticklabel'],fontsize=my_figconf['labelfontsize'])

     if np.size(my_figconf['xlabel']) > 0  :
         ax.set_xlabel(my_figconf['xlabel'],fontsize=my_figconf['labelfontsize'])
     if np.size(my_figconf['ylabel']) > 0  :
         ax.set_ylabel(my_figconf['ylabel'],fontsize=my_figconf['labelfontsize'])

     
     if np.size( my_figconf['axesrange'] ) > 0   :
         plt.ylim( my_figconf['axesrange'][2] , my_figconf['axesrange'][3] )
         plt.xlim( my_figconf['axesrange'][0] , my_figconf['axesrange'][1] )

     if my_figconf['save']   :
         plt.savefig( my_figconf['figpath'] + my_figconf['figname'] + '.' + my_figconf['figtype'] , dpi=my_figconf['dpi'])

     if my_figconf['show']    :
         plt.show()

     if my_figconf['close']   :
         plt.close()


def  plot_lines( x , y )   :

      import matplotlib.pyplot as plt

      #Line plot
      x=np.squeeze(x)
      y=np.squeeze(y)


      if not os.path.exists(figconf['figpath']):
         os.makedirs(figconf['figpath'])

      get_line_prop( x , y )

      if figconf['figure']  :

         fig=plt.figure(1,figsize=figconf['figsize'])

      ax = plt.axes(figconf['axessize'])

      if np.size( y.shape ) > 1  :
         nc = y.shape[1]
      else                       : 
         nc = 1

      if nc == 1         :
        p=ax.plot(x,y,figconf['lineprop'][0] )
      else               :
        for iv in range( 0 , nc )                  :
           p=ax.plot(x[:,iv],y[:,iv],figconf['lineprop'][iv] )


      if figconf['gridline']   :
           ax.grid(color=figconf['gridlinecolor'], linestyle=figconf['gridlinestyle']
                 , linewidth=figconf['gridlinewidth'] )

      if np.size(figconf['xtick']) > 0  :
           ax.set_xticks(figconf['xtick'])
           tmpl=ax.get_xticklabels(minor=False, which=None)
           ax.set_xticklabels(tmpl,fontsize=figconf['labelfontsize'])
      if np.size(figconf['ytick']) > 0  :
           ax.set_yticks(figconf['ytick'])
           tmpl=ax.get_yticklabels(minor=False, which=None)
           ax.set_yticklabels(tmpl,fontsize=figconf['labelfontsize'])

      if np.size(figconf['xlabel']) > 0  :
           ax.set_xlabel(figconf['xlabel'][0],fontsize=figconf['labelfontsize'])
      if np.size(figconf['ylabel']) > 0  :
           ax.set_ylabel(figconf['ylabel'][0],fontsize=figconf['labelfontsize'])

      if figconf['save']   :
          plt.savefig( figconf['figpath'] + figconf['figname'] + '.' + figconf['figtype'] , dpi=figconf['dpi'])

      if figconf['show']    :
          plt.show()

      if figconf['close']   :
          plt.close()


def  get_line_prop( x , y )      :

     #Obtain the line characteristics for a particular plot.
     #If multiple line properties are specified then combine them
    
     if np.size( y.shape ) > 1    : 
        nc = y.shape[1]  #Get number of variables to plot.
     else                         :
        nc = 1

     figconf['lineprop']=list()

     for iv in range( 0 , nc )  :

       tmplp=list()

       for prop in ['linestyle','linemarker','linecolor']  :

          if np.size( figconf[prop] ) == nc         : 
             tmplp.append(figconf[prop][iv])
          else                                      :
             tmplp.append(figconf[prop][0])      

       figconf['lineprop'].append(''.join(tmplp[:]))
       

def plot_bar( x , y )     :
      #To plot a histogram or a bar plot in general.
      import matplotlib.pyplot as plt



 
 


#Modified from http://wiki.scipy.org/Cookbook/Matplotlib/ColormapTransformations
   
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

