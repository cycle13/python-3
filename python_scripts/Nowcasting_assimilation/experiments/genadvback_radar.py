def genadvback_radar(nx, ny, u_motion1, v_motion1, ddt, code=-999):

    import numpy as np 
    import advection_back as advb
    from datetime import datetime
    from datetime import timedelta
    import os

    #PARAMETROS
    OBS = 'QC_ANGUIL'
    CASE = '20100111'
    DATA = '240'
    DBZ = 'QC_DATA'
    INPUT = 'motion_vectors_radar_prueba'
    DT = 'dt30'
    OUTPUT = 'advection_prueba'
    OUTCONF = 'advback_MSE_dt30'
    BASEDIR = '/home/aarruti'
    FILEDIR = BASEDIR + '/' + OBS + '/' + CASE + '/' + DATA + '/' + INPUT + '/' + DT
    BASEOUTDIR = '/home/aarruti'
    FILEOUTDIR = BASEOUTDIR + '/' + OBS + '/' + CASE + '/' + DATA + '/' + OUTPUT + '/' + OUTCONF
    if not os.path.exists(FILEOUTDIR):
        os.makedirs(FILEOUTDIR)

    #PARAMETROS
    year_start = 2010
    month_start = 01
    day_start = 11
    hour_start = 9
    minute_start = 00
    year_end = 2010
    month_end = 01
    day_end = 11
    hour_end = 20
    minute_end = 00
    formato_fecha = '%Y%m%d_%H%M'
    u_motions = u_motion1.reshape((nx,ny))
    v_motions = v_motion1.reshape((nx,ny))
    ddx = 2000
    ddt = 60 # Cada 60 segundos
    nmethod = 1 #Metodo para interpolar
    paso_t = 1 # Paso de tiempo de pronostico 1 minuto
    t_prono = 120 # Tiempo total de adveccion 2 horas
    dt_iniprono= 30 # Cada media hora se inicializa un pronostico
    UNDEF = -999

    inidate = datetime(year_start, month_start, day_start, hour_start, minute_start)
    enddate = datetime(year_end, month_end, day_end, hour_end, minute_end)

    vrange = [40]

    idate = inidate
    while idate <= enddate:
        idatename = idate.strftime(formato_fecha)
        path = FILEOUTDIR +'/'+ idatename
        if not os.path.exists(path):
            os.makedirs(path)
#    if not (os.path.exists(BASEDIR+'/'+ OBS +'/'+CASE+'/'+DATA+ '/' + DBZ + '/cappi_2km_resh2km_'+ idatename + '.npy') and os.path.exists(FILEDIR+'/motion_'+idatename+'_040_MSE       _uniform_promaroundmax.npz')):
        if not (os.path.exists(FILEDIR+'/motion_'+idatename+'_040_MSE_uniform_promaroundmax.npz')):
           print('No hay dato de reflectividad o de VM')
           idate = idate + timedelta(minutes=dt_iniprono)
           print idate
           continue
        ct1 = np.load(BASEDIR+'/'+ OBS +'/'+CASE+'/'+DATA+ '/' + DBZ + '/cappi_2km_resh2km_'+ idatename + '.npy')
#    print ct1
        for v in vrange:
#        path1 = path + '/%03d' %v
#        if not os.path.exists(path1):
#            os.makedirs(path1)
#        ct1 = np.load(BASEDIR+'/'+ OBS +'/'+CASE+'/'+DATA+ '/' + DBZ + '/cappi_2km_resh2km_'+ idatename + '.npy')
#        if not (os.path.exists(FILEDIR+'/motion_'+idatename+'_%03d_corrlineal_uniform_promaroundmax.npz' %v):            
#            idate = idate + timedelta(minutes=dt_iniprono)
#            continue
            data = np.load(FILEDIR+'/motion_'+idatename+'_%03d_MSE_uniform_promaroundmax.npz' %v)
            print(data.files)
            u_motion = data['u_motions']
            v_motion = data['v_motions']
            u_motion[u_motion<=UNDEF] = np.nan
            v_motion[v_motion<=UNDEF] = np.nan
            if np.isnan(u_motion).any() == True or np.isnan(v_motion).any() == True:
                print 'NaN'
            else:
                nx = nnx
                ny = nny
                dx = ddx
                dt_prono = ddt
                method = nmethod

                ct1[ct1 >= 1e10] = 0.0
                ct1[ct1 < 0.0] = 0.0
                campo = ct1

                iniforecast = idate + timedelta(minutes=paso_t)
                endforecast = iniforecast + timedelta(minutes=t_prono)

                iforecast = iniforecast
                iforecastname = iforecast.strftime(formato_fecha)

                ref_new = advb.advection_back(nx, ny, dx, dt_prono, u_motion, \
                v_motion, campo, method)

                np.save(FILEOUTDIR+'/'+idatename+'/'+iforecastname+'_%03d.npy' %v, ref_new)
                iforecast = iforecast+timedelta(minutes=paso_t)
                iforecastname = iforecast.strftime(formato_fecha)

                while iforecast <= endforecast:
                    campo = ref_new
                    ref_new = advb.advection_back(nx, ny, dx, dt_prono, u_motion, \
                    v_motion, campo, method)
                    np.save(FILEOUTDIR+'/'+idatename+'/'+iforecastname+'_%03d.npy' %v, ref_new)
                    iforecast = iforecast + timedelta(minutes=paso_t)
                    iforecastname = iforecast.strftime(formato_fecha)
        idate = idate + timedelta(minutes=dt_iniprono)






