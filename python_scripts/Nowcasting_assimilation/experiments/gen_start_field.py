def gen_start_field(field_start_0, field_start_1, fileoutdir, date_start, dt, dx, desp_max, sigma, sigmatrend, nx, ny, box_size, UNDEF=-999):
    import numpy as np
    from motion_vectors import motion_vectors as wmv
    #import w_motion_vectors_old as wmv

    #reload(wmv)

    FILEOUTDIR = fileoutdir

    field_t0 = field_start_0
    field_t1 = field_start_1

    r = field_t0.shape[0]/2
    H, W = field_t0.shape
    x, y = np.meshgrid(np.arange(H), np.arange(W))
    xc = 120
    yc = 120
    d = (x-xc)**2 + (y-yc)**2
    mask = 1-(d < r**2)
    field_t0[mask == 1] = UNDEF
    field_t1[mask == 1] = UNDEF

    field_t0[np.isnan(field_t0)] = UNDEF
    field_t1[np.isnan(field_t1)] = UNDEF
    field_t0[field_t0 > 1e3] = UNDEF
    field_t1[field_t1 > 1e3] = UNDEF
    field_t0[field_t0 < 5] = UNDEF
    field_t1[field_t1 < 5] = UNDEF

    print("start to calculate start field")
    min_box_fraction = (5*(2*box_size)*(2*box_size))/100
    u_motion, v_motion, correlation, reftrend, nref, local_correlation = \
    wmv.motion_vector(field_t0=field_t0, field_t1=field_t1, dt=dt, dx=dx, \
    box_size=box_size, sigma=sigma, min_box_fraction=min_box_fraction, \
    nx=nx, ny=ny, desp_max=desp_max, aux_inputi=130, aux_inputj=130, \
    motion_vector_option=1, motion_vector_weigth=0, \
    motion_vector_norm=2)

    #Remove outliers from u_motion and v_motion using a local mean filter.
    print("start the QA: remove outliers")
    u_motion1, v_motion1 = wmv.filter_outlier(fieldu=u_motion,fieldv=v_motion, nx=nx, ny=ny,threshold=10,box_size=60)

    u_motion_mean = wmv.undef_mean(u_motion1, nx, ny)
    v_motion_mean = wmv.undef_mean(v_motion1, nx, ny)

    #Replace undef u_motion and v_motion by the field mean
    print("Complete de field")
    u_motion2 = np.empty((nx, ny))
    v_motion2 = np.empty((nx, ny))
    u_motion2[:, :] = u_motion1
    v_motion2[:, :] = v_motion1
    u_motion2[u_motion2==UNDEF] = u_motion_mean
    v_motion2[v_motion2==UNDEF] = v_motion_mean

    #Smooth u_motion and v_motion using a gaussian filter
    print("Start the gaussian filter")
    u_motions = wmv.gaussian_filter(field0=u_motion2, dx=dx, sigma=sigma, \
    nx=nx, ny=ny)
    v_motions = wmv.gaussian_filter(field0=v_motion2, dx=dx, sigma=sigma, \
    nx=nx, ny=ny)

    print("Save the fields")
    np.savez(FILEOUTDIR+'/motion_' + \
    date_start.strftime('%Y%m%d%H%M%S') + '_040_MSE_uniform_promaroundmax.npz', \
    u_motion=u_motion, v_motion=v_motion, u_motion1=u_motion1, \
    v_motion1=v_motion1, u_motion2=u_motion2, v_motion2=v_motion2, u_motions=u_motions, v_motions=v_motions, \
    local_correlation=local_correlation, reftrend=reftrend)
    print("FINISH START FIELD")
    return u_motions.flatten(), v_motions.flatten()
