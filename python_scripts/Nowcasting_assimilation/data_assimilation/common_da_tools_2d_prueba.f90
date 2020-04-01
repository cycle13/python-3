MODULE common_da_tools_2d
!=======================================================================
!
! [PURPOSE:] Data assimilation tools for 1D models
!
!=======================================================================
!$USE OMP_LIB
  USE common_tools
  USE common_letkf

  IMPLICIT NONE

  PUBLIC

CONTAINS

!=======================================================================
!  LETKF DA for the 1D model
!=======================================================================
SUBROUTINE da_letkf(nx,nt,no,nens,nvar,xloc,tloc,xfens,xaens,obs,obsloc,ofens,Rdiag,loc_scale,&
           &inf_coefs,update_smooth_coef,dx,dy)
IMPLICIT NONE
INTEGER,INTENT(IN)         :: nx,nt,nvar,no              !State dimensions, space, time and variables
INTEGER,INTENT(IN)         :: nens                       !Number of ensemble members
REAL(r_size),INTENT(IN)    :: xloc(nx,2)                 !Location of state grid points (space)
REAL(r_size),INTENT(IN)    :: tloc(nt)                   !Location of state grid points (time)
REAL(r_size),INTENT(IN)    :: obsloc(no,3)               !Location of obs (space , time)
REAL(r_size),INTENT(IN)    :: xfens(nx,nens,nvar,nt)     !Forecast state ensemble  
REAL(r_size),INTENT(OUT)   :: xaens(nx,nens,nvar,nt)     !Analysis state ensemble 
REAL(r_size),INTENT(IN)    :: ofens(no,nens)             !Ensemble in observation space
REAL(r_size),INTENT(IN)    :: obs(no)                    !Observations 
REAL(r_size),INTENT(IN)    :: Rdiag(no)                  !Diagonal of observation error covariance matrix.
REAL(r_size),INTENT(IN)    :: inf_coefs(5)               !Mult inf, RTPS , RTPP , EPES, Additive inflation (State variables)
REAL(r_size),INTENT(IN)    :: update_smooth_coef         !Update smooth parameter.
REAL(r_size),INTENT(IN)    :: loc_scale(2)
REAL(r_size),INTENT(IN)    :: dx , dy 

REAL(r_size)               :: xfpert(nx,nens,nvar,nt)       !State and parameter forecast perturbations
REAL(r_size)               :: xapert(nx,nens,nvar,nt)       !State and parameter analysis perturbations
REAL(r_size)               :: xfmean(nx,nvar,nt)            !State and parameter ensemble mean (forecast)
REAL(r_size)               :: ofmean(no) , ofpert(no,nens)                         !Ensemble mean in obs space, ensemble perturbations in obs space (forecast)
REAL(r_size)               :: d(no)                                                !Observation departure


REAL(r_size)               :: ofpert_loc(no,nens) , Rdiag_loc(no)                  !Localized ensemble in obs space and diag of R
REAL(r_size)               :: Rwf_loc(no)                                          !Localization weigths.
REAL(r_size)               :: d_loc(no)                                            !Localized observation departure.
INTEGER                    :: no_loc                                               !Number of observations in the local domain.

REAL(r_size)               :: wa(nens,nens)                                        !Analysis weights
REAL(r_size)               :: wainf(nens,nens)                                     !Analysis weights after inflation.
REAL(r_size)               :: wamean(nens)                                         !Mean analysis weights
REAL(r_size)               :: pa(nens,nens)                                        !Analysis cov matrix in ensemble space)

REAL(r_size),PARAMETER     :: min_infl=1.0d0                                       !Minumn allowed multiplicative inflation.
REAL(r_size)               :: mult_inf
REAL(r_size)               :: grid_loc(3)
REAL(r_size)               :: work1d(nx)

INTEGER                    :: ix,ie,ke,it,iv,io


!Initialization
xfmean=0.0d0
xfpert=0.0d0
xaens =0.0d0
ofmean=0.0d0
ofpert=0.0d0
d     =0.0d0

wa=0.0d0
wamean=0.0d0
wainf =0.0d0
pa=0.0d0


!dx=xloc(2,1)-xloc(1,1)    !Assuming regular grid
!dy=xloc(nx+1,1)-xloc(1,1)

!Compute forecast ensemble mean and perturbations.


DO it = 1,nt

 DO ix = 1,nx

  DO iv = 1,nvar

   CALL com_mean( nens,xfens(ix,:,iv,it),xfmean(ix,iv,it) )
 
   xfpert(ix,:,iv,it) = xfens(ix,:,iv,it) - xfmean(ix,iv,it) 

  END DO

 END DO

END DO


!Compute mean departure and HXf


DO io = 1,no

    CALL com_mean( nens , ofens(io,:) , ofmean(io) )

    ofpert(io,:)=ofens(io,:) - ofmean(io) 

    d(io) = obs(io) - ofmean(io)   

ENDDO


!Main assimilation loop.

DO it = 1,nt

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ix, grid_loc, no_loc, ofpert_loc,d_loc,Rdiag_loc,Rwf_loc, & 
!$OMP &           wa,wamean,pa, wainf,mult_inf,     &
!$OMP &           iv, ie, ke) 
!$OMP DO 

  DO ix = 1,nx


   !Localize observations
   grid_loc(1)=xloc(ix,1)  !Set the grid point location in space
   grid_loc(2)=xloc(ix,2)
   grid_loc(3)=tloc(it)  !Set the grid point location in time


   CALL r_localization(no,no_loc,nens,ofpert,d,Rdiag,ofpert_loc,d_loc,Rdiag_loc,&
        & Rwf_loc,grid_loc,dx,obsloc,loc_scale)


   !WRITE(*,*)no_loc , no , nens , loc_scale , dx

   IF( no_loc > 0 )THEN   
    !We have observations for this grid point. Let's compute the analysis.

    !Set multiplicative inflation
    mult_inf=inf_coefs(1)  
  
    !Compute analysis weights
 
    CALL letkf_core( nens,no_loc,ofpert_loc(1:no_loc,:),Rdiag_loc(1:no_loc),   &
                    Rwf_loc(1:no_loc),d_loc(1:no_loc),mult_inf,wa,wamean,pa,min_infl )

    !Update state variables (apply RTPP and RTPS )
    IF( inf_coefs(2) /= 0.0d0) THEN                                                         !GYL - RTPP method (Zhang et al. 2005)
      CALL weight_RTPP(nens,inf_coefs(2),wa,wainf)                                          !GYL
    ELSE IF( inf_coefs(4) /= 0.0d0) THEN                                                    !EPES
      CALL weight_EPES(nens,inf_coefs(4),wa,wainf)
    ELSE IF( inf_coefs(5) /= 0.0d0) THEN                                                    !TODO: Code additive inflation
      write(*,*)"[Warning]: Additive inflation not implemented for LETKF yet"
    ELSE
      wainf = wa                                                                            !GYL
    END IF                                                                                  !GYL
   
    !Apply the weights and update the state variables. 
   
    DO iv=1,nvar
      !RTPS inflation is variable dependent, so we have to implement it here. 
      IF( inf_coefs(3) /= 0.0d0) THEN
        CALL weight_RTPS(nens,inf_coefs(3),wa,pa,xfens(ix,:,iv,it),wainf)
      ENDIF

      DO ie=1,nens
       xaens(ix,ie,iv,it) = xfmean(ix,iv,it)
       DO ke = 1,nens
          xaens(ix,ie,iv,it) = xaens(ix,ie,iv,it) &                                         !GYL - sum trans and transm here
              & + xfpert(ix,ke,iv,it) * (wainf(ke,ie) + wamean(ke))                         !GYL
       END DO
      END DO
    END DO

   ELSE  

    !We don't have observations for this grid point. We can do nothing :( 
    xaens(ix,:,:,it)=xfens(ix,:,:,it)

   ENDIF

  END DO

!$END OMP DO
!$OMP END PARALLEL


END DO  



IF( update_smooth_coef > 0.0d0 )THEN


! allocate( work1d( nx ) )

 DO it = 1,nt

  DO iv = 1,nvar

   DO ie = 1,nens

    !Smooth the update of each ensemble member.
    work1d = xaens(:,ie,iv,it) - xfens(:,ie,iv,it) !DA update
    CALL com_filter_lanczos( nx , update_smooth_coef , work1d )
    xaens(:,ie,iv,it) = xfens(:,ie,iv,it) + work1d 

   END DO

  END DO

 END DO

! deallocate( work1d )

ENDIF





END SUBROUTINE da_letkf


!=======================================================================
!  ETKF DA for 1D models
!=======================================================================

SUBROUTINE da_etkf(nx,nt,no,nens,nvar,xfens,xaens,obs,ofens,Rdiag,inf_coefs)

IMPLICIT NONE
INTEGER,INTENT(IN)         :: no , nt , nx , nvar        !Number of observations , number of times , number of state variables ,number of variables
INTEGER,INTENT(IN)         :: nens                       !Number of ensemble members
REAL(r_size),INTENT(IN)    :: xfens(nx,nens,nvar,nt)     !Forecast state ensemble  
REAL(r_size),INTENT(OUT)   :: xaens(nx,nens,nvar,nt)     !Analysis state ensemble 
REAL(r_size),INTENT(IN)    :: ofens(no,nens)             !Ensemble in observation space
REAL(r_size),INTENT(IN)    :: obs(no)                    !Observations 

REAL(r_size),INTENT(IN)    :: Rdiag(no)                  !Diagonal of observation error covariance matrix.

REAL(r_size),INTENT(INOUT) :: inf_coefs(5)               !Mult inf, RTPS , RTPP , EPES, Additive inflation (State variables)

REAL(r_size)               :: xfpert(nx,nens,nvar,nt)    !State and parameter forecast perturbations
REAL(r_size)               :: xfmean(nx,nvar,nt)         !State and parameter ensemble mean (forecast)
REAL(r_size)               :: ofmean(no) , ofpert(no,nens)                         !Ensemble mean in obs space, ensemble perturbations in obs space (forecast)
REAL(r_size)               :: d(no)                                                !Observation departure

REAL(r_size)               :: wa(nens,nens)                                        !Analysis weights
REAL(r_size)               :: wainf(nens,nens)                                     !Analysis weights after inflation.
REAL(r_size)               :: wamean(nens)                                         !Mean analysis weights
REAL(r_size)               :: pa(nens,nens)                                        !Analysis cov matrix in ensemble space)

REAL(r_size),PARAMETER     :: min_infl=1.0d0                                       !Minumn allowed multiplicative inflation.
REAL(r_size)               :: mult_inf
REAL(r_size)               :: work

REAL(r_size)               :: Rwf(no)                                              !Dummy variable.

INTEGER                    :: ie,ke,iv,io,ix,it

!Initialization
xfmean=0.0
xfpert=0.0

!Compute forecast ensemble mean and perturbations.

DO it = 1,nt

 DO ix = 1,nx

  DO iv = 1,nvar

   CALL com_mean( nens,xfens(ix,:,iv,it),xfmean(ix,iv,it) )
 
   xfpert(ix,:,iv,it) = xfens(ix,:,iv,it) - xfmean(ix,iv,it) 

  END DO

 END DO

END DO


!Compute mean departure and HXf

DO io = 1,no

    CALL com_mean( nens , ofens(io,:) , ofmean(io) )

    ofpert(io,:)=ofens(io,:) - ofmean(io) 

    d(io) = obs(io) - ofmean(io)   

ENDDO


!Set multiplicative inflation

mult_inf=inf_coefs(1)



!Main assimilation loop.


IF( no > 0 )THEN   
  !We have observations for this grid point. Let's compute the analysis.
   
  !Compute analysis weights
  Rwf=1.0d0 
  CALL letkf_core( nens,no,ofpert,Rdiag,Rwf,d,   &
                   mult_inf,wa,wamean,pa,min_infl )


  !Update state variables (apply RTPP and RTPS )
  IF( inf_coefs(2) /= 0.0d0) THEN                                                            !GYL - RTPP method (Zhang et al. 2005)
    CALL weight_RTPP(nens,inf_coefs(2),wa,wainf)                                             !GYL
  ELSE IF( inf_coefs(4) /= 0.0d0) THEN                                                       !EPES
    CALL weight_EPES(nens,inf_coefs(4),wa,wainf)
  ELSE IF( inf_coefs(5) /= 0.0d0) THEN                                                       !TODO: Code additive inflation
    write(*,*)"[Warning]: Additive inflation not implemented for ETKF yet"
  ELSE
    wainf = wa                                                                               !GYL
  END IF                                                                                     !GYL
   
  !Apply the weights and update the state variables.
  DO it=1,nt
   DO ix=1,nx 
    DO iv=1,nvar
     !RTPS inflation is variable dependent, so we have to implement it here. 
     IF( inf_coefs(3) /= 0.0d0) THEN
       CALL weight_RTPS(nens,inf_coefs(3),wa,pa,xfens(ix,:,iv,it),wainf)
     ENDIF

     DO ie=1,nens
       xaens(ix,ie,iv,it) = xfmean(ix,iv,it)
       DO ke = 1,nens
          xaens(ix,ie,iv,it) = xaens(ix,ie,iv,it) &                                          !GYL - sum trans and transm here
                        & + xfpert(ix,ke,iv,it) * ( wainf(ke,ie) + wamean(ke) )               !GYL
       END DO
     END DO

    END DO
   END DO
  END DO

ELSE  

  !We don't have observations for this grid point. We can do nothing :( 
  xaens(:,:,:,:)=xfens(:,:,:,:)

ENDIF



END SUBROUTINE da_etkf


!=======================================================================
!  R-LOCALIZATION for 2D models
!=======================================================================

SUBROUTINE  r_localization(no,no_loc,nens,ofpert,d,Rdiag,ofpert_loc,d_loc,Rdiag_loc    &
            &              ,Rwf_loc,grid_loc,dx,obsloc,loc_scale)
IMPLICIT NONE
INTEGER,INTENT(IN)         :: no , nens                  !Number of observations , number of ensemble members
REAL(r_size),INTENT(IN)    :: obsloc(no,3)               !Location of obs (space,time)
REAL(r_size),INTENT(IN)    :: grid_loc(3)                !Current grid point location (space,time)
!REAL(r_size),INTENT(IN)    :: x_min , x_max, y_min, y_max              !Space grid limits
REAL(r_size),INTENT(IN)    :: dx                         !Grid resolution
REAL(r_size),INTENT(IN)    :: ofpert(no,nens)            !Ensemble in observation space
REAL(r_size),INTENT(IN)    :: Rdiag(no)                  !Diagonal of observation error covariance matrix.
REAL(r_size)               :: d(no)                      !Observation departure
REAL(r_size)               :: loc_scale(2)               !Localization scale (space,time) , note that negative values
                                                         !will remove localization for that dimension.

INTEGER,INTENT(OUT)        :: no_loc                     !Number of observations in the local domain


REAL(r_size),INTENT(OUT)   :: ofpert_loc(no,nens) , Rdiag_loc(no)   !Localized ensemble in obs space and diag of R
REAL(r_size),INTENT(OUT)   :: Rwf_loc(no)                           !Localization weigthing function.
REAL(r_size),INTENT(OUT)   :: d_loc(no)                             !Localized observation departure.


INTEGER                    :: io
INTEGER                    :: if_loc(2)                        !1- means localize, 0- means no localization
REAL(r_size)               :: distance(2)                      !Distance between the current grid point and the observation (space,time)
REAL(r_size)               :: distance_tr(2)                   !Ignore observations farther than this threshold (space,time)


ofpert_loc=0.0d0
d_loc     =0.0d0
Rdiag_loc =0.0d0
Rwf_loc   =1.0d0
no_loc    =0



!If space localization will be applied.
if( loc_scale(1) <= 0)then
    !No localization in space.
    if_loc(1) = 0                                       !Ver bien estos parametros
    loc_scale(1)=1.0d0
else 
    if_loc(1) = 1
endif
!If temporal localization will be applied.
if( loc_scale(2) <= 0 )then
    !No localization in time.
    if_loc(2) = 0
    loc_scale(2)=1.0d0
else 
    if_loc(2) = 1
endif
!NOTE: To remove both space and time localizations is better to use the da_etkf routine.


!Compute distance threshold (space , time)
distance_tr = 2.0 * loc_scale * SQRT(10.0d0/3.0d0) 

!write(*,*) distance_tr , loc_scale

DO io = 1,no

   !Compute distance in space
!   CALL distance_x( grid_loc(1) , obsloc(io,1) , x_min , x_max , dx , distance(1) )
!   CALL distance_x( grid_loc(2) , obsloc(io,2) , y_min , y_max , dy , distance(2) )
   !write(*,*) grid_loc(1) , obsloc(io,1) , distance(1) 
   CALL com_distll_1(obsloc(io,2), obsloc(io,1) , grid_loc(2), grid_loc(1), distance(1))


   !write(*,*)obsloc(io,1), obsloc(io,2) , grid_loc(1), grid_loc(2), distance(1)
   !Compute distance in time 
   distance(2)=abs( grid_loc(3) - obsloc(io,3) )

   !write(*,*)distance , grid_loc , obsloc 


   !distance_tr check
   !write(*,*)distance(1),distance_tr(1),distance(2),distance_tr(2)
   IF( distance(1)*if_loc(1) < distance_tr(1) .and. distance(2)*if_loc(2) < distance_tr(2))THEN

     !This observation will be assimilated!
     no_loc = no_loc + 1

     ofpert_loc(no_loc,:) =ofpert(io,:)
     d_loc(no_loc)        =d(io)
     !Apply localization to the R matrix
     Rwf_loc(no_loc)      =  exp(-1.0* (( REAL( if_loc(1) , r_size ) * distance(1)/(2.0*loc_scale(1)))**2 + &
                                           ( REAL( if_loc(2) , r_size ) * distance(2)/(2.0*loc_scale(2)))**2))

     Rdiag_loc(no_loc)    =  Rdiag(io) / Rwf_loc( no_loc )

   ENDIF

ENDDO

END SUBROUTINE r_localization

!Compute distance between two grid points assuming cyclic boundary conditions.
SUBROUTINE distance_x( x1 , x2 , x_min , x_max , dx , distance)
IMPLICIT NONE
REAL(r_size), INTENT(IN)  :: x1,x2        !Possition
REAL(r_size), INTENT(IN)  :: dx           !Grid resolution
REAL(r_size), INTENT(IN)  :: x_min,x_max  !Grid limits
REAL(r_size), INTENT(OUT) :: distance     !Distance between x1 and x2


  distance=abs(x2 - x1) 

  if( distance > ( x_max - x1 ) + ( x2 - x_min ) + dx )then

    distance = ( x_max - x1 ) + ( x2 - x_min ) + dx

  endif

  if( distance > ( x_max - x2 ) + ( x1 - x_min ) + dx )then

    distance = ( x_max - x2 ) + ( x1 - x_min ) + dx 

  endif

END SUBROUTINE distance_x

!--------------------- --------------------------------------------------
! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
!-----------------------------------------------------------------------
SUBROUTINE com_distll_1(alon,alat,blon,blat,dist)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: alon
  REAL(r_size),INTENT(IN) :: alat
  REAL(r_size),INTENT(IN) :: blon
  REAL(r_size),INTENT(IN) :: blat
  REAL(r_size),INTENT(OUT) :: dist
  REAL(r_size),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2
  REAL(r_size) :: cosd

  lon1 = alon * pi * r180
  lon2 = blon * pi * r180
  lat1 = alat * pi * r180
  lat2 = blat * pi * r180

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)

  dist = ACOS( cosd ) * re

  RETURN
END SUBROUTINE com_distll_1



END MODULE common_da_tools_2d

