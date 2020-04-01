module covariance_matrix_tools
!=======================================================================
!
! [PURPOSE:] Common procedure for verification against grided data
!
!=======================================================================
!$USE OMP_LIB
  USE common_functions
! USE common_data
  USE common_random
!  USE common_verification
!  USE map_utils
!  USE common_smooth2d
!  USE ifport
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  !NAMELIST INPUT
!  INTEGER :: nbv=1  !Number of ensemble members.
!  integer, parameter :: max_npoints = 20 , max_vars = 20
!  integer :: npoints
!  real(r_size) :: plon(max_npoints)  , plat(max_npoints) , plev(max_npoints) 
!  character(50) :: pvarname(max_npoints)
!  real(r_size) :: dep(max_npoints) , error(max_npoints)
!  LOGICAL :: bootstrap=.true.
!  INTEGER :: bootstrap_samples=20
!  LOGICAL :: computehistogram=.true.
!  INTEGER :: max_histogram=50            !Number of bins

!  LOGICAL  :: computemoments=.true.
!  INTEGER  :: max_moments=4

!  LOGICAL      :: smoothcov=.false.
!  REAL(r_size) :: smoothcovlength=1.0d5  !Lanczos filter length scale in the same unit as dx.
!  REAL(r_size) :: smoothdx=1.0d3         !Data resolution 

!  LOGICAL :: computecovar=.true.
!  INTEGER :: skipx=1,skipy=1,skipz=1
!  INTEGER :: nignore=0
!  INTEGER , PARAMETER :: max_nignore = 20
!  character(50) :: ignorevarname(max_nignore)

!  LOGICAL :: computeindex=.true.       !Compute or not covariance and correlation index.
!  INTEGER :: delta = 1
!  INTEGER :: skip  = 1                 !Skip grid points in covariance strenght computation.

!  CHARACTER(clen) :: BASE_VARS(max_vars)  !This is a list of reference or observed variables.
!  CHARACTER(clen) :: COV_VARS(max_vars)   !This is the variables whose correlation / covariance with the observed variables will be measured.
!  INTEGER :: NBASE_VARS=0 , NCOV_VARS=0   !Number of base vars and cov vars.
!  REAL(r_sngl)    :: tr_rain = 5.0e-4 , tr_norain = 1.0e-5 !Condensate thresholds to separate rainy grid points.
!
!  character(50) :: inputendian='little_endian'
!  character(50) :: outputendian='big_endian'
!
!  CHARACTER(LEN=100) :: NAMELIST_FILE='./covariance_matrix.namelist'

!  !GRID INFO
!  TYPE(proj_info)  :: projection
!  TYPE(ctl_info)   :: ctl

!  integer :: nv !Total number of variables (vars times levs)

!  INTEGER,PARAMETER :: r_size=kind(0.0d0)
!  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)


CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

!SUBROUTINE read_namelist
!IMPLICIT NONE
!LOGICAL  :: file_exist
!INTEGER  :: ierr

!!In the current version PLEV represents not only the level but also the
!!variable. 

!NAMELIST / GENERAL / nbv , npoints , plon , plat , pvarname , plev,   &
!                     dep , error , bootstrap , bootstrap_samples  ,   &
!                     nignore , ignorevarname , skipx , skipy , skipz, &
!                     inputendian , outputendian , smoothcov ,         &
!                     smoothdx , smoothcovlength , computemoments ,   &
!                     max_moments , computeindex , delta , computecovar, &
!                     computehistogram , max_histogram , base_vars ,  &
!                     cov_vars , nbase_vars , ncov_vars , tr_rain , tr_norain

!plon=undef
!plat=undef
!plev=undef
!dep=undef
!error=undef

!INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)
!
!IF( .NOT. file_exist )THEN
!  WRITE(*,*)"ERROR: Could not find namelist"
!ENDIF
!
!OPEN(54,FILE=NAMELIST_FILE)
!READ(54,NML=GENERAL,IOSTAT=IERR)
!IF(IERR /=0)THEN
!WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
!WRITE(*,*)"Using default values"
!ENDIF
!REWIND(54)

!END SUBROUTINE read_namelist

SUBROUTINE compute_covar(nbv,nx,ny,nz,mv,ensemble,pensemble,bootstrap_samples,undef_mask,covar,covarmean,covarstd)
!Compute the covariance between the nx,ny,nz points of the array ensemble and
!the matrix M vectors pensemble.
IMPLICIT NONE
INTEGER     , INTENT(IN)  :: nbv,nx,ny,nz,mv
REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv)
LOGICAL, INTENT(IN)       :: undef_mask(nx,ny,nz)
REAL(r_sngl), INTENT(IN)  :: pensemble(nbv,mv)
INTEGER,      INTENT(IN)  :: bootstrap_samples
REAL(r_sngl), INTENT(OUT) :: covar(nx,ny,nz,mv) 
REAL(r_sngl), INTENT(OUT) :: covarmean(nx,ny,nz,mv)
REAL(r_sngl), INTENT(OUT) :: covarstd(nx,ny,nz,mv)
REAL(r_sngl)              :: cov , localcov(nx,ny,nz,mv)
INTEGER                   :: sampleindex(nbv)
INTEGER                   :: ii , jj , kk , im , is


!Ititialize output
covar    =0.0e0
covarmean=0.0e0
covarstd =0.0e0

DO im = 1,mv
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,cov)
  DO ii = 1,nx
   DO jj = 1,ny
    DO kk = 1,nz
     IF( undef_mask(ii,jj,kk) ) then
      CALL com_covar_sngl(nbv,ensemble(ii,jj,kk,:),pensemble(:,mv),covar(ii,jj,kk,im))
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDDO

!PERFORM BOOTSTRAP ESTIMATION 
IF( bootstrap_samples > 0 )THEN
     write(*,*)"Using bootstrap, with ",bootstrap_samples," samples"
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,im,is,localcov,sampleindex) 

    DO is = 1,bootstrap_samples
     CALL generate_sample(sampleindex,nbv)
     localcov=0.0d0
     DO im = 1,mv
      DO ii = 1,nx
       DO jj = 1,ny
        DO kk = 1,nz
         IF( undef_mask(ii,jj,kk) ) then
           CALL com_covar_sngl_sample(nbv,pensemble(:,mv),ensemble(ii,jj,kk,:),sampleindex,localcov(ii,jj,kk,im))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
     ENDDO
!$OMP CRITICAL     
    covarmean=covarmean+localcov
    covarstd =covarstd +localcov**2
!$OMP END CRITICAL
   ENDDO
!$OMP END PARALLEL DO
   covarmean=covarmean/REAL(bootstrap_samples,r_sngl)
   covarstd =sqrt( covarstd/REAL(bootstrap_samples,r_sngl)-( covarmean )**2 )

ENDIF

END SUBROUTINE compute_covar

SUBROUTINE covariance_strenght(var1_ens,var2_ens,mask,nmask,covst,corrst,corrstdist,cov_profile,corr_profile, & 
                               num_profile,nx,ny,nz,nbv,delta,skip_cov,undef_mask,undefbin)

IMPLICIT NONE
INTEGER, INTENT(IN) :: nx,ny,nz,nbv,delta,skip_cov
REAL(r_sngl), INTENT(IN) :: var1_ens(nx,ny,nz,nbv)  , undefbin  !Base variable
REAL(r_sngl), INTENT(IN) :: var2_ens(nx,ny,nz,nbv)              !Covariance variable
INTEGER,      INTENT(IN) :: mask(nx,ny,nz)                      !mask that allows a segemntation of the domain.
INTEGER,      INTENT(IN) :: nmask                               !maximum number of categories in mask.
LOGICAL ,     INTENT(IN) :: undef_mask(nx,ny,nz) 
REAL(r_sngl), INTENT(OUT):: covst(nx,ny,nz),corrst(nx,ny,nz),corrstdist(nx,ny,nz) !Covariance and correlation strhenght index.
REAL(r_sngl)             :: var1_std(nx,ny,nz) , var2_std(nx,ny,nz)
INTEGER :: imin , imax , jmin , jmax , ii , jj , kk , iii , jjj , contador ,imask
REAL(r_sngl)             :: pensemble(nbv) , tmpcov , tmpcorr , tmpdist
REAL(r_sngl)             :: tmp_cov_profile(delta) , tmp_corr_profile(delta) !, dist_profile(delta)
INTEGER                  :: tmp_num_profile(delta) , current_index
REAL(r_sngl),INTENT(OUT) :: cov_profile(delta,nz,nmask) , corr_profile(delta,nz,nmask)
INTEGER     ,INTENT(OUT) :: num_profile(delta,nz,nmask)

covst=0
corrst=0
corrstdist=0

cov_profile=0.0e0     !Note that if OPENMP is used this variable is initialized to 0
corr_profile=0.0e0    !Note that if OPENMP is used this variable is initialized to 0
num_profile=0.0e0     !Note that if OPENMP is used this variable is initialized to 0

!tr_norain=1.0e-5
!tr_rain  =1.0e-3

!$OMP PARALLEL DO PRIVATE(ii,jj,kk)
DO ii=1,nx
 DO jj=1,ny
  DO kk=1,nz
    CALL com_stdev_sngl(nbv,var1_ens(ii,jj,kk,:),var1_std(ii,jj,kk))
    CALL com_stdev_sngl(nbv,var2_ens(ii,jj,kk,:),var2_std(ii,jj,kk))
  ENDDO
 ENDDO
ENDDO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(imin,imax,jmin,jmax,ii,jj,kk,iii,jjj,pensemble,tmpcov,tmpcorr,contador,tmpdist, &
!$OMP tmp_cov_profile,tmp_corr_profile,tmp_num_profile,current_index,imask)                               &
!$OMP reduction(+:cov_profile,corr_profile,num_profile)

DO ii=1,nx
 DO jj=1,ny
  DO kk=1,nz
   IF( undef_mask(ii,jj,kk) )THEN

    imin=ii-delta
    imax=ii+delta
    if( imin < 1 )imin=1
    if( imax > nx)imax=nx
    jmin=jj-delta
    jmax=jj+delta
    if( jmin < 1 )jmin=1
    if( jmax > ny)jmax=ny

    pensemble=var1_ens(ii,jj,kk,:)

    tmp_cov_profile=0.0e0
    tmp_corr_profile=0.0e0
    !dist_profile=0.0e0
    tmp_num_profile=0

    DO iii=imin,imax,skip_cov
     DO jjj=jmin,jmax,skip_cov

      tmpdist= sqrt(real(iii-ii,r_sngl)**2 + real(jjj-jj,r_sngl)**2) !/real(delta,r_sngl)

      IF( undef_mask(iii,jjj,kk) .and. tmpdist <= real(delta,r_sngl) )THEN

       IF( var1_std(ii,jj,kk) > 0 .and. var2_std(iii,jjj,kk) > 0 )THEN

        CALL com_covar_sngl(nbv,var2_ens(iii,jjj,kk,:),pensemble,tmpcov)
  
        tmpcorr=tmpcov/( var1_std(ii,jj,kk) * var2_std(iii,jjj,kk) ) 

        current_index=FLOOR( tmpdist ) + 1

        IF( current_index >= delta) current_index=delta

        tmp_cov_profile(current_index)=tmp_cov_profile(current_index)+ABS(tmpcov) 
  
        tmp_corr_profile(current_index)=tmp_corr_profile(current_index)+ABS(tmpcorr)

        !dist_profile(current_index)=dist_profile(current_index)+tmpdist*ABS(tmpcorr)

        tmp_num_profile(current_index)=tmp_num_profile(current_index) + 1
    
       ENDIF
  
      ENDIF
       
     ENDDO
    ENDDO
   
    DO imask=1,nmask

      IF( mask(ii,jj,kk) == imask )THEN
        !Accumulate the profiles for the corresponding category.
        cov_profile(:,kk,imask) =cov_profile(:,kk,imask) + tmp_cov_profile 
        corr_profile(:,kk,imask)=corr_profile(:,kk,imask) + tmp_corr_profile
        num_profile(:,kk,imask)=num_profile(:,kk,imask) + tmp_num_profile
      ENDIF

    ENDDO
  
    contador=0 
    DO iii=1,delta
        IF(  tmp_num_profile(iii) >= 1 )THEN
          tmp_cov_profile(iii) = tmp_cov_profile(iii) / tmp_num_profile(iii)
          tmp_corr_profile(iii) = tmp_corr_profile(iii) / tmp_num_profile(iii)

          covst(ii,jj,kk) = covst(ii,jj,kk) + tmp_cov_profile(iii) 
          corrst(ii,jj,kk) = corrst(ii,jj,kk) + tmp_corr_profile(iii)
          corrstdist(ii,jj,kk) = corrstdist(ii,jj,kk) + tmp_corr_profile(iii) * iii
 
          contador=contador + 1
        ENDIF
    ENDDO

    IF( contador >= 0 )THEN
      corrstdist(ii,jj,kk)=corrstdist(ii,jj,kk)/corrst(ii,jj,kk)
      covst(ii,jj,kk)=covst(ii,jj,kk)/REAL(contador,r_sngl)
      corrst(ii,jj,kk)=corrst(ii,jj,kk)/REAL(contador,r_sngl)
    ELSE
      covst(ii,jj,kk)=undefbin
      corrst(ii,jj,kk)=undefbin
      corrstdist(ii,jj,kk)=undefbin
    ENDIF
   
   ELSE
 
    covst(ii,jj,kk)=undefbin
    corrst(ii,jj,kk)=undefbin
    corrstdist(ii,jj,kk)=undefbin

   ENDIF

  ENDDO
 ENDDO
ENDDO

!$OMP END PARALLEL DO

  DO ii=1,delta
   DO kk=1,nz
    DO jj=1,2
     IF( num_profile(ii,kk,jj) >= 1 )THEN
       cov_profile(ii,kk,jj)=cov_profile(ii,kk,jj) / real( num_profile(ii,kk,jj) , r_sngl )
       corr_profile(ii,kk,jj)=corr_profile(ii,kk,jj) / real( num_profile(ii,kk,jj) , r_sngl )
     ELSE
       cov_profile(ii,kk,jj)=undefbin
       corr_profile(ii,kk,jj)=undefbin
     ENDIF
    ENDDO
   ENDDO
  ENDDO

END SUBROUTINE covariance_strenght


SUBROUTINE generate_sample(sampledindex,n)
!Resample randomly picking with substitution.
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(OUT) :: sampledindex(n)
REAL(r_size) :: randomloc(n)
INTEGER :: intind(n) , ii

 CALL com_rand(n,randomloc)
 sampledindex=INT(1+randomloc*(n-1))

END SUBROUTINE generate_sample


SUBROUTINE read_ensemble(path,ensemble,undef_mask,nx,ny,nbv,selected_fields,n_selected_fields,undef,ie)
IMPLICIT NONE
INTEGER , INTENT(IN)       :: selected_fields(n_selected_fields) 
INTEGER , INTENT(IN)       :: n_selected_fields
INTEGER , INTENT(IN)       :: nx , ny , nz , nbv
REAL(r_sngl) , INTENT(IN)  :: undef
REAL(r_sngl) , INTENT(OUT) :: ensemble(nx,ny,n_selected_fields,nbv)
LOGICAL      , INTENT(OUT) :: undef_mask(nx,ny,n_selected_fields)
REAL(r_sngl)               :: bufr(nx,ny)
CHARACTER(*) , INTENT(IN)  :: path , ie                            !Data path and input endian
INTEGER                    :: ifield , ibv , iunit , reclength , ii , jj
CHARACTER(40)              :: filename

INQUIRE(IOLENGTH=reclength)reclength
reclength=nx*ny*reclength

undef_mask=.true.


filename='____.grd'

DO ibv = 1,nbv
  
   write(filename,'(I4.4)')nbv 

   OPEN(iunit, FILE=path // filename ,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=ie)

   !Read only the selected fields.
   DO ifield = 1,n_selected_fields
      READ(iunit,rec=selected_fields(ifield)) bufr 
      ensemble(:,:,ifield,ibv)=bufr 
      DO ii=1,nx
       DO jj=1,ny
         if( ensemble(ii,jj,ifield,ibv)==undef)then
           undef_mask=.false.
         endif
       ENDDO
      ENDDO 
   ENDDO

ENDDO 

END SUBROUTINE read_ensemble

END MODULE covariance_matrix_tools



