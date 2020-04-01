MODULE motion_vectors
!=======================================================================
!
! [PURPOSE:] Data assimilation tools for 1D models
!
!=======================================================================
!$USE OMP_LIB

  IMPLICIT NONE

  PUBLIC

CONTAINS


! THESE ROUTINES ARE FOR A MOTION VECTOR COMPUTATION USING PYTHON-FORTRAN 
! COMPILE USING: f2py -c -m test /home/jruiz/PYTHON/motion_vectors.f90
! FOR PARALLEL FUNCTIONALITY f2py -c -lgomp --f90flags="-fopenmp -lgomp" -m motion_vectors motion_vectors_parallel.f90 
! FOR DEBUG  FUNCTIONALITY f2py -c -lgomp --f90flags="-g -traceback" -m motion_vectors motion_vectors_parallel.f90 

SUBROUTINE MOTION_VECTOR(field_t0,field_t1,dt,dx,box_size,sigma,nx,ny,desp_max,min_box_fraction,u_motion,v_motion,max_corr &
       &      ,trend,nref,aux_output,aux_inputi,aux_inputj,motion_vector_option,motion_vector_weigth,motion_vector_norm)

IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , kkk , ini_i , ini_j , end_i , end_j , contador , n_points , n_points2
REAL(dp), PARAMETER :: UNDEF=-999.0d0 , min_dbz=-999.0d0
INTEGER , PARAMETER :: METRIC_OPTION = 1 ! Metric option for similarity quantification (1 = correlation , ...)
REAL(dp), INTENT(IN) :: field_t0(nx,ny),field_t1(nx,ny),dt,dx,sigma,min_box_fraction
INTEGER , INTENT(IN) :: nx,ny,desp_max,box_size
INTEGER , INTENT(IN) :: aux_inputi , aux_inputj
REAL(dp), INTENT(OUT) :: u_motion(nx,ny) , v_motion(nx,ny) , max_corr(nx,ny) , trend(nx,ny)
INTEGER , INTENT(OUT) :: nref(nx,ny) !Number of points with ref > min_dbz within the box.
REAL(dp), INTENT(OUT) :: aux_output(2*desp_max+ 1,2*desp_max+1)
INTEGER  :: tmp_dim(1)
REAL(dp),ALLOCATABLE :: corr_matrix(:,:) , vector_1(:) , vector_2(:) 
REAL(dp),ALLOCATABLE :: correlation_weigths(:,:) , vector_cw(:)   !Gaussian weights for weighted correlation.
REAL(dp) :: box_fraction , dist , tmp_corr  ,  despi ,despj,weigth
INTEGER  :: maxi , maxj  , i_point
REAL(dp) , PARAMETER :: correlation_threshold = 0.3  ! $
INTEGER , INTENT(IN) :: MOTION_VECTOR_OPTION  !0 - Maximum corr , 1 Average sourrounding max corr, 2 - Average over entire displacements.
INTEGER , INTENT(IN) :: MOTION_VECTOR_WEIGTH  !0 - Uniform weigth , 1 - Gaussian weigths.
INTEGER , INTENT(IN) :: MOTION_VECTOR_NORM    !0 - Linear correlation, 1 CDF Correlation , 2 -MSE
INTEGER , PARAMETER :: NT=6
REAL(dp)  :: THRESHOLDS(NT)
REAL(dp), ALLOCATABLE  :: field2_t0(:,:) , field2_t1(:,:)
INTEGER  :: nx2 , ny2

! Zero padding 
nx2= nx + 2*( desp_max + box_size )
ny2= ny + 2*( desp_max + box_size )

ALLOCATE( field2_t0(nx2,ny2) , field2_t1(nx2,ny2) )

field2_t0=min_dbz
field2_t1=min_dbz

!WRITE(*,*)field_t0(1,1)
!WRITE(*,*)field_t0(nx,1)
!WRITE(*,*)field_t0(1,ny)
!WRITE(*,*)field_t0(nx,ny)

field2_t0(desp_max+box_size+1:desp_max+box_size+nx,desp_max+box_size+1:desp_max+box_size+ny)=field_t0
field2_t1(desp_max+box_size+1:desp_max+box_size+nx,desp_max+box_size+1:desp_max+box_size+ny)=field_t1

where( field2_t0 < min_dbz) 
   field2_t0=min_dbz
end where
where( field2_t1 < min_dbz)
   field2_t1=min_dbz
end where 


DO ii=1,NT
   THRESHOLDS(ii)=ii*10
ENDDO 

!WRITE(*,*) " INSIDE W MOTION VECTOR ROUTINE "
!WRITE(*,*) " FIRST REF FIELD AT T 0 = ",field_t0(1,1)
!WRITE(*,*) " FIRST REF FIELD AT T 1 = ",field_t1(1,1)
!WRITE(*,*) " DELTA T                  ",dt
!WRITE(*,*) " DELTA X                  ",dx
!WRITE(*,*) " SIGMA                    ",box_size
!WRITE(*,*) " SIGMA_THRESHOLD          ",sigma
!WRITE(*,*) " NX ,  NY                 ",nx,ny
!WRITE(*,*) " DESP MAX                 ",desp_max
!WRITE(*,*) " MIN_BOX_FRACTION         ",min_box_fraction
!WRITE(*,*) " THRESHOLDS               ",THRESHOLDS

! First find the box size based on sigma and the threshold. In this implementation
! the box size is computed from sigma and sigma_threshold.

ALLOCATE(corr_matrix(desp_max*2+1,desp_max*2+1) , vector_1((box_size*2+1)**2) , vector_2((box_size*2+1)**2)  )

! Compute weigths for correlation.
ALLOCATE( correlation_weigths(box_size*2+1,box_size*2+1) , vector_cw((box_size*2+1)**2) )


SELECT CASE ( MOTION_VECTOR_WEIGTH )

CASE( 0 ) !Uniform weigths

   correlation_weigths=1.0d0

CASE( 1 ) !Gaussian weigths

  DO ii=-box_size,box_size
    DO jj=-box_size,box_size
      dist=( (ii*dx) ** 2 + (jj*dx) ** 2 )
      correlation_weigths( ii+box_size+1 , jj+box_size+1 ) = exp (-( dist / (sigma**2) ) )
    ENDDO
  ENDDO

END SELECT

!Initialize u_motion and v_motion with UNDEF values.
u_motion=undef
v_motion=undef
max_corr=-1.0d0
nref=0

n_points=(box_size * 2 + 1)**2
tmp_dim(1)=n_points

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,ini_i,end_i,ini_j,end_j,iii,jjj,kkk,corr_matrix & 
!$OMP&           ,vector_1,vector_2,vector_cw,contador,box_fraction      &
!$OMP&           ,despi,despj,tmp_corr,weigth,maxi,maxj,n_points2) 
!$OMP DO 


DO ii = desp_max+box_size+1,desp_max+box_size+nx

!  WRITE(*,*)"Processing row number ",ii
  DO jj = desp_max+box_size+1,desp_max+box_size+ny
 
        corr_matrix=-1.0d0

        ini_i=ii-box_size
        end_i=ii+box_size
        ini_j=jj-box_size
        end_j=jj+box_size

        vector_1=RESHAPE(field2_t0(ini_i:end_i,ini_j:end_j),tmp_dim)
        vector_cw=RESHAPE(correlation_weigths,tmp_dim)
        vector_cw=vector_cw/SUM(vector_cw)

        nref(ii-box_size-desp_max,jj-box_size-desp_max)=0.0d0
        DO iii=1,n_points
           if( vector_1(iii) > min_dbz)then
             nref(ii-box_size-desp_max,jj-box_size-desp_max)=nref(ii-box_size-desp_max,jj-box_size-desp_max)+1
           endif
        ENDDO

        if( nref(ii-box_size-desp_max,jj-box_size-desp_max) < min_box_fraction )CYCLE

        maxi=0
        maxj=0


        DO iii= -desp_max,desp_max
          DO jjj= -desp_max,desp_max

            vector_2=RESHAPE(field2_t1(ini_i+iii:end_i+iii,ini_j+jjj:end_j+jjj),tmp_dim)

            SELECT CASE( MOTION_VECTOR_NORM )
 
            CASE( 0 ) !Weigthed linear correlation. 
             CALL CORRELATION( vector_1,vector_2,vector_cw,n_points,tmp_corr)
            CASE( 1 ) !Weigthed CDF correlation
             CALL CDFCORRELATION(vector_1,vector_2,vector_cw,n_points,THRESHOLDS,NT,tmp_corr)
            CASE( 2 ) !Weighted MSE
             CALL MSE(vector_1,vector_2,vector_cw,n_points,tmp_corr)
             IF( tmp_corr .NE. UNDEF  )THEN
                tmp_corr = 1.0d0/(tmp_corr+1.0d0) !The index has to be possitively oriented.
             ENDIF
            END SELECT

            corr_matrix(iii+desp_max+1,jjj+desp_max+1)=tmp_corr
            IF( tmp_corr >= max_corr(ii-box_size-desp_max,jj-box_size-desp_max) )THEN
                max_corr(ii-box_size-desp_max,jj-box_size-desp_max)=tmp_corr
                maxj=jjj
                maxi=iii
                !Compute the mean reflectivity trend in the box.
                trend(ii-box_size-desp_max,jj-box_size-desp_max)=SUM( vector_2-vector_1, &
      &         vector_1>min_dbz .OR. vector_2 > min_dbz)/SUM(vector_1-vector_1+1.0d0,vector_1>min_dbz &
      &         .OR. vector_2 > min_dbz)
            ENDIF

          ENDDO
        ENDDO
 
     !!!!! FOR DEBUG
     IF( ii==aux_inputi+desp_max+box_size .AND. jj==aux_inputj+desp_max+box_size )THEN
         aux_output=corr_matrix
         
         !WRITE(*,*)maxj,maxi,max_corr(ii-box_size-desp_max,jj-box_size-desp_max)
         !vector_2=RESHAPE(field2_t1(ini_i+0:end_i+0,ini_j+1:end_j+1),tmp_dim)
         !DO iii=1,tmp_dim(1)
         !   WRITE(*,*)vector_2(iii),vector_1(iii)
         !ENDDO
       
     ENDIF 
     !!!!!!! $

     SELECT CASE ( MOTION_VECTOR_OPTION )

     CASE( 0 ) !Maximum correlation criteria. 

       despi=REAL(maxi,dp)
       despj=REAL(maxj,dp)
 
       IF( max_corr(ii-box_size-desp_max,jj-box_size-desp_max) > 0.0d0 )THEN
       u_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despj*dx/dt
       v_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despi*dx/dt
       ENDIF

     CASE( 1 ) !Mean correlation around the maximum value
         !Find the weighted location for the minimum.
         despi=0.0d0
         despj=0.0d0
         weigth=0.0d0
         DO iii=maxi-1,maxi+1
           DO jjj=maxj-1,maxj+1 
             IF( iii <= desp_max  .AND. iii >= -desp_max .AND. jjj <= desp_max .AND. jjj >= -desp_max )THEN
                IF( corr_matrix( iii+desp_max+1 , jjj+desp_max + 1) > 0)THEN
                    tmp_corr = corr_matrix( iii+desp_max+1 , jjj+desp_max + 1)**2
                    weigth= weigth + tmp_corr
                    despi = despi  + tmp_corr *iii
                    despj = despj  + tmp_corr *jjj
                ENDIF
             ENDIF                
           ENDDO
         ENDDO

         IF( weigth > 0.0)THEN
           u_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despj*dx/(dt*weigth)
           v_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despi*dx/(dt*weigth)
         ENDIF

     CASE( 2 ) !Mean correlation using the highest correlation points.
        weigth=0.0d0
        despi=0.0d0
        despj=0.0d0
        DO iii= -desp_max,desp_max 
           DO jjj= -desp_max,desp_max
             tmp_corr = corr_matrix( iii + desp_max + 1, jjj + desp_max + 1 )
             IF( max_corr(ii-box_size-desp_max,jj-box_size-desp_max) - tmp_corr < correlation_threshold .AND. tmp_corr > 0.0 )THEN
                 tmp_corr = tmp_corr ** 2
                 weigth=weigth+tmp_corr
                 despi=despi + iii * (tmp_corr)
                 despj=despj + jjj * (tmp_corr)
             ENDIF
           ENDDO
        ENDDO

         !Get the location of the maximum correlation within
         !corr_matrix
         IF( weigth > 0.0)THEN
           u_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despj*dx/(dt*weigth)
           v_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despi*dx/(dt*weigth)
         ENDIF

     END SELECT

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL

DEALLOCATE(corr_matrix , vector_1 , vector_2  )
DEALLOCATE( correlation_weigths , vector_cw )
DEALLOCATE( field2_t0 , field2_t1 )

!WRITE(*,*)"FINISH MOTION VECTOR COMPUTATION"
END SUBROUTINE MOTION_VECTOR


SUBROUTINE HORN_SCHUNCK(field_t0,field_t1,dt,dx,nx,ny,u_motion,v_motion,alpha,max_iter)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , ini_i , ini_j , end_i , end_j , contador , n_points
REAL(dp), PARAMETER :: UNDEF=-999.0d0 , min_box_fraction = 0.01d0 , min_dbz=0.0d0
INTEGER , PARAMETER :: METRIC_OPTION = 1 ! Metric option for similarity quantification (1 = correlation , ...)
REAL(dp), INTENT(IN) :: field_t0(nx,ny),field_t1(nx,ny),dt,dx
INTEGER , INTENT(IN) :: nx,ny
REAL(dp), INTENT(OUT) :: u_motion(nx,ny) , v_motion(nx,ny)
REAL(dp)              :: u_tmp(nx,ny)    , v_tmp(nx,ny)
REAL(dp)              :: Ix(nx,ny) , Iy(nx,ny) , It(nx,ny) , Ixsq(nx,ny) , Iysq(nx,ny)
INTEGER :: ITER
INTEGER , INTENT(IN) :: MAX_ITER
REAL(dp), INTENT(IN) :: ALPHA
REAL(dp), PARAMETER :: UPDATE_THRESHOLD = 1d-3
REAL(dp)            :: ALPHASQ
LOGICAL             :: CONVERGE
REAL(dp)            :: um,vm , tmp , UPDATE
REAL(dp)            :: mean_u , mean_v

WRITE(*,*) " INSIDE HORN SCHUNCK "
WRITE(*,*) " FIRST REF FIELD AT T 0 = ",field_t0(1,1)
WRITE(*,*) " FIRST REF FIELD AT T 1 = ",field_t1(1,1)
WRITE(*,*) " DELTA T                  ",dt
WRITE(*,*) " DELTA X                  ",dx
WRITE(*,*) " NX ,  NY                 ",nx,ny
! Initialize u_motion and v_motion with UNDEF values.
u_motion=0.0d0
v_motion=0.0d0
u_tmp=u_motion
v_tmp=v_motion
! First compute local trend and spatial derivative.
It = field_t1 - field_t0
Ix(2:nx-1,:)=field_t0(3:nx,:) - field_t0(1:nx-2,:)
Iy(:,2:ny-1)=field_t0(:,3:ny) - field_t0(:,1:ny-2)
Iysq=Iy**2
Ixsq=Ix**2
ALPHASQ=ALPHA**2
ITER=0
CONVERGE=.FALSE.

DO WHILE ( .NOT.CONVERGE .AND. ITER < MAX_ITER )

!Cost function minimization.
    UPDATE=0.0d0
    mean_u=0.0d0
    mean_v=0.0d0
    DO ii = 2 , nx - 1
     DO jj = 2 , ny - 1

       um=0.25*(u_motion(ii+1,jj)+u_motion(ii-1,jj)+u_motion(ii,jj+1)+u_motion(ii,jj-1))
       vm=0.25*(v_motion(ii+1,jj)+v_motion(ii-1,jj)+v_motion(ii,jj+1)+v_motion(ii,jj-1))

       tmp= -( Ix(ii,jj)*um + Iy(ii,jj)*vm + It(ii,jj) )/(ALPHASQ + Ixsq(ii,jj) + Iysq(ii,jj))

       u_tmp(ii,jj)= um + Ix(ii,jj)*tmp
       v_tmp(ii,jj)= vm + Iy(ii,jj)*tmp

       UPDATE=UPDATE + ABS(u_motion(ii,jj)-u_tmp(ii,jj)) + ABS(v_motion(ii,jj)-v_tmp(ii,jj))

       mean_u=mean_u+u_motion(ii,jj)
       mean_v=mean_v+v_motion(ii,jj)
     ENDDO
    ENDDO
    mean_u=mean_u/REAL(nx*ny,dp)
    mean_v=mean_v/REAL(nx*ny,dp)

    UPDATE=UPDATE/REAL(nx*ny,dp)
    IF( UPDATE <= UPDATE_THRESHOLD)THEN
      CONVERGE=.TRUE. !We reach convergence! We will stop iteration process.
    ENDIF

    !Boundary conditions (zero derivative at the boundaries.
    u_motion=u_tmp
    v_motion=v_tmp


    !u_motion(1,:)=u_motion(2,:)
    !u_motion(nx,:)=u_motion(nx-1,:)
    !u_motion(:,1)=u_motion(:,2)
    !u_motion(:,ny)=u_motion(:,ny-1)

    !v_motion(1,:)=v_motion(2,:)
    !v_motion(nx,:)=v_motion(nx-1,:)
    !v_motion(:,1)=v_motion(:,2)
    !v_motion(:,ny)=v_motion(:,ny-1)

    !Boundary conditions equal to mean u and v inside the domain.
    u_motion(1:2,:)=mean_u
    u_motion(nx-1:nx,:)=mean_u
    u_motion(:,1:2)=mean_u
    u_motion(:,ny-1:ny)=mean_u

    v_motion(1:2,:)=mean_v
    v_motion(nx-1:nx,:)=mean_v
    v_motion(:,1:2)=mean_v
    v_motion(:,ny-1:ny)=mean_v


ITER=ITER+1
WRITE(*,*)"ITERATION ",ITER," MOTION VECTOR UPDATE ",UPDATE
ENDDO


u_motion=u_motion*dx/dt
v_motion=v_motion*dx/dt


WRITE(*,*)"FINISH MOTION VECTOR COMPUTATION"
END SUBROUTINE HORN_SCHUNCK


SUBROUTINE CDFCORRELATION( A,B,W,L,THRESHOLDS,NT,CDFINDEX)
!This subroutine computes a correspondance between two vectors
!based on the comparison of the correspondent CDF of the two vectors.
!The idea is taking the probability of detection (POD) for different 
!thresholds and output a weigthted sum of the PODs corresponding to 
!different thresholds.
! CDFINDEX= (N1xPOD1 + ... NntxPODnt)/( N1 + ..... Nnt)
!Where nt is the number of thresholds (input), N is the number of points above
!the threshold in A and POD1 is the number of points which is above the threshold in A and B simultaneously.

INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER , INTENT(IN) :: L , NT  !Number of elements in inputs vectors and number of thresholds.
REAL(dp), INTENT(IN) :: A(L) , B(L) , THRESHOLDS(NT)  !Input vectors A is the reference vector.
REAL(dp), INTENT(IN) :: W(L) !Weigths.
REAL(dp), PARAMETER  :: UNDEF=-999.0d0
REAL(dp)             :: C(L) !Auxiliar arrays
REAL(dp), INTENT(OUT) :: CDFINDEX   !Output
INTEGER   :: ii,jj,N

    CDFINDEX=0.0d0
    N=0

DO ii=1,NT
       C=0.0d0
    WHERE( A >= THRESHOLDS(ii) .AND. B >= THRESHOLDS(ii) )
       C=1.0d0 * W
    ENDWHERE
    CDFINDEX= CDFINDEX + SUM(  C  )
       C=0.0d0
    WHERE( A >= THRESHOLDS(ii) )
       C=W
    ENDWHERE
    N=N + SUM( C )
ENDDO

!  IF( N .GT. 0)THEN
!  CDFINDEX=CDFINDEX/REAL(N,dp)
!  ELSE
!  CDFINDEX=0.0d0
!  ENDIF

END SUBROUTINE CDFCORRELATION

SUBROUTINE MSE( A,B,W,L,WMSE)
!This subroutine computes a weigthed MSE between two vectors.

INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER , INTENT(IN) :: L  !Number of elements in inputs vectors 
REAL(dp), INTENT(IN) :: A(L) , B(L)  !Input vectors A is the reference vector.
REAL(dp), INTENT(IN) :: W(L) !Weigths.
REAL(dp), PARAMETER  :: UNDEF=-999.0d0
REAL(dp)             :: TW       !Auxiliar arrays
REAL(dp), INTENT(OUT) :: WMSE   !Output
INTEGER   :: ii,jj,N

    WMSE=0.0d0
    TW=0.0d0

DO ii=1,L
    IF( A(ii) .NE. UNDEF .AND. B(ii) .NE. UNDEF )THEN
       WMSE = WMSE + W(ii)*( ( A(ii) - B(ii) )**2 )
       TW = TW + W(ii)
    ENDIF
ENDDO

IF( TW .GT. 0.0d0 )THEN
  WMSE= WMSE / TW
ELSE
  WMSE= UNDEF
ENDIF


END SUBROUTINE MSE

SUBROUTINE UNDEF_MEAN( field,nx,ny, mean)
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER , INTENT(IN) :: nx,ny
REAL(dp), INTENT(IN) :: field(nx,ny)
REAL(dp), PARAMETER  :: UNDEF=-999.0d0
REAL(dp), INTENT(OUT) :: mean
INTEGER   :: ii,jj,n

mean=0.0d0
n=0

DO ii=1,nx
  DO jj=1,ny
    IF( field(ii,jj) .NE. UNDEF)THEN
      mean = mean + field(ii,jj)    
      n=n+1
    ENDIF
  ENDDO
ENDDO

IF( n > 0)THEN
  mean = mean / REAL(n,dp)

ELSE
 
  mean = UNDEF

endif

END SUBROUTINE UNDEF_MEAN



SUBROUTINE COTREC_SOR(uo_field,vo_field,dx,nx,ny,alpha,max_iter,u_field,v_field,forcing,lambda)
!This subroutine is based on the algorithm described by Li et al 1995
!This smooths the motion vector field obtained with some other technique 
!by use a zero divergence constrained minimization.
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , ini_i , ini_j , end_i , end_j , contador , n_points
REAL(dp), PARAMETER :: UNDEF=-999.0d0 
REAL(dp), INTENT(IN) :: uo_field(nx,ny),vo_field(nx,ny),dx
REAL(dp), INTENT(OUT) :: u_field(nx,ny) , v_field(nx,ny)
REAL(dp), INTENT(OUT) :: lambda(nx,ny) , forcing(nx,ny)
REAL(dp)             :: up,um,vp,vm,RESIDUAL(nx,ny),dxsq,MAX_RESIDUAL
REAL(dp)             :: l1,l2,l3,l4 , lambdap,lambdam
INTEGER , INTENT(IN) :: nx,ny
INTEGER              :: ITER
INTEGER , INTENT(IN) :: MAX_ITER 
REAL(dp), INTENT(IN) :: ALPHA
REAL(dp), PARAMETER :: RESIDUAL_THRESHOLD = 1.0d-4
LOGICAL             :: CONVERGE

!Compute forcing. TREC WIND divergence.
!At the boundaries (which in this case can be quite complicated)
!a zero gradient condition is assumed.

DO ii = 1 , nx
   DO jj = 1 , ny

     IF( uo_field(ii,jj) .NE. UNDEF )THEN
      !Compute Ux
      IF(  ii .GT. 1 )THEN
       IF( vo_field(ii-1,jj) .EQ. UNDEF)vm=vo_field(ii,jj)
       IF( vo_field(ii-1,jj) .NE. UNDEF)vm=vo_field(ii-1,jj)
      ELSE
       vm=vo_field(ii,jj)
      ENDIF
      IF( ii .LT. nx )THEN
       IF( vo_field(ii+1,jj) .EQ. UNDEF)vp=vo_field(ii,jj)
       IF( vo_field(ii+1,jj) .NE. UNDEF)vp=vo_field(ii+1,jj)
      ELSE
       vp=vo_field(ii,jj)
      ENDIF

      IF(  jj .GT. 1 )THEN
       IF( uo_field(ii,jj-1) .EQ. UNDEF)um=uo_field(ii,jj)
       IF( uo_field(ii,jj-1) .NE. UNDEF)um=uo_field(ii,jj-1)
      ELSE
       um=uo_field(ii,jj)
      ENDIF
      IF( jj .LT. ny )THEN
       IF( uo_field(ii,jj+1) .EQ. UNDEF)up=uo_field(ii,jj)
       IF( uo_field(ii,jj+1) .NE. UNDEF)up=uo_field(ii,jj+1)
      ELSE
       up=uo_field(ii,jj)
      ENDIF

      forcing(ii,jj)=-2.0d0 * ( (up-um)/(2.0d0*dx) + (vp-vm)/(2.0d0*dx) )

     ELSE
      forcing(ii,jj)=UNDEF

     ENDIF
      !IF( forcing(ii,jj) .NE. UNDEF .AND. forcing(ii,jj) == 0.0d0 )WRITE(*,*)forcing(ii,jj)
   ENDDO
ENDDO

!Initialize  residual and Lambda
!RESIDUAL=0.0d0
lambda=0.0d0
u_field=0.0d0
v_field=0.0d0

!Start iterative procedure.
ITER=0
MAX_RESIDUAL=2.0d0*RESIDUAL_THRESHOLD
dxsq=dx**2
DO WHILE( MAX_RESIDUAL > RESIDUAL_THRESHOLD .AND. ITER < MAX_ITER )
!Compute residual.
MAX_RESIDUAL=0.0d0
RESIDUAL=0.0d0

  DO ii = 1 , nx
    DO jj = 1 , ny

       IF( forcing(ii,jj) .NE. UNDEF )THEN
       !The residual is computed only at the locations where the forcing is defined.

          IF( ii .LT. nx )THEN
             l1=lambda(ii+1,jj)
          ELSE
             l1=0.0d0
          ENDIF
          IF( ii .GT. 1 )THEN
             l2=lambda(ii-1,jj)
          ELSE
             l2=0.0d0
          ENDIF
          IF( jj .LT. ny )THEN
             l3=lambda(ii,jj+1)
          ELSE
             l3=0.0d0
          ENDIF
          IF( jj .GT. 1 )THEN
             l4=lambda(ii,jj-1)
          ELSE
             l4=0.0d0
          ENDIF

          !WRITE(*,*)forcing(ii,jj)     
          RESIDUAL(ii,jj)= l1 + l2 + l3 + l4 - 4*lambda(ii,jj) - dxsq*forcing(ii,jj) 

          IF( ABS(RESIDUAL(ii,jj)) .GT. MAX_RESIDUAL )THEN
              MAX_RESIDUAL=ABS(RESIDUAL(ii,jj))
          ENDIF

       ENDIF

    ENDDO
  ENDDO

  !UDDATE LAMBDA
  lambda = lambda + (ALPHA / 4.0d0) * RESIDUAL

  WRITE(*,*)"COTREC SOR, MAX_RESIDUAL ",MAX_RESIDUAL," ITERATION ",ITER

ITER=ITER+1
ENDDO


!Compute smooth U and smoooth V (equations 9 and 10 of Li et al 1995)
DO ii = 1 , nx
  DO jj = 1 , ny
     IF( uo_field(ii,jj) .NE. UNDEF) THEN
        IF( ii .LT. nx )THEN
          IF( vo_field(ii+1,jj) .NE. UNDEF)THEN 
             vp=vo_field(ii+1,jj)
             lambdap=lambda(ii+1,jj)
          ELSE
             vp=vo_field(ii,jj)
             lambdap=0.0d0
          ENDIF
        ELSE
          vp=vo_field(ii,jj)
          lambdap=0.0d0
        ENDIF
        IF( ii .GT. 1 )THEN
          IF( vo_field(ii-1,jj) .NE. UNDEF)THEN 
             vm=vo_field(ii-1,jj)
             lambdam=lambda(ii-1,jj)
          ELSE
             vm=vo_field(ii,jj)
             lambdam=0.0d0
          ENDIF
        ELSE
          vm=vo_field(ii,jj)
          lambdam=0.0d0
        ENDIF

        v_field(ii,jj) =0.25d0*(2.0d0*vo_field(ii,jj)+vp+vm) +0.5d0*(lambdap-lambdam)/(dx*2.0d0)


        IF( jj .LT. ny )THEN
          IF( uo_field(ii,jj+1) .NE. UNDEF)THEN
             up=uo_field(ii,jj+1)
             lambdap=lambda(ii,jj+1)
          ELSE
             up=uo_field(ii,jj)
             lambdap=0.0d0
          ENDIF
        ELSE
          up=uo_field(ii,jj)
          lambdap=0.0d0
        ENDIF
        IF( jj .GT. 1 )THEN
          IF( uo_field(ii,jj-1) .NE. UNDEF)THEN
             um=uo_field(ii,jj-1)
             lambdam=lambda(ii,jj-1)
          ELSE
             um=uo_field(ii,jj)
             lambdam=0.0d0
          ENDIF
        ELSE
          um=uo_field(ii,jj)
          lambdam=0.0d0
        ENDIF

        u_field(ii,jj) =0.25d0*(2.0d0*uo_field(ii,jj)+up+um) +0.5d0*(lambdap-lambdam)/(dx*2.0d0)

     ENDIF
  ENDDO
ENDDO

WHERE( uo_field == UNDEF )
     v_field = UNDEF 
     u_field = UNDEF
ENDWHERE


END SUBROUTINE COTREC_SOR

!-----------------------------------------------------------------------
! Correlation
!-----------------------------------------------------------------------
SUBROUTINE CORRELATION(A,B,W,N,COR)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)


INTEGER, INTENT(IN) :: N
REAL(dp) , INTENT(IN) :: A(N),B(N),W(N)
REAL(dp) , INTENT(OUT) :: COR

!Weighted correlation, A and B are two variables and W are the weights that 
!corresponds to each element of A and B.

REAL(dp) :: STDA , STDB , COV

  CALL STDEV(A,W,N,STDA)
  CALL STDEV(B,W,N,STDB)
  CALL COVAR(A,B,W,N,COV)

  COR = COV/STDA/STDB

END SUBROUTINE CORRELATION

!-----------------------------------------------------------------------
! Covariance
!-----------------------------------------------------------------------
SUBROUTINE COVAR(A,B,W,N,COV)
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = kind(1.0e0)
  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER,INTENT(IN) :: N
  REAL(dp),INTENT(IN) :: A(N)
  REAL(dp),INTENT(IN) :: B(N)
  REAL(dp),INTENT(IN) :: W(N)
  REAL(dp),INTENT(OUT) :: COV

  REAL(dp) :: MEANA,MEANB

  MEANA=SUM(W*A)/SUM(W)  !Weighted mean.
  MEANB=SUM(W*B)/SUM(W)  !Weighted mean

  COV = SUM( W*(A-MEANA)*(B-MEANB) )/SUM(W) !SUM( W*(A-MEANA)*(B-MEANB) )  !Weigthed covariance

  RETURN
END SUBROUTINE COVAR

!-----------------------------------------------------------------------
! Standard deviation
!-----------------------------------------------------------------------
SUBROUTINE STDEV(A,W,N,STD)
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = kind(1.0e0)
  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER,INTENT(IN) :: N
  REAL(dp),INTENT(IN) :: A(N)
  REAL(dp),INTENT(IN) :: W(N)
  REAL(dp),INTENT(OUT) :: STD

  !STD = SQRT( SUM ((A - SUM(A)/REAL(N,dp))**2)/REAL(N,dp) )

  STD = SQRT( SUM(W*(A(:) - SUM(W*A)/SUM(W) )**2)/SUM(W) )

  RETURN
END SUBROUTINE STDEV


SUBROUTINE BARNES(x,y,z,xo,yo,xr,yr,nx,ny,no,npassmax,zo)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER, INTENT(IN) :: nx , ny !Matrix dimensions.
INTEGER, INTENT(IN) :: no      !Number of observations
REAL(DP),INTENT(IN) :: xo(nx,ny),yo(nx,ny)
INTEGER,INTENT(IN) ::  npassmax !Number of filter passes.
REAL(DP),INTENT(IN) :: x(no),y(no),z(no) !Location and value of input observations.
REAL(DP),INTENT(IN) :: xr , yr
REAL(DP),INTENT(OUT) :: zo(nx,ny)
!Internal variables
REAL(DP) :: xp(nx*ny) , yp(nx*ny) , zp(nx*ny) , zpp(nx*ny) , zpa(nx*ny) ,zpb(nx*ny) , zpn(nx*ny)
REAL(DP) :: xr2 , yr2 , dx , dy , w  , wtot , wmax , suma
INTEGER  :: i , j , np , ni , npass , n
REAL(DP), PARAMETER :: cp=0.3

! DESCRIPTION
! $
! Barnes objective analysis:  a successive corrections scheme
! which is an easy-to-use cousin of kriging, classic objective
! analysis, or optimal interpolation. 
! $
! Returns matrix ZO containing elements corresponding to the elements
! of matrices XO and YO and determined by Barnes objective analysis.
! X and Y are vectors containing the input data coordinates, and Z
! is a vector with the input data.
! $
! XR and YR are the Barnes smoothing length scales in the x and
! y senses.  These remain fixed throughout the iterations (as
! recommended in the later Barnes papers).  Optional NPASS sets the 
! number of Barnes iterations (defaults to 3, which is usually
! sufficient).
! $
! HISTORY
! Translated to Fortran 90, 02-Nov-15, J. Ruiz
! Bug fix, 8-Sep-10, Carlos Carrillo 
! Original, 10-Aug-09, S. Pierce
! $
! REFERENCES
!$
! Barnes, S. L. (1994) Applications of the Barnes objective analysis
! scheme.  Part I:  effects of undersampling, wave position, and station
! randomness.  J. of Atmos. and Oceanic Tech., 11, 1433-1448.
! $
! Barnes, S. L. (1994) Applications of the Barnes objective analysis
! scheme.  Part II:  Improving derivative estimates.  J. of Atmos. and 
! Oceanic Tech., 11, 1449-1458.
! $
! Barnes, S. L. (1994) Applications of the Barnes objective analysis
! scheme.  Part III:  Tuning for minimum error.  J. of Atmos. and Oceanic 
! Tech., 11, 1459-1479.
! $
! Daley, R. (1991) Atmospheric data analysis, Cambridge Press, New York.
! Section 3.6.

!From matrix to row vector.
WRITE(*,*) " INPUT CHECK BARNES SUBROUTINE"
WRITE(*,*) " X , Y , Z ",X(1),Y(1),Z(1)
WRITE(*,*) " XO YO ",XO(1,1),YO(1,1)
WRITE(*,*) " XR , YR ",xr,yr
WRITE(*,*) " NX NY ",nx,ny
WRITE(*,*) " NO    ",no
WRITE(*,*) " NPASSMAX ",npassmax

DO i=1,nx
   DO j=1,ny
    xp((i-1)*nx+j) = xo(i,j)
    yp((i-1)*nx+j) = yo(i,j)
    zp((i-1)*nx+j) = 0.0d0
   ENDDO
ENDDO


xr2=xr**2
yr2=yr**2

!FIRST GUESS (WEIGTHED AVERAGE) 
!This step is the first pass of the filter.
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(dx,dy,w,np,wtot)
!$OMP DO
DO np=1,nx*ny
     zp(np)=0.0d0
     wtot=0.0d0
     DO i=1,no
       dx=(x(i)-xp(np))**2
       dy=(y(i)-yp(np))**2
       w=EXP(-dx/xr2-dy/yr2 )
       zp(np)=zp(np)+z(i)*w
       wtot=wtot+w
     ENDDO
     zp(np)=zp(np)/wtot
ENDDO
!$END OMP DO
!$OMP END PARALLEL


zpp=zp
zpb=0.0d0

!Update localization parameters.
!xr2=(xr*cp)**2
!yr2=(yr*cp)**2


DO npass=2,npassmax

!Update localization parameters.
xr2=(xr/npass)**2
yr2=(yr/npass)**2

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,dx,dy,w,np,ni,suma,wtot)
!$OMP DO
   DO ni=1,no
    IF( mod(npass,2) == 0)THEN
      wtot=0.0d0
      suma=0.0d0
      DO i=1,no
       dx=(x(i)-x(ni))**2
       dy=(y(i)-y(ni))**2
       w=exp(-dx/xr2-dy/yr2)
       suma=suma+(z(i)-zpb(i))*w
       wtot=wtot+w
      ENDDO
       zpa(ni)=zpb(ni)+suma/wtot
    ELSE
      wtot=0.0d0
      suma=0.0d0
      DO i=1,no
       dx=(x(i)-x(ni))**2
       dy=(y(i)-y(ni))**2
       w=exp(-dx/xr2-dy/yr2)
       suma=suma+(z(i)-zpa(i))*w
       wtot=wtot+w
      ENDDO
       zpb(ni)=zpa(ni)+suma/wtot

     ENDIF
  ENDDO
!$END OMP DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,dx,dy,w,np,ni,suma,wtot)
!$OMP DO

  DO np=1,nx*ny

      suma=0.0d0
      wtot=0.0d0
      wmax=0.0d0
      IF( mod(npass,2)==0 )THEN
        DO i=1,no
         dx=(x(i)-xp(np))**2
         dy=(y(i)-yp(np))**2
         w=exp(-dx/xr2-dy/yr2)
         suma=suma+(z(i)-zpa(i))*w
         wtot=wtot+w
         IF( w > wmax)THEN
           wmax=w
         ENDIF
        ENDDO
      ELSE
        DO i=1,no
         dx=(x(i)-xp(np))**2
         dy=(y(i)-yp(np))**2
         w=exp(-dx/xr2-dy/yr2)
         suma=suma+(z(i)-zpb(i))*w
         wtot=wtot+w
         IF( w > wmax)THEN
           wmax=w
         ENDIF
        ENDDO
       ENDIF
       !Modified Barnes Analysis. The update is greather where we have closer observation.
       zp(np)=zp(np)+wmax*suma/wtot

  ENDDO
!$END OMP DO
!$OMP END PARALLEL

    IF ( mod(npass,2)==0 )THEN
       zpn=zp
    ELSE
       zpp=zp
    ENDIF

WRITE(*,*)"rms z adjustment after pass ",npass ," is ",SUM((zpn-zpp)**2)/REAL(nx*ny,dp)   
ENDDO


DO i=1,nx
 DO j=1,ny
    zo(i,j) = zp((i-1)*nx+j) 
 ENDDO
ENDDO



END SUBROUTINE BARNES

SUBROUTINE GAUSSIAN_FILTER(field0,field1,dx,sigma,nx,ny)

IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj 
REAL(dp), PARAMETER :: UNDEF=-999.0d0 
REAL(dp), INTENT(IN) :: field0(nx,ny),dx,sigma
REAL(dp), INTENT(OUT):: field1(nx,ny)
INTEGER , INTENT(IN) :: nx,ny
REAL(dp)             :: tmp , suma , sumaw , dist

WRITE(*,*) " INSIDE GAUSSIAN FILTER "
WRITE(*,*) " DELTA X                  ",dx
WRITE(*,*) " SIGMA                    ",sigma
WRITE(*,*) " NX ,  NY                 ",nx,ny

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,dist,tmp,suma,sumaw) 
!$OMP DO 
DO ii = 1,nx
  DO jj = 1,ny
        !COMPUTE AVERAGED VALUE.
        suma=0.0d0
        sumaw=0.0d0
        DO iii=1,nx
          DO jjj=1,ny
               dist=( ((ii-iii)*dx) ** 2 + ((jj-jjj)*dx) ** 2 )/( sigma**2)
               IF( dist <= 10)THEN
                 tmp=exp(-dist)
                 suma=suma+field0(iii,jjj)*tmp
                 sumaw=sumaw+tmp
               ENDIF
          ENDDO
        ENDDO

  IF( sumaw > 0.0d0)THEN
   field1(ii,jj)=suma/sumaw
  ELSE
   field1(ii,jj)=UNDEF
  ENDIF

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL


WRITE(*,*)"FINISH GAUSSIAN FILTER"
END SUBROUTINE GAUSSIAN_FILTER

SUBROUTINE FILTER_OUTLIER(fieldu,fieldv,fieldoutu,fieldoutv,nx,ny,threshold,box_size)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj
REAL(dp), PARAMETER :: UNDEF=-999.0d0, pi=3.141592
REAL(dp), INTENT(IN) :: fieldu(nx,ny),fieldv(nx,ny),threshold
REAL(dp), INTENT(OUT):: fieldoutu(nx,ny), fieldoutv(nx,ny)
INTEGER , INTENT(IN) :: nx,ny,box_size
REAL(dp)             :: sumaw,v_mod,u_dif,v_dif,us,vs,dif_mod
INTEGER              :: ini_i,end_i,ini_j,end_j,np

fieldoutu=fieldu
fieldoutv=fieldv
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,sumaw, &
!$OMP&             ini_i,end_i,ini_j,end_j,np,dif_mod,us,vs,u_dif,v_dif,v_mod) 
!$OMP DO 
DO ii = 1,nx
  DO jj = 1,ny
        !COMPUTE LOCAL AVERAGE
        ini_i=ii-box_size
        end_i=ii+box_size
        ini_j=jj-box_size
        end_j=jj+box_size

        if(ini_i < 1)ini_i=1
        if(end_i > nx)end_i=nx
        if(ini_j < 1)ini_j=1
        if(end_j > ny)end_j=ny

        !Compute sample size
        np=(end_i-ini_i+1)*(end_j-ini_j+1)

        !sumasq=0.0d0
        sumaw=0.0d0
        us=0.0d0
        vs=0.0d0
        do iii=ini_i,end_i
          do jjj=ini_j,end_j
            if( fieldu(iii,jjj) .ne. UNDEF )then
             us=us+fieldu(iii,jjj)
             vs=vs+fieldv(iii,jjj)
             sumaw=sumaw+1.0d0
            endif
          enddo
        enddo
        if( sumaw == 0.0d0 )cycle
        us=us/real(sumaw)
        vs=vs/real(sumaw)
        u_dif=0.0d0
        v_dif=0.0d0
        v_mod=0.0d0
        u_dif=ABS(fieldu(ii,jj)-us)
        v_dif=ABS(fieldv(ii,jj)-vs)
        dif_mod=sqrt(u_dif**2+v_dif**2)
!        v_mod=sqrt(fieldu(ii,jj)**2+fieldv(ii,jj)**2)
        if(ABS(dif_mod)>threshold)then
            fieldoutu(ii,jj)=UNDEF
            fieldoutv(ii,jj)=UNDEF
        endif

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL
WRITE(*,*)"FINISH FILTER OUTLIER"
END SUBROUTINE FILTER_OUTLIER



!!!
!SUBROUTINE FILTER_OUTLIER(fieldu, fieldv, fieldoutu, fieldoutv, nx, ny, threshold, box_size, angle)
!IMPLICIT NONE
!INTEGER, PARAMETER :: sp = kind(1.0e0)
!INTEGER, PARAMETER :: dp = kind(1.0d0)
!INTEGER            :: ii , jj , iii , jjj
!REAL(dp), PARAMETER :: UNDEF=-999.0d0, pi=3.141592
!REAL(dp)             :: v_mod(,), vv_mod(,)
!REAL(dp), INTENT(IN) :: fieldu(nx,ny),fieldv(nx,ny),threshold
!REAL(dp), INTENT(OUT):: fieldoutu(nx,ny), fieldoutv(nx,ny)
!INTEGER , INTENT(IN) :: nx,ny,box_size, angle
!REAL(dp)             :: suma,sumaw,v_mod,vv_mod,local_angle,us,vs
!INTEGER              :: ini_i,end_i,ini_j,end_j,np

!fieldoutu=fieldu
!fieldoutv=fieldv
!!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,suma,sumaw, &
!!$OMP&                          ini_i,end_i,ini_j,end_j,np) 
!!$OMP DO 
!DO ii = 1,nx
!  DO jj = 1,ny
        !COMPUTE LOCAL AVERAGE
!        ini_i=ii-box_size
!        end_i=ii+box_size
!        ini_j=jj-box_size
!        end_j=jj+box_size

!        if(ini_i < 1)ini_i=1
!        if(end_i > nx)end_i=nx
!        if(ini_j < 1)ini_j=1
!        if(end_j > ny)end_j=ny

        !Compute sample size
!        np=(end_i-ini_i+1)*(end_j-ini_j+1)

!        suma=0.0d0
        !sumasq=0.0d0
!        sumaw=0.0d0
!        v_mod=0.0d0
!        us=0.0d0
!        vs=0.0d0
!        do iii=ini_i,end_i
!          do jjj=ini_j,end_j
!            if( fieldu(iii,jjj) .ne. UNDEF )then
!            v_mod=sqrt(fieldu(iii,jjj)**2+fieldv(iii,jjj)**2)
!            us=us+fieldu(iii,jjj)
!            vs=vs+fieldv(iii,jjj)
!            suma=suma+v_mod
!            sumaw=sumaw+1.0d0
!            endif
!          enddo
!        enddo
!        if( sumaw == 0.0d0 )cycle
!        suma=suma/real(sumaw)
!        us=us/real(sumaw)
!        vs=vs/real(sumaw)
!        vv_mod=0.0d0
!        vv_mod=sqrt(fieldu(ii,jj)**2+fieldv(ii,jj)**2)
!        local_angle=acos((fieldu(ii,jj)*us+fieldv(ii,jj)*vs)/real((sqrt(fieldu(ii,jj)**2+fieldv(ii,jj)**2)*(sqrt(us**2+vs**2)))))
!        local_angle=local_angle*180/real(3.1415926)
!        if(ABS(vv_mod - suma)>threshold .or. local_angle>angle)then
!            fieldoutu(ii,jj)=UNDEF
!            fieldoutv(ii,jj)=UNDEF
!        endif

!  ENDDO
!ENDDO

!!$END OMP DO
!!$OMP END PARALLEL


!WRITE(*,*)"FINISH FILTER OUTLIER"
!END SUBROUTINE FILTER_OUTLIER
!!!



!SUBROUTINE FILTER_OUTLIER(field0,field1,nx,ny,threshold,box_size)
!IMPLICIT NONE
!INTEGER, PARAMETER :: sp = kind(1.0e0)
!INTEGER, PARAMETER :: dp = kind(1.0d0)
!INTEGER            :: ii , jj , iii , jjj 
!REAL(dp), PARAMETER :: UNDEF=-999.0d0 
!REAL(dp), INTENT(IN) :: field0(nx,ny),threshold
!REAL(dp), INTENT(OUT):: field1(nx,ny)
!INTEGER , INTENT(IN) :: nx,ny,box_size
!REAL(dp)             :: suma,sumasq,sumaw
!INTEGER              :: ini_i,end_i,ini_j,end_j,np

!field1=field0

!!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,suma,sumasq,sumaw, &
!!$OMP&                          ini_i,end_i,ini_j,end_j,np) 
!!$OMP DO 
!!DO ii = 1,nx
!!  DO jj = 1,ny
!        !COMPUTE LOCAL AVERAGE
!!        ini_i=ii-box_size
!!        end_i=ii+box_size
!!        ini_j=jj-box_size
!!        end_j=jj+box_size
!
!!        if(ini_i < 1)ini_i=1
!!        if(end_i > nx)end_i=nx
!!        if(ini_j < 1)ini_j=1
!!        if(end_j > ny)end_j=ny
!
!        !Compute sample size
!!        np=(end_i-ini_i+1)*(end_j-ini_j+1)
!        
!!        suma=0.0d0
!!        sumaw=0.0d0
!!        do iii=ini_i,end_i
!!          do jjj=ini_j,end_j
!!            if( field0(iii,jjj) .ne. UNDEF )then
!!            suma=suma+field0(iii,jjj)
!!            sumaw=sumaw+1.0d0
!!            endif
!!          enddo
!!        enddo
!!        if( sumaw == 0.0d0 )cycle
!
!!        suma=suma/real(sumaw)
! 
!!        if( ABS(field0(ii,jj) - suma ) > threshold )then
!!            field1(ii,jj)=UNDEF
!!        endif
!
!
!!  ENDDO
!!ENDDO
!
!!$END OMP DO
!!$OMP END PARALLEL
!
!
!!WRITE(*,*)"FINISH FILTER OUTLIER"
!!END SUBROUTINE FILTER_OUTLIER
!
!
!!-----------------------------------------------------------------------
!! Damp lateral boundary
!!-----------------------------------------------------------------------
SUBROUTINE damp_latbnd(spec_bdy_width,nx,ny,var,vard)
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = kind(1.0e0)
  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER(dp),INTENT(in) :: spec_bdy_width
  INTEGER(dp),INTENT(in) :: nx,ny
  REAL(dp),INTENT(in)  :: var(nx,ny)
  REAL(dp),INTENT(out) :: vard(nx,ny)
  INTEGER :: i,j,k
  INTEGER :: ist,ied,jst,jed

  vard=var

  ist = 1
  ied = nx - 1
  jst = 1
  jed = ny - 1
    DO i = ist, spec_bdy_width
      vard(i,1:ny) = var(i,1:ny) * real(i-1,dp) / real(spec_bdy_width)
    END DO
    DO i = ied-spec_bdy_width+1, ied
      vard(i,1:ny) = var(i,1:ny) * real(ied-i,dp) / real(spec_bdy_width)
    END DO
    DO j = jst, spec_bdy_width
      vard(1:nx,j) = var(1:nx,j) * real(j-1,dp) / real(spec_bdy_width)
    END DO
    DO j = jed-spec_bdy_width+1, jed
      vard(1:nx,j) = var(1:nx,j) * real(jed-j,dp) / real(spec_bdy_width)
    END DO


  RETURN
END SUBROUTINE damp_latbnd




END MODULE motion_vectors



