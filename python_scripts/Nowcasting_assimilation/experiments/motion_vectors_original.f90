! THESE ROUTINES ARE FOR A MOTION VECTOR COMPUTATION USING
! PYTHON-FORTRAN 
! COMPILE USING: f2py -c -m test /home/jruiz/PYTHON/motion_vectors.f90
! FOR PARALLEL FUNCTIONALITY f2py -c -lgomp --f90flags="-fopenmp -lgomp"
! -m motion_vectors motion_vectors_parallel.f90 

SUBROUTINE
MOTION_VECTOR(field_t0,field_t1,dt,dx,sigma,sigma_threshold,nx,ny,desp_max,u_motion,v_motion,max_corr
&
       &      ,aux_output,aux_inputi,aux_inputj)


IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , ini_i , ini_j , end_i ,
end_j , contador , n_points

REAL(dp), PARAMETER :: UNDEF=-999.0d0 , min_box_fraction = 0.10d0 ,
min_dbz=0.0d0
INTEGER , PARAMETER :: METRIC_OPTION = 1 ! Metric option for similarity
quantification (1 = correlation , ...)

REAL(dp), INTENT(IN) ::
field_t0(nx,ny),field_t1(nx,ny),dt,dx,sigma,sigma_threshold
INTEGER , INTENT(IN) :: nx,ny,desp_max
INTEGER , INTENT(IN) :: aux_inputi , aux_inputj

REAL(dp), INTENT(OUT) :: u_motion(nx,ny) , v_motion(nx,ny) ,
max_corr(nx,ny) 
REAL(dp), INTENT(OUT) :: aux_output(2*desp_max+ 1,2*desp_max+1)
INTEGER  :: box_size , tmp_dim(1)

REAL(dp),ALLOCATABLE :: corr_matrix(:,:) , vector_1(:) , vector_2(:) 

REAL(dp),ALLOCATABLE :: correlation_weigths(:,:) , vector_cw(:)
!Gaussian weights for weighted correlation.
REAL(dp) :: box_fraction , dist , tmp_corr  ,  despi ,despj,weigth
INTEGER  :: maxi , maxj  , i_point

REAL(dp) , PARAMETER :: correlation_threshold = 0.3  !

INTEGER , PARAMETER :: MOTION_VECTOR_OPTION=1  !0 - Maximum corr , 1
Average sourrounding max corr, 2 - Average over entire displacements.
INTEGER , PARAMETER :: MOTION_VECTOR_WEIGTH=0  !0 - Uniform weigth , 1 -
Gaussian weigths.

WRITE(*,*) " INSIDE W MOTION VECTOR ROUTINE "
WRITE(*,*) " FIRST REF FIELD AT T 0 = ",field_t0(1,1)
WRITE(*,*) " FIRST REF FIELD AT T 1 = ",field_t1(1,1)
WRITE(*,*) " DELTA T                  ",dt
WRITE(*,*) " DELTA X                  ",dx
WRITE(*,*) " SIGMA                    ",sigma
WRITE(*,*) " SIGMA_THRESHOLD          ",sigma_threshold
WRITE(*,*) " NX ,  NY                 ",nx,ny
WRITE(*,*) " DESP MAX                 ",desp_max

!First find the box size based on sigma and the threshold. In this
!implementation
!the box size is computed from sigma and sigma_threshold.

box_size= CEILING ( sigma_threshold * sigma / dx )

ALLOCATE(corr_matrix(desp_max*2+1,desp_max*2+1) ,
vector_1((box_size*2+1)**2) , vector_2((box_size*2+1)**2)  )
!ALLOCATE(vector_1((box_size*2+1)**2) , vector_2((box_size*2+1)**2)  )

!Compute weigths for correlation.
ALLOCATE( correlation_weigths(box_size*2+1,box_size*2+1) ,
vector_cw((box_size*2+1)**2) )


SELECT CASE ( MOTION_VECTOR_WEIGTH )

CASE( 0 ) !Uniform weigths

   correlation_weigths=1.0d0

CASE( 1 ) !Gaussian weigths

   DO ii=-box_size,box_size
    DO jj=-box_size,box_size
      dist=( (ii*dx) ** 2 + (jj*dx) ** 2 )
      correlation_weigths( ii+box_size+1 , jj+box_size+1 ) = exp (-(
dist / (sigma**2) ) )
    ENDDO
  ENDDO

END SELECT


!Initialize u_motion and v_motion with UNDEF values.
u_motion=undef
v_motion=undef
max_corr=-1.0d0


!$OMP PARALLEL DEFAULT(SHARED)
!PRIVATE(ii,jj,ini_i,end_i,ini_j,end_j,iii,jjj,corr_matrix & 
!$OMP&
!,vector_1,vector_2,vector_cw,n_points,tmp_dim,contador,box_fraction
!&
!$OMP&           ,despi,despj,tmp_corr,weigth,maxi,maxj) 
!$OMP DO 
DO ii = 1,nx


  WRITE(*,*)"Processing row number ",ii
  DO jj = 1,ny


        corr_matrix=0.0d0

        ini_i=max(ii-box_size,1+desp_max)
        end_i=min(ii+box_size,nx-desp_max)
        ini_j=max(jj-box_size,1+desp_max)
        end_j=min(jj+box_size,ny-desp_max)

        n_points=(end_i-ini_i+1)*(end_j-ini_j+1)


        IF ( n_points   <= 1 )CYCLE   !No podemos calcular la velocidad
en este caso.
 
          contador=0
          DO iii=ini_i,end_i
            DO jjj=ini_j,end_j
               IF( field_t0(iii,jjj) > min_dbz )THEN
                 contador = contador + 1
               ENDIF
            ENDDO
          ENDDO


        box_fraction=REAL(contador,dp)/REAL(n_points,dp)

        IF ( box_fraction < min_box_fraction )CYCLE   !No podemos
calcular la velocidad en este caso.

        tmp_dim(1)=n_points

        vector_1(1:n_points)=RESHAPE(field_t0(ini_i:end_i,ini_j:end_j),tmp_dim)
        vector_cw(1:n_points)=RESHAPE(correlation_weigths(ini_i-ii+box_size+1:end_i-ii+box_size+1,ini_j-jj+box_size+1:
&
     &                                 end_j-jj+box_size+1),tmp_dim)
        vector_cw(1:n_points)=vector_cw(1:n_points)/SUM(vector_cw(1:n_points))

        maxi=0
        maxj=0


        DO iii= -desp_max,desp_max
          DO jjj= -desp_max,desp_max

            vector_2(1:n_points)=RESHAPE(field_t1(ini_i+iii:end_i+iii,ini_j+jjj:end_j+jjj),tmp_dim)
            !Compute the correlation between field_t0 and field_t1
            CALL CORRELATION( vector_1(1:n_points), vector_2(1:n_points)
            C, vector_cw(1:n_points) , n_points , tmp_corr )

            corr_matrix(iii+desp_max+1,jjj+desp_max+1)=tmp_corr
            IF( tmp_corr >= max_corr(ii,jj) )THEN
                max_corr(ii,jj)=tmp_corr
                maxj=jjj
                maxi=iii
            ENDIF

           ENDDO
         ENDDO
 
     !!!!! FOR DEBUG
     IF( ii==aux_inputi .AND. jj==aux_inputj )THEN
         aux_output=corr_matrix
     ENDIF 
     !!!!!!!


     SELECT CASE ( MOTION_VECTOR_OPTION )

     CASE( 0 ) !Maximum correlation criteria. 

       despi=REAL(maxi,dp)
       despj=REAL(maxj,dp)
 
       IF( max_corr(ii,jj) >= 0.0d0 )THEN
       u_motion(ii,jj)= despj*dx/dt
       v_motion(ii,jj)= despi*dx/dt
       ENDIF

     CASE( 1 ) !Mean correlation around the maximum value
         !Find the weighted location for the minimum.
         despi=0.0d0
         despj=0.0d0
         weigth=0.0d0
         DO iii=maxi-1,maxi+1
           DO jjj=maxj-1,maxj+1 
             IF( iii <= desp_max  .AND. iii >= -desp_max .AND. jjj <=
desp_max .AND. jjj >= -desp_max )THEN
                IF( corr_matrix( iii+desp_max+1 , jjj+desp_max + 1) >=
0)THEN
                    tmp_corr = corr_matrix( iii+desp_max+1 ,
jjj+desp_max + 1)**2
                    weigth=weigth + tmp_corr
                    despi = despi + tmp_corr *iii
                    despj = despj + tmp_corr *jjj
                ENDIF
             ENDIF                
           ENDDO
         ENDDO

         IF( weigth > 0.0)THEN
           u_motion(ii,jj)= despj*dx/(dt*weigth)
           v_motion(ii,jj)= despi*dx/(dt*weigth)
         ENDIF


     CASE( 2 ) !Mean correlation using the highest correlation points.
        weigth=0.0d0
        despi=0.0d0
        despj=0.0d0
        DO iii= -desp_max,desp_max 
           DO jjj= -desp_max,desp_max
             tmp_corr = corr_matrix( iii + desp_max + 1, jjj + desp_max + 1 )
             IF( max_corr(ii,jj) - tmp_corr < correlation_threshold .AND. tmp_corr > 0.0 )THEN
                 tmp_corr = tmp_corr ** 2
                 weigth=weigth+tmp_corr
                 despi=despi + iii * (tmp_corr )
                 despj=despj + jjj * (tmp_corr )
             ENDIF
           ENDDO
        ENDDO

         !Get the location of the maximum correlation within
         !corr_matrix
         IF( weigth > 0.0)THEN
           u_motion(ii,jj)= despj*dx/(dt*weigth)
           v_motion(ii,jj)= despi*dx/(dt*weigth)
         ENDIF

     END SELECT

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL

DEALLOCATE(corr_matrix , vector_1 , vector_2  )
DEALLOCATE( correlation_weigths , vector_cw )

WRITE(*,*)"FINISH MOTION VECTOR COMPUTATION"
END SUBROUTINE MOTION_VECTOR


SUBROUTINE
HORN_SCHUNCK(field_t0,field_t1,dt,dx,nx,ny,u_motion,v_motion,alpha,max_iter)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , ini_i , ini_j , end_i ,
end_j , contador , n_points
REAL(dp), PARAMETER :: UNDEF=-999.0d0 , min_box_fraction = 0.10d0 ,
min_dbz=0.0d0
INTEGER , PARAMETER :: METRIC_OPTION = 1 ! Metric option for similarity
quantification (1 = correlation , ...)
REAL(dp), INTENT(IN) :: field_t0(nx,ny),field_t1(nx,ny),dt,dx
INTEGER , INTENT(IN) :: nx,ny
REAL(dp), INTENT(OUT) :: u_motion(nx,ny) , v_motion(nx,ny)
REAL(dp)              :: u_tmp(nx,ny)    , v_tmp(nx,ny)
REAL(dp)              :: Ix(nx,ny) , Iy(nx,ny) , It(nx,ny) , Ixsq(nx,ny)
, Iysq(nx,ny)
INTEGER :: ITER
INTEGER , INTENT(IN) :: MAX_ITER
REAL(dp), INTENT(IN) :: ALPHA
REAL(dp), PARAMETER :: UPDATE_THRESHOLD = 0.1d0
REAL(dp)            :: ALPHASQ
LOGICAL             :: CONVERGE
REAL(dp)            :: um,vm , tmp , UPDATE

WRITE(*,*) " INSIDE HORN SCHUNCK "
WRITE(*,*) " FIRST REF FIELD AT T 0 = ",field_t0(1,1)
WRITE(*,*) " FIRST REF FIELD AT T 1 = ",field_t1(1,1)
WRITE(*,*) " DELTA T                  ",dt
WRITE(*,*) " DELTA X                  ",dx
WRITE(*,*) " NX ,  NY                 ",nx,ny
!Initialize u_motion and v_motion with UNDEF values.
u_motion=0.0d0
v_motion=0.0d0
u_tmp=u_motion
v_tmp=v_motion
!First compute local trend and spatial derivative.
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
    DO ii = 2 , nx - 1
     DO jj = 2 , ny - 1

       um=0.25*(u_motion(ii+1,jj)+u_motion(ii-1,jj)+u_motion(ii,jj+1)+u_motion(ii,jj-1))
       vm=0.25*(v_motion(ii+1,jj)+v_motion(ii-1,jj)+v_motion(ii,jj+1)+v_motion(ii,jj-1))

       tmp= -( Ix(ii,jj)*um + Iy(ii,jj)*vm + It(ii,jj) )/(ALPHASQ + Ixsq(ii,jj) + Iysq(ii,jj))

       u_tmp(ii,jj)= um + Ix(ii,jj)*tmp
       v_tmp(ii,jj)= vm + Iy(ii,jj)*tmp

       UPDATE=UPDATE + ABS(u_motion(ii,jj)-u_tmp(ii,jj)) + ABS(v_motion(ii,jj)-v_tmp(ii,jj))

     ENDDO
    ENDDO

    UPDATE=UPDATE/REAL(nx*ny,dp)
    IF( UPDATE <= UPDATE_THRESHOLD)THEN
!      CONVERGE=.TRUE. !We reach convergence! We will stop iteration
!      process.
    ENDIF

    !Boundary conditions (zero derivative at the boundaries.
    u_motion=u_tmp
    v_motion=v_tmp

    u_motion(1,:)=u_motion(2,:)
    u_motion(nx,:)=u_motion(nx-1,:)
    u_motion(:,1)=u_motion(:,2)
    u_motion(:,ny)=u_motion(:,ny-1)

    v_motion(1,:)=v_motion(2,:)
    v_motion(nx,:)=v_motion(nx-1,:)
    v_motion(:,1)=v_motion(:,2)
    v_motion(:,ny)=v_motion(:,ny-1)

ITER=ITER+1
WRITE(*,*)"ITERATION ",ITER," MOTION VECTOR UPDATE ",UPDATE
ENDDO


u_motion=u_motion*dx/dt
v_motion=v_motion*dx/dt


WRITE(*,*)"FINISH MOTION VECTOR COMPUTATION"
END SUBROUTINE HORN_SCHUNCK


SUBROUTINE
COTREC_SOR(uo_field,vo_field,dx,nx,ny,u_field,v_field,forcing,lambda)
!This subroutine is based on the algorithm described by Li et al 1995
!This smooths the motion vector field obtained with some other technique 
!by use a zero divergence constrained minimization.
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , ini_i , ini_j , end_i ,
end_j , contador , n_points
REAL(dp), PARAMETER :: UNDEF=-999.0d0 
REAL(dp), INTENT(IN) :: uo_field(nx,ny),vo_field(nx,ny),dx
REAL(dp), INTENT(OUT) :: u_field(nx,ny) , v_field(nx,ny)
REAL(dp), INTENT(OUT) :: lambda(nx,ny) , forcing(nx,ny)
REAL(dp)             :: up,um,vp,vm,RESIDUAL(nx,ny),dxsq,MAX_RESIDUAL
REAL(dp)             :: l1,l2,l3,l4 , lambdap,lambdam
INTEGER , INTENT(IN) :: nx,ny
INTEGER              :: ITER
INTEGER , PARAMETER :: MAX_ITER = 1000
REAL(dp), PARAMETER :: ALPHA=1.0d0
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
       IF( uo_field(ii-1,jj) .EQ. UNDEF)um=uo_field(ii,jj)
       IF( uo_field(ii-1,jj) .NE. UNDEF)um=uo_field(ii-1,jj)
      ELSE
       um=uo_field(ii,jj)
      ENDIF
      IF( ii .LT. nx )THEN
       IF( uo_field(ii+1,jj) .EQ. UNDEF)up=uo_field(ii,jj)
       IF( uo_field(ii+1,jj) .NE. UNDEF)up=uo_field(ii+1,jj)
      ELSE
       up=uo_field(ii,jj)
      ENDIF

      IF(  jj .GT. 1 )THEN
       IF( vo_field(ii,jj-1) .EQ. UNDEF)vm=vo_field(ii,jj)
       IF( vo_field(ii,jj-1) .NE. UNDEF)vm=vo_field(ii,jj-1)
      ELSE
       vm=vo_field(ii,jj)
      ENDIF
      IF( jj .LT. ny )THEN
       IF( vo_field(ii,jj+1) .EQ. UNDEF)vp=vo_field(ii,jj)
       IF( vo_field(ii,jj+1) .NE. UNDEF)vp=vo_field(ii,jj+1)
      ELSE
       vp=vo_field(ii,jj)
      ENDIF

      forcing(ii,jj)=-2.0d0 * ( (up-um)/(2.0d0*dx) + (vp-vm)/(2.0d0*dx)
)

     ELSE
      forcing(ii,jj)=UNDEF

     ENDIF
      !IF( forcing(ii,jj) .NE. UNDEF .AND. forcing(ii,jj) == 0.0d0
      !)WRITE(*,*)forcing(ii,jj)
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
       !The residual is computed only at the locations where the forcing
       !is defined.

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
          RESIDUAL(ii,jj)= l1 + l2 + l3 + l4 - 4*lambda(ii,jj) -
dxsq*forcing(ii,jj) 

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
          IF( uo_field(ii+1,jj) .NE. UNDEF)THEN 
             up=uo_field(ii+1,jj)
             lambdap=lambda(ii+1,jj)
          ELSE
             up=uo_field(ii,jj)
             lambdap=0.0d0
          ENDIF
        ELSE
          up=uo_field(ii,jj)
          lambdap=0.0d0
        ENDIF
        IF( ii .GT. 1 )THEN
          IF( uo_field(ii-1,jj) .NE. UNDEF)THEN 
             um=uo_field(ii-1,jj)
             lambdam=lambda(ii-1,jj)
          ELSE
             um=uo_field(ii,jj)
             lambdam=0.0d0
          ENDIF
        ELSE
          um=uo_field(ii,jj)
          lambdam=0.0d0
        ENDIF

        u_field(ii,jj) =0.25d0*(2.0d0*uo_field(ii,jj)+up+um)
+0.5d0*(lambdap-lambdam)/(dx*2.0d0)


        IF( jj .LT. ny )THEN
          IF( vo_field(ii,jj+1) .NE. UNDEF)THEN
             vp=vo_field(ii,jj+1)
             lambdap=lambda(ii,jj+1)
          ELSE
             vp=vo_field(ii,jj)
             lambdap=0.0d0
          ENDIF
        ELSE
          vp=vo_field(ii,jj)
          lambdap=0.0d0
        ENDIF
        IF( jj .GT. 1 )THEN
          IF( vo_field(ii,jj-1) .NE. UNDEF)THEN
             vm=vo_field(ii,jj-1)
             lambdam=lambda(ii,jj-1)
          ELSE
             vm=vo_field(ii,jj)
             lambdam=0.0d0
          ENDIF
        ELSE
          vm=vo_field(ii,jj)
          lambdam=0.0d0
        ENDIF

        v_field(ii,jj) =0.25d0*(2.0d0*vo_field(ii,jj)+vp+vm)
+0.5d0*(lambdap-lambdam)/(dx*2.0d0)

     ENDIF
  ENDDO
ENDDO

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

!Weighted correlation, A and B are two variables and W are the weights
!that 
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

  COV = SUM( W*(A-MEANA)*(B-MEANB) )/SUM(W) !SUM( W*(A-MEANA)*(B-MEANB)
  C)  !Weigthed covariance

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

!---------------------------------------------------------------------------



