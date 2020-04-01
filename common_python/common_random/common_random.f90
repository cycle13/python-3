MODULE common_random

!=======================================================================
!
! [PURPOSE:] General constants and procedures
!
! [ATTENTION:] This module calls 'SFMT.f90'
!
! [HISTORY:]
!   07/20/2004 Takemasa MIYOSHI  created
!   01/23/2009 Takemasa MIYOSHI  modified for SFMT
!
!=======================================================================
  IMPLICIT NONE
  PUBLIC
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)

  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.0d0
  REAL(r_size),PARAMETER :: cp=7.0d0 / 2.0d0 * rd
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0


!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
! RAND (random number with uniform distribution)
!-----------------------------------------------------------------------
SUBROUTINE com_rand(ndim,var)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(OUT) :: var(1:ndim)
  REAL(r_dble) :: genrand_res53
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    CALL init_gen_rand(iseed)
    first=.false.
  END IF

  DO i=1,ndim
    var(i) = genrand_res53()
  END DO

  RETURN
END SUBROUTINE com_rand
!-----------------------------------------------------------------------
! RANDN (random number with normal distribution)
!-----------------------------------------------------------------------
SUBROUTINE com_randn(ndim,var)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(OUT) :: var(1:ndim)
  REAL(r_size) :: rnd(2)
  REAL(r_dble) :: genrand_res53
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    CALL init_gen_rand(iseed)
    first=.false.
  END IF

  IF( MOD(ndim,2)==0 ) THEN
    DO i=1,ndim/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
  ELSE
    DO i=1,(ndim-1)/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
    rnd(1) = genrand_res53()
    rnd(2) = genrand_res53()
    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
  END IF

  RETURN
END SUBROUTINE com_randn

!-----------------------------------------------------------------------
! RANDN (random number with normal distribution) with auxiliar input to compute
! seed
! This routine is to generate different random sequences in parallel
! environments
! Aux can be the PID or a function of the PID.
!-----------------------------------------------------------------------
SUBROUTINE com_randn2(ndim,var,aux)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(OUT) :: var(1:ndim)
  REAL(r_size) :: rnd(2)
  INTEGER :: aux
  REAL(r_dble) :: genrand_res53
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000 + aux*1000
    CALL init_gen_rand(iseed)
    first=.false.
  END IF

  IF( MOD(ndim,2)==0 ) THEN
    DO i=1,ndim/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
  ELSE
    DO i=1,(ndim-1)/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
    rnd(1) = genrand_res53()
    rnd(2) = genrand_res53()
    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
  END IF

  RETURN
END SUBROUTINE com_randn2


END MODULE common_random
