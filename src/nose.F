!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/nose [1.0] *
!!
!!   NAME
!!     nose
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM 20-Feb-2001: Now npt_ifo is allocated to zero when not used
!!     CJM 11-apr-2001: adding routines to thermostat ao_type 
!!
!!   SOURCE
!******************************************************************************

MODULE nose
 
  USE ao_types, ONLY : ao_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE molecule_types, ONLY : particle_node_type
  USE message_passing, ONLY : mp_bcast, mp_sum
  USE termination, ONLY : stop_program
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : gasdev

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lnhc, initialize_npt, initialize_nhc
  PUBLIC :: lnhc_parameters_type, nhc_info_type, npt_info_type
  PUBLIC :: yoshida_coef
  PUBLIC :: particle_kinetic_type, particle_scale_type

  TYPE npt_info_type
     REAL ( dbl ) :: eps
     REAL ( dbl ) :: v
     REAL ( dbl ) :: f
     REAL ( dbl ) :: mass
  END TYPE npt_info_type

  TYPE nhc_info_type
     REAL ( dbl ) :: eta
     REAL ( dbl ) :: v
     REAL ( dbl ) :: f
     REAL ( dbl ) :: nkt
     REAL ( dbl ) :: mass
  END TYPE nhc_info_type

  TYPE particle_kinetic_type
     REAL ( dbl ), POINTER :: point
  END TYPE particle_kinetic_type

  TYPE particle_scale_type
     REAL ( dbl ), POINTER :: point
  END TYPE particle_scale_type

  TYPE lnhc_parameters_type
     CHARACTER ( LEN = 10 ) :: dis_type
     INTEGER :: nyosh, nc, nhc_len, num_nhc
     REAL ( dbl ), POINTER :: dt_yosh ( : )
     TYPE (nhc_info_type ), POINTER :: nvt ( :, : )
     TYPE (particle_kinetic_type ), POINTER :: p_kin ( :, : )
     TYPE (particle_scale_type ), POINTER :: p_scale ( :, : )
     REAL ( dbl ), POINTER :: v_scale ( : )
     REAL ( dbl ), POINTER :: s_kin ( : )
  END TYPE lnhc_parameters_type

  INTERFACE lnhc
     MODULE PROCEDURE lnhc_particles, lnhc_barostat, lnhc_ao
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS
!******************************************************************************
!!****** nose/lnhc_ao [1.0] *
!!
!!   NAME
!!     lnhc_ao
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!   CREATION DATE
!!     11-Apr-2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE lnhc_ao ( nhc, ao, group )
  IMPLICIT NONE
! Arguments
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc
  TYPE ( ao_type ), INTENT ( INOUT ) :: ao
  INTEGER, INTENT ( IN ) :: group

! Locals
  INTEGER :: i, j, nnodes, inc, iyosh, n, inhc, handle, nx1, nx2, ncoef
  REAL ( dbl ) :: scale, flops

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0._dbl

  nx1 = SIZE ( nhc % nvt ( :, 1 ) )
  nx2 = SIZE ( nhc % nvt ( 1, : ) )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dbl
  nhc % s_kin = 0.0_dbl
  ncoef  = SIZE ( ao % cr)
  DO i = 1, ncoef 
    nhc % p_kin(1,i) % point  = nhc % p_kin(1,i) % point + &
      ao % mass ( i ) * ao % cr ( i )  * ao % cr ( i ) 
  END DO

! force on the first bead in every thermostat chain (f1=sum(m*v**2)-ncoef*KbT)
  DO n = 1, nhc % num_nhc
    nhc % nvt ( 1, n ) % f = ( nhc % s_kin ( n ) - nhc % nvt ( 1, n ) % nkt )/ &
                               nhc % nvt ( 1, n ) % mass
  END DO

! perform multiple time stepping using Yoshida
NCLOOP: DO inc = 1, nhc % nc
YOSH:     DO iyosh = 1, nhc % nyosh
!
! update velocity on the last thermostat in the chain    ! O1
       nhc % nvt ( nhc % nhc_len, : ) % v = nhc % nvt (nhc % nhc_len, : ) % v + &
         nhc % nvt (nhc % nhc_len, : ) % f * 0.25_dbl * nhc % dt_yosh ( iyosh )
!
! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
      DO n = 1, nhc % num_nhc
        DO inhc = nhc % nhc_len - 1, 1, -1
          scale = exp (-0.125_dbl * nhc % nvt ( inhc+1, n )% v * &
                 nhc % dt_yosh ( iyosh ) )
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v * scale     ! scale
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v + &
            nhc % nvt ( inhc, n ) % f * 0.25_dbl * nhc % dt_yosh ( iyosh ) ! shift
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc , n ) % v * scale    ! scale
        END DO
      END DO

! the core of the operator ----- START------
! update nhc positions
      nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
        0.5_dbl * nhc % nvt ( :, : ) % v * nhc % dt_yosh ( iyosh )
! now accumulate the scale factor for particle velocities
      DO n = 1, nhc % num_nhc
        nhc % v_scale ( n ) = nhc % v_scale ( n ) * &
        exp ( -0.5_dbl * nhc % dt_yosh ( &
             iyosh ) * nhc % nvt ( 1, n ) % v )
      END DO
! the core of the operator ------ END ------

! update the force on first thermostat again (since coefficient velocities 
! have changed)
      DO n = 1, nhc % num_nhc
        nhc % nvt ( 1, n ) % f = ( nhc % s_kin ( n ) * nhc % v_scale ( n ) &
          * nhc % v_scale ( n ) - nhc % nvt ( 1, n ) % nkt ) &
          /nhc % nvt ( 1, n ) % mass
      END DO

! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
      DO inhc = 1, nhc % nhc_len - 1
        DO n = 1, nhc % num_nhc
          scale = exp( -0.125_dbl *nhc % nvt ( inhc+1, n ) % v * & 
            nhc % dt_yosh ( iyosh ) ) 
          nhc % nvt ( inhc, n ) % v = nhc % nvt( inhc, n ) % v * scale ! scale
          nhc % nvt ( inhc, n ) % v = nhc % nvt( inhc, n ) % v + &
            nhc % nvt ( inhc, n ) % f * 0.25_dbl * nhc % dt_yosh ( iyosh ) ! shift
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v * scale ! scale
        END DO
! updating the forces on all  the thermostats
        DO n = 1, nhc % num_nhc
          nhc % nvt ( inhc+1, n ) % f =  & 
	 (nhc % nvt ( inhc, n ) % mass * nhc % nvt ( inhc, n ) % v &
         *nhc % nvt ( inhc, n ) % v-nhc % nvt ( inhc + 1, n ) % nkt ) / &
          nhc % nvt ( inhc + 1, n ) % mass
        END DO
      END DO
! update velocity on last thermostat                                  ! O1
      nhc % nvt ( nhc % nhc_len, : ) % v = &
      nhc % nvt ( nhc % nhc_len, : ) % v + &
      nhc % nvt ( nhc % nhc_len, : ) % f * &
       0.25_dbl * nhc % dt_yosh ( iyosh )
    END DO YOSH
  END DO NCLOOP
! now scale the coefficients velocities
  DO i = 1, ncoef
    ao % cr ( i ) = ao % cr ( i ) * nhc % p_scale ( 1, i ) % point 
  END DO
!
  CALL timestop(flops,handle)
END SUBROUTINE lnhc_ao

!******************************************************************************
!!****** nose/lnhc_barostat [1.0] *
!!
!!   NAME
!!     lnhc_barostat
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!   CREATION DATE
!!     13-DEC-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE lnhc_barostat ( nhc, npt, group )

  IMPLICIT NONE

! Arguments
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc
  TYPE ( npt_info_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: npt
  INTEGER, INTENT ( IN ) :: group

! Locals
  INTEGER :: i, j, nnodes, inc, iyosh, n, inhc, handle, nx1, nx2, ncoef
  REAL ( dbl ) :: scale, flops

!------------------------------------------------------------------------------

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dbl

  nx1 = SIZE ( nhc % nvt, 1 )
  nx2 = SIZE ( nhc % nvt, 2 )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dbl
  nhc % s_kin = 0.0_dbl
  ncoef = 0
  DO i = 1, SIZE ( npt, 1 )
    DO j = 1, SIZE ( npt, 2 )
      ncoef = ncoef + 1
      nhc % p_kin( 1, ncoef ) % point = nhc % p_kin ( 1, ncoef ) % point &
          + npt ( i, j ) % mass * npt ( i, j ) % v * npt ( i, j ) % v
    END DO
  END DO
#if defined(__parallel)
! for replicated thermostats and distributed molecules, we have to
! sum s_kin (pointed to by p_kin%point) over all processors
  IF ( nhc % dis_type == 'INTER_REP' ) THEN
     CALL mp_sum ( nhc % s_kin, group )
  END IF
#endif
  flops = flops + REAL ( ncoef, dbl )

! force on the first bead in every thermostat chain
  DO n = 1, nhc % num_nhc
     nhc % nvt(1,n) % f = (nhc % s_kin(n)-nhc % nvt(1,n) % nkt)/ &
          nhc % nvt(1,n) % mass
  END DO
  flops = flops + REAL ( nhc % num_nhc * 2, dbl )

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhc % nc
     YOSH: DO iyosh = 1, nhc % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dbl*nhc % dt_yosh(iyosh)
        flops = flops + REAL ( 3 * nx2, dbl )

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhc % num_nhc
           DO inhc = nhc % nhc_len - 1, 1, -1
              scale = exp(-0.125_dbl*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dbl*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
        END DO
        flops = flops &
             + REAL ( nhc % num_nhc * (nhc % nhc_len - 1 ) * 8, dbl )

! the core of the operator ----- START------
! update nhc positions
        nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
             0.5_dbl*nhc % nvt ( :, : ) % v*nhc % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhc % num_nhc
           nhc % v_scale(n) = nhc % v_scale(n)*exp(-0.5_dbl*nhc % dt_yosh( &
                iyosh)*nhc % nvt(1,n) % v)
        END DO
        flops = flops + REAL(nhc % num_nhc*4)
! the core of the operator ------ END ------

! update the force on first thermostat again (since particle velocities
! have changed)
        DO n = 1, nhc % num_nhc
           nhc % nvt(1,n) % f = (nhc % s_kin(n)*nhc % v_scale(n)*nhc % v_scale(n &
                )-nhc % nvt(1,n) % nkt)/nhc % nvt(1,n) % mass
        END DO
        flops = flops + REAL(nhc % num_nhc*4)

! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
        DO inhc = 1, nhc % nhc_len - 1
           DO n = 1, nhc % num_nhc
              scale = exp(-0.125_dbl*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dbl*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
! updating the forces on all the thermostats
           DO n = 1, nhc % num_nhc
              nhc % nvt(inhc+1,n) % f = (nhc % nvt(inhc,n) % mass*nhc % nvt(inhc,n) % v &
                   *nhc % nvt(inhc,n) % v-nhc % nvt(inhc+1,n) % nkt)/ &
                   nhc % nvt(inhc+1,n) % mass
           END DO
        END DO
        flops = flops + REAL(nhc % nhc_len-1)*REAL(nhc % num_nhc*12)
! update velocity on last thermostat                                    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dbl*nhc % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

! now scale the particle velocities
  ncoef = 0
  DO i = 1,  SIZE ( npt, 1 )
    DO j = 1,  SIZE ( npt, 2 )
      ncoef = ncoef + 1
      npt ( i, j ) % v = npt ( i, j ) % v * nhc % p_scale( 1, ncoef ) % point
    END DO
  END DO
  flops = flops + REAL ( ncoef, dbl )

  flops = flops * 1.E-6_dbl
  CALL timestop(flops,handle)

END SUBROUTINE lnhc_barostat

!******************************************************************************
!!****** nose/lnhc_particles [1.0] *
!!
!!   NAME
!!     lnhc_particles
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE lnhc_particles ( nhc, pnode, group )

  IMPLICIT NONE

! Arguments
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  INTEGER, INTENT ( IN ) :: group

! Locals
  INTEGER :: i, nnodes, inc, iyosh, n, inhc, handle, nx1, nx2
  REAL ( dbl ) :: scale, flops

!------------------------------------------------------------------------------

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dbl

  nx1 = SIZE ( nhc % nvt, 1 )
  nx2 = SIZE ( nhc % nvt, 2 )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dbl
  nhc % s_kin = 0.0_dbl
  nnodes = SIZE ( pnode )
  DO i = 1, nnodes
     nhc % p_kin(1,i) % point = nhc % p_kin(1,i) % point &
          + pnode(i) % p % prop % mass*pnode(i) % p % v(1)*pnode(i) % p % v(1)
     nhc % p_kin(2,i) % point = nhc % p_kin(2,i) % point &
          + pnode(i) % p % prop % mass*pnode(i) % p % v(2)*pnode(i) % p % v(2)
     nhc % p_kin(3,i) % point = nhc % p_kin(3,i) % point &
          + pnode(i) % p % prop % mass*pnode(i) % p % v(3)*pnode(i) % p % v(3)
  END DO
#if defined(__parallel)
! for replicated thermostats and distributed molecules, we have to
! sum s_kin (pointed to by p_kin%point) over all processors
  IF ( nhc % dis_type == 'INTER_REP' ) THEN
     CALL mp_sum ( nhc % s_kin, group )
  END IF
#endif
  flops = flops + REAL ( nnodes * 9, dbl )

! force on the first bead in every thermostat chain
  DO n = 1, nhc % num_nhc
     nhc % nvt(1,n) % f = (nhc % s_kin(n)-nhc % nvt(1,n) % nkt)/ &
          nhc % nvt(1,n) % mass
  END DO
  flops = flops + REAL ( nhc % num_nhc * 2, dbl )

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhc % nc
     YOSH: DO iyosh = 1, nhc % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dbl*nhc % dt_yosh(iyosh)
        flops = flops + REAL ( 3 * nx2, dbl )

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhc % num_nhc
           DO inhc = nhc % nhc_len - 1, 1, -1
              scale = exp(-0.125_dbl*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dbl*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
        END DO
        flops = flops &
             + REAL ( nhc % num_nhc * (nhc % nhc_len - 1 ) * 8, dbl )

! the core of the operator ----- START------
! update nhc positions
        nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
             0.5_dbl*nhc % nvt ( :, : ) % v*nhc % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhc % num_nhc
           nhc % v_scale(n) = nhc % v_scale(n)*exp(-0.5_dbl*nhc % dt_yosh( &
                iyosh)*nhc % nvt(1,n) % v)
        END DO
        flops = flops + REAL(nhc % num_nhc*4)
! the core of the operator ------ END ------

! update the force on first thermostat again (since particle velocities
! have changed)
        DO n = 1, nhc % num_nhc
           nhc % nvt(1,n) % f = (nhc % s_kin(n)*nhc % v_scale(n)*nhc % v_scale(n &
                )-nhc % nvt(1,n) % nkt)/nhc % nvt(1,n) % mass
        END DO
        flops = flops + REAL(nhc % num_nhc*4)

! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
        DO inhc = 1, nhc % nhc_len - 1
           DO n = 1, nhc % num_nhc
              scale = exp(-0.125_dbl*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dbl*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
! updating the forces on all the thermostats
           DO n = 1, nhc % num_nhc
              nhc % nvt(inhc+1,n) % f = (nhc % nvt(inhc,n) % mass*nhc % nvt(inhc,n) % v &
                   *nhc % nvt(inhc,n) % v-nhc % nvt(inhc+1,n) % nkt)/ &
                   nhc % nvt(inhc+1,n) % mass
           END DO
        END DO
        flops = flops + REAL(nhc % nhc_len-1)*REAL(nhc % num_nhc*12)
! update velocity on last thermostat                                    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dbl*nhc % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

! now scale the particle velocities
  DO i = 1, nnodes
     pnode(i) % p % v(1) = pnode(i) % p % v(1)*nhc % p_scale(1,i) % point
     pnode(i) % p % v(2) = pnode(i) % p % v(2)*nhc % p_scale(2,i) % point
     pnode(i) % p % v(3) = pnode(i) % p % v(3)*nhc % p_scale(3,i) % point
  END DO
  flops = flops + REAL ( nnodes*3, dbl )

  flops = flops * 1.E-6_dbl
  CALL timestop(flops,handle)

END SUBROUTINE lnhc_particles

!******************************************************************************
!!****** nose/initialize_nhc [1.0] *
!!
!!   NAME
!!     initialize_nhc
!!
!!   FUNCTION
!!     Initializes the NHC velocities to the Maxwellian distribution
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_nhc(nhc,tau_nhc,temp_ext, globenv )
  IMPLICIT NONE
  TYPE (lnhc_parameters_type ), INTENT ( INOUT ) :: nhc
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  REAL ( dbl ), INTENT ( IN ) :: tau_nhc, temp_ext
  INTEGER :: i, j, number
  REAL ( dbl ) :: akin, temp, v

! first initializing the mass of the nhc variables
  nhc % nvt ( :, : ) % mass = nhc % nvt ( :, : ) % nkt * tau_nhc * tau_nhc
! initializing velocities
  DO i = 1, nhc % num_nhc
     DO j = 1, nhc % nhc_len
        v = gasdev(globenv % idum)
#if defined(__parallel)
        IF (nhc % dis_type=='INTER_REP') CALL mp_bcast(v,globenv % source, &
             globenv % group)
#endif
        nhc % nvt(j,i) % v = v
     END DO
  END DO
  akin = 0.0_dbl
  DO i = 1, nhc % num_nhc
     DO j = 1, nhc % nhc_len
        akin = akin + 0.5_dbl * ( nhc % nvt ( j , i ) % mass * &
                                  nhc % nvt ( j , i ) % v * &
                                  nhc % nvt ( j , i ) % v )
     END DO
  END DO
  number = nhc % num_nhc
#if defined(__parallel)
  IF (nhc % dis_type=='DIS_REP') CALL mp_sum(akin,globenv % group)
  IF (nhc % dis_type=='DIS_REP') CALL mp_sum(number,globenv % group)
#endif

! scale velocities to get the correct initial temperature
  temp = 2.0_dbl*akin/REAL(number)
  temp = sqrt(temp_ext/temp)
  DO i = 1, nhc % num_nhc
     DO j = 1, nhc % nhc_len
        nhc % nvt(j,i) % v = temp * nhc % nvt(j,i) % v
        nhc % nvt(j,i) % eta = 0.0_dbl
     END DO
  END DO

! initializing all of the forces on the thermostats
  DO i = 1, nhc % num_nhc
     DO j = 2, nhc % nhc_len
        nhc % nvt(j,i) % f = nhc % nvt(j-1,i) % mass*nhc % nvt(j-1,i) % v* &
             nhc % nvt(j-1,i) % v - nhc % nvt(j,i) % nkt
        nhc % nvt(j,i) % f = nhc % nvt(j,i) % f/nhc % nvt(j,i) % mass
     END DO
  END DO

END SUBROUTINE initialize_nhc

!!*****
!******************************************************************************
!!****** nose/initialize_npt [1.0] *
!!
!!   NAME
!!     initialize_npt
!!
!!   FUNCTION
!!     Initializes the NHC velocities to the Maxwellian distribution
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_npt ( npt, tau_cell, temp_ext, deth, nfree, ensemble, &
     globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( npt_info_type ), DIMENSION ( :, : ), INTENT ( OUT ) :: npt
  REAL ( dbl ), INTENT ( IN ) :: tau_cell, temp_ext, deth
  CHARACTER ( LEN = 20 ), INTENT ( IN ) :: ensemble
  INTEGER, INTENT ( IN ) :: nfree
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  INTEGER :: i, j, number, isos
  REAL ( dbl ) :: akin, temp, v

!------------------------------------------------------------------------------
! first initializing the mass of the nhc variables
  SELECT CASE ( ensemble )
  CASE ( 'NPT_I')
     npt ( :, : ) % mass = REAL ( nfree + 3, dbl ) &
          * temp_ext * tau_cell ** 2

  CASE ( 'NPT_F')
     npt ( :, : ) % mass = REAL ( nfree + 3, dbl ) &
          * temp_ext * tau_cell ** 2 / 3.0_dbl
  END SELECT

! initializing velocities
  DO i = 1, size ( npt,1)
     DO j = i, size ( npt,2)
        v = gasdev ( globenv % idum)
#if defined(__parallel)
        CALL mp_bcast ( v, globenv % source, globenv % group )
#endif
! symmetrizing the initial barostat velocities to ensure
! no rotation of the cell under NPT_F
        npt(j,i) % v = v
        npt(i,j) % v = v
     END DO
  END DO

  akin = 0.0_dbl
  DO i = 1, SIZE( npt,1 )
     DO j = 1, SIZE( npt,2 )
        akin = akin + 0.5_dbl*( npt(j,i) % mass* npt(j,i) % v* npt &
             (j,i) % v)
     END DO
  END DO

  number = SIZE ( npt, 1 ) * SIZE ( npt, 2 )

! scale velocities to get the correct initial temperature
  IF ( number /= 0 ) THEN
    temp = 2.0_dbl * akin / REAL ( number, dbl )
    temp = sqrt(temp_ext/temp)
  ENDIF
  DO i = 1, size(npt,1)
     DO j = i, size(npt,2)
        npt(j,i) % v = temp * npt(j,i) % v
        npt(i,j) % v = npt(j,i) % v
! dbg ISOTROPIC LIMIT
!         npt(j,i) % v = 0._dbl
!       write (*, *) 'INTIAL v_eps', npt(j,i) % v
! dbg ISOTROPIC LIMIT
     END DO
  END DO

! initialize the position
  IF ( ensemble == 'NPT_I' ) npt ( :, : ) % eps = LOG ( deth ) / 3.0_dbl

END SUBROUTINE initialize_npt

!!*****
!******************************************************************************
!!****** nose/yoshida_coef [1.0] *
!!
!!   NAME
!!     yoshida_coef
!!
!!   FUNCTION
!!     This lists the coefficients for the Yoshida method (higher
!!     order integrator used in NVT)
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE yoshida_coef ( nhc, dt )

  IMPLICIT NONE

! Arguments
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc
  REAL ( dbl ), INTENT ( IN ) :: dt

! Locals
  REAL ( dbl ), DIMENSION (nhc % nyosh) :: yosh_wt

!------------------------------------------------------------------------------

  SELECT CASE (nhc % nyosh)
  CASE DEFAULT
     CALL stop_program ( 'yoshida_coef', 'Value not available' )
  CASE (1)
     yosh_wt(1) = 1.0_dbl
  CASE (3)
     yosh_wt(1) = 1.0_dbl/(2.0_dbl-(2.0_dbl)**(1.0_dbl/3.0_dbl))
     yosh_wt(2) = 1.0_dbl - 2.0_dbl*yosh_wt(1)
     yosh_wt(3) = yosh_wt(1)
  CASE (5)
     yosh_wt(1) = 1.0_dbl/(4.0_dbl-(4.0_dbl)**(1.0_dbl/3.0_dbl))
     yosh_wt(2) = yosh_wt(1)
     yosh_wt(4) = yosh_wt(1)
     yosh_wt(5) = yosh_wt(1)
     yosh_wt(3) = 1.0_dbl - 4.0_dbl*yosh_wt(1)
  CASE (7)
     yosh_wt(1) = .78451361047756_dbl
     yosh_wt(2) = .235573213359357_dbl
     yosh_wt(3) = -1.17767998417887_dbl
     yosh_wt(4) = 1.0_dbl - 2.0_dbl*(yosh_wt(1)+yosh_wt(2)+yosh_wt(3))
     yosh_wt(5) = yosh_wt(3)
     yosh_wt(6) = yosh_wt(2)
     yosh_wt(7) = yosh_wt(1)
  CASE (9)
     yosh_wt(1) = 0.192_dbl
     yosh_wt(2) = 0.554910818409783619692725006662999_dbl
     yosh_wt(3) = 0.124659619941888644216504240951585_dbl
     yosh_wt(4) = -0.843182063596933505315033808282941_dbl
     yosh_wt(5) = 1.0_dbl - 2.0_dbl*(yosh_wt(1)+yosh_wt(2)+&
                  yosh_wt(3)+yosh_wt(4))
     yosh_wt(6) = yosh_wt(4)
     yosh_wt(7) = yosh_wt(3)
     yosh_wt(8) = yosh_wt(2)
     yosh_wt(9) = yosh_wt(1)
  CASE (15)
     yosh_wt(1) = 0.102799849391985_dbl
     yosh_wt(2) = -0.196061023297549e1_dbl
     yosh_wt(3) = 0.193813913762276e1_dbl
     yosh_wt(4) = -0.158240635368243_dbl
     yosh_wt(5) = -0.144485223686048e1_dbl
     yosh_wt(6) = 0.253693336566229_dbl
     yosh_wt(7) = 0.914844246229740_dbl
     yosh_wt(8) = 1.0_dbl - 2.0_dbl*(yosh_wt(1)+yosh_wt(2)+&
                 yosh_wt(3)+yosh_wt(4)+yosh_wt(5)+yosh_wt(6)+yosh_wt(7))
     yosh_wt(9) = yosh_wt(7)
     yosh_wt(10) = yosh_wt(6)
     yosh_wt(11) = yosh_wt(5)
     yosh_wt(12) = yosh_wt(4)
     yosh_wt(13) = yosh_wt(3)
     yosh_wt(14) = yosh_wt(2)
     yosh_wt(15) = yosh_wt(1)
  END SELECT
  nhc % dt_yosh = dt * yosh_wt / REAL ( nhc % nc, dbl )

END SUBROUTINE yoshida_coef

!!*****
!******************************************************************************

END MODULE nose

!******************************************************************************
