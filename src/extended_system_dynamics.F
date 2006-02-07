!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/extended_system_dynamics [1.0] *
!!
!!   NAME
!!     extended_system_dynamics
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM 20-Feb-2001: Now npt_ifo is allocated to zero when not used
!!     CJM 11-apr-2001: adding routines to thermostat ao_type 
!!     CJM 02-Aug-2003: renamed
!!
!!   SOURCE
!******************************************************************************

MODULE extended_system_dynamics
 
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE input_constants,                 ONLY: dyn_coeff_mass_global
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lnhc

  INTERFACE lnhc
     MODULE PROCEDURE lnhc_particles, lnhc_particle_set, lnhc_barostat, lnhc_coefs, lnhc_pw
  END INTERFACE

!!*****
!******************************************************************************
CONTAINS
!******************************************************************************
!!****** extended_system_dynamics/lnhc_coefs [1.0] *
!!
!!   NAME
!!     lnhc_coefs
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
!!     05-02-04 adapted to new coefs type
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE lnhc_coefs ( nhc, dyn_coeff_set, group )
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: handle, i, icoef, iel, ikind, &
                                                inc, inhc, iyosh, n, nkind, &
                                                nx1, nx2
    LOGICAL                                  :: gmass
    REAL(KIND=dp)                            :: flops, scale
    TYPE(dyn_coeff_type), POINTER            :: coeffs

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dp

  nx1 = SIZE ( nhc % nvt ,1 )
  nx2 = SIZE ( nhc % nvt ,2 )
  gmass=(dyn_coeff_set%mass_storage==dyn_coeff_mass_global)

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dp
  nhc % s_kin = 0.0_dp
  i  = 0               
  nkind = SIZE ( dyn_coeff_set%coeffs_of_kind)

  DO ikind=1,nkind
    IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
      coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
      DO iel = 1, coeffs%n_els
        DO icoef = 1, coeffs%ncoef_atom
          i = i + 1
          IF(gmass) THEN
            nhc % p_kin(1,i) % point  = nhc % p_kin(1,i) % point + &
            dyn_coeff_set%global_mass*coeffs%vel(iel,icoef)*coeffs%vel(iel,icoef)
          ELSE
            nhc % p_kin(1,i) % point  = nhc % p_kin(1,i) % point + &
            coeffs%masses(iel,icoef)*coeffs%vel(iel,icoef)*coeffs%vel(iel,icoef)
          END IF 
        END DO
      END DO
    END IF
  END DO
  IF ( nhc % dis_type == 'INTER_REP' )  CALL mp_sum ( nhc % s_kin, group )

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
         nhc % nvt (nhc % nhc_len, : ) % f * 0.25_dp * nhc % dt_yosh ( iyosh )
!
! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
      DO n = 1, nhc % num_nhc
        DO inhc = nhc % nhc_len - 1, 1, -1
          scale = EXP (-0.125_dp * nhc % nvt ( inhc+1, n )% v * &
                 nhc % dt_yosh ( iyosh ) )
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v * scale     ! scale
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v + &
            nhc % nvt ( inhc, n ) % f * 0.25_dp * nhc % dt_yosh ( iyosh ) ! shift
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc , n ) % v * scale    ! scale
        END DO
      END DO

! the core of the operator ----- START------
! update nhc positions
      nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
        0.5_dp * nhc % nvt ( :, : ) % v * nhc % dt_yosh ( iyosh )
! now accumulate the scale factor for particle velocities
      DO n = 1, nhc % num_nhc
        nhc % v_scale ( n ) = nhc % v_scale ( n ) * &
        EXP ( -0.5_dp * nhc % dt_yosh ( &
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
          scale = EXP( -0.125_dp *nhc % nvt ( inhc+1, n ) % v * & 
            nhc % dt_yosh ( iyosh ) ) 
          nhc % nvt ( inhc, n ) % v = nhc % nvt( inhc, n ) % v * scale ! scale
          nhc % nvt ( inhc, n ) % v = nhc % nvt( inhc, n ) % v + &
            nhc % nvt ( inhc, n ) % f * 0.25_dp * nhc % dt_yosh ( iyosh ) ! shift
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
       0.25_dp * nhc % dt_yosh ( iyosh )
    END DO YOSH
  END DO NCLOOP

  i=0
! now scale the coefficients velocities
  DO ikind=1,nkind
    IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
      coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
      DO iel = 1, coeffs%n_els
        DO icoef = 1, coeffs%ncoef_atom
          i = i + 1
          coeffs%vel(iel,icoef)=coeffs%vel(iel,icoef)*nhc%p_scale(1,i)%point
        END DO
      END DO
    END IF
  END DO
!
  CALL timestop(flops,handle)
END SUBROUTINE lnhc_coefs

!******************************************************************************
!!****** extended_system_dynamics/lnhc_pw [1.0] *
!!
!!   NAME
!!     lnhc_pw
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
SUBROUTINE lnhc_pw ( nhc, pw, group )
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(pw_type), INTENT(INOUT)             :: pw
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: handle, i, inc, inhc, iyosh, &
                                                n, ncoef, nx1, nx2
    REAL(KIND=dp)                            :: flops, scale

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dp

  nx1 = SIZE ( nhc % nvt ,1 )
  nx2 = SIZE ( nhc % nvt ,2 )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dp
  nhc % s_kin = 0.0_dp
  ncoef  = SIZE ( pw % cc )

  DO i = 1, ncoef
    nhc % p_kin ( 1, i ) % point  = nhc % p_kin ( 1, i ) % point + &
      pw % mass_cc ( i ) * REAL ( pw % cc ( i ),KIND=dp)  * REAL ( pw % cc ( i ),KIND=dp) 

    nhc % p_kin ( 2 , i ) % point  = nhc % p_kin ( 2, i ) % point + &
      pw % mass_cc ( i ) * AIMAG ( pw % cc ( i ) )  * AIMAG ( pw % cc ( i ) ) 
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
         nhc % nvt (nhc % nhc_len, : ) % f * 0.25_dp * nhc % dt_yosh ( iyosh )
!
! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
      DO n = 1, nhc % num_nhc
        DO inhc = nhc % nhc_len - 1, 1, -1
          scale = EXP (-0.125_dp * nhc % nvt ( inhc+1, n )% v * &
                 nhc % dt_yosh ( iyosh ) )
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v * scale     ! scale
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc, n ) % v + &
            nhc % nvt ( inhc, n ) % f * 0.25_dp * nhc % dt_yosh ( iyosh ) ! shift
          nhc % nvt ( inhc, n ) % v = nhc % nvt ( inhc , n ) % v * scale    ! scale
        END DO
      END DO

! the core of the operator ----- START------
! update nhc positions
      nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
        0.5_dp * nhc % nvt ( :, : ) % v * nhc % dt_yosh ( iyosh )
! now accumulate the scale factor for particle velocities
      DO n = 1, nhc % num_nhc
        nhc % v_scale ( n ) = nhc % v_scale ( n ) * &
        EXP ( -0.5_dp * nhc % dt_yosh ( &
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
          scale = EXP( -0.125_dp *nhc % nvt ( inhc+1, n ) % v * & 
            nhc % dt_yosh ( iyosh ) ) 
          nhc % nvt ( inhc, n ) % v = nhc % nvt( inhc, n ) % v * scale ! scale
          nhc % nvt ( inhc, n ) % v = nhc % nvt( inhc, n ) % v + &
            nhc % nvt ( inhc, n ) % f * 0.25_dp * nhc % dt_yosh ( iyosh ) ! shift
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
       0.25_dp * nhc % dt_yosh ( iyosh )
    END DO YOSH
  END DO NCLOOP
! now scale the coefficients velocities
  DO i = 1, ncoef
     pw % cc ( i ) =  CMPLX ( REAL ( pw % cc ( i ),KIND=dp) * nhc % p_scale ( 1, i ) % point, &
                      AIMAG ( pw % cc ( i ) ) * nhc % p_scale ( 2, i ) % point,KIND=dp )
  END DO
!
  CALL timestop(flops,handle)
END SUBROUTINE lnhc_pw

!******************************************************************************
!!****** extended_system_dynamics/lnhc_barostat [1.0] *
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


    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: handle, i, inc, inhc, iyosh, &
                                                j, n, ncoef, nx1, nx2
    REAL(KIND=dp)                            :: flops, scale

!------------------------------------------------------------------------------

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dp

  nx1 = SIZE ( nhc % nvt, 1 )
  nx2 = SIZE ( nhc % nvt, 2 )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dp
  nhc % s_kin = 0.0_dp
  ncoef = 0
  DO i = 1, SIZE ( npt, 1 )
    DO j = 1, SIZE ( npt, 2 )
      ncoef = ncoef + 1
      nhc % p_kin( 1, ncoef ) % point = nhc % p_kin ( 1, ncoef ) % point &
          + npt ( i, j ) % mass * npt ( i, j ) % v * npt ( i, j ) % v
    END DO
  END DO

  IF ( nhc % dis_type == 'INTER_REP' ) CALL mp_sum ( nhc % s_kin, group )

  flops = flops + REAL ( ncoef,KIND=dp)

! force on the first bead in every thermostat chain
  DO n = 1, nhc % num_nhc
     nhc % nvt(1,n) % f = (nhc % s_kin(n)-nhc % nvt(1,n) % nkt)/ &
          nhc % nvt(1,n) % mass
  END DO
  flops = flops + REAL ( nhc % num_nhc * 2,KIND=dp)

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhc % nc
     YOSH: DO iyosh = 1, nhc % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dp*nhc % dt_yosh(iyosh)
        flops = flops + REAL ( 3 * nx2,KIND=dp)

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhc % num_nhc
           DO inhc = nhc % nhc_len - 1, 1, -1
              scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
        END DO
        flops = flops &
             + REAL ( nhc % num_nhc * (nhc % nhc_len - 1 ) * 8,KIND=dp)

! the core of the operator ----- START------
! update nhc positions
        nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
             0.5_dp*nhc % nvt ( :, : ) % v*nhc % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhc % num_nhc
           nhc % v_scale(n) = nhc % v_scale(n)*EXP(-0.5_dp*nhc % dt_yosh( &
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
              scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh) ! shift
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
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dp*nhc % dt_yosh(iyosh)
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
  flops = flops + REAL ( ncoef,KIND=dp)

  flops = flops * 1.E-6_dp
  CALL timestop(flops,handle)

END SUBROUTINE lnhc_barostat

!******************************************************************************
!!****** extended_system_dynamics/lnhc_particle_set [1.0] *
!!
!!   NAME
!!     lnhc_particle_set
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
SUBROUTINE lnhc_particle_set ( nhc, molecule_kind_set, molecule_set,  &
                               atomic_kind_set, particle_set, local_molecules, &
                               group )


    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: group

    INTEGER :: first_atom, handle, i, ii, ikind, imol, imol_local, inc, inhc, &
      ipart, iyosh, last_atom, n, nmol_local, npart_local, nx1, nx2
    REAL(KIND=dp)                            :: flops, mass, scale
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule

!    TYPE(distribution_1d_type), POINTER      :: particles
!------------------------------------------------------------------------------

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dp

  nx1 = SIZE ( nhc % nvt, 1 )
  nx2 = SIZE ( nhc % nvt, 2 )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dp
  nhc % s_kin = 0.0_dp
  ii = 0
  DO ikind = 1, SIZE ( molecule_kind_set )
    nmol_local = local_molecules % n_el ( ikind )
    DO imol_local = 1, nmol_local
      imol = local_molecules % list ( ikind ) % array ( imol_local )
      molecule => molecule_set ( imol )
      CALL get_molecule ( molecule, first_atom = first_atom, &
                          last_atom = last_atom ) 
      DO ipart = first_atom, last_atom
        ii = ii + 1
        atomic_kind => particle_set ( ipart ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
        nhc % p_kin ( 1, ii ) % point = nhc % p_kin ( 1, ii ) % point &
            + mass * particle_set ( ipart ) % v ( 1 ) *            &
                     particle_set ( ipart ) % v ( 1 )
        nhc % p_kin ( 2, ii ) % point = nhc % p_kin ( 2, ii ) % point &
            + mass * particle_set ( ipart ) % v ( 2 ) *            &
                     particle_set ( ipart ) % v ( 2 )
        nhc % p_kin ( 3, ii ) % point = nhc % p_kin ( 3, ii ) % point &
            + mass * particle_set ( ipart ) % v ( 3 ) *            &
                     particle_set ( ipart ) % v ( 3 )
      END DO
    END DO
  END DO
  npart_local=ii

  IF ( nhc % dis_type == 'INTER_REP' )  CALL mp_sum ( nhc % s_kin, group )

  flops = flops + REAL ( npart_local * 9,KIND=dp)

! force on the first bead in every thermostat chain
  DO n = 1, nhc % num_nhc
     nhc % nvt(1,n) % f = (nhc % s_kin(n)-nhc % nvt(1,n) % nkt)/ &
          nhc % nvt(1,n) % mass
  END DO
  flops = flops + REAL ( nhc % num_nhc * 2,KIND=dp)

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhc % nc
     YOSH: DO iyosh = 1, nhc % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dp*nhc % dt_yosh(iyosh)
        flops = flops + REAL ( 3 * nx2,KIND=dp)

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhc % num_nhc
           DO inhc = nhc % nhc_len - 1, 1, -1
              scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
        END DO
        flops = flops &
             + REAL ( nhc % num_nhc * (nhc % nhc_len - 1 ) * 8,KIND=dp)

! the core of the operator ----- START------
! update nhc positions
        nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
             0.5_dp*nhc % nvt ( :, : ) % v*nhc % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhc % num_nhc
           nhc % v_scale(n) = nhc % v_scale(n)*EXP(-0.5_dp*nhc % dt_yosh( &
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
              scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh) ! shift
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
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dp*nhc % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

! now scale the particle velocities
  ii = 0
  DO ikind = 1, SIZE ( molecule_kind_set )
    nmol_local = local_molecules % n_el ( ikind )
    DO imol_local = 1, nmol_local
      imol = local_molecules % list ( ikind ) % array ( imol_local )
      molecule => molecule_set ( imol )
      CALL get_molecule ( molecule, first_atom = first_atom, &
                          last_atom = last_atom ) 
      DO ipart = first_atom, last_atom
        ii = ii + 1
        particle_set ( ipart ) % v ( 1 ) =  particle_set ( ipart ) % v ( 1 ) * &
                                            nhc % p_scale ( 1, ii ) % point
        particle_set ( ipart ) % v ( 2 ) =  particle_set ( ipart ) % v ( 2 ) * &
                                            nhc % p_scale ( 2, ii ) % point
        particle_set ( ipart ) % v ( 3 ) =  particle_set ( ipart ) % v ( 3 ) * &
                                            nhc % p_scale ( 3, ii ) % point
      END DO
    END DO
  END DO
  flops = flops + REAL ( npart_local*3,KIND=dp)

  flops = flops * 1.E-6_dp
  CALL timestop(flops,handle)

END SUBROUTINE lnhc_particle_set
!******************************************************************************
!!****** extended_system_dynamics/lnhc_particles[1.0] *
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
SUBROUTINE lnhc_particles ( nhc, molecule_kind_set, molecule_set,  &
                               atomic_kind_set, particle_set, local_molecules, &
                               vel, group )


    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    INTEGER, INTENT(IN)                      :: group

    INTEGER :: first_atom, handle, i, ii, ikind, imol, imol_local, inc, inhc, &
      ipart, iyosh, last_atom, n, nmol_local, npart_local, nx1, nx2
    REAL(KIND=dp)                            :: flops, mass, scale
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dp

  nx1 = SIZE ( nhc % nvt, 1 )
  nx2 = SIZE ( nhc % nvt, 2 )

! get force on first thermostat for all the chains in the system.
  nhc % v_scale = 1.0_dp
  nhc % s_kin = 0.0_dp
  ii = 0
  DO ikind = 1, SIZE ( molecule_kind_set )
    nmol_local = local_molecules % n_el ( ikind )
    DO imol_local = 1, nmol_local
      imol = local_molecules % list ( ikind ) % array ( imol_local )
      molecule => molecule_set ( imol )
      CALL get_molecule ( molecule, first_atom = first_atom, &
                          last_atom = last_atom ) 
      DO ipart = first_atom, last_atom
        ii = ii + 1
        atomic_kind => particle_set ( ipart ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
        nhc % p_kin ( 1, ii ) % point = nhc % p_kin ( 1, ii ) % point &
            + mass * vel ( 1, ipart ) * vel ( 1, ipart )
        nhc % p_kin ( 2, ii ) % point = nhc % p_kin ( 2, ii ) % point &
            + mass * vel ( 2, ipart ) * vel ( 2, ipart )
        nhc % p_kin ( 3, ii ) % point = nhc % p_kin ( 3, ii ) % point &
            + mass * vel ( 3, ipart ) * vel ( 3, ipart )
      END DO
    END DO
  END DO
  npart_local=ii

  IF ( nhc % dis_type == 'INTER_REP' )  CALL mp_sum ( nhc % s_kin, group )

  flops = flops + REAL ( npart_local * 9,KIND=dp)

! force on the first bead in every thermostat chain
  DO n = 1, nhc % num_nhc
     nhc % nvt(1,n) % f = (nhc % s_kin(n)-nhc % nvt(1,n) % nkt)/ &
          nhc % nvt(1,n) % mass
  END DO
  flops = flops + REAL ( nhc % num_nhc * 2,KIND=dp)

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhc % nc
     YOSH: DO iyosh = 1, nhc % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dp*nhc % dt_yosh(iyosh)
        flops = flops + REAL ( 3 * nx2,KIND=dp)

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhc % num_nhc
           DO inhc = nhc % nhc_len - 1, 1, -1
              scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh) ! shift
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
           END DO
        END DO
        flops = flops &
             + REAL ( nhc % num_nhc * (nhc % nhc_len - 1 ) * 8,KIND=dp)

! the core of the operator ----- START------
! update nhc positions
        nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
             0.5_dp*nhc % nvt ( :, : ) % v*nhc % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhc % num_nhc
           nhc % v_scale(n) = nhc % v_scale(n)*EXP(-0.5_dp*nhc % dt_yosh( &
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
              scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh) &
                   )
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
              nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                   nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh) ! shift
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
             nhc % nvt(nhc % nhc_len,:) % f*0.25_dp*nhc % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

! now scale the particle velocities
  ii = 0
  DO ikind = 1, SIZE ( molecule_kind_set )
    nmol_local = local_molecules % n_el ( ikind )
    DO imol_local = 1, nmol_local
      imol = local_molecules % list ( ikind ) % array ( imol_local )
      molecule => molecule_set ( imol )
      CALL get_molecule ( molecule, first_atom = first_atom, &
                          last_atom = last_atom ) 
      DO ipart = first_atom, last_atom
        ii = ii + 1
        vel ( 1, ipart ) =  vel ( 1, ipart ) * nhc % p_scale ( 1, ii ) % point
        vel ( 2, ipart ) =  vel ( 2, ipart ) * nhc % p_scale ( 2, ii ) % point
        vel ( 3, ipart ) =  vel ( 3, ipart ) * nhc % p_scale ( 3, ii ) % point
      END DO
    END DO
  END DO
  flops = flops + REAL ( npart_local*3,KIND=dp)

  flops = flops * 1.E-6_dp
  CALL timestop(flops,handle)

END SUBROUTINE lnhc_particles



!******************************************************************************

END MODULE extended_system_dynamics

!******************************************************************************
