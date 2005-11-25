!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/init_extended_system_variables [1.0] *
!!
!!   NAME
!!     init_extended_system_variables
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM 12-NOV-2003:  new
!!
!!   SOURCE
!******************************************************************************
MODULE init_extended_system_variables
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum, mp_allgather
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: gasdev
  USE input_constants,                 ONLY: nve_ensemble,&
                                             nvt_ensemble,&
                                             npt_i_ensemble,&
                                             npt_f_ensemble,&
                                             nph_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             isokin_ensemble
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_nhc_variables, init_barostat_variables, init_nhc_forces
  
  CONTAINS
!******************************************************************************
!!****** init_extended_system_variables/init_nhc_variables [1.0] *
!!
!!   NAME
!!     init_nhc_variables
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

SUBROUTINE init_nhc_variables ( nhc,tau_nhc,temp_ext, globenv )
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    REAL(KIND=dp), INTENT(IN)                :: tau_nhc, temp_ext
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

! locals
    INTEGER                                  :: i, j, number, isos, iproc, &
                                                rannumneed, tot_rn, &
                                                icount
    REAL(KIND=dp)                            :: akin, temp, v
    INTEGER, DIMENSION ( : ), ALLOCATABLE    :: rn_need_per_proc
    REAL (KIND=dp), DIMENSION ( : ), ALLOCATABLE    :: array_of_rn

! first initializing the mass of the nhc variables

  nhc % nvt ( :, : ) % mass = nhc % nvt ( :, : ) % nkt * tau_nhc * tau_nhc
  nhc % nvt ( :, : ) % eta = 0._dp
  nhc % nvt ( :, : ) % v = 0._dp
  nhc % nvt ( :, : ) % f = 0._dp


  SELECT CASE ( nhc % dis_type )
    CASE ( 'REP_REP' )  ! for NPT 
    CASE DEFAULT
      rannumneed = nhc%num_nhc*nhc%nhc_len
      ALLOCATE ( rn_need_per_proc( globenv % num_pe ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                                    'failed to allocate rn_need_per_proc')
      rn_need_per_proc ( : ) = 0

      CALL mp_allgather ( rannumneed, rn_need_per_proc, globenv % group )

      IF (nhc % dis_type=='INTER_REP') THEN
        tot_rn = rn_need_per_proc(1)
      ELSE
        tot_rn = SUM ( rn_need_per_proc )
      ENDIF

      DEALLOCATE ( rn_need_per_proc, STAT = isos )
      IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                                    'failed to deallocate rn_need_per_proc')

      ALLOCATE(array_of_rn( tot_rn), STAT = isos)
      IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                                    'failed to allocate array_of_rn')
      array_of_rn(:) = 0.0_dp
  END SELECT 

  SELECT CASE ( nhc % dis_type )
    CASE ( 'REP_REP' ) ! for NPT
      IF ( globenv % ionode ) THEN
! Map deterministically determined random number to nhc % v
        DO i = 1, nhc % num_nhc
           DO j = 1, nhc % nhc_len
              nhc % nvt(j,i) % v = gasdev ( globenv % idum )
           END DO
        END DO
  
        akin = 0.0_dp
        DO i = 1, nhc % num_nhc
           DO j = 1, nhc % nhc_len
              akin = akin + 0.5_dp * ( nhc % nvt ( j , i ) % mass * &
                                       nhc % nvt ( j , i ) % v * &
                                       nhc % nvt ( j , i ) % v )
           END DO
        END DO
        number = nhc % num_nhc

! scale velocities to get the correct initial temperature
        temp = 2.0_dp*akin/REAL(number)
        temp = sqrt(temp_ext/temp)
        DO i = 1, nhc % num_nhc
           DO j = 1, nhc % nhc_len
              nhc % nvt(j,i) % v = temp * nhc % nvt(j,i) % v
              nhc % nvt(j,i) % eta = 0.0_dp
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
      END IF
! Broadcast
      DO i = 1, SIZE ( nhc % nvt, 1 )
        DO j = 1, SIZE ( nhc % nvt, 2 )
          CALL mp_bcast ( nhc % nvt ( i, j ) % v, globenv % source, globenv % group )
          CALL mp_bcast ( nhc % nvt ( i, j ) % eta, globenv % source, globenv % group )
          CALL mp_bcast ( nhc % nvt ( i, j ) % f, globenv % source, globenv % group )
        END DO
      END DO

    CASE DEFAULT
      DO i=1,tot_rn
        array_of_rn(i) = gasdev ( globenv % idum )
      END DO
! Map deterministically determined random number to nhc % v
      DO i = 1, nhc % num_nhc
         icount = nhc % index(i)
         icount = (icount-1) * nhc % nhc_len
         DO j = 1, nhc % nhc_len
            icount = icount + 1
            nhc % nvt(j,i) % v = array_of_rn ( icount )
         END DO
      END DO

      DEALLOCATE ( array_of_rn, STAT = isos )
      IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                                'failed to deallocate array_of_rn')
  
      akin = 0.0_dp
      DO i = 1, nhc % num_nhc
         DO j = 1, nhc % nhc_len
            akin = akin + 0.5_dp * ( nhc % nvt ( j , i ) % mass * &
                                      nhc % nvt ( j , i ) % v * &
                                      nhc % nvt ( j , i ) % v )
         END DO
      END DO
      number = nhc % num_nhc

      IF (nhc % dis_type=='DIS_REP') CALL mp_sum(akin,globenv % group)
      IF (nhc % dis_type=='DIS_REP') CALL mp_sum(number,globenv % group)

! scale velocities to get the correct initial temperature
      temp = 2.0_dp*akin/REAL(number)
      temp = sqrt(temp_ext/temp)
      DO i = 1, nhc % num_nhc
         DO j = 1, nhc % nhc_len
            nhc % nvt(j,i) % v = temp * nhc % nvt(j,i) % v
            nhc % nvt(j,i) % eta = 0.0_dp
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

  END SELECT

END SUBROUTINE init_nhc_variables

!!*****
!******************************************************************************
!!****** init_extended_system_variables/init_barostat_variables [1.0] *
!!
!!   NAME
!!     init_barostat_variables
!!
!!   FUNCTION
!!     Initializes the barostat velocities to the Maxwellian distribution
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

SUBROUTINE init_barostat_variables ( npt, tau_cell, temp_ext, deth, nfree, ensemble, &
                                     cmass, globenv )


    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    REAL(KIND=dp), INTENT(IN)                :: tau_cell, temp_ext, deth, cmass
    INTEGER, INTENT(IN)                      :: nfree
    INTEGER, INTENT(IN)                      :: ensemble
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    INTEGER                                  :: i, j, number
    REAL(KIND=dp)                            :: akin, temp, v

!------------------------------------------------------------------------------
! first initializing the mass of the nhc variables

  SELECT CASE ( ensemble )
  CASE ( npt_i_ensemble )
     npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) &
          * temp_ext * tau_cell ** 2
  CASE ( npt_f_ensemble)
     npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) &
          * temp_ext * tau_cell ** 2 / 3.0_dp
  CASE ( nph_uniaxial_ensemble,nph_uniaxial_damped_ensemble )
     npt ( :, : ) % mass = cmass
          
!     npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) &
!          * temp_ext * tau_cell ** 2
  END SELECT

  IF ( globenv % ionode ) THEN
! initializing velocities
    DO i = 1, size ( npt,1)
       DO j = i, size ( npt,2)
          v = gasdev ( globenv % idum)
! symmetrizing the initial barostat velocities to ensure
! no rotation of the cell under NPT_F
          npt(j,i) % v = v
          npt(i,j) % v = v
       END DO
    END DO

    akin = 0.0_dp
    DO i = 1, SIZE( npt,1 )
       DO j = 1, SIZE( npt,2 )
          akin = akin + 0.5_dp*( npt(j,i) % mass* npt(j,i) % v* npt &
               (j,i) % v)
       END DO
    END DO

    number = SIZE ( npt, 1 ) * SIZE ( npt, 2 )

! scale velocities to get the correct initial temperature
    IF ( number /= 0 ) THEN
      temp = 2.0_dp * akin / REAL ( number,KIND=dp)
      temp = sqrt(temp_ext/temp)
    ENDIF
    DO i = 1, size(npt,1)
       DO j = i, size(npt,2)
          npt(j,i) % v = temp * npt(j,i) % v
          npt(i,j) % v = npt(j,i) % v
! dbg ISOTROPIC LIMIT
!       npt(j,i) % v = 0.0_dp
!       npt(i,j) % v = 0.0_dp
!       write (*, *) 'INTIAL v_eps', npt(j,i) % v
! dbg 
       END DO
    END DO

! Zero barostat velocities for nph_uniaxial
    SELECT CASE ( ensemble )
    CASE ( nph_uniaxial_ensemble, nph_uniaxial_damped_ensemble )
       npt(:,:) % v = 0.0_dp
    END SELECT
  ENDIF 

! Broadcast
  DO i = 1, SIZE ( npt, 1 )
    DO j = 1, SIZE ( npt, 2 )
       CALL mp_bcast ( npt ( i, j ) % v, globenv % source, globenv % group )
    END DO
  END DO

END SUBROUTINE init_barostat_variables
!!*****
!******************************************************************************
!!****** init_extended_system_variables/init_nhc_forces [1.0] *
!!
!!   NAME
!!     init_nhc_forces
!!
!!   FUNCTION
!!    Assigns extended parameters from the restart file.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE init_nhc_forces ( nhc_part, nhc_baro, nhc_coef )
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT), OPTIONAL                :: nhc_part, nhc_baro, nhc_coef

    INTEGER                                  :: i, j

! assign the forces

IF(PRESENT(nhc_part)) THEN
  DO i = 1, SIZE ( nhc_part % nvt, 2 )
     DO j = 2, SIZE ( nhc_part % nvt, 1 )
        nhc_part % nvt ( j, i ) % f = nhc_part % nvt ( j-1, i ) % mass * &
                                      nhc_part % nvt ( j-1, i ) % v *  &
                                      nhc_part % nvt ( j-1, i ) % v -  &
                                      nhc_part % nvt ( j, i ) % nkt
        nhc_part % nvt ( j, i ) % f = nhc_part % nvt ( j, i ) % f &
                                     /nhc_part % nvt ( j, i ) % mass
     END DO
  END DO
END IF

! assign the forces

IF(PRESENT(nhc_baro)) THEN
  DO i = 1, SIZE ( nhc_baro % nvt, 2 )
     DO j = 2, SIZE ( nhc_baro % nvt, 1 )
        nhc_baro % nvt ( j, i ) % f = nhc_baro % nvt ( j-1, i ) % mass * &
                                      nhc_baro % nvt ( j-1, i ) % v *  &
                                      nhc_baro % nvt ( j-1, i ) % v -  &
                                      nhc_baro % nvt ( j, i ) % nkt
        nhc_baro % nvt ( j, i ) % f = nhc_baro % nvt ( j, i ) % f &
                                     /nhc_baro % nvt ( j, i ) % mass
     END DO
  END DO
END IF

! assign the forces

IF(PRESENT(nhc_coef)) THEN
  DO i = 1, SIZE ( nhc_coef % nvt, 2 )
     DO j = 2, SIZE ( nhc_coef % nvt, 1 )
        nhc_coef % nvt ( j, i ) % f = nhc_coef % nvt ( j-1, i ) % mass * &
                                      nhc_coef % nvt ( j-1, i ) % v *  &
                                      nhc_coef % nvt ( j-1, i ) % v -  &
                                      nhc_coef % nvt ( j, i ) % nkt
        nhc_coef % nvt ( j, i ) % f = nhc_coef % nvt ( j, i ) % f &
                                     /nhc_coef % nvt ( j, i ) % mass
     END DO
  END DO
END IF

END SUBROUTINE init_nhc_forces

!!*****
!******************************************************************************
END MODULE init_extended_system_variables
!******************************************************************************
