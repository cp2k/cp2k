!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/extended_system_initialize [1.0] *
!!
!!   NAME
!!     extended_system_initialize
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
MODULE extended_system_initialize
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

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_npt, initialize_nhc
  PUBLIC :: set_yoshida_coef
  
  CONTAINS
!******************************************************************************
!!****** extended_system_initialize/initialize_nhc [1.0] *
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
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    REAL(KIND=dp), INTENT(IN)                    :: tau_nhc, temp_ext
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

! locals
    INTEGER                                  :: i, j, number, isos, iproc, &
                                                rannumneed, tot_rn, &
                                                icount
    REAL(KIND=dp)                                :: akin, temp, v
    INTEGER, DIMENSION ( : ), ALLOCATABLE    :: rn_need_per_proc
    REAL (KIND=dp), DIMENSION ( : ), ALLOCATABLE    :: array_of_rn

! first initializing the mass of the nhc variables

  nhc % nvt ( :, : ) % mass = nhc % nvt ( :, : ) % nkt * tau_nhc * tau_nhc

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

END SUBROUTINE initialize_nhc

!!*****
!******************************************************************************
!!****** extended_system_initialize/initialize_npt [1.0] *
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


    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    REAL(KIND=dp), INTENT(IN)                    :: tau_cell, temp_ext, deth
    INTEGER, INTENT(IN)                      :: nfree
    CHARACTER(LEN=20), INTENT(IN)            :: ensemble
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    INTEGER                                  :: i, j, number
    REAL(KIND=dp)                                :: akin, temp, v

!------------------------------------------------------------------------------
! first initializing the mass of the nhc variables

  SELECT CASE ( ensemble )
  CASE ( 'NPT_I')
     npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) &
          * temp_ext * tau_cell ** 2
  CASE ( 'NPT_F')
     npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) &
          * temp_ext * tau_cell ** 2 / 3.0_dp
  END SELECT

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

END SUBROUTINE initialize_npt

!!*****
!******************************************************************************
!!****** extended_system_initialize/set_yoshida_coef [1.0] *
!!
!!   NAME
!!     set_yoshida_coef
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

SUBROUTINE set_yoshida_coef ( nhc, dt )


    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    REAL(KIND=dp), INTENT(IN)                    :: dt

    REAL(KIND=dp), DIMENSION(nhc%nyosh)          :: yosh_wt

!------------------------------------------------------------------------------

  SELECT CASE (nhc % nyosh)
  CASE DEFAULT
     CALL stop_program ( 'yoshida_coef', 'Value not available' )
  CASE (1)
     yosh_wt(1) = 1.0_dp
  CASE (3)
     yosh_wt(1) = 1.0_dp/(2.0_dp-(2.0_dp)**(1.0_dp/3.0_dp))
     yosh_wt(2) = 1.0_dp - 2.0_dp*yosh_wt(1)
     yosh_wt(3) = yosh_wt(1)
  CASE (5)
     yosh_wt(1) = 1.0_dp/(4.0_dp-(4.0_dp)**(1.0_dp/3.0_dp))
     yosh_wt(2) = yosh_wt(1)
     yosh_wt(4) = yosh_wt(1)
     yosh_wt(5) = yosh_wt(1)
     yosh_wt(3) = 1.0_dp - 4.0_dp*yosh_wt(1)
  CASE (7)
     yosh_wt(1) = .78451361047756_dp
     yosh_wt(2) = .235573213359357_dp
     yosh_wt(3) = -1.17767998417887_dp
     yosh_wt(4) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+yosh_wt(3))
     yosh_wt(5) = yosh_wt(3)
     yosh_wt(6) = yosh_wt(2)
     yosh_wt(7) = yosh_wt(1)
  CASE (9)
     yosh_wt(1) = 0.192_dp
     yosh_wt(2) = 0.554910818409783619692725006662999_dp
     yosh_wt(3) = 0.124659619941888644216504240951585_dp
     yosh_wt(4) = -0.843182063596933505315033808282941_dp
     yosh_wt(5) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+&
                  yosh_wt(3)+yosh_wt(4))
     yosh_wt(6) = yosh_wt(4)
     yosh_wt(7) = yosh_wt(3)
     yosh_wt(8) = yosh_wt(2)
     yosh_wt(9) = yosh_wt(1)
  CASE (15)
     yosh_wt(1) = 0.102799849391985_dp
     yosh_wt(2) = -0.196061023297549e1_dp
     yosh_wt(3) = 0.193813913762276e1_dp
     yosh_wt(4) = -0.158240635368243_dp
     yosh_wt(5) = -0.144485223686048e1_dp
     yosh_wt(6) = 0.253693336566229_dp
     yosh_wt(7) = 0.914844246229740_dp
     yosh_wt(8) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+&
                 yosh_wt(3)+yosh_wt(4)+yosh_wt(5)+yosh_wt(6)+yosh_wt(7))
     yosh_wt(9) = yosh_wt(7)
     yosh_wt(10) = yosh_wt(6)
     yosh_wt(11) = yosh_wt(5)
     yosh_wt(12) = yosh_wt(4)
     yosh_wt(13) = yosh_wt(3)
     yosh_wt(14) = yosh_wt(2)
     yosh_wt(15) = yosh_wt(1)
  END SELECT
  nhc % dt_yosh = dt * yosh_wt / REAL ( nhc % nc,KIND=dp)

END SUBROUTINE set_yoshida_coef

!!*****
!******************************************************************************
END MODULE extended_system_initialize
!******************************************************************************
