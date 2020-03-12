!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods to apply the QTB thermostat to PI runs.
!>         Based on the PILE implementation from Felix Uhl (pint_pile.F)
!> \author Fabien Brieuc
!> \par History
!>      02.2018 created [Fabien Brieuc]
! **************************************************************************************************
MODULE pint_qtb
   USE cp_files,                        ONLY: open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level,&
                                              silent_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE fft_lib,                         ONLY: fft_1dm,&
                                              fft_create_plan_1dm,&
                                              fft_destroy_plan
   USE fft_plan,                        ONLY: fft_plan_type
   USE fft_tools,                       ONLY: FWFFT,&
                                              fft_plan_style,&
                                              fft_type
   USE input_constants,                 ONLY: propagator_rpmd
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_bcast
   USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                              rng_record_length,&
                                              rng_stream_type,&
                                              rng_stream_type_from_record
   USE pint_io,                         ONLY: pint_write_line
   USE pint_types,                      ONLY: normalmode_env_type,&
                                              pint_env_type,&
                                              qtb_therm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pint_qtb_step, &
             pint_qtb_init, &
             pint_qtb_release, &
             pint_calc_qtb_energy

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_qtb'

CONTAINS

   ! ***************************************************************************
   !> \brief initializes the data for a QTB run
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param qtb_therm ...
!> \param pint_env ...
!> \param normalmode_env ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_init(qtb_therm, pint_env, normalmode_env, section)
      TYPE(qtb_therm_type), POINTER                      :: qtb_therm
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(normalmode_env_type), POINTER                 :: normalmode_env
      TYPE(section_vals_type), POINTER                   :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_init', routineP = moduleN//':'//routineN

      CHARACTER(LEN=rng_record_length)                   :: rng_record
      INTEGER                                            :: i, j, p
      LOGICAL                                            :: restart
      REAL(KIND=dp)                                      :: dti2, ex
      REAL(KIND=dp), DIMENSION(3, 2)                     :: initial_seed
      TYPE(section_vals_type), POINTER                   :: rng_section

      IF (pint_env%propagator%prop_kind /= propagator_rpmd) THEN
         CPABORT("QTB is designed to work with the RPMD propagator only")
      ENDIF

      pint_env%e_qtb = 0.0_dp
      ALLOCATE (qtb_therm)
      qtb_therm%ref_count = 1
      qtb_therm%thermostat_energy = 0.0_dp

      !Get input parameters
      CALL section_vals_val_get(section, "TAU", r_val=qtb_therm%tau)
      CALL section_vals_val_get(section, "LAMBDA", r_val=qtb_therm%lamb)
      CALL section_vals_val_get(section, "TAUCUT", r_val=qtb_therm%taucut)
      CALL section_vals_val_get(section, "LAMBCUT", r_val=qtb_therm%lambcut)
      CALL section_vals_val_get(section, "FP", i_val=qtb_therm%fp)
      CALL section_vals_val_get(section, "NF", i_val=qtb_therm%nf)

      p = pint_env%p
      dti2 = 0.5_dp*pint_env%dt
      ALLOCATE (qtb_therm%c1(p))
      ALLOCATE (qtb_therm%c2(p))
      ALLOCATE (qtb_therm%g_fric(p))
      ALLOCATE (qtb_therm%massfact(p, pint_env%ndim))

      !Initialize everything
      qtb_therm%g_fric(1) = 1.0_dp/qtb_therm%tau
      DO i = 2, p
         qtb_therm%g_fric(i) = SQRT((1.d0/qtb_therm%tau)**2 + (qtb_therm%lamb)**2* &
                                    normalmode_env%lambda(i))
      END DO
      DO i = 1, p
         ex = -dti2*qtb_therm%g_fric(i)
         qtb_therm%c1(i) = EXP(ex)
         ex = qtb_therm%c1(i)*qtb_therm%c1(i)
         qtb_therm%c2(i) = SQRT(1.0_dp - ex)
      END DO
      DO j = 1, pint_env%ndim
         DO i = 1, pint_env%p
            qtb_therm%massfact(i, j) = SQRT(1.0_dp/pint_env%mass_fict(i, j))
         END DO
      END DO

      !prepare Random number generator
      NULLIFY (rng_section)
      rng_section => section_vals_get_subs_vals(section, &
                                                subsection_name="RNG_INIT")
      CALL section_vals_get(rng_section, explicit=restart)
      IF (restart) THEN
         CALL section_vals_val_get(rng_section, "_DEFAULT_KEYWORD_", &
                                   i_rep_val=1, c_val=rng_record)
         qtb_therm%gaussian_rng_stream = rng_stream_type_from_record(rng_record)
      ELSE
         initial_seed(:, :) = REAL(pint_env%thermostat_rng_seed, dp)
         qtb_therm%gaussian_rng_stream = rng_stream_type( &
                                         name="qtb_rng_gaussian", distribution_type=GAUSSIAN, &
                                         extended_precision=.TRUE., &
                                         seed=initial_seed)
      END IF

      !Initialization of the QTB random forces
      CALL pint_qtb_forces_init(pint_env, normalmode_env, qtb_therm, restart)

   END SUBROUTINE pint_qtb_init

! **************************************************************************************************
!> \brief ...
!> \param vold ...
!> \param vnew ...
!> \param p ...
!> \param ndim ...
!> \param masses ...
!> \param qtb_therm ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_step(vold, vnew, p, ndim, masses, qtb_therm)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vold, vnew
      INTEGER, INTENT(IN)                                :: p, ndim
      REAL(kind=dp), DIMENSION(:, :), INTENT(IN)         :: masses
      TYPE(qtb_therm_type), POINTER                      :: qtb_therm

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_step', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ibead, idim
      REAL(KIND=dp)                                      :: delta_ekin

      CALL timeset(routineN, handle)
      delta_ekin = 0.0_dp

      !update random forces
      DO ibead = 1, p
         qtb_therm%cpt(ibead) = qtb_therm%cpt(ibead) + 1
         !new random forces at every qtb_therm%step
         IF (qtb_therm%cpt(ibead) == 2*qtb_therm%step(ibead)) THEN
            IF (ibead == 1) THEN
               !update the rng status
               DO i = 1, qtb_therm%nf - 1
                  qtb_therm%rng_status(i) = qtb_therm%rng_status(i + 1)
               END DO
               CALL qtb_therm%gaussian_rng_stream%dump(qtb_therm%rng_status(qtb_therm%nf))
            END IF
            DO idim = 1, ndim
               !update random numbers
               DO i = 1, qtb_therm%nf - 1
                  qtb_therm%r(i, ibead, idim) = qtb_therm%r(i + 1, ibead, idim)
               END DO
               qtb_therm%r(qtb_therm%nf, ibead, idim) = qtb_therm%gaussian_rng_stream%next()
               !compute new random force through the convolution product
               qtb_therm%rf(ibead, idim) = 0.0_dp
               DO i = 1, qtb_therm%nf
                  qtb_therm%rf(ibead, idim) = qtb_therm%rf(ibead, idim) + &
                                              qtb_therm%h(i, ibead)*qtb_therm%r(i, ibead, idim)
               END DO
            END DO
            qtb_therm%cpt(ibead) = 0
         END IF
      END DO

      !perform MD step
      DO idim = 1, ndim
         DO ibead = 1, p
            vnew(ibead, idim) = qtb_therm%c1(ibead)*vold(ibead, idim) + &
                                qtb_therm%massfact(ibead, idim)*qtb_therm%c2(ibead)* &
                                qtb_therm%rf(ibead, idim)
            delta_ekin = delta_ekin + masses(ibead, idim)*( &
                         vnew(ibead, idim)*vnew(ibead, idim) - &
                         vold(ibead, idim)*vold(ibead, idim))
         END DO
      END DO

      qtb_therm%thermostat_energy = qtb_therm%thermostat_energy - 0.5_dp*delta_ekin

      CALL timestop(handle)
   END SUBROUTINE pint_qtb_step

   ! ***************************************************************************
   !> \brief releases the qtb environment
   !> \param qtb_therm qtb data to be released
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param qtb_therm ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_release(qtb_therm)

      TYPE(qtb_therm_type), POINTER                      :: qtb_therm

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(qtb_therm)) THEN
         qtb_therm%ref_count = qtb_therm%ref_count - 1
         IF (qtb_therm%ref_count == 0) THEN
            DEALLOCATE (qtb_therm%c1)
            DEALLOCATE (qtb_therm%c2)
            DEALLOCATE (qtb_therm%g_fric)
            DEALLOCATE (qtb_therm%massfact)
            DEALLOCATE (qtb_therm%rf)
            DEALLOCATE (qtb_therm%h)
            DEALLOCATE (qtb_therm%r)
            DEALLOCATE (qtb_therm%cpt)
            DEALLOCATE (qtb_therm%step)
            DEALLOCATE (qtb_therm%rng_status)
            DEALLOCATE (qtb_therm)
         END IF
      END IF
      NULLIFY (qtb_therm)

      RETURN
   END SUBROUTINE pint_qtb_release

   ! ***************************************************************************
   !> \brief returns the qtb kinetic energy contribution
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
! **************************************************************************************************
   SUBROUTINE pint_calc_qtb_energy(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_qtb_energy', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(pint_env%qtb_therm)) THEN
         pint_env%e_qtb = pint_env%qtb_therm%thermostat_energy
      END IF

      RETURN

   END SUBROUTINE pint_calc_qtb_energy

   ! ***************************************************************************
   !> \brief initialize the QTB random forces
   !> \author Fabien Brieuc
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param normalmode_env ...
!> \param qtb_therm ...
!> \param restart ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_forces_init(pint_env, normalmode_env, qtb_therm, restart)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(normalmode_env_type), POINTER                 :: normalmode_env
      TYPE(qtb_therm_type), POINTER                      :: qtb_therm
      LOGICAL                                            :: restart

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_forces_init', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: tmp1
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: filter
      INTEGER                                            :: handle, i, ibead, idim, log_unit, ndim, &
                                                            nf, p, print_level, step
      REAL(KIND=dp)                                      :: aa, bb, correct, dt, dw, fcut, h, kT, &
                                                            tmp, w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: fp1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(fft_plan_type)                                :: plan

      CALL timeset(routineN, handle)

      IF (fft_type /= 3) CALL cp_warn(__LOCATION__, "The FFT library in use cannot"// &
                                      " handle transformation of an arbitrary length.")

      p = pint_env%p
      ndim = pint_env%ndim
      dt = pint_env%dt
      IF (MOD(qtb_therm%nf, 2) /= 0) qtb_therm%nf = qtb_therm%nf + 1
      nf = qtb_therm%nf

      para_env => pint_env%logger%para_env

      ALLOCATE (qtb_therm%rng_status(nf))
      ALLOCATE (qtb_therm%h(nf, p))
      ALLOCATE (qtb_therm%step(p))

      !initialize random forces on ionode only
      IF (para_env%ionode) THEN

         NULLIFY (logger)
         logger => cp_get_default_logger()
         print_level = logger%iter_info%print_level

         !physical temperature (T) not the simulation one (TxP)
         kT = pint_env%kT*pint_env%propagator%temp_sim2phys

         ALLOCATE (fp(nf/2))
         ALLOCATE (filter(0:nf - 1))

         IF (print_level == debug_print_level) THEN
            !create log file if print_level is debug
            CALL open_file(file_name=TRIM(logger%iter_info%project_name)//".qtbLog", &
                           file_action="WRITE", file_status="UNKNOWN", unit_number=log_unit)
            WRITE (log_unit, '(A)') ' # Log file for the QTB random forces generation'
            WRITE (log_unit, '(A)') ' # ------------------------------------------------'
            WRITE (log_unit, '(A,I5)') ' # Number of beads P = ', p
            WRITE (log_unit, '(A,I6)') ' # Number of dimension 3*N = ', ndim
            WRITE (log_unit, '(A,I4)') ' # Number of filter parameters Nf=', nf
         END IF

         DO ibead = 1, p
            !fcut is adapted to the NM freq.
            !Note that lambda is the angular free ring freq. squared
            fcut = SQRT((1.d0/qtb_therm%taucut)**2 + (qtb_therm%lambcut)**2* &
                        normalmode_env%lambda(ibead))
            fcut = fcut/twopi
            !new random forces are drawn every step
            qtb_therm%step(ibead) = NINT(1.0_dp/(2.0_dp*fcut*dt))
            IF (qtb_therm%step(ibead) == 0) qtb_therm%step(ibead) = 1
            step = qtb_therm%step(ibead)
            !effective timestep h = step*dt = 1/(2*fcut)
            h = step*dt
            !angular freq. step - dw = 2*pi/(nf*h) = 2*wcut/nf
            dw = twopi/(nf*h)

            !generate f_P function
            IF (qtb_therm%fp == 0) THEN
               CALL pint_qtb_computefp0(pint_env, fp, fp1, dw, aa, bb, log_unit, ibead, print_level)
            ELSE
               CALL pint_qtb_computefp1(pint_env, fp, fp1, dw, aa, bb, log_unit, ibead, print_level)
            END IF
            fp = p*kT*fp ! fp is now in cp2k energy units

            IF (print_level == debug_print_level) THEN
               WRITE (log_unit, '(A,I4,A)') ' # --------  NM ', ibead, '  --------'
               WRITE (log_unit, '(A,I4,A)') ' # New random forces every ', step, ' MD steps'
               WRITE (log_unit, '(A,ES13.3,A)') ' # Angular cutoff freq. = ', twopi*fcut*4.1341e4_dp, ' rad/ps'
               WRITE (log_unit, '(A,ES13.3,A)') ' # Free ring polymer angular freq.= ', &
                  SQRT(normalmode_env%lambda(ibead))*4.1341e4_dp, ' rad/ps'
               WRITE (log_unit, '(A,ES13.3,A)') ' # Friction coeff. = ', qtb_therm%g_fric(ibead)*4.1341e4_dp, ' THz'
               WRITE (log_unit, '(A,ES13.3,A)') ' # Angular frequency step dw = ', dw*4.1341e4_dp, ' rad/ps'
            END IF

            !compute the filter in Fourier space
            IF (p == 1) THEN
               filter(0) = SQRT(kT)*(1.0_dp, 0.0_dp)
            ELSE IF (qtb_therm%fp == 1 .AND. ibead == 1) THEN
               filter(0) = SQRT(p*kT)*(1.0_dp, 0.0_dp)
            ELSE
               filter(0) = SQRT(p*kT*fp1(1))*(1.0_dp, 0.0_dp)
            ENDIF
            DO i = 1, nf/2
               w = i*dw
               tmp = 0.5_dp*w*h
               correct = SIN(tmp)/tmp
               filter(i) = SQRT(fp(i))/correct*(1.0_dp, 0.0_dp)
               filter(nf - i) = CONJG(filter(i))
            END DO

            !compute the filter in time space - FFT
            CALL pint_qtb_fft(filter, filter, plan, nf)
            !reordering + normalisation
            !normalisation : 1/nf comes from the DFT, 1/sqrt(step) is to
            !take into account the effective timestep h = step*dt and
            !1/sqrt(2.0_dp) is to take into account the fact that the
            !same random force is used for the two thermostat "half-steps"
            DO i = 0, nf/2 - 1
               tmp1 = filter(i)/(nf*SQRT(2.0_dp*step))
               filter(i) = filter(nf/2 + i)/(nf*SQRT(2.0_dp*step))
               filter(nf/2 + i) = tmp1
            END DO

            DO i = 0, nf - 1
               qtb_therm%h(i + 1, ibead) = REAL(filter(i), dp)
            END DO
         END DO

         DEALLOCATE (filter)
         DEALLOCATE (fp)
         IF (p > 1) DEALLOCATE (fp1)
      END IF

      CALL mp_bcast(qtb_therm%h, para_env%source, para_env%group)
      CALL mp_bcast(qtb_therm%step, para_env%source, para_env%group)

      ALLOCATE (qtb_therm%r(nf, p, ndim))
      ALLOCATE (qtb_therm%cpt(p))
      ALLOCATE (qtb_therm%rf(p, ndim))

      IF (restart) THEN
         CALL pint_qtb_restart(pint_env, qtb_therm)
      ELSE
         !update the rng status
         DO i = 1, qtb_therm%nf
            CALL qtb_therm%gaussian_rng_stream%dump(qtb_therm%rng_status(i))
         END DO
         !if no restart then initialize random numbers from scratch
         qtb_therm%cpt = 0
         DO idim = 1, ndim
            DO ibead = 1, p
               DO i = 1, nf
                  qtb_therm%r(i, ibead, idim) = qtb_therm%gaussian_rng_stream%next()
               END DO
            END DO
         END DO
      END IF

      !compute the first random forces
      DO idim = 1, ndim
         DO ibead = 1, p
            qtb_therm%rf(ibead, idim) = 0.0_dp
            DO i = 1, nf
               qtb_therm%rf(ibead, idim) = qtb_therm%rf(ibead, idim) + &
                                           qtb_therm%h(i, ibead)*qtb_therm%r(i, ibead, idim)
            END DO
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE pint_qtb_forces_init

   ! ***************************************************************************
   !> \brief control the generation of the first random forces in the case
   !> of a restart
   !> \author Fabien Brieuc
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param qtb_therm ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_restart(pint_env, qtb_therm)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(qtb_therm_type), POINTER                      :: qtb_therm

      INTEGER                                            :: begin, i, ibead, idim, istep

      begin = pint_env%first_step - MOD(pint_env%first_step, qtb_therm%step(1)) - &
              (qtb_therm%nf - 1)*qtb_therm%step(1)

      IF (begin <= 0) THEN
         qtb_therm%cpt = 0
         !update the rng status
         DO i = 1, qtb_therm%nf
            CALL qtb_therm%gaussian_rng_stream%dump(qtb_therm%rng_status(i))
         END DO
         !first random numbers initialized from scratch
         DO idim = 1, pint_env%ndim
            DO ibead = 1, pint_env%p
               DO i = 1, qtb_therm%nf
                  qtb_therm%r(i, ibead, idim) = qtb_therm%gaussian_rng_stream%next()
               END DO
            END DO
         END DO
         begin = 1
      ELSE
         qtb_therm%cpt(1) = 2*(qtb_therm%step(1) - 1)
         DO ibead = 2, pint_env%p
            qtb_therm%cpt(ibead) = 2*MOD(begin - 1, qtb_therm%step(ibead))
         END DO
      END IF

      !from istep = 1,2*(the last previous MD step - begin) because
      !the thermostat step is called two times per MD step
      !DO istep = 2*begin, 2*pint_env%first_step
      DO istep = 1, 2*(pint_env%first_step - begin + 1)
         DO ibead = 1, pint_env%p
            qtb_therm%cpt(ibead) = qtb_therm%cpt(ibead) + 1
            !new random forces at every qtb_therm%step
            IF (qtb_therm%cpt(ibead) == 2*qtb_therm%step(ibead)) THEN
               IF (ibead == 1) THEN
                  !update the rng status
                  DO i = 1, qtb_therm%nf - 1
                     qtb_therm%rng_status(i) = qtb_therm%rng_status(i + 1)
                  END DO
                  CALL qtb_therm%gaussian_rng_stream%dump(qtb_therm%rng_status(qtb_therm%nf))
               END IF
               DO idim = 1, pint_env%ndim
                  !update random numbers
                  DO i = 1, qtb_therm%nf - 1
                     qtb_therm%r(i, ibead, idim) = qtb_therm%r(i + 1, ibead, idim)
                  END DO
                  qtb_therm%r(qtb_therm%nf, ibead, idim) = qtb_therm%gaussian_rng_stream%next()
               END DO
               qtb_therm%cpt(ibead) = 0
            END IF
         END DO
      END DO

   END SUBROUTINE pint_qtb_restart

   ! ***************************************************************************
   !> \brief compute the f_P^(0) function necessary for coupling QTB with PIMD
   !> \param fp stores the computed function on the grid used for the generation
   !> of the filter h
   !> \param fp1 stores the computed function on an larger and finer grid
   !> \param dw angular frequency step
   !> \author Fabien Brieuc
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param fp ..
!> \param fp1 ..
!> \param dw ...
!> \param aa ...
!> \param bb ...
!> \param log_unit ...
!> \param ibead ...
!> \param print_level ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_computefp0(pint_env, fp, fp1, dw, aa, bb, log_unit, ibead, print_level)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(KIND=dp), DIMENSION(:)                        :: fp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: fp1
      REAL(KIND=dp)                                      :: dw, aa, bb
      INTEGER                                            :: log_unit, ibead, print_level

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_computefp0', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=200)                                 :: line
      INTEGER                                            :: i, j, k, n, niter, nx, p
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: kk
      REAL(KIND=dp)                                      :: dx, dx1, err, fprev, hbokT, malpha, op, &
                                                            r2, tmp, w, x1, xmax, xmin, xx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: h, x, x2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fpxk, xk, xk2

      n = SIZE(fp)
      p = pint_env%p

      !using the physical temperature (T) not the simulation one (TxP)
      hbokT = 1.0_dp/(pint_env%kT*pint_env%propagator%temp_sim2phys)

      !P = 1 : standard QTB
      !fp = theta(w, T) / kT
      IF (p == 1) THEN
         DO j = 1, n
            w = j*dw
            tmp = hbokT*w
            fp(j) = tmp*(0.5_dp + 1.0_dp/(EXP(tmp) - 1.0_dp))
         END DO

         IF (print_level == debug_print_level) THEN
            WRITE (log_unit, '(A)') ' # ------------------------------------------------'
            WRITE (log_unit, '(A)') ' # computed fp^(0) function'
            WRITE (log_unit, '(A)') ' # i, w(a.u.), fp'
            DO j = 1, n
               WRITE (log_unit, *) j, j*dw, j*0.5_dp*hbokt*dw, fp(j)
            END DO
         END IF
         ! P > 1: QTB-PIMD
      ELSE
         !**** initialization ****
         dx1 = 0.5_dp*hbokt*dw
         xmin = 1.0e-7_dp !these values allows for an acceptable
         dx = 0.05_dp !ratio between accuracy, computing time and
         xmax = 10000.0_dp !memory requirement - tested for P up to 1024
         nx = INT((xmax - xmin)/dx) + 1
         nx = nx + nx/5 !add 20% points to avoid any problems at the end
         !of the interval (probably unnecessary)
         IF (ibead == 1) THEN
            op = 1.0_dp/p
            malpha = op !mixing parameter alpha = 1/P
            niter = 30 !30 iterations are enough to converge

            IF (print_level == debug_print_level) THEN
               WRITE (log_unit, '(A)') ' # ------------------------------------------------'
               WRITE (log_unit, '(A)') ' # computing fp^(0) function'
               WRITE (log_unit, '(A)') ' # parameters used:'
               WRITE (log_unit, '(A,ES13.3)') ' # dx = ', dx
               WRITE (log_unit, '(A,ES13.3)') ' # xmin = ', xmin
               WRITE (log_unit, '(A,ES13.3)') ' # xmax = ', xmax
               WRITE (log_unit, '(A,I8,I8)') ' # nx, n = ', nx, n
            END IF

            ALLOCATE (x(nx))
            ALLOCATE (x2(nx))
            ALLOCATE (h(nx))
            ALLOCATE (fp1(nx))
            ALLOCATE (xk(p - 1, nx))
            ALLOCATE (xk2(p - 1, nx))
            ALLOCATE (kk(p - 1, nx))
            ALLOCATE (fpxk(p - 1, nx))

            ! initialize fp(x)
            ! fp1 = fp(x) = h(x/P)
            ! fpxk = fp(xk) = h(xk/P)
            DO j = 1, nx
               x(j) = xmin + (j - 1)*dx
               x2(j) = x(j)**2
               h(j) = x(j)/TANH(x(j))
               IF (x(j) <= 1.0e-10_dp) h(j) = 1.0_dp
               fp1(j) = op*x(j)/TANH(x(j)*op)
               IF (x(j)*op <= 1.0e-10_dp) fp1(j) = 1.0_dp
               DO k = 1, p - 1
                  xk2(k, j) = x2(j) + (p*SIN(k*pi*op))**2
                  xk(k, j) = SQRT(xk2(k, j))
                  kk(k, j) = NINT((xk(k, j) - xmin)/dx) + 1
                  fpxk(k, j) = xk(k, j)*op/TANH(xk(k, j)*op)
                  IF (xk(k, j)*op <= 1.0e-10_dp) fpxk(k, j) = 1.0_dp
               END DO
            END DO

            ! **** resolution ****
            ! compute fp(x)
            DO i = 1, niter
               err = 0.0_dp
               DO j = 1, nx
                  tmp = 0.0_dp
                  DO k = 1, p - 1
                     tmp = tmp + fpxk(k, j)*x2(j)/xk2(k, j)
                  END DO
                  fprev = fp1(j)
                  fp1(j) = malpha*(h(j) - tmp) + (1.0_dp - malpha)*fp1(j)
                  IF (j <= n) err = err + ABS(1.0_dp - fp1(j)/fprev) ! compute "errors"
               END DO
               err = err/n

               ! Linear regression on the last 20% of the F_P function
               CALL pint_qtb_linreg(fp1(8*nx/10:nx), x(8*nx/10:nx), aa, bb, r2, log_unit, print_level)

               ! compute the new F_P(xk*sqrt(P))
               ! through linear interpolation
               ! or linear extrapolation if outside of the range
               DO j = 1, nx
                  DO k = 1, p - 1
                     IF (kk(k, j) < nx) THEN
                        fpxk(k, j) = fp1(kk(k, j)) + (fp1(kk(k, j) + 1) - fp1(kk(k, j)))/dx* &
                                     (xk(k, j) - x(kk(k, j)))
                     ELSE
                        fpxk(k, j) = aa*xk(k, j) + bb
                     ENDIF
                  END DO
               END DO
            END DO

            IF (print_level == debug_print_level) THEN
               ! **** tests ****
               WRITE (log_unit, '(A,ES9.3)') ' # average error during computation: ', err
               WRITE (log_unit, '(A,ES9.3)') ' # slope of F_P at high freq. - theoretical: ', op
               WRITE (log_unit, '(A,ES9.3)') ' # slope of F_P at high freq. - calculated: ', aa
               WRITE (log_unit, '(A,F6.3)') ' # F_P at zero freq. - theoretical: ', 1.0_dp
               WRITE (log_unit, '(A,F6.3)') ' # F_P at zero freq. - calculated: ', fp1(1)
            ELSE IF (print_level > silent_print_level) THEN
               CALL pint_write_line("QTB| Initialization of random forces using fP0 function")
               CALL pint_write_line("QTB| Computation of fP0 function")
               WRITE (line, '(A,ES9.3)') 'QTB| average error  ', err
               CALL pint_write_line(TRIM(line))
               WRITE (line, '(A,ES9.3)') 'QTB| slope at high frequency - theoretical: ', op
               CALL pint_write_line(TRIM(line))
               WRITE (line, '(A,ES9.3)') 'QTB| slope at high frequency - calculated:  ', aa
               CALL pint_write_line(TRIM(line))
               WRITE (line, '(A,F6.3)') 'QTB| value at zero frequency - theoretical:  ', 1.0_dp
               CALL pint_write_line(TRIM(line))
               WRITE (line, '(A,F6.3)') 'QTB| value at zero frequency - calculated:  ', fp1(1)
               CALL pint_write_line(TRIM(line))
            END IF

            IF (print_level == debug_print_level) THEN
               ! **** write solution ****
               WRITE (log_unit, '(A)') ' # ------------------------------------------------'
               WRITE (log_unit, '(A)') ' # computed fp function'
               WRITE (log_unit, '(A)') ' # i, w(a.u.), x, fp'
               DO j = 1, nx
                  WRITE (log_unit, *) j, j*dw, xmin + (j - 1)*dx, fp1(j)
               END DO
            END IF

            DEALLOCATE (x)
            DEALLOCATE (x2)
            DEALLOCATE (h)
            DEALLOCATE (xk)
            DEALLOCATE (xk2)
            DEALLOCATE (kk)
            DEALLOCATE (fpxk)
         END IF

         ! compute values of fP on the grid points for the current NM
         ! through linear interpolation / regression
         DO j = 1, n
            x1 = j*dx1
            k = NINT((x1 - xmin)/dx) + 1
            IF (k > nx) THEN
               fp(j) = aa*x1 + bb
            ELSE IF (k <= 0) THEN
               CALL pint_write_line("QTB| error in fp computation x < xmin")
               CPABORT("Error in fp computation (x < xmin) in initialization of QTB random forces")
            ELSE
               xx = xmin + (k - 1)*dx
               IF (x1 > xx) THEN
                  fp(j) = fp1(k) + (fp1(k + 1) - fp1(k))/dx*(x1 - xx)
               ELSE
                  fp(j) = fp1(k) + (fp1(k) - fp1(k - 1))/dx*(x1 - xx)
               END IF
            END IF
         END DO

      END IF

   END SUBROUTINE pint_qtb_computefp0

   ! ***************************************************************************
   !> \brief compute the f_P^(1) function necessary for coupling QTB with PIMD
   !> \param fp stores the computed function on the grid used for the generation
   !> of the filter h
   !> \param fp1 stores the computed function on an larger and finer grid
   !> \param dw angular frequency step
   !> \author Fabien Brieuc
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param fp ..
!> \param fp1 ..
!> \param dw ...
!> \param aa ...
!> \param bb ...
!> \param log_unit ...
!> \param ibead ...
!> \param print_level ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_computefp1(pint_env, fp, fp1, dw, aa, bb, log_unit, ibead, print_level)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(KIND=dp), DIMENSION(:)                        :: fp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: fp1
      REAL(KIND=dp)                                      :: dw, aa, bb
      INTEGER                                            :: log_unit, ibead, print_level

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_computefp1', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=200)                                 :: line
      INTEGER                                            :: i, j, k, n, niter, nx, p
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: kk
      REAL(KIND=dp)                                      :: dx, dx1, err, fprev, hbokT, malpha, op, &
                                                            op1, r2, tmp, tmp1, xmax, xmin, xx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: h, x, x2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fpxk, xk, xk2

      n = SIZE(fp)
      p = pint_env%p

      !using the physical temperature (T) not the simulation one (TxP)
      hbokT = 1.0_dp/(pint_env%kT*pint_env%propagator%temp_sim2phys)

      !Centroid NM (ibead=1) : classical
      !fp = 1
      IF (ibead == 1) THEN
         DO j = 1, n
            fp(j) = 1.0_dp
         END DO
      ELSE
         !**** initialization ****
         dx1 = 0.5_dp*hbokt*dw
         xmin = 1.0e-3_dp !these values allows for an acceptable
         dx = 0.05_dp !ratio between accuracy, computing time and
         xmax = 10000.0_dp !memory requirement - tested for P up to 1024
         nx = INT((xmax - xmin)/dx) + 1
         nx = nx + nx/5 !add 20% points to avoid problem at the end
         !of the interval (probably unnecessary)
         op = 1.0_dp/p
         IF (ibead == 2) THEN
            op1 = 1.0_dp/(p - 1)
            malpha = op !mixing parameter alpha = 1/P
            niter = 40 !40 iterations are enough to converge

            IF (print_level == debug_print_level) THEN
               ! **** write solution ****
               WRITE (log_unit, '(A)') ' # ------------------------------------------------'
               WRITE (log_unit, '(A)') ' # computing fp^(1) function'
               WRITE (log_unit, '(A)') ' # parameters used:'
               WRITE (log_unit, '(A,ES13.3)') ' # dx = ', dx
               WRITE (log_unit, '(A,ES13.3)') ' # xmin = ', xmin
               WRITE (log_unit, '(A,ES13.3)') ' # xmax = ', xmax
               WRITE (log_unit, '(A,I8,I8)') ' # nx, n = ', nx, n
            END IF

            ALLOCATE (x(nx))
            ALLOCATE (x2(nx))
            ALLOCATE (h(nx))
            ALLOCATE (fp1(nx))
            ALLOCATE (xk(p - 1, nx))
            ALLOCATE (xk2(p - 1, nx))
            ALLOCATE (kk(p - 1, nx))
            ALLOCATE (fpxk(p - 1, nx))

            ! initialize F_P(x) = f_P(x_1)
            ! fp1 = fp(x) = h(x/(P-1))
            ! fpxk = fp(xk) = h(xk/(P-1))
            DO j = 1, nx
               x(j) = xmin + (j - 1)*dx
               x2(j) = x(j)**2
               h(j) = x(j)/TANH(x(j))
               IF (x(j) <= 1.0e-10_dp) h(j) = 1.0_dp
               fp1(j) = op1*x(j)/TANH(x(j)*op1)
               IF (x(j)*op1 <= 1.0e-10_dp) fp1(j) = 1.0_dp
               DO k = 1, p - 1
                  xk2(k, j) = x2(j) + (p*SIN(k*pi*op))**2
                  xk(k, j) = SQRT(xk2(k, j) - (p*SIN(pi*op))**2)
                  kk(k, j) = NINT((xk(k, j) - xmin)/dx) + 1
                  fpxk(k, j) = xk(k, j)*op1/TANH(xk(k, j)*op1)
                  IF (xk(k, j)*op1 <= 1.0e-10_dp) fpxk(k, j) = 1.0_dp
               END DO
            END DO

            ! **** resolution ****
            ! compute fp(x)
            DO i = 1, niter
               err = 0.0_dp
               DO j = 1, nx
                  tmp = 0.0_dp
                  DO k = 2, p - 1
                     tmp = tmp + fpxk(k, j)*x2(j)/xk2(k, j)
                  END DO
                  fprev = fp1(j)
                  tmp1 = 1.0_dp + (p*SIN(pi*op)/x(j))**2
                  fp1(j) = malpha*tmp1*(h(j) - 1.0_dp - tmp) + (1.0_dp - malpha)*fp1(j)
                  IF (j <= n) err = err + ABS(1.0_dp - fp1(j)/fprev) ! compute "errors"
               END DO
               err = err/n

               ! Linear regression on the last 20% of the F_P function
               CALL pint_qtb_linreg(fp1(8*nx/10:nx), x(8*nx/10:nx), aa, bb, r2, log_unit, print_level)

               ! compute the new F_P(xk*sqrt(P))
               ! through linear interpolation
               ! or linear extrapolation if outside of the range
               DO j = 1, nx
                  DO k = 1, p - 1
                     IF (kk(k, j) < nx) THEN
                        fpxk(k, j) = fp1(kk(k, j)) + (fp1(kk(k, j) + 1) - fp1(kk(k, j)))/dx* &
                                     (xk(k, j) - x(kk(k, j)))
                     ELSE
                        fpxk(k, j) = aa*xk(k, j) + bb
                     END IF
                  END DO
               END DO
            END DO

            IF (print_level == debug_print_level) THEN
               ! **** tests ****
               WRITE (log_unit, '(A,ES9.3)') ' # average error during computation: ', err
               WRITE (log_unit, '(A,ES9.3)') ' # slope of F_P at high freq. - theoretical: ', op1
               WRITE (log_unit, '(A,ES9.3)') ' # slope of F_P at high freq. - calculated: ', aa
            ELSE IF (print_level > silent_print_level) THEN
               CALL pint_write_line("QTB| Initialization of random forces using fP1 function")
               CALL pint_write_line("QTB| Computation of fP1 function")
               WRITE (line, '(A,ES9.3)') 'QTB| average error  ', err
               CALL pint_write_line(TRIM(line))
               WRITE (line, '(A,ES9.3)') 'QTB| slope at high frequency - theoretical: ', op1
               CALL pint_write_line(TRIM(line))
               WRITE (line, '(A,ES9.3)') 'QTB| slope at high frequency - calculated:  ', aa
               CALL pint_write_line(TRIM(line))
            END IF

            IF (print_level == debug_print_level) THEN
               ! **** write solution ****
               WRITE (log_unit, '(A)') ' # ------------------------------------------------'
               WRITE (log_unit, '(A)') ' # computed fp function'
               WRITE (log_unit, '(A)') ' # i, w(a.u.), x, fp'
               DO j = 1, nx
                  WRITE (log_unit, *) j, j*dw, xmin + (j - 1)*dx, fp1(j)
               END DO
            END IF

            DEALLOCATE (x2)
            DEALLOCATE (h)
            DEALLOCATE (xk)
            DEALLOCATE (xk2)
            DEALLOCATE (kk)
            DEALLOCATE (fpxk)
         END IF

         ! compute values of fP on the grid points for the current NM
         ! trough linear interpolation / regression
         DO j = 1, n
            tmp = (j*dx1)**2 - (p*SIN(pi*op))**2
            IF (tmp < 0.d0) THEN
               fp(j) = fp1(1)
            ELSE
               tmp = SQRT(tmp)
               k = NINT((tmp - xmin)/dx) + 1
               IF (k > nx) THEN
                  fp(j) = aa*tmp + bb
               ELSE IF (k <= 0) THEN
                  fp(j) = fp1(1)
               ELSE
                  xx = xmin + (k - 1)*dx
                  IF (tmp > xx) THEN
                     fp(j) = fp1(k) + (fp1(k + 1) - fp1(k))/dx*(tmp - xx)
                  ELSE
                     fp(j) = fp1(k) + (fp1(k) - fp1(k - 1))/dx*(tmp - xx)
                  END IF
               END IF
            END IF
         END DO

      END IF

   END SUBROUTINE pint_qtb_computefp1

   ! ***************************************************************************
   !> \brief perform a simple linear regression - y(x) = a*x + b
   !> \author Fabien Brieuc
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param y ...
!> \param x ...
!> \param a ...
!> \param b ...
!> \param r2 ...
!> \param log_unit ...
!> \param print_level ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_linreg(y, x, a, b, r2, log_unit, print_level)
      REAL(KIND=dp), DIMENSION(:)                        :: y, x
      REAL(KIND=dp)                                      :: a, b, r2
      INTEGER                                            :: log_unit, print_level

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_linreg', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=200)                                 :: line
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: xav, xvar, xycov, yav, yvar

      n = SIZE(y)

      xav = 0.0_dp
      yav = 0.0_dp
      xycov = 0.0_dp
      xvar = 0.0_dp
      yvar = 0.0_dp

      DO i = 1, n
         xav = xav + x(i)
         yav = yav + y(i)
         xycov = xycov + x(i)*y(i)
         xvar = xvar + x(i)**2
         yvar = yvar + y(i)**2
      END DO

      xav = xav/n
      yav = yav/n
      xycov = xycov/n
      xycov = xycov - xav*yav
      xvar = xvar/n
      xvar = xvar - xav**2
      yvar = yvar/n
      yvar = yvar - yav**2

      a = xycov/xvar
      b = yav - a*xav

      r2 = xycov/SQRT(xvar*yvar)

      IF (r2 < 0.9_dp) THEN
         IF (print_level == debug_print_level) THEN
            WRITE (log_unit, '(A, E10.3)') '# possible error during linear regression: r^2 = ', r2
         ELSE IF (print_level > silent_print_level) THEN
            WRITE (line, '(A,E10.3)') 'QTB| possible error during linear regression: r^2 = ', r2
            CALL pint_write_line(TRIM(line))
         END IF
      END IF

   END SUBROUTINE pint_qtb_linreg

! **************************************************************************************************
!> \brief ...
!> \param z_in ...
!> \param z_out ...
!> \param plan ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE pint_qtb_fft(z_in, z_out, plan, n)

      INTEGER                                            :: n
      TYPE(fft_plan_type)                                :: plan
      COMPLEX(KIND=dp), DIMENSION(n)                     :: z_out, z_in

      CHARACTER(len=*), PARAMETER :: routineN = 'pint_qtb_fft', routineP = moduleN//':'//routineN

      INTEGER                                            :: stat

      CALL fft_create_plan_1dm(plan, fft_type, FWFFT, .FALSE., n, 1, z_in, z_out, fft_plan_style)
      CALL fft_1dm(plan, z_in, z_out, 1.d0, stat)
      CALL fft_destroy_plan(plan)
   END SUBROUTINE pint_qtb_fft

END MODULE
