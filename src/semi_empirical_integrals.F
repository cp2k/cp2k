!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of wrappers for semi-empirical analytical/numerical Integrals
!>        routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
!> \par History
!>         05.2008 Teodoro Laino [tlaino] - University of Zurich - In core integrals
! *****************************************************************************
MODULE semi_empirical_integrals

  
  USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                             hfx_get_mult_cache_elements
  USE input_constants,                 ONLY: do_se_IS_slater
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE memory_utilities,                ONLY: reallocate
  USE semi_empirical_int_ana,          ONLY: corecore_ana,&
                                             corecore_el_ana,&
                                             rotint_ana,&
                                             rotnuc_ana
  USE semi_empirical_int_gks,          ONLY: corecore_gks,&
                                             drotint_gks,&
                                             drotnuc_gks,&
                                             rotint_gks,&
                                             rotnuc_gks
  USE semi_empirical_int_num,          ONLY: corecore_el_num,&
                                             corecore_num,&
                                             dcorecore_el_num,&
                                             dcorecore_num,&
                                             drotint_num,&
                                             drotnuc_num,&
                                             rotint_num,&
                                             rotnuc_num
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE semi_empirical_types,            ONLY: se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_integrals'
  PUBLIC :: rotint, rotnuc, corecore, corecore_el, drotint, drotnuc, dcorecore,&
            dcorecore_el

CONTAINS

! *****************************************************************************
!> \brief  wrapper for numerical/analytical 2 center 2 electrons integrals
!>         routines with possibility of incore storage/compression
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param w ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \param error ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE rotint(sepi,sepj,rij,w,anag,se_int_control,se_taper,store_int_env,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(2025), INTENT(OUT)   :: w
    LOGICAL                                  :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: buffer_left, buffer_size, &
                                                buffer_start, cache_size, &
                                                memory_usage, nbits, &
                                                new_size, nints
    INTEGER(KIND=int_8)                      :: mem_compression_counter
    LOGICAL                                  :: buffer_overflow, failure
    REAL(KIND=dp)                            :: eps_storage

    failure = .FALSE.
    w(:) = 0.0_dp
    IF (.NOT.store_int_env%memory_parameter%do_all_on_the_fly) THEN
       nints      = (sepi%natorb*(sepi%natorb+1)/2)*(sepj%natorb*(sepj%natorb+1)/2)
       cache_size = store_int_env%memory_parameter%cache_size
       eps_storage= store_int_env%memory_parameter%eps_storage_scaling
       IF (store_int_env%filling_containers) THEN
          mem_compression_counter = store_int_env%memory_parameter%actual_memory_usage * cache_size
          IF(mem_compression_counter > store_int_env%memory_parameter%max_compression_counter) THEN
             buffer_overflow = .TRUE.
             store_int_env%memory_parameter%ram_counter = store_int_env%nbuffer
          ELSE
             store_int_env%nbuffer = store_int_env%nbuffer+1
             buffer_overflow = .FALSE.
          END IF
          ! Compute Integrals
          IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
            CALL rotint_gks (sepi,sepj,rij,w,se_int_control=se_int_control,error=error)
          ELSE
            IF (anag) THEN
               CALL rotint_ana (sepi,sepj,rij,w,se_int_control=se_int_control,se_taper=se_taper,error=error)
            ELSE
               CALL rotint_num (sepi,sepj,rij,w,se_int_control=se_int_control,se_taper=se_taper,error=error)
            END IF
          END IF
          ! Store integrals if we did not go overflow
          IF (.NOT.buffer_overflow) THEN
             IF (store_int_env%compress) THEN
                ! Store integrals in the containers
                IF (store_int_env%nbuffer>SIZE(store_int_env%max_val_buffer)) THEN
                   new_size = store_int_env%nbuffer+1000
                   CALL reallocate(store_int_env%max_val_buffer, 1, new_size)
                END IF
                store_int_env%max_val_buffer(store_int_env%nbuffer) = MAXVAL(ABS(w(1:nints)))

                nbits        = EXPONENT(store_int_env%max_val_buffer(store_int_env%nbuffer)/eps_storage) + 1
                buffer_left  = nints
                buffer_start = 1
                DO WHILE (buffer_left > 0)
                   buffer_size = MIN(buffer_left, cache_size)
                   CALL hfx_add_mult_cache_elements(w(buffer_start:),&
                                                    buffer_size, nbits, &
                                                    store_int_env%integral_caches(nbits), &
                                                    store_int_env%integral_containers(nbits), &
                                                    eps_storage, 1.0_dp,&
                                                    store_int_env%memory_parameter%actual_memory_usage, &
                                                    .FALSE., error)
                   buffer_left  = buffer_left  - buffer_size
                   buffer_start = buffer_start + buffer_size
                END DO
             ELSE
                ! Skip compression
                memory_usage = store_int_env%memory_parameter%actual_memory_usage
                CPPrecondition((nints/1.2_dp)<=HUGE(0)-memory_usage,cp_failure_level,routineP,error,failure)
                IF (memory_usage+nints>SIZE(store_int_env%uncompressed_container)) THEN
                   new_size = INT((memory_usage+nints)*1.2_dp)
                   CALL reallocate(store_int_env%uncompressed_container, 1, new_size)
                END IF
                store_int_env%uncompressed_container(memory_usage:memory_usage+nints-1)=w(1:nints)
                store_int_env%memory_parameter%actual_memory_usage = memory_usage + nints
             END IF
          END IF
       ELSE
          ! Get integrals from the containers
          IF(store_int_env%memory_parameter%ram_counter == store_int_env%nbuffer) THEN
             buffer_overflow = .TRUE.
          ELSE
             store_int_env%nbuffer = store_int_env%nbuffer + 1
             buffer_overflow = .FALSE.
          END IF
          ! Get integrals from cache unless we overflowed
          IF (.NOT.buffer_overflow) THEN
             IF (store_int_env%compress) THEN
                ! Get Integrals from containers
                nbits        = EXPONENT(store_int_env%max_val_buffer(store_int_env%nbuffer)/eps_storage) + 1
                buffer_left  = nints
                buffer_start = 1
                DO WHILE (buffer_left > 0)
                   buffer_size = MIN(buffer_left, cache_size)
                   CALL hfx_get_mult_cache_elements(w(buffer_start:), &
                                                    buffer_size, nbits, &
                                                    store_int_env%integral_caches(nbits), &
                                                    store_int_env%integral_containers(nbits), &
                                                    eps_storage, 1.0_dp, &
                                                    store_int_env%memory_parameter%actual_memory_usage, &
                                                    .FALSE.)
                   buffer_left  = buffer_left  - buffer_size
                   buffer_start = buffer_start + buffer_size
                END DO
             ELSE
                ! Skip compression
                memory_usage = store_int_env%memory_parameter%actual_memory_usage
                w(1:nints)   = store_int_env%uncompressed_container(memory_usage:memory_usage+nints-1)
                store_int_env%memory_parameter%actual_memory_usage = memory_usage + nints
             END IF
          ELSE
             IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
               CALL rotint_gks (sepi,sepj,rij,w,se_int_control=se_int_control,error=error)
             ELSE
                IF (anag) THEN
                   CALL rotint_ana (sepi,sepj,rij,w,se_int_control=se_int_control,se_taper=se_taper,error=error)
                ELSE
                   CALL rotint_num (sepi,sepj,rij,w,se_int_control=se_int_control,se_taper=se_taper,error=error)
                END IF
             END IF
          END IF
       END IF
    ELSE
       IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
          CALL rotint_gks (sepi,sepj,rij,w,se_int_control=se_int_control,error=error)
       ELSE
          IF (anag) THEN
             CALL rotint_ana (sepi,sepj,rij,w,se_int_control=se_int_control,se_taper=se_taper,error=error)
          ELSE
             CALL rotint_num (sepi,sepj,rij,w,se_int_control=se_int_control,se_taper=se_taper,error=error)
          END IF
       END IF
    END IF
  END SUBROUTINE rotint

! *****************************************************************************
!> \brief wrapper for numerical/analytical 1 center 1 electron integrals
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param e1b ...
!> \param e2a ...
!> \param itype ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \param error ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE rotnuc (sepi,sepj,rij,e1b,e2a,itype,anag,se_int_control,se_taper,store_int_env,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(45), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), OPTIONAL, &
      POINTER                                :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotnuc', &
      routineP = moduleN//':'//routineN

    INTEGER :: buffer_left, buffer_size, buffer_start, cache_size, &
      memory_usage, nbits, new_size, nints, nints_1, nints_2
    INTEGER(KIND=int_8)                      :: mem_compression_counter
    LOGICAL                                  :: buffer_overflow, &
                                                do_all_on_the_fly, failure
    REAL(KIND=dp)                            :: eps_storage, w(90)

    failure           = .FALSE.
    do_all_on_the_fly = .TRUE.
    IF (PRESENT(e1b)) e1b(:) = 0.0_dp
    IF (PRESENT(e2a)) e2a(:) = 0.0_dp
    IF (PRESENT(store_int_env)) do_all_on_the_fly = store_int_env%memory_parameter%do_all_on_the_fly
    IF (.NOT.do_all_on_the_fly) THEN
       nints_1    = 0
       nints_2    = 0
       IF (PRESENT(e1b)) nints_1 = (sepi%natorb*(sepi%natorb+1)/2)
       IF (PRESENT(e2a)) nints_2 = (sepj%natorb*(sepj%natorb+1)/2)
       nints                     = nints_1 + nints_2
       ! This is the upper limit for an spd basis set
       CPPrecondition(nints<=90,cp_failure_level,routineP,error,failure)
       cache_size = store_int_env%memory_parameter%cache_size
       eps_storage= store_int_env%memory_parameter%eps_storage_scaling
       IF (store_int_env%filling_containers) THEN
          mem_compression_counter = store_int_env%memory_parameter%actual_memory_usage * cache_size
          IF(mem_compression_counter > store_int_env%memory_parameter%max_compression_counter) THEN
             buffer_overflow = .TRUE.
             store_int_env%memory_parameter%ram_counter = store_int_env%nbuffer
          ELSE
             store_int_env%nbuffer = store_int_env%nbuffer+1
             buffer_overflow = .FALSE.
          END IF
          ! Compute Integrals
          IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
             CALL rotnuc_gks (sepi,sepj,rij,e1b=e1b,e2a=e2a,&
                  se_int_control=se_int_control,error=error)
          ELSE
             IF (anag) THEN
                CALL rotnuc_ana (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
                     se_int_control=se_int_control, se_taper=se_taper, error=error)
             ELSE
                CALL rotnuc_num (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
                     se_int_control=se_int_control, se_taper=se_taper, error=error)
             END IF
          END IF
          ! Store integrals if we did not go overflow
          IF (.NOT.buffer_overflow) THEN
             IF (PRESENT(e1b)) w(        1:nints_1) = e1b(1:nints_1)
             IF (PRESENT(e2a)) w(nints_1+1:nints  ) = e2a(1:nints_2)

             IF (store_int_env%compress) THEN
                ! Store integrals in the containers
                IF (store_int_env%nbuffer>SIZE(store_int_env%max_val_buffer)) THEN
                   new_size = store_int_env%nbuffer+1000
                   CALL reallocate(store_int_env%max_val_buffer, 1, new_size)
                END IF
                store_int_env%max_val_buffer(store_int_env%nbuffer) = MAXVAL(ABS(w(1:nints)))

                nbits        = EXPONENT(store_int_env%max_val_buffer(store_int_env%nbuffer)/eps_storage) + 1
                buffer_left  = nints
                buffer_start = 1
                DO WHILE (buffer_left > 0)
                   buffer_size = MIN(buffer_left, cache_size)
                   CALL hfx_add_mult_cache_elements(w(buffer_start:),&
                                                    buffer_size, nbits, &
                                                    store_int_env%integral_caches(nbits), &
                                                    store_int_env%integral_containers(nbits), &
                                                    eps_storage, 1.0_dp,&
                                                    store_int_env%memory_parameter%actual_memory_usage, &
                                                    .FALSE., error)
                   buffer_left  = buffer_left  - buffer_size
                   buffer_start = buffer_start + buffer_size
                END DO
             ELSE
                ! Skip compression
                memory_usage = store_int_env%memory_parameter%actual_memory_usage
                CPPrecondition((nints/1.2_dp)<=HUGE(0)-memory_usage,cp_failure_level,routineP,error,failure)
                IF (memory_usage+nints>SIZE(store_int_env%uncompressed_container)) THEN
                   new_size = INT((memory_usage+nints)*1.2_dp)
                   CALL reallocate(store_int_env%uncompressed_container, 1, new_size)
                END IF
                store_int_env%uncompressed_container(memory_usage:memory_usage+nints-1)=w(1:nints)
                store_int_env%memory_parameter%actual_memory_usage = memory_usage + nints
             END IF
          END IF
       ELSE
          ! Get integrals from the containers
          IF(store_int_env%memory_parameter%ram_counter == store_int_env%nbuffer) THEN
             buffer_overflow = .TRUE.
          ELSE
             store_int_env%nbuffer = store_int_env%nbuffer + 1
             buffer_overflow = .FALSE.
          END IF
          ! Get integrals from cache unless we overflowed
          IF (.NOT.buffer_overflow) THEN
             IF (store_int_env%compress) THEN
                ! Get Integrals from containers
                nbits        = EXPONENT(store_int_env%max_val_buffer(store_int_env%nbuffer)/eps_storage) + 1
                buffer_left  = nints
                buffer_start = 1
                DO WHILE (buffer_left > 0 )
                   buffer_size = MIN(buffer_left, cache_size)
                   CALL hfx_get_mult_cache_elements(w(buffer_start:), &
                                                    buffer_size, nbits, &
                                                    store_int_env%integral_caches(nbits), &
                                                    store_int_env%integral_containers(nbits), &
                                                    eps_storage, 1.0_dp, &
                                                    store_int_env%memory_parameter%actual_memory_usage, &
                                                    .FALSE.)
                   buffer_left  = buffer_left  - buffer_size
                   buffer_start = buffer_start + buffer_size
                END DO
             ELSE
                ! Skip compression
                memory_usage = store_int_env%memory_parameter%actual_memory_usage
                w(1:nints)   = store_int_env%uncompressed_container(memory_usage:memory_usage+nints-1)
                store_int_env%memory_parameter%actual_memory_usage = memory_usage + nints
             END IF
             IF (PRESENT(e1b)) e1b(1:nints_1) = w(        1:nints_1)
             IF (PRESENT(e2a)) e2a(1:nints_2) = w(nints_1+1:nints  )
          ELSE
             IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
                CALL rotnuc_gks (sepi,sepj,rij,e1b=e1b,e2a=e2a,&
                     se_int_control=se_int_control, error=error)
             ELSE
                IF (anag) THEN
                   CALL rotnuc_ana (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
                        se_int_control=se_int_control, se_taper=se_taper, error=error)
                ELSE
                   CALL rotnuc_num (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
                        se_int_control=se_int_control, se_taper=se_taper, error=error)
                END IF
             END IF
          END IF
       END IF
    ELSE
       IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
          CALL rotnuc_gks (sepi,sepj,rij,e1b=e1b,e2a=e2a,&
               se_int_control=se_int_control, error=error)
       ELSE
          IF (anag) THEN
             CALL rotnuc_ana (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
                  se_int_control=se_int_control, se_taper=se_taper, error=error)
          ELSE
             CALL rotnuc_num (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
                  se_int_control=se_int_control, se_taper=se_taper, error=error)
          END IF
       END IF
    END IF

  END SUBROUTINE rotnuc

! *****************************************************************************
!> \brief  wrapper for numerical/analytical routines
!>         core-core integrals, since are evaluated only once do not need to be
!>         stored.
!>
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param enuc ...
!> \param itype ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \date   04.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE corecore (sepi,sepj,rij,enuc,itype,anag,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore', &
      routineP = moduleN//':'//routineN

    enuc = 0.0_dp
    IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
       CALL corecore_gks (sepi,sepj,rij,enuc=enuc,se_int_control=se_int_control,error=error)
    ELSE
       IF (anag) THEN
          CALL corecore_ana (sepi,sepj,rij,enuc=enuc,itype=itype,se_int_control=se_int_control,&
               se_taper=se_taper, error=error)
       ELSE
          CALL corecore_num (sepi,sepj,rij,enuc=enuc,itype=itype,se_int_control=se_int_control,&
               se_taper=se_taper, error=error)
       END IF
    END IF

  END SUBROUTINE corecore

! *****************************************************************************
!> \brief  wrapper for numerical/analytical routines
!>         core-core electrostatic (only) integrals
!>
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param enuc ...
!> \param itype ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \date   05.2009
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE corecore_el (sepi,sepj,rij,enuc,itype,anag,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore_el', &
      routineP = moduleN//':'//routineN

    enuc = 0.0_dp
    IF (anag) THEN
       CALL corecore_el_ana (sepi,sepj,rij,enuc=enuc,itype=itype,se_int_control=se_int_control,&
            se_taper=se_taper, error=error)
    ELSE
       CALL corecore_el_num (sepi,sepj,rij,enuc=enuc,itype=itype,se_int_control=se_int_control,&
            se_taper=se_taper, error=error)
    END IF

  END SUBROUTINE corecore_el

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param dw ...
!> \param delta ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \date   04.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE drotint (sepi,sepj,rij,dw,delta,anag,se_int_control,se_taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3, 2025), &
      INTENT(OUT)                            :: dw
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotint', &
      routineP = moduleN//':'//routineN

    dw(:,:) = 0.0_dp
    IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
       CALL  drotint_gks(sepi,sepj,rij,dw=dw,se_int_control=se_int_control,error=error)
    ELSE
       IF (anag) THEN
          CALL  rotint_ana(sepi,sepj,rij,dw=dw,   se_int_control=se_int_control,se_taper=se_taper,error=error)
       ELSE
          CALL drotint_num(sepi,sepj,rij,dw,delta,se_int_control=se_int_control,se_taper=se_taper,error=error)
       END IF
    END IF

  END SUBROUTINE drotint

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param de1b ...
!> \param de2a ...
!> \param itype ...
!> \param delta ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \date   04.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE drotnuc (sepi,sepj,rij,de1b,de2a,itype,delta,anag,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3, 45), &
      INTENT(OUT), OPTIONAL                  :: de1b, de2a
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotnuc', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(de1b)) de1b(:,:) = 0.0_dp
    IF (PRESENT(de2a)) de2a(:,:) = 0.0_dp
    IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
          CALL drotnuc_gks (sepi, sepj, rij, de1b=de1b, de2a=de2a,&
               se_int_control=se_int_control,error=error)
    ELSE
       IF (anag) THEN
          CALL rotnuc_ana  (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype,&
               se_int_control=se_int_control,se_taper=se_taper,error=error)
       ELSE
          CALL drotnuc_num (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype,&
               delta=delta,se_int_control=se_int_control,se_taper=se_taper,error=error)
       END IF
    END IF

  END SUBROUTINE drotnuc

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param denuc ...
!> \param itype ...
!> \param delta ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \date   04.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE dcorecore (sepi,sepj,rij,denuc,itype,delta,anag,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dcorecore', &
      routineP = moduleN//':'//routineN

    denuc = 0.0_dp
    IF (se_int_control%integral_screening == do_se_IS_slater ) THEN
       CALL corecore_gks (sepi,sepj,rij,denuc=denuc,se_int_control=se_int_control,error=error)
    ELSE
       IF (anag) THEN
          CALL corecore_ana (sepi,sepj,rij,denuc=denuc,itype=itype,se_int_control=se_int_control,&
               se_taper=se_taper, error=error)
       ELSE
          CALL dcorecore_num (sepi,sepj,rij,denuc=denuc,delta=delta,itype=itype,&
               se_int_control=se_int_control,se_taper=se_taper, error=error)
       END IF
    END IF

  END SUBROUTINE dcorecore

! *****************************************************************************
!> \brief  wrapper for numerical/analytical routines
!>         core-core electrostatic (only) integrals derivatives
!>
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param denuc ...
!> \param itype ...
!> \param delta ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \date   05.2009
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE dcorecore_el (sepi,sepj,rij,denuc,itype,delta,anag,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dcorecore_el', &
      routineP = moduleN//':'//routineN

    denuc = 0.0_dp
    IF (anag) THEN
       CALL corecore_el_ana (sepi,sepj,rij,denuc=denuc,itype=itype,se_int_control=se_int_control,&
            se_taper=se_taper, error=error)
    ELSE
       CALL dcorecore_el_num (sepi,sepj,rij,denuc=denuc,delta=delta,itype=itype,&
            se_int_control=se_int_control,se_taper=se_taper, error=error)
    END IF

  END SUBROUTINE dcorecore_el

END MODULE semi_empirical_integrals
