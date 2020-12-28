!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The type definitions for the PWDFT environment
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
MODULE pwdft_environment_types
   USE ISO_C_BINDING, ONLY: C_NULL_PTR, &
                            C_PTR
   USE cp_para_types, ONLY: cp_para_env_type
   USE cp_subsys_types, ONLY: cp_subsys_type
   USE input_section_types, ONLY: section_vals_type
   USE kinds, ONLY: dp
   USE qs_subsys_types, ONLY: qs_subsys_get, &
                              qs_subsys_release, &
                              qs_subsys_retain, &
                              qs_subsys_set, &
                              qs_subsys_type

#if defined(__SIRIUS)
   USE sirius, ONLY: sirius_free_handler
#endif

#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pwdft_environment_types'
   INTEGER, PRIVATE, SAVE :: last_pwdft_id = 0

   ! *** Public data types ***
   PUBLIC :: pwdft_environment_type, pwdft_energy_type

   ! *** Public subroutines ***
   PUBLIC :: pwdft_env_retain, &
             pwdft_env_release, &
             pwdft_env_set, &
             pwdft_env_get, &
             pwdft_env_create

! **************************************************************************************************
!> \brief The PWDFT energy type
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   TYPE pwdft_energy_type
      REAL(KIND=dp)                                    :: etotal
   END TYPE pwdft_energy_type

! **************************************************************************************************
!> \brief The PWDFT environment type
!> \param id_nr unique identifier
!> \param ref_count reference counter (see doc/ReferenceCounting.html)
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   TYPE pwdft_environment_type
      INTEGER                                          :: id_nr, &
                                                          ref_count
      TYPE(cp_para_env_type), POINTER                  :: para_env => NULL()
      TYPE(qs_subsys_type), POINTER                    :: qs_subsys => NULL()
      TYPE(section_vals_type), POINTER                 :: pwdft_input => NULL()
      TYPE(section_vals_type), POINTER                 :: force_env_input => NULL()
      TYPE(section_vals_type), POINTER                 :: xc_input => NULL()
      TYPE(pwdft_energy_type), POINTER                 :: energy => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER          :: forces => NULL()
      REAL(KIND=dp), DIMENSION(3, 3)                   :: stress
!     16 different functionals should be enough
      CHARACTER(len=80)                                :: xc_func(16)
      TYPE(C_PTR)                                      :: sctx = C_NULL_PTR
      TYPE(C_PTR)                                      :: gs_handler = C_NULL_PTR
      TYPE(C_PTR)                                      :: ks_handler = C_NULL_PTR
   END TYPE pwdft_environment_type

! **************************************************************************************************
!> \brief Type which allows the use of arrays of pointers
!> \param pwdft_env_ptr Pointer to the actual pwdft environment
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   TYPE pwdft_environment_p_type
      TYPE(pwdft_environment_type), POINTER :: pwdft_env_ptr
   END TYPE pwdft_environment_p_type

CONTAINS

! **************************************************************************************************
!> \brief Retains a pwdft environment (see doc/ReferenceCounting.html)
!> \param pwdft_env The pwdft environment to retain
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_env_retain(pwdft_env)

      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

!   ------------------------------------------------------------------------

      CPASSERT(ASSOCIATED(pwdft_env))
      CPASSERT(pwdft_env%ref_count > 0)
      pwdft_env%ref_count = pwdft_env%ref_count + 1
   END SUBROUTINE pwdft_env_retain

! **************************************************************************************************
!> \brief Releases the given pwdft environment (see doc/ReferenceCounting.html)
!> \param pwdft_env The pwdft environment to release
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_env_release(pwdft_env)

      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

!   ------------------------------------------------------------------------

      IF (ASSOCIATED(pwdft_env)) THEN
         CPASSERT(pwdft_env%ref_count > 0)
         pwdft_env%ref_count = pwdft_env%ref_count - 1
         IF (pwdft_env%ref_count < 1) THEN

#if defined(__SIRIUS)
            CALL sirius_free_handler(pwdft_env%gs_handler)
            CALL sirius_free_handler(pwdft_env%ks_handler)
            CALL sirius_free_handler(pwdft_env%sctx)
#endif

            IF (ASSOCIATED(pwdft_env%qs_subsys)) THEN
               CALL qs_subsys_release(pwdft_env%qs_subsys)
            END IF
            IF (ASSOCIATED(pwdft_env%energy)) THEN
               DEALLOCATE (pwdft_env%energy)
            END IF
            IF (ASSOCIATED(pwdft_env%forces)) THEN
               DEALLOCATE (pwdft_env%forces)
            END IF
            !
            DEALLOCATE (pwdft_env)
         END IF
      END IF
      NULLIFY (pwdft_env)
   END SUBROUTINE pwdft_env_release

! **************************************************************************************************
!> \brief Returns various attributes of the pwdft environment
!> \param pwdft_env The enquired pwdft environment
!> \param id_nr unique identifier
!> \param pwdft_input ...
!> \param force_env_input ...
!> \param xc_input ...
!> \param cp_subsys ...
!> \param qs_subsys ...
!> \param para_env ...
!> \param energy ...
!> \param forces ...
!> \param stress ...
!> \param sctx ...
!> \param gs_handler ...
!> \param ks_handler ...
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_env_get(pwdft_env, id_nr, pwdft_input, force_env_input, xc_input, &
                            cp_subsys, qs_subsys, para_env, energy, forces, stress, &
                            sctx, gs_handler, ks_handler)

      TYPE(pwdft_environment_type), POINTER              :: pwdft_env
      INTEGER, INTENT(OUT), OPTIONAL                     :: id_nr
      TYPE(section_vals_type), OPTIONAL, POINTER         :: pwdft_input, force_env_input, xc_input
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: cp_subsys
      TYPE(qs_subsys_type), OPTIONAL, POINTER            :: qs_subsys
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(pwdft_energy_type), OPTIONAL, POINTER         :: energy
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: forces
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL           :: stress
      TYPE(C_PTR), OPTIONAL                              :: sctx, gs_handler, ks_handler

!   ------------------------------------------------------------------------

      CPASSERT(ASSOCIATED(pwdft_env))
      CPASSERT(pwdft_env%ref_count > 0)

      IF (PRESENT(id_nr)) id_nr = pwdft_env%id_nr
      IF (PRESENT(pwdft_input)) pwdft_input => pwdft_env%pwdft_input
      IF (PRESENT(force_env_input)) force_env_input => pwdft_env%force_env_input
      IF (PRESENT(xc_input)) xc_input => pwdft_env%xc_input
      IF (PRESENT(qs_subsys)) qs_subsys => pwdft_env%qs_subsys
      IF (PRESENT(cp_subsys)) THEN
         CALL qs_subsys_get(pwdft_env%qs_subsys, cp_subsys=cp_subsys)
      END IF
      IF (PRESENT(para_env)) para_env => pwdft_env%para_env
      IF (PRESENT(energy)) energy => pwdft_env%energy
      IF (PRESENT(forces)) forces => pwdft_env%forces
      IF (PRESENT(stress)) stress(1:3, 1:3) = pwdft_env%stress(1:3, 1:3)
      IF (PRESENT(sctx)) sctx = pwdft_env%sctx
      IF (PRESENT(gs_handler)) gs_handler = pwdft_env%gs_handler
      IF (PRESENT(ks_handler)) ks_handler = pwdft_env%ks_handler

   END SUBROUTINE pwdft_env_get

! **************************************************************************************************
!> \brief Sets various attributes of the pwdft environment
!> \param pwdft_env The enquired pwdft environment
!> \param id_nr unique identifier
!> \param pwdft_input ...
!> \param force_env_input ...
!> \param xc_input ...
!> \param qs_subsys ...
!> \param cp_subsys ...
!> \param para_env ...
!> \param energy ...
!> \param forces ...
!> \param stress ...
!> \param sctx ...
!> \param gs_handler ...
!> \param ks_handler ...
!> \par History
!>      07.2018 initial create
!> \author JHU
!> \note
!>   For possible missing arguments see the attributes of pwdft_environment_type
! **************************************************************************************************
   SUBROUTINE pwdft_env_set(pwdft_env, id_nr, pwdft_input, force_env_input, xc_input, &
                            qs_subsys, cp_subsys, para_env, energy, forces, stress, &
                            sctx, gs_handler, ks_handler)

      TYPE(pwdft_environment_type), POINTER              :: pwdft_env
      INTEGER, INTENT(IN), OPTIONAL                      :: id_nr
      TYPE(section_vals_type), OPTIONAL, POINTER         :: pwdft_input, force_env_input, xc_input
      TYPE(qs_subsys_type), OPTIONAL, POINTER            :: qs_subsys
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: cp_subsys
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(pwdft_energy_type), OPTIONAL, POINTER         :: energy
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: forces
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL           :: stress
      TYPE(C_PTR), OPTIONAL                              :: sctx, gs_handler, ks_handler

!   ------------------------------------------------------------------------

      CPASSERT(ASSOCIATED(pwdft_env))
      CPASSERT(pwdft_env%ref_count > 0)

      IF (PRESENT(id_nr)) pwdft_env%id_nr = id_nr
      IF (PRESENT(para_env)) pwdft_env%para_env => para_env
      IF (PRESENT(pwdft_input)) pwdft_env%pwdft_input => pwdft_input
      IF (PRESENT(force_env_input)) pwdft_env%force_env_input => force_env_input
      IF (PRESENT(xc_input)) pwdft_env%xc_input => xc_input

      IF (PRESENT(qs_subsys)) THEN
         CALL qs_subsys_retain(qs_subsys)
         CALL qs_subsys_release(pwdft_env%qs_subsys)
         pwdft_env%qs_subsys => qs_subsys
      END IF
      IF (PRESENT(cp_subsys)) THEN
         CALL qs_subsys_set(pwdft_env%qs_subsys, cp_subsys=cp_subsys)
      END IF
      IF (PRESENT(energy)) pwdft_env%energy => energy
      IF (PRESENT(forces)) pwdft_env%forces => forces
      IF (PRESENT(stress)) pwdft_env%stress(1:3, 1:3) = stress(1:3, 1:3)
      IF (PRESENT(sctx)) pwdft_env%sctx = sctx
      IF (PRESENT(gs_handler)) pwdft_env%gs_handler = gs_handler
      IF (PRESENT(ks_handler)) pwdft_env%ks_handler = ks_handler

   END SUBROUTINE pwdft_env_set

! **************************************************************************************************
!> \brief Reinitializes the pwdft environment
!> \param pwdft_env The pwdft environment to be reinitialized
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_env_clear(pwdft_env)

      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

!   ------------------------------------------------------------------------

      CPASSERT(ASSOCIATED(pwdft_env))
      NULLIFY (pwdft_env%para_env)
      NULLIFY (pwdft_env%pwdft_input)
      NULLIFY (pwdft_env%force_env_input)
      IF (ASSOCIATED(pwdft_env%qs_subsys)) THEN
         CALL qs_subsys_release(pwdft_env%qs_subsys)
      END IF
      IF (ASSOCIATED(pwdft_env%energy)) THEN
         DEALLOCATE (pwdft_env%energy)
      END IF
      IF (ASSOCIATED(pwdft_env%forces)) THEN
         DEALLOCATE (pwdft_env%forces)
         NULLIFY (pwdft_env%forces)
      END IF
      pwdft_env%stress = 0.0_dp
      pwdft_env%sctx = C_NULL_PTR
      pwdft_env%gs_handler = C_NULL_PTR
      pwdft_env%ks_handler = C_NULL_PTR

   END SUBROUTINE pwdft_env_clear

! **************************************************************************************************
!> \brief Creates the pwdft environment
!> \param pwdft_env The pwdft environment to be created
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_env_create(pwdft_env)

      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

      ALLOCATE (pwdft_env)

      pwdft_env%ref_count = 1
      last_pwdft_id = last_pwdft_id + 1
      pwdft_env%id_nr = last_pwdft_id
      CALL pwdft_env_clear(pwdft_env)

   END SUBROUTINE pwdft_env_create

END MODULE pwdft_environment_types
