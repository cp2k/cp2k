!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE input_restart_rng
  USE cp_linked_list_val,              ONLY: cp_sll_val_create,&
                                             cp_sll_val_get_length,&
                                             cp_sll_val_type
  USE input_section_types,             ONLY: section_get_keyword_index,&
                                             section_type,&
                                             section_vals_add_values,&
                                             section_vals_type
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type
  USE parallel_rng_types,              ONLY: rng_record_length
  USE string_utilities,                ONLY: ascii_to_string
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_restart_rng'

  PUBLIC :: section_rng_val_set

CONTAINS

! *****************************************************************************
!> \brief routine to dump rngs.. fast implementation
!> \param rng_section ...
!> \param nsize ...
!> \param ascii ...
!> \par History
!>      02.2006 created [teo]
!>      - string dump (again) instead of integer ASCII code (07.03.06,MK)
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_rng_val_set(rng_section, nsize, ascii)

    TYPE(section_vals_type), POINTER         :: rng_section
    INTEGER, INTENT(IN)                      :: nsize
    INTEGER, DIMENSION(:, :)                 :: ascii

    CHARACTER(LEN=*), PARAMETER :: routineN = 'section_rng_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: ik, irk, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure = .FALSE.

    CPASSERT(ASSOCIATED(rng_section))
    CPASSERT(rng_section%ref_count > 0)

    NULLIFY (my_val,old_val,section,vals)

    section => rng_section%section

    ik = section_get_keyword_index(section,"_DEFAULT_KEYWORD_")

    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
                   "section "//TRIM(section%name)//" does not contain keyword "//&
                   "_DEFAULT_KEYWORD_",failure)

    DO
      IF (SIZE(rng_section%values,2)==1) EXIT
      CALL section_vals_add_values(rng_section)
    END DO

    vals => rng_section%values(ik,1)%list
    Nlist = 0

    IF (ASSOCIATED(vals)) THEN
      Nlist = cp_sll_val_get_length(vals)
    END IF

    DO irk=1,nsize

      CALL ascii_to_string(ascii(:,irk),rng_record)
      CALL val_create(val=my_val,lc_val=rng_record)

      IF (Nlist /= 0) THEN
        IF (irk == 1) THEN
          new_pos => vals
        ELSE
          new_pos => new_pos%rest
        END IF
        old_val => new_pos%first_el
        CALL val_release(old_val)
        new_pos%first_el => my_val
      ELSE
        IF (irk == 1) THEN
          NULLIFY (new_pos)
          CALL cp_sll_val_create(new_pos,first_el=my_val)
          vals => new_pos
        ELSE
          NULLIFY (new_pos%rest)
          CALL cp_sll_val_create(new_pos%rest,first_el=my_val)
          new_pos => new_pos%rest
        END IF
      END IF

      NULLIFY (my_val)

    END DO

    rng_section%values(ik,1)%list => vals

  END SUBROUTINE section_rng_val_set

END MODULE input_restart_rng
