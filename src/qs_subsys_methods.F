!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines that work on qs_subsys_type
!> \author Ole Schuett
! *****************************************************************************
MODULE qs_subsys_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_methods,                    ONLY: read_cell,&
                                             write_cell
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: cp_subsys_create
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE molecule_kind_types,             ONLY: num_ao_el_per_molecule
  USE qs_kind_types,                   ONLY: create_qs_kind_set,&
                                             qs_kind_type
  USE qs_subsys_types,                 ONLY: qs_subsys_set,&
                                             qs_subsys_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_subsys_methods'

  PUBLIC :: qs_subsys_create


CONTAINS

! *****************************************************************************
!> \brief Creates a qs_subsys. Optionally an existsing cp_subsys is used.
!> \param subsys ...
!> \param para_env ...
!> \param root_section ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param cp_subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_subsys_create(subsys, para_env, root_section,force_env_section, subsys_section, &
                              use_motion_section, cp_subsys, cell, cell_ref, error)
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section, &
                                                force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: cp_subsys
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_subsys_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: use_ref_cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: my_cell, my_cell_ref
    TYPE(cp_subsys_type), POINTER            :: my_cp_subsys
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: cell_section, kind_section

    NULLIFY(atomic_kind_set,qs_kind_set,cell_section,kind_section,my_cell,my_cell_ref,my_cp_subsys)

    IF(ASSOCIATED(subsys)) STOP "qs_subsys_create: subsys already associated"

    ! create cp_subsys
    IF (PRESENT(cp_subsys)) THEN
       my_cp_subsys => cp_subsys
    ELSE
       CALL cp_subsys_create(my_cp_subsys, para_env, root_section=root_section,&
                             force_env_section=force_env_section,&
                             subsys_section=subsys_section,&
                             use_motion_section=use_motion_section,&
                             error=error)
    END IF

    ! create cp_subsys%cell
    !TODO: moved to cp_subsys_create(), needs further disentanglement of cell_ref.
    IF (PRESENT(cell)) THEN
       my_cell => cell
       IF (PRESENT(cell_ref)) THEN
          my_cell_ref => cell_ref
          use_ref_cell = .TRUE.
       ELSE
          CALL cell_create(my_cell_ref, error=error)
          CALL cell_clone(my_cell, my_cell_ref, error)
          use_ref_cell = .FALSE.
       END IF
    ELSE
       cell_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
       CALL read_cell(my_cell, my_cell_ref, use_ref_cell=use_ref_cell, &
            cell_section=cell_section, para_env=para_env, error=error)
    END IF
    CALL cp_subsys_set(my_cp_subsys, cell=my_cell, error=error)
    CALL write_cell(my_cell,subsys_section,cell_ref=my_cell_ref,error=error)

    ! setup qs_kinds
    CALL cp_subsys_get(my_cp_subsys, atomic_kind_set=atomic_kind_set, error=error)
    kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)
    CALL create_qs_kind_set(qs_kind_set, atomic_kind_set, kind_section, &
                            para_env, force_env_section, error)

    CALL num_ao_el_per_molecule(my_cp_subsys%molecule_kinds_new%els, qs_kind_set)

    ! finally create qs_subsys
    ALLOCATE(subsys)
    CALL qs_subsys_set(subsys,&
                      cp_subsys=my_cp_subsys,&
                      cell_ref=my_cell_ref,&
                      use_ref_cell=use_ref_cell,&
                      qs_kind_set=qs_kind_set,&
                      error=error)

    IF (.NOT.PRESENT(cell))      CALL cell_release(my_cell,error=error)
    IF (.NOT.PRESENT(cell_ref))  CALL cell_release(my_cell_ref,error=error)
    IF (.NOT.PRESENT(cp_subsys)) CALL cp_subsys_release(my_cp_subsys,error=error)
  END SUBROUTINE qs_subsys_create
END MODULE qs_subsys_methods
