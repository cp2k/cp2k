!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A partial duplicate of f77_int_low.F for internal usage within CP2K.
!> \author Ole Schuett
! *****************************************************************************

! *****************************************************************************
SUBROUTINE cp_create_fenv_comm_int(new_env_id,input_file_path,input_declaration_cptr,&
   output_file_path,mpi_comm,ierr)
  USE f77_interface,                   ONLY: create_force_env
  USE input_section_types,             ONLY: section_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE iso_c_binding
  IMPLICIT NONE
  CHARACTER(len=*)                     :: input_file_path, output_file_path
  INTEGER                              :: new_env_id, mpi_comm, ierr
  TYPE(cp_error_type)                  :: error
  TYPE(c_ptr)                          :: input_declaration_cptr
  TYPE(section_type), POINTER          :: input_declaration

  NULLIFY(input_declaration)
  CALL c_f_pointer(input_declaration_cptr, input_declaration)
  CALL create_force_env(new_env_id,input_declaration,input_file_path,output_file_path,mpi_comm,ierr=ierr)
END SUBROUTINE cp_create_fenv_comm_int

! *****************************************************************************
RECURSIVE SUBROUTINE cp_destroy_fenv_int(env_id, ierr)
  USE f77_interface, ONLY: dfe => destroy_force_env
  IMPLICIT NONE
  INTEGER :: env_id, ierr

  CALL dfe(env_id,ierr)
END SUBROUTINE cp_destroy_fenv_int

! *****************************************************************************
SUBROUTINE cp_set_pos_int(env_id, new_pos, n_el, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: sp => set_pos
  IMPLICIT NONE
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: new_pos

  CALL sp(env_id,new_pos,n_el,ierr)
END SUBROUTINE cp_set_pos_int

! *****************************************************************************
SUBROUTINE cp_set_vel_int(env_id, new_vel, n_el, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: sv => set_vel
  IMPLICIT NONE
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: new_vel

  CALL sv(env_id,new_vel,n_el,ierr)
END SUBROUTINE cp_set_vel_int

! *****************************************************************************
SUBROUTINE cp_get_natom_int(env_id, natom, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: gna => get_natom
  IMPLICIT NONE
  INTEGER :: natom, env_id, ierr

  CALL gna(env_id,natom,ierr)
END SUBROUTINE cp_get_natom_int

! *****************************************************************************
SUBROUTINE cp_get_nparticle_int(env_id, nparticle, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: gnp => get_nparticle
  IMPLICIT NONE
  INTEGER :: nparticle, env_id, ierr

  CALL gnp(env_id,nparticle,ierr)
END SUBROUTINE cp_get_nparticle_int

! *****************************************************************************
SUBROUTINE cp_get_pos_int(env_id, pos, n_el, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: gp => get_pos
  IMPLICIT NONE
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: pos

  CALL gp(env_id,pos,n_el,ierr)
END SUBROUTINE cp_get_pos_int

! *****************************************************************************
SUBROUTINE cp_calc_force_int(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: cf => calc_force
  IMPLICIT NONE
  INTEGER :: env_id, ierr, n_el_pos, n_el_force
  REAL(dp) :: e_pot
  REAL(dp), DIMENSION(1:n_el_pos) :: pos
  REAL(dp), DIMENSION(1:n_el_force) :: force

  CALL cf(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)
END SUBROUTINE cp_calc_force_int

! *****************************************************************************
RECURSIVE SUBROUTINE cp_rep_init(rep_env_id,ierr)
  USE replica_methods, ONLY: rpin => rep_env_init_low
  IMPLICIT NONE
  INTEGER :: rep_env_id,ierr

  CALL rpin(rep_env_id,ierr)
END SUBROUTINE cp_rep_init

! *****************************************************************************
RECURSIVE SUBROUTINE cp_rep_destroy(rep_env_id,ierr)
  USE replica_methods, ONLY: rpin => rep_env_destroy_low
  IMPLICIT NONE
  INTEGER :: rep_env_id,ierr

  CALL rpin(rep_env_id,ierr)
END SUBROUTINE cp_rep_destroy

! *****************************************************************************
RECURSIVE SUBROUTINE cp_rep_calc_e_f(rep_env_id,calc_f,ierr)
  USE replica_methods, ONLY: rpef => rep_env_calc_e_f_low
  IMPLICIT NONE
  INTEGER :: rep_env_id,calc_f,ierr

  CALL rpef(rep_env_id,calc_f,ierr)
END SUBROUTINE cp_rep_calc_e_f

! *****************************************************************************
RECURSIVE SUBROUTINE cp_ep_init(ep_env_id,input_declaration_cptr,ierr)
  USE ep_methods,                      ONLY: ep_env_init_low
  USE input_section_types,             ONLY: section_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE iso_c_binding
  IMPLICIT NONE
  INTEGER :: ep_env_id, ierr
  TYPE(cp_error_type)                  :: error
  TYPE(c_ptr)                          :: input_declaration_cptr
  TYPE(section_type), POINTER          :: input_declaration
  NULLIFY(input_declaration)
  CALL c_f_pointer(input_declaration_cptr, input_declaration)
  CALL ep_env_init_low(ep_env_id, input_declaration, ierr)
END SUBROUTINE cp_ep_init

! *****************************************************************************
RECURSIVE SUBROUTINE cp_ep_calc_e_f(ep_env_id,calc_f,ierr)
  USE ep_methods, ONLY: epef => ep_env_calc_e_f_low
  IMPLICIT NONE
  INTEGER :: ep_env_id, calc_f, ierr

  CALL epef(ep_env_id, calc_f, ierr)
END SUBROUTINE cp_ep_calc_e_f

