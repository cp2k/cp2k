#include "gauxc/gauxc_config.f"

MODULE xc_gauxc

   USE iso_c_binding, ONLY: c_int, c_int64_t, c_double, c_bool, c_size_t, c_int32_t, c_null_ptr, c_ptr, c_f_pointer, c_associated, c_char
   USE kinds, ONLY: dp, default_string_length
   USE particle_types, ONLY: particle_type
   USE physcon, ONLY: bohr
   USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind, &
       & get_atomic_kind_set
   USE qs_integral_utils, ONLY: basis_set_list_setup
   USE basis_set_types, ONLY: gto_basis_set_type, gto_basis_set_p_type
   USE qs_kind_types, ONLY: qs_kind_type, get_qs_kind
   USE periodic_table, ONLY: get_ptable_info
   USE gauxc_status, ONLY: gauxc_status_type, gauxc_status_message
   USE gauxc_enums, ONLY: gauxc_executionspace, gauxc_radialquad, &
     & gauxc_atomicgridsizedefault, gauxc_pruningscheme
   USE gauxc_runtime_environment, ONLY: gauxc_runtime_environment_type, &
      & gauxc_runtime_environment_comm_rank, gauxc_runtime_environment_comm_size, &
      & gauxc_runtime_environment_new, gauxc_runtime_environment_delete
   USE gauxc_molecule, ONLY: gauxc_molecule_type, gauxc_molecule_new_from_atoms, gauxc_delete
   USE gauxc_atom, ONLY: gauxc_atom_type
   USE gauxc_basisset, ONLY: gauxc_basisset_type, gauxc_basisset_new, &
      & gauxc_basisset_new_from_shells, gauxc_delete
   USE gauxc_shell, ONLY: gauxc_shell_type
   USE gauxc_molgrid, ONLY: gauxc_molgrid_type, gauxc_molgrid_new_default, &
      & gauxc_delete
   USE gauxc_load_balancer, ONLY: gauxc_load_balancer_type, gauxc_load_balancer_factory_type, &
      & gauxc_load_balancer_factory_new, gauxc_load_balancer_factory_get_instance, gauxc_delete
   USE gauxc_molecular_weights, ONLY: gauxc_molecular_weights_type, &
      & gauxc_molecular_weights_factory_type, gauxc_molecular_weights_factory_new, &
      & gauxc_molecular_weights_modify_weights, gauxc_get_instance, gauxc_delete, &
      & gauxc_molecular_weights_settings
   USE gauxc_xc_functional, ONLY: gauxc_functional_type, gauxc_functional_from_string, gauxc_delete
   USE gauxc_integrator, ONLY: gauxc_integrator_type, &
       & gauxc_integrator_new, gauxc_integrator_eval_exc_vxc_rks, &
       & gauxc_integrator_eval_exc_vxc_uks, gauxc_delete
   USE gauxc_external_hdf5_read, ONLY: gauxc_read_hdf5_record
   USE cp_dbcsr_api, ONLY: dbcsr_p_type
   USE input_section_types, ONLY: section_vals_type, section_vals_val_get
#ifdef GAUXC_HAS_MPI
   USE mpi
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   INTERFACE
      SUBROUTINE print_gauxc_status_message(c_message) BIND(C, name="print_c_string")
         IMPORT :: c_ptr
         TYPE(c_ptr), VALUE :: c_message
      END SUBROUTINE print_gauxc_status_message
   END INTERFACE

   TYPE cp_gauxc_molecule_type
      TYPE(gauxc_molecule_type) :: molecule
   END TYPE

   TYPE cp_gauxc_basisset_type
      TYPE(gauxc_basisset_type) :: basis
   END TYPE

   TYPE cp_gauxc_grid_type
      TYPE(gauxc_molgrid_type) :: grid
      TYPE(gauxc_load_balancer_type) :: lb
      TYPE(gauxc_load_balancer_factory_type) :: lbf
      TYPE(gauxc_molecular_weights_type) :: mw
      TYPE(gauxc_molecular_weights_factory_type) :: mwf
   END TYPE cp_gauxc_grid_type

   TYPE cp_gauxc_integrator_type
      TYPE(gauxc_functional_type) :: func
      TYPE(gauxc_integrator_type) :: integrator
   END TYPE cp_gauxc_integrator_type

   TYPE cp_gauxc_xc_type
      REAL(c_double) :: exc
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: vxc_scalar, vxc_zeta
   END TYPE cp_gauxc_xc_type

   ! Copy paste from docu, move these wherever adequate
   TYPE(gauxc_runtime_environment_type) :: rt

   INTEGER(c_int) :: world_rank, world_size

   PUBLIC :: gauxc_init, gauxc_finalize, gauxc_create_molecule, &
          & gauxc_create_basisset, gauxc_destroy_molecule, gauxc_destroy_basisset, &
          & gauxc_create_grid, gauxc_create_integrator, gauxc_destroy_grid, &
    & gauxc_destroy_integrator, gauxc_compute_xc, dbcsr_to_dense, dense_to_dbcsr, &
           & read_execution_space, read_atomic_grid_size, read_radial_quad, &
           & read_pruning_scheme, &
           & cp_gauxc_grid_type, cp_gauxc_integrator_type, cp_gauxc_xc_type, cp_gauxc_molecule_type, cp_gauxc_basisset_type

CONTAINS

! **************************************************************************************************
!> \brief Checks gauxc status and prints error message before aborting
!> \param status the gauxc status to check
! **************************************************************************************************
   SUBROUTINE check_gauxc_status(status)
      TYPE(gauxc_status_type), INTENT(IN)                :: status

      IF (status%code /= 0) THEN
         WRITE (*, '(a,1x,i0)') "GauXC returned with status code", status%code
         IF (c_associated(status%message)) THEN
            WRITE (0, '(a)', ADVANCE='no') "GauxC status message: ["
            CALL print_gauxc_status_message(status%message)
            WRITE (0, '(a)') "]"
         ELSE
            WRITE (0, '(a)') "GauxC status message: [null]"
         END IF
         CPABORT("GauXC returned with non-zero status code")
      END IF
   END SUBROUTINE check_gauxc_status

! **************************************************************************************************
!> \brief ...
!> \param mpi_comm ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_init(mpi_comm, status)
      INTEGER, INTENT(IN)                                :: mpi_comm
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      rt = gauxc_runtime_environment_new(status, mpi_comm)
      CALL check_gauxc_status(status)

      world_rank = gauxc_runtime_environment_comm_rank(status, rt)
      CALL check_gauxc_status(status)

      world_size = gauxc_runtime_environment_comm_size(status, rt)

      IF (world_rank == 0) THEN
         IF (status%code /= 0) THEN
            PRINT '(a,1x,i0)', "GauXC returned with status code", status%code
            PRINT '(a)', gauxc_status_message(status)
         END IF
      END IF
      CALL check_gauxc_status(status)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_finalize(status)
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_runtime_environment_delete(status, rt)
      CALL check_gauxc_status(status)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param natoms ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_molecule(particle_set, natoms, status) RESULT(res)
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(cp_gauxc_molecule_type)                       :: res

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: atomic_number, i
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gauxc_atom_type), ALLOCATABLE, DIMENSION(:)   :: atoms

      ALLOCATE (atoms(natoms))

      DO i = 1, natoms
         atomic_kind => particle_set(i)%atomic_kind
         CALL get_atomic_kind(atomic_kind, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, number=atomic_number)
         atoms(i)%atomic_number = INT(atomic_number, c_int64_t)
         atoms(i)%x = REAL(particle_set(i)%r(1), c_double)
         atoms(i)%y = REAL(particle_set(i)%r(2), c_double)
         atoms(i)%z = REAL(particle_set(i)%r(3), c_double)
      END DO

      res%molecule = gauxc_molecule_new_from_atoms(status, atoms, INT(natoms, c_size_t))
      CALL check_gauxc_status(status)

      DEALLOCATE (atoms)

   END FUNCTION gauxc_create_molecule

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param natoms ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_basisset(qs_kind_set, particle_set, natoms, status) RESULT(res)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(cp_gauxc_basisset_type)                       :: res

      INTEGER                                            :: iatom, ikind, iprim, iset, ishell, lval, &
                                                            nkind, npgf, nset, nshell, &
                                                            nshell_total, shell_index
      REAL(c_double), DIMENSION(3)                       :: shell_origin
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gauxc_shell_type), ALLOCATABLE, DIMENSION(:)  :: shells
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis

      nkind = SIZE(qs_kind_set)

      ALLOCATE (basis_set_list(nkind))
      ! fixme: schreibe version dieser funktion, die als eingabe direkt die basis_set_list nimmt
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      nshell_total = 0
      DO iatom = 1, natoms
         atomic_kind => particle_set(iatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         gto_basis => basis_set_list(ikind)%gto_basis_set
         CPASSERT(ASSOCIATED(gto_basis))
         nshell_total = nshell_total + SUM(gto_basis%nshell)
      END DO

      ALLOCATE (shells(nshell_total))

      shell_index = 0
      DO iatom = 1, natoms ! for each atom
         atomic_kind => particle_set(iatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         gto_basis => basis_set_list(ikind)%gto_basis_set
         CPASSERT(ASSOCIATED(gto_basis))

         shell_origin(1) = REAL(particle_set(iatom)%r(1), c_double)
         shell_origin(2) = REAL(particle_set(iatom)%r(2), c_double)
         shell_origin(3) = REAL(particle_set(iatom)%r(3), c_double)

         nset = gto_basis%nset
         CPASSERT(nset == SIZE(gto_basis%nshell))
         DO iset = 1, nset ! for each shell group
            nshell = gto_basis%nshell(iset)
            npgf = gto_basis%npgf(iset) ! corresponds with nprim of gauxc

            DO ishell = 1, gto_basis%nshell(iset) ! for each shell within the shell group
               shell_index = shell_index + 1 ! global shell index, flattened over atoms and groups
               lval = gto_basis%l(ishell, iset)
               shells(shell_index)%l = INT(lval, c_int32_t)
               shells(shell_index)%pure = .TRUE._c_bool
               shells(shell_index)%nprim = INT(npgf, c_int32_t)
               shells(shell_index)%origin = shell_origin

               DO iprim = 1, npgf
                  shells(shell_index)%exponents(iprim) = REAL(gto_basis%zet(iprim, iset), c_double)
                  shells(shell_index)%coefficients(iprim) = REAL(gto_basis%gcc(iprim, ishell, iset), c_double)
               END DO
            END DO
         END DO
      END DO

      res%basis = gauxc_basisset_new_from_shells(status, shells, INT(nshell_total, c_size_t))
      CALL check_gauxc_status(status)

      DEALLOCATE (shells)
      DEALLOCATE (basis_set_list)

   END FUNCTION gauxc_create_basisset

! **************************************************************************************************
!> \brief ...
!> \param molecule ...
!> \param basis ...
!> \param gauxc_section ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_grid(molecule, basis, gauxc_section, status) RESULT(res)
      TYPE(cp_gauxc_molecule_type), INTENT(IN)           :: molecule
      TYPE(cp_gauxc_basisset_type), INTENT(in)           :: basis
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(cp_gauxc_grid_type)                           :: res

      CHARACTER(len=default_string_length)               :: grid_spec, lb_exec_space_str, &
                                                            prune_spec, rad_quad_spec
      INTEGER                                            :: batch_val, grid_val, lb_exec_val, &
                                                            prune_val, rad_quad_val
      INTEGER(c_int)                                     :: grid_type_local, int_exec_space_local, &
                                                            lb_exec_space_local, &
                                                            pruning_scheme_local, radial_quad_local
      INTEGER(c_int64_t)                                 :: batch_size

      CALL section_vals_val_get(gauxc_section, "GRID", i_val=grid_val)
      CALL section_vals_val_get(gauxc_section, "RADIAL_QUADRATURE", i_val=rad_quad_val)
      CALL section_vals_val_get(gauxc_section, "PRUNING_SCHEME", i_val=prune_val)
      CALL section_vals_val_get(gauxc_section, "BATCH_SIZE", i_val=batch_val)
      CALL section_vals_val_get(gauxc_section, "LB_EXECUTION_SPACE", i_val=lb_exec_val)

      SELECT CASE (grid_val)
      CASE (1)
         grid_spec = "fine"
      CASE (2)
         grid_spec = "ultrafine"
      CASE (3)
         grid_spec = "superfine"
      CASE (4)
         grid_spec = "gm3"
      CASE (5)
         grid_spec = "gm5"
      CASE DEFAULT
         grid_spec = "fine"
      END SELECT

      SELECT CASE (rad_quad_val)
      CASE (1)
         rad_quad_spec = "becke"
      CASE (2)
         rad_quad_spec = "muraknowles"
      CASE (3)
         rad_quad_spec = "treutlerahlrichs"
      CASE (4)
         rad_quad_spec = "murrayhandylaming"
      CASE DEFAULT
         rad_quad_spec = "muraknowles"
      END SELECT

      SELECT CASE (prune_val)
      CASE (1)
         prune_spec = "unpruned"
      CASE (2)
         prune_spec = "robust"
      CASE (3)
         prune_spec = "treutler"
      CASE DEFAULT
         prune_spec = "robust"
      END SELECT

      SELECT CASE (lb_exec_val)
      CASE (1)
         lb_exec_space_str = "host"
      CASE (2)
         lb_exec_space_str = "device"
      CASE DEFAULT
         lb_exec_space_str = "host"
      END SELECT

      batch_size = INT(batch_val, c_int64_t)

      grid_type_local = read_atomic_grid_size(grid_spec)
      radial_quad_local = read_radial_quad(rad_quad_spec)
      pruning_scheme_local = read_pruning_scheme(prune_spec)
      lb_exec_space_local = read_execution_space(lb_exec_space_str)
      int_exec_space_local = read_execution_space("host")

      res%grid = gauxc_molgrid_new_default(status, molecule%molecule, pruning_scheme_local, &
         & batch_size, radial_quad_local, grid_type_local)
      CALL check_gauxc_status(status)

      res%lbf = gauxc_load_balancer_factory_new(status, lb_exec_space_local)
      CALL check_gauxc_status(status)

      res%lb = gauxc_load_balancer_factory_get_instance(status, res%lbf, rt, molecule%molecule, res%grid, basis%basis)
      CALL check_gauxc_status(status)

      res%mwf = gauxc_molecular_weights_factory_new(status, int_exec_space_local)
      CALL check_gauxc_status(status)

      res%mw = gauxc_get_instance(status, res%mwf)
      CALL check_gauxc_status(status)

      CALL gauxc_molecular_weights_modify_weights(status, res%mw, res%lb)
      CALL check_gauxc_status(status)

   END FUNCTION gauxc_create_grid

! **************************************************************************************************
!> \brief ...
!> \param xc_functional_name ...
!> \param lb ...
!> \param gauxc_section ...
!> \param nspins ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_integrator(xc_functional_name, lb, gauxc_section, nspins, status) RESULT(res)
      CHARACTER(len=*), INTENT(IN)                       :: xc_functional_name
      TYPE(gauxc_load_balancer_type), INTENT(IN)         :: lb
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(cp_gauxc_integrator_type)                     :: res

      CHARACTER(len=default_string_length)               :: int_exec_space_str
      INTEGER(c_int)                                     :: int_exec_space_local, int_exec_val
      LOGICAL(c_bool)                                    :: polarized

      CALL section_vals_val_get(gauxc_section, "INT_EXECUTION_SPACE", i_val=int_exec_val)

      SELECT CASE (int_exec_val)
      CASE (1)
         int_exec_space_str = "host"
      CASE (2)
         int_exec_space_str = "device"
      CASE DEFAULT
         int_exec_space_str = "host"
      END SELECT

      polarized = (nspins == 2)
      res%func = gauxc_functional_from_string(status, xc_functional_name, polarized)
      CALL check_gauxc_status(status)

      int_exec_space_local = read_execution_space(int_exec_space_str)
      res%integrator = gauxc_integrator_new(status, res%func, lb, int_exec_space_local)
      CALL check_gauxc_status(status)

   END FUNCTION gauxc_create_integrator

! **************************************************************************************************
!> \brief ...
!> \param integrator ...
!> \param density_scalar ...
!> \param density_zeta ...
!> \param nspins ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_compute_xc(integrator, density_scalar, density_zeta, nspins, status) RESULT(res)
      TYPE(gauxc_integrator_type), INTENT(IN)            :: integrator
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: density_scalar
      REAL(c_double), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(cp_gauxc_xc_type)                             :: res

      ALLOCATE (res%vxc_scalar, mold=density_scalar)

      IF (nspins == 1) THEN
         CALL gauxc_integrator_eval_exc_vxc_rks(status, integrator, density_scalar, &
                                                res%exc, res%vxc_scalar)
      ELSE
         CPASSERT(PRESENT(density_zeta))
         ALLOCATE (res%vxc_zeta, mold=density_zeta)
         CALL gauxc_integrator_eval_exc_vxc_uks(status, integrator, density_scalar, &
                                                density_zeta, res%exc, res%vxc_scalar, res%vxc_zeta)
      END IF
      CALL check_gauxc_status(status)

   END FUNCTION gauxc_compute_xc

! **************************************************************************************************
!> \brief ...
!> \param molecule ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_molecule(molecule, status)
      TYPE(cp_gauxc_molecule_type), INTENT(INOUT)        :: molecule
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, molecule%molecule)
      CALL check_gauxc_status(status)
   END SUBROUTINE gauxc_destroy_molecule

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_basisset(basis, status)
      TYPE(cp_gauxc_basisset_type), INTENT(INOUT)        :: basis
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, basis%basis)
      CALL check_gauxc_status(status)
   END SUBROUTINE gauxc_destroy_basisset

! **************************************************************************************************
!> \brief ...
!> \param grid_result ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_grid(grid_result, status)
      TYPE(cp_gauxc_grid_type), INTENT(INOUT)            :: grid_result
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, grid_result%mw)
      CALL check_gauxc_status(status)
      CALL gauxc_delete(status, grid_result%mwf)
      CALL check_gauxc_status(status)
      CALL gauxc_delete(status, grid_result%lb)
      CALL check_gauxc_status(status)
      CALL gauxc_delete(status, grid_result%lbf)
      CALL check_gauxc_status(status)
      CALL gauxc_delete(status, grid_result%grid)
      CALL check_gauxc_status(status)
   END SUBROUTINE gauxc_destroy_grid

! **************************************************************************************************
!> \brief ...
!> \param integrator_result ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_integrator(integrator_result, status)
      TYPE(cp_gauxc_integrator_type), INTENT(INOUT)      :: integrator_result
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, integrator_result%integrator)
      CALL check_gauxc_status(status)
      CALL gauxc_delete(status, integrator_result%func)
      CALL check_gauxc_status(status)
   END SUBROUTINE gauxc_destroy_integrator

   ! **************************************************************************************************
   !> \brief ...
   !> \param dbcsr_mat ...
! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param dense_mat ...
!> \note UNKNOWN_COMMENT  **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense(dbcsr_mat, dense_mat)
      USE cp_dbcsr_api, ONLY: dbcsr_get_info, dbcsr_iterator_type, &
                              dbcsr_iterator_start, dbcsr_iterator_next_block, dbcsr_iterator_blocks_left, &
                              dbcsr_iterator_stop
      TYPE(dbcsr_p_type), INTENT(IN)                     :: dbcsr_mat
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: dense_mat

      INTEGER                                            :: icol, irow, ncols, nrows
      REAL(c_double), POINTER                            :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(dbcsr_mat%matrix, nfullrows_total=nrows, nfullcols_total=ncols)
      IF (.NOT. ALLOCATED(dense_mat)) THEN
         ALLOCATE (dense_mat(nrows, ncols))
         dense_mat = 0.0_c_double
      END IF

      CALL dbcsr_iterator_start(iter, dbcsr_mat%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block)
         dense_mat(irow:irow + SIZE(block, 1) - 1, icol:icol + SIZE(block, 2) - 1) = &
            dense_mat(irow:irow + SIZE(block, 1) - 1, icol:icol + SIZE(block, 2) - 1) + block(:, :)
      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE dbcsr_to_dense

   ! **************************************************************************************************
   !> \brief Convert a dense matrix to dbcsr format using a template matrix
   !> \param dense_mat Input dense matrix
   !> \param template_dbcsr Template dbcsr matrix for structure
   !> \param dbcsr_mat Output dbcsr matrix
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param dense_mat ...
!> \param template_dbcsr ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_to_dbcsr(dense_mat, template_dbcsr) RESULT(dbcsr_mat)
      USE cp_dbcsr_api, ONLY: dbcsr_get_info, dbcsr_iterator_type, &
                              dbcsr_iterator_start, dbcsr_iterator_next_block, &
                              dbcsr_iterator_blocks_left, dbcsr_iterator_stop, &
                              dbcsr_put_block, dbcsr_finalize, dbcsr_work_create, &
                              dbcsr_init_p, dbcsr_create, dbcsr_copy, dbcsr_release
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: dense_mat
      TYPE(dbcsr_p_type), INTENT(IN)                     :: template_dbcsr
      TYPE(dbcsr_p_type)                                 :: dbcsr_mat

      INTEGER                                            :: icol, irow, ncols, nrows
      REAL(c_double), POINTER                            :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(template_dbcsr%matrix, nfullrows_total=nrows, nfullcols_total=ncols)

      CALL dbcsr_init_p(dbcsr_mat%matrix)
      CALL dbcsr_create(dbcsr_mat%matrix, template=template_dbcsr%matrix, name="VXC from GauXC")
      CALL dbcsr_work_create(dbcsr_mat%matrix, work_mutable=.TRUE.)
      CALL dbcsr_finalize(dbcsr_mat%matrix)

      CALL dbcsr_iterator_start(iter, template_dbcsr%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block)
         CALL dbcsr_put_block(dbcsr_mat%matrix, irow, icol, &
                              dense_mat(irow:irow + SIZE(block, 1) - 1, &
                                        icol:icol + SIZE(block, 2) - 1))
      END DO
      CALL dbcsr_iterator_stop(iter)

   END FUNCTION dense_to_dbcsr

   ! **************************************************************************************************
   !> \brief ...
   !> \param spec ...
   !> \return ...
! **************************************************************************************************
!> \brief ...
!> \param spec ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION read_execution_space(spec) RESULT(val)
      CHARACTER(len=*), INTENT(IN)                       :: spec
      INTEGER(c_int)                                     :: val

      SELECT CASE (spec)
      CASE ("host")
         val = gauxc_executionspace%host
      CASE ("device")
         val = gauxc_executionspace%device
      CASE DEFAULT
         val = gauxc_executionspace%host
      END SELECT
   END FUNCTION read_execution_space

! **************************************************************************************************
!> \brief ...
!> \param spec ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION read_atomic_grid_size(spec) RESULT(val)
      CHARACTER(len=*), INTENT(IN)                       :: spec
      INTEGER(c_int)                                     :: val

      SELECT CASE (spec)
      CASE ("fine")
         val = gauxc_atomicgridsizedefault%finegrid
      CASE ("ultrafine")
         val = gauxc_atomicgridsizedefault%ultrafinegrid
      CASE ("superfine")
         val = gauxc_atomicgridsizedefault%superfinegrid
      CASE ("gm3")
         val = gauxc_atomicgridsizedefault%gm3
      CASE ("gm5")
         val = gauxc_atomicgridsizedefault%gm5
      CASE DEFAULT
         val = gauxc_atomicgridsizedefault%finegrid
      END SELECT
   END FUNCTION read_atomic_grid_size

! **************************************************************************************************
!> \brief ...
!> \param spec ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION read_radial_quad(spec) RESULT(val)
      CHARACTER(len=*), INTENT(IN)                       :: spec
      INTEGER(c_int)                                     :: val

      SELECT CASE (spec)
      CASE ("becke")
         val = gauxc_radialquad%becke
      CASE ("muraknowles")
         val = gauxc_radialquad%mura_knowles
      CASE ("treutlerahlrichs")
         val = gauxc_radialquad%treutler_ahlrichs
      CASE ("murrayhandylaming")
         val = gauxc_radialquad%murray_handy_laming
      CASE DEFAULT
         val = gauxc_radialquad%mura_knowles
      END SELECT
   END FUNCTION read_radial_quad

! **************************************************************************************************
!> \brief ...
!> \param spec ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION read_pruning_scheme(spec) RESULT(val)
      CHARACTER(len=*), INTENT(IN)                       :: spec
      INTEGER(c_int)                                     :: val

      SELECT CASE (spec)
      CASE ("unpruned")
         val = gauxc_pruningscheme%unpruned
      CASE ("robust")
         val = gauxc_pruningscheme%robust
      CASE ("treutler")
         val = gauxc_pruningscheme%treutler
      CASE DEFAULT
         val = gauxc_pruningscheme%robust
      END SELECT
   END FUNCTION read_pruning_scheme

END MODULE
