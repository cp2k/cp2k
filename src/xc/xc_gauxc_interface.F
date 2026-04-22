#include "gauxc/gauxc_config.f"

#define GAUXC_RETURN_IF_ERROR(status) IF (status%status%code /= 0) RETURN

MODULE xc_gauxc_interface

   USE iso_c_binding, ONLY: &
      c_associated, &
      c_bool, &
      c_char, &
      c_double, &
      c_f_pointer, &
      c_int, &
      c_int32_t, &
      c_int64_t, &
      c_null_ptr, &
      c_ptr, &
      c_size_t
   USE particle_types, ONLY: &
      particle_type
   USE qs_kind_types, ONLY: &
      get_qs_kind, &
      qs_kind_type
   USE input_section_types, ONLY: &
      section_vals_type, &
      section_vals_val_get
   USE cp_dbcsr_api, ONLY: &
      dbcsr_p_type

#ifdef __GAUXC

   USE kinds, ONLY: &
      default_string_length, &
      dp
   USE physcon, ONLY: &
      bohr
   USE atomic_kind_types, ONLY: &
      atomic_kind_type, &
      get_atomic_kind, &
      get_atomic_kind_set
   USE qs_integral_utils, ONLY: &
      basis_set_list_setup
   USE basis_set_types, ONLY: &
      gto_basis_set_p_type, &
      gto_basis_set_type, &
      write_gto_basis_set
   USE periodic_table, ONLY: &
      get_ptable_info
   USE gauxc_status, ONLY: &
      gauxc_status_message, &
      gauxc_status_type
   USE gauxc_enums, ONLY: &
      gauxc_atomicgridsizedefault, &
      gauxc_executionspace, &
      gauxc_pruningscheme, &
      gauxc_radialquad
   USE gauxc_runtime_environment, ONLY: &
      gauxc_runtime_environment_comm_rank, &
      gauxc_runtime_environment_comm_size, &
      gauxc_runtime_environment_delete, &
      gauxc_runtime_environment_new, &
      gauxc_runtime_environment_type
   USE gauxc_molecule, ONLY: &
      gauxc_delete, &
      gauxc_molecule_new_from_atoms, &
      gauxc_molecule_type
   USE gauxc_atom, ONLY: &
      gauxc_atom_type
   USE gauxc_basisset, ONLY: &
      gauxc_basisset_new, &
      gauxc_basisset_new_from_shells, &
      gauxc_basisset_type, &
      gauxc_delete
   USE gauxc_shell, ONLY: &
      gauxc_shell_type
   USE gauxc_molgrid, ONLY: &
      gauxc_delete, &
      gauxc_molgrid_new_default, &
      gauxc_molgrid_type
   USE gauxc_load_balancer, ONLY: &
      gauxc_delete, &
      gauxc_load_balancer_factory_get_instance, &
      gauxc_load_balancer_factory_new, &
      gauxc_load_balancer_factory_type, &
      gauxc_load_balancer_type
   USE gauxc_molecular_weights, ONLY: &
      gauxc_delete, &
      gauxc_get_instance, &
      gauxc_molecular_weights_factory_new, &
      gauxc_molecular_weights_factory_type, &
      gauxc_molecular_weights_modify_weights, &
      gauxc_molecular_weights_settings, &
      gauxc_molecular_weights_type
   USE gauxc_xc_functional, ONLY: &
      gauxc_delete, &
      gauxc_functional_from_string, &
      gauxc_functional_type
   USE gauxc_integrator, ONLY: &
      gauxc_delete, &
      gauxc_integrator_eval_exc_vxc_rks, &
      gauxc_integrator_eval_exc_vxc_uks, &
      gauxc_integrator_new, &
      gauxc_integrator_type
   USE gauxc_external_hdf5_read, ONLY: &
      gauxc_read_hdf5_record
   USE string_utilities, ONLY: &
      uppercase
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#ifndef __GAUXC

   ! The module still exists as an empty shell when compiling without GauXC.

   TYPE cp_gauxc_molecule_type
   END TYPE cp_gauxc_molecule_type

   TYPE cp_gauxc_basisset_type
   END TYPE cp_gauxc_basisset_type

   TYPE cp_gauxc_grid_type
   END TYPE cp_gauxc_grid_type

   TYPE cp_gauxc_integrator_type
   END TYPE cp_gauxc_integrator_type

   TYPE cp_gauxc_status_type
   END TYPE cp_gauxc_status_type

#else

   INTERFACE
      SUBROUTINE print_gauxc_status_message(c_message) BIND(C, name="print_c_string")
         IMPORT :: c_ptr
         TYPE(c_ptr), VALUE :: c_message
      END SUBROUTINE print_gauxc_status_message
   END INTERFACE

   ! TODO can we make the single fields private somehow?

   TYPE cp_gauxc_molecule_type
      TYPE(gauxc_molecule_type) :: molecule
   END TYPE cp_gauxc_molecule_type

   TYPE cp_gauxc_basisset_type
      TYPE(gauxc_basisset_type) :: basis
   END TYPE cp_gauxc_basisset_type

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

   TYPE cp_gauxc_status_type
      TYPE(gauxc_status_type) :: status
   END TYPE cp_gauxc_status_type

   TYPE(gauxc_runtime_environment_type) :: rt

#endif

   TYPE cp_gauxc_xc_type
      REAL(c_double) :: exc
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: vxc_scalar, vxc_zeta
   END TYPE cp_gauxc_xc_type

   CHARACTER(len=*), PARAMETER :: no_gauxc_message = "Compile CP2K with GauXC to use this functionality!"

   PUBLIC :: &
      cp_gauxc_basisset_type, &
      cp_gauxc_grid_type, &
      cp_gauxc_integrator_type, &
      cp_gauxc_molecule_type, &
      cp_gauxc_status_type, &
      cp_gauxc_xc_type, &
      gauxc_check_status, &
      gauxc_compute_xc, &
      gauxc_create_basisset, &
      gauxc_create_grid, &
      gauxc_create_integrator, &
      gauxc_create_molecule, &
      gauxc_destroy_basisset, &
      gauxc_destroy_grid, &
      gauxc_destroy_integrator, &
      gauxc_destroy_molecule, &
      gauxc_finalize, &
      gauxc_init
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mpi_comm ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_init(mpi_comm, status)
      INTEGER, INTENT(IN)                                :: mpi_comm
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status

#ifdef __GAUXC
      INTEGER(c_int)                                     :: world_rank

      rt = gauxc_runtime_environment_new(status%status, mpi_comm)
      GAUXC_RETURN_IF_ERROR(status)

      world_rank = gauxc_runtime_environment_comm_rank(status%status, rt)
      GAUXC_RETURN_IF_ERROR(status)
#else
      MARK_USED(mpi_comm)
      MARK_USED(status)
#endif
   END SUBROUTINE gauxc_init

! **************************************************************************************************
!> \brief ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_finalize(status)
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status

#ifdef __GAUXC
      CALL gauxc_runtime_environment_delete(status%status, rt)
      GAUXC_RETURN_IF_ERROR(status)
#else
      MARK_USED(status)
#endif
   END SUBROUTINE gauxc_finalize

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param natoms ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_molecule(particle_set, status) RESULT(res)
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status
      TYPE(cp_gauxc_molecule_type)                       :: res

#ifdef __GAUXC
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: atomic_number, i, natoms
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gauxc_atom_type), ALLOCATABLE, DIMENSION(:)   :: atoms

      natoms = SIZE(particle_set)
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

      res%molecule = gauxc_molecule_new_from_atoms(status%status, atoms, INT(natoms, c_size_t))
      GAUXC_RETURN_IF_ERROR(status)

      DEALLOCATE (atoms)
#else
      MARK_USED(particle_set)
      MARK_USED(natoms)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END FUNCTION gauxc_create_molecule

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param natoms ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_basisset(qs_kind_set, particle_set, status) RESULT(res)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status
      TYPE(cp_gauxc_basisset_type)                       :: res

#ifdef __GAUXC
      INTEGER                                            :: iatom, ikind, iprim, iset, ishell, lval, &
                                                            nkind, npgf, nset, nshell, &
                                                            nshell_total, shell_index, natoms
      REAL(c_double), DIMENSION(3)                       :: shell_origin
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gauxc_shell_type), ALLOCATABLE, DIMENSION(:)  :: shells
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis

      nkind = SIZE(qs_kind_set)
      natoms = SIZE(particle_set)

      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      nshell_total = 0
      DO iatom = 1, natoms
         atomic_kind => particle_set(iatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         gto_basis => basis_set_list(ikind)%gto_basis_set
         CPASSERT(ASSOCIATED(gto_basis))
         CALL write_gto_basis_set(gto_basis, 6, "–––––GAUXC-CREATE-BASISSET–––––")
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
               ! FIXME hardcoded true param
               ! pure=1: spherical Gaussians; pure=0: cartesian Gaussians
               shells(shell_index)%pure = .TRUE._c_bool
               shells(shell_index)%nprim = INT(npgf, c_int32_t)
               shells(shell_index)%origin = shell_origin

               DO iprim = 1, npgf
                  shells(shell_index)%exponents(iprim) = &
                     REAL(gto_basis%zet(iprim, iset), c_double)
                  shells(shell_index)%coefficients(iprim) = &
                     REAL(gto_basis%gcc(iprim, ishell, iset), c_double)
               END DO
            END DO
         END DO
      END DO

      res%basis = gauxc_basisset_new_from_shells(status%status, shells, INT(nshell_total, c_size_t))
      GAUXC_RETURN_IF_ERROR(status)

      DEALLOCATE (shells)
      DEALLOCATE (basis_set_list)

#else
      MARK_USED(qs_kind_set)
      MARK_USED(particle_set)
      MARK_USED(natoms)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
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
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status
      TYPE(cp_gauxc_grid_type)                           :: res

#ifdef __GAUXC
      CHARACTER(len=default_string_length)               :: grid_spec, lb_exec_space_str, &
                                                            prune_spec, rad_quad_spec
      INTEGER                                            :: batch_val
      INTEGER(c_int)                                     :: grid_type_local, int_exec_space_local, &
                                                            lb_exec_space_local, &
                                                            pruning_scheme_local, radial_quad_local
      INTEGER(c_int64_t)                                 :: batch_size

      CALL section_vals_val_get( &
         gauxc_section, &
         "GRID", &
         c_val=grid_spec)
      CALL section_vals_val_get( &
         gauxc_section, &
         "RADIAL_QUADRATURE", &
         c_val=rad_quad_spec)
      CALL section_vals_val_get( &
         gauxc_section, &
         "PRUNING_SCHEME", &
         c_val=prune_spec)
      CALL section_vals_val_get( &
         gauxc_section, &
         "BATCH_SIZE", &
         i_val=batch_val)
      CALL section_vals_val_get( &
         gauxc_section, &
         "LB_EXECUTION_SPACE", &
         c_val=lb_exec_space_str)

      batch_size = INT(batch_val, c_int64_t)

      grid_type_local = read_atomic_grid_size(grid_spec)
      radial_quad_local = read_radial_quad(rad_quad_spec)
      pruning_scheme_local = read_pruning_scheme(prune_spec)
      lb_exec_space_local = read_execution_space(lb_exec_space_str)
      int_exec_space_local = read_execution_space("host")

      res%grid = gauxc_molgrid_new_default( &
                 status%status, &
                 molecule%molecule, &
                 pruning_scheme_local, &
                 batch_size, &
                 radial_quad_local, &
                 grid_type_local)
      GAUXC_RETURN_IF_ERROR(status)

      res%lbf = gauxc_load_balancer_factory_new( &
                status%status, &
                lb_exec_space_local)
      GAUXC_RETURN_IF_ERROR(status)

      res%lb = gauxc_load_balancer_factory_get_instance( &
               status%status, &
               res%lbf, &
               rt, &
               molecule%molecule, &
               res%grid, &
               basis%basis)
      GAUXC_RETURN_IF_ERROR(status)

      res%mwf = gauxc_molecular_weights_factory_new( &
                status%status, &
                int_exec_space_local)
      GAUXC_RETURN_IF_ERROR(status)

      res%mw = gauxc_get_instance( &
               status%status, &
               res%mwf)
      GAUXC_RETURN_IF_ERROR(status)

      CALL gauxc_molecular_weights_modify_weights( &
         status%status, &
         res%mw, &
         res%lb)
      GAUXC_RETURN_IF_ERROR(status)

#else
      MARK_USED(molecule)
      MARK_USED(basis)
      MARK_USED(gauxc_section)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END FUNCTION gauxc_create_grid

! **************************************************************************************************
!> \brief ...
!> \param xc_functional_name ...
!> \param grid ...
!> \param gauxc_section ...
!> \param nspins ...
!> \param status ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_create_integrator(xc_functional_name, grid, gauxc_section, nspins, status) RESULT(res)
      CHARACTER(len=*), INTENT(IN)                       :: xc_functional_name
      TYPE(cp_gauxc_grid_type), INTENT(IN)         :: grid
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status
      TYPE(cp_gauxc_integrator_type)                     :: res

#ifdef __GAUXC
      CHARACTER(len=default_string_length)               :: int_exec_space_str
      INTEGER(c_int)                                     :: int_exec_space_local
      LOGICAL(c_bool)                                    :: polarized

      CALL section_vals_val_get( &
         gauxc_section, &
         "INT_EXECUTION_SPACE", &
         c_val=int_exec_space_str)

      polarized = (nspins == 2)
      res%func = gauxc_functional_from_string( &
                 status%status, &
                 xc_functional_name, &
                 polarized)
      GAUXC_RETURN_IF_ERROR(status)

      int_exec_space_local = read_execution_space(int_exec_space_str)
      res%integrator = gauxc_integrator_new( &
                       status%status, &
                       res%func, &
                       grid%lb, &
                       int_exec_space_local)
      GAUXC_RETURN_IF_ERROR(status)

#else
      MARK_USED(xc_functional_name)
      MARK_USED(grid)
      MARK_USED(gauxc_section)
      MARK_USED(nspins)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
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
      TYPE(cp_gauxc_integrator_type), INTENT(IN)            :: integrator
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: density_scalar
      REAL(c_double), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status
      TYPE(cp_gauxc_xc_type)                             :: res

#ifdef __GAUXC

      IF (.NOT. ALLOCATED(res%vxc_scalar)) THEN
         ALLOCATE (res%vxc_scalar, mold=density_scalar)
      ELSE
         CPASSERT(ALL(SHAPE(res%vxc_scalar) == SHAPE(density_scalar)))
      END IF
      res%vxc_scalar = 0._dp

      IF (nspins == 1) THEN
         CALL gauxc_integrator_eval_exc_vxc_rks( &
            status%status, &
            integrator%integrator, &
            density_scalar, &
            res%exc, &
            res%vxc_scalar)
      ELSE
         CPASSERT(PRESENT(density_zeta))

         IF (.NOT. ALLOCATED(res%vxc_zeta)) THEN
            ALLOCATE (res%vxc_zeta, mold=density_zeta)
         ELSE
            CPASSERT(ALL(SHAPE(res%vxc_zeta) == SHAPE(density_scalar)))
         END IF
         res%vxc_zeta = 0._dp

         CALL gauxc_integrator_eval_exc_vxc_uks( &
            status%status, &
            integrator%integrator, &
            density_scalar, &
            density_zeta, &
            res%exc, &
            res%vxc_scalar, &
            res%vxc_zeta)
      END IF
      GAUXC_RETURN_IF_ERROR(status)

#else
      MARK_USED(integrator)
      MARK_USED(density_scalar)
      MARK_USED(density_zeta)
      MARK_USED(nspins)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END FUNCTION gauxc_compute_xc

! **************************************************************************************************
!> \brief ...
!> \param molecule ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_molecule(molecule, status)
      TYPE(cp_gauxc_molecule_type), INTENT(INOUT)        :: molecule
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status

#ifdef __GAUXC
      CALL gauxc_delete(status%status, molecule%molecule)
      GAUXC_RETURN_IF_ERROR(status)
#else
      MARK_USED(molecule)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END SUBROUTINE gauxc_destroy_molecule

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_basisset(basis, status)
      TYPE(cp_gauxc_basisset_type), INTENT(INOUT)        :: basis
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status

#ifdef __GAUXC
      CALL gauxc_delete(status%status, basis%basis)
      GAUXC_RETURN_IF_ERROR(status)
#else
      MARK_USED(basis)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END SUBROUTINE gauxc_destroy_basisset

! **************************************************************************************************
!> \brief ...
!> \param grid_result ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_grid(grid_result, status)
      TYPE(cp_gauxc_grid_type), INTENT(INOUT)            :: grid_result
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status

#ifdef __GAUXC
      CALL gauxc_delete(status%status, grid_result%mw)
      GAUXC_RETURN_IF_ERROR(status)
      CALL gauxc_delete(status%status, grid_result%mwf)
      GAUXC_RETURN_IF_ERROR(status)
      CALL gauxc_delete(status%status, grid_result%lb)
      GAUXC_RETURN_IF_ERROR(status)
      CALL gauxc_delete(status%status, grid_result%lbf)
      GAUXC_RETURN_IF_ERROR(status)
      CALL gauxc_delete(status%status, grid_result%grid)
      GAUXC_RETURN_IF_ERROR(status)
#else
      MARK_USED(grid_result)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END SUBROUTINE gauxc_destroy_grid

! **************************************************************************************************
!> \brief ...
!> \param integrator_result ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_destroy_integrator(integrator_result, status)
      TYPE(cp_gauxc_integrator_type), INTENT(INOUT)      :: integrator_result
      TYPE(cp_gauxc_status_type), INTENT(OUT)            :: status

#ifdef __GAUXC
      CALL gauxc_delete(status%status, integrator_result%integrator)
      GAUXC_RETURN_IF_ERROR(status)
      CALL gauxc_delete(status%status, integrator_result%func)
      GAUXC_RETURN_IF_ERROR(status)
#else
      MARK_USED(integrator_result)
      MARK_USED(status)
      CPABORT(no_gauxc_message)
#endif
   END SUBROUTINE gauxc_destroy_integrator

! **************************************************************************************************
!> \brief Checks gauxc status and prints error message before aborting
!> \param status the gauxc status to check
! **************************************************************************************************
   SUBROUTINE gauxc_check_status(status)
      TYPE(cp_gauxc_status_type), INTENT(IN)                :: status

#ifdef __GAUXC
      IF (status%status%code /= 0) THEN
         WRITE (0, '(a,1x,i0)') "GauXC returned with status code", status%status%code
         IF (c_associated(status%status%message)) THEN
            WRITE (0, '(a)', ADVANCE='no') "GauXC status message: ["
            CALL print_gauxc_status_message(status%status%message)
            WRITE (0, '(a)') "]"
         ELSE
            WRITE (0, '(a)') "GauXC status message: [null]"
         END IF
         CPABORT("GauXC returned with non-zero status code")
      END IF
#else
      MARK_USED(status)
#endif
   END SUBROUTINE gauxc_check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From here on, it's private helpers !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef __GAUXC

! **************************************************************************************************
!> \brief ...
!> \param spec ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION read_execution_space(spec) RESULT(val)
      CHARACTER(len=*), INTENT(IN)                       :: spec
      INTEGER(c_int)                                     :: val

      CHARACTER(len=LEN(spec))                           :: spec_upper

      spec_upper = spec
      CALL uppercase(spec_upper)

      SELECT CASE (spec_upper)
      CASE ("HOST")
         val = gauxc_executionspace%host
      CASE ("DEVICE")
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

      CHARACTER(len=LEN(spec))                           :: spec_upper

      spec_upper = spec
      CALL uppercase(spec_upper)

      SELECT CASE (spec_upper)
      CASE ("FINE")
         val = gauxc_atomicgridsizedefault%finegrid
      CASE ("ULTRAFINE")
         val = gauxc_atomicgridsizedefault%ultrafinegrid
      CASE ("SUPERFINE")
         val = gauxc_atomicgridsizedefault%superfinegrid
      CASE ("GM3")
         val = gauxc_atomicgridsizedefault%gm3
      CASE ("GM5")
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

      CHARACTER(len=LEN(spec))                           :: spec_upper

      spec_upper = spec
      CALL uppercase(spec_upper)

      SELECT CASE (spec_upper)
      CASE ("BECKE")
         val = gauxc_radialquad%becke
      CASE ("MURAKNOWLES")
         val = gauxc_radialquad%mura_knowles
      CASE ("TREUTLERAHLRICHS")
         val = gauxc_radialquad%treutler_ahlrichs
      CASE ("MURRAYHANDYLAMING")
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

      CHARACTER(len=LEN(spec))                           :: spec_upper

      spec_upper = spec
      CALL uppercase(spec_upper)

      SELECT CASE (spec_upper)
      CASE ("UNPRUNED")
         val = gauxc_pruningscheme%unpruned
      CASE ("ROBUST")
         val = gauxc_pruningscheme%robust
      CASE ("TREUTLER")
         val = gauxc_pruningscheme%treutler
      CASE DEFAULT
         val = gauxc_pruningscheme%robust
      END SELECT
   END FUNCTION read_pruning_scheme

#endif

END MODULE xc_gauxc_interface
