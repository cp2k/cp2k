#include "gauxc/gauxc_config.f"

MODULE xc_gauxc

   USE iso_c_binding, ONLY: c_int, c_int64_t, c_double, c_bool, c_size_t, c_int32_t, c_null_ptr
   USE kinds, ONLY: dp
   USE particle_types, ONLY: particle_type
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
      & gauxc_integrator_new, gauxc_eval_exc_vxc, gauxc_delete
   USE gauxc_external_hdf5_read, ONLY: gauxc_read_hdf5_record
#ifdef GAUXC_HAS_MPI
   USE mpi
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE gauxc_grid_result_type
      TYPE(gauxc_molgrid_type) :: grid
      TYPE(gauxc_load_balancer_type) :: lb
      TYPE(gauxc_load_balancer_factory_type) :: lbf
      TYPE(gauxc_molecular_weights_type) :: mw
      TYPE(gauxc_molecular_weights_factory_type) :: mwf
   END TYPE gauxc_grid_result_type

   TYPE gauxc_integrator_result_type
      TYPE(gauxc_functional_type) :: func
      TYPE(gauxc_integrator_type) :: integrator
   END TYPE gauxc_integrator_result_type

   TYPE gauxc_xc_result_type
      REAL(c_double) :: exc
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: vxc_scalar, vxc_zeta
   END TYPE gauxc_xc_result_type

   ! Copy paste from docu, move these wherever adequate
   TYPE(gauxc_runtime_environment_type) :: rt

   CHARACTER(len=:), ALLOCATABLE :: grid_spec, rad_quad_spec, prune_spec, &
      & lb_exec_space_str, int_exec_space_str
   INTEGER(c_int) :: grid_type, radial_quad, pruning_scheme
   INTEGER(c_int) :: lb_exec_space, int_exec_space
   INTEGER(c_int) :: world_rank, world_size
   INTEGER(c_int64_t) :: batch_size

   PUBLIC :: gauxc_init, gauxc_finalize, gauxc_create_molecule, &
      & gauxc_create_basisset, gauxc_destroy_molecule, gauxc_destroy_basisset, &
      & gauxc_create_grid, gauxc_create_integrator, gauxc_destroy_grid, &
      & gauxc_destroy_integrator, gauxc_compute_xc, &
      & read_execution_space, read_atomic_grid_size, read_radial_quad, &
      & read_pruning_scheme, &
      & gauxc_grid_result_type, gauxc_integrator_result_type, gauxc_xc_result_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mpi_comm ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_init(mpi_comm, status)
      INTEGER, INTENT(IN)                                :: mpi_comm
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      grid_spec = "fine"
      rad_quad_spec = "muraknowles"
      prune_spec = "robust"
      lb_exec_space_str = "host"
      int_exec_space_str = "host"
      batch_size = 512_c_int64_t

      rt = gauxc_runtime_environment_new(status, mpi_comm)
      CPASSERT(status%code == 0)

      world_rank = gauxc_runtime_environment_comm_rank(status, rt)
      CPASSERT(status%code == 0)

      world_size = gauxc_runtime_environment_comm_size(status, rt)

      IF (world_rank == 0) THEN
         IF (status%code /= 0) THEN
            PRINT '(a,1x,i0)', "GauXC returned with status code", status%code
            PRINT '(a)', gauxc_status_message(status)
         END IF

         ! This output is prolly not supposed to come at init
         PRINT '(a)', &
           & "Configuration", &
           & "-> Grid              : "//grid_spec, &
           & "-> Radial quadrature : "//rad_quad_spec, &
           & "-> Pruning scheme    : "//prune_spec, &
           & ""
      END IF
      CPASSERT(status%code == 0)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_finalize(status)
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_runtime_environment_delete(status, rt)
      CPASSERT(status%code == 0)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param natoms ...
!> \param molecule ...
!> \param status ...
! **************************************************************************************************
   FUNCTION gauxc_create_molecule(particle_set, natoms, status) RESULT(res)
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(gauxc_molecule_type)                          :: res

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

      res = gauxc_molecule_new_from_atoms(status, atoms, INT(natoms, c_size_t))
      CPASSERT(status%code == 0)

      DEALLOCATE (atoms)

   END FUNCTION gauxc_create_molecule

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param natoms ...
!> \param basis ...
!> \param status ...
! **************************************************************************************************
   FUNCTION gauxc_create_basisset(qs_kind_set, particle_set, natoms, status) RESULT(res)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(gauxc_basisset_type)                          :: res

      CHARACTER(LEN=2)                                   :: element_symbol
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
         CALL get_atomic_kind(atomic_kind, element_symbol=element_symbol)
         NULLIFY (gto_basis)
         DO ikind = 1, nkind
            IF (qs_kind_set(ikind)%element_symbol .EQ. element_symbol) THEN
               gto_basis => basis_set_list(ikind)%gto_basis_set
               EXIT
            END IF
         END DO
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

      res = gauxc_basisset_new_from_shells(status, shells, INT(nshell_total, c_size_t))
      CPASSERT(status%code == 0)

      DEALLOCATE (shells)
      DEALLOCATE (basis_set_list)

   END FUNCTION gauxc_create_basisset

   FUNCTION gauxc_create_grid(molecule, basis, status) RESULT(res)
      TYPE(gauxc_molecule_type), INTENT(IN)              :: molecule
      TYPE(gauxc_basisset_type), INTENT(in)              :: basis
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(gauxc_grid_result_type)                       :: res

      INTEGER(c_int)                                     :: grid_type_local, int_exec_space_local, &
                                                            lb_exec_space_local, &
                                                            pruning_scheme_local, radial_quad_local

      grid_type_local = read_atomic_grid_size(grid_spec)
      radial_quad_local = read_radial_quad(rad_quad_spec)
      pruning_scheme_local = read_pruning_scheme(prune_spec)
      lb_exec_space_local = read_execution_space(lb_exec_space_str)
      int_exec_space_local = read_execution_space(int_exec_space_str)

      res%grid = gauxc_molgrid_new_default(status, molecule, pruning_scheme_local, &
         & batch_size, radial_quad_local, grid_type_local)
      CPASSERT(status%code == 0)

      res%lbf = gauxc_load_balancer_factory_new(status, lb_exec_space_local)
      CPASSERT(status%code == 0)

      res%lb = gauxc_load_balancer_factory_get_instance(status, res%lbf, rt, molecule, res%grid, basis)
      CPASSERT(status%code == 0)

      res%mwf = gauxc_molecular_weights_factory_new(status, int_exec_space_local)
      CPASSERT(status%code == 0)

      res%mw = gauxc_get_instance(status, res%mwf)
      CPASSERT(status%code == 0)

      CALL gauxc_molecular_weights_modify_weights(status, res%mw, res%lb)
      CPASSERT(status%code == 0)

   END FUNCTION gauxc_create_grid

   FUNCTION gauxc_create_integrator(xc_functional_name, status, lb) RESULT(res)
      CHARACTER(len=*), INTENT(IN)                       :: xc_functional_name
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(gauxc_load_balancer_type), INTENT(IN)         :: lb
      TYPE(gauxc_integrator_result_type)                 :: res

      LOGICAL(c_bool), PARAMETER                         :: hybrid = .TRUE._c_bool

      INTEGER(c_int)                                     :: int_exec_space_local

      res%func = gauxc_functional_from_string(status, xc_functional_name, hybrid)
      CPASSERT(status%code == 0)

      int_exec_space_local = read_execution_space(int_exec_space_str)
      res%integrator = gauxc_integrator_new(status, res%func, lb, int_exec_space_local)
      CPASSERT(status%code == 0)

   END FUNCTION gauxc_create_integrator

   FUNCTION gauxc_compute_xc(density_scalar, density_zeta, integrator, status) RESULT(res)
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: density_scalar, density_zeta
      TYPE(gauxc_integrator_type), INTENT(IN)            :: integrator
      TYPE(gauxc_status_type), INTENT(OUT)               :: status
      TYPE(gauxc_xc_result_type)                         :: res

      CALL gauxc_eval_exc_vxc(status, integrator, density_scalar, density_zeta, &
         & "PBE", res%exc, res%vxc_scalar, res%vxc_zeta)
      CPASSERT(status%code == 0)

   END FUNCTION gauxc_compute_xc

   SUBROUTINE gauxc_destroy_molecule(molecule, status)
      TYPE(gauxc_molecule_type), INTENT(INOUT)           :: molecule
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, molecule)
      CPASSERT(status%code == 0)
   END SUBROUTINE gauxc_destroy_molecule

   SUBROUTINE gauxc_destroy_basisset(basis, status)
      TYPE(gauxc_basisset_type), INTENT(INOUT)           :: basis
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, basis)
      CPASSERT(status%code == 0)
   END SUBROUTINE gauxc_destroy_basisset

   SUBROUTINE gauxc_destroy_grid(grid_result, status)
      TYPE(gauxc_grid_result_type), INTENT(INOUT)        :: grid_result
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, grid_result%mw)
      CPASSERT(status%code == 0)
      CALL gauxc_delete(status, grid_result%mwf)
      CPASSERT(status%code == 0)
      CALL gauxc_delete(status, grid_result%lb)
      CPASSERT(status%code == 0)
      CALL gauxc_delete(status, grid_result%lbf)
      CPASSERT(status%code == 0)
      CALL gauxc_delete(status, grid_result%grid)
      CPASSERT(status%code == 0)
   END SUBROUTINE gauxc_destroy_grid

   SUBROUTINE gauxc_destroy_integrator(integrator_result, status)
      TYPE(gauxc_integrator_result_type), INTENT(INOUT)  :: integrator_result
      TYPE(gauxc_status_type), INTENT(OUT)               :: status

      CALL gauxc_delete(status, integrator_result%integrator)
      CPASSERT(status%code == 0)
      CALL gauxc_delete(status, integrator_result%func)
      CPASSERT(status%code == 0)
   END SUBROUTINE gauxc_destroy_integrator

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
