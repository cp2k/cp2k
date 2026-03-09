#include "gauxc/gauxc_config.f"

MODULE xc_gauxc

   USE iso_c_binding, ONLY: c_int, c_int64_t, c_double, c_bool, c_size_t, c_int32_t
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
     & gauxc_runtime_environment_new, gauxc_delete
   USE gauxc_molecule, ONLY: gauxc_molecule_type, gauxc_molecule_new_from_atoms, gauxc_delete
   USE gauxc_atom, ONLY: gauxc_atom_type
   USE gauxc_basisset, ONLY: gauxc_basisset_type, gauxc_basisset_new, &
      & gauxc_basisset_new_from_shells, gauxc_delete
   USE gauxc_shell, ONLY: gauxc_shell_type
   USE gauxc_molgrid, ONLY: gauxc_molgrid_type, gauxc_molgrid_new_default, &
     & gauxc_delete
   USE gauxc_load_balancer, ONLY: gauxc_load_balancer_type, gauxc_load_balancer_factory_type, &
     & gauxc_load_balancer_factory_new, gauxc_delete
   USE gauxc_molecular_weights, ONLY: gauxc_molecular_weights_type, &
     & gauxc_molecular_weights_factory_type, gauxc_molecular_weights_factory_new, &
     & gauxc_molecular_weights_modify_weights, gauxc_get_instance, gauxc_delete, &
     & gauxc_molecular_weights_settings
   USE gauxc_xc_functional, ONLY: gauxc_functional_type, gauxc_functional_from_string, gauxc_delete
   USE gauxc_integrator, ONLY: gauxc_integrator_type, &
     & gauxc_eval_exc_vxc, gauxc_delete
   USE gauxc_external_hdf5_read, ONLY: gauxc_read_hdf5_record
#ifdef GAUXC_HAS_MPI
   USE mpi
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Copy paste from docu, move these wherever adequate
   TYPE(gauxc_status_type) :: status
   TYPE(gauxc_runtime_environment_type) :: rt
   TYPE(gauxc_molecule_type) :: mol
   TYPE(gauxc_basisset_type) :: basis
   TYPE(gauxc_molgrid_type) :: grid
   TYPE(gauxc_load_balancer_factory_type) :: lbf
   TYPE(gauxc_load_balancer_type) :: lb
   TYPE(gauxc_molecular_weights_factory_type) :: mwf
   TYPE(gauxc_molecular_weights_type) :: mw
   TYPE(gauxc_functional_type) :: func
   TYPE(gauxc_integrator_type) :: integrator

   ! Copy paste from docu, move these wherever adequate
   CHARACTER(len=:), ALLOCATABLE :: input_file, model, grid_spec, rad_quad_spec, prune_spec, &
  & lb_exec_space_str, int_exec_space_str
   INTEGER(c_int) :: grid_type, radial_quad, pruning_scheme
   INTEGER(c_int) :: lb_exec_space, int_exec_space
   INTEGER(c_int) :: world_rank, world_size
   INTEGER(c_int64_t) :: batch_size
   REAL(c_double) :: exc, t_exc, t_start, t_end
   INTEGER :: error

   PUBLIC :: gauxc_init, gauxc_finalize, gauxc_create_molecule, &
      & gauxc_create_basisset, gauxc_destroy_molecule, gauxc_destroy_basisset

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mpi_comm ...
! **************************************************************************************************
   SUBROUTINE gauxc_init(mpi_comm)
      INTEGER:: mpi_comm

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
      CPASSERT(status%code == 0)

      IF (world_rank == 0) THEN
         IF (status%code /= 0) THEN
            PRINT '(a,1x,i0)', "GauXC returned with status code", status%code
            PRINT '(a)', gauxc_status_message(status)
         END IF

         ! This output is prolly not supposed to come at init
         PRINT '(a)', &
           & "Configuration", &
           & "-> Input file        : "//input_file, &
           & "-> Model             : "//model, &
           & "-> Grid              : "//grid_spec, &
           & "-> Radial quadrature : "//rad_quad_spec, &
           & "-> Pruning scheme    : "//prune_spec, &
           & ""
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE gauxc_finalize()
      CALL gauxc_delete(status, rt)
      CPASSERT(status%code == 0)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param natoms ...
!> \param molecule ...
! **************************************************************************************************
   SUBROUTINE gauxc_create_molecule(particle_set, natoms, molecule)
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(gauxc_molecule_type), INTENT(OUT)             :: molecule

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

      molecule = gauxc_molecule_new_from_atoms(status, atoms, INT(natoms, c_size_t))
      CPASSERT(status%code == 0)

      DEALLOCATE (atoms)

   END SUBROUTINE gauxc_create_molecule

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param natoms ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE gauxc_create_basisset(qs_kind_set, particle_set, natoms, basis)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(gauxc_basisset_type), INTENT(OUT)             :: basis

      INTEGER                                            :: iatom, ikind, iprim, iset, &
                                                            istart, lval, nkind, npgf, nset, &
                                                            nshell, nshell_total, shell_index, &
                                                            local_ishell
      REAL(c_double), DIMENSION(3)                       :: shell_origin
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gauxc_shell_type), ALLOCATABLE, DIMENSION(:)  :: shells
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis

      nkind = SIZE(qs_kind_set)

      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

       nshell_total = 0
       DO iatom = 1, natoms
          atomic_kind => particle_set(iatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind, kind_number=ikind)
          gto_basis => basis_set_list(ikind)%gto_basis_set
          IF (ASSOCIATED(gto_basis)) THEN
             nshell_total = nshell_total + SUM(gto_basis%nshell)
          END IF
       END DO

      ALLOCATE (shells(nshell_total))

      shell_index = 0
      DO iatom = 1, natoms
         atomic_kind => particle_set(iatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)

         gto_basis => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(gto_basis)) CYCLE

         shell_origin(1) = REAL(particle_set(iatom)%r(1), c_double)
         shell_origin(2) = REAL(particle_set(iatom)%r(2), c_double)
         shell_origin(3) = REAL(particle_set(iatom)%r(3), c_double)

         nset = gto_basis%nset
         istart = 1
         DO iset = 1, nset
            nshell = gto_basis%nshell(iset)
            npgf = gto_basis%npgf(iset)

            DO local_ishell = istart, istart + nshell - 1
               shell_index = shell_index + 1
               lval = gto_basis%l(1, local_ishell)
               shells(shell_index)%l = INT(lval, c_int32_t)
               shells(shell_index)%pure = .TRUE._c_bool
               shells(shell_index)%nprim = INT(npgf, c_int32_t)
               shells(shell_index)%origin = shell_origin

               DO iprim = 1, npgf
                  shells(shell_index)%exponents(iprim) = REAL(gto_basis%zet(iprim, local_ishell), c_double)
                  shells(shell_index)%coefficients(iprim) = REAL(gto_basis%gcc(iprim, local_ishell, 1), c_double)
               END DO
            END DO

            istart = istart + nshell
         END DO
      END DO

      basis = gauxc_basisset_new_from_shells(status, shells, INT(nshell_total, c_size_t))
      CPASSERT(status%code == 0)

      DEALLOCATE (shells)
      DEALLOCATE (basis_set_list)

   END SUBROUTINE gauxc_create_basisset

   ! **************************************************************************************************
   !> \brief ...
   !> \param molecule ...
   ! **************************************************************************************************
   SUBROUTINE gauxc_destroy_molecule(molecule)
      TYPE(gauxc_molecule_type), INTENT(INOUT)           :: molecule

      CALL gauxc_delete(status, molecule)
      CPASSERT(status%code == 0)
   END SUBROUTINE gauxc_destroy_molecule

   ! **************************************************************************************************
   !> \brief ...
   !> \param basis ...
   ! **************************************************************************************************
   SUBROUTINE gauxc_destroy_basisset(basis)
      TYPE(gauxc_basisset_type), INTENT(INOUT)           :: basis

      CALL gauxc_delete(status, basis)
      CPASSERT(status%code == 0)
   END SUBROUTINE gauxc_destroy_basisset

END MODULE
