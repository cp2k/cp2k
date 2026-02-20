#include "gauxc/gauxc_config.f"

MODULE xc_gauxc

   USE iso_c_binding, ONLY: c_int, c_int64_t, c_double, c_bool
   ! use flap, only : command_line_interface
   USE gauxc_status, ONLY: gauxc_status_type, gauxc_status_message
   USE gauxc_enums, ONLY: gauxc_executionspace, gauxc_radialquad, &
     & gauxc_atomicgridsizedefault, gauxc_pruningscheme
   USE gauxc_runtime_environment, ONLY: gauxc_runtime_environment_type, &
     & gauxc_runtime_environment_comm_rank, gauxc_runtime_environment_comm_size, &
     & gauxc_runtime_environment_new, gauxc_delete
   USE gauxc_molecule, ONLY: gauxc_molecule_type, gauxc_molecule_new, gauxc_delete
   USE gauxc_basisset, ONLY: gauxc_basisset_type, gauxc_basisset_new, gauxc_delete
   USE gauxc_molgrid, ONLY: gauxc_molgrid_type, gauxc_molgrid_new_default, &
     & gauxc_delete
   USE gauxc_load_balancer, ONLY: gauxc_load_balancer_type, gauxc_load_balancer_factory_type, &
     & gauxc_load_balancer_factory_new, gauxc_get_shared_instance, gauxc_delete
   USE gauxc_molecular_weights, ONLY: gauxc_molecular_weights_type, &
     & gauxc_molecular_weights_factory_type, gauxc_molecular_weights_factory_new, &
     & gauxc_molecular_weights_modify_weights, gauxc_get_instance, gauxc_delete, &
     & gauxc_molecular_weights_settings
   USE gauxc_xc_functional, ONLY: gauxc_functional_type, gauxc_functional_from_string, gauxc_delete
   USE gauxc_integrator, ONLY: gauxc_integrator_type, gauxc_integrator_factory_type, &
     & gauxc_integrator_factory_new, gauxc_get_instance, gauxc_eval_exc_vxc, gauxc_delete
   USE gauxc_matrix, ONLY: gauxc_matrix_type, gauxc_matrix_empty, gauxc_matrix_data, gauxc_delete
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
   TYPE(gauxc_integrator_factory_type) :: intf
   TYPE(gauxc_integrator_type) :: integrator
   TYPE(gauxc_matrix_type) :: p_s, p_z, vxc_s, vxc_z

   ! Copy paste from docu, move these wherever adequate
   CHARACTER(len=:), ALLOCATABLE :: input_file, model, grid_spec, rad_quad_spec, prune_spec, &
  & lb_exec_space_str, int_exec_space_str
   INTEGER(c_int) :: grid_type, radial_quad, pruning_scheme
   INTEGER(c_int) :: lb_exec_space, int_exec_space
   INTEGER(c_int) :: world_rank, world_size
   INTEGER(c_int64_t) :: batch_size
   REAL(c_double) :: exc, t_exc, t_start, t_end
   INTEGER :: error

   PUBLIC :: gauxc_init, gauxc_finalize

CONTAINS

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

   SUBROUTINE gauxc_finalize()
      CALL gauxc_delete(status, rt)
      CPASSERT(status%code == 0)
   END SUBROUTINE

END MODULE
