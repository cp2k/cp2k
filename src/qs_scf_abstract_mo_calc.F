!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_abstract_mo_calc
   USE input_section_types,             ONLY: section_vals_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type

   IMPLICIT NONE

   TYPE, ABSTRACT :: AbstractMOCalc
   CONTAINS
      PROCEDURE(abstract_mo_calc_init), deferred, pass(self) :: init
      PROCEDURE(abstract_mo_calc_pre), deferred, pass(self)  :: pre
      PROCEDURE(abstract_mo_calc_run), deferred, pass(self)  :: run
   END TYPE

   INTERFACE
      SUBROUTINE abstract_mo_calc_init(self, qs_env, scf_env, scf_section)
         IMPORT AbstractMOCalc, qs_environment_type, qs_scf_env_type, section_vals_type

         CLASS(AbstractMOCalc), INTENT(INOUT)              :: self
         TYPE(qs_environment_type), POINTER, INTENT(IN)    :: qs_env
         TYPE(qs_scf_env_type), POINTER, INTENT(IN)        :: scf_env
         TYPE(section_vals_type), POINTER, INTENT(IN)      :: scf_section
      END SUBROUTINE abstract_mo_calc_init
   END INTERFACE

   INTERFACE
      SUBROUTINE abstract_mo_calc_pre(self)
         IMPORT AbstractMOCalc
         CLASS(AbstractMOCalc), INTENT(INOUT) :: self
      END SUBROUTINE abstract_mo_calc_pre
   END INTERFACE

   INTERFACE
      SUBROUTINE abstract_mo_calc_run(self, diis_step, energy_only)
         IMPORT AbstractMOCalc
         CLASS(AbstractMOCalc), INTENT(INOUT) :: self
         LOGICAL, INTENT(INOUT)               :: diis_step
         LOGICAL, INTENT(INOUT)               :: energy_only
      END SUBROUTINE abstract_mo_calc_run
   END INTERFACE
END MODULE qs_scf_abstract_mo_calc
