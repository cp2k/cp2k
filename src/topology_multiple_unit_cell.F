!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles the multiple unit cell option regarding atomic coordinates
!> \author Teodoro Laino [tlaino] - 05.2009
! *****************************************************************************
MODULE topology_multiple_unit_cell
  USE cell_types,                      ONLY: cell_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_multiple_unit_cell'

  PRIVATE

! *** Public parameters ***
  PUBLIC :: topology_muc

CONTAINS

! *****************************************************************************
!> \brief Handles the multiple_unit_cell for the atomic coordinates..
!> \author Teodoro Laino [tlaino] - 05.2009
! *****************************************************************************
  SUBROUTINE topology_muc(topology, subsys_section, error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_muc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: unit_str
    INTEGER                                  :: handle, i, ind, j, k, m, n, &
                                                natoms, nrep
    INTEGER, DIMENSION(:), POINTER           :: iwork, multiple_unit_cell
    LOGICAL                                  :: check, explicit, failure, &
                                                scale
    REAL(KIND=dp), DIMENSION(3)              :: trsl, trsl_i, trsl_j, trsl_k
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(multiple_unit_cell,iwork,cell)
    CALL section_vals_val_get(subsys_section, "TOPOLOGY%MULTIPLE_UNIT_CELL",&
                              i_vals=multiple_unit_cell, error=error)
    ! Fail is one of the value is set to zero..
    CALL cp_assert(ALL(multiple_unit_cell>0),cp_fatal_level,cp_assertion_failed,routineP,&
         "SUBSYS%TOPOLOGY%MULTIPLE_UNIT_CELL accepts only integer values larger than 0! "//&
         "A value of 0 or negative is meaningless!"//&
CPSourceFileRef)
    IF (ANY(multiple_unit_cell/=1)) THEN
       ! Check that the setup between CELL and TOPOLOGY is the same..
       CALL section_vals_val_get(subsys_section, "CELL%MULTIPLE_UNIT_CELL",&
                                 i_vals=iwork, error=error)
       CALL cp_assert(ALL(iwork==multiple_unit_cell),cp_fatal_level,cp_assertion_failed,routineP,&
         "SUBSYS%TOPOLOGY%MULTIPLE_UNIT_CELL and SUBSYS%CELL%MULTIPLE_UNIT_CELL have been "//&
         "setup to two different values!! Correct this error!!"//&
CPSourceFileRef)
       cell => topology%cell_muc
       natoms = topology%natoms*REAL(PRODUCT(multiple_unit_cell),KIND=dp)

       ! Check, if velocities are provided, that they are consistent in number with the atoms...
       work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
       CALL section_vals_get(work_section,explicit=explicit,error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(work_section,'_DEFAULT_KEYWORD_',n_rep_val=nrep,error=error)
          check = nrep==natoms
          CALL cp_assert(check,cp_fatal_level,cp_assertion_failed,routineP,&
               "Number of available entries in VELOCITY section is not compatible with the number "//&
               "of atoms!"//&
CPSourceFileRef)
       END IF

       CALL reallocate(topology%atom_info%id_molname, 1,    natoms)
       CALL reallocate(topology%atom_info%id_resname, 1,    natoms)
       CALL reallocate(topology%atom_info%resid,   1,    natoms)
       CALL reallocate(topology%atom_info%id_atmname, 1,    natoms)
       CALL reallocate(topology%atom_info%r,             1,3,1,natoms)
       CALL reallocate(topology%atom_info%atm_mass,      1,    natoms)
       CALL reallocate(topology%atom_info%atm_charge,    1,    natoms)
       CALL reallocate(topology%atom_info%occup,         1,    natoms)
       CALL reallocate(topology%atom_info%beta,          1,    natoms)
       CALL reallocate(topology%atom_info%id_element,       1,    natoms)
       ind = 0
       DO k = 1, multiple_unit_cell(3)
          trsl_k = cell%hmat(:,3)*REAL(k-1,KIND=dp)
          DO j = 1, multiple_unit_cell(2)
             trsl_j = cell%hmat(:,2)*REAL(j-1,KIND=dp)
             DO i = 1, multiple_unit_cell(1)
                trsl_i = cell%hmat(:,1)*REAL(i-1,KIND=dp)
                trsl   = trsl_i + trsl_j + trsl_k
                ind = ind + 1
                IF (ind==1) CYCLE

                ! loop over atoms
                n = (ind-1)*topology%natoms
                DO m = 1, topology%natoms
                    topology%atom_info%id_atmname(n+m)  =  topology%atom_info%id_atmname(m)
                    topology%atom_info%r(1,n+m)         =  topology%atom_info%r(1,m)+trsl(1)
                    topology%atom_info%r(2,n+m)         =  topology%atom_info%r(2,m)+trsl(2)
                    topology%atom_info%r(3,n+m)         =  topology%atom_info%r(3,m)+trsl(3)
                    topology%atom_info%id_molname(n+m)  =  topology%atom_info%id_molname(m)
                    topology%atom_info%id_resname(n+m)  =  topology%atom_info%id_resname(m)
                    topology%atom_info%resid(n+m)       =  topology%atom_info%resid(m)
                    topology%atom_info%id_element(n+m)  =  topology%atom_info%id_element(m)
                    topology%atom_info%atm_mass(n+m)    =  topology%atom_info%atm_mass(m)
                    topology%atom_info%atm_charge(n+m)  =  topology%atom_info%atm_charge(m)
                END DO
             END DO
          END DO
       END DO
       topology%natoms            = natoms

       ! Deallocate the coordinate section (will be rebuilt later with the whole atomic set)
       work_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
       CALL section_vals_get(work_section,explicit=explicit,error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,error=error)
          CALL section_vals_val_get(work_section,"SCALED",l_val=scale,error=error)
       END IF
       CALL section_vals_remove_values(work_section, error)
       IF (explicit) THEN
          CALL section_vals_val_set(work_section,"UNIT",c_val=unit_str,error=error)
          CALL section_vals_val_set(work_section,"SCALED",l_val=scale,error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE topology_muc

END MODULE topology_multiple_unit_cell
