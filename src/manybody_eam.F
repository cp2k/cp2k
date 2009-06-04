!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      EAM potential
!> \author CJM, I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
MODULE manybody_eam

  USE bibliography,                    ONLY: Foiles1986,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: ea_type,&
                                             eam_pot_type,&
                                             pair_potential_pp_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_force_eam, density_nonbond
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_eam'

CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author CJM
! *****************************************************************************
  SUBROUTINE density_nonbond (fist_nonbond_env, particle_set, cell, para_env, error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_nonbond', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iend, igrp, ikind, ilist, ipair, &
      iparticle, istart, jkind, kind_a, kind_b, nkinds, nparticle, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pot_types
    LOGICAL                                  :: check, do_eam, failure
    REAL(KIND=dp)                            :: fac, rab2, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rho
    TYPE(eam_pot_type), POINTER              :: eam_a, eam_b
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    failure = .FALSE.
    CALL timeset(routineN,handle)
    do_eam = .FALSE.
    CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, &
         potparm = potparm, r_last_update=r_last_update,&
         r_last_update_pbc=r_last_update_pbc,error=error)
    nkinds = SIZE ( potparm%pot, 1 )
    DO ikind = 1, nkinds
       DO jkind = ikind, nkinds
          IF (ANY(potparm%pot(ikind,jkind)%pot%type==ea_type)) do_eam = .TRUE.
       END DO
    END DO
    ! Only if EAM potential are present
    IF (do_eam) THEN
       ! Add citation
       CALL cite_reference(Foiles1986)
       ! Check that order of multiple defined potential in case we do EAM
       ALLOCATE(pot_types(SIZE(potparm%pot(1,1)%pot%type)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, SIZE(potparm%pot(1,1)%pot%type)
          pot_types(i) = potparm%pot(1,1)%pot%type(i)
       END DO

       NULLIFY(eam_a, eam_b)
       nparticle = SIZE (particle_set)
       ALLOCATE ( rho(nparticle), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       rho = 0._dp
       ! Starting the force loop
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          IF (neighbor_kind_pair%npairs ==0) CYCLE
          Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart = neighbor_kind_pair%grp_kind_start(igrp)
             iend   = neighbor_kind_pair%grp_kind_end(igrp)
             ikind  = neighbor_kind_pair%ij_kind(1,igrp)
             jkind  = neighbor_kind_pair%ij_kind(2,igrp)
             check  = (SIZE(pot_types)==SIZE(potparm%pot(ikind,jkind)%pot%type))
             CPPrecondition(check,cp_failure_level,routineP,error,failure)
             check  = ALL(pot_types==potparm%pot(ikind,jkind)%pot%type)
             CPPrecondition(check,cp_failure_level,routineP,error,failure)

             DO i = 1, SIZE(potparm%pot(ikind,jkind)%pot%type)
                IF (potparm%pot(ikind,jkind)%pot%type(i)/=ea_type) CYCLE
                rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
                CALL matvec_3x3(cell_v, cell%hmat,REAL(neighbor_kind_pair%cell_vector,KIND=dp))
                DO ipair = istart, iend
                   atom_a = neighbor_kind_pair%list(1,ipair)
                   atom_b = neighbor_kind_pair%list(2,ipair)
                   fac    = 1.0_dp
                   IF (atom_a==atom_b) fac = 0.5_dp
                   rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                   rab = rab+cell_v
                   rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                   IF (rab2 <= rab2_max) THEN
                      kind_a = particle_set(atom_a)%atomic_kind%kind_number
                      kind_b = particle_set(atom_b)%atomic_kind%kind_number
                      eam_a => potparm%pot(kind_a,kind_a)%pot%set(i)%eam
                      eam_b => potparm%pot(kind_b,kind_b)%pot%set(i)%eam
                      CALL get_rho_eam ( eam_a, eam_b, rab2, atom_a, atom_b, rho, fac)
                   END IF
                END DO
             END DO
          END DO Kind_Group_Loop
       END DO
       CALL mp_sum( rho, para_env%group )
       DO iparticle=1,nparticle
          particle_set(iparticle)%rho = rho(iparticle)
       END DO

       DEALLOCATE ( rho, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)

  END SUBROUTINE density_nonbond

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_rho_eam ( eam_a, eam_b, rab2, atom_a, atom_b, rho, fac)
    TYPE(eam_pot_type), POINTER              :: eam_a, eam_b
    REAL(dp), INTENT(IN)                     :: rab2
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(INOUT)                  :: rho(:)
    REAL(dp), INTENT(IN)                     :: fac

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_rho_eam', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: index
    REAL(dp)                                 :: qq, rab, rhoi, rhoj

    rab = SQRT(rab2)

    ! Particle A
    index = INT(rab/eam_b%drar) + 1
    IF ( index > eam_b%npoints) THEN
       index = eam_b%npoints
    ELSEIF (index < 1) THEN
       index = 1
    ENDIF
    qq = rab - eam_b%rval(index)
    rhoi = eam_b%rho(index) + qq * eam_b%rhop(index)

    ! Particle B
    index = INT(rab / eam_a%drar) + 1
    IF (index > eam_a%npoints) THEN
       index = eam_a%npoints
    ELSEIF (index < 1) THEN
       index = 1
    ENDIF
    qq = rab - eam_a%rval(index)
    rhoj = eam_a%rho(index) + qq * eam_a%rhop(index)

    rho(atom_a) = rho(atom_a) + rhoi*fac
    rho(atom_b) = rho(atom_b) + rhoj*fac
  END SUBROUTINE get_rho_eam

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_force_eam ( rab2, eam_a, eam_b, particle_set, atom_a, atom_b, f_eam )
    REAL(dp), INTENT(IN)                     :: rab2
    TYPE(eam_pot_type), POINTER              :: eam_a, eam_b
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(OUT)                    :: f_eam

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_force_eam', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: index
    REAL(KIND=dp)                            :: denspi, denspj, fcp, qq, rab

    rab  = SQRT(rab2)

    ! Particle A
    index = INT(rab/eam_a%drar) + 1
    IF (index > eam_a%npoints) THEN
       index = eam_a%npoints
    ELSEIF (index < 1) THEN
       index = 1
    ENDIF
    qq = rab - eam_a%rval(index)
    denspi = eam_a%rhop(index) +  qq*(eam_a%rhop(index+1) - eam_a%rhop(index))/eam_a%drar


    ! Particle B
    index = INT(rab/eam_b%drar) + 1
    IF (index > eam_b%npoints) THEN
       index = eam_b%npoints
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF
    qq = rab - eam_b%rval(index)
    denspj = eam_b%rhop(index) + qq*(eam_b%rhop(index+1) - eam_b%rhop(index))/eam_b%drar

    fcp = denspj * particle_set(atom_a)%f_embed + denspi * particle_set(atom_b)%f_embed
    f_eam =  fcp / rab
  END SUBROUTINE get_force_eam

END MODULE manybody_eam

