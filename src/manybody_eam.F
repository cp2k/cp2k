!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      EAM potential
!> \author CJM, I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
MODULE manybody_eam

  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: ea_type,&
                                             eam_pot_type,&
                                             pair_potential_pp_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_force_eam,density_nonbond
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_eam'

CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author CJM
! *****************************************************************************
  SUBROUTINE density_nonbond ( fist_nonbond_env, particle_set, cell, para_env, error )
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_nonbond', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iend, igrp, ikind, ilist, ipair, &
      iparticle, istart, istat, jkind, nkinds, nparticle
    LOGICAL                                  :: do_eam, failure
    REAL(KIND=dp)                            :: flops, rab2, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rho
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    failure = .FALSE.
    CALL timeset ( routineN, 'I', 'Mflops', handle )
    flops = 0.0_dp
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
       NULLIFY(eam)
       nparticle = SIZE ( particle_set )
       ALLOCATE ( rho ( nparticle ), STAT = istat )
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "rho",nparticle*dp_size)
       rho ( : ) = 0._dp
       ! Starting the force loop
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          IF (neighbor_kind_pair%npairs ==0) CYCLE
          Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart  = neighbor_kind_pair%grp_kind_start(igrp)
             iend    = neighbor_kind_pair%grp_kind_end(igrp)
             ikind   = neighbor_kind_pair%ij_kind(1,igrp)
             jkind   = neighbor_kind_pair%ij_kind(2,igrp)
             DO i = 1, SIZE(potparm%pot(ikind,jkind)%pot%type)
                IF (potparm%pot(ikind,jkind)%pot%type(i)/=ea_type) CYCLE
                rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
                CALL matvec_3x3(cell_v, cell%hmat,REAL(neighbor_kind_pair%cell_vector,KIND=dp))
                DO ipair = istart, iend
                   atom_a = neighbor_kind_pair%list(1,ipair)
                   atom_b = neighbor_kind_pair%list(2,ipair)
                   rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                   rab = rab+cell_v
                   rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                   IF (rab2 <= rab2_max) THEN
                      eam => potparm%pot(ikind,jkind)%pot%set(i)%eam
                      CALL get_rho_eam ( eam, rab2, atom_a, atom_b, rho  )
                   END IF
                END DO
             END DO
          END DO Kind_Group_Loop
       END DO
       CALL mp_sum( rho, para_env%group )
       DO iparticle=1,nparticle
          particle_set(iparticle)%rho = rho(iparticle)
       END DO

       DEALLOCATE ( rho, STAT = istat )
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "deallocate rho")
    END IF
    CALL timestop ( flops, handle )

  END SUBROUTINE density_nonbond

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_rho_eam ( eam, rab2, atom_a, atom_b, rho )
    TYPE(eam_pot_type), POINTER              :: eam
    REAL(dp), INTENT(IN)                     :: rab2
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(INOUT)                  :: rho( : )

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_rho_eam', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: index
    REAL(dp)                                 :: qq, rab, rhoi, rhoj

    rab = SQRT ( rab2 )
    index = INT ( rab / eam % drar ) + 1

    IF ( index > eam%npoints) THEN
       index = eam%npoints
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF

    qq = rab - eam % rval ( index )

    rhoi = eam % rho ( index ) + qq * eam % rhop ( index )
    rhoj = eam % rho ( index ) + qq * eam % rhop ( index )

    rho ( atom_a ) = rho ( atom_a ) + rhoi
    rho ( atom_b ) = rho ( atom_b ) + rhoj
  END SUBROUTINE get_rho_eam

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_force_eam ( rab2, eam, particle_set, atom_a, atom_b, f_eam )
    REAL(dp), INTENT(IN)                     :: rab2
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(OUT)                    :: f_eam

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_force_eam', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: index
    REAL(dp)                                 :: denspi, denspj, fcp, qq, rab

    f_eam = 0._dp
    rab = SQRT ( rab2 )
    index = INT ( rab / eam % drar ) + 1

    IF ( index > eam%npoints ) THEN
       index = eam%npoints
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF
    qq = rab - eam % rval ( index )

    denspi = eam % rhop ( index ) +  &
         qq * ( eam % rhop ( index + 1 ) - &
         eam % rhop ( index ) ) / eam % drar
    denspj = eam % rhop ( index ) + &
         qq *( eam % rhop ( index + 1 ) - &
         eam % rhop ( index ) ) / eam % drar

    fcp = denspj * particle_set ( atom_a ) % f_embed + &
         denspi * particle_set ( atom_b ) % f_embed
    f_eam =  fcp / rab
  END SUBROUTINE get_force_eam

END MODULE manybody_eam

