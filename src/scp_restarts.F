!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routine handling the collection of the restart information for the SCP
!> \author C. J. Mundy
!> \History -) Reorganization of the SCP RESTART [tlaino] - Teodoro Laino
!>             Zurich Uiversity 03.2008 
! *****************************************************************************
MODULE scp_restarts
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE particle_types,                  ONLY: particle_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_restarts'
  
  PUBLIC :: write_scp_coeff_set, read_aux_coeff_set

CONTAINS
      
! *****************************************************************************
  SUBROUTINE write_scp_coeff_set (ires, scp_env, atomic_kind_set, particle_set,  error)

    INTEGER, INTENT(IN)                      :: ires
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_scp_coeff_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, iel, ikind, &
                                                irep, n_el, natom, natom_tot, &
                                                nkind, nsgf, nsgf_max, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cindex
    INTEGER, POINTER                         :: mindex( : )
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbuffer, pindex
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: coeffs
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (para_env, mindex)

    CALL get_scp_env(scp_env=scp_env,aux_coeff_set=aux_coeff_set,&
         para_env=para_env,error=error)

    IF (ASSOCIATED(aux_coeff_set)) THEN          
       natom_tot=SIZE(particle_set)
       nsgf_max=0
       n_el=0
       DO ikind=1,SIZE(atomic_kind_set)
          atomic_kind=> atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               aux_basis_set=aux_basis_set,&
               natom=natom)
          IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
          CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
               nsgf=nsgf)
          nsgf_max=MAX(nsgf_max,nsgf)
          n_el = n_el + natom
       END DO

       ALLOCATE (pbuffer(n_el,nsgf_max),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       ALLOCATE (cindex(natom_tot),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       ALLOCATE (mindex(n_el),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       ALLOCATE (pindex(natom_tot,nsgf_max),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   

       pindex=0.0_dp
       cindex=0
       mindex=0

       nkind=SIZE(aux_coeff_set%coeffs_of_kind)
       DO ikind=1, nkind
          coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
          IF(.NOT.ASSOCIATED(coeffs)) CYCLE
          DO iel= 1,coeffs%n_els
             iatom=aux_coeff_set%distribution%list(ikind)%array(iel)
             cindex(iatom)=iatom
             pindex(iatom,:)=coeffs%c(iel,:)
          END DO
       END DO
       CALL mp_sum(pindex, para_env%group)
       CALL mp_sum(cindex, para_env%group)

       iel = 0
       DO iatom = 1, natom_tot
          IF ( cindex ( iatom ) == 0 ) THEN
             CYCLE 
          ELSE
             iel = iel + 1
             mindex(iel) = cindex ( iatom )
             pbuffer(iel,:) = pindex ( iatom, : )
          END IF
       END DO
       
       IF (ires>0) WRITE(ires)n_el, nsgf_max, natom_tot
       IF (ires>0) WRITE(ires)mindex
       DO irep = 1, SIZE(mindex)
          IF (ires>0) WRITE(ires) pbuffer(irep,:)
       END DO
       DEALLOCATE (mindex,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       DEALLOCATE (pindex,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       DEALLOCATE (cindex,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (pbuffer,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
    END IF

    CALL timestop(handle)

  END SUBROUTINE write_scp_coeff_set

! *****************************************************************************
  SUBROUTINE read_aux_coeff_set(rst_unit, aux_coeff_set, para_env, error)

    INTEGER, INTENT(IN)                      :: rst_unit
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_aux_coeff_set', &
      routineP = moduleN//':'//routineN

    INTEGER :: iatom, iel, ikind, irep, natom_tot, natom_tot_read, &
      ncoef_atom_max, ncoef_tot, nel_tot, nel_tot_read, nkind, nsgf_max_read, &
      stat
    INTEGER, POINTER                         :: mindex( : )
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbuffer, pindex
    TYPE(aux_coeff_type), POINTER            :: coeffs

    failure = .FALSE.
    nel_tot=aux_coeff_set%nel_tot
    ncoef_tot=aux_coeff_set%ncoef_tot
    ncoef_atom_max=aux_coeff_set%ncoef_atom_max
    nkind=SIZE(aux_coeff_set%coeffs_of_kind)

    IF (para_env%ionode) THEN
       READ(rst_unit)nel_tot_read, nsgf_max_read, natom_tot_read
    END IF
    CALL mp_bcast(nel_tot_read,para_env%source,para_env%group)
    CALL mp_bcast(nsgf_max_read, para_env%source,para_env%group)
    CALL mp_bcast(natom_tot_read, para_env%source,para_env%group)
    
    CPPrecondition(nel_tot_read==nel_tot,cp_failure_level,routineP,error,failure)
    CPPrecondition(nsgf_max_read==ncoef_atom_max,cp_failure_level,routineP,error,failure)

    natom_tot = natom_tot_read
    ALLOCATE (pbuffer(nel_tot,nsgf_max_read),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (mindex(nel_tot),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (pindex(natom_tot,nsgf_max_read),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (para_env%ionode) THEN
       READ(rst_unit)mindex
       DO irep = 1, nel_tot
          READ(rst_unit) pbuffer(irep,:)
       END DO       
    END IF
    CALL mp_bcast(mindex,para_env%source,para_env%group)
    CALL mp_bcast(pbuffer,para_env%source,para_env%group)

    DO iel = 1, nel_tot
       iatom = mindex(iel)
       pindex(iatom,:) = pbuffer(iel,:)
    END DO
    
    nkind=SIZE(aux_coeff_set%coeffs_of_kind)
    DO ikind=1, nkind
       coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
       IF(.NOT.ASSOCIATED(coeffs)) CYCLE
       DO iel= 1,coeffs%n_els
          iatom=aux_coeff_set%distribution%list(ikind)%array(iel)
          coeffs%c(iel,:)=pindex(iatom,:)
       END DO
    END DO

    DEALLOCATE (mindex,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
    DEALLOCATE (pindex,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
    DEALLOCATE (pbuffer,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   

  END SUBROUTINE read_aux_coeff_set

END MODULE scp_restarts
