!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of DFTB3 Terms
!> \author JGH
! *****************************************************************************
MODULE qs_dftb3_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb3_methods'

  PUBLIC :: build_dftb3_diagonal

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param mcharge ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_dftb3_diagonal(qs_env,ks_matrix,rho,mcharge,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(dp), DIMENSION(:)                   :: mcharge
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb3_diagonal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_i, atom_j, handle, i, &
                                                ia, iatom, ic, icol, ikind, &
                                                irow, is, istat, jatom, &
                                                jkind, natom, nimg
    INTEGER, DIMENSION(3)                    :: cellind
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, kind_of
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: failure, found, use_virial
    REAL(KIND=dp)                            :: dr, eb3, eloc, fi, gmij, ua, &
                                                ui, uj, zeff
    REAL(KIND=dp), DIMENSION(3)              :: fij, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsblock, ksblock, pblock, &
                                                sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: n_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)
    NULLIFY(atprop)

    ! Energy
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,atprop=atprop,error=error)

    eb3 = 0.0_dp
    CALL get_qs_env(qs_env=qs_env,local_particles=local_particles,error=error)
    DO ikind=1,SIZE(local_particles%n_el)
       CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind, error=error)
       CALL get_dftb_atom_param(dftb_kind,dudq=ua,zeff=zeff)
       DO ia=1,local_particles%n_el(ikind)
          iatom=local_particles%list(ikind)%array(ia)
          eloc = - 1.0_dp/6.0_dp * ua * mcharge(iatom)**3
          eb3 = eb3 + eloc
          IF(atprop%energy) THEN
             ! we have to add the part not covered by 0.5*Tr(FP)
             eloc = -0.5_dp*eloc - 0.25_dp * ua * zeff * mcharge(iatom)**2
             atprop%atecoul(iatom) = atprop%atecoul(iatom) + eloc
          END IF
       END DO
    END DO
    CALL get_qs_env(qs_env=qs_env,para_env=para_env,error=error)
    CALL mp_sum(eb3,para_env%group)
    energy%dftb3 = eb3

    ! Forces and Virial
    IF ( calculate_forces ) THEN
       CALL get_qs_env(qs_env=qs_env,matrix_s_kp=matrix_s,natom=natom,force=force,&
            cell=cell,virial=virial,particle_set=particle_set,error=error)
       ! virial
       use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

       ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                kind_of=kind_of,atom_of_kind=atom_of_kind)
       CALL qs_rho_get(rho, rho_ao_kp=matrix_p, error=error)
       IF ( SIZE(matrix_p,1) == 2) THEN
          DO ic=1,SIZE(matrix_p,2)
             CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                  alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
          END DO
       END IF
       !
       nimg = SIZE(matrix_p,2)
       NULLIFY(cell_to_index)
       IF (nimg>1) THEN
          NULLIFY(kpoints)
          CALL get_qs_env(qs_env=qs_env,kpoints=kpoints,error=error)
          CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index,error=error)
       END IF
       NULLIFY(n_list)
       CALL get_qs_env(qs_env=qs_env,sab_orb=n_list,error=error)
       CALL neighbor_list_iterator_create(nl_iterator,n_list)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
              iatom=iatom,jatom=jatom,r=rij,cell=cellind)

          dr = SQRT(SUM(rij**2))
          IF(iatom==jatom .AND. dr < 1.0e-6_dp) CYCLE
 
          icol = MAX(iatom,jatom)
          irow = MIN(iatom,jatom)
 
          IF(nimg==1) THEN
             ic = 1
          ELSE
             ic = cell_to_index(cellind(1),cellind(2),cellind(3))
             CPPostcondition(ic > 0,cp_failure_level,routineP,error,failure)
          END IF

          ikind  = kind_of(iatom)
          atom_i = atom_of_kind(iatom)
          CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind, error=error)
          CALL get_dftb_atom_param(dftb_kind,dudq=ui)
          jkind  = kind_of(jatom)
          atom_j = atom_of_kind(jatom)
          CALL get_qs_kind(qs_kind_set(jkind),dftb_parameter=dftb_kind, error=error)
          CALL get_dftb_atom_param(dftb_kind,dudq=uj)
          !
          gmij = -0.5_dp*(ui*mcharge(iatom)**2 + uj*mcharge(jatom)**2)
          !
          NULLIFY(pblock)
          CALL cp_dbcsr_get_block_p(matrix=matrix_p(1,ic)%matrix,&
               row=irow,col=icol,block=pblock,found=found)
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
          DO i=1,3
             NULLIFY(dsblock)
             CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i,ic)%matrix,&
                  row=irow,col=icol,block=dsblock,found=found)
             CPPostcondition(found,cp_failure_level,routineP,error,failure)
             IF(irow==iatom) THEN
                fi = -gmij*SUM(pblock*dsblock)
             ELSE
                fi =  gmij*SUM(pblock*dsblock)
             END IF
             force(ikind)%rho_elec(i,atom_i) = force(ikind)%rho_elec(i,atom_i) + fi
             force(jkind)%rho_elec(i,atom_j) = force(jkind)%rho_elec(i,atom_j) - fi
             fij(i) = fi
          END DO
          IF (use_virial) THEN
             CALL virial_pair_force ( virial%pv_virial, 1._dp, fij, rij, error)
             IF(atprop%stress) THEN
                CALL virial_pair_force (atprop%atstress(:,:,iatom), 0.5_dp, fij, rij, error)
                CALL virial_pair_force (atprop%atstress(:,:,jatom), 0.5_dp, fij, rij, error)
             END IF
          END IF

       END DO
       CALL neighbor_list_iterator_release(nl_iterator)
       !
       DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF ( SIZE(matrix_p,1) == 2 ) THEN
          DO ic=1,SIZE(matrix_p,2)
             CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                               alpha_scalar=1.0_dp,beta_scalar=-1.0_dp,error=error)
          END DO
       END IF
    END IF

    ! KS matrix
    IF ( .NOT. just_energy ) THEN
       CALL get_qs_env(qs_env=qs_env,matrix_s_kp=matrix_s,natom=natom,error=error)
       ALLOCATE (kind_of(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,kind_of=kind_of)
       !
       nimg = SIZE(ks_matrix,2)
       NULLIFY(cell_to_index)
       IF (nimg>1) THEN
          NULLIFY(kpoints)
          CALL get_qs_env(qs_env=qs_env,kpoints=kpoints,error=error)
          CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index,error=error)
       END IF

       NULLIFY(n_list)
       CALL get_qs_env(qs_env=qs_env,sab_orb=n_list,error=error)
       CALL neighbor_list_iterator_create(nl_iterator,n_list)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
              iatom=iatom,jatom=jatom,r=rij,cell=cellind)

          icol = MAX(iatom,jatom)
          irow = MIN(iatom,jatom)

          IF(nimg==1) THEN
             ic = 1
          ELSE
             ic = cell_to_index(cellind(1),cellind(2),cellind(3))
             CPPostcondition(ic > 0,cp_failure_level,routineP,error,failure)
          END IF

          ikind  = kind_of(iatom)
          CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind, error=error)
          CALL get_dftb_atom_param(dftb_kind,dudq=ui)
          jkind  = kind_of(jatom)
          CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind, error=error)
          CALL get_dftb_atom_param(dftb_kind,dudq=uj)
          gmij = -0.5_dp*(ui*mcharge(iatom)**2 + uj*mcharge(jatom)**2)
          !
          NULLIFY(sblock)
          CALL cp_dbcsr_get_block_p(matrix=matrix_s(1,ic)%matrix,&
               row=irow,col=icol,block=sblock,found=found)
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
          DO is=1,SIZE(ks_matrix,1)
             NULLIFY(ksblock)
             CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is,ic)%matrix,&
                  row=irow,col=icol,block=ksblock,found=found)
             CPPostcondition(found,cp_failure_level,routineP,error,failure)
             ksblock = ksblock - 0.5_dp*gmij*sblock
          END DO

       END DO
       CALL neighbor_list_iterator_release(nl_iterator)
       !
       DEALLOCATE (kind_of,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_dftb3_diagonal

END MODULE qs_dftb3_methods

