!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Hamiltonian matrices in SCPTB
!> \author JGH
! *****************************************************************************
MODULE scptb_ks_matrix
  USE ai_coulomb,                      ONLY: coulomb2_new
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type,&
                                             scptb_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_get_block_p, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_p_type, cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE gamma,                           ONLY: fgamma_0
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE linear_systems,                  ONLY: solve_system
  USE mathconstants,                   ONLY: dfac,&
                                             fourpi,&
                                             oorootpi,&
                                             pi
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_charges
  USE orbital_pointers,                ONLY: indso,&
                                             ncoset,&
                                             nsoset
  USE orbital_transformation_matrices, ONLY: c2s_tramat
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_scp_charge
  USE qs_dftb_coulomb,                 ONLY: dftb_ewald_overlap,&
                                             dftb_spme_evaluate
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential,          ONLY: integrate_scp_rspace
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE scptb_types,                     ONLY: &
       get_scptb_parameter, scp_vector_add, scp_vector_copy, &
       scp_vector_create, scp_vector_dot, scp_vector_mult, scp_vector_p_type, &
       scp_vector_release, scp_vector_set, scp_vector_type, &
       scptb_parameter_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_ks_matrix'

  PUBLIC :: build_scptb_ks_matrix

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
! *****************************************************************************
  SUBROUTINE build_scptb_ks_matrix(qs_env,calculate_forces,just_energy)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'build_scptb_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, handle, iatom, ic, ikind, ispin, istat, lmaxscp, &
      natom, natom_kind, nkind, nspins, output_unit
    LOGICAL                                  :: do_es, do_scp, failure
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mcharge
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, dmcharge
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: ks_matrixkp, matrixkp_h, &
                                                matrixkp_p, matrixkp_s
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind
    TYPE(section_vals_type), POINTER         :: scf_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(dft_control, logger, scf_section)
    NULLIFY(particle_set, ks_env, ks_matrixkp, rho, energy)
    CPASSERT(ASSOCIATED(qs_env))
    logger => cp_get_default_logger()

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    para_env=para_env,&
                    rho=rho,&
                    energy=energy,&
                    matrix_h_kp=matrixkp_h,&
                    matrix_ks_kp=ks_matrixkp)

    CPASSERT(SIZE(ks_matrixkp,1)>0)
    CPASSERT(SIZE(ks_matrixkp,2)>0)

    scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF")
    nspins=dft_control%nspins
    CPASSERT(ASSOCIATED(rho))

    do_es = .FALSE.
    do_scp = .FALSE.
    nkind = SIZE(atomic_kind_set)
    DO ikind=1,nkind
      CALL get_atomic_kind(atomic_kind_set(ikind),natom=natom_kind)
      CALL get_qs_kind(qs_kind_set(ikind), scptb_parameter=scptb_kind)
      CALL get_scptb_parameter(scptb_kind,lmaxscp=lmaxscp)
      IF(lmaxscp > 0) do_scp = .TRUE.
      IF(lmaxscp >= 0) do_es = .TRUE.
    END DO
    scptb_control => dft_control%qs_control%scptb_control
    IF(.NOT.scptb_control%do_scp) do_scp = .FALSE.
    IF(.NOT.scptb_control%do_scc) do_es  = .FALSE.
    ! SCP requires self-consistent charge method
    IF(do_scp) do_es = .TRUE.

    ! copy the core matrix into the fock matrix
    CPASSERT(ASSOCIATED(matrixkp_h))
    DO ispin=1,nspins
       DO ic=1,SIZE(matrixkp_h,2)
          CALL cp_dbcsr_copy(ks_matrixkp(ispin,ic)%matrix,matrixkp_h(1,ic)%matrix)
       END DO
    END DO

    IF (do_es) THEN
       ! Mulliken charges
       CALL get_qs_env(qs_env,&
                       particle_set=particle_set,&
                       matrix_s_kp=matrixkp_s)
       CALL qs_rho_get(rho,rho_ao_kp=matrixkp_p)
       natom=SIZE(particle_set)
       ALLOCATE(charges(natom,nspins),STAT=istat)
       CPASSERT(istat==0)
       ALLOCATE(mcharge(natom),STAT=istat)
       CPASSERT(istat==0)
       ALLOCATE(dmcharge(natom,3),STAT=istat)
       CPASSERT(istat==0)
       IF(calculate_forces) THEN
          CALL mulliken_charges(matrixkp_p,matrixkp_s,para_env,charges,dmcharge)
       ELSE
          charges = 0._dp
          CALL mulliken_charges(matrixkp_p,matrixkp_s,para_env,charges)
       END IF
       nkind = SIZE(atomic_kind_set)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind),natom=natom_kind)
          CALL get_qs_kind(qs_kind_set(ikind),scptb_parameter=scptb_kind)
          CALL get_scptb_parameter(scptb_kind,zeff=zeff)
          DO iatom=1,natom_kind
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            mcharge(atom_a) = zeff - SUM(charges(atom_a,1:nspins))
          END DO
       END DO
       DEALLOCATE(charges,STAT=istat)
       CPASSERT(istat==0)

       IF (do_scp) THEN
          CALL scp_coulomb(qs_env,ks_matrixkp,matrixkp_s,rho,mcharge,calculate_forces,just_energy)
       ELSE
          CALL tb_coulomb(qs_env,ks_matrixkp,matrixkp_s,rho,mcharge,calculate_forces,just_energy)
       END IF

       DEALLOCATE(mcharge,STAT=istat)
       CPASSERT(istat==0)
       DEALLOCATE(dmcharge,STAT=istat)
       CPASSERT(istat==0)
    ELSE
       energy%hartree = 0._dp
    END IF

    energy%total = energy%core + energy%hartree + energy%qmmm_el + energy%repulsive + &
                   energy%dispersion

    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
       extension=".scfLog")
    IF (output_unit>0) THEN
        WRITE (UNIT=output_unit,FMT="(/,(T9,A,T60,F20.10))")&
             "Repulsive pair potential energy:               ",energy%repulsive,&
             "Zeroth order Hamiltonian energy:               ",energy%core,&
             "Total kinetic energy:                          ",energy%kinetic,&
             "Charge fluctuation energy:                     ",energy%hartree,&
             "London dispersion energy:                      ",energy%dispersion
        IF (qs_env%qmmm) THEN
           WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                "QM/MM Electrostatic energy:                    ",energy%qmmm_el
        END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%DETAILED_ENERGY")
    ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
    IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
       CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array)
       DO ispin=1,SIZE(mo_derivs)
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff_b=mo_coeff, occupation_numbers=occupation_numbers )
          IF(.NOT.mo_array(ispin)%mo_set%use_mo_coeff_b) THEN
             CPASSERT(.FALSE.)
          ENDIF
          CALL cp_dbcsr_multiply('n','n',1.0_dp,ks_matrixkp(1,ispin)%matrix,mo_coeff,&
               0.0_dp,mo_derivs(ispin)%matrix)
       ENDDO
    ENDIF


    CALL timestop(handle)

  END SUBROUTINE build_scptb_ks_matrix

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param matrix_s ...
!> \param rho ...
!> \param mcharge ...
!> \param calculate_forces ...
!> \param just_energy ...
! *****************************************************************************
  SUBROUTINE tb_coulomb(qs_env,ks_matrix,matrix_s,rho,mcharge,calculate_forces,just_energy)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: ks_matrix, matrix_s
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(inout)                          :: mcharge
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'tb_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, ewald_type, i, ia, iat, iatom, ic, icol, ikind, irow, is, &
      istat, ja, jat, jatom, jkind, lmaxscp, nai, naj, natom, ndim, nimg, &
      nkind, nsize
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    INTEGER, DIMENSION(3)                    :: cellind, periodic
    INTEGER, DIMENSION(:), POINTER           :: atomi_list, atomj_list
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: defined, found, use_virial
    REAL(KIND=dp)                            :: alpha, beta, chabs, deth, dq, &
                                                dr, dr2, ehartree, eta, f0, &
                                                ff, fi, gmij, noa, nob, t
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: qcharge
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drcharge, gmcharge
    REAL(KIND=dp), DIMENSION(0:1)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: fij, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsblock, ksblock, pblock, &
                                                sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: n_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind
    TYPE(virial_type), POINTER               :: virial

    NULLIFY(virial, atprop)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,virial=virial,atprop=atprop,&
         energy=energy,para_env=para_env,particle_set=particle_set,cell=cell)
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
         local_particles=local_particles,natom=natom)

    scptb_control => dft_control%qs_control%scptb_control
    use_virial=.FALSE.
    ndim = 1
    IF (calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,force=force)
       ndim = 4
       use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer).AND.scptb_control%do_ewald
    END IF

    nimg = dft_control%nimages
    NULLIFY(cell_to_index)
    IF (nimg>1) THEN
       NULLIFY(kpoints)
       CALL get_qs_env(qs_env=qs_env,kpoints=kpoints)
       CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index)
    END IF

    energy%hartree = 0._dp

    chabs = SUM(ABS(mcharge))
    CALL mp_sum(chabs,para_env%group)
    IF ( chabs > 1.e-8_dp ) THEN
       IF(scptb_control%do_ewald) THEN
          ALLOCATE(gmcharge(natom,ndim),STAT=istat)
          CPASSERT(istat==0)
          gmcharge = 0._dp
          ! Ewald sum
          NULLIFY(ewald_env,ewald_pw)
          CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw)
          CALL get_cell(cell=cell,periodic=periodic,deth=deth)
          CALL ewald_env_get(ewald_env,alpha=alpha,ewald_type=ewald_type)
          CALL get_qs_env(qs_env=qs_env,sab_tbe=n_list)
          CALL dftb_ewald_overlap(gmcharge,mcharge,alpha,n_list,&
                                  virial,use_virial,atprop)
          SELECT CASE(ewald_type)
          CASE DEFAULT
            CALL cp_unimplemented_error(routineP,"Invalid Ewald type")
          CASE(do_ewald_none)
            CPABORT("Not allowed with SCPTB")
            CPASSERT(.FALSE.)
          CASE(do_ewald_ewald)
            CALL cp_unimplemented_error(routineP,"Standard Ewald not implemented in SCPTB")
          CASE(do_ewald_pme)
            CALL cp_unimplemented_error(routineP,"PME not implemented in SCPTB")
          CASE(do_ewald_spme)
            CALL dftb_spme_evaluate(ewald_env,ewald_pw,particle_set,cell,&
              gmcharge,mcharge,calculate_forces,virial,use_virial,atprop)
          END SELECT

          ! add self charge interaction and background charge contribution
          CALL mp_sum(gmcharge(:,1),para_env%group)
          gmcharge(:,1) = gmcharge(:,1) - 2._dp*alpha*oorootpi*mcharge(:)
          IF ( ANY (periodic(:)==1) ) THEN
            gmcharge(:,1) = gmcharge(:,1) - pi / alpha**2 / deth
          END IF
          energy%hartree = energy%hartree + 0.5_dp*SUM(mcharge(:)*gmcharge(:,1))
          IF ( calculate_forces ) THEN
            ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
            CPASSERT(istat==0)

            CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

            gmcharge(:,2)=gmcharge(:,2)*mcharge(:)
            gmcharge(:,3)=gmcharge(:,3)*mcharge(:)
            gmcharge(:,4)=gmcharge(:,4)*mcharge(:)
            DO iatom=1,natom
              ikind  = kind_of(iatom)
              iat = atom_of_kind(iatom)
              force(ikind)%rho_elec(1,iat) = force(ikind)%rho_elec(1,iat) - gmcharge(iatom,2)
              force(ikind)%rho_elec(2,iat) = force(ikind)%rho_elec(2,iat) - gmcharge(iatom,3)
              force(ikind)%rho_elec(3,iat) = force(ikind)%rho_elec(3,iat) - gmcharge(iatom,4)
            END DO
          END IF
          IF ( .NOT. just_energy ) THEN
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            IF ( calculate_forces .AND. SIZE(matrix_p,1) == 2) THEN
              DO ic=1,nimg
                CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                     alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
              END DO
            END IF

            nsize = SIZE(ks_matrix,1)
            NULLIFY(n_list)
            CALL get_qs_env(qs_env=qs_env,sab_orb=n_list)
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
                  CPASSERT(ic > 0)
               END IF

               gmij = 0.5_dp*(gmcharge(iatom,1)+gmcharge(jatom,1))

               CALL cp_dbcsr_get_block_p(matrix=matrix_s(1,ic)%matrix,&
                    row=irow,col=icol,block=sblock,found=found)
               CPASSERT(found)
               DO is=1,nsize
                  NULLIFY(ksblock)
                  CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is,ic)%matrix,&
                       row=irow,col=icol,block=ksblock,found=found)
                  CPASSERT(found)
                  ksblock = ksblock - gmij*sblock
               END DO
               IF ( calculate_forces .AND. iatom /= jatom ) THEN
                  ikind  = kind_of(iatom)
                  iat = atom_of_kind(iatom)
                  jkind  = kind_of(jatom)
                  jat = atom_of_kind(jatom)
                  NULLIFY(pblock)
                  CALL cp_dbcsr_get_block_p(matrix=matrix_p(1,ic)%matrix,&
                       row=irow,col=icol,block=pblock,found=found)
                  CPASSERT(found)
                  DO i=1,3
                     NULLIFY(dsblock)
                     CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i,ic)%matrix,&
                          row=irow,col=icol,block=dsblock,found=found)
                     CPASSERT(found)
                     fi = -gmij*SUM(pblock*dsblock)*2.0_dp
                     force(ikind)%rho_elec(i,iat) = force(ikind)%rho_elec(i,iat) + fi
                     force(jkind)%rho_elec(i,jat) = force(jkind)%rho_elec(i,jat) - fi
                     fij(i) = fi
                  END DO
                  IF (use_virial) THEN
                     CALL virial_pair_force ( virial%pv_virial, 1._dp, fij, rij)
                     IF (atprop%stress) THEN
                       CALL virial_pair_force (atprop%atstress(:,:,iatom), 0.5_dp, fij, rij)
                       CALL virial_pair_force (atprop%atstress(:,:,jatom), 0.5_dp, fij, rij)
                     END IF
                  END IF
               END IF

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)

            IF ( calculate_forces .AND. SIZE(matrix_p,1) == 2) THEN
              DO ic=1,nimg
                CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                              alpha_scalar=1.0_dp,beta_scalar=-1.0_dp)
              END DO
            END IF
          END IF

          ! direct sum correction for Gaussians
          ehartree = 0._dp
          IF(.NOT.just_energy .OR. calculate_forces) THEN
            ALLOCATE(qcharge(natom),STAT=istat)
            CPASSERT(istat==0)
            qcharge = 0._dp
          END IF
          IF (calculate_forces) THEN
            ALLOCATE(drcharge(3,natom),STAT=istat)
            CPASSERT(istat==0)
            drcharge = 0._dp
          END IF
          NULLIFY(n_list)
          CALL get_qs_env(qs_env=qs_env,sab_scp=n_list)
          CALL neighbor_list_iterator_create(nl_iterator,n_list)
          DO WHILE (neighbor_list_iterate(nl_iterator)==0)
             CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)

             dr2 = SUM(rij(:)**2)
             dr = SQRT(dr2)
             IF (dr > 1.e-10) THEN
                CALL get_qs_kind(qs_kind_set(ikind), scptb_parameter=scptb_kind)
                CALL get_scptb_parameter(scptb_kind,lmaxscp=lmaxscp,ag=alpha)
                IF (lmaxscp < 0) CYCLE
                noa = (alpha/pi)**1.5_dp
                CALL get_qs_kind(qs_kind_set(jkind), scptb_parameter=scptb_kind)
                CALL get_scptb_parameter(scptb_kind,lmaxscp=lmaxscp,ag=beta)
                IF (lmaxscp < 0) CYCLE
                eta = alpha*beta/(alpha+beta)
                f0 = 2.0_dp*SQRT(pi**5/(alpha+beta))/(alpha*beta)
                nob = (beta/pi)**1.5_dp
                t = eta*dr2
                CALL fgamma_0(0,t,f)
                ff = noa*nob*f(0)*f0
                ehartree = ehartree + mcharge(iatom)*mcharge(jatom)*ff
                ehartree = ehartree - mcharge(iatom)*mcharge(jatom)/dr
                IF (.NOT.just_energy .OR. calculate_forces) THEN
                   qcharge(iatom) = qcharge(iatom) + mcharge(jatom)*ff
                   qcharge(jatom) = qcharge(jatom) + mcharge(iatom)*ff
                   qcharge(iatom) = qcharge(iatom) - mcharge(jatom)/dr
                   qcharge(jatom) = qcharge(jatom) - mcharge(iatom)/dr
                END IF
                IF ( calculate_forces ) THEN
                   ff = noa*nob*f0*(EXP(-t)-f(0))/dr2
                   drcharge(:,iatom) = drcharge(:,iatom) + ff*rij(:)*mcharge(iatom)*mcharge(jatom)
                   drcharge(:,jatom) = drcharge(:,jatom) - ff*rij(:)*mcharge(iatom)*mcharge(jatom)
                   drcharge(:,iatom) = drcharge(:,iatom) + rij(:)/(dr*dr2)*mcharge(iatom)*mcharge(jatom)
                   drcharge(:,jatom) = drcharge(:,jatom) - rij(:)/(dr*dr2)*mcharge(iatom)*mcharge(jatom)
                   IF (use_virial) THEN
                      fij(:) = ff*rij(:)*mcharge(iatom)*mcharge(jatom) + &
                               rij(:)/(dr*dr2)*mcharge(iatom)*mcharge(jatom)
                      CALL virial_pair_force ( virial%pv_virial, -1._dp, fij, rij)
                   END IF
                END IF

             END IF

          END DO
          CALL neighbor_list_iterator_release(nl_iterator)
          CALL mp_sum(ehartree,para_env%group)
          energy%hartree = energy%hartree + ehartree

          IF(.NOT.just_energy.OR.calculate_forces) THEN
            IF ( calculate_forces ) THEN
               CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
               IF (SIZE(matrix_p,1) == 2) THEN
                  DO ic=1,nimg
                     CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                          alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
                  END DO
               END IF
            END IF
            ! calculate KS matrix
            CALL mp_sum(qcharge,para_env%group)

            NULLIFY(n_list)
            CALL get_qs_env(qs_env=qs_env,sab_orb=n_list)
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
                  CPASSERT(ic > 0)
               END IF

               dq = qcharge(iatom)+qcharge(jatom)
               IF(.NOT.just_energy) THEN
                  CALL cp_dbcsr_get_block_p(matrix=matrix_s(1,ic)%matrix,&
                       row=irow,col=icol,block=sblock,found=found)
                  CPASSERT(found)
                  DO is=1,nsize
                     NULLIFY(ksblock)
                     CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is,ic)%matrix,&
                          row=irow,col=icol,block=ksblock,found=found)
                     CPASSERT(found)
                     ksblock = ksblock - 0.5_dp*dq*sblock
                  END DO
               END IF
               IF ( calculate_forces .AND. iatom /= jatom ) THEN
                  ikind = kind_of(iatom)
                  iat   = atom_of_kind(iatom)
                  jkind = kind_of(jatom)
                  jat   = atom_of_kind(jatom)
                  NULLIFY(pblock)
                  CALL cp_dbcsr_get_block_p(matrix=matrix_p(1,ic)%matrix,&
                       row=irow,col=icol,block=pblock,found=found)
                  CPASSERT(found)
                  DO i=1,3
                     NULLIFY(dsblock)
                     CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i,ic)%matrix,&
                          row=irow,col=icol,block=dsblock,found=found)
                     CPASSERT(found)
                     fij(i) = dq*SUM(pblock*dsblock)
                     force(ikind)%rho_elec(i,iat) = force(ikind)%rho_elec(i,iat) - fij(i)
                     force(jkind)%rho_elec(i,jat) = force(jkind)%rho_elec(i,jat) + fij(i)
                  END DO
                  IF (use_virial) THEN
                     CALL virial_pair_force ( virial%pv_virial, -1._dp, fij, rij)
                  END IF
               END IF

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)

            IF ( calculate_forces ) THEN
               DO iatom=1,natom
                  ikind  = kind_of(iatom)
                  iat = atom_of_kind(iatom)
                  force(ikind)%rho_elec(1,iat) = force(ikind)%rho_elec(1,iat) - drcharge(1,iatom)
                  force(ikind)%rho_elec(2,iat) = force(ikind)%rho_elec(2,iat) - drcharge(2,iatom)
                  force(ikind)%rho_elec(3,iat) = force(ikind)%rho_elec(3,iat) - drcharge(3,iatom)
               END DO
               IF (SIZE(matrix_p,1) == 2) THEN
                  DO ic=1,nimg
                     CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                                       alpha_scalar=1.0_dp,beta_scalar=-1.0_dp)
                  END DO
               END IF
            END IF
          END IF

          DEALLOCATE(gmcharge,STAT=istat)
          CPASSERT(istat==0)
          IF(.NOT.just_energy .OR. calculate_forces) THEN
            DEALLOCATE(qcharge,STAT=istat)
            CPASSERT(istat==0)
          END IF
          IF ( calculate_forces ) THEN
            DEALLOCATE(drcharge,STAT=istat)
            CPASSERT(istat==0)
            DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
            CPASSERT(istat==0)
          END IF
       ELSE
         ! direct sum
         IF (nimg > 1) THEN
            CALL cp_unimplemented_error(fromWhere=routineP, &
                 message="SCPTB direct sum electrostatic cannot be used with k-points.", &
                 error_level=cp_failure_level)
         END IF
         IF(.NOT.just_energy .OR. calculate_forces) THEN
            ALLOCATE(qcharge(natom),STAT=istat)
            CPASSERT(istat==0)
            qcharge = 0._dp
         END IF
         IF (calculate_forces) THEN
            ALLOCATE(drcharge(3,natom),STAT=istat)
            CPASSERT(istat==0)
            drcharge = 0._dp
         END IF
         nkind = SIZE(atomic_kind_set)
         DO ikind=1,nkind
           CALL get_atomic_kind(atomic_kind_set(ikind),natom=nai,atom_list=atomi_list)
           CALL get_qs_kind(qs_kind_set(ikind),scptb_parameter=scptb_kind)
           CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=alpha)
           IF (.NOT.defined) CYCLE
           IF (lmaxscp < 0) CYCLE
           noa = (alpha/pi)**1.5_dp
           DO jkind=1,ikind
             CALL get_atomic_kind(atomic_kind_set(jkind),natom=naj,atom_list=atomj_list)
             CALL get_qs_kind(qs_kind_set(jkind), scptb_parameter=scptb_kind)
             CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=beta)
             eta = alpha*beta/(alpha+beta)
             f0 = 2.0_dp*SQRT(pi**5/(alpha+beta))/(alpha*beta)
             IF (.NOT.defined) CYCLE
             IF (lmaxscp < 0) CYCLE
             nob = (beta/pi)**1.5_dp
             DO ia=1,nai
               iatom=atomi_list(ia)
               DO ja=1,local_particles%n_el(jkind)
                 jatom=local_particles%list(jkind)%array(ja)
                 IF(ikind==jkind .AND. jatom >= iatom) CYCLE
                 rij = particle_set(iatom)%r - particle_set(jatom)%r
                 rij = pbc(rij,cell)
                 dr2 = SUM(rij(:)**2)
                 t = eta*dr2
                 CALL fgamma_0(0,t,f)
                 ff = noa*nob*f(0)*f0
                 energy%hartree = energy%hartree + mcharge(iatom)*mcharge(jatom)*ff
                 IF (.NOT.just_energy .OR. calculate_forces) THEN
                    qcharge(iatom) = qcharge(iatom) + mcharge(jatom)*ff
                    qcharge(jatom) = qcharge(jatom) + mcharge(iatom)*ff
                 END IF
                 IF ( calculate_forces ) THEN
                    ff = noa*nob*f0*(EXP(-t)-f(0))/dr2
                    drcharge(:,iatom) = drcharge(:,iatom) + ff*rij(:)*mcharge(iatom)*mcharge(jatom)
                    drcharge(:,jatom) = drcharge(:,jatom) - ff*rij(:)*mcharge(iatom)*mcharge(jatom)
                 END IF
               END DO
             END DO
           END DO
         END DO
         CALL mp_sum(energy%hartree,para_env%group)
         IF(.NOT.just_energy.OR.calculate_forces) THEN
            CALL mp_sum(qcharge,para_env%group)
            ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
            CPASSERT(istat==0)
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                     kind_of=kind_of, atom_of_kind=atom_of_kind)
            IF ( calculate_forces ) THEN
               CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
               IF (SIZE(matrix_p,1) == 2) THEN
                  DO ic=1,nimg
                     CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                          alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
                  END DO
               END IF
            END IF

            ! calculate KS matrix
            nsize = SIZE(ks_matrix)
            CALL cp_dbcsr_iterator_start(iter,matrix_s(1,1)%matrix)
            DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
               CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
               dq = qcharge(iatom)+qcharge(jatom)
               IF(.NOT.just_energy) THEN
                  DO is=1,nsize
                     NULLIFY(ksblock)
                     CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is,1)%matrix,&
                          row=iatom,col=jatom,block=ksblock,found=found)
                     CPASSERT(found)
                     ksblock = ksblock - 0.5_dp*dq*sblock
                  END DO
               END IF
               IF ( calculate_forces ) THEN
                  ikind = kind_of(iatom)
                  iat   = atom_of_kind(iatom)
                  jkind = kind_of(jatom)
                  jat   = atom_of_kind(jatom)
                  NULLIFY(pblock)
                  CALL cp_dbcsr_get_block_p(matrix=matrix_p(1,1)%matrix,&
                       row=iatom,col=jatom,block=pblock,found=found)
                  CPASSERT(found)
                  DO i=1,3
                     NULLIFY(dsblock)
                     CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i,1)%matrix,&
                          row=iatom,col=jatom,block=dsblock,found=found)
                     CPASSERT(found)
                     fi = dq*SUM(pblock*dsblock)
                     force(ikind)%rho_elec(i,iat) = force(ikind)%rho_elec(i,iat) - fi
                     force(jkind)%rho_elec(i,jat) = force(jkind)%rho_elec(i,jat) + fi
                  END DO
               END IF
            END DO
            CALL cp_dbcsr_iterator_stop(iter)

            IF ( calculate_forces ) THEN
               DO iatom=1,natom
                  ikind  = kind_of(iatom)
                  iat = atom_of_kind(iatom)
                  force(ikind)%rho_elec(1,iat) = force(ikind)%rho_elec(1,iat) + drcharge(1,iatom)
                  force(ikind)%rho_elec(2,iat) = force(ikind)%rho_elec(2,iat) + drcharge(2,iatom)
                  force(ikind)%rho_elec(3,iat) = force(ikind)%rho_elec(3,iat) + drcharge(3,iatom)
               END DO
               DEALLOCATE(drcharge,STAT=istat)
               CPASSERT(istat==0)
               IF (SIZE(matrix_p,1) == 2) THEN
                  DO ic=1,nimg
                     CALL cp_dbcsr_add(matrix_p(1,ic)%matrix,matrix_p(2,ic)%matrix,&
                                       alpha_scalar=1.0_dp,beta_scalar=-1.0_dp)
                  END DO
               END IF
            END IF
 
            DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
            CPASSERT(istat==0)
            DEALLOCATE(qcharge,STAT=istat)
            CPASSERT(istat==0)
         END IF

       END IF
    END IF

  END SUBROUTINE tb_coulomb

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrixkp ...
!> \param matrixkp_s ...
!> \param rho ...
!> \param mcharge ...
!> \param calculate_forces ...
!> \param just_energy ...
! *****************************************************************************
  SUBROUTINE scp_coulomb(qs_env,ks_matrixkp,matrixkp_s,rho,mcharge,&
             calculate_forces,just_energy)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: ks_matrixkp, matrixkp_s
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(inout)                          :: mcharge
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_coulomb', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: ndiis = 5
    REAL(KIND=dp), PARAMETER                 :: epsrel = 1.e-8_dp

    INTEGER :: blk, i, iat, iatom, ibas, ikind, is, istat, j, jat, jatom, &
      jkind, kint, kmax, l, lmaxscp, mdiis, mpos, nat, natom, nkind, norb, &
      nsize
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, lscp, &
                                                natoms, nbasis
    LOGICAL                                  :: defined, do_kpoints, found
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: adef
    REAL(KIND=dp)                            :: ak, alpha, dq, dqv, fi, rhom1
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: charges, zeta
    REAL(KIND=dp), DIMENSION(1:3)            :: pol
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsblock, ksblock, pblock, &
                                                sblock
    REAL(KIND=dp), DIMENSION(ndiis+1, 1)     :: vdiis
    REAL(KIND=dp), &
      DIMENSION(ndiis+1, ndiis+1)            :: bdiis
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p, matrix_s
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scp_vector_p_type), &
      DIMENSION(ndiis)                       :: cres, cval
    TYPE(scp_vector_type), POINTER           :: cp, cpol, cself, cvec, cw
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
         do_kpoints=do_kpoints,&
         energy=energy,force=force,para_env=para_env)

    IF(do_kpoints) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="SCPTB not implemented with k-points.", &
            error_level=cp_failure_level)
    END IF

    NULLIFY(ks_matrix,matrix_s)
    ks_matrix => ks_matrixkp(:,1)
    matrix_s => matrixkp_s(:,1)

    energy%hartree = 0._dp

    natom = SIZE(mcharge,1)
    ALLOCATE(charges(natom),stat=istat)
    CPASSERT(istat==0)
    dqv = 1._dp/SQRT(4._dp*pi)
    charges(:) = mcharge*dqv
    scptb_control => dft_control%qs_control%scptb_control

    NULLIFY(cvec,cpol,cself,cp,cw)

    ! Create scp vectors
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set)
    nkind = SIZE(atomic_kind_set)
    ALLOCATE(nbasis(nkind),natoms(nkind),lscp(nkind),zeta(nkind),adef(nkind),stat=istat)
    CPASSERT(istat==0)
    DO ikind=1,nkind
       CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
       CALL get_qs_kind(qs_kind_set(ikind), scptb_parameter=scptb_kind)
       CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=alpha)
       IF(defined) THEN
         norb = (lmaxscp+1)**2
       ELSE
         norb = 0
       END IF
       natoms(ikind) = nat
       nbasis(ikind) = norb
       lscp(ikind)   = lmaxscp
       zeta(ikind)   = alpha
       adef(ikind)   = defined
    END DO
    NULLIFY(cvec,cself,cpol,cp,cw)
    CALL scp_vector_create(cvec, nkind, natoms, nbasis)
    CALL scp_vector_set(cvec, 0._dp)
    CALL scp_vector_create(cself, nkind, natoms, nbasis)
    CALL scp_vector_set(cself, 0._dp)
    CALL scp_vector_create(cpol, nkind, natoms, nbasis)
    CALL scp_vector_set(cpol, 0._dp)
    ! set vector of multipoles
    CALL scp_set_charge(cvec,charges,adef,natoms,atomic_kind_set)
    DO ikind=1,nkind
       IF(adef(ikind)) THEN
         CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
         CALL get_qs_kind(qs_kind_set(ikind),scptb_parameter=scptb_kind)
         CALL get_scptb_parameter(scptb_kind,pol=pol)
         alpha = zeta(ikind)
         DO iat=1,natoms(ikind)
            DO ibas=1,nbasis(ikind)
               l = indso(1,ibas)
               IF(l > 0) THEN
                  cpol%vector(ikind)%vmat(ibas,iat) = 1._dp/pol(l)
               END IF
               cself%vector(ikind)%vmat(ibas,iat) = &
                  4._dp/(2*l+1)/dfac(2*l+1)*SQRT(2._dp*pi)*alpha**(l+0.5_dp)
            END DO
         END DO
       END IF
    END DO

    ! Calclulate scp vector
    DO iat=1,ndiis
       NULLIFY(cval(iat)%vec,cres(iat)%vec)
       CALL scp_vector_create(cval(iat)%vec, nkind, natoms, nbasis)
       CALL scp_vector_create(cres(iat)%vec, nkind, natoms, nbasis)
    END DO
    CALL scp_vector_create(cp, nkind, natoms, nbasis)
    CALL scp_vector_create(cw, nkind, natoms, nbasis)

    kmax = MIN(SUM(nbasis),100) * 10
    mdiis = 0
    DO kint=1,kmax
       mdiis = MIN(mdiis + 1, ndiis)
       mpos = MOD(kint-1,ndiis) + 1
       CALL scp_vector_set(cw, 0._dp)
       CALL apply_scp_vector(cvec,cw,cself,cpol,cp,qs_env,scptb_control%do_ewald,&
                             calculate_forces=.FALSE.)
       rhom1 = scp_dot(cw,cw)
       IF ( SQRT(rhom1) < epsrel ) EXIT
       CALL scp_vector_copy(cvec, cval(mpos)%vec)
       CALL scp_vector_copy(cw, cres(mpos)%vec)
       IF (mdiis > 100) THEN
          vdiis = 0._dp
          vdiis(mdiis+1,1) = 1._dp
          bdiis = 0._dp
          bdiis(mdiis+1,1:mdiis) = 1._dp
          bdiis(1:mdiis,mdiis+1) = 1._dp
          DO i=1,mdiis
             DO j=1,i
                bdiis(i,j) = scp_dot(cres(i)%vec,cres(j)%vec)
                bdiis(j,i) = bdiis(i,j)
             END DO
          END DO
          CALL solve_system ( bdiis, mdiis+1, vdiis )
          CALL scp_vector_set(cvec, 0._dp)
          DO i=1,mdiis
             CALL scp_vector_add(vdiis(i,1), cval(i)%vec, cvec)
          END DO
       ELSE
          CALL scp_vector_add(-0.1_dp, cw, cvec)
       END IF
       CALL scp_set_charge(cvec,charges,adef,natoms,atomic_kind_set)
    END DO
    DO iat=1,ndiis
       CALL scp_vector_release(cval(iat)%vec)
       CALL scp_vector_release(cres(iat)%vec)
    END DO
    ! SCP energy
    CALL scp_vector_set(cw, 0._dp)
    IF (calculate_forces) THEN
       CALL apply_scp_vector(cvec,cw,cself,cpol,cp,qs_env,scptb_control%do_ewald,&
                             calculate_forces=.TRUE.)
    ELSE
       CALL apply_scp_vector(cvec,cw,cself,cpol,cp,qs_env,scptb_control%do_ewald,&
                             calculate_forces=.FALSE.)
    END IF
    CALL scp_vector_dot(ak, cw, cvec)
    energy%hartree = 0.5_dp*ak

    IF (.NOT.just_energy) THEN

       ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
       CPASSERT(istat==0)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)
       IF ( calculate_forces ) THEN
          CALL qs_rho_get(rho, rho_ao=matrix_p)
          IF (SIZE(matrix_p) == 2) THEN
             CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
                  alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
          END IF
       END IF

       ! calculate KS matrix
       nsize = SIZE(ks_matrix)
       CALL cp_dbcsr_iterator_start(iter,matrix_s(1)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
          ikind = kind_of(iatom)
          iat   = atom_of_kind(iatom)
          jkind = kind_of(jatom)
          jat   = atom_of_kind(jatom)
          dq = 0.5_dp*(cw%vector(ikind)%vmat(1,iat) + cw%vector(jkind)%vmat(1,jat))/SQRT(fourpi)
          DO is=1,nsize
             NULLIFY(ksblock)
             CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is)%matrix,&
                  row=iatom,col=jatom,block=ksblock,found=found)
             ksblock = ksblock - dq*sblock
          END DO
          IF ( calculate_forces .AND. iatom /= jatom ) THEN
             NULLIFY(pblock)
             CALL cp_dbcsr_get_block_p(matrix=matrix_p(1)%matrix,&
                  row=iatom,col=jatom,block=pblock,found=found)
             CPASSERT(found)
             dq = (cw%vector(ikind)%vmat(1,iat) + &
                   cw%vector(ikind)%vmat(1,iat) * cself%vector(ikind)%vmat(1,iat) + &
                   cw%vector(jkind)%vmat(1,jat) + &
                   cw%vector(jkind)%vmat(1,jat) * cself%vector(jkind)%vmat(1,jat))/SQRT(fourpi)
             DO i=1,3
                NULLIFY(dsblock)
                CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i)%matrix,&
                     row=iatom,col=jatom,block=dsblock,found=found)
                CPASSERT(found)
                fi = dq*SUM(pblock*dsblock)
!deb what about these forces?
!               force(ikind)%rho_elec(i,iat) = force(ikind)%rho_elec(i,iat) + fi
!               force(jkind)%rho_elec(i,jat) = force(jkind)%rho_elec(i,jat) - fi
             END DO
          END IF
       END DO
       CALL cp_dbcsr_iterator_stop(iter)

       IF (calculate_forces) THEN
          IF (SIZE(matrix_p) == 2) THEN
             CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
                           alpha_scalar=1.0_dp,beta_scalar=-1.0_dp)
          END IF
       END IF

    END IF

    CALL scp_vector_release(cself)
    CALL scp_vector_release(cvec)
    CALL scp_vector_release(cpol)
    CALL scp_vector_release(cw)
    CALL scp_vector_release(cp)
    DEALLOCATE(nbasis,natoms,lscp,zeta,adef,charges,stat=istat)
    CPASSERT(istat==0)

  END SUBROUTINE scp_coulomb

! *****************************************************************************
!> \brief ...
!> \param vec1 ...
!> \param vec2 ...
!> \retval res ...
! *****************************************************************************
  FUNCTION scp_dot(vec1,vec2) RESULT(res)
    TYPE(scp_vector_type)                    :: vec1, vec2
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: i

! dot product without l=0 term

    res = 0._dp
    DO i=1,SIZE(vec1%vector)
      res = res + SUM(vec1%vector(i)%vmat(2:,:)*vec2%vector(i)%vmat(2:,:))
    END DO

  END FUNCTION scp_dot

! *****************************************************************************
!> \brief ...
!> \param vec ...
!> \param charge ...
!> \param adef ...
!> \param natoms ...
!> \param atomic_kind_set ...
! *****************************************************************************
  SUBROUTINE scp_set_charge(vec,charge,adef,natoms,atomic_kind_set)
    TYPE(scp_vector_type)                    :: vec
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: charge
    LOGICAL, DIMENSION(:), INTENT(IN)        :: adef
    INTEGER, DIMENSION(:), INTENT(IN)        :: natoms
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    INTEGER                                  :: iat, iatom, ikind, nkind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    nkind = SIZE(adef)
    DO ikind=1,nkind
       IF(adef(ikind)) THEN
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,atom_list=atom_list)
         DO iat=1,natoms(ikind)
            iatom = atom_list(iat)
            vec%vector(ikind)%vmat(1,iat) = charge(iatom)
         END DO
       END IF
    END DO

  END SUBROUTINE scp_set_charge

! *****************************************************************************
!> \brief ...
!> \param vec ...
!> \param message ...
! *****************************************************************************
  SUBROUTINE scp_print(vec,message)
    TYPE(scp_vector_type)                    :: vec
    CHARACTER(LEN=*), INTENT(IN)             :: message

    INTEGER                                  :: iat, iba, ikind, natom, nbas, &
                                                nkind

    WRITE(*,*) message
    nkind = SIZE(vec%vector)
    DO ikind=1,nkind
       nbas = SIZE(vec%vector(ikind)%vmat,1)
       natom = SIZE(vec%vector(ikind)%vmat,2)
       DO iat=1,natom
          DO iba=1,nbas
             WRITE(*,*) "Kind=",ikind,"  Atom=",iat," Pol=",vec%vector(ikind)%vmat(iba,iat) 
          END DO
       END DO
    END DO

  END SUBROUTINE scp_print
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param cout ...
!> \param cself ...
!> \param cpol ...
!> \param cdum ...
!> \param qs_env ...
!> \param do_ewald ...
!> \param calculate_forces ...
! *****************************************************************************
  SUBROUTINE apply_scp_vector(cin,cout,cself,cpol,cdum,qs_env,do_ewald,calculate_forces)

    TYPE(scp_vector_type), POINTER           :: cin, cout, cself, cpol, cdum
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: do_ewald, calculate_forces

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_scp_vector', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: dr = 0.001_dp

    INTEGER                                  :: i, ia, iatom, ikind, istat, &
                                                j, ja, jatom, jkind, li, lj, &
                                                lmaxscp, nai, naj, ni, nj, &
                                                nkind
    INTEGER, DIMENSION(:), POINTER           :: atomi_list, atomj_list
    LOGICAL                                  :: defined, failure, use_virial
    REAL(KIND=dp)                            :: alpha, beta, ehartree
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: aintegral, maint, paint
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: daintegral
    REAL(KIND=dp), DIMENSION(16)             :: normi, normj
    REAL(KIND=dp), DIMENSION(3)              :: drij, fij, rij
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cii, cij, coi, coj
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: scp_pot, scp_rho, scp_rho_g
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,force=force,virial=virial,para_env=para_env)
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    IF (do_ewald) THEN
       CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,&
                       pw_pools=pw_pools,poisson_env=poisson_env)
       CALL pw_pool_create_pw(auxbas_pw_pool,scp_rho%pw,use_data=REALDATA3D,in_space=REALSPACE)
       CALL pw_pool_create_pw(auxbas_pw_pool,scp_pot%pw,use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE)
       CALL pw_pool_create_pw(auxbas_pw_pool,scp_rho_g%pw,use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE)

       CALL calculate_scp_charge(scp_rho,qs_env,cin)
       CALL pw_transfer(scp_rho%pw, scp_rho_g%pw)

       ! Getting the Hartree energy and Hartree potential.  Also getting the stress tensor
       ! from the Hartree term if needed.
       IF (use_virial .AND. calculate_forces) THEN
          h_stress(:,:) = 0.0_dp
          CALL pw_poisson_solve(poisson_env,scp_rho_g%pw,ehartree,scp_pot%pw,h_stress=h_stress)
          virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe,dp)
       ELSE
          CALL pw_poisson_solve(poisson_env,scp_rho_g%pw,ehartree,scp_pot%pw)
       END IF
       CALL pw_pool_give_back_pw(auxbas_pw_pool,scp_rho_g%pw)

       CALL pw_transfer(scp_pot%pw,scp_rho%pw)
       CALL integrate_scp_rspace(scp_rho,qs_env,cout,calculate_forces)

       CALL pw_pool_give_back_pw(auxbas_pw_pool,scp_rho%pw)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,scp_pot%pw)
    ELSE
       ! direct sum
       NULLIFY(atomic_kind_set,qs_kind_set,particle_set,cell)
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
                       particle_set=particle_set,cell=cell)
       nkind = SIZE(atomic_kind_set)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind),natom=nai,atom_list=atomi_list)
          CALL get_qs_kind(qs_kind_set(ikind), scptb_parameter=scptb_kind)
          CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=alpha)
          IF (.NOT.defined) CYCLE
          IF (lmaxscp < 0) CYCLE
          li = lmaxscp
          ni = (li+1)**2
          coi => cout%vector(ikind)%vmat
          cii => cin%vector(ikind)%vmat
          CALL calc_norm(normi,li,alpha)
          DO jkind=1,ikind
             CALL get_atomic_kind(atomic_kind_set(jkind), natom=naj,atom_list=atomj_list)
             CALL get_qs_kind(qs_kind_set(jkind), scptb_parameter=scptb_kind)
             CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=beta)
             IF (.NOT.defined) CYCLE
             IF (lmaxscp < 0) CYCLE
             lj = lmaxscp
             nj = (lj+1)**2
             coj => cout%vector(jkind)%vmat
             cij => cin%vector(jkind)%vmat
             ALLOCATE (aintegral(ni,nj),STAT=istat)
             CPASSERT(istat==0)
             IF(calculate_forces) THEN
                ALLOCATE (daintegral(ni,nj,3),paint(ni,nj),maint(ni,nj),STAT=istat)
                CPASSERT(istat==0)
             END IF
             CALL calc_norm(normj,lj,beta)
             DO iatom=1,nai
                ia = atomic_kind_set(ikind)%atom_list(iatom)
                DO jatom=1,naj
                   ja = atomic_kind_set(jkind)%atom_list(jatom)
                   IF(ikind==jkind .AND. jatom > iatom) CYCLE
                   rij = particle_set(ia)%r - particle_set(ja)%r
                   rij = pbc(rij,cell)
                   CALL calc_int(aintegral,rij,alpha,li,beta,lj)
                   ! add normalization constants to integrals
                   DO j=1,nj
                      DO i=1,ni
                         aintegral(i,j) = aintegral(i,j)*normi(i)*normj(j)
                      END DO
                   END DO
                   ! add contribution to vector
                   coi(1:ni,iatom) = coi(1:ni,iatom) + MATMUL(aintegral(1:ni,1:nj),cij(1:nj,jatom))
                   IF(ikind==jkind .AND. jatom == iatom) CYCLE
                   coj(1:nj,jatom) = coj(1:nj,jatom) + MATMUL(cii(1:ni,iatom),aintegral(1:ni,1:nj))
                   IF(calculate_forces) THEN
                      DO i=1,3
                         drij = rij
                         drij(i) = drij(i) + dr
                         CALL calc_int(paint,drij,alpha,li,beta,lj)
                         drij(i) = drij(i) - 2.0_dp*dr
                         CALL calc_int(maint,drij,alpha,li,beta,lj)
                         daintegral(:,:,i) = (paint - maint)/(2.0_dp*dr)
                      END DO
                      DO j=1,nj
                         DO i=1,ni
                            daintegral(i,j,1) = daintegral(i,j,1)*normi(i)*normj(j)
                            daintegral(i,j,2) = daintegral(i,j,2)*normi(i)*normj(j)
                            daintegral(i,j,3) = daintegral(i,j,3)*normi(i)*normj(j)
                         END DO
                      END DO
                      fij(1) = SUM(cii(1:ni,iatom) * MATMUL(daintegral(1:ni,1:nj,1),cij(1:nj,jatom)))
                      fij(2) = SUM(cii(1:ni,iatom) * MATMUL(daintegral(1:ni,1:nj,2),cij(1:nj,jatom)))
                      fij(3) = SUM(cii(1:ni,iatom) * MATMUL(daintegral(1:ni,1:nj,3),cij(1:nj,jatom)))
                      force(ikind)%rho_elec(:,iatom) = force(ikind)%rho_elec(:,iatom) + fij
                      force(jkind)%rho_elec(:,jatom) = force(jkind)%rho_elec(:,jatom) - fij
                   END IF
                END DO
             END DO
             DEALLOCATE (aintegral,STAT=istat)
             CPASSERT(istat==0)
             IF(calculate_forces) THEN
                DEALLOCATE (daintegral,paint,maint,STAT=istat)
                CPASSERT(istat==0)
             END IF
          END DO
       END DO

    END IF

    ! remove self interaction
    CALL scp_vector_copy(cself, cdum)
    CALL scp_vector_mult(cin, cdum)
    CALL scp_vector_add(-1._dp, cdum, cout)
    ! add scp polarization
    CALL scp_vector_copy(cpol, cdum)
    CALL scp_vector_mult(cin, cdum)
    CALL scp_vector_add(1._dp, cdum, cout)

  END SUBROUTINE apply_scp_vector

! *****************************************************************************
!> \brief ...
!> \param aintegral ...
!> \param rij ...
!> \param alpha ...
!> \param li ...
!> \param beta ...
!> \param lj ...
! *****************************************************************************
  SUBROUTINE calc_int(aintegral,rij,alpha,li,beta,lj)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: aintegral
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), INTENT(IN)                :: alpha
    INTEGER, INTENT(IN)                      :: li
    REAL(KIND=dp), INTENT(IN)                :: beta
    INTEGER, INTENT(IN)                      :: lj

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, mi, mj, ni, nj, nl
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dr2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vij
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: v
    REAL(KIND=dp), DIMENSION(0:9)            :: f
    REAL(KIND=dp), DIMENSION(1)              :: zi, zj

    failure = .FALSE.
    dr2 = SUM(rij(:)**2)
    zi(1) = alpha
    zj(1) = beta
    nl = li + lj + 1
    ni = ncoset(li)
    nj = ncoset(lj)
    ALLOCATE (vij(ni,nj),v(ni,nj,nl),STAT=istat)
    CPASSERT(istat==0)

    ! cartesian integrals
    CALL coulomb2_new(li,1,zi,0,lj,1,zj,0,rij,dr2,vij,v,f)
    ! transform to spherical
    mi = nsoset(li)
    mj = nsoset(lj)
    v(1:ni,1:mj,1) = MATMUL(vij(1:ni,1:nj),TRANSPOSE(c2s_tramat(1:mj,1:nj)))
    aintegral(1:mi,1:mj) = MATMUL(c2s_tramat(1:mi,1:ni),v(1:ni,1:mj,1))

    DEALLOCATE (vij,v,STAT=istat)
    CPASSERT(istat==0)

  END SUBROUTINE calc_int

! *****************************************************************************
!> \brief ...
!> \param anorm ...
!> \param lm ...
!> \param alpha ...
! *****************************************************************************
  SUBROUTINE calc_norm(anorm,lm,alpha)

    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: anorm
    INTEGER, INTENT(IN)                      :: lm
    REAL(KIND=dp), INTENT(IN)                :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, il, l

    anorm = 0.0_dp
    i = 0
    DO l=0,lm
       DO il=1,2*l+1
          i = i + 1
          anorm(i) = 2.0_dp**(l+2)*oorootpi/dfac(2*l+1) * alpha**(l+1.5_dp)
       END DO
    END DO

  END SUBROUTINE calc_norm

! *****************************************************************************

END MODULE scptb_ks_matrix

