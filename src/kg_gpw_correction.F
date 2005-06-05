!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_gpw_correction
!!
!!   NAME
!!     kg_gpw_correction
!!
!!     
!!   AUTHOR
!!     MI (20.12.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

 
MODULE kg_gpw_correction

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_init,&
                                             coeff_zero,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type, & 
                                             coeff_copy
  USE cp_control_types,                ONLY: dft_control_type
  USE input_constants,                 ONLY: sic_none
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pool_create,&
                                             rs_pools_dealloc
  USE cube_utils,                      ONLY: destroy_cube_info,&
                                             init_cube_info,&
                                             cube_info_type
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gaussian_gridlevel,&
                                             init_gaussian_gridlevel
  USE kg_gpw_collocate_den,            ONLY: integrate_mol_potential
  USE kg_gpw_fm_mol_types,             ONLY: get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             fm_mol_blocks_type,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type
  USE kg_gpw_pw_env_types,             ONLY: get_molbox_env,&
                                             get_rho_mol_block,&
                                             kg_molbox_env_type,&
                                             kg_sub_pw_env_type,&
                                             rho_mol_blocks_type
  USE kinds,                           ONLY: dp, dp_size, int_size
  USE l_utils,                         ONLY: init_l_info,&
                                             l_info_retain,&
                                             l_info_release,&
                                             l_info_type
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grids,                        ONLY: 
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             transfer_matrix, &
                                             real_matrix_type, &
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_exc_calc,&
                                             xc_vxc_pw_create1

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: module_name='kg_gpw_correction'

  PUBLIC :: kg_gpw_ekin_mol

!***
!****************************************************************************

CONTAINS

!****f* kg_gpw_correction/kg_gpw_ekin_mol
!!
!!   NAME
!!     kg_gpw_ekin_mol
!!
!!   FUNCTION
!!     loop over all the molecules where :
!!     the molecular densityis passed to the calculation of the ekin potential
!!     the ekin potential is integrated in real space
!!     the ks matrix is updated aggordingly (block by block)
!!     the forces coming from this contribution are calculated if required
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs environment
!!     - molbox_env : cell, internal positions, rho, pw_env for each molecule
!!     - fm_mol_set : other info about the molecule, maybe not neede
!!     - h : ks-matrix in sparse form
!!     - ekin_mol : total energy from this term (to be added to the total energy)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!    
!!
!!*** **********************************************************************

  SUBROUTINE kg_gpw_ekin_mol(qs_env,molbox_env,fm_mol_set,ks_global,p_global,&
                             ekin_mol,calculate_forces,just_energy,error)

    TYPE(qs_environment_type), POINTER           :: qs_env
    TYPE(kg_molbox_env_type), DIMENSION(:), &
      POINTER                                    :: molbox_env
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                    :: fm_mol_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                    :: ks_global
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                    :: p_global
    REAL(dp), INTENT(OUT)                        :: ekin_mol
    LOGICAL, INTENT(IN), OPTIONAL                :: calculate_forces, just_energy
        TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'kg_gpw_ekin_mol', & 
      routineP = module_name//':'//routine_name

    INTEGER :: handle, iat_mol, iat_kind, iatom, ikind, imol, imolecule_kind,&
               ispin, istat, nat_mol, natom, nmol_local, nmolecule_kind, nspins 
    INTEGER, DIMENSION(:), POINTER               :: i_atom, i_kind
    INTEGER, DIMENSION(:), ALLOCATABLE           :: atom_of_kind
    LOGICAL                                      :: failure, &
                                                    my_calculate_forces, &
                                                    my_just_energy
    CHARACTER(LEN=40), DIMENSION(3)              :: functionals_tmp
    LOGICAL, DIMENSION(3)                        :: gradient_functionals_tmp,&
                                                    crossterms_tmp
    REAL(dp)                                     :: ekin_imol
    REAL(dp), DIMENSION(:,:), POINTER            :: r_mbox
    REAL(dp), DIMENSION(:,:), POINTER            :: forces_mol

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                    :: atomic_kind_set
    TYPE(cell_type), POINTER                     :: cell_mol
    TYPE(coeff_type), DIMENSION(:), POINTER      :: vxc_rho, vxc_tau
    TYPE(cp_para_env_type), POINTER              :: para_env
    TYPE(dft_control_type), POINTER              :: dft_control
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                    :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER            :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER            :: fm_mol
    TYPE(kg_molbox_env_type), POINTER            :: molbox
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                    :: particle_set
    TYPE(pw_env_type) , POINTER                  :: pw_env_mol
    TYPE(pw_p_type), DIMENSION(:), POINTER       :: my_vxc_rho, my_vxc_tau, &
                                                    rho_g, rho_r, tau
    TYPE(pw_pool_type), POINTER                  :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), POINTER   :: force_global
    TYPE(qs_rho_type), POINTER                   :: rho_global
    TYPE(qs_rho_type), POINTER                   :: rho_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER                                    :: rho_mol_blocks
    TYPE(rho_mol_blocks_type), POINTER           :: rho_block

!------------------------------------------------------------------------------

    failure=.FALSE.

    CALL timeset(routine_name,"I",'',handle)

    ! initialize to zero
    ekin_mol = 0.0_dp

    NULLIFY(atomic_kind_set,particle_set,para_env,vxc_tau,tau)
    NULLIFY(dft_control,rho_global,rho_g,rho_r,my_vxc_rho,my_vxc_tau)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,rho=rho_global,para_env=para_env)

    ! Modify temporarely dft_control to calculate only the kinetic energy functional
    functionals_tmp(1:3) = dft_control%xc_control%functionals(1:3)
    gradient_functionals_tmp(1:3) = dft_control%xc_control%gradient_functionals(1:3)
    crossterms_tmp(1:3) = dft_control%xc_control%crossterms(1:3)

    dft_control%xc_control%functionals(1) = "NONE"
    dft_control%xc_control%functionals(2) = "NONE"
    dft_control%xc_control%gradient_functionals(1) = .FALSE.
    dft_control%xc_control%gradient_functionals(2) = .FALSE.
    dft_control%xc_control%crossterms(1) = .FALSE.
    dft_control%xc_control%crossterms(2) = .FALSE.

    ! Sic is not implemented
    IF(dft_control%sic_method_id /= sic_none) &
       CALL  stop_program(routine_name,module_name,__LINE__,"KG_GPW with SIC  not implemented")

    nspins = dft_control%nspins

    ! initial allocations
    ALLOCATE(rho_r(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (rho_global%tau_r_valid) THEN
       ALLOCATE(tau(nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF 

    ! for gradient corrected functional the density in g space might be useful
    IF ( rho_global%rho_g_valid ) THEN
       ALLOCATE(rho_g(nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    CPPrecondition(ASSOCIATED(molbox_env),cp_failure_level,routineP,error,failure)
    nmolecule_kind = SIZE(molbox_env,1)
    CPPrecondition(ASSOCIATED(fm_mol_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(fm_mol_set,1)==nmolecule_kind,cp_failure_level,routineP,error,failure)

    my_just_energy = .FALSE.
    IF(PRESENT(just_energy)) my_just_energy = just_energy
    my_calculate_forces = .FALSE.
    IF(PRESENT(calculate_forces)) my_calculate_forces = calculate_forces

    IF (my_calculate_forces) THEN
    ! Array of index within the kind 
    ! to associate the forces to the right position in forces array

      CALL get_qs_env(qs_env=qs_env,particle_set=particle_set)
      natom = SIZE(particle_set)
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL get_qs_env(qs_env=qs_env, force=force_global)

    END IF

    IF(.NOT. failure) THEN

      ! Loop over the molecule kinds
      DO imolecule_kind = 1,nmolecule_kind

        NULLIFY(molbox,pw_env_mol,cell_mol,rho_mol_blocks)
        molbox => molbox_env(imolecule_kind)
        CALL get_molbox_env(molbox_env=molbox,natom=nat_mol,cell_mol=cell_mol,&
                            pw_env_mol=pw_env_mol, rho_mol_blocks=rho_mol_blocks)
        nmol_local = SIZE(rho_mol_blocks,1)
        CALL pw_env_get(pw_env=pw_env_mol,auxbas_pw_pool=auxbas_pw_pool,error=error)

        NULLIFY(fm_mol,fm_mol_blocks)
        fm_mol => fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol,&
                               fm_mol_blocks=fm_mol_blocks)

        ! Prepare array forces
        IF(my_calculate_forces) THEN
          ALLOCATE(forces_mol(3,nat_mol),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ELSE
          NULLIFY(forces_mol)
        END IF

        ! Loop over the molecules of one kind
        DO imol = 1,nmol_local 

          NULLIFY(rho_block, rho_mol, r_mbox)
          rho_block => rho_mol_blocks(imol)
          CALL get_rho_mol_block(rho_block=rho_block, rho_mol=rho_mol,&
                                 r_in_molbox=r_mbox)

          NULLIFY(mol_block,i_atom,i_kind)
          mol_block => fm_mol_blocks(imol)
          CALL get_fm_mol_block(fm_mol_block = mol_block,&
                                index_atom = i_atom,&
                                index_kind = i_kind)

          DO ispin=1,nspins
            rho_r(ispin)%pw => rho_mol%rho_r(ispin)%pw
            IF(rho_mol%tau_r_valid) tau(ispin)%pw => rho_mol%tau_r(ispin)%pw
            IF(rho_mol%rho_g_valid) rho_g(ispin)%pw => rho_mol%rho_g(ispin)%pw
          END DO

          ekin_imol = 0.0_dp
          IF (my_just_energy) THEN
            ekin_imol = xc_exc_calc(rho_r=rho_r,tau=tau,&
                     rho_g=rho_g, dft_control=dft_control,&
                     cell=cell_mol, pw_pool=auxbas_pw_pool,&
                     error=error)
          ELSE
            CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                                   rho_g=rho_g,tau=tau,exc=ekin_imol,&
                                   dft_control=dft_control,&
                                   cell=cell_mol, pw_pool=auxbas_pw_pool,&
                                   error=error)
          END IF

          ekin_mol = ekin_mol + ekin_imol

          ! we have pw data for the xc, here we transfer to coeff
          ! pw -> coeff
          IF (ASSOCIATED(my_vxc_rho)) THEN
             ALLOCATE(vxc_rho(nspins),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,nspins
               CALL coeff_init(vxc_rho(ispin),pw=my_vxc_rho(ispin)%pw,&
                               error=error)
             END DO
             DEALLOCATE(my_vxc_rho,STAT=istat)
          END IF
          IF (ASSOCIATED(my_vxc_tau)) THEN
             ALLOCATE(vxc_tau(nspins),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,nspins
                CALL coeff_init(vxc_tau(ispin),pw=my_vxc_tau(ispin)%pw,&
                                error=error)
             END DO
             DEALLOCATE(my_vxc_tau,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

          ! Integrate the potential
          IF (.NOT. just_energy) THEN
             CPPostcondition(ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
             DO ispin=1,nspins

               vxc_rho(ispin)%pw%cr3d  =  vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol

               IF(my_calculate_forces) THEN
                 ! Initialize the forces
                 forces_mol(1:3,1:nat_mol) = 0.0_dp
                 CALL integrate_mol_potential(qs_env=qs_env,vxc_mol=vxc_rho(ispin),&
                                             matrix_p=p_global(ispin),&
                                             matrix_h=ks_global(ispin),&
                                             pw_env=pw_env_mol,&
                                             atom=i_atom, kind=i_kind, ratom=r_mbox,&
                                             forces_mol=forces_mol,&
                                             error=error)
               ELSE
                 CALL integrate_mol_potential(qs_env=qs_env,vxc_mol=vxc_rho(ispin),&
                                             matrix_p=p_global(ispin),&
                                             matrix_h=ks_global(ispin),&
                                             pw_env=pw_env_mol,&
                                             atom=i_atom, kind=i_kind, ratom=r_mbox,&
                                             error=error)
               END IF

               CALL pw_pool_give_back_coeff(auxbas_pw_pool,vxc_rho(ispin),error=error)

             END DO
             DEALLOCATE(vxc_rho,STAT=istat)
             CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

             IF (ASSOCIATED(vxc_tau)) THEN
               DO ispin=1,nspins

                 vxc_tau(ispin)%pw%cr3d  =  vxc_tau(ispin)%pw%cr3d*vxc_tau(ispin)%pw%pw_grid%dvol

                 IF(my_calculate_forces) THEN
                   CALL integrate_mol_potential(qs_env=qs_env,vxc_mol=vxc_tau(ispin),&
                                               matrix_p=p_global(ispin),&
                                               matrix_h=ks_global(ispin),&
                                               pw_env=pw_env_mol,&
                                               atom=i_atom, kind=i_kind, ratom=r_mbox,&
                                               forces_mol=forces_mol,&
                                               compute_tau=.TRUE.,error=error)
                 ELSE
                   CALL integrate_mol_potential(qs_env=qs_env,vxc_mol=vxc_rho(ispin),&
                                               matrix_p=p_global(ispin),&
                                               matrix_h=ks_global(ispin),&
                                               pw_env=pw_env_mol,&
                                               atom=i_atom, kind=i_kind, ratom=r_mbox,&
                                               compute_tau=.TRUE.,error=error)
                 END IF
                 CALL pw_pool_give_back_coeff(auxbas_pw_pool,vxc_tau(ispin),error=error)

               END DO
               DEALLOCATE(vxc_tau,STAT=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

             END IF
          END IF

          ! Copy the forces in the global array
          IF(my_calculate_forces) THEN

             DO iat_mol = 1,nat_mol
               ikind = i_kind(iat_mol)
               iatom = i_atom(iat_mol)
               iat_kind = atom_of_kind(iatom)
               ! Change sign: it is the correction
               force_global(ikind)%kg_gpw_ekin_mol(1:3,iat_kind) = -forces_mol(1:3,iat_mol)
             END DO

          END IF

        END DO  ! imol

        IF(ASSOCIATED(forces_mol)) THEN
          DEALLOCATE(forces_mol,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

      END DO ! imolecule_kind

       CALL mp_sum(ekin_mol,para_env%group)
    END IF ! failure


    ! Deallocate
    DEALLOCATE(rho_r,STAT=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    IF (ASSOCIATED(rho_g)) THEN
       DEALLOCATE(rho_g,STAT=istat)
       CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    END IF
    IF (ASSOCIATED(tau)) THEN
       DEALLOCATE(tau,STAT=istat)
       CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    END IF

    ! Set back dft_control%xc_control to the correct parameters for the functionals
    dft_control%xc_control%functionals(1:3) = functionals_tmp(1:3)
    dft_control%xc_control%gradient_functionals(1:3) = gradient_functionals_tmp(1:3)
    dft_control%xc_control%crossterms(1:3) = crossterms_tmp(1:3)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_gpw_ekin_mol


!***************************************************************************

END MODULE kg_gpw_correction
