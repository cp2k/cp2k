!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A collection of methods to treat the QM/MM electrostatic coupling
!> \par History
!>      5.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_gpw_energy
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_spline_utils,                 ONLY: pw_prolongate_s3
  USE cube_utils,                      ONLY: cube_info_type
  USE f77_blas
  USE input_constants,                 ONLY: do_par_atom,&
                                             do_qmmm_coulomb,&
                                             do_qmmm_gauss,&
                                             do_qmmm_none,&
                                             do_qmmm_swave,&
                                             spline3_nopbc_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE mm_collocate_potential,          ONLY: collocate_gf_rspace_NoPBC
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_se_energy,                  ONLY: build_se_qmmm_matrix
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_per_pot_p_type,&
                                             qmmm_per_pot_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qmmm_util,                       ONLY: spherical_cutoff_factor
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_gpw_energy'

  PUBLIC :: qmmm_el_coupling
  PUBLIC :: qmmm_elec_with_gaussian,&
            qmmm_elec_with_gaussian_LR, qmmm_elec_with_gaussian_LG
!***
CONTAINS

! *****************************************************************************
!> \brief Main Driver to compute the QM/MM Electrostatic Coupling
!> \param qs_env ...
!> \param qmmm_env ...
!> \param mm_particles ...
!> \param mm_cell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE  qmmm_el_coupling(qs_env,qmmm_env,mm_particles,mm_cell,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_el_coupling', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, iw2
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: input_section, &
                                                interp_section, print_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    IF (.NOT.failure) THEN
       NULLIFY(ks_qmmm_env_loc, pw_pools, pw_env,input_section, dft_control)
       CALL get_qs_env(qs_env=qs_env,&
                       pw_env=pw_env,&
                       para_env=para_env,&
                       input=input_section,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       subsys=subsys,&
                       dft_control=dft_control,&
                       error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)


       CALL pw_env_get(pw_env=pw_env, pw_pools=pw_pools, error=error)
       print_section => section_vals_get_subs_vals(input_section,"QMMM%PRINT",error=error)
       iw = cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",&
            extension=".qmmmLog",error=error)
       IF (iw>0) &
            WRITE(iw,'(T2,"QMMM|",1X,A)')"Information on the QM/MM Electrostatic Potential:"
       !
       ! Initializing vectors:
       !        Zeroing v_qmmm_rspace
       CALL pw_zero(ks_qmmm_env_loc%v_qmmm_rspace%pw,error=error)
       IF (dft_control%qs_control%semi_empirical) THEN
          ! SEMIEMPIRICAL
          SELECT CASE(qmmm_env%qmmm_coupl_type)
          CASE(do_qmmm_coulomb,do_qmmm_none)
             CALL build_se_qmmm_matrix(qs_env,qmmm_env,mm_particles,mm_cell,para_env,error)
             IF( qmmm_env%qmmm_coupl_type==do_qmmm_none) THEN
                IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')&
                     "No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
             END IF
          CASE (do_qmmm_gauss,do_qmmm_swave)
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="GAUSS or SWAVE QM/MM electrostatic coupling not yet implemented for SE.", &
                  error=error, error_level=cp_failure_level)
          CASE DEFAULT
             IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')"Unknown Coupling..."
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       ELSEIF (dft_control%qs_control%dftb) THEN
          ! DFTB
          CALL cp_unimplemented_error(fromWhere=routineP, &
          message="QM/MM electrostatic coupling not yet implemented for DFTB.", &
          error=error, error_level=cp_failure_level)
       ELSE
          ! QS
          SELECT CASE(qmmm_env%qmmm_coupl_type)
          CASE(do_qmmm_coulomb)
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Coulomb QM/MM electrostatic coupling not implemented for GPW/GAPW.", &
                  error=error, error_level=cp_failure_level)
          CASE(do_qmmm_gauss,do_qmmm_swave)
             IF (iw>0) &
                  WRITE(iw,'(T2,"QMMM|",1X,A)')&
                  "QM/MM Coupling computed collocating the Gaussian Potential Functions."
             interp_section => section_vals_get_subs_vals(input_section,&
                  "QMMM%INTERPOLATOR",error=error)
             CALL    qmmm_elec_with_gaussian(qmmm_env=qmmm_env,&
                                             v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                             mm_particles=mm_particles,&
                                             aug_pools=qmmm_env%aug_pools,&
                                             para_env=para_env,&
                                             eps_mm_rspace=qmmm_env%eps_mm_rspace,&
                                             cube_info=ks_qmmm_env_loc%cube_info,&
                                             pw_pools=pw_pools,&
                                             auxbas_grid=qmmm_env%gridlevel_info%auxbas_grid,&
                                             coarser_grid=qmmm_env%gridlevel_info%coarser_grid,&
                                             interp_section=interp_section,&
                                             mm_cell=mm_cell,&
                                             error=error)
          CASE(do_qmmm_none)
             IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')&
                  "No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
          CASE DEFAULT
             IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')"Unknown Coupling..."
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
          ! Dump info on the electrostatic potential if requested
          IF (BTEST(cp_print_key_should_output(logger%iter_info,print_section,&
               "POTENTIAL",error=error),cp_p_file)) THEN
             iw2 = cp_print_key_unit_nr(logger,print_section,"POTENTIAL",&
                  extension=".qmmmLog",error=error)
             CALL cp_pw_to_cube(ks_qmmm_env_loc%v_qmmm_rspace%pw,iw2,&
                  particles=particles,&
                  stride=section_get_ivals(print_section,"POTENTIAL%STRIDE",error),&
                  title="QM/MM: MM ELECTROSTATIC POTENTIAL ", error=error)
             CALL cp_print_key_finished_output(iw2,logger,print_section,&
                  "POTENTIAL", error=error)
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,print_section,&
            "PROGRAM_RUN_INFO", error=error)
       CALL timestop(handle)
    END IF
  END SUBROUTINE qmmm_el_coupling

! *****************************************************************************
!> \brief Compute the QM/MM electrostatic Interaction collocating the gaussian
!>      Electrostatic Potential
!> \param qmmm_env ...
!> \param v_qmmm ...
!> \param mm_particles ...
!> \param aug_pools ...
!> \param cube_info ...
!> \param para_env ...
!> \param eps_mm_rspace ...
!> \param pw_pools ...
!> \param auxbas_grid ...
!> \param coarser_grid ...
!> \param interp_section ...
!> \param mm_cell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian(qmmm_env, v_qmmm, mm_particles,&
       aug_pools, cube_info, para_env, eps_mm_rspace, pw_pools,&
       auxbas_grid, coarser_grid, interp_section, mm_cell, error)
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(pw_p_type), INTENT(INOUT)           :: v_qmmm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    INTEGER, INTENT(IN)                      :: auxbas_grid, coarser_grid
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, igrid, &
                                                ilevel, kind_interp, lb(3), &
                                                ngrids, ub(3)
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_chrg),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_index),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(aug_pools),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pools),cp_failure_level,routineP,error,failure)
    !Statements
    CALL timeset(routineN,handle)
    ngrids=SIZE(pw_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw,error=error)
    END DO

    CALL qmmm_elec_with_gaussian_low( grids, mm_particles,&
         qmmm_env%mm_atom_chrg, qmmm_env%mm_el_pot_radius, qmmm_env%mm_atom_index,  &
         qmmm_env%num_mm_atoms, cube_info, para_env, eps_mm_rspace, qmmm_env%pgfs,  &
         auxbas_grid, coarser_grid, aug_pools, qmmm_env%potentials,       &
         mm_cell=mm_cell, dOmmOqm=qmmm_env%dOmmOqm, periodic=qmmm_env%periodic,     &
         per_potentials=qmmm_env%per_potentials, par_scheme=qmmm_env%par_scheme,    &
         qmmm_spherical_cutoff=qmmm_env%spherical_cutoff, error=error)

    IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
       CALL qmmm_elec_with_gaussian_low( grids, qmmm_env%added_charges%added_particles,   &
         qmmm_env%added_charges%mm_atom_chrg, qmmm_env%added_charges%mm_el_pot_radius,    &
         qmmm_env%added_charges%mm_atom_index, qmmm_env%added_charges%num_mm_atoms,       &
         cube_info, para_env, eps_mm_rspace, qmmm_env%added_charges%pgfs,auxbas_grid,     &
         coarser_grid, aug_pools, qmmm_env%added_charges%potentials,            &
         mm_cell=mm_cell, dOmmOqm=qmmm_env%dOmmOqm, periodic=qmmm_env%periodic,           &
         per_potentials=qmmm_env%per_potentials, par_scheme=qmmm_env%par_scheme,          &
         qmmm_spherical_cutoff=qmmm_env%spherical_cutoff, error=error)
    END IF
    ! Sumup all contributions according the parallelization scheme
    IF (qmmm_env%par_scheme==do_par_atom) THEN
       DO ilevel=1,SIZE(grids)
          CALL mp_sum(grids(ilevel)%pw%cr3d,para_env%group)
       END DO
    END IF
    ! RealSpace Interpolation
    CALL section_vals_val_get(interp_section,"kind", i_val=kind_interp, error=error)
    SELECT CASE(kind_interp)
    CASE(spline3_nopbc_interp, spline3_pbc_interp)
       ! Spline Iterpolator
       CALL mp_sync(para_env%group)
       CALL timeset(TRIM(routineN)//":spline3Int",handle2)
       DO Ilevel = coarser_grid, auxbas_grid+1, -1
          CALL pw_prolongate_s3(grids(Ilevel  )%pw,&
                                grids(Ilevel-1)%pw,&
                                aug_pools(Ilevel)%pool,&
                                param_section=interp_section,&
                                error=error)
       END DO
       CALL timestop(handle2)
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    lb = v_qmmm%pw%pw_grid%bounds_local(1,:)
    ub = v_qmmm%pw%pw_grid%bounds_local(2,:)

    v_qmmm%pw%cr3d = grids(auxbas_grid)%pw%cr3d (lb(1):ub(1),&
                                                 lb(2):ub(2),&
                                                 lb(3):ub(3) )

    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)

    CALL timestop(handle)
  END SUBROUTINE qmmm_elec_with_gaussian

! *****************************************************************************
!> \brief Compute the QM/MM electrostatic Interaction collocating the gaussian
!>      Electrostatic Potential - Low Level
!> \param tmp_grid ...
!> \param mm_particles ...
!> \param mm_charges ...
!> \param mm_el_pot_radius ...
!> \param mm_atom_index ...
!> \param num_mm_atoms ...
!> \param cube_info ...
!> \param para_env ...
!> \param eps_mm_rspace ...
!> \param pgfs ...
!> \param auxbas_grid ...
!> \param coarser_grid ...
!> \param aug_pools ...
!> \param potentials ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param periodic ...
!> \param per_potentials ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error variable to control error logging ,stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_low( tmp_grid, mm_particles, mm_charges,&
       mm_el_pot_radius, mm_atom_index, num_mm_atoms, cube_info, para_env,   &
       eps_mm_rspace, pgfs, auxbas_grid, coarser_grid, aug_pools,  &
       potentials, mm_cell, dOmmOqm, periodic, per_potentials, par_scheme,   &
       qmmm_spherical_cutoff, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: tmp_grid
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    INTEGER, INTENT(IN)                      :: auxbas_grid, coarser_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    LOGICAL, INTENT(IN)                      :: periodic
    TYPE(qmmm_per_pot_p_type), &
      DIMENSION(:), POINTER                  :: per_potentials
    INTEGER, INTENT(IN)                      :: par_scheme
    REAL(KIND=dp), INTENT(IN)                :: qmmm_spherical_cutoff(2)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_low', &
      routineNb = 'qmmm_elec_gaussian_low', routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, IGauss, &
                                                ilevel, Imm, IndMM, IRadTyp, &
                                                LIndMM, myind, n_rep_real(3), &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: bo2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, height, &
                                                sph_chrg_factor, W
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    TYPE(qmmm_gaussian_type), POINTER        :: pgf
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot
    TYPE(qmmm_pot_type), POINTER             :: pot

    NULLIFY(pgf, pot, per_pot, xdat, ydat, zdat)
    CALL timeset(routineN,handle)
    CALL timeset(routineNb//"_G",handle2)
    bo2 = tmp_grid(auxbas_grid)%pw%pw_grid%bounds
    ALLOCATE (xdat(bo2(1,1):bo2(2,1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(bo2(1,2):bo2(2,2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(bo2(1,3):bo2(2,3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (par_scheme==do_par_atom) myind = 0
    Radius: DO IRadTyp = 1, SIZE(pgfs)
       pgf => pgfs(IRadTyp)%pgf
       pot => potentials(IRadTyp)%pot
       n_rep_real = 0
       IF (periodic) THEN
          per_pot => per_potentials(IRadTyp)%pot
          n_rep_real = per_pot%n_rep_real
       END IF
       Gaussian: DO IGauss = 1, pgf%Number_of_Gaussians
          alpha     = 1.0_dp / pgf%Gk(IGauss)
          alpha     = alpha * alpha
          height    = pgf%Ak(IGauss)
          ilevel    = pgf%grid_level(IGauss)
          Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
             IF (par_scheme==do_par_atom) THEN
                myind = myind + 1
                IF (MOD(myind,para_env%num_pe)/=para_env%mepos) CYCLE Atoms
             END IF
             LIndMM    =   pot%mm_atom_index(Imm)
             IndMM     =   mm_atom_index(LIndMM)
             ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm, mm_cell)+dOmmOqm
             W         =   mm_charges(LIndMM) * height
             ! Possible Spherical Cutoff
             IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
                CALL spherical_cutoff_factor(qmmm_spherical_cutoff, ra, sph_chrg_factor, error)
                W = W * sph_chrg_factor
             END IF
             IF (ABS(W)<= EPSILON(0.0_dp)) CYCLE Atoms
             CALL collocate_gf_rspace_NoPBC(zetp=alpha,&
                                            rp=ra,&
                                            scale=-1.0_dp,&
                                            W=W,&
                                            pwgrid=tmp_grid(ilevel)%pw,&
                                            cube_info=cube_info(ilevel),&
                                            eps_mm_rspace=eps_mm_rspace,&
                                            xdat=xdat,&
                                            ydat=ydat,&
                                            zdat=zdat,&
                                            bo2=bo2,&
                                            n_rep_real=n_rep_real,&
                                            mm_cell=mm_cell)
          END DO Atoms
       END DO Gaussian
    END DO Radius
    IF (ASSOCIATED(xdat)) THEN
       DEALLOCATE (xdat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(ydat)) THEN
       DEALLOCATE (ydat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(zdat)) THEN
       DEALLOCATE (zdat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    CALL timestop(handle2)
    CALL timeset(routineNb//"_R",handle2)
    IF (periodic) THEN
       ! Long Range Part of the QM/MM Potential with Gaussians With Periodic Boundary Conditions
       CALL qmmm_elec_with_gaussian_LG    (pgfs=pgfs,&
                                           cgrid=tmp_grid(coarser_grid)%pw,&
                                           mm_charges=mm_charges,&
                                           mm_atom_index=mm_atom_index,&
                                           mm_particles=mm_particles,&
                                           para_env=para_env,&
                                           coarser_grid_level=coarser_grid,&
                                           per_potentials=per_potentials,&
                                           mm_cell=mm_cell,&
                                           dOmmOqm=dOmmOqm,&
                                           par_scheme=par_scheme,&
                                           qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                           error=error)
    ELSE
       ! Long Range Part of the QM/MM Potential with Gaussians
       CALL qmmm_elec_with_gaussian_LR    (pgfs=pgfs,&
                                           grid=tmp_grid(coarser_grid)%pw,&
                                           mm_charges=mm_charges,&
                                           mm_atom_index=mm_atom_index,&
                                           mm_particles=mm_particles,&
                                           mm_el_pot_radius=mm_el_pot_radius,&
                                           para_env=para_env,&
                                           coarser_grid_level=coarser_grid,&
                                           potentials=potentials,&
                                           mm_cell=mm_cell,&
                                           dOmmOqm=dOmmOqm,&
                                           par_scheme=par_scheme,&
                                           qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                           error=error)
    END IF
    CALL timestop(handle2)
    CALL timestop(handle)

  END SUBROUTINE qmmm_elec_with_gaussian_low

! *****************************************************************************
!> \brief Compute the QM/MM electrostatic Interaction collocating
!>      (1/R - Sum_NG Gaussians) on the coarser grid level in G-SPACE
!>      Long Range QM/MM Electrostatic Potential with Gaussian - Low Level
!>      PERIODIC BOUNDARY CONDITION VERSION
!> \param pgfs ...
!> \param cgrid ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_particles ...
!> \param para_env ...
!> \param coarser_grid_level ...
!> \param per_potentials ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2004 created [tlaino]
!> \author Teodoro Laino
!> \note
!>      This version includes the explicit code of Eval_Interp_Spl3_pbc
!>      in order to achieve better performance
! *****************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_LG(pgfs, cgrid, mm_charges, mm_atom_index,&
       mm_particles, para_env, coarser_grid_level, per_potentials,&
       mm_cell, dOmmOqm, par_scheme, qmmm_spherical_cutoff, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: cgrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: coarser_grid_level
    TYPE(qmmm_per_pot_p_type), &
      DIMENSION(:), POINTER                  :: per_potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, INTENT(IN)                      :: par_scheme
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: qmmm_spherical_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_LG', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ii1, ii2, ii3, ii4, ij1, ij2, ij3, ij4, ik1, ik2, &
      ik3, ik4, Imm, IndMM, IRadTyp, ivec(3), j, k, LIndMM, my_j, my_k, &
      myind, npts(3)
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: a1, a2, a3, abc_X(4,4), abc_X_Y(4), b1, b2, b3, c1, c2, &
      c3, d1, d2, d3, dr1, dr1c, dr2, dr2c, dr3, dr3c, e1, e2, e3, f1, f2, &
      f3, g1, g2, g3, h1, h2, h3, p1, p2, p3, q1, q2, q3, qt, r1, r2, r3, &
      rt1, rt2, rt3, rv1, rv2, rv3, s1, s2, s3, s4, sph_chrg_factor, t1, t2, &
      t3, t4, u1, u2, u3, v1, v2, v3, v4, val, xd1, xd2, xd3, xs1, xs2, xs3
    REAL(KIND=dp), DIMENSION(3)              :: ra, vec
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid, grid2
    TYPE(pw_type), POINTER                   :: pw
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(grid, pw)
    dr1c = cgrid%pw_grid%dr(1)
    dr2c = cgrid%pw_grid%dr(2)
    dr3c = cgrid%pw_grid%dr(3)
    gbo  = cgrid%pw_grid%bounds
    bo   = cgrid%pw_grid%bounds_local
    grid2=>cgrid%cr3d
    IF (par_scheme==do_par_atom) myind = 0
    Radius: DO IRadTyp = 1, SIZE(pgfs)
       per_pot => per_potentials(IRadTyp)%pot
       pw => per_pot%TabLR
       npts = pw%pw_grid%npts
       dr1  = pw%pw_grid%dr(1)
       dr2  = pw%pw_grid%dr(2)
       dr3  = pw%pw_grid%dr(3)
       grid => pw%cr3d(:,:,:)
       Atoms: DO Imm = 1, SIZE(per_pot%mm_atom_index)
          IF (par_scheme==do_par_atom) THEN
             myind = myind + 1
             IF (MOD(myind,para_env%num_pe)/=para_env%mepos) CYCLE
          END IF
          LIndMM    =   per_pot%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
          qt        =   mm_charges(LIndMM)
          ! Possible Spherical Cutoff
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, ra, sph_chrg_factor, error)
             qt = qt * sph_chrg_factor
          END IF
          IF (ABS(qt)<= EPSILON(0.0_dp)) CYCLE Atoms
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k=k-gbo(1,3)
             xs3  = REAL(my_k,dp)*dr3c
             my_j=bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2c
             rv3 = rt3 - xs3
             vec(3) = rv3
             ivec(3) = FLOOR(vec(3)/pw%pw_grid%dr(3))
                   xd3  = (vec(3)/dr3)-REAL(ivec(3),kind=dp)
                   ik1 = MODULO(ivec(3)-1,npts(3))+1
                   ik2 = MODULO(ivec(3)  ,npts(3))+1
                   ik3 = MODULO(ivec(3)+1,npts(3))+1
                   ik4 = MODULO(ivec(3)+2,npts(3))+1
             p1  = 3.0_dp + xd3
             p2  = p1*p1
             p3  = p2*p1
             q1  = 2.0_dp + xd3
             q2  = q1*q1
             q3  = q2*q1
             r1  = 1.0_dp + xd3
             r2  = r1*r1
             r3  = r2*r1
             u1  = xd3
             u2  = u1*u1
             u3  = u2*u1
             v1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*p1 + 12.0_dp*p2 - p3)
             v2 = -22.0_dp/3.0_dp + 10.0_dp*q1 - 4.0_dp*q2 + 0.5_dp*q3
             v3 =   2.0_dp/3.0_dp -  2.0_dp*r1 + 2.0_dp*r2 - 0.5_dp*r3
             v4 =   1.0_dp/6.0_dp*u3
             DO j =  bo(1,2), bo(2,2)
                xs1 = (bo(1,1)-gbo(1,1))*dr1c
                rv2 = rt2 - xs2
                vec(2) = rv2
                ivec(2) = FLOOR(vec(2)/pw%pw_grid%dr(2))
                xd2  = (vec(2)/dr2)-REAL(ivec(2),kind=dp)
                   ij1 = MODULO(ivec(2)-1,npts(2))+1
                   ij2 = MODULO(ivec(2)  ,npts(2))+1
                   ij3 = MODULO(ivec(2)+1,npts(2))+1
                   ij4 = MODULO(ivec(2)+2,npts(2))+1
                e1  = 3.0_dp + xd2
                e2  = e1*e1
                e3  = e2*e1
                f1  = 2.0_dp + xd2
                f2  = f1*f1
                f3  = f2*f1
                g1  = 1.0_dp + xd2
                g2  = g1*g1
                g3  = g2*g1
                h1  = xd2
                h2  = h1*h1
                h3  = h2*h1
                s1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*e1 + 12.0_dp*e2 - e3)
                s2 = -22.0_dp/3.0_dp + 10.0_dp*f1 - 4.0_dp*f2 + 0.5_dp*f3
                s3 =   2.0_dp/3.0_dp -  2.0_dp*g1 + 2.0_dp*g2 - 0.5_dp*g3
                s4 =   1.0_dp/6.0_dp*h3
                DO i =  bo(1,1), bo(2,1)
                   rv1  = rt1 - xs1
                   vec(1) = rv1
                   ivec(1) = FLOOR(vec(1)/pw%pw_grid%dr(1))
                   xd1  = (vec(1)/dr1)-REAL(ivec(1),kind=dp)
                   ii1 = MODULO(ivec(1)-1,npts(1))+1
                   ii2 = MODULO(ivec(1)  ,npts(1))+1
                   ii3 = MODULO(ivec(1)+1,npts(1))+1
                   ii4 = MODULO(ivec(1)+2,npts(1))+1
                   !
                   ! Spline Interpolation 
                   !

                   a1  = 3.0_dp + xd1
                   a2  = a1*a1
                   a3  = a2*a1
                   b1  = 2.0_dp + xd1
                   b2  = b1*b1
                   b3  = b2*b1
                   c1  = 1.0_dp + xd1
                   c2  = c1*c1
                   c3  = c2*c1
                   d1  = xd1
                   d2  = d1*d1
                   d3  = d2*d1
                   t1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*a1 + 12.0_dp*a2 - a3)
                   t2 = -22.0_dp/3.0_dp + 10.0_dp*b1 - 4.0_dp*b2 + 0.5_dp*b3
                   t3 =   2.0_dp/3.0_dp -  2.0_dp*c1 + 2.0_dp*c2 - 0.5_dp*c3
                   t4 =   1.0_dp/6.0_dp*d3

                   abc_X(1,1) = grid(ii1,ij1,ik1)*v1 + grid(ii1,ij1,ik2)*v2 + grid(ii1,ij1,ik3)*v3 + grid(ii1,ij1,ik4)*v4
                   abc_X(1,2) = grid(ii1,ij2,ik1)*v1 + grid(ii1,ij2,ik2)*v2 + grid(ii1,ij2,ik3)*v3 + grid(ii1,ij2,ik4)*v4
                   abc_X(1,3) = grid(ii1,ij3,ik1)*v1 + grid(ii1,ij3,ik2)*v2 + grid(ii1,ij3,ik3)*v3 + grid(ii1,ij3,ik4)*v4
                   abc_X(1,4) = grid(ii1,ij4,ik1)*v1 + grid(ii1,ij4,ik2)*v2 + grid(ii1,ij4,ik3)*v3 + grid(ii1,ij4,ik4)*v4
                   abc_X(2,1) = grid(ii2,ij1,ik1)*v1 + grid(ii2,ij1,ik2)*v2 + grid(ii2,ij1,ik3)*v3 + grid(ii2,ij1,ik4)*v4
                   abc_X(2,2) = grid(ii2,ij2,ik1)*v1 + grid(ii2,ij2,ik2)*v2 + grid(ii2,ij2,ik3)*v3 + grid(ii2,ij2,ik4)*v4
                   abc_X(2,3) = grid(ii2,ij3,ik1)*v1 + grid(ii2,ij3,ik2)*v2 + grid(ii2,ij3,ik3)*v3 + grid(ii2,ij3,ik4)*v4
                   abc_X(2,4) = grid(ii2,ij4,ik1)*v1 + grid(ii2,ij4,ik2)*v2 + grid(ii2,ij4,ik3)*v3 + grid(ii2,ij4,ik4)*v4
                   abc_X(3,1) = grid(ii3,ij1,ik1)*v1 + grid(ii3,ij1,ik2)*v2 + grid(ii3,ij1,ik3)*v3 + grid(ii3,ij1,ik4)*v4
                   abc_X(3,2) = grid(ii3,ij2,ik1)*v1 + grid(ii3,ij2,ik2)*v2 + grid(ii3,ij2,ik3)*v3 + grid(ii3,ij2,ik4)*v4
                   abc_X(3,3) = grid(ii3,ij3,ik1)*v1 + grid(ii3,ij3,ik2)*v2 + grid(ii3,ij3,ik3)*v3 + grid(ii3,ij3,ik4)*v4
                   abc_X(3,4) = grid(ii3,ij4,ik1)*v1 + grid(ii3,ij4,ik2)*v2 + grid(ii3,ij4,ik3)*v3 + grid(ii3,ij4,ik4)*v4
                   abc_X(4,1) = grid(ii4,ij1,ik1)*v1 + grid(ii4,ij1,ik2)*v2 + grid(ii4,ij1,ik3)*v3 + grid(ii4,ij1,ik4)*v4
                   abc_X(4,2) = grid(ii4,ij2,ik1)*v1 + grid(ii4,ij2,ik2)*v2 + grid(ii4,ij2,ik3)*v3 + grid(ii4,ij2,ik4)*v4
                   abc_X(4,3) = grid(ii4,ij3,ik1)*v1 + grid(ii4,ij3,ik2)*v2 + grid(ii4,ij3,ik3)*v3 + grid(ii4,ij3,ik4)*v4
                   abc_X(4,4) = grid(ii4,ij4,ik1)*v1 + grid(ii4,ij4,ik2)*v2 + grid(ii4,ij4,ik3)*v3 + grid(ii4,ij4,ik4)*v4

                   abc_X_Y(1) = abc_X(1,1)*t1 + abc_X(2,1)*t2 + abc_X(3,1)*t3 + abc_X(4,1)*t4
                   abc_X_Y(2) = abc_X(1,2)*t1 + abc_X(2,2)*t2 + abc_X(3,2)*t3 + abc_X(4,2)*t4
                   abc_X_Y(3) = abc_X(1,3)*t1 + abc_X(2,3)*t2 + abc_X(3,3)*t3 + abc_X(4,3)*t4
                   abc_X_Y(4) = abc_X(1,4)*t1 + abc_X(2,4)*t2 + abc_X(3,4)*t3 + abc_X(4,4)*t4

                   val = abc_X_Y(1)*s1 + abc_X_Y(2)*s2 + abc_X_Y(3)*s3 + abc_X_Y(4)*s4

                   grid2(i,j,k) = grid2(i,j,k) - val * qt
                   xs1 = xs1 + dr1c
                END DO
                xs2 = xs2 + dr2c
             END DO
          END DO LoopOnGrid
       END DO Atoms
    END DO Radius
    CALL timestop(handle)
  END SUBROUTINE qmmm_elec_with_gaussian_LG

! *****************************************************************************
!> \brief Compute the QM/MM electrostatic Interaction collocating
!>      (1/R - Sum_NG Gaussians) on the coarser grid level.
!>      Long Range QM/MM Electrostatic Potential with Gaussian - Low Level
!> \param pgfs ...
!> \param grid ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_particles ...
!> \param mm_el_pot_radius ...
!> \param para_env ...
!> \param coarser_grid_level ...
!> \param potentials ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_LR(pgfs, grid, mm_charges, mm_atom_index,&
       mm_particles, mm_el_pot_radius, para_env,coarser_grid_level, potentials,&
       mm_cell, dOmmOqm, par_scheme, qmmm_spherical_cutoff, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: grid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: coarser_grid_level
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, INTENT(IN)                      :: par_scheme
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: qmmm_spherical_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_LR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, Imm, IndMM, &
                                                IRadTyp, ix, j, k, LIndMM, &
                                                my_j, my_k, myind, n1, n2, n3
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dr1, dr2, dr3, dx, qt, r, r2, &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                rx, rx2, rx3, &
                                                sph_chrg_factor, Term, xs1, &
                                                xs2, xs3
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid2
    TYPE(qmmm_pot_type), POINTER             :: pot

    CALL timeset(routineN,handle)
    failure = .FALSE.
    n1   = grid%pw_grid%npts(1)
    n2   = grid%pw_grid%npts(2)
    n3   = grid%pw_grid%npts(3)
    dr1  = grid%pw_grid%dr(1)
    dr2  = grid%pw_grid%dr(2)
    dr3  = grid%pw_grid%dr(3)
    gbo  = grid%pw_grid%bounds
    bo   = grid%pw_grid%bounds_local
    grid2=>grid%cr3d
    IF (par_scheme==do_par_atom) myind=0
    Radius: DO IRadTyp = 1, SIZE(pgfs)
       pot => potentials(IRadTyp)%pot
       dx     =  Pot%dx
       pot0_2 => Pot%pot0_2
       Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
          IF (par_scheme==do_par_atom) THEN
             myind = myind + 1
             IF (MOD(myind,para_env%num_pe)/=para_env%mepos) CYCLE
          END IF
          LIndMM    =   pot%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
          qt        =   mm_charges(LIndMM)
          ! Possible Spherical Cutoff
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, ra, sph_chrg_factor, error)
             qt = qt * sph_chrg_factor
          END IF
          IF (ABS(qt)<= EPSILON(0.0_dp)) CYCLE Atoms
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k=k-gbo(1,3)
             xs3  = REAL(my_k,dp)*dr3
             my_j=bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                xs1 = (bo(1,1)-gbo(1,1))*dr1
                rv2 = rt2 - xs2
                DO i =  bo(1,1), bo(2,1)
                   rv1  = rt1 - xs1
                   r2   = rv1*rv1 + rv2*rv2 + rv3*rv3
                   r    = SQRT(r2)
                   ix  = FLOOR(r/dx)+1
                   rx  = (r-REAL(ix-1,dp)*dx)/dx
                   rx2 = rx*rx
                   rx3 = rx2*rx
                   Term = pot0_2(1,ix  )*(1.0_dp-3.0_dp*rx2+2.0_dp*rx3)  &
                         +pot0_2(2,ix  )*(rx-2.0_dp*rx2+rx3)             &
                         +pot0_2(1,ix+1)*(3.0_dp*rx2-2.0_dp*rx3)         &
                         +pot0_2(2,ix+1)*(-rx2+rx3)
                   grid2(i,j,k) = grid2(i,j,k) - Term * qt
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
       END DO Atoms
    END DO Radius
    CALL timestop(handle)
  END SUBROUTINE qmmm_elec_with_gaussian_LR

END MODULE qmmm_gpw_energy
