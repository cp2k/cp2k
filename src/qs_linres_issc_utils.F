!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Chemical shift calculation by dfpt
!>      Initialization of the issc_env, creation of the special neighbor lists
!>      Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!>      Write output
!>      Deallocate everything
!> \note
!>      The psi0 should be localized
!>      the Sebastiani method works within the assumption that the orbitals are
!>      completely contained in the simulation box
! *****************************************************************************
MODULE qs_linres_issc_utils

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_methods,               ONLY: linres_solver
  USE qs_linres_types,                 ONLY: get_issc_env,&
                                             issc_env_create,&
                                             issc_env_type,&
                                             linres_control_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: set_up_op_sm
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: issc_env_cleanup, issc_env_init, issc_response

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_issc_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the issc environment
! *****************************************************************************
  SUBROUTINE issc_response(issc_env,p_env,qs_env,iatom,error)
    !
    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_response', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ijdir, ispin, &
                                                istat, jdir, nao, nmo, &
                                                nspins, output_unit
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: chk, fro
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fc_psi0, h1_psi0, psi0_order, &
                                                psi1, psi1_fc
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: efg_psi0, psi1_efg, psi1_pso, &
                                                pso_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: issc_section, lr_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    !
    NULLIFY(dft_control, linres_control, lr_section, issc_section)
    NULLIFY(logger, mpools, psi1,h1_psi0, mo_coeff,para_env)
    NULLIFY(tmp_fm_struct, psi1_fc, psi1_efg, psi1_pso, pso_psi0, fc_psi0, efg_psi0, psi0_order)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    issc_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%SPINSPIN",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)
    IF(output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "*** Self consistent optimization of the response wavefunctions ***"
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    mpools=mpools,&
                    linres_control=linres_control,&
                    mos=mos,&
                    para_env=para_env,&
                    error=error)

    nspins = dft_control%nspins

    CALL get_issc_env(issc_env=issc_env, &
                      !nao=nao, &
                      !centers_set=centers_set, &
                      !nbr_center=nbr_center, &
                      !center_list=center_list, &
                      !list_cubes=list_cubes, &
                      psi1_efg=psi1_efg, &
                      psi1_pso=psi1_pso, &
                      psi1_fc=psi1_fc, &
                      efg_psi0=efg_psi0, &
                      pso_psi0=pso_psi0, &
                      fc_psi0=fc_psi0, &
                      !psi0_order=psi0_order,&
                      error=error)
    !
    ! allocate the vectors
    ALLOCATE(psi0_order(nspins))
    ALLOCATE(psi1(nspins),h1_psi0(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       psi0_order(ispin)%matrix => mo_coeff
       CALL cp_fm_get_info(mo_coeff,ncol_global=nmo,nrow_global=nao,error=error)
       NULLIFY(tmp_fm_struct,psi1(ispin)%matrix,h1_psi0(ispin)%matrix)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
            &                   ncol_global=nmo,&
            &                   context=mo_coeff%matrix_struct%context,&
            &                   error=error)
       CALL cp_fm_create(psi1(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(h1_psi0(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    ENDDO
    chk = 0.0_dp
    !
    ! operator efg
    ijdir = 0
    DO idir = 1,3
    DO jdir = idir,3
       ijdir = ijdir+1
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1_efg(ispin,ijdir)%matrix,0.0_dp,error=error)
       ENDDO
       !IF(.NOT.current_env%simple_done(idir)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator efg_"//ACHAR(idir+119)//ACHAR(jdir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(p_psi0(ispin,ijdir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(efg_psi0(ispin,ijdir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          !
          linres_control%converged = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,error=error)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_efg(ispin,ijdir)%matrix,error=error)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
             chk = chk + fro
          ENDDO
          !
          ! print response functions
          !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
          !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          !   ncubes = SIZE(list_cubes,1)
          !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
          !   DO ispin = 1,nspins
          !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
          !            centers_set(ispin)%array,print_key,'psi1_efg',&
          !            idir=ijdir,ispin=ispin,error=error)
          !   ENDDO ! ispin
          !ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          ! Write the result in the restart file
          !current_env%simple_done(ijdir) = .TRUE.
          !current_env%simple_converged(ijdir) = linres_control%converged
          !CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir,error=error)
       !ENDIF

    ENDDO ! jdir
    ENDDO ! idir
    !
    ! operator pso
    DO idir = 1,3
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1_pso(ispin,idir)%matrix,0.0_dp,error=error)
       ENDDO
       !IF(.NOT.current_env%simple_done(idir+6)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator pso_"//ACHAR(idir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(pso_psi0(ispin,idir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          ! 
          linres_control%converged = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,error=error)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_pso(ispin,idir)%matrix,error=error)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
             chk = chk + fro
          ENDDO
          !
          ! print response functions
          !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
          !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          !   ncubes = SIZE(list_cubes,1)
          !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
          !   DO ispin = 1,nspins
          !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
          !           centers_set(ispin)%array,print_key,'psi1_pso',&
          !           idir=idir,ispin=ispin,error=error)
          !   ENDDO ! ispin
          !ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          ! Write the result inthe restart file
          !current_env%simple_done(idir+6) = .TRUE.
          !current_env%simple_converged(idir+6) = linres_control%converged
          !CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir+3,error=error)
       !ENDIF
    ENDDO ! idir
    !
    ! operator fc
    DO ispin = 1,nspins
       CALL cp_fm_set_all(psi1_fc(ispin)%matrix,0.0_dp,error=error)
    ENDDO
    !IF(.NOT.current_env%simple_done(idir+9)) THEN
       IF(output_unit>0) THEN
          WRITE(output_unit,"(T10,A)") "Response to the perturbation operator fc"
       ENDIF
       !
       !Initial guess for psi1
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
          !CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
          !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
       ENDDO
       !
       !DO scf cycle to optimize psi1
       DO ispin = 1,nspins
          CALL cp_fm_to_fm(fc_psi0(ispin)%matrix,h1_psi0(ispin)%matrix,error=error)
       ENDDO
       !
       ! 
       linres_control%converged = .FALSE.
       CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,error=error)
       !
       !
       ! copy the response
       DO ispin=1,nspins
          CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_fc(ispin)%matrix,error=error)
          CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
          chk = chk + fro
       ENDDO
       !
       ! print response functions
       !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
       !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
       !   ncubes = SIZE(list_cubes,1)
       !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
       !   DO ispin = 1,nspins
       !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
       !           centers_set(ispin)%array,print_key,'psi1_pso',&
       !           idir=idir,ispin=ispin,error=error)
       !   ENDDO ! ispin
       !ENDIF ! print response functions
       !
       !
       IF(output_unit>0) THEN
          WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
       ENDIF
       ! Write the result inthe restart file
       !current_env%simple_done(idir+9) = .TRUE.
       !current_env%simple_converged(idir+9) = linres_control%converged
       !CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir+3,error=error)
    !ENDIF
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| response: CheckSum =',chk
    ENDIF
    !
    ! clean up
    DO ispin = 1,nspins
       CALL cp_fm_release(psi1(ispin)%matrix,error=error)
       CALL cp_fm_release(h1_psi0(ispin)%matrix,error=error)
    ENDDO
    DEALLOCATE(psi1,h1_psi0,psi0_order,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE issc_response

! *****************************************************************************
!> \brief Initialize the issc environment
! *****************************************************************************
  SUBROUTINE issc_env_init(issc_env,qs_env,error)
    !
    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, idir, ini, ir, &
                                                ispin, istat, m, n, n_rep, &
                                                nao, natom, nspins, &
                                                output_unit
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, gapw
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: issc_section, lr_section

!

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(linres_control)
    NULLIFY(logger,issc_section)
    NULLIFY(tmp_fm_struct)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         &                             extension=".linresLog",error=error)
    
    IF(issc_env%ref_count /= 0) THEN
       CALL issc_env_cleanup(issc_env,error=error)
    ENDIF

    IF(output_unit>0) THEN
       WRITE(output_unit,"(/,T20,A,/)") "*** Start indirect spin-spin coupling Calculation ***"
       WRITE(output_unit,"(T10,A,/)") "Inizialization of the ISSC environment"
    ENDIF

    CALL issc_env_create(issc_env,error=error)
    !
    issc_section => section_vals_get_subs_vals(qs_env%input, &
         &          "PROPERTIES%LINRES%SPINSPIN",error=error)
    !CALL section_vals_val_get(nmr_section,"INTERPOLATE_SHIFT",l_val=nmr_env%interpolate_shift,error=error)
    !CALL section_vals_val_get(nmr_section,"SHIFT_GAPW_RADIUS",r_val=nmr_env%shift_gapw_radius,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         &          dft_control=dft_control,&
         &          linres_control=linres_control,&
         &          particle_set=particle_set,&
         &          mos=mos,&
         &          error=error)
    !
    !
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure.AND.(output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get indirect spin-spin coupling parameters within PBC you need to localize zero order orbitals '
    ENDIF
    !
    ! read the list of atoms on which the issc need to be calculated
    CALL section_vals_val_get(issc_section,"ISSC_ON_ATOM_LIST",n_rep_val=n_rep,error=error)
    !
    !
    NULLIFY(issc_env%issc_on_atom_list)
    n = 0
    DO ir = 1,n_rep
       NULLIFY(list)
       CALL section_vals_val_get(issc_section,"ISSC_ON_ATOM_LIST",i_rep_val=ir,i_vals=list,error=error)
       IF(ASSOCIATED(list)) THEN
          CALL reallocate(issc_env%issc_on_atom_list,1,n+SIZE(list))
          DO ini = 1,SIZE(list)
             issc_env%issc_on_atom_list(ini+n) = list(ini)
          ENDDO
          n = n + SIZE(list)
       ENDIF
    ENDDO
    !
    IF(.NOT.ASSOCIATED(issc_env%issc_on_atom_list)) THEN
       ALLOCATE(issc_env%issc_on_atom_list(natom),STAT=istat)
       CPPrecondition(istat.EQ.0,cp_warning_level,routineP,error,failure)
       DO iatom = 1,natom
          issc_env%issc_on_atom_list(iatom) = iatom
       ENDDO
    ENDIF
    issc_env%issc_natms = SIZE(issc_env%issc_on_atom_list)
    !
    !
    ! Initialize the issc tensor
    ALLOCATE(issc_env%issc(3,3,issc_env%issc_natms,issc_env%issc_natms),&
         &   issc_env%issc_loc(3,3,issc_env%issc_natms,issc_env%issc_natms),&
         &   STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    issc_env%issc(:,:,:,:) = 0.0_dp
    issc_env%issc_loc(:,:,:,:) = 0.0_dp
    !
    ! allocation
    ALLOCATE(issc_env%efg_psi0(nspins,6),issc_env%pso_psi0(nspins,3),issc_env%fc_psi0(nspins),&
         &   issc_env%psi1_efg(nspins,6),issc_env%psi1_pso(nspins,3),issc_env%psi1_fc(nspins),&
         &   STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       !mo_coeff => current_env%psi0_order(ispin)%matrix
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,ncol_global=m,nrow_global=nao,error=error)

       NULLIFY(tmp_fm_struct)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
            &                   ncol_global=m,&
            &                   context=mo_coeff%matrix_struct%context,&
            &                   error=error)
       DO idir = 1,6
          NULLIFY(issc_env%psi1_efg(ispin,idir)%matrix,issc_env%efg_psi0(ispin,idir)%matrix)
          CALL cp_fm_create(issc_env%psi1_efg(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(issc_env%efg_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
       ENDDO
       DO idir = 1,3
          NULLIFY(issc_env%psi1_pso(ispin,idir)%matrix,issc_env%pso_psi0(ispin,idir)%matrix)
          CALL cp_fm_create(issc_env%psi1_pso(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(issc_env%pso_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
       ENDDO
       NULLIFY(issc_env%psi1_fc(ispin)%matrix,issc_env%fc_psi0(ispin)%matrix)
       CALL cp_fm_create(issc_env%psi1_fc(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(issc_env%fc_psi0(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    ENDDO
    !
    ! efg, pso and fc operators
    CALL allocate_matrix_set(issc_env%matrix_efg,6,error=error)
    CALL set_up_op_sm(issc_env%matrix_efg(1)%matrix,qs_env,symmetry="symmetric",&
         &            name="efg (3xx-rr)/3",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_efg(1)%matrix,issc_env%matrix_efg(2)%matrix,&
            'efg xy',target_symmetry="symmetric",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_efg(1)%matrix,issc_env%matrix_efg(3)%matrix,&
            'efg xz',target_symmetry="symmetric",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_efg(1)%matrix,issc_env%matrix_efg(4)%matrix,&
            'efg (3yy-rr)/3',target_symmetry="symmetric",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_efg(1)%matrix,issc_env%matrix_efg(5)%matrix,&
            'efg yz',target_symmetry="symmetric",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_efg(1)%matrix,issc_env%matrix_efg(6)%matrix,&
            'efg (3zz-rr)/3',target_symmetry="symmetric",error=error)
    CALL set_matrix(issc_env%matrix_efg(1)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_efg(2)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_efg(3)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_efg(4)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_efg(5)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_efg(6)%matrix,0.0_dp)
    !
    CALL allocate_matrix_set(issc_env%matrix_pso,3,error=error)
    CALL set_up_op_sm(issc_env%matrix_pso(1)%matrix,qs_env,symmetry="antisymmetric",&
         &            name="pso x",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_pso(1)%matrix,issc_env%matrix_pso(2)%matrix,&
         'pso y',target_symmetry="antisymmetric",error=error)
    CALL replicate_matrix_structure(issc_env%matrix_pso(1)%matrix,issc_env%matrix_pso(3)%matrix,&
         'pso z',target_symmetry="antisymmetric",error=error)
    CALL set_matrix(issc_env%matrix_pso(1)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_pso(2)%matrix,0.0_dp)
    CALL set_matrix(issc_env%matrix_pso(3)%matrix,0.0_dp)
    !
    CALL allocate_matrix_set(issc_env%matrix_fc,1,error=error)
    CALL set_up_op_sm(issc_env%matrix_fc(1)%matrix,qs_env,symmetry="symmetric",&
         &            name="fc",error=error)
    CALL set_matrix(issc_env%matrix_fc(1)%matrix,0.0_dp)
    !
    ! Conversion factors
    IF(output_unit>0) THEN
       WRITE(output_unit,"(T2,A,T60,I4,A)")&
            & " ISSC| spin-spin coupling computed for ",issc_env%issc_natms,' atoms'
    ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE issc_env_init

! *****************************************************************************
!> \brief Deallocate the issc environment
!> \par History
! *****************************************************************************
  SUBROUTINE issc_env_cleanup(issc_env,error)

    TYPE(issc_env_type)                      :: issc_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir, ispin, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    issc_env%ref_count = issc_env%ref_count - 1
    IF(issc_env%ref_count == 0 ) THEN
       IF(ASSOCIATED(issc_env%issc_on_atom_list)) THEN
          DEALLOCATE(issc_env%issc_on_atom_list,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(issc_env%issc)) THEN
          DEALLOCATE(issc_env%issc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(issc_env%issc_loc)) THEN
          DEALLOCATE(issc_env%issc_loc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !efg_psi0
       IF(ASSOCIATED(issc_env%efg_psi0)) THEN
          DO idir = 1,SIZE(issc_env%efg_psi0,2)
             DO ispin = 1,SIZE(issc_env%efg_psi0,1)
                CALL cp_fm_release(issc_env%efg_psi0(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%efg_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !pso_psi0
       IF(ASSOCIATED(issc_env%pso_psi0)) THEN
          DO idir = 1,SIZE(issc_env%pso_psi0,2)
             DO ispin = 1,SIZE(issc_env%pso_psi0,1)
                CALL cp_fm_release(issc_env%pso_psi0(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%pso_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !fc_psi0
       IF(ASSOCIATED(issc_env%fc_psi0)) THEN
          DO ispin = 1,SIZE(issc_env%fc_psi0,1)
             CALL cp_fm_release(issc_env%fc_psi0(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(issc_env%fc_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_efg
       IF(ASSOCIATED(issc_env%psi1_efg)) THEN
          DO idir = 1,SIZE(issc_env%psi1_efg,2)
             DO ispin = 1,SIZE(issc_env%psi1_efg,1)
                CALL cp_fm_release(issc_env%psi1_efg(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%psi1_efg,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_pso
       IF(ASSOCIATED(issc_env%psi1_pso)) THEN
          DO idir = 1,SIZE(issc_env%psi1_pso,2)
             DO ispin = 1,SIZE(issc_env%psi1_pso,1)
                CALL cp_fm_release(issc_env%psi1_pso(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%psi1_pso,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_fc
       IF(ASSOCIATED(issc_env%psi1_fc)) THEN
          DO ispin = 1,SIZE(issc_env%psi1_fc,1)
             CALL cp_fm_release(issc_env%psi1_fc(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(issc_env%psi1_fc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       ! cubes
       !IF(ASSOCIATED(issc_env%list_cubes)) THEN
       !   DEALLOCATE(issc_env%list_cubes,STAT=istat)
       !   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !ENDIF
       !
       !matrix_efg
       IF(ASSOCIATED(issc_env%matrix_efg)) THEN
          CALL deallocate_matrix_set(issc_env%matrix_efg,error=error)
       ENDIF
       !
       !matrix_pso
       IF(ASSOCIATED(issc_env%matrix_pso)) THEN
          CALL deallocate_matrix_set(issc_env%matrix_pso,error=error)
       ENDIF
       !
       !matrix_fc
       IF(ASSOCIATED(issc_env%matrix_fc)) THEN
          CALL deallocate_matrix_set(issc_env%matrix_fc,error=error)
       ENDIF

    ENDIF ! ref count

  END SUBROUTINE issc_env_cleanup

END MODULE qs_linres_issc_utils
