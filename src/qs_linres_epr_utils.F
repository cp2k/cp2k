!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_epr_utils *
!!
!!   NAME
!!     qs_linres_epr_utils
!!
!!   FUNCTION
!!     g tensor calculation by dfpt
!!     Initialization of the epr_env, creation of the special neighbor lists
!!     Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!!     Write output
!!     Deallocate everything
!!
!!   NOTE
!!     The psi0 should be localized
!!     the Sebastiani method works within the assumption that the orbitals are
!!     completely contained in the simulation box
!!
!!   AUTHOR
!!      MI
!!
!!   MODIFICATION HISTORY
!!      created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_epr_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_constants,                 ONLY: restart_guess,&
                                             epr_gauge_r,&
                                             epr_gauge_r_and_Bader,&
                                             epr_gauge_r_and_step_function
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE mathconstants,                   ONLY: fourpi,&
                                             twopi
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: a_fine,&
                                             a_bohr,&
                                             e_charge,&
                                             e_gfactor,&
                                             e_mass,&
                                             mu_perm
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_zero
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_epr_atom_current,      ONLY: list_3c_jrho_atom
  USE qs_linres_epr_current,           ONLY: epr_response_current
  USE qs_linres_epr_op,                ONLY: set_vecp
  USE qs_linres_methods,               ONLY: linres_scf
  USE qs_linres_types,                 ONLY: &
       deallocate_jrho_atom_set, deallocate_nablavks_atom_set, &
       epr_env_create, epr_env_type, get_epr_env, init_jrho_atom_set, &
       init_nablavks_atom_set, jrho_atom_type, linres_control_type, &
       nablavks_atom_type, set_epr_env
  USE qs_loc_methods,                  ONLY: qs_print_cubes
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: set_up_op_sm
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals
  USE qs_rho_atom_types,               ONLY: deallocate_rho_atom_set,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_do_epr, epr_env_cleanup, epr_env_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_utils'

CONTAINS

  SUBROUTINE epr_do_epr(epr_env,p_env,qs_env,error)
    !
    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_do_epr', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, idir, ii, iii, ispin, ist_true, istat, istate, &
      istate2, nao, ncubes, nmo, nstates(2), output_unit
    INTEGER, DIMENSION(:), POINTER           :: list_cubes
    INTEGER, DIMENSION(:, :, :), POINTER     :: statetrueindex
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: dk(3), dkl(3), dl(3)
    REAL(dp), DIMENSION(:), POINTER          :: dkl_vec_ii, dkl_vec_iii
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_dklxp0
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fm_work_ii, fm_work_iii, &
                                                h1_psi0, psi1
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: p_psi0, psi1_D, psi1_p, &
                                                psi1_rxp, rxp_psi0
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: epr_section, lr_section

!

    failure   = .FALSE.
    CALL timeset(routineN,"I"," ",handle)
    
    NULLIFY(ao_mo_fm_pools, cell, dft_control, linres_control, lr_section, epr_section)
    NULLIFY(logger, mpools, psi1,h1_psi0, vecbuf_dklxp0, mos, mo_coeff,para_env )
    NULLIFY(fm_work_ii, fm_work_iii, tmp_fm_struct, dkl_vec_ii,dkl_vec_iii )

    NULLIFY(list_cubes, statetrueindex, centers_set)
    NULLIFY(psi1_p, psi1_rxp, psi1_D, p_psi0, rxp_psi0)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    epr_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%EPR",error=error)
    CALL section_vals_val_get(epr_section,"GAUGE",i_val=epr_env%gauge,error=error)
    SELECT CASE(epr_env%gauge)
    CASE(epr_gauge_r)
       epr_env%gauge_name="R"
    CASE(epr_gauge_r_and_Bader)
       epr_env%gauge_name="R_AND_BADER"
    CASE(epr_gauge_r_and_step_function)
       epr_env%gauge_name="R_AND_STEP_FUNCTION"
    CASE DEFAULT
       CALL stop_program("epr_env_init","unknown gauge, try again...")
    END SELECT

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "*** Self consistent optimization of the response wavefunctions ***"
    END IF

    CALL get_qs_env(qs_env=qs_env,&
          dft_control= dft_control,&
          mpools=mpools,cell=cell,&
          linres_control=linres_control,&
          mos=mos,para_env=para_env,error=error)

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    CALL mpools_get(mpools, ao_mo_fm_pools=ao_mo_fm_pools,error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, psi1, name=routineP//":psi1",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, h1_psi0, name=routineP//":h1_psi0",error=error)

    !Check if restart and from where
    IF(epr_env%restart_epr) THEN
       !CALL epr_restart_epr(epr_env,qs_env,action="read",error=error)
    ENDIF

    CALL get_epr_env(epr_env=epr_env, nao=nao, nstates=nstates, &
                     centers_set=centers_set, list_cubes=list_cubes,&
                     statetrueindex=statetrueindex, &
                     psi1_p=psi1_p, psi1_rxp=psi1_rxp, psi1_D=psi1_D,&
                     p_psi0=p_psi0, rxp_psi0=rxp_psi0,error=error)
    !
    ! operator p
    DO idir = 1,3
       CALL cp_fm_set_all(psi1_p(1,idir)%matrix,0.0_dp,error=error)
    ENDDO
    DO idir = 1,3
       IF(.NOT.epr_env%simpleepr_done(idir)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator P_"//ACHAR(idir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,dft_control%nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(p_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,dft_control%nspins
             CALL cp_fm_to_fm(p_psi0(ispin,idir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO

          linres_control%converged = .FALSE.
          CALL linres_scf(p_env,qs_env,psi1,h1_psi0,error=error)

          IF(output_unit>0) THEN
             WRITE(output_unit,'(T10,A,E24.16)') 'Second order energy P_'//ACHAR(idir+119)//'=',&
                  & p_env%etotal
          ENDIF
          !
          ! the optimized wfns are copied in the  fm
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)")&
                  & "Store the psi1 for the calculation of the response current density "
          ENDIF
          DO ispin=1, dft_control%nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_p(ispin,idir)%matrix,error=error)
          ENDDO
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             DO ispin = 1,dft_control%nspins
                CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
                     &              centers_set(ispin)%array,ionode,'psi1_p',&
                     &              idir=idir,ispin=ispin,stride=section_get_ival(epr_section,&
                     &              "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE",error=error),&
                     &              error=error)
             ENDDO ! ispin
          ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          ! Write the result in the restart file
          epr_env%simpleepr_done(idir) = .TRUE.
          epr_env%simpleepr_converged(idir) = linres_control%converged
          !CALL epr_restart_epr(epr_env,qs_env,action="write",ivec=idir,error=error)
       ENDIF
    ENDDO ! idir
    !
    ! operator rxp
    DO idir = 1,3
       CALL cp_fm_set_all(psi1_rxp(1,idir)%matrix,0.0_dp,error=error)
    ENDDO
    
    DO idir = 1,3
       IF(.NOT.epr_env%simpleepr_done(idir+3)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator L_"//ACHAR(idir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,dft_control%nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,dft_control%nspins
             CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          ! 
          linres_control%converged = .FALSE.
          CALL linres_scf(p_env,qs_env,psi1,h1_psi0,error=error)
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,'(T10,A,E24.16)') 'Second order energy L_'//ACHAR(idir+119)//'=',&
                  & p_env%etotal
          END IF
          !
          !   the optimized wfns are copied in the  fm
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)")&
                  & "Store the psi1 for the calculation of the response current density "
          ENDIF
          DO ispin=1,dft_control%nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_rxp(ispin,idir)%matrix,error=error)
          ENDDO
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             DO ispin = 1,dft_control%nspins
                CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
                     &              centers_set(ispin)%array,ionode,'psi1_rxp',&
                     &              idir=idir,ispin=ispin,stride=section_get_ival(epr_section,&
                     &              "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE",error=error),&
                     &              error=error)
             ENDDO ! ispin
          ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          ! Write the result inthe restart file
          epr_env%simpleepr_done(idir+3) = .TRUE.
          !CALL epr_restart_epr(epr_env,qs_env,action="write",ivec=idir+3,error=error)
       ENDIF
    ENDDO ! idir
    !
    ! operator D
    IF(epr_env%full_epr) THEN
       DO ispin = 1,dft_control%nspins
          CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
          DO idir = 1,3
             CALL cp_fm_set_all(psi1_D(ispin,idir)%matrix,0.0_dp,error=error)
          ENDDO
       ENDDO
       !
       ! The correction is state depedent a loop over the states is necessary
       ALLOCATE(vecbuf_dklxp0(1,nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       vecbuf_dklxp0(1,nao) = 0.0_dp
       !
       ALLOCATE(fm_work_ii(dft_control%nspins),fm_work_iii(dft_control%nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,dft_control%nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo,mo_coeff=mo_coeff)
          NULLIFY(tmp_fm_struct,fm_work_ii(ispin)%matrix,fm_work_iii(ispin)%matrix)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
               &                   ncol_global=nmo,para_env=para_env,&
               &                   context=mo_coeff%matrix_struct%context,&
               &                   error=error)

          CALL cp_fm_create(fm_work_ii(ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_set_all(fm_work_ii(ispin)%matrix,0.0_dp,error=error)
          CALL cp_fm_create(fm_work_iii(ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_set_all(fm_work_iii(ispin)%matrix,0.0_dp,error=error)

          CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       ENDDO ! ispin
       !
       ALLOCATE(dkl_vec_ii(MAXVAL(nstates)),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dkl_vec_iii(MAXVAL(nstates)),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DO idir = 1,3
          !
          ! here we run over the max number of states, then we need 
          ! to be careful with overflow while doing uks calculations.
          DO istate = 1,MAXVAL(nstates) !nstates(ispin)
             !
             IF(output_unit>0) THEN
                WRITE(output_unit,"(T10,A,I4,A)")&
                     & "Response to the perturbation operator (dk-dl)xp for -state- ",&
                     & istate," in dir. "//ACHAR(idir+119)
             ENDIF
             !
             DO ispin = 1,dft_control%nspins
                !
                ! take care that no overflow can occur for uks
                IF(istate.GT.nstates(ispin)) THEN
                   !
                   ! set h1_psi0 and psi1 to zero to avoid problems in linres_scf
                   CALL cp_fm_set_all(h1_psi0(ispin)%matrix,0.0_dp,error=error)
                   CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
                   CYCLE
                ENDIF
                !
                dkl_vec_ii (:)=0.0_dp
                dkl_vec_iii(:)=0.0_dp
                !
                ist_true = statetrueindex(idir,istate,ispin)
                !
                ! the initial guess is the previous set of psi1, just optimized
                CALL set_vecp(idir,ii,iii)
                dk(1:3) = centers_set(ispin)%array(1:3,ist_true)
                !
                DO istate2 = 1,nstates(ispin)
                   dl(1:3) = centers_set(ispin)%array(1:3,istate2)
                   dkl = pbc(dl,dk,cell)
                   dkl_vec_ii(istate2)  = dkl(ii)
                   dkl_vec_iii(istate2) = dkl(iii)
                ENDDO
                !
                ! First term
                ! Rescale the ground state orbitals by (dk-dl)_ii
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo,mo_coeff=mo_coeff)
                CALL cp_fm_to_fm(mo_coeff,fm_work_ii(ispin)%matrix,error=error)
                CALL cp_fm_column_scale(fm_work_ii(ispin)%matrix,dkl_vec_ii(1:homo))
                !
                ! Apply the p_iii operator
                ! fm_work_iii = -p_iii * (dk-dl)_ii * C0
                CALL cp_sm_fm_multiply(epr_env%op_p_ao(iii)%matrix,fm_work_ii(ispin)%matrix,&
                     &                 fm_work_iii(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,&
                     &                 error=error)
                !
                ! Copy in h1_psi0
                ! h1_psi0_i = fm_work_iii
                CALL cp_fm_to_fm(fm_work_iii(ispin)%matrix,h1_psi0(ispin)%matrix,error=error)
                !
                ! Second term
                ! Rescale the ground state orbitals by (dk-dl)_iii
                CALL cp_fm_to_fm(mo_coeff,fm_work_iii(ispin)%matrix,error=error)
                CALL cp_fm_column_scale(fm_work_iii(ispin)%matrix,dkl_vec_iii(1:homo))
                !
                ! Apply the p_ii operator
                ! fm_work_ii = -p_ii * (dk-dl)_iii * C0
                CALL cp_sm_fm_multiply(epr_env%op_p_ao(ii)%matrix,fm_work_iii(ispin)%matrix,&
                     &                 fm_work_ii(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,&
                     &                 error=error)
                !
                ! Copy in h1_psi0
                ! h1_psi0_i = fm_work_iii - fm_work_ii
                CALL cp_fm_scale_and_add(1.0_dp,h1_psi0(ispin)%matrix,&
                     &                  -1.0_dp,fm_work_ii(ispin)%matrix,&
                     &                   error=error)
             ENDDO
             !
             ! Optimize the response wavefunctions
             CALL linres_scf(p_env,qs_env,psi1,h1_psi0,error=error)
             !
             IF(output_unit>0) THEN
                WRITE(output_unit,'(T10,A,I4,A,E24.16)') 'Second order energy D_'//ACHAR(idir+119)//&
                     & ' for -state- ',istate,'=',p_env%etotal
                WRITE(output_unit,"(T10,A,/)")&
                     & "Store the psi1 vector for the calculation of the response current density "
             ENDIF
             !
             DO ispin = 1,dft_control%nspins
                !
                ! take care that no overflow can occur for uks
                IF(istate.GT.nstates(ispin)) CYCLE
                !
                ! need to reset those guys
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo)
                ist_true = statetrueindex(idir,istate,ispin)
                !
                ! the optimized wfns are copied in the fm
                CALL cp_fm_get_submatrix(psi1(ispin)%matrix,vecbuf_dklxp0 ,&
                     &                   1,ist_true,nao,1,transpose=.TRUE.,&
                     &                   error=error)
                CALL cp_fm_set_submatrix(psi1_D(ispin,idir)%matrix,vecbuf_dklxp0,&
                     &                   1,ist_true,nao,1,transpose=.TRUE.,&
                     &                   error=error)
                epr_env%fullepr_done(idir*ispin,istate) = .TRUE.
             ENDDO ! ispin
             !
          ENDDO ! istate
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             DO ispin = 1,dft_control%nspins
                CALL qs_print_cubes(qs_env,psi1_D(ispin,idir)%matrix,&
                     &              ncubes,list_cubes,&
                     &              centers_set(ispin)%array,ionode,'psi1_D',&
                     &              idir=idir,ispin=ispin,stride=section_get_ival(epr_section,&
                     &              "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE",error=error),&
                     &              error=error)
             ENDDO
          ENDIF ! print response functions
          ! CALL epr_restart_epr(epr_env,qs_env,action="write",ivec=idir+6,error=error)
       ENDDO ! idir
       !
       ! clean up
       DO ispin = 1,dft_control%nspins
          CALL cp_fm_release(fm_work_ii(ispin)%matrix,error=error)
          CALL cp_fm_release(fm_work_iii(ispin)%matrix,error=error)
       ENDDO
       DEALLOCATE(fm_work_ii,fm_work_iii,dkl_vec_ii,dkl_vec_iii,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(vecbuf_dklxp0,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !DO ispin = 1,dft_control%nspins
       !   DO idir = 1,3
       !      CALL cp_fm_set_all(psi1_D(ispin,idir)%matrix,0.0_dp,error=error)
       !   ENDDO
       !ENDDO
    ENDIF
    !
    !If current density and shift are calculated here
    ! probably a lot of memory is required
    ! Better to free what is not anymore necessary in p_env

    CALL p_env_release(p_env,error=error)
    !
    !Calculate current density and induced field
    CALL epr_response_current(epr_env,qs_env,psi1,h1_psi0,error=error)
    !
    ! clean up
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, psi1,error=error)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, h1_psi0,error=error)
    NULLIFY(psi1,h1_psi0)
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    CALL timestop(0.0_dp, handle)
    !
  END SUBROUTINE epr_do_epr

! *****************************************************************************
!!****f* cp2k/qs_linres_epr_utils/epr_env_init
!!
!!   NAME
!!     epr_env_init
!!
!!   FUNCTION
!!     Initialize the epr environment
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [MI]
!!
!****************************************************************************

  SUBROUTINE epr_env_init(epr_env,qs_env,error)
    !
    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i_B, icount, idir, ini, ir, ispin, istat, &
      istate, istate2, istate_next, j, k, n_mo(2), n_rep, nao, natom, ncubes, &
      nmo, nmoloc, nspins, nstates, output_unit
    INTEGER, DIMENSION(:), POINTER           :: bounds, list
    LOGICAL                                  :: failure, gapw
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: state_done
    REAL(dp)                                 :: center(3), center2(3), dist, &
                                                mdist, rab(3)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: vks_atom_set
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: epr_section, lr_section

    CALL timeset("epr_env_init","I"," ",handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, cell, dft_control, linres_control, scf_control)
    NULLIFY(logger, matrix_s, mos, mpools, epr_section, particle_set)
    NULLIFY(auxbas_pw_pool, pw_env)
    NULLIFY(jrho1_atom_set)

    n_mo(1:2) = 0
    nao = 0
    nmoloc = 0

    logger => cp_error_get_logger(error)
    !ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         &                             extension=".linresLog",error=error)
    
    IF(epr_env%ref_count /= 0) THEN
      CALL epr_env_cleanup(epr_env,error=error)
    END IF

    IF(output_unit>0) THEN
      WRITE(output_unit,"(/,T20,A,/)") "*** Start EPR g tensor calculation ***"
      WRITE(output_unit,"(T10,A,/)") "Initialization of the EPR environment"
    ENDIF

    CALL epr_env_create(epr_env,error=error)
    !
    ! If surrent_density or full_epr different allocations are required
    epr_section => section_vals_get_subs_vals(qs_env%input, &
         &                                    "PROPERTIES%LINRES%EPR",error=error)
    CALL section_vals_val_get(epr_section,"FULL",l_val=epr_env%full_epr,error=error)
    CALL section_vals_val_get(epr_section,"RESTART_EPR",l_val=epr_env%restart_epr,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         &          atomic_kind_set=atomic_kind_set,&
         &          cell=cell,&
         &          dft_control=dft_control,&
         &          linres_control=linres_control,&
         &          matrix_s=matrix_s,&
         &          mos=mos,&
         &          mpools=mpools,&
         &          particle_set=particle_set,&
         &          pw_env=pw_env,&
         &          scf_control=scf_control,error=error)
    !
    ! Check if restat also psi0 should be restarted
    IF(epr_env%restart_epr .AND. scf_control%density_guess/=restart_guess)THEN
       CALL stop_program("epr_env_init","restart_epr requires density_guess=restart")
    ENDIF
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure .AND. (output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get EPR parameters within PBC you need localized zero order orbitals '
    ENDIF
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    ALLOCATE(epr_env%centers_set(nspins),STAT=istat)
    nstates = 0
    epr_env%nstates = 0
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=n_mo(ispin),homo=homo)
       nmoloc = SIZE(linres_control%localized_wfn_control%centers_set(ispin)%array,2)
       CPPrecondition(nmoloc==homo,cp_warning_level,routineP,error,failure)
       ALLOCATE(epr_env%centers_set(ispin)%array(3,nmoloc),STAT=istat)
       CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       ! point to the psi0 centers
       DO idir = 1,3
          epr_env%centers_set(ispin)%array(idir,:) = &
               & linres_control%localized_wfn_control%centers_set(ispin)%array(idir,:)
       ENDDO
       epr_env%nstates(ispin) = homo
       nstates = MAX(nstates,homo)
    ENDDO
    epr_env%nao =nao
    ALLOCATE(epr_env%statetrueindex(3,nstates,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    epr_env%statetrueindex = 0
    ALLOCATE(state_done(3,nstates),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspins
       state_done(:,:) = .FALSE.
       epr_env%statetrueindex(1,1,ispin) = 1
       center(1) = epr_env%centers_set(ispin)%array(1,1)
       center(2) = epr_env%centers_set(ispin)%array(2,1)
       center(3) = epr_env%centers_set(ispin)%array(3,1)
       state_done(1,1) = .TRUE.
       icount = 1
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo)
       !
       DO idir = 1,3
          ini = 1
          IF(idir == 1) ini = 2
          DO istate = ini,nmo
             mdist = 100.0_dp
             !
             DO istate2 = 1,nmo
                IF(.NOT.state_done(idir,istate2)) THEN
                   center2(1) = epr_env%centers_set(ispin)%array(1,istate2)
                   center2(2) = epr_env%centers_set(ispin)%array(2,istate2)
                   center2(3) = epr_env%centers_set(ispin)%array(3,istate2)
                   !
                   rab = pbc(center,center2,cell)
                   CALL set_vecp(idir,j,k)
                   dist = SQRT(rab(j)*rab(j)+rab(k)*rab(k))
                   !
                   IF(dist .LT. mdist) THEN
                      mdist =dist
                      istate_next = istate2
                   ENDIF
                ENDIF
             ENDDO ! istate2
             !
             icount = icount + 1
             state_done(idir,istate_next) = .TRUE.
             epr_env%statetrueindex(idir,icount,ispin) = istate_next
             center(1) = epr_env%centers_set(ispin)%array(1,istate_next)
             center(2) = epr_env%centers_set(ispin)%array(2,istate_next)
             center(3) = epr_env%centers_set(ispin)%array(3,istate_next)
          ENDDO ! istate
          icount = 0
       ENDDO ! idir
    ENDDO
    DEALLOCATE(state_done,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    ! Conversion factors
    ! Magical constant twopi/cell%deth just like in NMR shift (basically undo scale_fac in qs_linres_nmr_current.F)
    epr_env%g_free_factor = -1.0_dp * e_gfactor
    epr_env%g_zke_factor = e_gfactor * ( a_fine )**2
    epr_env%g_so_factor = ( a_fine )**2 * ( -1.0_dp * e_gfactor - 1.0_dp ) * twopi / cell%deth
    epr_env%g_so_factor_gapw =  ( a_fine )**2 * ( -1.0_dp * e_gfactor - 1.0_dp )
    ! * 2 because B_ind = 2 * B_beta
    epr_env%g_soo_factor = 4.0_dp * fourpi * ( a_fine )**2 * twopi / cell%deth
    ! 2 * 2 * 1/4 * e^2 / m * a_0^2 * 2/3 * mu_0 / (omega * 1e-30 )
    epr_env%g_soo_chicorr_factor = 2.0_dp / 3.0_dp * e_charge**2 / e_mass * a_bohr**2 * &
                                   mu_perm / (cell%deth /cell%unit_of_length**3*1.0E-30_dp)

    IF(output_unit>0) THEN
       WRITE(output_unit,"(T2,A,T60,A)") "  EPR|","Gauge used "//TRIM(epr_env%gauge_name)
       IF(epr_env%full_epr) THEN
          WRITE(output_unit,"(T2,A,T60,A)") "  EPR|","Full orbital dependent correction"
       ENDIF
    ENDIF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
         &    "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN

       NULLIFY(bounds,list)
       ncubes = 0
       CALL section_vals_val_get(epr_section,&
            &                    "PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LU_BOUNDS",&
            &                    i_vals=bounds,error=error)
       ncubes = bounds(2) - bounds(1)  + 1
       IF(ncubes > 0 ) THEN
          ALLOCATE( epr_env%list_cubes(ncubes),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ir = 1,ncubes
             epr_env%list_cubes(ir) = bounds(1) + (ir-1)
          ENDDO
       ENDIF
       IF(.NOT. ASSOCIATED(epr_env%list_cubes)) THEN
          CALL section_vals_val_get(epr_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
               &                    n_rep_val=n_rep,error=error)
          ncubes = 0
          DO ir = 1,n_rep
             NULLIFY(list)
             CALL section_vals_val_get(epr_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
                  &                    i_rep_val=ir,i_vals=list,error=error)
             IF(ASSOCIATED(list)) THEN
                CALL reallocate(epr_env%list_cubes,1,ncubes+ SIZE(list))
                DO ini = 1, SIZE(list)
                   epr_env%list_cubes(ini+ncubes) = list(ini)
                ENDDO
                ncubes = ncubes + SIZE(list)
             ENDIF
          ENDDO  ! ir
       ENDIF
       IF(.NOT. ASSOCIATED(epr_env%list_cubes)) THEN
          ALLOCATE( epr_env%list_cubes(nstates),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ir = 1,nstates
             epr_env%list_cubes(ir) = ir
          ENDDO
       ENDIF
    ENDIF

    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=epr_env%ao_mo_fm_pools, error=error)
    ! for the chemical shift we need 6 psi1, i.e. 6 optimization procedures
    ! They become 9 if full epr is calculated, i.e. with the correction term too
    ! All of them are required at the end of the optimization procedure
    ! if the current density and the induced fiels have to be calculated
    ! If instead only the shift is needed, only one psi1 should be enough, providing
    ! that after every optimization the corresponding shift contribution is calculated
    ! prepare the psi1

    ALLOCATE(epr_env%psi1_p(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(epr_env%psi1_rxp(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       DO idir = 1,3
          NULLIFY(epr_env%psi1_p(ispin,idir)%matrix)
          CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
               &                 epr_env%psi1_p(ispin,idir)%matrix,error=error)
          NULLIFY(epr_env%psi1_rxp(ispin,idir)%matrix)
          CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
               &                 epr_env%psi1_rxp(ispin,idir)%matrix,error=error)
       ENDDO
    ENDDO

    IF(epr_env%full_epr) THEN
       ALLOCATE(epr_env%psi1_D(nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
          DO idir = 1,3
             NULLIFY(epr_env%psi1_D(ispin,idir)%matrix)
             CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                  &                 epr_env%psi1_D(ispin,idir)%matrix,error=error)
          ENDDO
       ENDDO
    ENDIF
    !
    ! px py pz operator acting on the atomic orbitals
    CALL allocate_matrix_set(epr_env%op_p_ao,3,error=error)
    CALL set_up_op_sm(epr_env%op_p_ao(1)%matrix,qs_env,symmetry="antisymmetric",&
         &            name="OP_P",error=error)
    CALL set_matrix(epr_env%op_p_ao(1)%matrix,0.0_dp)
    DO idir = 2,3
       CALL replicate_matrix_structure(epr_env%op_p_ao(1)%matrix, &
            &                          epr_env%op_p_ao(idir)%matrix,&
            &                          "epr_env%op_p_ao"//"-"//TRIM(&
            &                          ADJUSTL(cp_to_string(idir))),&
            &                          target_symmetry="antisymmetric",&
            &                          error=error)
       CALL set_matrix(epr_env%op_p_ao(idir)%matrix,0.0_dp)
    ENDDO
    !
    ! for the rxp we cannot calculate it a priori because it is in facts (r-dk)xp
    ! where dk is the center of the orbital it is applied to. We would need nstate operators
    ! What we can store is (r-dk)xp|psi0k> for each k, which is a full matrix only
    ! Therefore we prepare here the full matrix p_psi0 and rxp_psi0
    ! We also need a temporary sparse matrix where to store the integrals during the calculation
    ALLOCATE(epr_env%p_psi0(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(epr_env%rxp_psi0(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       DO idir = 1,3
          NULLIFY(epr_env%p_psi0(ispin,idir)%matrix)
          CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
               &                 epr_env%p_psi0(ispin,idir)%matrix,error=error)
          NULLIFY(epr_env%rxp_psi0(ispin,idir)%matrix)
          CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
               &                 epr_env%rxp_psi0(ispin,idir)%matrix,error=error)
       ENDDO
    ENDDO

    CALL allocate_matrix_set(epr_env%op_rmd_ao,3,error=error)
    CALL set_up_op_sm(epr_env%op_rmd_ao(1)%matrix,qs_env,symmetry="none",&
         &            name="OP_RMD",error=error)
    CALL set_matrix(epr_env%op_rmd_ao(1)%matrix,0.0_dp)
    DO idir = 2,3
       CALL replicate_matrix_structure(epr_env%op_rmd_ao(1)%matrix, &
            &                          epr_env%op_rmd_ao(idir)%matrix,&
            &                          "epr_env%op_rmd_ao"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(idir))),&
            &                          target_symmetry="none",error=error)
       CALL set_matrix(epr_env%op_rmd_ao(idir)%matrix,0.0_dp)
    ENDDO
    !
    ! Current density matrix in x y and z:
    !
    !jp0_ao
    ALLOCATE(epr_env%jp0_ao(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       NULLIFY(epr_env%jp0_ao(ispin)%matrix)
       CALL replicate_matrix_structure(epr_env%op_rmd_ao(1)%matrix, &
            &                          epr_env%jp0_ao(ispin)%matrix,"epr_env%jp0_ao"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       CALL set_matrix(epr_env%jp0_ao(ispin)%matrix,0.0_dp)
    ENDDO
    !
    !jp1_ao
    ALLOCATE( epr_env%jp1_ao(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       NULLIFY(epr_env%jp1_ao(ispin)%matrix)
       CALL replicate_matrix_structure(epr_env%op_rmd_ao(1)%matrix, &
            &                          epr_env%jp1_ao(ispin)%matrix,"epr_env%jp1_ao"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       CALL set_matrix(epr_env%jp1_ao(ispin)%matrix,0.0_dp)
    ENDDO
    !
    !jp2_ao
    ALLOCATE( epr_env%jp2_ao(2,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,2
       DO ispin = 1,nspins
          NULLIFY(epr_env%jp2_ao(idir,ispin)%matrix)
          CALL replicate_matrix_structure(epr_env%op_rmd_ao(1)%matrix, &
               &                          epr_env%jp2_ao(idir,ispin)%matrix,"epr_env%jp2_ao"//&
               &                          "-"//TRIM(ADJUSTL(cp_to_string(idir)))//&
               &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
               &                          target_symmetry="none",error=error)
       ENDDO
    ENDDO
    !
    ! If the current density on the grid needs to be stored
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
      ALLOCATE(epr_env%jrho1_set(3,3),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i_B = 1,3
         DO idir = 1,3
            NULLIFY(epr_env%jrho1_set(idir,i_B)%rho)
            ALLOCATE(epr_env%jrho1_set(idir,i_B)%rho,STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            epr_env%jrho1_set(idir,i_B)%rho%ref_count = 1
            NULLIFY( epr_env%jrho1_set(idir,i_B)%rho%rho_r)
            NULLIFY( epr_env%jrho1_set(idir,i_B)%rho%rho_g)

            ALLOCATE(epr_env%jrho1_set(idir,i_B)%rho%rho_r(nspins),stat=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            ALLOCATE(epr_env%jrho1_set(idir,i_B)%rho%rho_g(nspins),stat=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            DO ispin=1,nspins
               CALL pw_pool_create_pw(auxbas_pw_pool,&
                    &                  epr_env%jrho1_set(idir,i_B)%rho%rho_r(ispin)%pw,&
                    &                  use_data=REALDATA3D,in_space=REALSPACE,error=error)
               CALL pw_zero(epr_env%jrho1_set(idir,i_B)%rho%rho_r(ispin)%pw)

               CALL pw_pool_create_pw(auxbas_pw_pool,&
                    &                  epr_env%jrho1_set(idir,i_B)%rho%rho_g(ispin)%pw,&
                    &                  use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                    &                  error=error)
               CALL pw_zero(epr_env%jrho1_set(idir,i_B)%rho%rho_g(ispin)%pw)
            ENDDO
            CALL qs_rho_retain(epr_env%jrho1_set(idir,i_B)%rho,error=error)
         ENDDO
      ENDDO
    !
    ! Initialize local current density if GAPW calculation
    IF(gapw) THEN
       !write(*,*) 'epr_env_init: enter list_3c_jrho_atom'
       CALL list_3c_jrho_atom(qs_env,error=error)
       !write(*,*) 'epr_env_init: enter init_jrho_atom_set'
       CALL init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,&
            &                  nspins,.TRUE.,&
            &                  error=error)
       CALL init_nablavks_atom_set(nablavks_atom_set,atomic_kind_set,&
                   nspins,error=error)
       !write(*,*) 'epr_env_init: enter set_epr_env'
       CALL set_epr_env(epr_env=epr_env,jrho1_atom_set=jrho1_atom_set,&
                        nablavks_atom_set=nablavks_atom_set,&
            &           error=error)
       !write(*,*) 'epr_env_init: done with gapw'
    ENDIF
    !

    ! Bind
    ALLOCATE(epr_env%bind_set(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i_B = 1,3
      DO idir = 1,3
        NULLIFY(epr_env%bind_set(idir,i_B)%rho)
        ALLOCATE(epr_env%bind_set(idir,i_B)%rho,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        epr_env%bind_set(idir,i_B)%rho%ref_count = 1
        NULLIFY( epr_env%bind_set(idir,i_B)%rho%rho_r)
        NULLIFY( epr_env%bind_set(idir,i_B)%rho%rho_g)

        ALLOCATE(epr_env%bind_set(idir,i_B)%rho%rho_r(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(epr_env%bind_set(idir,i_B)%rho%rho_g(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%bind_set(idir,i_B)%rho%rho_r(1)%pw,&
             use_data=REALDATA3D,in_space=REALSPACE,error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%bind_set(idir,i_B)%rho%rho_g(1)%pw,&
             use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
      END DO
    END DO

    ! Nabla_V_ks
    ALLOCATE(epr_env%nablavks_set(3,dft_control%nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
      DO ispin = 1,nspins
        NULLIFY(epr_env%nablavks_set(idir,ispin)%rho)
        ALLOCATE(epr_env%nablavks_set(idir,ispin)%rho,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        epr_env%nablavks_set(idir,ispin)%rho%ref_count = 1
        NULLIFY( epr_env%nablavks_set(idir,ispin)%rho%rho_r)
        NULLIFY( epr_env%nablavks_set(idir,ispin)%rho%rho_g)
        ALLOCATE(epr_env%nablavks_set(idir,ispin)%rho%rho_r(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(epr_env%nablavks_set(idir,ispin)%rho%rho_g(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%nablavks_set(idir,ispin)%rho%rho_r(1)%pw,&
             use_data=REALDATA3D,in_space=REALSPACE,error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%nablavks_set(idir,ispin)%rho%rho_g(1)%pw,&
             use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
      END DO
    END DO

    ! Initialize the g tensor components
    ALLOCATE(epr_env%g_total(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(epr_env%g_so(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(epr_env%g_soo(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    epr_env%g_total = 0.0_dp
    epr_env%g_zke = 0.0_dp
    epr_env%g_so = 0.0_dp
    epr_env%g_soo = 0.0_dp
    ! Initialize the Bind0 tensor components
    ALLOCATE(epr_env%bind0(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    epr_env%bind0 = 0.0_dp

    ALLOCATE(epr_env%basisfun_center(3,epr_env%nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    epr_env%basisfun_center = 0.0_dp
    epr_env%simpleepr_done(1:6) = .FALSE.
    ALLOCATE(epr_env%fullepr_done(3*nspins,nstates),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    epr_env%fullepr_done = .FALSE.
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE epr_env_init

! *****************************************************************************
!!****f* cp2k/qs_linres_epr_utils/epr_env_cleanup
!!
!!   NAME
!!     epr_env_cleanup
!!
!!   FUNCTION
!!     Deallocate the epr environment
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE epr_env_cleanup(epr_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_B, idir, ispin, istat, j
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(.NOT. failure) THEN
      epr_env%ref_count = epr_env%ref_count - 1
      IF(epr_env%ref_count == 0 ) THEN
        !psi1_p
        IF(ASSOCIATED(epr_env%psi1_p)) THEN
          DO idir = 1,SIZE(epr_env%psi1_p,2)
            DO ispin = 1,SIZE(epr_env%psi1_p,1)
              CALL fm_pool_give_back_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                   epr_env%psi1_p(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(epr_env%psi1_p, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !psi1_rxp
        IF(ASSOCIATED(epr_env%psi1_rxp)) THEN
          DO idir = 1,SIZE(epr_env%psi1_rxp,2)
            DO ispin = 1,SIZE(epr_env%psi1_rxp,1)
              CALL fm_pool_give_back_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                   epr_env%psi1_rxp(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(epr_env%psi1_rxp, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !psi1_D
        IF(ASSOCIATED(epr_env%psi1_D)) THEN
          DO idir = 1,SIZE(epr_env%psi1_D,2)
            DO ispin = 1,SIZE(epr_env%psi1_D,1)
              CALL fm_pool_give_back_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                   epr_env%psi1_D(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(epr_env%psi1_D, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !p_psi0
        IF(ASSOCIATED(epr_env%p_psi0)) THEN
          DO idir = 1,SIZE(epr_env%p_psi0,2)
            DO ispin = 1,SIZE(epr_env%p_psi0,1)
              CALL fm_pool_give_back_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                   epr_env%p_psi0(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(epr_env%p_psi0, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !rxp_psi0
        IF(ASSOCIATED(epr_env%rxp_psi0)) THEN
          DO idir = 1,SIZE(epr_env%rxp_psi0,2)
            DO ispin = 1,SIZE(epr_env%rxp_psi0,1)
              CALL fm_pool_give_back_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                   epr_env%rxp_psi0(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(epr_env%rxp_psi0, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        NULLIFY(epr_env%ao_mo_fm_pools)

        DO ispin = 1,SIZE(epr_env%centers_set,1)
          DEALLOCATE(epr_env%centers_set(ispin)%array,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO
        DEALLOCATE(epr_env%centers_set,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        IF(ASSOCIATED(epr_env%list_cubes)) THEN
          DEALLOCATE(epr_env%list_cubes)
        END IF
        !op_p_ao
        IF(ASSOCIATED(epr_env%op_p_ao)) THEN
          CALL deallocate_matrix_set(epr_env%op_p_ao,error=error)
        END IF
        !op_rmd_ao
        IF(ASSOCIATED(epr_env%op_rmd_ao)) THEN
          CALL deallocate_matrix_set(epr_env%op_rmd_ao,error=error)
        END IF
        !jp0_ao
        IF(ASSOCIATED(epr_env%jp0_ao)) THEN
          DO i = 1,SIZE(epr_env%jp0_ao,1)
              CALL deallocate_matrix(epr_env%jp0_ao(i)%matrix,error=error)
          END DO
          DEALLOCATE(epr_env%jp0_ao,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF        
        !jp1_ao
        IF(ASSOCIATED(epr_env%jp1_ao)) THEN
          DO i = 1,SIZE(epr_env%jp1_ao,1)
              CALL deallocate_matrix(epr_env%jp1_ao(i)%matrix,error=error)
          END DO
          DEALLOCATE(epr_env%jp1_ao,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !jp2_ao
        IF(ASSOCIATED(epr_env%jp2_ao)) THEN
           DO i = 1,SIZE(epr_env%jp2_ao,1)
              DO j = 1,SIZE(epr_env%jp2_ao,2)
                 CALL deallocate_matrix(epr_env%jp2_ao(i,j)%matrix,error=error)
              END DO
           END DO
           DEALLOCATE(epr_env%jp2_ao,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! Current density on the grid
        IF(ASSOCIATED(epr_env%jrho1_set)) THEN
          DO i_B = 1,3
            DO idir = 1,3
              CALL qs_rho_release(epr_env%jrho1_set(idir,i_B)%rho,error=error)
            END DO
          END DO
        END IF
        ! Local current density, atom by atom (only gapw)
        IF(ASSOCIATED(epr_env%jrho1_atom_set)) THEN
          CALL deallocate_jrho_atom_set(epr_env%jrho1_atom_set,error=error)
        END IF

        !fullepr_done
        IF(ASSOCIATED(epr_env%fullepr_done)) THEN
          DEALLOCATE(epr_env%fullepr_done,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(epr_env%basisfun_center)) THEN
          DEALLOCATE(epr_env%basisfun_center,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(epr_env%statetrueindex)) THEN
          DEALLOCATE(epr_env%statetrueindex,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        ! nablavks_set
        IF(ASSOCIATED(epr_env%nablavks_set)) THEN
          DEALLOCATE(epr_env%nablavks_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! nablavks_atom_set
        IF(ASSOCIATED(epr_env%nablavks_atom_set)) THEN
          CALL deallocate_nablavks_atom_set(epr_env%nablavks_atom_set,error=error)
        END IF
        ! vks_atom_set
        IF(ASSOCIATED(epr_env%vks_atom_set)) THEN
          CALL deallocate_rho_atom_set(epr_env%vks_atom_set)
        END IF
        ! bind_set
        IF(ASSOCIATED(epr_env%bind_set)) THEN
          DEALLOCATE(epr_env%bind_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! bind_atom_set
        IF(ASSOCIATED(epr_env%bind_atom_set)) THEN
          DEALLOCATE(epr_env%bind_atom_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_total
        IF(ASSOCIATED(epr_env%g_total)) THEN
          DEALLOCATE(epr_env%g_total,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_so
        IF(ASSOCIATED(epr_env%g_so)) THEN
          DEALLOCATE(epr_env%g_so,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_soo
        IF(ASSOCIATED(epr_env%g_soo)) THEN
          DEALLOCATE(epr_env%g_soo,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! bind0
        IF(ASSOCIATED(epr_env%bind0)) THEN
          DEALLOCATE(epr_env%bind0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

      END IF  ! ref count
    END IF ! failure

  END SUBROUTINE epr_env_cleanup


END MODULE  qs_linres_epr_utils
