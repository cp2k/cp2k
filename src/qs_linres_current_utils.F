!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Chemical shift calculation by dfpt
!>      Initialization of the nmr_env, creation of the special neighbor lists
!>      Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!>      Write output
!>      Deallocate everything
!> \note
!>      The psi0 should be localized
!>      the Sebastiani method works within the assumption that the orbitals are
!>      completely contained in the simulation box
!> \par History
!>       created 07-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_current_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_i_utils,                ONLY: cp_2d_i_p_type
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
  USE input_constants,                 ONLY: current_gauge_r,&
                                             current_gauge_r_and_Bader,&
                                             current_gauge_r_and_step_func,&
                                             current_orb_center_atom,&
                                             current_orb_center_box,&
                                             current_orb_center_common,&
                                             current_orb_center_wannier
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_release
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_atom_current,          ONLY: list_3c_jrho_atom
  USE qs_linres_current,               ONLY: chi_analytic,&
                                             current_response_current
  USE qs_linres_methods,               ONLY: linres_scf
  USE qs_linres_op,                    ONLY: set_vecp
  USE qs_linres_types,                 ONLY: current_env_create,&
                                             current_env_type,&
                                             deallocate_jrho_atom_set,&
                                             get_current_env,&
                                             init_jrho_atom_set,&
                                             jrho_atom_type,&
                                             linres_control_type,&
                                             set_current_env
  USE qs_loc_methods,                  ONLY: qs_print_cubes
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: set_up_op_sm
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_rho_types,                    ONLY: qs_rho_retain
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
  PUBLIC :: current_do_current, current_do_chi, current_env_cleanup, current_env_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_current_utils'

CONTAINS

! *****************************************************************************
  SUBROUTINE current_do_current(current_env,p_env,qs_env,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_do_current', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i, icenter, idir, ii, iii, ispin, ist_true, &
      istat, istate, j, jcenter, jstate, max_nbr_center, max_states, nao, &
      nbr_center(2), ncubes, nmo, nspins, nstates(2), output_unit
    INTEGER, DIMENSION(:), POINTER           :: list_cubes
    INTEGER, DIMENSION(:, :, :), POINTER     :: statetrueindex
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: dk(3), dkl(3), dl(3)
    REAL(dp), DIMENSION(:), POINTER          :: dkl_vec_ii, dkl_vec_iii
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_dklxp0
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
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
    TYPE(section_vals_type), POINTER         :: current_section, lr_section, &
                                                print_key

!
!

    failure   = .FALSE.
    CALL timeset(routineN,"I"," ",handle)
    !
    NULLIFY(ao_mo_fm_pools, cell, dft_control, linres_control, lr_section, current_section,&
         &  logger, mpools, psi1,h1_psi0, vecbuf_dklxp0, mos, mo_coeff,para_env, &
         &  fm_work_ii, fm_work_iii, tmp_fm_struct, dkl_vec_ii,dkl_vec_iii, &
         &  list_cubes, statetrueindex, centers_set, center_list, psi1_p, psi1_rxp, psi1_D, &
         &  p_psi0, rxp_psi0)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    current_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%CURRENT",error=error)

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

    nspins = dft_control%nspins

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    CALL mpools_get(mpools, ao_mo_fm_pools=ao_mo_fm_pools,error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, psi1, name=routineP//":psi1",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, h1_psi0, name=routineP//":h1_psi0",error=error)

    CALL get_current_env(current_env=current_env, nao=nao, nstates=nstates, &
                         centers_set=centers_set, nbr_center=nbr_center, center_list=center_list, &
                         list_cubes=list_cubes, statetrueindex=statetrueindex, &
                         psi1_p=psi1_p, psi1_rxp=psi1_rxp, psi1_D=psi1_D,&
                         p_psi0=p_psi0, rxp_psi0=rxp_psi0,error=error)
    !
    ! operator p
    DO idir = 1,3
       CALL cp_fm_set_all(psi1_p(1,idir)%matrix,0.0_dp,error=error)
    ENDDO
    DO idir = 1,3
       IF(.NOT.current_env%simple_done(idir)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator P_"//ACHAR(idir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(p_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(p_psi0(ispin,idir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          !
          linres_control%converged = .FALSE.
          CALL linres_scf(p_env,qs_env,psi1,h1_psi0,error=error)
          !
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
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_p(ispin,idir)%matrix,error=error)
          ENDDO
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             print_key => section_vals_get_subs_vals(current_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
             DO ispin = 1,nspins
                CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
                      centers_set(ispin)%array,print_key,'psi1_p',&
                      idir=idir,ispin=ispin,error=error)
             ENDDO ! ispin
          ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          ! Write the result in the restart file
          current_env%simple_done(idir) = .TRUE.
          current_env%simple_converged(idir) = linres_control%converged
          !CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir,error=error)
       ENDIF
    ENDDO ! idir
    !
    ! operator rxp
    DO idir = 1,3
       CALL cp_fm_set_all(psi1_rxp(1,idir)%matrix,0.0_dp,error=error)
    ENDDO
    
    DO idir = 1,3
       IF(.NOT.current_env%simple_done(idir+3)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator L_"//ACHAR(idir+119)
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
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_rxp(ispin,idir)%matrix,error=error)
          ENDDO
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             print_key => section_vals_get_subs_vals(current_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
             DO ispin = 1,nspins
                CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
                     centers_set(ispin)%array,print_key,'psi1_rxp',&
                     idir=idir,ispin=ispin,error=error)
             ENDDO ! ispin
          ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          ! Write the result inthe restart file
          current_env%simple_done(idir+3) = .TRUE.
          !CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir+3,error=error)
       ENDIF
    ENDDO ! idir
    !
    ! operator D
    IF(current_env%full) THEN
       DO ispin = 1,nspins
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
       ALLOCATE(fm_work_ii(nspins),fm_work_iii(nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
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
       max_nbr_center = MAXVAL(nbr_center(1:nspins))
       max_states = MAXVAL(nstates(1:nspins))
       !
       ALLOCATE(dkl_vec_ii(max_states),dkl_vec_iii(max_states),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !
       DO idir = 1,3
          !
          ! here we run over the max number of states, then we need 
          ! to be careful with overflow while doing uks calculations.
          DO icenter = 1,max_nbr_center
             !
             IF(output_unit>0) THEN
                WRITE(output_unit,"(T10,A,I4,A)")&
                     & "Response to the perturbation operator (dk-dl)xp for -state- ",&
                     & icenter," in dir. "//ACHAR(idir+119)
             ENDIF
             !
             DO ispin = 1,nspins
                !
                ! take care that no overflow can occur for uks
                IF(icenter.GT.nbr_center(ispin)) THEN
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
                ist_true = statetrueindex(idir,icenter,ispin)
                !
                ! the initial guess is the previous set of psi1, just optimized
                CALL set_vecp(idir,ii,iii)
                dk(1:3) = centers_set(ispin)%array(1:3,ist_true)
                !
                DO jcenter = 1,nbr_center(ispin)
                   dl(1:3) = centers_set(ispin)%array(1:3,jcenter)
                   dkl = pbc(dl,dk,cell)
                   DO j = center_list(ispin)%array(1,jcenter),center_list(ispin)%array(1,jcenter+1)-1
                      jstate = center_list(ispin)%array(2,j)
                      dkl_vec_ii(jstate)  = dkl(ii)
                      dkl_vec_iii(jstate) = dkl(iii)
                   ENDDO
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
                CALL cp_sm_fm_multiply(current_env%op_p_ao(iii)%matrix,fm_work_ii(ispin)%matrix,&
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
                CALL cp_fm_column_scale(fm_work_iii(ispin)%matrix,dkl_vec_iii(1:nmo))
                !
                ! Apply the p_ii operator
                ! fm_work_ii = -p_ii * (dk-dl)_iii * C0
                CALL cp_sm_fm_multiply(current_env%op_p_ao(ii)%matrix,fm_work_iii(ispin)%matrix,&
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
                     & ' for -state- ',icenter,'=',p_env%etotal
                WRITE(output_unit,"(T10,A,/)")&
                     & "Store the psi1 vector for the calculation of the response current density "
             ENDIF
             !
             DO ispin = 1,nspins
                !
                ! take care that no overflow can occur for uks
                IF(icenter.GT.nbr_center(ispin)) CYCLE
                !
                ! need to reset those guys
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo)
                ist_true = statetrueindex(idir,icenter,ispin)
                DO i = center_list(ispin)%array(1,ist_true),center_list(ispin)%array(1,ist_true+1)-1
                   istate = center_list(ispin)%array(2,i)
                   !
                   ! the optimized wfns are copied in the fm
                   CALL cp_fm_get_submatrix(psi1(ispin)%matrix,vecbuf_dklxp0 ,&
                        &                   1,istate,nao,1,transpose=.TRUE.,&
                        &                   error=error)
                   CALL cp_fm_set_submatrix(psi1_D(ispin,idir)%matrix,vecbuf_dklxp0,&
                        &                   1,istate,nao,1,transpose=.TRUE.,&
                        &                   error=error)
                ENDDO
                current_env%full_done(idir*ispin,icenter) = .TRUE.
             ENDDO ! ispin
             !
          ENDDO ! istate
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             print_key => section_vals_get_subs_vals(current_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
             DO ispin = 1,nspins
                CALL qs_print_cubes(qs_env,psi1_D(ispin,idir)%matrix,&
                     ncubes,list_cubes,centers_set(ispin)%array,print_key,'psi1_D',&
                      idir=idir,ispin=ispin,error=error)
             ENDDO
          ENDIF ! print response functions
          ! CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir+6,error=error)
       ENDDO ! idir
       !
       ! clean up
       DO ispin = 1,nspins
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
    CALL current_response_current(current_env,qs_env,psi1,h1_psi0,error=error)
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
  END SUBROUTINE current_do_current

  SUBROUTINE current_do_chi(current_env,qs_env,error)
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_do_chi', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iB, iiB, iiiB, ispin, &
                                                nspins
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control

    failure = .FALSE.
    NULLIFY(dft_control)
    !
    CALL timeset(routineN,"I"," ",handle)
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          dft_control=dft_control,&
         &          error=error)
    nspins = dft_control%nspins
    DO iB = 1,3
       CALL set_vecp(iB,iiB,iiiB)
       DO ispin = 1,nspins
          CALL chi_analytic(current_env,qs_env,ispin,iB,error=error)
       ENDDO
    ENDDO
    !
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE current_do_chi

  SUBROUTINE current_env_init(current_env,qs_env,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i, i_B, iao, iatom, ibox, icount, idir, ii, ini, &
      ir, ispin, istat, istate, istate2, istate_next, ix, iy, iz, j, k, &
      max_nbr_center, max_states, n_mo(2), n_rep, nao, natom, nbr_box, &
      ncubes, nmo, nmoloc, nspins, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buff, first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: bounds, list, nbox
    LOGICAL                                  :: failure, gapw, is0
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: state_done
    REAL(dp)                                 :: center(3), center2(3), dist, &
                                                mdist, r(3), rab(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rbuff
    REAL(dp), DIMENSION(:), POINTER          :: common_center
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
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: current_section, lr_section

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, cell, dft_control, linres_control, scf_control,&
         &  logger, matrix_s, mos, mpools, current_section, particle_set,&
         &  auxbas_pw_pool, pw_env, jrho1_atom_set, common_center)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         &                             extension=".linresLog",error=error)
    
    IF(current_env%ref_count /= 0) THEN
       CALL current_env_cleanup(current_env,error=error)
    END IF

    IF(output_unit>0) THEN
       WRITE(output_unit,"(/,T20,A,/)") "*** Start current Calculation ***"
       WRITE(output_unit,"(T10,A,/)") "Inizialization of the current environment"
    ENDIF

    CALL current_env_create(current_env,error=error)
    !
    current_env%chi_tensor(:,:,:)     = 0.0_dp
    current_env%chi_tensor_loc(:,:,:) = 0.0_dp
    !
    ! If current_density or full_nmr different allocations are required
    current_section => section_vals_get_subs_vals(qs_env%input, &
         &         "PROPERTIES%LINRES%CURRENT",error=error)
    CALL section_vals_val_get(current_section,"GAUGE",i_val=current_env%gauge,error=error)
    SELECT CASE(current_env%gauge)
    CASE(current_gauge_r)
       current_env%gauge_name="R"
    CASE(current_gauge_r_and_Bader)
       current_env%gauge_name="R_AND_BADER"
    CASE(current_gauge_r_and_step_func)
       current_env%gauge_name="R_AND_STEP_FUNCTION"
    CASE DEFAULT
       CALL stop_program("current_env_init","unknown gauge, try again...")
    END SELECT

    ! no need of this guy anymore
    !CALL section_vals_val_get(current_section,"FULL",l_val=current_env%full,error=error)
    current_env%full=.TRUE.

    CALL section_vals_val_get(current_section,"ORBITAL_CENTER",i_val=current_env%orb_center,error=error)

    SELECT CASE(current_env%orb_center)
    CASE(current_orb_center_wannier)
       current_env%orb_center_name = "WANNIER"
    CASE(current_orb_center_common)
       current_env%orb_center_name = "COMMON"
       current_env%full = .FALSE.
       ! Is there a user specified common_center?
       CALL section_vals_val_get(current_section,"COMMON_CENTER",r_vals=common_center,error=error)
    CASE(current_orb_center_atom)
       current_env%orb_center_name = "ATOM"
    CASE(current_orb_center_box)
       current_env%orb_center_name = "BOX"
       ! Is there a user specified nbox?
       CALL section_vals_val_get(current_section,"NBOX",i_vals=nbox,error=error)
    CASE DEFAULT
       CALL stop_program("current_env_init","unknown orbital center, try again...")
    END SELECT

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
         &          scf_control=scf_control, &
         &          error=error)
    !
    ! Check if restat also psi0 should be restarted
    !IF(current_env%restart_current .AND. scf_control%density_guess/=restart_guess)THEN
    !   CALL stop_program("current_env_init","restart_nmr requires density_guess=restart")
    !ENDIF
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure .AND. (output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get CURRENT parameters within PBC you need localized zero order orbitals '
    ENDIF
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    ALLOCATE(current_env%centers_set(nspins),current_env%center_list(nspins),STAT=istat)
    CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
    !
    !
    IF(current_env%do_qmmm) THEN
       DO ispin = 1,nspins
          DO istate = 1,SIZE(linres_control%localized_wfn_control%centers_set(ispin)%array,2)
             ! just to be sure...
             r(:) = pbc(linres_control%localized_wfn_control%centers_set(ispin)%array(:,istate),cell)
             IF(r(1).LT.0.0_dp) THEN
                linres_control%localized_wfn_control%centers_set(ispin)%array(1,istate) = &
                     r(1) + cell%hmat(1,1)
             ENDIF
             IF(r(2).LT.0.0_dp) THEN
                linres_control%localized_wfn_control%centers_set(ispin)%array(2,istate) = &
                     r(2) + cell%hmat(2,2)
             ENDIF
             IF(r(3).LT.0.0_dp) THEN
                linres_control%localized_wfn_control%centers_set(ispin)%array(3,istate) = &
                     r(3) + cell%hmat(3,3)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    max_states = 0
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=n_mo(ispin),homo=homo)
       nmoloc = SIZE(linres_control%localized_wfn_control%centers_set(ispin)%array,2)
       CPPrecondition(nmoloc==homo,cp_warning_level,routineP,error,failure)
       current_env%nstates(ispin) = homo
       max_states = MAX(max_states,homo)
       !
       ALLOCATE(current_env%center_list(ispin)%array(2,homo+1),&
            &   current_env%centers_set(ispin)%array(3,homo),STAT=istat)
       CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       current_env%center_list(ispin)%array(:,:) = 0
       current_env%centers_set(ispin)%array(:,:) = HUGE(0.0_dp)
       !
       ! point to the psi0 centers
       SELECT CASE(current_env%orb_center)
       CASE(current_orb_center_wannier)
          !
          ! use the wannier center as -center-
          current_env%nbr_center(ispin) = homo
          DO istate = 1,homo
             current_env%centers_set(ispin)%array(1:3,istate) = &
                  & linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate)
             current_env%center_list(ispin)%array(:,istate) = istate
          ENDDO
          current_env%center_list(ispin)%array(1,homo+1) = homo+1
          !
       CASE(current_orb_center_common)
          !
          ! use a common -center-
          current_env%centers_set(ispin)%array(:,1) = common_center(:)
          current_env%nbr_center(ispin) = 1
          current_env%center_list(ispin)%array(1,1) = 1
          current_env%center_list(ispin)%array(1,2) = homo+1
          DO istate = 1,homo
             current_env%center_list(ispin)%array(2,istate) = istate
          ENDDO
          !
       CASE(current_orb_center_atom)
          !
          ! use the atom as -center-
          ALLOCATE(buff(homo),STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
          buff(:) = 0
          !
          DO istate = 1,homo
             mdist = HUGE(0.0_dp)
             DO iatom = 1,natom
                r   = pbc(particle_set(iatom)%r(:),cell)
                rab = pbc(r,linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate),&
                     &cell)
                dist = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)
                IF(dist.LT.mdist) THEN
                   buff(istate) = iatom
                   mdist = dist
                ENDIF
             ENDDO
          ENDDO
          !
          i = 0
          ii = 1
          current_env%center_list(ispin)%array(1,1) = 1
          DO iatom = 1,natom
             j = 0
             is0 = .TRUE.
             DO istate = 1,homo
                IF(buff(istate).EQ.iatom) THEN
                   j = j + 1
                   i = i + 1
                   is0 = .FALSE.
                   current_env%center_list(ispin)%array(2,i) = istate
                ENDIF
             ENDDO
             IF(.NOT.is0) THEN
                IF(output_unit>0) THEN
                   WRITE(output_unit,'(A,I6,A,I6)') 'clustering ',j,' center(s) on atom ',iatom
                ENDIF
                current_env%center_list(ispin)%array(1,ii+1) = &
                     & current_env%center_list(ispin)%array(1,ii) + j
                current_env%centers_set(ispin)%array(:,ii) = &
                     &  pbc(particle_set(iatom)%r,cell)
                ii = ii + 1
             ENDIF
          ENDDO
          current_env%nbr_center(ispin) = ii-1
          !
          DEALLOCATE(buff,STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       CASE(current_orb_center_box)
          !
          ! use boxes as -center-
          nbr_box = nbox(1)*nbox(2)*nbox(3)
          ALLOCATE(rbuff(3,nbr_box),buff(homo),STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
          rbuff(:,:) = HUGE(0.0_dp)
          buff(:) = 0
          !
          ibox = 1
          DO iz = 1,nbox(3)
             DO iy = 1,nbox(2)
                DO ix = 1,nbox(1)
                   rbuff(1,ibox) = cell%hmat(1,1)*((REAL(ix,KIND=dp)-0.5_dp)/REAL(nbox(1),KIND=dp)-0.5_dp)
                   rbuff(2,ibox) = cell%hmat(2,2)*((REAL(iy,KIND=dp)-0.5_dp)/REAL(nbox(2),KIND=dp)-0.5_dp)
                   rbuff(3,ibox) = cell%hmat(3,3)*((REAL(iz,KIND=dp)-0.5_dp)/REAL(nbox(3),KIND=dp)-0.5_dp)
                   ibox = ibox+1
                ENDDO
             ENDDO
          ENDDO
          !
          DO istate = 1,homo
             mdist = HUGE(0.0_dp)
             DO ibox = 1,nbr_box
                rab(:) = pbc(rbuff(:,ibox), &
                     &       linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate),&
                     &       cell)
                dist = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)
                IF(dist.LT.mdist) THEN
                   buff(istate) = ibox
                   mdist = dist
                ENDIF
             ENDDO
          ENDDO
          !
          i = 0
          ii = 1
          current_env%center_list(ispin)%array(1,1) = 1
          DO ibox = 1,nbr_box
             j = 0
             is0 = .TRUE.
             DO istate = 1,homo
                IF(buff(istate).EQ.ibox) THEN
                   j = j + 1
                   i = i + 1
                   is0 = .FALSE.
                   current_env%center_list(ispin)%array(2,i) = istate
                ENDIF
             ENDDO
             IF(.NOT.is0) THEN
                IF(output_unit>0) THEN
                   WRITE(output_unit,'(A,I6,A,I6)') 'clustering ',j,' center(s) on box ',ibox
                ENDIF
                current_env%center_list(ispin)%array(1,ii+1) = &
                  & current_env%center_list(ispin)%array(1,ii) + j
                current_env%centers_set(ispin)%array(:,ii) = rbuff(:,ibox)
                ii = ii + 1
             ENDIF
          ENDDO
          current_env%nbr_center(ispin) = ii-1
          !
          DEALLOCATE(buff,rbuff,STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       CASE DEFAULT
          CALL stop_program("current_env_init","unknown orbital center, try again...")
       END SELECT
    ENDDO
    !
    !
    IF(current_env%do_qmmm.AND.&
         & (current_env%orb_center.EQ.current_orb_center_wannier.OR.&
         &  current_env%orb_center.EQ.current_orb_center_atom.OR.&
         &  current_env%orb_center.EQ.current_orb_center_box)) THEN
       IF(output_unit>0) THEN
          WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING '
          WRITE(output_unit,*) 'orbital center shifted to match the '
          WRITE(output_unit,*) 'center of the box (L/2 L/2 L/2) (superdirty...)'
          WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING '
       ENDIF
       DO ispin = 1,nspins
          DO istate = 1,current_env%nbr_center(ispin)
             IF(current_env%centers_set(ispin)%array(1,istate).LE.0.0_dp) THEN
                current_env%centers_set(ispin)%array(1,istate) = &
                     current_env%centers_set(ispin)%array(1,istate) + cell%hmat(1,1)
             ENDIF
             IF(current_env%centers_set(ispin)%array(2,istate).LE.0.0_dp) THEN
                current_env%centers_set(ispin)%array(2,istate) = &
                     current_env%centers_set(ispin)%array(2,istate) + cell%hmat(2,2)
             ENDIF
             IF(current_env%centers_set(ispin)%array(3,istate).LE.0.0_dp) THEN
                current_env%centers_set(ispin)%array(3,istate) = &
                  current_env%centers_set(ispin)%array(3,istate) + cell%hmat(3,3)
             ENDIF
          ENDDO
       ENDDO
       ! printing
       IF(output_unit>0) THEN
          DO ispin = 1,nspins
             WRITE(output_unit,'(/,T2,A,I2)') "WANNIER CENTERS for spin ",ispin
             WRITE(output_unit,'(/,T18,A)') "--------------- Shifted Centers --------------- "
             DO istate = 1,current_env%nbr_center(ispin)
                WRITE(output_unit,'(T5,A,I6,3F12.6)')&
                     'state ',istate,current_env%centers_set(ispin)%array(1:3,istate)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    max_nbr_center = MAXVAL(current_env%nbr_center(1:nspins))
    current_env%nao = nao
    ALLOCATE(current_env%statetrueindex(3,max_nbr_center,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    current_env%statetrueindex(:,:,:) = 0
    ALLOCATE(state_done(3,max_nbr_center),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspins
       state_done(:,:) = .FALSE.
       current_env%statetrueindex(1,1,ispin) = 1
       center(1) = current_env%centers_set(ispin)%array(1,1)
       center(2) = current_env%centers_set(ispin)%array(2,1)
       center(3) = current_env%centers_set(ispin)%array(3,1)
       state_done(1,1) = .TRUE.
       icount = 1
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo)
       !
       DO idir = 1,3
          ini = 1
          IF(idir == 1) ini = 2
          DO istate = ini,current_env%nbr_center(ispin)
             mdist = HUGE(0.0_dp)
             !
             DO istate2 = 1,current_env%nbr_center(ispin)
                IF(.NOT.state_done(idir,istate2)) THEN
                   center2(1) = current_env%centers_set(ispin)%array(1,istate2)
                   center2(2) = current_env%centers_set(ispin)%array(2,istate2)
                   center2(3) = current_env%centers_set(ispin)%array(3,istate2)
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
             current_env%statetrueindex(idir,icount,ispin) = istate_next
             center(1) = current_env%centers_set(ispin)%array(1,istate_next)
             center(2) = current_env%centers_set(ispin)%array(2,istate_next)
             center(3) = current_env%centers_set(ispin)%array(3,istate_next)
          ENDDO ! istate
          icount = 0
       ENDDO ! idir
    ENDDO
    DEALLOCATE(state_done,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,"(T2,A)") "CURRENT| Gauge used "//TRIM(current_env%gauge_name)
       WRITE(output_unit,"(T2,A)") "CURRENT| Orbital center used "//TRIM(current_env%orb_center_name)
       IF(current_env%orb_center.EQ.current_orb_center_common) THEN
          WRITE(output_unit,"(T2,A,3F10.6)") "CURRENT| Common center ",common_center(1:3)
       ELSEIF(current_env%orb_center.EQ.current_orb_center_box) THEN
          WRITE(output_unit,"(T2,A,3I5)") "CURRENT| nbr boxes in each direction ",nbox(1:3)
       ENDIF
       !IF(current_env%full) THEN
       !   WRITE(output_unit,"(T2,A,T2,A)") "CURRENT|","Full orbital dependent correction"
       !ENDIF
    ENDIF

    IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
         &    "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN

       NULLIFY(bounds,list)
       ncubes = 0
       CALL section_vals_val_get(current_section,&
            &                    "PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LU_BOUNDS",&
            &                    i_vals=bounds,error=error)
       ncubes = bounds(2) - bounds(1)  + 1
       IF(ncubes > 0 ) THEN
          ALLOCATE( current_env%list_cubes(ncubes),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ir = 1,ncubes
             current_env%list_cubes(ir) = bounds(1) + (ir-1)
          ENDDO
       ENDIF
       IF(.NOT. ASSOCIATED(current_env%list_cubes)) THEN
          CALL section_vals_val_get(current_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
               &                    n_rep_val=n_rep,error=error)
          ncubes = 0
          DO ir = 1,n_rep
             NULLIFY(list)
             CALL section_vals_val_get(current_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
                  &                    i_rep_val=ir,i_vals=list,error=error)
             IF(ASSOCIATED(list)) THEN
                CALL reallocate(current_env%list_cubes,1,ncubes+ SIZE(list))
                DO ini = 1, SIZE(list)
                   current_env%list_cubes(ini+ncubes) = list(ini)
                ENDDO
                ncubes = ncubes + SIZE(list)
             ENDIF
          ENDDO  ! ir
       ENDIF
       IF(.NOT. ASSOCIATED(current_env%list_cubes)) THEN
          ALLOCATE(current_env%list_cubes(max_states),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ir = 1,max_states
             current_env%list_cubes(ir) = ir
          ENDDO
       ENDIF
    ENDIF
    IF (BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
         &    "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
    ENDIF

    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=current_env%ao_mo_fm_pools, error=error)
    ! for the chemical shift we need 6 psi1, i.e. 6 optimization procedures
    ! They become 9 if full nmr is calculated, i.e. with the correction term too
    ! All of them are required at the end of the optimization procedure
    ! if the current density and the induced fiels have to be calculated
    ! If instead only the shift is needed, only one psi1 should be enough, providing
    ! that after every optimization the corresponding shift contribution is calculated
    ! prepare the psi1

    ALLOCATE(current_env%psi1_p(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(current_env%psi1_rxp(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       DO idir = 1,3
          NULLIFY(current_env%psi1_p(ispin,idir)%matrix)
          CALL fm_pool_create_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
               &                 current_env%psi1_p(ispin,idir)%matrix,error=error)
          NULLIFY(current_env%psi1_rxp(ispin,idir)%matrix)
          CALL fm_pool_create_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
               &                 current_env%psi1_rxp(ispin,idir)%matrix,error=error)
       ENDDO
    ENDDO

    IF(current_env%full) THEN
       ALLOCATE(current_env%psi1_D(nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
          DO idir = 1,3
             NULLIFY(current_env%psi1_D(ispin,idir)%matrix)
             CALL fm_pool_create_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
                  &                 current_env%psi1_D(ispin,idir)%matrix,error=error)
          ENDDO
       ENDDO
    ENDIF
    !
    ! px py pz operator acting on the atomic orbitals
    CALL allocate_matrix_set(current_env%op_p_ao,3,error=error)
    CALL set_up_op_sm(current_env%op_p_ao(1)%matrix,qs_env,symmetry="antisymmetric",&
         &            name="OP_P",error=error)
    CALL set_matrix(current_env%op_p_ao(1)%matrix,0.0_dp)
    DO idir = 2,3
       CALL replicate_matrix_structure(current_env%op_p_ao(1)%matrix, &
            &                          current_env%op_p_ao(idir)%matrix,&
            &                          "current_env%op_p_ao"//"-"//TRIM(&
            &                          ADJUSTL(cp_to_string(idir))),&
            &                          target_symmetry="antisymmetric",&
            &                          error=error)
       CALL set_matrix(current_env%op_p_ao(idir)%matrix,0.0_dp)
    ENDDO
    !
    ! for the rxp we cannot calculate it a priori because it is in facts (r-dk)xp
    ! where dk is the center of the orbital it is applied to. We would need nstate operators
    ! What we can store is (r-dk)xp|psi0k> for each k, which is a full matrix only
    ! Therefore we prepare here the full matrix p_psi0 and rxp_psi0
    ! We also need a temporary sparse matrix where to store the integrals during the calculation
    ALLOCATE(current_env%p_psi0(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(current_env%rxp_psi0(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       DO idir = 1,3
          NULLIFY(current_env%p_psi0(ispin,idir)%matrix)
          CALL fm_pool_create_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
               &                 current_env%p_psi0(ispin,idir)%matrix,error=error)
          NULLIFY(current_env%rxp_psi0(ispin,idir)%matrix)
          CALL fm_pool_create_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
               &                 current_env%rxp_psi0(ispin,idir)%matrix,error=error)
       ENDDO
    ENDDO

    CALL allocate_matrix_set(current_env%op_rmd_ao,3,error=error)
    CALL set_up_op_sm(current_env%op_rmd_ao(1)%matrix,qs_env,symmetry="none",&
         &            name="OP_RMD",error=error)
    CALL set_matrix(current_env%op_rmd_ao(1)%matrix,0.0_dp)
    DO idir = 2,3
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          current_env%op_rmd_ao(idir)%matrix,&
            &                          "current_env%op_rmd_ao"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(idir))),&
            &                          target_symmetry="none",error=error)
       CALL set_matrix(current_env%op_rmd_ao(idir)%matrix,0.0_dp)
    ENDDO
    !
    ! Current density matrix in x y and z:
    !
    !jp0_ao
    ALLOCATE(current_env%jp0_ao(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       NULLIFY(current_env%jp0_ao(ispin)%matrix)
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          current_env%jp0_ao(ispin)%matrix,"current_env%jp0_ao"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       CALL set_matrix(current_env%jp0_ao(ispin)%matrix,0.0_dp)
    ENDDO
    !
    !jp1_ao
    ALLOCATE(current_env%jp1_ao(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       NULLIFY(current_env%jp1_ao(ispin)%matrix)
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          current_env%jp1_ao(ispin)%matrix,"current_env%jp1_ao"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       CALL set_matrix(current_env%jp1_ao(ispin)%matrix,0.0_dp)
    ENDDO
    !
    !jp2_ao
    ALLOCATE(current_env%jp2_ao(2,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,2
       DO ispin = 1,nspins
          NULLIFY(current_env%jp2_ao(idir,ispin)%matrix)
          CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
               &                          current_env%jp2_ao(idir,ispin)%matrix,"current_env%jp2_ao"//&
               &                          "-"//TRIM(ADJUSTL(cp_to_string(idir)))//&
               &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
               &                          target_symmetry="none",error=error)
       ENDDO
    ENDDO
    !
    ! If the current density on the grid needs to be stored
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
    ALLOCATE(current_env%jrho1_set(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i_B = 1,3
       DO idir = 1,3
          NULLIFY(current_env%jrho1_set(idir,i_B)%rho)
          ALLOCATE(current_env%jrho1_set(idir,i_B)%rho,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          current_env%jrho1_set(idir,i_B)%rho%ref_count = 1
          NULLIFY( current_env%jrho1_set(idir,i_B)%rho%rho_r)
          NULLIFY( current_env%jrho1_set(idir,i_B)%rho%rho_g)
          
          ALLOCATE(current_env%jrho1_set(idir,i_B)%rho%rho_r(nspins),stat=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(current_env%jrho1_set(idir,i_B)%rho%rho_g(nspins),stat=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ispin=1,nspins
             CALL pw_pool_create_pw(auxbas_pw_pool,&
                  &                  current_env%jrho1_set(idir,i_B)%rho%rho_r(ispin)%pw,&
                  &                  use_data=REALDATA3D,in_space=REALSPACE,error=error)
             CALL pw_zero(current_env%jrho1_set(idir,i_B)%rho%rho_r(ispin)%pw)
             
             CALL pw_pool_create_pw(auxbas_pw_pool,&
                  &                  current_env%jrho1_set(idir,i_B)%rho%rho_g(ispin)%pw,&
                  &                  use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                  &                  error=error)
             CALL pw_zero(current_env%jrho1_set(idir,i_B)%rho%rho_g(ispin)%pw)
          ENDDO
          CALL qs_rho_retain(current_env%jrho1_set(idir,i_B)%rho,error=error)
       ENDDO
    ENDDO
    !
    ! Initialize local current density if GAPW calculation
    IF(gapw) THEN
       CALL list_3c_jrho_atom(qs_env,error=error)
       CALL init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,&
            &                  nspins,.TRUE.,& ! vw need to rm the T
            &                  error=error)
       CALL set_current_env(current_env=current_env,jrho1_atom_set=jrho1_atom_set,&
            &           error=error)
    ENDIF
    !
    ALLOCATE(current_env%basisfun_center(3,current_env%nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    current_env%basisfun_center = 0.0_dp


    ALLOCATE (first_sgf(natom),last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)
    !Build 3 arrays where for each contracted basis function
    !the x y and z coordinates of the center are given
    DO iatom = 1,natom
       DO iao = first_sgf(iatom),last_sgf(iatom)
          DO idir = 1,3
             current_env%basisfun_center(idir,iao) = particle_set(iatom)%r(idir)
          END DO
       END DO
    END DO

    DEALLOCATE (first_sgf,last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    
    current_env%simple_done(1:6) = .FALSE.

    ALLOCATE(current_env%full_done(3*nspins,max_nbr_center),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    current_env%full_done = .FALSE.

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)
    
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE current_env_init


  SUBROUTINE current_env_cleanup(current_env,error)

    TYPE(current_env_type)                   :: current_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_B, idir, ispin, istat, j
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(.NOT. failure) THEN
       current_env%ref_count = current_env%ref_count - 1
       IF(current_env%ref_count == 0 ) THEN
          !psi1_p
          IF(ASSOCIATED(current_env%psi1_p)) THEN
             DO idir = 1,SIZE(current_env%psi1_p,2)
                DO ispin = 1,SIZE(current_env%psi1_p,1)
                   CALL fm_pool_give_back_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
                        current_env%psi1_p(ispin,idir)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%psi1_p, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !psi1_rxp
          IF(ASSOCIATED(current_env%psi1_rxp)) THEN
             DO idir = 1,SIZE(current_env%psi1_rxp,2)
                DO ispin = 1,SIZE(current_env%psi1_rxp,1)
                   CALL fm_pool_give_back_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
                        current_env%psi1_rxp(ispin,idir)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%psi1_rxp, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !psi1_D
          IF(ASSOCIATED(current_env%psi1_D)) THEN
             DO idir = 1,SIZE(current_env%psi1_D,2)
                DO ispin = 1,SIZE(current_env%psi1_D,1)
                   CALL fm_pool_give_back_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
                        current_env%psi1_D(ispin,idir)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%psi1_D, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !p_psi0
          IF(ASSOCIATED(current_env%p_psi0)) THEN
             DO idir = 1,SIZE(current_env%p_psi0,2)
                DO ispin = 1,SIZE(current_env%p_psi0,1)
                   CALL fm_pool_give_back_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
                        current_env%p_psi0(ispin,idir)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%p_psi0, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !rxp_psi0
          IF(ASSOCIATED(current_env%rxp_psi0)) THEN
             DO idir = 1,SIZE(current_env%rxp_psi0,2)
                DO ispin = 1,SIZE(current_env%rxp_psi0,1)
                   CALL fm_pool_give_back_fm(current_env%ao_mo_fm_pools(ispin)%pool,&
                        current_env%rxp_psi0(ispin,idir)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%rxp_psi0, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          
          NULLIFY(current_env%ao_mo_fm_pools)
          
          DO ispin = 1,SIZE(current_env%centers_set,1)
             DEALLOCATE(current_env%centers_set(ispin)%array,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(current_env%centers_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          DO ispin = 1,SIZE(current_env%center_list,1)
             DEALLOCATE(current_env%center_list(ispin)%array,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(current_env%center_list,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          
          IF(ASSOCIATED(current_env%list_cubes)) THEN
             DEALLOCATE(current_env%list_cubes)
          END IF
          !op_p_ao
          IF(ASSOCIATED(current_env%op_p_ao)) THEN
             CALL deallocate_matrix_set(current_env%op_p_ao,error=error)
          END IF
          !op_rmd_ao
          IF(ASSOCIATED(current_env%op_rmd_ao)) THEN
             CALL deallocate_matrix_set(current_env%op_rmd_ao,error=error)
          END IF
          !jp0_ao
          IF(ASSOCIATED(current_env%jp0_ao)) THEN
             DO i = 1,SIZE(current_env%jp0_ao,1)
                CALL deallocate_matrix(current_env%jp0_ao(i)%matrix,error=error)
             END DO
             DEALLOCATE(current_env%jp0_ao,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !jp1_ao
          IF(ASSOCIATED(current_env%jp1_ao)) THEN
             DO i = 1,SIZE(current_env%jp1_ao,1)
                CALL deallocate_matrix(current_env%jp1_ao(i)%matrix,error=error)
             END DO
             DEALLOCATE(current_env%jp1_ao,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !jp2_ao
          IF(ASSOCIATED(current_env%jp2_ao)) THEN
             DO i = 1,SIZE(current_env%jp2_ao,1)
                DO j = 1,SIZE(current_env%jp2_ao,2)
                   CALL deallocate_matrix(current_env%jp2_ao(i,j)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%jp2_ao,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          ! Current density on the grid
          IF(ASSOCIATED(current_env%jrho1_set)) THEN
             DO i_B = 1,3
                DO idir = 1,3
                   DO i=1,SIZE(current_env%jrho1_set(idir,i_B)%rho%rho_r)
                      CALL pw_release(current_env%jrho1_set(idir,i_B)%rho%rho_r(i)%pw,error=error)
                   END DO
                   DO i=1,SIZE(current_env%jrho1_set(idir,i_B)%rho%rho_g)
                      CALL pw_release(current_env%jrho1_set(idir,i_B)%rho%rho_g(i)%pw,error=error)
                   END DO
                   DEALLOCATE(current_env%jrho1_set(idir,i_B)%rho%rho_r,stat=istat)
                   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                   DEALLOCATE(current_env%jrho1_set(idir,i_B)%rho%rho_g,stat=istat)
                   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                   DEALLOCATE(current_env%jrho1_set(idir,i_B)%rho,STAT=istat)
                   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                   !
                   ! here ref_count=2 and qs_rho_release doesnt relase the structure
                   !CALL qs_rho_release(current_env%jrho1_set(idir,i_B)%rho,error=error)
                END DO
             END DO
             DEALLOCATE(current_env%jrho1_set,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !IF(ASSOCIATED(current_env%jrho1_set)) THEN
          !   DO i_B = 1,3
          !      DO idir = 1,3
          !         CALL qs_rho_release(current_env%jrho1_set(idir,i_B)%rho,error=error)
          !      END DO
          !   END DO
          !END IF
          ! Local current density, atom by atom (only gapw)
          IF(ASSOCIATED(current_env%jrho1_atom_set)) THEN
             CALL deallocate_jrho_atom_set(current_env%jrho1_atom_set,error=error)
          END IF
          
          !fullnmr_done
          IF(ASSOCIATED(current_env%full_done)) THEN
             DEALLOCATE(current_env%full_done,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(current_env%basisfun_center)) THEN
             DEALLOCATE(current_env%basisfun_center,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(current_env%statetrueindex)) THEN
             DEALLOCATE(current_env%statetrueindex,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

       END IF  ! ref count
    END IF ! failure

  END SUBROUTINE current_env_cleanup

END MODULE  qs_linres_current_utils
