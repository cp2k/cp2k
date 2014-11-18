!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, &
       convert_offsets_to_sizes, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_set, dbcsr_distribution_obj, &
       dbcsr_type_antisymmetric
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: current_gauge_atom,&
                                             current_gauge_r,&
                                             current_gauge_r_and_step_func,&
                                             current_orb_center_atom,&
                                             current_orb_center_box,&
                                             current_orb_center_common,&
                                             current_orb_center_wannier,&
                                             ot_precond_full_all
  USE input_section_types,             ONLY: section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_linres_methods,               ONLY: linres_read_restart,&
                                             linres_solver,&
                                             linres_write_restart
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
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             localized_wfn_control_type,&
                                             qs_loc_env_new_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_types,                    ONLY: qs_rho_clear,&
                                             qs_rho_create,&
                                             qs_rho_retain,&
                                             qs_rho_set
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             rs_grid_release
  USE scf_control_types,               ONLY: scf_control_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: current_response, current_env_cleanup, current_env_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_current_utils'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param p_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_response(current_env,p_env,qs_env,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_response', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: my_pos
    INTEGER :: first_center, handle, i, icenter, idir, ii, iii, ispin, &
      ist_true, istat, istate, j, jcenter, jstate, max_nbr_center, &
      max_states, nao, natom, nbr_center(2), ncubes, nmo, nspins, nstates(2), &
      output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: list_cubes, rbs
    INTEGER, DIMENSION(:, :, :), POINTER     :: statetrueindex
    LOGICAL                                  :: append_cube, failure, ionode, &
                                                should_stop
    REAL(dp)                                 :: dk(3), dkl(3), dl(3)
    REAL(dp), DIMENSION(:), POINTER          :: dkl_vec_ii, dkl_vec_iii
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_dklxp0
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op_p_ao
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fm_work_ii, fm_work_iii, &
                                                h1_psi0, hpsi0_ptr, &
                                                psi0_order, psi1, psi1_ptr
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: p_psi0, psi1_D, psi1_p, &
                                                psi1_rxp, rxp_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: current_section, lr_section, &
                                                print_key

    failure = .FALSE.
    CALL timeset(routineN,handle)
    !
    NULLIFY(cell, dft_control, linres_control, lr_section, current_section,&
            logger, mpools, psi1,h1_psi0, vecbuf_dklxp0, mo_coeff,para_env, &
            fm_work_ii, fm_work_iii, tmp_fm_struct, dkl_vec_ii,dkl_vec_iii, &
            list_cubes, statetrueindex, centers_set, center_list, psi1_p, psi1_rxp, psi1_D, &
            p_psi0, rxp_psi0, psi0_order, op_p_ao, sab_orb, particle_set)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    current_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%CURRENT",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)
    IF(output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "*** Self consistent optimization of the response wavefunctions ***"
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control= dft_control,&
                    mpools=mpools,cell=cell,&
                    linres_control=linres_control,&
                    sab_orb=sab_orb,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set,&
                    dbcsr_dist=dbcsr_dist,&
                    para_env=para_env,error=error)

    nspins = dft_control%nspins

    CALL get_current_env(current_env=current_env, &
                         nao=nao, &
                         nstates=nstates, &
                         centers_set=centers_set, &
                         nbr_center=nbr_center, &
                         center_list=center_list, &
                         list_cubes=list_cubes, &
                         statetrueindex=statetrueindex, &
                         psi1_p=psi1_p, &
                         psi1_rxp=psi1_rxp, &
                         psi1_D=psi1_D,&
                         p_psi0=p_psi0, &
                         rxp_psi0=rxp_psi0, &
                         psi0_order=psi0_order,&
                         error=error)
    !
    ! allocate the vectors
    IF(current_env%full) THEN
       ALLOCATE(psi1(nspins),h1_psi0(nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
          mo_coeff => psi0_order(ispin)%matrix
          NULLIFY(tmp_fm_struct,psi1(ispin)%matrix,h1_psi0(ispin)%matrix)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                                   ncol_global=nstates(ispin),&
                                   context=mo_coeff%matrix_struct%context,&
                                   error=error)
          CALL cp_fm_create(psi1(ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(h1_psi0(ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       ENDDO
       !
       ! prepare for allocation
       natom = SIZE(particle_set,1)
       ALLOCATE (first_sgf(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (last_sgf(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL get_particle_set(particle_set, qs_kind_set,&
            first_sgf=first_sgf,&
            last_sgf=last_sgf,error=error)
       ALLOCATE (rbs(natom), STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
       CALL array_nullify (row_blk_sizes)
       CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
       DEALLOCATE (first_sgf,last_sgf,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !
       ! rebuild the linear momentum matrices
       CALL cp_dbcsr_allocate_matrix_set(op_p_ao,3,error=error)
       ALLOCATE(op_p_ao(1)%matrix)
       CALL cp_dbcsr_init(op_p_ao(1)%matrix,error=error)
       CALL cp_dbcsr_create(matrix=op_p_ao(1)%matrix, &
            name="OP_P", &
            dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
            row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
            nblks=0, nze=0, mutable_work=.TRUE., &
            error=error)
       CALL cp_dbcsr_alloc_block_from_nbl(op_p_ao(1)%matrix,sab_orb,error=error)
       !
       !
       CALL array_release (row_blk_sizes)
       !
       !
       CALL cp_dbcsr_set(op_p_ao(1)%matrix,0.0_dp,error=error)
       DO idir = 2,3
          ALLOCATE(op_p_ao(idir)%matrix)
          CALL cp_dbcsr_init(op_p_ao(idir)%matrix,error=error)
          CALL cp_dbcsr_copy(op_p_ao(idir)%matrix,op_p_ao(1)%matrix,&
               "current_env%op_p_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),&
               error=error)
          CALL cp_dbcsr_set(op_p_ao(idir)%matrix,0.0_dp,error=error)
       ENDDO
       !
       !CALL p_xyz_ao(op_p_ao,qs_env,minimum_image=.FALSE.,error=error)
       CALL build_lin_mom_matrix(qs_env,op_p_ao,error)
       !
    ENDIF
    !
    ! set response to zero
    !
    DO idir = 1,3
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1_p(ispin,idir)%matrix  ,0.0_dp,error=error)
          CALL cp_fm_set_all(psi1_rxp(ispin,idir)%matrix,0.0_dp,error=error)
          IF(current_env%full) CALL cp_fm_set_all(psi1_D(ispin,idir)%matrix,0.0_dp,error=error)
       ENDDO
    ENDDO
    !
    !
    !
    ! load restart file
    !
    first_center = 0
    IF(linres_control%linres_restart) THEN
       DO idir = 1,3
          ! operator p
          psi1_ptr => psi1_p(:,idir)
          CALL linres_read_restart(qs_env,lr_section,psi1_ptr,idir,"nmr_p",error=error)
          ! operator rxp
          psi1_ptr => psi1_rxp(:,idir)
          CALL linres_read_restart(qs_env,lr_section,psi1_ptr,idir,"nmr_rxp",error=error)
       ENDDO
    ENDIF
    !
    !
    !
    !
    should_stop =.FALSE.
    current_env%all_pert_op_done =.FALSE.
    ! operator p
    DO idir = 1,3
       IF(should_stop) EXIT
       IF(output_unit>0) THEN
          WRITE(output_unit,"(T10,A)") "Response to the perturbation operator P_"//ACHAR(idir+119)
       ENDIF
       !
       ! Initial guess for psi1
       psi1_ptr => psi1_p(:,idir)
       hpsi0_ptr => p_psi0(:,idir)
       !
       !
       linres_control%converged = .FALSE.
       CALL linres_solver(p_env,qs_env,psi1_ptr,hpsi0_ptr,psi0_order,should_stop,error=error)
       !
       !
       ! print response functions
       IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
            &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          ncubes = SIZE(list_cubes,1)
          print_key => section_vals_get_subs_vals(current_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
          append_cube = section_get_lval(current_section,"PRINT%RESPONSE_FUNCTION_CUBES%APPEND",error=error)
          my_pos = "REWIND"
          IF(append_cube) THEN
            my_pos = "APPEND"
          END IF
          !

          DO ispin = 1,nspins
             CALL qs_print_cubes(qs_env,psi1_ptr(ispin)%matrix,ncubes,list_cubes,&
                  centers_set(ispin)%array,print_key,'psi1_p',&
                  idir=idir,ispin=ispin,file_position=my_pos,error=error)
          ENDDO ! ispin
       ENDIF ! print response functions
       !
       ! write restart file
       CALL linres_write_restart(qs_env,lr_section,psi1_ptr,idir,"nmr_p",error=error)
    ENDDO ! idir
    !
    ! operator rxp
    DO idir = 1,3
       IF(should_stop) EXIT
       IF(output_unit>0) THEN
          WRITE(output_unit,"(T10,A)") "Response to the perturbation operator L_"//ACHAR(idir+119)
       ENDIF
       !
       ! Initial guess for psi1
       psi1_ptr => psi1_rxp(:,idir)
       hpsi0_ptr => rxp_psi0(:,idir)
       !
       !
       linres_control%converged = .FALSE.
       CALL linres_solver(p_env,qs_env,psi1_ptr,hpsi0_ptr,psi0_order,should_stop,error=error)
       !
       ! print response functions
       IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
            &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          ncubes = SIZE(list_cubes,1)
          print_key => section_vals_get_subs_vals(current_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)

          DO ispin = 1,nspins
             CALL qs_print_cubes(qs_env,psi1_ptr(ispin)%matrix,ncubes,list_cubes,&
                  centers_set(ispin)%array,print_key,'psi1_rxp',&
                  idir=idir,ispin=ispin,file_position=my_pos, error=error)
          ENDDO ! ispin
       ENDIF ! print response functions
       !
       ! write restart file
       CALL linres_write_restart(qs_env,lr_section,psi1_ptr,idir,"nmr_rxp",error=error)
    ENDDO ! idir
    IF(.NOT. should_stop) current_env%all_pert_op_done =.TRUE.
    !
    ! operator D
    IF(current_env%full) THEN
       current_env%all_pert_op_done =.FALSE.
       !
       !
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
       ENDDO
       !
       ! The correction is state depedent a loop over the states is necessary
       max_nbr_center = MAXVAL(nbr_center(1:nspins))
       max_states = MAXVAL(nstates(1:nspins))
       !
       ALLOCATE(vecbuf_dklxp0(1,nao),fm_work_ii(nspins),fm_work_iii(nspins),&
                dkl_vec_ii(max_states),dkl_vec_iii(max_states),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       vecbuf_dklxp0(1,nao) = 0.0_dp
       !
       DO ispin = 1,nspins
          nmo = nstates(ispin)
          mo_coeff => psi0_order(ispin)%matrix
          NULLIFY(tmp_fm_struct,fm_work_ii(ispin)%matrix,fm_work_iii(ispin)%matrix)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                                   ncol_global=nmo,para_env=para_env,&
                                   context=mo_coeff%matrix_struct%context,&
                                   error=error)

          CALL cp_fm_create(fm_work_ii(ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_set_all(fm_work_ii(ispin)%matrix,0.0_dp,error=error)
          CALL cp_fm_create(fm_work_iii(ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_set_all(fm_work_iii(ispin)%matrix,0.0_dp,error=error)
          CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       ENDDO ! ispin
       !
       DO idir = 1,3
          IF(should_stop) EXIT
          DO ispin = 1,nspins
            CALL cp_fm_set_all(psi1_D(ispin,idir)%matrix,0.0_dp,error=error)
          END DO
          first_center = 0
          IF(linres_control%linres_restart) THEN
            psi1_ptr => psi1_D(:,idir)
            CALL linres_read_restart(qs_env,lr_section,psi1_ptr,idir,"nmr_dxp-",ind=first_center,error=error)
          END IF 
          IF(first_center>0) THEN
            IF(output_unit>0) THEN
                  WRITE(output_unit,"(T10,A,I4,A)")&
                     & "Response to the perturbation operators (dk-dl)xp up to state ",&
                       first_center," have been read from restart"
            ENDIF
          END IF
          !
          ! here we run over the max number of states, then we need
          ! to be careful with overflow while doing uks calculations.
          DO icenter = 1,max_nbr_center
             !
            IF(should_stop) EXIT
            IF(icenter > first_center) THEN
             IF(output_unit>0) THEN
                WRITE(output_unit,"(T10,A,I4,A)")&
                     & "Response to the perturbation operator (dk-dl)xp for -state- ",&
                       icenter," in dir. "//ACHAR(idir+119)
             ENDIF
             !
             DO ispin = 1,nspins
                nmo = nstates(ispin)
                mo_coeff => psi0_order(ispin)%matrix
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
                CALL cp_fm_to_fm(mo_coeff,fm_work_ii(ispin)%matrix,error=error)
                CALL cp_fm_column_scale(fm_work_ii(ispin)%matrix,dkl_vec_ii(1:nmo))
                !
                ! Apply the p_iii operator
                ! fm_work_iii = -p_iii * (dk-dl)_ii * C0
                CALL cp_dbcsr_sm_fm_multiply(op_p_ao(iii)%matrix,fm_work_ii(ispin)%matrix,&
                                          fm_work_iii(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,&
                                          error=error)
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
                CALL cp_dbcsr_sm_fm_multiply(op_p_ao(ii)%matrix,fm_work_iii(ispin)%matrix,&
                                          fm_work_ii(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,&
                                          error=error)
                !
                ! Copy in h1_psi0
                ! h1_psi0_i = fm_work_iii - fm_work_ii
                CALL cp_fm_scale_and_add(1.0_dp,h1_psi0(ispin)%matrix,&
                     &                  -1.0_dp,fm_work_ii(ispin)%matrix,&
                                         error=error)
             ENDDO

             !
             ! Optimize the response wavefunctions
             CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,should_stop,error=error)
             !
             IF(output_unit>0) THEN
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
                ist_true = statetrueindex(idir,icenter,ispin)
                DO i = center_list(ispin)%array(1,ist_true),center_list(ispin)%array(1,ist_true+1)-1
                   istate = center_list(ispin)%array(2,i)
                   !
                   ! the optimized wfns are copied in the fm
                   CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_D(ispin,idir)%matrix,1,istate,istate)
                ENDDO
                current_env%full_done(idir*ispin,icenter) = .TRUE.
             ENDDO ! ispin
             !
           ELSE
             DO ispin = 1,nspins
                current_env%full_done(idir*ispin,icenter) = .TRUE.
             ENDDO ! ispin
           END IF
           psi1_ptr => psi1_D(:,idir)
           CALL  linres_write_restart(qs_env,lr_section,psi1_ptr,idir,"nmr_dxp-",ind=icenter,error=error)

          ENDDO ! center
          !
          ! print response functions
          IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
               &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
             ncubes = SIZE(list_cubes,1)
             print_key => section_vals_get_subs_vals(current_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
             DO ispin = 1,nspins
                CALL qs_print_cubes(qs_env,psi1_D(ispin,idir)%matrix,&
                     ncubes,list_cubes,centers_set(ispin)%array,print_key,'psi1_D',&
                      idir=idir,ispin=ispin, file_position=my_pos,error=error)
             ENDDO
          ENDIF ! print response functions
          !
       ENDDO ! idir
       IF(.NOT. should_stop) current_env%all_pert_op_done =.TRUE.
       !
       ! clean up
       DO ispin = 1,nspins
          CALL cp_fm_release(fm_work_ii(ispin)%matrix,error=error)
          CALL cp_fm_release(fm_work_iii(ispin)%matrix,error=error)
       ENDDO
       DEALLOCATE(fm_work_ii,fm_work_iii,dkl_vec_ii,dkl_vec_iii,vecbuf_dklxp0,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ENDIF
    !
    ! clean up
    IF(current_env%full) THEN
       CALL cp_dbcsr_deallocate_matrix_set(op_p_ao,error=error)
       DO ispin = 1,nspins
          CALL cp_fm_release(psi1(ispin)%matrix,error=error)
          CALL cp_fm_release(h1_psi0(ispin)%matrix,error=error)
       ENDDO
       DEALLOCATE(psi1,h1_psi0,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE current_response

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_env_init(current_env,qs_env,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i, iao, iatom, ibox, icenter, icount, idir, ii, &
      ini, ir, is, ispin, istat, istate, istate2, istate_next, ix, iy, iz, j, &
      jstate, k, max_nbr_center, max_states, n, n_rep, nao, natom, nbr_box, &
      ncubes, nmo, nspins, nstate, nstate_list(2), output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buff, first_sgf, last_sgf
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: state_list
    INTEGER, DIMENSION(:), POINTER           :: bounds, list, nbox, &
                                                selected_states_on_atom_list
    LOGICAL                                  :: failure, gapw, is0, &
                                                uniform_occupation
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: state_done
    REAL(dp)                                 :: center(3), center2(3), dist, &
                                                mdist, r(3), rab(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rbuff
    REAL(dp), DIMENSION(:), POINTER          :: common_center
    REAL(dp), DIMENSION(:, :), POINTER       :: center_array
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: current_section, lr_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, qs_kind_set, cell, dft_control, linres_control, scf_control,&
            logger, mos, mpools, current_section, particle_set, mo_coeff,&
            auxbas_pw_pool, pw_env, jrho1_atom_set, common_center, tmp_fm_struct,&
            para_env, qs_loc_env,localized_wfn_control, rho_g, rho_r)
    !
    !
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    mos=mos,&
                    mpools=mpools,&
                    particle_set=particle_set,&
                    pw_env=pw_env,&
                    scf_control=scf_control, &
                    para_env=para_env,&
                    error=error)
    !
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    CALL get_mo_set(mo_set=mos(1)%mo_set,nao=nao)
    !
    max_states = 0
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo)
       max_states = MAX(max_states,nmo)
    ENDDO
    !
    !
    !
    ! some checks
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo,&
                       uniform_occupation=uniform_occupation)
       !
       ! check that homo=nmo
       IF (nmo.NE.homo) CALL stop_program(routineN,moduleN,__LINE__,"nmo != homo")
       !
       ! check that the nbr of localized states is equal to the nbr of states
!       nmoloc = SIZE(linres_control%localized_wfn_control%centers_set(ispin)%array,2)
!       IF (nmoloc.NE.nmo) CALL stop_program(routineN,moduleN,__LINE__,&
!                                            "The nbr of localized functions "//&
!            &                               "is not equal to the nbr of states.")
       !
       ! check that the occupation is uniform
       IF (.NOT.uniform_occupation) CALL stop_program(routineN,moduleN,__LINE__,"nmo != homo")
    ENDDO
    !
    !
    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".linresLog",error=error)

    IF(output_unit>0) THEN
       WRITE(output_unit,"(/,T20,A,/)") "*** Start current Calculation ***"
       WRITE(output_unit,"(T10,A,/)") "Inizialization of the current environment"
    ENDIF

    IF(current_env%ref_count /= 0) CALL current_env_cleanup(current_env,qs_env,error=error)

    CALL current_env_create(current_env,error=error)
    !
    !
    current_env%chi_tensor(:,:,:)     = 0.0_dp
    current_env%chi_tensor_loc(:,:,:) = 0.0_dp
    current_env%nao                   = nao
    current_env%full                  = .TRUE.
    current_env%do_selected_states    = .FALSE.
    !
    ! If current_density or full_nmr different allocations are required
    current_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES%CURRENT",error=error)
    !
    ! Select the gauge
    CALL section_vals_val_get(current_section,"GAUGE",i_val=current_env%gauge,error=error)
    SELECT CASE(current_env%gauge)
    CASE(current_gauge_r)
       current_env%gauge_name="R"
    CASE(current_gauge_r_and_step_func)
       current_env%gauge_name="R_AND_STEP_FUNCTION"
    CASE(current_gauge_atom)
       current_env%gauge_name="ATOM"
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown gauge, try again...")
    END SELECT
    !
    ! maximal radius to build the atom gauge
    CALL section_vals_val_get(current_section,"GAUGE_ATOM_RADIUS",&
         r_val=current_env%gauge_atom_radius,error=error)
    ! use old gauge atom
    CALL section_vals_val_get(current_section,"USE_OLD_GAUGE_ATOM",&
         l_val=current_env%use_old_gauge_atom,error=error)
    ! chi for pbc
    CALL section_vals_val_get(current_section,"CHI_PBC",l_val=current_env%chi_pbc,error=error)
    !
    ! use old gauge atom
    CALL section_vals_val_get(current_section,"USE_OLD_GAUGE_ATOM",&
         l_val=current_env%use_old_gauge_atom,error=error)
    !
    ! chi for pbc
    CALL section_vals_val_get(current_section,"CHI_PBC",l_val=current_env%chi_pbc,error=error)
    !
    ! which center for the orbitals shall we use
    CALL section_vals_val_get(current_section,"ORBITAL_CENTER",i_val=current_env%orb_center,error=error)
    SELECT CASE(current_env%orb_center)
    CASE(current_orb_center_wannier)
       !
       current_env%orb_center_name = "WANNIER"
    CASE(current_orb_center_common)
       !
       current_env%orb_center_name = "COMMON"
       current_env%full = .FALSE.
       !
       ! Is there a user specified common_center?
       CALL section_vals_val_get(current_section,"COMMON_CENTER",r_vals=common_center,error=error)
    CASE(current_orb_center_atom)
       !
       current_env%orb_center_name = "ATOM"
    CASE(current_orb_center_box)
       !
       current_env%orb_center_name = "BOX"
       !
       ! Is there a user specified nbox?
       CALL section_vals_val_get(current_section,"NBOX",i_vals=nbox,error=error)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown orbital center, try again...")
    END SELECT
    !
    ! Check if restat also psi0 should be restarted
    !IF(current_env%restart_current .AND. scf_control%density_guess/=restart_guess)THEN
    !   CALL stop_program(routineN,moduleN,__LINE__,"restart_nmr requires density_guess=restart")
    !ENDIF
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure) THEN
       IF(output_unit>0) THEN
          WRITE(output_unit,'(A)') &
               ' To get CURRENT parameters within PBC you need localized zero order orbitals '
       ENDIF
    ENDIF
    qs_loc_env => linres_control%qs_loc_env
    CALL get_qs_loc_env(qs_loc_env,localized_wfn_control=localized_wfn_control,error=error)
    !
    !
    ALLOCATE(current_env%centers_set(nspins),current_env%center_list(nspins),&
             state_list(max_states,nspins),STAT=istat)
    CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
    state_list(:,:) = HUGE(0)
    nstate_list(:)  = HUGE(0)
    !
    !
    !
    ! if qmmm is selected change the definition of the center of the box 0 -> L/2
    IF(current_env%do_qmmm) THEN
       DO ispin = 1,nspins
          DO istate = 1,SIZE(localized_wfn_control%centers_set(ispin)%array,2)
             ! just to be sure...
             r(:) = pbc(localized_wfn_control%centers_set(ispin)%array(:,istate),cell)
             IF(r(1).LT.0.0_dp) THEN
                localized_wfn_control%centers_set(ispin)%array(1,istate) = &
                     r(1) + cell%hmat(1,1)
             ENDIF
             IF(r(2).LT.0.0_dp) THEN
                localized_wfn_control%centers_set(ispin)%array(2,istate) = &
                     r(2) + cell%hmat(2,2)
             ENDIF
             IF(r(3).LT.0.0_dp) THEN
                localized_wfn_control%centers_set(ispin)%array(3,istate) = &
                     r(3) + cell%hmat(3,3)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    !
    !
    ! if the user has requested to compute the response for a subset of the states
    ! we collect them here. it requies the states to be localized!
    CALL section_vals_val_get(current_section,"SELECTED_STATES_ATOM_RADIUS",&
                              r_val=current_env%selected_states_atom_radius,error=error)
    CALL section_vals_val_get(current_section,"SELECTED_STATES_ON_ATOM_LIST",n_rep_val=n_rep,error=error)
    !
    current_env%do_selected_states = n_rep.GT.0
    !
    ! for the moment selected states doesnt work with the preconditioner FULL_ALL
    IF(linres_control%preconditioner_type.EQ.ot_precond_full_all.AND.current_env%do_selected_states) &
         CALL stop_program(routineN,moduleN,__LINE__,&
                           "Selected states doesnt work with the preconditioner FULL_ALL")
    !
    !
    NULLIFY(current_env%selected_states_on_atom_list)
    n = 0
    DO ir = 1,n_rep
       NULLIFY(list)
       CALL section_vals_val_get(current_section,"SELECTED_STATES_ON_ATOM_LIST",&
                                 i_rep_val=ir,i_vals=list,error=error)
       IF(ASSOCIATED(list)) THEN
          CALL reallocate(current_env%selected_states_on_atom_list,1,n+SIZE(list))
          DO ini = 1,SIZE(list)
             current_env%selected_states_on_atom_list(ini+n) = list(ini)
          ENDDO
          n = n + SIZE(list)
       ENDIF
    ENDDO
    !
    ! build the subset
    IF(current_env%do_selected_states) THEN
       selected_states_on_atom_list => current_env%selected_states_on_atom_list
       DO ispin = 1,nspins
          center_array => localized_wfn_control%centers_set(ispin)%array
          nstate = 0
          DO istate = 1,SIZE(center_array,2)
             DO i = 1,SIZE(selected_states_on_atom_list,1)
                iatom = selected_states_on_atom_list(i)
                r(:)  = pbc(center_array(1:3,istate)-particle_set(iatom)%r(:),cell)
                IF(SQRT(DOT_PRODUCT(r,r)).LE.current_env%selected_states_atom_radius) THEN
                   !
                   ! add the state to the list
                   nstate = nstate + 1
                   state_list(nstate,ispin) = istate
                   EXIT
                ENDIF
             ENDDO
          ENDDO
          nstate_list(ispin) = nstate
       ENDDO
    ELSE
       DO ispin = 1,nspins
          center_array => localized_wfn_control%centers_set(ispin)%array
          nstate = 0
          DO istate = 1,SIZE(center_array,2)
             nstate = nstate + 1
             state_list(nstate,ispin) = istate
          ENDDO
          nstate_list(ispin) = nstate
       ENDDO
    ENDIF
    !
    !
    !
    ! clustering the states
    DO ispin = 1,nspins
       nstate = nstate_list(ispin)
       current_env%nstates(ispin) = nstate
       !
       ALLOCATE(current_env%center_list(ispin)%array(2,nstate+1),&
                current_env%centers_set(ispin)%array(3,nstate),STAT=istat)
       CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       current_env%center_list(ispin)%array(:,:) = HUGE(0)
       current_env%centers_set(ispin)%array(:,:) = HUGE(0.0_dp)
       !
       center_array => localized_wfn_control%centers_set(ispin)%array
       !
       ! point to the psi0 centers
       SELECT CASE(current_env%orb_center)
       CASE(current_orb_center_wannier)
          !
          ! use the wannier center as -center-
          current_env%nbr_center(ispin) = nstate
          DO is = 1,nstate
             istate = state_list(is,ispin)
             current_env%centers_set(ispin)%array(1:3,is) = center_array(1:3,istate)
             current_env%center_list(ispin)%array(1,is) = is
             current_env%center_list(ispin)%array(2,is) = istate
          ENDDO
          current_env%center_list(ispin)%array(1,nstate+1) = nstate+1
          !
       CASE(current_orb_center_common)
          !
          ! use a common -center-
          current_env%centers_set(ispin)%array(:,1) = common_center(:)
          current_env%nbr_center(ispin) = 1
          current_env%center_list(ispin)%array(1,1) = 1
          current_env%center_list(ispin)%array(1,2) = nstate+1
          DO is = 1,nstate
             istate = state_list(is,ispin)
             current_env%center_list(ispin)%array(2,is) = istate
          ENDDO
          !
       CASE(current_orb_center_atom)
          !
          ! use the atom as -center-
          ALLOCATE(buff(nstate_list(ispin)),STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
          buff(:) = 0
          !
          DO is = 1,nstate
             istate = state_list(is,ispin)
             mdist = HUGE(0.0_dp)
             DO iatom = 1,natom
                r   = pbc(particle_set(iatom)%r(:),cell)
                rab = pbc(r,center_array(1:3,istate),cell)
                dist = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)
                IF(dist.LT.mdist) THEN
                   buff(is) = iatom
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
             DO is = 1,nstate
                istate = state_list(is,ispin)
                IF(buff(is).EQ.iatom) THEN
                   j = j + 1
                   i = i + 1
                   is0 = .FALSE.
                   current_env%center_list(ispin)%array(2,i) = istate
                ENDIF
             ENDDO
             IF(.NOT.is0) THEN
                IF(output_unit>0) THEN
                   WRITE(output_unit,'(T2,A,I6,A,I6)') 'clustering ',j,' center(s) on atom ',iatom
                ENDIF
                current_env%center_list(ispin)%array(1,ii+1) = &
                       current_env%center_list(ispin)%array(1,ii) + j
                current_env%centers_set(ispin)%array(:,ii) = &
                        pbc(particle_set(iatom)%r,cell)
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
          ALLOCATE(rbuff(3,nbr_box),buff(nstate),STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
          rbuff(:,:) = HUGE(0.0_dp)
          buff(:) = 0
          !
          ibox = 1
          DO iz = 1,nbox(3)
             DO iy = 1,nbox(2)
                DO ix = 1,nbox(1)
                   rbuff(1,ibox) = cell%hmat(1,1)*((REAL(ix,dp)-0.5_dp)/REAL(nbox(1),dp)-0.5_dp)
                   rbuff(2,ibox) = cell%hmat(2,2)*((REAL(iy,dp)-0.5_dp)/REAL(nbox(2),dp)-0.5_dp)
                   rbuff(3,ibox) = cell%hmat(3,3)*((REAL(iz,dp)-0.5_dp)/REAL(nbox(3),dp)-0.5_dp)
                   ibox = ibox+1
                ENDDO
             ENDDO
          ENDDO
          !
          DO is = 1,nstate
             istate = state_list(is,ispin)
             mdist = HUGE(0.0_dp)
             DO ibox = 1,nbr_box
                rab(:) = pbc(rbuff(:,ibox),center_array(1:3,istate),cell)
                dist = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)
                IF(dist.LT.mdist) THEN
                   buff(is) = ibox
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
             DO is = 1,nstate
                istate = state_list(is,ispin)
                IF(buff(is).EQ.ibox) THEN
                   j = j + 1
                   i = i + 1
                   is0 = .FALSE.
                   current_env%center_list(ispin)%array(2,i) = istate
                ENDIF
             ENDDO
             IF(.NOT.is0) THEN
                IF(output_unit>0) THEN
                   WRITE(output_unit,'(T2,A,I6,A,I6)') 'clustering ',j,' center(s) on box ',ibox
                ENDIF
                current_env%center_list(ispin)%array(1,ii+1) = &
                    current_env%center_list(ispin)%array(1,ii) + j
                current_env%centers_set(ispin)%array(:,ii) = rbuff(:,ibox)
                ii = ii + 1
             ENDIF
          ENDDO
          current_env%nbr_center(ispin) = ii-1
          !
          DEALLOCATE(buff,rbuff,STAT=istat)
          CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Unknown orbital center, try again...")
       END SELECT
    ENDDO
    !
    !
    !
    ! block the states for each center
    ALLOCATE(current_env%psi0_order(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       NULLIFY(tmp_fm_struct)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                                ncol_global=nstate_list(ispin),para_env=para_env,&
                                context=mo_coeff%matrix_struct%context,&
                                error=error)
       CALL cp_fm_create(current_env%psi0_order(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       CALL cp_fm_set_all(current_env%psi0_order(ispin)%matrix,0.0_dp,error=error)
    ENDDO
    !
    !
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       jstate = 0
       DO icenter = 1,current_env%nbr_center(ispin)
          DO i = current_env%center_list(ispin)%array(1,icenter),&
               current_env%center_list(ispin)%array(1,icenter+1)-1
             !
             istate = current_env%center_list(ispin)%array(2,i)
             jstate = jstate + 1
             !
             IF(current_env%do_selected_states) THEN ! this should be removed. always reorder the states
                ! the blocking works (so far) with all the precond except FULL_ALL
                ! the center_list(ispin)%array(2,i) array is obsolete then
                current_env%center_list(ispin)%array(2,i) = jstate
                CALL cp_fm_to_fm(mo_coeff,current_env%psi0_order(ispin)%matrix,1,istate,jstate)
             ELSE
                ! for the moment we just copy them
                CALL cp_fm_to_fm(mo_coeff,current_env%psi0_order(ispin)%matrix,1,istate,istate)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
    !
    !
    !
    IF(current_env%do_qmmm.AND.&
         & (current_env%orb_center.EQ.current_orb_center_wannier.OR.&
            current_env%orb_center.EQ.current_orb_center_atom.OR.&
            current_env%orb_center.EQ.current_orb_center_box)) THEN
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
       DO ispin = 1,nspins
          IF(output_unit>0) THEN
             WRITE(output_unit,'(/,T2,A,I2)') "WANNIER CENTERS for spin ",ispin
             WRITE(output_unit,'(/,T18,A)') "--------------- Shifted Centers --------------- "
             DO istate = 1,current_env%nbr_center(ispin)
                WRITE(output_unit,'(T5,A,I6,3F12.6)')&
                     'state ',istate,current_env%centers_set(ispin)%array(1:3,istate)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    !
    !
    !
    ! reorder the center dependent responses
    max_nbr_center = MAXVAL(current_env%nbr_center(1:nspins))
    current_env%nao = nao
    ALLOCATE(current_env%statetrueindex(3,max_nbr_center,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    current_env%statetrueindex(:,:,:) = 0
    IF(.TRUE.) THEN
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
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo)
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
                      IF(dist.LT.mdist) THEN
                         mdist = dist
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
    ELSE
       DO ispin = 1,nspins
          DO icenter = 1,current_env%nbr_center(ispin)
             current_env%statetrueindex(:,icenter,ispin) = icenter
          ENDDO
       ENDDO
    ENDIF
    !
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,"(T2,A,T60,A)" ) "CURRENT| Gauge used",&
            REPEAT(' ',20-LEN_TRIM(current_env%gauge_name))//TRIM(current_env%gauge_name)
       WRITE(output_unit,"(T2,A,T79,L1)") "CURRENT| Use old gauge code",current_env%use_old_gauge_atom
       WRITE(output_unit,"(T2,A,T79,L1)") "CURRENT| Compute chi for PBC calculation",current_env%chi_pbc
       IF(current_env%gauge.EQ.current_gauge_atom) THEN
          WRITE(output_unit,"(T2,A,T68,E12.6)") "CURRENT| Radius for building the gauge",&
                 current_env%gauge_atom_radius
       ENDIF
       WRITE(output_unit,"(T2,A,T60,A)") "CURRENT| Orbital center used",&
            REPEAT(' ',20-LEN_TRIM(current_env%orb_center_name))//TRIM(current_env%orb_center_name)
       IF(current_env%orb_center.EQ.current_orb_center_common) THEN
          WRITE(output_unit,"(T2,A,T50,3F10.6)") "CURRENT| Common center",common_center(1:3)
       ELSEIF(current_env%orb_center.EQ.current_orb_center_box) THEN
          WRITE(output_unit,"(T2,A,T60,3I5)") "CURRENT| Nbr boxes in each direction",nbox(1:3)
       ENDIF
       !
       DO ispin = 1,nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo)
          WRITE(output_unit,'(T2,A,I6,A,I6,A,I2)') 'CURRENT| Compute', &
                  nstate_list(ispin),' selected response functions out of',&
                  nmo,' for spin ',ispin
          !
          WRITE(output_unit,'(T2,A,I6,A,I2)') 'CURRENT| There is a total of ', &
                 current_env%nbr_center(ispin),' (clustered) center(s) for spin ',ispin
       ENDDO
    ENDIF

    IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
         &    "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN

       NULLIFY(bounds,list)
       ncubes = 0
       CALL section_vals_val_get(current_section,&
            &                    "PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LU_BOUNDS",&
                                 i_vals=bounds,error=error)
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
                                    n_rep_val=n_rep,error=error)
          ncubes = 0
          DO ir = 1,n_rep
             NULLIFY(list)
             CALL section_vals_val_get(current_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
                                       i_rep_val=ir,i_vals=list,error=error)
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

    ! for the chemical shift we need 6 psi1, i.e. 6 optimization procedures
    ! They become 9 if full nmr is calculated, i.e. with the correction term too
    ! All of them are required at the end of the optimization procedure
    ! if the current density and the induced fiels have to be calculated
    ! If instead only the shift is needed, only one psi1 should be enough, providing
    ! that after every optimization the corresponding shift contribution is calculated
    ! prepare the psi1
    !
    ! for the rxp we cannot calculate it a priori because it is in facts (r-dk)xp
    ! where dk is the center of the orbital it is applied to. We would need nstate operators
    ! What we can store is (r-dk)xp|psi0k> for each k, which is a full matrix only
    ! Therefore we prepare here the full matrix p_psi0 and rxp_psi0
    ! We also need a temporary sparse matrix where to store the integrals during the calculation
    ALLOCATE(current_env%p_psi0(nspins,3),current_env%rxp_psi0(nspins,3),&
             current_env%psi1_p(nspins,3),current_env%psi1_rxp(nspins,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(current_env%full) THEN
       ALLOCATE(current_env%psi1_D(nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    DO ispin = 1,nspins
       mo_coeff => current_env%psi0_order(ispin)%matrix
       NULLIFY(tmp_fm_struct)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                                ncol_global=current_env%nstates(ispin),&
                                context=mo_coeff%matrix_struct%context,&
                                error=error)
       DO idir = 1,3
          NULLIFY(current_env%psi1_p(ispin,idir)%matrix,current_env%psi1_rxp(ispin,idir)%matrix)
          NULLIFY(current_env%p_psi0(ispin,idir)%matrix,current_env%rxp_psi0(ispin,idir)%matrix)
          CALL cp_fm_create(current_env%psi1_p(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(current_env%psi1_rxp(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(current_env%p_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(current_env%rxp_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
          IF(current_env%full) THEN
             NULLIFY(current_env%psi1_D(ispin,idir)%matrix)
             CALL cp_fm_create(current_env%psi1_D(ispin,idir)%matrix,tmp_fm_struct,error=error)
          ENDIF
       ENDDO
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    ENDDO
    !
    ! If the current density on the grid needs to be stored
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
    ALLOCATE(current_env%jrho1_set(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
       NULLIFY(rho_r, rho_g)
       ALLOCATE(rho_r(nspins), rho_g(nspins))
       DO ispin=1,nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_r(ispin)%pw,&
                                 use_data=REALDATA3D,in_space=REALSPACE,error=error)
          CALL pw_zero(rho_r(ispin)%pw, error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_g(ispin)%pw,&
                                 use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                 error=error)
          CALL pw_zero(rho_g(ispin)%pw, error=error)
       ENDDO
       NULLIFY(current_env%jrho1_set(idir)%rho)
       CALL qs_rho_create(current_env%jrho1_set(idir)%rho, error)
       CALL qs_rho_set(current_env%jrho1_set(idir)%rho,&
                       rho_r=rho_r, rho_g=rho_g, error=error)
       CALL qs_rho_retain(current_env%jrho1_set(idir)%rho,error=error)
    ENDDO
    !
    ! Initialize local current density if GAPW calculation
    IF(gapw) THEN
       CALL init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,nspins,error=error)
       CALL set_current_env(current_env=current_env,jrho1_atom_set=jrho1_atom_set,&
                            error=error)
    ENDIF
    !
    ALLOCATE(current_env%basisfun_center(3,current_env%nao),&
             first_sgf(natom),last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    current_env%basisfun_center = 0.0_dp
    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf,error=error)
    !Build 3 arrays where for each contracted basis function
    !the x y and z coordinates of the center are given
    DO iatom = 1,natom
       DO iao = first_sgf(iatom),last_sgf(iatom)
          DO idir = 1,3
             current_env%basisfun_center(idir,iao) = particle_set(iatom)%r(idir)
          END DO
       END DO
    END DO

    DEALLOCATE(first_sgf,last_sgf,state_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    current_env%simple_done(1:6) = .FALSE.

    ALLOCATE(current_env%full_done(3*nspins,max_nbr_center),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    current_env%full_done = .FALSE.
    !
    ! allocate pointer for the gauge
    ALLOCATE(current_env%rs_gauge(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(current_env%rs_gauge(1)%rs,current_env%rs_gauge(2)%rs,current_env%rs_gauge(3)%rs)

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,"PRINT%PROGRAM_RUN_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE current_env_init


! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_env_cleanup(current_env,qs_env,error)

    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idir, ispin, istat
    LOGICAL                                  :: failure
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs

    failure=.FALSE.

    NULLIFY(pw_env,rs_descs)

    IF(.NOT.failure) THEN
       current_env%ref_count = current_env%ref_count - 1
       IF(current_env%ref_count == 0 ) THEN
          !
          !psi0_order
          IF(ASSOCIATED(current_env%psi0_order)) THEN
             DO ispin = 1,SIZE(current_env%psi0_order,1)
                CALL cp_fm_release(current_env%psi0_order(ispin)%matrix,error=error)
             ENDDO
             DEALLOCATE(current_env%psi0_order,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          !
          !psi1_p
          IF(ASSOCIATED(current_env%psi1_p)) THEN
             DO idir = 1,SIZE(current_env%psi1_p,2)
                DO ispin = 1,SIZE(current_env%psi1_p,1)
                   CALL cp_fm_release(current_env%psi1_p(ispin,idir)%matrix,error=error)
                ENDDO
             ENDDO
             DEALLOCATE(current_env%psi1_p,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          !
          !psi1_rxp
          IF(ASSOCIATED(current_env%psi1_rxp)) THEN
             DO idir = 1,SIZE(current_env%psi1_rxp,2)
                DO ispin = 1,SIZE(current_env%psi1_rxp,1)
                   CALL cp_fm_release(current_env%psi1_rxp(ispin,idir)%matrix,error=error)
                ENDDO
             ENDDO
             DEALLOCATE(current_env%psi1_rxp,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          !
          !psi1_D
          IF(ASSOCIATED(current_env%psi1_D)) THEN
             DO idir = 1,SIZE(current_env%psi1_D,2)
                DO ispin = 1,SIZE(current_env%psi1_D,1)
                   CALL cp_fm_release(current_env%psi1_D(ispin,idir)%matrix,error=error)
                ENDDO
             ENDDO
             DEALLOCATE(current_env%psi1_D,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          !
          !p_psi0
          IF(ASSOCIATED(current_env%p_psi0)) THEN
             DO idir = 1,SIZE(current_env%p_psi0,2)
                DO ispin = 1,SIZE(current_env%p_psi0,1)
                   CALL cp_fm_release(current_env%p_psi0(ispin,idir)%matrix,error=error)
                ENDDO
             ENDDO
             DEALLOCATE(current_env%p_psi0,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          !
          !rxp_psi0
          IF(ASSOCIATED(current_env%rxp_psi0)) THEN
             DO idir = 1,SIZE(current_env%rxp_psi0,2)
                DO ispin = 1,SIZE(current_env%rxp_psi0,1)
                   CALL cp_fm_release(current_env%rxp_psi0(ispin,idir)%matrix,error=error)
                ENDDO
             ENDDO
             DEALLOCATE(current_env%rxp_psi0,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          !
          !
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
             DEALLOCATE(current_env%list_cubes,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          ! Current density on the grid
          IF(ASSOCIATED(current_env%jrho1_set)) THEN
             DO idir = 1,3
                CALL qs_rho_clear(current_env%jrho1_set(idir)%rho, error)
                DEALLOCATE(current_env%jrho1_set(idir)%rho)
             END DO
             DEALLOCATE(current_env%jrho1_set)
          END IF
          !
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
          IF(ASSOCIATED(current_env%selected_states_on_atom_list)) THEN
             DEALLOCATE(current_env%selected_states_on_atom_list,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          ! give back the gauge
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env,rs_descs=rs_descs,error=error)
          DO idir=1,3
             IF(ASSOCIATED(current_env%rs_gauge(idir)%rs)) THEN
                DO i=1,SIZE(current_env%rs_gauge(idir)%rs)
                  CALL rs_grid_release(current_env%rs_gauge(idir)%rs(i)%rs_grid, error=error)
                END DO
                DEALLOCATE(current_env%rs_gauge(idir)%rs,STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             ENDIF
          ENDDO
          !
          DEALLOCATE(current_env%rs_gauge,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          ! give back the buf
          IF(ASSOCIATED(current_env%rs_buf)) THEN
             DO i=1,SIZE(current_env%rs_buf)
               CALL rs_grid_release(current_env%rs_buf(i)%rs_grid, error=error)
             END DO
             DEALLOCATE(current_env%rs_buf,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDIF

       END IF  ! ref count
    END IF ! failure

  END SUBROUTINE current_env_cleanup

END MODULE qs_linres_current_utils
