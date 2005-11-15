!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/xas_tp_scf   [1.0] *
!!
!!   NAME
!!      xas_methods
!!
!!   FUNCTION
!!      xas_scf for the tp method
!!      It is repeaated for every atom that have to be excited
!!
!!   AUTHOR
!!     MI (05.2005)
!!
!!   MODIFICATION HISTORY
!!     created 05.2005
!!
!!   SOURCE
!******************************************************************************


MODULE xas_tp_scf

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                             cp_fm_get_info,&
                                             cp_fm_maxabsval,&
                                             cp_fm_p_type,&
                                             cp_fm_init_random,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type,cp_fm_to_fm
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: ot_precond_full_kinetic,&
                                             xas_tp_hh, xas_tp_fh,&
                                             xas_dscf, do_loc_none ,&
                                             xas_scf_general,xas_scf_ot
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE qs_ot_scf,                       ONLY: ot_scf_mini
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner,&
                                             preconditioner_type
  USE qs_diis,                         ONLY: qs_diis_b_step
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_loc_methods,                  ONLY: optimize_loc_berry  
  USE qs_loc_types,                    ONLY: get_qs_loc_env,& 
                                             qs_loc_env_new_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,set_mo_set,&
                                             mo_set_p_type,mo_set_type,&
                                             write_mo_set
  USE qs_ot_eigensolver,                ONLY: ot_eigensolver
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf,                          ONLY: eigensolver,&
                                             init_scf_run,&
                                             scf_env_density_mixing
  USE qs_scf_types,                    ONLY: ao_mixing_nr,&
                                             special_diag_method_nr,&
                                             general_diag_method_nr,&
                                             ot_method_nr,&
                                             qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type,pbc
  USE sparse_matrix_types,             ONLY: real_matrix_p_type, &
                                             real_matrix_type,&
                                             set_matrix,&
                                             transfer_matrix
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env, set_xas_env,&
                                             xas_environment_type

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: module_name='xas_tp_scf'  

! *** Public subroutines ***

  PUBLIC :: xas_do_tp_scf 
 
!****************************************************************************

CONTAINS


!!****f* xas_methods/xas_do_tp_scf [1.0] *
!!
!!   NAME
!!     xas_do_tp_scf
!!
!!   FUNCTION
!!     perform an scf loop to calculate the xas spectrum 
!!     given by the excitation of a inner state of a selected atom
!!     by using the transition potential method
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - xas_env: the environment for XAS  calculations 
!!     - scf_env: the scf_env where to perform the scf procedure
!!     - qs_env: the qs_env, the scf_env and xas_env live in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE xas_do_tp_scf (xas_control,xas_env,iatom,scf_env,qs_env,globenv,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(xas_environment_type), POINTER      :: xas_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "xas_tp_scf", &
      routineP = module_name//':'//routineN

    INTEGER                                  :: handle, handle2,  homo, ic_max, ir_max,ispin, &
                                                istat, iter_count, nao, nmo,nvirtual2, &
                                                nvirtual1,output_unit,i,istate,uno_iter
    LOGICAL :: diis_step, do_level_shift, energy_only, failure,& 
               gapw, ionode, should_stop, uniform_occupation, use_cholesky, use_jacobi, id_equal
    REAL(KIND=dp)                            :: a_max,b_max, diis_error, IP_energy, &
                                                orthonormality, uno_eps, &
                                                ra(3), rac(3), rc(3), t1, t2, tot1_h, tot1_s
    REAL(dp), DIMENSION(6)                   :: weights
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers, eigenvalues
    REAL(KIND=dp), DIMENSION(:), POINTER     :: all_evals
    REAL(KIND=dp), DIMENSION(:), POINTER     :: uno_evals
    REAL(dp), DIMENSION(:,:), POINTER        :: centers_wfn
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: vecbuffer,vecbuffer2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                     :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: scf_work1
    TYPE(cp_fm_type), POINTER                :: all_vectors
    TYPE(cp_fm_type), POINTER                :: excvec_coeff
    TYPE(cp_fm_type), POINTER                :: excvec_overlap
    TYPE(cp_fm_type), POINTER                :: ortho, scf_work2,fm_work
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_fm_type), POINTER                :: uno_orbs
    TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER  :: op_fm_set
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: moloc_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(preconditioner_type), POINTER       :: preconditioner
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE( qs_loc_env_new_type ), POINTER     :: qs_loc_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_p, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric
    INTEGER, DIMENSION(:), POINTER           :: col_indices
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      POINTER                                    :: op_sm_set
    INTEGER :: m, n

    NULLIFY(dft_control,matrix_h,matrix_s,matrix_ks,matrix_p,rho,energy,scf_control,logger,&
         ks_env,mos,atomic_kind_set,particle_set,vecbuffer,vecbuffer2, occupation_numbers,&
         scf_work1,ortho,scf_work2,excvec_coeff,excvec_overlap,mo_coeff,fm_work,&
         col_indices, eigenvalues)
    NULLIFY(centers_wfn,qs_loc_env,localized_wfn_control,op_sm_set,op_fm_set,moloc_coeff,cell)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(xas_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         dft_control=dft_control,&
         scf_control=scf_control)

    CALL write_checkpoint_information("entering "//routineN,globenv)

    CALL timeset(routineN,"I"," ",handle)

    ionode = logger%para_env%source==logger%para_env%mepos
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    energy_only = .FALSE.

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
            "XAS_TP_SCF WAVEFUNCTION OPTIMIZATION"
    END IF

!   True the GAPW method is used
    gapw = dft_control%qs_control%gapw

    CALL get_qs_env(qs_env=qs_env,&
         matrix_h=matrix_h,&
         matrix_s=matrix_s,energy=energy,&
         particle_set=particle_set,&
         ks_env=ks_env)

    scf_work1 => scf_env%scf_work1
    scf_work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    CALL get_qs_env(qs_env,matrix_ks=matrix_ks,rho=rho,mos=mos)
    matrix_p => rho%rho_ao
    nao = mos(1)%mo_set%nao
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp
    ALLOCATE(vecbuffer2(1,xas_control%nexc_search),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer2 = 0.0_dp

    ! do some assertions here on these matrices having the same structure, 
    ! as is currently required
    DO ispin=1,SIZE(matrix_s)
      id_equal=(matrix_s(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_p)
      id_equal=(matrix_p(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_h)
      id_equal=(matrix_h(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_ks)
      id_equal=(matrix_ks(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    ! end sparsity check

    ! Some elements from the xas_env
     CALL get_xas_env(xas_env=xas_env,excvec_coeff=excvec_coeff,nvirtual2=nvirtual2,&
          excvec_overlap=excvec_overlap,fm_work=fm_work,centers_wfn=centers_wfn)
     CPPrecondition(ASSOCIATED(excvec_overlap),cp_failure_level,routineP,error,failure)

    iter_count = 0
    diis_step = .FALSE.
    use_jacobi = .FALSE.

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,&
            FMT="(/,T3,A,T9,A,T36,A,T49,A,T68,A,/,T3,A)")&
            "Step","Update method","Time","Convergence","Total energy",&
            REPEAT("-",77)
    END IF
    
    IF (.NOT.logger%print_keys%timecumul) t1 = m_walltime()

!   *** SCF loop ***

    scf_loop: DO
       CALL timeset("scf_iter","I"," ",handle2)

       IF (ionode) CALL m_flush(output_unit)

       iter_count = iter_count + 1

! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date

       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
            error=error,&
            calculate_forces=.FALSE.,just_energy=energy_only)


       scf_env%p_mix = xas_env%p_mix
       scf_env%iter_param = 0.0_dp
       IF (scf_env%mixing_method==ao_mixing_nr) scf_env%iter_param = scf_env%p_mix

       SELECT CASE (xas_env%scf_method)
         CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown scf method method:"//&
                 cp_to_string(xas_env%scf_method),error,failure)

         CASE(xas_scf_general) ! diagonalisation (default)

            DO ispin=1,dft_control%nspins 
             CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_work1(ispin)%matrix)
            ENDDO

            IF (iter_count > 1) THEN
               CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_work1,&
                    scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                    xas_env%eps_diis,s_matrix=matrix_s,error=error)
            END IF

            do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
               ((scf_control%density_guess == "CORE").OR.(iter_count > 1)))

            IF (diis_step) THEN
               scf_env%p_mix = 1.0_dp
               scf_env%iter_param = diis_error
               IF (use_jacobi) THEN
                  scf_env%iter_method = "DIIS/Jacobi"
               ELSE
                  scf_env%iter_method = "DIIS/Diag."
               END IF
            ELSE
               IF (use_jacobi) THEN
                  scf_env%iter_method = "Mixing/Jacobi"
               ELSE
                  scf_env%iter_method = "Mixing/Diag."
               END IF
            END IF

            IF ((iter_count > 1).AND.(scf_env%iter_delta < scf_control%eps_jacobi)) THEN
               use_jacobi = .TRUE.
            ELSE
               use_cholesky = scf_control%use_cholesky
               use_jacobi = .FALSE.
            END IF

            scf_env%iter_delta=0.0_dp

            DO ispin=1,dft_control%nspins

               CALL eigensolver(scf_work1(ispin)%matrix, mos(ispin)%mo_set,ortho,scf_work2,&
                    do_level_shift,scf_control%level_shift,&
                    use_cholesky=use_cholesky,&
                    work_syevx=scf_control%work_syevx,&
                    use_jacobi=use_jacobi,&
                    jacobi_threshold=scf_control%jacobi_threshold,&
                    smear=scf_control%smear, &
                    error=error)

               IF(ispin == 1) THEN
! ** use the maximum overlap criterion to find the index of the excited orbital
                 CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
                      eigenvalues=eigenvalues,homo=homo)
!TEST : Calculate the centers of the wavefunctions to check the localization
                 IF(.FALSE.) THEN
                   CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
                   CALL get_qs_loc_env(qs_loc_env=qs_loc_env,&
                        localized_wfn_control=localized_wfn_control,&
                        moloc_coeff=moloc_coeff,&
                        op_sm_set=op_sm_set,op_fm_set=op_fm_set,&
                        cell=cell, weights=weights)
 
                   CALL cp_fm_to_fm(mo_coeff,moloc_coeff(1)%matrix,xas_control%nexc_search,1,1 )
                   CALL optimize_loc_berry(do_loc_none, localized_wfn_control, &
                        moloc_coeff(1)%matrix, op_sm_set, &
                        op_fm_set, cell, weights, 1, ionode, error=error)

                   ra(1:3) = particle_set(iatom)%r(1:3)
                   DO istate = 1,xas_control%nexc_search
                      centers_wfn(1,istate) = localized_wfn_control%centers_set(1)%array(1,istate)
                      centers_wfn(2,istate) = localized_wfn_control%centers_set(1)%array(2,istate)
                      centers_wfn(3,istate) = localized_wfn_control%centers_set(1)%array(3,istate)
                      rc(1:3) = centers_wfn(1:3,istate)
                      rac = pbc(ra,rc,cell)
                      IF(ionode) WRITE(*,'(I4,4f12.6 )' ) istate, rac, eigenvalues(istate)
                   END DO 
                 END IF
!TEST

                 CALL cp_sm_fm_multiply(matrix_s(1)%matrix,mo_coeff,fm_work,ncol=nmo)
                 CALL cp_fm_gemm("T","N",1,xas_control%nexc_search,nao,1.0_dp,excvec_coeff,&
                      fm_work,0.0_dp,excvec_overlap,b_first_col=1,error=error)
                 CALL cp_fm_get_info ( matrix=excvec_overlap, col_indices = col_indices,&
                        nrow_global = m, ncol_global = n )
!        CALL cp_fm_get_info ( matrix=excvec_coeff, nrow_global = m, ncol_global = n )
!        CALL cp_fm_get_info ( matrix=fm_work, nrow_global = m, ncol_global = n )
                 CALL cp_fm_get_submatrix(excvec_overlap,vecbuffer2,1,1,&
                      1,xas_control%nexc_search,transpose=.FALSE.,error=error)
!dbg
! if(ionode)  WRITE(*,*) globenv%para_env%mepos, vecbuffer2(1,1:SIZE(vecbuffer2,2))
! if(ionode) WRITE(*,*) MAXLOC(ABS(vecbuffer2(1,:))), MAXVAL(ABS(vecbuffer2(1,:)))
!dbg
                 uniform_occupation = .FALSE.
                 CALL get_mo_set(mos(1)%mo_set,&
                        occupation_numbers=occupation_numbers)

                 b_max = 0.0_dp 

                 DO i = 1,xas_control%nexc_search
                   a_max = ABS(vecbuffer2(1,i))
                   IF(a_max > b_max)THEN
                      ic_max = i 
                      b_max = a_max
                   ENDIF
                 END DO

                 IF(ic_max /= dft_control% xas_estate) THEN
                   ir_max = dft_control% xas_estate
                   dft_control% xas_estate = ic_max

                   IF(xas_control%xas_method == xas_tp_hh) THEN
                       occupation_numbers(dft_control%xas_estate) = 0.5_dp
                       occupation_numbers(ir_max) = 1.0_dp
                   ELSE IF(xas_control%xas_method == xas_tp_fh) THEN
                       occupation_numbers(dft_control%xas_estate) = 0.0_dp
                       occupation_numbers(ir_max) = 1.0_dp
                   ELSE IF(xas_control%xas_method == xas_dscf) THEN
                       occupation_numbers(dft_control%xas_estate) = 0.0_dp
                       occupation_numbers(ir_max) = 1.0_dp
                   END IF

                 END IF
                 ! IOnization Potential
                 IP_energy = eigenvalues(dft_control% xas_estate)
!dbg
!         IF(ionode) WRITE(*,'(A, I4, A,f5.2)') ' st ', dft_control% xas_estate,&
!                ' occ ',  occupation_numbers(dft_control% xas_estate)
!dbg
                 CALL set_mo_set(mos(1)%mo_set,&
                      uniform_occupation=uniform_occupation)
                 CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,dft_control%xas_estate,&
                      nao,1,transpose=.TRUE.,error=error)
                 CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
                      nao,1,transpose=.TRUE.,error=error)
               END IF 

               IF (logger%print_keys%each_scf_step) THEN
                  CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
               END IF

               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                    scf_env%p_mix_new(ispin)%matrix,error=error)

            ENDDO

         CASE(xas_scf_ot) ! orbital transforms
! in case of LSD the first spin qs_ot_env will drive the minimization
! in the case of a restricted calculation, it will make sure the spin orbitals are equal
            orthogonality_metric=>matrix_s(1)%matrix

            CALL ot_scf_mini(mos,qs_env%mo_derivs,orthogonality_metric, &
                             energy%total, energy_only,scf_env%iter_delta, &
                             scf_env%qs_ot_env)

            DO ispin=1,dft_control%nspins
              IF(ispin == 1) THEN
! ** use the maximum overlap criterion to find the index of the excited orbital
                 CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
                      eigenvalues=eigenvalues,homo=homo)

                 CALL cp_sm_fm_multiply(matrix_s(1)%matrix,mo_coeff,fm_work,ncol=nmo)

                 CALL cp_fm_gemm("T","N",1,xas_control%nexc_search,nao,1.0_dp,excvec_coeff,&
                      fm_work,0.0_dp,excvec_overlap,b_first_col=1,error=error)
                 CALL cp_fm_get_info ( matrix=excvec_overlap, col_indices = col_indices)

                 CALL cp_fm_get_submatrix(excvec_overlap,vecbuffer2,1,1,&
                      1,xas_control%nexc_search,transpose=.TRUE.,error=error)

                 uniform_occupation = .FALSE.
                 CALL get_mo_set(mos(1)%mo_set,&
                        occupation_numbers=occupation_numbers)

                 b_max = 0.0_dp
 
                 DO i = 1,xas_control%nexc_search
                   a_max = ABS(vecbuffer2(1,i))
                   IF(a_max > b_max)THEN
                      ic_max = i 
                      b_max = a_max
                   ENDIF
                 END DO

                 IF(ic_max /= dft_control% xas_estate) THEN
                   ir_max = dft_control% xas_estate
                   dft_control% xas_estate = ic_max

                   IF(xas_control%xas_method == xas_tp_hh) THEN
                       occupation_numbers(dft_control%xas_estate) = 0.5_dp
                       occupation_numbers(ir_max) = 1.0_dp
                   ELSE IF(xas_control%xas_method == xas_tp_fh) THEN
                       occupation_numbers(dft_control%xas_estate) = 0.0_dp
                       occupation_numbers(ir_max) = 1.0_dp
                   ELSE IF(xas_control%xas_method == xas_dscf) THEN
                       occupation_numbers(dft_control%xas_estate) = 0.0_dp
                       occupation_numbers(ir_max) = 1.0_dp
                   END IF

                 END IF

                 ! IOnization Potential
                 IP_energy = eigenvalues(dft_control% xas_estate)

                 CALL set_mo_set(mos(1)%mo_set,&
                      uniform_occupation=uniform_occupation)
                 CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,dft_control%xas_estate,&
                      nao,1,transpose=.TRUE.,error=error)
                 CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
                      nao,1,transpose=.TRUE.,error=error)
               END IF

               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                     rho%rho_ao(ispin)%matrix,error=error)
            ENDDO

            scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
            scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min

       END SELECT


       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
          CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
                                    scf_env%p_mix,scf_env%iter_delta,&
                                    qs_env=qs_env,error=error)
       ENDIF


       t2 = m_walltime()

       IF (ionode.AND.logger%print_keys%scf.and.scf_env%print_iter_line) THEN
          WRITE (UNIT=output_unit,&
               FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
               iter_count,TRIM(scf_env%iter_method),&
               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
       END IF

! ** convergence check
       CALL external_control(should_stop,"SCF",error)
       IF (scf_env%iter_delta < xas_control%eps_scf) THEN
          IF (ionode.AND.logger%print_keys%scf) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                  "*** SCF run converged in ",iter_count," steps ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       ELSE IF (should_stop.OR.&
            iter_count == xas_control%max_scf) THEN
          IF (ionode.AND.logger%print_keys%scf) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       END IF

       IF (.NOT.logger%print_keys%timecumul) t1 = m_walltime()

!   *** mixing methods have the new density matrix in p_mix_new
       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
           DO ispin=1,dft_control%nspins
              CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                                       rho%rho_ao(ispin)%matrix)
           END DO
       ENDIF


! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=gapw, error=error)
       IF(gapw) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
       ENDIF

       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)      

       CALL timestop(0.0_dp,handle2)

    END DO scf_loop
    IF (ionode.AND.logger%print_keys%scf) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
         IF(gapw) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
              "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
              "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
              "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
         END IF
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy:                                  ",energy%total
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Ionization potential of the excited atom:      ",IP_energy
       CALL m_flush(output_unit)
    END IF

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    DEALLOCATE(vecbuffer,vecbuffer2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    !Still not sure if this is really possible
!dbg
!     write(*,*) xas_env%scf_method, xas_scf_general
!dbg
    IF(.NOT.xas_env%scf_method==xas_scf_general) THEN
    ! Rotate the wfn to get the eugenstate of the KS hamiltonian
    ! Only ispin=1 should be needed
      DO ispin=1,dft_control%nspins
        CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
             eigenvalues=eigenvalues,homo=homo)
        CALL calculate_subspace_eigenvalues(mo_coeff,&
             matrix_ks(ispin)%matrix,eigenvalues, &
             para_env=qs_env%para_env, &
             do_rotation=.TRUE.)
      END DO  ! ispin
    END IF

    !Calculate the virtual states from the KS matrix matrix_ks(1)
!dbg
!     write(*,*) 'nv2 ', nvirtual2
!dbg
    IF(nvirtual2 .GT. 0) THEN
      NULLIFY(mo_coeff,uno_orbs,uno_evals,preconditioner)
      ALLOCATE(preconditioner,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL init_preconditioner(preconditioner,para_env=qs_env%para_env,&
           blacs_env=qs_env%blacs_env)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo)
      CALL make_preconditioner(preconditioner, &
                  precon_type=ot_precond_full_kinetic, &
                  matrix_h=matrix_ks(1)%matrix,matrix_s=matrix_s(1)%matrix,&
                  matrix_t=qs_env%kinetic(1)%matrix, &
                  mo_coeff=mo_coeff,energy_gap=0.2_dp)

      CALL get_xas_env(xas_env=xas_env,unoccupied_orbs=uno_orbs,&
           unoccupied_evals=uno_evals,unoccupied_eps=uno_eps,unoccupied_max_iter=uno_iter)
      CALL cp_fm_init_random(uno_orbs,nvirtual2)
      CALL ot_eigensolver(matrix_h=matrix_ks(1)%matrix,matrix_s=matrix_s(1)%matrix, &
           matrix_c=uno_orbs,matrix_orthogonal_space=mo_coeff,&
           preconditioner=preconditioner,eps_gradient=uno_eps,&
           iter_max=uno_iter,globenv=globenv,size_ortho_space=nmo)
      CALL calculate_subspace_eigenvalues(uno_orbs,matrix_ks(1)%matrix,&
           uno_evals,para_env=qs_env%para_env,do_rotation=.TRUE.)
      CALL destroy_preconditioner(preconditioner)
      DEALLOCATE(preconditioner,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! Prapare arrays for the calculation of the spectrum
    IF(.NOT. xas_control%xas_method == xas_dscf) THEN
      ! Copy the final vectors in the array
      NULLIFY(all_vectors,all_evals)
      CALL get_xas_env(xas_env=xas_env,all_vectors=all_vectors,&
           all_evals=all_evals,nvirtual1=nvirtual1)
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,mo_coeff=mo_coeff,&
           nao=nao,nmo=nmo,homo=homo)

      IF(xas_control%emission) THEN
        CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nvirtual1,&
             source_start=1,target_start=1)
      ELSE
        CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nvirtual1,&
             source_start=homo+1,target_start=1)
      END IF
      IF(nvirtual2 .GT. 0) THEN
        CALL cp_fm_to_fm(uno_orbs,all_vectors,ncol=nvirtual2,&
             source_start=1,target_start=1+nvirtual1)
      END IF

      IF(xas_control%emission) THEN
        DO istate = 1,nvirtual1
          all_evals(istate) = eigenvalues(istate)
        ENDDO
      ELSE
        DO istate = 1,nvirtual1
          all_evals(istate) = eigenvalues(homo+istate)
        ENDDO
      END IF
      DO istate = 1,nvirtual2
        all_evals(istate+nvirtual1) = uno_evals(istate)
      END DO 
    END IF

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv)


  END SUBROUTINE xas_do_tp_scf

END MODULE xas_tp_scf
