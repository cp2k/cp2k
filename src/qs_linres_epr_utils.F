!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
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
  USE coefficient_types,               ONLY: coeff_zero
  USE convert_units,                   ONLY: convert_to_cp2k_units
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
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_get_next_line,&
                                             parser_get_object,&
                                             parser_release
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_constants,                 ONLY: restart_guess
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: fourpi !R
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_methods,               ONLY: linres_scf
  USE qs_linres_epr_atom_current,      ONLY: list_3c_jrho_atom
  USE qs_linres_epr_current,           ONLY: epr_response_current
  USE qs_linres_epr_op,                ONLY: set_vecp
!RUSE qs_linres_epr_shift,             ONLY: epr_print_shift
  USE qs_linres_types,                 ONLY: deallocate_jrho_atom_set,&
                                             deallocate_nablavks_atom_set,&
                                             get_epr_env,&
                                             init_jrho_atom_set,&
                                             init_nablavks_atom_set,&
                                             jrho_atom_type,&
                                             nablavks_atom_type,&
                                             linres_control_type,&
                                             epr_env_create,&
                                             epr_env_type,&
                                             set_epr_env
  USE qs_loc_methods,                  ONLY: qs_print_cubes
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: set_up_op_sm
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain
  USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals
  USE qs_rho_atom_types,               ONLY: rho_atom_type,&
                                             deallocate_rho_atom_set !R
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

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_do_epr', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, idir, ii, iii, ispin, ist_true, istat, istate, &
      istate2, nao, ncubes, nmo, nstates(2), output_unit
    INTEGER, DIMENSION(:), POINTER           :: list_cubes
    INTEGER, DIMENSION(:, :), POINTER        :: statetrueindex
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: dk(3), dkl(3), dl(3)
    REAL(dp), DIMENSION(:), POINTER          :: dkl_vec_ii, dkl_vec_iii
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_dklxp0
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: h1_psi0, psi1
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: p_psi0, psi1_D, psi1_p, &
                                                psi1_rxp, rxp_psi0
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_work_ii, fm_work_iii, &
                                                mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: lr_section, epr_section

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
!       CALL epr_restart_epr(epr_env,qs_env,action="read",error=error)
    END IF

    CALL get_epr_env(epr_env=epr_env, nao=nao, nstates=nstates, &
                     centers_set=centers_set, list_cubes=list_cubes,&
                     statetrueindex=statetrueindex, &
                     psi1_p=psi1_p, psi1_rxp=psi1_rxp, psi1_D=psi1_D,&
                     p_psi0=p_psi0, rxp_psi0=rxp_psi0)

    ! operator p
    DO idir = 1,3
       CALL cp_fm_set_all(psi1_p(1,idir)%matrix,0.0_dp,error=error)
    END DO 
    DO idir = 1,3
      IF(.NOT. epr_env%simpleepr_done(idir)) THEN

         IF(output_unit>0) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,I4/)")&
           "Response to the perturbation operator p in dir. ", idir
         END IF

        !Initial guess for psi1
        !Maybe it could be better than this
      !   IF(idir==1) THEN
           DO ispin = 1,dft_control%nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
      !       CALL cp_fm_to_fm(p_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
      !       CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
           END DO 
      !   END IF

         !DO scf cycle to optimize psi1
         DO ispin = 1,dft_control%nspins
           CALL cp_fm_to_fm(p_psi0(ispin,idir)%matrix, h1_psi0(ispin)%matrix,error=error)
         END DO
         linres_control%converged = .FALSE.
         CALL linres_scf(p_env, qs_env, psi1, h1_psi0,error=error)

!         CALL test_func(qs_env,psi1,epr_env,error)

         IF(linres_control%converged) THEN

          ! the optimized wfns are copied in the  fm
           IF(output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
             "Store the psi1 for the calculation of the response current density "
           END IF
           DO ispin=1, dft_control%nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix, psi1_p(ispin,idir)%matrix,error=error)
           END DO

         ! print response functions
           IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
              "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
              ncubes = SIZE(list_cubes,1)
              DO ispin = 1,dft_control%nspins
                CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
                     centers_set(ispin)%array,ionode,'psi1_p',&
                     idir=idir,ispin=ispin,stride=section_get_ival(epr_section,&
                     "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE",error=error),error=error)
              END DO  ! ispin
           END IF  ! print response functions
         END IF

         IF(output_unit>0) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
            "Write the resulting psi1 in restart file "
         END IF
         ! Write the result inthe restart file
         epr_env%simpleepr_done(idir) = .TRUE.
         epr_env%simpleepr_converged(idir) = linres_control%converged
!        CALL epr_restart_epr(epr_env,qs_env,action="write",ivec=idir,error=error)

      END IF  
    END DO ! idir

    ! operator rxp
    DO idir = 1,3
       CALL cp_fm_set_all(psi1_rxp(1,idir)%matrix,0.0_dp,error=error)
    END DO 

    DO idir = 1,3 

      IF(.NOT. epr_env%simpleepr_done(idir+3)) THEN

         IF(output_unit>0) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,I4/)")&
           "Response to the perturbation operator rxp in dir. ", idir
         END IF

         !Initial guess for psi1
         !Maybe it could be better than this
         DO ispin = 1,dft_control%nspins
           CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
!             CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
!             CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
         END DO 

         !DO scf cycle to optimize psi1
         DO ispin = 1,dft_control%nspins
            CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, h1_psi0(ispin)%matrix,error=error)
!dpg
            ! print perturbation operator
!            ncubes = SIZE(list_cubes,1)
!            CALL qs_print_cubes(qs_env,h1_psi0(ispin)%matrix,ncubes,list_cubes,&
!                     centers_set(ispin)%array,ionode,'pert_rxp',&
!                     idir=idir,ispin=ispin, stride=section_get_ival(epr_section,&
!                     "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE"),error=error)
!dpg
         END DO

         ! Optimize the response wavefunctions 
         linres_control%converged = .FALSE.
         CALL linres_scf(p_env, qs_env, psi1, h1_psi0,error=error)

         IF(linres_control%converged) THEN

           !   the optimized wfns are copied in the  fm
           IF(output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
             "Store the psi1 for the calculation of the response current density "
           END IF
           DO ispin=1, dft_control%nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix, psi1_rxp(ispin,idir)%matrix,error=error)
           END DO
           ! print response functions
           IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
              "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
              ncubes = SIZE(list_cubes,1)
              DO ispin = 1,dft_control%nspins
                CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
                     centers_set(ispin)%array,ionode,'psi1_rxp',&
                     idir=idir,ispin=ispin,stride=section_get_ival(epr_section,&
                     "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE",error=error),error=error)
              END DO  ! ispin
           END IF  ! print response functions

         END IF

         IF(output_unit>0) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
            "Write the resulting psi1 in restart file "
         END IF
         ! Write the result inthe restart file
         epr_env%simpleepr_done(idir+3) = .TRUE.
!        CALL epr_restart_epr(epr_env,qs_env,action="write",ivec=idir+3,error=error)
        END IF
     END DO ! idir


     ! operator D
     IF(epr_env%full_epr) THEN
       DO ispin = 1,dft_control%nspins
         CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
       END DO 
      DO idir = 1,3
         CALL cp_fm_set_all(psi1_D(1,idir)%matrix,0.0_dp,error=error)
      END DO 


!      The correction is state depedent a loop over the states is necessary
       ALLOCATE(vecbuf_dklxp0(1,nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       vecbuf_dklxp0(1,nao) = 0.0_dp

       DO ispin=1, dft_control%nspins

         CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,homo=homo)
         NULLIFY(tmp_fm_struct)
         CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
               ncol_global=nmo,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
         CALL cp_fm_create (fm_work_ii, tmp_fm_struct,error=error )
         CALL cp_fm_set_all(fm_work_ii,0.0_dp,error=error)
         CALL cp_fm_create (fm_work_iii, tmp_fm_struct,error=error )
         CALL cp_fm_set_all(fm_work_iii,0.0_dp,error=error)
         CALL cp_fm_struct_release ( tmp_fm_struct,error=error )

         ALLOCATE(dkl_vec_ii(nstates(ispin)),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(dkl_vec_iii(nstates(ispin)),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO idir = 1,3
           DO istate = 1,nstates(ispin)
             ist_true = statetrueindex(idir,istate)
             !the initial guess is the previous set of psi1, just optimized
             CALL set_vecp(idir,ii,iii)
             dk(1:3) = centers_set(ispin)%array(1:3,ist_true)
              
             IF(output_unit>0) THEN
                WRITE (UNIT=output_unit,FMT="(T10,A,I4,A,I4/)")&
               "Response to the perturbation operator (dk-dl)xp for the k state ",&
                ist_true, " in dir. ", idir
             END IF

             DO istate2 = 1,nstates(ispin)
                dl(1:3)  = centers_set(ispin)%array(1:3,istate2)
                dkl = pbc(dl,dk,cell)
                dkl_vec_ii(istate2)  = dkl(ii)
                dkl_vec_iii(istate2) = dkl(iii)
             END DO

         ! First term
             ! Rescale the ground state orbitals by (dk-dl)_ii
             CALL cp_fm_to_fm(mo_coeff,fm_work_ii,error=error)
             CALL cp_fm_column_scale(fm_work_ii,dkl_vec_ii(1:homo)) 
             ! Apply the p_iii operator
             CALL cp_sm_fm_multiply(epr_env%op_p_ao(iii)%matrix,fm_work_ii,&
                   fm_work_iii,ncol=nmo,alpha=-1.0_dp,error=error)
             ! Copy in h1_psi1
             CALL cp_fm_to_fm(fm_work_iii,h1_psi0(ispin)%matrix,error=error)


         ! Second term
             ! Rescale the ground state orbitals by (dk-dl)_iii
             CALL cp_fm_to_fm(mo_coeff,fm_work_iii,error=error)
             CALL cp_fm_column_scale(fm_work_iii,dkl_vec_iii(1:homo)) 
             ! Apply the p_iii operator
             CALL cp_sm_fm_multiply(epr_env%op_p_ao(ii)%matrix,fm_work_iii,&
                   fm_work_ii,ncol=nmo,alpha=-1.0_dp,error=error)
             ! Copy in h1_psi1
             CALL cp_fm_scale_and_add(1.0_dp,h1_psi0(ispin)%matrix,-1.0_dp,fm_work_ii,error=error)

             ! Optimize the response wavefunctions 
             CALL linres_scf(p_env, qs_env, psi1, h1_psi0,error=error)

             IF(output_unit>0) THEN
               WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
               "Store the psi1 vector for the calculation of the response current density "
             END IF
             ! the optimized wfns are copied in the  fm
             CALL cp_fm_get_submatrix(psi1(ispin)%matrix,vecbuf_dklxp0 ,&
                    1,ist_true,nao,1,transpose=.TRUE.,error=error)
             CALL cp_fm_set_submatrix(psi1_D(ispin,idir)%matrix, vecbuf_dklxp0,&
                    1,ist_true,nao,1,transpose=.TRUE.,error=error)
             epr_env%fullepr_done(idir*ispin,istate) = .TRUE.

           END DO  ! istate
           ! print response functions
           IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
              "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
              ncubes = SIZE(list_cubes,1)
               CALL qs_print_cubes(qs_env,psi1_D(ispin,idir)%matrix,&
                    ncubes,list_cubes,&
                    centers_set(ispin)%array,ionode,'psi1_D',&
                    idir=idir,ispin=ispin,stride=section_get_ival(epr_section,&
                    "PRINT%RESPONSE_FUNCTION_CUBES%STRIDE",error=error),error=error)
           END IF  ! print response functions

!          CALL epr_restart_epr(epr_env,qs_env,action="write",ivec=idir+6,error=error)
          
         END DO  ! idir

         CALL cp_fm_release(fm_work_ii,error=error)
         CALL cp_fm_release(fm_work_iii,error=error)
         DEALLOCATE(dkl_vec_ii,dkl_vec_iii,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO  ! ispin
     END IF

    !If current density and shift are calculated here
    ! probably a lot of memory is required
    ! Better to free what is not anymore necessary in p_env
     CALL p_env_release(p_env,error=error)
 
    !Calculate current density and induced field
     CALL epr_response_current(epr_env,qs_env,psi1,h1_psi0,error=error)

    ! Print Results
!R   CALL epr_print_shift(epr_env,qs_env,epr_section,error=error)
!R STOP 'wait'
    !----------!
    ! clean up !
    !----------!
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, psi1,error=error)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, h1_psi0,error=error)
     NULLIFY(psi1,h1_psi0)
     CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
          "PRINT%PROGRAM_RUN_INFO",error=error)     

    CALL timestop(0.0_dp, handle)

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

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_env_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: label
!R  CHARACTER(LEN=default_string_length)     :: nics_file_name
    INTEGER :: handle, homo, i_B, icount, idir, ini, ir, ispin, istat, &
      istate, istate2, istate_next, j, k, n_mo(2), n_rep, nao, nat_print, &
      natom, ncubes, nmoloc, nspins, nstates, output_unit, i1 !R
    INTEGER, DIMENSION(:), POINTER           :: bounds, list
    LOGICAL                                  :: failure, gapw, ionode, my_end
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: state_done
    REAL(dp)                                 :: center(3), center2(3), dist, &
                                                mdist, rab(3)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(rho_atom_type), DIMENSION(:), & !R
      POINTER                                :: vks_atom_set !R
    TYPE(nablavks_atom_type), DIMENSION(:), & !R
      POINTER                                :: nablavks_atom_set !R
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
    TYPE(section_vals_type), POINTER         :: lr_section, epr_section

    CALL timeset("epr_env_init","I"," ",handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, cell, dft_control, linres_control, scf_control)
    NULLIFY(logger, matrix_s, mos, mpools, epr_section, particle_set)
    NULLIFY(auxbas_pw_pool,pw_env)

    n_mo(1:2) = 0
    nao = 0
    nmoloc = 0

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    IF(epr_env%ref_count /= 0) THEN
      CALL epr_env_cleanup(epr_env,error=error)
    END IF


    IF(output_unit>0) THEN

      WRITE (UNIT=output_unit,FMT="(/,T20,A,/)")&
         "*** Start EPR g tensor Calculation ***"
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "Initialization of the EPR environment"

    END IF

    CALL epr_env_create(epr_env,error=error)

    ! If surrent_density or full_epr different allocations are required
    epr_section => section_vals_get_subs_vals(qs_env%input, &
                             "PROPERTIES%LINRES%EPR",error=error)
    CALL section_vals_val_get(epr_section,"FULL",l_val=epr_env%full_epr,error=error)
    CALL section_vals_val_get(epr_section,"RESTART_EPR",l_val=epr_env%restart_epr,error=error)
!R  CALL section_vals_val_get(epr_section,"NICS",l_val=epr_env%do_nics,error=error)
!R  IF(epr_env%do_nics)THEN
!R    CALL section_vals_val_get(epr_section,"NICS_FILE_NAME",&
!R         c_val=nics_file_name,error=error)
!R    CALL parser_create(parser,nics_file_name,error=error)
!R    CALL parser_get_next_line(parser,1,error=error)
!R    CALL parser_get_object  (parser,epr_env%n_nics)
!R    ALLOCATE(epr_env%r_nics(3,epr_env%n_nics),STAT=istat)
!R    CALL parser_get_next_line(parser,2,error=error)
!R    DO j = 1,epr_env%n_nics
!R      CALL parser_get_object(parser,label)
!R      CALL parser_get_object(parser,epr_env%r_nics(1,j))
!R      CALL parser_get_object(parser,epr_env%r_nics(2,j))
!R      CALL parser_get_object(parser,epr_env%r_nics(3,j))
!R      CALL convert_to_cp2k_units("ANGSTROM",length=epr_env%r_nics(1,j))
!R      CALL convert_to_cp2k_units("ANGSTROM",length=epr_env%r_nics(2,j))
!R      CALL convert_to_cp2k_units("ANGSTROM",length=epr_env%r_nics(3,j))
!R      CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
!R      IF (my_end) EXIT
!R    END DO
!R    CALL parser_release(parser,error=error)
!R  END IF


    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    matrix_s=matrix_s,&
                    mos=mos,&
                    mpools=mpools,&
                    particle_set=particle_set,&
                    pw_env=pw_env,&
                    scf_control=scf_control,error=error)

    ! Check if restat also psi0 should be restarted
    IF(epr_env%restart_epr .AND. scf_control%density_guess/=restart_guess)THEN
      CALL stop_program("epr_env_init","restart_epr requires density_guess=restart")
    END IF

    ! check that the psi0 are localized and you have all the centers
     CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
     IF(failure .AND. (output_unit>0)) THEN
       WRITE (UNIT=output_unit,FMT='(A)') &
          ' To get EPR parameters within PBC you need localized zero order orbitals '
     END IF
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
           linres_control%localized_wfn_control%centers_set(ispin)%array(idir,:)
       END DO
       epr_env%nstates(ispin) = homo
       nstates = MAX(nstates,homo)
     END DO  
     epr_env%nao =nao
     ALLOCATE(epr_env%statetrueindex(3,nstates),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     epr_env%statetrueindex = 0
     ALLOCATE(state_done(3,nstates),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     state_done = .FALSE.

     epr_env%statetrueindex(1,1) = 1
     center(1)  = epr_env%centers_set(1)%array(1,1)
     center(2)  = epr_env%centers_set(1)%array(2,1)
     center(3)  = epr_env%centers_set(1)%array(3,1)
     state_done(1,1) = .TRUE. 
     icount = 1

     DO idir = 1,3
       ini = 1
       IF(idir == 1) ini = 2
       DO istate = ini, nstates
         mdist = 100.0_dp

         DO istate2 = 1,nstates
           IF(.NOT.state_done(idir,istate2)) THEN
             center2(1)  = epr_env%centers_set(1)%array(1,istate2)
             center2(2)  = epr_env%centers_set(1)%array(2,istate2)
             center2(3)  = epr_env%centers_set(1)%array(3,istate2)

             rab = pbc(center,center2,cell)
             CALL set_vecp(idir,j,k)
             dist = SQRT(rab(j)*rab(j)+rab(k)*rab(k))

             IF(dist .LT. mdist) THEN
                mdist =dist
                istate_next = istate2
             END IF
           END IF
         END DO  ! istate2

         icount = icount + 1
         state_done(idir,istate_next) = .TRUE. 
         epr_env%statetrueindex(idir,icount) = istate_next

         center(1)  = epr_env%centers_set(1)%array(1,istate_next)
         center(2)  = epr_env%centers_set(1)%array(2,istate_next)
         center(3)  = epr_env%centers_set(1)%array(3,istate_next)
       END DO  ! istate
       icount = 0
     END DO  ! idir
     DEALLOCATE(state_done,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


!    Conversion factors
     epr_env%g_free_factor = 2.0023192778_dp
     epr_env%g_zke_factor = -1.0_dp * 2.0023192778_dp * ( 1.0_dp / 137.03602_dp )**2
     epr_env%g_so_factor =  ( 1.0_dp / 137.03602_dp )**2 * ( 2.0023192778_dp - 1.0_dp )
     epr_env%g_soo_factor = 2.0_dp * 2.0_dp * ( 1.0_dp / 137.03602_dp ) * fourpi * ( 1.0_dp / 137.03602_dp )**2
!R!    factor for the CHEMICAL SHIFTS: alpha^2 *  ppm.
!R     epr_env%shift_factor = ( 1.0_dp / 137.03602_dp )**2   * 1.0E+6_dp / cell%deth
!R!    chi_factor =  1/4 * e^2/m * a_0 ^2
!R     epr_env%chi_factor = 1.9727566E-29_dp / 1.0E-30_dp ! -> displayed in 10^-30 J/T^2
!R!    Factor to convert 10^-30 J/T^2 into ppm cgs = ppm cm^3/mol
!R!      = 10^-30 * mu_0/4pi * N_A * 10^6 * 10^6  [one 10^6 for ppm, one for m^3 -> cm^3]
!R     epr_env%chi_SI2ppmcgs = 6.022045_dp/1.0E+2_dp 
!R!    Chi to Shift: 10^-30  *  2/3  mu_0 / Omega  * 1/ppm
!R     epr_env%chi_SI2shiftppm = 1.0E-30_dp * 8.37758041E-7_dp/ &
!R                      (cell%deth /cell%unit_of_length**3*1.0E-30_dp) * 1.0E+6_dp

     IF(output_unit>0) THEN
       IF(epr_env%full_epr) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "  EPR|","Full orbital dependent correction"
       END IF
!R     IF(epr_env%do_nics) THEN
!R       WRITE (UNIT=output_unit,FMT="(T2,A,T50,I5,A)")&
!R          "  EPR|   NICS computed in ",epr_env%n_nics," additional points"
!R       WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
!R          "  EPR|   NICS coordinates read on file ",nics_file_name
!R     END IF
!R     WRITE (UNIT=output_unit,FMT="(T2,A,T60,ES15.6)")&
!R          "  EPR|      Shift factor (ppm)", epr_env%shift_factor
!R     WRITE (UNIT=output_unit,FMT="(T2,A,T60,ES15.6)")&
!R          "  EPR|        Chi factor  (SI)", epr_env%chi_factor
!R     WRITE (UNIT=output_unit,FMT="(T2,A,T60,ES15.6)")&
!R          "  EPR| Conversion Chi (ppm/cgs)", epr_env%chi_SI2ppmcgs
!R     WRITE (UNIT=output_unit,FMT="(T2,A,T60,ES15.6)")&
!R          "  EPR| Conversion Chi to Shift", epr_env%chi_SI2shiftppm
     END IF

!R   epr_env%store_current = .FALSE.
     IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
                  "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN

       NULLIFY(bounds,list)
       ncubes = 0
       CALL section_vals_val_get(epr_section,&
          "PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LU_BOUNDS",&
           i_vals=bounds,error=error) 
       ncubes = bounds(2) - bounds(1)  + 1
       IF(ncubes > 0 ) THEN
         ALLOCATE( epr_env%list_cubes(ncubes),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO ir = 1,ncubes
           epr_env%list_cubes(ir) = bounds(1) + (ir-1)
         END DO
       END IF
       IF(.NOT. ASSOCIATED(epr_env%list_cubes)) THEN
         CALL section_vals_val_get(epr_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
             n_rep_val=n_rep,error=error)
         ncubes = 0
         DO ir = 1,n_rep
          NULLIFY(list)
          CALL section_vals_val_get(epr_section,"PRINT%RESPONSE_FUNCTION_CUBES%CUBES_LIST",&
             i_rep_val=ir,i_vals=list,error=error)
          IF(ASSOCIATED(list)) THEN
            CALL reallocate(epr_env%list_cubes,1,ncubes+ SIZE(list))
            DO ini = 1, SIZE(list)
              epr_env%list_cubes(ini+ncubes) = list(ini)
            END DO
            ncubes = ncubes + SIZE(list)
          END IF
         END DO  ! ir
       END IF
       IF(.NOT. ASSOCIATED(epr_env%list_cubes)) THEN
         ALLOCATE( epr_env%list_cubes(nstates),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO ir = 1,nstates
           epr_env%list_cubes(ir) = ir
         END DO
       END IF
     END IF
!R   IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
!R                "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
!R     epr_env%store_current = .TRUE.
!R   END IF

!R   ALLOCATE (epr_env%do_calc_cs_atom(natom), STAT=istat)
!R   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R   epr_env%do_calc_cs_atom = 0
!R
!R   IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
!R                "PRINT%SHIELDING_TENSOR",error=error),cp_p_file)) THEN
!R
!R     NULLIFY(bounds,list)
!R     nat_print = 0
!R     CALL section_vals_val_get(epr_section,&
!R        "PRINT%SHIELDING_TENSOR%ATOMS_LU_BOUNDS",&
!R         i_vals=bounds,error=error) 
!R     nat_print = bounds(2) - bounds(1)  + 1
!R     IF(nat_print > 0) THEN
!R       ALLOCATE(epr_env%cs_atom_list(nat_print),STAT=istat)
!R       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R       DO ir = 1,nat_print
!R         epr_env%cs_atom_list(ir) = bounds(1) + (ir-1)
!R         epr_env%do_calc_cs_atom(bounds(1) + (ir-1)) = 1
!R       END DO
!R     END IF
!R
!R     IF(.NOT. ASSOCIATED(epr_env%cs_atom_list)) THEN
!R       CALL section_vals_val_get(epr_section,"PRINT%SHIELDING_TENSOR%ATOMS_LIST",&
!R           n_rep_val=n_rep,error=error)
!R       nat_print = 0
!R       DO ir = 1,n_rep
!R        NULLIFY(list)
!R        CALL section_vals_val_get(epr_section,"PRINT%SHIELDING_TENSOR%ATOMS_LIST",&
!R           i_rep_val=ir,i_vals=list,error=error)
!R        IF(ASSOCIATED(list)) THEN
!R          CALL reallocate(epr_env%cs_atom_list,1,nat_print + SIZE(list))
!R          DO ini = 1, SIZE(list)
!R            epr_env%cs_atom_list(ini+nat_print) = list(ini)
!R            epr_env%do_calc_cs_atom(list(ini)) = 1
!R          END DO
!R          nat_print = nat_print + SIZE(list)
!R        END IF
!R       END DO  ! ir
!R     END IF
!R
!R     IF(.NOT. ASSOCIATED(epr_env%cs_atom_list)) THEN
!R       ALLOCATE(epr_env%cs_atom_list(natom),STAT=istat)
!R       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R       DO ir = 1,natom
!R         epr_env%cs_atom_list(ir) = ir
!R       END DO
!R       epr_env%do_calc_cs_atom = 1
!R     END IF
!R
!R   ELSE
!R      NULLIFY(epr_env%cs_atom_list)
!R   END IF

!R   IF(output_unit>0) THEN
!R     IF(ASSOCIATED(epr_env%cs_atom_list)) THEN
!R       WRITE (UNIT=output_unit,FMT="(T2,A,T40,I5,A)")&
!R          "  EPR|   Shielding tensor computed for ", SIZE(epr_env%cs_atom_list,1) ," atoms"
!R     ELSE
!R       WRITE (UNIT=output_unit,FMT="(T2,A,T50)")&
!R          "  EPR|   Shielding tensor not computed at the atomic positions"
!R
!R     END IF
!R   END IF

     CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=epr_env%ao_mo_fm_pools, error=error)
    ! for the g tensor we need 6 psi1, i.e. 6 optimization procedures
    ! They become 9 if full epr is calculated, i.e. with the correction term too
    ! All of them are required at the end of the optimization procedure
    ! if the current density and the induced fiels have to be calculated 

     ALLOCATE(epr_env%psi1_p(nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(epr_env%psi1_rxp(nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO ispin = 1,nspins
       DO idir = 1,3
         NULLIFY(epr_env%psi1_p(ispin,idir)%matrix)
         CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
              epr_env%psi1_p(ispin,idir)%matrix,error=error)
         NULLIFY(epr_env%psi1_rxp(ispin,idir)%matrix)
         CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
              epr_env%psi1_rxp(ispin,idir)%matrix,error=error)
       END DO 
     END DO 

     IF(epr_env%full_epr) THEN
       ALLOCATE(epr_env%psi1_D(nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
         DO idir = 1,3
           NULLIFY(epr_env%psi1_D(ispin,idir)%matrix)
           CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
                epr_env%psi1_D(ispin,idir)%matrix,error=error)
         END DO
       END DO
     END IF

     ! px py pz operator acting on the atomic orbitals
     CALL allocate_matrix_set(epr_env%op_p_ao,3,error=error)
     CALL set_up_op_sm(epr_env%op_p_ao(1)%matrix,qs_env,symmetry="antisymmetric",&
          name="OP_P",error=error)
     CALL set_matrix(epr_env%op_p_ao(1)%matrix,0.0_dp)
     DO idir = 2,3 
       CALL replicate_matrix_structure(epr_env%op_p_ao(1)%matrix, &
            epr_env%op_p_ao(idir)%matrix,&
            "epr_env%op_p_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),&
            target_symmetry="antisymmetric",error=error)
       CALL set_matrix(epr_env%op_p_ao(idir)%matrix,0.0_dp)
     END DO

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
              epr_env%p_psi0(ispin,idir)%matrix,error=error)
         NULLIFY(epr_env%rxp_psi0(ispin,idir)%matrix)
         CALL fm_pool_create_fm(epr_env%ao_mo_fm_pools(ispin)%pool,&
              epr_env%rxp_psi0(ispin,idir)%matrix,error=error)
       END DO 
     END DO 

     CALL allocate_matrix_set(epr_env%op_rmd_ao,3,error=error)
     CALL set_up_op_sm(epr_env%op_rmd_ao(1)%matrix,qs_env,symmetry="none",&
          name="OP_RMD",error=error)
     CALL set_matrix(epr_env%op_rmd_ao(1)%matrix,0.0_dp)
     DO idir = 2,3 
       CALL replicate_matrix_structure( epr_env%op_rmd_ao(1)%matrix, &
            epr_env%op_rmd_ao(idir)%matrix,"epr_env%op_rmd_ao"//&
            "-"//TRIM(ADJUSTL(cp_to_string(idir))),target_symmetry="none",error=error)
       CALL set_matrix(epr_env%op_rmd_ao(idir)%matrix,0.0_dp)
     END DO

    ! Current density matrix in x y and z:
     CALL replicate_matrix_structure( epr_env%op_rmd_ao(1)%matrix, &
           epr_env%jp1_ao,"epr_env%jp1_ao",target_symmetry="none",error=error)
     CALL set_matrix(epr_env%jp1_ao,0.0_dp)

     ALLOCATE( epr_env%jp2_ao(2),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO idir = 1,2
       NULLIFY(epr_env%jp2_ao(idir)%matrix)
       CALL replicate_matrix_structure( epr_env%op_rmd_ao(1)%matrix, &
          epr_env%jp2_ao(idir)%matrix,"epr_env%jp2_ao"//&
          "-"//TRIM(ADJUSTL(cp_to_string(idir))),target_symmetry="none",error=error)
     END DO

     ! If the current density on the grid needs to be stored
!R   IF(epr_env%store_current) THEN
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
               CALL pw_pool_init_coeff(auxbas_pw_pool,&
                    epr_env%jrho1_set(idir,i_B)%rho%rho_r(ispin),&
                    use_data=REALDATA3D,in_space=REALSPACE,error=error)
               CALL coeff_zero(epr_env%jrho1_set(idir,i_B)%rho%rho_r(ispin))

               CALL pw_pool_init_coeff(auxbas_pw_pool,&
                    epr_env%jrho1_set(idir,i_B)%rho%rho_g(ispin),&
                    use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
               CALL coeff_zero(epr_env%jrho1_set(idir,i_B)%rho%rho_g(ispin))
               
           END DO
           CALL qs_rho_retain(epr_env%jrho1_set(idir,i_B)%rho,error=error)
         END DO 
       END DO
         
!R   END IF

     ! Initialize local current density if GAPW calculation
     gapw = dft_control%qs_control%gapw
     IF(gapw) THEN
       CALL list_3c_jrho_atom(qs_env,error=error)
       CALL init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,&
            nspins,.TRUE.,error=error) !R
       CALL init_nablavks_atom_set(nablavks_atom_set,atomic_kind_set,&
            nspins,.TRUE.,error=error) !R

       CALL allocate_rho_atom_internals(qs_env, vks_atom_set,error=error) !R
       CALL set_epr_env(epr_env=epr_env, jrho1_atom_set=jrho1_atom_set, nablavks_atom_set=nablavks_atom_set, &
                        vks_atom_set=vks_atom_set, error=error) !R
     END IF

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
         CALL pw_pool_init_coeff(auxbas_pw_pool,&
              epr_env%bind_set(idir,i_B)%rho%rho_r(1),&
              use_data=REALDATA3D,in_space=REALSPACE,error=error)
         CALL pw_pool_init_coeff(auxbas_pw_pool,&
              epr_env%bind_set(idir,i_B)%rho%rho_g(1),&
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
         CALL pw_pool_init_coeff(auxbas_pw_pool,&
              epr_env%nablavks_set(idir,ispin)%rho%rho_r(1),&
              use_data=REALDATA3D,in_space=REALSPACE,error=error)
         CALL pw_pool_init_coeff(auxbas_pw_pool,&
              epr_env%nablavks_set(idir,ispin)%rho%rho_g(1),&
              use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
       END DO
     END DO

    ! Initialize the chemical shift tensor
!R  epr_env%chi_tensor = 0.0_dp
!R  ALLOCATE(epr_env%chemical_shift(3,3,natom),STAT=istat)
!R  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R  epr_env%chemical_shift = 0.0_dp
!R  ALLOCATE(epr_env%chemical_shift_loc(3,3,natom),STAT=istat)
!R  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R  epr_env%chemical_shift = 0.0_dp
!R  IF(epr_env%do_nics) THEN
!R    ALLOCATE(epr_env%chemical_shift_loc_nics(3,3,epr_env%n_nics),STAT=istat)
!R    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    epr_env%chemical_shift_loc_nics = 0.0_dp
!R    ALLOCATE(epr_env%chemical_shift_nics(3,3,epr_env%n_nics),STAT=istat)
!R    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    epr_env%chemical_shift_nics = 0.0_dp
!R  END IF

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
    DO i1 = 1,3
       epr_env%g_total(i1,i1) = epr_env%g_free_factor
    END DO
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
         "PRINT%PROGRAM_RUN_INFO",error=error)

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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_B, idir, ispin, istat
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
!R      IF(ASSOCIATED(epr_env%cs_atom_list)) THEN
!R        DEALLOCATE(epr_env%cs_atom_list)
!R      END IF
!R      IF(ASSOCIATED(epr_env%do_calc_cs_atom)) THEN
!R        DEALLOCATE(epr_env%do_calc_cs_atom)
!R      END IF
        !op_p_ao
        IF(ASSOCIATED(epr_env%op_p_ao)) THEN
          CALL deallocate_matrix_set(epr_env%op_p_ao,error=error)
        END IF
        !op_rmd_ao
        IF(ASSOCIATED(epr_env%op_rmd_ao)) THEN
          CALL deallocate_matrix_set(epr_env%op_rmd_ao,error=error)
        END IF
        !jp1_ao 
        IF(ASSOCIATED(epr_env%jp1_ao)) THEN
          CALL deallocate_matrix(epr_env%jp1_ao,error=error)
        END IF
        !jp2_ao 
        IF(ASSOCIATED(epr_env%jp2_ao)) THEN
          DO idir = 1,2
              CALL deallocate_matrix(epr_env%jp2_ao(idir)%matrix,error=error)
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
!R      !chemical_shift
!R      IF(ASSOCIATED(epr_env%chemical_shift)) THEN
!R        DEALLOCATE(epr_env%chemical_shift,STAT=istat)
!R        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R      END IF
!R      IF(ASSOCIATED(epr_env%chemical_shift_loc)) THEN
!R        DEALLOCATE(epr_env%chemical_shift_loc,STAT=istat)
!R        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R      END IF
        IF(ASSOCIATED(epr_env%basisfun_center)) THEN
          DEALLOCATE(epr_env%basisfun_center,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(epr_env%statetrueindex)) THEN
          DEALLOCATE(epr_env%statetrueindex,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
!R      ! nics
!R      IF(ASSOCIATED(epr_env%r_nics)) THEN
!R        DEALLOCATE(epr_env%r_nics,STAT=istat)
!R        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R      END IF
!R      IF(ASSOCIATED(epr_env%chemical_shift_nics)) THEN
!R        DEALLOCATE(epr_env%chemical_shift_nics,STAT=istat)
!R        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R      END IF
!R      IF(ASSOCIATED(epr_env%chemical_shift_loc_nics)) THEN
!R        DEALLOCATE(epr_env%chemical_shift_loc_nics,STAT=istat)
!R        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R      END IF
        ! nablavks_set !R
        IF(ASSOCIATED(epr_env%nablavks_set)) THEN
          DEALLOCATE(epr_env%nablavks_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! nablavks_atom_set !R
        IF(ASSOCIATED(epr_env%nablavks_atom_set)) THEN
          CALL deallocate_nablavks_atom_set(epr_env%nablavks_atom_set,error=error)
        END IF
        ! vks_atom_set !R
        IF(ASSOCIATED(epr_env%vks_atom_set)) THEN
          CALL deallocate_rho_atom_set(epr_env%vks_atom_set)
        END IF
        ! bind_set !R
        IF(ASSOCIATED(epr_env%bind_set)) THEN
          DEALLOCATE(epr_env%bind_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! bind_atom_set !R
        IF(ASSOCIATED(epr_env%bind_atom_set)) THEN
          DEALLOCATE(epr_env%bind_atom_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! matrix_v !R
        IF(ASSOCIATED(epr_env%matrix_v)) THEN
          DEALLOCATE(epr_env%matrix_v,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_total !R
        IF(ASSOCIATED(epr_env%g_total)) THEN
          DEALLOCATE(epr_env%g_total,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_so !R
        IF(ASSOCIATED(epr_env%g_so)) THEN
          DEALLOCATE(epr_env%g_so,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_soo !R
        IF(ASSOCIATED(epr_env%g_soo)) THEN
          DEALLOCATE(epr_env%g_soo,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! bind0 !R
        IF(ASSOCIATED(epr_env%bind0)) THEN
          DEALLOCATE(epr_env%bind0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
      END IF  ! ref count
    END IF ! failure
    
  END SUBROUTINE epr_env_cleanup


END MODULE  qs_linres_epr_utils
