!   CP2K: A general program to perform molecular dynamics simulations         !

!!****** cp2k/qs_ot_scf [1.0] *
!!
!!   NAME
!!     qs_ot_scf
!!
!!   FUNCTION
!! 
!!     basic functionality for using ot in the scf routines.
!!
!!   AUTHOR
!!     Joost VandeVondele (25.08.2002)
!!
!!   MODIFICATION HISTORY
!!     01.2003 : Joost VandeVondele : adapted for LSD
!!
!!   SOURCE
!******************************************************************************

MODULE qs_ot_scf
! *****************************************************************************
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: ot_ls_2pnt,&
                                             ot_ls_3pnt,&
                                             ot_ls_gold,&
                                             ot_ls_none,&
                                             ot_mini_cg,&
                                             ot_mini_diis,&
                                             ot_mini_sd,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_full_all,&
                                             ot_precond_none,&
                                             ot_precond_s_inverse,&
                                             ot_precond_sparse_diag
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict
  USE qs_ot,                           ONLY: qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_get_orbitals,&
                                             qs_ot_get_p,&
                                             qs_ot_init,&
                                             qs_ot_settings_init,&
                                             qs_ot_type
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE


! *** Public subroutines ***

  PUBLIC :: ot_scf_init
  PUBLIC :: ot_scf_mini
  PUBLIC :: ot_scf_destroy
  PUBLIC :: ot_scf_read_input

!!***
! *****************************************************************************

CONTAINS

  SUBROUTINE ot_scf_read_input(qs_ot_env,globenv)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(global_environment_type), &
         pointer                             :: globenv

    CHARACTER(LEN=40)                        :: keyword, value
    INTEGER                                  :: handle, ispin, ls_method, &
         nspin, ot_method, &
         preconditioner_method
    LOGICAL                                  :: explicit, found
    TYPE(section_vals_type), POINTER         :: ot_section

    ! executable

    CALL timeset("ot_scf_read_input","I","",handle)

    ! decide default settings
    CALL qs_ot_settings_init(qs_ot_env(1)%settings)

    ! use ot input new style
    ot_section=>section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%DFT%SCF%OT")
    CALL section_vals_get(ot_section,explicit=explicit)

       CALL section_vals_val_get(ot_section,"MINIMIZER",i_val=ot_method)
       ! compatibility
       SELECT CASE(ot_method)
       CASE (ot_mini_sd)
         qs_ot_env(1)%settings%ot_method="SD"
       CASE (ot_mini_cg)
         qs_ot_env(1)%settings%ot_method="CG"
       CASE (ot_mini_diis)
         qs_ot_env(1)%settings%ot_method="DIIS"
         CALL section_vals_val_get(ot_section,"N_DIIS",i_val=qs_ot_env(1)%settings%diis_m)
       CASE DEFAULT
         CALL stop_program("READ OTSCF MINIMIZER","VALUE UNKNOWN")
       END SELECT
       CALL section_vals_val_get(ot_section,"SAFER_DIIS",l_val=qs_ot_env(1)%settings%safer_diis)
       CALL section_vals_val_get(ot_section,"LINESEARCH",i_val=ls_method)
       SELECT CASE(ls_method)
       CASE (ot_ls_none)
         qs_ot_env(1)%settings%line_search_method="NONE"
       CASE (ot_ls_2pnt)
         qs_ot_env(1)%settings%line_search_method="2PNT"
       CASE (ot_ls_3pnt)
         qs_ot_env(1)%settings%line_search_method="3PNT"
       CASE (ot_ls_gold)
         qs_ot_env(1)%settings%line_search_method="GOLD"
         CALL section_vals_val_get(ot_section,"GOLD_TARGET",r_val=qs_ot_env(1)%settings%gold_target)
       CASE DEFAULT
         CALL stop_program("READ OTSCF LS","VALUE UNKNOWN")
       END SELECT
       CALL section_vals_val_get(ot_section,"PRECONDITIONER",i_val=qs_ot_env(1)%settings%preconditioner_type)
       SELECT CASE(qs_ot_env(1)%settings%preconditioner_type)
       CASE(ot_precond_none)
           qs_ot_env(1)%settings%preconditioner_name="NONE"
       CASE(ot_precond_full_single)
           qs_ot_env(1)%settings%preconditioner_name="FULL_SINGLE"
       CASE(ot_precond_full_all)
           qs_ot_env(1)%settings%preconditioner_name="FULL_ALL"
       CASE(ot_precond_full_kinetic)
           qs_ot_env(1)%settings%preconditioner_name="FULL_KINETIC"
       CASE(ot_precond_s_inverse)
           qs_ot_env(1)%settings%preconditioner_name="FULL_S_INVERSE"
       CASE(ot_precond_sparse_diag)
           qs_ot_env(1)%settings%preconditioner_name="SPARSE_DIAG"
       CASE DEFAULT
         CALL stop_program("READ OTSCF PRECONDITIONER","VALUE UNKNOWN")
       END SELECT
       CALL section_vals_val_get(ot_section,"STEPSIZE",r_val=qs_ot_env(1)%settings%ds_min)
       CALL section_vals_val_get(ot_section,"ENERGY_GAP",r_val=qs_ot_env(1)%settings%energy_gap)
       CALL section_vals_val_get(ot_section,"EPS_TAYLOR",r_val=qs_ot_env(1)%settings%eps_taylor)
       CALL section_vals_val_get(ot_section,"MAX_TAYLOR",i_val=qs_ot_env(1)%settings%max_taylor)
       CALL section_vals_val_get(ot_section,"ROTATION",l_val=qs_ot_env(1)%settings%do_rotation)

    ! write OT output

    IF (globenv%ionode) THEN
    WRITE(globenv%scr,'(A)') " "
    WRITE(globenv%scr,'(A)') "  ----------------------------------- OT --------------------------------------"
    WRITE(globenv%scr,'(A)') " "
    WRITE(globenv%scr,'(A)') "  An efficient orbital transformation method "
    WRITE(globenv%scr,'(A)') "  for electronic structure calculations "
    WRITE(globenv%scr,'(A)') "  Joost VandeVondele and Juerg Hutter"
    WRITE(globenv%scr,'(A)') "  J. Chem. Phys. 118, 4365, 2003"
    WRITE(globenv%scr,'(A)') " "
    WRITE(globenv%scr,'(A,L2)') "  Allowing for rotations: ",qs_ot_env(1)%settings%do_rotation
    SELECT CASE (qs_ot_env(1)%settings%OT_METHOD)
    CASE ("SD")
         WRITE(globenv%scr,'(A)') "  minimizer      : SD             : steepest descent"
    CASE ("CG")
         WRITE(globenv%scr,'(A)') "  minimizer      : CG             : conjugate gradient"
    CASE ("DIIS")
         WRITE(globenv%scr,'(A)') "  minimizer      : DIIS           : direct inversion in the iterative subspace"
         WRITE(globenv%scr,'(A,I3,A)') "                            using : - ",qs_ot_env(1)%settings%diis_m," diis vectors"
         IF (qs_ot_env(1)%settings%safer_diis) THEN
            WRITE(globenv%scr,'(A,I3,A)') "                                    - safer DIIS on"
         ELSE
            WRITE(globenv%scr,'(A,I3,A)') "                                    - safer DIIS off"
         ENDIF
    CASE DEFAULT
         WRITE(globenv%scr,'(3A)') "  minimizer      :      ",qs_ot_env(1)%settings%OT_METHOD," : UNKNOWN"
    END SELECT
    SELECT CASE (qs_ot_env(1)%settings%preconditioner_name)
    CASE ("FULL_SINGLE")
         WRITE(globenv%scr,'(A)') "  preconditioner : FULL_SINGLE    : diagonalization based "
    CASE ("FULL_ALL")
         WRITE(globenv%scr,'(A)') "  preconditioner : FULL_ALL       : diagonalization, state selective "
    CASE ("FULL_KINETIC")
         WRITE(globenv%scr,'(A)') "  preconditioner : FULL_KINETIC   : cholesky inversion of T + eS"
    CASE ("FULL_S_INVERSE")
         WRITE(globenv%scr,'(A)') "  preconditioner : FULL_S_INVERSE : cholesky inversion of S"
    CASE ("SPARSE_DIAG")
         WRITE(globenv%scr,'(A)')  &
           "  preconditioner : SPARSE_DIAG    : diagonal atomic block diagonalization"
    CASE ("NONE")
         WRITE(globenv%scr,'(A)') "  preconditioner : NONE  "
    CASE DEFAULT 
         WRITE(globenv%scr,'(3A)') "  preconditioner : ",qs_ot_env(1)%settings%preconditioner_name," : UNKNOWN"
    END SELECT
    IF (qs_ot_env(1)%settings%OT_METHOD .EQ."SD".OR.qs_ot_env(1)%settings%OT_METHOD.EQ."CG") THEN
       SELECT CASE (qs_ot_env(1)%settings%line_search_method)
         CASE ("2PNT")
                    WRITE(globenv%scr,'(A)') "  line search    : 2PNT           : 2 energies, one gradient "
         CASE ("3PNT")
                    WRITE(globenv%scr,'(A)') "  line search    : 3PNT           : 3 energies "
         CASE ("GOLD")
                    WRITE(globenv%scr,'(A)') "  line search    : GOLD           : bracketing and golden section search "
                    WRITE(globenv%scr,'(A,F14.8)') "                   target relative accuracy : ", &
                               qs_ot_env(1)%settings%gold_target
         CASE ("NONE")
                    WRITE(globenv%scr,'(A)') "  line search    : NONE "
         CASE DEFAULT
                    WRITE(globenv%scr,'(3A)') "  line search : ", &
                               qs_ot_env(1)%settings%line_search_method," : UNKNOWN"
       END SELECT
       WRITE(globenv%scr,'(A,F14.8)') "  stepsize       :",qs_ot_env(1)%settings%ds_min
       WRITE(globenv%scr,'(A,F14.8)') "  energy_gap     :",qs_ot_env(1)%settings%energy_gap
       WRITE(globenv%scr,'(A,E14.5)') "  eps_taylor     :",qs_ot_env(1)%settings%eps_taylor
       WRITE(globenv%scr,'(A,I14)') "  max_taylor     :",qs_ot_env(1)%settings%max_taylor
       WRITE(globenv%scr,'(A)') "  "
       WRITE(globenv%scr,'(A)') "  ----------------------------------- OT --------------------------------------"
    ENDIF
 END IF

    ! copy the ot settings type so it is identical 
    nspin=SIZE(qs_ot_env)
    DO ispin=2,nspin
       qs_ot_env(ispin)%settings=qs_ot_env(1)%settings
    ENDDO


    CALL timestop(0.0_dp,handle)

  END SUBROUTINE ot_scf_read_input
!
! performs the actual minimisation, needs only limited info
! updated for restricted calculations
! matrix_dedc is the derivative of the energy with respect to the orbitals
! a null pointer for matrix_s implies that matrix_s is the unit matrix
!
!
SUBROUTINE ot_scf_mini(mo_array,matrix_dedc,matrix_s,energy, &
                                     energy_only,delta,qs_ot_env)
 

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dedc
    TYPE(real_matrix_type), POINTER          :: matrix_s
    REAL(KIND=dp)                            :: energy
    LOGICAL                                  :: energy_only
    REAL(KIND=dp)                            :: delta
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    CHARACTER(LEN=*), PARAMETER              :: RoutineN = "ot_scf_mini"

    INTEGER                                  :: handle, ispin, k, n, nspin
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset("ot_scf_mini","I","",handle)

    ! we use the size of the ot_env to indicate the spin. 
    ! this needs correction in case of restricted calculations
    nspin=SIZE(qs_ot_env)  

    IF ( .NOT. energy_only) THEN
       IF (qs_ot_env(1)%settings%do_rotation) THEN
          DO ispin=1,nspin
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k)
             CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,mo_coeff,matrix_dedc(ispin)%matrix, &
                             0.0_dp,qs_ot_env(ispin)%rot_mat_chcf)
          ENDDO
       ENDIF
    ENDIF

    ! notice we use (1) for driving all output and the minimization in case of LSD
    qs_ot_env(1)%etotal=energy
    CALL ot_mini(qs_ot_env,matrix_dedc)
    delta           =qs_ot_env(1)%delta
    energy_only     =qs_ot_env(1)%energy_only

    ! generate the orbitals using the new matrix_x
    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k)
       IF (ASSOCIATED(matrix_s)) THEN
         CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_x, &
                                       qs_ot_env(ispin)%matrix_sx, k)
       ELSE
         CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx)
       ENDIF
       CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin))
       CALL qs_ot_get_orbitals(mo_coeff,qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin))
    ENDDO
    IF (qs_ot_env(1)%restricted) THEN
       CALL mo_set_restrict(mo_array)
    ENDIF
    CALL timestop(0.0_dp,handle)

END SUBROUTINE ot_scf_mini
!
! initialises qs_ot_env so that mo_coeff is the current point
! and that the mimizization can be started.
!
SUBROUTINE ot_scf_init(mo_array,matrix_s,qs_ot_env)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: handle, ispin, k, n, nspin
    TYPE(cp_fm_type), POINTER                :: mo_coeff

! executable

    CALL timeset("ot_scf_init","I","",handle)

    ! adapted for work with the restricted keyword
    nspin=SIZE(qs_ot_env)

    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
       
       CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k)
   
       ! allocate
       CALL qs_ot_allocate(qs_ot_env(ispin),mo_coeff%matrix_struct)
   
       ! set c0,sc0
       CALL cp_fm_to_fm(mo_coeff,qs_ot_env(ispin)%matrix_c0)
       IF (ASSOCIATED(matrix_s)) THEN
          CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_c0, &
                                 qs_ot_env(ispin)%matrix_sc0,k)
       ELSE
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_c0,qs_ot_env(ispin)%matrix_sc0)
       ENDIF

       ! init
       CALL qs_ot_init(qs_ot_env(ispin))
   
       ! set x
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp)
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_sx,0.0_dp)

       IF (qs_ot_env(ispin)%settings%do_rotation) THEN
           CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp)
       ENDIF
   
       ! get c
       CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin))

    END DO
    CALL timestop(0.0_dp,handle)
END SUBROUTINE

SUBROUTINE ot_scf_destroy(qs_ot_env)

    TYPE(qs_ot_type)                         :: qs_ot_env

   CALL qs_ot_destroy(qs_ot_env)

END SUBROUTINE

END MODULE

