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
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm, &
                                             cp_fm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_ot,                           ONLY: qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_get_orbitals,&
                                             qs_ot_get_p,&
                                             qs_ot_init,&
                                             qs_ot_type
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE qs_parser,                       ONLY: close_file,&
                                             finish_parser,&
                                             open_file,&
                                             read_object,&
                                             start_parser,&
                                             test_object
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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
    TYPE(global_environment_type), INTENT(IN)       :: globenv
    TYPE(qs_ot_type), POINTER, DIMENSION(:)         :: qs_ot_env

    ! locals 
    INTEGER :: handle
    CHARACTER(LEN=40) :: keyword
    CHARACTER(LEN=40) :: value
    LOGICAL :: found
    INTEGER :: nspin,ispin
    ! executable
    CALL timeset("ot_scf_read_input","I","",handle)
    ! decide default settings

    nspin=SIZE(qs_ot_env)

    qs_ot_env(1)%settings%ot_method="CG"
    qs_ot_env(1)%settings%diis_m=7 
    qs_ot_env(1)%settings%preconditioner_type="FULL_KINETIC"
    qs_ot_env(1)%settings%line_search_method="2PNT"
    qs_ot_env(1)%settings%ds_min=0.15_dp

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="OTSCF",section_found=found)

    IF (found) THEN 
    DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("MINIMIZER")
               CALL read_object(value,lower_to_upper=.TRUE.)
               SELECT CASE (TRIM(value))
               CASE ("SD")
                    qs_ot_env(1)%settings%ot_method="SD"
               CASE ("CG")
                    qs_ot_env(1)%settings%ot_method="CG"
               CASE ("DIIS")
                    qs_ot_env(1)%settings%ot_method="DIIS"
                    qs_ot_env(1)%settings%diis_m=7 
                    IF (test_object().EQ."INT") THEN
                       CALL read_object(qs_ot_env(1)%settings%diis_m)
                    ENDIF
               CASE DEFAULT
                    IF (globenv%ionode) WRITE(globenv%scr,*) value
                    CALL stop_program("READ OTSCF MINIMIZER","VALUE UNKNOWN")
               END SELECT
          CASE ("LINESEARCH")
               CALL read_object(value,lower_to_upper=.TRUE.)
               SELECT CASE (TRIM(value))
               CASE ("NONE")
                    qs_ot_env(1)%settings%line_search_method="NONE"
               CASE ("2PNT")
                    qs_ot_env(1)%settings%line_search_method="2PNT"
               CASE ("3PNT")
                    qs_ot_env(1)%settings%line_search_method="3PNT"
               CASE ("GOLD")
                    qs_ot_env(1)%settings%gold_target=0.01 ! 1 percent in the position 
                    qs_ot_env(1)%settings%line_search_method="GOLD"
                    IF (test_object().EQ."FLT") THEN
                       CALL read_object(qs_ot_env(1)%settings%gold_target)
                    ENDIF
               CASE DEFAULT
                    IF (globenv%ionode) WRITE(globenv%scr,*) value
                    CALL stop_program("READ OTSCF LINESEARCH","VALUE UNKNOWN")
               END SELECT
          CASE ("PRECONDITIONER")
               CALL read_object(value,lower_to_upper=.TRUE.)
               SELECT CASE (TRIM(value)) ! somehow in order of efficiency
               CASE ("FULL_SINGLE")
                    qs_ot_env(1)%settings%preconditioner_type="FULL_SINGLE"
               CASE ("FULL_KINETIC")
                    qs_ot_env(1)%settings%preconditioner_type="FULL_KINETIC"
               CASE ("FULL_S_INVERSE")
                    qs_ot_env(1)%settings%preconditioner_type="FULL_S_INVERSE"
               CASE ("SPARSE_DIAG")
                    qs_ot_env(1)%settings%preconditioner_type="SPARSE_DIAG"
               CASE ("NONE")
                    qs_ot_env(1)%settings%preconditioner_type="NONE"
               CASE DEFAULT ! others are too experimental to be here
                    IF (globenv%ionode) WRITE(globenv%scr,*) value
                    CALL stop_program("READ OTSCF PRECONDITIONER","VALUE UNKNOWN")
               END SELECT
          CASE ("STEPSIZE")
               CALL read_object(qs_ot_env(1)%settings%ds_min)
          CASE DEFAULT
               IF (globenv%ionode) WRITE(globenv%scr,*) keyword
               CALL stop_program("READ OTSCF","KEYWORD UNKNOWN")
          END SELECT 
    ENDDO
    ENDIF
    CALL finish_parser()
    IF (globenv%ionode) THEN
    WRITE(6,*) ""
    WRITE(6,*) " ----------------------------------- OT --------------------------------------"
    WRITE(6,*) ""
    WRITE(6,*) " An efficient orbital transformation method "
    WRITE(6,*) " for electronic structure calculations "
    WRITE(6,*) " Joost VandeVondele and Juerg Hutter"
    WRITE(6,*) " J. Chem. Phys. 118, 4365, 2003"
    WRITE(6,*) ""
    SELECT CASE (qs_ot_env(1)%settings%OT_METHOD)
    CASE ("SD")
         WRITE(6,*) " minimizer      : SD             : steepest descent"
    CASE ("CG")
         WRITE(6,*) " minimizer      : CG             : conjugate gradient"
    CASE ("DIIS")
         WRITE(6,*) " minimizer      : DIIS           : direct inversion in the iterative subspace"
         WRITE(6,*) "                           using : ",qs_ot_env(1)%settings%diis_m," diis vectors"
    CASE DEFAULT
         WRITE(6,*) " minimizer      :      ",qs_ot_env(1)%settings%OT_METHOD," : UNKNOWN"
    END SELECT
    SELECT CASE (qs_ot_env(1)%settings%preconditioner_type)
    CASE ("FULL_SINGLE")
         WRITE(6,*) " preconditioner : FULL_SINGLE    : diagonalization based "
    CASE ("FULL_KINETIC")
         WRITE(6,*) " preconditioner : FULL_KINETIC   : cholesky inversion of T + eS"
    CASE ("FULL_S_INVERSE")
         WRITE(6,*) " preconditioner : FULL_S_INVERSE : cholesky inversion of S"
    CASE ("SPARSE_DIAG")
         WRITE(6,*)  &
           " preconditioner : SPARSE_DIAG    : diagonal atomic block diagonalization"
    CASE ("NONE")
         WRITE(6,*) " preconditioner : NONE  "
    CASE DEFAULT 
         WRITE(6,*) " preconditioner : ",qs_ot_env(1)%settings%preconditioner_type," : UNKNOWN"
    END SELECT
    IF (qs_ot_env(1)%settings%OT_METHOD .EQ."SD".OR.qs_ot_env(1)%settings%OT_METHOD.EQ."CG") THEN
       SELECT CASE (qs_ot_env(1)%settings%line_search_method)
         CASE ("2PNT")
                    WRITE(6,*) " line search    : 2PNT           : 2 energies, one gradient "
         CASE ("3PNT")
                    WRITE(6,*) " line search    : 3PNT           : 3 energies "
         CASE ("GOLD")
                    WRITE(6,*) " line search    : GOLD           : bracketing and golden section search "
                    WRITE(6,*) "                  target relative accuracy : ", &
                               qs_ot_env(1)%settings%gold_target
         CASE ("NONE")
                    WRITE(6,*) " line search    : NONE "
         CASE DEFAULT
                    WRITE(6,*) " line search : ", &
                               qs_ot_env(1)%settings%line_search_method," : UNKNOWN"
       END SELECT
    ENDIF
    WRITE(6,*) " stepsize       :",qs_ot_env(1)%settings%ds_min
    WRITE(6,*) " "
    WRITE(6,*) " ----------------------------------- OT --------------------------------------"
    ENDIF
 
    DO ispin=2,nspin
       qs_ot_env(ispin)%settings=qs_ot_env(1)%settings
    ENDDO


    CALL timestop(0.0_dp,handle)

END SUBROUTINE
!
! performs the actual minimisation, needs only limited info
!
!
SUBROUTINE ot_scf_mini(mo_array,matrix_ks,matrix_s,energy, &
                                     energy_only,delta,qs_ot_env,globenv)
 

    TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mo_array
    TYPE(real_matrix_type), POINTER                 :: matrix_s
    TYPE(real_matrix_p_type), POINTER, DIMENSION(:) :: matrix_ks
    TYPE(global_environment_type), INTENT(IN)       :: globenv
    TYPE(qs_ot_type), POINTER, DIMENSION(:)         :: qs_ot_env
    LOGICAL                                         :: energy_only
    REAL(KIND = dp)                                        :: energy,delta

    INTEGER nspin,ispin,handle
    TYPE(cp_fm_type), POINTER        :: mo_eigenvectors

    CALL timeset("ot_scf_mini","I","",handle)

    nspin=SIZE(mo_array)

    IF ( .NOT. energy_only) THEN
       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,eigenvectors=mo_eigenvectors)
          CALL cp_sm_fm_multiply(matrix_ks(ispin)%matrix,mo_eigenvectors,&
                                       qs_ot_env(ispin)%matrix_hc,&
                                       qs_ot_env(ispin)%settings%k,&
                                       para_env=globenv%para_env)
          SELECT CASE (nspin)
          CASE (1)
               CALL cp_fm_scale_and_add(4.0_dp,qs_ot_env(ispin)%matrix_hc)
          CASE (2)            
               CALL cp_fm_scale_and_add(2.0_dp,qs_ot_env(ispin)%matrix_hc)
          CASE DEFAULT 
               CALL stop_program("OT_SCF_MINI","NSPIN NYI")
          END SELECT 
       ENDDO
    ENDIF

    ! notice we use (1) for driving all output and the minimization
    qs_ot_env(1)%etotal=energy
    CALL ot_mini(qs_ot_env)
    delta           =qs_ot_env(1)%delta
    energy_only     =qs_ot_env(1)%energy_only

    DO ispin=1,nspin
       CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_x, &
                                      qs_ot_env(ispin)%matrix_sx, &
                                      qs_ot_env(ispin)%settings%k,&
                                      para_env=globenv%para_env)
       CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin))
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,eigenvectors=mo_eigenvectors)
       CALL qs_ot_get_orbitals(mo_eigenvectors,qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin))
    ENDDO

    CALL timestop(0.0_dp,handle)

END SUBROUTINE
!
! initialises qs_ot_env so that mo_eigenvectors is the current point
! and that the mimizization can be started.
!
SUBROUTINE ot_scf_init(mo_array,matrix_s,qs_ot_env)

    TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mo_array
    TYPE(real_matrix_type), POINTER                 :: matrix_s
    TYPE(qs_ot_type), POINTER, DIMENSION(:)         :: qs_ot_env

    ! locals 
    INTEGER :: handle,nrow_block,ncol_block,n,k
    TYPE(cp_fm_type), POINTER              :: mo_eigenvectors
    INTEGER nspin,ispin
    ! executable
    CALL timeset("ot_scf_init","I","",handle)
    nspin=SIZE(mo_array)

    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,eigenvectors=mo_eigenvectors)
       CALL cp_fm_get_info(mo_eigenvectors,nrow_block=nrow_block, &
                                           ncol_block=ncol_block, &
                                                   nrow_global=n, &
                                                     ncol_global=k)
   

       qs_ot_env(ispin)%settings%n=n
       qs_ot_env(ispin)%settings%k=k
       qs_ot_env(ispin)%settings%ortho_k=k
       qs_ot_env(ispin)%settings%nrow_block=nrow_block
       qs_ot_env(ispin)%settings%ncol_block=ncol_block

       ! allocate
       CALL qs_ot_allocate(qs_ot_env(ispin),mo_eigenvectors%matrix_struct%para_env, &
                                            mo_eigenvectors%matrix_struct%context)
   
       ! set c0,sc0
       CALL cp_fm_to_fm(mo_eigenvectors,qs_ot_env(ispin)%matrix_c0)
       CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_c0, &
                              qs_ot_env(ispin)%matrix_sc0,qs_ot_env(ispin)%settings%k, &
                                                para_env=mo_eigenvectors%matrix_struct%para_env)

       ! init
       CALL qs_ot_init(qs_ot_env(ispin))
   
       ! set x
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp)
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_sx,0.0_dp)
   
       ! get c
       CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin))

    END DO
    CALL timestop(0.0_dp,handle)
END SUBROUTINE

SUBROUTINE ot_scf_destroy(qs_ot_env)

    TYPE(qs_ot_type)                                :: qs_ot_env

   CALL qs_ot_destroy(qs_ot_env)

END SUBROUTINE

END MODULE

