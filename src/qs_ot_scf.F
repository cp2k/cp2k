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
!!
!!   SOURCE
!******************************************************************************

MODULE qs_ot_scf

! *****************************************************************************

  USE kinds,                    ONLY: int_size, wp => dp
  USE global_types,             ONLY: global_environment_type
  USE termination,              ONLY: stop_memory,&
                                      stop_program
  USE timings,                  ONLY: timeset,&
                                      timestop

  USE qs_blacs,                 ONLY: cp_fm_set_all, &
                                      copy_blacs_to_blacs_matrix, &
                                      cp_sm_fm_multiply, &
                                      cp_fm_add
  USE cp_full_matrix,           ONLY: cp_full_matrix_type, cp_fm_get_info
  USE sparse_matrix_types,      ONLY: real_matrix_type
  USE qs_ot,                    ONLY: qs_ot_type,qs_ot_allocate, &
                                      qs_ot_destroy,qs_ot_init,qs_ot_get_p, &
                                      qs_ot_get_orbitals
  USE qs_ot_minimizer,               ONLY: ot_mini

  USE qs_parser,         ONLY: close_file,&
       open_file, finish_parser, read_object, start_parser,&
       test_object

  IMPLICIT NONE


! *** Public subroutines ***

  PUBLIC :: ot_scf_init
  PUBLIC :: ot_scf_mini
  PUBLIC :: ot_scf_destroy
  PUBLIC :: ot_scf_read_input

! *****************************************************************************

CONTAINS

SUBROUTINE ot_scf_read_input(qs_ot_env,globenv)
    TYPE(global_environment_type), INTENT(IN)       :: globenv
    TYPE(qs_ot_type)                                :: qs_ot_env

    ! locals 
    INTEGER :: handle
    CHARACTER(LEN=40) :: keyword
    CHARACTER(LEN=40) :: value
    LOGICAL :: found
    ! executable
    CALL timeset("ot_scf_read_input","I","",handle)
    ! decide default settings
    qs_ot_env%settings%ot_method="DIIS"
    qs_ot_env%settings%diis_m=7 
    qs_ot_env%settings%preconditioner_type="FULL_KINETIC"
    qs_ot_env%settings%line_search_method="2PNT"
    qs_ot_env%settings%ds_min=0.15_wp

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="OTSCF",section_found=found)
    if (found) then 
    DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("MINIMIZER")
               CALL read_object(value,lower_to_upper=.TRUE.)
               SELECT CASE (TRIM(value))
               CASE ("SD")
                    qs_ot_env%settings%ot_method="SD"
               CASE ("CG")
                    qs_ot_env%settings%ot_method="CG"
               CASE ("DIIS")
                    qs_ot_env%settings%ot_method="DIIS"
                    qs_ot_env%settings%diis_m=7 
                    IF (test_object().eq."INT") THEN
                       CALL read_object(qs_ot_env%settings%diis_m)
                    ENDIF
               CASE DEFAULT
                    IF (globenv%ionode) write(globenv%scr,*) value
                    CALL stop_program("READ OTSCF MINIMIZER","VALUE UNKNOWN")
               END SELECT
          CASE ("LINESEARCH")
               CALL read_object(value,lower_to_upper=.TRUE.)
               SELECT CASE (TRIM(value))
               CASE ("NONE")
                    qs_ot_env%settings%line_search_method="NONE"
               CASE ("2PNT")
                    qs_ot_env%settings%line_search_method="2PNT"
               CASE ("3PNT")
                    qs_ot_env%settings%line_search_method="3PNT"
               CASE ("GOLD")
                    qs_ot_env%settings%gold_target=0.01 ! 1 percent in the position 
                    qs_ot_env%settings%line_search_method="GOLD"
                    IF (test_object().eq."FLT") THEN
                       CALL read_object(qs_ot_env%settings%gold_target)
                    ENDIF
               CASE DEFAULT
                    IF (globenv%ionode) write(globenv%scr,*) value
                    CALL stop_program("READ OTSCF LINESEARCH","VALUE UNKNOWN")
               END SELECT
          CASE ("PRECONDITIONER")
               CALL read_object(value,lower_to_upper=.TRUE.)
               SELECT CASE (TRIM(value)) ! somehow in order of efficiency
               CASE ("FULL_SINGLE")
                    qs_ot_env%settings%preconditioner_type="FULL_SINGLE"
               CASE ("FULL_KINETIC")
                    qs_ot_env%settings%preconditioner_type="FULL_KINETIC"
               CASE ("FULL_S_INVERSE")
                    qs_ot_env%settings%preconditioner_type="FULL_S_INVERSE"
               CASE ("SPARSE_DIAG")
                    qs_ot_env%settings%preconditioner_type="SPARSE_DIAG"
               CASE ("NONE")
                    qs_ot_env%settings%preconditioner_type="NONE"
               CASE DEFAULT ! others are too experimental to be here
                    IF (globenv%ionode) write(globenv%scr,*) value
                    CALL stop_program("READ OTSCF PRECONDITIONER","VALUE UNKNOWN")
               END SELECT
          CASE ("STEPSIZE")
               CALL read_object(qs_ot_env%settings%ds_min)
          CASE DEFAULT
               IF (globenv%ionode) write(globenv%scr,*) keyword
               CALL stop_program("READ OTSCF","KEYWORD UNKNOWN")
          END SELECT 
    ENDDO
    endif

    CALL finish_parser()

    CALL timestop(0.0_wp,handle)

END SUBROUTINE
!
! performs the actual minimisation, needs only limited info
!
!
SUBROUTINE ot_scf_mini(mo_eigenvectors,matrix_ks,matrix_s,energy, &
                                     energy_only,delta,qs_ot_env,globenv)

    TYPE(cp_full_matrix_type), POINTER              :: mo_eigenvectors
    TYPE(real_matrix_type), POINTER                 :: matrix_ks,matrix_s
    TYPE(global_environment_type), INTENT(IN)       :: globenv
    TYPE(qs_ot_type)                                :: qs_ot_env
    LOGICAL                                         :: energy_only
    REAL(wp)                                        :: energy,delta

    IF ( .NOT. energy_only) THEN
       CALL cp_sm_fm_multiply(matrix_ks,mo_eigenvectors,&
                                    qs_ot_env%matrix_hc,&
                                   qs_ot_env%settings%k,&
                               para_env=globenv%para_env)
       CALL cp_fm_add(4.0_wp,qs_ot_env%matrix_hc,0.0_wp,qs_ot_env%matrix_hc)
    ENDIF

    qs_ot_env%etotal=energy
    CALL ot_mini(qs_ot_env,globenv)
    delta           =qs_ot_env%delta
    energy_only     =qs_ot_env%energy_only

    CALL cp_sm_fm_multiply(matrix_s,qs_ot_env%matrix_x, &
                                   qs_ot_env%matrix_sx, &
                                   qs_ot_env%settings%k,&
                                  para_env=globenv%para_env)
    CALL qs_ot_get_p(qs_ot_env%matrix_x,qs_ot_env%matrix_sx,qs_ot_env,globenv)
    CALL qs_ot_get_orbitals(mo_eigenvectors,qs_ot_env%matrix_x,qs_ot_env,globenv)

END SUBROUTINE
!
! initialises qs_ot_env so that mo_eigenvectors is the current point
! and that the mimizization can be started.
!
SUBROUTINE ot_scf_init(mo_eigenvectors,matrix_s,qs_ot_env,globenv)

    TYPE(cp_full_matrix_type), POINTER              :: mo_eigenvectors
    TYPE(real_matrix_type), POINTER                 :: matrix_s
    TYPE(global_environment_type), INTENT(IN)       :: globenv
    TYPE(qs_ot_type)                                :: qs_ot_env

    ! locals 
    INTEGER :: handle,nrow_block,ncol_block,n,k
    ! executable
    CALL timeset("ot_scf_init","I","",handle)

    CALL cp_fm_get_info(mo_eigenvectors,nrow_block=nrow_block, &
                                        ncol_block=ncol_block, &
                                                nrow_global=n, &
                                                  ncol_global=k)


    qs_ot_env%settings%n=n
    qs_ot_env%settings%k=k
    qs_ot_env%settings%ortho_k=k
    qs_ot_env%settings%nrow_block=nrow_block
    qs_ot_env%settings%ncol_block=ncol_block

    ! allocate
    CALL qs_ot_allocate(qs_ot_env,globenv)

    ! set c0,sc0
    CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,qs_ot_env%matrix_c0)
    CALL cp_sm_fm_multiply(matrix_s,qs_ot_env%matrix_c0, &
                           qs_ot_env%matrix_sc0,qs_ot_env%settings%k, &
                                             para_env=globenv%para_env)

    ! init
    CALL qs_ot_init(qs_ot_env,globenv)

    ! set x
    CALL cp_fm_set_all(qs_ot_env%matrix_x,0.0_wp)
    CALL cp_fm_set_all(qs_ot_env%matrix_sx,0.0_wp)

    ! get c
    CALL qs_ot_get_p(qs_ot_env%matrix_x,qs_ot_env%matrix_sx,qs_ot_env,globenv)

    CALL timestop(0.0_wp,handle)
END SUBROUTINE

SUBROUTINE ot_scf_destroy(qs_ot_env,globenv)

    TYPE(global_environment_type), INTENT(IN)       :: globenv
    TYPE(qs_ot_type)                                :: qs_ot_env

   CALL qs_ot_destroy(qs_ot_env,globenv)

END SUBROUTINE

END MODULE

