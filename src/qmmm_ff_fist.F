!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qmmm_ff_fist [1.0] *
!!
!!   NAME
!!     qmmm_ff_fist
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************


MODULE qmmm_ff_fist
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_describe_html, section_describe_index_html, section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "qmmm_ff_fist"

  PUBLIC :: qmmm_ff_precond_only_qm,&
            qmmm_coordinate_control

CONTAINS
  !!*************************************************************************
  !!****** qmmm_ff_fist/qmmm_coordinate_control [1.0] *
  !!
  !!   NAME
  !!     qmmm_coordinate_control
  !!
  !!   FUNCTION
  !!     Modifies the atom_info%label_atmname
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [tlaino]
  !!
  !***************************************************************************
  SUBROUTINE qmmm_coordinate_control(topology, globenv, error)
    IMPLICIT NONE
    
    TYPE(topology_parameters_type), &
         INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    LOGICAL :: explicit
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_coordinate_control', &
         routineP = moduleN//':'//routineN
    TYPE(section_vals_type), POINTER :: qmmm_section
    INTEGER                                  :: handle, iw, iatm,ikind,nkind
    INTEGER, DIMENSION(:), POINTER :: mm_indexes
    
    iw = globenv%scr
    CALL write_checkpoint_information("entering "//routineN,globenv)
    CALL timeset(routineN,'I','',handle)
    IF((globenv%ionode).AND.(globenv%print%level>0))THEN
       WRITE(iw,*) "  Entering qmmm_coordinate_control"
    END IF
    explicit = .FALSE.
    qmmm_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    qmmm_section => section_vals_get_subs_vals(qmmm_section,"QMMM",error=error)
    CALL section_vals_get(qmmm_section,explicit=explicit,error=error)
    IF (explicit) THEN
       qmmm_section => section_vals_get_subs_vals(qmmm_section,"QM_KIND",error=error)
       CALL section_vals_get(qmmm_section,n_repetition=nkind,error=error)
       DO ikind=1,nkind
          CALL section_vals_val_get(qmmm_section,"MM_INDEX",i_rep_section=ikind,i_vals=mm_indexes,&
               error=error)
          DO iatm = 1, SIZE(mm_indexes)
             topology%atom_info%label_atmname(mm_indexes(iatm)) = "_QM_"//&
                  TRIM(topology%atom_info%label_atmname(mm_indexes(iatm)))
             topology%atom_info%label_molname(mm_indexes(iatm)) = "_QM_"//&
                  TRIM(topology%atom_info%label_molname(mm_indexes(iatm)))
             topology%atom_info%label_resname(mm_indexes(iatm)) = "_QM_"//&
                  TRIM(topology%atom_info%label_resname(mm_indexes(iatm)))
          END DO
       END DO
    ELSE
       IF((globenv%ionode).AND.(globenv%print%level>0))THEN
          WRITE(iw,*) "  Exiting  qmmm_coordinate_control:: Skipping QMMM modifications! "
       END IF
    END IF

    IF((globenv%ionode).AND.(globenv%print%level>0))THEN
       WRITE(iw,*) "  Exiting  qmmm_coordinate_control"
    END IF
    
    CALL timestop(0.0_dp,handle)
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END SUBROUTINE qmmm_coordinate_control
  
  !!****f* qmmm_ff_fist/qmmm_ff_precond_only_qm [1.0] *
  !!
  !!   NAME
  !!     qmmm_ff_precond_only_qm
  !!
  !!   FUNCTION
  !!
  !!     This function handles the atom names and modifies
  !!     the "_QM_" prefix, in order to find the parameters
  !!     and return .TRUE. if all input labels belong to QM atoms
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [tlaino]
  !!
  !!*************************************************************************
  FUNCTION qmmm_ff_precond_only_qm(id1, id2, id3, id4) RESULT(only_qm)
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=default_string_length), INTENT(INOUT) :: id1
    CHARACTER(LEN=default_string_length), INTENT(INOUT), OPTIONAL :: id2, id3, id4
    LOGICAL :: only_qm
    ! Local Variables
    CHARACTER(LEN=default_string_length) :: tmp
    INTEGER :: j

    only_qm = .FALSE.
    IF (INDEX(id1,"_QM_") /= 0 ) THEN
       only_qm = .TRUE.
       tmp=TRIM(id1(5:))
       clean_string_1: DO j = 1, default_string_length
          id1(j:j) = " "
       END DO clean_string_1
       id1 = TRIM(tmp)
    END IF
    
    IF (PRESENT(id2)) THEN
       IF (INDEX(id2,"_QM_") /= 0 ) THEN
          tmp=TRIM(id2(5:))
          clean_string_2: DO j = 1, default_string_length
             id2(j:j) = " "
          END DO clean_string_2
          id2 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

    IF (PRESENT(id3)) THEN
       IF (INDEX(id3,"_QM_") /= 0 ) THEN
          tmp=TRIM(id3(5:))
          clean_string_3: DO j = 1, default_string_length
             id3(j:j) = " "
          END DO clean_string_3
          id3 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

    IF (PRESENT(id4)) THEN
       IF (INDEX(id4,"_QM_") /= 0 ) THEN
          tmp=TRIM(id4(5:))
          clean_string_4: DO j = 1, default_string_length
             id4(j:j) = " "
          END DO clean_string_4
          id4 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

   END FUNCTION qmmm_ff_precond_only_qm

END MODULE qmmm_ff_fist
