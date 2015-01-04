!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - taken out of input_cp2k_motion
!> \author Ole Schuett
! *****************************************************************************
MODULE input_cp2k_motion_print
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level,&
                                             silent_print_level
  USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_dcd_aligned_cell,&
                                             dump_pdb,&
                                             dump_xmol
  USE input_cp2k_subsys,               ONLY: create_structure_data_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_motion_print'

PUBLIC :: create_motion_print_section, add_format_keyword

CONTAINS

! *****************************************************************************
!> \brief creates the motion%print section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_motion_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_motion_print_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    NULLIFY(keyword, section, print_key)

    CALL section_create(section,name="print",&
         description="Controls the printing properties during an MD run",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
         error=error)

    CALL cp_print_key_section_create(print_key,"TRAJECTORY",&
         description="Controls the output of the trajectory",&
         print_level=low_print_level, common_iter_levels=1,&
         filename="",unit_str="angstrom",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
         description="Specifies the format of the output file for the trajectory.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"SHELL_TRAJECTORY",&
         description="Controls the output of the trajectory of shells when the shell-model is used ",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="angstrom",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
         description="Specifies the format of the output file for the trajectory of shells.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"CORE_TRAJECTORY",&
         description="Controls the output of the trajectory of cores when the shell-model is used ",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="angstrom",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
         description="Specifies the format of the output file for the trajectory of cores.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"CELL",&
         description="Controls the output of the simulation cell. "//&
         "For later analysis of the trajectory it is recommendable that the "//&
         "frequency of printing is the same as the one used for the trajectory file.",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"VELOCITIES",&
         description="Controls the output of the velocities",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="bohr*au_t^-1",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the velocities.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"SHELL_VELOCITIES",&
         description="Controls the output of the velocities of shells when the shell model is used",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="bohr*au_t^-1",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the velocities of shells.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"CORE_VELOCITIES",&
         description="controls the output of the velocities of cores when the shell model is used",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="bohr*au_t^-1",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the velocities of cores.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL create_structure_data_section(print_key, error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"FORCE_MIXING_LABELS",&
         description="Controls the output of the force mixing (FORCE_EVAL&QMMM&FORCE_MIXING) labels",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the force mixing labels.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"FORCES",&
         description="Controls the output of the forces",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="hartree*bohr^-1",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the forces.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"SHELL_FORCES",&
         description="controls the output of the forces on shells when shell-model is used",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="hartree*bohr^-1",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the forces on shells.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"CORE_FORCES",&
         description="controls the output of the forces on cores when shell-model is used",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",unit_str="hartree*bohr^-1",error=error)
    CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
         description="Specifies the format of the output file for the forces on cores.", error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"MIXED_ENERGIES",&
         description="Controls the output of the energies of the two"//&
         "regular FORCE_EVALS in the MIXED method"//&
         "printed is step,time,Etot,E_F1,E_F2,CONS_QNT",&
         print_level=low_print_level, common_iter_levels=1,&
         filename="",error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"STRESS",&
         description="Controls the output of the stress tensor",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"RESTART",&
         description="Controls the dumping of the restart file during runs. "//&
         "By default keeps a short history of three restarts. See also RESTART_HISTORY", &
         each_iter_names=s2a("MD"),each_iter_values=(/20/), &
         print_level=silent_print_level, common_iter_levels=1,  &
         add_last=add_last_numeric, filename="",error=error)

    CALL keyword_create(keyword, name="BACKUP_COPIES",&
         description="Specifies the maximum number of backup copies.",&
         usage="BACKUP_COPIES {int}",&
         default_i_val=3, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SPLIT_RESTART_FILE",&
                        description="If specified selected input sections, which are growing with the "//&
                                    "number of atoms in the system, are written to another restart file "//&
                                    "in binary format instead of the default restart file in human "//&
                                    "readable ASCII format. This split of the restart file may "//&
                                    "provide significant memory savings and an accelerated I/O for "//&
                                    "systems with a very large number of atoms",&
                        usage="SPLIT_RESTART_FILE yes",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"RESTART_HISTORY",&
         description="Dumps unique restart files during the run keeping all of them."//&
         "Most useful if recovery is needed at a later point.",&
         print_level=low_print_level, common_iter_levels=0, &
         each_iter_names=s2a("MD","GEO_OPT","ROT_OPT"),each_iter_values=(/500,500,500/), &
         filename="",error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)


    CALL cp_print_key_section_create(print_key,"TRANSLATION_VECTOR",&
         description="Dumps the translation vector applied along an MD (if any). Useful"//&
         " for postprocessing of QMMM trajectories in which the QM fragment is continuously"//&
         " centered in the QM box",&
         print_level=high_print_level, common_iter_levels=1,&
         filename="",error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

  END SUBROUTINE create_motion_print_section


! *****************************************************************************
!> \brief creates the FORMAT keyword
!> \param keyword ...
!> \param section will contain the pint section
!> \param pos ...
!> \param description ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino 10.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE add_format_keyword(keyword, section, pos, description, error)
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN)                      :: pos
    CHARACTER(LEN=*), INTENT(IN)             :: description
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_format_keyword', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)

    IF (pos) THEN

       CALL keyword_create(keyword, name="FORMAT",&
            description=description,usage="FORMAT (ATOMIC|DCD|PDB|XMOL|XYZ)",&
            default_i_val=dump_xmol,&
            enum_c_vals=s2a("ATOMIC","DCD","DCD_ALIGNED_CELL","PDB","XMOL","XYZ"),&
            enum_i_vals= (/dump_atomic,dump_dcd,dump_dcd_aligned_cell,dump_pdb,dump_xmol,dump_xmol/),&
            enum_desc=s2a("Write only the coordinates X,Y,Z without element symbols to a formatted file",&
                          "Write the coordinates (no element labels) and the cell information to a binary file",&
                          "Like DCD, but the dumped coordinates refer to an aligned cell following the common convention: "//&
                          "the cell vector <b>a</b> is aligned with the <i>x</i> axis and the cell vector <b>b</b> lies in "//&
                          "the <i>xy</i> plane. This allows the reconstruction of scaled coordinates from the DCD data only.",&
                          "Write the atomic information in PDB format to a formatted file",&
                          "Mostly known as XYZ format, provides in a formatted file: element_symbol X Y Z",&
                          "Alias name for XMOL"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_OCCUP",&
            variants=(/"CHARGE_O"/),&
            description="Write the MM charges to the OCCUP field of the PDB file",&
            usage="CHARGE_OCCUP logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_BETA",&
            variants=(/"CHARGE_B"/),&
            description="Write the MM charges to the BETA field of the PDB file",&
            usage="CHARGE_BETA logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_EXTENDED",&
            description="Write the MM charges to the very last field of the PDB file (starting from column 81)",&
            usage="CHARGE_EXTENDED logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    ELSE

       CALL keyword_create(keyword, name="FORMAT",&
            description=description,usage="FORMAT (ATOMIC|DCD|XMOL|XYZ)",&
            default_i_val=dump_xmol,&
            enum_c_vals=s2a("ATOMIC","DCD","XMOL","XYZ"),&
            enum_i_vals= (/dump_atomic,dump_dcd,dump_xmol,dump_xmol/),&
            enum_desc=s2a("Write only the coordinates X,Y,Z without element symbols to a formatted file",&
                          "Write the coordinates (no element labels) and the cell information to a binary file",&
                          "Mostly known as XYZ format, provides in a formatted file: element_symbol X Y Z",&
                          "Alias name for XMOL"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE add_format_keyword

END MODULE input_cp2k_motion_print
