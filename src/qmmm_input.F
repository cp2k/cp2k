!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_iput [1.0] *
!!
!!   NAME
!!     qmmm_input
!!
!!   FUNCTION
!!     Set the QM/MM environment defined in Namelist QMMM
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     6.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_input
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             read_atomic_kind_set, &
                                             get_atomic_kind,&
                                             set_atomic_kind
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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dbl
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_create,&
                                             qmmm_env_type
  USE qs_environment,                  ONLY: qs_init_subsys
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_retain,&
                                             cell_type,&
                                             init_qmmm_cell,&
                                             cell_clone
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: connectivity_control,&
                                             force_field_control,&
                                             coordinate_control
  USE topology_pdb,                    ONLY: write_coordinate_pdb
  USE topology_psf,                    ONLY: write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_constraint_pack,&
                                             topology_coordinate_pack
  USE qs_parser,                       ONLY: finish_parser,&
                                             open_file,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             allocate_potential,&
                                             set_potential
  USE cp_keywords,                     ONLY: keyword_type, &
                                             list_keywords, &
                                             find_keyword

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_input'

  PUBLIC :: qmmm_namelist_read
  PRIVATE :: init_keywords
  INTEGER, PARAMETER, PRIVATE :: key_qm_atoms       =1
  INTEGER, PARAMETER, PRIVATE :: key_a_qm           =2
  INTEGER, PARAMETER, PRIVATE :: key_b_qm           =3
  INTEGER, PARAMETER, PRIVATE :: key_c_qm           =4
  INTEGER, PARAMETER, PRIVATE :: key_abc_qm         =5
  INTEGER, PARAMETER, PRIVATE :: key_coupling_type  =6
  INTEGER, PARAMETER, PRIVATE :: key_radius         =7
  INTEGER, PARAMETER, PRIVATE :: key_eps_mm_rspace  =8
  INTEGER, PARAMETER, PRIVATE :: Nkeywords          =8 ! To be updated when New Keywords added !
  TYPE(keyword_type), PRIVATE, DIMENSION(1:Nkeywords) :: keywords 
!***
!****************************************************************************
CONTAINS
  !!****f* qmmm_input/init_keywords [1.0] *
  !!
  !!   NAME
  !!     init_keywords
  !!
  !!   FUNCTION
  !!
  !!     init the keywords of the qmmm section
  !!     define the constants as module parameters, 
  !!     and increase the lenght before adding here
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
  !!     06.2004 created [tlaino]
  !!
  !!*************************************************************************

  SUBROUTINE init_keywords()
    keywords(key_QM_ATOMS)=keyword_type("QM_ATOMS", &
         "QM_ATOMS 3 \n H 2 \n  1 2 \n O 1 \n  3",&
         "Set and specifies the number, the typo and the index of the QM Atoms",&
         "QMATOMS","","")
    keywords(key_A_QM)=keyword_type("A_QM",&
         "A_QM Real Real Real",&
         "Set the A vector of the QM Cell to the Reals value. The other values will be"//&
         "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
         "","","")
    keywords(key_B_QM)=keyword_type("B_QM",&
         "B_QM Real Real Real",&
         "Set the B vector of the QM Cell to the Reals value. The other values will be"//&
         "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
         "","","")
    keywords(key_C_QM)=keyword_type("C_QM",&
         "C_QM Real Real Real",&
         "Set the C vector of the QM Cell to the Reals value. The other values will be"//&
         "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
         "","","")
    keywords(key_ABC_QM)=keyword_type("ABC_QM",&
         "ABC_QM Real Real Real",&
         "Set the A,B,C vectors of  QM Cell to the Reals value."//&
         "It changes the overall shape of the QM cell.",&
         "","","")
    keywords(key_COUPLING_TYPE)=keyword_type("E_COUPL",&
         "E_COUPL [NONE,GRID,GAUSS]",&
         "Set the Coupling's type between QM and MM system. Possible Options:"//&
         "(1) NONE: Just Mechanical Coupling. (2) GRID: Computing the electrostatic"//&
         "potential analytically on the Grid. (3) GAUSS: Mapping the MM electrostatic"//&
         "potential on the Grid with Gaussian.",&
         "ECOUPL","","")
    keywords(key_RADIUS)=keyword_type("RADIUS",&
         "RADIUS \n O 0.7 \n H 0.4 \n END",&
         "Set the charge radius to evaluate the MM electrostatic potential in the QM/MM"//&
         "coupling.",&
         "","","")
    keywords(key_EPS_MM_RSPACE)=keyword_type("EPS_MM_RSPACE",&
         "EPS_MM_RSPACE  Real",&
         "Set the threshold g(r) - threshold = 0. Where g(r) is a gaussian function."//&
         "Standard value set to 1.0D-10",&
         "","","")
  END SUBROUTINE init_keywords

  !!****f* qmmm_input/qmmm_namelist_read [1.0] *
  !!
  !!   NAME
  !!     qmmm_namelist_read
  !!
  !!   FUNCTION
  !!     Read "QMMM" Namelist Arguments
  !!     Arguments:
  !!     &QMMM
  !!        QM_ATOMS Number_of_QM_atoms
  !!        Atom_LABEL Number_of_this_atomic_Kind
  !!        Index1 Index2 ... IndexN
  !!        Atom_LABEL Number_of_this_atomic_Kind
  !!        Index1 Index2 ... IndexL
  !!   
  !!        ABC_QM x y z
  !!        ECOUPL (E_COUPL) [GRID] [GAUSS]
  !!        EPS_MM_RSPACE Real
  !!               
  !!     &END
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
  !!     05.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE qmmm_namelist_read(globenv,&
                                qmmm_env,&
                                subsys_mm,&
                                qm_atom_type,&
                                qm_atom_index,&
                                mm_atom_index,&
                                qm_subsys_index,&
                                qm_cell_small,&
                                qmmm_coupl_type,&
                                eps_mm_rspace,&
                                error)
    IMPLICIT NONE
    ! Arguments
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_atom_index, qm_subsys_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    CHARACTER(len=default_string_length), INTENT(OUT) :: qmmm_coupl_type
    REAL(KIND=dbl), INTENT(OUT)              :: eps_mm_rspace
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_namelist_read',&
         routineP = moduleN//':'//routineN
    INTEGER :: i, j, icount, size_mm_system, stat, ikind
    LOGICAL :: is_mm
    INTEGER                                  :: num_qm_atoms, qm_num_atom_kind, qm_num_atom_tot
    LOGICAL                                  :: failure, section_found
    CHARACTER(len=default_string_length)     :: keyword, qm_atom_kind
    CHARACTER(len=3)                         :: test_result
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type),POINTER           :: atomic_kind
    CHARACTER(LEN=default_string_length)     :: atmname, tmp_kind
    REAL(dbl)                                :: tmp_radius
    TYPE(elp_potential_type), POINTER        :: elp_potential
    

    qmmm_coupl_type = "NONE"
    eps_mm_rspace   = 1.0E-10_dbl
    !
    ! Initialize core_radius_charge
    !
    tmp_radius=0.0_dbl
    CALL cp_subsys_get(subsys=subsys_mm(1)%subsys,&
                       atomic_kinds=atomic_kinds,&
                       error=error)
    DO Ikind = 1, SIZE(atomic_kinds%els)
       atomic_kind => atomic_kinds%els(Ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            elp_potential=elp_potential)
       CALL set_potential  (potential=elp_potential,&
                            core_charge_radius=tmp_radius)
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
                            elp_potential=elp_potential)
    END DO
    !
    ! Check for &QMMM Namelist
    !
    CALL init_keywords()
    CALL start_parser(file_name=globenv%input_file_name,&
         globenv=globenv,&
         start_section_label="QMMM",&
         section_found=section_found)
    
    IF (section_found) THEN
       !
       ! Read Info from QMMM Namelist
       !
       DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          test_result = test_object(newline=.TRUE.)
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (find_keyword(keywords,keyword))
          CASE(key_QM_ATOMS)
             qm_num_atom_tot = 0
             CALL read_object(num_qm_atoms)
             ALLOCATE(qm_atom_index(num_qm_atoms),& 
                  qm_atom_type(num_qm_atoms),&
                  qm_subsys_index(num_qm_atoms),&
                  stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ! Read the QM Atomic Indexes..
             DO WHILE (qm_num_atom_tot /= num_qm_atoms)
                CALL read_object(qm_atom_kind,newline=.TRUE.)
                CALL read_object(qm_num_atom_kind)
                IF (qm_num_atom_kind /= 0) THEN
                   i = qm_num_atom_tot + 1
                   CALL read_object(qm_atom_index(i),newline=.TRUE.)
                   qm_atom_type(i) = qm_atom_kind
                   DO i = qm_num_atom_tot + 2, qm_num_atom_kind + qm_num_atom_tot
                      CALL read_object(qm_atom_index(i))
                      qm_atom_type(i) = qm_atom_kind
                   END DO
                END IF
                qm_num_atom_tot = qm_num_atom_tot + qm_num_atom_kind
             END DO
             qm_subsys_index=1
          CASE (key_A_QM)
             CALL read_object(qm_cell_small%hmat(1,1))
             CALL read_object(qm_cell_small%hmat(2,1))
             CALL read_object(qm_cell_small%hmat(3,1))
             qm_cell_small%hmat=qm_cell_small%hmat(:,:)*qm_cell_small%unit_of_length
          CASE (key_ABC_QM)
             CALL read_object(qm_cell_small%hmat(1,1))
             CALL read_object(qm_cell_small%hmat(2,2))
             CALL read_object(qm_cell_small%hmat(3,3))
             qm_cell_small%hmat=qm_cell_small%hmat(:,:)*qm_cell_small%unit_of_length
          CASE (key_B_QM)
             CALL read_object(qm_cell_small%hmat(1,2))
             CALL read_object(qm_cell_small%hmat(2,2))
             CALL read_object(qm_cell_small%hmat(3,2))
             qm_cell_small%hmat=qm_cell_small%hmat(:,:)*qm_cell_small%unit_of_length
          CASE (key_C_QM)
             CALL read_object(qm_cell_small%hmat(1,3))
             CALL read_object(qm_cell_small%hmat(2,3))
             CALL read_object(qm_cell_small%hmat(3,3))
             qm_cell_small%hmat=qm_cell_small%hmat(:,:)*qm_cell_small%unit_of_length
          CASE(key_COUPLING_TYPE)
             CALL read_object(qmmm_coupl_type)
          CASE(key_RADIUS)
             DO
                CALL read_object(tmp_kind, lower_to_upper=.TRUE.,newline=.TRUE.)
                IF(tmp_kind=="END") EXIT
                WRITE(*,*) "Kind: ",tmp_kind
                CALL read_object(tmp_radius)
                tmp_radius = tmp_radius*qm_cell_small%unit_of_length
                WRITE(*,*) "RRRADIUS in atomic units",tmp_kind,"  ",tmp_radius
                !
                ! Modify the MM forcefield to store the information on the QM/MM 
                ! Electrostatic radius potential
                !
                DO Ikind = 1, SIZE(atomic_kinds%els)
                   atomic_kind => atomic_kinds%els(Ikind)
                   CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                        name=atmname)
                   IF (TRIM(tmp_kind)==atmname) THEN
                      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                           elp_potential=elp_potential)
                      CALL set_potential(potential=elp_potential,&
                                         core_charge_radius=tmp_radius)
                      CALL set_atomic_kind(atomic_kind=atomic_kind,&
                                           elp_potential=elp_potential)
                   END IF
                END DO
             END DO
          CASE(key_EPS_MM_RSPACE)
             CALL read_object(eps_mm_rspace)
          CASE DEFAULT
             CALL list_keywords(keywords,globenv%scr)
             CALL stop_parser(routine=routineN,message="INVALID_KEYWORD "&
                  //keyword)
          END SELECT
       END DO
    ELSE
       !
       ! Error: QM/MM calculation requested but no QMMM section found
       ! 
       CALL stop_parser(routineN,moduleN,__LINE__,"QMMM SECTION_NOT_FOUND")
    ENDIF
    ! Initialize cell parameters..
    CALL init_qmmm_cell(qm_cell_small)
    ! Build MM atoms list
    size_mm_system = SIZE(subsys_mm(1)%subsys%particles%els) -&
         SIZE(qm_atom_index)
    ALLOCATE(mm_atom_index(size_mm_system))
    icount = 0
    DO i = 1, SIZE(subsys_mm(1)%subsys%particles%els)
       is_mm = .TRUE.
       DO j = 1, SIZE(qm_atom_index)
          IF (i.EQ.qm_atom_index(j)) THEN
             is_mm = .FALSE.
          END IF
       END DO
       IF (is_mm) THEN
          icount = icount + 1
          mm_atom_index(icount) = i
       END IF
    END DO
    !
    CPPostcondition(icount==size_mm_system,cp_failure_level,routineP,error,failure)
    ! Terminate parser..
    CALL finish_parser()
  END SUBROUTINE qmmm_namelist_read

END MODULE qmmm_input

