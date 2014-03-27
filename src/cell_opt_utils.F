!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains a functional that calculates the energy and its derivatives
!>      for the geometry optimizer
!> \par History
!>      03.2008 - Teodoro Laino [tlaino] - University of Zurich - Cell Optimization
! *****************************************************************************
MODULE cell_opt_utils
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell_param,&
                                             set_cell_param
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_units_rad
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set_cell,&
                                             force_env_type
  USE input_constants,                 ONLY: &
       cell_sym_cubic, cell_sym_hexagonal, cell_sym_monoclinic, &
       cell_sym_orthorhombic, cell_sym_rhombohedral, cell_sym_tetragonal_ab, &
       cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, cell_sym_triclinic
  USE input_cp2k_global,               ONLY: create_global_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: sqrt3
  USE mathlib,                         ONLY: angle
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_utils'

  PUBLIC :: get_ut_cell_matrix, get_dg_dh, gopt_new_error_create,&
            gopt_new_error_release, read_external_press_tensor

CONTAINS
! *****************************************************************************
!> \param Transform a general cell matrix into an upper triangular one
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE get_ut_cell_matrix(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_ut_cell_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(3)                    :: periodic
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: cell_angle, cell_length
    TYPE(cell_type), POINTER                 :: cell

    failure  = .FALSE.
    CALL force_env_get(force_env, cell=cell, error=error)
    ! orthorhombic cells are diagonal, and should stay exactly like this
    IF (.NOT. cell%orthorhombic) THEN
      CALL get_cell_param(cell,cell_length,cell_angle,cp_units_rad,periodic=periodic,&
                          error=error)
      CALL set_cell_param(cell,cell_length,cell_angle,periodic=periodic,&
                          do_init_cell=.TRUE.,error=error)
    END IF
    CALL force_env_set_cell(force_env, cell=cell, error=error)

  END SUBROUTINE get_ut_cell_matrix

! *****************************************************************************
!> \param creates a new error used for cell optimization algorithm
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE gopt_new_error_create(new_error, root_section, para_env, project_name,&
       id_run, error)
    TYPE(cp_error_type), INTENT(out)         :: new_error
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(len=default_string_length), &
      INTENT(OUT)                            :: project_name
    INTEGER, INTENT(IN)                      :: id_run
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gopt_new_error_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: c_val, input_file_path, &
                                                output_file_path
    CHARACTER(len=default_string_length)     :: label
    INTEGER                                  :: i, lp, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: loc_logger, logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section

    failure = .FALSE.
    NULLIFY(loc_logger, logger, enum, keyword, section)
    logger => cp_error_get_logger(error)

    CALL create_global_section(section,error=error)
    keyword => section_get_keyword(section,"RUN_TYPE",error=error)
    CALL keyword_get(keyword,enum=enum,error=error)
    label = TRIM(enum_i2c(enum,id_run,error=error))
    CALL section_release(section, error)

    ! Redirecting output of the sub_calculation to a different file
    CALL section_vals_val_get(root_section,"GLOBAL%PROJECT_NAME",c_val=project_name,error=error)
    input_file_path = TRIM(project_name)
    lp=LEN_TRIM(input_file_path)
    i = logger%iter_info%iteration(logger%iter_info%n_rlevel)
    input_file_path(lp+1:LEN(input_file_path))="-"//TRIM(label)//"-"//ADJUSTL(cp_to_string(i))
    lp=LEN_TRIM(input_file_path)
    CALL section_vals_val_set(root_section,"GLOBAL%PROJECT_NAME",c_val=input_file_path(1:lp),error=error)
    CALL section_vals_val_set(root_section,"GLOBAL%RUN_TYPE",i_val=id_run,error=error)

    ! Redirecting output into a new file
    output_file_path=input_file_path(1:lp)//".out"
    IF (para_env%mepos==para_env%source) THEN
       CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
            file_action="WRITE", file_position="APPEND",unit_number=unit_nr)
    ELSE
       unit_nr=-1
    END IF
    CALL cp_logger_create(loc_logger, para_env=para_env, default_global_unit_nr=unit_nr, &
         close_global_unit_on_dealloc=.FALSE.)
    CALL section_vals_val_get(root_section,"GLOBAL%PROJECT",c_val=c_val,error=error)
    IF (c_val/="") THEN
       CALL cp_logger_set(loc_logger,local_filename=TRIM(c_val)//"_localLog")
    END IF
    loc_logger%iter_info%project_name=c_val
    CALL section_vals_val_get(root_section,"GLOBAL%PRINT_LEVEL",&
       i_val=loc_logger%iter_info%print_level,error=error)

    ! Initialize error
    CALL cp_error_init(new_error, stop_level=cp_fatal_level, logger=loc_logger)
    CALL cp_logger_release(loc_logger)

  END SUBROUTINE gopt_new_error_create

! *****************************************************************************
!> \param releases a new error used for cell optimization algorithm
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE gopt_new_error_release(new_error, root_section, para_env, project_name, id_run, error)
    TYPE(cp_error_type), INTENT(inout)       :: new_error
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(len=default_string_length), &
      INTENT(IN)                             :: project_name
    INTEGER, INTENT(IN)                      :: id_run
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gopt_new_error_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(new_error)
    IF (para_env%mepos==para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger)
       CALL close_file(unit_number=unit_nr)
    END IF
    CALL cp_error_dealloc_ref(new_error, error)
    CALL section_vals_val_set(root_section,"GLOBAL%RUN_TYPE",i_val=id_run,error=error)
    CALL section_vals_val_set(root_section,"GLOBAL%PROJECT_NAME",c_val=project_name,error=error)

  END SUBROUTINE gopt_new_error_release

! *****************************************************************************
!> \param Reads the external pressure tensor
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE read_external_press_tensor(geo_section, cell, pres_ext, mtrx, rot, error)
    TYPE(section_vals_type), POINTER         :: geo_section
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(OUT)               :: pres_ext
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: mtrx
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: rot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_external_press_tensor', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind, j
    LOGICAL                                  :: check, failure
    REAL(KIND=dp), DIMENSION(3, 3)           :: pres_ext_tens
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pvals

    failure = .FALSE.
    NULLIFY(pvals)
    mtrx          = 0.0_dp
    pres_ext_tens = 0.0_dp
    pres_ext      = 0.0_dp
    CALL section_vals_val_get(geo_section,"EXTERNAL_PRESSURE", r_vals=pvals, error=error)
    check = (SIZE(pvals)==1).OR.(SIZE(pvals)==9)
    CALL cp_assert(check,cp_failure_level,cp_assertion_failed,routineP,&
         "EXTERNAL_PRESSURE can have 1 or 9 components only! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    IF (SIZE(pvals)==9) THEN
       ind=0
       DO i = 1,3
          DO j = 1,3
             ind=ind+1
             pres_ext_tens(j,i) = pvals(ind)
          END DO
       END DO
       ! Also the pressure tensor must be oriented in the same canonical directions
       ! of the simulation cell
       pres_ext_tens = MATMUL(TRANSPOSE(rot),pres_ext_tens)
       DO i=1,3
          pres_ext=pres_ext+pres_ext_tens(i,i)
       ENDDO
       pres_ext=pres_ext/3.0_dp
       DO i=1,3
          pres_ext_tens(i,i)=pres_ext_tens(i,i)-pres_ext
       ENDDO
    ELSE
       pres_ext=pvals(1)
    END IF

    IF(ANY(pres_ext_tens > 1.0E-5_dp))THEN
       mtrx=cell%deth*MATMUL(cell%h_inv,MATMUL(pres_ext_tens,TRANSPOSE(cell%h_inv)))
    END IF

  END SUBROUTINE read_external_press_tensor

! *****************************************************************************
!> \param Computes the derivatives for the cell
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE get_dg_dh(gradient,av_ptens,pres_ext,cell,mtrx,keep_angles,&
                       keep_symmetry,pres_int,error)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: av_ptens
    REAL(KIND=dp), INTENT(IN)                :: pres_ext
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mtrx
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_angles, keep_symmetry
    REAL(KIND=dp), INTENT(OUT)               :: pres_int
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dg_dh', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_keep_angles, &
                                                my_keep_symmetry
    REAL(KIND=dp)                            :: a, cosa, cosah, g, norm, &
                                                norm_b, norm_c, sina, sinah
    REAL(KIND=dp), DIMENSION(3, 3)           :: correction, pten_hinv_old, &
                                                ptens

    failure = .FALSE.
    my_keep_angles = .FALSE.
    IF (PRESENT(keep_angles)) my_keep_angles =  keep_angles
    my_keep_symmetry = .FALSE.
    IF (PRESENT(keep_symmetry)) my_keep_symmetry =  keep_symmetry
    gradient = 0.0_dp
    ptens = av_ptens

    ! Evaluating the internal pressure
    pres_int = 0.0_dp
    DO i=1,3
       pres_int = pres_int + ptens(i,i)
    ENDDO
    pres_int = pres_int/3.0_dp

    ptens(1,1) = av_ptens(1,1) - pres_ext
    ptens(2,2) = av_ptens(2,2) - pres_ext
    ptens(3,3) = av_ptens(3,3) - pres_ext

    pten_hinv_old = cell%deth * MATMUL(cell%h_inv, ptens)
    correction = MATMUL(mtrx , cell%hmat)

    gradient(1) = pten_hinv_old(1,1) - correction(1,1)
    gradient(2) = pten_hinv_old(2,1) - correction(2,1)
    gradient(3) = pten_hinv_old(2,2) - correction(2,2)
    gradient(4) = pten_hinv_old(3,1) - correction(3,1)
    gradient(5) = pten_hinv_old(3,2) - correction(3,2)
    gradient(6) = pten_hinv_old(3,3) - correction(3,3)

    IF (my_keep_angles) THEN
       ! If we want to keep the angles constant we have to project out the
       ! components of the cell angles
       norm_b        = DOT_PRODUCT(cell%hmat(:,2),cell%hmat(:,2))
       norm          = DOT_PRODUCT(cell%hmat(1:2,2),gradient(2:3))
       gradient(2:3) = cell%hmat(1:2,2)/norm_b*norm
       norm_c        = DOT_PRODUCT(cell%hmat(:,3),cell%hmat(:,3))
       norm          = DOT_PRODUCT(cell%hmat(1:3,3),gradient(4:6))
       gradient(4:6) = cell%hmat(1:3,3)/norm_c*norm
       ! Retain an exact orthorhombic cell
       ! (off-diagonal elements must remain zero identically to keep QS fast)
       IF (cell%orthorhombic) THEN
          gradient(2) = 0.0_dp
          gradient(4) = 0.0_dp
          gradient(5) = 0.0_dp
       END IF
    END IF

    IF (my_keep_symmetry) THEN
       SELECT CASE (cell%symmetry_id)
       CASE (cell_sym_cubic,&
             cell_sym_tetragonal_ab,&
             cell_sym_tetragonal_ac,&
             cell_sym_tetragonal_bc,&
             cell_sym_orthorhombic)
          SELECT CASE (cell%symmetry_id)
          CASE (cell_sym_cubic)
             g = (gradient(1) + gradient(3) + gradient(6))/3.0_dp
             gradient(1) = g
             gradient(3) = g
             gradient(6) = g
          CASE (cell_sym_tetragonal_ab,&
                cell_sym_tetragonal_ac,&
                cell_sym_tetragonal_bc)
             SELECT CASE (cell%symmetry_id)
             CASE (cell_sym_tetragonal_ab)
                g = 0.5_dp*(gradient(1) + gradient(3))
                gradient(1) = g
                gradient(3) = g
             CASE (cell_sym_tetragonal_ac)
                g = 0.5_dp*(gradient(1) + gradient(6))
                gradient(1) = g
                gradient(6) = g
             CASE (cell_sym_tetragonal_bc)
                g = 0.5_dp*(gradient(3) + gradient(6))
                gradient(3) = g
                gradient(6) = g
             END SELECT
          CASE (cell_sym_orthorhombic)
             ! Nothing else to do
          END SELECT
          gradient(2) = 0.0_dp
          gradient(4) = 0.0_dp
          gradient(5) = 0.0_dp
       CASE (cell_sym_hexagonal)
          g = 0.5_dp*(gradient(1) + 0.5_dp*(gradient(2) + sqrt3*gradient(3)))
          gradient(1) = g
          gradient(2) = 0.5_dp*g
          gradient(3) = sqrt3*gradient(2)
          gradient(4) = 0.0_dp
          gradient(5) = 0.0_dp
       CASE (cell_sym_rhombohedral)
          a = (angle(cell%hmat(:,3),cell%hmat(:,2)) +&
               angle(cell%hmat(:,1),cell%hmat(:,3)) +&
               angle(cell%hmat(:,1),cell%hmat(:,2)))/3.0_dp
          cosa = COS(a)
          sina = SIN(a)
          cosah = COS(0.5_dp*a)
          sinah = SIN(0.5_dp*a)
          norm = cosa/cosah
          norm_c = SQRT(1.0_dp - norm*norm)
          g = (gradient(1) + gradient(2)*cosa + gradient(3)*sina +&
               gradient(4)*cosah*norm + gradient(5)*sinah*norm + gradient(6)*norm_c)/3.0_dp
          gradient(1) = g
          gradient(2) = g*cosa
          gradient(3) = g*sina
          gradient(4) = g*cosah*norm
          gradient(5) = g*sinah*norm
          gradient(6) = g*norm_c
       CASE (cell_sym_monoclinic)
          gradient(2) = 0.0_dp
          gradient(5) = 0.0_dp
       CASE (cell_sym_triclinic)
          ! Nothing to do
       END SELECT
    END IF

    gradient = -gradient

  END SUBROUTINE get_dg_dh

END MODULE cell_opt_utils
