!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A generic framework to calculate step lengths for 1D line search
!> \author Ole Schuett
! **************************************************************************************************
MODULE linesearch
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'linesearch'

   PUBLIC :: linesearch_type

   INTEGER, PARAMETER :: linesearch_method_adapt = 1
   INTEGER, PARAMETER :: linesearch_method_2pnt = 2
   INTEGER, PARAMETER :: linesearch_method_3pnt = 3
   INTEGER, PARAMETER :: linesearch_method_gold = 4
   INTEGER, PARAMETER :: linesearch_method_none = 5

   TYPE linesearch_2pnt_type
      REAL(KIND=dp), DIMENSION(2)                    :: energies = 0.0
      REAL(KIND=dp)                                  :: scan_step = 0.0
      REAL(KIND=dp)                                  :: last_step_size = 0.0
      REAL(KIND=dp)                                  :: max_step_size = 0.0
      INTEGER                                        :: count = 1
   END TYPE linesearch_2pnt_type

   TYPE linesearch_3pnt_type
      REAL(KIND=dp), DIMENSION(3)                    :: energies = 0.0
      REAL(KIND=dp), DIMENSION(3)                    :: scan_steps = 0.0
      REAL(KIND=dp)                                  :: last_step_size = 0.0
      REAL(KIND=dp)                                  :: max_step_size = 0.0
      REAL(KIND=dp)                                  :: tiny_step_size = 0.0
      INTEGER                                        :: count = 1
   END TYPE linesearch_3pnt_type

   TYPE linesearch_adapt_type
      REAL(KIND=dp)                                  :: last_step_size = 0.0
      REAL(KIND=dp)                                  :: left_x = 0.0
      REAL(KIND=dp)                                  :: middle_x = 0.0
      REAL(KIND=dp)                                  :: right_x = 0.0
      REAL(KIND=dp)                                  :: left_e = 0.0
      REAL(KIND=dp)                                  :: middle_e = 0.0
      REAL(KIND=dp)                                  :: right_e = 0.0
      LOGICAL                                        :: have_left = .FALSE.
      LOGICAL                                        :: have_middle = .FALSE.
      LOGICAL                                        :: have_right = .FALSE.
      INTEGER                                        :: count = 0
   END TYPE linesearch_adapt_type

   TYPE linesearch_gold_type
      REAL(KIND=dp)                                  :: scan_steps = 0.0
      REAL(KIND=dp)                                  :: last_step_size = 0.0
      REAL(KIND=dp)                                  :: eps_step_size = 0.0
      REAL(KIND=dp)                                  :: left_x = 0.0
      REAL(KIND=dp)                                  :: middle_x = 0.0
      REAL(KIND=dp)                                  :: right_x = 0.0
      REAL(KIND=dp)                                  :: left_e = 0.0
      REAL(KIND=dp)                                  :: middle_e = 0.0
      REAL(KIND=dp)                                  :: right_e = 0.0
      LOGICAL                                        :: have_left = .FALSE.
      LOGICAL                                        :: have_middle = .FALSE.
      LOGICAL                                        :: have_right = .FALSE.
      LOGICAL                                        :: gave_up = .FALSE.
   END TYPE linesearch_gold_type

   TYPE linesearch_type
      PRIVATE
      REAL(KIND=dp), PUBLIC                          :: step_size = 0.0_dp
      LOGICAL, PUBLIC                                :: starts = .FALSE.
      TYPE(linesearch_adapt_type), POINTER           :: state_adapt => Null()
      TYPE(linesearch_2pnt_type), POINTER            :: state_2pnt => Null()
      TYPE(linesearch_3pnt_type), POINTER            :: state_3pnt => Null()
      TYPE(linesearch_gold_type), POINTER            :: state_gold => Null()
      INTEGER                                        :: iw = -1
      INTEGER                                        :: method = -1
      CHARACTER(LEN=10)                              :: label = ""
      REAL(KIND=dp)                                  :: init_step_size = 0.0_dp
      REAL(dp)                                       :: eps_step_size = 0.0_dp
      REAL(dp)                                       :: max_step_size = 0.0_dp
      REAL(dp)                                       :: tiny_step_size = 0.0_dp
   END TYPE linesearch_type

   PUBLIC :: linesearch_create_section, linesearch_init, linesearch_finalize
   PUBLIC :: linesearch_step, linesearch_reset

CONTAINS

! **************************************************************************************************
!> \brief Declare the line search input section.
!> \param section ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_create_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_create_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_section, printkey

      NULLIFY (keyword, print_section, printkey)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "LINE_SEARCH", repeats=.FALSE., &
                          description="Detail settings or linesearch method.")

      CALL keyword_create(keyword, name="METHOD", &
                          description="Linesearch method.", &
                          default_i_val=linesearch_method_adapt, &
                          enum_c_vals=s2a("ADAPT", "3PNT", "2PNT", "GOLD", "NONE"), &
                          enum_desc=s2a("extrapolates usually based on 3 points, uses additional points on demand, very robust.", &
                                        "extrapolate based on 3 points", &
                                 "extrapolate based on 2 points and the slope (super fast, but might get stuck at saddle points)", &
                                        "perform 1D golden section search of the minimum (very expensive)", &
                                        "always take steps of fixed INITIAL_STEP_SIZE"), &
                          enum_i_vals=(/linesearch_method_adapt, linesearch_method_3pnt, linesearch_method_2pnt, &
                                        linesearch_method_gold, linesearch_method_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="INITIAL_STEP_SIZE", &
                          description="Initial step length", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MAX_STEP_SIZE", &
                          description="Maximum step length", &
                          default_r_val=3.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TINY_STEP_SIZE", &
                          description="Step length taken if negative step is suggested.", &
                          default_r_val=0.002_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="EPS_STEP_SIZE", &
                          description="Convergence criterion of GOLD method.", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(print_section, name="PRINT", &
                          description="Print section", &
                          n_keywords=0, n_subsections=1, repeats=.TRUE.)

      CALL cp_print_key_section_create(printkey, "RUN_INFO", &
                                       description="General run information", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")

      CALL section_add_subsection(print_section, printkey)
      CALL section_release(printkey)
      CALL section_add_subsection(section, print_section)
      CALL section_release(print_section)

   END SUBROUTINE linesearch_create_section

! **************************************************************************************************
!> \brief Initialize linesearch from given input section
!> \param this ...
!> \param section ...
!> \param label ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_init(this, section, label)
      TYPE(linesearch_type), INTENT(INOUT)               :: this
      TYPE(section_vals_type), POINTER                   :: section
      CHARACTER(LEN=*)                                   :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_init', &
         routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger

      CALL section_vals_val_get(section, "METHOD", i_val=this%method)
      CALL section_vals_val_get(section, "INITIAL_STEP_SIZE", r_val=this%init_step_size)
      CALL section_vals_val_get(section, "MAX_STEP_SIZE", r_val=this%max_step_size)
      CALL section_vals_val_get(section, "TINY_STEP_SIZE", r_val=this%tiny_step_size)
      CALL section_vals_val_get(section, "EPS_STEP_SIZE", r_val=this%eps_step_size)

      CPASSERT(LEN_TRIM(label) <= 10)
      this%label = label
      logger => cp_get_default_logger()
      this%iw = cp_print_key_unit_nr(logger, section, "PRINT%RUN_INFO", &
                                     extension=".linesearchlog")

      CALL linesearch_init_low(this)

   END SUBROUTINE linesearch_init

! **************************************************************************************************
!> \brief Helper routine to (re)-initialize line search machinery
!> \param this ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_init_low(this)
      TYPE(linesearch_type), INTENT(INOUT)               :: this

      this%step_size = 0.0_dp
      this%starts = .TRUE.

      SELECT CASE (this%method)
      CASE (linesearch_method_adapt)
         ALLOCATE (this%state_adapt)
         this%state_adapt%last_step_size = this%init_step_size
      CASE (linesearch_method_2pnt)
         ALLOCATE (this%state_2pnt)
         this%state_2pnt%max_step_size = this%max_step_size
         this%state_2pnt%last_step_size = this%init_step_size
      CASE (linesearch_method_3pnt)
         ALLOCATE (this%state_3pnt)
         this%state_3pnt%last_step_size = this%init_step_size
         this%state_3pnt%max_step_size = this%max_step_size
         this%state_3pnt%tiny_step_size = this%tiny_step_size
      CASE (linesearch_method_gold)
         ALLOCATE (this%state_gold)
         this%state_gold%last_step_size = this%init_step_size
         this%state_gold%eps_step_size = this%eps_step_size
      CASE (linesearch_method_none)
         ! nothing todo
      CASE DEFAULT
         CPABORT("unkown method")
      END SELECT

   END SUBROUTINE linesearch_init_low

! **************************************************************************************************
!> \brief Finzalize line search machinery
!> \param this ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_finalize(this)
      TYPE(linesearch_type), INTENT(INOUT)               :: this

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_finalize', &
         routineP = moduleN//':'//routineN

      SELECT CASE (this%method)
      CASE (linesearch_method_adapt)
         DEALLOCATE (this%state_adapt)
      CASE (linesearch_method_2pnt)
         DEALLOCATE (this%state_2pnt)
      CASE (linesearch_method_3pnt)
         DEALLOCATE (this%state_3pnt)
      CASE (linesearch_method_gold)
         DEALLOCATE (this%state_gold)
      CASE (linesearch_method_none)
         ! nothing todo
      CASE DEFAULT
         CPABORT("unkown method")
      END SELECT

      !TODO: should finish printkey, but don't have the section here
   END SUBROUTINE linesearch_finalize

! **************************************************************************************************
!> \brief Reset line search to initial state
!> \param this ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_reset(this)
      TYPE(linesearch_type), INTENT(INOUT)               :: this

      CALL linesearch_finalize(this)
      CALL linesearch_init_low(this)
   END SUBROUTINE linesearch_reset

! **************************************************************************************************
!> \brief Calculate step length of next line search step.
!> \param this ...
!> \param energy ...
!> \param slope ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_step(this, energy, slope)
      TYPE(linesearch_type), INTENT(INOUT)               :: this
      REAL(KIND=dp), INTENT(IN)                          :: energy, slope

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_step', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: is_done
      REAL(KIND=dp)                                      :: step_size

      SELECT CASE (this%method)
      CASE (linesearch_method_adapt)
         CALL linesearch_adapt(this%state_adapt, energy, step_size, is_done, this%iw, TRIM(this%label))
      CASE (linesearch_method_2pnt)
         CALL linesearch_2pnt(this%state_2pnt, energy, slope, step_size, is_done, this%iw, TRIM(this%label))
      CASE (linesearch_method_3pnt)
         CALL linesearch_3pnt(this%state_3pnt, energy, step_size, is_done, this%iw, TRIM(this%label))
      CASE (linesearch_method_gold)
         CALL linesearch_gold(this%state_gold, energy, step_size, is_done, this%iw, TRIM(this%label))
      CASE (linesearch_method_none)
         step_size = this%init_step_size ! take steps of fixed length
         is_done = .TRUE.
      CASE DEFAULT
         CPABORT("unkown method")
      END SELECT

      this%step_size = step_size
      this%starts = is_done
   END SUBROUTINE linesearch_step

! **************************************************************************************************
!> \brief Perform a 2pnt linesearch
!> \param this ...
!> \param energy ...
!> \param slope ...
!> \param step_size ...
!> \param is_done ...
!> \param unit_nr ...
!> \param label ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_2pnt(this, energy, slope, step_size, is_done, unit_nr, label)
      TYPE(linesearch_2pnt_type), INTENT(INOUT)          :: this
      REAL(KIND=dp), INTENT(IN)                          :: energy, slope
      REAL(KIND=dp), INTENT(OUT)                         :: step_size
      LOGICAL, INTENT(OUT)                               :: is_done
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN)                       :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_2pnt', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: a, b, c, pred_energy, x2

      this%energies(this%count) = energy
      is_done = .FALSE.

      SELECT CASE (this%count)
      CASE (1)
         step_size = 2.0_dp*this%last_step_size
         this%scan_step = step_size
         this%count = 2
      CASE (2)
         c = this%energies(1)
         b = -slope
         x2 = this%scan_step
         a = (this%energies(2)-b*x2-c)/(x2**2)

         IF (a < 0.0_dp) THEN
            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| had to quench curvature"
            a = 1.0E-15_dp
         ENDIF

         step_size = -b/(2.0_dp*a)
         pred_energy = a*step_size**2+b*step_size+c

         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| 2pnt suggested step_size: ", step_size
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| 2pnt predicted energy", pred_energy

         IF (pred_energy > this%energies(1) .OR. pred_energy > this%energies(2)) THEN
            CPABORT(label//"LS| predicted energy not below test points")
         ENDIF

         IF (step_size > this%max_step_size) THEN
            step_size = this%max_step_size
            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| limiting step_size to MAX_STEP_SIZE"
         ENDIF

         this%last_step_size = step_size
         this%count = 1
         is_done = .TRUE.
      CASE DEFAULT
         CPABORT("this should not happen")
      END SELECT

   END SUBROUTINE linesearch_2pnt

! **************************************************************************************************
!> \brief Perform a 3pnt linesearch
!> \param this ...
!> \param energy ...
!> \param step_size ...
!> \param is_done ...
!> \param unit_nr ...
!> \param label ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_3pnt(this, energy, step_size, is_done, unit_nr, label)
      TYPE(linesearch_3pnt_type), INTENT(INOUT)          :: this
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), INTENT(OUT)                         :: step_size
      LOGICAL, INTENT(OUT)                               :: is_done
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN)                       :: label

      REAL(KIND=dp)                                      :: a, b, c, denom, pred_energy, x1, x2, x3, &
                                                            y1, y2, y3

      this%energies(this%count) = energy
      is_done = .FALSE.

      SELECT CASE (this%count)
      CASE (1)
         step_size = (2.0_dp/3.0_dp)*this%last_step_size
         IF (step_size < this%tiny_step_size) THEN
            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| initial step size too small, using TINY_STEP_SIZE"
            step_size = this%tiny_step_size
         ENDIF
         this%scan_steps(1) = 0.0_dp
         this%scan_steps(2) = step_size
         this%count = 2
      CASE (2)
         IF (this%energies(1) > this%energies(2)) THEN
            step_size = 2.0_dp*this%scan_steps(2)
         ELSE
            step_size = 0.5_dp*this%scan_steps(2)
         END IF
         this%scan_steps(3) = step_size
         this%count = 3
      CASE (3)
         ! fitting y = a*x^2 + b*x + c
         y1 = this%energies(1); y2 = this%energies(2); y3 = this%energies(3)
         x1 = this%scan_steps(1); x2 = this%scan_steps(2); x3 = this%scan_steps(3)

         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| 3pnt scan_steps: ", this%scan_steps
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| 3pnt energies: ", this%energies

         ! Cramer's Rule
         denom = (x1-x2)*(x1-x3)*(x2-x3)
         a = (x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))/denom
         b = (x3**2*(y1-y2)+x2**2*(y3-y1)+x1**2*(y2-y3))/denom
         c = (x2*x3*(x2-x3)*y1+x3*x1*(x3-x1)*y2+x1*x2*(x1-x2)*y3)/denom

         step_size = -b/(2.0_dp*a)
         pred_energy = a*step_size**2+b*step_size+c
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| 3pnt suggested step_size: ", step_size
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| 3pnt predicted energy", pred_energy

         IF (a < 0) THEN
            step_size = -2.0_dp*step_size
            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| inverting step size"
         END IF

         IF (step_size < 0) THEN
            step_size = this%tiny_step_size
            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| makeing a step of size TINY_STEP_SIZE"
         ENDIF

         IF (step_size > this%max_step_size) THEN
            step_size = this%max_step_size
            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| limiting step_size to MAX_STEP_SIZE"
         ENDIF

         this%last_step_size = step_size
         this%count = 1
         is_done = .TRUE.
      CASE DEFAULT
         CPABORT("this should not happen")
      END SELECT

   END SUBROUTINE linesearch_3pnt

! **************************************************************************************************
!> \brief Perform an adaptive linesearch
!> \param this ...
!> \param energy ...
!> \param step_size ...
!> \param is_done ...
!> \param unit_nr ...
!> \param label ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_adapt(this, energy, step_size, is_done, unit_nr, label)
      TYPE(linesearch_adapt_type), INTENT(INOUT)         :: this
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), INTENT(OUT)                         :: step_size
      LOGICAL, INTENT(OUT)                               :: is_done
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN)                       :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_adapt', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: grow_factor = 2.0_dp, &
                                                            shrink_factor = 0.5_dp

      REAL(KIND=dp)                                      :: a, b, c, denom, pred_energy, x1, x2, x3, &
                                                            y1, y2, y3

      is_done = .FALSE.
      this%count = this%count+1

      IF (.NOT. this%have_left) THEN
         this%left_x = 0.0_dp
         this%left_e = energy
         this%have_left = .TRUE.
         step_size = this%last_step_size

      ELSE IF (.NOT. (this%have_middle .OR. this%have_right)) THEN
         IF (energy < this%left_e) THEN
            this%middle_e = energy
            this%middle_x = this%last_step_size
            this%have_middle = .TRUE.
            step_size = this%middle_x*grow_factor
         ELSE
            this%right_e = energy
            this%right_x = this%last_step_size
            this%have_right = .TRUE.
            step_size = this%right_x*shrink_factor
         ENDIF

      ELSE IF (.NOT. this%have_right) THEN
         IF (energy < this%middle_e) THEN
            this%middle_e = energy
            this%middle_x = this%last_step_size
            step_size = this%middle_x*grow_factor
         ELSE
            this%right_e = energy
            this%right_x = this%last_step_size
            this%have_right = .TRUE.
         ENDIF

      ELSE IF (.NOT. this%have_middle) THEN
         IF (energy > this%left_e) THEN
            this%right_e = energy
            this%right_x = this%last_step_size
            step_size = this%right_x*shrink_factor
         ELSE
            this%middle_e = energy
            this%middle_x = this%last_step_size
            this%have_middle = .TRUE.
         ENDIF
      ENDIF

      IF (this%count > 3) THEN
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| Need extra step"
      ENDIF
      IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| adapt: ", this%have_left, this%have_middle, this%have_right
      IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| adapt: scan_steps: ", this%left_x, this%middle_x, this%right_x
      IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| adapt: energies: ", this%left_e, this%middle_e, this%right_e

      IF (this%have_left .AND. this%have_middle .AND. this%have_right) THEN
         ! fitting y = a*x^2 + b*x + c
         y1 = this%left_e; y2 = this%middle_e; y3 = this%right_e
         x1 = this%left_x; x2 = this%middle_x; x3 = this%right_x

         ! Cramer's rule
         denom = (x1-x2)*(x1-x3)*(x2-x3)
         a = (x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))/denom
         b = (x3**2*(y1-y2)+x2**2*(y3-y1)+x1**2*(y2-y3))/denom
         c = (x2*x3*(x2-x3)*y1+x3*x1*(x3-x1)*y2+x1*x2*(x1-x2)*y3)/denom

         IF (ABS(a) /= 0.0_dp) THEN
            step_size = -b/(2.0_dp*a)
         ELSE
            step_size = 0.0_dp
         ENDIF

         CPASSERT(step_size >= 0.0_dp)
         pred_energy = a*step_size**2+b*step_size+c
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| adapt: suggested step_size: ", step_size
         IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| adapt: predicted energy", pred_energy

         ! reset
         is_done = .TRUE.
         this%count = 0
         this%have_left = .FALSE.
         this%have_middle = .FALSE.
         this%have_right = .FALSE.
         this%left_e = 0.0
         this%middle_e = 0.0
         this%right_e = 0.0
         this%left_x = 0.0
         this%middle_x = 0.0
         this%right_x = 0.0
      END IF

      this%last_step_size = step_size
   END SUBROUTINE linesearch_adapt

! **************************************************************************************************
!> \brief Perform a gold linesearch
!> \param this ...
!> \param energy ...
!> \param step_size ...
!> \param is_done ...
!> \param unit_nr ...
!> \param label ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE linesearch_gold(this, energy, step_size, is_done, unit_nr, label)
      TYPE(linesearch_gold_type), INTENT(INOUT)          :: this
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), INTENT(OUT)                         :: step_size
      LOGICAL, INTENT(OUT)                               :: is_done
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN)                       :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'linesearch_gold', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER :: phi = (1.0_dp+SQRT(5.0_dp))/2.0_dp

      REAL(KIND=dp)                                      :: a, b, d

      is_done = .FALSE.

      IF (this%gave_up) &
         CPABORT("had to give up, should not be called again")

      IF (.NOT. this%have_left) THEN
         this%left_x = 0.0_dp
         this%left_e = energy
         this%have_left = .TRUE.
         step_size = this%last_step_size

      ELSE IF (.NOT. (this%have_middle .OR. this%have_right)) THEN
         IF (energy < this%left_e) THEN
            this%middle_e = energy
            this%middle_x = this%scan_steps
            this%have_middle = .TRUE.
            step_size = this%middle_x*phi
         ELSE
            this%right_e = energy
            this%right_x = this%scan_steps
            this%have_right = .TRUE.
            step_size = this%right_x/phi
         ENDIF

      ELSE IF (.NOT. this%have_right) THEN
         IF (energy < this%middle_e) THEN
            this%middle_e = energy
            this%middle_x = this%scan_steps
            step_size = this%middle_x*phi
         ELSE
            this%right_e = energy
            this%right_x = this%scan_steps
            this%have_right = .TRUE.
         ENDIF

      ELSE IF (.NOT. this%have_middle) THEN
         IF (energy > this%left_e) THEN
            this%right_e = energy
            this%right_x = this%scan_steps
            step_size = this%right_x/phi
         ELSE
            this%middle_e = energy
            this%middle_x = this%scan_steps
            this%have_middle = .TRUE.
         ENDIF

      ELSE !up and running
         a = this%middle_x-this%left_x
         b = this%right_x-this%middle_x
         IF (energy < this%middle_e) THEN
            IF (a < b) THEN
               this%left_e = this%middle_e
               this%left_x = this%middle_x
            ELSE
               this%right_e = this%middle_e
               this%right_x = this%middle_x
            ENDIF
            this%middle_e = energy
            this%middle_x = this%scan_steps
         ELSE
            IF (a < b) THEN
               this%right_e = energy
               this%right_x = this%scan_steps
            ELSE
               this%left_e = energy
               this%left_x = this%scan_steps
            ENDIF
         ENDIF
      ENDIF

      IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| gold: ", this%have_left, this%have_middle, this%have_right
      IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| gold: ", this%left_x, this%middle_x, this%right_x
      IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| gold: ", this%left_e, this%middle_e, this%right_e

      IF (this%have_left .AND. this%have_middle .AND. this%have_right) THEN
         a = this%middle_x-this%left_x
         b = this%right_x-this%middle_x
         IF (ABS(MIN(a, b)*phi-MAX(a, b)) > 1.0E-10) &
            CPABORT("golden-ratio gone")

         IF (a < b) THEN
            step_size = this%middle_x+a/phi
         ELSE
            step_size = this%middle_x-b/phi
         END IF

         d = ABS(this%right_x-this%left_x)/(ABS(this%middle_x)+ABS(step_size))
         IF (d < this%eps_step_size) THEN
            step_size = this%middle_x
            this%last_step_size = step_size
            is_done = .TRUE.

            IF (unit_nr > 0) WRITE (unit_nr, *) label, "LS| gold done, step-size: ", step_size

            this%have_left = .FALSE.
            this%have_middle = .FALSE.
            this%have_right = .FALSE.
            this%left_e = 0.0
            this%middle_e = 0.0
            this%right_e = 0.0
            this%left_x = 0.0
            this%middle_x = 0.0
            this%right_x = 0.0
         ENDIF

      END IF

      IF (step_size < 1E-10) CPABORT("linesearch failed / done")

      this%scan_steps = step_size
   END SUBROUTINE linesearch_gold

END MODULE linesearch
