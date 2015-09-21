!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the Minima Crawling global optimization scheme
!> \author Ole Schuett
! *****************************************************************************
MODULE glbopt_mincrawl
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE glbopt_history,                  ONLY: history_add,&
                                             history_finalize,&
                                             history_fingerprint,&
                                             history_fingerprint_type,&
                                             history_init,&
                                             history_lookup,&
                                             history_type
  USE input_constants,                 ONLY: dump_xmol
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE particle_methods,                ONLY: write_particle_coordinates
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: kelvin
  USE swarm_message,                   ONLY: swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_type
#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_mincrawl'


 PUBLIC :: mincrawl_type
 PUBLIC :: mincrawl_init, mincrawl_finalize
 PUBLIC :: mincrawl_steer

 TYPE minima_type
   INTEGER                                             :: id = -1
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pos
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: escape_hist
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: tempdist
   REAL(KIND=dp)                                       :: Epot = -1.0
   TYPE(history_fingerprint_type)                      :: fp
   LOGICAL                                             :: disabled = .FALSE.
   INTEGER                                             :: n_active = 0
   INTEGER                                             :: n_sampled = 0
 END TYPE minima_type

 TYPE minima_p_type
   TYPE(minima_type), POINTER                          :: p => Null()
 END TYPE minima_p_type

 TYPE worker_state_type
   TYPE(minima_type), POINTER                          :: start_minima => Null()
   INTEGER                                             :: tempstep = 0
   INTEGER                                             :: iframe = 1
 END TYPE worker_state_type

 TYPE mincrawl_type
   PRIVATE
   TYPE(history_type)                                  :: history
   TYPE(worker_state_type), DIMENSION(:), ALLOCATABLE  :: workers
   TYPE(minima_p_type), DIMENSION(:), ALLOCATABLE      :: minimas
   REAL(KIND=dp)                                       :: tempstep_base = 0
   INTEGER                                             :: tempstep_max = 0
   REAL(KIND=dp)                                       :: tempdist_init_width = 0
   REAL(KIND=dp)                                       :: tempdist_update_width = 0
   REAL(KIND=dp)                                       :: tempdist_update_height = 0
   INTEGER                                             :: esc_hist_len = 0
   INTEGER                                             :: tempstep_init = 0
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: tempdist_init
   INTEGER                                             :: n_minima = 0
   INTEGER                                             :: n_workers = 0
   INTEGER                                             :: worker_per_min = 0
   INTEGER                                             :: iw = 0
   INTEGER                                             :: minima_traj_unit = 0
   TYPE(section_vals_type), POINTER                    :: mincrawl_section => Null()
   TYPE(rng_stream_type),POINTER                       :: rng_stream => Null()
   TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set => Null()
 END TYPE mincrawl_type

 CONTAINS


! *****************************************************************************
!> \brief Initializes master for Minima Crawling
!> \param this ...
!> \param glbopt_section ...
!> \param n_workers ...
!> \param iw ...
!> \param particle_set ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE mincrawl_init(this, glbopt_section, n_workers, iw, particle_set)
    TYPE(mincrawl_type)                      :: this
    TYPE(section_vals_type), POINTER         :: glbopt_section
    INTEGER, INTENT(IN)                      :: n_workers, iw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    INTEGER                                  :: i
    REAL(kind=dp)                            :: temp_in_kelvin
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: history_section

    NULLIFY(logger, history_section)

    ! read input
    this%mincrawl_section => section_vals_get_subs_vals(glbopt_section, "MINIMA_CRAWLING")
    CALL section_vals_val_get(this%mincrawl_section,"TEMPSTEP_BASE",         r_val=this%tempstep_base)
    CALL section_vals_val_get(this%mincrawl_section,"TEMPSTEP_MAX",          i_val=this%tempstep_max)
    CALL section_vals_val_get(this%mincrawl_section,"TEMPDIST_INIT_WIDTH",   r_val=this%tempdist_init_width)
    CALL section_vals_val_get(this%mincrawl_section,"TEMPDIST_UPDATE_WIDTH", r_val=this%tempdist_update_width)
    CALL section_vals_val_get(this%mincrawl_section,"TEMPDIST_UPDATE_HEIGHT",r_val=this%tempdist_update_height)
    CALL section_vals_val_get(this%mincrawl_section,"TEMPERATURE_INIT", r_val=temp_in_kelvin)
    this%tempstep_init = temp2tempstep(this, temp_in_kelvin / kelvin)
    CALL section_vals_val_get(this%mincrawl_section,"WORKER_PER_MINIMA", i_val=this%worker_per_min)
    CALL section_vals_val_get(this%mincrawl_section,"ESCAPE_HISTORY_LENGTH", i_val=this%esc_hist_len)

    !init minima trajectory
    logger => cp_get_default_logger()
    this%minima_traj_unit = cp_print_key_unit_nr(logger, &
      this%mincrawl_section,"MINIMA_TRAJECTORY", &
      middle_name="minima", extension=".xyz", &
      file_action="WRITE", file_position="REWIND")

    !init history
    history_section => section_vals_get_subs_vals(glbopt_section, "HISTORY")
    CALL history_init(this%history, history_section, iw=iw)

    !allocate data structures
    ALLOCATE(this%minimas(1000)) !will be grown if needed

    ALLOCATE(this%workers(n_workers))
    this%n_workers = n_workers
    this%iw = iw
    this%particle_set => particle_set

    ! call fermi-like stepfunction for initial temp-dist
    ALLOCATE(this%tempdist_init(this%tempstep_max))
    this%tempdist_init = 0.0
    DO i=1, this%tempstep_max
        this%tempdist_init(i) = 1.0 / (1.0 + EXP((this%tempstep_init - i)/this%tempdist_init_width))
    ENDDO

    CALL create_rng_stream(this%rng_stream, name="mincrawl")
 END SUBROUTINE mincrawl_init


! *****************************************************************************
!> \brief Central steering routine of Minima Crawling
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE mincrawl_steer(this, report, cmd)
    TYPE(mincrawl_type)                      :: this
    TYPE(swarm_message_type)                 :: report, cmd

    CHARACTER(len=default_string_length)     :: status
    INTEGER                                  :: wid
    TYPE(minima_type), POINTER               :: best_minima

    CALL swarm_message_get(report, "status", status)
    CALL swarm_message_get(report, "worker_id", wid)

    IF(TRIM(status) == "initial_hello") THEN
       this%workers(wid)%tempstep = this%tempstep_init
       CALL swarm_message_add(cmd, "command", "md_and_gopt")
       CALL swarm_message_add(cmd, "iframe", 1)
       CALL swarm_message_add(cmd, "temperature", tempstep2temp(this,this%workers(wid)%tempstep))
       RETURN
    ENDIF

    IF(TRIM(status) == "ok") &
      CALL mincrawl_register_minima(this, report)


    IF(.FALSE.) CALL print_tempdist(best_minima)

    best_minima => choose_promising_minima(this)

    IF(.NOT. ASSOCIATED(best_minima)) THEN ! no suitable minima found
      CALL swarm_message_add(cmd, "command", "wait")
      !WRITE(this%iw,*) " MINCRAWL| Waiting until new minima become available"
      RETURN
    ENDIF

    best_minima%n_active = best_minima%n_active + 1
    best_minima%n_sampled = best_minima%n_sampled + 1
    this%workers(wid)%start_minima => best_minima
    this%workers(wid)%tempstep = choose_tempstep(this, best_minima)

    CALL swarm_message_add(cmd, "command", "md_and_gopt")
    CALL swarm_message_add(cmd, "iframe",  this%workers(wid)%iframe)
    CALL swarm_message_add(cmd, "temperature", tempstep2temp(this, this%workers(wid)%tempstep))
    CALL swarm_message_add(cmd, "positions",  best_minima%pos)

    IF(this%iw > 0) THEN
      WRITE(this%iw,'(1X,A,T71,I10)') &
       "MINCRAWL| Total number of found minima", this%n_minima
      WRITE(this%iw,'(1X,A,T71,I10)') &
       "MINCRAWL| Sampling minima with id",best_minima%id
      WRITE(this%iw,'(1X,A,I10,A,A,T71,F10.3)')&
       "MINCRAWL| Temperature  (step ", this%workers(wid)%tempstep," ) ",&
       "[Kelvin]", kelvin * tempstep2temp(this, this%workers(wid)%tempstep)
    ENDIF

 END SUBROUTINE mincrawl_steer


! *****************************************************************************
!> \brief Helper routine for mincrawl_steer, choses minimum based on its score.
!> \param this ...
!> \retval minima ...
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION choose_promising_minima(this) RESULT(minima)
    TYPE(mincrawl_type)                      :: this
    TYPE(minima_type), POINTER               :: minima

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: score, score_best

    score_best = HUGE(1.0)
    NULLIFY(minima)

    DO i=1, this%n_minima
       IF(this%minimas(i)%p%disabled) CYCLE
       IF(this%minimas(i)%p%n_active > this%worker_per_min) CYCLE
       score = minima_score(this%minimas(i)%p)
!       WRITE (*,*) "Minima: ", i, " active: ",this%minimas(i)%active, " E_expect: ", E_expect
       IF(score < score_best) THEN
          score_best = score
          minima => this%minimas(i)%p
       ENDIF
    END DO
 END FUNCTION choose_promising_minima


! *****************************************************************************
!> \brief Helper routine for choose_promising_minima, calculates a minimum's score
!> \param minima ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION minima_score(minima) RESULT(res)
    TYPE(minima_type), POINTER               :: minima
    REAL(KIND=dp)                            :: res

    res = SUM(minima%escape_hist) / SIZE(minima%escape_hist)
 END FUNCTION minima_score


! *****************************************************************************
!> \brief Helper routine for mincrawl_steer, samples from a temp-dist.
!> \param this ...
!> \param minima ...
!> \retval step ...
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION choose_tempstep(this, minima) RESULT(step)
    TYPE(mincrawl_type)                      :: this
    TYPE(minima_type), POINTER               :: minima
    INTEGER                                  :: step

    REAL(KIND=dp)                            :: a, r

   DO
      r = next_random_number(this%rng_stream)
      step = INT(r*SIZE(minima%tempdist)) + 1
      a = 1.0 -  2.0*ABS(minima%tempdist(step)-0.5)
      r = next_random_number(this%rng_stream)
      IF(r < a) EXIT
   END DO

 END FUNCTION choose_tempstep


! *****************************************************************************
!> \brief Debugging routine, prints a minimum's temp-distribution.
!> \param minima ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE print_tempdist(minima)
    TYPE(minima_type), POINTER               :: minima

    INTEGER                                  :: i

!WRITE (*,*) "tempdist: ", SUM(minima%tempdist, DIM=1)

   DO i=1, SIZE(minima%tempdist)
      WRITE (*,*) "tempstep: ",i,  minima%tempdist(i)
   ENDDO
 END SUBROUTINE print_tempdist


! *****************************************************************************
!> \brief Helper routine, convertes a  discrete temp-step to a temperature.
!> \param this ...
!> \param step ...
!> \retval temp_in_au ...
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION tempstep2temp(this, step) RESULT(temp_in_au)
    TYPE(mincrawl_type)                      :: this
    INTEGER                                  :: step
    REAL(KIND=dp)                            :: temp_in_au

   temp_in_au = ( this%tempstep_base**step )/ kelvin
 END FUNCTION tempstep2temp


! *****************************************************************************
!> \brief Helper routine, convertes a temperature to a discrete temp-step.
!> \param this ...
!> \param temp_in_au ...
!> \retval step ...
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION temp2tempstep(this, temp_in_au) RESULT(step)
    TYPE(mincrawl_type)                      :: this
    REAL(KIND=dp)                            :: temp_in_au
    INTEGER                                  :: step

    CHARACTER(len=*), PARAMETER :: routineN = 'temp2tempstep', &
      routineP = moduleN//':'//routineN

   step = INT(LOG(temp_in_au*kelvin)/LOG(this%tempstep_base))
   !WRITE(*,*) "temp: ", temp_in_au*kelvin, this%tempstep_base
   !WRITE(*,*) "step: ", step
   IF(step > this%tempstep_max) CPABORT("temp2tempstep: step > tempstep_max")
 END FUNCTION temp2tempstep



! *****************************************************************************
!> \brief Helper routine for mincrawl_steer
!>        Incorporates information of new report into history.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE mincrawl_register_minima(this, report)
    TYPE(mincrawl_type)                      :: this
    TYPE(swarm_message_type)                 :: report

    CHARACTER(len=*), PARAMETER :: routineN = 'mincrawl_register_minima', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: new_mid, tempstep, wid
    LOGICAL                                  :: minima_known
    REAL(KIND=dp)                            :: report_Epot
    REAL(KIND=dp), DIMENSION(:), POINTER     :: report_positions
    TYPE(history_fingerprint_type)           :: report_fp
    TYPE(minima_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: minimas_tmp
    TYPE(minima_type), POINTER               :: new_minima, start_minima

    NULLIFY(start_minima, new_minima, report_positions)

    CALL swarm_message_get(report, "worker_id", wid)
    CALL swarm_message_get(report, "Epot", report_Epot)
    CALL swarm_message_get(report, "positions", report_positions)
    CALL swarm_message_get(report, "iframe", this%workers(wid)%iframe)

    start_minima => this%workers(wid)%start_minima
    tempstep = this%workers(wid)%tempstep

    report_fp = history_fingerprint(report_Epot, report_positions)
    CALL history_lookup(this%history, report_fp, minima_known)

    IF(ASSOCIATED(start_minima)) THEN
      start_minima%n_active = start_minima%n_active - 1
      IF(start_minima%n_active < 0) CPABORT("negative n_active")

      ! update tempdist and escape_hist
      IF(minima_known) THEN
         CALL update_tempdist(this, start_minima%tempdist, tempstep, -1)
      ELSE
         CALL update_tempdist(this, start_minima%tempdist, tempstep, +1)
         start_minima%escape_hist(:) = EOSHIFT(start_minima%escape_hist, 1)
         start_minima%escape_hist(1) = report_Epot
      ENDIF

    ENDIF

    IF(.NOT. minima_known) THEN
       this%n_minima = this%n_minima + 1
       IF(this%n_minima > SIZE(this%minimas)) THEN
          minimas_tmp(:) = this%minimas
          DEALLOCATE(this%minimas)
          ALLOCATE(this%minimas(SIZE(minimas_tmp) + 1000))
          this%minimas(:SIZE(minimas_tmp)) = minimas_tmp
       ENDIF

       new_mid = this%n_minima
       ALLOCATE(this%minimas(new_mid)%p)
       new_minima => this%minimas(new_mid)%p
       new_minima%id = new_mid
       ALLOCATE(new_minima%escape_hist(this%esc_hist_len))
       ALLOCATE(new_minima%tempdist(this%tempstep_max))

       new_minima%escape_hist(:) = report_Epot !init with Epot

       IF(ASSOCIATED(start_minima)) THEN
          new_minima%tempdist(:) = start_minima%tempdist ! inherit tempdist
       ELSE
          new_minima%tempdist(:) = this%tempdist_init
       ENDIF

       new_minima%Epot = report_Epot
       new_minima%fp   = report_fp
       ALLOCATE(new_minima%pos(SIZE(report_positions)))
       new_minima%pos(:)  = report_positions

       IF(ASSOCIATED(start_minima))THEN
         IF(report_Epot < start_minima%Epot) THEN
           start_minima%disabled = .TRUE.
           IF(this%iw > 0) WRITE(this%iw,'(1X,A,T71,I10)') &
             "MINCRAWL| Disabling minimum with id", start_minima%id
         ENDIF
       ENDIF

       IF(this%iw > 0) WRITE(this%iw,'(1X,A,T71,I10)') &
          "MINCRAWL| Adding new minima with id", new_mid

       CALL history_add(this%history, report_fp, id=new_mid)
       CALL write_minima_traj(this, wid, new_mid, report_Epot, report_positions)
    ENDIF
    DEALLOCATE(report_positions)
 END SUBROUTINE mincrawl_register_minima


! *****************************************************************************
!> \brief Helper routine for mincrawl_register_minima.
!>        Adds or substracts small Gaussian from a minimum's temp-distribution.
!> \param this ...
!> \param tempdist ...
!> \param center ...
!> \param direction ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE update_tempdist(this, tempdist, center, direction)
    TYPE(mincrawl_type)                      :: this
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: tempdist
    INTEGER                                  :: center, direction

    INTEGER                                  :: i

   DO i=1, SIZE(tempdist)
     tempdist(i) = tempdist(i) + this%tempdist_update_height &
          * REAL(direction,KIND=dp) * EXP(-((center-i)/this%tempdist_update_width)**2)
     tempdist(i) = MAX(0.0_dp, MIN(1.0_dp, tempdist(i)))
   ENDDO
 END SUBROUTINE update_tempdist


! *****************************************************************************
!> \brief Helper routine for mincrawl_register_minima, write trajectory.
!> \param this ...
!> \param worker_id ...
!> \param minimum_id ...
!> \param Epot ...
!> \param positions ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE write_minima_traj(this, worker_id, minimum_id, Epot, positions)
    TYPE(mincrawl_type), INTENT(INOUT)       :: this
    INTEGER, INTENT(IN)                      :: worker_id, minimum_id
    REAL(KIND=dp), INTENT(IN)                :: Epot
    REAL(KIND=dp), DIMENSION(:), POINTER     :: positions

    CHARACTER(len=default_string_length)     :: title, unit_str
    REAL(KIND=dp)                            :: unit_conv

    IF(this%minima_traj_unit <= 0) RETURN

    WRITE (title, '(A,I8,A,I5,A,F20.10)') 'minimum_id = ',minimum_id, &
      ' worker_id = ',worker_id, ' Epot = ', Epot

    !get the conversion factor for the length unit
    CALL section_vals_val_get(this%mincrawl_section,"MINIMA_TRAJECTORY%UNIT",&
            c_val=unit_str)
    unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))

    CALL write_particle_coordinates(this%particle_set,&
                                    iunit=this%minima_traj_unit, &
                                    output_format=dump_xmol, &
                                    content="POS", &
                                    title=TRIM(title), &
                                    array=positions,&
                                    unit_conv=unit_conv)
 END SUBROUTINE write_minima_traj


! *****************************************************************************
!> \brief Finalizes master for Minima Crawling
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE mincrawl_finalize(this)
    TYPE(mincrawl_type)                      :: this

    INTEGER                                  :: i
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)

    DO i=1, this%n_minima
       !WRITE (*,*) "Minima: ", i, " n_sampled: ",this%minimas(i)%n_sampled
       DEALLOCATE(this%minimas(i)%p)
    ENDDO

    logger => cp_get_default_logger()
    CALL cp_print_key_finished_output(this%minima_traj_unit, logger, &
      this%mincrawl_section, "MINIMA_TRAJECTORY")

    CALL history_finalize(this%history)
    CALL delete_rng_stream(this%rng_stream)
 END SUBROUTINE mincrawl_finalize


END MODULE glbopt_mincrawl

