!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the Minima Hopping global optimization scheme
!> \author Ole Schuett
! *****************************************************************************
MODULE glbopt_minhop
  USE bibliography,                    ONLY: Goedecker2004,&
                                             cite_reference
  USE glbopt_history,                  ONLY: history_add,&
                                             history_finalize,&
                                             history_fingerprint,&
                                             history_fingerprint_match,&
                                             history_fingerprint_type,&
                                             history_init,&
                                             history_lookup,&
                                             history_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE physcon,                         ONLY: kelvin
  USE swarm_message,                   ONLY: swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_type
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_minhop'


 PUBLIC :: minhop_type
 PUBLIC :: minhop_init, minhop_finalize
 PUBLIC :: minhop_steer

 TYPE worker_state_type
   REAL(KIND=dp)                                       :: Eaccept = -1.0
   REAL(KIND=dp)                                       :: temp = -1.0
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pos
   REAL(KIND=dp)                                       :: Epot = -1.0
   TYPE(history_fingerprint_type)                      :: fp
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pos_hop
   REAL(KIND=dp)                                       :: Epot_hop = HUGE(1.0)
   TYPE(history_fingerprint_type)                      :: fp_hop
   INTEGER                                             :: minima_id = -1
   INTEGER                                             :: iframe = 1
 END TYPE worker_state_type

 TYPE minima_state_type
   REAL(KIND=dp)                                       :: Eaccept = -1.0
   REAL(KIND=dp)                                       :: temp = -1.0
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pos
   REAL(KIND=dp)                                       :: Epot = -1.0
   TYPE(history_fingerprint_type)                      :: fp
   LOGICAL                                             :: disabled = .FALSE.
   INTEGER                                             :: n_active = 0
   INTEGER                                             :: n_sampled = 0
 END TYPE minima_state_type

 TYPE minhop_type
   PRIVATE
   TYPE(history_type), DIMENSION(:),ALLOCATABLE        :: history
   TYPE(worker_state_type), DIMENSION(:), ALLOCATABLE  :: worker_state
   TYPE(minima_state_type), DIMENSION(:), ALLOCATABLE  :: minima_state
   INTEGER                                             :: n_minima = 0
   REAL(KIND=dp)                                       :: beta1 = 0
   REAL(KIND=dp)                                       :: beta2 = 0
   REAL(KIND=dp)                                       :: beta3 = 0
   REAL(KIND=dp)                                       :: Eaccept0 = 0
   REAL(KIND=dp)                                       :: temp_init = 0
   REAL(KIND=dp)                                       :: temp_max = 0
   REAL(KIND=dp)                                       :: temp_min = 0
   REAL(KIND=dp)                                       :: alpha1 = 0
   REAL(KIND=dp)                                       :: alpha2 = 0
   INTEGER                                             :: n_accepted = 0
   INTEGER                                             :: n_rejected = 0
   INTEGER                                             :: iw = 0
   INTEGER                                             :: n_workers =0
   LOGICAL                                             :: share_history = .FALSE.
 END TYPE minhop_type

 CONTAINS


! *****************************************************************************
!> \brief Initializes master for Minima Hopping
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE minhop_init(this, glbopt_section, n_workers, iw, error)
    TYPE(minhop_type)                        :: this
    TYPE(section_vals_type), POINTER         :: glbopt_section
    INTEGER, INTENT(IN)                      :: n_workers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, n_histories
    REAL(kind=dp)                            :: temp_in_kelvin
    TYPE(section_vals_type), POINTER         :: history_section, &
                                                minhop_section

    CALL cite_reference(Goedecker2004)

    ! read input
    minhop_section => section_vals_get_subs_vals(glbopt_section, "MINIMA_HOPPING", error=error)
    CALL section_vals_val_get(minhop_section,"BETA_1",  r_val=this%beta1,error=error)
    CALL section_vals_val_get(minhop_section,"BETA_2",  r_val=this%beta2,error=error)
    CALL section_vals_val_get(minhop_section,"BETA_3",  r_val=this%beta3,error=error)
    CALL section_vals_val_get(minhop_section,"ALPHA_1", r_val=this%alpha1,error=error)
    CALL section_vals_val_get(minhop_section,"ALPHA_2", r_val=this%alpha2,error=error)
    CALL section_vals_val_get(minhop_section,"E_ACCEPT_INIT", r_val=this%Eaccept0,error=error)
    CALL section_vals_val_get(minhop_section,"TEMPERATURE_INIT", r_val=temp_in_kelvin,error=error)
    this%temp_init = temp_in_kelvin / kelvin
    CALL section_vals_val_get(minhop_section,"SHARE_HISTORY", l_val=this%share_history,error=error)

    ! allocate history / histories
    history_section => section_vals_get_subs_vals(glbopt_section, "HISTORY", error=error)
    n_histories = n_workers
    IF(this%share_history) n_histories = 1
    ALLOCATE(this%history(n_histories))

    !only the first history shall write to iw
    CALL history_init(this%history(1), history_section, iw=iw, error=error)
    DO i=2, n_histories
      CALL history_init(this%history(i), history_section, iw=-1, error=error)
    END DO

    ALLOCATE(this%worker_state(n_workers))
    this%n_workers = n_workers
    this%iw = iw

    IF(this%iw>0) THEN
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| beta_1", this%beta1
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| beta_2", this%beta2
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| beta_3", this%beta3
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| alpha_1", this%alpha1
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| alpha_2", this%alpha2
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| Initial acceptance energy [Hartree]", this%Eaccept0
      WRITE(this%iw, '(A,T71,F10.3)') " MINHOP| Initial temperature [Kelvin]", this%temp_init * kelvin
      WRITE(this%iw, '(A,T71,L10)') " MINHOP| All workers share a single history", this%share_history
    ENDIF
 END SUBROUTINE minhop_init


! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE minhop_steer(this, report, cmd)
    TYPE(minhop_type)                        :: this
    TYPE(swarm_message_type)                 :: report, cmd

    CHARACTER(len=default_string_length)     :: status
    INTEGER                                  :: hid, iframe, wid
    LOGICAL                                  :: minima_known
    REAL(KIND=dp)                            :: report_Epot
    REAL(KIND=dp), DIMENSION(:), POINTER     :: report_positions
    TYPE(history_fingerprint_type)           :: report_fp

    NULLIFY(report_positions)
    CALL swarm_message_get(report, "worker_id", wid)
    CALL swarm_message_get(report, "status", status)

    IF(TRIM(status) == "initial_hello") THEN
       this%worker_state(wid)%temp = this%temp_init
       this%worker_state(wid)%Eaccept = this%Eaccept0
       CALL swarm_message_add(cmd, "command", "md_and_gopt")
       CALL swarm_message_add(cmd, "iframe", 1)
       CALL swarm_message_add(cmd, "temperature", this%worker_state(wid)%temp)
       IF(this%iw > 0) WRITE(this%iw,'(1X,A,1X,I10,1X,A,7X,F10.3)') &
          "MINHOP| Sending worker", wid, &
          "initial temperature [Kelvin]",this%worker_state(wid)%temp * kelvin
       RETURN
    ENDIF

    hid = wid ! history_id = worker_id unless ....
    IF(this%share_history) hid = 1 !...there is ONE shared history.

    CALL swarm_message_get(report, "Epot", report_Epot)
    CALL swarm_message_get(report, "positions", report_positions)

    report_fp = history_fingerprint(this%history(hid), report_Epot, report_positions)

    IF(.NOT. ALLOCATED(this%worker_state(wid)%pos)) THEN
       !init (first real report)
       this%worker_state(wid)%Epot        = report_Epot
       ALLOCATE(this%worker_state(wid)%pos(SIZE(report_positions)))
       this%worker_state(wid)%pos(:)         = report_positions
       this%worker_state(wid)%fp          = report_fp
    END IF

    IF(history_fingerprint_match(this%history(hid), this%worker_state(wid)%fp, report_fp)) THEN
       ! not escaped
       IF(this%iw>0) WRITE(this%iw,'(A)') " MINHOP| Not escaped"
       this%worker_state(wid)%temp = this%worker_state(wid)%temp * this%beta1 !increasing temperature
    ELSE
       ! escaped
       CALL history_lookup(this%history(hid), report_fp, minima_known)
       IF(minima_known) THEN
          IF(this%iw>0) WRITE(this%iw,'(A)') " MINHOP| Escaped, old minima"
          this%worker_state(wid)%temp = this%worker_state(wid)%temp * this%beta2 !increasing temperature
       ELSE
          IF(this%iw>0) WRITE(this%iw,'(A)') " MINHOP| Escaped, new minima"
          this%worker_state(wid)%temp = this%worker_state(wid)%temp * this%beta3 !decreasing temperature
          CALL history_add(this%history(hid), report_fp)
       ENDIF

       IF(report_Epot < this%worker_state(wid)%Epot_hop) THEN
          ! new locally lowest
          IF(this%iw>0) WRITE(this%iw,'(A)') " MINHOP| New locally lowest"
          this%worker_state(wid)%Epot_hop = report_Epot
          IF(.NOT. ALLOCATED(this%worker_state(wid)%pos_hop))&
             ALLOCATE(this%worker_state(wid)%pos_hop(SIZE(report_positions)))
          this%worker_state(wid)%pos_hop(:)  = report_positions
          this%worker_state(wid)%fp_hop   = report_fp
       ENDIF

       IF(this%worker_state(wid)%Epot_hop - this%worker_state(wid)%Epot < this%worker_state(wid)%Eaccept) THEN
          ! accept
          IF(this%iw>0) WRITE(this%iw,'(A)') " MINHOP| Accept"
          this%worker_state(wid)%Epot = this%worker_state(wid)%Epot_hop
          this%worker_state(wid)%pos(:)  = this%worker_state(wid)%pos_hop
          this%worker_state(wid)%fp   = this%worker_state(wid)%fp_hop
          this%worker_state(wid)%Epot_hop = HUGE(1.0)

          this%worker_state(wid)%Eaccept = this%worker_state(wid)%Eaccept * this%alpha1 !decreasing Eaccept
          this%n_accepted = this%n_accepted + 1
       ELSE
          ! not accept
          IF(this%iw>0) WRITE(this%iw,'(A)') " MINHOP| Reject"
          this%worker_state(wid)%Eaccept = this%worker_state(wid)%Eaccept * this%alpha2 !increasing Eaccept
          this%n_rejected = this%n_rejected + 1
       ENDIF
    END IF

    IF(this%iw>0) THEN
        WRITE (this%iw,'(A,15X,E20.10)')&
         " MINHOP| Worker's acceptance Energy [Hartree]",this%worker_state(wid)%Eaccept
        WRITE (this%iw,'(A,22X,F20.3)')&
         " MINHOP| Worker's temperature [Kelvin]",this%worker_state(wid)%temp * kelvin
    END IF

    CALL swarm_message_get(report, "iframe", iframe)
    CALL swarm_message_add(cmd, "iframe", iframe)
    CALL swarm_message_add(cmd, "command", "md_and_gopt")
    CALL swarm_message_add(cmd, "positions",  this%worker_state(wid)%pos)
    CALL swarm_message_add(cmd, "temperature", this%worker_state(wid)%temp)

    IF(this%iw > 0) THEN
      WRITE(this%iw,'(A,30X,I10)')&
       " MINHOP| Total number of accepted minima",this%n_accepted
      WRITE(this%iw,'(A,30X,I10)')&
       " MINHOP| Total number of rejected minima",this%n_rejected
    ENDIF

    DEALLOCATE(report_positions)
 END SUBROUTINE minhop_steer


! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE minhop_finalize(this)
    TYPE(minhop_type)                        :: this

    INTEGER                                  :: i

    DO i=1, SIZE(this%history)
      CALL history_finalize(this%history(i))
    END DO
 END SUBROUTINE minhop_finalize


END MODULE glbopt_minhop

