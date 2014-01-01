!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_master
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE glbopt_mincrawl,                 ONLY: mincrawl_finalize,&
                                             mincrawl_init,&
                                             mincrawl_steer,&
                                             mincrawl_type
  USE glbopt_minhop,                   ONLY: minhop_finalize,&
                                             minhop_init,&
                                             minhop_steer,&
                                             minhop_type
  USE input_constants,                 ONLY: glbopt_do_mincrawl,&
                                             glbopt_do_minhop
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE swarm_message,                   ONLY: swarm_message_get,&
                                             swarm_message_type
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_master'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_steer


 TYPE glbopt_master_type
   PRIVATE
   REAL(KIND=dp)                                       :: E_lowest = HUGE(1.0_dp)
   REAL(KIND=dp)                                       :: E_target = TINY(1.0_dp)
   INTEGER                                             :: iw = 0
   INTEGER(int_8)                                      :: total_md_steps = 0
   INTEGER(int_8)                                      :: total_gopt_steps = 0
   INTEGER(int_8)                                      :: count_reports = 0
   INTEGER                                             :: method
   TYPE(minhop_type),POINTER                           :: minhop
   TYPE(mincrawl_type), POINTER                        :: mincrawl
   INTEGER                                             :: i_iteration = 0
 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(this, para_env, swarm_section, input_path, n_walkers, iw, error)
    TYPE(glbopt_master_type)                 :: this
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: swarm_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(section_vals_type), POINTER         :: glbopt_section

    this%iw = iw

    glbopt_section => section_vals_get_subs_vals(swarm_section, "GLOBAL_OPT", error=error)
    CALL section_vals_val_get(glbopt_section,"E_TARGET", r_val=this%E_target, error=error)
    CALL section_vals_val_get(glbopt_section,"METHOD", i_val=this%method, error=error)

    SELECT CASE (this%method)
       CASE(glbopt_do_minhop)
         ALLOCATE(this%minhop)
         CALL minhop_init(this%minhop, glbopt_section, n_walkers, iw, error)
       CASE(glbopt_do_mincrawl)
         ALLOCATE(this%mincrawl)
         CALL mincrawl_init(this%mincrawl, glbopt_section, n_walkers, iw, error)
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT
 END SUBROUTINE glbopt_master_init


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(this, report, cmd, should_stop)
    TYPE(glbopt_master_type)                 :: this
    TYPE(swarm_message_type)                 :: report, cmd
    LOGICAL, INTENT(INOUT)                   :: should_stop

    CHARACTER(len=default_string_length)     :: status
    INTEGER                                  :: gopt_steps, md_steps, &
                                                report_worker_id
    REAL(KIND=dp)                            :: report_Epot

    this%i_iteration = this%i_iteration + 1

    CALL swarm_message_get(report, "worker_id", report_worker_id)
    CALL swarm_message_get(report, "status", status)

    IF(TRIM(status) == "ok") THEN
       CALL swarm_message_get(report, "Epot", report_Epot)
       CALL swarm_message_get(report, "md_steps", md_steps)
       CALL swarm_message_get(report, "gopt_steps", gopt_steps)
       this%E_lowest         = MIN(this%E_lowest, report_Epot)
       this%total_md_steps   = this%total_md_steps + md_steps
       this%total_gopt_steps = this%total_gopt_steps + gopt_steps
       this%count_reports    = this%count_reports + 1

       IF(this%iw>0) THEN
          WRITE (this%iw,'(A,46X,I8)')&
             " GLBOPT| Reporting worker ", report_worker_id
          WRITE (this%iw,'(A,20X,E15.8)')&
             " GLBOPT| Reported potential energy [Hartree] ", report_Epot
          WRITE (this%iw,'(A,13X,E15.8)')&
             " GLBOPT| Lowest reported potential energy [Hartree] ", this%E_lowest
          WRITE (this%iw,'(A,T71,F10.1)')&
             " GLBOPT| Average number of MD steps", REAL(this%total_md_steps,KIND=dp)/this%count_reports
          WRITE (this%iw,'(A,T71,F10.1)')&
             " GLBOPT| Average number of Geo-Opt steps", REAL(this%total_gopt_steps,KIND=dp)/this%count_reports
       END IF
    END IF

    IF(this%E_lowest<this%E_target) THEN
       IF(this%iw>0)  WRITE (this%iw,"(A,I8,A)" )  &
          " GLBOPT| Reached E_pot < E_target after ", this%i_iteration," iterations. Quitting."
       should_stop = .TRUE.
    ELSE
       SELECT CASE (this%method)
         CASE(glbopt_do_minhop)
             CALL minhop_steer(this%minhop, report, cmd)
         CASE(glbopt_do_mincrawl)
             CALL mincrawl_steer(this%mincrawl, report, cmd)
         CASE DEFAULT
          STOP "Unkown glbopt_method"
      END SELECT
    ENDIF
 END SUBROUTINE glbopt_master_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(this)
    TYPE(glbopt_master_type)                 :: this

    SELECT CASE (this%method)
       CASE(glbopt_do_minhop)
         CALL minhop_finalize(this%minhop)
         DEALLOCATE(this%minhop)
       CASE(glbopt_do_mincrawl)
         CALL mincrawl_finalize(this%mincrawl)
         DEALLOCATE(this%mincrawl)
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT
 END SUBROUTINE glbopt_master_finalize


END MODULE glbopt_master

