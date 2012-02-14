!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for SCPTB 
!> \author JGH (12.2011)
! *****************************************************************************
MODULE scptb_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                             create_gto_from_sto_basis,&
                                             deallocate_sto_basis_set,&
                                             set_sto_basis_set,&
                                             sto_basis_set_type
  USE cp_control_types,                ONLY: scptb_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE erf_fn,                          ONLY: erfc
  USE external_potential_types,        ONLY: set_potential
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE scptb_parameters,                ONLY: scptb_default_parameter
  USE scptb_types,                     ONLY: allocate_scptb_parameter,&
                                             scptb_parameter_type,&
                                             set_scptb_parameter
  USE string_utilities,                ONLY: uppercase
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_utils'

  PUBLIC :: scptb_parameter_init

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE scptb_parameter_init(atomic_kind_set,scptb_control,print_section,para_env,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(scptb_control_type), INTENT(IN)     :: scptb_control
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scptb_parameter_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: symbol
    CHARACTER(LEN=default_string_length)     :: iname
    INTEGER                                  :: i, ikind, is, l, nkind, &
                                                nshell, output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: lq, nq
    LOGICAL                                  :: failure, print_info
    REAL(KIND=dp)                            :: a, b, c, d, eff, r, rmax, rmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    TYPE(sto_basis_set_type), POINTER        :: sto_basis

    failure = .FALSE.
    output_unit = -1
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    print_info = (BTEST(cp_print_key_should_output(logger%iter_info,print_section,&
         "KINDS/BASIS_SET",error=error),cp_p_file)) 
    IF (print_info) THEN
       output_unit = cp_print_key_unit_nr(logger,print_section,"KINDS",extension=".Log",error=error)
       IF ( output_unit > 0 ) THEN
         WRITE(output_unit,"(/,A)") " SCPTB| A set of SCPTB "//&
              "parameters for material sciences."
         WRITE(output_unit,"(A)") " SCPTB| J. Hutter, Y. Misteli, R. Koitz"
         WRITE(output_unit,"(A)") " SCPTB| University of Zurich, 2012"
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,print_section,"KINDS",error=error)
    END IF

    nkind  = SIZE(atomic_kind_set)
    DO ikind = 1, nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,name=iname)
      CALL uppercase(iname)
      NULLIFY(scptb_parameter)
      CALL allocate_scptb_parameter(scptb_parameter,error)
      CALL set_scptb_parameter(scptb_parameter,key=iname,error=error)

      CALL scptb_parameter_from_file(scptb_parameter,scptb_control,para_env,error)
      CALL scptb_default_parameter(scptb_parameter,error)
      IF (scptb_parameter%defined) THEN
         CALL set_potential(potential=atomic_kind%all_potential,zeff=scptb_parameter%zeff,&
              zeff_correction=0.0_dp)
         ! basis set
         NULLIFY(sto_basis)
         CALL allocate_sto_basis_set (sto_basis,error)
         nshell = SUM(scptb_parameter%norb)
         ALLOCATE (nq(nshell),lq(nshell),zet(nshell),STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
         nq=0
         lq=0
         zet=0._dp
         ALLOCATE (symbol(nshell),STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
         symbol=""
         is=0
         DO l=0,scptb_parameter%lmaxorb
            DO i=1,scptb_parameter%norb(l)
               is = is + 1
               zet(is) = scptb_parameter%zeta(i,l)
               lq(is) = l
               nq(is) = scptb_parameter%nqm(i,l)
               IF(l==0) WRITE(symbol(is),'(I1,A1)') nq(is),"S"
               IF(l==1) WRITE(symbol(is),'(I1,A1)') nq(is),"P"
               IF(l==2) WRITE(symbol(is),'(I1,A1)') nq(is),"D"
               IF(l==3) WRITE(symbol(is),'(I1,A1)') nq(is),"F"
            END DO
         END DO

         IF (nshell > 0) THEN
            CALL set_sto_basis_set(sto_basis,name=scptb_parameter%atomname,&
                 nshell=nshell,symbol=symbol,nq=nq,lq=lq,zet=zet)
            CALL create_gto_from_sto_basis(sto_basis,atomic_kind%orb_basis_set,&
                 scptb_control%sto_ng,error=error)
         END IF

         DEALLOCATE (nq,lq,zet,symbol,STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
         CALL deallocate_sto_basis_set (sto_basis,error)

         !set interaction radius
         a = SQRT(0.5_dp*scptb_parameter%a0)
         rmin = 1._dp
         rmax = 100._dp
         DO
            r = 0.5_dp*(rmax+rmin)
            eff = erfc(a*r)/r
            IF (eff < scptb_control%epspair) THEN
               rmax = r
            ELSE
               rmin = r
            END IF
            IF ((rmax-rmin) < 1.e-2) EXIT
         END DO
         scptb_parameter%repair = rmax

         a = SQRT(scptb_parameter%crep(1))
         b = SQRT(scptb_parameter%crep(2))
         c = SQRT(scptb_parameter%crep(3))
         d = SQRT(scptb_parameter%crep(4))
         rmin = 1._dp
         rmax = 100._dp
         DO
            r = 0.5_dp*(rmax+rmin)
            eff = (c+d*r)*EXP(-(a*r+b*r*r))
            IF (eff < scptb_control%epspair) THEN
               rmax = r
            ELSE
               rmin = r
            END IF
            IF ((rmax-rmin) < 1.e-2) EXIT
         END DO
         scptb_parameter%rcpair = rmax

      END IF

      atomic_kind%scptb_parameter => scptb_parameter

    END DO

  END SUBROUTINE scptb_parameter_init

  SUBROUTINE scptb_parameter_from_file(scptb_parameter,scptb_control,para_env,error)

    TYPE(scptb_parameter_type), &
      INTENT(INOUT)                          :: scptb_parameter
    TYPE(scptb_control_type), INTENT(IN)     :: scptb_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scptb_parameter_from_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: fname
    LOGICAL                                  :: at_end
    TYPE(cp_parser_type), POINTER            :: parser

    fname = scptb_control%parameter_file
    IF ( fname /= "" ) THEN
       CALL parser_create(parser,fname,para_env=para_env,error=error)
       DO
         at_end = .FALSE.
         CALL parser_get_next_line(parser,1,at_end,error=error)
         IF ( at_end ) EXIT
!           CALL parser_get_object(parser,name_a,lower_to_upper=.TRUE.,error=error)
!           CALL parser_get_object(parser,name_b,lower_to_upper=.TRUE.,error=error)
!           !Checking Names
!           IF ( (iname==name_a .AND. jname==name_b) ) THEN
!             CALL parser_get_object(parser,skfn,string_length=8,error=error)
!             sk_files(ikind,jkind) = TRIM(dftb_control%sk_file_path)//"/"//&
!                                     TRIM(skfn)
!             found = .TRUE.
!             EXIT
!           END IF
!           !Checking Element
!           IF ( (iel==name_a .AND. jel==name_b) ) THEN
!             CALL parser_get_object(parser,skfn,string_length=8,error=error)
!             sk_files(ikind,jkind) = TRIM(dftb_control%sk_file_path)//"/"//&
!                                     TRIM(skfn)
!             found = .TRUE.
!             EXIT
!           END IF
       END DO
       CALL parser_release(parser,error=error)
    END IF

  END SUBROUTINE scptb_parameter_from_file

END MODULE scptb_utils

