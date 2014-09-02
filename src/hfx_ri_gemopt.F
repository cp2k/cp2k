!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief optimizer for geminal basis sets
!> \author JGH October 2009
! *****************************************************************************
MODULE hfx_ri_gemopt
  USE basis_set_types,                 ONLY: geminal_basis_set_type,&
                                             get_geminal_basis_set,&
                                             write_geminal_basis_set
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hfx_ri_methods,                  ONLY: hfx_ri_energy,&
                                             hfx_ri_env_create
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast
  USE powell,                          ONLY: opt_state_type,&
                                             powell_optimize
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE ri_environment_types,            ONLY: get_ri_env,&
                                             ri_environment_type,&
                                             ri_vector_set,&
                                             ri_vector_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_gemopt'

  PUBLIC :: geminal_optimize

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief optimizer for geminal bsis sets
!>
!> \param qs_env ...
!> \param print_key ...
!> \param iunit ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE geminal_optimize (qs_env,print_key,iunit,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: print_key
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geminal_optimize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, n10, n100
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ehf
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(opt_state_type)                     :: ostate
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff
    TYPE(section_vals_type), POINTER         :: powell_section

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,energy=energy,para_env=para_env,error=error)
    ehf = energy%ex
    CALL hfx_ri_env_create (qs_env,error)

    NULLIFY(ri_env)
    CALL get_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,error=error)
    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,error=error)

    powell_section => section_vals_get_subs_vals(section_vals=print_key,&
                   subsection_name="POWELL",error=error)

    CALL get_geminfo(qs_env,x,error)
    ostate%nf = 0
    ostate%nvar = SIZE(x)
    CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
    CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
    CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)
    ostate%iprint = 1
    ostate%unit  = iunit

    ostate%state = 0
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")')
      WRITE(iunit,'(" POWELL| Number of parameters in optimization",T71,I10)') ostate%nvar
    END IF
    n10 = MAX(ostate%maxfun/10,1)
    n100 = MAX(ostate%maxfun/100,1)

    DO

      IF ( ostate%state == 2 ) THEN
        CALL set_geminfo (qs_env,x,error)
        CALL ri_vector_set(ri_coeff, 0._dp, error)
        CALL hfx_ri_energy (qs_env,error)
        ostate%f = energy%ex - ehf
      END IF

      IF ( ostate%state == -1 ) EXIT

      IF ( para_env%source == para_env%mepos ) THEN
         CALL powell_optimize (ostate%nvar, x, ostate)
      END IF

      CALL mp_bcast(x,para_env%source,para_env%group)
      CALL mp_bcast(ostate%state,para_env%source,para_env%group)
      CALL mp_bcast(ostate%fopt,para_env%source,para_env%group)
      CALL mp_bcast(ostate%nf,para_env%source,para_env%group)

      IF ( MOD(ostate%nf,n10) == 0 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
              INT(REAL(ostate%nf,dp)/REAL(ostate%maxfun,dp)*100._dp),ostate%fopt
        WRITE(iunit,'(T9,3F24.12)') x
      END IF
      IF ( MOD(ostate%nf,n100) == 0 .AND. iunit > 0 ) THEN
        CALL write_gemprog (ostate,x,error)
      END IF

    END DO

    ostate%state = 8
    IF ( para_env%source == para_env%mepos ) CALL powell_optimize (ostate%nvar, x, ostate)
    CALL mp_bcast(x,para_env%source,para_env%group)
    CALL set_geminfo (qs_env,x,error)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,G20.10)') ostate%fopt
      CALL write_geminfo (qs_env,iunit,error)
    END IF

    IF ( ASSOCIATED(x) ) THEN
       DEALLOCATE (x,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE geminal_optimize

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param x ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_geminfo (qs_env,x,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_geminfo', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: restriction
    INTEGER                                  :: ikind, iset, istat, nkind, &
                                                nopt, np, nset
    INTEGER, DIMENSION(:), POINTER           :: npgf
    LOGICAL                                  :: failure, noopt
    REAL(KIND=dp)                            :: a, b, c
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: zet
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    failure = .FALSE.
    NULLIFY(qs_kind_set)

    CALL get_qs_env(qs_env=qs_env,qs_kind_set=qs_kind_set,error=error)
    nkind = SIZE(qs_kind_set)

    nopt = 0
    DO ikind = 1,nkind
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind_set(ikind),no_optimize=noopt,&
                            geminal_basis_set=geminal_basis_set)
       IF ( noopt ) CYCLE
       IF ( .NOT.ASSOCIATED(geminal_basis_set) ) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,type_restriction=restriction,&
                                  nset=nset,npgf=npgf)
       DO iset = 1, nset
         np = npgf(iset)
         CPPrecondition(np==1,cp_failure_level,routineP,error,failure)
       END DO
       IF ( restriction=="NO" ) THEN
         nopt = nopt + 3*nset
       ELSE
         nopt = nopt + 2*nset
       END IF
    END DO
    ALLOCATE (x(nopt),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    nopt = 0
    DO ikind = 1,nkind
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind_set(ikind),no_optimize=noopt,&
                            geminal_basis_set=geminal_basis_set)
       IF ( noopt ) CYCLE
       IF ( .NOT.ASSOCIATED(geminal_basis_set) ) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,nset=nset,zet=zet)
       DO iset = 1, nset
         a = SQRT(zet(1,1,1,iset))
         c = zet(2,1,1,iset)/a
         b = SQRT(zet(2,2,1,iset)-c*c)
         SELECT CASE (restriction)
           CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
           CASE ("NO")
             nopt = nopt+1
             x(nopt) = a
             nopt = nopt+1
             x(nopt) = b
             nopt = nopt+1
             x(nopt) = c
           CASE ("RS")
             CPPrecondition(ABS(c) < 1.E-6,cp_failure_level,routineP,error,failure)
             nopt = nopt+1
             x(nopt) = a
             nopt = nopt+1
             x(nopt) = b
           CASE ("RR")
             CPPrecondition(ABS(a*a-b*b-c*c) < 1.E-6,cp_failure_level,routineP,error,failure)
             nopt = nopt+1
             x(nopt) = a
             nopt = nopt+1
             x(nopt) = b
         END SELECT
       END DO
    END DO

  END SUBROUTINE get_geminfo
! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param x ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_geminfo (qs_env,x,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_geminfo', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: restriction
    INTEGER                                  :: ikind, iset, nkind, nopt, nset
    LOGICAL                                  :: failure, noopt
    REAL(KIND=dp)                            :: a, b, c
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: zet, zeth
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    failure = .FALSE.
    NULLIFY(qs_kind_set)

    CALL get_qs_env(qs_env=qs_env,qs_kind_set=qs_kind_set,error=error)
    nkind = SIZE(qs_kind_set)

    nopt = 0
    DO ikind = 1,nkind
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind_set(ikind),no_optimize=noopt,&
                            geminal_basis_set=geminal_basis_set)
       IF ( noopt ) CYCLE
       IF ( .NOT.ASSOCIATED(geminal_basis_set) ) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,type_restriction=restriction,&
            nset=nset,zet=zet,zeth=zeth)
       DO iset = 1, nset
         SELECT CASE (restriction)
           CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
           CASE ("NO")
             nopt = nopt + 1
             a = x(nopt)
             nopt = nopt + 1
             b = x(nopt)
             nopt = nopt + 1
             c = x(nopt)
           CASE ("RS")
             nopt = nopt + 1
             a = x(nopt)
             nopt = nopt + 1
             b = x(nopt)
             c = 0._dp
           CASE ("RR")
             nopt = nopt + 1
             a = x(nopt)
             nopt = nopt + 1
             b = x(nopt)
             c = SQRT(ABS(a*a-b*b))
         END SELECT

         zet(1,1,1,iset) = a*a
         zet(2,1,1,iset) = a*c
         zet(1,2,1,iset) = a*c
         zet(2,2,1,iset) = b*b+c*c
         zeth(1,1,1,iset) = a*a
         zeth(2,1,1,iset) = -a*c
         zeth(1,2,1,iset) = -a*c
         zeth(2,2,1,iset) = b*b+c*c
       END DO

!       kind_radius = 0.0_dp
!       DO iset=1,nset
!         set_radius(iset) = 0.0_dp
!         DO ipgf=1,npgf(iset)
!           pgf_radius(ipgf,iset) = 0.0_dp
!           DO ishell=1,nshell(iset)
!             l = MAX(lshell(ishell,iset),ls(iset))
!             gcca = gcc(ipgf,ishell,iset)
!             zab = zetg(1,1,ipgf,iset) + zetg(2,2,ipgf,iset)
!             zd  = (zetg(1,1,ipgf,iset) - zetg(2,2,ipgf,iset))**2 + 4._dp*zetg(1,2,ipgf,iset)
!             zeta = 0.5_dp*MIN(zab+SQRT(zd),zab-SQRT(zd))
!             pgf_radius(ipgf,iset) = MAX(pgf_radius(ipgf,iset),&
!                                         exp_radius(l,zeta,qs_control%eps_pgf_orb,gcca))
!           END DO
!           set_radius(iset) = MAX(set_radius(iset),pgf_radius(ipgf,iset))
!         END DO
!         kind_radius = MAX(kind_radius,set_radius(iset))
!       END DO
!       CALL set_geminal_basis_set(geminal_basis_set=geminal_basis_set,&
!                                  pgf_radius=pgf_radius,&
!                                  set_radius=set_radius,&
!                                  kind_radius=kind_radius)

    END DO

  END SUBROUTINE set_geminfo
! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param iunit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_geminfo (qs_env,iunit,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_geminfo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind
    LOGICAL                                  :: noopt
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    NULLIFY(qs_kind_set)

    CALL get_qs_env(qs_env=qs_env,qs_kind_set=qs_kind_set,error=error)
    nkind = SIZE(qs_kind_set)

    DO ikind = 1,nkind
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind_set(ikind),no_optimize=noopt,&
                            geminal_basis_set=geminal_basis_set)
       IF ( noopt ) CYCLE
       IF ( .NOT.ASSOCIATED(geminal_basis_set) ) CYCLE

       CALL write_geminal_basis_set(geminal_basis_set,iunit,error)
    END DO

  END SUBROUTINE write_geminfo
! *****************************************************************************
!> \brief ...
!> \param ostate ...
!> \param x ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_gemprog (ostate,x,error)
    TYPE(opt_state_type), INTENT(IN)         :: ostate
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_gemprog', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iunit

    CALL open_file(file_name="powell.res",file_status="UNKNOWN",&
                   file_action="READWRITE",file_position="APPEND",unit_number=iunit)

    WRITE(iunit,*) ostate%nf, ostate%fopt, ostate%xopt

    CALL close_file(iunit)

  END SUBROUTINE write_gemprog
! *****************************************************************************

END MODULE hfx_ri_gemopt

