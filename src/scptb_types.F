!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the SCPTB parameter types.
!> \author JGH (12.2011)
! *****************************************************************************
MODULE scptb_types

  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_types'

! *****************************************************************************
  TYPE scptb_parameter_type
    ! PRIVATE
    CHARACTER(LEN=default_string_length)   :: key
    CHARACTER(LEN=default_string_length)   :: parameterization
    CHARACTER(LEN=default_string_length)   :: atomname
    LOGICAL                                :: defined
    INTEGER                                :: z          !atomic number
    REAL(KIND=dp)                          :: zeff       !effective core charge
    INTEGER                                :: natorb     !total number of orbitals
    INTEGER                                :: lmaxorb    !max angular momentum orbitals
    INTEGER                                :: lmaxscp    !max angular momentum scp density
    INTEGER, DIMENSION(0:3)                :: norb       !number of orbitals per l QN
    INTEGER, DIMENSION(10,0:3)             :: nqm        !principal quantum numbers
    REAL(KIND=dp), DIMENSION(10,0:3)       :: zeta       !orbital exponents
    REAL(KIND=dp), DIMENSION(10,0:3)       :: hcore      !core potential energy
    REAL(KIND=dp), DIMENSION(10,0:3)       :: occupation !free atom occupation
    REAL(KIND=dp)                          :: energy     !free atom energy
    REAL(KIND=dp), DIMENSION(3)            :: crep       !core repulsion parameters
    REAL(KIND=dp), DIMENSION(1:3)          :: pol        !polarization for SCP
    REAL(KIND=dp)                          :: ag         !width of SCP charge
    REAL(KIND=dp)                          :: rcpair     !cutoff radius for core pair potentials
  END TYPE scptb_parameter_type

  TYPE scptb_parameter_p_type
     TYPE(scptb_parameter_type), POINTER    :: scptb_param
  END TYPE scptb_parameter_p_type

! *****************************************************************************

  TYPE scp_kind_vector_type
     INTEGER                                           :: natom
     INTEGER                                           :: nbasis
     REAL(dp), DIMENSION(:,:), POINTER                 :: vmat
  END TYPE scp_kind_vector_type

  TYPE scp_vector_type
     INTEGER                                           :: total_length
     TYPE(scp_kind_vector_type), DIMENSION(:), POINTER :: vector
  END TYPE scp_vector_type

  TYPE scp_vector_p_type
     TYPE(scp_vector_type), POINTER    :: vec
  END TYPE scp_vector_p_type

! *****************************************************************************

  PUBLIC :: scptb_parameter_type, &
            scptb_parameter_p_type, &
            allocate_scptb_parameter,&
            deallocate_scptb_parameter,&
            get_scptb_parameter,&
            set_scptb_parameter,&
            write_scptb_parameter

  PUBLIC :: scp_kind_vector_type,&
            scp_vector_type,&
            scp_vector_p_type,&
            scp_vector_create,&
            scp_vector_release,&
            scp_vector_set,&
            scp_vector_add,&
            scp_vector_mult,&
            scp_vector_copy,&
            scp_vector_scale,&
            scp_vector_dot,&
            scp_vector_norm,&
            scp_vector_print,&
            scp_vector_sync

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE allocate_scptb_parameter(scptb_parameter,error)

    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_scptb_parameter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    IF (ASSOCIATED(scptb_parameter)) &
            CALL deallocate_scptb_parameter(scptb_parameter,error)

    ALLOCATE (scptb_parameter,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    scptb_parameter % key = ""
    scptb_parameter % parameterization = ""
    scptb_parameter % atomname = ""
    scptb_parameter % defined = .FALSE.
    scptb_parameter % z = -1
    scptb_parameter % zeff = -1.0_dp
    scptb_parameter % natorb = 0
    scptb_parameter % lmaxorb = -1
    scptb_parameter % lmaxscp = -1
    scptb_parameter % norb = 0
    scptb_parameter % nqm = 0
    scptb_parameter % zeta = 0.0_dp
    scptb_parameter % hcore = 0.0_dp
    scptb_parameter % occupation = 0.0_dp
    scptb_parameter % energy = 0.0_dp
    scptb_parameter % crep = 0.0_dp
    scptb_parameter % pol = 0.0_dp
    scptb_parameter % ag = 0.0_dp
    scptb_parameter % rcpair = 0.0_dp

  END SUBROUTINE allocate_scptb_parameter

! *****************************************************************************
  SUBROUTINE deallocate_scptb_parameter(scptb_parameter,error)

    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_scptb_parameter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(scptb_parameter),cp_failure_level,routineP,error,failure)
    DEALLOCATE (scptb_parameter,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_scptb_parameter

! *****************************************************************************
  SUBROUTINE get_scptb_parameter(scptb_parameter,key,atomname,parameterization,defined,z,zeff,natorb,&
    lmaxorb,lmaxscp,norb,nqm,zeta,hcore,occupation,energy,crep,pol,ag,rcpair,error)

    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: key, atomname, &
                                                parameterization
    LOGICAL, INTENT(OUT), OPTIONAL           :: defined
    INTEGER, INTENT(OUT), OPTIONAL           :: z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zeff
    INTEGER, INTENT(OUT), OPTIONAL           :: natorb, lmaxorb, lmaxscp
    INTEGER, DIMENSION(0:3), OPTIONAL        :: norb
    INTEGER, DIMENSION(10, 0:3), OPTIONAL    :: nqm
    REAL(KIND=dp), DIMENSION(10, 0:3), &
      OPTIONAL                               :: zeta, hcore, occupation
    REAL(KIND=dp), OPTIONAL                  :: energy
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: crep
    REAL(KIND=dp), DIMENSION(1:3), OPTIONAL  :: pol
    REAL(KIND=dp), OPTIONAL                  :: ag, rcpair
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_scptb_parameter', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(scptb_parameter),cp_failure_level,routineP,error,failure)

    IF (PRESENT(key)) key = scptb_parameter%key
    IF (PRESENT(atomname)) atomname = scptb_parameter%atomname
    IF (PRESENT(parameterization)) parameterization = scptb_parameter%parameterization
    IF (PRESENT(defined)) defined = scptb_parameter%defined
    IF (PRESENT(z)) z = scptb_parameter%z
    IF (PRESENT(zeff)) zeff = scptb_parameter%zeff
    IF (PRESENT(natorb)) natorb = scptb_parameter%natorb
    IF (PRESENT(lmaxorb)) lmaxorb = scptb_parameter%lmaxorb
    IF (PRESENT(lmaxscp)) lmaxscp = scptb_parameter%lmaxscp
    IF (PRESENT(norb)) norb = scptb_parameter%norb
    IF (PRESENT(nqm)) nqm = scptb_parameter%nqm
    IF (PRESENT(zeta)) zeta = scptb_parameter%zeta
    IF (PRESENT(hcore)) hcore = scptb_parameter%hcore
    IF (PRESENT(occupation)) occupation = scptb_parameter%occupation
    IF (PRESENT(energy)) energy = scptb_parameter%energy
    IF (PRESENT(crep)) crep = scptb_parameter%crep
    IF (PRESENT(pol)) pol = scptb_parameter%pol
    IF (PRESENT(ag)) ag = scptb_parameter%ag
    IF (PRESENT(rcpair)) rcpair = scptb_parameter%rcpair

  END SUBROUTINE get_scptb_parameter

! *****************************************************************************
  SUBROUTINE set_scptb_parameter(scptb_parameter,key,atomname,parameterization,defined,z,zeff,natorb,&
    lmaxorb,lmaxscp,norb,nqm,zeta,hcore,occupation,energy,crep,pol,ag,rcpair,error)

    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    CHARACTER(LEN=*), OPTIONAL               :: key, atomname, &
                                                parameterization
    LOGICAL, OPTIONAL                        :: defined
    INTEGER, OPTIONAL                        :: z
    REAL(KIND=dp), OPTIONAL                  :: zeff
    INTEGER, OPTIONAL                        :: natorb, lmaxorb, lmaxscp
    INTEGER, DIMENSION(0:3), OPTIONAL        :: norb
    INTEGER, DIMENSION(10, 0:3), OPTIONAL    :: nqm
    REAL(KIND=dp), DIMENSION(10, 0:3), &
      OPTIONAL                               :: zeta, hcore, occupation
    REAL(KIND=dp), OPTIONAL                  :: energy
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: crep
    REAL(KIND=dp), DIMENSION(1:3), OPTIONAL  :: pol
    REAL(KIND=dp), OPTIONAL                  :: ag, rcpair
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_scptb_parameter', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(scptb_parameter),cp_failure_level,routineP,error,failure)

    IF (PRESENT(key)) scptb_parameter%key = key
    IF (PRESENT(atomname)) scptb_parameter%atomname = atomname
    IF (PRESENT(parameterization)) scptb_parameter%parameterization = parameterization
    IF (PRESENT(defined)) scptb_parameter%defined = defined
    IF (PRESENT(z)) scptb_parameter%z = z
    IF (PRESENT(zeff)) scptb_parameter%zeff = zeff
    IF (PRESENT(natorb)) scptb_parameter%natorb = natorb
    IF (PRESENT(lmaxorb)) scptb_parameter%lmaxorb = lmaxorb
    IF (PRESENT(lmaxscp)) scptb_parameter%lmaxscp = lmaxscp
    IF (PRESENT(norb)) scptb_parameter%norb = norb
    IF (PRESENT(nqm)) scptb_parameter%nqm = nqm
    IF (PRESENT(zeta)) scptb_parameter%zeta = zeta
    IF (PRESENT(hcore)) scptb_parameter%hcore = hcore
    IF (PRESENT(occupation)) scptb_parameter%occupation = occupation
    IF (PRESENT(energy)) scptb_parameter%energy = energy
    IF (PRESENT(crep)) scptb_parameter%crep = crep
    IF (PRESENT(pol)) scptb_parameter%pol = pol
    IF (PRESENT(ag)) scptb_parameter%ag = ag
    IF (PRESENT(rcpair)) scptb_parameter%rcpair = rcpair

  END SUBROUTINE set_scptb_parameter

! *****************************************************************************
  SUBROUTINE write_scptb_parameter(scptb_parameter,section,error)

    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_scptb_parameter', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atomname, key, &
                                                parameterization
    INTEGER                                  :: i, l, lmaxorb, lmaxscp, &
                                                natorb, output_unit
    INTEGER, DIMENSION(0:3)                  :: norb
    INTEGER, DIMENSION(10, 0:3)              :: nqm
    LOGICAL                                  :: defined
    REAL(KIND=dp)                            :: ag, energy, rcpair, zeff
    REAL(KIND=dp), DIMENSION(10, 0:3)        :: hcore, zeta
    REAL(KIND=dp), DIMENSION(1:3)            :: pol
    REAL(KIND=dp), DIMENSION(3)              :: crep
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (ASSOCIATED(scptb_parameter).AND.&
        BTEST(cp_print_key_should_output(logger%iter_info,section,&
         "PRINT%KINDS/POTENTIAL",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,section,"PRINT%KINDS",&
            extension=".Log",error=error)

       IF (output_unit >0) THEN
          CALL get_scptb_parameter(scptb_parameter,key=key,atomname=atomname,&
               parameterization=parameterization,defined=defined,error=error)

          WRITE (UNIT=output_unit,FMT="(/,T10,A,T67,A14)") " SCPTB  parameters: ",TRIM(atomname)
          WRITE (UNIT=output_unit,FMT="(T67,A14)") TRIM(key)
          WRITE (UNIT=output_unit,FMT="(T67,A14)") TRIM(parameterization)
          IF ( defined ) THEN
            CALL get_scptb_parameter(scptb_parameter,&
                 zeff=zeff,natorb=natorb,lmaxorb=lmaxorb,lmaxscp=lmaxscp,&
                 norb=norb,nqm=nqm,zeta=zeta,hcore=hcore,&
                 energy=energy,crep=crep,pol=pol,ag=ag,rcpair=rcpair,&
                 error=error)
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)") "Effective core charge:",zeff
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)") "Total number of orbitals:",natorb
            WRITE (UNIT=output_unit,FMT="(T16,A,T30,A,T71,A)") "l-QM n-QM  i","Exponent","H core"
            DO l=0,lmaxorb
               DO i=1,norb(l)
                  WRITE (UNIT=output_unit,FMT="(T16,3I4,T30,F12.6,T71,F10.3)") &
                        l,i,nqm(i,l),zeta(i,l),hcore(i,l)
               END DO
            END DO
            WRITE (UNIT=output_unit,FMT="(T16,A,T61,F20.10)") "Energy of free atom [au]:",energy
            WRITE (UNIT=output_unit,FMT="(T16,A,T41,4F10.4)") "Core repulsion potential ",(crep(i),i=1,4)
            WRITE (UNIT=output_unit,FMT="(T16,A,T41,4F10.4)") "SCP Polarization ",(pol(l),l=1,lmaxscp)
            WRITE (UNIT=output_unit,FMT="(T16,A,T61,F20.5)") "Gaussian exponents for SCP ",ag
            WRITE (UNIT=output_unit,FMT="(T16,A,T61,F20.5)") "Pair potential cutoff (core) [bohr]:",rcpair
          ELSE
            WRITE (UNIT=output_unit,FMT="(T55,A)")&
                 "Parameters are not defined"
         END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,section,&
            "PRINT%KINDS",error=error)
    END IF

  END SUBROUTINE write_scptb_parameter

! *****************************************************************************
  SUBROUTINE scp_vector_create(scp_vec, nkind, natoms, nbasis, error)
    TYPE(scp_vector_type), POINTER           :: scp_vec
    INTEGER, INTENT(IN)                      :: nkind
    INTEGER, DIMENSION(:), INTENT(IN)        :: natoms, nbasis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, tl
    LOGICAL                                  :: failure

    CALL scp_vector_release(scp_vec, error)

    ALLOCATE(scp_vec,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(scp_vec%vector(nkind),stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    tl=0
    DO ikind=1,nkind
       scp_vec%vector(ikind)%natom = natoms(ikind)
       scp_vec%vector(ikind)%nbasis = nbasis(ikind)
       tl=tl+natoms(ikind)*nbasis(ikind)
       ALLOCATE(scp_vec%vector(ikind)%vmat(nbasis(ikind),natoms(ikind)),stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       scp_vec%vector(ikind)%vmat=0._dp
    END DO
    scp_vec%total_length = tl

  END SUBROUTINE scp_vector_create

  SUBROUTINE scp_vector_release(scp_vec, error)
    TYPE(scp_vector_type), POINTER           :: scp_vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat
    LOGICAL                                  :: failure

    IF(ASSOCIATED(scp_vec)) THEN
      IF(ASSOCIATED(scp_vec%vector)) THEN
        DO ikind=1,SIZE(scp_vec%vector)
           DEALLOCATE(scp_vec%vector(ikind)%vmat,stat=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO
        DEALLOCATE(scp_vec%vector,stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(scp_vec%vector)
      END IF
      DEALLOCATE(scp_vec,stat=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE scp_vector_release

  SUBROUTINE scp_vector_set(scp_vec, value, error)
    TYPE(scp_vector_type)                    :: scp_vec
    REAL(KIND=dp), INTENT(IN)                :: value
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(scp_vec%vector)
      scp_vec%vector(i)%vmat(:,:) = value
    END DO

  END SUBROUTINE scp_vector_set

  SUBROUTINE scp_vector_add(alpha, vecx, vecy, error)
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(scp_vector_type)                    :: vecx, vecy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(vecx%vector)
      vecy%vector(i)%vmat(:,:) = vecy%vector(i)%vmat(:,:) + alpha*vecx%vector(i)%vmat(:,:)
    END DO

  END SUBROUTINE scp_vector_add

  SUBROUTINE scp_vector_copy(vec_in, vec_out, error)
    TYPE(scp_vector_type)                    :: vec_in, vec_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(vec_in%vector)
      vec_out%vector(i)%vmat(:,:) = vec_in%vector(i)%vmat(:,:)
    END DO

  END SUBROUTINE scp_vector_copy

  SUBROUTINE scp_vector_scale(scale, vec, error)
    REAL(KIND=dp), INTENT(IN)                :: scale
    TYPE(scp_vector_type)                    :: vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(vec%vector)
      vec%vector(i)%vmat(:,:) = scale*vec%vector(i)%vmat(:,:)
    END DO

  END SUBROUTINE scp_vector_scale

  SUBROUTINE scp_vector_mult(vecx, vecy, error)
    TYPE(scp_vector_type)                    :: vecx, vecy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_mult', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(vecx%vector)
      vecy%vector(i)%vmat(:,:) = vecy%vector(i)%vmat(:,:) * vecx%vector(i)%vmat(:,:)
    END DO

  END SUBROUTINE scp_vector_mult

  SUBROUTINE scp_vector_dot(RESULT, vec1, vec2, error)
    REAL(KIND=dp), INTENT(OUT)               :: RESULT
    TYPE(scp_vector_type)                    :: vec1, vec2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_dot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

     RESULT = 0._dp
     DO i=1,SIZE(vec1%vector)
       RESULT = RESULT + SUM(vec1%vector(i)%vmat(:,:)*vec2%vector(i)%vmat(:,:))
     END DO

  END SUBROUTINE scp_vector_dot

  SUBROUTINE scp_vector_norm(RESULT, vec, error)
    REAL(KIND=dp), INTENT(OUT)               :: RESULT
    TYPE(scp_vector_type)                    :: vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

     RESULT = 0._dp
     DO i=1,SIZE(vec%vector)
       RESULT = RESULT + SUM(vec%vector(i)%vmat(:,:)*vec%vector(i)%vmat(:,:))
     END DO
     RESULT = SQRT ( RESULT/REAL(vec%total_length,KIND=dp) )

  END SUBROUTINE scp_vector_norm

  SUBROUTINE scp_vector_sync(vec, para_env, error)
    TYPE(scp_vector_type)                    :: vec
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_sync', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

     DO i=1,SIZE(vec%vector)
       CALL mp_sum(vec%vector(i)%vmat(:,:),para_env%group)
     END DO

  END SUBROUTINE scp_vector_sync

  SUBROUTINE scp_vector_print(vec, error)
    TYPE(scp_vector_type)                    :: vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_vector_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ia, ib, j, na, nb

    WRITE(*,*) "Index   Type   Atom    Basis       Value"
    j = 0
    DO i=1,SIZE(vec%vector)
       na = vec%vector(i)%natom
       nb = vec%vector(i)%nbasis
       DO ia = 1,na
          DO ib = 1,nb
             j = j+1
             WRITE(*,"(i5,i7,i7,i7,F20.10)") j,i,ia,ib,vec%vector(i)%vmat(ib,ia)
          END DO
       END DO
    END DO

  END SUBROUTINE scp_vector_print

! *****************************************************************************

END MODULE scptb_types

