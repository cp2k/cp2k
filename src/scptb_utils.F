!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for SCPTB
!> \author JGH (12.2011), RK (07.2012)
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
  USE external_potential_types,        ONLY: set_potential
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE qs_kind_types,                   ONLY: qs_kind_type,&
                                             set_qs_kind
  USE scptb_parameters,                ONLY: scptb_default_parameter
  USE scptb_types,                     ONLY: allocate_scptb_parameter,&
                                             scptb_parameter_type,&
                                             set_scptb_parameter
  USE string_utilities,                ONLY: uppercase
  USE xml_parser,                      ONLY: XML_PARSE,&
                                             xml_close,&
                                             xml_get,&
                                             xml_ok,&
                                             xml_open
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_utils'

  PUBLIC :: scptb_parameter_init

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param scptb_control ...
!> \param print_section ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_parameter_init(atomic_kind_set,qs_kind_set,scptb_control,print_section,para_env,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scptb_control_type), INTENT(INOUT)  :: scptb_control
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
    REAL(KIND=dp)                            :: a, c, d, eff, r, rmax, rmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet
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

    ! read global parameters from xml file
    IF(para_env%ionode) CALL scptb_global_parameter_from_file(scptb_control,error)
    CALL scptb_global_parameter_broadcast(scptb_control,para_env,error)

    nkind  = SIZE(atomic_kind_set)
    DO ikind = 1, nkind
      CALL get_atomic_kind(atomic_kind_set(ikind), name=iname)
      CALL uppercase(iname)
      NULLIFY(scptb_parameter)
      CALL allocate_scptb_parameter(scptb_parameter,error)
      CALL set_scptb_parameter(scptb_parameter,key=iname,defined=.FALSE., error=error)
    
      ! only for ionode
      IF(para_env%ionode) CALL scptb_parameter_from_file(scptb_parameter,scptb_control,error)
      CALL scptb_parameter_broadcast(scptb_parameter,para_env,error)
      
      IF (scptb_parameter%defined .EQV. .FALSE.) THEN
        output_unit = cp_print_key_unit_nr(logger,print_section,"KINDS",extension=".Log",error=error)
        IF ( output_unit > 0 ) THEN
           WRITE(output_unit, "(A,I3,2X,A,A)") " SCPTB| Kind ", ikind, TRIM(iname), &
            ": No parameter xml given or element undefined. Using defaults."
        END IF
        !defaults, but undefined for most elements atm
        CALL scptb_default_parameter(scptb_parameter,error) 
        CALL cp_print_key_finished_output(output_unit,logger,print_section,"KINDS",error=error)
      END IF

      !initialization

      IF (scptb_parameter%defined) THEN
         CALL set_potential(potential=qs_kind_set(ikind)%all_potential,zeff=scptb_parameter%zeff,&
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
         !Aggregate all the electron shell parameters into arrays
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
            CALL create_gto_from_sto_basis(sto_basis, qs_kind_set(ikind)%orb_basis_set,&
                                           scptb_control%sto_ng,error=error)
         END IF

         DEALLOCATE (nq,lq,zet,symbol,STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

         CALL deallocate_sto_basis_set (sto_basis,error)
         !throw away the STO basis set, since from now on we use the GTO

         a = scptb_parameter%crep(1)
         c = scptb_parameter%crep(2)
         d = scptb_parameter%crep(3)

         rmin = 1._dp
         rmax = 100._dp
         DO
            r = 0.5_dp*(rmax+rmin)
            eff = (c+d*r**2)*EXP(-(a*r))
            IF (eff < scptb_control%epspair) THEN
               rmax = r
            ELSE
               rmin = r
            END IF
            IF ((rmax-rmin) < 1.e-2) EXIT
         END DO
         scptb_parameter%rcpair = 0.5_dp*rmax

         !lmaxscp: make sure it is also consistent with pure SCC-TB
         IF(scptb_parameter%lmaxscp < 0) THEN
            IF(scptb_parameter%ag > 0.0_dp) THEN
               scptb_parameter%lmaxscp = 0
            END IF
         END IF

      END IF

      CALL set_qs_kind(qs_kind_set(ikind), scptb_parameter=scptb_parameter)

    END DO

  END SUBROUTINE scptb_parameter_init

! *****************************************************************************
!> \brief ...
!> \param scptb_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_global_parameter_from_file(scptb_control,error)

    TYPE(scptb_control_type), INTENT(INOUT)  :: scptb_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'scptb_global_parameter_from_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: attribs(2,10), author, &
                                                DATA(10), fname, found_key, &
                                                tag, vers
    INTEGER                                  :: no_attribs, no_data, &
                                                tag_count, unrecognized_tags
    LOGICAL :: endtag, has_A, has_ag, has_B, has_C, has_D, has_energy, &
      has_eta, has_hcore, has_nqm, has_occ, has_z, has_zeff, tag_over
    TYPE(XML_PARSE)                          :: xml_info

    DATA has_z, has_zeff, has_ag, has_energy, has_hcore, has_occ,   &
                    has_nqm, has_eta, has_A, has_B, has_C, has_D/12*.FALSE./

    fname = scptb_control%parameter_file
    IF ( fname /= "" ) THEN
        CALL xml_open(xml_info, fname, .TRUE.)
       
        tag_count = 0
        tag=""
        DATA=""
        attribs=""
        unrecognized_tags = 0
        tag_over = .FALSE.
        found_key = ""
        author = "unknown author"
        vers = " v 0.0"

        DO WHILE ( xml_ok(xml_info) )   !Until the file is read entirely
            CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
            CALL uppercase(tag)
            SELECT CASE (TRIM(tag)) !Decide on the first-level tag
            CASE ("")
                !just garbage from the parser
            CASE ("NM12_PARAMETERS")
                !just the header, do nothing
            CASE ("META_INFO")
                tag_over = .FALSE.
                DO WHILE (tag_over .EQV. .FALSE.)
                    CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                    CALL uppercase(tag)
                    SELECT CASE (tag)
                    CASE ("META_INFO")
                        IF (endtag .EQV. .TRUE.) tag_over = .TRUE.
                    CASE ("VERSION")
                        IF (endtag .EQV. .FALSE.) vers = DATA(1)
                    CASE ("AUTHOR")
                        IF (endtag .EQV. .FALSE.) author = DATA(1)
                    CASE ("REF")
                        !IF (endtag .eqv. .FALSE.)  !do nothing
                    CASE ("COMMENT")
                        !IF (endtag .eqv. .FALSE.)  !do nothing
                    CASE DEFAULT
                        unrecognized_tags = unrecognized_tags + 1
                    END SELECT
                END DO

            CASE ("GENERAL_PARAMETERS") !general params
                tag_over = .FALSE.
                DO WHILE (tag_over .EQV. .FALSE.)
                    CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                    CALL uppercase(tag)
                    SELECT CASE (tag)
                    CASE ("GENERAL_PARAMETERS")
                        IF (endtag .EQV. .TRUE.) tag_over = .TRUE.
                    CASE ("D3_S6")
                        IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_control%sd3(1) 
                    CASE ("D3_SR6")
                        IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_control%sd3(2) 
                    CASE ("D3_S8")
                        IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_control%sd3(3) 
                    CASE DEFAULT
                        unrecognized_tags = unrecognized_tags + 1
                    END SELECT
                END DO

            CASE ("ELEMENT") !now find the element we are looking for
                tag_count=tag_count+1
                tag_over = .FALSE.
                
           CASE DEFAULT
                unrecognized_tags = unrecognized_tags + 1
           END SELECT

       END DO

       CALL xml_close(xml_info) !That's it, close the file.
        
    END IF

  END SUBROUTINE scptb_global_parameter_from_file

! *****************************************************************************
!> \brief ...
!> \param scptb_parameter ...
!> \param scptb_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_parameter_from_file(scptb_parameter,scptb_control,error)

    TYPE(scptb_parameter_type), &
      INTENT(INOUT)                          :: scptb_parameter
    TYPE(scptb_control_type), INTENT(INOUT)  :: scptb_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scptb_parameter_from_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: attribs(2,10), author, &
                                                DATA(10), fname, found_key, &
                                                search_key, tag, vers
    INTEGER                                  :: j, l, n, no_attribs, no_data, &
                                                tag_count, unrecognized_tags
    LOGICAL :: endtag, has_A, has_ag, has_B, has_C, has_D, has_energy, &
      has_eta, has_hcore, has_nqm, has_occ, has_z, has_zeff, tag_over, &
      tag_over2
    TYPE(XML_PARSE)                          :: xml_info

    DATA has_z, has_zeff, has_ag, has_energy, has_hcore, has_occ,   &
                    has_nqm, has_eta, has_A, has_B, has_C, has_D/12*.FALSE./

    fname = scptb_control%parameter_file
    IF ( fname /= "" ) THEN
        CALL xml_open(xml_info, fname, .TRUE.)
       
        tag_count = 0
        tag=""
        DATA=""
        attribs=""
        unrecognized_tags = 0
        tag_over = .FALSE.
        found_key = ""
        author = "unknown author"
        vers = " v 0.0"

        search_key = TRIM(scptb_parameter%key)
        CALL uppercase(search_key)

        DO WHILE ( xml_ok(xml_info) )   !Until the file is read entirely
            CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
            CALL uppercase(tag)
            SELECT CASE (TRIM(tag)) !Decide on the first-level tag
            CASE ("")
                !just garbage from the parser
            CASE ("NM12_PARAMETERS")
                !just the header, do nothing
            CASE ("META_INFO")
                tag_over = .FALSE.
                DO WHILE (tag_over .EQV. .FALSE.)
                    CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                    CALL uppercase(tag)
                    SELECT CASE (tag)
                    CASE ("META_INFO")
                        IF (endtag .EQV. .TRUE.) tag_over = .TRUE.
                    CASE ("VERSION")
                        IF (endtag .EQV. .FALSE.) vers = DATA(1)
                    CASE ("AUTHOR")
                        IF (endtag .EQV. .FALSE.) author = DATA(1)
                    CASE ("REF")
                        !IF (endtag .eqv. .FALSE.)  !do nothing
                    CASE ("COMMENT")
                        !IF (endtag .eqv. .FALSE.)  !do nothing
                    CASE DEFAULT
                        unrecognized_tags = unrecognized_tags + 1
                    END SELECT
                END DO

                WRITE(scptb_parameter%parameterization,'(A,A,A)') TRIM(author), ", ", TRIM(vers)    

            CASE ("GENERAL_PARAMETERS") !general params
                tag_over = .FALSE.
                DO WHILE (tag_over .EQV. .FALSE.)
                    CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                    CALL uppercase(tag)
                    SELECT CASE (tag)
                    CASE ("GENERAL_PARAMETERS")
                        IF (endtag .EQV. .TRUE.) tag_over = .TRUE.
                    CASE ("D3_S6")
                        ! do nothing in this sweep of data
                    CASE ("D3_SR6")
                        ! do nothing in this sweep of data
                    CASE ("D3_S8")
                        ! do nothing in this sweep of data
                    CASE DEFAULT
                        unrecognized_tags = unrecognized_tags + 1
                    END SELECT
                END DO

            CASE ("ELEMENT") !now find the element we are looking for
                tag_count=tag_count+1
                tag_over = .FALSE.
                found_key = TRIM(attribs(2,1))
                CALL uppercase(found_key)
                
                IF (found_key == search_key) THEN
                    ! read all the values

                    scptb_parameter%atomname = search_key


                    DO WHILE (tag_over .EQV. .FALSE.)
                        CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                        CALL uppercase(tag)
                        SELECT CASE (tag)
                        CASE ("ELEMENT")
                            IF (endtag .EQV. .TRUE.) tag_over = .TRUE.
                        CASE ("Z")
                            IF (endtag .EQV. .FALSE.) READ(DATA(1),'(I4)') scptb_parameter%z ; has_z = .TRUE.
                        CASE ("ZEFF")
                            IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_parameter%zeff ; has_zeff = .TRUE.
                        CASE ("CHARGE_WIDTH")
                            IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_parameter%ag ; has_ag = .TRUE.
                        CASE ("ENERGY")
                            IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_parameter%energy ; has_energy = .TRUE.

                        CASE ("ELECTRON_SHELL")
                            !first deal with the attributes
                            DO j=1,no_attribs
                                IF (attribs(1,j) == "n") READ(attribs(2,j),*) n
                                IF (attribs(1,j) == "l") READ(attribs(2,j),*) l
                            END DO

                            !now deal with the sub tags
                            tag_over2 = .FALSE.
                            DO WHILE (tag_over2 .EQV. .FALSE.)
                                CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                                CALL uppercase(tag)
                                SELECT CASE (tag)
                                    CASE ("ELECTRON_SHELL")
                                        IF (endtag .EQV. .TRUE.) tag_over2 = .TRUE.
                                    CASE ("ZETA")
                                        !deal with the values
                                        IF (endtag .EQV. .FALSE.) THEN
                                            READ(DATA(1),*) scptb_parameter%zeta(n,l)
                                            scptb_parameter%natorb = scptb_parameter%natorb + 1
                                            scptb_parameter%norb(l) = scptb_parameter%norb(l) + 1
                                            IF (l > scptb_parameter%lmaxorb) scptb_parameter%lmaxorb = l
                                        END IF
                                    CASE ("H_CORE")    
                                        IF (endtag .EQV. .FALSE.) THEN
                                            READ(DATA(1),*) scptb_parameter%hcore(n,l) 
                                            has_hcore = .TRUE.
                                        END IF
                                    CASE ("OCC")
                                        IF (endtag .EQV. .FALSE.) &
                                           READ(DATA(1),*) scptb_parameter%occupation(n,l) ; has_occ = .TRUE.
                                    CASE ("N")
                                        IF (endtag .EQV. .FALSE.) &
                                           READ(DATA(1),*) scptb_parameter%nqm(n,l) ; has_nqm = .TRUE.
                                    CASE DEFAULT
                                        IF (endtag .EQV. .FALSE.) unrecognized_tags = unrecognized_tags + 1
                                END SELECT 
                            END DO

                        CASE ("SCP")
                            !deal with the sub tags
                            tag_over2 = .FALSE.
                            DO WHILE (tag_over2 .EQV. .FALSE.)
                                CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                                CALL uppercase(tag)
                                SELECT CASE (tag)
                                    CASE ("SCP")
                                        IF (endtag .EQV. .TRUE.) tag_over2 = .TRUE.
                                    CASE ("ETA")
                                        !first deal with the attributes
                                        DO j=1,no_attribs
                                            IF (attribs(1,j) == "l") READ(attribs(2,j),*) l ; has_eta = .TRUE.
                                        END DO
                                        IF (endtag .EQV. .FALSE.) THEN
                                            IF(l>0) READ(DATA(1),*) scptb_parameter%pol(l)
                                            IF (l > scptb_parameter%lmaxscp) scptb_parameter%lmaxscp = l
                                        END IF
                                    CASE DEFAULT
                                        IF (endtag .EQV. .FALSE.) unrecognized_tags = unrecognized_tags + 1
                                END SELECT 
                            END DO

                        CASE ("REPULSIVE")
                            !deal with the sub tags
                            tag_over2 = .FALSE.
                            DO WHILE (tag_over2 .EQV. .FALSE.)
                                CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                                CALL uppercase(tag)
                                SELECT CASE (tag)
                                    CASE ("REPULSIVE")
                                        IF (endtag .EQV. .TRUE.) tag_over2 = .TRUE.
                                    CASE ("A")
                                        IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_parameter%crep(1) ; has_A = .TRUE.
                                    CASE ("C")  
                                        IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_parameter%crep(2) ; has_C = .TRUE.
                                    CASE ("D")
                                        IF (endtag .EQV. .FALSE.) READ(DATA(1),*) scptb_parameter%crep(3) ; has_D = .TRUE.
                                    CASE DEFAULT
                                        IF (endtag .EQV. .FALSE.) unrecognized_tags = unrecognized_tags + 1
                                END SELECT 
                            END DO
                        CASE DEFAULT
                            IF (endtag .EQV. .FALSE.) unrecognized_tags = unrecognized_tags + 1
                        END SELECT
                    END DO
                ELSE
                    !when the current element is not the one we are looking for,
                    !just keep running through the file until the next element.
                    DO WHILE (tag_over .EQV. .FALSE.)
                        CALL xml_get(xml_info, tag, endtag, attribs, no_attribs, DATA, no_data)
                        CALL uppercase(tag)
                        IF (tag=="ELEMENT" .AND. endtag .EQV. .TRUE.) tag_over = .TRUE.
                    END DO
                END IF

           CASE DEFAULT
                unrecognized_tags = unrecognized_tags + 1
           END SELECT

       END DO

       CALL xml_close(xml_info) !That's it, close the file.
        
       !An element is only considered properly defined if all non-optional parameters were
       !specified. That doesn't safeguard against incomplete or meaningless input.
       IF ( has_z .AND. has_zeff .AND. has_ag .AND. has_hcore .AND.   &
            has_nqm .AND. has_A .AND. has_C .AND. has_D) THEN
           scptb_parameter%defined = .TRUE. 
       ELSE
           scptb_parameter%defined = .FALSE. 
       END IF

    END IF

  END SUBROUTINE scptb_parameter_from_file

! *****************************************************************************
!> \brief ...
!> \param scptb_control ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_global_parameter_broadcast(scptb_control,para_env,error)

    TYPE(scptb_control_type), INTENT(INOUT)  :: scptb_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'scptb_global_parameter_broadcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, source
    LOGICAL                                  :: failure

    failure = .FALSE.

    source = para_env%source
    group = para_env%group
    
    IF(para_env%ionode) THEN
       CPPostcondition(source==para_env%mepos,cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(source/=para_env%mepos,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_bcast(scptb_control%sd3,source,group)

  END SUBROUTINE scptb_global_parameter_broadcast

! *****************************************************************************
!> \brief ...
!> \param scptb_parameter ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_parameter_broadcast(scptb_parameter,para_env,error)

    TYPE(scptb_parameter_type), &
      INTENT(INOUT)                          :: scptb_parameter
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scptb_parameter_broadcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, source
    LOGICAL                                  :: failure

    failure = .FALSE.

    source = para_env%source
    group = para_env%group
    
    IF(para_env%ionode) THEN
       CPPostcondition(source==para_env%mepos,cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(source/=para_env%mepos,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_bcast(scptb_parameter%key,source,group)
    CALL mp_bcast(scptb_parameter%parameterization,source,group)
    CALL mp_bcast(scptb_parameter%atomname,source,group)
    CALL mp_bcast(scptb_parameter%defined,source,group)
    CALL mp_bcast(scptb_parameter%z,source,group)
    CALL mp_bcast(scptb_parameter%zeff,source,group)
    CALL mp_bcast(scptb_parameter%natorb,source,group)
    CALL mp_bcast(scptb_parameter%lmaxorb,source,group)
    CALL mp_bcast(scptb_parameter%lmaxscp,source,group)
    CALL mp_bcast(scptb_parameter%norb,source,group)
    CALL mp_bcast(scptb_parameter%nqm,source,group)
    CALL mp_bcast(scptb_parameter%zeta,source,group)
    CALL mp_bcast(scptb_parameter%hcore,source,group)
    CALL mp_bcast(scptb_parameter%occupation,source,group)
    CALL mp_bcast(scptb_parameter%energy,source,group)
    CALL mp_bcast(scptb_parameter%crep,source,group)
    CALL mp_bcast(scptb_parameter%pol,source,group)
    CALL mp_bcast(scptb_parameter%ag,source,group)
    CALL mp_bcast(scptb_parameter%rcpair,source,group)

  END SUBROUTINE scptb_parameter_broadcast

END MODULE scptb_utils

