!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Subroutine input_torsions changed (DG) 05-Dec-2000
!>      Output formats changed (DG) 05-Dec-2000
!>      JGH (26-01-2002) : force field parameters stored in tables, not in
!>        matrices. Input changed to have parameters labeled by the position
!>        and not atom pairs (triples etc)
!>      Teo (11.2005) : Moved all information on force field  pair_potential to
!>                      a much lighter memory structure
!>      Teo 09.2006   : Splitted all routines force_field I/O in a separate file
!>      10.2008 Teodoro Laino [tlaino] :  splitted from force_fields_input in order
!>              to collect in a unique module only the functions to read external FF
!> \author CJM
! *****************************************************************************
MODULE force_fields_ext
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object,&
                                             parser_search_string,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE force_field_kind_types,          ONLY: do_ff_g96
  USE force_field_types,               ONLY: amber_info_type,&
                                             charmm_info_type,&
                                             force_field_type,&
                                             gromos_info_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_amber,                  ONLY: rdparm_amber_8
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_ext'

  PRIVATE
  PUBLIC :: read_force_field_charmm,&
            read_force_field_amber,&
            read_force_field_gromos

CONTAINS

! *****************************************************************************
!> \brief Reads the GROMOS force_field
!> \param ff_type ...
!> \param para_env ...
!> \param mm_section ...
!> \param error ...
!> \author ikuo
! *****************************************************************************
  SUBROUTINE read_force_field_gromos ( ff_type , para_env,  mm_section, error )

    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_gromos', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: ekt = 2.5_dp

    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length), &
      DIMENSION(21)                          :: avail_section
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: namearray(:)
    INTEGER                                  :: handle, iatom, icon, itemp, &
                                                itype, iw, jatom, ncon, &
                                                ntype, offset, stat
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: cosphi0, cost2, csq, sdet
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(gromos_info_type), POINTER          :: gro_info

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger,parser)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    avail_section( 1) = "TITLE"
    avail_section( 2) = "TOPPHYSCON"
    avail_section( 3) = "TOPVERSION"
    avail_section( 4) = "ATOMTYPENAME"
    avail_section( 5) = "RESNAME"
    avail_section( 6) = "SOLUTEATOM"
    avail_section( 7) = "BONDTYPE"
    avail_section( 8) = "BONDH"
    avail_section( 9) = "BOND"
    avail_section(10) = "BONDANGLETYPE"
    avail_section(11) = "BONDANGLEH"
    avail_section(12) = "BONDANGLE"
    avail_section(13) = "IMPDIHEDRALTYPE"
    avail_section(14) = "IMPDIHEDRALH"
    avail_section(15) = "IMPDIHEDRAL"
    avail_section(16) = "DIHEDRALTYPE"
    avail_section(17) = "DIHEDRALH"
    avail_section(18) = "DIHEDRAL"
    avail_section(19) = "LJPARAMETERS"
    avail_section(20) = "SOLVENTATOM"
    avail_section(21) = "SOLVENTCONSTR"

    gro_info => ff_type%gro_info
    gro_info%ff_gromos_type = ff_type%ff_type
    NULLIFY(namearray)
    ! ATOMTYPENAME SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the ATOMTYPENAME section'
    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(4))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(namearray,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,namearray(itype),lower_to_upper=.TRUE.,error=error)
          IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(namearray(itype))
       END DO
    END IF
    CALL parser_release(parser,error=error)

    ! SOLVENTCONSTR SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the SOLVENTATOM section'
    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)

    label = TRIM(avail_section(21))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ncon,error=error)
       CALL reallocate(gro_info%solvent_k,1,ncon)
       CALL reallocate(gro_info%solvent_r0,1,ncon)
       DO icon=1,ncon
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          CALL parser_get_object(parser,gro_info%solvent_r0(icon),error=error)
          gro_info%solvent_k(icon)=0.0_dp
       END DO
    END IF
    CALL parser_release(parser,error=error)

    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
    ! BONDTYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the BONDTYPE section'
    label = TRIM(avail_section(7))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%bond_k,1,ntype)
       CALL reallocate(gro_info%bond_r0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%bond_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%bond_r0(itype),error=error)
          IF(ff_type%ff_type==do_ff_g96) THEN
             gro_info%bond_k(itype) = cp_unit_to_cp2k(gro_info%bond_k(itype),"kjmol*nm^-4",error=error)
          ELSE ! Assume its G87
             gro_info%bond_k(itype) = (2.0_dp) * gro_info%bond_k(itype) * gro_info%bond_r0(itype)**2
             gro_info%bond_k(itype) = cp_unit_to_cp2k(gro_info%bond_k(itype),"kjmol*nm^-2",error=error)
          END IF
          gro_info%bond_r0(itype)= cp_unit_to_cp2k(gro_info%bond_r0(itype),"nm",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT BONDTYPE INFO HERE!!!!"
       END DO
    END IF

    ! BONDANGLETYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the BONDANGLETYPE section'
    label = TRIM(avail_section(10))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%bend_k,1,ntype)
       CALL reallocate(gro_info%bend_theta0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%bend_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%bend_theta0(itype),error=error)
          gro_info%bend_theta0(itype) = cp_unit_to_cp2k(gro_info%bend_theta0(itype),"deg",error=error)
          IF(ff_type%ff_type==do_ff_g96) THEN
             gro_info%bend_theta0(itype) = COS(gro_info%bend_theta0(itype))
          ELSE ! Assume its G87
             cost2 = COS(gro_info%bend_theta0(itype))*COS(gro_info%bend_theta0(itype))
             sdet = cost2*cost2 - (2.0_dp*cost2-1.0_dp)*(1.0_dp-ekt/gro_info%bend_k(itype))
             csq = (cost2-SQRT(sdet))/(2.0_dp*cost2-1.0_dp)
             gro_info%bend_k(itype) = ekt/ACOS(csq)**2
          END IF
          gro_info%bend_k(itype) = cp_unit_to_cp2k(gro_info%bend_k(itype),"kjmol",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT BONDANGLETYPE INFO HERE!!!!"
       END DO
    END IF

    ! IMPDIHEDRALTYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the IMPDIHEDRALTYPE section'
    label = TRIM(avail_section(13))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%impr_k,1,ntype)
       CALL reallocate(gro_info%impr_phi0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%impr_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%impr_phi0(itype),error=error)
          gro_info%impr_phi0(itype) = cp_unit_to_cp2k(gro_info%impr_phi0(itype),"deg",error=error)
          gro_info%impr_k(itype) = cp_unit_to_cp2k(gro_info%impr_k(itype),"kjmol*deg^-2",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT IMPDIHEDRALTYPE INFO HERE!!!!"
       END DO
    END IF

    ! DIHEDRALTYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the DIHEDRALTYPE section'
    label = TRIM(avail_section(16))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%torsion_k,1,ntype)
       CALL reallocate(gro_info%torsion_m,1,ntype)
       CALL reallocate(gro_info%torsion_phi0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%torsion_k(itype),error=error)
          CALL parser_get_object(parser,cosphi0,error=error)
          CALL parser_get_object(parser,gro_info%torsion_m(itype),error=error)
          gro_info%torsion_phi0(itype) = ACOS(cosphi0)
          gro_info%torsion_k(itype) = cp_unit_to_cp2k(gro_info%torsion_k(itype),"kjmol",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT DIHEDRALTYPE INFO HERE!!!!"
       END DO
    END IF

    CALL parser_release(parser,error=error)

    ! LJPARAMETERS SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the LJPARAMETERS section'
    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(19))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(gro_info%nonbond_a)) offset = SIZE(gro_info%nonbond_a)
       ntype=SIZE(namearray)
       CALL reallocate(gro_info%nonbond_a,1,ntype)
       CALL reallocate(gro_info%nonbond_a_14,1,ntype)
       CALL reallocate(gro_info%nonbond_c6,1,ntype,1,ntype)
       CALL reallocate(gro_info%nonbond_c12,1,ntype,1,ntype)
       CALL reallocate(gro_info%nonbond_c6_14,1,ntype,1,ntype)
       CALL reallocate(gro_info%nonbond_c12_14,1,ntype,1,ntype)

       gro_info%nonbond_c12 = 0._dp
       gro_info%nonbond_c6  = 0._dp
       gro_info%nonbond_c12_14 = 0._dp
       gro_info%nonbond_c6_14  = 0._dp

       DO itype=1,ntype*(ntype+1)/2
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,iatom,error=error)
          CALL parser_get_object(parser,jatom,error=error)
          IF(iatom==jatom) THEN
             gro_info%nonbond_a(iatom) = namearray(iatom)
             gro_info%nonbond_a_14(iatom) = namearray(iatom)
          END IF
          CALL parser_get_object(parser,gro_info%nonbond_c12(iatom,jatom),error=error)
          CALL parser_get_object(parser,gro_info%nonbond_c6(iatom,jatom),error=error)
          CALL parser_get_object(parser,gro_info%nonbond_c12_14(iatom,jatom),error=error)
          CALL parser_get_object(parser,gro_info%nonbond_c6_14(iatom,jatom),error=error)
          gro_info%nonbond_c6(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c6(iatom,jatom),"kjmol*nm^6",error=error)
          gro_info%nonbond_c12(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c12(iatom,jatom),"kjmol*nm^12",error=error)
          gro_info%nonbond_c6_14(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c6_14(iatom,jatom),"kjmol*nm^6",error=error)
          gro_info%nonbond_c12_14(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c12_14(iatom,jatom),"kjmol*nm^12",error=error)

          gro_info%nonbond_c6_14(jatom,iatom)= gro_info%nonbond_c6_14(iatom,jatom)
          gro_info%nonbond_c12_14(jatom,iatom)= gro_info%nonbond_c12_14(iatom,jatom)
          gro_info%nonbond_c6(jatom,iatom)= gro_info%nonbond_c6(iatom,jatom)
          gro_info%nonbond_c12(jatom,iatom)= gro_info%nonbond_c12(iatom,jatom)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT LJPARAMETERS INFO HERE!!!!"
       END DO
    END IF
    CALL parser_release(parser,error=error)

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)

    DEALLOCATE(namearray,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE read_force_field_gromos

! *****************************************************************************
!> \brief Reads the charmm force_field
!> \param ff_type ...
!> \param para_env ...
!> \param mm_section ...
!> \param error ...
!> \author ikuo
! *****************************************************************************
  SUBROUTINE read_force_field_charmm ( ff_type, para_env, mm_section, error )

    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_charmm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: label, string, string2, &
                                                string3, string4
    CHARACTER(LEN=default_string_length), &
      DIMENSION(1)                           :: bond_section
    CHARACTER(LEN=default_string_length), &
      DIMENSION(19)                          :: avail_section
    CHARACTER(LEN=default_string_length), &
      DIMENSION(2)                           :: angl_section, impr_section, &
                                                nbon_section, thet_section
    INTEGER                                  :: dummy, handle, ilab, iw, &
                                                nbend, nbond, nimpr, &
                                                nnonbond, nonfo, ntorsion, nub
    LOGICAL                                  :: failure, found
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger,parser)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    avail_section(1) = "BOND"     ; bond_section(1) = avail_section(1)
    avail_section(11)= "BONDS"
    avail_section(2) = "ANGL"     ; angl_section(1) = avail_section(2)
    avail_section(3) = "THETA"    ; angl_section(2) = avail_section(3)
    avail_section(12)= "THETAS"
    avail_section(13)= "ANGLE"
    avail_section(14)= "ANGLES"
    avail_section(4) = "DIHEDRAL" ; thet_section(1) = avail_section(4)
    avail_section(15)= "DIHEDRALS"
    avail_section(5) = "PHI"      ; thet_section(2) = avail_section(5)
    avail_section(6) = "IMPROPER" ; impr_section(1) = avail_section(6)
    avail_section(7) = "IMPH"     ; impr_section(2) = avail_section(7)
    avail_section(16)= "IMPHI"
    avail_section(8) = "NONBONDED"; nbon_section(1) = avail_section(8)
    avail_section(9) = "NBOND"    ; nbon_section(2) = avail_section(9)
    avail_section(10)= "HBOND"
    avail_section(17)= "NBFIX"
    avail_section(18)= "CMAP"
    avail_section(19)= "END"

    chm_info => ff_type%chm_info
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Read in all the Bonds info from the param file here
    !      Vbond = Kb(b-b0)^2
    !      UNITS for Kb: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
    !-----------------------------------------------------------------------------
    nbond = 0
    DO ilab = 1, SIZE(bond_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(bond_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL uppercase ( string2 )
                nbond = nbond + 1
                CALL reallocate(chm_info%bond_a,1,nbond)
                CALL reallocate(chm_info%bond_b,1,nbond)
                CALL reallocate(chm_info%bond_k,1,nbond)
                CALL reallocate(chm_info%bond_r0,1,nbond)
                chm_info%bond_a(nbond) = string
                chm_info%bond_b(nbond) = string2
                CALL parser_get_object(parser,chm_info%bond_k(nbond),error=error)
                CALL parser_get_object(parser,chm_info%bond_r0(nbond),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM BOND ",nbond,&
                     TRIM(chm_info%bond_a(nbond))," ",&
                     TRIM(chm_info%bond_b(nbond))," ",&
                     chm_info%bond_k(nbond),&
                     chm_info%bond_r0(nbond)
                ! Do some units conversion into internal atomic units
                chm_info%bond_r0(nbond) = cp_unit_to_cp2k(chm_info%bond_r0(nbond),"angstrom",error=error)
                chm_info%bond_k(nbond) = cp_unit_to_cp2k(chm_info%bond_k(nbond),"kcalmol*angstrom^-2",error=error)
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in all the Bends and UB info from the param file here
    !      Vangle = Ktheta(theta-theta0)^2
    !      UNITS for Ktheta: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
    !      FACTOR of "2" rolled into Ktheta
    !      Vub = Kub(S-S0)^2
    !      UNITS for Kub: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
    !-----------------------------------------------------------------------------
    nbend = 0
    nub = 0
    DO ilab = 1, SIZE(angl_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(angl_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL parser_get_object(parser,string3,error=error)
                CALL uppercase ( string2 )
                CALL uppercase ( string3 )
                nbend = nbend + 1
                CALL reallocate(chm_info%bend_a,1,nbend)
                CALL reallocate(chm_info%bend_b,1,nbend)
                CALL reallocate(chm_info%bend_c,1,nbend)
                CALL reallocate(chm_info%bend_k,1,nbend)
                CALL reallocate(chm_info%bend_theta0,1,nbend)
                chm_info%bend_a(nbend) = string
                chm_info%bend_b(nbend) = string2
                chm_info%bend_c(nbend) = string3
                CALL parser_get_object(parser,chm_info%bend_k(nbend),error=error)
                CALL parser_get_object(parser,chm_info%bend_theta0(nbend),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM BEND ",nbend,&
                     TRIM(chm_info%bend_a(nbend))," ",&
                     TRIM(chm_info%bend_b(nbend))," ",&
                     TRIM(chm_info%bend_c(nbend))," ",&
                     chm_info%bend_k(nbend),&
                     chm_info%bend_theta0(nbend)
                ! Do some units conversion into internal atomic units
                chm_info%bend_theta0(nbend) = cp_unit_to_cp2k(chm_info%bend_theta0(nbend),"deg",error=error)
                chm_info%bend_k(nbend) = cp_unit_to_cp2k(chm_info%bend_k(nbend),"kcalmol*rad^-2",error=error)
                IF (parser_test_next_token(parser,error=error) == "FLT") THEN
                   nub = nub + 1
                   CALL reallocate(chm_info%ub_a,1,nub)
                   CALL reallocate(chm_info%ub_b,1,nub)
                   CALL reallocate(chm_info%ub_c,1,nub)
                   CALL reallocate(chm_info%ub_k,1,nub)
                   CALL reallocate(chm_info%ub_r0,1,nub)
                   chm_info%ub_a(nub) = string
                   chm_info%ub_b(nub) = string2
                   chm_info%ub_c(nub) = string3
                   CALL parser_get_object(parser,chm_info%ub_k(nub),error=error)
                   CALL parser_get_object(parser,chm_info%ub_r0(nub),error=error)
                   IF(iw>0) WRITE(iw,*) "    CHM UB ",nub,&
                        TRIM(chm_info%ub_a(nub))," ",&
                        TRIM(chm_info%ub_b(nub))," ",&
                        TRIM(chm_info%ub_c(nub))," ",&
                        chm_info%ub_k(nub),&
                        chm_info%ub_r0(nub)
                   ! Do some units conversion into internal atomic units
                   chm_info%ub_r0(nub) = cp_unit_to_cp2k(chm_info%ub_r0(nub),"angstrom",error=error)
                   chm_info%ub_k(nub) = cp_unit_to_cp2k(chm_info%ub_k(nub),"kcalmol*angstrom^-2",error=error)
                END IF
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Read in all the Dihedrals info from the param file here
    !      Vtorsion = Kphi(1+COS(n(phi)-delta))
    !      UNITS for Kphi: [(kcal/mol)] to [Eh]
    !-----------------------------------------------------------------------------
    ntorsion = 0
    DO ilab = 1, SIZE(thet_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(thet_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL parser_get_object(parser,string3,error=error)
                CALL parser_get_object(parser,string4,error=error)
                CALL uppercase ( string2 )
                CALL uppercase ( string3 )
                CALL uppercase ( string4 )
                ntorsion = ntorsion + 1
                CALL reallocate(chm_info%torsion_a,1,ntorsion)
                CALL reallocate(chm_info%torsion_b,1,ntorsion)
                CALL reallocate(chm_info%torsion_c,1,ntorsion)
                CALL reallocate(chm_info%torsion_d,1,ntorsion)
                CALL reallocate(chm_info%torsion_k,1,ntorsion)
                CALL reallocate(chm_info%torsion_m,1,ntorsion)
                CALL reallocate(chm_info%torsion_phi0,1,ntorsion)
                chm_info%torsion_a(ntorsion) = string
                chm_info%torsion_b(ntorsion) = string2
                chm_info%torsion_c(ntorsion) = string3
                chm_info%torsion_d(ntorsion) = string4
                CALL parser_get_object(parser,chm_info%torsion_k(ntorsion),error=error)
                CALL parser_get_object(parser,chm_info%torsion_m(ntorsion),error=error)
                CALL parser_get_object(parser,chm_info%torsion_phi0(ntorsion),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM TORSION ",ntorsion,&
                     TRIM(chm_info%torsion_a(ntorsion))," ",&
                     TRIM(chm_info%torsion_b(ntorsion))," ",&
                     TRIM(chm_info%torsion_c(ntorsion))," ",&
                     TRIM(chm_info%torsion_d(ntorsion))," ",&
                     chm_info%torsion_k(ntorsion),&
                     chm_info%torsion_m(ntorsion),&
                     chm_info%torsion_phi0(ntorsion)
                ! Do some units conversion into internal atomic units
                chm_info%torsion_phi0(ntorsion) = cp_unit_to_cp2k(chm_info%torsion_phi0(ntorsion),&
                     "deg",error=error)
                chm_info%torsion_k(ntorsion) = cp_unit_to_cp2k(chm_info%torsion_k(ntorsion),"kcalmol",&
                     error=error)
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Read in all the Improper info from the param file here
    !      Vimpr = Kpsi(psi-psi0)^2
    !      UNITS for Kpsi: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
    !-----------------------------------------------------------------------------
    nimpr = 0
    DO ilab = 1, SIZE(impr_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(impr_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL parser_get_object(parser,string3,error=error)
                CALL parser_get_object(parser,string4,error=error)
                CALL uppercase ( string2 )
                CALL uppercase ( string3 )
                CALL uppercase ( string4 )
                nimpr = nimpr + 1
                CALL reallocate(chm_info%impr_a,1,nimpr)
                CALL reallocate(chm_info%impr_b,1,nimpr)
                CALL reallocate(chm_info%impr_c,1,nimpr)
                CALL reallocate(chm_info%impr_d,1,nimpr)
                CALL reallocate(chm_info%impr_k,1,nimpr)
                CALL reallocate(chm_info%impr_phi0,1,nimpr)
                chm_info%impr_a(nimpr) = string
                chm_info%impr_b(nimpr) = string2
                chm_info%impr_c(nimpr) = string3
                chm_info%impr_d(nimpr) = string4
                CALL parser_get_object(parser,chm_info%impr_k(nimpr),error=error)
                CALL parser_get_object(parser,dummy,error=error)
                CALL parser_get_object(parser,chm_info%impr_phi0(nimpr),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM IMPROPERS ",nimpr,&
                     TRIM(chm_info%impr_a(nimpr))," ",&
                     TRIM(chm_info%impr_b(nimpr))," ",&
                     TRIM(chm_info%impr_c(nimpr))," ",&
                     TRIM(chm_info%impr_d(nimpr))," ",&
                     chm_info%impr_k(nimpr),&
                     chm_info%impr_phi0(nimpr)
                ! Do some units conversion into internal atomic units
                chm_info%impr_phi0(nimpr) = cp_unit_to_cp2k(chm_info%impr_phi0(nimpr),"deg",error=error)
                chm_info%impr_k(nimpr) = cp_unit_to_cp2k(chm_info%impr_k(nimpr),"kcalmol",error=error)
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Read in all the Nonbonded info from the param file here
    !-----------------------------------------------------------------------------
    nnonbond = 0
    nonfo = 0
    DO ilab = 1, SIZE(nbon_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(nbon_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                nnonbond = nnonbond + 1
                CALL reallocate(chm_info%nonbond_a,1,nnonbond)
                CALL reallocate(chm_info%nonbond_eps,1,nnonbond)
                CALL reallocate(chm_info%nonbond_rmin2,1,nnonbond)
                chm_info%nonbond_a(nnonbond) = string
                CALL parser_get_object(parser,chm_info%nonbond_eps(nnonbond),error=error)
                CALL parser_get_object(parser,chm_info%nonbond_eps(nnonbond),error=error)
                CALL parser_get_object(parser,chm_info%nonbond_rmin2(nnonbond),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM NONBOND ",nnonbond,&
                     TRIM(chm_info%nonbond_a(nnonbond))," ",&
                     chm_info%nonbond_eps(nnonbond),&
                     chm_info%nonbond_rmin2(nnonbond)
                chm_info%nonbond_rmin2(nnonbond) = cp_unit_to_cp2k(chm_info%nonbond_rmin2(nnonbond),&
                     "angstrom", error=error)
                chm_info%nonbond_eps(nnonbond) = cp_unit_to_cp2k(chm_info%nonbond_eps(nnonbond),&
                     "kcalmol",error=error)
                IF (parser_test_next_token(parser,error=error) == "FLT") THEN
                   nonfo = nonfo + 1
                   CALL reallocate(chm_info%nonbond_a_14,1,nonfo)
                   CALL reallocate(chm_info%nonbond_eps_14,1,nonfo)
                   CALL reallocate(chm_info%nonbond_rmin2_14,1,nonfo)
                   chm_info%nonbond_a_14(nonfo) = chm_info%nonbond_a(nnonbond)
                   CALL parser_get_object(parser,chm_info%nonbond_eps_14(nonfo),error=error)
                   CALL parser_get_object(parser,chm_info%nonbond_eps_14(nonfo),error=error)
                   CALL parser_get_object(parser,chm_info%nonbond_rmin2_14(nonfo),error=error)
                   IF(iw>0) WRITE(iw,*) "    CHM ONFO ",nonfo,&
                        TRIM(chm_info%nonbond_a_14(nonfo))," ",&
                        chm_info%nonbond_eps_14(nonfo),&
                        chm_info%nonbond_rmin2_14(nonfo)
                   chm_info%nonbond_rmin2_14(nonfo) = cp_unit_to_cp2k(chm_info%nonbond_rmin2_14(nonfo),&
                        "angstrom",error=error)
                   chm_info%nonbond_eps_14(nonfo) = cp_unit_to_cp2k(chm_info%nonbond_eps_14(nonfo),&
                        "kcalmol",error=error)
                END IF
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE read_force_field_charmm

! *****************************************************************************
!> \brief Reads the AMBER force_field
!> \param ff_type ...
!> \param para_env ...
!> \param mm_section ...
!> \param particle_set ...
!> \param error ...
!> \author Teodoro Laino [tlaino, teodoro.laino-AT-gmail.com] - 11.2008
! *****************************************************************************
  SUBROUTINE read_force_field_amber ( ff_type, para_env, mm_section, particle_set, error )

    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_amber', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw
    LOGICAL                                  :: failure
    TYPE(amber_info_type), POINTER           :: amb_info
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    amb_info => ff_type%amb_info

    ! Read the Amber topology file to gather the forcefield information
    CALL rdparm_amber_8(ff_type%ff_file_name, iw, para_env, do_connectivity=.FALSE.,&
         do_forcefield=.TRUE., particle_set=particle_set, amb_info=amb_info, error=error)

    !-----------------------------------------------------------------------------
    ! 1. Converts all the Bonds info from the param file here
    !      Vbond = Kb(b-b0)^2
    !      UNITS for Kb: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
    !-----------------------------------------------------------------------------
    DO i = 1, SIZE(amb_info%bond_a)
       IF(iw>0) WRITE(iw,*) "    AMB BOND ",i,&
            TRIM(amb_info%bond_a(i))," ",&
            TRIM(amb_info%bond_b(i))," ",&
            amb_info%bond_k(i),&
            amb_info%bond_r0(i)

       ! Do some units conversion into internal atomic units
       amb_info%bond_r0(i) = cp_unit_to_cp2k(amb_info%bond_r0(i),"angstrom",error=error)
       amb_info%bond_k(i) = cp_unit_to_cp2k(amb_info%bond_k(i),"kcalmol*angstrom^-2",error=error)
    END DO

    !-----------------------------------------------------------------------------
    ! 2. Converts all the Bends info from the param file here
    !      Vangle = Ktheta(theta-theta0)^2
    !      UNITS for Ktheta: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
    !      FACTOR of "2" rolled into Ktheta
    !      Vub = Kub(S-S0)^2
    !      UNITS for Kub: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
    !-----------------------------------------------------------------------------
    DO i = 1, SIZE(amb_info%bend_a)
       IF(iw>0) WRITE(iw,*) "    AMB BEND ",i,&
            TRIM(amb_info%bend_a(i))," ",&
            TRIM(amb_info%bend_b(i))," ",&
            TRIM(amb_info%bend_c(i))," ",&
            amb_info%bend_k(i),&
            amb_info%bend_theta0(i)

       ! Do some units conversion into internal atomic units
       amb_info%bend_theta0(i) = cp_unit_to_cp2k(amb_info%bend_theta0(i),"rad",error=error)
       amb_info%bend_k(i) = cp_unit_to_cp2k(amb_info%bend_k(i),"kcalmol*rad^-2",error=error)
    END DO

    !-----------------------------------------------------------------------------
    ! 3. Converts all the Dihedrals info from the param file here
    !      Vtorsion = Kphi(1+COS(n(phi)-delta))
    !      UNITS for Kphi: [(kcal/mol)] to [Eh]
    !-----------------------------------------------------------------------------
    DO i = 1, SIZE(amb_info%torsion_a)
       IF(iw>0) WRITE(iw,*) "    AMB TORSION ",i,&
            TRIM(amb_info%torsion_a(i))," ",&
            TRIM(amb_info%torsion_b(i))," ",&
            TRIM(amb_info%torsion_c(i))," ",&
            TRIM(amb_info%torsion_d(i))," ",&
            amb_info%torsion_k(i),&
            amb_info%torsion_m(i),&
            amb_info%torsion_phi0(i)

       ! Do some units conversion into internal atomic units
       amb_info%torsion_phi0(i) = cp_unit_to_cp2k(amb_info%torsion_phi0(i),"rad",error=error)
       amb_info%torsion_k(i) = cp_unit_to_cp2k(amb_info%torsion_k(i),"kcalmol",error=error)
    END DO

    !-----------------------------------------------------------------------------
    ! 4. Converts all the Nonbonded info from the param file here
    !-----------------------------------------------------------------------------
    DO i = 1, SIZE(amb_info%nonbond_eps)
       IF(iw>0) WRITE(iw,*) "    AMB NONBOND ",i,&
            TRIM(amb_info%nonbond_a(i))," ",&
            amb_info%nonbond_eps(i),&
            amb_info%nonbond_rmin2(i)

       ! Do some units conversion into internal atomic units
       amb_info%nonbond_rmin2(i) = cp_unit_to_cp2k(amb_info%nonbond_rmin2(i),"angstrom", error=error)
       amb_info%nonbond_eps(i) = cp_unit_to_cp2k(amb_info%nonbond_eps(i),"kcalmol",error=error)
    END DO
    CALL cp_print_key_finished_output(iw,logger,mm_section,"PRINT%FF_INFO",error=error)
    CALL timestop(handle)
  END SUBROUTINE read_force_field_amber

! *****************************************************************************
!> \brief This function is simply a wrap to the parser_get_next_line..
!>      Comments: This routine would not be necessary if the continuation
!>                char for CHARMM would not be the "-".. How can you choose this
!>                character in a file of numbers as a continuation char????
!>                This sounds simply crazy....
!> \param parser ...
!> \param nline ...
!> \param error ...
!> \author Teodoro Laino - Zurich University - 06.2007
! *****************************************************************************
  SUBROUTINE charmm_get_next_line(parser, nline, error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'charmm_get_next_line', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), PARAMETER              :: continuation_char = "-"

    INTEGER                                  :: i, len_line

    DO i = 1, nline
       len_line = LEN_TRIM(parser%input_line)
       DO WHILE (parser%input_line(len_line:len_line)==continuation_char)
          CALL parser_get_next_line(parser,1,error=error)
          len_line = LEN_TRIM(parser%input_line)
       END DO
       CALL parser_get_next_line(parser,1,error=error)
    END DO

  END SUBROUTINE charmm_get_next_line

END MODULE force_fields_ext
