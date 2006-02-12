!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
!-----------------------------------------------------------------------------!


!!****h* cp2k/leap2fist [1.0] *
!!
!!   NAME
!!     leap2fist
!!
!!   FUNCTION
!!     Program to convert from AMBER parm type to FIST and DL_POLY parm type
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2004 created [tlaino]
!!     teodoro.laino .at. gmail.com 
!!
!!
!!   SOURCE
!!     Please report any bug to the author
!!
!****************************************************************************

PROGRAM leap2fist
  USE kinds,                           ONLY: dp,&
                                             sp
  USE machine,                         ONLY: m_getlog,&
                                             m_hostnm,&
                                             m_getarg,&
                                             m_iargc

  IMPLICIT NONE
  ! 
  ! Parameters Data...
  !
  REAL(KIND=dp), PARAMETER :: CV_CHRG = 18.2223_dp, UNDEF=9999999
  INTEGER, PARAMETER  ::                            &
                        start_io_units        =  10,&
                        max_io_units          = 100,&
                        atom_name_length      =   4,&
                        title_sub_length      =   4,&
                        default_string_length =  80,&
                        title_length          =  80,&
                        date_length           =  25,&
                        input_line_length     = 100,&
                        NDIM                  =   3  ! Cell dimension
  CHARACTER(LEN=date_length)           :: Today
  CHARACTER(LEN=default_string_length) :: ThisIsMe , MyMachine
  REAL(KIND=dp)                        :: pi
  INTEGER :: NATOM        ! Total number of Atoms
  INTEGER :: NTYPES       ! Total number of distinct atom types
  INTEGER :: NBONH        ! Number of bonds containing hydrogens
  INTEGER :: MBONA        ! Number of bonds not containing hydrogens
  INTEGER :: NTHETH       ! Number of angles containing hydrogens
  INTEGER :: MTHETA       ! Number of angles not containing hydrogens
  INTEGER :: NPHIH        ! Number of dihedrals containing hydrogens
  INTEGER :: MPHIA        ! Number of dihedrals not containing hydrogens
  INTEGER :: NHPARM       !    currently NOT USED
  INTEGER :: NPARM        !    set to 1 if LES is used
  INTEGER :: NNB          !    number of excluded atoms
  INTEGER :: NRES         ! Number of residues
  INTEGER :: NBONA        !    MBONA  + number of constraint bonds     ( in v.8 NBONA=MBONA)
  INTEGER :: NTHETA       !    MTHETA + number of constraint angles    ( in v.8 NBONA=MBONA)
  INTEGER :: NPHIA        !    MPHIA  + number of constraint dihedrals ( in v.8 NBONA=MBONA)
  INTEGER :: NUMBND       ! Number of unique bond types
  INTEGER :: NUMANG       ! Number of unique angle types
  INTEGER :: NPTRA        ! Number of unique dihedral types
  INTEGER :: NATYP        ! Number of atom types in parameter file
  INTEGER :: NPHB         ! Number of distinct 10-12 hydrogen bond pair types
  INTEGER :: IFPERT,     &!    Variable not used in this converter...
             NBPER,      &!    Variable not used in this converter...
             NGPER,      &!    Variable not used in this converter...
             NDPER,      &!    Variable not used in this converter...
             MBPER,      &!    Variable not used in this converter...
             MGPER,      &!    Variable not used in this converter...
             MDPER,      &!    Variable not used in this converter...
             IFBOX,      &!    Variable not used in this converter...
             NMXRS,      &!    Variable not used in this converter...
             IFCAP,      &!    Variable not used in this converter...
             NUMEXTRA,   &!    Variable not used in this converter...
             NCOPY        !    Variable not used in this converter...
  CHARACTER (LEN=title_length)                            :: TITLE 
  CHARACTER (LEN=atom_name_length), POINTER, DIMENSION(:) :: IGRAPH,&
                                                             LABRES,&
                                                             ISYMBL
  INTEGER, DIMENSION(:), POINTER                          :: IAC,&
                                                             ICO,&
                                                             IPRES,&
                                                             IBH,JBH,ICBH,&
                                                             IB,JB,ICB,&
                                                             ITH,JTH,KTH,ICTH,&
                                                             IT,JT,KT,ICT,&
                                                             IPH,JPH,KPH,LPH,ICPH,&
                                                             IP,JP,KP,LP,ICP
  REAL (KIND=dp), DIMENSION(:), POINTER                   :: CHRG,&
                                                             AMASS,&
                                                             RK,&
                                                             REQ,&
                                                             TK,&
                                                             TEQ,&
                                                             PK,&
                                                             PN,&
                                                             PHASE,&
                                                             CN1,&
                                                             CN2,&
                                                             ASOL,&
                                                             BSOL,&
                                                             BOX
  REAL (KIND=dp)                                          :: BETA
  !
  ! Local Variables...
  !                    
  CHARACTER (LEN=default_string_length), POINTER, DIMENSION(:) :: command_line 
  CHARACTER (LEN=default_string_length) :: input_filename, output_filename, psf_filename
  INTEGER :: unit_o, unit_i, unit_p
  LOGICAL :: verbose, dlpoly, fist, psf_xplor, amber_impropers
  LOGICAL, POINTER, DIMENSION(:) :: io_units
  CHARACTER(LEN=5), POINTER, DIMENSION(:) :: molnames

  NULLIFY( command_line, IGRAPH, ISYMBL, LABRES, IAC, ICO, IPRES,IBH,JBH,ICBH,&
           IB,JB,ICB,ITH,JTH,KTH,ICTH,IT,JT,KT,ICT,IPH,JPH,KPH,LPH,ICPH,&
           IP,JP,KP,LP,ICP,CHRG,RK,REQ,TK,TEQ,PK,PN,PHASE,CN1,CN2,ASOL,BSOL,&
           BOX, AMASS, io_units, molnames)
  CALL get_command_line(command_line)
  !
  ! Start conversion...
  !
  CALL open_file(input_filename,  "old", "formatted", unit_i)
  CALL open_file(output_filename, "unknown", "formatted", unit_o)
  CALL rdparm_amber_8(unit_i)
  IF (.NOT.dlpoly) CALL conform_amber_2_charmm
  CALL initialize_var
  IF (fist) THEN
     CALL open_file(psf_filename, "unknown", "formatted", unit_p)
     !
     ! Generate the potential parameter file...
     !
     CALL create_potpar_file(unit_o)
     !
     ! Generate the PSF file...
     !
     CALL generate_molecules()
     CALL write_psf(unit_p)
     CALL close_file(unit_p)
  ELSEIF (dlpoly) THEN
     CALL write_dlpoly_field(unit_o)
  END IF
  CALL close_file(unit_o)
  CALL close_file(unit_i)
  CALL release_amber_8_structures

CONTAINS

  SUBROUTINE initialize_var
    IMPLICIT NONE
    CHARACTER (LEN=10)        :: time
    CHARACTER (LEN=5)         :: zone
    CHARACTER (LEN=8)         :: date
    INTEGER, DIMENSION(8)     :: values

    pi = ATAN(1.0_dp)*4.0_dp
    CALL DATE_AND_TIME(date=date, time=time, zone=zone, values=values)
    Today=TRIM(date)//" "//TRIM(time)//" "//TRIM(zone)
    CALL m_getlog(ThisIsMe)    
    CALL m_hostnm(MyMachine)
    
  END SUBROUTINE initialize_var
    
  REAL(KIND=dp) FUNCTION CONVERT_ANGLE(radiant, degree) RESULT(conversion)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: radiant, degree

    IF (PRESENT(radiant)) THEN
       conversion = radiant / pi * 180.0_dp
    ELSE IF (PRESENT(degree)) THEN
       conversion = degree  / 180.0_dp * pi
    END IF
  END FUNCTION CONVERT_ANGLE
    

  SUBROUTINE open_file(filename, status, form, my_unit)
    IMPLICIT NONE
    INTEGER, SAVE :: ifirst
    DATA ifirst /0/
    INTEGER, INTENT(OUT) :: my_unit
    CHARACTER (LEN=*), INTENT(IN) :: filename, status, form
    INTEGER :: I
    
    IF (ifirst == 0) THEN
       ALLOCATE(io_units(start_io_units:max_io_units))
       io_units = .FALSE.
       ifirst = ifirst + 1
    END IF
    DO I = start_io_units, max_io_units
       IF (.not.io_units(i)) EXIT
    END DO
    IF (i.GT.max_io_units) CALL stop_converter("Exceeded maximum number of contemporary opened units!")
    my_unit = i
    io_units(my_unit) = .TRUE.
    OPEN (UNIT=my_unit, file=filename, status=status, form=form)
  END SUBROUTINE open_file

  SUBROUTINE close_file(my_unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: my_unit

    io_units(my_unit) = .FALSE.
    CLOSE(my_unit)
  END SUBROUTINE close_file

  SUBROUTINE print_help_banner
    IMPLICIT NONE
    
    WRITE(*,'(A)')"Usage: leap2fist -top filename -out filename "&
                            //" [options]"
    WRITE(*,'(A)')"Options:"
    WRITE(*,'(2X,A,T20,A)')"-top","Topology file generated with leap from AMBER V.8"
    WRITE(*,'(2X,A,T20,A)')"-pot","Optional. Specify the file name of the potential file"
    WRITE(*,'(2X,A,T20,A)')"    ","Default is the top file's root name"
    WRITE(*,'(2X,A,T20,A)')"-psf","Optional. Specify the file name of the psf file"
    WRITE(*,'(2X,A,T20,A)')"    ","Default is the top file's root name"
    WRITE(*,'(2X,A,T20,A)')"-xplor","Optional. PSF will be dumped out according the X-PLOR format"
    WRITE(*,'(2X,A,T20,A)')"    ","This format is also readable by VMD and can be used within NAMD."
    WRITE(*,'(2X,A,T20,A)')"-impropers","Optional. Put impropers in a different section w.r.t. proper torsions."
    WRITE(*,'(2X,A,T20,A)')"    ","Default = .FALSE., i.e. impropers and propers are in the same section (AMBER style)."
    WRITE(*,'(2X,A,T20,A)')"-verbose","Print verbose information during executions" 
    WRITE(*,'(2X,A,T20,A)')"-dlpoly","Converts parameter file into dlpoly format." 
    WRITE(*,'(2X,A,T20,A)')"    ","    "
    WRITE(*,'(2X,A,T20,A)')"Comments:",&
         "This program converts an AMBER PARMTOP file into a CHARM PARM file to be used with Fist."

    CALL stop_converter(" ")
  END SUBROUTINE print_help_banner

  SUBROUTINE get_command_line(command_line)
    IMPLICIT NONE
    CHARACTER (LEN=default_string_length), POINTER, DIMENSION(:) :: command_line
    INTEGER :: narg, ind_root, i

    input_filename  = "NULL"
    output_filename = "NULL"
    psf_filename    = "NULL"
    verbose         = .FALSE.
    dlpoly          = .FALSE.
    fist            = .TRUE.
    psf_xplor       = .FALSE.
    amber_impropers = .TRUE.
    narg = m_iargc()
    IF (narg.LT.2) THEN
       CALL print_help_banner
    END IF
    ALLOCATE(command_line(narg))
    DO i=1, narg
       CALL m_getarg(i, command_line(i))
    END DO
    DO I = 1, SIZE(command_line)
       SELECT CASE(TRIM(command_line(i)))
       CASE ("-top")
          input_filename   = TRIM(command_line(i+1))
       CASE ("-pot") 
          output_filename  = TRIM(command_line(i+1))
       CASE ("-psf") 
          psf_filename  = TRIM(command_line(i+1))
       CASE ("-verbose")
          verbose = .TRUE.
       CASE ("-impropers")
          amber_impropers = .FALSE.
       CASE ("-xplor")
          psf_xplor = .TRUE.
       CASE ("-dlpoly")
          dlpoly = .TRUE.
          fist   = .FALSE.
       CASE DEFAULT
          ! just do nothing...
       END SELECT
    END DO
    !
    ind_root = INDEX(input_filename,".")
    IF (ind_root == 0) ind_root = LEN_TRIM(input_filename)
    IF(INDEX(output_filename,"NULL") /= 0 )&
         output_filename = TRIM(input_filename(1:ind_root)//"pot")
    IF(INDEX(psf_filename,"NULL") /= 0 )&
         psf_filename = TRIM(input_filename(1:ind_root)//"psf")
    IF (dlpoly) output_filename = "FIELD"
    IF (dlpoly) psf_filename    = "NULL"
  END SUBROUTINE get_command_line

  SUBROUTINE rdparm_amber_8(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: unit
    CHARACTER (LEN=default_string_length)  :: fmt, section    
    LOGICAL                                :: failure
    INTEGER                                :: i

    failure = check_amber_8_std(unit)
    IF (.NOT.failure) THEN
       DO WHILE (get_section_parmtop(unit, section, fmt))
          SELECT CASE (TRIM(section))
          CASE("TITLE")
             READ(unit,fmt)TITLE
          CASE("POINTERS")
             READ(unit,fmt)NATOM, NTYPES, NBONH,    MBONA, NTHETH, MTHETA, NPHIH,  &
                           MPHIA, NHPARM, NPARM,    NNB,   NRES,   NBONA,  NTHETA, &
                           NPHIA, NUMBND, NUMANG,   NPTRA, NATYP,  NPHB,   IFPERT, &
                           NBPER, NGPER,  NDPER,    MBPER, MGPER,  MDPER,  IFBOX,  &
                           NMXRS, IFCAP,  NUMEXTRA, NCOPY                           
             WRITE(*,1000) NATOM, NTYPES, NBONH,    MBONA, NTHETH, MTHETA, NPHIH,  &
                           MPHIA, NHPARM, NPARM,    NNB,   NRES,   NBONA,  NTHETA, &
                           NPHIA, NUMBND, NUMANG,   NPTRA, NATYP,  NPHB,   IFBOX,  &
                           NMXRS, IFCAP,  NUMEXTRA, NCOPY                   

             ALLOCATE( IGRAPH (NATOM              ),&
                       CHRG   (NATOM              ),&
                       AMASS  (NATOM              ),&
                       IAC    (NATOM              ),&
                       ICO    (NTYPES*NTYPES      ),&
                       LABRES (NRES               ),&
                       IPRES  (NRES               ),&
                       RK     (NUMBND             ),&
                       REQ    (NUMBND             ),&
                       TK     (NUMANG             ),&
                       TEQ    (NUMANG             ),&
                       PK     (NPTRA              ),&
                       PN     (NPTRA              ),&
                       PHASE  (NPTRA              ),&
                       CN1    (NTYPES*(NTYPES+1)/2),&
                       CN2    (NTYPES*(NTYPES+1)/2),&
                       ASOL   (NTYPES*(NTYPES+1)/2),&
                       BSOL   (NTYPES*(NTYPES+1)/2),&
                       IBH    (NBONH              ),&
                       JBH    (NBONH              ),&
                       ICBH   (NBONH              ),&
                       IB     (NBONA              ) )
             ALLOCATE( JB     (NBONA              ),&
                       ICB    (NBONA              ),&
                       ITH    (NTHETH             ),&
                       JTH    (NTHETH             ),&
                       KTH    (NTHETH             ),&
                       ICTH   (NTHETH             ),&
                       IT     (NTHETA             ),&
                       JT     (NTHETA             ),&
                       KT     (NTHETA             ),&
                       ICT    (NTHETA             ),&
                       IPH    (NPHIH              ),&
                       JPH    (NPHIH              ),&
                       KPH    (NPHIH              ),&
                       LPH    (NPHIH              ),&
                       ICPH   (NPHIH              ),&
                       IP     (NPHIA              ),&
                       JP     (NPHIA              ),&
                       KP     (NPHIA              ),&
                       LP     (NPHIA              ),&
                       ICP    (NPHIA              ),&
                       ISYMBL (NATOM              ),&
                       BOX    (NDIM               ) &
                     )  
          CASE("ATOM_NAME")
             READ(unit,fmt)(IGRAPH(i),                                i=1,NATOM)
          CASE("CHARGE")                                              
             READ(unit,fmt)(CHRG(i),                                  i=1,NATOM)
          CASE("MASS")                                              
             READ(unit,fmt)(AMASS(i),                                  i=1,NATOM)
          CASE("ATOM_TYPE_INDEX")                                     
             READ(unit,fmt)(IAC(i),                                   i=1,NATOM)
          CASE("NONBONDED_PARM_INDEX")                                
             READ(unit,fmt)(ICO(i),                                   i=1,NTYPES*NTYPES)
          CASE("RESIDUE_LABEL")                                       
             READ(unit,fmt)(LABRES(i),                                i=1,NRES)
          CASE("RESIDUE_POINTER")                                     
             READ(unit,fmt)(IPRES(i),                                 i=1,NRES)
          CASE("BOND_FORCE_CONSTANT")                                 
             READ(unit,fmt)(RK(i),                                    i=1,NUMBND)
          CASE("BOND_EQUIL_VALUE")                                    
             READ(unit,fmt)(REQ(i),                                   i=1,NUMBND)
          CASE("ANGLE_FORCE_CONSTANT")                                
             READ(unit,fmt)(TK(i),                                    i=1,NUMANG)
          CASE("ANGLE_EQUIL_VALUE")                                   
             READ(unit,fmt)(TEQ(i),                                   i=1,NUMANG)
          CASE("DIHEDRAL_FORCE_CONSTANT")                             
             READ(unit,fmt)(PK(i),                                    i=1,NPTRA)
          CASE("DIHEDRAL_PERIODICITY")                                
             READ(unit,fmt)(PN(i),                                    i=1,NPTRA)
          CASE("DIHEDRAL_PHASE")                                      
             READ(unit,fmt)(PHASE(i),                                 i=1,NPTRA)
          CASE("LENNARD_JONES_ACOEF")                                 
             READ(unit,fmt)(CN1(i),                                   i=1,NTYPES*(NTYPES+1)/2)
          CASE("LENNARD_JONES_BCOEF")                                 
             READ(unit,fmt)(CN2(i),                                   i=1,NTYPES*(NTYPES+1)/2)
          CASE("HBOND_ACOEF")
             READ(unit,fmt)(ASOL(i),                                  i=1,NPHB)
          CASE("HBOND_BCOEF")
             READ(unit,fmt)(BSOL(i),                                  i=1,NPHB)
          CASE("BONDS_INC_HYDROGEN")                                  
             READ(unit,fmt)(IBH(i),JBH(i),ICBH(i),                    i=1,NBONH)
          CASE("BONDS_WITHOUT_HYDROGEN")                              
             READ(unit,fmt)(IB(i),JB(i),ICB(i),                       i=1,NBONA)
          CASE("ANGLES_INC_HYDROGEN")
             READ(unit,fmt)(ITH(i),JTH(i),KTH(i),ICTH(i),             i=1,NTHETH)
          CASE("ANGLES_WITHOUT_HYDROGEN")
             READ(unit,fmt)(IT(i),JT(i),KT(i),ICT(i),                 i=1,NTHETA)
          CASE("DIHEDRALS_INC_HYDROGEN")
             READ(unit,fmt)(IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i),      i=1,NPHIH)
          CASE("DIHEDRALS_WITHOUT_HYDROGEN")
             READ(unit,fmt)(IP(i),JP(i),KP(i),LP(i),ICP(i),           i=1,NPHIA)
          CASE("AMBER_ATOM_TYPE")
             READ(unit,fmt)(ISYMBL(i),                                i=1,NATOM)
          CASE("BOX_DIMENSIONS")
             READ(unit,fmt)BETA,(BOX(i),                              i=1,NDIM)
          CASE DEFAULT
             ! Just Ignore other sections...
          END SELECT
       END DO
    END IF

    RETURN
    !
    ! FORMAT SECTION
    !
1000 FORMAT(T2, &
          /' NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7, &
          /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7, &
          /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7, &
          /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7, &
          /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7, &
          /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7, &
          /' NCOPY  = ',i7/)

  END SUBROUTINE rdparm_amber_8

  SUBROUTINE  conform_amber_2_charmm
    IMPLICIT NONE
    INTEGER :: i, num_types, num_types2, j, k, len_type, l
    CHARACTER(LEN=atom_name_length) :: my_type, my_form
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: my_charge
    INTEGER,  DIMENSION(:), ALLOCATABLE :: my_index, my_done
    REAL(KIND=dp) :: my_ch
    !
    ! Modify AMBER TYPES if two symbols have different charge
    ! This is not allowed in CHARMM !
    !
    ALLOCATE(my_done(NATOM))
    my_done = 0
    DO i = 1, NATOM
       my_type    = ISYMBL(i)
       len_type   = LEN_TRIM(my_type)
       num_types  = COUNT(ISYMBL==my_type)
       IF (num_types == 1.OR.my_done(i)==1) THEN 
          my_done(i) = 1
          CYCLE
       END IF
       my_done(i) = 1
       ALLOCATE(my_charge(num_types), my_index(num_types))
       k = 0
       ! Collect Information
       DO j = 1, NATOM
          IF (ISYMBL(j) /= my_type) CYCLE
          k = k + 1
          my_charge(k) = CHRG(j)
          my_index(k)  = j
       END DO
       ! Post-Processing of information       
       IF (ALL(my_charge==my_charge(1))) THEN 
          DEALLOCATE(my_charge, my_index)
          CYCLE
       END IF       
       k = 1
       DO j = 1, num_types
          my_ch = my_charge(j)
          num_types2 = COUNT(my_charge==my_ch)
          IF ( my_ch==-9999 ) THEN 
             CYCLE
          ELSEIF(num_types2 < 10) THEN
             IF (len_type>3) CALL stop_converter(" A type with more than 4 character should be created!")
          ELSEIF(num_types2 < 100) THEN
             IF (len_type>2) CALL stop_converter(" A type with more than 4 character should be created!")
          ELSEIF(num_types2 < 1000) THEN
             IF (len_type>1) CALL stop_converter(" A type with more than 4 character should be created!")
          ELSE
             CALL stop_converter(" A type with more than 4 character should be created!")
          END IF
          DO l = 1, num_types
             IF (my_charge(l) /= my_ch) CYCLE
             IF (k<10) THEN
                my_form = '(I1)'
             ELSEIF (K<100) THEN
                my_form = '(I2)'
             ELSEIF (K<1000) THEN
                my_form = '(I3)'
             ELSE
                CALL stop_converter(" Should never reach this point!")
             END IF
             WRITE(ISYMBL(my_index(l))(LEN_TRIM(ISYMBL(my_index(l)))+1:),my_form)k
             my_charge(l) = -9999
             my_done(my_index(l)) = 1
          END DO
          k = k + 1
       END DO
       IF (ANY(my_charge /= -9999)) CALL stop_converter(" Should never reach this point!")
       DEALLOCATE(my_charge, my_index)
    END DO
    IF (ANY(my_done/=1)) CALL stop_converter(" Should never reach this point!")
    DEALLOCATE(my_done)
  END SUBROUTINE conform_amber_2_charmm

  SUBROUTINE release_amber_8_structures
    IMPLICIT NONE
    IF (ASSOCIATED( command_line)) DEALLOCATE ( command_line)
    IF (ASSOCIATED( io_units)) DEALLOCATE( io_units)
    IF (ASSOCIATED( IGRAPH ))  DEALLOCATE( IGRAPH )
    IF (ASSOCIATED( CHRG   ))  DEALLOCATE( CHRG   )
    IF (ASSOCIATED( AMASS  ))  DEALLOCATE( AMASS  )
    IF (ASSOCIATED( IAC    ))  DEALLOCATE( IAC    )
    IF (ASSOCIATED( ICO    ))  DEALLOCATE( ICO    )
    IF (ASSOCIATED( LABRES ))  DEALLOCATE( LABRES )
    IF (ASSOCIATED( IPRES  ))  DEALLOCATE( IPRES  )
    IF (ASSOCIATED( RK     ))  DEALLOCATE( RK     )
    IF (ASSOCIATED( REQ    ))  DEALLOCATE( REQ    )
    IF (ASSOCIATED( TK     ))  DEALLOCATE( TK     )
    IF (ASSOCIATED( TEQ    ))  DEALLOCATE( TEQ    )
    IF (ASSOCIATED( PK     ))  DEALLOCATE( PK     )
    IF (ASSOCIATED( PN     ))  DEALLOCATE( PN     )
    IF (ASSOCIATED( PHASE  ))  DEALLOCATE( PHASE  )
    IF (ASSOCIATED( CN1    ))  DEALLOCATE( CN1    )
    IF (ASSOCIATED( CN2    ))  DEALLOCATE( CN2    )
    IF (ASSOCIATED( ASOL   ))  DEALLOCATE( ASOL   )
    IF (ASSOCIATED( BSOL   ))  DEALLOCATE( BSOL   )
    IF (ASSOCIATED( IBH    ))  DEALLOCATE( IBH    )
    IF (ASSOCIATED( JBH    ))  DEALLOCATE( JBH    )
    IF (ASSOCIATED( ICBH   ))  DEALLOCATE( ICBH   )
    IF (ASSOCIATED( IB     ))  DEALLOCATE( IB     )
    IF (ASSOCIATED( JB     ))  DEALLOCATE( JB     )
    IF (ASSOCIATED( ICB    ))  DEALLOCATE( ICB    )
    IF (ASSOCIATED( ITH    ))  DEALLOCATE( ITH    )
    IF (ASSOCIATED( JTH    ))  DEALLOCATE( JTH    )
    IF (ASSOCIATED( KTH    ))  DEALLOCATE( KTH    )
    IF (ASSOCIATED( ICTH   ))  DEALLOCATE( ICTH   )
    IF (ASSOCIATED( IT     ))  DEALLOCATE( IT     )
    IF (ASSOCIATED( JT     ))  DEALLOCATE( JT     )
    IF (ASSOCIATED( KT     ))  DEALLOCATE( KT     )
    IF (ASSOCIATED( ICT    ))  DEALLOCATE( ICT    )
    IF (ASSOCIATED( IPH    ))  DEALLOCATE( IPH    )
    IF (ASSOCIATED( JPH    ))  DEALLOCATE( JPH    )
    IF (ASSOCIATED( KPH    ))  DEALLOCATE( KPH    )
    IF (ASSOCIATED( LPH    ))  DEALLOCATE( LPH    )
    IF (ASSOCIATED( ICPH   ))  DEALLOCATE( ICPH   )
    IF (ASSOCIATED( IP     ))  DEALLOCATE( IP     )
    IF (ASSOCIATED( JP     ))  DEALLOCATE( JP     )
    IF (ASSOCIATED( KP     ))  DEALLOCATE( KP     )
    IF (ASSOCIATED( LP     ))  DEALLOCATE( LP     )
    IF (ASSOCIATED( ICP    ))  DEALLOCATE( ICP    )
    IF (ASSOCIATED( ISYMBL ))  DEALLOCATE( ISYMBL )
    IF (ASSOCIATED( BOX    ))  DEALLOCATE( BOX    )
  END SUBROUTINE release_amber_8_structures


  SUBROUTINE stop_converter(msg)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: msg

    CALL release_amber_8_structures
    IF (LEN_TRIM(msg) == 0) STOP
    WRITE(*,'(A)')msg
    WRITE(*,'(A)')"Execution aborted... stopping now!"
    STOP
  END SUBROUTINE stop_converter


  LOGICAL FUNCTION  get_section_parmtop(unit, section, fmt) RESULT(another_section)
    IMPLICIT NONE
    CHARACTER (LEN=default_string_length), INTENT(OUT) :: fmt, section
    INTEGER, INTENT(IN)                                :: unit
    CHARACTER (LEN=input_line_length)                  :: input_line
    INTEGER :: indflag, start_f, end_f
    
    ! section
    READ(unit,'(A)',END=100)input_line
    DO WHILE (INDEX(input_line,"%FLAG") == 0)
       READ(unit,'(A)',END=100)input_line
    END DO
    IF (verbose) WRITE(*,*)"Reading from PrmTop  ::"//input_line
    indflag = INDEX(input_line,"%FLAG")+LEN_TRIM("%FLAG")
    DO WHILE (INDEX(input_line(indflag:indflag)," ") /= 0)
       indflag = indflag + 1
    END DO
    section = TRIM(input_line(indflag:))
    IF (verbose) WRITE(*,*)"Reading section name ::"//TRIM(section)
    ! format
    READ(unit,'(A)')input_line
    IF (INDEX(input_line,"%FORMAT") == 0) THEN
       CALL stop_converter(" Expecting %FORMAT.. :: "//"FATAL ERROR !")
    END IF
    start_f = INDEX(input_line,"(")
    end_f   = INDEX(input_line,")")
    fmt     = input_line(start_f:end_f)
    IF (verbose) WRITE(*,*)"Format specified     :: "//TRIM(fmt)
    another_section = .TRUE.
    RETURN
100 another_section = .FALSE.
  END FUNCTION get_section_parmtop

  LOGICAL FUNCTION  check_amber_8_std(unit) RESULT(failure_AMBER_V8)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                :: unit
    CHARACTER (LEN=default_string_length)              :: line

    REWIND(unit)
    READ(unit,'(A)')line
    failure_AMBER_V8 = .TRUE.
    IF (INDEX(line,"%VERSION ") /= 0) THEN
       failure_AMBER_V8 = .FALSE.
       WRITE(*,'(A)')"Amber PrmTop V.8 :: "//TRIM(line)
    END IF
    IF (verbose) WRITE(*,*)" This is an AMBER V.8 PRMTOP format file :: ",.NOT.failure_AMBER_V8
  END FUNCTION check_amber_8_std


  SUBROUTINE generate_molecules
    IMPLICIT NONE
    INTEGER, DIMENSION(:), POINTER :: my_i, my_j
    INTEGER :: my_bonds, ibond, istart, jstart, imol, i, k
    CHARACTER (LEN=10), POINTER, DIMENSION(:) :: map_mol_name
    CHARACTER (LEN=10) :: MYNAME

    NULLIFY(map_mol_name)
    my_bonds = NBONH+MBONA
    ALLOCATE(my_i(my_bonds), my_j(my_bonds))
    my_i(1:NBONH) = IBH/3+1; my_i(NBONH+1:my_bonds) = IB/3+1
    my_j(1:NBONH) = JBH/3+1; my_j(NBONH+1:my_bonds) = JB/3+1    
    
    !Zero the arrays
    ALLOCATE(map_mol_name(NATOM))
    map_mol_name(:) = "MOL0000000"
    ibond = 0
    imol  = 0
    DO WHILE ((COUNT(MASK=my_i>0)+COUNT(MASK=my_j>0))>0)
       ibond = ibond + 1
       istart = my_i(ibond)
       jstart = my_j(ibond)
       if ((istart == 0).AND.(jstart==0)) CYCLE
       imol = imol + 1
       MYNAME = give_mol_name("MOL0000000",imol)
       my_i(ibond) = 0
       my_j(ibond) = 0
       map_mol_name(istart) = MYNAME
       map_mol_name(jstart) = MYNAME

       CALL build_mol_low(map_mol_name, MYNAME, istart, my_i, my_j)
       CALL build_mol_low(map_mol_name, MYNAME, jstart, my_i, my_j)
    END DO
    !
    !  Identify atoms with no bonds...
    !
    my_i(1:NBONH) = IBH/3+1; my_i(NBONH+1:my_bonds) = IB/3+1
    my_j(1:NBONH) = JBH/3+1; my_j(NBONH+1:my_bonds) = JB/3+1
    DO i = 1, NATOM
       IF ((COUNT(MASK=my_i==i)+COUNT(MASK=my_j==i)) == 0) THEN
          imol = imol + 1
          MYNAME = give_mol_name("MOL0000000",imol)
          map_mol_name(i) = MYNAME
       END IF
    END DO
    DEALLOCATE(my_i, my_j)
    CALL clean_mol_name(map_mol_name,imol)
    IF (verbose) WRITE(*,'(8A10)')map_mol_name
    ALLOCATE(molnames(NATOM))
    DO i=1, NATOM
       molnames(i) = "MOL00"
       DO k = 4, 10
          IF (map_mol_name(i)(k:k) /= "0") THEN
             IF (k==9) THEN
                molnames(i)(4:5) = map_mol_name(i)(k:10)
             ELSEIF (k == 10) THEN
                molnames(i)(5:5) = map_mol_name(i)(k:10)
             ELSE
                WRITE(*,'(A)')"More than 99 different molecules.. Wee need to implement a special case for this!"
                CALL stop_converter("Error in generate_molecules: More than 99 DIFFERENT molecules present!")
             END IF
             EXIT 
          END IF
       END DO
    END DO
    DEALLOCATE(map_mol_name)
  END SUBROUTINE generate_molecules

  RECURSIVE SUBROUTINE build_mol_low(map_mol_name, MYNAME, index, my_i, my_j)
    IMPLICIT NONE
    CHARACTER (LEN=10), POINTER, DIMENSION(:) :: map_mol_name
    CHARACTER (LEN=10), INTENT(IN) :: MYNAME
    INTEGER, DIMENSION(:), POINTER :: my_i, my_j
    INTEGER, INTENT(IN) :: index
    INTEGER :: istart, jstart, my_bond

    DO WHILE (COUNT(MASK=my_i==index)>0)
       my_bond = 1
       DO WHILE (my_bond <= SIZE(my_i))
          IF (my_i(my_bond) == index) THEN
             istart = my_i(my_bond)
             jstart = my_j(my_bond)
             my_i(my_bond)  = 0
             my_j(my_bond)  = 0
             map_mol_name(istart) = MYNAME
             map_mol_name(jstart) = MYNAME
             CALL build_mol_low(map_mol_name, MYNAME, istart, my_i, my_j)
             CALL build_mol_low(map_mol_name, MYNAME, jstart, my_i, my_j)          
          END IF
          my_bond = my_bond + 1
       END DO
    END DO

    DO WHILE (COUNT(MASK=my_j==index)>0)
       my_bond = 1
       DO WHILE (my_bond <= SIZE(my_j))
          IF (my_j(my_bond) == index) THEN
             istart = my_i(my_bond)
             jstart = my_j(my_bond)
             my_i(my_bond)  = 0
             my_j(my_bond)  = 0
             map_mol_name(istart) = MYNAME
             map_mol_name(jstart) = MYNAME
             CALL build_mol_low(map_mol_name, MYNAME, istart, my_i, my_j)
             CALL build_mol_low(map_mol_name, MYNAME, jstart, my_i, my_j)          
          END IF
          my_bond = my_bond + 1
       END DO
    END DO

  END SUBROUTINE build_mol_low

  FUNCTION give_mol_name(basename,imol) RESULT(my_name)
    IMPLICIT NONE
    CHARACTER (LEN=10), INTENT(IN)  :: basename
    CHARACTER (LEN=10) :: my_name
    INTEGER, INTENT(IN) :: imol
    integer :: I
   
    my_name = basename
    WRITE(my_name(4:10),'(I7)')imol
    DO I = 4, 10
       IF (my_name(I:I)==" ") my_name(I:I)="0"
    END DO
  END FUNCTION give_mol_name

  SUBROUTINE clean_mol_name(map_mol_name, nmols)
    IMPLICIT NONE
    CHARACTER (LEN=10), POINTER, DIMENSION(:) :: map_mol_name
    INTEGER, INTENT(IN) :: nmols
    INTEGER, POINTER, DIMENSION(:) ::  ifirst, ilast, idim
    INTEGER :: my_first, my_last, I, J, imol, my_dim, my_res
    CHARACTER (LEN=10) :: MY_NAME    
    CHARACTER (LEN=default_string_length), DIMENSION (:), POINTER :: MY_LABRES

    ALLOCATE(ifirst(nmols), ilast(nmols), idim(nmols), MY_LABRES(NATOM))
    !
    DO my_res = 1, NRES
       my_first = IPRES(my_res)
       IF (my_res == NRES) THEN 
          my_last = NATOM
       ELSE
          my_last = IPRES(my_res + 1) - 1
       END IF
       MY_LABRES(my_first:my_last) = LABRES(my_res)
    END DO
    !
    my_first  = 1
    ifirst(1) = 1
    my_first  = 1
    DO i = 1, nmols
       DO j = my_first+1, NATOM
          IF (map_mol_name(j) /= map_mol_name(my_first)) THEN
             my_last  = j - 1
             my_first = j
             ilast(i) = my_last
             idim(i)  = ilast(i) - ifirst(i) + 1
             IF (i /=nmols) ifirst(i+1) = my_first
             EXIT
          END IF
       END DO
       IF (j == NATOM+1) THEN
          IF (i /= nmols) CALL stop_converter("Error in clean_mol_name :: Something seems to be unlogical...")
          ilast(i) = NATOM
          idim(i)  = ilast(i) - ifirst(i) + 1
       END IF
    END DO
    IF (verbose) WRITE(*,'(A)')"IFIRST ::"
    IF (verbose) WRITE(*,*)ifirst
    IF (verbose) WRITE(*,'(A)')"ILAST  ::"
    IF (verbose) WRITE(*,*)ilast
    imol = 0
    DO j = 1, nmols
       IF ( idim(j) == 0 ) CYCLE
       imol = imol + 1
       MY_NAME = give_mol_name("MOL0000000",imol)
       map_mol_name(ifirst(j):ilast(j)) = MY_NAME
       my_dim  = idim(j)
       idim(j) = 0
       DO i = 1, nmols
          IF ( idim(i) == 0 ) CYCLE
          IF ( idim(i) == my_dim ) THEN
             ! check residues...
             IF (ALL(MY_LABRES(ifirst(i):ilast(i)) == MY_LABRES(ifirst(j):ilast(j)))) THEN
                ! Molecules are really the same
                IF (verbose) WRITE(*,'(A)')"Find similar molecules",MY_LABRES(ifirst(i):ilast(i))
                IF (verbose) WRITE(*,'(A)')"Find similar molecules",MY_LABRES(ifirst(j):ilast(j))
                map_mol_name(ifirst(i):ilast(i)) = MY_NAME
                idim(i) = 0
             END IF
          END IF
       END DO
    END DO

    DEALLOCATE(ifirst, ilast, MY_LABRES)
  END SUBROUTINE clean_mol_name

  SUBROUTINE write_psf(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER :: my_sec, my_res, my_bonds, local, my_phi, my_theta, my_atom
    INTEGER :: i, ind, ifirst, ilast
    INTEGER, DIMENSION(:), POINTER :: my_i, my_j, my_k, my_l
    CHARACTER (LEN=2) :: int_format

    NULLIFY( my_i, my_j, my_k, my_l)
    IF(verbose) WRITE(*,'(A)') "    Writing out PSF file "//TRIM(psf_filename)    
    !
    ! Section Title
    !
    my_sec = 1
    int_format="I7"
    if (psf_xplor) int_format="I8"
    WRITE(unit,'(A)') "PSF"
    WRITE(unit,'(A)') " "
    WRITE(unit,'(I7,A)') my_sec, " !NTITLE"
    WRITE(unit,'(A)')" Conversion from AMBER PARMTOP ::"//TRIM(TITLE)
    WRITE(unit,'(A)') " "
    
    !
    ! Section Atom
    !
    WRITE(unit,'('//int_format//',A)')NATOM," !NATOM"
    DO my_res = 1, NRES
       ifirst = IPRES(my_res)
       IF (my_res == NRES) THEN 
          ilast = NATOM
       ELSE
          ilast = IPRES(my_res + 1) - 1
       END IF
       DO my_atom = ifirst, ilast
          WRITE(unit,'(I7,1X,A5,I7,1X,A5,A5,A5,F10.5,F10.5,I7)') &
               my_atom,&
               molnames(my_atom),&
               my_res ,&
               LABRES(my_res) ,&
               IGRAPH(my_atom),&
               ISYMBL(my_atom),&
               CHRG(my_atom)/CV_CHRG ,&
               AMASS(my_atom) ,0
       END DO
    END DO

    WRITE(unit,*) " "
    !
    ! Section Bonds
    !
    my_bonds = NBONH+MBONA
    ALLOCATE(my_i(my_bonds), my_j(my_bonds))
    my_i(1:NBONH) = IBH/3+1; my_i(NBONH+1:my_bonds) = IB/3+1
    my_j(1:NBONH) = JBH/3+1; my_j(NBONH+1:my_bonds) = JB/3+1
    WRITE(unit,*) " ",my_bonds," !NBOND"
    WRITE(unit,'(8'//int_format//')')(my_i(local),my_j(local), local=1,my_bonds)
    DEALLOCATE(my_i, my_j)
    WRITE(unit,*) ""
    !
    ! Section Theta
    !
    my_theta = NTHETH + MTHETA
    ALLOCATE(my_i(my_theta), my_j(my_theta), my_k(my_theta))
    my_i(1:NTHETH) = ITH/3+1; my_i(NTHETH+1:my_theta) = IT/3+1
    my_j(1:NTHETH) = JTH/3+1; my_j(NTHETH+1:my_theta) = JT/3+1
    my_k(1:NTHETH) = KTH/3+1; my_k(NTHETH+1:my_theta) = KT/3+1
    WRITE(unit,*) " ",my_theta," !NTHETA"
    WRITE(unit,'(9'//int_format//')') (my_i(local),my_j(local),my_k(local), local=1,my_theta)
    WRITE(unit,*) ""
    DEALLOCATE(my_i, my_j, my_k)
    !
    ! Section PHI
    !
    IF (.NOT.amber_impropers) THEN
       my_phi = COUNT(MASK=LPH>=0) + COUNT(MASK=LP>=0)
       ALLOCATE(my_i(my_phi), my_j(my_phi), my_k(my_phi), my_l(my_phi))
       ind = 0
       DO I = 1, NPHIH
          IF (LPH(i) >= 0) THEN
             ind = ind + 1
             my_i(ind) = IPH(i)/3+1
             my_j(ind) = JPH(i)/3+1
             my_k(ind) = ABS(KPH(i))/3+1
             my_l(ind) = ABS(LPH(i))/3+1
          END IF
       END DO
       DO I = 1, NPHIA
          IF (LP(i) >= 0) THEN
             ind = ind + 1
             my_i(ind) = IP(i)/3+1
             my_j(ind) = JP(i)/3+1
             my_k(ind) = ABS(KP(i))/3+1
             my_l(ind) = ABS(LP(i))/3+1
          END IF
       END DO
       IF (ind .NE. my_phi) CALL stop_converter("Error in write_psf :: evaluation of proper torsion.")
       WRITE(unit,*) " ",my_phi," !NPHI"
       WRITE(unit,'(8'//int_format//')') (my_i(local),my_j(local),my_k(local),my_l(local), local=1,my_phi)
       WRITE(unit,*) ""
       DEALLOCATE(my_i, my_j, my_k, my_l)
       !
       ! Section IMPHI
       !
       my_phi = COUNT(MASK=LPH<0) + COUNT(MASK=LP<0)
       ALLOCATE(my_i(my_phi), my_j(my_phi), my_k(my_phi), my_l(my_phi))
       ind = 0
       DO I = 1, NPHIH
          IF (LPH(i) < 0) THEN
             ind = ind + 1
             my_i(ind) = IPH(i)/3+1
             my_j(ind) = JPH(i)/3+1
             my_k(ind) = ABS(KPH(i))/3+1
             my_l(ind) = ABS(LPH(i))/3+1
          END IF
       END DO
       DO I = 1, NPHIA
          IF (LP(i) < 0) THEN
             ind = ind + 1
             my_i(ind) = IP(i)/3+1
             my_j(ind) = JP(i)/3+1
             my_k(ind) = ABS(KP(i))/3+1
             my_l(ind) = ABS(LP(i))/3+1
          END IF
       END DO
       IF (ind .NE. my_phi) CALL stop_converter("Error in write_psf :: evaluation of improper torsion.")
       WRITE(unit,*) " ",my_phi," !NIMPHI"
       WRITE(unit,'(8'//int_format//')') (my_i(local),my_j(local),my_k(local),my_l(local), local=1,my_phi)
       WRITE(unit,*) ""
       DEALLOCATE(my_i, my_j, my_k, my_l)
    ELSE
       my_phi = NPHIH + NPHIA
       ALLOCATE(my_i(my_phi), my_j(my_phi), my_k(my_phi), my_l(my_phi))
       ind = 0
       DO I = 1, NPHIH
          ind = ind + 1
          my_i(ind) = IPH(i)/3+1
          my_j(ind) = JPH(i)/3+1
          my_k(ind) = ABS(KPH(i))/3+1
          my_l(ind) = ABS(LPH(i))/3+1
       END DO
       DO I = 1, NPHIA          
          ind = ind + 1
          my_i(ind) = IP(i)/3+1
          my_j(ind) = JP(i)/3+1
          my_k(ind) = ABS(KP(i))/3+1 
          my_l(ind) = ABS(LP(i))/3+1
       END DO
       IF (ind .NE. my_phi) CALL stop_converter("Error in write_psf :: evaluation of proper torsion.")
       WRITE(unit,*) " ",my_phi," !NPHI"
       WRITE(unit,'(8'//int_format//')') (my_i(local),my_j(local),my_k(local),my_l(local), local=1,my_phi)
       WRITE(unit,*) ""
       DEALLOCATE(my_i, my_j, my_k, my_l)
       WRITE(unit,*) " ",0," !NIMPHI"
       WRITE(unit,*) ""
    END IF
    !
    ! Section NDON
    !
    WRITE(unit,*) " 0 !NDON"
    WRITE(unit,*) ""
    !
    ! Section NACC
    !
    WRITE(unit,*) " 0 !NACC"
    WRITE(unit,*) ""
    !
    ! Section NNB
    !
    WRITE(unit,*) " 0 !NNB"
    WRITE(unit,*) ""    
    !
    ! Section NGRP 
    !
    WRITE(unit,*) " 0 !NGRP"
    WRITE(unit,*) ""    
    DEALLOCATE(molnames)
  END SUBROUTINE write_psf



  SUBROUTINE  create_potpar_file(unit)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: unit
    INTEGER :: my_bond, my_angle, my_torsion, my_improper, my_nonbond
    INTEGER :: my_count, my_atom
    INTEGER,       DIMENSION (:), POINTER :: my_i, my_j ,my_k, my_l ! atom indexs
    INTEGER,       DIMENSION (:), POINTER :: my_atom_index          ! index of atoms in residue
    INTEGER,       DIMENSION (:), POINTER :: dihedral_type          ! 0: Proper ; 1:Improper
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_kb, my_req          ! stretching
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_ka, my_theta        ! bending
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_a, my_delta, my_m   ! torsion
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_sigma, my_epsilon, my_charges
    CHARACTER(LEN=4), DIMENSION (:), POINTER :: my_type

    NULLIFY( my_i, my_j, my_k, my_l, my_atom_index, dihedral_type,&
             my_kb, my_ka, my_req, my_theta, my_a, my_delta, my_m,&
             my_sigma, my_epsilon, my_type, my_charges)
    ALLOCATE(my_atom_index(NATOM))
    DO my_atom = 1, NATOM
       my_atom_index(my_atom) = (my_atom-1)*3
    END DO
    !
    ! Header file...
    !
    WRITE(unit,1000)"AMBER FF Converted into CHARMM FF style"
    WRITE(unit,1000)"Generated on :: "//TRIM(Today)//" by :: "//TRIM(ThisIsMe)
    WRITE(unit,1000)TRIM(MyMachine)
    WRITE(unit,1000)"Leap Title :: "//TRIM(TITLE)
    WRITE(unit,1000)"Send all comments related to the FFs conversion to"
    WRITE(unit,1000)"teodoro.laino@gmail.com"
    WRITE(unit,*)""
    !
    ! Bonds
    !
    WRITE(unit,1001)
    my_count = &
         NUMBER_UNIQUE_RES_BONDS(my_atom_index, my_i, my_j, my_kb, my_req)
    IF (verbose) THEN
       DO my_bond = 1, my_count
          WRITE(*,'(2A10,2I5,2F12.6)')ISYMBL(my_i(my_bond)), ISYMBL(my_j(my_bond)),&
                                      my_i(my_bond), my_j(my_bond), my_kb(my_bond),&
                                      my_req(my_bond)
       END DO
    END IF
    my_count = &
         GIVE_BACK_UNIQUE_TYPES(my_i=my_i, my_j=my_j, A=my_kb, B=my_req)
    DO my_bond = 1, my_count
       WRITE(unit,2001)give_back_type(my_i(my_bond)),&
                       give_back_type(my_j(my_bond)),&
                       my_kb(my_bond), my_req(my_bond)
    END DO
    WRITE(unit,*)""
    !
    ! Angles
    !
    WRITE(unit,1002)
    my_count = &
         NUMBER_UNIQUE_RES_ANGLES(my_atom_index, my_i, my_j, my_k, my_ka, my_theta)
    my_count = &
         GIVE_BACK_UNIQUE_TYPES(my_i=my_i, my_j=my_j, my_k=my_k, A=my_ka, B=my_theta)
    DO my_angle = 1, my_count
       WRITE(unit,2002)give_back_type(my_i(my_angle)),&
                       give_back_type(my_j(my_angle)),&
                       give_back_type(my_k(my_angle)),&
                       my_ka(my_angle), CONVERT_ANGLE(radiant=my_theta(my_angle))
    END DO
    WRITE(unit,*)""
    !
    ! Dihedrals
    !
    WRITE(unit,1003)
    my_count = &
         NUMBER_UNIQUE_RES_DIHEDRALS(my_atom_index, my_i, my_j, my_k, my_l, my_a, my_delta,&
                                     my_m, dihedral_type)
    my_count = &
         GIVE_BACK_UNIQUE_TYPES(my_i=my_i, my_j=my_j, my_k=my_k, my_l=my_l,&
                                A=my_a, B=my_m, C=my_delta, D=dihedral_type)
    DO my_torsion = 1, my_count
       IF (dihedral_type(my_torsion) == 0) THEN
          WRITE(unit,2003)give_back_type(my_i(my_torsion)),&
                          give_back_type(my_j(my_torsion)),&
                          give_back_type(my_k(my_torsion)),&
                          give_back_type(my_l(my_torsion)),&
                          my_a(my_torsion),&
                          INT(my_m(my_torsion)),&
                          CONVERT_ANGLE(radiant=my_delta(my_torsion)) 
       END IF
    END DO
    DO my_improper = 1, my_count
       IF (dihedral_type(my_improper) == 1) THEN
          WRITE(unit,2003)give_back_type(my_i(my_improper)),&
                          give_back_type(my_j(my_improper)),&
                          give_back_type(my_k(my_improper)),&
                          give_back_type(my_l(my_improper)),&
                          my_a(my_improper),&
                          INT(my_m(my_improper)),&
                          CONVERT_ANGLE(radiant=my_delta(my_improper))
       END IF
    END DO
    WRITE(unit,*)""
    !
    ! Impropers.. The improper list is empty because AMBER treats impropers
    !             in the same way as torsion (same potential). So all 
    !             impropers have been added to the torsion list...
    !
    WRITE(unit,1004)
    WRITE(unit,*)""
    !
    ! Nonbonded
    !
    WRITE(unit,1005)
    my_count = NUMBER_UNIQUE_LJ(my_atom_index, type_i=my_type, vdw_a=my_sigma, vdw_b=my_epsilon,&
                                charges=my_charges)
    DO my_nonbond = 1, my_count
       WRITE(unit,2004)my_type(my_nonbond),&
                       my_epsilon(my_nonbond),&
                       my_sigma(my_nonbond)
    END DO
    WRITE(unit,*)""   
    !
    ! End of CHARMM format
    !
    WRITE(unit,'(A)')"END"
    !
    ! Print out unique CHARGES...
    !
    WRITE(unit,1006)
    DO my_nonbond = 1, my_count
       WRITE(unit,2005)my_type(my_nonbond),&
                       my_charges(my_nonbond)
    END DO
    WRITE(unit,1007)

    IF (ASSOCIATED(my_i         )) DEALLOCATE(my_i         )
    IF (ASSOCIATED(my_j         )) DEALLOCATE(my_j         )
    IF (ASSOCIATED(my_k         )) DEALLOCATE(my_k         )
    IF (ASSOCIATED(my_l         )) DEALLOCATE(my_l         )
    IF (ASSOCIATED(my_atom_index)) DEALLOCATE(my_atom_index)           
    IF (ASSOCIATED(dihedral_type)) DEALLOCATE(dihedral_type)        
    IF (ASSOCIATED(my_a         )) DEALLOCATE(my_a         )        
    IF (ASSOCIATED(my_m         )) DEALLOCATE(my_m         )        
    IF (ASSOCIATED(my_delta     )) DEALLOCATE(my_delta     )        
    IF (ASSOCIATED(my_epsilon   )) DEALLOCATE(my_epsilon   )
    IF (ASSOCIATED(my_sigma     )) DEALLOCATE(my_sigma     )
    RETURN
1000 FORMAT("*>>>>>>>",T12,A,T73,"<<<<<<<")
1001 FORMAT("BONDS",/,"!",/,"!V(bond) = Kb(b - b0)**2",/,"!",/,"!Kb: kcal/mole/A**2",/, &
            "!b0: A",/,"!",/,"!atom type Kb          b0",/,"!")
1002 FORMAT("ANGLES",/,"!",/,"!V(angle) = Ktheta(Theta - Theta0)**2",/,"!",/,           &
            "!V(Urey-Bradley) = Kub(S - S0)**2",/,"!",/,"!Ktheta: kcal/mole/rad**2",/,  &
            "!Theta0: degrees",/,"!Kub: kcal/mole/A**2 (Urey-Bradley)",/,"!S0: A",/,    &
            "!",/,"!atom types     Ktheta    Theta0   Kub     S0",/,"!")
1003 FORMAT("DIHEDRALS",/,"!",/,"!V(dihedral) = Kchi(1 + cos(n(chi) - delta))",/,       &
            "!",/,"!Kchi: kcal/mole",/,"!n: multiplicity",/,"!delta: degrees",/,        &
            "!",/,"!atom types             Kchi    n   delta",/,"!")
1004 FORMAT("IMPROPER",/,"!",/,"!V(improper) = Kpsi(psi - psi0)**2",/,"!",/,            &
            "!Kpsi: kcal/mole/rad**2",/,"!psi0: degrees",/,                             &
            "!note that the second column of numbers (0) is ignored",/,"!",/,           &
            "!atom types           Kpsi                   psi0",/,"!")
1005 FORMAT("NONBONDED",/,"!",/,                                                        &
            "!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]",/, &
            "!",/,"!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)",/,               &
            "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j",/,"!",/,                       &
            "!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4"&
            ,/,"!")
1006 FORMAT(/,"!",/,"! This Section can be cutted & pasted into the Fist input file..", &
            /"!",/,"CHARGES")
1007 FORMAT("END CHARGES")

2001 FORMAT(A4,1X,A4,1X,2F15.9)                     ! bond
2002 FORMAT(A4,1X,A4,1X,A4,1X,2F15.9)               ! angle
2003 FORMAT(A4,1X,A4,1X,A4,1X,A4,1X,F15.9,I5,F15.9) ! torsion
2004 FORMAT(A4,1X,"    0.000000000",2F15.9)         ! nonbond
2005 FORMAT(A4,1X,F15.9)                            ! charges
  END SUBROUTINE create_potpar_file


  SUBROUTINE write_dlpoly_field (unit)
    IMPLICIT NONE
    !!$ local
    INTEGER :: i, j
    INTEGER :: ires, my_nres
    INTEGER :: nrept=1
    INTEGER :: my_count
    
    INTEGER, INTENT(IN) :: unit ! output unit
    !!$ force field parameters
    INTEGER,       DIMENSION (:), POINTER :: my_i, my_j ,my_k, my_l ! atom indexs
    INTEGER,       DIMENSION (:), POINTER :: my_atom_index          ! index of atoms in residue
    INTEGER,       DIMENSION (:), POINTER :: dihedral_type          ! 0: Proper ; 1:Improper
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_kb, my_req          ! stretching
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_ka, my_theta        ! bending
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_a, my_delta, my_m   ! torsion
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_ki                  ! improper
    REAL(KIND=dp), DIMENSION (:), POINTER :: vdw_a, vdw_b           ! lj
    REAL(KIND=dp), DIMENSION (:), POINTER :: my_scale_el, my_scale_LJ  ! scale factors
    CHARACTER(LEN=4), DIMENSION (:), POINTER :: type_i, type_j         ! lj
    CHARACTER(LEN=4), DIMENSION (:), POINTER :: my_res_name
    
    NULLIFY(my_atom_index)
    NULLIFY(my_i         )
    NULLIFY(my_j         )
    NULLIFY(my_k         )
    NULLIFY(my_l         )
    NULLIFY(dihedral_type)
    NULLIFY(my_kb        )
    NULLIFY(my_req       )
    NULLIFY(my_ka        )
    NULLIFY(my_theta     )
    NULLIFY(my_a         )
    NULLIFY(my_delta     )
    NULLIFY(my_m         )
    NULLIFY(my_ki        )
    NULLIFY(vdw_a        )
    NULLIFY(vdw_b        )
    NULLIFY(my_scale_el  )
    NULLIFY(my_scale_LJ  )
    NULLIFY(type_i       )
    NULLIFY(type_j       )
    NULLIFY(my_res_name  )

    !!$ FILED formatting; integers :: i5, reals :: e12.4
    !!$ header :: a80
    WRITE(unit,'(a)')"FIELD FILE GENERATED FROM AMBER8 PARM FILE"
    
    !!$ units  :: a40
    WRITE(unit,'(a)')"UNITS kcal"
    
    !!$ Type of molecules
    my_nres=NUMBER_UNIQUE_RES(unique_res_name=my_res_name)  ! to be fixed for connected units...
    WRITE(unit,'(a9,2x,i5)')"MOLECULES",my_nres
    
    !!$ write the field for each residue type
    unique_residue : DO ires = 1, my_nres
       
       !!$ molecule name
       WRITE(unit,'(a4)') my_res_name(ires)
       
       !!$ number of molecules of type ires
       my_count= NUMBER_UNIQUE_RES(my_res_name=my_res_name(ires))
       WRITE(unit,'(a8,i5)')"NUMMOLS ", my_count
       
       !!$ list of atoms with mass and charge
       !!$ ifriz and igrp are not specified
       my_count = NUMBER_UNIQUE_RES_ATOMS(my_res_name(ires), my_atom_index)
       WRITE(unit,'(a,2x,i10)')"ATOMS", my_count
       atoms : DO j = 1, SIZE(my_atom_index)
          i = my_atom_index(j)/3+1
          WRITE(unit,'(a8,e12.4,e12.4,i5)')ISYMBL(i),AMASS(i),CHRG(i)/CV_CHRG,nrept
       END DO atoms
       
       !!$ write the bonds
       my_count = &
            NUMBER_UNIQUE_RES_BONDS(my_atom_index, my_i, my_j, my_kb, my_req)
       WRITE(unit,'(a,2x,i10)')"BONDS",my_count
       bonds : DO i=1,my_count
          WRITE(unit,'(a4,i5,i5,e12.4,e12.4)') &
               "harm", my_i(i), my_j(i), 2.0*my_kb(i), my_req(i)
       END DO bonds
       
       !!$ write the angles
       my_count = &
            NUMBER_UNIQUE_RES_ANGLES(my_atom_index, my_i, my_j, my_k, my_ka, my_theta)
       WRITE(unit,'(a,2x,i10)')"ANGLES",my_count
       angles : DO i=1,my_count
          WRITE(unit,'(a4,i5,i5,i5,e12.4,e12.4)') &
               "harm", my_i(i), my_j(i), my_k(i), 2.0*my_ka(i), CONVERT_ANGLE(radiant=my_theta(i))
       END DO angles
       
       !!$ write dihedrals
       my_count = &
            NUMBER_UNIQUE_RES_DIHEDRALS(my_atom_index, my_i, my_j, my_k, my_l, my_a, my_delta,&
                                        my_m, dihedral_type, my_scale_el, my_scale_LJ)
       WRITE(unit,'(a,2x,i10)')"DIHEDRALS",my_count
       dihedrals : DO i=1,my_count
          WRITE(unit,'(a4,4i5,5e12.4)') &
               "cos ", my_i(i), my_j(i), my_k(i), my_l(i), my_a(i),&
               CONVERT_ANGLE(radiant=my_delta(i)), my_m(i),my_scale_el(i), my_scale_LJ(i)
       END DO dihedrals
       
       WRITE(unit,'(a6)')"FINISH"
    END DO unique_residue
    
    !!$ write the van der Waals parameters
    my_count = &
         NUMBER_UNIQUE_LJ(my_atom_index, type_i, type_j, vdw_a, vdw_b)
    WRITE(unit,'(a,2x,i10)')"VDW",my_count  
    vdw : DO i=1,my_count
       WRITE(unit,'(a8,a8,a4,1x,e12.4,e12.4)')type_i(i), type_j(i),"12-6",vdw_a(i), vdw_b(i)
    END DO vdw
    
    WRITE(unit,'(a5)')"CLOSE"
    !
    ! Deallocate all structures...
    !
     IF (ASSOCIATED(my_i         )) DEALLOCATE(my_i         )
     IF (ASSOCIATED(my_j         )) DEALLOCATE(my_j         )
     IF (ASSOCIATED(my_k         )) DEALLOCATE(my_k         )
     IF (ASSOCIATED(my_l         )) DEALLOCATE(my_l         )
     IF (ASSOCIATED(my_atom_index)) DEALLOCATE(my_atom_index)           
     IF (ASSOCIATED(dihedral_type)) DEALLOCATE(dihedral_type)          
     IF (ASSOCIATED(my_kb        )) DEALLOCATE(my_kb        )
     IF (ASSOCIATED(my_req       )) DEALLOCATE(my_req       )   
     IF (ASSOCIATED(my_ka        )) DEALLOCATE(my_ka        )
     IF (ASSOCIATED(my_theta     )) DEALLOCATE(my_theta     )   
     IF (ASSOCIATED(my_a         )) DEALLOCATE(my_a         )
     IF (ASSOCIATED(my_delta     )) DEALLOCATE(my_delta     )
     IF (ASSOCIATED(my_m         )) DEALLOCATE(my_m         )
     IF (ASSOCIATED(my_ki        )) DEALLOCATE(my_ki        )          
     IF (ASSOCIATED(vdw_a        )) DEALLOCATE(vdw_a        )
     IF (ASSOCIATED(vdw_b        )) DEALLOCATE(vdw_b        )    
     IF (ASSOCIATED(my_scale_el  )) DEALLOCATE(my_scale_el  )
     IF (ASSOCIATED(my_scale_LJ  )) DEALLOCATE(my_scale_LJ  )
     IF (ASSOCIATED(type_i       )) DEALLOCATE(type_i       )
     IF (ASSOCIATED(type_j       )) DEALLOCATE(type_j       ) 
     IF (ASSOCIATED(my_res_name  )) DEALLOCATE(my_res_name  )
     
    RETURN
  END SUBROUTINE write_dlpoly_field
  
  INTEGER FUNCTION NUMBER_UNIQUE_RES(unique_res_name, my_res_name) RESULT(unique_res)
    IMPLICIT NONE
    CHARACTER (LEN=4), DIMENSION(:), POINTER, OPTIONAL :: unique_res_name
    CHARACTER (LEN=4), INTENT(IN), OPTIONAL  :: my_res_name
    INTEGER :: I, j
    CHARACTER (LEN=4), DIMENSION(:), POINTER :: local_res 

    IF (PRESENT(unique_res_name)) THEN
       !
       ! Determine the number of unique residues and give back their names...
       !
       ALLOCATE(local_res(SIZE(LABRES)))
       local_res  = LABRES
       unique_res = 0
       DO i = 1, NRES
          IF (TRIM(local_res(i)) == "NULL") CYCLE
          unique_res = unique_res + 1
          DO j = i+1, NRES
             IF (TRIM(local_res(j)) == TRIM(local_res(i))) local_res(j) = "NULL"
          END DO
       END DO
       IF(verbose) WRITE(*,'(A,I5)')"number_unique_res :: unique_res ::",unique_res 
       ALLOCATE(unique_res_name(unique_res))
       unique_res = 0
       DO i = 1, NRES
          IF (TRIM(local_res(i)) /= "NULL") THEN 
             unique_res = unique_res + 1
             unique_res_name(unique_res) = TRIM(local_res(i)) 
          END IF
       END DO
       IF (verbose) WRITE(*,'(A)')"number_unique_res :: unique_res_name :: ",&
            unique_res_name
       IF (unique_res /= SIZE(unique_res_name))&
            CALL stop_converter("Array of different size :: NUMBER_UNIQUE_RES")
       DEALLOCATE(local_res)
    ELSEIF (PRESENT(my_res_name)) THEN
       !
       ! Gives back the number of residues of my_res_name 
       !
       unique_res = 0
       DO i = 1, NRES
          IF (TRIM(LABRES(i)) == TRIM(my_res_name)) THEN 
             unique_res = unique_res + 1
          END IF
       END DO
       IF (verbose) WRITE(*,'(A,I7)')"number_unique_res :: my_resname ::"&
            //my_res_name//" :: ",unique_res
    END IF

  END FUNCTION NUMBER_UNIQUE_RES


  INTEGER FUNCTION NUMBER_UNIQUE_RES_ATOMS(my_res_name, my_atom_index) RESULT(number_of_atom)
    IMPLICIT NONE
    CHARACTER (LEN=4), INTENT(IN)            :: my_res_name
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index 
    INTEGER :: I, J, istart, iend

    IF (ASSOCIATED(my_atom_index)) DEALLOCATE(my_atom_index)
    DO I = 1, NRES
       IF (TRIM(LABRES(I)) == TRIM(my_res_name)) EXIT
    END DO
    IF (verbose) WRITE(*,'(A)')"RESNAME in number_unique_res_atoms :: "//my_res_name
    istart = IPRES(i)
    IF ( i.EQ.NRES ) THEN 
       iend = NATOM
    ELSE
       iend = ipres(i+1)-1
    END IF
    ALLOCATE(my_atom_index(iend-istart+1))
    number_of_atom = 0
    DO J = istart, iend
       number_of_atom = number_of_atom + 1
       my_atom_index(number_of_atom) = (J-1)*3
    END DO    
    IF (verbose) WRITE(*,'(A,I7)')"NUMBER OF ATOMS in number_unique_res_atoms :: ",number_of_atom
  END FUNCTION NUMBER_UNIQUE_RES_ATOMS


  INTEGER FUNCTION NUMBER_UNIQUE_RES_BONDS(my_atom_index, my_i, my_j, my_kb, my_req)&
       RESULT(number_of_bonds)
    IMPLICIT NONE
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, my_i, my_j
    REAL(KIND=dp),    DIMENSION (:), POINTER :: my_kb, my_req

    IF (ASSOCIATED(my_i))   DEALLOCATE(my_i)
    IF (ASSOCIATED(my_j))   DEALLOCATE(my_j)
    IF (ASSOCIATED(my_kb))  DEALLOCATE(my_kb)
    IF (ASSOCIATED(my_req)) DEALLOCATE(my_req)
    number_of_bonds = 0
    CALL loop_on_bonds(number_of_bonds, NBONH, IBH, JBH, my_atom_index)
    CALL loop_on_bonds(number_of_bonds, MBONA, IB,  JB,  my_atom_index)
    ALLOCATE(my_i(number_of_bonds),&
             my_j(number_of_bonds),&
            my_kb(number_of_bonds),&
           my_req(number_of_bonds) )
    number_of_bonds = 0
    CALL loop_on_bonds(number_of_bonds, NBONH, IBH, JBH, my_atom_index, ICBH, my_i, my_j, my_kb, my_req)    
    CALL loop_on_bonds(number_of_bonds, MBONA, IB,  JB,  my_atom_index, ICB,  my_i, my_j, my_kb, my_req)
    IF (verbose) WRITE(*,*)"Total number of bonds :: ",number_of_bonds
    IF (SIZE(my_i) /= number_of_bonds)&
         CALL stop_converter("Array of different size :: NUMBER_UNIQUE_RES_BONDS")
  END FUNCTION NUMBER_UNIQUE_RES_BONDS


  SUBROUTINE loop_on_bonds(number_of_bonds, ndim, ib_l, jb_l, my_atom_index,&
       icb_l, my_i, my_j, my_kb, my_req)
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: number_of_bonds
    INTEGER, INTENT(in)    :: ndim
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, ib_l, jb_l
    INTEGER,          DIMENSION (:), POINTER, OPTIONAL :: my_i, my_j, icb_l
    REAL(KIND=dp),    DIMENSION (:), POINTER, OPTIONAL :: my_kb, my_req
    INTEGER :: i

    bonds: DO i = 1, ndim
       IF (ANY(my_atom_index==ib_l(i)).AND.ANY(my_atom_index==ib_l(i))) THEN
          number_of_bonds = number_of_bonds + 1
          IF (PRESENT(my_req)) THEN
             my_i(number_of_bonds)   = (ib_l(i)-my_atom_index(1))/3+1
             my_j(number_of_bonds)   = (jb_l(i)-my_atom_index(1))/3+1
             my_kb(number_of_bonds)  = RK(icb_l(i))
             my_req(number_of_bonds) = REQ(icb_l(i))
          END IF
       END IF
    END DO bonds
    IF (verbose) WRITE(*,*)"Entering loop_on_bonds :: number_of_bonds ::",number_of_bonds
  END SUBROUTINE loop_on_bonds

  INTEGER FUNCTION NUMBER_UNIQUE_RES_ANGLES(my_atom_index, my_i, my_j, my_k, my_ka, my_theta)&
       RESULT(number_of_angles)
    IMPLICIT NONE
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, my_i, my_j, my_k
    REAL(KIND=dp),    DIMENSION (:), POINTER :: my_ka, my_theta

    IF (ASSOCIATED(my_i))     DEALLOCATE(my_i)
    IF (ASSOCIATED(my_j))     DEALLOCATE(my_j)
    IF (ASSOCIATED(my_k))     DEALLOCATE(my_k)    
    IF (ASSOCIATED(my_ka))    DEALLOCATE(my_ka)
    IF (ASSOCIATED(my_theta)) DEALLOCATE(my_theta)

    number_of_angles = 0
    CALL loop_on_angles(number_of_angles, NTHETH, ITH, JTH, KTH, my_atom_index)
    CALL loop_on_angles(number_of_angles, MTHETA, IT,  JT,  KT,  my_atom_index)
    ALLOCATE(my_i(number_of_angles),&
             my_j(number_of_angles),&
             my_k(number_of_angles),&
            my_ka(number_of_angles),&
         my_theta(number_of_angles) )
    number_of_angles = 0
    CALL loop_on_angles(number_of_angles, NTHETH, ITH, JTH, KTH, my_atom_index,&
         ICTH, my_i, my_j, my_k, my_ka, my_theta)
    CALL loop_on_angles(number_of_angles, MTHETA, IT,   JT,  KT, my_atom_index,&
         ICT , my_i, my_j, my_k, my_ka, my_theta)
    IF (verbose) WRITE(*,*)"Total number of angles :: ",number_of_angles
    IF (SIZE(my_i) /= number_of_angles)&
         CALL stop_converter("Array of different size :: NUMBER_UNIQUE_RES_ANGLES")
  END FUNCTION NUMBER_UNIQUE_RES_ANGLES


  SUBROUTINE loop_on_angles(number_of_angles, ndim, it_l, jt_l, kt_l, my_atom_index,&
       ict_l, my_i, my_j, my_k, my_ka, my_theta)
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: number_of_angles
    INTEGER, INTENT(in)    :: ndim
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, it_l, jt_l, kt_l
    INTEGER,          DIMENSION (:), POINTER, OPTIONAL :: my_i, my_j, my_k, ict_l
    REAL(KIND=dp),    DIMENSION (:), POINTER, OPTIONAL :: my_ka, my_theta
    INTEGER :: i

    angles: DO i = 1, ndim
       IF  (ANY(my_atom_index==it_l(i)).AND.ANY(my_atom_index==jt_l(i)).AND.&
            ANY(my_atom_index==kt_l(i))) THEN
          number_of_angles = number_of_angles + 1
          IF (PRESENT(my_theta)) THEN
             my_i(number_of_angles)     = (it_l(i)-my_atom_index(1))/3+1
             my_j(number_of_angles)     = (jt_l(i)-my_atom_index(1))/3+1
             my_k(number_of_angles)     = (kt_l(i)-my_atom_index(1))/3+1
             my_ka(number_of_angles)    = TK(ict_l(i))
             my_theta(number_of_angles) = TEQ(ict_l(i))
          END IF
       END IF
    END DO angles
    IF (verbose) WRITE(*,*)"Entering loop_on_angles :: number_of_angles ::",number_of_angles
  END SUBROUTINE loop_on_angles
  

  INTEGER FUNCTION NUMBER_UNIQUE_RES_DIHEDRALS(my_atom_index, my_i, my_j, my_k, my_l, my_a, my_delta,&
       my_m, dihedral_type, my_scale_el, my_scale_LJ) RESULT(number_of_dihedrals)
    IMPLICIT NONE
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, my_i, my_j, my_k, my_l
    INTEGER,          DIMENSION (:), POINTER, OPTIONAL :: dihedral_type
    REAL(KIND=dp),    DIMENSION (:), POINTER :: my_a, my_delta, my_m
    REAL(KIND=dp),    DIMENSION (:), POINTER, OPTIONAL :: my_scale_el, my_scale_LJ

    IF (ASSOCIATED(my_i))        DEALLOCATE(my_i)
    IF (ASSOCIATED(my_j))        DEALLOCATE(my_j)
    IF (ASSOCIATED(my_k))        DEALLOCATE(my_k)
    IF (ASSOCIATED(my_l))        DEALLOCATE(my_l)
    IF (ASSOCIATED(my_a))        DEALLOCATE(my_a)
    IF (ASSOCIATED(my_delta))    DEALLOCATE(my_delta)
    IF (ASSOCIATED(my_m))        DEALLOCATE(my_m)
    IF (PRESENT(my_scale_el)) THEN
       IF (ASSOCIATED(my_scale_el)) DEALLOCATE(my_scale_el)
    END IF
    IF (PRESENT(my_scale_LJ)) THEN 
       IF (ASSOCIATED(my_scale_LJ)) DEALLOCATE(my_scale_LJ)
    END IF
    IF (PRESENT(dihedral_type)) THEN 
       IF (ASSOCIATED(dihedral_type)) DEALLOCATE(dihedral_type)
    END IF
    number_of_dihedrals = 0
    CALL loop_on_dihedrals(number_of_dihedrals, NPHIH, IPH, JPH, KPH, LPH, my_atom_index)
    CALL loop_on_dihedrals(number_of_dihedrals, MPHIA, IP,  JP,  KP,   LP, my_atom_index)
    ALLOCATE(my_i(number_of_dihedrals),&
             my_j(number_of_dihedrals),&
             my_k(number_of_dihedrals),&
             my_l(number_of_dihedrals),&
             my_a(number_of_dihedrals),&
         my_delta(number_of_dihedrals),&
             my_m(number_of_dihedrals) )
    IF (PRESENT(my_scale_el)) THEN
       ALLOCATE(my_scale_el(number_of_dihedrals))
    END IF
    IF (PRESENT(my_scale_LJ)) THEN 
       ALLOCATE(my_scale_LJ(number_of_dihedrals))
    END IF
    IF (PRESENT(dihedral_type)) THEN 
       ALLOCATE(dihedral_type(number_of_dihedrals))
    END IF
    number_of_dihedrals = 0
    IF (PRESENT(my_scale_el).AND.PRESENT(my_scale_LJ)) THEN
       CALL loop_on_dihedrals(number_of_dihedrals, NPHIH, IPH, JPH, KPH, LPH, my_atom_index,&
            ICPH, my_i, my_j, my_k, my_l, my_a, my_delta, my_m, dihedral_type, my_scale_el, &
            my_scale_LJ )
       CALL loop_on_dihedrals(number_of_dihedrals, MPHIA, IP,  JP,  KP,   LP, my_atom_index,&
            ICP,  my_i, my_j, my_k, my_l, my_a, my_delta, my_m, dihedral_type, my_scale_el, &
            my_scale_LJ )
    ELSE
       CALL loop_on_dihedrals(number_of_dihedrals, NPHIH, IPH, JPH, KPH, LPH, my_atom_index,&
            ICPH, my_i, my_j, my_k, my_l, my_a, my_delta, my_m, dihedral_type)
       CALL loop_on_dihedrals(number_of_dihedrals, MPHIA, IP,  JP,  KP,   LP, my_atom_index,&
            ICP,  my_i, my_j, my_k, my_l, my_a, my_delta, my_m, dihedral_type)
    END IF
    IF (verbose) WRITE(*,*)"Total number of dihedrals :: ",number_of_dihedrals
    IF (SIZE(my_i) /= number_of_dihedrals)&
         CALL stop_converter("Array of different size :: NUMBER_UNIQUE_RES_DIHEDRALS")
  END FUNCTION NUMBER_UNIQUE_RES_DIHEDRALS


  SUBROUTINE loop_on_dihedrals(number_of_dihedrals, ndim, it_l, jt_l, kt_l, lt_l, my_atom_index,&
       ict_l, my_i, my_j, my_k, my_l, my_a, my_delta, my_m,  dihedral_type, my_scale_el, my_scale_LJ )
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: number_of_dihedrals
    INTEGER, INTENT(in)    :: ndim
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, it_l, jt_l, kt_l, lt_l
    INTEGER,          DIMENSION (:), POINTER, OPTIONAL :: my_i, my_j, my_k, ict_l, my_l, dihedral_type
    REAL(KIND=dp),    DIMENSION (:), POINTER, OPTIONAL :: my_a, my_delta, my_m
    REAL(KIND=dp),    DIMENSION (:), POINTER, OPTIONAL :: my_scale_el, my_scale_LJ
    INTEGER :: i

    dihedrals: DO i = 1, ndim
       IF  (ANY(my_atom_index==it_l(i)).AND.ANY(my_atom_index==jt_l(i)).AND.&
            ANY(my_atom_index==kt_l(i)).AND.ANY(my_atom_index==lt_l(i))) THEN
          number_of_dihedrals = number_of_dihedrals + 1
          IF (PRESENT(my_delta)) THEN
             my_i(number_of_dihedrals)          = (it_l(i)-my_atom_index(1))/3+1
             my_j(number_of_dihedrals)          = (jt_l(i)-my_atom_index(1))/3+1
             my_k(number_of_dihedrals)          = (ABS(kt_l(i))-my_atom_index(1))/3+1
             my_l(number_of_dihedrals)          = (ABS(lt_l(i))-my_atom_index(1))/3+1
             my_a(number_of_dihedrals)          = PK(ict_l(i))
             my_m(number_of_dihedrals)          = PN(ict_l(i))
             my_delta(number_of_dihedrals)      = PHASE(ict_l(i))
             IF (PRESENT(my_scale_el).AND.PRESENT(my_scale_LJ)) THEN
                my_scale_el(number_of_dihedrals)   = 1.0_dp/1.2_dp 
                my_scale_LJ(number_of_dihedrals)   = 0.5_dp
                IF (kt_l(i) < 0 .OR. lt_l(i) < 0) THEN
                   my_scale_el(number_of_dihedrals)   = 0.0_dp
                   my_scale_LJ(number_of_dihedrals)   = 0.0_dp
                END IF
             END IF
             IF (PRESENT(dihedral_type)) THEN
                dihedral_type(number_of_dihedrals) = 0
                IF (lt_l(i) < 0) dihedral_type(number_of_dihedrals) = 1
             END IF
          END IF
       END IF
    END DO dihedrals
    IF (verbose) WRITE(*,*)"Entering loop_on_dihedrals :: number_of_dihedrals ::",number_of_dihedrals
  END SUBROUTINE loop_on_dihedrals

  
  INTEGER FUNCTION  NUMBER_UNIQUE_LJ(my_atom_index, type_i, type_j, vdw_a, vdw_b, charges)&
       RESULT(number_of_LJ)
    IMPLICIT NONE
    INTEGER,          DIMENSION (:), POINTER :: my_atom_index, local_index
    CHARACTER (LEN=4),DIMENSION (:), POINTER :: type_i, local_type
    CHARACTER (LEN=4),DIMENSION (:), POINTER, OPTIONAL :: type_j
    REAL(KIND=dp),    DIMENSION (:), POINTER :: vdw_a, vdw_b
    REAL(KIND=dp),    DIMENSION (:), POINTER, OPTIONAL :: charges
    INTEGER :: i, j, my_lj_index, num_local_types, i_iac, j_iac
    REAL(KINd=dp) :: F12, F6, my_rij, my_rij6, my_epsij
    LOGICAL :: Isj, Isch

    Isj  = PRESENT(type_j)
    Isch = PRESENT(charges)
    IF (ASSOCIATED(type_i)) DEALLOCATE(type_i)
    IF (Isj) THEN 
       IF (ASSOCIATED(type_j)) DEALLOCATE(type_j)
    END IF
    IF (Isch) THEN
       IF (ASSOCIATED(charges)) DEALLOCATE(charges)
    END IF
    IF (ASSOCIATED(vdw_a )) DEALLOCATE(vdw_a )
    IF (ASSOCIATED(vdw_b )) DEALLOCATE(vdw_b )

    ALLOCATE(local_index(SIZE(my_atom_index)), local_type(SIZE(my_atom_index)))
    local_index = my_atom_index/3+1
    DO I = 1, SIZE(my_atom_index)
       j = local_index(i)
       local_type(i) = ISYMBL(j)
    END DO
    !
    DO i = 1, SIZE(local_type)
       IF (TRIM(local_type(i)) == "NULL") cycle
       DO j = i+1, SIZE(local_type)
          IF (TRIM(local_type(j)) == TRIM(local_type(i))) local_type(j)="NULL"
       END DO
    END DO
    num_local_types = COUNT(MASK=local_type.NE."NULL")
    IF (verbose) WRITE(*,*)"Number_unique_LJ :: num_local_types ::",num_local_types
    IF (Isj) THEN
       ALLOCATE(type_i(num_local_types*(num_local_types+1)/2),&
                type_j(num_local_types*(num_local_types+1)/2),&
                vdw_a (num_local_types*(num_local_types+1)/2),&
                vdw_b (num_local_types*(num_local_types+1)/2))
    ELSE
       ALLOCATE(type_i(num_local_types),&
                vdw_a (num_local_types),&
                vdw_b (num_local_types))       
    END IF
    number_of_LJ= 0
    DO i = 1, SIZE(local_type)
       IF (TRIM(local_type(i)) == "NULL") CYCLE
       IF (Isj) THEN
          DO j = i, SIZE(local_type)
             IF (TRIM(local_type(j)) == "NULL") CYCLE
             number_of_LJ= number_of_LJ+ 1
             type_i(number_of_LJ) = local_type(i)
             type_j(number_of_LJ) = local_type(j)
             i_iac = local_index(i)
             j_iac = local_index(j)
             my_lj_index = ICO ( NTYPES * (IAC(i_iac)-1) + IAC(j_iac))
             IF (my_lj_index .GT. 0) THEN
                vdw_a(number_of_LJ) = CN1(my_lj_index)
                vdw_b(number_of_LJ) = CN2(my_lj_index)
             ELSE
                CALL stop_converter("NUMBER_UNIQUE_LJ :: 10-12 not yet implemented!")
                vdw_a(number_of_LJ) = ASOL(ABS(my_lj_index))
                vdw_b(number_of_LJ) = BSOL(ABS(my_lj_index))            
             END IF
          END DO
       ELSE
          number_of_LJ= number_of_LJ+ 1
          type_i(number_of_LJ) = local_type(i)
          i_iac = local_index(i)
          my_lj_index = ICO ( NTYPES * (IAC(i_iac)-1) + IAC(i_iac))
          IF (my_lj_index .GT. 0) THEN
             F12 =   CN1(my_lj_index);  F6  =   CN2(my_lj_index)
             IF ( F6 == 0.0_dp ) THEN
                IF ( F12 /= 0.0_dp ) &
                     CALL stop_converter("NUMBER_UNIQUE_LJ :: Anomalies in LJ (0.0) coefficients!")
                my_rij   = 0.0_dp
                my_epsij = 0.0_dp
             ELSE
                my_rij6   = (2.0_dp*F12/F6)
                my_rij    = my_rij6**(1.0_dp/6.0_dp)
                my_epsij  = F6/(2.0_dp*my_rij6) 
             END IF
             vdw_a(number_of_LJ) = my_rij/2.0_dp
             vdw_b(number_of_LJ) = my_epsij                
          ELSE
             CALL stop_converter("NUMBER_UNIQUE_LJ :: 10-12 not yet implemented!")
             vdw_a(number_of_LJ) = ASOL(ABS(my_lj_index))
             vdw_b(number_of_LJ) = BSOL(ABS(my_lj_index))            
          END IF
       END IF
    END DO
    IF (Isch) THEN
       IF (Isj) THEN
          CALL stop_converter("NUMBER_UNIQUE_LJ :: Charges not implemented when handling LJ couples!!")
       END IF
       ALLOCATE(charges(number_of_LJ))
       DO i = 1, number_of_LJ
          DO j = 1, NATOM
             IF (TRIM(type_i(i)) == TRIM(ISYMBL(j))) EXIT
          END DO
          IF (j == NATOM+1) CALL stop_converter("NUMBER_UNIQUE_LJ :: Something wrong assigning charges..")
          charges(i) = CHRG(j)/CV_CHRG
       END DO
    END IF
    DEALLOCATE (local_index, local_type)
  END FUNCTION NUMBER_UNIQUE_LJ

  INTEGER FUNCTION GIVE_BACK_UNIQUE_TYPES(my_i, my_j, my_k, my_l, A, B, C, D)   &
       RESULT(unique_dim)
    IMPLICIT NONE
    INTEGER, POINTER, DIMENSION(:) :: my_i, my_j
    INTEGER, POINTER, DIMENSION(:), OPTIONAL ::  my_k, my_l, D
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL :: A, B, C
    INTEGER :: i0, i1, j0, j1, k0, k1, l0, l1, m, n
    INTEGER, POINTER, DIMENSION(:) :: loc_i, loc_j, loc_k, loc_l, loc_D
    REAL(KIND=dp), DIMENSION(:), POINTER :: loc_A, loc_B, loc_C
    LOGICAL :: IsA, IsB, IsC, IsD, Is_k, Is_l, Li, Lj, Lk, Ll, jump

    IsA = PRESENT(A)
    IsB = PRESENT(B)
    IsC = PRESENT(C)
    IsD = PRESENT(D)
    Is_k = PRESENT(my_k)
    Is_l = PRESENT(my_l)
    IF (verbose)          WRITE(*,'("Entering give_back_unique_types (my_i)::",T40,10I8)')my_i
    IF (verbose)          WRITE(*,'("Entering give_back_unique_types (my_j)::",T40,10I8)')my_j
    IF (verbose.AND.Is_k) WRITE(*,'("Entering give_back_unique_types (my_k)::",T40,10I8)')my_k
    IF (verbose.AND.Is_l) WRITE(*,'("Entering give_back_unique_types (my_l)::",T40,10I8)')my_l
    IF (verbose.AND.IsA)  WRITE(*,'("Entering vector A ::",T40,10F8.2)')A
    IF (verbose.AND.IsB)  WRITE(*,'("Entering vector B ::",T40,10F8.2)')B
    IF (verbose.AND.IsC)  WRITE(*,'("Entering vector C ::",T40,10F8.2)')C
    IF (verbose.AND.IsD)  WRITE(*,'("Entering vector D ::",T40,10I8)'  )D
    IF (verbose)          WRITE(*,*)""
    ALLOCATE(loc_i(SIZE(my_i)), loc_j(SIZE(my_j)))
    IF (Is_k) ALLOCATE(loc_k(SIZE(my_k)))
    IF (Is_l) ALLOCATE(loc_l(SIZE(my_l)))
    IF (IsA)  ALLOCATE(loc_A(SIZE(A)))
    IF (IsB)  ALLOCATE(loc_B(SIZE(B)))
    IF (IsC)  ALLOCATE(loc_C(SIZE(C)))
    IF (IsD)  ALLOCATE(loc_D(SIZE(D)))
    loc_i = (my_i-1)*3
    loc_j = (my_j-1)*3
    IF (Is_k) loc_k = (my_k-1)*3
    IF (Is_l) loc_l = (my_l-1)*3
    IF (IsA)  loc_A = A
    IF (IsB)  loc_B = B
    IF (IsC)  loc_C = C
    IF (IsD)  loc_D = D
    k0 = UNDEF
    l0 = UNDEF
    k1 = UNDEF
    l1 = UNDEF
    DO m = 1, SIZE(loc_i)
       i0 = loc_i(m)/3+1
       j0 = loc_j(m)/3+1
       IF (Is_k) k0 = loc_k(m)/3+1
       IF (Is_l) l0 = loc_l(m)/3+1
       IF ((i0 == UNDEF).AND.(j0 == UNDEF)) CYCLE
       IF (Is_k.AND.(k0 == UNDEF)) CYCLE
       IF (Is_l.AND.(l0 == UNDEF)) CYCLE
       DO n = m+1, SIZE(loc_i)             
          i1 = loc_i(n)/3+1
          j1 = loc_j(n)/3+1
          IF (Is_k) k1 = loc_k(n)/3+1
          IF (Is_l) l1 = loc_l(n)/3+1
          IF ((i1 == UNDEF).AND.(j1 == UNDEF)) CYCLE
          IF (Is_k.AND.(k1 == UNDEF)) CYCLE
          IF (Is_l.AND.(l1 == UNDEF)) CYCLE          
          Li = (TRIM(ISYMBL(i0)) == TRIM(ISYMBL(i1)))
          Lj = (TRIM(ISYMBL(j0)) == TRIM(ISYMBL(j1)))
          Lk = .TRUE.
          Ll = .TRUE.          
          IF (Is_k) Lk = (TRIM(ISYMBL(k0)) == TRIM(ISYMBL(k1)))
          IF (Is_l) Ll = (TRIM(ISYMBL(l0)) == TRIM(ISYMBL(l1)))
          !
          ! Invert the selection and check...
          !
          IF (.NOT.Is_k.AND.(.NOT.Is_l)) THEN
             IF ((Li.AND.Lj).OR.((TRIM(ISYMBL(i0)) == TRIM(ISYMBL(j1))).AND.&
                                 (TRIM(ISYMBL(j0)) == TRIM(ISYMBL(i1))))) THEN
                Lj = .TRUE.
                Li = .TRUE.
             END IF
          ELSE IF (.NOT.Is_l) THEN
             IF ((Li.AND.Lk).OR.((TRIM(ISYMBL(i0)) == TRIM(ISYMBL(k1))).AND.&
                                 (TRIM(ISYMBL(k0)) == TRIM(ISYMBL(i1))))) THEN
                 Li = .TRUE.
                 Lk = .TRUE.
              END IF
          ELSE 
             IF ((Li.AND.Lj.AND.Lk.AND.Ll).OR.((TRIM(ISYMBL(i0)) == TRIM(ISYMBL(l1))).AND.&
                                               (TRIM(ISYMBL(j0)) == TRIM(ISYMBL(k1))).AND.&
                                               (TRIM(ISYMBL(k0)) == TRIM(ISYMBL(j1))).AND.&
                                               (TRIM(ISYMBL(l0)) == TRIM(ISYMBL(i1))))) THEN
                 Li = .TRUE.
                 Lj = .TRUE.
                 Lk = .TRUE.
                 Ll = .TRUE.
              END IF
          END IF
          IF (Li .AND. Lj .AND. Lk .AND. Ll) THEN
             jump = .TRUE.
             IF (IsA) jump = jump .AND. (loc_A(m) == loc_A(n))
             IF (IsB) jump = jump .AND. (loc_B(m) == loc_B(n))
             IF (IsC) jump = jump .AND. (loc_C(m) == loc_C(n))
             IF (IsD) jump = jump .AND. (loc_D(m) == loc_D(n))
             IF (jump) THEN
                loc_i(n) = (UNDEF-1)*3
                loc_j(n) = (UNDEF-1)*3
                IF (Is_k) loc_k(n) = (UNDEF-1)*3
                IF (Is_l) loc_l(n) = (UNDEF-1)*3
                IF (IsA) loc_A(n)  = UNDEF
                IF (IsB) loc_B(n)  = UNDEF
                IF (IsC) loc_C(n)  = UNDEF
                IF (IsD) loc_D(n)  = UNDEF
             END IF
          END IF
       END DO
    END DO
    unique_dim = COUNT(MASK=loc_i /= (UNDEF-1)*3)
    IF (verbose) WRITE(*,*)"After first loop :: unique_dim == ", unique_dim
    DEALLOCATE(my_i, my_j)
    IF (Is_k) DEALLOCATE(my_k)
    IF (Is_l) DEALLOCATE(my_l)
    IF (IsA)  DEALLOCATE(A)
    IF (IsB)  DEALLOCATE(B)
    IF (IsC)  DEALLOCATE(C)
    IF (IsD)  DEALLOCATE(D)
    !
    ! Assign values and go back...
    !
    ALLOCATE(my_i(unique_dim), my_j(unique_dim))
    IF (Is_k) ALLOCATE(my_k(unique_dim))
    IF (Is_l) ALLOCATE(my_l(unique_dim))
    IF (IsA)  ALLOCATE(A(unique_dim))
    IF (IsB)  ALLOCATE(B(unique_dim))
    IF (IsC)  ALLOCATE(C(unique_dim))
    IF (IsD)  ALLOCATE(D(unique_dim))
    unique_dim = 0
    DO m = 1, SIZE(loc_i)
       i0 = loc_i(m)/3+1
       j0 = loc_j(m)/3+1
       IF (Is_k) k0 = loc_k(m)/3+1
       IF (Is_l) l0 = loc_l(m)/3+1
       IF ((i0 == UNDEF).AND.(j0 == UNDEF)) CYCLE
       IF (Is_k.AND.(k0 == UNDEF)) CYCLE
       IF (Is_l.AND.(l0 == UNDEF)) CYCLE
       unique_dim = unique_dim + 1
       my_i(unique_dim) = loc_i(m)
       my_j(unique_dim) = loc_j(m)
       IF (Is_k) my_k(unique_dim) = loc_k(m)
       IF (Is_l) my_l(unique_dim) = loc_l(m)
       IF (IsA) A(unique_dim) = loc_A(m)
       IF (IsB) B(unique_dim) = loc_B(m)
       IF (IsC) C(unique_dim) = loc_C(m)
       IF (IsD) D(unique_dim) = loc_D(m)
    END DO 
    IF (verbose)          WRITE(*,'("Exiting give_back_unique_types (my_i)::",T40,10I8)')my_i
    IF (verbose)          WRITE(*,'("Exiting give_back_unique_types (my_j)::",T40,10I8)')my_j
    IF (verbose.AND.Is_k) WRITE(*,'("Exiting give_back_unique_types (my_k)::",T40,10I8)')my_k
    IF (verbose.AND.Is_l) WRITE(*,'("Exiting give_back_unique_types (my_l)::",T40,10I8)')my_l
    IF (verbose.AND.IsA)  WRITE(*,'("Exiting vector A ::",T40,10F8.2)')A
    IF (verbose.AND.IsB)  WRITE(*,'("Exiting vector B ::",T40,10F8.2)')B
    IF (verbose.AND.IsC)  WRITE(*,'("Exiting vector C ::",T40,10F8.2)')C
    IF (verbose.AND.IsD)  WRITE(*,'("Exiting vector D ::",T40,10I8)'  )D
    IF (verbose)          WRITE(*,*)""
    DEALLOCATE(loc_i, loc_j)
    IF (Is_k) DEALLOCATE(loc_k)
    IF (Is_l) DEALLOCATE(loc_l)
    IF (IsA)  DEALLOCATE(loc_A)
    IF (IsB)  DEALLOCATE(loc_B)
    IF (IsC)  DEALLOCATE(loc_C)
    IF (IsD)  DEALLOCATE(loc_D)
  END FUNCTION GIVE_BACK_UNIQUE_TYPES

  CHARACTER(LEN=4) FUNCTION GIVE_BACK_TYPE(my_index) RESULT(my_type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: my_index

    my_type = ISYMBL(my_index/3+1)

  END FUNCTION GIVE_BACK_TYPE

END PROGRAM leap2fist

