!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/atoms_input [1.0] *
!!
!!   NAME
!!     atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!   cjm, Feb-20-2001 : added all the extended variables to
!!   system_type
!!   gt 23-09-2002 : major changes. Pointer part is allocated/deallocated
!!                   and initialized here. Atomic coordinates can now be 
!!                   read also from &COORD section in the input file.
!!                   If &COORD is not found, .dat file is read.  
!!                   If & coord is found and .NOT. 'INIT', parsing of the .dat
!!                   is performed to get the proper coords/vel/eta variables
!!    CJM 31-7-03  : Major rewrite.  No more atype
!!
!!   SOURCE
!******************************************************************************

MODULE atoms_input

  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             test_next
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             scaled_to_real, &
                                             get_hinv
  USE string_utilities,                ONLY: xstring
  USE termination,                     ONLY: stop_program
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             atom_info_type
  USE util,                            ONLY: get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atoms_input, read_md_restart, read_pimd_atoms_input

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** atoms_input/read_atoms_input [1.0] *
!!
!!   NAME
!!     read_atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!  OPTIONS
!!
!!  INPUTS
!!*** *************************************************************************

SUBROUTINE read_atoms_input ( topology, globenv )


    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    INTEGER, PARAMETER                       :: nblock_atom = 10000

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atom_info_type),POINTER             :: atom_info
    CHARACTER(LEN=6)                         :: label
    REAL(KIND=dp)                           :: unit_of_length
    INTEGER                                  :: atm_equal_mol, iatom, ierror, &
                                                ilen, iw, natom, ia, ie
    LOGICAL                                  :: scaled_coordinates

    LOGICAL                                  :: PM1, PM2, PM3, PM4

!   ---------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering read_atoms_input"


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. get cell and topolgoy%atom_info
  !-----------------------------------------------------------------------------
  atom_info => topology%atom_info
  cell => topology%cell
  CALL get_cell(cell=cell,&
                scaled_coordinates=scaled_coordinates,&
                unit_of_length=unit_of_length)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. NULLIFY and ALLOCATE all that is necessary in topology
  !-----------------------------------------------------------------------------
  CALL reallocate(atom_info%label_molname,1,nblock_atom)
  CALL reallocate(atom_info%label_resname,1,nblock_atom)
  CALL reallocate(atom_info%label_resid,1,nblock_atom)
  CALL reallocate(atom_info%label_atmname,1,nblock_atom)
  CALL reallocate(atom_info%element,1,nblock_atom)
  CALL reallocate(atom_info%r,1,3,1,nblock_atom)
  CALL reallocate(atom_info%atm_mass,1,nblock_atom)
  CALL reallocate(atom_info%atm_charge,1,nblock_atom)




  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Read in the coordinates from &COORD section in the input file
  !-----------------------------------------------------------------------------
  natom = 0
  label = '&COORD'
  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
      CALL stop_program ('read_atoms_input','No &COORD section found')
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')

       natom = natom + 1

       IF (natom > SIZE(atom_info%label_atmname)) THEN
        CALL reallocate(atom_info%label_molname,1,natom+nblock_atom)
        CALL reallocate(atom_info%label_resname,1,natom+nblock_atom)
        CALL reallocate(atom_info%label_resid,1,natom+nblock_atom)
        CALL reallocate(atom_info%label_atmname,1,natom+nblock_atom)
        CALL reallocate(atom_info%element,1,natom+nblock_atom)
        CALL reallocate(atom_info%r,1,3,1,natom+nblock_atom)
        CALL reallocate(atom_info%atm_mass,1,natom+nblock_atom)
        CALL reallocate(atom_info%atm_charge,1,natom+nblock_atom)
       END IF

       ilen = default_string_length
       atom_info%label_molname(natom) = ""
       atom_info%label_resname(natom) = ""
       atom_info%label_resid(natom) = 1
       atom_info%label_atmname(natom) = ""
       atom_info%element(natom) = ""

       ! For default, set atom name to residue name to molecule name
       CALL get_next(atom_info%label_atmname(natom),ilen)
       atom_info%label_molname(natom) = atom_info%label_atmname(natom)
       atom_info%label_resname(natom) = atom_info%label_atmname(natom)
       atom_info%element(natom) = atom_info%label_atmname(natom)

       ! *** Read the atomic coordinates ***
       CALL get_next(atom_info%r(1,natom))
       CALL get_next(atom_info%r(2,natom))
       CALL get_next(atom_info%r(3,natom))

       !Future user optional inputs
       IF (test_next()=='C') THEN
         CALL get_next(atom_info%label_molname(natom),ilen)
         atom_info%label_resname(natom) = atom_info%label_molname(natom)
       END IF
       atom_info%atm_mass(natom) = 0.0_dp 
       atom_info%atm_charge(natom) = 0.0_dp 

     call read_line
     END DO

    CALL parser_end
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Get rid of unnecessary spaces
  !-----------------------------------------------------------------------------
  topology%natoms = natom
  CALL reallocate(atom_info%label_molname,1,natom)
  CALL reallocate(atom_info%label_resname,1,natom)
  CALL reallocate(atom_info%label_resid,1,natom)
  CALL reallocate(atom_info%label_atmname,1,natom)
  CALL reallocate(atom_info%r,1,3,1,natom)
  CALL reallocate(atom_info%atm_mass,1,natom)
  CALL reallocate(atom_info%atm_charge,1,natom)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Convert coordinates into internal cp2k coordinates
  !-----------------------------------------------------------------------------
  DO iatom=1,natom
    IF (scaled_coordinates) THEN
      atom_info%r(:,iatom) = scaled_to_real(atom_info%r(:,iatom),cell)
    ELSE
      atom_info%r(:,iatom) = atom_info%r(:,iatom)*unit_of_length
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Post processing for special molecule name at the end of xyz
  !-----------------------------------------------------------------------------
  atm_equal_mol = 0
  DO iatom=1,natom
    IF(atom_info%label_molname(iatom)==atom_info%label_atmname(iatom))THEN
      atm_equal_mol = atm_equal_mol + 1
      topology%bondparm_factor=0.0_dp
    END IF
  END DO
  IF ((atm_equal_mol /= 0) .AND. (atm_equal_mol /= natom)) THEN
    CALL stop_program ('read_atoms_input','check coord optional molname')
  END IF


  IF(PM1) WRITE(iw,*) "  Exiting read_atoms_input"

END SUBROUTINE read_atoms_input

!******************************************************************************
!!****** atoms_input/read_md_restart [1.0] *
!!
!!   NAME
!!     read_md_restart
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!*** *************************************************************************
SUBROUTINE read_md_restart ( rtype, globenv, particle_set, cell, cell_ref,  &
                             nhc_part, nhc_baro, npt_info ) 
    CHARACTER(len=*), INTENT(IN)             :: rtype
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part, nhc_baro
    TYPE(npt_info_type), POINTER             :: npt_info( :, : )

    CHARACTER(len=20)                        :: filen
    INTEGER                                  :: i, ia, ie, iunit, iw, j, &
                                                size_pnhc1, size_pnhc2,  &
                                                size_bnhc1, size_bnhc2,  &
                                                size_baro1, size_baro2,  &
                                                icount, isos
    REAL (KIND=dp), POINTER, DIMENSION ( : )        :: eta_p, veta_p, fnhc_p, mnhc_p
    REAL (KIND=dp), POINTER, DIMENSION ( : )        :: eta_b, veta_b, fnhc_b, mnhc_b
    REAL (KIND=dp), POINTER, DIMENSION ( :,: )      :: veps, meps, h, href, r, v, f
    LOGICAL                                  :: exists

  iw = globenv % scr

  NULLIFY ( eta_p, veta_p, fnhc_p, mnhc_p )
  NULLIFY ( eta_b, veta_b, fnhc_b, mnhc_b )
  NULLIFY ( r, v, f, h, href )
  NULLIFY ( veps, meps )

  CALL xstring ( globenv % project_name, ia, ie )
  filen = globenv % project_name ( ia:ie ) // '.restart'
  ia = MIN ( LEN ( filen ), 20 )
  IF (globenv % ionode) &
    WRITE ( iw, '( A, T61, A )' ) ' READ MD| Read pos/vel/eta from file ', &
           ADJUSTR ( filen ( 1:ia ) )

  IF ( globenv % ionode ) THEN
    ALLOCATE ( r ( 3, SIZE ( particle_set ) ), STAT = isos )
    IF (isos/=0) CALL stop_program('atoms_input', &
    'failed to allocate r' )
    ALLOCATE ( v ( 3, SIZE ( particle_set ) ), STAT = isos )
    IF (isos/=0) CALL stop_program('atoms_input', &
    'failed to allocate v' )
    ALLOCATE ( f ( 3, SIZE ( particle_set ) ), STAT = isos )
    IF (isos/=0) CALL stop_program('atoms_input', &
    'failed to allocate f' )
    ALLOCATE ( h ( 3, 3 ), STAT = isos )
    IF (isos/=0) CALL stop_program('atoms_input', &
    'failed to allocate h' )
    ALLOCATE ( href ( 3, 3 ), STAT = isos )
    IF (isos/=0) CALL stop_program('atoms_input', &
    'failed to allocate href' )
    INQUIRE ( FILE = filen, EXIST = exists )
    IF ( exists ) THEN
       iunit = get_unit()
       OPEN ( iunit, file = filen )
       IF ( rtype == 'INIT' ) THEN
          DO i = 1, SIZE ( particle_set )
             READ ( iunit, * ) r ( :, i )
             f ( :, i ) = 0.0_dp
          END DO
          READ ( iunit, * ) h ( 1, 1:3 )
          READ ( iunit, * ) h ( 2, 1:3 )
          READ ( iunit, * ) h ( 3, 1:3 )

       ELSE IF ( rtype == 'POS' )THEN
          DO i = 1, SIZE ( particle_set )
             READ ( iunit, * ) r ( :, i )
             f ( :, i ) = 0.0_dp
          END DO
          READ ( iunit, * ) h( 1, 1:3 )
          READ ( iunit, * ) h( 2, 1:3 )
          READ ( iunit, * ) h( 3, 1:3 )

          READ ( iunit, * ) href( 1, 1:3 )
          READ ( iunit, * ) href( 2, 1:3 )
          READ ( iunit, * ) href( 3, 1:3 )
       ELSE IF ( rtype == 'POSVEL' ) THEN

          DO i = 1, SIZE ( particle_set )
             READ ( iunit, * ) r ( :, i )
             f ( :, i ) = 0.0_dp
          END DO
          READ ( iunit, * ) h( 1, 1:3 )
          READ ( iunit, * ) h( 2, 1:3 )
          READ ( iunit, * ) h( 3, 1:3 )

          READ ( iunit, * ) href( 1, 1:3 )
          READ ( iunit, * ) href( 2, 1:3 )
          READ ( iunit, * ) href( 3, 1:3 )

          DO i = 1, SIZE ( particle_set )
             READ ( iunit, * ) v ( :, i )
          END DO

       ELSE IF ( rtype == 'ALL' ) THEN

          DO i = 1, SIZE ( particle_set )
             READ ( iunit, * ) r ( :, i )
             f ( :, i ) = 0.0_dp
          END DO

          READ ( iunit, * ) h ( 1, 1:3 )
          READ ( iunit, * ) h ( 2, 1:3 )
          READ ( iunit, * ) h ( 3, 1:3 )

          READ ( iunit, * ) href ( 1, 1:3 )
          READ ( iunit, * ) href ( 2, 1:3 )
          READ ( iunit, * ) href ( 3, 1:3 )


          DO i = 1, SIZE ( particle_set )
             READ ( iunit, * ) v ( :, i )
          END DO

! particle thermostat positions and velocities
          READ ( iunit, * ) size_pnhc1, size_pnhc2

          ALLOCATE ( eta_p ( size_pnhc1*size_pnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate eta_p' )

          ALLOCATE ( veta_p ( size_pnhc1*size_pnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate veta_p' )

          ALLOCATE ( fnhc_p ( size_pnhc1*size_pnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate fnhc_p' )

          ALLOCATE ( mnhc_p ( size_pnhc1*size_pnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate mnhc_p' )

          DO i = 1, size_pnhc1*size_pnhc2
             READ ( iunit, * ) eta_p ( i )
          END DO
          DO i = 1, size_pnhc1*size_pnhc2
             READ ( iunit, * ) veta_p ( i )
          END DO
          DO i = 1, size_pnhc1*size_pnhc2
             READ ( iunit, * ) fnhc_p ( i )
          END DO
          DO i = 1, size_pnhc1*size_pnhc2
             READ ( iunit, * ) mnhc_p ( i )
          END DO
! barostat thermostat positions and velocities
          READ ( iunit, * ) size_bnhc1, size_bnhc2

          ALLOCATE ( eta_b ( size_bnhc1*size_bnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate eta_b' )

          ALLOCATE ( veta_b ( size_bnhc1*size_bnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate veta_b' )

          ALLOCATE ( fnhc_b ( size_bnhc1*size_bnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate fnhc_b' )

          ALLOCATE ( mnhc_b ( size_bnhc1*size_bnhc2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate mnhc_b' )

          DO i = 1, size_bnhc1*size_bnhc2
             READ ( iunit, * ) eta_b ( i )
          END DO
          DO i = 1, size_bnhc1*size_bnhc2
             READ ( iunit, * ) veta_b ( i )
          END DO
          DO i = 1, size_bnhc1*size_bnhc2
             READ ( iunit, * ) fnhc_b ( i )
          END DO
          DO i = 1, size_bnhc1*size_bnhc2
             READ ( iunit, * ) mnhc_b ( i )
          END DO
! barostat velocities
          READ ( iunit, * ) size_baro1, size_baro2

          ALLOCATE ( veps ( size_baro1, size_baro2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate veps' )

          ALLOCATE ( meps ( size_baro1, size_baro2 ), STAT = isos )
          IF (isos/=0) CALL stop_program('atoms_input', &
          'failed to allocate meps' )
          DO i = 1, size_baro1
            DO j = 1, size_baro2
              READ ( iunit, * ) veps ( i, j )
            ENDDO
          ENDDO
          DO i = 1, size_baro1
            DO j = 1, size_baro2
              READ ( iunit, * )  meps ( i, j )
            ENDDO
          ENDDO
       ELSE
          CALL stop_program ( 'read_md_restart', 'this RESTART type not programmed' )
       END IF

       CLOSE ( iunit )

    ELSE

        CALL stop_program ( 'read_md_restart', 'No MD RESTART file found ' )

    END IF
  END IF

!-----------------------------------------------------------------------------
! Broadcast data and assign to structures in parallel
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! particles, velocities, forces
!-----------------------------------------------------------------------------
  CALL reallocate ( r, 1, 3, 1, SIZE ( particle_set ) )
  CALL mp_bcast ( r, globenv % source, globenv % group )

  CALL reallocate ( v, 1, 3, 1, SIZE ( particle_set ) )
  CALL mp_bcast ( v, globenv % source, globenv % group )

  CALL reallocate ( f, 1, 3, 1, SIZE ( particle_set ) )
  CALL mp_bcast ( f, globenv % source, globenv % group )

  DO i = 1, SIZE ( particle_set )
    particle_set ( i ) % r ( : ) = r ( :, i )
    particle_set ( i ) % v ( : ) = v ( :, i )
    particle_set ( i ) % f ( : ) = f ( :, i )
  END DO

  DEALLOCATE ( r, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate r' )
  DEALLOCATE ( v, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate v' )
  DEALLOCATE ( f, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate f' )

!-----------------------------------------------------------------------------
! cells
!-----------------------------------------------------------------------------
  CALL reallocate ( h, 1, 3, 1, 3 )
  CALL mp_bcast ( h, globenv % source, globenv % group )

  CALL reallocate ( href, 1, 3, 1, 3 )
  CALL mp_bcast ( href, globenv % source, globenv % group )
  DO i = 1, 3
    DO j = 1, 3
      cell % hmat ( i, j ) = h ( i, j )
      cell_ref % hmat ( i, j ) = href ( i, j )
    END DO
  END DO

  CALL get_hinv ( cell_ref )
  CALL get_hinv ( cell )

  DEALLOCATE ( h, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate h' )
  DEALLOCATE ( href, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate href' )

!-----------------------------------------------------------------------------
! particle thermostat
!-----------------------------------------------------------------------------
  CALL mp_bcast ( size_pnhc1, globenv % source, globenv % group )
  CALL mp_bcast ( size_pnhc2, globenv % source, globenv % group )
  
  CALL reallocate ( eta_p, 1, size_pnhc1*size_pnhc2 )
  CALL mp_bcast ( eta_p, globenv % source, globenv % group )

  CALL reallocate ( veta_p, 1, size_pnhc1*size_pnhc2 )
  CALL mp_bcast ( veta_p, globenv % source, globenv % group )

  CALL reallocate ( fnhc_p, 1, size_pnhc1*size_pnhc2 )
  CALL mp_bcast ( fnhc_p, globenv % source, globenv % group )

  CALL reallocate ( mnhc_p, 1, size_pnhc1*size_pnhc2 )
  CALL mp_bcast ( mnhc_p, globenv % source, globenv % group )

  DO i = 1, nhc_part % num_nhc
     icount = nhc_part % index(i)
     icount = (icount-1) * nhc_part % nhc_len
     DO j = 1, nhc_part % nhc_len
        icount = icount + 1
        nhc_part % nvt(j,i) % eta = eta_p ( icount )
        nhc_part % nvt(j,i) % v = veta_p ( icount )
        nhc_part % nvt(j,i) % f = fnhc_p ( icount )
        nhc_part % nvt(j,i) % mass = mnhc_p ( icount )
     END DO
  END DO

  DEALLOCATE ( eta_p, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate eta_p' )
  DEALLOCATE ( veta_p, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate veta_p' )
  DEALLOCATE ( fnhc_p, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate fnhc_p' )
  DEALLOCATE ( mnhc_p, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate mnhc_p' )

!-----------------------------------------------------------------------------
! baroastat thermostat
!-----------------------------------------------------------------------------
  CALL mp_bcast ( size_bnhc1, globenv % source, globenv % group )
  CALL mp_bcast ( size_bnhc2, globenv % source, globenv % group )
  
  CALL reallocate ( eta_b, 1, size_bnhc1*size_bnhc2 )
  CALL mp_bcast ( eta_b, globenv % source, globenv % group )

  CALL reallocate ( veta_b, 1, size_bnhc1*size_bnhc2 )
  CALL mp_bcast ( veta_b, globenv % source, globenv % group )

  CALL reallocate ( fnhc_b, 1, size_bnhc1*size_bnhc2 )
  CALL mp_bcast ( fnhc_b, globenv % source, globenv % group )

  CALL reallocate ( mnhc_b, 1, size_bnhc1*size_bnhc2 )
  CALL mp_bcast ( mnhc_b, globenv % source, globenv % group )

  DO i = 1, nhc_baro % num_nhc
     icount = nhc_part % index(i)
     icount = (icount-1) * nhc_part % nhc_len
     DO j = 1, nhc_baro % nhc_len
        icount = icount + 1
        nhc_baro % nvt(j,i) % eta = eta_b ( icount )
        nhc_baro % nvt(j,i) % v = veta_b ( icount )
        nhc_baro % nvt(j,i) % f = fnhc_b ( icount )
        nhc_baro % nvt(j,i) % mass = mnhc_b ( icount )
     END DO
  END DO

  DEALLOCATE ( eta_b, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate eta_b' )
  DEALLOCATE ( veta_b, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate veta_b' )
  DEALLOCATE ( fnhc_b, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate fnhc_b' )
  DEALLOCATE ( mnhc_b, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate mnhc_b' )
!------------------------------------------------------------------------
!barostat
!------------------------------------------------------------------------
  CALL mp_bcast ( size_baro1, globenv % source, globenv % group )
  CALL mp_bcast ( size_baro2, globenv % source, globenv % group )
  
  CALL reallocate ( veps, 1, size_baro1, 1, size_baro2 )
  CALL mp_bcast ( veps, globenv % source, globenv % group )

  CALL reallocate ( meps, 1, size_baro1, 1, size_baro2 )
  CALL mp_bcast ( meps, globenv % source, globenv % group )

  IF ( size_baro1 /= SIZE ( npt_info, 1 ) ) &
  CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
  IF ( size_baro2 /= SIZE ( npt_info, 2 ) ) &
  CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )

  DO i = 1, size_baro1
    DO j = 1, size_baro2
      npt_info ( i, j ) % v = veps ( i, j )
    ENDDO
  ENDDO
  DO i = 1, size_baro1
    DO j = 1, size_baro2
      npt_info ( i, j ) % mass = meps ( i, j )
    ENDDO
  ENDDO

  DEALLOCATE ( veps, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate veps' )

  DEALLOCATE ( meps, STAT = isos )
  IF (isos/=0) CALL stop_program('atoms_input', &
  'failed to allocate meps' )
   


!..write some information to output
  IF ( globenv%ionode .AND. globenv%print_level>0 ) THEN
    WRITE ( iw, '( A )' ) ' RESTART MD | Box parameters [Angstrom]'
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 1, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 2, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 3, i ), i = 1, 3 )
    WRITE ( iw, '( A, T71, I10 )' ) &
        ' RESTART MD| Number of atoms read ', SIZE ( particle_set )
    IF ( globenv % print_level > 4 ) THEN
         IF ( rtype == 'POS' ) THEN
           CALL print_c ( iw, particle_set )
         ELSE IF ( rtype == 'POSVEL' ) THEN
            CALL print_cv ( iw, particle_set )
         END IF
    END IF
    WRITE ( iw, '()' )
  END IF

END SUBROUTINE read_md_restart

!******************************************************************************
!!****** atoms_input/read_pimd_atoms_input [1.0] *
!!
!!   NAME
!!     read_pimd_atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!*** *************************************************************************
SUBROUTINE read_pimd_atoms_input ( rtype, globenv, particle_set, cell,  &
                                  cell_ref, nhc_part, nhc_baro, npt_info ) 

    CHARACTER(len=*), INTENT(IN)             :: rtype
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part, nhc_baro
    TYPE(npt_info_type), POINTER             :: npt_info( :, : )

    CHARACTER(len=20)                        :: filen
    INTEGER                                  :: i, ia, ie, iunit, iw, j, &
                                                size1, size2
    LOGICAL                                  :: exists

  iw = globenv % scr

  CALL xstring ( globenv % project_name, ia, ie )
  filen = globenv % project_name ( ia:ie ) // '.dat'
  ia = MIN ( LEN ( filen ), 20 )
  WRITE ( iw, '( A, T61, A )' ) ' READ_PIMD| Read pos/vel/eta from file ', &
           ADJUSTR ( filen ( 1:ia ) )

  INQUIRE ( FILE = filen, EXIST = exists )
  IF ( exists ) THEN
     iunit = get_unit()
     OPEN ( iunit, file = filen )

     IF ( rtype == 'INIT' ) THEN
        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0.0_dp
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

     ELSE IF ( rtype == 'POS' )THEN
        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0.0_dp
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )
     ELSE IF ( rtype == 'POSVEL' ) THEN

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % v ( : )
        END DO

     ELSE IF ( rtype == 'ALL' ) THEN

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0.0_dp
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % v ( : )
        END DO

! particle thermostat positions and velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( nhc_part % nvt, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_part' )
        IF ( size2 /= SIZE ( nhc_part % nvt, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_part' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % eta 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % v 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % mass 
          ENDDO
        ENDDO
! barostat thermostat positions and velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( nhc_baro % nvt, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_baro' )
        IF ( size2 /= SIZE ( nhc_baro % nvt, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_baro' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % eta 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % v 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % mass 
          ENDDO
        ENDDO
! barostat velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( npt_info, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
        IF ( size2 /= SIZE ( npt_info, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) npt_info ( i, j ) % v
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * )  npt_info ( i, j ) % mass 
          ENDDO
        ENDDO
     ELSE
        CALL stop_program ( 'read_md_restart', 'this RESTART type not programmed' )
     END IF

     CLOSE ( iunit )

  ELSE

     CALL stop_program ( 'read_md_restart', 'No MD RESTART file found ' )

  END IF

!..write some information to output
 IF ( globenv%ionode .AND. globenv%print_level>4 ) THEN
    WRITE ( iw, '( A )' ) ' READ PIMD| Box parameters [Angstrom]'
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' READ PIMD| ', ( cell % hmat ( 1, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 2, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' READ PIMD| ', ( cell % hmat ( 3, i ), i = 1, 3 )
    WRITE ( iw, '( A, T71, I10 )' ) &
        ' READ PIMD| Number of atoms read ', SIZE ( particle_set )
    IF ( globenv % print_level > 4 ) THEN
         IF ( rtype == 'POS' ) THEN
           CALL print_c ( iw, particle_set )
         ELSE IF ( rtype == 'POSVEL' ) THEN
            CALL print_cv ( iw, particle_set )
         END IF
    END IF
    WRITE ( iw, '()' )
  END IF

END SUBROUTINE read_pimd_atoms_input

!******************************************************************************

SUBROUTINE print_c ( iw, particle_set )


    INTEGER, INTENT(IN)                      :: iw
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )

    INTEGER                                  :: i, n

!------------------------------------------------------------------------------

  n = SIZE ( particle_set )
  WRITE ( iw, '( A )' ) ' RESTART MD | Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T26, I10, 3F15.5 )' ) ' RESTART MD| ', i, particle_set ( i ) % r ( : )
  END DO

END SUBROUTINE print_c

!******************************************************************************

SUBROUTINE print_cv ( iw, particle_set )


    INTEGER, INTENT(IN)                      :: iw
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )

    INTEGER                                  :: i, n

!------------------------------------------------------------------------------

  n = SIZE ( particle_set )
  WRITE ( iw, '( A )' ) ' RESTART MD | Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T8, I8, 3F10.4, 5X, 3F10.4 )' ) &
          ' RESTART MD| ', i, particle_set ( i ) % r ( : ), &
                        particle_set ( i ) % v ( : )
  END DO

END SUBROUTINE print_cv

!******************************************************************************

END MODULE atoms_input
