!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/nhc [1.0] *
!!
!!   NAME
!!     nhc
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM, 20-Feb-01
!!     JGH (10-Mar-2001) 
!!     CJM (10-Apr-2001) 
!!
!!   SOURCE
!******************************************************************************

MODULE nhc

  USE coefficient_types, ONLY : coeff_type, PLANEWAVES, ATOMICORBITALS
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : molecule_structure_type, molecule_type
  USE message_passing, ONLY : mp_sum
  USE nose, ONLY : npt_info_type, lnhc_parameters_type, yoshida_coef
  USE util, ONLY : get_share
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nhc_particle_setup, nhc_barostat_setup, nhc_coef_setup

!!***
!******************************************************************************
CONTAINS
!******************************************************************************
!!****** nhc/nhc_barostat_setup [1.0] *
!!
!!   NAME
!!     nhc_barostat
!!
!!   FUNCTION
!!     Creates the thermostatting for the barostat
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM, 20-Feb-01  : nhc structure allocated to zero when not in use
!!     JGH (10-Mar-2001) : set nhc variables to zero when not in use
!!
!!*** *************************************************************************

SUBROUTINE nhc_barostat_setup ( simpar, nhc, globenv )
! creates the thermostatting maps use global_mod
  IMPLICIT NONE
! begin global declaratation
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( lnhc_parameters_type ), INTENT ( OUT ) :: nhc
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
! begin local definitions
  INTEGER :: i, iw
  INTEGER :: isos, number
!
  iw = globenv % scr
!
  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
  CASE  DEFAULT
    ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
       'failed to allocate nhc%nvt' )
    nhc % nyosh = 0
    nhc % nc = 0
    nhc % nhc_len = 0
    nhc % num_nhc = 0
  CASE ( 'NPT' )
    nhc % nyosh = simpar % nyosh
    nhc % nc = simpar % nc
    nhc % nhc_len = simpar % nhclen

! Counting the total number of thermostats ( 1 for both NPT_I and
! NPT_F )

    nhc % num_nhc = 1

    IF ( simpar % ensemble == 'NPT_F' ) THEN
      number = 9
    ELSE
      number = 1
    ENDIF

! Now we can allocate the target array s_kin and p_kin

    ALLOCATE ( nhc % s_kin ( nhc % num_nhc ), STAT = isos)
    IF ( isos /= 0 ) CALL stop_program ( 'nhc_barostat_setup', &
      'failed to allocate nhc%s_kin' )
    ALLOCATE ( nhc % v_scale ( nhc % num_nhc ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
      'failed to allocate nhc%v_scale' )
    ALLOCATE ( nhc % p_kin ( 1, number ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
      'failed to allocate nhc%p_kin' )
    ALLOCATE ( nhc % p_scale ( 1, number ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
      'failed to allocate nhc%p_scale')

! Begin the mapping loop

    DO i = 1, number
      nhc % p_kin ( 1, i ) % point => nhc % s_kin ( 1 )
      nhc % p_scale ( 1, i ) % point => nhc % v_scale ( 1 )
    END DO

    ALLOCATE ( nhc % nvt ( nhc % nhc_len, nhc % num_nhc ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
       'failed to allocate nhc%nvt' )
!
! Now that we know how many there are stick this into nhc % static % nkt
! (number of degrees of freedom times k_B T )

    DO i = 1, nhc % num_nhc
      nhc % nvt ( 1, i ) % nkt = simpar % temp_ext * number
! dbg ISOTROPIC LIMIT
!      nhc % nvt ( 1, i ) % nkt = simpar % temp_ext
! dbg
    END DO

! getting the number of degrees of freedom times k_B T for the rest of the chain

    DO i = 2, nhc % nhc_len
      nhc % nvt ( i, : ) %nkt = simpar % temp_ext
    END DO

! Set up the Yoshida weights

    ALLOCATE (nhc % dt_yosh (1 : nhc % nyosh ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
      'failed to allocate nhc%dt_yosh')
    CALL yoshida_coef( nhc, simpar % dt )
  END SELECT


END SUBROUTINE nhc_barostat_setup
!******************************************************************************
!!****** nhc/nhc_particle_setup [1.0] *
!!
!!   NAME
!!     nhc_particle_setup
!!
!!   FUNCTION
!!     Creates the thermostatting maps
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     29-Nov-00 (JGH) correct counting of DOF if constraints are off
!!     CJM, 20-Feb-01  : nhc structure allocated to zero when not in use
!!     JGH (10-Mar-2001) : set nhc variables to zero when not in use
!!
!!*** *************************************************************************

SUBROUTINE nhc_particle_setup ( simpar, molecule, mol_setup, nhc, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( INOUT ) :: mol_setup
  TYPE ( lnhc_parameters_type ), INTENT ( OUT ) :: nhc
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  CHARACTER ( LEN = 20 ) :: opt
  INTEGER :: i, j, k, nmol_type, number, iw, natoms
  INTEGER :: const_per_molecule, nmol
  INTEGER :: atm_offset, mol_offset, start, finish, iconst, iglobal
  INTEGER :: icust, ii, jj, kk, ieq, isos
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: mnum, tnum, deg_of_freedom
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: const_mol, tot_const
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: natoms_mol, mol_off
  INTEGER, DIMENSION ( :, : ), ALLOCATABLE :: point
  REAL ( dbl ), SAVE, TARGET :: scale
  LOGICAL :: match, noneflag, nointer

!------------------------------------------------------------------------------

  iw = globenv % scr

  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
    CASE ( 'NVE' )
     ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%nvt', 0 )
     nhc % nyosh = 0
     nhc % nc = 0
     nhc % nhc_len = 0
     nhc % num_nhc = 0
    CASE ( 'NVT', 'NPT' )
     nhc%nyosh = simpar%nyosh
     nhc%nc = simpar%nc
     nhc%nhc_len = simpar%nhclen

! allocating the temporary arrays mnum and tnum
     nmol_type = size(mol_setup)
     nmol = size(molecule)
     ALLOCATE (mnum(nmol_type ),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'mnum', nmol_type )
     ALLOCATE (tnum(nmol_type ),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tnum', nmol_type )
     ALLOCATE (const_mol(nmol),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol', nmol )
     ALLOCATE (tot_const(nmol),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const', nmol )
     ALLOCATE (natoms_mol(nmol),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'natoms_mol', nmol )
     ALLOCATE (point(2,nmol),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','point', 2 * nmol )
     ALLOCATE (mol_off(nmol_type+1),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','mol_off', nmol_type + 1 )

! Counting the total number of thermostats
     number = 0
     match = .FALSE.
     noneflag = .FALSE.
     nointer = .TRUE.

! 1) count the 'global'
! (In principle, every moleule type can have GLOBAL but
! there is only ONE thermostat.  So, only count it once)
    DO i = 1, nmol_type
       IF ((mol_setup(i)%nhcopt=="GLOBAL") .AND. ( .NOT. match)) THEN
          number = number + 1
          match = .TRUE.
          nointer = .FALSE.
       END IF
    END DO

! 2) count the user defined options
    DO i = 1, nmol_type
       opt = mol_setup(i) %nhcopt
       IF ((opt/="GLOBAL") .AND. (opt/="MASSIVE") .AND. &
            (opt/="MOLECULE") .AND. (opt/="NONE")) THEN
! checking to see if 2 or more molecules are on the same user
! defined thermostat.  Again, only one thermostat per
! unique user defined option
          match = .FALSE.
          CHECK: DO j = 1, i - 1
             IF (opt==mol_setup(j)%nhcopt) THEN
                match = .TRUE.
                EXIT CHECK
             END IF
          END DO CHECK
          IF ( .NOT. match) number = number + 1
          IF ( .NOT. match) nointer = .FALSE.
       END IF
    END DO

! 3) count the 'molecule'
    DO i = 1, nmol_type
       IF (mol_setup(i)%nhcopt=="MOLECULE") THEN
          number = number + mol_setup(i) %num_mol
       END IF
    END DO

! 4) count the 'massive'
    DO i = 1, nmol_type
       IF (mol_setup(i)%nhcopt=="MASSIVE") THEN
          number = number + 3*mol_setup(i) %num_mol*mol_setup(i) %molpar%natom
       END IF
    END DO
    IF ( number == 0 ) THEN
       CALL stop_program ( 'nhc_paticle_setup', "Number of NHC = 0 -> use NVE" )
    END IF

!
! here we decide on what parallel algorithm to use
! if there are only massive and molecule type thermostats we can use
! a local scheme, in all other cases we assume a coupling of
! degrees of freedom from different processors
!
    nhc%dis_type = 'NONE'
    IF (globenv%num_pe>1) THEN
     IF (nointer) THEN

! no interaction between thermostats on different processors
        nhc%dis_type = 'DIS_REP'

! we only count thermostats on this processor
        number = 0
        DO i = 1, size(molecule)
           j = molecule(i) %moltype
           IF (mol_setup(j)%nhcopt=="MOLECULE") THEN
              number = number + 1
           ELSE IF (mol_setup(j)%nhcopt=="MASSIVE") THEN
              number = number + 3*mol_setup(j) %molpar%natom
           ELSE
              CALL stop_program ( 'nhc_particle_setup', '' )
           END IF
        END DO
     ELSE
! interaction between thermostats on different processors
        nhc%dis_type = 'INTER_REP'
     END IF
    END IF

! some utility arrays for the setup of thermostats
! Define the local variable point(i,j) i=1,2 and j=1->nmol_type.
! This array stores the locations of the first and last atom in
! of a molecule type.  It will be used to map a specific atom
! belonging to a molecule onto a thermostat.  nconst(j) j=1->nmol_type
! stores the TOTAL number of constraints of each molecule type to
! compute the degrees of freedom.
    IF (nhc%dis_type=='NONE' .OR. nhc%dis_type=='DIS_REP' ) THEN
       atm_offset = 0
       DO i = 1, size(molecule)
          j = molecule(i) %moltype
          natoms_mol(i) = mol_setup(j) %molpar%natom
          point(1,i) = atm_offset + 1
          point(2,i) = atm_offset + natoms_mol ( i )
          const_mol(i) = mol_setup(j) %molpar%ndcon &
             + 3*mol_setup(j) %molpar%n3x3con
          tot_const(i) = const_mol(i)
          atm_offset = point(2,i)
       END DO

    ELSE IF (nhc%dis_type=='INTER_REP' ) THEN
       atm_offset = 0
       mol_off = 0
       DO i = 1, nmol_type
          j = molecule ( i ) % moltype
          natoms_mol(i) = mol_setup(j) %molpar%natom
          point(1,i) = atm_offset + 1
          point(2,i) = atm_offset + natoms_mol(i)*get_share(mol_setup(j)%num_mol, &
               globenv%num_pe,globenv%mepos)
          const_mol(i) = mol_setup(j) %molpar%ndcon + 3*mol_setup(j) %molpar%n3x3con
          tot_const(i) = const_mol(i)*mol_setup(j) %num_mol
          atm_offset = point(2,i)
       END DO
       DO i = 1, nmol_type
          mol_off(i+1) = mol_off(i) + get_share(mol_setup(i)%num_mol,globenv%num_pe, &
               globenv%mepos)
       END DO
    END IF

    IF ( .NOT. simpar % constraint ) THEN
      const_mol = 0._dbl
      tot_const = 0._dbl
    END IF

! count the number of atoms on this processor
    natoms = sum(molecule ( : )%natoms_mol)

! Now we can allocate the target array s_kin and p_kin
    SELECT CASE (simpar%ensemble(1:3))
    CASE ( 'NPT' )
       ALLOCATE (deg_of_freedom(number+1),STAT=isos)
       IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
          'deg_of_freedom', number + 1 )
    CASE DEFAULT
       ALLOCATE (deg_of_freedom(number),STAT=isos)
       IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
          'deg_of_freedom', number )
    END SELECT

    ALLOCATE (nhc%s_kin(number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%s_kin', number )
    ALLOCATE (nhc%v_scale(number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%v_scale', number )
    ALLOCATE (nhc%p_kin(3,natoms),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%p_kin', 3 * natoms )
    ALLOCATE (nhc%p_scale(3,natoms),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'nhc%p_scale', 3 * natoms )

! Check to see if there is at least one global thermostatting option defined
    iglobal = 0
    DO i = 1, nmol_type
       IF (mol_setup(i)%nhcopt=="GLOBAL") THEN
          iglobal = 1
          EXIT
       END IF
    END DO

! initialize variable 'scale' for the 'none' option
    scale = 1.0_dbl

! counter variable number is to determine total number of thermostats
! since the global is ALWAYS on thermostat one, if there is at least
! one GLOBAL option, the counter starts at 1 else, it starts at 0
    IF (iglobal==1) THEN
       number = 1
    ELSE
       number = 0
    END IF
    icust = 0

! initialize deg_of_freedom
    deg_of_freedom = 0
    IF (nhc%dis_type=="NONE" .OR. nhc%dis_type=="DIS_REP") THEN

!  Begin the mapping loop
     DO i = 1, size(molecule)
        j = molecule(i) %moltype
        opt = mol_setup(j) %nhcopt

! Search for user defined thermostatting options first
        IF ((opt/="GLOBAL") .AND. (opt/="MOLECULE") .AND. &
             (opt/="MASSIVE") .AND. (opt/="NONE")) THEN
           icust = icust + 1

! The first user defined option gets thrown on a thermostat, no questions asked
           IF (icust==1) THEN
              number = number + 1

! map the thermostat number to the molecule
              molecule(i) %nvt_num = number

! save the molecule number that is customized
              mnum(icust) = i

! save the thermostat number that is customized
              tnum(icust) = number
              DO ii = 1, 3
                 DO jj = point(1,i), point(2,i)
                    nhc%p_kin(ii,jj) %point => nhc%s_kin(number)
                    nhc%p_scale(ii,jj) %point => nhc%v_scale(number)
                 END DO
              END DO
              deg_of_freedom(number) = deg_of_freedom(number) + tot_const(i)
           ELSE

! if icust > 1 check to see if it is equal to the previous user-defined options
              ieq = 0
! The loop is to icust-1 because we know that it is equal to itself
              DO jj = 1, icust - 1
                 IF (mol_setup(j)%nhcopt==mol_setup(mnum(jj))%nhcopt) THEN
                    ieq = jj
                    EXIT
                 END IF
              END DO

! if the user defined option found a match, map everything !
!    on the same thermostat
              IF (ieq /= 0 ) THEN
! map the thermostat number to the molecule
                 molecule(i) %nvt_num = tnum(mnum(ieq))
                 DO ii = 1, 3
                    DO jj = point(1,i), point(2,i)
                       nhc%p_kin(ii,jj) %point => nhc%s_kin(tnum(mnum(ieq)))
                       nhc%p_scale(ii,jj) %point => nhc%v_scale(tnum(mnum( &
                            ieq)))
                    END DO
                 END DO
                 deg_of_freedom(tnum(mnum(ieq))) = deg_of_freedom(tnum(mnum( &
                      ieq))) + tot_const(i)

! if the user defined option does not find a match.  Map it to a new
!    thermostat, no questions asked
              ELSE
                 number = number + 1

! map the thermostat number to the molecule
                 molecule(i) %nvt_num = number
                 mnum(icust) = i
                 tnum(icust) = number
                 DO ii = 1, 3
                    DO jj = point(1,i), point(2,i)
                       nhc%p_kin(ii,jj) %point => nhc%s_kin(number)
                       nhc%p_scale(ii,jj) %point => nhc%v_scale(number)
                    END DO
                 END DO
                 deg_of_freedom ( number ) = deg_of_freedom ( number ) &
                      + tot_const ( i )
              END IF
           END IF

! doing global thermostat option
        ELSE IF (mol_setup(j)%nhcopt=="GLOBAL") THEN
! map the thermostat number to the molecule
           molecule(i) %nvt_num = 1
           DO ii = 1, 3
              DO jj = point(1,i), point(2,i)
                 nhc%p_kin(ii,jj) %point => nhc%s_kin(1)
                 nhc%p_scale(ii,jj) %point => nhc%v_scale(1)
              END DO
           END DO
           deg_of_freedom(1) = deg_of_freedom(1) + tot_const(i)

! doing molecular thermostat option
        ELSE IF (mol_setup(j)%nhcopt=="MOLECULE") THEN
           DO ii = point(1,i), point(2,i), natoms_mol(i)
              number = number + 1

! map the thermostat number to the molecule
              molecule(i) %nvt_num = number
              start = ii
              finish = ii + molecule(i) %natoms_mol - 1
              DO jj = 1, 3
                 DO kk = start, finish
                    nhc%p_kin(jj,kk) %point => nhc%s_kin(number)
                    nhc%p_scale(jj,kk) %point => nhc%v_scale(number)
                 END DO
              END DO
              deg_of_freedom(number) = const_mol(j)
           END DO

! doing massive thermostat option
        ELSE IF ( mol_setup ( j ) % nhcopt == "MASSIVE" ) THEN
           IF ( tot_const ( i ) > 0 ) THEN
              CALL stop_program ( 'nhc_particle_setup', &
                   "Massive thermostats with constraints not possible" )
           END IF
           DO ii = point(1,i), point(2,i)
              DO jj = 1, 3
                 number = number + 1
                 nhc%p_kin(jj,ii) %point => nhc%s_kin(number)
                 nhc%p_scale(jj,ii) %point => nhc%v_scale(number)
              END DO
           END DO

        ELSE IF (mol_setup(j)%nhcopt=="NONE") THEN
           noneflag = .TRUE.
! NOTE: the thermostat counter is not incremented here.
! vscale(i,j) points to scale=1.d0 for eternity.  Thus
! the the velocites of particles that are not connected
! to a thermostat will not get scaled in lnhc.f
           DO ii = 1, 3
              DO jj = point(1,i), point(2,i)
                 nhc%p_scale(ii,jj) %point => scale
              END DO
           END DO
        END IF
     END DO

    ELSE IF (nhc%dis_type=="INTER_REP") THEN

!  Begin the mapping loop
     DO i = 1, nmol_type
        opt = mol_setup(i) %nhcopt

! Search for user defined thermostatting options first
        IF ((opt/="GLOBAL") .AND. (opt/="MOLECULE") .AND. &
             (opt/="MASSIVE") .AND. (opt/="NONE")) THEN
           icust = icust + 1

! The first user defined option gets thrown on a thermostat, no questions asked
           IF (icust==1) THEN
              number = number + 1

! map the thermostat number to the molecule
              DO k = mol_off(i) + 1, mol_off(i+1)
                 molecule(k) %nvt_num = number
              END DO

! save the molecule number that is customized
              mnum(icust) = i

! save the thermostat number that is customized
              tnum(icust) = number
              DO ii = 1, 3
                 DO jj = point(1,i), point(2,i)
                    nhc%p_kin(ii,jj) %point => nhc%s_kin(number)
                    nhc%p_scale(ii,jj) %point => nhc%v_scale(number)
                 END DO
              END DO
              deg_of_freedom(number) = deg_of_freedom(number) + tot_const(i)
           ELSE

! if icust > 1 check to see if it is equal to the previous user-defined options
              ieq = 0

! The loop is to icust-1 because we know that it is equal to itself
              DO jj = 1, icust - 1
                 IF (mol_setup(i)%nhcopt==mol_setup(mnum(jj))%nhcopt) THEN
                    ieq = jj
                    EXIT
                 END IF
              END DO

! if the user defined option found a match, map everything
!    on the same thermostat
              IF (ieq /= 0 ) THEN

! map the thermostat number to the molecule
                 DO k = mol_off(i) + 1, mol_off(i+1)
                    molecule(k) %nvt_num = number
                 END DO
                 DO ii = 1, 3
                    DO jj = point(1,i), point(2,i)
                       nhc%p_kin(ii,jj) %point => nhc%s_kin(tnum(mnum(ieq)))
                       nhc%p_scale(ii,jj) %point => nhc%v_scale(tnum(mnum( &
                            ieq)))
                    END DO
                 END DO
                 deg_of_freedom(tnum(mnum(ieq))) = deg_of_freedom(tnum(mnum( &
                      ieq))) + tot_const(i)

! if the user defined option does not find a match.  Map it to a new
! thermostat, no questions asked
              ELSE
                 number = number + 1

! map the thermostat number to the molecule
                 DO k = mol_off(i) + 1, mol_off(i+1)
                    molecule(k) %nvt_num = number
                 END DO
                 mnum(icust) = i
                 tnum(icust) = number
                 DO ii = 1, 3
                    DO jj = point(1,i), point(2,i)
                       nhc%p_kin(ii,jj) %point => nhc%s_kin(number)
                       nhc%p_scale(ii,jj) %point => nhc%v_scale(number)
                    END DO
                 END DO
                 deg_of_freedom(number) = deg_of_freedom(number) + &
                      tot_const(i)
              END IF
           END IF

! doing global thermostat option
        ELSE IF (mol_setup(i)%nhcopt=="GLOBAL") THEN

! map the thermostat number to the molecule
           DO k = mol_off(i) + 1, mol_off(i+1)
              molecule(k) %nvt_num = number
           END DO
           DO ii = 1, 3
              DO jj = point(1,i), point(2,i)
                 nhc%p_kin(ii,jj) %point => nhc%s_kin(1)
                 nhc%p_scale(ii,jj) %point => nhc%v_scale(1)
              END DO
           END DO
           deg_of_freedom(1) = deg_of_freedom(1) + tot_const(i)

! doing molecular thermostat option
        ELSE IF (mol_setup(i)%nhcopt=="MOLECULE") THEN

! map the thermostat number to the molecule
           DO k = mol_off(i) + 1, mol_off(i+1)
              molecule(k) %nvt_num = number
           END DO
           DO ii = point(1,i), point(2,i), natoms_mol(i)
              number = number + 1
              start = ii
              finish = ii + molecule(i) %natoms_mol - 1
              DO jj = 1, 3
                 DO kk = start, finish
                    nhc%p_kin(jj,kk) %point => nhc%s_kin(number)
                    nhc%p_scale(jj,kk) %point => nhc%v_scale(number)
                 END DO
              END DO
              deg_of_freedom(number) = const_mol(i)
           END DO

! doing massive thermostat option
        ELSE IF ( mol_setup ( i ) % nhcopt == "MASSIVE" ) THEN
           IF ( tot_const ( i ) > 0 ) THEN
              CALL stop_program ( 'nhc_particle_setup', &
                   "Massive thermostats with constraints not possible" )
           END IF
           DO ii = point(1,i), point(2,i)
              DO jj = 1, 3
                 number = number + 1
                 nhc%p_kin(jj,ii) %point => nhc%s_kin(number)
                 nhc%p_scale(jj,ii) %point => nhc%v_scale(number)
              END DO
           END DO

        ELSE IF (mol_setup(i)%nhcopt=="NONE") THEN
           noneflag = .TRUE.
! NOTE: the thermostat counter is not incremented here.
! vscale(i,j) points to scale=1.d0 for eternity.  Thus
! the the velocites of particles that are not connected
! to a thermostat will not get scaled in lnhc.f
           DO ii = 1, 3
              DO jj = point(1,i), point(2,i)
                 nhc%p_scale(ii,jj) %point => scale
              END DO
           END DO
        END IF
     END DO
    END IF

! deallocating the temporary arrays mnum and tnum
    DEALLOCATE (mnum,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','mnum' )
    DEALLOCATE (tnum,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','tnum' )
    DEALLOCATE (const_mol,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol' )
    DEALLOCATE (tot_const,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const' )
    DEALLOCATE (natoms_mol,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'natoms_mol' )
    DEALLOCATE (mol_off,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'mol_off' )

    ALLOCATE (nhc%nvt(nhc%nhc_len,number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','nhc%nvt', 0 )

! Sum up the number of degrees of freedom on each thermostat.
! first: initialize the target
    natoms = sum ( molecule ( : ) % natoms_mol )
    nhc%s_kin = 0._dbl
    DO i = 1, 3
       DO j = 1, natoms
          nhc % p_kin ( i, j ) % point = nhc % p_kin ( i, j ) % point + 1
       END DO
    END DO

#if defined(__parallel)
! if thermostats are replicated but molecules distributed, we have to
!    sum s_kin over all processors
    IF ( nhc % dis_type == 'INTER_REP' ) THEN
       CALL mp_sum ( nhc % s_kin, globenv % group )
    END IF
#endif

!
! Now that we know how many there are stick this into nhc%static%nkt
! (number of degrees of freedom times k_B T for the first thermostat
!  on the chain)

    nhc % num_nhc = number

! We know the total number of system thermostats. For nvt it is number

    IF ( ( number == 1 ) .AND. ( .NOT. noneflag ) .AND. &
       ( .NOT. nhc % dis_type == 'DIS_REP' ) ) THEN
     nhc % nvt ( 1, 1 ) % nkt = simpar % temp_ext * ( nhc % s_kin ( 1 ) - &
       deg_of_freedom ( 1 ) - 3 )
    ELSE

! re-initializing simpar%nfree to zero because of multiple thermostats
     simpar % nfree = 0._dbl

     DO i = 1, number
        nhc % nvt ( 1, i ) % nkt = simpar % temp_ext * ( nhc % s_kin ( i ) - &
        deg_of_freedom ( i ) )
     END DO
    END IF

! getting the number of degrees of freedom times k_B T for the rest
! of the chain

    DO i = 2, nhc % nhc_len
       nhc % nvt ( i, : ) % nkt = simpar % temp_ext
    END DO

! Set up the Yoshida weights
    ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'nhc%dt_yosh', nhc % nyosh )
    CALL yoshida_coef ( nhc, simpar % dt )

  END SELECT


END SUBROUTINE nhc_particle_setup

!******************************************************************************
!!****** nhc/nhc_coef_setup [1.0] *
!!
!!   NAME
!!     nhc_coef_setup
!!
!!   FUNCTION
!!     Creates the thermostatting for polarization degrees of freedom
!!
!!   AUTHOR
!!     CJM 10-Apr-01
!!
!!   MODIFICATION HISTORY
!!     NONE  
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE nhc_coef_setup ( simpar, coeff, nhc, globenv )
! creates the thermostatting maps use global_mod
  IMPLICIT NONE
! begin global declaratation
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
  TYPE ( lnhc_parameters_type ), INTENT ( OUT ) :: nhc
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
! begin local definitions
  INTEGER :: i, j, iw
  INTEGER :: isos, number, ncoef
!
  iw = globenv % scr
!
  SELECT CASE ( simpar % ensemble_coef ( 1 : 3 ) )
  CASE  DEFAULT
    ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_coef_setup', &
       'failed to allocate nhc%nvt' )
    nhc % nyosh = 0
    nhc % nc = 0
    nhc % nhc_len = 0
    nhc % num_nhc = 0
  CASE ( 'NVT', 'NPT' )
    nhc % nyosh = simpar % nyosh_coef
    nhc % nc = simpar % nc_coef
    nhc % nhc_len = simpar % nhclen_coef

    SELECT CASE ( coeff % in_use )
    CASE ( ATOMICORBITALS )
      ncoef = SIZE ( coeff % ao % cr )
      IF ( simpar % nhcopt_coef == 'GLOBAL' )  THEN
        nhc % num_nhc = 1
      ELSEIF ( simpar % nhcopt_coef == 'MASSIVE' )  THEN
        nhc % num_nhc = ncoef
      ENDIF
    CASE ( PLANEWAVES )
      ncoef = SIZE ( coeff % pw % cc )
      IF ( simpar % nhcopt_coef == 'GLOBAL' )  THEN
        nhc % num_nhc = 1
      ELSEIF ( simpar % nhcopt_coef == 'MASSIVE' )  THEN
        nhc % num_nhc = 2 * ncoef
      ENDIF
    END SELECT    

! Now we can allocate the target array s_kin and p_kin

    ALLOCATE ( nhc % s_kin ( nhc % num_nhc ), STAT = isos)
    IF ( isos /= 0 ) CALL stop_program ( 'nhc_coef_setup', &
      'failed to allocate nhc%s_kin' )
    ALLOCATE ( nhc % v_scale ( nhc % num_nhc ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_coef_setup', &
      'failed to allocate nhc%v_scale' )
    IF ( coeff % in_use == ATOMICORBITALS ) THEN
      ALLOCATE ( nhc % p_kin ( 1, ncoef ), STAT = isos )
    ELSEIF ( coeff % in_use == PLANEWAVES ) THEN
      ALLOCATE ( nhc % p_kin ( 2, ncoef ), STAT = isos )
    ENDIF
    IF (isos/=0) CALL stop_program('nhc_coef_setup', &
      'failed to allocate nhc % p_kin' )
    IF ( coeff % in_use == ATOMICORBITALS ) THEN
      ALLOCATE ( nhc % p_scale ( 1, ncoef ), STAT = isos )
    ELSEIF ( coeff % in_use == PLANEWAVES ) THEN
      ALLOCATE ( nhc % p_scale ( 2, ncoef ), STAT = isos )
    ENDIF
    IF (isos/=0) CALL stop_program('nhc_coef_setup', &
      'failed to allocate nhc%p_scale')

! Begin the mapping loop

    IF ( simpar % nhcopt_coef == 'GLOBAL' )  THEN
      DO i = 1, ncoef
        DO j = 1, SIZE ( nhc % p_kin, 1 )
          nhc % p_kin ( j, i ) % point => nhc % s_kin ( 1 )
          nhc % p_scale ( j, i ) % point => nhc % v_scale ( 1 )
        END DO
      END DO
    ELSEIF ( simpar % nhcopt_coef == 'MASSIVE' )  THEN
      number = 0
      DO i = 1, ncoef
        DO j = 1, SIZE ( nhc % p_kin, 1 )
          number = number + 1
          nhc % p_kin ( j, i ) % point => nhc % s_kin ( number )
          nhc % p_scale ( j, i ) % point => nhc % v_scale ( number )
        END DO
      END DO
    ENDIF


    ALLOCATE ( nhc % nvt ( nhc % nhc_len, nhc % num_nhc ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_ao_setup', &
       'failed to allocate nhc%nvt' )

! Find the number of degrees of freedom on each thermostat
    nhc % s_kin = 0._dbl
    DO i = 1, ncoef
      DO j = 1, SIZE ( nhc % p_kin, 1 ) 
        nhc % p_kin ( j, i ) % point  = nhc % p_kin ( j, i ) % point  + 1._dbl
      END DO
    ENDDO
     
!
! Now that we know how many there are stick this into nhc % static % nkt
! (number of degrees of freedom times k_B T )

    DO i = 1, nhc % num_nhc
     nhc % nvt ( 1, i ) % nkt = simpar % temp_coef * nhc % s_kin ( i )
    END DO

! getting the number of degrees of freedom times k_B T for the rest of the chain

    DO i = 2, nhc % nhc_len
      nhc % nvt ( i, : ) % nkt = simpar % temp_coef
    END DO

! Set up the Yoshida weights
   
    ALLOCATE (nhc % dt_yosh (1 : nhc % nyosh ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_coef_setup', &
      'failed to allocate nhc%dt_yosh')
    CALL yoshida_coef( nhc, simpar % dt )
  END SELECT

END SUBROUTINE nhc_coef_setup

!******************************************************************************

END MODULE nhc

!******************************************************************************
