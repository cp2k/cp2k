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

  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             get_dyn_coeff
  USE extended_system_initialize,      ONLY: set_yoshida_coef
  USE extended_system_types,           ONLY: lnhc_parameters_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

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
! begin global declaratation
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(lnhc_parameters_type), INTENT(OUT)  :: nhc
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, isos, iw, number

! begin local definitions
!

  iw = globenv % scr
!
  nhc % dis_type = 'NONE'
  IF ( globenv % num_pe > 1 ) nhc % dis_type = 'DIS_REP'
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
    nhc % dis_type = 'NONE'
   
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
    CALL set_yoshida_coef( nhc, simpar % dt )
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
!!     CJM(10-NOV-2001) : New parallelization with new molecule structures
!!
!!*** *************************************************************************

SUBROUTINE nhc_particle_setup ( simpar, local_molecules, molecule, &
                                molecule_kind_set, nhc, globenv )


    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(lnhc_parameters_type), INTENT(OUT)  :: nhc
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=default_string_length)     :: nhcopt
    INTEGER :: atm_offset, finish, i, icount, ii, ikind, imol, isos, iw, j, &
      jj, kk, n3x3con, natom, natoms_local, ndcon, nkind, nmol_local, &
      nmol_per_kind, nmolecule, number, start
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: const_mol, deg_of_freedom, &
                                                natom_mol, &
                                                nmol_local_of_kind, tot_const
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: point
    LOGICAL                                  :: have_global, nointer, noneflag
    REAL(dbl), SAVE, TARGET                  :: scale
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

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

     nkind = SIZE ( molecule_kind_set )

! Counting the global number of thermostats
     number = 0
     noneflag = .FALSE.
     have_global = .FALSE.
     nointer = .TRUE.
    

! 1) count the 'global' to determine distribution scheme
    DO i = 1, nkind
       molecule_kind => molecule_kind_set ( i )
       CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt )
       IF ( ( nhcopt == "GLOBAL" ) ) THEN
          have_global = .TRUE.
          nointer = .FALSE.
          number = number + 1
       END IF
    END DO

! 2) count the 'molecule' 
    DO i = 1, nkind
       molecule_kind => molecule_kind_set ( i )
       CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt, nmolecule=nmolecule )
       IF ( ( nhcopt == "MOLECULE" ) ) THEN
          number = number + nmolecule
       END IF
    END DO

! 3) count the 'massive' 
    DO i = 1, nkind
       molecule_kind => molecule_kind_set ( i )
       CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt, natom=natom )
       IF ( ( nhcopt == "MASSIVE" ) ) THEN
          number = number + 3 * natom
       END IF
    END DO
!
! here we decide on what parallel algorithm to use
! if there are only massive and molecule type thermostats we can use
! a local scheme, in cases involving any combination with a
! global thermostat we assume a coupling of  degrees of freedom 
! from different processors
!
    nhc % dis_type = 'NONE'
    IF ( globenv % num_pe > 1 ) THEN
      IF ( nointer ) THEN
        nhc % dis_type = 'DIS_REP'
! we only count thermostats on this processor
        number = 0
        DO ikind = 1, nkind
           nmol_local = local_molecules % n_el ( ikind ) 
           DO imol = 1, nmol_local
             i = local_molecules % list ( ikind ) % array ( imol )
             molecule_kind => molecule ( i ) % molecule_kind
             CALL get_molecule_kind ( molecule_kind,  nhcopt=nhcopt, natom=natom )
             IF (nhcopt=="MOLECULE") THEN
                number = number + 1
             ELSE IF (nhcopt=="MASSIVE") THEN
                number = number + 3*natom
             ELSE
                CALL stop_program ( 'nhc_particle_setup', '' )
             END IF
           END DO
        END DO
      ELSE
        nhc % dis_type = 'INTER_REP' 
      END IF
    END IF

! compute the TOTAL number of molecules and atoms on THIS PROC and
! TOTAL number of molecules of IKIND on THIS PROC
    
    ALLOCATE ( nmol_local_of_kind ( nkind ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nmol_local_of_kind', nkind )
    natoms_local = 0
    DO ikind = 1, nkind
      nmol_per_kind =  local_molecules % n_el ( ikind )
      nmol_local_of_kind ( ikind ) = nmol_per_kind
      DO imol = 1, nmol_per_kind
        i = local_molecules % list ( ikind ) % array ( imol )
        molecule_kind => molecule ( i ) % molecule_kind
        CALL get_molecule_kind ( molecule_kind, natom=natom)
        natoms_local = natoms_local + natom
      END DO
    END DO
    nmol_local = SUM ( nmol_local_of_kind )

    IF ( nhc % dis_type == 'NONE' .OR. nhc % dis_type == 'DIS_REP' ) THEN
      ALLOCATE ( const_mol ( nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol', nmol_local) 
      ALLOCATE ( tot_const ( nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const', nmol_local) 
      ALLOCATE ( natom_mol ( nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'natom_mol', nmol_local) 
      ALLOCATE ( point ( 2, nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'point', nmol_local*2) 
      point ( :, : ) = 0
      icount = 0
      atm_offset = 0
      DO ikind = 1, nkind
         nmol_per_kind =  local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           icount = icount + 1
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule_kind => molecule ( i ) % molecule_kind
           CALL get_molecule_kind ( molecule_kind, natom=natom, ndist=ndcon, ng3x3=n3x3con )
           natom_mol ( icount ) = natom
           point ( 1, icount ) = atm_offset + 1 
           point ( 2, icount ) = atm_offset + natom 
           const_mol ( icount ) =  ndcon + 3 *  n3x3con
           tot_const ( icount ) = const_mol ( icount )
           atm_offset =  point ( 2, icount )
         END DO
      END DO

    ELSE IF ( nhc % dis_type == 'INTER_REP' ) THEN
      ALLOCATE ( const_mol ( nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol', nkind ) 
      ALLOCATE ( tot_const ( nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const', nkind ) 
      ALLOCATE ( natom_mol ( nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'natom_mol', nkind ) 
      ALLOCATE ( point ( 2, nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'point', nkind*2) 
      point ( :, : ) = 0
      atm_offset = 0
      DO ikind = 1, nkind
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, natom=natom, ndist = ndcon, &
                                 ng3x3 = n3x3con, nmolecule = nmolecule )
        const_mol ( ikind ) =  ndcon + 3 * n3x3con
        tot_const ( ikind ) = const_mol ( ikind ) * nmolecule
        natom_mol ( ikind ) =  natom
        point ( 1, ikind ) = atm_offset + 1 
        point ( 2, ikind ) = atm_offset + natom * nmol_local_of_kind ( ikind )
        atm_offset =  point ( 2, ikind )
      END DO
    ENDIF

    IF ( .NOT. simpar % constraint ) THEN
      const_mol = 0._dbl
      tot_const = 0._dbl
    END IF



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
    ALLOCATE (nhc%p_kin(3,natoms_local),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%p_kin', 3 * natoms_local )
    ALLOCATE (nhc%p_scale(3,natoms_local),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'nhc%p_scale', 3 * natoms_local )

! initialize variable 'scale' for the 'none' option
    scale = 1.0_dbl

! initialize deg_of_freedom
    deg_of_freedom = 0
    number = 0
    IF ( have_global ) number = 1

!  Begin the mapping loop
    IF ( nhc % dis_type == 'NONE' .OR. nhc % dis_type == 'DIS_REP' ) THEN
      icount = 0
      DO ikind = 1, nkind
        nmol_local =  local_molecules % n_el ( ikind ) 
        DO imol = 1, nmol_local
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule_kind => molecule ( i ) % molecule_kind
          icount = icount + 1
          CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt )
! doing Global thermostat option
          IF ( nhcopt=="GLOBAL") THEN
! map the thermostat number to the molecule
             DO ii = 1, 3
                DO jj = point ( 1, icount ), point ( 2, icount )
                   nhc%p_kin ( ii, jj ) % point => nhc % s_kin ( 1 )
                   nhc%p_scale ( ii, jj ) % point => nhc % v_scale ( 1 )
                END DO
             END DO
             deg_of_freedom ( 1 ) = deg_of_freedom ( 1 ) + tot_const ( icount )

! doing molecular thermostat option
          ELSE IF ( nhcopt == "MOLECULE" ) THEN
             DO ii = point ( 1, icount ), point ( 2, icount ), natom_mol ( icount )
                number = number + 1
                start = ii
                finish = ii + natom_mol ( icount )
! map the thermostat number to the molecule
                DO jj = 1, 3
                  DO kk = start, finish
                    nhc%p_kin(jj,kk) %point => nhc%s_kin(number)
                    nhc%p_scale(jj,kk) %point => nhc%v_scale(number)
                  END DO
                END DO
                deg_of_freedom ( number ) = const_mol ( icount )
             END DO
! doing massive thermostat option
          ELSE IF ( nhcopt == "MASSIVE" ) THEN
             IF ( tot_const ( icount ) > 0 ) CALL stop_program ( 'nhc_particle_setup', &
                                 "Massive thermostats with constraints are impossible" )
             DO ii =  point ( 1, icount ), point ( 2, icount )
                DO jj = 1, 3
                   number = number + 1
                   nhc%p_kin(jj,ii) %point => nhc%s_kin(number)
                   nhc%p_scale(jj,ii) %point => nhc%v_scale(number)
                END DO
             END DO

          ELSE IF (nhcopt=="NONE") THEN
             noneflag = .TRUE.
! NOTE: the thermostat counter is not incremented here.
! nhc % p_scale(i,j) points to scale=1.d0 for eternity.  Thus
! the the velocites of particles that are not connected
! to a thermostat will not get scaled in lnhc.f
             DO ii = 1, 3
                DO jj = point ( 1, icount ), point ( 2, icount )
                   nhc%p_scale(ii,jj) %point => scale
                END DO
             END DO

          END IF
        END DO
      END DO
    ELSEIF ( nhc % dis_type == 'INTER_REP' ) THEN    
      DO ikind = 1, nkind
! doing global thermostat option
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt = nhcopt )
        IF ( nhcopt=="GLOBAL") THEN
! map the thermostat number to the molecule
           DO ii = 1, 3
              DO jj = point ( 1, ikind ), point ( 2, ikind )
                 nhc%p_kin ( ii, jj ) % point => nhc % s_kin ( 1 )
                 nhc%p_scale ( ii, jj ) % point => nhc % v_scale ( 1 )
              END DO
           END DO
           deg_of_freedom ( 1 ) = deg_of_freedom ( 1 ) + tot_const ( ikind )

! doing molecular thermostat option
        ELSEIF ( nhcopt=="MOLECULE") THEN
! map the thermostat number to the molecule
          DO ii = point ( 1, ikind ), point ( 2, ikind ), natom_mol ( ikind )
            number = number + 1
            start = ii
            finish = ii + natom_mol ( ikind )
            DO jj = 1, 3
              DO kk = start, finish
                nhc%p_kin ( jj, kk ) % point => nhc % s_kin ( number )
                nhc%p_scale ( jj, kk ) % point => nhc % v_scale ( number )
              END DO
            END DO
            deg_of_freedom ( number ) = deg_of_freedom ( number ) + tot_const ( ikind )
          END DO
! doing  massive thermostat option
        ELSEIF ( nhcopt=="MASSIVE") THEN
          IF ( tot_const ( ikind ) > 0 ) CALL stop_program ( 'nhc_particle_setup', &
                              "Massive thermostats with constraints are impossible" )
          DO ii = point ( 1, ikind ), point ( 2, ikind )
            DO jj = 1, 3
              number = number + 1
              nhc%p_kin ( jj, ii ) % point => nhc % s_kin ( number )
              nhc%p_scale ( jj, ii ) % point => nhc % v_scale ( number )
            END DO
          END DO
! doing none option
        ELSEIF ( nhcopt=="NONE") THEN
          noneflag = .TRUE.
          DO ii = 1, 3
            DO jj = point ( 1, ikind ), point ( 2, ikind )
              nhc%p_scale ( ii, jj ) % point => scale
            END DO
          END DO
        END IF
      END DO
    ENDIF

    DEALLOCATE (const_mol,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol' )
    DEALLOCATE (tot_const,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const' )
    DEALLOCATE ( natom_mol, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'natom_mol' )
    DEALLOCATE ( point, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'point' )
    DEALLOCATE ( nmol_local_of_kind, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nmol_local_of_kind' )

    ALLOCATE (nhc%nvt(nhc%nhc_len,number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','nhc%nvt', 0 )

! Sum up the number of degrees of freedom on each thermostat.
! first: initialize the target


    nhc%s_kin = 0._dbl
    DO i = 1, 3
       DO j = 1, natoms_local
          nhc % p_kin ( i, j ) % point = nhc % p_kin ( i, j ) % point + 1
       END DO
    END DO

! if thermostats are replicated but molecules distributed, we have to
! sum s_kin over all processors
    IF ( nhc % dis_type == 'INTER_REP' ) CALL mp_sum ( nhc % s_kin, globenv % group )

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
    CALL set_yoshida_coef ( nhc, simpar % dt )

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
!!     gt 04-02-04 adapted to new kind of coefficents
!!
!!   SOURCE
!******************************************************************************
 SUBROUTINE nhc_coef_setup ( simpar, dyn_coeff_set, nhc, globenv )

    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(lnhc_parameters_type), INTENT(OUT)  :: nhc
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, ikind, isos, iw, j, &
                                                ncoef_local, ncoef_of_kind, &
                                                nkind, number
    LOGICAL                                  :: nointer

!------------------------------------------------------------------------------

  iw = globenv % scr

  SELECT CASE ( simpar % ensemble_coef ( 1 : 3 ) )
  CASE  DEFAULT
    ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', &
       'allocate nhc%nvt',0 )
    nhc % nyosh = 0
    nhc % nc = 0
    nhc % nhc_len = 0
    nhc % num_nhc = 0
  CASE ( 'NVT', 'NPT' )
    nhc % nyosh = simpar % nyosh_coef
    nhc % nc = simpar % nc_coef
    nhc % nhc_len = simpar % nhclen_coef

    number = 0
    ncoef_local = 0
    nointer = .TRUE.

    nkind = SIZE ( dyn_coeff_set%coeffs_of_kind)
    DO ikind=1,nkind
      IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN  
        CALL get_dyn_coeff(coeffs=dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                           ncoef_of_kind=ncoef_of_kind)
         ncoef_local= ncoef_local + ncoef_of_kind 
       ENDIF
    END DO

! Counting the global number of thermostats

    IF ( simpar % nhcopt_coef == 'GLOBAL' )  THEN
      nointer=.FALSE.
      number=1
    ELSEIF ( simpar % nhcopt_coef == 'MASSIVE' )  THEN
      number=ncoef_local 
    ENDIF

! coefficients are local. communications needed only with GLOBAL
    nhc % dis_type = 'NONE'
    IF ( globenv % num_pe > 1 ) THEN
      IF ( nointer ) THEN
        nhc % dis_type = 'DIS_REP'
      ELSE
        nhc % dis_type = 'INTER_REP'
      END IF
    END IF

! Now we can allocate the target array s_kin and p_kin

    ALLOCATE ( nhc % s_kin ( number ), STAT = isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_coef_setup', &
      'nhc%s_kin', number  )
    ALLOCATE ( nhc % v_scale ( number  ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', &
      'nhc%v_scale', number )
    ALLOCATE ( nhc % p_kin ( 1, ncoef_local ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', &
      'nhc % p_kin', ncoef_local )
    ALLOCATE ( nhc % p_scale ( 1, ncoef_local ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', &
      'nhc%p_scale', ncoef_local)
    nhc % num_nhc = number

! Begin the mapping loop

    i=0
    number = 0
    IF ( simpar % nhcopt_coef == 'GLOBAL' )  THEN
      DO i= 1, ncoef_local
        DO j = 1, SIZE ( nhc % p_kin, 1 )
          nhc % p_kin ( j, i ) % point => nhc % s_kin ( 1 )
          nhc % p_scale ( j, i ) % point => nhc % v_scale ( 1 )
        END DO
      END DO
    ELSEIF ( simpar % nhcopt_coef == 'MASSIVE' )  THEN
      DO i = 1, ncoef_local 
        DO j = 1, SIZE ( nhc % p_kin, 1 )
          number = number + 1
          nhc % p_kin ( j, i ) % point => nhc % s_kin ( number )
          nhc % p_scale ( j, i ) % point => nhc % v_scale ( number )
        END DO
      END DO
    ENDIF

    ALLOCATE ( nhc % nvt ( nhc % nhc_len, nhc % num_nhc ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', 'nhc%nvt', nhc%nhc_len*nhc%num_nhc)

! Find the number of degrees of freedom on each thermostat
    nhc % s_kin = 0._dbl
    DO i = 1, ncoef_local
      DO j = 1, SIZE ( nhc % p_kin, 1 ) 
        nhc % p_kin ( j, i ) % point  = nhc % p_kin ( j, i ) % point  + 1._dbl
      END DO
    ENDDO

! if thermostats are replicated but molecules distributed, we have to
! sum s_kin over all processors
    IF ( nhc % dis_type == 'INTER_REP' ) CALL mp_sum ( nhc % s_kin, globenv % group )

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
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', &
      'nhc%dt_yosh',nhc % nyosh )
    CALL set_yoshida_coef( nhc, simpar % dt )
  END SELECT

 END SUBROUTINE nhc_coef_setup

!******************************************************************************

END MODULE nhc
