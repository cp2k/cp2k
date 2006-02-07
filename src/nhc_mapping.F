!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/nhc_mapping [1.0] *
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

MODULE nhc_mapping

  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             get_dyn_coeff
  USE extended_system_types,           ONLY: lnhc_parameters_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_nhcopt_global, do_nhcopt_massive, do_nhcopt_molecule, &
       isokin_ensemble, langevin_ensemble, nhc_global, nhc_massive, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_ensemble, reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nhc_to_particle_mapping, nhc_to_barostat_mapping, nhc_to_coef_mapping

!!***
!******************************************************************************
CONTAINS
!******************************************************************************
!!****** nhc/nhc_to_barostat_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_barostat_mapping
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

SUBROUTINE nhc_to_barostat_mapping ( simpar, nhc, globenv )
! creates the thermostatting maps use global_mod
! begin global declaratation
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER                                  :: i, isos, iw, number

! begin local definitions
!

  iw = globenv % scr
!
  nhc % dis_type = 'NONE'
  IF ( globenv % num_pe > 1 ) nhc % dis_type = 'REP_REP'
  SELECT CASE (simpar%ensemble)
  CASE  DEFAULT
    ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
    IF (isos/=0) CALL stop_program('nhc_barostat_setup', &
       'failed to allocate nhc%nvt' )
    nhc % nyosh = 0
    nhc % nc = 0
    nhc % nhc_len = 0
    nhc % num_nhc = 0
  CASE ( npt_i_ensemble, npt_f_ensemble )
    nhc % nyosh = simpar % nyosh
    nhc % nc = simpar % nc
    nhc % nhc_len = simpar % nhclen

! Counting the total number of thermostats ( 1 for both NPT_I and
! NPT_F )

    nhc % num_nhc = 1
! allocate the index array
    ALLOCATE ( nhc % index ( 1 ), STAT = isos ) 
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_barostat_setup', &
      'index', number )
   
    IF ( simpar % ensemble == npt_f_ensemble ) THEN
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
    nhc % index ( 1 ) = 1

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


  END SELECT

END SUBROUTINE nhc_to_barostat_mapping
!******************************************************************************
!!****** nhc/nhc_to_particle_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_particle_mapping
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

SUBROUTINE nhc_to_particle_mapping ( simpar, local_molecules, molecule_set, &
                                     molecule_kind_set, nhc, globenv )


    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER :: atm_offset, atom_increment, finish, first_atom, i, icount, &
      ielement, ii, iii, ikind, imol, imol_local, ipart, isos, iw, j, jj, kk, &
      last_atom, ll, n3x3con, n4x6con, natom, natoms_local, ndcon, nhcopt, &
      nkind, nmol_local, nmol_per_kind, nmolecule, number, number_total, &
      start, sum_of_thermostats
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: const_mol, deg_of_freedom, &
                                                tot_const
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: point
    INTEGER, DIMENSION(:), POINTER           :: massive_atom_list
    LOGICAL                                  :: have_global, nointer, noneflag
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  iw = globenv % scr
  NULLIFY ( massive_atom_list ) 

  SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
      CALL stop_program ( 'nhc_particle_setup', 'unknown ensemble ' )
    CASE ( nve_ensemble,&
           isokin_ensemble,&
           nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble,&
           reftraj_ensemble,&
           langevin_ensemble )
     ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%nvt', 0 )
     nhc % nyosh = 0
     nhc % nc = 0
     nhc % nhc_len = 0
     nhc % num_nhc = 0
     nhc % dis_type = ''
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble )
     nhc%nyosh = simpar%nyosh
     nhc%nc = simpar%nc
     nhc%nhc_len = simpar%nhclen

     nkind = SIZE ( molecule_kind_set )

! Counting the global number of thermostats
     sum_of_thermostats = 0
     noneflag = .FALSE.
! variable to denote  a global thermostat
     have_global = .FALSE.
! variable to denote independent thermostats (no communication necessary)
     nointer = .TRUE.
    

! 1) count the 'global' to determine distribution scheme
    DO i = 1, nkind
       molecule_kind => molecule_kind_set ( i )
       CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt )
       IF ( ( nhcopt == do_nhcopt_global ) ) THEN
          have_global = .TRUE.
          nointer = .FALSE.
          sum_of_thermostats = 1
       END IF
    END DO

! 2) count the 'molecule' 
    DO i = 1, nkind
       molecule_kind => molecule_kind_set ( i )
       CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt, nmolecule=nmolecule )
       IF ( ( nhcopt == do_nhcopt_molecule ) ) THEN
          sum_of_thermostats = sum_of_thermostats + nmolecule
! If we have ONE kind and ONE molecule, then effectively we have a GLOBAL thermostat
! and the degrees of freedom will be computed correctly for this special case
          IF ( ( nmolecule == 1 ) .AND. ( nkind == 1 ) ) nointer = .FALSE.
       END IF
    END DO

! 3) count the 'massive' 
    DO i = 1, nkind
       molecule_kind => molecule_kind_set ( i )
       CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt, nmolecule=nmolecule, &
                                natom=natom )
       IF ( ( nhcopt == do_nhcopt_massive) ) THEN
          sum_of_thermostats = sum_of_thermostats + 3 * natom * nmolecule
       END IF
    END DO
    number = sum_of_thermostats
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
        nhc % dis_type = 'DIS_REP' ! DIStributed thermostats, no interaction
! we only count thermostats on this processor
        number = 0
        DO ikind = 1, nkind
           nmol_local = local_molecules % n_el ( ikind ) 
           molecule_kind => molecule_kind_set ( ikind )
           CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt, natom=natom )
           IF (nhcopt==do_nhcopt_molecule) THEN
             number = number + nmol_local
           ELSE IF (nhcopt==do_nhcopt_massive) THEN
             number = number + 3*nmol_local*natom
           ELSE
             CALL stop_program ( 'nhc_particle_setup', '' )
           END IF
        END DO
      ELSE
        nhc % dis_type = 'INTER_REP' ! REPlicated thermostats, INTERacting via communication
      END IF
    END IF
! Allocate
    ALLOCATE ( nhc % index ( number ), STAT = isos ) 
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
      'index', number )

! compute the TOTAL number of molecules and atoms on THIS PROC and
! TOTAL number of molecules of IKIND on THIS PROC
    
    natoms_local = 0
    nmol_local = 0
    DO ikind = 1, nkind
      molecule_kind => molecule_kind_set ( ikind )
      CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt, natom=natom )
      natoms_local = natoms_local + natom * local_molecules % n_el ( ikind )
      nmol_local = nmol_local + local_molecules % n_el ( ikind )
    END DO
    IF ( nhc % dis_type == 'NONE' .OR. nhc % dis_type == 'DIS_REP' ) THEN
      ALLOCATE ( const_mol ( nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol', nmol_local) 
      ALLOCATE ( tot_const ( nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const', nmol_local) 
      ALLOCATE ( point ( 2, nmol_local ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'point', nmol_local*2) 
      point ( :, : ) = 0
      atm_offset = 0
      icount = 0
      DO ikind = 1, nkind
         nmol_per_kind =  local_molecules % n_el ( ikind )
         molecule_kind => molecule_kind_set ( ikind )
         CALL get_molecule_kind ( molecule_kind, ndist = ndcon, &
                                  ng3x3 = n3x3con, ng4x6 = n4x6con, &
                                  natom = natom )
         DO imol = 1, nmol_per_kind
           icount = icount + 1
           point ( 1, icount ) = atm_offset + 1 
           point ( 2, icount ) = atm_offset + natom 
           const_mol ( icount ) =  ndcon + 3 *  n3x3con + 6 * n4x6con
           tot_const ( icount ) = const_mol ( icount )
           atm_offset =  point ( 2, icount )
         END DO
      END DO

    ELSE IF ( nhc % dis_type == 'INTER_REP' ) THEN
      ALLOCATE ( const_mol ( nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol', nkind ) 
      ALLOCATE ( tot_const ( nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const', nkind ) 
      ALLOCATE ( point ( 2, nkind ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'point', nkind*2) 
      point ( :, : ) = 0
      atm_offset = 0
      DO ikind = 1, nkind
        nmol_per_kind =  local_molecules % n_el ( ikind )
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, natom=natom, ndist = ndcon, &
                                 ng3x3 = n3x3con, ng4x6 = n4x6con, nmolecule=nmolecule )
        const_mol ( ikind ) =  ndcon + 3 * n3x3con + 6 * n4x6con
        tot_const ( ikind ) = const_mol ( ikind ) * nmolecule
        point ( 1, ikind ) = atm_offset + 1 
        point ( 2, ikind ) = atm_offset + natom * nmol_per_kind
        atm_offset =  point ( 2, ikind )
      END DO
    ENDIF

    IF ( .NOT. simpar % constraint ) THEN
      const_mol = 0.0_dp
      tot_const = 0.0_dp
    END IF



! Now we can allocate the target array s_kin and p_kin
    ALLOCATE (deg_of_freedom(number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
          'deg_of_freedom', number )
    ALLOCATE (nhc%s_kin(number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%s_kin', number )
    ALLOCATE (nhc%v_scale(number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%v_scale', number )
    ALLOCATE (nhc%p_kin(3,natoms_local),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'nhc%p_kin', 3 * natoms_local )
    ALLOCATE (nhc%p_scale(3,natoms_local),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'nhc%p_scale', 3 * natoms_local )

! initialize deg_of_freedom
    deg_of_freedom = 0
    number = 0
    IF ( have_global ) number = 1

    CALL massive_list_generate ( molecule_set, molecule_kind_set, &
                                 local_molecules, globenv, massive_atom_list ) 
! DEAL WITH DIS & GLOBAL CASE
    icount = 0
    IF ( nhc % dis_type == 'NONE' .OR. nhc % dis_type == 'DIS_REP' ) THEN
      DO ikind = 1, nkind
        nmol_local =  local_molecules % n_el ( ikind ) 
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt )
        DO imol_local = 1, nmol_local
          IF ( nhcopt==do_nhcopt_global) THEN
             icount = icount + 1
             DO ii = 1, 3
                DO jj = point ( 1, icount ), point ( 2, icount )
                   nhc%p_kin ( ii, jj ) % point => nhc % s_kin ( 1 )
                   nhc%p_scale ( ii, jj ) % point => nhc % v_scale ( 1 )
                END DO
             END DO
             deg_of_freedom ( 1 ) = deg_of_freedom ( 1 ) + tot_const ( icount )
             nhc % index ( 1 ) = 1
          END IF
        END DO
      END DO
    END IF

! DEAL WITH DIS & MOLECULE CASE
    IF ( nhc % dis_type == 'NONE' .OR. nhc % dis_type == 'DIS_REP' ) THEN
      DO ikind = 1, nkind
        nmol_local =  local_molecules % n_el ( ikind ) 
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt )
        DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          IF ( nhcopt==do_nhcopt_molecule) THEN
            icount = icount + 1
            number = number + 1
            nhc % index ( number ) = imol
            deg_of_freedom ( number ) = const_mol ( icount )
            DO jj = 1, 3
              DO kk = point ( 1, icount ), point ( 2, icount )
                nhc%p_kin(jj,kk) %point => nhc%s_kin(number)
                nhc%p_scale(jj,kk) %point => nhc%v_scale(number)
              END DO
            END DO
          END IF
        END DO
      END DO
    END IF

! Communicate the total_number thermostat
    number_total = number
    CALL mp_sum ( number_total, globenv % group )

! DEAL WITH DIS & MASSIVE CASE
    IF ( nhc % dis_type == 'NONE' .OR. nhc % dis_type == 'DIS_REP' ) THEN
      DO ikind = 1, nkind
        nmol_local =  local_molecules % n_el ( ikind ) 
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt=nhcopt )
        DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom = first_atom, &
                              last_atom = last_atom ) 
          IF ( nhcopt==do_nhcopt_massive) THEN
            icount = icount + 1
            IF ( tot_const ( icount ) > 0 ) CALL stop_program ( 'nhc_particle_setup', &
                                 "Massive thermostats with constraints are impossible" )
            atom_increment = 0
            DO ii =  point ( 1, icount ), point ( 2, icount )
               ipart = first_atom + atom_increment
               DO  iii = 1, SIZE ( massive_atom_list )
                 IF ( ipart == massive_atom_list ( iii ) ) THEN
                   ielement = iii
                   EXIT
                 END IF
               END DO
               atom_increment = atom_increment + 1
               DO jj = 1, 3
                  number = number + 1
                  nhc % index ( number ) = number_total + ( ielement - 1 ) * 3 + jj
                  nhc%p_kin(jj,ii) %point => nhc%s_kin(number)
                  nhc%p_scale(jj,ii) %point => nhc%v_scale(number)
               END DO
            END DO
            IF ( first_atom + atom_increment -1 /= last_atom ) &
            CALL stop_program ( 'nhc_particle_setup',          &
            "inconsistent mapping of particles" )
          END IF
        END DO
      END DO
    END IF

! DEAL WITH INTER_REP AND GLOBAL
    IF ( nhc % dis_type == 'INTER_REP' ) THEN    
      DO ikind = 1, nkind
! doing global thermostat option
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt = nhcopt )
        IF ( nhcopt==do_nhcopt_global) THEN
! map the thermostat number to the molecule
           DO ii = 1, 3
              DO jj = point ( 1, ikind ), point ( 2, ikind )
                 nhc%p_kin ( ii, jj ) % point => nhc % s_kin ( 1 )
                 nhc%p_scale ( ii, jj ) % point => nhc % v_scale ( 1 )
              END DO
           END DO
           deg_of_freedom ( 1 ) = deg_of_freedom ( 1 ) + tot_const ( ikind )
           nhc % index ( 1 ) = 1
        END IF
      END DO
    ENDIF 

! DEAL WITH INTER_REP AND MOLECULE
    IF ( nhc % dis_type == 'INTER_REP' ) THEN    
      DO ikind = 1, nkind
! doing global thermostat option
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt = nhcopt, natom = natom )
! doing molecular thermostat option
        IF ( nhcopt==do_nhcopt_molecule ) THEN
! map the thermostat number to the molecule
          number = number + 1
          nhc % index ( number ) = number
          deg_of_freedom ( number ) = deg_of_freedom ( number ) + tot_const ( ikind )
          DO jj = 1, 3
            DO kk = point ( 1, ikind ), point ( 2, ikind ), natom
              start = kk
              finish = kk + ( natom - 1 )
              DO ll = start, finish
                nhc%p_kin ( jj, ll ) % point => nhc % s_kin ( number )
                nhc%p_scale ( jj, ll ) % point => nhc % v_scale ( number )
              END DO
            END DO
          END DO
        ENDIF
      END DO
    END IF

! DEAL WITH INTER_REP AND MASSIVE
    IF ( nhc % dis_type == 'INTER_REP' ) THEN    
      DO ikind = 1, nkind
! doing global thermostat option
        molecule_kind => molecule_kind_set ( ikind )
        CALL get_molecule_kind ( molecule_kind, nhcopt = nhcopt )
! doing  massive thermostat option
        IF ( nhcopt==do_nhcopt_massive) THEN
          IF ( tot_const ( ikind ) > 0 ) CALL stop_program ( 'nhc_particle_setup', &
                              "Massive thermostats with constraints are impossible" )
          DO ii = point ( 1, ikind ), point ( 2, ikind )
            DO jj = 1, 3
              number = number + 1
              nhc % index ( number ) = number
              nhc%p_kin ( jj, ii ) % point => nhc % s_kin ( number )
              nhc%p_scale ( jj, ii ) % point => nhc % v_scale ( number )
            END DO
          END DO
        END IF
      END DO
    ENDIF


    DEALLOCATE (const_mol,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'const_mol' )
    DEALLOCATE (tot_const,STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'tot_const' )
    DEALLOCATE ( point, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'point' )

    ALLOCATE (nhc%nvt(nhc%nhc_len,number),STAT=isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup','nhc%nvt', 0 )

! Sum up the number of degrees of freedom on each thermostat.
! first: initialize the target


    nhc%s_kin = 0.0_dp
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

! We know the total number of system thermostats. 

    IF ( ( sum_of_thermostats == 1 ) .AND. ( .NOT. noneflag ) .AND. &
       ( .NOT. nhc % dis_type == 'DIS_REP' ) ) THEN
     nhc % nvt ( 1, 1 ) % nkt = simpar % temp_ext * ( nhc % s_kin ( 1 ) - &
       deg_of_freedom ( 1 ) - simpar % nfree )
    ELSE

! re-initializing simpar%nfree to zero because of multiple thermostats
     simpar % nfree = 0.0_dp

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

  END SELECT
  
  IF ( ASSOCIATED ( massive_atom_list ) ) THEN 
    DEALLOCATE ( massive_atom_list, STAT = isos ) 
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', 'massive_atom_list' )
  END IF

END SUBROUTINE nhc_to_particle_mapping

!******************************************************************************
!!****** nhc/nhc_to_coef_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_coef_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting for the
!!     polarization  coefficients'  degrees of freedom
!!
!!   AUTHOR
!!     CJM 10-Apr-01
!!
!!   MODIFICATION HISTORY
!!     gt 04-02-04 adapted to new kind of coefficents
!!     gt 11-02-05 adapted to new input                
!!
!!   SOURCE
!******************************************************************************
 SUBROUTINE nhc_to_coef_mapping ( dyn_coeff_set, nhc, globenv )

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER                                  :: i, ikind, isos, iw, j, &
                                                kg_coeff_ens, kg_coeff_nhc, &
                                                ncoef_local, ncoef_of_kind, &
                                                nkind, number
    LOGICAL                                  :: nointer
    REAL(KIND=dp)                            :: kg_coeff_temp
    TYPE(section_vals_type), POINTER         :: kg_coeff_section, nose_section

!------------------------------------------------------------------------------

  iw = globenv % scr
  
  kg_coeff_section => section_vals_get_subs_vals(globenv%input_file,&
       "MOTION%MD")
  kg_coeff_section => section_vals_get_subs_vals(kg_coeff_section,&
       "kg_coeff")
  nose_section => section_vals_get_subs_vals(kg_coeff_section,&
       "nose")

  CALL section_vals_val_get(kg_coeff_section,"ensemble",i_val=kg_coeff_ens)
  CALL section_vals_val_get(kg_coeff_section,"nhcopt",i_val=kg_coeff_nhc)
  CALL section_vals_val_get(kg_coeff_section,"temperature",r_val=kg_coeff_temp)
  
  SELECT CASE ( kg_coeff_ens )
  CASE  DEFAULT
    ALLOCATE ( nhc % nvt ( 0, 0 ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', &
       'allocate nhc%nvt',0 )
    nhc % nyosh = 0
    nhc % nc = 0
    nhc % nhc_len = 0
    nhc % num_nhc = 0
    nhc % dis_type = ''
  CASE ( nvt_ensemble )
    CALL section_vals_val_get(nose_section,"yoshida",i_val=nhc % nyosh )
    CALL section_vals_val_get(nose_section,"mts",i_val=nhc % nc )
    CALL section_vals_val_get(nose_section,"length",i_val=nhc % nhc_len )
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

    IF ( kg_coeff_nhc == nhc_global )  THEN
      nointer=.FALSE.
      number=1
    ELSEIF ( kg_coeff_nhc == nhc_massive )  THEN
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

    ALLOCATE ( nhc % index ( number ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_coef_setup', &
      'index', number )
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
    IF ( kg_coeff_nhc == nhc_global )  THEN
      DO i= 1, ncoef_local
        DO j = 1, SIZE ( nhc % p_kin, 1 )
          nhc % p_kin ( j, i ) % point => nhc % s_kin ( 1 )
          nhc % p_scale ( j, i ) % point => nhc % v_scale ( 1 )
        END DO
      END DO
      nhc % index ( 1 ) = 1
    ELSEIF (  kg_coeff_nhc == nhc_massive )  THEN
      DO i = 1, ncoef_local 
        DO j = 1, SIZE ( nhc % p_kin, 1 )
          number = number + 1
          nhc % p_kin ( j, i ) % point => nhc % s_kin ( number )
          nhc % p_scale ( j, i ) % point => nhc % v_scale ( number )
          nhc % index ( number ) = number
        END DO
      END DO
    ENDIF

    ALLOCATE ( nhc % nvt ( nhc % nhc_len, nhc % num_nhc ), STAT = isos )
    IF (isos/=0) CALL stop_memory('nhc_coef_setup', 'nhc%nvt', nhc%nhc_len*nhc%num_nhc)

! Find the number of degrees of freedom on each thermostat
    nhc % s_kin = 0.0_dp
    DO i = 1, ncoef_local
      DO j = 1, SIZE ( nhc % p_kin, 1 ) 
        nhc % p_kin ( j, i ) % point  = nhc % p_kin ( j, i ) % point  + 1.0_dp
      END DO
    ENDDO

! if thermostats are replicated but molecules distributed, we have to
! sum s_kin over all processors
    IF ( nhc % dis_type == 'INTER_REP' ) CALL mp_sum ( nhc % s_kin, globenv % group )

!
! Now that we know how many there are stick this into nhc % static % nkt
! (number of degrees of freedom times k_B T )

    DO i = 1, nhc % num_nhc
     nhc % nvt ( 1, i ) % nkt = kg_coeff_temp * nhc % s_kin ( i )
    END DO

! getting the number of degrees of freedom times k_B T for the rest of the chain

    DO i = 2, nhc % nhc_len
      nhc % nvt ( i, : ) % nkt = kg_coeff_temp
    END DO
  END SELECT

 END SUBROUTINE nhc_to_coef_mapping

!******************************************************************************
 SUBROUTINE massive_list_generate ( molecule_set, molecule_kind_set, &
                                    local_molecules, globenv, massive_atom_list ) 
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, POINTER                         :: massive_atom_list( : )

    INTEGER :: first_atom, i, ikind, imol, iproc, isos, j, k, natom, ncount, &
      nhcopt, nkind, nmol_per_kind, num_massive_atm, num_massive_atm_local, &
      offset
    INTEGER, DIMENSION(:), POINTER           :: array_num_massive_atm, &
                                                local_atm_list, work
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

  num_massive_atm_local = 0
  NULLIFY(local_atm_list)
  CALL reallocate(local_atm_list,1,num_massive_atm_local)

  nkind = SIZE ( molecule_kind_set )
  DO ikind = 1,  nkind
    nmol_per_kind =  local_molecules%n_el(ikind)
    DO imol = 1, nmol_per_kind
      i = local_molecules%list(ikind)%array(imol)
      molecule => molecule_set ( i )
      molecule_kind => molecule % molecule_kind
      CALL get_molecule_kind(molecule_kind,natom=natom,nhcopt=nhcopt)
      IF(nhcopt == do_nhcopt_massive) THEN
        num_massive_atm_local = num_massive_atm_local + natom
        CALL reallocate(local_atm_list,1,num_massive_atm_local)
        CALL get_molecule (molecule,first_atom=first_atom)
        DO j=1,natom
          local_atm_list(num_massive_atm_local-natom+j) = first_atom -1 + j
        END DO
      END IF
    END DO
  END DO

  ALLOCATE(array_num_massive_atm(globenv%num_pe), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'array_num_massive_atm', globenv%num_pe )
  CALL mp_allgather(num_massive_atm_local,array_num_massive_atm,globenv%group)

  num_massive_atm = SUM(array_num_massive_atm)
  ALLOCATE(massive_atom_list(num_massive_atm), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'massive_atom_list', num_massive_atm )

  offset = 0
  DO iproc=1,globenv%num_pe
    ncount = array_num_massive_atm(iproc)
    ALLOCATE(work(ncount), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'work', ncount )
    IF(globenv%mepos == (iproc-1)) THEN
      DO i=1,ncount
        work(i) = local_atm_list(i)
      END DO
    ELSE
      work(:) = 0
    END IF
    CALL mp_bcast(work,iproc-1,globenv%group)
    DO i=1,ncount
      massive_atom_list(offset+i) = work(i)
    END DO
    DEALLOCATE(work, STAT = isos)
    IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'work' )
    offset = offset + array_num_massive_atm(iproc)
  END DO

  DO i=1,num_massive_atm
    DO j=1,num_massive_atm
      IF(massive_atom_list(i)<massive_atom_list(j)) THEN
        k = massive_atom_list ( i )
        massive_atom_list ( i )  = massive_atom_list ( j )
        massive_atom_list ( j ) = k
      ENDIF
    END DO
  END DO

  DEALLOCATE(local_atm_list, STAT = isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
      'local_atm_list' )
  DEALLOCATE(array_num_massive_atm, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'nhc_particle_setup', &
       'array_num_massive_atm' )

  END SUBROUTINE massive_list_generate

!******************************************************************************

END MODULE nhc_mapping
