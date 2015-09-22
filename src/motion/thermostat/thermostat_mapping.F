!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! *****************************************************************************
MODULE thermostat_mapping

  USE cell_types,                      ONLY: use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: map_info_type
  USE input_constants,                 ONLY: do_region_defined,&
                                             do_region_global,&
                                             do_region_massive,&
                                             do_region_molecule,&
                                             do_thermo_communication,&
                                             do_thermo_no_communication
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate,&
                                             sort
#include "../../common/cp_common_uses.f90"

  IMPLICIT NONE

  PUBLIC :: thermostat_mapping_region,&
            adiabatic_mapping_region,&
            init_baro_map_info

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_mapping'

CONTAINS
! *****************************************************************************
!> \brief Main general setup for adiabatic thermostat regions (Nose only)
!> \param map_info ...
!> \param deg_of_freedom ...
!> \param massive_atom_list ...
!> \param molecule_kind_set ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param para_env ...
!> \param natoms_local ...
!> \param simpar ...
!> \param number ...
!> \param region ...
!> \param gci ...
!> \param shell ...
!> \param map_loc_thermo_gen ...
!> \param sum_of_thermostats ...
!> \author CJM - PNNL
! *****************************************************************************
  SUBROUTINE adiabatic_mapping_region ( map_info, deg_of_freedom, massive_atom_list,&
       molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, simpar, &
       number, region, gci, shell, map_loc_thermo_gen, sum_of_thermostats)

    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(OUT)                     :: natoms_local
    TYPE(simpar_type), POINTER               :: simpar
    INTEGER, INTENT(INOUT)                   :: number
    INTEGER, INTENT(IN)                      :: region
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: shell
    INTEGER, DIMENSION(:), POINTER           :: map_loc_thermo_gen
    INTEGER, INTENT(INOUT)                   :: sum_of_thermostats

    CHARACTER(LEN=*), PARAMETER :: routineN = 'adiabatic_mapping_region', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nkind, nmol_local, &
                                                nsize, number_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY ( const_mol, tot_const, point)
    CPASSERT(.NOT.ASSOCIATED(deg_of_freedom))
    CPASSERT(.NOT.ASSOCIATED(massive_atom_list))

    nkind = SIZE(molecule_kind_set)
    CALL adiabatic_region_evaluate(map_info%dis_type, natoms_local, nmol_local,&
         const_mol, tot_const, point, local_molecules, molecule_kind_set, molecule_set,&
         simpar, shell)

    ! Now we can allocate the target array s_kin and p_kin..
    SELECT CASE(region)
    CASE(do_region_global, do_region_molecule, do_region_massive)
       nsize = number
    CASE DEFAULT
       ! STOP PROGRAM
    END SELECT
    ALLOCATE (map_info%s_kin(nsize))
    ALLOCATE (map_info%v_scale(nsize))
    ALLOCATE (map_info%p_kin(3,natoms_local))
    ALLOCATE (map_info%p_scale(3,natoms_local))
! nullify thermostat pointers
    ! Allocate index array to 1
    ALLOCATE ( map_info%index(1))
    ALLOCATE ( map_info%map_index(1))
    ALLOCATE (deg_of_freedom(1))

    CALL massive_list_generate (molecule_set, molecule_kind_set, &
         local_molecules, para_env, massive_atom_list, region, shell)

    CALL adiabatic_mapping_region_low(region, map_info, nkind, point,&
         deg_of_freedom, local_molecules, const_mol, massive_atom_list,&
         tot_const, molecule_set, number_of_thermostats, shell, gci,&
         map_loc_thermo_gen)

     number = number_of_thermostats
     sum_of_thermostats=number
     CALL mp_sum ( sum_of_thermostats, para_env%group )

!    check = (number==number_of_thermostats)
!    CPPrecondition(check,cp_fatal_level,routineP,failure)
    DEALLOCATE (const_mol)
    DEALLOCATE (tot_const)
    DEALLOCATE (point)

    CALL timestop(handle)

  END SUBROUTINE adiabatic_mapping_region

! *****************************************************************************
!> \brief Performs the real mapping for the thermostat region
!> \param region ...
!> \param map_info ...
!> \param nkind ...
!> \param point ...
!> \param deg_of_freedom ...
!> \param local_molecules ...
!> \param const_mol ...
!> \param massive_atom_list ...
!> \param tot_const ...
!> \param molecule_set ...
!> \param ntherm ...
!> \param shell ...
!> \param gci ...
!> \param map_loc_thermo_gen ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
  SUBROUTINE adiabatic_mapping_region_low( region, map_info, nkind, point,&
       deg_of_freedom, local_molecules, const_mol, massive_atom_list, tot_const,&
       molecule_set, ntherm, shell, gci, map_loc_thermo_gen)

    INTEGER, INTENT(IN)                      :: region
    TYPE(map_info_type), POINTER             :: map_info
    INTEGER                                  :: nkind
    INTEGER, DIMENSION(:, :), POINTER        :: point
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, DIMENSION(:), POINTER           :: const_mol, massive_atom_list, &
                                                tot_const
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    INTEGER, INTENT(OUT)                     :: ntherm
    LOGICAL, INTENT(IN)                      :: shell
    TYPE(global_constraint_type), POINTER    :: gci
    INTEGER, DIMENSION(:), POINTER           :: map_loc_thermo_gen

    CHARACTER(LEN=*), PARAMETER :: routineN = 'adiabatic_mapping_region_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, first_shell, glob_therm_num, handle, icount, &
      ielement, ii, ikind, imol, imol_local, ipart, jj, k, kk, last_atom, &
      last_shell, nglob_cns, nmol_local, number
    LOGICAL                                  :: check, failure, &
                                                global_constraints, &
                                                have_thermostat
    REAL(dp), SAVE, TARGET                   :: unity
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)
    unity          = 1.0_dp
    failure        = .FALSE.
    global_constraints = ASSOCIATED(gci)
    deg_of_freedom = 0
    icount         = 0
    number         = 0
    ntherm         = 0
    nglob_cns      = 0
    IF (global_constraints) nglob_cns = gci%ntot-gci%nrestraint
    IF ( region==do_region_global) THEN
       ! Global Region
       check = ( map_info%dis_type == do_thermo_communication )
       CPASSERT(check)
       DO ikind = 1, nkind
          DO jj = point ( 1, ikind ), point ( 2, ikind )
             IF ( map_loc_thermo_gen ( jj ) /= HUGE ( 0 ) ) THEN
               DO ii = 1, 3
                  map_info%p_kin(ii,jj)%point => map_info%s_kin(1)
                  map_info%p_scale(ii,jj)%point => map_info%v_scale(1)
               END DO
             ELSE
               DO ii = 1, 3
                  NULLIFY ( map_info%p_kin(ii,jj)%point )
                  map_info%p_scale(ii,jj)%point => unity
               END DO
             ENDIF
          END DO
          deg_of_freedom(1) = deg_of_freedom(1) + tot_const(ikind)
          map_info%index(1)     = 1
          map_info%map_index(1) = 1
          number = 1
       END DO
       deg_of_freedom(1) = deg_of_freedom(1) + nglob_cns
    ELSE IF ( region==do_region_molecule) THEN
       ! Molecular Region
       IF   ( map_info%dis_type == do_thermo_no_communication ) THEN
          ! This is the standard case..
          DO ikind = 1, nkind
             nmol_local =  local_molecules % n_el (ikind)
             DO imol_local = 1, nmol_local
                imol = local_molecules%list(ikind)%array(imol_local)
                number = number + 1
                have_thermostat=.TRUE.
! determine if the local molecule belongs to a thermostat
                DO kk = point ( 1, number ), point ( 2, number )
                 !  WRITE ( *, * ) 'kk', size(map_loc_thermo_gen), kk, map_loc_thermo_gen ( kk )
                   IF ( map_loc_thermo_gen ( kk ) == HUGE ( 0 ) ) THEN
                     have_thermostat=.FALSE.
                     EXIT
                   ENDIF
                END DO
! If molecule has a thermostat, then map
                IF ( have_thermostat ) THEN
! glob_therm_num is the global thermostat number attached to the local molecule
! We can test to make sure all atoms in the local molecule belong to the same
! global thermostat as a way to detect errors.
                  glob_therm_num = map_loc_thermo_gen ( point ( 1, number ) )
                  ntherm = ntherm + 1
                  CALL reallocate ( map_info%index, 1, ntherm )
                  CALL reallocate ( map_info%map_index, 1, ntherm )
                  CALL reallocate ( deg_of_freedom, 1, ntherm )
                  map_info%index(ntherm)     = glob_therm_num
                  map_info%map_index(ntherm) = ntherm
                  deg_of_freedom ( ntherm )  = const_mol ( number )
                  DO kk = point ( 1, number ), point ( 2, number )
                     DO jj = 1, 3
                        map_info%p_kin(jj,kk) %point => map_info%s_kin(ntherm)
                        map_info%p_scale(jj,kk) %point => map_info%v_scale(ntherm)
                     END DO
                  END DO
! If molecule has no thermostat, then nullify pointers to that atom
                ELSE
                  DO kk = point ( 1, number ), point ( 2, number )
                    DO jj = 1, 3
                       NULLIFY ( map_info%p_kin(jj,kk) %point )
                       map_info%p_scale(jj,kk) %point => unity
                    END DO
                  END DO
                END IF
             END DO
          END DO
       ELSE IF ( map_info%dis_type == do_thermo_communication ) THEN
          ! This case is quite rare and happens only when we have one molecular
          ! kind and one molecule..
          CPASSERT(nkind==1)
          number = number + 1
          ntherm = ntherm + 1
          map_info%index(ntherm)     = ntherm
          map_info%map_index(ntherm) = ntherm
          deg_of_freedom ( ntherm )  = deg_of_freedom ( ntherm ) + tot_const ( nkind )
          DO kk = point ( 1, nkind ), point ( 2, nkind )
             IF ( map_loc_thermo_gen ( kk ) /= HUGE ( 0 ) ) THEN
               DO jj = 1, 3
                  map_info%p_kin ( jj, kk ) % point => map_info%s_kin(ntherm)
                  map_info%p_scale ( jj, kk ) % point => map_info%v_scale(ntherm)
               END DO
             ELSE
             END IF
               DO jj = 1, 3
                  NULLIFY ( map_info%p_kin ( jj, kk ) % point )
                  map_info%p_scale ( jj, kk ) % point => unity
               END DO
          END DO
       ELSE
          CPASSERT(.FALSE.)
       END IF
       IF (nglob_cns/=0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Molecular thermostats with global constraints are impossible!")
       END IF
    ELSE IF ( region==do_region_massive) THEN
       ! Massive Region
       check = ( map_info%dis_type == do_thermo_no_communication )
       CPASSERT(check)
       DO ikind = 1, nkind
          nmol_local =  local_molecules % n_el ( ikind )
          DO imol_local = 1, nmol_local
             icount = icount + 1
             imol = local_molecules % list ( ikind ) % array ( imol_local )
             molecule => molecule_set ( imol )
             CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom,&
                  first_shell=first_shell, last_shell=last_shell)
             IF (shell) THEN
                first_atom = first_shell
                last_atom  = last_shell
             ELSE
                IF ((tot_const(icount)>0).OR.(nglob_cns/=0)) THEN
                   CALL stop_program(routineN,moduleN,__LINE__,&
                                     "Massive thermostats with constraints are impossible!")
                END IF
             END IF
             k = 0

             have_thermostat=.TRUE.
             DO ii =  point ( 1, icount ), point ( 2, icount )
               IF ( map_loc_thermo_gen ( ii ) /= 1 ) THEN
                  have_thermostat = .FALSE.
                  EXIT
               END IF
             END DO

             IF ( have_thermostat ) THEN
               DO ii =  point ( 1, icount ), point ( 2, icount )
                  ipart = first_atom + k
                  ielement = locate(massive_atom_list,ipart)
                  k = k + 1
                  DO jj = 1, 3
                     ntherm = ntherm + 1
                     CALL reallocate ( map_info%index, 1, ntherm )
                     CALL reallocate ( map_info%map_index, 1, ntherm )
                     map_info%index(ntherm)     = (ielement - 1)*3 + jj
                     map_info%map_index(ntherm) = ntherm
                     map_info%p_kin(jj,ii) %point   => map_info%s_kin(ntherm)
                     map_info%p_scale(jj,ii) %point => map_info%v_scale(ntherm)
                  END DO
               END DO
             ELSE
               DO ii =  point ( 1, icount ), point ( 2, icount )
                 DO jj = 1, 3
                   NULLIFY ( map_info%p_kin(jj,ii) %point )
                   map_info%p_scale(jj,ii) %point => unity
                 END DO
               END DO
             END IF
             IF ( first_atom + k -1 /= last_atom ) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "Inconsistent mapping of particles")
             END IF
          END DO
       END DO
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"Invalid region!")
    END IF

    CALL timestop(handle)

  END SUBROUTINE adiabatic_mapping_region_low


! *****************************************************************************
!> \brief creates the mapping between the system and the thermostats
!> \param dis_type ...
!> \param natoms_local ...
!> \param nmol_local ...
!> \param const_mol ...
!> \param tot_const ...
!> \param point ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param simpar ...
!> \param shell ...
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE adiabatic_region_evaluate(dis_type, natoms_local, nmol_local, const_mol,&
       tot_const, point, local_molecules, molecule_kind_set, molecule_set, simpar, shell)
    INTEGER, INTENT(IN)                      :: dis_type
    INTEGER, INTENT(OUT)                     :: natoms_local, nmol_local
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    TYPE(simpar_type), POINTER               :: simpar
    LOGICAL, INTENT(IN)                      :: shell

    CHARACTER(LEN=*), PARAMETER :: routineN = 'adiabatic_region_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: atm_offset, first_atom, handle, icount, ikind, ilist, imol, &
      imol_local, katom, last_atom, natom, nc, nfixd, nkind, nmol_per_kind, &
      nmolecule, nshell
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)

    natoms_local = 0
    nmol_local   = 0
    nkind = SIZE ( molecule_kind_set )
    NULLIFY(fixd_list, molecule_kind, molecule)
    ! Compute the TOTAL number of molecules and atoms on THIS PROC and
    ! TOTAL number of molecules of IKIND on THIS PROC
    DO ikind = 1, nkind
       molecule_kind => molecule_kind_set ( ikind )
       CALL get_molecule_kind ( molecule_kind, natom=natom, nshell=nshell )
       IF (shell) THEN
          IF (nshell/=0) THEN
             natoms_local = natoms_local + nshell * local_molecules % n_el ( ikind )
             nmol_local   = nmol_local + local_molecules % n_el ( ikind )
          END IF
       ELSE
          natoms_local = natoms_local + natom * local_molecules % n_el ( ikind )
          nmol_local   = nmol_local + local_molecules % n_el ( ikind )
       END IF
    END DO

    CPASSERT(.NOT.ASSOCIATED(const_mol))
    CPASSERT(.NOT.ASSOCIATED(tot_const))
    CPASSERT(.NOT.ASSOCIATED(point))
    IF ( dis_type == do_thermo_no_communication ) THEN
       ALLOCATE ( const_mol (nmol_local))
       ALLOCATE ( tot_const (nmol_local))
       ALLOCATE ( point (2, nmol_local))

       point (:,:)= 0
       atm_offset = 0
       icount = 0
       DO ikind = 1, nkind
          nmol_per_kind =  local_molecules % n_el ( ikind )
          molecule_kind => molecule_kind_set ( ikind )
          CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
               fixd_list=fixd_list, nshell=nshell)
          IF (shell) natom = nshell
          DO imol_local = 1, nmol_per_kind
             icount = icount + 1
             point ( 1, icount )  = atm_offset + 1
             point ( 2, icount )  = atm_offset + natom
             IF (.NOT.shell) THEN
                ! nc keeps track of all constraints but not fixed ones..
                ! Let's identify fixed atoms for this molecule
                nfixd = 0
                imol = local_molecules%list(ikind)%array(imol_local)
                molecule => molecule_set(imol)
                CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom)
                IF (ASSOCIATED(fixd_list)) THEN
                   DO katom = first_atom, last_atom
                      DO ilist = 1, SIZE(fixd_list)
                         IF ( ( katom == fixd_list(ilist)%fixd        ) .AND. &
                              (.NOT. fixd_list(ilist)%restraint%active)) THEN
                            SELECT CASE(fixd_list(ilist)%itype)
                            CASE(use_perd_x,use_perd_y,use_perd_z)
                               nfixd=nfixd+1
                            CASE(use_perd_xy,use_perd_xz,use_perd_yz)
                               nfixd=nfixd+2
                            CASE(use_perd_xyz)
                               nfixd=nfixd+3
                            END SELECT
                         END IF
                      END DO
                   END DO
                END IF
                const_mol ( icount ) = nc + nfixd
                tot_const ( icount ) = const_mol ( icount )
             END IF
             atm_offset =  point ( 2, icount )
          END DO
       END DO
    ELSE IF ( dis_type == do_thermo_communication ) THEN
       ALLOCATE ( const_mol ( nkind ))
       ALLOCATE ( tot_const ( nkind ))
       ALLOCATE ( point ( 2, nkind ))
       point ( :, : ) = 0
       atm_offset     = 0
       ! nc keeps track of all constraints but not fixed ones..
       DO ikind = 1, nkind
          nmol_per_kind =  local_molecules % n_el ( ikind )
          molecule_kind => molecule_kind_set ( ikind )
          CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
               nmolecule=nmolecule, nconstraint_fixd=nfixd,nshell=nshell)
          IF (shell) natom = nshell
          IF (.NOT.shell) THEN
             const_mol ( ikind ) = nc
             ! Let's consider the fixed atoms only for the total number of constraints
             ! in case we are in REPLICATED/INTERACTING thermostats
             tot_const ( ikind ) = const_mol ( ikind ) * nmolecule + nfixd
          END IF
          point ( 1, ikind ) = atm_offset + 1
          point ( 2, ikind ) = atm_offset + natom * nmol_per_kind
          atm_offset =  point (2, ikind)
       END DO
    ENDIF
    IF (( .NOT. simpar % constraint ).OR.shell) THEN
       const_mol = 0
       tot_const = 0
    END IF

    CALL timestop(handle)

  END SUBROUTINE adiabatic_region_evaluate

! *****************************************************************************
!> \brief Main general setup thermostat regions (thermostat independent)
!> \param map_info ...
!> \param deg_of_freedom ...
!> \param massive_atom_list ...
!> \param molecule_kind_set ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param para_env ...
!> \param natoms_local ...
!> \param simpar ...
!> \param number ...
!> \param region ...
!> \param gci ...
!> \param shell ...
!> \param map_loc_thermo_gen ...
!> \param sum_of_thermostats ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
  SUBROUTINE thermostat_mapping_region ( map_info, deg_of_freedom, massive_atom_list,&
       molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, simpar, &
       number, region, gci, shell, map_loc_thermo_gen, sum_of_thermostats)

    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(OUT)                     :: natoms_local
    TYPE(simpar_type), POINTER               :: simpar
    INTEGER, INTENT(IN)                      :: number, region
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: shell
    INTEGER, DIMENSION(:), POINTER           :: map_loc_thermo_gen
    INTEGER, INTENT(IN)                      :: sum_of_thermostats

    CHARACTER(LEN=*), PARAMETER :: routineN = 'thermostat_mapping_region', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nkind, nmol_local, &
                                                nsize, number_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    LOGICAL                                  :: check, failure

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY ( const_mol, tot_const, point)
    CPASSERT(.NOT.ASSOCIATED(deg_of_freedom))
    CPASSERT(.NOT.ASSOCIATED(massive_atom_list))

    nkind = SIZE(molecule_kind_set)
    CALL mapping_region_evaluate(map_info%dis_type, natoms_local, nmol_local,&
         const_mol, tot_const, point, local_molecules, molecule_kind_set, molecule_set,&
         region, simpar, shell, map_loc_thermo_gen, sum_of_thermostats, para_env)

    ! Now we can allocate the target array s_kin and p_kin..
    SELECT CASE(region)
    CASE(do_region_global, do_region_molecule, do_region_massive)
       nsize = number
    CASE(do_region_defined)
       nsize = sum_of_thermostats
    END SELECT
    ALLOCATE (map_info%s_kin(nsize))
    ALLOCATE (map_info%v_scale(nsize))
    ALLOCATE (map_info%p_kin(3,natoms_local))
    ALLOCATE (map_info%p_scale(3,natoms_local))
    ! Allocate index array
    ALLOCATE ( map_info%index(number))
    ALLOCATE ( map_info%map_index(number))
    ALLOCATE (deg_of_freedom(number))

    CALL massive_list_generate (molecule_set, molecule_kind_set, &
         local_molecules, para_env, massive_atom_list, region, shell)

    CALL  thermostat_mapping_region_low(region, map_info, nkind, point,&
         deg_of_freedom, local_molecules, const_mol, massive_atom_list,&
         tot_const, molecule_set, number_of_thermostats, shell, gci,&
         map_loc_thermo_gen)

    check = (number==number_of_thermostats)
    CPASSERT(check)
    DEALLOCATE (const_mol)
    DEALLOCATE (tot_const)
    DEALLOCATE (point)

    CALL timestop(handle)

  END SUBROUTINE thermostat_mapping_region

! *****************************************************************************
!> \brief Performs the real mapping for the thermostat region
!> \param region ...
!> \param map_info ...
!> \param nkind ...
!> \param point ...
!> \param deg_of_freedom ...
!> \param local_molecules ...
!> \param const_mol ...
!> \param massive_atom_list ...
!> \param tot_const ...
!> \param molecule_set ...
!> \param number ...
!> \param shell ...
!> \param gci ...
!> \param map_loc_thermo_gen ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
  SUBROUTINE  thermostat_mapping_region_low(region, map_info, nkind, point,&
       deg_of_freedom, local_molecules, const_mol, massive_atom_list, tot_const,&
       molecule_set, number, shell, gci, map_loc_thermo_gen)

    INTEGER, INTENT(IN)                      :: region
    TYPE(map_info_type), POINTER             :: map_info
    INTEGER                                  :: nkind
    INTEGER, DIMENSION(:, :), POINTER        :: point
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, DIMENSION(:), POINTER           :: const_mol, massive_atom_list, &
                                                tot_const
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    INTEGER, INTENT(OUT)                     :: number
    LOGICAL, INTENT(IN)                      :: shell
    TYPE(global_constraint_type), POINTER    :: gci
    INTEGER, DIMENSION(:), POINTER           :: map_loc_thermo_gen

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'thermostat_mapping_region_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, first_shell, handle, i, icount, ielement, ii, &
      ikind, imap, imol, imol_local, ipart, itmp, jj, k, kk, last_atom, &
      last_shell, nglob_cns, nmol_local
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp, wrk
    LOGICAL                                  :: check, failure, &
                                                global_constraints
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)

    failure        = .FALSE.
    global_constraints = ASSOCIATED(gci)
    deg_of_freedom = 0
    icount         = 0
    number         = 0
    nglob_cns      = 0
    IF (global_constraints) nglob_cns = gci%ntot-gci%nrestraint
    IF ( region==do_region_global) THEN
       ! Global Region
       check = ( map_info%dis_type == do_thermo_communication )
       CPASSERT(check)
       DO ikind = 1, nkind
          DO jj = point ( 1, ikind ), point ( 2, ikind )
             DO ii = 1, 3
                map_info%p_kin(ii,jj)%point => map_info%s_kin(1)
                map_info%p_scale(ii,jj)%point => map_info%v_scale(1)
             END DO
          END DO
          deg_of_freedom(1) = deg_of_freedom(1) + tot_const(ikind)
          map_info%index(1)     = 1
          map_info%map_index(1) = 1
          number = 1
       END DO
       deg_of_freedom(1) = deg_of_freedom(1) + nglob_cns
    ELSE IF ( region==do_region_defined) THEN
       ! User defined Region to thermostat
       check = ( map_info%dis_type == do_thermo_communication )
       CPASSERT(check)
       ! Lets' identify the matching of the local thermostat w.r.t. the global one
       itmp = SIZE(map_loc_thermo_gen)
       ALLOCATE(tmp(itmp))
       ALLOCATE(wrk(itmp))
       tmp(:) = map_loc_thermo_gen
       CALL sort(tmp, itmp, wrk)
       number = 1
       map_info%index(number)     = tmp(1)
       map_info%map_index(number) = tmp(1)
       deg_of_freedom(number)     = tot_const(tmp(1))
       DO i = 2, itmp
          IF (tmp(i)/=tmp(i-1)) THEN
             number = number + 1
             map_info%index(number)     = tmp(i)
             map_info%map_index(number) = tmp(i)
             deg_of_freedom(number) = tot_const(tmp(i))
          END IF
       END DO
       DEALLOCATE(tmp)
       DEALLOCATE(wrk)
       DO jj = 1, SIZE(map_loc_thermo_gen)
          DO ii = 1, 3
             imap = map_loc_thermo_gen(jj)
             map_info%p_kin(ii,jj)%point => map_info%s_kin(imap)
             map_info%p_scale(ii,jj)%point => map_info%v_scale(imap)
          END DO
       END DO
       IF (nglob_cns/=0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "User Defined thermostats with global constraints not implemented!")
       END IF
    ELSE IF ( region==do_region_molecule) THEN
       ! Molecular Region
       IF      ( map_info%dis_type == do_thermo_no_communication ) THEN
          ! This is the standard case..
          DO ikind = 1, nkind
             nmol_local =  local_molecules % n_el (ikind)
             DO imol_local = 1, nmol_local
                imol = local_molecules%list(ikind)%array(imol_local)
                number = number + 1
                map_info%index(number)     = imol
                map_info%map_index(number) = number
                deg_of_freedom ( number )  = const_mol ( number )
                DO kk = point ( 1, number ), point ( 2, number )
                   DO jj = 1, 3
                      map_info%p_kin(jj,kk) %point => map_info%s_kin(number)
                      map_info%p_scale(jj,kk) %point => map_info%v_scale(number)
                   END DO
                END DO
             END DO
          END DO
       ELSE IF ( map_info%dis_type == do_thermo_communication ) THEN
          ! This case is quite rare and happens only when we have one molecular
          ! kind and one molecule..
          CPASSERT(nkind==1)
          number = number + 1
          map_info%index(number)     = number
          map_info%map_index(number) = number
          deg_of_freedom ( number )  = deg_of_freedom ( number ) + tot_const ( nkind )
          DO kk = point ( 1, nkind ), point ( 2, nkind )
             DO jj = 1, 3
                map_info%p_kin ( jj, kk ) % point => map_info%s_kin(number)
                map_info%p_scale ( jj, kk ) % point => map_info%v_scale(number)
             END DO
          END DO
       ELSE
          CPASSERT(.FALSE.)
       END IF
       IF (nglob_cns/=0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Molecular thermostats with global constraints are impossible!")
       END IF
    ELSE IF ( region==do_region_massive) THEN
       ! Massive Region
       check = ( map_info%dis_type == do_thermo_no_communication )
       CPASSERT(check)
       DO ikind = 1, nkind
          nmol_local =  local_molecules % n_el ( ikind )
          DO imol_local = 1, nmol_local
             icount = icount + 1
             imol = local_molecules % list ( ikind ) % array ( imol_local )
             molecule => molecule_set ( imol )
             CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom,&
                  first_shell=first_shell, last_shell=last_shell)
             IF (shell) THEN
                first_atom = first_shell
                last_atom  = last_shell
             ELSE
                IF ((tot_const(icount)>0).OR.(nglob_cns/=0)) THEN
                   CALL stop_program(routineN,moduleN,__LINE__,&
                                     "Massive thermostats with constraints are impossible!")
                END IF
             END IF
             k = 0
             DO ii =  point ( 1, icount ), point ( 2, icount )
                ipart = first_atom + k
                ielement = locate(massive_atom_list,ipart)
                k = k + 1
                DO jj = 1, 3
                   number = number + 1
                   map_info%index(number)     = (ielement - 1)*3 + jj
                   map_info%map_index(number) = number
                   map_info%p_kin(jj,ii) %point   => map_info%s_kin(number)
                   map_info%p_scale(jj,ii) %point => map_info%v_scale(number)
                END DO
             END DO
             IF ( first_atom + k -1 /= last_atom ) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "Inconsistent mapping of particles")
             END IF
          END DO
       END DO
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"Invalid region!")
    END IF

    CALL timestop(handle)

  END SUBROUTINE thermostat_mapping_region_low

! *****************************************************************************
!> \brief creates the mapping between the system and the thermostats
!> \param dis_type ...
!> \param natoms_local ...
!> \param nmol_local ...
!> \param const_mol ...
!> \param tot_const ...
!> \param point ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param region ...
!> \param simpar ...
!> \param shell ...
!> \param map_loc_thermo_gen ...
!> \param sum_of_thermostats ...
!> \param para_env ...
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE mapping_region_evaluate(dis_type, natoms_local, nmol_local, const_mol,&
       tot_const, point, local_molecules, molecule_kind_set, molecule_set, region,&
       simpar, shell, map_loc_thermo_gen, sum_of_thermostats, para_env)
    INTEGER, INTENT(IN)                      :: dis_type
    INTEGER, INTENT(OUT)                     :: natoms_local, nmol_local
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    INTEGER, INTENT(IN)                      :: region
    TYPE(simpar_type), POINTER               :: simpar
    LOGICAL, INTENT(IN)                      :: shell
    INTEGER, DIMENSION(:), POINTER           :: map_loc_thermo_gen
    INTEGER, INTENT(IN)                      :: sum_of_thermostats
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mapping_region_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: atm_offset, first_atom, handle, i, iatm, icount, id_region, &
      ikind, ilist, imol, imol_local, j, jatm, katom, last_atom, natom, nc, &
      nfixd, nkind, nmol_per_kind, nmolecule, nshell
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)

    natoms_local = 0
    nmol_local   = 0
    nkind = SIZE ( molecule_kind_set )
    NULLIFY(fixd_list, molecule_kind, molecule, colv_list, g3x3_list, g4x6_list)
    ! Compute the TOTAL number of molecules and atoms on THIS PROC and
    ! TOTAL number of molecules of IKIND on THIS PROC
    DO ikind = 1, nkind
       molecule_kind => molecule_kind_set ( ikind )
       CALL get_molecule_kind ( molecule_kind, natom=natom, nshell=nshell )
       IF (shell) THEN
          IF (nshell/=0) THEN
             natoms_local = natoms_local + nshell * local_molecules % n_el ( ikind )
             nmol_local   = nmol_local + local_molecules % n_el ( ikind )
          END IF
       ELSE
          natoms_local = natoms_local + natom * local_molecules % n_el ( ikind )
          nmol_local   = nmol_local + local_molecules % n_el ( ikind )
       END IF
    END DO

    CPASSERT(.NOT.ASSOCIATED(const_mol))
    CPASSERT(.NOT.ASSOCIATED(tot_const))
    CPASSERT(.NOT.ASSOCIATED(point))
    IF ( dis_type == do_thermo_no_communication ) THEN
       ALLOCATE ( const_mol (nmol_local))
       ALLOCATE ( tot_const (nmol_local))
       ALLOCATE ( point (2, nmol_local))

       point (:,:)= 0
       atm_offset = 0
       icount = 0
       DO ikind = 1, nkind
          nmol_per_kind =  local_molecules % n_el ( ikind )
          molecule_kind => molecule_kind_set ( ikind )
          CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
               fixd_list=fixd_list, nshell=nshell)
          IF (shell) natom = nshell
          DO imol_local = 1, nmol_per_kind
             icount = icount + 1
             point ( 1, icount )  = atm_offset + 1
             point ( 2, icount )  = atm_offset + natom
             IF (.NOT.shell) THEN
                ! nc keeps track of all constraints but not fixed ones..
                ! Let's identify fixed atoms for this molecule
                nfixd = 0
                imol = local_molecules%list(ikind)%array(imol_local)
                molecule => molecule_set(imol)
                CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom)
                IF (ASSOCIATED(fixd_list)) THEN
                   DO katom = first_atom, last_atom
                      DO ilist = 1, SIZE(fixd_list)
                         IF ( ( katom == fixd_list(ilist)%fixd        ) .AND. &
                              (.NOT. fixd_list(ilist)%restraint%active)) THEN
                            SELECT CASE(fixd_list(ilist)%itype)
                            CASE(use_perd_x,use_perd_y,use_perd_z)
                               nfixd=nfixd+1
                            CASE(use_perd_xy,use_perd_xz,use_perd_yz)
                               nfixd=nfixd+2
                            CASE(use_perd_xyz)
                               nfixd=nfixd+3
                            END SELECT
                         END IF
                      END DO
                   END DO
                END IF
                const_mol ( icount ) = nc + nfixd
                tot_const ( icount ) = const_mol ( icount )
             END IF
             atm_offset =  point ( 2, icount )
          END DO
       END DO
    ELSE IF ( dis_type == do_thermo_communication ) THEN
       IF (region==do_region_defined) THEN
          ! Setup of the arbitrary region
          ALLOCATE ( tot_const (sum_of_thermostats))
          ALLOCATE ( point (2, 0))
          ALLOCATE ( const_mol (0))
          atm_offset = 0
          tot_const  = 0
          const_mol  = 0
          point      = 0
          DO ikind = 1, nkind
             nmol_per_kind =  local_molecules % n_el ( ikind )
             molecule_kind => molecule_kind_set ( ikind )
             CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
                  fixd_list=fixd_list, colv_list=colv_list, g3x3_list=g3x3_list,&
                  g4x6_list=g4x6_list, nshell=nshell)
             IF (shell) natom = nshell
             DO imol_local = 1, nmol_per_kind
                IF (.NOT.shell) THEN
                   ! First if nc is not zero let's check if all atoms of a molecule
                   ! are in the same thermostatting region..
                   imol = local_molecules%list(ikind)%array(imol_local)
                   molecule => molecule_set(imol)
                   id_region = map_loc_thermo_gen(atm_offset+1)
                   IF (ALL(map_loc_thermo_gen(atm_offset+1:atm_offset+natom)==id_region)) THEN
                      ! All the atoms of a molecule are within the same thermostatting
                      ! region.. this is the easy case..
                      tot_const(id_region) = tot_const(id_region) + nc
                   ELSE
                      ! If not let's check the single constraints defined for this molecule
                      ! and continue only when atoms involved in the constraint belong to
                      ! the same thermostatting region
                      IF (ASSOCIATED(colv_list)) THEN
                         DO i = 1, SIZE(colv_list)
                            IF (.NOT.colv_list(i)%restraint%active) THEN
                               iatm = atm_offset + colv_list(i)%i_atoms(1)
                               DO j = 2, SIZE(colv_list(i)%i_atoms)
                                  jatm = atm_offset + colv_list(i)%i_atoms(j)
                                  IF (map_loc_thermo_gen(iatm)/=map_loc_thermo_gen(jatm)) THEN
                                     CALL stop_program(routineN,moduleN,__LINE__,&
                                          "User Defined Region: "//&
                                          "A constraint (COLV) was defined between two thermostatting regions! "//&
                                          "This is not allowed!")
                                  END IF
                               END DO
                               id_region = map_loc_thermo_gen(iatm)
                               tot_const(id_region) = tot_const(id_region) + 1
                            END IF
                         END DO
                      END IF
                      IF (ASSOCIATED(g3x3_list)) THEN
                         DO i = 1, SIZE(g3x3_list)
                            IF (.NOT.g3x3_list(i)%restraint%active) THEN
                               iatm = atm_offset + g3x3_list(i)%a
                               jatm = atm_offset + g3x3_list(i)%b
                               IF (map_loc_thermo_gen(iatm)/=map_loc_thermo_gen(jatm)) THEN
                                  CALL stop_program(routineN,moduleN,__LINE__,&
                                       "User Defined Region: "//&
                                       "A constraint (G3X3) was defined between two thermostatting regions! "//&
                                       "This is not allowed!")
                               END IF
                               jatm = atm_offset + g3x3_list(i)%c
                               IF (map_loc_thermo_gen(iatm)/=map_loc_thermo_gen(jatm)) THEN
                                  CALL stop_program(routineN,moduleN,__LINE__,&
                                       "User Defined Region: "//&
                                       "A constraint (G3X3) was defined between two thermostatting regions! "//&
                                       "This is not allowed!")
                               END IF
                            END IF
                            id_region = map_loc_thermo_gen(iatm)
                            tot_const(id_region) = tot_const(id_region) + 3
                         END DO
                      END IF
                      IF (ASSOCIATED(g4x6_list)) THEN
                         DO i = 1, SIZE(g4x6_list)
                            IF (.NOT.g4x6_list(i)%restraint%active) THEN
                               iatm = atm_offset + g4x6_list(i)%a
                               jatm = atm_offset + g4x6_list(i)%b
                               IF (map_loc_thermo_gen(iatm)/=map_loc_thermo_gen(jatm)) THEN
                                  CALL stop_program(routineN,moduleN,__LINE__,&
                                       " User Defined Region: "//&
                                       "A constraint (G4X6) was defined between two thermostatting regions! "//&
                                       "This is not allowed!")
                               END IF
                               jatm = atm_offset + g4x6_list(i)%c
                               IF (map_loc_thermo_gen(iatm)/=map_loc_thermo_gen(jatm)) THEN
                                  CALL stop_program(routineN,moduleN,__LINE__,&
                                       " User Defined Region: "//&
                                       "A constraint (G4X6) was defined between two thermostatting regions! "//&
                                       "This is not allowed!")
                               END IF
                               jatm = atm_offset + g4x6_list(i)%d
                               IF (map_loc_thermo_gen(iatm)/=map_loc_thermo_gen(jatm)) THEN
                                  CALL stop_program(routineN,moduleN,__LINE__,&
                                       " User Defined Region: "//&
                                       "A constraint (G4X6) was defined between two thermostatting regions! "//&
                                       "This is not allowed!")
                               END IF
                            END IF
                            id_region = map_loc_thermo_gen(iatm)
                            tot_const(id_region) = tot_const(id_region) + 6
                         END DO
                      END IF
                   END IF
                   ! Here we handle possibly fixed atoms
                   IF (ASSOCIATED(fixd_list)) THEN
                      CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom)
                      iatm = 0
                      DO katom = first_atom, last_atom
                         iatm = iatm + 1
                         DO ilist = 1, SIZE(fixd_list)
                            IF ( ( katom == fixd_list(ilist)%fixd        ) .AND. &
                                 (.NOT. fixd_list(ilist)%restraint%active)) THEN
                               id_region = map_loc_thermo_gen(atm_offset+iatm)
                               SELECT CASE(fixd_list(ilist)%itype)
                               CASE(use_perd_x, use_perd_y, use_perd_z)
                                  tot_const(id_region) = tot_const(id_region) + 1
                               CASE(use_perd_xy, use_perd_xz, use_perd_yz)
                                  tot_const(id_region) = tot_const(id_region) + 2
                               CASE(use_perd_xyz)
                                  tot_const(id_region) = tot_const(id_region) + 3
                               END SELECT
                            END IF
                         END DO
                      END DO
                   END IF
                END IF
                atm_offset =  atm_offset + natom
             END DO
          END DO
          CALL mp_sum(tot_const, para_env%group)
       ELSE
          ALLOCATE ( const_mol ( nkind ))
          ALLOCATE ( tot_const ( nkind ))
          ALLOCATE ( point ( 2, nkind ))
          point ( :, : ) = 0
          atm_offset     = 0
          ! nc keeps track of all constraints but not fixed ones..
          DO ikind = 1, nkind
             nmol_per_kind =  local_molecules % n_el ( ikind )
             molecule_kind => molecule_kind_set ( ikind )
             CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
                  nmolecule=nmolecule, nconstraint_fixd=nfixd,nshell=nshell)
             IF (shell) natom = nshell
             IF (.NOT.shell) THEN
                const_mol ( ikind ) = nc
                ! Let's consider the fixed atoms only for the total number of constraints
                ! in case we are in REPLICATED/INTERACTING thermostats
                tot_const ( ikind ) = const_mol ( ikind ) * nmolecule + nfixd
             END IF
             point ( 1, ikind ) = atm_offset + 1
             point ( 2, ikind ) = atm_offset + natom * nmol_per_kind
             atm_offset =  point (2, ikind)
          END DO
       END IF
    ENDIF
    IF (( .NOT. simpar % constraint ).OR.shell) THEN
       const_mol = 0
       tot_const = 0
    END IF

    CALL timestop(handle)

  END SUBROUTINE mapping_region_evaluate

! *****************************************************************************
!> \brief ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param local_molecules ...
!> \param para_env ...
!> \param massive_atom_list ...
!> \param region ...
!> \param shell ...
! *****************************************************************************
  SUBROUTINE massive_list_generate ( molecule_set, molecule_kind_set, &
       local_molecules, para_env, massive_atom_list, region, shell)

    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, POINTER                         :: massive_atom_list( : )
    INTEGER, INTENT(IN)                      :: region
    LOGICAL, INTENT(IN)                      :: shell

    CHARACTER(LEN=*), PARAMETER :: routineN = 'massive_list_generate', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, first_shell, handle, i, ikind, imol, iproc, j, &
      natom, ncount, nkind, nmol_per_kind, nshell, num_massive_atm, &
      num_massive_atm_local, offset, stat
    INTEGER, DIMENSION(:), POINTER           :: array_num_massive_atm, &
                                                local_atm_list, work
    LOGICAL                                  :: failure
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)

    failure = .FALSE.
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
          CALL get_molecule_kind(molecule_kind,natom=natom,nshell=nshell)
          IF(region == do_region_massive) THEN
             IF (shell) THEN
                natom = nshell
             END IF
             num_massive_atm_local = num_massive_atm_local + natom
             CALL reallocate(local_atm_list,1,num_massive_atm_local)
             CALL get_molecule (molecule,first_atom=first_atom,first_shell=first_shell)
             IF (shell) THEN
                first_atom = first_shell
             END IF
             DO j=1,natom
                local_atm_list(num_massive_atm_local-natom+j) = first_atom -1 + j
             END DO
          END IF
       END DO
    END DO

    ALLOCATE(array_num_massive_atm(para_env%num_pe))
    CALL mp_allgather(num_massive_atm_local,array_num_massive_atm,para_env%group)

    num_massive_atm = SUM(array_num_massive_atm)
    ALLOCATE(massive_atom_list(num_massive_atm))

    offset = 0
    DO iproc=1,para_env%num_pe
       ncount = array_num_massive_atm(iproc)
       ALLOCATE(work(ncount))
       IF(para_env%mepos == (iproc-1)) THEN
          DO i=1,ncount
             work(i) = local_atm_list(i)
          END DO
       ELSE
          work(:) = 0
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       DO i=1,ncount
          massive_atom_list(offset+i) = work(i)
       END DO
       DEALLOCATE(work)
       offset = offset + array_num_massive_atm(iproc)
    END DO

    ! Sort atom list
    ALLOCATE (work(num_massive_atm),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "work",int_size*num_massive_atm)
    CALL sort(massive_atom_list,num_massive_atm,work)
    DEALLOCATE (work)

    DEALLOCATE(local_atm_list)
    DEALLOCATE(array_num_massive_atm)

    CALL timestop(handle)

  END SUBROUTINE massive_list_generate

! *****************************************************************************
!> \brief Initialize the map_info for barostat thermostat
!> \param map_info ...
!> \param ndeg ...
!> \param num_thermo ...
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE init_baro_map_info(map_info, ndeg, num_thermo)

    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, INTENT(IN)                      :: ndeg, num_thermo

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_baro_map_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i

    CALL timeset(routineN,handle)

    ALLOCATE (map_info%s_kin(num_thermo))
    ALLOCATE (map_info%v_scale(num_thermo))
    ALLOCATE (map_info%p_kin(1,ndeg))
    ALLOCATE (map_info%p_scale(1,ndeg))
    ! Allocate the index array
    ALLOCATE (map_info%index(1))
    ALLOCATE (map_info%map_index(1))

    ! Begin the mapping loop
    DO i = 1, ndeg
       map_info%p_kin(1,i)%point   => map_info%s_kin(1)
       map_info%p_scale(1,i)%point => map_info%v_scale(1)
    END DO
    map_info%index(1)     = 1
    map_info%map_index(1) = 1

    CALL timestop(handle)

  END SUBROUTINE init_baro_map_info

END MODULE thermostat_mapping

