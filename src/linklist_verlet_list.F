!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linklist_verlet_list [1.0] *
!!
!!   NAME
!!     linklist_verlet_list
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_verlet_list

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp_size,&
                                             int_size
  USE linklist_types,                  ONLY: linklist_internal_data_type,&
                                             subcell_data_type
  USE linklist_utilities,              ONLY: get_match
  USE molecule_types,                  ONLY: linklist_exclusion,&
                                             particle_node_type
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_list_types,          ONLY: add_neighbor_list,&
                                             add_neighbor_node,&
                                             allocate_neighbor_list_set,&
                                             deallocate_neighbor_list_set,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_node,&
                                             init_neighbor_list,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: include_list,&
                                             sort

  IMPLICIT NONE

  INTEGER :: output_unit
  LOGICAL :: print_cell_parameters
  CHARACTER(LEN=*), PARAMETER :: module_name = "linklist_verlet_list"


! *** Global types of the module ***

  TYPE atoms_type
    INTEGER, DIMENSION(:), POINTER :: list
  END TYPE atoms_type

  TYPE pbc_coord_type
    REAL(dbl), DIMENSION(:,:), POINTER :: s
  END TYPE pbc_coord_type

  PRIVATE
  PUBLIC :: build_verlet_lists, rebuild_verlet_lists

!!***
!******************************************************************************

CONTAINS

!******************************************************************************

! to get the Verlet neighbor list for nonbonded interactions

  SUBROUTINE build_verlet_lists ( atomic_kind_set, part, cell, ll_data, pnode )

!   Purpose: Build all the required neighbor lists FIST

!   History: - Creation (28.08.2000,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), POINTER                 :: cell
    TYPE(linklist_internal_data_type), &
      INTENT(IN)                             :: ll_data
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode

    CHARACTER(LEN=*), PARAMETER :: module_name = "linklist_verlet_list", &
      routine_name = "build_verlet_lists", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER                                  :: atom_a, group, handle, ikind, &
                                                inode, istat, mype, natoms, &
                                                nkind, nnodes, npe, &
                                                periodic( 3 )
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom
    LOGICAL                                  :: ionode
    REAL(dbl)                                :: lat_vec( 3 ), ra( 3 ), &
                                                sa( 3 ), subcells, &
                                                unit_of_length
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atoms
    TYPE(global_environment_type)            :: globenv
    TYPE(pbc_coord_type)                     :: pbc_coord

!   ---------------------------------------------------------------------------

    globenv = ll_data % globenv

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)
  

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    IF ((ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > LOW)) THEN
      print_cell_parameters = .TRUE.
    ELSE
      print_cell_parameters = .FALSE.
    END IF

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Allocate work storage ***

    nkind = SIZE ( atomic_kind_set )

    ALLOCATE ( natom ( nkind ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( routine, "natom", nkind * int_size )

    ALLOCATE ( atoms ( nkind ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( routine, "atoms", nkind * int_size )

!   *** fill atomic kind info ***

    DO ikind = 1, nkind

      atomic_kind => atomic_kind_set ( ikind )

      NULLIFY ( atoms ( ikind ) % list )

      CALL get_atomic_kind ( atomic_kind = atomic_kind, &
                             natom = natom ( ikind ), &
                             atom_list = atoms ( ikind ) % list )
    END DO

!   *** Calculate PBC coordinates ***
    natoms = SIZE ( part )
    ALLOCATE ( pbc_coord % s ( 3, natoms ), STAT = istat )
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s",&
                                     3*natoms*dp_size)
    nnodes = SIZE ( pnode )
    DO inode=1,nnodes
       
      ALLOCATE ( pnode ( inode ) % nl, STAT = istat )
      IF ( istat /= 0 ) CALL stop_memory ( routine, "node", 1 )
      atom_a = pnode ( inode ) % p % iatom
      ra = pnode ( inode ) % p % r
      sa  = real_to_scaled ( ra ( : ), cell )
      lat_vec ( 1 ) = periodic(1)*INT(sa(1)+SIGN(0.5_dbl,sa(1)))
      lat_vec ( 2 ) = periodic(2)*INT(sa(2)+SIGN(0.5_dbl,sa(2)))
      lat_vec ( 3 ) = periodic(3)*INT(sa(3)+SIGN(0.5_dbl,sa(3)))
      pbc_coord % s ( :, atom_a ) = sa ( : ) - lat_vec ( : )

    END DO

!   *** Build the verlet neighbor lists ***

    CALL build_verlet_neighbor_lists( ll_data % subcells, ll_data % rlist_cut,  &
                                      ll_data % rlist_cutsq, pbc_coord, atoms, &
                                      natom, ll_data % sab_ppnl, cell,         &
                                      ll_data % ll_cell, part, pnode )

    IF (ionode.AND.globenv%print%sab_orb_neighbor_lists) THEN
      CALL write_neighbor_lists(part, cell, ll_data % sab_ppnl,unit_of_length_name, &
                                unit_of_length,"VERLET NEIGHBOR LISTS")
    END IF

!   *** Release work storage ***

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      NULLIFY (atoms(ikind)%list)
    ENDDO

    DEALLOCATE (pbc_coord%s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s")

    DEALLOCATE (atoms,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atoms")

    CALL timestop(0.0_dbl,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_verlet_lists

! *****************************************************************************

  SUBROUTINE build_verlet_neighbor_lists( subcells, rlist_cut, rlist_cutsq,  &
                                          pbc_coord, atoms, natom, sab_orb,  &
                                          cell, ll_cell, part, pnode )

!   Purpose: Build the neighbor lists for calculation of the orbital-orbital
!            overlap matrix.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************
    REAL(dbl), INTENT(IN)                    :: subcells, rlist_cut( :, : ), &
                                                rlist_cutsq( :, : )
    TYPE(pbc_coord_type), INTENT(IN)         :: pbc_coord
    TYPE(atoms_type), INTENT(IN)             :: atoms( : )
    INTEGER, INTENT(IN)                      :: natom( : )
    TYPE(neighbor_list_set_p_type), POINTER  :: sab_orb( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(subcell_data_type), POINTER         :: ll_cell( :, : )
    TYPE(particle_type), INTENT(IN)          :: part( : )
    TYPE(particle_node_type), INTENT(INOUT)  :: pnode( : )

    CHARACTER(LEN=*), PARAMETER :: &
      routine_name = "build_verlet_neighbor_lists", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"
    INTEGER, PARAMETER                       :: els = 50

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: ab, atom_a, atom_b, cell_a( 3 ), handle, iatom, id, iexclude, &
      igrid, iijk, ikind, inode, istat, jatom, jgrid, jkind, kgrid, nexcl, &
      nkinds, nnode, periodic( 3 )
    INTEGER, DIMENSION(:), POINTER           :: elist, ind, ncell, ngrid
    INTEGER, POINTER                         :: ijk( :, :, :, : ), &
                                                nijk( :, :, : )
    LOGICAL                                  :: wrap_x, wrap_y, wrap_z
    REAL(dbl) :: a_max, a_min, abc( 3 ), b_max, b_min, c_max, c_min, r( 3 ), &
      ra( 3 ), rab2_max, rab_max, rb( 3 ), sa_max( 3 ), sa_max_pbc( 3 ), &
      sa_min( 3 ), sa_min_pbc( 3 ), sa_pbc( 3 ), sb_pbc( 3 ), unit_of_length
    REAL(dbl), POINTER                       :: grid_max( :, :, :, : ), &
                                                grid_min( :, :, :, : ), &
                                                sab_max( : )
    TYPE(linklist_exclusion), POINTER        :: ll_exclude

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  abc=abc, &
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)
    nkinds = SIZE ( rlist_cut, 1 )

! Check for ncells and ngrids
    ALLOCATE ( ll_cell ( nkinds, nkinds ), STAT = istat )
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"ll_cell",nkinds*nkinds*dp_size)
    END IF
    DO ikind = 1, nkinds
      DO jkind = 1, nkinds
        NULLIFY ( ll_cell ( ikind, jkind ) % grid_min )
        NULLIFY ( ll_cell ( ikind, jkind ) % grid_max )
        NULLIFY ( ll_cell ( ikind, jkind ) % nijk )
        NULLIFY ( ll_cell ( ikind, jkind ) % ijk )
        rab2_max = rlist_cutsq ( ikind, jkind )
        rab_max = rlist_cut ( ikind, jkind )
        r(:) = rab_max
        ll_cell ( ikind, jkind ) % sab_max(:) = real_to_scaled ( r ( : ), cell )
        sab_max => ll_cell ( ikind, jkind ) % sab_max

        DO id = 1, 3
          IF ( sab_max ( id ) <=0.5_dbl ) THEN
            ll_cell ( ikind, jkind ) % ncell ( id ) = 0
          ELSE
            ll_cell ( ikind, jkind ) % ncell ( id ) = periodic ( id ) *  &
                                           CEILING ( sab_max ( id ) )
          END IF
        END DO

        ll_cell(ikind,jkind)%ngrid(:) = MAX(1,NINT(0.5_dbl*subcells/sab_max(:)))

        ncell => ll_cell ( ikind, jkind ) % ncell
        ngrid => ll_cell ( ikind, jkind ) % ngrid
        IF ( ikind <= jkind ) THEN
          WRITE ( output_unit, '(A,T50,2(I2),T60,3(I3))' )  &
                 ' BUILD_VERLET_LISTS| NCELL for interaction', ikind, jkind, ncell ( : )
          WRITE ( output_unit, '(A,T50,2(I2),T60,3(I3))' )  &
                 ' BUILD_VERLET_LISTS| NGRID for interaction', ikind, jkind, ngrid ( : )
        END IF
        IF ( SUM ( ngrid ) /= 3 ) THEN
          ALLOCATE(ll_cell(ikind,jkind)%grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*dp_size)
          END IF
          grid_min => ll_cell ( ikind, jkind ) % grid_min

          ALLOCATE(ll_cell(ikind,jkind)%grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*dp_size)
          END IF
          grid_max => ll_cell ( ikind, jkind ) % grid_max

          ALLOCATE(ll_cell(ikind,jkind)%nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk => ll_cell ( ikind, jkind ) % nijk
          nijk(:,:,:) = 0

          ALLOCATE(ll_cell(ikind,jkind)%ijk(natom(jkind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",&
                           natom(jkind)*PRODUCT(ngrid)*int_size)
          END IF 
          ijk => ll_cell ( ikind, jkind ) % ijk

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,dbl)/REAL(ngrid(1),dbl) - 0.5_dbl
            a_max = REAL(igrid,dbl)/REAL(ngrid(1),dbl) - 0.5_dbl
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,dbl)/REAL(ngrid(2),dbl) - 0.5_dbl
              b_max = REAL(jgrid,dbl)/REAL(ngrid(2),dbl) - 0.5_dbl
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,dbl)/REAL(ngrid(3),dbl) - 0.5_dbl
                c_max = REAL(kgrid,dbl)/REAL(ngrid(3),dbl) - 0.5_dbl
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO jatom=1,natom(jkind)
            atom_b = atoms(jkind)%list(jatom)
            sb_pbc(:) = pbc_coord % s ( :, atom_b )
            igrid = MAX(1,CEILING((sb_pbc(1) + 0.5_dbl)*ngrid(1)))
            jgrid = MAX(1,CEILING((sb_pbc(2) + 0.5_dbl)*ngrid(2)))
            kgrid = MAX(1,CEILING((sb_pbc(3) + 0.5_dbl)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = atom_b
          END DO
        END IF
      END DO
    END DO

    ALLOCATE ( elist ( els ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'elist', els )

    ALLOCATE ( ind ( els ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'ind', els )

    IF (ASSOCIATED(sab_orb)) THEN
      DO ab=1,SIZE(sab_orb)
        CALL deallocate_neighbor_list_set(sab_orb(ab)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sab_orb")
    END IF

    nnode = SIZE ( pnode )

    ALLOCATE ( sab_orb ( nnode ), STAT = istat )
    IF (istat /= 0) CALL stop_memory ( routine, "sab_orb", nnode )

    DO ab = 1, SIZE ( sab_orb )
      NULLIFY ( sab_orb ( ab ) % neighbor_list_set )
    END DO

!   *** Print headline ***

    IF (print_cell_parameters) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE VERLET NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all nodes ***
    DO inode = 1, nnode
      CALL get_atomic_kind ( pnode ( inode ) % p % atomic_kind, &
                             kind_number = ikind )
      atom_a = pnode ( inode ) % p % iatom
      ra = part ( atom_a ) % r
      pnode ( inode ) % nneighbor = 0
      nexcl = pnode ( inode ) % nexcl
      IF ( nexcl > SIZE ( elist ) ) STOP 'elist'
      ll_exclude => pnode ( inode ) % ex
! list number of exclusion
      DO iexclude = 1, nexcl
         elist ( iexclude ) = ll_exclude % p % iatom
         ll_exclude => ll_exclude % next
      END DO
! sort elist
      CALL sort ( elist, nexcl, ind )
! allocate neighbor list
      cell_a = (/0,0,0/)
      CALL allocate_neighbor_list_set ( neighbor_list_set = &
                                        sab_orb ( inode ) % neighbor_list_set, &
                                        r_max = 0._dbl )
                                        
      CALL add_neighbor_list( neighbor_list_set =  &
                              sab_orb ( inode ) % neighbor_list_set, &
                              atom = atom_a,&
                              cell = cell_a,&
                              neighbor_list = pnode ( inode ) % nl )
! loop over kinds
      DO jkind = 1, nkinds
        rab2_max = rlist_cutsq ( ikind, jkind )
        sab_max => ll_cell ( ikind, jkind ) % sab_max
        ncell => ll_cell ( ikind, jkind ) % ncell
        ngrid => ll_cell ( ikind, jkind ) % ngrid
        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,ngrid,&
            scaled_to_real(1.0_dbl/REAL(ngrid(:),dbl),cell)/unit_of_length
        END IF

!   *** Check, if we have to consider a subcell grid ***
        IF (SUM(ngrid) == 3) THEN
          DO jatom=1,natom(jkind)
            atom_b = atoms(jkind)%list(jatom)
            rb = part ( atom_b ) % r
            IF ( .NOT. include_list ( atom_a, atom_b ) ) CYCLE
! look for a match in the exclusion list 
            SELECT CASE ( nexcl )
            CASE DEFAULT
              IF ( .NOT. get_match ( atom_b, elist, nexcl ) ) THEN
                  CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                              rab2_max, pnode ( inode ) % nneighbor, &
                                              pnode ( inode ) % nl )
              END IF
            CASE ( 0 )
              CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                          rab2_max, pnode ( inode ) % nneighbor, &
                                          pnode ( inode ) % nl )
            CASE ( 1 )
              IF ( atom_b /= elist ( 1 ) ) THEN
                CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                            rab2_max, pnode ( inode ) % nneighbor, &
                                            pnode ( inode ) % nl )
              END IF
            END SELECT
          ENDDO

        ELSE

          sa_pbc(:) = pbc_coord % s ( :, atom_a )
          sa_min(:) = sa_pbc(:) - sab_max(:)
          sa_max(:) = sa_pbc(:) + sab_max(:)

          grid_min => ll_cell ( ikind, jkind ) % grid_min
          grid_max => ll_cell ( ikind, jkind ) % grid_max
          nijk => ll_cell ( ikind, jkind ) % nijk
          ijk => ll_cell ( ikind, jkind ) % ijk

          DO igrid = 1, ngrid ( 1 )
            DO jgrid = 1, ngrid ( 2 )
              DO kgrid = 1, ngrid ( 3 )

                wrap_x = .FALSE.
                wrap_y = .FALSE.
                wrap_z = .FALSE.

                IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                IF (sa_max(1) > 0.5_dbl) THEN
                  sa_max_pbc(1)=sa_max(1)-periodic(1)*INT(sa_max(1)+SIGN(0.5_dbl,sa_max(1)))
                  IF (sa_max_pbc(1) > grid_min(1, igrid, jgrid, kgrid )) &
                  wrap_x = .TRUE.
                ELSEIF (sa_min(1) < -0.5_dbl) THEN
                  sa_min_pbc(1)=sa_min(1)-periodic(1)*INT(sa_min(1)+SIGN(0.5_dbl,sa_min(1)))
                  IF (sa_min_pbc(1) <= grid_max(1, igrid, jgrid, kgrid )) & 
                  wrap_x=.TRUE.
                ENDIF

                IF (.NOT.wrap_x) THEN
                  IF (sa_max(1) < grid_min(1, igrid, jgrid, kgrid )) CYCLE
                  IF (sa_min(1) >= grid_max(1, igrid, jgrid, kgrid )) CYCLE
                ENDIF

                IF (sa_max(2) > 0.5_dbl) THEN
                  sa_max_pbc(2)=sa_max(2)-periodic(2)*INT(sa_max(2)+SIGN(0.5_dbl,sa_max(2)))
                  IF (sa_max_pbc(2) > grid_min(2, igrid, jgrid, kgrid )) &
                  wrap_y=.TRUE.
                ELSEIF (sa_min(2) < -0.5_dbl) THEN
                  sa_min_pbc(2)=sa_min(2)-periodic(2)*INT(sa_min(2)+SIGN(0.5_dbl,sa_min(2)))
                  IF (sa_min_pbc(2) <= grid_max(2, igrid, jgrid, kgrid )) &
                  wrap_y=.TRUE.
                ENDIF

                IF (.NOT.wrap_y) THEN
                  IF (sa_max(2) < grid_min(2, igrid, jgrid, kgrid )) CYCLE
                  IF (sa_min(2) >= grid_max(2, igrid, jgrid, kgrid )) CYCLE
                ENDIF

                IF (sa_max(3) > 0.5_dbl) THEN
                  sa_max_pbc(3)=sa_max(3)-periodic(3)*INT(sa_max(3)+SIGN(0.5_dbl,sa_max(3)))
                  IF (sa_max_pbc(3) > grid_min(3, igrid, jgrid, kgrid )) &
                  wrap_z=.TRUE.
                ELSEIF (sa_min(3) < -0.5_dbl) THEN
                  sa_min_pbc(3)=sa_min(3)-periodic(3)*INT(sa_min(3)+SIGN(0.5_dbl,sa_min(3)))
                  IF (sa_min_pbc(3) <= grid_max(3, igrid, jgrid, kgrid )) &
                  wrap_z=.TRUE.
                ENDIF

                IF (.NOT.wrap_z) THEN
                  IF (sa_max(3) < grid_min(3, igrid, jgrid, kgrid )) CYCLE
                  IF (sa_min(3) >= grid_max(3, igrid, jgrid, kgrid )) CYCLE
                ENDIF


                DO iijk = 1, nijk ( igrid, jgrid, kgrid )
                  atom_b = ijk ( iijk, igrid, jgrid, kgrid )
                  rb = part ( atom_b ) % r
                  IF ( .NOT. include_list ( atom_a, atom_b ) ) CYCLE
! look for a match in the exclusion list 
                  SELECT CASE ( nexcl )
                  CASE DEFAULT
                    IF ( .NOT. get_match ( atom_b, elist, nexcl ) ) THEN
                      CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                                  rab2_max, pnode ( inode ) % nneighbor, &
                                                  pnode ( inode ) % nl )
                    END IF
                  CASE ( 0 )
                    CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                                rab2_max, pnode ( inode ) % nneighbor, &
                                                pnode ( inode ) % nl )
                  CASE ( 1 )
                    IF ( atom_b /= elist ( 1 ) ) THEN
                      CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                                  rab2_max, pnode ( inode ) % nneighbor, &
                                                  pnode ( inode ) % nl )
                    END IF
                  END SELECT
                END DO

              END DO
            END DO
          END DO

        END IF

      END DO
    END DO

    DEALLOCATE ( elist, STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'elist' )

    DEALLOCATE (ind,STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'ind' )

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE build_verlet_neighbor_lists

!******************************************************************************

! to get the Verlet neighbor list for nonbonded interactions

  SUBROUTINE rebuild_verlet_lists ( atomic_kind_set, part, cell, ll_data, pnode )

!   Purpose: Build all the required neighbor lists FIST

!   History: - Creation (28.08.2000,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), POINTER                 :: cell
    TYPE(linklist_internal_data_type), &
      INTENT(IN)                             :: ll_data
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode

    CHARACTER(LEN=*), PARAMETER :: module_name = "linklist_verlet_list", &
      routine_name = "rebuild_verlet_lists", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER :: atom_a, group, handle, ikind, inode, istat, mype, natoms, &
      nkind, nnodes, npe, output_unit, periodic( 3 )
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom
    LOGICAL                                  :: ionode
    REAL(dbl)                                :: lat_vec( 3 ), ra( 3 ), sa( 3 )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atoms
    TYPE(global_environment_type)            :: globenv
    TYPE(pbc_coord_type)                     :: pbc_coord

!   ---------------------------------------------------------------------------

    globenv = ll_data % globenv

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)
  

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    IF ((ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > LOW)) THEN
      print_cell_parameters = .TRUE.
    ELSE
      print_cell_parameters = .FALSE.
    END IF

    CALL get_cell(cell=cell, &
                  periodic=periodic)

!   *** Allocate work storage ***

    nkind = SIZE ( atomic_kind_set )

    ALLOCATE ( natom ( nkind ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( routine, "natom", nkind * int_size )

    ALLOCATE ( atoms ( nkind ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( routine, "atoms", nkind * int_size )

!   *** fill atomic kind info ***

    DO ikind = 1, nkind

      atomic_kind => atomic_kind_set ( ikind )

      NULLIFY ( atoms ( ikind ) % list )

      CALL get_atomic_kind ( atomic_kind = atomic_kind, &
                             natom = natom ( ikind ), &
                             atom_list = atoms ( ikind ) % list )
    END DO

!   *** Calculate PBC coordinates ***
    natoms = SIZE ( part )
    ALLOCATE ( pbc_coord % s ( 3, natoms ), STAT = istat )
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s",&
                                     3*natoms*dp_size)
    nnodes = SIZE ( pnode )
    DO inode=1,nnodes
      atom_a = pnode ( inode ) % p % iatom
      ra = pnode ( inode ) % p % r
      sa  = real_to_scaled ( ra ( : ), cell )
      lat_vec ( 1 ) = periodic(1)*INT(sa(1)+SIGN(0.5_dbl,sa(1)))
      lat_vec ( 2 ) = periodic(2)*INT(sa(2)+SIGN(0.5_dbl,sa(2)))
      lat_vec ( 3 ) = periodic(3)*INT(sa(3)+SIGN(0.5_dbl,sa(3)))
      pbc_coord % s ( :, atom_a ) = sa ( : ) - lat_vec ( : )
    END DO

!   *** Rebuild the verlet neighbor lists ***

    CALL rebuild_verlet_neighbor_lists( ll_data % rlist_cutsq, pbc_coord, atoms, &
                                        natom, ll_data % sab_ppnl, cell,         &
                                        ll_data % ll_cell, part, pnode )

!   *** Release work storage ***

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      NULLIFY (atoms(ikind)%list)
    ENDDO

    DEALLOCATE (pbc_coord%s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s")

    DEALLOCATE (atoms,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atoms")

    CALL timestop(0.0_dbl,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE rebuild_verlet_lists

! *****************************************************************************

  SUBROUTINE rebuild_verlet_neighbor_lists( rlist_cutsq, pbc_coord, atoms,  &
                                            natom, sab_orb, cell, ll_cell, part, pnode )

!   Purpose: Build the neighbor lists for calculation of the orbital-orbital
!            overlap matrix.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************
    REAL(dbl), INTENT(IN)                    :: rlist_cutsq( :, : )
    TYPE(pbc_coord_type), INTENT(IN)         :: pbc_coord
    TYPE(atoms_type), INTENT(IN)             :: atoms( : )
    INTEGER, INTENT(IN)                      :: natom( : )
    TYPE(neighbor_list_set_p_type), POINTER  :: sab_orb( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(subcell_data_type), POINTER         :: ll_cell( :, : )
    TYPE(particle_type), INTENT(IN)          :: part( : )
    TYPE(particle_node_type), INTENT(INOUT)  :: pnode( : )

    CHARACTER(LEN=*), PARAMETER :: &
      routine_name = "rebuild_verlet_neighbor_lists", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"
    INTEGER, PARAMETER                       :: els = 50

    INTEGER :: atom_a, atom_b, cell_a( 3 ), handle, iatom, iexclude, igrid, &
      iijk, ikind, inode, istat, jatom, jgrid, jkind, kgrid, nexcl, nkinds, &
      nnode, periodic( 3 )
    INTEGER, DIMENSION(:), POINTER           :: elist, ind, ncell, ngrid
    INTEGER, POINTER                         :: ijk( :, :, :, : ), &
                                                nijk( :, :, : )
    LOGICAL                                  :: wrap_x, wrap_y, wrap_z
    REAL(dbl) :: r( 3 ), ra( 3 ), rab2_max, rb( 3 ), sa_max( 3 ), &
      sa_max_pbc( 3 ), sa_min( 3 ), sa_min_pbc( 3 ), sa_pbc( 3 ), sb_pbc( 3 )
    REAL(dbl), POINTER                       :: grid_max( :, :, :, : ), &
                                                grid_min( :, :, :, : ), &
                                                sab_max( : )
    TYPE(linklist_exclusion), POINTER        :: ll_exclude

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    CALL get_cell(cell=cell, periodic=periodic )
    nkinds = SIZE ( ll_cell, 1 )
! Check for ncells and ngrids
    DO ikind = 1, nkinds
      DO jkind = 1, nkinds
        ngrid => ll_cell ( ikind, jkind ) % ngrid
        IF ( SUM ( ngrid ) /= 3 ) THEN
          nijk => ll_cell ( ikind, jkind ) % nijk
          nijk = 0
          ijk => ll_cell ( ikind, jkind ) % ijk
          ijk = 0
          DO jatom=1,natom(jkind)
            atom_b = atoms(jkind)%list(jatom)
            sb_pbc(:) = pbc_coord % s ( :, atom_b ) 
            igrid = MAX(1,CEILING((sb_pbc(1) + 0.5_dbl)*ngrid(1)))
            jgrid = MAX(1,CEILING((sb_pbc(2) + 0.5_dbl)*ngrid(2)))
            kgrid = MAX(1,CEILING((sb_pbc(3) + 0.5_dbl)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = atom_b
          END DO
        END IF 
      END DO
    END DO

    ALLOCATE ( elist ( els ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'elist', els )

    ALLOCATE ( ind ( els ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'ind', els )

    nnode = SIZE ( pnode )

!   *** Loop over all nodes ***
    DO inode = 1, nnode
      CALL get_atomic_kind ( pnode ( inode ) % p % atomic_kind, &
                             kind_number = ikind )
      atom_a = pnode ( inode ) % p % iatom
      ra = part ( atom_a ) % r
      pnode ( inode ) % nneighbor = 0
      nexcl = pnode ( inode ) % nexcl
      IF ( nexcl > SIZE ( elist ) ) STOP 'elist'
      ll_exclude => pnode ( inode ) % ex
! list number of exclusion
      DO iexclude = 1, nexcl
         elist ( iexclude ) = ll_exclude % p % iatom
         ll_exclude => ll_exclude % next
      END DO
! sort elist
      CALL sort ( elist, nexcl, ind )
! allocate neighbor list
      cell_a = (/0,0,0/)
      CALL get_neighbor_list(neighbor_list=pnode(inode)%nl,&
                             atom=atom_a, &
                             cell=cell_a )
      CALL init_neighbor_list(pnode(inode)%nl)
                                        
! loop over kinds
      DO jkind = 1, nkinds
        rab2_max = rlist_cutsq ( ikind, jkind )
        sab_max => ll_cell ( ikind, jkind ) % sab_max
        ncell => ll_cell ( ikind, jkind ) % ncell
        ngrid => ll_cell ( ikind, jkind ) % ngrid

!   *** Check, if we have to consider a subcell grid ***
        IF (SUM(ngrid) == 3) THEN
          DO jatom=1,natom(jkind)
            atom_b = atoms(jkind)%list(jatom)
            rb = part ( atom_b ) % r
            IF ( .NOT. include_list ( atom_a, atom_b ) ) CYCLE
! look for a match in the exclusion list 
            SELECT CASE ( nexcl )
            CASE DEFAULT
              IF ( .NOT. get_match ( atom_b, elist, nexcl ) ) THEN
                CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                            rab2_max, pnode ( inode ) % nneighbor, &
                                            pnode ( inode ) % nl )
              END IF
            CASE ( 0 )
                CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                            rab2_max, pnode ( inode ) % nneighbor, &
                                            pnode ( inode ) % nl )
            CASE ( 1 )
              IF ( atom_b /= elist ( 1 ) ) THEN
              END IF
                CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                            rab2_max, pnode ( inode ) % nneighbor, &
                                            pnode ( inode ) % nl )
            END SELECT
          ENDDO

        ELSE

          sa_pbc(:) = pbc_coord % s ( :, atom_a ) 
          sa_min(:) = sa_pbc(:) - sab_max(:)
          sa_max(:) = sa_pbc(:) + sab_max(:)

          grid_min => ll_cell ( ikind, jkind ) % grid_min
          grid_max => ll_cell ( ikind, jkind ) % grid_max
          nijk => ll_cell ( ikind, jkind ) % nijk
          ijk => ll_cell ( ikind, jkind ) % ijk

          DO igrid = 1, ngrid ( 1 )
            DO jgrid = 1, ngrid ( 2 )
              DO kgrid = 1, ngrid ( 3 )

                wrap_x = .FALSE.
                wrap_y = .FALSE.
                wrap_z = .FALSE.

                IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                IF (sa_max(1) > 0.5_dbl) THEN
                  sa_max_pbc(1)=sa_max(1)-periodic(1)*INT(sa_max(1)+SIGN(0.5_dbl,sa_max(1)))
                  IF (sa_max_pbc(1) > grid_min(1, igrid, jgrid, kgrid )) &
                  wrap_x = .TRUE.
                ELSEIF (sa_min(1) < -0.5_dbl) THEN
                  sa_min_pbc(1)=sa_min(1)-periodic(1)*INT(sa_min(1)+SIGN(0.5_dbl,sa_min(1)))
                  IF (sa_min_pbc(1) <= grid_max(1, igrid, jgrid, kgrid )) & 
                  wrap_x=.TRUE.
                ENDIF

                IF (.NOT.wrap_x) THEN
                  IF (sa_max(1) < grid_min(1, igrid, jgrid, kgrid )) CYCLE
                  IF (sa_min(1) >= grid_max(1, igrid, jgrid, kgrid )) CYCLE
                ENDIF

                IF (sa_max(2) > 0.5_dbl) THEN
                  sa_max_pbc(2)=sa_max(2)-periodic(2)*INT(sa_max(2)+SIGN(0.5_dbl,sa_max(2)))
                  IF (sa_max_pbc(2) > grid_min(2, igrid, jgrid, kgrid )) &
                  wrap_y=.TRUE.
                ELSEIF (sa_min(2) < -0.5_dbl) THEN
                  sa_min_pbc(2)=sa_min(2)-periodic(2)*INT(sa_min(2)+SIGN(0.5_dbl,sa_min(2)))
                  IF (sa_min_pbc(2) <= grid_max(2, igrid, jgrid, kgrid )) &
                  wrap_y=.TRUE.
                ENDIF

                IF (.NOT.wrap_y) THEN
                  IF (sa_max(2) < grid_min(2, igrid, jgrid, kgrid )) CYCLE
                  IF (sa_min(2) >= grid_max(2, igrid, jgrid, kgrid )) CYCLE
                ENDIF

                IF (sa_max(3) > 0.5_dbl) THEN
                  sa_max_pbc(3)=sa_max(3)-periodic(3)*INT(sa_max(3)+SIGN(0.5_dbl,sa_max(3)))
                  IF (sa_max_pbc(3) > grid_min(3, igrid, jgrid, kgrid )) &
                  wrap_z=.TRUE.
                ELSEIF (sa_min(3) < -0.5_dbl) THEN
                  sa_min_pbc(3)=sa_min(3)-periodic(3)*INT(sa_min(3)+SIGN(0.5_dbl,sa_min(3)))
                  IF (sa_min_pbc(3) <= grid_max(3, igrid, jgrid, kgrid )) &
                  wrap_z=.TRUE.
                ENDIF

                IF (.NOT.wrap_z) THEN
                  IF (sa_max(3) < grid_min(3, igrid, jgrid, kgrid )) CYCLE
                  IF (sa_min(3) >= grid_max(3, igrid, jgrid, kgrid )) CYCLE
                ENDIF


                DO iijk = 1, nijk ( igrid, jgrid, kgrid )
                  atom_b = ijk ( iijk, igrid, jgrid, kgrid )
                  rb = part ( atom_b ) % r
                  IF ( .NOT. include_list ( atom_a, atom_b ) ) CYCLE
! look for a match in the exclusion list 
                  SELECT CASE ( nexcl )
                  CASE DEFAULT
                    IF ( .NOT. get_match ( atom_b, elist, nexcl ) ) THEN
                      CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                                  rab2_max, pnode ( inode ) % nneighbor, &
                                                  pnode ( inode ) % nl )
                    END IF
                  CASE ( 0 )
                      CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                                  rab2_max, pnode ( inode ) % nneighbor, &
                                                  pnode ( inode ) % nl )
                  CASE ( 1 )
                    IF ( atom_b /= elist ( 1 ) ) THEN
                      CALL update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                                  rab2_max, pnode ( inode ) % nneighbor, &
                                                  pnode ( inode ) % nl )
                    END IF
                  END SELECT
                END DO

              END DO
            END DO
          END DO

        END IF

      END DO
    END DO

    DEALLOCATE ( elist, STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'elist' )

    DEALLOCATE (ind,STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ( 'verlet_list', 'ind' )

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE rebuild_verlet_neighbor_lists

! *****************************************************************************
  SUBROUTINE update_neighbor_node ( cell, ra, rb, atom_b, ncell, periodic, &
                                    rab2_max, nneighbor, neighbor_list )
    TYPE(cell_type), POINTER                 :: cell
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: ra, rb
    INTEGER, INTENT(IN)                      :: atom_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: ncell, periodic
    REAL(dbl), INTENT(IN)                    :: rab2_max
    INTEGER, INTENT(INOUT)                   :: nneighbor
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    INTEGER                                  :: cell_b( 3 ), icell, imax, &
                                                imin, jcell, jmax, jmin, &
                                                kcell, kmax, kmin, &
                                                lat_vec( 3 )
    REAL(dbl)                                :: rab( 3 ), rab2, sab( 3 ), &
                                                sab_pbc( 3 )

! compute distance vector

    rab ( : ) = rb( : ) - ra ( : )
    sab ( : ) = real_to_scaled(rab(:),cell)
! compute lattice translation vector
    lat_vec ( 1 ) = periodic(1)*INT(sab(1)+SIGN(0.5_dbl,sab(1)))
    lat_vec ( 2 ) = periodic(2)*INT(sab(2)+SIGN(0.5_dbl,sab(2)))
    lat_vec ( 3 ) = periodic(3)*INT(sab(3)+SIGN(0.5_dbl,sab(3)))

! now do lattice translations of  minimum image
    imax = ncell ( 1 ) - lat_vec ( 1 )
    imin = -ncell ( 1 ) - lat_vec ( 1 )
    jmax = ncell ( 2 ) - lat_vec ( 2 )
    jmin = -ncell ( 2 ) - lat_vec ( 2 )
    kmax = ncell ( 3 ) - lat_vec ( 3 )
    kmin = -ncell ( 3 ) - lat_vec ( 3 )

    DO icell=imin, imax
      cell_b(1) = icell
      DO jcell=jmin, jmax
        cell_b(2) = jcell
        DO kcell=kmin, kmax
          cell_b(3) = kcell

          sab_pbc(:) = sab(:) + REAL(cell_b(:),dbl)
          rab(:) = scaled_to_real(sab_pbc(:),cell)
          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

          IF (rab2 < rab2_max .AND. rab2 > 1.e-6_dbl ) THEN
            nneighbor = nneighbor + 1
            CALL add_neighbor_node(&
              neighbor_list=neighbor_list,&
              neighbor = atom_b,&
              cell = cell_b,&
              r = rab ( : ) )
          END IF

        END DO
      END DO
    END DO
  END SUBROUTINE update_neighbor_node
! *****************************************************************************

  SUBROUTINE write_neighbor_lists(part,cell,neighbor_lists,unit_of_length_name, &
                                  unit_of_length,name)

!   Purpose: Write a set of neighbor lists to the output unit.

!   History: - Creation (04.03.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), POINTER                 :: cell
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: neighbor_lists
    CHARACTER(LEN=*), INTENT(IN)             :: unit_of_length_name
    REAL(dbl), INTENT(IN)                    :: unit_of_length
    CHARACTER(LEN=*), INTENT(IN)             :: name

    INTEGER                                  :: ab, atom_a, atom_b, i, &
                                                nneighbor
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b
    REAL(dbl), DIMENSION(3)                  :: ra, rab, rb
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node

!   ---------------------------------------------------------------------------
!   *** Write headline ***

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,7X,A,2(11X,A),10X,A)")&
      TRIM(name)//" IN "//TRIM(unit_of_length_name),&
      "Atom  Neighbor  Cell(i,j,k)","X","Y","Z","Distance"

    DO ab=1,SIZE(neighbor_lists)

      neighbor_list_set => neighbor_lists(ab)%neighbor_list_set

      IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

!     *** Loop over all atoms and their corresponding neighbor lists ***

      neighbor_list => first_list(neighbor_list_set)

      DO WHILE (ASSOCIATED(neighbor_list))

        CALL get_neighbor_list(neighbor_list=neighbor_list,&
                               atom=atom_a,&
                               cell=cell_a,&
                               nnode=nneighbor)

        ra(:) = pbc(part(atom_a)%r,cell,cell_a)

        WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,I6,2X,3I4,3F12.6)")&
          atom_a,nneighbor,cell_a(:),ra(:)/unit_of_length

!       *** Direct the work pointer to the  ***
!       *** start point of the current list ***

        neighbor_node => first_node(neighbor_list)

!       *** Traverse the neighbor list of the current ***
!       *** atom and print the stored information     ***

        DO WHILE (ASSOCIATED(neighbor_node))

          CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                 neighbor=atom_b,&
                                 cell=cell_b,&
                                 r=rab)

          rb(:) = ra(:) + rab(:)

          WRITE (UNIT=output_unit,FMT="(T10,I6,2X,3I4,3F12.6,2X,F12.6)")&
            atom_b,cell_b,rb(:)/unit_of_length,&
            SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))/unit_of_length

          neighbor_node => next(neighbor_node)

        END DO

        neighbor_list => next(neighbor_list)

      END DO

    END DO

  END SUBROUTINE write_neighbor_lists

! *****************************************************************************
END MODULE linklist_verlet_list
