!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_neighbor_lists [1.0] *
!!
!!   NAME
!!     qs_neighbor_lists
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists for Quickstep (QS).
!!
!!   AUTHOR
!!     Matthias Krack (08.10.1999,26.03.2002)
!!
!!   MODIFICATION HISTORY
!!     - List update for sab_orb neighbor list (10.09.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_neighbor_lists

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE basis_set_types,          ONLY: get_gto_basis_set,&
                                      gto_basis_set_type
  USE external_potential_types, ONLY: get_potential,&
                                      gth_potential_type
  USE global_types,             ONLY: LOW,global_environment_type
  USE message_passing,          ONLY: mp_sync
  USE particle_types,           ONLY: particle_type
  USE qs_environment_types,     ONLY: get_qs_env,&
                                      qs_environment_type,&
                                      set_qs_env
  USE qs_neighbor_list_types,   ONLY: add_neighbor_list,&
                                      add_neighbor_node,&
                                      allocate_neighbor_list_set,&
                                      clean_neighbor_list_set,&
                                      deallocate_neighbor_list_set,&
                                      find_neighbor_list,&
                                      first_list,&
                                      first_node,&
                                      get_neighbor_list,&
                                      get_neighbor_list_set,&
                                      get_neighbor_node,&
                                      init_neighbor_list,&
                                      init_neighbor_list_set,&
                                      neighbor_list_p_type,&
                                      neighbor_list_set_p_type,&
                                      neighbor_list_set_type,&
                                      neighbor_list_type,&
                                      neighbor_node_type,&
                                      next
  USE simulation_cell,          ONLY: cell_type,&
                                      get_cell,&
                                      pbc,&
                                      real_to_scaled,&
                                      scaled_to_real
  USE string_utilities,         ONLY: compress
  USE termination,              ONLY: stop_memory,&
                                      stop_program
  USE timings,                  ONLY: timeset,&
                                      timestop

  IMPLICIT NONE

  PRIVATE

! *** Global types of the module ***

  TYPE atoms_type
    INTEGER, DIMENSION(:), POINTER :: list
  END TYPE atoms_type

  TYPE pbc_coord_type
    REAL(wp), DIMENSION(:,:), POINTER :: r,s
  END TYPE pbc_coord_type

! *** Global variables of the module ***

  TYPE(cell_type), POINTER              :: cell
  TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set,&
                                           sca_neighbor_list_set,&
                                           scb_neighbor_list_set
  TYPE(neighbor_list_type), POINTER     :: neighbor_list,&
                                           sca_neighbor_list,&
                                           scb_neighbor_list
  TYPE(neighbor_node_type), POINTER     :: neighbor_node,&
                                           sca_neighbor_node,&
                                           scb_neighbor_node
  CHARACTER(LEN=40)                     :: string
  CHARACTER(LEN=8)                      :: unit_of_length_name
  REAL(wp)                              :: a_max,a_min,&
                                           b_max,b_min,&
                                           c_max,c_min,&
                                           rab_max,rab2,rab2_max,&
                                           rac_max,rac2,rac2_max,&
                                           rca_max,rca2,rca2_max,&
                                           rcb_max,rcb2,rcb2_max,&
                                           subcells,unit_of_length
  INTEGER                               :: atom_a,atom_b,atom_c,group,&
                                           iatom,icell,igrid,iijk,ikind,&
                                           ineighbor,ipe,istat,&
                                           jatom,jcell,jgrid,jkind,&
                                           katom,kcell,kgrid,kkind,&
                                           maxatom,mype,&
                                           nkind,nneighbor,npe,output_unit
  LOGICAL                               :: gth_potential_present,ionode,&
                                           print_cell_parameters

  REAL(wp), DIMENSION(3) :: r,ra_pbc,rab,rac,rb,rc,rc_pbc,s,s_max,s_min,sab,&
                            sab_max,sac_max,sb,sb_max,sb_min,sb_pbc,sc,sc_max,&
                            sc_min,sa_pbc,sc_pbc,sca_max,scb_max
  INTEGER, DIMENSION(3)  :: cell_a,cell_b,cell_c,ncell,ngrid,periodic

  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_ppl,&
                                                           sca_ppnl,scb_ppnl
  TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

  TYPE(atoms_type), DIMENSION(:), ALLOCATABLE           :: atoms
  TYPE(neighbor_list_p_type), DIMENSION(:), ALLOCATABLE :: kind_a
  TYPE(pbc_coord_type), DIMENSION(:), ALLOCATABLE       :: pbc_coord
  REAL(wp), DIMENSION(:), ALLOCATABLE                   :: orb_radius,&
                                                           ppl_radius,&
                                                           ppnl_radius
  INTEGER, DIMENSION(:), ALLOCATABLE                    :: atom_of_kind,&
                                                           natom,neighbors
  LOGICAL, DIMENSION(:), ALLOCATABLE                    :: orb_present,&
                                                           ppl_present,&
                                                           ppnl_present
  REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE             :: grid_max,grid_min
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE                :: nijk
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE              :: ijk

! *** Public subroutines ***

  PUBLIC :: build_neighbor_lists,&
            update_neighbor_lists

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_neighbor_lists(qs_env,globenv)

!   Purpose: Build all the required neighbor lists for Quickstep.

!   History: - Creation (28.08.2000,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    INTEGER                           :: handle

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    CALL timeset("build_neighbor_lists","I","",handle)

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

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sca_ppnl=sca_ppnl,&
                    scb_ppnl=scb_ppnl)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxatom=maxatom)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Allocate work storage ***

    ALLOCATE (kind_a(maxatom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a",maxatom*int_size)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (natom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom",nkind*int_size)

    ALLOCATE (atoms(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atoms",nkind*int_size)

    ALLOCATE (pbc_coord(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord",nkind*int_size)

    ALLOCATE (orb_present(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"orb_present",nkind*int_size)

    ALLOCATE (orb_radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"orb_radius",nkind*wp_size)
    orb_radius = 0.0_wp

    IF (gth_potential_present) THEN

      ALLOCATE (ppl_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_present",nkind*int_size)

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius",nkind*wp_size)
      ppl_radius = 0.0_wp

      ALLOCATE (ppnl_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppnl_present",nkind*int_size)

      ALLOCATE (ppnl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppnl_radius",nkind*wp_size)
      ppnl_radius = 0.0_wp

    END IF

!   *** Calculate PBC coordinates ***

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atoms(ikind)%list)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           gth_potential=gth_potential,&
                           orb_basis_set=orb_basis_set,&
                           natom=natom(ikind),&
                           atom_list=atoms(ikind)%list)

      IF (ASSOCIATED(gth_potential)) THEN
        CALL get_potential(potential=gth_potential,&
                           ppl_present=ppl_present(ikind),&
                           ppl_radius=ppl_radius(ikind),&
                           ppnl_present=ppnl_present(ikind),&
                           ppnl_radius=ppnl_radius(ikind))
      END IF

      IF (ASSOCIATED(orb_basis_set)) THEN

        orb_present(ikind) = .TRUE.

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=orb_radius(ikind))

        ALLOCATE (pbc_coord(ikind)%r(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r",&
                                         3*natom(ikind)*wp_size)
        ALLOCATE (pbc_coord(ikind)%s(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s",&
                                         3*natom(ikind)*wp_size)

        DO iatom=1,natom(ikind)
          atom_a = atoms(ikind)%list(iatom)
          ra_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          pbc_coord(ikind)%r(:,iatom) = ra_pbc(:)
          pbc_coord(ikind)%s(:,iatom) = real_to_scaled(ra_pbc(:),cell)
        END DO

      ELSE

        orb_present(ikind) = .FALSE.

        NULLIFY (pbc_coord(ikind)%r)
        NULLIFY (pbc_coord(ikind)%s)

      END IF

    END DO

!   *** Build the orbital-orbital overlap neighbor lists ***

    CALL build_sab_orb_neighbor_lists()

!   *** Put the generated neighbor lists into the QS environment ***

    CALL set_qs_env(qs_env=qs_env,&
                    sab_orb=sab_orb)

    IF (ionode.AND.globenv%print%sab_orb_neighbor_lists) THEN
      CALL write_neighbor_lists(sab_orb,"ORBITAL ORBITAL NEIGHBOR LISTS")
    END IF

!   *** Build orbital GTH-PPL operator overlap list ***

    IF (gth_potential_present) THEN

      IF (ANY(ppl_present)) THEN

        CALL build_gth_ppl_neighbor_lists()

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_qs_env(qs_env=qs_env,&
                        sac_ppl=sac_ppl)

        IF (ionode.AND.globenv%print%sac_ppl_neighbor_lists) THEN
          CALL write_neighbor_lists(sac_ppl,"ORBITAL GTH-PPL NEIGHBOR LISTS")
        END IF

      END IF

      IF (ANY(ppnl_present)) THEN

        ALLOCATE (neighbors(maxatom),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbors",maxatom*int_size)

        CALL build_gth_ppnl_neighbor_lists()

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_qs_env(qs_env=qs_env,&
                        sca_ppnl=sca_ppnl,&
                        scb_ppnl=scb_ppnl)

        IF (globenv%print%sca_ppnl_neighbor_lists) THEN
          IF (npe == 1) THEN
            CALL write_neighbor_lists(sca_ppnl,&
                                      "GTH-PPNL ORBITAL(A) NEIGHBOR LISTS")
          ELSE
            DO ipe=0,npe-1
              IF (ipe == mype) THEN
                WRITE (string,"(A,I5)")&
                  "SCA_PPNL NEIGHBOR LISTS OF PROCESS",ipe
                CALL compress(string)
                CALL write_neighbor_lists(sca_ppnl,string)
              END IF
              CALL mp_sync(group)
            END DO
          END IF
        END IF

        IF (globenv%print%scb_ppnl_neighbor_lists) THEN
          IF (npe == 1) THEN
            CALL write_neighbor_lists(scb_ppnl,&
                                      "GTH-PPNL ORBITAL(B) NEIGHBOR LISTS")
          ELSE
            DO ipe=0,npe-1
              IF (ipe == mype) THEN
                WRITE (string,"(A,I5)")&
                  "SCB_PPNL NEIGHBOR LISTS OF PROCESS",ipe
                CALL compress(string)
                CALL write_neighbor_lists(scb_ppnl,string)
              END IF
              CALL mp_sync(group)
            END DO
          END IF
        END IF

      END IF

    END IF

!   *** Release work storage ***

    DEALLOCATE (kind_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a")

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      NULLIFY (atoms(ikind)%list)
      IF (ASSOCIATED(pbc_coord(ikind)%r)) THEN
        DEALLOCATE (pbc_coord(ikind)%r,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r")
      END IF
      IF (ASSOCIATED(pbc_coord(ikind)%s)) THEN
        DEALLOCATE (pbc_coord(ikind)%s,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s")
      END IF
    END DO

    DEALLOCATE (atoms,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atoms")

    DEALLOCATE (pbc_coord,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord")

    DEALLOCATE (orb_present,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"orb_present")

    DEALLOCATE (orb_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"orb_radius")

    IF (gth_potential_present) THEN

      DEALLOCATE (ppl_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_present")

      DEALLOCATE (ppl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius")

      DEALLOCATE (ppnl_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppnl_present")

      DEALLOCATE (ppnl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppnl_radius")

      IF (ALLOCATED(neighbors)) THEN
        DEALLOCATE (neighbors,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbors")
      END IF

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_sab_orb_neighbor_lists()

!   Purpose: Build the neighbor lists for calculation of the orbital-orbital
!            overlap matrix.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_sab_orb_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    INTEGER :: ab,handle,natom_a
    LOGICAL :: equal_kinds

!   ---------------------------------------------------------------------------

    CALL timeset("build_sab_orb_neighbor_lists","I","",handle)

    IF (ASSOCIATED(sab_orb)) THEN
      DO ab=1,SIZE(sab_orb)
        CALL deallocate_neighbor_list_set(sab_orb(ab)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sab_orb")
    END IF

    ALLOCATE (sab_orb(nkind*(nkind + 1)/2),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sab_orb",0)

    DO ab=1,SIZE(sab_orb)
      NULLIFY (sab_orb(ab)%neighbor_list_set)
    END DO

!   *** Print headline ***

    IF (print_cell_parameters) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE ORBITAL ORBITAL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.orb_present(ikind)) CYCLE

      DO jkind=ikind,nkind

        IF (.NOT.orb_present(jkind)) CYCLE

        ab = ikind + jkind*(jkind - 1)/2

        equal_kinds = (ikind == jkind)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = orb_radius(ikind) + orb_radius(jkind)
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,ngrid,&
            scaled_to_real(1.0_wp/REAL(ngrid(:),wp),cell)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          sab_orb(ab)%neighbor_list_set,&
                                        r_max=rab_max)
        neighbor_list_set => sab_orb(ab)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom)%neighbor_list)
          END DO

          DO jatom=1,natom(jkind)

            atom_b = atoms(jkind)%list(jatom)
            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),wp)
                  rb(:) = scaled_to_real(sb(:),cell)

                  IF (equal_kinds) THEN
                    natom_a = jatom
                  ELSE
                    natom_a = natom(ikind)
                  END IF

                  DO iatom=1,natom_a
                    rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                    IF (rab2 < rab2_max) THEN
                      CALL add_neighbor_node(&
                        neighbor_list=kind_a(iatom)%neighbor_list,&
                        neighbor=atom_b,&
                        cell=cell_b,&
                        r=rab(:))
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(ikind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",&
                             natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
              b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = iatom
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom)%neighbor_list)
          END DO

          DO jatom=1,natom(jkind)

            atom_b = atoms(jkind)%list(jatom)
            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),wp)
                  rb(:) = scaled_to_real(sb(:),cell)
                  sb_min(:) = sb(:) - sab_max(:)
                  sb_max(:) = sb(:) + sab_max(:)

                  IF (sb_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sb_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sb_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sb_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sb_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sb_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sb_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sb_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          iatom = ijk(iijk,igrid,jgrid,kgrid)
                          IF (equal_kinds) THEN
                            IF (jatom < iatom) CYCLE
                          END IF
                          rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                          IF (rab2 < rab2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_a(iatom)%neighbor_list,&
                              neighbor=atom_b,&
                              cell=cell_b,&
                              r=rab(:))
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

      END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_sab_orb_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_gth_ppl_neighbor_lists()

!   Purpose: Build the neighbor lists for the calculation of the local part of
!            the GTH pseudo potential (PPL).

!   History: - Creation (18.03.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_gth_ppl_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    INTEGER :: ac,handle

!   ---------------------------------------------------------------------------

    CALL timeset("build_gth_ppl_neighbor_lists","I","",handle)

    IF (ASSOCIATED(sac_ppl)) THEN
      DO ac=1,SIZE(sac_ppl)
        CALL deallocate_neighbor_list_set(sac_ppl(ac)%neighbor_list_set)
      END DO
      DEALLOCATE (sac_ppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sac_ppl")
    END IF

    ALLOCATE (sac_ppl(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sac_ppl",0)

    DO ac=1,SIZE(sac_ppl)
      NULLIFY (sac_ppl(ac)%neighbor_list_set)
    END DO

!   *** Print headline ***

    IF (print_cell_parameters) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE ORBITAL GTH-PPL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.orb_present(ikind)) CYCLE

      DO kkind=1,nkind

        ac = ikind + (kkind - 1)*nkind

        IF (.NOT.ppl_present(kkind)) CYCLE

!       *** Calculate the square of the maximum interaction distance ***

        rac_max = orb_radius(ikind) + ppl_radius(kkind)
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/sac_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,kkind,ngrid,&
            scaled_to_real(1.0_wp/REAL(ngrid(:),wp),cell)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          sac_ppl(ac)%neighbor_list_set,&
                                        r_max=rac_max)
        neighbor_list_set => sac_ppl(ac)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom)%neighbor_list)
          END DO

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)
            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)

                  DO iatom=1,natom(ikind)
                    rac(:) = rc(:) - pbc_coord(ikind)%r(:,iatom)
                    rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                    IF (rac2 < rac2_max) THEN
                      CALL add_neighbor_node(&
                        neighbor_list=kind_a(iatom)%neighbor_list,&
                        neighbor=atom_c,&
                        cell=cell_c,&
                        r=rac(:))
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(ikind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
              b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = iatom
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom)%neighbor_list)
          END DO

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)
            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)
                  sc_min(:) = sc(:) - sac_max(:)
                  sc_max(:) = sc(:) + sac_max(:)

                  IF (sc_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sc_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sc_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sc_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          iatom = ijk(iijk,igrid,jgrid,kgrid)
                          rac(:) = rc(:) - pbc_coord(ikind)%r(:,iatom)
                          rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                          IF (rac2 < rac2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_a(iatom)%neighbor_list,&
                              neighbor=atom_c,&
                              cell=cell_c,&
                              r=rac(:))
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

      END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_gth_ppl_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_gth_ppnl_neighbor_lists()

!   Purpose: Build the neighbor lists for the calculation of the non-local part
!            of the GTH pseudo potential (PPNL). Atom c is a PPNL operator
!            atom, atom a is an atom in the cell(0,0,0), and atom b is an atom
!            is an atom in any cell(i,j,k).

!   History: - Creation (18.03.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_gth_ppnl_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    INTEGER :: cb,handle,icellb,jcellb,kcellb

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rcb

!   ---------------------------------------------------------------------------

    CALL timeset("build_gth_ppnl_neighbor_lists","I","",handle)

    IF (ASSOCIATED(sca_ppnl)) THEN
      DO cb=1,SIZE(sca_ppnl)
        CALL deallocate_neighbor_list_set(sca_ppnl(cb)%neighbor_list_set)
      END DO
      DEALLOCATE (sca_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sca_ppnl")
    END IF

    ALLOCATE (sca_ppnl(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sca_ppnl",0)

    DO cb=1,SIZE(sca_ppnl)
      NULLIFY (sca_ppnl(cb)%neighbor_list_set)
    END DO

    IF (ASSOCIATED(scb_ppnl)) THEN
      DO cb=1,SIZE(scb_ppnl)
        CALL deallocate_neighbor_list_set(scb_ppnl(cb)%neighbor_list_set)
      END DO
      DEALLOCATE (scb_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"scb_ppnl")
    END IF

    ALLOCATE (scb_ppnl(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scb_ppnl",0)

    DO cb=1,SIZE(scb_ppnl)
      NULLIFY (scb_ppnl(cb)%neighbor_list_set)
    END DO

!   *** Allocate work storage ***

    ALLOCATE (rcb(3,maxatom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rcb",3*maxatom*wp_size)

!   *** Print headline ***

    IF (print_cell_parameters) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE GTH-PPNL ORBTIAL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO kkind=1,nkind

      IF (.NOT.ppnl_present(kkind)) CYCLE

      DO jkind=1,nkind

        IF (.NOT.orb_present(jkind)) CYCLE

        cb = kkind + (jkind - 1)*nkind

!       *** Calculate the square of the maximum interaction distance ***

        rcb_max = ppnl_radius(kkind) + orb_radius(jkind)
        rcb2_max = rcb_max*rcb_max

        r(:) = rcb_max
        scb_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(scb_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/scb_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            kkind,jkind,ngrid,&
            scaled_to_real(1.0_wp/REAL(ngrid(:),wp),cell)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          sca_ppnl(cb)%neighbor_list_set,&
                                        r_max=rcb_max)
        sca_neighbor_list_set => sca_ppnl(kkind)%neighbor_list_set

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          scb_ppnl(cb)%neighbor_list_set,&
                                        r_max=rcb_max)
        scb_neighbor_list_set => scb_ppnl(kkind)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)

!           *** Distribute the GTH-PPNL operator atoms ***

            IF (mype /= MODULO(atom_c,npe)) CYCLE

            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)

                  sca_neighbor_list =>&
                    find_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c)

                  IF (.NOT.ASSOCIATED(sca_neighbor_list)) THEN
                    CALL add_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=sca_neighbor_list)
                    CALL add_neighbor_list(&
                      neighbor_list_set=scb_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=scb_neighbor_list)
                  ELSE
                    scb_neighbor_list =>&
                      find_neighbor_list(&
                        neighbor_list_set=scb_neighbor_list_set,&
                        atom=atom_c,&
                        cell=cell_c)
                  END IF

                  nneighbor = 1

                  DO jatom=1,natom(jkind)
                    rcb(:,nneighbor) = pbc_coord(jkind)%r(:,jatom) - rc(:)
                    rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                           rcb(2,nneighbor)*rcb(2,nneighbor) +&
                           rcb(3,nneighbor)*rcb(3,nneighbor)
                    IF (rcb2 < rcb2_max) THEN
                      neighbors(nneighbor) = atoms(jkind)%list(jatom)
                      nneighbor =  nneighbor + 1
                    END IF
                  END DO

                  IF (nneighbor > 1) THEN
                    DO ineighbor=1,nneighbor-1
                      CALL add_neighbor_node(&
                        neighbor_list=sca_neighbor_list,&
                        neighbor=neighbors(ineighbor),&
                        cell=cell_a,&
                        r=rcb(:,ineighbor))
                    END DO
                  END IF

                  DO icellb=cell_c(1)-ncell(1),cell_c(1)+ncell(1)
                    cell_b(1) = icellb
                    DO jcellb=cell_c(2)-ncell(2),cell_c(2)+ncell(2)
                      cell_b(2) = jcellb
                      DO kcellb=cell_c(3)-ncell(3),cell_c(3)+ncell(3)
                        cell_b(3) = kcellb

                        IF ((icellb == 0).AND.&
                            (jcellb == 0).AND.&
                            (kcellb == 0)) CYCLE

                        s(:) = REAL(cell_b(:),wp)
                        r(:) = scaled_to_real(s(:),cell) - rc(:)

                        nneighbor = 1

                        DO jatom=1,natom(jkind)
                          rcb(:,nneighbor) =&
                            pbc_coord(jkind)%r(:,jatom) + r(:)
                          rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                                 rcb(2,nneighbor)*rcb(2,nneighbor) +&
                                 rcb(3,nneighbor)*rcb(3,nneighbor)
                          IF (rcb2 < rcb2_max) THEN
                            neighbors(nneighbor) = atoms(jkind)%list(jatom)
                            nneighbor =  nneighbor + 1
                          END IF
                        END DO

                        IF (nneighbor > 1) THEN
                          DO ineighbor=1,nneighbor-1
                            CALL add_neighbor_node(&
                              neighbor_list=scb_neighbor_list,&
                              neighbor=neighbors(ineighbor),&
                              cell=cell_b,&
                              r=rcb(:,ineighbor))
                          END DO
                        END IF

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(jkind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",natom(jkind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
              b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO jatom=1,natom(jkind)
            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)
            igrid = MAX(1,CEILING((sb_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sb_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sb_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = jatom
          END DO

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)

!           *** Distribute the GTH-PPNL operator atoms ***

            IF (mype /= MODULO(atom_c,npe)) CYCLE

            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)

                  sc_min(:) = sc(:) - scb_max(:)
                  sc_max(:) = sc(:) + scb_max(:)

                  IF (sc_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sc_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sc_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  sca_neighbor_list =>&
                    find_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c)

                  IF (.NOT.ASSOCIATED(sca_neighbor_list)) THEN
                    CALL add_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=sca_neighbor_list)
                    CALL add_neighbor_list(&
                      neighbor_list_set=scb_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=scb_neighbor_list)
                  ELSE
                    scb_neighbor_list =>&
                      find_neighbor_list(&
                        neighbor_list_set=scb_neighbor_list_set,&
                        atom=atom_c,&
                        cell=cell_c)
                  END IF

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sc_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        nneighbor = 1

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          jatom = ijk(iijk,igrid,jgrid,kgrid)
                          rcb(:,nneighbor) =&
                            pbc_coord(jkind)%r(:,jatom) - rc(:)
                          rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                                 rcb(2,nneighbor)*rcb(2,nneighbor) +&
                                 rcb(3,nneighbor)*rcb(3,nneighbor)
                          IF (rcb2 < rcb2_max) THEN
                            neighbors(nneighbor) =&
                              atoms(jkind)%list(jatom)
                            nneighbor =  nneighbor + 1
                          END IF
                        END DO

                        IF (nneighbor > 1) THEN
                          DO ineighbor=1,nneighbor-1
                            CALL add_neighbor_node(&
                              neighbor_list=sca_neighbor_list,&
                              neighbor=neighbors(ineighbor),&
                              cell=cell_a,&
                              r=rcb(:,ineighbor))
                          END DO
                        END IF

                      END DO
                    END DO
                  END DO

                  DO icellb=cell_c(1)-ncell(1),cell_c(1)+ncell(1)
                    cell_b(1) = icellb
                    DO jcellb=cell_c(2)-ncell(2),cell_c(2)+ncell(2)
                      cell_b(2) = jcellb
                      DO kcellb=cell_c(3)-ncell(3),cell_c(3)+ncell(3)
                        cell_b(3) = kcellb

                        IF ((icellb == 0).AND.&
                            (jcellb == 0).AND.&
                            (kcellb == 0)) CYCLE

                        s(:) = REAL(cell_b(:),wp)
                        r(:) = scaled_to_real(s(:),cell) - rc(:)

                        s_max(:) = sc_max(:) - cell_b(:)
                        s_min(:) = sc_min(:) - cell_b(:)

                        IF (s_max(1) < grid_min(1,1,1,1)) CYCLE
                        IF (s_max(2) < grid_min(2,1,1,1)) CYCLE
                        IF (s_max(3) < grid_min(3,1,1,1)) CYCLE

                        igrid = ngrid(1)
                        jgrid = ngrid(2)
                        kgrid = ngrid(3)

                        IF (s_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (s_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (s_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO igrid=1,ngrid(1)
                          DO jgrid=1,ngrid(2)
                            DO kgrid=1,ngrid(3)

                              IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                              IF (s_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                              IF (s_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                              IF (s_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                              IF (s_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                              IF (s_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                              IF (s_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                              nneighbor = 1

                              DO iijk=1,nijk(igrid,jgrid,kgrid)
                                jatom = ijk(iijk,igrid,jgrid,kgrid)
                                rcb(:,nneighbor) =&
                                  pbc_coord(jkind)%r(:,jatom) + r(:)
                                  rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                                         rcb(2,nneighbor)*rcb(2,nneighbor) +&
                                         rcb(3,nneighbor)*rcb(3,nneighbor)
                                IF (rcb2 < rcb2_max) THEN
                                  neighbors(nneighbor) =&
                                    atoms(jkind)%list(jatom)
                                  nneighbor =  nneighbor + 1
                                END IF
                              END DO

                              IF (nneighbor > 1) THEN
                                DO ineighbor=1,nneighbor-1
                                  CALL add_neighbor_node(&
                                    neighbor_list=scb_neighbor_list,&
                                    neighbor=neighbors(ineighbor),&
                                    cell=cell_b,&
                                    r=rcb(:,ineighbor))
                                END DO
                              END IF

                            END DO
                          END DO
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (rcb,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rcb")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_gth_ppnl_neighbor_lists

! *****************************************************************************

  SUBROUTINE update_neighbor_lists(qs_env,globenv)

!   Purpose: Update all the required neighbor lists for Quickstep.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE update_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    INTEGER                           :: handle,nparticle

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    CALL timeset("update_neighbor_lists","I","",handle)

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sca_ppnl=sca_ppnl,&
                    scb_ppnl=scb_ppnl)

    nparticle = SIZE(particle_set)

    ALLOCATE (atom_of_kind(nparticle),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"atom_of_kind",nparticle*int_size)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             gth_potential_present=gth_potential_present,&
                             maxatom=maxatom)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Allocate work storage ***

    ALLOCATE (kind_a(maxatom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a",maxatom*int_size)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (natom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom",nkind*int_size)

    ALLOCATE (atoms(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atoms",nkind*int_size)

    ALLOCATE (pbc_coord(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord",nkind*int_size)

!   *** Calculate PBC coordinates ***

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atoms(ikind)%list)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           gth_potential=gth_potential,&
                           orb_basis_set=orb_basis_set,&
                           natom=natom(ikind),&
                           atom_list=atoms(ikind)%list)

      IF (ASSOCIATED(orb_basis_set)) THEN

        ALLOCATE (pbc_coord(ikind)%r(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r",&
                                         3*natom(ikind)*wp_size)
        ALLOCATE (pbc_coord(ikind)%s(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s",&
                                         3*natom(ikind)*wp_size)

        DO iatom=1,natom(ikind)
          atom_a = atoms(ikind)%list(iatom)
          ra_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          pbc_coord(ikind)%r(:,iatom) = ra_pbc(:)
          pbc_coord(ikind)%s(:,iatom) = real_to_scaled(ra_pbc(:),cell)
        END DO

      ELSE

        NULLIFY (pbc_coord(ikind)%r)
        NULLIFY (pbc_coord(ikind)%s)

      END IF

    END DO

!   *** Update the orbital-orbital overlap neighbor lists ***

    CALL update_sab_orb_neighbor_lists()

    IF (ionode.AND.globenv%print%sab_orb_neighbor_lists) THEN
      CALL write_neighbor_lists(sab_orb,"ORBITAL ORBITAL NEIGHBOR LISTS")
    END IF

!   *** Update the orbital GTH-PPL operator overlap list ***

    IF (gth_potential_present) THEN

!      IF (ANY(ppnl_present)) THEN

        CALL update_gth_ppl_neighbor_lists()

        IF (ionode.AND.globenv%print%sac_ppl_neighbor_lists) THEN
          CALL write_neighbor_lists(sac_ppl,"ORBITAL GTH-PPL NEIGHBOR LISTS")
        END IF

!      END IF

!      IF (ANY(ppnl_present)) THEN

        ALLOCATE (neighbors(maxatom),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbors",maxatom*int_size)

        CALL update_gth_ppnl_neighbor_lists()

        IF (globenv%print%sca_ppnl_neighbor_lists) THEN
          IF (npe == 1) THEN
            CALL write_neighbor_lists(sca_ppnl,&
                                      "GTH-PPNL ORBITAL(A) NEIGHBOR LISTS")
          ELSE
            DO ipe=0,npe-1
              IF (ipe == mype) THEN
                WRITE (string,"(A,I5)")&
                  "SCA_PPNL NEIGHBOR LISTS OF PROCESS",ipe
                CALL compress(string)
                CALL write_neighbor_lists(sca_ppnl,string)
              END IF
              CALL mp_sync(group)
            END DO
          END IF
        END IF

        IF (globenv%print%scb_ppnl_neighbor_lists) THEN
          IF (npe == 1) THEN
            CALL write_neighbor_lists(scb_ppnl,&
                                      "GTH-PPNL ORBITAL(B) NEIGHBOR LISTS")
          ELSE
            DO ipe=0,npe-1
              IF (ipe == mype) THEN
                WRITE (string,"(A,I5)")&
                  "SCB_PPNL NEIGHBOR LISTS OF PROCESS",ipe
                CALL compress(string)
                CALL write_neighbor_lists(scb_ppnl,string)
              END IF
              CALL mp_sync(group)
            END DO
          END IF
        END IF

!      END IF

    END IF

!   *** Release work storage ***

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    DEALLOCATE (kind_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a")

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      NULLIFY (atoms(ikind)%list)
      IF (ASSOCIATED(pbc_coord(ikind)%r)) THEN
        DEALLOCATE (pbc_coord(ikind)%r,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r")
      END IF
      IF (ASSOCIATED(pbc_coord(ikind)%s)) THEN
        DEALLOCATE (pbc_coord(ikind)%s,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s")
      END IF
    END DO

    DEALLOCATE (atoms,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atoms")

    DEALLOCATE (pbc_coord,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord")

    IF (gth_potential_present) THEN

!      DEALLOCATE (ppnl_present,STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routine,"ppnl_present")

      IF (ALLOCATED(neighbors)) THEN
        DEALLOCATE (neighbors,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbors")
      END IF

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE update_neighbor_lists

! *****************************************************************************

  SUBROUTINE update_sab_orb_neighbor_lists()

!   Purpose: Update the neighbor lists for the calculation of the orbital-
!            orbital overlap matrix.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE update_sab_orb_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    INTEGER :: ab,handle,natom_a
    LOGICAL :: equal_kinds

!   ---------------------------------------------------------------------------

    CALL timeset("update_sab_orb_neighbor_lists","I","",handle)

    IF (.NOT.ASSOCIATED(sab_orb)) THEN
      CALL stop_program(routine,&
                        "sab_orb is not associated and cannot be updated")
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind
      DO jkind=ikind,nkind

        ab = ikind + jkind*(jkind - 1)/2

        neighbor_list_set => sab_orb(ab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        equal_kinds = (ikind == jkind)

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   r_max=rab_max)
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/sab_max(:)))

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          neighbor_list => first_list(neighbor_list_set)

          DO WHILE (ASSOCIATED(neighbor_list))
            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a)
            CALL init_neighbor_list(neighbor_list)
            iatom = atom_of_kind(atom_a)
            kind_a(iatom)%neighbor_list => neighbor_list
            neighbor_list => next(neighbor_list)
          END DO

          DO jatom=1,natom(jkind)

            atom_b = atoms(jkind)%list(jatom)
            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),wp)
                  rb(:) = scaled_to_real(sb(:),cell)

                  IF (equal_kinds) THEN
                    natom_a = jatom
                  ELSE
                    natom_a = natom(ikind)
                  END IF

                  DO iatom=1,natom_a
                    rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                    IF (rab2 < rab2_max) THEN
                      CALL add_neighbor_node(&
                        neighbor_list=kind_a(iatom)%neighbor_list,&
                        neighbor=atom_b,&
                        cell=cell_b,&
                        r=rab(:))
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(ikind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",&
                             natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
              b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          neighbor_list => first_list(neighbor_list_set)

          DO WHILE (ASSOCIATED(neighbor_list))
            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a)
            CALL init_neighbor_list(neighbor_list)
            iatom = atom_of_kind(atom_a)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = iatom
            kind_a(iatom)%neighbor_list => neighbor_list
            neighbor_list => next(neighbor_list)
          END DO

          DO jatom=1,natom(jkind)

            atom_b = atoms(jkind)%list(jatom)
            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),wp)
                  rb(:) = scaled_to_real(sb(:),cell)
                  sb_min(:) = sb(:) - sab_max(:)
                  sb_max(:) = sb(:) + sab_max(:)

                  IF (sb_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sb_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sb_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sb_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sb_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sb_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sb_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sb_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          iatom = ijk(iijk,igrid,jgrid,kgrid)
                          IF (equal_kinds) THEN
                            IF (jatom < iatom) CYCLE
                          END IF
                          rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                          IF (rab2 < rab2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_a(iatom)%neighbor_list,&
                              neighbor=atom_b,&
                              cell=cell_b,&
                              r=rab(:))
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          CALL clean_neighbor_list_set(neighbor_list_set)

!         *** Release work storage ***

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

      END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE update_sab_orb_neighbor_lists

! *****************************************************************************

  SUBROUTINE update_gth_ppl_neighbor_lists()

!   Purpose: Update the neighbor lists for the calculation of the local part of
!            the GTH pseudo potential (PPL).

!   History: - Creation (17.09.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE update_gth_ppl_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    INTEGER :: ac,handle

!   ---------------------------------------------------------------------------

    CALL timeset("update_gth_ppl_neighbor_lists","I","",handle)

    IF (.NOT.ASSOCIATED(sac_ppl)) THEN
      CALL stop_program(routine,&
                        "sac_ppl is not associated and cannot be updated")
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind
      DO kkind=1,nkind

        ac = ikind + (kkind - 1)*nkind

        neighbor_list_set => sac_ppl(ac)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   r_max=rac_max)
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/sac_max(:)))

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          neighbor_list => first_list(neighbor_list_set)

          DO WHILE (ASSOCIATED(neighbor_list))
            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a)
            CALL init_neighbor_list(neighbor_list)
            iatom = atom_of_kind(atom_a)
            kind_a(iatom)%neighbor_list => neighbor_list
            neighbor_list => next(neighbor_list)
          END DO

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)
            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)

                  DO iatom=1,natom(ikind)
                    rac(:) = rc(:) - pbc_coord(ikind)%r(:,iatom)
                    rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                    IF (rac2 < rac2_max) THEN
                      CALL add_neighbor_node(&
                        neighbor_list=kind_a(iatom)%neighbor_list,&
                        neighbor=atom_c,&
                        cell=cell_c,&
                        r=rac(:))
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(ikind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
              b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          neighbor_list => first_list(neighbor_list_set)

          DO WHILE (ASSOCIATED(neighbor_list))
            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a)
            CALL init_neighbor_list(neighbor_list)
            iatom = atom_of_kind(atom_a)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = iatom
            kind_a(iatom)%neighbor_list => neighbor_list
            neighbor_list => next(neighbor_list)
          END DO

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)
            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)
                  sc_min(:) = sc(:) - sac_max(:)
                  sc_max(:) = sc(:) + sac_max(:)

                  IF (sc_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sc_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sc_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sc_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          iatom = ijk(iijk,igrid,jgrid,kgrid)
                          rac(:) = rc(:) - pbc_coord(ikind)%r(:,iatom)
                          rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                          IF (rac2 < rac2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_a(iatom)%neighbor_list,&
                              neighbor=atom_c,&
                              cell=cell_c,&
                              r=rac(:))
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          CALL clean_neighbor_list_set(neighbor_list_set)

!         *** Release work storage ***

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

      END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE update_gth_ppl_neighbor_lists

! *****************************************************************************

  SUBROUTINE update_gth_ppnl_neighbor_lists()

!   Purpose: Update the neighbor lists for the calculation of the non-local
!            part of the GTH pseudo potential (PPNL). Atom c is a PPNL operator
!            atom, atom a is an atom in the cell(0,0,0), and atom b is an atom
!            is an atom in any cell(i,j,k).

!   History: - Creation (17.09.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE update_gth_ppnl_neighbor_lists (MODULE qs_neighbor_lists)"

!   *** Local variables ***

    INTEGER :: cb,handle,icellb,jcellb,kcellb

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rcb

!   ---------------------------------------------------------------------------

    CALL timeset("update_gth_ppnl_neighbor_lists","I","",handle)

    IF (.NOT.ASSOCIATED(sca_ppnl)) THEN
      CALL stop_program(routine,&
                        "sca_ppnl is not associated and cannot be updated")
    END IF

    IF (.NOT.ASSOCIATED(scb_ppnl)) THEN
      CALL stop_program(routine,&
                        "scb_ppnl is not associated and cannot be updated")
    END IF

!   *** Allocate work storage ***

    ALLOCATE (rcb(3,maxatom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rcb",3*maxatom*wp_size)

!   *** Loop over all atomic kind pairs ***

    DO kkind=1,nkind
      DO jkind=1,nkind

        cb = kkind + (jkind - 1)*nkind

        neighbor_list_set => sca_ppnl(cb)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        IF (jkind == 1) THEN
          sca_neighbor_list_set => sca_ppnl(kkind)%neighbor_list_set
          scb_neighbor_list_set => scb_ppnl(kkind)%neighbor_list_set
          CALL init_neighbor_list_set(sca_neighbor_list_set)
          CALL init_neighbor_list_set(scb_neighbor_list_set)
        END IF

!       *** Calculate the square of the maximum interaction distance ***

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   r_max=rcb_max)
        rcb2_max = rcb_max*rcb_max

        r(:) = rcb_max
        scb_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(scb_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/scb_max(:)))

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)

!           *** Distribute the GTH-PPNL operator atoms ***

            IF (mype /= MODULO(atom_c,npe)) CYCLE

            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)

                  sca_neighbor_list =>&
                    find_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c)

                  IF (.NOT.ASSOCIATED(sca_neighbor_list)) THEN
                    CALL add_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=sca_neighbor_list)
                    CALL add_neighbor_list(&
                      neighbor_list_set=scb_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=scb_neighbor_list)
                  ELSE
                    scb_neighbor_list =>&
                      find_neighbor_list(&
                        neighbor_list_set=scb_neighbor_list_set,&
                        atom=atom_c,&
                        cell=cell_c)
                  END IF

                  nneighbor = 1

                  DO jatom=1,natom(jkind)
                    rcb(:,nneighbor) = pbc_coord(jkind)%r(:,jatom) - rc(:)
                    rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                           rcb(2,nneighbor)*rcb(2,nneighbor) +&
                           rcb(3,nneighbor)*rcb(3,nneighbor)
                    IF (rcb2 < rcb2_max) THEN
                      neighbors(nneighbor) = atoms(jkind)%list(jatom)
                      nneighbor =  nneighbor + 1
                    END IF
                  END DO

                  IF (nneighbor > 1) THEN
                    DO ineighbor=1,nneighbor-1
                      CALL add_neighbor_node(&
                        neighbor_list=sca_neighbor_list,&
                        neighbor=neighbors(ineighbor),&
                        cell=cell_a,&
                        r=rcb(:,ineighbor))
                    END DO
                  END IF

                  DO icellb=cell_c(1)-ncell(1),cell_c(1)+ncell(1)
                    cell_b(1) = icellb
                    DO jcellb=cell_c(2)-ncell(2),cell_c(2)+ncell(2)
                      cell_b(2) = jcellb
                      DO kcellb=cell_c(3)-ncell(3),cell_c(3)+ncell(3)
                        cell_b(3) = kcellb

                        IF ((icellb == 0).AND.&
                            (jcellb == 0).AND.&
                            (kcellb == 0)) CYCLE

                        s(:) = REAL(cell_b(:),wp)
                        r(:) = scaled_to_real(s(:),cell) - rc(:)

                        nneighbor = 1

                        DO jatom=1,natom(jkind)
                          rcb(:,nneighbor) =&
                            pbc_coord(jkind)%r(:,jatom) + r(:)
                          rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                                 rcb(2,nneighbor)*rcb(2,nneighbor) +&
                                 rcb(3,nneighbor)*rcb(3,nneighbor)
                          IF (rcb2 < rcb2_max) THEN
                            neighbors(nneighbor) = atoms(jkind)%list(jatom)
                            nneighbor =  nneighbor + 1
                          END IF
                        END DO

                        IF (nneighbor > 1) THEN
                          DO ineighbor=1,nneighbor-1
                            CALL add_neighbor_node(&
                              neighbor_list=scb_neighbor_list,&
                              neighbor=neighbors(ineighbor),&
                              cell=cell_b,&
                              r=rcb(:,ineighbor))
                          END DO
                        END IF

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(jkind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",natom(jkind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
              b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO jatom=1,natom(jkind)
            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)
            igrid = MAX(1,CEILING((sb_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sb_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sb_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = jatom
          END DO

          DO katom=1,natom(kkind)

            atom_c = atoms(kkind)%list(katom)

!           *** Distribute the GTH-PPNL operator atoms ***

            IF (mype /= MODULO(atom_c,npe)) CYCLE

            sc_pbc(:) = pbc_coord(kkind)%s(:,katom)

            DO icell=-ncell(1),ncell(1)
              cell_c(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_c(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_c(3) = kcell

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),wp)
                  rc(:) = scaled_to_real(sc(:),cell)

                  sc_min(:) = sc(:) - scb_max(:)
                  sc_max(:) = sc(:) + scb_max(:)

                  IF (sc_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sc_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sc_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  sca_neighbor_list =>&
                    find_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c)

                  IF (.NOT.ASSOCIATED(sca_neighbor_list)) THEN
                    CALL add_neighbor_list(&
                      neighbor_list_set=sca_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=sca_neighbor_list)
                    CALL add_neighbor_list(&
                      neighbor_list_set=scb_neighbor_list_set,&
                      atom=atom_c,&
                      cell=cell_c,&
                      neighbor_list=scb_neighbor_list)
                  ELSE
                    scb_neighbor_list =>&
                      find_neighbor_list(&
                        neighbor_list_set=scb_neighbor_list_set,&
                        atom=atom_c,&
                        cell=cell_c)
                  END IF

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sc_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sc_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sc_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        nneighbor = 1

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          jatom = ijk(iijk,igrid,jgrid,kgrid)
                          rcb(:,nneighbor) =&
                            pbc_coord(jkind)%r(:,jatom) - rc(:)
                          rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                                 rcb(2,nneighbor)*rcb(2,nneighbor) +&
                                 rcb(3,nneighbor)*rcb(3,nneighbor)
                          IF (rcb2 < rcb2_max) THEN
                            neighbors(nneighbor) =&
                              atoms(jkind)%list(jatom)
                            nneighbor =  nneighbor + 1
                          END IF
                        END DO

                        IF (nneighbor > 1) THEN
                          DO ineighbor=1,nneighbor-1
                            CALL add_neighbor_node(&
                              neighbor_list=sca_neighbor_list,&
                              neighbor=neighbors(ineighbor),&
                              r=rcb(:,ineighbor),&
                              cell=cell_a)
                          END DO
                        END IF

                      END DO
                    END DO
                  END DO

                  DO icellb=cell_c(1)-ncell(1),cell_c(1)+ncell(1)
                    cell_b(1) = icellb
                    DO jcellb=cell_c(2)-ncell(2),cell_c(2)+ncell(2)
                      cell_b(2) = jcellb
                      DO kcellb=cell_c(3)-ncell(3),cell_c(3)+ncell(3)
                        cell_b(3) = kcellb

                        IF ((icellb == 0).AND.&
                            (jcellb == 0).AND.&
                            (kcellb == 0)) CYCLE

                        s(:) = REAL(cell_b(:),wp)
                        r(:) = scaled_to_real(s(:),cell) - rc(:)

                        s_max(:) = sc_max(:) - cell_b(:)
                        s_min(:) = sc_min(:) - cell_b(:)

                        IF (s_max(1) < grid_min(1,1,1,1)) CYCLE
                        IF (s_max(2) < grid_min(2,1,1,1)) CYCLE
                        IF (s_max(3) < grid_min(3,1,1,1)) CYCLE

                        igrid = ngrid(1)
                        jgrid = ngrid(2)
                        kgrid = ngrid(3)

                        IF (s_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (s_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (s_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO igrid=1,ngrid(1)
                          DO jgrid=1,ngrid(2)
                            DO kgrid=1,ngrid(3)

                              IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                              IF (s_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                              IF (s_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                              IF (s_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                              IF (s_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                              IF (s_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                              IF (s_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                              nneighbor = 1

                              DO iijk=1,nijk(igrid,jgrid,kgrid)
                                jatom = ijk(iijk,igrid,jgrid,kgrid)
                                rcb(:,nneighbor) =&
                                  pbc_coord(jkind)%r(:,jatom) + r(:)
                                  rcb2 = rcb(1,nneighbor)*rcb(1,nneighbor) +&
                                         rcb(2,nneighbor)*rcb(2,nneighbor) +&
                                         rcb(3,nneighbor)*rcb(3,nneighbor)
                                IF (rcb2 < rcb2_max) THEN
                                  neighbors(nneighbor) =&
                                    atoms(jkind)%list(jatom)
                                  nneighbor =  nneighbor + 1
                                END IF
                              END DO

                              IF (nneighbor > 1) THEN
                                DO ineighbor=1,nneighbor-1
                                  CALL add_neighbor_node(&
                                    neighbor_list=scb_neighbor_list,&
                                    neighbor=neighbors(ineighbor),&
                                    r=rcb(:,ineighbor),&
                                    cell=cell_b)
                                END DO
                              END IF

                            END DO
                          END DO
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

        IF (jkind == nkind) CALL clean_neighbor_list_set(neighbor_list_set)

      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (rcb,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rcb")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE update_gth_ppnl_neighbor_lists

! *****************************************************************************

  SUBROUTINE write_neighbor_lists(neighbor_lists,name)

!   Purpose: Write a set of neighbor lists to the output unit.

!   History: - Creation (04.03.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: neighbor_lists
    CHARACTER(LEN=*), INTENT(IN)                          :: name

!   *** Local variables ***

    INTEGER :: ab,atom_a,atom_b,i,nneighbor

    REAL(wp), DIMENSION(3) :: ra,rab,rb
    INTEGER, DIMENSION(3)  :: cell_a,cell_b

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

        ra(:) = pbc(particle_set(atom_a)%r,cell,cell_a)

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

END MODULE qs_neighbor_lists
