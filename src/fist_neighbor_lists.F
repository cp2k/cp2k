!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_neighbor_lists [1.0] *
!!
!!   NAME
!!     fist_neighbor_lists
!!
!!   FUNCTION
!!     Build all neighbor lists for FIST.
!!
!!   AUTHOR
!!     Matthias Krack (19.11.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fist_neighbor_lists

  USE kinds, ONLY: int_size,&
                   dp,&
                   dp_size

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: LOW,global_environment_type
  USE molecule_types,                  ONLY: particle_node_type,&
                                             linklist_exclusion
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_list_types,          ONLY: add_neighbor_list,&
                                             add_neighbor_node,&
                                             allocate_neighbor_list_set,&
                                             clean_neighbor_list_set,&
                                             deallocate_neighbor_list_set,&
                                             find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_node,&
                                             init_neighbor_list,&
                                             init_neighbor_list_set,&
                                             neighbor_list_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next,&
                                             set_neighbor_node
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global types of the module ***

  TYPE atoms_type
    INTEGER, DIMENSION(:), POINTER :: list
  END TYPE atoms_type

  TYPE excl_type
    INTEGER, DIMENSION(:,:), POINTER :: list
  END TYPE excl_type

  TYPE pbc_coord_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: r,s
  END TYPE pbc_coord_type

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "fist_neighbor_lists"

! *** Global variables of the module ***

  TYPE(linklist_exclusion), POINTER     :: excl_node
  TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
  TYPE(neighbor_list_type), POINTER     :: neighbor_list
  TYPE(neighbor_node_type), POINTER     :: neighbor_node
  CHARACTER(LEN=40)                     :: string
  CHARACTER(LEN=8)                      :: unit_of_length_name
  REAL(KIND = dp)                              :: a_max,a_min,&
                                           b_max,b_min,&
                                           c_max,c_min,&
                                           rab_max,rab2,rab2_max,&
                                           subcells,unit_of_length
  INTEGER                               :: atom_a,atom_b,group,&
                                           iatom,icell,iexcl,igrid,iijk,ikind,&
                                           ineighbor,inode,ipe,istat,&
                                           jatom,jcell,jgrid,jkind,&
                                           kcell,kgrid,maxatom,maxexcl,mype,&
                                           nkind,nneighbor,nnode,npe,&
                                           output_unit
  LOGICAL                               :: ionode,print_cell_parameters

  REAL(KIND = dp), DIMENSION(3) :: r,ra_pbc,rab,rb,s,s_max,s_min,sab,sab_max,sb,&
                            sb_max,sb_min,sb_pbc,sa_pbc
  INTEGER, DIMENSION(3)  :: cell_a,cell_b,ncell,ngrid,periodic

  TYPE(atoms_type), DIMENSION(:), ALLOCATABLE           :: atoms
  TYPE(excl_type), DIMENSION(:), ALLOCATABLE            :: excl_a
  TYPE(neighbor_list_p_type), DIMENSION(:), ALLOCATABLE :: kind_a
  TYPE(pbc_coord_type), DIMENSION(:), ALLOCATABLE       :: pbc_coord
  REAL(KIND = dp), DIMENSION(:), ALLOCATABLE                   :: kind_radius
  INTEGER, DIMENSION(:), ALLOCATABLE                    :: atom_of_kind,&
                                                           kind_of,&
                                                           natom,neighbors
  REAL(KIND = dp), DIMENSION(:,:,:,:), ALLOCATABLE             :: grid_max,grid_min
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE                :: nijk
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE              :: ijk

! *** Public subroutines ***

  PUBLIC :: build_fist_neighbor_lists,&
            rebuild_fist_neighbor_lists

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_fist_neighbor_lists(atomic_kind_set,particle_set,pnode,&
                                       cell,r_cut,nonbonded,globenv)

!   Purpose: Build all the required neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                              :: cell
    TYPE(global_environment_type), INTENT(IN)             :: globenv
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nonbonded
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(particle_node_type), DIMENSION(:), POINTER       :: pnode
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)                  :: r_cut

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_fist_neighbor_lists"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    INTEGER                           :: handle,nparticle

!   ---------------------------------------------------------------------------

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

    nparticle = SIZE(particle_set)

    ALLOCATE (atom_of_kind(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",nparticle*int_size)
    ALLOCATE (kind_of(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of",nparticle*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of,&
                             maxatom=maxatom)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!MK *** pnode and neighbor list structure do not fit very well ***

    nnode = SIZE(pnode)
    maxexcl = 0
    DO inode=1,nnode
      maxexcl = MAX(maxexcl,pnode(inode)%nexcl)
    END DO

!   *** Allocate work storage ***

    ALLOCATE (kind_a(maxatom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a",maxatom*int_size)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (excl_a(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"excl_a",nkind*int_size)

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
                           natom=natom(ikind),&
                           atom_list=atoms(ikind)%list)

      ALLOCATE (excl_a(ikind)%list(maxexcl,natom(ikind)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"excl_a(ikind)%list",&
                                       maxexcl*natom(ikind)*dp_size)
      excl_a(ikind)%list(:,:) = 0

      ALLOCATE (pbc_coord(ikind)%r(3,natom(ikind)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r",&
                                       3*natom(ikind)*dp_size)
      ALLOCATE (pbc_coord(ikind)%s(3,natom(ikind)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s",&
                                       3*natom(ikind)*dp_size)

      DO iatom=1,natom(ikind)
        atom_a = atoms(ikind)%list(iatom)
        ra_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
        pbc_coord(ikind)%r(:,iatom) = ra_pbc(:)
        pbc_coord(ikind)%s(:,iatom) = real_to_scaled(ra_pbc(:),cell)
      END DO

    END DO

    DO inode=1,nnode
      atom_a = pnode(inode)%p%iatom
      ikind = kind_of(atom_a)
      iatom = atom_of_kind(atom_a)
      excl_node => pnode(inode)%ex
      DO iexcl=1,pnode(inode)%nexcl
        excl_a(ikind)%list(iexcl,iatom) = excl_node%p%iatom
        excl_node => excl_node%next
      END DO
    END DO

!   *** Build the nonbonded neighbor lists ***

    CALL build_nonbonded(atomic_kind_set,particle_set,cell,r_cut,nonbonded)

   IF (ionode.AND.globenv%print%sab_orb_neighbor_lists) THEN
!    IF (.TRUE.) THEN
      CALL write_neighbor_lists(nonbonded,"NONBONDED NEIGHBOR LISTS",&
                                particle_set,cell)
    END IF

!   *** Release work storage ***

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of")

    DEALLOCATE (kind_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a")

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      NULLIFY (atoms(ikind)%list)
      IF (ASSOCIATED(excl_a(ikind)%list)) THEN
        DEALLOCATE (excl_a(ikind)%list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"excl_a(ikind)%list")
      END IF
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

    DEALLOCATE (excl_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"excl_a")

    DEALLOCATE (pbc_coord,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord")

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_fist_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_nonbonded(atomic_kind_set,particle_set,cell,r_cut,nonbonded)

!   Purpose: Build the nonbonded neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                              :: cell
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nonbonded
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)                  :: r_cut

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_nonbonded"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: ab,handle,natom_a
    LOGICAL :: equal_kinds

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(nonbonded)) THEN
      DO ab=1,SIZE(nonbonded)
        CALL deallocate_neighbor_list_set(nonbonded(ab)%neighbor_list_set)
      END DO
      DEALLOCATE (nonbonded,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nonbonded")
    END IF

    ALLOCATE (nonbonded(nkind*(nkind + 1)/2),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nonbonded",0)

    DO ab=1,SIZE(nonbonded)
      NULLIFY (nonbonded(ab)%neighbor_list_set)
    END DO

!   *** Print headline ***

    IF (print_cell_parameters) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE NONBONDED NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind
      DO jkind=ikind,nkind

        ab = ikind + jkind*(jkind - 1)/2

        equal_kinds = (ikind == jkind)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = r_cut(ikind,jkind)
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_dp*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,ngrid,&
            scaled_to_real(1.0_dp/REAL(ngrid(:),dp),cell)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          nonbonded(ab)%neighbor_list_set,&
                                        r_max=rab_max)
        neighbor_list_set => nonbonded(ab)%neighbor_list_set

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

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),dp)
                  rb(:) = scaled_to_real(sb(:),cell)

                  IF (equal_kinds) THEN
                    natom_a = jatom
                  ELSE
                    natom_a = natom(ikind)
                  END IF

                  DO iatom=1,natom_a
                    rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                    sab ( : ) = real_to_scaled ( rab ( : ), cell )
                    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                    IF (rab2 < rab2_max) THEN
                      IF (rab2 > 1.0E-6_dp) THEN
                        CALL add_neighbor_node(&
                          neighbor_list=kind_a(iatom)%neighbor_list,&
                          neighbor=atom_b,&
                          cell=cell_b,&
                          r=sab(:),&
                          exclusion_list=excl_a(ikind)%list(:,iatom))
                      END IF
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*dp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*dp_size)
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
            a_min = REAL(igrid-1,dp)/REAL(ngrid(1),dp) - 0.5_dp
            a_max = REAL(igrid,dp)/REAL(ngrid(1),dp) - 0.5_dp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,dp)/REAL(ngrid(2),dp) - 0.5_dp
              b_max = REAL(jgrid,dp)/REAL(ngrid(2),dp) - 0.5_dp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,dp)/REAL(ngrid(3),dp) - 0.5_dp
                c_max = REAL(kgrid,dp)/REAL(ngrid(3),dp) - 0.5_dp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_dp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_dp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_dp)*ngrid(3)))
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

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),dp)
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
                          sab ( : ) = real_to_scaled ( rab ( : ), cell )
                          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                          IF (rab2 < rab2_max) THEN
                            IF (rab2 > 1.0E-6_dp) THEN
                              CALL add_neighbor_node(&
                                neighbor_list=kind_a(iatom)%neighbor_list,&
                                neighbor=atom_b,&
                                cell=cell_b,&
                                r=sab(:),&
                                exclusion_list=excl_a(ikind)%list(:,iatom))
                            END IF
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

  END SUBROUTINE build_nonbonded

! *****************************************************************************

  SUBROUTINE rebuild_fist_neighbor_lists(atomic_kind_set,particle_set,pnode,&
                                         cell,r_cut,nonbonded,globenv)

!   Purpose: Rebuild all the required neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                              :: cell
    TYPE(global_environment_type), INTENT(IN)             :: globenv
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nonbonded
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(particle_node_type), DIMENSION(:), POINTER       :: pnode
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)                  :: r_cut

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "rebuild_fist_neighbor_lists"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind

    INTEGER :: handle,nparticle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    nparticle = SIZE(particle_set)

    ALLOCATE (atom_of_kind(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",nparticle*int_size)
    ALLOCATE (kind_of(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of",nparticle*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of,&
                             maxatom=maxatom)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!MK *** pnode and neighbor list structure do not fit very well ***

    nnode = SIZE(pnode)
    maxexcl = 0
    DO inode=1,nnode
      maxexcl = MAX(maxexcl,pnode(inode)%nexcl)
    END DO

!   *** Allocate work storage ***

    ALLOCATE (kind_a(maxatom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a",maxatom*int_size)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (excl_a(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"excl_a",nkind*int_size)

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
                           natom=natom(ikind),&
                           atom_list=atoms(ikind)%list)

      ALLOCATE (excl_a(ikind)%list(maxexcl,natom(ikind)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"excl_a(ikind)%list",&
                                       maxexcl*natom(ikind)*dp_size)
      excl_a(ikind)%list(:,:) = 0

      ALLOCATE (pbc_coord(ikind)%r(3,natom(ikind)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r",&
                                       3*natom(ikind)*dp_size)
      ALLOCATE (pbc_coord(ikind)%s(3,natom(ikind)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s",&
                                       3*natom(ikind)*dp_size)

      DO iatom=1,natom(ikind)
        atom_a = atoms(ikind)%list(iatom)
        ra_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
        pbc_coord(ikind)%r(:,iatom) = ra_pbc(:)
        pbc_coord(ikind)%s(:,iatom) = real_to_scaled(ra_pbc(:),cell)
      END DO

    END DO

    DO inode=1,nnode
      atom_a = pnode(inode)%p%iatom
      ikind = kind_of(atom_a)
      iatom = atom_of_kind(atom_a)
      excl_node => pnode(inode)%ex
      DO iexcl=1,pnode(inode)%nexcl
        excl_a(ikind)%list(iexcl,iatom) = excl_node%p%iatom
        excl_node => excl_node%next
      END DO
    END DO

!   *** Rebuild the nonbonded neighbor lists ***

    CALL rebuild_nonbonded(atomic_kind_set,particle_set,cell,r_cut,nonbonded)

    IF (ionode.AND.globenv%print%sab_orb_neighbor_lists) THEN
      CALL write_neighbor_lists(nonbonded,"NONBONDED NEIGHBOR LISTS",&
                                particle_set,cell)
    END IF

!   *** Release work storage ***

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of")

    DEALLOCATE (kind_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_a")

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      NULLIFY (atoms(ikind)%list)
      IF (ASSOCIATED(excl_a(ikind)%list)) THEN
        DEALLOCATE (excl_a(ikind)%list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"excl_a(ikind)%list")
      END IF
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

    DEALLOCATE (excl_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"excl_a")

    DEALLOCATE (pbc_coord,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord")

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE rebuild_fist_neighbor_lists

! *****************************************************************************

  SUBROUTINE rebuild_nonbonded(atomic_kind_set,particle_set,cell,r_cut,&
                               nonbonded)

!   Purpose: Rebuild the nonbonded neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                              :: cell
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nonbonded
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)                  :: r_cut

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "rebuild_nonbonded"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: ab,handle,natom_a
    LOGICAL :: equal_kinds

!   ---------------------------------------------------------------------------

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind
      DO jkind=ikind,nkind

        ab = ikind + jkind*(jkind - 1)/2

        neighbor_list_set => nonbonded(ab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        equal_kinds = (ikind == jkind)

        rab_max = r_cut(ikind,jkind)
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_dp*subcells/sab_max(:)))

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

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),dp)
                  rb(:) = scaled_to_real(sb(:),cell)

                  IF (equal_kinds) THEN
                    natom_a = jatom
                  ELSE
                    natom_a = natom(ikind)
                  END IF

                  DO iatom=1,natom_a
                    rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                    sab ( : ) = real_to_scaled ( rab ( : ), cell )
                    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                    IF (rab2 < rab2_max) THEN
                      IF (rab2 > 1.0E-6_dp) THEN
                        CALL add_neighbor_node(&
                          neighbor_list=kind_a(iatom)%neighbor_list,&
                          neighbor=atom_b,&
                          cell=cell_b,&
                          r=sab(:),&
                          exclusion_list=excl_a(ikind)%list(:,iatom))
                      END IF
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*dp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*dp_size)
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
            a_min = REAL(igrid-1,dp)/REAL(ngrid(1),dp) - 0.5_dp
            a_max = REAL(igrid,dp)/REAL(ngrid(1),dp) - 0.5_dp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,dp)/REAL(ngrid(2),dp) - 0.5_dp
              b_max = REAL(jgrid,dp)/REAL(ngrid(2),dp) - 0.5_dp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,dp)/REAL(ngrid(3),dp) - 0.5_dp
                c_max = REAL(kgrid,dp)/REAL(ngrid(3),dp) - 0.5_dp
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
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_dp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_dp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_dp)*ngrid(3)))
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

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),dp)
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
                          sab ( : ) = real_to_scaled ( rab ( : ), cell )
                          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                          IF (rab2 < rab2_max) THEN
                            IF (rab2 > 1.0E-6_dp) THEN
                              CALL add_neighbor_node(&
                                neighbor_list=kind_a(iatom)%neighbor_list,&
                                neighbor=atom_b,&
                                cell=cell_b,&
                                r=sab(:),&
                                exclusion_list=excl_a(ikind)%list(:,iatom))
                            END IF
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

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

!MK     CALL clean_neighbor_list_set(neighbor_list_set)

      END DO
    END DO

  END SUBROUTINE rebuild_nonbonded

! *****************************************************************************

  SUBROUTINE write_neighbor_lists(neighbor_lists,name,particle_set,cell)

!   Purpose: Write a set of neighbor lists to the output unit.

!   History: - Creation (04.03.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                              :: cell
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: neighbor_lists
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    CHARACTER(LEN=*), INTENT(IN)                          :: name

!   *** Local variables ***

    INTEGER :: ab,atom_a,atom_b,i,nneighbor

    REAL(KIND = dp), DIMENSION(3) :: ra,rab,rb
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

END MODULE fist_neighbor_lists
