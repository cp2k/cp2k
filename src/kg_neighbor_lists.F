!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_neighbor_lists [1.0] *
!!
!!   NAME
!!     kg_neighbor_lists
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists for KG 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_neighbor_lists

  USE kinds, ONLY: int_size,&
                   dbl,&      
                   dp_size

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE basis_set_types,          ONLY: get_gto_basis_set,&
                                      gto_basis_set_type
  USE checkpoint_handler,       ONLY: write_checkpoint_information
  USE dft_types,                ONLY: dft_control_type
  USE external_potential_types, ONLY: get_potential,&
                                      gth_potential_type
  USE global_types,             ONLY: LOW,global_environment_type
  USE message_passing,          ONLY: mp_sync
  USE particle_types,           ONLY: particle_type
  USE kg_environment_types,     ONLY: get_kg_env,&
                                      kg_environment_type,&
                                      set_kg_env
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
                                      next,&
                                      set_neighbor_node
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
    REAL(dbl), DIMENSION(:,:), POINTER :: r,s
  END TYPE pbc_coord_type

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_neighbor_lists"

! *** Global variables of the module ***

  TYPE(cell_type), POINTER              :: cell
  TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
  TYPE(neighbor_list_type), POINTER     :: neighbor_list
  TYPE(neighbor_node_type), POINTER     :: neighbor_node
  CHARACTER(LEN=40)                     :: string
  CHARACTER(LEN=8)                      :: unit_of_length_name
  REAL(dbl)                              :: a_max,a_min,&
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

  REAL(dbl), DIMENSION(3) :: r,ra_pbc,rab,rac,rb,rc,rc_pbc,s,s_max,s_min,&
                            sac_max,sb,sb_max,sb_min,sb_pbc,sc,sc_max,&
                            sc_min,sa_pbc,sc_pbc,sca_max,scb_max
  INTEGER, DIMENSION(3)  :: cell_a,cell_b,cell_c,ncell,ngrid,periodic

  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sac_ppl
  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sac_ppl_aux
  TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

  TYPE(atoms_type), DIMENSION(:), ALLOCATABLE           :: atoms
  TYPE(neighbor_list_p_type), DIMENSION(:), ALLOCATABLE :: kind_a,kind_b
  TYPE(pbc_coord_type), DIMENSION(:), ALLOCATABLE       :: pbc_coord
  REAL(dbl), DIMENSION(:), ALLOCATABLE                   :: aux_radius,&
                                                           orb_radius,&
                                                           ppl_radius
  INTEGER, DIMENSION(:), ALLOCATABLE                    :: atom_of_kind,&
                                                           natom,neighbors
  LOGICAL, DIMENSION(:), ALLOCATABLE                    :: aux_present,&
                                                           orb_present,&
                                                           ppl_present
  REAL(dbl), DIMENSION(:,:,:,:), ALLOCATABLE             :: grid_max,grid_min
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE                :: nijk
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE              :: ijk
  TYPE(dft_control_type), POINTER :: dft_control

! *** Public subroutines ***

  PUBLIC :: build_neighbor_lists,&
            rebuild_neighbor_lists

!!***
! *****************************************************************************

CONTAINS

!*****************************************************************************

  SUBROUTINE build_neighbor_lists(kg_env,globenv)

!   Purpose: Build all the required neighbor lists for Kim-Gordon MD.

!   History: - Creation (28.06.2003,gt)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(kg_environment_type), INTENT(INOUT)  :: kg_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_neighbor_lists"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set
    REAL(dbl)                          :: r_skin
    INTEGER                           :: handle

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

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

!   *** Set Verlet skin ***

    r_skin = 0.0_dbl

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set,&
                    sac_ppl=sac_ppl,& 
                    sac_ppl_aux=sac_ppl_aux)

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
    IF (istat /= 0) CALL stop_memory(routine,"orb_radius",nkind*dp_size)
    orb_radius = 0.0_dbl

    ALLOCATE (aux_present(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"aux_present",nkind*int_size)

    ALLOCATE (aux_radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"aux_radius",nkind*dp_size)
    aux_radius = 0.0_dbl

    IF (gth_potential_present) THEN

      ALLOCATE (ppl_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_present",nkind*int_size)

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius",nkind*dp_size)
      ppl_radius = 0.0_dbl

    END IF

!   *** Calculate PBC coordinates ***

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atoms(ikind)%list)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           gth_potential=gth_potential,&
                           orb_basis_set=orb_basis_set,&
                           aux_basis_set=aux_basis_set,&
                           natom=natom(ikind),&
                           atom_list=atoms(ikind)%list)

      IF (ASSOCIATED(gth_potential)) THEN
        CALL get_potential(potential=gth_potential,&
                           ppl_present=ppl_present(ikind),&
                           ppl_radius=ppl_radius(ikind))  
      END IF


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
      IF(ASSOCIATED(orb_basis_set).OR.ASSOCIATED(aux_basis_set)) THEN

        IF (ASSOCIATED(orb_basis_set)) THEN

          orb_present(ikind) = .TRUE.

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 kind_radius=orb_radius(ikind))

        END IF

        IF (ASSOCIATED(aux_basis_set)) THEN

          aux_present(ikind) = .TRUE.

          CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                                 kind_radius=aux_radius(ikind))

        END IF

      ELSE

        orb_present(ikind) = .FALSE.
        aux_present(ikind) = .FALSE.

!        NULLIFY (pbc_coord(ikind)%r)
!        NULLIFY (pbc_coord(ikind)%s)

      END IF

    END DO

!   *** Build orbital GTH-PPL operator overlap list ***

    IF (gth_potential_present) THEN

      IF (ANY(ppl_present)) THEN

        CALL build_gth_neighbor_lists(r_skin)

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_kg_env(kg_env=kg_env,&
                        sac_ppl=sac_ppl)

!        IF (ionode.AND.globenv%print%sac_ppl_neighbor_lists) THEN
!         CALL write_neighbor_lists(sac_ppl,"ORBITAL GTH-PPL NEIGHBOR LISTS")
!       END IF
       
        CALL build_gth_neighbor_lists_aux(r_skin)

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_kg_env(kg_env=kg_env,&
                         sac_ppl_aux=sac_ppl_aux)

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

    IF (ALLOCATED(aux_present)) THEN
      DEALLOCATE (aux_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"aux_present")
    END IF

    IF (ALLOCATED(aux_radius)) THEN
      DEALLOCATE (aux_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"aux_radius")
    END IF

    IF (gth_potential_present) THEN

      DEALLOCATE (ppl_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_present")

      DEALLOCATE (ppl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius")

      IF (ALLOCATED(kind_b)) THEN
        DEALLOCATE (kind_b,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"kind_b")
      END IF

      IF (ALLOCATED(neighbors)) THEN
        DEALLOCATE (neighbors,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbors")
      END IF

    END IF

    CALL timestop(0.0_dbl,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_gth_neighbor_lists(r_skin)

!   Purpose: Build the neighbor lists for the calculation of the local part of
!            the GTH pseudo potential (PPL).

!   History: - Creation (18.03.2002,MK)

!   ***************************************************************************

    REAL(dbl), INTENT(IN) :: r_skin

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name =&
      "build_gth_neighbor_lists"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: ac,handle

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

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

        rac_max = orb_radius(ikind) + ppl_radius(kkind) + r_skin
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_dbl*subcells/sac_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,kkind,ngrid,&
            scaled_to_real(1.0_dbl/REAL(ngrid(:),dbl),cell)/unit_of_length
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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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
            CALL stop_memory(routine,"ijk",natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

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

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_dbl)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_dbl)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_dbl)*ngrid(3)))
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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE build_gth_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_gth_neighbor_lists_aux(r_skin)

!   Purpose: Build the neighbor lists for the calculation of the local part of
!            the GTH pseudo potential (PPL).

!   History: - Creation (18.03.2002,MK)

!   ***************************************************************************

    REAL(dbl), INTENT(IN) :: r_skin

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name =&
      "build_gth_neighbor_lists_aux"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: ac,handle

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    IF (ASSOCIATED(sac_ppl_aux)) THEN
      DO ac=1,SIZE(sac_ppl_aux)
        CALL deallocate_neighbor_list_set(sac_ppl_aux(ac)%neighbor_list_set)
      END DO
      DEALLOCATE (sac_ppl_aux,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sac_ppl_aux")
    END IF

    ALLOCATE (sac_ppl_aux(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sac_ppl_aux",0)

    DO ac=1,SIZE(sac_ppl_aux)
      NULLIFY (sac_ppl_aux(ac)%neighbor_list_set)
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

      IF (.NOT.aux_present(ikind)) CYCLE

      DO kkind=1,nkind

        ac = ikind + (kkind - 1)*nkind

        IF (.NOT.ppl_present(kkind)) CYCLE

!       *** Calculate the square of the maximum interaction distance ***

        rac_max = aux_radius(ikind) + ppl_radius(kkind) + r_skin
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_dbl*subcells/sac_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_cell_parameters) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,kkind,ngrid,&
            scaled_to_real(1.0_dbl/REAL(ngrid(:),dbl),cell)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          sac_ppl_aux(ac)%neighbor_list_set,&
                                        r_max=rac_max)
        neighbor_list_set => sac_ppl_aux(ac)%neighbor_list_set

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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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
            CALL stop_memory(routine,"ijk",natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

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

          DO iatom=1,natom(ikind)
            atom_a = atoms(ikind)%list(iatom)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_dbl)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_dbl)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_dbl)*ngrid(3)))
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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE build_gth_neighbor_lists_aux

! *****************************************************************************

  SUBROUTINE rebuild_neighbor_lists(kg_env,globenv)

!   Purpose: Rebuild all the required neighbor lists for Quickstep.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(kg_environment_type), INTENT(INOUT)  :: kg_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "rebuild_neighbor_lists"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set
    INTEGER                           :: handle,nparticle

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    group = globenv%group
    ionode = globenv%ionode
!   mype = globenv%mepos
!   npe = globenv%num_pe
    output_unit = globenv%scr

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set,&
                    sac_ppl_aux=sac_ppl_aux,&
                    sac_ppl=sac_ppl) 

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
                           orb_basis_set=orb_basis_set,&
                           aux_basis_set=aux_basis_set,&
                           natom=natom(ikind),&
                           atom_list=atoms(ikind)%list)

!      IF (ASSOCIATED(orb_basis_set).OR.ASSOCIATED(aux_basis_set)) THEN

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

!      ELSE

!        NULLIFY (pbc_coord(ikind)%r)
!        NULLIFY (pbc_coord(ikind)%s)

!      END IF

    END DO

!   *** Rebuild the orbital GTH-PPL operator overlap list ***

    IF (gth_potential_present) THEN

      IF (ASSOCIATED(sac_ppl)) THEN

        CALL rebuild_gth_neighbor_lists()

!        IF (ionode.AND.globenv%print%sac_ppl_neighbor_lists) THEN
!          CALL write_neighbor_lists(sac_ppl,"ORBITAL GTH-PPL NEIGHBOR LISTS")
!        END IF

      END IF

      IF (ASSOCIATED(sac_ppl_aux)) THEN
        CALL rebuild_gth_neighbor_lists_aux()
      END IF

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
      IF (ALLOCATED(neighbors)) THEN
        DEALLOCATE (neighbors,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbors")
      END IF
    END IF

    CALL timestop(0.0_dbl,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE rebuild_neighbor_lists

! *****************************************************************************

  SUBROUTINE rebuild_gth_neighbor_lists()

!   Purpose: Rebuild the neighbor lists for the calculation of the local part
!            of the GTH pseudo potential (PPL).

!   History: - Creation (17.09.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name =&
      "rebuild_gth_neighbor_lists"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: ac,handle

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

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
        ngrid(:) = MAX(1,NINT(0.5_dbl*subcells/sac_max(:)))

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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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
            CALL stop_memory(routine,"ijk",natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

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

          neighbor_list => first_list(neighbor_list_set)

          DO WHILE (ASSOCIATED(neighbor_list))
            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a)
            CALL init_neighbor_list(neighbor_list)
            iatom = atom_of_kind(atom_a)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_dbl)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_dbl)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_dbl)*ngrid(3)))
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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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

!MK          CALL clean_neighbor_list_set(neighbor_list_set)

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

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE rebuild_gth_neighbor_lists

! *****************************************************************************

  SUBROUTINE rebuild_gth_neighbor_lists_aux()

!   Purpose: Rebuild the neighbor lists for the calculation of the local part
!            of the GTH pseudo potential (PPL).

!   History: - Creation (17.09.2002,MK)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name =&
      "rebuild_gth_neighbor_lists_aux"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: ac,handle

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind
      DO kkind=1,nkind

        ac = ikind + (kkind - 1)*nkind

        neighbor_list_set => sac_ppl_aux(ac)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   r_max=rac_max)
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_dbl*subcells/sac_max(:)))

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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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
            CALL stop_memory(routine,"ijk",natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

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

          neighbor_list => first_list(neighbor_list_set)

          DO WHILE (ASSOCIATED(neighbor_list))
            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a)
            CALL init_neighbor_list(neighbor_list)
            iatom = atom_of_kind(atom_a)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_dbl)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_dbl)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_dbl)*ngrid(3)))
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

                  sc(:) = sc_pbc(:) + REAL(cell_c(:),dbl)
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

!MK          CALL clean_neighbor_list_set(neighbor_list_set)

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

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE rebuild_gth_neighbor_lists_aux

! *****************************************************************************

  SUBROUTINE write_neighbor_lists(neighbor_lists,name)

!   Purpose: Write a set of neighbor lists to the output unit.

!   History: - Creation (04.03.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: neighbor_lists
    CHARACTER(LEN=*), INTENT(IN)                          :: name

!   *** Local variables ***

    INTEGER :: ab,atom_a,atom_b,i,nneighbor

    REAL(dbl), DIMENSION(3) :: ra,rab,rb
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

END MODULE kg_neighbor_lists
