!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/nl [1.0] *
!!
!!   NAME
!!     nl
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists
!!
!!   AUTHOR
!!     Matthias Krack (08.10.1999)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE nl

  USE kinds, ONLY: wp => dp

  USE atoms,                 ONLY: atom_info,natom
  USE atomic_kinds,          ONLY: all_potential_present,&
                                   gth_potential_present,&
                                   kind_info,&
                                   nkind
  USE global_types,          ONLY: global_environment_type
  USE memory_utilities,      ONLY: reallocate
  USE method_specifications, ONLY: allchem,&
                                   gpw
  USE neighbor_list_types,   ONLY: allocate_neighbor_list_set,&
                                   deallocate_neighbor_list_set,&
                                   first_neighbor_list,&
                                   first_neighbor_node,&
                                   get_neighbor_node,&
                                   neighbor_list_set_type,&
                                   neighbor_list_type,&
                                   neighbor_node_type,&
                                   next_neighbor_list,&
                                   next_neighbor_node
  USE nl_utilities,          ONLY: build_neighbor_list,&
                                   write_neighbor_list_set
  USE termination,           ONLY: stop_memory
  USE timings,               ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

! *** Global variables ***

  TYPE(neighbor_list_set_type), POINTER :: aux_neighbor_list_set,&
                                           orb_neighbor_list_set,&
                                           verfc_neighbor_list_set,&
                                           vppl_neighbor_list_set,&
                                           vppnl_ac_neighbor_list_set,&
                                           vppnl_ca_neighbor_list_set
  LOGICAL                               :: init_pointers_done = .FALSE.

! *** Public variables ***

  PUBLIC :: aux_neighbor_list_set,&
            orb_neighbor_list_set,&
            verfc_neighbor_list_set,&
            vppl_neighbor_list_set,&
            vppnl_ca_neighbor_list_set

! *** Public subroutines ***

  PUBLIC :: build_neighbor_list_sets

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_neighbor_list_sets(globenv)

!   Purpose: Build all the required neighbor lists.

!   History: - Creation (28.08.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type) :: globenv

!   *** Local variables ***

    INTEGER :: handle,output_unit
    LOGICAL :: ionode

!   ---------------------------------------------------------------------------

    CALL timeset("build_neighbor_list_sets","I","",handle)

!   *** Initialize the neighbor list set pointers ***

    IF (.NOT.init_pointers_done) THEN
      NULLIFY (aux_neighbor_list_set)
      NULLIFY (orb_neighbor_list_set)
      NULLIFY (verfc_neighbor_list_set)
      NULLIFY (vppl_neighbor_list_set)
      NULLIFY (vppnl_ac_neighbor_list_set)
      NULLIFY (vppnl_ca_neighbor_list_set)
      init_pointers_done = .TRUE.
    END IF

    output_unit = globenv%scr
    ionode = globenv%ionode

!   *** Build the overlap neighbor list set ***

    CALL build_nl_overlap()

    IF (globenv%print%nl_overlap) THEN
      IF (ionode) THEN
        CALL write_neighbor_list_set(orb_neighbor_list_set,output_unit)
        IF (allchem) THEN
          CALL write_neighbor_list_set(aux_neighbor_list_set,output_unit)
        END IF
      END IF
    END IF

!   *** Build the core Hamiltonian neighbor list sets ***

    CALL build_nl_core_hamiltonian(globenv)

    IF (all_potential_present) THEN
      IF (ionode) THEN
        IF (gpw.AND.globenv%print%nl_verfc) THEN
          CALL write_neighbor_list_set(verfc_neighbor_list_set,output_unit)
        END IF
      END IF
    END IF

    IF (gth_potential_present) THEN
      IF (ionode) THEN
        IF (globenv%print%nl_vppl) THEN
          CALL write_neighbor_list_set(vppl_neighbor_list_set,output_unit)
        END IF
        IF (globenv%print%nl_vppnl) THEN
          CALL write_neighbor_list_set(vppnl_ac_neighbor_list_set,output_unit)
          CALL write_neighbor_list_set(vppnl_ca_neighbor_list_set,output_unit)
        END IF
      END IF
      CALL deallocate_neighbor_list_set(vppnl_ac_neighbor_list_set)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_neighbor_list_sets

! *****************************************************************************

  SUBROUTINE build_nl_overlap()

!   Purpose: Build a neighbor list based on the overlap of the Gaussian-type
!            basis functions.

!   History: - Creation (23.05.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER  :: atom_a,iatom,handle,ikind,jkind
    REAL(wp) :: radius_a,radius_b

!   ---------------------------------------------------------------------------

    CALL timeset("build_nl_overlap","I","",handle)

!   *** Overlap of the orbital basis functions ***
!   *** with the orbital basis functions       ***

    CALL allocate_neighbor_list_set(orb_neighbor_list_set,"ORB-OVERLAP")

    DO ikind=1,nkind
      radius_a = kind_info(ikind)%orb_basis_set%kind_radius
      DO iatom=1,kind_info(ikind)%natom
        atom_a = kind_info(ikind)%atom_list(iatom)
        DO jkind=1,nkind
          radius_b = kind_info(jkind)%orb_basis_set%kind_radius
          CALL build_neighbor_list(orb_neighbor_list_set,&
                                   atom_a,(/0,0,0/),radius_a,&
                                   kind_info(jkind),radius_b)
        END DO
      END DO
    END DO

!   *** Overlap of the orbital basis functions ***
!   *** with the auxiliary basis functions     ***

    IF (allchem) THEN
      CALL allocate_neighbor_list_set(aux_neighbor_list_set,"AUX-OVERLAP")
      DO ikind=1,nkind
        radius_a = kind_info(ikind)%orb_basis_set%kind_radius
        DO iatom=1,kind_info(ikind)%natom
          atom_a = kind_info(ikind)%atom_list(iatom)
          DO jkind=1,nkind
            radius_b = kind_info(jkind)%aux_basis_set%kind_radius
            CALL build_neighbor_list(aux_neighbor_list_set,&
                                     atom_a,(/0,0,0/),radius_a,&
                                     kind_info(jkind),radius_b)
          END DO
        END DO
      END DO
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_nl_overlap

! *****************************************************************************

  SUBROUTINE build_nl_core_hamiltonian(globenv)

!   Purpose: Build a neighbor list based on the overlap of the Gaussian-type
!            basis functions with the atomic potential operator.

!   History: - Creation (23.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_nl_core_hamiltonian (MODULE nl)"
    INTEGER, PARAMETER          :: nblock = 100

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list
    TYPE(neighbor_node_type), POINTER :: neighbor_node

    REAL(wp) :: radius_a,radius_c
    INTEGER  :: atom_a,atom_c,handle,iatom,iatom_prj_ppnl,ikind,istat,kkind,&
                natom_prj_ppnl
    LOGICAL  :: new_atom_prj_ppnl

    INTEGER, DIMENSION(3) :: cell

    INTEGER, DIMENSION(:,:), POINTER  :: atom_prj_ppnl

!   ---------------------------------------------------------------------------

    CALL timeset("build_nl_core_hamiltonian","I","",handle)

!   *** Allocate the required neighbor list data structures ***

    IF (gpw.AND.all_potential_present) THEN
      CALL allocate_neighbor_list_set(verfc_neighbor_list_set,"VERFC")
    END IF

    IF (gth_potential_present) THEN
      CALL allocate_neighbor_list_set(vppl_neighbor_list_set,"VPPL")
      CALL allocate_neighbor_list_set(vppnl_ac_neighbor_list_set,"VPPNL_AC")
      CALL allocate_neighbor_list_set(vppnl_ca_neighbor_list_set,"VPPNL_CA")
    END IF

!   *** Create the neighbor lists for all atoms of each atomic kind ***

    cell(:) = (/0,0,0/)

    DO ikind=1,nkind
      radius_a = kind_info(ikind)%orb_basis_set%kind_radius
      DO iatom=1,kind_info(ikind)%natom
        atom_a = kind_info(ikind)%atom_list(iatom)
        DO kkind=1,nkind
          IF (ASSOCIATED(kind_info(kkind)%all_potential)) THEN
            IF (gpw) THEN
              radius_c = kind_info(kkind)%core_charge_radius
              CALL build_neighbor_list(verfc_neighbor_list_set,&
                                       atom_a,cell(:),radius_a,&
                                       kind_info(kkind),radius_c)
            END IF
          ELSE IF (ASSOCIATED(kind_info(kkind)%gth_potential)) THEN
            IF ((kind_info(kkind)%gth_potential%nexp_ppl > 0).OR.allchem) THEN
              radius_c = kind_info(kkind)%gth_potential%ppl_radius
              CALL build_neighbor_list(vppl_neighbor_list_set,&
                                       atom_a,cell(:),radius_a,&
                                       kind_info(kkind),radius_c)
            END IF
            IF (kind_info(kkind)%gth_potential%nppnl > 0) THEN
              radius_c = kind_info(kkind)%gth_potential%ppnl_radius
              CALL build_neighbor_list(vppnl_ac_neighbor_list_set,&
                                       atom_a,cell(:),radius_a,&
                                       kind_info(kkind),radius_c)
            END IF
          END IF
        END DO
      END DO
    END DO

!   *** Create the neighbor list set "VPPNL_CA" ***
!   *** using the neighbor list set "VPPNL_AC"  ***

    IF (gth_potential_present) THEN

      NULLIFY (atom_prj_ppnl)

      atom_prj_ppnl => reallocate(atom_prj_ppnl,1,4,1,nblock)

      natom_prj_ppnl = 0

!     *** Build a list of all operator atoms "c" ***

      neighbor_list => first_neighbor_list(vppnl_ac_neighbor_list_set)

      DO WHILE (ASSOCIATED(neighbor_list))

        neighbor_node => first_neighbor_node(neighbor_list)

        DO WHILE (ASSOCIATED(neighbor_node))

          CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                 neighbor=atom_c,&
                                 cell=cell(:))

!         *** Check, if the operator atom is already in the list ***

          new_atom_prj_ppnl = .TRUE.

          DO iatom_prj_ppnl=1,natom_prj_ppnl
            IF ((atom_prj_ppnl(1,iatom_prj_ppnl) == atom_c).AND.&
                (atom_prj_ppnl(2,iatom_prj_ppnl) == cell(1)).AND.&
                (atom_prj_ppnl(3,iatom_prj_ppnl) == cell(2)).AND.&
                (atom_prj_ppnl(4,iatom_prj_ppnl) == cell(3))) THEN
              new_atom_prj_ppnl = .FALSE.
              EXIT
            END IF
          END DO

!         *** Add the operator atom, if it was not found in the list ***

          IF (new_atom_prj_ppnl) THEN
            natom_prj_ppnl = natom_prj_ppnl + 1
            IF (natom_prj_ppnl > SIZE(atom_prj_ppnl,2)) THEN
              atom_prj_ppnl =>&
                reallocate(atom_prj_ppnl,1,4,1,natom_prj_ppnl+nblock)
            END IF
            atom_prj_ppnl(1,natom_prj_ppnl) = atom_c
            atom_prj_ppnl(2:4,natom_prj_ppnl) = cell(:)
          END IF

          neighbor_node => next_neighbor_node(neighbor_node)

        END DO

        neighbor_list => next_neighbor_list(neighbor_list)

      END DO

      atom_prj_ppnl => reallocate(atom_prj_ppnl,1,4,1,natom_prj_ppnl)

!     *** Build the neighbor list set "VPPNL_CA" ***

      DO iatom_prj_ppnl=1,natom_prj_ppnl
        IF (globenv%mepos /= MODULO(iatom_prj_ppnl,globenv%num_pe)) CYCLE
        atom_c = atom_prj_ppnl(1,iatom_prj_ppnl)
        cell(:) = atom_prj_ppnl(2:4,iatom_prj_ppnl)
        kkind = atom_info(atom_c)%kind
        radius_c = kind_info(kkind)%gth_potential%ppnl_radius
        DO ikind=1,nkind
          radius_a = kind_info(ikind)%orb_basis_set%kind_radius
          CALL build_neighbor_list(vppnl_ca_neighbor_list_set,&
                                   atom_c,cell,radius_c,&
                                   kind_info(ikind),radius_a)
        END DO
      END DO

!     *** Release work storage ***

      DEALLOCATE (atom_prj_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_prj_ppnl")

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_nl_core_hamiltonian

! *****************************************************************************

END MODULE nl
