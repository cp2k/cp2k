!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_neighbor_lists [1.0] *
!!
!!   NAME
!!     kg_neighbor_lists
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists.
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_neighbor_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc,&
                                             real_to_scaled
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type,&
                                             qs_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_neighbor_lists,               ONLY: build_gth_ppl_neighbor_lists,&
                                             build_qs_neighbor_lists,&
                                             local_atoms_type,&
                                             write_neighbor_lists
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"  

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_neighbor_lists'

! *** Public subroutines ***

  PUBLIC :: build_kg_neighbor_lists

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_kg_neighbor_lists(kg_env,globenv,error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_kg_neighbor_lists'

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_control_type), POINTER           :: qs_control

    CALL write_checkpoint_information("entering "//routineN,globenv,error)

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control)
    qs_control => dft_control%qs_control

    SELECT CASE (qs_control%method)
    CASE ("KG_POL")
      CALL build_kg_neighbor_lists1(kg_env,globenv,kg_env%input,error)
    CASE ("KG_NOPOL")
      CALL build_kg_neighbor_lists1(kg_env,globenv,kg_env%input,error)
    CASE ("KG_GPW")
      CALL build_kg_neighbor_lists2(kg_env,globenv,kg_env%input,error)
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"method unknown")
    END SELECT

    CALL write_checkpoint_information("leaving "//routineN,globenv,error)

  END SUBROUTINE build_kg_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_kg_neighbor_lists1(kg_env,globenv,force_env_section,error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_kg_neighbor_lists1'

    CHARACTER(LEN=default_string_length)     :: print_key_path
    INTEGER :: atom_a, handle, iatom, iatom_local, ikind, istat, iw, maxatom, &
      maxatom_local, natom_a, natom_local_a, nkind, output_unit
    LOGICAL                                  :: elp_potential_present, &
                                                gth_potential_present, &
                                                print_subcell_grid
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: aux_present, orb_present, &
                                                ppl_present
    REAL(dp)                                 :: r_skin
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: aux_radius, orb_radius, &
                                                ppl_radius
    REAL(dp), DIMENSION(3)                   :: r_pbc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, orb_basis_set
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sac_ppl, sac_ppl_aux
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL write_checkpoint_information("entering "//routineN,globenv,error)
    CALL timeset(routineN,"I","",handle)

    NULLIFY (atomic_kind)
    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (distribution_1d)
    NULLIFY (gth_potential)
    NULLIFY (elp_potential)
    NULLIFY (orb_basis_set)
    NULLIFY (aux_basis_set)
    NULLIFY (particle_set)
    NULLIFY (sac_ppl)
    NULLIFY (sac_ppl_aux)
    print_subcell_grid = .FALSE.
    output_unit = cp_print_key_unit_nr(logger,force_env_section,"SUBSYS%PRINT%SUBCELL",&
         extension=".Log",error=error)
    IF (output_unit>0)  print_subcell_grid = .TRUE.

!   *** Set Verlet skin ***

    r_skin = 0.0_dp

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    local_particles=distribution_1d,&
                    particle_set=particle_set,&
                    sac_ppl=sac_ppl,&
                    sac_ppl_aux=sac_ppl_aux)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             elp_potential_present=elp_potential_present,&
                             gth_potential_present=gth_potential_present,&
                             maxatom=maxatom)

!   *** Allocate work storage ***

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom",nkind*int_size)

    ALLOCATE (orb_present(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "orb_present",nkind*int_size)

    ALLOCATE (orb_radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "orb_radius",nkind*dp_size)
    orb_radius(:) = 0.0_dp

    ALLOCATE (aux_present(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "aux_present",nkind*int_size)

    ALLOCATE (aux_radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "aux_radius",nkind*dp_size)
    aux_radius(:) = 0.0_dp

    IF (gth_potential_present.OR.elp_potential_present) THEN

      ALLOCATE (ppl_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "ppl_present",nkind*int_size)

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "ppl_radius",nkind*dp_size)
      ppl_radius = 0.0_dp

    END IF

!   *** Initialize the local data structures ***

    maxatom_local = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom(ikind)%list)
      NULLIFY (atom(ikind)%list_local_a_index)
      NULLIFY (atom(ikind)%list_local_b_index)
      NULLIFY (atom(ikind)%r_pbc)
      NULLIFY (atom(ikind)%s_pbc)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom(ikind)%list,&
                           gth_potential=gth_potential,&
                           elp_potential=elp_potential,&
                           aux_basis_set=aux_basis_set,&
                           orb_basis_set=orb_basis_set)

      natom_a = SIZE(atom(ikind)%list)
      natom_local_a = distribution_1d%n_el(ikind)

      maxatom_local = MAX(maxatom_local,natom_local_a)

      IF (ASSOCIATED(orb_basis_set)) THEN

        orb_present(ikind) = .TRUE.

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=orb_radius(ikind))

      ELSE

        orb_present(ikind) = .FALSE.

      END IF

      IF (ASSOCIATED(aux_basis_set)) THEN

        aux_present(ikind) = .TRUE.

        CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                               kind_radius=aux_radius(ikind))

      ELSE

        aux_present(ikind) = .FALSE.

      END IF
 
      IF (orb_present(ikind).OR.aux_present(ikind)) THEN 

        IF (natom_local_a > 0) THEN

          ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "atom(ikind)%list_local_a_index",&
                                           natom_local_a*int_size)

!         *** Build index vector for mapping ***

          DO iatom_local=1,natom_local_a
            atom_a = distribution_1d%list(ikind)%array(iatom_local)
            DO iatom=1,natom_a
              IF (atom_a == atom(ikind)%list(iatom)) THEN
                atom(ikind)%list_local_a_index(iatom_local) = iatom
                EXIT
              END IF
            END DO
          END DO

        END IF

        ALLOCATE (atom(ikind)%r_pbc(3,natom_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%r_pbc",3*natom_a*dp_size)

        ALLOCATE (atom(ikind)%s_pbc(3,natom_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%s_pbc",3*natom_a*dp_size)

!       *** Calculate PBC coordinates ***

        DO iatom=1,natom_a
          atom_a = atom(ikind)%list(iatom)
          r_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
          atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:),cell)
        END DO

      END IF

      IF (ASSOCIATED(gth_potential)) THEN
        CALL get_potential(potential=gth_potential,&
                           ppl_present=ppl_present(ikind),&
                           ppl_radius=ppl_radius(ikind)) 
      ELSE IF (ASSOCIATED(elp_potential)) THEN
        CALL get_potential(potential=elp_potential,&
                           ppl_present=ppl_present(ikind),&
                           ppl_radius=ppl_radius(ikind))  
      END IF

    END DO

!   *** Build orbital PPL operator overlap list ***

    IF (gth_potential_present.OR.elp_potential_present) THEN

      IF (ANY(ppl_present).AND.ANY(orb_present)) THEN

        CALL build_gth_ppl_neighbor_lists(sac_ppl,atom,cell,r_skin,maxatom_local,&
                                          print_subcell_grid,output_unit,&
                                          orb_present,ppl_present,&
                                          orb_radius,ppl_radius,&
                                          name="ORBITAL GTH-PPL",&
                                          particle_set=particle_set)

        CALL set_kg_env(kg_env=kg_env,sac_ppl=sac_ppl)

        print_key_path = "DFT%PRINT%NEIGHBOR_LISTS/SAC_PPL"

        IF (BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
                                             print_key_path,error=error),&
                  cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger=logger,&
                                    basis_section=force_env_section,&
                                    print_key_path=print_key_path,&
                                    extension=".out",&
                                    middle_name="orb_sac_ppl",&
                                    local=.TRUE.,&
                                    log_filename=.FALSE.,&
                                    file_position="REWIND",&
                                    error=error)
          CALL write_neighbor_lists(sac_ppl,particle_set,cell,globenv,iw,&
                                    "ORBITAL GTH-PPL",error)
          CALL cp_print_key_finished_output(unit_nr=iw,&
                                            logger=logger,&
                                            basis_section=force_env_section,&
                                            print_key_path=print_key_path,&
                                            local=.TRUE.,&
                                            error=error)
        END IF

      END IF

      IF (ANY(ppl_present).AND.ANY(aux_present)) THEN

        CALL build_gth_ppl_neighbor_lists(sac_ppl_aux,atom,cell,r_skin,maxatom_local,&
                                          print_subcell_grid,output_unit,&
                                          aux_present,ppl_present,&
                                          aux_radius,ppl_radius,&
                                          name="AUXILIARY ORBITAL GTH-PPL")

        CALL set_kg_env(kg_env=kg_env,sac_ppl_aux=sac_ppl_aux)

        print_key_path = "DFT%PRINT%NEIGHBOR_LISTS/SAC_PPL"

        IF (BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
                                             print_key_path,error=error),&
                  cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger=logger,&
                                    basis_section=force_env_section,&
                                    print_key_path=print_key_path,&
                                    extension=".out",&
                                    middle_name="aux_sac_ppl",&
                                    local=.TRUE.,&
                                    log_filename=.FALSE.,&
                                    file_position="REWIND",&
                                    error=error)
          CALL write_neighbor_lists(sac_ppl,particle_set,cell,globenv,iw,&
                                    "AUXILIARY ORBITAL GTH-PPL",error)
          CALL cp_print_key_finished_output(unit_nr=iw,&
                                            logger=logger,&
                                            basis_section=force_env_section,&
                                            print_key_path=print_key_path,&
                                            local=.TRUE.,&
                                            error=error)
        END IF

      END IF

    END IF

!   *** Release work storage ***

    DO ikind=1,nkind
      NULLIFY (atom(ikind)%list)
      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_a_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%list_local_a_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%list_local_b_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_b_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%list_local_b_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%r_pbc)) THEN
        DEALLOCATE (atom(ikind)%r_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%r_pbc")
      END IF
      IF (ASSOCIATED(atom(ikind)%s_pbc)) THEN
        DEALLOCATE (atom(ikind)%s_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%s_pbc")
      END IF
    END DO

    DEALLOCATE (atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom")

    DEALLOCATE (orb_present,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "orb_present")

    DEALLOCATE (orb_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "orb_radius")

    DEALLOCATE (aux_present,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "aux_present")

    DEALLOCATE (aux_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "aux_radius")

    IF (gth_potential_present.OR.elp_potential_present) THEN

      DEALLOCATE (ppl_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "ppl_present")

      DEALLOCATE (ppl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "ppl_radius")

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
                                      "SUBSYS%PRINT%SUBCELL",error=error)    
    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv,error)

  END SUBROUTINE build_kg_neighbor_lists1

! *****************************************************************************

  SUBROUTINE build_kg_neighbor_lists2(kg_env,globenv,force_env_section,error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_kg_neighbor_lists2'

    INTEGER                                  :: handle
    TYPE(qs_environment_type), POINTER       :: qs_env

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routineN,globenv,error)

    CALL timeset(routineN,"I","",handle)

    NULLIFY (qs_env)

    CALL get_kg_env(kg_env=kg_env,&
                    sub_qs_env=qs_env)

    IF(ASSOCIATED(qs_env)) THEN
      CALL build_qs_neighbor_lists(qs_env,globenv,.TRUE.,force_env_section,error=error)
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "sub_qs_env not defined")
    END IF

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv,error)

  END SUBROUTINE build_kg_neighbor_lists2

! *****************************************************************************

END MODULE kg_neighbor_lists

! *****************************************************************************
