!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****** cp2k/topology_util [1.0] *
!!
!!   NAME
!!     topology_util
!!
!!   FUNCTION
!!     Collection of subroutine needed for topology related things
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!    jgh (23-05-2004) Last atom of molecule information added
!!
!!   SOURCE
!******************************************************************************

MODULE topology_util
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE checkpoint_handler,              ONLY: Pwarning
  USE colvar_types,                    ONLY: colvar_create,&
                                             colvar_setup,&
                                             dist_colvar_id,&
                                             colvar_counters,&
                                             colvar_clone,&
                                             colvar_p_reallocate
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE external_potential_types,        ONLY: allocate_potential,&
                                             elp_potential_type,&
                                             get_potential,&
                                             set_potential
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_bondparm_covalent, do_bondparm_vdw, do_conn_off, do_fist, &
       do_qmmm_none, do_skip_ei_12, do_skip_ei_13, do_skip_ei_14, &
       do_skip_vdw_12, do_skip_vdw_13, do_skip_vdw_14
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       allocate_molecule_kind_set, atom_type, bend_type, bond_type, &
       colvar_constraint_type, fixd_constraint_type, &
       g3x3_constraint_type, g4x6_constraint_type, get_molecule_kind, &
       impr_type, molecule_kind_type, onfo_type, set_molecule_kind, &
       torsion_type, ub_type,&
       setup_colvar_counters
  USE molecule_types_new,              ONLY: &
       allocate_molecule_set, get_molecule, local_colvar_constraint_type, &
       local_constraint_type, local_g3x3_constraint_type, &
       local_g4x6_constraint_type, local_molecule_type, molecule_type, &
       set_molecule, set_molecule_set
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             generate_particle_distances,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: massunit
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_utilities,                ONLY: integer_to_string,&
                                             str_comp,&
                                             str_search,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             constraint_info_type,&
                                             topology_parameters_type
  USE util,                            ONLY: find_boundary
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_util'


  PRIVATE
  PUBLIC :: topology_set_atm_mass,&
            topology_reorder_atoms,&
            topology_connectivity_pack,&
            topology_constraint_pack,&
            topology_coordinate_pack,&
            topology_generate_bend,&
            topology_generate_bond,&
            topology_generate_dihe,&
            topology_generate_impr,&
            topology_generate_molecule,&
            topology_generate_onfo,&
            topology_generate_ub,&
            topology_molecules_check,&
            check_subsys_element

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_reorder_atoms
!!
!!
!!   FUNCTION
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_reorder_atoms(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_reorder_atoms', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      POINTER                                :: telement(:), &
                                                tlabel_atmname(:), &
                                                tlabel_resname(:)
    INTEGER                                  :: handle, iatm, iatm_count, &
                                                imol_count, inum, ityp, iw, &
                                                location, natom, nnum, ntyp, &
                                                stat
    INTEGER, POINTER                         :: mnum(:), mtyp(:), &
                                                new_position(:), &
                                                tmap_mol_num(:), &
                                                tmap_mol_typ(:)
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tr
    REAL(KIND=dp), POINTER                   :: tatm_charge(:), tatm_mass(:)
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    atom_info => topology%atom_info
    conn_info => topology%conn_info
    
    natom = topology%natoms
    
    NULLIFY(mtyp,mnum,new_position)
    NULLIFY(tlabel_resname,tlabel_atmname,telement)
    NULLIFY(tr,tatm_charge,tatm_mass)
    NULLIFY(tmap_mol_typ,tmap_mol_num)
    
    !ALLOCATE all the temporary arrays needed for this routine
    ALLOCATE(mtyp(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(mnum(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(new_position(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tlabel_resname(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tlabel_atmname(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tr(3,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tatm_charge(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tatm_mass(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmap_mol_typ(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmap_mol_num(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(telement(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iatm=1,natom
       mtyp(iatm) = atom_info%map_mol_typ(iatm)
       mnum(iatm) = atom_info%map_mol_num(iatm)
    END DO
    ntyp = MAXVAL(mtyp)
    nnum = MAXVAL(mnum)
    
    imol_count = 0
    iatm_count = 0
    
    !Let's determine the new order that we want everything to be in
    DO ityp=1,ntyp
       DO inum=1,nnum
          imol_count = imol_count + 1
          DO iatm=1,natom
             IF ( (ityp==mtyp(iatm)) .AND. (inum==mnum(iatm)) ) THEN
                iatm_count = iatm_count + 1
                new_position(iatm_count) = iatm
                CALL integer_to_string(imol_count,atom_info%label_resname(iatm))
                atom_info%label_resname(iatm) = "R"//atom_info%label_resname(iatm)
             END IF
          END DO
       END DO
    END DO
    
    !Lets swap the atoms now
    DO iatm=1,natom
       location = new_position(iatm)
       tlabel_resname(iatm)       = atom_info%label_resname(location)
       tlabel_atmname(iatm)       = atom_info%label_atmname(location)
       telement(iatm)             = atom_info%element(location)
       tr(1,iatm)                 = atom_info%r(1,location)
       tr(2,iatm)                 = atom_info%r(2,location)
       tr(3,iatm)                 = atom_info%r(3,location)
       tatm_charge(iatm)          = atom_info%atm_charge(location)
       tatm_mass(iatm)            = atom_info%atm_mass(location)
       tmap_mol_typ(iatm)         = atom_info%map_mol_typ(location)
       tmap_mol_num(iatm)         = atom_info%map_mol_num(location)
    END DO
    DO iatm=1,natom
       atom_info%label_resname(iatm)   = tlabel_resname(iatm)
       atom_info%label_atmname(iatm)   = tlabel_atmname(iatm)
       atom_info%element(iatm)         = telement(iatm)
       atom_info%r(1,iatm)             = tr(1,iatm)
       atom_info%r(2,iatm)             = tr(2,iatm)
       atom_info%r(3,iatm)             = tr(3,iatm)
       atom_info%atm_charge(iatm)      = tatm_charge(iatm)
       atom_info%atm_mass(iatm)        = tatm_mass(iatm)
       atom_info%map_mol_typ(iatm)     = tmap_mol_typ(iatm)
       atom_info%map_mol_num(iatm)     = tmap_mol_num(iatm)
    END DO
    
    !DEALLOCATE all the temporary arrays needed for this routine
    DEALLOCATE(mtyp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mnum,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(new_position,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tlabel_resname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tlabel_atmname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(telement,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tr,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tatm_charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !DEALLOCATE all the stuff in topology structure so it can be re-runned
    DEALLOCATE(atom_info%map_mol_typ,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(atom_info%map_mol_num,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%bond_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%bond_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%theta_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%theta_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%theta_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%phi_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%phi_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%phi_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%phi_d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_reorder_atoms

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_molecule
!!
!!
!!   FUNCTION
!!     Use information from bond list to generate molecule. (ie clustering)
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_molecule(topology,qmmm,qmmm_env,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_in_kind, atom_in_mol, first, handle, i, iatm, iatom, &
      ifirst, iiatom, iii, ilast, inum, itype, iw, j, jtype, jump1, jump2, &
      last, max_mol_num, myind, natom, nbond, nmol, nnum, ntype, stat
    INTEGER, DIMENSION(:), POINTER           :: checkme, qm_atom_index
    LOGICAL                                  :: do_again, failure, my_qmmm
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    NULLIFY(checkme, qm_atom_index)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    
    ntype = 1
    itype = 0
    nnum = 1
    inum = 0
    !
    ! QM/MM coordinate_control
    ! 
    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm
    
    natom = topology%natoms
    nbond = 0
    IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
    
    ALLOCATE(checkme(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF (ASSOCIATED(atom_info%map_mol_typ)) DEALLOCATE(atom_info%map_mol_typ,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (atom_info%map_mol_typ(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF (ASSOCIATED(atom_info%map_mol_num)) DEALLOCATE(atom_info%map_mol_num,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (atom_info%map_mol_num(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF (ASSOCIATED(atom_info%map_mol_res)) DEALLOCATE(atom_info%map_mol_res,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (atom_info%map_mol_res(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !Zero the arrays
    checkme(:) = 0
    atom_info%map_mol_typ(:) = 0
    atom_info%map_mol_num(:) = 0
    atom_info%map_mol_res(:) = 1
    
    !Parse the atoms list to see how many different molecule types there are
    ntype = 1
    atom_info%map_mol_typ(1)=1
    DO iatom=2,natom
       IF (topology%conn_type == do_conn_off ) THEN
          ntype=ntype+1
          atom_info%map_mol_typ(iatom)=ntype        
       ELSE
          IF(atom_info%label_molname(iatom-1)==atom_info%label_molname(iatom)) THEN
             atom_info%map_mol_typ(iatom)=ntype
          ELSE
             ntype=ntype+1
             atom_info%map_mol_typ(iatom)=ntype
          END IF
       END IF
    END DO
    
    DO itype=1,ntype
       !search for first atom in the type list and last to cut down time
       CALL find_boundary(atom_info%map_mol_typ,natom,first,last,itype)
       IF(iw>0) WRITE(iw,*) "      itype ",itype,first,last
       
       inum = 0

       IF (first /= last) THEN
          DO iatom=first,last
             checkme(:) = 0             
             IF(atom_info%map_mol_num(iatom)==0) THEN
                inum=inum+1
                checkme(iatom)=1
                DO WHILE (SUM(checkme)>0)
                   DO iiatom=first,last
                      IF((checkme(iiatom)==1)) THEN
                         checkme(iiatom)=0
                         atom_info%map_mol_num(iiatom) = inum
                         DO iii=1,nbond
                            IF((conn_info%bond_a(iii)==iiatom).AND.&
                                 (atom_info%map_mol_num(conn_info%bond_b(iii))==0)) THEN 
                               checkme(conn_info%bond_b(iii))=1
                            END IF
                            IF((conn_info%bond_b(iii)==iiatom).AND.&
                                 (atom_info%map_mol_num(conn_info%bond_a(iii))==0)) THEN
                               checkme(conn_info%bond_a(iii))=1
                            END IF
                         END DO
                      END IF
                   END DO
                END DO
             END IF             
          END DO
       ELSE          
          atom_info%map_mol_num(first) = 1
       END IF
    END DO
    
    DEALLOCATE(checkme,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF(topology%para_res) THEN
       DO itype=1,ntype
          CALL find_boundary(atom_info%map_mol_typ,natom,first,last,itype)
          nmol = atom_info%map_mol_num(last)
          DO jtype = 1,nmol
             CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,&
                  natom,first,last,itype,jtype)
             inum = 1
             atom_info%map_mol_res(first) = inum
             DO iatom=first+1,last
                IF((atom_info%label_resid(iatom-1)/=atom_info%label_resid(iatom)).OR.&
                     (atom_info%label_resname(iatom-1)/=atom_info%label_resname(iatom)))THEN
                   inum = inum + 1
                END IF
                atom_info%map_mol_res(iatom) = inum
             END DO
          END DO
       END DO
    END IF
    
    IF(iw>0) THEN
       DO iatom=1,natom
          WRITE(iw,*) "      iatom,map_mol_typ,map_mol_num ",iatom,&
               atom_info%map_mol_typ(iatom),atom_info%map_mol_num(iatom)
       END DO
    END IF

    IF (my_qmmm) THEN
       do_again = .FALSE.
       IF (iw>0) WRITE(iw,*)"MAP_MOL_NUM ",atom_info%map_mol_num
       IF (iw>0) WRITE(iw,*)"MAP_MOL_TYP ",atom_info%map_mol_typ
       ALLOCATE(qm_atom_index(SIZE(qmmm_env%qm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       qm_atom_index = qmmm_env%qm_atom_index
       CPPostcondition(ALL(qm_atom_index /= 0),cp_failure_level,routineP,error,failure)
       DO myind = 1, SIZE(qm_atom_index)
          IF (qm_atom_index(myind) ==0) CYCLE
          CALL find_boundary(atom_info%map_mol_typ,natom,ifirst,ilast,atom_info%map_mol_typ(qm_atom_index(myind)))
          CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,natom,ifirst,ilast,&
               atom_info%map_mol_typ(qm_atom_index(myind)),atom_info%map_mol_num(qm_atom_index(myind)))
          IF (iw>0) WRITE(iw,*)"qm fragment:: ifirst, ilast",ifirst,ilast
          CPPostcondition(((ifirst/=0).OR.(ilast/=natom)),cp_failure_level,routineP,error,failure)
          DO iatm = ifirst, ilast
             atom_info%label_molname(iatm) = "_QM_"//&
                  TRIM(atom_info%label_molname(iatm))
             IF (iw>0) WRITE(iw,*)"QM Molecule name :: ", atom_info%label_molname(iatm)
             WHERE (qm_atom_index == iatm) qm_atom_index = 0
          END DO
          DO iatm = 1, ifirst-1
             IF ( ANY(qm_atom_index == iatm) ) do_again = .TRUE.
          END DO
          DO iatm = ilast+1, natom
             IF ( ANY(qm_atom_index == iatm) ) do_again = .TRUE.
          END DO
          IF (iw>0) WRITE(iw,*)" Another QM fragment? :: ",do_again

          IF (ifirst /= 1) THEN
             jump1 = atom_info%map_mol_typ(ifirst)  - atom_info%map_mol_typ(ifirst-1)
             CPPostcondition(jump1<=1.AND.jump1>=0,cp_failure_level,routineP,error,failure)
             jump1 = ABS(jump1-1)
          ELSE
             jump1 = 0
          END IF
          IF (ilast /= natom) THEN
             jump2 = atom_info%map_mol_typ(ilast+1) - atom_info%map_mol_typ(ilast)
             CPPostcondition(jump2<=1.AND.jump2>=0,cp_failure_level,routineP,error,failure)
             jump2 = ABS(jump2-1)
          ELSE
             jump2 = 0
          END IF

          ! Changing mol_type consistently
          DO iatm = ifirst, natom
             atom_info%map_mol_typ(iatm) = atom_info%map_mol_typ(iatm) + jump1
          END DO
          DO iatm = ilast+1, natom
             atom_info%map_mol_typ(iatm) = atom_info%map_mol_typ(iatm) + jump2
          END DO

          IF (jump1 == 1) THEN
             DO iatm = ifirst, ilast
                atom_info%map_mol_num(iatm) = 1
             END DO
          END IF
          
          IF (jump2 == 1) THEN
             CALL find_boundary(atom_info%map_mol_typ,natom,first,last,atom_info%map_mol_typ(ilast+1))
             CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,natom,ifirst,ilast,&
                  atom_info%map_mol_typ(ilast+1),atom_info%map_mol_num(ilast+1))
             atom_in_mol = ilast-ifirst+1
             inum = 1
             DO iatm = first, last, atom_in_mol
                atom_info%map_mol_num(iatm:iatm+atom_in_mol-1) = inum
                inum = inum + 1
             END DO
          END IF

          IF (.NOT.do_again) EXIT
       END DO
       DEALLOCATE(qm_atom_index,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       IF(iw>0) THEN
          WRITE(iw,*)"After the QM/MM Setup:"
          DO iatom=1,natom
             WRITE(iw,*) "      iatom,map_mol_typ,map_mol_num ",iatom,&
                  atom_info%map_mol_typ(iatom),atom_info%map_mol_num(iatom)
          END DO
       END IF
    END IF
    !
    ! Further check : see if the number of atoms belonging to same molecule kinds
    !                 are equal
    ntype = MAXVAL(atom_info%map_mol_typ)
    IF (iw>0) WRITE(*,*)"SUMMARY:: Number of molecules type found:",ntype
    DO i = 1, ntype
       atom_in_kind   = COUNT(atom_info%map_mol_typ == i)
       IF (iw>0) WRITE(*,*)"Molecule type: ",i,"contains ",atom_in_kind,' atoms'
       IF (atom_in_kind <= 1) CYCLE
       CALL find_boundary(atom_info%map_mol_typ,natom,first,last,i)
       IF (iw>0) WRITE(*,*)"Boundary atoms: ",first, last
       CPPostcondition(last-first+1==atom_in_kind,cp_failure_level,routineP,error,failure)
       max_mol_num    = MAXVAL(atom_info%map_mol_num(first:last))
       IF (iw>0) WRITE(*,*)"Number of molecules of type ",i," is ::",max_mol_num       
       atom_in_mol    = atom_in_kind / max_mol_num
       IF (iw>0) WRITE(*,*)"Number of atoms per each molecule: ",atom_in_mol
       IF (iw>0) WRITE(*,*)"MAP_MOL_TYP:: ",atom_info%map_mol_typ(first:last)
       IF (iw>0) WRITE(*,*)"MAP_MOL_NUM:: ",atom_info%map_mol_num(first:last)
       !
       DO j = 1, max_mol_num
          IF (COUNT(atom_info%map_mol_num(first:last)==j) /= atom_in_mol) THEN
             IF(iw>0) THEN
                WRITE(iw,*)"molecule type:",i,"molecule num:",j," has ",COUNT(atom_info%map_mol_num(first:last)==j),&
                     " atoms instead of ",atom_in_mol," ."
             END IF
             CALL stop_program ("topology_generate_molecule","Two molecules of the same kind"//&
                  " have been created with different numbers of atoms!")
          END IF
       END DO
    END DO
    
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE topology_generate_molecule

!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_set_atm_mass
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Use info from periodic table and set atm_mass
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_set_atm_mass(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_set_atm_mass', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    CHARACTER(LEN=2)                         :: upper_sym_1, upper_sym_2
    CHARACTER(LEN=2), POINTER                :: atm_symbol(:)
    INTEGER                                  :: handle, i, iatom, ielem, &
                                                ielem_found, iw, natom, &
                                                nelem, stat
    LOGICAL                                  :: element_found, failure
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    atom_info => topology%atom_info
    
    natom = topology%natoms
    
    !Get the temporary array that we need
    ALLOCATE(atm_symbol(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !------------------------------------------------------------------------------
    DO iatom=1,natom
       atm_symbol(iatom) = ''
       IF(ASSOCIATED(atom_info%element)) THEN
          atm_symbol(iatom) = atom_info%element(iatom)
       END IF
    END DO
    !------------------------------------------------------------------------------
    !Set up the arrays that we need for this routine
    DO iatom=1,natom
       nelem = UBOUND(ptable,1)
       element_found=.FALSE.
       ielem_found=0
       upper_sym_1=atm_symbol(iatom)
       CALL uppercase(upper_sym_1)
       DO ielem=1,nelem 
          upper_sym_2=ptable(ielem)%symbol
          CALL uppercase(upper_sym_2)
          IF(upper_sym_1 == upper_sym_2)THEN
             element_found=.TRUE.
             ielem_found=ielem
          END IF
       END DO
       IF (.NOT. element_found) THEN
          IF (iw>0) WRITE(iw,'(A)') "WARNING : topology_set_atm_mass: The element " &
               //atm_symbol(iatom)//" was not found in the periodic table"
       ENDIF
       
       atom_info%atm_mass(iatom) = ptable(ielem_found)% amass
       IF(iw>0) WRITE(iw,'(7X,A,A5)') "In topology_set_atm_mass :: iatom = ",ptable(ielem_found)%symbol
    END DO

    DEALLOCATE(atm_symbol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    
  END SUBROUTINE topology_set_atm_mass
!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_bond
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Use info from periodic table and assumptions to generate bonds
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_bond(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bond', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: upper_sym_1, upper_sym_2
    CHARACTER(LEN=2), POINTER                :: atm_symbol(:)
    INTEGER                                  :: cbond, handle, i, iatom, &
                                                ibond, ielem, ielem_found, &
                                                iw, jatom, natom, nelem, stat
    LOGICAL                                  :: element_found, failure
    REAL(KIND=dp), POINTER                   :: radius(:)
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    
    cbond = 0
    
    natom = topology%natoms
    
    NULLIFY(atm_symbol)
    NULLIFY(radius)
    
    !Get the temporary array that we need
    ALLOCATE(atm_symbol(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(radius(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !------------------------------------------------------------------------------
    DO iatom=1,natom
       atm_symbol(iatom) = ''
       IF(ASSOCIATED(atom_info%element)) THEN
          atm_symbol(iatom) = atom_info%element(iatom)
       END IF
    END DO
    !------------------------------------------------------------------------------
    !Set up the arrays that we need for this routine
    DO iatom=1,natom
       nelem = UBOUND(ptable,1)
       element_found=.FALSE.
       ielem_found=0
       upper_sym_1=atm_symbol(iatom)
       CALL uppercase(upper_sym_1)
       DO ielem=1,nelem 
          upper_sym_2=ptable(ielem)%symbol
          CALL uppercase(upper_sym_2)
          IF(upper_sym_1 == upper_sym_2)THEN
             element_found=.TRUE.
             ielem_found=ielem
          END IF
       END DO
       IF (.NOT. element_found) THEN
          IF (iw>0) WRITE(iw,'(A)') "WARNING : topology_generate_bond: The element " &
               //atm_symbol(iatom)//" was not found in the periodic table"
       ENDIF
       
       IF(topology%bondparm_type==do_bondparm_covalent) THEN
          radius(iatom) = ptable(ielem_found)%covalent_radius
       ELSE IF(topology%bondparm_type==do_bondparm_vdw) THEN
          radius(iatom) = ptable(ielem_found)%vdw_radius
       ELSE
          CALL stop_program ("topology_generate_bond","illegal bondparm_type")
       END IF
       CALL convert_to_cp2k_units("ANGSTROM",length=radius(iatom))
       
       atom_info%atm_mass(iatom) = ptable(ielem_found)% amass
       IF(iw>0) WRITE(iw,'(7X,A,A5)') "In topology_generate_bond :: iatom = ",ptable(ielem_found)%symbol
    END DO
    !------------------------------------------------------------------------------
    IF( (topology%bondparm_type==do_bondparm_covalent).OR.&
         (topology%bondparm_type==do_bondparm_vdw) ) THEN
       CALL generate_particle_distances(atom_info%r,topology%cell,radius,&
            topology%bondparm_factor,conn_info%bond_a,conn_info%bond_b,&
            topology%bondparm_type,globenv)
       IF (iw>0) THEN
          WRITE(iw,'(A,I6)')"Total number of generated bonds :: ",SIZE(conn_info%bond_a)
          ! Dump info on the computed bonds
          DO ibond = 1, SIZE(conn_info%bond_a)
             WRITE(iw,'("BOND Nr.:",I6," between atoms:",2I6)')ibond,conn_info%bond_a(ibond),&
                  conn_info%bond_b(ibond)
          END DO
       END IF
    ELSE 
       CALL stop_program ("topology_generate_bond","illegal bondparm_type")
    END IF
    
    IF(topology%para_res) THEN
       DO ibond=1,SIZE(conn_info%bond_a)
          iatom = conn_info%bond_a(ibond)
          jatom = conn_info%bond_b(ibond)
          IF((atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)).OR.&
               (atom_info%label_resid(iatom)/=atom_info%label_resid(jatom)).OR.&
               (atom_info%label_resname(iatom)/=atom_info%label_resname(jatom)))THEN
             IF(iw>0) WRITE(iw,*) "      PARA_RES, bond between molecules atom ",&
                  iatom,jatom
             cbond = cbond + 1
             CALL reallocate(conn_info%c_bond_a,1,cbond+1)
             CALL reallocate(conn_info%c_bond_b,1,cbond+1)
             conn_info%c_bond_a(cbond) = iatom
             conn_info%c_bond_b(cbond) = jatom
          ELSE
             IF(atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)) THEN
                CALL stop_program ("topology_generate_bond",&
                     "bonds between different molecule types??? ")
             END IF
          END IF
       END DO
    END IF
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    DEALLOCATE(atm_symbol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)

  END SUBROUTINE topology_generate_bond

!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_bend
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of bonds, generate a list of bends
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_bend(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bend', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    INTEGER                                  :: handle, ibond, itheta, iw, &
                                                jbond, jtheta, nbond, ntheta
    LOGICAL                                  :: failure
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    conn_info => topology%conn_info
    nbond = 0
    !IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
    IF(ASSOCIATED(conn_info%bond_a)) THEN
      nbond = SIZE(conn_info%bond_a)
    ELSE
      CALL reallocate(conn_info%bond_a,1,nbond)
      CALL reallocate(conn_info%bond_b,1,nbond)
    END IF
    ntheta = 0
    CALL reallocate(conn_info%theta_a,1,ntheta+nblock)
    CALL reallocate(conn_info%theta_b,1,ntheta+nblock)
    CALL reallocate(conn_info%theta_c,1,ntheta+nblock)
    !------------------------------------------------------------------------------
    ntheta=0
    DO ibond=1,nbond
       DO jbond=ibond+1,nbond
          !Test case A
          IF( (conn_info%bond_b(ibond)==conn_info%bond_a(jbond)) ) THEN
             ntheta = ntheta + 1
             IF(ntheta > SIZE(conn_info%theta_a)) THEN
                CALL reallocate(conn_info%theta_a,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_b,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_c,1,ntheta+nblock)
             END IF
             conn_info%theta_a(ntheta) =  conn_info%bond_a(ibond)
             conn_info%theta_b(ntheta) =  conn_info%bond_b(ibond)
             conn_info%theta_c(ntheta) =  conn_info%bond_b(jbond)
             CYCLE
          END IF
          !Test case B
          IF( (conn_info%bond_b(ibond)==conn_info%bond_b(jbond)) ) THEN
             ntheta = ntheta + 1
             IF(ntheta > SIZE(conn_info%theta_a)) THEN
                CALL reallocate(conn_info%theta_a,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_b,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_c,1,ntheta+nblock)
             END IF
             conn_info%theta_a(ntheta) =  conn_info%bond_a(ibond)
             conn_info%theta_b(ntheta) =  conn_info%bond_b(ibond)
             conn_info%theta_c(ntheta) =  conn_info%bond_a(jbond)
             CYCLE
          END IF
          !Test case c
          IF( (conn_info%bond_a(ibond)==conn_info%bond_a(jbond)) ) THEN
             ntheta = ntheta + 1
             IF(ntheta > SIZE(conn_info%theta_a)) THEN
                CALL reallocate(conn_info%theta_a,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_b,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_c,1,ntheta+nblock)
             END IF
             conn_info%theta_a(ntheta) =  conn_info%bond_b(ibond)
             conn_info%theta_b(ntheta) =  conn_info%bond_a(ibond)
             conn_info%theta_c(ntheta) =  conn_info%bond_b(jbond)
             CYCLE
          END IF
          !Test case A
          IF( (conn_info%bond_a(ibond)==conn_info%bond_b(jbond)) ) THEN
             ntheta = ntheta + 1
             IF(ntheta > SIZE(conn_info%theta_a)) THEN
                CALL reallocate(conn_info%theta_a,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_b,1,ntheta+nblock)
                CALL reallocate(conn_info%theta_c,1,ntheta+nblock)
             END IF
             conn_info%theta_a(ntheta) =  conn_info%bond_b(ibond)
             conn_info%theta_b(ntheta) =  conn_info%bond_a(ibond)
             conn_info%theta_c(ntheta) =  conn_info%bond_a(jbond)
             CYCLE
          END IF
       END DO
    END DO
    !------------------------------------------------------------------------------
    
    CALL reallocate(conn_info%theta_a,1,ntheta)
    CALL reallocate(conn_info%theta_b,1,ntheta)
    CALL reallocate(conn_info%theta_c,1,ntheta)
    
    !check for duplicates
    DO itheta=1,ntheta
       DO jtheta=itheta+1,ntheta
          IF( (conn_info%theta_a(itheta) == conn_info%theta_a(jtheta)) .AND. &
               (conn_info%theta_b(itheta) == conn_info%theta_b(jtheta)) .AND. &
               (conn_info%theta_c(itheta) == conn_info%theta_c(jtheta)) ) THEN
             IF(iw>0) WRITE(iw,*) "    Duplicate BEND ",itheta,&
                  conn_info%theta_a(itheta),&
                  conn_info%theta_b(itheta),&
                  conn_info%theta_c(itheta)
             IF(iw>0) WRITE(iw,*) "    Duplicate BEND ",jtheta,&
                  conn_info%theta_a(jtheta),&
                  conn_info%theta_b(jtheta),&
                  conn_info%theta_c(jtheta)
             CALL stop_program("topology_generate_bend","duplicate bend")
          END IF
       END DO
    END DO

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_generate_bend

!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_ub
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of bend, generate a list of Urey-Bradley
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_ub(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_ub', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, itheta, iw, ntheta
    LOGICAL                                  :: failure
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    conn_info => topology%conn_info
    
    ntheta = SIZE(conn_info%theta_a)
    CALL reallocate(conn_info%ub_a,1,ntheta)
    CALL reallocate(conn_info%ub_b,1,ntheta)
    CALL reallocate(conn_info%ub_c,1,ntheta)
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    DO itheta=1,ntheta
       conn_info%ub_a(itheta) = conn_info%theta_a(itheta)
       conn_info%ub_b(itheta) = conn_info%theta_b(itheta)
       conn_info%ub_c(itheta) = conn_info%theta_c(itheta)
    END DO

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)

  END SUBROUTINE topology_generate_ub

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_dihe
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of bend, generate a list of dihe
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_dihe(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_dihe', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    INTEGER                                  :: handle, iphi, itheta, iw, &
                                                jphi, jtheta, nphi, ntheta
    LOGICAL                                  :: failure
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    conn_info => topology%conn_info
    ntheta = SIZE(conn_info%theta_a)
    nphi = 0
    CALL reallocate(conn_info%phi_a,1,nphi+nblock)
    CALL reallocate(conn_info%phi_b,1,nphi+nblock)
    CALL reallocate(conn_info%phi_c,1,nphi+nblock)
    CALL reallocate(conn_info%phi_d,1,nphi+nblock)
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    DO itheta=1,ntheta
       DO jtheta=itheta+1,ntheta
          !Test case A
          IF( (conn_info%theta_a(itheta)==conn_info%theta_b(jtheta)) .AND. &
               (conn_info%theta_b(itheta)==conn_info%theta_a(jtheta)) ) THEN
             nphi = nphi + 1
             IF(nphi > SIZE(conn_info%phi_a)) THEN
                CALL reallocate(conn_info%phi_a,1,nphi+nblock)
                CALL reallocate(conn_info%phi_b,1,nphi+nblock)
                CALL reallocate(conn_info%phi_c,1,nphi+nblock)
                CALL reallocate(conn_info%phi_d,1,nphi+nblock)
             END IF
             conn_info%phi_a(nphi) = conn_info%theta_c(itheta)
             conn_info%phi_b(nphi) = conn_info%theta_b(itheta)
             conn_info%phi_c(nphi) = conn_info%theta_a(itheta)
             conn_info%phi_d(nphi) = conn_info%theta_c(jtheta)
             CYCLE
          END IF
          !Test case B
          IF( (conn_info%theta_a(itheta)==conn_info%theta_b(jtheta)) .AND. &
               (conn_info%theta_b(itheta)==conn_info%theta_c(jtheta)) ) THEN
             nphi = nphi + 1
             IF(nphi > SIZE(conn_info%phi_a)) THEN
                CALL reallocate(conn_info%phi_a,1,nphi+nblock)
                CALL reallocate(conn_info%phi_b,1,nphi+nblock)
                CALL reallocate(conn_info%phi_c,1,nphi+nblock)
                CALL reallocate(conn_info%phi_d,1,nphi+nblock)
             END IF
             conn_info%phi_a(nphi) = conn_info%theta_c(itheta)
             conn_info%phi_b(nphi) = conn_info%theta_b(itheta)
             conn_info%phi_c(nphi) = conn_info%theta_a(itheta)
             conn_info%phi_d(nphi) = conn_info%theta_a(jtheta)
             CYCLE
          END IF
          !Test case C
          IF( (conn_info%theta_b(itheta)==conn_info%theta_a(jtheta)) .AND. &
               (conn_info%theta_c(itheta)==conn_info%theta_b(jtheta)) ) THEN
             nphi = nphi + 1
             IF(nphi > SIZE(conn_info%phi_a)) THEN
                CALL reallocate(conn_info%phi_a,1,nphi+nblock)
                CALL reallocate(conn_info%phi_b,1,nphi+nblock)
                CALL reallocate(conn_info%phi_c,1,nphi+nblock)
                CALL reallocate(conn_info%phi_d,1,nphi+nblock)
             END IF
             conn_info%phi_a(nphi) = conn_info%theta_a(itheta)
             conn_info%phi_b(nphi) = conn_info%theta_b(itheta)
             conn_info%phi_c(nphi) = conn_info%theta_c(itheta)
             conn_info%phi_d(nphi) = conn_info%theta_c(jtheta)
             CYCLE
          END IF
          !Test case D
          IF( (conn_info%theta_b(itheta)==conn_info%theta_c(jtheta)) .AND. &
               (conn_info%theta_c(itheta)==conn_info%theta_b(jtheta)) ) THEN
             nphi = nphi + 1
             IF(nphi > SIZE(conn_info%phi_a)) THEN
                CALL reallocate(conn_info%phi_a,1,nphi+nblock)
                CALL reallocate(conn_info%phi_b,1,nphi+nblock)
                CALL reallocate(conn_info%phi_c,1,nphi+nblock)
                CALL reallocate(conn_info%phi_d,1,nphi+nblock)
             END IF
             conn_info%phi_a(nphi) = conn_info%theta_a(itheta)
             conn_info%phi_b(nphi) = conn_info%theta_b(itheta)
             conn_info%phi_c(nphi) = conn_info%theta_c(itheta)
             conn_info%phi_d(nphi) = conn_info%theta_a(jtheta)
             CYCLE
          END IF
       END DO
    END DO
    !------------------------------------------------------------------------------
    CALL reallocate(conn_info%phi_a,1,nphi)
    CALL reallocate(conn_info%phi_b,1,nphi)
    CALL reallocate(conn_info%phi_c,1,nphi)
    CALL reallocate(conn_info%phi_d,1,nphi)
    
    ! check for duplicate dihe. Can yank this out once the subroutine is DEBUGED
    DO iphi=1,nphi
       DO jphi=iphi+1,nphi
          IF( (conn_info%phi_a(iphi)==conn_info%phi_a(jphi)) .AND. &
               (conn_info%phi_b(iphi)==conn_info%phi_b(jphi)) .AND. &
               (conn_info%phi_c(iphi)==conn_info%phi_c(jphi)) .AND. &
               (conn_info%phi_d(iphi)==conn_info%phi_d(jphi)) ) THEN
             IF(iw>0) THEN
                WRITE(iw,*) "      Dihe ",iphi,&
                     conn_info%phi_a(iphi),conn_info%phi_b(iphi),&  
                     conn_info%phi_c(iphi),conn_info%phi_d(iphi)
                WRITE(iw,*) "      Dihe ",jphi,&
                     conn_info%phi_a(jphi),conn_info%phi_b(jphi),&  
                     conn_info%phi_c(jphi),conn_info%phi_d(jphi)
             END IF
             CALL stop_program("topology_generate_dihe","duplicate dihe")
          END IF
       END DO
    END DO

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    
  END SUBROUTINE topology_generate_dihe
 
!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_impr
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of bends, generate a list of impr
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_impr(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_impr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, handle, iatom, &
                                                ibond, iw, natom, nbond, &
                                                nimpr, stat
    INTEGER, DIMENSION(:), POINTER           :: checkme
    LOGICAL                                  :: failure
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    
    natom  = SIZE(atom_info%label_atmname)
    nbond = 0
    IF ( ASSOCIATED ( conn_info % bond_a ) ) &
         nbond = SIZE(conn_info%bond_a)
    nimpr  = 0
    
    CALL reallocate(conn_info%impr_a,1,nimpr)
    CALL reallocate(conn_info%impr_b,1,nimpr)
    CALL reallocate(conn_info%impr_c,1,nimpr)
    CALL reallocate(conn_info%impr_d,1,nimpr)

    DO iatom=1,natom
       NULLIFY(checkme)
       counter=0
       DO ibond=1,nbond
          IF(iatom==conn_info%bond_a(ibond)) THEN
             counter=counter+1
             CALL reallocate(checkme,1,counter)
             checkme(counter)=conn_info%bond_b(ibond)
          END IF
          IF(iatom==conn_info%bond_b(ibond)) THEN
             counter=counter+1
             CALL reallocate(checkme,1,counter)
             checkme(counter)=conn_info%bond_a(ibond)
          END IF
       END DO
       
       IF(counter==3) THEN
          nimpr=nimpr+1
          CALL reallocate(conn_info%impr_a,1,nimpr)
          CALL reallocate(conn_info%impr_b,1,nimpr)
          CALL reallocate(conn_info%impr_c,1,nimpr)
          CALL reallocate(conn_info%impr_d,1,nimpr)
          conn_info%impr_a(nimpr) = iatom
          conn_info%impr_b(nimpr) = checkme(1)
          conn_info%impr_c(nimpr) = checkme(2)
          conn_info%impr_d(nimpr) = checkme(3)
       END IF
       
       IF(counter/=0) THEN 
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)

  END SUBROUTINE topology_generate_impr

!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_onfo
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of torsion, generate a list of onfo
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_onfo(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_onfo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iphi, iw, nphi
    LOGICAL                                  :: failure
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    conn_info => topology%conn_info    
    nphi = SIZE(conn_info%phi_a)
    CALL reallocate(conn_info%onfo_a,1,nphi)
    CALL reallocate(conn_info%onfo_b,1,nphi)
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    DO iphi=1,nphi
       conn_info%onfo_a(iphi) = conn_info%phi_a(iphi)
       conn_info%onfo_b(iphi) = conn_info%phi_d(iphi)
    END DO
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    
  END SUBROUTINE topology_generate_onfo

!!*****
!******************************************************************************
!!****** topology/topology_connectivity_pack [1.0] *
!!
!!   NAME
!!     topology_connectivity_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     topology connectivity pack
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE topology_connectivity_pack(molecule_kind_set,molecule_set,&
       topology,globenv,subsys_section,error)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_connectivity_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: atm_a, atm_b, atm_c, atm_d, c_intra, counter, first, handle, &
      i, ibond, ikind, imol, inum, ires, ityp, iw, j, k, last, natom, nbond, &
      nelectron, nhcopt, nimpr, nmol, nonfo, nphi, nsgf, ntheta, nub, stat
    INTEGER, DIMENSION(:), POINTER           :: checkme, first_list, &
                                                last_list, molecule_list
    LOGICAL                                  :: failure, found
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list
    TYPE(local_molecule_type), POINTER       :: lmi
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(onfo_type), DIMENSION(:), POINTER   :: onfo_list
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list
    TYPE(ub_type), DIMENSION(:), POINTER     :: ub_list

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Set the topology%[nmol_type,nmol,nmol_conn]
    !-----------------------------------------------------------------------------
    natom  = topology%natoms
    topology%nmol      = 1
    topology%nmol_type = 1
    topology%nmol_conn = 0
    DO i=1,natom-1
       IF((atom_info%map_mol_typ(i+1)/=atom_info%map_mol_typ(i)) .OR. &
            (atom_info%map_mol_res(i+1)/=atom_info%map_mol_res(i))) &
            topology%nmol_type = topology%nmol_type + 1
       IF((atom_info%map_mol_typ(i+1)/=atom_info%map_mol_typ(i)) .OR. &
            (atom_info%map_mol_num(i+1)/=atom_info%map_mol_num(i)) .OR. &
            (atom_info%map_mol_res(i+1)/=atom_info%map_mol_res(i))) &
            topology%nmol = topology%nmol + 1
       IF((atom_info%map_mol_typ(i+1)==atom_info%map_mol_typ(i)) .AND. &
            (atom_info%map_mol_num(i+1)==atom_info%map_mol_num(i)) .AND. &
            (atom_info%map_mol_res(i+1)/=atom_info%map_mol_res(i))) &
            topology%nmol_conn = topology%nmol_conn + 1
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Allocate the molecule_kind_set
    !-----------------------------------------------------------------------------
    IF(topology%nmol_type <= 0) THEN
       CALL stop_program("topology_connectivity_pack","no molecule kind defined")
    ELSE
       NULLIFY(molecule_kind_set)
       i = topology%nmol_type
       CALL allocate_molecule_kind_set(molecule_kind_set,i)
       IF(iw>0) WRITE(iw,*) "    Allocated molecule_kind_set, Dimenstion of ",&
            SIZE(molecule_kind_set)
    END IF
        
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Allocate the molecule_set
    !-----------------------------------------------------------------------------
    IF(topology%nmol <= 0) THEN
       CALL stop_program("topology_connectivity_pack","no molecule defined")
    ELSE
       NULLIFY(molecule_set)
       i = topology%nmol
       CALL allocate_molecule_set(molecule_set,i)
       IF(iw>0) WRITE(iw,*) "    Allocated molecule_set, dimenstion of ",&
            topology%nmol
    END IF
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the molecule_kind_set%[kind_number,name,nhcopt,nsgf,nelectron]
    !-----------------------------------------------------------------------------
    counter=0
    nhcopt = topology % nhcopt
    natom = topology%natoms
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       nsgf=0
       nelectron=0
       i=0
       found = .FALSE.
       DO j=1,natom-1
          IF((atom_info%map_mol_typ(j+1)/=atom_info%map_mol_typ(j)) .OR. &
               (atom_info%map_mol_res(j+1)/=atom_info%map_mol_res(j))) THEN
             i = i + 1
             IF(ikind==i) THEN
                found = .TRUE.
                EXIT
             END IF
          END IF
       END DO
       name = TRIM(atom_info%label_molname(j))
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=ikind,&
            molname_generated=topology%molname_generated,&
            name=TRIM(name),&
            nhcopt=nhcopt,&
            nsgf=nsgf,&
            nelectron=nelectron)
    END DO



    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Set the molecule_list for molecule_kind in molecule_kind_set
    !-----------------------------------------------------------------------------
    counter=1
    natom = topology%natoms
    DO ikind=1,SIZE(molecule_kind_set)
       i=-1
       imol=0
       ityp=0
       inum=0
       ires=0
       found = .FALSE.
       DO j=1,natom
          IF(topology%para_res) THEN
             IF((atom_info%map_mol_typ(j) /= ityp) .OR. &
                  (atom_info%map_mol_num(j) /= inum) .OR. &
                  (atom_info%map_mol_res(j) /= ires)) imol = imol + 1
             IF((atom_info%map_mol_typ(j) /= ityp) .OR. &
                  (atom_info%map_mol_res(j) /= ires)) i = i + 1
             ityp = atom_info%map_mol_typ(j)
             inum = atom_info%map_mol_num(j)
             ires = atom_info%map_mol_res(j)
             IF(ikind==i) THEN
                found = .TRUE.
                EXIT
             END IF
          ELSE
             IF((atom_info%map_mol_typ(j) /= ityp) .OR. &
                  (atom_info%map_mol_num(j) /= inum)) imol = imol + 1
             IF((atom_info%map_mol_typ(j) /= ityp)) i = i + 1
             ityp = atom_info%map_mol_typ(j)
             inum = atom_info%map_mol_num(j)
             ires = atom_info%map_mol_res(j)
             IF(ikind==i) THEN
                found = .TRUE.
                EXIT
             END IF
          END IF
       END DO
       IF(ikind==topology%nmol_type) imol=topology%nmol+1
       ALLOCATE(molecule_list(imol-counter),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i=1,SIZE(molecule_list)
          molecule_list(i)=i+counter-1
       END DO
       molecule_kind => molecule_kind_set(ikind)
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list)
       IF(iw>0) WRITE(iw,*) "      molecule_list",ikind,molecule_list(:)
       counter=imol
    END DO
    DO ikind=1,SIZE(molecule_kind_set)-topology%nmol_type
       ALLOCATE(molecule_list(1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       molecule_list(1)=imol+ikind
       molecule_kind => molecule_kind_set(ikind+topology%nmol_type)
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list)
       IF(iw>0) WRITE(iw,*) "      molecule_list",ikind+topology%nmol_type,molecule_list(:)
    END DO



    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the molecule_set(imol)%molecule_kind via set_molecule
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list)
       DO i=1,SIZE(molecule_list)
          molecule => molecule_set(molecule_list(i))
          CALL set_molecule(molecule,molecule_kind=molecule_kind)
       END DO
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Set the molecule_set(imol)%[first_atom,last_atom] via set_molecule_set
    !-----------------------------------------------------------------------------
    ALLOCATE(first_list(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(last_list(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    first_list(:) = 0
    last_list(:) = 0
    DO ikind=1,topology%nmol
       i=0
       imol=0
       ityp=0
       inum=0
       ires=0
       found = .FALSE.
       counter=0
       DO j=1,natom
          IF((atom_info%map_mol_typ(j) /= ityp) .OR. &
               (atom_info%map_mol_num(j) /= inum) .OR. &
               (atom_info%map_mol_res(j) /= ires)) THEN
             ityp = atom_info%map_mol_typ(j)
             inum = atom_info%map_mol_num(j)
             ires = atom_info%map_mol_res(j)
             i = i + 1
             imol = imol + 1
             IF(ikind==i) THEN
                counter=imol
                found = .TRUE.
                EXIT
             END IF
          END IF
       END DO
       IF(ikind==topology%nmol_type) imol=topology%nmol
       first_list(ikind)=j
    END DO
    DO ikind=1,topology%nmol-1
       last_list(ikind)=first_list(ikind+1)-1
    END DO
    last_list(topology%nmol)=topology%natoms
    CALL set_molecule_set(molecule_set,first_list,last_list)
    DEALLOCATE(first_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Set and NULLIFY the molecule_set(imol)%lmi via set_molecule_set
    !-----------------------------------------------------------------------------
    DO i=1,SIZE(molecule_set)
       molecule => molecule_set(i)
       NULLIFY(lmi)
       ALLOCATE(lmi,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lmi%nstates=0
       NULLIFY(lmi%states)
       CALL set_molecule(molecule,lmi=lmi)
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Set the atom_list for molecule_kind in molecule_kind_set (PART 1)
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,SIZE(molecule_set)
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          CALL get_molecule(molecule=molecule,&
               first_atom=first,last_atom=last)
          natom = 0
          IF(first/=0 .AND. last/=0) natom = last-first+1
          ALLOCATE(atom_list(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,natom
             !Atomic kind information will be filled in (PART 2)
             NULLIFY(atom_list(i)%atomic_kind)
             atom_list(i)%name=atom_info%label_atmname(i+first-1)
             IF(iw>0) WRITE(iw,'(5X,A,I5,I5,1X,A5)') "atom_list ",ikind,i,atom_list(i)%name
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,atom_list=atom_list)
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Set the molecule_kind%[nbond,bond_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,topology%nmol
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          c_intra = 0      !number of extra connector intra-mol interactions
          counter = i
          nbond = 0
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
          ALLOCATE(checkme(nbond),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          checkme(:)=0
          DO j=1,nbond
             atm_a = conn_info%bond_a(j)
             atm_b = conn_info%bond_b(j)
             IF((first<=atm_a .AND. atm_a <=last) .AND.&
                  (first<=atm_b .AND. atm_b <=last)) checkme(j)=1
          END DO
          IF(ASSOCIATED(conn_info%c_bond_a)) THEN
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) c_intra = c_intra + 1
             END DO
          END IF
          IF(iw>0) THEN 
             WRITE(iw,*) "    Total number bonds for molecule type ",SUM(checkme)+c_intra
             WRITE(iw,*) "    checkme (bonds inside  molecules) :: ",checkme
             WRITE(iw,*) "    c_intra (bonds between molecules) :: ",c_intra
          END IF
          ibond = SUM(checkme) + c_intra
          ALLOCATE(bond_list(ibond),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ibond = 0
          DO j=1,nbond
             IF(checkme(j)==1) THEN
                ibond = ibond + 1
                bond_list(ibond)%a=conn_info%bond_a(j)-first+1
                bond_list(ibond)%b=conn_info%bond_b(j)-first+1
                bond_list(ibond)%name="BOND"
                !point this to the right bond_kind_type if using force field
                NULLIFY(bond_list(ibond)%bond_kind)
                IF(iw>0) THEN
                   WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",&
                        ikind,"  bond",&
                        conn_info%bond_a(j),&
                        conn_info%bond_b(j),&
                        "offset number at",&
                        conn_info%bond_a(j)-first+1,&
                        conn_info%bond_b(j)-first+1
                END IF
             END IF
          END DO
          IF(c_intra/=0) THEN
             ibond = SIZE(bond_list) - c_intra + 1
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) THEN
                   bond_list(ibond)%a=conn_info%c_bond_a(j)-first+1
                   bond_list(ibond)%b=conn_info%c_bond_b(j)-first+1
                   bond_list(ibond)%name="BOND"
                   !point this to the right bond_kind_type if using force field
                   NULLIFY(bond_list(ibond)%bond_kind)
                   IF(iw>0) THEN
                      WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",&
                           ikind,"c_bond",&
                           conn_info%c_bond_a(j),&
                           conn_info%c_bond_b(j),&
                           "offset number at",&
                           conn_info%c_bond_a(j)-first+1,&
                           conn_info%c_bond_b(j)-first+1
                   END IF
                   ibond = ibond + 1
                END IF
             END DO
          END IF
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               nbond=SIZE(bond_list),bond_list=bond_list)
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Set the molecule_kind%[nbend,bend_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,topology%nmol
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          ntheta = 0
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          IF(ASSOCIATED(conn_info%theta_a)) ntheta = SIZE(conn_info%theta_a)
          ALLOCATE(checkme(ntheta),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          checkme(:)=0
          DO j=1,ntheta
             atm_a = conn_info%theta_a(j)
             atm_b = conn_info%theta_b(j)
             atm_c = conn_info%theta_c(j)
             IF((first<=atm_a .AND. atm_a <=last) .AND.&
                  (first<=atm_b .AND. atm_b <=last) .AND.&
                  (first<=atm_c .AND. atm_c <=last)) checkme(j)=1
          END DO
          IF(ASSOCIATED(conn_info%c_bond_a)) THEN
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) THEN
                   DO k=1,ntheta
                      IF( ( (atm_a==conn_info%theta_a(k)) .AND.&
                           (atm_b==conn_info%theta_b(k)) ) .OR.&
                           ( (atm_a==conn_info%theta_b(k)) .AND.&
                           (atm_b==conn_info%theta_a(k)) ) .OR.&
                           ( (atm_a==conn_info%theta_b(k)) .AND.&
                           (atm_b==conn_info%theta_c(k)) ) .OR.&
                           ( (atm_a==conn_info%theta_c(k)) .AND.&
                           (atm_b==conn_info%theta_b(k)) ) ) THEN
                         checkme(k)=1
                      END IF
                   END DO
                END IF
             END DO
          END IF

          IF(iw>0) WRITE(iw,*) "    Total number bends for molecule type ",&
               SUM(checkme)
          i = SUM(checkme)
          ALLOCATE(bend_list(i),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          i = 0
          DO j=1,ntheta
             IF(checkme(j)==1) THEN
                i = i + 1
                bend_list(i)%a=conn_info%theta_a(j)-first+1
                bend_list(i)%b=conn_info%theta_b(j)-first+1
                bend_list(i)%c=conn_info%theta_c(j)-first+1
                bend_list(i)%name="BEND"
                !point this to the right bend_kind_type if using force field
                NULLIFY(bend_list(i)%bend_kind)
                IF(iw>0) THEN
                   WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                        "molecule_kind",ikind,"bend",&
                        conn_info%theta_a(j),&
                        conn_info%theta_b(j),&
                        conn_info%theta_c(j),&
                        "offset number at",&
                        conn_info%theta_a(j)-first+1,&
                        conn_info%theta_b(j)-first+1,&
                        conn_info%theta_c(j)-first+1
                END IF
             END IF
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               nbend=SIZE(bend_list),bend_list=bend_list)
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 12. Set the molecule_kind%[nub,ub_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,topology%nmol
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          nub = 0
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          IF(ASSOCIATED(conn_info%ub_a)) nub = SIZE(conn_info%ub_a)
          ALLOCATE(checkme(nub),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          checkme(:)=0
          DO j=1,nub
             atm_a = conn_info%ub_a(j)
             atm_b = conn_info%ub_b(j)
             atm_c = conn_info%ub_c(j)
             IF((first<=atm_a .AND. atm_a <=last) .AND.&
                  (first<=atm_b .AND. atm_b <=last) .AND.&
                  (first<=atm_c .AND. atm_c <=last)) checkme(j)=1
          END DO
          IF(ASSOCIATED(conn_info%c_bond_a)) THEN
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) THEN
                   DO k=1,nub
                      IF( ( (atm_a==conn_info%ub_a(k)) .AND.&
                           (atm_b==conn_info%ub_b(k)) ) .OR.&
                           ( (atm_a==conn_info%ub_b(k)) .AND.&
                           (atm_b==conn_info%ub_a(k)) ) .OR.&
                           ( (atm_a==conn_info%ub_b(k)) .AND.&
                           (atm_b==conn_info%ub_c(k)) ) .OR.&
                           ( (atm_a==conn_info%ub_c(k)) .AND.&
                           (atm_b==conn_info%ub_b(k)) ) ) THEN
                         checkme(k)=1
                      END IF
                   END DO
                END IF
             END DO
          END IF

          IF(iw>0) WRITE(iw,*) "    Total number UB for molecule type ",&
               SUM(checkme)
          i = SUM(checkme)
          ALLOCATE(ub_list(i),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          i = 0
          DO j=1,nub
             IF(checkme(j)==1) THEN
                i = i + 1
                ub_list(i)%a=conn_info%ub_a(j)-first+1
                ub_list(i)%b=conn_info%ub_b(j)-first+1
                ub_list(i)%c=conn_info%ub_c(j)-first+1
                ub_list(i)%name="UB"
                !point this to the right ub_kind_type if using force field
                NULLIFY(ub_list(i)%ub_kind)
                IF(iw>0) THEN
                   WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                        "molecule_kind",ikind,"UB",&
                        conn_info%ub_a(j),&
                        conn_info%ub_b(j),&
                        conn_info%ub_c(j),&
                        "offset number at",&
                        conn_info%ub_a(j)-first+1,&
                        conn_info%ub_b(j)-first+1,&
                        conn_info%ub_c(j)-first+1
                END IF
             END IF
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               nub=SIZE(ub_list),ub_list=ub_list)
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 13. Set the molecule_kind%[ntorsion,torsion_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,topology%nmol
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          nphi = 0
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          IF(ASSOCIATED(conn_info%phi_a)) nphi = SIZE(conn_info%phi_a)
          ALLOCATE(checkme(nphi),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          checkme(:)=0
          DO j=1,nphi
             atm_a = conn_info%phi_a(j)
             atm_b = conn_info%phi_b(j)
             atm_c = conn_info%phi_c(j)
             atm_d = conn_info%phi_d(j)
             IF((first<=atm_a .AND. atm_a <=last) .AND.&
                  (first<=atm_b .AND. atm_b <=last) .AND.&
                  (first<=atm_c .AND. atm_c <=last) .AND.&
                  (first<=atm_d .AND. atm_d <=last)) checkme(j)=1
          END DO
          IF(ASSOCIATED(conn_info%c_bond_a)) THEN
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) THEN
                   DO k=1,nphi
                      IF( ( (atm_a==conn_info%phi_a(k)) .AND.&
                           (atm_b==conn_info%phi_b(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_b(k)) .AND.&
                           (atm_b==conn_info%phi_c(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_c(k)) .AND.&
                           (atm_b==conn_info%phi_d(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_d(k)) .AND.&
                           (atm_b==conn_info%phi_c(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_c(k)) .AND.&
                           (atm_b==conn_info%phi_b(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_b(k)) .AND.&
                           (atm_b==conn_info%phi_a(k)) ) ) THEN
                         checkme(k)=1
                      END IF
                   END DO
                END IF
             END DO
          END IF

          IF(iw>0) WRITE(iw,*) "    Total number TOR for molecule type ",&
               SUM(checkme)
          i = SUM(checkme)
          ALLOCATE(torsion_list(i),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          i = 0
          DO j=1,nphi
             IF(checkme(j)==1) THEN
                i = i + 1
                torsion_list(i)%a=conn_info%phi_a(j)-first+1
                torsion_list(i)%b=conn_info%phi_b(j)-first+1
                torsion_list(i)%c=conn_info%phi_c(j)-first+1
                torsion_list(i)%d=conn_info%phi_d(j)-first+1
                torsion_list(i)%name="TORSION"
                !point this to the right torsion_kind_type if using force field
                NULLIFY(torsion_list(i)%torsion_kind)
                IF(iw>0) THEN
                   WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                        "molecule_kind",ikind,"TOR",&
                        conn_info%phi_a(j),&
                        conn_info%phi_b(j),&
                        conn_info%phi_c(j),&
                        conn_info%phi_d(j),&
                        "offset number at",&
                        conn_info%phi_a(j)-first+1,&
                        conn_info%phi_b(j)-first+1,&
                        conn_info%phi_c(j)-first+1,&
                        conn_info%phi_d(j)-first+1
                END IF
             END IF
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               ntorsion=SIZE(torsion_list),torsion_list=torsion_list)
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 14. Set the molecule_kind%[nimpr,impr_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,topology%nmol
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          nimpr = 0
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          IF(ASSOCIATED(conn_info%impr_a)) nimpr = SIZE(conn_info%impr_a)
          ALLOCATE(checkme(nimpr),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          checkme(:)=0
          DO j=1,nimpr
             atm_a = conn_info%impr_a(j)
             atm_b = conn_info%impr_b(j)
             atm_c = conn_info%impr_c(j)
             atm_d = conn_info%impr_d(j)
             IF((first<=atm_a .AND. atm_a <=last) .AND.&
                  (first<=atm_b .AND. atm_b <=last) .AND.&
                  (first<=atm_c .AND. atm_c <=last) .AND.&
                  (first<=atm_d .AND. atm_d <=last)) checkme(j)=1
          END DO
          IF(ASSOCIATED(conn_info%c_bond_a)) THEN
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) THEN
                   DO k=1,nimpr
                      IF( ( (atm_a==conn_info%impr_a(k)) .AND.&
                           (atm_b==conn_info%impr_b(k)) ) .OR.&
                           ( (atm_a==conn_info%impr_a(k)) .AND.&
                           (atm_b==conn_info%impr_c(k)) ) .OR.&
                           ( (atm_a==conn_info%impr_a(k)) .AND.&
                           (atm_b==conn_info%impr_d(k)) ) .OR.&
                           ( (atm_a==conn_info%impr_b(k)) .AND.&
                           (atm_b==conn_info%impr_a(k)) ) .OR.&
                           ( (atm_a==conn_info%impr_c(k)) .AND.&
                           (atm_b==conn_info%impr_a(k)) ) .OR.&
                           ( (atm_a==conn_info%impr_d(k)) .AND.&
                           (atm_b==conn_info%impr_a(k)) ) ) THEN
                         checkme(k)=1
                      END IF
                   END DO
                END IF
             END DO
          END IF

          IF(iw>0) WRITE(iw,*) "    Total number IMPR for molecule type ",&
               SUM(checkme)
          i = SUM(checkme)
          ALLOCATE(impr_list(i),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          i = 0
          DO j=1,nimpr
             IF(checkme(j)==1) THEN
                i = i + 1
                impr_list(i)%a=conn_info%impr_a(j)-first+1
                impr_list(i)%b=conn_info%impr_b(j)-first+1
                impr_list(i)%c=conn_info%impr_c(j)-first+1
                impr_list(i)%d=conn_info%impr_d(j)-first+1
                impr_list(i)%name="IMPROPER"
                !point this to the right impr_kind_type if using force field
                NULLIFY(impr_list(i)%impr_kind)
                IF(iw>0) THEN
                   WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                        "molecule_kind",ikind,"IMPR",&
                        conn_info%impr_a(j),&
                        conn_info%impr_b(j),&
                        conn_info%impr_c(j),&
                        conn_info%impr_d(j),&
                        "offset number at",&
                        conn_info%impr_a(j)-first+1,&
                        conn_info%impr_b(j)-first+1,&
                        conn_info%impr_c(j)-first+1,&
                        conn_info%impr_d(j)-first+1
                END IF
             END IF
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               nimpr=SIZE(impr_list),impr_list=impr_list)
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 15. Set the molecule_kind%[nonfo,onfo_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    counter = 0
    DO ikind=1,topology%nmol
       molecule => molecule_set(ikind)
       molecule_kind => molecule_set(ikind)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          nonfo = 0
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          IF(ASSOCIATED(conn_info%onfo_a)) nonfo = SIZE(conn_info%onfo_a)
          ALLOCATE(checkme(nonfo),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          checkme(:)=0
          DO j=1,nonfo
             atm_a = conn_info%onfo_a(j)
             atm_b = conn_info%onfo_b(j)
             IF((first<=atm_a .AND. atm_a <=last) .AND.&
                  (first<=atm_b .AND. atm_b <=last)) checkme(j)=1
          END DO
          IF(ASSOCIATED(conn_info%c_bond_a)) THEN
             DO j=1,SIZE(conn_info%c_bond_a)
                atm_a = conn_info%c_bond_a(j)
                atm_b = conn_info%c_bond_b(j)
                IF(first<=atm_a .AND. atm_a <=last) THEN
                   DO k=1,nonfo
                      IF( ( (atm_a==conn_info%phi_a(k)) .AND.&
                           (atm_b==conn_info%phi_b(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_b(k)) .AND.&
                           (atm_b==conn_info%phi_c(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_c(k)) .AND.&
                           (atm_b==conn_info%phi_d(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_d(k)) .AND.&
                           (atm_b==conn_info%phi_c(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_c(k)) .AND.&
                           (atm_b==conn_info%phi_b(k)) ) .OR.&
                           ( (atm_a==conn_info%phi_b(k)) .AND.&
                           (atm_b==conn_info%phi_a(k)) ) ) THEN
                         checkme(k)=1
                      END IF
                   END DO
                END IF
             END DO
          END IF

          IF(iw>0) WRITE(iw,*) "    Total number ONFO for molecule type ",&
               SUM(checkme)
          i = SUM(checkme)
          ALLOCATE(onfo_list(i),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          i = 0
          DO j=1,nonfo
             IF(checkme(j)==1) THEN
                i = i + 1
                onfo_list(i)%a=conn_info%onfo_a(j)-first+1
                onfo_list(i)%b=conn_info%onfo_b(j)-first+1
                onfo_list(i)%name="ONFO"
                NULLIFY(onfo_list(i)%onfo_kind)
                IF(iw>0) THEN
                   WRITE(iw,'(7X,A,I3,1X,A,I4,I4,1X,A,I4,I4)') &
                        "molecule_kind",ikind,"ONFO",&
                        conn_info%onfo_a(j),&
                        conn_info%onfo_b(j),&
                        "offset number at",&
                        conn_info%onfo_a(j)-first+1,&
                        conn_info%onfo_b(j)-first+1
                END IF
             END IF
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               nonfo=SIZE(onfo_list),onfo_list=onfo_list)
          DEALLOCATE(checkme,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_connectivity_pack

!******************************************************************************
!!****** topology_util/topology_coordinate_pack [1.0] *
!!
!!   NAME
!!     topology_coordinate_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Take info readin from different file format and stuff it into
!!     compatible data structure in cp2k
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino - modified in order to optimize the list of molecules 
!!                     to build the exclusion lists
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

  SUBROUTINE topology_coordinate_pack (particle_set,atomic_kind_set,&
       molecule_kind_set,molecule_set,topology,globenv,qmmm,qmmm_env,&
       subsys_section,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_coordinate_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname, &
                                                my_elp_description(2), &
                                                my_elp_name
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: element
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: work(:)
    INTEGER :: counter, first, first_atom, handle, i, iatom, ibend, ibond, &
      ikind, imol_now, imolecule, itorsion, iw, j, k, handle3, kk,&
      last, last_atom, natom, nbend, nbond, ntorsion, position, stat, handle2
    INTEGER, DIMENSION(:), POINTER :: atomkind, iatomlist, first_mol, kind_of, &
      last_mol, list, list2, molecule_list, natom_of_kind
    LOGICAL                                  :: failure, found, my_qmmm
    LOGICAL, POINTER, DIMENSION(:)           :: atomlist, atomlist_ei, excl_molecule_list
    REAL(KIND=dp)                            :: alpha_core_charge, qeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge, mass
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list, ex_bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list, ex_bond_list
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: ex_torsion_list, torsion_list

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
        
    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm=qmmm
    atom_info => topology%atom_info

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Determine topology%[natom_type,atom_names] and save mass(natom_type) 
    !    and element(natom_type) 
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_1',handle2)
    counter=0
    NULLIFY(work,mass,element,charge)
    ALLOCATE(work(topology%natoms),mass(topology%natoms),element(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(charge(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    work(:)=""
    IF (iw>0) WRITE(iw,*)"molecule_kind_set ::",SIZE(molecule_kind_set)
    DO i = 1, SIZE(molecule_kind_set)
       DO j=1,SIZE(molecule_set)
          molecule => molecule_set(j)
          molecule_kind => molecule_set(j)%molecule_kind
          CALL get_molecule_kind(molecule_kind=molecule_kind,kind_number=ikind)
          IF(ikind==i) EXIT
       END DO
       IF (iw>0) WRITE(iw,*)"molecule number ::",j," has molecule kind number ::",i
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            natom=natom,atom_list=atom_list)
       CALL get_molecule(molecule=molecule,&
            first_atom=first,last_atom=last)
       IF (iw>0) WRITE(iw,*)"boundaries of molecules (first, last) ::",first,last
       DO j = 1, natom
          position=str_search(work,counter,atom_list(j)%name)
          IF(position==0)THEN
             counter = counter + 1
             work(counter)   = atom_list(j)%name
             mass(counter)   = atom_info%atm_mass(first+j-1)
             element(counter)= atom_info%element(first+j-1)
             charge(counter) = atom_info%atm_charge(first+j-1)
             IF(iw>0) WRITE(iw,'(7X,A,1X,A5,F10.5,5X,A2,5X,F10.5)') &
                  "NEW ATOMIC KIND",work(counter),mass(counter),element(counter),charge(counter)
          ELSE
             found = .FALSE.
             DO k = 1, counter
                IF( (str_comp(work(k),atom_list(j)%name)) .AND. &
                     (charge(k)==atom_info%atm_charge(first+j-1)) ) THEN
                   found = .TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter = counter + 1
                work(counter)   = atom_list(j)%name
                mass(counter)   = atom_info%atm_mass(first+j-1)
                element(counter)= atom_info%element(first+j-1)
                charge(counter) = atom_info%atm_charge(first+j-1)
                IF(iw>0) WRITE(iw,'(7X,A,1X,A5,F10.5,5X,A2,5X,F10.5)') &
                     "NEW ATOMIC KIND",work(counter),mass(counter),element(counter),charge(counter)
             END IF
          END IF
       END DO
    END DO
    topology%natom_type = counter
    ALLOCATE(atom_info%atom_names(topology%natom_type),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    atom_info%atom_names(1:counter)=work(1:counter)
    DEALLOCATE(work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL reallocate(mass,1,counter)
    CALL reallocate(element,1,counter)
    CALL reallocate(charge,1,counter)
    IF(iw>0) &
         WRITE(iw,'(5X,A,I3)') "Total Number of Atomic Kinds = ",topology%natom_type
    CALL timestop(0.0_dp,handle2)


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Allocate the data structure for the atomic kind information 
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_2',handle2)
    NULLIFY(atomic_kind_set)
    CALL allocate_atomic_kind_set(atomic_kind_set,topology%natom_type)
    CALL timestop(0.0_dp,handle2)


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3.  Allocate the data structure for the atomic information
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_3',handle2)
    NULLIFY(particle_set)
    CALL allocate_particle_set(particle_set,topology%natoms)
    CALL timestop(0.0_dp,handle2)


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the atomic_kind_set(ikind)%[name,kind_number,mass]
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_4',handle2)
    DO i=1,topology%natom_type
       atomic_kind => atomic_kind_set(i)
       mass(i) = mass(i) * massunit
       CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=i,&
            name=atom_info%atom_names(i),element_symbol=element(i),&
            mass=mass(i))
       IF (iw>0) THEN
          WRITE(iw,'(A,I5,A,I5,4A)')"Atomic Kind n.:",i," out of:",topology%natom_type,&
               " name:   ",TRIM(atom_info%atom_names(i)),"   element:   ",&
               TRIM(element(i))
       END IF
    END DO
    DEALLOCATE(mass,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(element,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Determine number of atom of each kind (ie natom_of_kind and kind_of)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_5',handle2)
    ALLOCATE(kind_of(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(natom_of_kind(topology%natom_type),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    kind_of(:) = 0
    natom_of_kind(:) = 0
    DO i=1,topology%natom_type
       DO j=1,topology%natoms
          IF( (atom_info%atom_names(i)==atom_info%label_atmname(j)) .AND. &
               (charge(i)==atom_info%atm_charge(j)) ) THEN
             natom_of_kind(i)=natom_of_kind(i)+1
             IF(kind_of(j)==0) kind_of(j)=i
          END IF
       END DO
    END DO
    IF (ANY(kind_of==0)) THEN
      WRITE(*,*)"Two molecules have been defined as identical molecules but atoms mismatch charges!!"
      CPPostcondition(.false.,cp_failure_level,routineP,error,failure)
    END IF 
    CALL timestop(0.0_dp,handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the atom_kind_set(ikind)%[natom,atom_list]
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_6',handle2)
    DO i=1,topology%natom_type
       atomic_kind => atomic_kind_set(i)
       NULLIFY (iatomlist)
       ALLOCATE(iatomlist(natom_of_kind(i)),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       counter=0
       DO j=1,topology%natoms
          IF(kind_of(j)==i) THEN
             counter=counter+1
             iatomlist(counter)=j
          END IF
       END DO
       IF(iw>0) THEN 
          WRITE(iw,'(A,I6,A)') "      Atomic kind ",i," contains particles"
          DO J=1,SIZE(iatomlist)  
             IF (MOD(J,5).EQ.0) THEN ! split long lines
                WRITE(iw,'(I12)') iatomlist(J)
             ELSE
                WRITE(iw,'(I12)',ADVANCE="NO") iatomlist(J)
             ENDIF
          ENDDO
          WRITE(iw,*)
       ENDIF
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
            natom=natom_of_kind(i),atom_list=iatomlist)
       DEALLOCATE(iatomlist,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(natom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Fill in some more infor for particle_set and lpi (local part info)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_7',handle2)
    DO i=1,topology%natoms
       ikind = kind_of(i)
       IF (iw>0) THEN
          WRITE(iw,*)"atom number :: ",i,"kind number ::",ikind
       END IF
       particle_set(i)%atomic_kind => atomic_kind_set(ikind)
       particle_set(i)%r(:) = atom_info%r(:,i)
       NULLIFY(particle_set(i)%lpi)
       ALLOCATE(particle_set(i)%lpi)
       NULLIFY(particle_set(i)%lpi%exclusion)
       ALLOCATE(particle_set(i)%lpi%exclusion(SIZE(atomic_kind_set)))
       DO k=1,SIZE(atomic_kind_set)
          NULLIFY(particle_set(i)%lpi%exclusion(k)%list_vdw)
          NULLIFY(particle_set(i)%lpi%exclusion(k)%list_ei)
       END DO
    END DO
    DEALLOCATE(kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Fill in the particle_set%lpi%exclusion%list_vdw
    ! 9. Fill in the particle_set%lpi%exclusion%list_ei (if required)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_89',handle2)
    ALLOCATE(atomkind(SIZE(particle_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(first_mol(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(last_mol(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,SIZE(particle_set)
       atomic_kind => particle_set(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind,kind_number=j)
       atomkind(i) = j
    END DO
    DO i=1,SIZE(molecule_set)
       molecule => molecule_set(i)
       CALL get_molecule(molecule=molecule,&
            molecule_kind=molecule_kind,&
            first_atom=first_mol(i),&
            last_atom=last_mol(i))
    END DO
    ALLOCATE(atomlist(SIZE(particle_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(topology%exclude_vdw/=topology%exclude_ei) THEN
       ALLOCATE(atomlist_ei(SIZE(particle_set)),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    ALLOCATE(excl_molecule_list(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO imol_now = 1, SIZE(molecule_set)
       stat = 1
       CALL get_excluded_molecule(imol_now,excl_molecule_list,nbond,&
            stat,first_mol,last_mol,molecule_set,error)
       ALLOCATE(ex_bond_list(nbond),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ibond = 1
       DO i=1,SIZE(excl_molecule_list)
          IF (excl_molecule_list(i))   THEN
             molecule => molecule_set(i)
             CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind,&
                  bond_list=bond_list)
             ex_bond_list(ibond:ibond+SIZE(bond_list)-1)%a = &
                  bond_list(:)%a+first_mol(i)-1
             ex_bond_list(ibond:ibond+SIZE(bond_list)-1)%b = &
                  bond_list(:)%b+first_mol(i)-1
             ibond = ibond + SIZE(bond_list)
          END IF
       END DO
       stat = 2
       CALL get_excluded_molecule(imol_now,excl_molecule_list,nbend,&
            stat,first_mol,last_mol,molecule_set,error)
       ALLOCATE(ex_bend_list(nbend),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ibend = 1
       DO i=1,SIZE(excl_molecule_list)
          IF (excl_molecule_list(i))   THEN
             molecule => molecule_set(i)
             CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind,&
                  bend_list=bend_list)
             ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%a = &
                  bend_list(:)%a+first_mol(i)-1
             ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%b = &
                  bend_list(:)%b+first_mol(i)-1
             ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%c = &
                  bend_list(:)%c+first_mol(i)-1
             ibend = ibend + SIZE(bend_list)
          END IF
       END DO
       stat = 3
       CALL get_excluded_molecule(imol_now,excl_molecule_list,ntorsion,&
            stat,first_mol,last_mol,molecule_set,error)
       ALLOCATE(ex_torsion_list(ntorsion),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       itorsion = 1
       DO i=1,SIZE(excl_molecule_list)
          IF (excl_molecule_list(i))   THEN
             molecule => molecule_set(i)
             CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind,&
                  torsion_list=torsion_list)
             ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%a = &
                  torsion_list(:)%a+first_mol(i)-1
             ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%b = &
                  torsion_list(:)%b+first_mol(i)-1
             ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%c = &
                  torsion_list(:)%c+first_mol(i)-1
             ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%d = &
                  torsion_list(:)%d+first_mol(i)-1
             itorsion = itorsion + SIZE(torsion_list)
          END IF
       END DO
       DO iatom = first_mol(imol_now), last_mol(imol_now)
          atomlist        = .TRUE. 
          atomlist(iatom) = .FALSE.
       
          IF  (topology%exclude_vdw==do_skip_vdw_12 .OR. &
               topology%exclude_vdw==do_skip_vdw_13 .OR. &
               topology%exclude_vdw==do_skip_vdw_14) THEN
             DO j=1,nbond
                IF(  (ex_bond_list(j)%a==iatom).OR. &
                     (ex_bond_list(j)%b==iatom) ) THEN
                   atomlist(ex_bond_list(j)%a) = .FALSE.
                   atomlist(ex_bond_list(j)%b) = .FALSE.
                END IF
             END DO
          END IF
       
          IF  (topology%exclude_vdw==do_skip_vdw_13 .OR. &
               topology%exclude_vdw==do_skip_vdw_14) THEN
             DO j=1,nbend
                IF(  (ex_bend_list(j)%a==iatom).OR. &
                     (ex_bend_list(j)%c==iatom) ) THEN
                   atomlist(ex_bend_list(j)%a) = .FALSE.
                   atomlist(ex_bend_list(j)%c) = .FALSE.
                END IF
             END DO
          END IF
          
          IF  (topology%exclude_vdw==do_skip_vdw_14) THEN
             DO j=1,ntorsion
                IF(  (ex_torsion_list(j)%a==iatom).OR. &
                     (ex_torsion_list(j)%d==iatom) ) THEN
                   atomlist(ex_torsion_list(j)%a) = .FALSE.
                   atomlist(ex_torsion_list(j)%d) = .FALSE.
                END IF
             END DO
          END IF
          
          IF(topology%exclude_vdw/=topology%exclude_ei) THEN
             atomlist_ei        = .TRUE. 
             atomlist_ei(iatom) = .FALSE.
             IF  (topology%exclude_ei==do_skip_ei_12 .OR. &
                  topology%exclude_ei==do_skip_ei_13 .OR. &
                  topology%exclude_ei==do_skip_ei_14) THEN
                DO j=1,nbond
                   IF(  (ex_bond_list(j)%a==iatom).OR. &
                        (ex_bond_list(j)%b==iatom) ) THEN
                      atomlist_ei(ex_bond_list(j)%a) = .FALSE.
                      atomlist_ei(ex_bond_list(j)%b) = .FALSE.
                   END IF
                END DO
             END IF
             
             IF  (topology%exclude_ei==do_skip_ei_13 .OR. &
                  topology%exclude_ei==do_skip_ei_14) THEN
                DO j=1,nbend
                   IF(  (ex_bend_list(j)%a==iatom).OR. &
                        (ex_bend_list(j)%c==iatom) ) THEN
                      atomlist_ei(ex_bend_list(j)%a) = .FALSE.
                      atomlist_ei(ex_bend_list(j)%c) = .FALSE.
                   END IF
                END DO
             END IF
             
             IF  (topology%exclude_ei==do_skip_ei_14) THEN
                DO j=1,ntorsion
                   IF(  (ex_torsion_list(j)%a==iatom).OR. &
                        (ex_torsion_list(j)%d==iatom) ) THEN
                      atomlist_ei(ex_torsion_list(j)%a) = .FALSE.
                      atomlist_ei(ex_torsion_list(j)%d) = .FALSE.
                   END IF
                END DO
             END IF
          END IF
          
          IF (my_qmmm) THEN 
             DO j=1,SIZE(molecule_set)
                molecule => molecule_set(j)
                CALL get_molecule(molecule=molecule,&
                     first_atom=first_atom,last_atom=last_atom)
                IF ((iatom >= first_atom).AND.(iatom <= last_atom)) EXIT
             END DO
             imolecule=j
             !
             ! exclude all intramolecular interactions between QM atoms 
             !
             IF ((qmmm_env%qmmm_coupl_type == do_qmmm_none).AND.&
                  ANY(qmmm_env%qm_molecule_index == imolecule)) THEN
                DO j = 1, SIZE(qmmm_env%qm_atom_index)
                   counter = qmmm_env%qm_atom_index(j)
                   IF ((counter >= first_atom).AND.(counter <= last_atom)) THEN
                      atomlist(counter) = .FALSE.
                      IF(topology%exclude_vdw/=topology%exclude_ei) THEN
                         atomlist_ei(counter) = .FALSE.
                      END IF
                   END IF
                END DO
             END IF
          END IF

          DO j=1,SIZE(atomic_kind_set)
             counter = 0
             NULLIFY(list)
             DO k=1,SIZE(atomic_kind_set(j)%atom_list)
                kk = atomic_kind_set(j)%atom_list(k)
                IF(.NOT.atomlist(kk)) THEN
                   counter = counter + 1
                END IF
             END DO
             IF (counter>0) THEN
                ALLOCATE(list(counter),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter = 0
                DO k=1,SIZE(atomic_kind_set(j)%atom_list)
                   kk = atomic_kind_set(j)%atom_list(k)
                   IF(.NOT.atomlist(kk)) THEN
                      counter = counter + 1
                      list(counter) = kk
                   END IF
                END DO
                particle_set(iatom)%lpi%exclusion(j)%list_vdw => list
             END IF
             NULLIFY(list2)
             IF(topology%exclude_vdw==topology%exclude_ei) THEN
                IF(counter > 0) THEN
                   ALLOCATE(list2(counter),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   list2(:) = list(:)
                END IF
             ELSE
                counter = 0
                DO k=1,SIZE(atomic_kind_set(j)%atom_list)
                   kk = atomic_kind_set(j)%atom_list(k)
                   IF(.NOT.atomlist_ei(kk)) THEN
                      counter = counter + 1
                   END IF
                END DO
                IF (counter>0) THEN
                   ALLOCATE(list2(counter),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   counter = 0
                   DO k=1,SIZE(atomic_kind_set(j)%atom_list)
                      kk = atomic_kind_set(j)%atom_list(k)
                      IF(.NOT.atomlist_ei(kk)) THEN
                         counter = counter + 1
                         list2(counter) = kk
                      END IF
                   END DO
                END IF
             END IF
             IF (counter>0) THEN
                particle_set(iatom)%lpi%exclusion(j)%list_ei => list2
             END IF
             IF (iw>0) THEN 
                IF (ASSOCIATED(list)) &
                     WRITE(iw,*)"exclusion list_vdw :: ",&
                     "atom num :",iatom,"kind num. :",j,"exclusion list ::",&
                     list
                IF(topology%exclude_vdw/=topology%exclude_ei) THEN
                   IF (ASSOCIATED(list2)) &
                        WRITE(iw,*)"exclusion list_ei :: ",&
                        "atom num :",iatom,"kind num. :",j,"exclusion list ::",&
                        list2
                END IF
             END IF
          END DO
       END DO
       DEALLOCATE(ex_bond_list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ex_bend_list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ex_torsion_list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(atomlist,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(topology%exclude_vdw/=topology%exclude_ei) THEN
       DEALLOCATE(atomlist_ei,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(excl_molecule_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(atomkind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(first_mol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_mol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Set the atomic_kind_set()%elp_potentail%[qeff] (PART 1)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_10',handle2)
    IF(globenv%prog_name_id == do_fist) THEN
       DO i=1,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(i)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
          qeff = -10000.0_dp
          alpha_core_charge = -10000.0_dp
          qeff = charge(i)
          my_elp_name = TRIM("FIST_"//TRIM(atmname))
          my_elp_description(1) = "FIST Electrostatic Potential"
          my_elp_description(2) = "                            "
          NULLIFY(elp_potential)
          CALL allocate_potential(elp_potential)
          CALL set_potential(potential=elp_potential,&
               qeff=qeff,&
               alpha_core_charge=alpha_core_charge,&
               name=my_elp_name,&
               description=my_elp_description)
          CALL set_atomic_kind(atomic_kind=atomic_kind,elp_potential=elp_potential)
       END DO
    END IF
    DEALLOCATE(charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Set the atom_list for molecule_kind in molecule_kind_set (PART 2)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_11',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            natom=natom,molecule_list=molecule_list,&
            atom_list=atom_list)
       molecule => molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,&
            first_atom=first,last_atom=last)
       DO j=1,natom
          DO k=1,SIZE(atomic_kind_set)
             atomic_kind => atomic_kind_set(k)
             CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
             IF(globenv%prog_name_id == do_fist) THEN
                CALL get_atomic_kind(atomic_kind=atomic_kind,elp_potential=elp_potential)
                CALL get_potential(potential=elp_potential,qeff=qeff)
                IF( (atom_list(j)%name == atmname) .AND. &
                     (qeff==atom_info%atm_charge(first+j-1)) ) THEN
                   atom_list(j)%atomic_kind => atomic_kind_set(k)
                   EXIT
                END IF
             ELSE
                IF(atom_list(j)%name == atmname) THEN
                   atom_list(j)%atomic_kind => atomic_kind_set(k)
                   EXIT
                END IF
             END IF
          END DO
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,atom_list=atom_list)
    END DO
    CALL timestop(0.0_dp,handle2)

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_coordinate_pack
 
!******************************************************************************
!!****** topology_util/get_excluded_molecule [1.0] *
!!
!!   NAME
!!     get_excluded_molecule
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Gives back a list of molecules used to build up the exclusion lists
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!*** *************************************************************************
  SUBROUTINE get_excluded_molecule(imol,excl_molecule_list,ntarg,&
       stat,first_mol,last_mol,molecule_set,error)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN)            :: imol
    LOGICAL, POINTER, DIMENSION(:) :: excl_molecule_list
    INTEGER, POINTER, DIMENSION(:) :: first_mol, last_mol
    INTEGER, INTENT(OUT)           :: ntarg
    INTEGER, INTENT(INOUT)         :: stat
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'get_excluded_molecule', &
         routineP = moduleN//':'//routineN    
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list
    INTEGER :: i, j, my_mol, il
    INTEGER :: first_atom0, last_atom0
    INTEGER :: jatom, a, b, c, d, nloc, offset
    LOGICAL :: failure

    failure = .FALSE.
    excl_molecule_list = .FALSE.
    ntarg = 0
    excl_molecule_list(imol) = .TRUE.
    first_atom0 = first_mol(imol)
    last_atom0  = last_mol(imol)
    offset = first_atom0 - 1
    molecule => molecule_set(imol)
    CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
    SELECT CASE (stat)
    CASE(1)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            bond_list=bond_list, nbond=nloc)
       ntarg = ntarg + nloc
       DO j = 1, SIZE(bond_list)
          a = offset+bond_list(j)%a
          b = offset+bond_list(j)%b
          ! check boundaries for a
          IF ((a<first_atom0).OR.(a>last_atom0)) THEN
             jatom = a
             DO i = 1, SIZE(first_mol)
                IF ((first_mol(i) <= jatom).AND.(last_mol(i) >= jatom)) EXIT
             END DO
             CPPostcondition(i<=SIZE(first_mol),cp_failure_level,routineP,error,failure)
             IF (.NOT. excl_molecule_list(i)) THEN
                excl_molecule_list(i) = .TRUE.
                molecule => molecule_set(i)
                CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind,&
                     nbond=nloc)
                ntarg = ntarg + nloc
             END IF
          END IF
          ! check boundaries for b
          IF ((b<first_atom0).OR.(b>last_atom0)) THEN
             jatom = b
             DO i = 1, SIZE(first_mol)
                IF ((first_mol(i) <= jatom).AND.(last_mol(i) >= jatom)) EXIT
             END DO
             CPPostcondition(i<=SIZE(first_mol),cp_failure_level,routineP,error,failure)
             IF (.NOT. excl_molecule_list(i)) THEN
                excl_molecule_list(i) = .TRUE.
                molecule => molecule_set(i)
                CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind,&
                     nbond=nloc)
                ntarg = ntarg + nloc
             END IF
          END IF
       END DO
       ! Let's look for external bonds
       DO i = 1, SIZE(molecule_set)
          molecule => molecule_set(i)
          CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               bond_list=bond_list)
          il = last_mol(i)-first_mol(i)+1
          IF (ALL(bond_list(:)%a<=il).AND.ALL(bond_list(:)%b<=il)) CYCLE
          DO j = 1, SIZE(bond_list)
             a = first_mol(i)-1+bond_list(j)%a
             b = first_mol(i)-1+bond_list(j)%b
             IF ((a>=first_atom0.AND.a<=last_atom0).OR.&
                 (b>=first_atom0.AND.b<=last_atom0)) THEN
                IF (.NOT. excl_molecule_list(i)) THEN
                   excl_molecule_list(i) = .TRUE.
                   molecule => molecule_set(i)
                   CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                   CALL get_molecule_kind(molecule_kind=molecule_kind,&
                        nbond=nloc)
                   ntarg = ntarg + nloc
                   EXIT
                END IF
             END IF
          END DO
       END DO
    CASE(2)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            bend_list=bend_list, nbend=nloc)
       ntarg = ntarg + nloc
       DO j = 1, SIZE(bend_list)
          a = offset + bend_list(j)%a
          c = offset + bend_list(j)%c
          ! check boundaries for a
          IF ((a<first_atom0).OR.(a>last_atom0)) THEN
             jatom = a
             DO i = 1, SIZE(first_mol)
                IF ((first_mol(i) <= jatom).AND.(last_mol(i) >= jatom)) EXIT
             END DO
             CPPostcondition(i<=SIZE(first_mol),cp_failure_level,routineP,error,failure)
             IF (.NOT. excl_molecule_list(i)) THEN
                excl_molecule_list(i) = .TRUE.
                molecule => molecule_set(i)
                CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind,&
                     nbend=nloc)
                ntarg = ntarg + nloc
             END IF
          END IF
          ! check boundaries for c
          IF ((c<first_atom0).OR.(c>last_atom0)) THEN
             jatom = c
             DO i = 1, SIZE(first_mol)
                IF ((first_mol(i) <= jatom).AND.(last_mol(i) >= jatom)) EXIT
             END DO
             CPPostcondition(i<=SIZE(first_mol),cp_failure_level,routineP,error,failure)
             IF (.NOT. excl_molecule_list(i)) THEN
                excl_molecule_list(i) = .TRUE.
                molecule => molecule_set(i)
                CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind,&
                     nbend=nloc)
                ntarg = ntarg + nloc
             END IF
          END IF
       END DO
       ! Let's look for external bends
       DO i = 1, SIZE(molecule_set)
          molecule => molecule_set(i)
          CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               bend_list=bend_list)
          il = last_mol(i)-first_mol(i)+1
          IF (ALL(bend_list(:)%a<=il).AND.ALL(bend_list(:)%c<=il)) CYCLE
          DO j = 1, SIZE(bend_list)
             a = first_mol(i)-1+bend_list(j)%a
             c = first_mol(i)-1+bend_list(j)%c
             IF ((a>=first_atom0.AND.a<=last_atom0).OR.&
                 (c>=first_atom0.AND.c<=last_atom0)) THEN
                IF (.NOT. excl_molecule_list(i)) THEN
                   excl_molecule_list(i) = .TRUE.
                   molecule => molecule_set(i)
                   CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                   CALL get_molecule_kind(molecule_kind=molecule_kind,&
                        nbend=nloc)
                   ntarg = ntarg + nloc
                   EXIT
                END IF
             END IF
          END DO
       END DO
    CASE(3)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            torsion_list=torsion_list, ntorsion=nloc)
       ntarg = ntarg + nloc
       DO j = 1, SIZE(torsion_list)
          a = offset + torsion_list(j)%a
          d = offset + torsion_list(j)%d
          ! check boundaries for a
          IF ((a<first_atom0).OR.(a>last_atom0)) THEN
             jatom = a
             DO i = 1, SIZE(first_mol)
                IF ((first_mol(i) <= jatom).AND.(last_mol(i) >= jatom)) EXIT
             END DO
             CPPostcondition(i<=SIZE(first_mol),cp_failure_level,routineP,error,failure)
             IF (.NOT. excl_molecule_list(i)) THEN
                excl_molecule_list(i) = .TRUE.
                molecule => molecule_set(i)
                CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind,&
                     ntorsion=nloc)
                ntarg = ntarg + nloc
             END IF
          END IF
          ! check boundaries for d
          IF ((d<first_atom0).OR.(d>last_atom0)) THEN
             jatom = d
             DO i = 1, SIZE(first_mol)
                IF ((first_mol(i) <= jatom).AND.(last_mol(i) >= jatom)) EXIT
             END DO
             CPPostcondition(i<=SIZE(first_mol),cp_failure_level,routineP,error,failure)
             IF (.NOT. excl_molecule_list(i)) THEN
                excl_molecule_list(i) = .TRUE.
                molecule => molecule_set(i)
                CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind,&
                     ntorsion=nloc)
                ntarg = ntarg + nloc
             END IF
          END IF
       END DO
       ! Let's look for external torsions
       DO i = 1, SIZE(molecule_set)
          molecule => molecule_set(i)
          CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               torsion_list=torsion_list)
          il = last_mol(i)-first_mol(i)+1
          IF (ALL(torsion_list(:)%a<=il).AND.ALL(torsion_list(:)%d<=il)) CYCLE
          DO j = 1, SIZE(torsion_list)
             a = first_mol(i)-1+torsion_list(j)%a
             d = first_mol(i)-1+torsion_list(j)%d
             IF ((a>=first_atom0.AND.a<=last_atom0).OR.&
                 (d>=first_atom0.AND.d<=last_atom0)) THEN
                IF (.NOT. excl_molecule_list(i)) THEN
                   excl_molecule_list(i) = .TRUE.
                   molecule => molecule_set(i)
                   CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
                   CALL get_molecule_kind(molecule_kind=molecule_kind,&
                        ntorsion=nloc)
                   ntarg = ntarg + nloc
                   EXIT
                END IF
             END IF
          END DO
       END DO
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)       
    END SELECT

  END SUBROUTINE get_excluded_molecule
 
!******************************************************************************
!!****** topology_util/topology_constraint_pack [1.0] *
!!
!!   NAME
!!     topology_constraint_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in all the information needed for the constraints
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

  SUBROUTINE topology_constraint_pack ( molecule_kind_set,molecule_set,&
                                      topology,qmmm_env,particle_set,globenv,subsys_section,error)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_constraint_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: molname, name
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_typeh
    INTEGER :: dum, first, first_atom, handle, i, offset, iw, j, k, kk, last, &
      m, nbond, nhdist, nfixed_atoms, ng3x3, ng4x6, nmolecule, nrep, search1, &
      search2, stat, n_start_colv, ncolv_tot, ncolv_mol, kdim, ii, natom, &
      last_atom
    INTEGER, DIMENSION(:), POINTER           :: molecule_list,&
                                                colv_mol_num,&
                                                colv_mol_map
    LOGICAL                                  :: failure, fix_atom_molname, &
                                                fix_atom_qmmm
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: missed_molname
    REAL(KIND=dp)                            :: rmod, rvec(3), xdum
    REAL(KIND=dp), DIMENSION(:), POINTER     :: hdist
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(colvar_counters)                    :: ncolv
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    TYPE(local_colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: lcolv
    TYPE(local_constraint_type), POINTER     :: lci
    TYPE(local_g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: lg3x3
    TYPE(local_g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: lg4x6
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(section_vals_type), POINTER         :: hbonds_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    cons_info => topology%cons_info
    hbonds_section => section_vals_get_subs_vals(subsys_section,"CONSTRAINT%HBONDS",error=error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. NULLIFY the molecule_set(imol)%lci via set_molecule_set
    !-----------------------------------------------------------------------------
    DO i=1,topology%nmol
       molecule => molecule_set(i)
       NULLIFY(lci)
       ALLOCATE(lci,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(lci%lcolv)
       NULLIFY(lci%lg3x3)
       NULLIFY(lci%lg4x6)
       CALL set_molecule(molecule,lci=lci)
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Add more stuff to cons_info%const_dist_* if constraint hydrogen is on
    !-----------------------------------------------------------------------------
    IF(topology%const_hydr) THEN
       topology%const_colv = .TRUE.
       NULLIFY(atom_typeh, hdist)
       CALL section_vals_val_get(hbonds_section,"ATOM_TYPE",n_rep_val=nrep,error=error)
       IF (nrep/=0) &
            CALL section_vals_val_get(hbonds_section,"ATOM_TYPE",c_vals=atom_typeh,error=error)
       CALL section_vals_val_get(hbonds_section,"TARGETS",n_rep_val=nrep,error=error)
       IF (nrep/=0) &
            CALL section_vals_val_get(hbonds_section,"TARGETS",r_vals=hdist,error=error)
       IF (ASSOCIATED(hdist)) THEN
          CPPostcondition(SIZE(hdist)==SIZE(atom_typeh),cp_failure_level,routineP,error,failure)
       END IF
       nhdist = 0
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               bond_list=bond_list, nbond=nbond, atom_list=atom_list,&
               molecule_list=molecule_list)
          ! Let's tag all requested atoms involving Hydrogen 
          ! on the first molecule of this kind
          molecule => molecule_set(molecule_list(1))
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
          natom = last_atom - first_atom +1
          DO j =  1, natom
             atomic_kind => atom_list(j)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name,element_symbol=element_symbol) 
             IF((name(1:1)=="H").OR.(TRIM(element_symbol)=="H")) THEN
                DO k=1,nbond
                   IF ((bond_list(k)%a==j).OR.(bond_list(k)%b==j))  THEN
                      nhdist = nhdist + 1
                   END IF
                END DO
             END IF
          END DO
       END DO
       n_start_colv = cons_info%nconst_colv
       cons_info%nconst_colv = nhdist + n_start_colv
       CALL reallocate(cons_info%const_colv_mol,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%const_colv_target,1,cons_info%nconst_colv)
       CALL colvar_p_reallocate(cons_info%colvar_set,1,cons_info%nconst_colv)
       !
       nhdist = 0
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               bond_list=bond_list, nbond=nbond, atom_list=atom_list,&
               molecule_list=molecule_list)
          molecule => molecule_set(molecule_list(1))
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
          natom = last_atom  - first_atom + 1
          offset = first_atom - 1
          DO j =  1, natom
             atomic_kind => atom_list(j)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name,element_symbol=element_symbol) 
             IF((name(1:1)=="H").OR.(TRIM(element_symbol)=="H")) THEN
                DO k=1,nbond
                   IF ((bond_list(k)%a==j).OR.(bond_list(k)%b==j))  THEN
                      rvec = particle_set(offset+bond_list(k)%a)%r-particle_set(offset+bond_list(k)%b)%r
                      rmod = SQRT(DOT_PRODUCT(rvec,rvec))
                      IF (ASSOCIATED(hdist)) THEN
                         IF (SIZE(hdist)>0) THEN
                            IF (bond_list(k)%a==j) atomic_kind => atom_list(bond_list(k)%b)%atomic_kind
                            IF (bond_list(k)%b==j) atomic_kind => atom_list(bond_list(k)%a)%atomic_kind
                            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                 name=name,element_symbol=element_symbol) 
                            DO m = 1, SIZE(hdist)
                               IF (TRIM(name)==TRIM(atom_typeh(m))) EXIT
                               IF (TRIM(element_symbol)==TRIM(atom_typeh(m))) EXIT
                            END DO
                            IF (m<=SIZE(hdist)) THEN
                               rmod = hdist(m)
                            END IF
                         END IF
                      END IF
                      nhdist = nhdist + 1
                      cons_info%const_colv_mol(nhdist+n_start_colv)     = i
                      cons_info%const_colv_target(nhdist+n_start_colv)  = rmod
                      CALL colvar_create(cons_info%colvar_set(nhdist+n_start_colv)%colvar,&
                           dist_colvar_id, error)
                      cons_info%colvar_set(nhdist+n_start_colv)%colvar%dist_param%i_at = bond_list(k)%a  
                      cons_info%colvar_set(nhdist+n_start_colv)%colvar%dist_param%j_at = bond_list(k)%b
                      CALL colvar_setup(cons_info%colvar_set(nhdist+n_start_colv)%colvar, error)
                   END IF
                END DO
             END IF
          END DO
       END DO      
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Set the distance constraint molecule_kind_set(ikind)%colv_list
    !-----------------------------------------------------------------------------
    IF(topology%const_colv) THEN
       DO i = 1, SIZE(cons_info%colvar_set)
          ii = cons_info%const_colv_mol(i)
          IF (ii==0) CYCLE
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ncolv_mol = COUNT(cons_info%const_colv_mol==ii)
          ALLOCATE(colv_list(ncolv_mol),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ncolv_mol = 0
          DO j = 1, SIZE(cons_info%colvar_set)
             IF (cons_info%const_colv_mol(j)==ii) THEN
                ncolv_mol = ncolv_mol + 1
                kdim = SIZE(cons_info%colvar_set(j)%colvar%i_atom)
                ALLOCATE(colv_list(ncolv_mol)%i_atoms(kdim),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                colv_list(ncolv_mol)%type_id        = cons_info%colvar_set(j)%colvar%type_id
                colv_list(ncolv_mol)%i_atoms        = cons_info%colvar_set(j)%colvar%i_atom
                IF (cons_info%const_colv_target(j)==-HUGE(0.0_dp)) THEN
                   CALL colvar_eval_mol_f(cons_info%colvar_set(j)%colvar, topology%cell,&
                        particle_set, error=error)
                   colv_list(ncolv_mol)%expected_value = cons_info%colvar_set(j)%colvar%ss
                ELSE
                   colv_list(ncolv_mol)%expected_value = cons_info%const_colv_target(j)
                END IF
             END IF
          END DO
          CALL setup_colvar_counters(colv_list,ncolv)
          CALL set_molecule_kind(molecule_kind,colv_list=colv_list,ncolv=ncolv)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lcolv(ncolv_mol),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ncolv_mol = 0
             DO k=1,SIZE(cons_info%colvar_set)
                IF (cons_info%const_colv_mol(k)==ii) THEN
                   ncolv_mol = ncolv_mol + 1
                   lcolv(ncolv_mol)%lambda    = 0.0_dp
                   lcolv(ncolv_mol)%sigma     = 0.0_dp
                   ! Set Up colvar variable
                   NULLIFY(lcolv(ncolv_mol)%colvar, lcolv(ncolv_mol)%colvar_old)
                   ! Colvar
                   CALL colvar_clone(lcolv(ncolv_mol)%colvar, cons_info%colvar_set(k)%colvar,&
                        i_atom_offset=first_atom-1,error=error)
                   ! Colvar_old
                   CALL colvar_clone(lcolv(ncolv_mol)%colvar_old,lcolv(ncolv_mol)%colvar,error=error)
                   ! Check for consistency in the constraint definition
                   IF  (ANY(lcolv(ncolv_mol)%colvar%i_atom > last_atom).OR.&
                        ANY(lcolv(ncolv_mol)%colvar%i_atom < first_atom)) THEN
                      WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
                      WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
                           " but the atoms specified in the constraint and the atoms defined for",&
                           " the molecule DO NOT match!",&
                           "This could be very probable due to a wrong connectivity, or an error",&
                           " in the constraint specification in the input file.",&
                           " Please check it carefully!"
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)                   
                   END IF
                END IF
             END DO
             CALL set_molecule(molecule=molecule,lcolv=lcolv)
          END DO
          WHERE (cons_info%const_colv_mol==ii) cons_info%const_colv_mol=0
       END DO
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the group 3x3 constraint g3x3_list
    !-----------------------------------------------------------------------------
    IF(topology%const_33) THEN
       DO i=1,cons_info%nconst_g33
          ii = cons_info%const_g33_mol(i)
          IF (ii==0) CYCLE
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ng3x3 = COUNT(cons_info%const_g33_mol==ii)
          ALLOCATE(g3x3_list(ng3x3),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ng3x3 = 0
          DO j=1,SIZE(cons_info%const_g33_mol)
             IF (cons_info%const_g33_mol(j)==ii) THEN
                ng3x3 = ng3x3 + 1
                g3x3_list(ng3x3)%a   = cons_info%const_g33_a(j)
                g3x3_list(ng3x3)%b   = cons_info%const_g33_b(j)
                g3x3_list(ng3x3)%c   = cons_info%const_g33_c(j)
                g3x3_list(ng3x3)%dab = cons_info%const_g33_dab(j)
                g3x3_list(ng3x3)%dac = cons_info%const_g33_dac(j)
                g3x3_list(ng3x3)%dbc = cons_info%const_g33_dbc(j)
             END IF
          END DO     
          CALL set_molecule_kind(molecule_kind,ng3x3=ng3x3,g3x3_list=g3x3_list)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lg3x3(ng3x3),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO k=1,ng3x3
                lg3x3(k)%scale     = 0.0_dp
                lg3x3(k)%scale_old = 0.0_dp
                lg3x3(k)%fa        = 0.0_dp
                lg3x3(k)%fb        = 0.0_dp
                lg3x3(k)%fc        = 0.0_dp
                lg3x3(k)%ra_old    = 0.0_dp
                lg3x3(k)%rb_old    = 0.0_dp
                lg3x3(k)%rc_old    = 0.0_dp
                lg3x3(k)%va        = 0.0_dp
                lg3x3(k)%vb        = 0.0_dp
                lg3x3(k)%vc        = 0.0_dp
                lg3x3(k)%lambda    = 0.0_dp
                IF  ((g3x3_list(k)%a+first_atom-1 < first_atom).OR.&
                     (g3x3_list(k)%b+first_atom-1 < first_atom).OR.& 
                     (g3x3_list(k)%c+first_atom-1 < first_atom).OR.& 
                     (g3x3_list(k)%a+first_atom-1 > last_atom ).OR.&
                     (g3x3_list(k)%b+first_atom-1 > last_atom ).OR.& 
                     (g3x3_list(k)%c+first_atom-1 > last_atom )) THEN
                   WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
                   WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
                        " but the atoms specified in the constraint and the atoms defined for",&
                        " the molecule DO NOT match!",&
                        "This could be very probable due to a wrong connectivity, or an error",&
                        " in the constraint specification in the input file.",&
                        " Please check it carefully!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)                   
                END IF
             END DO
             CALL set_molecule(molecule=molecule,lg3x3=lg3x3)
          END DO
          WHERE (cons_info%const_g33_mol==ii) cons_info%const_g33_mol=0
       END DO
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Set the group 4x6 constraint g4x6_list
    !-----------------------------------------------------------------------------
    IF(topology%const_46) THEN
       DO i=1,cons_info%nconst_g46
          ii = cons_info%const_g46_mol(i)
          IF (ii==0) CYCLE
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ng4x6 = COUNT(cons_info%const_g46_mol==ii)
          ALLOCATE(g4x6_list(ng4x6),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ng4x6 = 0
          DO j=1,SIZE(cons_info%const_g46_a)
             IF (cons_info%const_g46_mol(j)==ii) THEN
                ng4x6 = ng4x6 + 1
                g4x6_list(ng4x6)%a   = cons_info%const_g46_a(j)
                g4x6_list(ng4x6)%b   = cons_info%const_g46_b(j)
                g4x6_list(ng4x6)%c   = cons_info%const_g46_c(j)
                g4x6_list(ng4x6)%d   = cons_info%const_g46_d(j)
                g4x6_list(ng4x6)%dab = cons_info%const_g46_dab(j)
                g4x6_list(ng4x6)%dac = cons_info%const_g46_dac(j)
                g4x6_list(ng4x6)%dbc = cons_info%const_g46_dbc(j)
                g4x6_list(ng4x6)%dad = cons_info%const_g46_dad(j)
                g4x6_list(ng4x6)%dbd = cons_info%const_g46_dbd(j)
                g4x6_list(ng4x6)%dcd = cons_info%const_g46_dcd(j)
             END IF
          END DO
          CALL set_molecule_kind(molecule_kind,ng4x6=ng4x6,g4x6_list=g4x6_list)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lg4x6(ng4x6),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO k=1,ng4x6
                lg4x6(k)%scale     = 0.0_dp
                lg4x6(k)%scale_old = 0.0_dp
                lg4x6(k)%fa        = 0.0_dp
                lg4x6(k)%fb        = 0.0_dp
                lg4x6(k)%fc        = 0.0_dp
                lg4x6(k)%fd        = 0.0_dp
                lg4x6(k)%fe        = 0.0_dp
                lg4x6(k)%ff        = 0.0_dp
                lg4x6(k)%ra_old    = 0.0_dp
                lg4x6(k)%rb_old    = 0.0_dp
                lg4x6(k)%rc_old    = 0.0_dp
                lg4x6(k)%rd_old    = 0.0_dp
                lg4x6(k)%re_old    = 0.0_dp
                lg4x6(k)%rf_old    = 0.0_dp
                lg4x6(k)%va        = 0.0_dp
                lg4x6(k)%vb        = 0.0_dp
                lg4x6(k)%vc        = 0.0_dp
                lg4x6(k)%vd        = 0.0_dp
                lg4x6(k)%ve        = 0.0_dp
                lg4x6(k)%vf        = 0.0_dp
                lg4x6(k)%lambda    = 0.0_dp
                IF  ((g4x6_list(k)%a+first_atom-1 < first_atom).OR.&
                     (g4x6_list(k)%b+first_atom-1 < first_atom).OR.& 
                     (g4x6_list(k)%c+first_atom-1 < first_atom).OR.& 
                     (g4x6_list(k)%d+first_atom-1 < first_atom).OR.& 
                     (g4x6_list(k)%a+first_atom-1 > last_atom ).OR.&
                     (g4x6_list(k)%b+first_atom-1 > last_atom ).OR.& 
                     (g4x6_list(k)%c+first_atom-1 > last_atom ).OR.& 
                     (g4x6_list(k)%d+first_atom-1 > last_atom )) THEN
                   WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
                   WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
                        " but the atoms specified in the constraint and the atoms defined for",&
                        " the molecule DO NOT match!",&
                        "This could be very probable due to a wrong connectivity, or an error",&
                        " in the constraint specification in the input file.",&
                        " Please check it carefully!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)                   
                END IF
             END DO
             CALL set_molecule(molecule=molecule,lg4x6=lg4x6)
          END DO
          WHERE (cons_info%const_g46_mol==ii) cons_info%const_g46_mol=0
       END DO
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the group fixed_atom constraint fixd_list
    !-----------------------------------------------------------------------------
    IF(topology%const_atom) THEN
       ALLOCATE(missed_molname(SIZE(cons_info%fixed_molnames,1)),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       missed_molname=.TRUE.
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list, name=molname)
          WHERE(molname.EQ.cons_info%fixed_molnames)
            missed_molname=.FALSE.
          END WHERE
          ! Try to figure out how many atoms of the list belong to this molecule_kind
          nfixed_atoms = 0
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first, last_atom = last )
             fix_atom_molname=.FALSE.
             IF (ASSOCIATED(cons_info%fixed_molnames)) THEN
                IF (ANY(cons_info%fixed_molnames.EQ.molname)) fix_atom_molname=.TRUE.
             ENDIF
             fix_atom_qmmm=.FALSE.
             IF (PRESENT(qmmm_env)) THEN
                IF  ((cons_info%freeze_qm.AND.ANY(qmmm_env%qm_molecule_index == molecule_list(j))).OR.&
                     (cons_info%freeze_mm.AND.ALL(qmmm_env%qm_molecule_index /= molecule_list(j)))) THEN
                  fix_atom_qmmm=.TRUE.
                ENDIF
             ENDIF
             DO k = first, last
                IF  (ANY(cons_info%fixed_atoms == k) .OR. fix_atom_qmmm .OR. fix_atom_molname) THEN
                        nfixed_atoms = nfixed_atoms + 1 
                END IF
             END DO
          END DO
          ALLOCATE(fixd_list(nfixed_atoms),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          kk = 0
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first, last_atom = last )
             fix_atom_molname=.FALSE.
             IF (ASSOCIATED(cons_info%fixed_molnames)) THEN
                IF (ANY(cons_info%fixed_molnames.EQ.molname)) fix_atom_molname=.TRUE.
             ENDIF
             fix_atom_qmmm=.FALSE.
             IF (PRESENT(qmmm_env)) THEN
                IF  ((cons_info%freeze_qm.AND.ANY(qmmm_env%qm_molecule_index == molecule_list(j))).OR.&
                     (cons_info%freeze_mm.AND.ALL(qmmm_env%qm_molecule_index /= molecule_list(j)))) THEN
                  fix_atom_qmmm=.TRUE.
                ENDIF
             ENDIF
             DO k = first, last
                IF  (ANY(cons_info%fixed_atoms == k) .OR. fix_atom_qmmm .OR. fix_atom_molname) THEN                     
                      kk = kk + 1
                      fixd_list(kk)%fixd = k
                      particle_set(k)%is_fixed = .TRUE.
                END IF
             END DO
          END DO
          IF (iw>0) WRITE(iw,*)"FIXD_LIST ::",fixd_list
          CALL set_molecule_kind(molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list)
       END DO
       CPPostcondition(COUNT(missed_molname)==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(missed_molname,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE topology_constraint_pack
!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_molecules_check
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Check and verify that all molecules of the same kind are bonded the same
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_molecules_check(topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_molecules_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, first, first2, &
                                                handle, iatom, imol, &
                                                imol_typ, iw, last, last2, &
                                                natom, nmol, nmol_typ
    LOGICAL                                  :: failure
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1.
    !-----------------------------------------------------------------------------
    natom = topology%natoms
    nmol_typ = atom_info%map_mol_typ(natom)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Check each molecule have the same number of atoms
    !-----------------------------------------------------------------------------
    DO imol_typ=1,nmol_typ
       CALL find_boundary(atom_info%map_mol_typ,natom,first,last,imol_typ)
       !get back the first and last atom for this molecule type
       nmol = atom_info%map_mol_num(last)
       IF(nmol > 1) THEN
          CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,&
               natom,first,last,imol_typ,1)
          counter = last - first + 1
          DO imol=2,nmol
             CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,&
                  natom,first,last,imol_typ,imol)
             !get back the first and last atom of this molecule
             iatom = last - first + 1
             IF(iatom /= counter) THEN
                CALL stop_program ("topology_molecules_check",&
                     "different number of atoms for same molecule kind"//&
                     " imol_typ= "//cp_to_string(imol_typ)//&
                     " imol= "//cp_to_string(imol))
             END IF
          END DO
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Check that each molecule has the same atom sequences
    !-----------------------------------------------------------------------------
    DO imol_typ=1,nmol_typ
       CALL find_boundary(atom_info%map_mol_typ,natom,first,last,imol_typ)
       !get back the first and last atom for this molecule type
       nmol = atom_info%map_mol_num(last)
       IF(nmol > 1) THEN
          CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,&
               natom,first,last,imol_typ,1)
          counter = last - first + 1
          DO imol=2,nmol
             CALL find_boundary(atom_info%map_mol_typ,atom_info%map_mol_num,&
                  natom,first2,last2,imol_typ,imol)
             DO iatom=1,counter
                IF( atom_info%label_atmname(first+iatom-1) /= &
                     atom_info%label_atmname(first2+iatom-1) ) THEN
                   CALL stop_program ("topology_molecules_check",&
                        "different atom name for same molecule kind"&
                        //" imol_typ= "//TRIM(cp_to_string(imol_typ))&
                        //" imol= "//TRIM(cp_to_string(imol))&
                        //" "//TRIM(atom_info%label_atmname(first+iatom-1))&
                       //" "//TRIM(atom_info%label_atmname(first2+iatom-1)))
                END IF
             END DO
          END DO
       END IF
    END DO


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Check that each molecule have the same bond sequences
    !-----------------------------------------------------------------------------
    ! Need to put in later, this is not necesary for simple molecules as the two 
    ! above check will catch all the problems. This will be needed in the 
    ! instance where you have a bunch of ring per se and one of the bond is 
    ! streatched.


    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)

  END SUBROUTINE topology_molecules_check


!!****f* cp_subsystem_methods/check_subsys_element [1.0] *
!!
!!   NAME
!!     check_subsys_element
!!
!!   FUNCTION
!!     Check and returns the ELEMENT label
!!
!!   NOTES
!!
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE check_subsys_element(element_in, element_out, globenv, force_env_section, error)
    CHARACTER(len=*), INTENT(IN)             :: element_in
    CHARACTER(len=default_string_length), &
      INTENT(OUT)                            :: element_out
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'check_subsys_element', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: current_symbol, &
                                                element_symbol, keyword
    INTEGER                                  :: handle, i, i_rep, n_rep
    INTEGER, SAVE                            :: iwarn = 20
    LOGICAL                                  :: failure, found
    TYPE(section_vals_type), POINTER         :: kind_section

    found = .FALSE.
    failure = .FALSE.
    element_symbol = element_in
    element_out = ""
    CALL uppercase(element_symbol)
    DO i=LBOUND(ptable,1),UBOUND(ptable,1)
       current_symbol = ptable(i)%symbol
       CALL uppercase(current_symbol)
       IF (TRIM(current_symbol) == TRIM(element_symbol)) THEN
          element_out = TRIM(element_symbol)
          found = .TRUE.
          EXIT
       END IF
    END DO
    IF (.NOT.found) THEN
       kind_section => section_vals_get_subs_vals(force_env_section,&
            "KIND",error=error)
       CALL section_vals_get(kind_section,n_repetition=n_rep,error=error)
       DO i_rep=1,n_rep
          CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
               c_val=keyword,i_rep_section=i_rep,error=error)
          CALL uppercase(keyword)
          IF (TRIM(keyword)==TRIM(element_symbol)) THEN
             CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
                  keyword_name="ELEMENT",n_rep_val=i,error=error)
             IF (i>0) THEN
                CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
                     keyword_name="ELEMENT",c_val=element_symbol,error=error)
                element_out = TRIM(element_symbol)
                found = .TRUE.                
                EXIT
             END IF
          END IF
       END DO
    END IF
    IF (.NOT.found) THEN
       IF (iwarn > 0) THEN
          CALL Pwarning("No element specification provided for type:"//element_in//".",globenv, error)
          iwarn = iwarn - 1
          IF (iwarn == 0) THEN
             CALL Pwarning("Further Output regarding unknown element types will be suppressed!",&
                  globenv,error)
          END IF
       END IF
    END IF

  END SUBROUTINE check_subsys_element

END MODULE topology_util
