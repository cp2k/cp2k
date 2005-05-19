!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             allocate_potential,&
                                             set_potential
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             bend_type,&
                                             bond_type,&
                                             dist_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             fixd_constraint_type,&
                                             get_molecule_kind,&
                                             impr_type,&
                                             molecule_kind_type,&
                                             onfo_type,&
                                             set_molecule_kind,&
                                             torsion_type,&
                                             ub_type
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             get_molecule,&
                                             local_constraint_type,&
                                             local_dist_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type,&
                                             local_molecule_type,&
                                             molecule_type,&
                                             set_molecule,&
                                             set_molecule_set
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: massunit
  USE simulation_cell,                 ONLY: pbc
  USE string_utilities,                ONLY: integer_to_string,&
                                             uppercase,&
                                             str_search
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             atom_info_type,&
                                             connectivity_info_type,&
                                             constraint_info_type
  USE qmmm_ff_fist, ONLY: qmmm_ff_precond_only_qm
  USE qmmm_types,   ONLY: qmmm_env_mm_type
  USE input_constants,   ONLY: do_qmmm_none
  

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_util"


  PRIVATE
  PUBLIC :: find_boundary,&
            find_boundary2,&
            find_boundary3,&
            topology_set_atm_mass,&
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
            topology_generate_ub

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
SUBROUTINE topology_reorder_atoms(topology,globenv)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_reorder_atoms', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: PM1, PM2, PM3, PM4

  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

    CHARACTER (LEN=default_string_length), POINTER   :: tlabel_resname(:)
    CHARACTER (LEN=default_string_length), POINTER   :: tlabel_atmname(:)
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: tr
    INTEGER, POINTER :: tmap_mol_typ(:)
    INTEGER, POINTER :: tmap_mol_num(:)
    REAL(KIND=dp), POINTER :: tatm_charge(:)
    REAL(KIND=dp), POINTER :: tatm_mass(:)
    CHARACTER (LEN=default_string_length), POINTER   :: telement(:)

    INTEGER      :: istat,iw
    INTEGER      :: location,iatm_count,imol_count
    INTEGER      :: ityp,inum,iatm,ntyp,nnum,natom
    INTEGER, POINTER :: new_position(:)
    INTEGER, POINTER :: mnum(:)
    INTEGER, POINTER :: mtyp(:)

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1)  WRITE(iw,*) "    Entering topology_reorder_atoms"

  atom_info => topology%atom_info
  conn_info => topology%conn_info

  natom = topology%natoms

  NULLIFY(mtyp,mnum,new_position)
  NULLIFY(tlabel_resname,tlabel_atmname,telement)
  NULLIFY(tr,tatm_charge,tatm_mass)
  NULLIFY(tmap_mol_typ,tmap_mol_num)

  !ALLOCATE all the temporary arrays needed for this routine
  ALLOCATE(mtyp(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mtyp')
  ALLOCATE(mnum(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mnum')
  ALLOCATE(new_position(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','new_position')
  ALLOCATE(tlabel_resname(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_resname')
  ALLOCATE(tlabel_atmname(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_atmname')
  ALLOCATE(tr(3,natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tr')
  ALLOCATE(tatm_charge(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tatm_charge')
  ALLOCATE(tatm_mass(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tatm_mass')
  ALLOCATE(tmap_mol_typ(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tmap_mol_typ')
  ALLOCATE(tmap_mol_num(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tmap_mol_num')
  ALLOCATE(telement(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','telement')

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
  DEALLOCATE(mtyp,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mtyp')
  DEALLOCATE(mnum,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mnum')
  DEALLOCATE(new_position,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','new_position')
  DEALLOCATE(tlabel_resname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_resname')
  DEALLOCATE(tlabel_atmname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_atmname')
  DEALLOCATE(telement,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','telement')
  DEALLOCATE(tr,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tr')
  DEALLOCATE(tatm_charge,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tatm_charge')

  !DEALLOCATE all the stuff in topology structure so it can be re-runned
  DEALLOCATE(atom_info%map_mol_typ,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','map_mol_typ')
  DEALLOCATE(atom_info%map_mol_num,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','map_mol_num')
  DEALLOCATE(conn_info%bond_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','bond_a')
  DEALLOCATE(conn_info%bond_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','bond_b')
  DEALLOCATE(conn_info%theta_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','theta_a')
  DEALLOCATE(conn_info%theta_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','theta_b')
  DEALLOCATE(conn_info%theta_c,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','theta_c')
  DEALLOCATE(conn_info%phi_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_a')
  DEALLOCATE(conn_info%phi_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_b')
  DEALLOCATE(conn_info%phi_c,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_c')
  DEALLOCATE(conn_info%phi_d,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_d')


  IF(PM1) WRITE(iw,*) "    Exiting  topology_reorder_atoms"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF


!------------------------------------------------------------------------------
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
SUBROUTINE topology_generate_molecule(topology,qmmm,qmmm_env,globenv,error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv
    TYPE(qmmm_env_mm_type), POINTER, OPTIONAL :: qmmm_env
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error    

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_molecule', &
      routineP = moduleN//':'//routineN

  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

  INTEGER                                  :: first, iatom, iiatom, iii, &
                                              inum, istat, itype, iw, last, &
                                              natom, nbond, nnum, ntype
  INTEGER                                  :: map_mol_typ_orig, jump
  INTEGER                                  :: ifirst, ilast, myind, iatm, stat, jatm
  INTEGER, POINTER, DIMENSION(:)           :: checkme, qm_atom_index
  LOGICAL                                  :: PM1, PM2, PM3, PM4, my_qmmm, do_again
  LOGICAL                                  :: failure
  INTEGER                                  :: jtype,nmol

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  failure = .FALSE.
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
  NULLIFY(checkme, qm_atom_index)
  IF(PM1)  WRITE(iw,*) "    Entering topology_generate_molecule"

!------------------------------------------------------------------------------

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
 
  ALLOCATE(checkme(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('generate_molecule','checkme')

  IF (ASSOCIATED(atom_info%map_mol_typ)) &
    DEALLOCATE(atom_info%map_mol_typ,STAT=istat)
  IF (istat /= 0) CALL stop_memory('generate_molecule','dealloc :: map_mol_typ')
  ALLOCATE (atom_info%map_mol_typ(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('generate_molecule','alloc :: map_mol_typ')

  IF (ASSOCIATED(atom_info%map_mol_num)) &
    DEALLOCATE(atom_info%map_mol_num,STAT=istat)
  IF (istat /= 0) CALL stop_memory('generate_molecule','dealloc :: map_mol_num')
  ALLOCATE (atom_info%map_mol_num(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('generate_molecule','alloc :: map_mol_num')

  IF (ASSOCIATED(atom_info%map_mol_res)) &
    DEALLOCATE(atom_info%map_mol_res,STAT=istat)
  IF (istat /= 0) CALL stop_memory('generate_molecule','dealloc :: map_mol_res')
  ALLOCATE (atom_info%map_mol_res(natom),STAT=istat)

  !Zero the arrays
  checkme(:) = 0
  atom_info%map_mol_typ(:) = 0
  atom_info%map_mol_num(:) = 0
  atom_info%map_mol_res(:) = 1

  !Parse the atoms list to see how many different molecule types there are
  ntype = 1
  atom_info%map_mol_typ(1)=1
  DO iatom=2,natom
    IF(atom_info%label_molname(iatom-1)==atom_info%label_molname(iatom)) THEN
      atom_info%map_mol_typ(iatom)=ntype
    ELSE
      ntype=ntype+1
      atom_info%map_mol_typ(iatom)=ntype
    END IF
  END DO

  DO itype=1,ntype
    !search for first atom in the type list and last to cut down time
    CALL find_boundary(atom_info%map_mol_typ,natom,first,last,itype)
    IF(PM2) WRITE(iw,*) "      itype ",itype,first,last

    inum = 0

    DO iatom=first,last
      checkme(:) = 0

      IF(atom_info%map_mol_num(iatom)==0) THEN
        inum=inum+1
        checkme(iatom)=1
        DO WHILE (SUM(checkme)>0)
          DO iiatom=first,last
            IF((checkme(iiatom)==1)) THEN
              checkme(iiatom)=0
              atom_info%map_mol_num(iiatom)=inum
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
  END DO

  DEALLOCATE(checkme,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','checkme')

  IF(topology%para_res) THEN
    DO itype=1,ntype
      CALL find_boundary(atom_info%map_mol_typ,natom,first,last,itype)
      nmol = atom_info%map_mol_num(last)
      DO jtype = 1,nmol
        CALL find_boundary2(atom_info%map_mol_typ,atom_info%map_mol_num,&
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

  IF (my_qmmm) THEN
     do_again = .FALSE.
     IF (PM4) WRITE(iw,*)"MAP_MOL_NUM ",atom_info%map_mol_num
     IF (PM4) WRITE(iw,*)"MAP_MOL_TYP ",atom_info%map_mol_typ
     ALLOCATE(qm_atom_index(SIZE(qmmm_env%qm_atom_index)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     qm_atom_index = qmmm_env%qm_atom_index
     CPPostcondition(ALL(qm_atom_index /= 0),cp_failure_level,routineP,error,failure)
     DO myind = 1, SIZE(qm_atom_index)
        ifirst = qm_atom_index(myind)
        IF ( ifirst == 0 ) CYCLE
        ilast  = ifirst
        DO WHILE ((atom_info%map_mol_num(ifirst) == atom_info%map_mol_num(ilast)).AND.&
                  (atom_info%map_mol_typ(ifirst) == atom_info%map_mol_typ(ilast)).AND.&
                  ifirst > 1)
           ifirst = ifirst - 1
        END DO
        CPPostcondition(ifirst>=1,cp_failure_level,routineP,error,failure)
        IF ( ifirst > 1 ) ifirst = ifirst + 1
        DO WHILE ((atom_info%map_mol_num(ifirst) == atom_info%map_mol_num(ilast)).AND.&
                  (atom_info%map_mol_typ(ifirst) == atom_info%map_mol_typ(ilast)).AND.&
                  ilast < natom )
           ilast = ilast + 1
        END DO
        CPPostcondition(ilast<=natom,cp_failure_level,routineP,error,failure)
        IF ( ilast < natom ) ilast = ilast - 1
        IF (PM2) WRITE(iw,*)"qm fragment:: ifirst, ilast",ifirst,ilast
        DO iatm = ifirst, ilast
           atom_info%label_molname(iatm) = "_QM_"//&
                TRIM(atom_info%label_molname(iatm))
           IF (PM2) WRITE(iw,*)"QM Molecule name :: ", atom_info%label_molname(iatm)
           WHERE (qm_atom_index == iatm) qm_atom_index = 0
        END DO
        DO iatm = 1, ifirst-1
           IF ( ANY(qm_atom_index == iatm) ) do_again = .TRUE.
        END DO
        DO iatm = ilast+1, natom
           IF ( ANY(qm_atom_index == iatm) ) do_again = .TRUE.
        END DO
        IF (PM2) WRITE(iw,*)" Another QM fragment? :: ",do_again
        DO iatm = ifirst, ilast
           atom_info%map_mol_num(iatm) = 1           
        END DO
        map_mol_typ_orig = atom_info%map_mol_typ(ifirst)
        IF (ifirst /= 1) THEN
           IF (atom_info%map_mol_typ(ifirst-1) == map_mol_typ_orig) THEN
              DO iatm = ifirst, ilast
                 atom_info%map_mol_typ(iatm) = map_mol_typ_orig + 1
              END DO
           END IF
        END IF
        IF (ilast /= natom) THEN
           inum = atom_info%map_mol_num(ilast+1)
           DO iatm = ilast+1, natom
              IF (atom_info%map_mol_num(iatm) == 1) EXIT
              atom_info%map_mol_num(iatm) = atom_info%map_mol_num(iatm) - inum +1
           END DO
           map_mol_typ_orig = atom_info%map_mol_typ(ifirst)
           IF (atom_info%map_mol_typ(ilast+1) < map_mol_typ_orig+1) THEN
              jump = 2
              IF (ifirst == 1) jump = 1
              DO iatm = ilast+1, natom
                 atom_info%map_mol_typ(iatm) =  atom_info%map_mol_typ(iatm) + jump
              END DO
           ELSE IF (atom_info%map_mol_typ(ilast+1) > map_mol_typ_orig+1) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
           END IF
        END IF
        IF (.NOT.do_again) EXIT
     END DO
     DEALLOCATE(qm_atom_index,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF

  IF(PM4) THEN
    DO iatom=1,natom
      WRITE(iw,*) "      iatom,map_mol_typ,map_mol_num ",iatom,&
                 atom_info%map_mol_typ(iatom),atom_info%map_mol_num(iatom)
    END DO
  END IF

!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_molecule"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_molecule

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary(num_array,ntot,first,last,search)


    INTEGER, POINTER                         :: num_array(:)
    INTEGER                                  :: ntot, first, last, search

    INTEGER                                  :: i
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

  found=.FALSE.
  first=0
  last=ntot

  DO i=1,ntot
    IF(num_array(i)==search) THEN
      IF(.NOT.found) THEN
        first = i
      END IF
      found=.TRUE.
    ELSE
      IF(found) THEN
        last = i-1
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary2
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary2(num_array1,num_array2,ntot,first,last,search1,search2)


    INTEGER, POINTER                         :: num_array1(:), num_array2(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

  found=.FALSE.
  first=0
  last=ntot

  CALL find_boundary(num_array1,ntot,tfirst,tlast,search1)

  last=tlast

  DO i=tfirst,tlast
    IF(num_array2(i)==search2) THEN
      IF(.NOT.found) THEN
        first = i
      END IF
      found=.TRUE.
    ELSE
      IF(found) THEN
        last = i-1
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary2

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary3
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary3(num_array1,num_array2,num_array3,ntot,first,last,&
                          search1,search2,search3)


    INTEGER, POINTER                         :: num_array1(:), num_array2(:), &
                                                num_array3(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2, search3

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

  found=.FALSE.
  first=0
  last=ntot

  CALL find_boundary2(num_array1,num_array2,ntot,tfirst,tlast,search1,search2)

  last=tlast

  DO i=tfirst,tlast
    IF(num_array3(i)==search3) THEN
      IF(.NOT.found) THEN
        first = i
      END IF
      found=.TRUE.
    ELSE
      IF(found) THEN
        last = i-1
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary3


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
SUBROUTINE topology_set_atm_mass(topology,globenv)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_set_atm_mass', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000
  TYPE(atom_info_type),POINTER     :: atom_info

  CHARACTER(len=default_string_length)       :: label_atmname
  CHARACTER(LEN=2), POINTER                :: atm_symbol(:)
  CHARACTER(LEN=2)                         :: upper_sym_1,upper_sym_2
  LOGICAL                                  :: element_found
  INTEGER                                  :: ielem_found
  INTEGER                                  :: first, i, iatom, ibond, &
       ielem, istat, iw, j, jatom, &
       k, last, natom, nbond, nelem, &
       search
  INTEGER, POINTER                         :: tbond_a(:), tbond_b(:)
  LOGICAL                                  :: PM1, PM2, PM3, PM4, only_qm
  REAL(KIND=dp)                                :: dr(3), r2, rbond, rbond2, &
       rr(3)
  REAL(KIND=dp), POINTER                       :: covalent_radius(:), &
       vdw_radius(:)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_set_atm_mass"

  atom_info => topology%atom_info

  natom = topology%natoms

  !Get the temporary array that we need
  ALLOCATE(atm_symbol(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('set_atm_mass','allocate:: atm_symbol',natom)
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
        WRITE(iw,'(A)') "WARNING : topology_set_atm_mass: The element " &
                        //atm_symbol(iatom)//" was not found in the periodic table"
    ENDIF

    atom_info%atm_mass(iatom) = ptable(ielem_found)% amass
    IF(PM2) WRITE(iw,'(7X,A,A5)') "iatom = ",ptable(ielem_found)%symbol
  END DO

  DEALLOCATE(atm_symbol,STAT=istat)
  IF(istat/=0) CALL stop_memory('set_atm_mass','deallocate:: atm_symbol',natom)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Exiting  topology_set_atm_mass"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
SUBROUTINE topology_generate_bond(topology,globenv)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bond', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000
  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

  CHARACTER(len=default_string_length)       :: label_atmname
    CHARACTER(LEN=2), POINTER                :: atm_symbol(:)
    CHARACTER(LEN=2)                         :: upper_sym_1,upper_sym_2
    LOGICAL                                  :: element_found
    INTEGER                                  :: ielem_found
    INTEGER                                  :: first, i, iatom, ibond, &
                                                ielem, istat, iw, j, jatom, &
                                                k, last, natom, nbond, nelem, &
                                                search
    INTEGER, POINTER                         :: tbond_a(:), tbond_b(:)
    LOGICAL                                  :: PM1, PM2, PM3, PM4, only_qm
    REAL(KIND=dp)                                :: dr(3), r2, rbond, rbond2, &
                                                rr(3)
    REAL(KIND=dp), POINTER                       :: covalent_radius(:), &
                                                vdw_radius(:)
    INTEGER :: cbond

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_bond"

  atom_info => topology%atom_info
  conn_info => topology%conn_info
  
  cbond = 0

  natom = topology%natoms

  !Get the temporary array that we need
  ALLOCATE(tbond_a(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_a',nblock)
  ALLOCATE(tbond_b(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_b',nblock)
  ALLOCATE(atm_symbol(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','atm_symbol',natom)
  ALLOCATE(covalent_radius(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','covalent_radius',natom)
  ALLOCATE(vdw_radius(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','vdw_radius',natom)
!------------------------------------------------------------------------------
  DO iatom=1,natom
    atm_symbol(iatom) = ''
    IF(ASSOCIATED(atom_info%element)) THEN
      atm_symbol(iatom) = atom_info%element(iatom)
    END IF
    !
    ! In general %element should always be set... I added this stuff for qmmm 
    ! calculation because if you read coordinates from cp2k input %element is not
    ! setted up!! We should change this...
    !
  !  IF(atm_symbol(iatom) == '') THEN
  !     label_atmname = atom_info%label_atmname(iatom)
  !     only_qm = qmmm_ff_precond_only_qm(id1=label_atmname)
  !     atm_symbol(iatom) = label_atmname
  !  END IF
  END DO
!------------------------------------------------------------------------------
  !Set up the arrays that we need for this routine
  tbond_a(:) = 0
  tbond_b(:) = 0

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
        !EXIT
      END IF
    END DO
    IF (.NOT. element_found) THEN
        WRITE(iw,'(A)') "WARNING : topology_generate_bond: The element " &
                        //atm_symbol(iatom)//" was not found in the periodic table"
    ENDIF
    covalent_radius(iatom) = ptable(ielem_found)%covalent_radius
    CALL convert_to_cp2k_units("ANGSTROM",length=covalent_radius(iatom))

    vdw_radius(iatom) = ptable(ielem_found)%vdw_radius
    CALL convert_to_cp2k_units("ANGSTROM",length=vdw_radius(iatom))

    atom_info%atm_mass(iatom) = ptable(ielem_found)% amass
    IF(PM2) WRITE(iw,'(7X,A,A5)') "iatom = ",ptable(ielem_found)%symbol
  END DO
!------------------------------------------------------------------------------
  !Find bonds according to periodic table
  nbond = 0
  DO iatom=1,natom
    DO jatom=iatom+1,natom
      dr(1) = atom_info%r(1,iatom) - atom_info%r(1,jatom)
      dr(2) = atom_info%r(2,iatom) - atom_info%r(2,jatom)
      dr(3) = atom_info%r(3,iatom) - atom_info%r(3,jatom)

      rr = pbc ( dr, topology%cell )
      r2 = DOT_PRODUCT(rr,rr)

      IF(r2 <= 0.01_dp) THEN
        CALL stop_program ("topology_generate_bond",&
                           "bond distance between atoms less then 0.1")
      END IF

      IF(topology%bondparm_type=="COVALENT") THEN
        rbond = covalent_radius(iatom)+covalent_radius(jatom)
      ELSE IF(topology%bondparm_type=="VDW") THEN
        rbond = MAX(vdw_radius(iatom),vdw_radius(jatom))
      ELSE
        CALL stop_program ("topology_generate_bond",&
                           "illegal bondparm_type")
      END IF
      rbond2 = rbond*rbond

      !Look for the right factor to use...
      rbond2 = rbond2*(topology%bondparm_factor)

      !Test the distance to the sum of the covalent radius
      IF(r2 < rbond2) THEN
        nbond=nbond+1
        IF(nbond > SIZE(tbond_a)) THEN
          CALL reallocate(tbond_a,1,nbond+nblock)
          CALL reallocate(tbond_b,1,nbond+nblock)
        END IF
        tbond_a(nbond)=iatom
        tbond_b(nbond)=jatom
        IF(PM2) WRITE(iw,'(7X,A,I5,I5,I5)') "bond ",nbond,iatom,jatom

        IF(topology%para_res) THEN
        IF((atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)).OR.&
           (atom_info%label_resid(iatom)/=atom_info%label_resid(jatom)).OR.&
           (atom_info%label_resname(iatom)/=atom_info%label_resname(jatom)))THEN
            IF(PM1) WRITE(iw,*) "      PARA_RES, bond between molecules atom ",&
                                iatom,jatom
            cbond = cbond + 1
            CALL reallocate(conn_info%c_bond_a,1,cbond+1)
            CALL reallocate(conn_info%c_bond_b,1,cbond+1)
            conn_info%c_bond_a(cbond) = iatom
            conn_info%c_bond_b(cbond) = jatom
        END IF
        ELSE
        IF(atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)) THEN
            CALL stop_program ("topology_generate_bond",&
                               "bonds between different molecule types??? ")
        END IF
        END IF

      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Copy the stuff from temporary array into the more permanent one
  ALLOCATE(conn_info%bond_a(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_a',nbond)
  ALLOCATE(conn_info%bond_b(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_b',nbond)
  DO ibond=1,nbond
    conn_info%bond_a(ibond) = tbond_a(ibond)
    conn_info%bond_b(ibond) = tbond_b(ibond)
  END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  DEALLOCATE(tbond_a,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_a')
  DEALLOCATE(tbond_b,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_b')
  DEALLOCATE(atm_symbol,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','atm_symbol')
  DEALLOCATE(covalent_radius,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','covalent_radius')
  DEALLOCATE(vdw_radius,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','vdw_radius')

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_bond"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
SUBROUTINE topology_generate_bend(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bend', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000
  TYPE(connectivity_info_type),POINTER     :: conn_info

    INTEGER                                  :: nbond,ibond,jbond,&
                                                ntheta,itheta,jtheta,&
                                                iw,istat
    LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_bend"

  conn_info => topology%conn_info

  nbond = 0
  IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
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
        IF(PM1) WRITE(iw,*) "    Duplicate BEND ",itheta,&
                            conn_info%theta_a(itheta),&
                            conn_info%theta_b(itheta),&
                            conn_info%theta_c(itheta)
        IF(PM1) WRITE(iw,*) "    Duplicate BEND ",jtheta,&
                            conn_info%theta_a(jtheta),&
                            conn_info%theta_b(jtheta),&
                            conn_info%theta_c(jtheta)
        CALL stop_program("topology_generate_bend","duplicate bend")
      END IF
    END DO
  END DO
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_bend"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

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
SUBROUTINE topology_generate_ub(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_ub', &
      routineP = moduleN//':'//routineN

    TYPE(connectivity_info_type),POINTER     :: conn_info

    INTEGER                                  :: itheta, ntheta
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER                                  :: iw

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_ub"

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
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_ub"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_ub
 
 

!
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
SUBROUTINE topology_generate_dihe(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_dihe', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

  TYPE(connectivity_info_type),POINTER     :: conn_info

    INTEGER                                  :: iphi, jphi, iw, nphi, ntheta
    INTEGER, POINTER                         :: tmp_a(:), tmp_b(:), tmp_c(:), &
                                                tmp_d(:)
    LOGICAL                                  :: PM1, PM2, PM3, PM4

    INTEGER                                  :: itheta,jtheta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_dihe"

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
        IF(PM1) THEN
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
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_dihe"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
SUBROUTINE topology_generate_impr(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_impr', &
      routineP = moduleN//':'//routineN

    TYPE(atom_info_type),POINTER             :: atom_info
    TYPE(connectivity_info_type),POINTER     :: conn_info

    INTEGER                                  :: iatom, natom
    INTEGER                                  :: ibond, nbond
    INTEGER                                  :: nimpr
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER                                  :: iw,istat

    INTEGER                                  :: counter,atm_num(3)
    INTEGER, DIMENSION(:), POINTER           :: checkme

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_impr"

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

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
      DEALLOCATE(checkme,STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_generate_impr','checkme')
    END IF
  END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_impr"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
SUBROUTINE topology_generate_onfo(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_onfo', &
      routineP = moduleN//':'//routineN

    TYPE(connectivity_info_type),POINTER     :: conn_info

    INTEGER                                  :: iphi, nphi
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER                                  :: iw

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_onfo"

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
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_onfo"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
                                      topology,globenv)


    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_connectivity_pack', &
      routineP = moduleN//':'//routineN

  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

    CHARACTER(LEN=default_string_length)     :: nhcopt,name
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                counter, first, i, istat, iw, &
                                                j, k,last, natom, nbond, nmol, &
                                                nphi, ntheta, nub, nimpr, nonfo
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, DIMENSION(:), POINTER           :: checkme, first_list, last_list
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(ub_type), DIMENSION(:), POINTER     :: ub_list
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list
    TYPE(onfo_type), DIMENSION(:), POINTER   :: onfo_list
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list
    TYPE(local_molecule_type), POINTER       :: lmi
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    INTEGER                                  :: ikind,jkind,ityp,inum,ires,imol
    INTEGER                                  :: ibond,c_intra
    LOGICAL                                  :: found
    INTEGER                                  :: nsgf,nelectron,handle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  CALL timeset(routineN,'I','',handle)

  IF(PM1) WRITE(iw,*) "  Entering topology_connectivity_pack"

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
    IF(PM2) WRITE(iw,*) "    Allocated molecule_kind_set, Dimenstion of ",&
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
    IF(PM2) WRITE(iw,*) "    Allocated molecule_set, dimenstion of ",&
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
    name = TRIM(atom_info%label_molname(j))//"-"//&
           TRIM(atom_info%label_resname(j))
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
         kind_number=ikind,&
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
    ALLOCATE(molecule_list(imol-counter),STAT=istat)
    IF(istat/=0) &
      CALL stop_memory('connectivity_pack','molecule_list',imol-counter+1)
    DO i=1,SIZE(molecule_list)
      molecule_list(i)=i+counter-1
    END DO
    molecule_kind => molecule_kind_set(ikind)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list)
    IF(PM2) WRITE(iw,*) "      molecule_list",ikind,molecule_list(:)
    counter=imol
  END DO
  DO ikind=1,SIZE(molecule_kind_set)-topology%nmol_type
    ALLOCATE(molecule_list(1),STAT=istat)
    IF(istat/=0) &
      CALL stop_memory('connectivity_pack','molecule_list',imol-counter+1)
    molecule_list(1)=imol+ikind
    molecule_kind => molecule_kind_set(ikind+topology%nmol_type)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list)
    IF(PM2) WRITE(iw,*) "      molecule_list",ikind+topology%nmol_type,molecule_list(:)
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
  ALLOCATE(first_list(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('connectivity_pack','frist_list',topology%nmol)
  ALLOCATE(last_list(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('connectivity_pack','last_list',topology%nmol)
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
  DEALLOCATE(first_list,STAT=istat)
  IF(istat/=0) CALL stop_memory('connectivity_pack','frist_list')
  DEALLOCATE(last_list,STAT=istat)
  IF(istat/=0) CALL stop_memory('connectivity_pack','last_list')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Set and NULLIFY the molecule_set(imol)%lmi via set_molecule_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_set)
    molecule => molecule_set(i)
    NULLIFY(lmi)
    ALLOCATE(lmi,STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','lmi',1)
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
      ALLOCATE(atom_list(natom),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','atom_list',natom)
      DO i=1,natom
        !Atomic kind information will be filled in (PART 2)
        NULLIFY(atom_list(i)%atomic_kind)
        atom_list(i)%name=atom_info%label_atmname(i+first-1)
        IF(PM2) WRITE(iw,'(5X,A,I5,I5,1X,A5)') "atom_list ",ikind,i,atom_list(i)%name
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
      ALLOCATE(checkme(nbond),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nbond)
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
      IF(PM2) WRITE(iw,*) "    Total number bonds for molecule type ",&
                          SUM(checkme)+c_intra
      
      ibond = SUM(checkme) + c_intra
      ALLOCATE(bond_list(ibond),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','bond_list',ibond)
      ibond = 0
      DO j=1,nbond
        IF(checkme(j)==1) THEN
          ibond = ibond + 1
          bond_list(ibond)%a=conn_info%bond_a(j)-first+1
          bond_list(ibond)%b=conn_info%bond_b(j)-first+1
          bond_list(ibond)%name="BOND"
          !point this to the right bond_kind_type if using force field
          NULLIFY(bond_list(ibond)%bond_kind)
          IF(PM2) THEN
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
            IF(PM2) THEN
              WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",&
                          ikind,"c_bond",&
                          conn_info%c_bond_a(j),&
                          conn_info%c_bond_b(j),&
                          "offset number at",&
                          conn_info%c_bond_a(j)-first+1,&
                          conn_info%c_bond_b(j)-first+1
            ibond = ibond + 1
            END IF
          END IF
        END DO
      END IF
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             nbond=SIZE(bond_list),bond_list=bond_list)
      DEALLOCATE(checkme,STAT=istat)
      IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
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
      c_intra = 0      !number of extra connector intra-mol interactions
      counter = i
      ntheta = 0
      CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
      IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
      IF(ASSOCIATED(conn_info%theta_a)) ntheta = SIZE(conn_info%theta_a)
      ALLOCATE(checkme(ntheta),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nbond)
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
                    (atm_b==conn_info%theta_c(k)) ) ) THEN
                checkme(k)=1
              END IF
            END DO
          END IF
        END DO
      END IF

      IF(PM2) WRITE(iw,*) "    Total number bends for molecule type ",&
                          SUM(checkme)+c_intra
      i = SUM(checkme) + c_intra
      ALLOCATE(bend_list(i),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','bend_list',i)
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
          IF(PM2) THEN
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
      DEALLOCATE(checkme,STAT=istat)
      IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
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
      c_intra = 0      !number of extra connector intra-mol interactions
      counter = i
      nub = 0
      CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
      IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
      IF(ASSOCIATED(conn_info%ub_a)) nub = SIZE(conn_info%ub_a)
      ALLOCATE(checkme(nub),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nub)
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
                    (atm_b==conn_info%ub_c(k)) ) ) THEN
                checkme(k)=1
              END IF
            END DO
          END IF
        END DO
      END IF

      IF(PM2) WRITE(iw,*) "    Total number UB for molecule type ",&
                          SUM(checkme)+c_intra
      i = SUM(checkme) + c_intra
      ALLOCATE(ub_list(i),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','ub_list',i)
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
          IF(PM2) THEN
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
      DEALLOCATE(checkme,STAT=istat)
      IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
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
      c_intra = 0      !number of extra connector intra-mol interactions
      counter = i
      nphi = 0
      CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
      IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
      IF(ASSOCIATED(conn_info%phi_a)) nphi = SIZE(conn_info%phi_a)
      ALLOCATE(checkme(nphi),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nbond)
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
                    (atm_b==conn_info%phi_d(k)) ) ) THEN
                checkme(k)=1
              END IF
            END DO
          END IF
        END DO
      END IF

      IF(PM2) WRITE(iw,*) "    Total number TOR for molecule type ",&
                          SUM(checkme)+c_intra
      i = SUM(checkme) + c_intra
      ALLOCATE(torsion_list(i),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','phi_list',i)
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
          IF(PM2) THEN
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
      DEALLOCATE(checkme,STAT=istat)
      IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
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
      c_intra = 0      !number of extra connector intra-mol interactions
      counter = i
      nimpr = 0
      CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
      IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
      IF(ASSOCIATED(conn_info%impr_a)) nimpr = SIZE(conn_info%impr_a)
      ALLOCATE(checkme(nimpr),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nimpr)
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

      IF(PM2) WRITE(iw,*) "    Total number IMPR for molecule type ",&
                          SUM(checkme)+c_intra
      i = SUM(checkme) + c_intra
      ALLOCATE(impr_list(i),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','impr_list',i)
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
          IF(PM2) THEN
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
      DEALLOCATE(checkme,STAT=istat)
      IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
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
      c_intra = 0      !number of extra connector intra-mol interactions
      counter = i
      nonfo = 0
      CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
      IF(ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
      IF(ASSOCIATED(conn_info%onfo_a)) nonfo = SIZE(conn_info%onfo_a)
      ALLOCATE(checkme(nonfo),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nonfo)
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
                    (atm_b==conn_info%phi_d(k)) ) ) THEN
                checkme(k)=1
              END IF
            END DO
          END IF
        END DO
      END IF

      IF(PM2) WRITE(iw,*) "    Total number ONFO for molecule type ",&
                          SUM(checkme)+c_intra
      i = SUM(checkme) + c_intra
      ALLOCATE(onfo_list(i),STAT=istat)
      IF(istat/=0) CALL stop_memory('connectivity_pack','onfo_list',i)
      i = 0
      DO j=1,nonfo
        IF(checkme(j)==1) THEN
          i = i + 1
          onfo_list(i)%a=conn_info%onfo_a(j)-first+1
          onfo_list(i)%b=conn_info%onfo_b(j)-first+1
          onfo_list(i)%name="ONFO"
          !point this to the right ONFO if using force field
          NULLIFY(onfo_list(i)%onfo_kind)
          IF(PM2) THEN
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
      DEALLOCATE(checkme,STAT=istat)
      IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
    END IF
  END DO


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_connectivity_pack"

  CALL timestop(0.0_dp,handle)

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
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
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

SUBROUTINE topology_coordinate_pack (particle_set,atomic_kind_set,&
                              molecule_kind_set,molecule_set,&
                              topology,globenv,qmmm,qmmm_env,error)


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
      INTENT(IN), OPTIONAL                   :: globenv
    TYPE(qmmm_env_mm_type), POINTER, OPTIONAL        :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    LOGICAL, OPTIONAL, INTENT(IN) :: qmmm

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_coordinate_pack', &
      routineP = moduleN//':'//routineN

  TYPE(atom_info_type),POINTER             :: atom_info

    CHARACTER(LEN=default_string_length), &
      POINTER                                :: work(:)
    CHARACTER(LEN=default_string_length)     :: my_elp_name, my_elp_description(2)
    INTEGER                                  :: counter, first, i, ikind, &
                                                istat, iw, j, k, last, natom
    INTEGER, DIMENSION(:), POINTER           :: atomkind,atomlist, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    REAL(KIND=dp), DIMENSION(:), POINTER         :: charge, mass
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list,ex_bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list,ex_bond_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list,ex_torsion_list
    INTEGER, DIMENSION(:), POINTER :: list,list2
    INTEGER                                  :: imolecule,first_atom,last_atom,&
                                                nbond,nbend,ntorsion

    CHARACTER(LEN=default_string_length)     :: atmname
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: qeff, alpha_core_charge
    TYPE(elp_potential_type), POINTER        :: elp_potential
    LOGICAL :: my_qmmm
    
    INTEGER :: iatom,handle

    INTEGER, DIMENSION(:), POINTER :: first_mol,last_mol,natom_mol,nbond_mol,nbend_mol,ntorsion_mol
    INTEGER :: imol_now,imol_minus,imol_plus
    INTEGER :: ibond,ibend,itorsion,offset
    LOGICAL :: mol_change,mol_change2


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  CALL timeset(routineN,'I','',handle)

  IF(PM1) WRITE(iw,*) "  Entering topology_coordinate_pack"
  my_qmmm = .FALSE.
  IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm=qmmm

  atom_info => topology%atom_info


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Determine topology%[natom_type,atom_names] and save mass(natom_type) 
  !-----------------------------------------------------------------------------
  counter=0
  ALLOCATE(work(topology%natoms),mass(topology%natoms),STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','work,mass',topology%natoms)
  work(:)=""
  DO i = 1, SIZE(molecule_kind_set)
    DO j=1,SIZE(molecule_set)
      molecule => molecule_set(j)
      molecule_kind => molecule_set(j)%molecule_kind
      CALL get_molecule_kind(molecule_kind=molecule_kind,kind_number=ikind)
      IF(ikind==i) EXIT
    END DO
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list)
    CALL get_molecule(molecule=molecule,&
                           first_atom=first,last_atom=last)
    DO j = 1, natom
      IF(str_search(work,counter,atom_list(j)%name)==0)THEN
         counter = counter + 1
         work(counter) = atom_list(j)%name
         mass(counter) = atom_info%atm_mass(first+j-1)
        IF(PM2) WRITE(iw,'(7X,A,1X,A5,F10.5)') &
                "NEW ATOMIC KIND",work(counter),mass(counter)
      END IF
    END DO
  END DO
  topology%natom_type = counter
  ALLOCATE(atom_info%atom_names(topology%natom_type),STAT=istat)
  IF(istat/=0) CALL stop_memory ('coordinate_pack','atom_names')
  atom_info%atom_names(1:counter)=work(1:counter)
  DEALLOCATE(work,STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','work')
  CALL reallocate(mass,1,counter)
  IF(PM2) &
    WRITE(iw,'(5X,A,I3)') "Total Number of Atomic Kinds = ",topology%natom_type


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Allocate the data structure for the atomic kind information 
  !-----------------------------------------------------------------------------
  NULLIFY(atomic_kind_set)
  CALL allocate_atomic_kind_set(atomic_kind_set,topology%natom_type)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3.  Allocate the data structure for the atomic information
  !-----------------------------------------------------------------------------
  NULLIFY(particle_set)
  CALL allocate_particle_set(particle_set,topology%natoms)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Set the atomic_kind_set(ikind)%[name,kind_number,mass]
  !-----------------------------------------------------------------------------
  DO i=1,topology%natom_type
    atomic_kind => atomic_kind_set(i)
    mass(i) = mass(i) * massunit
    CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=i,&
                         name=atom_info%atom_names(i),mass=mass(i)) 
  END DO
  DEALLOCATE(mass,STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','mass')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Determine number of atom of each kind (ie natom_of_kind and kind_of)
  !-----------------------------------------------------------------------------
  ALLOCATE(kind_of(topology%natoms),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_of')
  ALLOCATE(natom_of_kind(topology%natom_type),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_of_kind')
  kind_of(:) = 0
  natom_of_kind(:) = 0
  DO i=1,topology%natom_type
    DO j=1,topology%natoms
      IF(atom_info%atom_names(i)==atom_info%label_atmname(j)) THEN
        natom_of_kind(i)=natom_of_kind(i)+1
        IF(kind_of(j)==0) kind_of(j)=i
      END IF
    END DO
  END DO

  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Set the atom_kind_set(ikind)%[natom,atom_list]
  !-----------------------------------------------------------------------------
  DO i=1,topology%natom_type
    atomic_kind => atomic_kind_set(i)
    NULLIFY (atomlist)
    ALLOCATE(atomlist(natom_of_kind(i)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomlist')
    counter=0
    DO j=1,topology%natoms
      IF(kind_of(j)==i) THEN
        counter=counter+1
        atomlist(counter)=j
      END IF
    END DO
    IF(PM2) THEN 
      WRITE(iw,'(A,I6,A)') "      Atomic kind ",i," contains particles"
      DO J=1,SIZE(atomlist)  
         IF (MOD(J,5).EQ.0) THEN ! split long lines
            write(iw,'(I12)') atomlist(J)
         ELSE
            write(iw,'(I12)',ADVANCE="NO") atomlist(J)
         ENDIF
      ENDDO
      WRITE(iw,*)
    ENDIF
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         natom=natom_of_kind(i),atom_list=atomlist)
    DEALLOCATE(atomlist,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomlist')
  END DO
  DEALLOCATE(natom_of_kind,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_of_kind')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Fill in some more infor for particle_set and lpi (local part info)
  !-----------------------------------------------------------------------------
  DO i=1,topology%natoms
     ikind = kind_of(i)
     IF (PM2) THEN
        WRITE(iw,*)"atom number :: ",i,"kind number ::",ikind
     END IF
     particle_set(i)%atomic_kind => atomic_kind_set(ikind)
     particle_set(i)%r(:) = atom_info%r(:,i)
     NULLIFY(particle_set(i)%lpi)
     ALLOCATE(particle_set(i)%lpi)
     NULLIFY(particle_set(i)%lpi%exclusion)
     ALLOCATE(particle_set(i)%lpi%exclusion(SIZE(atomic_kind_set)))
  END DO
  DEALLOCATE(kind_of,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_of')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Fill in the particle_set%lpi%exclusion%list_vdw
  !-----------------------------------------------------------------------------
  ALLOCATE(atomkind(SIZE(particle_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','atomkind')
  ALLOCATE(natom_mol(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','natom_mol')
  ALLOCATE(nbond_mol(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','nbond_mol')
  ALLOCATE(nbend_mol(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','nbend_mol')
  ALLOCATE(ntorsion_mol(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ntorsion_mol')
  ALLOCATE(first_mol(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','first_mol')
  ALLOCATE(last_mol(SIZE(molecule_set)),STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','last_mol')
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
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom_mol(i),&
                           nbond=nbond_mol(i),&
                           nbend=nbend_mol(i),&
                           ntorsion=ntorsion_mol(i))
  END DO
  imol_now=1
  mol_change=.FALSE.
  mol_change2=.FALSE.
  DO iatom=1,SIZE(particle_set)
    IF(iatom==first_mol(imol_now)) mol_change =.TRUE.
    IF(iatom==last_mol(imol_now)) mol_change2=.TRUE.
    IF(mol_change) THEN
      imol_minus = MAX(1,imol_now-topology%exclude_check)
      imol_plus  = MIN(SIZE(molecule_set),imol_now+topology%exclude_check)
      natom = SUM(natom_mol(imol_minus:imol_plus))
      nbond = SUM(nbond_mol(imol_minus:imol_plus))
      nbend = SUM(nbend_mol(imol_minus:imol_plus))
      ntorsion = SUM(ntorsion_mol(imol_minus:imol_plus))

      ALLOCATE(atomlist(natom),STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','atomlist')
      ALLOCATE(ex_bond_list(nbond),STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bond_list')
      ALLOCATE(ex_bend_list(nbend),STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bend_list')
      ALLOCATE(ex_torsion_list(ntorsion),STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_torsion_list')
      ibond = 1
      ibend = 1
      itorsion = 1
      DO i=imol_minus,imol_plus
        molecule => molecule_set(i)
        CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               bond_list=bond_list,&
                               bend_list=bend_list,&
                               torsion_list=torsion_list)
        ex_bond_list(ibond:ibond+SIZE(bond_list)-1)%a = &
                    bond_list(:)%a+first_mol(i)-1
        ex_bond_list(ibond:ibond+SIZE(bond_list)-1)%b = &
                    bond_list(:)%b+first_mol(i)-1
        ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%a = &
                    bend_list(:)%a+first_mol(i)-1
        ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%b = &
                    bend_list(:)%b+first_mol(i)-1
        ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%c = &
                    bend_list(:)%c+first_mol(i)-1
        ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%a = &
                       torsion_list(:)%a+first_mol(i)-1
        ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%b = &
                       torsion_list(:)%b+first_mol(i)-1
        ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%c = &
                       torsion_list(:)%c+first_mol(i)-1
        ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%d = &
                       torsion_list(:)%d+first_mol(i)-1
        ibond = ibond + SIZE(bond_list)
        ibend = ibend + SIZE(bend_list)
        itorsion = itorsion + SIZE(torsion_list)
      END DO
      offset = first_mol(imol_minus)
      mol_change=.FALSE.
    END IF

    atomlist(:) = 1 
    atomlist(iatom-offset+1) = 0    !exclude yourself???

    IF(topology%exclude_vdw=="1-2" .OR. &
       topology%exclude_vdw=="1-3" .OR. &
       topology%exclude_vdw=="1-4") THEN
      DO j=1,nbond
         IF (PM2) WRITE(iw,*)"bond_list ::",&
                             ex_bond_list(j)%a,&
                             ex_bond_list(j)%b
         IF(  (ex_bond_list(j)%a==iatom).OR. &
              (ex_bond_list(j)%b==iatom) ) THEN
            atomlist(ex_bond_list(j)%a-offset+1) = 0
            atomlist(ex_bond_list(j)%b-offset+1) = 0
         END IF
      END DO
    END IF

    IF(topology%exclude_vdw=="1-3" .OR. &
       topology%exclude_vdw=="1-4") THEN
      DO j=1,nbend
         IF (PM2) WRITE(iw,*)"bend_list ::",&
                             ex_bend_list(j)%a,&
                             ex_bend_list(j)%b,&
                             ex_bend_list(j)%c
         IF(  (ex_bend_list(j)%a==iatom).OR. &
              (ex_bend_list(j)%c==iatom) ) THEN
            atomlist(ex_bend_list(j)%a-offset+1) = 0
            atomlist(ex_bend_list(j)%c-offset+1) = 0
         END IF
      END DO
    END IF

    IF(topology%exclude_vdw=="1-4") THEN
      DO j=1,ntorsion
         IF (PM2) WRITE(iw,*)"tors_list ::",&
                             ex_torsion_list(j)%a,&
                             ex_torsion_list(j)%b,&
                             ex_torsion_list(j)%c,&
                             ex_torsion_list(j)%d
         IF(  (ex_torsion_list(j)%a==iatom).OR. &
              (ex_torsion_list(j)%d==iatom) ) THEN
            atomlist(ex_torsion_list(j)%a-offset+1) = 0
            atomlist(ex_torsion_list(j)%d-offset+1) = 0
         END IF
      END DO
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
             IF ((counter >= first_atom).AND.(counter <= last_atom)) &
                  atomlist(counter-offset+1) = 0
          END DO
       END IF
    END IF

    DO j=1,SIZE(atomic_kind_set)
      counter = 0
      NULLIFY(list)
      DO k=1,natom
        IF((atomkind(k+offset-1)==j).AND.(atomlist(k)==0)) THEN
          counter = counter + 1
          CALL reallocate(list,1,counter)
          list(counter) = k+offset-1
        END IF
      END DO
      particle_set(iatom)%lpi%exclusion(j)%list_vdw => list
      IF(topology%exclude_vdw==topology%exclude_ei) THEN
        NULLIFY(list2)
        IF(counter > 0) THEN
          CALL reallocate(list2,1,counter)
          list2(:) = list(:)
        END IF
        particle_set(iatom)%lpi%exclusion(j)%list_ei => list2
      END IF
      IF (PM2) THEN 
         IF (ASSOCIATED(list)) &
         WRITE(iw,*)"exclusion list_vdw :: ",&
                   "atom num :",iatom,"kind num. :",j,"exclusion list ::",&
                   list
      END IF
    END DO
    IF(mol_change2) THEN
      imol_now = imol_now + 1
      DEALLOCATE(atomlist,STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','atomlist')
      DEALLOCATE(ex_bond_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bond_list')
      DEALLOCATE(ex_bend_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bend_list')
      DEALLOCATE(ex_torsion_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('topology_coordinate_pack','x_torsion_list')
      mol_change2=.FALSE.
    END IF
  END DO
  DEALLOCATE(atomkind,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomkind')
  DEALLOCATE(natom_mol,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_mol')
  DEALLOCATE(nbond_mol,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','nbond_mol')
  DEALLOCATE(nbend_mol,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','nbend_mol')
  DEALLOCATE(ntorsion_mol,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','ntorsion_mol')
  DEALLOCATE(first_mol,STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','first_mol')
  DEALLOCATE(last_mol,STAT=istat)
  IF(istat/=0) CALL stop_memory('topology_coordinate_pack','last_mol')

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. Fill in the particle_set%lpi%exclusion%list_ei
  !-----------------------------------------------------------------------------
  IF(topology%exclude_vdw/=topology%exclude_ei) THEN
    ALLOCATE(atomkind(SIZE(particle_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','atomkind')
    ALLOCATE(natom_mol(SIZE(molecule_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','natom_mol')
    ALLOCATE(nbond_mol(SIZE(molecule_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','nbond_mol')
    ALLOCATE(nbend_mol(SIZE(molecule_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','nbend_mol')
    ALLOCATE(ntorsion_mol(SIZE(molecule_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ntorsion_mol')
    ALLOCATE(first_mol(SIZE(molecule_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','first_mol')
    ALLOCATE(last_mol(SIZE(molecule_set)),STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','last_mol')
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
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom_mol(i),&
                             nbond=nbond_mol(i),&
                             nbend=nbend_mol(i),&
                             ntorsion=ntorsion_mol(i))
    END DO
    imol_now=1
    mol_change=.FALSE.
    mol_change2=.FALSE.
    DO iatom=1,SIZE(particle_set)
  
      IF(iatom==first_mol(imol_now)) mol_change =.TRUE.
      IF(iatom==last_mol(imol_now-1)) mol_change2=.TRUE.
      IF(mol_change) THEN
        imol_minus = MAX(1,imol_now-topology%exclude_check)
        imol_plus  = MIN(SIZE(molecule_set),imol_now+topology%exclude_check)
        natom = SUM(natom_mol(imol_minus:imol_plus))
        nbond = SUM(nbond_mol(imol_minus:imol_plus))
        nbend = SUM(nbend_mol(imol_minus:imol_plus))
        ntorsion = SUM(ntorsion_mol(imol_minus:imol_plus))
  
        ALLOCATE(atomlist(natom),STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','atomlist')
        ALLOCATE(ex_bond_list(nbond),STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bond_list')
        ALLOCATE(ex_bend_list(nbend),STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bend_list')
        ALLOCATE(ex_torsion_list(ntorsion),STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_torsion_list')
        ibond = 1
        ibend = 1
        itorsion = 1
        DO i=imol_minus,imol_plus
          molecule => molecule_set(i)
          CALL get_molecule(molecule=molecule,molecule_kind=molecule_kind)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
                                 bond_list=bond_list,&
                                 bend_list=bend_list,&
                                 torsion_list=torsion_list)
          ex_bond_list(ibond:ibond+SIZE(bond_list)-1)%a = &
                      bond_list(:)%a+first_mol(i)-1
          ex_bond_list(ibond:ibond+SIZE(bond_list)-1)%b = &
                      bond_list(:)%b+first_mol(i)-1
          ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%a = &
                      bend_list(:)%a+first_mol(i)-1
          ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%b = &
                      bend_list(:)%b+first_mol(i)-1
          ex_bend_list(ibend:ibend+SIZE(bend_list)-1)%c = &
                      bend_list(:)%c+first_mol(i)-1
          ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%a = &
                         torsion_list(:)%a+first_mol(i)-1
          ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%b = &
                         torsion_list(:)%b+first_mol(i)-1
          ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%c = &
                         torsion_list(:)%c+first_mol(i)-1
          ex_torsion_list(itorsion:itorsion+SIZE(torsion_list)-1)%d = &
                         torsion_list(:)%d+first_mol(i)-1
          ibond = ibond + SIZE(bond_list)
          ibend = ibend + SIZE(bend_list)
          itorsion = itorsion + SIZE(torsion_list)
        END DO
        offset = first_mol(imol_minus)
        mol_change=.FALSE.
      END IF
  
      atomlist(:) = 1 
      atomlist(iatom-offset+1) = 0    !exclude yourself???
  
      IF(topology%exclude_ei=="1-2" .OR. &
         topology%exclude_ei=="1-3" .OR. &
         topology%exclude_ei=="1-4") THEN
        DO j=1,nbond
           IF (PM2) WRITE(iw,*)"bond_list ::",&
                               ex_bond_list(j)%a,&
                               ex_bond_list(j)%b
           IF(  (ex_bond_list(j)%a==iatom).OR. &
                (ex_bond_list(j)%b==iatom) ) THEN
              atomlist(ex_bond_list(j)%a-offset+1) = 0
              atomlist(ex_bond_list(j)%b-offset+1) = 0
           END IF
        END DO
      END IF
  
      IF(topology%exclude_ei=="1-3" .OR. &
         topology%exclude_ei=="1-4") THEN
        DO j=1,nbend
           IF (PM2) WRITE(iw,*)"bend_list ::",&
                               ex_bend_list(j)%a,&
                               ex_bend_list(j)%b,&
                               ex_bend_list(j)%c
           IF(  (ex_bend_list(j)%a==iatom).OR. &
                (ex_bend_list(j)%c==iatom) ) THEN
              atomlist(ex_bend_list(j)%a-offset+1) = 0
              atomlist(ex_bend_list(j)%c-offset+1) = 0
           END IF
        END DO
      END IF
  
      IF(topology%exclude_ei=="1-4") THEN
        DO j=1,ntorsion
           IF (PM2) WRITE(iw,*)"tors_list ::",&
                               ex_torsion_list(j)%a,&
                               ex_torsion_list(j)%b,&
                               ex_torsion_list(j)%c,&
                               ex_torsion_list(j)%d
           IF(  (ex_torsion_list(j)%a==iatom).OR. &
                (ex_torsion_list(j)%d==iatom) ) THEN
              atomlist(ex_torsion_list(j)%a-offset+1) = 0
              atomlist(ex_torsion_list(j)%d-offset+1) = 0
           END IF
        END DO
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
               IF ((counter >= first_atom).AND.(counter <= last_atom)) &
                    atomlist(counter)= 0 
            END DO
         END IF
      END IF
  
      DO j=1,SIZE(atomic_kind_set)
        counter = 0
        NULLIFY(list)
        DO k=1,natom
          IF((atomkind(k+offset-1)==j).AND.(atomlist(k)==0)) THEN
            counter = counter + 1
            CALL reallocate(list,1,counter)
            list(counter) = k+offset-1
          END IF
        END DO
        particle_set(iatom)%lpi%exclusion(j)%list_ei => list
        IF (PM2) THEN 
           IF (ASSOCIATED(list)) &
           WRITE(iw,*)"exclusion list_ei :: ",&
                     "atom num :",iatom,"kind num. :",j,"exclusion list ::",&
                     list
        END IF
      END DO
      IF(mol_change2) THEN
        imol_now = imol_now + 1
        DEALLOCATE(atomlist,STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','atomlist')
        DEALLOCATE(ex_bond_list,STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bond_list')
        DEALLOCATE(ex_bend_list,STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','ex_bend_list')
        DEALLOCATE(ex_torsion_list,STAT=istat)
        IF(istat/=0) CALL stop_memory('topology_coordinate_pack','x_torsion_list')
        mol_change2=.FALSE.
      END IF
    END DO
    DEALLOCATE(atomkind,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomkind')
    DEALLOCATE(natom_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_mol')
    DEALLOCATE(nbond_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','nbond_mol')
    DEALLOCATE(nbend_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','nbend_mol')
    DEALLOCATE(ntorsion_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','ntorsion_mol')
    DEALLOCATE(first_mol,STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','first_mol')
    DEALLOCATE(last_mol,STAT=istat)
    IF(istat/=0) CALL stop_memory('topology_coordinate_pack','last_mol')    
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 10. Set the atomic_kind_set()%elp_potentail%[qeff] (PART 1)
  !-----------------------------------------------------------------------------
  IF(globenv%program_name == "FIST") THEN
    DO i=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(i)
      CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
      qeff = -10000.0_dp
      alpha_core_charge = -10000.0_dp
      DO j=1,SIZE(atom_info%label_atmname)
        IF(atmname == atom_info%label_atmname(j)) THEN
          
          qeff = atom_info%atm_charge(j)
          EXIT
        END IF
      END DO
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


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 11. Set the atom_list for molecule_kind in molecule_kind_set (PART 2)
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list)
    DO j=1,natom
      DO k=1,topology%natom_type
        IF(atom_list(j)%name == atom_info%atom_names(k)) THEN
          atom_list(j)%atomic_kind => atomic_kind_set(k)
        END IF
      END DO
    END DO
  END DO


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_coordinate_pack"
  CALL timestop(0.0_dp,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF

END SUBROUTINE topology_coordinate_pack


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
                                      topology,qmmm_env,particle_set,globenv)

  TYPE(molecule_kind_type), DIMENSION(:), &
    POINTER                                :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), &
    POINTER                                :: molecule_set
  TYPE(topology_parameters_type), &
    INTENT(INOUT)                          :: topology
  TYPE(qmmm_env_mm_type), POINTER, OPTIONAL        :: qmmm_env
    TYPE(particle_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: particle_set
  TYPE(global_environment_type), &
    INTENT(IN), OPTIONAL                   :: globenv

  CHARACTER(len=*), PARAMETER :: routineN = 'topology_constraint_pack', &
    routineP = moduleN//':'//routineN

  TYPE(constraint_info_type),POINTER     :: cons_info

  INTEGER                                  :: istat, iw, i, j, k, natom, dum, kk
  INTEGER                                  :: counter,nmolecule,ndist,nbond
  INTEGER                                  :: ng3x3,ng4x6,nfixed_atoms
  INTEGER                                  :: first, last,search1,search2
  LOGICAL                                  :: PM1, PM2, PM3, PM4
  REAL(KIND=dp)                            :: xdum
  INTEGER, DIMENSION(:), POINTER           :: molecule_list
  CHARACTER(LEN=default_string_length)  :: name,element_symbol
  TYPE(atomic_kind_type), POINTER        :: atomic_kind
  TYPE(molecule_kind_type), POINTER        :: molecule_kind
  TYPE(molecule_type), POINTER        :: molecule

  TYPE(dist_constraint_type), DIMENSION(:), POINTER :: dist_list
  TYPE(g3x3_constraint_type), DIMENSION(:), POINTER :: g3x3_list
  TYPE(g4x6_constraint_type), DIMENSION(:), POINTER :: g4x6_list
  TYPE(fixd_constraint_type), DIMENSION(:), POINTER :: fixd_list
  TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
  TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
  TYPE(local_constraint_type), POINTER       :: lci
  TYPE ( local_dist_constraint_type ), DIMENSION ( : ) , POINTER :: ldist
  TYPE ( local_g3x3_constraint_type ), DIMENSION ( : ) , POINTER :: lg3x3
  TYPE ( local_g4x6_constraint_type ), DIMENSION ( : ) , POINTER :: lg4x6
  TYPE(atom_info_type),POINTER     :: atom_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering topology_constraint_pack"

  cons_info => topology%cons_info
  atom_info => topology%atom_info
  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. NULLIFY the molecule_set(imol)%lci via set_molecule_set
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol
    molecule => molecule_set(i)
    NULLIFY(lci)
    ALLOCATE(lci,STAT=istat)
    IF(istat/=0) CALL stop_memory('constraint_pack','lci',1)
    NULLIFY(lci%ldist)
    NULLIFY(lci%lg3x3)
    NULLIFY(lci%lg4x6)
    CALL set_molecule(molecule,lci=lci)
  END DO

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Add more stuff to cons_info%const_dist_* if constraint hydrogen is on
  !-----------------------------------------------------------------------------
  IF(topology%const_hydr) THEN
    ndist = cons_info%nconst_dist
    DO i=1,SIZE(molecule_kind_set)
      molecule_kind => molecule_kind_set(i)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             atom_list=atom_list,&
                             nbond=nbond,bond_list=bond_list)
      DO j=1,SIZE(atom_list)
        atomic_kind => atom_list(j)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,element_symbol=element_symbol) 
        IF((TRIM(name)=="H").OR.(TRIM(element_symbol)=="H")) THEN
          DO k=1,nbond
            IF((bond_list(k)%a==j).OR.(bond_list(k)%b==j)) THEN
              ndist = ndist + 1
              CALL reallocate(cons_info%const_dist_mol,1,ndist)
              CALL reallocate(cons_info%const_dist_a,1,ndist)
              CALL reallocate(cons_info%const_dist_b,1,ndist)
              CALL reallocate(cons_info%const_dist_dab,1,ndist)
              cons_info%const_dist_mol(ndist)  = i
              cons_info%const_dist_a(ndist)    = bond_list(k)%a
              cons_info%const_dist_b(ndist)    = bond_list(k)%b
              cons_info%const_dist_dab(ndist)  = bond_list(k)%bond_kind%r0
            END IF
          END DO
        END IF
      END DO
    END DO
    cons_info%nconst_dist = ndist
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. sort and find duplicates cons_info%const_dist_*
  !-----------------------------------------------------------------------------
  IF(topology%const_dist) THEN
    ndist = cons_info%nconst_dist
    !Sort the first column
    DO i=1,ndist
      DO j=1,ndist
        IF(cons_info%const_dist_mol(i)<cons_info%const_dist_mol(j)) THEN
           !
           dum = cons_info%const_dist_mol(i)
           cons_info%const_dist_mol(i) = cons_info%const_dist_mol(j)
           cons_info%const_dist_mol(j) = dum
           !
           dum = cons_info%const_dist_a(i)
           cons_info%const_dist_a(i) = cons_info%const_dist_a(j)
           cons_info%const_dist_a(j) = dum
           !
           dum = cons_info%const_dist_b(i)
           cons_info%const_dist_b(i) = cons_info%const_dist_b(j)
           cons_info%const_dist_b(j) = dum
           !
           xdum = cons_info%const_dist_dab(i)
           cons_info%const_dist_dab(i) = cons_info%const_dist_dab(j)
           cons_info%const_dist_dab(j) = xdum
           !
        END IF
      END DO
    END DO
    !Sort the second column
    i=1
    DO WHILE (i<ndist)
      search1 = cons_info%const_dist_mol(i)
      CALL find_boundary(cons_info%const_dist_mol,ndist,first,last,search1)
      IF(.NOT.(first==last)) THEN
        DO j=first,last
          DO k=first,last
            IF(cons_info%const_dist_a(j)<cons_info%const_dist_a(k)) THEN
               !
               dum = cons_info%const_dist_mol(j)
               cons_info%const_dist_mol(j) = cons_info%const_dist_mol(k)
               cons_info%const_dist_mol(k) = dum
               !
               dum = cons_info%const_dist_a(j)
               cons_info%const_dist_a(j) = cons_info%const_dist_a(k)
               cons_info%const_dist_a(k) = dum
               !
               dum = cons_info%const_dist_b(j)
               cons_info%const_dist_b(j) = cons_info%const_dist_b(k)
               cons_info%const_dist_b(k) = dum
               !
               xdum = cons_info%const_dist_dab(j)
               cons_info%const_dist_dab(j) = cons_info%const_dist_dab(k)
               cons_info%const_dist_dab(k) = xdum
            END IF
          END DO
        END DO
        i=last+1
      ELSE
        i=i+1
      END IF
    END DO
    !Sort the third column
    i=1
    DO WHILE (i<ndist)
      search1 = cons_info%const_dist_mol(i)
      search2 = cons_info%const_dist_a(i)
      CALL find_boundary2(cons_info%const_dist_mol,cons_info%const_dist_a,&
                          ndist,first,last,search1,search2)
      IF(.NOT.(first==last)) THEN
        DO j=first,last
          DO k=first,last
            IF(cons_info%const_dist_b(j)<cons_info%const_dist_b(k)) THEN
               !
               dum = cons_info%const_dist_mol(j)
               cons_info%const_dist_mol(j) = cons_info%const_dist_mol(k)
               cons_info%const_dist_mol(k) = dum
               !
               dum = cons_info%const_dist_a(j)
               cons_info%const_dist_a(j) = cons_info%const_dist_a(k)
               cons_info%const_dist_a(k) = dum
               !
               dum = cons_info%const_dist_b(j)
               cons_info%const_dist_b(j) = cons_info%const_dist_b(k)
               cons_info%const_dist_b(k) = dum
               !
               xdum = cons_info%const_dist_dab(j)
               cons_info%const_dist_dab(j) = cons_info%const_dist_dab(k)
               cons_info%const_dist_dab(k) = xdum
            END IF
          END DO
        END DO
        i=last+1
      ELSE
        i=i+1
      END IF
    END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Set the distance constraint molecule_kind_set(ikind)%dist_list
  !-----------------------------------------------------------------------------
  IF(topology%const_dist) THEN
    DO i=1,SIZE(molecule_kind_set)
      molecule_kind => molecule_kind_set(i)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             nmolecule=nmolecule,molecule_list=molecule_list)
      CALL find_boundary(cons_info%const_dist_mol,cons_info%nconst_dist,&
                         first,last,i)
      ndist = last-first+1
      ALLOCATE(dist_list(ndist),STAT=istat)
      IF (istat /= 0) CALL stop_memory('topology_constraint_pack','dist_list')
      DO j=1,nmolecule
        molecule => molecule_set(molecule_list(j))
        ALLOCATE(ldist(ndist),STAT=istat)
        IF (istat /= 0) CALL stop_memory('topology_constraint_pack','ldist')
        DO k=1,ndist
          ldist(k)%del_lam = 0.0_dp
          ldist(k)%lambda = 0.0_dp
          ldist(k)%scale = 0.0_dp
          ldist(k)%scale_old = 0.0_dp
          ldist(k)%sigma = 0.0_dp
          ldist(k)%f(:) = 0.0_dp
          ldist(k)%ra_old(:) = 0.0_dp
          ldist(k)%rb_old(:) = 0.0_dp
        END DO
        CALL set_molecule(molecule=molecule,ldist=ldist)
      END DO
      ndist = 0
      DO j=first,last
        ndist = ndist + 1
        dist_list(ndist)%a   = cons_info%const_dist_a(j)
        dist_list(ndist)%b   = cons_info%const_dist_b(j)
        dist_list(ndist)%dab = cons_info%const_dist_dab(j)
      END DO
      CALL set_molecule_kind(molecule_kind,dist_list=dist_list,ndist=ndist)
    END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Set the group 3x3 constraint g3x3_list
  !-----------------------------------------------------------------------------
  IF(topology%const_33) THEN
    DO i=1,cons_info%nconst_g33
      j = cons_info%const_g33_mol(i)
      molecule_kind => molecule_kind_set(j)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             nmolecule=nmolecule,molecule_list=molecule_list)
      CALL find_boundary(cons_info%const_g33_mol,cons_info%nconst_g33,&
                         first,last,i)
      ng3x3 = last-first+1
      ALLOCATE(g3x3_list(ng3x3),STAT=istat)
      IF (istat /= 0) CALL stop_memory('topology_constraint_pack','g3x3_list')
      DO j=1,nmolecule
        molecule => molecule_set(molecule_list(j))
        ALLOCATE(lg3x3(ng3x3),STAT=istat)
        IF (istat /= 0) CALL stop_memory('topology_constraint_pack','lg3x3')
        DO k=1,ng3x3
          lg3x3(k)%scale = 0.0_dp
          lg3x3(k)%scale_old = 0.0_dp
          lg3x3(k)%fa = 0.0_dp
          lg3x3(k)%fb = 0.0_dp
          lg3x3(k)%fc = 0.0_dp
          lg3x3(k)%ra_old = 0.0_dp
          lg3x3(k)%rb_old = 0.0_dp
          lg3x3(k)%rc_old = 0.0_dp
          lg3x3(k)%va = 0.0_dp
          lg3x3(k)%vb = 0.0_dp
          lg3x3(k)%vc = 0.0_dp
          lg3x3(k)%va_old = 0.0_dp
          lg3x3(k)%vb_old = 0.0_dp
          lg3x3(k)%vc_old = 0.0_dp
          lg3x3(k)%del_lam = 0.0_dp
          lg3x3(k)%lambda = 0.0_dp
        END DO
        CALL set_molecule(molecule=molecule,lg3x3=lg3x3)
      END DO
      ng3x3 = 0
      DO j=1,SIZE(cons_info%const_g33_a)
        ng3x3 = ng3x3 + 1
        g3x3_list(ng3x3)%a   = cons_info%const_g33_a(j)
        g3x3_list(ng3x3)%b   = cons_info%const_g33_b(j)
        g3x3_list(ng3x3)%c   = cons_info%const_g33_c(j)
        g3x3_list(ng3x3)%dab = cons_info%const_g33_dab(j)
        g3x3_list(ng3x3)%dac = cons_info%const_g33_dac(j)
        g3x3_list(ng3x3)%dbc = cons_info%const_g33_dbc(j)
      END DO
      CALL set_molecule_kind(molecule_kind,ng3x3=ng3x3,g3x3_list=g3x3_list)
    END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Set the group 4x6 constraint g4x6_list
  !-----------------------------------------------------------------------------
  IF(topology%const_46) THEN
    DO i=1,cons_info%nconst_g46
      j = cons_info%const_g46_mol(i)
      molecule_kind => molecule_kind_set(j)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             nmolecule=nmolecule,molecule_list=molecule_list)
      CALL find_boundary(cons_info%const_g46_mol,cons_info%nconst_g46,&
                         first,last,i)
      ng4x6 = last-first+1
      ALLOCATE(g4x6_list(ng4x6),STAT=istat)
      IF (istat /= 0) CALL stop_memory('topology_constraint_pack','g4x6_list')
      DO j=1,nmolecule
        molecule => molecule_set(molecule_list(j))
        ALLOCATE(lg4x6(ng4x6),STAT=istat)
        IF (istat /= 0) CALL stop_memory('topology_constraint_pack','lg4x6')
        DO k=1,ng4x6
          lg4x6(k)%scale = 0.0_dp
          lg4x6(k)%scale_old = 0.0_dp
          lg4x6(k)%fa = 0.0_dp
          lg4x6(k)%fb = 0.0_dp
          lg4x6(k)%fc = 0.0_dp
          lg4x6(k)%fd = 0.0_dp
          lg4x6(k)%fe = 0.0_dp
          lg4x6(k)%ff = 0.0_dp
          lg4x6(k)%ra_old = 0.0_dp
          lg4x6(k)%rb_old = 0.0_dp
          lg4x6(k)%rc_old = 0.0_dp
          lg4x6(k)%rd_old = 0.0_dp
          lg4x6(k)%re_old = 0.0_dp
          lg4x6(k)%rf_old = 0.0_dp
          lg4x6(k)%va = 0.0_dp
          lg4x6(k)%vb = 0.0_dp
          lg4x6(k)%vc = 0.0_dp
          lg4x6(k)%vd = 0.0_dp
          lg4x6(k)%ve = 0.0_dp
          lg4x6(k)%vf = 0.0_dp
          lg4x6(k)%va_old = 0.0_dp
          lg4x6(k)%vb_old = 0.0_dp
          lg4x6(k)%vc_old = 0.0_dp
          lg4x6(k)%vd_old = 0.0_dp
          lg4x6(k)%ve_old = 0.0_dp
          lg4x6(k)%vf_old = 0.0_dp
          lg4x6(k)%del_lam = 0.0_dp
          lg4x6(k)%lambda = 0.0_dp
        END DO
        CALL set_molecule(molecule=molecule,lg4x6=lg4x6)
      END DO
      ng4x6 = 0
      DO j=1,SIZE(cons_info%const_g46_a)
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
      END DO
      CALL set_molecule_kind(molecule_kind,ng4x6=ng4x6,g4x6_list=g4x6_list)
    END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Set the group fixed_atom constraint fixd_list
  !-----------------------------------------------------------------------------
  IF(topology%const_atom) THEN
     DO i=1,SIZE(molecule_kind_set)
        molecule_kind => molecule_kind_set(i)
        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               nmolecule=nmolecule,molecule_list=molecule_list)
        ! Try to figure out how many atoms of the list belong to this molecule_kind
        nfixed_atoms = 0
        DO j=1,nmolecule
           molecule => molecule_set(molecule_list(j))
           CALL get_molecule ( molecule, first_atom = first, last_atom = last )
           DO k = first, last
              IF (PRESENT(qmmm_env)) THEN
                 IF  (ANY(cons_info%fixed_atoms == k).OR.&
                      (cons_info%freeze_qm.AND.ANY(qmmm_env%qm_molecule_index == molecule_list(j))).OR.&
                      (cons_info%freeze_mm.AND.ALL(qmmm_env%qm_molecule_index /= molecule_list(j))))&
                      nfixed_atoms = nfixed_atoms + 1 
              ELSE
                 IF  (ANY(cons_info%fixed_atoms == k))&
                      nfixed_atoms = nfixed_atoms + 1 
              END IF
           END DO
        END DO
        ALLOCATE(fixd_list(nfixed_atoms),STAT=istat)
        IF (istat /= 0) CALL stop_memory('topology_constraint_pack','fixd_list')
        kk = 0
        DO j=1,nmolecule
           molecule => molecule_set(molecule_list(j))
           CALL get_molecule ( molecule, first_atom = first, last_atom = last )
           DO k = first, last
              IF (PRESENT(qmmm_env)) THEN
                 IF  (ANY(cons_info%fixed_atoms == k).OR.&
                      (cons_info%freeze_qm.AND.ANY(qmmm_env%qm_molecule_index == molecule_list(j))).OR.&
                      (cons_info%freeze_mm.AND.ALL(qmmm_env%qm_molecule_index /= molecule_list(j)))) THEN
                    kk = kk + 1
                    fixd_list(kk)%fixd = k
                    particle_set(k)%is_fixed = .TRUE.
                 END IF
              ELSE
                 IF  (ANY(cons_info%fixed_atoms == k)) THEN                     
                    kk = kk + 1
                    fixd_list(kk)%fixd = k
                    particle_set(k)%is_fixed = .TRUE.
                 END IF
              END IF
           END DO
        END DO
        IF (PM1) WRITE(*,*)"FIXD_LIST ::",fixd_list
        CALL set_molecule_kind(molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list)
     END DO
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_constraint_pack"
  CALL write_checkpoint_information("leaving "//routineN,globenv)

END SUBROUTINE topology_constraint_pack


END MODULE topology_util
