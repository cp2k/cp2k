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
  USE checkpoint_handler,              ONLY: Pwarning
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: ptable
  USE string_utilities,                ONLY: integer_to_string,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_util'

  TYPE atom_bond_list_type
     INTEGER, POINTER, DIMENSION(:) :: bonds => NULL()
  END TYPE atom_bond_list_type

  PRIVATE
  PUBLIC :: topology_set_atm_mass,&
            topology_reorder_atoms,&
            topology_molecules_check,&
            check_subsys_element,&
            reorder_structure,&
            find_molecule,&
            atom_bond_list_type
  
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
  SUBROUTINE topology_reorder_atoms(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
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
!!    reorder_structure 
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     each atom will be assigned a molecule number based on bonded fragments
!!
!!   AUTHOR
!!     Joost 05.2006
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reorder_structure(work, list1, list2, N, error)
    TYPE(atom_bond_list_type), &
      DIMENSION(:), INTENT(INOUT)            :: work
    INTEGER, DIMENSION(:), INTENT(IN)        :: list1, list2
    INTEGER, INTENT(IN)                      :: N
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorder_structure', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, I, Nbond, stat
    INTEGER, DIMENSION(:), POINTER           :: wrk_tmp
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO I=1,N
       atom_a=list1(I)
       atom_b=list2(I)

       wrk_tmp => work(atom_a)%bonds
       Nbond=SIZE(wrk_tmp)
       ALLOCATE(work(atom_a)%bonds(Nbond+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(atom_a)%bonds(1:Nbond) = wrk_tmp
       work(atom_a)%bonds(Nbond+1) = atom_b
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       wrk_tmp => work(atom_b)%bonds
       Nbond=SIZE(wrk_tmp)
       ALLOCATE(work(atom_b)%bonds(Nbond+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(atom_b)%bonds(1:Nbond) = wrk_tmp
       work(atom_b)%bonds(Nbond+1) = atom_a
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO    

  END SUBROUTINE reorder_structure
    
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_molecule
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     each atom will be assigned a molecule number based on bonded fragments
!!
!!   AUTHOR
!!     Joost 05.2006
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE find_molecule(atom_bond_list,mol_info)
    TYPE(atom_bond_list_type), &
      DIMENSION(:), INTENT(IN)               :: atom_bond_list
    INTEGER, DIMENSION(:), POINTER           :: mol_info

    INTEGER                                  :: I, N, nmol

    N=SIZE(atom_bond_list)
    nmol=0
    DO I=1,N
       IF (mol_info(I)==-1) THEN
          nmol=nmol+1
          CALL spread_mol(atom_bond_list,mol_info,i,nmol)
       ENDIF
    ENDDO
  END SUBROUTINE find_molecule

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     spread_mol
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     spreads the molnumber over the bonded list
!!
!!   AUTHOR
!!     Joost 05.2006
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  RECURSIVE SUBROUTINE spread_mol(atom_bond_list,mol_info,iatom,imol)
    TYPE(atom_bond_list_type), &
      DIMENSION(:), INTENT(IN)               :: atom_bond_list
    INTEGER, DIMENSION(:), POINTER           :: mol_info
    INTEGER, INTENT(IN)                      :: iatom, imol

    INTEGER                                  :: atom_b, i

    mol_info(iatom)=imol
    DO I=1,SIZE(atom_bond_list(iatom)%bonds)
       atom_b=atom_bond_list(iatom)%bonds(I)
       IF (mol_info(atom_b)==-1) &
            CALL spread_mol(atom_bond_list,mol_info,atom_b,imol)
       IF (mol_info(atom_b)/=imol) STOP "internal error"
    ENDDO
  END SUBROUTINE spread_mol
 
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
  SUBROUTINE topology_set_atm_mass(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
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
  SUBROUTINE topology_molecules_check(topology,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_molecules_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, first, first_loc, &
                                                handle, i, iatom, iw, k, &
                                                loc_counter, mol_num, &
                                                mol_typ, n, natom, stat
    LOGICAL                                  :: failure, icheck_num, &
                                                icheck_typ
    TYPE(atom_bond_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom_bond_list
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

    IF (iw>0) WRITE(iw,'(A)')"Start of Molecule_Check",&
         "  Checking consistency between the generated molecules"

    ALLOCATE(atom_bond_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       ALLOCATE(atom_bond_list(I)%bonds(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%bond_a)) N = SIZE(conn_info%bond_a)
    CALL reorder_structure(atom_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)

    mol_typ = atom_info%map_mol_typ(1)
    mol_num = atom_info%map_mol_num(1)
    counter     = 1
    loc_counter = 1
    first       = 1
    first_loc   = 1
    DO iatom = 2, natom
       icheck_num = (atom_info%map_mol_num(iatom) == mol_num)
       icheck_typ = (atom_info%map_mol_typ(iatom) == mol_typ)
       IF ((icheck_typ.AND.(.NOT.icheck_num)).OR.(.NOT.icheck_typ)) THEN
          !-----------------------------------------------------------------------------
          !-----------------------------------------------------------------------------
          ! 1. Check each molecule have the same number of atoms
          !-----------------------------------------------------------------------------
          IF (counter /= loc_counter) THEN
             CALL stop_program ("topology_molecules_check",&
                  "different number of atoms for same molecule kind"//&
                  " molecule type  = "//cp_to_string(mol_typ)//&
                  " molecule number= "//cp_to_string(mol_num)//&
                  " expected number of atoms="//cp_to_string(counter)//" found="//&
                  cp_to_string(loc_counter))            
          END IF
       END IF
       IF (.NOT.icheck_typ) THEN
          first       = iatom
          first_loc   = iatom
          counter     = 1
          loc_counter = 1
          mol_typ     = atom_info%map_mol_typ(iatom)
       END IF
       IF (icheck_num) THEN 
          IF (icheck_typ) loc_counter = loc_counter + 1
          !-----------------------------------------------------------------------------
          !-----------------------------------------------------------------------------
          ! 2. Check that each molecule has the same atom sequences
          !-----------------------------------------------------------------------------
          IF( TRIM(atom_info%label_atmname(iatom)) /= &
              TRIM(atom_info%label_atmname(first+loc_counter-1)) ) THEN
             CALL stop_program ("topology_molecules_check",&
                  "different atom name for same molecule kind"//&
                  " molecule type  = "//cp_to_string(mol_typ)//&
                  " molecule number= "//cp_to_string(mol_num)//&
                  " expected atom name="//TRIM(atom_info%label_atmname(first+loc_counter-1))//&
                  " found="//TRIM(atom_info%label_atmname(iatom)))
          END IF
          !-----------------------------------------------------------------------------
          !-----------------------------------------------------------------------------
          ! 3. Check that each molecule have the same bond sequences
          !-----------------------------------------------------------------------------
          IF (SIZE(atom_bond_list(iatom)%bonds)/=SIZE(atom_bond_list(iatom)%bonds)) failure = .TRUE.
          IF (.NOT.failure) THEN
             DO k = 1,SIZE(atom_bond_list(iatom)%bonds)
                IF (ALL(atom_bond_list(first+loc_counter-1)%bonds-first/=&
                     atom_bond_list(iatom)%bonds(k)-first_loc)) THEN
                   CALL stop_program ("topology_molecules_check",&
                        "different bonds for same molecule kind"//&
                        " molecule type  = "//cp_to_string(mol_typ)//&
                        " molecule number= "//cp_to_string(mol_num)//&
                        " expected bonds="//cp_to_string(iatom-first_loc+1)//" - "//&
                        cp_to_string(atom_bond_list(iatom)%bonds(k)-first_loc)//&
                        " NOT FOUND! Check the connectivity of your system.")
                END IF
             END DO
          END IF
       ELSE
          mol_num = atom_info%map_mol_num(iatom)
          loc_counter = 1
          first_loc   = iatom
       END IF
       IF (mol_num==1.AND.icheck_typ) counter = counter + 1
    END DO
    IF (iw>0) WRITE(iw,'(A)')"End of Molecule_Check"

    DO I=1,natom
       DEALLOCATE(atom_bond_list(I)%bonds,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(atom_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
  SUBROUTINE check_subsys_element(element_in, element_out, root_section, force_env_section, error)
    CHARACTER(len=*), INTENT(IN)             :: element_in
    CHARACTER(len=default_string_length), &
      INTENT(OUT)                            :: element_out
    TYPE(section_vals_type), POINTER         :: root_section, &
                                                force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'check_subsys_element', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: current_symbol, &
                                                element_symbol, keyword
    INTEGER                                  :: i, i_rep, n_rep
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
          CALL Pwarning("No element specification provided for type:"//element_in//".",root_section, error)
          iwarn = iwarn - 1
          IF (iwarn == 0) THEN
             CALL Pwarning("Further Output regarding unknown element types will be suppressed!",&
                  root_section,error)
          END IF
       END IF
    END IF

  END SUBROUTINE check_subsys_element

END MODULE topology_util
