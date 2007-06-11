!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/topology_generate_util [1.0] *
!!
!!   NAME
!!     topology_generate_util
!!
!!   FUNCTION
!!     Collection of subroutine needed for topology related things
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!    Teodor Laino 09.2006 - Major rewriting with linear scaling routines
!!
!!   SOURCE
!******************************************************************************

MODULE topology_generate_util
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             deallocate_atomic_kind_set,&
                                             set_atomic_kind
  USE cell_types,                      ONLY: pbc
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_deallocate,&
                                             fist_neighbor_type
  USE fist_neighbor_lists,             ONLY: build_fist_neighbor_lists
  USE input_constants,                 ONLY: do_add,&
                                             do_bondparm_covalent,&
                                             do_bondparm_vdw,&
                                             do_conn_off,&
                                             do_remove
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             deallocate_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE string_utilities,                ONLY: integer_to_string,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: array1_list_type,&
                                             array2_list_type,&
                                             find_molecule,&
                                             give_back_molecule,&
                                             reorder_list_array,&
                                             reorder_structure
  USE util,                            ONLY: find_boundary,&
                                             sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_generate_util'

  PRIVATE
  PUBLIC :: topology_generate_bend,&
            topology_generate_bond,&
            topology_generate_dihe,&
            topology_generate_impr,&
            topology_generate_onfo,&
            topology_generate_ub,&
            topology_generate_molecule
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_generate_util [1.0] *
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
  SUBROUTINE topology_generate_molecule(topology,qmmm,qmmm_env,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_in_kind, atom_in_mol, first, handle, handle2, i, iatm, &
      iatom, iend, ifirst, ilast, inum, istart, iw, j, jump1, jump2, last, &
      max_mol_num, mol_num, mol_res, mol_typ, myind, N, natom, nlocl, ntype, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, wrk1, wrk2
    LOGICAL                                  :: do_again, failure, my_qmmm
    TYPE(array1_list_type), ALLOCATABLE, &
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
    NULLIFY(qm_atom_index)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    !
    ! QM/MM coordinate_control
    ! 
    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm
    
    natom = topology%natoms
    stat  = 0
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
    atom_info%map_mol_typ(:) =  0
    atom_info%map_mol_num(:) = -1
    atom_info%map_mol_res(:) =  1
    
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
    IF(iw>0) WRITE(iw,'(/,A)') "Start of molecule generation"

    ! convert a simple list of bonds to a list of bonds per atom
    ! (each bond is present in the forward and backward direction
    ALLOCATE(atom_bond_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       ALLOCATE(atom_bond_list(I)%array1(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%bond_a)) N = SIZE(conn_info%bond_a)
    CALL reorder_structure(atom_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)
    CALL find_molecule(atom_bond_list,atom_info%map_mol_num, atom_info%label_molname)   
    DO I=1,natom
       DEALLOCATE(atom_bond_list(I)%array1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(atom_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    ! Modify according map_mol_typ the array map_mol_num
    !
    IF(iw>0) WRITE(iw,'(/,A)') "End of molecule generation"
    IF(iw>0) WRITE(iw,'(/,A)') "Checking for non-continuous generated molecules"
    ! Check molecule number
    ALLOCATE(wrk1(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wrk2(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    wrk1 = atom_info%map_mol_num
    CALL sort(wrk1,natom,wrk2)
    istart  = 1
    mol_typ = wrk1(istart)
    DO i = 2, natom
       IF (mol_typ /= wrk1(i)) THEN
          iend  = i-1
          first = MINVAL(wrk2(istart:iend))
          last  = MAXVAL(wrk2(istart:iend))
          nlocl = last-first+1
          IF (iend-istart+1 /= nlocl) THEN
             CALL stop_program ("topology_generate_molecule",&
                  "Detected a non contiguous molecule!! "//&
                  "Too late at this stage everything should be already ordered! "//&
                  "Try with the keyword reordering.. It may help!")
          END IF
          istart = i
       END IF
    END DO
    iend  = i-1
    first = MINVAL(wrk2(istart:iend))
    last  = MAXVAL(wrk2(istart:iend))
    nlocl = last-first+1
    IF (iend-istart+1 /= nlocl) THEN
       CALL stop_program ("topology_generate_molecule",&
            "Detected a non contiguous molecule!! "//&
            "Too late at this stage everything should be already ordered! "//&
            "Try with the keyword reordering.. It may help!")
    END IF
    DEALLOCATE(wrk1, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wrk2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
    IF(iw>0) WRITE(iw,'(/,A)') "End of check"
    IF(iw>0) WRITE(iw,'(/,A)') "Start of renumbering molecules"
    mol_typ = atom_info%map_mol_typ(1)
    mol_num = atom_info%map_mol_num(1)
    DO i = 2, natom
       IF (atom_info%map_mol_typ(i)/=mol_typ) THEN
          myind   = atom_info%map_mol_num(i) - mol_num + 1
          CPPostcondition(myind/=atom_info%map_mol_num(i-1),cp_failure_level,routineP,error,failure)
          mol_typ = atom_info%map_mol_typ(i)
          mol_num = atom_info%map_mol_num(i)
       END IF
       atom_info%map_mol_num(i) = atom_info%map_mol_num(i) - mol_num + 1
    END DO
    IF(iw>0) WRITE(iw,'(/,A)') "End of renumbering molecules"
    !
    ! Residues
    !
    CALL timeset(routineN//"_PARA_RES",'I','',handle2)
    IF(iw>0) WRITE(iw,'(/,A,L2)') "Starting PARA_RES: ",topology%para_res
    IF(topology%para_res) THEN
       mol_res = 1
       mol_typ = atom_info%map_mol_typ(1)
       mol_num = atom_info%map_mol_num(1)
       atom_info%map_mol_res(1) = mol_res
       DO i = 2, natom
          IF((atom_info%label_resid(i-1)/=atom_info%label_resid(i)).OR.&
             (TRIM(atom_info%label_resname(i-1))/=TRIM(atom_info%label_resname(i))))THEN
             mol_res = mol_res + 1
          END IF
          IF (atom_info%map_mol_typ(i)/=mol_typ.OR.atom_info%map_mol_num(i)/=mol_num) THEN
             mol_typ = atom_info%map_mol_typ(i)
             mol_num = atom_info%map_mol_num(i)
             mol_res = 1
          END IF
          atom_info%map_mol_res(i) = mol_res
       END DO
    END IF
    IF(iw>0) WRITE(iw,'(/,A)') "End of PARA_RES"
    CALL timestop(0.0_dp,handle2)
    IF(iw>0) THEN
       DO iatom=1,natom
          WRITE(iw,'(4(1X,A,":",I0),2(1X,A,1X,A))') "iatom",iatom,&
                                       "map_mol_typ",atom_info%map_mol_typ(iatom),&
                                       "map_mol_num",atom_info%map_mol_num(iatom),&
                                       "map_mol_res",atom_info%map_mol_res(iatom),&
                                       "mol_name:",TRIM(atom_info%label_molname(iatom)),&
                                       "res_name:",TRIM(atom_info%label_resname(iatom))
       END DO
    END IF
    IF (my_qmmm) THEN
       do_again = .FALSE.
       IF (iw>0) WRITE(iw,*)"MAP_MOL_NUM ",atom_info%map_mol_num
       IF (iw>0) WRITE(iw,*)"MAP_MOL_TYP ",atom_info%map_mol_typ
       IF (iw>0) WRITE(iw,*)"MAP_MOL_RES ",atom_info%map_mol_res
       ALLOCATE(qm_atom_index(SIZE(qmmm_env%qm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       qm_atom_index = qmmm_env%qm_atom_index
       CPPostcondition(ALL(qm_atom_index /= 0),cp_failure_level,routineP,error,failure)
       DO myind = 1, SIZE(qm_atom_index)
          IF (qm_atom_index(myind) ==0) CYCLE
          CALL find_boundary(atom_info%map_mol_typ,natom,ifirst,ilast,&
               atom_info%map_mol_typ(qm_atom_index(myind)))
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
    IF (iw>0) THEN 
       WRITE(iw,*)"SUMMARY:: Number of molecules type found:",ntype
       ntype = MAXVAL(atom_info%map_mol_typ)
       DO i = 1, ntype
          atom_in_kind   = COUNT(atom_info%map_mol_typ == i)
          WRITE(iw,*)"Molecule type: ",i,"contains ",atom_in_kind,' atoms'
          IF (atom_in_kind <= 1) CYCLE
          CALL find_boundary(atom_info%map_mol_typ,natom,first,last,i)
          WRITE(iw,*)"Boundary atoms: ",first, last
          CPPostcondition(last-first+1==atom_in_kind,cp_failure_level,routineP,error,failure)
          max_mol_num    = MAXVAL(atom_info%map_mol_num(first:last))
          WRITE(iw,*)"Number of molecules of type ",i," is ::",max_mol_num       
          atom_in_mol    = atom_in_kind / max_mol_num
          WRITE(iw,*)"Number of atoms per each molecule: ",atom_in_mol
          WRITE(iw,*)"MAP_MOL_TYP:: ",atom_info%map_mol_typ(first:last)
          WRITE(iw,*)"MAP_MOL_NUM:: ",atom_info%map_mol_num(first:last)
          WRITE(iw,*)"MAP_MOL_RES:: ",atom_info%map_mol_res(first:last)
          !
          DO j = 1, max_mol_num
             IF (COUNT(atom_info%map_mol_num(first:last)==j) /= atom_in_mol) THEN
                WRITE(iw,*)"molecule type:",i,"molecule num:",j," has ",&
                     COUNT(atom_info%map_mol_num(first:last)==j),&
                     " atoms instead of ",atom_in_mol," ."
                CALL stop_program ("topology_generate_molecule","Two molecules of the same kind"//&
                     " have been created with different numbers of atoms!")
             END IF
          END DO
       END DO
    END IF
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE topology_generate_molecule
  
!******************************************************************************
!!****** topology_generate_util [1.0] *
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
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_bond(topology,para_env,subsys_section,&
       force_env_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section, &
                                                force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bond', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: upper_sym_1, upper_sym_2
    CHARACTER(LEN=2), POINTER                :: atm_symbol(:)
    INTEGER :: cbond, handle, handle2, i, iatm1, iatm2, iatom, ibond, idim, &
      ielem, ielem_found, iw, j, jatom, k, n_bonds, n_heavy_bonds, &
      n_hydr_bonds, natom, nelem, npairs, output_unit, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bond_a, bond_b, list, map_nb
    LOGICAL                                  :: connectivity_ok, &
                                                element_found, failure, &
                                                ionode, print_info
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: h_list
    REAL(KIND=dp)                            :: bondparm_factor, cell_v(3), &
                                                dr(3), ksign, my_maxrad, r2, &
                                                r2_min, rbond, rbond2, tmp
    REAL(KIND=dp), DIMENSION(1, 1)           :: r_max, r_minsq
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pbc_coord
    TYPE(array2_list_type), DIMENSION(:), &
      POINTER                                :: bond_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: bond_section

    failure = .FALSE.    
    NULLIFY(logger, particle_set, atomic_kind_set, nonbonded)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,'I','',handle)
    IF (.NOT.failure) THEN
       atom_info => topology%atom_info
       conn_info => topology%conn_info
       bondparm_factor = topology%bondparm_factor
       cbond = 0    
       natom = topology%natoms    
       NULLIFY(atm_symbol,radius)
       ! Allocate temporary arrays
       ALLOCATE(atm_symbol(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(radius(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(list(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(h_list(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(pbc_coord(3,natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       h_list = .FALSE.
       CALL timeset(TRIM(routineN)//"_1",'I','',handle2)
       DO iatom=1,natom
          atm_symbol(iatom) = ''
          IF(ASSOCIATED(atom_info%element)) THEN
             atm_symbol(iatom) = atom_info%element(iatom)
          END IF
       END DO
       nelem = UBOUND(ptable,1)
       DO iatom=1,natom
          list(iatom) = iatom
          element_found=.FALSE.
          ielem_found  =0
          upper_sym_1  =atm_symbol(iatom)
          CALL uppercase(upper_sym_1)
          DO ielem=1,nelem 
             upper_sym_2=ptable(ielem)%symbol
             CALL uppercase(upper_sym_2)
             IF(upper_sym_1 == upper_sym_2)THEN
                element_found=.TRUE.
                ielem_found=ielem
                EXIT
             END IF
          END DO
          IF (upper_sym_2=="H ") h_list(iatom) = .TRUE.
          IF (.NOT. element_found) THEN
             failure = .TRUE.
             IF (output_unit>0) WRITE(output_unit,'(A)') "ERROR : topology_generate_bond: The element " &
                  //atm_symbol(iatom)//" was not found in the periodic table.",&
                  " Specify a valid element in order to use generate topology"
          ENDIF
          
          IF(topology%bondparm_type==do_bondparm_covalent) THEN
             radius(iatom) = ptable(ielem_found)%covalent_radius
          ELSE IF(topology%bondparm_type==do_bondparm_vdw) THEN
             radius(iatom) = ptable(ielem_found)%vdw_radius
          ELSE
             CALL stop_program ("topology_generate_bond","illegal bondparm_type")
          END IF
          CALL convert_to_cp2k_units("ANGSTROM",length=radius(iatom))
          IF(iw>0) WRITE(iw,'(T2,"GENERATE|",5X,A,T77,A5)') &
               "In topology_generate_bond :: iatom = ",ptable(ielem_found)%symbol
       END DO
       CPPostcondition(.NOT.failure,cp_failure_level,routineP,error,failure)
       CALL timestop(0.0_dp,handle2)
       CALL timeset(TRIM(routineN)//"_2",'I','',handle2)
       ! Initialize fake particle_set and atomic_kinds to generate the bond list
       ! using the neighboring list routine
       CALL allocate_atomic_kind_set(atomic_kind_set,1,error)
       CALL allocate_particle_set(particle_set,natom)
       !
       my_maxrad = MAXVAL(radius)*2.0_dp
       atomic_kind => atomic_kind_set(1)
       CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=1,&
            name="XXX",element_symbol="XXX",mass=0.0_dp,atom_list=list)
       CALL section_vals_val_get(subsys_section,"TOPOLOGY%GENERATE%BONDLENGTH_MAX",r_val=tmp,&
            error=error)
       r_max     = tmp
       IF (my_maxrad*bondparm_factor>r_max(1,1).AND.(.NOT.topology%molname_generated)) THEN
          IF (output_unit>0) THEN
             WRITE(output_unit,'(T2,"GENERATE|",A)')&
                  " ERROR in connectivity generation!",&
                  " The THRESHOLD to select possible bonds is bigger than the MAX bondlength",&
                  " used to build the neighbors lists. Increase the BONDLENGTH_MAX patameter"
             WRITE(output_unit,'(T2,"GENERATE|",2(A,F11.6),A)')&
                  " Present THRESHOLD (",my_maxrad*bondparm_factor," ). "//&
                  " Present BONDLENGTH_MAX (",r_max(1,1)," )"
          END IF
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       DO i=1,natom
          particle_set(i)%atomic_kind => atomic_kind_set(1)
          particle_set(i)%r(1) = atom_info%r(1,i)
          particle_set(i)%r(2) = atom_info%r(2,i)
          particle_set(i)%r(3) = atom_info%r(3,i)
          pbc_coord(:,i) = pbc(atom_info%r(:,i),topology%cell)
          NULLIFY(particle_set(i)%list_exclude_vdw)
          NULLIFY(particle_set(i)%list_exclude_ei)
          ALLOCATE(particle_set(i)%list_exclude_vdw(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(particle_set(i)%list_exclude_ei(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
       END DO
       CALL section_vals_val_get(subsys_section,"TOPOLOGY%GENERATE%BONDLENGTH_MIN",r_val=tmp,&
            error=error)
       r_minsq = tmp*tmp
       CALL timestop(0.0_dp,handle2)
       CALL timeset(TRIM(routineN)//"_3",'I','',handle2)
       CALL build_fist_neighbor_lists(atomic_kind_set,particle_set,&
            cell=topology%cell, r_max=r_max, r_minsq=r_minsq, nonbonded=nonbonded,&
            para_env=para_env, build_from_scratch=.TRUE.,&
            force_env_section=force_env_section,error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T2,"GENERATE| Number of prescreened bonds (neighbors):",T71,I10)')&
               nonbonded%neighbor_kind_pairs(1)%npairs
       END IF
       npairs = 0
       DO i = 1, SIZE(nonbonded%neighbor_kind_pairs)
          npairs = npairs + nonbonded%neighbor_kind_pairs(i)%npairs
       END DO
       ALLOCATE(bond_a(npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(bond_b(npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(map_nb(npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       idim = 0
       DO j = 1, SIZE(nonbonded%neighbor_kind_pairs)
          DO i = 1, nonbonded%neighbor_kind_pairs(j)%npairs
             idim = idim + 1
             bond_a(idim) = nonbonded%neighbor_kind_pairs(j)%list(1,i)
             bond_b(idim) = nonbonded%neighbor_kind_pairs(j)%list(2,i)
             map_nb(idim) = j
          END DO
       END DO
       CALL timestop(0.0_dp,handle2)
       CALL timeset(TRIM(routineN)//"_4",'I','',handle2)
       ! We have a list of neighbors let's order the list w.r.t. the particle number
       ALLOCATE(bond_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I = 1, natom
          ALLOCATE(bond_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(bond_list(I)%array2(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       CALL reorder_structure(bond_list, bond_a, bond_b, map_nb, SIZE(bond_a),&
            error=error)
       DEALLOCATE(bond_a,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(bond_b,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(map_nb,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Find the Real bonds in the system
       ! Let's start with heavy atoms.. hydrogens will be treated only later on...
       ! Heavy atoms loop
       CALL reallocate(conn_info%bond_a,1,1)
       CALL reallocate(conn_info%bond_b,1,1)
       connectivity_ok =.FALSE.
       ! No need to check consistency between provided molecule name and
       ! generated connectivity since we overrided the molecule definition.
       IF (topology%create_molecules) THEN
          atom_info%label_molname = "TO_DEFINE_LATER"
          ! A real name assignment will then be performed in the reorder module..
       END IF
       ! It may happen that the connectivity created is fault for the missing
       ! of one bond.. this external loop ensures that everything was created
       ! fits exactly with the definition of molecules..
       DO WHILE (.NOT.connectivity_ok)
          n_heavy_bonds = 0
          n_bonds       = 0
          DO iatm1 = 1, natom
             IF (h_list(iatm1)) CYCLE
             DO j = 1,SIZE(bond_list(iatm1)%array1)
                iatm2 = bond_list(iatm1)%array1(j)
                IF (atom_info%label_molname(iatm1)/=atom_info%label_molname(iatm2)) CYCLE
                IF(h_list(iatm2).OR.(iatm2<=iatm1)) CYCLE
                k     = bond_list(iatm1)%array2(j)
                ksign = SIGN(1.0_dp,REAL(k,KIND=dp))
                k     = ABS(k)
                CALL matvec_3x3(cell_v,topology%cell%hmat,&
                     REAL(nonbonded%neighbor_kind_pairs(k)%cell_vector,KIND=dp))
                dr = pbc_coord(:,iatm1) - pbc_coord(:,iatm2) - ksign*cell_v
                r2 = DOT_PRODUCT(dr,dr)
                IF(r2 <= 0.01_dp) THEN
                   CALL stop_program ("topology_generate_bond",&
                        "bond distance between atoms less then 0.01")
                END IF
                
                ! Screen neighbors
                IF(topology%bondparm_type==do_bondparm_covalent) THEN
                   rbond = radius(iatm1) + radius(iatm2)
                ELSE IF(topology%bondparm_type==do_bondparm_vdw) THEN
                   rbond = MAX(radius(iatm1),radius(iatm2))
                END IF
                rbond2 = rbond*rbond
                rbond2 = rbond2*(bondparm_factor)**2
                !Test the distance to the sum of the covalent radius
                IF(r2 <= rbond2) THEN
                   n_heavy_bonds = n_heavy_bonds + 1
                   CALL add_bonds_list(conn_info, iatm1, iatm2, n_heavy_bonds)
                END IF
             END DO
          END DO
          n_hydr_bonds  = 0
          n_bonds       = n_heavy_bonds
          ! Now check bonds formed by hydrogens...
          ! The hydrogen valence is 1 so we can choose the closest atom..
          DO iatm1 = 1, natom
             IF (.NOT.h_list(iatm1)) CYCLE
             r2_min = HUGE(0.0_dp)
             ibond  = -1
             print_info = .TRUE.
             DO j = 1,SIZE(bond_list(iatm1)%array1)
                iatm2 = bond_list(iatm1)%array1(j)
                print_info = .FALSE.
                IF (atom_info%label_molname(iatm1)/=atom_info%label_molname(iatm2)) CYCLE
                IF (h_list(iatm2).AND.(iatm2<=iatm1)) CYCLE
                k     = bond_list(iatm1)%array2(j)
                ksign = SIGN(1.0_dp,REAL(k,KIND=dp))
                k     = ABS(k)
                CALL matvec_3x3(cell_v,topology%cell%hmat,&
                     REAL(nonbonded%neighbor_kind_pairs(k)%cell_vector,KIND=dp))
                dr = pbc_coord(:,iatm1) - pbc_coord(:,iatm2) - ksign*cell_v
                r2 = DOT_PRODUCT(dr,dr)
                IF(r2 <= 0.01_dp) THEN
                   CALL stop_program ("topology_generate_bond",&
                        "bond distance between atoms less then 0.01")
                END IF
                IF (r2 <= r2_min) THEN
                   r2_min = r2
                   ibond  = iatm2
                END IF
             END DO
             IF (ibond == -1) THEN
                IF (output_unit>0.AND.print_info) THEN
                   WRITE(output_unit,'(T2,"GENERATE|",1X,A,I10,A)')&
                        "WARNING:: No connections detected for Hydrogen - Atom Nr:",iatm1," !"
                END IF
             ELSE
                n_hydr_bonds = n_hydr_bonds + 1
                n_bonds      = n_bonds + 1
                CALL add_bonds_list(conn_info, MIN(iatm1, ibond), MAX(iatm1, ibond), n_bonds)
             END IF
          END DO          
          IF (output_unit>0) THEN
             WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')&
                  " Preliminary Number of Bonds generated:",n_bonds
          END IF
          ! External defined bonds (useful for complex connectivity)
          bond_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%GENERATE%BOND",&
               error=error)
          CALL connectivity_external_control(section=bond_section,&
                                             Iarray1=conn_info%bond_a,&
                                             Iarray2=conn_info%bond_b,&
                                             nvar=n_bonds,&
                                             topology=topology,&
                                             output_unit=output_unit,&
                                             error=error)
          ! Resize arrays to their proper size..
          CALL reallocate(conn_info%bond_a,1,n_bonds)
          CALL reallocate(conn_info%bond_b,1,n_bonds)
          IF (topology%create_molecules) THEN
             ! Since we create molecule names we're not sure that all atoms are contiguous
             ! so we need to reorder them on the basis of the generated name
             IF (.NOT.topology%reorder_atom) THEN
                topology%reorder_atom = .TRUE.
                WRITE(output_unit,'(T2,"GENERATE|",A)')&
                     " Molecules names have been generated. Now reordering particle set in order to have ",&
                     " atoms belonging to the same molecule in a sequential order."
             END IF
             connectivity_ok = .TRUE.
          ELSE
             ! Check created connectivity and possibly give the OK to proceed
             connectivity_ok = check_generate_mol(conn_info%bond_a,conn_info%bond_b,&
                  atom_info, bondparm_factor, output_unit, error)
          END IF
          IF (my_maxrad*bondparm_factor>r_max(1,1).AND.(.NOT.topology%molname_generated)) THEN
             IF (output_unit>0) THEN
                WRITE(output_unit,'(T2,"GENERATE|",A)')&
                     " ERROR in connectivity generation!",&
                     " The THRESHOLD to select possible bonds is bigger than the MAX bondlength",&
                     " used to build the neighbors lists. Increase the BONDLENGTH_MAX patameter"
                WRITE(output_unit,'(T2,"GENERATE|",2(A,F11.6),A)')&
                     " Present THRESHOLD (",my_maxrad*bondparm_factor," ). "//&
                     " Present BONDLENGTH_MAX (",r_max(1,1)," )"
             END IF
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       END DO
       IF (connectivity_ok.AND.(output_unit>0)) THEN
          WRITE(output_unit,'(T2,"GENERATE|",A)')&
               "  Achieved consistency in connectivity generation."
       END IF
       CALL fist_neighbor_deallocate(nonbonded,error)
       CALL timestop(0.0_dp,handle2)
       CALL timeset(TRIM(routineN)//"_6",'I','',handle2)
       ! Deallocate temporary working arrays
       DO I = 1, natom
          DEALLOCATE(bond_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(bond_list(I)%array2,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(bond_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(pbc_coord,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(atm_symbol,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(radius,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL deallocate_particle_set(particle_set)
       CALL deallocate_atomic_kind_set(atomic_kind_set,error)
       !
       CALL timestop(0.0_dp,handle2)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Number of Bonds generated:",&
               n_bonds
       END IF
       CALL timeset(TRIM(routineN)//"_7",'I','',handle2)
       ! If PARA_RES then activate RESIDUES
       CALL reallocate(conn_info%c_bond_a,1,0)
       CALL reallocate(conn_info%c_bond_b,1,0)
       IF(topology%para_res) THEN
          DO ibond=1,SIZE(conn_info%bond_a)
             iatom = conn_info%bond_a(ibond)
             jatom = conn_info%bond_b(ibond)
             IF((atom_info%label_molname(iatom)  /=atom_info%label_molname(jatom)).OR.&
                  (atom_info%label_resid(iatom)  /=atom_info%label_resid(jatom)).OR.&
                  (atom_info%label_resname(iatom)/=atom_info%label_resname(jatom)))THEN
                IF(iw>0) WRITE(iw,*) "      PARA_RES, bond between molecules atom ",&
                     iatom,jatom
                cbond = cbond + 1
                CALL reallocate(conn_info%c_bond_a,1,cbond)
                CALL reallocate(conn_info%c_bond_b,1,cbond)
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
       CALL timestop(0.0_dp,handle2)
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)
  END SUBROUTINE topology_generate_bond

!******************************************************************************
!!****** check_generate_mol [1.0] *
!!
!!   NAME
!!     check_generate_mol
!!
!!   SYNOPSIS
!!
!!
!!   FUNCTION
!!     Performs a check on the generated connectivity
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  FUNCTION check_generate_mol(bond_a,bond_b,atom_info,bondparm_factor,output_unit,error) &
       RESULT(conn_ok)
    INTEGER, DIMENSION(:), POINTER           :: bond_a, bond_b
    TYPE(atom_info_type), POINTER            :: atom_info
    REAL(KIND=dp), INTENT(INOUT)             :: bondparm_factor
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: conn_ok

    CHARACTER(len=*), PARAMETER :: routineN = 'check_generate_mol', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: ctmp1, ctmp2, ctmp3
    INTEGER                                  :: handle, i, idim, itype, j, &
                                                mol_natom, natom, nsize, stat
    INTEGER, DIMENSION(:), POINTER           :: mol_map, mol_map_o, wrk
    INTEGER, DIMENSION(:, :), POINTER        :: mol_info
    LOGICAL                                  :: failure, found
    LOGICAL, DIMENSION(:), POINTER           :: icheck
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: bond_list

    CALL timeset(routineN,'I','',handle)
    conn_ok = .TRUE.
    failure = .FALSE.
    IF (.NOT.failure) THEN
       natom   = SIZE(atom_info%label_atmname)
       ALLOCATE(bond_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I = 1, natom
          ALLOCATE(bond_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       CALL reorder_structure(bond_list, bond_a, bond_b, SIZE(bond_a),&
            error=error)    
       ALLOCATE(mol_map(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mol_map_o(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(wrk(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, natom
          mol_map(i) = str2id(s2s(TRIM(atom_info%label_molname(i))))
       END DO
       mol_map_o = mol_map
       CALL sort(mol_map,natom,wrk)
       !
       ! mol(i,1) : stores id of the molecule
       ! mol(i,2) : stores the total number of atoms forming that kind of molecule
       ! mol(i,3) : contains the number of molecules generated for that kind
       ! mol(i,4) : contains the number of atoms forming one molecule of that kind
       ! Connectivity will be considered correct only if for each i:
       !
       !               mol(i,2) = mol(i,3)*mol(i,4)
       !
       ! If not, very probably, a bond is missing increase bondparm by 10% and let's
       ! check if the newest connectivity is bug free..
       !
       ALLOCATE(mol_info(1,4),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       itype = mol_map(1)
       idim  = 1
       mol_info(1,1) = itype
       DO i = 2, natom
          IF (mol_map(i)/=itype) THEN
             nsize = SIZE(mol_info,1)+1
             itype = mol_map(i)
             CALL reallocate(mol_info,1,nsize,1,SIZE(mol_info,2))
             mol_info(nsize,  1) = itype
             mol_info(nsize-1,2) = idim
             idim  = 1             
          ELSE
             idim  = idim + 1
          END IF
       END DO
       nsize = SIZE(mol_info,1)
       mol_info(nsize,2) = idim
       DO i =1, nsize
          mol_info(i,3) = 0
          mol_info(i,4) = 0
       END DO
       ! 
       ALLOCATE(icheck(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       icheck = .FALSE.
       DO i = 1, natom
          IF (icheck(i)) CYCLE
          itype     = mol_map_o(i)
          mol_natom = 0
          CALL give_back_molecule(icheck,bond_list,i,mol_natom,mol_map_o,mol_map_o(i))
          DO j = 1, SIZE(mol_info)
             IF (itype==mol_info(j,1)) EXIT
          END DO
          mol_info(j,3) = mol_info(j,3) + 1
          IF (mol_info(j,4)==0) mol_info(j,4) = mol_natom
          IF (mol_info(j,4)/=mol_natom) THEN
             ! Two same molecules have been found with different number
             ! of atoms. This usually indicates a missing bond in the
             ! generated connectivity
             ! Set connectivity to .false. EXIT and increase bondparm_factor by 1.05
             conn_ok = .FALSE.
             bondparm_factor = bondparm_factor * 1.05_dp
             IF (output_unit<0) EXIT
             WRITE(output_unit,'(/,T2,"GENERATE|",A)')" WARNING in connectivity generation!"
             WRITE(output_unit,'(T2,"GENERATE|",A)')&
                  ' Two molecules/residues named ('//TRIM(id2str(itype))//') have different '//&
                  ' number of atoms.'
             CALL integer_to_string(i, ctmp1)
             CALL integer_to_string(mol_natom, ctmp2)
             CALL integer_to_string(mol_info(j,4), ctmp3)
             WRITE(output_unit,'(T2,"GENERATE|",A)')' Molecule starting at position ('//&
                  TRIM(ctmp1)//') has Nr. <'//TRIM(ctmp2)//&
                  '> of atoms.',' while the other same molecules have Nr. <'//&
                  TRIM(ctmp3)//'> of atoms!' 
             WRITE(output_unit,'(T2,"GENERATE|",A)')&
                  ' Increasing bondparm_factor by 1.05.. An error was found in the generated',&
                  ' connectivity. Retry...'
             WRITE(output_unit,'(T2,"GENERATE|",A,F11.6,A,/)')&
                  " Present value of BONDPARM_FACTOR (",bondparm_factor," )."
             EXIT
          END IF
       END DO
       DEALLOCATE(icheck,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(mol_info,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(mol_map,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(mol_map_o,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(wrk,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I = 1, natom
          DEALLOCATE(bond_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(bond_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL timestop(0.0_dp,handle)
    END IF
  END FUNCTION check_generate_mol

!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     connectivity_external_control
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Add/Remove a bond to the generated list 
!!     Particularly useful for system with complex connectivity
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE connectivity_external_control(section, Iarray1, Iarray2, Iarray3, Iarray4, nvar,&
       topology, output_unit, is_impr, error)
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, DIMENSION(:), POINTER           :: Iarray1, Iarray2
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Iarray3, Iarray4
    INTEGER, INTENT(INOUT)                   :: nvar
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN), OPTIONAL            :: is_impr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'connectivity_external_control', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=8)                         :: fmt
    INTEGER                                  :: do_action, do_it, i, j, k, &
                                                n_rep, n_rep_val, natom, &
                                                new_size, nsize, stat
    INTEGER, DIMENSION(:), POINTER           :: atlist, Ilist1, Ilist2, &
                                                Ilist3, Ilist4
    LOGICAL                                  :: explicit, failure, ip3, ip4

    failure = .FALSE.
    natom   = topology%natoms
    ! Preliminary sort of arrays
    ip3   = PRESENT(Iarray3)
    ip4   = PRESENT(Iarray4)
    nsize = 2
    IF (ip3)         nsize = nsize + 1
    IF (ip3.AND.ip4) nsize = nsize + 1
    ! Put the lists always in the canonical order
    CALL reorder_list_array(Iarray1, Iarray2, Iarray3, Iarray4, nsize, nvar, error)
    ! Go on with external control
    CALL section_vals_get(section,explicit=explicit,n_repetition=n_rep,error=error)
    IF (explicit) THEN
       NULLIFY(Ilist1, Ilist2, Ilist3, Ilist4, atlist)
       ALLOCATE(Ilist1(nvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Ilist2(nvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       Ilist1 = Iarray1(1:nvar)
       Ilist2 = Iarray2(1:nvar)
       SELECT CASE(nsize)
       CASE(2) !do nothing
       CASE(3)
          ALLOCATE(Ilist3(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Ilist3 = Iarray3(1:nvar)
       CASE(4)
          ALLOCATE(Ilist3(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Ilist4(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Ilist3 = Iarray3(1:nvar)
          Ilist4 = Iarray4(1:nvar)
       CASE DEFAULT
          ! Should never reach this point
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       CALL list_canonical_order(Ilist1, Ilist2, Ilist3, Ilist4, nsize, is_impr)
       !
       DO i = 1, n_rep
          CALL section_vals_val_get(section,"ATOMS",i_rep_section=i,n_rep_val=n_rep_val,&
               error=error)
          CALL section_vals_val_get(section,"_SECTION_PARAMETERS_",i_rep_section=i,&
               i_val=do_action, error=error)
          !
          DO j = 1, n_rep_val
             CALL section_vals_val_get(section,"ATOMS",i_rep_section=i,i_rep_val=j,&
                  i_vals=atlist,error=error)
 CPPostcondition(SIZE(atlist)==nsize,cp_failure_level,routineP,error,failure)
             CALL integer_to_string(nsize-1,fmt)
             CALL check_element_list(do_it,do_action,atlist,Ilist1,Ilist2,Ilist3,Ilist4,&
                  is_impr)
             IF (do_action==do_add) THEN
                ! Add to the element to the list
                IF (do_it>0) THEN
                   nvar = nvar + 1
                   IF (output_unit>0) THEN
                      WRITE(output_unit,'(T2,"ADD|",1X,A,I6,'//TRIM(fmt)//'(A,I6),A,T64,A,I6)')&
                           "element (",&
                           atlist(1),(",",atlist(k),k=2,nsize),") added."," NEW size::",nvar
                   END IF
                   IF (nvar > SIZE(Iarray1)) THEN
                      new_size = INT(5+1.2*nvar)
                      CALL reallocate(Iarray1,1,new_size)
                      CALL reallocate(Iarray2,1,new_size)
                      SELECT CASE(nsize)
                      CASE(3)
                         CALL reallocate(Iarray3,1,new_size)
                      CASE(4)
                         CALL reallocate(Iarray3,1,new_size)
                         CALL reallocate(Iarray4,1,new_size)
                      END SELECT
                   END IF
                   ! Using Ilist instead of atlist the canonical order is preserved..
                   Iarray1(do_it+1:nvar) = Iarray1(do_it:nvar-1)
                   Iarray2(do_it+1:nvar) = Iarray2(do_it:nvar-1)
                   Iarray1(do_it) = Ilist1(do_it)
                   Iarray2(do_it) = Ilist2(do_it)
                   SELECT CASE(nsize)
                   CASE(3)
                      Iarray3(do_it+1:nvar) = Iarray3(do_it:nvar-1)
                      Iarray3(do_it) = Ilist3(do_it)
                   CASE(4)
                      Iarray3(do_it+1:nvar) = Iarray3(do_it:nvar-1)
                      Iarray4(do_it+1:nvar) = Iarray4(do_it:nvar-1)
                      Iarray3(do_it) = Ilist3(do_it)
                      Iarray4(do_it) = Ilist4(do_it)
                   END SELECT
                ELSE
                   IF (output_unit>0) THEN
                      WRITE(output_unit,'(T2,"ADD|",1X,A,I6,'//TRIM(fmt)//'(A,I6),A,T80,A)')&
                           "element (",&
                           atlist(1),(",",atlist(k),k=2,nsize),") already found.","X"
                   END IF
                END IF
             ELSE
                ! Remove element from the list
                IF (do_it>0) THEN 
                   nvar = nvar - 1
                   IF (output_unit>0) THEN
                      WRITE(output_unit,'(T2,"RMV|",1X,A,I6,'//TRIM(fmt)//'(A,I6),A,T64,A,I6)')&
                           "element (",&
                           atlist(1),(",",atlist(k),k=2,nsize),") removed."," NEW size::",nvar
                   END IF
                   Iarray1(do_it:nvar)=Iarray1(do_it+1:nvar+1)
                   Iarray2(do_it:nvar)=Iarray2(do_it+1:nvar+1)
                   Iarray1(nvar+1)=-HUGE(0)
                   Iarray2(nvar+1)=-HUGE(0)
                   SELECT CASE(nsize)
                   CASE(3)
                      Iarray3(do_it:nvar)=Iarray3(do_it+1:nvar+1)
                      Iarray3(nvar+1)=-HUGE(0)
                   CASE(4)
                      Iarray3(do_it:nvar)=Iarray3(do_it+1:nvar+1)
                      Iarray4(do_it:nvar)=Iarray4(do_it+1:nvar+1)
                      Iarray3(nvar+1)=-HUGE(0)
                      Iarray4(nvar+1)=-HUGE(0)
                   END SELECT
                ELSE
                   IF (output_unit>0) THEN
                      WRITE(output_unit,'(T2,"RMV|",1X,A,I6,'//TRIM(fmt)//'(A,I6),A,T80,A)')&
                           "element (",&
                           atlist(1),(",",atlist(k),k=2,nsize),") not found.","X"
                   END IF
                END IF
             END IF

          END DO
       END DO
       DEALLOCATE(Ilist1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(Ilist2,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE(nsize)
       CASE(2) ! do nothing
       CASE(3)
          DEALLOCATE(Ilist3,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CASE(4)
          DEALLOCATE(Ilist3,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Ilist4,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CASE DEFAULT
          ! Should never reach this point
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF
  END SUBROUTINE connectivity_external_control

!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     list_canonical_order
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Orders list in the canonical order: the extrema of the list are such
!!     that the first extrema is always smaller or equal to the last extrema.
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE list_canonical_order(Ilist1, Ilist2, Ilist3, Ilist4, nsize, is_impr)
    INTEGER, DIMENSION(:), POINTER           :: Ilist1, Ilist2
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Ilist3, Ilist4
    INTEGER, INTENT(IN)                      :: nsize
    LOGICAL, INTENT(IN), OPTIONAL            :: is_impr

    INTEGER                                  :: i, ss(3), tmp1, tmp2, tmp3, &
                                                tt(3)
    LOGICAL                                  :: do_impr

    do_impr= .FALSE.
    IF (PRESENT(is_impr)) do_impr = is_impr    
    SELECT CASE(nsize)
    CASE(2)
       DO i = 1, SIZE(Ilist1)
          tmp1 = Ilist1(i)
          tmp2 = Ilist2(i)
          Ilist1(i) = MIN(tmp1,tmp2)
          Ilist2(i) = MAX(tmp1,tmp2)
       END DO
    CASE(3)
       DO i = 1, SIZE(Ilist1)
          tmp1 = Ilist1(i)
          tmp2 = Ilist3(i)
          Ilist1(i) = MIN(tmp1,tmp2)
          Ilist3(i) = MAX(tmp1,tmp2)
       END DO
    CASE(4)
       DO i = 1, SIZE(Ilist1)
          IF (.NOT.do_impr) THEN
             tmp1 = Ilist1(i)
             tmp2 = Ilist4(i)
             Ilist1(i) = MIN(tmp1,tmp2)
             IF (Ilist1(i)==tmp2) THEN
                tmp3 = Ilist3(i)
                Ilist3(i) = Ilist2(i)
                Ilist2(i) = tmp3
             END IF
             Ilist4(i) = MAX(tmp1,tmp2)
          ELSE
             tt(1) = Ilist2(i)
             tt(2) = Ilist3(i)
             tt(3) = Ilist4(i)
             CALL sort(tt,3,ss)
             Ilist2(i) = tt(1)
             Ilist3(i) = tt(2)
             Ilist4(i) = tt(3)
          END IF
       END DO
    END SELECT

  END SUBROUTINE list_canonical_order

!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     check_element_list
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     finds an element in the ordered list
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE check_element_list(do_it,do_action,atlist,Ilist1,Ilist2,Ilist3,Ilist4,&
       is_impr)
    INTEGER, INTENT(OUT)                     :: do_it
    INTEGER, INTENT(IN)                      :: do_action
    INTEGER, DIMENSION(:), POINTER           :: atlist, Ilist1, Ilist2
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Ilist3, Ilist4
    LOGICAL, INTENT(IN), OPTIONAL            :: is_impr

    INTEGER                                  :: i, iend, istart, ndim, &
                                                new_size, nsize, ss(3), tmp1, &
                                                tmp2, tmp3, tt(3)
    INTEGER, DIMENSION(4)                    :: tmp
    LOGICAL                                  :: do_impr, found

    do_impr= .FALSE.
    IF (PRESENT(is_impr)) do_impr = is_impr
    found = .FALSE.
    nsize = SIZE(atlist)
    ndim  = SIZE(Ilist1)
    DO i = 1, nsize
       tmp(i)=atlist(i)
    END DO
    SELECT CASE(nsize)
    CASE(2)
       tmp1 = tmp(1)
       tmp2 = tmp(2)
       tmp(1) = MIN(tmp1,tmp2)
       tmp(2) = MAX(tmp1,tmp2)
    CASE(3)
       tmp1 = tmp(1)
       tmp2 = tmp(3)
       tmp(1) = MIN(tmp1,tmp2)
       tmp(3) = MAX(tmp1,tmp2)
    CASE(4)
       IF (.NOT.do_impr) THEN
          tmp1 = tmp(1)
          tmp2 = tmp(4)
          tmp(1) = MIN(tmp1,tmp2)
          IF (tmp(1)==tmp2) THEN
             tmp3   = tmp(3)
             tmp(3) = tmp(2)
             tmp(2) = tmp3
          END IF
          tmp(4) = MAX(tmp1,tmp2)
       ELSE
          tt(1) = tmp(2)
          tt(2) = tmp(3)
          tt(3) = tmp(4)
          CALL sort(tt,3,ss)
          tmp(2) = tt(1)
          tmp(3) = tt(2)
          tmp(4) = tt(3)
       END IF
    END SELECT
    ! boundary to search
    DO istart = 1, ndim
       IF (Ilist1(istart)>=tmp(1)) EXIT
    END DO
    IF (Ilist1(istart)>tmp(1).AND.(istart/=1)) istart=istart-1
    DO iend = istart, ndim
       IF (Ilist1(iend)/=tmp(1)) EXIT
    END DO
    IF (iend==ndim+1)  iend = ndim
    ! Final search in array
    SELECT CASE(nsize)
    CASE(2)
       DO i = istart, iend
          IF ((Ilist1(i)>tmp(1)).OR.(Ilist2(i)>tmp(2))) EXIT
          IF ((Ilist1(i)==tmp(1)).AND.(Ilist2(i)==tmp(2))) THEN
             found = .TRUE.
             EXIT
          END IF
       END DO
    CASE(3)
       DO i = istart, iend
          IF ((Ilist1(i)>tmp(1)).OR.(Ilist2(i)>tmp(2)).OR.(Ilist3(i)>tmp(3))) EXIT
          IF ((Ilist1(i)==tmp(1)).AND.(Ilist2(i)==tmp(2)).AND.(Ilist3(i)==tmp(3))) THEN
             found = .TRUE.
             EXIT
          END IF
       END DO       
    CASE(4)
       DO i = istart, iend
          IF ((Ilist1(i)>tmp(1)).OR.(Ilist2(i)>tmp(2)).OR.(Ilist3(i)>tmp(3)).OR.(Ilist4(i)>tmp(4))) EXIT
          IF ((Ilist1(i)==tmp(1)).AND.(Ilist2(i)==tmp(2))&
               .AND.(Ilist3(i)==tmp(3)).AND.(Ilist4(i)==tmp(4))) THEN
             found = .TRUE.
             EXIT
          END IF
       END DO
    END SELECT
    SELECT CASE (do_action)
    CASE (do_add)
       IF (found) THEN
          do_it = -i
          ! Nothing to modify. Element already present
          ! in this case ABS(do_it) gives the exact location of the element
          ! in the list
       ELSE
          ! Let's add the element in the right place of the list.. so that we can keep the
          ! canonical order
          ! in this case do_it gives the index of the list with indexes bigger than
          ! the one we're searching for
          ! At the end do_it gives the exact location of the element in the canonical list
          do_it =  i
          new_size = ndim+1
          CALL reallocate(Ilist1,1,new_size)
          CALL reallocate(Ilist2,1,new_size)
          Ilist1(i+1:new_size)=Ilist1(i:ndim)
          Ilist2(i+1:new_size)=Ilist2(i:ndim)
          Ilist1(i) = tmp(1)
          Ilist2(i) = tmp(2)
          SELECT CASE(nsize)
          CASE(3)
             CALL reallocate(Ilist3,1,new_size)
             Ilist3(i+1:new_size)=Ilist3(i:ndim)
             Ilist3(i) = tmp(3)
          CASE(4)
             CALL reallocate(Ilist3,1,new_size)
             CALL reallocate(Ilist4,1,new_size)
             Ilist3(i+1:new_size)=Ilist3(i:ndim)
             Ilist4(i+1:new_size)=Ilist4(i:ndim)
             Ilist3(i) = tmp(3)
             Ilist4(i) = tmp(4)
          END SELECT
       END IF
    CASE(do_remove)
       IF (found) THEN
          do_it =   i
          ! Let's delete the element in position do_it
          new_size = ndim-1
          Ilist1(i:new_size)=Ilist1(i+1:ndim)
          Ilist2(i:new_size)=Ilist2(i+1:ndim)
          CALL reallocate(Ilist1,1,new_size)
          CALL reallocate(Ilist2,1,new_size)
          SELECT CASE(nsize)
          CASE(3)
             Ilist3(i:new_size)=Ilist3(i+1:ndim)
             CALL reallocate(Ilist3,1,new_size)
          CASE(4)
             Ilist3(i:new_size)=Ilist3(i+1:ndim)
             Ilist4(i:new_size)=Ilist4(i+1:ndim)
             CALL reallocate(Ilist3,1,new_size)
             CALL reallocate(Ilist4,1,new_size)
          END SELECT
       ELSE
          do_it =  -i
          ! Nothing to modify. Element not present in the list
          ! in this case ABS(do_it) gives the exact location of the element
          ! in the list
       END IF
    END SELECT
  END SUBROUTINE check_element_list

!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     add_bonds_list
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Adds a bond to the generated bond list 
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE add_bonds_list(conn_info, atm1, atm2, n_bonds)
    TYPE(connectivity_info_type), POINTER    :: conn_info
    INTEGER, INTENT(IN)                      :: atm1, atm2, n_bonds

    INTEGER                                  :: new_size, old_size

    old_size      = SIZE(conn_info%bond_a)
    IF(n_bonds > old_size) THEN
       new_size=INT(5+1.2*old_size)
       CALL reallocate(conn_info%bond_a,1,new_size)
       CALL reallocate(conn_info%bond_b,1,new_size)
    END IF
    conn_info%bond_a(n_bonds) = atm1
    conn_info%bond_b(n_bonds) = atm2
  END SUBROUTINE add_bonds_list

!******************************************************************************
!!****** topology_generate_util [1.0] *
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
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_bend(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bend', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, i, iw, &
                                                natom, nbond, nsize, ntheta, &
                                                output_unit, stat
    LOGICAL                                  :: failure, ionode
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: bond_list
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: bend_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,'I','',handle)
    IF (.NOT.failure) THEN
       ionode = logger%para_env%mepos==logger%para_env%source
       IF (ionode) THEN
          output_unit= cp_logger_get_default_unit_nr(logger)
       ELSE
          output_unit = -1
       END IF
       conn_info => topology%conn_info
       nbond = 0
       natom = topology%natoms
       ! This call is for connectivity off
       IF(ASSOCIATED(conn_info%bond_a)) THEN
          nbond = SIZE(conn_info%bond_a)
       ELSE
          CALL reallocate(conn_info%bond_a,1,nbond)
          CALL reallocate(conn_info%bond_b,1,nbond)
       END IF
       ntheta = 0
       nsize  = INT(5+1.2*ntheta)
       CALL reallocate(conn_info%theta_a,1,nsize)
       CALL reallocate(conn_info%theta_b,1,nsize)
       CALL reallocate(conn_info%theta_c,1,nsize)
       ! Get list of bonds to pre-process theta
       ALLOCATE(bond_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I = 1, natom
          ALLOCATE(bond_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       CALL reorder_structure(bond_list, conn_info%bond_a, conn_info%bond_b, nbond, error)
       ! All the dirty job is handled by this routine.. for bends it_levl is equal 3
       CALL timeset(routineN//"_1",'I','',handle2)
       CALL match_iterative_path(Iarray1=bond_list,&
                                 Iarray2=bond_list,&
                                 max_levl=3,&
                                 nvar=ntheta,&
                                 Oarray1 =conn_info%theta_a,&
                                 Oarray2 =conn_info%theta_b,&
                                 Oarray3 =conn_info%theta_c,&
                                 error=error)
       CALL timestop(0.0_dp,handle2)
       DO I = 1, natom
          DEALLOCATE(bond_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(bond_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Preliminary Number of Bends generated:",&
               ntheta
       END IF
       ! External defined bends (useful for complex connectivity)
       bend_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%GENERATE%ANGLE",error=error)
       CALL connectivity_external_control(section=bend_section,&
                                          Iarray1=conn_info%theta_a,&
                                          Iarray2=conn_info%theta_b,&
                                          Iarray3=conn_info%theta_c,&
                                          nvar=ntheta,&
                                          topology=topology,&
                                          output_unit=output_unit,&
                                          error=error)
       ! Resize arrays to their proper size..
       CALL reallocate(conn_info%theta_a,1,ntheta)
       CALL reallocate(conn_info%theta_b,1,ntheta)
       CALL reallocate(conn_info%theta_c,1,ntheta)
       
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Number of Bends generated:",&
               ntheta
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)
  END SUBROUTINE topology_generate_bend

!
!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     match_iterative_path
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Routine matching iteratively along a graph 
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  RECURSIVE SUBROUTINE match_iterative_path(Iarray1, Iarray2, Iarray3,&
       max_levl, Oarray1, Oarray2, Oarray3, Oarray4, Ilist, it_levl, nvar,&
       error)
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: Iarray1
    TYPE(array1_list_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: Iarray2, Iarray3
    INTEGER, INTENT(IN)                      :: max_levl
    INTEGER, DIMENSION(:), POINTER           :: Oarray1, Oarray2
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Oarray3, Oarray4
    INTEGER, DIMENSION(:), INTENT(INOUT), &
      OPTIONAL                               :: Ilist
    INTEGER, INTENT(IN), OPTIONAL            :: it_levl
    INTEGER, INTENT(INOUT)                   :: nvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'match_iterative_path', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind, j, my_levl, natom, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: my_list
    LOGICAL                                  :: check, failure
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: wrk

    failure = .FALSE.
    check   = max_levl>=2.AND.max_levl<=4
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    IF (.NOT.PRESENT(Ilist)) THEN
       SELECT CASE(max_levl)
       CASE(2)
    CPPostcondition(.NOT.PRESENT(Iarray2),cp_failure_level,routineP,error,failure) 
    CPPostcondition(.NOT.PRESENT(Iarray3),cp_failure_level,routineP,error,failure) 
    CPPostcondition(.NOT.PRESENT(Oarray3),cp_failure_level,routineP,error,failure) 
    CPPostcondition(.NOT.PRESENT(Oarray4),cp_failure_level,routineP,error,failure) 
       CASE(3)
    CPPostcondition(PRESENT(Iarray2),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.PRESENT(Iarray3),cp_failure_level,routineP,error,failure) 
    CPPostcondition(PRESENT(Oarray3),cp_failure_level,routineP,error,failure) 
    CPPostcondition(.NOT.PRESENT(Oarray4),cp_failure_level,routineP,error,failure) 
       CASE(4)
    CPPostcondition(PRESENT(Iarray2),cp_failure_level,routineP,error,failure)
    CPPostcondition(PRESENT(Iarray3),cp_failure_level,routineP,error,failure)
    CPPostcondition(PRESENT(Oarray3),cp_failure_level,routineP,error,failure) 
    CPPostcondition(PRESENT(Oarray4),cp_failure_level,routineP,error,failure) 
       END SELECT
    END IF
    natom   = SIZE(Iarray1)
    IF (.NOT.failure) THEN
       IF (.NOT.PRESENT(Ilist)) THEN
          ! Start a new loop.. Only the first time the routine is called
          ALLOCATE(my_list(max_levl),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, natom
             my_levl =  1
             my_list = -1
             my_list(my_levl) = i
             CALL match_iterative_path(Iarray1=Iarray1,&
                                       Iarray2=Iarray2,&
                                       Iarray3=Iarray3,&
                                       it_levl=my_levl+1,&
                                       max_levl=max_levl,&
                                       Oarray1=Oarray1,&
                                       Oarray2=Oarray2,&
                                       Oarray3=Oarray3,&
                                       Oarray4=Oarray4,&
                                       nvar   =nvar,&
                                       Ilist  =my_list,&
                                       error  =error)
          END DO
          DEALLOCATE(my_list,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          SELECT CASE(it_levl)
          CASE(2)
             wrk => Iarray1
          CASE(3)
             wrk => Iarray2
          CASE(4)
             wrk => Iarray3
          END SELECT
          i = Ilist(it_levl-1)
          DO j = 1, SIZE(Iarray1(i)%array1)
             ind = wrk(i)%array1(j)
             IF (ANY(Ilist==ind)) CYCLE
             IF (it_levl<max_levl) THEN
                Ilist(it_levl) = ind
                CALL match_iterative_path(Iarray1=Iarray1,&
                                          Iarray2=Iarray2,&
                                          Iarray3=Iarray3,&
                                          it_levl=it_levl+1,&
                                          max_levl=max_levl,&
                                          Oarray1=Oarray1,&
                                          Oarray2=Oarray2,&
                                          Oarray3=Oarray3,&
                                          Oarray4=Oarray4,&
                                          nvar   =nvar,&
                                          Ilist  =Ilist,&
                                          error  =error)
                Ilist(it_levl) = -1
             ELSEIF (it_levl==max_levl) THEN
                IF (Ilist(1)>ind) CYCLE
                Ilist(it_levl) = ind
                nvar = nvar + 1
                SELECT CASE(it_levl)
                CASE(2)
                   IF(nvar > SIZE(Oarray1)) THEN
                      CALL reallocate(Oarray1,1,INT(5+1.2*nvar))
                      CALL reallocate(Oarray2,1,INT(5+1.2*nvar))
                   END IF
                   Oarray1(nvar) = Ilist(1)
                   Oarray2(nvar) = Ilist(2)
                CASE(3)
                   IF(nvar > SIZE(Oarray1)) THEN
                      CALL reallocate(Oarray1,1,INT(5+1.2*nvar))
                      CALL reallocate(Oarray2,1,INT(5+1.2*nvar))
                      CALL reallocate(Oarray3,1,INT(5+1.2*nvar))
                   END IF
                   Oarray1(nvar) = Ilist(1)
                   Oarray2(nvar) = Ilist(2)
                   Oarray3(nvar) = Ilist(3)
                CASE(4)
                   IF(nvar > SIZE(Oarray1)) THEN
                      CALL reallocate(Oarray1,1,INT(5+1.2*nvar))
                      CALL reallocate(Oarray2,1,INT(5+1.2*nvar))
                      CALL reallocate(Oarray3,1,INT(5+1.2*nvar))
                      CALL reallocate(Oarray4,1,INT(5+1.2*nvar))
                   END IF
                   Oarray1(nvar) = Ilist(1)
                   Oarray2(nvar) = Ilist(2)
                   Oarray3(nvar) = Ilist(3)
                   Oarray4(nvar) = Ilist(4)
                CASE DEFAULT
                   !should never reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END SELECT
                Ilist(it_levl) = -1
             ELSE
                !should never reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE match_iterative_path
    
!
!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     topology_generate_ub
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     The list of Urey-Bradley is equal to the list of bends
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_ub(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_ub', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, itheta, iw, ntheta, &
                                                output_unit
    LOGICAL                                  :: failure, ionode
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error)    
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
       output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
       output_unit = -1
    END IF
    CALL timeset(routineN,'I','',handle)
    conn_info => topology%conn_info    
    ntheta = SIZE(conn_info%theta_a)
    CALL reallocate(conn_info%ub_a,1,ntheta)
    CALL reallocate(conn_info%ub_b,1,ntheta)
    CALL reallocate(conn_info%ub_c,1,ntheta)

    DO itheta=1,ntheta
       conn_info%ub_a(itheta) = conn_info%theta_a(itheta)
       conn_info%ub_b(itheta) = conn_info%theta_b(itheta)
       conn_info%ub_c(itheta) = conn_info%theta_c(itheta)
    END DO
    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Number of UB generated:",&
            ntheta
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)

  END SUBROUTINE topology_generate_ub

!******************************************************************************
!!****** topology_generate_util [1.0] *
!!
!!   NAME
!!     topology_generate_dihe
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Generate a list of torsions from bonds    
!! 
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_dihe(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_dihe', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    INTEGER                                  :: handle, i, iw, natom, nbond, &
                                                nphi, nsize, output_unit, stat
    LOGICAL                                  :: failure, ionode
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: bond_list
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: torsion_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
       output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
       output_unit = -1
    END IF
    CALL timeset(routineN,'I','',handle)
    IF (.NOT.failure) THEN
       conn_info => topology%conn_info
       nphi   = 0
       nsize  = INT(5+1.2*nphi)
       CALL reallocate(conn_info%phi_a,1,nsize)
       CALL reallocate(conn_info%phi_b,1,nsize)
       CALL reallocate(conn_info%phi_c,1,nsize)
       CALL reallocate(conn_info%phi_d,1,nsize)
       ! Get list of bonds to pre-process phi
       natom = topology%natoms
       ALLOCATE(bond_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I = 1, natom
          ALLOCATE(bond_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       nbond = SIZE(conn_info%bond_a)
       CALL reorder_structure(bond_list, conn_info%bond_a, conn_info%bond_b, nbond, error)
       ! All the dirty job is handled by this routine.. for torsions it_levl is equal 4
       CALL match_iterative_path(Iarray1=bond_list,&
                                 Iarray2=bond_list,&
                                 Iarray3=bond_list,&
                                 max_levl=4,&
                                 nvar=nphi,&
                                 Oarray1 =conn_info%phi_a,&
                                 Oarray2 =conn_info%phi_b,&
                                 Oarray3 =conn_info%phi_c,&
                                 Oarray4 =conn_info%phi_d,&
                                 error=error)
       DO I = 1, natom
          DEALLOCATE(bond_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(bond_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Preliminary Number of Torsions generated:",&
               nphi
       END IF
       ! External defined torsions (useful for complex connectivity)
       torsion_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%GENERATE%TORSION",error=error)
       CALL connectivity_external_control(section=torsion_section,&
                                          Iarray1=conn_info%phi_a,&
                                          Iarray2=conn_info%phi_b,&
                                          Iarray3=conn_info%phi_c,&
                                          Iarray4=conn_info%phi_d,&
                                          nvar=nphi,&
                                          topology=topology,&
                                          output_unit=output_unit,&
                                          error=error)
       ! Resize arrays to their proper size..
       CALL reallocate(conn_info%phi_a,1,nphi)
       CALL reallocate(conn_info%phi_b,1,nphi)
       CALL reallocate(conn_info%phi_c,1,nphi)
       CALL reallocate(conn_info%phi_d,1,nphi)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Number of Torsions generated:",&
               nphi
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)
    
  END SUBROUTINE topology_generate_dihe

!******************************************************************************
!!****** topology_generate_util [1.0] *
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
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_impr(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_impr', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: atm_symbol
    INTEGER                                  :: handle, i, ind, iw, j, natom, &
                                                nbond, nimpr, nsize, &
                                                output_unit, stat
    LOGICAL                                  :: accept_impr, failure, ionode
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: bond_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: impr_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
       output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
       output_unit = -1
    END IF
    CALL timeset(routineN,'I','',handle)
    IF (.NOT.failure) THEN
       atom_info => topology%atom_info
       conn_info => topology%conn_info       
       natom  = topology%natoms
       nbond = SIZE(conn_info%bond_a)
       nimpr  = 0 
       nsize  = INT(5+1.2*nimpr)
       CALL reallocate(conn_info%impr_a,1,nsize)
       CALL reallocate(conn_info%impr_b,1,nsize)
       CALL reallocate(conn_info%impr_c,1,nsize)
       CALL reallocate(conn_info%impr_d,1,nsize)
       ! Get list of bonds to pre-process phi
       ALLOCATE(bond_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I = 1, natom
          ALLOCATE(bond_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       nbond = SIZE(conn_info%bond_a)
       CALL reorder_structure(bond_list, conn_info%bond_a, conn_info%bond_b, nbond, error)
       DO I = 1, natom
          ! Count all atoms with three bonds
          IF (SIZE(bond_list(I)%array1)==3) THEN
             ! Problematic cases::
             ! Nitrogen
             accept_impr = .TRUE.
             atm_symbol  = atom_info%element(i)
             CALL uppercase(atm_symbol)
             IF (atm_symbol=="N ") THEN
                accept_impr = .FALSE.
                ! Impropers on Nitrogen only when there is another atom close to it
                ! with other 3 bonds 
                DO j = 1, 3
                   ind = bond_list(I)%array1(j)
                   IF (SIZE(bond_list(ind)%array1)==3) accept_impr = .TRUE.
                END DO
             END IF
             IF (.NOT.accept_impr) CYCLE
             nimpr = nimpr + 1
             IF (nimpr>SIZE(conn_info%impr_a)) THEN
                nsize  = INT(5+1.2*nimpr)
                CALL reallocate(conn_info%impr_a,1,nsize)
                CALL reallocate(conn_info%impr_b,1,nsize)
                CALL reallocate(conn_info%impr_c,1,nsize)
                CALL reallocate(conn_info%impr_d,1,nsize)
             END IF
             conn_info%impr_a(nimpr) = i
             conn_info%impr_b(nimpr) = bond_list(I)%array1(1)
             conn_info%impr_c(nimpr) = bond_list(I)%array1(2)
             conn_info%impr_d(nimpr) = bond_list(I)%array1(3)
          END IF
       END DO
       DO I = 1, natom
          DEALLOCATE(bond_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(bond_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       ! External defined impropers (useful for complex connectivity)
       impr_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%GENERATE%IMPROPER",error=error)
       CALL connectivity_external_control(section=impr_section,&
                                          Iarray1=conn_info%impr_a,&
                                          Iarray2=conn_info%impr_b,&
                                          Iarray3=conn_info%impr_c,&
                                          Iarray4=conn_info%impr_d,&
                                          nvar=nimpr,&
                                          topology=topology,&
                                          output_unit=output_unit,&
                                          is_impr=.TRUE.,&
                                          error=error)
       ! Resize arrays to their proper size..
       CALL reallocate(conn_info%impr_a,1,nimpr)
       CALL reallocate(conn_info%impr_b,1,nimpr)
       CALL reallocate(conn_info%impr_c,1,nimpr)
       CALL reallocate(conn_info%impr_d,1,nimpr)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Number of Impropers generated:",&
               nimpr
       END IF       
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)

  END SUBROUTINE topology_generate_impr

!******************************************************************************
!!****** topology_generate_util [1.0] *
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
  SUBROUTINE topology_generate_onfo(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_onfo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iphi, iw, nphi, &
                                                output_unit
    LOGICAL                                  :: failure, ionode
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error) 
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
       output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
       output_unit = -1
    END IF
    CALL timeset(routineN,'I','',handle)
    
    conn_info => topology%conn_info    
    nphi = SIZE(conn_info%phi_a)
    CALL reallocate(conn_info%onfo_a,1,nphi)
    CALL reallocate(conn_info%onfo_b,1,nphi)

    DO iphi=1,nphi
       conn_info%onfo_a(iphi) = conn_info%phi_a(iphi)
       conn_info%onfo_b(iphi) = conn_info%phi_d(iphi)
    END DO
    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,"GENERATE|",1X,A,T71,I10)')" Number of 1-4 interactions generated:",&
            nphi
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)
    
  END SUBROUTINE topology_generate_onfo    


END MODULE topology_generate_util
