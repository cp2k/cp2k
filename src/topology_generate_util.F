!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
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
!!    jgh (23-05-2004) Last atom of molecule information added
!!
!!   SOURCE
!******************************************************************************

MODULE topology_generate_util
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_bondparm_covalent,&
                                             do_bondparm_vdw,&
                                             do_conn_off
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: generate_particle_distances
  USE periodic_table,                  ONLY: ptable
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: atom_bond_list_type,&
                                             find_molecule,&
                                             reorder_structure
  USE util,                            ONLY: find_boundary
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_in_kind, atom_in_mol, first, handle, handle2, i, iatm, &
      iatom, ifirst, ilast, inum, iw, j, jump1, jump2, last, max_mol_num, &
      mol_num, mol_res, mol_typ, myind, N, natom, ntype, stat
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: do_again, failure, my_qmmm
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
       ALLOCATE(atom_bond_list(I)%bonds(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%bond_a)) N = SIZE(conn_info%bond_a)
    CALL reorder_structure(atom_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)
    CALL find_molecule(atom_bond_list,atom_info%map_mol_num)   
    DO I=1,natom
       DEALLOCATE(atom_bond_list(I)%bonds,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(atom_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    ! Modify according map_mol_typ the array map_mol_num
    !
    IF(iw>0) WRITE(iw,'(/,A)') "End of molecule generation"
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
          WRITE(iw,'(4(A," : ",I7))') "iatom",iatom,&
                                      "map_mol_typ",atom_info%map_mol_typ(iatom),&
                                      "map_mol_num",atom_info%map_mol_num(iatom),&
                                      "map_mol_res",atom_info%map_mol_res(iatom)
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
          !
          DO j = 1, max_mol_num
             IF (COUNT(atom_info%map_mol_num(first:last)==j) /= atom_in_mol) THEN
                WRITE(iw,*)"molecule type:",i,"molecule num:",j," has ",COUNT(atom_info%map_mol_num(first:last)==j),&
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
!!     
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE topology_generate_bond(topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
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
       CALL generate_particle_distances(atom_info%r,&
            atom_info%label_molname,&
            topology%cell,radius,&
            topology%bondparm_factor,&
            conn_info%bond_a,&
            conn_info%bond_b,&
            topology%bondparm_type,&
            para_env)
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
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    DEALLOCATE(atm_symbol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)

  END SUBROUTINE topology_generate_bond

!
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
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
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)
  END SUBROUTINE topology_generate_bend

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
!!     Using a list of bend, generate a list of Urey-Bradley
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
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
!!     Using a list of bend, generate a list of dihe
!!
!!   AUTHOR
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
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
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)

  END SUBROUTINE topology_generate_impr

!
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GENERATE_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    conn_info => topology%conn_info    
    nphi = SIZE(conn_info%phi_a)
    CALL reallocate(conn_info%onfo_a,1,nphi)
    CALL reallocate(conn_info%onfo_b,1,nphi)

    DO iphi=1,nphi
       conn_info%onfo_a(iphi) = conn_info%phi_a(iphi)
       conn_info%onfo_b(iphi) = conn_info%phi_d(iphi)
    END DO
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GENERATE_INFO",error=error)
    
  END SUBROUTINE topology_generate_onfo    


END MODULE topology_generate_util
