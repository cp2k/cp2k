!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the molecule kind structure types and the corresponding
!>      functionality
!> \par History
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
!> \author MK (22.08.2003)
! *****************************************************************************
MODULE molecule_kind_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE colvar_types,                    ONLY: &
       angle_colvar_id, colvar_counters, combine_colvar_id, coord_colvar_id, &
       dfunct_colvar_id, dist_colvar_id, gyration_colvar_id, &
       hydronium_colvar_id, plane_distance_colvar_id, &
       plane_plane_angle_colvar_id, population_colvar_id, qparm_colvar_id, &
       reaction_path_colvar_id, rotation_colvar_id, torsion_colvar_id, &
       xyz_diag_colvar_id, xyz_outerdiag_colvar_id
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE force_field_kind_types,          ONLY: &
       bend_kind_type, bond_kind_type, impr_kind_dealloc_ref, impr_kind_type, &
       opbend_kind_type, torsion_kind_dealloc_ref, torsion_kind_type, &
       ub_kind_dealloc_ref, ub_kind_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE shell_potential_types,           ONLY: shell_kind_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molecule_kind_types'

! *** Define the derived structure types ***

! *****************************************************************************
  TYPE atom_type
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER :: id_name
  END TYPE atom_type

! *****************************************************************************
  TYPE shell_type
    INTEGER :: a
    CHARACTER(LEN=default_string_length)  :: name
    TYPE(shell_kind_type), POINTER :: shell_kind
  END TYPE shell_type

! *****************************************************************************
  TYPE bond_type
    INTEGER :: a,b
    INTEGER :: id_type, itype
    TYPE(bond_kind_type), POINTER :: bond_kind
  END TYPE bond_type

! *****************************************************************************
  TYPE bend_type
    INTEGER :: a,b,c
    INTEGER :: id_type, itype
    TYPE(bend_kind_type), POINTER :: bend_kind
  END TYPE bend_type

! *****************************************************************************
  TYPE ub_type
    INTEGER :: a,b,c
    INTEGER :: id_type, itype
    TYPE(ub_kind_type), POINTER :: ub_kind
  END TYPE ub_type

! *****************************************************************************
  TYPE torsion_type
    INTEGER :: a,b,c,d
    INTEGER :: id_type, itype
    TYPE(torsion_kind_type), POINTER :: torsion_kind
  END TYPE torsion_type

! *****************************************************************************
  TYPE impr_type
    INTEGER :: a,b,c,d
    INTEGER :: id_type, itype
    TYPE(impr_kind_type), POINTER :: impr_kind
  END TYPE impr_type

! *****************************************************************************
  TYPE opbend_type
    INTEGER :: a,b,c,d
    INTEGER :: id_type, itype
    TYPE(opbend_kind_type), POINTER :: opbend_kind
  END TYPE opbend_type

! *****************************************************************************
  TYPE restraint_type
     LOGICAL       :: active
     REAL(KIND=dp) :: k0
  END TYPE restraint_type

! *****************************************************************************
  TYPE colvar_constraint_type
     INTEGER                        :: type_id
     INTEGER                        :: inp_seq_num
     LOGICAL                        :: use_points
     REAL(KIND = dp)                :: expected_value
     REAL(KIND = dp)                :: expected_value_growth_speed
     INTEGER, POINTER, DIMENSION(:) :: i_atoms
     TYPE(restraint_type)           :: restraint
  END TYPE colvar_constraint_type

! *****************************************************************************
  TYPE g3x3_constraint_type
     INTEGER                        :: a,b,c
     REAL(KIND = dp)                :: dab,dac,dbc
     TYPE(restraint_type)           :: restraint
  END TYPE g3x3_constraint_type

! *****************************************************************************
  TYPE g4x6_constraint_type
     INTEGER                        :: a,b,c,d
     REAL(KIND = dp)                :: dab,dac,dbc,dad,dbd,dcd
     TYPE(restraint_type)           :: restraint
  END TYPE g4x6_constraint_type

! *****************************************************************************
  TYPE vsite_constraint_type
     INTEGER                        :: a,b,c,d
     REAL(KIND = dp)                :: wbc, wdc
     TYPE(restraint_type)           :: restraint
  END TYPE vsite_constraint_type

! *****************************************************************************
  TYPE fixd_constraint_type
     TYPE(restraint_type)           :: restraint
     INTEGER                        :: fixd, itype
     REAL(KIND=dp), DIMENSION(3)    :: coord
  END TYPE fixd_constraint_type

! *****************************************************************************
  TYPE local_fixd_constraint_type
     INTEGER                        :: ifixd_index, ikind
  END TYPE local_fixd_constraint_type

! *****************************************************************************
  TYPE molecule_kind_type
    TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
    TYPE(bond_kind_type), DIMENSION(:), POINTER       :: bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
    TYPE(bend_kind_type), DIMENSION(:), POINTER       :: bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
    TYPE(ub_kind_type), DIMENSION(:), POINTER         :: ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER              :: ub_list
    TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
    TYPE(impr_kind_type), DIMENSION(:), POINTER       :: impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER            :: impr_list
    TYPE(opbend_kind_type), DIMENSION(:), POINTER     :: opbend_kind_set
    TYPE(opbend_type), DIMENSION(:), POINTER          :: opbend_list
    TYPE(colvar_constraint_type), DIMENSION(:), POINTER :: colv_list
    TYPE(g3x3_constraint_type), DIMENSION(:), POINTER   :: g3x3_list
    TYPE(g4x6_constraint_type), DIMENSION(:), POINTER   :: g4x6_list
    TYPE(vsite_constraint_type), DIMENSION(:), POINTER  :: vsite_list
    TYPE(fixd_constraint_type), DIMENSION(:), POINTER   :: fixd_list
    TYPE(shell_type), DIMENSION(:), POINTER           :: shell_list
    CHARACTER(LEN=default_string_length)              :: name
    REAL(KIND = dp)                                   :: charge,&
                                                         mass
    INTEGER                                           :: kind_number,&
                                                         natom,&
                                                         nbond,&
                                                         nbend,&
                                                         nimpr,&
                                                         nopbend,&
                                                         ntorsion,&
                                                         nub,&
                                                         ng3x3, ng3x3_restraint,&
                                                         ng4x6, ng4x6_restraint,&
                                                         nvsite, nvsite_restraint,&
                                                         nfixd, nfixd_restraint,&
                                                         nmolecule, nshell
    TYPE(colvar_counters)                             :: ncolv
    INTEGER                                           :: nsgf, nelectron,&
                                                         nelectron_alpha,&
                                                         nelectron_beta
    INTEGER, DIMENSION(:), POINTER                    :: molecule_list
    LOGICAL                                           :: molname_generated
  END TYPE molecule_kind_type

  ! *** Public subroutines ***
  PUBLIC :: allocate_molecule_kind_set,&
            deallocate_molecule_kind_set,&
            get_molecule_kind,&
            get_molecule_kind_set,&
            set_molecule_kind,&
            write_molecule_kind_set,&
            setup_colvar_counters

  ! *** Public data types ***
  PUBLIC :: atom_type,&
            bend_type,&
            bond_type,&
            ub_type,&
            torsion_type,&
            impr_type,&
            opbend_type,&
            colvar_constraint_type,&
            g3x3_constraint_type,&
            g4x6_constraint_type,&
            vsite_constraint_type,&
            fixd_constraint_type,&
            local_fixd_constraint_type,&
            molecule_kind_type,&
            shell_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param colv_list ...
!> \param ncolv ...
! *****************************************************************************
  SUBROUTINE setup_colvar_counters(colv_list, ncolv)
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(colvar_counters)                    :: ncolv

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_colvar_counters', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k
    LOGICAL                                  :: failure

    failure=.FALSE.
    ncolv%ndist          = 0
    ncolv%nangle         = 0
    ncolv%ndfunct        = 0
    ncolv%ntorsion       = 0
    ncolv%ncoord         = 0
    ncolv%nplane_dist    = 0
    ncolv%nplane_angle   = 0
    ncolv%nrot           = 0
    ncolv%nqparm         = 0
    ncolv%nxyz_diag      = 0
    ncolv%nxyz_outerdiag = 0
    ncolv%nhydronium     = 0
    ncolv%nreactionpath  = 0
    ncolv%ncombinecvs    = 0
    ncolv%nrestraint     = 0
    ncolv%npopulation    = 0
    ncolv%ngyration      = 0

    IF (ASSOCIATED(colv_list)) THEN
       DO k = 1, SIZE(colv_list)
          IF (colv_list(k)%restraint%active) ncolv%nrestraint  = ncolv%nrestraint + 1
          SELECT CASE(colv_list(k)%type_id)
          CASE(angle_colvar_id)
             ncolv%nangle = ncolv%nangle + 1
          CASE(coord_colvar_id)
             ncolv%ncoord = ncolv%ncoord + 1
          CASE(population_colvar_id)
             ncolv%npopulation = ncolv%npopulation + 1
          CASE(gyration_colvar_id)
             ncolv%ngyration = ncolv%ngyration + 1
          CASE(rotation_colvar_id)
             ncolv%nrot = ncolv%nrot + 1
          CASE(dist_colvar_id)
             ncolv%ndist = ncolv%ndist + 1
          CASE(dfunct_colvar_id)
             ncolv%ndfunct = ncolv%ndfunct + 1
          CASE(plane_distance_colvar_id)
             ncolv%nplane_dist = ncolv%nplane_dist + 1
          CASE(plane_plane_angle_colvar_id)
             ncolv%nplane_angle = ncolv%nplane_angle + 1
          CASE(torsion_colvar_id)
             ncolv%ntorsion = ncolv%ntorsion + 1
          CASE(qparm_colvar_id)
             ncolv%nqparm = ncolv%nqparm + 1
          CASE(xyz_diag_colvar_id)
             ncolv%nxyz_diag = ncolv%nxyz_diag + 1
          CASE(xyz_outerdiag_colvar_id)
             ncolv%nxyz_outerdiag = ncolv%nxyz_outerdiag + 1
          CASE(hydronium_colvar_id)
             ncolv%nhydronium = ncolv%nhydronium + 1
          CASE(reaction_path_colvar_id)
             ncolv%nreactionpath = ncolv%nreactionpath + 1
          CASE(combine_colvar_id)
             ncolv%ncombinecvs = ncolv%ncombinecvs + 1
          CASE DEFAULT
             CPASSERT(.FALSE.)
          END SELECT
       END DO
    END IF
    ncolv%ntot = ncolv%ndist           +&
                 ncolv%nangle          +&
                 ncolv%ntorsion        +&
                 ncolv%ncoord          +&
                 ncolv%nplane_dist     +&
                 ncolv%nplane_angle    +&
                 ncolv%ndfunct         +&
                 ncolv%nrot            +&
                 ncolv%nqparm          +&
                 ncolv%nxyz_diag       +&
                 ncolv%nxyz_outerdiag  +&
                 ncolv%nhydronium      +&
                 ncolv%nreactionpath   +&
                 ncolv%ncombinecvs     +&
                 ncolv%npopulation     +&
                 ncolv%ngyration


  END SUBROUTINE setup_colvar_counters

! *****************************************************************************
!> \brief   Allocate and initialize a molecule kind set.
!> \param molecule_kind_set ...
!> \param nmolecule_kind ...
!> \date    22.08.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_molecule_kind_set(molecule_kind_set,nmolecule_kind)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    INTEGER, INTENT(IN)                      :: nmolecule_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_molecule_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imolecule_kind
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(molecule_kind_set)) THEN
      CALL deallocate_molecule_kind_set(molecule_kind_set)
    END IF

    ALLOCATE (molecule_kind_set(nmolecule_kind))

    DO imolecule_kind=1,nmolecule_kind
      NULLIFY (molecule_kind_set(imolecule_kind)%atom_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bond_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bend_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%colv_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%ub_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%ub_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%impr_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%impr_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%opbend_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%opbend_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%g3x3_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%g4x6_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%vsite_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%fixd_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%shell_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%torsion_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bond_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%bend_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%torsion_kind_set)
      molecule_kind_set(imolecule_kind)%charge = 0.0_dp
      molecule_kind_set(imolecule_kind)%mass = 0.0_dp
      molecule_kind_set(imolecule_kind)%name = ""
      molecule_kind_set(imolecule_kind)%molname_generated = .FALSE.
      molecule_kind_set(imolecule_kind)%kind_number = imolecule_kind
      molecule_kind_set(imolecule_kind)%natom = 0
      molecule_kind_set(imolecule_kind)%nbend = 0
      molecule_kind_set(imolecule_kind)%nbond = 0
      molecule_kind_set(imolecule_kind)%nimpr = 0
      molecule_kind_set(imolecule_kind)%nopbend = 0
      molecule_kind_set(imolecule_kind)%nub = 0
      CALL setup_colvar_counters(molecule_kind_set(imolecule_kind)%colv_list,&
        molecule_kind_set(imolecule_kind)%ncolv)
      molecule_kind_set(imolecule_kind)%ng3x3 = 0
      molecule_kind_set(imolecule_kind)%ng4x6 = 0
      molecule_kind_set(imolecule_kind)%nvsite = 0
      molecule_kind_set(imolecule_kind)%nfixd = 0
      molecule_kind_set(imolecule_kind)%ng3x3_restraint = 0
      molecule_kind_set(imolecule_kind)%ng4x6_restraint = 0
      molecule_kind_set(imolecule_kind)%nvsite_restraint = 0
      molecule_kind_set(imolecule_kind)%nfixd_restraint = 0
      molecule_kind_set(imolecule_kind)%nmolecule = 0
      molecule_kind_set(imolecule_kind)%ntorsion = 0
      molecule_kind_set(imolecule_kind)%nshell = 0
      NULLIFY (molecule_kind_set(imolecule_kind)%molecule_list)
    END DO

  END SUBROUTINE allocate_molecule_kind_set

! *****************************************************************************
!> \brief   Deallocate a molecule kind set.
!> \param molecule_kind_set ...
!> \date    22.08.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_molecule_kind_set(molecule_kind_set)

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_molecule_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imolecule_kind, j, &
                                                nmolecule_kind
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(molecule_kind_set)) THEN

      nmolecule_kind = SIZE(molecule_kind_set)

      DO imolecule_kind=1,nmolecule_kind

        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%atom_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%atom_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bend_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bend_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bend_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bend_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%ub_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%ub_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%ub_kind_set)) THEN
           CALL ub_kind_dealloc_ref(molecule_kind_set(imolecule_kind)%ub_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%impr_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%impr_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%impr_kind_set)) THEN
           DO i = 1, SIZE(molecule_kind_set(imolecule_kind)%impr_kind_set)
              CALL impr_kind_dealloc_ref()  !This Subroutine doesn't deallocate anything, maybe needs to be implemented
           END DO
           DEALLOCATE (molecule_kind_set(imolecule_kind)%impr_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%opbend_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%opbend_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%opbend_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%opbend_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bond_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bond_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bond_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bond_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%colv_list)) THEN
           DO j = 1, SIZE(molecule_kind_set(imolecule_kind)%colv_list)
              DEALLOCATE (molecule_kind_set(imolecule_kind)%colv_list(j)%i_atoms)
           END DO
           DEALLOCATE (molecule_kind_set(imolecule_kind)%colv_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%g3x3_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%g3x3_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%g4x6_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%g4x6_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%vsite_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%vsite_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%fixd_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%fixd_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%torsion_kind_set)) THEN
           DO i = 1, SIZE(molecule_kind_set(imolecule_kind)%torsion_kind_set)
              CALL torsion_kind_dealloc_ref(molecule_kind_set(imolecule_kind)%torsion_kind_set(i))
           END DO
           DEALLOCATE (molecule_kind_set(imolecule_kind)%torsion_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%shell_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%shell_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%torsion_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%torsion_list)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%molecule_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%molecule_list)
        ENDIF
      END DO

      DEALLOCATE (molecule_kind_set)
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
           "The pointer molecule_kind_set is not associated and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_molecule_kind_set

! *****************************************************************************
!> \brief   Get informations about a molecule kind.
!> \param molecule_kind ...
!> \param atom_list ...
!> \param bond_list ...
!> \param bend_list ...
!> \param ub_list ...
!> \param impr_list ...
!> \param opbend_list ...
!> \param colv_list ...
!> \param fixd_list ...
!> \param g3x3_list ...
!> \param g4x6_list ...
!> \param vsite_list ...
!> \param torsion_list ...
!> \param shell_list ...
!> \param name ...
!> \param mass ...
!> \param charge ...
!> \param kind_number ...
!> \param natom ...
!> \param nbend ...
!> \param nbond ...
!> \param nub ...
!> \param nimpr ...
!> \param nopbend ...
!> \param nconstraint ...
!> \param nconstraint_fixd ...
!> \param nfixd ...
!> \param ncolv ...
!> \param ng3x3 ...
!> \param ng4x6 ...
!> \param nvsite ...
!> \param nfixd_restraint ...
!> \param ng3x3_restraint ...
!> \param ng4x6_restraint ...
!> \param nvsite_restraint ...
!> \param nrestraints ...
!> \param nmolecule ...
!> \param nsgf ...
!> \param nshell ...
!> \param ntorsion ...
!> \param molecule_list ...
!> \param nelectron ...
!> \param nelectron_alpha ...
!> \param nelectron_beta ...
!> \param bond_kind_set ...
!> \param bend_kind_set ...
!> \param ub_kind_set ...
!> \param impr_kind_set ...
!> \param opbend_kind_set ...
!> \param torsion_kind_set ...
!> \param molname_generated ...
!> \date    27.08.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_molecule_kind(molecule_kind,atom_list,bond_list,bend_list,&
                               ub_list,impr_list,opbend_list,colv_list,fixd_list,&
                               g3x3_list,g4x6_list,vsite_list,torsion_list,shell_list,&
                               name,mass,charge,kind_number,natom,nbend,nbond,nub,&
                               nimpr,nopbend,nconstraint,nconstraint_fixd,nfixd,ncolv,ng3x3,ng4x6,&
                               nvsite,nfixd_restraint,ng3x3_restraint, ng4x6_restraint,&
                               nvsite_restraint,nrestraints,nmolecule,nsgf,nshell,ntorsion,&
                               molecule_list,nelectron, nelectron_alpha,nelectron_beta, &
                               bond_kind_set,bend_kind_set,&
                               ub_kind_set,impr_kind_set,opbend_kind_set,torsion_kind_set,&
                               molname_generated)

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atom_list
    TYPE(bond_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_list
    TYPE(bend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_list
    TYPE(ub_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: ub_list
    TYPE(impr_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_list
    TYPE(opbend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: opbend_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: colv_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: fixd_list
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g4x6_list
    TYPE(vsite_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: vsite_list
    TYPE(torsion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_list
    TYPE(shell_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: shell_list
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), OPTIONAL                  :: mass, charge
    INTEGER, INTENT(OUT), OPTIONAL           :: kind_number, natom, nbend, &
                                                nbond, nub, nimpr, nopbend, &
                                                nconstraint, &
                                                nconstraint_fixd, nfixd
    TYPE(colvar_counters), INTENT(out), &
      OPTIONAL                               :: ncolv
    INTEGER, INTENT(OUT), OPTIONAL :: ng3x3, ng4x6, nvsite, nfixd_restraint, &
      ng3x3_restraint, ng4x6_restraint, nvsite_restraint, nrestraints, &
      nmolecule, nsgf, nshell, ntorsion
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: molecule_list
    INTEGER, INTENT(OUT), OPTIONAL           :: nelectron, nelectron_alpha, &
                                                nelectron_beta
    TYPE(bond_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_kind_set
    TYPE(bend_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_kind_set
    TYPE(ub_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ub_kind_set
    TYPE(impr_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_kind_set
    TYPE(opbend_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: opbend_kind_set
    TYPE(torsion_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_kind_set
    LOGICAL, INTENT(OUT), OPTIONAL           :: molname_generated

    CHARACTER(len=*), PARAMETER :: routineN = 'get_molecule_kind', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(atom_list)) atom_list => molecule_kind%atom_list
      IF (PRESENT(bend_list)) bend_list => molecule_kind%bend_list
      IF (PRESENT(bond_list)) bond_list => molecule_kind%bond_list
      IF (PRESENT(impr_list)) impr_list => molecule_kind%impr_list
      IF (PRESENT(opbend_list)) opbend_list => molecule_kind%opbend_list
      IF (PRESENT(ub_list)) ub_list => molecule_kind%ub_list
      IF (PRESENT(bond_kind_set)) bond_kind_set => molecule_kind%bond_kind_set
      IF (PRESENT(bend_kind_set)) bend_kind_set => molecule_kind%bend_kind_set
      IF (PRESENT(ub_kind_set)) ub_kind_set => molecule_kind%ub_kind_set
      IF (PRESENT(impr_kind_set)) impr_kind_set => molecule_kind%impr_kind_set
      IF (PRESENT(opbend_kind_set)) opbend_kind_set => molecule_kind%opbend_kind_set
      IF (PRESENT(torsion_kind_set)) torsion_kind_set => molecule_kind%torsion_kind_set
      IF (PRESENT(colv_list)) colv_list => molecule_kind%colv_list
      IF (PRESENT(g3x3_list)) g3x3_list => molecule_kind%g3x3_list
      IF (PRESENT(g4x6_list)) g4x6_list => molecule_kind%g4x6_list
      IF (PRESENT(vsite_list)) vsite_list => molecule_kind%vsite_list
      IF (PRESENT(fixd_list)) fixd_list => molecule_kind%fixd_list
      IF (PRESENT(torsion_list)) torsion_list => molecule_kind%torsion_list
      IF (PRESENT(shell_list)) shell_list => molecule_kind%shell_list
      IF (PRESENT(name)) name = molecule_kind%name
      IF (PRESENT(molname_generated)) molname_generated = molecule_kind%molname_generated
      IF (PRESENT(mass)) mass = molecule_kind%mass
      IF (PRESENT(charge)) charge = molecule_kind%charge
      IF (PRESENT(kind_number)) kind_number = molecule_kind%kind_number
      IF (PRESENT(natom)) natom = molecule_kind%natom
      IF (PRESENT(nbend)) nbend = molecule_kind%nbend
      IF (PRESENT(nbond)) nbond = molecule_kind%nbond
      IF (PRESENT(nub)) nub = molecule_kind%nub
      IF (PRESENT(nimpr)) nimpr = molecule_kind%nimpr
      IF (PRESENT(nopbend)) nopbend = molecule_kind%nopbend
      IF (PRESENT(nconstraint)) nconstraint = (molecule_kind%ncolv%ntot - molecule_kind%ncolv%nrestraint) +&
                                              3*(molecule_kind%ng3x3-molecule_kind%ng3x3_restraint) +&
                                              6*(molecule_kind%ng4x6-molecule_kind%ng4x6_restraint)+&
                                              3*(molecule_kind%nvsite-molecule_kind%nvsite_restraint)
      IF (PRESENT(ncolv)) ncolv = molecule_kind%ncolv
      IF (PRESENT(ng3x3)) ng3x3= molecule_kind%ng3x3
      IF (PRESENT(ng4x6)) ng4x6= molecule_kind%ng4x6
      IF (PRESENT(nvsite)) nvsite= molecule_kind%nvsite
      ! Number of atoms that have one or more components fixed
      IF (PRESENT(nfixd)) nfixd= molecule_kind%nfixd
      ! Number of degrees of freedom fixed
      IF (PRESENT(nconstraint_fixd)) THEN
         nconstraint_fixd = 0
         IF (molecule_kind%nfixd/=0) THEN
            DO i = 1, SIZE(molecule_kind%fixd_list)
               IF (molecule_kind%fixd_list(i)%restraint%active) CYCLE
               SELECT CASE(molecule_kind%fixd_list(i)%itype)
               CASE(use_perd_x, use_perd_y, use_perd_z)
                  nconstraint_fixd=nconstraint_fixd+1
               CASE(use_perd_xy, use_perd_xz, use_perd_yz)
                  nconstraint_fixd=nconstraint_fixd+2
               CASE(use_perd_xyz)
                  nconstraint_fixd=nconstraint_fixd+3
               END SELECT
            END DO
         END IF
      END IF
      IF (PRESENT(ng3x3_restraint)) ng3x3_restraint= molecule_kind%ng3x3_restraint
      IF (PRESENT(ng4x6_restraint)) ng4x6_restraint= molecule_kind%ng4x6_restraint
      IF (PRESENT(nvsite_restraint)) nvsite_restraint= molecule_kind%nvsite_restraint
      IF (PRESENT(nfixd_restraint)) nfixd_restraint= molecule_kind%nfixd_restraint
      IF (PRESENT(nrestraints)) nrestraints = molecule_kind%ncolv%nrestraint +&
                                              molecule_kind%ng3x3_restraint +&
                                              molecule_kind%ng4x6_restraint+&
                                              molecule_kind%nvsite_restraint
      IF (PRESENT(nmolecule)) nmolecule = molecule_kind%nmolecule
      IF (PRESENT(nshell)) nshell= molecule_kind%nshell
      IF (PRESENT(ntorsion)) ntorsion= molecule_kind%ntorsion
      IF (PRESENT(nsgf)) nsgf = molecule_kind%nsgf
      IF (PRESENT(nelectron)) nelectron = molecule_kind%nelectron
      IF (PRESENT(nelectron_alpha)) nelectron_alpha = molecule_kind%nelectron_beta
      IF (PRESENT(nelectron_beta)) nelectron_beta = molecule_kind%nelectron_alpha
      IF (PRESENT(molecule_list)) molecule_list => molecule_kind%molecule_list

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE get_molecule_kind

! *****************************************************************************
!> \brief   Get informations about a molecule kind set.
!> \param molecule_kind_set ...
!> \param maxatom ...
!> \param natom ...
!> \param nbond ...
!> \param nbend ...
!> \param nub ...
!> \param ntorsion ...
!> \param nimpr ...
!> \param nopbend ...
!> \param nconstraint ...
!> \param nconstraint_fixd ...
!> \param nmolecule ...
!> \param nrestraints ...
!> \date    27.08.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_molecule_kind_set(molecule_kind_set,maxatom,natom,&
                                   nbond,nbend,nub,ntorsion,nimpr,nopbend,&
                                   nconstraint,nconstraint_fixd,nmolecule,&
                                   nrestraints)

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    INTEGER, INTENT(OUT), OPTIONAL :: maxatom, natom, nbond, nbend, nub, &
      ntorsion, nimpr, nopbend, nconstraint, nconstraint_fixd, nmolecule, &
      nrestraints

    CHARACTER(len=*), PARAMETER :: routineN = 'get_molecule_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER :: ibend, ibond, iimpr, imolecule_kind, iopbend, itorsion, iub, &
      na, nc, nc_fixd, nfixd_restraint, nm, nmolecule_kind, nrestraints_tot
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    IF (ASSOCIATED(molecule_kind_set)) THEN

      IF (PRESENT(maxatom)) maxatom = 0
      IF (PRESENT(natom)) natom = 0
      IF (PRESENT(nbond)) nbond = 0
      IF (PRESENT(nbend)) nbend = 0
      IF (PRESENT(nub)) nub = 0
      IF (PRESENT(ntorsion)) ntorsion = 0
      IF (PRESENT(nimpr)) nimpr = 0
      IF (PRESENT(nopbend)) nopbend = 0
      IF (PRESENT(nconstraint)) nconstraint = 0
      IF (PRESENT(nconstraint_fixd)) nconstraint_fixd = 0
      IF (PRESENT(nrestraints)) nrestraints = 0
      IF (PRESENT(nmolecule)) nmolecule = 0

      nmolecule_kind = SIZE(molecule_kind_set)

      DO imolecule_kind=1,nmolecule_kind

        molecule_kind => molecule_kind_set(imolecule_kind)

        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               natom=na,&
                               nbond=ibond,&
                               nbend=ibend,&
                               nub=iub,&
                               ntorsion=itorsion,&
                               nimpr=iimpr,&
                               nopbend=iopbend,&
                               nconstraint=nc,&
                               nconstraint_fixd=nc_fixd,&
                               nfixd_restraint=nfixd_restraint,&
                               nrestraints=nrestraints_tot,&
                               nmolecule=nm)
        IF (PRESENT(maxatom)) maxatom = MAX(maxatom,na)
        IF (PRESENT(natom)) natom = natom + na*nm
        IF (PRESENT(nbond)) nbond = nbond + ibond*nm
        IF (PRESENT(nbend)) nbend = nbend + ibend*nm
        IF (PRESENT(nub)) nub = nub + iub*nm
        IF (PRESENT(ntorsion)) ntorsion = ntorsion + itorsion*nm
        IF (PRESENT(nimpr)) nimpr = nimpr + iimpr*nm
        IF (PRESENT(nopbend)) nopbend = nopbend + iopbend*nm
        IF (PRESENT(nconstraint)) nconstraint = nconstraint + nc*nm + nc_fixd
        IF (PRESENT(nconstraint_fixd)) nconstraint_fixd = nconstraint_fixd + nc_fixd
        IF (PRESENT(nmolecule)) nmolecule = nmolecule + nm
        IF (PRESENT(nrestraints)) nrestraints = nrestraints + nm*nrestraints_tot + nfixd_restraint

      END DO

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind_set is not associated")

    END IF

  END SUBROUTINE get_molecule_kind_set

! *****************************************************************************
!> \brief   Set the components of a molecule kind.
!> \param molecule_kind ...
!> \param name ...
!> \param mass ...
!> \param charge ...
!> \param kind_number ...
!> \param molecule_list ...
!> \param atom_list ...
!> \param nbond ...
!> \param bond_list ...
!> \param nbend ...
!> \param bend_list ...
!> \param nub ...
!> \param ub_list ...
!> \param nimpr ...
!> \param impr_list ...
!> \param nopbend ...
!> \param opbend_list ...
!> \param ntorsion ...
!> \param torsion_list ...
!> \param fixd_list ...
!> \param ncolv ...
!> \param colv_list ...
!> \param ng3x3 ...
!> \param g3x3_list ...
!> \param ng4x6 ...
!> \param nfixd ...
!> \param g4x6_list ...
!> \param nvsite ...
!> \param vsite_list ...
!> \param ng3x3_restraint ...
!> \param ng4x6_restraint ...
!> \param nfixd_restraint ...
!> \param nshell ...
!> \param shell_list ...
!> \param nvsite_restraint ...
!> \param bond_kind_set ...
!> \param bend_kind_set ...
!> \param ub_kind_set ...
!> \param torsion_kind_set ...
!> \param impr_kind_set ...
!> \param opbend_kind_set ...
!> \param nelectron ...
!> \param nsgf ...
!> \param molname_generated ...
!> \date    27.08.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
 SUBROUTINE set_molecule_kind(molecule_kind,name,mass,charge,kind_number,&
                               molecule_list,atom_list,nbond,bond_list,&
                               nbend,bend_list,nub,ub_list,nimpr,impr_list,&
                               nopbend, opbend_list, ntorsion,&
                               torsion_list,fixd_list,ncolv,colv_list,ng3x3,&
                               g3x3_list,ng4x6,nfixd,g4x6_list,nvsite,&
                               vsite_list,ng3x3_restraint,ng4x6_restraint,&
                               nfixd_restraint,nshell, shell_list,&
                               nvsite_restraint,bond_kind_set,bend_kind_set,&
                               ub_kind_set,torsion_kind_set,impr_kind_set,&
                               opbend_kind_set, nelectron,nsgf,&
                               molname_generated)

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    REAL(KIND=dp), OPTIONAL                  :: mass, charge
    INTEGER, INTENT(IN), OPTIONAL            :: kind_number
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: molecule_list
    TYPE(atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atom_list
    INTEGER, INTENT(IN), OPTIONAL            :: nbond
    TYPE(bond_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_list
    INTEGER, INTENT(IN), OPTIONAL            :: nbend
    TYPE(bend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_list
    INTEGER, INTENT(IN), OPTIONAL            :: nub
    TYPE(ub_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: ub_list
    INTEGER, INTENT(IN), OPTIONAL            :: nimpr
    TYPE(impr_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_list
    INTEGER, INTENT(IN), OPTIONAL            :: nopbend
    TYPE(opbend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: opbend_list
    INTEGER, INTENT(IN), OPTIONAL            :: ntorsion
    TYPE(torsion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: fixd_list
    TYPE(colvar_counters), INTENT(IN), &
      OPTIONAL                               :: ncolv
    TYPE(colvar_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: colv_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng3x3
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g3x3_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng4x6, nfixd
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g4x6_list
    INTEGER, INTENT(IN), OPTIONAL            :: nvsite
    TYPE(vsite_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: vsite_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng3x3_restraint, &
                                                ng4x6_restraint, &
                                                nfixd_restraint, nshell
    TYPE(shell_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: shell_list
    INTEGER, INTENT(IN), OPTIONAL            :: nvsite_restraint
    TYPE(bond_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_kind_set
    TYPE(bend_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_kind_set
    TYPE(ub_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ub_kind_set
    TYPE(torsion_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_kind_set
    TYPE(impr_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_kind_set
    TYPE(opbend_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: opbend_kind_set
    INTEGER, INTENT(IN), OPTIONAL            :: nelectron, nsgf
    LOGICAL, INTENT(IN), OPTIONAL            :: molname_generated

    CHARACTER(len=*), PARAMETER :: routineN = 'set_molecule_kind', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(atom_list)) THEN
        n = SIZE(atom_list)
        molecule_kind%natom = n
        molecule_kind%atom_list => atom_list
      END IF
      IF (PRESENT(molname_generated)) molecule_kind%molname_generated = molname_generated
      IF (PRESENT(name)) molecule_kind%name = name
      IF (PRESENT(mass)) molecule_kind%mass = mass
      IF (PRESENT(charge)) molecule_kind%charge = charge
      IF (PRESENT(kind_number)) molecule_kind%kind_number = kind_number
      IF (PRESENT(nbond)) molecule_kind%nbond = nbond
      IF (PRESENT(bond_list)) molecule_kind%bond_list => bond_list
      IF (PRESENT(nbend)) molecule_kind%nbend = nbend
      IF (PRESENT(nelectron)) molecule_kind%nelectron = nelectron
      IF (PRESENT(nsgf)) molecule_kind%nsgf = nsgf
      IF (PRESENT(bend_list)) molecule_kind%bend_list => bend_list
      IF (PRESENT(nub)) molecule_kind%nub = nub
      IF (PRESENT(ub_list)) molecule_kind%ub_list => ub_list
      IF (PRESENT(ntorsion)) molecule_kind%ntorsion = ntorsion
      IF (PRESENT(torsion_list)) molecule_kind%torsion_list => torsion_list
      IF (PRESENT(nimpr)) molecule_kind%nimpr = nimpr
      IF (PRESENT(impr_list)) molecule_kind%impr_list => impr_list
      IF (PRESENT(nopbend)) molecule_kind%nopbend = nopbend
      IF (PRESENT(opbend_list)) molecule_kind%opbend_list => opbend_list
      IF (PRESENT(ncolv)) molecule_kind%ncolv = ncolv
      IF (PRESENT(colv_list)) molecule_kind%colv_list => colv_list
      IF (PRESENT(ng3x3)) molecule_kind%ng3x3 = ng3x3
      IF (PRESENT(g3x3_list)) molecule_kind%g3x3_list => g3x3_list
      IF (PRESENT(ng4x6)) molecule_kind%ng4x6 = ng4x6
      IF (PRESENT(nvsite)) molecule_kind%nvsite = nvsite
      IF (PRESENT(nfixd)) molecule_kind%nfixd = nfixd
      IF (PRESENT(nfixd_restraint)) molecule_kind%nfixd_restraint = nfixd_restraint
      IF (PRESENT(ng3x3_restraint)) molecule_kind%ng3x3_restraint = ng3x3_restraint
      IF (PRESENT(ng4x6_restraint)) molecule_kind%ng4x6_restraint = ng4x6_restraint
      IF (PRESENT(nvsite_restraint)) molecule_kind%nvsite_restraint = nvsite_restraint
      IF (PRESENT(g4x6_list)) molecule_kind%g4x6_list => g4x6_list
      IF (PRESENT(vsite_list)) molecule_kind%vsite_list => vsite_list
      IF (PRESENT(fixd_list)) molecule_kind%fixd_list => fixd_list
      IF (PRESENT(bond_kind_set)) molecule_kind%bond_kind_set => bond_kind_set
      IF (PRESENT(bend_kind_set)) molecule_kind%bend_kind_set => bend_kind_set
      IF (PRESENT(ub_kind_set)) molecule_kind%ub_kind_set => ub_kind_set
      IF (PRESENT(torsion_kind_set)) molecule_kind%torsion_kind_set => torsion_kind_set
      IF (PRESENT(impr_kind_set)) molecule_kind%impr_kind_set => impr_kind_set
      IF (PRESENT(opbend_kind_set)) molecule_kind%opbend_kind_set => opbend_kind_set
      IF (PRESENT(nshell)) molecule_kind%nshell = nshell
      IF (PRESENT(shell_list)) molecule_kind%shell_list => shell_list
      IF (PRESENT(molecule_list)) THEN
        n = SIZE(molecule_list)
        molecule_kind%nmolecule = n
        molecule_kind%molecule_list => molecule_list
      END IF

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE set_molecule_kind

! *****************************************************************************
!> \brief   Write a molecule kind data set to the output unit.
!> \param molecule_kind ...
!> \param output_unit ...
!> \date    24.09.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_molecule_kind(molecule_kind,output_unit)
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    INTEGER, INTENT(in)                      :: output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'write_molecule_kind', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: iatom, imolecule, natom, &
                                                nmolecule
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure = .FALSE.
    IF (output_unit>0) THEN
      IF (ASSOCIATED(molecule_kind)) THEN
        natom = SIZE(molecule_kind%atom_list)
        nmolecule = SIZE(molecule_kind%molecule_list)

        IF ( natom == 1 ) THEN
          atomic_kind => molecule_kind%atom_list(1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
          WRITE (UNIT=output_unit,FMT="(/,T2,I5,A,T36,A,A,T64,A)")&
            molecule_kind%kind_number,&
            ". Molecule kind: "//TRIM(molecule_kind%name),&
            "Atomic kind name:   ",TRIM(name)
          WRITE (UNIT=output_unit,FMT="(T9,A,L1,T55,A,T75,I6)")&
             "Automatic name: ",molecule_kind%molname_generated,&
             "Number of molecules:",nmolecule
        ELSE
          WRITE (UNIT=output_unit,FMT="(/,T2,I5,A,T50,A,T75,I6,/,T22,A)")&
            molecule_kind%kind_number,&
            ". Molecule kind: "//TRIM(molecule_kind%name),&
            "Number of atoms:    ",natom,&
            "Atom         Atomic kind name"
          DO iatom=1,natom
            atomic_kind => molecule_kind%atom_list(iatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
            WRITE (UNIT=output_unit,FMT="(T20,I6,(7X,A18))")&
              iatom,TRIM(name)
          END DO
          WRITE (UNIT=output_unit,FMT="(/,T9,A,L1)")&
             "The name was automatically generated: ",&
             molecule_kind%molname_generated
          WRITE (UNIT=output_unit,FMT="(T9,A,I6,/,T9,A,(T30,5I10))")&
            "Number of molecules: ",nmolecule,"Molecule list:",&
            (molecule_kind%molecule_list(imolecule),imolecule=1,nmolecule)
          IF ( molecule_kind%nbond > 0 ) &
            WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
            "Number of bonds:       ",molecule_kind%nbond
          IF ( molecule_kind%nbend > 0 ) &
            WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
            "Number of bends:       ",molecule_kind%nbend
          IF ( molecule_kind%nub > 0 ) &
            WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
            "Number of Urey-Bradley:",molecule_kind%nub
          IF ( molecule_kind%ntorsion > 0 ) &
            WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
            "Number of torsions:    ",molecule_kind%ntorsion
          IF ( molecule_kind%nimpr > 0 ) &
            WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
            "Number of improper:    ",molecule_kind%nimpr
          IF ( molecule_kind%nopbend > 0 ) &
            WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
            "Number of out opbends:    ",molecule_kind%nopbend
        END IF

      ELSE
          CPASSERT(.FALSE.)
      END IF

    END IF
  END SUBROUTINE write_molecule_kind

! *****************************************************************************
!> \brief   Write a moleculeatomic kind set data set to the output unit.
!> \param molecule_kind_set ...
!> \param subsys_section ...
!> \date    24.09.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_molecule_kind_set(molecule_kind_set,subsys_section)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(len=*), PARAMETER :: routineN = 'write_molecule_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, imolecule_kind, natom, nbend, nbond, nimpr, nmolecule, &
      nmolecule_kind, nopbend, ntors, ntotal, nub, output_unit
    LOGICAL                                  :: all_single_atoms, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_get_default_logger()
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%MOLECULES",extension=".Log")
    IF (output_unit>0) THEN
      IF (ASSOCIATED(molecule_kind_set)) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "MOLECULE KIND INFORMATION"

        nmolecule_kind = SIZE(molecule_kind_set)

        all_single_atoms=.TRUE.
        DO imolecule_kind=1,nmolecule_kind
           natom = SIZE(molecule_kind_set(imolecule_kind)%atom_list)
           nmolecule = SIZE(molecule_kind_set(imolecule_kind)%molecule_list)
           IF (natom*nmolecule>1) all_single_atoms=.FALSE.
        ENDDO

        IF (all_single_atoms) THEN
           WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") &
              "All atoms are their own molecule, skipping detailed information"
        ELSE
           DO imolecule_kind=1,nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL write_molecule_kind(molecule_kind,output_unit)
           END DO
        ENDIF

        CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
                                   nbond=nbond,&
                                   nbend=nbend,&
                                   nub=nub,&
                                   ntorsion=ntors,&
                                   nimpr=nimpr,&
                                   nopbend=nopbend)
        ntotal = nbond+nbend+nub+ntors+nimpr+nopbend
        IF ( ntotal > 0 ) THEN
           WRITE(UNIT=output_unit,FMT="(/,/,T2,A,T45,A30,I6)") &
             "MOLECULE KIND SET INFORMATION", &
             "Total Number of bonds:       ",nbond
           WRITE (UNIT=output_unit,FMT="(T45,A30,I6)")&
             "Total Number of bends:       ",nbend
           WRITE (UNIT=output_unit,FMT="(T45,A30,I6)")&
             "Total Number of Urey-Bradley:",nub
           WRITE (UNIT=output_unit,FMT="(T45,A30,I6)")&
             "Total Number of torsions:    ",ntors
           WRITE (UNIT=output_unit,FMT="(T45,A30,I6)")&
             "Total Number of improper:    ",nimpr
           WRITE (UNIT=output_unit,FMT="(T45,A30,I6)")&
             "Total Number of opbends:    ",nopbend
        END IF
      ELSE
         CPASSERT(.FALSE.)
      END IF

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%MOLECULES")

    CALL timestop(handle)

  END SUBROUTINE write_molecule_kind_set

END MODULE molecule_kind_types
