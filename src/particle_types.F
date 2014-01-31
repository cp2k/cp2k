!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the data structure for the particle information.
!> \par History
!>      - Atomic kind added in particle_type (MK,08.01.2002)
!>      - Functionality for particle_type added (MK,14.01.2002)
!>      - Allow for general coordinate input (MK,13.09.2003)
!>      - Molecule concept introduced (MK,26.09.2003)
!>      - Last atom information added (jgh,23.05.2004)
!>      - particle_type cleaned (MK,03.02.2005)
!> \author CJM, MK
! *****************************************************************************
MODULE particle_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             set_cell_param
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE external_potential_types,        ONLY: fist_potential_type,&
                                             get_potential
  USE f77_blas
  USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_dcd_aligned_cell,&
                                             dump_pdb,&
                                             dump_xmol
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             sp,&
                                             sp_size
  USE mathconstants,                   ONLY: degree
  USE mathlib,                         ONLY: angle,&
                                             dihedral_angle
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE physcon,                         ONLY: massunit
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE shell_potential_types,           ONLY: get_shell,&
                                             shell_kind_type
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort,&
                                             sort_unique
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'particle_types'

  ! Data types
! *****************************************************************************
  TYPE particle_type
     TYPE(atomic_kind_type), POINTER       :: atomic_kind ! atomic kind information
     REAL(KIND = dp), DIMENSION(3)         :: f,&         ! force
                                              r,&         ! position
                                              v           ! velocity
     ! Particle dependent terms for shell-model
     INTEGER                               :: atom_index, t_region_index, shell_index
  END TYPE particle_type

  ! Public data types

  PUBLIC :: particle_type

  ! Public subroutines

  PUBLIC :: allocate_particle_set,&
            deallocate_particle_set,&
            get_particle_set,&
            write_fist_particle_coordinates,&
            write_qs_particle_coordinates,&
            write_particle_distances,&
            write_particle_coordinates,&
            write_structure_data,&
            write_particle_matrix,&
            update_particle_set,&
            update_particle_pos_or_vel

  PUBLIC :: get_particle_pos_or_vel

CONTAINS

! *****************************************************************************
!> \brief   Allocate a particle set.
!> \author  MK
!> \date    14.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_particle_set(particle_set,nparticle,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: nparticle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(particle_set)) THEN
       CALL deallocate_particle_set(particle_set,error=error)
    END IF
    ALLOCATE (particle_set(nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iparticle=1,nparticle
      NULLIFY (particle_set(iparticle)%atomic_kind)
      particle_set(iparticle)%f(:) = 0.0_dp
      particle_set(iparticle)%r(:) = 0.0_dp
      particle_set(iparticle)%v(:) = 0.0_dp
      particle_set(iparticle)%shell_index = 0
      particle_set(iparticle)%atom_index = 0
      particle_set(iparticle)%t_region_index = 0
    END DO

  END SUBROUTINE allocate_particle_set

! *****************************************************************************
!> \brief   Deallocate a particle set.
!> \author  MK
!> \date    14.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_particle_set(particle_set,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(particle_set)) THEN
      DEALLOCATE (particle_set,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer particle_set is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_particle_set

! *****************************************************************************
!> \brief   Get the components of a particle set.
!> \author  MK
!> \date    14.01.2002
!> \par History
!>      - particle type cleaned (13.10.2003,MK)
!>      - refactoring and add besis set option (17.08.2010,jhu)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_particle_set(particle_set,first_sgf,last_sgf,nsgf,&
                              basis_set_id,basis,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, DIMENSION(:), INTENT(INOUT), &
      OPTIONAL                               :: first_sgf, last_sgf, nsgf
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, iparticle, isgf, &
                                                nparticle, ns
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)

    nparticle = SIZE(particle_set)
    IF (PRESENT(first_sgf)) THEN
      CPPrecondition(SIZE(first_sgf) >= nparticle,cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(last_sgf)) THEN
      CPPrecondition(SIZE(last_sgf) >= nparticle,cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(nsgf)) THEN
      CPPrecondition(SIZE(nsgf) >= nparticle,cp_failure_level,routineP,error,failure)
    END IF

    IF (PRESENT(first_sgf).OR.PRESENT(last_sgf).OR.PRESENT(nsgf)) THEN
      isgf = 0
      DO iparticle=1,nparticle
        IF (PRESENT(basis)) THEN
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,kind_number=ikind)
          CALL get_gto_basis_set(gto_basis_set=basis(ikind)%gto_basis_set,nsgf=ns)
        ELSEIF (PRESENT(basis_set_id)) THEN
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,nsgf=ns,&
                             basis_set_id=basis_set_id)
        ELSE
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,nsgf=ns)
        END IF
        IF (PRESENT(nsgf)) nsgf(iparticle) = ns
        IF (PRESENT(first_sgf)) first_sgf(iparticle) = isgf + 1
        isgf = isgf + ns
        IF (PRESENT(last_sgf)) last_sgf(iparticle) = isgf
      END DO
    END IF
    IF (PRESENT(first_sgf)) THEN
      IF (SIZE(first_sgf) > nparticle) first_sgf(nparticle+1) = isgf + 1
    END IF

  END SUBROUTINE get_particle_set

! *****************************************************************************
!> \brief   Should be able to write a few formats e.g. xmol, and some binary
!>          format (dcd) some format can be used for x,v,f
!>
!>          FORMAT   CONTENT                                    UNITS x, v, f
!>          XMOL     POS, VEL, FORCE, POS_VEL, POS_VEL_FORCE    Angstrom, a.u., a.u.
!>
!> \author  MK
!> \date    14.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_particle_coordinates(particle_set,iunit,output_format,&
                                        content,title,cell,array,unit_conv,&
                                        charge_occup,charge_beta,&
                                        charge_extended,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: iunit, output_format
    CHARACTER(LEN=*)                         :: content, title
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: array
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: unit_conv
    LOGICAL, INTENT(IN), OPTIONAL            :: charge_occup, charge_beta, &
                                                charge_extended
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_particle_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=120)                       :: line
    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=4)                         :: name
    CHARACTER(LEN=default_string_length)     :: atm_name, my_format
    INTEGER                                  :: handle, iatom, natom, stat
    LOGICAL                                  :: dummy, failure, &
                                                my_charge_beta, &
                                                my_charge_extended, &
                                                my_charge_occup
    REAL(KIND=dp)                            :: angle_alpha, angle_beta, &
                                                angle_gamma, factor, qeff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: arr
    REAL(KIND=dp), DIMENSION(3)              :: abc, angles, f, r, v
    REAL(KIND=dp), DIMENSION(3, 3)           :: h
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:) :: x4, y4, z4
    TYPE(cell_type), POINTER                 :: cell_dcd
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(shell_kind_type), POINTER           :: shell

    CALL timeset(routineN,handle)

    failure = .FALSE.
    natom = SIZE(particle_set)
    IF (PRESENT(array)) THEN
       SELECT CASE(TRIM(content))
       CASE ("POS_VEL","POS_VEL_FORCE")
          CALL stop_program(routineN,moduleN,__LINE__,"Illegal usage")
       END SELECT
    END IF
    factor = 1.0_dp
    IF (PRESENT(unit_conv)) THEN
       factor = unit_conv
    END IF
    SELECT CASE (output_format)
    CASE (dump_xmol)
       WRITE(iunit,"(I8)") natom
       WRITE(iunit,"(A)")  TRIM(title)
       DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
               element_symbol=element_symbol)
          IF (LEN_TRIM(element_symbol) == 0)  THEN
             CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                  name=atm_name)
             dummy = qmmm_ff_precond_only_qm(id1=atm_name)
             my_format="(A4,"
             name = TRIM(atm_name)
          ELSE
             my_format="(T2,A2,"
             name = TRIM(element_symbol)
          END IF
          SELECT CASE (TRIM(content))
          CASE ("POS")
             IF (PRESENT(array)) THEN
                r(1:3) = get_particle_pos_or_vel(iatom,particle_set,array)
             ELSE
                r(:) = particle_set(iatom)%r(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),r(1:3)*factor
          CASE ("VEL")
             IF (PRESENT(array)) THEN
                v(1:3) = get_particle_pos_or_vel(iatom,particle_set,array)
             ELSE
                v(:) = particle_set(iatom)%v(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),v(1:3)*factor
          CASE ("FORCE")
             IF (PRESENT(array)) THEN
                f(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                f(:) = particle_set(iatom)%f(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),f(1:3)*factor
          CASE ("FORCE_MIXING_LABELS")
             IF (PRESENT(array)) THEN
                f(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                f(:) = particle_set(iatom)%f(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),f(1:3)*factor
          END SELECT
       END DO
    CASE (dump_atomic)
       DO iatom=1,natom
          SELECT CASE (TRIM(content))
          CASE ("POS")
             IF (PRESENT(array)) THEN
                r(1:3) = get_particle_pos_or_vel(iatom,particle_set,array)
             ELSE
                r(:) = particle_set(iatom)%r(:)
             END IF
             WRITE (iunit,"(3F20.10)")r(1:3)*factor
          CASE ("VEL")
             IF (PRESENT(array)) THEN
                v(1:3) = get_particle_pos_or_vel(iatom,particle_set,array)
             ELSE
                v(:) = particle_set(iatom)%v(:)
             END IF
             WRITE (iunit,"(3F20.10)")v(1:3)*factor
          CASE ("FORCE")
             IF (PRESENT(array)) THEN
                f(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                f(:) = particle_set(iatom)%f(:)
             END IF
             WRITE (iunit,"(3F20.10)")f(1:3)*factor
          CASE ("FORCE_MIXING_LABELS")
             IF (PRESENT(array)) THEN
                f(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                f(:) = particle_set(iatom)%f(:)
             END IF
             WRITE (iunit,"(3F20.10)")f(1:3)*factor
          END SELECT
       END DO
    CASE (dump_dcd,dump_dcd_aligned_cell)
       IF (.NOT.(PRESENT(cell))) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Cell is not present! Report this bug!")
       END IF
       CALL get_cell(cell,alpha=angle_alpha,beta=angle_beta,gamma=angle_gamma,&
                     abc=abc)
       IF (.NOT.cell%orthorhombic.AND.(output_format == dump_dcd_aligned_cell)) THEN
          ! In the case of a non-orthorhombic cell adopt a common convention
          ! for the orientation of the cell with respect to the Cartesian axes:
          ! Cell vector a is aligned with the x axis and the cell vector b lies
          ! in the xy plane.
          NULLIFY (cell_dcd)
          CALL cell_create(cell_dcd,error=error)
          CALL cell_clone(cell,cell_dcd,error=error)
          angles(1) = angle_alpha/degree
          angles(2) = angle_beta/degree
          angles(3) = angle_gamma/degree
          CALL set_cell_param(cell_dcd,abc,angles,&
                              do_init_cell=.TRUE.,error=error)
          h(1:3,1:3) = MATMUL(cell_dcd%hmat(1:3,1:3),cell%h_inv(1:3,1:3))
          CALL cell_release(cell_dcd,error=error)
       END IF
       ALLOCATE (arr(3,natom),STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "arr",3*natom*dp_size)
       IF (PRESENT(array)) THEN
          arr(1:3,1:natom) = RESHAPE(array,(/3,natom/))
       ELSE
          SELECT CASE (TRIM(content))
          CASE ("POS")
             DO iatom=1,natom
                arr(1:3,iatom) = particle_set(iatom)%r(1:3)
             END DO
          CASE ("VEL")
             DO iatom=1,natom
                arr(1:3,iatom) = particle_set(iatom)%v(1:3)
             END DO
          CASE ("FORCE")
             DO iatom=1,natom
                arr(1:3,iatom) = particle_set(iatom)%f(1:3)
             END DO
          CASE DEFAULT
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "Illegal DCD dump type")
          END SELECT
       END IF
       ALLOCATE (x4(natom),STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "x4",natom*sp_size)
       ALLOCATE (y4(natom),STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "y4",natom*sp_size)
       ALLOCATE (z4(natom),STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "z4",natom*sp_size)
       IF (.NOT.cell%orthorhombic.AND.(output_format == dump_dcd_aligned_cell)) THEN
          x4(1:natom) = MATMUL(h(1,1:3),arr(1:3,1:natom))
          y4(1:natom) = MATMUL(h(2,1:3),arr(1:3,1:natom))
          z4(1:natom) = MATMUL(h(3,1:3),arr(1:3,1:natom))
       ELSE
          x4(1:natom) = arr(1,1:natom)
          y4(1:natom) = arr(2,1:natom)
          z4(1:natom) = arr(3,1:natom)
       END IF
       WRITE (iunit) abc(1)*factor,angle_gamma,abc(2)*factor,&
                     angle_beta,angle_alpha,abc(3)*factor
       WRITE (iunit) x4*REAL(factor,KIND=sp)
       WRITE (iunit) y4*REAL(factor,KIND=sp)
       WRITE (iunit) z4*REAL(factor,KIND=sp)
       ! Release work storage
       DEALLOCATE (arr,STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"arr")
       DEALLOCATE (x4,STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"x4")
       DEALLOCATE (y4,STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"y4")
       DEALLOCATE (z4,STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"z4")
    CASE (dump_pdb)
       my_charge_occup = .FALSE.
       IF (PRESENT(charge_occup)) my_charge_occup = charge_occup
       my_charge_beta = .FALSE.
       IF (PRESENT(charge_beta)) my_charge_beta = charge_beta
       my_charge_extended = .FALSE.
       IF (PRESENT(charge_extended)) my_charge_extended = charge_extended
       IF (LEN_TRIM(title) > 0) THEN
          WRITE (UNIT=iunit,FMT="(A6,T11,A)")&
           "REMARK",TRIM(title)
       END IF
       CALL get_cell(cell,alpha=angle_alpha,beta=angle_beta,gamma=angle_gamma,abc=abc)
       ! COLUMNS       DATA TYPE      CONTENTS
       ! --------------------------------------------------
       !  1 -  6       Record name    "CRYST1"
       !  7 - 15       Real(9.3)      a (Angstroms)
       ! 16 - 24       Real(9.3)      b (Angstroms)
       ! 25 - 33       Real(9.3)      c (Angstroms)
       ! 34 - 40       Real(7.2)      alpha (degrees)
       ! 41 - 47       Real(7.2)      beta (degrees)
       ! 48 - 54       Real(7.2)      gamma (degrees)
       ! 56 - 66       LString        Space group
       ! 67 - 70       Integer        Z value
       WRITE (UNIT=iunit,FMT="(A6,3F9.3,3F7.2)")&
        "CRYST1",abc(1:3)*factor,angle_alpha,angle_beta,angle_gamma
       WRITE (UNIT=line(1:6),FMT="(A6)") "ATOM  "
       DO iatom=1,natom
          line = ""
          ! COLUMNS        DATA TYPE       CONTENTS
          !  1 -  6        Record name     "ATOM  "
          !  7 - 11        Integer         Atom serial number
          ! 13 - 16        Atom            Atom name
          ! 17             Character       Alternate location indicator
          ! 18 - 20        Residue name    Residue name
          ! 22             Character       Chain identifier
          ! 23 - 26        Integer         Residue sequence number
          ! 27             AChar           Code for insertion of residues
          ! 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstrom
          ! 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstrom
          ! 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstrom
          ! 55 - 60        Real(6.2)       Occupancy
          ! 61 - 66        Real(6.2)       Temperature factor (Default = 0.0)
          ! 73 - 76        LString(4)      Segment identifier, left-justified
          ! 77 - 78        LString(2)      Element symbol, right-justified
          ! 79 - 80        LString(2)      Charge on the atom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,name=atm_name,&
                               fist_potential=fist_potential,shell=shell)
          IF (LEN_TRIM(element_symbol) == 0)  THEN
             dummy = qmmm_ff_precond_only_qm(id1=atm_name)
          END IF
          name = TRIM(atm_name)
          IF (ASSOCIATED(fist_potential)) THEN
             CALL get_potential(potential=fist_potential,qeff=qeff)
          ELSE
             qeff = 0.0_dp
          END IF
          IF (ASSOCIATED(shell)) CALL get_shell(shell=shell,charge=qeff,error=error)
          WRITE (UNIT=line( 1: 6),FMT="(A6)") "ATOM  "
          WRITE (UNIT=line( 7:11),FMT="(I5)") MODULO(iatom,100000)
          WRITE (UNIT=line(13:16),FMT="(A4)") ADJUSTL(name)
          ! WRITE (UNIT=line(18:20),FMT="(A3)") TRIM(resname)
          ! WRITE (UNIT=line(23:26),FMT="(I4)") MODULO(idres,10000)
          SELECT CASE (TRIM(content))
          CASE ("POS")
             IF (PRESENT(array)) THEN
                r(1:3) = get_particle_pos_or_vel(iatom,particle_set,array)
             ELSE
                r(:) = particle_set(iatom)%r(:)
             END IF
             WRITE (UNIT=line(31:54),FMT="(3F8.3)") r(1:3)*factor
          CASE DEFAULT
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "PDB dump only for trajectory available")
          END SELECT
          IF (my_charge_occup) THEN
             WRITE (UNIT=line(55:60),FMT="(F6.2)") qeff
          ELSE
             WRITE (UNIT=line(55:60),FMT="(F6.2)") 0.0_dp
          END IF
          IF (my_charge_beta) THEN
             WRITE (UNIT=line(61:66),FMT="(F6.2)") qeff
          ELSE
             WRITE (UNIT=line(61:66),FMT="(F6.2)") 0.0_dp
          END IF
          ! WRITE (UNIT=line(73:76),FMT="(A4)") ADJUSTL(TRIM(molname))
          WRITE (UNIT=line(77:78),FMT="(A2)") ADJUSTR(TRIM(element_symbol))
          IF (my_charge_extended) THEN
             WRITE (UNIT=line(81:),FMT="(SP,F0.8)") qeff
          END IF
          WRITE (UNIT=iunit,FMT="(A)") TRIM(line)
       END DO
       WRITE (UNIT=iunit,FMT="(A)") "END"
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Illegal dump type")
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE write_particle_coordinates

! *****************************************************************************
!> \brief   Write the atomic coordinates to the output unit.
!> \author  MK
!> \date    05.06.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_fist_particle_coordinates(particle_set,subsys_section,&
             charges, error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'write_fist_particle_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: iatom, ikind, iw, natom
    REAL(KIND=dp)                            :: conv, mass, qeff
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%ATOMIC_COORDINATES",extension=".coordLog",error=error)

    CALL section_vals_val_get(subsys_section,"PRINT%ATOMIC_COORDINATES%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (iw>0) THEN
      WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
        "MODULE FIST:  ATOMIC COORDINATES IN "//TRIM(unit_str)
      WRITE (UNIT=iw,&
           FMT="(/,T3,A,7X,2(A1,11X),A1,8X,A8,5X,A6,/)")&
        "Atom  Kind  ATM_TYP","X","Y","Z","  q(eff)","  Mass"

      natom = SIZE(particle_set)
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             kind_number=ikind,&
                             name=name,&
                             mass=mass,&
                             qeff=qeff)
        IF (ASSOCIATED(charges)) qeff=charges(iatom)
        WRITE (UNIT=iw,&
               FMT="(T2,I5,1X,I4,3X,A4,3X,3F12.6,4X,F6.2,2X,F11.4)")&
          iatom,ikind,name,&
          particle_set(iatom)%r(1:3)*conv,qeff,mass/massunit
      END DO
      WRITE (iw,'(/)')
    END IF

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%ATOMIC_COORDINATES", error=error)

  END SUBROUTINE write_fist_particle_coordinates

! *****************************************************************************
!> \brief   Write the atomic coordinates to the output unit.
!> \author  MK
!> \date    05.06.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_qs_particle_coordinates(particle_set,subsys_section,label,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    CHARACTER(LEN=*), INTENT(IN)             :: label
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'write_qs_particle_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: unit_str
    INTEGER                                  :: handle, iatom, ikind, iw, &
                                                natom, z
    REAL(KIND=dp)                            :: conv, mass, zeff
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%ATOMIC_COORDINATES",extension=".coordLog",error=error)

    CALL section_vals_val_get(subsys_section,"PRINT%ATOMIC_COORDINATES%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (iw>0) THEN
      WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
        "MODULE "//TRIM(label)//":  ATOMIC COORDINATES IN "//TRIM(unit_str)
      WRITE (UNIT=iw,&
           FMT="(/,T3,A,7X,2(A1,11X),A1,8X,A8,5X,A6,/)")&
        "Atom  Kind  Element","X","Y","Z","  Z(eff)","  Mass"

      natom = SIZE(particle_set)
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             kind_number=ikind,&
                             element_symbol=element_symbol,&
                             mass=mass,&
                             z=z,&
                             zeff=zeff)
        WRITE (UNIT=iw,&
               FMT="(T2,I7,1X,I5,1X,A2,1X,I3,3F12.6,4X,F6.2,2X,F11.4)")&
          iatom,ikind,element_symbol,z,&
          particle_set(iatom)%r(1:3)*conv,zeff,mass/massunit
      END DO
      WRITE (iw,'(/)')
    END IF

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%ATOMIC_COORDINATES", error=error)

    CALL timestop(handle)

  END SUBROUTINE write_qs_particle_coordinates

! *****************************************************************************
!> \brief   Write the matrix of the particle distances to the output unit.
!> \author  Matthias Krack
!> \date    06.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_particle_distances(particle_set,cell,subsys_section,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_particle_distances', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: unit_str
    INTEGER                                  :: handle, iatom, iw, jatom, &
                                                natom, stat
    INTEGER, DIMENSION(3)                    :: periodic
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: conv, dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: distance_matrix
    REAL(KIND=dp), DIMENSION(3)              :: rab
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%INTERATOMIC_DISTANCES",extension=".distLog",error=error)

    CALL section_vals_val_get(subsys_section,"PRINT%INTERATOMIC_DISTANCES%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (iw>0) THEN
       CALL get_cell(cell=cell, periodic=periodic)
       natom = SIZE(particle_set)
       ALLOCATE (distance_matrix(natom,natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       distance_matrix(:,:) = 0.0_dp
       DO iatom=1,natom
          DO jatom=iatom+1,natom
             rab(:) = pbc(particle_set(iatom)%r(:),&
                  particle_set(jatom)%r(:),cell)
             dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
             distance_matrix(iatom,jatom) = dab*conv
             distance_matrix(jatom,iatom) = distance_matrix(iatom,jatom)
          END DO
       END DO

       ! Print the distance matrix
       WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
            "INTERATOMIC DISTANCES IN "//TRIM(unit_str)

       CALL write_particle_matrix(distance_matrix,particle_set,iw,error=error)
    END IF

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%INTERATOMIC_DISTANCES", error=error)

    CALL timestop(handle)

  END SUBROUTINE write_particle_distances

! *****************************************************************************
  SUBROUTINE write_particle_matrix(matrix,particle_set,iw,el_per_part,Ilist,error)
    REAL(KIND=dp), DIMENSION(:, :)           :: matrix
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: iw
    INTEGER, INTENT(IN), OPTIONAL            :: el_per_part
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Ilist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: from, i, iatom, icol, jatom, &
                                                katom, my_el_per_part, natom, &
                                                to
    INTEGER, DIMENSION(:), POINTER           :: my_list
    LOGICAL                                  :: failure

    failure = .FALSE.
    my_el_per_part = 1
    IF (PRESENT(el_per_part)) my_el_per_part = el_per_part
    IF (PRESENT(Ilist)) THEN
       natom = SIZE(Ilist)
    ELSE
       natom = SIZE(particle_set)
    END IF
    ALLOCATE(my_list(natom))
    IF (PRESENT(Ilist)) THEN
       my_list = Ilist
    ELSE
       DO i = 1, natom
          my_list(i) = i
       END DO
    END IF
    natom = natom*my_el_per_part
    DO jatom=1,natom,5
       from = jatom
       to = MIN(from+4,natom)
       WRITE (UNIT=iw,FMT="(/,T2,11X,5(4X,I5,4X))")&
            (icol,icol=from,to)
       DO iatom=1,natom
          katom = iatom/my_el_per_part
          IF (MOD(iatom,my_el_per_part)/=0) katom = katom+1
          CALL get_atomic_kind(atomic_kind=particle_set(my_list(katom))%atomic_kind,&
               element_symbol=element_symbol)
          WRITE (UNIT=iw,FMT="(T2,I5,2X,A2,2X,5(1X,F12.6))")&
               iatom,element_symbol,&
               (matrix(iatom,icol),icol=from,to)
       END DO
    END DO
    DEALLOCATE(my_list)
  END SUBROUTINE write_particle_matrix

! *****************************************************************************
!> \brief   Write structure data requested by a separate structure data input
!>          section to the output unit.
!>          input_section can be either motion_section or subsys_section.
!>
!> \author  MK
!> \date    11.03.04
!> \version 1.0
!> \par History
!>          Recovered (23.03.06,MK)
! *****************************************************************************
  SUBROUTINE write_structure_data(particle_set,cell,input_section,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: input_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_structure_data', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string, unit_str
    INTEGER                                  :: handle, i, i_rep, istat, iw, &
                                                n, n_rep, n_vals, natom, &
                                                new_size, old_size, wrk2(2), &
                                                wrk3(3), wrk4(4)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: work
    INTEGER, DIMENSION(:), POINTER           :: atomic_indices, index_list
    LOGICAL                                  :: failure, unique
    REAL(KIND=dp)                            :: conv, dab
    REAL(KIND=dp), DIMENSION(3)              :: r, rab, rbc, rcd, s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY (atomic_indices)
    NULLIFY (index_list)
    NULLIFY (logger)
    NULLIFY (section)
    string = ""

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger=logger,&
                              basis_section=input_section,&
                              print_key_path="PRINT%STRUCTURE_DATA",&
                              extension=".coordLog",&
                              error=error)

    CALL section_vals_val_get(input_section,"PRINT%STRUCTURE_DATA%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    CALL uppercase(unit_str)
    IF (iw > 0) THEN
       natom = SIZE(particle_set)
       section => section_vals_get_subs_vals(section_vals=input_section,&
                                             subsection_name="PRINT%STRUCTURE_DATA",&
                                             error=error)

      WRITE (UNIT=iw,FMT="(/,T2,A)") "REQUESTED STRUCTURE DATA"
      ! Print the requested atomic position vectors
      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="POSITION",&
                                n_rep_val=n_rep,&
                                error=error)
      IF (n_rep > 0) THEN
         WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Position vectors r(i) of the atoms i in "//TRIM(unit_str)
         old_size = 0
         DO i_rep=1,n_rep
            CALL section_vals_val_get(section_vals=section,&
                                      keyword_name="POSITION",&
                                      i_rep_val=i_rep,&
                                      i_vals=atomic_indices,&
                                      error=error)
            n_vals = SIZE(atomic_indices)
            new_size = old_size + n_vals
            CALL reallocate(index_list,1,new_size)
            index_list(old_size+1:new_size) = atomic_indices(1:n_vals)
            old_size = new_size
         END DO
         ALLOCATE (work(new_size),STAT=istat)
         CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
         CALL sort(index_list,new_size,work)
         DEALLOCATE (work,STAT=istat)
         CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
         DO i=1,new_size
            WRITE (UNIT=string,FMT="(A,I0,A)") "(",index_list(i),")"
            IF ((index_list(i) < 1).OR.(index_list(i) > natom)) THEN
               WRITE (UNIT=iw,FMT="(T3,A)")&
                "Invalid atomic index "//TRIM(string)//" specified. Print request is ignored."
               CYCLE
            END IF
            IF (i > 1) THEN
               ! Skip redundant indices
               IF (index_list(i) == index_list(i-1)) CYCLE
            END IF
            WRITE (UNIT=iw,FMT="(T3,A,T20,A,3F13.6)")&
             "r"//TRIM(string),"=",pbc(particle_set(index_list(i))%r(1:3),cell)*conv
         END DO
         DEALLOCATE (index_list,STAT=istat)
         CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
      END IF

      ! Print the requested atomic position vectors in scaled coordinates
      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="POSITION_SCALED",&
                                n_rep_val=n_rep,&
                                error=error)
      IF (n_rep > 0) THEN
         WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Position vectors s(i) of the atoms i in scaled coordinates"
         old_size = 0
         DO i_rep=1,n_rep
            CALL section_vals_val_get(section_vals=section,&
                                      keyword_name="POSITION_SCALED",&
                                      i_rep_val=i_rep,&
                                      i_vals=atomic_indices,&
                                      error=error)
            n_vals = SIZE(atomic_indices)
            new_size = old_size + n_vals
            CALL reallocate(index_list,1,new_size)
            index_list(old_size+1:new_size) = atomic_indices(1:n_vals)
            old_size = new_size
         END DO
         ALLOCATE (work(new_size),STAT=istat)
         CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
         CALL sort(index_list,new_size,work)
         DEALLOCATE (work,STAT=istat)
         CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
         DO i=1,new_size
            WRITE (UNIT=string,FMT="(A,I0,A)") "(",index_list(i),")"
            IF ((index_list(i) < 1).OR.(index_list(i) > natom)) THEN
               WRITE (UNIT=iw,FMT="(T3,A)")&
                "Invalid atomic index "//TRIM(string)//" specified. Print request is ignored."
               CYCLE
            END IF
            IF (i > 1) THEN
               ! Skip redundant indices
               IF (index_list(i) == index_list(i-1)) CYCLE
            END IF
            r(1:3) = pbc(particle_set(index_list(i))%r(1:3),cell)
            CALL real_to_scaled(s,r,cell)
            WRITE (UNIT=iw,FMT="(T3,A,T20,A,3F13.6)")&
             "s"//TRIM(string),"=",s(1:3)
         END DO
         DEALLOCATE (index_list,STAT=istat)
         CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
      END IF

      ! Print the requested distances
      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="DISTANCE",&
                                n_rep_val=n,&
                                error=error)
      IF (n > 0) THEN
        WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Distance vector r(i,j) between the atom i and j in "//&
          TRIM(unit_str)
        DO i=1,n
          CALL section_vals_val_get(section_vals=section,&
                                    keyword_name="DISTANCE",&
                                    i_rep_val=i,&
                                    i_vals=atomic_indices,&
                                    error=error)
          string = ""
          WRITE (UNIT=string,FMT="(A,2(I0,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),")"
          wrk2=atomic_indices
          CALL sort_unique(wrk2, unique)
          IF (((wrk2(1)>=1).AND.(wrk2(SIZE(wrk2))<=natom)).AND.unique) THEN
            rab(:) = pbc(particle_set(atomic_indices(1))%r(:),&
                         particle_set(atomic_indices(2))%r(:),cell)
            dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
            WRITE (UNIT=iw,FMT="(T3,A,T20,A,3F13.6,3X,A,F13.6)")&
              "r"//TRIM(string),"=",rab(:)*conv,&
              "|r| =",dab*conv
          ELSE
            WRITE (UNIT=iw,FMT="(T3,A)")&
              "Invalid atomic indices "//TRIM(string)//" specified. Print request is ignored."
          END IF
        END DO
      END IF

      ! Print the requested angles
      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="ANGLE",&
                                n_rep_val=n,&
                                error=error)
      IF (n > 0) THEN
        WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Angle a(i,j,k) between the atomic distance vectors r(j,i) and "//&
          "r(j,k) in DEGREE"
        DO i=1,n
          CALL section_vals_val_get(section_vals=section,&
                                    keyword_name="ANGLE",&
                                    i_rep_val=i,&
                                    i_vals=atomic_indices,&
                                    error=error)
          string = ""
          WRITE (UNIT=string,FMT="(A,3(I0,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),",",atomic_indices(3),")"
          wrk3=atomic_indices
          CALL sort_unique(wrk3, unique)
          IF (((wrk3(1)>=1).AND.(wrk3(SIZE(wrk3))<=natom)).AND.unique) THEN
            rab(:) = pbc(particle_set(atomic_indices(1))%r(:),&
                         particle_set(atomic_indices(2))%r(:),cell)
            rbc(:) = pbc(particle_set(atomic_indices(2))%r(:),&
                         particle_set(atomic_indices(3))%r(:),cell)
            WRITE (UNIT=iw,FMT="(T3,A,T26,A,F9.3)")&
              "a"//TRIM(string),"=",angle(-rab,rbc)*degree
          ELSE
            WRITE (UNIT=iw,FMT="(T3,A)")&
              "Invalid atomic indices "//TRIM(string)//" specified. Print request is ignored."
          END IF
        END DO
      END IF

      ! Print the requested dihedral angles
      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="DIHEDRAL_ANGLE",&
                                n_rep_val=n,&
                                error=error)
      IF (n > 0) THEN
        WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Dihedral angle d(i,j,k,l) between the planes (i,j,k) and (j,k,l) "//&
          "in DEGREE"
        DO i=1,n
          CALL section_vals_val_get(section_vals=section,&
                                    keyword_name="DIHEDRAL_ANGLE",&
                                    i_rep_val=i,&
                                    i_vals=atomic_indices,&
                                    error=error)
          string = ""
          WRITE (UNIT=string,FMT="(A,4(I0,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),",",&
                atomic_indices(3),",",atomic_indices(4),")"
          wrk4=atomic_indices
          CALL sort_unique(wrk4, unique)
          IF (((wrk4(1)>=1).AND.(wrk4(SIZE(wrk4))<=natom)).AND.unique) THEN
            rab(:) = pbc(particle_set(atomic_indices(1))%r(:),&
                         particle_set(atomic_indices(2))%r(:),cell)
            rbc(:) = pbc(particle_set(atomic_indices(2))%r(:),&
                         particle_set(atomic_indices(3))%r(:),cell)
            rcd(:) = pbc(particle_set(atomic_indices(3))%r(:),&
                         particle_set(atomic_indices(4))%r(:),cell)
            WRITE (UNIT=iw,FMT="(T3,A,T26,A,F9.3)")&
              "d"//TRIM(string),"=",dihedral_angle(rab,rbc,rcd)*degree
          ELSE
            WRITE (UNIT=iw,FMT="(T3,A)")&
              "Invalid atomic indices "//TRIM(string)//" specified. Print request is ignored."
          END IF
        END DO
      END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,input_section,&
         "PRINT%STRUCTURE_DATA", error=error)

    CALL timestop(handle)

  END SUBROUTINE write_structure_data

! *****************************************************************************
  SUBROUTINE update_particle_set ( particle_set, int_group,  pos, vel, for, add, error )

    TYPE(particle_type), POINTER             :: particle_set( : )
    INTEGER, INTENT(IN)                      :: int_group
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: pos( :, : ), vel( :, : ), &
                                                for(:,:)
    LOGICAL, INTENT(IN), OPTIONAL            :: add
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iparticle, nparticle
    LOGICAL                                  :: failure, my_add, update_for, &
                                                update_pos, update_vel

    CALL timeset(routineN,handle)

    failure = .FALSE.
    nparticle  = SIZE ( particle_set )
    update_pos = PRESENT(pos)
    update_vel = PRESENT(vel)
    update_for = PRESENT(for)
    my_add     = .FALSE.
    IF (PRESENT(add)) my_add = add

    IF (update_pos) THEN
       CALL mp_sum( pos,int_group)
       IF (my_add) THEN
          DO iparticle=1,nparticle
             particle_set(iparticle)%r(:) = particle_set(iparticle)%r(:) + pos(:,iparticle)
          END DO
       ELSE
          DO iparticle=1,nparticle
             particle_set(iparticle)%r(:) = pos(:,iparticle)
          END DO
       END IF
    END IF
    IF (update_vel) THEN
       CALL mp_sum( vel,int_group)
       IF (my_add) THEN
          DO iparticle=1,nparticle
             particle_set(iparticle)%v(:) = particle_set(iparticle)%v(:) + vel(:,iparticle)
          END DO
       ELSE
          DO iparticle=1,nparticle
             particle_set(iparticle)%v(:) = vel(:,iparticle)
          END DO
       END IF
    END IF
    IF (update_for) THEN
       CALL mp_sum( for,int_group)
       IF (my_add) THEN
          DO iparticle=1,nparticle
             particle_set(iparticle)%f(:) = particle_set(iparticle)%f(:) + for(:,iparticle)
          END DO
       ELSE
          DO iparticle=1,nparticle
             particle_set(iparticle)%f(:) = for(:,iparticle)
          END DO
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE update_particle_set

! *****************************************************************************
!> \brief   Return the atomic position or velocity of atom iatom in x from a
!>          packed vector even if core-shell particles are present
!> \author  Matthias Krack
!> \date    25.11.2010
!> \version 1.0
! *****************************************************************************
  FUNCTION get_particle_pos_or_vel(iatom,particle_set,vector) RESULT(x)

    INTEGER, INTENT(IN)                      :: iatom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector
    REAL(KIND=dp), DIMENSION(3)              :: x

    INTEGER                                  :: ic, is
    REAL(KIND=dp)                            :: fc, fs, mass

    ic = 3*(iatom - 1)
    IF (particle_set(iatom)%shell_index == 0) THEN
       x(1:3) = vector(ic+1:ic+3)
    ELSE
       is = 3*(SIZE(particle_set) + particle_set(iatom)%shell_index - 1)
       mass = particle_set(iatom)%atomic_kind%mass
       fc = particle_set(iatom)%atomic_kind%shell%mass_core/mass
       fs = particle_set(iatom)%atomic_kind%shell%mass_shell/mass
       x(1:3) = fc*vector(ic+1:ic+3) + fs*vector(is+1:is+3)
    END IF

  END FUNCTION get_particle_pos_or_vel

! *****************************************************************************
!> \brief   Update the atomic position or velocity by x and return the updated
!>          atomic position or velocity in x even if core-shell particles are
!>          present
!> \author  Matthias Krack
!> \date    26.11.2010
!> \version 1.0
!> \note    particle-set is not changed, only the positions or velocities in
!>          the packed vector are updated
! *****************************************************************************
  SUBROUTINE update_particle_pos_or_vel(iatom,particle_set,x,vector)

    INTEGER, INTENT(IN)                      :: iatom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: x
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: vector

    INTEGER                                  :: ic, is
    REAL(KIND=dp)                            :: fc, fs, mass

    ic = 3*(iatom - 1)
    IF (particle_set(iatom)%shell_index == 0) THEN
       vector(ic+1:ic+3) = vector(ic+1:ic+3) + x(1:3)
       x(1:3) = vector(ic+1:ic+3)
    ELSE
       is = 3*(SIZE(particle_set) + particle_set(iatom)%shell_index - 1)
       mass = particle_set(iatom)%atomic_kind%mass
       fc = particle_set(iatom)%atomic_kind%shell%mass_core/mass
       fs = particle_set(iatom)%atomic_kind%shell%mass_shell/mass
       vector(ic+1:ic+3) = vector(ic+1:ic+3) + x(1:3)
       vector(is+1:is+3) = vector(is+1:is+3) + x(1:3)
       x(1:3) = fc*vector(ic+1:ic+3) + fs*vector(is+1:is+3)
    END IF

  END SUBROUTINE update_particle_pos_or_vel

END MODULE particle_types
