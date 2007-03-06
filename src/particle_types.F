!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/particle_types [1.0] *
!!
!!   NAME
!!     particle_types
!!
!!   FUNCTION
!!     Define the data structure for the particle information.
!!
!!   AUTHOR
!!     CJM, MK
!!
!!   MODIFICATION HISTORY
!!     - Atomic kind added in particle_type (MK,08.01.2002)
!!     - Functionality for particle_type added (MK,14.01.2002)
!!     - Allow for general coordinate input (MK,13.09.2003)
!!     - Molecule concept introduced (MK,26.09.2003)
!!     - Last atom information added (jgh,23.05.2004)
!!     - particle_type cleaned (MK,03.02.2005)
!!
!!   SOURCE
!******************************************************************************

MODULE particle_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_xmol
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size,&
                                             sp
  USE mathconstants,                   ONLY: degree
  USE mathlib,                         ONLY: angle,&
                                             dihedral_angle
  USE message_passing,                 ONLY: mp_sum
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             massunit
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'particle_types'

  ! Data types
  TYPE particle_type
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_vdw
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_ei
     TYPE(atomic_kind_type), POINTER       :: atomic_kind ! atomic kind information
     REAL(KIND = dp), DIMENSION(3)         :: f,&         ! force
                                              r,&         ! position
                                              v           ! velocity
     ! particle dependent terms for EAM
     REAL(KIND = dp)                       :: rho, f_embed
     ! particle dependent terms for shell-model
     INTEGER                               :: atom_index, shell_index
  END TYPE particle_type

  ! Public data types

  PUBLIC :: particle_type

  ! Public subroutines

  PUBLIC :: allocate_particle_set,&
            clone_particle,&
            deallocate_particle_set,&
            get_particle_set,&
            write_fist_particle_coordinates,&
            write_qs_particle_coordinates,&
            write_particle_distances,&
            write_particle_coordinates,&
            write_structure_data,&
            write_particle_matrix,&
            update_particle_set

!!***
! *****************************************************************************

CONTAINS

!!****f* particle_types/clone_particle [1.0] *
!!
!!   NAME
!!     clone_particle
!!
!!   FUNCTION
!!     copies all the attributes of a particle to another
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_in: the particle to copy
!!     - p_out: the place where to copy
!!     - nullify_non_owned_pointers: if the pointers to shared structures
!!       should be nullified. If false, copies the pointer (but does not
!!       try to perform any memory handling).
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE clone_particle(p_in,p_out,nullify_non_owned_pointers,error)
    TYPE(particle_type), INTENT(in)          :: p_in
    TYPE(particle_type), INTENT(out)         :: p_out
    LOGICAL, INTENT(in)                      :: nullify_non_owned_pointers
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'clone_particle', &
      routineP = moduleN//':'//routineN

  NULLIFY(p_out%atomic_kind)
  NULLIFY(p_out%list_exclude_vdw)
  NULLIFY(p_out%list_exclude_ei)
  IF (.NOT.nullify_non_owned_pointers) THEN
     p_out%atomic_kind => p_in%atomic_kind
     p_out%list_exclude_vdw => p_in%list_exclude_vdw
     p_out%list_exclude_ei => p_in%list_exclude_ei
  END IF
  p_out%r=p_in%r
  p_out%f=p_in%f
  p_out%v=p_in%v
  p_out%shell_index=p_in%shell_index
  p_out%atom_index=p_in%atom_index

END SUBROUTINE clone_particle
!***************************************************************************

! *****************************************************************************

  SUBROUTINE allocate_particle_set(particle_set,nparticle)

!   Purpose: Allocate a particle set.

!   History: - Creation (14.01.2002,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: nparticle

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_particle_set'

    INTEGER                                  :: iparticle, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) CALL deallocate_particle_set(particle_set)

    ALLOCATE (particle_set(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "particle_set",nparticle*int_size)

    DO iparticle=1,nparticle
      NULLIFY (particle_set(iparticle)%atomic_kind)
      NULLIFY (particle_set(iparticle)%list_exclude_vdw)
      NULLIFY (particle_set(iparticle)%list_exclude_ei)
      particle_set(iparticle)%f(:) = 0.0_dp
      particle_set(iparticle)%r(:) = 0.0_dp
      particle_set(iparticle)%v(:) = 0.0_dp
      particle_set(iparticle)%shell_index = 0
      particle_set(iparticle)%atom_index = 0
    END DO

  END SUBROUTINE allocate_particle_set

! *****************************************************************************

  SUBROUTINE deallocate_particle_set(particle_set)

!   Purpose: Deallocate a particle set.

!   History: - Creation (14.01.2002,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_particle_set'

    INTEGER                                  :: iparticle, istat

!   ---------------------------------------------------------------------------


    IF (ASSOCIATED(particle_set)) THEN

      DO iparticle=1,SIZE(particle_set)
         IF (ASSOCIATED(particle_set(iparticle)%list_exclude_vdw,&
                        particle_set(iparticle)%list_exclude_ei)) THEN
            DEALLOCATE(particle_set(iparticle)%list_exclude_vdw,STAT=istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                 "particle_set(iparticle)%list_exclude_vdw")
            NULLIFY(particle_set(iparticle)%list_exclude_vdw,&
                    particle_set(iparticle)%list_exclude_ei)
         ELSE
            IF (ASSOCIATED(particle_set(iparticle)%list_exclude_vdw)) THEN
               DEALLOCATE(particle_set(iparticle)%list_exclude_vdw,STAT=istat)
               IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                    "particle_set(iparticle)%list_exclude_vdw")
               NULLIFY(particle_set(iparticle)%list_exclude_vdw)
            END IF
            IF (ASSOCIATED(particle_set(iparticle)%list_exclude_ei)) THEN
               DEALLOCATE(particle_set(iparticle)%list_exclude_ei,STAT=istat)
               IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                    "particle_set(iparticle)%list_exclude_ei")
               NULLIFY(particle_set(iparticle)%list_exclude_ei)
            END IF
         END IF
         NULLIFY(particle_set(iparticle)%atomic_kind)
      END DO
      DEALLOCATE (particle_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "particle_set")
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer particle_set is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_particle_set

! *****************************************************************************

  SUBROUTINE get_particle_set(particle_set,first_cgf,first_nco,first_nso,&
                              first_sgf,last_cgf,last_nco,last_nso,&
                              last_sgf,ncgf,ncotot,nsotot,nsgf,r)

!   Purpose: Get the components of a particle set.

!   History: - Creation (14.01.2002,MK)
!            - particle type cleaned (13.10.2003,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: first_cgf, first_nco, &
                                                first_nso, first_sgf, &
                                                last_cgf, last_nco, last_nso, &
                                                last_sgf, ncgf
    INTEGER, INTENT(OUT), OPTIONAL           :: ncotot, nsotot
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: nsgf
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: r

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'get_particle_set'

    INTEGER                                  :: icgf, ico, iparticle, isgf, &
                                                iso, maxco, maxso, nc, &
                                                nparticle, ns, nset
    TYPE(gto_basis_set_type), POINTER        :: orb_basis

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) THEN

      nparticle = SIZE(particle_set)

      IF (PRESENT(first_cgf)) THEN
        IF (SIZE(first_cgf) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <first_cgf>")
        END IF
      END IF

      IF (PRESENT(first_nco)) THEN
        IF (SIZE(first_nco) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <first_nco>")
        END IF
      END IF

      IF (PRESENT(first_nso)) THEN
        IF (SIZE(first_nso) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <first_nso>")
        END IF
      END IF


      IF (PRESENT(first_sgf)) THEN
        IF (SIZE(first_sgf) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <first_sgf>")
        END IF
      END IF

      IF (PRESENT(last_cgf)) THEN
        IF (SIZE(last_cgf) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <last_cgf>")
        END IF
      END IF

      IF (PRESENT(last_nco)) THEN
        IF (SIZE(last_nco) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <last_nco>")
        END IF
      END IF

      IF (PRESENT(last_sgf)) THEN
        IF (SIZE(last_sgf) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <last_sgf>")
        END IF
      END IF

      IF (PRESENT(ncgf)) THEN
        IF (SIZE(ncgf) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <ncgf>")
        END IF
      END IF

      IF (PRESENT(nsgf)) THEN
        IF (SIZE(nsgf) < nparticle) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <nsgf>")
        END IF
      END IF

      IF (PRESENT(first_cgf).OR.&
          PRESENT(last_cgf).OR.&
          PRESENT(first_sgf).OR.&
          PRESENT(last_sgf).OR.&
          PRESENT(ncgf).OR.&
          PRESENT(nsgf)) THEN
        icgf = 0
        isgf = 0
        DO iparticle=1,nparticle
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                               ncgf=nc,&
                               nsgf=ns)
          IF (PRESENT(ncgf)) ncgf(iparticle) = nc
          IF (PRESENT(first_cgf)) first_cgf(iparticle) = icgf + 1
          icgf = icgf + nc
          IF (PRESENT(last_cgf)) last_cgf(iparticle) = icgf
          IF (PRESENT(nsgf)) nsgf(iparticle) = ns
          IF (PRESENT(first_sgf)) first_sgf(iparticle) = isgf + 1
          isgf = isgf + ns
          IF (PRESENT(last_sgf)) last_sgf(iparticle) = isgf
        END DO
        IF (PRESENT(first_cgf)) THEN
          IF (SIZE(first_cgf) > nparticle) first_cgf(nparticle+1) = icgf + 1
        END IF
        IF (PRESENT(first_sgf)) THEN
          IF (SIZE(first_sgf) > nparticle) first_sgf(nparticle+1) = isgf + 1
        END IF
      END IF

      IF (PRESENT(first_nco).OR.&
          PRESENT(last_nco).OR.&
          PRESENT(ncotot)) THEN
        ico = 0
        DO iparticle = 1,nparticle
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                               orb_basis_set=orb_basis)
          CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nset,maxco=maxco)
          nc = maxco*nset
          IF(PRESENT(first_nco)) first_nco(iparticle) = ico + 1
          IF(PRESENT(last_nco))  last_nco(iparticle)  = ico + nc
          ico = ico + nc
        ENDDO
        IF(PRESENT(ncotot)) ncotot = ico

      ENDIF

      IF (PRESENT(first_nso).OR.&
          PRESENT(last_nso).OR.&
          PRESENT(nsotot)) THEN
        iso = 0
        DO iparticle = 1,nparticle
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                               orb_basis_set=orb_basis)
          CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nset,maxso=maxso)
          ns = maxso*nset
          IF(PRESENT(first_nso)) first_nso(iparticle) = iso + 1
          IF(PRESENT(last_nso))  last_nso(iparticle)  = iso + ns
          iso = iso + ns
        ENDDO
        IF(PRESENT(nsotot)) nsotot = iso
      ENDIF

      IF (PRESENT(r)) THEN
        DO iparticle=1,nparticle
          r(:,iparticle) = particle_set(iparticle)%r(:)
        END DO
      END IF

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer particle_set is not associated")

    END IF

  END SUBROUTINE get_particle_set

! *****************************************************************************
!  should be able to write a few formats e.g. xmol, and some binary format (dcd)
!  some format can be used for x,v,f
!  should work for fist, qs, kg
!  still needs some rethinking of the units
!  FORMAT   CONTENT                                    UNITS x,v,f
!  XMOL     POS, VEL, FORCE, POS_VEL, POS_VEL_FORCE    Angstrom, a.u., a.u. !hmmmm
! *****************************************************************************

  SUBROUTINE write_particle_coordinates(particle_set,iunit,FORMAT,&
       content,title,cell,array)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: iunit, FORMAT
    CHARACTER(LEN=*)                         :: content, title
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: array

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=4)                         :: name
    CHARACTER(LEN=default_string_length)     :: atm_name, my_format
    INTEGER                                  :: i, iatom, istat, j, natom, &
                                                stat
    LOGICAL                                  :: dummy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: arr
    REAL(KIND=dp), DIMENSION(3)              :: f, r, v
    REAL(KIND=sp), DIMENSION(:), POINTER     :: x4, y4, z4

    natom = SIZE(particle_set)
    IF (PRESENT(array)) THEN
       SELECT CASE(TRIM(content))
       CASE ("POS_VEL","POS_VEL_FORCE")
          CALL stop_program ("particle_types","illegal usage of write_particle_coordinates")
       END SELECT
    END IF
    SELECT CASE (FORMAT)
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
                r(:) = array((iatom-1)*3+1:(iatom-1)*3+3)*angstrom
             ELSE
                r(:) = particle_set(iatom)%r(:)*angstrom
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),r(1:3)
          CASE ("VEL")
             IF (PRESENT(array)) THEN
                v(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                v(:) = particle_set(iatom)%v(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),v(1:3)
          CASE ("FORCE")
             IF (PRESENT(array)) THEN
                f(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                f(:) = particle_set(iatom)%f(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),f(1:3)
          CASE ("POS_VEL")
             r(:) = particle_set(iatom)%r(:)*angstrom
             WRITE (iunit,TRIM(my_format)//"1X,6F20.10)") TRIM(name),r(1:3),&
                  particle_set(iatom)%v(1:3)
          CASE ("POS_VEL_FORCE")
             r(:) = particle_set(iatom)%r(:)*angstrom
             WRITE (iunit,TRIM(my_format)//"1X,9F20.10)") TRIM(name),r(1:3),&
                  particle_set(iatom)%v(1:3),&
                  particle_set(iatom)%f(1:3)
          END SELECT
       END DO
    CASE (dump_atomic)
       DO iatom=1,natom
          SELECT CASE (TRIM(content))
          CASE ("POS")
             IF (PRESENT(array)) THEN
                r(:) = array((iatom-1)*3+1:(iatom-1)*3+3)*angstrom
             ELSE
                r(:) = particle_set(iatom)%r(:)*angstrom
             END IF
             WRITE (iunit,"(3F20.10)") r(1:3)
          CASE ("VEL")
             IF (PRESENT(array)) THEN
                v(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                v(:) = particle_set(iatom)%v(:)
             END IF
             WRITE (iunit,"(3F20.10)")v(1:3)
          CASE ("FORCE")
             IF (PRESENT(array)) THEN
                f(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                f(:) = particle_set(iatom)%f(:)
             END IF
             WRITE (iunit,"(3F20.10)") f(1:3)
          CASE ("POS_VEL")
             WRITE (iunit,"(6F20.10)")particle_set(iatom)%r(1:3),&
                  particle_set(iatom)%v(1:3)
          CASE ("POS_VEL_FORCE")
             r(:) = particle_set(iatom)%r(:)
             WRITE (iunit,"(9F20.10)")r(1:3),&
                  particle_set(iatom)%v(1:3),&
                  particle_set(iatom)%f(1:3)
          END SELECT
       END DO
    CASE (dump_dcd)
       IF (.NOT.(PRESENT(cell))) CALL stop_memory('particle_types','cell not present')
       IF (PRESENT(array)) THEN
          ALLOCATE(arr(3,natom),stat=stat)
          IF (stat /=0) CALL stop_program ("particle_types","arr")
          arr = RESHAPE(array,(/3,natom/))
       END IF
       NULLIFY(x4,y4,z4)
       ALLOCATE(x4(natom),y4(natom),z4(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory('particle_types','x4,y4,z4')
       WRITE(iunit) ((cell%hmat(i,j)*angstrom,i=1,j),j=1,3)
       SELECT CASE (TRIM(content))
       CASE ("POS")
          IF (PRESENT(array)) THEN
             x4(:) = arr(1,:)*angstrom
             y4(:) = arr(2,:)*angstrom
             z4(:) = arr(3,:)*angstrom
          ELSE
             x4(:) = particle_set(:)%r(1)*angstrom
             y4(:) = particle_set(:)%r(2)*angstrom
             z4(:) = particle_set(:)%r(3)*angstrom
          END IF
          WRITE(iunit) x4
          WRITE(iunit) y4
          WRITE(iunit) z4
       CASE ("VEL")
          IF (PRESENT(array)) THEN
             x4(:) = arr(1,:)*angstrom*femtoseconds
             y4(:) = arr(2,:)*angstrom*femtoseconds
             z4(:) = arr(3,:)*angstrom*femtoseconds
          ELSE
             x4(:) = particle_set(:)%v(1)*angstrom*femtoseconds
             y4(:) = particle_set(:)%v(2)*angstrom*femtoseconds
             z4(:) = particle_set(:)%v(3)*angstrom*femtoseconds
          END IF
          WRITE(iunit) x4
          WRITE(iunit) y4
          WRITE(iunit) z4
       CASE ("FORCE")
          IF (PRESENT(array)) THEN
             x4(:) = arr(1,:)
             y4(:) = arr(2,:)
             z4(:) = arr(3,:)
          ELSE
             x4(:) = particle_set(:)%f(1)
             y4(:) = particle_set(:)%f(2)
             z4(:) = particle_set(:)%f(3)
          END IF
          WRITE(iunit) x4
          WRITE(iunit) y4
          WRITE(iunit) z4
       CASE DEFAULT
          CALL stop_program ("particle_types","illegal dcd dump type")
       END SELECT
       IF (PRESENT(array)) THEN
          DEALLOCATE(arr,STAT=istat)
          IF (istat /= 0) CALL stop_memory('particle_types','arr')
       END IF
       DEALLOCATE(x4,y4,z4,STAT=istat)
       IF (istat /= 0) CALL stop_memory('particle_types','x4,y4,z4')
    CASE DEFAULT
       CALL stop_program ("particle_types","illegal dump type")
    END SELECT

  END SUBROUTINE write_particle_coordinates

! *****************************************************************************

  SUBROUTINE write_fist_particle_coordinates(particle_set,cell,force_env_section,error)

!   Purpose: Write the atomic coordinates to the output unit.

!   History: - Creation (05.06.2000,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: iatom, ikind, iw, natom
    REAL(KIND=dp)                            :: mass, qeff, unit_of_length
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%ATOMIC_COORDINATES",extension=".coordLog",error=error)

    IF (iw>0) THEN

      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)

      WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
        "ATOMIC COORDINATES IN "//TRIM(unit_of_length_name)
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
        WRITE (UNIT=iw,&
               FMT="(T2,I5,1X,I4,3X,A4,3X,3F12.6,4X,F6.2,2X,F11.4)")&
          iatom,ikind,name,&
          particle_set(iatom)%r(1:3)/unit_of_length,qeff,mass/massunit
      END DO
      WRITE (iw,'(/)')
    END IF

    CALL cp_print_key_finished_output(iw,logger,force_env_section,&
         "SUBSYS%PRINT%ATOMIC_COORDINATES", error=error)

  END SUBROUTINE write_fist_particle_coordinates

  SUBROUTINE write_qs_particle_coordinates(particle_set,cell,force_env_section,error)

!   Purpose: Write the atomic coordinates to the output unit.

!   History: - Creation (05.06.2000,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER                                  :: iatom, ikind, iw, natom, z
    REAL(KIND=dp)                            :: mass, unit_of_length, zeff
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%ATOMIC_COORDINATES",extension=".coordLog",error=error)

    IF (iw>0) THEN

      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)

      WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
        "ATOMIC COORDINATES IN "//TRIM(unit_of_length_name)
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
               FMT="(T2,I5,1X,I4,3X,A2,2X,I3,3F12.6,4X,F6.2,2X,F11.4)")&
          iatom,ikind,element_symbol,z,&
          particle_set(iatom)%r(1:3)/unit_of_length,zeff,mass/massunit
      END DO
      WRITE (iw,'(/)')
    END IF

    CALL cp_print_key_finished_output(iw,logger,force_env_section,&
         "SUBSYS%PRINT%ATOMIC_COORDINATES", error=error)

  END SUBROUTINE write_qs_particle_coordinates


  SUBROUTINE write_particle_distances(particle_set,cell,force_env_section,error)

!   Purpose: Write the matrix of the particle distances to the output unit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_particle_distances'

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER                                  :: iatom, istat, iw, jatom, natom
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(KIND=dp)                            :: dab, unit_of_length
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: distance_matrix
    REAL(KIND=dp), DIMENSION(3)              :: rab
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%INTERATOMIC_DISTANCES",extension=".distLog",error=error)

    IF (iw>0) THEN

      CALL get_cell(cell=cell,&
                    periodic=periodic,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)


      natom = SIZE(particle_set)

      ALLOCATE (distance_matrix(natom,natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "distance_matrix",dp_size*natom*natom)

      distance_matrix(:,:) = 0.0_dp

      DO iatom=1,natom
        DO jatom=iatom+1,natom
          rab(:) = pbc(particle_set(iatom)%r(:),&
                       particle_set(jatom)%r(:),cell)
          dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
          distance_matrix(iatom,jatom) = dab/unit_of_length
          distance_matrix(jatom,iatom) = distance_matrix(iatom,jatom)
        END DO
      END DO

!     *** Print the distance matrix ***

      WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
        "INTERATOMIC DISTANCES IN "//TRIM(unit_of_length_name)

      CALL write_particle_matrix(distance_matrix,particle_set,iw,error=error)

    END IF

    CALL cp_print_key_finished_output(iw,logger,force_env_section,&
         "SUBSYS%PRINT%INTERATOMIC_DISTANCES", error=error)
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

  SUBROUTINE write_structure_data(particle_set,cell,force_env_section,error)

    ! Write structure data requested by a separate structure data input
    ! section to the output unit.

    ! Creation (11.03.04,MK)
    ! Recovered (23.03.06,MK)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_structure_data', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: handle, i, iw, n, natom
    INTEGER, DIMENSION(:), POINTER           :: atomic_indices
    REAL(KIND=dp)                            :: dab, unit_of_length
    REAL(KIND=dp), DIMENSION(3)              :: rab, rbc, rcd
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: section

! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY (atomic_indices)
    NULLIFY (logger)
    NULLIFY (section)

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger=logger,&
                              basis_section=force_env_section,&
                              print_key_path="SUBSYS%PRINT%STRUCTURE_DATA",&
                              extension=".coordLog",&
                              error=error)

    IF (iw > 0) THEN

      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)

      natom = SIZE(particle_set)

      section => section_vals_get_subs_vals(section_vals=force_env_section,&
                                            subsection_name="SUBSYS%PRINT%STRUCTURE_DATA",&
                                            error=error)

      WRITE (UNIT=iw,FMT="(/,T2,A)") "REQUESTED STRUCTURE DATA"

      ! Print the requested atomic position vectors

      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="POSITION",&
                                n_rep_val=n,&
                                error=error)
      IF (n > 0) THEN
        WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Position vector r(i) of atom i in "//TRIM(unit_of_length_name)
        DO i=1,n
          CALL section_vals_val_get(section_vals=section,&
                                    keyword_name="POSITION",&
                                    i_rep_val=i,&
                                    i_vals=atomic_indices,&
                                    error=error)
          string = ""
          WRITE (UNIT=string,FMT="(A,I6,A)")&
            "(",atomic_indices(1),")"
          CALL compress(string,full=.TRUE.)
          IF (valid_indices(atomic_indices,1,1,natom)) THEN
            WRITE (UNIT=iw,FMT="(T3,A,T20,A,3F13.6)")&
              "r"//TRIM(string),"=",&
              pbc(particle_set(atomic_indices(1))%r(1:3),cell)/unit_of_length
          ELSE
            WRITE (UNIT=iw,FMT="(T3,A)")&
              "Invalid atomic index "//TRIM(string)//" specified. Print request is ignored."
          END IF
        END DO
      END IF

      ! Print the requested distances

      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="DISTANCE",&
                                n_rep_val=n,&
                                error=error)
      IF (n > 0) THEN
        WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Distance vector r(i,j) between the atom i and j in "//&
          TRIM(unit_of_length_name)
        DO i=1,n
          CALL section_vals_val_get(section_vals=section,&
                                    keyword_name="DISTANCE",&
                                    i_rep_val=i,&
                                    i_vals=atomic_indices,&
                                    error=error)
          string = ""
          WRITE (UNIT=string,FMT="(A,2(I6,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),")"
          CALL compress(string,full=.TRUE.)
          IF (valid_indices(atomic_indices,2,1,natom)) THEN
            rab(:) = pbc(particle_set(atomic_indices(1))%r(:),&
                         particle_set(atomic_indices(2))%r(:),cell)
            dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
            WRITE (UNIT=iw,FMT="(T3,A,T20,A,3F13.6,3X,A,F13.6)")&
              "r"//TRIM(string),"=",rab(:)/unit_of_length,&
              "|r| =",dab/unit_of_length
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
          WRITE (UNIT=string,FMT="(A,3(I6,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),",",atomic_indices(3),")"
          CALL compress(string,full=.TRUE.)
          IF (valid_indices(atomic_indices,3,1,natom)) THEN
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
          WRITE (UNIT=string,FMT="(A,4(I6,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),",",&
                atomic_indices(3),",",atomic_indices(4),")"
          CALL compress(string,full=.TRUE.)
          IF (valid_indices(atomic_indices,4,1,natom)) THEN
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

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE write_structure_data

  ! ***************************************************************************

  FUNCTION valid_indices(indices,n,min_index,max_index) RESULT(indices_ok)

    ! Check for a set of valid indices.
    ! Creation (23.03.06,MK)

    INTEGER, DIMENSION(:), POINTER           :: indices
    INTEGER, INTENT(IN)                      :: n, min_index, max_index
    LOGICAL                                  :: indices_ok

    INTEGER                                  :: i, j

! -------------------------------------------------------------------------

    indices_ok = .TRUE.

    DO i=1,n
      IF ((indices(i) < min_index).OR.&
          (indices(i) > max_index)) THEN
        indices_ok = .FALSE.
        RETURN
      END IF
      DO j=i+1,n
        IF (indices(j) == indices(i)) THEN
          indices_ok = .FALSE.
          RETURN
        END IF
      END DO
    END DO

  END FUNCTION valid_indices

  ! ***************************************************************************

  !******************************************************************************
  SUBROUTINE update_particle_set ( particle_set, int_group,  pos, vel, for, add )
    
    TYPE(particle_type), POINTER             :: particle_set( : )
    INTEGER, INTENT(IN)                      :: int_group
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: pos( :, : ), vel( :, : ), &
                                                for(:,:)
    LOGICAL, INTENT(IN), OPTIONAL            :: add

    INTEGER                                  :: iparticle, nparticle
    LOGICAL                                  :: my_add, update_for, &
                                                update_pos, update_vel

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
    ENDIF
    
  END SUBROUTINE update_particle_set

END MODULE particle_types
