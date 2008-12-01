!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_xmol,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             sp
  USE mathconstants,                   ONLY: degree
  USE mathlib,                         ONLY: angle,&
                                             dihedral_angle
  USE message_passing,                 ONLY: mp_sum
  USE physcon,                         ONLY: massunit
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort_unique
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'particle_types'

  ! Data types
! *****************************************************************************
  TYPE particle_type
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_vdw
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_ei
     TYPE(atomic_kind_type), POINTER       :: atomic_kind ! atomic kind information
     REAL(KIND = dp), DIMENSION(3)         :: f,&         ! force
                                              r,&         ! position
                                              v           ! velocity
     ! Particle dependent terms for EAM
     REAL(KIND=dp)                         :: rho, f_embed
     ! Particle dependent terms for shell-model
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

CONTAINS

! *****************************************************************************
!> \brief copies all the attributes of a particle to another
!> \param p_in the particle to copy
!> \param p_out the place where to copy
!> \param nullify_non_owned_pointers if the pointers to shared structures
!>        should be nullified. If false, copies the pointer (but does not
!>        try to perform any memory handling).
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
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

    INTEGER                                  :: iparticle, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(particle_set)) THEN
       DO iparticle=1,SIZE(particle_set)
          IF (ASSOCIATED(particle_set(iparticle)%list_exclude_vdw,&
               particle_set(iparticle)%list_exclude_ei)) THEN
             DEALLOCATE(particle_set(iparticle)%list_exclude_vdw,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             NULLIFY(particle_set(iparticle)%list_exclude_ei)
          ELSE
             IF (ASSOCIATED(particle_set(iparticle)%list_exclude_vdw)) THEN
                DEALLOCATE(particle_set(iparticle)%list_exclude_vdw,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(particle_set(iparticle)%list_exclude_ei)) THEN
                DEALLOCATE(particle_set(iparticle)%list_exclude_ei,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
         END IF
         NULLIFY(particle_set(iparticle)%atomic_kind)
      END DO
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
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_particle_set(particle_set,first_cgf,first_nco,first_nso,&
                              first_sgf,last_cgf,last_nco,last_nso,&
                              last_sgf,ncgf,ncotot,nsotot,nsgf,r,basis_set_id,first_pgf,last_pgf,npgf,error)    !JT

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
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: first_pgf, last_pgf, npgf
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icgf, ico, iparticle, ipgf, &
                                                isgf, iso, maxco, maxso, &
                                                my_basis_set_id, nc, np, &
                                                nparticle, ns, nset
    TYPE(gto_basis_set_type), POINTER        :: orb_basis

    IF( PRESENT(basis_set_id) ) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

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
          PRESENT(first_pgf).OR.&  !JT
          PRESENT(last_pgf).OR.&   !JT
          PRESENT(ncgf).OR.&
          PRESENT(nsgf).OR.&
          PRESENT(npgf)) THEN          !JT
        ipgf = 0                       !JT
        icgf = 0
        isgf = 0
        DO iparticle=1,nparticle
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                               ncgf=nc,&
                               nsgf=ns,&
                               npgf=np,&                      !JT
                               basis_set_id=my_basis_set_id)
          IF (PRESENT(ncgf)) ncgf(iparticle) = nc
          IF (PRESENT(first_cgf)) first_cgf(iparticle) = icgf + 1
          icgf = icgf + nc
          IF (PRESENT(last_cgf)) last_cgf(iparticle) = icgf
          IF (PRESENT(nsgf)) nsgf(iparticle) = ns
          IF (PRESENT(first_sgf)) first_sgf(iparticle) = isgf + 1
          isgf = isgf + ns
          IF (PRESENT(last_sgf)) last_sgf(iparticle) = isgf
!JTs
          IF (PRESENT(npgf)) npgf(iparticle) = np
          IF (PRESENT(first_pgf)) first_pgf(iparticle) = ipgf + 1
          ipgf = ipgf + np
          IF (PRESENT(last_pgf)) last_pgf(iparticle) = ipgf
!JTe
        END DO
        IF (PRESENT(first_cgf)) THEN
          IF (SIZE(first_cgf) > nparticle) first_cgf(nparticle+1) = icgf + 1
        END IF
        IF (PRESENT(first_sgf)) THEN
          IF (SIZE(first_sgf) > nparticle) first_sgf(nparticle+1) = isgf + 1
        END IF
!JTs
        IF (PRESENT(first_pgf)) THEN
          IF (SIZE(first_pgf) > nparticle) first_pgf(nparticle+1) = ipgf + 1
        END IF
!JTe
      END IF

      IF (PRESENT(first_nco).OR.&
          PRESENT(last_nco).OR.&
          PRESENT(ncotot)) THEN
        ico = 0
        DO iparticle = 1,nparticle
          SELECT CASE (my_basis_set_id)
          CASE (use_orb_basis_set)
            CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                                 orb_basis_set=orb_basis)
            CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nset,maxco=maxco)
          CASE (use_aux_fit_basis_set)
            CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                                 aux_fit_basis_set=orb_basis)
            CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nset,maxco=maxco)
          END SELECT
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
          SELECT CASE (my_basis_set_id)
          CASE (use_orb_basis_set)
            CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                                 orb_basis_set=orb_basis)
            CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nset,maxso=maxso)
          CASE (use_aux_fit_basis_set)
            CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                                 aux_fit_basis_set=orb_basis)
            CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nset,maxso=maxso)
          END SELECT
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
       content,title,cell,array,unit_conv,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: iunit, output_format
    CHARACTER(LEN=*)                         :: content, title
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: array
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: unit_conv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_particle_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=4)                         :: name
    CHARACTER(LEN=default_string_length)     :: atm_name, my_format
    INTEGER                                  :: iatom, natom, stat
    LOGICAL                                  :: dummy, failure
    REAL(KIND=dp)                            :: alpha, beta, factor, gamma
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: arr
    REAL(KIND=dp), DIMENSION(3)              :: abc, f, r, v
    REAL(KIND=sp), DIMENSION(:), POINTER     :: x4, y4, z4

    failure = .FALSE.
    natom = SIZE(particle_set)
    IF (PRESENT(array)) THEN
       SELECT CASE(TRIM(content))
       CASE ("POS_VEL","POS_VEL_FORCE")
          CALL stop_program ("particle_types","illegal usage of write_particle_coordinates")
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
                r(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                r(:) = particle_set(iatom)%r(:)
             END IF
             WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),r(1:3)*factor
          CASE ("VEL")
             IF (PRESENT(array)) THEN
                v(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
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
          END SELECT
       END DO
    CASE (dump_atomic)
       DO iatom=1,natom
          SELECT CASE (TRIM(content))
          CASE ("POS")
             IF (PRESENT(array)) THEN
                r(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
             ELSE
                r(:) = particle_set(iatom)%r(:)
             END IF
             WRITE (iunit,"(3F20.10)")r(1:3)*factor
          CASE ("VEL")
             IF (PRESENT(array)) THEN
                v(:) = array((iatom-1)*3+1:(iatom-1)*3+3)
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
          END SELECT
       END DO
    CASE (dump_dcd)
       IF (.NOT.(PRESENT(cell))) THEN 
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Cell is not present! Report this bug!")
       END IF
       IF (PRESENT(array)) THEN
          ALLOCATE(arr(3,natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          arr = RESHAPE(array,(/3,natom/))
       END IF
       NULLIFY(x4,y4,z4)
       ALLOCATE(x4(natom),y4(natom),z4(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_cell(cell, alpha=alpha, beta=beta, gamma=gamma, abc=abc)
       WRITE(iunit) abc(1)*factor, gamma, abc(2)*factor, beta, alpha, abc(3)*factor
       SELECT CASE (TRIM(content))
       CASE ("POS")
          IF (PRESENT(array)) THEN
             x4(:) = arr(1,:)
             y4(:) = arr(2,:)
             z4(:) = arr(3,:)
          ELSE
             x4(:) = particle_set(:)%r(1)
             y4(:) = particle_set(:)%r(2)
             z4(:) = particle_set(:)%r(3)
          END IF
          WRITE(iunit) x4*REAL(factor,KIND=sp)
          WRITE(iunit) y4*REAL(factor,KIND=sp)
          WRITE(iunit) z4*REAL(factor,KIND=sp)
       CASE ("VEL")
          IF (PRESENT(array)) THEN
             x4(:) = arr(1,:)
             y4(:) = arr(2,:)
             z4(:) = arr(3,:)
          ELSE
             x4(:) = particle_set(:)%v(1)
             y4(:) = particle_set(:)%v(2)
             z4(:) = particle_set(:)%v(3)
          END IF
          WRITE(iunit) x4*REAL(factor,KIND=sp)
          WRITE(iunit) y4*REAL(factor,KIND=sp)
          WRITE(iunit) z4*REAL(factor,KIND=sp)
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
          WRITE(iunit) x4*REAL(factor,KIND=sp)
          WRITE(iunit) y4*REAL(factor,KIND=sp)
          WRITE(iunit) z4*REAL(factor,KIND=sp)
       CASE DEFAULT
          CALL stop_program ("particle_types","illegal dcd dump type")
       END SELECT
       IF (PRESENT(array)) THEN
          DEALLOCATE(arr,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(x4,y4,z4,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE DEFAULT
       CALL stop_program ("particle_types","illegal dump type")
    END SELECT

  END SUBROUTINE write_particle_coordinates

! *****************************************************************************
!> \brief   Write the atomic coordinates to the output unit.
!> \author  MK
!> \date    05.06.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_fist_particle_coordinates(particle_set,cell,subsys_section,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: subsys_section
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
  SUBROUTINE write_qs_particle_coordinates(particle_set,cell,subsys_section,label,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: subsys_section
    CHARACTER(LEN=*), INTENT(IN)             :: label
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: unit_str
    INTEGER                                  :: iatom, ikind, iw, natom, z
    REAL(KIND=dp)                            :: conv, mass, zeff
    TYPE(cp_logger_type), POINTER            :: logger

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
               FMT="(T2,I5,1X,I4,3X,A2,2X,I3,3F12.6,4X,F6.2,2X,F11.4)")&
          iatom,ikind,element_symbol,z,&
          particle_set(iatom)%r(1:3)*conv,zeff,mass/massunit
      END DO
      WRITE (iw,'(/)')
    END IF

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%ATOMIC_COORDINATES", error=error)

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
    INTEGER                                  :: iatom, iw, jatom, natom, stat
    INTEGER, DIMENSION(3)                    :: periodic
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: conv, dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: distance_matrix
    REAL(KIND=dp), DIMENSION(3)              :: rab
    TYPE(cp_logger_type), POINTER            :: logger

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

       !     *** Print the distance matrix ***
       WRITE (UNIT=iw,FMT="(/,/,T2,A)")&
            "INTERATOMIC DISTANCES IN "//TRIM(unit_str)

       CALL write_particle_matrix(distance_matrix,particle_set,iw,error=error)
    END IF

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%INTERATOMIC_DISTANCES", error=error)
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
    INTEGER                                  :: handle, i, iw, n, natom, &
                                                wrk1(1), wrk2(2), wrk3(3), &
                                                wrk4(4)
    INTEGER, DIMENSION(:), POINTER           :: atomic_indices
    LOGICAL                                  :: failure, unique
    REAL(KIND=dp)                            :: conv, dab
    REAL(KIND=dp), DIMENSION(3)              :: rab, rbc, rcd
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY (atomic_indices)
    NULLIFY (logger)
    NULLIFY (section)

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger=logger,&
                              basis_section=input_section,&
                              print_key_path="PRINT%STRUCTURE_DATA",&
                              extension=".coordLog",&
                              error=error)

    CALL section_vals_val_get(input_section,"PRINT%STRUCTURE_DATA%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (iw > 0) THEN
      natom = SIZE(particle_set)
      section => section_vals_get_subs_vals(section_vals=input_section,&
                                            subsection_name="PRINT%STRUCTURE_DATA",&
                                            error=error)

      WRITE (UNIT=iw,FMT="(/,T2,A)") "REQUESTED STRUCTURE DATA"
      ! Print the requested atomic position vectors
      CALL section_vals_val_get(section_vals=section,&
                                keyword_name="POSITION",&
                                n_rep_val=n,&
                                error=error)
      IF (n > 0) THEN
        WRITE (UNIT=iw,FMT="(/,T3,A,/)")&
          "Position vector r(i) of atom i in "//TRIM(unit_str)
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
          wrk1=atomic_indices
          CALL sort_unique(wrk1, unique)
          IF (((wrk1(1)>=1).AND.(wrk1(SIZE(wrk1))<=natom)).AND.unique) THEN
            WRITE (UNIT=iw,FMT="(T3,A,T20,A,3F13.6)")&
              "r"//TRIM(string),"=",&
              pbc(particle_set(atomic_indices(1))%r(1:3),cell)*conv
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
          TRIM(unit_str)
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
          WRITE (UNIT=string,FMT="(A,3(I6,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),",",atomic_indices(3),")"
          CALL compress(string,full=.TRUE.)
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
          WRITE (UNIT=string,FMT="(A,4(I6,A))")&
            "(",atomic_indices(1),",",atomic_indices(2),",",&
                atomic_indices(3),",",atomic_indices(4),")"
          CALL compress(string,full=.TRUE.)
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

    INTEGER                                  :: iparticle, nparticle
    LOGICAL                                  :: failure, my_add, update_for, &
                                                update_pos, update_vel

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
    ENDIF
    
  END SUBROUTINE update_particle_set

END MODULE particle_types
