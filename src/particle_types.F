!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
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
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_bondparm_covalent,&
                                             do_bondparm_vdw,&
                                             dump_atomic,&
                                             dump_dcd,&
                                             dump_xmol
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size,&
                                             sp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             massunit
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: find_boundary
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'particle_types'

! *** Data types ***

  TYPE exclusion_list_type
    INTEGER, DIMENSION(:), POINTER :: list_vdw
    INTEGER, DIMENSION(:), POINTER :: list_ei
  END TYPE exclusion_list_type

  TYPE local_particle_type
    TYPE(exclusion_list_type), DIMENSION(:), POINTER :: exclusion
  END TYPE local_particle_type

  TYPE particle_type
     TYPE(atomic_kind_type), POINTER       :: atomic_kind ! atomic kind info
     TYPE(local_particle_type), POINTER    :: lpi         ! local particle info
     REAL(KIND = dp), DIMENSION(3)         :: f,&         ! force
                                              r,&         ! position
                                              v           ! velocity
! particle dependent terms for EAM
     REAL(KIND = dp)                       :: rho, f_embed
     LOGICAL                               :: is_fixed    ! force is set to zero
  END TYPE particle_type

! *** Public data types ***

  PUBLIC :: exclusion_list_type,&
            local_particle_type,&
            particle_type

! *** Public subroutines ***

  PUBLIC :: allocate_particle_set,&
            deallocate_particle_set,&
            generate_particle_distances,&
            get_particle_set,&
            write_fist_particle_coordinates,&
            write_qs_particle_coordinates,&
            write_particle_distances, &
            write_particle_coordinates,&
            clone_particle

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'clone_particle', &
      routineP = moduleN//':'//routineN

  NULLIFY(p_out%atomic_kind)
  NULLIFY(p_out%lpi)
  IF (.NOT.nullify_non_owned_pointers) THEN
     p_out%atomic_kind => p_in%atomic_kind
     p_out%lpi => p_in%lpi
  END IF
  p_out%r=p_in%r
  p_out%f=p_in%f
  p_out%v=p_in%v
  p_out%is_fixed=p_in%is_fixed
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
      NULLIFY (particle_set(iparticle)%lpi)
      particle_set(iparticle)%f(:) = 0.0_dp
      particle_set(iparticle)%r(:) = 0.0_dp
      particle_set(iparticle)%v(:) = 0.0_dp
      particle_set(iparticle)%is_fixed = .FALSE.
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

    INTEGER                                  :: i, iparticle, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) THEN

      DO iparticle=1,SIZE(particle_set)

        IF (ASSOCIATED(particle_set(iparticle)%lpi)) THEN
          IF (ASSOCIATED(particle_set(iparticle)%lpi%exclusion)) THEN
            DO i=1,SIZE(particle_set(iparticle)%lpi%exclusion)
              IF (ASSOCIATED(particle_set(iparticle)%lpi%exclusion(i)%list_vdw)) THEN
                DEALLOCATE (particle_set(iparticle)%lpi%exclusion(i)%list_vdw,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                 "particle_set(iparticle)%exclusion(i)%list_vdw")
              END IF
              IF (ASSOCIATED(particle_set(iparticle)%lpi%exclusion(i)%list_ei)) THEN
                DEALLOCATE (particle_set(iparticle)%lpi%exclusion(i)%list_ei,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                 "particle_set(iparticle)%exclusion(i)%list_ei")
              END IF
            END DO
            DEALLOCATE (particle_set(iparticle)%lpi%exclusion,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routineN,moduleN,__LINE__,&
                               "particle_set(iparticle)%lpi%exclusion")
            END IF
          END IF
          DEALLOCATE (particle_set(iparticle)%lpi,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routineN,moduleN,__LINE__,&
                             "particle_set(iparticle)%lpi")
          END IF
        END IF
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

  SUBROUTINE write_particle_coordinates(particle_set,cell,iunit,FORMAT,&
                                        content,title)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: iunit, FORMAT
    CHARACTER(LEN=*)                         :: content, title

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=4)                         :: name
    CHARACTER(LEN=default_string_length)     :: atm_name, my_format
    INTEGER                                  :: iatom, istat, natom
    LOGICAL                                  :: dummy
    REAL(KIND=dp), DIMENSION(3)              :: r
    REAL(KIND=sp), DIMENSION(:), POINTER     :: x4, y4, z4

    natom = SIZE(particle_set)

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
          r(:) = particle_set(iatom)%r(:)*angstrom
          WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),r(1:3)
        CASE ("VEL")
          WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),&
                                            particle_set(iatom)%v(1:3)
        CASE ("FORCE")
          WRITE (iunit,TRIM(my_format)//"1X,3F20.10)") TRIM(name),&
                                            particle_set(iatom)%f(1:3)
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
         r(:) = particle_set(iatom)%r(:)
         WRITE (iunit,"(3F20.10)") r(1:3)
        CASE ("VEL")
          WRITE (iunit,"(3F20.10)")particle_set(iatom)%v(1:3)
        CASE ("FORCE")
          WRITE (iunit,"(3F20.10)") particle_set(iatom)%f(1:3)
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
       NULLIFY(x4,y4,z4)
       ALLOCATE(x4(natom),y4(natom),z4(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory('particle_types','x4,y4,z4')
       SELECT CASE (TRIM(content))
       CASE ("POS")
         x4(:) = particle_set(:)%r(1)*angstrom
         y4(:) = particle_set(:)%r(2)*angstrom
         z4(:) = particle_set(:)%r(3)*angstrom
         WRITE(iunit) x4
         WRITE(iunit) y4
         WRITE(iunit) z4
       CASE ("VEL")
         x4(:) = particle_set(:)%v(1)*angstrom*femtoseconds
         y4(:) = particle_set(:)%v(2)*angstrom*femtoseconds
         z4(:) = particle_set(:)%v(3)*angstrom*femtoseconds
         WRITE(iunit) x4
         WRITE(iunit) y4
         WRITE(iunit) z4
        CASE ("FORCE")
         x4(:) = particle_set(:)%f(1)
         y4(:) = particle_set(:)%f(2)
         z4(:) = particle_set(:)%f(3)
         WRITE(iunit) x4
         WRITE(iunit) y4
         WRITE(iunit) z4
       CASE DEFAULT
         CALL stop_program ("particle_types","illegal dcd dump type")
       END SELECT
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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

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
        IF (particle_set(iatom)%is_fixed) THEN
          WRITE (iw,&
               "(T2,I5,1X,I4,3X,A4,3X,3F12.6,4X,F6.2,2X,F11.4,1X,A1)")&
            iatom,ikind,name,&
            particle_set(iatom)%r(1:3)/unit_of_length,qeff,mass/massunit,"!"
        ELSE
          WRITE (iw,&
               "(T2,I5,1X,I4,3X,A4,3X,3F12.6,4X,F6.2,2X,F11.4)")&
            iatom,ikind,name,&
            particle_set(iatom)%r(1:3)/unit_of_length,qeff,mass/massunit
        END IF
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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

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
        IF (particle_set(iatom)%is_fixed) THEN
          WRITE (iw,&
               "(T2,I5,1X,I4,3X,A2,2X,I3,3F12.6,4X,F6.2,2X,F11.4,1X,A1)")&
            iatom,ikind,element_symbol,z,&
            particle_set(iatom)%r(1:3)/unit_of_length,zeff,mass/massunit,"!"
        ELSE
          WRITE (iw,&
               "(T2,I5,1X,I4,3X,A2,2X,I3,3F12.6,4X,F6.2,2X,F11.4)")&
            iatom,ikind,element_symbol,z,&
            particle_set(iatom)%r(1:3)/unit_of_length,zeff,mass/massunit
        END IF
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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_particle_distances'

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER                                  :: from, iatom, icol, ikind, &
                                                istat, iw, jatom, natom, to
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

      DO jatom=1,natom,5
        from = jatom
        to = MIN(from+4,natom)
        WRITE (UNIT=iw,FMT="(/,T2,11X,5(4X,I5,4X))")&
          (icol,icol=from,to)
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               kind_number=ikind,&
                               element_symbol=element_symbol)
          WRITE (UNIT=iw,FMT="(T2,I5,2X,A2,2X,5(1X,F12.6))")&
            iatom,element_symbol,&
            (distance_matrix(iatom,icol),icol=from,to)
        END DO
      END DO

    END IF

    CALL cp_print_key_finished_output(iw,logger,force_env_section,&
         "SUBSYS%PRINT%INTERATOMIC_DISTANCES", error=error)
  END SUBROUTINE write_particle_distances

! *****************************************************************************

SUBROUTINE generate_particle_distances(r,cell,radius,bondparm_factor,&
                                       bond_a,bond_b,bondparm_type,globenv)

    REAL(KIND=dp), POINTER                   :: r(:,:)
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), POINTER                   :: radius(:)
    REAL(KIND=dp)                            :: bondparm_factor
    INTEGER, POINTER                         :: bond_a(:), bond_b(:)
    INTEGER                                  :: bondparm_type
    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER, PARAMETER                       :: nblock = 1000

    INTEGER                                  :: first, i, iatom, iproc, &
                                                istat, j, jatom, k, last, &
                                                natom, nbond, workload
    INTEGER, POINTER                         :: nbond_loc(:)
    REAL(KIND=dp)                            :: dr(3), r2, rbond, rbond2, &
                                                rr(3)
    REAL(KIND=dp), POINTER                   :: tbond_a(:), tbond_b(:), &
                                                ttbond_a(:), ttbond_b(:)

!------------------------------------------------------------------------------

  NULLIFY(bond_a,bond_b,tbond_a,tbond_b,ttbond_a,ttbond_b)
  NULLIFY(nbond_loc)
!------------------------------------------------------------------------------
  CALL reallocate(tbond_a,1,nblock)
  CALL reallocate(tbond_b,1,nblock)

  workload = 0
  nbond = 0
  natom = SIZE(r,2)
  DO iatom=1,natom
    DO jatom=iatom+1+globenv%mepos,natom,globenv%num_pe
      workload = workload + 1
      dr(1) = r(1,iatom) - r(1,jatom)
      dr(2) = r(2,iatom) - r(2,jatom)
      dr(3) = r(3,iatom) - r(3,jatom)

      rr = pbc ( dr, cell )
      r2 = DOT_PRODUCT(rr,rr)

      IF(r2 <= 0.01_dp) THEN
        CALL stop_program ("topology_generate_bond",&
                           "bond distance between atoms less then 0.1")
      END IF

      IF(bondparm_type==do_bondparm_covalent) THEN
        rbond = radius(iatom) + radius(jatom)
      ELSE IF(bondparm_type==do_bondparm_vdw) THEN
        rbond = MAX(radius(iatom),radius(jatom))
      END IF

      rbond2 = rbond*rbond

      !Look for the right factor to use...
      rbond2 = rbond2*(bondparm_factor)

      !Test the distance to the sum of the covalent radius
      IF(r2 < rbond2) THEN
        nbond=nbond+1
        IF(nbond > SIZE(tbond_a)) THEN
          CALL reallocate(tbond_a,1,nbond+nblock)
          CALL reallocate(tbond_b,1,nbond+nblock)
        END IF
        tbond_a(nbond)=iatom
        tbond_b(nbond)=jatom
      END IF
    END DO
  END DO
  CALL reallocate(tbond_a,1,nbond)
  CALL reallocate(tbond_b,1,nbond)
!------------------------------------------------------------------------------
  !Put together the difference pieces for the bond
  
  !Determine the total number of bonds
  ALLOCATE(nbond_loc(globenv%num_pe),STAT=istat)
  IF(istat/=0) CALL stop_program('topology_generate_bond','nbond_loc')
  nbond_loc(:)=0
  CALL mp_allgather(nbond,nbond_loc,globenv%group)
  nbond = SUM(nbond_loc)

  !Reallocate those memory
  ALLOCATE(bond_a(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_a',nbond)
  ALLOCATE(bond_b(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_b',nbond)
  bond_a(:)=0
  bond_b(:)=0

  workload=0
  !Gather everything
  DO iproc=0,globenv%num_pe-1
    CALL reallocate(ttbond_a,1,nbond_loc(iproc+1))
    CALL reallocate(ttbond_b,1,nbond_loc(iproc+1))
    IF(iproc==globenv%mepos) THEN
      ttbond_a(:) = tbond_a(:)
      ttbond_b(:) = tbond_b(:)
    ELSE
      ttbond_a(:) = 0
      ttbond_b(:) = 0
    ENDIF
    CALL mp_bcast(ttbond_a,iproc,globenv%group)
    CALL mp_bcast(ttbond_b,iproc,globenv%group)
    bond_a(workload+1:workload+nbond_loc(iproc+1)) = ttbond_a(1:nbond_loc(iproc+1))
    bond_b(workload+1:workload+nbond_loc(iproc+1)) = ttbond_b(1:nbond_loc(iproc+1))
    workload=workload+nbond_loc(iproc+1)
  ENDDO

  !Sort the array
  !Sort the first column
  DO i=1,nbond
    DO j=1,nbond
      IF(bond_a(i) < bond_a(j)) THEN
        workload  = bond_a(i)
        bond_a(i) = bond_a(j)
        bond_a(j) = workload
        workload  = bond_b(i)
        bond_b(i) = bond_b(j)
        bond_b(j) = workload
      END IF
    END DO
  END DO
  !Sort the second column
  i=1
  DO WHILE (i<nbond)
    CALL find_boundary(bond_a,nbond,first,last,bond_a(i))
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(bond_b(j) < bond_b(k)) THEN
            workload = bond_a(j)
            bond_a(j) = bond_a(k)
            bond_a(k) = workload
            workload = bond_b(j)
            bond_b(j) = bond_b(k)
            bond_b(k) = workload
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO
!------------------------------------------------------------------------------
  DEALLOCATE(nbond_loc,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','nbond_loc')
  DEALLOCATE(tbond_a,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_a')
  DEALLOCATE(tbond_b,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_b')
  DEALLOCATE(ttbond_a,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','ttbond_a')
  DEALLOCATE(ttbond_b,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','ttbond_b')
!------------------------------------------------------------------------------
END SUBROUTINE generate_particle_distances

END MODULE particle_types
