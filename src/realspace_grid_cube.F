!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****s* cp2k/realspace_grid_cube [1.0] *
!!
!!   NAME
!!     realspace_grid_cube
!!
!!   FUNCTION
!!     Generate Gaussian cube files
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!
!!***
!******************************************************************************

MODULE realspace_grid_cube
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_types,                        ONLY: COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             pw_type
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_create,&
                                             rs_grid_release,&
                                             rs_pw_transfer,&
                                             pw2rs
  USE timings,                         ONLY: timeset,&
                                             timestop

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rs_grid_to_cube, &
            rs_pw_to_cube

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_cube'

!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
SUBROUTINE rs_grid_to_cube ( rs, iunit, ionode, title, particles, stride, error )
    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, INTENT(IN)                      :: iunit
    LOGICAL, INTENT(IN)                      :: ionode
    CHARACTER(*), OPTIONAL                   :: title
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    INTEGER, INTENT(in), OPTIONAL            :: stride
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_to_cube', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I1, I2, I3, iat, L1, L2, L3, &
                                                my_stride, np, U1, U2, U3, z
    LOGICAL                                  :: failure
    TYPE(particle_list_type), POINTER        :: my_particles

!-----------------------------------------------------------------------------!
! not working for distributed grids !

   failure=.FALSE.
   my_stride=1
   IF (PRESENT(stride)) my_stride=stride
   CPPrecondition(my_stride>0,cp_failure_level,routineP,error,failure)

   ! only dh diagonal is being used
   CPPrecondition(rs%orthorhombic,cp_failure_level,routineP,error,failure)
   CPPrecondition(rs%ngpts == rs%ngpts_local,cp_failure_level,routineP,error,failure)
   IF (ionode) THEN
     ! this format seems to work for e.g. molekel and gOpenmol
     WRITE(iunit,'(a11)') "-Quickstep-"
     IF (PRESENT(title)) THEN
        WRITE(iunit,*) TRIM(title)
     ELSE
        WRITE(iunit,*) "No Title"
     ENDIF

     np=0
     NULLIFY(my_particles)
     IF (PRESENT(particles)) my_particles=>particles
     IF (ASSOCIATED(my_particles)) np=my_particles%n_els

     WRITE(iunit,'(I5,3f12.6)') np,0.0_dp,0._dp,0._dp !start of cube

     WRITE(iunit,'(I5,3f12.6)') (rs%npts(1)+my_stride-1)/my_stride,&
          rs%dh(1,1)*REAL(my_stride,dp),0.0_dp,0.0_dp
     WRITE(iunit,'(I5,3f12.6)') (rs%npts(2)+my_stride-1)/my_stride,&
          0.0_dp,rs%dh(2,2)*REAL(my_stride,dp),0.0_dp
     WRITE(iunit,'(I5,3f12.6)') (rs%npts(3)+my_stride-1)/my_stride,&
          0.0_dp,0.0_dp,rs%dh(3,3)*REAL(my_stride,dp)

     IF (ASSOCIATED(my_particles)) THEN
        DO iat=1,np
           CALL get_atomic_kind(my_particles%els(iat)%atomic_kind,z=z)

           WRITE(iunit,'(I5,4f12.6)') z,0._dp,my_particles%els(iat)%r
        END DO
     END IF

     ! write(iunit,'(I5,4f12.6)') 1, 0.0,0.0,0.0,0.0
     L1=LBOUND(rs%r,1)
     L2=LBOUND(rs%r,2)
     L3=LBOUND(rs%r,3)
     U1=UBOUND(rs%r,1)
     U2=UBOUND(rs%r,2)
     U3=UBOUND(rs%r,3)
     DO I1=L1,U1,my_stride
       DO I2=L2,U2,my_stride
        WRITE(iunit,'(6E13.5)') (rs%r(MODULO(I1-L1,U1-L1+1)+L1,MODULO(I2-L2,U2-L2+1)+L2, &
                                      MODULO(I3-L3,U3-L3+1)+L3),I3=L3,U3,my_stride)
       ENDDO
     ENDDO
   ENDIF
END SUBROUTINE rs_grid_to_cube
!******************************************************************************

!******************************************************************************
! optionally pass a preallocated rs grid
!******************************************************************************
SUBROUTINE rs_pw_to_cube ( pw, iunit, ionode, title, rs, particles, stride,&
     error )
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: iunit
    LOGICAL, INTENT(IN)                      :: ionode
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    TYPE(realspace_grid_type), OPTIONAL, &
      POINTER                                :: rs
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    INTEGER, INTENT(in), OPTIONAL            :: stride
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: handle, nsmax
    TYPE(realspace_grid_type), POINTER       :: rs_local

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

   CALL timeset("rs_pw_to_cube",handle)
   IF (PRESENT(rs)) THEN
     rs_local=>rs
   ELSE
     nsmax=-1
     NULLIFY(rs_local)
     CALL rs_grid_create(rs_local,pw%pw_grid,nsmax,error=error)
   ENDIF

   CALL rs_pw_transfer(rs_local,pw,pw2rs)
   IF (PRESENT(particles)) THEN
      CALL rs_grid_to_cube(rs_local,iunit,ionode,title,particles=particles,&
           stride=stride,error=error)
   ELSE
      CALL rs_grid_to_cube(rs_local,iunit,ionode,title,stride=stride, error=error)
   END IF

   IF (.NOT. PRESENT(rs)) THEN
       CALL rs_grid_release(rs_local,error=error)
   ENDIF
   CALL timestop(handle)

END SUBROUTINE rs_pw_to_cube

!******************************************************************************

END MODULE realspace_grid_cube

!******************************************************************************
