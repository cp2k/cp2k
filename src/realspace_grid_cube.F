!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Generate Gaussian cube files
! *****************************************************************************
MODULE realspace_grid_cube
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_maxloc,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pw_to_cube

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_cube'

CONTAINS

  SUBROUTINE pw_to_cube ( pw, unit_nr, title, particles, stride, zero_tails, error )
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: unit_nr
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: stride
    LOGICAL, INTENT(IN), OPTIONAL            :: zero_tails
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_cube', &
      routineP = moduleN//':'//routineN

    INTEGER :: checksum, dest, gid, handle, i, I1, I2, I3, iat, ip, L1, L2, &
      L3, my_rank, my_stride(3), np, num_pe, rank(2), source, stat, tag, U1, &
      U2, U3, z
    LOGICAL                                  :: failure, my_zero_tails
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buf
    TYPE(particle_list_type), POINTER        :: my_particles

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_zero_tails = .FALSE.
    IF(PRESENT(zero_tails)) my_zero_tails = zero_tails
    my_stride=1
    IF (PRESENT(stride)) THEN
       CALL cp_assert(SIZE(stride)==1.OR.SIZE(stride)==3,cp_fatal_level,cp_assertion_failed,routineP,&
            "STRIDE keyword can accept only 1 (the same for X,Y,Z) or 3 values. Correct your input file."//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
       IF (SIZE(stride)==1) THEN
          DO i = 1,3
             my_stride(i) = stride(1)
          END DO
       ELSE
          my_stride = stride(1:3)
       END IF
       CPPrecondition(my_stride(1)>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(my_stride(2)>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(my_stride(3)>0,cp_failure_level,routineP,error,failure)
    END IF

    IF (unit_nr>0) THEN
      ! this format seems to work for e.g. molekel and gOpenmol
      ! latest version of VMD can read non orthorhombic cells
      WRITE(unit_nr,'(a11)') "-Quickstep-"
      IF (PRESENT(title)) THEN
         WRITE(unit_nr,*) TRIM(title)
      ELSE
         WRITE(unit_nr,*) "No Title"
      ENDIF

      np=0
      NULLIFY(my_particles)
      my_particles=>particles
      ! cube files can only be written for 99999 particles due to a format limitation (I5)
      ! so we limit the number of particles written.
      IF (ASSOCIATED(my_particles)) np=MIN(99999,my_particles%n_els)

      WRITE(unit_nr,'(I5,3f12.6)') np,0.0_dp,0._dp,0._dp !start of cube

      WRITE(unit_nr,'(I5,3f12.6)') (pw%pw_grid%npts(1)+my_stride(1)-1)/my_stride(1),&
           pw%pw_grid%dh(1,1)*REAL(my_stride(1),dp),pw%pw_grid%dh(2,1)*REAL(my_stride(1),dp),&
           pw%pw_grid%dh(3,1)*REAL(my_stride(1),dp)
      WRITE(unit_nr,'(I5,3f12.6)') (pw%pw_grid%npts(2)+my_stride(2)-1)/my_stride(2),&
           pw%pw_grid%dh(1,2)*REAL(my_stride(2),dp),pw%pw_grid%dh(2,2)*REAL(my_stride(2),dp),&
           pw%pw_grid%dh(3,2)*REAL(my_stride(2),dp)
      WRITE(unit_nr,'(I5,3f12.6)') (pw%pw_grid%npts(3)+my_stride(3)-1)/my_stride(3),&
           pw%pw_grid%dh(1,3)*REAL(my_stride(3),dp),pw%pw_grid%dh(2,3)*REAL(my_stride(3),dp),&
           pw%pw_grid%dh(3,3)*REAL(my_stride(3),dp)

      IF (ASSOCIATED(my_particles)) THEN
         DO iat=1,np
            CALL get_atomic_kind(my_particles%els(iat)%atomic_kind,z=z)

            WRITE(unit_nr,'(I5,4f12.6)') z,0._dp,my_particles%els(iat)%r
         END DO
      END IF
    ENDIF

    ! shortcut
    L1=pw%pw_grid%bounds(1,1)
    L2=pw%pw_grid%bounds(1,2)
    L3=pw%pw_grid%bounds(1,3)
    U1=pw%pw_grid%bounds(2,1)
    U2=pw%pw_grid%bounds(2,2)
    U3=pw%pw_grid%bounds(2,3)

    ALLOCATE(buf(L3:U3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    my_rank     = pw%pw_grid%para%my_pos
    gid         = pw%pw_grid%para%group
    num_pe      = pw%pw_grid%para%group_size
    tag         = 1

    rank(1)=unit_nr
    rank(2)=my_rank
    checksum=0
    IF (unit_nr>0) checksum=1

    CALL mp_sum(checksum,gid)
    CPPostcondition(checksum==1,cp_failure_level,routineP,error,failure)

    CALL mp_maxloc(rank,gid)
    CPPostcondition(rank(1)>0,cp_failure_level,routineP,error,failure)

    dest = rank(2)
    DO I1=L1,U1,my_stride(1)
      DO I2=L2,U2,my_stride(2)

       ! cycling through the CPUs, check if the current ray (I1,I2) is local to that CPU
       IF (pw%pw_grid%para%mode .NE. PW_MODE_LOCAL ) THEN
         DO ip=0,num_pe-1
            IF (pw%pw_grid%para%bo (1,1,ip,1)<=I1-L1+1 .AND. pw%pw_grid%para%bo (2,1,ip,1)>=I1-L1+1 .AND. &
                pw%pw_grid%para%bo (1,2,ip,1)<=I2-L2+1 .AND. pw%pw_grid%para%bo (2,2,ip,1)>=I2-L2+1 ) THEN
                source = ip
            ENDIF
         ENDDO
       ELSE
         source = dest
       ENDIF

       IF (source==dest) THEN
          IF (my_rank==source) THEN
             buf(:)=pw%cr3d(I1,I2,:)
          ENDIF
       ELSE
          IF(my_rank==source) THEN
             buf(:)=pw%cr3d(I1,I2,:)
             CALL mp_send(buf,dest,tag,gid)
          ENDIF
          IF(my_rank==dest) THEN
             CALL mp_recv(buf,source,tag,gid)
          ENDIF
       ENDIF

       IF (unit_nr>0) THEN
         IF(my_zero_tails) THEN
           DO I3 = L3,U3
             IF(buf(I3)<1.E-7_dp) buf(I3) = 0.0_dp
           END DO
         END IF
         WRITE(unit_nr,'(6E13.5)') (buf(I3),I3=L3,U3,my_stride(3))
       ENDIF

       ! this double loop generates so many messages that it can overload
       ! the message passing system, e.g. on XT3
       ! we therefore put a barrier here that limits the amount of message
       ! that flies around at any given time.
       ! if ever this routine becomes a bottleneck, we should go for a
       ! more complicated rewrite
       CALL mp_sync(gid)

      ENDDO
    ENDDO

    DEALLOCATE(buf)

    CALL timestop(handle)

  END SUBROUTINE pw_to_cube

END MODULE realspace_grid_cube
