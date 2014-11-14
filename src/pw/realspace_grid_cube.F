!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Generate Gaussian cube files
! *****************************************************************************
MODULE realspace_grid_cube
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_maxloc,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pw_to_cube, cube_to_pw

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_cube'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param stride ...
!> \param zero_tails ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pw_to_cube ( pw, unit_nr, title, particles_r, particles_z, stride, zero_tails, error )
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: unit_nr
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: particles_r
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: particles_z
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: stride
    LOGICAL, INTENT(IN), OPTIONAL            :: zero_tails
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_cube', &
      routineP = moduleN//':'//routineN

    INTEGER :: checksum, dest, gid, handle, i, I1, I2, I3, iat, ip, L1, L2, &
      L3, my_rank, my_stride(3), np, num_pe, rank(2), source, stat, tag, U1, &
      U2, U3
    LOGICAL                                  :: failure, my_zero_tails
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buf

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

      CPPrecondition(PRESENT(particles_z) .EQV. PRESENT(particles_r),cp_failure_level,routineP,error,failure)
      np=0
      IF(PRESENT(particles_z)) THEN
         CPPrecondition(SIZE(particles_z) == SIZE(particles_r, dim=2),cp_failure_level,routineP,error,failure)
         ! cube files can only be written for 99999 particles due to a format limitation (I5)
         ! so we limit the number of particles written.
         np=MIN(99999,SIZE(particles_z))
      END IF

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

      IF(PRESENT(particles_z)) THEN
         DO iat=1,np
            WRITE(unit_nr,'(I5,4f12.6)') particles_z(iat),0._dp,particles_r(:,iat)
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
  
  
! *****************************************************************************
!> \brief  Computes the external density on the grid 
!>         hacked from external_read_density
!> \param grid     pw to read from cube file
!> \param filename name of cube file
!> \param scaling  scale values before storing
!> \param error ...
!> \date   01.2014
!> \author M.Watkins
! *****************************************************************************
  SUBROUTINE cube_to_pw(grid, filename, scaling, error)

    TYPE(pw_type), POINTER                   :: grid
    CHARACTER(len=*), INTENT(in)             :: filename
    REAL(kind=dp), INTENT(in)                :: scaling
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_to_pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extunit, gid, handle, i, ip, &
                                                j, k, master, my_rank, nat, &
                                                ndum, num_pe, tag
    INTEGER, DIMENSION(3)                    :: lbounds, lbounds_local, &
                                                npoints, npoints_local, &
                                                ubounds, ubounds_local
    LOGICAL                                  :: failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: buffer
    REAL(kind=dp), DIMENSION(3)              :: dr, rdum

    failure=.FALSE.
    CALL timeset(routineN,handle)

    !get rs grids and parallel envnment
    gid     = grid%pw_grid%para%group
    my_rank = grid%pw_grid%para%my_pos
    num_pe  = grid%pw_grid%para%group_size   
    tag     = 1

    DO i =1,3
       dr(i) = grid%pw_grid%dh(i,i)  
    END DO
    
   npoints = grid%pw_grid%npts
   lbounds = grid%pw_grid%bounds(1,:)
   ubounds = grid%pw_grid%bounds(2,:)

   npoints_local = grid%pw_grid%npts_local  
   lbounds_local = grid%pw_grid%bounds_local(1,:)
   ubounds_local = grid%pw_grid%bounds_local(2,:)

   !pw grids at most pencils - all processors have a full set of z data for x,y 
   ALLOCATE(buffer(lbounds(3):ubounds(3)))

   IF ( my_rank == 0 ) THEN
      WRITE(*,FMT="(/,(T3,A,T61,A30))")"Reading the cube file:     ",filename

      CALL open_file(file_name=filename,&
           file_status="OLD",&
           file_form="FORMATTED",&
           file_action="READ", &
           unit_number=extunit)
      
      !skip header comments
      DO i =1,2
         READ(extunit,*) 
      ENDDO
      READ(extunit,*) nat,rdum
      DO i=1,3
         READ(extunit,*) ndum,rdum
         IF (ndum/=npoints(i).OR.(ABS(rdum(i)-dr(i)) > 1e-4 )) THEN   
            WRITE(*,*) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
            WRITE(*,*) "Restart from density | ",ndum," DIFFERS FROM ", npoints(i) 
            WRITE(*,*) "Restart from density | ",rdum," DIFFERS FROM ", dr(i)
         ENDIF
      ENDDO
      !ignore atomic postion data - read from coord or topology instead
      DO i =1,nat
         READ(extunit,*)
      ENDDO
   ENDIF

   !master sends all data to everyone
   DO i=lbounds(1),ubounds(1)
      DO j=lbounds(2),ubounds(2) 
         IF ( my_rank .EQ. 0 ) THEN 
            READ(extunit,*) (buffer(k),k=lbounds(3),ubounds(3))   
            IF ( num_pe .GT. 1 ) THEN 
               DO ip=1,num_pe-1
                  CALL mp_send(buffer(lbounds(3):ubounds(3)),ip,tag,gid)
               ENDDO
            ENDIF
         ELSE
            master=0  
            CALL mp_recv(buffer(lbounds(3):ubounds(3)),master,tag,gid)
         ENDIF
         CALL mp_sync(gid)
         
         !only use data that is local to me - i.e. in slice of pencil I own
         IF ( (lbounds_local(1).le.i).AND.(i.le.ubounds_local(1)).AND.(lbounds_local(2).le.j)&
              .AND.(j.le.ubounds_local(2)) ) THEN
            !allow scaling of external potential values by factor 'scaling' (SCALING_FACTOR in input file)
            grid%cr3d(i,j,lbounds(3):ubounds(3)) = buffer(lbounds(3):ubounds(3))*scaling
         ENDIF

      ENDDO
   ENDDO

   IF (my_rank == 0) CALL close_file(unit_number=extunit)
       
   CALL mp_sync(gid)

   CALL timestop(handle)
    
  END SUBROUTINE cube_to_pw

END MODULE realspace_grid_cube
