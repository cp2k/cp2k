!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (30-Nov-2000): ESSL FFT Library added
!>      JGH (05-Jan-2001): Added SGI library FFT
!>      JGH (14-Jan-2001): Added parallel 3d FFT
!>      JGH (10-Feb-2006): New interface type
!>      JGH (31-Mar-2008): Remove local allocates and reshapes (performance)
!>                         Possible problems can be related with setting arrays
!>                         not to zero
!>                         Some interfaces could be further simplified by avoiding
!>                         an initial copy. However, this assumes contiguous arrays
!>      IAB (15-Oct-2008): Moved mp_cart_sub calls out of cube_tranpose_* and into
!>                         fft_scratch type, reducing number of calls dramatically
!>      IAB (05-Dec-2008): Moved all other non-essential MPI calls into scratch type
!>      IAB (09-Jan-2009): Added fft_plan_type to store FFT data, including cached FFTW plans
!>      IAB (13-Feb-2009): Extended plan caching to serial 3D FFT (fft3d_s)
!>      IAB (09-Oct-2009): Added OpenMP directives to parallel 3D FFT
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2008-2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************
MODULE fft_tools

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fast,                            ONLY: zero_c
  USE fft_lib,                         ONLY: &
       fft_1dm, fft_3d, fft_create_plan_1dm, fft_create_plan_3d, &
       fft_destroy_plan, fft_do_cleanup, fft_do_init, fft_get_lengths, &
       fft_library
  USE fft_plan,                        ONLY: fft_plan_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size,&
                                             sp
  USE message_passing,                 ONLY: &
       mp_alltoall, mp_cart_coords, mp_cart_rank, mp_cart_sub, &
       mp_comm_compare, mp_comm_free, mp_comm_null, mp_environ, mp_irecv, &
       mp_isend, mp_rank_compare, mp_sum, mp_sync, mp_waitall
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB

#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fft_tools'

  ! Types for the pool of scratch data needed in FFT routines
  ! keep the subroutine "is_equal" up-to-date
  ! needs a default initialization
  TYPE fft_scratch_sizes
     INTEGER                              :: nx=0, ny=0, nz=0
     INTEGER                              :: lmax=0, mmax=0, nmax=0
     INTEGER                              :: mx1=0, mx2=0, mx3=0
     INTEGER                              :: my1=0, my2=0, my3=0
     INTEGER                              :: mz1=0, mz2=0, mz3=0
     INTEGER                              :: mcz1=0, mcz2=0, mcy3=0, mcx2=0
     INTEGER                              :: lg=0, mg=0
     INTEGER                              :: nbx=0, nbz=0
     INTEGER                              :: nmray=0, nyzray=0
     INTEGER                              :: gs_group=0, rs_group=0
     INTEGER, DIMENSION(2)                :: g_pos=0, r_pos=0, r_dim=0
     INTEGER                              :: numtask=0
  END TYPE fft_scratch_sizes

  TYPE fft_scratch_type
     INTEGER                              :: fft_scratch_id
     INTEGER                              :: tf_type
     LOGICAL                              :: in_use
     INTEGER                              :: group
     INTEGER, DIMENSION(3)                :: nfft
     ! to be used in cube_transpose_* routines
     INTEGER, DIMENSION(2)                :: cart_sub_comm=mp_comm_null, dim, pos
     ! to be used in fft3d_s
     COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER &
                                          :: ziptr, zoptr
     ! to be used in fft3d_ps : block distribution
     COMPLEX(KIND=dp), DIMENSION(:, :), POINTER &
            :: p1buf, p2buf, p3buf, p4buf, p5buf, p6buf, p7buf
     ! to be used in fft3d_ps : plane distribution
     COMPLEX(KIND=dp), DIMENSION(:, :), POINTER &
                                          :: r1buf, r2buf
     COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER &
                                          :: tbuf
     ! to be used in fft3d_pb
     COMPLEX(KIND=dp), DIMENSION(:, :), POINTER &
            :: a1buf, a2buf, a3buf, a4buf, a5buf, a6buf
     ! to be used in communication routines
     INTEGER, DIMENSION(:), POINTER       :: scount, rcount, sdispl, rdispl
     INTEGER, DIMENSION(:,:), POINTER     :: pgcube
     INTEGER, DIMENSION(:), POINTER       :: xzcount, yzcount, xzdispl, yzdispl
     INTEGER                              :: in, mip
     REAL(KIND=dp)                        :: rsratio
     COMPLEX(KIND=dp), DIMENSION(:), POINTER &
                                          :: xzbuf, yzbuf
     COMPLEX(KIND=sp), DIMENSION(:), POINTER &
                                          :: xzbuf_sgl, yzbuf_sgl
     COMPLEX(KIND=dp), DIMENSION(:, :), POINTER &
            :: rbuf1, rbuf2, rbuf3, rbuf4, rbuf5, rbuf6, rr
     COMPLEX(KIND=sp), DIMENSION(:, :), POINTER &
                                          :: ss, tt
     INTEGER, DIMENSION(:,:), POINTER     :: pgrid
     INTEGER, DIMENSION(:), POINTER       :: xcor, zcor, pzcoord
     TYPE(fft_scratch_sizes)              :: sizes
     TYPE(fft_plan_type), DIMENSION (6)   :: fft_plan
     INTEGER                              :: last_tick
  END TYPE fft_scratch_type

  TYPE fft_scratch_pool_type
     TYPE(fft_scratch_type),POINTER       :: fft_scratch
     TYPE(fft_scratch_pool_type),POINTER  :: fft_scratch_next
  END TYPE fft_scratch_pool_type

  INTEGER, SAVE                           :: init_fft_pool = 0
  ! the clock for fft pool. Allows to identify the least recently used scratch
  INTEGER, SAVE                           :: tick_fft_pool = 0
  ! limit the number of scratch pools to fft_pool_scratch_limit.
  INTEGER, SAVE                           :: fft_pool_scratch_limit=15
  TYPE(fft_scratch_pool_type),POINTER,SAVE:: fft_scratch_first
  ! END of types for the pool of scratch data needed in FFT routines

  PRIVATE
  PUBLIC :: init_fft, fft3d, finalize_fft
  PUBLIC :: fft_radix_operations
  PUBLIC :: FWFFT, BWFFT
  PUBLIC :: FFT_RADIX_CLOSEST, FFT_RADIX_NEXT
  PUBLIC :: FFT_RADIX_NEXT_ODD

  INTEGER, PARAMETER :: FWFFT = +1, BWFFT = -1
  INTEGER, PARAMETER :: FFT_RADIX_CLOSEST = 493, FFT_RADIX_NEXT = 494
  INTEGER, PARAMETER :: FFT_RADIX_ALLOWED = 495, FFT_RADIX_DISALLOWED = 496
  INTEGER, PARAMETER :: FFT_RADIX_NEXT_ODD = 497

  REAL(KIND=dp), PARAMETER :: ratio_sparse_alltoall = 0.5_dp

  ! these saved variables are FFT globals
  INTEGER, SAVE :: fft_type = 0
  LOGICAL, SAVE :: alltoall_sgl = .FALSE.
  LOGICAL, SAVE :: use_fftsg_sizes = .TRUE.
  INTEGER, SAVE :: fft_plan_style = 1

  ! these are only needed for pw_methods_cuda (-D__PW_CUDA)
  PUBLIC :: get_fft_scratch, release_fft_scratch
  PUBLIC :: cube_transpose_1, cube_transpose_2
  PUBLIC :: yz_to_x, x_to_yz, xz_to_yz, yz_to_xz
  PUBLIC :: fft_scratch_sizes, fft_scratch_type

  INTERFACE fft3d
     MODULE PROCEDURE fft3d_s, fft3d_ps, fft3d_pb
  END INTERFACE

#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
    INTERFACE
      INTEGER(KIND=C_INT) FUNCTION cudaFreeHost(buffer)  BIND(C, name="cudaFreeHost")
        USE ISO_C_BINDING
        IMPLICIT NONE
        TYPE(C_PTR), VALUE :: buffer
      END FUNCTION cudaFreeHost
    END INTERFACE

     INTERFACE
      INTEGER(KIND=C_INT) FUNCTION cudaHostAlloc(buffer, length, flag) BIND(C, name="cudaHostAlloc")
        USE ISO_C_BINDING
        IMPLICIT NONE
        TYPE(C_PTR)                   :: buffer
        INTEGER(KIND=C_SIZE_T), VALUE :: length
        INTEGER(KIND=C_INT), VALUE    :: flag
      END FUNCTION cudaHostAlloc
    END INTERFACE
#endif

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param fftlib ...
!> \param alltoall ...
!> \param fftsg_sizes ...
!> \param pool_limit ...
!> \param wisdom_file ...
!> \param plan_style ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE init_fft ( fftlib, alltoall, fftsg_sizes, pool_limit, wisdom_file,&
       plan_style, error )

    CHARACTER(LEN=*), INTENT(IN)             :: fftlib
    LOGICAL, INTENT(IN)                      :: alltoall, fftsg_sizes
    INTEGER, INTENT(IN)                      :: pool_limit
    CHARACTER(LEN=*), INTENT(IN)             :: wisdom_file
    INTEGER, INTENT(IN)                      :: plan_style
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_fft', &
      routineP = moduleN//':'//routineN

    use_fftsg_sizes = fftsg_sizes
    alltoall_sgl = alltoall
    fft_pool_scratch_limit = pool_limit
    fft_type = fft_library ( fftlib )
    fft_plan_style = plan_style

    IF ( fft_type <= 0 ) CALL stop_program(routineN,moduleN,__LINE__,&
                                           "Unknown FFT library: "//TRIM(fftlib))

    CALL fft_do_init(fft_type,wisdom_file)

    ! setup the FFT scratch pool, if one is associated, clear first
    CALL release_fft_scratch_pool(error)
    CALL init_fft_scratch_pool(error)

  END SUBROUTINE init_fft

! *****************************************************************************
!> \brief does whatever is needed to finalize the current fft setup
!> \param para_env ...
!> \param wisdom_file ...
!> \param error ...
!> \par History
!>      10.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE finalize_fft( para_env, wisdom_file, error )
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: wisdom_file
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'finalize_fft', &
      routineP = moduleN//':'//routineN

! release the FFT scratch pool

    CALL release_fft_scratch_pool(error)

    ! finalize fft libs

    CALL fft_do_cleanup(fft_type, wisdom_file, para_env%ionode)

  END SUBROUTINE finalize_fft

! *****************************************************************************
!> \brief Determine the allowed lengths of FFT's   '''
!> \param radix_in ...
!> \param radix_out ...
!> \param operation ...
!> \par History
!>      new library structure (JGH)
!> \author Ari Seitsonen
! *****************************************************************************
  SUBROUTINE fft_radix_operations ( radix_in, radix_out, operation )

    INTEGER, INTENT(IN)                      :: radix_in
    INTEGER, INTENT(OUT)                     :: radix_out
    INTEGER, INTENT(IN)                      :: operation

    CHARACTER(len=*), PARAMETER :: routineN = 'fft_radix_operations', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: fft_type_sg = 1

    INTEGER                                  :: i, iloc, ldata
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: DATA

    ldata = 1024
    ALLOCATE ( DATA(ldata) )
    DATA = -1

    ! if the user wants to use fftsg sizes go for it
    IF (use_fftsg_sizes) THEN
       CALL fft_get_lengths ( fft_type_sg, DATA, ldata )
    ELSE
       CALL fft_get_lengths ( fft_type, DATA, ldata )
    ENDIF

    iloc = 0
    DO i = 1, ldata
       IF ( DATA ( i ) == radix_in ) THEN
          iloc = i
          EXIT
       ELSE
          IF ( OPERATION == FFT_RADIX_ALLOWED ) THEN
             CYCLE
          ELSE IF ( DATA ( i ) > radix_in ) THEN
             iloc = i
             EXIT
          END IF
       END IF
    END DO

    IF ( iloc == 0 ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Index to radix array not found.")
    END IF

    IF ( OPERATION == FFT_RADIX_ALLOWED ) THEN
       IF ( DATA ( iloc ) == radix_in ) THEN
          radix_out = FFT_RADIX_ALLOWED
       ELSE
          radix_out = FFT_RADIX_DISALLOWED
       END IF

    ELSE IF ( OPERATION == FFT_RADIX_CLOSEST ) THEN
       IF ( DATA ( iloc ) == radix_in ) THEN
          radix_out = DATA ( iloc )
       ELSE
          IF ( ABS ( DATA ( iloc - 1 ) - radix_in ) <= &
               ABS ( DATA ( iloc ) - radix_in ) ) THEN
             radix_out = DATA ( iloc - 1 )
          ELSE
             radix_out = DATA ( iloc )
          END IF
       END IF

    ELSE IF ( OPERATION == FFT_RADIX_NEXT ) THEN
       radix_out = DATA ( iloc )

    ELSE IF ( OPERATION == FFT_RADIX_NEXT_ODD ) THEN
       DO i = iloc, ldata
          IF ( MOD ( DATA ( i ), 2 ) == 1 ) THEN
             radix_out = DATA ( i )
             EXIT
          END IF
       END DO
       IF ( MOD ( radix_out, 2 ) == 0 ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "No odd radix found.")
       END IF

    ELSE
       CALL stop_program (routineN,moduleN,__LINE__,&
                         "Disallowed radix operation.")
    END IF

    DEALLOCATE(DATA)

  END SUBROUTINE fft_radix_operations

! *****************************************************************************
!> \brief Calls the 3D-FFT function from the initialized library
!> \param fsign ...
!> \param n ...
!> \param zin ...
!> \param zout ...
!> \param scale ...
!> \param status ...
!> \param debug ...
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE fft3d_s ( fsign, n, zin, zout, scale, status, debug )

    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL, TARGET        :: zout
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    INTEGER, INTENT(OUT), OPTIONAL           :: status
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    CHARACTER(len=*), PARAMETER :: routineN = 'fft3d_s', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(1, 1, 1), &
      TARGET                                 :: zdum
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: zoptr
    INTEGER                                  :: handle, ld(3), lo(3), sign, &
                                                stat
    LOGICAL                                  :: fft_in_place, test
    REAL(KIND=dp)                            :: in_sum, norm, out_sum
    TYPE(cp_error_type)                      :: error
    TYPE(fft_scratch_type), POINTER          :: fft_scratch

    CALL timeset(routineN,handle)

    IF ( PRESENT ( scale ) ) THEN
       norm = scale
    ELSE
       norm = 1.0_dp
    END IF

    IF ( PRESENT ( debug ) ) THEN
       test = debug
    ELSE
       test = .FALSE.
    END IF

    IF ( PRESENT ( zout ) ) THEN
       fft_in_place = .FALSE.
    ELSE
       fft_in_place = .TRUE.
    END IF

    IF ( test ) THEN
       in_sum = SUM ( ABS ( zin ) )
    ENDIF

    ld ( 1 ) = SIZE ( zin ,1 )
    ld ( 2 ) = SIZE ( zin ,2 )
    ld ( 3 ) = SIZE ( zin ,3 )

    IF ( n(1) /= ld(1) .OR. n(2) /= ld(2) .OR. n(3) /= ld(3) ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Size and dimension (zin) have to be the same.")
    END IF

    sign = fsign
    CALL get_fft_scratch(fft_scratch,tf_type=400,n=n,error=error)

    IF ( fft_in_place ) THEN
       zoptr => zdum
       IF ( fsign == FWFFT ) THEN
          CALL fft_3d( fft_scratch%fft_plan(1), norm, zin, zoptr, stat )
       ELSE
          CALL fft_3d( fft_scratch%fft_plan(2), norm, zin, zoptr, stat )
       ENDIF
    ELSE
       IF ( fsign == FWFFT ) THEN
          CALL fft_3d( fft_scratch%fft_plan(3), norm, zin, zout, stat )
       ELSE
          CALL fft_3d( fft_scratch%fft_plan(4), norm, zin, zout, stat )
       ENDIF
    END IF

    CALL release_fft_scratch(fft_scratch,error)

    IF ( PRESENT ( zout ) ) THEN
       lo ( 1 ) = SIZE ( zout ,1 )
       lo ( 2 ) = SIZE ( zout ,2 )
       lo ( 3 ) = SIZE ( zout ,3 )
       IF ( n(1) /= lo(1) .OR. n(2) /= lo(2) .OR. n(3) /= lo(3) ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Size and dimension (zout) have to be the same.")
       END IF
    END IF

    IF ( PRESENT ( status ) ) THEN
       status = stat
    END IF

    IF ( test ) THEN
       IF ( PRESENT ( zout ) ) THEN
          out_sum = SUM ( ABS ( zout ) )
          WRITE ( *, '(A)') "  Out of place 3D FFT (local)  : fft3d_s"
          WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
          WRITE ( *, '(A,T60,3I7)') "     Input array dimensions ",ld
          WRITE ( *, '(A,T60,3I7)') "     Output array dimensions ",lo
          WRITE ( *, '(A,T61,E20.14)') "     Sum of input data ",in_sum
          WRITE ( *, '(A,T61,E20.14)') "     Sum of output data ",out_sum
       ELSE
          out_sum = SUM ( ABS ( zin ) )
          WRITE ( *, '(A)') "  In place 3D FFT (local)  : fft3d_s"
          WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
          WRITE ( *, '(A,T60,3I7)') "     Input/output array dimensions ",ld
          WRITE ( *, '(A,T61,E20.14)') "     Sum of input data ",in_sum
          WRITE ( *, '(A,T61,E20.14)') "     Sum of output data ",out_sum
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE fft3d_s

! *****************************************************************************
!> \brief ...
!> \param fsign ...
!> \param n ...
!> \param cin ...
!> \param gin ...
!> \param gs_group ...
!> \param rs_group ...
!> \param yzp ...
!> \param nyzray ...
!> \param bo ...
!> \param scale ...
!> \param status ...
!> \param debug ...
! *****************************************************************************
  SUBROUTINE fft3d_ps ( fsign, n, cin, gin, gs_group, rs_group, yzp, nyzray, &
       bo, scale, status, debug )

    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: cin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: gin
    INTEGER, INTENT(IN)                      :: gs_group, rs_group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nyzray
    INTEGER, DIMENSION(:, :, 0:, :), &
      INTENT(IN)                             :: bo
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    INTEGER, INTENT(OUT), OPTIONAL           :: status
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    CHARACTER(len=*), PARAMETER :: routineN = 'fft3d_ps', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: pbuf, qbuf, rbuf, sbuf
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: tbuf
    INTEGER :: g_pos, handle, ierr, iout, lg, lmax, mcx2, mcz1, mcz2, mg, &
      mmax, mx1, mx2, my1, mz2, n1, n2, nmax, numtask, numtask_g, numtask_r, &
      nx, ny, nz, r_dim(2), r_pos(2), rp, sign, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: p2p
    LOGICAL                                  :: test
    REAL(KIND=dp)                            :: norm, sum_data
    TYPE(cp_error_type)                      :: error
    TYPE(fft_scratch_sizes)                  :: fft_scratch_size
    TYPE(fft_scratch_type), POINTER          :: fft_scratch

    CALL timeset(routineN,handle)

    IF ( PRESENT ( debug ) ) THEN
       test = debug
    ELSE
       test = .FALSE.
    END IF

    CALL mp_environ ( numtask_g, g_pos, gs_group )
    CALL mp_environ ( numtask_r, r_dim, r_pos, rs_group )
    IF ( numtask_g /= numtask_r ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Real space and G space groups are different.")
    END IF
    numtask = numtask_r
    CALL mp_comm_compare ( rs_group, gs_group, iout )
    IF ( iout >3 ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Real space and G space groups are different.")
    END IF

    IF ( PRESENT ( scale ) ) THEN
       norm = scale
    ELSE
       norm = 1.0_dp
    END IF

    sign = fsign

    lg = SIZE ( gin ,1 )
    mg = SIZE ( gin ,2 )

    nx = SIZE ( cin ,1 )
    ny = SIZE ( cin ,2 )
    nz = SIZE ( cin ,3 )

    IF ( mg == 0 ) THEN
       mmax = 1
    ELSE
       mmax = mg
    END IF
    lmax = MAX ( lg, (nx*ny*nz)/mmax + 1 )

    ALLOCATE ( p2p ( 0 : numtask - 1 ), STAT = ierr )
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "p2p",int_size*numtask)

    CALL mp_rank_compare ( gs_group, rs_group, p2p )

    rp = p2p ( g_pos )
    mx1 = bo ( 2, 1, rp, 1 ) - bo ( 1, 1, rp, 1 ) + 1
    my1 = bo ( 2, 2, rp, 1 ) - bo ( 1, 2, rp, 1 ) + 1
    mx2 = bo ( 2, 1, rp, 2 ) - bo ( 1, 1, rp, 2 ) + 1
    mz2 = bo ( 2, 3, rp, 2 ) - bo ( 1, 3, rp, 2 ) + 1

    n1 = MAXVAL ( bo(2,1,:,1) - bo(1,1,:,1) + 1 )
    n2 = MAXVAL ( bo(2,2,:,1) - bo(1,2,:,1) + 1 )
    nmax = MAX( (2*n2)/numtask, 2 ) * mx2*mz2
    nmax = MAX( nmax, n1*MAXVAL(nyzray) )
    n1 = MAXVAL ( bo(2,1,:,2) )
    n2 = MAXVAL ( bo(2,3,:,2) )

    fft_scratch_size%nx  = nx
    fft_scratch_size%ny  = ny
    fft_scratch_size%nz  = nz
    fft_scratch_size%lmax= lmax
    fft_scratch_size%mmax= mmax
    fft_scratch_size%mx1 = mx1
    fft_scratch_size%mx2 = mx2
    fft_scratch_size%my1 = my1
    fft_scratch_size%mz2 = mz2
    fft_scratch_size%lg  = lg
    fft_scratch_size%mg  = mg
    fft_scratch_size%nbx = n1
    fft_scratch_size%nbz = n2
    mcz1 = MAXVAL ( bo(2,3,:,1) - bo(1,3,:,1) + 1 )
    mcx2 = MAXVAL ( bo(2,1,:,2) - bo(1,1,:,2) + 1 )
    mcz2 = MAXVAL ( bo(2,3,:,2) - bo(1,3,:,2) + 1 )
    fft_scratch_size%mcz1 = mcz1
    fft_scratch_size%mcx2 = mcx2
    fft_scratch_size%mcz2 = mcz2
    fft_scratch_size%nmax      = nmax
    fft_scratch_size%nmray     = MAXVAL(nyzray)
    fft_scratch_size%nyzray    = nyzray(g_pos)
    fft_scratch_size%gs_group  = gs_group
    fft_scratch_size%rs_group  = rs_group
    fft_scratch_size%g_pos     = g_pos
    fft_scratch_size%r_pos     = r_pos
    fft_scratch_size%r_dim     = r_dim
    fft_scratch_size%numtask   = numtask

    IF ( test ) THEN
       IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A)') "  Parallel 3D FFT : fft3d_ps"
          WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
          WRITE ( *, '(A,T67,2I7)') "     Array dimensions (gin) ",lg,mg
          WRITE ( *, '(A,T60,3I7)') "     Array dimensions (cin) ",nx,ny,nz
       END IF
    END IF

    IF ( r_dim ( 2 ) > 1 ) THEN

       !
       ! real space is distributed over x and y coordinate
       ! we have two stages of communication
       !

       IF ( r_dim ( 1 ) == 1 ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "This processor distribution is not supported.")
       END IF
       CALL get_fft_scratch(fft_scratch,tf_type=300,n=n,fft_sizes=fft_scratch_size,error=error)

       IF ( sign == FWFFT ) THEN
          ! cin -> gin

          IF ( test ) THEN
             sum_data = ABS ( SUM ( cin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A)') "  Two step communication algorithm "
                WRITE ( *, '(A,T60,3I7)') "     Transform Z ",n(3),mx1*my1
                WRITE ( *, '(A,T60,3I7)') "     Transform Y ",n(2),mx2*mz2
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          pbuf => fft_scratch%p1buf
          qbuf => fft_scratch%p2buf

          ! FFT along z
          CALL fft_1dm ( fft_scratch%fft_plan(1), cin, qbuf, norm , stat)

          rbuf => fft_scratch%p3buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( qbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
             END IF
          END IF

          ! Exchange data ( transpose of matrix )
          CALL cube_transpose_2 ( qbuf, rs_group, bo(:,:,:,1), bo(:,:,:,2), rbuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( rbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) T",sum_data
             END IF
          END IF

          pbuf => fft_scratch%p4buf

          ! FFT along y
          CALL fft_1dm ( fft_scratch%fft_plan(2), rbuf, pbuf, 1.0_dp , stat)

          qbuf => fft_scratch%p5buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( pbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) TS",sum_data
             END IF
          END IF

          ! Exchange data ( transpose of matrix ) and sort
          CALL xz_to_yz ( pbuf, rs_group, r_dim, g_pos, p2p, yzp, nyzray, &
               bo ( :, : , : , 2 ), qbuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( qbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) TS",sum_data
             END IF
          END IF

          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(3), qbuf, gin, 1.0_dp , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
             END IF
          END IF

       ELSE IF ( sign == BWFFT ) THEN
          ! gin -> cin

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A)') "  Two step communication algorithm "
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
                WRITE ( *, '(A,T60,3I7)') "     Transform Y ",n(2),mx2*mz2
                WRITE ( *, '(A,T60,3I7)') "     Transform Z ",n(3),mx1*my1
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          pbuf => fft_scratch%p7buf

          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(4), gin, pbuf, norm , stat)

          qbuf => fft_scratch%p4buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( pbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) TS",sum_data
             END IF
          END IF

          ! Exchange data ( transpose of matrix ) and sort
          CALL yz_to_xz ( pbuf, rs_group, r_dim, g_pos, p2p, yzp, nyzray, &
               bo ( :, : , : , 2 ), qbuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( qbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) TS",sum_data
             END IF
          END IF

          rbuf => fft_scratch%p3buf

          ! FFT along y
          CALL fft_1dm ( fft_scratch%fft_plan(5), qbuf, rbuf, 1.0_dp , stat)

          pbuf => fft_scratch%p2buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( rbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) T",sum_data
             END IF
          END IF

          ! Exchange data ( transpose of matrix )
          CALL cube_transpose_1 ( rbuf, rs_group, bo(:,:,:,2), bo(:,:,:,1), pbuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( pbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) T",sum_data
             END IF
          END IF

          qbuf => fft_scratch%p1buf

          ! FFT along z
          CALL fft_1dm ( fft_scratch%fft_plan(6), pbuf, cin, 1.0_dp , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( cin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
             END IF
          END IF

       ELSE

          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Illegal fsign parameter.")

       ENDIF


       CALL release_fft_scratch(fft_scratch,error)

    ELSE

       !
       ! real space is only distributed over x coordinate
       ! we have one stage of communication, after the transform of
       ! direction x
       !

       CALL get_fft_scratch(fft_scratch,tf_type=200,n=n,fft_sizes=fft_scratch_size,error=error)

       sbuf => fft_scratch%r1buf
       tbuf => fft_scratch%tbuf

       CALL zero_c(sbuf)
       CALL zero_c(tbuf)

       IF ( sign == FWFFT ) THEN
          ! cin -> gin

          IF ( test ) THEN
             sum_data = ABS ( SUM ( cin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A)') "     One step communication algorithm "
                WRITE ( *, '(A,T60,3I7)') "     Transform YZ ",n(2),n(3),nx
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          ! FFT along y and z
          CALL fft_1dm ( fft_scratch%fft_plan(1), cin, sbuf, 1._dp , stat)
          CALL fft_1dm ( fft_scratch%fft_plan(2), sbuf, tbuf, 1._dp , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( tbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) TS",sum_data
             END IF
          END IF

          ! Exchange data ( transpose of matrix ) and sort
          CALL yz_to_x ( tbuf, gs_group, g_pos, p2p, yzp, nyzray, &
               bo ( :, :, :, 2 ), sbuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( sbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) TS",sum_data
             END IF
          END IF
          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(3), sbuf, gin, norm , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
             END IF
          END IF

       ELSE IF ( sign == BWFFT ) THEN
          ! gin -> cin

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A)') "  One step communication algorithm "
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
                WRITE ( *, '(A,T60,3I7)') "     Transform YZ ",n(2),n(3),nx
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
            END IF
          END IF

          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(4), gin, sbuf, norm , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( sbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) TS",sum_data
             END IF
          END IF

          ! Exchange data ( transpose of matrix ) and sort
          CALL x_to_yz ( sbuf, gs_group, g_pos, p2p, yzp, nyzray, &
               bo ( :, :, :, 2 ), tbuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( tbuf ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) TS",sum_data
             END IF
          END IF

          ! FFT along y and z
          CALL fft_1dm ( fft_scratch%fft_plan(5), tbuf, sbuf, 1._dp , stat)
          CALL fft_1dm ( fft_scratch%fft_plan(6), sbuf, cin, 1._dp , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( cin ) )
             CALL mp_sum ( sum_data, gs_group )
             IF ( g_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
             END IF
          END IF
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Illegal fsign parameter.")
       ENDIF

       CALL release_fft_scratch(fft_scratch,error)

    ENDIF

    DEALLOCATE ( p2p, STAT = ierr )
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"p2p")

    IF ( PRESENT ( status ) ) THEN
       status = stat
    END IF
    CALL timestop(handle)


  END SUBROUTINE fft3d_ps

! *****************************************************************************
!> \brief ...
!> \param fsign ...
!> \param n ...
!> \param zin ...
!> \param gin ...
!> \param group ...
!> \param bo ...
!> \param scale ...
!> \param status ...
!> \param debug ...
! *****************************************************************************
  SUBROUTINE fft3d_pb ( fsign, n, zin, gin, group, bo, scale, status, debug )

    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(3), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: gin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:, :), &
      INTENT(IN)                             :: bo
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    INTEGER, INTENT(OUT), OPTIONAL           :: status
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    CHARACTER(len=*), PARAMETER :: routineN = 'fft3d_pb', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: abuf, bbuf
    INTEGER                                  :: handle, lg(2), lz(3), mcx2, &
                                                mcy3, mcz1, mcz2, mx1, mx2, &
                                                mx3, my1, my2, my3, my_pos, &
                                                mz1, mz2, mz3, np, sign, stat
    INTEGER, DIMENSION(2)                    :: dim, pos
    LOGICAL                                  :: failure, test
    REAL(KIND=dp)                            :: norm, sum_data
    TYPE(cp_error_type)                      :: error
    TYPE(fft_scratch_sizes)                  :: fft_scratch_size
    TYPE(fft_scratch_type), POINTER          :: fft_scratch

!------------------------------------------------------------------------------
! "Real Space"  1) xyZ      or      1) xYZ
!               2) xYz      or         not used
! "G Space"     3) Xyz      or      3) XYz
!
! There is one communicator (2-dimensional) for all distributions
! np = n1 * n2, where np is the total number of processors
! If n2 = 1, we have the second case and only one transpose step is needed
!
! Assignment of dimensions to axis for different steps
! First case: 1) n1=x; n2=y
!             2) n1=x; n2=z
!             3) n1=y; n2=z
! Second case 1) n1=x
!             3) n1=z
!
! The more general case with two communicators for the initial and final
! distribution is not covered.
!------------------------------------------------------------------------------

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CALL mp_environ ( np, dim, pos, group )
    CALL mp_cart_rank ( group, pos, my_pos )

    IF ( PRESENT ( debug ) ) THEN
       test = debug
    ELSE
       test = .FALSE.
    END IF

    IF ( PRESENT ( scale ) ) THEN
       norm = scale
    ELSE
       norm = 1.0_dp
    END IF

    sign = fsign

    IF ( test ) THEN
       lg ( 1 ) = SIZE ( gin ,1 )
       lg ( 2 ) = SIZE ( gin ,2 )
       lz ( 1 ) = SIZE ( zin ,1 )
       lz ( 2 ) = SIZE ( zin ,2 )
       lz ( 3 ) = SIZE ( zin ,3 )
       IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A)') "  Parallel 3D FFT : fft3d_pb"
          WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
          WRITE ( *, '(A,T67,2I7)') "     Array dimensions (gin) ",lg
          WRITE ( *, '(A,T60,3I7)') "     Array dimensions (cin) ",lz
       END IF
    END IF

    mx1 = bo ( 2, 1, my_pos, 1 ) - bo ( 1, 1, my_pos, 1 ) + 1
    my1 = bo ( 2, 2, my_pos, 1 ) - bo ( 1, 2, my_pos, 1 ) + 1
    mz1 = bo ( 2, 3, my_pos, 1 ) - bo ( 1, 3, my_pos, 1 ) + 1
    mx2 = bo ( 2, 1, my_pos, 2 ) - bo ( 1, 1, my_pos, 2 ) + 1
    my2 = bo ( 2, 2, my_pos, 2 ) - bo ( 1, 2, my_pos, 2 ) + 1
    mz2 = bo ( 2, 3, my_pos, 2 ) - bo ( 1, 3, my_pos, 2 ) + 1
    mx3 = bo ( 2, 1, my_pos, 3 ) - bo ( 1, 1, my_pos, 3 ) + 1
    my3 = bo ( 2, 2, my_pos, 3 ) - bo ( 1, 2, my_pos, 3 ) + 1
    mz3 = bo ( 2, 3, my_pos, 3 ) - bo ( 1, 3, my_pos, 3 ) + 1
    fft_scratch_size%mx1 = mx1
    fft_scratch_size%mx2 = mx2
    fft_scratch_size%mx3 = mx3
    fft_scratch_size%my1 = my1
    fft_scratch_size%my2 = my2
    fft_scratch_size%my3 = my3
    fft_scratch_size%mz1 = mz1
    fft_scratch_size%mz2 = mz2
    fft_scratch_size%mz3 = mz3
    mcz1 = MAXVAL ( bo ( 2, 3, :, 1 ) - bo ( 1, 3, :, 1 ) + 1 )
    mcx2 = MAXVAL ( bo ( 2, 1, :, 2 ) - bo ( 1, 1, :, 2 ) + 1 )
    mcz2 = MAXVAL ( bo ( 2, 3, :, 2 ) - bo ( 1, 3, :, 2 ) + 1 )
    mcy3 = MAXVAL ( bo ( 2, 2, :, 3 ) - bo ( 1, 2, :, 3 ) + 1 )
    fft_scratch_size%mcz1 = mcz1
    fft_scratch_size%mcx2 = mcx2
    fft_scratch_size%mcz2 = mcz2
    fft_scratch_size%mcy3 = mcy3
    fft_scratch_size%gs_group  = group
    fft_scratch_size%rs_group  = group
    fft_scratch_size%g_pos     = my_pos
    fft_scratch_size%numtask   = DIM(1)*DIM(2)

    IF ( DIM ( 1 ) > 1 .AND. DIM ( 2 ) > 1 ) THEN

       !
       ! First case; two stages of communication
       !

       CALL get_fft_scratch(fft_scratch,tf_type=100,n=n,fft_sizes=fft_scratch_size,error=error)

       IF ( sign == FWFFT ) THEN
          ! Stage 1 -> 3

          bbuf => fft_scratch%a2buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( zin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A)') "  Two step communication algorithm "
                WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx1*my1
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          ! FFT along z
          CALL fft_1dm ( fft_scratch%fft_plan(1), zin, bbuf, norm , stat)

          abuf => fft_scratch%a3buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( bbuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
             END IF
          END IF

          CALL cube_transpose_2 ( bbuf, group, bo(:,:,:,1), bo(:,:,:,2), abuf, fft_scratch, error )

          bbuf => fft_scratch%a4buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( abuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx2*mz2
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
             END IF
          END IF

          ! FFT along y
          CALL fft_1dm ( fft_scratch%fft_plan(2), abuf, bbuf, 1.0_dp , stat)

          abuf => fft_scratch%a5buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( bbuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) T",sum_data
             END IF
          END IF

          CALL cube_transpose_4 ( bbuf, group, bo(:,:,:,2), bo(:,:,:,3), abuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( abuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my3*mz3
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) ",sum_data
             END IF
          END IF


          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(3), abuf, gin, 1.0_dp , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
             END IF
          END IF

       ELSEIF ( sign == BWFFT ) THEN
          ! Stage 3 -> 1

          bbuf => fft_scratch%a5buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A)') "  Two step communication algorithm "
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my3*mz3
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(4), gin, bbuf, 1.0_dp , stat)

          abuf => fft_scratch%a4buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( bbuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
             END IF
          END IF

          CALL cube_transpose_3 ( bbuf, group, bo(:,:,:,3), bo(:,:,:,2), abuf, fft_scratch, error )

          bbuf => fft_scratch%a3buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( abuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx2*mz2
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
             END IF
          END IF

          ! FFT along y
          CALL fft_1dm ( fft_scratch%fft_plan(5), abuf, bbuf, 1.0_dp , stat)

          abuf => fft_scratch%a2buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( bbuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) T",sum_data
             END IF
          END IF

          CALL cube_transpose_1 ( bbuf, group, bo(:,:,:,2), bo(:,:,:,1), abuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( abuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx1*my1
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) ",sum_data
             END IF
          END IF

          ! FFT along z
          CALL fft_1dm ( fft_scratch%fft_plan(6), abuf, zin, norm , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( zin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
             END IF
          END IF

       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Illegal fsign parameter.")
       END IF

       CALL release_fft_scratch(fft_scratch,error)

    ELSEIF ( DIM ( 2 ) == 1 ) THEN

       !
       ! Second case; one stage of communication
       !

       CALL get_fft_scratch(fft_scratch,tf_type=101,n=n,fft_sizes=fft_scratch_size,error=error)

       IF ( sign == FWFFT ) THEN
          ! Stage 1 -> 3

          IF ( test ) THEN
             sum_data = ABS ( SUM ( zin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A)') "  one step communication algorithm "
                WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx1*my1
                WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx1*mz1
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          abuf => fft_scratch%a3buf
          bbuf => fft_scratch%a4buf
          ! FFT along z and y
          CALL fft_1dm ( fft_scratch%fft_plan(1), zin, abuf, norm , stat)
          CALL fft_1dm ( fft_scratch%fft_plan(2), abuf, bbuf, 1.0_dp , stat)

          abuf => fft_scratch%a5buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( bbuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
             END IF
          END IF

          CALL cube_transpose_6 ( bbuf, group, bo(:,:,:,1), bo(:,:,:,3), abuf, fft_scratch, error )

          IF ( test ) THEN
             sum_data = ABS ( SUM ( abuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my3*mz3
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
             END IF
          END IF

          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(3), abuf, gin, 1.0_dp , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
             END IF
          END IF

       ELSEIF ( sign == BWFFT ) THEN
          ! Stage 3 -> 1

          IF ( test ) THEN
             sum_data = ABS ( SUM ( gin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A)') "  one step communication algorithm "
                WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my3*mz3
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
             END IF
          END IF

          bbuf => fft_scratch%a5buf

          ! FFT along x
          CALL fft_1dm ( fft_scratch%fft_plan(4), gin, bbuf, 1.0_dp , stat)

          abuf => fft_scratch%a4buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( bbuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
             END IF
          END IF

          CALL cube_transpose_5 ( bbuf, group, bo(:,:,:,3), bo(:,:,:,1), abuf, fft_scratch, error )

          bbuf => fft_scratch%a3buf

          IF ( test ) THEN
             sum_data = ABS ( SUM ( abuf ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx1*mz1
                WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx1*my1
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
             END IF
          END IF

          ! FFT along y
          CALL fft_1dm ( fft_scratch%fft_plan(5), abuf, bbuf, 1.0_dp , stat)

          ! FFT along z
          CALL fft_1dm ( fft_scratch%fft_plan(6), bbuf, zin, norm , stat)

          IF ( test ) THEN
             sum_data = ABS ( SUM ( zin ) )
             CALL mp_sum ( sum_data, group )
             IF ( my_pos == 0 ) THEN
                WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
             END IF
          END IF

       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Illegal fsign parameter.")
       ENDIF

       CALL release_fft_scratch(fft_scratch,error)

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,&
                         "This partition not implemented.")

    END IF

    IF ( PRESENT ( status ) ) THEN
       status = stat
    END IF

    CALL timestop(handle)

  END SUBROUTINE fft3d_pb

! *****************************************************************************
!> \brief ...
!> \param sb ...
!> \param group ...
!> \param my_pos ...
!> \param p2p ...
!> \param yzp ...
!> \param nray ...
!> \param bo ...
!> \param tb ...
!> \param fft_scratch ...
!> \param error ...
!> \par History
!>      15. Feb. 2006 : single precision all_to_all
!> \author JGH (14-Jan-2001)
! *****************************************************************************
  SUBROUTINE x_to_yz ( sb, group, my_pos, p2p, yzp, nray, bo, tb, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, INTENT(IN)                      :: group, my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: tb
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'x_to_yz', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rr
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: ss, tt
    INTEGER                                  :: handle, ip, ir, ix, ixx, iy, &
                                                iz, mpr, nm, np, nr, nx
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl

    CALL timeset(routineN,handle)

    np = SIZE ( p2p )
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl

    IF ( alltoall_sgl ) THEN
       ss => fft_scratch%ss
       tt => fft_scratch%tt
       ss(:,:) = CMPLX(sb(:,:),KIND=sp)
       tt(:,:) = 0._sp
    ELSE
       rr => fft_scratch%rr
    END IF

    mpr = p2p ( my_pos )
    nm = MAXVAL ( nray ( 0 : np - 1 ) )
    nr = nray ( my_pos )
!$omp parallel do default(none), &
!$omp             private(ix,nx), &
!$omp             shared(np,p2p,bo,nr,scount,sdispl)
    DO ip = 0, np - 1
       ix = p2p ( ip )
       nx = bo ( 2, 1, ix ) - bo ( 1, 1, ix ) + 1
       scount ( ip ) = nr * nx
       sdispl ( ip ) = nr * ( bo ( 1, 1, ix ) - 1 )
    END DO
!$omp end parallel do
    nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
!$omp parallel do default(none), &
!$omp             private(nr), &
!$omp             shared(np,nray,nx,rcount,rdispl,nm)
    DO ip = 0, np - 1
       nr = nray ( ip )
       rcount ( ip ) = nr * nx
       rdispl ( ip ) = nm * nx * ip
    END DO
!$omp end parallel do
    IF ( alltoall_sgl ) THEN
       CALL mp_alltoall ( ss, scount, sdispl, tt, rcount, rdispl, group )
    ELSE
       CALL mp_alltoall ( sb, scount, sdispl, rr, rcount, rdispl, group )
    END IF

    nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
!$omp parallel do default(none) collapse(2) &
!$omp             private(ixx,ir,iy,iz,ix) &
!$omp             shared(np,nray,nx,alltoall_sgl,yzp,tt,rr,tb)
    DO ip = 0, np - 1
       DO ix = 1, nx
          ixx = nray(ip) * ( ix - 1 )
          IF ( alltoall_sgl ) THEN
             DO ir = 1, nray ( ip )
                iy = yzp ( 1, ir, ip )
                iz = yzp ( 2, ir, ip )
                tb ( iy, iz, ix ) = tt ( ir + ixx, ip )
             END DO
          ELSE
             DO ir = 1, nray ( ip )
                iy = yzp ( 1, ir, ip )
                iz = yzp ( 2, ir, ip )
                tb ( iy, iz, ix ) = rr ( ir + ixx, ip )
             END DO
          END IF
       END DO
    END DO
!$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE x_to_yz

! *****************************************************************************
!> \brief ...
!> \param tb ...
!> \param group ...
!> \param my_pos ...
!> \param p2p ...
!> \param yzp ...
!> \param nray ...
!> \param bo ...
!> \param sb ...
!> \param fft_scratch ...
!> \param error ...
!> \par History
!>      15. Feb. 2006 : single precision all_to_all
!> \author JGH (14-Jan-2001)
! *****************************************************************************
  SUBROUTINE yz_to_x ( tb, group, my_pos, p2p, yzp, nray, bo, sb, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: tb
    INTEGER, INTENT(IN)                      :: group, my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: sb
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'yz_to_x', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rr
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: ss, tt
    INTEGER                                  :: handle, ip, ir, ix, ixx, iy, &
                                                iz, mpr, nm, np, nr, nx
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl

    CALL timeset(routineN,handle)

    np = SIZE ( p2p )
    mpr = p2p ( my_pos )
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl

    IF ( alltoall_sgl ) THEN
       ss => fft_scratch%ss
       tt => fft_scratch%tt
       ss=0._sp
       tt=0._sp
    ELSE
       rr => fft_scratch%rr
    END IF

    nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
!$omp parallel do default(none) collapse(2) &
!$omp             private(ip, ixx, ir, iy, iz, ix) &
!$omp             shared(np,nray,nx,alltoall_sgl,yzp,tb,tt,rr)
    DO ip = 0, np - 1
       DO ix = 1, nx
          ixx = nray(ip) * ( ix - 1 )
          IF ( alltoall_sgl ) THEN
             DO ir = 1, nray ( ip )
                iy = yzp ( 1, ir, ip )
                iz = yzp ( 2, ir, ip )
                tt ( ir + ixx, ip ) = CMPLX(tb ( iy, iz, ix ),KIND=sp)
             END DO
          ELSE
             DO ir = 1, nray ( ip )
                iy = yzp ( 1, ir, ip )
                iz = yzp ( 2, ir, ip )
                rr ( ir + ixx, ip ) = tb ( iy, iz, ix )
             END DO
          END IF
       END DO
    END DO
!$omp end parallel do
    nm = MAXVAL ( nray ( 0 : np - 1 ) )
    nr = nray ( my_pos )
!$omp parallel do default(none), &
!$omp             private(ix,nx), &
!$omp             shared(np,p2p,bo,rcount,rdispl,nr)
    DO ip = 0, np - 1
       ix = p2p ( ip )
       nx = bo ( 2, 1, ix ) - bo ( 1, 1, ix ) + 1
       rcount ( ip ) = nr * nx
       rdispl ( ip ) = nr * ( bo ( 1, 1, ix ) - 1 )
    END DO
!$omp end parallel do
    nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
!$omp parallel do default(none), &
!$omp             private(nr), &
!$omp             shared(np,nray,scount,sdispl,nx,nm)
    DO ip = 0, np - 1
       nr = nray ( ip )
       scount ( ip ) = nr * nx
       sdispl ( ip ) = nm * nx * ip
    END DO
!$omp end parallel do

    IF ( alltoall_sgl ) THEN
       CALL mp_alltoall ( tt, scount, sdispl, ss, rcount, rdispl, group )
       sb = ss
    ELSE
       CALL mp_alltoall ( rr, scount, sdispl, sb, rcount, rdispl, group )
    END IF

    CALL timestop(handle)

  END SUBROUTINE yz_to_x

! *****************************************************************************
!> \brief ...
!> \param sb ...
!> \param group ...
!> \param dims ...
!> \param my_pos ...
!> \param p2p ...
!> \param yzp ...
!> \param nray ...
!> \param bo ...
!> \param tb ...
!> \param fft_scratch ...
!> \param error ...
!> \par History
!>      15. Feb. 2006 : single precision all_to_all
!> \author JGH (18-Jan-2001)
! *****************************************************************************
  SUBROUTINE yz_to_xz ( sb, group, dims, my_pos, p2p, yzp, nray, bo, tb, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(2), INTENT(IN)        :: dims
    INTEGER, INTENT(IN)                      :: my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: tb
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'yz_to_xz', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: xzbuf, yzbuf
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: xzbuf_sgl, yzbuf_sgl
    INTEGER                                  :: handle, icrs, ip, ipl, ipr, &
                                                ir, ix, iz, jj, jx, jy, jz, &
                                                myx, myz, np, npx, npz, nx, &
                                                nz, rs_pos
    INTEGER, DIMENSION(:), POINTER           :: pzcoord, rcount, rdispl, &
                                                scount, sdispl, xcor, zcor
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    CALL timeset(routineN,handle)

    np = SIZE ( p2p )

    rs_pos = p2p ( my_pos )

    IF ( alltoall_sgl ) THEN
       yzbuf_sgl => fft_scratch%yzbuf_sgl
       xzbuf_sgl => fft_scratch%xzbuf_sgl
    ELSE
       yzbuf => fft_scratch%yzbuf
       xzbuf => fft_scratch%xzbuf
    END IF
    npx = dims ( 1 )
    npz = dims ( 2 )
    pgrid  => fft_scratch%pgrid
    xcor   => fft_scratch%xcor
    zcor   => fft_scratch%zcor
    pzcoord=> fft_scratch%pzcoord
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl

    nx = SIZE ( sb ,2 )

! If the send and recv counts are not already cached, then
! calculate and store them
    IF ( fft_scratch%in == 0 ) THEN

      scount = 0

      DO ix = 0, npx - 1
         ip = pgrid ( ix, 0 )
         xcor ( bo ( 1, 1, ip ) : bo ( 2, 1, ip ) ) = ix
      END DO
      DO iz = 0, npz - 1
         ip = pgrid ( 0, iz )
         zcor ( bo ( 1, 3, ip ) : bo ( 2, 3, ip ) ) = iz
      END DO
      DO jx = 1, nx
         IF ( alltoall_sgl ) THEN
            DO ir = 1, nray ( my_pos )
               jy = yzp ( 1, ir, my_pos )
               jz = yzp ( 2, ir, my_pos )
               ip = pgrid ( xcor ( jx ), zcor ( jz ) )
               scount ( ip ) = scount ( ip ) + 1
            END DO
         ELSE
            DO ir = 1, nray ( my_pos )
               jy = yzp ( 1, ir, my_pos )
               jz = yzp ( 2, ir, my_pos )
               ip = pgrid ( xcor ( jx ), zcor ( jz ) )
               scount ( ip ) = scount ( ip ) + 1
            END DO
         END IF
      END DO

      CALL mp_alltoall ( scount, rcount, 1, group )
      fft_scratch%yzcount = scount
      fft_scratch%xzcount = rcount

      ! Work out the correct displacements in the buffers
      sdispl(0) = 0
      rdispl(0) = 0
      DO ip = 1, np - 1
         sdispl ( ip ) = sdispl (ip-1) + scount (ip-1)
         rdispl ( ip ) = rdispl (ip-1) + rcount (ip-1)
      END DO

      fft_scratch%yzdispl = sdispl
      fft_scratch%xzdispl = rdispl

      icrs = 0
      DO ip = 0, np - 1
         IF ( scount(ip) /= 0 ) icrs = icrs +1
         IF ( rcount(ip) /= 0 ) icrs = icrs +1
      END DO
      CALL mp_sum(icrs,group)
      fft_scratch%rsratio = REAL(icrs,KIND=dp)/(REAL(2*np,KIND=dp)*REAL(np,KIND=dp))

      fft_scratch%in = 1
    ELSE
      scount = fft_scratch%yzcount
      rcount = fft_scratch%xzcount
      sdispl = fft_scratch%yzdispl
      rdispl = fft_scratch%xzdispl
    END IF

! Do the actual packing
!$omp parallel do default(none), &
!$omp             private(ipl,jj,nx,ir,jx,jy,jz),&
!$omp             shared(np,p2p,pzcoord,bo,nray,yzp,zcor),&
!$omp             shared(yzbuf,sb,scount,sdispl,my_pos),&
!$omp             shared(yzbuf_sgl,alltoall_sgl)
    DO ip = 0, np - 1
       IF (scount(ip) == 0) CYCLE
       ipl = p2p(ip)
       jj = 0
       nx = bo (2, 1, ipl) - bo (1, 1, ipl) + 1
       DO ir = 1, nray(my_pos)
         jz = yzp( 2, ir, my_pos )
         IF ( zcor ( jz ) == pzcoord(ipl) ) THEN
           jj = jj + 1
           jy = yzp (1, ir, my_pos )
           IF( alltoall_sgl ) THEN
             DO jx = 0, nx - 1
              yzbuf_sgl ( sdispl (ip) + jj + jx * scount(ip) / nx ) = CMPLX(sb( ir, jx + bo (1, 1, ipl) ),KIND=sp)
             END DO
           ELSE
             DO jx = 0, nx - 1
              yzbuf ( sdispl (ip) + jj + jx * scount(ip) / nx ) = sb( ir, jx + bo (1, 1, ipl) )
             END DO
           END IF
         END IF
       END DO
    END DO
!$omp end parallel do

    IF ( alltoall_sgl ) THEN
       CALL mp_alltoall ( yzbuf_sgl, scount, sdispl, xzbuf_sgl, rcount, rdispl, group )
    ELSE
       IF ( fft_scratch%rsratio < ratio_sparse_alltoall  ) THEN
          CALL sparse_alltoall ( yzbuf, scount, sdispl, xzbuf, rcount, rdispl, group )
       ELSE
          CALL mp_alltoall ( yzbuf, scount, sdispl, xzbuf, rcount, rdispl, group )
       END IF
    END IF

    myx = fft_scratch%sizes%r_pos ( 1 )
    myz = fft_scratch%sizes%r_pos ( 2 )
    nz = bo ( 2, 3, rs_pos ) - bo ( 1, 3, rs_pos ) + 1

!$omp parallel do default(none), &
!$omp             private(ipr,jj,ir,jx,jy,jz),&
!$omp             shared(tb,np,p2p,bo,rs_pos,nray),&
!$omp             shared(yzp,alltoall_sgl,zcor,myz),&
!$omp             shared(xzbuf,xzbuf_sgl,nz,rdispl)
    DO ip = 0, np - 1
       ipr = p2p ( ip )
       jj = 0
       DO jx = 0, bo ( 2, 1, rs_pos ) - bo ( 1, 1, rs_pos )
          DO ir = 1, nray ( ip )
             jz = yzp ( 2, ir, ip )
             IF ( alltoall_sgl ) THEN
                IF ( zcor ( jz ) == myz ) THEN
                   jj = jj + 1
                   jy = yzp ( 1, ir, ip )
                   jz = jz - bo ( 1, 3, rs_pos ) + 1
                   tb ( jy, jz + jx * nz ) = xzbuf_sgl ( jj + rdispl(ipr) )
                END IF
             ELSE
                IF ( zcor ( jz ) == myz ) THEN
                   jj = jj + 1
                   jy = yzp ( 1, ir, ip )
                   jz = jz - bo ( 1, 3, rs_pos ) + 1
                   tb ( jy, jz + jx * nz ) = xzbuf ( jj + rdispl(ipr) )
                END IF
             END IF
          END DO
       END DO
    END DO
!$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE yz_to_xz

! *****************************************************************************
!> \brief ...
!> \param sb ...
!> \param group ...
!> \param dims ...
!> \param my_pos ...
!> \param p2p ...
!> \param yzp ...
!> \param nray ...
!> \param bo ...
!> \param tb ...
!> \param fft_scratch ...
!> \param error ...
!> \par History
!>      15. Feb. 2006 : single precision all_to_all
!> \author JGH (19-Jan-2001)
! *****************************************************************************
  SUBROUTINE xz_to_yz ( sb, group, dims, my_pos, p2p, yzp, nray, bo, tb, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(2), INTENT(IN)        :: dims
    INTEGER, INTENT(IN)                      :: my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: tb
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xz_to_yz', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: xzbuf, yzbuf
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: xzbuf_sgl, yzbuf_sgl
    INTEGER                                  :: handle, icrs, ip, ipl, ir, &
                                                ix, ixx, iz, jj, jx, jy, jz, &
                                                mp, myx, myz, np, npx, npz, &
                                                nx, nz
    INTEGER, DIMENSION(:), POINTER           :: pzcoord, rcount, rdispl, &
                                                scount, sdispl, xcor, zcor
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    CALL timeset(routineN,handle)

    np = SIZE ( p2p )

    IF ( alltoall_sgl ) THEN
       yzbuf_sgl => fft_scratch%yzbuf_sgl
       xzbuf_sgl => fft_scratch%xzbuf_sgl
    ELSE
       yzbuf => fft_scratch%yzbuf
       xzbuf => fft_scratch%xzbuf
    END IF
    npx = dims ( 1 )
    npz = dims ( 2 )
    pgrid  => fft_scratch%pgrid
    xcor   => fft_scratch%xcor
    zcor   => fft_scratch%zcor
    pzcoord=> fft_scratch%pzcoord
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl

! If the send and recv counts are not already cached, then
! calculate and store them
    IF ( fft_scratch%in == 0 ) THEN

      rcount = 0
      nx = MAXVAL ( bo ( 2, 1, : ) )

      DO ix = 0, npx - 1
         ip = pgrid ( ix, 0 )
         xcor ( bo ( 1, 1, ip ) : bo ( 2, 1, ip ) ) = ix
      END DO
      DO iz = 0, npz - 1
         ip = pgrid ( 0, iz )
         zcor ( bo ( 1, 3, ip ) : bo ( 2, 3, ip ) ) = iz
      END DO
      DO jx = 1, nx
         DO ir = 1, nray ( my_pos )
            jy = yzp ( 1, ir, my_pos )
            jz = yzp ( 2, ir, my_pos )
            ip = pgrid ( xcor ( jx ), zcor ( jz ) )
            rcount ( ip ) = rcount ( ip ) + 1
         END DO
      END DO

      CALL mp_alltoall ( rcount, scount, 1, group )
      fft_scratch%xzcount = scount
      fft_scratch%yzcount = rcount

      ! Work out the correct displacements in the buffers
      sdispl(0) = 0
      rdispl(0) = 0
      DO ip = 1, np - 1
         sdispl ( ip ) = sdispl (ip-1) + scount (ip-1)
         rdispl ( ip ) = rdispl (ip-1) + rcount (ip-1)
      END DO

      fft_scratch%xzdispl = sdispl
      fft_scratch%yzdispl = rdispl

       icrs = 0
       DO ip = 0, np - 1
          IF ( scount(ip) /= 0 ) icrs = icrs +1
          IF ( rcount(ip) /= 0 ) icrs = icrs +1
       END DO
       CALL mp_sum(icrs,group)
       fft_scratch%rsratio = REAL(icrs,KIND=dp)/(REAL(2*np,KIND=dp)*REAL(np,KIND=dp))

       fft_scratch%in = 1
    ELSE
       scount = fft_scratch%xzcount
       rcount = fft_scratch%yzcount
       sdispl = fft_scratch%xzdispl
       rdispl = fft_scratch%yzdispl
    END IF


! Now do the actual packing
    myx = fft_scratch%sizes%r_pos ( 1 )
    myz = fft_scratch%sizes%r_pos ( 2 )
    mp = p2p ( my_pos )
    nz = bo ( 2, 3, mp ) - bo ( 1, 3, mp ) + 1
    nx = bo ( 2, 1, mp ) - bo ( 1, 1, mp ) + 1

!$omp parallel do default(none), &
!$omp             private(jj,ipl,ir,jx,jy,jz,ixx),&
!$omp             shared(np,p2p,nray,yzp,zcor,myz,bo,mp),&
!$omp             shared(alltoall_sgl,nx,scount,sdispl),&
!$omp             shared(xzbuf,xzbuf_sgl,sb,nz)
    DO ip = 0, np - 1
       jj = 0
       ipl = p2p ( ip )
       DO ir = 1, nray ( ip )
          jz = yzp ( 2, ir, ip )
          IF ( zcor ( jz ) == myz ) THEN
             jj = jj + 1
             jy = yzp ( 1, ir, ip )
             jz = yzp ( 2, ir, ip ) - bo ( 1, 3, mp ) + 1
             IF ( alltoall_sgl ) THEN
                DO jx = 0, nx - 1
                   ixx = jj + jx * scount ( ipl )/nx
                   xzbuf_sgl ( ixx + sdispl(ipl) ) = CMPLX(sb ( jy, jz + jx * nz ),KIND=sp)
                END DO
             ELSE
                DO jx = 0, nx - 1
                   ixx = jj + jx * scount ( ipl )/nx
                   xzbuf ( ixx + sdispl(ipl) ) = sb ( jy, jz + jx * nz )
                END DO
             END IF
          END IF
       END DO
    END DO
!$omp end parallel do

    IF ( alltoall_sgl ) THEN
       CALL mp_alltoall ( xzbuf_sgl, scount, sdispl, yzbuf_sgl, rcount, rdispl, group )
    ELSE
       IF ( fft_scratch%rsratio < ratio_sparse_alltoall  ) THEN
          CALL sparse_alltoall ( xzbuf, scount, sdispl, yzbuf, rcount, rdispl, group )
       ELSE
          CALL mp_alltoall ( xzbuf, scount, sdispl, yzbuf, rcount, rdispl, group )
       END IF
    END IF

!$omp parallel do default(none), &
!$omp             private(ipl,jj,nx,ir,jx,jy,jz),&
!$omp             shared(p2p,pzcoord,bo,nray,my_pos,yzp),&
!$omp             shared(rcount,rdispl,tb,yzbuf,zcor),&
!$omp             shared(yzbuf_sgl,alltoall_sgl,np)
    DO ip = 0, np - 1
       IF (rcount(ip) == 0) CYCLE
       ipl = p2p(ip)
       jj = 0
       nx = bo (2, 1, ipl) - bo (1, 1, ipl) + 1
       DO ir = 1, nray(my_pos)
         jz = yzp( 2, ir, my_pos )
         IF ( zcor ( jz ) == pzcoord(ipl)) THEN
           jj = jj + 1
           jy = yzp (1, ir, my_pos )
           IF ( alltoall_sgl ) THEN
             DO jx = 0, nx - 1
               tb( ir, jx + bo (1, 1, ipl) ) = yzbuf_sgl ( rdispl (ip) + jj + jx * rcount(ip) / nx )
             END DO
           ELSE
             DO jx = 0, nx - 1
               tb( ir, jx + bo (1, 1, ipl) ) = yzbuf ( rdispl (ip) + jj + jx * rcount(ip) / nx )
             END DO
           END IF
         END IF
       END DO
    END DO
!$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE xz_to_yz

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param group ...
!> \param boin ...
!> \param boout ...
!> \param sout ...
!> \param fft_scratch ...
!> \param error ...
!> \par History
!>      none
!> \author JGH (20-Jan-2001)
! *****************************************************************************
  SUBROUTINE cube_transpose_1 ( cin, group, boin, boout, sout, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: cin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_transpose_1', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rbuf
    INTEGER                                  :: handle, ip, ipl, ir, is, ixy, &
                                                iz, mip, mz, np, nx, ny, nz, &
                                                sub_group
    INTEGER, DIMENSION(2)                    :: dim, pos
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    CALL timeset(routineN,handle)

    sub_group = fft_scratch%cart_sub_comm(2)
    mip    =  fft_scratch%mip
    dim    =  fft_scratch%dim
    pos    =  fft_scratch%pos
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl
    pgrid  => fft_scratch%pgcube
    np = DIM ( 2 )

    nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
    nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1

!$omp parallel do default(none), &
!$omp             private(ipl,ny), &
!$omp             shared(np,pgrid,boout,scount,sdispl,nx,nz)
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 2 )
       ny = boout ( 2, 2, ipl ) - boout ( 1, 2, ipl ) + 1
       scount ( ip ) = nx * nz * ny
       sdispl ( ip ) = nx * nz * ( boout ( 1, 2, ipl ) - 1 )
    END DO
!$omp end parallel do
    ny = boout ( 2, 2, mip ) - boout ( 1, 2, mip ) + 1
    mz = MAXVAL ( boin ( 2, 3, : ) - boin ( 1, 3, : ) + 1 )
!$omp parallel do default(none), &
!$omp             private(ipl,nz), &
!$omp             shared(np,pgrid,boin,nx,ny,rcount,rdispl,mz)
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 2 )
       nz = boin ( 2, 3, ipl ) - boin ( 1, 3, ipl ) + 1
       rcount ( ip ) = nx * nz * ny
       rdispl ( ip ) = nx * ny * mz * ip
    END DO
!$omp end parallel do

    rbuf => fft_scratch%rbuf1

    CALL mp_alltoall ( cin, scount, sdispl, rbuf, rcount, rdispl, sub_group )

!$omp parallel do default(none) collapse(2) &
!$omp             private(ip,ipl,nz,iz,is,ir) &
!$omp             shared(nx,ny,np,pgrid,boin,sout,rbuf)
    DO ixy = 1, nx * ny
       DO ip = 0, np - 1
          ipl = pgrid ( ip, 2 )
          nz = boin ( 2, 3, ipl ) - boin ( 1, 3, ipl ) + 1
          DO iz = 1, nz
             is = boin ( 1, 3, ipl ) + iz - 1
             ir = iz + nz * ( ixy - 1 )
             sout ( is, ixy ) = rbuf ( ir, ip )
          END DO
       END DO
    END DO
!$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE cube_transpose_1

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param group ...
!> \param boin ...
!> \param boout ...
!> \param sout ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cube_transpose_2 ( cin, group, boin, boout, sout, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: cin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_transpose_2', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rbuf
    INTEGER                                  :: handle, ip, ipl, ir, ixy, iz, &
                                                mip, mz, np, nx, ny, nz, &
                                                sub_group
    INTEGER, DIMENSION(2)                    :: dim, pos
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    CALL timeset(routineN,handle)

    sub_group = fft_scratch%cart_sub_comm(2)
    mip    =  fft_scratch%mip
    dim    =  fft_scratch%dim
    pos    =  fft_scratch%pos
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl
    pgrid  => fft_scratch%pgcube
    np = DIM ( 2 )

    nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
    ny = boin ( 2, 2, mip ) - boin ( 1, 2, mip ) + 1
    mz = MAXVAL ( boout ( 2, 3, : ) - boout ( 1, 3, : ) + 1 )

    rbuf => fft_scratch%rbuf2

!$omp parallel default(none), &
!$omp          private(ip,ipl,nz,iz,ir), &
!$omp          shared(nx,ny,np,pgrid,boout,rbuf,cin,scount,sdispl,mz)
!$omp do collapse(2)
    DO ixy = 1, nx * ny
       DO ip = 0, np - 1
          ipl = pgrid ( ip, 2 )
          nz = boout ( 2, 3, ipl ) - boout ( 1, 3, ipl ) + 1
          DO iz = boout ( 1, 3, ipl ), boout ( 2, 3, ipl )
             ir = iz - boout ( 1, 3, ipl ) + 1 + ( ixy - 1 ) * nz
             rbuf ( ir, ip ) = cin ( iz, ixy )
          END DO
       END DO
    END DO
!$omp end do
!$omp do
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 2 )
       nz = boout ( 2, 3, ipl ) - boout ( 1, 3, ipl ) + 1
       scount ( ip ) = nx * ny * nz
       sdispl ( ip ) = nx * ny * mz * ip
    END DO
!$omp end do
!$omp end parallel
    nz = boout ( 2, 3, mip ) - boout ( 1, 3, mip ) + 1
!$omp parallel do default(none), &
!$omp             private(ipl,ny), &
!$omp             shared(np,pgrid,boin,nx,nz,rcount,rdispl)
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 2 )
       ny = boin ( 2, 2, ipl ) - boin ( 1, 2, ipl ) + 1
       rcount ( ip ) = nx * ny * nz
       rdispl ( ip ) = nx * nz * ( boin ( 1, 2, ipl ) - 1 )
    END DO
!$omp end parallel do

    CALL mp_alltoall ( rbuf, scount, sdispl, sout, rcount, rdispl, sub_group )

    CALL timestop(handle)

  END SUBROUTINE cube_transpose_2

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param group ...
!> \param boin ...
!> \param boout ...
!> \param sout ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cube_transpose_3 ( cin, group, boin, boout, sout, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: cin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_transpose_3', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rbuf
    INTEGER                                  :: handle, ip, ipl, ir, is, ixz, &
                                                iy, lb, mip, my, my_id, np, &
                                                num_threads, nx, ny, nz, &
                                                sub_group, ub
    INTEGER, DIMENSION(2)                    :: dim, pos
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    CALL timeset(routineN,handle)

    sub_group = fft_scratch%cart_sub_comm(1)
    mip    =  fft_scratch%mip
    dim    =  fft_scratch%dim
    pos    =  fft_scratch%pos
    np = DIM ( 1 )
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl
    pgrid  => fft_scratch%pgcube

    ny = boin ( 2, 2, mip ) - boin ( 1, 2, mip ) + 1
    nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
!$omp parallel do default(none), &
!$omp             private(ipl, nx), &
!$omp             shared(np,pgrid,boout,ny,nz,scount,sdispl)
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 1 )
       nx = boout ( 2, 1, ipl ) - boout ( 1, 1, ipl ) + 1
       scount ( ip ) = nx * nz * ny
       sdispl ( ip ) = ny * nz * ( boout ( 1, 1, ipl ) - 1 )
    END DO
!$omp end parallel do
    nx = boout ( 2, 1, mip ) - boout ( 1, 1, mip ) + 1
    my = MAXVAL ( boin ( 2, 2, : ) - boin ( 1, 2, : ) + 1 )
!$omp parallel do default(none), &
!$omp             private(ipl, ny), &
!$omp             shared(np,pgrid,boin,nx,nz,my,rcount,rdispl)
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 1 )
       ny = boin ( 2, 2, ipl ) - boin ( 1, 2, ipl ) + 1
       rcount ( ip ) = nx * nz * ny
       rdispl ( ip ) = nx * my * nz * ip
    END DO
!$omp end parallel do

    rbuf => fft_scratch%rbuf3
    num_threads = 1
    my_id = 0
!$omp parallel default(none), &
!$omp          private(num_threads, my_id, lb, ub), &
!$omp          shared(rbuf)
!$  num_threads = MIN(omp_get_max_threads(), SIZE(rbuf,2))
!$  my_id = omp_get_thread_num()
    IF (my_id < num_threads) THEN
      lb = (SIZE(rbuf,2)*my_id)/num_threads
      ub = (SIZE(rbuf,2)*(my_id+1))/num_threads - 1
      rbuf(:,lb:ub) = 0.0_dp
    END IF
!$omp end parallel

    CALL mp_alltoall ( cin, scount, sdispl, rbuf, rcount, rdispl, sub_group )

!$omp parallel do default(none) collapse(2) &
!$omp             private(ip,ipl,ny,iy,is,ir) &
!$omp             shared(nx,nz,np,pgrid,boin,rbuf,sout)
    DO ixz = 1, nx * nz
       DO ip = 0, np - 1
          ipl = pgrid ( ip, 1 )
          ny = boin ( 2, 2, ipl ) - boin ( 1, 2, ipl ) + 1
          DO iy = 1, ny
             is = boin ( 1, 2, ipl ) + iy - 1
             ir = iy + ny * ( ixz - 1 )
             sout ( is, ixz ) = rbuf ( ir, ip )
          END DO
       END DO
    END DO
!$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE cube_transpose_3

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param group ...
!> \param boin ...
!> \param boout ...
!> \param sout ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cube_transpose_4 ( cin, group, boin, boout, sout, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: cin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_transpose_4', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rbuf
    INTEGER                                  :: handle, ip, ipl, ir, iy, izx, &
                                                lb, mip, my, my_id, np, &
                                                num_threads, nx, ny, nz, &
                                                sub_group, ub
    INTEGER, DIMENSION(2)                    :: dim, pos
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    CALL timeset(routineN,handle)

    sub_group = fft_scratch%cart_sub_comm(1)
    mip    =  fft_scratch%mip
    dim    =  fft_scratch%dim
    pos    =  fft_scratch%pos
    np = DIM ( 1 )
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl
    pgrid  => fft_scratch%pgcube

    nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
    nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
    my = MAXVAL ( boout ( 2, 2, : ) - boout ( 1, 2, : ) + 1 )

    rbuf => fft_scratch%rbuf4
    num_threads = 1
    my_id = 0
!$omp parallel default(none), &
!$omp          private(num_threads,my_id,lb,ub,ip,ipl,ny,iy,ir), &
!$omp          shared(rbuf,nz,nx,np,pgrid,boout,cin,my,scount,sdispl)
!$  num_threads = MIN(omp_get_max_threads(), SIZE(rbuf,2))
!$  my_id = omp_get_thread_num()
    IF (my_id < num_threads) THEN
      lb = (SIZE(rbuf,2)*my_id)/num_threads
      ub = (SIZE(rbuf,2)*(my_id+1))/num_threads - 1
      rbuf(:,lb:ub) = 0.0_dp
    END IF
!$omp barrier

!$omp do collapse(2)
    DO izx = 1, nz * nx
       DO ip = 0, np - 1
          ipl = pgrid ( ip, 1 )
          ny = boout ( 2, 2, ipl ) - boout ( 1, 2, ipl ) + 1
          DO iy = boout ( 1, 2, ipl ), boout ( 2, 2, ipl )
             ir = iy - boout ( 1, 2, ipl ) + 1 + ( izx - 1 ) * ny
             rbuf ( ir, ip ) = cin ( iy, izx )
          END DO
       END DO
    END DO
!$omp end do
!$omp do
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 1 )
       ny = boout ( 2, 2, ipl ) - boout ( 1, 2, ipl ) + 1
       scount ( ip ) = nx * ny * nz
       sdispl ( ip ) = nx * nz * my * ip
    END DO
!$omp end do
!$omp end parallel
    ny = boout ( 2, 2, mip ) - boout ( 1, 2, mip ) + 1
!$omp parallel do default(none), &
!$omp             private(ipl,nx), &
!$omp             shared(np,pgrid,boin,rcount,rdispl,ny,nz)
    DO ip = 0, np - 1
       ipl = pgrid ( ip, 1 )
       nx = boin ( 2, 1, ipl ) - boin ( 1, 1, ipl ) + 1
       rcount ( ip ) = nx * ny * nz
       rdispl ( ip ) = ny * nz * ( boin ( 1, 1, ipl ) - 1 )
    END DO
!$omp end parallel do

    CALL mp_alltoall ( rbuf, scount, sdispl, sout, rcount, rdispl, sub_group )

    CALL timestop(handle)

  END SUBROUTINE cube_transpose_4

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param group ...
!> \param boin ...
!> \param boout ...
!> \param sout ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cube_transpose_5 ( cin, group, boin, boout, sout, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: cin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_transpose_5', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rbuf
    INTEGER                                  :: handle, ip, ir, is, ixz, iy, &
                                                lb, mip, my, my_id, np, &
                                                num_threads, nx, ny, nz, ub
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl

    CALL timeset(routineN,handle)

    np     =  fft_scratch%sizes%numtask
    mip    =  fft_scratch%mip
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl

    ny = boin ( 2, 2, mip ) - boin ( 1, 2, mip ) + 1
    nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
!$omp parallel do default(none), &
!$omp             private(nx), &
!$omp             shared(np,boout,ny,nz,scount,sdispl)
    DO ip = 0, np - 1
       nx = boout ( 2, 1, ip ) - boout ( 1, 1, ip ) + 1
       scount ( ip ) = nx * nz * ny
       sdispl ( ip ) = ny * nz * ( boout ( 1, 1, ip ) - 1 )
    END DO
!$omp end parallel do
    nx = boout ( 2, 1, mip ) - boout ( 1, 1, mip ) + 1
    my = MAXVAL ( boin ( 2, 2, : ) - boin ( 1, 2, : ) + 1 )
!$omp parallel do default(none), &
!$omp             private(ny), &
!$omp             shared(np,boin,nx,nz,rcount,rdispl,my)
    DO ip = 0, np - 1
       ny = boin ( 2, 2, ip ) - boin ( 1, 2, ip ) + 1
       rcount ( ip ) = nx * nz * ny
       rdispl ( ip ) = nx * my * nz * ip
    END DO
!$omp end parallel do

    rbuf => fft_scratch%rbuf5
    num_threads = 1
    my_id = 0
!$omp parallel default(none), &
!$omp          private(num_threads, my_id, lb, ub), &
!$omp          shared(rbuf)
!$  num_threads = MIN(omp_get_max_threads(), SIZE(rbuf,2))
!$  my_id = omp_get_thread_num()
    IF (my_id < num_threads) THEN
      lb = (SIZE(rbuf,2)*my_id)/num_threads
      ub = (SIZE(rbuf,2)*(my_id+1))/num_threads - 1
      rbuf(:,lb:ub) = 0.0_dp
    END IF
!$omp end parallel


    CALL mp_alltoall ( cin, scount, sdispl, rbuf, rcount, rdispl, group )

!$omp parallel do default(none) collapse(2) &
!$omp             private(ip,ny,iy,is,ir) &
!$omp             shared(nx,nz,np,boin,sout,rbuf)
    DO ixz = 1, nx * nz
       DO ip = 0, np - 1
          ny = boin ( 2, 2, ip ) - boin ( 1, 2, ip ) + 1
          DO iy = 1, ny
             is = boin ( 1, 2, ip ) + iy - 1
             ir = iy + ny * ( ixz - 1 )
             sout ( is, ixz ) = rbuf ( ir, ip )
          END DO
       END DO
    END DO
!$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE cube_transpose_5

! *****************************************************************************
!> \brief ...
!> \param cin ...
!> \param group ...
!> \param boin ...
!> \param boout ...
!> \param sout ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cube_transpose_6 ( cin, group, boin, boout, sout, fft_scratch, error )

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: cin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cube_transpose_6', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: rbuf
    INTEGER                                  :: handle, ip, ir, iy, izx, lb, &
                                                mip, my, my_id, np, &
                                                num_threads, nx, ny, nz, ub
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl, scount, sdispl

    CALL timeset(routineN,handle)

    np     =  fft_scratch%sizes%numtask
    mip    =  fft_scratch%mip
    scount => fft_scratch%scount
    rcount => fft_scratch%rcount
    sdispl => fft_scratch%sdispl
    rdispl => fft_scratch%rdispl

    nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
    nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
    my = MAXVAL ( boout ( 2, 2, : ) - boout ( 1, 2, : ) + 1 )

    rbuf => fft_scratch%rbuf5
    num_threads = 1
    my_id = 0
!$omp parallel default(none), &
!$omp          private(num_threads,my_id,lb,ub,ip,ny,iy,ir), &
!$omp          shared(rbuf,nx,nz,np,boout,cin,my,scount,sdispl)
!$  num_threads = MIN(omp_get_max_threads(), SIZE(rbuf,2))
!$  my_id = omp_get_thread_num()
    IF (my_id < num_threads) THEN
      lb = (SIZE(rbuf,2)*my_id)/num_threads
      ub = (SIZE(rbuf,2)*(my_id+1))/num_threads - 1
      rbuf(:,lb:ub) = 0.0_dp
    END IF
!$omp barrier

!$omp do collapse(2)
    DO izx = 1, nz * nx
       DO ip = 0, np - 1
          ny = boout ( 2, 2, ip ) - boout ( 1, 2, ip ) + 1
          DO iy = boout ( 1, 2, ip ), boout ( 2, 2, ip )
             ir = iy - boout ( 1, 2, ip ) + 1 + ( izx - 1 ) * ny
             rbuf ( ir, ip ) = cin ( iy, izx )
          END DO
       END DO
    END DO
!$omp end do
!$omp do
    DO ip = 0, np - 1
       ny = boout ( 2, 2, ip ) - boout ( 1, 2, ip ) + 1
       scount ( ip ) = nx * ny * nz
       sdispl ( ip ) = nx * nz * my * ip
    END DO
!$omp end do
!$omp end parallel
    ny = boout ( 2, 2, mip ) - boout ( 1, 2, mip ) + 1
!$omp parallel do default(none), &
!$omp             private(nx), &
!$omp             shared(np,boin,rcount,rdispl,nz,ny)
    DO ip = 0, np - 1
       nx = boin ( 2, 1, ip ) - boin ( 1, 1, ip ) + 1
       rcount ( ip ) = nx * ny * nz
       rdispl ( ip ) = ny * nz * ( boin ( 1, 1, ip ) - 1 )
    END DO
!$omp end parallel do

    CALL mp_alltoall ( rbuf, scount, sdispl, sout, rcount, rdispl, group )

    CALL timestop(handle)

  END SUBROUTINE cube_transpose_6

! *****************************************************************************
!> \brief ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_fft_scratch_pool(error)

    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_fft_scratch_pool', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL release_fft_scratch_pool(error)

    ! Allocate first scratch and mark it as used
    ALLOCATE (fft_scratch_first,STAT=ierr)
    CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (fft_scratch_first%fft_scratch,STAT=ierr)
    CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    NULLIFY(fft_scratch_first%fft_scratch_next)
    fft_scratch_first%fft_scratch%fft_scratch_id=0
    fft_scratch_first%fft_scratch%in_use=.TRUE.
    fft_scratch_first%fft_scratch%group=0
    NULLIFY(fft_scratch_first%fft_scratch%ziptr)
    NULLIFY(fft_scratch_first%fft_scratch%zoptr)
    NULLIFY(fft_scratch_first%fft_scratch%p1buf)
    NULLIFY(fft_scratch_first%fft_scratch%p2buf)
    NULLIFY(fft_scratch_first%fft_scratch%p3buf)
    NULLIFY(fft_scratch_first%fft_scratch%p4buf)
    NULLIFY(fft_scratch_first%fft_scratch%p5buf)
    NULLIFY(fft_scratch_first%fft_scratch%p6buf)
    NULLIFY(fft_scratch_first%fft_scratch%p7buf)
    NULLIFY(fft_scratch_first%fft_scratch%r1buf)
    NULLIFY(fft_scratch_first%fft_scratch%r2buf)
    NULLIFY(fft_scratch_first%fft_scratch%tbuf)
    NULLIFY(fft_scratch_first%fft_scratch%a1buf)
    NULLIFY(fft_scratch_first%fft_scratch%a2buf)
    NULLIFY(fft_scratch_first%fft_scratch%a3buf)
    NULLIFY(fft_scratch_first%fft_scratch%a4buf)
    NULLIFY(fft_scratch_first%fft_scratch%a5buf)
    NULLIFY(fft_scratch_first%fft_scratch%a6buf)
    NULLIFY(fft_scratch_first%fft_scratch%scount,fft_scratch_first%fft_scratch%rcount,&
         fft_scratch_first%fft_scratch%sdispl,fft_scratch_first%fft_scratch%rdispl)
    NULLIFY(fft_scratch_first%fft_scratch%rr,&
         fft_scratch_first%fft_scratch%ss, fft_scratch_first%fft_scratch%tt)
    NULLIFY(fft_scratch_first%fft_scratch%xzbuf,fft_scratch_first%fft_scratch%yzbuf,&
         fft_scratch_first%fft_scratch%xzbuf_sgl,fft_scratch_first%fft_scratch%yzbuf_sgl)
    NULLIFY(fft_scratch_first%fft_scratch%pgrid)
    NULLIFY(fft_scratch_first%fft_scratch%pgcube)
    NULLIFY(fft_scratch_first%fft_scratch%xcor, fft_scratch_first%fft_scratch%zcor)
    NULLIFY(fft_scratch_first%fft_scratch%pzcoord)
    NULLIFY(fft_scratch_first%fft_scratch%xzcount, fft_scratch_first%fft_scratch%yzcount,&
         fft_scratch_first%fft_scratch%xzdispl, fft_scratch_first%fft_scratch%yzdispl)
    NULLIFY(fft_scratch_first%fft_scratch%rbuf1,fft_scratch_first%fft_scratch%rbuf2,&
         fft_scratch_first%fft_scratch%rbuf3,fft_scratch_first%fft_scratch%rbuf4)
    NULLIFY(fft_scratch_first%fft_scratch%rbuf5,fft_scratch_first%fft_scratch%rbuf6)
    fft_scratch_first%fft_scratch%in = 0
    fft_scratch_first%fft_scratch%rsratio = 1._dp
    DO i=1,6
       fft_scratch_first%fft_scratch%fft_plan(i)%valid = .FALSE.
    END DO
    ! this is a very special scratch, it seems, we always keep it 'most - recent' so we will never delete it
    fft_scratch_first%fft_scratch%last_tick=HUGE(fft_scratch_first%fft_scratch%last_tick)

    init_fft_pool = init_fft_pool + 1

  END SUBROUTINE init_fft_scratch_pool

! *****************************************************************************
!> \brief ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE deallocate_fft_scratch_type(fft_scratch,error)
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
    USE ISO_C_BINDING
#endif
    TYPE(fft_scratch_type), INTENT(INOUT)    :: fft_scratch
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_fft_scratch_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
    COMPLEX(KIND=dp), POINTER :: dummy_ptr_z
#endif

    failure = .FALSE.
    ! deallocate structures
    IF(ASSOCIATED(fft_scratch%ziptr)) THEN
       DEALLOCATE(fft_scratch%ziptr,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%zoptr)) THEN
       DEALLOCATE(fft_scratch%zoptr,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p1buf)) THEN
       DEALLOCATE(fft_scratch%p1buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p2buf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%p2buf(1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%p2buf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p3buf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%p3buf(1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%p3buf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p4buf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%p4buf(1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%p4buf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p5buf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%p5buf(1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%p5buf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p6buf)) THEN
       DEALLOCATE(fft_scratch%p6buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%p7buf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%p7buf(1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%p7buf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%r1buf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%r1buf(1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%r1buf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%r2buf)) THEN
       DEALLOCATE(fft_scratch%r2buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%tbuf)) THEN
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
       dummy_ptr_z => fft_scratch%tbuf(1,1,1)
       ierr = cudaFreeHost(c_loc(dummy_ptr_z))
#else
       DEALLOCATE(fft_scratch%tbuf,STAT=ierr)
#endif
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%a1buf)) THEN
       DEALLOCATE(fft_scratch%a1buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%a2buf)) THEN
       DEALLOCATE(fft_scratch%a2buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%a3buf)) THEN
       DEALLOCATE(fft_scratch%a3buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%a4buf)) THEN
       DEALLOCATE(fft_scratch%a4buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%a5buf)) THEN
       DEALLOCATE(fft_scratch%a5buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%a6buf)) THEN
       DEALLOCATE(fft_scratch%a6buf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%scount)) THEN
       DEALLOCATE(fft_scratch%scount,fft_scratch%rcount,&
            fft_scratch%sdispl,fft_scratch%rdispl,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%rr)) THEN
       DEALLOCATE(fft_scratch%rr,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%xzbuf)) THEN
       DEALLOCATE(fft_scratch%xzbuf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%yzbuf)) THEN
       DEALLOCATE(fft_scratch%yzbuf,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%xzbuf_sgl)) THEN
       DEALLOCATE(fft_scratch%xzbuf_sgl,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%yzbuf_sgl)) THEN
       DEALLOCATE(fft_scratch%yzbuf_sgl,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%ss)) THEN
       DEALLOCATE(fft_scratch%ss,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%tt)) THEN
       DEALLOCATE(fft_scratch%tt,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%pgrid)) THEN
       DEALLOCATE(fft_scratch%pgrid,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%pgcube)) THEN
       DEALLOCATE(fft_scratch%pgcube,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%xcor)) THEN
       DEALLOCATE(fft_scratch%xcor,fft_scratch%zcor,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%pzcoord)) THEN
       DEALLOCATE(fft_scratch%pzcoord,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%xzcount)) THEN
       DEALLOCATE(fft_scratch%xzcount,fft_scratch%yzcount,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(fft_scratch%xzdispl,fft_scratch%yzdispl,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
       fft_scratch%in=0
       fft_scratch%rsratio=1._dp
    END IF
    IF(ASSOCIATED(fft_scratch%rbuf1)) THEN
       DEALLOCATE(fft_scratch%rbuf1,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%rbuf2)) THEN
       DEALLOCATE(fft_scratch%rbuf2,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%rbuf3)) THEN
       DEALLOCATE(fft_scratch%rbuf3,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%rbuf4)) THEN
       DEALLOCATE(fft_scratch%rbuf4,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%rbuf5)) THEN
       DEALLOCATE(fft_scratch%rbuf5,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(fft_scratch%rbuf6)) THEN
       DEALLOCATE(fft_scratch%rbuf6,STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (fft_scratch%cart_sub_comm(1) .NE. mp_comm_null) THEN
       CALL mp_comm_free(fft_scratch%cart_sub_comm(1))
    END IF
    IF (fft_scratch%cart_sub_comm(2) .NE. mp_comm_null) THEN
       CALL mp_comm_free(fft_scratch%cart_sub_comm(2))
    END IF
    fft_scratch%cart_sub_comm = mp_comm_null

    CALL fft_destroy_plan(fft_scratch%fft_plan(1))
    CALL fft_destroy_plan(fft_scratch%fft_plan(2))
    CALL fft_destroy_plan(fft_scratch%fft_plan(3))
    CALL fft_destroy_plan(fft_scratch%fft_plan(4))
    CALL fft_destroy_plan(fft_scratch%fft_plan(5))
    CALL fft_destroy_plan(fft_scratch%fft_plan(6))

  END SUBROUTINE deallocate_fft_scratch_type

! *****************************************************************************
!> \brief ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_fft_scratch_pool(error)

    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_fft_scratch_pool', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure
    TYPE(fft_scratch_pool_type), POINTER     :: fft_scratch, &
                                                fft_scratch_current

    failure = .FALSE.

    IF (init_fft_pool == 0) NULLIFY(fft_scratch_first)

    fft_scratch => fft_scratch_first
    DO
       IF (ASSOCIATED(fft_scratch)) THEN
          fft_scratch_current => fft_scratch
          fft_scratch => fft_scratch_current%fft_scratch_next
          NULLIFY(fft_scratch_current%fft_scratch_next)

          CALL deallocate_fft_scratch_type(fft_scratch_current%fft_scratch,error)

          DEALLOCATE(fft_scratch_current%fft_scratch,STAT=ierr)
          CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(fft_scratch_current,STAT=ierr)
          CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
       ELSE
          EXIT
       END IF
    END DO

    init_fft_pool = 0

  END SUBROUTINE release_fft_scratch_pool

! *****************************************************************************
!> \brief ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE resize_fft_scratch_pool(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'resize_fft_scratch_pool', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, last_tick, nscratch
    LOGICAL                                  :: failure
    TYPE(fft_scratch_pool_type), POINTER     :: fft_scratch_current, &
                                                fft_scratch_old

    failure=.FALSE.
    nscratch=0

    last_tick=HUGE(last_tick)
    NULLIFY(fft_scratch_old)

    ! start at the global pool, count, and find a deletion candidate
    fft_scratch_current => fft_scratch_first
    DO
       IF (ASSOCIATED(fft_scratch_current)) THEN
          nscratch=nscratch+1
          ! is this a candidate for deletion (i.e. least recently used, and not in use)
          IF (.NOT.fft_scratch_current%fft_scratch%in_use) THEN
             IF (fft_scratch_current%fft_scratch%last_tick<last_tick) THEN
                last_tick=fft_scratch_current%fft_scratch%last_tick
                fft_scratch_old=>fft_scratch_current
             ENDIF
          ENDIF
          fft_scratch_current=>fft_scratch_current%fft_scratch_next
       ELSE
          EXIT
       ENDIF
    ENDDO

    ! we should delete a scratch
    IF (nscratch>fft_pool_scratch_limit) THEN
       ! note that we never deallocate the first (special) element of the list
       IF (ASSOCIATED(fft_scratch_old)) THEN
          fft_scratch_current => fft_scratch_first
          DO
             IF (ASSOCIATED(fft_scratch_current)) THEN
                ! should we delete the next in the list?
                IF (ASSOCIATED(fft_scratch_current%fft_scratch_next,fft_scratch_old)) THEN
                   ! fix the linked list
                   fft_scratch_current%fft_scratch_next=>fft_scratch_old%fft_scratch_next

                   ! deallocate the element
                   CALL deallocate_fft_scratch_type(fft_scratch_old%fft_scratch,error)
                   DEALLOCATE(fft_scratch_old%fft_scratch,STAT=ierr)
                   CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
                   DEALLOCATE(fft_scratch_old,STAT=ierr)
                   CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

                ELSE
                   fft_scratch_current=>fft_scratch_current%fft_scratch_next
                ENDIF
             ELSE
                EXIT
             ENDIF
          ENDDO

       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "The number of the scratches exceeded the limit, but none could be deallocated")
       ENDIF
    ENDIF

  END SUBROUTINE resize_fft_scratch_pool

! *****************************************************************************
!> \brief ...
!> \param fft_scratch ...
!> \param tf_type ...
!> \param n ...
!> \param fft_sizes ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_fft_scratch(fft_scratch,tf_type,n,fft_sizes,error)
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
    USE ISO_C_BINDING
#endif

    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    INTEGER, INTENT(IN)                      :: tf_type
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    TYPE(fft_scratch_sizes), INTENT(IN), &
      OPTIONAL                               :: fft_sizes
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_fft_scratch', &
      routineP = moduleN//':'//routineN

    INTEGER :: coord(2), DIM(2), handle, i, ierr, ix, iz, lg, lmax, m1, m2, &
      mcx2, mcy3, mcz1, mcz2, mg, mmax, mx1, mx2, my1, my3, mz1, mz2, mz3, &
      nbx, nbz, nm, nmax, nmray, nn, np, nx, ny, nyzray, nz, pos(2)
    INTEGER, DIMENSION(3)                    :: pcoord
    LOGICAL                                  :: equal, failure
    LOGICAL, DIMENSION(2)                    :: dims
    TYPE(fft_scratch_pool_type), POINTER     :: fft_scratch_current, &
                                                fft_scratch_last, &
                                                fft_scratch_new
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
    INTEGER(KIND=C_SIZE_T)    :: length
    TYPE(C_PTR)               :: cptr_r1buf, cptr_tbuf, &
                                 cptr_p2buf, cptr_p3buf, cptr_p4buf, cptr_p5buf, cptr_p7buf
    INTEGER(KIND=C_INT), PARAMETER :: cudaHostAllocDefault=0, &
                                      cudaHostAllocPortable=1, &
                                      cudaHostAllocMapped=2
#endif
    CALL timeset(routineN,handle)

    failure = .FALSE.

    ! this is the place to check that the scratch_pool does not grow without limits
    ! before we add a new scratch check the size of the pool and release some of the list if needed
    CALL resize_fft_scratch_pool(error)

    ! get the required scratch
    tick_fft_pool=tick_fft_pool+1
    fft_scratch_current => fft_scratch_first
    DO
       IF (ASSOCIATED(fft_scratch_current)) THEN
          IF(fft_scratch_current%fft_scratch%in_use) THEN
             fft_scratch_last => fft_scratch_current
             fft_scratch_current => fft_scratch_current%fft_scratch_next
             CYCLE
          END IF
          IF(tf_type /= fft_scratch_current%fft_scratch%tf_type) THEN
             fft_scratch_last => fft_scratch_current
             fft_scratch_current => fft_scratch_current%fft_scratch_next
             CYCLE
          END IF
          IF(.NOT.ALL(n==fft_scratch_current%fft_scratch%nfft))THEN
             fft_scratch_last => fft_scratch_current
             fft_scratch_current => fft_scratch_current%fft_scratch_next
             CYCLE
          END IF
          IF(PRESENT(fft_sizes)) THEN
             IF(fft_sizes%gs_group /= fft_scratch_current%fft_scratch%group) THEN
                fft_scratch_last => fft_scratch_current
                fft_scratch_current => fft_scratch_current%fft_scratch_next
                CYCLE
             END IF
             CALL is_equal(fft_sizes,fft_scratch_current%fft_scratch%sizes,equal)
             IF(.NOT. equal) THEN
                fft_scratch_last => fft_scratch_current
                fft_scratch_current => fft_scratch_current%fft_scratch_next
                CYCLE
             ENDIF
          END IF
          ! Success
          fft_scratch => fft_scratch_current%fft_scratch
          fft_scratch_current%fft_scratch%in_use = .TRUE.
          EXIT
       ELSE
          ! We cannot find the scratch type in this pool
          ! Generate a new scratch set
          ALLOCATE (fft_scratch_new,STAT=ierr)
          CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (fft_scratch_new%fft_scratch,STAT=ierr)
          CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
          fft_scratch_new%fft_scratch%group = 0
          NULLIFY(fft_scratch_new%fft_scratch%ziptr)
          NULLIFY(fft_scratch_new%fft_scratch%zoptr)
          NULLIFY(fft_scratch_new%fft_scratch%p1buf)
          NULLIFY(fft_scratch_new%fft_scratch%p2buf)
          NULLIFY(fft_scratch_new%fft_scratch%p3buf)
          NULLIFY(fft_scratch_new%fft_scratch%p4buf)
          NULLIFY(fft_scratch_new%fft_scratch%p5buf)
          NULLIFY(fft_scratch_new%fft_scratch%p6buf)
          NULLIFY(fft_scratch_new%fft_scratch%p7buf)
          NULLIFY(fft_scratch_new%fft_scratch%r1buf)
          NULLIFY(fft_scratch_new%fft_scratch%r2buf)
          NULLIFY(fft_scratch_new%fft_scratch%tbuf)
          NULLIFY(fft_scratch_new%fft_scratch%a1buf)
          NULLIFY(fft_scratch_new%fft_scratch%a2buf)
          NULLIFY(fft_scratch_new%fft_scratch%a3buf)
          NULLIFY(fft_scratch_new%fft_scratch%a4buf)
          NULLIFY(fft_scratch_new%fft_scratch%a5buf)
          NULLIFY(fft_scratch_new%fft_scratch%a6buf)
          NULLIFY(fft_scratch_new%fft_scratch%scount,fft_scratch_new%fft_scratch%rcount,&
               fft_scratch_new%fft_scratch%sdispl,fft_scratch_new%fft_scratch%rdispl)
          NULLIFY(fft_scratch_new%fft_scratch%rr,&
               fft_scratch_new%fft_scratch%ss, fft_scratch_new%fft_scratch%tt)
          NULLIFY(fft_scratch_new%fft_scratch%xzbuf, fft_scratch_new%fft_scratch%yzbuf,&
               fft_scratch_new%fft_scratch%xzbuf_sgl, fft_scratch_new%fft_scratch%yzbuf_sgl)
          NULLIFY(fft_scratch_new%fft_scratch%pgcube)
          NULLIFY(fft_scratch_new%fft_scratch%pgrid)
          NULLIFY(fft_scratch_new%fft_scratch%xcor, fft_scratch_new%fft_scratch%zcor)
          NULLIFY(fft_scratch_new%fft_scratch%pzcoord)
          NULLIFY(fft_scratch_new%fft_scratch%xzcount, fft_scratch_new%fft_scratch%yzcount)
          NULLIFY(fft_scratch_new%fft_scratch%xzdispl, fft_scratch_new%fft_scratch%yzdispl)
          NULLIFY(fft_scratch_new%fft_scratch%rbuf1, fft_scratch_new%fft_scratch%rbuf2,&
               fft_scratch_new%fft_scratch%rbuf3, fft_scratch_new%fft_scratch%rbuf4)
          NULLIFY(fft_scratch_new%fft_scratch%rbuf5, fft_scratch_new%fft_scratch%rbuf6)
          fft_scratch_new%fft_scratch%in=0
          fft_scratch_new%fft_scratch%rsratio=1._dp
          DO i=1,6
             fft_scratch_new%fft_scratch%fft_plan(i)%valid = .FALSE.
          END DO

          fft_scratch_new%fft_scratch%cart_sub_comm=mp_comm_null

          IF ( tf_type .NE. 400 ) THEN
             fft_scratch_new%fft_scratch%sizes=fft_sizes
             np = fft_sizes%numtask
             ALLOCATE ( fft_scratch_new%fft_scratch%scount(0:np-1), fft_scratch_new%fft_scratch%rcount(0:np-1),&
                  fft_scratch_new%fft_scratch%sdispl(0:np-1), fft_scratch_new%fft_scratch%rdispl(0:np-1),&
                  fft_scratch_new%fft_scratch%pgcube(0:np-1, 2),STAT=ierr)
          END IF

          SELECT CASE (tf_type)
          CASE DEFAULT
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          CASE (100)    ! fft3d_pb: full cube distribution
             mx1   = fft_sizes%mx1
             my1   = fft_sizes%my1
             mx2   = fft_sizes%mx2
             mz2   = fft_sizes%mz2
             my3   = fft_sizes%my3
             mz3   = fft_sizes%mz3
             ALLOCATE ( fft_scratch_new%fft_scratch%a1buf(mx1*my1,n(3)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a2buf(n(3),mx1*my1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a3buf(mx2*mz2,n(2)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a4buf(n(2),mx2*mz2),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a5buf(my3*mz3,n(1)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a6buf(n(1),my3*mz3),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             fft_scratch_new%fft_scratch%group = fft_sizes%gs_group

             CALL mp_environ ( nn, dim, pos, fft_sizes%rs_group )
             CALL mp_cart_rank ( fft_sizes%rs_group, pos, fft_scratch_new%fft_scratch%mip )
             fft_scratch_new%fft_scratch%dim = dim
             fft_scratch_new%fft_scratch%pos = pos
             mcz1   = fft_sizes%mcz1
             mcx2   = fft_sizes%mcx2
             mcz2   = fft_sizes%mcz2
             mcy3   = fft_sizes%mcy3
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf1(mx2*my1*mcz2,0:DIM(2)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf2(mx1*my1*mcz2,0:DIM(2)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf3(mx2*mz3*mcy3,0:DIM(1)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf4(mx2*mz2*mcy3,0:DIM(1)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

             dims = (/.TRUE.,.FALSE./)
             CALL mp_cart_sub( fft_sizes%rs_group, dims, fft_scratch_new%fft_scratch%cart_sub_comm(1))
             dims = (/.FALSE.,.TRUE./)
             CALL mp_cart_sub( fft_sizes%rs_group, dims, fft_scratch_new%fft_scratch%cart_sub_comm(2))

             !initialise pgcube
             DO i = 0, DIM(1) - 1
                coord = (/i, pos(2)/)
                CALL mp_cart_rank ( fft_sizes%rs_group, coord, fft_scratch_new%fft_scratch%pgcube (i, 1) )
             END DO
             DO i = 0, DIM(2) - 1
                coord = (/pos(1), i/)
                CALL mp_cart_rank ( fft_sizes%rs_group, coord, fft_scratch_new%fft_scratch%pgcube (i, 2) )
             END DO

             !set up fft plans
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(1), fft_type, FWFFT, .TRUE., n(3), mx1*my1, &
                  fft_scratch_new%fft_scratch%a1buf, fft_scratch_new%fft_scratch%a2buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(2), fft_type, FWFFT, .TRUE., n(2), mx2*mz2, &
                  fft_scratch_new%fft_scratch%a3buf, fft_scratch_new%fft_scratch%a4buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(3), fft_type, FWFFT, .TRUE., n(1), my3*mz3, &
                  fft_scratch_new%fft_scratch%a5buf, fft_scratch_new%fft_scratch%a6buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(4), fft_type, BWFFT, .TRUE., n(1), my3*mz3, &
                  fft_scratch_new%fft_scratch%a6buf, fft_scratch_new%fft_scratch%a5buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(5), fft_type, BWFFT, .TRUE., n(2), mx2*mz2, &
                  fft_scratch_new%fft_scratch%a4buf, fft_scratch_new%fft_scratch%a3buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(6), fft_type, BWFFT, .TRUE., n(3), mx1*my1, &
                  fft_scratch_new%fft_scratch%a2buf, fft_scratch_new%fft_scratch%a1buf, fft_plan_style)


          CASE (101)    ! fft3d_pb: full cube distribution (dim 1)
             mx1   = fft_sizes%mx1
             my1   = fft_sizes%my1
             mz1   = fft_sizes%mz1
             my3   = fft_sizes%my3
             mz3   = fft_sizes%mz3
             ALLOCATE ( fft_scratch_new%fft_scratch%a1buf(mx1*my1,n(3)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a2buf(n(3),mx1*my1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             fft_scratch_new%fft_scratch%group = fft_sizes%gs_group
             ALLOCATE ( fft_scratch_new%fft_scratch%a3buf(mx1*mz1,n(2)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a4buf(n(2),mx1*mz1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a5buf(my3*mz3,n(1)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%a6buf(n(1),my3*mz3),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

             CALL mp_environ ( nn, dim, pos, fft_sizes%rs_group )
             CALL mp_cart_rank ( fft_sizes%rs_group, pos, fft_scratch_new%fft_scratch%mip )
             fft_scratch_new%fft_scratch%dim = dim
             fft_scratch_new%fft_scratch%pos = pos
             mcy3   = fft_sizes%mcy3
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf5(mx1*mz3*mcy3,0:DIM(1)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf6(mx1*mz1*mcy3,0:DIM(1)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

             !set up fft plans
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(1), fft_type, FWFFT, .TRUE., n(3), mx1*my1, &
                  fft_scratch_new%fft_scratch%a1buf, fft_scratch_new%fft_scratch%a3buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(2), fft_type, FWFFT, .TRUE., n(2), mx1*mz1, &
                  fft_scratch_new%fft_scratch%a3buf, fft_scratch_new%fft_scratch%a4buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(3), fft_type, FWFFT, .TRUE., n(1), my3*mz3, &
                  fft_scratch_new%fft_scratch%a5buf, fft_scratch_new%fft_scratch%a6buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(4), fft_type, BWFFT, .TRUE., n(1), my3*mz3, &
                  fft_scratch_new%fft_scratch%a6buf, fft_scratch_new%fft_scratch%a5buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(5), fft_type, BWFFT, .TRUE., n(2), mx1*mz1, &
                  fft_scratch_new%fft_scratch%a4buf, fft_scratch_new%fft_scratch%a3buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(6), fft_type, BWFFT, .TRUE., n(3), mx1*my1, &
                  fft_scratch_new%fft_scratch%a3buf, fft_scratch_new%fft_scratch%a1buf, fft_plan_style)


          CASE (200)    ! fft3d_ps: plane distribution
             nx    = fft_sizes%nx
             ny    = fft_sizes%ny
             nz    = fft_sizes%nz
             mx2   = fft_sizes%mx2
             lmax  = fft_sizes%lmax
             mmax  = fft_sizes%mmax
             lg    = fft_sizes%lg
             mg    = fft_sizes%mg
             np    = fft_sizes%numtask
             nmray = fft_sizes%nmray
             nyzray= fft_sizes%nyzray
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
             length = INT(2 * dp_size * MAX(mmax,1) * MAX(lmax,1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_r1buf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_r1buf, fft_scratch_new%fft_scratch%r1buf, (/MAX(mmax,1),MAX(lmax,1)/))
             length = INT(2 * dp_size * MAX(ny,1) * MAX(nz,1) * MAX(nx,1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_tbuf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_tbuf, fft_scratch_new%fft_scratch%tbuf, (/MAX(ny,1),MAX(nz,1),MAX(nx,1)/))
#else
             ALLOCATE ( fft_scratch_new%fft_scratch%r1buf(mmax,lmax),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%tbuf(ny,nz,nx),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
#endif
             fft_scratch_new%fft_scratch%group = fft_sizes%gs_group
             ALLOCATE ( fft_scratch_new%fft_scratch%r2buf(lg,mg),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             nm = nmray*mx2
             IF ( alltoall_sgl ) THEN
                ALLOCATE ( fft_scratch_new%fft_scratch%ss(mmax,lmax),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
                ALLOCATE ( fft_scratch_new%fft_scratch%tt(nm,0:np-1),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ELSE
                ALLOCATE ( fft_scratch_new%fft_scratch%rr(nm,0:np-1),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             END IF

             !set up fft plans
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(1), fft_type, FWFFT, .TRUE., nz, nx*ny, &
                  fft_scratch_new%fft_scratch%tbuf, fft_scratch_new%fft_scratch%r1buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(2), fft_type, FWFFT, .TRUE., ny, nx*nz, &
                  fft_scratch_new%fft_scratch%r1buf, fft_scratch_new%fft_scratch%tbuf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(3), fft_type, FWFFT, .TRUE., n(1), nyzray, &
                  fft_scratch_new%fft_scratch%r1buf, fft_scratch_new%fft_scratch%r2buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(4), fft_type, BWFFT, .TRUE., n(1), nyzray, &
                  fft_scratch_new%fft_scratch%r2buf, fft_scratch_new%fft_scratch%r1buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(5), fft_type, BWFFT, .TRUE., ny, nx*nz, &
                  fft_scratch_new%fft_scratch%tbuf, fft_scratch_new%fft_scratch%r1buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(6), fft_type, BWFFT, .TRUE., nz, nx*ny, &
                  fft_scratch_new%fft_scratch%r1buf, fft_scratch_new%fft_scratch%tbuf, fft_plan_style)

          CASE (300)    ! fft3d_ps: block distribution
             mx1   = fft_sizes%mx1
             mx2   = fft_sizes%mx2
             my1   = fft_sizes%my1
             mz2   = fft_sizes%mz2
             mcx2  = fft_sizes%mcx2
             lg    = fft_sizes%lg
             mg    = fft_sizes%mg
             nmax  = fft_sizes%nmax
             nmray = fft_sizes%nmray
             nyzray= fft_sizes%nyzray
             m1    = fft_sizes%r_dim(1)
             m2    = fft_sizes%r_dim(2)
             nbx   = fft_sizes%nbx
             nbz   = fft_sizes%nbz
             ALLOCATE ( fft_scratch_new%fft_scratch%p1buf(mx1*my1,n(3)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%p6buf(lg,mg),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
#if defined ( __PW_CUDA ) && !defined ( __PW_CUDA_NO_HOSTALLOC )
             length = INT(2 * dp_size * MAX(n(3),1) * MAX(mx1*my1,1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_p2buf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_p2buf, fft_scratch_new%fft_scratch%p2buf, (/MAX(n(3),1),MAX(mx1*my1,1)/))
             length = INT(2 * dp_size * MAX(mx2*mz2,1) * MAX(n(2),1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_p3buf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_p3buf, fft_scratch_new%fft_scratch%p3buf, (/MAX(mx2*mz2,1),MAX(n(2),1)/))
             length = INT(2 * dp_size * MAX(n(2),1) * MAX(mx2*mz2,1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_p4buf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_p4buf, fft_scratch_new%fft_scratch%p4buf, (/MAX(n(2),1),MAX(mx2*mz2,1)/))
             length = INT(2 * dp_size * MAX(nyzray,1) * MAX(n(1),1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_p5buf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_p5buf, fft_scratch_new%fft_scratch%p5buf, (/MAX(nyzray,1),MAX(n(1),1)/))
             length = INT(2 * dp_size * MAX(mg,1) * MAX(lg,1), KIND=C_SIZE_T)
             ierr = cudaHostAlloc(cptr_p7buf, length, cudaHostAllocDefault)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             CALL c_f_pointer(cptr_p7buf, fft_scratch_new%fft_scratch%p7buf, (/MAX(mg,1),MAX(lg,1)/))
#else
             ALLOCATE ( fft_scratch_new%fft_scratch%p2buf(n(3),mx1*my1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%p3buf(mx2*mz2,n(2)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%p4buf(n(2),mx2*mz2),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%p5buf(nyzray,n(1)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%p7buf(mg,lg),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
#endif
             IF ( alltoall_sgl ) THEN
                ALLOCATE ( fft_scratch_new%fft_scratch%yzbuf_sgl(mg*lg),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
                ALLOCATE ( fft_scratch_new%fft_scratch%xzbuf_sgl(n(2)*mx2*mz2),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ELSE
                ALLOCATE ( fft_scratch_new%fft_scratch%yzbuf(mg*lg),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
                ALLOCATE ( fft_scratch_new%fft_scratch%xzbuf(n(2)*mx2*mz2),STAT=ierr)
                CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             END IF
             ALLOCATE ( fft_scratch_new%fft_scratch%pgrid(0:m1-1,0:m2-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%xcor(nbx),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%zcor(nbz),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%pzcoord(0:np-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%xzcount(0:np-1), &
                  fft_scratch_new%fft_scratch%yzcount(0:np-1) )
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%xzdispl(0:np-1), &
                  fft_scratch_new%fft_scratch%yzdispl(0:np-1) )
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             fft_scratch_new%fft_scratch%in=0
             fft_scratch_new%fft_scratch%rsratio=1._dp
             fft_scratch_new%fft_scratch%group = fft_sizes%gs_group

             CALL mp_environ ( nn, dim, pos, fft_sizes%rs_group )
             CALL mp_cart_rank ( fft_sizes%rs_group, pos, fft_scratch_new%fft_scratch%mip )
             fft_scratch_new%fft_scratch%dim = dim
             fft_scratch_new%fft_scratch%pos = pos
             mcz1   = fft_sizes%mcz1
             mcz2   = fft_sizes%mcz2
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf1(mx2*my1*mcz2,0:DIM(2)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%rbuf2(mx1*my1*mcz2,0:DIM(2)-1),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

             dims = (/.FALSE.,.TRUE./)
             CALL mp_cart_sub( fft_sizes%rs_group, dims, fft_scratch_new%fft_scratch%cart_sub_comm(2))

             !initialise pgcube
             DO i = 0, DIM(2) - 1
                coord = (/pos(1), i/)
                CALL mp_cart_rank ( fft_sizes%rs_group, coord, fft_scratch_new%fft_scratch%pgcube (i, 2) )
             END DO

             !initialise pgrid
             DO ix = 0, m1 - 1
                DO iz = 0, m2 - 1
                   coord = (/ix, iz/)
                   CALL mp_cart_rank ( fft_sizes%rs_group, coord, fft_scratch_new%fft_scratch%pgrid ( ix, iz ) )
                END DO
             END DO

             !initialise pzcoord
             DO i = 0, np -1
                CALL mp_cart_coords( fft_sizes%rs_group, i, pcoord)
                fft_scratch_new%fft_scratch%pzcoord(i) = pcoord(2)
             END DO

             !set up fft plans
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(1), fft_type, FWFFT, .TRUE., n(3), mx1*my1, &
                  fft_scratch_new%fft_scratch%p1buf, fft_scratch_new%fft_scratch%p2buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(2), fft_type, FWFFT, .TRUE., n(2), mx2*mz2, &
                  fft_scratch_new%fft_scratch%p3buf, fft_scratch_new%fft_scratch%p4buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(3), fft_type, FWFFT, .TRUE., n(1), nyzray, &
                  fft_scratch_new%fft_scratch%p5buf, fft_scratch_new%fft_scratch%p6buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(4), fft_type, BWFFT, .TRUE., n(1), nyzray, &
                  fft_scratch_new%fft_scratch%p6buf, fft_scratch_new%fft_scratch%p7buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(5), fft_type, BWFFT, .TRUE., n(2), mx2*mz2, &
                  fft_scratch_new%fft_scratch%p4buf, fft_scratch_new%fft_scratch%p3buf, fft_plan_style)
             CALL fft_create_plan_1dm(fft_scratch_new%fft_scratch%fft_plan(6), fft_type, BWFFT, .TRUE., n(3), mx1*my1, &
                  fft_scratch_new%fft_scratch%p3buf, fft_scratch_new%fft_scratch%p1buf, fft_plan_style)


          CASE (400)    ! serial FFT
             np = 0
             ALLOCATE ( fft_scratch_new%fft_scratch%ziptr(n(1),n(2),n(3)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
             ALLOCATE ( fft_scratch_new%fft_scratch%zoptr(n(1),n(2),n(3)),STAT=ierr)
             CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

             !in place plans
             CALL fft_create_plan_3d(fft_scratch_new%fft_scratch%fft_plan(1), fft_type, .TRUE., FWFFT, n, &
                  fft_scratch_new%fft_scratch%ziptr, fft_scratch_new%fft_scratch%ziptr, fft_plan_style)
             CALL fft_create_plan_3d(fft_scratch_new%fft_scratch%fft_plan(2), fft_type, .TRUE., BWFFT, n, &
                  fft_scratch_new%fft_scratch%ziptr, fft_scratch_new%fft_scratch%ziptr, fft_plan_style)
             ! out of place plans
             CALL fft_create_plan_3d(fft_scratch_new%fft_scratch%fft_plan(3), fft_type, .FALSE., FWFFT, n, &
                  fft_scratch_new%fft_scratch%ziptr, fft_scratch_new%fft_scratch%zoptr, fft_plan_style)
             CALL fft_create_plan_3d(fft_scratch_new%fft_scratch%fft_plan(4), fft_type, .FALSE., BWFFT, n, &
                  fft_scratch_new%fft_scratch%ziptr, fft_scratch_new%fft_scratch%zoptr, fft_plan_style)

          END SELECT

          NULLIFY(fft_scratch_new%fft_scratch_next)
          fft_scratch_new%fft_scratch%fft_scratch_id = &
               fft_scratch_last%fft_scratch%fft_scratch_id + 1
          fft_scratch_new%fft_scratch%in_use=.TRUE.
          fft_scratch_new%fft_scratch%nfft = n
          fft_scratch_last%fft_scratch_next => fft_scratch_new
          fft_scratch_new%fft_scratch%tf_type = tf_type
          fft_scratch => fft_scratch_new%fft_scratch
          EXIT

       END IF
    END DO

    fft_scratch%last_tick=tick_fft_pool

    CALL timestop(handle)

  END SUBROUTINE get_fft_scratch

! *****************************************************************************
!> \brief ...
!> \param fft_scratch ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_fft_scratch(fft_scratch,error)

    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_fft_scratch', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: scratch_id
    LOGICAL                                  :: failure
    TYPE(fft_scratch_pool_type), POINTER     :: fft_scratch_current

    FAILURE = .FALSE.

    scratch_id = fft_scratch%fft_scratch_id

    fft_scratch_current => fft_scratch_first
    DO
       IF (ASSOCIATED(fft_scratch_current)) THEN
          IF(scratch_id==fft_scratch_current%fft_scratch%fft_scratch_id) THEN
             fft_scratch%in_use = .FALSE.
             NULLIFY(fft_scratch)
             EXIT
          END IF
          fft_scratch_current => fft_scratch_current%fft_scratch_next
       ELSE
          ! We cannot find the scratch type in this pool
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          EXIT
       END IF
    END DO

  END SUBROUTINE release_fft_scratch

! *****************************************************************************
!> \brief ...
!> \param rs ...
!> \param scount ...
!> \param sdispl ...
!> \param rq ...
!> \param rcount ...
!> \param rdispl ...
!> \param group ...
! *****************************************************************************
 SUBROUTINE sparse_alltoall ( rs, scount, sdispl, rq, rcount, rdispl, group )
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: rs
    INTEGER, DIMENSION(:), POINTER           :: scount, sdispl
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: rq
    INTEGER, DIMENSION(:), POINTER           :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: msgin, msgout
    INTEGER                                  :: ip, n, nr, ns, pos, rn, sn
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rreq, sreq

    CALL mp_sync ( group )
    CALL mp_environ ( n, pos, group )
    ALLOCATE(sreq(0:n-1))
    ALLOCATE(rreq(0:n-1))
    nr = 0
    DO ip = 0,n-1
       IF ( rcount(ip) == 0 ) CYCLE
       IF ( ip == pos ) CYCLE
       msgout => rq(rdispl(ip)+1:rdispl(ip)+rcount(ip))
       CALL mp_irecv(msgout,ip,group,rn)
       nr = nr+1
       rreq(nr-1) = rn
    END DO
    ns = 0
    DO ip = 0,n-1
       IF ( scount(ip) == 0 ) CYCLE
       IF ( ip == pos ) CYCLE
       msgin => rs(sdispl(ip)+1:sdispl(ip)+scount(ip))
       CALL mp_isend(msgin,ip,group,sn)
       ns = ns+1
       sreq(ns-1) = sn
    END DO
    IF ( rcount(pos) /= 0 ) THEN
       IF( rcount(pos) /= scount(pos) ) STOP
       rq(rdispl(pos)+1:rdispl(pos)+rcount(pos)) = rs(sdispl(pos)+1:sdispl(pos)+scount(pos))
    END IF
    CALL mp_waitall(sreq(0:ns-1))
    CALL mp_waitall(rreq(0:nr-1))
    DEALLOCATE(sreq)
    DEALLOCATE(rreq)
    CALL mp_sync ( group )

  END SUBROUTINE sparse_alltoall


! *****************************************************************************
!> \brief  test data structures for equality. It is assumed that if they are
!>         different for one mpi task they are different for all (??)
!> \param fft_size_1 ...
!> \param fft_size_2 ...
!> \param equal ...
! *****************************************************************************
  SUBROUTINE is_equal(fft_size_1,fft_size_2,equal)
    TYPE(fft_scratch_sizes)                  :: fft_size_1, fft_size_2
    LOGICAL                                  :: equal

    equal=.TRUE.

    equal=equal.AND.fft_size_1%nx==fft_size_2%nx
    equal=equal.AND.fft_size_1%ny==fft_size_2%ny
    equal=equal.AND.fft_size_1%nz==fft_size_2%nz

    equal=equal.AND.fft_size_1%lmax==fft_size_2%lmax
    equal=equal.AND.fft_size_1%mmax==fft_size_2%mmax
    equal=equal.AND.fft_size_1%nmax==fft_size_2%nmax

    equal=equal.AND.fft_size_1%mx1==fft_size_2%mx1
    equal=equal.AND.fft_size_1%mx2==fft_size_2%mx2
    equal=equal.AND.fft_size_1%mx3==fft_size_2%mx3

    equal=equal.AND.fft_size_1%my1==fft_size_2%my1
    equal=equal.AND.fft_size_1%my2==fft_size_2%my2
    equal=equal.AND.fft_size_1%my3==fft_size_2%my3

    equal=equal.AND.fft_size_1%mcz1==fft_size_2%mcz1
    equal=equal.AND.fft_size_1%mcx2==fft_size_2%mcx2
    equal=equal.AND.fft_size_1%mcz2==fft_size_2%mcz2
    equal=equal.AND.fft_size_1%mcy3==fft_size_2%mcy3

    equal=equal.AND.fft_size_1%lg==fft_size_2%lg
    equal=equal.AND.fft_size_1%mg==fft_size_2%mg

    equal=equal.AND.fft_size_1%nbx==fft_size_2%nbx
    equal=equal.AND.fft_size_1%nbz==fft_size_2%nbz

    equal=equal.AND.fft_size_1%nmray==fft_size_2%nmray
    equal=equal.AND.fft_size_1%nyzray==fft_size_2%nyzray

    equal=equal.AND.fft_size_1%gs_group==fft_size_2%gs_group
    equal=equal.AND.fft_size_1%rs_group==fft_size_2%rs_group

    equal=equal.AND.ALL(fft_size_1%g_pos==fft_size_2%g_pos)
    equal=equal.AND.ALL(fft_size_1%r_pos==fft_size_2%r_pos)
    equal=equal.AND.ALL(fft_size_1%r_dim==fft_size_2%r_dim)

    equal=equal.AND.fft_size_1%numtask==fft_size_2%numtask

  END SUBROUTINE is_equal

END MODULE fft_tools

