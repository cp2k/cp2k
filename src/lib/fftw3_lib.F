!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH 23-Jan-2006 : initial version
!>      Adapted for new interface
!>      IAB 09-Jan-2009 : Modified to cache plans in fft_plan_type
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!>      IAB 09-Oct-2009 : Added OpenMP directives to 1D FFT, and planning routines
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************
SUBROUTINE fftw3_get_lengths ( DATA, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: DATA
  INTEGER, INTENT(INOUT)                   :: max_length

  INTEGER, DIMENSION(:), ALLOCATABLE       :: dlocal, idx
  INTEGER                                  :: h, i, j, k, m, maxn, &
                                              maxn_elevens, maxn_fives, &
                                              maxn_sevens, maxn_thirteens, &
                                              maxn_threes, maxn_twos, &
                                              ndata, nmax, number

!------------------------------------------------------------------------------
! compute ndata
!! FFTW can do arbitrary(?) lenghts, maybe you want to limit them to some
!!    powers of small prime numbers though...

  maxn_twos = 15
  maxn_threes = 3
  maxn_fives = 2
  maxn_sevens = 1
  maxn_elevens = 1
  maxn_thirteens = 0
  maxn = 37748736

  ndata = 0
  DO h = 0, maxn_twos
     nmax = HUGE(0) / 2**h
     DO i = 0, maxn_threes
        DO j = 0, maxn_fives
           DO k = 0, maxn_sevens
              DO m = 0, maxn_elevens
                 number = (3**i) * (5**j) * (7**k) * (11**m)

                 IF ( number > nmax ) CYCLE

                 number = number * 2 ** h
                 IF ( number >= maxn ) CYCLE

                 ndata = ndata + 1
              END DO
           END DO
        END DO
     END DO
  END DO

  ALLOCATE ( dlocal ( ndata ), idx ( ndata ) )

  ndata = 0
  dlocal ( : ) = 0
  DO h = 0, maxn_twos
     nmax = HUGE(0) / 2**h
     DO i = 0, maxn_threes
        DO j = 0, maxn_fives
           DO k = 0, maxn_sevens
              DO m = 0, maxn_elevens
                 number = (3**i) * (5**j) * (7**k) * (11**m)

                 IF ( number > nmax ) CYCLE

                 number = number * 2 ** h
                 IF ( number >= maxn ) CYCLE

                 ndata = ndata + 1
                 dlocal ( ndata ) = number
              END DO
           END DO
        END DO
     END DO
  END DO

  CALL sortint ( dlocal, ndata, idx )
  ndata = MIN ( ndata, max_length )
  DATA(1:ndata) = dlocal(1:ndata)
  max_length = ndata

  DEALLOCATE ( dlocal, idx )

END SUBROUTINE fftw3_get_lengths

! *****************************************************************************

SUBROUTINE fftw3_create_plan_3d(plan, zin, zout, plan_style, aligned)

  USE fft_plan,                           ONLY : fft_plan_type, dp

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zout
  INTEGER                                            :: plan_style
  LOGICAL                                            :: aligned
  INTEGER                                            :: n1,n2,n3

#if defined ( __FFTW3 )
#include "fftw3.f"
  INTEGER :: fftw_plan_type
  SELECT CASE(plan_style)
  CASE(1)
         fftw_plan_type = FFTW_ESTIMATE 
  CASE(2)
         fftw_plan_type = FFTW_MEASURE
  CASE(3)
         fftw_plan_type = FFTW_PATIENT
  CASE(4)
         fftw_plan_type = FFTW_EXHAUSTIVE
  CASE DEFAULT
   STOP "fftw3_create_plan_3d"
  END SELECT

  IF ( .NOT. aligned ) THEN
    fftw_plan_type = fftw_plan_type + FFTW_UNALIGNED
  END IF

  n1 = plan%n_3d(1)
  n2 = plan%n_3d(2)
  n3 = plan%n_3d(3)

#if defined(__SGL) || defined(__FFTSGL)
  IF ( plan%fft_in_place .AND. plan%fsign == +1 ) THEN
    CALL sfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zin,FFTW_FORWARD,fftw_plan_type)
  ELSE IF ( plan%fft_in_place .AND. plan%fsign == -1 ) THEN
    CALL sfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zin,FFTW_BACKWARD,fftw_plan_type)
  ELSE IF ( .NOT. plan%fft_in_place .AND. plan%fsign == +1 ) THEN
    CALL sfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zout,FFTW_FORWARD,fftw_plan_type)
  ELSE IF ( .NOT. plan%fft_in_place .AND. plan%fsign == -1 ) THEN
    CALL sfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zout,FFTW_BACKWARD,fftw_plan_type)
  END IF
#else
  IF ( plan%fft_in_place .AND. plan%fsign == +1 ) THEN
    CALL dfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zin,FFTW_FORWARD,fftw_plan_type)
  ELSE IF ( plan%fft_in_place .AND. plan%fsign == -1 ) THEN
    CALL dfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zin,FFTW_BACKWARD,fftw_plan_type)
  ELSE IF ( .NOT. plan%fft_in_place .AND. plan%fsign == +1 ) THEN
    CALL dfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zout,FFTW_FORWARD,fftw_plan_type)
  ELSE IF ( .NOT. plan%fft_in_place .AND. plan%fsign == -1 ) THEN
    CALL dfftw_plan_dft_3d(plan%fftw_plan,n1,n2,n3,zin,zout,FFTW_BACKWARD,fftw_plan_type)
  END IF
#endif

#endif

END SUBROUTINE fftw3_create_plan_3d

SUBROUTINE fftw33d ( plan, scale, zin, zout, stat )

  USE fft_plan,                           ONLY : fft_plan_type, dp

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT(IN)                      :: plan
  REAL(KIND=dp), INTENT(IN)                            :: scale
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zout
  INTEGER, INTENT(OUT)                                 :: stat

  INTEGER                                  :: n1, n2, n3

!------------------------------------------------------------------------------

  n1 = plan%n_3d(1)
  n2 = plan%n_3d(2)
  n3 = plan%n_3d(3)

#if defined ( __FFTW3 )
  stat = 1

#if defined(__SGL) || defined(__FFTSGL)
  IF ( plan%fft_in_place ) THEN
    CALL sfftw_execute_dft(plan%fftw_plan, zin, zin)
  ELSE
    CALL sfftw_execute_dft(plan%fftw_plan, zin, zout)
  END IF
#else
  IF ( plan%fft_in_place ) THEN
    CALL dfftw_execute_dft(plan%fftw_plan, zin, zin)
  ELSE
    CALL dfftw_execute_dft(plan%fftw_plan, zin, zout)
  END IF
#endif

  IF ( scale /= 1.0_dp ) THEN
#if defined(__SGL) || defined(__FFTSGL)
    IF ( plan%fft_in_place ) THEN
      CALL csscal(n1*n2*n3,scale,zin,1)
    ELSE
      CALL csscal(n1*n2*n3,scale,zout,1)
    END IF
#else
    IF ( plan%fft_in_place ) THEN
      CALL zdscal(n1*n2*n3,scale,zin,1)
    ELSE
      CALL zdscal(n1*n2*n3,scale,zout,1)
    END IF
#endif
  END IF

#else

  stat = 0

#endif

END SUBROUTINE fftw33d

! *****************************************************************************

SUBROUTINE fftw3_create_plan_1dm(plan, zin, zout, plan_style, aligned)

  USE fft_plan,                           ONLY : fft_plan_type, dp

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zout
  INTEGER                                            :: plan_style
  LOGICAL                                            :: aligned
  INTEGER                                            :: ii,di,io,DO, num_threads, num_rows

!$  INTEGER :: omp_get_num_threads, omp_get_thread_num

#if defined ( __FFTW3 )
#include "fftw3.f"
  INTEGER :: fftw_plan_type 
  SELECT CASE(plan_style)
  CASE(1)
         fftw_plan_type = FFTW_ESTIMATE
  CASE(2)
         fftw_plan_type = FFTW_MEASURE
  CASE(3)
         fftw_plan_type = FFTW_PATIENT
  CASE(4)
         fftw_plan_type = FFTW_EXHAUSTIVE
  CASE DEFAULT
   STOP "fftw3_create_plan_1dm"
  END SELECT

  IF ( .NOT. aligned ) THEN
    fftw_plan_type = fftw_plan_type + FFTW_UNALIGNED
  END IF

num_threads = 1
!$omp parallel default(none), &
!$omp          shared(num_threads)
!$ num_threads = omp_get_num_threads()
!$omp end parallel

num_rows = plan%m / num_threads
!$   plan%num_threads_needed = num_threads

! Check for number of rows less than num_threads
!$ IF (plan%m < num_threads) THEN
!$   num_rows = 1
!$   plan%num_threads_needed = plan%m
!$ ENDIF
#if defined(__SGL) || defined(__FFTSGL)
! Check for odd numbers of elements per thread (for SSE alignment issues)
!$ IF (MOD(num_rows*plan%n, 2) .EQ. 1 .AND. plan%num_threads_needed .NE. 1) THEN
!$   num_rows = num_rows + 1
!$   plan%num_threads_needed = (plan%m + num_rows - 1) / num_rows
!$ ENDIF
#endif
! Check for total number of rows not divisible by num_threads
!$ IF (num_rows*plan%num_threads_needed .NE. plan%m) THEN
!$   plan%need_alt_plan = .TRUE.
!$ ENDIF

!$  plan%num_rows = num_rows
  ii = 1
  di = plan%n
  io = 1
  DO = plan%n
  IF ( plan%fsign == +1 .AND. plan%trans ) THEN
    ii = plan%m
    di = 1
  ELSEIF ( plan%fsign == -1 .AND. plan%trans ) THEN
    io = plan%m
    DO = 1
  END IF

#if defined(__SGL) || defined(__FFTSGL)
  IF ( plan%fsign == +1 ) THEN
    CALL sfftw_plan_many_dft(plan%fftw_plan,1,plan%n,num_rows,zin,0,ii,di,&
              zout,0,io,DO,FFTW_FORWARD,fftw_plan_type)
  ELSE
    CALL sfftw_plan_many_dft(plan%fftw_plan,1,plan%n,num_rows,zin,0,ii,di,&
              zout,0,io,DO,FFTW_BACKWARD,fftw_plan_type)
  END IF
#else
  IF ( plan%fsign == +1 ) THEN
    CALL dfftw_plan_many_dft(plan%fftw_plan,1,plan%n,num_rows,zin,0,ii,di,&
              zout,0,io,DO,FFTW_FORWARD,fftw_plan_type)
  ELSE
    CALL dfftw_plan_many_dft(plan%fftw_plan,1,plan%n,num_rows,zin,0,ii,di,&
              zout,0,io,DO,FFTW_BACKWARD,fftw_plan_type)
  END IF
#endif

!$ IF (plan%need_alt_plan) THEN
!$  plan%alt_num_rows = plan%m - (plan%num_threads_needed - 1)*num_rows
#if defined(__SGL) || defined(__FFTSGL)
!$  IF ( plan%fsign == +1 ) THEN
!$    CALL sfftw_plan_many_dft(plan%alt_fftw_plan,1,plan%n,plan%alt_num_rows,zin,0,ii,di,&
!$              zout,0,io,DO,FFTW_FORWARD,fftw_plan_type)
!$  ELSE
!$    CALL sfftw_plan_many_dft(plan%alt_fftw_plan,1,plan%n,plan%alt_num_rows,zin,0,ii,di,&
!$              zout,0,io,DO,FFTW_BACKWARD,fftw_plan_type)
!$  END IF
#else
!$  IF ( plan%fsign == +1 ) THEN
!$    CALL dfftw_plan_many_dft(plan%alt_fftw_plan,1,plan%n,plan%alt_num_rows,zin,0,ii,di,&
!$              zout,0,io,DO,FFTW_FORWARD,fftw_plan_type)
!$  ELSE
!$    CALL dfftw_plan_many_dft(plan%alt_fftw_plan,1,plan%n,plan%alt_num_rows,zin,0,ii,di,&
!$              zout,0,io,DO,FFTW_BACKWARD,fftw_plan_type)
!$  END IF
#endif
!$ END IF

#endif

END SUBROUTINE fftw3_create_plan_1dm

SUBROUTINE fftw3_destroy_plan ( plan )

  USE fft_plan,                           ONLY: fft_plan_type

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT (INOUT)   :: plan

#if defined ( __FFTW3 )
#if defined(__SGL) || defined(__FFTSGL)
  CALL sfftw_destroy_plan(plan%fftw_plan)
!$  IF (plan%need_alt_plan) THEN
!$    CALL sfftw_destroy_plan(plan%alt_fftw_plan)
!$  END IF
#else
  CALL dfftw_destroy_plan(plan%fftw_plan)
!$  IF (plan%need_alt_plan) THEN
!$    CALL dfftw_destroy_plan(plan%alt_fftw_plan)
!$  END IF
#endif
#endif

END SUBROUTINE fftw3_destroy_plan

SUBROUTINE fftw31dm ( plan, zin, zout, scale, stat )

  USE fft_plan,                 ONLY : fft_plan_type, dp

  IMPLICIT NONE
#if defined (__NAG)
  INTEGER*8          :: dummy_64_bit
  INTEGER, PARAMETER :: integer8_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer8_kind = 8
#endif

  TYPE (fft_plan_type), INTENT (IN)                  :: plan
  COMPLEX(KIND=dp), DIMENSION(*), TARGET, INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), TARGET, INTENT(INOUT)      :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale
  INTEGER, INTENT(OUT)                               :: stat
  INTEGER :: num_rows, my_id, in_offset, out_offset, scal_offset
  COMPLEX(KIND=dp), POINTER :: zin_ptr, zout_ptr, zscal_ptr
  INTEGER (KIND=integer8_kind)        :: fftw_plan

!$  INTEGER :: omp_get_num_threads, omp_get_thread_num

!------------------------------------------------------------------------------

my_id = 0
num_rows = plan%m

!$omp parallel default(none), &
!$omp          private(my_id,num_rows,zin_ptr,zout_ptr,zscal_ptr,in_offset,out_offset,scal_offset,fftw_plan), &
!$omp          shared(zin,zout), &
!$omp          shared(plan,scale,stat)
!$ my_id = omp_get_thread_num()

!$ if (my_id < plan%num_threads_needed) then

fftw_plan = plan%fftw_plan

in_offset = 1
out_offset = 1
scal_offset = 1

!$ in_offset = 1 + plan%num_rows * my_id * plan%n
!$ out_offset = 1 + plan%num_rows * my_id * plan%n
!$ IF ( plan%fsign == +1 .AND. plan%trans ) THEN
!$  in_offset = 1 + plan%num_rows*my_id
!$ ELSEIF ( plan%fsign == -1 .AND. plan%trans ) THEN
!$  out_offset = 1 + plan%num_rows*my_id
!$ ENDIF
!$ scal_offset = 1 + plan%n*plan%num_rows*my_id 
!$ IF ( plan%need_alt_plan .AND. my_id .EQ. omp_get_num_threads() - 1 ) THEN
!$   num_rows = plan%alt_num_rows
!$   fftw_plan = plan%alt_fftw_plan
!$ ELSE
!$   num_rows = plan%num_rows
!$ ENDIF

zin_ptr => zin(in_offset)
zout_ptr => zout(out_offset)
zscal_ptr => zout(scal_offset)

#if defined ( __FFTW3 )
  stat=1
#if defined(__SGL) || defined(__FFTSGL)
  CALL sfftw_execute_dft(fftw_plan, zin_ptr, zout_ptr)
!$ endif
! all theads need to meet at this barrier
!$omp barrier
!$ if (my_id < plan%num_threads_needed) then
  IF ( scale /= 1.0_dp ) CALL csscal(plan%n*num_rows,scale,zscal_ptr,1)
!$ endif
#else
  CALL dfftw_execute_dft(fftw_plan, zin_ptr, zout_ptr)
!$ endif
! all theads need to meet at this barrier
!$omp barrier
!$ if (my_id < plan%num_threads_needed) then
  IF ( scale /= 1.0_dp ) CALL zdscal(plan%n*num_rows,scale,zscal_ptr,1)
!$ endif
#endif

#else

  stat=0

!$ else
!$ end if

#endif

!$omp end parallel

END SUBROUTINE fftw31dm

!     Copyright (c) 2003, 2006 Matteo Frigo
!     Copyright (c) 2003, 2006 Massachusetts Institute of Technology
!     
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     
!     This is an example implementation of Fortran wisdom export/import
!     to/from a Fortran unit (file), exploiting the generic
!     dfftw_export_wisdom/dfftw_import_wisdom functions.
!     
!     We cannot compile this file into the FFTW library itself, lest all
!     FFTW-calling programs be required to link to the Fortran I/O
!     libraries.
!
!     adapted to become more standard Fortran 90 [2007-10] Joost VandeVondele
!     and added some namespacing
!     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! *****************************************************************************
      SUBROUTINE fftw_write_char(c, iunit)
         IMPLICIT NONE
         CHARACTER c
         INTEGER iunit
         WRITE(iunit,'(a)',ADVANCE="NO") c
      END      

! *****************************************************************************
      SUBROUTINE fftw_export_wisdom_to_file(iunit)
         IMPLICIT NONE
         INTEGER iunit
         EXTERNAL fftw_write_char
#if defined ( __FFTW3 )
#if defined (__SGL) || defined(__FFTSGL)
         CALL sfftw_export_wisdom(fftw_write_char, iunit)
#else
         CALL dfftw_export_wisdom(fftw_write_char, iunit)
#endif
#endif
      END

!     Fortran 77 does not have any portable way to read an arbitrary
!     file one character at a time [needs to wait for stream IO of F2003].
!     The best alternative seems to be to
!     read a whole line into a buffer, since for fftw-exported wisdom we
!     can bound the line length.  (If the file contains longer lines,
!     then the lines will be truncated and the wisdom import should
!     simply fail.)  Ugh (and not thread safe).

! *****************************************************************************
      SUBROUTINE fftw_read_char(ic, iunit)
         IMPLICIT NONE
         INTEGER ic
         INTEGER iunit
         CHARACTER*256 buf
         SAVE buf
         INTEGER ibuf
         DATA ibuf/257/
         SAVE ibuf
         IF (ibuf .LT. 257) THEN
            ic = ICHAR(buf(ibuf:ibuf))
            ibuf = ibuf + 1
            RETURN
         ENDIF
         READ(iunit,123,END=666) buf
         ic = ICHAR(buf(1:1))
         ibuf = 2
         RETURN
 666     ic = -1
         ibuf = 257
 123     FORMAT(a256)
      END
      
! *****************************************************************************
      SUBROUTINE fftw_import_wisdom_from_file(isuccess, iunit)
         IMPLICIT NONE
         INTEGER isuccess
         INTEGER iunit
         EXTERNAL fftw_read_char
         isuccess=0
#if defined ( __FFTW3 )
#if defined (__SGL) || defined(__FFTSGL)
         CALL sfftw_import_wisdom(isuccess, fftw_read_char, iunit)
#else
         CALL dfftw_import_wisdom(isuccess, fftw_read_char, iunit)
#endif
#endif
      END
