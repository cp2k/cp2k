!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****h* cp2k/orbital_transformation_matrices [1.0] *
!!
!!   NAME
!!     orbital_transformation_matrices
!!
!!   FUNCTION
!!     Calculation of the spherical harmonics and the corresponding orbital
!!     transformation matrices.
!!
!!   AUTHOR
!!     Matthias Krack (08.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE orbital_transformation_matrices

! *****************************************************************************

! Index:

! SUBROUTINE init_spherical_harmonics(maxl,globenv)
! SUBROUTINE write_matrix(matrix,l,lunit,headline)

! *****************************************************************************

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             pi
  USE mathlib
!,          ONLY: OPERATOR(.over.)

  USE orbital_pointers,                ONLY: co,&
                                             nco,&
                                             nso
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  TYPE orbtramat_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: c2s,slm,s2c
  END TYPE orbtramat_type

  TYPE(orbtramat_type), DIMENSION(:), POINTER :: orbtramat

  INTEGER, SAVE :: current_maxl = -1

! *** Public variables ***

  PUBLIC :: orbtramat

! *** Public data types ***

  PUBLIC :: orbtramat_type

! *** Public subroutines ***

  PUBLIC :: init_spherical_harmonics

!!***
! *****************************************************************************

CONTAINS

!!*****
!!****f* cp2k/init_spherical_harmonics [1.0] *
!!
!!   NAME
!!     init_spherical_harmonics
!!
!!   FUNCTION
!!     Build up the orbital transformation matrices for the
!!     transformation and for the backtransformation of orbitals from
!!     the Cartesian representation to the spherical representation.
!!
!!   NOTES
!!     Literature: H. B. Schlegel, M. J. Frisch,
!!                 Int. J. Quantum Chem. 54, 83 (1995)
!!
!!   HISTORY
!!     Creation (09.07.1999, Matthias Krack)
!!
!!   SOURCE
!!
  SUBROUTINE init_spherical_harmonics(maxl,globenv)

    INTEGER, INTENT(IN)                      :: maxl
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE init_spherical_harmonics (MODULE spherical_harmonics)"

    CHARACTER(LEN=78)                        :: headline
    INTEGER                                  :: expo, i, ic, ic1, ic2, is, &
                                                istat, j, k, l, lx, lx1, lx2, &
                                                ly, ly1, ly2, lz, lz1, lz2, &
                                                m, ma, nc, ns, output_unit
    REAL(KIND=dp)                            :: s, s1, s2

!   ---------------------------------------------------------------------------
!   *** Quick return, if the current initialization is sufficient ***

    IF (maxl <= current_maxl) RETURN

!   *** Deallocate the old orbital transformation matrices ***

    IF (current_maxl > -1) THEN
      DO l=0,SIZE(orbtramat,1)-1
        DEALLOCATE (orbtramat(l)%c2s,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"orbtramat(l)%c2s")
        DEALLOCATE (orbtramat(l)%slm,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"orbtramat(l)%slm")
        DEALLOCATE (orbtramat(l)%s2c,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"orbtramat(l)%s2c")
      END DO
      DEALLOCATE (orbtramat,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"orbtramat")
    END IF

!   *** Allocate and initialize the orbital transformation matrices ***

    ALLOCATE (orbtramat(0:maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"orbtramat",0)

    DO l=0,maxl
      nc = nco(l)
      ns = nso(l)
      ALLOCATE (orbtramat(l)%c2s(ns,nc),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"orbtramat(l)%s2c",ns*nc*dp_size)
      END IF
      orbtramat(l)%c2s(:,:) = 0.0_dp
      ALLOCATE (orbtramat(l)%slm(ns,nc),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"orbtramat(l)%s2c",ns*nc*dp_size)
      END IF
      orbtramat(l)%slm(:,:) = 0.0_dp
      ALLOCATE (orbtramat(l)%s2c(ns,nc),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"orbtramat(l)%s2c",ns*nc*dp_size)
      END IF
      orbtramat(l)%s2c(:,:) = 0.0_dp
    END DO

!   *** Loop over all angular momentum quantum numbers l ***

    DO l=0,maxl

!     *** Build the orbital transformation matrix for the     ***
!     *** transformation from Cartesian to spherical orbitals ***
!     *** (c2s, formula 15)                                   ***

      DO lx=0,l
        DO ly=0,l-lx
          lz = l - lx - ly
          ic = co(lx,ly,lz)
          DO m=-l,l
            is = l + m + 1
            ma = ABS(m)
            j = lx + ly - ma
            IF ((j >= 0).AND.(MODULO(j,2) == 0)) THEN
              j = j/2
              s1 = 0.0_dp
              DO i=0,(l-ma)/2
                s2 = 0.0_dp
                DO k=0,j
                  IF (((m < 0).AND.(MODULO(ABS(ma-lx),2) == 1)).OR.&
                      ((m > 0).AND.(MODULO(ABS(ma-lx),2) == 0))) THEN
                    expo = (ma - lx + 2*k)/2
                    s = (-1.0_dp)**expo*SQRT(2.0_dp)
                  ELSE IF ((m == 0).AND.(MODULO(lx,2) == 0)) THEN
                    expo = k - lx/2
                    s = (-1.0_dp)**expo
                  ELSE
                    s = 0.0_dp
                  END IF
                  s2 = s2 + (j.over.k)*(ma.over.(lx-2*k))*s
                END DO
                s1 = s1 + (l.over.i)*(i.over.j)*(-1.0_dp)**i*fac(2*l-2*i)/&
                          fac(l-ma-2*i)*s2
              END DO
              orbtramat(l)%c2s(is,ic) =&
                SQRT((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l-ma))/&
                     (fac(lx)*fac(ly)*fac(lz)*fac(2*l)*fac(l+ma)))*s1/&
                (2.0_dp**l*fac(l))
            ELSE
              orbtramat(l)%c2s(is,ic) = 0.0_dp
            END IF
          END DO
        END DO
      END DO

!     *** Build the corresponding transformation matrix for the ***
!     *** transformation from spherical to Cartesian orbitals   ***
!     *** (s2c = s*TRANSPOSE(c2s), formulas 18 and 19)          ***

      DO lx1=0,l
        DO ly1=0,l-lx1
          lz1 = l - lx1 - ly1
          ic1 = co(lx1,ly1,lz1)
          s1 = SQRT((fac(lx1)*fac(ly1)*fac(lz1))/&
                    (fac(2*lx1)*fac(2*ly1)*fac(2*lz1)))
          DO lx2=0,l
            DO ly2=0,l-lx2
              lz2 = l - lx2 - ly2
              ic2 = co(lx2,ly2,lz2)
              lx = lx1 + lx2
              ly = ly1 + ly2
              lz = lz1 + lz2
              IF ((MODULO(lx,2) == 0).AND.&
                  (MODULO(ly,2) == 0).AND.&
                  (MODULO(lz,2) == 0)) THEN
                s2 = SQRT((fac(lx2)*fac(ly2)*fac(lz2))/&
                          (fac(2*lx2)*fac(2*ly2)*fac(2*lz2)))
                s = fac(lx)*fac(ly)*fac(lz)*s1*s2/&
                    (fac(lx/2)*fac(ly/2)*fac(lz/2))
                DO is=1,nso(l)
                  orbtramat(l)%s2c(is,ic1) = orbtramat(l)%s2c(is,ic1) +&
                                             s*orbtramat(l)%c2s(is,ic2)
                END DO
              END IF
            END DO
          END DO
        END DO
      END DO

!     *** Build up the real spherical harmonics ***

      s = SQRT(0.25_dp*dfac(2*l+1)/pi)

      DO lx=0,l
        DO ly=0,l-lx
          lz = l - lx - ly
          ic = co(lx,ly,lz)
          DO m=-l,l
            is = l + m + 1
!MK         s2 = (-1.0_dp)**m*s !-> alternative S(lm) definition
            orbtramat(l)%slm(is,ic) =&
              s*orbtramat(l)%c2s(is,ic)/&
!MK              s2*orbtramat(l)%c2s(is,ic)/&
              SQRT(dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))
          END DO
        END DO
      END DO

    END DO

!   *** Print the spherical harmonics and   ***
!   *** the orbital transformation matrices ***

    IF (globenv%ionode.AND.globenv%print%spherical_harmonics) THEN

      output_unit = globenv%scr

      DO l=0,maxl

        nc = nco(l)
        ns = nso(l)

        headline = "SPHERICAL HARMONICS"
        CALL write_matrix(orbtramat(l)%slm,l,output_unit,headline)

        headline = "CARTESIAN ORBITAL TO SPHERICAL ORBITAL "//&
                   "TRANSFORMATION MATRIX"
        CALL write_matrix(orbtramat(l)%c2s,l,output_unit,headline)

        headline = "SPHERICAL ORBITAL TO CARTESIAN ORBITAL "//&
                   "TRANSFORMATION MATRIX"
        CALL write_matrix(orbtramat(l)%s2c,l,output_unit,headline)

      END DO

      WRITE (UNIT=output_unit,FMT="(A)") ""

    END IF

!   *** Save the current initialization status ***

    current_maxl = maxl

  END SUBROUTINE init_spherical_harmonics

! *****************************************************************************

  SUBROUTINE write_matrix(matrix,l,lunit,headline)

!   Purpose: Print a matrix to the logical unit "lunit".

!   History: - Creation (07.06.2000, Matthias Krack)

!   ***************************************************************************

!   List of variables:

!   headline: Headline of the matrix.
!   lunit   : Logical unit number.
!   matrix  : Matrix to be printed.

!   ---------------------------------------------------------------------------

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    INTEGER, INTENT(IN)                      :: l, lunit
    CHARACTER(LEN=*), INTENT(IN)             :: headline

    CHARACTER(12)                            :: symbol
    CHARACTER(LEN=78)                        :: string
    INTEGER                                  :: from, i, ic, is, jc, lx, ly, &
                                                lz, m, nc, to

!   ---------------------------------------------------------------------------
!   *** Write headline ***

    WRITE (UNIT=lunit,FMT="(/,/,T2,A)") TRIM(headline)

!   *** Write the matrix in the defined format ***

    nc = nco(l)

    DO ic=1,nc,6
      from = ic
      to = MIN(nc,from+5)
      i = 1
      DO lx=l,0,-1
        DO ly=l-lx,0,-1
          lz = l - lx - ly
          jc = co(lx,ly,lz)
          IF ((jc >= from).AND.(jc <= to)) THEN
            symbol = cgf_symbol(1,(/lx,ly,lz/))
            WRITE (UNIT=string(i:),FMT="(A)") TRIM(symbol(3:12))
            i = i + 11
          END IF
        END DO
      END DO
      WRITE (UNIT=lunit,FMT="(/,T13,A)") TRIM(string)
      symbol = ""
      DO m=-l,l
        is = l + m + 1
        symbol = sgf_symbol(1,l,m)
        WRITE (UNIT=lunit,FMT="(T4,A4,6(1X,F10.6))")&
          symbol(3:6),(matrix(is,jc),jc=from,to)
      END DO
    END DO

  END SUBROUTINE write_matrix

! *****************************************************************************

END MODULE orbital_transformation_matrices
!!*****
