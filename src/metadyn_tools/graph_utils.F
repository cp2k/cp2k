!-----------------------------------------------------------------------------!
!   FES: a fast and general program to map metadynamics on grids              !
!   Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013 !
!                 Teodoro Laino                                               !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Module containing utils for mapping FESs
!> \author Teodoro Laino [tlaino] - 06.2009
!> \par History
!>     06.2009 created [tlaino]
!>     teodoro.laino .at. gmail.com 
!>
!> \par Note
!>     Please report any bug to the author
! *****************************************************************************
MODULE graph_utils
  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  TYPE mep_input_data_type
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: minima
     INTEGER                                :: max_iter
     INTEGER                                :: nreplica
     REAL(KIND=dp)                          :: kb
  END TYPE mep_input_data_type

  PUBLIC :: get_val_res,&
            mep_input_data_type,&
            point_pbc,&
            point_no_pbc,&
            derivative,&
            pbc

CONTAINS

! *****************************************************************************
!> \brief computes the derivative of the FES w.r.t CVs
!> \param fes ...
!> \param pos0 ...
!> \param iperd ...
!> \param ndim ...
!> \param ngrid ...
!> \param dp_grid ...
!> \retval der ...
!> \par History
!>      06.2009 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  FUNCTION derivative(fes, pos0, iperd, ndim, ngrid, dp_grid) RESULT(der)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fes
    INTEGER, DIMENSION(:), INTENT(IN)        :: pos0, iperd
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), INTENT(IN)        :: ngrid
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: dp_grid
    REAL(KIND=dp), DIMENSION(ndim)           :: der

    INTEGER                                  :: i, j, pnt
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pos

    ALLOCATE(pos(ndim))
    pos(:) = pos0
    DO i = 1, ndim
       der(i) = 0.0_dp
       DO j = 1,-1,-2
          pos(i) = pos0(i) + j
          pnt    = point_pbc(pos,iperd,ngrid,ndim)
          der(i) = der(i) + REAL(j,KIND=dp)*(-fes(pnt))
       END DO
       pos(i) = pos0(i)
       der(i) = der(i) / (2.0_dp*dp_grid(i))
    END DO
    DEALLOCATE(pos)

  END FUNCTION derivative

! *****************************************************************************
!> \brief Computes the pointer to the 1D array given the n-dimensional position
!>        PBC version
!> \param pos ...
!> \param iperd ...
!> \param ngrid ...
!> \param ndim ...
!> \retval pnt ...
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  FUNCTION point_pbc (pos,iperd,ngrid,ndim) RESULT(pnt)
    INTEGER, DIMENSION(:), INTENT(IN)        :: pos, iperd, ngrid
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER                                  :: pnt

    INTEGER                                  :: idim, lpnt

    idim=1
    pnt=pos(idim)
    IF (iperd(idim) == 1) THEN
       lpnt = pos(idim)
       lpnt = 1000*ngrid(idim)+lpnt
       lpnt = MOD(lpnt,ngrid(idim))
       IF (lpnt==0) lpnt=ngrid(idim)
       pnt = lpnt
    END IF
    DO idim=2,ndim
       lpnt = pos(idim)
       IF (iperd(idim) == 1) THEN
          lpnt = 1000*ngrid(idim)+lpnt
          lpnt = MOD(lpnt,ngrid(idim))
          IF (lpnt==0) lpnt=ngrid(idim)
       END IF
       pnt=pnt+(lpnt-1) * PRODUCT(ngrid(1:idim-1))
    END DO

  END FUNCTION point_pbc

! *****************************************************************************
!> \brief Computes the pointer to the 1D array given the n-dimensional position
!>        PBC version
!> \param pos ...
!> \param iperd ...
!> \param ngrid ...
!> \param ndim ...
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE pbc (pos,iperd,ngrid,ndim)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: pos
    INTEGER, DIMENSION(:), INTENT(IN)        :: iperd, ngrid
    INTEGER, INTENT(IN)                      :: ndim

    INTEGER                                  :: idim, lpnt

    DO idim=1,ndim
       IF (iperd(idim) == 1) THEN
          lpnt = pos(idim)
          lpnt = 1000*ngrid(idim)+lpnt
          lpnt = MOD(lpnt,ngrid(idim))
          IF (lpnt==0) lpnt=ngrid(idim)
          pos(idim) = lpnt
       END IF
    END DO
  END SUBROUTINE pbc

! *****************************************************************************
!> \brief Computes the pointer to the 1D array given the n-dimensional position
!>        non-PBC version
!> \param pos ...
!> \param ngrid ...
!> \param ndim ...
!> \retval pnt ...
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  FUNCTION point_no_pbc (pos, ngrid, ndim) RESULT(pnt)
    INTEGER, DIMENSION(:), INTENT(IN)        :: pos, ngrid
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER                                  :: pnt

    INTEGER                                  :: i

    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO

  END FUNCTION point_no_pbc

! *****************************************************************************
!> \brief Parser informations from the cp2k input/restart
!> \param unit ...
!> \param section ...
!> \param keyword ...
!> \param subsection ...
!> \param i_val ...
!> \param r_val ...
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE get_val_res(unit, section, keyword, subsection, i_val, r_val)
    INTEGER, INTENT(IN)                      :: unit
    CHARACTER(len=*)                         :: section
    CHARACTER(len=*), OPTIONAL               :: keyword, subsection
    INTEGER, INTENT(OUT), OPTIONAL           :: i_val
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: r_val

    CHARACTER(len=512)                       :: line
    INTEGER                                  :: my_ind, stat

    REWIND(unit)
    CALL search(unit,TRIM(section),line,stat=stat)

    IF (stat/=0) THEN 
       WRITE(*,*)"Pattern: "//TRIM(section)//" not found in input file!"
       STOP "Search failed!"
    END IF

    IF (PRESENT(keyword)) THEN
       CALL search(unit,TRIM(keyword),line,stat)
       IF (stat/=0) THEN
         ! if the keyword is not found, let's give back values that will trigger a problem..
         IF (PRESENT(i_val)) i_val=-HUGE(1)
         IF (PRESENT(r_val)) r_val=-HUGE(0.0_dp)
       ELSE
         ! Otherwise read the value
         my_ind=INDEX(line,TRIM(keyword))+LEN_TRIM(keyword)+1
         IF (PRESENT(i_val)) READ(line(my_ind:),*)i_val
         IF (PRESENT(r_val)) READ(line(my_ind:),*)r_val
       END IF
    END IF

    IF (PRESENT(subsection)) THEN
       CALL search(unit,TRIM(subsection),line,stat)
    END IF

  END SUBROUTINE get_val_res

  ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param unit ...
!> \param key ...
!> \param line ...
!> \param stat ...
! *****************************************************************************
  SUBROUTINE search(unit, key, line, stat)
    INTEGER, INTENT(in)                      :: unit
    CHARACTER(LEN=*), INTENT(IN)             :: key
    CHARACTER(LEN=512), INTENT(OUT)          :: line
    INTEGER, INTENT(out)                     :: stat

    stat = 99
    DO WHILE (.TRUE.)
       READ(unit,'(A)',ERR=100,END=100) line
       IF (INDEX(line,TRIM(key)) /= 0) THEN 
          stat = 0
          EXIT 
       END IF
    END DO
100 CONTINUE
  END SUBROUTINE search

END MODULE graph_utils
