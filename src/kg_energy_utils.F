!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!


MODULE kg_energy_utils

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE linear_systems,                  ONLY: solve_system
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             test_next
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_energy_section, energy_min_type, min_info_type
  PUBLIC :: print_convergence, print_energy_info, diis                 

  TYPE min_info_type
    REAL ( dbl ) :: e, de, tol, e_hartree, e_xc, e_pol, err, e_ppl, e_ppl2
  END TYPE min_info_type

  TYPE energy_min_type
    CHARACTER ( len=8 ) :: mintype
    INTEGER :: iter
    INTEGER :: ndiis
    REAL ( dbl ) :: max_tol
    REAL ( dbl ) :: sd_delta
    TYPE ( min_info_type ) :: min_info
  END TYPE energy_min_type

CONTAINS

!******************************************************************************

!!>---------------------------------------------------------------------------!
!! SECTION: &energy ... &end                                                  !
!!                                                                            !
!!  mintype:  [DIIS, SD]                 METHOD for minimization              !
!!  sd_delta: delta                      delta t for steepest descent         !
!!  max_tol:  max_tol                    maximum tolerance for convergence    !
!!  iter:     max_iter                   maximum number of iterations         !
!!  ndiis:    ndiis                      number of DIIS vectors               !
!!                                                                            !
!!<---------------------------------------------------------------------------!

SUBROUTINE read_energy_section ( energy_min, globenv )


    TYPE(energy_min_type), INTENT(OUT)       :: energy_min
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=20)                        :: string
    CHARACTER(len=5)                         :: label
    INTEGER                                  :: ierror, ilen, iw

!------------------------------------------------------------------------------

  iw = globenv % scr

!..defaults
  energy_min % max_tol = 1.0E-6_dbl
  energy_min % iter = 10000
  energy_min % ndiis  = 5
  energy_min % sd_delta = 1._dbl
  energy_min % mintype = 'SD'

! initialize min variables
  energy_min % min_info % e = 0._dbl
  energy_min % min_info % de = 0._dbl
  energy_min % min_info % err = 1.0E+6_dbl
  energy_min % min_info % tol = 1.0E+6_dbl
  energy_min % min_info % e_xc = 0._dbl
  energy_min % min_info % e_hartree = 0._dbl
  energy_min % min_info % e_ppl = 0._dbl
  energy_min % min_info % e_ppl2 = 0._dbl

!..parse the input section

  label = '&ENERGY'

  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )

  IF ( ierror /= 0 ) THEN

     IF ( globenv % ionode ) THEN
        CALL stop_program ( 'read_energy_section', 'no input section &ENERGY found' )
     END IF

  ELSE

     CALL read_line()

     DO WHILE ( test_next ( ) /= 'X' )

        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string ) 

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error ( )
           CALL stop_program ( 'read_energy_section', 'unknown option' )

        CASE ( 'MINTYPE' )
           ilen = 8
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           energy_min % mintype = string ( 1:ilen )

        CASE ( 'TOL' )
           CALL get_next ( energy_min % max_tol )

        CASE ( 'ITER' )
           CALL get_next ( energy_min % iter )

        CASE ( 'NDIIS' )
           CALL get_next ( energy_min % ndiis )

        CASE ( 'SD_DELTA' )
           CALL get_next ( energy_min % sd_delta )

        END SELECT

! Get the next line
        CALL read_line

     END DO

  END IF
  CALL parser_end

!..end of parsing the input section
!..write some information to output
  IF ( globenv % print_level >= 0 ) THEN
     WRITE ( iw, '( A )' ) ' ENERGY MINIMIZATION'
     WRITE ( iw, '( A, T71, A )' ) &
          ' ENERGY| Minimization type ', ADJUSTR ( energy_min % mintype )
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' ENERGY| Max Number of iteration ', energy_min % iter
     WRITE ( iw, '( A, T71, G10.4 )' ) &
          ' ENERGY| Max tolerance ', energy_min % max_tol
     WRITE ( iw, '( A, T71, G10.4 )' ) &
          ' ENERGY| delta         ', energy_min % sd_delta
     IF (energy_min % mintype == 'DIIS' ) THEN
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' ENERGY| number of diis vectors ', energy_min % ndiis
     END IF
  END IF

END SUBROUTINE read_energy_section

!******************************************************************************
!
  SUBROUTINE print_energy_info ( iter, min_info, min_type, globenv )

    INTEGER, INTENT(IN)                      :: iter
    TYPE(min_info_type), INTENT(IN)          :: min_info
    CHARACTER(len=*), INTENT(IN)             :: min_type
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    integer                                  :: iw

    iw = globenv % scr
    
    WRITE (iw,*)
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
    WRITE (iw,'(A,T71,I10)') ' ITERATION NUMBER        = ', iter
    WRITE (iw,'(A,T66,E15.7)') ' MAX FORCE ON COEFFICENTS = ', min_info%tol
    IF ( min_type == 'DIIS' ) &
    WRITE (iw,'(A,T66,E15.7)') ' ODIIS ERROR              = ', min_info%err
    WRITE (iw,'(A,T66,E15.7)') ' TOTAL HARTREE ENERGY  = ', min_info%e_hartree 
    WRITE (iw,'(A,T66,E15.7)') ' KINETIC+XC ENERGY  = ', min_info%e_xc
    WRITE (iw,'(A,T66,E15.7)') ' PPL        ENERGY  = ', min_info%e_ppl
    WRITE (iw,'(A,T66,E15.7)') ' PPL  LR    ENERGY  = ', min_info%e_ppl2
    WRITE (iw,'(A,T66,E15.7)') ' TOTAL ENERGY       = ', min_info%e
    WRITE (iw,'(A,T66,E15.7)') ' DELTA E           = ', min_info%de
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
!
  END SUBROUTINE print_energy_info
!
!******************************************************************************

  SUBROUTINE print_convergence ( conv_flag, globenv )

    LOGICAL, INTENT(IN)                      :: conv_flag
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    integer                                  :: iw

    iw = globenv % scr

     IF (conv_flag) THEN
       WRITE (iw,*)
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,'(T23,A)') '*    CONVERGENCE ACHIEVED    *'
       WRITE (iw,'(T23,A)') '*    IN COEF MINIMIZATION    *'
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,*)
     ELSE
       WRITE (iw,*)
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,'(T23,A)') '*  WARNING: NO CONVERGENCE!  *'
       WRITE (iw,'(T23,A)') '*  RESTART WITH LAST COEF    *'
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,*)
     END IF

  END SUBROUTINE print_convergence 

!------------------------------------------------------------------------------!
  SUBROUTINE DIIS ( c_pos, c_force, dtw, fconv, hess, &
                   e_hist, c_hist, nhist, ihist)
!-----
!  c_pos: basis functions coefficents
!  c_force: forces on the basis functions coefficents
!  dtw: is the time-step, dt
!  fconv: maximum gradient       
!  hess: is the Hessian  
!  e_hist: is a history vector
!  c_hist: is the history vector of the coefficients 
!  nhist: number of histories to be stored
!  ihist: history counter
!--------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: c_pos, c_force
    REAL(dbl), INTENT(IN)                    :: dtw
    REAL(dbl), INTENT(INOUT)                 :: fconv
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: hess
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: e_hist, c_hist
    INTEGER, INTENT(IN)                      :: nhist
    INTEGER, INTENT(INOUT)                   :: ihist

    INTEGER                                  :: icoef, j, ncoef
    REAL(dbl), ALLOCATABLE                   :: dcg(:)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: cg, eg

! allocating the arrays

   ncoef = size(c_pos)
   ALLOCATE(cg(ncoef)) 
   ALLOCATE(eg(ncoef))
   ALLOCATE(dcg(nhist+1))

! increment the counters for the histories and the iterations
  ihist=ihist+1

! forming the trial vector from the natural force 
  DO icoef = 1, ncoef
    cg (icoef) = c_pos (icoef) + dtw/hess(icoef) * c_force (icoef)
  END DO

! store the histories
  CALL store_hist(c_force, dtw, cg, e_hist, c_hist, hess, nhist, ihist)

! solve for the coefficients
  CALL solve(dcg,e_hist,ihist,ncoef)

! the new coefficients!
  DO icoef = 1, ncoef
    c_pos (icoef) = 0._dbl
  END DO
  eg(:)=0._dbl
  DO j=1, ihist
    DO icoef = 1, ncoef
      c_pos (icoef) = c_pos (icoef) + dcg(j) * c_hist(j,icoef)
      eg (icoef) = eg(icoef) + dcg(j) * e_hist(j,icoef)
    ENDDO
  ENDDO
  fconv=maxval(abs(eg))
! deallocating the arrays
  DEALLOCATE(cg)
  DEALLOCATE(eg)
  DEALLOCATE(dcg)
  return
  END SUBROUTINE diis
!------------------------------------------------------------------------------!
  SUBROUTINE STORE_HIST(c_forces,dtw,cg,e_hist,c_hist,hess,nhist,ihist)
!------------------------------------------------------------------------------!

    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: c_forces
    REAL(dbl), INTENT(IN)                    :: dtw
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: cg
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: e_hist, c_hist
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: hess
    INTEGER, INTENT(IN)                      :: nhist
    INTEGER, INTENT(INOUT)                   :: ihist

    INTEGER                                  :: icoef, j, ncoef

  ncoef = size(c_forces)

! constructing the error vector and the coefficient vector over histories
  IF (nhist==1) THEN
    ihist=nhist
    goto 2000
  ELSE
    IF (ihist>nhist) THEN
      ihist=nhist
      DO icoef=1,ncoef
        DO j=2, nhist
          e_hist(j-1,icoef)=e_hist(j,icoef)
          c_hist(j-1,icoef)=c_hist(j,icoef)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
2000  continue
  DO icoef = 1, ncoef
    c_hist (ihist,icoef) = cg (icoef)
    e_hist (ihist,icoef) = dtw * c_forces (icoef) / hess(icoef)
  ENDDO
  return
  END SUBROUTINE store_hist
!------------------------------------------------------------------------------!
  SUBROUTINE SOLVE(dcg,e_hist,ihist,ncoef)
!------------------------------------------------------------------------------!
    REAL(dbl), INTENT(OUT)                   :: dcg(:)
    REAL(dbl), DIMENSION(:, :), intent(IN)   :: e_hist
    INTEGER, INTENT(IN)                      :: IHIST, ncoef

    INTEGER                                  :: icoef, ios, j, k, mysize
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: bcg, eigenvec

! constructing the solution vector

  do j=1, ihist
   dcg(j)=0._dbl
  enddo
  dcg(ihist+1)=-1._dbl

! constructing the overlap matrix
  ALLOCATE (bcg(ihist+1,ihist+1), STAT = ios)
  IF ( ios /= 0 )  CALL stop_memory ( 'diis', 'bcg', ihist +1 )
  ALLOCATE (eigenvec(ihist+1,1), STAT = ios)
  IF ( ios /= 0 )  CALL stop_memory ( 'diis', 'eigenvec', ihist +1 )

  bcg=0._dbl
  do j=1,ihist
    do k=1,ihist
      do icoef = 1, ncoef
        bcg(j,k) = bcg(j,k) + e_hist(j,icoef) * e_hist(k,icoef)
      enddo
    enddo
  enddo
  do j=1, ihist
    bcg(ihist+1,j) = -1._dbl
    bcg(j,ihist+1) = -1._dbl
  enddo
  do j=1, ihist + 1
    eigenvec(j,1) = dcg(j)
  end do

  bcg (ihist+1,ihist+1) = 0._dbl
  mysize = ihist + 1

! get the solution vector

  call solve_system ( bcg, mysize, eigenvec)
 
  do j=1, ihist
     dcg(j) =  eigenvec(j,1) 
  end do

  DEALLOCATE(bcg)
  IF ( ios /= 0 )  CALL stop_memory ( 'diis', 'deall bcg' )
  DEALLOCATE(eigenvec)
  IF ( ios /= 0 )  CALL stop_memory ( 'diis', 'deall eigenvec' )
  RETURN
  END SUBROUTINE solve
!******************************************************************************

END MODULE kg_energy_utils

!******************************************************************************
