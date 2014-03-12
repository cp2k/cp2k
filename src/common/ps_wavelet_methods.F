!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition and initialisation of the ps_wavelet data type.
!> \history 01.2014 Renamed from ps_wavelet_types to disentangle dependencies (Ole Schuett)
!> \author Florian Schiffmann (09.2007,fschiff)
! *****************************************************************************
MODULE ps_wavelet_methods

  USE bibliography,                    ONLY: Genovese2006,&
                                             Genovese2007,&
                                             cite_reference
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_alltoall,&
                                             mp_cart_rank
  USE ps_wavelet_kernel,               ONLY: createKernel
  USE ps_wavelet_types,                ONLY: WAVELET0D,&
                                             ps_wavelet_type
  USE ps_wavelet_util,                 ONLY: F_FFT_dimensions,&
                                             PSolver,&
                                             P_FFT_dimensions,&
                                             S_FFT_dimensions
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_poisson_types,                ONLY: pw_poisson_parameter_type
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_methods'

! *** Public data types ***

  PUBLIC :: ps_wavelet_create,&
            RS_z_slice_distribution,&
            cp2k_distribution_to_z_slices,&
            z_slices_to_cp2k_distribution,&
            ps_wavelet_solve

CONTAINS

! *****************************************************************************
!> \brief creates the ps_wavelet_type which is needed for the link to
!>      the Poisson Solver of Luigi Genovese
!> \param wavlet wavelet to create
!> \param pw_grid the grid that is used to create the wavelet kernel
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Flroian Schiffmann
! *****************************************************************************
  SUBROUTINE ps_wavelet_create(poisson_params,wavelet,pw_grid, error)
    TYPE(pw_poisson_parameter_type), &
      INTENT(IN)                             :: poisson_params
    TYPE(ps_wavelet_type), POINTER           :: wavelet
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ps_wavelet_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iproc, nproc, nx, ny, &
                                                nz, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: hx, hy, hz

    CALL timeset(routineN,handle)

    CALL cite_reference(Genovese2006)
    CALL cite_reference(Genovese2007)

    failure=.FALSE.

    ALLOCATE(wavelet, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    nx = pw_grid % npts ( 1 )
    ny = pw_grid % npts ( 2 )
    nz = pw_grid % npts ( 3 )

    hx = pw_grid % dr ( 1 )
    hy = pw_grid % dr ( 2 )
    hz = pw_grid % dr ( 3 )

    nproc=PRODUCT(pw_grid % para % rs_dims)

    iproc=pw_grid % para % rs_mpo

    NULLIFY(wavelet%karray,wavelet%rho_z_sliced)

    wavelet%geocode           = poisson_params%wavelet_geocode
    wavelet%method            = poisson_params%wavelet_method
    wavelet%special_dimension = poisson_params%wavelet_special_dimension
    wavelet%itype_scf         = poisson_params%wavelet_scf_type
    wavelet%datacode          = "D"

    IF (poisson_params%wavelet_method==WAVELET0D) THEN
       IF (hx.NE.hy) &
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Poisson solver for non cubic cells not yet implemented",&
                error=error,failure=failure)
        IF(hz.NE.hy) &
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                "Poisson solver for non cubic cells not yet implemented",&
                 error=error,failure=failure)
    END IF

    CALL RS_z_slice_distribution(wavelet,pw_grid, error)

    CALL timestop(handle)
  END SUBROUTINE ps_wavelet_create

! *****************************************************************************
  SUBROUTINE RS_z_slice_distribution(wavelet,pw_grid, error)

    TYPE(ps_wavelet_type), POINTER           :: wavelet
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'RS_z_slice_distribution', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: geocode
    INTEGER                                  :: handle, iproc, m1, m2, m3, &
                                                md1, md2, md3, n1, n2, n3, &
                                                nd1, nd2, nd3, nproc, nx, ny, &
                                                nz, stat, z_dim
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: hx, hy, hz

    failure=.FALSE.

    CALL timeset(routineN,handle)
    nproc=PRODUCT(pw_grid % para % rs_dims)
    iproc= pw_grid % para % rs_mpo
    geocode=wavelet%geocode
    nx = pw_grid % npts ( 1 )
    ny = pw_grid % npts ( 2 )
    nz = pw_grid % npts ( 3 )
    hx = pw_grid % dr ( 1 )
    hy = pw_grid % dr ( 2 )
    hz = pw_grid % dr ( 3 )

    !calculate Dimensions for the z-distributed density and for the kernel

    IF (geocode == 'P') THEN
       CALL P_FFT_dimensions(nx,ny,nz,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
    ELSE IF (geocode == 'S') THEN
       CALL S_FFT_dimensions(nx,ny,nz,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
    ELSE IF (geocode == 'F') THEN
       CALL F_FFT_dimensions(nx,ny,nz,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
    END IF

    wavelet%PS_grid(1)=md1
    wavelet%PS_grid(2)=md3
    wavelet%PS_grid(3)=md2
    z_dim=md2/nproc
    !!!!!!!!!      indicies y and z are interchanged    !!!!!!!
    ALLOCATE(wavelet%rho_z_sliced(md1,md3,z_dim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL createKernel(geocode,nx,ny,nz,hx,hy,hz,wavelet%itype_scf,iproc,nproc,wavelet%karray,&
                       pw_grid % para % rs_group ,error)

    CALL timestop(handle)
  END SUBROUTINE RS_z_slice_distribution

! *****************************************************************************
  SUBROUTINE cp2k_distribution_to_z_slices (density , wavelet, pw_grid, error)

    TYPE(pw_type), POINTER                   :: density
    TYPE(ps_wavelet_type), POINTER           :: wavelet
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp2k_distribution_to_z_slices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dest, handle, i, ierr, ii, &
                                                iproc, j, k, l, local_z_dim, &
                                                loz, m, m2, md2, nproc
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, &
                                                sdispl, tmp
    INTEGER, DIMENSION(2)                    :: cart_pos, lox, loy
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: max_val_low, max_val_up
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rbuf, sbuf

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(wavelet),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN

       nproc=PRODUCT(pw_grid % para % rs_dims)
       iproc=pw_grid % para % rs_mpo
       md2=wavelet%PS_grid(3)
       m2= pw_grid % npts ( 3 )
       lb(:)=pw_grid % bounds_local ( 1, : )
       ub(:)=pw_grid % bounds_local ( 2, : )
       local_z_dim=MAX((md2/nproc),1)

       ALLOCATE(sbuf(PRODUCT(pw_grid % npts_local)),stat=ierr)
       ALLOCATE(rbuf(PRODUCT(wavelet%PS_grid)/nproc))
       ALLOCATE(scount(nproc),sdispl(nproc),rcount(nproc),rdispl(nproc),tmp(nproc))

       rbuf=0.0_dp
       ii=1
       DO k=lb(3),ub(3)
          DO j=lb(2),ub(2)
             DO i=lb(1),ub(1)
                sbuf(ii)=density % cr3d(i,j,k)
                ii=ii+1
             END DO
          END DO
       END DO


       IF(wavelet%geocode=='S'.OR.wavelet%geocode=='F')THEN
          max_val_low=0._dp
          max_val_up=0._dp
          IF(lb(2)==pw_grid%bounds(1,2))max_val_low=MAXVAL(ABS(density % cr3d(:,lb(2),:)))
          IF(ub(2)==pw_grid%bounds(2,2))max_val_up=MAXVAL(ABS(density % cr3d(:,ub(2),:)))
          CALL cp_assert((max_val_low.LT.0.0001_dp),cp_warning_level,cp_assertion_failed,routineP,&
               "Density hits the lower boundary of the system in XZ-plane"//&
               CPSourceFileRef,&
               only_ionode=.TRUE.)
          CALL cp_assert((max_val_up.LT.0.0001_dp),cp_warning_level,cp_assertion_failed,routineP,&
               "Density hits the upper boundary of the system in XZ-plane"//&
               CPSourceFileRef,&
               only_ionode=.TRUE.)
          IF(wavelet%geocode=='F')THEN
             max_val_low=0._dp
             max_val_up=0._dp
             IF(lb(1)==pw_grid%bounds(1,1))max_val_low=MAXVAL(ABS(density % cr3d(lb(1),:,:)))
             IF(ub(1)==pw_grid%bounds(2,1))max_val_up=MAXVAL(ABS(density % cr3d(ub(1),:,:)))
             CALL cp_assert((max_val_low.LT.0.0001_dp),cp_warning_level,cp_assertion_failed,routineP,&
                  "Density hits the lower boundary of the system in YZ-plane"//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
             CALL cp_assert((max_val_up.LT.0.0001_dp),cp_warning_level,cp_assertion_failed,routineP,&
                  "Density hits the upper boundary of the system in YZ-plane"//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
             max_val_low=0._dp
             max_val_up=0._dp
             IF(lb(3)==pw_grid%bounds(1,3))max_val_low=MAXVAL(ABS(density % cr3d(:,:,lb(3))))
             IF(ub(3)==pw_grid%bounds(2,3))max_val_up=MAXVAL(ABS(density % cr3d(:,:,ub(3))))
             CALL cp_assert((max_val_low.LT.0.0001_dp),cp_warning_level,cp_assertion_failed,routineP,&
                  "Density hits the lower boundary of the system in XY-plane"//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
             CALL cp_assert((max_val_up.LT.0.0001_dp),cp_warning_level,cp_assertion_failed,routineP,&
                  "Density hits the upper boundary of the system in XY-plane"//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
       END IF

       DO i = 0,pw_grid % para % rs_dims(1)-1
          DO j= 0, pw_grid % para % rs_dims(2)-1
             cart_pos=(/i,j/)
             CALL mp_cart_rank ( pw_grid % para % rs_group, &
                                 cart_pos, &
                                 dest )
             IF((ub(1).GE.lb(1)).AND.(ub(2).GE.lb(2)))THEN
                IF(dest*local_z_dim.LE.m2)THEN
                   IF((dest+1)*local_z_dim.LE.m2)THEN
                      scount(dest+1)=ABS((ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*local_z_dim)
                   ELSE
                      scount(dest+1)=ABS((ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*MOD(m2,local_z_dim))
                   END IF
                ELSE
                   scount(dest+1)=0
                END IF
             ELSE
                scount(dest+1)=0
             END IF
             lox = get_limit ( pw_grid % npts ( 1 ), pw_grid % para % rs_dims ( 1 ), i )
             loy = get_limit ( pw_grid % npts ( 2 ), pw_grid % para % rs_dims ( 2 ), j )
             IF((lox(2).GE.lox(1)).AND.(loy(2).GE.loy(1)))THEN
                IF(iproc*local_z_dim.LE.m2)THEN
                   IF((iproc+1)*local_z_dim.LE.m2)THEN
                      rcount(dest+1)=ABS((lox(2)-lox(1)+1)*(loy(2)-loy(1)+1)*local_z_dim)
                   ELSE
                      rcount(dest+1)=ABS((lox(2)-lox(1)+1)*(loy(2)-loy(1)+1)*MOD(m2,local_z_dim))
                   END IF
                ELSE
                   rcount(dest+1)=0
                END IF
             ELSE
                rcount(dest+1)=0
             END IF

          END DO
       END DO
       sdispl(1)=0
       rdispl(1)=0
       DO i = 2,nproc
          sdispl(i)=sdispl(i-1)+scount(i-1)
          rdispl(i)=rdispl(i-1)+rcount(i-1)
       END DO
       CALL mp_alltoall(sbuf,scount,sdispl,rbuf,rcount,rdispl,pw_grid%para%rs_group)
       !!!! and now, how to put the right cubes to the right position!!!!!!

       wavelet%rho_z_sliced=0.0_dp

       DO i = 0,pw_grid % para % rs_dims(1)-1
          DO j= 0, pw_grid % para % rs_dims(2)-1
             cart_pos=(/i,j/)
             CALL mp_cart_rank ( pw_grid % para % rs_group, &
                                 cart_pos, &
                                 dest )

             lox = get_limit ( pw_grid % npts ( 1 ), pw_grid % para % rs_dims ( 1 ), i )
             loy = get_limit ( pw_grid % npts ( 2 ), pw_grid % para % rs_dims ( 2 ), j )
             IF(iproc*local_z_dim.LE.m2)THEN
                IF((iproc+1)*local_z_dim.LE.m2)THEN
                   loz=local_z_dim
                ELSE
                   loz=MOD(m2,local_z_dim)
                END IF
                ii=1
                DO k=1,loz
                   DO l = loy(1),loy(2)
                      DO m = lox(1),lox(2)
                         wavelet%rho_z_sliced(m,l,k)=rbuf(ii+rdispl(dest+1))
                         ii=ii+1
                      END DO
                   END DO
                END DO
             END IF
          END DO
       END DO

       DEALLOCATE(sbuf,rbuf,scount,sdispl,rcount,rdispl,tmp)

    END IF
    CALL timestop(handle)

  END SUBROUTINE cp2k_distribution_to_z_slices

! *****************************************************************************
  SUBROUTINE z_slices_to_cp2k_distribution(density , wavelet, pw_grid, error)

    TYPE(pw_type), POINTER                   :: density
    TYPE(ps_wavelet_type), POINTER           :: wavelet
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'z_slices_to_cp2k_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dest, i, ierr, ii, iproc, j, &
                                                k, l, local_z_dim, loz, m, &
                                                m2, md2, nproc
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, &
                                                sdispl, tmp
    INTEGER, DIMENSION(2)                    :: cart_pos, lox, loy, min_x, &
                                                min_y
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rbuf, sbuf

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(wavelet),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       nproc=PRODUCT(pw_grid % para % rs_dims)
       iproc=pw_grid % para % rs_mpo
       md2=wavelet%PS_grid(3)
       m2= pw_grid % npts ( 3 )

       lb(:)=pw_grid % bounds_local ( 1, : )
       ub(:)=pw_grid % bounds_local ( 2, : )

       local_z_dim=MAX((md2/nproc),1)

       ALLOCATE(rbuf(PRODUCT(pw_grid % npts_local)),stat=ierr)
       ALLOCATE(sbuf(PRODUCT(wavelet%PS_grid)/nproc))
       ALLOCATE(scount(nproc),sdispl(nproc),rcount(nproc),rdispl(nproc),tmp(nproc))
       scount=0
       rcount=0
       rbuf=0.0_dp
       ii=1
       IF(iproc*local_z_dim.LE.m2)THEN
          IF((iproc+1)*local_z_dim.LE.m2)THEN
             loz=local_z_dim
          ELSE
             loz=MOD(m2,local_z_dim)
          END IF
       ELSE
          loz=0
       END IF

       min_x=get_limit ( pw_grid % npts ( 1 ), pw_grid % para % rs_dims ( 1 ), 0 )
       min_y=get_limit ( pw_grid % npts ( 2 ) , pw_grid % para % rs_dims ( 2 ), 0 )
       DO i = 0,pw_grid % para % rs_dims(1)-1
          DO j= 0, pw_grid % para % rs_dims(2)-1
             cart_pos=(/i,j/)
             CALL mp_cart_rank ( pw_grid % para % rs_group, &
                                 cart_pos, &
                                 dest )
             IF((ub(1).GE.lb(1)).AND.(ub(2).GE.lb(2)))THEN
                IF(dest*local_z_dim.LE.m2)THEN
                   IF((dest+1)*local_z_dim.LE.m2)THEN
                      rcount(dest+1)=ABS((ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*local_z_dim)
                   ELSE
                      rcount(dest+1)=ABS((ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*MOD(m2,local_z_dim))
                   END IF
                ELSE
                   rcount(dest+1)=0
                END IF
             ELSE
                rcount(dest+1)=0
             END IF
             lox = get_limit ( pw_grid % npts ( 1 ), pw_grid % para % rs_dims ( 1 ), i )
             loy = get_limit ( pw_grid % npts ( 2 ) , pw_grid % para % rs_dims ( 2 ), j )
             IF((lox(2).GE.lox(1)).AND.(loy(2).GE.loy(1)))THEN
                scount(dest+1)=ABS((lox(2)-lox(1)+1)*(loy(2)-loy(1)+1)*loz)
                DO k=lox(1)-min_x(1)+1,lox(2)-min_x(1)+1
                   DO l=loy(1)-min_y(1)+1,loy(2)-min_y(1)+1
                      DO m=1,loz
                         sbuf(ii)= wavelet%rho_z_sliced(k,l,m)
                         ii=ii+1
                      END DO
                   END DO
                END DO
             ELSE
                scount(dest+1)=0
             END IF
          END DO
       END DO
       sdispl(1)=0
       rdispl(1)=0
       DO i = 2,nproc
          sdispl(i)=sdispl(i-1)+scount(i-1)
          rdispl(i)=rdispl(i-1)+rcount(i-1)
       END DO
       CALL mp_alltoall(sbuf,scount,sdispl,rbuf,rcount,rdispl,pw_grid % para % rs_group)

       !!!! and now, how to put the right cubes to the right position!!!!!!

       DO i = 0,pw_grid % para % rs_dims(1)-1
          DO j= 0, pw_grid % para % rs_dims(2)-1
             cart_pos=(/i,j/)
             CALL mp_cart_rank ( pw_grid % para % rs_group, &
                                 cart_pos, &
                                 dest )
             IF(dest*local_z_dim.LE.m2)THEN
                IF((dest+1)*local_z_dim.LE.m2)THEN
                   loz=local_z_dim
                ELSE
                   loz=MOD(m2,local_z_dim)
                END IF
                ii=1
                IF(lb(3)+(dest*local_z_dim).LE.ub(3))THEN
                   DO m = lb(1),ub(1)
                      DO l = lb(2),ub(2)
                         DO k= lb(3)+(dest*local_z_dim),lb(3)+(dest*local_z_dim)+loz-1
                            density%cr3d(m,l,k)=rbuf(ii+rdispl(dest+1))
                            ii=ii+1
                         END DO
                      END DO
                   END DO
                END IF
             END IF
          END DO
       END DO
       DEALLOCATE(sbuf,rbuf,scount,sdispl,rcount,rdispl,tmp)
    END IF

  END SUBROUTINE z_slices_to_cp2k_distribution

! *****************************************************************************
    SUBROUTINE ps_wavelet_solve(wavelet,pw_grid,eh,error)

    TYPE(ps_wavelet_type), POINTER           :: wavelet
    TYPE(pw_grid_type), POINTER              :: pw_grid
    REAL(KIND=dp)                            :: eh
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ps_wavelet_solve', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: datacode, geocode
    INTEGER                                  :: handle, iproc, nproc, nx, ny, &
                                                nz
    REAL(KIND=dp)                            :: hx, hy, hz

    CALL timeset(routineN,handle)
    nproc=PRODUCT(pw_grid % para % rs_dims)
    iproc= pw_grid % para % rs_mpo
    geocode=wavelet%geocode
    datacode=wavelet%datacode
    nx = pw_grid % npts ( 1 )
    ny = pw_grid % npts ( 2 )
    nz = pw_grid % npts ( 3 )
    hx = pw_grid % dr ( 1 )
    hy = pw_grid % dr ( 2 )
    hz = pw_grid % dr ( 3 )

    CALL PSolver(geocode,datacode,iproc,nproc,nx,ny,nz,hx,hy,hz,&
         wavelet%rho_z_sliced,wavelet%karray,pw_grid,error)
    CALL timestop(handle)
  END SUBROUTINE ps_wavelet_solve

END MODULE ps_wavelet_methods
