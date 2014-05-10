!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that collects GA functionality for semi-empirical
!>
!> \author CJM and Bruce Palmer (05.2011)
!> \par History
!>
! *****************************************************************************
MODULE se_ga_tools
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ga_environment_types,            ONLY: ga_environment_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_ga_tools'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.


! *****************************************************************************

! Public subroutine

  PUBLIC :: se_deallocate_local_ga, se_ga_ks_accumulate, &
            se_allocate_local_ga, se_ga_initialize, &
            se_ga_allocate_local_info, se_ga_deallocate_local_info, &
            se_ga_get_nbin, se_ga_release, se_ga_diag_add, &
            se_ga_put_pmatrix, se_ga_pair_list_init

  INTERFACE se_ga_allocate_local_info
     MODULE PROCEDURE se_ga_allocate_local_2, se_ga_allocate_local_1
  END INTERFACE

  INTERFACE se_ga_deallocate_local_info
     MODULE PROCEDURE se_ga_deallocate_local_2, se_ga_deallocate_local_1
  END INTERFACE

  INTERFACE se_ga_ks_accumulate
     MODULE PROCEDURE se_ga_ks_accumulate_1, se_ga_ks_accumulate_2
  END INTERFACE


#ifdef _USE_GA
#include "mafdecls.fh"
#include "global.fh"
#endif

CONTAINS
! *****************************************************************************
!> \brief Module that collects GA functionality
!>
!> \author CJM and Bruce Palmer (05.2011)
!> \par History
!>
! *****************************************************************************
  SUBROUTINE se_ga_release ( qs_env, error )
    IMPLICIT NONE
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
   INTEGER :: stat
   TYPE ( ga_environment_type ), POINTER :: ga_env
   LOGICAL :: status, failure
   CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_release', &
     routineP = moduleN//':'//routineN

   CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
   NULLIFY ( ga_env )

   CALL get_qs_env(qs_env, ga_env=ga_env, error=error)

   CPPrecondition(ASSOCIATED(ga_env),cp_failure_level,routineP,error,failure)

! write timing
!   WRITE ( *, * ) ga_nodeid (), ' TIME_PUT', ga_env % time_put, ga_env % iput
!   WRITE ( *, * ) ga_nodeid (), ' TIME_GET', ga_env % time_get, ga_env % iget
!   WRITE ( *, * ) ga_nodeid (), ' TIME_ACC', ga_env % time_acc, ga_env % iacc
!
!   call ga_igop(1, ga_env % iput, 1, '+')
!   call ga_igop(1, ga_env % iget, 1, '+')
!   call ga_igop(1, ga_env % iacc, 1, '+')
!   WRITE ( *,'(i4,a,i16)') ga_nodeid(),' Total puts: ', ga_env % iput
!   WRITE ( *,'(i4,a,i16)') ga_nodeid(),' Total gets: ', ga_env % iget
!   WRITE ( *,'(i4,a,i16)') ga_nodeid(),' Total accs: ', ga_env % iacc

!   call ga_igop(1, ga_env % ielm_get, 1, '+')
!   call ga_igop(1, ga_env % ielm_put, 1, '+')
!   call ga_igop(1, ga_env % ielm_acc, 1, '+')
!   WRITE ( *,'(i4,a,i16)') ga_nodeid(),' Total elm puts: ', ga_env % ielm_put
!   WRITE ( *,'(i4,a,i16)') ga_nodeid(),' Total elm gets: ', ga_env % ielm_get
!   WRITE ( *,'(i4,a,i16)') ga_nodeid(),' Total elm accs: ', ga_env % ielm_acc

   status=ga_destroy( ga_env % g_offset )
   status=ga_destroy( ga_env % g_ks )
   status=ga_destroy( ga_env % g_atoms )
   status=ga_destroy( ga_env % g_d_info )
   status=ga_destroy( ga_env % g_cellnum )
   status=ga_destroy( ga_env % g_count )
   status=ga_destroy( ga_env % g_cnt )

#endif
  END SUBROUTINE se_ga_release
! *****************************************************************************
!> \brief Specific GA tool interfaces for the semi-empirical code
!> \author CJM
! *****************************************************************************
  SUBROUTINE se_ga_pair_list_init ( qs_env, error )
    IMPLICIT NONE
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error
#ifdef _USE_GA
    TYPE ( cp_para_env_type ), POINTER      :: para_env
    TYPE ( ga_environment_type ), POINTER   :: ga_env
    TYPE ( cell_type ), POINTER             :: cell
    TYPE ( dft_control_type ), POINTER      :: dft_control
    TYPE ( semi_empirical_control_type ), POINTER      :: se_control
    REAL ( dp ) ::  xbox, ybox, zbox
    REAL ( dp ) ::  xcell, ycell, zcell, cx, cy, cz, rc, rcut2
    INTEGER ::  kxmin,kymin,kzmin,kxmax,kymax,kzmax,kc,kx,ky,kz
    INTEGER ::  nbins, nc, ix, iy, iz, itmp
    INTEGER ::  jx, jy, jz, idx
    INTEGER ::  ierr, ibeg, iend, me, nproc, i, j, icnt
    INTEGER ::  g_cnt, one
    INTEGER ::  lo(2), hi(2), dims(2),chunk(2), ndim, ld
    INTEGER ::  pair(2)
    INTEGER ::  pdx, pdy, pdz, number
    INTEGER ::  kcmax, nix(2000), niy(2000), niz(2000)
    LOGICAL ::  status, wrap

    NULLIFY ( cell, ga_env, dft_control, se_control )

    CALL get_qs_env(qs_env, dft_control=dft_control, &
                    para_env=para_env,ga_env=ga_env, &
                    cell = cell, error=error )


    se_control => dft_control % qs_control % se_control
    ga_env % rcoul = 2.0_dp * se_control % cutoff_cou
    ga_env % ncells = se_control % ga_ncells

    WRITE ( *, * ) 'AAAA', ga_env % ncells, ga_env % rcoul

    xbox=cell%hmat(1,1)
    ybox=cell%hmat(2,2)
    zbox=cell%hmat(3,3)

    IF ( ga_env % ncells == 0 ) THEN
! use cut-off
      pdx=INT ( xbox/ga_env%rcoul )
      IF ( ( xbox/REAL ( pdx, dp ) < ga_env % rcoul ) .AND. ( pdx > 1 ) ) THEN
        pdx = pdx - 1
      ELSE IF ( pdx == 0 ) THEN
        pdx = 1
      END IF
      pdy=INT ( ybox/ga_env%rcoul )
      IF ( ( ybox/REAL ( pdy, dp ) < ga_env % rcoul ) .AND. ( pdy > 1 ) ) THEN
        pdy = pdy - 1
      ELSE IF ( pdy == 0 ) THEN
        pdy = 1
      END IF
      pdz=INT ( zbox/ga_env%rcoul )
      IF ( ( zbox/REAL ( pdz, dp ) < ga_env % rcoul ) .AND. ( pdz > 1 ) ) THEN
        pdz = pdz - 1
      ELSE IF ( pdz == 0 ) THEN
        pdz = 1
      END IF
      ga_env % ncells = pdx*pdy*pdz
      WRITE ( *, * ) 'pds ', pdx,pdy,pdz
    ELSE
! use input driven number and compute pdx, pdy, pdz
      number = ga_env % ncells
      CALL grid_factor ( number, cell%hmat(1,1), cell%hmat(2,2), cell%hmat(3,3), pdx, pdy, pdz )
    END IF
! Assign value for use later
    ga_env % pdx = pdx
    ga_env % pdy = pdy
    ga_env % pdz = pdz
!
!  This  finds the neighboring cells
!
   xcell = xbox / REAL ( pdx, dp )
   kxmax = INT(ga_env % rcoul/xcell)+1
!
   ycell = ybox / REAL ( pdy, dp )
   kymax = INT(ga_env % rcoul/ycell)+1
!
   zcell = zbox / REAL ( pdz, dp )
   kzmax = INT(ga_env % rcoul/zcell)+1
!
   WRITE ( *, * ) 'ks ', kxmax, kymax, kzmax
!
  rcut2 = ga_env % rcoul**2
!
  kc = 1
  nix(kc) = 0
  niy(kc) = 0
  niz(kc) = 0
!
  kxmin = 0
  DO kx = kxmin, kxmax
    kymin = -kymax
    IF (kx==0) kymin = 0
    DO ky = kymin, kymax
      kzmin = -kzmax
      IF (kx==0.AND.ky==0) kzmin = 1
      DO kz = kzmin, kzmax
!        write(6,'(a,6i3)') 'klims: ',kxmin,kxmax,kymin,kymax,kzmin,kzmax
        IF (kx/=0) THEN
          cx = xcell*REAL(iabs(kx)-1, dp)
        ELSE
          cx = 0.0_dp
        ENDIF
        IF (ky/=0) THEN
          cy = ycell*REAL(iabs(ky)-1, dp)
        ELSE
          cy = 0.0_dp
        ENDIF
        IF (kz/=0) THEN
          cz = zcell*REAL(iabs(kz)-1, dp)
        ELSE
          cz = 0.0_dp
        ENDIF
        rc = cx**2+cy**2+cz**2
        IF (rc<rcut2) THEN
          kc = kc + 1
          nix(kc) = kx
          niy(kc) = ky
          niz(kc) = kz
        ENDIF
      END DO
    END DO
  END DO
  kcmax = kc
  WRITE ( *, * ) 'kcmax ', kcmax
! We now have unique neighbors of each cell

!
!  Initialize system
!
  me = ga_nodeid()
  nproc = ga_nnodes()
!
  rcut2 = ga_env % rcoul**2
  xcell = xbox/REAL(pdx, dp)
  ycell = ybox/REAL(pdy, dp)
  zcell = zbox/REAL(pdz, dp)
  WRITE ( *, * ) 'boxs ', xbox, ybox, zbox
  WRITE ( *, * ) 'cells ', xcell, ycell, zcell
!
!  Find total number of atom bins
!
  nbins = pdx*pdy*pdz
  WRITE ( *, * ) 'nbins ', nbins
  one = 1
!
!  Find number of unique neighbors
!
!  if (me.eq.0) write(6,'(a,i8)') 'Number of neighbors found: ',kcmax
!  if (me.eq.0) then
!    do i = 1, kcmax
!      write(6,'(i4,a,3i4)') i,' nix, niy, niz: ',nix(i),niy(i),niz(i)
!    end do
!  endif
!
!  Set up distributed pair list
!
  ibeg = INT(REAL(nbins,dp)*REAL(me,dp)/REAL(nproc,dp)) + 1
  IF (me.lt.nproc-1) THEN
    iend = INT(REAL(nbins,dp)*REAL(me+1,dp)/REAL(nproc,dp))
  ELSE
    iend = nbins
  ENDIF
!
!  Add all pairs to the distributed pair list
!
  icnt = 0
  DO i = ibeg, iend
    itmp = i - 1
    ix = MOD(itmp,pdx)
    itmp = (itmp-ix)/pdy
    iy = MOD(itmp,pdy)
    iz = (itmp-iy)/pdy
    ix = ix + 1
    iy = iy + 1
    iz = iz + 1
    DO nc = 1, kcmax
      wrap = .FALSE.
      jx = ix + nix(nc)
      IF (jx>pdx) THEN
        jx = MOD(jx,pdx)
        wrap = .TRUE.
      ENDIF
      IF (jx<1) THEN
        jx = pdx+jx
        wrap = .TRUE.
      ENDIF
      jy = iy + niy(nc)
      IF (jy>pdy) THEN
        jy = MOD(jy,pdy)
        wrap = .TRUE.
      ENDIF
      IF (jy<1) THEN
        jy = pdy+jy
        wrap = .TRUE.
      ENDIF
      jz = iz + niz(nc)
      IF (jz>pdz) THEN
        jz = MOD(jz,pdz)
        wrap = .TRUE.
      ENDIF
      IF (jz<1) THEN
        jz = pdz+jz
        wrap = .TRUE.
      ENDIF
      j = (jz-1)*pdx*pdy + (jy-1)*pdx + jx
      IF (wrap) THEN
        IF (iabs(ix-jx).ne.0) THEN
          cx = xcell*(iabs(ix-jx)-1)
        ELSE
          cx = 0.0_dp
        ENDIF
        IF (iabs(iy-jy)/=0) THEN
          cy = ycell*(iabs(iy-jy)-1)
        ELSE
          cy = 0.0_dp
        ENDIF
        IF (iabs(iz-jz)/=0) THEN
          cz = zcell*(iabs(iz-jz)-1)
        ELSE
          cz = 0.0_dp
        ENDIF
        rc = cx**2 + cy**2 + cz**2
        IF (rc<rcut2) CYCLE
      ENDIF
      icnt = icnt + 1
    END DO
  END DO
!
!  Calculated all pair interactions on this processor
!
! Hold on to this
  ga_env % npairs = icnt
  WRITE ( *, * ) 'NPAIRS', icnt
  CALL ga_igop(1,ga_env%npairs,1,'+')
!  if (me.eq.0) write(6,'(a,i8)') 'Total number of pairs found: ',npairs
!
!  Create global counter and distributed pair list
!
  ndim = 1
  dims(1) = 1
  g_cnt = ga_create_handle()
  CALL ga_set_data(g_cnt, ndim, dims, MT_INT)
  status = ga_allocate(g_cnt)
!
  ndim = 2
  dims(1) = 2
  dims(2) = ga_env % npairs
  chunk(1) = 2
  chunk(2) = -1
! hold on g_pair_list
  ga_env % g_pair_list = ga_create_handle()
  CALL ga_set_data(ga_env % g_pair_list, ndim, dims, MT_INT)
  CALL ga_set_chunk(ga_env % g_pair_list, chunk)
  status = ga_allocate(ga_env % g_pair_list)
!
  CALL ga_zero(g_cnt)
  lo(1) = 1
  hi(1) = 2
  ld = 2
  DO i = ibeg, iend
    itmp = i - 1
    ix = MOD(itmp,pdx)
    itmp = (itmp-ix)/pdy
    iy = MOD(itmp,pdy)
    iz = (itmp-iy)/pdy
    ix = ix + 1
    iy = iy + 1
    iz = iz + 1
    pair(1) = i
    DO nc = 1, kcmax
      wrap = .FALSE.
      jx = ix + nix(nc)
      IF (jx>pdx) THEN
        jx = MOD(jx,pdx)
        wrap = .TRUE.
      ENDIF
      IF (jx<1) THEN
        jx = pdx+jx
        wrap = .TRUE.
      ENDIF
      jy = iy + niy(nc)
      IF (jy>pdy) THEN
        jy = MOD(jy,pdy)
        wrap = .TRUE.
      ENDIF
      IF (jy<1) THEN
        jy = pdy+jy
        wrap = .TRUE.
      ENDIF
      jz = iz + niz(nc)
      IF (jz>pdz) THEN
        jz = MOD(jz,pdz)
        wrap = .TRUE.
      ENDIF
      IF (jz<1) THEN
        jz = pdz+jz
        wrap = .TRUE.
      ENDIF
      j = (jz-1)*pdx*pdy + (jy-1)*pdx + jx
      IF (wrap) THEN
        IF (iabs(ix-jx)/=0) THEN
          cx = xcell*(iabs(ix-jx)-1)
        ELSE
          cx = 0.0_dp
        ENDIF
        IF (iabs(iy-jy)/=0) THEN
          cy = ycell*(iabs(iy-jy)-1)
        ELSE
          cy = 0.0_dp
        ENDIF
        IF (iabs(iz-jz)/=0) THEN
          cz = zcell*(iabs(iz-jz)-1)
        ELSE
          cz = 0.0_dp
        ENDIF
        rc = cx**2 + cy**2 + cz**2
        IF (rc<rcut2) CYCLE
      ENDIF
      idx = nga_read_inc(g_cnt, one, one)
!
      lo(2) = idx+1
      hi(2) = idx+1
      pair(2) = j
      CALL nga_put(ga_env % g_pair_list,lo,hi,pair,ld)
    END DO
  END DO
  CALL ga_sync
!
!  Evaluate interactions between all pairs
!
  status = ga_destroy(g_cnt)
#endif
  END SUBROUTINE se_ga_pair_list_init
! *****************************************************************************
!> \brief Specific GA tool interfaces for the semi-empirical code
!> \author CJM
! *****************************************************************************
  SUBROUTINE se_ga_initialize ( qs_env, error )
    IMPLICIT NONE
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error
#ifdef _USE_GA
    TYPE ( cp_para_env_type ), POINTER      :: para_env
    TYPE ( ga_environment_type ), POINTER   :: ga_env
    TYPE ( distribution_1d_type ), POINTER  :: local_particles
    TYPE ( particle_type ), DIMENSION(:), &
      POINTER                                 :: particle_set
    TYPE ( cell_type ), POINTER               :: cell
    TYPE ( semi_empirical_type ), POINTER     :: se_kind
    TYPE ( atomic_kind_type ), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE ( atomic_kind_type ), POINTER       :: atomic_kind
    INTEGER :: nx, ny, nz, number,  xbin, ybin, zbin, nbin
    INTEGER :: nparticle_kind, iparticle_kind, nparticle_local, iparticle_local
    INTEGER :: iparticle, k, j, icnt, natoms, natorb
    REAL ( dp ) :: dx, dy, dz, r_pbc ( 3 )
    INTEGER :: stat, chunk ( 3 ), dims ( 3 )
    INTEGER :: handle, lo ( 2 ), hi ( 2 )
    INTEGER :: ioff, iloc
    INTEGER, PARAMETER :: one = 1
    INTEGER :: ii, jj
    INTEGER, SAVE :: max_rows
    LOGICAL, SAVE :: first_time=.TRUE.
    REAL ( dp ), POINTER :: block ( :, : )
    REAL ( dp ), DIMENSION ( : ), POINTER :: tmp_info
    LOGICAL :: failure, status, check, defined, found
    CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_intialize', &
      routineP = moduleN//':'//routineN

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    NULLIFY ( para_env, cell, ga_env )
    NULLIFY ( atomic_kind_set, particle_set, local_particles )
    NULLIFY ( tmp_info )

    CALL get_qs_env(qs_env, para_env=para_env,  &
                    particle_set=particle_set, &
                    atomic_kind_set=atomic_kind_set, &
                    local_particles=local_particles, ga_env=ga_env, &
                    cell = cell, error=error )


    natoms = SIZE ( particle_set )



    CPPrecondition ( ASSOCIATED ( ga_env ), cp_failure_level,routineP,error,failure)


    number = ga_env % ncells
    nx = ga_env % pdx
    ny = ga_env % pdy
    nz = ga_env % pdz

!    CALL grid_factor ( number, cell%hmat(1,1), cell%hmat(2,2), cell%hmat(3,3), nx, ny, nz )

    dx=cell%hmat(1,1)/REAL ( nx, dp )
    dy=cell%hmat(2,2)/REAL ( ny, dp )
    dz=cell%hmat(3,3)/REAL ( nz, dp )

    ga_env % g_cnt=ga_create_handle( )
    CALL ga_set_data ( ga_env % g_cnt, one, one, MT_F_INT )
    status=ga_allocate ( ga_env % g_cnt )
    CALL ga_zero ( ga_env % g_cnt )

    ga_env % g_cellnum=ga_create_handle( )
    CALL ga_set_data ( ga_env % g_cellnum, one, number, MT_F_INT )
    status=ga_allocate ( ga_env % g_cellnum )
    CALL ga_zero ( ga_env % g_cellnum )

    ga_env % g_count=ga_create_handle ( )
    CALL ga_set_data ( ga_env % g_count, one, number, MT_F_INT )
    status=ga_allocate ( ga_env % g_count )
    CALL ga_zero ( ga_env % g_count )

    ga_env % g_offset=ga_create_handle ( )
    CALL ga_set_data ( ga_env % g_offset, one, number, MT_F_INT )
    status=ga_allocate ( ga_env % g_offset )
    CALL ga_zero ( ga_env % g_offset )

    IF ( para_env%mepos == 0 ) &
    CALL nga_put ( ga_env % g_count, one, one, one, one )

    ga_env % g_atoms=ga_create_handle ( )
    CALL ga_set_data ( ga_env % g_atoms, 1, natoms, MT_F_INT )
    status=ga_allocate ( ga_env % g_atoms )
    CALL ga_zero ( ga_env % g_atoms )

    IF ( first_time ) max_rows = 0
    nparticle_kind = SIZE(atomic_kind_set)
    DO iparticle_kind=1,nparticle_kind
      atomic_kind => atomic_kind_set(iparticle_kind)
      CALL get_atomic_kind(atomic_kind=atomic_kind, se_parameter=se_kind )
      CALL get_se_param ( se_kind, defined=defined, natorb=natorb )
      IF (.NOT. defined .AND. natorb < 1) CYCLE
      nparticle_local = local_particles%n_el(iparticle_kind)
      DO iparticle_local=1,nparticle_local
        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
        IF ( first_time ) &
        max_rows = MAX ( natorb, max_rows )
        r_pbc ( : )  = pbc ( particle_set ( iparticle ) % r ( : ), cell )
! should be between 0 -> L
        CALL pbc_wrap ( r_pbc, cell%hmat(1,1), cell%hmat(2,2), cell%hmat (3,3) )
        xbin = INT ( r_pbc  ( 1 ) / dx )
        ybin = INT ( r_pbc  ( 2 ) / dy )
        zbin = INT ( r_pbc  ( 3 ) / dz )
        nbin = zbin * nx * ny + ybin* nx + xbin  +  1
        IF ( ( nbin > ga_env % ncells ) .OR. ( nbin < 1 ) ) THEN
           WRITE ( *, * ) "WARNING", iparticle, nbin, ga_env % ncells
        END IF
        CALL nga_acc ( ga_env % g_cellnum, nbin, nbin, one, one, one )
      END DO
    END DO

    IF ( first_time ) THEN
!      CALL mp_max ( max_rows, para_env%group )
      first_time = .FALSE.
    END IF
    !WRITE ( *, * ) 'MAXROWS', max_rows
    max_rows=4

! set se_dims
    ga_env % se_dims ( 1 ) = max_rows
    ga_env % se_dims ( 2 ) = max_rows
    ga_env % se_dims ( 3 ) = natoms

! compute offset array
    CALL ga_sync ( )
    CALL ga_scan_add(ga_env % g_cellnum, ga_env % g_offset, ga_env % g_count, 1, number, 1)
    CALL ga_zero ( ga_env % g_count )

    dims(1) = ga_env % se_dims(1)*ga_env % se_dims(1) + 5
    dims(2) = ga_env % se_dims(3)
    chunk(1) = dims(1)
    chunk(2) = -1
    ALLOCATE(tmp_info(dims(1)))
    tmp_info = 0.0_dp

    ga_env % g_d_info = ga_create_handle()
    CALL ga_set_data(ga_env % g_d_info,2,dims,MT_F_DBL)
    CALL ga_set_chunk(ga_env % g_d_info,chunk)
    status = ga_allocate(ga_env % g_d_info)

    dims(1) = ga_env % se_dims(1)
    dims(2) = ga_env % se_dims(1)
    dims(3) = ga_env % se_dims(3)

    chunk(1) = dims(1)
    chunk(2) = dims(1)
    chunk(3) = -1

    ga_env % g_ks = ga_create_handle()
    CALL ga_set_data(ga_env % g_ks,3,dims,MT_F_DBL)
    CALL ga_set_chunk(ga_env % g_ks,chunk)
    status = ga_allocate(ga_env % g_ks)


    nparticle_kind = SIZE(atomic_kind_set)
    DO iparticle_kind=1,nparticle_kind
      atomic_kind => atomic_kind_set(iparticle_kind)
      CALL get_atomic_kind(atomic_kind=atomic_kind, se_parameter=se_kind )
      CALL get_se_param ( se_kind, defined=defined, natorb=natorb )
     ! WRITE ( *, * ) 'NATORB', natorb
      nparticle_local = local_particles%n_el(iparticle_kind)
      DO iparticle_local=1,nparticle_local
        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
        r_pbc ( : )  = pbc ( particle_set ( iparticle ) % r ( : ), cell )
! should be between 0 -> L
        CALL pbc_wrap ( r_pbc, cell%hmat(1,1), cell%hmat(2,2), cell%hmat (3,3) )
        xbin = INT ( r_pbc  ( 1 ) / dx )
        ybin = INT ( r_pbc  ( 2 ) / dy )
        zbin = INT ( r_pbc  ( 3 ) / dz )
        nbin = zbin * nx * ny + ybin* nx + xbin  +  1
        iloc = nga_read_inc(ga_env % g_count, nbin, 1)
        CALL nga_get(ga_env % g_offset, nbin, nbin, ioff, 1)
        CALL nga_put(ga_env % g_atoms, iparticle, iparticle, ioff+iloc+1,1)

        tmp_info(1) = REAL ( iparticle_kind, dp )
        tmp_info(2) = REAL ( iparticle, dp )
        tmp_info(3) = particle_set ( iparticle ) % r ( 1 )
        tmp_info(4) = particle_set ( iparticle ) % r ( 2 )
        tmp_info(5) = particle_set ( iparticle ) % r ( 3 )
       ! WRITE ( *, * ) 'TMP ', iparticle, SNGL ( tmp_info ( 3:5 ))
        icnt = 0
        DO k = 1, natorb
          DO j = 1, natorb
            icnt = icnt + 1
            tmp_info(5+icnt) = 0.0_dp
          END DO
        END DO
        lo(1) = 1
        lo(2) = ioff+iloc+1
        hi(1) = ga_env%se_dims ( 1 )**2 + 5
        hi(2) = ioff+iloc+1
        CALL nga_put(ga_env % g_d_info, lo, hi, tmp_info, dims(1))
      END DO
    END DO
    DEALLOCATE(tmp_info)
    CALL ga_sync ()
    CALL timestop(handle)
#endif
  END SUBROUTINE se_ga_initialize
! *****************************************************************************
!> \brief Module that collects GA functionality
!>
!> \author CJM and Bruce Palmer (05.2011)
!> \par History
!>
! *****************************************************************************
  SUBROUTINE se_allocate_local_ga ( buff, max_rows, error )
    REAL ( dp ), POINTER :: buff ( :, : )
    INTEGER, INTENT ( IN ) :: max_rows
    TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
   INTEGER :: stat
   LOGICAL :: failure
   CHARACTER(len=*), PARAMETER :: routineN = 'se_allocate_local_ga', &
     routineP = moduleN//':'//routineN
   ALLOCATE ( buff ( max_rows, max_rows ), STAT=stat )
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
  END SUBROUTINE se_allocate_local_ga
! *****************************************************************************
!> \brief Module that collects GA functionality
!>
!> \author CJM and Bruce Palmer (05.2011)
!> \par History
!>
! *****************************************************************************
  SUBROUTINE se_deallocate_local_ga ( buff, error )
    REAL ( dp ), POINTER :: buff ( :, : )
    TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
   INTEGER :: stat
   LOGICAL :: failure, status
   CHARACTER(len=*), PARAMETER :: routineN = 'se_deallocate_local_ga', &
     routineP = moduleN//':'//routineN
   DEALLOCATE(buff,STAT=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   NULLIFY ( buff )
#endif
  END SUBROUTINE se_deallocate_local_ga
! *****************************************************************************
    SUBROUTINE se_ga_put_pmatrix ( matrix, ga_env, error )
    IMPLICIT NONE
    TYPE ( cp_dbcsr_type ), POINTER :: matrix
    TYPE ( ga_environment_type ), POINTER :: ga_env
    TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
   TYPE ( cp_dbcsr_iterator ) :: iter
   LOGICAL :: failure
   INTEGER :: handle, stat
   INTEGER :: iblock_row, iblock_col, rows_block
   INTEGER ::  lo ( 2 ), hi ( 2 ), ld
   INTEGER :: blk, idx
   INTEGER, PARAMETER :: one = 1
   REAL(kind=dp), DIMENSION(:, :), POINTER  ::   block
   CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_put_pmatrix', &
     routineP = moduleN//':'//routineN

   CALL timeset(routineN,handle)
   CALL ga_sync
   CALL cp_dbcsr_iterator_start (iter, matrix )
   DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
      CALL cp_dbcsr_iterator_next_block (iter, iblock_row,iblock_col,  block ,blk)
      IF ( iblock_row == iblock_col) THEN
         rows_block=SIZE(block,1)
         CALL nga_get ( ga_env % g_atoms, iblock_row, iblock_row, idx, one )
         lo(1) = 6
         lo(2) = idx
         hi(1) = 5+rows_block**2
         hi(2) = idx
         ld = 1
         CALL nga_put(ga_env % g_d_info, lo, hi, block, ld)
      END IF
   END DO
   CALL cp_dbcsr_iterator_stop ( iter )
   CALL timestop(handle)
#endif
  END SUBROUTINE se_ga_put_pmatrix
! *****************************************************************************
    SUBROUTINE se_ga_diag_add ( matrix, ga_env, error )
    IMPLICIT NONE
    TYPE ( cp_dbcsr_type ), POINTER :: matrix
    TYPE ( ga_environment_type ), POINTER :: ga_env
    TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
   TYPE ( cp_dbcsr_iterator ) :: iter
   LOGICAL :: failure
   INTEGER :: handle, stat
   INTEGER :: iblock_row, iblock_col, rows_block
   INTEGER ::  lo ( 3 ), hi ( 3 ), ld ( 2 )
   INTEGER :: blk, idx
   INTEGER, PARAMETER :: one = 1
   REAL(kind=dp), DIMENSION(:, :), POINTER  ::  dum, block
   CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_diag_add', &
     routineP = moduleN//':'//routineN

   CALL timeset(routineN,handle)
   CALL ga_sync
   CALL cp_dbcsr_iterator_start (iter, matrix )
   DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
      CALL cp_dbcsr_iterator_next_block (iter, iblock_row,iblock_col,  block ,blk)
      IF ( iblock_row == iblock_col) THEN
         rows_block=SIZE(block,1)
         CALL nga_get ( ga_env % g_atoms, iblock_row, iblock_row, idx, one )
         lo(1) = 1
         lo(2) = 1
         lo(3) = idx
         hi(1) = rows_block
         hi(2) = rows_block
         hi(3) = idx
         ld(1) = rows_block
         ld(2) = rows_block
         ALLOCATE ( dum ( rows_block, rows_block ),STAT=stat )
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         CALL nga_get(ga_env % g_ks, lo, hi, dum, ld)
         block = block + dum
         DEALLOCATE ( dum , STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
   END DO
   CALL cp_dbcsr_iterator_stop ( iter )
   CALL timestop(handle)
#endif
  END SUBROUTINE se_ga_diag_add
! *****************************************************************************
    SUBROUTINE se_ga_ks_accumulate_2 ( ga_env, ks_i, ks_j, ioff, joff, isize, jsize )
      TYPE ( ga_environment_type ), POINTER :: ga_env
      REAL ( dp ), DIMENSION ( :, :, : ), POINTER :: ks_i, ks_j
      INTEGER, INTENT ( IN ) :: ioff, joff, isize, jsize
#ifdef _USE_GA
      INTEGER :: ilo ( 3 ), jlo ( 3 ), ihi ( 3 ), jhi ( 3 )
      REAL ( dp ), PARAMETER :: one=1.0_dp
      ilo(1) = 1
      ilo(2) = 1
      ilo(3) = ioff + 1
      jlo(1) = 1
      jlo(2) = 1
      jlo(3) = joff + 1
      ihi(1) = ga_env % se_dims(1)
      ihi(2) = ga_env % se_dims(2)
      ihi(3) = ioff + isize
      jhi(1) = ga_env % se_dims(1)
      jhi(2) = ga_env % se_dims(2)
      jhi(3) = joff + jsize
      CALL nga_acc ( ga_env % g_ks, ilo, ihi, ks_i, ga_env % se_dims ( 1 ), one )
      CALL nga_acc ( ga_env % g_ks, jlo, jhi, ks_j, ga_env % se_dims ( 1 ), one )
#endif
    END SUBROUTINE se_ga_ks_accumulate_2
! *****************************************************************************
    SUBROUTINE se_ga_ks_accumulate_1 ( ga_env, ks_i, ioff, isize )
      TYPE ( ga_environment_type ), POINTER :: ga_env
      REAL ( dp ), DIMENSION ( :, :, : ), POINTER :: ks_i
      INTEGER, INTENT ( IN ) :: ioff, isize
#ifdef _USE_GA
      INTEGER :: ilo ( 3 ), ihi ( 3 )
      REAL ( dp ), PARAMETER :: one=1.0_dp
      ilo(1) = 1
      ilo(2) = 1
      ilo(3) = ioff + 1
      ihi(1) = ga_env % se_dims(1)
      ihi(2) = ga_env % se_dims(2)
      ihi(3) = ioff + isize
      CALL nga_acc ( ga_env % g_ks, ilo, ihi, ks_i, ga_env % se_dims ( 1 ), one )
#endif
    END SUBROUTINE se_ga_ks_accumulate_1
! *****************************************************************************
    SUBROUTINE se_ga_deallocate_local_2 ( dens_i_info, dens_j_info, ks_i, ks_j, error )
     REAL ( dp ), POINTER :: dens_i_info ( :, : ), dens_j_info ( :, : )
     REAL ( dp ), POINTER :: ks_i ( :, :, : ), ks_j ( :, :, : )
     TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
     INTEGER :: stat
     CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_deallocate_local_2', &
       routineP = moduleN//':'//routineN
     LOGICAL :: failure

     IF ( ASSOCIATED ( dens_i_info ) )  DEALLOCATE ( dens_i_info, STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF ( ASSOCIATED ( dens_j_info ) )  DEALLOCATE ( dens_j_info, STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF ( ASSOCIATED ( ks_i ) )  DEALLOCATE ( ks_i, STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF ( ASSOCIATED ( ks_j ) )  DEALLOCATE ( ks_j, STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     NULLIFY ( dens_i_info )
     NULLIFY ( dens_j_info )
     NULLIFY ( ks_i )
     NULLIFY ( ks_j )
#endif
    END SUBROUTINE se_ga_deallocate_local_2
! *****************************************************************************
    SUBROUTINE se_ga_deallocate_local_1 ( dens_i_info, ks_i, error )
     REAL ( dp ), POINTER, OPTIONAL :: dens_i_info ( :, : )
     REAL ( dp ), POINTER, OPTIONAL:: ks_i ( :, :, : )
     TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
     INTEGER :: stat
     CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_deallocate_local_1', &
       routineP = moduleN//':'//routineN
     LOGICAL :: failure

     IF ( PRESENT ( dens_i_info ) ) THEN
       IF ( ASSOCIATED ( dens_i_info ) )  DEALLOCATE ( dens_i_info, STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY ( dens_i_info )
     END IF
     IF ( PRESENT ( ks_i ) ) THEN
       IF ( ASSOCIATED ( ks_i ) )  DEALLOCATE ( ks_i, STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY ( ks_i )
     END IF
#endif
    END SUBROUTINE se_ga_deallocate_local_1
! *****************************************************************************
    SUBROUTINE se_ga_allocate_local_2 ( ga_env, dens_i_info, dens_j_info, ks_i, ks_j, &
                                          ibin, jbin, isize, jsize, ioff, joff, error )
     TYPE ( ga_environment_type ), POINTER :: ga_env
     INTEGER, INTENT ( IN ) :: ibin, jbin
     INTEGER, INTENT ( OUT ) :: isize, jsize, ioff, joff
        REAL ( dp ), POINTER :: dens_i_info ( :, : ), dens_j_info ( :, : )
REAL ( dp ), POINTER :: ks_i ( :, :, : ), ks_j ( :, :, : )
        TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
        CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_allocate_local_2', &
        routineP = moduleN//':'//routineN
        LOGICAL ::  failure
     INTEGER, PARAMETER  :: one=1
     INTEGER :: ilo ( 2 ), jlo ( 2 ), ihi ( 2 ), jhi ( 2 ), stat
     CALL nga_get(ga_env%g_offset,ibin,ibin,ioff,one)
     CALL nga_get(ga_env%g_cellnum,ibin,ibin,isize,one)
     CALL nga_get(ga_env%g_offset,jbin,jbin,joff,one)
     CALL nga_get(ga_env%g_cellnum,jbin,jbin,jsize,one)
     IF ( ( isize == 0 ) .OR. ( jsize == 0 ) ) RETURN
     ilo(1) = 1
     ilo(2) = ioff + 1
     jlo(1) = 1
     jlo(2) = joff + 1
     ihi(1) = ga_env % se_dims(1)*ga_env % se_dims(1) + 5
     ihi(2) = ioff + isize
     jhi(1) = ga_env % se_dims(1)*ga_env % se_dims(1) + 5
     jhi(2) = joff + jsize
     ALLOCATE ( dens_i_info ( ihi ( 1 ), isize ), STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE ( dens_j_info ( jhi ( 1 ), jsize ), STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE ( ks_i ( ga_env % se_dims ( 1 ), ga_env % se_dims ( 1 ), isize ), STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE ( ks_j ( ga_env % se_dims ( 1 ), ga_env % se_dims ( 1 ), jsize ), STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ks_i = 0.0_dp
     ks_j = 0.0_dp
     CALL nga_get(ga_env % g_d_info,ilo,ihi,dens_i_info,ihi(1))
     CALL nga_get(ga_env % g_d_info,jlo,jhi,dens_j_info,jhi(1))
#else
     isize = 1
     jsize = 1
     ioff = 0
     joff = 0
#endif
    END SUBROUTINE se_ga_allocate_local_2
! *****************************************************************************
    SUBROUTINE se_ga_allocate_local_1 ( ga_env, dens_i_info, ks_i,&
                                       ibin, isize, ioff,  error )
     TYPE ( ga_environment_type ), POINTER :: ga_env
     INTEGER, INTENT ( IN ) :: ibin
     INTEGER, INTENT ( OUT ):: isize
     INTEGER, INTENT ( OUT ), OPTIONAL :: ioff
     REAL ( dp ), POINTER :: dens_i_info ( :, : )
     REAL ( dp ), POINTER, OPTIONAL :: ks_i ( :, :, : )
     TYPE(cp_error_type), INTENT(inout) :: error
#ifdef _USE_GA
     CHARACTER(len=*), PARAMETER :: routineN = 'se_ga_allocate_local_1', &
       routineP = moduleN//':'//routineN
     LOGICAL ::  failure
     INTEGER :: myioff
     INTEGER, PARAMETER  :: one=1
     INTEGER :: ilo ( 2 ), ihi ( 2 ), stat
     CALL nga_get(ga_env%g_offset,ibin,ibin,myioff,one)
     CALL nga_get(ga_env%g_cellnum,ibin,ibin,isize,one)
     IF ( isize == 0 ) RETURN
     IF ( PRESENT ( ioff ) )   ioff =  myioff
     ilo(1) = 1
     ilo(2) = myioff + 1
     ihi(1) = ga_env % se_dims(1) * ga_env % se_dims(1) + 5
     ihi(2) = myioff + isize
     ALLOCATE ( dens_i_info ( ihi ( 1 ), isize ), STAT=stat )
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     CALL nga_get(ga_env % g_d_info,ilo,ihi,dens_i_info,ihi(1))
     IF ( PRESENT ( ks_i ) ) THEN
      ALLOCATE ( ks_i ( ga_env % se_dims ( 1 ), ga_env % se_dims ( 1 ), isize ), STAT=stat )
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ks_i = 0.0_dp
     END IF
#else
     isize = 1
     IF ( PRESENT ( ioff ) )   ioff = 0
#endif
    END SUBROUTINE se_ga_allocate_local_1

! *****************************************************************************
    SUBROUTINE se_ga_get_nbin ( ga_handle, nbin, lzero )
     INTEGER, INTENT ( IN ) :: ga_handle
     INTEGER, INTENT ( OUT ) :: nbin
     LOGICAL, INTENT ( IN ), OPTIONAL :: lzero
#ifdef _USE_GA
     LOGICAL :: my_lzero
     INTEGER, PARAMETER  :: one=1
     my_lzero = .FALSE.
     IF ( PRESENT ( lzero ) ) my_lzero = lzero
     IF ( my_lzero ) CALL ga_zero ( ga_handle )
     nbin = nga_read_inc ( ga_handle, one, one )
#else
     nbin = 1
#endif
    END SUBROUTINE se_ga_get_nbin

! *****************************************************************************
    SUBROUTINE grid_factor(p,xdim,ydim,zdim,idx,idy,idz)
    INTEGER, INTENT(IN)                      :: p
    REAL(dp), INTENT(IN)                     :: xdim, ydim, zdim
    INTEGER, INTENT(OUT)                     :: idx, idy, idz

    INTEGER, PARAMETER                       :: MAX_FACTOR = 10000

    INTEGER                                  :: fac(MAX_FACTOR), i, ifac, ip, &
                                                ix, iy, iz, j, pmax, &
                                                prime(MAX_FACTOR)

!

       i = 1
!
!    factor p completely
!    first, find all prime numbers, besides 1, less than or equal to
!    the square root of p
!
      ip = INT(SQRT(DBLE(p)))+1
      pmax = 0
      DO i = 2, ip
        DO j = 1, pmax
          IF (MOD(i,prime(j))==0) go to 100
        END DO
        pmax = pmax + 1
        IF (pmax>MAX_FACTOR) WRITE(6,*) 'Overflow in grid_factor'
        prime(pmax) = i
        100   CONTINUE
      END DO
!
!    find all prime factors of p
!
      ip = p
      ifac = 0
      DO i = 1, pmax
        200   IF (MOD(ip,prime(i))==0) THEN
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        ENDIF
      END DO
!
!    determine three factors of p of approximately the
!    same size
!
      idx = 1
      idy = 1
      idz = 1
      DO i = ifac, 1, -1
        ix = xdim / idx
        iy = ydim / idy
        iz = zdim / idz
        IF (ix>=iy.and.ix>=iz.and.ix>1) THEN
          idx = fac(i)*idx
        ELSEIF (iy>=ix.and.iy>=iz.and.iy>1) THEN
          idy = fac(i)*idy
        ELSEIF (iz>=ix.and.iz>=iy.and.iz>1) THEN
          idz = fac(i)*idz
        ELSE
          WRITE(6,*) 'Too many processors in grid factoring routine'
        ENDIF
      END DO
    END SUBROUTINE grid_factor
! *****************************************************************************
    SUBROUTINE pbc_wrap ( r, x, y, z )
    REAL(kind=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r
    REAL(kind=dp), INTENT(IN)                :: x, y, z

      DO WHILE ( r ( 1 ) > x )
        r ( 1 ) = r ( 1 ) - x
      END DO
      DO WHILE ( r ( 2 ) > y )
        r ( 2 ) = r ( 2 ) - y
      END DO
      DO WHILE ( r ( 3 ) > z )
        r ( 3 ) = r ( 3 ) - z
      END DO
      DO WHILE ( r ( 1 ) < 0.0_dp )
        r ( 1 ) = r ( 1 ) + x
      END DO
      DO WHILE ( r ( 2 ) < 0.0_dp )
        r ( 2 ) = r ( 2 ) + y
      END DO
      DO WHILE ( r ( 3 ) < 0.0_dp )
        r ( 3 ) = r ( 3 ) + z
      END DO
    END SUBROUTINE pbc_wrap
! *****************************************************************************
END MODULE se_ga_tools
