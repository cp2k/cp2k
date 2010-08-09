!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Cluster methods for atoms, particles, or molecules
!> \par History
!> \author VW (2010)
! *****************************************************************************
MODULE cluster_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cluster_types,                   ONLY: box_type,&
                                             cluster_type
  USE input_constants,                 ONLY: clustering_box,&
                                             clustering_none
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: radians
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cluster_methods'

  LOGICAL, PARAMETER, PRIVATE          :: dbg = .FALSE.

! *** Public subroutines ***

  PUBLIC :: clustering_particles, cluster_deallocate

CONTAINS

! *****************************************************************************
!> \brief Clustering particles
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE clustering_particles(clusters, clustering_section, cell, &
       particle_kind_set, particle_set, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(section_vals_type), POINTER         :: clustering_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: particle_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'clustering_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, method
    LOGICAL                                  :: failure
    TYPE(box_type)                           :: boxes

    failure = .FALSE.

    CALL timeset(routineN,handle)

    CALL section_vals_val_get(clustering_section,"METHOD",i_val=method,error=error)
    !
    ! box my particles
    CALL boxing_particles(cell, particle_kind_set, particle_set, boxes, error=error)
    !
    ! cluster my particles
    SELECT CASE(method)
    CASE(clustering_none) ! ATOM
       CALL cluster_none(clusters, boxes, error)
    CASE(clustering_box) ! BOX
       CALL cluster_box(clusters, boxes, error)
    !CASE(clustering_h_on_heavier)
    ! call cluster_h_on_heavier(clusters, boxes, error)
    !CASE(clustering_smallest_on_larger)
    !call cluster_smallest_on_larger(clusters, boxes, particle_kind_set, particle_set, error)
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    IF(dbg) CALL cluster_print(clusters, error)

    !
    ! cleanup
    CALL box_deallocate(boxes, error=error)

    CALL timestop(handle)

  END SUBROUTINE clustering_particles

! *****************************************************************************
!> \brief Clustering the particles based on boxes
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE cluster_box(clusters, boxes, error)

    TYPE(cluster_type)                       :: clusters
    TYPE(box_type)                           :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_box', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibox, jbox, kbox, &
                                                nclusters, npart, nparticles
    INTEGER, DIMENSION(:), POINTER           :: part_i
    INTEGER, DIMENSION(:, :, :), POINTER     :: box_p

!

    nparticles = boxes%npart
    nclusters = boxes%nbox_filled

    CALL cluster_allocate(clusters, nclusters, nparticles, error)

    nclusters = 1
    nparticles = 0
    clusters%cluster_ptr(nclusters) = 1
    box_p => boxes%box_p
    part_i => boxes%part_i
    DO kbox = 0,boxes%nbox(3)-1
    DO jbox = 0,boxes%nbox(2)-1
    DO ibox = 0,boxes%nbox(1)-1
       npart = box_p(ibox+1,jbox,kbox)-box_p(ibox,jbox,kbox)
       IF(npart.GT.0) THEN
          DO i = box_p(ibox,jbox,kbox),box_p(ibox+1,jbox,kbox)-1
             nparticles = nparticles + 1
             clusters%particle_idx(nparticles) = part_i(i)
          ENDDO
          nclusters = nclusters + 1
          clusters%cluster_ptr(nclusters) = clusters%cluster_ptr(nclusters-1) + npart
       ENDIF
    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE cluster_box

! *****************************************************************************
!> \brief Clustering the smallest particle block on larger one based on boxes
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE cluster_smallest_on_larger(clusters, boxes, error)

    TYPE(cluster_type)                       :: clusters
    TYPE(box_type)                           :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_smallest_on_larger', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibox, jbox, kbox, &
                                                nclusters, npart, nparticles
    INTEGER, DIMENSION(:), POINTER           :: part_i
    INTEGER, DIMENSION(:, :, :), POINTER     :: box_p
    LOGICAL                                  :: failure

    failure = .FALSE.

!
! not yet tere...

    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)

    nparticles = boxes%npart
    nclusters = boxes%npart ! worst case

    CALL cluster_allocate(clusters, nclusters, nparticles, error)

    nclusters = 1
    nparticles = 0
    clusters%cluster_ptr(nclusters) = 1
    box_p => boxes%box_p
    part_i => boxes%part_i
    DO kbox = 0,boxes%nbox(3)-1
    DO jbox = 0,boxes%nbox(2)-1
    DO ibox = 0,boxes%nbox(1)-1
       npart = box_p(ibox+1,jbox,kbox)-box_p(ibox,jbox,kbox)
       IF(npart.GT.0) THEN
          DO i = box_p(ibox,jbox,kbox),box_p(ibox+1,jbox,kbox)-1
             nparticles = nparticles + 1
             clusters%particle_idx(nparticles) = part_i(i)
          ENDDO
          nclusters = nclusters + 1
          clusters%cluster_ptr(nclusters) = clusters%cluster_ptr(nclusters-1) + npart
       ENDIF
    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE cluster_smallest_on_larger

! *****************************************************************************
!> \brief Clustering none 
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE cluster_none(clusters, boxes, error)

    TYPE(cluster_type)                       :: clusters
    TYPE(box_type)                           :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibox, jbox, kbox, &
                                                nclusters, npart, nparticles
    INTEGER, DIMENSION(:), POINTER           :: part_i
    INTEGER, DIMENSION(:, :, :), POINTER     :: box_p

!

    nparticles = boxes%npart
    nclusters = boxes%npart

    CALL cluster_allocate(clusters, nclusters, nparticles, error)

    nclusters = 1
    nparticles = 0
    clusters%cluster_ptr(nclusters) = 1
    box_p => boxes%box_p
    part_i => boxes%part_i
    DO kbox = 0,boxes%nbox(3)-1
    DO jbox = 0,boxes%nbox(2)-1
    DO ibox = 0,boxes%nbox(1)-1
       npart = box_p(ibox+1,jbox,kbox)-box_p(ibox,jbox,kbox)
       IF(npart.GT.0) THEN
          DO i = box_p(ibox,jbox,kbox),box_p(ibox+1,jbox,kbox)-1
             nparticles = nparticles + 1
             clusters%particle_idx(nparticles) = part_i(i)
             nclusters = nclusters + 1
             clusters%cluster_ptr(nclusters) = clusters%cluster_ptr(nclusters-1) + 1
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE cluster_none

! *****************************************************************************
!> \brief Box my particles
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE boxing_particles(cell, particle_kind_set,particle_set,&
       boxes, error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: particle_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(box_type)                           :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'boxing_particles', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: box_secure = 1.0_dp, &
                                                box_size_guess = 8.0_dp

    INTEGER                                  :: handle, i, iatom, ibox, &
                                                istat, ixyz, jbox, kbox, &
                                                natms, nparticles
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cnt
    INTEGER, DIMENSION(3)                    :: nbox, periodic
    INTEGER, DIMENSION(:), POINTER           :: part_i
    INTEGER, DIMENSION(:, :, :), POINTER     :: box_p
    LOGICAL                                  :: failure, orthorhombic
    REAL(dp)                                 :: alpha, beta, gamma
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: xyz, xyz_tmp
    REAL(dp), DIMENSION(3)                   :: abc, cell_center, dbox, lbox
    REAL(dp), DIMENSION(3, 3)                :: g, h, h_inv, u
    REAL(dp), DIMENSION(:, :), POINTER       :: xyz_p
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    natms = SIZE(particle_set,1)
    !
    !
    CALL get_cell(cell, &
                  orthorhombic=orthorhombic, &
                  periodic=periodic, &
                  h=h, h_inv=h_inv, &
                  alpha=alpha, &
                  beta=beta, &
                  gamma=gamma, &
                  abc=abc)
    !
    !
    ALLOCATE(xyz(3,natms),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    ! g = canonical h
    g(:,:) = 0.0_dp
    g(1,1) = 1.0_dp; g(2,2) = 1.0_dp; g(3,3) = 1.0_dp;
    SELECT CASE(SUM(periodic))
    CASE(0)
    CASE(1)
       DO ixyz = 1,3
          IF(periodic(ixyz).EQ.1) g(ixyz,ixyz) = abc(ixyz)
       ENDDO
    CASE(2)
       ! some more work here
       IF(periodic(3).EQ.0) THEN
          g(1,1) = abc(1)
          g(1,2) = abc(2) * COS(radians * abc(3))
          g(2,2) = abc(2) * SIN(radians * abc(3))
       ENDIF
       IF(periodic(2).EQ.0) THEN
          g(1,1) = abc(1)
          g(1,3) = abc(2) * COS(radians * abc(3))
          g(3,3) = abc(2) * SIN(radians * abc(3))
       ENDIF
       IF(periodic(1).EQ.0) THEN
          g(2,2) = abc(1)
          g(2,3) = abc(2) * COS(radians * abc(3))
          g(3,3) = abc(2) * SIN(radians * abc(3))
       ENDIF
    CASE(3)
       g(1,1) = abc(1)
       g(1,2) = abc(2) * COS(radians * gamma)
       g(2,2) = abc(3) * SIN(radians * gamma)
       g(1,3) = abc(3) * COS(radians * beta)
       g(2,3) = (abc(2) * abc(3) * COS(radians * alpha) - g(1,2) * g(1,3)) / g(2,2)
       g(3,3) = SQRT(abc(3)**2 - g(1,3)**2 - g(2,3)**2)
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    !
    ! compute u matrix
    u = MATMUL(g,h_inv)
    !
    ! wrap the atoms and center them
    cell_center(:) = 0.5_dp*( h(:,1) + h(:,2) + h(:,3) )
    DO iatom = 1,natms
       xyz(:,iatom) = pbc(particle_set(iatom)%r(:), cell) + cell_center(:)
    ENDDO
    !
    ! center for non periodic
    DO ixyz = 1,3
       IF(periodic(ixyz).NE.1) xyz(ixyz,:) = xyz(ixyz,:) - MINVAL(xyz(ixyz,:)) + box_secure / 2.0_dp
    ENDDO
    !
    ! rotate coordinate if needed
    IF(.NOT.orthorhombic) THEN
       ALLOCATE(xyz_tmp(3,natms),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL dcopy(3*natms, xyz, 1, xyz_tmp, 1)
       CALL dgemm('n', 'n', 3, natms, 3, 1.0_dp, u, 3, xyz_tmp, 3 , 0.0_dp, xyz, 3)
       DEALLOCATE(xyz_tmp,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    !
    ! compute the size of the box
    DO ixyz = 1,3
       IF(periodic(ixyz).EQ.1) THEN
          lbox(ixyz) = g(ixyz,ixyz)
       ELSE
          lbox(ixyz) = MAXVAL(xyz(ixyz,:)) + box_secure
       ENDIF
    ENDDO
    !
    ! box the atoms
    nbox(1) = CEILING(lbox(1) / box_size_guess)
    nbox(2) = CEILING(lbox(2) / box_size_guess)
    nbox(3) = CEILING(lbox(3) / box_size_guess)
    dbox(1) = lbox(1) / REAL(nbox(1),dp)
    dbox(2) = lbox(2) / REAL(nbox(2),dp)
    dbox(3) = lbox(3) / REAL(nbox(3),dp)

    !
    !
    IF(dbg) THEN
       WRITE(*,*) routineN//' natms=',natms
       WRITE(*,*) routineN//' periodic=',periodic
       WRITE(*,*) routineN//' h(:,1)=',h(:,1)
       WRITE(*,*) routineN//' h(:,2)=',h(:,2)
       WRITE(*,*) routineN//' h(:,3)=',h(:,3)
       WRITE(*,*) routineN//' g(:,1)=',g(:,1)
       WRITE(*,*) routineN//' g(:,2)=',g(:,2)
       WRITE(*,*) routineN//' g(:,3)=',g(:,3)
       WRITE(*,*) routineN//' u(:,1)=',u(:,1)
       WRITE(*,*) routineN//' u(:,2)=',u(:,2)
       WRITE(*,*) routineN//' u(:,3)=',u(:,3)
       WRITE(*,*) routineN//'  alpha=',alpha
       WRITE(*,*) routineN//'   beta=',beta
       WRITE(*,*) routineN//'  gamma=',gamma
       WRITE(*,*) routineN//'    abc=',abc
       WRITE(*,*) routineN//'lbox(:)=',lbox(:)
       WRITE(*,*) routineN//'   nbox=',nbox
       WRITE(*,*) routineN//'   dbox=',dbox
       WRITE(*,*) routineN//' box_size_guess=',box_size_guess
    ENDIF

    CALL box_allocate(boxes, nbox, natms, error=error)

    !
    ! count particles per box
    ALLOCATE(cnt(0:nbox(1)-1,0:nbox(2)-1,0:nbox(3)-1))
    cnt(:,:,:) = 0
    DO iatom = 1,natms
       ibox = FLOOR(xyz(1,iatom) / dbox(1))
       jbox = FLOOR(xyz(2,iatom) / dbox(2))
       kbox = FLOOR(xyz(3,iatom) / dbox(3))
       IF(ibox.EQ.-1.OR.ibox.EQ.nbox(1).OR.&
          jbox.EQ.-1.OR.jbox.EQ.nbox(2).OR.&
          kbox.EQ.-1.OR.kbox.EQ.nbox(3)) THEN
          WRITE(*,*) 'ibox,jbox,kbox',ibox,jbox,kbox,xyz(:,iatom)
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       cnt(ibox,jbox,kbox) = cnt(ibox,jbox,kbox) + 1
    ENDDO

    !
    ! set the table
    box_p => boxes%box_p
    part_i => boxes%part_i
    xyz_p => boxes%xyz
    i = 1
    boxes%nbox_filled = 0
    DO kbox = 0,boxes%nbox(3)-1
    DO jbox = 0,boxes%nbox(2)-1
       box_p(0,jbox,kbox) = i
    DO ibox = 0,boxes%nbox(1)-1
       nparticles = cnt(ibox,jbox,kbox)
       IF(nparticles.GT.0) boxes%nbox_filled = boxes%nbox_filled + 1
       box_p(ibox+1,jbox,kbox) = box_p(ibox,jbox,kbox) + nparticles
       i = i + nparticles
    ENDDO
    ENDDO
    ENDDO

    !
    ! fill the table
    DO iatom = 1,natms
       ibox = FLOOR(xyz(1,iatom) / dbox(1))
       jbox = FLOOR(xyz(2,iatom) / dbox(2))
       kbox = FLOOR(xyz(3,iatom) / dbox(3))
       IF(ibox.EQ.-1.OR.ibox.EQ.nbox(1).OR.&
          jbox.EQ.-1.OR.jbox.EQ.nbox(2).OR.&
          kbox.EQ.-1.OR.kbox.EQ.nbox(3)) THEN
          WRITE(*,*) 'ibox,jbox,kbox',ibox,jbox,kbox,xyz(:,iatom)
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       i = box_p(ibox,jbox,kbox) + cnt(ibox,jbox,kbox) - 1
       part_i(i) = iatom
       xyz_p(:,i) = xyz(:,iatom)
       cnt(ibox,jbox,kbox) = cnt(ibox,jbox,kbox) - 1
    ENDDO

    !
    ! cleanup
    DEALLOCATE(xyz,cnt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE boxing_particles

  SUBROUTINE box_allocate(boxes, nbox, npart, error)
    TYPE(box_type)                           :: boxes
    INTEGER, DIMENSION(3), INTENT(in)        :: nbox
    INTEGER, INTENT(in)                      :: npart
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'box_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    boxes%nbox(:) = nbox(:)
    boxes%npart = npart
    ALLOCATE(boxes%box_p(0:nbox(1),0:nbox(2)-1,0:nbox(3)-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes%part_i(npart),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes%xyz(3,npart),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE box_allocate

  SUBROUTINE box_deallocate(boxes, error)
    TYPE(box_type)                           :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'box_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    boxes%nbox(:) = 0
    boxes%npart = 0
    IF(ASSOCIATED(boxes%box_p)) THEN
       DEALLOCATE(boxes%box_p,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF(ASSOCIATED(boxes%part_i)) THEN
       DEALLOCATE(boxes%part_i,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF(ASSOCIATED(boxes%xyz)) THEN
       DEALLOCATE(boxes%xyz,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
  END SUBROUTINE box_deallocate

  SUBROUTINE cluster_allocate(clusters, nclusters, nparticles, error)
    TYPE(cluster_type)                       :: clusters
    INTEGER, INTENT(in)                      :: nclusters, nparticles
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    clusters%nclusters = nclusters
    clusters%nparticles = nparticles
    ALLOCATE(clusters%cluster_ptr(nclusters+1), clusters%particle_idx(nparticles),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE cluster_allocate

  SUBROUTINE cluster_deallocate(clusters, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    DEALLOCATE(clusters%cluster_ptr, clusters%particle_idx,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE cluster_deallocate

  SUBROUTINE cluster_print(clusters, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icluster, ip, iparticle

    WRITE(*,'(/,/,T2,A)') 'CLUSTERING THE PARTICLES'
    WRITE(*,'(T5,A,I6,A)') 'Found ',clusters%nclusters,' clusters'
    WRITE(*,'(T5,A)') ' cluster   number of particles   particles'

    DO icluster = 1,clusters%nclusters
       WRITE(*,'(T5,I6,T20,I6)') icluster,&
            clusters%cluster_ptr(icluster+1)-clusters%cluster_ptr(icluster)
       IF(.TRUE.) THEN
          DO ip = clusters%cluster_ptr(icluster),clusters%cluster_ptr(icluster+1)-1
             iparticle = clusters%particle_idx(ip)
             WRITE(*,'(T38,I6)') iparticle
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE cluster_print

END MODULE cluster_methods
