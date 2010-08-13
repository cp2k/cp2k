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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cluster_types,                   ONLY: box_type,&
                                             cluster_type
  USE input_constants,                 ONLY: clustering_box,&
                                             clustering_h_on_heavier,&
                                             clustering_none
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: radians
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cluster_methods'

  LOGICAL, PARAMETER, PRIVATE          :: dbg = .FALSE., careful_mod = .TRUE.

! *** Public subroutines ***

  PUBLIC :: cluster_allocate, cluster_deallocate
  PUBLIC :: clustering_particles
  PUBLIC :: cluster_part2cluster, cluster_member_particles

CONTAINS

! *****************************************************************************
!> \brief Clustering particles
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE clustering_particles(clusters, clustering_section, cell, &
       particle_set, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(section_vals_type), POINTER         :: clustering_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'clustering_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, method
    LOGICAL                                  :: failure
    REAL(dp)                                 :: box_size
    TYPE(box_type)                           :: boxes

    failure = .FALSE.

    CALL timeset(routineN,handle)

    CALL section_vals_val_get(clustering_section,"METHOD",i_val=method,error=error)
    CALL section_vals_val_get(clustering_section,"BOX_SIZE",r_val=box_size,error=error)
    !
    ! box my particles
    CALL boxing_particles(cell, particle_set, boxes, box_size, error=error)
    !
    ! cluster my particles
    SELECT CASE(method)
    CASE(clustering_none) ! ATOM
       CALL cluster_none(clusters, boxes, error)
    CASE(clustering_box) ! BOX
       CALL cluster_box(clusters, boxes, error)
    CASE(clustering_h_on_heavier)
       CALL cluster_h_on_heavier(clusters, boxes, cell, particle_set, error)
    !CASE(clustering_smallest_on_larger)
    !call cluster_smallest_on_larger(clusters, boxes, particle_set, error)
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    IF(dbg) CALL cluster_print(clusters, cell, particle_set, error)

    !
    ! check if we didnt loose any particles...
    CALL cluster_check(clusters, error)

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
!> \brief Clustering the hydrogen particle block on larger one based on boxes
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE cluster_h_on_heavier(clusters, boxes, cell, particle_set, error)

    TYPE(cluster_type)                       :: clusters
    TYPE(box_type)                           :: boxes
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_h_on_heavier', &
      routineP = moduleN//':'//routineN
    REAL(dp), PARAMETER                      :: radius = 2.0_dp

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER :: i, ibox, ipart, istat, j, jbox, jpart, kbox, lbox, mbox, &
      mpart, n_part, nbox, nclusters, nnbox(3), npart, nparticles
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: clustered
    INTEGER, DIMENSION(:), POINTER           :: cluster_ptr, part_i, &
                                                particle_idx
    INTEGER, DIMENSION(:, :, :), POINTER     :: box_p
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dist, r(3)
    TYPE(atomic_kind_type), POINTER          :: atom_kind

    failure = .FALSE.

    nparticles = boxes%npart
    nclusters = boxes%npart ! worst case

    CALL cluster_allocate(clusters, nclusters, nparticles, error)

    ALLOCATE(clustered(nparticles),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    clustered(:) = 0

    DO ipart = 1,nparticles
       atom_kind => particle_set(ipart)%atomic_kind
       CALL get_atomic_kind(atom_kind,element_symbol=element_symbol)
       IF(TRIM(element_symbol).EQ."H") clustered(ipart) = 1
    ENDDO

    nclusters = 1
    nparticles = 0
    particle_idx => clusters%particle_idx
    cluster_ptr => clusters%cluster_ptr
    cluster_ptr(nclusters) = 1
    box_p => boxes%box_p
    part_i => boxes%part_i
    nnbox(:) = boxes%nbox(:)
    DO kbox = 0,nnbox(3)-1
    DO jbox = 0,nnbox(2)-1
    DO ibox = 0,nnbox(1)-1
       mpart = box_p(ibox+1,jbox,kbox)-box_p(ibox,jbox,kbox)
       IF(mpart.EQ.0) CYCLE
       !
       ! cluster h on heavier
       DO i = box_p(ibox,jbox,kbox),box_p(ibox+1,jbox,kbox)-1
          ipart = part_i(i)
          !
          ! find a none hydrogen atom
          IF(clustered(ipart).NE.0) CYCLE
          !
          ! add it to the cluster
          n_part = 1
          nparticles = nparticles + 1
          particle_idx(nparticles) = ipart
          !
          ! find the hydrogen atoms that are not clustered yet
          DO nbox = MAX(kbox-1,0),MIN(kbox+1,nnbox(3)-1)
          DO mbox = MAX(jbox-1,0),MIN(jbox+1,nnbox(2)-1)
          DO lbox = MAX(ibox-1,0),MIN(ibox+1,nnbox(1)-1)
             npart = box_p(lbox+1,mbox,nbox)-box_p(lbox,mbox,nbox)
             IF(npart.EQ.0) CYCLE
             DO j = box_p(lbox,mbox,nbox),box_p(lbox+1,mbox,nbox)-1
                jpart = part_i(j)
                IF(clustered(jpart).NE.1) CYCLE
                !
                ! if possible add it to the cluster
                r = pbc(particle_set(ipart)%r(:), particle_set(jpart)%r(:), cell)
                dist = SQRT(DOT_PRODUCT(r,r))
                IF(dist.LE.radius) THEN
                   clustered(jpart) = -1
                   n_part = n_part + 1
                   nparticles = nparticles + 1
                   particle_idx(nparticles) = jpart
                ENDIF
             ENDDO
          ENDDO !lbox
          ENDDO !mbox
          ENDDO !nbox
          !
          ! increament the cluster
          nclusters = nclusters + 1
          cluster_ptr(nclusters) = cluster_ptr(nclusters-1) + n_part
       ENDDO !i
       !
       !
    ENDDO !ibox
    ENDDO !jbox
    ENDDO !kbox
    !
    ! add the remaining hydrogens to the cluster 
    IF(nparticles.NE.boxes%npart) THEN
       DO kbox = 0,nnbox(3)-1
       DO jbox = 0,nnbox(2)-1
       DO ibox = 0,nnbox(1)-1
          mpart = box_p(ibox+1,jbox,kbox)-box_p(ibox,jbox,kbox)
          IF(mpart.EQ.0) CYCLE
          n_part = 0
          !
          ! box remaining h
          DO i = box_p(ibox,jbox,kbox),box_p(ibox+1,jbox,kbox)-1
             ipart = part_i(i)
             !
             ! find a none clustered hydrogen atom
             IF(clustered(ipart).NE.1) CYCLE
             n_part = n_part + 1
             clustered(ipart) = -1
             nparticles = nparticles + 1
             particle_idx(nparticles) = ipart
          ENDDO
          !
          ! increment the cluster if found any particles
          IF(n_part.GT.0) THEN
             nclusters = nclusters + 1
             cluster_ptr(nclusters) = cluster_ptr(nclusters-1) + n_part
          ENDIF
       ENDDO !ibox
       ENDDO !jbox
       ENDDO !kbox
    ENDIF
    !
    ! reset to the correct value
    clusters%nclusters = nclusters - 1

    DEALLOCATE(clustered,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE cluster_h_on_heavier

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
  SUBROUTINE boxing_particles(cell, particle_set, boxes, box_size_guess, error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(box_type)                           :: boxes
    REAL(kind=dp), INTENT(in)                :: box_size_guess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'boxing_particles', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: box_secure = 1.0e-6_dp

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
    ! wrap the atoms and center them
    cell_center(:) = 0.5_dp*( h(:,1) + h(:,2) + h(:,3) )
    DO iatom = 1,natms
       xyz(:,iatom) = pbc(particle_set(iatom)%r(:), cell) + cell_center(:)
    ENDDO

    !
    ! g = canonical h
    g(:,:) = 0.0_dp
    g(1,1) = abc(1)
    g(1,2) = abc(2) * COS(radians * gamma)
    g(2,2) = abc(3) * SIN(radians * gamma)
    g(1,3) = abc(3) * COS(radians * beta)
    g(2,3) = ( abc(2) * abc(3) * COS(radians * alpha) - g(1,2) * g(1,3) ) / g(2,2)
    g(3,3) = SQRT( abc(3)**2 - g(1,3)**2 - g(2,3)**2 )

    !
    ! compute u matrix
    u = MATMUL(g,h_inv)

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
    ! shift the xyz and set the size of the box and secure it
    xyz(:,:) = xyz(:,:) + box_secure / 2.0_dp
    DO ixyz = 1,3
       lbox(ixyz) = g(ixyz,ixyz) + box_secure
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
       WRITE(*,'(T2,A,I6)') routineN//' natms=',natms
       WRITE(*,'(T2,A,3I2)') routineN//' periodic=',periodic
       WRITE(*,'(T2,A,L1)') routineN//' orthorhombic=',orthorhombic
       WRITE(*,'(T2,A,3F12.4)') routineN//' h(:,1)=',h(:,1)
       WRITE(*,'(T2,A,3F12.4)') routineN//' h(:,2)=',h(:,2)
       WRITE(*,'(T2,A,3F12.4)') routineN//' h(:,3)=',h(:,3)
       WRITE(*,'(T2,A,3F12.4)') routineN//' g(:,1)=',g(:,1)
       WRITE(*,'(T2,A,3F12.4)') routineN//' g(:,2)=',g(:,2)
       WRITE(*,'(T2,A,3F12.4)') routineN//' g(:,3)=',g(:,3)
       WRITE(*,'(T2,A,3F12.4)') routineN//' u(:,1)=',u(:,1)
       WRITE(*,'(T2,A,3F12.4)') routineN//' u(:,2)=',u(:,2)
       WRITE(*,'(T2,A,3F12.4)') routineN//' u(:,3)=',u(:,3)
       WRITE(*,'(T2,A,F12.4)') routineN//'  alpha=',alpha
       WRITE(*,'(T2,A,F12.4)') routineN//'   beta=',beta
       WRITE(*,'(T2,A,F12.4)') routineN//'  gamma=',gamma
       WRITE(*,'(T2,A,3F12.4)') routineN//'    abc=',abc
       WRITE(*,'(T2,A,3F12.4)') routineN//'lbox(:)=',lbox(:)
       WRITE(*,'(T2,A,3I4)') routineN//'   nbox=',nbox
       WRITE(*,'(T2,A,3F12.4)') routineN//'   dbox=',dbox
       WRITE(*,'(T2,A,F12.4)') routineN//' box_size_guess=',box_size_guess
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

  SUBROUTINE cluster_check(clusters, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icluster, ip, iparticle, istat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: icheck
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE(icheck(clusters%nparticles),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    icheck(:) = 0

    DO icluster = 1,clusters%nclusters
       DO ip = clusters%cluster_ptr(icluster),clusters%cluster_ptr(icluster+1)-1
          iparticle = clusters%particle_idx(ip)
          IF(icheck(iparticle).NE.0) THEN
             WRITE(*,*) ' ERROR: The particle ',iparticle,' is present more than once...'
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ELSE
             icheck(iparticle) =  1
          ENDIF
       ENDDO
    ENDDO

    IF(SUM(icheck).NE.clusters%nparticles) THEN
       WRITE(*,*) ' ERROR: The number of particles is not correct...'
       WRITE(*,*) ' clusters%nparticles', clusters%nparticles,' SUM(icheck)',SUM(icheck)
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    DEALLOCATE(icheck,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE cluster_check

  SUBROUTINE cluster_print(clusters, cell, particle_set, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: icluster, id, ip, iparticle
    TYPE(atomic_kind_type), POINTER          :: atom_kind

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
 
    IF(.TRUE.) THEN
       id = 126
       OPEN(id)
       DO icluster = 1,clusters%nclusters
          WRITE(id,'(I6)') clusters%cluster_ptr(icluster+1)-clusters%cluster_ptr(icluster)
          WRITE(id,'(A,I6)') 'icluster',icluster
          DO ip = clusters%cluster_ptr(icluster),clusters%cluster_ptr(icluster+1)-1
             iparticle = clusters%particle_idx(ip)
             atom_kind => particle_set(iparticle)%atomic_kind
             CALL get_atomic_kind(atom_kind,element_symbol=element_symbol)
             WRITE(id,'(A,3F15.4)') element_symbol,pbc(particle_set(iparticle)%r(:),cell)*angstrom
          ENDDO
       ENDDO
       CLOSE(id)
    ENDIF

  END SUBROUTINE cluster_print


  FUNCTION cluster_part2cluster (clusters, particle, error) RESULT (cluster)
    TYPE(cluster_type), INTENT(IN)           :: clusters
    INTEGER, INTENT(IN)                      :: particle
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: cluster

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_part2cluster', &
      routineP = moduleN//':'//routineN

    IF (careful_mod) THEN
       CALL cp_assert (ASSOCIATED (clusters%particle_to_cluster),&
            cp_fatal_level, cp_caller_error, routineN,&
            "Clusters not setup.", error=error)
       CALL cp_assert (SIZE (clusters%particle_to_cluster) .EQ. clusters%nparticles,&
            cp_warning_level, cp_internal_error, routineN,&
            "Mismatch between declared and actual number of particles.",&
            error=error)
       CALL cp_assert (particle .LE. SIZE (clusters%particle_to_cluster),&
            cp_fatal_level, cp_caller_error, routineN,&
            "Invalid particle number.", error=error)
    ENDIF
    cluster = clusters%particle_to_cluster(particle)
    IF (careful_mod) THEN
       CALL cp_assert (cluster .LE. clusters%nclusters, cp_warning_level,&
            cp_internal_error, routineN, "Obtained cluster nonexistent.",&
            error=error)
    ENDIF
  END FUNCTION cluster_part2cluster

  FUNCTION cluster_member_particles (clusters, cluster, error) RESULT (members)
    TYPE(cluster_type), INTENT(IN)           :: clusters
    INTEGER, INTENT(IN)                      :: cluster
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, DIMENSION(:), POINTER           :: members

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_member_particles', &
      routineP = moduleN//':'//routineN

    IF (careful_mod) THEN
       CALL cp_assert (ASSOCIATED (clusters%cluster_ptr), cp_fatal_level,&
            cp_caller_error, routineN, "Clusters not setup.", error=error)
       CALL cp_assert (ASSOCIATED (clusters%particle_idx), cp_fatal_level,&
            cp_caller_error, routineN, "Clusters not setup.", error=error)
       CALL cp_assert (SIZE (clusters%cluster_ptr) .EQ. clusters%nclusters+1,&
            cp_warning_level, cp_internal_error, routineN,&
            "Mismatch between declared and actual number of particles.",&
            error=error)
       CALL cp_assert (cluster .LE. SIZE (clusters%cluster_ptr),&
            cp_fatal_level, cp_caller_error, routineN,&
            "Invalid cluster number.", error=error)
    ENDIF
    members => clusters%particle_idx(clusters%cluster_ptr(cluster)&
                                    :clusters%cluster_ptr(cluster+1)-1)
    IF (careful_mod) THEN
       CALL cp_assert (cluster .LE. clusters%nclusters, cp_warning_level,&
            cp_internal_error, routineN, "Obtained cluster nonexistent.",&
            error=error)
    ENDIF
  END FUNCTION cluster_member_particles


END MODULE cluster_methods
