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
  USE f77_blas
  USE input_constants,                 ONLY: clustering_box,&
                                             clustering_none
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
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
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: boxes

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
    !CASE("H_ON_HEAVIER")
    ! call cluster_h_on_heavier(clusters, boxes, error)
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
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_box', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ibox, jbox, kbox, &
                                                nclusters, nparticles

!
! First round, how many boxes do we have
!

    nclusters = 0
    nparticles = 0
    !WRITE(*,*) routineN//' bound',LBOUND(boxes),UBOUND(boxes)
    DO kbox = 0,UBOUND(boxes,3)-1
    DO jbox = 0,UBOUND(boxes,2)-1
    DO ibox = 0,UBOUND(boxes,1)-1
       IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
          nclusters = nclusters + 1
          nparticles = nparticles + boxes(ibox,jbox,kbox)%nparticles
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !
    !WRITE(*,*) 'nclusters',nclusters,' nparticles',nparticles
    CALL cluster_allocate(clusters, nclusters, nparticles, error)
    !
    ! Secound round, 
    !
    nclusters = 1
    nparticles = 0
    clusters%cluster_ptr(nclusters) = 1
    DO kbox = 0,UBOUND(boxes,3)-1
    DO jbox = 0,UBOUND(boxes,2)-1
    DO ibox = 0,UBOUND(boxes,1)-1
       IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
          DO iatom = 1,boxes(ibox,jbox,kbox)%nparticles
             nparticles = nparticles + 1
             clusters%particle_idx(nparticles) = boxes(ibox,jbox,kbox)%atom_list(iatom)
          ENDDO
          nclusters = nclusters + 1
          clusters%cluster_ptr(nclusters) = clusters%cluster_ptr(nclusters-1) + boxes(ibox,jbox,kbox)%nparticles
       ENDIF
    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE cluster_box

! *****************************************************************************
!> \brief Clustering none 
!> \par History
!>      none
!> \author VW (2010)
! *****************************************************************************
  SUBROUTINE cluster_none(clusters, boxes, error)

    TYPE(cluster_type)                       :: clusters
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ibox, jbox, kbox, &
                                                nclusters, nparticles

!

    nparticles = 0
    DO kbox = 0,UBOUND(boxes,3)-1
    DO jbox = 0,UBOUND(boxes,2)-1
    DO ibox = 0,UBOUND(boxes,1)-1
       IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
          nparticles = nparticles + boxes(ibox,jbox,kbox)%nparticles
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    nclusters = nparticles
    !
    CALL cluster_allocate(clusters, nclusters, nparticles, error)
    !
    !
    nclusters = 1
    nparticles = 0
    clusters%cluster_ptr(nclusters) = 1
    DO kbox = 0,UBOUND(boxes,3)-1
    DO jbox = 0,UBOUND(boxes,2)-1
    DO ibox = 0,UBOUND(boxes,1)-1
       IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
          DO iatom = 1,boxes(ibox,jbox,kbox)%nparticles
             nparticles = nparticles + 1
             clusters%particle_idx(nparticles) = boxes(ibox,jbox,kbox)%atom_list(iatom)
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
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'boxing_particles', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: box_size_guess = 8.0_dp

    INTEGER                                  :: handle, iatom, ibox, im, ip, &
                                                jbox, jm, jp, kbox, km, kp, &
                                                natms, nparticles, ntot, &
                                                ntot_in_cell
    INTEGER, DIMENSION(3)                    :: nbox, periodic
    LOGICAL                                  :: failure, orthorhombic, skin
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: xyz
    REAL(dp), DIMENSION(3)                   :: cell_center, dbox, r
    REAL(dp), DIMENSION(3, 3)                :: h
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    natms = SIZE(particle_set,1)
    !
    !
    !get_cell(cell,alpha,beta,gamma,deth,orthorhombic,abc,periodic,h,h_inv,id_nr)
    CALL get_cell(cell, orthorhombic=orthorhombic, periodic=periodic, h=h)
    !
    ! so far only for ortho and 3d perodic
    CPPostcondition(orthorhombic,cp_failure_level,routineP,error,failure)
    CPPostcondition(COUNT(periodic.EQ.1).EQ.3,cp_failure_level,routineP,error,failure)
    !
    ! box the atoms
    nbox(1) = CEILING(h(1,1) / box_size_guess)
    nbox(2) = CEILING(h(2,2) / box_size_guess)
    nbox(3) = CEILING(h(3,3) / box_size_guess)
    dbox(1) = h(1,1) / REAL(nbox(1),dp)
    dbox(2) = h(2,2) / REAL(nbox(2),dp)
    dbox(3) = h(3,3) / REAL(nbox(3),dp)

    IF(dbg) THEN
       WRITE(*,*) routineN//' natms=',natms
       WRITE(*,*) routineN//' h(:,1)=',h(:,1)
       WRITE(*,*) routineN//' h(:,2)=',h(:,2)
       WRITE(*,*) routineN//' h(:,3)=',h(:,3)
       WRITE(*,*) routineN//' nbox=',nbox
       WRITE(*,*) routineN//' dbox=',dbox
       WRITE(*,*) routineN//' box_size_guess=',box_size_guess
    ENDIF

    ALLOCATE(xyz(3,natms))
    cell_center(:) = 0.5_dp*( h(:,1) + h(:,2) + h(:,3) )
    DO iatom = 1,natms
       xyz(:,iatom) = pbc(particle_set(iatom)%r(:), cell) + cell_center(:)
       !write(*,*) routineN//' xyz=',iatom, xyz(:,iatom)
    ENDDO

    ALLOCATE(boxes(-1:nbox(1),-1:nbox(2),-1:nbox(3)))

    !
    ! First pass
    DO iatom = 1,natms
       ibox = FLOOR(xyz(1,iatom)/dbox(1))
       jbox = FLOOR(xyz(2,iatom)/dbox(2))
       kbox = FLOOR(xyz(3,iatom)/dbox(3))
       IF(ibox.EQ.-1.OR.ibox.EQ.nbox(1).OR.&
          jbox.EQ.-1.OR.jbox.EQ.nbox(2).OR.&
          kbox.EQ.-1.OR.kbox.EQ.nbox(3)) THEN
          WRITE(*,*) 'ibox,jbox,kbox',ibox,jbox,kbox,xyz(:,iatom)
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       !write(*,*) 'ibox,jbox,kbox',ibox,jbox,kbox,xyz(:,iatom)
       boxes(ibox,jbox,kbox)%nparticles = boxes(ibox,jbox,kbox)%nparticles + 1
       !
       ! The skin
       !ortho only
       im = FLOOR((xyz(1,iatom)-h(1,1))/dbox(1))
       ip = FLOOR((xyz(1,iatom)+h(1,1))/dbox(1))
       jm = FLOOR((xyz(2,iatom)-h(2,2))/dbox(2))
       jp = FLOOR((xyz(2,iatom)+h(2,2))/dbox(2))
       km = FLOOR((xyz(3,iatom)-h(3,3))/dbox(3))
       kp = FLOOR((xyz(3,iatom)+h(3,3))/dbox(3))

       !write(*,*) 'im,ip,jm,jp,km,kp',im,ip,jm,jp,km,kp
       !
       !!
       skin = .FALSE.
       IF(im.EQ.     -1) THEN
          ibox = im; skin = .TRUE.
       ELSEIF(ip.EQ.nbox(1)) THEN
          ibox = ip; skin = .TRUE.
       ENDIF
       IF(jm.EQ.     -1) THEN
          jbox = jm; skin = .TRUE.
       ELSEIF(jp.EQ.nbox(2)) THEN
          jbox = jp; skin = .TRUE.
       ENDIF
       IF(km.EQ.     -1) THEN
          kbox = km; skin = .TRUE.
       ELSEIF(kp.EQ.nbox(3)) THEN
          kbox = kp; skin = .TRUE.
       ENDIF
       !!
       IF(skin) THEN
          boxes(ibox,jbox,kbox)%nparticles = boxes(ibox,jbox,kbox)%nparticles + 1
       ENDIF
    ENDDO


    IF(.FALSE.) THEN
       WRITE(*,*) 'FIRST CHECK'
       DO kbox = -1,nbox(3)
       DO jbox = -1,nbox(2)
       DO ibox = -1,nbox(1)
          IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
             WRITE(*,*) 'box=',ibox,jbox,kbox
             WRITE(*,*) 'nbr atom=',boxes(ibox,jbox,kbox)%nparticles
             !DO iatom = 1,boxes(ibox,jbox,kbox)%nparticles
             !   WRITE(*,'(A,I3,3E14.6)') 'iatom=',boxes(ibox,jbox,kbox)%atom_list(iatom)
             !ENDDO
          ENDIF
       ENDDO
       ENDDO
       ENDDO
    ENDIF
    !
    ! allocate
    DO kbox = -1,nbox(3)
    DO jbox = -1,nbox(2)
    DO ibox = -1,nbox(1)
       nparticles = boxes(ibox,jbox,kbox)%nparticles
       IF(nparticles.GE.0) THEN
          ALLOCATE(boxes(ibox,jbox,kbox)%atom_list(nparticles))
          ALLOCATE(boxes(ibox,jbox,kbox)%xyz(3,nparticles))
       ENDIF
       ! reset to 0 for next round
       boxes(ibox,jbox,kbox)%nparticles = 0
    ENDDO
    ENDDO
    ENDDO
    !
    ! Second pass
    DO iatom = 1,natms
       ibox = FLOOR(xyz(1,iatom)/dbox(1))
       jbox = FLOOR(xyz(2,iatom)/dbox(2))
       kbox = FLOOR(xyz(3,iatom)/dbox(3))
       nparticles = boxes(ibox,jbox,kbox)%nparticles
       nparticles = nparticles + 1
       boxes(ibox,jbox,kbox)%nparticles = nparticles
       boxes(ibox,jbox,kbox)%atom_list(nparticles) = iatom
       boxes(ibox,jbox,kbox)%xyz(:,nparticles) = xyz(:,iatom)
       !
       ! The skin
       !ortho only
       im = FLOOR((xyz(1,iatom)-h(1,1))/dbox(1))
       ip = FLOOR((xyz(1,iatom)+h(1,1))/dbox(1))
       jm = FLOOR((xyz(2,iatom)-h(2,2))/dbox(2))
       jp = FLOOR((xyz(2,iatom)+h(2,2))/dbox(2))
       km = FLOOR((xyz(3,iatom)-h(3,3))/dbox(3))
       kp = FLOOR((xyz(3,iatom)+h(3,3))/dbox(3))
       !
       !!
       r(:) = 0.0_dp
       skin = .FALSE.
       IF(im.EQ.     -1) THEN
          ibox = im; r(1)=-h(1,1); skin = .TRUE.
       ELSEIF(ip.EQ.nbox(1)) THEN
          ibox = ip; r(1)=h(1,1); skin = .TRUE.
       ENDIF
       IF(jm.EQ.     -1) THEN
          jbox = jm; r(2)=-h(2,2); skin = .TRUE.
       ELSEIF(jp.EQ.nbox(2)) THEN
          jbox = jp; r(2)=h(2,2); skin = .TRUE.
       ENDIF
       IF(km.EQ.     -1) THEN
          kbox = km; r(3)=-h(3,3); skin = .TRUE.
       ELSEIF(kp.EQ.nbox(3)) THEN
          kbox = kp; r(3)=h(3,3); skin = .TRUE.
       ENDIF
       !!
       IF(skin) THEN
          nparticles = boxes(ibox,jbox,kbox)%nparticles
          nparticles = nparticles + 1
          boxes(ibox,jbox,kbox)%nparticles = nparticles
          boxes(ibox,jbox,kbox)%atom_list(nparticles) = iatom
          !write(*,*) 'r',r
          boxes(ibox,jbox,kbox)%xyz(:,nparticles) = xyz(:,iatom) + r(:)
       ENDIF
    ENDDO


    IF(.FALSE.) THEN
       WRITE(*,*) 'SECOND CHECK'
       ntot = 0
       ntot_in_cell = 0
       DO kbox = -1,nbox(3)
       DO jbox = -1,nbox(2)
       DO ibox = -1,nbox(1)
          IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
             WRITE(*,*) 'box=',ibox,jbox,kbox
             WRITE(*,*) 'nbr atom=',boxes(ibox,jbox,kbox)%nparticles
             ntot = ntot + boxes(ibox,jbox,kbox)%nparticles
             IF(ibox.GT.-1.AND.ibox.LT.nbox(1).AND.&
                jbox.GT.-1.AND.jbox.LT.nbox(2).AND.&
                kbox.GT.-1.AND.kbox.LT.nbox(3)) THEN
                ntot_in_cell = ntot_in_cell + boxes(ibox,jbox,kbox)%nparticles
             ENDIF
             DO iatom = 1,boxes(ibox,jbox,kbox)%nparticles
                WRITE(*,'(A,I3)') 'iatom=',boxes(ibox,jbox,kbox)%atom_list(iatom)
             ENDDO
          ENDIF
       ENDDO
       ENDDO
       ENDDO
       !
       WRITE(*,*) ntot_in_cell!ntot
       WRITE(*,*) 
       DO kbox = -1,nbox(3)
       DO jbox = -1,nbox(2)
       DO ibox = -1,nbox(1)
          IF(boxes(ibox,jbox,kbox)%nparticles.GT.0) THEN
             !WRITE(*,*) 'box=',ibox,jbox,kbox
             !WRITE(*,*) 'nbr atom=',boxes(ibox,jbox,kbox)%nparticles
             DO iatom = 1,boxes(ibox,jbox,kbox)%nparticles
                IF(ibox.EQ.-1.OR.ibox.EQ.nbox(1).OR.&
                   jbox.EQ.-1.OR.jbox.EQ.nbox(2).OR.&
                   kbox.EQ.-1.OR.kbox.EQ.nbox(3)) THEN
                   WRITE(*,'(A,3E14.6)') 'H',boxes(ibox,jbox,kbox)%xyz(:,iatom)*0.5_dp
                ELSE
                   IF(.FALSE.)WRITE(*,'(A,3E14.6)') 'X',boxes(ibox,jbox,kbox)%xyz(:,iatom)*0.5_dp
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       ENDDO
       ENDDO

    ENDIF

    DEALLOCATE(xyz)

    CALL timestop(handle)

  END SUBROUTINE boxing_particles

  SUBROUTINE box_deallocate(boxes, error)
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: boxes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: ibox, jbox, kbox

    IF(ASSOCIATED(boxes)) THEN
       DO kbox = LBOUND(boxes,3),UBOUND(boxes,3)
       DO jbox = LBOUND(boxes,2),UBOUND(boxes,2)
       DO ibox = LBOUND(boxes,1),UBOUND(boxes,1)
          IF(ASSOCIATED(boxes(ibox,jbox,kbox)%atom_list)) THEN
             DEALLOCATE(boxes(ibox,jbox,kbox)%atom_list)
          ENDIF
          IF(ASSOCIATED(boxes(ibox,jbox,kbox)%xyz)) THEN
             DEALLOCATE(boxes(ibox,jbox,kbox)%xyz)
          ENDIF
       ENDDO
       ENDDO
       ENDDO
       DEALLOCATE(boxes)
    ENDIF
  END SUBROUTINE box_deallocate

  SUBROUTINE cluster_allocate(clusters, nclusters, nparticles, error)
    TYPE(cluster_type)                       :: clusters
    INTEGER, INTENT(in)                      :: nclusters, nparticles
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    clusters%nclusters = nclusters
    clusters%nparticles = nparticles
    ALLOCATE(clusters%cluster_ptr(nclusters+1), clusters%particle_idx(nparticles))
  END SUBROUTINE cluster_allocate

  SUBROUTINE cluster_deallocate(clusters, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    DEALLOCATE(clusters%cluster_ptr, clusters%particle_idx)
  END SUBROUTINE cluster_deallocate

  SUBROUTINE cluster_print(clusters, error)
    TYPE(cluster_type)                       :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: icluster, ip, iparticle

    WRITE(*,'(/,/,T2,A)') 'CLUSTERING THE PARTICLES'
    WRITE(*,'(T5,A,I6,A)') 'Found ',clusters%nclusters,' clusters'
    WRITE(*,'(T5,A)') ' cluster   number of particles'

    DO icluster = 1,clusters%nclusters
       WRITE(*,'(T5,I6,T20,I6)') icluster,&
            clusters%cluster_ptr(icluster+1)-clusters%cluster_ptr(icluster)
       IF(.FALSE.) THEN
          DO ip = clusters%cluster_ptr(icluster),clusters%cluster_ptr(icluster+1)-1
             iparticle = clusters%particle_idx(ip)
             WRITE(*,*) 'ip',ip,' iparticle',iparticle
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE cluster_print

END MODULE cluster_methods
