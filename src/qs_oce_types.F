!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_oce_types
!!
!!   NAME
!!     qs_oce_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE qs_oce_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_node,&
                                             qlist_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE sparse_matrix_types,             ONLY: allocate_matrix,&
                                             allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_type, &
                                             real_matrix_p_type,&
                                             replicate_matrix_structure
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: locate,&
                                             sort

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_oce_types"

! *** Define a oce matrix type ***

  TYPE vtriple_type
    TYPE(qlist_type), DIMENSION(:,:,:), POINTER :: neighbor
  END TYPE vtriple_type


  TYPE oce_matrix_type
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER  :: hardo,softo
    TYPE(vtriple_type), DIMENSION(:,:), POINTER      :: voce 
    LOGICAL, DIMENSION(:,:), POINTER                 :: block_added
  END TYPE

! *** Public data types ***

  PUBLIC :: oce_matrix_type

! *** Public subroutines ***

  PUBLIC :: allocate_oce_set, &
            deallocate_oce_set, &
            deallocate_vtriple,&
            prepare_oce_coeff,&
            retrieve_sac_list,&
            vtriple_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_oce_set(oce_set,natom,nkind,ndim)

!   Purpose: Allocate and initialize the matrix set of oce coefficients.

!   History:

!   ***************************************************************************

    TYPE(oce_matrix_type), POINTER           :: oce_set
    INTEGER, INTENT(IN)                      :: natom, ndim, nkind

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE allocate_oce_set (MODULE oce_types)"

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ho,so
    INTEGER                                  :: ia, ik, iset, istat

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(oce_set)) CALL deallocate_oce_set(oce_set)

    ALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"oce_set",0)
 
    NULLIFY(ho,so)
!    CALL allocate_matrix_set(oce_set%hardo,ndim)
!    CALL allocate_matrix_set(oce_set%softo,ndim)
    CALL allocate_matrix_set(ho,ndim)
    CALL allocate_matrix_set(so,ndim)

    oce_set%hardo => ho
    oce_set%softo => so

    ALLOCATE (oce_set%voce(nkind,natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"voce",nkind*natom)
    DO ia = 1,natom
      DO ik = 1,nkind
        NULLIFY(oce_set%voce(ik,ia)%neighbor)
      ENDDO
    ENDDO

    ALLOCATE(oce_set%block_added(natom,natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"block_added",natom*natom)
    oce_set%block_added = .FALSE.
  
  END SUBROUTINE allocate_oce_set

! *****************************************************************************

  SUBROUTINE deallocate_oce_set(oce_set)

!   Purpose: Deallocate the matrix set of oce coefficients

! *****************************************************************************

    TYPE(oce_matrix_type), POINTER         :: oce_set

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE deallocate_oce_set (MODULE oce_types)"

    INTEGER                                :: istat
!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(oce_set)) RETURN

    IF(ASSOCIATED(oce_set%hardo)) CALL deallocate_matrix_set(oce_set%hardo)
    IF(ASSOCIATED(oce_set%softo)) CALL deallocate_matrix_set(oce_set%softo)
    IF(ASSOCIATED(oce_set%voce))  CALL deallocate_vtriple(oce_set%voce)
    IF(ASSOCIATED(oce_set%block_added)) DEALLOCATE(oce_set%block_added,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"deallocate_block_added",0)
    DEALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"deallocate_oce_set",0)

  END SUBROUTINE deallocate_oce_set

! *****************************************************************************

  SUBROUTINE deallocate_vtriple(vtriple)

    TYPE(vtriple_type), DIMENSION(:,:), POINTER :: vtriple

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_vtriple"

    TYPE(qlist_type), POINTER                   :: sac_oce_neighbor
    INTEGER  iat, icell, istat, jcell, kcell, kkind, natom, nkind

! *****************************************************************************

    IF(ASSOCIATED(vtriple)) THEN
      nkind = SIZE(vtriple,1)
      natom = SIZE(vtriple,2)
      DO iat = 1,natom
         DO kkind=1,nkind
           IF (ASSOCIATED(vtriple(kkind,iat)%neighbor)) THEN
             DO kcell=LBOUND(vtriple(kkind,iat)%neighbor,3),&
                      UBOUND(vtriple(kkind,iat)%neighbor,3)
               DO jcell=LBOUND(vtriple(kkind,iat)%neighbor,2),&
                        UBOUND(vtriple(kkind,iat)%neighbor,2)
                 DO icell=LBOUND(vtriple(kkind,iat)%neighbor,1),&
                          UBOUND(vtriple(kkind,iat)%neighbor,1)
                   sac_oce_neighbor => vtriple(kkind,iat)%neighbor(icell,&
                                                               jcell,&
                                                               kcell)
                   IF (sac_oce_neighbor%n > 0) THEN
                     DEALLOCATE (sac_oce_neighbor%r2,STAT=istat)
                     IF (istat /= 0) THEN
                       CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "sac_oce_neighbor%r2")
                     END IF
                     DEALLOCATE (sac_oce_neighbor%r,STAT=istat)
                     IF (istat /= 0) THEN
                       CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "sac_oce_neighbor%r")
                     END IF
                     DEALLOCATE (sac_oce_neighbor%list,STAT=istat)
                     IF (istat /= 0) THEN
                       CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "sac_oce_neighbor%list")
                     END IF
                     DEALLOCATE (sac_oce_neighbor%index_list,STAT=istat)
                     IF (istat /= 0) THEN
                       CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "sac_oce_neighbor%index_list")
                     END IF
                   END IF
                 END DO
               END DO
             END DO
             DEALLOCATE (vtriple(kkind,iat)%neighbor,STAT=istat)
             IF (istat /= 0) THEN
               CALL stop_memory(routine_name,module_name,__LINE__,&
                               "vtriple(kkind,iat)%neighbor")
             END IF
           END IF
         END DO
      END DO
      DEALLOCATE(vtriple,STAT=istat)
      IF (istat /= 0) THEN
                       CALL stop_memory(routine_name,module_name,__LINE__,&
                       "vtriple")
      END IF
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer vtriple is not associated")
    ENDIF

  END SUBROUTINE deallocate_vtriple
! *****************************************************************************

  SUBROUTINE prepare_oce_coeff(particle_set,oce,natom,nsgf, & 
                               first_sgf,last_sgf,distribution_2d)

!   Purpose: Create the structure for the matrix for oce coefficients

!   History: - Creation (5-02-04 MI)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(particle_type), DIMENSION(:), &
                         POINTER             :: particle_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    INTEGER, INTENT(IN)                      :: natom,nsgf
    INTEGER, DIMENSION(natom)                :: last_sgf , &
                                                first_sgf
                                     
    CHARACTER(LEN=*), PARAMETER :: routine_name = "prepare_oce_coeff"

    INTEGER, DIMENSION(:), ALLOCATABLE      :: last_nso,&
                                               first_nso
    CHARACTER(LEN=12)           :: cgfsym
    CHARACTER(LEN=80)           :: name
    INTEGER                     :: i,ikind, istat, maxder, nsotot

!   ---------------------------------------------------------------------------

    ALLOCATE (first_nso(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "first_nso",natom*int_size)
    ALLOCATE (last_nso(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "last_nso",natom*int_size)


    CALL get_particle_set(particle_set=particle_set, &
                          first_nso=first_nso,last_nso=last_nso,&
                          nsotot=nsotot)

    CALL allocate_matrix(matrix=oce%hardo(1)%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsotot,&
                         ncol=nsgf,&
                         first_row=first_nso,&
                         last_row=last_nso,&
                         first_col=first_sgf,&
                         last_col=last_sgf,&
                         matrix_name="HARD OCE COEFF",&
                         matrix_symmetry="non symmetric",&
                         distribution_2d=distribution_2d)

    DEALLOCATE (first_nso,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "first_nso")
    DEALLOCATE (last_nso,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "last_nso")

    name = "SOFT OCE COEFF"
    CALL compress(name)
    CALL uppercase(name)

    CALL replicate_matrix_structure(source=oce%hardo(1)%matrix,&
                                    target=oce%softo(1)%matrix,&
                                    target_name=TRIM(name),&
                                    target_symmetry="non symmetric")

    maxder = SIZE(oce%hardo,1)

    DO i=2,maxder
      cgfsym = cgf_symbol(1,indco(1:3,i))
      name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE HARD OCE COEFF "//&
             "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)

      CALL replicate_matrix_structure(source=oce%hardo(1)%matrix,&
                                      target=oce%hardo(i)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="non symmetric")
      
      name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE SOFT OCE COEFF "//&
             "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL replicate_matrix_structure(source=oce%softo(1)%matrix,&
                                      target=oce%softo(i)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="non symmetric")
    END DO

!   Initialize the control variable block_added
    oce%block_added(1:natom,1:natom) = .FALSE.

  END SUBROUTINE prepare_oce_coeff

! *****************************************************************************
  SUBROUTINE retrieve_sac_list(vlist,sac_list,iatom,ikind,nkind,globenv)

    TYPE(vtriple_type), DIMENSION(:,:), &
                     POINTER                 :: vlist
    TYPE(neighbor_list_set_p_type),  &
                  DIMENSION(:), POINTER      :: sac_list
    INTEGER, INTENT(IN)                      :: iatom,ikind,nkind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "retrieve_sac_list"

!   *** Local variables ***

    TYPE(qlist_type), POINTER            :: sac_list_neighbor
    TYPE(neighbor_list_type), POINTER    :: sac_list_neighbor_list
    TYPE(neighbor_node_type), POINTER    :: sac_list_neighbor_node

    REAL(dp)                             :: rac2
    REAL(dp), DIMENSION(3)               :: rac
    INTEGER, DIMENSION(3)                :: cell_c,cell_c_max,cell_c_min
    INTEGER :: i, iac, icell, istat, jcell, katom, kcell, kkind, kneighbor, n, nneighbor

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(sac_list)) THEN
    DO kkind=1,nkind
      iac = ikind + nkind*(kkind - 1)
      IF (.NOT.ASSOCIATED(sac_list(iac)%neighbor_list_set)) CYCLE

      sac_list_neighbor_list =>&
          find_neighbor_list(neighbor_list_set=&
                            sac_list(iac)%neighbor_list_set,&
                            atom=iatom)
      CALL get_neighbor_list(neighbor_list=sac_list_neighbor_list,&
                                     nnode=nneighbor)
!     *** Find the proper cell index ranges ***
      cell_c_max(:) = 0
      cell_c_min(:) = 0

      sac_list_neighbor_node => first_node(sac_list_neighbor_list)

      DO kneighbor=1,nneighbor
        CALL get_neighbor_node(neighbor_node=sac_list_neighbor_node,&
                               cell=cell_c)
        DO i=1,3
           cell_c_max(i) = MAX(cell_c_max(i),cell_c(i))
           cell_c_min(i) = MIN(cell_c_min(i),cell_c(i))
        END DO
        sac_list_neighbor_node => next(sac_list_neighbor_node)
      END DO


      ALLOCATE (vlist(kkind,iatom)%neighbor(cell_c_min(1):cell_c_max(1),&
                                  cell_c_min(2):cell_c_max(2),&
                                  cell_c_min(3):cell_c_max(3)),&
                                  STAT=istat)
      IF (istat /= 0 ) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,&
                         "vlist(kkind,iatom)%neighbor",&
                         (cell_c_max(1) - cell_c_min(1) + 1)*&
                         (cell_c_max(2) - cell_c_min(2) + 1)*&
                         (cell_c_max(3) - cell_c_min(3) + 1)*&
                          int_size)
       END IF
       vlist(kkind,iatom)%neighbor(:,:,:)%n = 0

!     *** Allocate and initialize the sac_list neighbor lists ***
!      *** Find proper array size for each cell ***

       sac_list_neighbor_node => first_node(sac_list_neighbor_list)

       DO kneighbor=1,nneighbor
          CALL get_neighbor_node(neighbor_node=sac_list_neighbor_node,&
                                cell=cell_c)
          vlist(kkind,iatom)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n =&
              vlist(kkind,iatom)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n + 1
              sac_list_neighbor_node => next(sac_list_neighbor_node)
       END DO

!      *** Allocate the sac_list neighbor list ***

       DO kcell=cell_c_min(3),cell_c_max(3)
          DO jcell=cell_c_min(2),cell_c_max(2)
             DO icell=cell_c_min(1),cell_c_max(1)
               sac_list_neighbor => vlist(kkind,iatom)%neighbor(icell,jcell,kcell)
               NULLIFY (sac_list_neighbor%sac)
               NULLIFY (sac_list_neighbor%r2)
               NULLIFY (sac_list_neighbor%r)
               NULLIFY (sac_list_neighbor%index_list)
               NULLIFY (sac_list_neighbor%list)
               n = sac_list_neighbor%n
               IF (n > 0) THEN
                 CALL reallocate(sac_list_neighbor%r2,1,n)
                 CALL reallocate(sac_list_neighbor%r,1,3,1,n)
                 CALL reallocate(sac_list_neighbor%index_list,1,n)
                 CALL reallocate(sac_list_neighbor%list,1,n)
               END IF
            END DO
         END DO
       END DO

!      *** Fill sac_list neighbor lists ***

      vlist(kkind,iatom)%neighbor(:,:,:)%n = 0

      sac_list_neighbor_node => first_node(sac_list_neighbor_list)

      DO kneighbor=1,nneighbor
        CALL get_neighbor_node(neighbor_node=sac_list_neighbor_node,&
                               neighbor=katom,&
                               cell=cell_c,&
                               r=rac)
         sac_list_neighbor => vlist(kkind,iatom)%neighbor(cell_c(1),&
                                                    cell_c(2),&
                                                    cell_c(3))
         sac_list_neighbor%n = sac_list_neighbor%n + 1
         sac_list_neighbor%list(sac_list_neighbor%n) = katom
         rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
         sac_list_neighbor%r2(sac_list_neighbor%n) = rac2
         sac_list_neighbor%r(:,sac_list_neighbor%n) = rac(:)
         sac_list_neighbor_node => next(sac_list_neighbor_node)

       END DO
       DO kcell=cell_c_min(3),cell_c_max(3)
         DO jcell=cell_c_min(2),cell_c_max(2)
            DO icell=cell_c_min(1),cell_c_max(1)

               sac_list_neighbor => vlist(kkind,iatom)%neighbor(icell,jcell,kcell)

!              *** Sort sac_list neighbor lists ***
               IF (sac_list_neighbor%n > 0) THEN
                  CALL sort(sac_list_neighbor%list,&
                          sac_list_neighbor%n,&
                          sac_list_neighbor%index_list)
               ELSE
                  CYCLE
               END IF
            END DO ! icell
         END DO ! jcell
       END DO ! kcell

    ENDDO !kkind
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
               "The pointer sac_list is not associated",&
                globenv)
    END IF
  END SUBROUTINE  retrieve_sac_list
! *****************************************************************************

END MODULE qs_oce_types
