!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
                                             list_type,&
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

  TYPE voce_type
    TYPE(list_type), DIMENSION(:,:,:), POINTER :: neighbor
  END TYPE voce_type


  TYPE oce_matrix_type
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER  :: hardo,softo
    TYPE(voce_type), DIMENSION(:,:), POINTER         :: voce 
  END TYPE

! *** Public data types ***

  PUBLIC :: oce_matrix_type

! *** Public subroutines ***

  PUBLIC :: allocate_oce_set, &
            deallocate_oce_set, &
            prepare_oce_coeff,&
            retrieve_sac_oce,&
            voce_type

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
    IF (istat /= 0) CALL stop_memory(routine,"voce",0)
    DO ia = 1,natom
      DO ik = 1,nkind
        NULLIFY(oce_set%voce(ik,ia)%neighbor)
      ENDDO
    ENDDO
  
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
    IF(ASSOCIATED(oce_set%voce))  CALL deallocate_voce(oce_set%voce)
    DEALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"deallocate_oce_set",0)

  END SUBROUTINE deallocate_oce_set

! *****************************************************************************

  SUBROUTINE deallocate_voce(voce)

    TYPE(voce_type), DIMENSION(:,:), POINTER :: voce

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_voce"

    TYPE(list_type), POINTER    :: sac_oce_neighbor
    INTEGER  iat, icell, istat, jcell, kcell, kkind, natom, nkind

! *****************************************************************************

    IF(ASSOCIATED(voce)) THEN
      nkind = SIZE(voce,1)
      natom = SIZE(voce,2)
      DO iat = 1,natom
         DO kkind=1,nkind
           IF (ASSOCIATED(voce(kkind,iat)%neighbor)) THEN
             DO kcell=LBOUND(voce(kkind,iat)%neighbor,3),&
                      UBOUND(voce(kkind,iat)%neighbor,3)
               DO jcell=LBOUND(voce(kkind,iat)%neighbor,2),&
                        UBOUND(voce(kkind,iat)%neighbor,2)
                 DO icell=LBOUND(voce(kkind,iat)%neighbor,1),&
                          UBOUND(voce(kkind,iat)%neighbor,1)
                   sac_oce_neighbor => voce(kkind,iat)%neighbor(icell,&
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
             DEALLOCATE (voce(kkind,iat)%neighbor,STAT=istat)
             IF (istat /= 0) THEN
               CALL stop_memory(routine_name,module_name,__LINE__,&
                               "voce(kkind,iat)%neighbor")
             END IF
           END IF
         END DO
       END DO
    ELSE

    ENDIF

  END SUBROUTINE deallocate_voce
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

  END SUBROUTINE prepare_oce_coeff

! *****************************************************************************
  SUBROUTINE retrieve_sac_oce(voce,sac_oce,iatom,ikind,nkind)

    TYPE(voce_type), DIMENSION(:,:), &
                     POINTER            :: voce
    TYPE(neighbor_list_set_p_type),  &
                  DIMENSION(:), POINTER :: sac_oce
    INTEGER, INTENT(IN)                  :: iatom,ikind,nkind

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "retrieve_sac_oce"

!   *** Local variables ***

    TYPE(list_type), POINTER            :: sac_oce_neighbor
    TYPE(neighbor_list_type), POINTER   :: sac_oce_neighbor_list
    TYPE(neighbor_node_type), POINTER   :: sac_oce_neighbor_node

    REAL(dp)                            :: rac2
    REAL(dp), DIMENSION(3)              :: rac
    INTEGER, DIMENSION(3)               :: cell_c,cell_c_max,cell_c_min
    INTEGER :: i, iac, icell, istat, jcell, katom, kcell, kkind, kneighbor, n, nneighbor

!   ---------------------------------------------------------------------------

    DO kkind=1,nkind
      iac = ikind + nkind*(kkind - 1)
      IF (.NOT.ASSOCIATED(sac_oce(iac)%neighbor_list_set)) CYCLE

      sac_oce_neighbor_list =>&
          find_neighbor_list(neighbor_list_set=&
                            sac_oce(iac)%neighbor_list_set,&
                            atom=iatom)
      CALL get_neighbor_list(neighbor_list=sac_oce_neighbor_list,&
                                     nnode=nneighbor)
!     *** Find the proper cell index ranges ***
      cell_c_max(:) = 0
      cell_c_min(:) = 0

      sac_oce_neighbor_node => first_node(sac_oce_neighbor_list)

      DO kneighbor=1,nneighbor
        CALL get_neighbor_node(neighbor_node=sac_oce_neighbor_node,&
                               cell=cell_c)
        DO i=1,3
           cell_c_max(i) = MAX(cell_c_max(i),cell_c(i))
           cell_c_min(i) = MIN(cell_c_min(i),cell_c(i))
        END DO
        sac_oce_neighbor_node => next(sac_oce_neighbor_node)
      END DO


      ALLOCATE (voce(kkind,iatom)%neighbor(cell_c_min(1):cell_c_max(1),&
                                  cell_c_min(2):cell_c_max(2),&
                                  cell_c_min(3):cell_c_max(3)),&
                                  STAT=istat)
      IF (istat /= 0 ) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,&
                         "voce(kkind,iatom)%neighbor",&
                         (cell_c_max(1) - cell_c_min(1) + 1)*&
                         (cell_c_max(2) - cell_c_min(2) + 1)*&
                         (cell_c_max(3) - cell_c_min(3) + 1)*&
                          int_size)
       END IF
       voce(kkind,iatom)%neighbor(:,:,:)%n = 0

!     *** Allocate and initialize the sac_oce neighbor lists ***
!      *** Find proper array size for each cell ***

       sac_oce_neighbor_node => first_node(sac_oce_neighbor_list)

       DO kneighbor=1,nneighbor
          CALL get_neighbor_node(neighbor_node=sac_oce_neighbor_node,&
                                cell=cell_c)
          voce(kkind,iatom)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n =&
              voce(kkind,iatom)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n + 1
              sac_oce_neighbor_node => next(sac_oce_neighbor_node)
       END DO

!      *** Allocate the sac_oce neighbor list ***

       DO kcell=cell_c_min(3),cell_c_max(3)
          DO jcell=cell_c_min(2),cell_c_max(2)
             DO icell=cell_c_min(1),cell_c_max(1)
               sac_oce_neighbor => voce(kkind,iatom)%neighbor(icell,jcell,kcell)
               NULLIFY (sac_oce_neighbor%sac)
               NULLIFY (sac_oce_neighbor%r2)
               NULLIFY (sac_oce_neighbor%r)
               NULLIFY (sac_oce_neighbor%index_list)
               NULLIFY (sac_oce_neighbor%list)
               n = sac_oce_neighbor%n
               IF (n > 0) THEN
                 CALL reallocate(sac_oce_neighbor%r2,1,n)
                 CALL reallocate(sac_oce_neighbor%r,1,3,1,n)
                 CALL reallocate(sac_oce_neighbor%index_list,1,n)
                 CALL reallocate(sac_oce_neighbor%list,1,n)
               END IF
            END DO
         END DO
       END DO

!      *** Fill sac_oce neighbor lists ***

      voce(kkind,iatom)%neighbor(:,:,:)%n = 0

      sac_oce_neighbor_node => first_node(sac_oce_neighbor_list)

      DO kneighbor=1,nneighbor
        CALL get_neighbor_node(neighbor_node=sac_oce_neighbor_node,&
                               neighbor=katom,&
                               cell=cell_c,&
                               r=rac)
         sac_oce_neighbor => voce(kkind,iatom)%neighbor(cell_c(1),&
                                                    cell_c(2),&
                                                    cell_c(3))
         sac_oce_neighbor%n = sac_oce_neighbor%n + 1
         sac_oce_neighbor%list(sac_oce_neighbor%n) = katom
         rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
         sac_oce_neighbor%r2(sac_oce_neighbor%n) = rac2
         sac_oce_neighbor%r(:,sac_oce_neighbor%n) = rac(:)
         sac_oce_neighbor_node => next(sac_oce_neighbor_node)

       END DO
       DO kcell=cell_c_min(3),cell_c_max(3)
         DO jcell=cell_c_min(2),cell_c_max(2)
            DO icell=cell_c_min(1),cell_c_max(1)

               sac_oce_neighbor => voce(kkind,iatom)%neighbor(icell,jcell,kcell)

!              *** Sort sac_oce neighbor lists ***
               IF (sac_oce_neighbor%n > 0) THEN
                  CALL sort(sac_oce_neighbor%list,&
                          sac_oce_neighbor%n,&
                          sac_oce_neighbor%index_list)
               ELSE
                  CYCLE
               END IF
            END DO ! icell
         END DO ! jcell
       END DO ! kcell

    ENDDO !kkind
  END SUBROUTINE  retrieve_sac_oce
! *****************************************************************************

END MODULE qs_oce_types
