!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_oce_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE qs_neighbor_list_types,          ONLY: &
       find_neighbor_list, first_node, get_neighbor_list, get_neighbor_node, &
       neighbor_list_set_p_type, neighbor_list_type, neighbor_node_type, &
       next, qlist_type, reduced_3c_list_type
  USE sap_kind_types,                  ONLY: release_sap_int,&
                                             sap_int_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate,&
                                             sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oce_types'

! *** Define a oce matrix type ***

! *****************************************************************************
  TYPE vtriple_type
    TYPE(qlist_type), DIMENSION(:,:,:), POINTER :: neighbor
  END TYPE vtriple_type

! *****************************************************************************
  TYPE oce_matrix_type
    TYPE(sap_int_type), DIMENSION(:), POINTER        :: intac
  END TYPE

! *** Public data types ***

!  PUBLIC :: dist_list_4oce_type,oce_couple,oce_matrix_type,vtriple_type
  PUBLIC :: oce_matrix_type,vtriple_type

! *** Public subroutines ***

  PUBLIC :: allocate_oce_set, &
            allocate_vtriple,&
            deallocate_vtriple,&
            create_oce_set,&
            deallocate_oce_set,&
            build_reduced_3c_lists,&
            retrieve_sac_list

CONTAINS

! *****************************************************************************
!> \brief   Allocate and initialize the matrix set of oce coefficients.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_oce_set(oce_set,natom,nkind,ndim,error)
    TYPE(oce_matrix_type), POINTER           :: oce_set
    INTEGER, INTENT(IN)                      :: natom, nkind, ndim
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_oce_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat
    LOGICAL                                  :: failure = .FALSE.

    ALLOCATE(oce_set%intac(nkind*nkind),STAT=istat)
    CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
    DO i=1,nkind*nkind
      NULLIFY(oce_set%intac(i)%alist)
    END DO

  END SUBROUTINE allocate_oce_set

! *****************************************************************************
  SUBROUTINE allocate_vtriple(vtriple,nkind,natom,error)

    TYPE(vtriple_type), DIMENSION(:, :), &
      POINTER                                :: vtriple
    INTEGER, INTENT(IN)                      :: nkind, natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_vtriple', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, ikind, istat

    IF(ASSOCIATED(vtriple)) CALL deallocate_vtriple(vtriple,error=error)

    ALLOCATE (vtriple(nkind,natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"vtriple",nkind*natom)
    DO iat = 1,natom
      DO ikind = 1,nkind
        NULLIFY(vtriple(ikind,iat)%neighbor)
      ENDDO
    ENDDO

  END SUBROUTINE allocate_vtriple

! *****************************************************************************
  SUBROUTINE create_oce_set(oce_set,error)

    TYPE(oce_matrix_type), POINTER           :: oce_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_oce_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF(ASSOCIATED(oce_set)) CALL deallocate_oce_set(oce_set,error=error)

    ALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"oce_set",0)

    NULLIFY(oce_set%intac)

  END SUBROUTINE create_oce_set

! *****************************************************************************
!> \brief  Deallocate the matrix set of oce coefficients 
!> \author  
!> \date    
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_oce_set(oce_set,error)
    TYPE(oce_matrix_type), POINTER           :: oce_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_oce_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (.NOT.ASSOCIATED(oce_set)) RETURN

    IF(ASSOCIATED(oce_set%intac)) CALL release_sap_int(oce_set%intac,error=error)

    DEALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"deallocate_oce_set",0)

  END SUBROUTINE deallocate_oce_set

! *****************************************************************************
  SUBROUTINE deallocate_vtriple(vtriple,error)

    TYPE(vtriple_type), DIMENSION(:, :), &
      POINTER                                :: vtriple
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_vtriple', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, icell, istat, jcell, &
                                                kcell, kkind, natom, nkind
    LOGICAL                                  :: failure
    TYPE(qlist_type), POINTER                :: sac_oce_neighbor

    failure = .FALSE.

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
                     CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
                     DEALLOCATE (sac_oce_neighbor%r,STAT=istat)
                     CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
                     DEALLOCATE (sac_oce_neighbor%list,STAT=istat)
                     CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
                     DEALLOCATE (sac_oce_neighbor%index_list,STAT=istat)
                     CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
                   END IF
                 END DO
               END DO
             END DO
             DEALLOCATE (vtriple(kkind,iat)%neighbor,STAT=istat)
             CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
           END IF
         END DO
      END DO
      DEALLOCATE(vtriple,STAT=istat)
      CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
    ELSE
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE deallocate_vtriple

! *****************************************************************************
  SUBROUTINE build_reduced_3c_lists(atomic_kind_set,iatom,jatom,ikind,jkind, &
                         sbc_list, vlist,reduced_3c_list,pippo,eps_rho,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: iatom, jatom, ikind, jkind
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sbc_list
    TYPE(vtriple_type), DIMENSION(:, :), &
      POINTER                                :: vlist
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_list
    LOGICAL                                  :: pippo
    REAL(dp), INTENT(IN)                     :: eps_rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_reduced_3c_lists', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ibc, iblock, istat, katom, kkind, kkneighbor, &
      kneighbor, nblock_kkind, nkind, nneighbor, nneighbor_old
    INTEGER, DIMENSION(3)                    :: cell_c
    INTEGER, DIMENSION(:), POINTER           :: index_atom
    LOGICAL                                  :: add_node, failure, paw_atom
    REAL(KIND=dp)                            :: dac2, dbc2, maxab
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rac2_tmp, rbc2_tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rac_tmp, rbc_tmp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(neighbor_list_type), POINTER        :: sbc_list_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sbc_list_neighbor_node
    TYPE(qlist_type), DIMENSION(:, :, :), &
      POINTER                                :: tmp_neigh
    TYPE(qlist_type), POINTER                :: sac_list_neighbor

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY(index_atom, rac_tmp, rac2_tmp, rbc_tmp, rbc2_tmp)
    nneighbor_old = 0
    nkind = SIZE(atomic_kind_set,1)

    DO kkind=1,nkind
       nblock_kkind = 0

       IF(.NOT.ASSOCIATED(vlist(kkind,iatom)%neighbor)) CYCLE

       ibc = jkind + nkind*(kkind - 1)

       IF (.NOT.ASSOCIATED(sbc_list(ibc)%neighbor_list_set)) CYCLE
       NULLIFY(sbc_list_neighbor_list)

       sbc_list_neighbor_list =>&
             find_neighbor_list(neighbor_list_set=&
             sbc_list(ibc)%neighbor_list_set,&
             atom=jatom)

       IF (.NOT.ASSOCIATED(sbc_list_neighbor_list)) CYCLE

       atomic_kind => atomic_kind_set(kkind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            paw_atom=paw_atom)

       IF(.NOT. paw_atom .AND. pippo)  CYCLE

       CALL get_neighbor_list(neighbor_list=sbc_list_neighbor_list,&
                                     nnode=nneighbor)

       sbc_list_neighbor_node => first_node(sbc_list_neighbor_list)

       ! Allocate temporary arrays
       IF(nneighbor > 0) THEN
         ALLOCATE(rac_tmp(3,nneighbor),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(rac2_tmp(nneighbor),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(rbc_tmp(3,nneighbor),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(rbc2_tmp(nneighbor),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(index_atom(nneighbor),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       knodeloop: DO kneighbor=1,nneighbor
          CALL get_neighbor_node(neighbor_node=sbc_list_neighbor_node,&
                            neighbor=katom,&
                            cell=cell_c,&
                            r=rbc)
          dbc2 = rbc(1)*rbc(1)+ rbc(2)*rbc(2) + rbc(3)*rbc(3)

          tmp_neigh => vlist(kkind,iatom)%neighbor
          DO i=1, 3
            IF ( (cell_c(i) < LBOUND(tmp_neigh,i)) .OR. &
                 (cell_c(i) > UBOUND(tmp_neigh,i))) THEN
              sbc_list_neighbor_node => next(sbc_list_neighbor_node)
              CYCLE knodeloop
            END IF
          END DO

          sac_list_neighbor => vlist(kkind,iatom)%neighbor(cell_c(1),&
                                                           cell_c(2),&
                                                           cell_c(3))

          IF (sac_list_neighbor%n == 0) THEN
              sbc_list_neighbor_node => next(sbc_list_neighbor_node)
              CYCLE
          END IF

!         *** Locate operator atom in the sac_oce neighbor list ***

           kkneighbor = locate(sac_list_neighbor%list,katom)

           IF (kkneighbor == 0) THEN
               sbc_list_neighbor_node => next(sbc_list_neighbor_node)
               CYCLE
           END IF

           rac(:) = sac_list_neighbor%r(:,kkneighbor)
           dac2 = rac(1)*rac(1)+ rac(2)*rac(2) + rac(3)*rac(3)

!  Check whether this node has to be really added to the list
           add_node = .TRUE.
           IF(pippo) THEN
              CALL check_oce_overlap(atomic_kind_set,ikind,rac,dac2,&
                   jkind,rbc,dbc2,kkind,maxab)
              IF(maxab .LT. eps_rho) THEN
                add_node = .FALSE.
              END IF
           END IF
           IF(.NOT. add_node) THEN
             sbc_list_neighbor_node => next(sbc_list_neighbor_node)
             CYCLE
           END IF

           ! add one block to the reduced list and fill in data
           nblock_kkind = nblock_kkind + 1

           rac_tmp(1:3,nblock_kkind) = rac(1:3)
           rac2_tmp(nblock_kkind) = dac2
           rbc_tmp(1:3,nblock_kkind) = rbc(1:3)
           rbc2_tmp(nblock_kkind) = dbc2
           index_atom(nblock_kkind) = katom

           sbc_list_neighbor_node => next(sbc_list_neighbor_node)

       END DO  knodeloop ! kneighbor

       ! Copy the reduced list data in the reduced list array
       reduced_3c_list(kkind)%nnode = nblock_kkind
       IF(nblock_kkind .GT. 0) THEN
         ALLOCATE(reduced_3c_list(kkind)%rac(3,nblock_kkind),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(reduced_3c_list(kkind)%rac2(nblock_kkind),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(reduced_3c_list(kkind)%rbc(3,nblock_kkind),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(reduced_3c_list(kkind)%rbc2(nblock_kkind),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(reduced_3c_list(kkind)%index_atom(nblock_kkind),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         DO iblock = 1,nblock_kkind
           reduced_3c_list(kkind)%rac(1:3,iblock)=rac_tmp(1:3,iblock)
           reduced_3c_list(kkind)%rac2(iblock)=rac2_tmp(iblock)
           reduced_3c_list(kkind)%rbc(1:3,iblock)=rbc_tmp(1:3,iblock)
           reduced_3c_list(kkind)%rbc2(iblock)=rbc2_tmp(iblock)
           reduced_3c_list(kkind)%index_atom(iblock)=index_atom(iblock)
         END DO
       END IF

       IF(nneighbor > 0) THEN
         DEALLOCATE(index_atom, rac_tmp, rac2_tmp, rbc_tmp, rbc2_tmp,STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
    ENDDO  ! kkind

    CALL timestop(handle)

  END SUBROUTINE build_reduced_3c_lists

! *****************************************************************************
  SUBROUTINE check_oce_overlap(atomic_kind_set,ikind,rac,dac2,jkind,rbc,dbc2,&
             kkind,maxab)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: ikind
    REAL(dp), INTENT(IN)                     :: rac(3)
    REAL(dp)                                 :: dac2
    INTEGER, INTENT(IN)                      :: jkind
    REAL(dp), INTENT(IN)                     :: rbc(3)
    REAL(dp)                                 :: dbc2
    INTEGER, INTENT(IN)                      :: kkind
    REAL(dp), INTENT(INOUT)                  :: maxab

    CHARACTER(len=*), PARAMETER :: routineN = 'check_oce_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iexp, iexp_zmin, iset, &
                                                ishell, jexp_zmin, jset, &
                                                jshell, lc_max, lshella, &
                                                lshellb, nseta, nsetb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nprjc, &
                                                nshella, nshellb
    INTEGER, DIMENSION(:, :), POINTER        :: lsha, lshb
    REAL(dp) :: gauss_value_a, gauss_value_b, gcca_zmin, gccb_zmin, &
      max_overlap_ab, overlap_ab, rcprj, zetamin, zetbmin
    REAL(dp), DIMENSION(:, :), POINTER       :: zeta, zetb, zetc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcca, gccb
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_c

    NULLIFY(atomic_kind,orb_basis_set,paw_proj_c)
    NULLIFY(gcca,gccb,la_max,lb_max,npgfa,npgfb,nprjc,zeta,zetb,zetc)

! ikind
    atomic_kind => atomic_kind_set(ikind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
         orb_basis_set=orb_basis_set)
    CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
         lmax=la_max,lmin=la_min,l=lsha,npgf=npgfa,&
         nshell=nshella,nset=nseta,gcc=gcca,zet=zeta)

! jkind
    atomic_kind => atomic_kind_set(jkind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
         orb_basis_set=orb_basis_set)
    CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
         lmax=lb_max,lmin=lb_min,l=lshb,nshell=nshellb,&
         npgf=npgfb,nset=nsetb,gcc=gccb,zet=zetb)

! kkind
    atomic_kind => atomic_kind_set(kkind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
         paw_proj_set=paw_proj_c)
    CALL get_paw_proj_set(paw_proj_set=paw_proj_c,&
         maxl=lc_max, nprj=nprjc, zetprj=zetc,rcprj=rcprj)

    max_overlap_ab  = 0.0_dp
    DO iset = 1,nseta
       zetamin = HUGE(1.0_dp)
       gcca_zmin = 1.0_dp
       DO iexp = 1,npgfa(iset)
         IF(zeta(iexp,iset) .LT. zetamin) THEN
           zetamin = MIN(zetamin,zeta(iexp,iset))
           iexp_zmin = iexp
         END IF
       END DO
       DO ishell = 1,nshella(iset)
         lshella = lsha(ishell,iset)
         gcca_zmin = gcca(iexp_zmin,ishell,iset)
!         gauss_value_a = (sqrt(dac2)- rcprj)**(lshella)*exp(-zetamin*(sqrt(dac2)- rcprj)**2)
         IF(dac2<1.0E-8_dp ) THEN
           gauss_value_a = 1.0_dp
         ELSE
           gauss_value_a = gcca_zmin*(SQRT(dac2))**(lshella)*EXP(-zetamin*(SQRT(dac2))**2)
         END IF
!      write(*,*) ' AC ',sqrt(dac2), iset,zetamin,lshella,gcca_zmin,gauss_value_a
         DO jset = 1,nsetb
           zetbmin = HUGE(1.0_dp)
           gccb_zmin = 1.0_dp
           DO iexp = 1,npgfb(jset)
             IF(zetb(iexp,jset) .LT. zetbmin) THEN
               zetbmin = zetb(iexp,jset)
               jexp_zmin = iexp
             END IF
           END DO
           DO jshell = 1,nshellb(jset)
             lshellb = lshb(jshell,jset)
             gccb_zmin = gccb(jexp_zmin,jshell,jset)

!             gauss_value_b = (sqrt(dbc2)- rcprj)**(lshellb)*exp(-zetbmin*(sqrt(dbc2)- rcprj)**2)
             IF(dbc2<1.0E-8_dp ) THEN
               gauss_value_b = 1.0_dp
             ELSE
               gauss_value_b = gccb_zmin*(SQRT(dbc2))**(lshellb)*EXP(-zetbmin*(SQRT(dbc2))**2)
             END IF

!      write(*,*) ' BC ',sqrt(dbc2), jset,zetbmin,lshellb,gccb_zmin,gauss_value_b

             overlap_ab = gauss_value_b*gauss_value_a
             max_overlap_ab = MAX(overlap_ab,max_overlap_ab)

!             DO kset = 0,lc_max
!               zetcmin= 100.0_dp
!               DO ip = 1,nprjc(kset)
!                  zetcmin = MIN(zetcmin,zetc(ip,kset))
!               END DO
!               gauss_value_c = (rcprj**kset)*exp(-zetcmin*rcprj**2)
!
!               overlap_abc = overlap_ab*gauss_value_c
!               max_overlap_abc = MAX(overlap_abc,max_overlap_abc)
!
!             END DO ! kset

           END DO  ! lshellb
         END DO  ! jset
       END DO  ! lshella
    END DO  ! iset

    maxab = max_overlap_ab

  END SUBROUTINE
!
! *****************************************************************************
  SUBROUTINE retrieve_sac_list(vlist,sac_list,iatom,ikind,nkind,error)

    TYPE(vtriple_type), DIMENSION(:, :), &
      POINTER                                :: vlist
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sac_list
    INTEGER, INTENT(IN)                      :: iatom, ikind, nkind
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'retrieve_sac_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iac, icell, istat, &
                                                jcell, katom, kcell, kkind, &
                                                kneighbor, n, nneighbor
    INTEGER, DIMENSION(3)                    :: cell_c, cell_c_max, cell_c_min
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rac2
    REAL(dp), DIMENSION(3)                   :: rac
    TYPE(neighbor_list_type), POINTER        :: sac_list_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sac_list_neighbor_node
    TYPE(qlist_type), POINTER                :: sac_list_neighbor

    CALL timeset(routineN,handle)

    failure = .FALSE.

    IF (ASSOCIATED(sac_list)) THEN
      DO kkind=1,nkind
         iac = ikind + nkind*(kkind - 1)
         IF (.NOT.ASSOCIATED(sac_list(iac)%neighbor_list_set)) CYCLE

         sac_list_neighbor_list =>&
            find_neighbor_list(neighbor_list_set=&
                              sac_list(iac)%neighbor_list_set,&
                              atom=iatom)
         IF (.NOT.ASSOCIATED(sac_list_neighbor_list)) CYCLE
         CALL get_neighbor_list(neighbor_list=sac_list_neighbor_list,&
                                     nnode=nneighbor)
!        *** Find the proper cell index ranges ***
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
           CALL stop_memory(routineN,moduleN,__LINE__,&
                           "vlist(kkind,iatom)%neighbor",&
                           (cell_c_max(1) - cell_c_min(1) + 1)*&
                           (cell_c_max(2) - cell_c_min(2) + 1)*&
                           (cell_c_max(3) - cell_c_min(3) + 1)*&
                            int_size)
         END IF
         vlist(kkind,iatom)%neighbor(:,:,:)%n = 0

!        *** Allocate and initialize the sac_list neighbor lists ***
!        *** Find proper array size for each cell ***

         sac_list_neighbor_node => first_node(sac_list_neighbor_list)

         DO kneighbor=1,nneighbor
           CALL get_neighbor_node(neighbor_node=sac_list_neighbor_node,&
                                 cell=cell_c)
           vlist(kkind,iatom)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n =&
              vlist(kkind,iatom)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n + 1
              sac_list_neighbor_node => next(sac_list_neighbor_node)
         END DO

!        *** Allocate the sac_list neighbor list ***

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

!        *** Fill sac_list neighbor lists ***

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

     CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE  retrieve_sac_list

END MODULE qs_oce_types
