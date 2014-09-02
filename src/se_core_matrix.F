!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Hamiltonian integral matrix <a|H|b> for
!>      semi-empirical methods
!> \author JGH
! *****************************************************************************
MODULE se_core_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_distribute, &
       cp_dbcsr_get_block_diag, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_replicate_all, cp_dbcsr_set, &
       cp_dbcsr_sum_replicated, cp_dbcsr_type
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_mndod, do_method_pdg, &
       do_method_pm3, do_method_pm6, do_method_pnnl, do_method_rm1, &
       use_orb_basis_set
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_int_arrays,       ONLY: rij_threshold
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_core_matrix'

  PUBLIC :: build_se_core_matrix

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param calculate_forces ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_se_core_matrix(qs_env,para_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_se_core_matrix', &
      routineP = moduleN//':'//routineN
    INTEGER, DIMENSION(4), PARAMETER         :: ix = (/1,3,4,2/)

    INTEGER :: atom_a, atom_b, handle, i, iatom, icol, icor, ikind, inode, &
      irow, istat, itype, iw, j, jatom, jkind, natom, natorb_a, nkind, nr_a, &
      nra, nrb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nrt
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, found, &
                                                use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: se_defined
    REAL(KIND=dp)                            :: delta, dr, econst, eheat, &
                                                eisol, kh, udd, uff, upp, &
                                                uss, ZPa, ZPb, ZSa, ZSb
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ZPt, ZSt
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: hmt, umt
    REAL(KIND=dp), DIMENSION(16)             :: ha, hb, ua
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(:), POINTER     :: beta_a, sto_exponents_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsmat, h_block, h_blocka, &
                                                pabmat, pamat, s_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: diagmat_h, diagmat_p
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(semi_empirical_type), POINTER       :: se_kind_a
    TYPE(virial_type), POINTER               :: virial

!   REAL(KIND=dp), DIMENSION(3)              :: R

    CALL timeset(routineN,handle)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

     NULLIFY ( rho,force,atomic_kind_set,qs_kind_set,sab_orb,&
               diagmat_h,diagmat_p,particle_set, matrix_p )

     ! calculate overlap matrix
     IF(calculate_forces) THEN
        CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, error=error)
        CALL build_overlap_matrix(qs_env,nderivative=1, matrix_s=matrix_s,&
                              matrix_name="OVERLAP",&
                              basis_set_id_a=use_orb_basis_set,&
                              basis_set_id_b=use_orb_basis_set, &
                              sab_nl=qs_env%sab_orb,&
                              error=error)
        CALL set_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
        CALL get_qs_env(qs_env=qs_env, virial=virial, error=error)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
     ELSE
        CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, error=error)
        CALL build_overlap_matrix(qs_env,matrix_s=matrix_s,&
                              matrix_name="OVERLAP",&
                              basis_set_id_a=use_orb_basis_set,&
                              basis_set_id_b=use_orb_basis_set, &
                              sab_nl=qs_env%sab_orb,&
                              error=error)
        CALL set_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
        use_virial = .FALSE.
     END IF

     CALL get_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h,&
                     matrix_s=matrix_s,&
                     atomic_kind_set=atomic_kind_set,&
                     qs_kind_set=qs_kind_set,&
                     dft_control=dft_control,error=error)

     IF(calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       particle_set=particle_set,&
                       rho=rho,&
                       force=force,error=error)
       matrix_p => rho%rho_ao

       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, alpha_scalar=1.0_dp, beta_scalar=1.0_dp, error=error)
       END IF
       natom = SIZE (particle_set)
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       delta = dft_control%qs_control%se_control%delta
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind)
       ALLOCATE(diagmat_p,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_dbcsr_init (diagmat_p, error=error)
       CALL cp_dbcsr_get_block_diag(matrix_p(1)%matrix, diagmat_p, error=error)
       CALL cp_dbcsr_replicate_all(diagmat_p,error=error)
     END IF

     ! Allocate the core Hamiltonian matrix
     CALL cp_dbcsr_allocate_matrix_set(matrix_h,1,error=error)
     ALLOCATE(matrix_h(1)%matrix,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL cp_dbcsr_init(matrix_h(1)%matrix, error=error)
     CALL cp_dbcsr_copy(matrix_h(1)%matrix,matrix_s(1)%matrix,"CORE HAMILTONIAN MATRIX",error=error)
     CALL cp_dbcsr_set(matrix_h(1)%matrix,0.0_dp,error=error)

     ! Allocate a diagonal block matrix
     ALLOCATE(diagmat_h,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL cp_dbcsr_init(diagmat_h, error=error)
     CALL cp_dbcsr_get_block_diag(matrix_s(1)%matrix, diagmat_h, error=error)
     CALL cp_dbcsr_set(diagmat_h, 0.0_dp, error=error)
     CALL cp_dbcsr_replicate_all(diagmat_h,error=error)

     CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb,error=error)

     ! kh might be set in qs_control
     itype = get_se_type(dft_control%qs_control%method_id)
     kh = 0.5_dp

     nkind = SIZE(atomic_kind_set)

     ALLOCATE (se_defined(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (hmt(16,nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (umt(16,nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     ALLOCATE (ZSt(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (ZPt(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (nrt(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     econst   = 0.0_dp
     DO ikind=1,nkind
        CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
        CALL get_qs_kind(qs_kind_set(ikind), se_parameter=se_kind_a)
        CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a,&
             beta=beta_a,uss=uss,upp=upp,udd=udd,uff=uff,eisol=eisol,eheat=eheat,&
             nr=nr_a,sto_exponents=sto_exponents_a)
        econst=econst-(eisol-eheat)*REAL(natom,dp)
        se_defined(ikind) = (defined .AND. natorb_a >= 1)
        hmt( 1   ,ikind) = beta_a(0)
        hmt( 2: 4,ikind) = beta_a(1)
        hmt( 5: 9,ikind) = beta_a(2)
        hmt(10:16,ikind) = beta_a(3)
        umt( 1   ,ikind) = uss
        umt( 2: 4,ikind) = upp
        umt( 5: 9,ikind) = udd
        umt(10:16,ikind) = uff

        ZSt(ikind) = sto_exponents_a(0)
        ZPt(ikind) = sto_exponents_a(1)
        nrt(ikind) = nr_a

     END DO
     qs_env%energy%core_self    = econst

     CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,inode=inode,r=rij)
        IF (.NOT.se_defined(ikind)) CYCLE
        IF (.NOT.se_defined(jkind)) CYCLE
        ha(1:16) = hmt(1:16,ikind)
        ua(1:16) = umt(1:16,ikind)
        hb(1:16) = hmt(1:16,jkind)

        nra = nrt(ikind)
        nrb = nrt(jkind)
        ZSa = ZSt(ikind)
        ZSb = ZSt(jkind)
        ZPa = ZPt(ikind)
        ZPb = ZPt(jkind)

        IF(inode==1) THEN
          SELECT CASE (dft_control%qs_control%method_id)
            CASE (do_method_am1,do_method_rm1,do_method_mndo,do_method_pdg,&
                  do_method_pm3,do_method_pm6,do_method_mndod, do_method_pnnl)
              NULLIFY(h_blocka)
              CALL cp_dbcsr_get_block_p(diagmat_h,iatom,iatom,h_blocka,found)
              CPPostcondition(ASSOCIATED(h_blocka),cp_failure_level,routineP,error,failure)
              IF(calculate_forces) THEN
                CALL cp_dbcsr_get_block_p(diagmat_p,iatom,iatom,pamat,found)
                CPPostcondition(ASSOCIATED(pamat),cp_failure_level,routineP,error,failure)
              END IF
          END SELECT
        END IF
        dr = SUM(rij(:)**2)
        IF (iatom == jatom .AND. dr < rij_threshold) THEN

          SELECT CASE (dft_control%qs_control%method_id)
             CASE DEFAULT
               CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             CASE (do_method_am1,do_method_rm1,do_method_mndo,do_method_pdg,&
                   do_method_pm3,do_method_pm6,do_method_mndod, do_method_pnnl)
               DO i=1,SIZE(h_blocka,1)
                 h_blocka(i,i)=h_blocka(i,i)+ua(i)
               END DO
          END SELECT

        ELSE
          IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
          ELSE
            irow = jatom
            icol = iatom
          END IF
          NULLIFY (h_block)
          CALL cp_dbcsr_get_block_p(matrix_h(1)%matrix,&
               irow,icol,h_block,found)
          CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)
          ! two-centre one-electron term
          NULLIFY(s_block)

!          CALL cp_dbcsr_get_block_p(matrix_s(1)%matrix,&
!               irow,icol,s_block,found)
!          CPPostcondition(ASSOCIATED(s_block),cp_failure_level,routineP,error,failure)
!
!          if( irow == iatom )then
!            R= -rij
!            call makeS(R,nra,nrb,ZSa,ZSb,ZPa,ZPb,S)
!          else
!            R= rij
!            call makeS(R,nrb,nra,ZSb,ZSa,ZPb,ZPa,S)
!          endif
!
!          do i=1,4
!            do j=1,4
!              s_block(i,j)=S(ix(i),ix(j))
!            enddo
!          enddo

          CALL cp_dbcsr_get_block_p(matrix_s(1)%matrix,&
               irow,icol,s_block,found)
          CPPostcondition(ASSOCIATED(s_block),cp_failure_level,routineP,error,failure)
          IF ( irow == iatom ) THEN
            DO i=1,SIZE(h_block,1)
              DO j=1,SIZE(h_block,2)
                h_block(i,j)=h_block(i,j)+kh*(ha(i)+hb(j))*s_block(i,j)
              END DO
            END DO
          ELSE
            DO i=1,SIZE(h_block,1)
              DO j=1,SIZE(h_block,2)
                h_block(i,j)=h_block(i,j)+kh*(ha(j)+hb(i))*s_block(i,j)
              END DO
            END DO
          END IF
          IF(calculate_forces) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)

!            if( irow == iatom )then
!              R= -rij
!              call makedS(R,nra,nrb,ZSa,ZSb,ZPa,ZPb,dS)
!            else
!              R= rij
!              call makedS(R,nrb,nra,ZSb,ZSa,ZPb,ZPa,dS)
!            endif

            CALL cp_dbcsr_get_block_p(matrix_p(1)%matrix,irow,icol,pabmat,found)
            CPPostcondition(ASSOCIATED(pabmat),cp_failure_level,routineP,error,failure)
            DO icor=1,3
              force_ab(icor) = 0._dp

!              CALL cp_dbcsr_get_block_p(matrix_s(icor+1)%matrix,irow,icol,dsmat,found)
!              CPPostcondition(ASSOCIATED(dsmat),cp_failure_level,routineP,error,failure)
!
!              do i=1,4
!                do j=1,4
!                  dsmat(i,j)=dS(ix(i),ix(j),icor)
!                enddo
!              enddo

              CALL cp_dbcsr_get_block_p(matrix_s(icor+1)%matrix,irow,icol,dsmat,found)
              CPPostcondition(ASSOCIATED(dsmat),cp_failure_level,routineP,error,failure)
              dsmat=2._dp*kh*dsmat*pabmat
              IF ( irow == iatom ) THEN
                 DO i=1,SIZE(h_block,1)
                   DO j=1,SIZE(h_block,2)
                     force_ab(icor)=force_ab(icor)+(ha(i)+hb(j))*dsmat(i,j)
                   END DO
                 END DO
               ELSE
                 DO i=1,SIZE(h_block,1)
                   DO j=1,SIZE(h_block,2)
                     force_ab(icor)=force_ab(icor)+(ha(j)+hb(i))*dsmat(i,j)
                   END DO
                 END DO
               END IF
            END DO
          END IF

        END IF

        IF(calculate_forces .AND. (iatom/=jatom .OR. dr > rij_threshold)) THEN
          IF ( irow == iatom ) force_ab = -force_ab
          force(ikind)%all_potential(:,atom_a) =&
              force(ikind)%all_potential(:,atom_a) - force_ab(:)
          force(jkind)%all_potential(:,atom_b) =&
              force(jkind)%all_potential(:,atom_b) + force_ab(:)
          IF (use_virial) THEN
            CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
          END IF
        END IF

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     DEALLOCATE (se_defined,hmt,umt,ZSt,ZPt,nrt,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL cp_dbcsr_sum_replicated(diagmat_h, error=error)
     CALL cp_dbcsr_distribute(diagmat_h,error=error)
     CALL cp_dbcsr_add(matrix_h(1)%matrix, diagmat_h,1.0_dp,1.0_dp,error=error)
     CALL set_qs_env(qs_env=qs_env,matrix_h=matrix_h,error=error)

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log",error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,para_env,&
             scale=evolt,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", error=error)
     END IF

     IF(calculate_forces) THEN
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,alpha_scalar=1.0_dp,beta_scalar=-1.0_dp,error=error)
       END IF
       DEALLOCATE(atom_of_kind,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_dbcsr_deallocate_matrix(diagmat_p,error=error)
     END IF

     CALL cp_dbcsr_deallocate_matrix(diagmat_h,error=error)

     CALL timestop(handle)

  END SUBROUTINE build_se_core_matrix

! *****************************************************************************
!> \brief ...
!> \param R ...
!> \param nra ...
!> \param nrb ...
!> \param ZSA ...
!> \param ZSB ...
!> \param ZPA ...
!> \param ZPB ...
!> \param S ...
! *****************************************************************************
SUBROUTINE makeS(R,nra,nrb,ZSA,ZSB,ZPA,ZPB,S)

    REAL(kind=dp), DIMENSION(3)              :: R
    INTEGER                                  :: nra, nrb
    REAL(kind=dp)                            :: ZSA, ZSB, ZPA, ZPB
    REAL(kind=dp), DIMENSION(4, 4)           :: S

    INTEGER, DIMENSION(4, 4), PARAMETER :: nc1 = RESHAPE( (/  2,  4,  4,  6,  &
      4,  3,  6,  7,  4,  6,  4,  8,  6,  7,  8,  5 /) , (/ 4, 4/) ), nc2 = &
      RESHAPE( (/  4,  4,  8,  8,  6,  8,  6, 12,  8,  8, 12,  8, 10, 12, 14, &
      16 /) , (/ 4, 4/) ), nc3 = RESHAPE( (/  4,  6,  8, 10,  4,  8,  8, 12,  &
      8,  6, 12, 14,  8, 12,  8, 16 /) , (/ 4, 4/) ), nc4 = RESHAPE( (/  4,  8&
      , 11, 14,  8,  6, 12, 14, 11, 12, 10, 20, 14, 14, 20, 12 /) , (/ 4, 4/) &
      ), nc5 = RESHAPE( (/  2,  4,  6,  8,  4,  4,  8,  8,  6,  8,  6, 12,  8,&
      8, 12,  8 /) , (/ 4, 4/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c1 = RESHAPE( (/  1,  1,  1,  1&
      ,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  2,  3, -1, -2&
      ,  1,  2, -2, -1, -3,  1, -3, -2, -1, -4,  0, -1, -2,  2, -1,  1, -2, -1&
      ,  2, -2,  3, -3,  2, -1, -3,  6,  0, -1, -1, -2,  1,  0, -2, -4, -1,  2&
      , -1, -3,  2,  4,  3, -4,  0,  0,  0, -3,  0,  0,  1, -1,  0,  1,  0,  3&
      , -3, -1,  3,  1,  0,  0,  0, -1,  0,  0,  1,  2,  0, -1,  0,  3,  1, -2&
      , -3,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0,  1, -1,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c2 = RESHAPE( (/  1,  1,  1,  1&
      ,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  1,  2, -2, -1&
      ,  1,  1, -3, -2, -1,  1, -4, -3, -2, -1,  1, -1,  1,  1,  1,  1, -2,  1&
      ,  1,  1,  1, -3,  1,  1,  1,  1, -1, -1, -1,  2,  1, -1, -2, -2,  3, -2&
      , -2, -3,  6,  2, -1, -3,  0,  0,  1, -2, -2, -1,  1,  1, -3,  2, -1,  3&
      , -4, -3, -2, -1,  0,  0, -1, -1,  1,  1,  1, -2, -1, -1,  2,  3, -4,  2&
      ,  4,  3,  0,  0, -1, -2,  0, -1,  0, -2,  3,  2, -2, -1,  6,  2, -1, -3&
      ,  0,  0, -1, -1,  0,  1,  0,  1, -1, -1,  1, -1,  1, -3, -1,  3,  0,  0&
      ,  0,  0,  0,  0,  0, -2,  0,  0,  2,  0, -4,  2,  4,  3,  0,  0,  0,  0&
      ,  0,  0,  0,  1,  0,  0, -1,  0,  1,  1, -2, -3,  0,  0,  0,  0,  0,  0&
      ,  0,  1,  0,  0,  1,  0,  0, -3, -1,  3,  0,  0,  0,  0,  0,  0,  0,  1&
      ,  0,  0, -1,  0,  0,  1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0, -2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c3 = RESHAPE( (/ -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -3, -4,  1, -1&
      , -2, -3,  1,  1, -1, -2,  2,  1,  1, -1,  1,  1,  1,  1,  1,  1,  1,  1&
      ,  1,  2,  1,  1,  1,  1,  3,  1,  1, -1, -3, -6, -1,  1,  2, -2,  1, -2&
      ,  2,  1, -2,  2, -3,  3,  0,  2,  3,  4,  0,  1,  2,  3, -1, -1,  1,  2&
      , -2, -1, -3,  1,  0,  1, -1, -4,  0,  1,  1,  2, -1,  1,  2,  4,  1, -2&
      ,  3,  3,  0,  0,  3,  6,  0, -1, -2,  2, -1,  0, -2, -1,  2, -2,  1, -3&
      ,  0,  0,  1, -1,  0, -1, -1,  3,  1,  0, -1,  1, -1, -1, -1, -3,  0,  0&
      ,  0,  4,  0,  0,  0, -2,  0,  0, -2, -4,  0,  2,  0, -3,  0,  0,  0,  1&
      ,  0,  0,  0,  1,  0,  0, -1, -2,  0,  1,  0, -3,  0,  0,  0,  0,  0,  0&
      ,  0, -3,  0,  0,  1, -1,  0,  1,  0,  3,  0,  0,  0,  0,  0,  0,  0, -1&
      ,  0,  0,  1, -1,  0, -1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  2,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1&
      ,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c4 = RESHAPE( (/ -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1&
      ,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -2, -3,  1,  1, -1, -2&
      ,  2,  1,  2, -1,  3,  2,  1,  3, -1,  1,  2,  3, -1, -1,  1,  2, -2, -1&
      , -1,  1, -3, -2, -1, -2,  0,  1, -1, -3,  1, -1,  1,  1, -1,  1, -1,  2&
      , -3,  1,  2, -1,  0, -1,  2,  4, -1,  1, -1, -1,  2, -1, -1, -1,  4, -1&
      , -1, -3,  0,  1, -1, -1, -1,  0,  1,  2, -1, -1, -1, -1, -1, -2, -1,  3&
      ,  0, -1,  2, -1,  1,  0, -1, -2, -2,  1,  2,  2,  1,  2, -2,  1,  0,  0&
      , -2,  4,  0,  0, -1,  1,  2, -1,  1, -1, -4,  1,  1,  2,  0,  0,  1, -3&
      ,  0,  0,  1, -1,  1,  1, -1, -1,  3, -1,  1, -3,  0,  0, -1,  3,  0,  0&
      , -1, -2, -1,  1,  0, -1,  3,  2, -1, -1,  0,  0,  0, -3,  0,  0,  1,  2&
      ,  0, -1,  0, -1, -3, -2, -1,  1,  0,  0,  0,  1,  0,  0,  0, -1,  0,  0&
      ,  0,  2, -1, -1,  2,  0,  0,  0,  0, -1,  0,  0,  0,  1,  0,  0,  0, -1&
      ,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0&
      ,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c5 = RESHAPE( (/ -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1, -1, -2, -3,  1,  1&
      , -1, -2,  2,  1,  2, -1,  3,  2,  1,  3,  0,  1, -1, -3,  1,  1,  1,  1&
      , -1,  1,  1,  2, -3,  1,  2,  1,  0,  1,  1,  1, -1, -1,  1,  2,  1,  1&
      , -1,  1,  1, -2,  1, -3,  0,  0,  2, -1,  0,  0,  1,  2, -2, -1, -2,  2&
      ,  1, -2, -2, -3,  0,  0,  1,  3,  0,  0,  1,  1,  1, -1,  1,  1, -3,  1&
      , -1,  1,  0,  0,  0,  3,  0,  0, -1, -2,  0, -1,  0, -1,  3,  2, -1,  3&
      ,  0,  0,  0,  1,  0,  0, -1, -1,  0,  1,  0, -2, -1, -1, -2, -1,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma1 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  0,  2,  3,  4,  2,  &
      2,  4,  5,  3,  4,  4,  6,  4,  5,  6,  6,  0,  1,  1,  3,  1,  0,  3,  &
      4,  1,  3,  2,  5,  3,  4,  5,  4,  0,  0,  0,  2,  0,  0,  2,  3,  0,  &
      2,  0,  4,  2,  3,  4,  2,  0,  0,  0,  1,  0,  0,  1,  2,  0,  1,  0,  &
      3,  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,  0,  &
      1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma2 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  1,  2,  3,  4,  2,  &
      3,  4,  5,  3,  4,  5,  6,  4,  5,  6,  7,  1,  1,  3,  4,  2,  3,  3,  &
      5,  3,  4,  5,  5,  4,  5,  6,  7,  0,  0,  2,  3,  1,  2,  2,  4,  2,  &
      3,  4,  4,  3,  4,  5,  6,  0,  0,  2,  2,  1,  2,  1,  4,  2,  2,  4,  &
      3,  3,  4,  5,  6,  0,  0,  1,  1,  0,  1,  0,  3,  1,  1,  3,  2,  2,  &
      3,  4,  5,  0,  0,  1,  1,  0,  1,  0,  3,  1,  1,  3,  1,  2,  3,  4,  &
      5,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  1,  2,  3,  4,  0,  &
      0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  1,  2,  3,  4,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  1,  2,  3,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  1,  0,  0,  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma3 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  1,  2,  3,  4,  2,  &
      3,  4,  5,  3,  4,  5,  6,  4,  5,  6,  7,  1,  2,  3,  4,  1,  3,  4,  &
      5,  3,  3,  5,  6,  4,  5,  5,  7,  0,  1,  2,  3,  0,  2,  3,  4,  2,  &
      2,  4,  5,  3,  4,  4,  6,  0,  1,  2,  3,  0,  2,  2,  4,  2,  1,  4,  &
      5,  2,  4,  3,  6,  0,  0,  1,  2,  0,  1,  1,  3,  1,  0,  3,  4,  1,  &
      3,  2,  5,  0,  0,  1,  2,  0,  1,  1,  3,  1,  0,  3,  4,  1,  3,  1,  &
      5,  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  2,  3,  0,  2,  0,  4,  0,  &
      0,  0,  1,  0,  0,  0,  2,  0,  0,  2,  3,  0,  2,  0,  4,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  1,  2,  0,  1,  0,  3,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  1,  2,  0,  1,  0,  3,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma4 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  2,  3,  4,  5,  3,  &
      4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  0,  2,  3,  4,  2,  2,  4,  &
      5,  3,  4,  4,  6,  4,  5,  6,  6,  0,  2,  3,  4,  2,  2,  4,  5,  3,  &
      4,  4,  6,  4,  5,  6,  6,  0,  1,  2,  3,  1,  0,  3,  4,  2,  3,  4,  &
      5,  3,  4,  5,  6,  0,  1,  2,  3,  1,  0,  3,  4,  2,  3,  2,  5,  3,  &
      4,  5,  4,  0,  0,  2,  3,  0,  0,  2,  3,  2,  2,  2,  5,  3,  3,  5,  &
      4,  0,  0,  1,  2,  0,  0,  2,  3,  1,  2,  2,  4,  2,  3,  4,  2,  0,  &
      0,  1,  2,  0,  0,  1,  2,  1,  1,  0,  4,  2,  2,  4,  2,  0,  0,  0,  &
      2,  0,  0,  1,  2,  0,  1,  0,  4,  2,  2,  4,  2,  0,  0,  0,  1,  0,  &
      0,  0,  1,  0,  0,  0,  3,  1,  1,  3,  0,  0,  0,  0,  1,  0,  0,  0,  &
      1,  0,  0,  0,  3,  1,  1,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  3,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      2,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  &
      0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma5 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  0,  2,  3,  4,  2,  &
      2,  4,  5,  3,  4,  4,  6,  4,  5,  6,  6,  0,  1,  2,  3,  1,  2,  3,  &
      4,  2,  3,  4,  5,  3,  4,  5,  6,  0,  0,  2,  3,  0,  0,  3,  3,  2,  &
      3,  2,  5,  3,  3,  5,  4,  0,  0,  1,  2,  0,  0,  2,  3,  1,  2,  2,  &
      4,  2,  3,  4,  4,  0,  0,  0,  2,  0,  0,  2,  2,  0,  2,  0,  4,  2,  &
      2,  4,  2,  0,  0,  0,  1,  0,  0,  1,  1,  0,  1,  0,  3,  1,  1,  3,  &
      2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb1 = RESHAPE( (/  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1,  1,  1,  1,  &
      2,  1,  1,  1,  1,  2,  1,  1,  1,  1,  2,  0,  2,  3,  2,  2,  4,  2,  &
      2,  3,  2,  4,  2,  2,  2,  2,  4,  0,  3,  4,  3,  3,  0,  3,  3,  4,  &
      3,  6,  3,  3,  3,  3,  6,  0,  0,  0,  4,  0,  0,  4,  4,  0,  4,  0,  &
      4,  4,  4,  4,  8,  0,  0,  0,  5,  0,  0,  5,  5,  0,  5,  0,  5,  5,  &
      5,  5,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  6,  0,  6,  6,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  7,  0,  0,  7,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb2 = RESHAPE( (/  1,  1,  1,  &
      1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  0,  2,  2,  2,  &
      2,  0,  2,  2,  2,  2,  0,  2,  2,  2,  2,  0,  3,  0,  0,  0,  0,  3,  &
      0,  0,  0,  0,  3,  0,  0,  0,  0,  1,  2,  3,  1,  3,  3,  2,  3,  3,  &
      1,  3,  2,  3,  3,  3,  3,  0,  0,  1,  4,  1,  1,  5,  1,  1,  4,  1,  &
      5,  1,  1,  1,  1,  0,  0,  4,  5,  2,  4,  4,  4,  4,  5,  4,  4,  4,  &
      4,  4,  4,  0,  0,  2,  3,  0,  2,  0,  2,  2,  3,  2,  7,  2,  2,  2,  &
      2,  0,  0,  3,  4,  0,  3,  0,  5,  3,  4,  5,  6,  5,  5,  5,  5,  0,  &
      0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  0,  3,  3,  3,  3,  0,  0,  0,  &
      0,  0,  0,  0,  6,  0,  0,  6,  0,  4,  6,  6,  6,  0,  0,  0,  0,  0,  &
      0,  0,  4,  0,  0,  4,  0,  0,  4,  4,  4,  0,  0,  0,  0,  0,  0,  0,  &
      5,  0,  0,  5,  0,  0,  5,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  5,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  6,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb3 = RESHAPE( (/  1,  1,  1,  &
      1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  0,  &
      2,  2,  2,  2,  0,  2,  2,  2,  2,  0,  2,  0,  0,  0,  0,  3,  0,  0,  &
      0,  0,  3,  0,  0,  0,  0,  3,  0,  1,  3,  3,  3,  2,  3,  1,  3,  3,  &
      2,  3,  3,  1,  3,  2,  3,  0,  1,  1,  1,  0,  1,  4,  1,  1,  5,  1,  &
      1,  4,  1,  5,  1,  0,  2,  4,  4,  0,  4,  5,  4,  4,  4,  4,  4,  5,  &
      4,  4,  4,  0,  0,  2,  2,  0,  2,  3,  2,  2,  0,  2,  2,  3,  2,  7,  &
      2,  0,  0,  3,  5,  0,  3,  4,  5,  3,  0,  5,  5,  4,  5,  6,  5,  0,  &
      0,  0,  3,  0,  0,  0,  3,  0,  0,  3,  3,  0,  3,  0,  3,  0,  0,  0,  &
      4,  0,  0,  0,  6,  0,  0,  6,  6,  0,  6,  0,  6,  0,  0,  0,  0,  0,  &
      0,  0,  4,  0,  0,  4,  4,  0,  4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  &
      5,  0,  0,  5,  7,  0,  5,  0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  5,  0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      6,  0,  0,  0,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb4 = RESHAPE( (/  2,  2,  2,  &
      2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  3,  3,  3,  3,  4,  3,  &
      3,  3,  3,  4,  3,  3,  3,  3,  4,  0,  1,  1,  1,  1,  0,  1,  1,  1,  &
      1,  2,  1,  1,  1,  1,  2,  0,  2,  4,  4,  2,  4,  4,  2,  4,  4,  0,  &
      4,  4,  2,  4,  0,  0,  0,  2,  2,  0,  2,  0,  0,  2,  0,  6,  2,  2,  &
      0,  2,  6,  0,  3,  0,  0,  3,  0,  5,  5,  0,  5,  4,  0,  0,  5,  0,  &
      2,  0,  1,  3,  5,  1,  0,  1,  1,  3,  1,  2,  5,  5,  1,  5,  8,  0,  &
      0,  1,  3,  0,  0,  4,  6,  1,  4,  6,  3,  3,  6,  3,  6,  0,  0,  4,  &
      1,  0,  0,  2,  4,  4,  2,  4,  1,  1,  4,  1,  4,  0,  0,  2,  4,  0,  &
      0,  5,  5,  2,  5,  0,  6,  4,  5,  6,  8,  0,  0,  0,  2,  0,  0,  3,  &
      3,  0,  3,  0,  4,  2,  3,  4,  6,  0,  0,  0,  5,  0,  0,  0,  6,  0,  &
      0,  0,  2,  5,  6,  2,  0,  0,  0,  0,  3,  0,  0,  0,  4,  0,  0,  0,  &
      7,  3,  4,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  &
      0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  6,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  4,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  7,  0,  0,  7,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  5,  0,  0,  5,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb5 = RESHAPE( (/  2,  2,  2,  &
      2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  0,  3,  3,  3,  3,  4&
      ,  3,  3,  3,  3,  4,  3,  3,  3,  3,  4,  0,  0,  4,  4,  0,  0,  4,  0&
      ,  4,  4,  0,  4,  4,  0,  4,  0,  0,  1,  0,  0,  1,  2,  0,  5,  0,  0&
      ,  6,  0,  0,  5,  0,  6,  0,  0,  1,  5,  0,  0,  5,  1,  1,  5,  2,  5&
      ,  5,  1,  5,  2,  0,  0,  2,  1,  0,  0,  1,  6,  2,  1,  4,  1,  1,  6&
      ,  1,  8,  0,  0,  0,  2,  0,  0,  2,  3,  0,  2,  0,  6,  2,  3,  6,  4&
      ,  0,  0,  0,  3,  0,  0,  3,  4,  0,  3,  0,  2,  3,  4,  2,  6,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  7,  0,  0,  7,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  4,  0,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  5,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )

    INTEGER                                  :: k, k1, k2, mu
    REAL(kind=dp)                            :: cp, ct, fac1, fac2, J, Jc, &
                                                Jcc, Jss, rr, sp, st, xx, yy, &
                                                za, zb
    REAL(kind=dp), DIMENSION(3)              :: v
    REAL(kind=dp), DIMENSION(3, 3)           :: Arot

S(:,:)=  0.0_dp

      v(:)=R(:)
      rr=SQRT(DOT_PRODUCT(v,v))

IF( rr < 1.0e-20_dp )THEN

      DO mu=1,4
        S(mu,mu)=1.0_dp
      ENDDO

ELSE

      fac1=1.0_dp
      IF(nra == 1)THEN
        fac1=fac1*2.0_dp
      ELSE
        IF(nra == 2)THEN
          fac1=fac1*SQRT(4.0_dp/3.0_dp)
        ELSE
          IF(nra == 3)THEN
            fac1=fac1*SQRT(8.0_dp/45.0_dp)
          ELSE
            IF(nra == 4)THEN
              fac1=fac1*SQRT(4.0_dp/315.0_dp)
            ELSE
              WRITE(*,*) 'nra= ',nra
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(nrb == 1)THEN
        fac1=fac1*2.0_dp
      ELSE
        IF(nrb == 2)THEN
          fac1=fac1*SQRT(4.0_dp/3.0_dp)
        ELSE
          IF(nrb == 3)THEN
            fac1=fac1*SQRT(8.0_dp/45.0_dp)
          ELSE
            IF(nrb == 4)THEN
              fac1=fac1*SQRT(4.0_dp/315.0_dp)
            ELSE
              WRITE(*,*) 'nrb= ',nrb
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      ct= -v(3)/rr
      IF(ABS(ct) < 1.0_dp)THEN
        st= SQRT(1.0_dp-ct**2)
        cp= -v(1)/(rr*st)
        sp= -v(2)/(rr*st)
        Arot(1,1)=  ct*cp
        Arot(1,2)=    -sp
        Arot(1,3)=  st*cp
        Arot(2,1)=  ct*sp
        Arot(2,2)=     cp
        Arot(2,3)=  st*sp
        Arot(3,1)= -st
        Arot(3,2)=  0.0_dp
        Arot(3,3)=  ct
      ELSE
        Arot(1,1)=  ct
        Arot(1,2)=  0.0_dp
        Arot(1,3)=  0.0_dp
        Arot(2,1)=  0.0_dp
        Arot(2,2)=  1.0_dp
        Arot(2,3)=  0.0_dp
        Arot(3,1)=  0.0_dp
        Arot(3,2)=  0.0_dp
        Arot(3,3)=  ct
      ENDIF

      za=ZSA
      zb=ZSB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)

      J=0.0_dp
      DO k=1,nc1(nra,nrb)
        J=J+REAL(c1(nra,nrb,k),dp)*AA(ma1(nra,nrb,k),xx)*BB(mb1(nra,nrb,k),yy)
      ENDDO
      J=J*rr**(nra+nrb+1)
      J=J/2.0_dp**(nra+nrb+2)

      S(1,1)=S(1,1)+fac1*fac2*J

      za=ZPA
      zb=ZSB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)

      Jc=0.0_dp
      DO k=1,nc2(nra,nrb)
        Jc=Jc+REAL(c2(nra,nrb,k),dp)*AA(ma2(nra,nrb,k),xx)*BB(mb2(nra,nrb,k),yy)
      ENDDO
      Jc=Jc*rr**(nra+nrb+1)
      Jc=Jc/2.0_dp**(nra+nrb+2)

      DO k1=1,3
        S(k1+1,1)=S(k1+1,1) &
     &          +SQRT(3.0_dp)*Arot(k1,3)*fac1*fac2*Jc
      ENDDO

      za=ZSA
      zb=ZPB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)

      Jc=0.0_dp
      DO k=1,nc3(nra,nrb)
        Jc=Jc+REAL(c3(nra,nrb,k),dp)*AA(ma3(nra,nrb,k),xx)*BB(mb3(nra,nrb,k),yy)
      ENDDO
      Jc=Jc*rr**(nra+nrb+1)
      Jc=Jc/2.0_dp**(nra+nrb+2)

      DO k1=1,3
        S(1,k1+1)=S(1,k1+1) &
     &          -SQRT(3.0_dp)*Arot(k1,3)*fac1*fac2*Jc
      ENDDO

      za=ZPA
      zb=ZPB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)

      Jss=0.0_dp
      DO k=1,nc4(nra,nrb)
        Jss=Jss+REAL(c4(nra,nrb,k),dp)*AA(ma4(nra,nrb,k),xx)*BB(mb4(nra,nrb,k),yy)
      ENDDO
      Jss=Jss*rr**(nra+nrb+1)
      Jss=Jss/2.0_dp**(nra+nrb+2)

      Jcc=0.0_dp
      DO k=1,nc5(nra,nrb)
        Jcc=Jcc+REAL(c5(nra,nrb,k),dp)*AA(ma5(nra,nrb,k),xx)*BB(mb5(nra,nrb,k),yy)
      ENDDO
      Jcc=Jcc*rr**(nra+nrb+1)
      Jcc=Jcc/2.0_dp**(nra+nrb+2)

      DO k1=1,3
        DO k2=1,3
          S(k1+1,k2+1)=S(k1+1,k2+1) &
     &            +1.5_dp*Arot(k1,1)*Arot(k2,1)*fac1*fac2*Jss &
     &            +1.5_dp*Arot(k1,2)*Arot(k2,2)*fac1*fac2*Jss &
     &            -3.0_dp*Arot(k1,3)*Arot(k2,3)*fac1*fac2*Jcc
        ENDDO
      ENDDO

ENDIF

END SUBROUTINE makeS

! *****************************************************************************
!> \brief ...
!> \param R ...
!> \param nra ...
!> \param nrb ...
!> \param ZSA ...
!> \param ZSB ...
!> \param ZPA ...
!> \param ZPB ...
!> \param dS ...
! *****************************************************************************
SUBROUTINE makedS(R,nra,nrb,ZSA,ZSB,ZPA,ZPB,dS)

    REAL(kind=dp), DIMENSION(3)              :: R
    INTEGER                                  :: nra, nrb
    REAL(kind=dp)                            :: ZSA, ZSB, ZPA, ZPB
    REAL(kind=dp), DIMENSION(4, 4, 3)        :: dS

    INTEGER, DIMENSION(4, 4), PARAMETER :: nc1 = RESHAPE( (/  2,  4,  4,  6,  &
      4,  3,  6,  7,  4,  6,  4,  8,  6,  7,  8,  5 /) , (/ 4, 4/) ), nc2 = &
      RESHAPE( (/  4,  4,  8,  8,  6,  8,  6, 12,  8,  8, 12,  8, 10, 12, 14, &
      16 /) , (/ 4, 4/) ), nc3 = RESHAPE( (/  4,  6,  8, 10,  4,  8,  8, 12,  &
      8,  6, 12, 14,  8, 12,  8, 16 /) , (/ 4, 4/) ), nc4 = RESHAPE( (/  4,  8&
      , 11, 14,  8,  6, 12, 14, 11, 12, 10, 20, 14, 14, 20, 12 /) , (/ 4, 4/) &
      ), nc5 = RESHAPE( (/  2,  4,  6,  8,  4,  4,  8,  8,  6,  8,  6, 12,  8,&
      8, 12,  8 /) , (/ 4, 4/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c1 = RESHAPE( (/  1,  1,  1,  1&
      ,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  2,  3, -1, -2&
      ,  1,  2, -2, -1, -3,  1, -3, -2, -1, -4,  0, -1, -2,  2, -1,  1, -2, -1&
      ,  2, -2,  3, -3,  2, -1, -3,  6,  0, -1, -1, -2,  1,  0, -2, -4, -1,  2&
      , -1, -3,  2,  4,  3, -4,  0,  0,  0, -3,  0,  0,  1, -1,  0,  1,  0,  3&
      , -3, -1,  3,  1,  0,  0,  0, -1,  0,  0,  1,  2,  0, -1,  0,  3,  1, -2&
      , -3,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0,  1, -1,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c2 = RESHAPE( (/  1,  1,  1,  1&
      ,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1,  1,  1,  2, -2, -1&
      ,  1,  1, -3, -2, -1,  1, -4, -3, -2, -1,  1, -1,  1,  1,  1,  1, -2,  1&
      ,  1,  1,  1, -3,  1,  1,  1,  1, -1, -1, -1,  2,  1, -1, -2, -2,  3, -2&
      , -2, -3,  6,  2, -1, -3,  0,  0,  1, -2, -2, -1,  1,  1, -3,  2, -1,  3&
      , -4, -3, -2, -1,  0,  0, -1, -1,  1,  1,  1, -2, -1, -1,  2,  3, -4,  2&
      ,  4,  3,  0,  0, -1, -2,  0, -1,  0, -2,  3,  2, -2, -1,  6,  2, -1, -3&
      ,  0,  0, -1, -1,  0,  1,  0,  1, -1, -1,  1, -1,  1, -3, -1,  3,  0,  0&
      ,  0,  0,  0,  0,  0, -2,  0,  0,  2,  0, -4,  2,  4,  3,  0,  0,  0,  0&
      ,  0,  0,  0,  1,  0,  0, -1,  0,  1,  1, -2, -3,  0,  0,  0,  0,  0,  0&
      ,  0,  1,  0,  0,  1,  0,  0, -3, -1,  3,  0,  0,  0,  0,  0,  0,  0,  1&
      ,  0,  0, -1,  0,  0,  1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0, -2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c3 = RESHAPE( (/ -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -3, -4,  1, -1&
      , -2, -3,  1,  1, -1, -2,  2,  1,  1, -1,  1,  1,  1,  1,  1,  1,  1,  1&
      ,  1,  2,  1,  1,  1,  1,  3,  1,  1, -1, -3, -6, -1,  1,  2, -2,  1, -2&
      ,  2,  1, -2,  2, -3,  3,  0,  2,  3,  4,  0,  1,  2,  3, -1, -1,  1,  2&
      , -2, -1, -3,  1,  0,  1, -1, -4,  0,  1,  1,  2, -1,  1,  2,  4,  1, -2&
      ,  3,  3,  0,  0,  3,  6,  0, -1, -2,  2, -1,  0, -2, -1,  2, -2,  1, -3&
      ,  0,  0,  1, -1,  0, -1, -1,  3,  1,  0, -1,  1, -1, -1, -1, -3,  0,  0&
      ,  0,  4,  0,  0,  0, -2,  0,  0, -2, -4,  0,  2,  0, -3,  0,  0,  0,  1&
      ,  0,  0,  0,  1,  0,  0, -1, -2,  0,  1,  0, -3,  0,  0,  0,  0,  0,  0&
      ,  0, -3,  0,  0,  1, -1,  0,  1,  0,  3,  0,  0,  0,  0,  0,  0,  0, -1&
      ,  0,  0,  1, -1,  0, -1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  2,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1&
      ,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c4 = RESHAPE( (/ -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1&
      ,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -2, -3,  1,  1, -1, -2&
      ,  2,  1,  2, -1,  3,  2,  1,  3, -1,  1,  2,  3, -1, -1,  1,  2, -2, -1&
      , -1,  1, -3, -2, -1, -2,  0,  1, -1, -3,  1, -1,  1,  1, -1,  1, -1,  2&
      , -3,  1,  2, -1,  0, -1,  2,  4, -1,  1, -1, -1,  2, -1, -1, -1,  4, -1&
      , -1, -3,  0,  1, -1, -1, -1,  0,  1,  2, -1, -1, -1, -1, -1, -2, -1,  3&
      ,  0, -1,  2, -1,  1,  0, -1, -2, -2,  1,  2,  2,  1,  2, -2,  1,  0,  0&
      , -2,  4,  0,  0, -1,  1,  2, -1,  1, -1, -4,  1,  1,  2,  0,  0,  1, -3&
      ,  0,  0,  1, -1,  1,  1, -1, -1,  3, -1,  1, -3,  0,  0, -1,  3,  0,  0&
      , -1, -2, -1,  1,  0, -1,  3,  2, -1, -1,  0,  0,  0, -3,  0,  0,  1,  2&
      ,  0, -1,  0, -1, -3, -2, -1,  1,  0,  0,  0,  1,  0,  0,  0, -1,  0,  0&
      ,  0,  2, -1, -1,  2,  0,  0,  0,  0, -1,  0,  0,  0,  1,  0,  0,  0, -1&
      ,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0&
      ,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: c5 = RESHAPE( (/ -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1, -1, -2, -3,  1,  1&
      , -1, -2,  2,  1,  2, -1,  3,  2,  1,  3,  0,  1, -1, -3,  1,  1,  1,  1&
      , -1,  1,  1,  2, -3,  1,  2,  1,  0,  1,  1,  1, -1, -1,  1,  2,  1,  1&
      , -1,  1,  1, -2,  1, -3,  0,  0,  2, -1,  0,  0,  1,  2, -2, -1, -2,  2&
      ,  1, -2, -2, -3,  0,  0,  1,  3,  0,  0,  1,  1,  1, -1,  1,  1, -3,  1&
      , -1,  1,  0,  0,  0,  3,  0,  0, -1, -2,  0, -1,  0, -1,  3,  2, -1,  3&
      ,  0,  0,  0,  1,  0,  0, -1, -1,  0,  1,  0, -2, -1, -1, -2, -1,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma1 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  0,  2,  3,  4,  2,  &
      2,  4,  5,  3,  4,  4,  6,  4,  5,  6,  6,  0,  1,  1,  3,  1,  0,  3,  &
      4,  1,  3,  2,  5,  3,  4,  5,  4,  0,  0,  0,  2,  0,  0,  2,  3,  0,  &
      2,  0,  4,  2,  3,  4,  2,  0,  0,  0,  1,  0,  0,  1,  2,  0,  1,  0,  &
      3,  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,  0,  &
      1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma2 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  1,  2,  3,  4,  2,  &
      3,  4,  5,  3,  4,  5,  6,  4,  5,  6,  7,  1,  1,  3,  4,  2,  3,  3,  &
      5,  3,  4,  5,  5,  4,  5,  6,  7,  0,  0,  2,  3,  1,  2,  2,  4,  2,  &
      3,  4,  4,  3,  4,  5,  6,  0,  0,  2,  2,  1,  2,  1,  4,  2,  2,  4,  &
      3,  3,  4,  5,  6,  0,  0,  1,  1,  0,  1,  0,  3,  1,  1,  3,  2,  2,  &
      3,  4,  5,  0,  0,  1,  1,  0,  1,  0,  3,  1,  1,  3,  1,  2,  3,  4,  &
      5,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  1,  2,  3,  4,  0,  &
      0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  1,  2,  3,  4,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  1,  2,  3,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  1,  0,  0,  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma3 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  1,  2,  3,  4,  2,  &
      3,  4,  5,  3,  4,  5,  6,  4,  5,  6,  7,  1,  2,  3,  4,  1,  3,  4,  &
      5,  3,  3,  5,  6,  4,  5,  5,  7,  0,  1,  2,  3,  0,  2,  3,  4,  2,  &
      2,  4,  5,  3,  4,  4,  6,  0,  1,  2,  3,  0,  2,  2,  4,  2,  1,  4,  &
      5,  2,  4,  3,  6,  0,  0,  1,  2,  0,  1,  1,  3,  1,  0,  3,  4,  1,  &
      3,  2,  5,  0,  0,  1,  2,  0,  1,  1,  3,  1,  0,  3,  4,  1,  3,  1,  &
      5,  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  2,  3,  0,  2,  0,  4,  0,  &
      0,  0,  1,  0,  0,  0,  2,  0,  0,  2,  3,  0,  2,  0,  4,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  1,  2,  0,  1,  0,  3,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  1,  2,  0,  1,  0,  3,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma4 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  2,  3,  4,  5,  3,  &
      4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  0,  2,  3,  4,  2,  2,  4,  &
      5,  3,  4,  4,  6,  4,  5,  6,  6,  0,  2,  3,  4,  2,  2,  4,  5,  3,  &
      4,  4,  6,  4,  5,  6,  6,  0,  1,  2,  3,  1,  0,  3,  4,  2,  3,  4,  &
      5,  3,  4,  5,  6,  0,  1,  2,  3,  1,  0,  3,  4,  2,  3,  2,  5,  3,  &
      4,  5,  4,  0,  0,  2,  3,  0,  0,  2,  3,  2,  2,  2,  5,  3,  3,  5,  &
      4,  0,  0,  1,  2,  0,  0,  2,  3,  1,  2,  2,  4,  2,  3,  4,  2,  0,  &
      0,  1,  2,  0,  0,  1,  2,  1,  1,  0,  4,  2,  2,  4,  2,  0,  0,  0,  &
      2,  0,  0,  1,  2,  0,  1,  0,  4,  2,  2,  4,  2,  0,  0,  0,  1,  0,  &
      0,  0,  1,  0,  0,  0,  3,  1,  1,  3,  0,  0,  0,  0,  1,  0,  0,  0,  &
      1,  0,  0,  0,  3,  1,  1,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  3,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      2,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  &
      0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: ma5 = RESHAPE( (/  2,  3,  4,  &
      5,  3,  4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,  0,  2,  3,  4,  2,  &
      2,  4,  5,  3,  4,  4,  6,  4,  5,  6,  6,  0,  1,  2,  3,  1,  2,  3,  &
      4,  2,  3,  4,  5,  3,  4,  5,  6,  0,  0,  2,  3,  0,  0,  3,  3,  2,  &
      3,  2,  5,  3,  3,  5,  4,  0,  0,  1,  2,  0,  0,  2,  3,  1,  2,  2,  &
      4,  2,  3,  4,  4,  0,  0,  0,  2,  0,  0,  2,  2,  0,  2,  0,  4,  2,  &
      2,  4,  2,  0,  0,  0,  1,  0,  0,  1,  1,  0,  1,  0,  3,  1,  1,  3,  &
      2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb1 = RESHAPE( (/  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1,  1,  1,  1,  &
      2,  1,  1,  1,  1,  2,  1,  1,  1,  1,  2,  0,  2,  3,  2,  2,  4,  2,  &
      2,  3,  2,  4,  2,  2,  2,  2,  4,  0,  3,  4,  3,  3,  0,  3,  3,  4,  &
      3,  6,  3,  3,  3,  3,  6,  0,  0,  0,  4,  0,  0,  4,  4,  0,  4,  0,  &
      4,  4,  4,  4,  8,  0,  0,  0,  5,  0,  0,  5,  5,  0,  5,  0,  5,  5,  &
      5,  5,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  6,  0,  6,  6,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  7,  0,  0,  7,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb2 = RESHAPE( (/  1,  1,  1,  &
      1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  0,  2,  2,  2,  &
      2,  0,  2,  2,  2,  2,  0,  2,  2,  2,  2,  0,  3,  0,  0,  0,  0,  3,  &
      0,  0,  0,  0,  3,  0,  0,  0,  0,  1,  2,  3,  1,  3,  3,  2,  3,  3,  &
      1,  3,  2,  3,  3,  3,  3,  0,  0,  1,  4,  1,  1,  5,  1,  1,  4,  1,  &
      5,  1,  1,  1,  1,  0,  0,  4,  5,  2,  4,  4,  4,  4,  5,  4,  4,  4,  &
      4,  4,  4,  0,  0,  2,  3,  0,  2,  0,  2,  2,  3,  2,  7,  2,  2,  2,  &
      2,  0,  0,  3,  4,  0,  3,  0,  5,  3,  4,  5,  6,  5,  5,  5,  5,  0,  &
      0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  0,  3,  3,  3,  3,  0,  0,  0,  &
      0,  0,  0,  0,  6,  0,  0,  6,  0,  4,  6,  6,  6,  0,  0,  0,  0,  0,  &
      0,  0,  4,  0,  0,  4,  0,  0,  4,  4,  4,  0,  0,  0,  0,  0,  0,  0,  &
      5,  0,  0,  5,  0,  0,  5,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  5,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  6,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb3 = RESHAPE( (/  1,  1,  1,  &
      1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  0,  &
      2,  2,  2,  2,  0,  2,  2,  2,  2,  0,  2,  0,  0,  0,  0,  3,  0,  0,  &
      0,  0,  3,  0,  0,  0,  0,  3,  0,  1,  3,  3,  3,  2,  3,  1,  3,  3,  &
      2,  3,  3,  1,  3,  2,  3,  0,  1,  1,  1,  0,  1,  4,  1,  1,  5,  1,  &
      1,  4,  1,  5,  1,  0,  2,  4,  4,  0,  4,  5,  4,  4,  4,  4,  4,  5,  &
      4,  4,  4,  0,  0,  2,  2,  0,  2,  3,  2,  2,  0,  2,  2,  3,  2,  7,  &
      2,  0,  0,  3,  5,  0,  3,  4,  5,  3,  0,  5,  5,  4,  5,  6,  5,  0,  &
      0,  0,  3,  0,  0,  0,  3,  0,  0,  3,  3,  0,  3,  0,  3,  0,  0,  0,  &
      4,  0,  0,  0,  6,  0,  0,  6,  6,  0,  6,  0,  6,  0,  0,  0,  0,  0,  &
      0,  0,  4,  0,  0,  4,  4,  0,  4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  &
      5,  0,  0,  5,  7,  0,  5,  0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  5,  0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      6,  0,  0,  0,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb4 = RESHAPE( (/  2,  2,  2,  &
      2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  3,  3,  3,  3,  4,  3,  &
      3,  3,  3,  4,  3,  3,  3,  3,  4,  0,  1,  1,  1,  1,  0,  1,  1,  1,  &
      1,  2,  1,  1,  1,  1,  2,  0,  2,  4,  4,  2,  4,  4,  2,  4,  4,  0,  &
      4,  4,  2,  4,  0,  0,  0,  2,  2,  0,  2,  0,  0,  2,  0,  6,  2,  2,  &
      0,  2,  6,  0,  3,  0,  0,  3,  0,  5,  5,  0,  5,  4,  0,  0,  5,  0,  &
      2,  0,  1,  3,  5,  1,  0,  1,  1,  3,  1,  2,  5,  5,  1,  5,  8,  0,  &
      0,  1,  3,  0,  0,  4,  6,  1,  4,  6,  3,  3,  6,  3,  6,  0,  0,  4,  &
      1,  0,  0,  2,  4,  4,  2,  4,  1,  1,  4,  1,  4,  0,  0,  2,  4,  0,  &
      0,  5,  5,  2,  5,  0,  6,  4,  5,  6,  8,  0,  0,  0,  2,  0,  0,  3,  &
      3,  0,  3,  0,  4,  2,  3,  4,  6,  0,  0,  0,  5,  0,  0,  0,  6,  0,  &
      0,  0,  2,  5,  6,  2,  0,  0,  0,  0,  3,  0,  0,  0,  4,  0,  0,  0,  &
      7,  3,  4,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  &
      0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  6,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  4,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  0,  0,  7,  0,  0,  7,  0,  0,  0,  0,  0,  0,  &
      0,  0,  0,  0,  0,  0,  5,  0,  0,  5,  0 /) , (/ 4, 4,20/) )
    INTEGER, DIMENSION(4, 4, 20), PARAMETER :: mb5 = RESHAPE( (/  2,  2,  2,  &
      2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  0,  3,  3,  3,  3,  4&
      ,  3,  3,  3,  3,  4,  3,  3,  3,  3,  4,  0,  0,  4,  4,  0,  0,  4,  0&
      ,  4,  4,  0,  4,  4,  0,  4,  0,  0,  1,  0,  0,  1,  2,  0,  5,  0,  0&
      ,  6,  0,  0,  5,  0,  6,  0,  0,  1,  5,  0,  0,  5,  1,  1,  5,  2,  5&
      ,  5,  1,  5,  2,  0,  0,  2,  1,  0,  0,  1,  6,  2,  1,  4,  1,  1,  6&
      ,  1,  8,  0,  0,  0,  2,  0,  0,  2,  3,  0,  2,  0,  6,  2,  3,  6,  4&
      ,  0,  0,  0,  3,  0,  0,  3,  4,  0,  3,  0,  2,  3,  4,  2,  6,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  7,  0,  0,  7,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  4,  0,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  5,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0&
      ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) , (/ 4, 4,20/) )

    INTEGER                                  :: k, k1, k2, mu
    REAL(kind=dp)                            :: cp, ct, dJ, dJc, dJcc, dJss, &
                                                dxx, dyy, f, fac1, fac2, J, &
                                                Jc, Jcc, Jss, rr, sp, st, w, &
                                                w1, w2, xx, yy, za, zb
    REAL(kind=dp), DIMENSION(3)              :: dcp, dct, dsp, dst, v
    REAL(kind=dp), DIMENSION(3, 3)           :: Arot
    REAL(kind=dp), DIMENSION(3, 3, 3)        :: dArot

dS(:,:,:)=  0.0_dp

      v(:)=R(:)
      rr=SQRT(DOT_PRODUCT(v,v))

IF( rr < 1.0e-20_dp )THEN

      DO mu=1,4
        dS(mu,mu,:)=0.0_dp
      ENDDO

ELSE

      fac1=1.0_dp
      IF(nra == 1)THEN
        fac1=fac1*2.0_dp
      ELSE
        IF(nra == 2)THEN
          fac1=fac1*SQRT(4.0_dp/3.0_dp)
        ELSE
          IF(nra == 3)THEN
            fac1=fac1*SQRT(8.0_dp/45.0_dp)
          ELSE
            IF(nra == 4)THEN
              fac1=fac1*SQRT(4.0_dp/315.0_dp)
            ELSE
              WRITE(*,*) 'nra= ',nra
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(nrb == 1)THEN
        fac1=fac1*2.0_dp
      ELSE
        IF(nrb == 2)THEN
          fac1=fac1*SQRT(4.0_dp/3.0_dp)
        ELSE
          IF(nrb == 3)THEN
            fac1=fac1*SQRT(8.0_dp/45.0_dp)
          ELSE
            IF(nrb == 4)THEN
              fac1=fac1*SQRT(4.0_dp/315.0_dp)
            ELSE
              WRITE(*,*) 'nrb= ',nrb
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      ct= -v(3)/rr
      IF(ABS(ct) >= 1.0_dp)THEN

        dct(:)=v(:)*v(3)/rr**3
        dct(3)=dct(3)-1.0_dp/rr

        Arot(1,1)=  ct
        Arot(1,2)=  0.0_dp
        Arot(1,3)=  0.0_dp
        Arot(2,1)=  0.0_dp
        Arot(2,2)=  1.0_dp
        Arot(2,3)=  0.0_dp
        Arot(3,1)=  0.0_dp
        Arot(3,2)=  0.0_dp
        Arot(3,3)=  ct

        dArot(1,1,:)= dct(:)
        dArot(1,2,:)= 0.0_dp
        dArot(1,3,:)= 0.0_dp
        dArot(2,1,:)= 0.0_dp
        dArot(2,2,:)= 0.0_dp
        dArot(2,3,:)= 0.0_dp
        dArot(3,1,:)= 0.0_dp
        dArot(3,2,:)= 0.0_dp
        dArot(3,3,:)= dct(:)

      ELSE

        xx= SQRT(v(1)**2+v(2)**2)
        st= xx/rr
        cp= -v(1)/xx
        sp= -v(2)/xx

        dct(:)=v(:)*v(3)/rr**3
        dct(3)=dct(3)-1.0_dp/rr
        dst(:)= -ct*dct(:)/st
        dcp(:)=v(:)*v(1)/(rr**3*st)
        dcp(:)=dcp(:)+v(1)*dst(:)/(rr*st**2)
        dcp(1)=dcp(1)-1.0_dp/(rr*st)
        dsp(:)=v(:)*v(2)/(rr**3*st)
        dsp(:)=dsp(:)+v(2)*dst(:)/(rr*st**2)
        dsp(2)=dsp(2)-1.0_dp/(rr*st)

        Arot(1,1)=  ct*cp
        Arot(1,2)=    -sp
        Arot(1,3)=  st*cp
        Arot(2,1)=  ct*sp
        Arot(2,2)=     cp
        Arot(2,3)=  st*sp
        Arot(3,1)= -st
        Arot(3,2)=  0.0_dp
        Arot(3,3)=  ct

        dArot(1,1,:)=  dct(:)*cp+ct*dcp(:)
        dArot(1,2,:)= -dsp(:)
        dArot(1,3,:)=  dst(:)*cp+st*dcp(:)
        dArot(2,1,:)=  dct(:)*sp+ct*dsp(:)
        dArot(2,2,:)=  dcp(:)
        dArot(2,3,:)=  dst(:)*sp+st*dsp(:)
        dArot(3,1,:)= -dst(:)
        dArot(3,2,:)=  0.0_dp
        dArot(3,3,:)=  dct(:)

      ENDIF

      za=ZSA
      zb=ZSB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)
      dxx=0.5_dp*(za+zb)
      dyy=0.5_dp*(za-zb)

      w= 0.0_dp
      w1=0.0_dp
      w2=0.0_dp
      f=rr**(nra+nrb+1)/2.0_dp**(nra+nrb+2)
      DO k=1,nc1(nra,nrb)
        w=w+REAL(c1(nra,nrb,k),dp)*AA(ma1(nra,nrb,k),xx)*BB(mb1(nra,nrb,k),yy)
        w1=w1+REAL(c1(nra,nrb,k),dp)*AA(ma1(nra,nrb,k)+1,xx)*BB(mb1(nra,nrb,k),yy)
        w2=w2+REAL(c1(nra,nrb,k),dp)*AA(ma1(nra,nrb,k),xx)*BB(mb1(nra,nrb,k)+1,yy)
      ENDDO
      J=  f*w
      dJ= f*REAL(nra+nrb+1,dp)*w/rr
      dJ=dJ-dxx*f*w1
      dJ=dJ-dyy*f*w2

      dS(1,1,:)=dS(1,1,:)+fac1*fac2*dJ*v(:)/rr

      za=ZPA
      zb=ZSB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)
      dxx=0.5_dp*(za+zb)
      dyy=0.5_dp*(za-zb)

      w= 0.0_dp
      w1=0.0_dp
      w2=0.0_dp
      f=rr**(nra+nrb+1)/2.0_dp**(nra+nrb+2)
      DO k=1,nc2(nra,nrb)
        w=w+REAL(c2(nra,nrb,k),dp)*AA(ma2(nra,nrb,k),xx)*BB(mb2(nra,nrb,k),yy)
        w1=w1+REAL(c2(nra,nrb,k),dp)*AA(ma2(nra,nrb,k)+1,xx)*BB(mb2(nra,nrb,k),yy)
        w2=w2+REAL(c2(nra,nrb,k),dp)*AA(ma2(nra,nrb,k),xx)*BB(mb2(nra,nrb,k)+1,yy)
      ENDDO
      Jc=  f*w
      dJc= f*REAL(nra+nrb+1,dp)*w/rr
      dJc=dJc-dxx*f*w1
      dJc=dJc-dyy*f*w2

      DO k1=1,3
        dS(k1+1,1,:)=dS(k1+1,1,:) &
     &          +SQRT(3.0_dp)*Arot(k1,3)*fac1*fac2*dJc*v(:)/rr &
     &          +SQRT(3.0_dp)*dArot(k1,3,:)*fac1*fac2*Jc
      ENDDO

      za=ZSA
      zb=ZPB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)
      dxx=0.5_dp*(za+zb)
      dyy=0.5_dp*(za-zb)

      w= 0.0_dp
      w1=0.0_dp
      w2=0.0_dp
      f=rr**(nra+nrb+1)/2.0_dp**(nra+nrb+2)
      DO k=1,nc3(nra,nrb)
        w=w+REAL(c3(nra,nrb,k),dp)*AA(ma3(nra,nrb,k),xx)*BB(mb3(nra,nrb,k),yy)
        w1=w1+REAL(c3(nra,nrb,k),dp)*AA(ma3(nra,nrb,k)+1,xx)*BB(mb3(nra,nrb,k),yy)
        w2=w2+REAL(c3(nra,nrb,k),dp)*AA(ma3(nra,nrb,k),xx)*BB(mb3(nra,nrb,k)+1,yy)
      ENDDO
      Jc=  f*w
      dJc= f*REAL(nra+nrb+1,dp)*w/rr
      dJc=dJc-dxx*f*w1
      dJc=dJc-dyy*f*w2

      DO k1=1,3
        dS(1,k1+1,:)=dS(1,k1+1,:) &
     &          -SQRT(3.0_dp)*Arot(k1,3)*fac1*fac2*dJc*v(:)/rr &
     &          -SQRT(3.0_dp)*dArot(k1,3,:)*fac1*fac2*Jc
      ENDDO

      za=ZPA
      zb=ZPB
      fac2=SQRT(za**(2*nra+1)*zb**(2*nrb+1))
      xx=0.5_dp*rr*(za+zb)
      yy=0.5_dp*rr*(za-zb)
      dxx=0.5_dp*(za+zb)
      dyy=0.5_dp*(za-zb)

      w= 0.0_dp
      w1=0.0_dp
      w2=0.0_dp
      f=rr**(nra+nrb+1)/2.0_dp**(nra+nrb+2)
      DO k=1,nc4(nra,nrb)
        w=w+REAL(c4(nra,nrb,k),dp)*AA(ma4(nra,nrb,k),xx)*BB(mb4(nra,nrb,k),yy)
        w1=w1+REAL(c4(nra,nrb,k),dp)*AA(ma4(nra,nrb,k)+1,xx)*BB(mb4(nra,nrb,k),yy)
        w2=w2+REAL(c4(nra,nrb,k),dp)*AA(ma4(nra,nrb,k),xx)*BB(mb4(nra,nrb,k)+1,yy)
      ENDDO
      Jss=  f*w
      dJss= f*REAL(nra+nrb+1,dp)*w/rr
      dJss=dJss-dxx*f*w1
      dJss=dJss-dyy*f*w2

      w= 0.0_dp
      w1=0.0_dp
      w2=0.0_dp
      f=rr**(nra+nrb+1)/2.0_dp**(nra+nrb+2)
      DO k=1,nc5(nra,nrb)
        w=w+REAL(c5(nra,nrb,k),dp)*AA(ma5(nra,nrb,k),xx)*BB(mb5(nra,nrb,k),yy)
        w1=w1+REAL(c5(nra,nrb,k),dp)*AA(ma5(nra,nrb,k)+1,xx)*BB(mb5(nra,nrb,k),yy)
        w2=w2+REAL(c5(nra,nrb,k),dp)*AA(ma5(nra,nrb,k),xx)*BB(mb5(nra,nrb,k)+1,yy)
      ENDDO
      Jcc=  f*w
      dJcc= f*REAL(nra+nrb+1,dp)*w/rr
      dJcc=dJcc-dxx*f*w1
      dJcc=dJcc-dyy*f*w2

      DO k1=1,3
        DO k2=1,3
          dS(k1+1,k2+1,:)=dS(k1+1,k2+1,:) &
     &            +1.5_dp*Arot(k1,1)*Arot(k2,1)*fac1*fac2*dJss*v(:)/rr &
     &            +1.5_dp*dArot(k1,1,:)*Arot(k2,1)*fac1*fac2*Jss &
     &            +1.5_dp*Arot(k1,1)*dArot(k2,1,:)*fac1*fac2*Jss &
     &            +1.5_dp*Arot(k1,2)*Arot(k2,2)*fac1*fac2*dJss*v(:)/rr &
     &            +1.5_dp*dArot(k1,2,:)*Arot(k2,2)*fac1*fac2*Jss &
     &            +1.5_dp*Arot(k1,2)*dArot(k2,2,:)*fac1*fac2*Jss &
     &            -3.0_dp*Arot(k1,3)*Arot(k2,3)*fac1*fac2*dJcc*v(:)/rr &
     &            -3.0_dp*dArot(k1,3,:)*Arot(k2,3)*fac1*fac2*Jcc &
     &            -3.0_dp*Arot(k1,3)*dArot(k2,3,:)*fac1*fac2*Jcc
        ENDDO
      ENDDO

ENDIF

END SUBROUTINE makedS

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
! *****************************************************************************
FUNCTION AA(n,x)

    INTEGER                                  :: n
    REAL(kind=dp)                            :: x, AA

    REAL(kind=dp)                            :: p

IF(n == 0)THEN
  p=     1.0_dp
ELSE
IF(n == 1)THEN
  p=     1.0_dp+x
ELSE
IF(n == 2)THEN
  p=     2.0_dp+x*(&
         2.0_dp+x)
ELSE
IF(n == 3)THEN
  p=     6.0_dp+x*(&
         6.0_dp+x*(&
         3.0_dp+x))
ELSE
IF(n == 4)THEN
  p=    24.0_dp+x*(&
        24.0_dp+x*(&
        12.0_dp+x*(&
         4.0_dp+x)))
ELSE
IF(n == 5)THEN
  p=   120.0_dp+x*(&
       120.0_dp+x*(&
        60.0_dp+x*(&
        20.0_dp+x*(&
         5.0_dp+x))))
ELSE
IF(n == 6)THEN
  p=   720.0_dp+x*(&
       720.0_dp+x*(&
       360.0_dp+x*(&
       120.0_dp+x*(&
        30.0_dp+x*(&
         6.0_dp+x)))))
ELSE
IF(n == 7)THEN
  p=  5040.0_dp+x*(&
      5040.0_dp+x*(&
      2520.0_dp+x*(&
       840.0_dp+x*(&
       210.0_dp+x*(&
        42.0_dp+x*(&
         7.0_dp+x))))))
ELSE
IF(n == 8)THEN
  p= 40320.0_dp+x*(&
     40320.0_dp+x*(&
     20160.0_dp+x*(&
      6720.0_dp+x*(&
      1680.0_dp+x*(&
       336.0_dp+x*(&
        56.0_dp+x*(&
         8.0_dp+x)))))))
ELSE
IF(n == 9)THEN
  p= 362880.0_dp+x*(&
     362880.0_dp+x*(&
     181440.0_dp+x*(&
      60480.0_dp+x*(&
      15120.0_dp+x*(&
       3024.0_dp+x*(&
        504.0_dp+x*(&
         72.0_dp+x*(&
          9.0_dp+x))))))))
ELSE
IF(n == 10)THEN
  p= 3628800.0_dp+x*(&
     3628800.0_dp+x*(&
     1814400.0_dp+x*(&
      604800.0_dp+x*(&
      151200.0_dp+x*(&
       30240.0_dp+x*(&
        5040.0_dp+x*(&
         720.0_dp+x*(&
          90.0_dp+x*(&
          10.0_dp+x)))))))))
ELSE
  p=1.0_dp
  WRITE(*,*) ' n= ',n,' in AA(n,x) '
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF
ENDIF

AA=EXP(-x)*p/x**(n+1)

END FUNCTION AA

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param y ...
! *****************************************************************************
FUNCTION BB(n,y)

    INTEGER                                  :: n
    REAL(kind=dp)                            :: y, BB

IF( ABS(y) > 1.0e-20_dp )THEN
  BB= REAL((-1)**(n+1),dp)*AA(n,-y)-AA(n,y)
ELSE
  IF( MOD(n,2) == 0 )THEN
    BB=    2.0_dp/REAL(n+1,dp)
  ELSE
    BB= -y*2.0_dp/REAL(n+2,dp)
  ENDIF
ENDIF

END FUNCTION BB

END MODULE se_core_matrix

