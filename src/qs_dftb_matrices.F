!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_dftb_matrices [1.1] *
!!
!!   FUNCTION
!!     Calculation of Overlap and Hamiltonian matrices in DFTB
!!
!!   AUTHOR
!!     JGH
!!
!!***
!******************************************************************************

MODULE qs_dftb_matrices

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type,&
                                             dftb_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_charges
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE particle_types,                  ONLY: particle_type
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                             qs_dftb_pairpot_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_dftb_coulomb,                 ONLY: build_dftb_coulomb
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             set_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_output,            ONLY: write_sparse_matrix
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix_set,&
                                             allocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             sparse_times_local,&
                                             cp_sm_scale_and_add,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  INTEGER,DIMENSION(16),PARAMETER        :: orbptr = (/ 0, 1, 1, 1, &
                                   2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3 /)

  ! Maximum number of points used for interpolation
  INTEGER, PARAMETER                     :: max_inter = 5
  ! Maximum number of points used for extrapolation
  INTEGER, PARAMETER                     :: max_extra = 9
  ! see also qs_dftb_parameters
  REAL(dp), PARAMETER                    :: slako_d0 = 1._dp
  ! pointer to skab
  INTEGER, DIMENSION(0:3,0:3,0:3,0:3,0:3):: iptr
  ! screening for gamma function
  REAL(dp), PARAMETER                    :: tol_gamma = 1.e-4_dp
  ! small real number
  REAL(dp), PARAMETER                    :: rtiny = TINY(1._dp)

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_matrices'

  PUBLIC :: build_dftb_matrices, build_dftb_ks_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_dftb_matrices(qs_env,para_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dftb_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nkind, ikind, jkind, &
      natom, natorb_a, natorb_b, iab, nlist, ilist, iatom, jatom, nnode, nmat, &
      inode, iw, istat, ik, jk, ngrd, ngrdcut, llm, last_jatom, i, icol, irow, &
      lmaxi, lmaxj, l1, l2, m, la, lb, nrow, ncol, neighbor_list_id, n_urpoly, &
      spdim, atom_a, atom_b, n1, n2
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, felem, lelem
    LOGICAL                                  :: defined, failure, new_atom_b
    REAL(KIND=dp)                            :: dr, dgrd, ga, gb, urep_cut, s_cut,&
                                                erep, ddr, drp, drm
    REAL(KIND=dp), DIMENSION(0:3)            :: skself, eta_a, eta_b
    REAL(KIND=dp), DIMENSION(2)              :: surr
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, force_w,rij, srep
    REAL(KIND=dp), DIMENSION(10)             :: urep
    REAL(KIND=dp), DIMENSION(20)             :: skab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fmat, smat, pblock, wblock,&
                                                sblock, fblock, gblock, dgblock, dsblock, &
                                                dfblock, dsblock_i, spxr, scoeff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind_a, dftb_kind_b
    TYPE(qs_dftb_pairpot_type), POINTER      :: dftb_param_ij
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s,&
                                                matrix_w, gamma_matrix
    TYPE(distribution_2d_type), POINTER      :: distribution_2d

!   ----------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    ! set pointers
    iptr = 0
    DO la=0,3
      DO lb=0,3
        llm=0
        DO l1=0,MAX(la,lb)
          DO l2=0,MIN(l1,la,lb)
            DO m=0,l2
              llm=llm+1
              iptr(l1,l2,m,la,lb)=llm
            END DO
          END DO
        END DO
      END DO
    END DO

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! Allocate the overlap and Hamiltonian matrix
    CALL setup_matrices(qs_env,calculate_forces,error)

    NULLIFY (matrix_h,matrix_s,matrix_p,atomic_kind_set,sab_orb)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    matrix_h=matrix_h,matrix_s=matrix_s,&
                    dft_control=dft_control,error=error)

    dftb_control => dft_control%qs_control%dftb_control

    ! gamma matrix allocation
    IF ( dftb_control%self_consistent ) THEN
      IF(calculate_forces) THEN
        nmat=4
      ELSE
        nmat=1
      END IF
      NULLIFY(gamma_matrix)
      CALL allocate_matrix_set(gamma_matrix,nmat,error=error)
      NULLIFY(distribution_2d)
      CALL get_qs_env(qs_env=qs_env,&
          particle_set=particle_set,&
          neighbor_list_id=neighbor_list_id,&
          distribution_2d=distribution_2d,error=error)
      natom = SIZE(particle_set)
      nrow = natom
      ncol = natom
      ALLOCATE (felem(natom),lelem(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO iatom = 1, natom
        felem(iatom) = iatom
        lelem(iatom) = iatom
      ENDDO
      CALL allocate_matrix(matrix=gamma_matrix(1)%matrix,&
                nblock_row=natom,nblock_col=natom,&
                nrow=nrow,ncol=ncol,&
                first_row=felem,last_row=lelem,&
                first_col=felem,last_col=lelem,&
                matrix_name="GAMMA MATRIX",&
                matrix_symmetry="symmetric",&
                sparsity_id=neighbor_list_id, &
                distribution_2d=distribution_2d,error=error)
      DO i=2,nmat
        CALL replicate_matrix_structure(source=gamma_matrix(1)%matrix,&
                                    target=gamma_matrix(i)%matrix,&
                                    target_name="DERIVATIVE GAMMA MATRIX",&
                                    target_symmetry="symmetric",error=error)
      END DO
      dftb_control%gamma_matrix => gamma_matrix
      DEALLOCATE (felem,lelem,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF(calculate_forces) THEN
      NULLIFY (rho,rho,force,particle_set,matrix_w)
      CALL get_qs_env(qs_env=qs_env,&
                      particle_set=particle_set,&
                      rho=rho,matrix_w=matrix_w,&
                      force=force,error=error)
      matrix_p => rho%rho_ao
      IF (SIZE(matrix_p) == 2) THEN
         CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,&
              matrix_p(2)%matrix,1.0_dp,error=error)
         CALL cp_sm_scale_and_add(matrix_w(1)%matrix, 1.0_dp,&
              matrix_w(2)%matrix,1.0_dp,error=error)
      END IF
      natom = SIZE (particle_set)
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
    END IF

    erep = 0._dp

    CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            natom=natom,&
                            dftb_parameter=dftb_kind_a)
       CALL get_dftb_atom_param(dftb_kind_a,&
              defined=defined,lmax=lmaxi,skself=skself,&
              eta=eta_a,natorb=natorb_a)

       IF (.NOT.defined .OR. natorb_a < 1) CYCLE

       DO jkind=1,nkind
           atomic_kind => atomic_kind_set(jkind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                dftb_parameter=dftb_kind_b)
           CALL get_dftb_atom_param(dftb_kind_b,&
                  defined=defined,lmax=lmaxj,eta=eta_b,natorb=natorb_b)

           IF (.NOT.defined .OR. natorb_b < 1) CYCLE

           ! retrieve information on F and S matrix
           ik = MAX(ikind,jkind)
           jk = MIN(ikind,jkind)
           dftb_param_ij => dftb_control%dftb_potential(ik,jk)
           ngrd = dftb_param_ij%ngrd
           ngrdcut = dftb_param_ij%ngrdcut
           dgrd = dftb_param_ij%dgrd
           ddr = dgrd*0.1_dp
           llm = dftb_param_ij%llm
           fmat => dftb_param_ij%fmat
           smat => dftb_param_ij%smat
           ! repulsive pair potential
           n_urpoly = dftb_param_ij%n_urpoly
           urep_cut = dftb_param_ij%urep_cut
           urep = dftb_param_ij%urep
           spxr => dftb_param_ij%spxr
           scoeff => dftb_param_ij%scoeff
           spdim = dftb_param_ij%spdim
           s_cut = dftb_param_ij%s_cut
           srep = dftb_param_ij%srep
           surr = dftb_param_ij%surr

           iab = ikind + nkind*(jkind - 1)
           IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_orb(iab)%neighbor_list_set,&
                                      nlist=nlist)
           NULLIFY(sab_orb_neighbor_list)

           DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
               sab_orb_neighbor_list => &
                          first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
               sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                    atom=iatom,&
                                    nnode=nnode)
             last_jatom = 0

             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

             DO inode=1,nnode
               CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                      neighbor=jatom,r=rij)

               dr = SQRT(SUM(rij(:)**2))
               IF (NINT(dr/dgrd) <= ngrdcut) THEN

                 IF (jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                 ELSE
                   new_atom_b = .FALSE.
                 END IF
                 IF (new_atom_b) THEN
                   icol = MAX(iatom,jatom)
                   irow = MIN(iatom,jatom)
                   NULLIFY (sblock)
                   CALL add_block_node(matrix=matrix_s(1)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=sblock,error=error)
                   NULLIFY (fblock)
                   CALL add_block_node(matrix=matrix_h(1)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=fblock,error=error)
                   IF ( dftb_control%self_consistent ) THEN
                     NULLIFY (gblock)
                     CALL add_block_node(matrix=gamma_matrix(1)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=gblock,error=error)
                   END IF
                   IF (calculate_forces) THEN
                     NULLIFY (pblock)
                     CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=pblock)
                     CPPrecondition(ASSOCIATED(pblock),cp_failure_level,routineP,error,failure)
                     NULLIFY (wblock)
                     CALL get_block_node(matrix=matrix_w(1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=wblock)
                     CPPrecondition(ASSOCIATED(wblock),cp_failure_level,routineP,error,failure)
                     force_ab = 0._dp
                     force_w  = 0._dp
                     IF ( dftb_control%self_consistent ) THEN
                       DO i=2,4
                         NULLIFY (dsblock)
                         CALL add_block_node(matrix=matrix_s(i)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=dsblock,error=error)
                         NULLIFY (dgblock)
                         CALL add_block_node(matrix=gamma_matrix(i)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=dgblock,error=error)
                       END DO
                     END IF
                   END IF
                 END IF

                 IF (iatom == jatom .AND. dr < 0.001_dp) THEN
                   ! diagonal block
                   DO i=1,natorb_a
                     sblock(i,i) = sblock(i,i) + 1._dp
                     fblock(i,i) = fblock(i,i) + skself(orbptr(i))
                   END DO
                 ELSE
                   ! off-diagonal block
                   IF ( irow == iatom ) THEN
                     CALL getskz(smat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                     CALL turnsk(sblock,skab,rij,dr,lmaxi,lmaxj)
                     CALL getskz(fmat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                     CALL turnsk(fblock,skab,rij,dr,lmaxi,lmaxj)
                   ELSE
                     CALL getskz(smat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                     CALL turnsk(sblock,skab,-rij,dr,lmaxj,lmaxi)
                     CALL getskz(fmat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                     CALL turnsk(fblock,skab,-rij,dr,lmaxj,lmaxi)
                   END IF
                   IF(calculate_forces) THEN
                     n1 = SIZE(fblock,1)
                     n2 = SIZE(fblock,2)
                     ALLOCATE (dfblock(n1,n2),dsblock(n1,n2),STAT=istat)
                     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                     dfblock=0._dp
                     dsblock=0._dp
                     DO i=1,3
                       rij(i) = rij(i) + ddr
                       dr = SQRT(SUM(rij(:)**2))
                       IF ( irow == iatom ) THEN
                         CALL getskz(smat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dsblock,skab,rij,dr,lmaxi,lmaxj)
                         CALL getskz(fmat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dfblock,skab,rij,dr,lmaxi,lmaxj)
                       ELSE
                         CALL getskz(smat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dsblock,skab,-rij,dr,lmaxj,lmaxi)
                         CALL getskz(fmat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dfblock,skab,-rij,dr,lmaxj,lmaxi)
                       END IF
                       rij(i) = rij(i) - 2._dp*ddr
                       dr = SQRT(SUM(rij(:)**2))
                       dsblock = -dsblock
                       dfblock = -dfblock
                       IF ( irow == iatom ) THEN
                         CALL getskz(smat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dsblock,skab,rij,dr,lmaxi,lmaxj)
                         CALL getskz(fmat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dfblock,skab,rij,dr,lmaxi,lmaxj)
                       ELSE
                         CALL getskz(smat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dsblock,skab,-rij,dr,lmaxj,lmaxi)
                         CALL getskz(fmat,skab,dr,ngrd,ngrdcut,dgrd,llm)
                         CALL turnsk(dfblock,skab,-rij,dr,lmaxj,lmaxi)
                       END IF
                       rij(i) = rij(i) + ddr
                       dr = SQRT(SUM(rij(:)**2))
                       dfblock = -0.5_dp*dfblock/ddr
                       dsblock = -0.5_dp*dsblock/ddr
                       force_ab(i) = force_ab(i) + SUM(dfblock*pblock)
                       force_w(i) = force_w(i) - SUM(dsblock*wblock)
                       IF ( dftb_control%self_consistent ) THEN
                         NULLIFY (dsblock_i)
                         CALL get_block_node(matrix=matrix_s(1+i)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=dsblock_i)
                         CPPrecondition(ASSOCIATED(dsblock),cp_failure_level,routineP,error,failure)
                         dsblock_i = dsblock_i + dsblock
                       END IF
                     ENDDO
                     DEALLOCATE (dfblock,dsblock,STAT=istat)
                     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                 END IF

                 IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
                   atom_a = atom_of_kind(iatom)
                   atom_b = atom_of_kind(jatom)
                   IF ( irow == iatom ) force_ab = -force_ab
                   IF ( irow == iatom ) force_w = -force_w
                   force(ikind)%all_potential(:,atom_a) =&
                       force(ikind)%all_potential(:,atom_a) - force_ab(:)
                   force(jkind)%all_potential(:,atom_b) =&
                       force(jkind)%all_potential(:,atom_b) + force_ab(:)
                   force(ikind)%overlap(:,atom_a) =&
                       force(ikind)%overlap(:,atom_a) - force_w(:)
                   force(jkind)%overlap(:,atom_b) =&
                       force(jkind)%overlap(:,atom_b) + force_w(:)
                 END IF

                 ! gamma matrix
                 IF ( dftb_control%self_consistent ) THEN
                   ga = eta_a(0)
                   gb = eta_b(0)
                   gblock(1,1)= gblock(1,1) + gamma_rab_sr(dr,ga,gb)
                   IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
                     DO i=1,3
                       rij(i) = rij(i) + ddr
                       drp = SQRT(SUM(rij(:)**2))
                       rij(i) = rij(i) - 2._dp*ddr
                       drm = SQRT(SUM(rij(:)**2))
                       rij(i) = rij(i) + ddr
                       NULLIFY (dgblock)
                       CALL get_block_node(matrix=gamma_matrix(1+i)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=dgblock)
                       CPPrecondition(ASSOCIATED(dgblock),cp_failure_level,routineP,error,failure)
                       dgblock(1,1)= dgblock(1,1) + &
                         0.5_dp*(gamma_rab_sr(drp,ga,gb)-gamma_rab_sr(drm,ga,gb))/ddr
                     END DO
                   END IF
                 END IF

               END IF

               ! repulsive potential
               IF (dr <= urep_cut .OR. spdim > 0) THEN
                  CALL urep_egr(rij,dr,erep,force_ab,&
                      n_urpoly,urep,spdim,s_cut,srep,spxr,scoeff,surr,calculate_forces)
                  IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
                    atom_a = atom_of_kind(iatom)
                    atom_b = atom_of_kind(jatom)
                    force(ikind)%repulsive(:,atom_a) =&
                        force(ikind)%repulsive(:,atom_a) - force_ab(:)
                    force(jkind)%repulsive(:,atom_b) =&
                        force(jkind)%repulsive(:,atom_b) + force_ab(:)
                 END IF
               END IF

               sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO ! inode => jatom(atom B)

           END DO ! ilist => iatom(atom A)

        END DO ! jkind

     END DO ! ikind

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,para_env,&
             output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", error=error)
     END IF

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,para_env,&
             output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/OVERLAP", error=error)
     END IF

     IF(calculate_forces) THEN
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix,1.0_dp,&
               matrix_p(2)%matrix,-1.0_dp,error=error)
       END IF
       DEALLOCATE(atom_of_kind)
     END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_dftb_matrices

! *****************************************************************************

  SUBROUTINE build_dftb_ks_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins,&
      output_unit, ncol_global, natom, istat, nkind, ikind, iatom, atom_a
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers, &
                                                scaling_factor, mcharge
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: charges
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s
    TYPE(section_vals_type), POINTER         :: scf_section

!   ----------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(matrix_h, dft_control, logger, scf_section)
    NULLIFY(matrix_s, matrix_p, particle_set)

    energy%hartree = 0._dp

    logger => cp_error_get_logger(error)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            atomic_kind_set=atomic_kind_set,&
            matrix_h=matrix_h,&
            para_env=para_env,error=error)

       scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF",error=error)
       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

       DO ispin=1,nspins
         ! copy the core matrix into the fock matrix
         CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix,error=error)
       END DO

       ! Mulliken charges
       CALL get_qs_env(qs_env=qs_env,&
            particle_set=particle_set,&
            matrix_s=matrix_s,error=error)
       matrix_p => rho%rho_ao
       natom=SIZE(particle_set)
       ALLOCATE(charges(natom,nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !
       CALL mulliken_charges(matrix_p,matrix_s(1)%matrix,para_env,charges)
       !
       ALLOCATE(mcharge(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nkind = SIZE(atomic_kind_set)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             natom=natom,dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind,zeff=zeff)
          DO iatom=1,natom
            atom_a = atomic_kind%atom_list(iatom)
            mcharge(atom_a) = zeff - SUM(charges(atom_a,1:nspins))
         END DO
       END DO
       DEALLOCATE(charges,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       CALL build_dftb_coulomb(ks_env,qs_env,ks_matrix,rho,mcharge,energy,&
               calculate_forces,just_energy,error)

       DEALLOCATE(mcharge,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       energy%hartree = energy%hartree + energy%core_overlap + energy%core_self
       energy%total = energy%core + energy%hartree + energy%qmmm_el 

       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T9,A,T60,F20.10))")&
               "Zeroth order Hamiltonian energy:               ",energy%core,&
               "Charge fluctuation energy:                     ",energy%hartree
          IF (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DETAILED_ENERGY", error=error)
       ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,ncol_global=ncol_global,error=error)
             CALL cp_sm_fm_multiply(ks_matrix(ispin)%matrix,mo_coeff,&
                                          mo_derivs(ispin)%matrix, ncol_global,error=error)
             ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
             scaling_factor=2.0_dp*occupation_numbers
             CALL cp_fm_column_scale(mo_derivs(ispin)%matrix,scaling_factor)
             DEALLOCATE(scaling_factor)
          ENDDO
       ENDIF

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_dftb_ks_matrix

! *****************************************************************************

  SUBROUTINE setup_matrices(qs_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'setup_matrices', routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: natom, nsgf, neighbor_list_id,&
      nkind, istat, nmat, i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_h

!   ----------------------------------------------------------------------------

    logger => cp_error_get_logger(error)
    NULLIFY(logger,matrix_s,matrix_h,particle_set,sab_orb,distribution_2d,&
            atomic_kind_set)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    neighbor_list_id=neighbor_list_id, &
                    distribution_2d=distribution_2d,error=error)

    dftb_control => dft_control%qs_control%dftb_control

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)

    IF ( dftb_control%self_consistent .AND. calculate_forces ) THEN
      ! we have to store the derivative overlap matrices
      nmat = 4
    ELSE
      nmat = 1
    END IF

    CALL allocate_matrix_set(matrix_s,nmat,error=error)
    CALL allocate_matrix_set(matrix_h,1,error=error)

    CALL allocate_matrix(matrix=matrix_s(1)%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=first_sgf,&
                         last_row=last_sgf,&
                         first_col=first_sgf,&
                         last_col=last_sgf,&
                         matrix_name="OVERLAP MATRIX",&
                         matrix_symmetry="symmetric",&
                         sparsity_id=neighbor_list_id, &
                         distribution_2d=distribution_2d,error=error)

    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                    target=matrix_h(1)%matrix,&
                                    target_name="CORE HAMILTONIAN MATRIX",&
                                    target_symmetry="symmetric",error=error)
    IF ( nmat > 1 ) THEN
      DO i=2,nmat
        CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                        target=matrix_s(i)%matrix,&
                                        target_name="OVERLAP DERIVATIVE MATRIX",&
                                        target_symmetry="symmetric",error=error)
      END DO
    END IF

    CALL set_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
    CALL set_qs_env(qs_env=qs_env,matrix_h=matrix_h,error=error)

  END SUBROUTINE setup_matrices

! *****************************************************************************
! module dftb_sktools
!
! Creation: 07. Feb. 2004, TH
!
  subroutine getskz(slakotab,skpar,dx,ngrd,ngrdcut,dgrd,llm)
    !
    ! Gets matrix elements on z axis, as they are stored in the tables
    !
    integer, intent(in)   :: llm, ngrd, ngrdcut
    real(dp), intent(out) :: skpar(llm)
    real(dp), intent(in)  :: dx,dgrd,slakotab(:,:)
    !
    integer               :: clgp
    !
    skpar = 0._dp
    !
    ! Determine closest grid point
    !
    clgp = nint(dx/dgrd)
    !
    ! Screen elements which are too far away
    !
    if (clgp > ngrdcut) return
    !
    ! The grid point is either contained in the table --> matrix element 
    ! can be interpolated, or it is outside the table --> matrix element 
    ! needs to be extrapolated.
    !
    if (clgp > ngrd) then
      ! 
      ! Extrapolate external matrix elements if table does not finish with zero
      !
      call extrapol(slakotab,skpar,dx,ngrd,dgrd,llm)
    else
      !
      ! Interpolate tabulated matrix elements
      !
      call interpol(slakotab,skpar,dx,ngrd,dgrd,llm,clgp)
    end if
  end subroutine getskz

  subroutine interpol(slakotab,skpar,dx,ngrd,dgrd,llm,clgp)
    integer, intent(in)    :: clgp                 ! closest grid point
    integer, intent(in)    :: ngrd                 ! number of grid points
    real(dp), intent(in)   :: dgrd                 ! grid spacing
    integer, intent(in)    :: llm                  ! dimension of skpar
    real(dp), intent(out)  :: skpar(llm)           ! interpolated Slater-Koster values
    real(dp), intent(in)   :: slakotab(:,:)        ! Slater-Koster table
    real(dp), intent(in)   :: dx                   ! diatomic distance
    !
    integer                :: l,k !loop variables
    real(dp)               :: xa(max_inter),ya(max_inter),error
    integer                :: fgpm,lgpm ! first and last grid points -1
    !
    ! Get grid points for interpolation
    !
    lgpm = min(clgp+max_inter/2,ngrd) - 1
    fgpm = lgpm - max_inter + 1
    do k = 1,max_inter
      xa(k) = (fgpm+k)*dgrd
    end do
    !
    ! Interpolate matrix elements for all orbitals
    !
    do l = 1, llm
      !
      ! Read SK parameters from table
      !
      ya(1:max_inter) = slakotab(fgpm:lgpm,l)
      call polint(xa,ya,max_inter,dx,skpar(l),error)
    end do
  end subroutine interpol

  subroutine extrapol(slakotab,skpar,dx,ngrd,dgrd,llm)
    integer, intent(in)    :: llm ! dimension of Slako vector
    real(dp), intent(out)  :: skpar(llm) ! extrapolated Slater-Koster value
    real(dp), intent(in)   :: slakotab(:,:) ! Slater-Koster table
    integer, intent(in)    :: ngrd ! number of grid points
    real(dp), intent(in)   :: dgrd ! grid spacing
    real(dp), intent(in)   :: dx ! diatomic distance
    !
    integer               :: l,k !loop variables
    real(dp)              :: xa(max_extra),ya(max_extra)
    real(dp)              :: error
    integer               :: fgp,lgp,nzero,ntable
    !
    ! Get the last three grid points of the table
    !
    nzero = max_extra/3
    ntable = max_extra-nzero
    !
    ! Get the three last distances from the table
    !
    do k = 1,ntable
      xa(k) = (ngrd-(max_extra-3)+k)*dgrd
    end do
    do k = 1,nzero
      xa(ntable+k) = (ngrd+k-1)*dgrd + slako_d0
      ya(ntable+k) = 0.0
    end do
    !
    ! Extrapolate matrix elements for all orbitals
    !
    do l = 1,llm
      !
      ! Read SK parameters from table
      !
      fgp = ngrd + 1 - (max_extra-3)
      lgp = ngrd
      ya(1:max_extra-3) = slakotab(fgp:lgp,l)
      call polint(xa,ya,max_extra,dx,skpar(l),error)
    end do
  end subroutine extrapol  

  subroutine turnsk(mat,skab,dxv,dx,lmaxa,lmaxb)
    !
    ! Purpose: Turn matrix element from z-axis to orientation of dxv
    !
    ! History: Creation 13. Jan 2004, TH
    !
    implicit none
    !
    real(dp), intent(in)    :: dxv(3),dx    ! distance vector and distance
    real(dp), intent(in)    :: skab(:)      ! Integrals on z-axis
    real(dp), intent(inout) :: mat(:,:)     ! final integrals
    integer, intent(in)     :: lmaxa, lmaxb ! max angular momentum
    !
    ! These routines are taken from an old TB code (unknown to TH). 
    ! They are highly optimised and taken because they are time critical. 
    ! They are explicit, so not recursive, and work up to d functions.
    ! 
    ! Set variables necessary for rotation of matrix elements
    !
    real(dp) :: rr2(6) ! r_i^2/r, replicated in rr2(4:6) for index convenience later
    real(dp) :: rr(6)  ! r_i/r, direction vector, rr(4:6) are replicated from 1:3 
    real(dp) :: rinv
    integer  :: lmaxab ! lmax of A and B
    integer  :: minlmaxab ! min of lmaxa and lmaxb
    !
    ! Determine l quantum limits. 
    !
    lmaxab = max(lmaxa,lmaxb)
    if (lmaxab.gt.2)  stop 'lmax=2'
    minlmaxab = min(lmaxa,lmaxb)
    !
    ! s-s interaction
    !
    call skss(skab,mat)
    !
    if (lmaxab.le.0) return
    !
    rr2(1:3) = dxv(1:3)**2
    rr(1:3) = dxv(1:3)
    rinv = 1.0_dp/dx
    !
    rr(1:3) = rr(1:3)*rinv
    rr(4:6) = rr(1:3)
    rr2(1:3) = rr2(1:3)*rinv**2
    rr2(4:6) = rr2(1:3)
    !
    ! s-p, p-s and p-p interaction
    !
    if (minlmaxab.ge.1) then
      call skpp(skab,mat,iptr(:,:,:,lmaxa,lmaxb))
      call sksp(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
      call sksp(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
    else
      if (lmaxb.ge.1) then
        call sksp(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
      else
        call sksp(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
      end if
    end if
    ! 
    ! If there is only s-p interaction we have finished
    !
    if (lmaxab.le.1) return
    !
    ! at least one atom has d functions
    !
    if (minlmaxab.eq.2) then
      !
      ! in case both atoms have d functions
      !
      call skdd(skab,mat,iptr(:,:,:,lmaxa,lmaxb))
      call sksd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
      call sksd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
      call skpd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
      call skpd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
    else
      !
      ! One atom has d functions, the other has s or s and p functions
      !
      if (lmaxa.eq.0) then
        !
        ! atom b has d, the atom a only s functions
        !
        call sksd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
      else if (lmaxa.eq.1) then
        !
        ! atom b has d, the atom a s and p functions
        !
        call sksd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
        call skpd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.true.)
      else
        !
        ! atom a has d functions
        !
        if (lmaxb.eq.0) then
          !
          ! atom a has d, atom b has only s functions
          !
          call sksd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
        else
          !
          ! atom a has d, atom b has s and p functions
          !
          call sksd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
          call skpd(skab,mat,iptr(:,:,:,lmaxa,lmaxb),.false.)
        end if
      end if
    end if
    !
  contains
    !
    ! The subroutines to turn the matrix elements are taken as internal subroutines
    ! as it is beneficial to inline them.
    !
    ! They are both turning the matrix elements and placing them appropriately 
    ! into the matrix block
    !
    subroutine skss(skpar,mat)
      !
      ! Purpose: s-s interaction (no rotation necessary)
      !
      implicit none
      !
      real(dp), intent(inout) :: mat(:,:)
      real(dp), intent(in) :: skpar(:)
      !
      mat(1,1) = mat(1,1) + skpar(1)
      !
    end subroutine skss
    
    subroutine sksp(skpar,mat,ind,transposed)
      !
      ! s-p interaction (simple rotation)
      !
      implicit none
      !
      logical, intent(in)     :: transposed
      integer, intent(in)     :: ind(0:,0:,0:)
      real(dp), intent(inout) :: mat(:,:)
      real(dp), intent(in)    :: skpar(:)
      !
      integer                 :: l
      real(dp)                :: skp
      !
      skp = skpar(ind(1,0,0))
      if (transposed) then
        do l = 1,3
          mat(1,l+1) = mat(1,l+1) + rr(l)*skp
        end do
      else
        do l = 1,3
          mat(l+1,1) = mat(l+1,1) - rr(l)*skp
        end do
      end if
      !
    end subroutine sksp
    !
    !
    subroutine skpp(skpar,mat,ind)
      !
      implicit none
      !
      real(dp), intent(in)    :: skpar(:)
      real(dp), intent(inout) :: mat(:,:)
      integer, intent(in)     :: ind(0:,0:,0:)
      !
      real(dp)                :: matel(6),epp(6)
      integer                 :: l,ir,is,ii,k
      real(dp)                :: skppp,skpps
      !
      epp(1:3) = rr2(1:3)
      do l = 1,3
        epp(l+3) = rr(l)*rr(l+1)
      end do
      skppp = skpar(ind(1,1,1))
      skpps = skpar(ind(1,1,0))
      !
      do l = 1,3
        matel(l) = epp(l)*skpps + (1._dp-epp(l))*skppp
      end do
      do l = 4,6
        matel(l) = epp(l)*(skpps - skppp)
      end do
      !
      do ir = 1,3
        do is = 1,ir
          ii = ir - is
          k = 3*ii-(ii*(ii-1))/2+is
          mat(is+1,ir+1) = mat(is+1,ir+1) + matel(k)
          mat(ir+1,is+1) = mat(ir+1,is+1) + matel(k)
        end do
      end do
    end subroutine skpp
    !
    !
    subroutine sksd(skpar,mat,ind,transposed)
      !
      implicit none
      !
      logical, intent(in)     :: transposed
      integer, intent(in)     :: ind(0:,0:,0:)
      real(dp), intent(inout) :: mat(:,:)
      real(dp), intent(in)    :: skpar(:)
      !
      integer                 :: l
      real(dp)                :: es(5),d4,d5,r3
      real(dp)                :: sksds
      !
      sksds = skpar(ind(2,0,0))
      r3 = sqrt(3._dp)
      d4 = rr2(3) - 0.5_dp*(rr2(1)+rr2(2))
      d5 = rr2(1) - rr2(2)
      !
      do l = 1,3
        es(l) = r3*rr(l)*rr(l+1)
      end do
      es(4) = 0.5_dp*r3*d5
      es(5) = d4
      !
      if (transposed) then
        do l = 1,5
          mat(1,l+4) = mat(1,l+4) + es(l)*sksds
        end do
      else
        do l = 1,5
          mat(l+4,1) = mat(l+4,1) + es(l)*sksds
        end do
      end if
    end subroutine sksd
    !
    !
    subroutine skpd(skpar,mat,ind,transposed)
      !
      implicit none
      !
      logical, intent(in)     :: transposed
      integer, intent(in)     :: ind(0:,0:,0:)
      real(dp), intent(inout) :: mat(:,:)
      real(dp), intent(in)    :: skpar(:)
      !
      real(dp)                :: epd(13,2),dm(15)
      real(dp)                :: r3,d3,d4,d5,d6,sktmp
      integer                 :: l,ir,is,k,m
      !
      r3 = sqrt(3.0_dp)
      d3 = rr2(1) + rr2(2)
      d4 = rr2(3) - 0.5_dp*d3
      d5 = rr2(1) - rr2(2)
      d6 = rr(1)*rr(2)*rr(3)
      do l = 1,3
        epd(l,1) = r3*rr2(l)*rr(l+1)
        epd(l,2) = rr(l+1)*(1.0_dp-2._dp*rr2(l))
        epd(l+4,1) = r3*rr2(l)*rr(l+2)
        epd(l+4,2) = rr(l+2)*(1.0_dp-2*rr2(l))
        epd(l+7,1) = 0.5_dp*r3*rr(l)*d5
        epd(l+10,1) = rr(l)*d4
      end do
      !
      epd(4,1) = r3*d6
      epd(4,2) = -2._dp*d6
      epd(8,2) = rr(1)*(1.0_dp-d5)
      epd(9,2) = -rr(2)*(1.0_dp+d5)
      epd(10,2) = -rr(3)*d5
      epd(11,2) = -r3*rr(1)*rr2(3)
      epd(12,2) = -r3*rr(2)*rr2(3)
      epd(13,2) = r3*rr(3)*d3
      !
      dm(1:15) = 0.0_dp
      !
      do m = 1,2
        sktmp = skpar(ind(2,1,m-1))
        dm(1)=dm(1)+epd(1,m)*sktmp
        dm(2)=dm(2)+epd(6,m)*sktmp
        dm(3)=dm(3)+epd(4,m)*sktmp
        dm(5)=dm(5)+epd(2,m)*sktmp
        dm(6)=dm(6)+epd(7,m)*sktmp
        dm(7)=dm(7)+epd(5,m)*sktmp
        dm(9)=dm(9)+epd(3,m)*sktmp
        do l = 8,13
          dm(l+2) = dm(l+2)+epd(l,m)*sktmp
        end do
      end do
      !
      dm(4) = dm(3)
      dm(8) = dm(3)
      !
      if (transposed) then
        do ir = 1,5
          do is = 1,3
            k=3*(ir-1)+is
            mat(is+1,ir+4) = mat(is+1,ir+4) + dm(k)
          end do
        end do
      else
        do ir = 1,5
          do is = 1,3
            k=3*(ir-1)+is
            mat(ir+4,is+1) = mat(ir+4,is+1) - dm(k)
          end do
        end do
      end if
      !
    end subroutine skpd
    !
    !
    subroutine skdd(skpar,mat,ind)
      !
      implicit none
      !
      integer, intent(in)     :: ind(0:,0:,0:)
      real(dp), intent(inout) :: mat(:,:)
      real(dp), intent(in)    :: skpar(:)
      !
      real(dp)                :: e(15,3),dm(15),dd(3),r3,d3,d4,d5
      integer                 :: ii,ir,is,l,m,k
      !
      r3 = sqrt(3._dp)
      d3 = rr2(1) + rr2(2)
      d4 = rr2(3) - 0.5_dp*d3
      d5 = rr2(1) - rr2(2)
      do l = 1,3
        e(l,1) = rr2(l)*rr2(l+1)
        e(l,2) = rr2(l) + rr2(l+1) - 4._dp*e(l,1)
        e(l,3) = rr2(l+2) + e(l,1)
        e(l,1) = 3._dp*e(l,1)
      end do
      e(4,1) = d5**2
      e(4,2) = d3 - e(4,1)
      e(4,3) = rr2(3) + 0.25_dp*e(4,1)
      e(4,1) = 0.75_dp*e(4,1)
      e(5,1) = d4**2
      e(5,2) = 3._dp*rr2(3)*d3
      e(5,3) = 0.75_dp*d3**2
      dd(1) = rr(1)*rr(3)
      dd(2) = rr(2)*rr(1)
      dd(3) = rr(3)*rr(2)
      do l = 1,2
        e(l+5,1) = 3._dp*rr2(l+1)*dd(l)
        e(l+5,2) = dd(l)*(1._dp-4._dp*rr2(l+1))
        e(l+5,3) = dd(l)*(rr2(l+1)-1._dp)
      end do
      e(8,1) = dd(1)*d5*1.5_dp
      e(8,2) = dd(1)*(1.0_dp-2.0_dp*d5)
      e(8,3) = dd(1)*(0.5_dp*d5-1.0_dp)
      e(9,1) = d5*0.5_dp*d4*r3
      e(9,2) = -d5*rr2(3)*r3
      e(9,3) = d5*0.25_dp*(1.0_dp+rr2(3))*r3
      e(10,1) = rr2(1)*dd(3)*3.0_dp
      e(10,2) = (0.25_dp-rr2(1))*dd(3)*4.0_dp
      e(10,3) = dd(3)*(rr2(1)-1.0_dp)
      e(11,1) = 1.5_dp*dd(3)*d5
      e(11,2) = -dd(3)*(1.0_dp+2.0_dp*d5)
      e(11,3) = dd(3)*(1.0_dp+0.5_dp*d5)
      e(13,3) = 0.5_dp*d5*dd(2)
      e(13,2) = -2.0_dp*dd(2)*d5
      e(13,1) = e(13,3)*3.0_dp
      e(12,1) = d4*dd(1)*r3
      e(14,1) = d4*dd(3)*r3
      e(15,1) = d4*dd(2)*r3
      e(15,2) = -2.0_dp*r3*dd(2)*rr2(3)
      e(15,3) = 0.5_dp*r3*(1.0_dp+rr2(3))*dd(2)
      e(14,2) = r3*dd(3)*(d3-rr2(3))
      e(14,3) = -r3*0.5_dp*dd(3)*d3
      e(12,2) = r3*dd(1)*(d3-rr2(3))
      e(12,3) = -r3*0.5_dp*dd(1)*d3
      !
      dm(1:15) = 0._dp
      do l = 1,15
        do m = 1,3
          dm(l) = dm(l)+e(l,m)*skpar(ind(2,2,m-1))
        end do
      end do
      !
      do ir = 1,5
        do is = 1,ir
          ii = ir-is
          k = 5*ii-(ii*(ii-1))/2+is
          mat(ir+4,is+4) = mat(ir+4,is+4) + dm(k)
          mat(is+4,ir+4) = mat(is+4,ir+4) + dm(k)
        end do
      end do
    end subroutine skdd
    !
  end subroutine turnsk

  subroutine polint(xa,ya,n,x,y,dy)
    integer, intent(in)   :: n
    real(dp), intent(in)  :: xa(n),ya(n),x
    real(dp), intent(out) :: y,dy
    !
    real(dp)              :: c(n),d(n)
    real(dp)              :: den,dif,dift,ho,hp,w
    integer               :: i,m,ns
    !
    ns=1

    dif=abs(x-xa(1))
    do i = 1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
        ns=i
        dif=dift
      endif
      c(i)=ya(i)
      d(i)=ya(i)
    end do
    !
    y=ya(ns)
    ns=ns-1
    do m = 1,n-1
      do i = 1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) stop 'POLINT'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
      end do
      if (2*ns.lt.n-m)then
        dy=c(ns+1)
      else
        dy=d(ns)
        ns=ns-1
      endif
      y=y+dy
    end do
!
    return
  end subroutine polint

! *****************************************************************************
  subroutine urep_egr(rv,r,erep,derep,&
                      n_urpoly,urep,spdim,s_cut,srep,spxr,scoeff,surr,dograd)

    real(dp), intent(in)              :: rv(3)    ! interatomic distance vector
    real(dp), intent(in)              :: r        ! interatomic distance
    real(dp), intent(inout)           :: erep 
    real(dp), intent(inout)           :: derep(3)
    integer, intent(in)               :: n_urpoly,spdim
    real(dp), intent(in)              :: urep(:),s_cut,srep(3),surr(2)
    real(dp), intent(in)              :: spxr(:,:),scoeff(:,:)
    logical, intent(in)               :: dograd

    integer     :: ic, nsp, jsp, isp
    real(dp)    :: rz,de_z 
    !
    derep=0._dp
    de_z = 0._dp
    if (n_urpoly > 0) then
      !
      ! polynomial part
      !
      rz = urep(1) - r
      if (rz <= rtiny) RETURN
      do ic = 2,n_urpoly
        erep = erep + urep(ic) * rz**(ic)
      end do
      if (dograd) then
        do ic = 2,n_urpoly
          de_z = de_z - ic*urep(ic) * rz**(ic-1)
        end do
      end if
    else if (spdim > 0) then
      !
      ! spline part
      !
      ! This part is kind of proprietary Paderborn code and I won't reverse-engeneer
      ! everything in detail. What is obvious is documented.
      !
      ! This part has 4 regions:
      ! a) very long range is screened
      ! b) short-range is extrapolated with e-functions
      ! ca) normal range is approximated with a spline
      ! cb) longer range is extrapolated with an higher degree spline
      !
      if (r > s_cut) RETURN  ! screening (condition a)
      !
      if (r < s_cut) then
        ! a) short range
        erep = erep + exp(-srep(1)*r + srep(2)) + srep(3)
        if (dograd) de_z = de_z -srep(1)*exp(-srep(1)*r + srep(2))
      else
        !
        ! condition c). First determine between which places the spline is located:
        !
        ispg: do isp = 1,spdim ! condition ca)
          if (r <  spxr(isp,1)) cycle ispg ! distance is smaller than this spline range
          if (r >= spxr(isp,2)) cycle ispg ! distance is larger than this spline range
          ! at this point we have found the correct spline interval
          rz = r - spxr(isp,1)
          if (isp /= spdim) then
            nsp = 3 ! condition ca
            do jsp = 0,nsp
              erep = erep + scoeff(isp,jsp+1)*rz**(jsp)
            end do
            if (dograd) then
              do jsp = 1,nsp
                de_z = de_z + jsp*scoeff(isp,jsp+1)*rz**(jsp-1)
              end do
            end if
          else
            nsp = 5 ! condition cb
            do jsp = 0,nsp
              if( jsp <= 3 ) then
                erep = erep + scoeff(isp,jsp+1)*rz**(jsp)
              else
                erep = erep + surr(jsp-3)*rz**(jsp)
              endif
            end do
            if (dograd) then
              do jsp = 1,nsp
                if( jsp <= 3 ) then
                  de_z = de_z + jsp*scoeff(isp,jsp+1)*rz**(jsp-1)
                else
                  de_z = de_z + jsp*surr(jsp-3)*rz**(jsp-1)
                endif
              end do
            end if
          end if
          exit ispg
        end do ispg
      end if
    end if
    !
    if (dograd) then
       if ( r > 0.e-12_dp ) derep(1:3) = (de_z/r)*rv(1:3)
    end if

  end subroutine urep_egr

! gamma************************************************************************

  function gamma_rab_sr(r,ga,gb) result(gamma)
    !
    ! Purpose: Computes the short-range gamma parameter from exact Coulomb 
    !          interaction of normalized exp(-a*r) charge distribution - 1/r
    !
    !          Documentation refers to the Elstner et al, PRB 58 (1998) 7260 paper
    !
    real(dp), intent(in)    :: r,ga,gb ! distance and onsite hardness values for elements a and b
    real(dp) :: gamma
    real(dp) :: a,b         ! local work parameters based on ga and gb
    real(dp) :: g_sum,fac
    !
    gamma = 0.0_dp
    a = 3.2_dp*ga ! 3.2 = 16/5 in Eq. 18 and ff.
    b = 3.2_dp*gb
    g_sum = a + b
    if (g_sum.lt.tol_gamma) return ! hardness screening
    if (r < rtiny) then     ! This is for short distances but non-onsite terms
      ! This gives also correct diagonal elements (a=b, r=0)
      gamma = 0.5_dp*(a*b/g_sum + (a*b)**2/g_sum**3) 
      return
    end if
    !
    ! distinguish two cases: Gamma's are very close, e.g. for the same atom type, 
    !                        and Gamma's are different
    !
    if (abs(a-b) < rtiny) then
      fac = 1.6_dp*r*a*b/g_sum*(1.0_dp + a*b/g_sum**2)
      gamma = -(48.0_dp + 33._dp*fac + (9.0_dp + fac)*fac**2)*exp(-fac)/(48._dp*r)
    else
      gamma = -exp(-a*r)*(0.5_dp*a*b**4/(a**2-b**2)**2 - &
                (b**6 - 3._dp*a**2*b**4)/(r*(a**2-b**2)**3)) - & ! a-> b
               exp(-b*r)*(0.5_dp*b*a**4/(b**2-a**2)**2 - &
                (a**6 - 3._dp*b**2*a**4)/(r*(b**2-a**2)**3)) ! b-> a
    end if
    !
  end function gamma_rab_sr

! *****************************************************************************

END MODULE qs_dftb_matrices

! *****************************************************************************

