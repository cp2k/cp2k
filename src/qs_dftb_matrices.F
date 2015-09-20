!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Overlap and Hamiltonian matrices in DFTB
!> \author JGH
! *****************************************************************************
MODULE qs_dftb_matrices
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             is_hydrogen
  USE atprop_types,                    ONLY: atprop_array_init,&
                                             atprop_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type,&
                                             dftb_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       convert_offsets_to_sizes, cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_finalize, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, dbcsr_distribution_obj, &
       dbcsr_type_antisymmetric, dbcsr_type_symmetric
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_create,&
                                             ewald_env_get,&
                                             ewald_env_release,&
                                             ewald_env_set,&
                                             ewald_environment_type,&
                                             read_ewald_section
  USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                             ewald_pw_release,&
                                             ewald_pw_type
  USE input_constants,                 ONLY: do_fist_pol_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_charges
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             deallocate_particle_set,&
                                             particle_type
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE qmmm_types_low,                  ONLY: qmmm_env_qm_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qmmm_util,                       ONLY: spherical_cutoff_factor
  USE qs_dftb_coulomb,                 ONLY: build_dftb_coulomb
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                             qs_dftb_pairpot_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type,&
                                             set_ks_env
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE spme,                            ONLY: spme_forces,&
                                             spme_potential
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

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
  REAL(dp), PARAMETER                    :: rtiny = 1.e-10_dp
  ! eta(0) for mm atoms and non-scc qm atoms
  REAL(dp), PARAMETER                    :: eta_mm = 0.47_dp
  ! step size for qmmm finite difference
  REAL(dp), PARAMETER                    :: ddrmm = 0.0001_dp

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_matrices'

  PUBLIC :: build_dftb_matrices, build_dftb_ks_matrix, build_dftb_overlap
  PUBLIC :: build_dftb_qmmm_matrix, build_dftb_qmmm_matrix_zero, &
            build_dftb_qmmm_matrix_pc, deriv_dftb_qmmm_matrix, &
            deriv_dftb_qmmm_matrix_pc

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param calculate_forces ...
! *****************************************************************************
  SUBROUTINE build_dftb_matrices(qs_env,para_env,calculate_forces)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dftb_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER :: after, atom_a, atom_b, handle, i, iatom, ic, icol, ikind, img, &
      inode, irow, istat, iw, jatom, jkind, l1, l2, la, lb, llm, lmaxi, &
      lmaxj, m, n1, n2, n_urpoly, natom, natorb_a, natorb_b, nderivatives, &
      ngrd, ngrdcut, nimg, nkind, nmat, spdim
    INTEGER, DIMENSION(3)                    :: cell
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: defined, found, hb_sr_damp, &
                                                use_virial
    REAL(KIND=dp)                            :: ddr, dgam, dgrd, dr, drm, &
                                                drp, erep, erepij, f0, f1, &
                                                foab, fow, ga, gb, hb_para, &
                                                s_cut, urep_cut
    REAL(KIND=dp), DIMENSION(0:3)            :: eta_a, eta_b, skself
    REAL(KIND=dp), DIMENSION(10)             :: urep
    REAL(KIND=dp), DIMENSION(2)              :: surr
    REAL(KIND=dp), DIMENSION(3)              :: drij, force_ab, force_rr, &
                                                force_w, rij, srep
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: dfblock, dsblock, fblock, &
      fmatij, fmatji, gblock, pblock, sblock, scoeff, smatij, smatji, spxr, &
      wblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(block_p_type), DIMENSION(2:4)       :: dgblocks, dsblocks
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: gamma_matrix
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_h, matrix_p, matrix_s, &
                                                matrix_w
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind_a, dftb_kind_b
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), POINTER               :: dftb_potential
    TYPE(qs_dftb_pairpot_type), POINTER      :: dftb_param_ij, dftb_param_ji
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

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

    NULLIFY(logger, virial, atprop)
    logger => cp_get_default_logger()

    NULLIFY (matrix_h,matrix_s,matrix_p,matrix_w,gamma_matrix,atomic_kind_set,&
             qs_kind_set,sab_orb,ks_env)

    CALL get_qs_env(qs_env=qs_env,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    matrix_h_kp=matrix_h,&
                    matrix_s_kp=matrix_s,&
                    atprop=atprop,&
                    dft_control=dft_control,&
                    ks_env=ks_env)
    
    dftb_control => dft_control%qs_control%dftb_control
    nimg = dft_control%nimages
    ! Allocate the overlap and Hamiltonian matrix
    CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb)
    nderivatives=0
    IF(dftb_control%self_consistent .AND. calculate_forces)nderivatives=1
    CALL setup_matrices2(qs_env,nderivatives,nimg,matrix_s,"OVERLAP",sab_orb)
    CALL setup_matrices2(qs_env,0,nimg,matrix_h,"CORE HAMILTONIAN",sab_orb)
    CALL set_ks_env(ks_env,matrix_s_kp=matrix_s)
    CALL set_ks_env(ks_env,matrix_h_kp=matrix_h)

    NULLIFY (dftb_potential)
    CALL get_qs_env(qs_env=qs_env,dftb_potential=dftb_potential)
    NULLIFY (particle_set)
    CALL get_qs_env(qs_env=qs_env,particle_set=particle_set)

    ! gamma matrix allocation
    IF ( dftb_control%self_consistent ) THEN
       IF(calculate_forces) THEN
          nmat=4
       ELSE
          nmat=1
       END IF
       CALL get_qs_env(qs_env=qs_env,&
                       gamma_matrix=gamma_matrix)
       CALL setup_gamma(qs_env,nmat,gamma_matrix,sab_orb)
    END IF

    IF(calculate_forces) THEN
       NULLIFY (rho,force,matrix_w)
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       matrix_w_kp=matrix_w,&
                       virial=virial,&
                       force=force)
       CALL qs_rho_get(rho,rho_ao_kp=matrix_p)

       IF (SIZE(matrix_p,1) == 2) THEN
          DO img=1,nimg
             CALL cp_dbcsr_add(matrix_p(1,img)%matrix,matrix_p(2,img)%matrix,&
                  alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
             CALL cp_dbcsr_add(matrix_w(1,img)%matrix,matrix_w(2,img)%matrix,&
                  alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
          END DO
       END IF
       natom = SIZE(particle_set)
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       CPASSERT(istat==0)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind)
       use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    END IF
    ! atomic energy decomposition
    IF (atprop%energy) THEN
      natom = SIZE(particle_set)
      CALL atprop_array_init(atprop%atecc,natom)
    END IF

    NULLIFY(cell_to_index)
    IF (nimg>1) THEN
       CALL get_ks_env(ks_env=ks_env,kpoints=kpoints)
       CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index)
    END IF

    erep = 0._dp

    nkind = SIZE(atomic_kind_set)

    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
            iatom=iatom,jatom=jatom,inode=inode,r=rij,cell=cell)
       CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
       CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind_a)
       CALL get_dftb_atom_param(dftb_kind_a,&
              defined=defined,lmax=lmaxi,skself=skself,&
              eta=eta_a,natorb=natorb_a)

       IF (.NOT.defined .OR. natorb_a < 1) CYCLE
       CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind_b)
       CALL get_dftb_atom_param(dftb_kind_b,&
              defined=defined,lmax=lmaxj,eta=eta_b,natorb=natorb_b)

       IF (.NOT.defined .OR. natorb_b < 1) CYCLE

       ! retrieve information on F and S matrix
       dftb_param_ij => dftb_potential(ikind,jkind)
       dftb_param_ji => dftb_potential(jkind,ikind)
       ! assume table size and type is symmetric
       ngrd = dftb_param_ij%ngrd
       ngrdcut = dftb_param_ij%ngrdcut
       dgrd = dftb_param_ij%dgrd
       ddr = dgrd*0.1_dp
       CPASSERT(dftb_param_ij%llm==dftb_param_ji%llm)
       llm = dftb_param_ij%llm
       fmatij => dftb_param_ij%fmat
       smatij => dftb_param_ij%smat
       fmatji => dftb_param_ji%fmat
       smatji => dftb_param_ji%smat
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

       dr = SQRT(SUM(rij(:)**2))
       IF (NINT(dr/dgrd) <= ngrdcut) THEN

         IF(nimg==1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1),cell(2),cell(3))
            CPASSERT(ic > 0)
         END IF

         icol = MAX(iatom,jatom)
         irow = MIN(iatom,jatom)
         NULLIFY(sblock,fblock,gblock)
         CALL cp_dbcsr_get_block_p(matrix=matrix_s(1,ic)%matrix,&
                row=irow,col=icol,BLOCK=sblock,found=found)
         CPASSERT(found)
         CALL cp_dbcsr_get_block_p(matrix=matrix_h(1,ic)%matrix,&
                row=irow,col=icol,BLOCK=fblock,found=found)
         CPASSERT(found)
         IF ( dftb_control%self_consistent ) THEN
               CALL cp_dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix,&
                   row=irow,col=icol,BLOCK=gblock,found=found)
            CPASSERT(found)
         END IF

   
         IF (calculate_forces) THEN
            NULLIFY (pblock)
            CALL cp_dbcsr_get_block_p(matrix=matrix_p(1,ic)%matrix,&
                 row=irow,col=icol,block=pblock,found=found)
            CPASSERT(ASSOCIATED(pblock))
            NULLIFY (wblock)
            CALL cp_dbcsr_get_block_p(matrix=matrix_w(1,ic)%matrix,&
                 row=irow,col=icol,block=wblock,found=found)
            CPASSERT(ASSOCIATED(wblock))
            IF ( dftb_control%self_consistent ) THEN
               DO i=2,4
                  NULLIFY(dsblocks(i)%block)
                  CALL cp_dbcsr_get_block_p(matrix=matrix_s(i,ic)%matrix,&
                         row=irow,col=icol,BLOCK=dsblocks(i)%block,found=found)
                  CPASSERT(found)
                  NULLIFY (dgblocks(i)%block)
                  CALL cp_dbcsr_get_block_p(matrix=gamma_matrix(i)%matrix,&
                         row=irow,col=icol,BLOCK=dgblocks(i)%block,found=found)
                  CPASSERT(found)
               END DO
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
           CALL compute_block_sk(sblock,smatij,smatji,rij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)
           CALL compute_block_sk(fblock,fmatij,fmatji,rij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)
           IF(calculate_forces) THEN
             force_ab = 0._dp
             force_w  = 0._dp
             n1 = SIZE(fblock,1)
             n2 = SIZE(fblock,2)
             ! make sure that displacement is in the correct direction depending on the position
             ! of the block (upper or lower triangle)
             f0=1.0_dp
             IF ( irow == iatom ) f0=-1.0_dp

             ALLOCATE (dfblock(n1,n2),dsblock(n1,n2),STAT=istat)
             CPASSERT(istat==0)

             DO i=1,3 
               drij = rij 
               dfblock=0._dp; dsblock=0._dp

               drij(i) = rij(i) - ddr * f0
               CALL compute_block_sk(dsblock,smatij,smatji,drij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)
               CALL compute_block_sk(dfblock,fmatij,fmatji,drij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)             
 
               dsblock=-dsblock 
               dfblock=-dfblock 
 
               drij(i) = rij(i) + ddr * f0
               CALL compute_block_sk(dsblock,smatij,smatji,drij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)
               CALL compute_block_sk(dfblock,fmatij,fmatji,drij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)

               dfblock = dfblock/ (2.0_dp*ddr)
               dsblock = dsblock/ (2.0_dp*ddr)

               foab = 2.0_dp * SUM(dfblock*pblock)
               fow = -2.0_dp * SUM(dsblock*wblock)

               force_ab(i) = force_ab(i) + foab
               force_w(i) = force_w(i) + fow
               IF ( dftb_control%self_consistent ) THEN
                  CPASSERT(ASSOCIATED(dsblocks(i+1)%block))
                  dsblocks(i+1)%block = dsblocks(i+1)%block + dsblock
               END IF
             ENDDO
             IF ( use_virial ) THEN
               CALL virial_pair_force ( virial%pv_virial, -f0, force_ab, rij)
               CALL virial_pair_force ( virial%pv_virial, -f0, force_w, rij)
               IF (atprop%stress) THEN
                  f1 = 0.5_dp*f0
                  CALL virial_pair_force (atprop%atstress(:,:,iatom),-f1,force_ab,rij)
                  CALL virial_pair_force (atprop%atstress(:,:,iatom),-f1,force_w,rij)
                  CALL virial_pair_force (atprop%atstress(:,:,jatom),-f1,force_ab,rij)
                  CALL virial_pair_force (atprop%atstress(:,:,jatom),-f1,force_w,rij)
               END IF
             END IF
             DEALLOCATE (dfblock,dsblock,STAT=istat)
             CPASSERT(istat==0)
           END IF
         END IF

         IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            IF ( irow == iatom ) force_ab = -force_ab
            IF ( irow == iatom ) force_w = -force_w
            force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) - force_ab(:)
            force(jkind)%all_potential(:,atom_b) = force(jkind)%all_potential(:,atom_b) + force_ab(:)
            force(ikind)%overlap(:,atom_a) = force(ikind)%overlap(:,atom_a) - force_w(:)
            force(jkind)%overlap(:,atom_b) = force(jkind)%overlap(:,atom_b) + force_w(:)
         END IF

         ! gamma matrix
         IF ( dftb_control%self_consistent ) THEN
            hb_sr_damp = dftb_control%hb_sr_damp
            IF (hb_sr_damp) THEN
               ! short range correction enabled only when iatom XOR jatom are hydrogens
               hb_sr_damp = is_hydrogen(particle_set(iatom)%atomic_kind).NEQV.&
                            is_hydrogen(particle_set(jatom)%atomic_kind)
            END IF
            IF (hb_sr_damp) THEN
               hb_para = dftb_control%hb_sr_para
            ELSE
               hb_para = 0.0_dp
            END IF
            ga = eta_a(0)
            gb = eta_b(0)
            gblock(1,1)= gblock(1,1) + gamma_rab_sr(dr,ga,gb,hb_para)
            IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
               drp = dr + ddr
               drm = dr - ddr
               dgam = 0.5_dp*(gamma_rab_sr(drp,ga,gb,hb_para)-gamma_rab_sr(drm,ga,gb,hb_para))/ddr
               DO i=1,3
                  CPASSERT(ASSOCIATED(dgblocks(i+1)%block))
                  IF ( irow == iatom ) THEN
                     dgblocks(i+1)%block(1,1)= dgblocks(i+1)%block(1,1) + dgam*rij(i)/dr
                  ELSE
                     dgblocks(i+1)%block(1,1)= dgblocks(i+1)%block(1,1) - dgam*rij(i)/dr
                  END IF
               END DO
            END IF
         END IF

       END IF

       ! repulsive potential
       IF ((dr <= urep_cut .OR. spdim > 0) .AND. dr > 0.001_dp) THEN
          erepij = 0._dp
          CALL urep_egr(rij,dr,erepij,force_rr,&
              n_urpoly,urep,spdim,s_cut,srep,spxr,scoeff,surr,calculate_forces)
          erep = erep + erepij
          IF(atprop%energy) THEN
             atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*erepij
             atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*erepij
          END IF
          IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            force(ikind)%repulsive(:,atom_a) =&
                force(ikind)%repulsive(:,atom_a) - force_rr(:)
            force(jkind)%repulsive(:,atom_b) =&
                force(jkind)%repulsive(:,atom_b) + force_rr(:)
            IF ( use_virial ) THEN
              CALL virial_pair_force ( virial%pv_virial, -1._dp, force_rr, rij)
              IF(atprop%stress) THEN
                CALL virial_pair_force(atprop%atstress(:,:,iatom),-0.5_dp,force_rr,rij)
                CALL virial_pair_force(atprop%atstress(:,:,jatom),-0.5_dp,force_rr,rij)
              END IF
            END IF
          END IF
       END IF

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     IF ( dftb_control%self_consistent ) THEN
        DO i=1,SIZE(gamma_matrix)
           CALL cp_dbcsr_finalize(gamma_matrix(i)%matrix)
        ENDDO
        CALL set_ks_env(ks_env,gamma_matrix=gamma_matrix)
     ENDIF
     DO i=1,SIZE(matrix_s,1)
        DO img=1,nimg
           CALL cp_dbcsr_finalize(matrix_s(i,img)%matrix)
        END DO
     ENDDO
     DO i=1,SIZE(matrix_h,1)
        DO img=1,nimg
           CALL cp_dbcsr_finalize(matrix_h(i,img)%matrix)
        END DO
     ENDDO

     ! set repulsive energy
     CALL mp_sum(erep,para_env%group)
     energy%repulsive = erep

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN"),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log")
       CALL section_vals_val_get(qs_env%input,"DFT%PRINT%AO_MATRICES%NDIGITS",i_val=after)
       after = MIN(MAX(after,1),16)
       DO img=1,nimg
          CALL cp_dbcsr_write_sparse_matrix(matrix_h(1,img)%matrix,4,after,qs_env,para_env,&
               output_unit=iw)
       END DO

       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN")
     END IF

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP"),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",&
            extension=".Log")
       CALL section_vals_val_get(qs_env%input,"DFT%PRINT%AO_MATRICES%NDIGITS",i_val=after)
       after = MIN(MAX(after,1),16)
       DO img=1,nimg
          CALL cp_dbcsr_write_sparse_matrix(matrix_s(1,img)%matrix,4,after,qs_env,para_env,&
             output_unit=iw)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
                 qs_env%input,"DFT%PRINT%AO_MATRICES/DERIVATIVES"),cp_p_file)) THEN
            DO i=2,SIZE(matrix_s,1)
              CALL cp_dbcsr_write_sparse_matrix(matrix_s(i,img)%matrix,4,after,qs_env,para_env,&
                   output_unit=iw)
            END DO
          END IF
       END DO

       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/OVERLAP")
     END IF

     IF (calculate_forces) THEN
       IF (SIZE(matrix_p,1) == 2) THEN
          DO img=1,nimg
             CALL cp_dbcsr_add(matrix_p(1,img)%matrix,matrix_p(2,img)%matrix,alpha_scalar=1.0_dp,&
                  beta_scalar=-1.0_dp)
          END DO
       END IF
       DEALLOCATE(atom_of_kind,STAT=istat)
       CPASSERT(istat==0)
     END IF

    CALL timestop(handle)

  END SUBROUTINE build_dftb_matrices

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
! *****************************************************************************
  SUBROUTINE build_dftb_ks_matrix(qs_env,calculate_forces,just_energy)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ikind, &
                                                img, ispin, istat, natom, &
                                                nkind, nspins, output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pc_ener, qmmm_el, zeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mcharge, occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p1, mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: ks_matrix, matrix_h, &
                                                matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: scf_section

    CALL timeset(routineN,handle)
    NULLIFY(dft_control, logger, scf_section,matrix_p, particle_set, ks_env,&
            ks_matrix, rho, energy)
    logger => cp_get_default_logger()
    failure=.FALSE.
    CPASSERT(ASSOCIATED(qs_env))

    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    matrix_h_kp=matrix_h,&
                    para_env=para_env,&
                    ks_env=ks_env,&
                    matrix_ks_kp=ks_matrix,&
                    rho=rho,&
                    energy=energy)

    energy%hartree = 0.0_dp
    energy%qmmm_el = 0.0_dp

    scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF")
    nspins=dft_control%nspins
    CPASSERT(ASSOCIATED(matrix_h))
    CPASSERT(ASSOCIATED(rho))
    CPASSERT(SIZE(ks_matrix)>0)

    DO ispin=1,nspins
       DO img=1,SIZE(ks_matrix,2)
          ! copy the core matrix into the fock matrix
          CALL cp_dbcsr_copy(ks_matrix(ispin,img)%matrix,matrix_h(1,img)%matrix)
       END DO
    END DO

    IF ( dft_control%qs_control%dftb_control%self_consistent ) THEN
      ! Mulliken charges
      CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,&
           matrix_s_kp=matrix_s)
      CALL qs_rho_get(rho,rho_ao_kp=matrix_p)
      natom=SIZE(particle_set)
      ALLOCATE(charges(natom,nspins),STAT=istat)
      CPASSERT(istat==0)
      !
      CALL mulliken_charges(matrix_p,matrix_s,para_env,charges)
      !
      ALLOCATE(mcharge(natom),STAT=istat)
      CPASSERT(istat==0)
      nkind = SIZE(atomic_kind_set)
      DO ikind=1,nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
         CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
         CALL get_dftb_atom_param(dftb_kind,zeff=zeff)
         DO iatom=1,natom
           atom_a = atomic_kind_set(ikind)%atom_list(iatom)
           mcharge(atom_a) = zeff - SUM(charges(atom_a,1:nspins))
         END DO
      END DO
      DEALLOCATE(charges,STAT=istat)
      CPASSERT(istat==0)

      CALL build_dftb_coulomb(qs_env,ks_matrix,rho,mcharge,energy,&
            calculate_forces,just_energy)

      DEALLOCATE(mcharge,STAT=istat)
      CPASSERT(istat==0)

    END IF

    IF (qs_env%qmmm) THEN
       CPASSERT(SIZE(ks_matrix,2)==1)
       DO ispin = 1, nspins
          ! If QM/MM sumup the 1el Hamiltonian
          CALL cp_dbcsr_add(ks_matrix(ispin,1)%matrix,qs_env%ks_qmmm_env%matrix_h(1)%matrix,&
               1.0_dp,1.0_dp)
          CALL qs_rho_get(rho,rho_ao=matrix_p1)
          ! Compute QM/MM Energy
          CALL cp_dbcsr_trace(qs_env%ks_qmmm_env%matrix_h(1)%matrix,&
               matrix_p1(ispin)%matrix,trace=qmmm_el)
          energy%qmmm_el = energy%qmmm_el + qmmm_el
       END DO
       pc_ener = qs_env%ks_qmmm_env%pc_ener
       energy%qmmm_el = energy%qmmm_el + pc_ener
    END IF

    energy%total = energy%core + energy%hartree + energy%qmmm_el + &
                   energy%repulsive + energy%dispersion + energy%dftb3

    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
       extension=".scfLog")
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T9,A,T60,F20.10))")&
            "Repulsive pair potential energy:               ",energy%repulsive,&
            "Zeroth order Hamiltonian energy:               ",energy%core,&
            "Charge fluctuation energy:                     ",energy%hartree,&
            "London dispersion energy:                      ",energy%dispersion
       IF ( dft_control%qs_control%dftb_control%dftb3_diagonal ) THEN
          WRITE (UNIT=output_unit,FMT="(T9,A,T60,F20.10)")&
               "DFTB3 3rd Order Energy Correction              ",energy%dftb3
       END IF
       IF (qs_env%qmmm) THEN
          WRITE (UNIT=output_unit,FMT="(T9,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%DETAILED_ENERGY")
    ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
    IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
       CPASSERT(SIZE(ks_matrix,2)==1)
       CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array)
       DO ispin=1,SIZE(mo_derivs)
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff_b=mo_coeff, occupation_numbers=occupation_numbers )
          IF(.NOT.mo_array(ispin)%mo_set%use_mo_coeff_b) THEN
             CPASSERT(.FALSE.)
          ENDIF
          CALL cp_dbcsr_multiply('n','n',1.0_dp,ks_matrix(ispin,1)%matrix,mo_coeff,&
               0.0_dp,mo_derivs(ispin)%matrix)
       ENDDO
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE build_dftb_ks_matrix

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nderivative ...
!> \param matrix_s ...
! *****************************************************************************
  SUBROUTINE build_dftb_overlap(qs_env,nderivative,matrix_s)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dftb_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, indder, inode, irow, istat, j, &
      jatom, jkind, l1, l2, la, lb, llm, lmaxi, lmaxj, m, n1, n2, natom, &
      natorb_a, natorb_b, ngrd, ngrdcut, nkind
    LOGICAL                                  :: defined, found
    REAL(KIND=dp)                            :: ddr, dgrd, dr, f0
    REAL(KIND=dp), DIMENSION(0:3)            :: skself
    REAL(KIND=dp), DIMENSION(3)              :: drij, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsblock, dsblockm, sblock, &
                                                smatij, smatji
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dsblock1
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(block_p_type), DIMENSION(2:10)      :: dsblocks
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind_a, dftb_kind_b
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), POINTER               :: dftb_potential
    TYPE(qs_dftb_pairpot_type), POINTER      :: dftb_param_ij, dftb_param_ji
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)

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
    logger => cp_get_default_logger()

    NULLIFY ( atomic_kind_set, qs_kind_set, sab_orb)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
                    dft_control=dft_control)

    dftb_control => dft_control%qs_control%dftb_control

    NULLIFY (dftb_potential)
    CALL get_qs_env(qs_env=qs_env,&
                    dftb_potential=dftb_potential)

    nkind = SIZE(atomic_kind_set)

    ! Allocate the overlap matrix
    CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb)
    CALL setup_matrices1(qs_env,nderivative,matrix_s,'OVERLAP',sab_orb)

    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
            iatom=iatom,jatom=jatom,inode=inode,r=rij)

       CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
       CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind_a)
       CALL get_dftb_atom_param(dftb_kind_a,&
              defined=defined,lmax=lmaxi,skself=skself,&
              natorb=natorb_a)

       IF (.NOT.defined .OR. natorb_a < 1) CYCLE

       CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind_b)
       CALL get_dftb_atom_param(dftb_kind_b,&
              defined=defined,lmax=lmaxj,natorb=natorb_b)

       IF (.NOT.defined .OR. natorb_b < 1) CYCLE

       ! retrieve information on F and S matrix
       dftb_param_ij => dftb_potential(ikind,jkind)
       dftb_param_ji => dftb_potential(jkind,ikind)
       ! assume table size and type is symmetric
       ngrd = dftb_param_ij%ngrd
       ngrdcut = dftb_param_ij%ngrdcut
       dgrd = dftb_param_ij%dgrd
       ddr = dgrd*0.1_dp
       CPASSERT(dftb_param_ij%llm==dftb_param_ji%llm)
       llm = dftb_param_ij%llm
       smatij => dftb_param_ij%smat
       smatji => dftb_param_ji%smat

       dr = SQRT(SUM(rij(:)**2))
       IF (NINT(dr/dgrd) <= ngrdcut) THEN

          icol = MAX(iatom,jatom); irow = MIN(iatom,jatom)
           
          NULLIFY(sblock)
          CALL cp_dbcsr_get_block_p(matrix=matrix_s(1)%matrix,&
               row=irow,col=icol,BLOCK=sblock,found=found)
          CPASSERT(found)

          IF (nderivative.gt.0) THEN
             DO i=2,SIZE(matrix_s,1)
                NULLIFY(dsblocks(i)%block)
                CALL cp_dbcsr_get_block_p(matrix=matrix_s(i)%matrix,&
                     row=irow,col=icol,BLOCK=dsblocks(i)%block,found=found)
             END DO
          END IF

          IF (iatom == jatom .AND. dr < 0.001_dp) THEN
             ! diagonal block
             DO i=1,natorb_a
                sblock(i,i) = sblock(i,i) + 1._dp
             END DO
          ELSE
             ! off-diagonal block
             CALL compute_block_sk(sblock,smatij,smatji,rij,ngrd,ngrdcut,dgrd,&
                                   llm,lmaxi,lmaxj,irow,iatom)

             IF(nderivative.ge.1) THEN
                n1 = SIZE(sblock,1); n2 = SIZE(sblock,2)
                indder=1 ! used to put the 2nd derivatives in the correct matric (5=xx,8=yy,10=zz)

                ALLOCATE (dsblock1(n1,n2,3),dsblock(n1,n2),dsblockm(n1,n2),STAT=istat)
                CPASSERT(istat==0)
                dsblock1=0.0_dp
                DO i=1,3
                   dsblock=0._dp; dsblockm=0.0_dp
                   drij=rij
                   f0=1.0_dp; IF ( irow == iatom ) f0=-1.0_dp
 
                   drij(i) = rij(i) - ddr * f0
                   CALL compute_block_sk(dsblockm,smatij,smatji,drij,ngrd,ngrdcut,dgrd,&
                                         llm,lmaxi,lmaxj,irow,iatom)
 
                   drij(i) = rij(i) + ddr * f0
                   CALL compute_block_sk(dsblock,smatij,smatji,drij,ngrd,ngrdcut,dgrd,&
                                         llm,lmaxi,lmaxj,irow,iatom)
    
                   dsblock1(:,:,i)= (dsblock+dsblockm)
                   dsblock = dsblock-dsblockm
                   dsblock = dsblock/(2.0_dp*ddr)
    
                   CPASSERT(ASSOCIATED(dsblocks(i+1)%block))
                   dsblocks(i+1)%block = dsblocks(i+1)%block + dsblock
                   IF(nderivative.gt.1) THEN
                      indder=indder+5-i
                      dsblocks(indder)%block=0.0_dp
                      dsblocks(indder)%block=dsblocks(indder)%block+&
                              (dsblock1(:,:,i)-2.0_dp*sblock)/ddr**2
                   END IF
                ENDDO

                IF (nderivative.gt.1 ) THEN
                   DO i=1,2
                      DO j=i+1,3
                         dsblock=0._dp; dsblockm=0.0_dp
                         drij=rij
                         f0=1.0_dp; IF ( irow == iatom ) f0=-1.0_dp
   
                         drij(i) = rij(i) - ddr * f0; drij(j) = rij(j) - ddr * f0
                         CALL compute_block_sk(dsblockm,smatij,smatji,drij,ngrd,ngrdcut,dgrd,&
                                               llm,lmaxi,lmaxj,irow,iatom)
                      
                         drij(i) = rij(i) + ddr * f0; drij(j) = rij(j) + ddr * f0
                         CALL compute_block_sk(dsblock,smatij,smatji,drij,ngrd,ngrdcut,dgrd,&
                                               llm,lmaxi,lmaxj,irow,iatom)
                      
                         indder=2+2*i+j 
                         dsblocks(indder)%block=0.0_dp
                         dsblocks(indder)%block=dsblocks(indder)%block+(dsblock+dsblockm-&
                                  dsblock1(:,:,i)-dsblock1(:,:,j)+2.0_dp*sblock)/(2.0_dp*ddr**2)
                      END DO
                   END DO   
                END IF

                DEALLOCATE (dsblock1,dsblock,dsblockm,STAT=istat)
                CPASSERT(istat==0)
             END IF
          END IF
       END IF
     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     DO i=1,SIZE(matrix_s,1)
        CALL cp_dbcsr_finalize(matrix_s(i)%matrix)
     ENDDO

    CALL timestop(handle)

  END SUBROUTINE build_dftb_overlap

! *****************************************************************************
!> \brief Constructs the 1-el DFTB hamiltonian
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \author JGH 10.2014 [created]
! *****************************************************************************
  SUBROUTINE build_dftb_qmmm_matrix(qs_env,qmmm_env,particles_mm,mm_cell,para_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_qmmm_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iatom, ikind, &
                                                istat, jatom, natom, natorb, &
                                                nkind
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: defined, failure, found
    REAL(KIND=dp)                            :: pc_ener, zeff
    REAL(KIND=dp), DIMENSION(0:3)            :: eta_a
    REAL(KIND=dp), DIMENSION(:), POINTER     :: qpot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particles_qm,&
                    qs_kind_set=qs_kind_set,&
                    rho=rho,&
                    natom=natom)
    dftb_control => dft_control%qs_control%dftb_control

    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=qs_env%input)
    NULLIFY(matrix_s)
    CALL build_dftb_overlap(qs_env,0,matrix_s)

    ALLOCATE(qpot(natom),STAT=istat)
    CPASSERT(istat==0)
    qpot = 0.0_dp
    pc_ener = 0.0_dp

    nkind = SIZE(atomic_kind_set)
    DO ikind=1,nkind
       CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list)
       NULLIFY(dftb_kind)
       CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
       CALL get_dftb_atom_param(dftb_kind,zeff=zeff,&
            defined=defined,eta=eta_a,natorb=natorb)
       ! use mm charge smearing for non-scc cases
       IF(.NOT.dftb_control%self_consistent) eta_a(0) = eta_mm
       IF (.NOT.defined .OR. natorb < 1) CYCLE
       DO i = 1, SIZE(list)
          iatom = list(i)
          CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%Potentials,particles_mm,&
               qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,&
               qmmm_env%spherical_cutoff,particles_qm)
          ! Possibly added charges
          IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
             CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%added_charges%potentials,&
                  qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                  qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,qmmm_env%spherical_cutoff,&
                  particles_qm)
          END IF
          pc_ener = pc_ener + qpot(iatom)*zeff
       END DO
    END DO

    ! Allocate the core Hamiltonian matrix
    CALL get_qs_env(qs_env=qs_env,ks_qmmm_env=ks_qmmm_env_loc)
    matrix_h => ks_qmmm_env_loc%matrix_h
    CALL cp_dbcsr_allocate_matrix_set(matrix_h,1)
    ALLOCATE(matrix_h(1)%matrix,STAT=istat)
    CPASSERT(istat==0)
    CALL cp_dbcsr_init(matrix_h(1)%matrix)
    CALL cp_dbcsr_copy(matrix_h(1)%matrix,matrix_s(1)%matrix,&
                       name="QMMM HAMILTONIAN MATRIX")
    CALL cp_dbcsr_set(matrix_h(1)%matrix,0.0_dp)

    CALL cp_dbcsr_iterator_start(iter, matrix_s(1)%matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
       NULLIFY(hblock)
       CALL cp_dbcsr_get_block_p(matrix=matrix_h(1)%matrix,&
            row=iatom,col=jatom,block=hblock,found=found)
       CPASSERT(found)
       hblock = hblock - 0.5_dp*sblock*(qpot(iatom)+qpot(jatom))
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    ks_qmmm_env_loc%matrix_h => matrix_h
    ks_qmmm_env_loc%pc_ener = pc_ener

    DEALLOCATE(qpot,STAT=istat)
    CPASSERT(istat==0)

    CALL cp_dbcsr_deallocate_matrix_set ( matrix_s)

    CALL timestop(handle)

  END SUBROUTINE build_dftb_qmmm_matrix

! *****************************************************************************
!> \brief Constructs an empty 1-el DFTB hamiltonian
!> \param qs_env ...
!> \param para_env ...
!> \author JGH 10.2014 [created]
! *****************************************************************************
  SUBROUTINE build_dftb_qmmm_matrix_zero(qs_env,para_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_qmmm_matrix_zero', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_s
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=qs_env%input)
    NULLIFY(matrix_s)
    CALL build_dftb_overlap(qs_env,0,matrix_s)

    ! Allocate the core Hamiltonian matrix
    CALL get_qs_env(qs_env=qs_env,ks_qmmm_env=ks_qmmm_env_loc)
    matrix_h => ks_qmmm_env_loc%matrix_h
    CALL cp_dbcsr_allocate_matrix_set(matrix_h,1)
    ALLOCATE(matrix_h(1)%matrix,STAT=istat)
    CPASSERT(istat==0)
    CALL cp_dbcsr_init(matrix_h(1)%matrix)
    CALL cp_dbcsr_copy(matrix_h(1)%matrix,matrix_s(1)%matrix,&
                       name="QMMM HAMILTONIAN MATRIX")
    CALL cp_dbcsr_set(matrix_h(1)%matrix,0.0_dp)
    ks_qmmm_env_loc%matrix_h => matrix_h
    ks_qmmm_env_loc%pc_ener = 0.0_dp

    CALL cp_dbcsr_deallocate_matrix_set ( matrix_s)

    CALL timestop(handle)

  END SUBROUTINE build_dftb_qmmm_matrix_zero

! *****************************************************************************
!> \brief Constructs the 1-el DFTB hamiltonian
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \author JGH 10.2014 [created]
! *****************************************************************************
  SUBROUTINE build_dftb_qmmm_matrix_pc(qs_env,qmmm_env,particles_mm,mm_cell,para_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_qmmm_matrix_pc', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, do_ipol, ewald_type, handle, i, iatom, ikind, imm, imp, &
      indmm, ipot, istat, jatom, natom, natorb, nkind, nmm
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: defined, do_multipoles, &
                                                failure, found
    REAL(KIND=dp)                            :: alpha, pc_ener, zeff
    REAL(KIND=dp), DIMENSION(0:3)            :: eta_a
    REAL(KIND=dp), DIMENSION(2)              :: rcutoff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges_mm, qpot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: atoms_mm, particles_qm
    TYPE(qmmm_pot_type), POINTER             :: Pot
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: ewald_section, &
                                                poisson_section, print_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particles_qm,&
                    qs_kind_set=qs_kind_set,&
                    rho=rho,&
                    natom=natom)
    dftb_control => dft_control%qs_control%dftb_control

    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=qs_env%input)
    NULLIFY(matrix_s)
    CALL build_dftb_overlap(qs_env,0,matrix_s)

    ALLOCATE(qpot(natom),STAT=istat)
    CPASSERT(istat==0)
    qpot = 0.0_dp
    pc_ener = 0.0_dp

    ! Create Ewald environments
    poisson_section => section_vals_get_subs_vals(qs_env%input,"MM%POISSON")
    CALL ewald_env_create(ewald_env,para_env)
    CALL ewald_env_set(ewald_env,poisson_section=poisson_section)
    ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD")
    CALL read_ewald_section(ewald_env,ewald_section)
    print_section => section_vals_get_subs_vals(qs_env%input,"PRINT%GRID_INFORMATION")
    CALL ewald_pw_create(ewald_pw,ewald_env,mm_cell,mm_cell,print_section=print_section)

    CALL ewald_env_get(ewald_env,ewald_type=ewald_type,do_multipoles=do_multipoles,do_ipol=do_ipol)
    IF(do_multipoles) CALL cp_unimplemented_error(fromWhere=routineP, &
       message="No multipole force fields allowed in DFTB QM/MM",error_level=cp_fatal_level)
    IF(do_ipol /= do_fist_pol_none) CALL cp_unimplemented_error(fromWhere=routineP, &
       message="No polarizable force fields allowed in DFTB QM/MM",error_level=cp_fatal_level)

    SELECT CASE(ewald_type)
       CASE(do_ewald_pme)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="PME Ewald type not implemented for DFTB/QMMM",error_level=cp_fatal_level)
       CASE(do_ewald_ewald,do_ewald_spme)
          DO ipot = 1,SIZE(qmmm_env%Potentials)
             Pot => qmmm_env%Potentials(ipot)%Pot
             nmm = SIZE(Pot%mm_atom_index)
             ! get a 'clean' mm particle set
             NULLIFY(atoms_mm)
             CALL allocate_particle_set(atoms_mm,nmm)
             ALLOCATE(charges_mm(nmm),STAT=istat)
             CPASSERT(istat==0)
             DO Imp=1,nmm
                Imm = Pot%mm_atom_index(Imp)
                IndMM = qmmm_env%mm_atom_index(Imm)
                atoms_mm(Imp)%r = particles_mm(IndMM)%r
                atoms_mm(Imp)%atomic_kind => particles_mm(IndMM)%atomic_kind
                charges_mm(Imp) = qmmm_env%mm_atom_chrg(Imm)
             END DO
             IF(ewald_type == do_ewald_ewald) THEN
                CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="Ewald not implemented for DFTB/QMMM",error_level=cp_fatal_level)
             ELSE IF(ewald_type == do_ewald_spme) THEN
                ! spme electrostatic potential
                CALL spme_potential(ewald_env,ewald_pw,mm_cell,atoms_mm,charges_mm,particles_qm,qpot)
             END IF
             CALL deallocate_particle_set(atoms_mm)
             DEALLOCATE(charges_mm,STAT=istat)
             CPASSERT(istat==0)
          END DO
          IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
             DO ipot = 1,SIZE(qmmm_env%added_charges%Potentials)
                Pot => qmmm_env%added_charges%Potentials(ipot)%Pot
                nmm = SIZE(Pot%mm_atom_index)
                ! get a 'clean' mm particle set
                NULLIFY(atoms_mm)
                CALL allocate_particle_set(atoms_mm,nmm)
                ALLOCATE(charges_mm(nmm),STAT=istat)
                CPASSERT(istat==0)
                DO Imp=1,nmm
                   Imm = Pot%mm_atom_index(Imp)
                   IndMM = qmmm_env%added_charges%mm_atom_index(Imm)
                   atoms_mm(Imp)%r = qmmm_env%added_charges%added_particles(IndMM)%r
                   atoms_mm(Imp)%atomic_kind => qmmm_env%added_charges%added_particles(IndMM)%atomic_kind
                   charges_mm(Imp) = qmmm_env%added_charges%mm_atom_chrg(Imm)
                END DO
                IF(ewald_type == do_ewald_ewald) THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Ewald not implemented for DFTB/QMMM",error_level=cp_fatal_level)
                ELSE IF(ewald_type == do_ewald_spme) THEN
                   ! spme electrostatic potential
                   CALL spme_potential(ewald_env,ewald_pw,mm_cell,atoms_mm,charges_mm,particles_qm,qpot)
                END IF
                CALL deallocate_particle_set(atoms_mm)
                DEALLOCATE(charges_mm,STAT=istat)
                CPASSERT(istat==0)
             END DO
          END IF
          CALL mp_sum(qpot,para_env%group)
          ! Add Ewald and DFTB short range corrections
          ! This is effectively using a minimum image convention!
          ! Set rcutoff to values compatible with alpha Ewald
          CALL ewald_env_get(ewald_env,rcut=rcutoff(1),alpha=alpha)
          rcutoff(2) = 0.025_dp*rcutoff(1)
          rcutoff(1) = 2.0_dp*rcutoff(1)
          nkind = SIZE(atomic_kind_set)
          DO ikind=1,nkind
             CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list)
             NULLIFY(dftb_kind)
             CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
             CALL get_dftb_atom_param(dftb_kind,zeff=zeff,&
                  defined=defined,eta=eta_a,natorb=natorb)
             ! use mm charge smearing for non-scc cases
             IF(.NOT.dftb_control%self_consistent) eta_a(0) = eta_mm
             IF (.NOT.defined .OR. natorb < 1) CYCLE
             DO i = 1, SIZE(list)
                iatom = list(i)
                CALL build_mm_pot(qpot(iatom),1,eta_a(0),qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,rcutoff,&
                     particles_qm)
                CALL build_mm_pot(qpot(iatom),2,alpha,qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,rcutoff,&
                     particles_qm)
                ! Possibly added charges
                IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                   CALL build_mm_pot(qpot(iatom),1,eta_a(0),qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,rcutoff,particles_qm)
                   CALL build_mm_pot(qpot(iatom),2,alpha,qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,rcutoff,particles_qm)
                END IF
                pc_ener = pc_ener + qpot(iatom)*zeff
             END DO
          END DO
       CASE(do_ewald_none)
          ! Simply summing up charges with 1/R (DFTB corrected)
          nkind = SIZE(atomic_kind_set)
          DO ikind=1,nkind
             CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list)
             NULLIFY(dftb_kind)
             CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
             CALL get_dftb_atom_param(dftb_kind,zeff=zeff,&
                  defined=defined,eta=eta_a,natorb=natorb)
             ! use mm charge smearing for non-scc cases
             IF(.NOT.dftb_control%self_consistent) eta_a(0) = eta_mm
             IF (.NOT.defined .OR. natorb < 1) CYCLE
             DO i = 1, SIZE(list)
                iatom = list(i)
                CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,&
                     qmmm_env%spherical_cutoff,particles_qm)
                ! Possibly added charges
                IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                   CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,qmmm_env%spherical_cutoff,&
                        particles_qm)
                END IF
                pc_ener = pc_ener + qpot(iatom)*zeff
             END DO
          END DO
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Unknown Ewald type!",error_level=cp_fatal_level)
    END SELECT

    ! Allocate the core Hamiltonian matrix
    CALL get_qs_env(qs_env=qs_env,ks_qmmm_env=ks_qmmm_env_loc)
    matrix_h => ks_qmmm_env_loc%matrix_h 
    CALL cp_dbcsr_allocate_matrix_set(matrix_h,1)
    ALLOCATE(matrix_h(1)%matrix,STAT=istat)
    CPASSERT(istat==0)
    CALL cp_dbcsr_init(matrix_h(1)%matrix)
    CALL cp_dbcsr_copy(matrix_h(1)%matrix,matrix_s(1)%matrix,&
                       name="QMMM HAMILTONIAN MATRIX")
    CALL cp_dbcsr_set(matrix_h(1)%matrix,0.0_dp)

    CALL cp_dbcsr_iterator_start(iter, matrix_s(1)%matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
       NULLIFY(hblock)
       CALL cp_dbcsr_get_block_p(matrix=matrix_h(1)%matrix,&
            row=iatom,col=jatom,block=hblock,found=found)
       CPASSERT(found)
       hblock = hblock - 0.5_dp*sblock*(qpot(iatom)+qpot(jatom))
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    ks_qmmm_env_loc%matrix_h => matrix_h
    ks_qmmm_env_loc%pc_ener = pc_ener

    DEALLOCATE(qpot,STAT=istat)
    CPASSERT(istat==0)

    ! Release Ewald environment
    CALL ewald_env_release(ewald_env)
    CALL ewald_pw_release(ewald_pw)

    CALL cp_dbcsr_deallocate_matrix_set ( matrix_s)

    CALL timestop(handle)

  END SUBROUTINE build_dftb_qmmm_matrix_pc

! *****************************************************************************
!> \brief Constructs the derivative w.r.t. 1-el DFTB hamiltonian QMMM terms
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \param calc_force ...
!> \param Forces ...
!> \param Forces_added_charges ...
!> \author JGH 10.2014 [created]
! *****************************************************************************
  SUBROUTINE deriv_dftb_qmmm_matrix(qs_env,qmmm_env,particles_mm,mm_cell,para_env,&
       calc_force, Forces, Forces_added_charges)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'deriv_dftb_qmmm_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, blk, handle, i, iatom, ikind, iqm, istat, jatom, &
      natom, natorb, nkind, nspins, number_qm_atoms
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: defined, failure, found
    REAL(KIND=dp)                            :: fi, gmij, zeff
    REAL(KIND=dp), DIMENSION(0:3)            :: eta_a
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mcharge, qpot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, dsblock, Forces_QM, &
                                                pblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (calc_force) THEN
       NULLIFY (rho, atomic_kind_set, qs_kind_set, particles_qm)
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       atomic_kind_set=atomic_kind_set,&
                       qs_kind_set=qs_kind_set,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       dft_control=dft_control,&
                       particle_set=particles_qm,&
                       natom=number_qm_atoms)
       dftb_control => dft_control%qs_control%dftb_control

       NULLIFY(matrix_s)
       CALL build_dftb_overlap(qs_env,1,matrix_s)
       CALL qs_rho_get(rho, rho_ao=matrix_p)

       nspins = dft_control%nspins
       nkind  = SIZE(atomic_kind_set)
       ! Mulliken charges
       ALLOCATE(charges(number_qm_atoms,nspins),STAT=istat)
       CPASSERT(istat==0)
       !
       CALL mulliken_charges(matrix_p,matrix_s(1)%matrix,para_env,charges)
       !
       ALLOCATE(mcharge(number_qm_atoms),STAT=istat)
       CPASSERT(istat==0)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
          DO iatom=1,natom
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            mcharge(atom_a) = zeff - SUM(charges(atom_a,1:nspins))
          END DO
       END DO
       DEALLOCATE(charges,STAT=istat)
       CPASSERT(istat==0)

       ALLOCATE(qpot(number_qm_atoms),STAT=istat)
       CPASSERT(istat==0)
       qpot = 0.0_dp
       ALLOCATE(Forces_QM(3,number_qm_atoms),stat=istat)
       CPASSERT(istat==0)
       Forces_QM= 0.0_dp

       ! calculate potential and forces from classical charges
       iqm = 0
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list)
          NULLIFY(dftb_kind)
          CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind,&
               defined=defined,eta=eta_a,natorb=natorb)
          ! use mm charge smearing for non-scc cases
          IF(.NOT.dftb_control%self_consistent) eta_a(0) = eta_mm
          IF (.NOT.defined .OR. natorb < 1) CYCLE
          DO i = 1, SIZE(list)
             iatom = list(i)
             iqm = iqm + 1
             CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%Potentials,particles_mm,&
               qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,&
               qmmm_env%spherical_cutoff,particles_qm)
             CALL build_mm_dpot(mcharge(iatom),0,eta_a(0),qmmm_env%Potentials,particles_mm,&
               qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,Forces,Forces_QM(:,iqm),&
               qmmm_env%spherical_cutoff,particles_qm)
             ! Possibly added charges
             IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%added_charges%potentials,&
                  qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                  qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,qmmm_env%spherical_cutoff,&
                  particles_qm)
                CALL build_mm_dpot(mcharge(iatom),0,eta_a(0),qmmm_env%added_charges%potentials,&
                  qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                  qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,Forces_added_charges,&
                  Forces_QM(:,iqm),qmmm_env%spherical_cutoff,particles_qm)
             END IF
          END DO
       END DO

       ! Transfer QM gradients to the QM particles..
       iqm = 0
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
          NULLIFY(dftb_kind)
          CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind,defined=defined,natorb=natorb)
          IF (.NOT.defined .OR. natorb < 1) CYCLE
          DO i = 1, SIZE(list)
             iqm = iqm + 1
             iatom = qmmm_env%qm_atom_index(list(i))
             particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:,iqm)
          END DO
       END DO

       ! derivatives from qm charges
       Forces_QM = 0.0_dp
       IF ( SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
               alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
       END IF
       !
       CALL cp_dbcsr_iterator_start(iter, matrix_s(1)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
          !
          IF(iatom==jatom) CYCLE
          !
          gmij = -0.5_dp*(qpot(iatom) + qpot(jatom))
          NULLIFY(pblock)
          CALL cp_dbcsr_get_block_p(matrix=matrix_p(1)%matrix,&
               row=iatom,col=jatom,block=pblock,found=found)
          CPASSERT(found)
          DO i=1,3
             NULLIFY(dsblock)
             CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i)%matrix,&
                  row=iatom,col=jatom,block=dsblock,found=found)
             CPASSERT(found)
             fi = -2.0_dp*gmij*SUM(pblock*dsblock)
             Forces_QM(i,iatom) = Forces_QM(i,iatom) + fi
             Forces_QM(i,jatom) = Forces_QM(i,jatom) - fi
          END DO
       END DO
       CALL cp_dbcsr_iterator_stop(iter)
       !
       IF ( SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
                            alpha_scalar=1.0_dp,beta_scalar=-1.0_dp)
       END IF
       !
       ! Transfer QM gradients to the QM particles..
       CALL mp_sum(Forces_QM, para_env%group)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
          DO i = 1, SIZE(list)
             iqm = list(i)
             iatom = qmmm_env%qm_atom_index(iqm)
             particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:,iqm)
          END DO
       END DO
       !
       DEALLOCATE(mcharge,STAT=istat)
       CPASSERT(istat==0)
       !
       ! MM forces will be handled directly from the QMMM module in the same way
       ! as for GPW/GAPW methods
       DEALLOCATE(Forces_QM,stat=istat)
       CPASSERT(istat==0)
       DEALLOCATE(qpot,STAT=istat)
       CPASSERT(istat==0)

       CALL cp_dbcsr_deallocate_matrix_set ( matrix_s)

    END IF
    CALL timestop(handle)

  END SUBROUTINE deriv_dftb_qmmm_matrix

! *****************************************************************************
!> \brief Constructs the derivative w.r.t. 1-el DFTB hamiltonian QMMM terms
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \param calc_force ...
!> \param Forces ...
!> \param Forces_added_charges ...
!> \author JGH 10.2014 [created]
! *****************************************************************************
  SUBROUTINE deriv_dftb_qmmm_matrix_pc(qs_env,qmmm_env,particles_mm,mm_cell,para_env,&
       calc_force, Forces, Forces_added_charges)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'deriv_dftb_qmmm_matrix_pc', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, blk, do_ipol, ewald_type, handle, i, iatom, ikind, &
      imm, imp, indmm, ipot, iqm, istat, jatom, natom, natorb, nkind, nmm, &
      nspins, number_qm_atoms
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: defined, do_multipoles, &
                                                failure, found
    REAL(KIND=dp)                            :: alpha, fi, gmij, zeff
    REAL(KIND=dp), DIMENSION(0:3)            :: eta_a
    REAL(KIND=dp), DIMENSION(2)              :: rcutoff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges_mm, mcharge, qpot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, dsblock, Forces_MM, &
                                                Forces_QM, pblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: atoms_mm, particles_qm
    TYPE(qmmm_pot_type), POINTER             :: Pot
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: ewald_section, &
                                                poisson_section, print_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (calc_force) THEN
       NULLIFY (rho, atomic_kind_set, qs_kind_set, particles_qm)
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       atomic_kind_set=atomic_kind_set,&
                       qs_kind_set=qs_kind_set,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       dft_control=dft_control,&
                       particle_set=particles_qm,&
                       natom=number_qm_atoms)
       dftb_control => dft_control%qs_control%dftb_control

       NULLIFY(matrix_s)
       CALL build_dftb_overlap(qs_env,1,matrix_s)
       CALL qs_rho_get(rho, rho_ao=matrix_p)

       nspins = dft_control%nspins
       nkind  = SIZE(atomic_kind_set)
      ! Mulliken charges
       ALLOCATE(charges(number_qm_atoms,nspins),STAT=istat)
       CPASSERT(istat==0)
       !
       CALL mulliken_charges(matrix_p,matrix_s(1)%matrix,para_env,charges)
       !
       ALLOCATE(mcharge(number_qm_atoms),STAT=istat)
       CPASSERT(istat==0)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
          DO iatom=1,natom
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            mcharge(atom_a) = zeff - SUM(charges(atom_a,1:nspins))
          END DO
       END DO
       DEALLOCATE(charges,STAT=istat)
       CPASSERT(istat==0)

       ALLOCATE(qpot(number_qm_atoms),STAT=istat)
       CPASSERT(istat==0)
       qpot = 0.0_dp
       ALLOCATE(Forces_QM(3,number_qm_atoms),stat=istat)
       CPASSERT(istat==0)
       Forces_QM= 0.0_dp

       ! Create Ewald environments
       poisson_section => section_vals_get_subs_vals(qs_env%input,"MM%POISSON")
       CALL ewald_env_create(ewald_env,para_env)
       CALL ewald_env_set(ewald_env,poisson_section=poisson_section)
       ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD")
       CALL read_ewald_section(ewald_env,ewald_section)
       print_section => section_vals_get_subs_vals(qs_env%input,"PRINT%GRID_INFORMATION")
       CALL ewald_pw_create(ewald_pw,ewald_env,mm_cell,mm_cell,print_section=print_section)

       CALL ewald_env_get(ewald_env,ewald_type=ewald_type,do_multipoles=do_multipoles,do_ipol=do_ipol)
       IF(do_multipoles) CALL cp_unimplemented_error(fromWhere=routineP, &
          message="No multipole force fields allowed in DFTB QM/MM",error_level=cp_fatal_level)
       IF(do_ipol /= do_fist_pol_none) CALL cp_unimplemented_error(fromWhere=routineP, &
          message="No polarizable force fields allowed in DFTB QM/MM",error_level=cp_fatal_level)

       SELECT CASE(ewald_type)
          CASE(do_ewald_pme)
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="PME Ewald type not implemented for DFTB/QMMM",error_level=cp_fatal_level)
          CASE(do_ewald_ewald,do_ewald_spme)
             DO ipot = 1,SIZE(qmmm_env%Potentials)
                Pot => qmmm_env%Potentials(ipot)%Pot
                nmm = SIZE(Pot%mm_atom_index)
                ! get a 'clean' mm particle set
                NULLIFY(atoms_mm)
                CALL allocate_particle_set(atoms_mm,nmm)
                ALLOCATE(charges_mm(nmm),STAT=istat)
                CPASSERT(istat==0)
                DO Imp=1,nmm
                   Imm = Pot%mm_atom_index(Imp)
                   IndMM = qmmm_env%mm_atom_index(Imm)
                   atoms_mm(Imp)%r = particles_mm(IndMM)%r
                   atoms_mm(Imp)%atomic_kind => particles_mm(IndMM)%atomic_kind
                   charges_mm(Imp) = qmmm_env%mm_atom_chrg(Imm)
                END DO
                ! force array for mm atoms
                ALLOCATE(Forces_MM(3,nmm),stat=istat)
                CPASSERT(istat==0)
                Forces_MM= 0.0_dp
                IF(ewald_type == do_ewald_ewald) THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Ewald not implemented for DFTB/QMMM",error_level=cp_fatal_level)
                ELSE IF(ewald_type == do_ewald_spme) THEN
                   ! spme electrostatic potential
                   CALL spme_potential(ewald_env,ewald_pw,mm_cell,atoms_mm,charges_mm, &
                                       particles_qm,qpot)
                   ! forces QM
                   CALL spme_forces(ewald_env,ewald_pw,mm_cell,atoms_mm,charges_mm, &
                                    particles_qm,mcharge,Forces_QM)
                   ! forces MM
                   CALL spme_forces(ewald_env,ewald_pw,mm_cell,particles_qm,mcharge,&
                                    atoms_mm,charges_mm,Forces_MM)
                END IF
                CALL deallocate_particle_set(atoms_mm)
                DEALLOCATE(charges_mm,STAT=istat)
                CPASSERT(istat==0)
                ! transfer MM forces
                CALL mp_sum(Forces_MM,para_env%group)
                DO Imp=1,nmm
                   Imm = Pot%mm_atom_index(Imp)
                   Forces(:,Imm) = Forces(:,Imm) - Forces_MM(:,Imp)
                END DO
                DEALLOCATE(Forces_MM,stat=istat)
                CPASSERT(istat==0)
             END DO

             IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                DO ipot = 1,SIZE(qmmm_env%added_charges%Potentials)
                   Pot => qmmm_env%added_charges%Potentials(ipot)%Pot
                   nmm = SIZE(Pot%mm_atom_index)
                   ! get a 'clean' mm particle set
                   NULLIFY(atoms_mm)
                   CALL allocate_particle_set(atoms_mm,nmm)
                   ALLOCATE(charges_mm(nmm),STAT=istat)
                   CPASSERT(istat==0)
                   DO Imp=1,nmm
                      Imm = Pot%mm_atom_index(Imp)
                      IndMM = qmmm_env%added_charges%mm_atom_index(Imm)
                      atoms_mm(Imp)%r = qmmm_env%added_charges%added_particles(IndMM)%r
                      atoms_mm(Imp)%atomic_kind => qmmm_env%added_charges%added_particles(IndMM)%atomic_kind
                      charges_mm(Imp) = qmmm_env%added_charges%mm_atom_chrg(Imm)
                   END DO
                   ! force array for mm atoms
                   ALLOCATE(Forces_MM(3,nmm),stat=istat)
                   CPASSERT(istat==0)
                   Forces_MM= 0.0_dp
                   IF(ewald_type == do_ewald_ewald) THEN
                      CALL cp_unimplemented_error(fromWhere=routineP, &
                           message="Ewald not implemented for DFTB/QMMM",error_level=cp_fatal_level)
                   ELSE IF(ewald_type == do_ewald_spme) THEN
                      ! spme electrostatic potential
                      CALL spme_potential(ewald_env,ewald_pw,mm_cell,atoms_mm,&
                                          charges_mm,particles_qm,qpot)
                      ! forces QM
                      CALL spme_forces(ewald_env,ewald_pw,mm_cell,atoms_mm,charges_mm,&
                                       particles_qm,mcharge,Forces_QM)
                      ! forces MM
                      CALL spme_forces(ewald_env,ewald_pw,mm_cell,particles_qm,mcharge,&
                                       atoms_mm,charges_mm,Forces_MM)
                   END IF
                   CALL deallocate_particle_set(atoms_mm)
                   ! transfer MM forces
                   CALL mp_sum(Forces_MM,para_env%group)
                   DO Imp=1,nmm
                      Imm = Pot%mm_atom_index(Imp)
                      Forces_added_charges(:,Imm) = Forces_added_charges(:,Imm) - Forces_MM(:,Imp)
                   END DO
                   DEALLOCATE(Forces_MM,stat=istat)
                   CPASSERT(istat==0)
                END DO
             END IF
             CALL mp_sum(qpot,para_env%group)
             CALL mp_sum(Forces_QM, para_env%group)
             CPASSERT(istat==0)
             ! Add Ewald and DFTB short range corrections
             ! This is effectively using a minimum image convention!
             ! Set rcutoff to values compatible with alpha Ewald
             CALL ewald_env_get(ewald_env,rcut=rcutoff(1),alpha=alpha)
             rcutoff(2) = 0.025_dp*rcutoff(1)
             rcutoff(1) = 2.0_dp*rcutoff(1)
             nkind = SIZE(atomic_kind_set)
             iqm = 0
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list)
                NULLIFY(dftb_kind)
                CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
                CALL get_dftb_atom_param(dftb_kind,&
                     defined=defined,eta=eta_a,natorb=natorb)
                ! use mm charge smearing for non-scc cases
                IF(.NOT.dftb_control%self_consistent) eta_a(0) = eta_mm
                IF (.NOT.defined .OR. natorb < 1) CYCLE
                DO i = 1, SIZE(list)
                   iatom = list(i)
                   iqm = iqm + 1
                   CALL build_mm_pot(qpot(iatom),1,eta_a(0),qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,rcutoff,particles_qm)
                   CALL build_mm_dpot(mcharge(iatom),1,eta_a(0),qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,Forces,Forces_QM(:,iqm),&
                     rcutoff,particles_qm)
                   CALL build_mm_pot(qpot(iatom),2,alpha,qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,rcutoff,particles_qm)
                   CALL build_mm_dpot(mcharge(iatom),2,alpha,qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,Forces,Forces_QM(:,iqm),&
                     rcutoff,particles_qm)
                   ! Possibly added charges
                   IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                      CALL build_mm_pot(qpot(iatom),1,eta_a(0),qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,rcutoff,particles_qm)
                      CALL build_mm_dpot(mcharge(iatom),1,eta_a(0),qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,Forces_added_charges,Forces_QM(:,iqm),&
                        rcutoff,particles_qm)
                      CALL build_mm_pot(qpot(iatom),2,alpha,qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,rcutoff,particles_qm)
                      CALL build_mm_dpot(mcharge(iatom),2,alpha,qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,Forces_added_charges,Forces_QM(:,iqm),&
                        rcutoff,particles_qm)
                   END IF
                END DO
             END DO

          CASE(do_ewald_none)
             ! Simply summing up charges with 1/R (DFTB corrected)
             ! calculate potential and forces from classical charges
             iqm = 0
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list)
                NULLIFY(dftb_kind)
                CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
                CALL get_dftb_atom_param(dftb_kind,&
                     defined=defined,eta=eta_a,natorb=natorb)
                ! use mm charge smearing for non-scc cases
                IF(.NOT.dftb_control%self_consistent) eta_a(0) = eta_mm
                IF (.NOT.defined .OR. natorb < 1) CYCLE
                DO i = 1, SIZE(list)
                   iatom = list(i)
                   iqm = iqm + 1
                   CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,&
                     qmmm_env%spherical_cutoff,particles_qm)
                   CALL build_mm_dpot(mcharge(iatom),0,eta_a(0),qmmm_env%Potentials,particles_mm,&
                     qmmm_env%mm_atom_chrg,qmmm_env%mm_atom_index,mm_cell,iatom,Forces,Forces_QM(:,iqm),&
                     qmmm_env%spherical_cutoff,particles_qm)
                   ! Possibly added charges
                   IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                      CALL build_mm_pot(qpot(iatom),0,eta_a(0),qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,qmmm_env%spherical_cutoff,&
                        particles_qm)
                      CALL build_mm_dpot(mcharge(iatom),0,eta_a(0),qmmm_env%added_charges%potentials,&
                        qmmm_env%added_charges%added_particles,qmmm_env%added_charges%mm_atom_chrg,&
                        qmmm_env%added_charges%mm_atom_index,mm_cell,iatom,Forces_added_charges,&
                        Forces_QM(:,iqm),qmmm_env%spherical_cutoff,particles_qm)
                   END IF
                END DO
             END DO
          CASE DEFAULT
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Unknown Ewald type!",error_level=cp_fatal_level)
       END SELECT

       ! Transfer QM gradients to the QM particles..
       iqm = 0
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
          NULLIFY(dftb_kind)
          CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind,defined=defined,natorb=natorb)
          IF (.NOT.defined .OR. natorb < 1) CYCLE
          DO i = 1, SIZE(list)
             iqm = iqm + 1
             iatom = qmmm_env%qm_atom_index(list(i))
             particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:,iqm)
          END DO
       END DO

       ! derivatives from qm charges
       Forces_QM = 0.0_dp
       IF ( SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
               alpha_scalar=1.0_dp,beta_scalar=1.0_dp)
       END IF
       !
       CALL cp_dbcsr_iterator_start(iter, matrix_s(1)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
          !
          IF(iatom==jatom) CYCLE
          !
          gmij = -0.5_dp*(qpot(iatom) + qpot(jatom))
          NULLIFY(pblock)
          CALL cp_dbcsr_get_block_p(matrix=matrix_p(1)%matrix,&
               row=iatom,col=jatom,block=pblock,found=found)
          CPASSERT(found)
          DO i=1,3
             NULLIFY(dsblock)
             CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i)%matrix,&
                  row=iatom,col=jatom,block=dsblock,found=found)
             CPASSERT(found)
             fi = -2.0_dp*gmij*SUM(pblock*dsblock)
             Forces_QM(i,iatom) = Forces_QM(i,iatom) + fi
             Forces_QM(i,jatom) = Forces_QM(i,jatom) - fi
          END DO
       END DO
       CALL cp_dbcsr_iterator_stop(iter)
       !
       IF ( SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
                            alpha_scalar=1.0_dp,beta_scalar=-1.0_dp)
       END IF
       !
       ! Transfer QM gradients to the QM particles..
       CALL mp_sum(Forces_QM, para_env%group)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
          DO i = 1, SIZE(list)
             iqm = list(i)
             iatom = qmmm_env%qm_atom_index(iqm)
             particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:,iqm)
          END DO
       END DO
       !
       DEALLOCATE(mcharge,STAT=istat)
       CPASSERT(istat==0)
       !
       ! MM forces will be handled directly from the QMMM module in the same way
       ! as for GPW/GAPW methods
       DEALLOCATE(Forces_QM,stat=istat)
       CPASSERT(istat==0)
       DEALLOCATE(qpot,STAT=istat)
       CPASSERT(istat==0)

       ! Release Ewald environment
       CALL ewald_env_release(ewald_env)
       CALL ewald_pw_release(ewald_pw)

       CALL cp_dbcsr_deallocate_matrix_set ( matrix_s)

    END IF

    CALL timestop(handle)

  END SUBROUTINE deriv_dftb_qmmm_matrix_pc

! *****************************************************************************
!> \brief ...
!> \param qpot ...
!> \param pot_type ...
!> \param qm_alpha ...
!> \param potentials ...
!> \param particles_mm ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_cell ...
!> \param IndQM ...
!> \param qmmm_spherical_cutoff ...
!> \param particles_qm ...
! *****************************************************************************
  SUBROUTINE build_mm_pot(qpot, pot_type, qm_alpha, potentials,&
       particles_mm, mm_charges, mm_atom_index, mm_cell, IndQM, &
       qmmm_spherical_cutoff, particles_qm)

    REAL(KIND=dp), INTENT(INOUT)             :: qpot
    INTEGER, INTENT(IN)                      :: pot_type
    REAL(KIND=dp), INTENT(IN)                :: qm_alpha
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, INTENT(IN)                      :: IndQM
    REAL(KIND=dp), INTENT(IN)                :: qmmm_spherical_cutoff(2)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm

    CHARACTER(len=*), PARAMETER :: routineN = 'build_mm_pot', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: qsmall = 1.0e-15_dp

    INTEGER                                  :: handle, Imm, Imp, IndMM, Ipot
    REAL(KIND=dp)                            :: dr, qeff, rt1, rt2, rt3, &
                                                sph_chrg_factor, sr
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc, rij
    TYPE(qmmm_pot_type), POINTER             :: Pot

    CALL timeset(routineN,handle)
    ! Loop Over MM atoms
    ! Loop over Pot stores atoms with the same charge
    MainLoopPot: DO Ipot = 1, SIZE(Potentials)
       Pot    => Potentials(Ipot)%Pot
       ! Loop over atoms belonging to this type
       LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
          Imm = Pot%mm_atom_index(Imp)
          IndMM = mm_atom_index(Imm)
          r_pbc=pbc(particles_mm(IndMM)%r-particles_qm(IndQM)%r, mm_cell)
          rt1= r_pbc(1)
          rt2= r_pbc(2)
          rt3= r_pbc(3)
          rij = (/rt1,rt2,rt3/)
          dr = SQRT(SUM(rij**2))
          qeff = mm_charges(Imm)
          ! Computes the screening factor for the spherical cutoff (if defined)
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, rij, sph_chrg_factor)
             qeff = qeff * sph_chrg_factor
          END IF
          IF (ABS(qeff)<=qsmall) CYCLE
          IF(dr > rtiny) THEN
             IF (pot_type==0) THEN
                sr = gamma_rab_sr(dr,qm_alpha,eta_mm,0.0_dp)
                qpot = qpot + qeff*(1.0_dp/dr - sr)
             ELSE IF (pot_type==1) THEN
                sr = gamma_rab_sr(dr,qm_alpha,eta_mm,0.0_dp)
                qpot = qpot - qeff*sr
             ELSE IF (pot_type==2) THEN
                sr = erfc(qm_alpha*dr)/dr
                qpot = qpot + qeff*sr
             ELSE 
                STOP
             END IF
          END IF
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(handle)
  END SUBROUTINE build_mm_pot

! *****************************************************************************
!> \brief ...
!> \param qcharge ...
!> \param pot_type ...
!> \param qm_alpha ...
!> \param potentials ...
!> \param particles_mm ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_cell ...
!> \param IndQM ...
!> \param forces ...
!> \param forces_qm ...
!> \param qmmm_spherical_cutoff ...
!> \param particles_qm ...
! *****************************************************************************
  SUBROUTINE build_mm_dpot(qcharge, pot_type, qm_alpha, potentials,&
       particles_mm, mm_charges, mm_atom_index, mm_cell, IndQM, &
       forces, forces_qm, qmmm_spherical_cutoff, particles_qm)

    REAL(KIND=dp), INTENT(IN)                :: qcharge
    INTEGER, INTENT(IN)                      :: pot_type
    REAL(KIND=dp), INTENT(IN)                :: qm_alpha
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, INTENT(IN)                      :: IndQM
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: forces_qm
    REAL(KIND=dp), INTENT(IN)                :: qmmm_spherical_cutoff(2)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm

    CHARACTER(len=*), PARAMETER :: routineN = 'build_mm_dpot', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: qsmall = 1.0e-15_dp

    INTEGER                                  :: handle, Imm, Imp, IndMM, Ipot
    REAL(KIND=dp)                            :: dr, drm, drp, dsr, fsr, qeff, &
                                                rt1, rt2, rt3, sph_chrg_factor
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, r_pbc, rij
    TYPE(qmmm_pot_type), POINTER             :: Pot

    CALL timeset(routineN,handle)
    ! Loop Over MM atoms
    ! Loop over Pot stores atoms with the same charge
    MainLoopPot: DO Ipot = 1, SIZE(Potentials)
       Pot    => Potentials(Ipot)%Pot
       ! Loop over atoms belonging to this type
       LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
          Imm = Pot%mm_atom_index(Imp)
          IndMM = mm_atom_index(Imm)
          r_pbc=pbc(particles_mm(IndMM)%r-particles_qm(IndQM)%r, mm_cell)
          rt1= r_pbc(1)
          rt2= r_pbc(2)
          rt3= r_pbc(3)
          rij = (/rt1,rt2,rt3/)
          dr = SQRT(SUM(rij**2))
          qeff = mm_charges(Imm)
          ! Computes the screening factor for the spherical cutoff (if defined)
          ! We neglect derivative of cutoff function for gradients!!!
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, rij, sph_chrg_factor)
             qeff = qeff * sph_chrg_factor
          END IF
          IF (ABS(qeff)<=qsmall) CYCLE
          IF(dr > rtiny) THEN
             drp = dr + ddrmm
             drm = dr - ddrmm
             IF (pot_type==0) THEN
                dsr = 0.5_dp*(gamma_rab_sr(drp,qm_alpha,eta_mm,0.0_dp)-&
                              gamma_rab_sr(drm,qm_alpha,eta_mm,0.0_dp))/ddrmm
                fsr = qeff*qcharge*(-1.0_dp/(dr*dr) - dsr)
             ELSE IF (pot_type==1) THEN
                dsr = 0.5_dp*(gamma_rab_sr(drp,qm_alpha,eta_mm,0.0_dp)-&
                              gamma_rab_sr(drm,qm_alpha,eta_mm,0.0_dp))/ddrmm
                fsr = -qeff*qcharge*dsr
             ELSE IF (pot_type==2) THEN
                dsr = 0.5_dp*(erfc(qm_alpha*drp)/drp-erfc(qm_alpha*drm)/drm)/ddrmm
                fsr = qeff*qcharge*dsr
             ELSE
                STOP
             END IF
             force_ab = -fsr * rij/dr
          ELSE
             force_ab = 0.0_dp
          END IF
          ! The array of QM forces are really the forces
          forces_qm(:)  = forces_qm(:)  - force_ab
          ! The one of MM atoms are instead gradients
          forces(:,Imm) = forces(:,Imm) - force_ab
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(handle)
  END SUBROUTINE build_mm_dpot

! *****************************************************************************
!> \brief ...
!> \param block ...
!> \param smatij ...
!> \param smatji ...
!> \param rij ...
!> \param ngrd ...
!> \param ngrdcut ...
!> \param dgrd ...
!> \param llm ...
!> \param lmaxi ...
!> \param lmaxj ...
!> \param irow ...
!> \param iatom ...
! *****************************************************************************
  SUBROUTINE compute_block_sk(block,smatij,smatji,rij,ngrd,ngrdcut,dgrd,llm,lmaxi,lmaxj,irow,iatom)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block, smatij, smatji
    REAL(KIND=dp), DIMENSION(3)              :: rij
    INTEGER                                  :: ngrd, ngrdcut
    REAL(KIND=dp)                            :: dgrd
    INTEGER                                  :: llm, lmaxi, lmaxj, irow, iatom

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_block_sk', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr
    REAL(KIND=dp), DIMENSION(20)             :: skabij, skabji

     dr = SQRT(SUM(rij(:)**2))
     CALL getskz(smatij,skabij,dr,ngrd,ngrdcut,dgrd,llm)
     CALL getskz(smatji,skabji,dr,ngrd,ngrdcut,dgrd,llm)
     IF ( irow == iatom ) THEN
       CALL turnsk(block,skabji,skabij,rij,dr,lmaxi,lmaxj)
     ELSE
       CALL turnsk(block,skabij,skabji,-rij,dr,lmaxj,lmaxi)
     END IF

  END SUBROUTINE compute_block_sk

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nderivative ...
!> \param matrices ...
!> \param mnames ...
!> \param sab_nl ...
! *****************************************************************************
  SUBROUTINE setup_matrices1(qs_env,nderivative,matrices,mnames,sab_nl)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrices
    CHARACTER(LEN=*)                         :: mnames
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_matrices1', &
      routineP = moduleN//':'//routineN

    CHARACTER(1)                             :: symmetry_type
    CHARACTER(LEN=default_string_length)     :: matnames
    INTEGER                                  :: i, istat, natom, &
                                                neighbor_list_id, nkind, &
                                                nmat, nsgf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    NULLIFY(particle_set, atomic_kind_set)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    dbcsr_dist=dbcsr_dist,&
                    neighbor_list_id=neighbor_list_id)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    CALL get_qs_kind_set(qs_kind_set,nsgf=nsgf)

    ALLOCATE (first_sgf(natom),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPASSERT(istat==0)

    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)

    nmat = 0
    IF(nderivative==0)nmat=1
    IF(nderivative==1)nmat=4
    IF(nderivative==2)nmat=10
    CPASSERT(nmat>0)

    ALLOCATE (row_blk_sizes(natom), STAT=istat)
    CPASSERT(istat==0)
    CALL convert_offsets_to_sizes (first_sgf, row_blk_sizes, last_sgf)

    CALL cp_dbcsr_allocate_matrix_set(matrices,nmat)

    ! Up to 2nd derivative take care to get the symmetries correct
    DO i=1,nmat
       IF(i.gt.1)THEN
          matnames=TRIM(mnames)//" DERIVATIVE MATRIX DFTB"
          symmetry_type=dbcsr_type_antisymmetric
          IF(i.gt.4)symmetry_type=dbcsr_type_symmetric
       ELSE
          symmetry_type=dbcsr_type_symmetric
          matnames=TRIM(mnames)//" MATRIX DFTB"
       END IF
       ALLOCATE(matrices(i)%matrix)
       CALL cp_dbcsr_init(matrices(i)%matrix)
       CALL cp_dbcsr_create(matrix=matrices(i)%matrix, &
            name=TRIM(matnames), &
            dist=dbcsr_dist, matrix_type=symmetry_type,&
            row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
            nze=0, mutable_work=.TRUE.)
       CALL cp_dbcsr_alloc_block_from_nbl(matrices(i)%matrix,sab_nl)
    END DO

    DEALLOCATE (first_sgf,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE (last_sgf,STAT=istat)
    CPASSERT(istat==0)

    DEALLOCATE(row_blk_sizes)

  END SUBROUTINE setup_matrices1

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nderivative ...
!> \param nimg ...
!> \param matrices ...
!> \param mnames ...
!> \param sab_nl ...
! *****************************************************************************
  SUBROUTINE setup_matrices2(qs_env,nderivative,nimg,matrices,mnames,sab_nl)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: nderivative, nimg
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrices
    CHARACTER(LEN=*)                         :: mnames
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_matrices2', &
      routineP = moduleN//':'//routineN

    CHARACTER(1)                             :: symmetry_type
    CHARACTER(LEN=default_string_length)     :: matnames
    INTEGER                                  :: i, img, istat, natom, &
                                                neighbor_list_id, nkind, &
                                                nmat, nsgf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    NULLIFY(particle_set, atomic_kind_set)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    dbcsr_dist=dbcsr_dist,&
                    neighbor_list_id=neighbor_list_id)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    CALL get_qs_kind_set(qs_kind_set,nsgf=nsgf)

    ALLOCATE (first_sgf(natom),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPASSERT(istat==0)

    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)

    nmat = 0
    IF(nderivative==0)nmat=1
    IF(nderivative==1)nmat=4
    IF(nderivative==2)nmat=10
    CPASSERT(nmat>0)

    ALLOCATE (row_blk_sizes(natom), STAT=istat)
    CPASSERT(istat==0)
    CALL convert_offsets_to_sizes (first_sgf, row_blk_sizes, last_sgf)

    CALL cp_dbcsr_allocate_matrix_set(matrices,nmat,nimg)

    ! Up to 2nd derivative take care to get the symmetries correct
    DO img=1,nimg
       DO i=1,nmat
          IF(i.gt.1)THEN
             matnames=TRIM(mnames)//" DERIVATIVE MATRIX DFTB"
             symmetry_type=dbcsr_type_antisymmetric
             IF(i.gt.4)symmetry_type=dbcsr_type_symmetric
          ELSE
             symmetry_type=dbcsr_type_symmetric
             matnames=TRIM(mnames)//" MATRIX DFTB"
          END IF
          ALLOCATE(matrices(i,img)%matrix)
          CALL cp_dbcsr_init(matrices(i,img)%matrix)
          CALL cp_dbcsr_create(matrix=matrices(i,img)%matrix, &
               name=TRIM(matnames), &
               dist=dbcsr_dist, matrix_type=symmetry_type,&
               row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
               nze=0, mutable_work=.TRUE.)
          CALL cp_dbcsr_alloc_block_from_nbl(matrices(i,img)%matrix,sab_nl)
       END DO
    END DO

    DEALLOCATE (first_sgf,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE (last_sgf,STAT=istat)
    CPASSERT(istat==0)

    DEALLOCATE(row_blk_sizes)

  END SUBROUTINE setup_matrices2

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nmat ...
!> \param gammat ...
!> \param sab_nl ...
! *****************************************************************************
  SUBROUTINE setup_gamma(qs_env,nmat,gammat,sab_nl)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: nmat
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: gammat
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_gamma', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, istat, natom, ncol, &
                                                neighbor_list_id, nrow
    INTEGER, DIMENSION(:), POINTER           :: felem, lelem, row_blk_sizes
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

       CALL get_qs_env(qs_env=qs_env,&
                       particle_set=particle_set,&
                       neighbor_list_id=neighbor_list_id,&
                       dbcsr_dist=dbcsr_dist)
       natom = SIZE(particle_set)
       nrow = natom
       ncol = natom
       ALLOCATE (felem(natom),lelem(natom),STAT=istat)
       CPASSERT(istat==0)
       DO iatom = 1, natom
          felem(iatom) = iatom
          lelem(iatom) = iatom
       ENDDO

       ALLOCATE (row_blk_sizes(natom), STAT=istat)
       CPASSERT(istat==0)
       CALL convert_offsets_to_sizes (felem, row_blk_sizes, lelem)
 
       CALL cp_dbcsr_allocate_matrix_set(gammat,nmat)
       ALLOCATE(gammat(1)%matrix,STAT=istat)
       CPASSERT(istat==0)
       CALL cp_dbcsr_init(gammat(1)%matrix)
 
       CALL cp_dbcsr_create(matrix=gammat(1)%matrix, &
            name="GAMMA MATRIX", &
            dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric,&
            row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
            nze=0, mutable_work=.TRUE.)
 
       DO i=2,nmat
          ALLOCATE(gammat(i)%matrix)
          CALL cp_dbcsr_init(gammat(i)%matrix)
 
          CALL cp_dbcsr_create(matrix=gammat(i)%matrix, &
               name="DERIVATIVE GAMMA MATRIX", &
               dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
               row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
               nze=0, mutable_work=.TRUE.)
       END DO
 
       DEALLOCATE (row_blk_sizes,felem,lelem,STAT=istat)
       CPASSERT(istat==0)

       ! setup the matrices using the neighbor list
       DO i=1,nmat
          CALL cp_dbcsr_alloc_block_from_nbl(gammat(i)%matrix,sab_nl)
       END DO

  END SUBROUTINE setup_gamma

! *****************************************************************************
!> \brief Gets matrix elements on z axis, as they are stored in the tables
!> \param slakotab ...
!> \param skpar ...
!> \param dx ...
!> \param ngrd ...
!> \param ngrdcut ...
!> \param dgrd ...
!> \param llm ...
!> \author 07. Feb. 2004, TH
! *****************************************************************************
  SUBROUTINE getskz(slakotab,skpar,dx,ngrd,ngrdcut,dgrd,llm)
    REAL(dp), INTENT(in)                     :: slakotab(:,:), dx
    INTEGER, INTENT(in)                      :: ngrd, ngrdcut
    REAL(dp), INTENT(in)                     :: dgrd
    INTEGER, INTENT(in)                      :: llm
    REAL(dp), INTENT(out)                    :: skpar(llm)

    INTEGER                                  :: clgp

    skpar = 0._dp
    !
    ! Determine closest grid point
    !
    clgp = NINT(dx/dgrd)
    !
    ! Screen elements which are too far away
    !
    IF (clgp > ngrdcut) RETURN
    !
    ! The grid point is either contained in the table --> matrix element
    ! can be interpolated, or it is outside the table --> matrix element
    ! needs to be extrapolated.
    !
    IF (clgp > ngrd) THEN
      !
      ! Extrapolate external matrix elements if table does not finish with zero
      !
      CALL extrapol(slakotab,skpar,dx,ngrd,dgrd,llm)
    ELSE
      !
      ! Interpolate tabulated matrix elements
      !
      CALL interpol(slakotab,skpar,dx,ngrd,dgrd,llm,clgp)
    END IF
  END SUBROUTINE getskz

! *****************************************************************************
!> \brief ...
!> \param slakotab ...
!> \param skpar ...
!> \param dx ...
!> \param ngrd ...
!> \param dgrd ...
!> \param llm ...
!> \param clgp ...
! *****************************************************************************
  SUBROUTINE interpol(slakotab,skpar,dx,ngrd,dgrd,llm,clgp)
    REAL(dp), INTENT(in)                     :: slakotab(:,:), dx
    INTEGER, INTENT(in)                      :: ngrd
    REAL(dp), INTENT(in)                     :: dgrd
    INTEGER, INTENT(in)                      :: llm
    REAL(dp), INTENT(out)                    :: skpar(llm)
    INTEGER, INTENT(in)                      :: clgp

    INTEGER                                  :: fgpm, k, l, lgpm
    REAL(dp)                                 :: error, xa(max_inter), &
                                                ya(max_inter)

    lgpm = MIN(clgp+max_inter/2,ngrd)
    fgpm = lgpm - max_inter + 1
    DO k = 0,max_inter-1
      xa(k+1) = (fgpm+k)*dgrd
    END DO
    !
    ! Interpolate matrix elements for all orbitals
    !
    DO l = 1, llm
      !
      ! Read SK parameters from table
      !
      ya(1:max_inter) = slakotab(fgpm:lgpm,l)
      CALL polint(xa,ya,max_inter,dx,skpar(l),error)
    END DO
  END SUBROUTINE interpol

! *****************************************************************************
!> \brief ...
!> \param slakotab ...
!> \param skpar ...
!> \param dx ...
!> \param ngrd ...
!> \param dgrd ...
!> \param llm ...
! *****************************************************************************
  SUBROUTINE extrapol(slakotab,skpar,dx,ngrd,dgrd,llm)
    REAL(dp), INTENT(in)                     :: slakotab(:,:), dx
    INTEGER, INTENT(in)                      :: ngrd
    REAL(dp), INTENT(in)                     :: dgrd
    INTEGER, INTENT(in)                      :: llm
    REAL(dp), INTENT(out)                    :: skpar(llm)

    INTEGER                                  :: fgp, k, l, lgp, ntable, nzero
    REAL(dp)                                 :: error, xa(max_extra), &
                                                ya(max_extra)

    nzero = max_extra/3
    ntable = max_extra-nzero
    !
    ! Get the three last distances from the table
    !
    DO k = 1,ntable
      xa(k) = (ngrd-(max_extra-3)+k)*dgrd
    END DO
    DO k = 1,nzero
      xa(ntable+k) = (ngrd+k-1)*dgrd + slako_d0
      ya(ntable+k) = 0.0
    END DO
    !
    ! Extrapolate matrix elements for all orbitals
    !
    DO l = 1,llm
      !
      ! Read SK parameters from table
      !
      fgp = ngrd + 1 - (max_extra-3)
      lgp = ngrd
      ya(1:max_extra-3) = slakotab(fgp:lgp,l)
      CALL polint(xa,ya,max_extra,dx,skpar(l),error)
    END DO
  END SUBROUTINE extrapol

! *****************************************************************************
!> \brief   Turn matrix element from z-axis to orientation of dxv
!> \param mat ...
!> \param skab1 ...
!> \param skab2 ...
!> \param dxv ...
!> \param dx ...
!> \param lmaxa ...
!> \param lmaxb ...
!> \date    13. Jan 2004
!> \par Notes
!>          These routines are taken from an old TB code (unknown to TH).
!>          They are highly optimised and taken because they are time critical.
!>          They are explicit, so not recursive, and work up to d functions.
!>
!>          Set variables necessary for rotation of matrix elements
!>
!>          r_i^2/r, replicated in rr2(4:6) for index convenience later
!>          r_i/r, direction vector, rr(4:6) are replicated from 1:3
!>          lmax of A and B
!> \author  TH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE turnsk(mat,skab1,skab2,dxv,dx,lmaxa,lmaxb)
    REAL(dp), INTENT(inout)                  :: mat(:,:)
    REAL(dp), INTENT(in)                     :: skab1(:), skab2(:), dxv(3), dx
    INTEGER, INTENT(in)                      :: lmaxa, lmaxb

    INTEGER                                  :: lmaxab, minlmaxab
    REAL(dp)                                 :: rinv, rr(6), rr2(6)

    lmaxab = MAX(lmaxa,lmaxb)
    ! Determine l quantum limits.
    IF (lmaxab.gt.2)  STOP 'lmax=2'
    minlmaxab = MIN(lmaxa,lmaxb)
    !
    ! s-s interaction
    !
    CALL skss(skab1,mat)
    !
    IF (lmaxab.le.0) RETURN
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
    IF (minlmaxab.ge.1) THEN
      CALL skpp(skab1,mat,iptr(:,:,:,lmaxa,lmaxb))
      CALL sksp(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
      CALL sksp(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
    ELSE
      IF (lmaxb.ge.1) THEN
        CALL sksp(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
      ELSE
        CALL sksp(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
      END IF
    END IF
    !
    ! If there is only s-p interaction we have finished
    !
    IF (lmaxab.le.1) RETURN
    !
    ! at least one atom has d functions
    !
    IF (minlmaxab.eq.2) THEN
      !
      ! in case both atoms have d functions
      !
      CALL skdd(skab1,mat,iptr(:,:,:,lmaxa,lmaxb))
      CALL sksd(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
      CALL sksd(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
      CALL skpd(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
      CALL skpd(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
    ELSE
      !
      ! One atom has d functions, the other has s or s and p functions
      !
      IF (lmaxa.eq.0) THEN
        !
        ! atom b has d, the atom a only s functions
        !
        CALL sksd(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
      ELSE IF (lmaxa.eq.1) THEN
        !
        ! atom b has d, the atom a s and p functions
        !
        CALL sksd(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
        CALL skpd(skab2,mat,iptr(:,:,:,lmaxa,lmaxb),.TRUE.)
      ELSE
        !
        ! atom a has d functions
        !
        IF (lmaxb.eq.0) THEN
          !
          ! atom a has d, atom b has only s functions
          !
          CALL sksd(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
        ELSE
          !
          ! atom a has d, atom b has s and p functions
          !
          CALL sksd(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
          CALL skpd(skab1,mat,iptr(:,:,:,lmaxa,lmaxb),.FALSE.)
        END IF
      END IF
    END IF
    !
  CONTAINS
    !
    ! The subroutines to turn the matrix elements are taken as internal subroutines
    ! as it is beneficial to inline them.
    !
    ! They are both turning the matrix elements and placing them appropriately
    ! into the matrix block
    !
! *****************************************************************************
!> \brief   s-s interaction (no rotation necessary)
!> \param skpar ...
!> \param mat ...
!> \version 1.0
! *****************************************************************************
    SUBROUTINE skss(skpar,mat)
    REAL(dp), INTENT(in)                     :: skpar(:)
    REAL(dp), INTENT(inout)                  :: mat(:,:)

      mat(1,1) = mat(1,1) + skpar(1)
      !
    END SUBROUTINE skss

! *****************************************************************************
!> \brief  s-p interaction (simple rotation)
!> \param skpar ...
!> \param mat ...
!> \param ind ...
!> \param transposed ...
!> \version 1.0
! *****************************************************************************
    SUBROUTINE sksp(skpar,mat,ind,transposed)
    REAL(dp), INTENT(in)                     :: skpar(:)
    REAL(dp), INTENT(inout)                  :: mat(:,:)
    INTEGER, INTENT(in)                      :: ind(0:,0:,0:)
    LOGICAL, INTENT(in)                      :: transposed

    INTEGER                                  :: l
    REAL(dp)                                 :: skp

      skp = skpar(ind(1,0,0))
      IF (transposed) THEN
        DO l = 1,3
          mat(1,l+1) = mat(1,l+1) + rr(l)*skp
        END DO
      ELSE
        DO l = 1,3
          mat(l+1,1) = mat(l+1,1) - rr(l)*skp
        END DO
      END IF
      !
    END SUBROUTINE sksp

! *****************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
! *****************************************************************************
    SUBROUTINE skpp(skpar,mat,ind)
    REAL(dp), INTENT(in)                     :: skpar(:)
    REAL(dp), INTENT(inout)                  :: mat(:,:)
    INTEGER, INTENT(in)                      :: ind(0:,0:,0:)

    INTEGER                                  :: ii, ir, is, k, l
    REAL(dp)                                 :: epp(6), matel(6), skppp, skpps

      epp(1:3) = rr2(1:3)
      DO l = 1,3
        epp(l+3) = rr(l)*rr(l+1)
      END DO
      skppp = skpar(ind(1,1,1))
      skpps = skpar(ind(1,1,0))
      !
      DO l = 1,3
        matel(l) = epp(l)*skpps + (1._dp-epp(l))*skppp
      END DO
      DO l = 4,6
        matel(l) = epp(l)*(skpps - skppp)
      END DO
      !
      DO ir = 1,3
        DO is = 1,ir-1
          ii = ir - is
          k = 3*ii-(ii*(ii-1))/2+is
          mat(is+1,ir+1) = mat(is+1,ir+1) + matel(k)
          mat(ir+1,is+1) = mat(ir+1,is+1) + matel(k)
        END DO
        mat(ir+1,ir+1) = mat(ir+1,ir+1) + matel(ir)
      END DO
    END SUBROUTINE skpp

! *****************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
!> \param transposed ...
! *****************************************************************************
    SUBROUTINE sksd(skpar,mat,ind,transposed)
    REAL(dp), INTENT(in)                     :: skpar(:)
    REAL(dp), INTENT(inout)                  :: mat(:,:)
    INTEGER, INTENT(in)                      :: ind(0:,0:,0:)
    LOGICAL, INTENT(in)                      :: transposed

    INTEGER                                  :: l
    REAL(dp)                                 :: d4, d5, es(5), r3, sksds

      sksds = skpar(ind(2,0,0))
      r3 = SQRT(3._dp)
      d4 = rr2(3) - 0.5_dp*(rr2(1)+rr2(2))
      d5 = rr2(1) - rr2(2)
      !
      DO l = 1,3
        es(l) = r3*rr(l)*rr(l+1)
      END DO
      es(4) = 0.5_dp*r3*d5
      es(5) = d4
      !
      IF (transposed) THEN
        DO l = 1,5
          mat(1,l+4) = mat(1,l+4) + es(l)*sksds
        END DO
      ELSE
        DO l = 1,5
          mat(l+4,1) = mat(l+4,1) + es(l)*sksds
        END DO
      END IF
    END SUBROUTINE sksd

! *****************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
!> \param transposed ...
! *****************************************************************************
    SUBROUTINE skpd(skpar,mat,ind,transposed)
    REAL(dp), INTENT(in)                     :: skpar(:)
    REAL(dp), INTENT(inout)                  :: mat(:,:)
    INTEGER, INTENT(in)                      :: ind(0:,0:,0:)
    LOGICAL, INTENT(in)                      :: transposed

    INTEGER                                  :: ir, is, k, l, m
    REAL(dp)                                 :: d3, d4, d5, d6, dm(15), &
                                                epd(13,2), r3, sktmp

      r3 = SQRT(3.0_dp)
      d3 = rr2(1) + rr2(2)
      d4 = rr2(3) - 0.5_dp*d3
      d5 = rr2(1) - rr2(2)
      d6 = rr(1)*rr(2)*rr(3)
      DO l = 1,3
        epd(l,1) = r3*rr2(l)*rr(l+1)
        epd(l,2) = rr(l+1)*(1.0_dp-2._dp*rr2(l))
        epd(l+4,1) = r3*rr2(l)*rr(l+2)
        epd(l+4,2) = rr(l+2)*(1.0_dp-2*rr2(l))
        epd(l+7,1) = 0.5_dp*r3*rr(l)*d5
        epd(l+10,1) = rr(l)*d4
      END DO
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
      DO m = 1,2
        sktmp = skpar(ind(2,1,m-1))
        dm(1)=dm(1)+epd(1,m)*sktmp
        dm(2)=dm(2)+epd(6,m)*sktmp
        dm(3)=dm(3)+epd(4,m)*sktmp
        dm(5)=dm(5)+epd(2,m)*sktmp
        dm(6)=dm(6)+epd(7,m)*sktmp
        dm(7)=dm(7)+epd(5,m)*sktmp
        dm(9)=dm(9)+epd(3,m)*sktmp
        DO l = 8,13
          dm(l+2) = dm(l+2)+epd(l,m)*sktmp
        END DO
      END DO
      !
      dm(4) = dm(3)
      dm(8) = dm(3)
      !
      IF (transposed) THEN
        DO ir = 1,5
          DO is = 1,3
            k=3*(ir-1)+is
            mat(is+1,ir+4) = mat(is+1,ir+4) + dm(k)
          END DO
        END DO
      ELSE
        DO ir = 1,5
          DO is = 1,3
            k=3*(ir-1)+is
            mat(ir+4,is+1) = mat(ir+4,is+1) - dm(k)
          END DO
        END DO
      END IF
      !
    END SUBROUTINE skpd

! *****************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
! *****************************************************************************
    SUBROUTINE skdd(skpar,mat,ind)
    REAL(dp), INTENT(in)                     :: skpar(:)
    REAL(dp), INTENT(inout)                  :: mat(:,:)
    INTEGER, INTENT(in)                      :: ind(0:,0:,0:)

    INTEGER                                  :: ii, ir, is, k, l, m
    REAL(dp)                                 :: d3, d4, d5, dd(3), dm(15), &
                                                e(15,3), r3

      r3 = SQRT(3._dp)
      d3 = rr2(1) + rr2(2)
      d4 = rr2(3) - 0.5_dp*d3
      d5 = rr2(1) - rr2(2)
      DO l = 1,3
        e(l,1) = rr2(l)*rr2(l+1)
        e(l,2) = rr2(l) + rr2(l+1) - 4._dp*e(l,1)
        e(l,3) = rr2(l+2) + e(l,1)
        e(l,1) = 3._dp*e(l,1)
      END DO
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
      DO l = 1,2
        e(l+5,1) = 3._dp*rr2(l+1)*dd(l)
        e(l+5,2) = dd(l)*(1._dp-4._dp*rr2(l+1))
        e(l+5,3) = dd(l)*(rr2(l+1)-1._dp)
      END DO
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
      DO l = 1,15
        DO m = 1,3
          dm(l) = dm(l)+e(l,m)*skpar(ind(2,2,m-1))
        END DO
      END DO
      !
      DO ir = 1,5
        DO is = 1,ir-1
          ii = ir-is
          k = 5*ii-(ii*(ii-1))/2+is
          mat(ir+4,is+4) = mat(ir+4,is+4) + dm(k)
          mat(is+4,ir+4) = mat(is+4,ir+4) + dm(k)
        END DO
        mat(ir+4,ir+4) = mat(ir+4,ir+4) + dm(ir)
      END DO
    END SUBROUTINE skdd
    !
  END SUBROUTINE turnsk

! *****************************************************************************
!> \brief ...
!> \param xa ...
!> \param ya ...
!> \param n ...
!> \param x ...
!> \param y ...
!> \param dy ...
! *****************************************************************************
  SUBROUTINE polint(xa,ya,n,x,y,dy)
    INTEGER, INTENT(in)                      :: n
    REAL(dp), INTENT(in)                     :: ya(n), xa(n), x
    REAL(dp), INTENT(out)                    :: y, dy

    INTEGER                                  :: i, m, ns
    REAL(dp)                                 :: c(n), d(n), den, dif, dift, &
                                                ho, hp, w

!
!

    ns=1

    dif=ABS(x-xa(1))
    DO i = 1,n
      dift=ABS(x-xa(i))
      IF (dift.lt.dif) THEN
        ns=i
        dif=dift
      ENDIF
      c(i)=ya(i)
      d(i)=ya(i)
    END DO
    !
    y=ya(ns)
    ns=ns-1
    DO m = 1,n-1
      DO i = 1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        IF(den.eq.0.) STOP 'POLINT'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
      END DO
      IF (2*ns.lt.n-m)THEN
        dy=c(ns+1)
      ELSE
        dy=d(ns)
        ns=ns-1
      ENDIF
      y=y+dy
    END DO
!
    RETURN
  END SUBROUTINE polint

! *****************************************************************************
!> \brief ...
!> \param rv ...
!> \param r ...
!> \param erep ...
!> \param derep ...
!> \param n_urpoly ...
!> \param urep ...
!> \param spdim ...
!> \param s_cut ...
!> \param srep ...
!> \param spxr ...
!> \param scoeff ...
!> \param surr ...
!> \param dograd ...
! *****************************************************************************
  SUBROUTINE urep_egr(rv,r,erep,derep,&
                      n_urpoly,urep,spdim,s_cut,srep,spxr,scoeff,surr,dograd)

    REAL(dp), INTENT(in)                     :: rv(3), r
    REAL(dp), INTENT(inout)                  :: erep, derep(3)
    INTEGER, INTENT(in)                      :: n_urpoly
    REAL(dp), INTENT(in)                     :: urep(:)
    INTEGER, INTENT(in)                      :: spdim
    REAL(dp), INTENT(in)                     :: s_cut, srep(3)
    REAL(dp), POINTER                        :: spxr(:,:), scoeff(:,:)
    REAL(dp), INTENT(in)                     :: surr(2)
    LOGICAL, INTENT(in)                      :: dograd

    INTEGER                                  :: ic, isp, jsp, nsp
    REAL(dp)                                 :: de_z, rz

    derep=0._dp
    de_z = 0._dp
    IF (n_urpoly > 0) THEN
      !
      ! polynomial part
      !
      rz = urep(1) - r
      IF (rz <= rtiny) RETURN
      DO ic = 2,n_urpoly
        erep = erep + urep(ic) * rz**(ic)
      END DO
      IF (dograd) THEN
        DO ic = 2,n_urpoly
          de_z = de_z - ic*urep(ic) * rz**(ic-1)
        END DO
      END IF
    ELSE IF (spdim > 0) THEN
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
      IF (r > s_cut) RETURN  ! screening (condition a)
      !
      IF (r < spxr(1,1)) THEN
        ! a) short range
        erep = erep + EXP(-srep(1)*r + srep(2)) + srep(3)
        IF (dograd) de_z = de_z -srep(1)*EXP(-srep(1)*r + srep(2))
      ELSE
        !
        ! condition c). First determine between which places the spline is located:
        !
        ispg: DO isp = 1,spdim ! condition ca)
          IF (r <  spxr(isp,1)) CYCLE ispg ! distance is smaller than this spline range
          IF (r >= spxr(isp,2)) CYCLE ispg ! distance is larger than this spline range
          ! at this point we have found the correct spline interval
          rz = r - spxr(isp,1)
          IF (isp /= spdim) THEN
            nsp = 3 ! condition ca
            DO jsp = 0,nsp
              erep = erep + scoeff(isp,jsp+1)*rz**(jsp)
            END DO
            IF (dograd) THEN
              DO jsp = 1,nsp
                de_z = de_z + jsp*scoeff(isp,jsp+1)*rz**(jsp-1)
              END DO
            END IF
          ELSE
            nsp = 5 ! condition cb
            DO jsp = 0,nsp
              IF( jsp <= 3 ) THEN
                erep = erep + scoeff(isp,jsp+1)*rz**(jsp)
              ELSE
                erep = erep + surr(jsp-3)*rz**(jsp)
              ENDIF
            END DO
            IF (dograd) THEN
              DO jsp = 1,nsp
                IF( jsp <= 3 ) THEN
                  de_z = de_z + jsp*scoeff(isp,jsp+1)*rz**(jsp-1)
                ELSE
                  de_z = de_z + jsp*surr(jsp-3)*rz**(jsp-1)
                ENDIF
              END DO
            END IF
          END IF
          EXIT ispg
        END DO ispg
      END IF
    END IF
    !
    IF (dograd) THEN
       IF ( r > 1.e-12_dp ) derep(1:3) = (de_z/r)*rv(1:3)
    END IF

  END SUBROUTINE urep_egr

! *****************************************************************************
!> \brief  Computes the short-range gamma parameter from exact Coulomb
!>         interaction of normalized exp(-a*r) charge distribution - 1/r
!> \param r ...
!> \param ga ...
!> \param gb ...
!> \param hb_para ...
!> \retval gamma ...
!> \par Literature
!>         Elstner et al, PRB 58 (1998) 7260
!> \par History
!>      10.2008 Axel Kohlmeyer - adding sr_damp
!>      08.2014 JGH - adding flexible exponent for damping
!> \version 1.1
! *****************************************************************************
  FUNCTION gamma_rab_sr(r,ga,gb,hb_para) RESULT(gamma)
    REAL(dp), INTENT(in)                     :: r, ga, gb, hb_para
    REAL(dp)                                 :: gamma

    REAL(dp)                                 :: a, b, fac, g_sum

    gamma = 0.0_dp
    a = 3.2_dp*ga ! 3.2 = 16/5 in Eq. 18 and ff.
    b = 3.2_dp*gb
    g_sum = a + b
    IF (g_sum < tol_gamma) RETURN ! hardness screening
    IF (r < rtiny) THEN     ! This is for short distances but non-onsite terms
      ! This gives also correct diagonal elements (a=b, r=0)
      gamma = 0.5_dp*(a*b/g_sum + (a*b)**2/g_sum**3)
      RETURN
    END IF
    !
    ! distinguish two cases: Gamma's are very close, e.g. for the same atom type,
    !                        and Gamma's are different
    !
    IF (ABS(a-b) < rtiny) THEN
      fac = 1.6_dp*r*a*b/g_sum*(1.0_dp + a*b/g_sum**2)
      gamma = -(48.0_dp + 33._dp*fac + (9.0_dp + fac)*fac**2)*EXP(-fac)/(48._dp*r)
    ELSE
      gamma = -EXP(-a*r)*(0.5_dp*a*b**4/(a**2-b**2)**2 - &
                (b**6 - 3._dp*a**2*b**4)/(r*(a**2-b**2)**3)) - & ! a-> b
               EXP(-b*r)*(0.5_dp*b*a**4/(b**2-a**2)**2 - &
                (a**6 - 3._dp*b**2*a**4)/(r*(b**2-a**2)**3)) ! b-> a
    END IF
    !
    ! damping function for better short range hydrogen bonds.
    ! functional form from Hu H. et al., J. Phys. Chem. A 2007, 111, 5685-5691
    ! according to Elstner M, Theor. Chem. Acc. 2006, 116, 316-325,
    ! this should only be applied to a-b pairs involving hydrogen.
    IF (hb_para > 0.0_dp) THEN
      gamma = gamma * EXP(-(0.5_dp*(ga+gb))**hb_para *r*r)
    END IF
  END FUNCTION gamma_rab_sr

END MODULE qs_dftb_matrices

