!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the moment integrals <a|r^m|b> and <a|r x d/dr|b>
!> \par History
!>      added angular moments (JGH 11.2012)
!> \author JGH (20.07.2006)
! *****************************************************************************
MODULE qs_moments

  USE ai_angmom,                       ONLY: angmom
  USE ai_moments,                      ONLY: contract_cossin,&
                                             cossin,&
                                             moment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_get_info,&
                                             cp_cfm_p_type,&
                                             cp_cfm_release
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_trace
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_double,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             put_results
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE moments_utils,                   ONLY: get_reference_point
  USE orbital_pointers,                ONLY: current_maxl,&
                                             indco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: bohr,&
                                             debye
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_moments'

! *** Public subroutines ***

  PUBLIC :: build_berry_moment_matrix, build_local_moment_matrix
  PUBLIC :: qs_moment_berry_phase, qs_moment_locop

CONTAINS

! *****************************************************************************
  SUBROUTINE build_local_moment_matrix(qs_env,moments,nmoments,ref_point,ref_points,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: moments
    INTEGER, INTENT(IN)                      :: nmoments
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: ref_point
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: ref_points
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_moment_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, istat, &
      jatom, jkind, jset, last_jatom, maxco, maxsgf, natom, ncoa, ncob, &
      nkind, nm, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                npgfa, npgfb, nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: mab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc, rc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: mint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    IF ( nmoments < 1 ) RETURN

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atomic_kind_set,cell,particle_set,sab_orb)
    NULLIFY (matrix_s)

    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s,error=error)

    nm = (6 + 11*nmoments + 6*nmoments**2 + nmoments**3)/6 - 1
    CALL cp_dbcsr_allocate_matrix_set ( moments, nm ,error=error)

    DO i=1,nm
       ALLOCATE(moments(i)%matrix)
       CALL cp_dbcsr_init (moments(i)%matrix, error=error)
       CALL cp_dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments", error)
       CALL cp_dbcsr_set(moments(i)%matrix,0.0_dp,error=error)
    END DO

    NULLIFY (atomic_kind_set,particle_set,sab_orb,cell)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,cell=cell,&
                    sab_orb=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf)

    ALLOCATE (mab(maxco,maxco,nm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mab(:,:,:) = 0.0_dp

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    work(:,:) = 0.0_dp

    ALLOCATE (mint(nm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nm
      NULLIFY (mint(i)%block)
    END DO

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       IF(inode==1) last_jatom = 0

       ! this guarentees minimum image convention
       ! anything else would not make sense
       IF (jatom == last_jatom) THEN
         CYCLE
       END IF

       last_jatom = jatom

       IF (iatom <= jatom) THEN
         irow = iatom
         icol = jatom
       ELSE
         irow = jatom
         icol = iatom
       END IF

       DO i=1,nm
         NULLIFY (mint(i)%block)
         CALL cp_dbcsr_get_block_p(matrix=moments(i)%matrix,&
              row=irow,col=icol,BLOCK=mint(i)%block,found=found)
         mint(i)%block = 0._dp
       END DO

       ! fold atomic position back into unit cell
       IF ( PRESENT (ref_points) ) THEN
          rc(:) = 0.5_dp * (ref_points(:,iatom) + ref_points(:,jatom))
       ELSE IF ( PRESENT (ref_point) ) THEN
          rc(:) = ref_point(:)
       ELSE
          rc(:) = 0._dp
       END IF
       ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
       ! by folding around the center, such screwing can be avoided for a proper choice of center.
       ra(:) = pbc(particle_set(iatom)%r(:)-rc,cell)+rc
       rb(:) = pbc(particle_set(jatom)%r(:)-rc,cell)+rc
       ! we dont use PBC at this point
       rab(:) = ra(:) - rb(:)
       rac(:) = ra(:) - rc(:)
       rbc(:) = rb(:) - rc(:)
       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO jset=1,nsetb

           IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

           ncob = npgfb(jset)*ncoset(lb_max(jset))
           sgfb = first_sgfb(1,jset)

!          *** Calculate the primitive integrals ***

           CALL moment(la_max(iset),npgfa(iset),zeta(:,iset),&
                       rpgfa(:,iset),la_min(iset),&
                       lb_max(jset),npgfb(jset),zetb(:,jset),&
                       rpgfb(:,jset),nmoments,rac,rbc,mab)

!          *** Contraction step ***

           DO i=1,nm

             CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,mab(1,1,i),SIZE(mab,1),&
                        sphi_b(1,sgfb),SIZE(sphi_b,1),&
                        0.0_dp,work(1,1),SIZE(work,1))

             IF (iatom <= jatom) THEN

               CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          work(1,1),SIZE(work,1),&
                          1.0_dp,mint(i)%block(sgfa,sgfb),&
                          SIZE(mint(i)%block,1))

             ELSE

               CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                          1.0_dp,work(1,1),SIZE(work,1),&
                          sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          1.0_dp,mint(i)%block(sgfb,sgfa),&
                          SIZE(mint(i)%block,1))

             END IF

           END DO

         END DO
       END DO

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

!   *** Release work storage ***

    DEALLOCATE (mab,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nm
      NULLIFY (mint(i)%block)
    END DO
    DEALLOCATE (mint,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE build_local_moment_matrix

! *****************************************************************************
  SUBROUTINE build_local_magmom_matrix(qs_env,magmom,nmoments,ref_point,ref_points,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: magmom
    INTEGER, INTENT(IN)                      :: nmoments
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: ref_point
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: ref_points
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_magmom_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, istat, &
      jatom, jkind, jset, last_jatom, maxco, maxsgf, natom, ncoa, ncob, &
      nkind, nm, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                npgfa, npgfb, nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: mab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc, rc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: mint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    IF ( nmoments < 1 ) RETURN

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atomic_kind_set,cell,particle_set,sab_orb)
    NULLIFY (matrix_s)

    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s,error=error)

    ! magnetic dipoles/angular moments only
    nm = 3
    CALL cp_dbcsr_allocate_matrix_set ( magmom, nm ,error=error)

    DO i=1,nm
       ALLOCATE(magmom(i)%matrix)
       CALL cp_dbcsr_init (magmom(i)%matrix, error=error)
       CALL cp_dbcsr_copy(magmom(i)%matrix, matrix_s(1)%matrix, "Moments", error)
       CALL cp_dbcsr_set(magmom(i)%matrix,0.0_dp,error=error)
    END DO

    NULLIFY (atomic_kind_set,particle_set,sab_orb,cell)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,cell=cell,&
                    sab_orb=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf)

    ALLOCATE (mab(maxco,maxco,nm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mab(:,:,:) = 0.0_dp

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    work(:,:) = 0.0_dp

    ALLOCATE (mint(nm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nm
      NULLIFY (mint(i)%block)
    END DO

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       IF(inode==1) last_jatom = 0

       ! this guarentees minimum image convention
       ! anything else would not make sense
       IF (jatom == last_jatom) THEN
         CYCLE
       END IF

       last_jatom = jatom

       IF (iatom <= jatom) THEN
         irow = iatom
         icol = jatom
       ELSE
         irow = jatom
         icol = iatom
       END IF

       DO i=1,nm
         NULLIFY (mint(i)%block)
         CALL cp_dbcsr_get_block_p(matrix=magmom(i)%matrix,&
              row=irow,col=icol,BLOCK=mint(i)%block,found=found)
         mint(i)%block = 0._dp
       END DO

       ! fold atomic position back into unit cell
       IF ( PRESENT (ref_points) ) THEN
          rc(:) = 0.5_dp * (ref_points(:,iatom) + ref_points(:,jatom))
       ELSE IF ( PRESENT (ref_point) ) THEN
          rc(:) = ref_point(:)
       ELSE
          rc(:) = 0._dp
       END IF
       ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
       ! by folding around the center, such screwing can be avoided for a proper choice of center.
       ra(:) = pbc(particle_set(iatom)%r(:)-rc,cell)+rc
       rb(:) = pbc(particle_set(jatom)%r(:)-rc,cell)+rc
       ! we dont use PBC at this point
       rab(:) = ra(:) - rb(:)
       rac(:) = ra(:) - rc(:)
       rbc(:) = rb(:) - rc(:)
       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO jset=1,nsetb

           IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

           ncob = npgfb(jset)*ncoset(lb_max(jset))
           sgfb = first_sgfb(1,jset)

!          *** Calculate the primitive integrals ***

           CALL angmom(la_max(iset),npgfa(iset),zeta(:,iset),&
                       rpgfa(:,iset),la_min(iset),&
                       lb_max(jset),npgfb(jset),zetb(:,jset),&
                       rpgfb(:,jset),rac,rbc,mab)

!          *** Contraction step ***

           DO i=1,nm

             CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,mab(1,1,i),SIZE(mab,1),&
                        sphi_b(1,sgfb),SIZE(sphi_b,1),&
                        0.0_dp,work(1,1),SIZE(work,1))

             IF (iatom <= jatom) THEN

               CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          work(1,1),SIZE(work,1),&
                          1.0_dp,mint(i)%block(sgfa,sgfb),&
                          SIZE(mint(i)%block,1))

             ELSE

               CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                          1.0_dp,work(1,1),SIZE(work,1),&
                          sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          1.0_dp,mint(i)%block(sgfb,sgfa),&
                          SIZE(mint(i)%block,1))

             END IF

           END DO

         END DO
       END DO

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

!   *** Release work storage ***

    DEALLOCATE (mab,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nm
      NULLIFY (mint(i)%block)
    END DO
    DEALLOCATE (mint,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE build_local_magmom_matrix

! *****************************************************************************
  SUBROUTINE build_berry_moment_matrix(qs_env,cosmat,sinmat,kvec,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: kvec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_moment_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, istat, jatom, &
      jkind, jset, ldab, ldsa, ldsb, ldwork, natom, ncoa, ncob, nkind, nseta, &
      nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found
    REAL(dp), DIMENSION(:, :), POINTER       :: cblock, cosab, sblock, sinab, &
                                                work
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atomic_kind_set,particle_set,sab_orb,cell)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,cell=cell,&
                    sab_orb=sab_orb,error=error)

    CALL cp_dbcsr_set(sinmat,0.0_dp,error=error)
    CALL cp_dbcsr_set(cosmat,0.0_dp,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxco=ldwork )
    ldab = ldwork
    ALLOCATE(cosab(ldab,ldab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(sinab(ldab,ldab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(ldwork,ldwork),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       lb_min       =>  basis_set_b%lmin
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       ldsa = SIZE(sphi_a,1)
       ldsb = SIZE(sphi_b,1)

       IF (iatom <= jatom) THEN
         irow = iatom
         icol = jatom
       ELSE
         irow = jatom
         icol = iatom
       END IF

       NULLIFY (cblock)
       CALL cp_dbcsr_get_block_p(matrix=cosmat,&
            row=irow,col=icol,BLOCK=cblock,found=found)
       NULLIFY (sblock)
       CALL cp_dbcsr_get_block_p(matrix=sinmat,&
            row=irow,col=icol,BLOCK=sblock,found=found)
       IF(ASSOCIATED(cblock).AND..NOT.ASSOCIATED(sblock).OR.&
          .NOT.ASSOCIATED(cblock).AND.ASSOCIATED(sblock)) THEN
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(cblock).AND.ASSOCIATED(sblock)) THEN

          ra(:) = pbc(particle_set(iatom)%r(:),cell)
          rb(:) = ra+rab
          dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

          DO iset=1,nseta

             ncoa = npgfa(iset)*ncoset(la_max(iset))
             sgfa = first_sgfa(1,iset)

             DO jset=1,nsetb

                IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

!          *** Calculate the primitive integrals ***
                CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                     lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                     ra,rb,kvec,cosab,sinab)
                CALL contract_cossin(cblock,sblock,&
                     iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                     jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                     cosab,sinab,ldab,work,ldwork)

             END DO
          END DO

       ENDIF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE(cosab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(sinab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE build_berry_moment_matrix

! *****************************************************************************
  SUBROUTINE qs_moment_berry_phase (qs_env,magnetic,nmoments,reference,ref_point,unit_number,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: magnetic
    INTEGER, INTENT(IN)                      :: nmoments, reference
    REAL(dp), DIMENSION(:), POINTER          :: ref_point
    INTEGER, INTENT(IN)                      :: unit_number
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_phase', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=8), ALLOCATABLE, &
      DIMENSION(:)                           :: rlab
    CHARACTER(LEN=default_string_length)     :: description
    COMPLEX(dp)                              :: xphase(3), zdet, zdeta, &
                                                zi(3), zij(3,3), zijk(3,3,3), &
                                                zijkl(3,3,3,3), zphase(3), zz
    INTEGER                                  :: handle, i, ia, idim, ispin, &
                                                istat, ix, iy, iz, j, k, l, &
                                                nao, nm, nmo, nmom, nmotot, &
                                                tmp_dim
    LOGICAL                                  :: failure, ghost, uniform
    REAL(dp)                                 :: charge, ci(3), cij(3,3), dd, &
                                                occ, trace
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: mmom
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rmom
    REAL(dp), DIMENSION(3)                   :: kvec, rcc, ria
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: eigrmat
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, opvec
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: op_fm_set
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp

     failure = .FALSE.
     CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(qs_env%ls_scf_env)) THEN
      IF(unit_number > 0) WRITE(unit_number,*) "Periodic moment calculation not implemented in linear scaling code"
      RETURN
    ENDIF

     IF (.NOT.failure) THEN

       CALL timeset(routineN,handle)

       ! restrict maximum moment available
       nmom = MIN(nmoments,2)

       nm = (6 + 11*nmom + 6*nmom**2 + nmom**3)/6 - 1
       ! rmom(:,1)=electronic
       ! rmom(:,2)=nuclear
       ! rmom(:,1)=total
       ALLOCATE (rmom(nm+1,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (rlab(nm+1),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       rmom = 0.0_dp
       rlab = ""
       IF(magnetic) THEN
          nm=3
          ALLOCATE (mmom(nm),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          mmom = 0._dp
       END IF

       para_env => qs_env%para_env
       NULLIFY ( dft_control, rho, cell, particle_set )
       CALL get_qs_env ( qs_env, dft_control=dft_control, rho=rho, cell=cell,&
                         particle_set=particle_set, error=error)
       NULLIFY ( local_particles )
       CALL get_qs_env ( qs_env, local_particles=local_particles, error=error)
       NULLIFY ( matrix_s, mos )
       CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, mos=mos, error=error)

       NULLIFY ( cosmat, sinmat )
       ALLOCATE(cosmat, sinmat)
       CALL cp_dbcsr_init(cosmat, error=error)
       CALL cp_dbcsr_init(sinmat, error=error)
       CALL cp_dbcsr_copy(cosmat,matrix_s(1)%matrix,'COS MOM',error=error)
       CALL cp_dbcsr_copy(sinmat,matrix_s(1)%matrix,'SIN MOM',error=error)
       CALL cp_dbcsr_set(cosmat,0.0_dp,error=error)
       CALL cp_dbcsr_set(sinmat,0.0_dp,error=error)

       ALLOCATE ( op_fm_set( 2, dft_control%nspins ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( opvec( dft_control%nspins ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( eigrmat( dft_control%nspins ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nmotot = 0
       DO ispin = 1, dft_control%nspins
         NULLIFY(tmp_fm_struct,mo_coeff)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
         nmotot = nmotot + nmo
         CALL cp_fm_create (opvec(ispin)%matrix , mo_coeff%matrix_struct ,error=error)
         CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,&
              ncol_global=nmo,para_env=para_env,context=mo_coeff%matrix_struct%context,&
              error=error)
         DO i = 1, SIZE ( op_fm_set, 1 )
            NULLIFY(op_fm_set(i,ispin)%matrix)
            CALL cp_fm_create (op_fm_set(i,ispin)%matrix , tmp_fm_struct ,error=error)
         END DO
         CALL cp_cfm_create ( eigrmat(ispin)%matrix, op_fm_set(1,ispin)%matrix%matrix_struct ,&
         error=error)
         CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       END DO

       ! occupation
       DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=occ,uniform_occupation=uniform)
         IF (.NOT.uniform) THEN
            CALL cp_unimplemented_error(fromWhere=routineP, &
                 message="Berry phase moments for non uniform MOs' occupation numbers not implemented", &
                 error=error, error_level=cp_failure_level)
         END IF
       END DO

       ! reference point
       CALL get_reference_point(rcc,qs_env=qs_env,reference=reference,ref_point=ref_point,error=error)
       rcc = pbc(rcc,cell)

       ! label
       DO l = 1, nm
          ix = indco(1,l+1)
          iy = indco(2,l+1)
          iz = indco(3,l+1)
          CALL set_label (rlab(l+1),ix,iy,iz)
       END DO

       ! nuclear contribution
       DO ia = 1,SIZE(particle_set)
         atomic_kind => particle_set(ia)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=charge,ghost=ghost)
         IF (.NOT. ghost) THEN
            rmom(1,2) = rmom(1,2) - charge
         ENDIF
       END DO
       ria = twopi * MATMUL ( cell%h_inv, rcc )
       zphase = CMPLX ( COS(ria), SIN(ria), dp )**rmom(1,2)

       zi = 0._dp
       zij = 0._dp
       zijk = 0._dp
       zijkl = 0._dp

       DO l = 1, nmom
         SELECT CASE (l)
           CASE (1)
             ! Dipole
             zi(:) = CMPLX ( 1._dp, 0._dp, dp )
             DO ia = 1,SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=charge,ghost=ghost)
               IF (.NOT. ghost) THEN
                  ria = particle_set(ia)%r
                  ria = pbc(ria,cell)
                  DO i = 1, 3
                    kvec(:) = twopi*cell%h_inv(i,:)
                    dd = SUM ( kvec(:) * ria(:) )
                    zdeta = CMPLX(COS(dd),SIN(dd),KIND=dp)**charge
                    zi(i) = zi(i) * zdeta
                  END DO
               ENDIF
             END DO
             zi = zi * zphase
             ci = AIMAG(LOG(zi)) / twopi
             rmom(2:4,2) = MATMUL ( cell%hmat, ci )
           CASE (2)
             ! Quadrupole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 1 not implemented", &
                  error=error, error_level=cp_failure_level)
             zij(:,:) = CMPLX ( 1._dp, 0._dp, dp )
             DO ia = 1,SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=charge)
               ria = particle_set(ia)%r
               ria = pbc(ria,cell)
               DO i = 1, 3
                 DO j = i, 3
                   kvec(:) = twopi*(cell%h_inv(i,:)+cell%h_inv(j,:))
                   dd = SUM ( kvec(:) * ria(:) )
                   zdeta = CMPLX(COS(dd),SIN(dd),KIND=dp)**charge
                   zij(i,j) = zij(i,j) * zdeta
                 END DO
                 zij(j,i) = zij(i,j)
               END DO
             END DO
             DO i=1,3
               DO j=1,3
                 zij(i,j) = zij(i,j)*zphase(i)*zphase(j)
                 zz = zij(i,j)/zi(i)/zi(j)
                 cij(i,j) = AIMAG(LOG(zz))/twopi
               END DO
             END DO
             cij = 0.5_dp*cij/twopi/twopi
             cij = MATMUL ( MATMUL ( cell%hmat, cij ), TRANSPOSE(cell%hmat) )
             DO k = 4, 9
               ix = indco(1,k+1)
               iy = indco(2,k+1)
               iz = indco(3,k+1)
               IF ( ix == 0 ) THEN
                 rmom(k+1,2) = cij(iy,iz)
               ELSE IF ( iy == 0 ) THEN
                 rmom(k+1,2) = cij(ix,iz)
               ELSE IF ( iz == 0 ) THEN
                 rmom(k+1,2) = cij(ix,iy)
               END IF
             END DO
           CASE (3)
             ! Octapole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 2 not implemented", &
                  error=error, error_level=cp_failure_level)
           CASE (4)
             ! Hexadecapole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 3 not implemented", &
                  error=error, error_level=cp_failure_level)
           CASE DEFAULT
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 4 not implemented", &
                  error=error, error_level=cp_failure_level)
         END SELECT
       END DO

       ! electronic contribution

       ria = twopi * REAL(nmotot,dp) * occ * MATMUL ( cell%h_inv, rcc )
       xphase = CMPLX ( COS(ria), SIN(ria), dp )

       ! charge
       trace = 0.0_dp
       DO ispin=1,dft_control%nspins
         CALL cp_dbcsr_trace(rho%rho_ao(ispin)%matrix,matrix_s(1)%matrix,trace,error=error)
         rmom(1,1) = rmom(1,1) + trace
       END DO

       zi = 0._dp
       zij = 0._dp
       zijk = 0._dp
       zijkl = 0._dp

       DO l = 1, nmom
         SELECT CASE (l)
           CASE (1)
             ! Dipole
             DO i = 1, 3
                kvec(:) = twopi*cell%h_inv(i,:)
                CALL build_berry_moment_matrix(qs_env,cosmat,sinmat,kvec,error)
                IF(qs_env%run_rtp)THEN
                   CALL get_qs_env(qs_env,rtp=rtp,error=error)
                   CALL get_rtp(rtp,mos_new=mos_new,error=error)
                   CALL op_orbbas_rtp (cosmat,sinmat,mos,op_fm_set,opvec,mos_new,error)
                ELSE
                   CALL op_orbbas (cosmat,sinmat,mos,op_fm_set,opvec,error)
                END IF
                zdet = CMPLX ( 1._dp, 0._dp, dp )
                DO ispin = 1, dft_control%nspins
                  CALL cp_cfm_get_info(eigrmat(ispin)%matrix,ncol_local=tmp_dim,error=error)
                  DO idim=1,tmp_dim
                     eigrmat(ispin)%matrix%local_data(:,idim) = &
                          CMPLX (op_fm_set(1,ispin)%matrix%local_data(:,idim), &
                                 -op_fm_set(2,ispin)%matrix%local_data(:,idim),dp)
                  END DO
                  CALL cp_cfm_lu_decompose ( eigrmat(ispin)%matrix, zdeta ,error)
                  zdet = zdet * zdeta
                  IF (dft_control%nspins==1) THEN
                     zdet = zdet * zdeta
                  ENDIF
                END DO
                zi(i) = zdet
             END DO
             zi = zi * xphase
           CASE (2)
             ! Quadrupole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Berry phase moments bigger than 1 not implemented", &
                error=error, error_level=cp_failure_level)
             DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i,:)+cell%h_inv(j,:))
                  CALL build_berry_moment_matrix(qs_env,cosmat,sinmat,kvec,error)
                  IF(qs_env%run_rtp)THEN
                     CALL get_qs_env(qs_env,rtp=rtp,error=error)
                     CALL get_rtp(rtp,mos_new=mos_new,error=error)
                     CALL op_orbbas_rtp (cosmat,sinmat,mos,op_fm_set,opvec,mos_new,error)
                  ELSE
                     CALL op_orbbas (cosmat,sinmat,mos,op_fm_set,opvec,error)
                  END IF
                  zdet = CMPLX ( 1._dp, 0._dp, dp )
                  DO ispin = 1, dft_control%nspins
                    CALL cp_cfm_get_info(eigrmat(ispin)%matrix,ncol_local=tmp_dim,error=error)
                    DO idim=1,tmp_dim
                    eigrmat(ispin)%matrix%local_data(:,idim) = &
                          CMPLX (op_fm_set(1,ispin)%matrix%local_data(:,idim), &
                                 -op_fm_set(2,ispin)%matrix%local_data(:,idim),dp)
                    END DO
                    CALL cp_cfm_lu_decompose ( eigrmat(ispin)%matrix, zdeta ,error)
                    zdet = zdet * zdeta
                    IF (dft_control%nspins==1) THEN
                       zdet = zdet * zdeta
                    ENDIF
                  END DO
                  zij(i,j) = zdet*xphase(i)*xphase(j)
                  zij(j,i) = zdet*xphase(i)*xphase(j)
               END DO
             END DO
           CASE (3)
             ! Octapole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 2 not implemented", &
                  error=error, error_level=cp_failure_level)
           CASE (4)
             ! Hexadecapole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 3 not implemented", &
                  error=error, error_level=cp_failure_level)
           CASE DEFAULT
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 4 not implemented", &
                  error=error, error_level=cp_failure_level)
         END SELECT
       END DO
       DO l = 1, nmom
         SELECT CASE (l)
           CASE (1)
             ! Dipole
             ci = AIMAG(LOG(zi))/twopi
             rmom(2:4,1) = MATMUL ( cell%hmat, ci )
           CASE (2)
             ! Quadrupole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Berry phase moments bigger than 1 not implemented", &
                error=error, error_level=cp_failure_level)
             DO i=1,3
               DO j=1,3
                 zz = zij(i,j)/zi(i)/zi(j)
                 cij(i,j) =AIMAG(LOG(zz))/twopi
               END DO
             END DO
             cij = 0.5_dp*cij/twopi/twopi
             cij = MATMUL ( MATMUL ( cell%hmat, cij ), TRANSPOSE(cell%hmat) )
             DO k = 4, 9
               ix = indco(1,k+1)
               iy = indco(2,k+1)
               iz = indco(3,k+1)
               IF ( ix == 0 ) THEN
                 rmom(k+1,1) = cij(iy,iz)
               ELSE IF ( iy == 0 ) THEN
                 rmom(k+1,1) = cij(ix,iz)
               ELSE IF ( iz == 0 ) THEN
                 rmom(k+1,1) = cij(ix,iy)
               END IF
             END DO
           CASE (3)
             ! Octapole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 2 not implemented", &
                  error=error, error_level=cp_failure_level)
           CASE (4)
             ! Hexadecapole
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 3 not implemented", &
                  error=error, error_level=cp_failure_level)
           CASE DEFAULT
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Berry phase moments bigger than 4 not implemented", &
                  error=error, error_level=cp_failure_level)
         END SELECT
       END DO

       rmom(:,3) = rmom(:,1) + rmom(:,2)
       description="[DIPOLE]"
       CALL cp_results_erase(results=qs_env%results,description=description,error=error)
       CALL put_results(results=qs_env%results,description=description,&
                        values=rmom(2:4,3),error=error)
       IF(magnetic) THEN
          CALL print_moments(unit_number,nmom,rmom,rlab,rcc,cell,periodic=.TRUE.,mmom=mmom)
       ELSE
          CALL print_moments(unit_number,nmom,rmom,rlab,rcc,cell,periodic=.TRUE.)
       END IF

       DEALLOCATE (rmom,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (rlab,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF(magnetic) THEN
          DEALLOCATE ( mmom, STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       CALL cp_dbcsr_deallocate_matrix ( cosmat, error )
       CALL cp_dbcsr_deallocate_matrix ( sinmat, error )

       DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix,error=error)
         CALL cp_cfm_release(eigrmat(ispin)%matrix,error=error)
         DO i = 1, SIZE ( op_fm_set, 1 )
           CALL cp_fm_release(op_fm_set(i,ispin)%matrix,error=error)
         END DO
       END DO
       DEALLOCATE ( op_fm_set, STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( opvec, STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( eigrmat, STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       CALL timestop(handle)

     END IF

  END SUBROUTINE qs_moment_berry_phase

! *****************************************************************************
  SUBROUTINE op_orbbas (cosmat,sinmat,mos,op_fm_set,opvec,error)

    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: op_fm_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: opvec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: i, nao, nmo
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    DO i=1,SIZE ( op_fm_set, 2 )   ! spin
       CALL get_mo_set(mo_set=mos(i)%mo_set,nao=nao,mo_coeff=mo_coeff,nmo=nmo)
       CALL cp_dbcsr_sm_fm_multiply(cosmat, mo_coeff, opvec(i)%matrix, ncol=nmo, error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff,opvec(i)%matrix,0.0_dp,&
            op_fm_set(1,i)%matrix,error=error)
       CALL cp_dbcsr_sm_fm_multiply(sinmat, mo_coeff, opvec(i)%matrix, ncol=nmo, error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff,opvec(i)%matrix,0.0_dp,&
            op_fm_set(2,i)%matrix,error=error)
    ENDDO

  END SUBROUTINE op_orbbas

! *****************************************************************************
  SUBROUTINE op_orbbas_rtp (cosmat,sinmat,mos,op_fm_set,opvec,mos_new,error)

    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: op_fm_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: opvec, mos_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: i, icol, lcol, nao, newdim, &
                                                nmo
    LOGICAL                                  :: double_col, double_row
    TYPE(cp_fm_struct_type), POINTER         :: newstruct, newstruct1
    TYPE(cp_fm_type), POINTER                :: mo_coeff, work, work1, work2

    DO i=1,SIZE ( op_fm_set, 2 )   ! spin
       CALL get_mo_set(mo_set=mos(i)%mo_set,nao=nao,mo_coeff=mo_coeff,nmo=nmo)
       CALL cp_fm_get_info(mos_new(2*i)%matrix,ncol_local=lcol,ncol_global=nmo,error=error)
       double_col=.TRUE.
       double_row=.FALSE.
       CALL cp_fm_struct_double(newstruct,&
                                mos_new(2*i)%matrix%matrix_struct,&
                                mos_new(2*i)%matrix%matrix_struct%context,&
                                double_col,&
                                double_row,&
                                error)

       CALL cp_fm_create(work,matrix_struct=newstruct,error=error)
       CALL cp_fm_create(work1,matrix_struct=newstruct,error=error)
       CALL cp_fm_create(work2,matrix_struct=newstruct,error=error)
       CALL cp_fm_get_info(work,ncol_global=newdim,error=error)

       CALL cp_fm_set_all(work,0.0_dp,0.0_dp,error)
       DO icol=1,lcol
          work%local_data(:,icol)=mos_new(2*i-1)%matrix%local_data(:,icol)
          work%local_data(:,icol+lcol)=mos_new(2*i)%matrix%local_data(:,icol)
       END DO

       CALL cp_dbcsr_sm_fm_multiply(cosmat, work, work1, ncol=newdim, error=error)
       CALL cp_dbcsr_sm_fm_multiply(sinmat, work, work2, ncol=newdim, error=error)

       DO icol=1,lcol
          work%local_data(:,icol)=work1%local_data(:,icol)-work2%local_data(:,icol+lcol)
          work%local_data(:,icol+lcol)=work1%local_data(:,icol+lcol)+work2%local_data(:,icol)
       END DO

       CALL cp_fm_release(work1,error)
       CALL cp_fm_release(work2,error)

       CALL cp_fm_struct_double(newstruct1,&
                                op_fm_set(1,i)%matrix%matrix_struct,&
                                op_fm_set(1,i)%matrix%matrix_struct%context,&
                                double_col,&
                                double_row,&
                                error)

       CALL cp_fm_create(work1,matrix_struct=newstruct1,error=error)

       CALL cp_fm_gemm("T","N",nmo,newdim,nao,1.0_dp,mos_new(2*i-1)%matrix,&
            work,0.0_dp,work1,error=error)

       DO icol=1,lcol
          op_fm_set(1,i)%matrix%local_data(:,icol)=work1%local_data(:,icol)
          op_fm_set(2,i)%matrix%local_data(:,icol)=work1%local_data(:,icol+lcol)
       END DO


       CALL cp_fm_gemm("T","N",nmo,newdim,nao,1.0_dp,mos_new(2*i)%matrix,&
            work,0.0_dp,work1,error=error)

       DO icol=1,lcol
          op_fm_set(1,i)%matrix%local_data(:,icol)=&
               op_fm_set(1,i)%matrix%local_data(:,icol)+work1%local_data(:,icol+lcol)
          op_fm_set(2,i)%matrix%local_data(:,icol)=&
               op_fm_set(2,i)%matrix%local_data(:,icol)-work1%local_data(:,icol)
       END DO

       CALL cp_fm_release(work,error)
       CALL cp_fm_release(work1,error)
       CALL cp_fm_struct_release(newstruct,error=error)
       CALL cp_fm_struct_release(newstruct1,error=error)

    ENDDO

  END SUBROUTINE op_orbbas_rtp

! *****************************************************************************
  SUBROUTINE qs_moment_locop (qs_env,magnetic,nmoments,reference,ref_point,unit_number,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: magnetic
    INTEGER, INTENT(IN)                      :: nmoments, reference
    REAL(dp), DIMENSION(:), POINTER          :: ref_point
    INTEGER, INTENT(IN)                      :: unit_number
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_locop', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=8), ALLOCATABLE, &
      DIMENSION(:)                           :: rlab
    CHARACTER(LEN=default_string_length)     :: description
    INTEGER                                  :: handle, i, ia, iatom, ikind, &
                                                ispin, istat, ix, iy, iz, l, &
                                                nm, nmom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: charge, dd, strace, trace
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: mmom
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rmom
    REAL(dp), DIMENSION(3)                   :: rcc, ria
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: magmom, matrix_s, moments
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_rho_type), POINTER               :: rho

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

       CALL timeset(routineN,handle)

       ! reference point
       CALL get_reference_point(rcc,qs_env=qs_env,reference=reference,ref_point=ref_point,error=error)

       ! only allow for moments up to maxl set by basis
       nmom = MIN (nmoments,current_maxl)
       ! electronic contribution
       NULLIFY ( moments )
       CALL build_local_moment_matrix(qs_env,moments,nmom,ref_point=rcc,error=error)
       NULLIFY ( dft_control, rho, cell, particle_set )
       CALL get_qs_env ( qs_env, dft_control=dft_control, rho=rho, cell=cell,&
                         particle_set=particle_set, error=error)
       para_env => qs_env%para_env
       NULLIFY ( matrix_s )
       CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s,error=error)

       nm = SIZE(moments)
       ALLOCATE (rmom(nm+1,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (rlab(nm+1),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       rmom=0.0_dp
       rlab = ""

       trace = 0.0_dp
       DO ispin=1,dft_control%nspins
         CALL cp_dbcsr_trace(rho%rho_ao(ispin)%matrix,matrix_s(1)%matrix,trace,error=error)
         rmom(1,1) = rmom(1,1) + trace
       END DO

       DO i = 1, SIZE(moments)
          strace = 0._dp
          DO ispin=1,dft_control%nspins
            CALL cp_dbcsr_trace(rho%rho_ao(ispin)%matrix,moments(i)%matrix,trace,error=error)
            strace = strace + trace
          END DO
          rmom(i+1,1) = strace
       END DO

       CALL cp_dbcsr_deallocate_matrix_set ( moments, error )

       ! magnetic moments
       IF (magnetic) THEN
          NULLIFY ( magmom )
          CALL build_local_magmom_matrix(qs_env,magmom,nmom,ref_point=rcc,error=error)
          nm = SIZE(magmom)
          ALLOCATE (mmom(nm),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, SIZE(magmom)
             strace = 0._dp
             DO ispin=1,dft_control%nspins
               CALL cp_dbcsr_trace(rho%rho_ao(ispin)%matrix,magmom(i)%matrix,trace,error=error)
               strace = strace + trace
             END DO
             mmom(i) = strace
          END DO
       END IF

       ! nuclear contribution
       CALL get_qs_env(qs_env=qs_env,&
                      local_particles=local_particles,error=error)
       DO ikind = 1,SIZE(local_particles%n_el)
          DO ia = 1,local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            ! fold atomic positions back into unit cell
            ria = pbc(particle_set(iatom)%r-rcc,cell)+rcc
            ria = ria - rcc
            atomic_kind => particle_set(iatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=charge)
            rmom(1,2) = rmom(1,2) - charge
            DO l = 1, nm
              ix = indco(1,l+1)
              iy = indco(2,l+1)
              iz = indco(3,l+1)
              dd = 1._dp
              IF (ix > 0) dd = dd * ria(1)**ix
              IF (iy > 0) dd = dd * ria(2)**iy
              IF (iz > 0) dd = dd * ria(3)**iz
              rmom(l+1,2) = rmom(l+1,2) - charge * dd
              CALL set_label (rlab(l+1),ix,iy,iz)
            END DO
          END DO
       END DO
       CALL mp_sum(rmom(:,2),para_env%group)
       rmom(:,:) = -rmom(:,:)
       rmom(:,3) = rmom(:,1) + rmom(:,2)

       description="[DIPOLE]"
       CALL cp_results_erase(results=qs_env%results,description=description,error=error)
       CALL put_results(results=qs_env%results,description=description,&
                        values=rmom(2:4,3),error=error)
       IF(magnetic) THEN
          CALL print_moments(unit_number,nmom,rmom,rlab,rcc,cell,periodic=.FALSE.,mmom=mmom)
       ELSE
          CALL print_moments(unit_number,nmom,rmom,rlab,rcc,cell,periodic=.FALSE.)
       END IF

       DEALLOCATE (rmom,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (rlab,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF (magnetic) THEN
          DEALLOCATE (mmom,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       CALL timestop(handle)

     END IF

  END SUBROUTINE qs_moment_locop

! *****************************************************************************
  SUBROUTINE set_label(label,ix,iy,iz)
    CHARACTER(LEN=*), INTENT(OUT)            :: label
    INTEGER, INTENT(IN)                      :: ix, iy, iz

    INTEGER                                  :: i

    label = ""
    DO i=1,ix
      WRITE(label(i:),"(A1)") "X"
    END DO
    DO i=ix+1,ix+iy
      WRITE(label(i:),"(A1)") "Y"
    END DO
    DO i=ix+iy+1,ix+iy+iz
      WRITE(label(i:),"(A1)") "Z"
    END DO

  END SUBROUTINE set_label

! *****************************************************************************
  SUBROUTINE print_moments(unit_number,nmom,rmom,rlab,rcc,cell,periodic,mmom)
    INTEGER, INTENT(IN)                      :: unit_number, nmom
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: rmom
    CHARACTER(LEN=8), DIMENSION(:)           :: rlab
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rcc
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL                                  :: periodic
    REAL(dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: mmom

    INTEGER                                  :: i, i0, i1, j, l
    REAL(dp)                                 :: dd

    IF (unit_number > 0) THEN
      DO l=0,nmom
        SELECT CASE (l)
          CASE (0)
            WRITE(unit_number,"(T3,A,T32,3F14.8)") "Reference Point [Bohr]",rcc
            WRITE(unit_number,"(T3,A)") "Charges"
            WRITE(unit_number,"(T5,A,T18,F14.8,T36,A,T42,F14.8,T60,A,T68,F14.8)") &
             "Electronic=",rmom(1,1),"Core=",rmom(1,2),"Total=",rmom(1,3)
          CASE (1)
            IF (periodic) THEN
                WRITE(unit_number,"(T3,A)") "Dipole vectors are based on the periodic (Berry phase) operator."
                WRITE(unit_number,"(T3,A)") "They are defined modulo integer multiples of the cell matrix [Debye]."
                WRITE(unit_number,"(T3,A,3(F14.8,1X),A)") "[X] [",cell%hmat(1,:)*debye,"] [i]"
                WRITE(unit_number,"(T3,A,3(F14.8,1X),A)") "[Y]=[",cell%hmat(2,:)*debye,"]*[j]"
                WRITE(unit_number,"(T3,A,3(F14.8,1X),A)") "[Z] [",cell%hmat(3,:)*debye,"] [k]"
            ELSE
                WRITE(unit_number,"(T3,A)") "Dipoles are based on the traditional operator."
            ENDIF
            dd=SQRT(SUM(rmom(2:4,3)**2))*debye
            WRITE(unit_number,"(T3,A)") "Dipole moment [Debye]"
            WRITE(unit_number,"(T5,3(A,A,F14.8,1X),T60,A,T68,F14.8)") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye,i=2,4),"Total=",dd
          CASE (2)
            WRITE(unit_number,"(T3,A)") "Quadrupole moment [Debye*Angstrom]"
            WRITE(unit_number,"(T17,3(A,A,F14.8,9X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr,i=5,7)
            WRITE(unit_number,"(T17,3(A,A,F14.8,9X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr,i=8,10)
          CASE (3)
            WRITE(unit_number,"(T3,A)") "Octapole moment [Debye*Angstrom**2]"
            WRITE(unit_number,"(T7,4(A,A,F14.8,3X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr,i=11,14)
            WRITE(unit_number,"(T7,4(A,A,F14.8,3X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr,i=15,18)
            WRITE(unit_number,"(T7,4(A,A,F14.8,3X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr,i=19,20)
          CASE (4)
            WRITE(unit_number,"(T3,A)") "Hexadecapole moment [Debye*Angstrom**3]"
            WRITE(unit_number,"(T6,4(A,A,F14.8,2X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr/bohr,i=21,24)
            WRITE(unit_number,"(T6,4(A,A,F14.8,2X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr/bohr,i=25,28)
            WRITE(unit_number,"(T6,4(A,A,F14.8,2X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr/bohr,i=29,32)
            WRITE(unit_number,"(T6,4(A,A,F14.8,2X))") &
             (TRIM(rlab(i)),"=",rmom(i,3)*debye/bohr/bohr/bohr,i=32,35)
          CASE DEFAULT
            WRITE(unit_number,"(T3,A,A,I2)") "Higher moment [Debye*Angstrom**(L-1)]",&
               "  L=",l
            i0 = (6 + 11*(l-1) + 6*(l-1)**2 + (l-1)**3)/6
            i1 = (6 + 11*l + 6*l**2 + l**3)/6 - 1
            dd = debye/(bohr)**(l-1)
            DO i=i0,i1,3
              WRITE(unit_number,"(T18,3(A,A,F14.8,4X))") &
                (TRIM(rlab(j+1)),"=",rmom(j+1,3)*dd,j=i,MIN(i1,i+2))
            END DO
        END SELECT
      END DO
      IF(PRESENT(mmom)) THEN
        IF(nmom>=1) THEN
          dd=SQRT(SUM(mmom(1:3)**2))
          WRITE(unit_number,"(T3,A)") "Magnetic Dipole Moment [WRONG]"
          WRITE(unit_number,"(T5,3(A,A,F14.8,1X),T60,A,T68,F14.8)") &
             (TRIM(rlab(i+1)),"=",mmom(i),i=1,3),"Total=",dd
        END IF
      END IF
    END IF

  END SUBROUTINE print_moments

END MODULE qs_moments

