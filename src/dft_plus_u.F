!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief   Add the DFT+U contribution to the Hamiltonian matrix
!> \details The implemented methods refers to:\n
!>          S. L. Dudarev, D. Nguyen Manh, and A. P. Sutton,
!>          Philos. Mag. B \b 75, 613 (1997)\n
!>          S. L. Dudarev et al.,
!>          Phys. Rev. B \b 57, 1505 (1998)
!> \author  Matthias Krack (MK)
!> \date    14.01.2008
!> \version 1.0
! *****************************************************************************
MODULE dft_plus_u
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE bibliography,                    ONLY: Dudarev1997,&
                                             Dudarev1998,&
                                             cite_reference
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_get_block_diag, &
       cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
       cp_dbcsr_set, cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix,&
                                             write_fm_with_basis_info
  USE cp_fm_basic_linalg,              ONLY: cp_fm_transpose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             low_print_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: plus_u_lowdin,&
                                             plus_u_mulliken,&
                                             plus_u_mulliken_charges
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: jacobi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_symbols,                 ONLY: sgf_symbol
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type,&
                                             set_qs_kind
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dft_plus_u'

  PUBLIC :: plus_u

CONTAINS
! *****************************************************************************
!> \brief         Add the DFT+U contribution to the Hamiltonian matrix.\n
!>                Wrapper routine for all "+U" methods
!> \param[in]     qs_env Quickstep environment
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param[in,out] matrix_w Energy weighted density matrices for each spin
!> \date          14.01.2008
!> \author        Matthias Krack (MK)
!> \version       1.0
! *****************************************************************************
  SUBROUTINE plus_u(qs_env,matrix_h,matrix_w)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_w

    CHARACTER(LEN=*), PARAMETER :: routineN = 'plus_u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit, &
                                                print_level
    LOGICAL                                  :: failure, orthonormal_basis, &
                                                should_output
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPASSERT(ASSOCIATED(qs_env))

    NULLIFY (input, dft_control)

    logger => cp_get_default_logger()

    CALL get_qs_env(qs_env=qs_env,&
                    input=input,&
                    dft_control=dft_control)

    CALL cite_reference(Dudarev1997)
    CALL cite_reference(Dudarev1998)

    ! Later we could save here some time, if the method in use has this property
    ! which then has to be figured out here

    orthonormal_basis = .FALSE.

    ! Setup print control

    print_level = logger%iter_info%print_level
    should_output = (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                          "DFT%PRINT%PLUS_U"),cp_p_file).AND.&
                     (.NOT.PRESENT(matrix_w)))
    output_unit = cp_print_key_unit_nr(logger,input,"DFT%PRINT%PLUS_U",&
                                       extension=".plus_u",&
                                       ignore_should_output=should_output,&
                                       log_filename=.FALSE.)

    ! Select DFT+U method

    SELECT CASE (dft_control%plus_u_method_id)
    CASE (plus_u_lowdin)
      IF (orthonormal_basis) THEN
        ! For an orthonormal basis the Lowdin method and the Mulliken method
        ! are equivalent
        CALL mulliken(qs_env,orthonormal_basis,matrix_h,&
                      should_output,output_unit,print_level)
      ELSE
        CALL lowdin(qs_env,matrix_h,matrix_w,&
                    should_output,output_unit,print_level)
      END IF
    CASE (plus_u_mulliken)
      CALL mulliken(qs_env,orthonormal_basis,matrix_h,&
                    should_output,output_unit,print_level)
    CASE (plus_u_mulliken_charges)
      CALL mulliken_charges(qs_env,orthonormal_basis,matrix_h,matrix_w,&
                            should_output,output_unit,print_level)
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Invalid DFT+U method requested",only_ionode=.TRUE.)
    END SELECT

    CALL cp_print_key_finished_output(output_unit,logger,input,"DFT%PRINT%PLUS_U",&
                                      ignore_should_output=should_output)

    CALL timestop(handle)

  END SUBROUTINE plus_u

! *****************************************************************************
!> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
!>                using a method based on Lowdin charges
!>                \f[Q = S^{1/2} P S^{1/2}\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in]     qs_env Quickstep environment
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param[in,out] matrix_w Energy weighted density matrices for each spin
!> \param should_output ...
!> \param output_unit ...
!> \param print_level ...
!> \date          02.07.2008
!> \par
!>  \f{eqnarray*}{
!>   E^{\rm DFT+U} & = & E^{\rm DFT} +  E^{\rm U}\\\
!>                 & = & E^{\rm DFT} +  \frac{1}{2}(U - J)\sum_\mu (q_\mu - q_\mu^2)\\[1ex]
!>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
!>                          & = & \frac{\partial E^{\rm DFT}}
!>                                     {\partial P_{\mu\nu}} +
!>                                \frac{\partial E^{\rm U}}
!>                                     {\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} +
!>                                \frac{\partial E^{\rm U}}{\partial q_\mu}
!>                                \frac{\partial q_\mu}{\partial P_{\mu\nu}}\\\
!>  \f}
!> \author        Matthias Krack (MK)
!> \version       1.0
! *****************************************************************************
  SUBROUTINE lowdin(qs_env,matrix_h,matrix_w,should_output,output_unit,&
                    print_level)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_w
    LOGICAL, INTENT(IN)                      :: should_output
    INTEGER, INTENT(IN)                      :: output_unit, print_level

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lowdin', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: spin_info
    CHARACTER(LEN=6), ALLOCATABLE, &
      DIMENSION(:)                           :: symbol
    CHARACTER(LEN=default_string_length)     :: atomic_kind_name
    INTEGER :: atom_a, handle, i, i0, iatom, ikind, iorb, isb, iset, isgf, &
      ishell, ispin, j, jsb, jset, jsgf, jshell, lu, m, max_scf, n, natom, &
      natom_of_kind, nkind, norb, nsb, nsbsize, nset, nsgf, nsgf_kind, nspin, &
      stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf_atom
    INTEGER, DIMENSION(1)                    :: iloc
    INTEGER, DIMENSION(:), POINTER           :: atom_list, nshell, orbitals
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgf, l, last_sgf
    LOGICAL :: debug, dft_plus_u_atom, failure, fm_work1_local_alloc, &
      fm_work2_local_alloc, found, just_energy, smear
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: orb_occ
    REAL(KIND=dp)                            :: eps_scf, eps_u_ramping, &
                                                fspin, occ, trq, trq2, &
                                                u_minus_j, u_minus_j_target, &
                                                u_ramping
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: q_eigval
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: q_eigvec, q_matrix, q_work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: q_block, v_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: sm_h, sm_p, sm_q, sm_s, sm_v, &
                                                sm_w
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    TYPE(cp_fm_type), POINTER                :: fm_s_half, fm_work1, fm_work2
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_scf_env_type), POINTER           :: scf_env

    CALL timeset(routineN,handle)

    debug = .FALSE. ! Set to .TRUE. to print debug information

    NULLIFY (atom_list)
    NULLIFY (atomic_kind_set)
    NULLIFY (qs_kind_set)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (first_sgf)
    NULLIFY (fm_s_half)
    NULLIFY (fm_work1)
    NULLIFY (fm_work2)
    NULLIFY (fmstruct)
    NULLIFY (matrix_p)
    NULLIFY (matrix_s)
    NULLIFY (l)
    NULLIFY (last_sgf)
    NULLIFY (nshell)
    NULLIFY (orb_basis_set)
    NULLIFY (orbitals)
    NULLIFY (particle_set)
    NULLIFY (q_block)
    NULLIFY (rho)
    NULLIFY (scf_env)
    NULLIFY (sm_h)
    NULLIFY (sm_p)
    NULLIFY (sm_q)
    NULLIFY (sm_s)
    NULLIFY (sm_v)
    NULLIFY (sm_w)
    NULLIFY (v_block)
    NULLIFY (para_env)
    NULLIFY (blacs_env)

    failure = .FALSE.

    smear = .FALSE.
    max_scf = -1
    eps_scf = 1.0E30_dp

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    dft_control=dft_control,&
                    energy=energy,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    rho=rho,&
                    scf_env=scf_env,&
                    para_env=para_env,&
                    blacs_env=blacs_env)

    CPASSERT(ASSOCIATED(atomic_kind_set))
    CPASSERT(ASSOCIATED(dft_control))
    CPASSERT(ASSOCIATED(energy))
    CPASSERT(ASSOCIATED(matrix_s))
    CPASSERT(ASSOCIATED(particle_set))
    CPASSERT(ASSOCIATED(rho))

    sm_s => matrix_s(1)%matrix ! Overlap matrix in sparse format
    CALL qs_rho_get(rho, rho_ao=matrix_p)     ! Density matrices in sparse format

    energy%dft_plus_u = 0.0_dp

    nspin = dft_control%nspins

    IF (nspin == 2) THEN
      fspin = 1.0_dp
    ELSE
      fspin = 0.5_dp
    END IF

    ! Get the total number of atoms, contracted spherical Gaussian basis
    ! functions, and atomic kinds

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (first_sgf_atom(natom),STAT=stat)
    CPASSERT(stat == 0)
    first_sgf_atom(:) = 0

    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf_atom)

    IF (PRESENT(matrix_h).OR.PRESENT(matrix_w)) THEN
      just_energy  = .FALSE.
    ELSE
      just_energy = .TRUE.
    END IF

    ! Retrieve S^(1/2) from the SCF environment

    fm_s_half => scf_env%s_half
    CPASSERT(ASSOCIATED(fm_s_half))

    ! Try to retrieve (full) work matrices from the SCF environment and reuse
    ! them, if available

    fm_work1_local_alloc = .FALSE.
    fm_work2_local_alloc = .FALSE.

    IF (ASSOCIATED(scf_env%scf_work1)) THEN
      IF (ASSOCIATED(scf_env%scf_work1(1)%matrix)) THEN
        fm_work1 => scf_env%scf_work1(1)%matrix
      END IF
    END IF

    fm_work2 => scf_env%scf_work2

    IF ((.NOT.ASSOCIATED(fm_work1)).OR.&
        (.NOT.ASSOCIATED(fm_work2))) THEN
      CALL cp_fm_struct_create(fmstruct=fmstruct,&
                               para_env=para_env,&
                               context=blacs_env,&
                               nrow_global=nsgf,&
                               ncol_global=nsgf)
      IF (.NOT.ASSOCIATED(fm_work1)) THEN
        CALL cp_fm_create(matrix=fm_work1,&
                          matrix_struct=fmstruct,&
                          name="FULL WORK MATRIX 1")
        fm_work1_local_alloc = .TRUE.
      END IF
      IF (.NOT.ASSOCIATED(fm_work2)) THEN
        CALL cp_fm_create(matrix=fm_work2,&
                          matrix_struct=fmstruct,&
                          name="FULL WORK MATRIX 2")
        fm_work2_local_alloc = .TRUE.
      END IF
      CALL cp_fm_struct_release(fmstruct=fmstruct)
    END IF

    ! Create local block diagonal matrices

    ALLOCATE (sm_q,STAT=stat)
    CPASSERT(stat == 0)
    CALL cp_dbcsr_init(sm_q)
    CALL cp_dbcsr_get_block_diag(sm_s,sm_q)

    ALLOCATE (sm_v,STAT=stat)
    CPASSERT(stat == 0)
    CALL cp_dbcsr_init(sm_v)
    CALL cp_dbcsr_get_block_diag(sm_s,sm_v)

    ! Loop over all spins

    DO ispin=1,nspin

      IF (PRESENT(matrix_h)) THEN
        ! Hamiltonian matrix for spin ispin in sparse format
        sm_h => matrix_h(ispin)%matrix
      ELSE
        NULLIFY (sm_h)
      END IF

      IF (PRESENT(matrix_w)) THEN
        ! Energy weighted density matrix for spin ispin in sparse format
        sm_w => matrix_w(ispin)%matrix
      ELSE
        NULLIFY (sm_w)
      END IF

      sm_p => matrix_p(ispin)%matrix ! Density matrix for spin ispin in sparse format

      CALL cp_dbcsr_set(sm_q,0.0_dp)
      CALL cp_dbcsr_set(sm_v,0.0_dp)

      ! Calculate S^(1/2)*P*S^(1/2) as a full matrix (Lowdin)

      CALL cp_dbcsr_sm_fm_multiply(sm_p,fm_s_half,fm_work1,nsgf)
      CALL cp_gemm(transa="N",&
                      transb="N",&
                      m=nsgf,&
                      n=nsgf,&
                      k=nsgf,&
                      alpha=1.0_dp,&
                      matrix_a=fm_s_half,&
                      matrix_b=fm_work1,&
                      beta=0.0_dp,&
                      matrix_c=fm_work2)
      IF (debug) THEN
        CALL cp_dbcsr_write_sparse_matrix(sm_p,4,6,qs_env,para_env,&
                                       output_unit=output_unit)
        CALL write_fm_with_basis_info(fm_s_half,4,6,qs_env,para_env,&
                                      output_unit=output_unit)
        CALL write_fm_with_basis_info(fm_work2,4,6,qs_env,para_env,&
                                      output_unit=output_unit)
      END IF ! debug

      ! Copy occupation matrix to sparse matrix format, finally we are only
      ! interested in the diagonal (atomic) blocks, i.e. the previous full
      ! matrix product is not the most efficient choice, anyway.

      CALL copy_fm_to_dbcsr(fm_work2,sm_q,keep_sparsity=.TRUE.)

      ! E[DFT+U] = E[DFT] + E[U]
      !          = E[DFT] + (U - J)*(Tr(q) - Tr(q*q))/2

      ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
      !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
      !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

      ! Loop over all atomic kinds

      DO ikind=1,nkind

        ! Load the required atomic kind data

        CALL get_atomic_kind(atomic_kind_set(ikind),&
                             atom_list=atom_list,&
                             name=atomic_kind_name,&
                             natom=natom_of_kind)

        CALL get_qs_kind(qs_kind_set(ikind),&
                             dft_plus_u_atom=dft_plus_u_atom,&
                             l_of_dft_plus_u=lu,&
                             nsgf=nsgf_kind,&
                             basis_set=orb_basis_set,&
                             u_minus_j=u_minus_j,&
                             u_minus_j_target=u_minus_j_target,&
                             u_ramping=u_ramping,&
                             eps_u_ramping=eps_u_ramping,&
                             orbitals=orbitals,&
                             eps_scf=eps_scf,&
                             max_scf=max_scf,&
                             smear=smear)


        ! Check, if the atoms of this atomic kind need a DFT+U correction

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
        IF (.NOT.dft_plus_u_atom) CYCLE
        IF (lu < 0) CYCLE

        ! Apply U ramping if requested

        IF ((ispin == 1).AND.(u_ramping > 0.0_dp)) THEN
          IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
            u_minus_j = MIN(u_minus_j + u_ramping,u_minus_j_target)
            CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
          END IF
          IF (should_output.AND.(output_unit > 0)) THEN
            WRITE (UNIT=output_unit,FMT="(T3,A,3X,A,F0.3,A)")&
              "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)),&
              "U(eff) = ",u_minus_j*evolt," eV"
          END IF
        END IF

        IF (u_minus_j == 0.0_dp) CYCLE

        ! Load the required Gaussian basis set data

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgf,&
                               l=l,&
                               last_sgf=last_sgf,&
                               nset=nset,&
                               nshell=nshell)

        ! Count the relevant shell blocks of this atomic kind

        nsb = 0
        DO iset=1,nset
          DO ishell=1,nshell(iset)
            IF (l(ishell,iset) == lu) nsb = nsb + 1
          END DO
        END DO

        nsbsize = (2*lu + 1)
        n = nsb*nsbsize

        ALLOCATE (q_matrix(n,n),STAT=stat)
        CPASSERT(stat == 0)
        q_matrix(:,:) = 0.0_dp

        ! Print headline if requested

        IF (should_output.AND.(print_level > low_print_level)) THEN
          IF (output_unit > 0) THEN
            ALLOCATE (symbol(nsbsize),STAT=stat)
            CPASSERT(stat == 0)
            DO m=-lu,lu
              symbol(lu+m+1) = sgf_symbol(0,lu,m)
            END DO
            IF (nspin > 1) THEN
              WRITE (UNIT=spin_info,FMT="(A8,I2)") " of spin",ispin
            ELSE
              spin_info = ""
            END IF
            WRITE (UNIT=output_unit,FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))")&
              "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ",ikind,&
              ": "//TRIM(atomic_kind_name),&
              "Atom   Shell  ",(ADJUSTR(symbol(i)),i=1,nsbsize)," Trace"
            DEALLOCATE (symbol,STAT=stat)
            CPASSERT(stat == 0)
          END IF
        END IF

        ! Loop over all atoms of the current atomic kind

        DO iatom=1,natom_of_kind

          atom_a = atom_list(iatom)

          q_matrix(:,:) = 0.0_dp

          ! Get diagonal block

          CALL cp_dbcsr_get_block_p(matrix=sm_q,&
                                    row=atom_a,&
                                    col=atom_a,&
                                    block=q_block,&
                                    found=found)

          IF (ASSOCIATED(q_block)) THEN

            ! Calculate energy contribution to E(U)

            i = 0
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                IF (l(ishell,iset) /= lu) CYCLE
                DO isgf=first_sgf(ishell,iset),last_sgf(ishell,iset)
                  i = i + 1
                  j = 0
                  DO jset=1,nset
                    DO jshell=1,nshell(jset)
                      IF (l(jshell,jset) /= lu) CYCLE
                      DO jsgf=first_sgf(jshell,jset),last_sgf(jshell,jset)
                        j = j + 1
                        q_matrix(i,j) = q_block(isgf,jsgf)
                      END DO ! next contracted spherical Gaussian function "jsgf"
                    END DO ! next shell "jshell"
                  END DO ! next shell set "jset"
                END DO ! next contracted spherical Gaussian function "isgf"
              END DO ! next shell "ishell"
            END DO ! next shell set "iset"

            ! Perform the requested manipulations of the (initial) orbital occupations

            IF (ASSOCIATED(orbitals)) THEN
              IF ((qs_env%scf_env%iter_delta >= eps_scf).OR.&
                  ((qs_env%scf_env%outer_scf%iter_count == 0).AND.&
                   (qs_env%scf_env%iter_count <= max_scf))) THEN
                ALLOCATE (orb_occ(nsbsize),STAT=stat)
                CPASSERT(stat == 0)
                ALLOCATE (q_eigval(n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigval(:) = 0.0_dp
                ALLOCATE (q_eigvec(n,n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigvec(:,:) = 0.0_dp
                norb = SIZE(orbitals)
                CALL jacobi(q_matrix,q_eigval,q_eigvec)
                q_matrix(:,:) = 0.0_dp
                DO isb=1,nsb
                  trq = 0.0_dp
                  DO i=(isb-1)*nsbsize+1,isb*nsbsize
                    trq = trq + q_eigval(i)
                  END DO
                  IF (smear) THEN
                    occ = trq/REAL(norb,KIND=dp)
                  ELSE
                    occ = 1.0_dp/fspin
                  END IF
                  orb_occ(:) = .FALSE.
                  iloc = MAXLOC(q_eigvec(:,isb*nsbsize))
                  jsb = INT((iloc(1)-1)/nsbsize) + 1
                  i = 0
                  i0 = (jsb-1)*nsbsize + 1
                  iorb = -1000
                  DO j=i0,jsb*nsbsize
                    i = i + 1
                    IF (i > norb) THEN
                      DO m=-lu,lu
                        IF (.NOT.orb_occ(lu+m+1)) THEN
                          iorb = i0 + lu + m
                          orb_occ(lu+m+1) = .TRUE.
                        END IF
                      END DO
                    ELSE
                      iorb = i0 + lu + orbitals(i)
                      orb_occ(lu+orbitals(i)+1) = .TRUE.
                    END IF
                    CPASSERT(iorb /= -1000)
                    iloc = MAXLOC(q_eigvec(iorb,:))
                    q_eigval(iloc(1)) = MIN(occ,trq)
                    q_matrix(:,iloc(1)) = q_eigval(iloc(1))*q_eigvec(:,iloc(1)) ! backtransform left
                    trq = trq - q_eigval(iloc(1))
                  END DO
                END DO
                q_matrix(:,:) = MATMUL(q_matrix,TRANSPOSE(q_eigvec)) ! backtransform right
                DEALLOCATE (orb_occ,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigval,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigvec,STAT=stat)
                CPASSERT(stat == 0)
              END IF
            END IF ! orbitals associated

            trq = 0.0_dp
            trq2 = 0.0_dp

            DO i=1,n
              trq = trq + q_matrix(i,i)
              DO j=1,n
                trq2 = trq2 + q_matrix(i,j)*q_matrix(j,i)
              END DO
            END DO

            trq = fspin*trq
            trq2 = fspin*fspin*trq2

            ! Calculate energy contribution to E(U)

            energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*u_minus_j*(trq - trq2)/fspin

            ! Calculate potential V(U) = dE(U)/dq

            IF (.NOT.just_energy) THEN

              CALL cp_dbcsr_get_block_p(matrix=sm_v,&
                                        row=atom_a,&
                                        col=atom_a,&
                                        block=v_block,&
                                        found=found)
              CPASSERT(ASSOCIATED(v_block))

              i = 0
              DO iset=1,nset
                DO ishell=1,nshell(iset)
                  IF (l(ishell,iset) /= lu) CYCLE
                  DO isgf=first_sgf(ishell,iset),last_sgf(ishell,iset)
                    i = i + 1
                    j = 0
                    DO jset=1,nset
                      DO jshell=1,nshell(jset)
                        IF (l(jshell,jset) /= lu) CYCLE
                        DO jsgf=first_sgf(jshell,jset),last_sgf(jshell,jset)
                          j = j + 1
                          IF (isgf == jsgf) THEN
                            v_block(isgf,isgf) = u_minus_j*(0.5_dp - fspin*q_matrix(i,i))
                          ELSE
                            v_block(isgf,jsgf) = -u_minus_j*fspin*q_matrix(j,i)
                          END IF
                        END DO ! next contracted spherical Gaussian function "jsgf"
                      END DO ! next shell "jshell"
                    END DO ! next shell set "jset"
                  END DO ! next contracted spherical Gaussian function "isgf"
                END DO ! next shell "ishell"
              END DO ! next shell set "iset"

            END IF ! not just energy

          END IF ! q_block associated

          ! Consider print requests

          IF (should_output.AND.(print_level > low_print_level)) THEN
            CALL mp_sum(q_matrix, para_env%group)
            IF (output_unit > 0) THEN
              ALLOCATE (q_work(nsb,nsbsize),STAT=stat)
              CPASSERT(stat == 0)
              q_work(:,:) = 0.0_dp
              DO isb=1,nsb
                j = 0
                DO i=(isb-1)*nsbsize+1,isb*nsbsize
                  j = j + 1
                  q_work(isb,j) = q_matrix(i,i)
                END DO
              END DO
              DO isb=1,nsb
                WRITE (UNIT=output_unit,FMT="(T3,I6,2X,I6,2X,10F8.3)")&
                  atom_a,isb,q_work(isb,:),SUM(q_work(isb,:))
              END DO
              WRITE (UNIT=output_unit,FMT="(T12,A,2X,10F8.3)")&
                "Total",(SUM(q_work(:,i)),i=1,nsbsize),SUM(q_work)
              WRITE (UNIT=output_unit,FMT="(A)") ""
              DEALLOCATE (q_work,STAT=stat)
              CPASSERT(stat == 0)
              IF (debug) THEN
                ! Print the DFT+U occupation matrix
                WRITE (UNIT=output_unit,FMT="(T9,70I10)") (i,i=1,n)
                DO i=1,n
                  WRITE (UNIT=output_unit,FMT="(T3,I6,70F10.6)") i,q_matrix(i,:)
                END DO
                ! Print the eigenvalues and eigenvectors of the occupation matrix
                ALLOCATE (q_eigval(n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigval(:) = 0.0_dp
                ALLOCATE (q_eigvec(n,n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigvec(:,:) = 0.0_dp
                CALL jacobi(q_matrix,q_eigval,q_eigvec)
                WRITE (UNIT=output_unit,FMT="(/,T9,70I10)") (i,i=1,n)
                WRITE (UNIT=output_unit,FMT="(T9,71F10.6)") (q_eigval(i),i=1,n),&
                                                            SUM(q_eigval(1:n))
                DO i=1,n
                  WRITE (UNIT=output_unit,FMT="(T3,I6,70F10.6)") i,q_eigvec(i,:)
                END DO
                DEALLOCATE (q_eigval,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigvec,STAT=stat)
                CPASSERT(stat == 0)
              END IF ! debug
            END IF
            IF (debug) THEN
              ! Print the full atomic occupation matrix block
              ALLOCATE (q_work(nsgf_kind,nsgf_kind),STAT=stat)
              CPASSERT(stat == 0)
              q_work(:,:) = 0.0_dp
              IF (ASSOCIATED(q_block)) q_work(:,:) = q_block(:,:)
              CALL mp_sum(q_work, para_env%group)
              IF (output_unit > 0) THEN
                norb = SIZE(q_work,1)
                WRITE (UNIT=output_unit,FMT="(/,T9,200I10)") (i,i=1,norb)
                DO i=1,norb
                  WRITE (UNIT=output_unit,FMT="(T3,I6,200F10.6)") i,q_work(i,:)
                END DO
                ALLOCATE (q_eigval(norb),STAT=stat)
                CPASSERT(stat == 0)
                 q_eigval(:) = 0.0_dp
                ALLOCATE (q_eigvec(norb,norb),STAT=stat)
                CPASSERT(stat == 0)
                q_eigvec(:,:) = 0.0_dp
                CALL jacobi(q_work,q_eigval,q_eigvec)
                WRITE (UNIT=output_unit,FMT="(/,T9,200I10)") (i,i=1,norb)
                WRITE (UNIT=output_unit,FMT="(T9,201F10.6)") (q_eigval(i),i=1,norb),&
                                                             SUM(q_eigval(1:norb))
                DO i=1,norb
                  WRITE (UNIT=output_unit,FMT="(T3,I6,200F10.6)") i,q_eigvec(i,:)
                END DO
                DEALLOCATE (q_eigval,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigvec,STAT=stat)
                CPASSERT(stat == 0)
              END IF
              DEALLOCATE (q_work,STAT=stat)
              CPASSERT(stat == 0)
            END IF ! debug
          END IF ! should output

        END DO ! next atom "iatom" of atomic kind "ikind"

        IF (ALLOCATED(q_matrix)) THEN
          DEALLOCATE (q_matrix,STAT=stat)
          CPASSERT(stat == 0)
        END IF
      END DO ! next atomic kind "ikind"

      ! Add V(i,j)[U] to V(i,j)[DFT]

      IF (ASSOCIATED(sm_h)) THEN

        CALL cp_dbcsr_sm_fm_multiply(sm_v,fm_s_half,fm_work1,nsgf)
        CALL cp_fm_transpose(fm_work1,fm_work2)
        CALL cp_dbcsr_plus_fm_fm_t(sm_h,fm_s_half,fm_work2,nsgf)

      END IF ! An update of the Hamiltonian matrix is requested

      ! Calculate the contribution (non-Pulay part) to the derivatives
      ! w.r.t. the nuclear positions

      IF (ASSOCIATED(sm_w)) THEN

        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                       "Forces are not yet fully implemented for DFT+U method LOWDIN",&
                       only_ionode=.TRUE.)

      END IF ! W matrix update requested

    END DO ! next spin "ispin"

    ! Collect the energy contributions from all processes

    CALL mp_sum(energy%dft_plus_u, para_env%group)

    CALL cp_assert(energy%dft_plus_u >= 0.0_dp,cp_warning_level,&
                   cp_assertion_failed,routineP,&
                   "DFT+U energy contibution is negative possibly due "//&
                   "to unphysical Lowdin charges. Check your input, "//&
                   "if this warning persists or try a different method!",&
                   only_ionode=.TRUE.)

    ! Release (local) full matrices

    NULLIFY (fm_s_half)
    IF (fm_work1_local_alloc) THEN
      CALL cp_fm_release(matrix=fm_work1)
    END IF
    IF (fm_work2_local_alloc) THEN
      CALL cp_fm_release(matrix=fm_work2)
    END IF

    ! Release (local) sparse matrices

    CALL cp_dbcsr_deallocate_matrix(sm_q)
    CALL cp_dbcsr_deallocate_matrix(sm_v)

    CALL timestop(handle)

  END SUBROUTINE lowdin

! *****************************************************************************
!> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
!>                using a method based on the Mulliken population analysis
!>                \f[q_{\mu\nu} = \frac{1}{2} (P_{\mu\nu} S_{\nu\mu} +
!>                                             S_{\mu\nu} P_{\nu\mu})\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in]     qs_env Quickstep environment
!> \param orthonormal_basis ...
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param should_output ...
!> \param output_unit ...
!> \param print_level ...
!> \date          03.07.2008
!> \par
!>  \f{eqnarray*}{
!>   E^{\rm DFT+U} & = & E^{\rm DFT} + E^{\rm U}\\\
!>                 & = & E^{\rm DFT} + \frac{1}{2}\sum_A(U_A - J_A)(Tr(q_A) - Tr(q^2_A))\\[1ex]
!>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
!>                          & = & \frac{\partial E^{\rm DFT}}
!>                                     {\partial P_{\mu\nu}} +
!>                                \frac{\partial E^{\rm U}}
!>                                     {\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} + \sum_A
!>                                \frac{\partial E^{\rm U}}{\partial q_A}
!>                                \frac{\partial q_A}{\partial P_{\mu\nu}}\\\
!>  \f}
!> \author        Matthias Krack (MK)
!> \version       1.0
! *****************************************************************************
  SUBROUTINE mulliken(qs_env,orthonormal_basis,matrix_h,should_output,&
                      output_unit,print_level)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: orthonormal_basis
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h
    LOGICAL, INTENT(IN)                      :: should_output
    INTEGER, INTENT(IN)                      :: output_unit, print_level

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mulliken', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: spin_info
    CHARACTER(LEN=6), ALLOCATABLE, &
      DIMENSION(:)                           :: symbol
    CHARACTER(LEN=default_string_length)     :: atomic_kind_name
    INTEGER :: atom_a, handle, i, i0, iatom, ikind, iorb, isb, iset, isgf, &
      ishell, ispin, j, jsb, jset, jsgf, jshell, lu, m, max_scf, n, natom, &
      natom_of_kind, nkind, norb, nsb, nsbsize, nset, nsgf_kind, nspin, stat
    INTEGER, DIMENSION(1)                    :: iloc
    INTEGER, DIMENSION(:), POINTER           :: atom_list, nshell, orbitals
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgf, l, last_sgf
    LOGICAL                                  :: debug, dft_plus_u_atom, &
                                                failure, found, just_energy, &
                                                smear
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: is_plus_u_kind, orb_occ
    REAL(KIND=dp)                            :: eps_scf, eps_u_ramping, &
                                                fspin, occ, trq, trq2, &
                                                u_minus_j, u_minus_j_target, &
                                                u_ramping
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: q_eigval
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: q_eigvec, q_matrix, q_work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block, q_block, &
                                                s_block, v_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: kind_a
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: sm_h, sm_p, sm_q, sm_s, sm_v
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    debug = .FALSE. ! Set to .TRUE. to print debug information

    NULLIFY (atom_list)
    NULLIFY (atomic_kind_set)
    NULLIFY (qs_kind_set)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (first_sgf)
    NULLIFY (h_block)
    NULLIFY (matrix_p)
    NULLIFY (matrix_s)
    NULLIFY (l)
    NULLIFY (last_sgf)
    NULLIFY (nshell)
    NULLIFY (orb_basis_set)
    NULLIFY (p_block)
    NULLIFY (particle_set)
    NULLIFY (q_block)
    NULLIFY (rho)
    NULLIFY (s_block)
    NULLIFY (orbitals)
    NULLIFY (sm_h)
    NULLIFY (sm_p)
    NULLIFY (sm_q)
    NULLIFY (sm_s)
    NULLIFY (sm_v)
    NULLIFY (v_block)
    NULLIFY (para_env)

    failure = .FALSE.

    smear = .FALSE.
    max_scf = -1
    eps_scf = 1.0E30_dp

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    dft_control=dft_control,&
                    energy=energy,&
                    particle_set=particle_set,&
                    rho=rho,&
                    para_env=para_env)

    CPASSERT(ASSOCIATED(atomic_kind_set))
    CPASSERT(ASSOCIATED(dft_control))
    CPASSERT(ASSOCIATED(energy))
    CPASSERT(ASSOCIATED(particle_set))
    CPASSERT(ASSOCIATED(rho))

    IF (orthonormal_basis) THEN
      NULLIFY (sm_s)
    ELSE
      ! Get overlap matrix in sparse format
      CALL get_qs_env(qs_env=qs_env,&
                      matrix_s=matrix_s)
      CPASSERT(ASSOCIATED(matrix_s))
      sm_s => matrix_s(1)%matrix
    END IF

    ! Get density matrices in sparse format

    CALL qs_rho_get(rho, rho_ao=matrix_p)

    energy%dft_plus_u = 0.0_dp

    nspin = dft_control%nspins

    IF (nspin == 2) THEN
      fspin = 1.0_dp
    ELSE
      fspin = 0.5_dp
    END IF

    ! Get the total number of atoms, contracted spherical Gaussian basis
    ! functions, and atomic kinds

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (is_plus_u_kind(nkind),STAT=stat)
    CPASSERT(stat == 0)
    is_plus_u_kind(:) = .FALSE.

    IF (PRESENT(matrix_h)) THEN
      just_energy = .FALSE.
    ELSE
      just_energy = .TRUE.
    END IF

    ! Loop over all spins

    DO ispin=1,nspin

      IF (PRESENT(matrix_h)) THEN
        ! Hamiltonian matrix for spin ispin in sparse format
        sm_h => matrix_h(ispin)%matrix
      ELSE
        NULLIFY (sm_h)
      END IF

      ! Get density matrix for spin ispin in sparse format

      sm_p => matrix_p(ispin)%matrix

      IF (.NOT.ASSOCIATED(sm_q)) THEN
        ALLOCATE (sm_q,STAT=stat)
        CPASSERT(stat == 0)
        CALL cp_dbcsr_init(sm_q)
        CALL cp_dbcsr_get_block_diag(sm_p,sm_q)
      END IF
      CALL cp_dbcsr_set(sm_q,0.0_dp)

      IF (.NOT.ASSOCIATED(sm_v)) THEN
        ALLOCATE (sm_v,STAT=stat)
        CPASSERT(stat == 0)
        CALL cp_dbcsr_init(sm_v)
        CALL cp_dbcsr_get_block_diag(sm_p,sm_v)
      END IF
      CALL cp_dbcsr_set(sm_v,0.0_dp)

      DO iatom=1,natom

        CALL cp_dbcsr_get_block_p(matrix=sm_p,&
                                  row=iatom,&
                                  col=iatom,&
                                  block=p_block,&
                                  found=found)

        IF (.NOT.ASSOCIATED(p_block)) CYCLE

        CALL cp_dbcsr_get_block_p(matrix=sm_q,&
                                  row=iatom,&
                                  col=iatom,&
                                  block=q_block,&
                                  found=found)
        CPASSERT(ASSOCIATED(q_block))

        IF (orthonormal_basis) THEN
          ! S is the unit matrix
          DO isgf=1,SIZE(q_block,1)
            q_block(isgf,isgf) = p_block(isgf,isgf)
          END DO
        ELSE
          CALL cp_dbcsr_get_block_p(matrix=sm_s,&
                                    row=iatom,&
                                    col=iatom,&
                                    block=s_block,&
                                    found=found)
          CPASSERT(ASSOCIATED(s_block))
          ! Exploit that P and S are symmetric
          DO jsgf=1,SIZE(p_block,2)
            DO isgf=1,SIZE(p_block,1)
              q_block(isgf,jsgf) = p_block(isgf,jsgf)*s_block(isgf,jsgf)
            END DO
          END DO
        END IF ! orthonormal basis set

      END DO ! next atom "iatom"

      ! E[DFT+U] = E[DFT] + E[U]
      !          = E[DFT] + (U - J)*(Tr(q) - Tr(q*q))/2

      ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
      !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
      !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

      ! Loop over all atomic kinds

      DO ikind=1,nkind

        ! Load the required atomic kind data

        CALL get_atomic_kind(atomic_kind_set(ikind),&
                             atom_list=atom_list,&
                             name=atomic_kind_name,&
                             natom=natom_of_kind)

        CALL get_qs_kind(qs_kind_set(ikind),&
                         dft_plus_u_atom=dft_plus_u_atom,&
                         l_of_dft_plus_u=lu,&
                         nsgf=nsgf_kind,&
                         basis_set=orb_basis_set,&
                         u_minus_j=u_minus_j,&
                         u_minus_j_target=u_minus_j_target,&
                         u_ramping=u_ramping,&
                         eps_u_ramping=eps_u_ramping,&
                         orbitals=orbitals,&
                         eps_scf=eps_scf,&
                         max_scf=max_scf,&
                         smear=smear)

        ! Check, if the atoms of this atomic kind need a DFT+U correction

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
        IF (.NOT.dft_plus_u_atom) CYCLE
        IF (lu < 0) CYCLE

        ! Apply U ramping if requested

        IF ((ispin == 1).AND.(u_ramping > 0.0_dp)) THEN
          IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
            u_minus_j = MIN(u_minus_j + u_ramping,u_minus_j_target)
            CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
          END IF
          IF (should_output.AND.(output_unit > 0)) THEN
            WRITE (UNIT=output_unit,FMT="(T3,A,3X,A,F0.3,A)")&
              "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)),&
              "U(eff) = ",u_minus_j*evolt," eV"
          END IF
        END IF

        IF (u_minus_j == 0.0_dp) CYCLE

        is_plus_u_kind(ikind) = .TRUE.

        ! Load the required Gaussian basis set data

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgf,&
                               l=l,&
                               last_sgf=last_sgf,&
                               nset=nset,&
                               nshell=nshell)

        ! Count the relevant shell blocks of this atomic kind

        nsb = 0
        DO iset=1,nset
          DO ishell=1,nshell(iset)
            IF (l(ishell,iset) == lu) nsb = nsb + 1
          END DO
        END DO

        nsbsize = (2*lu + 1)
        n = nsb*nsbsize

        ALLOCATE (q_matrix(n,n),STAT=stat)
        CPASSERT(stat == 0)
        q_matrix(:,:) = 0.0_dp

        ! Print headline if requested

        IF (should_output.AND.(print_level > low_print_level)) THEN
          IF (output_unit > 0) THEN
            ALLOCATE (symbol(nsbsize),STAT=stat)
            CPASSERT(stat == 0)
            DO m=-lu,lu
              symbol(lu+m+1) = sgf_symbol(0,lu,m)
            END DO
            IF (nspin > 1) THEN
              WRITE (UNIT=spin_info,FMT="(A8,I2)") " of spin",ispin
            ELSE
              spin_info = ""
            END IF
            WRITE (UNIT=output_unit,FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))")&
              "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ",ikind,&
              ": "//TRIM(atomic_kind_name),&
              "Atom   Shell  ",(ADJUSTR(symbol(i)),i=1,nsbsize)," Trace"
            DEALLOCATE (symbol,STAT=stat)
            CPASSERT(stat == 0)
          END IF
        END IF

        ! Loop over all atoms of the current atomic kind

        DO iatom=1,natom_of_kind

          atom_a = atom_list(iatom)

          q_matrix(:,:) = 0.0_dp

          ! Get diagonal block

          CALL cp_dbcsr_get_block_p(matrix=sm_q,&
                                    row=atom_a,&
                                    col=atom_a,&
                                    block=q_block,&
                                    found=found)

          ! Calculate energy contribution to E(U)

          IF (ASSOCIATED(q_block)) THEN

            i = 0
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                IF (l(ishell,iset) /= lu) CYCLE
                DO isgf=first_sgf(ishell,iset),last_sgf(ishell,iset)
                  i = i + 1
                  j = 0
                  DO jset=1,nset
                    DO jshell=1,nshell(jset)
                      IF (l(jshell,jset) /= lu) CYCLE
                      DO jsgf=first_sgf(jshell,jset),last_sgf(jshell,jset)
                        j = j + 1
                        q_matrix(i,j) = q_block(isgf,jsgf)
                      END DO ! next contracted spherical Gaussian function "jsgf"
                    END DO ! next shell "jshell"
                  END DO ! next shell set "jset"
                END DO ! next contracted spherical Gaussian function "isgf"
              END DO ! next shell "ishell"
            END DO ! next shell set "iset"

            ! Perform the requested manipulations of the (initial) orbital occupations

            IF (ASSOCIATED(orbitals)) THEN
              IF ((qs_env%scf_env%iter_delta >= eps_scf).OR.&
                  ((qs_env%scf_env%outer_scf%iter_count == 0).AND.&
                   (qs_env%scf_env%iter_count <= max_scf))) THEN
                ALLOCATE (orb_occ(nsbsize),STAT=stat)
                CPASSERT(stat == 0)
                ALLOCATE (q_eigval(n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigval(:) = 0.0_dp
                ALLOCATE (q_eigvec(n,n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigvec(:,:) = 0.0_dp
                norb = SIZE(orbitals)
                CALL jacobi(q_matrix,q_eigval,q_eigvec)
                q_matrix(:,:) = 0.0_dp
                DO isb=1,nsb
                  trq = 0.0_dp
                  DO i=(isb-1)*nsbsize+1,isb*nsbsize
                    trq = trq + q_eigval(i)
                  END DO
                  IF (smear) THEN
                    occ = trq/REAL(norb,KIND=dp)
                  ELSE
                    occ = 1.0_dp/fspin
                  END IF
                  orb_occ(:) = .FALSE.
                  iloc = MAXLOC(q_eigvec(:,isb*nsbsize))
                  jsb = INT((iloc(1)-1)/nsbsize) + 1
                  i = 0
                  i0 = (jsb-1)*nsbsize + 1
                  iorb = -1000
                  DO j=i0,jsb*nsbsize
                    i = i + 1
                    IF (i > norb) THEN
                      DO m=-lu,lu
                        IF (.NOT.orb_occ(lu+m+1)) THEN
                          iorb = i0 + lu + m
                          orb_occ(lu+m+1) = .TRUE.
                        END IF
                      END DO
                    ELSE
                      iorb = i0 + lu + orbitals(i)
                      orb_occ(lu+orbitals(i)+1) = .TRUE.
                    END IF
                    CPASSERT(iorb /= -1000)
                    iloc = MAXLOC(q_eigvec(iorb,:))
                    q_eigval(iloc(1)) = MIN(occ,trq)
                    q_matrix(:,iloc(1)) = q_eigval(iloc(1))*q_eigvec(:,iloc(1)) ! backtransform left
                    trq = trq - q_eigval(iloc(1))
                  END DO
                END DO
                q_matrix(:,:) = MATMUL(q_matrix,TRANSPOSE(q_eigvec)) ! backtransform right
                DEALLOCATE (orb_occ,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigval,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigvec,STAT=stat)
                CPASSERT(stat == 0)
              END IF
            END IF ! orbitals associated

            trq = 0.0_dp
            trq2 = 0.0_dp

            DO i=1,n
              trq = trq + q_matrix(i,i)
              DO j=1,n
                trq2 = trq2 + q_matrix(i,j)*q_matrix(j,i)
              END DO
            END DO

            trq = fspin*trq
            trq2 = fspin*fspin*trq2

            ! Calculate energy contribution to E(U)

            energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*u_minus_j*(trq - trq2)/fspin

            ! Calculate potential V(U) = dE(U)/dq

            IF (.NOT.just_energy) THEN

              CALL cp_dbcsr_get_block_p(matrix=sm_v,&
                                        row=atom_a,&
                                        col=atom_a,&
                                        block=v_block,&
                                        found=found)
              CPASSERT(ASSOCIATED(v_block))

              i = 0
              DO iset=1,nset
                DO ishell=1,nshell(iset)
                  IF (l(ishell,iset) /= lu) CYCLE
                  DO isgf=first_sgf(ishell,iset),last_sgf(ishell,iset)
                    i = i + 1
                    j = 0
                    DO jset=1,nset
                      DO jshell=1,nshell(jset)
                        IF (l(jshell,jset) /= lu) CYCLE
                        DO jsgf=first_sgf(jshell,jset),last_sgf(jshell,jset)
                          j = j + 1
                          IF (isgf == jsgf) THEN
                            v_block(isgf,isgf) = u_minus_j*(0.5_dp - fspin*q_matrix(i,i))
                          ELSE
                            v_block(isgf,jsgf) = -u_minus_j*fspin*q_matrix(j,i)
                          END IF
                        END DO ! next contracted spherical Gaussian function "jsgf"
                      END DO ! next shell "jshell"
                    END DO ! next shell set "jset"
                  END DO ! next contracted spherical Gaussian function "isgf"
                END DO ! next shell "ishell"
              END DO ! next shell set "iset"

            END IF ! not just energy

          END IF ! q_block associated

          ! Consider print requests

          IF (should_output.AND.(print_level > low_print_level)) THEN
            CALL mp_sum(q_matrix, para_env%group)
            IF (output_unit > 0) THEN
              ALLOCATE (q_work(nsb,nsbsize),STAT=stat)
              CPASSERT(stat == 0)
              q_work(:,:) = 0.0_dp
              DO isb=1,nsb
                j = 0
                DO i=(isb-1)*nsbsize+1,isb*nsbsize
                  j = j + 1
                  q_work(isb,j) = q_matrix(i,i)
                END DO
              END DO
              DO isb=1,nsb
                WRITE (UNIT=output_unit,FMT="(T3,I6,2X,I6,2X,10F8.3)")&
                  atom_a,isb,q_work(isb,:),SUM(q_work(isb,:))
              END DO
              WRITE (UNIT=output_unit,FMT="(T12,A,2X,10F8.3)")&
                "Total",(SUM(q_work(:,i)),i=1,nsbsize),SUM(q_work)
              WRITE (UNIT=output_unit,FMT="(A)") ""
              DEALLOCATE (q_work,STAT=stat)
              CPASSERT(stat == 0)
              IF (debug) THEN
                ! Print the DFT+U occupation matrix
                WRITE (UNIT=output_unit,FMT="(T9,70I10)") (i,i=1,n)
                DO i=1,n
                  WRITE (UNIT=output_unit,FMT="(T3,I6,70F10.6)") i,q_matrix(i,:)
                END DO
                ! Print the eigenvalues and eigenvectors of the occupation matrix
                ALLOCATE (q_eigval(n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigval(:) = 0.0_dp
                ALLOCATE (q_eigvec(n,n),STAT=stat)
                CPASSERT(stat == 0)
                q_eigvec(:,:) = 0.0_dp
                CALL jacobi(q_matrix,q_eigval,q_eigvec)
                WRITE (UNIT=output_unit,FMT="(/,T9,70I10)") (i,i=1,n)
                WRITE (UNIT=output_unit,FMT="(T9,71F10.6)") (q_eigval(i),i=1,n),&
                                                            SUM(q_eigval(1:n))
                DO i=1,n
                  WRITE (UNIT=output_unit,FMT="(T3,I6,70F10.6)") i,q_eigvec(i,:)
                END DO
                DEALLOCATE (q_eigval,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigvec,STAT=stat)
                CPASSERT(stat == 0)
              END IF ! debug
            END IF
            IF (debug) THEN
              ! Print the full atomic occupation matrix block
              ALLOCATE (q_work(nsgf_kind,nsgf_kind),STAT=stat)
              CPASSERT(stat == 0)
              q_work(:,:) = 0.0_dp
              IF (ASSOCIATED(q_block)) q_work(:,:) = q_block(:,:)
              CALL mp_sum(q_work, para_env%group)
              IF (output_unit > 0) THEN
                norb = SIZE(q_work,1)
                WRITE (UNIT=output_unit,FMT="(/,T9,200I10)") (i,i=1,norb)
                DO i=1,norb
                  WRITE (UNIT=output_unit,FMT="(T3,I6,200F10.6)") i,q_work(i,:)
                END DO
                ALLOCATE (q_eigval(norb),STAT=stat)
                CPASSERT(stat == 0)
                 q_eigval(:) = 0.0_dp
                ALLOCATE (q_eigvec(norb,norb),STAT=stat)
                CPASSERT(stat == 0)
                q_eigvec(:,:) = 0.0_dp
                CALL jacobi(q_work,q_eigval,q_eigvec)
                WRITE (UNIT=output_unit,FMT="(/,T9,200I10)") (i,i=1,norb)
                WRITE (UNIT=output_unit,FMT="(T9,201F10.6)") (q_eigval(i),i=1,norb),&
                                                             SUM(q_eigval(1:norb))
                DO i=1,norb
                  WRITE (UNIT=output_unit,FMT="(T3,I6,200F10.6)") i,q_eigvec(i,:)
                END DO
                DEALLOCATE (q_eigval,STAT=stat)
                CPASSERT(stat == 0)
                DEALLOCATE (q_eigvec,STAT=stat)
                CPASSERT(stat == 0)
              END IF
              DEALLOCATE (q_work,STAT=stat)
              CPASSERT(stat == 0)
            END IF ! debug
          END IF ! should output

        END DO ! next atom "iatom" of atomic kind "ikind"

        IF (ALLOCATED(q_matrix)) THEN
          DEALLOCATE (q_matrix,STAT=stat)
          CPASSERT(stat == 0)
        END IF

      END DO ! next atomic kind "ikind"

      ! Add V(i,j)[U] to V(i,j)[DFT]

      IF (ASSOCIATED(sm_h)) THEN

        DO ikind=1,nkind

          IF (.NOT.is_plus_u_kind(ikind)) CYCLE

          kind_a => atomic_kind_set(ikind)

          CALL get_atomic_kind(atomic_kind=kind_a,&
                               atom_list=atom_list,&
                               natom=natom_of_kind)

          DO iatom=1,natom_of_kind

            atom_a = atom_list(iatom)

            CALL cp_dbcsr_get_block_p(matrix=sm_h,&
                                      row=atom_a,&
                                      col=atom_a,&
                                      block=h_block,&
                                      found=found)

            IF (.NOT.ASSOCIATED(h_block)) CYCLE

            CALL cp_dbcsr_get_block_p(matrix=sm_v,&
                                      row=atom_a,&
                                      col=atom_a,&
                                      block=v_block,&
                                      found=found)
            CPASSERT(ASSOCIATED(v_block))

            IF (orthonormal_basis) THEN
              DO isgf=1,SIZE(h_block,1)
                h_block(isgf,isgf) = h_block(isgf,isgf) + v_block(isgf,isgf)
              END DO
            ELSE
               CALL cp_dbcsr_get_block_p(matrix=sm_s,&
                                         row=atom_a,&
                                         col=atom_a,&
                                         block=s_block,&
                                         found=found)
              CPASSERT(ASSOCIATED(s_block))
              DO jsgf=1,SIZE(h_block,2)
                DO isgf=1,SIZE(h_block,1)
                  h_block(isgf,jsgf) = h_block(isgf,jsgf) + v_block(isgf,jsgf)*s_block(isgf,jsgf)
                END DO
              END DO
            END IF ! orthonormal basis set

          END DO ! next atom "iatom" of atomic kind "ikind"

        END DO ! Next atomic kind "ikind"

      END IF ! An update of the Hamiltonian matrix is requested

    END DO ! next spin "ispin"

    ! Collect the energy contributions from all processes

    CALL mp_sum(energy%dft_plus_u, para_env%group)

    CALL cp_assert(energy%dft_plus_u >= 0.0_dp,cp_warning_level,&
                   cp_assertion_failed,routineP,&
                   "DFT+U energy contibution is negative possibly due "//&
                   "to unphysical Mulliken charges. Check your input, "//&
                   "if this warning persists or try a different method!",&
                   only_ionode=.TRUE.)

    CALL cp_dbcsr_deallocate_matrix(sm_q)
    CALL cp_dbcsr_deallocate_matrix(sm_v)

    CALL timestop(handle)

  END SUBROUTINE mulliken

! *****************************************************************************
!> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
!>                using a method based on Mulliken charges
!>                \f[q_\mu = \sum_\nu \frac{1}{2}(P_{\mu\nu} S_{\nu\mu} +
!>                                                S_{\mu\nu} P_{\nu\mu})
!>                         = \sum_\nu P_{\mu\nu} S_{\nu\mu}\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in]     qs_env Quickstep environment
!> \param orthonormal_basis ...
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param[in,out] matrix_w Energy weighted density matrices for each spin
!> \param should_output ...
!> \param output_unit ...
!> \param print_level ...
!> \date          11.01.2008
!> \par
!>  \f{eqnarray*}{
!>   E^{\rm DFT+U} & = & E^{\rm DFT} +  E^{\rm U}\\\
!>                 & = & E^{\rm DFT} +  \frac{1}{2}(U - J)\sum_\mu (q_\mu - q_\mu^2)\\[1ex]
!>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
!>                          & = & \frac{\partial E^{\rm DFT}}
!>                                     {\partial P_{\mu\nu}} +
!>                                \frac{\partial E^{\rm U}}
!>                                     {\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} +
!>                                \frac{\partial E^{\rm U}}{\partial q_\mu}
!>                                \frac{\partial q_\mu}{\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} +
!>                                \frac{1}{2}(U - J)(1 - q_\mu - q_\nu) S_{\mu\nu}\\\
!>  \f}
!> \author        Matthias Krack (MK)
!> \version       1.0
!> \note          The use of any full matrices was avoided. Thus no ScaLAPACK
!>                calls are performed
! *****************************************************************************
  SUBROUTINE mulliken_charges(qs_env,orthonormal_basis,matrix_h,matrix_w,&
                              should_output,output_unit,print_level)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: orthonormal_basis
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_w
    LOGICAL, INTENT(IN)                      :: should_output
    INTEGER, INTENT(IN)                      :: output_unit, print_level

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mulliken_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: spin_info
    CHARACTER(LEN=6), ALLOCATABLE, &
      DIMENSION(:)                           :: symbol
    CHARACTER(LEN=default_string_length)     :: atomic_kind_name
    INTEGER :: atom_a, blk, handle, i, iatom, ikind, isb, iset, isgf, ishell, &
      ispin, jatom, jsgf, lu, m, natom, natom_of_kind, nkind, nsb, nset, &
      nsgf, nspin, sgf, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf_atom
    INTEGER, DIMENSION(:), POINTER           :: atom_list, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgf, l, last_sgf
    LOGICAL                                  :: dft_plus_u_atom, failure, &
                                                found, just_energy
    REAL(KIND=dp)                            :: eps_u_ramping, fspin, q, &
                                                u_minus_j, u_minus_j_target, &
                                                u_ramping, v
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dEdq, trps
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: q_ii
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block, s_block, &
                                                w_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: sm_h, sm_p, sm_s, sm_w
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    NULLIFY (atom_list)
    NULLIFY (atomic_kind_set)
    NULLIFY (qs_kind_set)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (first_sgf)
    NULLIFY (h_block)
    NULLIFY (matrix_p)
    NULLIFY (matrix_s)
    NULLIFY (l)
    NULLIFY (last_sgf)
    NULLIFY (nshell)
    NULLIFY (orb_basis_set)
    NULLIFY (p_block)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (s_block)
    NULLIFY (sm_h)
    NULLIFY (sm_p)
    NULLIFY (sm_s)
    NULLIFY (w_block)
    NULLIFY (para_env)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    dft_control=dft_control,&
                    energy=energy,&
                    particle_set=particle_set,&
                    rho=rho,&
                    para_env=para_env)

    CPASSERT(ASSOCIATED(atomic_kind_set))
    CPASSERT(ASSOCIATED(dft_control))
    CPASSERT(ASSOCIATED(energy))
    CPASSERT(ASSOCIATED(particle_set))
    CPASSERT(ASSOCIATED(rho))

    IF (orthonormal_basis) THEN
      NULLIFY (sm_s)
    ELSE
      ! Get overlap matrix in sparse format
      CALL get_qs_env(qs_env=qs_env,&
                      matrix_s=matrix_s)
      CPASSERT(ASSOCIATED(matrix_s))
      sm_s => matrix_s(1)%matrix
    END IF

    ! Get density matrices in sparse format

    CALL qs_rho_get(rho, rho_ao=matrix_p)

    energy%dft_plus_u = 0.0_dp

    nspin = dft_control%nspins

    IF (nspin == 2) THEN
      fspin = 1.0_dp
    ELSE
      fspin = 0.5_dp
    END IF

    ! Get the total number of atoms, contracted spherical Gaussian basis
    ! functions, and atomic kinds

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (first_sgf_atom(natom),STAT=stat)
    CPASSERT(stat==0)
    first_sgf_atom(:) = 0

    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf_atom)

    ALLOCATE (trps(nsgf),STAT=stat)
    CPASSERT(stat==0)
    trps(:) = 0.0_dp

    IF (PRESENT(matrix_h).OR.PRESENT(matrix_w)) THEN
      ALLOCATE (dEdq(nsgf),STAT=stat)
      CPASSERT(stat == 0)
      just_energy  = .FALSE.
    ELSE
      just_energy = .TRUE.
    END IF

    ! Loop over all spins

    DO ispin=1,nspin

      IF (PRESENT(matrix_h)) THEN
        ! Hamiltonian matrix for spin ispin in sparse format
        sm_h => matrix_h(ispin)%matrix
        dEdq(:) = 0.0_dp
      ELSE
        NULLIFY (sm_h)
      END IF

      IF (PRESENT(matrix_w)) THEN
        ! Energy weighted density matrix for spin ispin in sparse format
        sm_w => matrix_w(ispin)%matrix
        dEdq(:) = 0.0_dp
      ELSE
        NULLIFY (sm_w)
      END IF

      sm_p => matrix_p(ispin)%matrix ! Density matrix for spin ispin in sparse format

      ! Calculate Trace(P*S) assuming symmetric matrices

      trps(:) = 0.0_dp

      CALL cp_dbcsr_iterator_start(iter,sm_p)

      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

        CALL cp_dbcsr_iterator_next_block(iter,iatom,jatom,p_block,blk)

        IF (orthonormal_basis) THEN

          IF (iatom /= jatom) CYCLE

          IF (ASSOCIATED(p_block)) THEN
            sgf = first_sgf_atom(iatom)
            DO isgf=1,SIZE(p_block,1)
              trps(sgf) = trps(sgf) + p_block(isgf,isgf)
              sgf = sgf + 1
            END DO
          END IF

        ELSE

          CALL cp_dbcsr_get_block_p(matrix=sm_s,&
                                    row=iatom,&
                                    col=jatom,&
                                    block=s_block,&
                                    found=found)
          CPASSERT(ASSOCIATED(s_block))

          sgf = first_sgf_atom(jatom)
          DO jsgf=1,SIZE(p_block,2)
            DO isgf=1,SIZE(p_block,1)
              trps(sgf) = trps(sgf) + p_block(isgf,jsgf)*s_block(isgf,jsgf)
            END DO
            sgf = sgf + 1
          END DO

          IF (iatom /= jatom) THEN
            sgf = first_sgf_atom(iatom)
            DO isgf=1,SIZE(p_block,1)
              DO jsgf=1,SIZE(p_block,2)
                trps(sgf) = trps(sgf) + p_block(isgf,jsgf)*s_block(isgf,jsgf)
              END DO
              sgf = sgf + 1
            END DO
          END IF

        END IF ! orthonormal basis set

      END DO ! next atom "iatom"

      CALL cp_dbcsr_iterator_stop(iter)

      CALL mp_sum(trps, para_env%group)

      ! q <- Trace(PS)

      ! E[DFT+U] = E[DFT] + E[U]
      !          = E[DFT] + (U - J)*(q - q**2))/2

      ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
      !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
      !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

      ! Loop over all atomic kinds

      DO ikind=1,nkind

        ! Load the required atomic kind data
        CALL get_atomic_kind(atomic_kind_set(ikind),&
                             atom_list=atom_list,&
                             name=atomic_kind_name,&
                             natom=natom_of_kind)

        CALL get_qs_kind(qs_kind_set(ikind),&
                         dft_plus_u_atom=dft_plus_u_atom,&
                         l_of_dft_plus_u=lu,&
                         basis_set=orb_basis_set,&
                         u_minus_j=u_minus_j,&
                         u_minus_j_target=u_minus_j_target,&
                         u_ramping=u_ramping,&
                         eps_u_ramping=eps_u_ramping)

        ! Check, if this atom needs a DFT+U correction

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
        IF (.NOT.dft_plus_u_atom) CYCLE
        IF (lu < 0) CYCLE

        ! Apply U ramping if requested

        IF ((ispin == 1).AND.(u_ramping > 0.0_dp)) THEN
          IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
            u_minus_j = MIN(u_minus_j + u_ramping,u_minus_j_target)
            CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
          END IF
          IF (should_output.AND.(output_unit > 0)) THEN
            WRITE (UNIT=output_unit,FMT="(T3,A,3X,A,F0.3,A)")&
              "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)),&
              "U(eff) = ",u_minus_j*evolt," eV"
          END IF
        END IF

        IF (u_minus_j == 0.0_dp) CYCLE

        ! Load the required Gaussian basis set data

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgf,&
                               l=l,&
                               last_sgf=last_sgf,&
                               nset=nset,&
                               nshell=nshell)

        ! Count the relevant shell blocks of this atomic kind

        nsb = 0
        DO iset=1,nset
          DO ishell=1,nshell(iset)
            IF (l(ishell,iset) == lu) nsb = nsb + 1
          END DO
        END DO

        ALLOCATE (q_ii(nsb,2*lu+1),STAT=stat)
        CPASSERT(stat == 0)

        ! Print headline if requested

        IF (should_output.AND.(print_level > low_print_level)) THEN
          IF (output_unit > 0) THEN
            ALLOCATE (symbol(2*lu+1),STAT=stat)
            CPASSERT(stat == 0)
            DO m=-lu,lu
              symbol(lu+m+1) = sgf_symbol(0,lu,m)
            END DO
            IF (nspin > 1) THEN
              WRITE (UNIT=spin_info,FMT="(A8,I2)") " of spin",ispin
            ELSE
              spin_info = ""
            END IF
            WRITE (UNIT=output_unit,FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))")&
              "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ",ikind,&
              ": "//TRIM(atomic_kind_name),&
              "Atom   Shell  ",(ADJUSTR(symbol(i)),i=1,2*lu+1)," Trace"
            DEALLOCATE (symbol,STAT=stat)
            CPASSERT(stat == 0)
          END IF
        END IF

        ! Loop over all atoms of the current atomic kind

        DO iatom=1,natom_of_kind

          atom_a = atom_list(iatom)

          q_ii(:,:) = 0.0_dp

          ! Get diagonal block

          CALL cp_dbcsr_get_block_p(matrix=sm_p,&
                                    row=atom_a,&
                                    col=atom_a,&
                                    block=p_block,&
                                    found=found)

          ! Calculate E(U) and dE(U)/dq

          IF (ASSOCIATED(p_block)) THEN

            sgf = first_sgf_atom(atom_a)

            isb = 0
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                IF (l(ishell,iset) == lu) THEN
                  isb = isb + 1
                  i = 0
                  DO isgf=first_sgf(ishell,iset),last_sgf(ishell,iset)
                    q = fspin*trps(sgf)
                    i = i + 1
                    q_ii(isb,i) = q
                    energy%dft_plus_u = energy%dft_plus_u +&
                                        0.5_dp*u_minus_j*(q - q**2)/fspin
                    IF (.NOT.just_energy) THEN
                      dEdq(sgf) = dEdq(sgf) + u_minus_j*(0.5_dp - q)
                    END IF
                    sgf = sgf + 1
                  END DO ! next contracted spherical Gaussian function "isgf"
                ELSE
                  sgf = sgf + last_sgf(ishell,iset) - first_sgf(ishell,iset) + 1
                END IF ! angular momentum requested for DFT+U correction
              END DO ! next shell "ishell"
            END DO ! next shell set "iset"

          END IF ! this process is the owner of the sparse matrix block?

          ! Consider print requests

          IF (should_output.AND.(print_level > low_print_level)) THEN
            CALL mp_sum(q_ii, para_env%group)
            IF (output_unit > 0) THEN
              DO isb=1,nsb
                WRITE (UNIT=output_unit,FMT="(T3,I6,2X,I6,2X,10F8.3)")&
                  atom_a,isb,q_ii(isb,:),SUM(q_ii(isb,:))
              END DO
              WRITE (UNIT=output_unit,FMT="(T12,A,2X,10F8.3)")&
                "Total",(SUM(q_ii(:,i)),i=1,2*lu+1),SUM(q_ii)
              WRITE (UNIT=output_unit,FMT="(A)") ""
            END IF
          END IF ! should output

        END DO ! next atom "iatom" of atomic kind "ikind"

        IF (ALLOCATED(q_ii)) THEN
          DEALLOCATE (q_ii,STAT=stat)
          CPASSERT(stat == 0)
        END IF

      END DO ! next atomic kind "ikind"

      IF (.NOT.just_energy) THEN
        CALL mp_sum(dEdq, para_env%group)
      END IF

      ! Add V(i,j)[U] to V(i,j)[DFT]

      IF (ASSOCIATED(sm_h)) THEN

        CALL cp_dbcsr_iterator_start(iter,sm_h)

        DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

          CALL cp_dbcsr_iterator_next_block(iter,iatom,jatom,h_block,blk)

          IF (orthonormal_basis) THEN

            IF (iatom /= jatom) CYCLE

            IF (ASSOCIATED(h_block)) THEN
              sgf = first_sgf_atom(iatom)
              DO isgf=1,SIZE(h_block,1)
                h_block(isgf,isgf) = h_block(isgf,isgf) + dEdq(sgf)
                sgf = sgf + 1
              END DO
            END IF

          ELSE

            ! Request katom just to check for consistent sparse matrix pattern

            CALL cp_dbcsr_get_block_p(matrix=sm_s,&
                                      row=iatom,&
                                      col=jatom,&
                                      block=s_block,&
                                      found=found)
            CPASSERT(ASSOCIATED(s_block))

            ! Consider the symmetric form 1/2*(P*S + S*P) for the calculation

            sgf = first_sgf_atom(iatom)

            DO isgf=1,SIZE(h_block,1)
              IF (dEdq(sgf) /= 0.0_dp) THEN
                v = 0.5_dp*dEdq(sgf)
                DO jsgf=1,SIZE(h_block,2)
                  h_block(isgf,jsgf) = h_block(isgf,jsgf) + v*s_block(isgf,jsgf)
                END DO
              END IF
              sgf = sgf + 1
            END DO

            sgf = first_sgf_atom(jatom)

            DO jsgf=1,SIZE(h_block,2)
              IF (dEdq(sgf) /= 0.0_dp) THEN
                v = 0.5_dp*dEdq(sgf)
                DO isgf=1,SIZE(h_block,1)
                  h_block(isgf,jsgf) = h_block(isgf,jsgf) + v*s_block(isgf,jsgf)
                END DO
              END IF
              sgf = sgf + 1
            END DO

          END IF ! orthonormal basis set

        END DO ! Next atom "iatom"

        CALL cp_dbcsr_iterator_stop(iter)

      END IF ! An update of the Hamiltonian matrix is requested

      ! Calculate the contribution (non-Pulay part) to the derivatives
      ! w.r.t. the nuclear positions, which requires an update of the
      ! energy weighted density W.

      IF (ASSOCIATED(sm_w).AND.(.NOT.orthonormal_basis)) THEN

        CALL cp_dbcsr_iterator_start(iter,sm_p)

        DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

          CALL cp_dbcsr_iterator_next_block(iter,iatom,jatom,p_block,blk)

          ! Skip the diagonal blocks of the W matrix

          IF (iatom == jatom) CYCLE

          ! Request katom just to check for consistent sparse matrix patterns

          CALL cp_dbcsr_get_block_p(matrix=sm_w,&
                                    row=iatom,&
                                    col=jatom,&
                                    block=w_block,&
                                    found=found)
          CPASSERT(ASSOCIATED(w_block))

          ! Consider the symmetric form 1/2*(P*S + S*P) for the calculation

          sgf = first_sgf_atom(iatom)

          DO isgf=1,SIZE(w_block,1)
            IF (dEdq(sgf) /= 0.0_dp) THEN
              v = -0.5_dp*dEdq(sgf)
              DO jsgf=1,SIZE(w_block,2)
                w_block(isgf,jsgf) = w_block(isgf,jsgf) + v*p_block(isgf,jsgf)
              END DO
            END IF
            sgf = sgf + 1
          END DO

          sgf = first_sgf_atom(jatom)

          DO jsgf=1,SIZE(w_block,2)
            IF (dEdq(sgf) /= 0.0_dp) THEN
              v = -0.5_dp*dEdq(sgf)
              DO isgf=1,SIZE(w_block,1)
                w_block(isgf,jsgf) = w_block(isgf,jsgf) + v*p_block(isgf,jsgf)
              END DO
            END IF
            sgf = sgf + 1
          END DO

        END DO ! next block node "jatom"

        CALL cp_dbcsr_iterator_stop(iter)

      END IF ! W matrix update requested

   END DO ! next spin "ispin"

    ! Collect the energy contributions from all processes

    CALL mp_sum(energy%dft_plus_u, para_env%group)

    CALL cp_assert(energy%dft_plus_u >= 0.0_dp,cp_warning_level,&
                   cp_assertion_failed,routineP,&
                   "DFT+U energy contibution is negative possibly due "//&
                   "to unphysical Mulliken charges. Check your input, "//&
                   "if this warning persists or try a different method!",&
                   only_ionode=.TRUE.)

    ! Release local work storage

    IF (ALLOCATED(first_sgf_atom)) THEN
      DEALLOCATE (first_sgf_atom,STAT=stat)
      CPASSERT(stat==0)
    END IF

    IF (ALLOCATED(trps)) THEN
      DEALLOCATE (trps,STAT=stat)
      CPASSERT(stat==0)
    END IF

    IF (ALLOCATED(dEdq)) THEN
      DEALLOCATE (dEdq,STAT=stat)
      CPASSERT(stat==0)
    END IF

   CALL timestop(handle)

  END SUBROUTINE mulliken_charges

END MODULE dft_plus_u
