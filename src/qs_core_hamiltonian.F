!i----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the core Hamiltonian integral matrix <a|H|b> over
!>      Cartesian Gaussian-type functions.
!>
!>      <a|H|b> = <a|T|b> + <a|V|b>
!>
!>      Kinetic energy:
!>
!>      <a|T|b> = <a|-nabla**2/2|b>
!>                \_______________/
!>                        |
!>                     kinetic
!>
!>      Nuclear potential energy:
!>
!>      a) Allelectron calculation:
!>
!>                          erfc(r)
!>         <a|V|b> = -Z*<a|---------|b>
!>                             r
!>
!>                          1 - erf(r)
!>                 = -Z*<a|------------|b>
!>                              r
!>
!>                           1           erf(r)
!>                 = -Z*(<a|---|b> - <a|--------|b>)
!>                           r             r
!>
!>                           1
!>                 = -Z*(<a|---|b> - N*<ab||c>)
!>                           r
!>
!>                      -Z
!>                 = <a|---|b> + Z*N*<ab||c>
!>                       r
!>                   \_______/       \_____/
!>                       |              |
!>                    nuclear        coulomb
!>
!>      b) Pseudopotential calculation (Goedecker, Teter and Hutter; GTH):
!>
!>         <a|V|b> = <a|(V(local) + V(non-local))|b>
!>
!>                 = <a|(V(local)|b> + <a|V(non-local))|b>
!>
!>         <a|V(local)|b> = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!>                             (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!>                              C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!>
!>         <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!> \par Literature
!>      S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!>      C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!>      M. Krack and M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Joost VandeVondele (April 2003) : added LSD forces
!>      - Non-redundant calculation of the non-local part of the GTH PP
!>        (22.05.2003,MK)
!>      - New parallelization scheme (27.06.2003,MK)
!>      - OpenMP version (07.12.2003,JGH)
!>      - Binary search loop for VPPNL operators (09.01.2004,JGH,MK)
!>      - Refactoring of pseudopotential and nuclear attraction integrals (25.02.2009,JGH)
!>      - General refactoring (01.10.2010,JGH)
!>      - Refactoring related to the new kinetic energy and overlap routines (07.2014,JGH)
!> \author Matthias Krack (14.09.2000,21.03.02)
! *****************************************************************************
MODULE qs_core_hamiltonian
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE core_ae,                         ONLY: build_core_ae
  USE core_ppl,                        ONLY: build_core_ppl
  USE core_ppnl,                       ONLY: build_core_ppnl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_deallocate_matrix_set, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
       dbcsr_distribution_obj
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_matrix_dist,&
                                             cp_dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE ep_qs_types,                     ONLY: ep_qs_type
  USE input_constants,                 ONLY: do_ppl_analytic,&
                                             kg_tnadd_atomic,&
                                             rel_none,&
                                             rel_trans_atom,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kg_tnadd_mat,                    ONLY: build_tnadd_mat
  USE kinds,                           ONLY: dp
  USE lri_environment_methods,         ONLY: build_lri_matrices
  USE lri_environment_types,           ONLY: lri_environment_type
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_kinetic,                      ONLY: build_kinetic_matrix
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_oce_methods,                  ONLY: build_oce_matrices
  USE qs_oce_types,                    ONLY: allocate_oce_set,&
                                             create_oce_set,&
                                             oce_matrix_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_core_hamiltonian'

  PUBLIC :: build_core_hamiltonian_matrix
  PUBLIC :: dump_info_core_hamiltonian

CONTAINS

! *****************************************************************************
!> \brief Cosntruction of the QS Core Hamiltonian Matrix
!> \param qs_env ...
!> \param calculate_forces ...
!> \param error ...
!> \author Creation (11.03.2002,MK)
!>      Non-redundant calculation of the non-local part of the GTH PP (22.05.2003,MK)
!>      New parallelization scheme (27.06.2003,MK)
! *****************************************************************************
  SUBROUTINE build_core_hamiltonian_matrix(qs_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'build_core_hamiltonian_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, maxder, natom, nder, &
                                                nders, neighbor_list_id, &
                                                nkind, stat
    LOGICAL :: all_potential_present, do_ep, failure, gth_potential_present, &
      kgpot, lrigpw, ofdft, ppl_present, use_virial
    REAL(KIND=dp)                            :: eps_filter, eps_fit, eps_ppnl
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_h, matrix_p, &
      matrix_s, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_t, matrix_w
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux_fit, &
                                                sab_aux_fit_vs_orb, sab_orb, &
                                                sac_ae, sac_ppl, sap_oce, &
                                                sap_ppnl
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//"_forces",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY (atomic_kind_set)
    NULLIFY (qs_kind_set)
    NULLIFY (force)
    NULLIFY (matrix_h)
    NULLIFY (matrix_s)
    NULLIFY (matrix_t)
    NULLIFY (oce)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (sab_orb)
    NULLIFY (sac_ae)
    NULLIFY (sac_ppl)
    NULLIFY (sap_ppnl)
    NULLIFY (sap_oce)
    NULLIFY (sab_aux_fit)
    NULLIFY (sab_aux_fit_vs_orb)
    NULLIFY (scf_control)
    NULLIFY (matrix_w)
    NULLIFY (matrix_p)
    NULLIFY (para_env)

    ! is this a orbital-free method calculation
    ofdft = qs_env%dft_control%qs_control%ofgpw
    ! is this a LRI method calculation
    lrigpw = qs_env%dft_control%qs_control%lrigpw
    ! is this a KG method calculation using atomic TNADD potentials
    IF (qs_env%dft_control%qs_control%do_kg) THEN
       kgpot = (qs_env%kg_env%tnadd_method==kg_tnadd_atomic)
    ELSE
       kgpot = .FALSE.
    END IF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    matrix_h=matrix_h,&
                    kinetic=matrix_t,&
                    matrix_s=matrix_s,&
                    oce=oce, &
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    sac_ae=sac_ae,&
                    sac_ppl=sac_ppl,&
                    sap_ppnl=sap_ppnl,&
                    sap_oce=sap_oce,&
                    sab_aux_fit=sab_aux_fit,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    neighbor_list_id=neighbor_list_id,&
                    scf_control=scf_control,&
                    virial=virial,&
                    para_env=para_env,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    nders = 0
    IF (calculate_forces) THEN
       nder = 1

       CALL get_qs_env(qs_env=qs_env,force=force,ep_qs_env=ep_qs_env,matrix_w=matrix_w,error=error)

       do_ep=ASSOCIATED(ep_qs_env)
       IF (do_ep) do_ep=ep_qs_env%ep_active.AND.calculate_forces

       CALL get_qs_env(qs_env=qs_env, rho=rho,error=error)
       matrix_p => rho%rho_ao

       IF (do_ep) THEN
          IF (ASSOCIATED(ep_qs_env%dH_coeffs)) THEN
             matrix_p => ep_qs_env%dH_coeffs
          END IF
          IF (ASSOCIATED(ep_qs_env%dS_coeffs)) THEN
             matrix_w => ep_qs_env%dS_coeffs

          END IF
       END IF

       !     *** If LSD, then combine alpha density and beta density to
       !     *** total density: alpha <- alpha + beta   and
       !     *** spin density:   beta <- alpha - beta
       !     (since all things can be computed based on the sum of these matrices anyway)
       !     (matrix_p is restored at the end of the run, matrix_w is left in its modified state
       !     (as it should not be needed afterwards)
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                         alpha_scalar= 1.0_dp, beta_scalar=1.0_dp,error=error)
          CALL cp_dbcsr_add(matrix_p(2)%matrix, matrix_p(1)%matrix, &
                         alpha_scalar=-2.0_dp, beta_scalar=1.0_dp,error=error)
          CALL cp_dbcsr_add(matrix_w(1)%matrix, matrix_w(2)%matrix, &
                         alpha_scalar= 1.0_dp, beta_scalar=1.0_dp,error=error)
       END IF
    ELSE
       IF (cp_print_key_should_output(logger%iter_info,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error)/=0.OR.&
            BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
            "DFT%PRINT%OPTICAL_CONDUCTIVITY",error=error),cp_p_file)) THEN
          nder = 1
       ELSE
          nder = 0
       END IF
    END IF

    IF ((cp_print_key_should_output(logger%iter_info,qs_env%input,&
         "DFT%PRINT%AO_MATRICES/OVERLAP",error=error)/=0.AND.&
         BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
         "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error),cp_p_file)).OR.&
         BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
         "DFT%PRINT%OPTICAL_CONDUCTIVITY",error=error),cp_p_file)) THEN
       nders = 1
    END IF

    ! the delta pulse in the periodic case needs the momentum operator,
    ! which is equivalent to the derivative of the overlap matrix
    IF (ASSOCIATED(qs_env%dft_control%rtp_control)) THEN
       IF (qs_env%dft_control%rtp_control%apply_delta_pulse .AND. &
           qs_env%dft_control%rtp_control%periodic) THEN
           nders = 1
       ENDIF
    ENDIF

    ! filter for new matrices
    eps_filter = qs_env%dft_control%qs_control%eps_filter_matrix

    maxder = ncoset(nder)

    IF (calculate_forces) THEN
       ! S matrix
       CALL build_overlap_matrix(qs_env,nderivative=nders,matrix_s=matrix_s,&
                                 matrix_name="OVERLAP MATRIX",&
                                 basis_set_id_a=use_orb_basis_set,&
                                 basis_set_id_b=use_orb_basis_set, &
                                 sab_nl=sab_orb,calculate_forces=.TRUE.,&
                                 matrix_p=matrix_w(1)%matrix,&
                                 error=error)
       ! T matrix
       IF (.NOT.ofdft) &
          CALL build_kinetic_matrix(qs_env,matrix_t=matrix_t,&
                                 matrix_name="KINETIC ENERGY MATRIX",&
                                 basis_set_id=use_orb_basis_set,&
                                 sab_nl=sab_orb,calculate_forces=.TRUE.,&
                                 matrix_p=matrix_p(1)%matrix,&
                                 eps_filter=eps_filter,&
                                 error=error)
    ELSE
       ! S matrix
       CALL build_overlap_matrix(qs_env,nderivative=nders,matrix_s=matrix_s,&
                                 matrix_name="OVERLAP MATRIX",&
                                 basis_set_id_a=use_orb_basis_set,&
                                 basis_set_id_b=use_orb_basis_set, &
                                 sab_nl=sab_orb,&
                                 error=error)
       ! T matrix
       IF (.NOT.ofdft) &
          CALL build_kinetic_matrix(qs_env,matrix_t=matrix_t,&
                                 matrix_name="KINETIC ENERGY MATRIX",&
                                 basis_set_id=use_orb_basis_set,&
                                 sab_nl=sab_orb,&
                                 eps_filter=eps_filter,&
                                 error=error)
       IF( qs_env%dft_control%do_admm) THEN
         NULLIFY(matrix_s_aux_fit,matrix_s_aux_fit_vs_orb)
         CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s_aux_fit, error=error)
         CALL build_overlap_matrix(qs_env,matrix_s=matrix_s_aux_fit,&
                                   matrix_name="AUX_FIT_OVERLAP",&
                                   basis_set_id_a=use_aux_fit_basis_set,&
                                   basis_set_id_b=use_aux_fit_basis_set, &
                                   sab_nl=sab_aux_fit,&
                                   error=error)
         CALL set_qs_env(qs_env=qs_env,matrix_s_aux_fit=matrix_s_aux_fit,error=error)
         CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, error=error)
         CALL build_overlap_matrix(qs_env,matrix_s=matrix_s_aux_fit_vs_orb,&
                                   matrix_name="MIXED_OVERLAP",&
                                   basis_set_id_a=use_aux_fit_basis_set,&
                                   basis_set_id_b=use_orb_basis_set, &
                                   sab_nl=sab_aux_fit_vs_orb,&
                                   error=error)
         CALL set_qs_env(qs_env=qs_env,matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,error=error)
       END IF
    END IF

    ! initialize H matrix
    CALL cp_dbcsr_allocate_matrix_set(matrix_h,1,error)
    ALLOCATE(matrix_h(1)%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(matrix_h(1)%matrix, error=error)
    IF (.NOT.ofdft) THEN
       CALL cp_dbcsr_copy(matrix_h(1)%matrix,matrix_t(1)%matrix,&
           name="CORE HAMILTONIAN MATRIX",error=error)
       ! relativistic atomic correction to kinetic energy
       IF(qs_env%rel_control%rel_method /= rel_none)THEN
          IF(qs_env%rel_control%rel_transformation == rel_trans_atom)THEN
             CALL build_atomic_relmat(matrix_h, atomic_kind_set, qs_kind_set, particle_set, error)
          END IF
       END IF
    END IF

    !   *** Allocate the matrix of coefficients for one center expansions
    IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
       CALL create_oce_set(oce,error=error)
       CALL allocate_oce_set(oce,natom,nkind,maxder,error=error)
       ! force analytic ppl calcuation for GAPW methods
       qs_env%dft_control%qs_control%do_ppl_method=do_ppl_analytic
    ENDIF

    CALL get_qs_kind_set(qs_kind_set,&
                         gth_potential_present=gth_potential_present,&
                         all_potential_present=all_potential_present)
    ppl_present = ASSOCIATED(sac_ppl)

    IF (calculate_forces) THEN
       ! *** If LSD, then recover alpha density and beta density     ***
       ! *** from the total density (1) and the spin density (2)     ***
       ! *** The W matrix is neglected, since it will be destroyed   ***
       ! *** in the calling force routine after leaving this routine ***
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                         alpha_scalar= 0.5_dp, beta_scalar=0.5_dp,error=error)
          CALL cp_dbcsr_add(matrix_p(2)%matrix, matrix_p(1)%matrix, &
                         alpha_scalar=-1.0_dp, beta_scalar=1.0_dp,error=error)
       END IF
    END IF

    ! *** compute the ppl contribution to the core hamiltonian ***
    IF (ppl_present) THEN
      IF(qs_env%dft_control%qs_control%do_ppl_method==do_ppl_analytic) THEN
        CALL build_core_ppl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, error)
      END IF
    END IF

    ! *** compute the nuclear attraction contribution to the core hamiltonian ***
    IF (all_potential_present) THEN
      CALL build_core_ae(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ae, error)
    END IF

    ! *** compute the ppnl contribution to the core hamiltonian ***
    eps_ppnl = qs_env%dft_control%qs_control%eps_ppnl
    CALL build_core_ppnl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, error)

    ! *** GAPW one-center-expansion (oce) matrices
    IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
      eps_fit  = qs_env%dft_control%qs_control%gapw_control%eps_fit
      IF (ASSOCIATED(sap_oce)) &
         CALL build_oce_matrices(oce%intac,calculate_forces,nder,qs_kind_set,particle_set,&
              sap_oce,eps_fit,error)
    END IF

    ! *** LRIGPW  matrices
    IF(lrigpw) THEN
       CALL get_qs_env(qs_env=qs_env,lri_env=lri_env,error=error)
       CALL build_lri_matrices(lri_env,qs_env,calculate_forces,error)
    END IF

    ! *** KG atomic potentials for nonadditive kinetic energy
    IF (kgpot) THEN
       CALL get_qs_env(qs_env=qs_env,kg_env=kg_env,dbcsr_dist=dbcsr_dist,error=error)
       CALL build_tnadd_mat(kg_env, matrix_p, force, virial, calculate_forces, use_virial,&
                    qs_kind_set, atomic_kind_set, particle_set, sab_orb, dbcsr_dist, error)
    END IF

    ! *** Put the core Hamiltonian matrix in the QS environment ***
    CALL set_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,&
                    kinetic=matrix_t,&
                    matrix_h=matrix_h,&
                    oce=oce,error=error)

    IF(qs_env%rel_control%rel_method /= rel_none)THEN
       IF(qs_env%rel_control%rel_transformation /= rel_trans_atom)THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Relativistic corrections of this type are currently not implemented", &
               error=error, error_level=cp_failure_level)
       END IF
    END IF

    ! Print matrices if requested
    CALL dump_info_core_hamiltonian(matrix_s, matrix_t, matrix_h, &
         qs_env, calculate_forces, para_env, error)
    CALL timestop(handle)

  END SUBROUTINE build_core_hamiltonian_matrix

! *****************************************************************************
!> \brief Adds atomic blocks of relativistic correction for the kinetic energy
!> \param matrix_h ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_atomic_relmat(matrix_h, atomic_kind_set, qs_kind_set, particle_set, error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_atomic_relmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, iatom, ikind, istat, &
                                                jatom, natom
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hblock, reltmat
    TYPE(cp_dbcsr_iterator)                  :: iter

    failure = .FALSE.
    natom = SIZE(particle_set)
    ALLOCATE (kind_of(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set,kind_of=kind_of)

    CALL cp_dbcsr_iterator_start(iter,matrix_h(1)%matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, hblock, blk)
       IF (iatom==jatom) THEN
          ikind = kind_of(iatom)
          CALL get_qs_kind(qs_kind_set(ikind),reltmat=reltmat)
          hblock = hblock + reltmat
       END IF
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    DEALLOCATE (kind_of,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE build_atomic_relmat

! *****************************************************************************
!> \brief Possibly prints matrices after the construction of the Core
!>     Hamiltonian Matrix
!> \param matrix_s ...
!> \param matrix_t ...
!> \param matrix_h ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dump_info_core_hamiltonian(matrix_s, matrix_t, matrix_h, qs_env, &
       calculate_forces, para_env, error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_t, matrix_h
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_info_core_hamiltonian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, output_unit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_v
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    NULLIFY(logger, matrix_v)
    logger => cp_error_get_logger(error)

    ! *** Print the distribution of the overlap matrix blocks ***
    ! *** this duplicates causes duplicate printing at the force calc ***
    IF (.NOT. calculate_forces) THEN
       IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"PRINT%DISTRIBUTION",error=error),cp_p_file)) THEN
           output_unit = cp_print_key_unit_nr(logger,qs_env%input,"PRINT%DISTRIBUTION",&
                                              extension=".distribution",error=error)
           CALL cp_dbcsr_write_matrix_dist(matrix_s(1)%matrix,output_unit,para_env,error)
           CALL cp_print_key_finished_output(output_unit,logger,qs_env%input,"PRINT%DISTRIBUTION",&
                error=error)
       ENDIF
    ENDIF

    ! *** Print the overlap integral matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",&
            extension=".Log",error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_s(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
           "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error),cp_p_file)) THEN
          DO i=2,SIZE(matrix_s)
             CALL cp_dbcsr_write_sparse_matrix(matrix_s(i)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
          END DO
       END IF
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
           "DFT%PRINT%AO_MATRICES/OVERLAP", error=error)
    END IF

    ! *** Print the kinetic energy integral matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/KINETIC_ENERGY",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/KINETIC_ENERGY",&
            extension=".Log",error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_t(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/KINETIC_ENERGY", error=error)
    END IF

    ! *** Print the potential energy matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY",&
            extension=".Log",error=error)
       CALL cp_dbcsr_allocate_matrix_set(matrix_v,1,error=error)
       ALLOCATE(matrix_v(1)%matrix)
       CALL cp_dbcsr_init(matrix_v(1)%matrix, error=error)
       CALL cp_dbcsr_copy(matrix_v(1)%matrix,matrix_h(1)%matrix,&
            name="POTENTIAL ENERGY MATRIX",error=error)
       CALL cp_dbcsr_add(matrix_v(1)%matrix,matrix_t(1)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=-1.0_dp,error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_v(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_dbcsr_deallocate_matrix_set(matrix_v,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY", error=error)
    END IF

    ! *** Print the core Hamiltonian matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log",error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE dump_info_core_hamiltonian

END MODULE qs_core_hamiltonian
