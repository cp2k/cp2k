!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
!> \author Matthias Krack (14.09.2000,21.03.02)
! *****************************************************************************
MODULE qs_core_hamiltonian

  USE ai_kinetic,                      ONLY: kinetic
  USE ai_overlap_new,                  ONLY: overlap
  USE ai_overlap_ppl,                  ONLY: overlap_ppl
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE core_ppnl,                       ONLY: build_core_ppnl
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE ep_qs_types,                     ONLY: ep_qs_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE harris_env_types,                ONLY: harris_env_type
  USE input_constants,                 ONLY: diag_ot,&
                                             linear_response_run,&
                                             lr_precond_kinetic,&
                                             rel_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_all_potential,                ONLY: all_integrals,&
                                             vall_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       find_neighbor_list, first_list, first_node, get_neighbor_list, &
       get_neighbor_list_set, get_neighbor_node, neighbor_list_set_p_type, &
       neighbor_list_type, neighbor_node_type, next, reduced_3c_list_type
  USE qs_oce_methods,                  ONLY: set_up_oce
  USE qs_oce_types,                    ONLY: &
       allocate_oce_set, allocate_vtriple, build_reduced_3c_lists, &
       create_oce_set, deallocate_vtriple, dist_list_4oce_type, &
       oce_matrix_type, prepare_oce_coeff, retrieve_sac_list, vtriple_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_output,            ONLY: write_matrix_distribution,&
                                             write_sparse_matrix
  USE sparse_matrix_types,             ONLY: &
       add_block_node, add_matrices, allocate_matrix, allocate_matrix_set, &
       cp_sm_scale_and_add, deallocate_matrix_set, get_block_node, &
       real_matrix_p_type, replicate_matrix_structure
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_core_hamiltonian'

  ! *** Local types ***
! *****************************************************************************
  TYPE vppl_type
     REAL(KIND = dp)                          :: alpha_ppl,ppl_radius
     REAL(KIND = dp), DIMENSION(:), POINTER   :: cexp_ppl,r2
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: r
     INTEGER, DIMENSION(:), POINTER           :: neighbor
  END TYPE vppl_type
  
  PUBLIC :: build_core_hamiltonian_matrix

CONTAINS

! *****************************************************************************
!> \brief Cosntruction of the QS Core Hamiltonian Matrix
!> \author Creation (11.03.2002,MK)
!>      Non-redundant calculation of the non-local part of the GTH PP (22.05.2003,MK)
!>      New parallelization scheme (27.06.2003,MK) 
! *****************************************************************************
  SUBROUTINE build_core_hamiltonian_matrix(qs_env,globenv,calculate_forces,harris,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: harris
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'build_core_hamiltonian_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=80)                        :: name
    INTEGER :: atom_a, atom_b, atom_c, first_col, handle, i, iab, iac, iat, &
      iatom, ico, icol, ikind, ilist, inode, irow, iset, ithread, j, jatom, &
      jco, jkind, jset, katom, kkind, kneighbor, knode, last_jatom, ldai, &
      ldsab, lr_precond, maxblock, maxco, maxdco, maxder, maxl, maxlgto, &
      maxlppl, maxsgf, n, nat, natom, ncoa, ncob, ncol, nder, &
      neighbor_list_id, nkind, nlist, nneighbor, nnode, nnode_c, nrow, nseta, &
      nsetb, nsgf, nthread, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, first_sgf, &
                                                last_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                lb_max, lb_min, npgfa, npgfb, &
                                                nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: all_potential_present, build_kinetic_energy_matrix, do_ep, &
      epr_calculation, failure, gth_potential_present, harris_flag, &
      new_atom_b, oce_present, ppl_present, return_s_derivatives, &
      return_t_derivatives, rho0_present, use_virial, xas_calculation
    REAL(KIND=dp)                            :: dab, dac, dbc, f, f0, rab2, &
                                                rac2, rbc2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, rab, rac, &
                                                rbc
    REAL(KIND=dp), DIMENSION(3, 3)           :: virial_ab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, hab, p_block, pab, &
                                                rpgfa, rpgfb, sab, sphi_a, &
                                                sphi_b, w_block, work, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: habt, pabt, sabt, sdab, workt
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: sdabt
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces, sint, tint
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: oceht, ocest, sintt, tintt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dist_list_4oce_type), &
      DIMENSION(:), POINTER                  :: dist_list
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_3c, sac_oce, &
                                                sac_ppl, sac_ppnl, sbc_3c, &
                                                sbc_oce, sbc_ppnl
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local, &
                                                sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node, &
                                                sac_ppl_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s, &
                                                matrix_t, matrix_w
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce, &
                                                reduced_3c_rho0
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: lr_section
    TYPE(vall_type), DIMENSION(:), POINTER   :: vall
    TYPE(vall_type), DIMENSION(:, :), &
      POINTER                                :: vallt
    TYPE(virial_type), POINTER               :: virial
    TYPE(vppl_type), DIMENSION(:), POINTER   :: vppl
    TYPE(vppl_type), DIMENSION(:, :), &
      POINTER                                :: vpplt
    TYPE(vtriple_type), DIMENSION(:, :), &
      POINTER                                :: VH_3c_list

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//" (forces)",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY (atomic_kind_set)
    NULLIFY (distribution_2d)
    NULLIFY (force)
    NULLIFY (matrix_h)
    NULLIFY (habt)
    NULLIFY (oce,dist_list)
    NULLIFY (VH_3c_list)
    NULLIFY (matrix_p)
    NULLIFY (pabt)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (matrix_s)
    NULLIFY (sab_orb)
    NULLIFY (sabt)
    NULLIFY (sac_ppl)
    NULLIFY (sac_ppnl)
    NULLIFY (sbc_ppnl)
    NULLIFY (sac_oce)
    NULLIFY (sbc_oce)
    NULLIFY (sac_3c)
    NULLIFY (sbc_3c)
    NULLIFY (scf_control)
    NULLIFY (sdabt)
    NULLIFY (matrix_t)
    NULLIFY (matrix_w)
    NULLIFY (workt)
    NULLIFY (harris_env)
    NULLIFY (para_env)
    NULLIFY (atom_list)

    ! Default
    oce_present = .FALSE.
    rho0_present = .FALSE.

    IF (PRESENT(harris)) THEN
      harris_flag = harris
    ELSE
      harris_flag = .FALSE.
    END IF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    matrix_h=matrix_h,&
                    kinetic=matrix_t,&
                    oce=oce, &
                    particle_set=particle_set,&
                    matrix_s=matrix_s,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sac_ppnl=sac_ppnl,&
                    sbc_ppnl=sbc_ppnl,&
                    sac_oce=sac_oce,&
                    sbc_oce=sbc_oce,&
                    VH_3c_list=VH_3c_list,&
                    sac_3c=sac_3c,&
                    sbc_3c=sbc_3c,&
                    neighbor_list_id=neighbor_list_id,&
                    scf_control=scf_control,&
                    distribution_2d=distribution_2d,&
                    virial=virial,&
                    para_env=para_env,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    
    IF (calculate_forces) THEN
       nder = 1
       
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

       CALL get_qs_env(qs_env=qs_env,force=force,ep_qs_env=ep_qs_env,matrix_w=matrix_w,error=error)
       do_ep=ASSOCIATED(ep_qs_env)
       IF (do_ep) do_ep=ep_qs_env%ep_active.AND.calculate_forces
       IF (harris_flag) THEN
          CALL get_qs_env(qs_env=qs_env, harris_env=harris_env,error=error)
          matrix_p => harris_env%rho%rho_ao
       ELSE
          CALL get_qs_env(qs_env=qs_env, rho=rho,error=error)
          matrix_p => rho%rho_ao
       END IF
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
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,matrix_p(2)%matrix,1.0_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-2.0_dp,matrix_p(1)%matrix,1.0_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_w(1)%matrix, 1.0_dp,matrix_w(2)%matrix,1.0_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_w(2)%matrix,-2.0_dp,matrix_w(1)%matrix,1.0_dp,error=error)
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
       return_s_derivatives = .TRUE.
    ELSE
       return_s_derivatives = .FALSE.
    END IF
    
    maxder = ncoset(nder)

    CALL allocate_matrix_set(matrix_s,maxder,error=error)
    
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)
    
    ALLOCATE (first_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,last_sgf=last_sgf,&
         error=error)

    ! assert here that the S-matrix is of the type distribution_2d if this one is non-null
    ! we assert that the sparsity of s is determined by the neighbor list that is used to build it
    ! is that right ? i.e. this implies that even if forces are computed the sparsity of s remains the same.
    ! this seems necessary anyway.
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

!   *** Allocate the matrix of coefficients for one center expansions
    IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
       CALL create_oce_set(oce,error=error)
       CALL allocate_oce_set(oce,natom,nkind,maxder,error=error)
       CALL prepare_oce_coeff(particle_set,oce,natom,nsgf,&
            first_sgf,last_sgf,distribution_2d,error=error)
       ALLOCATE(dist_list(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1,natom
          dist_list(i)%ndist = 0
          NULLIFY(dist_list(i)%rdist)
       END DO
    ENDIF
    
    DEALLOCATE (first_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=2,maxder
       cgfsym = cgf_symbol(1,indco(1:3,i))
       name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE OVERLAP MATRIX "//&
            "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                      TARGET=matrix_s(i)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="antisymmetric",error=error)
   END DO
   
   !   *** Allocate the kinetic energy integral ***
   !   *** matrix (only needed for printing or  ***
   !   *** the OT preconditioner)               ***
   lr_precond = 0
   IF(globenv%run_type_id==linear_response_run) THEN
      lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
      CALL section_vals_val_get(lr_section,"PRECONDITIONER",i_val=lr_precond,error=error)
   END IF
   
   ! The transition potential method may need the kinetic matrix
   xas_calculation = qs_env%dft_control%do_xas_calculation
   
   ! The g tensor needs the kinetic matrix 
   CALL section_vals_val_get(qs_env%input,& 
        "PROPERTIES%LINRES%EPR%_SECTION_PARAMETERS_",l_val=epr_calculation,error=error)

   build_kinetic_energy_matrix = (scf_control%use_ot.OR.&
        cp_print_key_should_output(logger%iter_info,&
        qs_env%input,"DFT%PRINT%AO_MATRICES/KINETIC_ENERGY",error=error)/=0.OR.&
        cp_print_key_should_output(logger%iter_info,&
        qs_env%input,"DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY",error=error)/=0.OR.&
        (scf_control%use_diag .AND. scf_control%diagonalization%method == diag_ot ) .OR.&
        lr_precond == lr_precond_kinetic .OR. xas_calculation .OR. epr_calculation)
   ! Assume false by default
   return_t_derivatives = .FALSE.
  
   IF(qs_env%rel_control%rel_method/=rel_none)THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
              message="DKH method not yet available", &
              error=error, error_level=cp_failure_level) 
   END IF
 
   IF (build_kinetic_energy_matrix) THEN
      IF (cp_print_key_should_output(logger%iter_info,qs_env%input,&
           "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error)/=0) THEN
         CALL allocate_matrix_set(matrix_t,maxder,error=error)
         return_t_derivatives = .TRUE.
      ELSE
         CALL allocate_matrix_set(matrix_t,1,error=error)
         return_t_derivatives = .FALSE.
      END IF
      CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                      TARGET=matrix_t(1)%matrix,&
                                      target_name="KINETIC ENERGY MATRIX",error=error)
      DO i=2,SIZE(matrix_t)
         cgfsym = cgf_symbol(1,indco(1:3,i))
         name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE KINETIC ENERGY MATRIX "//&
              "W.R.T. THE NUCLEAR COORDINATES"
         CALL compress(name)
         CALL uppercase(name)
         CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                         TARGET=matrix_t(i)%matrix,&
                                         target_name=name,&
                                         target_symmetry="antisymmetric",error=error)
      END DO
   END IF

   ! *** Allocate the core Hamiltonian matrix ***
   CALL allocate_matrix_set(matrix_h,1,error=error)
   CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                   TARGET=matrix_h(1)%matrix,&
                                   target_name="CORE HAMILTONIAN MATRIX",error=error)

   ! *** Allocate work storage ***
   CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                            gth_potential_present=gth_potential_present,&
                            all_potential_present=all_potential_present, &
                            maxco=maxco,&
                            maxlgto=maxlgto,&
                            maxlppl=maxlppl,&
                            maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl)
    CALL init_orbital_pointers(maxl+nder+1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxdco,&
                             maxder=MAX(1,nder))

    ldsab = MAX(maxco,maxsgf)
    maxblock = MAX(SIZE(matrix_s),maxder)
    nthread = 1
    ldai = ncoset(maxl+nder+1)
    
    CALL reallocate(habt,1,ldsab,1,ldsab*maxder,0,nthread-1)
    CALL reallocate(sabt,1,ldsab,1,ldsab*maxblock,0,nthread-1)
    CALL reallocate(sdabt,1,maxdco,1,maxco,1,4,0,nthread-1)
    CALL reallocate(workt,1,ldsab,1,ldsab*maxder,0,nthread-1)
    IF (calculate_forces) THEN
       CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    END IF
    
    ALLOCATE (sintt(SIZE(matrix_s),0:nthread-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,SIZE(sintt,1)
       DO j=0,nthread-1
          NULLIFY (sintt(i,j)%block)
       END DO
    END DO
    
    IF (build_kinetic_energy_matrix) THEN
       ALLOCATE (tintt(SIZE(matrix_t),0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i=1,SIZE(tintt,1)
          DO j=0,nthread-1
             NULLIFY (tintt(i,j)%block)
          END DO
       END DO
    END IF

    ! ***  Initialize blocks for the one center expansione ***
    IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
       oce_present = ASSOCIATED(sac_oce)
       rho0_present = ASSOCIATED(sac_3c)
       IF(rho0_present) THEN
          CALL allocate_vtriple(VH_3c_list,nkind,natom,error=error)
       END IF
       ALLOCATE (oceht(SIZE(matrix_s),0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (ocest(SIZE(matrix_s),0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       DO i=1,SIZE(oceht,1)
          DO j=0,nthread-1
             NULLIFY (oceht(i,j)%block)
             NULLIFY (ocest(i,j)%block)
          END DO
       END DO
    ENDIF
    
    ! *** Load GTH pseudo potential data (local part -> PPL) ***
    ppl_present = ASSOCIATED(sac_ppl)

    IF (ppl_present) THEN
       ALLOCATE (vpplt(nkind,0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO kkind=1,nkind
          atomic_kind => atomic_kind_set(kkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,gth_potential=gth_potential)

          DO j=0,nthread-1
             NULLIFY (vpplt(kkind,j)%cexp_ppl)
             IF (ASSOCIATED(gth_potential)) THEN
                CALL get_potential(potential=gth_potential,&
                                   alpha_ppl=vpplt(kkind,j)%alpha_ppl,&
                                   cexp_ppl=vpplt(kkind,j)%cexp_ppl,&
                                   ppl_radius=vpplt(kkind,j)%ppl_radius)
             END IF
             NULLIFY (vpplt(kkind,j)%r2)
             NULLIFY (vpplt(kkind,j)%neighbor)
             NULLIFY (vpplt(kkind,j)%r)
          END DO
       END DO
    END IF

    IF (all_potential_present) THEN
       ALLOCATE(vallt(nkind,0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO kkind=1,nkind
          atomic_kind => atomic_kind_set(kkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               all_potential=all_potential)
          DO j=0,nthread-1
             NULLIFY(vallt(kkind,j)%neighbor)
             NULLIFY(vallt(kkind,j)%force_c)
             IF (ASSOCIATED(all_potential)) THEN
                CALL get_potential(potential=all_potential,&
                                   alpha_core_charge=vallt(kkind,j)%alpha_c,&
                                   zeff=vallt(kkind,j)%zeta_c,&
                                   ccore_charge=vallt(kkind,j)%core_charge,&
                                   core_charge_radius=vallt(kkind,j)%core_radius)
             ELSE
                vallt(kkind,j)%alpha_c = 0.0_dp
                vallt(kkind,j)%zeta_c = 0.0_dp
                vallt(kkind,j)%core_charge = 0.0_dp
                vallt(kkind,j)%core_radius = 0.0_dp
             END IF
          END DO
       END DO
    END IF

    ! *** Retrieve data of SAC_OCE and of SAC_3C neighbors of atom "iatom" ***
    DO ikind = 1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,atom_list=atom_list)
       DO iat = 1,nat
          iatom = atom_list(iat)
          IF(oce_present) THEN
             CALL retrieve_sac_list(oce%voce,sac_oce,iatom,ikind,nkind, error=error)
          END IF
          IF(rho0_present) THEN
             CALL retrieve_sac_list(VH_3c_list,sac_3c,iatom,ikind,nkind,error=error)
          ENDIF
       END DO
    END DO

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
       
       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                              first_sgf=first_sgfa,&
                              lmax=la_max,&
                              lmin=la_min,&
                              npgf=npgfa,&
                              nset=nseta,&
                              nsgf_set=nsgfa,&
                              pgf_radius=rpgfa,&
                              set_radius=set_radius_a,&
                              sphi=sphi_a,&
                              zet=zeta)
       
       DO jkind=1,nkind
          atomic_kind => atomic_kind_set(jkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
          
          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 first_sgf=first_sgfb,&
                                 lmax=lb_max,&
                                 lmin=lb_min,&
                                 npgf=npgfb,&
                                 nset=nsetb,&
                                 nsgf_set=nsgfb,&
                                 pgf_radius=rpgfb,&
                                 set_radius=set_radius_b,&
                                 sphi=sphi_b,&
                                 zet=zetb)
          
          iab = ikind + nkind*(jkind - 1)
          
          IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
          CALL get_neighbor_list_set(neighbor_list_set=sab_orb(iab)%neighbor_list_set,nlist=nlist)
          
          NULLIFY ( sab_orb_neighbor_list )
          ithread = 0
          hab => habt(:,:,ithread)
          sab => sabt(:,:,ithread)
          sdab => sdabt(:,:,:,ithread)
          work => workt(:,:,ithread)
          IF (calculate_forces) pab => pabt(:,:,ithread)
          sint => sintt(:,ithread)
          IF (build_kinetic_energy_matrix) THEN
             tint => tintt(:,ithread)
          END IF
          IF (ppl_present) THEN
             vppl => vpplt(:,ithread)
          END IF
          IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
             oceh => oceht(:,ithread)
             oces => ocest(:,ithread)
          END IF
          IF (all_potential_present) THEN
             vall => vallt(:,ithread)
          END IF
          ALLOCATE (ai_work(ldai,ldai,MAX(1,ncoset(maxlppl)),ncoset(nder+1)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
                sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,atom=iatom,nnode=nnode)

             IF (calculate_forces) atom_a = atom_of_kind(iatom)

             ! *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***
             IF (ppl_present) THEN
                DO kkind=1,nkind
                   iac = ikind + nkind*(kkind - 1)
                   IF (.NOT.ASSOCIATED(sac_ppl(iac)%neighbor_list_set)) CYCLE
                   sac_ppl_neighbor_list => find_neighbor_list(sac_ppl(iac)%neighbor_list_set,atom=iatom)
                   
                   CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,nnode=nneighbor)
                   
                   ALLOCATE (vppl(kkind)%r2(nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE (vppl(kkind)%neighbor(nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE (vppl(kkind)%r(3,nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   
                   sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)
                   
                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                                             neighbor=vppl(kkind)%neighbor(kneighbor),&
                                             r=vppl(kkind)%r(:,kneighbor))
                      vppl(kkind)%r2(kneighbor) =&
                           vppl(kkind)%r(1,kneighbor)*vppl(kkind)%r(1,kneighbor) +&
                           vppl(kkind)%r(2,kneighbor)*vppl(kkind)%r(2,kneighbor) +&
                           vppl(kkind)%r(3,kneighbor)*vppl(kkind)%r(3,kneighbor)
                      sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
                   END DO
                END DO
             END IF
             
             IF(all_potential_present) THEN
                DO kkind = 1,nkind
                   IF(ASSOCIATED(VH_3c_list(kkind,iatom)%neighbor)) &
                        vall(kkind)%neighbor =>  VH_3c_list(kkind,iatom)%neighbor
                END DO
             END IF
             
             last_jatom = 0
             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)
             
             DO inode=1,nnode
                CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)
                
                IF (jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                ELSE
                   new_atom_b = .FALSE.
                END IF
                
                IF(oce_present) THEN
                   CALL get_neighbor_node(sab_orb_neighbor_node,reduced_3c_oce=reduced_3c_oce)
                   
                   ! *** Build reduced 3-center lists
                   IF(.NOT. calculate_forces) &
                        CALL build_reduced_3c_lists(atomic_kind_set, iatom, jatom,&
                        ikind, jkind, sbc_oce, oce%voce, reduced_3c_oce,&
                        .TRUE.,qs_env%dft_control%qs_control%eps_rho_rspace,error)
                   
                   ! *** Calculate the coefficients of the one center expansions ***
                   CALL set_up_oce(atomic_kind_set,reduced_3c_oce,oce,ikind,iatom,jkind,jatom,&
                        rab,oceh,oces,nder,dist_list,error=error)
                END IF
                
                IF(rho0_present) THEN
                   CALL get_neighbor_node(sab_orb_neighbor_node,reduced_3c_rho0=reduced_3c_rho0)
                   ! *** Build reduced 3-center lists
                   IF(.NOT. calculate_forces) &
                        CALL build_reduced_3c_lists(atomic_kind_set, iatom, jatom, &
                        ikind, jkind, sbc_3c, VH_3c_list, reduced_3c_rho0,&
                        .FALSE.,qs_env%dft_control%qs_control%eps_rho_rspace,error)
                END IF
                
                IF (calculate_forces) atom_b = atom_of_kind(jatom)
                
                ! *** Use the symmetry of the first derivatives ***
                IF (iatom == jatom) THEN
                   f0 = 1.0_dp
                ELSE
                   f0 = 2.0_dp
                END IF
                
                ! *** Create matrix blocks for a new matrix block column ***
                IF (new_atom_b) THEN
                   IF (iatom <= jatom) THEN
                      irow = iatom
                      icol = jatom
                   ELSE
                      irow = jatom
                      icol = iatom
                   END IF
                   DO i=1,SIZE(matrix_s)
                      NULLIFY (sint(i)%block)
                      CALL add_block_node(matrix=matrix_s(i)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=sint(i)%block,error=error)
                   END DO
                   IF (build_kinetic_energy_matrix) THEN
                      DO i=1,SIZE(matrix_t)
                         NULLIFY (tint(i)%block)
                         CALL add_block_node(matrix=matrix_t(i)%matrix,&
                                             block_row=irow,&
                                             block_col=icol,&
                                             BLOCK=tint(i)%block,error=error)
                      END DO
                   END IF
                   NULLIFY (h_block)
                   CALL add_block_node(matrix=matrix_h(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=h_block,error=error)
                   IF (calculate_forces) THEN
                      CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=p_block)
                      IF (.NOT.ASSOCIATED(p_block)) THEN
                         CALL stop_program(routineN,moduleN,__LINE__,&
                              "A density matrix block is missing",para_env)
                      END IF
                      CALL get_block_node(matrix=matrix_w(1)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=w_block)
                      IF (.NOT.ASSOCIATED(w_block)) THEN
                         CALL stop_program(routineN,moduleN,__LINE__,&
                              "A weighted density matrix block is missing",para_env)
                      END IF
                   END IF
                END IF
                
                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)
                nrow = 1
                DO iset=1,nseta
                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)
                   DO jset=1,nsetb
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN
                         IF (calculate_forces.AND.((iatom /= jatom).OR.use_virial)) THEN
                            ! *** Decontract W matrix block ***
                            IF (iatom <= jatom) THEN
                               CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    w_block(sgfa,sgfb),SIZE(w_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                            ELSE
                               CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    w_block(sgfb,sgfa),SIZE(w_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                            END IF
                            CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                                 1.0_dp,work(1,1),SIZE(work,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,pab(1,1),SIZE(pab,1))
                            
                            ! *** Calculate the primitive overlap integrals ***
                            ! *** and the corresponding force contribution  ***
                            CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sab,nder,return_s_derivatives,&
                                 ai_work,ldai,sdab,pab,force_a)
                            force(ikind)%overlap(:,atom_a)=force(ikind)%overlap(:,atom_a) - 2.0_dp*force_a(:)
                            force(jkind)%overlap(:,atom_b)=force(jkind)%overlap(:,atom_b) + 2.0_dp*force_a(:)
                            IF (use_virial) THEN
                               CALL virial_pair_force ( virial%pv_virial, -f0, force_a, rab, error)
                            END IF
                         ELSE
                            ! *** Calculate the primitive overlap integrals ***
                            CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sab,nder,return_s_derivatives,&
                                 ai_work,ldai,sdab)
                         END IF
                         ! *** Contraction step (overlap matrix and its derivatives) ***
                         DO i=1,SIZE(matrix_s)
                            IF ((i > 1).AND.(.NOT.return_s_derivatives)) CYCLE
                            first_col = (i - 1)*SIZE(sab,1) + 1
                            CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                 1.0_dp,sab(1,first_col),SIZE(sab,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                            IF (iatom <= jatom) THEN
                               CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    work(1,1),SIZE(work,1),&
                                    1.0_dp,sint(i)%block(sgfa,sgfb),&
                                    SIZE(sint(i)%block,1))
                            ELSE
                               IF (i > 1) THEN
                                  f = -1.0_dp
                               ELSE
                                  f = 1.0_dp
                               END IF
                               CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                    f,work(1,1),SIZE(work,1),&
                                    sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    1.0_dp,sint(i)%block(sgfb,sgfa),&
                                    SIZE(sint(i)%block,1))
                            END IF
                         END DO
                         ! *** Decontract density matrix block ***
                         IF (calculate_forces) THEN
                            IF (iatom <= jatom) THEN
                               CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    p_block(sgfa,sgfb),SIZE(p_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                            ELSE
                               CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    p_block(sgfb,sgfa),SIZE(p_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                            END IF
                            CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                                 1.0_dp,work(1,1),SIZE(work,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,pab(1,1),SIZE(pab,1))
                         END IF
                         
                         IF (calculate_forces.AND.((iatom /= jatom).OR.use_virial)) THEN
                            ! *** Calculate the primitive kinetic energy integrals ***
                            ! *** and the corresponding force contribution         ***
                            CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sdab,hab,nder,return_t_derivatives,&
                                 ai_work,ldai,pab,force_a)
                            force(ikind)%kinetic(:,atom_a) =force(ikind)%kinetic(:,atom_a) + 2.0_dp*force_a(:)
                            force(jkind)%kinetic(:,atom_b) =force(jkind)%kinetic(:,atom_b) - 2.0_dp*force_a(:)
                            IF (use_virial) THEN
                               CALL virial_pair_force ( virial%pv_virial, f0, force_a, rab, error)
                            END IF
                         ELSE
                            ! *** Calculate the primitive kinetic energy integrals ***
                            CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sdab,hab,nder,return_t_derivatives,&
                                 ai_work,ldai)
                         END IF
                         
                         ! *** Contraction step (kinetic energy integral ***
                         ! *** matrix and its derivatives)               ***
                         IF (build_kinetic_energy_matrix) THEN
                            DO i=1,SIZE(matrix_t)
                               IF ((i > 1).AND.(.NOT.return_t_derivatives)) CYCLE
                               first_col = (i - 1)*SIZE(hab,1) + 1
                               CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                    1.0_dp,hab(1,first_col),SIZE(hab,1),&
                                    sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                               IF (iatom <= jatom) THEN
                                  CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                       1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                       work(1,1),SIZE(work,1),&
                                       1.0_dp,tint(i)%block(sgfa,sgfb),&
                                       SIZE(tint(i)%block,1))
                               ELSE
                                  IF (i > 1) THEN
                                     f = -1.0_dp
                                  ELSE
                                     f = 1.0_dp
                                  END IF
                                  CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                       f,work(1,1),SIZE(work,1),&
                                       sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                       1.0_dp,tint(i)%block(sgfb,sgfa),&
                                       SIZE(tint(i)%block,1))
                               END IF
                            END DO
                         END IF
                         
                         ! *** Loop over the GTH pseudo potential atoms (local part) ***
                         IF (ppl_present) THEN
                            DO kkind=1,nkind
                               IF (.NOT.ASSOCIATED(vppl(kkind)%neighbor)) CYCLE
                               DO kneighbor=1,SIZE(vppl(kkind)%neighbor)
                                  katom = vppl(kkind)%neighbor(kneighbor)
                                  rac(:) = vppl(kkind)%r(:,kneighbor)
                                  rac2 = vppl(kkind)%r2(kneighbor)
                                  dac = SQRT(rac2)
                                  
                                  IF (set_radius_a(iset) + vppl(kkind)%ppl_radius < dac) CYCLE
                                  rbc(:) = rac(:) - rab(:)
                                  rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                                  dbc = SQRT(rbc2)
                                  IF (set_radius_b(jset) + vppl(kkind)%ppl_radius < dbc) CYCLE
                                  
                                  ! *** Calculate the GTH pseudo potential forces ***
                                  IF (calculate_forces) THEN
                                     CALL overlap_ppl(&
                                          la_max(iset),la_min(iset),npgfa(iset),&
                                          rpgfa(:,iset),zeta(:,iset),&
                                          lb_max(jset),lb_min(jset),npgfb(jset),&
                                          rpgfb(:,jset),zetb(:,jset),&
                                          vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                                          vppl(kkind)%ppl_radius,&
                                          rab,dab,rac,dac,rbc,dbc,&
                                          hab,nder,nder,.FALSE.,ai_work,&
                                          pab,force_a,force_b)
                                     
                                     ! *** The derivatives w.r.t. atomic center c are    ***
                                     ! *** calculated using the translational invariance ***
                                     ! *** of the first derivatives                      ***
                                     atom_c = atom_of_kind(katom)
                                     force(ikind)%gth_ppl(:,atom_a) =force(ikind)%gth_ppl(:,atom_a) + f0*force_a(:)
                                     force(kkind)%gth_ppl(:,atom_c) =force(kkind)%gth_ppl(:,atom_c) - f0*force_a(:)
                                     
                                     force(jkind)%gth_ppl(:,atom_b) =force(jkind)%gth_ppl(:,atom_b) + f0*force_b(:)
                                     force(kkind)%gth_ppl(:,atom_c) =force(kkind)%gth_ppl(:,atom_c) - f0*force_b(:)
                                     
                                     IF (use_virial) THEN
                                        CALL virial_pair_force ( virial%pv_virial, f0, force_a, rac, error)
                                        CALL virial_pair_force ( virial%pv_virial, f0, force_b, rbc, error)
                                     END IF
                                  ELSE
                                     CALL overlap_ppl(&
                                          la_max(iset),la_min(iset),npgfa(iset),&
                                          rpgfa(:,iset),zeta(:,iset),&
                                          lb_max(jset),lb_min(jset),npgfb(jset),&
                                          rpgfb(:,jset),zetb(:,jset),&
                                          vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                                          vppl(kkind)%ppl_radius,&
                                          rab,dab,rac,dac,rbc,dbc,&
                                          hab,0,0,.FALSE.,ai_work)
                                  END IF
                               END DO
                            END DO
                         END IF
                      ELSE
                         DO jco=1,ncob
                            DO ico=1,ncoa
                               hab(ico,jco) = 0.0_dp
                            END DO
                         END DO
                         IF (calculate_forces) THEN
                            ! *** Decontract density matrix block ***
                            IF (iatom <= jatom) THEN
                               CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    p_block(sgfa,sgfb),SIZE(p_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                            ELSE
                               CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    p_block(sgfb,sgfa),SIZE(p_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                            END IF
                            CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                                 1.0_dp,work(1,1),SIZE(work,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,pab(1,1),SIZE(pab,1))
                         END IF
                      END IF

                      IF (all_potential_present ) THEN
                         IF(calculate_forces) THEN
                            CALL all_integrals(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,jkind,nkind,sbc_3c,&
                                 reduced_3c_rho0,vall,&
                                 hab,nder,pab,force_a,force_b,virial_ab,error=error)
                            force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) + f0*force_a(:)
                            force(jkind)%all_potential(:,atom_b) = force(jkind)%all_potential(:,atom_b) + f0*force_b(:)
                            IF (use_virial) THEN
                               virial%pv_virial = virial%pv_virial + f0*virial_ab
                            END IF
                            DO kkind = 1,nkind
                               IF (.NOT. ASSOCIATED(vall(kkind)%force_c)) CYCLE
                               nnode_c = reduced_3c_rho0(kkind)%nnode
                               DO knode = 1,nnode_c
                                  katom = reduced_3c_rho0(kkind)%index_atom(knode)
                                  atom_c = atom_of_kind(katom)
                                  force(kkind)%all_potential(1:3,atom_c) = force(kkind)%all_potential(1:3,atom_c) + &
                                                                           f0*vall(kkind)%force_c(1:3,knode)
                               END DO
                               DEALLOCATE(vall(kkind)%force_c,STAT=stat)
                               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                            END DO
                         ELSE
                            CALL all_integrals(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,jkind,nkind,sbc_3c,&
                                 reduced_3c_rho0,vall,&
                                 hab,nder,error=error)
                            DO kkind = 1,nkind
                               IF (.NOT. ASSOCIATED(vall(kkind)%force_c)) CYCLE
                               DEALLOCATE(vall(kkind)%force_c,STAT=stat)
                               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                            END DO
                         END IF
                      END IF
                      
                      ! *** Contraction step (core Hamiltonian matrix) ***
                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                      IF (iatom <= jatom) THEN
                         CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              work(1,1),SIZE(work,1),&
                              1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
                      ELSE
                         CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                              1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
                      END IF
                   END DO
                   nrow = nrow + ncoa
                END DO
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             END DO
             IF (ppl_present) THEN
                DO kkind=1,nkind
                   IF (ASSOCIATED(vppl(kkind)%r2)) THEN
                      DEALLOCATE (vppl(kkind)%r2,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      DEALLOCATE (vppl(kkind)%neighbor,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      DEALLOCATE (vppl(kkind)%r,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END DO
             END IF
             
             IF (all_potential_present) THEN
                DO kkind=1,nkind
                   IF (ASSOCIATED(vall(kkind)%neighbor)) THEN
                      NULLIFY (vall(kkind)%neighbor)
                      NULLIFY (vall(kkind)%force_c)
                   END IF
                END DO
             END IF
          END DO
          
          IF (ALLOCATED(ai_work)) THEN
             DEALLOCATE (ai_work,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END DO
    END DO

    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (pabt,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! *** If LSD, then recover alpha density and beta density     ***
       ! *** from the total density (1) and the spin density (2)     ***
       ! *** The W matrix is neglected, since it will be destroyed   ***
       ! *** in the calling force routine after leaving this routine ***
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 0.5_dp,matrix_p(2)%matrix,0.5_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-1.0_dp,matrix_p(1)%matrix,1.0_dp,error=error)
       END IF
    END IF

    ! *** compute the ppnl contribution to the core hamiltonian ***
    CALL build_core_ppnl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    atomic_kind_set, particle_set, sab_orb, sac_ppnl, sbc_ppnl, error)

    ! *** Put the core Hamiltonian matrix in the QS environment ***
    CALL set_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,&
                    kinetic=matrix_t,&
                    matrix_h=matrix_h,&
                    oce=oce,error=error)

    ! *** Release work storage ***
    IF (ppl_present) THEN
       DEALLOCATE (vpplt,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    IF (all_potential_present) THEN
       DEALLOCATE (vallt,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    IF(ASSOCIATED(oce))THEN
       CALL deallocate_vtriple(oce%voce,error=error)
    END IF
    IF(ASSOCIATED(VH_3c_list))THEN
       CALL deallocate_vtriple(VH_3c_list,error=error)
    END IF
    
    DEALLOCATE (habt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sdabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (workt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    
    DEALLOCATE (sintt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF (build_kinetic_energy_matrix) THEN
       DEALLOCATE (tintt,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
       DEALLOCATE (oceht,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (ocest,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i=1,SIZE(dist_list,1)
          IF(ASSOCIATED(dist_list(i)%rdist)) THEN
              DEALLOCATE(dist_list(i)%rdist,STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
       END DO
       DEALLOCATE(dist_list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    
    ! Print matrices if requested
    CALL dump_info_core_hamiltonian(matrix_s, matrix_t, matrix_h, oce, qs_env, para_env, error)
    CALL timestop(handle)
    
  END SUBROUTINE build_core_hamiltonian_matrix
  
! *****************************************************************************
!> \brief Possibly prints matrices after the construction of the Core 
!>     Hamiltonian Matrix
! *****************************************************************************
  SUBROUTINE dump_info_core_hamiltonian(matrix_s, matrix_t, matrix_h, oce, qs_env, para_env, error)
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_t, matrix_h
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_info_core_hamiltonian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_v

    CALL timeset(routineN,handle)

    NULLIFY(logger, matrix_v)
    logger => cp_error_get_logger(error)

    ! *** Print the distribution of the overlap matrix blocks ***
    CALL write_matrix_distribution(matrix_s(1)%matrix,qs_env%input,para_env,error)

    ! *** Print the overlap integral matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",&
            extension=".Log",error=error)
      CALL write_sparse_matrix(matrix_s(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
      IF (BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
           "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error),cp_p_file)) THEN
         DO i=2,SIZE(matrix_s)
            CALL write_sparse_matrix(matrix_s(i)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
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
       CALL write_sparse_matrix(matrix_t(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error),cp_p_file)) THEN
          DO i=2,SIZE(matrix_t)
             CALL write_sparse_matrix(matrix_t(i)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
          END DO
       END IF
      CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/KINETIC_ENERGY", error=error)
    END IF

    ! *** Print the potential energy matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY",&
            extension=".Log",error=error)
      CALL allocate_matrix_set(matrix_v,1,error=error)
      CALL replicate_matrix_structure(source=matrix_h(1)%matrix,&
                                      TARGET=matrix_v(1)%matrix,&
                                      target_name="POTENTIAL ENERGY MATRIX",error=error)
      CALL add_matrices(matrix_v(1)%matrix,1.0_dp,matrix_h(1)%matrix,-1.0_dp,matrix_t(1)%matrix)
      CALL write_sparse_matrix(matrix_v(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
      CALL deallocate_matrix_set(matrix_v,error=error)
      CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/POTENTIAL_ENERGY", error=error)
    END IF

    ! *** Print the core Hamiltonian matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", error=error)
    END IF

    ! *** Print the OCE HARD matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/OCE_HARD",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OCE_HARD",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(oce%hardo(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/OCE_HARD", error=error)
    END IF

    ! *** Print the OCE SOFT matrix, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/OCE_SOFT",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OCE_SOFT",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(oce%softo(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/OCE_SOFT", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE dump_info_core_hamiltonian

END MODULE qs_core_hamiltonian
