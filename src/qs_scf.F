!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_scf [1.0] *
!!
!!   NAME
!!     qs_scf
!!
!!   FUNCTION
!!     Routines for the Quickstep SCF run.
!!
!!   AUTHOR
!!     Matthias Krack (30.04.2001)
!!
!!   MODIFICATION HISTORY
!!     - Joost VandeVondele (02.2002)
!!          added code for: incremental (pab and gvg) update
!!                           initialisation (init_cube, l_info)
!!     - Joost VandeVondele (02.2002)
!!          called the poisson code of the classical part
!!          this takes into account the spherical cutoff and allows for
!!          isolated systems
!!     - Joost VandeVondele (02.2002)
!!          added multiple grid feature
!!          changed to spherical cutoff consistently (?)
!!          therefore removed the gradient correct functionals
!!     - updated with the new QS data structures (10.04.02,MK)
!!     - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!!     - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!!     - set_mo_occupation for smearing of the MO occupation numbers
!!       (17.04.02,MK)
!!     - MO level shifting added (22.04.02,MK)
!!     - Usage of TYPE mo_set_p_type
!!     - Joost VandeVondele (05.2002)
!!           added cholesky based diagonalisation
!!     - 05.2002 added pao method [fawzi]
!!     - parallel FFT (JGH 22.05.2002)
!!     - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE basis_set_types,          ONLY: get_gto_basis_set,&
                                      gto_basis_set_type
  USE coefficient_types,        ONLY: coeff_allocate,&
                                      coeff_copy,&
                                      coeff_transform_space,&
                                      coeff_type,&
                                      coeff_zero,&
                                      coeff_sumup
  USE qs_core_hamiltonian,      ONLY: build_core_hamiltonian_matrix
  USE cube_utils,               ONLY: cube_info_type,&
                                      destroy_cube_info,&
                                      init_cube_info
  USE dft_types,                ONLY: dft_control_type
  USE external_potential_types, ONLY: all_potential_type,&
                                      get_potential,&
                                      gth_potential_type
  USE global_types,             ONLY: DEBUG,&
                                      global_environment_type
  USE hartree,                  ONLY: calculate_hartree
  USE l_utils,                  ONLY: destroy_l_info,&
                                      init_l_info,&
                                      l_info_type
  USE mathconstants,            ONLY: fourpi,twopi
  USE memory_utilities,         ONLY: reallocate
  USE message_passing,          ONLY: mp_max,mp_sum,mp_sync,&
                                      mp_range,mp_allgather
  USE particle_types,           ONLY: particle_type
  USE qs_arpack,                ONLY: arpack_diag
  USE qs_blacs,                 ONLY: allocate_blacs_matrix,&
                                      blacs_add_to_element,&
                                      blacs_gemm,&
                                      blacs_get_element,&
                                      blacs_matrix_type,&
                                      blacs_set_all,&
                                      blacs_set_element,&
                                      blacs_syevx,&
                                      blacs_symm,&
                                      copy_blacs_to_blacs_matrix,&
                                      copy_blacs_to_sparse_matrix,&
                                      copy_sparse_to_blacs_matrix,&
                                      finish_blacs,&
                                      get_blacs_matrix_info,&
                                      power_blacs_matrix,&
                                      replicate_blacs_matrix,&
                                      start_blacs,&
                                      symmetrise_blacs_matrix, &
                                      blacs_cholesky_decompose, &
                                      blacs_cholesky_reduce, &
                                      blacs_cholesky_restore
  USE qs_diis,                  ONLY: eps_diis,max_diis,scf_diis
  USE qs_environment_types,     ONLY: get_qs_env,&
                                      qs_environment_type,&
                                      set_qs_env
  USE qs_force_types,           ONLY: qs_force_type,&
                                      write_qs_force
  USE qs_integrate_potential,   ONLY: integrate_v_rspace, &
                                      integrate_v_core_rspace
  USE qs_mo_types,              ONLY: allocate_mo_set,&
                                      calculate_density_matrix,&
                                      calculate_w_matrix,&
                                      correct_mo_eigenvalues,&
                                      get_mo_set,&
                                      mo_set_p_type,&
                                      mo_set_type,&
                                      read_mo_set,&
                                      set_mo_occupation,&
                                      write_mo_set
  USE qs_overlap,               ONLY: write_blacs_matrix,&
                                      write_sparse_matrix
  USE qs_parser,                ONLY: finish_parser,&
                                      read_object,&
                                      start_parser,&
                                      test_object
  USE simulation_cell,          ONLY: cell_type,&
                                      get_cell
  USE sparse_matrix_types,      ONLY: add_matrices,&
                                      allocate_matrix,&
                                      deallocate_matrix,&
                                      first_block_node,&
                                      get_block_node,&
                                      get_matrix_info,&
                                      next_block_node,&
                                      real_block_node_type,&
                                      real_matrix_p_type,&
                                      replicate_matrix_structure,&
                                      set_matrix,&
                                      symmetrise_diagonal_blocks,&
                                      transfer_matrix
  USE termination,              ONLY: stop_memory,&
                                      stop_program
  USE timings,                  ONLY: timeset,&
                                      timestop
  USE timesl,                   ONLY: cputime
  use pao_obj_function,         only: pao_gerd_functional_type,&
       pao_gerd_f_init, pao_gerd_f_dealloc_ref
  use cp_lbfgs_optimizer_gerd_f,only: cp_lbfgs_opt_gerd_f_type, &
       cp_opt_gerd_f_next,&
       cp_opt_gerd_f_init, cp_opt_gerd_f_dealloc_ref
  use pao_types,                only: pao_env_type, pao_env_get, &
       pao_glob_angles_type
  use cp_sparse_matrix,         only: cp_sparse_matrix_type, &
       cp_sparse_matrix_p_type
  use pao_proj_methods,         only: sm_matrix_transf_to_f,&
       sm_matrix_transf_to_m
  use pao_qs_env_methods,       only: pao_qs_env_did_change
  use cp_error_handling,        only: cp_error_type, cp_assert, &
       cp_error_message, cp_error_init, cp_error_dealloc_ref, cp_debug,&
       cp_internal_error
  use cp_log_handling,          only: cp_failure_level, cp_warning_level,&
       cp_to_string
  use pao_glob_angles_methods,  only: pao_g_ang_get,pao_g_ang_set
  use qs_build_KS_matrix,       only: qs_ks_env_type, qs_ks_init, &
       qs_ks_dealloc_ref, build_kohn_sham_matrix
  USE qs_collocate_density,     ONLY: calculate_rho_core,&
                                      calculate_rho_elec,&
                                      calculate_total_rho
  USE qs_core_energies,         ONLY: calculate_ecore,&
                                      calculate_ecore_overlap,&
                                      calculate_ecore_self

  IMPLICIT NONE

  PRIVATE

  TYPE(blacs_matrix_type), POINTER  :: ortho,scf_work1,scf_work2
  TYPE(cell_type), POINTER          :: cell
  TYPE(dft_control_type), POINTER   :: dft_control
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s

  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE(qs_force_type), DIMENSION(:), POINTER    :: force

  TYPE(real_matrix_p_type) :: w

  CHARACTER(LEN=10) :: density_guess
  REAL(wp)          :: eps_eigval,eps_scf,level_shift,p_mix,smear,&
       work_syevx
!!$       etotal,ex,total_rho_core_rspace,total_rho_elec_rspace,&
!!$                       total_rho_gspace,total_rho_rspace
  INTEGER           :: charge,max_scf,nelectron,nrebuild,nrow_block,&
                       ncol_block,nprow,npcol,maxl,nkind
  LOGICAL           :: gradient_functional,gth_potential_present,&
                       rebuild,use_cholesky,use_arpack

  character(len=*), private, parameter :: moduleN="qs_scf"

! *** Public variables ***

  PUBLIC :: nelectron

! *** Public subroutines ***

  PUBLIC :: read_scf_parameters,&
            scf,&
            write_scf_parameters

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE scf(qs_env,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE scf (MODULE qs_scf)"

!   *** Local variables ***

    TYPE(mo_set_type), POINTER :: mo_set

    REAL(wp) :: delta,diis_error,t1,t2
    INTEGER  :: context,handle,homo,igrid_level,imo,iscf,istat,nmo,&
                output_unit
    LOGICAL  :: calculate_forces,diis_step,do_level_shift,ionode

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors

    logical :: pao, failure
    type(pao_gerd_functional_type), pointer :: pao_obj_f
    type(cp_lbfgs_opt_gerd_f_type), pointer :: pao_optimizer
    type(cp_error_type) :: error
    type(pao_env_type),pointer :: pao_env
    type(cp_sparse_matrix_p_type),dimension(:),pointer :: min_density_m
    real(kind=wp)  ::pao_f,pao_last_f,pao_actual_grad
    
    type(qs_ks_env_type) :: ks_env
    type(real_matrix_p_type), pointer :: ks, p
    type(real_matrix_p_type), dimension(:), pointer :: ks_rmpv, p_rmpv,&
         p_old_rmpv
    real(kind=wp) :: ex,ec,ecore,ecore_overlap,&
         ecore_self,ehartree,&
         total_rho_gspace, total_rho_elec_rspace, etotal, &
         total_rho_core_rspace, total_rho_rspace
   

!   ---------------------------------------------------------------------------



!   *** Quick return, if no SCF iteration is requested ***

    IF (max_scf < 1) RETURN

    CALL timeset("scf","I","",handle)

    call cp_error_init(error)

    ionode = globenv%ionode
    output_unit = globenv%scr
    
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    h=h,&
                    s=s,&
                    particle_set=particle_set,&
                    pao_env=pao_env)

    calculate_forces = dft_control%forces
    pao=dft_control%qs_control%pao.and.associated(pao_env)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxl=maxl,&
                             nelectron=nelectron)

    nelectron = nelectron - charge
    call set_qs_env(qs_env, n_electrons=nelectron)

    CALL get_cell(cell=cell)

    nkind = SIZE(atomic_kind_set)

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "SCF WAVEFUNCTION OPTIMIZATION"
    END IF

    CALL init_scf_run(mo_set,context,qs_env,globenv,ecore_self=ecore_self, &
       ecore_overlap=ecore_overlap,calculate_forces=calculate_forces)

    call qs_ks_init(ks_env,qs_env=qs_env,global_env=globenv, error=error)
    ks_env%ecore_self=ecore_self ! ugly
    ks_env%ecore_overlap=ecore_overlap ! ugly

    call pao_opt_create(pao_optimizer,pao_obj_f,qs_env,globenv,error)

    CALL get_qs_env(qs_env,p=p_rmpv, k=ks_rmpv, p_old=p_old_rmpv)
    p => p_rmpv(1)
    ks => ks_rmpv(1)

    !!FM not nice! mv to qs_build_KS_matrix?
    CALL calculate_rho_core(ks_env%rho_core,ks_env%cube_info(1),&
         ks_env%l_info,&
         total_rho_core_rspace,qs_env,globenv)
    ks_env%total_rho_core_rspace=total_rho_core_rspace
    
    etotal = 0.0_wp
    iscf = 0
    diis_step = .FALSE.
    rebuild = .TRUE.

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,T3,A,T9,A,T34,A,T49,A,T68,A,/,T3,A)")&
        "Step","Update method","Time","Convergence","Total energy",&
        REPEAT("-",77)
    END IF

!   *** SCF loop ***

    scf_loop: DO

      iscf = iscf + 1

      t1 = cputime()

      IF (globenv%print%density_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_sparse_matrix(p%matrix,4,6,qs_env,globenv)
      END IF

      call build_kohn_sham_matrix(ks_env,qs_env=qs_env,globenv=globenv,&
           ks=ks_rmpv,p=p_rmpv,p_old=p_old_rmpv,rebuild=rebuild,&
           nrebuild=nrebuild,ex=ex,ec=ec,ecore=ecore,&
           ecore_overlap=ecore_overlap,&
           ecore_self=ecore_self,ehartree=ehartree,&
           total_rho_gspace=total_rho_gspace, &
           total_rho_elec_rspace=total_rho_elec_rspace, etotal=etotal, &
           total_rho_core_rspace=total_rho_core_rspace, &
           total_rho_rspace=total_rho_rspace,error=error,&
           calculate_forces=.FALSE.)

      rebuild = (MODULO(iscf,nrebuild) == 0)

      paoIf: if (pao) then
         call pao_qs_env_did_change(qs_env,global_env=globenv,h_changed=.true.)
         if (.not.cp_opt_gerd_f_next(pao_optimizer,f=pao_f,last_f=pao_last_f,&
              projected_gradient=pao_actual_grad, error=error)) then
            call cp_assert(.false.,cp_warning_level,&
                 error_nr=cp_internal_error,fromWhere=routine,&
                 message="pao optimizer reached convergence",&
                 error=error)
         end if
         call pao_env_get(pao_env, min_density_m=min_density_m,error=error)
         call sm_matrix_transf_to_f(source_m=min_density_m(1)%matrix%matrix,&
              target_m=ks%matrix, glob_angles=pao_env%angles_att,&
              qs_env=qs_env, global_env=globenv,error=error)
      else paoIf
      IF (globenv%print%kohn_sham_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_sparse_matrix(ks%matrix,4,6,qs_env,globenv)
      END IF

      CALL copy_sparse_to_blacs_matrix(ks%matrix,scf_work1,context,globenv)

      IF (iscf > 1) THEN
        CALL scf_diis(mo_set,scf_work1,scf_work2,delta,diis_error,diis_step,&
                      context,qs_env,globenv)
      END IF

      do_level_shift = ((level_shift /= 0.0_wp).and.&
           ((density_guess == "CORE").or.(iscf > 1)))
      CALL eigensolver(scf_work1,mo_set,ortho,scf_work2, &
                       do_level_shift,level_shift,context,globenv)
      
!     *** Update the MO occupation numbers, if smearing is used ***

      IF (smear /= 0.0_wp) CALL set_mo_occupation(mo_set,smear)

      IF (globenv%print%each_scf_step) THEN
        CALL write_mo_set(mo_set,atomic_kind_set,particle_set,4,6,context,&
                          globenv)
      END IF

      CALL calculate_density_matrix(mo_set,scf_work1,context,globenv)

      CALL copy_blacs_to_sparse_matrix(scf_work1,ks%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(ks%matrix)
      
   end if paoIf

      t2 = cputime()

      IF (diis_step) THEN
        CALL density_mixing(ks,p,1.0_wp,delta,globenv)
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T15,E10.2,T30,F8.2,T40,2F20.10)")&
            iscf,"DIIS",diis_error,t2 - t1,delta,etotal
        END IF
      ELSE
        CALL density_mixing(ks,p,p_mix,delta,globenv)
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T15,F6.2,T30,F8.2,T40,2F20.10)")&
            iscf,"Mixing",p_mix,t2 - t1,delta,etotal
        END IF
      END IF

      IF (delta < eps_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run converged ***"
        END IF
        EXIT scf_loop
      ELSE IF (iscf == max_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run NOT converged ***"
        END IF
        EXIT scf_loop
      END IF

    END DO scf_loop

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
        "Total electronic density (r-space): ",&
        total_rho_elec_rspace,total_rho_elec_rspace + REAL(nelectron,wp),&
        "Total core charge density (r-space):",&
        total_rho_core_rspace,total_rho_core_rspace - REAL(nelectron,wp)
      WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
        "Total charge density (r-space):     ",total_rho_rspace,&
        "Total charge density (g-space):     ",total_rho_gspace
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self,&
        "Core Hamiltonian energy:                       ",ecore,&
        "Hartree energy:                                ",ehartree,&
        "Exchange-correlation energy:                   ",ex + ec,&
        "Total energy:                                  ",etotal
    END IF

!   *** Undo mixing of the density matrix (restore original density) ***

    IF (.NOT.diis_step) CALL transfer_matrix(ks%matrix,p%matrix)

    CALL write_mo_set(mo_set,atomic_kind_set,particle_set,4,6,context,globenv)

    IF (globenv%print%density_matrix) THEN
      CALL write_sparse_matrix(p%matrix,4,6,qs_env,globenv)
    END IF

    IF (globenv%print%kohn_sham_matrix) THEN
      call build_kohn_sham_matrix(ks_env,qs_env=qs_env,globenv=globenv,&
           ks=ks_rmpv,p=p_rmpv,p_old=p_old_rmpv,rebuild=rebuild,&
           nrebuild=nrebuild,ex=ex,ec=ec,ecore=ecore,&
           ecore_overlap=ecore_overlap,&
           ecore_self=ecore_self,ehartree=ehartree,&
           total_rho_gspace=total_rho_gspace, &
           total_rho_elec_rspace=total_rho_elec_rspace, etotal=etotal, &
           total_rho_core_rspace=total_rho_core_rspace, &
           total_rho_rspace=total_rho_rspace,error=error,&
           calculate_forces=.FALSE.)
      CALL write_sparse_matrix(ks%matrix,4,6,qs_env,globenv)
    END IF

!   *** Calculate forces ***

    IF (calculate_forces) THEN
      NULLIFY (w%matrix)
      CALL replicate_matrix_structure(h(1)%matrix,w%matrix,"W MATRIX")
      CALL calculate_w_matrix(mo_set,scf_work1,context,globenv)
      CALL copy_blacs_to_sparse_matrix(scf_work1,w%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(w%matrix)
      IF (globenv%ionode.AND.globenv%print%w_matrix) THEN
        CALL write_sparse_matrix(w%matrix,4,6,qs_env,globenv)
      END IF
      CALL build_core_hamiltonian_matrix(qs_env,globenv,.TRUE.,p%matrix,&
                                         w%matrix)
      CALL build_kohn_sham_matrix(ks_env=ks_env,&
                                  qs_env=qs_env,&
                                  globenv=globenv,&
                                  ks=ks_rmpv,&
                                  p=p_rmpv,&
                                  p_old=p_old_rmpv,&
                                  rebuild=rebuild,&
                                  nrebuild=nrebuild,&
                                  ex=ex,&
                                  ec=ec,&
                                  ecore=ecore,&
                                  ecore_overlap=ecore_overlap,&
                                  ecore_self=ecore_self,&
                                  ehartree=ehartree,&
                                  total_rho_gspace=total_rho_gspace,&
                                  total_rho_elec_rspace=total_rho_elec_rspace,&
                                  etotal=etotal,&
                                  total_rho_core_rspace=total_rho_core_rspace,&
                                  total_rho_rspace=total_rho_rspace,&
                                  error=error,&
                                  calculate_forces=.TRUE.)
      CALL get_qs_env(qs_env=qs_env,force=force)
      CALL write_qs_force(force,globenv)
    END IF

!   *** Write restart file ***

    CALL write_mo_set(mo_set,context,globenv)

    call pao_opt_dealloc(pao_optimizer,pao_obj_f,qs_env,globenv,error)
    
    call qs_ks_dealloc_ref(ks_env, qs_env=qs_env, global_env=globenv, error=error)

    CALL finish_blacs(context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf

! *****************************************************************************

  SUBROUTINE init_scf_run(mo_set,context,qs_env,globenv, ecore_self, &
       ecore_overlap,calculate_forces)

!   Purpose: Initialise a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    INTEGER, INTENT(OUT)                      :: context
    LOGICAL, INTENT(IN)                       :: calculate_forces
    real(kind=wp), intent(out)                :: ecore_self, ecore_overlap

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_scf_run (MODULE qs_scf)"

!   *** Local variables ***

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors

    INTEGER :: handle,homo,ikind,istat,nao,ndep,nmo,output_unit
    LOGICAL :: ionode, pao

    INTEGER, DIMENSION(0:globenv%num_pe-1) :: nblock_pe,nelement_pe
    type(real_matrix_p_type), dimension(:), pointer :: k_rmpv, p_rmpv, &
         p_old_rmpv

!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Initialise BLACS ***

    CALL start_blacs(nprow,npcol,context,globenv)

!   *** Get the dimension of the full SCF matrices, ***
!   *** i.e. the total number of atomic orbitals    ***

    CALL get_matrix_info(matrix=h(1)%matrix,nrow=nao)

    IF (globenv%print%mo_eigenvalues.OR.&
        (level_shift /= 0.0_wp).OR.&
        (smear /= 0.0_wp)) THEN
      nmo = nao
    ELSE
      nmo = nelectron/2
    END IF

!   *** Allocate a MO set ***

    ALLOCATE (c(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"c",int_size)
    NULLIFY (c(1)%mo_set)

!   *** Allocate the distributed MO eigenvectors ***

    CALL allocate_mo_set(c(1)%mo_set,nao,nmo,nelectron,2.0_wp,nrow_block,&
                         ncol_block,context,globenv)

!   *** Put the MO set in the QS environment ***

    CALL set_qs_env(qs_env=qs_env,c=c)

!   *** Get the HOMO and the MO eigenvectors ***

    mo_set => c(1)%mo_set

    CALL get_mo_set(mo_set=mo_set,&
                    homo=homo,&
                    eigenvectors=mo_eigenvectors)

!   *** Get BLACS block size of the MO eigenvector matrix      ***
!   *** which has to fit to the other distributed SCF matrices ***

    CALL get_blacs_matrix_info(matrix=mo_eigenvectors,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

!   *** Allocate the distributed SCF matrices ***

    CALL allocate_blacs_matrix(new_matrix=ortho,&
                               nrow_global=nao,&
                               ncol_global=nao,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="ORTHOGONALIZATION MATRIX",&
                               context=context,&
                               globenv=globenv)

    CALL replicate_blacs_matrix(prototype_matrix=ortho,&
                                new_matrix=scf_work1,&
                                name="SCF WORK MATRIX 1")

    CALL replicate_blacs_matrix(prototype_matrix=ortho,&
                                new_matrix=scf_work2,&
                                name="SCF WORK MATRIX 2")

    CALL copy_sparse_to_blacs_matrix(h(1)%matrix,scf_work1,context,globenv)

!   *** Redistribute the core Hamiltonian matrix in ***
!   *** order to eliminate redundant atomic blocks  ***

    IF ((globenv%num_pe > 1).AND.gth_potential_present) THEN
      CALL deallocate_matrix(h(1)%matrix)
      CALL replicate_matrix_structure(s(1)%matrix,h(1)%matrix,&
                                      "CORE HAMILTONIAN MATRIX")
      CALL copy_blacs_to_sparse_matrix(scf_work1,h(1)%matrix,context,globenv)
      call set_qs_env(qs_env,h=h)
    END IF

    ! allocate k and p  and puth them in the qs_environment
    allocate(k_rmpv(1), p_rmpv(1), p_old_rmpv(1))
    call set_qs_env(qs_env, p=p_rmpv, k=k_rmpv, p_old=p_old_rmpv)
  
    NULLIFY (k_rmpv(1)%matrix)
    CALL replicate_matrix_structure(h(1)%matrix,k_rmpv(1)%matrix,"KOHN-SHAM MATRIX")

    NULLIFY (p_rmpv(1)%matrix)
    CALL replicate_matrix_structure(h(1)%matrix,p_rmpv(1)%matrix,"DENSITY MATRIX")

    NULLIFY (p_old_rmpv(1)%matrix)
    CALL replicate_matrix_structure(h(1)%matrix,p_old_rmpv(1)%matrix,"OLD DENSITY MATRIX")

    CALL calculate_ecore_self(ecore_self,qs_env)
    CALL calculate_ecore_overlap(ecore_overlap,qs_env,globenv,calculate_forces)

    IF (ionode.AND.globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self
    END IF

    if (use_cholesky .or. use_arpack) then
       ndep=0
       CALL copy_sparse_to_blacs_matrix(s(1)%matrix,ortho,context,globenv)
       CALL blacs_cholesky_decompose(ortho,context,globenv)
    else
       CALL calculate_ortho_matrix(ortho,scf_work2,ndep,context,globenv)
    endif

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,I10))")&
        "Number of electrons:                    ",nelectron,&
        "Number of occupied orbitals:            ",homo,&
        "Number of orbital functions:            ",nao,&
        "Number of independent orbital functions:",nao - ndep
    END IF

    IF (globenv%print%ortho_matrix) THEN
      CALL write_blacs_matrix(ortho,4,6,context,qs_env,globenv)
    END IF

    CALL calculate_first_density_matrix(ortho,mo_set,p_rmpv(1),scf_work1,scf_work2,&
                                        context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,work,ndep,context,globenv)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(blacs_matrix_type), POINTER          :: ortho,work
    INTEGER, INTENT(IN)                       :: context
    INTEGER, INTENT(OUT)                      :: ndep

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL copy_sparse_to_blacs_matrix(s(1)%matrix,ortho,context,globenv)
    CALL power_blacs_matrix(ortho,work,-0.5_wp,eps_eigval,ndep,work_syevx,&
                            context,globenv)
    CALL symmetrise_blacs_matrix(ortho,work,context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ortho_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work, &
    do_level_shift,level_shift,context,globenv)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.
!   ks will be modified

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    INTEGER, INTENT(IN)                       :: context
    LOGICAL, INTENT(IN)                       :: do_level_shift
    REAL(WP), INTENT(IN)                      :: level_shift

!   *** Local variables ***

    INTEGER :: handle,imo,nao,nmo,homo,ncv,handle2

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors
    REAL(wp), DIMENSION(:), POINTER  :: mo_eigenvalues

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    homo=homo, &
                    eigenvalues=mo_eigenvalues,&
                    eigenvectors=mo_eigenvectors)

     IF (do_level_shift) THEN
        IF (use_cholesky .or. use_arpack) THEN
           call stop_program("eigensolver","level shift not implemented")
        ENDIF
     ENDIF


    IF (use_cholesky) THEN

       CALL blacs_cholesky_reduce(ks,ortho,context,globenv)
       CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,work_syevx,context,globenv)
       CALL blacs_cholesky_restore(work,nmo,ortho,mo_eigenvectors,"SOLVE",context,globenv)

    ELSEIF (use_arpack) THEN ! not really recommended, since rather slow, except for huge basis sets

       ncv = nmo + 64
       CALL blacs_cholesky_reduce(ks,ortho,context,globenv)
       CALL arpack_diag(ks,work,mo_eigenvalues,nmo,ncv,.false.,context,globenv)
       CALL blacs_cholesky_restore(work,nmo,ortho,mo_eigenvectors,"SOLVE",context,globenv)

    ELSE

       CALL get_blacs_matrix_info(matrix=ks,nrow_global=nao)
       CALL blacs_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp, & 
                        work,context,globenv)
       CALL blacs_gemm("T","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp, &
                        ks,context,globenv)
 
       IF (do_level_shift) THEN
         DO imo=homo+1,nmo
            CALL blacs_add_to_element(ks,imo,imo,level_shift,context,&
                                    globenv)
         END DO
       ENDIF

       CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,work_syevx,context,globenv)
       CALL blacs_gemm("N","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
                        mo_eigenvectors,context,globenv)

       IF (do_level_shift) THEN
          CALL correct_mo_eigenvalues(mo_set,level_shift)
          CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
       END IF

    ENDIF
    ! write(6,*) "evals syevx",mo_eigenvalues(1:nmo)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************

  SUBROUTINE density_mixing(new_density,old_density,p_mix,delta,globenv)

!   Purpose: Perform a density mixing of the old (last SCF iteration) and the
!            new (current) density matrix.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_p_type), INTENT(INOUT)   :: new_density,old_density
    REAL(wp), INTENT(IN)                      :: p_mix
    REAL(wp), INTENT(OUT)                     :: delta

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    REAL(wp) :: r_mix
    INTEGER  :: handle,i,iblock_col,iblock_row,j,nblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: p_new_block,p_old_block

!   ---------------------------------------------------------------------------

    CALL timeset("density_mixing","I","",handle)

    delta = 0.0_wp

    r_mix = 1.0_wp - p_mix

    CALL get_matrix_info(matrix=old_density%matrix,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=new_density%matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=p_new_block)

        CALL get_block_node(matrix=old_density%matrix,&
                            block_row=iblock_row,&
                            block_col=iblock_col,&
                            block=p_old_block)

        DO j=1,SIZE(p_new_block,2)
          DO i=1,SIZE(p_new_block,1)
            p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
            delta = MAX(delta,ABS(p_new_block(i,j)))
            p_old_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
            p_new_block(i,j) = p_old_block(i,j) + r_mix*p_new_block(i,j)
          END DO
        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

    CALL mp_max(delta,globenv%group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE density_mixing

! *****************************************************************************

  SUBROUTINE calculate_first_density_matrix(ortho,mo_set,p,work1,work2,&
                                            context,globenv)

    TYPE(blacs_matrix_type), POINTER          :: ortho,work1,work2
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(real_matrix_p_type), INTENT(INOUT)   :: p
    INTEGER, INTENT(IN)                       :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_first_density_matrix (MODULE qs_scf)"

!   *** Local variables ***

    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(blacs_matrix_type), POINTER  :: mo_eigenvectors
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    REAL(wp) :: paa,saa,sum
    INTEGER  :: atom_a,first_sgfa,handle,iatom,ikind,iset,isgf,isgfa,ishell,&
                istat,la,last_sgfa,maxl,natom,nelec,nset,restart_unit

    INTEGER, DIMENSION(:), ALLOCATABLE :: econf

    INTEGER, DIMENSION(:), POINTER   :: atom_list,elec_conf,nshell
    INTEGER, DIMENSION(:,:), POINTER :: first_sgf,l,last_sgf

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_first_density_matrix","I","",handle)

    IF (density_guess == "RESTART") THEN

      CALL read_mo_set(mo_set,context,globenv)
      IF (level_shift /= 0.0_wp) THEN
        CALL get_mo_set(mo_set=mo_set,&
                        eigenvectors=mo_eigenvectors)
        CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
      END IF
      CALL calculate_density_matrix(mo_set,work1,context,globenv)
      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(p%matrix)

     ELSE IF (density_guess == "CORE") THEN

      CALL eigensolver(work1,mo_set,ortho,work2,.false.,0.0_wp,context,globenv)
      IF (level_shift /= 0.0_wp) THEN
        CALL get_mo_set(mo_set=mo_set,&
                        eigenvectors=mo_eigenvectors)
        CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
      END IF
      CALL calculate_density_matrix(mo_set,work1,context,globenv)
      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(p%matrix)

    ELSE IF (density_guess == "ATOMIC") THEN

      CALL blacs_set_all(work1,0.0_wp,context,globenv)

      sum = 0.0_wp

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxl=maxl)

      ALLOCATE (econf(0:maxl),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"econf",(maxl + 1)*int_size)

      CALL copy_sparse_to_blacs_matrix(s(1)%matrix,work2,context,globenv)

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             natom=natom,&
                             atom_list=atom_list,&
                             all_potential=all_potential,&
                             gth_potential=gth_potential,&
                             orb_basis_set=orb_basis_set)

        IF (ASSOCIATED(all_potential)) THEN
          CALL get_potential(potential=all_potential,elec_conf=elec_conf)
        ELSE IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,elec_conf=elec_conf)
        ELSE
          CYCLE
        END IF

        maxl = SIZE(elec_conf) - 1
        econf(:) = 0
        econf(0:maxl) = elec_conf(0:maxl)

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l,&
                               first_sgf=first_sgf,&
                               last_sgf=last_sgf)

        DO iset=1,nset
          DO ishell=1,nshell(iset)
            la = l(ishell,iset)
            nelec = 4*la + 2
            IF (econf(la) > 0) THEN
              IF (econf(la) >= nelec) THEN
                paa = 2.0_wp
                econf(la) = econf(la) - nelec
              ELSE
                paa = REAL(2*econf(la),wp)/REAL(nelec,wp)
                econf(la) = 0
              END IF
              first_sgfa = first_sgf(ishell,iset)
              last_sgfa = last_sgf(ishell,iset)
              DO isgfa=first_sgfa,last_sgfa
                DO iatom=1,natom
                  atom_a = atom_list(iatom)
                  isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                  CALL blacs_set_element(work1,isgf,isgf,paa,context,globenv)
                  CALL blacs_get_element(work2,isgf,isgf,saa,context,globenv)
                  sum = sum + paa*saa
                END DO
              END DO
            END IF
          END DO
        END DO

      END DO

      DEALLOCATE (econf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"econf")

      CALL mp_sum(sum,globenv%group)

      saa = REAL(nelectron,wp)/sum

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             natom=natom,&
                             atom_list=atom_list,&
                             orb_basis_set=orb_basis_set)

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l,&
                               first_sgf=first_sgf,&
                               last_sgf=last_sgf)

        DO iset=1,nset
          DO ishell=1,nshell(iset)
            first_sgfa = first_sgf(ishell,iset)
            last_sgfa = last_sgf(ishell,iset)
            DO isgfa=first_sgfa,last_sgfa
              DO iatom=1,natom
                atom_a = atom_list(iatom)
                isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                CALL blacs_get_element(work1,isgf,isgf,paa,context,globenv)
                paa = paa*saa
                CALL blacs_set_element(work1,isgf,isgf,paa,context,globenv)
              END DO
            END DO
          END DO
        END DO
      END DO

      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv)

    ELSE

      CALL stop_program("SUBROUTINE calculate_first_density_matrix "//&
                        "(MODULE qs_scf)",&
                        "An invalid keyword for the initial density "//&
                        "guess was specified")

    END IF

  END SUBROUTINE calculate_first_density_matrix

! *****************************************************************************

  SUBROUTINE read_scf_parameters(globenv)

!   Purpose: Read the parameters for the SCF run.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    charge = 0
    density_guess = "ATOMIC"
    eps_eigval = 1.0E-5_wp
    eps_scf = 1.0E-5_wp
    eps_diis = 0.1_wp
    level_shift = 0.0_wp
    max_diis = 4
    max_scf = 50
    nprow = 0
    npcol = 0
    nrebuild = 1
    nrow_block = 32
    ncol_block = 32
    p_mix = 0.4_wp
    smear = 0.0_wp
    work_syevx = 1.0_wp
    use_cholesky = .true.
    use_arpack = .false.

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="SCF")

    DO WHILE (test_object(newline=.TRUE.) /= "EOS")
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("CHARGE")
        CALL read_object(charge)
      CASE ("DENSITY_GUESS","SCF_GUESS","GUESS")
        CALL read_object(density_guess,lower_to_upper=.TRUE.)
      CASE ("DENSITY_MIXING","MIXING")
        CALL read_object(p_mix)
      CASE ("EPS_DIIS")
        CALL read_object(eps_diis)
      CASE ("EPS_EIGVAL")
        CALL read_object(eps_eigval)
        use_cholesky = .false.
      CASE ("CHOLESKY_ON")
        use_cholesky = .true.
      CASE ("CHOLESKY_OFF")
        use_cholesky = .false.
      CASE ("ARPACK_ON") 
        use_cholesky = .false.
        use_arpack = .true.
      CASE ("EPS_SCF")
        CALL read_object(eps_scf)
      CASE ("LEVEL_SHIFT","LSHIFT")
        CALL read_object(level_shift)
      CASE ("MAX_DIIS")
        CALL read_object(max_diis)
      CASE ("MAX_SCF")
        CALL read_object(max_scf)
      CASE ("NREBUILD")
        CALL read_object(nrebuild)
        nrebuild = MAX(1,nrebuild)
      CASE ("BLOCKSIZE")
        CALL read_object(nrow_block)
        ncol_block = nrow_block
      CASE ("PROCESS_GRID")
        CALL read_object(nprow)
        CALL read_object(npcol)
      CASE ("SMEAR")
        CALL read_object(smear)
      CASE ("WORK_SYEVX")
        CALL read_object(work_syevx)
        work_syevx = MIN(MAX(0.0_wp,work_syevx),1.0_wp)
      END SELECT
    END DO

    CALL finish_parser()

  END SUBROUTINE read_scf_parameters

! *****************************************************************************

  SUBROUTINE write_scf_parameters(globenv)

!   Purpose: Write the parameters for the SCF run to the output unit.

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (max_scf < 1) RETURN

    IF (globenv%print%scf) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/)") "SCF PARAMETERS"
      WRITE (UNIT=output_unit,&
             FMT="(T3,A,/,4(/,T3,A,I5),3(/,T3,A,ES9.2),4(/,T3,A,F5.2))")&
        "Density guess: "//TRIM(density_guess),&
        "charge:            ",charge,&
        "nrebuild:          ",nrebuild,&
        "max_scf:           ",max_scf,&
        "max_diis:          ",max_diis,&
        "eps_scf:           ",eps_scf,&
        "eps_diis:          ",eps_diis,&
        "eps_eigval:        ",eps_eigval,&
        "p_mix:             ",p_mix,&
        "work_syevx:        ",work_syevx,&
        "level_shift [a.u.]:",level_shift,&
        "smear [a.u.]:      ",smear
    END IF

  END SUBROUTINE write_scf_parameters

! *****************************************************************************

!!****f* qs_scf/pao_opt_create [1.0] *
!!
!!   NAME
!!     pao_opt_create
!!
!!   FUNCTION
!!    allocates and initializes the optimizer and the objective function
!!    
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine pao_opt_create(pao_optimizer,pao_obj_f,qs_env,global_env,&
     error)
  type(cp_lbfgs_opt_gerd_f_type), pointer :: pao_optimizer
  type(pao_gerd_functional_type), pointer :: pao_obj_f
  type(qs_environment_type), intent(inout), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure,pao
  character(len=*), parameter :: routineN='pao_opt_create',&
        routineP=moduleN//':'//routineN
  type(dft_control_type), pointer :: dft_control
  type(pao_env_type), pointer :: pao_env
  real(kind=wp), dimension(:), pointer :: angles
  type(real_matrix_p_type), dimension(:), pointer :: p
  failure=.false.
  
  nullify(pao_optimizer,pao_obj_f)

  ! checks if pao is active
  call get_qs_env(qs_env,dft_control=dft_control, pao_env=pao_env, p=p)
  pao=dft_control%qs_control%pao
  if (pao) then
     CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
     if (.not.failure) then
        allocate(pao_optimizer, pao_obj_f)
        call pao_gerd_f_init(pao_obj_f,qs_env=qs_env, global_env=global_env,&
             scf_work1=scf_work1, scf_work2=scf_work2, scf_work3=ortho,&
             error=error)
        call pao_g_ang_get(pao_env%angles_att,qs_env=qs_env,&
             global_env=global_env, angles=angles,error=error)
        if (.not.associated(angles)) then
           call pao_g_ang_set(pao_env%angles_att,qs_env=qs_env,&
                global_env=global_env,angles_val=0.0_wp)
           call pao_g_ang_get(pao_env%angles_att,qs_env=qs_env,&
                global_env=global_env, angles=angles,error=error)
        end if
        call cp_opt_gerd_f_init(pao_optimizer, global_env=global_env,&
             obj_funct=pao_obj_f, x0=angles,wanted_relative_f_delta=0.0_wp,&
             wanted_projected_gradient=0.0_wp,max_f_per_iter=5,error=error)

        call pao_qs_env_did_change(qs_env,s_changed=.true.,h_changed=.true.,&
             full_reset=.false., global_env=global_env)
        
        ! map first_density to min_bas (should be calculated in the min_bas)
        call sm_matrix_transf_to_m(source_m=p(1)%matrix, &
             target_m=pao_env%min_density_m(1)%matrix%matrix,&
             glob_angles=pao_env%angles_att,&
             qs_env=qs_env, global_env=global_env, error=error)
     end if     
  end if
end subroutine pao_opt_create
!***************************************************************************

!!****f* qs_scf/pao_opt_dealloc [1.0] *
!!
!!   NAME
!!     pao_opt_dealloc
!!
!!   FUNCTION
!!     deallocates the types needed for the optimization in the pao process:
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine pao_opt_dealloc(pao_optimizer,pao_obj_f,qs_env,global_env,&
     error)
  type(cp_lbfgs_opt_gerd_f_type), pointer :: pao_optimizer
  type(pao_gerd_functional_type), pointer :: pao_obj_f
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_opt_dealloc',&
        routineP=moduleN//':'//routineN
  integer :: stat
  failure=.false.
  
  if (associated(pao_optimizer).and.associated(pao_obj_f)) then
     call cp_opt_gerd_f_dealloc_ref(pao_optimizer)
     call pao_gerd_f_dealloc_ref(pao_obj_f)
     deallocate(pao_optimizer, pao_obj_f, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  end if
end subroutine pao_opt_dealloc
!***************************************************************************

END MODULE qs_scf
