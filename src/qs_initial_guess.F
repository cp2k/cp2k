!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2006  CP2K developers group
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_initial_guess [1.0] *
!!
!!   NAME
!!     qs_initial_guess
!!
!!   FUNCTION
!!     Routines to somehow generate an intial guess
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!      2006.03 Moved here from qs_scf.F [Joost VandeVondele]
!!
!!   SOURCE
!******************************************************************************

MODULE qs_initial_guess
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: atomic_guess,&
                                             core_guess,&
                                             densities_guess,&
                                             random_guess,&
                                             restart_guess
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_collocate_density,            ONLY: collocate_atomic_charge_density
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             make_basis
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict,&
                                             read_mo_set
  USE qs_scf_methods,                  ONLY: eigensolver,&
                                             simple_eigensolver
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: get_matrix_diagonal,&
                                             real_matrix_p_type,&
                                             set_matrix_diagonal
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_initial_guess'

  PUBLIC ::  calculate_first_density_matrix

CONTAINS
  
!!****f* qs_initial_guess/calculate_first_density_matrix *
!!
!!   NAME
!!     calculate_first_density_matrix
!!
!!   FUNCTION
!!     can use a variety of methods to come up with an initial
!!     density matrix and optionally an initial wavefunction
!!
!!   NOTES
!!     badly needs to be split in subroutines each doing one of the possible
!!     schemes
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 moved here from qs_scf [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE calculate_first_density_matrix(scf_env,qs_env,globenv,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calculate_first_density_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, group, handle, homo, iatom, ikind, iset, isgf, isgfa, &
      ishell, ispin, istat, la, maxl, maxll, nao, natom, ncount, nelectron, &
      nmo, nset, nspin, qs_env_id
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, l, last_sgfa
    LOGICAL                                  :: failure, id_equal, ortho_basis
    REAL(KIND=dp)                            :: maxocc, nelec, paa, scale, &
                                                total_rho, trps1, trps2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: econf, pdiag, sdiag
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work1
    TYPE(cp_fm_type), POINTER                :: mo_coeff, ortho, sv, work2
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: h_core_sparse, p_rmpv, &
                                                s_sparse
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input

!   ---------------------------------------------------------------------------

    NULLIFY(all_potential,atomic_kind, mo_coeff, sv,&
         gth_potential, orb_basis_set, atomic_kind_set, particle_set,&
         ortho,work2,work1,mo_array,s_sparse,p_rmpv,scf_control, &
         dft_control,h_core_sparse)
    NULLIFY(dft_section, input)

    CALL timeset("calculate_first_density_matrix","I"," ",handle)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,mos=mo_array, matrix_s=s_sparse,&
         matrix_h=h_core_sparse,input=input,&
         scf_control=scf_control, id_nr=qs_env_id, dft_control=dft_control,&
         error=error)

    nspin=dft_control%nspins

    p_rmpv => qs_env%rho%rho_ao
    work1 => scf_env%scf_work1
    work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)
    ortho_basis = .FALSE.
    IF ( dft_control%qs_control%semi_empirical ) THEN
       IF ( dft_control%qs_control%se_control%orthogonal_basis ) &
          ortho_basis = .TRUE.
    END IF

    IF (scf_control%use_ot.AND.&
        (.NOT.((scf_control%density_guess == random_guess).OR.&
               (scf_control%density_guess == atomic_guess).OR.&
               ((scf_control%density_guess == restart_guess).AND.&
                (scf_control%level_shift == 0.0_dp))))) THEN
       CALL stop_program("calculate_first_density_matrix",&
            "OT needs GUESS ATOMIC / RESTART : other options NYI")
    END IF

    IF (scf_control%density_guess == restart_guess) THEN

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,globenv,&
            id_nr=1,dft_section=dft_section,error=error)

       DO ispin=1,nspin
          IF (scf_control%level_shift /= 0.0_dp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  mo_coeff=mo_coeff)
             CALL cp_fm_to_fm(mo_coeff,ortho)
          END IF

          ! make all nmo vectors present orthonormal
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo, homo=homo)

          IF(ortho_basis) THEN
            CALL make_basis(mo_coeff,nmo)
          ELSE
            ! ortho so that one can restart for different positions (basis sets?)
            CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
            CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
            CALL make_basis(mo_coeff,homo, &
                 ortho=sv,otype="SV")
            CALL cp_fm_release(sv)
          ENDIF
          ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == random_guess) THEN

       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo)
          CALL cp_fm_init_random(mo_coeff,nmo)
          IF(ortho_basis) THEN
            CALL make_basis(mo_coeff,nmo)
          ELSE
            ! ortho so that one can restart for different positions (basis sets?)
            CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
            CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
            CALL make_basis(mo_coeff,nmo, &
                 ortho=sv,otype="SV")
            CALL cp_fm_release(sv)
          ENDIF
          ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == core_guess) THEN
 
       IF (dft_control%restricted) CALL stop_program("calculate_first_density_matrix","Option not yet implemented")

       ! work1(1) should contain the core hamiltonian
       ! we need a copy the core matrix for every spin (the dumbest solution,
       ! we should copy the mos)..
       
       CALL copy_sm_to_fm(h_core_sparse(1)%matrix,work1(1)%matrix)

       IF (dft_control%nspins.eq.2) THEN
          IF (.NOT.scf_control%use_ot) THEN
             CALL cp_fm_to_fm(work1(1)%matrix,&
                  work1(2)%matrix)
          ENDIF
       ENDIF
       DO ispin=1,nspin
          IF(ortho_basis) THEN
            CALL simple_eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,&
                 work2,.FALSE.,0.0_dp,work_syevx=scf_control%work_syevx,&
                 smear=0.0_dp,error=error)
          ELSE
            CALL eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,ortho,work2, &
               .FALSE.,0.0_dp,&
               use_cholesky=scf_control%use_cholesky, &
               work_syevx=scf_control%work_syevx,&
               use_jacobi=.FALSE.,&
               jacobi_threshold=scf_control%jacobi_threshold,&
               smear=0.0_dp,&
               error=error)
          END IF
          IF (scf_control%level_shift /= 0.0_dp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  mo_coeff=mo_coeff)
             CALL cp_fm_to_fm(mo_coeff,ortho)
          END IF
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == atomic_guess) THEN

      IF (dft_control%restricted) CALL stop_program("calculate_first_density_matrix","Option not yet implemented")

      group = qs_env%para_env%group

      natom = SIZE(particle_set)
      ALLOCATE (first_sgf(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"first_sgf",natom*int_size)
      CALL get_particle_set(particle_set=particle_set,&
                            first_sgf=first_sgf)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxl)
      ALLOCATE (econf(0:maxl),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"econf",(maxl + 1)*int_size)

      CALL get_mo_set(mo_array(1)%mo_set,nao=nao)
      ALLOCATE (pdiag(nao),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"pdiag",nao*dp_size)
      pdiag(:) = 0.0_dp

      ALLOCATE (sdiag(nao),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"sdiag",nao*dp_size)
      IF (ortho_basis) THEN
        sdiag=1._dp
      ELSE
        CALL get_matrix_diagonal(s_sparse(1)%matrix,sdiag)
        CALL mp_sum(sdiag,group)
      END IF

      DO ispin=1,nspin

        CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                        maxocc=maxocc,&
                        nelectron=nelectron)

        ncount = 0
        trps1 = 0.0_dp
        trps2 = 0.0_dp
        pdiag(:) = 0.0_dp

        IF (nelectron /= 0) THEN 

           DO ikind=1,SIZE(atomic_kind_set)

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

              maxll = MIN(SIZE(elec_conf) - 1,maxl)
              econf(:) = 0.0_dp
              econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)

              CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                     nset=nset,&
                                     nshell=nshell,&
                                     l=l,&
                                     first_sgf=first_sgfa,&
                                     last_sgf=last_sgfa)

              DO iset=1,nset
                 DO ishell=1,nshell(iset)
                    la = l(ishell,iset)
                    nelec = maxocc*REAL(2*la + 1,dp)
                    IF (econf(la) > 0.0_dp) THEN
                       IF (econf(la) >= nelec) THEN
                          paa = maxocc
                          econf(la) = econf(la) - nelec
                       ELSE
                          paa = maxocc*econf(la)/nelec
                          econf(la) = 0.0_dp
                          ncount = ncount + NINT(nelec/maxocc)
                       END IF
                       DO isgfa=first_sgfa(ishell,iset),last_sgfa(ishell,iset)
                          DO iatom=1,natom
                             atom_a = atom_list(iatom)
                             isgf = first_sgf(atom_a) + isgfa - 1
                             pdiag(isgf) = paa
                             IF (paa == maxocc) THEN
                                trps1 = trps1 + paa*sdiag(isgf)
                             ELSE
                                trps2 = trps2 + paa*sdiag(isgf)
                             END IF
                          END DO
                       END DO
                    END IF
                 END DO
              END DO

           END DO

           IF (trps2 == 0.0_dp) THEN
              DO isgf=1,nao
                 IF (sdiag(isgf) > 0.0_dp) pdiag(isgf) = pdiag(isgf)/sdiag(isgf)
              END DO
           ELSE
              scale = (REAL(nelectron,dp) - trps1)/trps2
              DO isgf=1,nao 
                 IF (pdiag(isgf) < maxocc) pdiag(isgf) = scale*pdiag(isgf)
              END DO
           END IF
        END IF

        CALL set_matrix_diagonal(p_rmpv(ispin)%matrix,pdiag)

        ! do we require an initial wavefunction ?
        ! assume p is a projector (it is not ..) and generate an initial  
        ! wavefunction and the corresponding dm
        IF (scf_control%use_ot) THEN
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo)
          CALL cp_fm_init_random(mo_coeff,nmo)
          CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")

          ! multiply times PS
          IF (ortho_basis) THEN
             CALL cp_fm_to_fm(mo_coeff,sv)
          ELSE
             CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
          ENDIF
          ! here we could easily multiply with the diag that we actually have replicated already 
          CALL cp_sm_fm_multiply(p_rmpv(ispin)%matrix,sv,mo_coeff,nmo)

          ! and ortho the result
          IF (ortho_basis) THEN
             CALL make_basis(mo_coeff,nmo)
          ELSE
             CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
             CALL make_basis(mo_coeff,nmo, ortho=sv,otype="SV")
          ENDIF
          CALL cp_fm_release(sv)

          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
        ENDIF

      END DO

      DEALLOCATE (econf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"econf")

      DEALLOCATE (first_sgf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"first_sgf")

      DEALLOCATE (pdiag,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"pdiag")

      DEALLOCATE (sdiag,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"sdiag")

    ELSE IF (scf_control%density_guess == densities_guess) THEN

      ! Collocation of the density into the PW-grid
      CALL collocate_atomic_charge_density(total_rho=total_rho, qs_env=qs_env, error=error)

      ! do some assertions here on these matrices having the same structure, 
      ! as is currently required
      DO ispin=1,SIZE(qs_env%matrix_s)
        id_equal=(qs_env%matrix_s(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      ! matrix_p = rho%rho_ao
      DO ispin=1,SIZE(qs_env%rho%rho_ao)
        id_equal=(qs_env%rho%rho_ao(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      DO ispin=1,SIZE(qs_env%matrix_h)
        id_equal=(qs_env%matrix_h(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      DO ispin=1,SIZE(qs_env%matrix_ks)
        id_equal=(qs_env%matrix_ks(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      ! end sparsity check

      qs_env%scf_env%iter_count = 1

      ! qs_ks_build_kohn_sham_matrix (Integration)
      !CALL qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env, qs_env=qs_env, &
      !                                  ks_matrix=qs_env%matrix_ks, rho=qs_env%rho, &
      !                                  energy=qs_env%energy, calculate_forces=.TRUE., &
      !                                  just_energy=.FALSE., kg_gpw=.FALSE., error=error)
      CALL qs_ks_did_change(ks_env=qs_env%ks_env, rho_changed=.TRUE., error=error)

      CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, qs_env=qs_env, &
                               calculate_forces=.TRUE., just_energy=.FALSE., error=error)

      ! diagonalization
      DO ispin = 1,qs_env%dft_control%nspins
        CALL copy_sm_to_fm(qs_env%matrix_ks(ispin)%matrix, qs_env%scf_env%scf_work1(ispin)%matrix)
      END DO

      qs_env%scf_env%iter_method = "Mixing/Diag"
      qs_env%scf_env%iter_delta = 0.0_dp

      DO ispin = 1,qs_env%dft_control%nspins
        CALL eigensolver(matrix_ks=qs_env%scf_env%scf_work1(ispin)%matrix, &
                         mo_set=qs_env%mos(ispin)%mo_set, ortho=qs_env%scf_env%ortho, &
                         work=qs_env%scf_env%scf_work2, do_level_shift=.FALSE., &
                         level_shift=scf_control%level_shift, &
                         use_cholesky=scf_control%use_cholesky, work_syevx=scf_control%work_syevx, &
                         use_jacobi=.FALSE., jacobi_threshold=scf_control%jacobi_threshold, &
                         smear=scf_control%smear, error=error)

        ! calculate_density_matrix
        CALL calculate_density_matrix(qs_env%mos(ispin)%mo_set, qs_env%rho%rho_ao(ispin)%matrix, &
                                      error=error)
      END DO

    ELSE

      CALL stop_program(routineP,&
                        "An invalid keyword for the initial density "//&
                        "guess was specified")

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_first_density_matrix

END MODULE qs_initial_guess
