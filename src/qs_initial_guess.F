!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to somehow generate an intial guess
!> \par History
!>       2006.03 Moved here from qs_scf.F [Joost VandeVondele]
! *****************************************************************************
MODULE qs_initial_guess

  USE atom_kind_orbitals,              ONLY: calculate_atomic_orbitals
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_check
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE external_potential_types,        ONLY: all_potential_type,&
                                             gth_potential_type
  USE f77_blas
  USE input_constants,                 ONLY: atomic_guess,&
                                             core_guess,&
                                             densities_guess,&
                                             history_guess,&
                                             mopac_guess,&
                                             no_guess,&
                                             random_guess,&
                                             restart_guess
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_collocate_density,            ONLY: collocate_atomic_charge_density
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             make_basis_simple,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict,&
                                             read_mo_set,&
                                             set_mo_occupation,&
                                             wfn_restart_file_name
  USE qs_scf_methods,                  ONLY: eigensolver,&
                                             simple_eigensolver
  USE qs_scf_types,                    ONLY: block_krylov_diag_method_nr,&
                                             ot_diag_method_nr,&
                                             qs_scf_env_type
  USE qs_wf_history_methods,           ONLY: wfi_update
  USE scf_control_types,               ONLY: scf_control_type
  USE scp_environment_types,           ONLY: scp_environment_type
  USE sparse_matrix_types,             ONLY: &
       cp_sm_scale_and_add, cp_sm_sm_trace, first_block_node, get_block_node, &
       get_matrix_diagonal, next_block_node, put_block_node, &
       real_block_node_type, real_matrix_p_type, set_matrix_diagonal
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_initial_guess'

  PUBLIC ::  calculate_first_density_matrix

CONTAINS

! *****************************************************************************
!> \brief can use a variety of methods to come up with an initial
!>      density matrix and optionally an initial wavefunction
!> \note
!>      badly needs to be split in subroutines each doing one of the possible
!>      schemes
!> \par History
!>      03.2006 moved here from qs_scf [Joost VandeVondele]
!>      06.2007 allow to skip the initial guess [jgh]
! *****************************************************************************
  SUBROUTINE calculate_first_density_matrix(scf_env,qs_env,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calculate_first_density_matrix', &
      routineP = moduleN//':'//routineN

    TYPE atom_matrix_type
      REAL(KIND=dp), DIMENSION(:,:), POINTER   :: mat
    END TYPE atom_matrix_type

    CHARACTER(LEN=default_path_length)       :: file_name, filename
    INTEGER :: atom_a, density_guess, group, handle, homo, i, iatom, id_nr, &
      ikind, iset, isgf, isgfa, ishell, ispin, istat, j, la, last_read, maxl, &
      maxll, nao, natom, ncount, nelectron, nmo, not_read, nset, nsgf, nspin, &
      nvec, qs_env_id, z, irow, icol
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, kind_of
    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, l, last_sgfa
    LOGICAL                                  :: exist, failure, id_equal, &
                                                ortho_basis, scp
    LOGICAL                                  :: ionode
    INTEGER                                  :: output_unit
    REAL(KIND=dp)                            :: maxocc, nelec, paa, scale, &
                                                sum, total_rho, trps1, trps2, &
                                                yy, zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: econf, pdiag, sdiag
    REAL(KIND=dp), DIMENSION(0:3)            :: edftb
    TYPE(atom_matrix_type), DIMENSION(:), &
      POINTER                                :: pmat
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work1
    TYPE(cp_fm_type), POINTER                :: mo_coeff, moa, mob, ortho, &
                                                sv, work2
    TYPE(cp_logger_type), POINTER            :: logger
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
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: dft_section, input

    logger => cp_error_get_logger(error)
    failure = .FALSE.
    NULLIFY(all_potential,atomic_kind, mo_coeff, sv,&
         gth_potential, orb_basis_set, atomic_kind_set, particle_set,&
         ortho,work2,work1,mo_array,s_sparse,p_rmpv,scf_control, &
         dft_control,h_core_sparse)
    NULLIFY(dft_section, input)
    NULLIFY(moa,mob)

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,mos=mo_array, matrix_s=s_sparse,&
         matrix_h=h_core_sparse,input=input,scp_env=scp_env,&
         scf_control=scf_control, id_nr=qs_env_id, dft_control=dft_control,&
         error=error)

    nspin=dft_control%nspins
    scp = dft_control%scp
    p_rmpv => qs_env%rho%rho_ao
    work1 => scf_env%scf_work1
    work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)

    density_guess=scf_control%density_guess

    ortho_basis = .FALSE.
    IF ( dft_control%qs_control%semi_empirical ) THEN
       IF ( dft_control%qs_control%se_control%orthogonal_basis ) &
          ortho_basis = .TRUE.
       ! if non-explicit basis, use Mopac Atomic guess
       IF (density_guess == atomic_guess) density_guess = mopac_guess
    END IF
    IF ( dft_control%qs_control%dftb ) THEN
       IF ( dft_control%qs_control%dftb_control%orthogonal_basis ) &
          ortho_basis = .TRUE.
       ! if non-explicit basis, use Mopac Atomic guess
       IF (density_guess == atomic_guess) density_guess = mopac_guess
    END IF

    IF (scf_control%use_ot.AND.&
        (.NOT.((density_guess == random_guess).OR.&
               (density_guess == atomic_guess).OR.&
               (density_guess == mopac_guess).OR.&
               (((density_guess == restart_guess).OR.&
                (density_guess == history_guess)).AND.&
                (scf_control%level_shift == 0.0_dp))))) THEN
       CALL stop_program("calculate_first_density_matrix",&
            "OT needs GUESS ATOMIC / RESTART / HISTORY RESTART: other options NYI")
    END IF

    ! if a restart was requested, check that the file exists,
    ! if not we fall back to an atomic guess. No kidding, the file name should remain
    ! in sync with read_mo_set_from_restart
    id_nr=0
    IF (density_guess == restart_guess) THEN
        CALL wfn_restart_file_name(file_name,exist,dft_section,logger,error=error)
        IF (.NOT.exist) THEN
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                "User requested to restart the wavefunction from the file named: "//&
                TRIM(file_name)//". This file does not exist. Please check the existence of"//&
                " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME."//&
                " Calculation continues using ATOMIC GUESS. "//&
CPSourceFileRef,&
                only_ionode=.TRUE.)
           density_guess = atomic_guess
        END IF
    ELSE IF (density_guess == history_guess) THEN
       CALL wfn_restart_file_name(file_name,exist,dft_section,logger,error=error)
       nvec = qs_env%wf_history%memory_depth
       not_read = nvec+1
       ! At this level we read the saved backup RESTART files..
       DO i=1,nvec
          j = i - 1
          filename = TRIM(file_name)
          IF (j/=0) filename = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(j))
          INQUIRE(FILE=filename,exist=exist)
          IF ((.NOT. exist) .AND. (i < not_read)) THEN
             not_read = i
          END IF
       END DO
       IF (not_read == 1) THEN
          density_guess = restart_guess
          filename = TRIM(file_name)
          INQUIRE(FILE=filename,exist=exist)
          IF (.NOT. exist) THEN
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "User requested to restart the wavefunction from a series of restart files named: "//&
                  TRIM(file_name)//" with extensions (.bak-n). These files do not exist."//&
                  " Even trying to switch to a plain restart wave-function failes because the"//&
                  " file named: "//TRIM(file_name)//" does not exist. Please check the existence of"//&
                  " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME. "//&
                  " Calculation continues using ATOMIC GUESS. "//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             density_guess = atomic_guess
          END IF
       END IF
       last_read = not_read - 1
    END IF

    IF (density_guess == restart_guess) THEN

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr=id_nr,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
            scp=scp, scp_env=scp_env, error=error)

       DO ispin=1,nspin
          IF (scf_control%level_shift /= 0.0_dp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                             mo_coeff=mo_coeff)
             CALL cp_fm_to_fm(mo_coeff,ortho,error=error)
          END IF

          ! make all nmo vectors present orthonormal
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo, homo=homo)

          IF(ortho_basis) THEN
            CALL make_basis_simple(mo_coeff,nmo,error=error)
          ELSE
            ! ortho so that one can restart for different positions (basis sets?)
            CALL make_basis_sm(mo_coeff,homo,s_sparse(1)%matrix,error=error)
          ENDIF
         ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array,error=error)

       DO ispin=1,nspin
          CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                 smear=qs_env%scf_control%smear,&
                                 error=error)
       ENDDO

       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (density_guess == history_guess) THEN
       IF (not_read > 1) THEN
          DO i=1, last_read
             j = last_read -i
             CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
             CALL read_mo_set(mo_array,atomic_kind_set,particle_set,qs_env%para_env,&
                  id_nr=j,multiplicity=dft_control%multiplicity,scp=scp, scp_env=scp_env,&
                  dft_section=dft_section, error=error)
          
             DO ispin=1,nspin
                IF (scf_control%level_shift /= 0.0_dp) THEN
                   CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                   mo_coeff=mo_coeff)
                   CALL cp_fm_to_fm(mo_coeff,ortho,error=error)
                END IF
            
                ! make all nmo vectors present orthonormal
                CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, nmo=nmo, homo=homo)
          
                IF(ortho_basis) THEN
                   CALL make_basis_simple(mo_coeff,nmo,error=error)
                ELSE
                   ! ortho so that one can restart for different positions (basis sets?)
                   CALL make_basis_sm(mo_coeff,homo,s_sparse(1)%matrix,error=error)
                ENDIF
                ! only alpha spin is kept for restricted
                IF (dft_control%restricted) EXIT
             END DO
             IF (dft_control%restricted) CALL mo_set_restrict(mo_array,error=error)

             DO ispin=1,nspin
                CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                       smear=qs_env%scf_control%smear,&
                                       error=error)
             ENDDO

             DO ispin=1,nspin
                CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                              p_rmpv(ispin)%matrix,error=error)
             ENDDO
          
             ! Write to extrapolation pipeline
             CALL wfi_update(wf_history=qs_env%wf_history, qs_env=qs_env, dt=1.0_dp, error=error)
          END DO
       END IF

    ELSE IF (density_guess == random_guess) THEN

       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo)
          CALL cp_fm_init_random(mo_coeff,nmo,error=error)
          IF(ortho_basis) THEN
            CALL make_basis_simple(mo_coeff,nmo,error=error)
          ELSE
            CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
          ENDIF
          ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array,error=error)

       DO ispin=1,nspin
          CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                 smear=qs_env%scf_control%smear,&
                                 error=error)
       ENDDO

       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (density_guess == core_guess) THEN

       ispin = 1

       ! Load core Hamiltonian into work matrix

       CALL copy_sm_to_fm(h_core_sparse(1)%matrix,work1(ispin)%matrix,error=error)

       ! Diagonalize the core Hamiltonian matrix and retrieve a first set of
       ! molecular orbitals (MOs)

       IF (ortho_basis) THEN
          CALL simple_eigensolver(matrix_ks=work1(ispin)%matrix,&
                                  mo_set=mo_array(ispin)%mo_set,&
                                  work=work2,&
                                  do_level_shift=.FALSE.,&
                                  level_shift=0.0_dp,&
                                  use_jacobi=.FALSE.,jacobi_threshold=0._dp,&
                                  error=error)
       ELSE
          CALL eigensolver(matrix_ks=work1(ispin)%matrix,&
                           mo_set=mo_array(ispin)%mo_set,&
                           ortho=ortho,&
                           work=work2,&
                           do_level_shift=.FALSE.,&
                           level_shift=0.0_dp,&
                           cholesky_method=scf_env%cholesky_method,&
                           use_jacobi=.FALSE.,&
                           jacobi_threshold=scf_control%diagonalization%jacobi_threshold,&
                           error=error)
       END IF

       ! Open shell case: copy alpha MOs to beta MOs

       IF (nspin == 2) THEN
          CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=moa)
          CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mob,nmo=nmo)
          CALL cp_fm_to_fm(moa,mob,nmo)
       END IF

       ! Build an initial density matrix (for each spin in the case of
       ! an open shell calculation) from the first MOs set

       DO ispin=1,nspin
         CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)
         CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                       p_rmpv(ispin)%matrix,&
                                       error=error)
       END DO

    ELSE IF (density_guess == atomic_guess) THEN

       group = qs_env%para_env%group
       natom = SIZE(particle_set)
       ALLOCATE (kind_of(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,kind_of=kind_of)
       
       ALLOCATE (pmat(SIZE(atomic_kind_set)),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(ikind)
          NULLIFY(pmat(ikind)%mat)
          CALL calculate_atomic_orbitals(atomic_kind,pmat=pmat(ikind)%mat,error=error)
          NULLIFY(atomic_kind)
       END DO

       scale = 1._dp
       IF (nspin==2) scale=0.5_dp
       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                          maxocc=maxocc,&
                          nelectron=nelectron)
          DO irow=1,p_rmpv(ispin)%matrix%nblock_row
            block_node => first_block_node(p_rmpv(ispin)%matrix,irow)
            DO WHILE (ASSOCIATED(block_node))
              CALL get_block_node(block_node=block_node,block_col=icol)
              IF (irow == icol) THEN
                ikind = kind_of(irow)
                pmat(ikind)%mat=pmat(ikind)%mat*scale
                CALL put_block_node(block_node,p_rmpv(ispin)%matrix,irow,icol,pmat(ikind)%mat)
                pmat(ikind)%mat=pmat(ikind)%mat/scale
              END IF
              block_node => next_block_node(block_node)
            END DO
          END DO
          CALL cp_sm_sm_trace(p_rmpv(ispin)%matrix,s_sparse(1)%matrix,trps1,qs_env%para_env,error)
          scale=REAL(nelectron,dp)/trps1
          CALL cp_sm_scale_and_add(p_rmpv(ispin)%matrix,alpha=scale,error=error)

          ! The orbital transformation method (OT) requires not only an
          ! initial density matrix, but also an initial wavefunction (MO set)
          IF (scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr .OR. &
               scf_env%method==block_krylov_diag_method_nr ) THEN
             IF (dft_control%restricted.AND.(ispin == 2)) THEN
                CALL mo_set_restrict(mo_array,error=error)
             ELSE
                CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                mo_coeff=mo_coeff,&
                                nmo=nmo, nao=nao, homo=homo)
                mo_coeff%local_data = 0.0_dp
                CALL cp_fm_init_random(mo_coeff,nmo,error=error)
                CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV",error=error)
                ! multiply times PS
                IF (ortho_basis) THEN
                   CALL cp_fm_to_fm(mo_coeff,sv,error=error)
                ELSE
                   ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
                   CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo,error=error)
                END IF
                CALL cp_sm_fm_multiply(p_rmpv(ispin)%matrix,sv,mo_coeff,homo,error=error)
                CALL cp_fm_release(sv,error=error)
                ! and ortho the result
                IF (ortho_basis) THEN
                   CALL make_basis_simple(mo_coeff,nmo,error=error)
                ELSE
                   CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
                END IF
             END IF

             CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                    smear=qs_env%scf_control%smear,&
                                    error=error)

             CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                           p_rmpv(ispin)%matrix,&
                                           error=error)
          END IF
       END DO

       DO ikind=1,SIZE(atomic_kind_set)
          IF(ASSOCIATED(pmat(ikind)%mat)) THEN
            DEALLOCATE (pmat(ikind)%mat,STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
       END DO
       DEALLOCATE (pmat,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       DEALLOCATE (kind_of,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ELSE IF (density_guess == mopac_guess) THEN

       group = qs_env%para_env%group
       natom = SIZE(particle_set)
       ALLOCATE (first_sgf(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,error=error)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxl)
       ALLOCATE (econf(0:maxl),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       CALL get_mo_set(mo_array(1)%mo_set,nao=nao)
       ALLOCATE (pdiag(nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       pdiag(:) = 0.0_dp
       
       ALLOCATE (sdiag(nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF (ortho_basis) THEN
          sdiag(:) = 1.0_dp
       ELSE
          CALL get_matrix_diagonal(s_sparse(1)%matrix,sdiag)
          CALL mp_sum(sdiag,group)
       END IF

       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                          maxocc=maxocc,&
                          nelectron=nelectron)
          
          ncount   = 0
          trps1    = 0.0_dp
          trps2    = 0.0_dp
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

                IF ( dft_control%qs_control%dftb ) THEN
                   CALL get_dftb_atom_param(dftb_parameter=atomic_kind%dftb_parameter,&
                        lmax=maxll,occupation=edftb)
                   maxll = MIN(maxll,maxl)
                   econf(0:maxl)=edftb(0:maxl)
                ELSEIF (ASSOCIATED(all_potential)) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,elec_conf=elec_conf, nsgf=nsgf, z=z, zeff=zeff)
                   maxll = MIN(SIZE(elec_conf) - 1,maxl)
                   econf(:) = 0.0_dp
                   econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)
                ELSE IF (ASSOCIATED(gth_potential)) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,elec_conf=elec_conf, nsgf=nsgf, z=z, zeff=zeff)
                   maxll = MIN(SIZE(elec_conf) - 1,maxl)
                   econf(:) = 0.0_dp
                   econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)
                ELSE
                   CYCLE
                END IF
          
                ! MOPAC TYEP GUESS
                IF (dft_control%qs_control%dftb) THEN
                   DO iatom=1,natom
                      atom_a = atom_list(iatom)
                      isgfa=first_sgf(atom_a)
                      DO la=0,maxll
                         SELECT CASE (la)
                         CASE (0)
                            pdiag(isgfa) = econf(0)
                         CASE (1)
                            pdiag(isgfa+1) = econf(1)/3._dp
                            pdiag(isgfa+2) = econf(1)/3._dp
                            pdiag(isgfa+3) = econf(1)/3._dp
                         CASE (2)
                            pdiag(isgfa+4) = econf(2)/5._dp
                            pdiag(isgfa+5) = econf(2)/5._dp
                            pdiag(isgfa+6) = econf(2)/5._dp
                            pdiag(isgfa+7) = econf(2)/5._dp
                            pdiag(isgfa+8) = econf(2)/5._dp
                         CASE (3)
                            pdiag(isgfa+ 9) = econf(3)/7._dp
                            pdiag(isgfa+10) = econf(3)/7._dp
                            pdiag(isgfa+11) = econf(3)/7._dp
                            pdiag(isgfa+12) = econf(3)/7._dp
                            pdiag(isgfa+13) = econf(3)/7._dp
                            pdiag(isgfa+14) = econf(3)/7._dp
                            pdiag(isgfa+15) = econf(3)/7._dp
                         CASE DEFAULT
                            CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END SELECT
                      END DO
                   END DO
                ELSEIF (dft_control%qs_control%semi_empirical) THEN
                  yy = REAL(dft_control%charge,KIND=dp)/REAL(nao,KIND=dp)
                  DO iatom=1,natom
                    atom_a = atom_list(iatom)
                    isgfa=first_sgf(atom_a)
                    SELECT CASE (nsgf)
                    CASE (1) ! s-basis
                      pdiag(isgfa   ) = (zeff         - yy)*0.5_dp*maxocc
                    CASE (4) ! sp-basis
                      IF (z == 1 ) THEN
                        ! special case: hydrogen with sp basis
                        pdiag(isgfa   ) = (zeff         - yy)*0.5_dp*maxocc
                        pdiag(isgfa+ 1) = 0._dp
                        pdiag(isgfa+ 2) = 0._dp
                        pdiag(isgfa+ 3) = 0._dp
                      ELSE
                        pdiag(isgfa   ) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa+ 1) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa+ 2) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa+ 3) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                      END IF
                    CASE (9) ! spd-basis
                      IF (z < 21 .OR. z > 30 .AND. z < 39 .OR. z > 48 .AND. z < 57) THEN
                         !   Main Group Element:  The "d" shell is formally empty.
                         pdiag(isgfa   ) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 1) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 2) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 3) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 4) = (             - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 5) = (             - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 6) = (             - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 7) = (             - yy)*0.5_dp*maxocc
                         pdiag(isgfa+ 8) = (             - yy)*0.5_dp*maxocc
                      ELSE IF (z < 99) THEN
                         sum = zeff - 9.0_dp*yy
                         !   First, put 2 electrons in the 's' shell
                         pdiag(isgfa   ) = (MAX(0.0_dp, MIN(sum, 2.0_dp)))*0.5_dp*maxocc
                         sum = sum - 2.0_dp
                         IF (sum > 0.0_dp) THEN
                            !   Now put as many electrons as possible into the 'd' shell
                            pdiag(isgfa+ 4) = (MAX(0.0_dp, MIN(sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                            pdiag(isgfa+ 5) = (MAX(0.0_dp, MIN(sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                            pdiag(isgfa+ 6) = (MAX(0.0_dp, MIN(sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                            pdiag(isgfa+ 7) = (MAX(0.0_dp, MIN(sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                            pdiag(isgfa+ 8) = (MAX(0.0_dp, MIN(sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                            sum = MAX(0.0_dp, sum-10.0_dp)
                            !   Put the remaining electrons in the 'p' shell
                            pdiag(isgfa+ 1) = (sum/3.0_dp)*0.5_dp*maxocc
                            pdiag(isgfa+ 2) = (sum/3.0_dp)*0.5_dp*maxocc
                            pdiag(isgfa+ 3) = (sum/3.0_dp)*0.5_dp*maxocc
                         END IF
                      END IF
                    CASE DEFAULT
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                    END SELECT
                  END DO
                ELSE
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
                     END DO ! ishell
                  END DO ! iset
                END IF
             END DO ! ikind
          
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

          ! if we have many added MOS, the atomic guess as we have it now seriously breaks down (more than usual that is)
          ! i.e. a number of states are being projected on the null space of the density matrix, which causes failure on reortho
          ! as temporary fix (to be removed when we have a real atomic guess), we just shift the diagonal a tiny bit.
          ! a better initial guess for a finite temperature would take the atomic orbital energies into account, and would do some finite
          ! temperature approach. This hack can be removed as soon as we have real atomic wavefunctions and energies here.
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,nmo=nmo)
          IF (COUNT(pdiag.NE.0.0_dp)<nmo .AND. scf_control%use_ot) THEN
             pdiag=pdiag+0.01_dp
          ENDIF

          CALL set_matrix_diagonal(p_rmpv(ispin)%matrix,pdiag)
          
          ! The orbital transformation method (OT) requires not only an
          ! initial density matrix, but also an initial wavefunction (MO set)
          IF (scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr .OR. &
               scf_env%method==block_krylov_diag_method_nr ) THEN
             IF (dft_control%restricted.AND.(ispin == 2)) THEN
                CALL mo_set_restrict(mo_array,error=error)
             ELSE
                CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                mo_coeff=mo_coeff,&
                                nmo=nmo, homo=homo)
                CALL cp_fm_init_random(mo_coeff,nmo,error=error)
                CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV",error=error)
                ! multiply times PS
                IF (ortho_basis) THEN
                   CALL cp_fm_to_fm(mo_coeff,sv,error=error)
                ELSE
                   CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo,error=error)
                END IF
                ! here we could easily multiply with the diag that we actually have replicated already
                CALL cp_sm_fm_multiply(p_rmpv(ispin)%matrix,sv,mo_coeff,homo,error=error)
                CALL cp_fm_release(sv,error=error)
                ! and ortho the result
                IF (ortho_basis) THEN
                   CALL make_basis_simple(mo_coeff,nmo,error=error)
                ELSE
                   CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
                END IF
             END IF

             CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                    smear=qs_env%scf_control%smear,&
                                    error=error)

             CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                           p_rmpv(ispin)%matrix,&
                                           error=error)
          END IF
       END DO

       DEALLOCATE (econf,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       DEALLOCATE (first_sgf,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       DEALLOCATE (pdiag,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       DEALLOCATE (sdiag,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
    ELSE IF (density_guess == densities_guess) THEN
       
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
       
       CALL qs_ks_did_change(ks_env=qs_env%ks_env, rho_changed=.TRUE., error=error)

       CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, qs_env=qs_env, calculate_forces=.TRUE.,&
            just_energy=.FALSE., error=error)

       ! diagonalization
       DO ispin = 1,qs_env%dft_control%nspins
          CALL copy_sm_to_fm(qs_env%matrix_ks(ispin)%matrix,&
                             qs_env%scf_env%scf_work1(ispin)%matrix,&
                             error=error)
       END DO

       qs_env%scf_env%iter_method = "Mixing/Diag"
       qs_env%scf_env%iter_delta = 0.0_dp

       DO ispin = 1,qs_env%dft_control%nspins
         CALL eigensolver(matrix_ks=qs_env%scf_env%scf_work1(ispin)%matrix, &
                          mo_set=mo_array(ispin)%mo_set,&
                          ortho=qs_env%scf_env%ortho, &
                          work=qs_env%scf_env%scf_work2,&
                          do_level_shift=.FALSE., &
                          level_shift=scf_control%level_shift, &
                          cholesky_method=scf_env%cholesky_method,&
                          use_jacobi=.FALSE.,&
                          jacobi_threshold=scf_control%diagonalization%jacobi_threshold, &
                          error=error)
         CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)
         CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                       qs_env%rho%rho_ao(ispin)%matrix,&
                                       error=error)
       END DO

    ELSE IF (density_guess == no_guess) THEN
      
    ELSE
       CALL stop_program(routineP,"An invalid keyword for the initial density guess was specified")
    END IF

    !
    ! once a run, we check the new bcsr code on real data, if things go wrong the user can
    ! disable this in the &QS section using the keyword CHECK_BCSR_CODE
    !
    IF (scf_control%use_ot) THEN
       IF (qs_env%dft_control%qs_control%check_bcsr_code) THEN
          ionode = logger%para_env%mepos==logger%para_env%source
          output_unit = -1
          IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)
          IF (output_unit>0) THEN
             WRITE(output_unit,'(T2,A)') ""
             WRITE(output_unit,'(T2,A)') "Testing the BCSR code ... " 
             WRITE(output_unit,'(T2,A)') "Report errors, and use the input keyword CHECK_BCSR_CODE to circumvent them"
          ENDIF
          CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
          CALL dbcsr_check(s_sparse(1)%matrix, mo_coeff,&
               qs_env%dft_control%qs_control%bcsr_code, error)
          qs_env%dft_control%qs_control%check_bcsr_code=.FALSE.
       ENDIF
    ENDIF


    CALL timestop(handle)

  END SUBROUTINE calculate_first_density_matrix

END MODULE qs_initial_guess
