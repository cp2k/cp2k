!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Provide various population analyses and print the requested output
!>        information
!>
!> \author  Matthias Krack (MK)
!> \date    09.07.2010
!> \version 1.0
! *****************************************************************************

MODULE population_analyses
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_copy, cp_dbcsr_deallocate_matrix, cp_dbcsr_get_block_p, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_p_type, cp_dbcsr_set, cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix,&
                                             write_fm_with_basis_info
  USE cp_fm_diag,                      ONLY: cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_diag,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: nso
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'population_analyses'

  PUBLIC :: lowdin_population_analysis,&
            mulliken_population_analysis

CONTAINS

! *****************************************************************************
!> \brief Perform a Lowdin population analysis based on a symmetric
!>        orthogonalisation of the density matrix using S^(1/2)
!>
!> \param qs_env ...
!> \param output_unit ...
!> \param print_level ...
!> \param error ...
!> \date    06.07.2010
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE lowdin_population_analysis(qs_env,output_unit,print_level,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit, print_level
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lowdin_population_analysis', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: headline
    INTEGER                                  :: handle, ispin, ndep, nsgf, &
                                                nspin, stat
    LOGICAL                                  :: failure, print_gop
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: orbpop
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: sm_p, sm_s
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    TYPE(cp_fm_type), POINTER                :: fm_s_half, fm_work1, fm_work2
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atomic_kind_set)
    NULLIFY (qs_kind_set)
    NULLIFY (fmstruct)
    NULLIFY (fm_s_half)
    NULLIFY (fm_work1)
    NULLIFY (fm_work2)
    NULLIFY (matrix_p)
    NULLIFY (matrix_s)
    NULLIFY (orbpop)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (scf_control)
    NULLIFY (sm_p)
    NULLIFY (sm_s)
    NULLIFY (orbpop)
    NULLIFY (para_env)
    NULLIFY (blacs_env)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    rho=rho,&
                    scf_control=scf_control,&
                    para_env=para_env,&
                    blacs_env=blacs_env,&
                    error=error)

    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)

    sm_s => matrix_s(1)%matrix ! Overlap matrix in sparse format
    matrix_p => rho%rho_ao     ! Density matrices in sparse format

    nspin = SIZE(matrix_p)

    ! Get the total number of contracted spherical Gaussian basis functions
    CALL get_qs_kind_set(qs_kind_set,nsgf=nsgf)

    ! Provide an array to store the orbital populations for each spin
    ALLOCATE (orbpop(nsgf,nspin),STAT=stat)
    CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
    orbpop(:,:) = 0.0_dp

    ! Write headline
    IF (output_unit > 0) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "LOWDIN POPULATION ANALYSIS"
    END IF

    ! Provide full size work matrices
    CALL cp_fm_struct_create(fmstruct=fmstruct,&
                             para_env=para_env,&
                             context=blacs_env,&
                             nrow_global=nsgf,&
                             ncol_global=nsgf,&
                             error=error)
    CALL cp_fm_create(matrix=fm_s_half,&
                      matrix_struct=fmstruct,&
                      name="S^(1/2) MATRIX",&
                      error=error)
    CALL cp_fm_create(matrix=fm_work1,&
                      matrix_struct=fmstruct,&
                      name="FULL WORK MATRIX 1",&
                      error=error)
    headline = "SYMMETRICALLY ORTHOGONALISED DENSITY MATRIX"
    CALL cp_fm_create(matrix=fm_work2,&
                      matrix_struct=fmstruct,&
                      name=TRIM(headline),&
                      error=error)
    CALL cp_fm_struct_release(fmstruct=fmstruct,&
                              error=error)

    ! Build full S^(1/2) matrix (computationally expensive)
    CALL copy_dbcsr_to_fm(sm_s,fm_s_half,error=error)
    CALL cp_fm_power(fm_s_half,fm_work1,0.5_dp,scf_control%eps_eigval,ndep,error=error)
    CALL cp_assert((ndep == 0),cp_warning_level,cp_assertion_failed,routineP,&
                   "Overlap matrix exhibits linear dependencies. At least some "//&
                   "eigenvalues have been quenched.",only_ionode=.TRUE.)

    ! Build Lowdin population matrix for each spin
    DO ispin=1,nspin
      sm_p => matrix_p(ispin)%matrix ! Density matrix for spin ispin in sparse format
      ! Calculate S^(1/2)*P*S^(1/2) as a full matrix (Lowdin)
      CALL cp_dbcsr_sm_fm_multiply(sm_p,fm_s_half,fm_work1,nsgf,error=error)
      CALL cp_gemm(transa="N",&
                      transb="N",&
                      m=nsgf,&
                      n=nsgf,&
                      k=nsgf,&
                      alpha=1.0_dp,&
                      matrix_a=fm_s_half,&
                      matrix_b=fm_work1,&
                      beta=0.0_dp,&
                      matrix_c=fm_work2,&
                      error=error)
      IF (print_level > 2) THEN
        ! Write the full Lowdin population matrix
        IF (nspin > 1) THEN
          IF (ispin == 1) THEN
            fm_work2%name = TRIM(headline)//" FOR ALPHA SPIN"
          ELSE
            fm_work2%name = TRIM(headline)//" FOR BETA SPIN"
          END IF
        END IF
        CALL write_fm_with_basis_info(fm_work2,4,6,qs_env,para_env,&
                                      output_unit=output_unit,error=error)
      END IF
      CALL cp_fm_get_diag(fm_work2,orbpop(:,ispin),error)
    END DO ! next spin ispin

    ! Write atomic populations and charges
    IF (output_unit > 0) THEN
      print_gop = (print_level > 1) ! Print also orbital populations
      CALL write_orbpop(orbpop,atomic_kind_set,qs_kind_set,particle_set,output_unit,print_gop,error)
    END IF

    ! Release local working storage
    CALL cp_fm_release(matrix=fm_s_half,error=error)
    CALL cp_fm_release(matrix=fm_work1,error=error)
    CALL cp_fm_release(matrix=fm_work2,error=error)
    IF (ASSOCIATED(orbpop)) THEN
      DEALLOCATE (orbpop,STAT=stat)
      CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE lowdin_population_analysis

! *****************************************************************************
!> \brief Perform a Mulliken population analysis
!>
!> \param qs_env ...
!> \param output_unit ...
!> \param print_level ...
!> \param error ...
!> \date    10.07.2010
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE mulliken_population_analysis(qs_env,output_unit,print_level,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit, print_level
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mulliken_population_analysis', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: headline
    INTEGER                                  :: blk, handle, iatom, isgf, &
                                                ispin, jatom, jsgf, natom, &
                                                nsgf, nspin, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf_atom
    LOGICAL                                  :: failure, found, print_gop
    REAL(KIND=dp)                            :: ps
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: orbpop, p_block, ps_block, &
                                                s_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: sm_p, sm_ps, sm_s
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atomic_kind_set)
    NULLIFY (qs_kind_set)
    NULLIFY (matrix_p)
    NULLIFY (matrix_s)
    NULLIFY (orbpop)
    NULLIFY (particle_set)
    NULLIFY (ps_block)
    NULLIFY (p_block)
    NULLIFY (rho)
    NULLIFY (sm_p)
    NULLIFY (sm_ps)
    NULLIFY (sm_s)
    NULLIFY (s_block)
    NULLIFY (para_env)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    rho=rho,&
                    para_env=para_env,&
                    error=error)

    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)

    sm_s => matrix_s(1)%matrix ! Overlap matrix in sparse format
    matrix_p => rho%rho_ao     ! Density matrices in sparse format

    nspin = SIZE(matrix_p)

    ! Get the total number of contracted spherical Gaussian basis functions
    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
    ALLOCATE (first_sgf_atom(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    first_sgf_atom(:) = 0

    CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_atom, error=error)

    ! Provide an array to store the orbital populations for each spin
    ALLOCATE (orbpop(nsgf,nspin),STAT=stat)
    CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
    orbpop(:,:) = 0.0_dp

    ! Write headline
    IF (output_unit > 0) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "MULLIKEN POPULATION ANALYSIS"
    END IF

    ! Create a DBCSR work matrix, if needed
    IF (print_level > 2) THEN
      ALLOCATE (sm_ps,STAT=stat)
      CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
      headline = "MULLIKEN NET ATOMIC ORBITAL AND OVERLAP POPULATION MATRIX"
      CALL cp_dbcsr_init(sm_ps,error=error)
      CALL cp_dbcsr_copy(matrix_b=sm_ps,&
                         matrix_a=sm_s,&
                         name=TRIM(headline),&
                         error=error)
      CALL cp_dbcsr_set(sm_ps,0.0_dp,error=error)
    END IF

    ! Build Mulliken population matrix for each spin
    DO ispin=1,nspin
      sm_p => matrix_p(ispin)%matrix ! Density matrix for spin ispin in sparse format
      ! Calculate Hadamard product of P and S as sparse matrix (Mulliken)
      ! CALL cp_dbcsr_hadamard_product(sm_p,sm_s,sm_ps,error=error)
      CALL cp_dbcsr_iterator_start(iter,sm_s)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
        CALL cp_dbcsr_iterator_next_block(iter,iatom,jatom,s_block,blk)
        IF (.NOT.(ASSOCIATED(s_block))) CYCLE
        CALL cp_dbcsr_get_block_p(matrix=sm_p,&
                               row=iatom,&
                               col=jatom,&
                               block=p_block,&
                               found=found)
        IF (print_level > 2) THEN
          CALL cp_dbcsr_get_block_p(matrix=sm_ps,&
                                 row=iatom,&
                                 col=jatom,&
                                 block=ps_block,&
                                 found=found)
          CPPostcondition(ASSOCIATED(ps_block),cp_failure_level,routineP,error,failure)
        END IF

        sgfb = first_sgf_atom(jatom)
        DO jsgf=1,SIZE(s_block,2)
          DO isgf=1,SIZE(s_block,1)
            ps = p_block(isgf,jsgf)*s_block(isgf,jsgf)
            IF (ASSOCIATED(ps_block)) ps_block(isgf,jsgf) = ps
            orbpop(sgfb,ispin) = orbpop(sgfb,ispin) + ps
          END DO
          sgfb = sgfb + 1
        END DO
        IF (iatom /= jatom) THEN
          sgfa = first_sgf_atom(iatom)
          DO isgf=1,SIZE(s_block,1)
            DO jsgf=1,SIZE(s_block,2)
              ps = p_block(isgf,jsgf)*s_block(isgf,jsgf)
              orbpop(sgfa,ispin) = orbpop(sgfa,ispin) + ps
            END DO
            sgfa = sgfa + 1
          END DO
        END IF
      END DO
      CALL cp_dbcsr_iterator_stop(iter)

      IF (print_level > 2) THEN
        ! Write the full Mulliken net AO and overlap population matrix
        IF (nspin > 1) THEN
          IF (ispin == 1) THEN
            sm_ps%matrix%m%name = TRIM(headline)//" FOR ALPHA SPIN"
          ELSE
            sm_ps%matrix%m%name = TRIM(headline)//" FOR BETA SPIN"
          END IF
        END IF
        CALL cp_dbcsr_write_sparse_matrix(sm_ps,4,6,qs_env,para_env,&
                                       output_unit=output_unit,&
                                       error=error)
      END IF
    END DO

    CALL mp_sum(orbpop,para_env%group)

    ! Write atomic populations and charges
    IF (output_unit > 0) THEN
      print_gop = (print_level > 1) ! Print also orbital populations
      CALL write_orbpop(orbpop,atomic_kind_set, qs_kind_set,particle_set,output_unit,print_gop,error)
    END IF

    ! Release local working storage
    IF (ASSOCIATED(sm_ps)) CALL cp_dbcsr_deallocate_matrix(sm_ps,error=error)
    IF (ASSOCIATED(orbpop)) THEN
      DEALLOCATE (orbpop,STAT=stat)
      CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
    END IF
    IF (ALLOCATED(first_sgf_atom)) THEN
      DEALLOCATE (first_sgf_atom,STAT=stat)
      CPPostconditionNoFail((stat == 0),cp_warning_level,routineP,error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE mulliken_population_analysis

! *****************************************************************************
!> \brief Write atomic orbital populations and net atomic charges
!>
!> \param orbpop ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param output_unit ...
!> \param print_orbital_contributions ...
!> \param error ...
!> \date    07.07.2010
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_orbpop(orbpop,atomic_kind_set,qs_kind_set,particle_set,output_unit,&
                          print_orbital_contributions,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: orbpop
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: print_orbital_contributions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_orbpop', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    INTEGER                                  :: handle, iao, iatom, ikind, &
                                                iset, isgf, ishell, iso, l, &
                                                natom, nset, nsgf, nspin
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), DIMENSION(3)              :: sumorbpop, totsumorbpop
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (lshell)
    NULLIFY (nshell)
    NULLIFY (orb_basis_set)
    NULLIFY (sgf_symbol)

    CPPrecondition(ASSOCIATED(orbpop),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)

    nspin = SIZE(orbpop,2)

    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

    ! Select and write headline
    IF (nspin == 1) THEN
      IF (print_orbital_contributions) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "# Orbital  AO symbol  Orbital population                Net charge"
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "#  Atom  Element  Kind  Atomic population                Net charge"
      END IF
    ELSE
      IF (print_orbital_contributions) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "# Orbital  AO symbol  Orbital population (alpha,beta)   Net charge  Spin moment"
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "#  Atom  Element  Kind  Atomic population (alpha,beta)   Net charge  Spin moment"
      END IF
    END IF

    totsumorbpop(:) = 0.0_dp

    iao = 1
    DO iatom=1,natom
      sumorbpop(:) = 0.0_dp
      NULLIFY (orb_basis_set)
      CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                           element_symbol=element_symbol,&
                           kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), orb_basis_set=orb_basis_set, zeff=zeff)
      IF (ASSOCIATED(orb_basis_set)) THEN
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=lshell,&
                               sgf_symbol=sgf_symbol)
        isgf = 1
        DO iset=1,nset
          DO ishell=1,nshell(iset)
            l = lshell(ishell,iset)
            DO iso=1,nso(l)
              IF (nspin == 1) THEN
                sumorbpop(1) = sumorbpop(1) + orbpop(iao,1)
                IF (print_orbital_contributions) THEN
                  IF (isgf == 1) WRITE (UNIT=output_unit,FMT="(A)") ""
                  WRITE (UNIT=output_unit,&
                         FMT="(T2,I9,2X,A2,1X,A,T30,F12.6)")&
                    iao,element_symbol,sgf_symbol(isgf),orbpop(iao,1)
                END IF
              ELSE
                sumorbpop(1:2) = sumorbpop(1:2) + orbpop(iao,1:2)
                sumorbpop(3) = sumorbpop(3) + orbpop(iao,1) - orbpop(iao,2)
                IF (print_orbital_contributions) THEN
                  IF (isgf == 1) WRITE (UNIT=output_unit,FMT="(A)") ""
                  WRITE (UNIT=output_unit,&
                         FMT="(T2,I9,2X,A2,1X,A,T29,2(1X,F12.6),T69,F12.6)")&
                    iao,element_symbol,sgf_symbol(isgf),orbpop(iao,1:2),&
                    orbpop(iao,1)-orbpop(iao,2)
                END IF
              END IF
              isgf = isgf + 1
              iao = iao + 1
            END DO
          END DO
        END DO
        IF (nspin == 1) THEN
          totsumorbpop(1) = totsumorbpop(1) + sumorbpop(1)
          totsumorbpop(3) = totsumorbpop(3) + zeff - sumorbpop(1)
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I7,5X,A2,2X,I6,T30,F12.6,T56,F12.6)")&
            iatom,element_symbol,ikind,sumorbpop(1),zeff-sumorbpop(1)
        ELSE
          totsumorbpop(1:2) = totsumorbpop(1:2) + sumorbpop(1:2)
          totsumorbpop(3) = totsumorbpop(3) + zeff - sumorbpop(1) - sumorbpop(2)
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I7,5X,A2,2X,I6,T29,4(1X,F12.6))")&
            iatom,element_symbol,ikind,sumorbpop(1:2),&
            zeff-sumorbpop(1)-sumorbpop(2),sumorbpop(3)
        END IF
      END IF ! atom has an orbital basis
    END DO ! next atom iatom

    ! Write total sums
    IF (print_orbital_contributions) WRITE (UNIT=output_unit,FMT="(A)") ""
    IF (nspin == 1) THEN
      WRITE (UNIT=output_unit,&
             FMT="(T2,A,T30,F12.6,T56,F12.6,/)")&
        "# Total charge",totsumorbpop(1),totsumorbpop(3)
    ELSE
      WRITE (UNIT=output_unit,&
             FMT="(T2,A,T29,4(1X,F12.6),/)")&
        "# Total charge and spin",totsumorbpop(1:2),totsumorbpop(3),&
        totsumorbpop(1)-totsumorbpop(2)
    END IF

    IF (output_unit > 0) CALL m_flush(output_unit)

    CALL timestop(handle)

  END SUBROUTINE write_orbpop

END MODULE population_analyses
