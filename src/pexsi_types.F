!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Environment storing all data that is needed in order to call the DFT 
!>        driver of the PEXSI library with data from the linear scaling quickstep 
!>        SCF environment, mainly parameters and intermediate data for the matrix
!>        conversion between DBCSR and CSR format.
!> \par History
!>       2014.11 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************

MODULE pexsi_types
  USE bibliography,                    ONLY: Lin2009,&
                                             Lin2013,&
                                             cite_reference
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_scale,&
                                             cp_dbcsr_type,&
                                             csr_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_dims_create,&
                                             mp_environ
  USE pexsi_interface,                 ONLY: cp_pexsi_get_options,&
                                             cp_pexsi_options,&
                                             cp_pexsi_plan_finalize,&
                                             cp_pexsi_plan_initialize,&
                                             cp_pexsi_set_options
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pexsi_types'

  LOGICAL, PARAMETER, PRIVATE          :: careful_mod = .FALSE.

  INTEGER, PARAMETER, PUBLIC           :: cp2k_to_pexsi = 1, pexsi_to_cp2k = 2

  PUBLIC :: lib_pexsi_env, lib_pexsi_init, lib_pexsi_finalize, &
            convert_nspin_cp2k_pexsi

! *****************************************************************************
!> \brief All PEXSI related data
!> \param options                       PEXSI options
!> \param plan                          PEXSI plan
!> \param mp_group                      message-passing group ID
!> \param mp_dims                       dimensions of the MPI cartesian grid used 
!>                                      for PEXSI
!> \param num_ranks_per_pole            number of MPI ranks per pole in PEXSI
!> \param kTS                           entropic energy contribution
!> \param matrix_w                      energy-weighted density matrix as needed 
!>                                      for the forces
!> \param csr_mat*                      intermediate matrices in CSR format
!> \param dbcsr_template_matrix_sym     Symmetric template matrix fixing DBCSR
!>                                      sparsity pattern
!> \param dbcsr_template_matrix_nonsym  Nonsymmetric template matrix fixing
!>                                      DBCSR sparsity pattern
!> \param csr_sparsity                  DBCSR matrix defining CSR sparsity
!> \param csr_screening                 whether distance screening should be
!>                                      applied to CSR matrices
!> \param max_ev_vector                 eigenvector corresponding to the largest 
!>                                      energy eigenvalue,
!>                                      returned by the Arnoldi method used to 
!>                                      determine the spectral radius deltaE
!> \param nspin                         number of spins
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  TYPE lib_pexsi_env
    TYPE(cp_dbcsr_type)                      :: dbcsr_template_matrix_sym, &
                                                dbcsr_template_matrix_nonsym
    TYPE(csr_type)                           :: csr_mat_p, csr_mat_ks, csr_mat_s, &
                                                csr_mat_E, csr_mat_F
    TYPE(cp_pexsi_options)                   :: options
    REAL(KIND=dp), DIMENSION(:), POINTER     :: kTS => NULL()

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_w => NULL()

    INTEGER(kind=int_8)                      :: plan
    INTEGER                                  :: nspin, num_ranks_per_pole, mp_group

    TYPE(cp_dbcsr_type), DIMENSION(:), &
      POINTER                                :: max_ev_vector
    TYPE(cp_dbcsr_type)                      :: csr_sparsity
    INTEGER, DIMENSION(2)                    :: mp_dims

    LOGICAL                                  :: csr_screening
  END TYPE lib_pexsi_env

CONTAINS

! *****************************************************************************
!> \brief Initialize PEXSI
!> \param pexsi_env All data needed by PEXSI 
!> \param mp_group message-passing group ID
!> \param nspin number of spins
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE lib_pexsi_init(pexsi_env, mp_group, nspin, error)
    TYPE(lib_pexsi_env), INTENT(INOUT)       :: pexsi_env
    INTEGER, INTENT(IN)                      :: mp_group, nspin
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lib_pexsi_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, mynode, &
                                                npSymbFact, numnodes, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL timeset(routineN,handle)

    CALL cite_reference(Lin2009)
    CALL cite_reference(Lin2013)

    pexsi_env%mp_group = mp_group
    CALL mp_environ(numnodes, mynode, pexsi_env%mp_group)

    ! Use all nodes available per pole by default or if the user tries to use 
    ! more nodes than MPI size
    IF ((pexsi_env%num_ranks_per_pole .GT. numnodes) &
       .OR. (pexsi_env%num_ranks_per_pole .EQ. 0)) THEN
         pexsi_env%num_ranks_per_pole = numnodes
    ENDIF

    ! Find num_ranks_per_pole from user input MIN_RANKS_PER_POLE s.t. num_ranks_per_pole
    ! is the smallest number greater or equal to min_ranks_per_pole that divides 
    ! MPI size without remainder.
    DO WHILE(MOD(numnodes, pexsi_env%num_ranks_per_pole) .NE. 0)
      pexsi_env%num_ranks_per_pole = pexsi_env%num_ranks_per_pole + 1
    ENDDO

    CALL cp_pexsi_get_options(pexsi_env%options,npSymbFact=npSymbFact,error=error)
    IF ((npSymbFact .GT. pexsi_env%num_ranks_per_pole) .OR. (npSymbFact .EQ. 0)) THEN
      ! Use maximum possible number of ranks for symbolic factorization
      CALL cp_pexsi_set_options(pexsi_env%options, npSymbFact=pexsi_env%num_ranks_per_pole, &
                                error=error)
    ENDIF

    ! Create dimensions for MPI cartesian grid for PEXSI
    pexsi_env%mp_dims=0
    CALL mp_dims_create(pexsi_env%num_ranks_per_pole, pexsi_env%mp_dims)

    ! allocations with nspin
    pexsi_env%nspin=nspin
    ALLOCATE(pexsi_env%kTS(nspin))
    ALLOCATE(pexsi_env%max_ev_vector(nspin))
    ALLOCATE(pexsi_env%matrix_w(nspin))
    DO ispin = 1, pexsi_env%nspin 
      ALLOCATE(pexsi_env%matrix_w(ispin)%matrix)
    ENDDO

    ! Initialize PEXSI
    pexsi_env%plan = cp_pexsi_plan_initialize(pexsi_env%mp_group,pexsi_env%mp_dims(1),&
                                              pexsi_env%mp_dims(2),mynode,error)

    ! Print PEXSI infos
    IF(unit_nr>0) CALL print_pexsi_info(pexsi_env, unit_nr, error)

    CALL timestop(handle)
  END SUBROUTINE lib_pexsi_init

! *****************************************************************************
!> \brief Release all PEXSI data
!> \param pexsi_env ...
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE lib_pexsi_finalize(pexsi_env, error)
    TYPE(lib_pexsi_env), INTENT(INOUT)       :: pexsi_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lib_pexsi_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin

    CALL timeset(routineN,handle)
    CALL cp_pexsi_plan_finalize(pexsi_env%plan, error)
    DEALLOCATE(pexsi_env%kTS)
    DEALLOCATE(pexsi_env%max_ev_vector)
    DO ispin = 1, pexsi_env%nspin
      DEALLOCATE(pexsi_env%matrix_w(ispin)%matrix)
    ENDDO
    DEALLOCATE(pexsi_env%matrix_w) 
    CALL timestop(handle)
  END SUBROUTINE lib_pexsi_finalize

! *****************************************************************************
!> \brief Scale various quantities with factors of 2. This converts spin restricted
!>        DFT calculations (PEXSI) to the unrestricted case (as is the case where
!>        the density matrix method is called in the linear scaling code).
!> \param[in] direction ...
!> \param[in] nspin ...
!> \param[in,out] numElectron ...
!> \param matrix_p ...
!> \param matrix_w ...
!> \param kTS ...
!> \param error ...
!> \par History
!>      01.2015 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE convert_nspin_cp2k_pexsi(direction, nspin, numElectron, matrix_p, matrix_w, kTS, error)
    INTEGER, INTENT(IN)                      :: direction, nspin
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: numElectron
    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      OPTIONAL                               :: matrix_p
    TYPE(cp_dbcsr_p_type), INTENT(INOUT), &
      OPTIONAL                               :: matrix_w
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: kTS
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'convert_nspin_cp2k_pexsi', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: scaling

    CALL timeset(routineN,handle)

    SELECT CASE(direction)
    CASE(cp2k_to_pexsi)
      scaling = 2.0_dp
    CASE(pexsi_to_cp2k)
      scaling = 0.5_dp
    END SELECT

    IF (PRESENT(numElectron)) numElectron = scaling*numElectron
    IF (PRESENT(matrix_p)) CALL cp_dbcsr_scale(matrix_p, scaling, error=error)
    IF (PRESENT(matrix_w)) CALL cp_dbcsr_scale(matrix_w%matrix, scaling, error=error)
    IF (PRESENT(kTS)) kTS = scaling*kTS

    CALL timestop(handle)
  END SUBROUTINE convert_nspin_cp2k_pexsi

! *****************************************************************************
!> \brief Print relevant options of PEXSI
!> \param pexsi_env ...
!> \param unit_nr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE print_pexsi_info(pexsi_env,unit_nr, error)
    TYPE(lib_pexsi_env), INTENT(IN)          :: pexsi_env
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: mynode, npSymbFact, numnodes, &
                                                numPole, ordering
    REAL(KIND=dp) :: gap, muInertiaExpansion, muInertiaTolerance, muMax0, &
      muMin0, muPEXSISafeGuard, numElectronPEXSITolerance, temperature

    CALL mp_environ(numnodes, mynode, pexsi_env%mp_group)

    CALL cp_pexsi_get_options(pexsi_env%options, temperature=temperature, gap=gap, &
           numPole=numPole, muMin0=muMin0, muMax0=muMax0, muInertiaTolerance=&
           muInertiaTolerance, muInertiaExpansion=muInertiaExpansion, &
           muPEXSISafeGuard=muPEXSISafeGuard, numElectronPEXSITolerance =&
           numElectronPEXSITolerance, ordering=ordering, npSymbFact=npSymbFact, &
           error = error)

    WRITE(unit_nr,'(/A)') " PEXSI| Initialized with parameters"
    WRITE(unit_nr,'(A,T61,E20.3)') " PEXSI|   Electronic temperature", temperature
    WRITE(unit_nr,'(A,T61,E20.3)') " PEXSI|   Spectral gap", gap
    WRITE(unit_nr,'(A,T61,I20)') " PEXSI|   Number of poles", numPole
    WRITE(unit_nr,'(A,T61,E20.3)') " PEXSI|   Convergence tolerance in number of electrons", &
                                   numElectronPEXSITolerance
    WRITE(unit_nr,'(A,T61,E20.3)') " PEXSI|   Convergence tolerance for inertia counting in mu", &
                                   muInertiaTolerance
    WRITE(unit_nr,'(A,T61,E20.3)') " PEXSI|   Restart tolerance for inertia counting in mu", &
                                   muPEXSISafeGuard
    WRITE(unit_nr,'(A,T61,E20.3)') " PEXSI|   Expansion of mu interval for inertia counting", &
                                   muInertiaExpansion

    WRITE(unit_nr,'(/A)') " PEXSI| Parallelization scheme"
    WRITE(unit_nr,'(A,T61,I20)') " PEXSI|   Number of poles processed in parallel", &
                                 numnodes/pexsi_env%num_ranks_per_pole
    WRITE(unit_nr,'(A,T61,I20)') " PEXSI|   Number of processors per pole", &
                                 pexsi_env%num_ranks_per_pole
    WRITE(unit_nr,'(A,T71,I5,T76,I5)') " PEXSI|   Process grid dimensions", &
                                       pexsi_env%mp_dims(1), pexsi_env%mp_dims(2)
    SELECT CASE (ordering)
      CASE(0) 
        WRITE(unit_nr,'(A,T61,A20)') " PEXSI|   Ordering strategy", "parallel"
      CASE(1) 
        WRITE(unit_nr,'(A,T61,A20)') " PEXSI|   Ordering strategy", "sequential"
      CASE(2) 
        WRITE(unit_nr,'(A,T61,A20)') " PEXSI|   Ordering strategy", "multiple minimum degree"
    END SELECT
    IF(ordering .EQ. 0) WRITE(unit_nr,'(A,T61,I20/)') &
                        " PEXSI|   Number of processors for symbolic factorization", npSymbFact

  END SUBROUTINE print_pexsi_info
END MODULE pexsi_types
