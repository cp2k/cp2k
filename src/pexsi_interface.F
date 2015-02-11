!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the PEXSI library v0.7.3, providing wrappers for all PEXSI 
!>        routines that are called inside CP2K.
!> \par History
!>       2014.12 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
MODULE pexsi_interface

#ifdef __LIBPEXSI
  USE f_ppexsi_interface,              ONLY: f_ppexsi_dft_driver,&
                                             f_ppexsi_load_real_symmetric_hs_matrix,&
                                             f_ppexsi_options,&
                                             f_ppexsi_plan_finalize,&
                                             f_ppexsi_plan_initialize,&
                                             f_ppexsi_retrieve_real_symmetric_dft_matrix,&
                                             f_ppexsi_set_default_options
#endif
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pexsi_interface'

  PUBLIC :: cp_pexsi_options, cp_pexsi_plan_initialize, & 
            cp_pexsi_load_real_symmetric_hs_matrix, cp_pexsi_dft_driver, &
            cp_pexsi_retrieve_real_symmetric_dft_matrix, cp_pexsi_plan_finalize,&
            cp_pexsi_set_options, cp_pexsi_get_options, cp_pexsi_set_default_options

  TYPE cp_pexsi_options
      PRIVATE
#ifdef __LIBPEXSI
      TYPE(f_ppexsi_options) :: options
#else
      INTEGER :: unused = -1
#endif
  END TYPE cp_pexsi_options

CONTAINS

! *****************************************************************************
!> \brief Set PEXSI internal options
!> \param pexsi_options ...
!> \param temperature ...
!> \param gap ...
!> \param deltaE ...
!> \param numPole ...
!> \param isInertiaCount ...
!> \param maxPEXSIIter ...
!> \param muMin0 ...
!> \param muMax0 ...
!> \param mu0 ...
!> \param muInertiaTolerance ...
!> \param muInertiaExpansion ...
!> \param muPEXSISafeGuard ...
!> \param numElectronPEXSITolerance ...
!> \param matrixType ...
!> \param isSymbolicFactorize ...
!> \param ordering ...
!> \param npSymbFact ...
!> \param verbosity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_set_options(pexsi_options, temperature, gap, deltaE, numPole, & 
                                  isInertiaCount, maxPEXSIIter, muMin0, muMax0, mu0, & 
                                  muInertiaTolerance, muInertiaExpansion, &
                                  muPEXSISafeGuard, numElectronPEXSITolerance, &
                                  matrixType, isSymbolicFactorize, ordering, & 
                                  npSymbFact, verbosity, error)

    TYPE(cp_pexsi_options), INTENT(INOUT)    :: pexsi_options
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: temperature, gap, deltaE
    INTEGER, INTENT(IN), OPTIONAL            :: numPole, isInertiaCount, &
                                                maxPEXSIIter
    REAL(KIND=real_8), INTENT(IN), OPTIONAL :: muMin0, muMax0, mu0, &
      muInertiaTolerance, muInertiaExpansion, muPEXSISafeGuard, &
      numElectronPEXSITolerance
    INTEGER, INTENT(IN), OPTIONAL            :: matrixType, &
                                                isSymbolicFactorize, &
                                                ordering, npSymbFact, &
                                                verbosity
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_set_options', &
      routineP = moduleN//':'//routineN

#ifdef __LIBPEXSI
    IF(PRESENT(temperature)) pexsi_options%options%temperature = temperature
    IF(PRESENT(gap)) pexsi_options%options%gap = gap
    IF(PRESENT(deltaE)) pexsi_options%options%deltaE = deltaE
    IF(PRESENT(numPole)) pexsi_options%options%numPole = numPole
    IF(PRESENT(isInertiaCount)) pexsi_options%options%isInertiaCount = isInertiaCount
    IF(PRESENT(maxPEXSIIter)) pexsi_options%options%maxPEXSIIter = maxPEXSIIter
    IF(PRESENT(muMin0)) pexsi_options%options%muMin0 = muMin0
    IF(PRESENT(muMax0)) pexsi_options%options%muMax0 = muMax0
    IF(PRESENT(mu0)) pexsi_options%options%mu0 = mu0
    IF(PRESENT(muInertiaTolerance)) & 
      pexsi_options%options%muInertiaTolerance = muInertiaTolerance
    IF(PRESENT(muInertiaExpansion)) & 
      pexsi_options%options%muInertiaExpansion = muInertiaExpansion
    IF(PRESENT(muPEXSISafeGuard)) & 
      pexsi_options%options%muPEXSISafeGuard = muPEXSISafeGuard
    IF(PRESENT(numElectronPEXSITolerance)) & 
      pexsi_options%options%numElectronPEXSITolerance = numElectronPEXSITolerance
    IF(PRESENT(matrixType)) pexsi_options%options%matrixType = matrixType
    IF(PRESENT(isSymbolicFactorize)) & 
      pexsi_options%options%isSymbolicFactorize = isSymbolicFactorize
    IF(PRESENT(ordering)) pexsi_options%options%ordering = ordering
    IF(PRESENT(npSymbFact)) pexsi_options%options%npSymbFact = npSymbFact
    IF(PRESENT(verbosity)) pexsi_options%options%verbosity = verbosity
#else
    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE cp_pexsi_set_options

! *****************************************************************************
!> \brief Access PEXSI internal options
!> \param pexsi_options ...
!> \param temperature ...
!> \param gap ...
!> \param deltaE ...
!> \param numPole ...
!> \param isInertiaCount ...
!> \param maxPEXSIIter ...
!> \param muMin0 ...
!> \param muMax0 ...
!> \param mu0 ...
!> \param muInertiaTolerance ...
!> \param muInertiaExpansion ...
!> \param muPEXSISafeGuard ...
!> \param numElectronPEXSITolerance ...
!> \param matrixType ...
!> \param isSymbolicFactorize ...
!> \param ordering ...
!> \param npSymbFact ...
!> \param verbosity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_get_options(pexsi_options, temperature, gap, deltaE, numPole, &
                                  isInertiaCount, maxPEXSIIter, muMin0, muMax0, mu0, &
                                  muInertiaTolerance, muInertiaExpansion, &
                                  muPEXSISafeGuard, numElectronPEXSITolerance, &
                                  matrixType, isSymbolicFactorize, ordering, & 
                                  npSymbFact, verbosity, error)
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL :: temperature, gap, deltaE
    INTEGER, INTENT(OUT), OPTIONAL           :: numPole, isInertiaCount, &
                                                maxPEXSIIter
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL :: muMin0, muMax0, mu0, &
      muInertiaTolerance, muInertiaExpansion, muPEXSISafeGuard, &
      numElectronPEXSITolerance
    INTEGER, INTENT(OUT), OPTIONAL           :: matrixType, &
                                                isSymbolicFactorize, &
                                                ordering, npSymbFact, &
                                                verbosity
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_get_options', &
      routineP = moduleN//':'//routineN

#ifdef __LIBPEXSI
    IF(PRESENT(temperature)) temperature = pexsi_options%options%temperature
    IF(PRESENT(gap)) gap = pexsi_options%options%gap
    IF(PRESENT(deltaE)) deltaE = pexsi_options%options%deltaE
    IF(PRESENT(numPole)) numPole = pexsi_options%options%numPole
    IF(PRESENT(isInertiaCount)) isInertiaCount = pexsi_options%options%isInertiaCount
    IF(PRESENT(maxPEXSIIter)) maxPEXSIIter = pexsi_options%options%maxPEXSIIter
    IF(PRESENT(muMin0)) muMin0 = pexsi_options%options%muMin0
    IF(PRESENT(muMax0)) muMax0 = pexsi_options%options%muMax0
    IF(PRESENT(mu0)) mu0 = pexsi_options%options%mu0
    IF(PRESENT(muInertiaTolerance)) &
      muInertiaTolerance = pexsi_options%options%muInertiaTolerance
    IF(PRESENT(muInertiaExpansion)) &
      muInertiaExpansion = pexsi_options%options%muInertiaExpansion
    IF(PRESENT(muPEXSISafeGuard)) &
      muPEXSISafeGuard = pexsi_options%options%muPEXSISafeGuard
    IF(PRESENT(numElectronPEXSITolerance)) &
      numElectronPEXSITolerance = pexsi_options%options%numElectronPEXSITolerance
    IF(PRESENT(matrixType)) matrixType = pexsi_options%options%matrixType
    IF(PRESENT(isSymbolicFactorize)) & 
      isSymbolicFactorize = pexsi_options%options%isSymbolicFactorize
    IF(PRESENT(ordering)) ordering = pexsi_options%options%ordering
    IF(PRESENT(npSymbFact)) npSymbFact = pexsi_options%options%npSymbFact
    IF(PRESENT(verbosity)) verbosity = pexsi_options%options%verbosity
#else
   ! assign intent-out arguments to silence compiler warnings
   IF(PRESENT(temperature)) temperature = 0.0_real_8
   IF(PRESENT(gap)) gap = 0.0_real_8
   IF(PRESENT(deltaE)) deltaE = 0.0_real_8
   IF(PRESENT(numPole)) numPole = -1
   IF(PRESENT(isInertiaCount)) isInertiaCount = -1
   IF(PRESENT(maxPEXSIIter)) maxPEXSIIter = -1
   IF(PRESENT(muMin0)) muMin0 = 0.0_real_8
   IF(PRESENT(muMax0)) muMax0 = 0.0_real_8
   IF(PRESENT(mu0)) mu0 = 0.0_real_8
   IF(PRESENT(muInertiaTolerance)) muInertiaTolerance = 0.0_real_8
   IF(PRESENT(muInertiaExpansion)) muInertiaExpansion = 0.0_real_8
   IF(PRESENT(muPEXSISafeGuard)) muPEXSISafeGuard = 0.0_real_8
   IF(PRESENT(numElectronPEXSITolerance)) numElectronPEXSITolerance = 0.0_real_8
   IF(PRESENT(matrixType)) matrixType = -1
   IF(PRESENT(isSymbolicFactorize)) isSymbolicFactorize = -1
   IF(PRESENT(ordering)) ordering = -1
   IF(PRESENT(npSymbFact)) npSymbFact = -1
   IF(PRESENT(verbosity)) verbosity = -1
   CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE cp_pexsi_get_options

! *****************************************************************************
!> \brief ...
!> \param pexsi_options ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_set_default_options(pexsi_options, error)
    TYPE(cp_pexsi_options), INTENT(OUT)      :: pexsi_options
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_set_default_options', &
      routineP = moduleN//':'//routineN

#ifdef  __LIBPEXSI
    CALL f_ppexsi_set_default_options(pexsi_options%options)
#else
    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE cp_pexsi_set_default_options

! *****************************************************************************
!> \brief ...
!> \param comm ...
!> \param numProcRow ...
!> \param numProcCol ...
!> \param outputFileIndex ...
!> \param error ...
!> \retval cp_pexsi_plan_initialize ...
! *****************************************************************************
  FUNCTION cp_pexsi_plan_initialize(comm, numProcRow, numProcCol, outputFileIndex, error)
    INTEGER, INTENT(IN)                      :: comm, numProcRow, numProcCol, &
                                                outputFileIndex
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8)                      :: cp_pexsi_plan_initialize

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_plan_initialize', &
      routineP = moduleN//':'//routineN

#ifdef  __LIBPEXSI
    INTEGER                                  :: info, handle
    CALL timeset(routineN,handle)
    cp_pexsi_plan_initialize = f_ppexsi_plan_initialize(comm, numProcRow, &
                                         numProcCol, outputFileIndex, info)
    CALL cp_assert(info .EQ. 0, cp_fatal_level, cp_internal_error, routineP,&
                   "Pexsi returned an error. Consider logPEXSI0 for details.",&
                   only_ionode = .TRUE., error=error)
    CALL timestop(handle)
#else
    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END FUNCTION cp_pexsi_plan_initialize

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param pexsi_options ...
!> \param nrows ...
!> \param nnz ...
!> \param nnzLocal ...
!> \param numColLocal ...
!> \param colptrLocal ...
!> \param rowindLocal ...
!> \param HnzvalLocal ...
!> \param isSIdentity ...
!> \param SnzvalLocal ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_load_real_symmetric_hs_matrix(plan,pexsi_options,nrows,nnz, &
                                                    nnzLocal,numColLocal,colptrLocal, &
                                                    rowindLocal,HnzvalLocal,isSIdentity, &
                                                    SnzvalLocal,error)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    INTEGER, INTENT(IN)                      :: nrows, nnz, nnzLocal, &
                                                numColLocal, colptrLocal(*), &
                                                rowindLocal(*)
    REAL(KIND=real_8), INTENT(IN)            :: HnzvalLocal(*)
    INTEGER, INTENT(IN)                      :: isSIdentity
    REAL(KIND=real_8), INTENT(IN)            :: SnzvalLocal(*)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_pexsi_load_real_symmetric_hs_matrix', &
      routineP = moduleN//':'//routineN


#ifdef  __LIBPEXSI
    INTEGER                                  :: handle, info

    CALL timeset(routineN,handle)
    CALL f_ppexsi_load_real_symmetric_hs_matrix(plan,pexsi_options%options,nrows,nnz,nnzLocal, &
                                                numColLocal,colptrLocal,rowindLocal, &
                                                HnzvalLocal,isSIdentity,SnzvalLocal,info)
    CALL cp_assert(info .EQ. 0, cp_fatal_level, cp_internal_error, routineP, &
                   "Pexsi returned an error. Consider logPEXSI0 for details.", &
                   only_ionode = .TRUE., error=error)
    CALL timestop(handle)
#else
    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE cp_pexsi_load_real_symmetric_hs_matrix

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param pexsi_options ...
!> \param numElectronExact ...
!> \param muPEXSI ...
!> \param numElectronPEXSI ...
!> \param muMinInertia ...
!> \param muMaxInertia ...
!> \param numTotalInertiaIter ...
!> \param numTotalPEXSIIter ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_dft_driver(plan,pexsi_options,numElectronExact,muPEXSI, &
                                 numElectronPEXSI,muMinInertia,muMaxInertia, &
                                 numTotalInertiaIter,numTotalPEXSIIter,error)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    REAL(KIND=real_8), INTENT(IN)            :: numElectronExact
    REAL(KIND=real_8), INTENT(out)           :: muPEXSI, numElectronPEXSI, &
                                                muMinInertia, muMaxInertia
    INTEGER, INTENT(out)                     :: numTotalInertiaIter, &
                                                numTotalPEXSIIter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_dft_driver', &
      routineP = moduleN//':'//routineN

#ifdef  __LIBPEXSI
    INTEGER                                  :: handle, info

    CALL timeset(routineN,handle)
    CALL f_ppexsi_dft_driver(plan,pexsi_options%options,numElectronExact,muPEXSI, &
                             numElectronPEXSI,muMinInertia,muMaxInertia, &
                             numTotalInertiaIter,numTotalPEXSIIter,info)
    CALL cp_assert(info .EQ. 0, cp_fatal_level, cp_internal_error, routineP,&
                   "Pexsi returned an error. Consider logPEXSI0 for details.",&
                   only_ionode = .TRUE., error=error)
    CALL timestop(handle)
#else
    ! assign intent-out arguments to silence compiler warnings
    muPEXSI = 0.0_real_8
    numElectronPEXSI = 0.0_real_8
    muMinInertia = 0.0_real_8
    muMaxInertia = 0.0_real_8
    numTotalInertiaIter = -1
    numTotalPEXSIIter = -1
    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE cp_pexsi_dft_driver

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param DMnzvalLocal ...
!> \param EDMnzvalLocal ...
!> \param FDMnzvalLocal ...
!> \param totalEnergyH ...
!> \param totalEnergyS ...
!> \param totalFreeEnergy ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_retrieve_real_symmetric_dft_matrix(plan,DMnzvalLocal,EDMnzvalLocal, &
                                                         FDMnzvalLocal,totalEnergyH, &
                                                         totalEnergyS,totalFreeEnergy,error)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    REAL(KIND=real_8), INTENT(out) :: DMnzvalLocal(*), EDMnzvalLocal(*), &
      FDMnzvalLocal(*), totalEnergyH, totalEnergyS, totalFreeEnergy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_pexsi_retrieve_real_symmetric_dft_matrix', &
      routineP = moduleN//':'//routineN


#ifdef  __LIBPEXSI
    INTEGER                                  :: handle, info

    CALL timeset(routineN,handle)
    CALL f_ppexsi_retrieve_real_symmetric_dft_matrix(plan,DMnzvalLocal,EDMnzvalLocal, &
                                                     FDMnzvalLocal,totalEnergyH, &
                                                     totalEnergyS,totalFreeEnergy,info)
    CALL cp_assert(info .EQ. 0, cp_fatal_level, cp_internal_error, routineP,&
                   "Pexsi returned an error. Consider logPEXSI0 for details.",&
                   only_ionode = .TRUE., error=error)
    CALL timestop(handle)
#else
    ! assign intent-out arguments to silence compiler warnings
    DMnzvalLocal(1) = 0.0_real_8
    EDMnzvalLocal(1) = 0.0_real_8
    FDMnzvalLocal(1) = 0.0_real_8
    totalEnergyH = 0.0_real_8
    totalEnergyS = 0.0_real_8
    totalFreeEnergy = 0.0_real_8

    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE cp_pexsi_retrieve_real_symmetric_dft_matrix

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_plan_finalize(plan,error)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_plan_finalize', &
      routineP = moduleN//':'//routineN


#ifdef  __LIBPEXSI
    INTEGER                                  :: info, handle
    CALL timeset(routineN,handle)
    CALL f_ppexsi_plan_finalize(plan,info)
    CALL cp_assert(info .EQ. 0, cp_fatal_level, cp_internal_error, routineP,&
                   "Pexsi returned an error. Consider logPEXSI0 for details.",&
                   only_ionode = .TRUE., error=error)
    CALL timestop(handle)
#else
    CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error, routineN, & 
                    "Requires linking to the PEXSI library.",&
                    only_ionode = .TRUE., error=error)
#endif
  END SUBROUTINE

END MODULE pexsi_interface
