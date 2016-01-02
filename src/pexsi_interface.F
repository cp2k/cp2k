!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the PEXSI library, providing wrappers for all PEXSI 
!>        routines that are called inside CP2K. Tested for PEXSI versions 
!>        v0.7.3 - v0.9.0
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
#include "./base/base_uses.f90"

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
! *****************************************************************************
  SUBROUTINE cp_pexsi_set_options(pexsi_options, temperature, gap, deltaE, numPole, & 
                                  isInertiaCount, maxPEXSIIter, muMin0, muMax0, mu0, & 
                                  muInertiaTolerance, muInertiaExpansion, &
                                  muPEXSISafeGuard, numElectronPEXSITolerance, &
                                  matrixType, isSymbolicFactorize, ordering, & 
                                  npSymbFact, verbosity)

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
    MARK_USED(pexsi_options)
    MARK_USED(temperature)
    MARK_USED(gap)
    MARK_USED(deltaE)
    MARK_USED(numPole)
    MARK_USED(isInertiaCount)
    MARK_USED(maxPEXSIIter)
    MARK_USED(muMin0)
    MARK_USED(muMax0)
    MARK_USED(mu0)
    MARK_USED(muInertiaTolerance)
    MARK_USED(muInertiaExpansion)
    MARK_USED(muPEXSISafeGuard)
    MARK_USED(numElectronPEXSITolerance)
    MARK_USED(matrixType)
    MARK_USED(isSymbolicFactorize)
    MARK_USED(ordering)
    MARK_USED(npSymbFact)
    MARK_USED(verbosity)
    CPABORT("Requires linking to the PEXSI library.")
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
! *****************************************************************************
  SUBROUTINE cp_pexsi_get_options(pexsi_options, temperature, gap, deltaE, numPole, &
                                  isInertiaCount, maxPEXSIIter, muMin0, muMax0, mu0, &
                                  muInertiaTolerance, muInertiaExpansion, &
                                  muPEXSISafeGuard, numElectronPEXSITolerance, &
                                  matrixType, isSymbolicFactorize, ordering, & 
                                  npSymbFact, verbosity)
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
   MARK_USED(pexsi_options)
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
   CPABORT("Requires linking to the PEXSI library.")
#endif
  END SUBROUTINE cp_pexsi_get_options

! *****************************************************************************
!> \brief ...
!> \param pexsi_options ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_set_default_options(pexsi_options)
    TYPE(cp_pexsi_options), INTENT(OUT)      :: pexsi_options

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_set_default_options', &
      routineP = moduleN//':'//routineN

#ifdef  __LIBPEXSI
    CALL f_ppexsi_set_default_options(pexsi_options%options)
#else
    CPABORT("Requires linking to the PEXSI library.")
#endif
  END SUBROUTINE cp_pexsi_set_default_options

! *****************************************************************************
!> \brief ...
!> \param comm ...
!> \param numProcRow ...
!> \param numProcCol ...
!> \param outputFileIndex ...
!> \retval cp_pexsi_plan_initialize ...
! *****************************************************************************
  FUNCTION cp_pexsi_plan_initialize(comm, numProcRow, numProcCol, outputFileIndex)
    INTEGER, INTENT(IN)                      :: comm, numProcRow, numProcCol, &
                                                outputFileIndex
    INTEGER(KIND=int_8)                      :: cp_pexsi_plan_initialize

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_plan_initialize', &
      routineP = moduleN//':'//routineN

#ifdef  __LIBPEXSI
    INTEGER                                  :: info, handle
    CALL timeset(routineN,handle)
    cp_pexsi_plan_initialize = f_ppexsi_plan_initialize(comm, numProcRow, &
                                         numProcCol, outputFileIndex, info)
    IF(info .NE. 0)&
       CPABORT("Pexsi returned an error. Consider logPEXSI0 for details.")
    CALL timestop(handle)
#else
    MARK_USED(comm)
    MARK_USED(numProcRow)
    MARK_USED(numProcCol)
    MARK_USED(outputFileIndex)
    CPABORT("Requires linking to the PEXSI library.")
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
! *****************************************************************************
  SUBROUTINE cp_pexsi_load_real_symmetric_hs_matrix(plan,pexsi_options,nrows,nnz, &
                                                    nnzLocal,numColLocal,colptrLocal, &
                                                    rowindLocal,HnzvalLocal,isSIdentity, &
                                                    SnzvalLocal)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    INTEGER, INTENT(IN)                      :: nrows, nnz, nnzLocal, &
                                                numColLocal, colptrLocal(*), &
                                                rowindLocal(*)
    REAL(KIND=real_8), INTENT(IN)            :: HnzvalLocal(*)
    INTEGER, INTENT(IN)                      :: isSIdentity
    REAL(KIND=real_8), INTENT(IN)            :: SnzvalLocal(*)

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_pexsi_load_real_symmetric_hs_matrix', &
      routineP = moduleN//':'//routineN


#ifdef  __LIBPEXSI
    INTEGER                                  :: handle, info

    CALL timeset(routineN,handle)
    CALL f_ppexsi_load_real_symmetric_hs_matrix(plan,pexsi_options%options,nrows,nnz,nnzLocal, &
                                                numColLocal,colptrLocal,rowindLocal, &
                                                HnzvalLocal,isSIdentity,SnzvalLocal,info)
    IF(info .NE. 0)&
       CPABORT("Pexsi returned an error. Consider logPEXSI0 for details.")
    CALL timestop(handle)
#else
    MARK_USED(plan)
    MARK_USED(pexsi_options)
    MARK_USED(nrows)
    MARK_USED(nnz)
    MARK_USED(nnzLocal)
    MARK_USED(numColLocal)
    MARK_USED(isSIdentity)
    CPABORT("Requires linking to the PEXSI library.")

    ! MARK_USED macro does not work on assumed shape variables
    IF(.FALSE.) THEN; DO
         IF(colptrLocal(1)>rowindLocal(1).OR.HnzvalLocal(1)>SnzvalLocal(1)) EXIT
    ENDDO; ENDIF
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
! *****************************************************************************
  SUBROUTINE cp_pexsi_dft_driver(plan,pexsi_options,numElectronExact,muPEXSI, &
                                 numElectronPEXSI,muMinInertia,muMaxInertia, &
                                 numTotalInertiaIter,numTotalPEXSIIter)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    REAL(KIND=real_8), INTENT(IN)            :: numElectronExact
    REAL(KIND=real_8), INTENT(out)           :: muPEXSI, numElectronPEXSI, &
                                                muMinInertia, muMaxInertia
    INTEGER, INTENT(out)                     :: numTotalInertiaIter, &
                                                numTotalPEXSIIter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_dft_driver', &
      routineP = moduleN//':'//routineN

#ifdef  __LIBPEXSI
    INTEGER                                  :: handle, info

    CALL timeset(routineN,handle)
    CALL f_ppexsi_dft_driver(plan,pexsi_options%options,numElectronExact,muPEXSI, &
                             numElectronPEXSI,muMinInertia,muMaxInertia, &
                             numTotalInertiaIter,numTotalPEXSIIter,info)
    IF(info .NE. 0)&
       CPABORT("Pexsi returned an error. Consider logPEXSI0 for details.")
    CALL timestop(handle)
#else
    MARK_USED(plan)
    MARK_USED(numelectronexact)
    MARK_USED(pexsi_options)
    ! assign intent-out arguments to silence compiler warnings
    muPEXSI = 0.0_real_8
    numElectronPEXSI = 0.0_real_8
    muMinInertia = 0.0_real_8
    muMaxInertia = 0.0_real_8
    numTotalInertiaIter = -1
    numTotalPEXSIIter = -1
    CPABORT("Requires linking to the PEXSI library.")
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
! *****************************************************************************
  SUBROUTINE cp_pexsi_retrieve_real_symmetric_dft_matrix(plan,DMnzvalLocal,EDMnzvalLocal, &
                                                         FDMnzvalLocal,totalEnergyH, &
                                                         totalEnergyS,totalFreeEnergy)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    REAL(KIND=real_8), INTENT(out) :: DMnzvalLocal(*), EDMnzvalLocal(*), &
      FDMnzvalLocal(*), totalEnergyH, totalEnergyS, totalFreeEnergy

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_pexsi_retrieve_real_symmetric_dft_matrix', &
      routineP = moduleN//':'//routineN


#ifdef  __LIBPEXSI
    INTEGER                                  :: handle, info

    CALL timeset(routineN,handle)
    CALL f_ppexsi_retrieve_real_symmetric_dft_matrix(plan,DMnzvalLocal,EDMnzvalLocal, &
                                                     FDMnzvalLocal,totalEnergyH, &
                                                     totalEnergyS,totalFreeEnergy,info)
    IF(info .NE. 0)&
       CPABORT("Pexsi returned an error. Consider logPEXSI0 for details.")
    CALL timestop(handle)
#else
    MARK_USED(plan)
    ! assign intent-out arguments to silence compiler warnings
    DMnzvalLocal(1) = 0.0_real_8
    EDMnzvalLocal(1) = 0.0_real_8
    FDMnzvalLocal(1) = 0.0_real_8
    totalEnergyH = 0.0_real_8
    totalEnergyS = 0.0_real_8
    totalFreeEnergy = 0.0_real_8

    CPABORT("Requires linking to the PEXSI library.")
#endif
  END SUBROUTINE cp_pexsi_retrieve_real_symmetric_dft_matrix

! *****************************************************************************
!> \brief ...
!> \param plan ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_plan_finalize(plan)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_plan_finalize', &
      routineP = moduleN//':'//routineN


#ifdef  __LIBPEXSI
    INTEGER                                  :: info, handle
    CALL timeset(routineN,handle)
    CALL f_ppexsi_plan_finalize(plan,info)
    IF(info .NE. 0)&
       CPABORT("Pexsi returned an error. Consider logPEXSI0 for details.")
    CALL timestop(handle)
#else
    MARK_USED(plan)
    CPABORT("Requires linking to the PEXSI library.")
#endif
  END SUBROUTINE

END MODULE pexsi_interface
