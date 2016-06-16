!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gaussian Process implementation
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_gaussprocess
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat_symm
   USE pao_types,                       ONLY: pao_env_type,&
                                              training_matrix_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_gaussprocess'

   PUBLIC ::pao_ml_gp_train, pao_ml_gp_predict

CONTAINS

! **************************************************************************************************
!> \brief Builds the covariance matrix
!> \param pao ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE pao_ml_gp_train(pao, training_matrices)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(training_matrix_type), DIMENSION(:), TARGET   :: training_matrices

      INTEGER                                            :: i, ikind, j, npoints
      REAL(dp), DIMENSION(:), POINTER                    :: idescr, jdescr
      TYPE(training_matrix_type), POINTER                :: training_matrix

      ! TODO this could be parallelized over ranks
      DO ikind = 1, SIZE(training_matrices)
         training_matrix => training_matrices(ikind)
         npoints = SIZE(training_matrix%inputs, 1) ! number of points
         CPASSERT(SIZE(training_matrix%outputs, 1) == npoints)
         CPASSERT(npoints > 0)

         !TODO proper output
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO|ML| Building covariance matrix for kind: ", &
            TRIM(training_matrix%kindname), " from ", npoints, "training points."

         ! build co-variance matrix A
         ALLOCATE (training_matrix%GP(npoints, npoints))
         DO i = 1, npoints
         DO j = i, npoints
            idescr => training_matrix%inputs(i, :)
            jdescr => training_matrix%inputs(j, :)
            training_matrix%GP(i, j) = kernel(pao%gp_scale, idescr, jdescr)
            training_matrix%GP(j, i) = training_matrix%GP(i, j)
         ENDDO
         ENDDO

         ! add noise of training data
         DO i = 1, npoints
            training_matrix%GP(i, i) = training_matrix%GP(i, i)+pao%gp_noise_var**2
         ENDDO

         ! invert covariance matrix
         ! TODO : cholesky decompose only, using a solve in predict when needed
         CALL invmat_symm(training_matrix%GP)

         ! TODO print min eigenvalue to estimate redundancy of training data and estimate gp_noise?
      ENDDO

   END SUBROUTINE pao_ml_gp_train

! **************************************************************************************************
!> \brief Uses covariance matrix to make prediction
!> \param pao ...
!> \param training_matrix ...
!> \param descriptor ...
!> \param output ...
!> \param variance ...
! **************************************************************************************************
   SUBROUTINE pao_ml_gp_predict(pao, training_matrix, descriptor, output, variance)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(training_matrix_type), INTENT(IN)             :: training_matrix
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: output
      REAL(dp), INTENT(OUT)                              :: variance

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_gp_predict', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, npoints
      REAL(dp)                                           :: cov
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: covs, weights

      npoints = SIZE(training_matrix%outputs, 1)

      ! calculate weights
      ALLOCATE (weights(npoints), covs(npoints))
      weights = 0.0_dp
      covs = 0.0_dp
      DO i = 1, npoints
         cov = kernel(pao%gp_scale, descriptor, training_matrix%inputs(i, :))
         covs(i) = cov
         weights(:) = weights(:)+cov*training_matrix%GP(i, :)
      ENDDO

      ! calculate variance
      variance = kernel(pao%gp_scale, descriptor, descriptor)
      DO i = 1, npoints
         cov = kernel(pao%gp_scale, descriptor, training_matrix%inputs(i, :))
         variance = variance-cov*weights(i)
      ENDDO

      ! build output
      output = 0.0_dp
      DO i = 1, npoints
         output(:) = output+weights(i)*training_matrix%outputs(i, :)
      ENDDO

   END SUBROUTINE pao_ml_gp_predict

! **************************************************************************************************
!> \brief Gaussian kernel used to measure covariance between two descriptors.
!> \param scale ...
!> \param descr1 ...
!> \param descr2 ...
!> \retval cov ...
! **************************************************************************************************
   FUNCTION kernel(scale, descr1, descr2) RESULT(dist)
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descr1, descr2
      REAL(dp)                                           :: dist

      REAL(dp)                                           :: fdist2

      fdist2 = SUM(((descr1-descr2)/scale)**2)
      dist = EXP(-fdist2/2.0_dp)
   END FUNCTION kernel

END MODULE pao_ml_gaussprocess
