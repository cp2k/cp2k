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
   USE pao_types,                       ONLY: pao_env_type,&
                                              training_matrix_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_gaussprocess'

   PUBLIC ::pao_ml_gp_train, pao_ml_gp_predict, pao_ml_gp_gradient

CONTAINS

! **************************************************************************************************
!> \brief Builds the covariance matrix
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_ml_gp_train(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      INTEGER                                            :: i, ikind, info, j, npoints
      REAL(dp), DIMENSION(:), POINTER                    :: idescr, jdescr
      TYPE(training_matrix_type), POINTER                :: training_matrix

      ! TODO this could be parallelized over ranks
      DO ikind = 1, SIZE(pao%ml_training_matrices)
         training_matrix => pao%ml_training_matrices(ikind)
         npoints = SIZE(training_matrix%inputs, 2) ! number of points
         CPASSERT(SIZE(training_matrix%outputs, 2) == npoints)
         IF (npoints == 0) CYCLE ! have no training data

         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO|ML| Building covariance matrix for kind: ", &
            TRIM(training_matrix%kindname), " from ", npoints, "training points."

         ! build covariance matrix
         ALLOCATE (training_matrix%GP(npoints, npoints))
         DO i = 1, npoints
         DO j = i, npoints
            idescr => training_matrix%inputs(:, i)
            jdescr => training_matrix%inputs(:, j)
            training_matrix%GP(i, j) = kernel(pao%gp_scale, idescr, jdescr)
            training_matrix%GP(j, i) = training_matrix%GP(i, j)
         ENDDO
         ENDDO

         ! add noise of training data
         DO i = 1, npoints
            training_matrix%GP(i, i) = training_matrix%GP(i, i)+pao%gp_noise_var**2
         ENDDO

         ! compute cholesky decomposition of covariance matrix
         CALL dpotrf("U", npoints, training_matrix%GP, npoints, info)
         CPASSERT(info == 0)
      ENDDO

   END SUBROUTINE pao_ml_gp_train

! **************************************************************************************************
!> \brief Uses covariance matrix to make prediction
!> \param pao ...
!> \param ikind ...
!> \param descriptor ...
!> \param output ...
!> \param variance ...
! **************************************************************************************************
   SUBROUTINE pao_ml_gp_predict(pao, ikind, descriptor, output, variance)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: ikind
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: output
      REAL(dp), INTENT(OUT)                              :: variance

      INTEGER                                            :: i, info, npoints
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: cov, weights
      TYPE(training_matrix_type), POINTER                :: training_matrix

      training_matrix => pao%ml_training_matrices(ikind)
      npoints = SIZE(training_matrix%outputs, 2)

      ! calculate covariances between descriptor and training-points
      ALLOCATE (cov(npoints))
      DO i = 1, npoints
         cov(i) = kernel(pao%gp_scale, descriptor, training_matrix%inputs(:, i))
      ENDDO

      ! calculate weights
      ALLOCATE (weights(npoints))
      weights(:) = cov(:)
      CALL DPOTRS("U", npoints, 1, training_matrix%GP, npoints, weights, npoints, info)
      CPASSERT(info == 0)

      ! calculate predicted output
      output = 0.0_dp
      DO i = 1, npoints
         output(:) = output+weights(i)*training_matrix%outputs(:, i)
      ENDDO

      ! calculate prediction's variance
      variance = kernel(pao%gp_scale, descriptor, descriptor)-DOT_PRODUCT(weights, cov)

      IF (variance < 0.0_dp) &
         CPABORT("PAO gaussian process found negative variance")

      DEALLOCATE (cov, weights)
   END SUBROUTINE pao_ml_gp_predict

! **************************************************************************************************
!> \brief Calculate gradient of Gaussian process
!> \param pao ...
!> \param ikind ...
!> \param descriptor ...
!> \param outer_deriv ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE pao_ml_gp_gradient(pao, ikind, descriptor, outer_deriv, gradient)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: ikind
      REAL(dp), DIMENSION(:), INTENT(IN), TARGET         :: descriptor
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: outer_deriv
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: gradient

      INTEGER                                            :: i, info, npoints
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: cov_deriv, weights_deriv
      REAL(dp), DIMENSION(SIZE(descriptor))              :: kg
      TYPE(training_matrix_type), POINTER                :: training_matrix

      training_matrix => pao%ml_training_matrices(ikind)
      npoints = SIZE(training_matrix%outputs, 2)

      ! calculate derivative of weights
      ALLOCATE (weights_deriv(npoints))
      DO i = 1, npoints
         weights_deriv(i) = SUM(outer_deriv*training_matrix%outputs(:, i))
      ENDDO

      ! calculate derivative of covariances
      ALLOCATE (cov_deriv(npoints))
      cov_deriv(:) = weights_deriv(:)
      CALL DPOTRS("U", npoints, 1, training_matrix%GP, npoints, cov_deriv, npoints, info)
      CPASSERT(info == 0)

      ! calculate derivative of kernel
      gradient(:) = 0.0_dp
      DO i = 1, npoints
         kg = kernel_grad(pao%gp_scale, descriptor, training_matrix%inputs(:, i))
         gradient(:) = gradient(:)+kg(:)*cov_deriv(i)
      ENDDO

      DEALLOCATE (cov_deriv, weights_deriv)
   END SUBROUTINE pao_ml_gp_gradient

! **************************************************************************************************
!> \brief Gaussian kernel used to measure covariance between two descriptors.
!> \param scale ...
!> \param descr1 ...
!> \param descr2 ...
!> \retval cov ...
! **************************************************************************************************
   PURE FUNCTION kernel(scale, descr1, descr2) RESULT(cov)
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descr1, descr2
      REAL(dp)                                           :: cov

      REAL(dp)                                           :: fdist2
      REAL(dp), DIMENSION(SIZE(descr1))                  :: diff

      diff = descr1-descr2
      fdist2 = SUM((diff/scale)**2)
      cov = EXP(-fdist2/2.0_dp)
   END FUNCTION kernel

! **************************************************************************************************
!> \brief Gradient of Gaussian kernel wrt descr1
!> \param scale ...
!> \param descr1 ...
!> \param descr2 ...
!> \retval cov ...
! **************************************************************************************************
   PURE FUNCTION kernel_grad(scale, descr1, descr2) RESULT(grad)
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descr1, descr2
      REAL(dp), DIMENSION(SIZE(descr1))                  :: grad

      REAL(dp)                                           :: cov, fdist2
      REAL(dp), DIMENSION(SIZE(descr1))                  :: diff

      diff = descr1-descr2
      fdist2 = SUM((diff/scale)**2)
      cov = EXP(-fdist2/2.0_dp)
      grad(:) = cov*(-diff/scale**2)

   END FUNCTION kernel_grad

END MODULE pao_ml_gaussprocess
