!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Neural Network implementation
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_neuralnet
   USE kinds,                           ONLY: dp
   USE pao_types,                       ONLY: pao_env_type,&
                                              training_matrix_type
   USE parallel_rng_types,              ONLY: create_rng_stream,&
                                              delete_rng_stream,&
                                              next_random_number,&
                                              rng_stream_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_neuralnet'

   PUBLIC ::pao_ml_nn_train, pao_ml_nn_predict, pao_ml_nn_gradient

   ! TODO turn these into input parameters
   REAL(dp), PARAMETER   :: step_size = 0.001_dp
   INTEGER, PARAMETER    :: nlayers = 3
   REAL(dp), PARAMETER   :: convergence_eps = 1e-7_dp
   INTEGER, PARAMETER    :: max_training_cycles = 50000

CONTAINS

! **************************************************************************************************
!> \brief Uses neural network to make a prediction
!> \param pao ...
!> \param ikind ...
!> \param descriptor ...
!> \param output ...
!> \param variance ...
! **************************************************************************************************
   SUBROUTINE pao_ml_nn_predict(pao, ikind, descriptor, output, variance)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: ikind
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: output
      REAL(dp), INTENT(OUT)                              :: variance

      TYPE(training_matrix_type), POINTER                :: training_matrix

      training_matrix => pao%ml_training_matrices(ikind)

      CALL nn_eval(training_matrix%NN, input=descriptor, prediction=output)

      variance = 0.0_dp ! Neural Networks don't provide a variance
   END SUBROUTINE pao_ml_nn_predict

! **************************************************************************************************
!> \brief Calculate gradient of neural network
!> \param pao ...
!> \param ikind ...
!> \param descriptor ...
!> \param outer_deriv ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE pao_ml_nn_gradient(pao, ikind, descriptor, outer_deriv, gradient)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: ikind
      REAL(dp), DIMENSION(:), INTENT(IN), TARGET         :: descriptor
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: outer_deriv
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: gradient

      INTEGER                                            :: i, ilayer, j, nlayers, width, width_in, &
                                                            width_out
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: backward, forward
      REAL(dp), DIMENSION(:, :, :), POINTER              :: A

      A => pao%ml_training_matrices(ikind)%NN

      nlayers = SIZE(A, 1)
      width = SIZE(A, 2); CPASSERT(SIZE(A, 2) == SIZE(A, 3))
      width_in = SIZE(descriptor)
      width_out = SIZE(outer_deriv)

      ALLOCATE (forward(0:nlayers, width), backward(0:nlayers, width))

      forward = 0.0_dp
      forward(0, 1:width_in) = descriptor

      DO ilayer = 1, nlayers
      DO i = 1, width
      DO j = 1, width
         forward(ilayer, i) = forward(ilayer, i)+A(ilayer, i, j)*TANH(forward(ilayer-1, j))
      ENDDO
      ENDDO
      ENDDO

      ! Turning Point ------------------------------------------------------------------------------
      backward = 0.0_dp
      backward(nlayers, 1:width_out) = outer_deriv(:)

      DO ilayer = nlayers, 1, -1
      DO i = 1, width
      DO j = 1, width
         backward(ilayer-1, j) = backward(ilayer-1, j)+backward(ilayer, i)*A(ilayer, i, j)*(1.0_dp-TANH(forward(ilayer-1, j))**2)
      ENDDO
      ENDDO
      ENDDO

      gradient(:) = backward(0, 1:width_in)

      DEALLOCATE (forward, backward)
   END SUBROUTINE pao_ml_nn_gradient

! **************************************************************************************************
!> \brief Trains the neural network on given training points
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_ml_nn_train(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      INTEGER                                            :: i, icycle, ikind, ilayer, ipoint, j, &
                                                            npoints, width, width_in, width_out
      REAL(dp)                                           :: bak, eps, error, error1, error2, num_grad
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: prediction
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: gradient
      TYPE(rng_stream_type), POINTER                     :: rng_stream
      TYPE(training_matrix_type), POINTER                :: training_matrix

      NULLIFY (rng_stream)

      ! TODO this could be parallelized over ranks
      DO ikind = 1, SIZE(pao%ml_training_matrices)
         training_matrix => pao%ml_training_matrices(ikind)

         npoints = SIZE(training_matrix%inputs, 2) ! number of points
         CPASSERT(SIZE(training_matrix%outputs, 2) == npoints)
         IF (npoints == 0) CYCLE

         !TODO proper output
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO|ML| Training neural network for kind: ", &
            TRIM(training_matrix%kindname), " from ", npoints, "training points."

         ! determine network width and allocate it
         width_in = SIZE(training_matrix%inputs, 1)
         width_out = SIZE(training_matrix%outputs, 1)
         width = MAX(width_in, width_out)
         ALLOCATE (training_matrix%NN(nlayers, width, width))

         ! initialize network with random numbers from -1.0 ... +1.0
         CALL create_rng_stream(rng_stream, name="pao_nn")
         DO ilayer = 1, nlayers
         DO i = 1, width
         DO j = 1, width
            training_matrix%NN(ilayer, i, j) = -1.0_dp+2.0_dp*next_random_number(rng_stream)
         ENDDO
         ENDDO
         ENDDO
         CALL delete_rng_stream(rng_stream)

         ! train the network using backpropagation
         ALLOCATE (gradient(nlayers, width, width))
         DO icycle = 1, max_training_cycles
            error = 0.0_dp
            gradient = 0.0_dp
            DO ipoint = 1, npoints
               CALL nn_backpropagate(training_matrix%NN, &
                                     input=training_matrix%inputs(:, ipoint), &
                                     goal=training_matrix%outputs(:, ipoint), &
                                     gradient=gradient, &
                                     error=error)
            ENDDO
            training_matrix%NN(:, :, :) = training_matrix%NN-step_size*gradient

            IF (pao%iw > 0 .AND. MOD(icycle, 100) == 0) WRITE (pao%iw, *) &
               "PAO|ML| ", TRIM(training_matrix%kindname), &
               " training-cycle:", icycle, "SQRT(error):", SQRT(error), "grad:", SUM(gradient**2)

            IF (SUM(gradient**2) < convergence_eps) EXIT
         ENDDO

         ! numeric gradient for debugging ----------------------------------------------------------
         IF (.FALSE.) THEN
            eps = 1e-4_dp
            ilayer = 1
            ipoint = 1
            error = 0.0_dp
            gradient = 0.0_dp
            CALL nn_backpropagate(training_matrix%NN, &
                                  input=training_matrix%inputs(:, ipoint), &
                                  goal=training_matrix%outputs(:, ipoint), &
                                  gradient=gradient, &
                                  error=error)

            ALLOCATE (prediction(width_out))
            DO i = 1, width
            DO j = 1, width
               bak = training_matrix%NN(ilayer, i, j)

               training_matrix%NN(ilayer, i, j) = bak+eps
               CALL nn_eval(training_matrix%NN, &
                            input=training_matrix%inputs(:, ipoint), &
                            prediction=prediction)
               error1 = SUM((training_matrix%outputs(:, ipoint)-prediction)**2)

               training_matrix%NN(ilayer, i, j) = bak-eps
               CALL nn_eval(training_matrix%NN, &
                            input=training_matrix%inputs(:, ipoint), &
                            prediction=prediction)
               error2 = SUM((training_matrix%outputs(:, ipoint)-prediction)**2)

               training_matrix%NN(ilayer, i, j) = bak
               num_grad = (error1-error2)/(2.0_dp*eps)
               IF (pao%iw > 0) WRITE (pao%iw, *) "PAO|ML| Numeric gradient:", i, j, gradient(ilayer, i, j), num_grad

            ENDDO
            ENDDO
            DEALLOCATE (prediction)
         ENDIF
         !------------------------------------------------------------------------------------------

         DEALLOCATE (gradient)

         ! test training points individually
         ALLOCATE (prediction(width_out))
         DO ipoint = 1, npoints
            CALL nn_eval(training_matrix%NN, &
                         input=training_matrix%inputs(:, ipoint), &
                         prediction=prediction)
            error = MAXVAL(ABS(training_matrix%outputs(:, ipoint)-prediction))
            IF (pao%iw > 0) WRITE (pao%iw, *) "PAO|ML| ", TRIM(training_matrix%kindname), &
               " verify training-point:", ipoint, "SQRT(error):", SQRT(error)
         ENDDO
         DEALLOCATE (prediction)

      ENDDO

   END SUBROUTINE pao_ml_nn_train

! **************************************************************************************************
!> \brief Evaluates the neural network for a given input
!> \param A ...
!> \param input ...
!> \param prediction ...
! **************************************************************************************************
   SUBROUTINE nn_eval(A, input, prediction)
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: A
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: input
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: prediction

      INTEGER                                            :: i, ilayer, j, nlayers, width, width_in, &
                                                            width_out
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: forward

      nlayers = SIZE(A, 1)
      width = SIZE(A, 2); CPASSERT(SIZE(A, 2) == SIZE(A, 3))
      width_in = SIZE(input)
      width_out = SIZE(prediction)

      ALLOCATE (forward(0:nlayers, width))

      forward = 0.0_dp
      forward(0, 1:width_in) = input(:)

      DO ilayer = 1, nlayers
      DO i = 1, width
      DO j = 1, width
         forward(ilayer, i) = forward(ilayer, i)+A(ilayer, i, j)*TANH(forward(ilayer-1, j))
      ENDDO
      ENDDO
      ENDDO

      prediction(:) = forward(nlayers, 1:width_out)

   END SUBROUTINE nn_eval

! **************************************************************************************************
!> \brief Uses backpropagation to calculate the gradient for a given training point
!> \param A ...
!> \param input ...
!> \param goal ...
!> \param error ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE nn_backpropagate(A, input, goal, error, gradient)
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: A
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: input, goal
      REAL(dp), INTENT(INOUT)                            :: error
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: gradient

      INTEGER                                            :: i, ilayer, j, nlayers, width, width_in, &
                                                            width_out
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: prediction
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: backward, forward

      nlayers = SIZE(A, 1)
      width = SIZE(A, 2); CPASSERT(SIZE(A, 2) == SIZE(A, 3))
      width_in = SIZE(input)
      width_out = SIZE(goal)

      ALLOCATE (forward(0:nlayers, width), prediction(width_out), backward(0:nlayers, width))

      forward = 0.0_dp
      forward(0, 1:width_in) = input

      DO ilayer = 1, nlayers
      DO i = 1, width
      DO j = 1, width
         forward(ilayer, i) = forward(ilayer, i)+A(ilayer, i, j)*TANH(forward(ilayer-1, j))
      ENDDO
      ENDDO
      ENDDO

      prediction(:) = forward(nlayers, 1:width_out)

      error = error+SUM((prediction-goal)**2)

      ! Turning Point ------------------------------------------------------------------------------
      backward = 0.0_dp
      backward(nlayers, 1:width_out) = prediction-goal

      DO ilayer = nlayers, 1, -1
      DO i = 1, width
      DO j = 1, width
         gradient(ilayer, i, j) = gradient(ilayer, i, j)+2.0_dp*backward(ilayer, i)*TANH(forward(ilayer-1, j))
         backward(ilayer-1, j) = backward(ilayer-1, j)+backward(ilayer, i)*A(ilayer, i, j)*(1.0_dp-TANH(forward(ilayer-1, j))**2)
      ENDDO
      ENDDO
      ENDDO

      DEALLOCATE (forward, backward, prediction)
   END SUBROUTINE nn_backpropagate

END MODULE pao_ml_neuralnet
