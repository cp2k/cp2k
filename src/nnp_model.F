!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with core routines for artificial neural networks
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_model

   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              max_line_length
   USE nnp_environment_types,           ONLY: &
        nnp_actfnct_cos, nnp_actfnct_exp, nnp_actfnct_gaus, nnp_actfnct_invsig, nnp_actfnct_lin, &
        nnp_actfnct_quad, nnp_actfnct_sig, nnp_actfnct_softplus, nnp_actfnct_tanh, nnp_arc_type, &
        nnp_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_model'

   PUBLIC :: nnp_write_arc, &
             nnp_predict, &
             nnp_gradients

CONTAINS

! ***************************************************************************
!> \brief Write neural network architecture information
!> \param nnp ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE nnp_write_arc(nnp, para_env)
      TYPE(nnp_type), POINTER                            :: nnp
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_write_arc', routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: my_label
      INTEGER                                            :: i, j, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      my_label = "NNP| "
      IF (para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         DO i = 1, nnp%n_ele
            WRITE (unit_nr, *) TRIM(my_label)//" Neural network specification for element "// &
               nnp%ele(i)//":"
            DO j = 1, nnp%n_layer
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A,1X,I2)') TRIM(my_label), &
                  nnp%arc(i)%n_nodes(j), "nodes in layer", j
            END DO
         END DO
      END IF

      RETURN

   END SUBROUTINE nnp_write_arc

! ***************************************************************************
!> \brief Predict energy by evaluating neural network
!> \param arc ...
!> \param nnp ...
!> \param i_com ...
! **************************************************************************************************
   SUBROUTINE nnp_predict(arc, nnp, i_com)
      TYPE(nnp_arc_type), INTENT(INOUT)                  :: arc
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: i_com

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_predict', routineP = moduleN//':'//routineN

      CHARACTER(len=max_line_length)                     :: msg_str
      INTEGER                                            :: handle, i, j
      REAL(KIND=dp)                                      :: norm

      CALL timeset(routineN, handle)

      DO i = 2, nnp%n_layer
         ! Calculate node(i)
         arc%layer(i)%node(:) = 0.0_dp
         !Perform matrix-vector product
         !y := alpha*A*x + beta*y
         !with A = layer(i)*weights
         !and  x = layer(i-1)%node
         !and  y = layer(i)%node
         CALL DGEMV('T', & !transpose matrix A
                    arc%n_nodes(i - 1), & !number of rows of A
                    arc%n_nodes(i), & !number of columns of A
                    1.0_dp, & ! alpha
                    arc%layer(i)%weights(:, :, i_com), & !matrix A
                    arc%n_nodes(i - 1), & !leading dimension of A
                    arc%layer(i - 1)%node, & !vector x
                    1, & !increment for the elements of x
                    1.0_dp, & !beta
                    arc%layer(i)%node, & !vector y
                    1) !increment for the elements of y

         ! Add bias weight
         DO j = 1, arc%n_nodes(i)
            arc%layer(i)%node(j) = arc%layer(i)%node(j) + arc%layer(i)%bweights(j, i_com)
         END DO

         ! Normalize by number of nodes in previous layer if requested
         IF (nnp%normnodes) THEN
            norm = 1.0_dp/REAL(arc%n_nodes(i - 1), dp)
            DO j = 1, arc%n_nodes(i)
               arc%layer(i)%node(j) = arc%layer(i)%node(j)*norm
            END DO
         END IF

         ! Store node values before application of activation function
         ! (needed for derivatives)
         DO j = 1, arc%n_nodes(i)
            arc%layer(i)%node_grad(j) = arc%layer(i)%node(j)
         END DO

         ! Apply activation function:
         SELECT CASE (nnp%actfnct(i - 1))
         CASE (nnp_actfnct_tanh)
            arc%layer(i)%node(:) = TANH(arc%layer(i)%node(:))
         CASE (nnp_actfnct_gaus)
            arc%layer(i)%node(:) = EXP(-0.5_dp*arc%layer(i)%node(:)**2)
         CASE (nnp_actfnct_lin)
            CONTINUE
         CASE (nnp_actfnct_cos)
            arc%layer(i)%node(:) = COS(arc%layer(i)%node(:))
         CASE (nnp_actfnct_sig)
            arc%layer(i)%node(:) = 1.0_dp/(1.0_dp + EXP(-1.0_dp*arc%layer(i)%node(:)))
         CASE (nnp_actfnct_invsig)
            arc%layer(i)%node(:) = 1.0_dp - 1.0_dp/(1.0_dp + EXP(-1.0_dp*arc%layer(i)%node(:)))
         CASE (nnp_actfnct_exp)
            arc%layer(i)%node(:) = EXP(-1.0_dp*arc%layer(i)%node(:))
         CASE (nnp_actfnct_softplus)
            arc%layer(i)%node(:) = LOG(EXP(arc%layer(i)%node(:)) + 1.0_dp)
         CASE (nnp_actfnct_quad)
            arc%layer(i)%node(:) = arc%layer(i)%node(:)**2
         CASE DEFAULT
            msg_str = "NNP| Error: Unknown activation function"
            CPABORT(msg_str)
         END SELECT
      END DO

      CALL timestop(handle)

      RETURN

   END SUBROUTINE nnp_predict

! ***************************************************************************
!> \brief Calculate gradients of neural network
!> \param arc ...
!> \param nnp ...
!> \param i_com ...
!> \param denergydsym ...
! **************************************************************************************************
   SUBROUTINE nnp_gradients(arc, nnp, i_com, denergydsym)
      TYPE(nnp_arc_type)                                 :: arc
      TYPE(nnp_type)                                     :: nnp
      INTEGER, INTENT(IN)                                :: i_com
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: denergydsym

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_gradients', routineP = moduleN//':'//routineN

      CHARACTER(len=max_line_length)                     :: msg_str
      INTEGER                                            :: handle, i, j, k
      REAL(KIND=dp)                                      :: norm

      CALL timeset(routineN, handle)

      norm = 1.0_dp

      DO i = 2, nnp%n_layer

         ! Apply activation function:
         SELECT CASE (nnp%actfnct(i - 1))
         CASE (nnp_actfnct_tanh)
            arc%layer(i)%node_grad(:) = 1.0_dp - arc%layer(i)%node(:)**2 !tanh(x)'=1-tanh(x)**2
         CASE (nnp_actfnct_gaus)
            arc%layer(i)%node_grad(j) = -1.0_dp*arc%layer(i)%node(j)*arc%layer(i)%node_grad(j)
         CASE (nnp_actfnct_lin)
            arc%layer(i)%node_grad(:) = 1.0_dp
         CASE (nnp_actfnct_cos)
            arc%layer(i)%node_grad(j) = -SIN(arc%layer(i)%node_grad(j))
         CASE (nnp_actfnct_sig)
            arc%layer(i)%node_grad(j) = EXP(-arc%layer(i)%node_grad(j))/ &
                                        (1.0_dp + EXP(-1.0_dp*arc%layer(i)%node_grad(j)))**2
         CASE (nnp_actfnct_invsig)
            arc%layer(i)%node_grad(j) = -1.0_dp*EXP(-1.0_dp*arc%layer(i)%node_grad(j))/ &
                                        (1.0_dp + EXP(-1.0_dp*arc%layer(i)%node_grad(j)))**2
         CASE (nnp_actfnct_exp)
            arc%layer(i)%node_grad(j) = -1.0_dp*arc%layer(i)%node(j)
         CASE (nnp_actfnct_softplus)
            arc%layer(i)%node_grad(j) = (EXP(arc%layer(i)%node(j)) + 1.0_dp)/ &
                                        EXP(arc%layer(i)%node(j))
         CASE (nnp_actfnct_quad)
            arc%layer(i)%node_grad(j) = 2.0_dp*arc%layer(i)%node_grad(j)
         CASE DEFAULT
            msg_str = "NNP| Error: Unknown activation function"
            CPABORT(msg_str)
         END SELECT
         ! Normalize by number of nodes in previous layer if requested
         IF (nnp%normnodes) THEN
            norm = 1.0_dp/REAL(arc%n_nodes(i - 1), dp)
            arc%layer(i)%node_grad(:) = norm*arc%layer(i)%node_grad(:)
         END IF

      END DO

      ! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
      DO j = 1, arc%n_nodes(2)
         DO i = 1, arc%n_nodes(1)
            arc%layer(2)%tmp_der(i, j) = arc%layer(2)%node_grad(j)*arc%layer(2)%weights(i, j, i_com)
         END DO
      END DO

      DO k = 3, nnp%n_layer
         ! Reset tmp_der:
         arc%layer(k)%tmp_der(:, :) = 0.0_dp
         !Perform matrix-matrix product
         !C := alpha*A*B + beta*C
         !with A = layer(k-1)%tmp_der
         !and  B = layer(k)%weights
         !and  C = tmp
         CALL DGEMM('N', & !don't transpose matrix A
                    'N', & !don't transpose matrix B
                    arc%n_nodes(1), & !number of rows of A
                    arc%n_nodes(k), & !number of columns of B
                    arc%n_nodes(k - 1), & !number of col of A and nb of rows of B
                    1.0_dp, & !alpha
                    arc%layer(k - 1)%tmp_der, & !matrix A
                    arc%n_nodes(1), & !leading dimension of A
                    arc%layer(k)%weights(:, :, i_com), & !matrix B
                    arc%n_nodes(k - 1), & !leading dimension of B
                    1.0_dp, & !beta
                    arc%layer(k)%tmp_der, & !matrix C
                    arc%n_nodes(1)) !leading dimension of C

         ! sum over all nodes in the target layer
         DO j = 1, arc%n_nodes(k)
            ! sum over input layer
            DO i = 1, arc%n_nodes(1)
               arc%layer(k)%tmp_der(i, j) = arc%layer(k)%node_grad(j)* &
                                            arc%layer(k)%tmp_der(i, j)
            END DO
         END DO
      END DO

      DO i = 1, arc%n_nodes(1)
         denergydsym(i) = arc%layer(nnp%n_layer)%tmp_der(i, 1)
      END DO

      CALL timestop(handle)

      RETURN

   END SUBROUTINE nnp_gradients

END MODULE nnp_model
