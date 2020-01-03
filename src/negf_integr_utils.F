!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Helper functions for integration routines.
!> \par History
!>   * 06.2017 created [Sergey Chulkov]
! **************************************************************************************************
MODULE negf_integr_utils
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "./base/base_uses.f90"
#:include 'negf_integr_utils.fypp'
IMPLICIT NONE
PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_integr_utils'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: equidistant_nodes_a_b, rescale_normalised_nodes
   PUBLIC :: get_arc_radius, get_arc_smallest_angle
   PUBLIC :: rescale_nodes_arc, rescale_nodes_cos, rescale_nodes_linear, rescale_nodes_pi_phi

   INTEGER, PARAMETER, PUBLIC :: contour_shape_linear = 0, &
                                 contour_shape_arc = 1

   INTERFACE equidistant_nodes_a_b
#:for nametype1, type1 in inst_params
      MODULE PROCEDURE equidistant_${nametype1}$nodes_a_b
#:endfor
   END INTERFACE

CONTAINS

#:for nametype1, type1 in inst_params
! **************************************************************************************************
!> \brief Compute equidistant nodes on an interval [a, b], where a and b are complex numbers.
!> \param a       lower bound
!> \param b       upper bound
!> \param nnodes  number of nodes
!> \param xnodes  array to store the nodes
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE equidistant_${nametype1}$nodes_a_b(a, b, nnodes, xnodes)
      ${type1}$, INTENT(in)                              :: a, b
      INTEGER, INTENT(in)                                :: nnodes
      ${type1}$, DIMENSION(nnodes), INTENT(out)          :: xnodes

      CHARACTER(len=*), PARAMETER :: routineN = 'equidistant_${nametype1}$nodes_a_b', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i
      ${type1}$                                          :: rscale

      CPASSERT(nnodes >= 1)

      rscale = (b-a)/REAL(nnodes-1, kind=dp)
      DO i = 1, nnodes
         xnodes(i) = a+rscale*REAL(i-1, kind=dp)
      END DO
   END SUBROUTINE equidistant_${nametype1}$nodes_a_b
#:endfor

   SUBROUTINE rescale_normalised_nodes(nnodes, tnodes, a, b, shape_id, xnodes, weights)
      INTEGER, INTENT(in)                                :: nnodes
      REAL(kind=dp), DIMENSION(nnodes), INTENT(in)       :: tnodes
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      INTEGER, INTENT(in)                                :: shape_id
      COMPLEX(kind=dp), DIMENSION(nnodes), INTENT(out), &
         OPTIONAL                                        :: xnodes, weights

      CHARACTER(len=*), PARAMETER :: routineN = 'rescale_normalised_nodes', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i
      REAL(kind=dp)                                      :: rscale
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: tnodes_angle

      CALL timeset(routineN, handle)

      SELECT CASE (shape_id)
      CASE (contour_shape_linear)
         IF (PRESENT(xnodes)) &
            CALL rescale_nodes_linear(nnodes, tnodes, a, b, xnodes)

         IF (PRESENT(weights)) &
            weights(:) = b-a

      CASE (contour_shape_arc)
         ALLOCATE (tnodes_angle(nnodes))

         tnodes_angle(:) = tnodes(:)
         CALL rescale_nodes_pi_phi(a, b, nnodes, tnodes_angle)

         IF (PRESENT(xnodes)) &
            CALL rescale_nodes_arc(nnodes, tnodes_angle, a, b, xnodes)

         IF (PRESENT(weights)) THEN
            rscale = (pi-get_arc_smallest_angle(a, b))*get_arc_radius(a, b)

            DO i = 1,nnodes
               weights(i) = rscale*CMPLX(SIN(tnodes_angle(i)), -COS(tnodes_angle(i)), kind=dp)
            END DO
         END IF

         DEALLOCATE (tnodes_angle)
      CASE DEFAULT
         CPABORT("Unimplemented integration shape")
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE rescale_normalised_nodes

! **************************************************************************************************
!> \brief Compute arc radius.
!> \param a       lower bound
!> \param b       upper bound
!> \return radius
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
!> \note Assuming Re(a) < Re(b) and Im(a) < Im(b)
!             c    *
!          r   *       B-------+------
!        a  *         /   .    |
!         *        r /      .  | delta
!        *          /  phi   . |
!        A---------*-----------+------
!        <--- r --><-l->
!                  <--- r --->
! **************************************************************************************************
   PURE FUNCTION get_arc_radius(a, b) RESULT(radius)
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      REAL(kind=dp)                                      :: radius

      COMPLEX(kind=dp) :: b_minus_a

      b_minus_a = b - a

      ! l = REAL(B - A); delta = AIMAG(B - A)
      ! radius = (l^2 + delta^2) / (2 * l)
      radius = 0.5_dp * REAL(b_minus_a * CONJG(b_minus_a), kind=dp) / REAL(b_minus_a, kind=dp)
   END FUNCTION get_arc_radius

! **************************************************************************************************
!> \brief Compute the angle phi.
!> \param a       lower bound
!> \param b       upper bound
!> \return angle
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
!> \note Assuming Re(a) < Re(b) and Im(a) < Im(b)
!             c    *
!          r   *       B-------+------
!        a  *         /   .    |
!         *        r /      .  | delta
!        *          /  phi   . |
!        A---------*-----------+------
!        <--- r --><-l->
!                  <--- r --->
! **************************************************************************************************
   PURE FUNCTION get_arc_smallest_angle(a, b) RESULT(phi)
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      REAL(kind=dp)                                      :: phi

      CHARACTER(len=*), PARAMETER :: routineN = 'get_arc_smallest_angle', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp) :: b_minus_a
      REAL(kind=dp)    :: delta2, l2

      b_minus_a = b - a

      ! l = REAL(B - A); delta = AIMAG(B - A)
      ! phi = arccos((l - radius)/radius) = arccos((l^2 - delta^2) / (l^2 + delta^2))
      l2 = REAL(b_minus_a, dp)
      l2 = l2 * l2
      delta2 = AIMAG(b_minus_a)
      delta2 = delta2 * delta2

      phi = ACOS((l2 - delta2) / (l2 + delta2))
   END FUNCTION get_arc_smallest_angle

   PURE FUNCTION get_axis_rotation_angle(a, b) RESULT(phi)
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      REAL(kind=dp)                                      :: phi

      CHARACTER(len=*), PARAMETER :: routineN = 'get_axis_rotation_angle', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp) :: b_minus_a

      b_minus_a = b - a
      phi = ACOS(REAL(b_minus_a, dp) / ABS(b_minus_a))
   END FUNCTION get_axis_rotation_angle

! **************************************************************************************************
!> \brief Rescale nodes [pi, phi] -> arc[a, b] .
!> \param nnodes        number of nodes
!> \param tnodes_angle  parametrically-defined nodes to rescale 
!> \param a             lower bound
!> \param b             upper bound
!> \param xnodes        rescaled nodes (initialised on exit)
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
!> \note Assuming Re(a) < Re(b) and Im(a) < Im(b)
! **************************************************************************************************
   SUBROUTINE rescale_nodes_arc(nnodes, tnodes_angle, a, b, xnodes)
      INTEGER, INTENT(in)                                :: nnodes
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: tnodes_angle
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      COMPLEX(kind=dp), DIMENSION(:), INTENT(out)        :: xnodes

      CHARACTER(len=*), PARAMETER :: routineN = 'rescale_nodes_arc', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp)                                   :: origin
      INTEGER                                            :: i
      REAL(kind=dp)                                      :: radius

      radius = get_arc_radius(a, b)
      origin = a + CMPLX(radius, 0.0_dp, kind=dp)

      DO i = 1,nnodes
         xnodes(i) = origin + radius * CMPLX(COS(tnodes_angle(i)), SIN(tnodes_angle(i)), kind=dp)
      END DO
   END SUBROUTINE rescale_nodes_arc

! **************************************************************************************************
!> \brief Rescale nodes tnodes(i) = cos(pi/2 * (1-tnodes(i))); tnodes \in [-1 .. 1] .
!> \param tnodes parametrically-defined nodes to rescale / rescaled nodes (modified on exit)
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
!> \note Assuming Re(a) < Re(b) and Im(a) < Im(b)
! **************************************************************************************************
   SUBROUTINE rescale_nodes_cos(nnodes, tnodes)
      INTEGER, INTENT(in)                                :: nnodes
      REAL(kind=dp), DIMENSION(nnodes), INTENT(inout)    :: tnodes

      CHARACTER(len=*), PARAMETER :: routineN = 'rescale_nodes_cos', &
         routineP = moduleN//':'//routineN

      tnodes(:) = COS(0.5_dp*pi*(1.0_dp-tnodes(:)))
   END SUBROUTINE rescale_nodes_cos

! **************************************************************************************************
!> \brief Rescale nodes [-1, 1] -> [a, b] .
!> \param nnodes        number of nodes
!> \param tnodes        parametrically-defined nodes to rescale
!> \param a             lower bound
!> \param b             upper bound
!> \param xnodes        rescaled nodes (initialised on exit)
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE rescale_nodes_linear(nnodes, tnodes, a, b, xnodes)
      INTEGER, INTENT(in)                                :: nnodes
      REAL(kind=dp), DIMENSION(nnodes), INTENT(in)       :: tnodes
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      COMPLEX(kind=dp), DIMENSION(nnodes), INTENT(out)   :: xnodes

      CHARACTER(len=*), PARAMETER :: routineN = 'rescale_nodes_linear', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp)                                   :: half_len, median

      median = 0.5_dp*(b+a)
      half_len = 0.5_dp*(b-a)

      xnodes(:) = median+half_len*tnodes(:)
   END SUBROUTINE rescale_nodes_linear

! **************************************************************************************************
!> \brief Rescale nodes [-1, 1] -> [pi, phi] .
!> \param nnodes        number of nodes
!> \param a             lower bound
!> \param b             upper bound
!> \param tnodes        parametrically-defined nodes to rescale / rescaled nodes (modified on exit)
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
!> \note Assuming Re(a) < Re(b) and Im(a) < Im(b)
! **************************************************************************************************
   SUBROUTINE rescale_nodes_pi_phi(a, b, nnodes, tnodes)
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      INTEGER, INTENT(in)                                :: nnodes
      REAL(kind=dp), DIMENSION(nnodes), INTENT(inout)    :: tnodes

      CHARACTER(len=*), PARAMETER :: routineN = 'rescale_nodes_pi_phi', &
         routineP = moduleN//':'//routineN

      REAL(kind=dp)                                      :: half_pi_minus_phi, phi

      phi = get_arc_smallest_angle(a, b)
      half_pi_minus_phi = 0.5_dp*(pi-phi)

      tnodes(:) = phi+half_pi_minus_phi*(1.0_dp-tnodes(:))
   END SUBROUTINE rescale_nodes_pi_phi
END MODULE negf_integr_utils
