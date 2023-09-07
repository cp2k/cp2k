!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      If parallel mode is distributed certain combination of
!>      "in_use" and "in_space" can not be used.
!>      For performance reasons it would be better to have the loops
!>      over g-vectros in the gather/scatter routines in new subprograms
!>      with the actual arrays (also the addressing) in the parameter list
!> \par History
!>      JGH (29-Dec-2000) : Changes for parallel use
!>      JGH (13-Mar-2001) : added timing calls
!>      JGH (26-Feb-2003) : OpenMP enabled
!>      JGH (17-Nov-2007) : Removed mass arrays
!>      JGH (01-Dec-2007) : Removed and renamed routines
!>      JGH (04-Jul-2019) : added pw_multiply routine
!>      03.2008 [tlaino] : Splitting pw_types into pw_types and pw_methods
!>      08.2023 (F.Stein) : Major refactoring
!> \author apsi
! **************************************************************************************************
MODULE pw_methods

   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit, &
                              cp_to_string
   USE fft_tools, ONLY: BWFFT, &
                        FWFFT, &
                        fft3d
   USE kahan_sum, ONLY: accurate_sum
   USE kinds, ONLY: dp
   USE machine, ONLY: m_memory
   USE mathconstants, ONLY: z_zero
   USE pw_copy_all, ONLY: pw_copy_match
   USE pw_fpga, ONLY: pw_fpga_c1dr3d_3d_dp, &
                      pw_fpga_c1dr3d_3d_sp, &
                      pw_fpga_init_bitstream, &
                      pw_fpga_r3dc1d_3d_dp, &
                      pw_fpga_r3dc1d_3d_sp
   USE pw_gpu, ONLY: pw_gpu_c1dr3d_3d, &
                     pw_gpu_c1dr3d_3d_ps, &
                     pw_gpu_r3dc1d_3d, &
                     pw_gpu_r3dc1d_3d_ps
   USE pw_grid_types, ONLY: HALFSPACE, &
                            PW_MODE_DISTRIBUTED, &
                            PW_MODE_LOCAL, &
                            pw_grid_type
   USE pw_types, ONLY: pw_r1d_type, pw_r3d_type, pw_c1d_type, pw_c3d_type, REALSPACE, RECIPROCALSPACE, pw_create, pw_release
#include "../base/base_uses.f90"

   IMPLICIT NONE

   #:include 'pw_types.fypp'

   PRIVATE

   PUBLIC :: pw_zero, pw_structure_factor, pw_smoothing
   PUBLIC :: pw_copy, pw_copy_to_array, pw_copy_from_array, pw_axpy, pw_transfer, pw_scale
   PUBLIC :: pw_gauss_damp, pw_compl_gauss_damp, pw_derive, pw_laplace, pw_dr2, pw_write, pw_multiply
   PUBLIC :: pw_log_deriv_gauss, pw_log_deriv_compl_gauss, pw_log_deriv_mix_cl, pw_log_deriv_trunc
   PUBLIC :: pw_gauss_damp_mix, pw_multiply_by
   PUBLIC :: pw_divide, pw_divide_by
   PUBLIC :: pw_integral_ab, pw_integral_a2b
   PUBLIC :: pw_dr2_gg, pw_integrate_function
   PUBLIC :: pw_set, pw_truncated

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_methods'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   INTEGER, PARAMETER, PUBLIC  ::  do_accurate_sum = 0, &
                                  do_standard_sum = 1

   INTERFACE pw_zero
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_zero_${kind}$
      #:endfor
   END INTERFACE pw_zero

   INTERFACE pw_copy
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            #:if kind[1:] == kind2[1:]
               MODULE PROCEDURE pw_copy_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endfor
   END INTERFACE pw_copy

   INTERFACE pw_copy_to_array
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            #:if kind[1:] == kind2[1:]
               MODULE PROCEDURE pw_copy_to_array_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endfor
   END INTERFACE pw_copy_to_array

   INTERFACE pw_copy_from_array
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            #:if kind[1:] == kind2[1:]
               MODULE PROCEDURE pw_copy_from_array_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endfor
   END INTERFACE pw_copy_from_array

   INTERFACE pw_scale
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_scale_${kind}$
      #:endfor
   END INTERFACE pw_scale

   INTERFACE pw_axpy
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_axpy_${kind}$
      #:endfor
   END INTERFACE pw_axpy

   INTERFACE pw_write
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_write_${kind}$
      #:endfor
   END INTERFACE pw_write

   INTERFACE pw_multiply
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            #:if kind[1:]==kind2[1:]
               MODULE PROCEDURE pw_multiply_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endfor
   END INTERFACE pw_multiply

   INTERFACE pw_divide
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_divide_${kind}$
      #:endfor
   END INTERFACE pw_divide

   INTERFACE pw_divide_by
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_divide_by_${kind}$
      #:endfor
   END INTERFACE pw_divide_by

   INTERFACE pw_integral_ab
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            #:if kind[1:] == kind2[1:]
               MODULE PROCEDURE pw_integral_ab_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endfor
   END INTERFACE pw_integral_ab

   INTERFACE pw_integral_a2b
      #:for kind in pw_kinds
         #:if kind.endswith('1d')
            MODULE PROCEDURE pw_integral_a2b_${kind}$
         #:endif
      #:endfor
   END INTERFACE pw_integral_a2b

   INTERFACE pw_multiply_by
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            #:if kind[1:]==kind2[1:]
               MODULE PROCEDURE pw_multiply_by_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endfor
   END INTERFACE pw_multiply_by

   INTERFACE pw_transfer
      #:for kind in pw_kinds
         #:for kind2 in pw_kinds
            MODULE PROCEDURE pw_transfer_${kind}$_${kind2}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE fft_wrap_pw1pw2
      #:for kind in pw_kinds
         #:if kind != 'r1d'
            #:for kind2 in pw_kinds
               #:if kind2 != 'r1d' and (kind[0] != 'r' or kind2[0] != 'r')
                  MODULE PROCEDURE fft_wrap_pw1pw2_${kind}$_${kind2}$
               #:endif
            #:endfor
         #:endif
      #:endfor
   END INTERFACE fft_wrap_pw1pw2

   INTERFACE pw_gather
      MODULE PROCEDURE pw_gather_s_c1d, pw_gather_p_c1d
   END INTERFACE

   INTERFACE pw_scatter
      MODULE PROCEDURE pw_scatter_s_c1d, pw_scatter_p_c1d
   END INTERFACE

   INTERFACE pw_set
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_set_value_${kind}$, pw_set_array_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_integrate_function
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_integrate_function_${kind}$
      #:endfor
   END INTERFACE

CONTAINS

   #:for kind, type in zip(pw_kinds, pw_types)
! **************************************************************************************************
!> \brief Set values of a pw type to zero
!> \param pw ...
!> \par History
!>      none
!> \author apsi
! **************************************************************************************************
      SUBROUTINE pw_zero_${kind}$ (pw)
         TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_zero_${kind}$'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
         pw%array = 0.0_dp
!$OMP END PARALLEL WORKSHARE

         CALL timestop(handle)

      END SUBROUTINE pw_zero_${kind}$

      #:for kind2, type2 in zip(pw_kinds, pw_types)
         #:if kind[1] == kind2[1]

! **************************************************************************************************
!> \brief copy a pw type variable
!> \param pw1 ...
!> \param pw2 ...
! **************************************************************************************************
            SUBROUTINE pw_copy_${kind}$_${kind2}$ (pw1, pw2)

               TYPE(pw_${kind}$_type), INTENT(IN)                 :: pw1
               TYPE(pw_${kind2}$_type), INTENT(INOUT)              :: pw2

               CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_copy_${kind}$_${kind2}$'

               INTEGER                                            :: handle
               #:if kind.endswith('1d') and kind2.endswith('1d')
                  INTEGER :: i, j, ng, ng1, ng2
               #:endif

               CALL timeset(routineN, handle)

               IF (pw1%in_space /= pw2%in_space) &
                  CPABORT("Copy requires both PW grids to be in the same space!")

               IF (ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
                  pw2%array = ${type2type('pw1%array', type, type2)}$
!$OMP END PARALLEL WORKSHARE
               ELSE IF (pw1%pw_grid%spherical .AND. pw2%pw_grid%spherical) THEN
                  IF (.NOT. pw_compatible(pw1%pw_grid, pw2%pw_grid)) &
                     CPABORT("Grids must be compatible if they are spherical!")

                  #:if kind.endswith('1d') and kind2.endswith('1d')
                     ng1 = SIZE(pw1%array)
                     ng2 = SIZE(pw2%array)
                     ng = MIN(ng1, ng2)
!$OMP PARALLEL DO PRIVATE(i) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                     DO i = 1, ng
                        pw2%array(i) = ${type2type('pw1%array(i)', type, type2)}$
                     END DO
!$OMP END PARALLEL DO
                     IF (ng2 > ng) THEN
!$OMP PARALLEL DO PRIVATE(i) DEFAULT(NONE) SHARED(ng, ng2, pw2)
                        DO i = ng + 1, ng2
                           pw2%array(i) = ${type2type('0.0_dp', type, type2)}$
                        END DO
!$OMP END PARALLEL DO
                     END IF
                  #:else
                     CPABORT("Spherical grids are only implemented for 1D grids!")
                  #:endif

               ELSE IF (.NOT. (pw1%pw_grid%spherical .OR. pw2%pw_grid%spherical)) THEN
                  #:if kind.endswith('1d') and kind2.endswith('1d')
                     ng1 = SIZE(pw1%array)
                     ng2 = SIZE(pw2%array)
                     ng = MIN(ng1, ng2)
                     IF (ng2 > ng1) CALL pw_zero(pw2)

                     IF (pw1%pw_grid%id_nr == pw2%pw_grid%reference) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng,pw1,pw2)
                           DO i = 1, ng
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(i) = ${type2type('pw1%array(j)', type, type2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng,pw1,pw2)
                           DO i = 1, ng
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(j) = ${type2type('pw1%array(i)', type, type2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE IF (pw2%pw_grid%id_nr == pw1%pw_grid%reference) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng,pw1,pw2)
                           DO i = 1, ng
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(i) = ${type2type('pw1%array(j)', type, type2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng,pw1,pw2)
                           DO i = 1, ng
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(j) = ${type2type('pw1%array(i)', type, type2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE
                        #:if kind == 'c1d' and kind2 == 'c1d'
                           CALL pw_copy_match(pw1, pw2)
                        #:else
                           CALL cp_abort(__LOCATION__, "Copying of PW grids of noncompatible grids "// &
                                         "not implemented for ${kind}$/${kind2}$!")
                        #:endif
                     END IF
                     #! 1D/3D grids
                  #:elif kind[1:]=='3d' and kind[1:]=='3d'
                     IF (ALL(SHAPE(pw1%array) == SHAPE(pw2%array))) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
                        pw2%array = ${type2type('pw1%array', type, type2)}$
!$OMP END PARALLEL WORKSHARE
                     ELSE
                        CPABORT("Arrays must have the same shape for a copy!")
                     END IF
                  #:else
                     CPABORT("Copy not available for mixed 1D/3D grids!")
                  #:endif
               ELSE
                  CPABORT("Incompatible grids")
               END IF

               CALL timestop(handle)

            END SUBROUTINE pw_copy_${kind}$_${kind2}$

! **************************************************************************************************
!> \brief copy a pw type into an array
!> \param pw1 ...
!> \param array ...
! **************************************************************************************************
            SUBROUTINE pw_copy_to_array_${kind}$_${kind2}$ (pw1, array)

               TYPE(pw_${kind}$_type), INTENT(IN)                 :: pw1
               ${type2}$, INTENT(INOUT)              :: array

               CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_copy_to_array_${kind}$_${kind2}$'

               INTEGER                                            :: handle

               CALL timeset(routineN, handle)

               IF (ANY(SHAPE(pw1%array) /= SHAPE(array))) &
                  CPABORT("Copy requires PW grid and array to have the same shape!")

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,array)
               array = ${type2type("pw1%array", type, type2)}$
!$OMP END PARALLEL WORKSHARE

               CALL timestop(handle)

            END SUBROUTINE pw_copy_to_array_${kind}$_${kind2}$

! **************************************************************************************************
!> \brief copy a pw type into an array
!> \param pw1 ...
!> \param array ...
! **************************************************************************************************
            SUBROUTINE pw_copy_from_array_${kind}$_${kind2}$ (pw1, array)

               TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw1
               ${type2}$, INTENT(IN)                              :: array

               CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_copy_from_array_${kind}$_${kind2}$'

               INTEGER                                            :: handle

               CALL timeset(routineN, handle)

               IF (ANY(SHAPE(pw1%array) /= SHAPE(array))) &
                  CPABORT("Copy requires PW grid and array to have the same shape!")

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,array)
               pw1%array = ${type2type("array", type2, type)}$
!$OMP END PARALLEL WORKSHARE

               CALL timestop(handle)

            END SUBROUTINE pw_copy_from_array_${kind}$_${kind2}$

         #:endif
      #:endfor

! **************************************************************************************************
!> \brief multiplies pw coeffs with a number
!> \param pw ...
!> \param a ...
! **************************************************************************************************
      SUBROUTINE pw_scale_${kind}$ (pw, a)

         TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw
         REAL(KIND=dp), INTENT(IN)                          :: a

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_scale_${kind}$'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(A, pw)
         pw%array = a*pw%array
!$OMP END PARALLEL WORKSHARE

         CALL timestop(handle)
      END SUBROUTINE pw_scale_${kind}$

! **************************************************************************************************
!> \brief pw2 = alpha*pw1 + beta*pw2
!>      alpha defaults to 1, beta defaults to 1
!> \param pw1 ...
!> \param pw2 ...
!> \param alpha ...
!> \param beta ...
!> \param allow_noncompatible_grids ...
!> \par History
!>      JGH (21-Feb-2003) : added reference grid functionality
!>      JGH (01-Dec-2007) : rename and remove complex alpha
!> \author apsi
!> \note
!>      Currently only summing up of respective types allowed,
!>      in order to avoid errors
! **************************************************************************************************
      SUBROUTINE pw_axpy_${kind}$ (pw1, pw2, alpha, beta, allow_noncompatible_grids)

         TYPE(pw_${kind}$_type), INTENT(IN)                 :: pw1
         TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw2
         REAL(KIND=dp), INTENT(in), OPTIONAL                :: alpha, beta
         LOGICAL, INTENT(IN), OPTIONAL                      :: allow_noncompatible_grids

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_axpy_${kind}$'

         INTEGER                                            :: handle
         LOGICAL                                            :: my_allow_noncompatible_grids
         REAL(KIND=dp)                                      :: my_alpha, my_beta
         #:if kind[1:]=='1d'
            INTEGER :: ng1, ng2, ng, i, j
         #:endif

         CALL timeset(routineN, handle)

         IF (pw1%in_space /= pw2%in_space) &
            CPABORT("Both pw grids must be in the same space to add them!")

         my_alpha = 1.0_dp
         IF (PRESENT(alpha)) my_alpha = alpha

         my_beta = 1.0_dp
         IF (PRESENT(beta)) my_beta = beta

         my_allow_noncompatible_grids = .FALSE.
         IF (PRESENT(allow_noncompatible_grids)) my_allow_noncompatible_grids = allow_noncompatible_grids

         IF (.NOT. (ASSOCIATED(pw1%array) .AND. ASSOCIATED(pw2%array))) &
            CPABORT("Both arrays must be associated!")

         IF (my_beta /= 1.0_dp) THEN
            IF (my_beta == 0.0_dp) THEN
               CALL pw_zero(pw2)
            ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw2,my_beta)
               pw2%array = pw2%array*my_beta
!$OMP END PARALLEL WORKSHARE
            END IF
         END IF

         IF (ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN

            IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1, pw2)
               pw2%array = pw2%array + pw1%array
!$OMP END PARALLEL WORKSHARE
            ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw2,pw1,my_alpha)
               pw2%array = pw2%array + my_alpha*pw1%array
!$OMP END PARALLEL WORKSHARE
            END IF

         ELSE IF (pw_compatible(pw1%pw_grid, pw2%pw_grid) .OR. my_allow_noncompatible_grids) THEN
            IF (pw1%pw_grid%spherical .AND. pw2%pw_grid%spherical) THEN
               #:if kind[1:]=='1d'
                  IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                     DO i = 1, ng
                        pw2%array(i) = pw2%array(i) + pw1%array(i)
                     END DO
!$OMP END PARALLEL DO
                  ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i) DEFAULT(NONE) SHARED(pw2,pw1,my_alpha,ng)
                     DO i = 1, ng
                        pw2%array(i) = pw2%array(i) + my_alpha*pw1%array(i)
                     END DO
!$OMP END PARALLEL DO
                  END IF
               #:else
                  CPABORT("Spherical grids are only supported with 1D grids!")
               #:endif
            ELSE IF (.NOT. (pw1%pw_grid%spherical .OR. pw2%pw_grid%spherical)) THEN

               #:if kind[1:]=='1d'
                  ng1 = SIZE(pw1%array)
                  ng2 = SIZE(pw2%array)
                  ng = MIN(ng1, ng2)
               #:endif

               IF (pw1%pw_grid%id_nr == pw2%pw_grid%reference) THEN
                  #:if kind[1:]!='1d'
                     CPABORT("Different reference grids only implemented with 1D grids!")
                  #:else
                     IF (my_alpha == 1.0_dp) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng2, pw1, pw2)
                           DO i = 1, ng2
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + pw1%array(j)
                           END DO
!$OMP END PARALLEL DO
                        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng1, pw1, pw2)
                           DO i = 1, ng1
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + pw1%array(i)
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE IF (my_alpha /= 0.0_dp) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng2, pw1, pw2)
                           DO i = 1, ng2
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + my_alpha*pw1%array(j)
                           END DO
!$OMP END PARALLEL DO
                        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng1, pw1, pw2)
                           DO i = 1, ng1
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + my_alpha*pw1%array(i)
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     END IF
                  #:endif
               ELSE IF (pw2%pw_grid%id_nr == pw1%pw_grid%reference) THEN
                  #:if kind[1:]!='1d'
                     CPABORT("Different reference grids only implemented with 1D grids!")
                  #:else

                     IF (my_alpha == 1.0_dp) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng2, pw1, pw2)
                           DO i = 1, ng2
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + pw1%array(j)
                           END DO
!$OMP END PARALLEL DO
                        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng1, pw1, pw2)
                           DO i = 1, ng1
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + pw1%array(i)
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE IF (my_alpha /= 0.0_dp) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng2, pw1, pw2)
                           DO i = 1, ng2
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + my_alpha*pw1%array(j)
                           END DO
!$OMP END PARALLEL DO
                        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng1, pw1, pw2)
                           DO i = 1, ng1
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + my_alpha*pw1%array(i)
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     END IF
                  #:endif
               ELSE
                  BLOCK
                     #:if kind[1:]=='1d'
                        INTEGER :: minsize, lbound1, lbound2, ubound1, ubound2

                        minsize = MIN(SIZE(pw1%array), SIZE(pw2%array))
                        lbound1 = LBOUND(pw1%array, 1)
                        ubound1 = lbound1 + minsize - 1
                        lbound2 = LBOUND(pw2%array, 1)
                        ubound2 = lbound2 + minsize - 1

                        IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2,lbound1,lbound2,ubound1,ubound2)
                           pw2%array(lbound2:ubound2) = pw2%array(lbound2:ubound2) + pw1%array(lbound1:ubound1)
!$OMP END PARALLEL WORKSHARE
                        ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw2,pw1,my_alpha,lbound1,lbound2,ubound1,ubound2)
                           pw2%array(lbound2:ubound2) = pw2%array(lbound2:ubound2) + my_alpha*pw1%array(lbound1:ubound1)
!$OMP END PARALLEL WORKSHARE
                        END IF
                     #:else
                        INTEGER, DIMENSION(3) :: minsize, lbound1, lbound2, ubound1, ubound2

                        minsize = MIN(SHAPE(pw1%array), SHAPE(pw2%array))
                        lbound1 = LBOUND(pw1%array)
                        ubound1 = lbound1 + minsize - 1
                        lbound2 = LBOUND(pw2%array)
                        ubound2 = lbound2 + minsize - 1

                        IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2,lbound1,lbound2,ubound1,ubound2)
                           pw2%array(lbound2(1):ubound2(1), lbound2(2):ubound2(2), lbound2(3):ubound2(3)) = &
                              pw2%array(lbound2(1):ubound2(1), lbound2(2):ubound2(2), lbound2(3):ubound2(3)) + &
                              pw1%array(lbound1(1):ubound1(1), lbound1(2):ubound1(2), lbound1(3):ubound1(3))
!$OMP END PARALLEL WORKSHARE
                        ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw2,pw1,my_alpha,lbound1,lbound2,ubound1,ubound2)
                           pw2%array(lbound2(1):ubound2(1), lbound2(2):ubound2(2), lbound2(3):ubound2(3)) = &
                              pw2%array(lbound2(1):ubound2(1), lbound2(2):ubound2(2), lbound2(3):ubound2(3)) + &
                              my_alpha*pw1%array(lbound1(1):ubound1(1), lbound1(2):ubound1(2), lbound1(3):ubound1(3))
!$OMP END PARALLEL WORKSHARE
                        END IF
                     #:endif
                  END BLOCK
               END IF
            ELSE
               CPABORT("Either both grids are spherical or non-spherical!")
            END IF
         ELSE
            CPABORT("Incompatible grids!")
         END IF
         CALL timestop(handle)
      END SUBROUTINE pw_axpy_${kind}$

      #:for kind2 in pw_kinds
         #:if kind[1:]==kind2[1:]
! **************************************************************************************************
!> \brief pw_out = pw_out + alpha * pw1 * pw2
!>      alpha defaults to 1
!> \param pw_out ...
!> \param pw1 ...
!> \param pw2 ...
!> \param alpha ...
! **************************************************************************************************
            SUBROUTINE pw_multiply_${kind}$_${kind2}$ (pw_out, pw1, pw2, alpha, beta)

               TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw_out, pw1
               TYPE(pw_${kind2}$_type), INTENT(IN)                 :: pw2
               REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha, beta

               CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_multiply'

               INTEGER                                            :: handle
               REAL(KIND=dp)                                      :: my_alpha

               CALL timeset(routineN, handle)

               my_alpha = 1.0_dp
               IF (PRESENT(alpha)) my_alpha = alpha

               IF (PRESENT(beta)) CALL pw_scale(pw_out, beta)

               IF (ASSOCIATED(pw_out%pw_grid, pw1%pw_grid) .AND. &
                   ASSOCIATED(pw_out%pw_grid, pw2%pw_grid)) THEN

                  IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw_out, pw1, pw2)
                     #:if kind==kind2 or kind[0]=='c'
                        pw_out%array = pw_out%array + pw1%array*pw2%array
                     #:else
                        pw_out%array = pw_out%array + pw1%array*REAL(pw2%array, KIND=dp)
                     #:endif
!$OMP END PARALLEL WORKSHARE
                  ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(my_alpha, pw_out, pw1, pw2)
                     #:if kind==kind2 or kind[0]=='c'
                        pw_out%array = pw_out%array + my_alpha*pw1%array*pw2%array
                     #:else
                        pw_out%array = pw_out%array + my_alpha*pw1%array*REAL(pw2%array, KIND=dp)
                     #:endif
!$OMP END PARALLEL WORKSHARE
                  END IF
               ELSE
                  CPABORT("Incompatible grids!")
               END IF

               CALL timestop(handle)
            END SUBROUTINE pw_multiply_${kind}$_${kind2}$
         #:endif
      #:endfor

! **************************************************************************************************
!> \brief pw_out = pw_out + alpha * pw1 / pw2
!>      alpha defaults to 1, threshold prevents division by zero
!> \param pw_out ...
!> \param pw1 ...
!> \param pw2 ...
!> \param alpha ...
! **************************************************************************************************
      SUBROUTINE pw_divide_${kind}$ (pw_out, pw1, pw2, alpha, beta, threshold)

         TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw_out
         TYPE(pw_${kind}$_type), INTENT(IN)                 :: pw1, pw2
         REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha, beta, threshold

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_divide_${kind}$'

         INTEGER                                            :: handle
         REAL(KIND=dp)                                      :: my_alpha, my_beta, my_threshold

         CALL timeset(routineN, handle)

         my_alpha = 1.0_dp
         IF (PRESENT(alpha)) my_alpha = alpha

         my_beta = 1.0_dp
         IF (PRESENT(beta)) my_beta = beta

         IF (my_beta /= 1.0_dp) CALL pw_scale(pw_out, my_beta)

         IF (ASSOCIATED(pw_out%pw_grid, pw1%pw_grid) .AND. &
             ASSOCIATED(pw_out%pw_grid, pw2%pw_grid)) THEN

            IF (PRESENT(threshold)) THEN
               my_threshold = MAX(threshold, 0.0_dp)
               IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw_out,pw1,pw2,my_threshold)
                  #:if kind[0]=='r'
                     pw_out%array = pw_out%array + pw1%array/SIGN(MAX(ABS(pw2%array), my_threshold), pw2%array)
                  #:else
                     pw_out%array = pw_out%array + pw1%array/MAX(ABS(pw2%array), my_threshold) &
                                    *EXP(CMPLX(0.0_dp, -ATAN2(AIMAG(pw2%array), REAL(pw2%array, KIND=dp)), KIND=dp))
                  #:endif
!$OMP END PARALLEL WORKSHARE
               ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(my_alpha,pw_out,pw1,pw2,my_threshold)
                  #:if kind[0]=='r'
                     pw_out%array = pw_out%array + my_alpha*pw1%array/SIGN(MAX(ABS(pw2%array), my_threshold), pw2%array)
                  #:else
                     pw_out%array = pw_out%array + my_alpha*pw1%array/MAX(ABS(pw2%array), my_threshold) &
                                    *EXP(CMPLX(0.0_dp, -ATAN2(AIMAG(pw2%array), REAL(pw2%array, KIND=dp)), KIND=dp))
                  #:endif
!$OMP END PARALLEL WORKSHARE
               END IF
            ELSE
               IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw_out, pw1, pw2)
                  pw_out%array = pw_out%array + pw1%array/pw2%array
!$OMP END PARALLEL WORKSHARE
               ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(my_alpha, pw_out, pw1, pw2)
                  pw_out%array = pw_out%array + my_alpha*pw1%array/pw2%array
!$OMP END PARALLEL WORKSHARE
               END IF
            END IF
         ELSE
            CPABORT("Incompatible grids!")
         END IF

         CALL timestop(handle)
      END SUBROUTINE pw_divide_${kind}$

! **************************************************************************************************
!> \brief writes a small description of the actual grid
!>      (change to output the data as cube file, maybe with an
!>      optional long_description arg?)
!> \param pw the pw data to output
!> \param unit_nr the unit to output to
!> \author Fawzi Mohamed, Frederick Stein
! **************************************************************************************************
      SUBROUTINE pw_write_${kind}$ (pw, unit_nr)

         TYPE(pw_${kind}$_type), INTENT(in)                 :: pw
         INTEGER, INTENT(in)                                :: unit_nr

         INTEGER                                            :: iostatus

         WRITE (unit=unit_nr, fmt="('<pw>:{ ')", iostat=iostatus)

         WRITE (unit=unit_nr, fmt="(a)", iostat=iostatus) " in_use=${kind}$"
         IF (ASSOCIATED(pw%array)) THEN
            #:if kind.endswith('1d')
               WRITE (unit=unit_nr, fmt="(' c${kind}$=<${type.split(", ")[0]}$(',i8,':',i8,')>')") &
                  LBOUND(pw%array, 1), UBOUND(pw%array, 1)
            #:elif kind.endswith('3d')
               WRITE (unit=unit_nr, fmt="(' c${kind}$=<${type.split(", ")[0]}$(',i8,':',i8,',',i8,':',i8,',',i8,':',i8,')>')") &
                  LBOUND(pw%array, 1), UBOUND(pw%array, 1), LBOUND(pw%array, 2), UBOUND(pw%array, 2), &
                  LBOUND(pw%array, 3), UBOUND(pw%array, 3)
            #:else
               #:stop 'Unknown pw kind:', kind
            #:endif
         ELSE
            WRITE (unit=unit_nr, fmt="(' c${kind}$=*null*')")
         END IF

         SELECT CASE (pw%in_space)
         CASE (REALSPACE)
            WRITE (unit=unit_nr, fmt="(a)", iostat=iostatus) " in_space=REALSPACE"
         CASE (RECIPROCALSPACE)
            WRITE (unit=unit_nr, fmt="(a)", iostat=iostatus) " in_space=RECIPROCALSPACE"
         CASE default
            WRITE (unit=unit_nr, fmt="(' in_space=',i8,',')", iostat=iostatus) &
               pw%in_space
         END SELECT

         WRITE (unit=unit_nr, fmt="(' pw_grid%id_nr=',i8,/,' }')", iostat=iostatus) &
            pw%pw_grid%id_nr

      END SUBROUTINE pw_write_${kind}$

      #:for kind2, type2 in zip(pw_kinds, pw_types)
         #! Mixed 1D/3D case is not meaningful
         #:if kind[1:] == kind2[1:]
! **************************************************************************************************
!> \brief Calculate integral over unit cell for functions in plane wave basis
!>      only returns the real part of it ......
!> \param pw1 ...
!> \param pw2 ...
!> \param sumtype ...
!> \param just_sum ...
!> \param local_only ...
!> \return ...
! **************************************************************************************************
            FUNCTION pw_integral_ab_${kind}$_${kind2}$ (pw1, pw2, sumtype, just_sum, local_only) RESULT(integral_value)

               TYPE(pw_${kind}$_type), INTENT(IN)                 :: pw1
               TYPE(pw_${kind2}$_type), INTENT(IN)                :: pw2
               INTEGER, INTENT(IN), OPTIONAL                      :: sumtype
               LOGICAL, INTENT(IN), OPTIONAL                      :: just_sum, local_only
               REAL(KIND=dp)                                      :: integral_value

               CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_integral_ab_${kind}$_${kind2}$'

               INTEGER                                            :: handle, loc_sumtype
               LOGICAL                                            :: my_just_sum, my_local_only

               CALL timeset(routineN, handle)

               loc_sumtype = do_accurate_sum
               IF (PRESENT(sumtype)) loc_sumtype = sumtype

               my_local_only = .FALSE.
               IF (PRESENT(local_only)) my_local_only = local_only

               IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                  CPABORT("Grids incompatible")
               END IF

               IF (pw1%in_space /= pw2%in_space) &
                  CPABORT("Both grids must be in the same space for integration!")

               my_just_sum = .FALSE.
               IF (PRESENT(just_sum)) my_just_sum = just_sum

               ! since we are only interested into the real part, we drop the imaginary part
               ! the dot product is assumed to be hermitian implying a complex conjugation of the first grid
               IF (loc_sumtype == do_standard_sum) THEN

                  ! Do standard sum
                  #:if kind[0] == 'r'
                     #:if kind2[0] == 'r'
                        integral_value = SUM(pw1%array*pw2%array)
                     #:else
                        integral_value = SUM(pw1%array*REAL(pw2%array, KIND=dp))
                     #:endif
                  #:else
                     #:if kind2[0] == 'r'
                        integral_value = SUM(REAL(pw1%array, KIND=dp)*pw2%array)
                     #:else
                        integral_value = SUM(REAL(CONJG(pw1%array)*pw2%array, KIND=dp))
                     #:endif
                  #:endif

               ELSE

                  ! Do accurate sum

                  #:if kind[0] == 'r'
                     #:if kind2[0] == 'r'
                        integral_value = accurate_sum(pw1%array*pw2%array)
                     #:else
                        integral_value = accurate_sum(pw1%array*REAL(pw2%array, KIND=dp))
                     #:endif
                  #:else
                     #:if kind2[0] == 'r'
                        integral_value = accurate_sum(REAL(pw1%array, KIND=dp)*pw2%array)
                     #:else
                        integral_value = accurate_sum(REAL(CONJG(pw1%array)*pw2%array, KIND=dp))
                     #:endif
                  #:endif

               END IF

               IF (.NOT. my_just_sum) THEN
                  #:if kind.endswith('3d')
                     integral_value = integral_value*pw1%pw_grid%dvol
                  #:else
                     integral_value = integral_value*pw1%pw_grid%vol
                  #:endif
               END IF

               #:if kind=='c1d'
                  IF (pw1%pw_grid%grid_span == HALFSPACE) THEN
                     integral_value = 2.0_dp*integral_value
                     IF (pw1%pw_grid%have_g0) integral_value = integral_value - &
                                                               REAL(CONJG(pw1%array(1))*pw2%array(1), KIND=dp)
                  END IF
               #:endif

               IF (.NOT. my_local_only .AND. pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED) &
                  CALL pw1%pw_grid%para%group%sum(integral_value)
               CALL timestop(handle)

            END FUNCTION pw_integral_ab_${kind}$_${kind2}$
         #:endif
      #:endfor

      #:if kind.endswith('1d')
! **************************************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param pw2 ...
!> \return ...
! **************************************************************************************************
         FUNCTION pw_integral_a2b_${kind}$ (pw1, pw2) RESULT(integral_value)

            TYPE(pw_${kind}$_type), INTENT(IN)                          :: pw1, pw2
            REAL(KIND=dp)                                      :: integral_value

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_integral_a2b_${kind}$'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)

            IF (pw1%in_space /= RECIPROCALSPACE .OR. pw2%in_space /= RECIPROCALSPACE) &
               CPABORT("Both grids must be in reciprocal space for this kind of integration!")

            IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
               CPABORT("Grids incompatible")
            END IF

            #:if kind[0] == 'r'
               integral_value = accurate_sum(pw1%array(:)*pw2%array(:) &
                                             *pw1%pw_grid%gsq(:))
            #:else
               integral_value = accurate_sum(REAL(CONJG(pw1%array(:)) &
                                                  *pw2%array(:), KIND=dp)*pw1%pw_grid%gsq(:))
            #:endif

            IF (pw1%pw_grid%grid_span == HALFSPACE) THEN
               integral_value = 2.0_dp*integral_value
            END IF

            integral_value = integral_value*pw1%pw_grid%vol

            IF (pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED) &
               CALL pw1%pw_grid%para%group%sum(integral_value)
            CALL timestop(handle)

         END FUNCTION pw_integral_a2b_${kind}$
      #:endif

      #:for kind2 in pw_kinds
         #:if kind[1:] == kind2[1:]
! **************************************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param pw2 ...
! **************************************************************************************************
            SUBROUTINE pw_multiply_by_${kind}$_${kind2}$ (pw1, pw2, alpha)
               TYPE(pw_${kind}$_type), INTENT(INOUT)                       :: pw1
               TYPE(pw_${kind2}$_type), INTENT(IN)                          :: pw2
               REAL(KIND=dp), INTENT(IN), OPTIONAL :: alpha

               CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pw_multiply_by_${kind}$_${kind2}$'

               INTEGER                                            :: handle
               REAL(KIND=dp) :: my_alpha

               CALL timeset(routineN, handle)

               IF (pw1%in_space /= pw2%in_space) &
                  CPABORT("Invalid space combination!")
               IF (.NOT. pw_compatible(pw1%pw_grid, pw2%pw_grid)) &
                  CPABORT("Incompatible grids!")

               my_alpha = 1.0_dp
               IF (PRESENT(alpha)) my_alpha = alpha

               IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
                  #:if kind==kind2 or kind[0]=='c'
                     pw1%array = pw1%array*pw2%array
                  #:else
                     pw1%array = pw1%array*REAL(pw2%array, KIND=dp)
                  #:endif
!$OMP END PARALLEL WORKSHARE
               ELSE IF (my_alpha == 0.0_dp) THEN
                  CALL pw_zero(pw1)
               ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2,my_alpha)
                  #:if kind==kind2 or kind[0]=='c'
                     pw1%array = my_alpha*pw1%array*pw2%array
                  #:else
                     pw1%array = my_alpha*pw1%array*REAL(pw2%array, KIND=dp)
                  #:endif
!$OMP END PARALLEL WORKSHARE
               END IF

               CALL timestop(handle)
            END SUBROUTINE pw_multiply_by_${kind}$_${kind2}$
         #:endif
      #:endfor

! **************************************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param pw2 ...
! **************************************************************************************************
      SUBROUTINE pw_divide_by_${kind}$ (pw1, pw2, threshold)
         TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw1
         TYPE(pw_${kind}$_type), INTENT(IN)                 :: pw2
         REAL(KIND=dp), INTENT(IN), OPTIONAL                :: threshold

         CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pw_divide_by_${kind}$'

         INTEGER                                            :: handle
         REAL(KIND=dp)                                      :: my_threshold

         CALL timeset(routineN, handle)

         IF (pw1%in_space /= pw2%in_space) &
            CPABORT("Invalid space combination!")
         IF (.NOT. pw_compatible(pw1%pw_grid, pw2%pw_grid)) &
            CPABORT("Incompatible grids!")

         IF (PRESENT(threshold)) THEN
            my_threshold = MAX(threshold, 0.0_dp)
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2,my_threshold)
            #:if kind[0] == 'r'
               pw1%array = pw1%array/SIGN(MAX(ABS(pw2%array), my_threshold), pw2%array)
            #:else
               pw1%array = pw1%array/MAX(ABS(pw2%array), my_threshold) &
                           *EXP(CMPLX(0.0_dp, -ATAN2(AIMAG(pw2%array), REAL(pw2%array, KIND=dp)), KIND=dp))
            #:endif
!$OMP END PARALLEL WORKSHARE
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
            pw1%array = pw1%array/pw2%array
!$OMP END PARALLEL WORKSHARE
         END IF

         CALL timestop(handle)
      END SUBROUTINE pw_divide_by_${kind}$
   #:endfor

! **************************************************************************************************
!> \brief Multiply all data points with a Gaussian damping factor
!>        Needed for longrange Coulomb potential
!>        V(\vec r)=erf(omega*r)/r
!>        V(\vec g)=\frac{4*\pi}{g**2}*exp(-g**2/omega**2)
!> \param pw ...
!> \param omega ...
!> \par History
!>      Frederick Stein (12-04-2019) created
!> \author Frederick Stein (12-Apr-2019)
!> \note
!>      Performs a Gaussian damping
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_gauss_damp(pw, omega)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gauss_damp'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: omega_2

      CALL timeset(routineN, handle)
      CPASSERT(omega >= 0)

      omega_2 = omega*omega
      omega_2 = 0.25_dp/omega_2

      IF (pw%in_space == RECIPROCALSPACE) THEN

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,omega_2)
         pw%array(:) = pw%array(:)*EXP(-pw%pw_grid%gsq(:)*omega_2)
!$OMP END PARALLEL WORKSHARE
      ELSE
         CPABORT("No suitable data field")
      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_gauss_damp

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of a Gaussian
!> \param pw ...
!> \param omega ...
!> \note
!>      Performs a Gaussian damping
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_log_deriv_gauss(pw, omega)

      TYPE(pw_c1d_type), INTENT(IN)                          :: pw
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_gauss'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: omega_2

      CALL timeset(routineN, handle)
      CPASSERT(omega >= 0)

      omega_2 = omega*omega
      omega_2 = 0.25_dp/omega_2

      IF (pw%in_space == RECIPROCALSPACE) THEN

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,omega_2)
         pw%array(:) = pw%array(:)*(1.0_dp + omega_2*pw%pw_grid%gsq(:))
!$OMP END PARALLEL WORKSHARE

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_log_deriv_gauss

! **************************************************************************************************
!> \brief Multiply all data points with a Gaussian damping factor
!>        Needed for longrange Coulomb potential
!>        V(\vec r)=erf(omega*r)/r
!>        V(\vec g)=\frac{4*\pi}{g**2}*exp(-g**2/omega**2)
!> \param pw ...
!> \param omega ...
!> \par History
!>      Frederick Stein (12-04-2019) created
!> \author Frederick Stein (12-Apr-2019)
!> \note
!>      Performs a Gaussian damping
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_compl_gauss_damp(pw, omega)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_compl_gauss_damp'

      INTEGER                                            :: cnt, handle, i
      REAL(KIND=dp)                                      :: omega_2, tmp

      CALL timeset(routineN, handle)

      omega_2 = omega*omega
      omega_2 = 0.25_dp/omega_2

      IF (pw%in_space == RECIPROCALSPACE) THEN
         cnt = SIZE(pw%array)

!$OMP PARALLEL DO PRIVATE(i, tmp) DEFAULT(NONE) SHARED(cnt, pw,omega_2)
         DO i = 1, cnt
            tmp = -omega_2*pw%pw_grid%gsq(i)
            IF (ABS(tmp) > 1.0E-5_dp) THEN
               pw%array(i) = pw%array(i)*(1.0_dp - EXP(tmp))
            ELSE
               pw%array(i) = pw%array(i)*(tmp + 0.5_dp*tmp*(tmp + (1.0_dp/3.0_dp)*tmp**2))
            END IF
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_compl_gauss_damp

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of the complementary Gaussian damping factor
!> \param pw ...
!> \param omega ...
!> \note
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_log_deriv_compl_gauss(pw, omega)

      TYPE(pw_c1d_type), INTENT(IN)                          :: pw
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_compl_gauss'

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: omega_2, tmp

      CALL timeset(routineN, handle)

      omega_2 = omega*omega
      omega_2 = 0.25_dp/omega_2

      IF (pw%in_space == RECIPROCALSPACE) THEN

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,omega_2)
         DO i = 1, SIZE(pw%array)
            tmp = omega_2*pw%pw_grid%gsq(i)
            ! For too small arguments, use the Taylor polynomial to prevent division by zero
            IF (ABS(tmp) >= 0.003_dp) THEN
               pw%array(i) = pw%array(i)*(1.0_dp - tmp*EXP(-tmp)/(1.0_dp - EXP(-tmp)))
            ELSE
               pw%array(i) = pw%array(i)*(0.5_dp*tmp - tmp**2/12.0_dp)
            END IF
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_log_deriv_compl_gauss

! **************************************************************************************************
!> \brief Multiply all data points with a Gaussian damping factor and mixes it with the original function
!>        Needed for mixed longrange/Coulomb potential
!>        V(\vec r)=(a+b*erf(omega*r))/r
!>        V(\vec g)=\frac{4*\pi}{g**2}*(a+b*exp(-g**2/omega**2))
!> \param pw ...
!> \param omega ...
!> \param scale_coul ...
!> \param scale_long ...
!> \par History
!>      Frederick Stein (16-Dec-2021) created
!> \author Frederick Stein (16-Dec-2021)
!> \note
!>      Performs a Gaussian damping
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_gauss_damp_mix(pw, omega, scale_coul, scale_long)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      REAL(KIND=dp), INTENT(IN)                          :: omega, scale_coul, scale_long

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gauss_damp_mix'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: omega_2

      CALL timeset(routineN, handle)

      omega_2 = omega*omega
      omega_2 = 0.25_dp/omega_2

      IF (pw%in_space == RECIPROCALSPACE) THEN

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw, omega_2, scale_coul, scale_long)
         pw%array(:) = pw%array(:)*(scale_coul + scale_long*EXP(-pw%pw_grid%gsq(:)*omega_2))
!$OMP END PARALLEL WORKSHARE

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_gauss_damp_mix

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of the mixed longrange/Coulomb potential
!>        Needed for mixed longrange/Coulomb potential
!> \param pw ...
!> \param omega ...
!> \param scale_coul ...
!> \param scale_long ...
!> \note
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_log_deriv_mix_cl(pw, omega, scale_coul, scale_long)

      TYPE(pw_c1d_type), INTENT(IN)                          :: pw
      REAL(KIND=dp), INTENT(IN)                          :: omega, scale_coul, scale_long

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_mix_cl'

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: omega_2, tmp

      CALL timeset(routineN, handle)

      omega_2 = omega*omega
      omega_2 = 0.25_dp/omega_2

      IF (pw%in_space == RECIPROCALSPACE) THEN

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,omega_2,scale_long,scale_coul)
         DO i = 1, SIZE(pw%array)
            tmp = omega_2*pw%pw_grid%gsq(i)
            pw%array(i) = pw%array(i)*(1.0_dp + scale_long*tmp*EXP(-tmp)/(scale_coul + scale_long*EXP(-tmp)))
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_log_deriv_mix_cl

! **************************************************************************************************
!> \brief Multiply all data points with a complementary cosine
!>        Needed for truncated Coulomb potential
!>        V(\vec r)=1/r if r<rc, 0 else
!>        V(\vec g)=\frac{4*\pi}{g**2}*(1-cos(g*rc))
!> \param pw ...
!> \param rcutoff ...
!> \par History
!>      Frederick Stein (07-06-2021) created
!> \author Frederick Stein (07-Jun-2021)
!> \note
!>      Multiplies by complementary cosine
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_truncated(pw, rcutoff)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      REAL(KIND=dp), INTENT(IN)                          :: rcutoff

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_truncated'

      INTEGER                                            :: cnt, handle, i
      REAL(KIND=dp)                                      :: tmp

      CALL timeset(routineN, handle)
      CPASSERT(rcutoff >= 0)

      IF (pw%in_space == RECIPROCALSPACE) THEN

         cnt = SIZE(pw%array)

!$OMP PARALLEL DO PRIVATE(i,tmp) DEFAULT(NONE) SHARED(cnt, pw, rcutoff)
         DO i = 1, cnt
            tmp = SQRT(pw%pw_grid%gsq(i))*rcutoff
            IF (tmp >= 0.005_dp) THEN
               pw%array(i) = pw%array(i)*(1.0_dp - COS(tmp))
            ELSE
               pw%array(i) = pw%array(i)*tmp**2/2.0_dp*(1.0 - tmp**2/12.0_dp)
            END IF
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_truncated

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of the complementary cosine
!>        This function is needed for virials using truncated Coulomb potentials
!> \param pw ...
!> \param rcutoff ...
!> \note
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_log_deriv_trunc(pw, rcutoff)

      TYPE(pw_c1d_type), INTENT(IN)                          :: pw
      REAL(KIND=dp), INTENT(IN)                          :: rcutoff

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_trunc'

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: rchalf, tmp

      CALL timeset(routineN, handle)
      CPASSERT(rcutoff >= 0)

      IF (pw%in_space == RECIPROCALSPACE) THEN

         rchalf = 0.5_dp*rcutoff
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,rchalf)
         DO i = 1, SIZE(pw%array)
            tmp = rchalf*SQRT(pw%pw_grid%gsq(i))
            ! For too small arguments, use the Taylor polynomial to prevent division by zero
            IF (ABS(tmp) >= 0.0001_dp) THEN
               pw%array(i) = pw%array(i)*(1.0_dp - tmp/TAN(tmp))
            ELSE
               pw%array(i) = pw%array(i)*tmp**2*(1.0_dp + tmp**2/15.0_dp)/3.0_dp
            END IF
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_log_deriv_trunc

! **************************************************************************************************
!> \brief Calculate the derivative of a plane wave vector
!> \param pw ...
!> \param n ...
!> \par History
!>      JGH (06-10-2002) allow only for inplace derivatives
!> \author JGH (25-Feb-2001)
!> \note
!>      Calculate the derivative dx^n(1) dy^n(2) dz^n(3) PW
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_derive(pw, n)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      INTEGER, DIMENSION(3), INTENT(IN)                  :: n

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_derive'

      COMPLEX(KIND=dp)                                   :: im
      INTEGER                                            :: handle, m

      CALL timeset(routineN, handle)

      CPASSERT(ALL(n >= 0))

      m = SUM(n)
      im = CMPLX(0.0_dp, 1.0_dp, KIND=dp)**m

      IF (pw%in_space == RECIPROCALSPACE) THEN

         IF (n(1) == 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
            pw%array(:) = pw%array(:)*pw%pw_grid%g(1, :)
!$OMP END PARALLEL WORKSHARE
         ELSE IF (n(1) > 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(n, pw)
            pw%array(:) = pw%array(:)*(pw%pw_grid%g(1, :)**n(1))
!$OMP END PARALLEL WORKSHARE
         END IF

         IF (n(2) == 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
            pw%array(:) = pw%array(:)*pw%pw_grid%g(2, :)
!$OMP END PARALLEL WORKSHARE
         ELSE IF (n(2) > 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(n, pw)
            pw%array(:) = pw%array(:)*(pw%pw_grid%g(2, :)**n(2))
!$OMP END PARALLEL WORKSHARE
         END IF

         IF (n(3) == 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
            pw%array(:) = pw%array(:)*pw%pw_grid%g(3, :)
!$OMP END PARALLEL WORKSHARE
         ELSE IF (n(3) > 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(n, pw)
            pw%array(:) = pw%array(:)*(pw%pw_grid%g(3, :)**n(3))
!$OMP END PARALLEL WORKSHARE
         END IF

         ! im can take the values 1, -1, i, -i
         ! skip this if im == 1
         IF (ABS(REAL(im, KIND=dp) - 1.0_dp) > 1.0E-10_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(im, pw)
            pw%array(:) = im*pw%array(:)
!$OMP END PARALLEL WORKSHARE
         END IF

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_derive

! **************************************************************************************************
!> \brief Calculate the Laplacian of a plane wave vector
!> \param pw ...
!> \par History
!>      Frederick Stein (01-02-2022) created
!> \author JGH (25-Feb-2001)
!> \note
!>      Calculate the derivative DELTA PW
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_laplace(pw)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_laplace'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (pw%in_space == RECIPROCALSPACE) THEN

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
         pw%array(:) = -pw%array(:)*pw%pw_grid%gsq(:)
!$OMP END PARALLEL WORKSHARE

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_laplace

! **************************************************************************************************
!> \brief Calculate the tensorial 2nd derivative of a plane wave vector
!> \param pw ...
!> \param pwdr2 ...
!> \param i ...
!> \param j ...
!> \par History
!>      none
!> \author JGH (05-May-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_dr2(pw, pwdr2, i, j)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw, pwdr2
      INTEGER, INTENT(IN)                                :: i, j

      REAL(KIND=dp), PARAMETER :: o3 = 0.33333333333333333333333_dp
      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_dr2'

      INTEGER                                            :: cnt, handle, ig
      REAL(KIND=dp)                                      :: gg

      CALL timeset(routineN, handle)

      IF (pw%in_space == RECIPROCALSPACE) THEN

         cnt = SIZE(pw%array)

         IF (i == j) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ig,gg) SHARED(cnt, i, pw, pwdr2)
            DO ig = 1, cnt
               gg = pw%pw_grid%g(i, ig)**2 - o3*pw%pw_grid%gsq(ig)
               pwdr2%array(ig) = gg*pw%array(ig)
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO PRIVATE (ig) DEFAULT(NONE) SHARED(cnt, i, j, pw, pwdr2)
            DO ig = 1, cnt
               pwdr2%array(ig) = pw%array(ig)*(pw%pw_grid%g(i, ig)*pw%pw_grid%g(j, ig))
            END DO
!$OMP END PARALLEL DO
         END IF

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_dr2

! **************************************************************************************************
!> \brief Calculate the tensorial 2nd derivative of a plane wave vector
!>      and divides by |G|^2. pwdr2_gg(G=0) is put to zero.
!> \param pw ...
!> \param pwdr2_gg ...
!> \param i ...
!> \param j ...
!> \par History
!>      none
!> \author RD (20-Nov-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE
!>      Adapted from pw_dr2
! **************************************************************************************************
   SUBROUTINE pw_dr2_gg(pw, pwdr2_gg, i, j)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw, pwdr2_gg
      INTEGER, INTENT(IN)                                :: i, j

      REAL(KIND=dp), PARAMETER :: o3 = 0.3333333333333333333333333_dp
      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_dr2_gg'

      INTEGER                                            :: cnt, handle, ig
      REAL(KIND=dp)                                      :: gg

      CALL timeset(routineN, handle)

      IF (pw%in_space == RECIPROCALSPACE) THEN

         cnt = SIZE(pw%array)

         IF (i == j) THEN
!$OMP PARALLEL DO PRIVATE (ig) DEFAULT(NONE) PRIVATE(gg) SHARED(cnt, i, pw, pwdr2_gg)
            DO ig = pw%pw_grid%first_gne0, cnt
               gg = pw%pw_grid%g(i, ig)**2 - o3*pw%pw_grid%gsq(ig)
               pwdr2_gg%array(ig) = gg*pw%array(ig)/pw%pw_grid%gsq(ig)
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO PRIVATE (ig) DEFAULT(NONE) SHARED(cnt, i, j, pw, pwdr2_gg)
            DO ig = pw%pw_grid%first_gne0, cnt
               pwdr2_gg%array(ig) = pw%array(ig)*(pw%pw_grid%g(i, ig)*pw%pw_grid%g(j, ig)) &
                                    /pw%pw_grid%gsq(ig)
            END DO
!$OMP END PARALLEL DO
         END IF

         IF (pw%pw_grid%have_g0) pwdr2_gg%array(1) = 0.0_dp

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_dr2_gg

! **************************************************************************************************
!> \brief Multiplies a G-space function with a smoothing factor of the form
!>      f(|G|) = exp((ecut - G^2)/sigma)/(1+exp((ecut - G^2)/sigma))
!> \param pw ...
!> \param ecut ...
!> \param sigma ...
!> \par History
!>      none
!> \author JGH (09-June-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_smoothing(pw, ecut, sigma)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      REAL(KIND=dp), INTENT(IN)                          :: ecut, sigma

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_smoothing'

      INTEGER                                            :: cnt, handle, ig
      REAL(KIND=dp)                                      :: arg, f

      CALL timeset(routineN, handle)

      IF (pw%in_space == RECIPROCALSPACE) THEN

         cnt = SIZE(pw%array)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ig, arg, f) SHARED(cnt, ecut, pw, sigma)
         DO ig = 1, cnt
            arg = (ecut - pw%pw_grid%gsq(ig))/sigma
            f = EXP(arg)/(1 + EXP(arg))
            pw%array(ig) = f*pw%array(ig)
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_smoothing

   #:for kind in pw_kinds
      #:for kind2 in pw_kinds
! **************************************************************************************************
!> \brief Generalize copy of pw types
!> \param pw1 ...
!> \param pw2 ...
!> \param debug ...
!> \par History
!>      JGH (13-Mar-2001) : added gather/scatter cases
!> \author JGH (25-Feb-2001)
!> \note
!>      Copy routine that allows for in_space changes
! **************************************************************************************************
         SUBROUTINE pw_transfer_${kind}$_${kind2}$ (pw1, pw2, debug)

            TYPE(pw_${kind}$_type), INTENT(IN)                          :: pw1
            TYPE(pw_${kind2}$_type), INTENT(INOUT)                       :: pw2
            LOGICAL, INTENT(IN), OPTIONAL                      :: debug

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_transfer_${kind}$_${kind2}$'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)
            !sample peak memory
            CALL m_memory()

            IF (pw1%in_space == REALSPACE .AND. pw2%in_space == REALSPACE) THEN

               #:if kind[1] == kind2[1]
                  CALL pw_copy(pw1, pw2)
               #:elif kind == 'c1d' and kind2 == 'c3d'
                  CALL pw_scatter(pw1, pw2%array)
               #:elif kind == 'c3d' and kind2 == 'c1d'
                  CALL pw_gather(pw2, pw1%array)
               #:else
                  CPABORT("No transfer operation implemented for the pair ${kind}$/${kind2}$!")
               #:endif

            ELSE IF (pw1%in_space == RECIPROCALSPACE .AND. &
                     pw2%in_space == RECIPROCALSPACE) THEN

               #:if kind[1] == kind2[1]
                  CALL pw_copy(pw1, pw2)
               #:elif kind == 'c1d' and kind2 == 'c3d'
                  CALL pw_scatter(pw1, pw2%array)
               #:elif kind == 'c3d' and kind2 == 'c1d'
                  CALL pw_gather(pw2, pw1%array)
               #:else
                  CPABORT("No transfer operation implemented for the pair ${kind}$/${kind2}$!")
               #:endif

            ELSE ! Both grids are in different space, so perform a FFT

               #:if kind != 'r1d' and kind2 != 'r1d' and (kind[0] != 'r' or kind2[0] != 'r')
                  ! FFT needed, all further tests done in fft_wrap_pw1pw2
                  CALL fft_wrap_pw1pw2_${kind}$_${kind2}$ (pw1, pw2, debug)
               #:else
                  MARK_USED(pw1)
                  MARK_USED(pw2)
                  MARK_USED(debug)
                  CPABORT("FFTs not implemented for ${kind}$/${kind2}$!")
               #:endif

            END IF

            CALL timestop(handle)

         END SUBROUTINE pw_transfer_${kind}$_${kind2}$
      #:endfor
   #:endfor

! **************************************************************************************************
!> \brief Gathers the pw vector from a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
   SUBROUTINE pw_gather_s_c1d(pw, c, scale)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: pw
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gather_s'

      INTEGER                                            :: gpt, handle, l, m, n

      CALL timeset(routineN, handle)

      ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                 ngpts => SIZE(pw%pw_grid%gsq), ghat => pw%pw_grid%g_hat)

         IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw, scale)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               pw%array(gpt) = scale*c(l, m, n)
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               pw%array(gpt) = c(l, m, n)
            END DO
!$OMP END PARALLEL DO
         END IF

      END ASSOCIATE

      CALL timestop(handle)

   END SUBROUTINE pw_gather_s_c1d

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE pw_gather_p_c1d(pw, c, scale)

      TYPE(pw_c1d_type), INTENT(INOUT)                   :: pw
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gather_p_c1d'

      INTEGER                                            :: gpt, handle, l, m, mn, n

      CALL timeset(routineN, handle)

      IF (pw%pw_grid%para%mode /= PW_MODE_DISTRIBUTED) THEN
         CPABORT("This grid type is not distributed")
      END IF

      ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                 ngpts => SIZE(pw%pw_grid%gsq), ghat => pw%pw_grid%g_hat, yzq => pw%pw_grid%para%yzq)

         IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw, scale)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               mn = yzq(m, n)
               pw%array(gpt) = scale*c(l, mn)
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               mn = yzq(m, n)
               pw%array(gpt) = c(l, mn)
            END DO
!$OMP END PARALLEL DO
         END IF

      END ASSOCIATE

      CALL timestop(handle)

   END SUBROUTINE pw_gather_p_c1d

! **************************************************************************************************
!> \brief Scatters a pw vector to a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
   SUBROUTINE pw_scatter_s_c1d(pw, c, scale)

      TYPE(pw_c1d_type), INTENT(IN)                      :: pw
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_scatter_s_c1d'

      INTEGER                                            :: gpt, handle, l, m, n

      CALL timeset(routineN, handle)

      IF (pw%in_space /= RECIPROCALSPACE) THEN
         CPABORT("Data has to be in RECIPROCALSPACE")
      END IF

      ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                 ghat => pw%pw_grid%g_hat, ngpts => SIZE(pw%pw_grid%gsq))

         ! should only zero the unused bits (but the zero is needed)
         IF (.NOT. PRESENT(scale)) c = 0.0_dp

         IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw, scale)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               c(l, m, n) = scale*pw%array(gpt)
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               c(l, m, n) = pw%array(gpt)
            END DO
!$OMP END PARALLEL DO
         END IF

      END ASSOCIATE

      IF (pw%pw_grid%grid_span == HALFSPACE) THEN

         ASSOCIATE (mapl => pw%pw_grid%mapl%neg, mapm => pw%pw_grid%mapm%neg, mapn => pw%pw_grid%mapn%neg, &
                    ghat => pw%pw_grid%g_hat, ngpts => SIZE(pw%pw_grid%gsq))

            IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw, scale)
               DO gpt = 1, ngpts
                  l = mapl(ghat(1, gpt)) + 1
                  m = mapm(ghat(2, gpt)) + 1
                  n = mapn(ghat(3, gpt)) + 1
                  c(l, m, n) = scale*CONJG(pw%array(gpt))
               END DO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw)
               DO gpt = 1, ngpts
                  l = mapl(ghat(1, gpt)) + 1
                  m = mapm(ghat(2, gpt)) + 1
                  n = mapn(ghat(3, gpt)) + 1
                  c(l, m, n) = CONJG(pw%array(gpt))
               END DO
!$OMP END PARALLEL DO
            END IF

         END ASSOCIATE

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_scatter_s_c1d

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE pw_scatter_p_c1d(pw, c, scale)

      TYPE(pw_c1d_type), INTENT(IN)                          :: pw
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_scatter_p_c1d'

      INTEGER                                            :: gpt, handle, l, m, mn, n

      CALL timeset(routineN, handle)

      IF (pw%in_space /= RECIPROCALSPACE) THEN
         CPABORT("Data has to be in RECIPROCALSPACE")
      END IF

      IF (pw%pw_grid%para%mode /= PW_MODE_DISTRIBUTED) THEN
         CPABORT("This grid type is not distributed")
      END IF

      ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                 ghat => pw%pw_grid%g_hat, yzq => pw%pw_grid%para%yzq, ngpts => SIZE(pw%pw_grid%gsq))

         IF (.NOT. PRESENT(scale)) c = z_zero

         IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw, scale)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               mn = yzq(m, n)
               c(l, mn) = scale*pw%array(gpt)
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw)
            DO gpt = 1, ngpts
               l = mapl(ghat(1, gpt)) + 1
               m = mapm(ghat(2, gpt)) + 1
               n = mapn(ghat(3, gpt)) + 1
               mn = yzq(m, n)
               c(l, mn) = pw%array(gpt)
            END DO
!$OMP END PARALLEL DO
         END IF

      END ASSOCIATE

      IF (pw%pw_grid%grid_span == HALFSPACE) THEN

         ASSOCIATE (mapm => pw%pw_grid%mapm%neg, mapn => pw%pw_grid%mapn%neg, mapl => pw%pw_grid%mapl%neg, &
                    ghat => pw%pw_grid%g_hat, ngpts => SIZE(pw%pw_grid%gsq), yzq => pw%pw_grid%para%yzq)

            IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw, scale)
               DO gpt = 1, ngpts
                  l = mapl(ghat(1, gpt)) + 1
                  m = mapm(ghat(2, gpt)) + 1
                  n = mapn(ghat(3, gpt)) + 1
                  mn = yzq(m, n)
                  c(l, mn) = scale*CONJG(pw%array(gpt))
               END DO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n) &
!$OMP             SHARED(c, pw)
               DO gpt = 1, ngpts
                  l = mapl(ghat(1, gpt)) + 1
                  m = mapm(ghat(2, gpt)) + 1
                  n = mapn(ghat(3, gpt)) + 1
                  mn = yzq(m, n)
                  c(l, mn) = CONJG(pw%array(gpt))
               END DO
!$OMP END PARALLEL DO
            END IF
         END ASSOCIATE
      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_scatter_p_c1d

   #! r1d is not supported (would require real scatter and copy to complex field, not yet implemented)
   #:for kind in pw_kinds
      #:if kind != 'r1d'
         #:for kind2 in pw_kinds
            #! real to real transforms are very special Fourier transforms not covered here
            #:if kind2 != 'r1d' and (kind[0] != 'r' or kind2[0] != 'r')
! **************************************************************************************************
!> \brief Generic function for 3d FFT of a pw grid
!> \param pw1 ...
!> \param pw2 ...
!> \param debug ...
!> \par History
!>      JGH (30-12-2000): New setup of functions and adaptation to parallelism
!>      JGH (04-01-2001): Moved routine from pws to this module, only covers
!>                        pw_types, no more coefficient types
!> \author apsi
!> \note
!>       fft_wrap_pw1pw2
! **************************************************************************************************
               SUBROUTINE fft_wrap_pw1pw2_${kind}$_${kind2}$ (pw1, pw2, debug)

                  TYPE(pw_${kind}$_type), INTENT(IN)                  :: pw1
                  TYPE(pw_${kind2}$_type), INTENT(INOUT)               :: pw2
                  LOGICAL, INTENT(IN), OPTIONAL                      :: debug

                  CHARACTER(len=*), PARAMETER                        :: routineN = 'fft_wrap_pw1pw2_${kind}$_${kind2}$'

                  COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, POINTER         :: grays
                  COMPLEX(KIND=dp), DIMENSION(:, :, :), CONTIGUOUS, POINTER      :: c_in, c_out
                  INTEGER                                            :: dir, handle, handle2, my_pos, nrays, &
                                                                        out_space, out_unit
                  INTEGER, DIMENSION(:), POINTER                     :: n
                  LOGICAL                                            :: test
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                  LOGICAL                                            :: use_pw_gpu
#endif
                  REAL(KIND=dp)                                      :: norm
                  INTEGER, DIMENSION(3) :: nloc

                  CALL timeset(routineN, handle2)
                  out_unit = cp_logger_get_default_io_unit()
                  CALL timeset(routineN//"_"//TRIM(ADJUSTL(cp_to_string( &
                                                           CEILING(pw1%pw_grid%cutoff/10)*10))), handle)

                  NULLIFY (c_in)
                  NULLIFY (c_out)

                  IF (PRESENT(debug)) THEN
                     test = debug
                  ELSE
                     test = .FALSE.
                  END IF

                  !..check if grids are compatible
                  IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                     IF (pw1%pw_grid%dvol /= pw2%pw_grid%dvol) THEN
                        CPABORT("PW grids not compatible")
                     END IF
                     IF (pw1%pw_grid%para%group /= pw2%pw_grid%para%group) THEN
                        CPABORT("PW grids have not compatible MPI groups")
                     END IF
                  END IF

                  !..prepare input
                  IF (pw1%in_space == REALSPACE) THEN
                     dir = FWFFT
                     norm = 1.0_dp/pw1%pw_grid%ngpts
                     out_space = RECIPROCALSPACE
                  ELSE IF (pw1%in_space == RECIPROCALSPACE) THEN
                     dir = BWFFT
                     norm = 1.0_dp
                     out_space = REALSPACE
                  ELSE
                     CPABORT("Error in space tag")
                  END IF

                  n => pw1%pw_grid%npts

                  IF (pw1%pw_grid%para%mode == PW_MODE_LOCAL) THEN

                     !
                     !..replicated data, use local FFT
                     !

                     IF (test .AND. out_unit > 0) THEN
                        WRITE (out_unit, '(A)') " FFT Protocol "
                        IF (dir == FWFFT) WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "FWFFT"
                        IF (dir == BWFFT) WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "BWFFT"
                        IF (pw1%in_space == REALSPACE) &
                           WRITE (out_unit, '(A,T72,A)') "  in space ", "REALSPACE"
                        IF (pw1%in_space == RECIPROCALSPACE) &
                           WRITE (out_unit, '(A,T66,A)') "  in space ", "RECIPROCALSPACE"
                        IF (out_space == REALSPACE) &
                           WRITE (out_unit, '(A,T72,A)') "  out space ", "REALSPACE"
                        IF (out_space == RECIPROCALSPACE) &
                           WRITE (out_unit, '(A,T66,A)') "  out space ", "RECIPROCALSPACE"
                        WRITE (out_unit, '(A,T66,E15.6)') "  scale factor", norm
                     END IF

                     IF (pw1%in_space == REALSPACE) THEN
                        #:if kind=='r3d' and kind2=='c1d'
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                           CALL pw_gpu_r3dc1d_3d(pw1, pw2, scale=norm)
#elif defined (__PW_FPGA)
                           ALLOCATE (c_out(n(1), n(2), n(3)))
                           ! check if bitstream for the fft size is present
                           ! if not, perform fft3d in CPU
                           IF (pw_fpga_init_bitstream(n) == 1) THEN
                              CALL pw_copy_to_array(pw1, c_out)
#if (__PW_FPGA_SP && __PW_FPGA)
                              CALL pw_fpga_r3dc1d_3d_sp(n, c_out)
#else
                              CALL pw_fpga_r3dc1d_3d_dp(n, c_out)
#endif
                              CALL zdscal(n(1)*n(2)*n(3), norm, c_out, 1)
                              CALL pw_gather(pw2, c_out)
                           ELSE
                              CALL pw_copy_to_array(pw1, c_out)
                              CALL fft3d(dir, n, c_out, scale=norm, debug=test)
                              CALL pw_gather(pw2, c_out)
                           END IF
                           DEALLOCATE (c_out)
#else
                           ALLOCATE (c_out(n(1), n(2), n(3)))
                           CALL pw_copy_to_array_r3d_c3d(pw1, c_out)
                           CALL fft3d(dir, n, c_out, scale=norm, debug=test)
                           CALL pw_gather(pw2, c_out)
                           DEALLOCATE (c_out)
#endif
                        #:else
                           #! If we start/end from/with a complex 3D array, use this one as a buffer
                           #! all other input grids have to be redistributed first into the complex 3D format stored in c_out
                           #:if kind2 == 'c3d'
                              c_out => pw2%array
                           #:else
                              ALLOCATE (c_out(n(1), n(2), n(3)))
                           #:endif
                           #:if kind == 'c3d'
                              c_in => pw1%array
                           #:elif kind == 'c1d'
                              CALL pw_scatter(pw1, c_out)
                           #:elif kind == 'r3d'
                              CALL pw_copy_to_array_r3d_c3d(pw1, c_out)
                           #:else
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif

                           #! Perform the actual Fourier Transform
                           #:if kind=='c3d'
                              CALL fft3d(dir, n, c_in, c_out, scale=norm, debug=test)
                           #:else
                              CALL fft3d(dir, n, c_out, scale=norm, debug=test)
                           #:endif

                           #! Reorder data if necessary
                           #:if kind2 == 'r3d'
                              CALL pw_copy_from_array_r3d_c3d(pw2, c_out)
                              DEALLOCATE (c_out)
                           #:elif kind2 == 'c1d'
                              CALL pw_gather(pw2, c_out)
                              DEALLOCATE (c_out)
                           #:elif kind2 != 'c3d'
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif
                        #:endif
                     ELSE ! Backward transformation
                        #:if kind == 'c1d' and kind2 == 'r3d'
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                           CALL pw_gpu_c1dr3d_3d(pw1, pw2, scale=norm)
#elif defined (__PW_FPGA)
                           ALLOCATE (c_out(n(1), n(2), n(3)))
                           ! check if bitstream for the fft size is present
                           ! if not, perform fft3d in CPU
                           IF (pw_fpga_init_bitstream(n) == 1) THEN
                              CALL pw_scatter(pw1, c_out)
                              ! transform using FPGA
#if (__PW_FPGA_SP && __PW_FPGA)
                              CALL pw_fpga_c1dr3d_3d_sp(n, c_out)
#else
                              CALL pw_fpga_c1dr3d_3d_dp(n, c_out)
#endif
                              CALL zdscal(n(1)*n(2)*n(3), norm, c_out, 1)
                              ! use real part only
                              CALL pw_copy_from_array(pw2, c_out)
                           ELSE
                              IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  PW_SCATTER : 3d -> 1d "
                              CALL pw_scatter(pw1, c_out)
                              ! transform
                              CALL fft3d(dir, n, c_out, scale=norm, debug=test)
                              ! use real part only
                              IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  REAL part "
                              CALL pw_copy_from_array(pw2, c_out)
                           END IF
                           DEALLOCATE (c_out)
#else
                           ALLOCATE (c_out(n(1), n(2), n(3)))
                           IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  PW_SCATTER : 3d -> 1d "
                           CALL pw_scatter(pw1, c_out)
                           ! transform
                           CALL fft3d(dir, n, c_out, scale=norm, debug=test)
                           ! use real part only
                           IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  REAL part "
                           CALL pw_copy_from_array_r3d_c3d(pw2, c_out)
                           DEALLOCATE (c_out)
#endif
                        #:else
                           #! If we start/end from/with a complex 3D array, use this one as a buffer
                           #! all other input grids have to be redistributed first into the complex 3D format stored in c_out
                           #:if kind2 == 'c3d'
                              c_out => pw2%array
                           #:else
                              ALLOCATE (c_out(n(1), n(2), n(3)))
                           #:endif
                           #:if kind == 'c3d'
                              c_in => pw1%array
                           #:elif kind == 'c3d'
                              c_in => pw1%array
                           #:elif kind == 'c1d'
                              CALL pw_scatter(pw1, c_out)
                           #:elif kind == 'r3d'
                              CALL pw_copy_to_array_r3d_c3d(pw1, c_out)
                           #:else
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif

                           #! Perform the actual Fourier Transform
                           #:if kind=='c3d'
                              CALL fft3d(dir, n, c_in, c_out, scale=norm, debug=test)
                           #:else
                              CALL fft3d(dir, n, c_out, scale=norm, debug=test)
                           #:endif

                           #! Reorder data if necessary
                           #:if kind2=='r3d'
                              CALL pw_copy_from_array_r3d_c3d(pw2, c_out)
                              DEALLOCATE (c_out)
                           #:elif kind2=='c1d'
                              CALL pw_gather(pw2, c_out)
                              DEALLOCATE (c_out)
                           #:elif kind2!= 'c3d'
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif
                        #:endif
                     END IF

                     IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') " End of FFT Protocol "

                  ELSE

                     !
                     !..parallel FFT
                     !

                     IF (test .AND. pw1%pw_grid%para%group_head .AND. out_unit > 0) THEN
                        WRITE (out_unit, '(A)') " FFT Protocol "
                        IF (dir == FWFFT) WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "FWFFT"
                        IF (dir == BWFFT) WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "BWFFT"
                        IF (pw1%in_space == REALSPACE) &
                           WRITE (out_unit, '(A,T72,A)') "  in space ", "REALSPACE"
                        IF (pw1%in_space == RECIPROCALSPACE) &
                           WRITE (out_unit, '(A,T66,A)') "  in space ", "RECIPROCALSPACE"
                        IF (out_space == REALSPACE) &
                           WRITE (out_unit, '(A,T72,A)') "  out space ", "REALSPACE"
                        IF (out_space == RECIPROCALSPACE) &
                           WRITE (out_unit, '(A,T66,A)') "  out space ", "RECIPROCALSPACE"
                        WRITE (out_unit, '(A,T66,E15.6)') "  scale factor", norm
                     END IF

                     my_pos = pw1%pw_grid%para%my_pos
                     nrays = pw1%pw_grid%para%nyzray(my_pos)
                     grays => pw1%pw_grid%grays
                     CPASSERT(SIZE(grays, 1) == n(1))
                     CPASSERT(SIZE(grays, 2) == nrays)

                     IF (pw1%in_space == REALSPACE) THEN
                        #:if kind == 'r3d' and kind2 == 'c1d'
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                           ! (no ray dist. is not efficient in CUDA)
                           use_pw_gpu = pw1%pw_grid%para%ray_distribution
                           IF (use_pw_gpu) THEN
                              CALL pw_gpu_r3dc1d_3d_ps(pw1, pw2, scale=norm)
                           ELSE
#endif
!..   prepare input
                              nloc = pw1%pw_grid%npts_local
                              ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                              CALL pw_copy_to_array_r3d_c3d(pw1, c_in)
                              grays = z_zero
                              !..transform
                              IF (pw1%pw_grid%para%ray_distribution) THEN
                                 CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%group, &
                                            pw1%pw_grid%para%rs_group, &
                                            pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                            pw1%pw_grid%para%bo, scale=norm, debug=test)
                              ELSE
                                 CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                                            pw1%pw_grid%para%bo, scale=norm, debug=test)
                              END IF
                              !..prepare output
                              IF (test .AND. pw1%pw_grid%para%group_head .AND. out_unit > 0) &
                                 WRITE (out_unit, '(A)') "  PW_GATHER : 2d -> 1d "
                              CALL pw_gather(pw2, grays)
                              DEALLOCATE (c_in)

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                           END IF
#endif
                        #:else
                           #! If we start from a complex 3D array, use this one as the input buffer
                           #! all other input grids have to be redistributed first into the complex 3D format
                           #! in parallel mode the called routine makes use of its own buffers
                           #:if kind == 'c3d'
                              c_in => pw1%array
                           #:else
                              ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                              #:if kind == 'c3d'
                                 c_in => pw1%array
                              #:elif kind == 'c1d'
                                 CALL pw_scatter(pw1, c_in)
                              #:elif kind == 'r3d'
                                 CALL pw_copy_to_array_r3d_c3d(pw1, c_in)
                              #:else
                                 #:stop 'Unknown kind {}'.format(kind)
                              #:endif
                           #:endif
                           grays = z_zero

                           #! Perform the actual Fourier Transform
                           IF (pw1%pw_grid%para%ray_distribution) THEN
                              CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%group, &
                                         pw1%pw_grid%para%rs_group, &
                                         pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                         pw1%pw_grid%para%bo, scale=norm, debug=test)
                           ELSE
                              CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                                         pw1%pw_grid%para%bo, scale=norm, debug=test)
                           END IF

                           #:if kind != 'c3d'
                              DEALLOCATE (c_in)
                           #:endif

                           #! Reorder data if necessary
                           #:if kind2=='r3d'
                              #:if kind == 'c3d'
                                 #! c_in will be overwritten, so we have to allocate it explicitly to prevent changes in the input array
                                 ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                              #:endif
                              BLOCK
                                 TYPE(pw_c1d_type) :: pw1_c1d
                                 CALL pw_create(pw1_c1d, pw1%pw_grid, RECIPROCALSPACE)
                                 CALL pw_gather(pw1_c1d, grays)
                                 CALL pw_scatter(pw1_c1d, c_in)
                                 CALL pw_release(pw1_c1d)
                              END BLOCK
                              CALL pw_copy_from_array(pw2, c_in)
                              #:if kind == 'c3d'
                                 DEALLOCATE (c_in)
                              #:endif
                           #:elif kind2=='c1d'
                              CALL pw_gather(pw2, grays)
                           #:elif kind2!= 'c3d'
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif
                        #:endif
                     ELSE
                        #:if kind == 'c1d' and kind2 == 'r3d'
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                           ! (no ray dist. is not efficient in CUDA)
                           use_pw_gpu = pw1%pw_grid%para%ray_distribution
                           IF (use_pw_gpu) THEN
                              CALL pw_gpu_c1dr3d_3d_ps(pw1, pw2, scale=norm)
                           ELSE
#endif
                              !..   prepare input
                              IF (test .AND. pw1%pw_grid%para%group_head .AND. out_unit > 0) &
                                 WRITE (out_unit, '(A)') "  PW_SCATTER : 2d -> 1d "
                              grays = z_zero
                              CALL pw_scatter(pw1, grays)
                              nloc = pw2%pw_grid%npts_local
                              ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                              !..transform
                              IF (pw1%pw_grid%para%ray_distribution) THEN
                                 CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%group, &
                                            pw1%pw_grid%para%rs_group, &
                                            pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                            pw1%pw_grid%para%bo, scale=norm, debug=test)
                              ELSE
                                 CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                                            pw1%pw_grid%para%bo, scale=norm, debug=test)
                              END IF
                              !..prepare output
                              IF (test .AND. pw1%pw_grid%para%group_head .AND. out_unit > 0) &
                                 WRITE (out_unit, '(A)') "  Real part "
                              CALL pw_copy_from_array_r3d_c3d(pw2, c_in)
                              DEALLOCATE (c_in)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                           END IF
#endif
                        #:else
                           #! If we start from a complex 3D array, use this one as the input buffer
                           #! all other input grids have to be redistributed first into the complex 3D format
                           #! in parallel mode the called routine makes use of its own buffers
                           #! grays is the grid in reciprocal space
                           grays = z_zero
                           #! c_in becomes the output array here (real space)
                           nloc = pw2%pw_grid%npts_local
                           ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                           #:if kind == 'c3d'
                              BLOCK
                                 TYPE(pw_c1d_type) :: pw1_c1d
                                 CALL pw_create(pw1_c1d, pw1%pw_grid, RECIPROCALSPACE)
                                 CALL pw_gather(pw1_c1d, pw1%array)
                                 CALL pw_scatter(pw1_c1d, grays)
                                 CALL pw_release(pw1_c1d)
                              END BLOCK
                           #:elif kind == 'c1d'
                              CALL pw_scatter(pw1, grays)
                           #:elif kind == 'r3d'
                              CALL pw_copy_to_array_r3d_c3d(pw1, c_in)
                              BLOCK
                                 TYPE(pw_c1d_type) :: pw1_c1d
                                 CALL pw_create(pw1_c1d, pw1%pw_grid, RECIPROCALSPACE)
                                 CALL pw_gather(pw1_c1d, c_in)
                                 CALL pw_scatter(pw1_c1d, grays)
                                 CALL pw_release(pw1_c1d)
                              END BLOCK
                           #:else
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif
                           #:if kind2 == 'c3d'
                              IF (ASSOCIATED(c_in)) DEALLOCATE (c_in)
                              c_in => pw2%array
                           #:endif

                           #! Perform the actual Fourier Transform
                           IF (pw1%pw_grid%para%ray_distribution) THEN
                              CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%group, &
                                         pw1%pw_grid%para%rs_group, &
                                         pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                         pw1%pw_grid%para%bo, scale=norm, debug=test)
                           ELSE
                              CALL fft3d(dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                                         pw1%pw_grid%para%bo, scale=norm, debug=test)
                           END IF

                           #! Reorder data if necessary
                           #:if kind2=='r3d'
                              CALL pw_copy_from_array_r3d_c3d(pw2, c_in)
                              DEALLOCATE (c_in)
                           #:elif kind2=='c1d'
                              CALL pw_gather(pw2, c_in)
                              DEALLOCATE (c_in)
                           #:elif kind2!= 'c3d'
                              #:stop 'Unknown kind {}'.format(kind)
                           #:endif
                        #:endif
                     END IF
                  END IF

                  IF (test .AND. pw1%pw_grid%para%group_head .AND. out_unit > 0) THEN
                     WRITE (out_unit, '(A)') " End of FFT Protocol "
                  END IF

                  CALL timestop(handle)
                  CALL timestop(handle2)

               END SUBROUTINE fft_wrap_pw1pw2_${kind}$_${kind2}$
            #:endif
         #:endfor
      #:endif
   #:endfor

! **************************************************************************************************
!> \brief ...
!> \param grida ...
!> \param gridb ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION pw_compatible(grida, gridb) RESULT(compat)

      TYPE(pw_grid_type), INTENT(IN)                     :: grida, gridb
      LOGICAL                                            :: compat

      compat = .FALSE.

      IF (grida%id_nr == gridb%id_nr) THEN
         compat = .TRUE.
      ELSE IF (grida%reference == gridb%id_nr) THEN
         compat = .TRUE.
      ELSE IF (gridb%reference == grida%id_nr) THEN
         compat = .TRUE.
      END IF

   END FUNCTION pw_compatible

! **************************************************************************************************
!> \brief Calculate the structure factor for point r
!> \param sf ...
!> \param r ...
!> \par History
!>      none
!> \author JGH (05-May-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE
! **************************************************************************************************
   SUBROUTINE pw_structure_factor(sf, r)

      TYPE(pw_c1d_type), INTENT(INOUT)                       :: sf
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_structure_factor'

      INTEGER                                            :: cnt, handle, ig
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)

      IF (sf%in_space == RECIPROCALSPACE) THEN

         cnt = SIZE(sf%array)

!$OMP PARALLEL DO PRIVATE (ig, arg) DEFAULT(NONE) SHARED(cnt, r, sf)
         DO ig = 1, cnt
            arg = DOT_PRODUCT(sf%pw_grid%g(:, ig), r)
            sf%array(ig) = CMPLX(COS(arg), -SIN(arg), KIND=dp)
         END DO
!$OMP END PARALLEL DO

      ELSE

         CPABORT("No suitable data field")

      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_structure_factor

   #:for kind, type in zip(pw_kinds, pw_types)
! **************************************************************************************************
!> \brief ...
!> \param fun ...
!> \param isign ...
!> \param oprt ...
!> \return ...
! **************************************************************************************************
      FUNCTION pw_integrate_function_${kind}$ (fun, isign, oprt, just_sum, local_sum) RESULT(total_fun)

         TYPE(pw_${kind}$_type), INTENT(IN)                          :: fun
         INTEGER, INTENT(IN), OPTIONAL                      :: isign
         CHARACTER(len=*), INTENT(IN), OPTIONAL             :: oprt
         LOGICAL, INTENT(IN), OPTIONAL :: just_sum, local_sum
         REAL(KIND=dp)                                      :: total_fun

         INTEGER                                            :: iop
         LOGICAL :: my_just_sum, my_local_sum

         iop = 0

         IF (PRESENT(oprt)) THEN
            SELECT CASE (oprt)
            CASE ("ABS", "abs")
               iop = 1
            CASE DEFAULT
               CPABORT("Unknown operator")
            END SELECT
         END IF

         my_just_sum = .FALSE.
         IF (PRESENT(just_sum)) my_just_sum = just_sum

         my_local_sum = .FALSE.
         IF (PRESENT(local_sum)) my_local_sum = local_sum

         total_fun = 0.0_dp

         IF (fun%in_space == REALSPACE) THEN
            ! do reduction using maximum accuracy
            IF (iop == 1) THEN
               total_fun = accurate_sum(ABS(fun%array))
            ELSE
               total_fun = accurate_sum(REAL(fun%array, KIND=dp))
            END IF
         ELSE IF (fun%in_space == RECIPROCALSPACE) THEN
            IF (iop == 1) &
               CPABORT("Operator ABS not implemented")
            #:if kind == 'c1d'
               IF (fun%pw_grid%have_g0) total_fun = REAL(fun%array(1), KIND=dp)
            #:else
               CPABORT("Integration in Reciprocalspace with pw kind ${kind}$ not implemented")
            #:endif
         ELSE
            CPABORT("No space defined")
         END IF

         IF (.NOT. my_just_sum) THEN
            IF (fun%in_space == REALSPACE) THEN
               total_fun = total_fun*fun%pw_grid%dvol
            ELSE IF (fun%in_space == RECIPROCALSPACE) THEN
               total_fun = total_fun*fun%pw_grid%vol
            ELSE
               CPABORT("Invalid pw space!")
            END IF
         END IF

         IF (.NOT. my_local_sum .AND. fun%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
            CALL fun%pw_grid%para%group%sum(total_fun)
         END IF

         IF (PRESENT(isign)) THEN
            total_fun = total_fun*SIGN(1._dp, REAL(isign, dp))
         END IF

      END FUNCTION pw_integrate_function_${kind}$

! **************************************************************************************************
!> \brief Initialize pw values using values from a real array data, which might be defined on a
!>        smaller grid than pw but which is contained in it
!> \param pw the pw to initialize
!> \param values the array holding the input data
!> \par History
!>      Created 12.2016
!> \author Nico Holmberg
! **************************************************************************************************
      SUBROUTINE pw_set_array_${kind}$ (pw, values)

         TYPE(pw_${kind}$_type), INTENT(INOUT)                       :: pw
         ${type}$, INTENT(IN)                                        :: values

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_set_array'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

         IF (ANY(SHAPE(pw%array) /= SHAPE(values))) &
            CPABORT("Arrays must have the same shape")

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw, values)
         pw%array = values
!$OMP END PARALLEL WORKSHARE

         CALL timestop(handle)

      END SUBROUTINE pw_set_array_${kind}$

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param value ...
! **************************************************************************************************
      SUBROUTINE pw_set_value_${kind}$ (pw, value)
         TYPE(pw_${kind}$_type), INTENT(IN)                          :: pw
         REAL(KIND=dp), INTENT(IN)                          :: value

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_set_value'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,value)
         pw%array = value
!$OMP END PARALLEL WORKSHARE

         CALL timestop(handle)

      END SUBROUTINE pw_set_value_${kind}$
   #:endfor

END MODULE pw_methods
