!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Subroutines to compute Green functions
! **************************************************************************************************
MODULE negf_green_methods
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_lu_invert,&
                                              cp_cfm_norm,&
                                              cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi,&
                                              z_mone,&
                                              z_one,&
                                              z_zero
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_green_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: sancho_work_matrices_type, sancho_work_matrices_create, sancho_work_matrices_release
   PUBLIC :: do_sancho, negf_contact_self_energy, negf_contact_broadening_matrix, negf_retarded_green_function

   TYPE sancho_work_matrices_type
      ! A_{n+1} = A_n + D_n + E_n
      TYPE(cp_cfm_type), POINTER                         :: a
      ! A0_{n+1} = A0_n + D_n
      TYPE(cp_cfm_type), POINTER                         :: a0
      ! A_inv = A^-1
      TYPE(cp_cfm_type), POINTER                         :: a_inv
      ! B_{n+1} = - B_n A_n^-1 B_n \equiv B_n A_n^-1 B_n
      TYPE(cp_cfm_type), POINTER                         :: b
      ! C_{n+1} = - C_n A_n^-1 C_n \equiv C_n A_n^-1 C_n
      TYPE(cp_cfm_type), POINTER                         :: c
      ! D_n = - B_n A_n^-1 C_n
      TYPE(cp_cfm_type), POINTER                         :: d
      ! E_n = - C_n A_n^-1 B_n
      TYPE(cp_cfm_type), POINTER                         :: e
      ! a scratch area for matrix multiplication
      TYPE(cp_cfm_type), POINTER                         :: scratch
   END TYPE sancho_work_matrices_type

CONTAINS
! **************************************************************************************************
!> \brief Create work matrices required for the Lopez-Sancho algorithm.
!> \param work      work matrices to create (allocated and initialised on exit)
!> \param fm_struct dense matrix structure
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE sancho_work_matrices_create(work, fm_struct)
      TYPE(sancho_work_matrices_type), INTENT(inout)     :: work
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CHARACTER(len=*), PARAMETER :: routineN = 'sancho_work_matrices_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncols, nrows

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(fm_struct))

      CALL cp_fm_struct_get(fm_struct, nrow_global=nrows, ncol_global=ncols)
      CPASSERT(nrows == ncols)

      NULLIFY (work%a, work%a0, work%a_inv, work%b, work%c, work%d, work%e, work%scratch)
      CALL cp_cfm_create(work%a, fm_struct)
      CALL cp_cfm_create(work%a0, fm_struct)
      CALL cp_cfm_create(work%a_inv, fm_struct)
      CALL cp_cfm_create(work%b, fm_struct)
      CALL cp_cfm_create(work%c, fm_struct)
      CALL cp_cfm_create(work%d, fm_struct)
      CALL cp_cfm_create(work%e, fm_struct)
      CALL cp_cfm_create(work%scratch, fm_struct)

      CALL timestop(handle)
   END SUBROUTINE sancho_work_matrices_create

! **************************************************************************************************
!> \brief Release work matrices.
!> \param work   work matrices to release
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE sancho_work_matrices_release(work)
      TYPE(sancho_work_matrices_type), INTENT(inout)     :: work

      CHARACTER(len=*), PARAMETER :: routineN = 'sancho_work_matrices_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(work%a)) CALL cp_cfm_release(work%a)
      IF (ASSOCIATED(work%a0)) CALL cp_cfm_release(work%a0)
      IF (ASSOCIATED(work%a_inv)) CALL cp_cfm_release(work%a_inv)
      IF (ASSOCIATED(work%b)) CALL cp_cfm_release(work%b)
      IF (ASSOCIATED(work%c)) CALL cp_cfm_release(work%c)
      IF (ASSOCIATED(work%d)) CALL cp_cfm_release(work%d)
      IF (ASSOCIATED(work%e)) CALL cp_cfm_release(work%e)
      IF (ASSOCIATED(work%scratch)) CALL cp_cfm_release(work%scratch)

      CALL timestop(handle)
   END SUBROUTINE sancho_work_matrices_release

! **************************************************************************************************
!> \brief Iterative method to compute a retarded surface Green's function at the point omega.
!> \param g_surf   computed retarded surface Green's function (initialised on exit)
!> \param omega    argument of the Green's function
!> \param h0       diagonal block of the Kohn-Sham matrix (must be Hermitian)
!> \param s0       diagonal block of the overlap matrix (must be Hermitian)
!> \param h1       off-fiagonal block of the Kohn-Sham matrix
!> \param s1       off-fiagonal block of the overlap matrix
!> \param conv     convergence threshold
!> \param transp   flag which indicates that the matrices h1 and s1 matrices should be transposed
!> \param work     a set of work matrices needed to compute the surface Green's function
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE do_sancho(g_surf, omega, h0, s0, h1, s1, conv, transp, work)
      TYPE(cp_cfm_type), POINTER                         :: g_surf
      COMPLEX(kind=dp), INTENT(in)                       :: omega
      TYPE(cp_fm_type), POINTER                          :: h0, s0, h1, s1
      REAL(kind=dp), INTENT(in)                          :: conv
      LOGICAL, INTENT(in)                                :: transp
      TYPE(sancho_work_matrices_type), INTENT(in)        :: work

      CHARACTER(len=*), PARAMETER :: routineN = 'do_sancho', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncols, nrows

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(g_surf, nrow_global=nrows, ncol_global=ncols)

      IF (debug_this_module) THEN
         CPASSERT(nrows == ncols)
      END IF

      ! A^{ret.}_0 = omega * s0 - h0
      CALL cp_fm_to_cfm(msourcer=s0, mtarget=work%a)
      CALL cp_cfm_scale_and_add_fm(omega, work%a, z_mone, h0)
      ! A0^{ret.}_0 = A^{ret.}_0
      CALL cp_cfm_to_cfm(work%a, work%a0)

      ! B^{ret.}_0 = omega * s1 - h1
      ! C^{ret.}_0 = B^{ret.}_0^T
      IF (transp) THEN
         ! beta = omega * s1 - h1
         CALL cp_fm_to_cfm(msourcer=s1, mtarget=work%c)
         CALL cp_cfm_scale_and_add_fm(omega, work%c, z_mone, h1)

         ! alpha = omega * s1' - h1'
         CALL cp_cfm_transpose(matrix=work%c, trans='T', matrixt=work%b)
      ELSE
         ! alpha  = omega * s1 - h1
         CALL cp_fm_to_cfm(msourcer=s1, mtarget=work%b)
         CALL cp_cfm_scale_and_add_fm(omega, work%b, z_mone, h1)

         ! beta = omega * s1' - h1'
         CALL cp_cfm_transpose(matrix=work%b, trans='T', matrixt=work%c)
      END IF

      ! actually compute the Green's function
      DO WHILE (cp_cfm_norm(work%b, 'M') + cp_cfm_norm(work%c, 'M') > conv)
         ! A_n^-1
         CALL cp_cfm_to_cfm(work%a, work%a_inv)
         CALL cp_cfm_lu_invert(work%a_inv)

         ! scratch <- A_n^-1 * B_n
         CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, z_one, work%a_inv, work%b, z_zero, work%scratch)
         ! E_n = - C_n A_n^-1 B_n
         CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, z_mone, work%c, work%scratch, z_zero, work%e)
         ! g_surf <- B_{n+1} = B_n A_n^-1 B_n
         ! keep B_n, as we need it to compute the matrix product B_n A_n^-1 C_n
         CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, z_one, work%b, work%scratch, z_zero, g_surf)

         ! scratch <- A_n^-1 * C_n
         CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, z_one, work%a_inv, work%c, z_zero, work%scratch)
         ! D_n = - B_n A_n^-1 C_n
         CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, z_mone, work%b, work%scratch, z_zero, work%d)
         ! we do not need B_n any longer, so the matrix B now holds the B_{n+1} matrix
         CALL cp_cfm_to_cfm(g_surf, work%b)
         ! C_{n+1} = C_n A_n^-1 C_n
         CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, z_one, work%c, work%scratch, z_zero, g_surf)
         CALL cp_cfm_to_cfm(g_surf, work%c)

         ! A0_{n+1} = A0_n + D_n = A0_n - B_n A_n^-1 C_n
         CALL cp_cfm_scale_and_add(z_one, work%a0, z_one, work%d)

         ! A_{n+1} = A0_n + D_n + E_n = A_n - B_n A_n^-1 C_n - C_n A_n^-1 B_n
         CALL cp_cfm_scale_and_add(z_one, work%a, z_one, work%d)
         CALL cp_cfm_scale_and_add(z_one, work%a, z_one, work%e)
      END DO

      ! g_surf = A0_n^-1
      CALL cp_cfm_to_cfm(work%a0, g_surf)
      CALL cp_cfm_lu_invert(g_surf)

      CALL timestop(handle)
   END SUBROUTINE do_sancho

! **************************************************************************************************
!> \brief Compute the contact self energy at point 'omega' as
!>   self_energy_C = [omega * S_SC0 - KS_SC0] * g_surf_c(omega - v_C) * [omega * S_SC0 - KS_SC0]^T,
!>   where C stands for the left (L) / right (R) contact.
!> \param self_energy_c contact self energy (initialised on exit)
!> \param omega         energy point where the contact self energy needs to be computed
!> \param g_surf_c      contact surface Green's function
!> \param h_sc0         scattering region -- contact off-diagonal block of the Kohn-Sham matrix
!> \param s_sc0         scattering region -- contact off-diagonal block of the overlap matrix
!> \param zwork1        complex work matrix of the same shape as s_sc0
!> \param zwork2        another complex work matrix of the same shape as s_sc0
!> \param transp        flag which indicates that transposed matrices (KS_C0S and S_C0S)
!>                      were actually passed
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_contact_self_energy(self_energy_c, omega, g_surf_c, h_sc0, s_sc0, zwork1, zwork2, transp)
      TYPE(cp_cfm_type), POINTER                         :: self_energy_c
      COMPLEX(kind=dp), INTENT(in)                       :: omega
      TYPE(cp_cfm_type), POINTER                         :: g_surf_c
      TYPE(cp_fm_type), POINTER                          :: h_sc0, s_sc0
      TYPE(cp_cfm_type), POINTER                         :: zwork1, zwork2
      LOGICAL, INTENT(in)                                :: transp

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_contact_self_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao_contact, nao_scattering

      CALL timeset(routineN, handle)

      ! zwork1 = omega * S_SC0   - KS_SC0     if transp = .FALSE., or
      !        = omega * S_SC0^T - KS_SC0^T   if transp = .TRUE.
      CALL cp_fm_to_cfm(msourcer=s_sc0, mtarget=zwork1)
      CALL cp_cfm_scale_and_add_fm(omega, zwork1, z_mone, h_sc0)

      IF (transp) THEN
         CALL cp_fm_get_info(s_sc0, nrow_global=nao_contact, ncol_global=nao_scattering)

         ! zwork2 = g_surf_c * [omega * S_SC0^T - KS_SC0^T]
         CALL cp_cfm_gemm('N', 'N', nao_contact, nao_scattering, nao_contact, z_one, g_surf_c, zwork1, z_zero, zwork2)
         ! [omega * S_SC0^T - KS_SC0^T]^T * g_surf_c * [omega * S_SC0^T - KS_SC0^T]
         CALL cp_cfm_gemm('T', 'N', nao_scattering, nao_scattering, nao_contact, z_one, zwork1, zwork2, z_zero, self_energy_c)
      ELSE
         CALL cp_fm_get_info(s_sc0, nrow_global=nao_scattering, ncol_global=nao_contact)

         ! zwork2 = [omega * S_SC0 - KS_SC0] * g_surf_c
         CALL cp_cfm_gemm('N', 'N', nao_scattering, nao_contact, nao_contact, z_one, zwork1, g_surf_c, z_zero, zwork2)
         ! [omega * S_SC0 - KS_SC0] * g_surf_c * [omega * S_SC0 - KS_SC0]^T
         CALL cp_cfm_gemm('N', 'T', nao_scattering, nao_scattering, nao_contact, z_one, zwork2, zwork1, z_zero, self_energy_c)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_contact_self_energy

! **************************************************************************************************
!> \brief Compute contact broadening matrix as
!>   gamma_C = i (self_energy_c^{ret.} - (self_energy_c^{ret.})^C)
!> \param gamma_c            broadening matrix (initialised on exit)
!> \param self_energy_c      retarded contact self-energy
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_contact_broadening_matrix(gamma_c, self_energy_c)
      TYPE(cp_cfm_type), POINTER                         :: gamma_c, self_energy_c

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_contact_broadening_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! gamma_contact = i (self_energy_contact^{ret.} - self_energy_contact^{adv.}) .
      !
      ! With no k-points, the Hamiltonian matrix is real-values, so
      ! self_energy_contact^{adv.} = self_energy_contact^{ret.}^C ,
      ! The above identity allows us to use a simplified expression for the broadening matrix:
      ! gamma_contact = i [self_energy_contact - self_energy_contact^C] .

      CALL cp_cfm_transpose(self_energy_c, 'C', gamma_c)
      CALL cp_cfm_scale_and_add(z_mone, gamma_c, z_one, self_energy_c)
      CALL cp_cfm_scale(gaussi, gamma_c)

      CALL timestop(handle)
   END SUBROUTINE negf_contact_broadening_matrix

! **************************************************************************************************
!> \brief Compute the retarded Green's function at point 'omega' as
!>   G_S^{ret.} = [ omega * S_S - KS_S - \sum_{contact} self_energy_{contact}^{ret.}]^{-1}.
!> \param g_ret_s            complex matrix to store the computed retarded Green's function
!> \param omega              energy point where the retarded Green's function needs to be computed
!> \param self_energy_ret_sum sum of contact self-energies
!> \param h_s                Kohn-Sham matrix block of the scattering region
!> \param s_s                overlap matrix block of the scattering region
!> \param v_hartree_s        contribution to the Kohn-Sham matrix from the external potential
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_retarded_green_function(g_ret_s, omega, self_energy_ret_sum, h_s, s_s, v_hartree_s)
      TYPE(cp_cfm_type), POINTER                         :: g_ret_s
      COMPLEX(kind=dp), INTENT(in)                       :: omega
      TYPE(cp_cfm_type), POINTER                         :: self_energy_ret_sum
      TYPE(cp_fm_type), POINTER                          :: h_s, s_s, v_hartree_s

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_retarded_green_function', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! g_ret_s = [ omega * S_S - H_S - self_energy_left - self_energy_right]^{-1}
      !
      ! omega * S_S - H_S - V_Hartree
      CALL cp_fm_to_cfm(msourcer=s_s, mtarget=g_ret_s)
      CALL cp_cfm_scale_and_add_fm(omega, g_ret_s, z_mone, h_s)
      IF (ASSOCIATED(v_hartree_s)) &
         CALL cp_cfm_scale_and_add_fm(z_one, g_ret_s, z_one, v_hartree_s)

      ! g_ret_s = [omega * S_S - H_S - \sum_{contact} self_energy_{contact}^{ret.} ]^-1
      CALL cp_cfm_scale_and_add(z_one, g_ret_s, z_mone, self_energy_ret_sum)

      CALL cp_cfm_lu_invert(g_ret_s)

      CALL timestop(handle)
   END SUBROUTINE negf_retarded_green_function
END MODULE negf_green_methods

