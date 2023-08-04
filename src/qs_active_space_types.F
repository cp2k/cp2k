!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The types needed for the calculation of active space Hamiltonians
!> \par History
!>      04.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_active_space_types

   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_csr_destroy,&
                                              dbcsr_csr_p_type,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_active_space_types'

   PUBLIC :: active_space_type, eri_type, eri_type_eri_element_func
   PUBLIC :: create_active_space_type, release_active_space_type
   PUBLIC :: csr_idx_to_combined, csr_idx_from_combined, get_irange_csr

! **************************************************************************************************
!> \brief Quantities needed for AS determination
!> \author JGH
! **************************************************************************************************
   TYPE eri_gpw_type
      LOGICAL                       :: redo_poisson = .FALSE.
      REAL(KIND=dp)                 :: cutoff = 0.0_dp
      REAL(KIND=dp)                 :: rel_cutoff = 0.0_dp
      REAL(KIND=dp)                 :: eps_grid = 0.0_dp
      INTEGER                       :: print_level = 0
      LOGICAL                       :: store_wfn = .FALSE.
   END TYPE eri_gpw_type

   TYPE eri_type
      INTEGER                       :: method = 0
      INTEGER                       :: OPERATOR = 0
      REAL(KIND=dp)                 :: operator_parameter = 0.0_dp
      INTEGER, DIMENSION(3)         :: periodicity = 0
      REAL(KIND=dp)                 :: cutoff_radius = 0.0_dp
      REAL(KIND=dp)                 :: eps_integral = 0.0_dp
      TYPE(eri_gpw_type)            :: eri_gpw = eri_gpw_type()
      TYPE(dbcsr_csr_p_type), &
         DIMENSION(:), POINTER      :: eri => NULL()
      INTEGER                       :: norb = 0

   CONTAINS
      PROCEDURE :: eri_foreach => eri_type_eri_foreach
   END TYPE eri_type

! **************************************************************************************************
!> \brief Abstract function object for the `eri_type_eri_foreach` method
! **************************************************************************************************
   TYPE, ABSTRACT :: eri_type_eri_element_func
   CONTAINS
      PROCEDURE(eri_type_eri_element_func_interface), DEFERRED :: func
   END TYPE eri_type_eri_element_func

   TYPE active_space_type
      INTEGER                       :: nelec_active = 0
      INTEGER, POINTER, DIMENSION(:, :) :: active_orbitals => NULL()
      INTEGER, POINTER, DIMENSION(:, :) :: inactive_orbitals => NULL()
      INTEGER                       :: nelec_inactive = 0
      INTEGER, DIMENSION(2)         :: nelec_inactive_spinwise = 0
      INTEGER, DIMENSION(2)         :: nelec_total = 0
      INTEGER                       :: multiplicity = 0
      INTEGER                       :: nspins = 0
      LOGICAL                       :: molecule = .FALSE.
      INTEGER                       :: model = 0
      REAL(KIND=dp)                 :: energy_total = 0.0_dp
      REAL(KIND=dp)                 :: energy_ref = 0.0_dp
      REAL(KIND=dp)                 :: energy_inactive = 0.0_dp
      REAL(KIND=dp)                 :: energy_active = 0.0_dp
      LOGICAL                       :: read_p_act = .FALSE.
      CHARACTER(LEN=default_path_length)        :: p_act_filename = ""
      TYPE(mo_set_type), DIMENSION(:), POINTER  :: mos_active => NULL()
      TYPE(mo_set_type), DIMENSION(:), POINTER  :: mos_inactive => NULL()
      TYPE(eri_type)                            :: eri = eri_type()
      TYPE(cp_fm_type), DIMENSION(:), POINTER   :: p_active => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER   :: ks_sub => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER   :: vxc_sub => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER   :: h_sub => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER   :: fock_sub => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: pmat_inactive => NULL()
   END TYPE active_space_type

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief The function signature to be implemented by a child of `eri_type_eri_element_func`
!> \param this object reference
!> \param i i-index
!> \param j j-index
!> \param k k-index
!> \param l l-index
!> \param val value of the integral at (i,j,k.l)
!> \return True if the ERI foreach loop should continue, false, if not
! **************************************************************************************************
      LOGICAL FUNCTION eri_type_eri_element_func_interface(this, i, j, k, l, val)
         IMPORT :: eri_type_eri_element_func, dp
         CLASS(eri_type_eri_element_func), INTENT(inout) :: this
         INTEGER, INTENT(in)                             :: i, j, k, l
         REAL(KIND=dp), INTENT(in)                       :: val
      END FUNCTION eri_type_eri_element_func_interface
   END INTERFACE

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Creates an active space environment type, nullifying all quantities.
!> \param active_space_env the active space environment to be initialized
! **************************************************************************************************
   SUBROUTINE create_active_space_type(active_space_env)
      TYPE(active_space_type), POINTER                   :: active_space_env

      IF (ASSOCIATED(active_space_env)) THEN
         CALL release_active_space_type(active_space_env)
      END IF

      ALLOCATE (active_space_env)
      NULLIFY (active_space_env%active_orbitals, active_space_env%inactive_orbitals)
      NULLIFY (active_space_env%mos_active, active_space_env%mos_inactive)
      NULLIFY (active_space_env%ks_sub, active_space_env%p_active)
      NULLIFY (active_space_env%vxc_sub, active_space_env%h_sub)
      NULLIFY (active_space_env%fock_sub, active_space_env%pmat_inactive)

   END SUBROUTINE create_active_space_type

! **************************************************************************************************
!> \brief Releases all quantities in the active space environment.
!> \param active_space_env the active space environment to be released
! **************************************************************************************************
   SUBROUTINE release_active_space_type(active_space_env)
      TYPE(active_space_type), POINTER                   :: active_space_env

      INTEGER                                            :: imo

      IF (ASSOCIATED(active_space_env)) THEN

         IF (ASSOCIATED(active_space_env%active_orbitals)) THEN
            DEALLOCATE (active_space_env%active_orbitals)
         END IF

         IF (ASSOCIATED(active_space_env%inactive_orbitals)) THEN
            DEALLOCATE (active_space_env%inactive_orbitals)
         END IF

         IF (ASSOCIATED(active_space_env%mos_active)) THEN
            DO imo = 1, SIZE(active_space_env%mos_active)
               CALL deallocate_mo_set(active_space_env%mos_active(imo))
            END DO
            DEALLOCATE (active_space_env%mos_active)
         END IF

         IF (ASSOCIATED(active_space_env%mos_inactive)) THEN
            DO imo = 1, SIZE(active_space_env%mos_inactive)
               CALL deallocate_mo_set(active_space_env%mos_inactive(imo))
            END DO
            DEALLOCATE (active_space_env%mos_inactive)
         END IF

         CALL release_eri_type(active_space_env%eri)

         CALL cp_fm_release(active_space_env%p_active)
         CALL cp_fm_release(active_space_env%ks_sub)
         CALL cp_fm_release(active_space_env%vxc_sub)
         CALL cp_fm_release(active_space_env%h_sub)
         CALL cp_fm_release(active_space_env%fock_sub)

         IF (ASSOCIATED(active_space_env%pmat_inactive)) &
            CALL dbcsr_deallocate_matrix_set(active_space_env%pmat_inactive)

         DEALLOCATE (active_space_env)
      END IF

   END SUBROUTINE release_active_space_type

! **************************************************************************************************
!> \brief Releases the ERI environment type.
!> \param eri_env the ERI environment to be released
! **************************************************************************************************
   SUBROUTINE release_eri_type(eri_env)
      TYPE(eri_type)                                     :: eri_env

      INTEGER                                            :: i

      IF (ASSOCIATED(eri_env%eri)) THEN

         DO i = 1, SIZE(eri_env%eri)
            CALL dbcsr_csr_destroy(eri_env%eri(i)%csr_mat)
            DEALLOCATE (eri_env%eri(i)%csr_mat)
         END DO
         DEALLOCATE (eri_env%eri)

      END IF

   END SUBROUTINE release_eri_type

! **************************************************************************************************
!> \brief calculates combined index (ij)
!> \param i Index j
!> \param j Index i
!> \param n Dimension in i or j direction
!> \returns The combined index
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   INTEGER FUNCTION csr_idx_to_combined(i, j, n) RESULT(ij)
      INTEGER, INTENT(IN)                                :: i, j, n

      CPASSERT(i <= j)
      CPASSERT(i <= n)
      CPASSERT(j <= n)

      ij = (i - 1)*n - ((i - 1)*(i - 2))/2 + (j - i + 1)

      CPASSERT(ij <= (n*(n + 1))/2)

   END FUNCTION csr_idx_to_combined

! **************************************************************************************************
!> \brief extracts indices i and j from combined index ij
!> \param ij The combined index
!> \param n Dimension in i or j direction
!> \param i Resulting i index
!> \param j Resulting j index
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE csr_idx_from_combined(ij, n, i, j)
      INTEGER, INTENT(IN)                                :: ij, n
      INTEGER, INTENT(OUT)                               :: i, j

      INTEGER                                            :: m, m0

      m = MAX(ij/n, 1)
      DO i = m, n
         m0 = (i - 1)*n - ((i - 1)*(i - 2))/2
         j = ij - m0 + i - 1
         IF (j <= n) EXIT
      END DO

      CPASSERT(i > 0 .AND. i <= n)
      CPASSERT(j > 0 .AND. j <= n)
      CPASSERT(i <= j)

   END SUBROUTINE csr_idx_from_combined

! **************************************************************************************************
!> \brief calculates index range for processor in group mp_group
!> \param nindex the number of indices
!> \param mp_group message-passing group ID
!> \return a range tuple
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   FUNCTION get_irange_csr(nindex, mp_group) RESULT(irange)
      USE message_passing, ONLY: mp_comm_type
      INTEGER, INTENT(IN)                                :: nindex
      TYPE(mp_comm_type), INTENT(IN)                     :: mp_group
      INTEGER, DIMENSION(2)                              :: irange

      REAL(KIND=dp)                                      :: rat

      ASSOCIATE (numtask => mp_group%num_pe, taskid => mp_group%mepos)

         IF (numtask == 1 .AND. taskid == 0) THEN
            irange(1) = 1
            irange(2) = nindex
         ELSEIF (numtask >= nindex) THEN
            IF (taskid >= nindex) THEN
               irange(1) = 1
               irange(2) = 0
            ELSE
               irange(1) = taskid + 1
               irange(2) = taskid + 1
            END IF
         ELSE
            rat = REAL(nindex, KIND=dp)/REAL(numtask, KIND=dp)
            irange(1) = NINT(rat*taskid) + 1
            irange(2) = NINT(rat*taskid + rat)
         END IF
      END ASSOCIATE

   END FUNCTION get_irange_csr

! **************************************************************************************************
!> \brief Calls the provided function for each element in the ERI
!> \param this object reference
!> \param nspin The spin number
!> \param active_orbitals the active orbital indices
!> \param fobj The function object from which to call `func(i, j, k, l, val)`
!> \param spin1 the first spin value
!> \param spin2 the second spin value
!> \par History
!>      04.2016 created [JHU]
!>      06.2016 factored out from qs_a_s_methods:fcidump [TMU]
!> \note Calls MPI, must be executed on all ranks.
! **************************************************************************************************
   SUBROUTINE eri_type_eri_foreach(this, nspin, active_orbitals, fobj, spin1, spin2)
      CLASS(eri_type), INTENT(in)              :: this
      CLASS(eri_type_eri_element_func)         :: fobj
      INTEGER, DIMENSION(:, :), INTENT(IN)     :: active_orbitals
      INTEGER, OPTIONAL                        :: spin1, spin2
      INTEGER                                  :: i1, i12, i12l, i2, i3, i34, i34l, i4, m1, m2, m3, m4, &
                                                  irange(2), irptr, nspin, nindex, nmo
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: colind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: erival
      REAL(KIND=dp)                            :: erint
      TYPE(mp_comm_type) :: mp_group

      IF (.NOT. PRESENT(spin1)) THEN
         spin1 = nspin
      END IF
      IF (.NOT. PRESENT(spin2)) THEN
         spin2 = nspin
      END IF

      ASSOCIATE (eri => this%eri(nspin)%csr_mat, norb => this%norb)
         nindex = (norb*(norb + 1))/2
         CALL mp_group%set_handle(eri%mp_group%get_handle())
         irange = get_irange_csr(nindex, mp_group)
         ALLOCATE (erival(nindex), colind(nindex))

         nmo = SIZE(active_orbitals, 1)

         DO m1 = 1, nmo
            i1 = active_orbitals(m1, spin1)
            DO m2 = m1, nmo
               i2 = active_orbitals(m2, spin1)
               i12 = csr_idx_to_combined(i1, i2, norb)

               IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
                  i12l = i12 - irange(1) + 1
                  irptr = eri%rowptr_local(i12l)
                  nindex = eri%nzerow_local(i12l)
                  colind(1:nindex) = eri%colind_local(irptr:irptr + nindex - 1)
                  erival(1:nindex) = eri%nzval_local%r_dp(irptr:irptr + nindex - 1)
               ELSE
                  erival = 0.0_dp
                  colind = 0
                  nindex = 0
               END IF

               CALL mp_group%sum(nindex)
               CALL mp_group%sum(colind(1:nindex))
               CALL mp_group%sum(erival(1:nindex))
               CALL mp_group%sync()

               DO i34l = 1, nindex
                  i34 = colind(i34l)
                  erint = erival(i34l)
                  CALL csr_idx_from_combined(i34, norb, i3, i4)

                  DO m3 = 1, nmo
                     IF (active_orbitals(m3, spin2) == i3) THEN
                        EXIT
                     END IF
                  END DO

                  DO m4 = 1, nmo
                     IF (active_orbitals(m4, spin2) == i4) THEN
                        EXIT
                     END IF
                  END DO

                  ! terminate the loop prematurely if the function returns false
                  IF (.NOT. fobj%func(m1, m2, m3, m4, erint)) RETURN
               END DO

            END DO
         END DO
      END ASSOCIATE
   END SUBROUTINE eri_type_eri_foreach

END MODULE qs_active_space_types
