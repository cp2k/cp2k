!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2017 created [Jan Wilhelm]
! **************************************************************************************************
MODULE bse
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE group_dist_types,                ONLY: get_group_dist,&
                                              group_dist_d1_type
   USE input_constants,                 ONLY: bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse'

   PUBLIC :: mult_B_with_W, fill_local_3c_arrays, do_subspace_iterations, deallocate_matrices_BSE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param B_iaQ_bse_local ...
!> \param homo ...
!> \param virtual ...
!> \param bse_spin_config ...
!> \param unit_nr ...
!> \param Eigenval ...
!> \param para_env ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo, virtual, bse_spin_config, unit_nr, &
                                     Eigenval, para_env, mp2_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                            B_bar_iaQ_bse_local, B_iaQ_bse_local
      INTEGER                                            :: homo, virtual, bse_spin_config, unit_nr
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(mp2_type)                                     :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'do_subspace_iterations'

      CHARACTER(LEN=10)                                  :: bse_davidson_abort_cond_string, &
                                                            success_abort_string
      CHARACTER(LEN=100)                                 :: char_a, en_exact_filename, en_filename, &
                                                            res_filename
      INTEGER :: bse_davidson_abort_cond, davidson_converged, fac_max_z_space, handle, i_iter, &
         j_print, local_RI_size, num_add_start_z_space, num_davidson_iter, &
         num_en_unconverged = -1, num_exact_en_unconverged = -1, num_exc_en, num_max_z_space, &
         num_new_t, num_res_unconverged = -1, num_Z_vectors, num_Z_vectors_init
      LOGICAL                                            :: bse_full_diag_debug
      REAL(kind=dp)                                      :: eps_exc_en, eps_res, max_en_diff, &
                                                            max_res_norm, z_space_energy_cutoff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: En_diffs, En_diffs_exact, Full_exc_spectrum, &
         Res_norms, Subspace_full_eigenval, Subspace_new_eigenval, Subspace_prev_eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: AZ_reshaped, M_ia_tmp, M_ji_tmp, RI_vector, &
         Subspace_new_eigenvec, Subspace_residuals_reshaped, Z_vectors_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, BZ, Subspace_add_dir, &
                                                            Subspace_ritzvec, W_vectors, Z_vectors

!MG: debugging

      CALL timeset(routineN, handle)

      !MG to del
      !Debug flag for exact diagonalization (only using lapack!!!)
      bse_full_diag_debug = .TRUE.

      bse_davidson_abort_cond = mp2_env%ri_g0w0%davidson_abort_cond
      num_exc_en = mp2_env%ri_g0w0%num_exc_en
      num_add_start_z_space = mp2_env%ri_g0w0%num_add_start_z_space
      fac_max_z_space = mp2_env%ri_g0w0%fac_max_z_space
      num_new_t = mp2_env%ri_g0w0%num_new_t
      num_davidson_iter = mp2_env%ri_g0w0%num_davidson_iter
      eps_res = mp2_env%ri_g0w0%eps_res
      eps_exc_en = mp2_env%ri_g0w0%eps_exc_en
      z_space_energy_cutoff = mp2_env%ri_g0w0%z_space_energy_cutoff

      num_Z_vectors_init = num_exc_en + num_add_start_z_space

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) "bse_spin_config", bse_spin_config
         WRITE (unit_nr, *) "num_exc_en", num_exc_en
         WRITE (unit_nr, *) "num_add_start_z_space", num_add_start_z_space
         WRITE (unit_nr, *) "num_Z_vectors_init", num_Z_vectors_init
         WRITE (unit_nr, *) "fac_max_z_space", fac_max_z_space
         WRITE (unit_nr, *) "num_new_t", num_new_t
         WRITE (unit_nr, *) "eps_res", eps_res
         WRITE (unit_nr, *) "num_davidson_iter", num_davidson_iter
         WRITE (unit_nr, *) "eps_exc_en", eps_exc_en
         WRITE (unit_nr, *) "bse_davidson_abort_cond", bse_davidson_abort_cond
         WRITE (unit_nr, *) "z_space_energy_cutoff", z_space_energy_cutoff
      END IF

      !MG to del
      ! print *,"Inside bse.F, before starting iterations"
      ! print *, "Printing B_bar_iaQ_bse_local of shape", SHAPE(B_bar_iaQ_bse_local)
      ! do i_print=1,SIZE(B_bar_iaQ_bse_local,3)
      !   print *, B_bar_iaQ_bse_local(:,:,i_print)
      ! end do
      ! print *, "Printing B_iaQ_bse_local of shape", SHAPE(B_iaQ_bse_local)
      ! do i_print=1,SIZE(B_iaQ_bse_local,3)
      !   print *, B_iaQ_bse_local(:,:,i_print)
      ! end do
      ! print *, "Printing B_abQ_bse_local of shape", SHAPE(B_abQ_bse_local)
      ! do k_print=1,SIZE(B_abQ_bse_local,3)
      !    do j_print=1,SIZE(B_abQ_bse_local,1)
      !       do i_print=1,SIZE(B_abQ_bse_local,2)
      !          print *,"(",i_print,j_print,k_print,")",B_abQ_bse_local(i_print,j_print,k_print)
      !          print *,"(",j_print,i_print,k_print,")",B_abQ_bse_local(j_print,i_print,k_print)
      !       end do
      !    end do
      ! end do
      ! print *, "Printing B_bar_ijQ_bse_local of shape", SHAPE(B_bar_ijQ_bse_local)
      ! do k_print=1,SIZE(B_bar_ijQ_bse_local,3)
      !    do j_print=1,SIZE(B_bar_ijQ_bse_local,1)
      !      do i_print=1,SIZE(B_bar_ijQ_bse_local,2)
      !         print *,"(",i_print,j_print,k_print,")",B_bar_ijQ_bse_local(i_print,j_print,k_print)
      !         print *,"(",j_print,i_print,k_print,")",B_bar_ijQ_bse_local(j_print,i_print,k_print)
      !     end do
      !   end do
      ! end do

      local_RI_size = SIZE(B_iaQ_bse_local, 3)

      num_Z_vectors = num_Z_vectors_init
      num_max_z_space = num_Z_vectors_init*fac_max_z_space

      !Check input parameters and correct them if necessary
      IF (num_new_t > num_Z_vectors_init) THEN
         num_new_t = num_Z_vectors_init
         IF (unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, "Number of added directions has to be smaller/equals than "// &
                         "initial dimension. Corrected num_new_t accordingly.")
         END IF
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) "Between BSE correction Warnings"
      END IF
      !If initial number is too big, already the first iteration causes trouble in LAPACK diagonal. (DORGQR)
      IF (2*num_Z_vectors_init > homo*virtual) THEN
         CALL cp_abort(__LOCATION__, "Initial dimension was too large and could not be corrected. "// &
                       "Choose another num_exc_en and num_add_start_z_space or adapt your basis set.")
      END IF
      IF (num_max_z_space .GE. homo*virtual) THEN
         !WRITE (unit_nr, *) "Incoming fac=",fac_max_z_space," and max=",num_max_z_space
         !WRITE (unit_nr, *) "homo*virtual",homo*virtual," and init=",num_Z_vectors_init
         fac_max_z_space = homo*virtual/num_Z_vectors_init
         num_max_z_space = num_Z_vectors_init*fac_max_z_space
         !WRITE (unit_nr, *) "Outgoing fac=",fac_max_z_space," and max=",num_max_z_space
         IF (fac_max_z_space == 0) THEN
            CALL cp_abort(__LOCATION__, "Maximal dimension was too large and could not be corrected. "// &
                          "Choose another fac_max_z_space and num_Z_vectors_init or adapt your basis set.")
         ELSE
            IF (unit_nr > 0) THEN
               CALL cp_warn(__LOCATION__, "Maximal dimension of Z space has to be smaller than homo*virtual. "// &
                            "Corrected fac_max_z_space accordingly.")
            END IF
         END IF
      END IF

      DO i_iter = 1, num_davidson_iter
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "Allocating Z_vec,AZ,BZ with dimensions (homo,virt,num_Z)", homo, virtual, num_Z_vectors
            WRITE (unit_nr, *) 'ProcNr', para_env%mepos, 'you really enter here for i_iter', i_iter
         END IF
         ALLOCATE (Z_vectors(homo, virtual, num_Z_vectors))
         Z_vectors = 0.0_dp
         !IF (unit_nr > 0) THEN
         !   print *, 'ProcNr',para_env%mepos,"Allocated Z_vec"
         !END IF
         !Dellocation procedures are a bit intricate, W_/Z_vectors and eigenvalues are needed for the next iteration,
         !  therefore we have to deallocate them separately from the other quantities
         IF (i_iter == 1) THEN
            CALL initial_guess_Z_vectors(Z_vectors, Eigenval, num_Z_vectors, homo, virtual)
            ALLOCATE (Subspace_prev_eigenval(num_exc_en))
            Subspace_prev_eigenval = 0.0_dp
         ELSE
            Z_vectors = W_vectors
            DEALLOCATE (W_vectors)
         END IF
         IF (unit_nr > 0) THEN
            PRINT *, 'ProcNr', para_env%mepos, "Allocated/rewritten Z arrays"
         END IF
         !IF (unit_nr > 0) THEN
         !   print *, 'ProcNr',para_env%mepos,"Allocated W_vectors"
         !END IF
         !WRITE (unit_nr, *) "Init - Checkpoint 1 for", i_iter
         !CALL mp_sync(para_env%group)
         CALL create_bse_work_arrays(AZ, Z_vectors_reshaped, AZ_reshaped, BZ, M_ia_tmp, M_ji_tmp, &
                                     RI_vector, Subspace_new_eigenval, Subspace_full_eigenval, Subspace_new_eigenvec, &
                                     Subspace_residuals_reshaped, Subspace_ritzvec, Subspace_add_dir, W_vectors, &
                                     homo, virtual, num_Z_vectors, local_RI_size, num_new_t, num_exc_en, &
                                     para_env)
         IF (unit_nr > 0) THEN
            PRINT *, 'ProcNr', para_env%mepos, "Allocated Work arrays"
         END IF
         !CALL mp_sync(para_env%group)
         !Check allocation of matrices
         ! if (.NOT. allocated(AZ)) print *, 'ProcNr',para_env%mepos,"Allocated AZ"
         ! if (.NOT. allocated(Z_vectors)) print *, 'ProcNr',para_env%mepos,"Allocated Z_vectors"
         ! if (.NOT. allocated(B_iaQ_bse_local)) print *, 'ProcNr',para_env%mepos,"Allocated B_iaQ_bse_local"
         ! if (.NOT. allocated(B_bar_ijQ_bse_local)) print *, 'ProcNr',para_env%mepos,"Allocated B_bar_ijQ_bse_local"
         ! if (.NOT. allocated(B_abQ_bse_local)) print *, 'ProcNr',para_env%mepos,"Allocated B_abQ_bse_local"
         ! if (.NOT. allocated(M_ia_tmp)) print *, 'ProcNr',para_env%mepos,"Allocated M_ia_tmp"
         ! if (.NOT. allocated(RI_vector)) print *, 'ProcNr',para_env%mepos,"Allocated RI_vector"
         ! CALL mp_sync(para_env%group)
         !print *, 'ProcNr',para_env%mepos,"compute AZ with",SHAPE(AZ), SHAPE(Z_vectors), SHAPE(B_iaQ_bse_local), &
         !SHAPE(B_bar_ijQ_bse_local), SHAPE(B_abQ_bse_local), &
         !SHAPE(M_ia_tmp), SHAPE(RI_vector), SHAPE(Eigenval), homo, virtual, num_Z_vectors, local_RI_size, bse_spin_config
         !CALL mp_sync(para_env%group)
         !print *, 'ProcNr',para_env%mepos,"Entering AZ"
         CALL compute_AZ(AZ, Z_vectors, B_iaQ_bse_local, B_bar_ijQ_bse_local, B_abQ_bse_local, &
                         M_ia_tmp, RI_vector, Eigenval, homo, virtual, num_Z_vectors, local_RI_size, &
                         para_env, bse_spin_config, z_space_energy_cutoff, i_iter, bse_full_diag_debug, &
                         Full_exc_spectrum, unit_nr)

         !MG: functionality of BZ not checked (issue with fm_mat_Q_static_bse_gemm in rpa_util needs to be checked!)
         !CALL compute_BZ(BZ, Z_vectors, B_iaQ_bse_local, B_bar_iaQ_bse_local, &
         !                M_ji_tmp, RI_vector, homo, virtual, num_Z_vectors, local_RI_size, &
         !                para_env)

         !WRITE (unit_nr, *) 'ProcNr',para_env%mepos,'finished AZ/BZ computation'!, AZ=',AZ
         IF (unit_nr > 0) THEN
            PRINT *, 'ProcNr', para_env%mepos, "Computed AZ"
         END IF

         !MG to check: Reshaping correct?
         AZ_reshaped = RESHAPE(AZ, (/homo*virtual, num_Z_vectors/))
         Z_vectors_reshaped = RESHAPE(Z_vectors, (/homo*virtual, num_Z_vectors/))

         ! Diagonalize M and extract smallest eigenvalues/corresponding eigenvector
         CALL compute_diagonalize_ZAZ(AZ_reshaped, Z_vectors_reshaped, num_Z_vectors, homo, virtual, Subspace_new_eigenval, &
                                      Subspace_new_eigenvec, num_new_t, Subspace_full_eigenval, para_env, unit_nr)
         IF (unit_nr > 0) THEN
            PRINT *, "Eigenval (eV) in iter=", i_iter, " is:", Subspace_new_eigenval(:6)*evolt
            !print *,"Eigenvec",Subspace_new_eigenvec
         END IF

         ! Threshold in energies
         CALL check_en_convergence(Subspace_full_eigenval, Subspace_prev_eigenval, eps_exc_en, num_en_unconverged, &
                                   num_exc_en, max_en_diff, En_diffs)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "MG: Largest change of desired exc ens =", max_en_diff
         END IF
         ! Compute residuals
         CALL compute_residuals(AZ_reshaped, Z_vectors_reshaped, Subspace_new_eigenval, Subspace_new_eigenvec, &
                                Subspace_residuals_reshaped, homo, virtual, num_new_t, num_Z_vectors, Subspace_ritzvec)

         ! print *,"Residuals are"
         ! DO k_print=1,homo*virtual
         !    print *,"k_print=",k_print
         !    print *,Subspace_residuals_reshaped(k_print,:)
         ! END DO

         !Abort, if residuals are small enough w.r.t threshold
         CALL check_res_convergence(Subspace_residuals_reshaped, num_new_t, eps_res, num_res_unconverged, &
                                    i_iter, max_res_norm, unit_nr, Res_norms)

         davidson_converged = -1
         IF (num_res_unconverged == 0 .AND. bse_davidson_abort_cond /= 0) THEN
            davidson_converged = 1
            success_abort_string = "RESIDUALS"
         ELSE IF (num_en_unconverged == 0 .AND. (bse_davidson_abort_cond /= 1)) THEN
            davidson_converged = 1
            success_abort_string = "ENERGIES"
         ELSE IF (i_iter == num_davidson_iter) THEN
            davidson_converged = -100
            success_abort_string = "-----"
         ELSE
            davidson_converged = -1
         END IF

         IF (bse_davidson_abort_cond == 0) THEN
            bse_davidson_abort_cond_string = "ENERGY"
         ELSE IF (bse_davidson_abort_cond == 1) THEN
            bse_davidson_abort_cond_string = "RESIDUAL"
         ELSE
            bse_davidson_abort_cond_string = "EITHER"
         END IF

         IF (davidson_converged == 1) THEN
            CALL success_message(Subspace_full_eigenval, homo, virtual, num_new_t, eps_res, num_res_unconverged, &
                                 bse_spin_config, unit_nr, num_exc_en, num_add_start_z_space, num_Z_vectors_init, &
                                 fac_max_z_space, num_davidson_iter, i_iter, num_Z_vectors, num_max_z_space, max_res_norm, &
                                 davidson_converged, max_en_diff, num_en_unconverged, bse_davidson_abort_cond_string, &
                                 eps_exc_en, success_abort_string, z_space_energy_cutoff)

            !Deallocate matrices, which are otherwise not cleared due to exiting the loop
            DEALLOCATE (AZ, & !BZ,
                        Z_vectors, M_ia_tmp, M_ji_tmp, RI_vector, Subspace_prev_eigenval, &
                        Subspace_new_eigenval, Subspace_new_eigenvec, Subspace_residuals_reshaped, &
                        Subspace_add_dir, AZ_reshaped, Z_vectors_reshaped, Subspace_ritzvec, Subspace_full_eigenval)

            EXIT
         ELSE IF (davidson_converged < -1) THEN
            CALL print_davidson_parameter(i_iter, num_davidson_iter, num_Z_vectors, num_res_unconverged, max_res_norm, &
                                          eps_res, num_en_unconverged, max_en_diff, eps_exc_en, num_exc_en, &
                                          num_Z_vectors_init, num_max_z_space, num_new_t, unit_nr, &
                                          success_abort_string, bse_davidson_abort_cond_string, z_space_energy_cutoff)

            CALL cp_abort(__LOCATION__, "BSE/TDA-Davidson did not converge using "// &
                          bse_davidson_abort_cond_string//" threshold condition!")
         END IF

         ! Calculate and add next orthonormal vector and update num_Z_vectors
         CALL compute_new_directions(homo, virtual, Subspace_residuals_reshaped, Subspace_new_eigenval, Eigenval, &
                                     num_new_t, Subspace_add_dir)

         ! print *,"Subspace_add_dir"
         ! print *, Subspace_add_dir
         ! print *,"Z_vectors"
         ! DO k_print=1,num_Z_vectors
         !    print *,k_print
         !    print *,Z_vectors(:,:,k_print)
         ! END DO

         !If exact-diag: compute difference to exact eigenvalues
         IF (bse_full_diag_debug) THEN
            ALLOCATE (En_diffs_exact(num_exc_en))
            num_exact_en_unconverged = 0
            DO j_print = 1, num_exc_en
               En_diffs_exact(j_print) = ABS(Subspace_full_eigenval(j_print) - Full_exc_spectrum(j_print))
               IF (En_diffs_exact(j_print) > eps_exc_en) num_exact_en_unconverged = num_exact_en_unconverged + 1
            END DO
         END IF

         ! res_filename = "residuals_init="//trim(char_a)//"_new="//trim(char_b)//"_max="//trim(char_c)//".dat"

         !Declare formatters
         WRITE (char_a, "(A9,I2,A6)") "(I5, I5, ", num_exc_en, "E11.4)"

         WRITE (res_filename, '(A15, I4, A5, I4, A5, I4, A4)') &
            "residuals_init=", num_Z_vectors_init, "_new=", num_new_t, "_max=", num_max_z_space, ".dat"
         WRITE (en_filename, '(A12, I4, A5, I4, A5, I4, A4)') &
            "energy_init=", num_Z_vectors_init, "_new=", num_new_t, "_max=", num_max_z_space, ".dat"
         WRITE (en_exact_filename, '(A18, I4, A5, I4, A5, I4, A4)') &
            "energy_exact_init=", num_Z_vectors_init, "_new=", num_new_t, "_max=", num_max_z_space, ".dat"

         ! IF (unit_nr > 0) THEN
         !    IF (i_iter == 1) THEN
         !       OPEN(10, file=res_filename, status="new")
         !       OPEN(11, file=en_filename, status="new")
         !    ELSE
         !       OPEN(10, file=res_filename, status="old", access="sequential", position="append")
         !       OPEN(11, file=en_filename, status="old", access="sequential", position="append")
         !    END IF

         !    WRITE(10,'(I4, I3, E11.4)') i_iter, num_res_unconverged, max_res_norm
         !    WRITE(11,char_a) i_iter, num_en_unconverged, En_diffs
         !    CLOSE(10)
         !    CLOSE(11)
         ! END IF

         ! IF (bse_full_diag_debug) THEN
         !    IF (unit_nr > 0) THEN
         !       IF (i_iter == 1) THEN
         !          OPEN(10, file=en_exact_filename, status="new")
         !       ELSE
         !          OPEN(10, file=en_exact_filename, status="old", access="sequential", position="append")
         !       END IF
         !       WRITE(10,char_a) i_iter, num_exact_en_unconverged ,En_diffs_exact
         !       CLOSE(10)
         !    END IF
         ! END IF

         !Check dimensions and orthonormalize vector system, depending on dimensionality
         CALL check_Z_space_dimension(W_vectors, Z_vectors, Subspace_add_dir, Subspace_ritzvec, &
                               num_Z_vectors, num_new_t, num_max_z_space, homo, virtual, i_iter, unit_nr, res_filename, en_filename)

         ! print *,"W_vectors"
         ! DO k_print=1,num_Z_vectors
         !    print *,k_print
         !    print *,W_vectors(:,:,k_print)
         ! END DO

         !Copy eigenvalues for threshold
         Subspace_prev_eigenval(:) = Subspace_full_eigenval(:num_exc_en)

         DEALLOCATE (AZ, & !BZ,
                     Z_vectors, M_ia_tmp, M_ji_tmp, RI_vector, &
                     Subspace_new_eigenval, Subspace_new_eigenvec, Subspace_residuals_reshaped, &
                     Subspace_add_dir, AZ_reshaped, Z_vectors_reshaped, Subspace_ritzvec, Subspace_full_eigenval, &
                     Res_norms, En_diffs)

         IF (bse_full_diag_debug) THEN
            DEALLOCATE (En_diffs_exact)
         END IF

         !Orthonorm:
         CALL orthonormalize_W(W_vectors, num_Z_vectors, homo, virtual)

         ! print *,"W_vectors"
         ! DO k_print=1,num_Z_vectors
         !    print *,k_print
         !    print *,W_vectors(:,:,k_print)
         ! END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param W_vectors ...
!> \param Z_vectors ...
!> \param Subspace_add_dir ...
!> \param Subspace_ritzvec ...
!> \param num_Z_vectors ...
!> \param num_new_t ...
!> \param num_max_z_space ...
!> \param homo ...
!> \param virtual ...
!> \param i_iter ...
!> \param unit_nr ...
!> \param res_filename ...
!> \param en_filename ...
! **************************************************************************************************
   SUBROUTINE check_Z_space_dimension(W_vectors, Z_vectors, Subspace_add_dir, Subspace_ritzvec, &
                                      num_Z_vectors, num_new_t, num_max_z_space, homo, virtual, i_iter, unit_nr, &
                                      res_filename, en_filename)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: W_vectors, Z_vectors, Subspace_add_dir, &
                                                            Subspace_ritzvec
      INTEGER                                            :: num_Z_vectors, num_new_t, &
                                                            num_max_z_space, homo, virtual, &
                                                            i_iter, unit_nr
      CHARACTER(LEN=100)                                 :: res_filename, en_filename

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_Z_space_dimension'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (num_Z_vectors + num_new_t .LE. num_max_z_space) THEN
         W_vectors(:, :, :num_Z_vectors) = Z_vectors(:, :, :)
         W_vectors(:, :, num_Z_vectors + 1:) = Subspace_add_dir
         num_Z_vectors = num_Z_vectors + num_new_t
      ELSE
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "MG: Resetting dimension in i_iter=", i_iter
         END IF
         ! IF (unit_nr > 0) THEN
         !    OPEN(10, file=res_filename, status="old", access="sequential", position="append")
         !       WRITE(10,*) "Resetting", i_iter
         !    CLOSE(10)
         !    OPEN(10, file=en_filename, status="old", access="sequential", position="append")
         !       WRITE(10,*) "Resetting", i_iter
         !    CLOSE(10)
         ! END IF
         DEALLOCATE (W_vectors)
         ALLOCATE (W_vectors(homo, virtual, 2*num_new_t))
         W_vectors(:, :, :num_new_t) = Subspace_ritzvec(:, :, :)
         W_vectors(:, :, num_new_t + 1:) = Subspace_add_dir
         num_Z_vectors = 2*num_new_t
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors_reshaped ...
!> \param AZ_reshaped ...
!> \param BZ ...
!> \param M_ia_tmp ...
!> \param M_ji_tmp ...
!> \param RI_vector ...
!> \param Subspace_new_eigenval ...
!> \param Subspace_full_eigenval ...
!> \param Subspace_new_eigenvec ...
!> \param Subspace_residuals_reshaped ...
!> \param Subspace_ritzvec ...
!> \param Subspace_add_dir ...
!> \param W_vectors ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param num_new_t ...
!> \param num_exc_en ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE create_bse_work_arrays(AZ, Z_vectors_reshaped, AZ_reshaped, BZ, M_ia_tmp, M_ji_tmp, &
                                     RI_vector, Subspace_new_eigenval, Subspace_full_eigenval, Subspace_new_eigenvec, &
                                     Subspace_residuals_reshaped, Subspace_ritzvec, Subspace_add_dir, W_vectors, &
                                     homo, virtual, num_Z_vectors, local_RI_size, num_new_t, num_exc_en, &
                                     para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Z_vectors_reshaped, AZ_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BZ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_ia_tmp, M_ji_tmp, RI_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_new_eigenval, &
                                                            Subspace_full_eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_new_eigenvec, &
                                                            Subspace_residuals_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Subspace_ritzvec, Subspace_add_dir, &
                                                            W_vectors
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size, num_new_t, num_exc_en
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_bse_work_arrays'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      !print *,"Procnr", para_env%mepos, "local_RI_size=",local_RI_size

      ALLOCATE (AZ(homo, virtual, num_Z_vectors))
      AZ = 0.0_dp

      ALLOCATE (Z_vectors_reshaped(homo*virtual, num_Z_vectors))
      Z_vectors_reshaped = 0.0_dp

      ALLOCATE (AZ_reshaped(homo*virtual, num_Z_vectors))
      AZ_reshaped = 0.0_dp

      !ALLOCATE (BZ(homo, virtual, num_Z_vectors))
      !BZ = 0.0_dp

      ALLOCATE (M_ia_tmp(homo, virtual))
      M_ia_tmp = 0.0_dp

      ALLOCATE (M_ji_tmp(homo, homo))
      M_ji_tmp = 0.0_dp

      ALLOCATE (RI_vector(local_RI_size, num_Z_vectors))
      RI_vector = 0.0_dp

      ALLOCATE (Subspace_new_eigenval(num_new_t))
      Subspace_new_eigenval = 0.0_dp

      ALLOCATE (Subspace_full_eigenval(num_Z_vectors))
      Subspace_full_eigenval = 0.0_dp

      ALLOCATE (Subspace_new_eigenvec(num_Z_vectors, num_new_t))
      Subspace_new_eigenvec = 0.0_dp

      ALLOCATE (subspace_residuals_reshaped(homo*virtual, num_new_t))
      subspace_residuals_reshaped = 0.0_dp

      ALLOCATE (Subspace_ritzvec(homo, virtual, num_new_t))
      Subspace_ritzvec = 0.0_dp

      ALLOCATE (Subspace_add_dir(homo, virtual, num_new_t))
      Subspace_add_dir = 0.0_dp

      ALLOCATE (W_vectors(homo, virtual, num_Z_vectors + num_new_t))
      W_vectors = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Subspace_full_eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param num_new_t ...
!> \param eps_res ...
!> \param num_res_unconverged ...
!> \param bse_spin_config ...
!> \param unit_nr ...
!> \param num_exc_en ...
!> \param num_add_start_z_space ...
!> \param num_Z_vectors_init ...
!> \param fac_max_z_space ...
!> \param num_davidson_iter ...
!> \param i_iter ...
!> \param num_Z_vectors ...
!> \param num_max_z_space ...
!> \param max_res_norm ...
!> \param davidson_converged ...
!> \param max_en_diff ...
!> \param num_en_unconverged ...
!> \param bse_davidson_abort_cond_string ...
!> \param eps_exc_en ...
!> \param success_abort_string ...
!> \param z_space_energy_cutoff ...
! **************************************************************************************************
   SUBROUTINE success_message(Subspace_full_eigenval, homo, virtual, num_new_t, eps_res, num_res_unconverged, &
                              bse_spin_config, unit_nr, num_exc_en, num_add_start_z_space, num_Z_vectors_init, &
                              fac_max_z_space, num_davidson_iter, i_iter, num_Z_vectors, num_max_z_space, max_res_norm, &
                              davidson_converged, max_en_diff, num_en_unconverged, bse_davidson_abort_cond_string, &
                              eps_exc_en, success_abort_string, z_space_energy_cutoff)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_full_eigenval
      INTEGER                                            :: homo, virtual, num_new_t
      REAL(KIND=dp)                                      :: eps_res
      INTEGER :: num_res_unconverged, bse_spin_config, unit_nr, num_exc_en, num_add_start_z_space, &
         num_Z_vectors_init, fac_max_z_space, num_davidson_iter, i_iter, num_Z_vectors, &
         num_max_z_space
      REAL(KIND=dp)                                      :: max_res_norm
      INTEGER                                            :: davidson_converged
      REAL(KIND=dp)                                      :: max_en_diff
      INTEGER                                            :: num_en_unconverged
      CHARACTER(LEN=10)                                  :: bse_davidson_abort_cond_string
      REAL(KIND=dp)                                      :: eps_exc_en
      CHARACTER(LEN=10)                                  :: success_abort_string
      REAL(KIND=dp)                                      :: z_space_energy_cutoff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'success_message'

      CHARACTER(LEN=10)                                  :: multiplet
      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: alpha

      CALL timeset(routineN, handle)

      !Prepare variables for printing
      SELECT CASE (bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
         multiplet = "Singlet"
      CASE (bse_triplet)
         alpha = 0.0_dp
         multiplet = "Triplet"
      END SELECT
      !Print statement
      !print *,"unit_nr=",unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T3,A)') '******************************************************************************'
         WRITE (unit_nr, '(T3,A)') '**                                                                          **'
         WRITE (unit_nr, '(T3,A)') '**                        BSE-TDA EXCITONIC ENERGIES                        **'
         WRITE (unit_nr, '(T3,A)') '**                                                                          **'
         WRITE (unit_nr, '(T3,A)') '******************************************************************************'
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') ' The excitation energies are calculated by iteratively diagonalizing: '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') '    A_iajb   =  (E_a-E_i) delta_ij delta_ab   +  alpha * v_iajb   -  W_ijab   '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A48,A7,A12,F3.1)') &
            ' The spin-dependent factor for the requested ', multiplet, " is alpha = ", alpha
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A16,T50,A22)') &
            ' Excitonic level', 'Excitation energy (eV)'
         !prints actual energies values
         DO i = 1, num_exc_en
            WRITE (unit_nr, '(T3,I16,T50,F22.3)') i, Subspace_full_eigenval(i)*evolt
         END DO

         WRITE (unit_nr, '(T3,A)') ' '

         !prints parameters of Davidson algorithm
         CALL print_davidson_parameter(i_iter, num_davidson_iter, num_Z_vectors, num_res_unconverged, max_res_norm, &
                                       eps_res, num_en_unconverged, max_en_diff, eps_exc_en, num_exc_en, &
                                       num_Z_vectors_init, num_max_z_space, num_new_t, unit_nr, &
                                       success_abort_string, bse_davidson_abort_cond_string, z_space_energy_cutoff)

         !Insert warning if energies are not converged (could probably be the case if one uses residual threshold)
         IF (num_en_unconverged > 0) THEN
            !CPABORT("UNCONVERED ENERGIES ARE PRINTED OUT - THRESHOLDS SHOULD BE CHOSEN MORE CAREFULLY!")
            WRITE (unit_nr, '(T3,A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE (unit_nr, '(T3,A2,T79,A2)') '!!', "!!"
            WRITE (unit_nr, '(T3,A2,T8,A65,T79,A2)') '!!', "THERE ARE UNCONVERGED ENERGIES PRINTED OUT, SOMETHING WENT WRONG!", "!!"
            WRITE (unit_nr, '(T3,A2,T79,A2)') '!!', "!!"
            WRITE (unit_nr, '(T3,A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param i_iter ...
!> \param num_davidson_iter ...
!> \param num_Z_vectors ...
!> \param num_res_unconverged ...
!> \param max_res_norm ...
!> \param eps_res ...
!> \param num_en_unconverged ...
!> \param max_en_diff ...
!> \param eps_exc_en ...
!> \param num_exc_en ...
!> \param num_Z_vectors_init ...
!> \param num_max_z_space ...
!> \param num_new_t ...
!> \param unit_nr ...
!> \param success_abort_string ...
!> \param bse_davidson_abort_cond_string ...
!> \param z_space_energy_cutoff ...
! **************************************************************************************************
   SUBROUTINE print_davidson_parameter(i_iter, num_davidson_iter, num_Z_vectors, num_res_unconverged, max_res_norm, &
                                       eps_res, num_en_unconverged, max_en_diff, eps_exc_en, num_exc_en, &
                                       num_Z_vectors_init, num_max_z_space, num_new_t, unit_nr, &
                                       success_abort_string, bse_davidson_abort_cond_string, z_space_energy_cutoff)

      INTEGER                                            :: i_iter, num_davidson_iter, &
                                                            num_Z_vectors, num_res_unconverged
      REAL(KIND=dp)                                      :: max_res_norm, eps_res
      INTEGER                                            :: num_en_unconverged
      REAL(KIND=dp)                                      :: max_en_diff, eps_exc_en
      INTEGER                                            :: num_exc_en, num_Z_vectors_init, &
                                                            num_max_z_space, num_new_t, unit_nr
      CHARACTER(LEN=10)                                  :: success_abort_string, &
                                                            bse_davidson_abort_cond_string
      REAL(KIND=dp)                                      :: z_space_energy_cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_davidson_parameter'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      WRITE (unit_nr, '(T3,A)') '******************************************************************************'
      WRITE (unit_nr, '(T3,A2,T15,A49,T79,A2)') &
         '**', "Parameters of the BSE-Davidson solver:", "**"
      WRITE (unit_nr, '(T3,A2,T79,A2)') &
         '**', "**"
      WRITE (unit_nr, '(T3,A2,T79,A2)') &
         '**', "**"
      WRITE (unit_nr, '(T3,A2,T10,A16,I5,A12,I5,A8,T79,A2)') &
         '**', "Converged after ", i_iter, " of maximal ", num_davidson_iter, " cycles,", "**"
      WRITE (unit_nr, '(T3,A2,T20,A11,A9,A7,A8,A20,T79,A2)') &
         '**', "because of ", success_abort_string, " using ", &
         bse_davidson_abort_cond_string, " threshold condition", "**"
      WRITE (unit_nr, '(T3,A2,T79,A2)') &
         '**', "**"
      WRITE (unit_nr, '(T3,A2,T10,A32,T65,I11,T79,A2)') &
         '**', "The Z space has at the end dim. ", num_Z_vectors, "**"
      WRITE (unit_nr, '(T3,A2,T10,A45,T65,I11,T79,A2)') &
         '**', "Number of unconverged residuals in subspace: ", num_res_unconverged, "**"
      WRITE (unit_nr, '(T3,A2,T10,A35,T65,E11.4,T79,A2)') &
         '**', "largest unconverged residual (eV): ", max_res_norm*evolt, "**"
      WRITE (unit_nr, '(T3,A2,T10,A45,T65,E11.4,T79,A2)') &
         '**', "threshold for convergence of residuals (eV): ", eps_res*evolt, "**"
      WRITE (unit_nr, '(T3,A2,T10,A45,T65,I11,T79,A2)') &
         '**', "Number of desired, but unconverged energies: ", num_en_unconverged, "**"
      WRITE (unit_nr, '(T3,A2,T10,A44,T65,E11.4,T79,A2)') &
         '**', "largest unconverged energy difference (eV): ", max_en_diff*evolt, "**"
      WRITE (unit_nr, '(T3,A2,T10,A44,T65,E11.4,T79,A2)') &
         '**', "threshold for convergence of energies (eV): ", eps_exc_en*evolt, "**"
      WRITE (unit_nr, '(T3,A2,T10,A40,T65,I11,T79,A2)') &
         '**', "number of computed excitation energies: ", num_exc_en, "**"

      IF (z_space_energy_cutoff > 0) THEN
         WRITE (unit_nr, '(T3,A2,T10,A37,T65,E11.4,T79,A2)') &
            '**', "cutoff for excitation energies (eV): ", z_space_energy_cutoff*evolt, "**"
      END IF

      WRITE (unit_nr, '(T3,A2,T10,A36,T65,I11,T79,A2)') &
         '**', "number of Z space at the beginning: ", num_Z_vectors_init, "**"
      WRITE (unit_nr, '(T3,A2,T10,A30,T65,I11,T79,A2)') &
         '**', "maximal dimension of Z space: ", num_max_z_space, "**"
      WRITE (unit_nr, '(T3,A2,T10,A31,T65,I11,T79,A2)') &
         '**', "added directions per iteration: ", num_new_t, "**"
      WRITE (unit_nr, '(T3,A2,T79,A2)') &
         '**', "**"
      WRITE (unit_nr, '(T3,A2,T79,A2)') &
         '**', "**"
      WRITE (unit_nr, '(T3,A)') '******************************************************************************'
      WRITE (unit_nr, '(T3,A)') ' '

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Subspace_full_eigenval ...
!> \param Subspace_prev_eigenval ...
!> \param eps_exc_en ...
!> \param num_en_unconverged ...
!> \param num_exc_en ...
!> \param max_en_diff ...
!> \param En_diffs ...
! **************************************************************************************************
   SUBROUTINE check_en_convergence(Subspace_full_eigenval, Subspace_prev_eigenval, eps_exc_en, num_en_unconverged, &
                                   num_exc_en, max_en_diff, En_diffs)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_full_eigenval, &
                                                            Subspace_prev_eigenval
      REAL(KIND=dp)                                      :: eps_exc_en
      INTEGER                                            :: num_en_unconverged, num_exc_en
      REAL(KIND=dp)                                      :: max_en_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: En_diffs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_en_convergence'

      INTEGER                                            :: handle, mu_l

      CALL timeset(routineN, handle)

      num_en_unconverged = 0
      ALLOCATE (En_diffs(num_exc_en))
      DO mu_l = 1, num_exc_en
         En_diffs(mu_l) = ABS(Subspace_full_eigenval(mu_l) - Subspace_prev_eigenval(mu_l))
         IF (En_diffs(mu_l) > eps_exc_en) num_en_unconverged = num_en_unconverged + 1
      END DO
      max_en_diff = MAXVAL(En_diffs)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Subspace_residuals_reshaped ...
!> \param num_new_t ...
!> \param eps_res ...
!> \param num_res_unconverged ...
!> \param i_iter ...
!> \param max_res_norm ...
!> \param unit_nr ...
!> \param Res_norms ...
! **************************************************************************************************
   SUBROUTINE check_res_convergence(Subspace_residuals_reshaped, num_new_t, eps_res, num_res_unconverged, &
                                    i_iter, max_res_norm, unit_nr, Res_norms)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_residuals_reshaped
      INTEGER                                            :: num_new_t
      REAL(KIND=dp)                                      :: eps_res
      INTEGER                                            :: num_res_unconverged, i_iter
      REAL(KIND=dp)                                      :: max_res_norm
      INTEGER                                            :: unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Res_norms

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_res_convergence'

      INTEGER                                            :: handle, mu_l

      CALL timeset(routineN, handle)

      num_res_unconverged = 0
      ALLOCATE (Res_norms(num_new_t))
      DO mu_l = 1, num_new_t
         Res_norms(mu_l) = NORM2(Subspace_residuals_reshaped(:, mu_l))
         IF (Res_norms(mu_l) > eps_res) THEN
            num_res_unconverged = num_res_unconverged + 1
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "MG: Unconverged res in i_iter=", i_iter, "is:", Res_norms(mu_l)
            END IF
         END IF
      END DO
      max_res_norm = MAXVAL(Res_norms)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) "MG: Maximal unconverged res (of ", num_res_unconverged, &
            " unconverged res in this step) in i_iter=", i_iter, "is:", max_res_norm
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param W_vectors ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE orthonormalize_W(W_vectors, num_Z_vectors, homo, virtual)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: W_vectors
      INTEGER                                            :: num_Z_vectors, homo, virtual

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'orthonormalize_W'

      INTEGER                                            :: handle, info_dor, info_orth, LWORK_dor, &
                                                            LWORK_W
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Tau_W, WORK_W, WORK_W_dor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: W_vectors_reshaped

      CALL timeset(routineN, handle)

      ALLOCATE (W_vectors_reshaped(homo*virtual, num_Z_vectors))
      W_vectors_reshaped = RESHAPE(W_vectors, (/homo*virtual, num_Z_vectors/))

      ALLOCATE (Tau_W(MIN(homo*virtual, num_Z_vectors)))
      Tau_W = 0.0_dp

      ALLOCATE (WORK_W(1))
      WORK_W = 0.0_dp

      ALLOCATE (WORK_W_dor(1))
      WORK_W_dor = 0.0_dp

      CALL DGEQRF(homo*virtual, num_Z_vectors, W_vectors_reshaped, homo*virtual, Tau_W, WORK_W, -1, info_orth)
      LWORK_W = INT(WORK_W(1))
      DEALLOCATE (WORK_W)
      ALLOCATE (WORK_W(LWORK_W))
      WORK_W = 0.0_dp
      CALL DGEQRF(homo*virtual, num_Z_vectors, W_vectors_reshaped, homo*virtual, Tau_W, WORK_W, LWORK_W, info_orth)
      IF (info_orth /= 0) THEN
         CPABORT("QR Decomp Step 1 doesnt work")
      END IF
      CALL DORGQR(homo*virtual, num_Z_vectors, MIN(homo*virtual, num_Z_vectors), W_vectors_reshaped, homo*virtual, &
                  Tau_W, WORK_W_dor, -1, info_dor)
      LWORK_dor = INT(WORK_W_dor(1))
      DEALLOCATE (WORK_W_dor)
      ALLOCATE (WORK_W_dor(LWORK_dor))
      WORK_W_dor = 0.0_dp
      CALL DORGQR(homo*virtual, num_Z_vectors, MIN(homo*virtual, num_Z_vectors), W_vectors_reshaped, homo*virtual, &
                  Tau_W, WORK_W_dor, LWORK_dor, info_dor)
      IF (info_orth /= 0) THEN
         CPABORT("QR Decomp Step 2 doesnt work")
      END IF

      W_vectors = RESHAPE(W_vectors_reshaped, (/homo, virtual, num_Z_vectors/))

      DEALLOCATE (WORK_W, WORK_W_dor, Tau_W, W_vectors_reshaped)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param homo ...
!> \param virtual ...
!> \param Subspace_residuals_reshaped ...
!> \param Subspace_new_eigenval ...
!> \param Eigenval ...
!> \param num_new_t ...
!> \param Subspace_add_dir ...
! **************************************************************************************************
   SUBROUTINE compute_new_directions(homo, virtual, Subspace_residuals_reshaped, Subspace_new_eigenval, Eigenval, &
                                     num_new_t, Subspace_add_dir)

      INTEGER                                            :: homo, virtual
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_residuals_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_new_eigenval
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: num_new_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Subspace_add_dir

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_new_directions'

      INTEGER                                            :: a_virt, handle, i_occ, mu_subspace, &
                                                            prec_neg
      REAL(KIND=dp)                                      :: prec_scalar
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_add_dir_reshaped

      CALL timeset(routineN, handle)

      ALLOCATE (Subspace_add_dir_reshaped(homo*virtual, num_new_t))

      prec_neg = 0
      DO mu_subspace = 1, num_new_t
         DO i_occ = 1, homo
            DO a_virt = 1, virtual
               !MG to check: Indexorder and range of indices
               prec_scalar = -1/(Subspace_new_eigenval(mu_subspace) - (Eigenval(a_virt + homo) - Eigenval(i_occ)))
               IF (prec_scalar < 0) THEN
                  prec_neg = prec_neg + 1
                  !prec_scalar = - prec_scalar
               END IF
               Subspace_add_dir_reshaped((i_occ - 1)*virtual + a_virt, mu_subspace) = prec_scalar* &
                                                              Subspace_residuals_reshaped((i_occ - 1)*virtual + a_virt, mu_subspace)
            END DO
         END DO
      END DO
      !if (prec_neg > 0) then
      !   print *, "MG: Prec has ",prec_neg," negative entries of overall ",homo*virtual*num_new_t," entries"
      !end if
      !print *,Subspace_add_dir_reshaped

      Subspace_add_dir = RESHAPE(Subspace_add_dir_reshaped, (/homo, virtual, num_new_t/))

      DEALLOCATE (Subspace_add_dir_reshaped)
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ_reshaped ...
!> \param Z_vectors_reshaped ...
!> \param Subspace_new_eigenval ...
!> \param Subspace_new_eigenvec ...
!> \param Subspace_residuals_reshaped ...
!> \param homo ...
!> \param virtual ...
!> \param num_new_t ...
!> \param num_Z_vectors ...
!> \param Subspace_ritzvec ...
! **************************************************************************************************
   SUBROUTINE compute_residuals(AZ_reshaped, Z_vectors_reshaped, Subspace_new_eigenval, Subspace_new_eigenvec, &
                                Subspace_residuals_reshaped, homo, virtual, num_new_t, num_Z_vectors, Subspace_ritzvec)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: AZ_reshaped, Z_vectors_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_new_eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_new_eigenvec, &
                                                            Subspace_residuals_reshaped
      INTEGER                                            :: homo, virtual, num_new_t, num_Z_vectors
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Subspace_ritzvec

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_residuals'

      INTEGER                                            :: handle, mu_subspace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_res_A, Subspace_res_ev

      CALL timeset(routineN, handle)

      ALLOCATE (Subspace_res_ev(homo*virtual, num_new_t))
      Subspace_res_ev = 0.0_dp

      ALLOCATE (Subspace_res_A(homo*virtual, num_new_t))
      Subspace_res_A = 0.0_dp

      !Compute all residuals in one loop, iterating over number of new/added t per iteration
      DO mu_subspace = 1, num_new_t

         CALL DGEMM("N", "N", homo*virtual, 1, num_Z_vectors, 1.0_dp, Z_vectors_reshaped, homo*virtual, &
                    Subspace_new_eigenvec(:, mu_subspace), num_Z_vectors, 0.0_dp, Subspace_res_ev(:, mu_subspace), homo*virtual)

         CALL DGEMM("N", "N", homo*virtual, 1, num_Z_vectors, 1.0_dp, AZ_reshaped, homo*virtual, &
                    Subspace_new_eigenvec(:, mu_subspace), num_Z_vectors, 0.0_dp, Subspace_res_A(:, mu_subspace), homo*virtual)

         Subspace_residuals_reshaped(:, mu_subspace) = Subspace_new_eigenval(mu_subspace)*Subspace_res_ev(:, mu_subspace) &
                                                       - Subspace_res_A(:, mu_subspace)

      END DO
      Subspace_ritzvec = RESHAPE(Subspace_res_ev, (/homo, virtual, num_new_t/))
      DEALLOCATE (Subspace_res_ev, Subspace_res_A)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ_reshaped ...
!> \param Z_vectors_reshaped ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
!> \param Subspace_new_eigenval ...
!> \param Subspace_new_eigenvec ...
!> \param num_new_t ...
!> \param Subspace_full_eigenval ...
!> \param para_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE compute_diagonalize_ZAZ(AZ_reshaped, Z_vectors_reshaped, num_Z_vectors, homo, virtual, Subspace_new_eigenval, &
                                      Subspace_new_eigenvec, num_new_t, Subspace_full_eigenval, para_env, unit_nr)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: AZ_reshaped, Z_vectors_reshaped
      INTEGER, INTENT(in)                                :: num_Z_vectors, homo, virtual
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_new_eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_new_eigenvec
      INTEGER, INTENT(in)                                :: num_new_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_full_eigenval
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(in)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_diagonalize_ZAZ'

      INTEGER                                            :: handle, i_Z_vector, j_Z_vector, LWORK, &
                                                            ZAZ_diag_info
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: WORK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ZAZ

      CALL timeset(routineN, handle)

      ALLOCATE (ZAZ(num_Z_vectors, num_Z_vectors))
      ZAZ(:, :) = 0.0_dp

      !Flatten AZ and Z matrices of a certain j_Z_vector w.r.t. occ and virt indices
      !Multiply for each j_Z_vec and write into matrix of dim (num_Z_vec, num_Z_vec)
      DO i_Z_vector = 1, num_Z_vectors
         DO j_Z_vector = 1, num_Z_vectors
            ZAZ(j_Z_vector, i_Z_vector) = DOT_PRODUCT(Z_vectors_reshaped(:, j_Z_vector), AZ_reshaped(:, i_Z_vector))
         END DO
      END DO
      IF (unit_nr > 0) THEN
         PRINT *, 'ProcNr', para_env%mepos, "Before Diag"
      END IF

      !MG to do: Check for symmetry of ZAZ!
      ALLOCATE (WORK(1))
      WORK = 0.0_dp
      CALL DSYEV("V", "U", num_Z_vectors, ZAZ, num_Z_vectors, Subspace_full_eigenval, WORK, -1, ZAZ_diag_info)
      LWORK = INT(WORK(1))
      DEALLOCATE (WORK)
      ALLOCATE (WORK(LWORK))
      WORK = 0.0_dp
      !MG to check: Usage of symmetric routine okay? (Correct LWORK?)
      CALL DSYEV("V", "U", num_Z_vectors, ZAZ, num_Z_vectors, Subspace_full_eigenval, WORK, LWORK, ZAZ_diag_info)

      IF (ZAZ_diag_info /= 0) THEN
         CPABORT("ZAZ could not be diagonalized successfully.")
      ELSE
         !print *,"Eigenvalues are", Subspace_full_eigenval
         !print *,"Eigenvectors are", ZAZ
         !print *,"WORK(1)=",WORK(1)
         !print *,"ZAZ_diag_info",ZAZ_diag_info
      END IF

      IF (unit_nr > 0) THEN
         PRINT *, 'ProcNr', para_env%mepos, "After Diag"
      END IF

      !print *,"Shape of Subspace_new_eigenval",SHAPE(Subspace_new_eigenval)
      !print *,"Shape of Subspace_full_eigenval",SHAPE(Subspace_full_eigenval)
      !print *, "num_new_t",num_new_t
      !print *, Subspace_full_eigenval(1:num_new_t)
      !print *, Subspace_new_eigenval
      !print *, Subspace_new_eigenval(:)
      !print *, Subspace_new_eigenval(:num_new_t)
      Subspace_new_eigenval = Subspace_full_eigenval(1:num_new_t)
      !print *,"Between eigenvalues and eigenvectors"
      Subspace_new_eigenvec(:, 1:num_new_t) = ZAZ(:, 1:num_new_t)
      !print *,"After eigenvalues and eigenvectors"
      DEALLOCATE (WORK)
      DEALLOCATE (ZAZ)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param BZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param M_ji_tmp ...
!> \param RI_vector ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_BZ(BZ, Z_vectors, B_iaQ_bse_local, B_bar_iaQ_bse_local, &
                         M_ji_tmp, RI_vector, homo, virtual, num_Z_vectors, local_RI_size, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BZ, Z_vectors, B_iaQ_bse_local, &
                                                            B_bar_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: M_ji_tmp, RI_vector
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i_Z_vector, LLL

      BZ(:, :, :) = 0.0_dp

      !CALL compute_v_ia_jb_part(BZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
      !                          num_Z_vectors, homo, virtual)

      DO i_Z_vector = 1, num_Z_vectors

         DO LLL = 1, local_RI_size

            ! M_ji^P = sum_b Z_jb*B_bi^P
            CALL DGEMM("N", "T", homo, homo, virtual, 1.0_dp, Z_vectors(:, :, i_Z_vector), homo, &
                       B_iaQ_bse_local(:, :, LLL), homo, 0.0_dp, M_ji_tmp, homo)
            ! (BZ)_ia = sum_jP M_ij^P*B^bar_ja^P
            CALL DGEMM("T", "N", homo, virtual, homo, 1.0_dp, M_ji_tmp, homo, &
                       B_bar_iaQ_bse_local, homo, 1.0_dp, BZ(:, :, i_Z_vector), homo)

         END DO

      END DO

      ! we make the sum to sum over all RI basis functions
      CALL para_env%sum(BZ)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param M_ia_tmp ...
!> \param RI_vector ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param para_env ...
!> \param bse_spin_config ...
!> \param z_space_energy_cutoff ...
!> \param i_iter ...
!> \param bse_full_diag_debug ...
!> \param Full_exc_spectrum ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE compute_AZ(AZ, Z_vectors, B_iaQ_bse_local, B_bar_ijQ_bse_local, B_abQ_bse_local, M_ia_tmp, &
                         RI_vector, Eigenval, homo, virtual, num_Z_vectors, local_RI_size, &
                         para_env, bse_spin_config, z_space_energy_cutoff, i_iter, bse_full_diag_debug, &
                         Full_exc_spectrum, unit_nr)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, Z_vectors, B_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_abQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: M_ia_tmp, RI_vector
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER                                            :: bse_spin_config
      REAL(KIND=dp)                                      :: z_space_energy_cutoff
      INTEGER                                            :: i_iter
      LOGICAL                                            :: bse_full_diag_debug
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Full_exc_spectrum
      INTEGER                                            :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_AZ'

      INTEGER                                            :: a, a_virt, b, diag_info, handle, i, &
                                                            i_occ, i_Z_vector, j, LLL, LWORK, m, n
      REAL(KIND=dp)                                      :: eigen_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: WORK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: A_full_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: A_full, v_iajb, W_ijab

!MG: i_iter only for Debugging
!MG to del

      CALL timeset(routineN, handle)
      !print *, 'ProcNr',para_env%mepos,"Entering AZ comp"
      AZ(:, :, :) = 0.0_dp
      !print *, 'ProcNr',para_env%mepos,"Comuting V"
      !MG to del
      ! print *,"3center Coulomb - B_iaQ, LLL=",local_RI_size
      ! do i=1,local_RI_size
      !    print *,"B_iaQ in LLL=",i," is: ",B_iaQ_bse_local(:,:,i)
      ! end do

      IF (i_iter == 1 .AND. bse_full_diag_debug) THEN
         ALLOCATE (W_ijab(homo, homo, virtual, virtual))
         ALLOCATE (A_full(homo, virtual, homo, virtual))
         ALLOCATE (A_full_reshaped(homo*virtual, homo*virtual))
         ALLOCATE (Full_exc_spectrum(homo*virtual))
         W_ijab = 0.0_dp
         A_full = 0.0_dp
         A_full_reshaped = 0.0_dp
         Full_exc_spectrum = 0.0_dp
      END IF

      CALL compute_v_ia_jb_part(AZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                num_Z_vectors, homo, virtual, bse_spin_config, v_iajb, bse_full_diag_debug, i_iter, &
                                unit_nr, para_env)

      !print *, 'ProcNr',para_env%mepos,"Computed V"
      IF (unit_nr > 0 .AND. i_iter == 1) THEN
         PRINT *, "MG: printing v_iajb"
         PRINT *, v_iajb
      END IF

      DO i_Z_vector = 1, num_Z_vectors

         ! JW TO DO: OMP PARALLELIZATION
         DO LLL = 1, local_RI_size

            ! M_ja^P = sum_b Z_jb*B_ba^P
            CALL DGEMM("N", "N", homo, virtual, virtual, 1.0_dp, Z_vectors(:, :, i_Z_vector), homo, &
                       B_abQ_bse_local(:, :, LLL), virtual, 0.0_dp, M_ia_tmp, homo)

            ! (AZ)_ia = sum_jP B_bar_ij^P*M_ja^P
            CALL DGEMM("N", "N", homo, virtual, homo, -1.0_dp, B_bar_ijQ_bse_local(:, :, LLL), homo, &
                       M_ia_tmp, homo, 1.0_dp, AZ(:, :, i_Z_vector), homo)

         END DO
      END DO

      !print *, 'ProcNr',para_env%mepos,"Computed W"
      IF (i_iter == 1 .AND. bse_full_diag_debug) THEN
         W_ijab = 0.0_dp
         !MG to del: Print screened 4c integrals for debugging purposes
         DO LLL = 1, local_RI_size
            DO i = 1, homo
               DO j = 1, homo
                  DO a = 1, virtual
                     DO b = 1, virtual
                        W_ijab(i, j, a, b) = W_ijab(i, j, a, b) + B_bar_ijQ_bse_local(i, j, LLL)*B_abQ_bse_local(a, b, LLL)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         ! we make the mp_sum to sum over all RI basis functions
         CALL para_env%sum(W_ijab)
         IF (unit_nr > 0) THEN
            PRINT *, "MG: printing W_ijab"
            PRINT *, W_ijab
         END IF
      END IF

      ! we make the mp_sum to sum over all RI basis functions
      CALL para_env%sum(AZ)

      !print *, 'ProcNr',para_env%mepos,"Summed over group"

      !MG to del
      !print *, "For comparison, first the array of eigenvalues:"
      !print *, Eigenval
      !print *, "Diagonal of A"
      ! add (e_a-e_i)*Z_ia
      DO i_occ = 1, homo
         DO a_virt = 1, virtual

            eigen_diff = Eigenval(a_virt + homo) - Eigenval(i_occ)
            IF (unit_nr > 0 .AND. i_iter == 1) THEN
               PRINT *, "Ediff at (i_occ,a_virt)=", i_occ, a_virt, " is: ", eigen_diff
            END IF

            AZ(i_occ, a_virt, :) = AZ(i_occ, a_virt, :) + Z_vectors(i_occ, a_virt, :)*eigen_diff

         END DO
      END DO

      !MG to refine: cut off contributions, which are too high in the excitation spectrum
      IF (z_space_energy_cutoff > 0) THEN
         DO i_occ = 1, homo
            DO a_virt = 1, virtual

               IF (Eigenval(a_virt + homo) > z_space_energy_cutoff .OR. -Eigenval(i_occ) > z_space_energy_cutoff) THEN
                  AZ(i_occ, a_virt, :) = 0
               END IF

            END DO
         END DO
      END IF

      !Debugging purposes: full diagonalization of A
      IF (i_iter == 1 .AND. bse_full_diag_debug) THEN
         n = 0
         DO i = 1, homo
            DO a = 1, virtual
               n = n + 1
               m = 0
               DO j = 1, homo
                  DO b = 1, virtual
                     m = m + 1
                     IF (a == b .AND. i == j) THEN
                        eigen_diff = Eigenval(a + homo) - Eigenval(i)
                     ELSE
                        eigen_diff = 0
                     END IF
                     A_full_reshaped(n, m) = eigen_diff + 2*v_iajb(i, a, j, b) - W_ijab(i, j, a, b)
                     A_full(i, a, j, b) = eigen_diff + 2*v_iajb(i, a, j, b) - W_ijab(i, j, a, b)
                     IF (a == b .AND. i == j .AND. unit_nr > 0) THEN
            PRINT *, "MG: i,a, A on diag, Ediff, v, W", i, a, A_full(i, a, j, b), eigen_diff, v_iajb(i, a, j, b), W_ijab(i, j, a, b)
                     END IF
                  END DO
               END DO
            END DO
         END DO

         IF (unit_nr > 0) THEN
            PRINT *, "MG: Full A:"
            PRINT *, A_full
            PRINT *, "MG: Full A reshaped:"
            PRINT *, A_full_reshaped
         END IF

         !MG to do: Check for symmetry of ZAZ!
         ALLOCATE (WORK(1))
         WORK = 0.0_dp
         CALL DSYEV("N", "U", homo*virtual, A_full_reshaped, homo*virtual, Full_exc_spectrum, WORK, -1, diag_info)
         LWORK = INT(WORK(1))
         DEALLOCATE (WORK)
         ALLOCATE (WORK(LWORK))
         WORK = 0.0_dp
         !MG to check: Usage of symmetric routine okay? (Correct LWORK?)
         CALL DSYEV("N", "U", homo*virtual, A_full_reshaped, homo*virtual, Full_exc_spectrum, WORK, LWORK, diag_info)

         DEALLOCATE (WORK)

         IF (unit_nr > 0) THEN
            PRINT *, "MG: Full spectrum of A:"
            PRINT *, Full_exc_spectrum*evolt
         END IF
         DEALLOCATE (W_ijab, v_iajb, A_full, A_full_reshaped)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param RI_vector ...
!> \param local_RI_size ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
!> \param bse_spin_config ...
!> \param v_iajb ...
!> \param bse_full_diag_debug ...
!> \param i_iter ...
!> \param unit_nr ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_v_ia_jb_part(AZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                   num_Z_vectors, homo, virtual, bse_spin_config, v_iajb, bse_full_diag_debug, i_iter, &
                                   unit_nr, para_env)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: AZ, Z_vectors, B_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: RI_vector
      INTEGER, INTENT(IN)                                :: local_RI_size, num_Z_vectors, homo, &
                                                            virtual, bse_spin_config
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: v_iajb
      LOGICAL                                            :: bse_full_diag_debug
      INTEGER, INTENT(IN)                                :: i_iter, unit_nr
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_v_ia_jb_part'

      INTEGER                                            :: a, a_virt, b, handle, i, i_occ, &
                                                            i_Z_vector, j, LLL
      REAL(KIND=dp)                                      :: alpha

!debugging:

      CALL timeset(routineN, handle)

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      RI_vector = 0.0_dp

      ! v_P = sum_jb B_jb^P Z_jb
      DO LLL = 1, local_RI_size
         DO i_Z_vector = 1, num_Z_vectors
            DO i_occ = 1, homo
               DO a_virt = 1, virtual

                  RI_vector(LLL, i_Z_vector) = RI_vector(LLL, i_Z_vector) + &
                                               Z_vectors(i_occ, a_virt, i_Z_vector)* &
                                               B_iaQ_bse_local(i_occ, a_virt, LLL)

               END DO
            END DO
         END DO
      END DO

      ! AZ = sum_P B_ia^P*v_P + ...
      DO LLL = 1, local_RI_size
         DO i_Z_vector = 1, num_Z_vectors
            DO i_occ = 1, homo
               DO a_virt = 1, virtual
                  !MG to check: Minus sign at v oder W? Factor for triplet/singlet
                  AZ(i_occ, a_virt, i_Z_vector) = AZ(i_occ, a_virt, i_Z_vector) + &
                                                  alpha*RI_vector(LLL, i_Z_vector)* &
                                                  B_iaQ_bse_local(i_occ, a_virt, LLL)

               END DO
            END DO
         END DO
      END DO
      IF (i_iter == 1 .AND. bse_full_diag_debug) THEN
         ALLOCATE (v_iajb(homo, virtual, homo, virtual))
         v_iajb = 0.0_dp
         !MG to del: Print 4c integrals for debugging purposes
         !print *,"MG: ProcNr and localRIsize",para_env%mepos,local_RI_size
         DO LLL = 1, local_RI_size
            DO i = 1, homo
               DO j = 1, homo
                  DO a = 1, virtual
                     DO b = 1, virtual
                        v_iajb(i, a, j, b) = v_iajb(i, a, j, b) + B_iaQ_bse_local(i, a, LLL)*B_iaQ_bse_local(j, b, LLL)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         ! we make the mp_sum to sum over all RI basis functions
         CALL para_env%sum(v_iajb)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...Eigenval
!> \param Z_vectors ...
!> \param Eigenval ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE initial_guess_Z_vectors(Z_vectors, Eigenval, num_Z_vectors, homo, virtual)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Z_vectors
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: num_Z_vectors, homo, virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initial_guess_Z_vectors'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, &
                                                            min_loc(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigen_diff_ia

      CALL timeset(routineN, handle)

      ALLOCATE (eigen_diff_ia(homo, virtual))

      DO i_occ = 1, homo
         DO a_virt = 1, virtual
            eigen_diff_ia(i_occ, a_virt) = Eigenval(a_virt + homo) - Eigenval(i_occ)
         END DO
      END DO

      DO i_Z_vector = 1, num_Z_vectors

         min_loc = MINLOC(eigen_diff_ia)

         Z_vectors(min_loc(1), min_loc(2), i_Z_vector) = 1.0_dp

         eigen_diff_ia(min_loc(1), min_loc(2)) = 1.0E20_dp

      END DO

      DEALLOCATE (eigen_diff_ia)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_ia ...
! **************************************************************************************************
   SUBROUTINE mult_B_with_W(fm_mat_S_ij_bse, fm_mat_S_ab_bse, fm_mat_S, &
                            fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                            fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, &
                            dimen_RI, homo, virtual, dimen_ia)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                                            fm_mat_S
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q_static_bse, &
                                                            fm_mat_Q_static_bse_gemm
      INTEGER, INTENT(IN)                                :: dimen_RI, homo, virtual, dimen_ia

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_B_with_W'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_type)                                   :: fm_mat_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_S_bar_ia_bse, fm_mat_S%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_S, fm_mat_S_bar_ia_bse)
      CALL cp_fm_set_all(fm_mat_S_bar_ia_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_S_bar_ij_bse, fm_mat_S_ij_bse%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_S_ij_bse, fm_mat_S_bar_ij_bse)
      CALL cp_fm_set_all(fm_mat_S_bar_ij_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_work, fm_mat_Q_static_bse_gemm%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_Q_static_bse_gemm, fm_mat_work)
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      ! get info of fm_mat_Q_static_bse and compute ((1+Q(0))^-1-1)
      CALL cp_fm_get_info(matrix=fm_mat_Q_static_bse_gemm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               fm_mat_Q_static_bse_gemm%local_data(iiB, jjB) = fm_mat_Q_static_bse_gemm%local_data(iiB, jjB) + 1.0_dp
            END IF
         END DO
      END DO

      ! calculate Trace(Log(Matrix)) as Log(DET(Matrix)) via cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_mat_Q_static_bse_gemm, n=dimen_RI, info_out=info_chol)

      CPASSERT(info_chol == 0)

      ! calculate [1+Q(i0)]^-1
      CALL cp_fm_cholesky_invert(fm_mat_Q_static_bse_gemm)

      ! symmetrize the result
      CALL cp_fm_upper_to_full(fm_mat_Q_static_bse_gemm, fm_mat_work)

      CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=homo**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_Q_static_bse, matrix_b=fm_mat_S_ij_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_S_bar_ij_bse)
      ! MG: Remove dimen_ia and use homo*virtual?
      ! fm_mat_S_bar_ia_bse has a different blacs_env as fm_mat_S_ij_bse since we take
      ! fm_mat_S from RPA. Therefore, we also need a different fm_mat_Q_static_bse_gemm
      CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_ia, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_Q_static_bse_gemm, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_S_bar_ia_bse)

      CALL cp_fm_release(fm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE

   ! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param B_iaQ_bse_local ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param gd_array ...
!> \param color_sub ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fill_local_3c_arrays(fm_mat_S_ab_bse, fm_mat_S, &
                                   fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                   B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                   B_iaQ_bse_local, dimen_RI, homo, virtual, &
                                   gd_array, color_sub, para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ab_bse, fm_mat_S, &
                                                            fm_mat_S_bar_ia_bse, &
                                                            fm_mat_S_bar_ij_bse
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                            B_bar_iaQ_bse_local, B_iaQ_bse_local
      INTEGER, INTENT(IN)                                :: dimen_RI, homo, virtual
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_local_3c_arrays'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL allocate_and_fill_local_array(B_iaQ_bse_local, fm_mat_S, gd_array, color_sub, homo, virtual, dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_bar_iaQ_bse_local, fm_mat_S_bar_ia_bse, gd_array, color_sub, homo, virtual, &
                                         dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_bar_ijQ_bse_local, fm_mat_S_bar_ij_bse, gd_array, color_sub, homo, homo, &
                                         dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_abQ_bse_local, fm_mat_S_ab_bse, gd_array, color_sub, virtual, virtual, &
                                         dimen_RI, para_env)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param B_local ...
!> \param fm_mat_S ...
!> \param gd_array ...
!> \param color_sub ...
!> \param small_size ...
!> \param big_size ...
!> \param dimen_RI ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_local_array(B_local, fm_mat_S, gd_array, &
                                            color_sub, small_size, big_size, dimen_RI, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: B_local
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub, small_size, big_size, dimen_RI
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_fill_local_array'

      INTEGER :: combi_index, end_RI, handle, handle1, i_comm, i_entry, iiB, imepos, jjB, &
         level_big_size, level_small_size, ncol_local, nrow_local, num_comm_cycles, RI_index, &
         size_RI, start_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, mepos_from_RI_index, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: matrix_el
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req_array

      CALL timeset(routineN, handle)

      ALLOCATE (mepos_from_RI_index(dimen_RI))
      mepos_from_RI_index = 0

      DO imepos = 0, para_env%num_pe - 1

         CALL get_group_dist(gd_array, pos=imepos, starts=start_RI, ends=end_RI)

         mepos_from_RI_index(start_RI:end_RI) = imepos

      END DO

      ! color_sub is automatically the number of the process since every subgroup has only one MPI rank
      CALL get_group_dist(gd_array, color_sub, start_RI, end_RI, size_RI)

      ALLOCATE (B_local(small_size, big_size, 1:size_RI))

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))

      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      num_comm_cycles = 10

      ! communicate not all due to huge memory overhead, since for every number in fm_mat_S, we store
      ! three additional ones (RI index, first MO index, second MO index!!)
      DO i_comm = 0, num_comm_cycles - 1

         num_entries_send = 0
         num_entries_rec = 0

         ! loop over RI index to get the number of sent entries
         DO jjB = 1, nrow_local

            RI_index = row_indices(jjB)

            IF (MODULO(RI_index, num_comm_cycles) /= i_comm) CYCLE

            imepos = mepos_from_RI_index(RI_index)

            num_entries_send(imepos) = num_entries_send(imepos) + ncol_local

         END DO

         CALL para_env%alltoall(num_entries_send, num_entries_rec, 1)

         ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
         ALLOCATE (buffer_send(0:para_env%num_pe - 1))

         ! allocate data message and corresponding indices
         DO imepos = 0, para_env%num_pe - 1

            ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
            buffer_rec(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
            buffer_send(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
            buffer_rec(imepos)%indx = 0

            ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
            buffer_send(imepos)%indx = 0

         END DO

         entry_counter(:) = 0

         ! loop over RI index for filling the send-buffer
         DO jjB = 1, nrow_local

            RI_index = row_indices(jjB)

            IF (MODULO(RI_index, num_comm_cycles) /= i_comm) CYCLE

            imepos = mepos_from_RI_index(RI_index)

            DO iiB = 1, ncol_local

               combi_index = col_indices(iiB)
               level_small_size = MAX(1, combi_index - 1)/MAX(big_size, 2) + 1
               level_big_size = combi_index - (level_small_size - 1)*big_size

               entry_counter(imepos) = entry_counter(imepos) + 1

               buffer_send(imepos)%msg(entry_counter(imepos)) = fm_mat_S%local_data(jjB, iiB)

               buffer_send(imepos)%indx(entry_counter(imepos), 1) = RI_index
               buffer_send(imepos)%indx(entry_counter(imepos), 2) = level_small_size
               buffer_send(imepos)%indx(entry_counter(imepos), 3) = level_big_size

            END DO

         END DO

         CALL timeset("BSE_comm_data", handle1)

         CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

         CALL timestop(handle1)

         ! fill B_local
         DO imepos = 0, para_env%num_pe - 1

            DO i_entry = 1, num_entries_rec(imepos)

               RI_index = buffer_rec(imepos)%indx(i_entry, 1) - start_RI + 1
               level_small_size = buffer_rec(imepos)%indx(i_entry, 2)
               level_big_size = buffer_rec(imepos)%indx(i_entry, 3)

               matrix_el = buffer_rec(imepos)%msg(i_entry)
               !MG to del
               !print *,"level_small_size",level_small_size
               !print *,"level_big_size",level_big_size
               !print *,"RI_index",RI_index

               B_local(level_small_size, level_big_size, RI_index) = matrix_el

            END DO

         END DO

         DO imepos = 0, para_env%num_pe - 1
            DEALLOCATE (buffer_send(imepos)%msg)
            DEALLOCATE (buffer_send(imepos)%indx)
            DEALLOCATE (buffer_rec(imepos)%msg)
            DEALLOCATE (buffer_rec(imepos)%indx)
         END DO

         DEALLOCATE (buffer_rec, buffer_send)

      END DO

      DEALLOCATE (num_entries_send, num_entries_rec)

      DEALLOCATE (mepos_from_RI_index)

      DEALLOCATE (entry_counter, req_array)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param &!fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_bse( &!fm_mat_S_ij_bse,fm_mat_S_ab_bse, &
      fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_Q_static_bse, &
                                                            fm_mat_Q_static_bse_gemm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_bse'

      INTEGER                                            :: handle

!fm_mat_S_ij_bse,fm_mat_S_ab_bse, &

      CALL timeset(routineN, handle)

      !CALL cp_fm_release(fm_mat_S_ij_bse(1))
      !CALL cp_fm_release(fm_mat_S_ab_bse(1))
      CALL cp_fm_release(fm_mat_Q_static_bse)
      CALL cp_fm_release(fm_mat_Q_static_bse_gemm)

      CALL timestop(handle)
   END SUBROUTINE deallocate_matrices_bse

END MODULE bse
