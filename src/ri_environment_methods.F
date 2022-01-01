!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates integral matrices for RIGPW method
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) heavily debugged
!>      updated for RI JGH [08.2017]
!> \authors JGH
!>          Dorothea Golze
! **************************************************************************************************
MODULE ri_environment_methods
   USE arnoldi_api,                     ONLY: arnoldi_conjugate_gradient
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_dot, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale_by_vector, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: allocate_lri_coefs,&
                                              lri_density_create,&
                                              lri_density_release,&
                                              lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE message_passing,                 ONLY: mp_sum
   USE pw_types,                        ONLY: pw_p_type
   USE qs_collocate_density,            ONLY: calculate_lri_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_o3c_methods,                  ONLY: calculate_o3c_integrals,&
                                              contract12_o3c
   USE qs_o3c_types,                    ONLY: get_o3c_iterator_info,&
                                              init_o3c_container,&
                                              o3c_container_type,&
                                              o3c_iterate,&
                                              o3c_iterator_create,&
                                              o3c_iterator_release,&
                                              o3c_iterator_type,&
                                              release_o3c_container
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE util,                            ONLY: get_limit

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_environment_methods'

   PUBLIC :: build_ri_matrices, calculate_ri_densities, ri_metric_solver

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_ri_matrices(lri_env, qs_env, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CALL calculate_ri_integrals(lri_env, qs_env, calculate_forces)

   END SUBROUTINE build_ri_matrices

! **************************************************************************************************
!> \brief calculates integrals needed for the RI density fitting,
!>        integrals are calculated once, before the SCF starts
!> \param lri_env ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE calculate_ri_integrals(lri_env, qs_env, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_ri_integrals'
      REAL(KIND=dp), DIMENSION(2), PARAMETER             :: fx = (/0.0_dp, 1.0_dp/)

      INTEGER                                            :: handle, i, i1, i2, iatom, ispin, izero, &
                                                            j, j1, j2, jatom, m, n, nbas
      INTEGER, DIMENSION(:, :), POINTER                  :: bas_ptr
      REAL(KIND=dp)                                      :: fpre
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eval
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: avec, fblk, fout
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p
      TYPE(dbcsr_type)                                   :: emat
      TYPE(dbcsr_type), POINTER                          :: fmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(lri_env))
      CPASSERT(ASSOCIATED(qs_env))

      ! overlap matrices
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      NULLIFY (rho, matrix_p)
      ! orbital overlap matrix (needed for N constraints and forces)
      ! we replicate this in order to not directly interact qith qs_env
      IF (calculate_forces) THEN
         ! charge constraint forces are calculated here
         CALL get_qs_env(qs_env=qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao=matrix_p)
         ALLOCATE (fmat)
         CALL dbcsr_create(fmat, template=matrix_p(1)%matrix)
         DO ispin = 1, dft_control%nspins
            fpre = lri_env%ri_fit%ftrm1n(ispin)/lri_env%ri_fit%ntrm1n
            CALL dbcsr_add(fmat, matrix_p(ispin)%matrix, fx(ispin), -fpre)
         END DO
         CALL build_overlap_matrix(ks_env=ks_env, matrix_s=lri_env%ob_smat, nderivative=0, &
                                   basis_type_a="ORB", basis_type_b="ORB", sab_nl=lri_env%soo_list, &
                                   calculate_forces=.TRUE., matrix_p=fmat)
         CALL dbcsr_release(fmat)
         DEALLOCATE (fmat)
      ELSE
         CALL build_overlap_matrix(ks_env=ks_env, matrix_s=lri_env%ob_smat, nderivative=0, &
                                   basis_type_a="ORB", basis_type_b="ORB", sab_nl=lri_env%soo_list)
      END IF
      ! RI overlap matrix
      CALL build_overlap_matrix(ks_env=ks_env, matrix_s=lri_env%ri_smat, nderivative=0, &
                                basis_type_a="RI_HXC", basis_type_b="RI_HXC", sab_nl=lri_env%saa_list)
      IF (calculate_forces) THEN
         ! calculate f*a(T) pseudo density matrix for forces
         bas_ptr => lri_env%ri_fit%bas_ptr
         avec => lri_env%ri_fit%avec
         fout => lri_env%ri_fit%fout
         ALLOCATE (fmat)
         CALL dbcsr_create(fmat, template=lri_env%ri_smat(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(fmat, lri_env%saa_list)
         CALL dbcsr_set(fmat, 0.0_dp)
         CALL dbcsr_iterator_start(iter, fmat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, fblk)
            i1 = bas_ptr(1, iatom)
            i2 = bas_ptr(2, iatom)
            j1 = bas_ptr(1, jatom)
            j2 = bas_ptr(2, jatom)
            IF (iatom <= jatom) THEN
               DO ispin = 1, dft_control%nspins
                  DO j = j1, j2
                     m = j - j1 + 1
                     DO i = i1, i2
                        n = i - i1 + 1
                        fblk(n, m) = fblk(n, m) + fout(i, ispin)*avec(j, ispin) + avec(i, ispin)*fout(j, ispin)
                     END DO
                  END DO
               END DO
            ELSE
               DO ispin = 1, dft_control%nspins
                  DO i = i1, i2
                     n = i - i1 + 1
                     DO j = j1, j2
                        m = j - j1 + 1
                        fblk(m, n) = fblk(m, n) + fout(i, ispin)*avec(j, ispin) + avec(i, ispin)*fout(j, ispin)
                     END DO
                  END DO
               END DO
            END IF
            IF (iatom == jatom) THEN
               fblk(:, :) = 0.25_dp*fblk(:, :)
            ELSE
               fblk(:, :) = 0.5_dp*fblk(:, :)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
         !
         CALL build_overlap_matrix(ks_env=ks_env, matrix_s=lri_env%ri_smat, nderivative=0, &
                                   basis_type_a="RI_HXC", basis_type_b="RI_HXC", sab_nl=lri_env%saa_list, &
                                   calculate_forces=.TRUE., matrix_p=fmat)
         CALL dbcsr_release(fmat)
         DEALLOCATE (fmat)
      END IF
      ! approximation (preconditioner) or exact inverse of RI overlap
      CALL dbcsr_allocate_matrix_set(lri_env%ri_sinv, 1)
      ALLOCATE (lri_env%ri_sinv(1)%matrix)
      SELECT CASE (lri_env%ri_sinv_app)
      CASE ("NONE")
         !
      CASE ("INVS")
         CALL dbcsr_create(lri_env%ri_sinv(1)%matrix, template=lri_env%ri_smat(1)%matrix)
         CALL invert_Hotelling(lri_env%ri_sinv(1)%matrix, lri_env%ri_smat(1)%matrix, &
                               threshold=1.e-10_dp, use_inv_as_guess=.FALSE., &
                               norm_convergence=1.e-10_dp, filter_eps=1.e-12_dp, silent=.FALSE.)
      CASE ("INVF")
         CALL dbcsr_create(emat, matrix_type=dbcsr_type_no_symmetry, template=lri_env%ri_smat(1)%matrix)
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL dbcsr_get_info(lri_env%ri_smat(1)%matrix, nfullrows_total=nbas)
         ALLOCATE (eval(nbas))
         CALL cp_dbcsr_syevd(lri_env%ri_smat(1)%matrix, emat, eval, para_env, blacs_env)
         izero = 0
         DO i = 1, nbas
            IF (eval(i) < 1.0e-10_dp) THEN
               eval(i) = 0.0_dp
               izero = izero + 1
            ELSE
               eval(i) = SQRT(1.0_dp/eval(i))
            END IF
         END DO
         CALL dbcsr_scale_by_vector(emat, eval, side='right')
         CALL dbcsr_create(lri_env%ri_sinv(1)%matrix, template=lri_env%ri_smat(1)%matrix)
         CALL dbcsr_multiply("N", "T", 1.0_dp, emat, emat, 0.0_dp, lri_env%ri_sinv(1)%matrix)
         DEALLOCATE (eval)
         CALL dbcsr_release(emat)
      CASE ("AINV")
         CALL dbcsr_create(lri_env%ri_sinv(1)%matrix, template=lri_env%ri_smat(1)%matrix)
         CALL invert_Hotelling(lri_env%ri_sinv(1)%matrix, lri_env%ri_smat(1)%matrix, &
                               threshold=1.e-5_dp, use_inv_as_guess=.FALSE., &
                               norm_convergence=1.e-4_dp, filter_eps=1.e-4_dp, silent=.FALSE.)
      CASE DEFAULT
         CPABORT("Unknown RI_SINV type")
      END SELECT

      ! solve Rx=n
      CALL ri_metric_solver(mat=lri_env%ri_smat(1)%matrix, &
                            vecr=lri_env%ri_fit%nvec, &
                            vecx=lri_env%ri_fit%rm1n, &
                            matp=lri_env%ri_sinv(1)%matrix, &
                            solver=lri_env%ri_sinv_app, &
                            ptr=lri_env%ri_fit%bas_ptr)

      ! calculate n(t)x
      lri_env%ri_fit%ntrm1n = SUM(lri_env%ri_fit%nvec(:)*lri_env%ri_fit%rm1n(:))

      ! calculate 3c-overlap integrals
      IF (ASSOCIATED(lri_env%o3c)) THEN
         CALL release_o3c_container(lri_env%o3c)
      ELSE
         ALLOCATE (lri_env%o3c)
      END IF
      CALL init_o3c_container(lri_env%o3c, dft_control%nspins, &
                              lri_env%orb_basis, lri_env%orb_basis, lri_env%ri_basis, &
                              lri_env%soo_list, lri_env%soa_list)

      CALL calculate_o3c_integrals(lri_env%o3c, calculate_forces, matrix_p)

      CALL timestop(handle)

   END SUBROUTINE calculate_ri_integrals
! **************************************************************************************************
!> \brief solver for RI systems (R*x=n)
!> \param mat ...
!> \param vecr ...
!> \param vecx ...
!> \param matp ...
!> \param solver ...
!> \param ptr ...
! **************************************************************************************************
   SUBROUTINE ri_metric_solver(mat, vecr, vecx, matp, solver, ptr)

      TYPE(dbcsr_type)                                   :: mat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vecr
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: vecx
      TYPE(dbcsr_type)                                   :: matp
      CHARACTER(LEN=*), INTENT(IN)                       :: solver
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ptr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ri_metric_solver'

      INTEGER                                            :: handle, max_iter, n
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: rerror, threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vect

      CALL timeset(routineN, handle)

      threshold = 1.e-10_dp
      max_iter = 100

      vecx(:) = vecr(:)

      SELECT CASE (solver)
      CASE ("NONE")
         CALL arnoldi_conjugate_gradient(mat, vecx, &
                                         converged=converged, threshold=threshold, max_iter=max_iter)
      CASE ("INVS", "INVF")
         converged = .FALSE.
         CALL ri_matvec(matp, vecr, vecx, ptr)
      CASE ("AINV")
         CALL arnoldi_conjugate_gradient(mat, vecx, matrix_p=matp, &
                                         converged=converged, threshold=threshold, max_iter=max_iter)
      CASE DEFAULT
         CPABORT("Unknown RI solver")
      END SELECT

      IF (.NOT. converged) THEN
         ! get error
         rerror = 0.0_dp
         n = SIZE(vecr)
         ALLOCATE (vect(n))
         CALL ri_matvec(mat, vecx, vect, ptr)
         vect(:) = vect(:) - vecr(:)
         rerror = MAXVAL(ABS(vect(:)))
         DEALLOCATE (vect)
         IF (rerror > threshold) THEN
            CPWARN("RI solver: CG did not converge properly")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ri_metric_solver

! **************************************************************************************************
!> \brief performs the fitting of the density and distributes the fitted
!>        density on the grid
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!>        lri_rho_struct where the fitted density is stored
!> \param qs_env ...
!> \param pmatrix ...
!> \param lri_rho_struct ...
!> \param atomic_kind_set ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE calculate_ri_densities(lri_env, qs_env, pmatrix, &
                                     lri_rho_struct, atomic_kind_set, para_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix
      TYPE(qs_rho_type), POINTER                         :: lri_rho_struct
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_ri_densities'

      INTEGER                                            :: atom_a, handle, i1, i2, iatom, ikind, &
                                                            ispin, n, natom, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :), POINTER                  :: bas_ptr
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: avec
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r

      CALL timeset(routineN, handle)

      nspin = SIZE(pmatrix, 1)
      CALL contract12_o3c(lri_env%o3c, pmatrix)
      CALL calculate_tvec_ri(lri_env, lri_env%o3c, para_env)
      CALL calculate_avec_ri(lri_env, pmatrix)
      !
      CALL get_qs_env(qs_env, lri_density=lri_density, natom=natom)
      CALL lri_density_release(lri_density)
      CALL lri_density_create(lri_density)
      lri_density%nspin = nspin
      ! allocate the arrays to hold RI expansion coefficients lri_coefs
      CALL allocate_lri_coefs(lri_env, lri_density, atomic_kind_set)
      ! reassign avec
      avec => lri_env%ri_fit%avec
      bas_ptr => lri_env%ri_fit%bas_ptr
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      DO ispin = 1, nspin
         lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atom_a = atom_of_kind(iatom)
            i1 = bas_ptr(1, iatom)
            i2 = bas_ptr(2, iatom)
            n = i2 - i1 + 1
            lri_coef(ikind)%acoef(atom_a, 1:n) = avec(i1:i2, ispin)
         END DO
      END DO
      CALL set_qs_env(qs_env, lri_density=lri_density)
      DEALLOCATE (atom_of_kind, kind_of)
      !
      CALL qs_rho_get(lri_rho_struct, rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r)
      DO ispin = 1, nspin
         lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
         CALL calculate_lri_rho_elec(rho_g(ispin), rho_r(ispin), qs_env, &
                                     lri_coef, tot_rho_r(ispin), "RI_HXC", .FALSE.)
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_ri_densities

! **************************************************************************************************
!> \brief assembles the global vector t=<P.T>
!> \param lri_env the lri environment
!> \param o3c ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE calculate_tvec_ri(lri_env, o3c, para_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calculate_tvec_ri'
      INTEGER, PARAMETER                                 :: msweep = 32

      INTEGER                                            :: handle, i1, i2, ibl, ibu, il, ispin, &
                                                            isweep, it, iu, katom, m, ma, mba, &
                                                            mepos, natom, nspin, nsweep, nthread
      INTEGER, DIMENSION(2, msweep)                      :: nlimit
      INTEGER, DIMENSION(:, :), POINTER                  :: bas_ptr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rm1t, tvec, tvl
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      nspin = SIZE(lri_env%ri_fit%tvec, 2)
      bas_ptr => lri_env%ri_fit%bas_ptr
      tvec => lri_env%ri_fit%tvec
      tvec = 0.0_dp
      natom = SIZE(bas_ptr, 2)
      nthread = 1
!$    nthread = omp_get_max_threads()

      IF (natom < 1000) THEN
         nsweep = 1
      ELSE
         nsweep = MIN(nthread, msweep)
      END IF

      nlimit = 0
      DO isweep = 1, nsweep
         nlimit(1:2, isweep) = get_limit(natom, nsweep, isweep - 1)
      END DO

      DO ispin = 1, nspin
         DO isweep = 1, nsweep
            il = nlimit(1, isweep)
            iu = nlimit(2, isweep)
            ma = iu - il + 1
            IF (ma < 1) CYCLE
            ibl = bas_ptr(1, il)
            ibu = bas_ptr(2, iu)
            mba = ibu - ibl + 1
            ALLOCATE (ta(mba, nthread))
            ta = 0.0_dp

            CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,o3c_iterator,ispin,il,iu,ibl,ta,bas_ptr)&
!$OMP PRIVATE (mepos,katom,tvl,i1,i2,m)

            mepos = 0
!$          mepos = omp_get_thread_num()

            DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
               CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, katom=katom, tvec=tvl)
               IF (katom < il .OR. katom > iu) CYCLE
               i1 = bas_ptr(1, katom) - ibl + 1
               i2 = bas_ptr(2, katom) - ibl + 1
               m = i2 - i1 + 1
               ta(i1:i2, mepos + 1) = ta(i1:i2, mepos + 1) + tvl(1:m, ispin)
            END DO
!$OMP END PARALLEL
            CALL o3c_iterator_release(o3c_iterator)

            ! sum over threads
            DO it = 1, nthread
               tvec(ibl:ibu, ispin) = tvec(ibl:ibu, ispin) + ta(1:mba, it)
            END DO
            DEALLOCATE (ta)
         END DO
      END DO

      ! global summation
      CALL mp_sum(tvec, para_env%group)

      rm1t => lri_env%ri_fit%rm1t
      DO ispin = 1, nspin
         ! solve Rx=t
         CALL ri_metric_solver(mat=lri_env%ri_smat(1)%matrix, &
                               vecr=tvec(:, ispin), &
                               vecx=rm1t(:, ispin), &
                               matp=lri_env%ri_sinv(1)%matrix, &
                               solver=lri_env%ri_sinv_app, &
                               ptr=bas_ptr)
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_tvec_ri
! **************************************************************************************************
!> \brief performs the fitting of the density in the RI method
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!> \param pmatrix ...
! **************************************************************************************************
   SUBROUTINE calculate_avec_ri(lri_env, pmatrix)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calculate_avec_ri'

      INTEGER                                            :: handle, ispin, nspin
      REAL(KIND=dp)                                      :: etr, nelec, nrm1t
      REAL(KIND=dp), DIMENSION(:), POINTER               :: nvec, rm1n
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: avec, rm1t, tvec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: smatrix

      CALL timeset(routineN, handle)

      nspin = SIZE(pmatrix)
      ! number of electrons
      smatrix => lri_env%ob_smat
      DO ispin = 1, nspin
         etr = 0.0_dp
         CALL dbcsr_dot(smatrix(1)%matrix, pmatrix(ispin)%matrix, etr)
         lri_env%ri_fit%echarge(ispin) = etr
      END DO

      tvec => lri_env%ri_fit%tvec
      rm1t => lri_env%ri_fit%rm1t
      nvec => lri_env%ri_fit%nvec
      rm1n => lri_env%ri_fit%rm1n

      ! calculate lambda
      DO ispin = 1, nspin
         nelec = lri_env%ri_fit%echarge(ispin)
         nrm1t = SUM(nvec(:)*rm1t(:, ispin))
         lri_env%ri_fit%lambda(ispin) = 2.0_dp*(nrm1t - nelec)/lri_env%ri_fit%ntrm1n
      END DO

      ! calculate avec = rm1t - lambda/2 * rm1n
      avec => lri_env%ri_fit%avec
      DO ispin = 1, nspin
         avec(:, ispin) = rm1t(:, ispin) - 0.5_dp*lri_env%ri_fit%lambda(ispin)*rm1n(:)
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_avec_ri

! **************************************************************************************************
!> \brief Mutiplies a replicated vector with a DBCSR matrix: vo = mat*vi
!> \param mat ...
!> \param vi ...
!> \param vo ...
!> \param ptr ...
! **************************************************************************************************
   SUBROUTINE ri_matvec(mat, vi, vo, ptr)

      TYPE(dbcsr_type)                                   :: mat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vi
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: vo
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ptr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ri_matvec'

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: group, handle, iatom, jatom, m1, m2, mb, &
                                                            n1, n2, nb
      LOGICAL                                            :: symm
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat, matrix_type=matrix_type, group=group)

      SELECT CASE (matrix_type)
      CASE (dbcsr_type_no_symmetry)
         symm = .FALSE.
      CASE (dbcsr_type_symmetric)
         symm = .TRUE.
      CASE (dbcsr_type_antisymmetric)
         CPABORT("NYI, antisymmetric matrix not permitted")
      CASE DEFAULT
         CPABORT("Unknown matrix type, ...")
      END SELECT

      vo(:) = 0.0_dp
      CALL dbcsr_iterator_start(iter, mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, block)
         n1 = ptr(1, iatom)
         n2 = ptr(2, iatom)
         nb = n2 - n1 + 1
         m1 = ptr(1, jatom)
         m2 = ptr(2, jatom)
         mb = m2 - m1 + 1
         CPASSERT(nb == SIZE(block, 1))
         CPASSERT(mb == SIZE(block, 2))
         vo(n1:n2) = vo(n1:n2) + MATMUL(block, vi(m1:m2))
         IF (symm .AND. (iatom /= jatom)) THEN
            vo(m1:m2) = vo(m1:m2) + MATMUL(TRANSPOSE(block), vi(n1:n2))
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL mp_sum(vo, group)

      CALL timestop(handle)

   END SUBROUTINE ri_matvec

END MODULE ri_environment_methods
