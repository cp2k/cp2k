!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Build up the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      Joost VandeVondele (02.2002)
!>            1) rewrote collocate_pgf for increased accuracy and speed
!>            2) collocate_core hack for PGI compiler
!>            3) added multiple grid feature
!>            4) new way to go over the grid
!>      Joost VandeVondele (05.2002)
!>            1) prelim. introduction of the real space grid type
!>      JGH [30.08.02] multigrid arrays independent from potential
!>      JGH [17.07.03] distributed real space code
!>      JGH [23.11.03] refactoring and new loop ordering
!>      JGH [04.12.03] OpneMP parallelization of main loops
!>      Joost VandeVondele (12.2003)
!>           1) modified to compute tau
!>      Joost removed incremental build feature
!>      Joost introduced map consistent
!>      Rewrote grid integration/collocation routines, [Joost VandeVondele,03.2007]
!>      JGH [26.06.15] modification to allow for k-points
!> \author Matthias Krack (03.04.2001)
! **************************************************************************************************
MODULE qs_integrate_potential_product
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type
   USE grid_api,                        ONLY: grid_integrate_task_list
   USE input_constants,                 ONLY: do_admm_exch_scaling_merlot
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                              rs_grid_release,&
                                              rs_grid_retain
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE task_list_methods,               ONLY: rs_copy_to_buffer,&
                                              rs_copy_to_matrices,&
                                              rs_gather_matrices,&
                                              rs_scatter_matrices
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_product'

! *** Public subroutines ***
! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

   PUBLIC :: integrate_v_rspace

CONTAINS

! **************************************************************************************************
!> \brief computes matrix elements corresponding to a given potential
!> \param v_rspace ...
!> \param hmat ...
!> \param hmat_kp ...
!> \param pmat ...
!> \param pmat_kp ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param force_adm whether force of in aux. dens. matrix is calculated
!> \param ispin ...
!> \param compute_tau ...
!> \param gapw ...
!> \param basis_type ...
!> \param pw_env_external ...
!> \param task_list_external ...
!> \par History
!>      IAB (29-Apr-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!>      Some refactoring, get priorities for options correct (JGH, 04.2014)
!>      Added options to allow for k-points
!>      For a smooth transition we allow for old and new (vector) matrices (hmat, pmat) (JGH, 06.2015)
!> \note
!>     integrates a given potential (or other object on a real
!>     space grid) = v_rspace using a multi grid technique (mgrid_*)
!>     over the basis set producing a number for every element of h
!>     (should have the same sparsity structure of S)
!>     additional screening is available using the magnitude of the
!>     elements in p (? I'm not sure this is a very good idea)
!>     this argument is optional
!>     derivatives of these matrix elements with respect to the ionic
!>     coordinates can be computed as well
! **************************************************************************************************
   SUBROUTINE integrate_v_rspace(v_rspace, hmat, hmat_kp, pmat, pmat_kp, &
                                 qs_env, calculate_forces, force_adm, ispin, &
                                 compute_tau, gapw, basis_type, pw_env_external, task_list_external)

      TYPE(pw_p_type)                                    :: v_rspace
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hmat
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: hmat_kp
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: pmat
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: pmat_kp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_adm
      INTEGER, INTENT(IN), OPTIONAL                      :: ispin
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau, gapw
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace'

      CHARACTER(len=default_string_length)               :: my_basis_type
      INTEGER                                            :: atom_a, group, handle, i, iatom, &
                                                            igrid_level, ikind, img, maxco, &
                                                            maxsgf_set, natoms, nimages, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: calculate_virial, distributed_grids, &
                                                            do_kp, my_compute_tau, my_force_adm, &
                                                            my_gapw, pab_required
      REAL(KIND=dp)                                      :: admm_scal_fac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: forces_array
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_matrix
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap, dhmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! we test here if the provided operator matrices are consistent
      CPASSERT(PRESENT(hmat) .OR. PRESENT(hmat_kp))
      do_kp = .FALSE.
      IF (PRESENT(hmat_kp)) do_kp = .TRUE.
      IF (PRESENT(pmat)) THEN
         CPASSERT(PRESENT(hmat))
      ELSE IF (PRESENT(pmat_kp)) THEN
         CPASSERT(PRESENT(hmat_kp))
      END IF

      NULLIFY (pw_env, admm_env)

      ! this routine works in two modes:
      ! normal mode : <a| V | b>
      ! tau mode    : < nabla a| V | nabla b>
      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau

      my_force_adm = .FALSE.
      IF (PRESENT(force_adm)) my_force_adm = force_adm

      ! this sets the basis set to be used. GAPW(==soft basis) overwrites basis_type
      ! default is "ORB"
      my_gapw = .FALSE.
      IF (PRESENT(gapw)) my_gapw = gapw
      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      ! get the task lists
      ! task lists have to be in sync with basis sets
      ! there is an option to provide the task list from outside (not through qs_env)
      ! outside option has highest priority
      CALL get_qs_env(qs_env=qs_env, &
                      task_list=task_list, &
                      task_list_soft=task_list_soft)
      IF (.NOT. my_basis_type == "ORB") THEN
         CPASSERT(PRESENT(task_list_external))
      END IF
      IF (my_gapw) task_list => task_list_soft
      IF (PRESENT(task_list_external)) task_list => task_list_external
      CPASSERT(ASSOCIATED(task_list))

      ! the information on the grids is provided through pw_env
      ! pw_env has to be the parent env for the potential grid (input)
      ! there is an option to provide an external grid
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      IF (PRESENT(pw_env_external)) pw_env => pw_env_external

      ! get all the general information on the system we are working on
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      natom=natoms, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      force=force, &
                      virial=virial)

      admm_scal_fac = 1.0_dp
      IF (my_force_adm) THEN
         CALL get_qs_env(qs_env=qs_env, admm_env=admm_env)
         ! Calculate bare scaling of force according to Merlot, 1. IF: ADMMP, 2. IF: ADMMS,
         IF ((.NOT. admm_env%charge_constrain) .AND. &
             (admm_env%scaling_model == do_admm_exch_scaling_merlot)) THEN
            admm_scal_fac = admm_env%gsi(ispin)**2
         ELSE IF (admm_env%charge_constrain .AND. &
                  (admm_env%scaling_model == do_admm_exch_scaling_merlot)) THEN
            admm_scal_fac = (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)
         END IF
      END IF

      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
      END DO

      ! get mpi group from rs_v
      group = rs_v(1)%rs_grid%desc%group

      ! assign from pw_env
      gridlevel_info => pw_env%gridlevel_info
      cube_info => pw_env%cube_info

      ! transform the potential on the rs_multigrids
      CALL potential_pw2rs(rs_v, v_rspace, pw_env)

      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         CPASSERT(do_kp)
      END IF
      nkind = SIZE(qs_kind_set)
      calculate_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer) .AND. calculate_forces
      pab_required = (PRESENT(pmat) .OR. PRESENT(pmat_kp)) .AND. calculate_forces

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type=my_basis_type)

      distributed_grids = .FALSE.
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         IF (rs_v(igrid_level)%rs_grid%desc%distributed) THEN
            distributed_grids = .TRUE.
         ENDIF
      ENDDO

      ALLOCATE (forces_array(3, natoms))

      IF (pab_required) THEN
         ! initialize the working pmat structures
         ALLOCATE (deltap(nimages))
         IF (do_kp) THEN
            DO img = 1, nimages
               deltap(img)%matrix => pmat_kp(img)%matrix
            END DO
         ELSE
            deltap(1)%matrix => pmat%matrix
         END IF

         ! Distribute matrix blocks.
         IF (distributed_grids) THEN
            CALL rs_scatter_matrices(deltap, task_list%pab_buffer, task_list, group)
         ELSE
            CALL rs_copy_to_buffer(deltap, task_list%pab_buffer, task_list)
         ENDIF
         DEALLOCATE (deltap)
      END IF

      ! Map all tasks from the grids
      CALL grid_integrate_task_list(task_list=task_list%grid_task_list, &
                                    compute_tau=my_compute_tau, &
                                    calculate_forces=calculate_forces, &
                                    calculate_virial=calculate_virial, &
                                    pab_blocks=task_list%pab_buffer, &
                                    rs_grids=rs_v, &
                                    hab_blocks=task_list%hab_buffer, &
                                    forces=forces_array, &
                                    virial=virial_matrix)

      IF (calculate_forces) THEN
         ALLOCATE (atom_of_kind(natoms), kind_of(natoms))
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
!$OMP PARALLEL DO DEFAULT(NONE)  PRIVATE(atom_a, ikind) &
!$OMP             SHARED(natoms, force, forces_array, atom_of_kind, kind_of, admm_scal_fac)
         DO iatom = 1, natoms
            atom_a = atom_of_kind(iatom)
            ikind = kind_of(iatom)
            force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) + admm_scal_fac*forces_array(:, iatom)
         END DO
!$OMP END PARALLEL DO
         DEALLOCATE (atom_of_kind, kind_of)
      END IF

      IF (calculate_virial) THEN
         virial%pv_virial = virial%pv_virial + admm_scal_fac*virial_matrix
      END IF

      ! Gather all matrix images into a single array.
      ALLOCATE (dhmat(nimages))
      IF (PRESENT(hmat_kp)) THEN
         CPASSERT(.NOT. PRESENT(hmat))
         DO img = 1, nimages
            dhmat(img)%matrix => hmat_kp(img)%matrix
         END DO
      ELSE
         CPASSERT(PRESENT(hmat) .AND. nimages == 1)
         dhmat(1)%matrix => hmat%matrix
      END IF

      ! Distribute matrix blocks.
      IF (distributed_grids) THEN
         CALL rs_gather_matrices(task_list%hab_buffer, dhmat, task_list, group)
      ELSE
         CALL rs_copy_to_matrices(task_list%hab_buffer, dhmat, task_list)
      ENDIF
      DEALLOCATE (dhmat)

      IF (ASSOCIATED(rs_v)) THEN
         DO i = 1, SIZE(rs_v)
            CALL rs_grid_release(rs_v(i)%rs_grid)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE integrate_v_rspace

END MODULE qs_integrate_potential_product
