!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines saving to files various data related to POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_read

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              mo_set_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_read'

   PUBLIC :: pod_read_restart_data

CONTAINS

! **************************************************************************************************
!> \brief read transformation matrix from the restart file
!> \param mat storage for the transformation matrix
!> \param fm_s full-matrix structure used for allocation of matrix
!> \param desc description of the transformation matrix
!> \param filename name of the transformation-matrix restart file
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_read_restart_tmat(mat, fm_s, desc, filename)

      ! Routine arguments
      TYPE(cp_fm_type), POINTER                          :: mat
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      CHARACTER(LEN=*)                                   :: desc
      CHARACTER(LEN=default_path_length)                 :: filename

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_read_restart_tmat', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr

! Routine name for debug purposes
! Local variables

      ! Memory allocation
      CALL cp_fm_create(matrix=mat, matrix_struct=fm_s, &
                        name=desc//' TRANSFORMATION MATRIX')

      ! Open the file
      CALL open_file(file_name=filename, &
                     file_action="READ", &
                     file_form="UNFORMATTED", &
                     file_status="OLD", &
                     unit_number=unit_nr)

      ! Read the matrix elements
      CALL cp_fm_read_unformatted(mat, unit_nr)

      ! Close the file
      IF (unit_nr > 0) &
         CALL close_file(unit_nr)

   END SUBROUTINE pod_read_restart_tmat

! **************************************************************************************************
!> \brief read coupling elements and localized wavefunctions from restart files
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param ns number of spin components
!> \param file_hab name of the Hamiltonian restart file
!> \param file_wfn name of the wavefunction restart file
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_read_restart_cpl(qs_env, ec, ns, file_hab, file_wfn)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: ns
      CHARACTER(LEN=default_path_length)                 :: file_hab, file_wfn

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_read_restart_cpl', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, id, j, k, n_ao, n_blocks, n_el, n_mo, &
                                                            n_spins, unit_hab, unit_wfn
      INTEGER, DIMENSION(:), POINTER                     :: vi
      REAL(KIND=dp)                                      :: mx_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_cntrl

! Routine name for debug purposes
! Local variables

      ! Initialization
      NULLIFY (blacs_env)
      NULLIFY (fm_s)
      NULLIFY (mo)
      NULLIFY (para_env)
      NULLIFY (qs_kind_set)
      NULLIFY (scf_cntrl)
      NULLIFY (vi)

      ! Number of basis functions
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Open the files
      CALL open_file(file_name=file_hab, &
                     file_action="READ", &
                     file_form="UNFORMATTED", &
                     file_status="OLD", &
                     unit_number=unit_hab)
      CALL open_file(file_name=file_wfn, &
                     file_action="READ", &
                     file_form="UNFORMATTED", &
                     file_status="OLD", &
                     unit_number=unit_wfn)

      ! Number of spin components
      READ (UNIT=unit_hab) n_spins
      IF (ns /= n_spins) &
         CPABORT('inconsistent number of spin components in Hamiltonian restart file')
      READ (UNIT=unit_wfn) n_spins
      IF (ns /= n_spins) &
         CPABORT('inconsistent number of spin components in wavefunction restart file')

      ! Maximal occupancy
      mx_occ = 2.0_dp
      IF (n_spins > 1) &
         mx_occ = 1.0_dp

      ! Number of blocks
      READ (UNIT=unit_hab) n_blocks
      IF (ec%n_blocks /= n_blocks) &
         CPABORT('inconsistent number of POD blocks in Hamiltonian restart file')
      READ (UNIT=unit_wfn) n_blocks
      IF (ec%n_blocks /= n_blocks) &
         CPABORT('inconsistent number of POD blocks in wavefunction restart file')

      ! Block dimensions
      ALLOCATE (vi(n_blocks))
      READ (UNIT=unit_hab) (vi(i), i=1, n_blocks)
      DO i = 1, n_blocks
         IF (vi(i) /= ec%block(i)%n_ao) &
            CPABORT('inconsistent POD block dimensions in Hamiltonian restart file')
      END DO
      READ (UNIT=unit_wfn) (vi(i), i=1, n_blocks)
      DO i = 1, n_blocks
         IF (vi(i) /= ec%block(i)%n_ao) &
            CPABORT('inconsistent POD block dimensions in wavefunction restart file')
      END DO
      DEALLOCATE (vi)

      ! Memory allocation
      DO i = 1, n_blocks
         ! Memory for states
         CPASSERT(.NOT. ASSOCIATED(ec%block(i)%mo))
         ALLOCATE (ec%block(i)%mo(n_spins))
         CPASSERT(ASSOCIATED(ec%block(i)%mo))
         DO j = 1, n_spins
            NULLIFY (ec%block(i)%mo(j)%mo_set)
         END DO
         ! Memory for couplings
         ec%block(i)%n_hab_blocks = ec%n_blocks - i
         CPASSERT(.NOT. ASSOCIATED(ec%block(i)%hab))
         IF (ec%block(i)%n_hab_blocks > 0) THEN
            ALLOCATE (ec%block(i)%hab(n_spins, ec%block(i)%n_hab_blocks))
            CPASSERT(ASSOCIATED(ec%block(i)%hab))
            DO j = 1, n_spins
               DO k = 1, ec%block(i)%n_hab_blocks
                  NULLIFY (ec%block(i)%hab(j, k)%matrix)
               END DO
            END DO
         END IF
      END DO

      ! Spin components
      DO i = 1, n_spins
         ! Diagonal blocks
         DO j = 1, n_blocks

            ! Number of electrons
            n_el = ec%block(j)%n_electrons
            IF (n_spins > 1) THEN
               n_el = n_el/2
               IF (MOD(ec%block(j)%n_electrons, 2) == 1) THEN
                  IF (i == 1) &
                     n_el = n_el + 1
               END IF
            END IF

            ! State energies
            n_mo = ec%block(j)%n_ao

            NULLIFY (ec%block(j)%mo(i)%mo_set)
            CALL allocate_mo_set(ec%block(j)%mo(i)%mo_set, n_ao, n_mo, n_el, &
                                 REAL(n_el, dp), mx_occ, 0.0_dp)
            mo => ec%block(j)%mo(i)%mo_set

            NULLIFY (mo%eigenvalues)
            ALLOCATE (mo%eigenvalues(n_mo))
            CPASSERT(ASSOCIATED(mo%eigenvalues))

            READ (UNIT=unit_hab) (mo%eigenvalues(k), k=1, n_mo)

            ! MO coefficients
            CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
            CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                                     nrow_global=n_ao, ncol_global=n_mo)
            CALL cp_fm_create(matrix=mo%mo_coeff, matrix_struct=fm_s, name='BLOCK STATES')

            CALL cp_fm_read_unformatted(mo%mo_coeff, unit_wfn)

            ! Occupancies
            NULLIFY (mo%occupation_numbers)
            ALLOCATE (mo%occupation_numbers(n_mo))
            CPASSERT(ASSOCIATED(mo%occupation_numbers))
            mo%occupation_numbers = 0.0_dp

            IF (n_el > 0) THEN
               CALL get_qs_env(qs_env, scf_control=scf_cntrl)
               CALL set_mo_occupation(mo_set=mo, smear=scf_cntrl%smear)
            END IF

            ! Clean memory
            CALL cp_fm_struct_release(fmstruct=fm_s)
            NULLIFY (fm_s)

         END DO

         ! Off-diagonal blocks
         DO j = 1, n_blocks
            id = 1
            DO k = 1, n_blocks

               ! Coupling elements
               IF (k > j) THEN

                  ! Memory allocation
                  CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                                           nrow_global=ec%block(j)%n_ao, ncol_global=ec%block(k)%n_ao)
                  CALL cp_fm_create(matrix=ec%block(j)%hab(i, id)%matrix, matrix_struct=fm_s, &
                                    name='H_KS OFF-DIAGONAL BLOCK')

                  CALL cp_fm_read_unformatted(ec%block(j)%hab(i, id)%matrix, unit_hab)

                  id = id + 1

                  ! Clean memory
                  CALL cp_fm_struct_release(fmstruct=fm_s)
                  NULLIFY (fm_s)

               END IF

            END DO

         END DO
      END DO

      ! Close file
      IF (unit_hab > 0) &
         CALL close_file(unit_hab)
      IF (unit_wfn > 0) &
         CALL close_file(unit_wfn)

   END SUBROUTINE pod_read_restart_cpl

! **************************************************************************************************
!> \brief read POD Hamiltonian and MOs from binary restart file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin components
!> \param fm_s full-matrix structure used for allocation of transformation matrices
!> \param file_tmf name of the forward-transformation-matrix restart file
!> \param file_tmr name of the reverse-transformation-matrix restart file
!> \param file_hab name of the Hamiltonian restart file
!> \param file_wfn name of the wavefunction restart file
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_read_restart_data(qs_env, ec, n_spins, fm_s, file_tmf, file_tmr, file_hab, file_wfn)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      CHARACTER(LEN=default_path_length)                 :: file_tmf, file_tmr, file_hab, file_wfn

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_read_restart_data', &
         routineP = moduleN//':'//routineN

! Routine name for debug purposes

      ! Forward transformation matrix
      CALL pod_read_restart_tmat(ec%tmatf, fm_s, 'S^(-1/2)', file_tmf)
      ! Reverse transformation matrix
      CALL pod_read_restart_tmat(ec%tmatr, fm_s, 'S^(+1/2)', file_tmr)

      ! Coupling elements and wavefunctions
      CALL pod_read_restart_cpl(qs_env, ec, n_spins, file_hab, file_wfn)

   END SUBROUTINE pod_read_restart_data

END MODULE et_coupling_pod_read
