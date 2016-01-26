!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods for deltaSCF calculations
! *****************************************************************************
MODULE qs_mom_methods

  USE bibliography,                    ONLY: Gilbert2008,&
                                             cite_reference
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type,&
                                             cp_fm_vectorsnorm,&
                                             cp_fm_vectorssum
  USE input_constants,                 ONLY: momproj_norm,&
                                             momproj_sum
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             duplicate_mo_set,&
                                             get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             set_mo_set
  USE qs_scf_diagonalization,          ONLY: general_eigenproblem
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE string_utilities,                ONLY: integer_to_string
  USE util,                            ONLY: sort,&
                                             sort_unique
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mom_methods'

  PUBLIC  :: do_mom_guess, do_mom_diag
  PRIVATE :: mom_is_unique_orbital_indices, mom_reoccupy_orbitals

CONTAINS

! *****************************************************************************
!> \brief check that every molecular orbital index appears only once in each
!>        (de-)occupation list supplied by user. Check that all the indices
!>        are positive integers and abort if it is not the case.
!> \param  iarr      list of molecular orbital indices to be checked
!> \retval is_unique .true. if all the elements are unique or the list contains
!>         exactly one 0 element (meaning no excitation)
!> \par History
!>      01.2016 created [Sergey Chulkov]
! *****************************************************************************
  FUNCTION mom_is_unique_orbital_indices(iarr) RESULT(is_unique)
    INTEGER, DIMENSION(:), POINTER           :: iarr
    LOGICAL                                  :: is_unique

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'mom_is_unique_orbital_indices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, norbs
    INTEGER, DIMENSION(:), POINTER           :: tmp_iarr

    CALL timeset(routineN, handle)

    CPASSERT(ASSOCIATED(iarr))
    norbs = SIZE(iarr)

    IF (norbs > 0) THEN
       ALLOCATE(tmp_iarr(norbs))

       tmp_iarr(:) = iarr(:)
       CALL sort_unique(tmp_iarr, is_unique)

       ! Ensure that all orbital indices are positive integers.
       ! A special value '0' means 'disabled keyword',
       ! it must appear once to be interpreted in such a way
       IF (tmp_iarr(1) < 0 .OR. (tmp_iarr(1) == 0 .AND. norbs > 1)) &
          CPABORT("MOM: all molecular orbital indices must be positive integer numbers")

       DEALLOCATE(tmp_iarr)
    ENDIF

    is_unique = .TRUE.

    CALL timestop(handle)

  END FUNCTION mom_is_unique_orbital_indices

! *****************************************************************************
!> \brief swap occupation numbers between molecular orbitals
!>        from occupation and de-occupation lists
!> \param mo_set        set of molecular orbitals
!> \param deocc_orb_set list of de-occupied orbital indices
!> \param occ_orb_set   list of newly occupied orbital indices
!> \param spin          spin component of the molecular orbitals;
!>                      to be used for diagnostic messages
!> \par History
!>      01.2016 created [Sergey Chulkov]
! *****************************************************************************
  SUBROUTINE mom_reoccupy_orbitals(mo_set, deocc_orb_set, occ_orb_set, spin)
    TYPE(mo_set_type), POINTER               :: mo_set
    INTEGER, DIMENSION(:), POINTER           :: deocc_orb_set, occ_orb_set
    CHARACTER(len=*), INTENT(in)             :: spin

    CHARACTER(len=*), PARAMETER :: routineN = 'mom_reoccupy_orbitals', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=10)                        :: str_iorb, str_norbs
    CHARACTER(len=3)                         :: str_prefix
    INTEGER                                  :: handle, homo, iorb, lfomo, &
                                                nao, nmo, norbs
    REAL(KIND=dp)                            :: maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occ_nums

    CALL timeset(routineN, handle)

    ! MOM electron excitation should preserve both the number of electrons and
    ! multiplicity of the electronic system thus ensuring the following constraint :
    ! norbs = SIZE(deocc_orb_set) == SIZE(occ_orb_set)
    norbs = SIZE(deocc_orb_set)

    ! the following assertion should never raise an exception
    CPASSERT(SIZE(deocc_orb_set) == SIZE(occ_orb_set))

    IF (deocc_orb_set(1) /= 0 .AND. occ_orb_set(1) /= 0) THEN

       ! MOM does not follow aufbau principle producing non-uniform orbitals occupation
       CALL set_mo_set(mo_set=mo_set, uniform_occupation=.FALSE.)
       CALL get_mo_set(mo_set=mo_set, maxocc=maxocc, &
          nao=nao, nmo=nmo, occupation_numbers=occ_nums)

       IF (deocc_orb_set(norbs) > nao .OR. occ_orb_set(norbs) > nao) THEN
          ! STOP: one of the molecular orbital index exceeds the number of atomic basis functions available
          CALL integer_to_string(nao, str_norbs)

          IF (deocc_orb_set(norbs) >= occ_orb_set(norbs)) THEN
             iorb = deocc_orb_set(norbs)
             str_prefix = 'de-'
          ELSE
             iorb = occ_orb_set(norbs)
             str_prefix = ''
          ENDIF
          CALL integer_to_string(iorb, str_iorb)

          CALL cp_abort(__LOCATION__, "Unable to "//TRIM(str_prefix)//"occupy "//&
               TRIM(spin)//" orbital No. "//TRIM(str_iorb)//&
               " since its index exceeds the number of atomic orbital functions available ("//&
               TRIM(str_norbs)//"). Please consider using a larger basis set.")
       ENDIF

       IF (deocc_orb_set(norbs) > nmo .OR. occ_orb_set(norbs) > nmo) THEN
          ! STOP: one of the molecular orbital index exceeds the number of constructed molecular orbitals
          IF (deocc_orb_set(norbs) >= occ_orb_set(norbs)) THEN
             iorb = deocc_orb_set(norbs)
          ELSE
             iorb = occ_orb_set(norbs)
          ENDIF

          IF (iorb-nmo > 1) THEN
             CALL integer_to_string(iorb-nmo, str_iorb)
             str_prefix = 's'
          ELSE
             str_iorb = 'an'
             str_prefix = ''
          ENDIF

          CALL integer_to_string(nmo, str_norbs)

          CALL cp_abort(__LOCATION__, "The number of molecular orbitals ("//TRIM(str_norbs)//&
               ") is not enough to perform MOM calculation. Please add "//&
               TRIM(str_iorb)//" extra orbital"//TRIM(str_prefix)//&
               " using the ADDED_MOS keyword in the SCF section of your input file.")
       ENDIF

       DO iorb = 1,norbs
          ! swap occupation numbers between two adjoint molecular orbitals
          IF (occ_nums(deocc_orb_set(iorb)) <= 0.0_dp) THEN
             CALL integer_to_string(deocc_orb_set(iorb), str_iorb)

             CALL cp_abort(__LOCATION__, "The "//TRIM(spin)//" orbital No. "//&
                  TRIM(str_iorb)//" is not occupied thus it cannot be deoccupied.")
          ENDIF

          IF (occ_nums(occ_orb_set(iorb)) > 0.0_dp) THEN
             CALL integer_to_string(occ_orb_set(iorb), str_iorb)

             CALL cp_abort(__LOCATION__, "The "//TRIM(spin)//" orbital No. "//&
                  TRIM(str_iorb)//" is already occupied thus it cannot be reoccupied.")
          ENDIF

          occ_nums(occ_orb_set(iorb)) = occ_nums(deocc_orb_set(iorb))
          occ_nums(deocc_orb_set(iorb)) = 0.0_dp
       ENDDO

       ! locate the lowest non-maxocc occupied orbital
       DO lfomo = 1,nmo
          IF (occ_nums(lfomo) /= maxocc) EXIT
       ENDDO

       ! locate the highest occupied orbital
       DO homo = nmo,1,-1
          IF (occ_nums(homo) > 0.0_dp) EXIT
       ENDDO

       CALL set_mo_set(mo_set=mo_set, homo=homo, lfomo=lfomo)

    ELSE IF (deocc_orb_set(1) /= 0 .OR. occ_orb_set(1) /= 0) THEN
       CALL cp_abort(__LOCATION__, &
            "Incorrect multiplicity of the MOM reference electronic state")
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE mom_reoccupy_orbitals

! *****************************************************************************
!> \brief initial guess for the maximum overlap method
!> \param nspins      number of spin components
!> \param mos         array of molecular orbitals
!> \param scf_control SCF control variables
!> \param p_rmpv      density matrix to be computed
!> \param scf_iter    default SCF iteration to activate the MOM
!> \par History
!>      01.2016 created [Sergey Chulkov]
! *****************************************************************************
  SUBROUTINE do_mom_guess(nspins, mos, scf_control, p_rmpv, scf_iter)
    INTEGER, INTENT(in)                      :: nspins
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_rmpv
    INTEGER, INTENT(in)                      :: scf_iter

    CHARACTER(len=*), PARAMETER :: routineN = 'do_mom_guess', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=10)                        :: str_iter
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    ! we are about to initialise the maximum overlap method,
    ! so cite the relevant reference first
    CALL cite_reference(Gilbert2008)

    ! ensure we do not have duplicated orbital indices
    IF (.NOT. &
        (mom_is_unique_orbital_indices(scf_control%diagonalization%mom_deoccA) .AND. &
         mom_is_unique_orbital_indices(scf_control%diagonalization%mom_deoccB) .AND. &
         mom_is_unique_orbital_indices(scf_control%diagonalization%mom_occA)   .AND. &
         mom_is_unique_orbital_indices(scf_control%diagonalization%mom_occB))) &
      CALL cp_abort(__LOCATION__, &
        "Duplicate orbital indices were found in the MOM section")

    ! ignore beta orbitals for spin-unpolarized calculations
    IF (nspins == 1 .AND. (scf_control%diagonalization%mom_deoccB(1) /= 0 &
                 .OR. scf_control%diagonalization%mom_occB(1) /= 0)) THEN

       CALL cp_warn(__LOCATION__,"Maximum overlap method will"//&
          " ignore beta orbitals since neither UKS nor ROKS calculation is performed")
    ENDIF

    ! compute the change in multiplicity and number of electrons
    IF (SIZE(scf_control%diagonalization%mom_deoccA) /= &
        SIZE(scf_control%diagonalization%mom_occA) .OR. &
        (nspins > 1 .AND. &
        SIZE(scf_control%diagonalization%mom_deoccB) /= &
        SIZE(scf_control%diagonalization%mom_occB))) THEN

       CALL cp_abort(__LOCATION__, "Incorrect multiplicity of the MOM reference"//&
            " electronic state or inconsistent number of electrons")
    ENDIF

    ! procees alpha orbitals
    IF (nspins >= 1) THEN
       CALL mom_reoccupy_orbitals(mos(1)%mo_set, &
            scf_control%diagonalization%mom_deoccA, &
            scf_control%diagonalization%mom_occA, 'alpha')

       CALL calculate_density_matrix(mos(1)%mo_set, p_rmpv(1)%matrix)
    ENDIF

    ! procees beta orbitals (if any)
    IF (nspins >= 2) THEN
       CALL mom_reoccupy_orbitals(mos(2)%mo_set, &
            scf_control%diagonalization%mom_deoccB, &
            scf_control%diagonalization%mom_occB, 'beta')

       CALL calculate_density_matrix(mos(2)%mo_set, p_rmpv(2)%matrix)
    ENDIF

    ! adjust the SCF iteration to activate the method if it has been set inappropriately
    IF (scf_control%diagonalization%mom_start < scf_iter) THEN
       IF (scf_control%diagonalization%mom_start > 0) THEN
          ! the inappropriately iteration number has been provided through the input file
          CALL integer_to_string(scf_iter, str_iter)
          CALL cp_warn(__LOCATION__, &
            "MOM procedure cannot be activated prior to SCF iteration no. "//&
            TRIM(str_iter)//" due to the initial SCF guess method used.")
       ENDIF
       scf_control%diagonalization%mom_start = scf_iter
    ENDIF

    ! MOM is now initialised properly
    scf_control%diagonalization%mom_didguess = .TRUE.

    CALL timestop(handle)

  END SUBROUTINE do_mom_guess

! *****************************************************************************
!> \brief do an SCF iteration, then compute occupation numbers of the new
!>  molecular orbitals according to their overlap with the previous ones
!> \param scf_env     SCF environment information
!> \param mos         array of molecular orbitals
!> \param matrix_ks   sparse Kohn-Sham matrix
!> \param matrix_s    sparse overlap matrix
!> \param scf_control SCF control variables
!> \param scf_section SCF input section
!> \param diis_step   have we done a DIIS step
!> \par History
!>      07.2014 created [Matt Watkins]
!>      01.2016 release version [Sergey Chulkov]
! *****************************************************************************
  SUBROUTINE do_mom_diag(scf_env,mos,matrix_ks,&
                             matrix_s,scf_control,scf_section,&
                             diis_step)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step

    CHARACTER(len=*), PARAMETER :: routineN = 'do_mom_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, iproj, ispin, &
                                                lfomo, nao, nmo, nspins
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: inds
    REAL(kind=dp)                            :: maxocc
    REAL(kind=dp), DIMENSION(:), POINTER     :: occ_nums, proj, tmp_occ_nums
    TYPE(cp_fm_struct_type), POINTER         :: mo_mo_fmstruct
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_orig, &
                                                overlap, svec
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos_orig

    CALL timeset(routineN, handle)

    IF (.NOT. scf_control%diagonalization%mom_didguess) &
      CALL cp_abort(__LOCATION__, &
       "The current implementation of the maximum overlap method is incompatible with the initial SCF guess")

    ! copy old molecular orbitals
    IF (scf_env%iter_count >= scf_control%diagonalization%mom_start) THEN
       ALLOCATE(mos_orig(SIZE(mos)))
       DO ispin=1,SIZE(mos)
          CALL duplicate_mo_set(mos_orig(ispin)%mo_set,mos(ispin)%mo_set)
       END DO
    END IF

    ! solve the eigenproblem
    CALL general_eigenproblem(scf_env, mos, matrix_ks, matrix_s, scf_control, scf_section, diis_step)

    ! number of spins == dft_control%nspins
    nspins = SIZE(matrix_ks)

    IF (scf_env%iter_count >= scf_control%diagonalization%mom_start) THEN
      DO ispin=1,nspins

         CALL get_mo_set(mo_set=mos(ispin)%mo_set, maxocc=maxocc, &
              nao=nao, nmo=nmo, occupation_numbers=occ_nums)

         ! TO DO: sparse-matrix variant; check if use_mo_coeff_b is set, and if yes use mo_coeff_b instead
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL get_mo_set(mo_set=mos_orig(ispin)%mo_set, mo_coeff=mo_coeff_orig)

         NULLIFY(svec, overlap, mo_mo_fmstruct)

         ! NOTE: whereas it is possible to store the result of the following
         ! matrix multiplication into work matrices scf_env%scf_work2(1:nspins),
         ! we create the our own work storage, since these mentioned works matrices
         ! are allocated only if neither OT nor davidson diagonalization is in use.
         CALL cp_fm_create(svec, mo_coeff%matrix_struct)

         ! S * C
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mo_coeff, svec, nmo)

         ! a temporary overlap matrix
         CALL cp_fm_struct_create(mo_mo_fmstruct, nrow_global=nmo, ncol_global=nmo, &
              para_env=mo_coeff%matrix_struct%para_env, context=mo_coeff%matrix_struct%context)
         CALL cp_fm_create(overlap, mo_mo_fmstruct)
         CALL cp_fm_struct_release(mo_mo_fmstruct)

         ! C(old occupied)^T * S * C(new)
         CALL cp_fm_column_scale(mo_coeff_orig, occ_nums)
         CALL cp_fm_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, mo_coeff_orig, svec, 0.0_dp, overlap)
         CALL cp_fm_release(svec)

         !CALL cp_fm_vectorsnorm(overlap,mos(ispin)%mo_set%occupation_numbers)

         ALLOCATE(proj(nmo))
         ALLOCATE(inds(nmo))
         ALLOCATE(tmp_occ_nums(nmo))

         ! project the new molecular orbitals into the space of the occupied old ones
         SELECT CASE(scf_control%diagonalization%mom_proj_formula)
         CASE (momproj_sum)
              ! proj_j = abs( \sum_i overlap(i, j) )
              CALL cp_fm_vectorssum(overlap, proj)

              DO iproj = 1,nmo
                 proj(iproj) = ABS(proj(iproj))
              END DO

         CASE (momproj_norm)
              ! proj_j = (\sum_i overlap(i, j)**2) ** 0.5
              CALL cp_fm_vectorsnorm(overlap, proj)

         CASE DEFAULT
              CPABORT("Unimplemented projection formula")
         END SELECT

         CALL cp_fm_release(overlap)


         tmp_occ_nums(:) = occ_nums(:)
         ! sort occupation numbers in ascending order
         CALL sort(tmp_occ_nums, nmo, inds)
         ! sort overlap projection in ascending order
         CALL sort(proj, nmo, inds)

         ! reorder occupation numbers according to overlap projections
         DO iproj = 1,nmo
            occ_nums(inds(iproj)) = tmp_occ_nums(iproj)
         ENDDO

         DEALLOCATE(tmp_occ_nums)
         DEALLOCATE(inds)
         DEALLOCATE(proj)

         ! locate the lowest non-fully occupied orbital
         DO lfomo = 1,nmo
            IF (occ_nums(lfomo) /= maxocc) EXIT
         ENDDO

         ! locate the highest occupied orbital
         DO homo = nmo,1,-1
            IF (occ_nums(homo) > 0.0_dp) EXIT
         ENDDO

         CALL set_mo_set(mo_set=mos(ispin)%mo_set, homo=homo, lfomo=lfomo)

      END DO

      DO ispin=1,SIZE(mos_orig)
         CALL deallocate_mo_set(mos_orig(ispin)%mo_set)
      ENDDO
      DEALLOCATE(mos_orig)
    ENDIF

    ! recompute density matrix
    DO ispin=1,nspins
       CALL calculate_density_matrix(mos(ispin)%mo_set, scf_env%p_mix_new(ispin, 1)%matrix)
    END DO

    CALL timestop(handle)

  END SUBROUTINE do_mom_diag

END MODULE qs_mom_methods
