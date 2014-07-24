!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the neighbor list data types and the corresponding functionality
! *****************************************************************************
MODULE fist_neighbor_list_types

  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE exclusion_types,                 ONLY: exclusion_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_neighbor_list_types'

! *****************************************************************************
  TYPE neighbor_kind_pairs_type
     INTEGER, POINTER, DIMENSION(:,:)  :: list, ij_kind
     INTEGER, POINTER, DIMENSION(:)    :: id_kind
     INTEGER, POINTER, DIMENSION(:)    :: grp_kind_start, grp_kind_end
     INTEGER                           :: cell_vector(3), npairs
     INTEGER                           :: ngrp_kind
     REAL(KIND=dp)                     :: rmax
     ! The *_scale arrays are scaling factors for the corresponding nonbonding
     ! interaction energies and forces for the pairs in 'list'. To keep the size
     ! of these arrays small, pairs whose interaction must be scaled are moved
     ! to beginning of the array 'list'. nscale is the number of elements in
     ! *_scale that are effectively used. This way one does not have to
     ! reallocate the *_scale arrays for every new scaled pair interaction.
     ! The field is_info is only used to switch between the regular nonbonded
     ! and the nonbonded14 splines for the van der waals interactions.
     REAL(KIND=dp), POINTER, DIMENSION(:)    :: ei_scale
     REAL(KIND=dp), POINTER, DIMENSION(:)    :: vdw_scale
     LOGICAL, POINTER, DIMENSION(:)          :: is_onfo
     INTEGER                                 :: nscale
  END TYPE neighbor_kind_pairs_type

! *****************************************************************************
  TYPE fist_neighbor_type
     TYPE(neighbor_kind_pairs_type), DIMENSION(:), POINTER :: neighbor_kind_pairs
     INTEGER                                               :: nlists
  END TYPE fist_neighbor_type

  PUBLIC :: neighbor_kind_pairs_type,&
            fist_neighbor_type,&
            fist_neighbor_init,&
            fist_neighbor_deallocate,&
            fist_neighbor_add

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param fist_neighbor ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fist_neighbor_deallocate(fist_neighbor,error)
    TYPE(fist_neighbor_type), POINTER        :: fist_neighbor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_neighbor_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(fist_neighbor)) THEN
       ! deallocate neighbor_kind_pairs
       IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs)) THEN
          DO i = 1, SIZE(fist_neighbor%neighbor_kind_pairs)
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%list)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%list, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%id_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%ij_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ij_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%ei_scale)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ei_scale, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%vdw_scale)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%vdw_scale, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%is_onfo)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%is_onfo, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
          DEALLOCATE(fist_neighbor%neighbor_kind_pairs, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(fist_neighbor, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE fist_neighbor_deallocate

! *****************************************************************************
!> \brief ...
!> \param fist_neighbor ...
!> \param ncell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fist_neighbor_init(fist_neighbor, ncell, error)
    TYPE(fist_neighbor_type), POINTER        :: fist_neighbor
    INTEGER, INTENT(IN)                      :: ncell(3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_neighbor_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, list_size, &
                                                nlistmin, stat
    LOGICAL                                  :: failure
    TYPE(neighbor_kind_pairs_type), &
      DIMENSION(:), POINTER                  :: new_pairs

    failure=.FALSE.
    CALL timeset ( routineN, handle )
    IF (.NOT.ASSOCIATED(fist_neighbor)) THEN
       ALLOCATE(fist_neighbor,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(fist_neighbor%neighbor_kind_pairs)
    ENDIF

    nlistmin=(2*MAXVAL(ncell)+1)**3
    IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs)) THEN
       IF (SIZE(fist_neighbor%neighbor_kind_pairs)<nlistmin) THEN
          ALLOCATE(new_pairs(nlistmin),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,SIZE(fist_neighbor%neighbor_kind_pairs)
             new_pairs(i)%list   =>fist_neighbor%neighbor_kind_pairs(i)%list
             list_size = SIZE(new_pairs(i)%list)
             ALLOCATE(new_pairs(i)%id_kind(list_size),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(new_pairs(i)%ei_scale(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(new_pairs(i)%vdw_scale(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(new_pairs(i)%is_onfo(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             NULLIFY(new_pairs(i)%ij_kind,&
                     new_pairs(i)%grp_kind_start,&
                     new_pairs(i)%grp_kind_end)
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%ij_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ij_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%id_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%ei_scale)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ei_scale, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%vdw_scale)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%vdw_scale, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%is_onfo)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%is_onfo, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          ENDDO
          DO i=SIZE(fist_neighbor%neighbor_kind_pairs)+1,nlistmin
             ALLOCATE(new_pairs(i)%list(2,0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(new_pairs(i)%id_kind(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             NULLIFY(new_pairs(i)%ij_kind,&
                     new_pairs(i)%grp_kind_start,&
                     new_pairs(i)%grp_kind_end)
             NULLIFY(new_pairs(i)%ei_scale,new_pairs(i)%vdw_scale,new_pairs(i)%is_onfo)
          ENDDO
          DEALLOCATE(fist_neighbor%neighbor_kind_pairs,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          fist_neighbor%neighbor_kind_pairs=>new_pairs
       ELSE
          DO i=1,SIZE(fist_neighbor%neighbor_kind_pairs)
             list_size = SIZE(fist_neighbor%neighbor_kind_pairs(i)%list)
             CALL reallocate(fist_neighbor%neighbor_kind_pairs(i)%id_kind,1,list_size)
          END DO
       ENDIF
    ELSE
       ALLOCATE(fist_neighbor%neighbor_kind_pairs(nlistmin),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i=1,nlistmin
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%list(2,0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ei_scale(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%vdw_scale(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%is_onfo(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          NULLIFY(fist_neighbor%neighbor_kind_pairs(i)%ij_kind,&
                  fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start,&
                  fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end)
       ENDDO
    ENDIF

    fist_neighbor%nlists=nlistmin
    DO i=1,nlistmin
       fist_neighbor%neighbor_kind_pairs(i)%npairs=0
       fist_neighbor%neighbor_kind_pairs(i)%list=HUGE(0)
       fist_neighbor%neighbor_kind_pairs(i)%id_kind=HUGE(0)
       fist_neighbor%neighbor_kind_pairs(i)%cell_vector=HUGE(0)
       fist_neighbor%neighbor_kind_pairs(i)%nscale = 0
    ENDDO
    CALL timestop ( handle )
  END SUBROUTINE fist_neighbor_init

! *****************************************************************************
!> \brief ...
!> \param neighbor_kind_pair ...
!> \param atom_a ...
!> \param atom_b ...
!> \param rab ...
!> \param check_spline ...
!> \param id_kind ...
!> \param skip ...
!> \param cell ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param exclusions ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fist_neighbor_add(neighbor_kind_pair, atom_a, atom_b, &
       rab, check_spline, id_kind, skip, cell, &
       ei_scale14, vdw_scale14, exclusions, error)
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(KIND=dp), DIMENSION(3)              :: rab
    LOGICAL, INTENT(OUT)                     :: check_spline
    INTEGER, INTENT(IN)                      :: id_kind
    LOGICAL, INTENT(IN)                      :: skip
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(IN)                :: ei_scale14, vdw_scale14
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL                               :: exclusions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_neighbor_add', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: eps_default = EPSILON(0.0_dp)*1.0E4_dp

    INTEGER                                  :: new_npairs, npairs, nscale, &
                                                old_npairs, stat
    INTEGER, DIMENSION(:), POINTER           :: new_id_kind
    INTEGER, DIMENSION(:, :), POINTER        :: new_list
    LOGICAL                                  :: ex_ei, ex_vdw, failure, &
                                                is_onfo
    REAL(KIND=dp), DIMENSION(3)              :: rabc

    IF (.NOT.PRESENT(exclusions)) THEN
       ex_ei = .FALSE.
       ex_vdw = .FALSE.
       is_onfo = .FALSE.
    ELSE
       ex_ei = ANY(exclusions(atom_a)%list_exclude_ei==atom_b)
       ex_vdw = ANY(exclusions(atom_a)%list_exclude_vdw==atom_b)
       is_onfo = ANY(exclusions(atom_a)%list_onfo==atom_b)
       IF (ex_ei.OR.ex_vdw.OR.is_onfo) THEN
          ! Check if this pair could correspond to a local interaction (bond, bend,
          ! or torsion) to which the exclusion lists and 14 potentials apply.
          !
          ! rab is the relative vector that may include some cell vectors. rabc is
          ! the 'shortest' possible relative vector, i.e. cell vectors are
          ! subtracted. When they are not the same, rab corresponds to a non-local
          ! interaction and the exclusion lists do not apply.
          rabc = pbc(rab, cell)
          IF ((ANY(ABS(rab-rabc)>eps_default))) THEN
             ex_ei = .FALSE.
             ex_vdw = .FALSE.
             is_onfo = .FALSE.
          END IF
       END IF
    ENDIF

    ! The skip option is .TRUE. for QM-QM pairs in an QM/MM run. In case these
    ! interactions have an ex_ei option, we store it in the neighbor list to
    ! do a proper bonded correction for the ewald summation. If there is no
    ! exclusion, the pair can be neglected.
    IF (skip.AND.(.NOT.ex_ei)) THEN
       ! If the pair is not present, checking is obviously not need.
       check_spline = .FALSE.
       RETURN
    END IF

    ! The check_spline is set to .TRUE. when the van derwaals is not excluded.
    ! Electrostatic interactions do not matter here as they are not evaluated
    ! with splines.
    check_spline = (.NOT.ex_vdw)

    ! If both types of interactions are excluded, the corresponding potentials
    ! will never be evaluated. At first sight such a pair would not need to be
    ! added to the neighborlists at all. However, they are still needed for
    ! proper corrections on interactions between the screening charges of bonded
    ! atoms when the ewald summation is used for the electrostatic interactions.

    ! If an interaction is excluded or scaled, store scale. If the interaction
    ! is an onfo, also store that property.
    IF (ex_ei.OR.ex_vdw.OR.is_onfo) THEN
       ! Allocate more memory for the scalings if necessary.
       nscale = neighbor_kind_pair%nscale
       IF (nscale==SIZE(neighbor_kind_pair%ei_scale)) THEN
          CALL reallocate(neighbor_kind_pair%ei_scale, 1, INT(5+1.2*nscale))
          CALL reallocate(neighbor_kind_pair%vdw_scale, 1, INT(5+1.2*nscale))
          CALL reallocate(neighbor_kind_pair%is_onfo, 1, INT(5+1.2*nscale))
       END IF
       nscale = nscale + 1
       IF (ex_ei) THEN
          neighbor_kind_pair%ei_scale(nscale) = 0.0_dp
       ELSE IF (is_onfo) THEN
          neighbor_kind_pair%ei_scale(nscale) = ei_scale14
       ELSE
          neighbor_kind_pair%ei_scale(nscale) = 1.0_dp
       END IF
       IF (ex_vdw) THEN
          neighbor_kind_pair%vdw_scale(nscale) = 0.0_dp
       ELSE IF (is_onfo) THEN
          neighbor_kind_pair%vdw_scale(nscale) = vdw_scale14
       ELSE
          neighbor_kind_pair%vdw_scale(nscale) = 1.0_dp
       END IF
       neighbor_kind_pair%is_onfo(nscale) = is_onfo
       neighbor_kind_pair%nscale = nscale
    ELSE
       nscale = HUGE(0)
    END IF

    ! Allocate more memory for the pair list if necessary.
    old_npairs=SIZE(neighbor_kind_pair%list,2)
    IF (old_npairs==neighbor_kind_pair%npairs) THEN
       failure = .FALSE.
       ! just a choice that will also grow for zero size arrays:
       new_npairs = INT(5+1.2*old_npairs)
       ! Pair Atoms Info
       ALLOCATE(new_list(2,new_npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       new_list(1:2,1:old_npairs)=neighbor_kind_pair%list(1:2,1:old_npairs)
       DEALLOCATE(neighbor_kind_pair%list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       neighbor_kind_pair%list=>new_list
       ! Kind Info
       ALLOCATE(new_id_kind(new_npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       new_id_kind(1:old_npairs)=neighbor_kind_pair%id_kind(1:old_npairs)
       DEALLOCATE(neighbor_kind_pair%id_kind,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       neighbor_kind_pair%id_kind=>new_id_kind
    ENDIF

    ! Store the pair ...
    npairs = neighbor_kind_pair%npairs + 1
    IF ((ex_ei.OR.ex_vdw.OR.is_onfo).AND.(npairs>nscale)) THEN
       ! ... after the previous pair that had scaling factors.
       neighbor_kind_pair%list(1,npairs) = neighbor_kind_pair%list(1,nscale)
       neighbor_kind_pair%list(2,npairs) = neighbor_kind_pair%list(2,nscale)
       neighbor_kind_pair%id_kind(npairs) = neighbor_kind_pair%id_kind(nscale)
       neighbor_kind_pair%list(1,nscale) = atom_a
       neighbor_kind_pair%list(2,nscale) = atom_b
       neighbor_kind_pair%id_kind(nscale) = id_kind
    ELSE
       ! ... at the end of the list.
       neighbor_kind_pair%list(1,npairs) = atom_a
       neighbor_kind_pair%list(2,npairs) = atom_b
       neighbor_kind_pair%id_kind(npairs) = id_kind
    END IF
    neighbor_kind_pair%npairs = npairs
  END SUBROUTINE fist_neighbor_add

END MODULE fist_neighbor_list_types
