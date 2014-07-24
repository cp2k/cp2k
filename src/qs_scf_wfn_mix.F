!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for GPW/GAPW
!> \par History
!>      Started as a copy from the relevant part of qs_scf
!> \author Joost VandeVondele (10.2003)
! *****************************************************************************
MODULE qs_scf_wfn_mix
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_mo_methods,                   ONLY: calculate_orthonormality
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             duplicate_mo_set,&
                                             mo_set_p_type,&
                                             write_mo_set
  USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                             special_diag_method_nr
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_wfn_mix'
  PUBLIC :: wfn_mix

CONTAINS

  ! *****************************************************************************
  !> \brief writes a new 'mixed' set of mos to restart file, without touching the current MOs
  ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param mos ...
!> \param particle_set ...
!> \param dft_section ...
!> \param atomic_kind_set ...
!> \param unoccupied_orbs ...
!> \param scf_env ...
!> \param matrix_s ...
!> \param output_unit ...
!> \param marked_states ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE wfn_mix(mos, particle_set, dft_section, atomic_kind_set, &
                    unoccupied_orbs, scf_env, matrix_s, output_unit, marked_states,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    INTEGER                                  :: output_unit
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: marked_states
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfn_mix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i_rep, ispin, mark_ind, mark_number, n_rep, &
      orig_mo_index, orig_spin_index, result_mo_index, result_spin_index
    LOGICAL                                  :: explicit, failure, &
                                                orig_is_virtual, overwrite_mos
    REAL(KIND=dp)                            :: orig_scale, orthonormality, &
                                                result_scale
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_vector
    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_y
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(section_vals_type), POINTER         :: update_section, &
                                                wfn_mix_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    wfn_mix_section => section_vals_get_subs_vals(dft_section,"PRINT%WFN_MIX",error=error)
    CALL section_vals_get(wfn_mix_section,explicit=explicit, error=error)

    ! only perform action if explicitly required
    IF (explicit) THEN

       IF (output_unit>0) THEN
          WRITE(output_unit,'()')
          WRITE(output_unit,'(T2,A)') "Performing wfn mixing"
          WRITE(output_unit,'(T2,A)') "====================="
       ENDIF

       ALLOCATE(mos_new(SIZE(mos)))
       DO ispin=1,SIZE(mos)
          CALL duplicate_mo_set(mos_new(ispin)%mo_set,mos(ispin)%mo_set,error)
       ENDDO

       ! a single vector matrix structure
       NULLIFY(fm_struct_vector)
       CALL cp_fm_struct_create(fm_struct_vector,template_fmstruct=mos(1)%mo_set%mo_coeff%matrix_struct, &
                                ncol_global=1, error=error)
       CALL cp_fm_create(matrix_x,fm_struct_vector,name="x",error=error)
       CALL cp_fm_create(matrix_y,fm_struct_vector,name="y",error=error)
       CALL cp_fm_struct_release(fm_struct_vector,error)

       update_section=>section_vals_get_subs_vals(wfn_mix_section,"UPDATE",error=error)
       CALL section_vals_get(update_section,n_repetition=n_rep,error=error)
       CALL section_vals_get(update_section,explicit=explicit, error=error)
       IF (.NOT. explicit) n_rep=0

       DO i_rep=1,n_rep

          CALL section_vals_val_get(update_section,"RESULT_MO_INDEX",i_rep_section=i_rep,i_val=result_mo_index,error=error)
          CALL section_vals_val_get(update_section,"RESULT_MARKED_STATE",i_rep_section=i_rep,i_val=mark_number,error=error)
          CALL section_vals_val_get(update_section,"RESULT_SPIN_INDEX",i_rep_section=i_rep,i_val=result_spin_index,error=error)
          CALL section_vals_val_get(update_section,"RESULT_SCALE",i_rep_section=i_rep,r_val=result_scale,error=error)

          mark_ind=1
          IF(mark_number.GT.0)result_mo_index=marked_states(mark_number,result_spin_index,mark_ind)

          CALL section_vals_val_get(update_section,"ORIG_MO_INDEX",i_rep_section=i_rep,i_val=orig_mo_index,error=error)
          CALL section_vals_val_get(update_section,"ORIG_MARKED_STATE",i_rep_section=i_rep,i_val=mark_number,error=error)
          CALL section_vals_val_get(update_section,"ORIG_SPIN_INDEX",i_rep_section=i_rep,i_val=orig_spin_index,error=error)
          CALL section_vals_val_get(update_section,"ORIG_SCALE",i_rep_section=i_rep,r_val=orig_scale,error=error)
          CALL section_vals_val_get(update_section,"ORIG_IS_VIRTUAL",i_rep_section=i_rep,l_val=orig_is_virtual,error=error)

          IF(orig_is_virtual)mark_ind=2
          IF(mark_number.GT.0)orig_mo_index=marked_states(mark_number,orig_spin_index,mark_ind)

          CALL section_vals_val_get(wfn_mix_section,"OVERWRITE_MOS",l_val=overwrite_mos,error=error)

          ! first get a copy of the proper orig
          IF (.NOT. ORIG_IS_VIRTUAL) THEN
             CALL cp_fm_to_fm(mos(orig_spin_index)%mo_set%mo_coeff,matrix_x,1,mos(orig_spin_index)%mo_set%nmo-orig_mo_index+1,1)
          ELSE
             CALL cp_fm_to_fm(unoccupied_orbs(orig_spin_index)%matrix,matrix_x,1,orig_mo_index,1)
          ENDIF

          ! now get a copy of the target
          CALL cp_fm_to_fm(mos_new(result_spin_index)%mo_set%mo_coeff,matrix_y, &
                           1,mos_new(result_spin_index)%mo_set%nmo-result_mo_index+1,1)

          ! properly mix
          CALL cp_fm_scale_and_add(result_scale,matrix_y,orig_scale,matrix_x,error)

          ! and copy back in the result mos
          CALL cp_fm_to_fm(matrix_y,mos_new(result_spin_index)%mo_set%mo_coeff, &
                           1,1,mos_new(result_spin_index)%mo_set%nmo-result_mo_index+1)

       ENDDO

       CALL cp_fm_release(matrix_x,error)
       CALL cp_fm_release(matrix_y,error)

       IF(scf_env%method==special_diag_method_nr) THEN
         CALL calculate_orthonormality(orthonormality,mos,error=error)
       ELSE
         CALL calculate_orthonormality(orthonormality,mos,matrix_s(1)%matrix,error=error)
       END IF

       IF (output_unit>0) THEN
          WRITE(output_unit,'()')
          WRITE(output_unit,'(T2,A,T61,E20.4)')  &
                 "Maximum deviation from MO S-orthonormality",orthonormality
          WRITE(output_unit,'(T2,A)') "Writing new MOs to file"
       ENDIF

       ! *** Write WaveFunction restart file ***

       DO ispin=1,SIZE(mos_new)
          IF(overwrite_mos)THEN
             CALL cp_fm_to_fm(mos_new(ispin)%mo_set%mo_coeff,mos(ispin)%mo_set%mo_coeff,error)
             IF(mos_new(1)%mo_set%use_mo_coeff_b)&
                CALL copy_fm_to_dbcsr(mos_new(ispin)%mo_set%mo_coeff,mos_new(ispin)%mo_set%mo_coeff_b,&
                                      error=error)
          END IF 
          IF(mos(1)%mo_set%use_mo_coeff_b)&
             CALL copy_fm_to_dbcsr(mos_new(ispin)%mo_set%mo_coeff,mos(ispin)%mo_set%mo_coeff_b,&
                                   error=error)
       END DO
       CALL write_mo_set(mos_new,particle_set,dft_section=dft_section, &
                         atomic_kind_set=atomic_kind_set,error=error)

       DO ispin=1,SIZE(mos_new)
          CALL deallocate_mo_set(mos_new(ispin)%mo_set,error)
       ENDDO
       DEALLOCATE(mos_new)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE wfn_mix

END MODULE qs_scf_wfn_mix
