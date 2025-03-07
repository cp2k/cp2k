!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input control types for NEGF/SMEAGOL transport calculations.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************

MODULE smeagol_control_types
   USE cp_units, ONLY: cp_unit_from_cp2k
   USE input_constants, ONLY: smeagol_bulklead_leftright, &
                              smeagol_gridmethod_traditional, &
                              smeagol_integraltype_gauss_legendre, &
                              smeagol_runtype_bulktransport, &
                              smeagol_runtype_emtransport
   USE input_section_types, ONLY: section_vals_get, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: default_string_length, &
                    dp
   USE physcon, ONLY: kelvin
   USE string_utilities, ONLY: integer_to_string
   USE util, ONLY: sort
#include "./base/base_uses.f90"
   #:include 'input_cp2k_smeagol.fypp'

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smeagol_control_types'

   PUBLIC :: smeagol_control_type
   PUBLIC :: smeagol_control_create, smeagol_control_release, read_smeagol_control

! **************************************************************************************************
!> \brief SMEAGOL-related input parameters
! **************************************************************************************************
   TYPE smeagol_control_type
      LOGICAL                                            :: smeagol_enabled = .FALSE.

      !> type of calculation
      INTEGER                                            :: run_type = smeagol_runtype_bulktransport

      !> regression test mode
      LOGICAL                                            :: do_regtest = .FALSE.

      !> current-induced forces. It is set automatically based on GLOBAL/run_type
      LOGICAL                                            :: emforces = .FALSE.

      !> scale factor to convert from CP2K (Hartree) to SMEAGOL (Rydberg) default energy unit
      REAL(kind=dp)                                      :: to_smeagol_energy_units = 2.0_dp

      !> number of cell images along i and j cell vectors
      INTEGER, DIMENSION(2)                              :: n_cell_images = (/1, 1/)

      !> what lead (bulk transport calculation)
      INTEGER                                            :: lead_label = smeagol_bulklead_leftright

      !> The length of the SMEAGOL project name is limited by 20 characters (hardcoded in SMEAGOL)
      CHARACTER(len=20)                                  :: project_name = "PROJECT_NAME"

      TYPE(smeagol_aux_control_type), POINTER            :: aux => NULL()

   END TYPE smeagol_control_type

! **************************************************************************************************
!> \brief SMEAGOL-related auxiliary input parameters. They remain unallocated when
!>        SMEAGOL support is disabled.
! **************************************************************************************************
   TYPE smeagol_aux_control_type
      !> MD or GEO_OPT iteration. In contrast with other components of smeagol_control_type that are read from the input file,
      !> these variables are initialised at the first NEGF iteration.
      !>
      !> index of MD or GEO_OPT iteration level.
      !> 0 if there is neither MD nor GEO_OPT iteration level (e.g. single-point energy calculation).
      !> -1 if it is not initialised by run_smeagol_emtrans() subroutine
      INTEGER                                            :: md_iter_level = -1
      !> The starting step value for MD / GEO_OPT iterations. The default value 0 can be overrited via STEP_START_VAL input keyword.
      INTEGER                                            :: md_first_step = 0

      !> BS.SubSystemsDelta(1:BS.Subsystems)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: deltabss_bs
      !> BS.SubSystemsBoundaries(1:BS.Subsystems, 1:2)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nebss_bs
      !>
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_bs

      REAL(kind=dp)                                      :: temperature = 300.0_dp/kelvin

      ! reademtr()
      #:for name1, keyword1, val1 in reademtr_local_llist
         LOGICAL                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_local_explicit_plist
         LOGICAL                                            :: isexplicit_${name1}$ = .FALSE.
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfmod_llist
         LOGICAL                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfcoop_llist
         LOGICAL                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in reademtr_local_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfmod_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfcoop_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in reademtr_local_rlist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfmod_rlist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_local_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_local_explicit_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_negfmod_rydberg_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_negfmod_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      ! ReadOptionsNEGF_DFT()
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_llist
         LOGICAL                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in readoptsnegf_negfmod_explicit_plist
         LOGICAL                                            :: isexplicit_${name1}$ = .FALSE.
      #:endfor

      #:for name1, keyword1, val1 in readoptsnegf_negfmod_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in readoptsnegf_negfmod_rlist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in readoptsnegf_negfmod_explicit_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in readoptsnegf_bfield_rydberg_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      ! emtrans_options()
      INTEGER                                            :: gridmethod = smeagol_gridmethod_traditional
      INTEGER                                            :: integraltype = smeagol_integraltype_gauss_legendre

      #:for name1, keyword1, val1 in emtoptions_negfmod_llist
         LOGICAL                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in emtoptions_local_explicit_ilist
         LOGICAL                                            :: isexplicit_${name1}$ = .FALSE.
      #:endfor

      #:for name1, keyword1, val1, unit1 in emtoptions_negfmod_explicit_ilist
         LOGICAL                                            :: isexplicit_${name1}$ = .FALSE.
      #:endfor

      #:for name1, keyword1, val1, unit1 in emtoptions_negfmod_explicit_rlist
         LOGICAL                                            :: isexplicit_${name1}$ = .FALSE.
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_local_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_local_explicit_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_sigma_ilist
         INTEGER                                            :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_rlist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_rlist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      #:for name1, keyword1, val1, unit1 in emtoptions_negfmod_rydberg_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor

      ! interface options
      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor
      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_explicit_plist
         REAL(kind=dp)                                      :: ${name1}$ = ${val1}$
      #:endfor
      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_explicit_plist
         LOGICAL                                            :: isexplicit_${name1}$ = .FALSE.
      #:endfor

   END TYPE smeagol_aux_control_type
CONTAINS

! **************************************************************************************************
!> \brief allocate control options for SMEAGOL calculation
!> \param smeagol_control an object to create
! **************************************************************************************************
   SUBROUTINE smeagol_control_create(smeagol_control)
      TYPE(smeagol_control_type), POINTER                :: smeagol_control

      CHARACTER(len=*), PARAMETER :: routineN = 'smeagol_control_create'

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(smeagol_control))
      CALL timeset(routineN, handle)

      ALLOCATE (smeagol_control)
#if defined(__SMEAGOL)
      ALLOCATE (smeagol_control%aux)
#else
      NULLIFY (smeagol_control%aux)
#endif

      CALL timestop(handle)
   END SUBROUTINE smeagol_control_create

! **************************************************************************************************
!> \brief release SMEAGOL control object
!> \param smeagol_control an object to release
! **************************************************************************************************
   SUBROUTINE smeagol_control_release(smeagol_control)
      TYPE(smeagol_control_type), POINTER                :: smeagol_control

      CHARACTER(len=*), PARAMETER :: routineN = 'smeagol_control_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(smeagol_control)) THEN
         IF (ASSOCIATED(smeagol_control%aux)) THEN
            IF (ALLOCATED(smeagol_control%aux%nebss_bs)) DEALLOCATE (smeagol_control%aux%nebss_bs)
            IF (ALLOCATED(smeagol_control%aux%deltabss_bs)) DEALLOCATE (smeagol_control%aux%deltabss_bs)
            IF (ALLOCATED(smeagol_control%aux%atomlist_bs)) DEALLOCATE (smeagol_control%aux%atomlist_bs)

            DEALLOCATE (smeagol_control%aux)
         END IF
         DEALLOCATE (smeagol_control)
      END IF

      CALL timestop(handle)
   END SUBROUTINE smeagol_control_release

! **************************************************************************************************
!> \brief Read SMEAGOL-related input parameters.
!> \param smeagol_control SMEAGOL control parameters
!> \param smeagol_section SMEAGOL input section
! **************************************************************************************************
   SUBROUTINE read_smeagol_control(smeagol_control, smeagol_section)
      TYPE(smeagol_control_type), POINTER                :: smeagol_control
      TYPE(section_vals_type), POINTER                   :: smeagol_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_smeagol_control'

      CHARACTER(len=default_string_length)               :: project_name
      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: n_cell_images_ptr

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(smeagol_section, "_SECTION_PARAMETERS_", l_val=smeagol_control%smeagol_enabled)

      ! SMEAGOL project name
      CALL section_vals_val_get(smeagol_section, "PROJECT_NAME", c_val=project_name)
      smeagol_control%project_name = project_name
      IF (smeagol_control%smeagol_enabled .AND. smeagol_control%project_name /= project_name) THEN
         CALL cp_warn(__LOCATION__, &
                      "SMEAGOL limits the length of the project name by 20 characters. "// &
                      "The project name is truncated to '"//TRIM(smeagol_control%project_name)//"'.")
      END IF

      ! reademtr() enum keywords
      CALL section_vals_val_get(smeagol_section, "RUN_TYPE", i_val=smeagol_control%run_type)

#if !defined(__SMEAGOL)
      IF (smeagol_control%run_type == smeagol_runtype_emtransport) THEN
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support. SMEAGOL RUN_TYPE EMTransport is not available")
      END IF
#endif

      CALL section_vals_val_get(smeagol_section, "REGRESSION_TEST", l_val=smeagol_control%do_regtest)

      CALL section_vals_val_get(smeagol_section, "BulkLead", i_val=smeagol_control%lead_label)

      ! NOTE: keyword NIMAGES_IJ seems to be identical to ignored SMEAGOL keywords BulkTransvCellSizeX and BulkTransvCellSizeY
      NULLIFY (n_cell_images_ptr)
      CALL section_vals_val_get(smeagol_section, "NIMAGES_IJ", i_vals=n_cell_images_ptr)
      smeagol_control%n_cell_images(1:2) = n_cell_images_ptr(1:2)

      ! Hartree > Rydberg scaling factor
      smeagol_control%to_smeagol_energy_units = cp_unit_from_cp2k(1.0_dp, "RY")

      IF (ASSOCIATED(smeagol_control%aux)) CALL read_smeagol_aux_control(smeagol_control%aux, smeagol_section)

      CALL timestop(handle)
   END SUBROUTINE read_smeagol_control

! **************************************************************************************************
!> \brief Read SMEAGOL-related auxiliary input parameters.
!> \param smeagol_control SMEAGOL control parameters
!> \param smeagol_section SMEAGOL input section
! **************************************************************************************************
   SUBROUTINE read_smeagol_aux_control(smeagol_control, smeagol_section)
      TYPE(smeagol_aux_control_type), POINTER            :: smeagol_control
      TYPE(section_vals_type), POINTER                   :: smeagol_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_smeagol_aux_control'

      CHARACTER(len=default_string_length)               :: nvals_str
      INTEGER                                            :: handle, i, n, nrep, n_unique
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ivec, indices
      INTEGER, DIMENSION(:), POINTER                     :: iptr
      LOGICAL                                            :: is_explicit, do_abort
      REAL(kind=dp), DIMENSION(:), POINTER               :: rptr
      TYPE(section_vals_type), POINTER                   :: subsection

      CALL timeset(routineN, handle)

      smeagol_control%md_iter_level = -1
      smeagol_control%md_first_step = 0

      CALL section_vals_val_get(smeagol_section, "TEMPERATURE", r_val=smeagol_control%temperature)

      ! reademtr() logical keywords
      ! The following keywords (variables) that are part of the reademtr_local_llist list are read but unused
      ! CB.WriteComplexBands (WriteComplexBands = .FALSE.) unimplemented.
      ! Presumably, some of these keywords should be removed.
      #:for name1, keyword1, val1 in reademtr_local_llist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", l_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1 in reademtr_negfmod_llist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", l_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1 in reademtr_negfcoop_llist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", l_val=smeagol_control%${name1}$)
      #:endfor

      ! reademtr() integer keywords
      #:for name1, keyword1, val1 in reademtr_local_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1 in reademtr_negfmod_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1 in reademtr_negfcoop_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor

      ! reademtr() real-valued keywords
      #:for name1, keyword1, val1 in reademtr_local_rlist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1 in reademtr_negfmod_rlist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor

      ! reademtr() physical-valued keywords
      #:for name1, keyword1, val1, unit1 in reademtr_local_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1, unit1 in reademtr_local_explicit_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", &
                                   r_val=smeagol_control%${name1}$, &
                                   explicit=smeagol_control%isexplicit_${name1}$)
      #:endfor
      #:for name1, keyword1, val1, unit1 in reademtr_negfmod_rydberg_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1, unit1 in reademtr_negfmod_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor

      ! ReadOptionsNEGF_DFT() logical keywords
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_llist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", l_val=smeagol_control%${name1}$)
      #:endfor

      ! ReadOptionsNEGF_DFT() integer keywords
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor

      ! ReadOptionsNEGF_DFT() real-valued keywords
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_rlist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor

      ! ReadOptionsNEGF_DFT() physical-valued keywords
      #:for name1, keyword1, val1, unit1 in readoptsnegf_negfmod_explicit_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", &
                                   r_val=smeagol_control%${name1}$, &
                                   explicit=smeagol_control%isexplicit_${name1}$)
      #:endfor

      #:for name1, keyword1, val1, unit1 in readoptsnegf_bfield_rydberg_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor

      ! emtrans_options() enum keywords
      CALL section_vals_val_get(smeagol_section, "EnergyGridType", i_val=smeagol_control%gridmethod)
      CALL section_vals_val_get(smeagol_section, "TypeOfIntegral", i_val=smeagol_control%integraltype)

      ! emtrans_options() logical keywords
      #:for name1, keyword1, val1 in emtoptions_negfmod_llist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", l_val=smeagol_control%${name1}$)
      #:endfor

      ! emtrans_options() integer keywords
      #:for name1, keyword1, val1 in emtoptions_local_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_local_explicit_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", &
                                   i_val=smeagol_control%${name1}$, &
                                   explicit=smeagol_control%isexplicit_${name1}$)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", &
                                   i_val=smeagol_control%${name1}$, &
                                   explicit=smeagol_control%isexplicit_${name1}$)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_sigma_ilist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", i_val=smeagol_control%${name1}$)
      #:endfor

      ! emtrans_options() real-valued keywords
      #:for name1, keyword1, val1 in emtoptions_negfmod_rlist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_rlist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", &
                                   r_val=smeagol_control%${name1}$, &
                                   explicit=smeagol_control%isexplicit_${name1}$)
      #:endfor

      ! emtrans_options() physical-valued keywords
      #:for name1, keyword1, val1, unit1 in emtoptions_negfmod_rydberg_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor

      ! interface physical-valued keywords
      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", r_val=smeagol_control%${name1}$)
      #:endfor
      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_explicit_plist
         CALL section_vals_val_get(smeagol_section, "${keyword1}$", &
                                   r_val=smeagol_control%${name1}$, &
                                   explicit=smeagol_control%isexplicit_${name1}$)
      #:endfor

      ! Bound states
      IF (smeagol_control%nbss <= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "BS.Subsystems keyword should have a positive integer value.")
      END IF

      subsection => section_vals_get_subs_vals(smeagol_section, "BS.SubSystemsBoundaries")
      CALL section_vals_get(subsection, explicit=is_explicit)
      IF (is_explicit) THEN
         CALL section_vals_val_get(subsection, "_DEFAULT_KEYWORD_", n_rep_val=nrep, explicit=is_explicit)
         IF (.NOT. is_explicit) nrep = 0

         do_abort = .FALSE.
         IF (smeagol_control%nbss == nrep) THEN
            IF (is_explicit) THEN
               ALLOCATE (smeagol_control%nebss_bs(nrep, 2))
               DO i = 1, nrep
                  CALL section_vals_val_get(subsection, "_DEFAULT_KEYWORD_", i_rep_val=i, i_vals=iptr)
                  IF (SIZE(iptr) == 2) THEN
                     IF (MINVAL(iptr) > 0 .AND. iptr(1) <= iptr(2)) THEN
                        smeagol_control%nebss_bs(i, 1:2) = iptr(1:2)
                     ELSE
                        do_abort = .TRUE.
                        EXIT
                     END IF
                  ELSE
                     do_abort = .TRUE.
                     EXIT
                  END IF
               END DO
            END IF
         ELSE
            do_abort = .TRUE.
         END IF

         IF (do_abort) THEN
            CALL integer_to_string(smeagol_control%nbss, nvals_str)
            CALL cp_abort(__LOCATION__, &
                          "BS.SubSystemsBoundaries section is expected to have BS.Subsystems ("//TRIM(nvals_str)// &
                          ") lines with two positive integer numbers on each line.")
         END IF

      END IF

      CALL section_vals_val_get(smeagol_section, "BS.SubSystemsDelta", explicit=is_explicit)
      IF (is_explicit) THEN
         CALL section_vals_val_get(smeagol_section, "BS.SubSystemsDelta", r_vals=rptr)
         IF (smeagol_control%nbss /= SIZE(rptr)) THEN ! do negative deltas make sense ? .OR. MINVAL(rptr) < 0
            CALL integer_to_string(smeagol_control%nbss, nvals_str)
            CALL cp_abort(__LOCATION__, &
                          "BS.SubSystemsDelta keyword is expected to have BS.Subsystems ("//TRIM(nvals_str)//") real numbers.")
         END IF

         n = SIZE(rptr)
         ALLOCATE (smeagol_control%deltabss_bs(n))
         smeagol_control%deltabss_bs(1:n) = rptr(1:n)
         !ELSE allocate and set deltabss_bs() to deltamin
      END IF

      CALL section_vals_val_get(smeagol_section, "AM.AtomListBS", explicit=is_explicit)
      IF (is_explicit) THEN
         CALL section_vals_val_get(smeagol_section, "AM.AtomListBS", i_vals=iptr)

         ! subsys is not available at this stage, so we cannot check that MINVAL(iptr) <= natoms so far
         IF (SIZE(iptr) == 0 .OR. MINVAL(iptr) <= 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "All atomic indices in AM.AtomListBS should be positive integer numbers.")
         END IF

         n = SIZE(iptr)
         ALLOCATE (ivec(n), indices(n))
         ivec(1:n) = iptr(1:n)
         CALL sort(ivec, n, indices)

         n_unique = 1
         DO i = 2, n
            IF (ivec(i) > ivec(i - 1)) n_unique = n_unique + 1
         END DO

         ALLOCATE (smeagol_control%atomlist_bs(n_unique))
         n_unique = 1
         smeagol_control%atomlist_bs(1) = ivec(1)

         DO i = 2, n
            IF (ivec(i) > ivec(i - 1)) THEN
               n_unique = n_unique + 1
               smeagol_control%atomlist_bs(n_unique) = ivec(i)
            END IF
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE read_smeagol_aux_control
END MODULE smeagol_control_types
