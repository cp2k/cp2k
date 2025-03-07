!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K+SMEAGOL interface.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************
MODULE smeagol_emtoptions
   #:include 'input_cp2k_smeagol.fypp'
   USE cell_types, ONLY: cell_type, &
                         scaled_to_real
   USE cp_dbcsr_api, ONLY: dbcsr_get_info, &
                           dbcsr_type
#if defined(__SMEAGOL)
   USE global_meshvar, ONLY: smeagolglobal_orbital_BS => orbital_BS
#endif
   USE input_constants, ONLY: smeagol_gridmethod_adaptive, &
                              smeagol_gridmethod_traditional, &
                              smeagol_integraltype_gauss_chebyshev, &
                              smeagol_integraltype_gauss_legendre
   USE kinds, ONLY: default_string_length, &
                    dp
#if defined(__SMEAGOL)
   USE mbfield, ONLY: smeagolglobal_ZeemanBx => ZeemanBx, &
                      smeagolglobal_ZeemanBy => ZeemanBy, &
                      smeagolglobal_ZeemanBz => ZeemanBz
#endif
   USE message_passing, ONLY: mp_para_env_type
#if defined(__SMEAGOL)
   USE mselfenergies, ONLY: SetOptionsSelfEnergies
   USE negfcoop, ONLY: coopinfo
   USE negfmod, ONLY: &
#:for name1, keyword1, val1 in reademtr_negfmod_llist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in reademtr_negfmod_ilist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in reademtr_negfmod_rlist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1, unit1 in reademtr_negfmod_rydberg_plist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1, unit1 in reademtr_negfmod_plist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in readoptsnegf_negfmod_llist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in readoptsnegf_negfmod_ilist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in readoptsnegf_negfmod_rlist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1, unit1 in readoptsnegf_negfmod_explicit_plist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in emtoptions_negfmod_llist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in emtoptions_negfmod_ilist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in emtoptions_negfmod_explicit_ilist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in emtoptions_negfmod_rlist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1 in emtoptions_negfmod_explicit_rlist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
#:for name1, keyword1, val1, unit1 in emtoptions_negfmod_rydberg_plist
      smeagolglobal_${name1}$ => ${name1}$, &
#:endfor
      smeagolglobal_deltabss_bs => deltabss_bs, &
      smeagolglobal_gamma_negf => gamma_negf, &
      smeagolglobal_emforces => emforces, &
      smeagolglobal_emSTT => emSTT, &
      smeagolglobal_emSTTLin => emSTTLin, &
      smeagolglobal_gridmethod => gridmethod, &
      smeagolglobal_integraltype => integraltype, &
      smeagolglobal_ndivxy => ndivxy, &
      smeagolglobal_ndivxyNL => ndivxyNL, &
      smeagolglobal_negf_base_comm => negf_base_comm, &
      smeagolglobal_nebss_bs => nebss_bs, &
      smeagolglobal_nprocs_hs => nprocs_hs
   USE sigma, ONLY: &
#:for name1, keyword1, val1 in emtoptions_sigma_ilist
      smeagolglobal_${name1}$ => ${name1}$
      #:endfor
#endif
      USE smeagol_control_types, ONLY: smeagol_control_type
      USE string_utilities, ONLY: integer_to_string
#include "./base/base_uses.f90"

      IMPLICIT NONE
      PRIVATE

      CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smeagol_emtoptions'

      PUBLIC :: reademtr, ReadOptionsNEGF_DFT, emtrans_options, emtrans_deallocate_global_arrays

   CONTAINS

      SUBROUTINE reademtr(smeagol_control, natoms, gamma_negf)
         TYPE(smeagol_control_type), POINTER                :: smeagol_control
         INTEGER, INTENT(in)                                :: natoms
         LOGICAL, INTENT(in)                                :: gamma_negf

         CHARACTER(LEN=*), PARAMETER :: routineN = 'reademtr'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

#if defined(__SMEAGOL)
         CPASSERT(ASSOCIATED(smeagol_control%aux))

         smeagolglobal_gamma_negf = gamma_negf

         IF (smeagol_control%aux%AtmRVCte <= 0) smeagol_control%aux%AtmRVCte = natoms

         #:for name1, keyword1, val1 in reademtr_negfmod_llist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         #:for name1, keyword1, val1 in reademtr_negfmod_ilist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         #:for name1, keyword1, val1 in reademtr_negfmod_rlist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         #:for name1, keyword1, val1, unit1 in reademtr_negfmod_rydberg_plist
            smeagolglobal_${name1}$ = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%${name1}$
         #:endfor
         #:for name1, keyword1, val1, unit1 in reademtr_negfmod_plist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         coopinfo%ccoop = smeagol_control%aux%COOPCalculate
         coopinfo%nbond = smeagol_control%aux%COOPNumberOfBonds
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
         MARK_USED(smeagol_control)
         MARK_USED(natoms)
         MARK_USED(gamma_negf)
#endif

         CALL timestop(handle)
      END SUBROUTINE reademtr

      SUBROUTINE ReadOptionsNEGF_DFT(smeagol_control, ucell, torqueflag, torquelin)
         TYPE(smeagol_control_type), POINTER                :: smeagol_control
         TYPE(cell_type), POINTER                           :: ucell
         LOGICAL, INTENT(in)                                :: torqueflag, torquelin

         CHARACTER(LEN=*), PARAMETER :: routineN = 'ReadOptionsNEGF_DFT'

         INTEGER                                            :: handle
         REAL(kind=dp), DIMENSION(3)                        :: coord_real, coord_scaled

         CALL timeset(routineN, handle)

#if defined(__SMEAGOL)
         CPASSERT(ASSOCIATED(smeagol_control%aux))

         smeagolglobal_emSTT = torqueflag
         smeagolglobal_emSTTLin = torquelin

         ! In case of the original SIESTA+SMEAGOL, 'TimeReversal' keyword is enabled by default, therefore 'EM.TimeReversal' is also enabled.
         ! In case of this CP2K+SMEAGOL interface, the default value of 'timereversal' variable is .FALSE.
         IF (smeagol_control%aux%timereversal) THEN
            CALL cp_warn(__LOCATION__, &
                         "The SMEAGOL keyword 'EM.TimeReversal' is not supported yet.")
         END IF

         #:for name1, keyword1, val1 in readoptsnegf_negfmod_llist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         #:for name1, keyword1, val1 in readoptsnegf_negfmod_ilist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         #:for name1, keyword1, val1 in readoptsnegf_negfmod_rlist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         ! readoptsnegf_negfmod_plist
         IF (.NOT. smeagol_control%aux%isexplicit_RGlxmin) smeagol_control%aux%RGlxmin = 0.0_dp
         IF (.NOT. smeagol_control%aux%isexplicit_RGlymin) smeagol_control%aux%RGlymin = 0.0_dp
         IF (.NOT. smeagol_control%aux%isexplicit_RGlzmin) smeagol_control%aux%RGlzmin = 0.0_dp
         IF (.NOT. smeagol_control%aux%isexplicit_RGlxmax) THEN
            coord_scaled(:) = (/1.0_dp, 0.0_dp, 0.0_dp/)
            CALL scaled_to_real(coord_real, coord_scaled, ucell)
            smeagol_control%aux%RGlxmax = coord_real(1)
         END IF
         IF (.NOT. smeagol_control%aux%isexplicit_RGlymax) THEN
            coord_scaled(:) = (/0.0_dp, 1.0_dp, 0.0_dp/)
            CALL scaled_to_real(coord_real, coord_scaled, ucell)
            smeagol_control%aux%RGlymax = coord_real(2)
         END IF
         IF (.NOT. smeagol_control%aux%isexplicit_RGlzmax) THEN
            coord_scaled(:) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            CALL scaled_to_real(coord_real, coord_scaled, ucell)
            smeagol_control%aux%RGlzmax = coord_real(3)
         END IF
         #:for name1, keyword1, val1, unit1 in readoptsnegf_negfmod_explicit_plist
            smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
         #:endfor

         ! options to add a Zeeman term to the Hamiltonian
         #:for name1, keyword1, val1, unit1 in readoptsnegf_bfield_rydberg_plist
            smeagolglobal_${name1}$ = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%${name1}$
         #:endfor
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
         MARK_USED(smeagol_control)
         MARK_USED(ucell)
         MARK_USED(torqueflag)
         MARK_USED(torquelin)
         ! local variables
         MARK_USED(coord_real)
         MARK_USED(coord_scaled)
#endif

         CALL timestop(handle)
      END SUBROUTINE ReadOptionsNEGF_DFT

      SUBROUTINE emtrans_options(smeagol_control, matrix_s, para_env, iter, istep, inicoor, iv, delta, nk)
         TYPE(smeagol_control_type), POINTER                :: smeagol_control
         TYPE(dbcsr_type), INTENT(in), POINTER              :: matrix_s
         TYPE(mp_para_env_type), POINTER                    :: para_env
         INTEGER, INTENT(in)                                :: iter, istep, inicoor, iv
         REAL(kind=dp), INTENT(in)                          :: delta
         INTEGER, INTENT(in)                                :: nk

         CHARACTER(LEN=*), PARAMETER :: routineN = 'emtrans_options'

         CHARACTER(len=default_string_length)               :: actual_val_str, expected_val_str
         INTEGER                                            :: GetRhoSingleLeadDefault, handle, i, iatom, n1, nblkcols_total, &
                                                               NParallelK, NParallelKbuf
         INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size

         CALL timeset(routineN, handle)

#if defined(__SMEAGOL)
         CPASSERT(ASSOCIATED(smeagol_control%aux))

         CALL dbcsr_get_info(matrix=matrix_s, nblkcols_total=nblkcols_total, &
                             col_blk_size=col_blk_size, col_blk_offset=col_blk_offset)
         ! number of atomic orbitals.
         ! The name n1 is meaningless, but it is used in the original SIESTA's version of emtrans_options() subroutine
         n1 = SUM(col_blk_size(1:nblkcols_total))

         IF ((iter == 1) .AND. (istep == inicoor) .AND. (iv == 0)) THEN

            IF (smeagol_control%aux%gridmethod == smeagol_gridmethod_traditional) THEN
               smeagolglobal_gridmethod = 'Traditional'
            ELSE IF (smeagol_control%aux%gridmethod == smeagol_gridmethod_adaptive) THEN
               smeagolglobal_gridmethod = 'Adaptivegrid'
            ELSE
               smeagolglobal_gridmethod = 'UNKNOWN'
            END IF

            IF (smeagol_control%aux%integraltype == smeagol_integraltype_gauss_legendre) THEN
               smeagolglobal_integraltype = 'gauss-legendre'
            ELSE IF (smeagol_control%aux%integraltype == smeagol_integraltype_gauss_chebyshev) THEN
               smeagolglobal_integraltype = 'gauss-chebyshev'
            ELSE
               smeagolglobal_integraltype = 'UNKNOWN'
            END IF

            smeagolglobal_negf_base_comm = para_env%get_handle()

            IF (MOD(smeagol_control%aux%ndivisions, 2) == 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "AM.NumberDivisions value must be odd.")
            END IF

            ! It seems that 'sigmatodisk' parameter is specific to the SIESTA interface.
            ! Consider setting smeagol_control%aux%sigmatodisk = .FALSE. and remove smeagol_control%aux%storesigma input keyword.
            IF (smeagol_control%aux%storesigma == 2) THEN
               smeagol_control%aux%sigmatodisk = .TRUE.
            ELSE IF (smeagol_control%aux%storesigma == 1) THEN
               smeagol_control%aux%sigmatodisk = .FALSE.
            END IF

            ! Bound states
            !     options:
            !       bs_add    : true => bound states are added
            !                   false=> normal smeagol
            !       bs_method : 0    => calculate bound states with effective Hamiltonian
            !                   1    => calculate bound states by adding a small
            !                           imaginary part to the selfenergies
            IF (smeagol_control%aux%bs_nmid == 0) smeagol_control%aux%bs_nmid = n1/2

            ! SC: Bad practice: real number comparison
            IF (smeagol_control%aux%SigmaWideBand /= 0.0_dp) smeagol_control%aux%m_svdtolzi = 0.0_dp

            IF (smeagol_control%aux%leadspdos) smeagol_control%aux%leadsdos = .TRUE.
            IF (smeagol_control%aux%curr_distKEne) smeagol_control%aux%curr_distK = .TRUE.
            IF (smeagol_control%aux%curr_distK) smeagol_control%aux%curr_dist = .TRUE.

            IF (smeagolglobal_emSTT .AND. smeagolglobal_emSTTLin .OR. smeagol_control%aux%curr_dist) THEN
               smeagol_control%aux%emldos2 = .TRUE.
               IF (smeagol_control%aux%curr_dist) THEN
                  GetRhoSingleLeadDefault = 3
               ELSE
                  GetRhoSingleLeadDefault = -3
               END IF
            ELSE
               GetRhoSingleLeadDefault = 0
            END IF

            ! current-induced forces
            ! The value of 'smeagol_control%emforces' is set in qs_energies().
            ! Calculation of forces is enabled automatically for certain run_types (energy_force, geo_opt, md) and disabled otherwise.
            IF (smeagol_control%aux%curr_dist) THEN
               smeagol_control%emforces = .TRUE.
            END IF

            IF (.NOT. smeagol_control%aux%isexplicit_nprocs_hs) smeagol_control%aux%nprocs_hs = smeagol_control%aux%nprocs_inverse
            smeagolglobal_nprocs_hs = smeagol_control%aux%nprocs_hs
           IF (.NOT. smeagol_control%aux%isexplicit_GetRhoSingleLead) smeagol_control%aux%GetRhoSingleLead = GetRhoSingleLeadDefault

            IF (smeagol_control%aux%MinChannelIndex < 1) smeagol_control%aux%MinChannelIndex = 1
            IF (smeagol_control%aux%MaxChannelIndex < 1) &
               smeagol_control%aux%MaxChannelIndex = smeagol_control%aux%MinChannelIndex + 4

            IF (smeagolglobal_emSTT .AND. smeagolglobal_emSTTLin .AND. smeagol_control%aux%GetRhoSingleLead /= -3) THEN
               CALL cp_warn(__LOCATION__, &
                            "EM.LDOSLeadsProjection should be set to -3. "// &
                            "If SpinTorque and STLinResp are T, otherwise the output "// &
                            "results for the spin transfer torque are incorrect.")
            END IF

            ! NParallelK
            NParallelK = smeagol_control%aux%NParallelK
            IF (MOD(para_env%num_pe, NParallelK) /= 0) then
               CALL cp_warn(__LOCATION__, &
                            "EM.ParallelOverKNum must be a divisor of the total number of "// &
                            "MPI processes used in a run; resetting the value of NParallelK.")
               NParallelKbuf = NParallelK
               DO NParallelK = NParallelKbuf, 1, -1
                  IF (MOD(para_env%num_pe, NParallelK) == 0) EXIT
               END DO
            END IF

            IF (NParallelK > para_env%num_pe) THEN
               CALL cp_warn(__LOCATION__, &
                            "EM.ParallelOverKNum can not be larger than the total number of "// &
                            "MPI processes used in a run; resetting the value of NParallelK.")
               NParallelK = -1
            END IF

            IF (NParallelK > nk) THEN
               CALL cp_warn(__LOCATION__, &
                            "EM.ParallelOverKNum can not be larger than the total number of "// &
                            "k-points used in a run; resetting the value of NParallelK.")
               NParallelK = -1
            END IF

            IF (NParallelK == -1) THEN
               DO NParallelK = nk, 1, -1
                  IF (MOD(para_env%num_pe, NParallelK) == 0) EXIT
               END DO
            END IF

            smeagol_control%aux%NParallelK = NParallelK

            IF (smeagol_control%aux%empdosk) smeagol_control%aux%empdos = .TRUE.
            IF (smeagol_control%aux%emldos2) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%TransmissionChannels) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%TransmissionMatrix) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%curr_dist) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%empdos) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%m_skipsvd < 2.0_dp) smeagol_control%aux%m_skipsvd = 10.0_dp

            IF (smeagol_control%aux%CallImpuritySolver) smeagol_control%aux%ComputeImpurityGfMatsubara = .TRUE.

            #:for name1, keyword1, val1 in emtoptions_negfmod_llist
               smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
            #:endfor
            smeagolglobal_emforces = smeagol_control%emforces

            ! emtoptions_local_ilist, emtoptions_local_explicit_ilist
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NxLeft) smeagol_control%aux%Sigma_NxLeft = smeagol_control%aux%Sigma_Nx
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NyLeft) smeagol_control%aux%Sigma_NyLeft = smeagol_control%aux%Sigma_Ny
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NxRight) smeagol_control%aux%Sigma_NxRight = smeagol_control%aux%Sigma_Nx
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NyRight) smeagol_control%aux%Sigma_NyRight = smeagol_control%aux%Sigma_Ny

            smeagolglobal_ndivxy(1) = smeagol_control%aux%Sigma_Nx
            smeagolglobal_ndivxy(2) = smeagol_control%aux%Sigma_Ny
            smeagolglobal_ndivxyNL(1, 1) = smeagol_control%aux%Sigma_NxLeft
            smeagolglobal_ndivxyNL(1, 2) = smeagol_control%aux%Sigma_NyLeft
            smeagolglobal_ndivxyNL(2, 1) = smeagol_control%aux%Sigma_NxRight
            smeagolglobal_ndivxyNL(2, 2) = smeagol_control%aux%Sigma_NyRight

            CALL SetOptionsSelfEnergies(smeagolglobal_ndivxyNL, 2)

            #:for name1, keyword1, val1 in emtoptions_negfmod_ilist
               smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
            #:endfor

            #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_ilist
               smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
            #:endfor

            #:for name1, keyword1, val1 in emtoptions_sigma_ilist
               smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
            #:endfor

            #:for name1, keyword1, val1 in emtoptions_negfmod_rlist
               smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
            #:endfor

            ! +++ bound states
            IF (.NOT. smeagol_control%aux%isexplicit_deltamin) smeagol_control%aux%deltamin = delta

            #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_rlist
               smeagolglobal_${name1}$ = smeagol_control%aux%${name1}$
            #:endfor

            #:for name1, keyword1, val1 in emtoptions_negfmod_rydberg_plist
               smeagolglobal_${name1}$ = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%${name1}$
            #:endfor

            ! At present read_options_ImpuritySolver() is a dummy libsmeagol.a subroutine,
            ! so we do not call it
            !IF (smeagol_control%ComputeImpurityGfMatsubara) THEN
            !   CALL read_options_ImpuritySolver(CallImpuritySolver,n1)
            !END IF
         END IF

         ! *** bound states (re-implemented ReadSpeciesBS() and ReadBSSubSystemBoundaries() subroutines)

         ! bound- state-related global allocatable arrays
         ALLOCATE (smeagolglobal_deltabss_bs(smeagol_control%aux%nbss))
         ALLOCATE (smeagolglobal_nebss_bs(smeagol_control%aux%nbss, 2))
         ALLOCATE (smeagolglobal_orbital_BS(n1))

         ! BS.SubSystemsDelta
         IF (ALLOCATED(smeagol_control%aux%deltabss_bs)) THEN
            smeagolglobal_deltabss_bs(:) = smeagol_control%aux%deltabss_bs(:)
         ELSE
            smeagolglobal_deltabss_bs(:) = smeagol_control%aux%deltamin
         END IF

         ! BS.SubSystemsBoundaries
         IF (ALLOCATED(smeagol_control%aux%nebss_bs)) THEN
            IF (MAXVAL(smeagol_control%aux%nebss_bs) > n1) THEN
               CALL integer_to_string(MAXVAL(smeagol_control%aux%nebss_bs), actual_val_str)
               CALL integer_to_string(n1, expected_val_str)
               CALL cp_abort(__LOCATION__, &
                             "The largest index in BS.SubSystemsBoundaries section ("//TRIM(actual_val_str)// &
                             ") exceeds the number of atomic orbitals ("//TRIM(expected_val_str)//").")
            END IF
            smeagolglobal_nebss_bs(:, :) = smeagol_control%aux%nebss_bs(:, :)
         ELSE
            smeagolglobal_nebss_bs(1, 1) = 1
            smeagolglobal_nebss_bs(1, 2) = n1/smeagol_control%aux%nbss
            DO i = 2, smeagol_control%aux%nbss
               smeagolglobal_nebss_bs(i, 1) = smeagolglobal_nebss_bs(i - 1, 2) + 1
               smeagolglobal_nebss_bs(i, 2) = i*n1/smeagol_control%aux%nbss
            END DO
            smeagolglobal_nebss_bs(smeagol_control%aux%nbss, 2) = n1
         END IF

         ! AM.AtomListBS
         IF (ALLOCATED(smeagol_control%aux%atomlist_bs)) THEN
            IF (MAXVAL(smeagol_control%aux%atomlist_bs) > nblkcols_total) THEN
               CALL integer_to_string(MAXVAL(smeagol_control%aux%atomlist_bs), actual_val_str)
               CALL integer_to_string(nblkcols_total, expected_val_str)
               CALL cp_abort(__LOCATION__, &
                             "The largest atomic index in AM.AtomListBS keyword ("//TRIM(actual_val_str)// &
                             ") exceeds the number of atoms ("//TRIM(expected_val_str)//") in the extended molecule.")
            END IF

            smeagolglobal_orbital_BS(:) = .FALSE.
            DO i = 1, SIZE(smeagol_control%aux%atomlist_bs)
               iatom = smeagol_control%aux%atomlist_bs(i)
               smeagolglobal_orbital_BS(col_blk_offset(iatom):col_blk_offset(iatom) + col_blk_size(iatom) - 1) = .TRUE.
            END DO
         ELSE
            smeagolglobal_orbital_BS(:) = .TRUE.
         END IF
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
         MARK_USED(smeagol_control)
         MARK_USED(matrix_s)
         MARK_USED(para_env)
         MARK_USED(iter)
         MARK_USED(istep)
         MARK_USED(inicoor)
         MARK_USED(iv)
         MARK_USED(delta)
         MARK_USED(nk)
         ! local variables
         MARK_USED(actual_val_str)
         MARK_USED(col_blk_offset)
         MARK_USED(col_blk_size)
         MARK_USED(expected_val_str)
         MARK_USED(GetRhoSingleLeadDefault)
         MARK_USED(i)
         MARK_USED(iatom)
         MARK_USED(n1)
         MARK_USED(nblkcols_total)
         MARK_USED(NParallelK)
         MARK_USED(NParallelKbuf)
#endif

         CALL timestop(handle)
      END SUBROUTINE emtrans_options

      SUBROUTINE emtrans_deallocate_global_arrays()

#if defined(__SMEAGOL)
         IF (ALLOCATED(smeagolglobal_deltabss_bs)) DEALLOCATE (smeagolglobal_deltabss_bs)
         IF (ALLOCATED(smeagolglobal_nebss_bs)) DEALLOCATE (smeagolglobal_nebss_bs)
         IF (ALLOCATED(smeagolglobal_orbital_BS)) DEALLOCATE (smeagolglobal_orbital_BS)
#endif

      END SUBROUTINE emtrans_deallocate_global_arrays

   END MODULE smeagol_emtoptions

