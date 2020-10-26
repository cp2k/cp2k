!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to perform the RTP in the velocity gauge
! **************************************************************************************************

MODULE rt_propagation_velocity_gauge
   USE ai_moments,                      ONLY: cossin
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_print,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry
   USE external_potential_types,        ONLY: gth_potential_p_type,&
                                              gth_potential_type,&
                                              sgp_potential_p_type,&
                                              sgp_potential_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: one,&
                                              twopi,&
                                              zero
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: c_light_au
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE sap_kind_types,                  ONLY: clist_type,&
                                              release_sap_int,&
                                              sap_int_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_velocity_gauge'

   PUBLIC :: velocity_gauge_ks_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE velocity_gauge_ks_matrix(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'velocity_gauge_ks_matrix'

      INTEGER                                            :: handle, idir, ispin
      LOGICAL                                            :: do_hfx
      REAL(KIND=dp), DIMENSION(3)                        :: vec_pot
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_im, matrix_s, &
                                                            momentum, nl_term
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, matrix_s, sab_orb, sab_all, matrix_ks, cell)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      sab_all=sab_all, &
                      sab_orb=sab_orb, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      cell=cell)

      !get vector potential
      CALL get_vector_potential(dft_control, cell, vec_pot)

      IF (dft_control%qs_control%gapw) THEN

         NULLIFY (input, matrix_ks_im)
         CALL get_qs_env(qs_env, &
                         input=input, &
                         matrix_ks_im=matrix_ks_im)

         ! allocate and build matrices for linear momentum term
         NULLIFY (momentum)
         CALL dbcsr_allocate_matrix_set(momentum, 3)
         DO idir = 1, 3
            CALL dbcsr_init_p(momentum(idir)%matrix)
            CALL dbcsr_create(momentum(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(momentum(idir)%matrix, sab_orb)
            CALL dbcsr_set(momentum(idir)%matrix, zero)
         ENDDO
         CALL build_lin_mom_matrix(qs_env, momentum)

         ! set imaginary part of KS matrix to zero
         ! if not done for hfx in subroutine hfx_ks_matrix()
         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
         CALL section_vals_get(hfx_sections, explicit=do_hfx)
         IF (.NOT. do_hfx) THEN
            DO ispin = 1, dft_control%nspins
               CALL dbcsr_set(matrix_ks_im(ispin)%matrix, zero)
            ENDDO
         ENDIF

         ! add linear term in vector potential to imaginary part of KS-matrix
         DO ispin = 1, dft_control%nspins
            DO idir = 1, 3
               CALL dbcsr_add(matrix_ks_im(ispin)%matrix, momentum(idir)%matrix, one, -2.0_dp*vec_pot(idir))
            ENDDO
         ENDDO

         CALL dbcsr_deallocate_matrix_set(momentum)

         ! add quadratic term to real part of KS matrix
         DO ispin = 1, dft_control%nspins
            DO idir = 1, 3
               CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_s(1)%matrix, one, vec_pot(idir)**2/c_light_au)
            ENDDO
         ENDDO

      ELSE
         NULLIFY (nl_term)
         CALL dbcsr_allocate_matrix_set(nl_term, dft_control%nspins)
         DO ispin = 1, dft_control%nspins
            CALL dbcsr_init_p(nl_term(ispin)%matrix)
            CALL dbcsr_create(nl_term(ispin)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(nl_term(ispin)%matrix, sab_all)
            CALL dbcsr_set(nl_term(ispin)%matrix, zero)
         ENDDO

         CALL velocity_gauge_nl_term(qs_env, nl_term, vec_pot)

         CALL dbcsr_deallocate_matrix_set(nl_term)
         ! CPABORT(" Not implemented for GPW yet!")
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param cell ...
!> \param vec_pot ...
! **************************************************************************************************
   PURE SUBROUTINE get_vector_potential(dft_control, cell, vec_pot)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: vec_pot

      CHARACTER(len=*), PARAMETER :: routineN = 'get_vector_potential'

      ! fieldstrength
      ! for a delta pulse: vec_pot = - c * kvec
      ! the speed of light cancels in the linear term in the Hamiltonian
      vec_pot(:) = cell%h_inv(1, :)*dft_control%rtp_control%delta_pulse_direction(1) + &
                   cell%h_inv(2, :)*dft_control%rtp_control%delta_pulse_direction(2) + &
                   cell%h_inv(3, :)*dft_control%rtp_control%delta_pulse_direction(3)
      vec_pot = vec_pot*twopi*dft_control%rtp_control%delta_pulse_scale

   END SUBROUTINE get_vector_potential

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nl_term ...
!> \param vec_pot ...
! **************************************************************************************************
   SUBROUTINE velocity_gauge_nl_term(qs_env, nl_term, vec_pot)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: nl_term
      REAL(KIND=dp), DIMENSION(3), INTENT(in)            :: vec_pot

      CHARACTER(len=*), PARAMETER :: routiuneN = "velocity_gauge_nl_term"

      INTEGER                                            :: handle, i, maxl, maxlgto, maxlppnl, nkind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int_cos, sap_int_sin

      CALL timeset(routiuneN, handle)

      NULLIFY (sap_ppnl)
      CALL get_qs_env(qs_env, &
                      sap_ppnl=sap_ppnl)

      IF (ASSOCIATED(sap_ppnl)) THEN
         NULLIFY (qs_kind_set, particle_set, cell)
         CALL get_qs_env(qs_env, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, &
                         cell=cell)

         nkind = SIZE(qs_kind_set)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxlgto=maxlgto, &
                              maxlppnl=maxlppnl)

         maxl = MAX(maxlppnl, maxlgto)
         CALL init_orbital_pointers(maxl + 1)

         ! initalize sab_int types to store the integrals
         NULLIFY (sap_int_cos, sap_int_sin)
         ALLOCATE (sap_int_cos(nkind*nkind), sap_int_sin(nkind*nkind))
         DO i = 1, SIZE(sap_int_cos)
            NULLIFY (sap_int_cos(i)%alist, sap_int_cos(i)%asort, sap_int_cos(i)%aindex)
            sap_int_cos(i)%nalist = 0
            NULLIFY (sap_int_sin(i)%alist, sap_int_sin(i)%asort, sap_int_sin(i)%aindex)
            sap_int_sin(i)%nalist = 0
         ENDDO

         ! calculate exponential integrals
         CALL build_sap_exp_ints(sap_int_cos, sap_int_sin, sap_ppnl, qs_kind_set, particle_set, &
                                 cell, kvec=vec_pot)

         CALL release_sap_int(sap_int_cos)
         CALL release_sap_int(sap_int_sin)
      ENDIF
      CALL dbcsr_print(nl_term(1)%matrix)
      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_nl_term

! **************************************************************************************************
!> \brief calculate <a|sin/cos|p> integrals and store in sap_int_type
!>        adapted from build_sap_ints
!>        Do this on each MPI task as the integrals need to be available globally.
!>        Might be faster than communicating as the integrals are obtained analytically.
!> \param sap_int_cos ...
!> \param sap_int_sin ...
!> \param sap_ppnl ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param kvec ...
! **************************************************************************************************
   SUBROUTINE build_sap_exp_ints(sap_int_cos, sap_int_sin, sap_ppnl, qs_kind_set, particle_set, cell, kvec)
      TYPE(sap_int_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: sap_int_cos, sap_int_sin
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sap_ppnl
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(in)            :: kvec

      CHARACTER(len=*), PARAMETER :: routiuneN = "build_sap_exp_ints"

      INTEGER :: handle, i, iac, iatom, ikind, ilist, iset, jneighbor, katom, kkind, l, lc_max, &
         lc_min, ldai, ldints, lppnl, maxco, maxl, maxlgto, maxlppnl, maxsgf, na, nb, ncoa, ncoc, &
         nkind, nlist, nneighbor, np, nppnl, nprjc, nseta, nsgfa, prjc, sgfa, slot
      INTEGER, DIMENSION(3)                              :: cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: dogth
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ai_work_cos, ai_work_sin, &
                                                            work_cos, work_sin
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rac, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, rpgfa, sphi_a, vprj_ppnl, zeta
      TYPE(clist_type), POINTER                          :: clist, clist_sin
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routiuneN, handle)

      nkind = SIZE(qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, &
                           maxco=maxco, &
                           maxlppnl=maxlppnl, &
                           maxsgf=maxsgf, &
                           maxlgto=maxlgto)

      ! maximum dimensions for allocations
      maxl = MAX(maxlppnl, maxlgto)
      ldints = MAX(maxco, ncoset(maxlppnl), maxsgf, maxlppnl)
      ldai = ncoset(maxl + 1)

      !set up direct access to basis and potential
      NULLIFY (gpotential, spotential)
      ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, sgp_potential=sgp_potential)
         NULLIFY (gpotential(ikind)%gth_potential)
         NULLIFY (spotential(ikind)%sgp_potential)
         IF (ASSOCIATED(gth_potential)) THEN
            gpotential(ikind)%gth_potential => gth_potential
         ELSE IF (ASSOCIATED(sgp_potential)) THEN
            spotential(ikind)%sgp_potential => sgp_potential
         END IF
      END DO

      !allocate sap int
      NULLIFY (clist)
      DO slot = 1, sap_ppnl(1)%nl_size

         ikind = sap_ppnl(1)%nlist_task(slot)%ikind
         kkind = sap_ppnl(1)%nlist_task(slot)%jkind
         iatom = sap_ppnl(1)%nlist_task(slot)%iatom
         katom = sap_ppnl(1)%nlist_task(slot)%jatom
         nlist = sap_ppnl(1)%nlist_task(slot)%nlist
         ilist = sap_ppnl(1)%nlist_task(slot)%ilist
         nneighbor = sap_ppnl(1)%nlist_task(slot)%nnode

         iac = ikind + nkind*(kkind - 1)
         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(gpotential(kkind)%gth_potential) .AND. &
             .NOT. ASSOCIATED(spotential(kkind)%sgp_potential)) CYCLE
         IF (.NOT. ASSOCIATED(sap_int_cos(iac)%alist)) THEN
            sap_int_cos(iac)%a_kind = ikind
            sap_int_cos(iac)%p_kind = kkind
            sap_int_cos(iac)%nalist = nlist
            ALLOCATE (sap_int_cos(iac)%alist(nlist))
            DO i = 1, nlist
               NULLIFY (sap_int_cos(iac)%alist(i)%clist)
               sap_int_cos(iac)%alist(i)%aatom = 0
               sap_int_cos(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int_cos(iac)%alist(ilist)%clist)) THEN
            sap_int_cos(iac)%alist(ilist)%aatom = iatom
            sap_int_cos(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE (sap_int_cos(iac)%alist(ilist)%clist(nneighbor))
            DO i = 1, nneighbor
               clist => sap_int_cos(iac)%alist(ilist)%clist(i)
               clist%catom = 0
               NULLIFY (clist%acint)
               NULLIFY (clist%achint)
               NULLIFY (clist%sgf_list)
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int_sin(iac)%alist)) THEN
            sap_int_sin(iac)%a_kind = ikind
            sap_int_sin(iac)%p_kind = kkind
            sap_int_sin(iac)%nalist = nlist
            ALLOCATE (sap_int_sin(iac)%alist(nlist))
            DO i = 1, nlist
               NULLIFY (sap_int_sin(iac)%alist(i)%clist)
               sap_int_sin(iac)%alist(i)%aatom = 0
               sap_int_sin(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int_sin(iac)%alist(ilist)%clist)) THEN
            sap_int_sin(iac)%alist(ilist)%aatom = iatom
            sap_int_sin(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE (sap_int_sin(iac)%alist(ilist)%clist(nneighbor))
            DO i = 1, nneighbor
               clist => sap_int_sin(iac)%alist(ilist)%clist(i)
               clist%catom = 0
               NULLIFY (clist%acint)
               NULLIFY (clist%achint)
               NULLIFY (clist%sgf_list)
            END DO
         END IF
      END DO

      ! actual calculation of the integrals <a|cos|p> and <a|sin|p>
      ! allocate temporary storage using maximum dimensions

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (basis_set, gpotential, ncoset, sap_ppnl, sap_int_cos, sap_int_sin, nkind, &
!$OMP         ldints, ldai, nco, cell, particle_set, kvec) &
!$OMP PRIVATE (slot, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, dac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta,&
!$OMP          rpgfa, set_radius_a, sphi_a, zeta, alpha_ppnl, cprj, lppnl, nppnl, nprj_ppnl,&
!$OMP          ppnl_radius, vprj_ppnl, clist, clist_sin, ra, rc, ncoa, sgfa, prjc, work_cos, work_sin,&
!$OMP          nprjc, rprjc, lc_max, lc_min, zetc, ncoc, ai_work_sin, ai_work_cos, na, nb, np, dogth)

      ALLOCATE (work_cos(ldints, ldints), work_sin(ldints, ldints))
      ALLOCATE (ai_work_cos(ldai, ldai), ai_work_sin(ldai, ldai))
      work_cos = 0.0_dp
      work_sin = 0.0_dp
      ai_work_cos = 0.0_dp
      ai_work_sin = 0.0_dp
      dogth = .FALSE.

      NULLIFY (first_sgfa, la_max, la_min, npgfa, nsgf_seta, rpgfa, set_radius_a, sphi_a, zeta)
      NULLIFY (alpha_ppnl, cprj, nprj_ppnl, vprj_ppnl)
      NULLIFY (clist, clist_sin)

!$OMP DO SCHEDULE(GUIDED)      
      DO slot = 1, sap_ppnl(1)%nl_size
         ikind = sap_ppnl(1)%nlist_task(slot)%ikind
         kkind = sap_ppnl(1)%nlist_task(slot)%jkind
         iatom = sap_ppnl(1)%nlist_task(slot)%iatom
         katom = sap_ppnl(1)%nlist_task(slot)%jatom
         nlist = sap_ppnl(1)%nlist_task(slot)%nlist
         ilist = sap_ppnl(1)%nlist_task(slot)%ilist
         nneighbor = sap_ppnl(1)%nlist_task(slot)%nnode
         jneighbor = sap_ppnl(1)%nlist_task(slot)%inode
         cell_c(:) = sap_ppnl(1)%nlist_task(slot)%cell(:)
         rac(1:3) = sap_ppnl(1)%nlist_task(slot)%r(1:3)
         dac = NORM2(rac)

         iac = ikind + nkind*(kkind - 1)
         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         ! get definition of gto basis set
         first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
         la_max => basis_set(ikind)%gto_basis_set%lmax
         la_min => basis_set(ikind)%gto_basis_set%lmin
         npgfa => basis_set(ikind)%gto_basis_set%npgf
         nseta = basis_set(ikind)%gto_basis_set%nset
         nsgfa = basis_set(ikind)%gto_basis_set%nsgf
         nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
         rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
         set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
         sphi_a => basis_set(ikind)%gto_basis_set%sphi
         zeta => basis_set(ikind)%gto_basis_set%zet

         IF (ASSOCIATED(gpotential(kkind)%gth_potential)) THEN
            ! GTH potential
            dogth = .TRUE.
            alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
            cprj => gpotential(kkind)%gth_potential%cprj
            lppnl = gpotential(kkind)%gth_potential%lppnl
            nppnl = gpotential(kkind)%gth_potential%nppnl
            nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
            ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
            vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
         ELSE
            CYCLE
         ENDIF

         clist => sap_int_cos(iac)%alist(ilist)%clist(jneighbor)
         clist_sin => sap_int_sin(iac)%alist(ilist)%clist(jneighbor)

         clist%catom = katom
         clist%cell = cell_c
         clist%rac = rac
         clist_sin%catom = katom
         clist_sin%cell = cell_c
         clist_sin%rac = rac

         ALLOCATE (clist%acint(nsgfa, nppnl, 1), clist%achint(nsgfa, nppnl, 1))
         clist%acint = 0.0_dp
         clist%achint = 0.0_dp
         clist%nsgf_cnt = 0

         ALLOCATE (clist_sin%acint(nsgfa, nppnl, 1), clist_sin%achint(nsgfa, nppnl, 1))
         clist_sin%acint = 0.0_dp
         clist_sin%achint = 0.0_dp
         clist_sin%nsgf_cnt = 0

         ra(:) = pbc(particle_set(iatom)%r(:), cell)
         rc(:) = ra + rac

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            IF (dogth) THEN
               prjc = 1
               work_cos = 0.0_dp
               work_sin = 0.0_dp
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE
                  lc_max = l + 2*(nprj_ppnl(l) - 1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)

                  CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lc_max, 1, zetc, rprjc, lc_min, ra, rc, kvec, ai_work_cos, ai_work_sin)

                  ! projector functions: Cartesian -> spherical
                  work_cos(1:ncoa, :prjc + nprjc - 1) = &
                     MATMUL(ai_work_cos(1:ncoa, 1:ncoc), cprj(1:ncoc, prjc:prjc + nprjc - 1))
                  work_sin(1:ncoa, :prjc + nprjc - 1) = &
                     MATMUL(ai_work_sin(1:ncoa, 1:ncoc), cprj(1:ncoc, prjc:prjc + nprjc - 1))

                  prjc = prjc + nprjc
               ENDDO

               ! contract gto basis set into acint
               na = nsgf_seta(iset)
               nb = nppnl
               np = ncoa
               clist%acint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work_cos(1:np, 1:nb - 1))
               clist_sin%acint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work_sin(1:np, 1:nb - 1))

               ! multiply with interaction matrix h_ij of the nl pp
               clist%achint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(clist%acint(sgfa:sgfa + na - 1, 1:nb, 1), vprj_ppnl(1:nb, 1:nb))
               clist_sin%achint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(clist_sin%acint(sgfa:sgfa + na - 1, 1:nb, 1), vprj_ppnl(1:nb, 1:nb))

            ENDIF

         ENDDO
         clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
         clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         clist_sin%maxac = MAXVAL(ABS(clist_sin%acint(:, :, 1)))
         clist_sin%maxach = MAXVAL(ABS(clist_sin%achint(:, :, 1)))
      ENDDO

      DEALLOCATE (work_cos, work_sin, ai_work_cos, ai_work_sin)

!$OMP END PARALLEL

      DEALLOCATE (basis_set, gpotential, spotential)

      CALL timestop(handle)

   END SUBROUTINE build_sap_exp_ints

END MODULE rt_propagation_velocity_gauge
