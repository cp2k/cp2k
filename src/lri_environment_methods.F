!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates integral matrices for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) heavily debugged
!> \authors JGH
!>          Dorothea Golze
! **************************************************************************************************
MODULE lri_environment_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE generic_os_integrals,            ONLY: int_overlap_aabb_os,&
                                              int_overlap_ab_os,&
                                              int_overlap_aba_os,&
                                              int_overlap_abb_os
   USE generic_shg_integrals,           ONLY: get_abb_same_kind,&
                                              int_overlap_ab_shg_low,&
                                              int_overlap_aba_shg_low,&
                                              int_overlap_abb_shg_low,&
                                              lri_precalc_angular_shg_part
   USE input_constants,                 ONLY: do_lri_inv,&
                                              do_lri_inv_auto,&
                                              do_lri_pseudoinv_diag,&
                                              do_lri_pseudoinv_svd
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: &
        allocate_lri_coefs, allocate_lri_ints, allocate_lri_ints_rho, allocate_lri_rhos, &
        deallocate_lri_ints, deallocate_lri_ints_rho, lri_density_create, lri_density_release, &
        lri_density_type, lri_environment_type, lri_int_rho_type, lri_int_type, lri_kind_type, &
        lri_list_type, lri_rhoab_type
   USE mathlib,                         ONLY: get_pseudo_inverse_diag,&
                                              get_pseudo_inverse_svd,&
                                              invmat_symm
   USE message_passing,                 ONLY: mp_max,&
                                              mp_sum
   USE particle_types,                  ONLY: particle_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_collocate_density,            ONLY: calculate_lri_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_methods'

   PUBLIC :: build_lri_matrices, calculate_lri_densities, &
             calculate_lri_integrals, calculate_lri_overlap_aabb, calculate_avec_lri, &
             lri_int_os, lri_int_shg

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_lri_matrices(lri_env, qs_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_lri_matrices', &
         routineP = moduleN//':'//routineN

      ! calculate the integrals needed to do the local (2-center) expansion
      ! of the (pair) densities
      CALL calculate_lri_integrals(lri_env, qs_env)

   END SUBROUTINE build_lri_matrices

! **************************************************************************************************
!> \brief calculates integrals needed for the LRI density fitting,
!>        integrals are calculated once, before the SCF starts
!> \param lri_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_lri_integrals(lri_env, qs_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iac, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, mepos, nba, nbb, &
                                                            nfa, nfb, nkind, nlist, nn, nneighbor, &
                                                            nthread
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: fbasa, fbasb, obasa, obasb
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(lri_list_type), POINTER                       :: lri_ints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, fbasa, fbasb, lrii, lri_ints, nl_iterator, &
               obasa, obasb, particle_set, soo_list, virial)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         CALL get_qs_env(qs_env=qs_env, cell=cell, dft_control=dft_control, &
                         nkind=nkind, particle_set=particle_set, virial=virial)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
         nthread = 1
!$       nthread = omp_get_max_threads()

         IF (ASSOCIATED(lri_env%lri_ints)) THEN
            CALL deallocate_lri_ints(lri_env%lri_ints)
         END IF

         ! allocate matrices storing the LRI integrals
         CALL allocate_lri_ints(lri_env, lri_env%lri_ints, nkind)
         lri_ints => lri_env%lri_ints

         CALL neighbor_list_iterator_create(nl_iterator, soo_list, nthread=nthread)
!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,nl_iterator,lri_env,lri_ints,nkind,use_virial)&
!$OMP PRIVATE (mepos,ikind,jkind,iatom,jatom,nlist,ilist,nneighbor,jneighbor,rab,dab,&
!$OMP          iac,obasa,obasb,fbasa,fbasb,lrii,nba,nbb,nfa,nfb,nn)

         mepos = 0
!$       mepos = omp_get_thread_num()

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)

            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                   nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                   iatom=iatom, jatom=jatom, r=rab)
            iac = ikind+nkind*(jkind-1)
            dab = SQRT(SUM(rab*rab))

            obasa => lri_env%orb_basis(ikind)%gto_basis_set
            obasb => lri_env%orb_basis(jkind)%gto_basis_set
            fbasa => lri_env%ri_basis(ikind)%gto_basis_set
            fbasb => lri_env%ri_basis(jkind)%gto_basis_set

            IF (.NOT. ASSOCIATED(obasa)) CYCLE
            IF (.NOT. ASSOCIATED(obasb)) CYCLE

            lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

            ! forces: not every pair is giving contributions
            ! no forces for self-pair aa
            IF (iatom == jatom .AND. dab < lri_env%delta) THEN
               lrii%calc_force_pair = .FALSE.
            ELSE
               ! forces for periodic self-pair aa' required for virial
               IF (iatom == jatom .AND. .NOT. use_virial) THEN
                  lrii%calc_force_pair = .FALSE.
               ELSE
                  lrii%calc_force_pair = .TRUE.
               ENDIF
            ENDIF

            nba = obasa%nsgf
            nbb = obasb%nsgf
            nfa = fbasa%nsgf
            nfb = fbasb%nsgf

            lrii%nba = nba
            lrii%nbb = nbb
            lrii%nfa = nfa
            lrii%nfb = nfb

            ! full method is used
            ! calculate integrals (fa,fb), (a,b), (a,b,fa) and (a,b,fb)
            IF (lri_env%use_shg_integrals) THEN
               CALL lri_int_shg(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                                iatom, jatom, ikind, jkind, .FALSE.)
            ELSE
               CALL lri_int_os(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                               iatom, jatom, ikind, .FALSE.)
            ENDIF

            ! Full LRI method
            IF (lrii%lrisr) THEN
               ! construct and invert S matrix
               IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                  lrii%sinv(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp_inv(1:nfa, 1:nfa)
               ELSE
                  nn = nfa+nfb
                  lrii%sinv(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa, 1:nfa)
                  lrii%sinv(1:nfa, nfa+1:nn) = lrii%sab(1:nfa, 1:nfb)
                  lrii%sinv(nfa+1:nn, 1:nfa) = TRANSPOSE(lrii%sab(1:nfa, 1:nfb))
                  lrii%sinv(nfa+1:nn, nfa+1:nn) = lri_env%bas_prop(jkind)%ri_ovlp(1:nfb, 1:nfb)
                  CALL inverse_lri_overlap(lri_env, lrii%sinv)
               ENDIF

               ! calculate Sinv*n and n*Sinv*n
               lrii%n(1:nfa) = lri_env%bas_prop(ikind)%int_fbas(1:nfa)
               IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                  CALL dgemv("N", nfa, nfa, 1.0_dp, lrii%sinv(1, 1), nfa, &
                             lrii%n(1), 1, 0.0_dp, lrii%sn, 1)
                  lrii%nsn = SUM(lrii%sn(1:nfa)*lrii%n(1:nfa))
               ELSE
                  lrii%n(nfa+1:nn) = lri_env%bas_prop(jkind)%int_fbas(1:nfb)
                  CALL dgemv("N", nn, nn, 1.0_dp, lrii%sinv(1, 1), nn, &
                             lrii%n(1), 1, 0.0_dp, lrii%sn, 1)
                  lrii%nsn = SUM(lrii%sn(1:nn)*lrii%n(1:nn))
               ENDIF
            END IF

            ! Distant Pair methods
            IF (lrii%lriff) THEN
               ! calculate Sinv*n and n*Sinv*n for A and B centers
               lrii%na(1:nfa) = lri_env%bas_prop(ikind)%int_fbas(1:nfa)
               lrii%nb(1:nfb) = lri_env%bas_prop(jkind)%int_fbas(1:nfb)
               CALL dgemv("N", nfa, nfa, 1.0_dp, lrii%asinv(1, 1), nfa, &
                          lrii%na(1), 1, 0.0_dp, lrii%sna, 1)
               lrii%nsna = SUM(lrii%sna(1:nfa)*lrii%na(1:nfa))
               CALL dgemv("N", nfb, nfb, 1.0_dp, lrii%bsinv(1, 1), nfb, &
                          lrii%nb(1), 1, 0.0_dp, lrii%snb, 1)
               lrii%nsnb = SUM(lrii%snb(1:nfb)*lrii%nb(1:nfb))
            END IF

         END DO
!$OMP END PARALLEL

         CALL neighbor_list_iterator_release(nl_iterator)

         IF (lri_env%debug) THEN
            CALL output_debug_info(lri_env, qs_env, lri_ints, soo_list)
         ENDIF

      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_lri_integrals

! **************************************************************************************************
!> \brief calculates overlap integrals (aabb) of the orbital basis set,
!>        reguired for LRI basis set optimization
!> \param lri_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_lri_overlap_aabb(lri_env, qs_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_overlap_aabb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iac, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, nba, nbb, nkind, &
                                                            nlist, nneighbor
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb
      TYPE(lri_int_rho_type), POINTER                    :: lriir
      TYPE(lri_list_type), POINTER                       :: lri_ints_rho
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      NULLIFY (cell, lriir, lri_ints_rho, nl_iterator, obasa, obasb, &
               particle_set, soo_list)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         CALL get_qs_env(qs_env=qs_env, nkind=nkind, particle_set=particle_set, &
                         cell=cell)

         IF (ASSOCIATED(lri_env%lri_ints_rho)) THEN
            CALL deallocate_lri_ints_rho(lri_env%lri_ints_rho)
         END IF

         CALL allocate_lri_ints_rho(lri_env, lri_env%lri_ints_rho, nkind)
         lri_ints_rho => lri_env%lri_ints_rho

         CALL neighbor_list_iterator_create(nl_iterator, soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                   iatom=iatom, jatom=jatom, r=rab)

            iac = ikind+nkind*(jkind-1)
            dab = SQRT(SUM(rab*rab))

            obasa => lri_env%orb_basis(ikind)%gto_basis_set
            obasb => lri_env%orb_basis(jkind)%gto_basis_set
            IF (.NOT. ASSOCIATED(obasa)) CYCLE
            IF (.NOT. ASSOCIATED(obasb)) CYCLE

            lriir => lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)

            nba = obasa%nsgf
            nbb = obasb%nsgf

            ! calculate integrals (aa,bb)
            CALL int_overlap_aabb_os(lriir%soaabb, obasa, obasb, rab, lri_env%debug, &
                                     lriir%dmax_aabb)

         END DO

         CALL neighbor_list_iterator_release(nl_iterator)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE calculate_lri_overlap_aabb

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
   SUBROUTINE calculate_lri_densities(lri_env, qs_env, pmatrix, &
                                      lri_rho_struct, atomic_kind_set, para_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix
      TYPE(qs_rho_type), POINTER                         :: lri_rho_struct
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_densities', &
         routineP = moduleN//':'//routineN

      TYPE(lri_density_type), POINTER                    :: lri_density

      CALL get_qs_env(qs_env, lri_density=lri_density)
      CALL calculate_avec_lri(lri_env, lri_density, pmatrix)
      CALL set_qs_env(qs_env, lri_density=lri_density)

      CALL distribute_lri_density_on_the_grid(lri_env, lri_density, qs_env, &
                                              lri_rho_struct, atomic_kind_set, para_env)

   END SUBROUTINE calculate_lri_densities

! **************************************************************************************************
!> \brief performs the fitting of the density; solves the linear system of
!>        equations; yield the expansion coefficients avec
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!> \param lri_density ...
!> \param pmatrix ...
! **************************************************************************************************
   SUBROUTINE calculate_avec_lri(lri_env, lri_density, pmatrix)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_avec_lri', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iac, iatom, ikind, ilist, ispin, jatom, jkind, jneighbor, mepos, nba, &
         nbb, nfa, nfb, nkind, nlist, nn, nneighbor, nspin, nthread
      LOGICAL                                            :: found, trans
      REAL(KIND=dp)                                      :: dab, rab(3), threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: m
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: paa, pab, pbb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pbij
      TYPE(dbcsr_type), POINTER                          :: pmat
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(lri_list_type), POINTER                       :: lri_rho
      TYPE(lri_rhoab_type), POINTER                      :: lrho
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      CALL timeset(routineN, handle)
      NULLIFY (lrii, lri_rho, nl_iterator, pbij, pmat, soo_list)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nspin = SIZE(pmatrix)
         nkind = lri_env%lri_ints%nkind
         nthread = 1
!$       nthread = omp_get_max_threads()

         CALL lri_density_release(lri_density)
         CALL lri_density_create(lri_density)
         lri_density%nspin = nspin

         ! allocate structure lri_rhos and vectors tvec and avec
         CALL allocate_lri_rhos(lri_env, lri_density%lri_rhos, nspin, nkind)

         DO ispin = 1, nspin
            pmat => pmatrix(ispin)%matrix
            lri_rho => lri_density%lri_rhos(ispin)%lri_list

            CALL neighbor_list_iterator_create(nl_iterator, soo_list, nthread=nthread)
!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,nl_iterator,lri_env,lri_rho,pmat,nkind)&
!$OMP PRIVATE (mepos,ikind,jkind,iatom,jatom,nlist,ilist,nneighbor,jneighbor,rab,dab,iac,&
!$OMP          trans,found,pbij,pab,paa,pbb,lrho,lrii,nba,nbb,nfa,nfb,nn,threshold,i,m)

            mepos = 0
!$          mepos = omp_get_thread_num()
            DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
               CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                      jatom=jatom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                      r=rab)

               iac = ikind+nkind*(jkind-1)
               dab = SQRT(SUM(rab*rab))

               IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

               ! get the density matrix Pab
               NULLIFY (pbij)
               IF (iatom <= jatom) THEN
                  CALL dbcsr_get_block_p(matrix=pmat, row=iatom, col=jatom, block=pbij, found=found)
                  trans = .FALSE.
               ELSE
                  CALL dbcsr_get_block_p(matrix=pmat, row=jatom, col=iatom, block=pbij, found=found)
                  trans = .TRUE.
               END IF
               CPASSERT(found)

               lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
               lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

               nba = lrii%nba
               nbb = lrii%nbb
               nfa = lrii%nfa
               nfb = lrii%nfb

               nn = nfa+nfb

               ! compute tvec = SUM_ab Pab *(a,b,x) and charge constraint
               ALLOCATE (pab(nba, nbb))
               IF (trans) THEN
                  pab(1:nba, 1:nbb) = TRANSPOSE(pbij(1:nbb, 1:nba))
               ELSE
                  pab(1:nba, 1:nbb) = pbij(1:nba, 1:nbb)
               END IF

               IF (lrii%lrisr) THEN
                  ! full LRI method
                  threshold = lri_env%eps_o3_int/MAX(SUM(ABS(pab(1:nba, 1:nbb))), 1.0e-14_dp)
                  lrho%charge = SUM(pab(1:nba, 1:nbb)*lrii%soo(1:nba, 1:nbb))
                  DO i = 1, nfa
                     IF (lrii%abascr(i) > threshold) THEN
                        lrho%tvec(i) = SUM(pab(1:nba, 1:nbb)*lrii%abaint(1:nba, 1:nbb, i))
                     ELSE
                        lrho%tvec(i) = 0.0_dp
                     END IF
                  END DO
                  IF (dab > lri_env%delta) THEN
                     DO i = 1, nfb
                        IF (lrii%abbscr(i) > threshold) THEN
                           lrho%tvec(nfa+i) = SUM(pab(1:nba, 1:nbb)*lrii%abbint(1:nba, 1:nbb, i))
                        ELSE
                           lrho%tvec(nfa+i) = 0.0_dp
                        END IF
                     END DO
                  ENDIF
                  !
                  IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                     lrho%nst = SUM(lrho%tvec(1:nfa)*lrii%sn(1:nfa))
                  ELSE
                     lrho%nst = SUM(lrho%tvec(1:nn)*lrii%sn(1:nn))
                  ENDIF
                  lrho%lambda = (lrho%charge-lrho%nst)/lrii%nsn
                  !
                  ! solve the linear system of equations
                  ALLOCATE (m(nn))
                  m = 0._dp
                  IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                     m(1:nfa) = lrho%tvec(1:nfa)+lrho%lambda*lrii%n(1:nfa)
                     CALL dgemv("N", nfa, nfa, 1.0_dp, lrii%sinv(1, 1), nfa, &
                                m(1), 1, 0.0_dp, lrho%avec, 1)
                  ELSE
                     m(1:nn) = lrho%tvec(1:nn)+lrho%lambda*lrii%n(1:nn)
                     CALL dgemv("N", nn, nn, 1.0_dp, lrii%sinv(1, 1), nn, &
                                m(1), 1, 0.0_dp, lrho%avec, 1)
                  ENDIF
                  DEALLOCATE (m)
               ENDIF

               IF (lrii%lriff) THEN
                  ! distant pair approximations
                  ALLOCATE (paa(nba, nbb), pbb(nba, nbb))
                  paa(1:nba, 1:nbb) = pab(1:nba, 1:nbb)
                  pbb(1:nba, 1:nbb) = pab(1:nba, 1:nbb)
                  !
                  SELECT CASE (lri_env%distant_pair_method)
                  CASE ("EW")
                     paa(1:nba, 1:nbb) = 0.5_dp*paa(1:nba, 1:nbb)
                     pbb(1:nba, 1:nbb) = 0.5_dp*pbb(1:nba, 1:nbb)
                  CASE DEFAULT
                     CPABORT("Unknown DISTANT_PAIR_METHOD in LRI")
                  END SELECT
                  !
                  threshold = lri_env%eps_o3_int/MAX(SUM(ABS(paa(1:nba, 1:nbb))), 1.0e-14_dp)
                  lrho%chargea = SUM(paa(1:nba, 1:nbb)*lrii%soo(1:nba, 1:nbb))
                  DO i = 1, nfa
                     IF (lrii%abascr(i) > threshold) THEN
                        lrho%tveca(i) = SUM(paa(1:nba, 1:nbb)*lrii%abaint(1:nba, 1:nbb, i))
                     ELSE
                        lrho%tveca(i) = 0.0_dp
                     END IF
                  END DO
                  threshold = lri_env%eps_o3_int/MAX(SUM(ABS(pbb(1:nba, 1:nbb))), 1.0e-14_dp)
                  lrho%chargeb = SUM(pbb(1:nba, 1:nbb)*lrii%soo(1:nba, 1:nbb))
                  DO i = 1, nfb
                     IF (lrii%abbscr(i) > threshold) THEN
                        lrho%tvecb(i) = SUM(pbb(1:nba, 1:nbb)*lrii%abbint(1:nba, 1:nbb, i))
                     ELSE
                        lrho%tvecb(i) = 0.0_dp
                     END IF
                  END DO
                  !
                  lrho%nsta = SUM(lrho%tveca(1:nfa)*lrii%sna(1:nfa))
                  lrho%nstb = SUM(lrho%tvecb(1:nfb)*lrii%snb(1:nfb))
                  lrho%lambdaa = (lrho%chargea-lrho%nsta)/lrii%nsna
                  lrho%lambdab = (lrho%chargeb-lrho%nstb)/lrii%nsnb
                  ! solve the linear system of equations
                  ALLOCATE (m(nfa))
                  m(1:nfa) = lrho%tveca(1:nfa)+lrho%lambdaa*lrii%na(1:nfa)
                  CALL dgemv("N", nfa, nfa, 1.0_dp, lrii%asinv(1, 1), nfa, &
                             m(1), 1, 0.0_dp, lrho%aveca, 1)
                  DEALLOCATE (m)
                  ALLOCATE (m(nfb))
                  m(1:nfb) = lrho%tvecb(1:nfb)+lrho%lambdab*lrii%nb(1:nfb)
                  CALL dgemv("N", nfb, nfb, 1.0_dp, lrii%bsinv(1, 1), nfb, &
                             m(1), 1, 0.0_dp, lrho%avecb, 1)
                  DEALLOCATE (m)
                  !
                  DEALLOCATE (paa, pbb)
               ENDIF

               DEALLOCATE (pab)

            END DO
!$OMP END PARALLEL
            CALL neighbor_list_iterator_release(nl_iterator)

         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_avec_lri

! **************************************************************************************************
!> \brief sums up avec and  distributes the fitted density on the grid
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!>        lri_rho_struct where the fitted density is stored
!> \param lri_density ...
!> \param qs_env ...
!> \param lri_rho_struct ...
!> \param atomic_kind_set ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE distribute_lri_density_on_the_grid(lri_env, lri_density, qs_env, &
                                                 lri_rho_struct, atomic_kind_set, para_env)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: lri_rho_struct
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'distribute_lri_density_on_the_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, atom_b, handle, iac, iatom, &
                                                            ikind, ilist, ispin, jatom, jkind, &
                                                            jneighbor, nat, natom, nfa, nfb, &
                                                            nkind, nspin
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      REAL(KIND=dp)                                      :: dab, fw, rab(3)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aci, acj, tot_rho_r
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(lri_list_type), POINTER                       :: lri_rho
      TYPE(lri_rhoab_type), POINTER                      :: lrho
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r

      CALL timeset(routineN, handle)
      NULLIFY (aci, acj, atomic_kind, atom_of_kind, lri_coef, lri_rho, &
               nl_iterator, soo_list, rho_r, rho_g, tot_rho_r)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nspin = lri_density%nspin
         nkind = lri_env%lri_ints%nkind

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=nat)
         ALLOCATE (atom_of_kind(nat))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)

         ! allocate the arrays to hold RI expansion coefficients lri_coefs
         CALL allocate_lri_coefs(lri_env, lri_density, atomic_kind_set)
         DO ispin = 1, nspin

            lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
            lri_rho => lri_density%lri_rhos(ispin)%lri_list

            ! sum up expansion coefficients
            CALL neighbor_list_iterator_create(nl_iterator, soo_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, ilist=ilist, inode=jneighbor, r=rab)
               dab = SQRT(SUM(rab*rab))
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               aci => lri_coef(ikind)%acoef(atom_a, :)
               acj => lri_coef(jkind)%acoef(atom_b, :)
               iac = ikind+nkind*(jkind-1)
               lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
               lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
               nfa = lrho%nfa
               nfb = lrho%nfb

               IF (lrii%lrisr) THEN
                  IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                     !self pair aa
                     aci(1:nfa) = aci(1:nfa)+lrho%avec(1:nfa)
                  ELSE
                     IF (iatom == jatom) THEN
                        !periodic self pair aa'
                        fw = lrii%wsr
                     ELSE
                        !pairs ab
                        fw = 2.0_dp*lrii%wsr
                     ENDIF
                     aci(1:nfa) = aci(1:nfa)+fw*lrho%avec(1:nfa)
                     acj(1:nfb) = acj(1:nfb)+fw*lrho%avec(nfa+1:nfa+nfb)
                  ENDIF
               ENDIF
               !
               IF (lrii%lriff) THEN
                  IF (iatom == jatom) THEN
                     fw = lrii%wff
                  ELSE
                     fw = 2.0_dp*lrii%wff
                  ENDIF
                  aci(1:nfa) = aci(1:nfa)+fw*lrho%aveca(1:nfa)
                  acj(1:nfb) = acj(1:nfb)+fw*lrho%avecb(1:nfb)
               ENDIF

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)

            ! replicate the acoef infomation
            DO ikind = 1, nkind
               atomic_kind => atomic_kind_set(ikind)
               CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
               DO iatom = 1, natom
                  aci => lri_coef(ikind)%acoef(iatom, :)
                  CALL mp_sum(aci, para_env%group)
               END DO
            END DO

         END DO

         !distribute fitted density on the grid
         CALL qs_rho_get(lri_rho_struct, rho_r=rho_r, rho_g=rho_g, &
                         tot_rho_r=tot_rho_r)
         DO ispin = 1, nspin
            CALL calculate_lri_rho_elec(rho_g(ispin), &
                                        rho_r(ispin), qs_env, &
                                        lri_density%lri_coefs(ispin)%lri_kinds, &
                                        tot_rho_r(ispin), "LRI")
         ENDDO

         CALL set_qs_env(qs_env, lri_density=lri_density)

         DEALLOCATE (atom_of_kind)

      END IF

      CALL timestop(handle)

   END SUBROUTINE distribute_lri_density_on_the_grid

! **************************************************************************************************
!> \brief calculates the lri intergrals according to the Obara-Saika (OS)
!>        scheme
!> \param lri_env ...
!> \param lrii ...
!> \param rab distance vector
!> \param obasa orb basis on center A
!> \param obasb orb basis on center B
!> \param fbasa aux basis on center A
!> \param fbasb aux basis on center B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_int_os(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                         iatom, jatom, ikind, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_os', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nba, nbb, nfa, nfb
      LOGICAL                                            :: do_force
      REAL(KIND=dp)                                      :: dab

      CALL timeset(routineN, handle)

      dab = SQRT(SUM(rab*rab))
      nba = obasa%nsgf
      nbb = obasb%nsgf
      nfa = fbasa%nsgf
      nfb = fbasb%nsgf

      !*** calculate overlap integrals; for iatom=jatom this is the self-overlap
      IF (iatom == jatom .AND. dab < lri_env%delta) THEN
         !*** (fa,fa)
         lrii%sab(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa, 1:nfa)
         !*** (a,a)
         lrii%soo(1:nba, 1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba, 1:nba)
         !*** (a,a,fa)
         CALL int_overlap_aba_os(lrii%abaint, rab=rab, oba=obasa, obb=obasb, &
                                 fba=fbasa, calculate_forces=.FALSE., debug=lri_env%debug, &
                                 dmax=lrii%dmax_aba, screen=lrii%abascr)
         IF (calculate_forces) THEN
            lrii%dsab = 0._dp
            lrii%dsoo = 0._dp
            lrii%dabdaint = 0.0_dp
            lrii%dabbint = 0.0_dp
         END IF
      ELSE
         IF (calculate_forces) THEN
            do_force = lrii%calc_force_pair
         ELSE
            do_force = .FALSE.
         END IF
         !*** (fa,fb)
         CALL int_overlap_ab_os(lrii%sab, lrii%dsab, rab, fbasa, fbasb, &
                                do_force, lri_env%debug, lrii%dmax_ab)
         !*** (a,b)
         CALL int_overlap_ab_os(lrii%soo, lrii%dsoo, rab, obasa, obasb, &
                                do_force, lri_env%debug, lrii%dmax_oo)
         !*** (a,b,fa)
         CALL int_overlap_aba_os(lrii%abaint, lrii%dabdaint, rab, obasa, obasb, fbasa, &
                                 do_force, lri_env%debug, lrii%dmax_aba, screen=lrii%abascr)
         !*** (a,b,fb)
         CALL int_overlap_abb_os(lrii%abbint, lrii%dabbint, rab, obasa, obasb, fbasb, &
                                 do_force, lri_env%debug, lrii%dmax_abb, screen=lrii%abbscr)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE lri_int_os

! **************************************************************************************************
!> \brief calcuates the lri integrals using solid harmonic Gaussians
!> \param lri_env ...
!> \param lrii ...
!> \param rab distance vector
!> \param obasa orb basis on A
!> \param obasb orb basis on B
!> \param fbasa aux basis on A
!> \param fbasb aux basis on B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param jkind kind atom B
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_int_shg(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                          iatom, jatom, ikind, jkind, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_shg', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nba, nbb, nfa, nfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: fba_index, fbb_index, oba_index, &
                                                            obb_index
      LOGICAL                                            :: do_force
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_fba, scon_fbb, scon_oba, scon_obb
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: scona_mix, sconb_mix

      CALL timeset(routineN, handle)
      NULLIFY (scon_oba, scon_obb, scon_fba, scon_fbb, scona_mix, sconb_mix, &
               oba_index, obb_index, fba_index, fbb_index)
      dab = SQRT(SUM(rab*rab))
      nba = obasa%nsgf
      nbb = obasb%nsgf
      nfa = fbasa%nsgf
      nfb = fbasb%nsgf

      !*** calculate overlap integrals; for iatom=jatom this is the self-overlap
      IF (iatom == jatom .AND. dab < lri_env%delta) THEN
         !*** (fa,fa)
         lrii%sab(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa, 1:nfa)
         !*** (a,a)
         lrii%soo(1:nba, 1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba, 1:nba)
         !*** (a,a,fa)
         lrii%abaint(1:nba, 1:nba, 1:nfa) = lri_env%bas_prop(ikind)%ovlp3
         lrii%abascr(1:nfa) = 1.0_dp
         IF (calculate_forces) THEN
            lrii%dsab = 0._dp
            lrii%dsoo = 0._dp
            lrii%dabdaint = 0.0_dp
            lrii%dabbint = 0.0_dp
         END IF
      ELSE
         IF (calculate_forces) THEN
            do_force = lrii%calc_force_pair
         ELSE
            do_force = .FALSE.
         END IF
         scon_oba => lri_env%bas_prop(ikind)%scon_orb
         scon_obb => lri_env%bas_prop(jkind)%scon_orb
         scon_fba => lri_env%bas_prop(ikind)%scon_ri
         scon_fbb => lri_env%bas_prop(jkind)%scon_ri
         scona_mix => lri_env%bas_prop(ikind)%scon_mix
         sconb_mix => lri_env%bas_prop(jkind)%scon_mix
         oba_index => lri_env%bas_prop(ikind)%orb_index
         fba_index => lri_env%bas_prop(ikind)%ri_index
         obb_index => lri_env%bas_prop(jkind)%orb_index
         fbb_index => lri_env%bas_prop(jkind)%ri_index
         CALL lri_precalc_angular_shg_part(obasa, obasb, fbasa, fbasb, rab, Waux_mat, dWaux_mat, &
                                           do_force)
         !*** (fa,fb)
         CALL int_overlap_ab_shg_low(lrii%sab, lrii%dsab, rab, fbasa, fbasb, scon_fba, scon_fbb, &
                                     Waux_mat, dWaux_mat, do_force, contraction_high=.FALSE.)
         !*** (a,b)
         CALL int_overlap_ab_shg_low(lrii%soo, lrii%dsoo, rab, obasa, obasb, scon_oba, scon_obb, &
                                     Waux_mat, dWaux_mat, do_force, contraction_high=.TRUE.)
         !*** (a,b,fa)
         CALL int_overlap_aba_shg_low(lrii%abaint, lrii%dabdaint, rab, obasa, obasb, fbasa, &
                                      scon_obb, scona_mix, oba_index, fba_index, &
                                      lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                      lri_env%cg_shg%ncg_none0, &
                                      Waux_mat, dWaux_mat, do_force, lrii%abascr)
         !*** (a,b,fb)
         IF (ikind == jkind) THEN
            CALL get_abb_same_kind(lrii%abbint, lrii%dabbint, lrii%abaint, lrii%dabdaint, &
                                   rab, obasa, fbasa, do_force)
            lrii%abbscr(:) = lrii%abascr(:)
         ELSE
            CALL int_overlap_abb_shg_low(lrii%abbint, lrii%dabbint, rab, obasa, obasb, fbasb, &
                                         scon_oba, sconb_mix, obb_index, fbb_index, &
                                         lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                         lri_env%cg_shg%ncg_none0, &
                                         Waux_mat, dWaux_mat, do_force, lrii%abbscr)
         ENDIF

         DEALLOCATE (Waux_mat, dWaux_mat)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE lri_int_shg

! **************************************************************************************************
!> \brief get inverse or pseudoinverse of lri overlap matrix for aux basis set
!> \param lri_env lri environment
!> \param sinv on entry overlap matrix, on exit its inverse
! **************************************************************************************************
   SUBROUTINE inverse_lri_overlap(lri_env, sinv)

      TYPE(lri_environment_type)                         :: lri_env
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sinv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'inverse_lri_overlap', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, info, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: anorm, rcond, rskip
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: s
      REAL(KIND=dp), EXTERNAL                            :: dlange

      CALL timeset(routineN, handle)

      NULLIFY (s)

      rskip = 1.E-8_dp ! parameter for pseudo inverse
      n = SIZE(sinv, 1)

      SELECT CASE (lri_env%lri_overlap_inv)
      CASE (do_lri_inv)
         CALL invmat_symm(sinv)
      CASE (do_lri_pseudoinv_svd)
         ALLOCATE (s(n, n))
         s(:, :) = sinv
         CALL get_pseudo_inverse_svd(s, sinv, rskip)
         DEALLOCATE (s)
      CASE (do_lri_pseudoinv_diag)
         ALLOCATE (s(n, n))
         s(:, :) = sinv
         CALL get_pseudo_inverse_diag(s, sinv, rskip)
         DEALLOCATE (s)
      CASE (do_lri_inv_auto)
         ! decide whether to calculate inverse or pseudoinverse
         ALLOCATE (s(n, n))
         s(:, :) = sinv
         ALLOCATE (work(3*n), iwork(n))
         ! norm of matrix
         anorm = dlange('1', n, n, sinv, n, work)
         ! Cholesky factorization (fails if matrix not positive definite)
         CALL dpotrf('U', n, sinv, n, info)
         IF (info == 0) THEN
            ! condition number
            CALL dpocon('U', n, sinv, n, anorm, rcond, work, iwork, info)
            IF (info /= 0) THEN
               CPABORT("DPOCON failed")
            ENDIF
            IF (LOG(1._dp/rcond) > lri_env%cond_max) THEN
               CALL get_pseudo_inverse_svd(s, sinv, rskip)
            ELSE
               CALL invmat_symm(sinv, "U")
            ENDIF
         ELSE
            CALL get_pseudo_inverse_svd(s, sinv, rskip)
         END IF
         DEALLOCATE (s, work, iwork)
      CASE DEFAULT
         CPABORT("No initialization for LRI overlap available?")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE inverse_lri_overlap

! **************************************************************************************************
!> \brief debug output
!> \param lri_env ...
!> \param qs_env ...
!> \param lri_ints ...
!> \param soo_list ...
! **************************************************************************************************
   SUBROUTINE output_debug_info(lri_env, qs_env, lri_ints, soo_list)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(lri_list_type), POINTER                       :: lri_ints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      CHARACTER(LEN=*), PARAMETER :: routineN = 'output_debug_info', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iac, ikind, ilist, iunit, jkind, &
                                                            jneighbor, nkind
      REAL(KIND=dp)                                      :: dmax_aabb, dmax_ab, dmax_aba, dmax_abb, &
                                                            dmax_oo
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(lri_int_rho_type), POINTER                    :: lriir
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      NULLIFY (input, logger, lrii, lriir, nl_iterator, para_env)
      CALL get_qs_env(qs_env, dft_control=dft_control, input=input, nkind=nkind, &
                      para_env=para_env)
      dmax_ab = 0._dp
      dmax_oo = 0._dp
      dmax_aba = 0._dp
      dmax_abb = 0._dp
      dmax_aabb = 0._dp

      CALL neighbor_list_iterator_create(nl_iterator, soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                ilist=ilist, inode=jneighbor)

         iac = ikind+nkind*(jkind-1)
         lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

         dmax_ab = MAX(dmax_ab, lrii%dmax_ab)
         dmax_oo = MAX(dmax_oo, lrii%dmax_oo)
         dmax_aba = MAX(dmax_aba, lrii%dmax_aba)
         dmax_abb = MAX(dmax_abb, lrii%dmax_abb)

         IF (dft_control%qs_control%lri_optbas) THEN
            lriir => lri_env%lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)
            dmax_aabb = MAX(dmax_aabb, lriir%dmax_aabb)
         ENDIF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL mp_max(dmax_ab, para_env%group)
      CALL mp_max(dmax_oo, para_env%group)
      CALL mp_max(dmax_aba, para_env%group)
      CALL mp_max(dmax_abb, para_env%group)
      CALL mp_max(dmax_aabb, para_env%group)

      logger => cp_get_default_logger()
      iunit = cp_print_key_unit_nr(logger, input, "PRINT%PROGRAM_RUN_INFO", &
                                   extension=".lridebug")

      IF (iunit > 0) THEN
         WRITE (iunit, FMT="(/,T2,A)") "DEBUG INFO FOR LRI INTEGRALS"
         WRITE (iunit, FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "// &
            "[ai|bi]; fit basis", dmax_ab
         WRITE (iunit, FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "// &
            "[a|b]; orbital basis", dmax_oo
         WRITE (iunit, FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "// &
            "[a|b|ai]", dmax_aba
         WRITE (iunit, FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "// &
            "[a|b|bi]", dmax_abb
         IF (dft_control%qs_control%lri_optbas) THEN
            WRITE (iunit, FMT="(T2,A,T69,ES12.5,/)") "Maximal deviation of integrals "// &
               "[aa|bb]; orbital basis", &
               dmax_aabb
         ENDIF
      ENDIF

      CALL cp_print_key_finished_output(iunit, logger, input, &
                                        "PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE output_debug_info

END MODULE lri_environment_methods
