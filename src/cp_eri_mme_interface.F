!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Minimax-Ewald method for periodic ERI's to be used in CP2K.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE cp_eri_mme_interface
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_create,&
                                              cell_release,&
                                              cell_type,&
                                              init_cell,&
                                              pbc
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                              cp_dbcsr_p_type
   USE cp_eri_mme_types,                ONLY: cp_eri_mme_finalize,&
                                              cp_eri_mme_init_read_input,&
                                              cp_eri_mme_param,&
                                              cp_eri_mme_set_params
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate,&
                                              eri_mme_2c_integrate_set,&
                                              eri_mme_3c_integrate_set
   USE eri_mme_test,                    ONLY: eri_mme_2c_perf_acc_test,&
                                              eri_mme_3c_perf_acc_test
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_eri_mme_interface'

   PUBLIC :: cp_eri_mme_2c_integrate, &
             cp_eri_mme_2c_integrate_s, &
             cp_eri_mme_3c_integrate, &
             cp_eri_mme_perf_acc_test

CONTAINS

! **************************************************************************************************
!> \brief high-level integration routine for 2c integrals over CP2K basis sets.
!>        Contiguous column-wise distribution and parallelization over pairs of sets.
!> \param param ...
!> \param para_env mpi environment for local columns
!> \param qs_env ...
!> \param basis_type_a ...
!> \param basis_type_b ...
!> \param hab columns of ERI matrix
!> \param first_b first column of hab
!> \param last_b last column of hab
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_2c_integrate(param, para_env, qs_env, basis_type_a, basis_type_b, hab, first_b, last_b)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_a, basis_type_b
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      INTEGER, INTENT(IN)                                :: first_b, last_b

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_eri_mme_2c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_end, atom_start, first_set, G_count, handle, iatom, ikind, iset, jatom, &
         jkind, jset, jset_end, jset_start, last_set, my_setpair, n_setpair, natom, nkind, nseta, &
         nseta_total, nsetb, nsetb_total, offset_a_end, offset_a_start, offset_b_end, &
         offset_b_start, R_count, set_end, set_offset_end, set_offset_start, set_start, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: eri_offsets
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: map_it_here
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: hab_full_set
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sphi_a, sphi_b, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      G_count = 0; R_count = 0; 
      ! get mapping between ERIs and atoms, sets, set offsets
      CALL get_eri_offsets(qs_env, basis_type_b, eri_offsets)

      atom_start = eri_offsets(first_b, 1)
      set_start = eri_offsets(first_b, 2)
      set_offset_start = eri_offsets(first_b, 3)

      atom_end = eri_offsets(last_b, 1)
      set_end = eri_offsets(last_b, 2)
      set_offset_end = eri_offsets(last_b, 3)

      ! get QS stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set)
      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)
      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      ! get total number of local set pairs to integrate
      nseta_total = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type_a)
         nseta_total = nseta_total+basis_set_a%nset
      ENDDO

      nsetb_total = 0
      DO jatom = atom_start, atom_end
         jkind = kind_of(jatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type=basis_type_b)
         nsetb_total = nsetb_total+basis_set_b%nset
      ENDDO

      n_setpair = nseta_total*nsetb_total

      my_setpair = 0

      offset_a_end = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type_a)

         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         npgfa => basis_set_a%npgf

         ra(:) = pbc(particle_set(iatom)%r, cell)

         DO iset = 1, nseta
            offset_a_start = offset_a_end
            offset_a_end = offset_a_end+nsgfa(iset)
            sgfa = first_sgfa(1, iset)

            offset_b_end = 0
            DO jatom = atom_start, atom_end
               jkind = kind_of(jatom)
               CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type=basis_type_b)

               first_sgfb => basis_set_b%first_sgf
               lb_max => basis_set_b%lmax
               lb_min => basis_set_b%lmin
               nsetb = basis_set_b%nset
               nsgfb => basis_set_b%nsgf_set
               sphi_b => basis_set_b%sphi
               zetb => basis_set_b%zet
               npgfb => basis_set_b%npgf

               rb(:) = pbc(particle_set(jatom)%r, cell)
               rab(:) = ra(:)-rb(:) ! pbc not needed?
               dab = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)

               jset_start = 1; jset_end = nsetb
               IF (jatom == atom_start) jset_start = set_start
               IF (jatom == atom_end) jset_end = set_end

               DO jset = jset_start, jset_end
                  first_set = 1; last_set = nsgfb(jset)
                  IF (jset == jset_start .AND. jatom == atom_start) first_set = set_offset_start
                  IF (jset == jset_end .AND. jatom == atom_end) last_set = set_offset_end

                  offset_b_start = offset_b_end
                  offset_b_end = offset_b_end+last_set+1-first_set
                  sgfb = first_sgfb(1, jset)
                  my_setpair = my_setpair+1
                  map_it_here = MODULO(my_setpair, para_env%num_pe) == para_env%mepos

                  IF (map_it_here) THEN
                     ALLOCATE (hab_full_set(nsgfa(iset), nsgfb(jset)))
                     CALL eri_mme_2c_integrate_set(param%par, &
                                                   la_min(iset), la_max(iset), lb_min(jset), lb_max(jset), &
                                                   npgfa(iset), zeta(:, iset), npgfb(jset), zetb(:, jset), &
                                                   ra, rb, &
                                                   hab_full_set, 1, 1, &
                                                   sphi_a, sgfa, nsgfa(iset), &
                                                   sphi_b, sgfb, nsgfb(jset), G_count, R_count)

                     hab(offset_a_start+1:offset_a_end, offset_b_start+1:offset_b_end) = &
                        hab_full_set(:, first_set:last_set)
                     DEALLOCATE (hab_full_set)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL update_local_counts(param, para_env, G_count_2c=G_count, R_count_2c=R_count)

      CALL mp_sum(hab, para_env%group)

      CALL timestop(handle)
   END SUBROUTINE cp_eri_mme_2c_integrate

! **************************************************************************************************
!> \brief high-level integration routine for 2c integrals over s-type functions.
!>        Parallelization over pairs of functions.
!> \param param ...
!> \param zeta ...
!> \param zetb ...
!> \param ra ...
!> \param rb ...
!> \param hab ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_2c_integrate_s(param, zeta, zetb, ra, rb, hab, para_env)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_eri_mme_2c_integrate_s', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: G_count, handle, ipgf, ipgf_prod, jpgf, &
                                                            npgf_prod, npgfa, npgfb, R_count
      INTEGER, DIMENSION(2)                              :: limits
      REAL(KIND=dp), DIMENSION(3)                        :: rab

      CALL timeset(routineN, handle)
      G_count = 0; R_count = 0

      hab(:, :) = 0.0_dp

      npgfa = SIZE(zeta)
      npgfb = SIZE(zetb)
      npgf_prod = npgfa*npgfb ! total number of integrals

      limits = get_limit(npgf_prod, para_env%num_pe, para_env%mepos)

      DO ipgf_prod = limits(1), limits(2)
         ipgf = (ipgf_prod-1)/npgfb+1
         jpgf = MOD(ipgf_prod-1, npgfb)+1
         rab(:) = ra(:, ipgf)-rb(:, jpgf)
         CALL eri_mme_2c_integrate(param%par, 0, 0, 0, 0, zeta(ipgf), &
                                   zetb(jpgf), rab, hab, ipgf-1, jpgf-1, G_count, R_count)
      ENDDO

      CALL update_local_counts(param, para_env, G_count_2c=G_count, R_count_2c=R_count)
      CALL mp_sum(hab, para_env%group)
      CALL timestop(handle)

   END SUBROUTINE cp_eri_mme_2c_integrate_s

! **************************************************************************************************
!> \brief high-level integration routine for 3c integrals (ab|c) over CP2K basis sets.
!>        For each local function of c, (ab|c) is written to a DBCSR matrix mat_ab.
!> \param param ...
!> \param para_env ...
!> \param qs_env ...
!> \param first_c start index of local range of c
!> \param last_c end index of local range of c
!> \param mat_ab DBCSR matrices for each c
!> \param basis_type_a ...
!> \param basis_type_b ...
!> \param basis_type_c ...
!> \param sab_nl neighbor list for a, b
!> \param transp return (c|ab) instead of (ab|c)
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_3c_integrate(param, para_env, qs_env, &
                                      first_c, last_c, mat_ab, &
                                      basis_type_a, basis_type_b, basis_type_c, &
                                      sab_nl, transp)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      INTEGER, INTENT(IN)                                :: first_c, last_c
      TYPE(cp_dbcsr_p_type), &
         DIMENSION(last_c-first_c+1), INTENT(INOUT)      :: mat_ab
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type_a, basis_type_b, basis_type_c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      LOGICAL, INTENT(IN), OPTIONAL                      :: transp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_eri_mme_3c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, atom_b, atom_end, atom_start, GG_count, GR_count, handle, iatom, ic, &
         ic_dbcsr, icol, ikind, inode, irow, iset, jatom, jkind, jset, katom, kkind, kset, &
         kset_end, kset_start, last_jatom, mepos, na, natom, nb, nc, nkind, nseta, nsetb, nsetc, &
         nthread, offset_a_end, offset_a_start, offset_b_end, offset_b_start, offset_c_end, &
         offset_c_start, RR_count, set_end, set_offset_end, set_offset_start, set_start, sgfa, &
         sgfb, sgfc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: eri_offsets
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lc_max, &
                                                            lc_min, npgfa, npgfb, npgfc, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfc
      LOGICAL                                            :: do_symmetric, found, to_be_asserted
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: habc
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: munu_block, rpgfb, sphi_a, sphi_b, &
                                                            sphi_c, zeta, zetb, zetc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_c
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      GG_count = 0; GR_count = 0; RR_count = 0

      nthread = 1

      ! get mapping between ERIs and atoms, sets, set offsets
      CALL get_eri_offsets(qs_env, basis_type_c, eri_offsets)
      atom_start = eri_offsets(first_c, 1)
      set_start = eri_offsets(first_c, 2)
      set_offset_start = eri_offsets(first_c, 3)

      atom_end = eri_offsets(last_c, 1)
      set_end = eri_offsets(last_c, 2)
      set_offset_end = eri_offsets(last_c, 3)

      ! get QS stuff
      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      natom=natom, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      cell=cell)

      nkind = SIZE(qs_kind_set)
      ALLOCATE (atom_of_kind(natom))
      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)

      ! get total number of local c functions to integrate (full sets)
      nc = 0
      DO katom = atom_start, atom_end
         kkind = kind_of(katom)
         CALL get_qs_kind(qs_kind=qs_kind_set(kkind), basis_set=basis_set_c, basis_type=basis_type_c)
         kset_start = 1; kset_end = basis_set_c%nset
         IF (katom == atom_start) kset_start = set_start
         IF (katom == atom_end) kset_end = set_end
         nc = nc+SUM(basis_set_c%nsgf_set(kset_start:kset_end))
      ENDDO

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)
      IF (do_symmetric) THEN
         CPASSERT(basis_type_a == basis_type_b)
      END IF

      ALLOCATE (basis_set_list_a(nkind), basis_set_list_b(nkind))
      CALL basis_set_list_setup(basis_set_list_a, basis_type_a, qs_kind_set)
      CALL basis_set_list_setup(basis_set_list_b, basis_type_b, qs_kind_set)

      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

      mepos = 0

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)

         ! exclude periodic images because method is periodic intrinsically
         IF (inode == 1) last_jatom = 0

         IF (jatom /= last_jatom) THEN
            last_jatom = jatom
         ELSE
            CYCLE
         END IF

         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)

         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         na = SUM(nsgfa)

         ra(:) = pbc(particle_set(iatom)%r, cell)

         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         nb = SUM(nsgfb)

         rb(:) = pbc(particle_set(jatom)%r, cell)

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         ALLOCATE (habc(na, nb, nc))
         habc(:, :, :) = 0.0_dp ! needs to be initialized due to screening

         rab(:) = pbc(rab, cell)
         dab = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)

         offset_a_end = 0
         DO iset = 1, nseta
            offset_a_start = offset_a_end
            offset_a_end = offset_a_end+nsgfa(iset)
            sgfa = first_sgfa(1, iset)

            offset_b_end = 0
            DO jset = 1, nsetb
               offset_b_start = offset_b_end
               offset_b_end = offset_b_end+nsgfb(jset)

               sgfb = first_sgfb(1, jset)

               ! Screening
               IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

               offset_c_end = 0
               DO katom = atom_start, atom_end

                  kkind = kind_of(katom)
                  CALL get_qs_kind(qs_kind=qs_kind_set(kkind), basis_set=basis_set_c, basis_type=basis_type_c)
                  first_sgfc => basis_set_c%first_sgf
                  lc_max => basis_set_c%lmax
                  lc_min => basis_set_c%lmin
                  nsetc = basis_set_c%nset
                  nsgfc => basis_set_c%nsgf_set
                  sphi_c => basis_set_c%sphi
                  zetc => basis_set_c%zet
                  npgfc => basis_set_c%npgf

                  rc(:) = pbc(particle_set(katom)%r, cell)

                  kset_start = 1; kset_end = nsetc
                  IF (katom == atom_start) kset_start = set_start
                  IF (katom == atom_end) kset_end = set_end

                  DO kset = kset_start, kset_end
                     offset_c_start = offset_c_end
                     offset_c_end = offset_c_end+nsgfc(kset)
                     sgfc = first_sgfc(1, kset)

                     CALL eri_mme_3c_integrate_set(param%par, &
                                                   la_min(iset), la_max(iset), &
                                                   lb_min(jset), lb_max(jset), &
                                                   lc_min(kset), lc_max(kset), &
                                                   npgfa(iset), npgfb(jset), npgfc(kset), &
                                                   zeta(:, iset), zetb(:, jset), zetc(:, kset), &
                                                   ra, rb, rc, &
                                                   habc, &
                                                   offset_a_start+1, &
                                                   offset_b_start+1, &
                                                   offset_c_start+1, &
                                                   sphi_a, sphi_b, sphi_c, &
                                                   sgfa, sgfb, sgfc, &
                                                   nsgfa(iset), nsgfb(jset), nsgfc(kset), &
                                                   GG_count, GR_count, RR_count, &
                                                   transp=transp)

                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         ic_dbcsr = 0
         DO ic = set_offset_start, set_offset_start+last_c-first_c
            ic_dbcsr = ic_dbcsr+1
            NULLIFY (munu_block)
            CALL cp_dbcsr_get_block_p(matrix=mat_ab(ic_dbcsr)%matrix, &
                                      row=irow, col=icol, block=munu_block, found=found)
            CPASSERT(found)
            munu_block(:, :) = 0.0_dp
            IF (irow .EQ. iatom) THEN
               to_be_asserted = SIZE(munu_block, 1) .EQ. SIZE(habc, 1) .AND. SIZE(munu_block, 2) .EQ. SIZE(habc, 2)
               CPASSERT(to_be_asserted)
               munu_block(:, :) = habc(:, :, ic)
            ELSE
               to_be_asserted = SIZE(munu_block, 2) .EQ. SIZE(habc, 1) .AND. SIZE(munu_block, 1) .EQ. SIZE(habc, 2)
               CPASSERT(to_be_asserted)
               munu_block(:, :) = TRANSPOSE(habc(:, :, ic))
            ENDIF
         ENDDO
         DEALLOCATE (habc)
      ENDDO

      DEALLOCATE (basis_set_list_a, basis_set_list_b)
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL update_local_counts(param, para_env, GG_count_3c=GG_count, GR_count_3c=GR_count, RR_count_3c=RR_count)

      CALL timestop(handle)
   END SUBROUTINE cp_eri_mme_3c_integrate

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param iw ...
!> \param eri_mme_test_section ...
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_perf_acc_test(para_env, iw, eri_mme_test_section)
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: iw
      TYPE(section_vals_type), POINTER                   :: eri_mme_test_section

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_eri_mme_perf_acc_test', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: count_r, G_count, GG_count, GR_count, i, &
                                                            ix, iy, iz, l_max, min_nR, nR, nR_xyz, &
                                                            nrep, nsample, nzet, R_count, RR_count
      LOGICAL                                            :: test_2c, test_3c, test_accuracy
      REAL(KIND=dp)                                      :: zet_fac, zetmax, zetmin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rabc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cell_par
      TYPE(cell_type), POINTER                           :: box
      TYPE(cp_eri_mme_param)                             :: param
      TYPE(section_vals_type), POINTER                   :: eri_mme_section

      NULLIFY (box, eri_mme_section, cell_par)

      eri_mme_section => section_vals_get_subs_vals(eri_mme_test_section, "ERI_MME")
      CALL cp_eri_mme_init_read_input(eri_mme_section, param)
      CALL section_vals_val_get(eri_mme_test_section, "TEST_3C", l_val=test_3c)
      CALL section_vals_val_get(eri_mme_test_section, "TEST_2C", l_val=test_2c)

      CALL section_vals_val_get(eri_mme_test_section, "ABC", r_vals=cell_par)
      CALL section_vals_val_get(eri_mme_test_section, "MIN_NPOS", i_val=min_nR)
      CALL section_vals_val_get(eri_mme_test_section, "NREP", i_val=nrep)
      CALL section_vals_val_get(eri_mme_test_section, "CHECK_2C_ACCURACY", l_val=test_accuracy)
      CALL section_vals_val_get(eri_mme_test_section, "LMAX", i_val=l_max)
      CALL section_vals_val_get(eri_mme_test_section, "ZET_MIN", r_val=zetmin)
      CALL section_vals_val_get(eri_mme_test_section, "ZET_MAX", r_val=zetmax)
      CALL section_vals_val_get(eri_mme_test_section, "NZET", i_val=nzet)
      CALL section_vals_val_get(eri_mme_test_section, "NSAMPLE_3C", i_val=nsample)

      IF (nzet .LE. 0) &
         CPABORT("Number of exponents NZET must be greater than 0.")

      CALL init_orbital_pointers(l_max)

      ! Create ranges of zet to be tested
      ALLOCATE (zet(nzet))

      zet(1) = zetmin
      IF (nzet .GT. 1) THEN
         zet_fac = (zetmax/zetmin)**(1.0_dp/(nzet-1))
         DO i = 1, nzet-1
            zet(i+1) = zet(i)*zet_fac
         ENDDO
      ENDIF

      ! initialize cell
      CALL cell_create(box)
      box%hmat = 0.0_dp
      box%hmat(1, 1) = cell_par(1)
      box%hmat(2, 2) = cell_par(2)
      box%hmat(3, 3) = cell_par(3)
      CALL init_cell(box)

      ! Create range of rab (atomic distances) to be tested
      nR_xyz = CEILING(REAL(min_nR, KIND=dp)**(1.0_dp/3.0_dp)-1.0E-06)
      nR = nR_xyz**3

      ALLOCATE (rabc(3, nR))
      count_r = 0
      DO ix = 1, nR_xyz
      DO iy = 1, nR_xyz
      DO iz = 1, nR_xyz
         count_r = count_r+1
         ! adding 10% of cell size to positions to avoid atoms exactly at boundary or center of a cell
         rabc(:, count_r) = pbc([ix*ABS(cell_par(1)), &
                                 iy*ABS(cell_par(2)), &
                                 iz*ABS(cell_par(3))]/nR_xyz+ &
                                0.1_dp*ABS(cell_par(:)), box)
      ENDDO
      ENDDO
      ENDDO

      ! initialize MME method
      CALL cp_eri_mme_set_params(param, box%hmat, box%orthorhombic, MINVAL(zet), [MAXVAL(zet)], [l_max], l_max, para_env)

      IF (iw > 0) WRITE (iw, '(T2, A, T61, I20)') "ERI_MME| Number of atomic distances:", nR

      G_count = 0; R_count = 0
      GG_count = 0; GR_count = 0; RR_count = 0

      IF (test_2c) CALL eri_mme_2c_perf_acc_test(param%par, l_max, zet, rabc, nrep, test_accuracy, para_env, iw, G_count, R_count)
      IF (test_3c) CALL eri_mme_3c_perf_acc_test(param%par, l_max, zet, rabc, nrep, nsample, &
                                                 para_env, iw, GG_count, GR_count, RR_count)
      CALL update_local_counts(param, para_env, G_count, R_count, GG_count, GR_count, RR_count)
      CALL cp_eri_mme_finalize(param)
      CALL cell_release(box)
   END SUBROUTINE cp_eri_mme_perf_acc_test

! **************************************************************************************************
!> \brief Update local counters to gather statistics on different paths taken in MME algorithm
!> (each Ewald sum can be performed over direct or reciprocal lattice vectors)
!> \param param ...
!> \param para_env ...
!> \param G_count_2c ...
!> \param R_count_2c ...
!> \param GG_count_3c ...
!> \param GR_count_3c ...
!> \param RR_count_3c ...
! **************************************************************************************************
   SUBROUTINE update_local_counts(param, para_env, G_count_2c, R_count_2c, GG_count_3c, GR_count_3c, RR_count_3c)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(INOUT), OPTIONAL                   :: G_count_2c, R_count_2c, GG_count_3c, &
                                                            GR_count_3c, RR_count_3c

      IF (PRESENT(G_count_2c)) THEN
         CALL mp_sum(G_count_2c, para_env%group)
         param%G_count_2c = param%G_count_2c+G_count_2c
      ENDIF

      IF (PRESENT(R_count_2c)) THEN
         CALL mp_sum(R_count_2c, para_env%group)
         param%R_count_2c = param%R_count_2c+R_count_2c
      ENDIF

      IF (PRESENT(GG_count_3c)) THEN
         CALL mp_sum(GG_count_3c, para_env%group)
         param%GG_count_3c = param%GG_count_3c+GG_count_3c
      ENDIF

      IF (PRESENT(GR_count_3c)) THEN
         CALL mp_sum(GR_count_3c, para_env%group)
         param%GR_count_3c = param%GR_count_3c+GR_count_3c
      ENDIF

      IF (PRESENT(RR_count_3c)) THEN
         CALL mp_sum(RR_count_3c, para_env%group)
         param%RR_count_3c = param%RR_count_3c+RR_count_3c
      ENDIF

   END SUBROUTINE update_local_counts

! **************************************************************************************************
!> \brief get pointer to atom, pointer to set and offset in a set for each spherical orbital of a
!>        basis.
!> \param qs_env ...
!> \param basis_type ...
!> \param eri_offsets (:,1) atom numbers
!>                    (:,2) set numbers
!>                    (:,3) set offsets
! **************************************************************************************************
   SUBROUTINE get_eri_offsets(qs_env, basis_type, eri_offsets)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: eri_offsets

      INTEGER                                            :: dimen_basis, iatom, ikind, iset, isgf, &
                                                            natom, nkind, nset, nsgf, offset, &
                                                            set_offset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(:), POINTER                     :: nsgf_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)
      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      dimen_basis = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type=basis_type)
         dimen_basis = dimen_basis+nsgf
      END DO

      ALLOCATE (eri_offsets(dimen_basis, 3))

      offset = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set, basis_type=basis_type)
         nset = basis_set%nset
         nsgf_set => basis_set%nsgf_set
         DO iset = 1, nset
            set_offset = 0
            DO isgf = 1, nsgf_set(iset)
               set_offset = set_offset+1
               eri_offsets(offset+set_offset, :) = [iatom, iset, set_offset]
            ENDDO
            offset = offset+nsgf_set(iset)
         ENDDO
      ENDDO
   END SUBROUTINE get_eri_offsets

END MODULE cp_eri_mme_interface
