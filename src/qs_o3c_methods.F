!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Methods used with 3-center overlap type integrals containers
!> \par History
!>      - none
!>      - 11.2018 fixed OMP race condition in contract3_o3c routine (A.Bussy)
!>      - 05.2019 Added a routine to compute 3-center integrals with libint (A.Bussy)
! **************************************************************************************************
MODULE qs_o3c_methods
   USE ai_contraction_sphi,             ONLY: abc_contract
   USE ai_overlap3,                     ONLY: overlap3
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE gamma,                           ONLY: init_md_ftable
   USE hfx_libint_wrapper,              ONLY: cp_libint_cleanup_3eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE input_constants,                 ONLY: do_potential_truncated
   USE kinds,                           ONLY: dp
   USE libint_3center,                  ONLY: eri_3center
   USE orbital_pointers,                ONLY: ncoset
   USE qs_o3c_types,                    ONLY: &
        get_o3c_container, get_o3c_iterator_info, get_o3c_vec, o3c_container_type, o3c_iterate, &
        o3c_iterator_create, o3c_iterator_release, o3c_iterator_type, o3c_vec_type, &
        set_o3c_container
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_o3c_methods'

   PUBLIC :: calculate_o3c_integrals, contract12_o3c, contract3_o3c, &
             calculate_o3c_libint_integrals

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param o3c ...
!> \param calculate_forces ...
!> \param matrix_p ...
! **************************************************************************************************
   SUBROUTINE calculate_o3c_integrals(o3c, calculate_forces, matrix_p)
      TYPE(o3c_container_type), POINTER                  :: o3c
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_p

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_o3c_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER :: egfa, egfb, egfc, handle, i, iatom, icol, ikind, irow, iset, ispin, j, jatom, &
         jkind, jset, katom, kkind, kset, mepos, ncoa, ncob, ncoc, ni, nj, nk, nseta, nsetb, &
         nsetc, nspin, nthread, sgfa, sgfb, sgfc
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lc_max, &
                                                            lc_min, npgfa, npgfb, npgfc, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfc
      LOGICAL                                            :: do_force, found, trans
      REAL(KIND=dp)                                      :: dij, dik, djk, fpre
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pmat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: sabc, sabc_contr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: iabdc, iadbc, idabc, sabdc, sdabc
      REAL(KIND=dp), DIMENSION(3)                        :: rij, rik, rjk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fi, fj, fk, pblock, rpgfa, rpgfb, rpgfc, &
                                                            sphi_a, sphi_b, sphi_c, tvec, zeta, &
                                                            zetb, zetc
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: iabc
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b, &
                                                            basis_set_list_c
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_c
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      do_force = .FALSE.
      IF (PRESENT(calculate_forces)) do_force = calculate_forces
      CALL get_o3c_container(o3c, nspin=nspin)

      ! basis sets
      CALL get_o3c_container(o3c, basis_set_list_a=basis_set_list_a, &
                             basis_set_list_b=basis_set_list_b, basis_set_list_c=basis_set_list_c)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,ncoset,nspin,basis_set_list_a,basis_set_list_b,&
!$OMP         basis_set_list_c,do_force,matrix_p)&
!$OMP PRIVATE (mepos,ikind,jkind,kkind,basis_set_a,basis_set_b,basis_set_c,rij,rik,rjk,&
!$OMP          first_sgfa,la_max,la_min,npgfa,nseta,nsgfa,rpgfa,set_radius_a,sphi_a,zeta,&
!$OMP          first_sgfb,lb_max,lb_min,npgfb,nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,&
!$OMP          first_sgfc,lc_max,lc_min,npgfc,nsetc,nsgfc,rpgfc,set_radius_c,sphi_c,zetc,&
!$OMP          iset,jset,kset,dij,dik,djk,ni,nj,nk,iabc,idabc,iadbc,iabdc,tvec,fi,fj,fk,ncoa,&
!$OMP          ncob,ncoc,sabc,sabc_contr,sdabc,sabdc,sgfa,sgfb,sgfc,egfa,egfb,egfc,i,j,&
!$OMP          pblock,pmat,ispin,iatom,jatom,katom,irow,icol,found,trans,fpre)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, &
                                    ikind=ikind, jkind=jkind, kkind=kkind, rij=rij, rik=rik, &
                                    integral=iabc, tvec=tvec, force_i=fi, force_j=fj, force_k=fk)
         CPASSERT(.NOT. ASSOCIATED(iabc))
         CPASSERT(.NOT. ASSOCIATED(tvec))
         CPASSERT(.NOT. ASSOCIATED(fi))
         CPASSERT(.NOT. ASSOCIATED(fj))
         CPASSERT(.NOT. ASSOCIATED(fk))
         ! basis
         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         basis_set_c => basis_set_list_c(kkind)%gto_basis_set
         ! center A
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! center B
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
         ! center C (RI)
         first_sgfc => basis_set_c%first_sgf
         lc_max => basis_set_c%lmax
         lc_min => basis_set_c%lmin
         npgfc => basis_set_c%npgf
         nsetc = basis_set_c%nset
         nsgfc => basis_set_c%nsgf_set
         rpgfc => basis_set_c%pgf_radius
         set_radius_c => basis_set_c%set_radius
         sphi_c => basis_set_c%sphi
         zetc => basis_set_c%zet

         ni = SUM(nsgfa)
         nj = SUM(nsgfb)
         nk = SUM(nsgfc)

         ALLOCATE (iabc(ni, nj, nk))
         iabc(:, :, :) = 0.0_dp
         IF (do_force) THEN
            ALLOCATE (fi(nk, 3), fj(nk, 3), fk(nk, 3))
            fi(:, :) = 0.0_dp
            fj(:, :) = 0.0_dp
            fk(:, :) = 0.0_dp
            ALLOCATE (idabc(ni, nj, nk, 3))
            idabc(:, :, :, :) = 0.0_dp
            ALLOCATE (iadbc(ni, nj, nk, 3))
            iadbc(:, :, :, :) = 0.0_dp
            ALLOCATE (iabdc(ni, nj, nk, 3))
            iabdc(:, :, :, :) = 0.0_dp
         ELSE
            NULLIFY (fi, fj, fk)
         END IF
         ALLOCATE (tvec(nk, nspin))
         tvec(:, :) = 0.0_dp

         rjk(1:3) = rik(1:3)-rij(1:3)
         dij = SQRT(SUM(rij**2))
         dik = SQRT(SUM(rik**2))
         djk = SQRT(SUM(rjk**2))

         DO iset = 1, nseta
            DO jset = 1, nsetb
               IF (set_radius_a(iset)+set_radius_b(jset) < dij) CYCLE
               DO kset = 1, nsetc
                  IF (set_radius_a(iset)+set_radius_c(kset) < dik) CYCLE
                  IF (set_radius_b(jset)+set_radius_c(kset) < djk) CYCLE

                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))

                  sgfa = first_sgfa(1, iset)
                  sgfb = first_sgfb(1, jset)
                  sgfc = first_sgfc(1, kset)

                  egfa = sgfa+nsgfa(iset)-1
                  egfb = sgfb+nsgfb(jset)-1
                  egfc = sgfc+nsgfc(kset)-1

                  IF (ncoa*ncob*ncoc > 0) THEN
                     ALLOCATE (sabc(ncoa, ncob, ncoc))
                     sabc(:, :, :) = 0.0_dp
                     IF (do_force) THEN
                        ALLOCATE (sdabc(ncoa, ncob, ncoc, 3))
                        sdabc(:, :, :, :) = 0.0_dp
                        ALLOCATE (sabdc(ncoa, ncob, ncoc, 3))
                        sabdc(:, :, :, :) = 0.0_dp
                        CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                      lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                      lc_max(kset), npgfc(kset), zetc(:, kset), rpgfc(:, kset), lc_min(kset), &
                                      rij, dij, rik, dik, rjk, djk, sabc, sdabc, sabdc)
                     ELSE
                        CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                      lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                      lc_max(kset), npgfc(kset), zetc(:, kset), rpgfc(:, kset), lc_min(kset), &
                                      rij, dij, rik, dik, rjk, djk, sabc)
                     END IF
                     ALLOCATE (sabc_contr(nsgfa(iset), nsgfb(jset), nsgfc(kset)))

                     CALL abc_contract(sabc_contr, sabc, &
                                       sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                       ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))
                     iabc(sgfa:egfa, sgfb:egfb, sgfc:egfc) = &
                        sabc_contr(1:nsgfa(iset), 1:nsgfb(jset), 1:nsgfc(kset))
                     IF (do_force) THEN
                        DO i = 1, 3
                           CALL abc_contract(sabc_contr, sdabc(:, :, :, i), &
                                             sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                             ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))
                           idabc(sgfa:egfa, sgfb:egfb, sgfc:egfc, i) = &
                              sabc_contr(1:nsgfa(iset), 1:nsgfb(jset), 1:nsgfc(kset))
                           CALL abc_contract(sabc_contr, sabdc(:, :, :, i), &
                                             sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                             ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))
                           iabdc(sgfa:egfa, sgfb:egfb, sgfc:egfc, i) = &
                              sabc_contr(1:nsgfa(iset), 1:nsgfb(jset), 1:nsgfc(kset))
                        END DO
                     END IF

                     DEALLOCATE (sabc_contr)
                     DEALLOCATE (sabc)
                  END IF
                  IF (do_force) THEN
                     DEALLOCATE (sdabc, sabdc)
                  END IF
               END DO
            END DO
         END DO
         IF (do_force) THEN
            ! translational invariance
            iadbc(:, :, :, :) = -idabc(:, :, :, :)-iabdc(:, :, :, :)
            !
            ! get the atom indices
            CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, &
                                       iatom=iatom, jatom=jatom, katom=katom)
            !
            ! contract over i and j to get forces
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
               trans = .FALSE.
            ELSE
               irow = jatom
               icol = iatom
               trans = .TRUE.
            END IF
            IF (iatom == jatom) THEN
               fpre = 1.0_dp
            ELSE
               fpre = 2.0_dp
            END IF
            ALLOCATE (pmat(ni, nj))
            pmat(:, :) = 0.0_dp
            DO ispin = 1, nspin
               CALL dbcsr_get_block_p(matrix=matrix_p(ispin)%matrix, &
                                      row=irow, col=icol, BLOCK=pblock, found=found)
               IF (found) THEN
                  IF (trans) THEN
                     pmat(:, :) = pmat(:, :)+TRANSPOSE(pblock(:, :))
                  ELSE
                     pmat(:, :) = pmat(:, :)+pblock(:, :)
                  END IF
               END IF
            END DO
            DO i = 1, 3
               DO j = 1, nk
                  fi(j, i) = fpre*SUM(pmat(:, :)*idabc(:, :, j, i))
                  fj(j, i) = fpre*SUM(pmat(:, :)*iadbc(:, :, j, i))
                  fk(j, i) = fpre*SUM(pmat(:, :)*iabdc(:, :, j, i))
               END DO
            END DO
            DEALLOCATE (pmat)
            !
            DEALLOCATE (idabc, iadbc, iabdc)
         END IF
         !
         CALL set_o3c_container(o3c_iterator, mepos=mepos, &
                                integral=iabc, tvec=tvec, force_i=fi, force_j=fj, force_k=fk)

      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE calculate_o3c_integrals

! **************************************************************************************************
!> \brief Computes the 3-center integrals of the o3c container based on libint for the given operator
!> \param o3c the 3-center integrals container
!> \param op the operator for the integrals
!> \param t_c_filename the name of the file with the truncated Coulomb data
!> \param para_env ...
!> \param r_cutoff the cutoff radius for the truncated coulomb operator
!> \param omega the range parameter for the erfc short range operator
!> \param eps_screen the screening threshold for sicarding integrals before contraction
!> \note The static initialization of the libint library needs to be done before hand
!>       In case the truncated coulomb operator is used, the potential parameter file must be read
!>       in advance too
! **************************************************************************************************
   SUBROUTINE calculate_o3c_libint_integrals(o3c, op, t_c_filename, para_env, r_cutoff, omega, &
                                             eps_screen)

      TYPE(o3c_container_type), POINTER                  :: o3c
      INTEGER, INTENT(IN)                                :: op
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: t_c_filename
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env
      REAL(dp), INTENT(IN), OPTIONAL                     :: r_cutoff, omega, eps_screen

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_o3c_libint_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER :: egfa, egfb, egfc, handle, ibasis, ikind, imax, iset, jkind, jset, kkind, kset, &
         m_max, maxli, maxlj, maxlk, mepos, ncoa, ncob, ncoc, ni, nj, nk, nseta, nsetb, nsetc, &
         nspin, nthread, sgfa, sgfb, sgfc, unit_id
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lc_max, &
                                                            lc_min, npgfa, npgfb, npgfc, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfc
      LOGICAL                                            :: do_screen
      REAL(dp)                                           :: my_eps_screen, my_omega, my_r_cutoff
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sabc
      REAL(dp), DIMENSION(3)                             :: ri, rij, rik, rj, rk
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi_a, sphi_b, sphi_c, tvec, zeta, &
                                                            zetb, zetc
      REAL(dp), DIMENSION(:, :, :), POINTER              :: iabc
      TYPE(cp_libint_t)                                  :: lib
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b, &
                                                            basis_set_list_c
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_c
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      NULLIFY (basis_set_list_a, basis_set_list_b, basis_set_list_c, basis_set_a, basis_set_b)
      NULLIFY (basis_set_c, iabc, tvec, first_sgfa, first_sgfb, first_sgfc, la_max, la_min, lb_max)
      NULLIFY (lb_min, lc_max, lc_min, npgfa, npgfb, npgfc, nsgfa, nsgfb, nsgfc)

      CALL timeset(routineN, handle)

      CALL get_o3c_container(o3c, nspin=nspin, basis_set_list_a=basis_set_list_a, &
                             basis_set_list_b=basis_set_list_b, basis_set_list_c=basis_set_list_c)

      !Need the max l for each basis for libint
      maxli = 0
      DO ibasis = 1, SIZE(basis_set_list_a)
         CALL get_gto_basis_set(gto_basis_set=basis_set_list_a(ibasis)%gto_basis_set, maxl=imax)
         maxli = MAX(maxli, imax)
      END DO
      maxlj = 0
      DO ibasis = 1, SIZE(basis_set_list_b)
         CALL get_gto_basis_set(gto_basis_set=basis_set_list_b(ibasis)%gto_basis_set, maxl=imax)
         maxlj = MAX(maxlj, imax)
      END DO
      maxlk = 0
      DO ibasis = 1, SIZE(basis_set_list_c)
         CALL get_gto_basis_set(gto_basis_set=basis_set_list_c(ibasis)%gto_basis_set, maxl=imax)
         maxlk = MAX(maxlk, imax)
      END DO
      m_max = maxli+maxlj+maxlk

      !Screening
      do_screen = .FALSE.
      IF (PRESENT(eps_screen)) THEN
         do_screen = .TRUE.
         my_eps_screen = eps_screen
      END IF

      !Short range parameters
      my_omega = 0.0_dp
      IF (PRESENT(omega)) my_omega = omega
      my_r_cutoff = 0.0_dp
      IF (PRESENT(r_cutoff)) my_r_cutoff = r_cutoff

      !Init the truncated Coulomb operator
      IF (op == do_potential_truncated) THEN
         CPASSERT(PRESENT(t_c_filename))
         CPASSERT(PRESENT(para_env))

         !open the file only if necessary
         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=t_c_filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      !Inint the initial gamma function before the OMP region as it is not thread safe
      CALL init_md_ftable(nmax=m_max)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,nspin,basis_set_list_a, basis_set_list_b,basis_set_list_c, &
!$OMP         maxli,maxlj,maxlk,my_r_cutoff,my_omega,ncoset,op,do_screen,my_eps_screen) &
!$OMP PRIVATE (basis_set_a,basis_set_b,basis_set_c,mepos,ikind,jkind,kkind,iabc,tvec,rij,rik, &
!$OMP          first_sgfa,la_max,la_min,npgfa,nseta,nsgfa,zeta,iset,ni,ri,ncoa,sgfa,egfa,sphi_a,&
!$OMP          first_sgfb,lb_max,lb_min,npgfb,nsetb,nsgfb,zetb,jset,nj,rj,ncob,sgfb,egfb,sphi_b,&
!$OMP          first_sgfc,lc_max,lc_min,npgfc,nsetc,nsgfc,zetc,kset,nk,rk,ncoc,sgfc,egfc,sphi_c,&
!$OMP          sabc,lib)

      mepos = 0
!$    mepos = omp_get_thread_num()

      !each thread needs its own lib because internal parameters could change at different rates
      CALL cp_libint_init_3eri(lib, MAX(maxli, maxlj, maxlk))
      CALL cp_libint_set_contrdepth(lib, 1)

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)

         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                    kkind=kkind, rij=rij, rik=rik, integral=iabc, tvec=tvec)

         !basis
         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         basis_set_c => basis_set_list_c(kkind)%gto_basis_set
         ! center A
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! center B
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         ! center C
         first_sgfc => basis_set_c%first_sgf
         lc_max => basis_set_c%lmax
         lc_min => basis_set_c%lmin
         npgfc => basis_set_c%npgf
         nsetc = basis_set_c%nset
         nsgfc => basis_set_c%nsgf_set
         sphi_c => basis_set_c%sphi
         zetc => basis_set_c%zet

         ni = SUM(nsgfa)
         nj = SUM(nsgfb)
         nk = SUM(nsgfc)

         ALLOCATE (iabc(ni, nj, nk))
         iabc(:, :, :) = 0.0_dp

         ALLOCATE (tvec(nk, nspin))
         tvec(:, :) = 0.0_dp

         !need positions for libint. Only relative positions are needed => set ri to 0.0
         ri = 0.0_dp
         rj = rij ! ri + rij
         rk = rik ! ri + rik

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            egfa = sgfa+nsgfa(iset)-1

            DO jset = 1, nsetb
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)
               egfb = sgfb+nsgfb(jset)-1

               DO kset = 1, nsetc
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))
                  sgfc = first_sgfc(1, kset)
                  egfc = sgfc+nsgfc(kset)-1

                  ALLOCATE (sabc(ncoa, ncob, ncoc))
                  sabc = 0.0_dp

                  CALL eri_3center(sabc, la_min(iset), la_max(iset), npgfa(iset), zeta(:, iset), &
                                   ri, lb_min(jset), lb_max(jset), npgfb(jset), zetb(:, jset), &
                                   rj, lc_min(kset), lc_max(kset), npgfc(kset), zetc(:, kset), &
                                   rk, lib, op, omega=my_omega, r_cutoff=my_r_cutoff)

                  IF (do_screen) THEN
                     IF (MAXVAL(ABS(sabc)) < my_eps_screen) THEN
                        DEALLOCATE (sabc)
                        CYCLE
                     END IF
                  END IF

                  CALL abc_contract(iabc(sgfa:egfa, sgfb:egfb, sgfc:egfc), sabc, &
                                    sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                    ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))

                  DEALLOCATE (sabc)

               END DO !kset
            END DO !jset
         END DO !iset

         CALL set_o3c_container(o3c_iterator, mepos=mepos, integral=iabc, tvec=tvec)

      END DO !o3c_iterator
      CALL cp_libint_cleanup_3eri(lib)
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE calculate_o3c_libint_integrals

! **************************************************************************************************
!> \brief Contraction of 3-tensor over indices 1 and 2 (assuming symmetry)
!>        t(k) = sum_ij (ijk)*p(ij)
!> \param o3c ...
!> \param matrix_p ...
! **************************************************************************************************
   SUBROUTINE contract12_o3c(o3c, matrix_p)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract12_o3c', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, icol, ik, irow, ispin, &
                                                            jatom, mepos, nk, nspin, nthread
      LOGICAL                                            :: found, ijsymmetric, trans
      REAL(KIND=dp)                                      :: fpre
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, tvec
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: iabc
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      nspin = SIZE(matrix_p, 1)
      CALL get_o3c_container(o3c, ijsymmetric=ijsymmetric)
      CPASSERT(ijsymmetric)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,matrix_p,nspin)&
!$OMP PRIVATE (mepos,ispin,iatom,jatom,ik,nk,irow,icol,iabc,tvec,found,pblock,trans,fpre)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, iatom=iatom, jatom=jatom, &
                                    integral=iabc, tvec=tvec)
         nk = SIZE(tvec, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF
         IF (iatom == jatom) THEN
            fpre = 1.0_dp
         ELSE
            fpre = 2.0_dp
         END IF

         DO ispin = 1, nspin
            CALL dbcsr_get_block_p(matrix=matrix_p(ispin)%matrix, &
                                   row=irow, col=icol, BLOCK=pblock, found=found)
            IF (found) THEN
               IF (trans) THEN
                  DO ik = 1, nk
                     tvec(ik, ispin) = fpre*SUM(TRANSPOSE(pblock(:, :))*iabc(:, :, ik))
                  END DO
               ELSE
                  DO ik = 1, nk
                     tvec(ik, ispin) = fpre*SUM(pblock(:, :)*iabc(:, :, ik))
                  END DO
               END IF
            END IF
         END DO

      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE contract12_o3c

! **************************************************************************************************
!> \brief Contraction of 3-tensor over index 3
!>        h(ij) = h(ij) + sum_k (ijk)*v(k)
!> \param o3c ...
!> \param vec ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE contract3_o3c(o3c, vec, matrix)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: vec
      TYPE(dbcsr_type), POINTER                          :: matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract3_o3c', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, icol, ik, irow, jatom, &
                                                            katom, mepos, nk, nthread, s1, s2
      LOGICAL                                            :: found, ijsymmetric, trans
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: iabc
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      CALL get_o3c_container(o3c, ijsymmetric=ijsymmetric)
      CPASSERT(ijsymmetric)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,vec,matrix)&
!$OMP PRIVATE (mepos,iabc,iatom,jatom,katom,irow,icol,trans,pblock,v,found,ik,nk,work,s1,s2)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, iatom=iatom, jatom=jatom, katom=katom, &
                                    integral=iabc)

         CALL get_o3c_vec(vec, katom, v)
         nk = SIZE(v)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF

         CALL dbcsr_get_block_p(matrix=matrix, row=irow, col=icol, BLOCK=pblock, found=found)

         IF (found) THEN
            s1 = SIZE(pblock, 1); s2 = SIZE(pblock, 2)
            ALLOCATE (work(s1, s2))
            work(:, :) = 0.0_dp

            IF (trans) THEN
               DO ik = 1, nk
                  CALL daxpy(s1*s2, v(ik), TRANSPOSE(iabc(:, :, ik)), 1, work(:, :), 1)
               END DO
            ELSE
               DO ik = 1, nk
                  CALL daxpy(s1*s2, v(ik), iabc(:, :, ik), 1, work(:, :), 1)
               END DO
            END IF

            ! Mulitple threads with same irow, icol but different katom (same even in PBCs) can try
            ! to access the dbcsr block at the same time. Prevent that by CRITICAL section but keep
            ! computations before hand in order to retain speed

!$OMP CRITICAL
            CALL dbcsr_get_block_p(matrix=matrix, row=irow, col=icol, BLOCK=pblock, found=found)
            CALL daxpy(s1*s2, 1.0_dp, work(:, :), 1, pblock(:, :), 1)
!$OMP END CRITICAL

            DEALLOCATE (work)
         END IF

      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE contract3_o3c

END MODULE qs_o3c_methods
