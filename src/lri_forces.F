!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates forces for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created Dorothea Golze [03.2014]
!> \authors Dorothea Golze
! **************************************************************************************************
MODULE lri_forces

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE lri_environment_methods,         ONLY: lri_int_os,&
                                              lri_int_shg
   USE lri_environment_types,           ONLY: &
        allocate_lri_force_components, deallocate_lri_force_components, lri_density_type, &
        lri_environment_type, lri_force_type, lri_int_type, lri_kind_type, lri_list_type, &
        lri_rhoab_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_o3c_types,                    ONLY: get_o3c_iterator_info,&
                                              o3c_iterate,&
                                              o3c_iterator_create,&
                                              o3c_iterator_release,&
                                              o3c_iterator_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_forces'

   PUBLIC :: calculate_lri_forces, calculate_ri_forces

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief calculates the lri forces
!> \param lri_env ...
!> \param lri_density ...
!> \param qs_env ...
!> \param pmatrix density matrix
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE calculate_lri_forces(lri_env, lri_density, qs_env, pmatrix, atomic_kind_set)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_forces', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, ispin, natom, &
                                                            nkind, nspin
      LOGICAL                                            :: use_virial
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v_dadr, v_dfdr
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (atomic_kind, force, lri_coef, v_dadr, v_dfdr, virial)

      IF (ASSOCIATED(lri_env%soo_list)) THEN

         nkind = lri_env%lri_ints%nkind
         nspin = SIZE(pmatrix)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         CALL get_qs_env(qs_env=qs_env, force=force, virial=virial)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         !calculate SUM_i integral(V*fbas_i)*davec/dR
         CALL calculate_v_dadr(lri_env, lri_density, pmatrix, atomic_kind_set, &
                               use_virial, virial)

         DO ispin = 1, nspin

            lri_coef => lri_density%lri_coefs(ispin)%lri_kinds

            DO ikind = 1, nkind
               atomic_kind => atomic_kind_set(ikind)
               CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
               DO iatom = 1, natom
                  v_dadr => lri_coef(ikind)%v_dadr(iatom, :)
                  v_dfdr => lri_coef(ikind)%v_dfdr(iatom, :)

                  force(ikind)%rho_lri_elec(:, iatom) = force(ikind)%rho_lri_elec(:, iatom) &
                                                        +v_dfdr(:)+v_dadr(:)

               END DO
            END DO
         ENDDO

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE calculate_lri_forces

! **************************************************************************************************
!> \brief calculates second term of derivative with respect to R, i.e.
!>        SUM_i integral(V * fbas_i)*davec/dR
!> \param lri_env ...
!> \param lri_density ...
!> \param pmatrix ...
!> \param atomic_kind_set ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_v_dadr(lri_env, lri_density, pmatrix, atomic_kind_set, &
                               use_virial, virial)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_v_dadr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, atom_b, handle, i, iac, iatom, &
                                                            ikind, ilist, ispin, jatom, jkind, &
                                                            jneighbor, k, natom, nba, nbb, nfa, &
                                                            nfb, nkind, nlist, nn, nspin
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      LOGICAL                                            :: found, trans
      REAL(KIND=dp), DIMENSION(3)                        :: dcharge, dlambda, force_a, force_b, &
                                                            nsdssn, nsdsst, nsdt, rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: st, v_dadra, v_dadrb, v_inta, v_intb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: davec, dssn, dsst, dtvec, pbij, sdssn, &
                                                            sdsst, sdt
      TYPE(dbcsr_type), POINTER                          :: pmat
      TYPE(gto_basis_set_type), POINTER                  :: fbasa, fbasb, obasa, obasb
      TYPE(lri_force_type), POINTER                      :: lri_force
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(lri_list_type), POINTER                       :: lri_rho
      TYPE(lri_rhoab_type), POINTER                      :: lrho
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      CALL timeset(routineN, handle)
      NULLIFY (atom_of_kind, lri_coef, lri_force, lrii, lri_rho, lrho, &
               nl_iterator, pbij, pmat, soo_list, v_dadra, v_dadrb, v_inta, &
               v_intb)
      NULLIFY (davec, dssn, dsst, dtvec, sdssn, sdsst, sdt, st)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nkind = lri_env%lri_ints%nkind
         nspin = SIZE(pmatrix)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)

         DO ispin = 1, nspin

            pmat => pmatrix(ispin)%matrix
            lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
            lri_rho => lri_density%lri_rhos(ispin)%lri_list

            CALL neighbor_list_iterator_create(nl_iterator, soo_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, &
                                      nlist=nlist, ilist=ilist, inode=jneighbor, r=rab)

               iac = ikind+nkind*(jkind-1)

               IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

               lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
               lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)

               ! zero contribution for pairs aa or bb and periodc pairs aa' or bb'
               ! calculate forces for periodic pairs aa' and bb' only for virial
               IF (.NOT. lrii%calc_force_pair) CYCLE

               nfa = lrii%nfa
               nfb = lrii%nfb
               nba = lrii%nba
               nbb = lrii%nbb
               nn = nfa+nfb

               CALL allocate_lri_force_components(lri_density%lri_force, nfa, nfb)
               lri_force => lri_density%lri_force
               st => lri_force%st
               dsst => lri_force%dsst
               sdsst => lri_force%sdsst
               dssn => lri_force%dssn
               sdssn => lri_force%sdssn
               sdt => lri_force%sdt
               dtvec => lri_force%dtvec
               davec => lri_force%davec

               ! get the density matrix Pab
               NULLIFY (pbij)
               IF (iatom <= jatom) THEN
                  CALL dbcsr_get_block_p(matrix=pmat, row=iatom, col=jatom, block=pbij, found=found)
                  trans = .FALSE.
               ELSE
                  CALL dbcsr_get_block_p(matrix=pmat, row=jatom, col=iatom, block=pbij, found=found)
                  trans = .TRUE.
               END IF

               ! Calculate derivatives of overlap integrals (a,b), (fa,fb), (a,fa,b), (a,b,fb)
               IF (ASSOCIATED(lrii%dsab)) DEALLOCATE (lrii%dsab)
               IF (ASSOCIATED(lrii%dsoo)) DEALLOCATE (lrii%dsoo)
               IF (ASSOCIATED(lrii%dabbint)) DEALLOCATE (lrii%dabbint)
               IF (ASSOCIATED(lrii%dabdaint)) DEALLOCATE (lrii%dabdaint)
               ALLOCATE (lrii%dsab(nfa, nfb, 3))
               lrii%dsab = 0._dp
               ALLOCATE (lrii%dsoo(nba, nbb, 3))
               lrii%dsoo = 0._dp
               ALLOCATE (lrii%dabbint(nba, nbb, nfb, 3))
               lrii%dabbint = 0._dp
               ALLOCATE (lrii%dabdaint(nba, nbb, nfa, 3))
               lrii%dabdaint = 0._dp
               ! calculate integrals
               obasa => lri_env%orb_basis(ikind)%gto_basis_set
               obasb => lri_env%orb_basis(jkind)%gto_basis_set
               fbasa => lri_env%ri_basis(ikind)%gto_basis_set
               fbasb => lri_env%ri_basis(jkind)%gto_basis_set
               IF (lri_env%use_shg_integrals) THEN
                  CALL lri_int_shg(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                                   iatom, jatom, ikind, jkind, .TRUE.)
               ELSE
                  CALL lri_int_os(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                                  iatom, jatom, ikind, .TRUE.)
               ENDIF

               ! compute dtvec/dRa = SUM_ab Pab *d(a,b,x)/dRa
               DO k = 1, 3
                  IF (trans) THEN
                     dcharge(k) = SUM(TRANSPOSE(pbij(1:nbb, 1:nba))*lrii%dsoo(1:nba, 1:nbb, k))
                     DO i = 1, nfa
                        dtvec(i, k) = SUM(TRANSPOSE(pbij(1:nbb, 1:nba))*lrii%dabdaint(1:nba, 1:nbb, i, k))
                     END DO
                     DO i = 1, nfb
                        dtvec(nfa+i, k) = SUM(TRANSPOSE(pbij(1:nbb, 1:nba))*lrii%dabbint(1:nba, 1:nbb, i, k))
                     END DO
                  ELSE
                     dcharge(k) = SUM(pbij(1:nba, 1:nbb)*lrii%dsoo(1:nba, 1:nbb, k))
                     DO i = 1, nfa
                        dtvec(i, k) = SUM(pbij(1:nba, 1:nbb)*lrii%dabdaint(1:nba, 1:nbb, i, k))
                     END DO
                     DO i = 1, nfb
                        dtvec(nfa+i, k) = SUM(pbij(1:nba, 1:nbb)*lrii%dabbint(1:nba, 1:nbb, i, k))
                     END DO
                  END IF
               ENDDO

               ! derivative of S (overlap) matrix dS
               !dS: dsaa and dsbb are zero, only work with ab blocks in following

               ! n*Sinv*dS*Sinv*t
               CALL dgemm("N", "N", nn, 1, nn, 1.0_dp, lrii%sinv, nn, lrho%tvec, nn, &
                          0.0_dp, st, nn)
               DO k = 1, 3
                  CALL dgemm("N", "N", nfa, 1, nfb, 1.0_dp, lrii%dsab(1:nfa, 1:nfb, k), nfa, &
                             st(nfa+1:nn), nn, 0.0_dp, dsst(1:nfa, k), nn)
               ENDDO
               CALL dgemm("T", "N", 3*nfb, 1, nfa, 1.0_dp, lrii%dsab, nfa, st(1:nfa), nn, &
                          0.0_dp, dsst(nfa+1:nn, 1:3), nn*3)
               CALL dgemm("N", "N", nn, 3, nn, 1.0_dp, lrii%sinv, nn, dsst, nn, &
                          0.0_dp, sdsst, nn)

               DO k = 1, 3
                  nsdsst(k) = SUM(sdsst(1:nn, k)*lrii%n(1:nn))
               ENDDO

               ! n*Sinv*dS*Sinv*n
               DO k = 1, 3
                  CALL dgemm("N", "N", nfa, 1, nfb, 1.0_dp, lrii%dsab(1:nfa, 1:nfb, k), nfa, &
                             lrii%sn(nfa+1:nn), nn, 0.0_dp, dssn(1:nfa, k), nn)
               ENDDO
               CALL dgemm("T", "N", 3*nfb, 1, nfa, 1.0_dp, lrii%dsab, nfa, lrii%sn(1:nfa), nn, &
                          0.0_dp, dssn(nfa+1:nn, 1:3), nn*3)
               CALL dgemm("N", "N", nn, 3, nn, 1.0_dp, lrii%sinv, nn, dssn, nn, &
                          0.0_dp, sdssn, nn)
               DO k = 1, 3
                  nsdssn(k) = SUM(sdssn(1:nn, k)*lrii%n(1:nn))
               ENDDO

               ! Sinv*dtvec/dR_a
               CALL dgemm("N", "N", nn, 3, nn, 1.0_dp, lrii%sinv, nn, dtvec, nn, &
                          0.0_dp, sdt, nn)
               DO k = 1, 3
                  nsdt(k) = SUM(sdt(1:nn, k)*lrii%n(1:nn))
               ENDDO

               ! dlambda/dRa
               DO k = 1, 3
                  dlambda(k) = (dcharge(k)+nsdsst(k)-nsdt(k))/lrii%nsn &
                               +(lrho%charge-lrho%nst)*nsdssn(k)/(lrii%nsn*lrii%nsn)
               ENDDO

               ! davec/dRa, contains: davec_a/dRa (1..nfa) and davec_b/dRa (nfa+1..nn)
               DO k = 1, 3
                  davec(1:nn, k) = -sdsst(1:nn, k)+sdt(1:nn, k) &
                                   +dlambda(k)*lrii%sn(1:nn)-lrho%lambda*sdssn(1:nn, k)
               ENDDO

               ! deallocate derivative integrals
               IF (ASSOCIATED(lrii%dsab)) DEALLOCATE (lrii%dsab)
               IF (ASSOCIATED(lrii%dsoo)) DEALLOCATE (lrii%dsoo)
               IF (ASSOCIATED(lrii%dabbint)) DEALLOCATE (lrii%dabbint)
               IF (ASSOCIATED(lrii%dabdaint)) DEALLOCATE (lrii%dabdaint)
               NULLIFY (lrii%dsab, lrii%dsoo, lrii%dabbint, lrii%dabdaint)
               WRITE (6, *) "lrii%dabdaint   ", ASSOCIATED(lrii%dabdaint)

               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)

               v_dadra => lri_coef(ikind)%v_dadr(atom_a, :)
               v_dadrb => lri_coef(jkind)%v_dadr(atom_b, :)
               v_inta => lri_coef(ikind)%v_int(atom_a, :)
               v_intb => lri_coef(jkind)%v_int(atom_b, :)

               ! sum over atom pairs
               DO k = 1, 3
                  force_a(k) = 2.0_dp*SUM(davec(1:nfa, k)*v_inta(1:nfa)) &
                               +2.0_dp*SUM(davec(nfa+1:nn, k)*v_intb(1:nfb))
                  v_dadra(k) = v_dadra(k)+force_a(k)
                  !minus: using that davec/dRa = - davec/dRb
                  force_b(k) = -2.0_dp*SUM(davec(nfa+1:nn, k)*v_intb(1:nfb)) &
                               -2.0_dp*SUM(davec(1:nfa, k)*v_inta(1:nfa))
                  v_dadrb(k) = v_dadrb(k)+force_b(k)
               ENDDO

               ! contribution to virial
               IF (use_virial) THEN
                  !periodic self-pairs aa' contribute only with factor 0.5
                  IF (iatom == jatom) THEN
                     CALL virial_pair_force(virial%pv_virial, 0.5_dp, force_a, rab)
                  ELSE
                     CALL virial_pair_force(virial%pv_virial, 1.0_dp, force_a, rab)
                  ENDIF
               ENDIF

               CALL deallocate_lri_force_components(lri_density%lri_force)

            END DO

            CALL neighbor_list_iterator_release(nl_iterator)

         ENDDO

         DEALLOCATE (atom_of_kind)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE calculate_v_dadr

! **************************************************************************************************
!> \brief calculates the ri forces
!> \param lri_env ...
!> \param lri_density ...
!> \param qs_env ...
!> \param pmatrix density matrix
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE calculate_ri_forces(lri_env, lri_density, qs_env, pmatrix, atomic_kind_set)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_ri_forces', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, ispin, natom, &
                                                            nkind, nspin
      LOGICAL                                            :: use_virial
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v_dadr, v_dfdr
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (atomic_kind, force, lri_coef, v_dadr, v_dfdr, virial)

      nkind = SIZE(atomic_kind_set)
      nspin = SIZE(pmatrix)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
      CALL get_qs_env(qs_env=qs_env, force=force, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      !calculate SUM_i integral(V*fbas_i)*davec/dR
      CALL calculate_v_dadr_ri(lri_env, lri_density, pmatrix, atomic_kind_set, &
                               use_virial, virial)

      DO ispin = 1, nspin

         lri_coef => lri_density%lri_coefs(ispin)%lri_kinds

         DO ikind = 1, nkind
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
            DO iatom = 1, natom
               v_dadr => lri_coef(ikind)%v_dadr(iatom, :)
               v_dfdr => lri_coef(ikind)%v_dfdr(iatom, :)

               force(ikind)%rho_lri_elec(:, iatom) = force(ikind)%rho_lri_elec(:, iatom) &
                                                     +v_dfdr(:)+v_dadr(:)

            END DO
         END DO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE calculate_ri_forces

! **************************************************************************************************
!> \brief calculates second term of derivative with respect to R, i.e.
!>        SUM_i integral(V * fbas_i)*davec/dR
!>        However contributions from charge constraint and derivative of overlap matrix have been
!>        calculated in calculate_ri_integrals
!> \param lri_env ...
!> \param lri_density ...
!> \param pmatrix ...
!> \param atomic_kind_set ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_v_dadr_ri(lri_env, lri_density, pmatrix, atomic_kind_set, &
                                  use_virial, virial)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pmatrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_v_dadr_ri', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, atom_b, atom_c, handle, i, i1, &
                                                            i2, iatom, ikind, ispin, jatom, jkind, &
                                                            katom, kkind, m, mepos, natom, nkind, &
                                                            nspin, nthread
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      INTEGER, DIMENSION(:, :), POINTER                  :: bas_ptr
      REAL(KIND=dp)                                      :: fscal
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, force_c, rij, rik, rjk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v_dadra, v_dadrb, v_dadrc
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fi, fj, fk, fo
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      nspin = SIZE(pmatrix)
      nkind = SIZE(atomic_kind_set)
      DO ispin = 1, nspin
         lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
         DO ikind = 1, nkind
            lri_coef(ikind)%v_dadr(:, :) = 0.0_dp
         END DO
      END DO

      bas_ptr => lri_env%ri_fit%bas_ptr

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
      ALLOCATE (atom_of_kind(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      ! contribution from 3-center overlap integral derivatives
      fo => lri_env%ri_fit%fout
      nthread = 1
!$    nthread = omp_get_max_threads()

      CALL o3c_iterator_create(lri_env%o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,o3c_iterator,bas_ptr,fo,nspin,atom_of_kind,lri_coef,virial,use_virial)&
!$OMP PRIVATE (mepos,ikind,jkind,kkind,iatom,jatom,katom,fi,fj,fk,fscal,&
!$OMP          force_a,force_b,force_c,i1,i2,m,i,ispin,atom_a,atom_b,atom_c,&
!$OMP          v_dadra,v_dadrb,v_dadrc,rij,rik,rjk)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, &
                                    ikind=ikind, jkind=jkind, kkind=kkind, &
                                    iatom=iatom, jatom=jatom, katom=katom, &
                                    rij=rij, rik=rik, force_i=fi, force_j=fj, force_k=fk)
         i1 = bas_ptr(1, katom)
         i2 = bas_ptr(2, katom)
         m = i2-i1+1
         DO i = 1, 3
            force_a(i) = 0.0_dp
            force_b(i) = 0.0_dp
            force_c(i) = 0.0_dp
            DO ispin = 1, nspin
               force_a(i) = force_a(i)+SUM(fi(1:m, i)*fo(i1:i2, ispin))
               force_b(i) = force_b(i)+SUM(fj(1:m, i)*fo(i1:i2, ispin))
               force_c(i) = force_c(i)+SUM(fk(1:m, i)*fo(i1:i2, ispin))
            END DO
         END DO
         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)
         atom_c = atom_of_kind(katom)
         !
         v_dadra => lri_coef(ikind)%v_dadr(atom_a, :)
         v_dadrb => lri_coef(jkind)%v_dadr(atom_b, :)
         v_dadrc => lri_coef(kkind)%v_dadr(atom_c, :)
         !
!$OMP CRITICAL(addforce)
         v_dadra(1:3) = v_dadra(1:3)+force_a(1:3)
         v_dadrb(1:3) = v_dadrb(1:3)+force_b(1:3)
         v_dadrc(1:3) = v_dadrc(1:3)+force_c(1:3)
         !
         IF (use_virial) THEN
            rjk(1:3) = rik(1:3)-rij(1:3)
            ! to be debugged
            fscal = 1.0_dp
            IF (iatom == jatom) fscal = 1.0_dp
            CALL virial_pair_force(virial%pv_virial, 1.0_dp, force_a, rik)
            CALL virial_pair_force(virial%pv_virial, 1.0_dp, force_b, rjk)
         ENDIF
!$OMP END CRITICAL(addforce)
      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      DEALLOCATE (atom_of_kind)

      CALL timestop(handle)

   END SUBROUTINE calculate_v_dadr_ri

END MODULE lri_forces
