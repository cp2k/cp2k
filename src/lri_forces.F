!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates forces for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created Dorothea Golze [03.2014]
!>      refactoring and distant pair approximation JGH [08.2017]
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
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lri_environment_types,           ONLY: &
        allocate_lri_force_components, deallocate_lri_force_components, lri_density_type, &
        lri_environment_type, lri_force_type, lri_int_type, lri_kind_type, lri_list_type, &
        lri_rhoab_type
   USE lri_integrals,                   ONLY: allocate_int_type,&
                                              deallocate_int_type,&
                                              dint_type,&
                                              lri_dint2
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
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
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: pmatrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_forces', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, ispin, natom, &
                                                            nkind, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_virial
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v_dadr, v_dfdr
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (atomic_kind, force, lri_coef, v_dadr, v_dfdr, virial)

      IF (ASSOCIATED(lri_env%soo_list)) THEN

         nkind = lri_env%lri_ints%nkind
         nspin = SIZE(pmatrix, 1)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         CALL get_qs_env(qs_env=qs_env, force=force, virial=virial)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

         !calculate SUM_i integral(V*fbas_i)*davec/dR
         CALL calculate_v_dadr_sr(lri_env, lri_density, pmatrix, cell_to_index, atomic_kind_set, &
                                  use_virial, virial)
         IF (lri_env%distant_pair_approximation) THEN
            CALL calculate_v_dadr_ff(lri_env, lri_density, pmatrix, cell_to_index, atomic_kind_set, &
                                     use_virial, virial)
         END IF

         DO ispin = 1, nspin

            lri_coef => lri_density%lri_coefs(ispin)%lri_kinds

            DO ikind = 1, nkind
               atomic_kind => atomic_kind_set(ikind)
               CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
               DO iatom = 1, natom
                  v_dadr => lri_coef(ikind)%v_dadr(iatom, :)
                  v_dfdr => lri_coef(ikind)%v_dfdr(iatom, :)

                  force(ikind)%rho_lri_elec(:, iatom) = force(ikind)%rho_lri_elec(:, iatom) &
                                                        + v_dfdr(:) + v_dadr(:)

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
!> \param cell_to_index ...
!> \param atomic_kind_set ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_v_dadr_sr(lri_env, lri_density, pmatrix, cell_to_index, atomic_kind_set, &
                                  use_virial, virial)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: pmatrix
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_v_dadr_sr', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, atom_b, handle, i, iac, iatom, ic, ikind, ilist, ispin, jatom, jkind, &
         jneighbor, k, mepos, natom, nba, nbb, nfa, nfb, nkind, nlist, nn, nspin, nthread
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      LOGICAL                                            :: found, use_cell_mapping
      REAL(KIND=dp)                                      :: ai, dab, dfw, fw, isn, threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pab
      REAL(KIND=dp), DIMENSION(3)                        :: dcharge, dlambda, force_a, force_b, &
                                                            idav, nsdssn, nsdsst, nsdt, rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: st, v_dadra, v_dadrb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dssn, dsst, dtvec, pbij, sdssn, sdsst, &
                                                            sdt
      TYPE(dbcsr_type), POINTER                          :: pmat
      TYPE(dint_type)                                    :: lridint
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
               nl_iterator, pbij, pmat, soo_list, v_dadra, v_dadrb)
      NULLIFY (dssn, dsst, dtvec, sdssn, sdsst, sdt, st)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nkind = lri_env%lri_ints%nkind
         nspin = SIZE(pmatrix, 1)
         nthread = 1
!$       nthread = omp_get_max_threads()
         use_cell_mapping = (SIZE(pmatrix, 2) > 1)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)

         DO ispin = 1, nspin

            lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
            lri_rho => lri_density%lri_rhos(ispin)%lri_list

            CALL neighbor_list_iterator_create(nl_iterator, soo_list, nthread=nthread)
!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,nl_iterator,pmatrix,nkind,lri_env,lri_rho,lri_coef,&
!$OMP         atom_of_kind,virial,use_virial,cell_to_index,use_cell_mapping,ispin)&
!$OMP PRIVATE (mepos,ikind,jkind,iatom,jatom,nlist,ilist,jneighbor,rab,iac,lrii,lrho,&
!$OMP          lri_force,nfa,nfb,nba,nbb,nn,st,nsdsst,nsdssn,dsst,sdsst,dssn,sdssn,sdt,&
!$OMP          nsdt,dtvec,pbij,found,obasa,obasb,fbasa,fbasb,i,k,threshold,&
!$OMP          dcharge,dlambda,atom_a,atom_b,v_dadra,v_dadrb,force_a,force_b,&
!$OMP          lridint,pab,fw,dfw,dab,ai,vint,isn,idav,cell,ic,pmat)

            mepos = 0
!$          mepos = omp_get_thread_num()

            DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
               CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, nlist=nlist, ilist=ilist, &
                                      inode=jneighbor, r=rab, cell=cell)

               iac = ikind + nkind*(jkind - 1)

               IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

               lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
               lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)

               ! only proceed if short-range exact LRI is needed
               IF (.NOT. lrii%lrisr) CYCLE
               fw = lrii%wsr
               dfw = lrii%dwsr

               ! zero contribution for pairs aa or bb and periodc pairs aa' or bb'
               ! calculate forces for periodic pairs aa' and bb' only for virial
               IF (.NOT. lrii%calc_force_pair) CYCLE

               nfa = lrii%nfa
               nfb = lrii%nfb
               nba = lrii%nba
               nbb = lrii%nbb
               nn = nfa + nfb

               IF (use_cell_mapping) THEN
                  ic = cell_to_index(cell(1), cell(2), cell(3))
                  CPASSERT(ic > 0)
               ELSE
                  ic = 1
               END IF
               pmat => pmatrix(ispin, ic)%matrix

               ! get the density matrix Pab
               ALLOCATE (pab(nba, nbb))
               NULLIFY (pbij)
               IF (iatom <= jatom) THEN
                  CALL dbcsr_get_block_p(matrix=pmat, row=iatom, col=jatom, block=pbij, found=found)
                  CPASSERT(found)
                  pab(1:nba, 1:nbb) = pbij(1:nba, 1:nbb)
               ELSE
                  CALL dbcsr_get_block_p(matrix=pmat, row=jatom, col=iatom, block=pbij, found=found)
                  CPASSERT(found)
                  pab(1:nba, 1:nbb) = TRANSPOSE(pbij(1:nbb, 1:nba))
               END IF

               obasa => lri_env%orb_basis(ikind)%gto_basis_set
               obasb => lri_env%orb_basis(jkind)%gto_basis_set
               fbasa => lri_env%ri_basis(ikind)%gto_basis_set
               fbasb => lri_env%ri_basis(jkind)%gto_basis_set

               ! Calculate derivatives of overlap integrals (a,b), (fa,fb), (a,fa,b), (a,b,fb)
               CALL allocate_int_type(lridint=lridint, nba=nba, nbb=nbb, nfa=nfa, nfb=nfb)
               CALL lri_dint2(lri_env, lrii, lridint, rab, obasa, obasb, fbasa, fbasb, &
                              iatom, jatom, ikind, jkind)

               NULLIFY (lri_force)
               CALL allocate_lri_force_components(lri_force, nfa, nfb)
               st => lri_force%st
               dsst => lri_force%dsst
               dssn => lri_force%dssn
               dtvec => lri_force%dtvec

               ! compute dtvec/dRa = SUM_ab Pab *d(a,b,x)/dRa
               DO k = 1, 3
                  threshold = 0.01_dp*lri_env%eps_o3_int/MAX(SUM(ABS(pab(1:nba, 1:nbb))), 1.0e-14_dp)
                  dcharge(k) = SUM(pab(1:nba, 1:nbb)*lridint%dsooint(1:nba, 1:nbb, k))
                  DO i = 1, nfa
                     IF (lrii%abascr(i) > threshold) THEN
                        dtvec(i, k) = SUM(pab(1:nba, 1:nbb)*lridint%dabaint(1:nba, 1:nbb, i, k))
                     END IF
                  END DO
                  DO i = 1, nfb
                     IF (lrii%abbscr(i) > threshold) THEN
                        dtvec(nfa + i, k) = SUM(pab(1:nba, 1:nbb)*lridint%dabbint(1:nba, 1:nbb, i, k))
                     END IF
                  END DO
               ENDDO

               DEALLOCATE (pab)

               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               ALLOCATE (vint(nn))
               vint(1:nfa) = lri_coef(ikind)%v_int(atom_a, 1:nfa)
               vint(nfa + 1:nn) = lri_coef(jkind)%v_int(atom_b, 1:nfb)

               isn = SUM(vint(1:nn)*lrii%sn(1:nn))
               DO k = 1, 3
                  ai = isn/lrii%nsn*dcharge(k)
                  force_a(k) = 2.0_dp*fw*ai
                  force_b(k) = -2.0_dp*fw*ai
               ENDDO

               ! derivative of S (overlap) matrix dS
               !dS: dsaa and dsbb are zero, only work with ab blocks in following
               st(1:nn) = MATMUL(lrii%sinv(1:nn, 1:nn), lrho%tvec(1:nn))
               DO k = 1, 3
                  dsst(1:nfa, k) = MATMUL(lridint%dsabint(1:nfa, 1:nfb, k), st(nfa + 1:nn))
                  dsst(nfa + 1:nn, k) = MATMUL(st(1:nfa), lridint%dsabint(1:nfa, 1:nfb, k))
                  nsdsst(k) = SUM(lrii%sn(1:nn)*dsst(1:nn, k))
                  dssn(1:nfa, k) = MATMUL(lridint%dsabint(1:nfa, 1:nfb, k), lrii%sn(nfa + 1:nn))
                  dssn(nfa + 1:nn, k) = MATMUL(lrii%sn(1:nfa), lridint%dsabint(1:nfa, 1:nfb, k))
                  nsdssn(k) = SUM(lrii%sn(1:nn)*dssn(1:nn, k))
                  nsdt(k) = SUM(dtvec(1:nn, k)*lrii%sn(1:nn))
               ENDDO
               ! dlambda/dRa
               DO k = 1, 3
                  dlambda(k) = (nsdsst(k) - nsdt(k))/lrii%nsn &
                               + (lrho%charge - lrho%nst)*nsdssn(k)/(lrii%nsn*lrii%nsn)
               ENDDO
               DO k = 1, 3
                  force_a(k) = force_a(k) + 2.0_dp*fw*isn*dlambda(k)
                  force_b(k) = force_b(k) - 2.0_dp*fw*isn*dlambda(k)
               ENDDO
               DO k = 1, 3
                  st(1:nn) = dtvec(1:nn, k) - dsst(1:nn, k) - lrho%lambda*dssn(1:nn, k)
                  idav(k) = SUM(vint(1:nn)*MATMUL(lrii%sinv(1:nn, 1:nn), st(1:nn)))
               ENDDO

               ! deallocate derivative integrals
               CALL deallocate_int_type(lridint=lridint)

               ! sum over atom pairs
               DO k = 1, 3
                  ai = 2.0_dp*fw*idav(k)
                  force_a(k) = force_a(k) + ai
                  force_b(k) = force_b(k) - ai
               ENDDO
               IF (ABS(dfw) > 0.0_dp) THEN
                  dab = SQRT(SUM(rab(1:3)*rab(1:3)))
                  ai = 2.0_dp*dfw/dab*SUM(lrho%avec(1:nn)*vint(1:nn))
                  DO k = 1, 3
                     force_a(k) = force_a(k) - ai*rab(k)
                     force_b(k) = force_b(k) + ai*rab(k)
                  ENDDO
               END IF

               DEALLOCATE (vint)

               v_dadra => lri_coef(ikind)%v_dadr(atom_a, :)
               v_dadrb => lri_coef(jkind)%v_dadr(atom_b, :)
!$OMP CRITICAL(addforces)
               DO k = 1, 3
                  v_dadra(k) = v_dadra(k) + force_a(k)
                  v_dadrb(k) = v_dadrb(k) + force_b(k)
               ENDDO
!$OMP END CRITICAL(addforces)

               ! contribution to virial
               IF (use_virial) THEN
                  !periodic self-pairs aa' contribute only with factor 0.5
!$OMP CRITICAL(addvirial)
                  IF (iatom == jatom) THEN
                     CALL virial_pair_force(virial%pv_virial, 0.5_dp, force_a, rab)
                  ELSE
                     CALL virial_pair_force(virial%pv_virial, 1.0_dp, force_a, rab)
                  ENDIF
!$OMP END CRITICAL(addvirial)
               ENDIF

               CALL deallocate_lri_force_components(lri_force)

            END DO
!$OMP END PARALLEL

            CALL neighbor_list_iterator_release(nl_iterator)

         ENDDO

         DEALLOCATE (atom_of_kind)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE calculate_v_dadr_sr

! **************************************************************************************************
!> \brief calculates second term of derivative with respect to R, i.e.
!>        SUM_i integral(V * fbas_i)*davec/dR for distant pair approximation
!> \param lri_env ...
!> \param lri_density ...
!> \param pmatrix ...
!> \param cell_to_index ...
!> \param atomic_kind_set ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_v_dadr_ff(lri_env, lri_density, pmatrix, cell_to_index, atomic_kind_set, &
                                  use_virial, virial)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: pmatrix
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_v_dadr_ff', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, atom_b, handle, i, iac, iatom, ic, ikind, ilist, ispin, jatom, jkind, &
         jneighbor, k, mepos, natom, nba, nbb, nfa, nfb, nkind, nlist, nn, nspin, nthread
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      LOGICAL                                            :: found, use_cell_mapping
      REAL(KIND=dp)                                      :: ai, dab, dfw, fw, isna, isnb, threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pab, wab, wbb
      REAL(KIND=dp), DIMENSION(3)                        :: dchargea, dchargeb, force_a, force_b, &
                                                            idava, idavb, rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: sta, stb, v_dadra, v_dadrb, vinta, vintb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dtveca, dtvecb, pbij
      TYPE(dbcsr_type), POINTER                          :: pmat
      TYPE(dint_type)                                    :: lridint
      TYPE(gto_basis_set_type), POINTER                  :: fbasa, fbasb, obasa, obasb
      TYPE(lri_force_type), POINTER                      :: lri_force_a, lri_force_b
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(lri_list_type), POINTER                       :: lri_rho
      TYPE(lri_rhoab_type), POINTER                      :: lrho
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      CALL timeset(routineN, handle)

      NULLIFY (atom_of_kind, lri_coef, lrii, lri_rho, lrho, &
               nl_iterator, pbij, pmat, soo_list, v_dadra, v_dadrb)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nkind = lri_env%lri_ints%nkind
         nspin = SIZE(pmatrix, 1)
         nthread = 1
!$       nthread = omp_get_max_threads()
         use_cell_mapping = (SIZE(pmatrix, 2) > 1)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)

         DO ispin = 1, nspin

            lri_coef => lri_density%lri_coefs(ispin)%lri_kinds
            lri_rho => lri_density%lri_rhos(ispin)%lri_list

            CALL neighbor_list_iterator_create(nl_iterator, soo_list, nthread=nthread)
!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,nl_iterator,pmatrix,nkind,lri_env,lri_rho,lri_coef,&
!$OMP         atom_of_kind,virial,use_virial,cell_to_index,use_cell_mapping,ispin)&
!$OMP PRIVATE (mepos,ikind,jkind,iatom,jatom,nlist,ilist,jneighbor,rab,iac,lrii,lrho,&
!$OMP          lri_force_a,lri_force_b,nfa,nfb,nba,nbb,nn,sta,stb,&
!$OMP          dtveca,dtvecb,pbij,found,obasa,obasb,fbasa,fbasb,i,k,threshold,&
!$OMP          dchargea,dchargeb,atom_a,atom_b,v_dadra,v_dadrb,force_a,force_b,&
!$OMP          lridint,pab,wab,wbb,fw,dfw,dab,ai,vinta,vintb,isna,isnb,idava,idavb,cell,ic,pmat)

            mepos = 0
!$          mepos = omp_get_thread_num()

            DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
               CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, nlist=nlist, ilist=ilist, &
                                      inode=jneighbor, r=rab, cell=cell)

               iac = ikind + nkind*(jkind - 1)

               IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

               lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
               lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)

               ! only proceed if short-range exact LRI is needed
               IF (.NOT. lrii%lriff) CYCLE
               fw = lrii%wff
               dfw = lrii%dwff

               ! zero contribution for pairs aa or bb and periodc pairs aa' or bb'
               ! calculate forces for periodic pairs aa' and bb' only for virial
               IF (.NOT. lrii%calc_force_pair) CYCLE

               nfa = lrii%nfa
               nfb = lrii%nfb
               nba = lrii%nba
               nbb = lrii%nbb
               nn = nfa + nfb

               IF (use_cell_mapping) THEN
                  ic = cell_to_index(cell(1), cell(2), cell(3))
                  CPASSERT(ic > 0)
               ELSE
                  ic = 1
               END IF
               pmat => pmatrix(ispin, ic)%matrix

               ! get the density matrix Pab
               ALLOCATE (pab(nba, nbb))
               NULLIFY (pbij)
               IF (iatom <= jatom) THEN
                  CALL dbcsr_get_block_p(matrix=pmat, row=iatom, col=jatom, block=pbij, found=found)
                  CPASSERT(found)
                  pab(1:nba, 1:nbb) = pbij(1:nba, 1:nbb)
               ELSE
                  CALL dbcsr_get_block_p(matrix=pmat, row=jatom, col=iatom, block=pbij, found=found)
                  CPASSERT(found)
                  pab(1:nba, 1:nbb) = TRANSPOSE(pbij(1:nbb, 1:nba))
               END IF

               ALLOCATE (wab(nba, nbb), wbb(nba, nbb))
               wab(1:nba, 1:nbb) = lri_env%wmat(ikind, jkind)%mat(1:nba, 1:nbb)
               wbb(1:nba, 1:nbb) = 1.0_dp - lri_env%wmat(ikind, jkind)%mat(1:nba, 1:nbb)

               obasa => lri_env%orb_basis(ikind)%gto_basis_set
               obasb => lri_env%orb_basis(jkind)%gto_basis_set
               fbasa => lri_env%ri_basis(ikind)%gto_basis_set
               fbasb => lri_env%ri_basis(jkind)%gto_basis_set

               ! Calculate derivatives of overlap integrals (a,b), (fa,fb), (a,fa,b), (a,b,fb)
               CALL allocate_int_type(lridint=lridint, nba=nba, nbb=nbb, nfa=nfa, nfb=nfb, &
                                      skip_dsab=.TRUE.)
               CALL lri_dint2(lri_env, lrii, lridint, rab, obasa, obasb, fbasa, fbasb, &
                              iatom, jatom, ikind, jkind)

               NULLIFY (lri_force_a, lri_force_b)
               CALL allocate_lri_force_components(lri_force_a, nfa, 0)
               CALL allocate_lri_force_components(lri_force_b, 0, nfb)
               dtveca => lri_force_a%dtvec
               dtvecb => lri_force_b%dtvec
               sta => lri_force_a%st
               stb => lri_force_b%st

               ! compute dtvec/dRa = SUM_ab Pab *d(a,b,x)/dRa
               threshold = 0.01_dp*lri_env%eps_o3_int/MAX(SUM(ABS(pab(1:nba, 1:nbb))), 1.0e-14_dp)
               DO k = 1, 3
                  dchargea(k) = SUM(wab(1:nba, 1:nbb)*pab(1:nba, 1:nbb)*lridint%dsooint(1:nba, 1:nbb, k))
                  DO i = 1, nfa
                     IF (lrii%abascr(i) > threshold) THEN
                        dtveca(i, k) = SUM(wab(1:nba, 1:nbb)*pab(1:nba, 1:nbb)*lridint%dabaint(1:nba, 1:nbb, i, k))
                     END IF
                  END DO
                  dchargeb(k) = SUM(wbb(1:nba, 1:nbb)*pab(1:nba, 1:nbb)*lridint%dsooint(1:nba, 1:nbb, k))
                  DO i = 1, nfb
                     IF (lrii%abbscr(i) > threshold) THEN
                        dtvecb(i, k) = SUM(wbb(1:nba, 1:nbb)*pab(1:nba, 1:nbb)*lridint%dabbint(1:nba, 1:nbb, i, k))
                     END IF
                  END DO
               ENDDO

               DEALLOCATE (pab, wab, wbb)

               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               vinta => lri_coef(ikind)%v_int(atom_a, :)
               vintb => lri_coef(jkind)%v_int(atom_b, :)

               isna = SUM(vinta(1:nfa)*lrii%sna(1:nfa))
               isnb = SUM(vintb(1:nfb)*lrii%snb(1:nfb))
               DO k = 1, 3
                  ai = isna/lrii%nsna*dchargea(k) + isnb/lrii%nsnb*dchargeb(k)
                  force_a(k) = 2.0_dp*fw*ai
                  force_b(k) = -2.0_dp*fw*ai
               ENDDO

               DO k = 1, 3
                  sta(1:nfa) = MATMUL(lrii%asinv(1:nfa, 1:nfa), dtveca(1:nfa, k))
                  idava(k) = SUM(vinta(1:nfa)*sta(1:nfa)) - isna/lrii%nsna*SUM(lrii%na(1:nfa)*sta(1:nfa))
                  stb(1:nfb) = MATMUL(lrii%bsinv(1:nfb, 1:nfb), dtvecb(1:nfb, k))
                  idavb(k) = SUM(vintb(1:nfb)*stb(1:nfb)) - isnb/lrii%nsnb*SUM(lrii%nb(1:nfb)*stb(1:nfb))
               ENDDO

               ! deallocate derivative integrals
               CALL deallocate_int_type(lridint=lridint)

               ! sum over atom pairs
               DO k = 1, 3
                  ai = 2.0_dp*fw*(idava(k) + idavb(k))
                  force_a(k) = force_a(k) + ai
                  force_b(k) = force_b(k) - ai
               ENDDO
               IF (ABS(dfw) > 0.0_dp) THEN
                  dab = SQRT(SUM(rab(1:3)*rab(1:3)))
                  ai = 2.0_dp*dfw/dab* &
                       (SUM(lrho%aveca(1:nfa)*vinta(1:nfa)) + &
                        SUM(lrho%avecb(1:nfb)*vintb(1:nfb)))
                  DO k = 1, 3
                     force_a(k) = force_a(k) - ai*rab(k)
                     force_b(k) = force_b(k) + ai*rab(k)
                  ENDDO
               END IF
               v_dadra => lri_coef(ikind)%v_dadr(atom_a, :)
               v_dadrb => lri_coef(jkind)%v_dadr(atom_b, :)

!$OMP CRITICAL(addforces)
               DO k = 1, 3
                  v_dadra(k) = v_dadra(k) + force_a(k)
                  v_dadrb(k) = v_dadrb(k) + force_b(k)
               ENDDO
!$OMP END CRITICAL(addforces)

               ! contribution to virial
               IF (use_virial) THEN
                  !periodic self-pairs aa' contribute only with factor 0.5
!$OMP CRITICAL(addvirial)
                  IF (iatom == jatom) THEN
                     CALL virial_pair_force(virial%pv_virial, 0.5_dp, force_a, rab)
                  ELSE
                     CALL virial_pair_force(virial%pv_virial, 1.0_dp, force_a, rab)
                  ENDIF
!$OMP END CRITICAL(addvirial)
               ENDIF

               CALL deallocate_lri_force_components(lri_force_a)
               CALL deallocate_lri_force_components(lri_force_b)

            END DO
!$OMP END PARALLEL

            CALL neighbor_list_iterator_release(nl_iterator)

         ENDDO

         DEALLOCATE (atom_of_kind)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE calculate_v_dadr_ff

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
                                                     + v_dfdr(:) + v_dadr(:)

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
         m = i2 - i1 + 1
         DO i = 1, 3
            force_a(i) = 0.0_dp
            force_b(i) = 0.0_dp
            force_c(i) = 0.0_dp
            DO ispin = 1, nspin
               force_a(i) = force_a(i) + SUM(fi(1:m, i)*fo(i1:i2, ispin))
               force_b(i) = force_b(i) + SUM(fj(1:m, i)*fo(i1:i2, ispin))
               force_c(i) = force_c(i) + SUM(fk(1:m, i)*fo(i1:i2, ispin))
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
         v_dadra(1:3) = v_dadra(1:3) + force_a(1:3)
         v_dadrb(1:3) = v_dadrb(1:3) + force_b(1:3)
         v_dadrc(1:3) = v_dadrc(1:3) + force_c(1:3)
         !
         IF (use_virial) THEN
            rjk(1:3) = rik(1:3) - rij(1:3)
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
