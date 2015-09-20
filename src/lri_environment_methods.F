!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates integral matrices for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) heavily debugged
!> \authors JGH 
!>          Dorothea Golze  
! *****************************************************************************
MODULE lri_environment_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_lri_inv,&
                                             do_lri_inv_auto,&
                                             do_lri_pseudoinv_diag,&
                                             do_lri_pseudoinv_svd
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: &
       allocate_lri_coefs, allocate_lri_ints, allocate_lri_ints_rho, &
       allocate_lri_rhos, deallocate_lri_ints, deallocate_lri_ints_rho, &
       lri_density_create, lri_density_release, lri_density_type, &
       lri_environment_type, lri_int_rho_type, lri_int_type, lri_kind_type, &
       lri_list_type, lri_rhoab_type
  USE lri_os_integrals,                ONLY: lri_int_aabb,&
                                             lri_int_ab,&
                                             lri_int_aba,&
                                             lri_int_abb
  USE lri_shg_integrals,               ONLY: lri_int_ab_shg,&
                                             lri_int_aba_shg,&
                                             lri_int_abb_shg,&
                                             precalc_angular_shg_part
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
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_methods'

  PUBLIC :: build_lri_matrices, calculate_lri_densities, calculate_lri_integrals,&
            calculate_lri_overlap_aabb, calculate_avec

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!> \param qs_env ...
!> \param calculate_forces ...
! *****************************************************************************
  SUBROUTINE build_lri_matrices(lri_env,qs_env,calculate_forces)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_lri_matrices', &
      routineP = moduleN//':'//routineN

! calculate the integrals needed to do the local (2-center) expansion
! of the (pair) densities

    CALL calculate_lri_integrals(lri_env,qs_env,calculate_forces)

  END SUBROUTINE build_lri_matrices

! *****************************************************************************
!> \brief calculates integrals needed for the LRI density fitting,
!>        integrals are calculated once, before the SCF starts
!> \param lri_env ...
!> \param qs_env ...
!> \param calculate_forces ...
! *****************************************************************************
  SUBROUTINE calculate_lri_integrals(lri_env,qs_env,calculate_forces)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iac, iatom, ikind, ilist, jatom, jkind, jneighbor, &
      nba, nbb, nfa, nfb, nkind, nlist, nn, nneighbor
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: fbasa, fbasb, obasa, obasb
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(cell, dft_control, fbasa, fbasb, lrii, lri_ints, nl_iterator, &
         obasa, obasb, particle_set, soo_list, virial)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
                      nkind=nkind,particle_set=particle_set,virial=virial)

      IF ( ASSOCIATED(lri_env%lri_ints) ) THEN
        CALL deallocate_lri_ints (lri_env%lri_ints)
      END IF

      ! allocate matrices storing the LRI integrals
      CALL allocate_lri_ints(lri_env,lri_env%lri_ints,nkind,&
                             calculate_forces,virial)
      lri_ints => lri_env%lri_ints

      CALL neighbor_list_iterator_create(nl_iterator,soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator)==0)

         CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
              nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
              iatom=iatom,jatom=jatom,r=rab)
        
         iac = ikind + nkind*(jkind - 1)
         dab = SQRT(SUM(rab*rab))

         obasa => lri_env%orb_basis(ikind)%gto_basis_set
         obasb => lri_env%orb_basis(jkind)%gto_basis_set
         fbasa => lri_env%ri_basis(ikind)%gto_basis_set
         fbasb => lri_env%ri_basis(jkind)%gto_basis_set

         IF (.NOT.ASSOCIATED(obasa)) CYCLE
         IF (.NOT.ASSOCIATED(obasb)) CYCLE

         lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

         nba = obasa%nsgf
         nbb = obasb%nsgf
         nfa = fbasa%nsgf
         nfb = fbasb%nsgf

         lrii%nba = nba
         lrii%nbb = nbb
         lrii%nfa = nfa
         lrii%nfb = nfb

         ra(:) = pbc(particle_set(iatom)%r, cell) 
         rb(:) = pbc(particle_set(jatom)%r, cell)

         ! calculate integrals (fa,fb), (a,b), (a,b,fa) and (a,b,fb)
         IF(lri_env%use_shg_integrals) THEN
           CALL lri_int_shg(lri_env,lrii,ra,rb,rab,obasa,obasb,fbasa,fbasb,&
                            iatom,jatom,ikind,jkind)
         ELSE
           CALL lri_int_os(lri_env,lrii,ra,rb,rab,obasa,obasb,fbasa,fbasb,&
                           iatom,jatom,ikind)
         ENDIF

         ! construct and invert S matrix
         lrii%sinv(1:nfa,1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa,1:nfa)
         IF(dab > lri_env%delta) THEN
           nn = nfa+nfb
           lrii%sinv(1:nfa,nfa+1:nn) = lrii%sab(1:nfa,1:nfb)
           lrii%sinv(nfa+1:nn,1:nfa) = TRANSPOSE(lrii%sab(1:nfa,1:nfb))
           lrii%sinv(nfa+1:nn,nfa+1:nn) = lri_env%bas_prop(jkind)%ri_ovlp(1:nfb,1:nfb)
         ENDIF
         CALL inverse_lri_overlap(lri_env,lrii%sinv)

         ! calculate Sinv*n and n*Sinv*n
         lrii%n(1:nfa) = lri_env%bas_prop(ikind)%int_fbas(1:nfa)
         IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           lrii%sn(1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),lrii%n(1:nfa))
           lrii%nsn = SUM(lrii%sn(1:nfa)*lrii%n(1:nfa))
         ELSE
           lrii%n(nfa+1:nn) = lri_env%bas_prop(jkind)%int_fbas(1:nfb)
           lrii%sn(1:nn) = MATMUL(lrii%sinv(1:nn,1:nn),lrii%n(1:nn))
           lrii%nsn = SUM(lrii%sn(1:nn)*lrii%n(1:nn))
         ENDIF

         ! calculate derivative of fit coefficients, needed for update of KS matrix
         IF(.NOT.dft_control%qs_control%lri_optbas) THEN
           CALL lri_calculate_derivative_acoef(lri_env,lrii,iatom,jatom,nba,nbb,&
                                               nfa,nfb,dab)
         ENDIF
 
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      IF(lri_env%debug) THEN
        CALL output_debug_info(lri_env,qs_env,lri_ints,soo_list)
      ENDIF

    END IF
 
    CALL timestop(handle)

  END SUBROUTINE calculate_lri_integrals

! *****************************************************************************
!> \brief calculates overlap integrals (aabb) of the orbital basis set,
!>        reguired for LRI basis set optimization
!> \param lri_env ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE calculate_lri_overlap_aabb(lri_env,qs_env)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_overlap_aabb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iac, iatom, ikind, &
                                                ilist, jatom, jkind, &
                                                jneighbor, nba, nbb, nkind, &
                                                nlist, nneighbor
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: obasa, obasb
    TYPE(lri_int_rho_type), POINTER          :: lriir
    TYPE(lri_list_type), POINTER             :: lri_ints_rho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(cell, lriir, lri_ints_rho, nl_iterator, obasa, obasb,&
            particle_set,soo_list)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      CALL get_qs_env(qs_env=qs_env,nkind=nkind,particle_set=particle_set,&
                      cell=cell)

      IF ( ASSOCIATED(lri_env%lri_ints_rho) ) THEN
        CALL deallocate_lri_ints_rho (lri_env%lri_ints_rho)
      END IF

      CALL allocate_lri_ints_rho(lri_env,lri_env%lri_ints_rho,nkind)
      lri_ints_rho => lri_env%lri_ints_rho

      CALL neighbor_list_iterator_create(nl_iterator,soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator)==0)

         CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
              nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
              iatom=iatom,jatom=jatom,r=rab)
        
         iac = ikind + nkind*(jkind - 1)
         dab = SQRT(SUM(rab*rab))

         obasa => lri_env%orb_basis(ikind)%gto_basis_set
         obasb => lri_env%orb_basis(jkind)%gto_basis_set
         IF (.NOT.ASSOCIATED(obasa)) CYCLE
         IF (.NOT.ASSOCIATED(obasb)) CYCLE

         lriir => lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)

         nba = obasa%nsgf
         nbb = obasb%nsgf
         ra(:) = pbc(particle_set(iatom)%r, cell) 
         rb(:) = pbc(particle_set(jatom)%r, cell) 
 
         ! calculate integrals (aa,bb) 
         CALL lri_int_aabb(lriir%soaabb,obasa,obasb,rab,ra,rb,lri_env%debug,&
                           lriir%dmax_aabb)

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_overlap_aabb

! *****************************************************************************
!> \brief performs the fitting of the density and distributes the fitted 
!>        density on the grid
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!>        lri_rho_struct where the fitted density is stored 
!> \param lri_density ...
!> \param qs_env ...
!> \param pmatrix ...
!> \param lri_rho_struct ...
!> \param atomic_kind_set ...
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE calculate_lri_densities(lri_env,lri_density,qs_env,pmatrix,&
                                    lri_rho_struct,atomic_kind_set,para_env)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(qs_rho_type), POINTER               :: lri_rho_struct
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_densities', &
      routineP = moduleN//':'//routineN

    CALL calculate_avec(lri_env,lri_density,qs_env,pmatrix)

    CALL distribute_lri_density_on_the_grid(lri_env,lri_density,qs_env,& 
                          lri_rho_struct,atomic_kind_set,para_env)   
 
   END SUBROUTINE calculate_lri_densities
 
! *****************************************************************************
!> \brief performs the fitting of the density; solves the linear system of  
!>        equations; yield the expansion coefficients avec
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!> \param lri_density ...
!> \param qs_env ...
!> \param pmatrix ...
! *****************************************************************************
  SUBROUTINE calculate_avec(lri_env,lri_density,qs_env,pmatrix)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_avec', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iac, iatom, ikind, ilist, ispin, jatom, jkind, &
      jneighbor, nba, nbb, nfa, nfb, nkind, nlist, nn, nneighbor, nspin, stat
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp)                            :: dab, rab(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pbij
    TYPE(cp_dbcsr_type), POINTER             :: pmat
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_list_type), POINTER             :: lri_rho
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(lrii, lri_rho, nl_iterator, pbij, pmat, soo_list)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      nspin = SIZE(pmatrix)
      nkind = lri_env%lri_ints%nkind

      CALL lri_density_release(lri_density)
      CALL lri_density_create(lri_density)
      lri_density%nspin = nspin

      ! allocate structure lri_rhos and vectors tvec and avec
      CALL allocate_lri_rhos(lri_env,lri_density%lri_rhos,nspin,nkind)

      DO ispin = 1, nspin
         pmat => pmatrix(ispin)%matrix
         lri_rho => lri_density%lri_rhos(ispin)%lri_list
 
         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,&
                 jatom=jatom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
                 r=rab)
  
            iac = ikind + nkind*(jkind - 1)
            dab = SQRT(SUM(rab*rab))

            IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

            ! get the density matrix Pab
            NULLIFY (pbij)
            IF (iatom <= jatom) THEN
              CALL cp_dbcsr_get_block_p(matrix=pmat,row=iatom,col=jatom,block=pbij,found=found)
              trans = .FALSE.
            ELSE
              CALL cp_dbcsr_get_block_p(matrix=pmat,row=jatom,col=iatom,block=pbij,found=found)
              trans = .TRUE.
            END IF
            CPASSERT(found)

            lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
  
            nba = lrii%nba
            nbb = lrii%nbb
            nfa = lrii%nfa
            nfb = lrii%nfb

            nn = nfa + nfb

            ! compute tvec = SUM_ab Pab *(a,b,x) and charge contraint 
            IF (trans) THEN
               lrho%charge = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%soo(1:nba,1:nbb)) 
               DO i=1,nfa
                  lrho%tvec(i) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%abaint(1:nba,1:nbb,i))
               END DO
               IF(dab > lri_env%delta) THEN
                 DO i=1,nfb
                    lrho%tvec(nfa+i) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%abbint(1:nba,1:nbb,i))
                 END DO
               ENDIF
            ELSE
               lrho%charge = SUM(pbij(1:nba,1:nbb)*lrii%soo(1:nba,1:nbb))
               DO i=1,nfa
                  lrho%tvec(i) = SUM(pbij(1:nba,1:nbb)*lrii%abaint(1:nba,1:nbb,i))
               END DO
               IF(dab > lri_env%delta) THEN
                 DO i=1,nfb
                    lrho%tvec(nfa+i) = SUM(pbij(1:nba,1:nbb)*lrii%abbint(1:nba,1:nbb,i))
                 END DO
               ENDIF
            END IF

            IF(iatom == jatom.AND.dab < lri_env%delta) THEN
              lrho%nst = SUM(lrho%tvec(1:nfa) * lrii%sn(1:nfa))
            ELSE
              lrho%nst = SUM(lrho%tvec(1:nn) * lrii%sn(1:nn))
            ENDIF
            lrho%lambda = (lrho%charge - lrho%nst)/lrii%nsn

            ! solve the linear system of equations
            ALLOCATE(m(nn),STAT=stat)
            CPASSERT(stat==0)
            m = 0._dp
            IF(iatom == jatom.AND.dab < lri_env%delta) THEN
              m(1:nfa) = lrho%tvec(1:nfa) + lrho%lambda * lrii%n(1:nfa)
              CALL dgemv("N",nfa,nfa,1.0_dp,lrii%sinv(1,1),nfa,&
                         m(1),1,0.0_dp,lrho%avec,1)
            ELSE
              m(1:nn) = lrho%tvec(1:nn) + lrho%lambda * lrii%n(1:nn)
              CALL dgemv("N",nn,nn,1.0_dp,lrii%sinv(1,1),nn,&
                         m(1),1,0.0_dp,lrho%avec,1)
            ENDIF
            DEALLOCATE(m,STAT=stat)
            CPASSERT(stat==0)
           
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

      END DO

      CALL set_qs_env(qs_env, lri_density=lri_density)

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_avec

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE distribute_lri_density_on_the_grid(lri_env,lri_density,qs_env,&
                                    lri_rho_struct,atomic_kind_set,para_env)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: lri_rho_struct
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'distribute_lri_density_on_the_grid', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iac, iatom, ikind, ilist, ispin, &
      jatom, jkind, jneighbor, nat, natom, nfa, nfb, nkind, nspin, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, rab(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: aci, acj, tot_rho_r
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: lri_coef
    TYPE(lri_list_type), POINTER             :: lri_rho
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(aci, acj, atomic_kind, atom_of_kind, lri_coef, lri_rho, &
         nl_iterator, soo_list, rho_r, rho_g, tot_rho_r)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      nspin = lri_density%nspin
      nkind = lri_env%lri_ints%nkind

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=nat)
      ALLOCATE(atom_of_kind(nat),STAT=stat)
      CPASSERT(stat==0)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      ! allocate the arrays to hold RI expansion coefficients lri_coefs
      CALL allocate_lri_coefs(lri_env,lri_density,atomic_kind_set)
      DO ispin = 1, nspin

         lri_coef => lri_density%lri_coefs(ispin)%lri_kinds  
         lri_rho  => lri_density%lri_rhos(ispin)%lri_list

         ! sum up expansion coefficients
         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
                 iatom=iatom,jatom=jatom,ilist=ilist,inode=jneighbor,r=rab)
            dab = SQRT(SUM(rab*rab))
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            aci => lri_coef(ikind)%acoef(atom_a,:)
            acj => lri_coef(jkind)%acoef(atom_b,:)
            iac = ikind + nkind*(jkind - 1)
            lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            nfa = lrho%nfa
            nfb = lrho%nfb
            IF(iatom == jatom.AND.dab < lri_env%delta) THEN
              !self pair aa
              aci(1:nfa) = aci(1:nfa) + lrho%avec(1:nfa)
            ELSE
              IF(iatom == jatom) THEN
                !periodic self pair aa' 
                aci(1:nfa) = aci(1:nfa) + lrho%avec(1:nfa)
                acj(1:nfb) = acj(1:nfb) + lrho%avec(nfa+1:nfa+nfb)
              ELSE
                !pairs ab
                aci(1:nfa) = aci(1:nfa) + 2.0_dp*lrho%avec(1:nfa)
                acj(1:nfb) = acj(1:nfb) + 2.0_dp*lrho%avec(nfa+1:nfa+nfb)
              ENDIF
            ENDIF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         ! replicate the acoef infomation 
         DO ikind=1,nkind
           atomic_kind => atomic_kind_set(ikind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
           DO iatom=1,natom
              aci => lri_coef(ikind)%acoef(iatom,:)
              CALL mp_sum(aci,para_env%group)
           END DO
         END DO

      END DO

      !distribute fitted density on the grid
      CALL qs_rho_get(lri_rho_struct, rho_r=rho_r, rho_g=rho_g,&
                      tot_rho_r=tot_rho_r)
      DO ispin=1,nspin
       CALL calculate_lri_rho_elec(rho_g(ispin),&
                                   rho_r(ispin), qs_env, &
                                   lri_density%lri_coefs(ispin)%lri_kinds,&
                                   tot_rho_r(ispin))
      ENDDO

      CALL set_qs_env(qs_env, lri_density=lri_density)

      DEALLOCATE(atom_of_kind,STAT=stat)
      CPASSERT(stat==0)

    END IF

    CALL timestop(handle)

  END SUBROUTINE distribute_lri_density_on_the_grid

! *****************************************************************************
!> \brief calculates the lri intergrals according to the Obara-Saika (OS)
!>        scheme
!> \param lri_env ...
!> \param lrii ...
!> \param ra position atom A
!> \param rb position atom B
!> \param rab distance vector
!> \param obasa orb basis on center A
!> \param obasb orb basis on center B
!> \param fbasa aux basis on center A
!> \param fbasb aux basis on center B
!> \param iatom index atom A 
!> \param jatom index atom B
!> \param ikind kind atom A
! *****************************************************************************
  SUBROUTINE lri_int_os(lri_env,lrii,ra,rb,rab,obasa,obasb,fbasa,fbasb,&
                        iatom,jatom,ikind)
    
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_int_type), POINTER              :: lrii
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: obasa, obasb, fbasa, fbasb
    INTEGER, INTENT(IN)                      :: iatom, jatom, ikind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_os', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nba, nbb, nfa, nfb
    REAL(KIND=dp)                            :: dab

    CALL timeset(routineN,handle)

    dab = SQRT(SUM(rab*rab))
    nba = obasa%nsgf
    nbb = obasb%nsgf
    nfa = fbasa%nsgf
    nfb = fbasb%nsgf

    !*** calculate overlap integrals; for iatom=jatom this is the self-overlap
    IF(iatom == jatom.AND.dab < lri_env%delta) THEN
      !*** (fa,fa)
      lrii%sab(1:nfa,1:nfa)=lri_env%bas_prop(ikind)%ri_ovlp(1:nfa,1:nfa)
      lrii%dsab = 0._dp
      !*** (a,a)
      lrii%soo(1:nba,1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba,1:nba)
      lrii%dsoo = 0._dp
      !*** (a,a,fa)
      CALL lri_int_aba(lrii%abaint,ra=ra,rb=rb,rab=rab,oba=obasa,obb=obasb,&
                       fba=fbasa, calculate_forces=.FALSE.,debug=lri_env%debug,&
                       dmax=lrii%dmax_aba)
      lrii%dabdaint = 0.0_dp
      lrii%dabbint  = 0.0_dp
    ELSE
      !*** (fa,fb)
      CALL lri_int_ab(lrii%sab,lrii%dsab,ra,rb,rab,fbasa,fbasb,&
                      lrii%calc_force_pair,lri_env%debug,&
                      lrii%dmax_ab)
      !*** (a,b)
      CALL lri_int_ab(lrii%soo,lrii%dsoo,ra,rb,rab,obasa,obasb,&
                      lrii%calc_force_pair,lri_env%debug,&
                      lrii%dmax_oo)
      !*** (a,b,fa)
      CALL lri_int_aba(lrii%abaint,lrii%dabdaint,ra,rb,rab,obasa,obasb,fbasa,&
                       lrii%calc_force_pair,lri_env%debug,lrii%dmax_aba)
      !*** (a,b,fb)
      CALL lri_int_abb(lrii%abbint,lrii%dabbint,ra,rb,rab,obasa,obasb,fbasb,&
                       lrii%calc_force_pair,lri_env%debug,lrii%dmax_abb)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE lri_int_os
 
! *****************************************************************************
!> \brief calcuates the lri integrals using solid harmonic Gaussians
!> \param lri_env ...
!> \param lrii ...
!> \param ra position of center A
!> \param rb position of center B
!> \param rab distance vector
!> \param obasa orb basis on A
!> \param obasb orb basis on B
!> \param fbasa aux basis on A
!> \param fbasb aux basis on B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param jkind kind atom B
! *****************************************************************************
  SUBROUTINE lri_int_shg(lri_env,lrii,ra,rb,rab,obasa,obasb,fbasa,fbasb,&
                        iatom,jatom,ikind,jkind)
    
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_int_type), POINTER              :: lrii
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: obasa, obasb, fbasa, fbasb
    INTEGER, INTENT(IN)                      :: iatom, jatom, ikind, jkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_shg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nba, nbb, nfa, nfb, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: scon_orba, scon_orbb, &
                                                scon_ria, scon_rib
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: dWaux_mat

    CALL timeset(routineN,handle)
    NULLIFY(scon_orba,scon_orbb,scon_ria,scon_rib,Waux_mat,dWaux_mat)
    failure = .FALSE.

    dab = SQRT(SUM(rab*rab))
    nba = obasa%nsgf
    nbb = obasb%nsgf
    nfa = fbasa%nsgf
    nfb = fbasb%nsgf

    !*** calculate overlap integrals; for iatom=jatom this is the self-overlap
    IF(iatom == jatom.AND.dab < lri_env%delta) THEN
      !*** (fa,fa)
      lrii%sab(1:nfa,1:nfa)=lri_env%bas_prop(ikind)%ri_ovlp(1:nfa,1:nfa)
      lrii%dsab = 0._dp
      !*** (a,a)
      lrii%soo(1:nba,1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba,1:nba)
      lrii%dsoo = 0._dp
      !*** (a,a,fa)
      CALL lri_int_aba(lrii%abaint,ra=ra,rb=rb,rab=rab,oba=obasa,obb=obasb,&
                       fba=fbasa, calculate_forces=.FALSE.,debug=lri_env%debug,&
                       dmax=lrii%dmax_aba)
      lrii%dabdaint = 0.0_dp
      lrii%dabbint  = 0.0_dp
    ELSE
      scon_orba   => lri_env%bas_prop(ikind)%scon_orb       
      scon_orbb   => lri_env%bas_prop(jkind)%scon_orb       
      scon_ria    => lri_env%bas_prop(ikind)%scon_ri       
      scon_rib    => lri_env%bas_prop(jkind)%scon_ri       
      CALL precalc_angular_shg_part(obasa,obasb,fbasa,fbasb,rab,Waux_mat,dWaux_mat,&
                                    lrii%calc_force_pair) 
      !*** (fa,fb)
      CALL lri_int_ab_shg(lrii%sab,lrii%dsab,rab,fbasa,fbasb,scon_ria,scon_rib,&
                          Waux_mat,dWaux_mat,lrii%calc_force_pair)
      !*** (a,b)
      CALL lri_int_ab_shg(lrii%soo,lrii%dsoo,rab,obasa,obasb,scon_orba,scon_orbb,&
                          Waux_mat,dWaux_mat,lrii%calc_force_pair)
      !*** (a,b,fa)
       CALL lri_int_aba_shg(lri_env,lrii%abaint,lrii%dabdaint,rab,obasa,obasb,fbasa,&
                            scon_orba,scon_orbb,scon_ria,Waux_mat,dWaux_mat,&
                            lrii%calc_force_pair)
      !*** (a,b,fb)
      CALL lri_int_abb_shg(lri_env,lrii%abbint,lrii%dabbint,rab,obasa,obasb,fbasb,&
                           scon_orba,scon_orbb,scon_rib,Waux_mat,dWaux_mat,&
                           lrii%calc_force_pair)
      DEALLOCATE(Waux_mat,dWaux_mat,STAT=stat)
      CPASSERT(stat==0)
    ENDIF

    CALL timestop(handle)

   END SUBROUTINE lri_int_shg

! *****************************************************************************
!> \brief calculate derivate of fit coefficients acoef with respect to 
!>        density matrix pmatrix 
!>        R = (a,b)/nsn  - SUM_i (a,b,ai)*sn(i)/nsn 
!>        Q = SUM_i sinv*(a,b,ai)
!>        derviate_aci = R + Q
!> \param lri_env ...
!> \param lrii ...
!> \param iatom ...
!> \param jatom ...
!> \param nba number of primary basis functions on a
!> \param nbb number of primary basis functions on b
!> \param nfa number of ri basis functions on a
!> \param nfb number of ri basis functions on b
!> \param dab distance between center a and b
! *****************************************************************************
  SUBROUTINE lri_calculate_derivative_acoef(lri_env,lrii,iatom,jatom,nba,nbb,&
                                            nfa,nfb,dab)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_int_type), POINTER              :: lrii
    INTEGER, INTENT(IN)                      :: iatom, jatom, nba, nbb, nfa, &
                                                nfb
    REAL(KIND=dp), INTENT(IN)                :: dab

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'lri_calculate_derivative_acoef', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nn, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dlambdadp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: abxint_work, Q, R

    failure = .FALSE.
    CALL timeset(routineN,handle)

    nn = nfa + nfb
 
    ALLOCATE(R(nba,nbb,nn),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(dlambdadp(nba,nbb),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(Q(nba,nbb,nn),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(abxint_work(nba,nbb,nn),STAT=stat)
    CPASSERT(stat==0)

    R(:,:,:)  = 0._dp
    dlambdadp = 0._dp
    Q(:,:,:)  = 0._dp
    abxint_work(:,:,:) = 0._dp

    abxint_work(1:nba,1:nbb,1:nfa)=lrii%abaint(1:nba,1:nbb,1:nfa)    !abaint
    abxint_work(1:nba,1:nbb,nfa+1:nn)=lrii%abbint(1:nba,1:nbb,1:nfb) !abbint
    IF(iatom == jatom.AND.dab < lri_env%delta) THEN
      ! R = Sinv*abx_int 
      CALL dgemm("N","N",nba*nbb,nfa,nfa,1.0_dp,abxint_work,nba*nbb,&
                 lrii%sinv,nfa,0.0_dp,R,nba*nbb)
      ! dlambdadp = soo/nsn - sn*abx_int/nsn
      dlambdadp(1:nba,1:nbb) = lrii%soo(1:nba,1:nbb)
      CALL dgemm("N","N",nba*nbb,1,nfa,-1.0_dp/lrii%nsn,abxint_work,nba*nbb,&
                 lrii%sn,nfa,1.0_dp/lrii%nsn,dlambdadp,nba*nbb)
      ! Q = dlambdadp * sn
      CALL dgemm("N","T",nba*nbb,nfa,1,1.0_dp,dlambdadp,nba*nbb,&
                 lrii%sn,nfa,0.0_dp,Q,nba*nbb)
      lrii%dacoef(1:nba,1:nbb,1:nfa)= R(1:nba,1:nbb,1:nfa)+Q(1:nba,1:nbb,1:nfa)
    ELSE
      ! R = Sinv*abx_int 
      CALL dgemm("N","N",nba*nbb,nn,nn,1.0_dp,abxint_work,nba*nbb,&
                 lrii%sinv,nn,0.0_dp,R,nba*nbb)
      ! dlambdadp = soo/nsn - sn*abx_int/nsn
      dlambdadp(1:nba,1:nbb) = lrii%soo(1:nba,1:nbb)
      CALL dgemm("N","N",nba*nbb,1,nn,-1.0_dp/lrii%nsn,abxint_work,nba*nbb,&
                 lrii%sn,nn,1.0_dp/lrii%nsn,dlambdadp,nba*nbb)
      ! Q = dlambdadp * sn
      CALL dgemm("N","T",nba*nbb,nn,1,1.0_dp,dlambdadp,nba*nbb,&
                 lrii%sn,nn,0.0_dp,Q,nba*nbb)
      lrii%dacoef(1:nba,1:nbb,1:nn)= R(1:nba,1:nbb,1:nn)+Q(1:nba,1:nbb,1:nn)
    ENDIF
    
    DEALLOCATE(abxint_work,STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(R,STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(dlambdadp,STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(Q,STAT=stat)
    CPASSERT(stat==0)

    CALL timestop(handle)
  
  END SUBROUTINE lri_calculate_derivative_acoef

! *****************************************************************************
!> \brief get inverse or pseudoinverse of lri overlap matrix for aux basis set
!> \param lri_env lri environment
!> \param sinv on entry overlap matrix, on exit its inverse
! *****************************************************************************
  SUBROUTINE inverse_lri_overlap(lri_env,sinv)

    TYPE(lri_environment_type)               :: lri_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sinv

    CHARACTER(LEN=*), PARAMETER :: routineN = 'inverse_lri_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, n, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: anorm, rcond, rskip
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: s
    REAL(KIND=dp), EXTERNAL                  :: dlange

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(s)

    rskip = 1.E-8_dp ! parameter for pseudo inverse
    n=SIZE(sinv,1)

    SELECT CASE(lri_env%lri_overlap_inv)
    CASE(do_lri_inv)
       CALL invmat_symm(sinv)
    CASE(do_lri_pseudoinv_svd)
       ALLOCATE(s(n,n),STAT=stat)
       CPASSERT(stat==0)
       s(:,:) = sinv
       CALL get_pseudo_inverse_svd(s,sinv,rskip)
       DEALLOCATE(s,STAT=stat)
       CPASSERT(stat==0)
    CASE(do_lri_pseudoinv_diag)
       ALLOCATE(s(n,n),STAT=stat)
       CPASSERT(stat==0)
       s(:,:) = sinv
       CALL get_pseudo_inverse_diag(s,sinv,rskip)
       DEALLOCATE(s,STAT=stat)
       CPASSERT(stat==0)
    CASE(do_lri_inv_auto)
       ! decide whether to calculate inverse or pseudoinverse
       ALLOCATE(s(n,n),STAT=stat)
       CPASSERT(stat==0)
       s(:,:) = sinv
       ALLOCATE(work(3*n),iwork(n),STAT=stat)
       CPASSERT(stat==0)
       ! norm of matrix
       anorm = dlange('1',n,n,sinv,n,work)
       ! Cholesky factorization (fails if matrix not positive definite)
       CALL dpotrf('U',n,sinv,n,info)
       IF (info == 0) THEN
           ! condition number
           CALL dpocon('U',n,sinv,n,anorm,rcond,work,iwork,info)
           IF (info /= 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,"DPOCON failed")
           ENDIF
           IF(LOG(1._dp/rcond) > lri_env%cond_max) THEN
             CALL get_pseudo_inverse_svd(s,sinv,rskip)
           ELSE
             CALL invmat_symm(sinv,"U")
           ENDIF
       ELSE
           CALL get_pseudo_inverse_svd(s,sinv,rskip)
       END IF
       DEALLOCATE(s,work,iwork,STAT=stat)
       CPASSERT(stat==0)
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,&
                 routineP,"No initialization for LRI overlap available?",&
                 only_ionode=.TRUE.)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE inverse_lri_overlap

! *****************************************************************************
!> \brief debug output 
!> \param lri_env ...
!> \param qs_env ...
!> \param lri_ints ...
!> \param soo_list ...
! *****************************************************************************
  SUBROUTINE output_debug_info(lri_env,qs_env,lri_ints,soo_list)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    CHARACTER(LEN=*), PARAMETER :: routineN = 'output_debug_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iac, ikind, ilist, &
                                                iunit, jkind, jneighbor, nkind
    REAL(KIND=dp)                            :: dmax_aabb, dmax_ab, dmax_aba, &
                                                dmax_abb, dmax_oo
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(lri_int_rho_type), POINTER          :: lriir
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
    NULLIFY(input, logger, lrii, lriir, nl_iterator, para_env)
    CALL get_qs_env(qs_env,dft_control=dft_control,input=input,nkind=nkind,&
                    para_env=para_env)
    dmax_ab   = 0._dp
    dmax_oo   = 0._dp
    dmax_aba  = 0._dp
    dmax_abb  = 0._dp
    dmax_aabb = 0._dp

    CALL neighbor_list_iterator_create(nl_iterator,soo_list)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)

       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
            ilist=ilist,inode=jneighbor)

       iac = ikind + nkind*(jkind - 1)
       lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

       dmax_ab  = MAX(dmax_ab,lrii%dmax_ab)
       dmax_oo  = MAX(dmax_oo,lrii%dmax_oo)
       dmax_aba = MAX(dmax_aba,lrii%dmax_aba)
       dmax_abb = MAX(dmax_abb,lrii%dmax_abb)

       IF(dft_control%qs_control%lri_optbas) THEN
         lriir => lri_env%lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)
         dmax_aabb = MAX(dmax_aabb,lriir%dmax_aabb)
       ENDIF

    END DO

    CALL neighbor_list_iterator_release(nl_iterator)
    CALL mp_max(dmax_ab,para_env%group)
    CALL mp_max(dmax_oo,para_env%group)
    CALL mp_max(dmax_aba,para_env%group)
    CALL mp_max(dmax_abb,para_env%group)
    CALL mp_max(dmax_aabb,para_env%group)

    logger => cp_get_default_logger()
    iunit=cp_print_key_unit_nr(logger,input,"PRINT%PROGRAM_RUN_INFO",&
                               extension=".lridebug")

    IF (iunit > 0) THEN
       WRITE(iunit,FMT="(/,T2,A)") "DEBUG INFO FOR LRI INTEGRALS"
       WRITE(iunit,FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "//&
                                            "[ai|bi]; fit basis", dmax_ab
       WRITE(iunit,FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "//&
                                            "[a|b]; orbital basis", dmax_oo
       WRITE(iunit,FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "//&
                                            "[a|b|ai]", dmax_aba
       WRITE(iunit,FMT="(T2,A,T69,ES12.5)") "Maximal deviation of integrals "//&
                                            "[a|b|bi]", dmax_abb
       IF(dft_control%qs_control%lri_optbas) THEN
          WRITE(iunit,FMT="(T2,A,T69,ES12.5,/)") "Maximal deviation of integrals "//&
                                               "[aa|bb]; orbital basis",&
                                               dmax_aabb
       ENDIF
    ENDIF

    CALL cp_print_key_finished_output(iunit,logger,input,&
                    "PRINT%PROGRAM_RUN_INFO")
    CALL timestop(handle)
  
  END SUBROUTINE output_debug_info

END MODULE lri_environment_methods
