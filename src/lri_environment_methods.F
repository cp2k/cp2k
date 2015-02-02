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
  USE ai_overlap,                      ONLY: overlap
  USE ai_overlap3,                     ONLY: overlap3
  USE ai_overlap_aabb,                 ONLY: overlap_aabb
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE lri_debug_integrals,             ONLY: overlap_aabb_test,&
                                             overlap_ab_test,&
                                             overlap_abc_test
  USE lri_environment_types,           ONLY: &
       allocate_lri_coefs, allocate_lri_ints, allocate_lri_ints_rho, &
       allocate_lri_rhos, deallocate_lri_ints, deallocate_lri_ints_rho, &
       lri_density_create, lri_density_release, lri_density_type, &
       lri_environment_type, lri_int_rho_type, lri_int_type, lri_kind_type, &
       lri_list_type, lri_rhoab_type
  USE mathlib,                         ONLY: invmat
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE orbital_pointers,                ONLY: ncoset
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_lri_matrices(lri_env,qs_env,calculate_forces,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_lri_matrices', &
      routineP = moduleN//':'//routineN

! calculate the integrals needed to do the local (2-center) expansion
! of the (pair) densities

    CALL calculate_lri_integrals(lri_env,qs_env,calculate_forces,error)

  END SUBROUTINE build_lri_matrices

! *****************************************************************************
!> \brief calculates integrals needed for the LRI density fitting,
!>        integrals are calculated once, before the SCF starts
!> \param lri_env ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_lri_integrals(lri_env,qs_env,calculate_forces,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iac, iatom, ikind, ilist, jatom, jkind, jneighbor, &
      nba, nbb, nfa, nfb, nkind, nlist, nn, nneighbor, stat
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
            obasa, obasb, particle_set,soo_list, virial)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
                      nkind=nkind,particle_set=particle_set,virial=virial,&
                      error=error)

      IF ( ASSOCIATED(lri_env%lri_ints) ) THEN
        CALL deallocate_lri_ints (lri_env%lri_ints,error)
      END IF

      ! allocate matrices storing the LRI integrals
      CALL allocate_lri_ints(lri_env,lri_env%lri_ints,nkind,&
                             calculate_forces,virial,error)
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

         ! calculate integrals (a,b,fa) and (a,b,fb)
         IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           CALL lri_int_aba(lrii%abaint,ra=ra,rb=rb,rab=rab,oba=obasa,obb=obasb,&
                            fba=fbasa, calculate_forces=.FALSE.,debug=lri_env%debug,&
                            dmax=lrii%dmax_aba,error=error)
           lrii%dabdaint = 0.0_dp
           lrii%dabbint  = 0.0_dp
         ELSE
           CALL lri_int_aba(lrii%abaint,lrii%dabdaint,ra,rb,rab,obasa,obasb,fbasa,&
                            lrii%calc_force_pair,lri_env%debug,lrii%dmax_aba,&
                            error=error)
           CALL lri_int_abb(lrii%abbint,lrii%dabbint,ra,rb,rab,obasa,obasb,fbasb,&
                            lrii%calc_force_pair,lri_env%debug,lrii%dmax_abb,&
                            error=error)
         ENDIF

         ! calculate integrals (fa,fb); for iatom=jatom this is the self-overlap
         IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           lrii%sab(1:nfa,1:nfa)=lri_env%bas_ovlp(ikind)%ri_ovlp(1:nfa,1:nfa)
           lrii%dsab = 0._dp
         ELSE
           CALL lri_int_ab(lrii%sab,lrii%dsab,ra,rb,rab,fbasa,fbasb,&
                           lrii%calc_force_pair,lri_env%debug,&
                           lrii%dmax_ab,error)
         ENDIF

         ! construct and invert S matrix
         lrii%sinv(1:nfa,1:nfa) = lri_env%bas_ovlp(ikind)%ri_ovlp(1:nfa,1:nfa)
         IF(dab > lri_env%delta) THEN
           nn = nfa+nfb
           lrii%sinv(1:nfa,nfa+1:nn) = lrii%sab(1:nfa,1:nfb)
           lrii%sinv(nfa+1:nn,1:nfa) = TRANSPOSE(lrii%sab(1:nfa,1:nfb))
           lrii%sinv(nfa+1:nn,nfa+1:nn) = lri_env%bas_ovlp(jkind)%ri_ovlp(1:nfb,1:nfb)
         ENDIF
         CALL invmat(lrii%sinv,stat,error)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ! calculate Sinv*n and n*Sinv*n
         lrii%n(1:nfa) = lri_env%bas_int(ikind)%int_fbas(1:nfa)
         IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           lrii%sn(1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),lrii%n(1:nfa))
           lrii%nsn = SUM(lrii%sn(1:nfa)*lrii%n(1:nfa))
         ELSE
           lrii%n(nfa+1:nn) = lri_env%bas_int(jkind)%int_fbas(1:nfb)
           lrii%sn(1:nn) = MATMUL(lrii%sinv(1:nn,1:nn),lrii%n(1:nn))
           lrii%nsn = SUM(lrii%sn(1:nn)*lrii%n(1:nn))
         ENDIF

         ! calculate integrals (a,b), overlap of primary basis
         IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           lrii%soo(1:nba,1:nba) = lri_env%bas_ovlp(ikind)%orb_ovlp(1:nba,1:nba)
           lrii%dsoo = 0._dp
         ELSE
           CALL lri_int_ab(lrii%soo,lrii%dsoo,ra,rb,rab,obasa,obasb,&
                           lrii%calc_force_pair,lri_env%debug,&
                           lrii%dmax_oo,error=error)
         ENDIF

         ! calculate derivative of fit coefficients, needed for update of KS matrix
         IF(.NOT.dft_control%qs_control%lri_optbas) THEN
           CALL lri_calculate_derivative_acoef(lri_env,lrii,iatom,jatom,nba,nbb,&
                                               nfa,nfb,dab,error)
         ENDIF
 
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      IF(lri_env%debug) THEN
        CALL output_debug_info(lri_env,qs_env,lri_ints,soo_list,error)
      ENDIF

    END IF
 
    CALL timestop(handle)

  END SUBROUTINE calculate_lri_integrals

! *****************************************************************************
!> \brief calculates overlap integrals (aabb) of the orbital basis set,
!>        reguired for LRI basis set optimization
!> \param lri_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_lri_overlap_aabb(lri_env,qs_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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
                      cell=cell,error=error)

      IF ( ASSOCIATED(lri_env%lri_ints_rho) ) THEN
        CALL deallocate_lri_ints_rho (lri_env%lri_ints_rho,error)
      END IF

      CALL allocate_lri_ints_rho(lri_env,lri_env%lri_ints_rho,nkind,error)
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
                           lriir%dmax_aabb,error)

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
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_lri_densities(lri_env,lri_density,qs_env,pmatrix,&
                                    lri_rho_struct,atomic_kind_set,para_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(qs_rho_type), POINTER               :: lri_rho_struct
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_densities', &
      routineP = moduleN//':'//routineN

    CALL calculate_avec(lri_env,lri_density,qs_env,pmatrix,error)

    CALL distribute_lri_density_on_the_grid(lri_env,lri_density,qs_env,& 
                          lri_rho_struct,atomic_kind_set,para_env,error)   
 
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_avec(lri_env,lri_density,qs_env,pmatrix,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

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

      CALL lri_density_release(lri_density,error)
      CALL lri_density_create(lri_density,error)
      lri_density%nspin = nspin

      ! allocate structure lri_rhos and vectors tvec and avec
      CALL allocate_lri_rhos(lri_env,lri_density%lri_rhos,nspin,nkind,error)

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
            CPPostcondition(found,cp_failure_level,routineP,error,failure)

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
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            m = 0._dp
            IF(iatom == jatom.AND.dab < lri_env%delta) THEN
              m(1:nfa) = lrho%tvec(1:nfa) + lrho%lambda * lrii%n(1:nfa)
              lrho%avec(1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),m(1:nfa))
            ELSE
              m(1:nn) = lrho%tvec(1:nn) + lrho%lambda * lrii%n(1:nn)
              lrho%avec(1:nn)  = MATMUL(lrii%sinv(1:nn,1:nn),m(1:nn))
            ENDIF
            DEALLOCATE(m,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

      END DO

      CALL set_qs_env(qs_env, lri_density=lri_density, error=error)

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
!> \param error ...
! *****************************************************************************
  SUBROUTINE distribute_lri_density_on_the_grid(lri_env,lri_density,qs_env,&
                                    lri_rho_struct,atomic_kind_set,para_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: lri_rho_struct
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      ! allocate the arrays to hold RI expansion coefficients lri_coefs
      CALL allocate_lri_coefs(lri_env,lri_density,atomic_kind_set,error)
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
      CALL qs_rho_get(lri_rho_struct, rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r, error=error)
      DO ispin=1,nspin
       CALL calculate_lri_rho_elec(rho_g(ispin),&
                                   rho_r(ispin), qs_env, lri_env,&
                                   lri_density%lri_coefs(ispin)%lri_kinds,&
                                   tot_rho_r(ispin), error)
      ENDDO

      CALL set_qs_env(qs_env, lri_density=lri_density, error=error)

      DEALLOCATE(atom_of_kind,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE distribute_lri_density_on_the_grid

! *****************************************************************************
!> \brief calculate integrals (a,b,fa)
!> \param abaint integral (a,b,fa)
!> \param dabdaint derivative of abaint with respect to A
!> \param ra ...
!> \param rb ...
!> \param rab ...
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fba auxiliary basis set at center A
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_aba(abaint,dabdaint,ra,rb,rab,oba,obb,fba,&
                         calculate_forces,debug,dmax,error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abaint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: dabdaint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fba
    LOGICAL, INTENT(IN)                      :: calculate_forces, debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_aba', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iset, jset, kaset, &
                                                m1, m2, m3, ncoa, ncob, ncoc, &
                                                nseta, nsetb, nsetca, sgfa, &
                                                sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, lca_max, lca_min, &
                                                npgfa, npgfb, npgfca, nsgfa, &
                                                nsgfb, nsgfca
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfca
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: saba
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sabda, sdaba, sdabda
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, rpgfca, sphi_a, &
                                                sphi_b, sphi_ca, zeta, zetb, &
                                                zetca

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, lca_max, lca_min, npgfa, npgfb,&
             npgfca, nsgfa, nsgfb, nsgfca) 
    NULLIFY(first_sgfa, first_sgfb, first_sgfca, set_radius_a, set_radius_b,&
            rpgfa, rpgfb, rpgfca, sphi_a, sphi_b, sphi_ca, zeta, zetb, zetca) 

    ! basis ikind
    first_sgfa   =>  oba%first_sgf
    la_max       =>  oba%lmax
    la_min       =>  oba%lmin
    npgfa        =>  oba%npgf
    nseta        =   oba%nset
    nsgfa        =>  oba%nsgf_set
    rpgfa        =>  oba%pgf_radius
    set_radius_a =>  oba%set_radius
    sphi_a       =>  oba%sphi
    zeta         =>  oba%zet
    ! basis jkind
    first_sgfb   =>  obb%first_sgf
    lb_max       =>  obb%lmax
    lb_min       =>  obb%lmin
    npgfb        =>  obb%npgf
    nsetb        =   obb%nset
    nsgfb        =>  obb%nsgf_set
    rpgfb        =>  obb%pgf_radius
    set_radius_b =>  obb%set_radius
    sphi_b       =>  obb%sphi
    zetb         =>  obb%zet

    ! basis RI A
    first_sgfca  =>  fba%first_sgf
    lca_max      =>  fba%lmax
    lca_min      =>  fba%lmin
    npgfca       =>  fba%npgf
    nsetca       =   fba%nset
    nsgfca       =>  fba%nsgf_set
    rpgfca       =>  fba%pgf_radius
    sphi_ca      =>  fba%sphi
    zetca        =>  fba%zet

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1

          ! calculate integrals abaint and derivative [d(a,b,a)/dA] dabdaint if requested
          rac = 0._dp
          dac = 0._dp
          rbc = -rab
          dbc = dab
          DO kaset=1,nsetca
             ncoc = npgfca(kaset)*ncoset(lca_max(kaset))
             sgfc = first_sgfca(1,kaset)
             m3=sgfc+nsgfca(kaset)-1
             IF(ncoa*ncob*ncoc > 0) THEN
                ALLOCATE(saba(ncoa,ncob,ncoc),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                saba(:,:,:)    = 0._dp
                ! integrals
                IF(calculate_forces) THEN
                   ALLOCATE(sdaba(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE(sabda(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE(sdabda(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sdaba(:,:,:,:) = 0._dp
                   sabda(:,:,:,:) = 0._dp
                   sdabda(:,:,:,:) = 0._dp
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                 lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                 lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                                 rab,dab,rac,dac,rbc,dbc,saba,sdaba,sabda,error=error)
                   !d(a,b,a)/dA = (da/dA,b,a) + (a,b,da/dA)
                   sdabda(:,:,:,:) = sdaba + sabda

                   DO i=1,3
                    CALL abc_contract(dabdaint(sgfa:m1,sgfb:m2,sgfc:m3,i),sdabda(:,:,:,i),& 
                         sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                         ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset),error)
                   ENDDO

                   DEALLOCATE(sdaba,sabda,sdabda,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ELSE
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                 lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                 lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                                 rab,dab,rac,dac,rbc,dbc,saba,error=error)
                ENDIF
                ! debug if requested
                IF(debug) THEN
                   CALL overlap_abc_test(la_max(iset),npgfa(iset),zeta(:,iset),la_min(iset),&
                                         lb_max(jset),npgfb(jset),zetb(:,jset),lb_min(jset),&
                                         lca_max(kaset),npgfca(kaset),zetca(:,kaset),lca_min(kaset),&
                                         ra,rb,ra,saba,dmax,error)
                ENDIF
                CALL abc_contract(abaint(sgfa:m1,sgfb:m2,sgfc:m3),saba,& 
                     sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset),error)
                DEALLOCATE(saba,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE lri_int_aba

! *****************************************************************************
!> \brief calculate integrals (a,b,fb)
!> \param abbint integral (a,b,fb)
!> \param dabbint derivative of abbint with respect to A
!> \param ra ...
!> \param rb ...
!> \param rab ...
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fbb auxiliary basis set at center B
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_abb(abbint,dabbint,ra,rb,rab,oba,obb,fbb,calculate_forces,&
                         debug,dmax,error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abbint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: dabbint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fbb
    LOGICAL, INTENT(IN)                      :: calculate_forces, debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_abb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iset, jset, kbset, &
                                                m1, m2, m3, ncoa, ncob, ncoc, &
                                                nseta, nsetb, nsetcb, sgfa, &
                                                sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, lcb_max, lcb_min, &
                                                npgfa, npgfb, npgfcb, nsgfa, &
                                                nsgfb, nsgfcb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfcb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: sabb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sdabb
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, rpgfcb, sphi_a, &
                                                sphi_b, sphi_cb, zeta, zetb, &
                                                zetcb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, lcb_max, lcb_min, npgfa, npgfb,&
            npgfcb, nsgfa, nsgfb, nsgfcb)
    NULLIFY(first_sgfa, first_sgfb, first_sgfcb, set_radius_a, set_radius_b,&
            rpgfa, rpgfb, rpgfcb, sphi_a, sphi_b, sphi_cb, zeta, zetb, zetcb) 

    ! basis ikind
    first_sgfa   =>  oba%first_sgf
    la_max       =>  oba%lmax
    la_min       =>  oba%lmin
    npgfa        =>  oba%npgf
    nseta        =   oba%nset
    nsgfa        =>  oba%nsgf_set
    rpgfa        =>  oba%pgf_radius
    set_radius_a =>  oba%set_radius
    sphi_a       =>  oba%sphi
    zeta         =>  oba%zet
    ! basis jkind
    first_sgfb   =>  obb%first_sgf
    lb_max       =>  obb%lmax
    lb_min       =>  obb%lmin
    npgfb        =>  obb%npgf
    nsetb        =   obb%nset
    nsgfb        =>  obb%nsgf_set
    rpgfb        =>  obb%pgf_radius
    set_radius_b =>  obb%set_radius
    sphi_b       =>  obb%sphi
    zetb         =>  obb%zet

    ! basis RI B
    first_sgfcb  =>  fbb%first_sgf
    lcb_max      =>  fbb%lmax
    lcb_min      =>  fbb%lmin
    npgfcb       =>  fbb%npgf
    nsetcb       =   fbb%nset
    nsgfcb       =>  fbb%nsgf_set
    rpgfcb       =>  fbb%pgf_radius
    sphi_cb      =>  fbb%sphi
    zetcb        =>  fbb%zet

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1

          ! calculate integrals abbint and derivative [d(a,b,b)/dA] dabbint if requested
          rac = rab
          dac = dab
          rbc = 0._dp
          dbc = 0._dp
          DO kbset=1,nsetcb
             ncoc = npgfcb(kbset)*ncoset(lcb_max(kbset))
             sgfc = first_sgfcb(1,kbset)
             m3=sgfc+nsgfcb(kbset)-1
             IF(ncoa*ncob*ncoc > 0) THEN
                ALLOCATE(sabb(ncoa,ncob,ncoc),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                sabb(:,:,:) = 0._dp
                IF(calculate_forces) THEN
                   ALLOCATE(sdabb(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sdabb(:,:,:,:) = 0._dp
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                                rab,dab,rac,dac,rbc,dbc,sabb,sdabc=sdabb,error=error)
                   DO i=1,3
                    CALL abc_contract(dabbint(sgfa:m1,sgfb:m2,sgfc:m3,i),sdabb(:,:,:,i),&
                         sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                         ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset),error)
                   ENDDO
                   DEALLOCATE(sdabb,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ELSE
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                                rab,dab,rac,dac,rbc,dbc,sabb,error=error)
                ENDIF
                ! debug if requested
                IF(debug) THEN
                   CALL overlap_abc_test(la_max(iset),npgfa(iset),zeta(:,iset),la_min(iset),&
                                         lb_max(jset),npgfb(jset),zetb(:,jset),lb_min(jset),&
                                         lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),lcb_min(kbset),&
                                         ra,rb,rb,sabb,dmax,error)
                ENDIF
                CALL abc_contract(abbint(sgfa:m1,sgfb:m2,sgfc:m3),sabb,&
                     sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset),error)
                DEALLOCATE(sabb,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
          END DO

       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE lri_int_abb

! *****************************************************************************
!> \brief calculate overlap integrals (a,b)
!> \param sab integral (a,b)
!> \param dsab derivative of sab with respect to A
!> \param ra ...
!> \param rb ...
!> \param rab ...
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_ab(sab,dsab,ra,rb,rab,fba,fbb,calculate_forces,debug,&
                        dmax,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dsab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: fba, fbb
    LOGICAL, INTENT(IN)                      :: calculate_forces, debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_ab', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iset, jset, lds, m1, m2, maxco, maxcoa, maxcob, &
      maxl, maxla, maxlb, ncoa, ncob, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sint
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: devab, swork
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb,&
            first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb,&
            sphi_a, sphi_b, zeta, zetb)

    ! basis ikind
    first_sgfa   =>  fba%first_sgf
    la_max       =>  fba%lmax
    la_min       =>  fba%lmin
    npgfa        =>  fba%npgf
    nseta        =   fba%nset
    nsgfa        =>  fba%nsgf_set
    rpgfa        =>  fba%pgf_radius
    set_radius_a =>  fba%set_radius
    sphi_a       =>  fba%sphi
    zeta         =>  fba%zet
    ! basis jkind
    first_sgfb   =>  fbb%first_sgf
    lb_max       =>  fbb%lmax
    lb_min       =>  fbb%lmin
    npgfb        =>  fbb%npgf
    nsetb        =   fbb%nset
    nsgfb        =>  fbb%nsgf_set
    rpgfb        =>  fbb%pgf_radius
    set_radius_b =>  fbb%set_radius
    sphi_b       =>  fbb%sphi
    zetb         =>  fbb%zet

    CALL get_gto_basis_set(fba,maxco=maxcoa,maxl=maxla)
    CALL get_gto_basis_set(fbb,maxco=maxcob,maxl=maxlb)
    maxco = MAX(maxcoa,maxcob)
    IF(calculate_forces) THEN
     maxl = MAX(maxla+1,maxlb)
    ELSE
     maxl = MAX(maxla,maxlb)
    ENDIF
    lds = ncoset(maxl)
    ALLOCATE(sint(maxco,maxco),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(calculate_forces) THEN
     !derivatives will be stored in devab(:,:,2:4)
     ALLOCATE(swork(lds,lds,4),devab(maxco,maxco,4),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     devab = 0._dp
    ELSE
     ALLOCATE(swork(lds,lds,1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    sint  = 0._dp
    swork = 0._dp

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1
          sint=0._dp
          swork = 0._dp

          ! calculate integrals
          IF(calculate_forces) THEN
             devab = 0._dp 
             CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,dab,sint,0,.FALSE.,swork,lds,sdab=devab)
             DO i=1,3
              !NOTE: devab(:,:,2:4) contains all derivatives for lmin=0 to lmax=lmax 
              !      correct after contraction (multiply with zero for elements l < lmin) 
              CALL ab_contract(dsab(sgfa:m1,sgfb:m2,i),devab(:,:,i+1),sphi_a(:,sgfa:),&
                   sphi_b(:,sgfb:),ncoa,ncob,nsgfa(iset),nsgfb(jset),error)
             ENDDO

          ELSE
             CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,dab,sint,0,.FALSE.,swork,lds)
          ENDIF
          ! debug if requested
          IF(debug) THEN
             CALL overlap_ab_test(la_max(iset),la_min(iset),npgfa(iset),zeta(:,iset),&
                                  lb_max(jset),lb_min(jset),npgfb(jset),zetb(:,jset),&
                                  ra,rb,sint,dmax,error)
          ENDIF

          CALL ab_contract(sab(sgfa:m1,sgfb:m2),sint,sphi_a(:,sgfa:),sphi_b(:,sgfb:),&
               ncoa,ncob,nsgfa(iset),nsgfb(jset),error)
       END DO
    END DO

    IF(calculate_forces) THEN
     DEALLOCATE(devab,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    DEALLOCATE(sint,swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_ab

! *****************************************************************************
!> \brief calculate overlap integrals (aa,bb)
!> \param saabb integral (aa,bb)
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param rab ...
!> \param ra ...
!> \param rb ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_aabb(saabb,oba,obb,rab,ra,rb,debug,dmax,error)

    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: saabb
    TYPE(gto_basis_set_type), POINTER        :: oba, obb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab, ra, rb
    LOGICAL, INTENT(IN)                      :: debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_aabb', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iset, isgfa1, jset, jsgfa2, kset, ksgfb1, lds, lset, &
      lsgfb2, m1, m2, m3, m4, maxco, maxcoa, maxcob, maxl, maxla, maxlb, &
      ncoa1, ncoa2, ncob1, ncob2, nseta, nsetb, sgfa1, sgfa2, sgfb1, sgfb2, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: asets_equal, bsets_equal, &
                                                failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: swork
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sint
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb,&
            first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb,&
            sphi_a, sphi_b, zeta, zetb)

    ! basis ikind
    first_sgfa   =>  oba%first_sgf
    la_max       =>  oba%lmax
    la_min       =>  oba%lmin
    npgfa        =>  oba%npgf
    nseta        =   oba%nset
    nsgfa        =>  oba%nsgf_set
    rpgfa        =>  oba%pgf_radius
    set_radius_a =>  oba%set_radius
    sphi_a       =>  oba%sphi
    zeta         =>  oba%zet
    ! basis jkind
    first_sgfb   =>  obb%first_sgf
    lb_max       =>  obb%lmax
    lb_min       =>  obb%lmin
    npgfb        =>  obb%npgf
    nsetb        =   obb%nset
    nsgfb        =>  obb%nsgf_set
    rpgfb        =>  obb%pgf_radius
    set_radius_b =>  obb%set_radius
    sphi_b       =>  obb%sphi
    zetb         =>  obb%zet

    CALL get_gto_basis_set(oba,maxco=maxcoa,maxl=maxla)
    CALL get_gto_basis_set(obb,maxco=maxcob,maxl=maxlb)
    maxco = MAX(maxcoa,maxcob)
    maxla = 2*maxla
    maxlb = 2*maxlb
    maxl = MAX(maxla,maxlb)
    lds = ncoset(maxl)
    ALLOCATE(sint(maxco,maxco,maxco,maxco),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(swork(lds,lds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sint  = 0._dp
    swork = 0._dp

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa1 = npgfa(iset)*ncoset(la_max(iset))
       sgfa1 = first_sgfa(1,iset)
       m1    = sgfa1 + nsgfa(iset)-1

       DO jset=iset,nseta

          ncoa2 = npgfa(jset)*ncoset(la_max(jset))
          sgfa2 = first_sgfa(1,jset)
          m2    = sgfa2 + nsgfa(jset)-1

          DO kset=1,nsetb

             ncob1 = npgfb(kset)*ncoset(lb_max(kset))
             sgfb1 = first_sgfb(1,kset)
             m3    = sgfb1 + nsgfb(kset)-1

             DO lset=kset,nsetb

                ncob2 = npgfb(lset)*ncoset(lb_max(lset))
                sgfb2 = first_sgfb(1,lset)
                m4=sgfb2+nsgfb(lset)-1

                ! check if sets are identical to spare some integral evaluation
                asets_equal = .FALSE.
                IF(iset == jset) asets_equal = .TRUE.
                bsets_equal = .FALSE.
                IF(kset == lset) bsets_equal = .TRUE.
                ! calculate integrals
                CALL overlap_aabb(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                                  la_max(jset),la_min(jset),npgfa(jset),rpgfa(:,jset),zeta(:,jset),&
                                  lb_max(kset),lb_min(kset),npgfb(kset),rpgfb(:,kset),zetb(:,kset),&
                                  lb_max(lset),lb_min(lset),npgfb(lset),rpgfb(:,lset),zetb(:,lset),&
                                  asets_equal,bsets_equal,rab,dab,sint,swork,lds)
                ! debug if requested
                IF (debug) THEN
                   CALL overlap_aabb_test(la_max(iset),la_min(iset),npgfa(iset),zeta(:,iset),&
                                          la_max(jset),la_min(jset),npgfa(jset),zeta(:,jset),&
                                          lb_max(kset),lb_min(kset),npgfb(kset),zetb(:,kset),&
                                          lb_max(lset),lb_min(lset),npgfb(lset),zetb(:,lset),&
                                          ra,rb,sint,dmax,error)
                ENDIF

                CALL abcd_contract(saabb(sgfa1:m1,sgfa2:m2,sgfb1:m3,sgfb2:m4),sint,sphi_a(:,sgfa1:),&
                                   sphi_a(:,sgfa2:),sphi_b(:,sgfb1:),sphi_b(:,sgfb2:),ncoa1,ncoa2,&
                                   ncob1,ncob2,nsgfa(iset),nsgfa(jset),nsgfb(kset),nsgfb(lset),error)

               ! account for the fact that some integrals are alike
                DO isgfa1 = sgfa1,m1
                  DO jsgfa2 = sgfa2,m2
                    DO ksgfb1 = sgfb1,m3
                      DO lsgfb2 = sgfb2,m4
                         saabb(jsgfa2,isgfa1,ksgfb1,lsgfb2) = saabb(isgfa1,jsgfa2,ksgfb1,lsgfb2)
                         saabb(isgfa1,jsgfa2,lsgfb2,ksgfb1) = saabb(isgfa1,jsgfa2,ksgfb1,lsgfb2)
                         saabb(jsgfa2,isgfa1,lsgfb2,ksgfb1) = saabb(isgfa1,jsgfa2,ksgfb1,lsgfb2)
                      END DO
                    END DO
                  END DO
                END DO
               
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE(sint,swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_aabb

! *****************************************************************************
!> \brief contract overlap integrals (a,b) and transfer to spherical Gaussians
!> \param abint ...
!> \param sab ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param ncoa ...
!> \param ncob ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ab_contract(abint,sab,sphi_a,sphi_b,ncoa,ncob,nsgfa,nsgfb,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: abint
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sab, sphi_a, sphi_b
    INTEGER, INTENT(IN)                      :: ncoa, ncob, nsgfa, nsgfb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ab_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m1, m2, msphia, msphib, nn, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cpp

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)

    m1 = SIZE(sab,1)
    m2 = SIZE(sab,2)

    nn = SIZE(abint,1)

    ALLOCATE(cpp(nsgfa,m2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL dgemm("T","N",nsgfa,m2,ncoa,1._dp,sphi_a,msphia,sab,m1,0.0_dp,cpp,nsgfa)
    CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,cpp,nsgfa,sphi_b,msphib,0.0_dp,&
               abint,nn)

    DEALLOCATE(cpp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE ab_contract

! *****************************************************************************
!> \brief contract three-center overlap integrals (a,b,c) and transfer
!>        to spherical Gaussians
!> \param abcint ...
!> \param sabc ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE abc_contract(abcint,sabc,sphi_a,sphi_b,sphi_c,ncoa,ncob,ncoc,&
                          nsgfa,nsgfb,nsgfc,error)

    REAL(KIND=dp), DIMENSION(:, :, :)        :: abcint, sabc
    REAL(KIND=dp), DIMENSION(:, :)           :: sphi_a, sphi_b, sphi_c
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, nsgfa, &
                                                nsgfb, nsgfc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'abc_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, isgfc, m1, m2, m3, &
                                                msphia, msphib, msphic, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: temp_ccc, work_cpc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: cpc, cpp

    CALL timeset(routineN,handle)

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)
    msphic = SIZE(sphi_c,1)

    m1 = SIZE(sabc,1)
    m2 = SIZE(sabc,2)
    m3 = SIZE(sabc,3)

    ALLOCATE(cpp(nsgfa,m2,m3),cpc(nsgfa,m2,nsgfc),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cpp = 0._dp
    cpc = 0._dp
    ALLOCATE(work_cpc(nsgfa,m2),temp_ccc(nsgfa,nsgfb),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    work_cpc(:,:) = 0._dp
    temp_ccc(:,:) = 0._dp

    CALL dgemm("T","N",nsgfa,m2*m3,ncoa,1._dp,sphi_a,msphia,sabc,m1,0.0_dp,cpp,nsgfa)
    CALL dgemm("N","N",nsgfa*m2,nsgfc,ncoc,1._dp,cpp,nsgfa*m2,sphi_c,msphic,0.0_dp,&
               cpc,nsgfa*m2)

    DO isgfc=1,nsgfc
     work_cpc(:,:)=cpc(:,:,isgfc)
     CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,work_cpc,nsgfa,sphi_b,msphib,&
                0.0_dp,temp_ccc,nsgfa)
     abcint(:,:,isgfc)=temp_ccc(:,:)
    END DO

    DEALLOCATE(cpp,cpc,work_cpc,temp_ccc,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE abc_contract

! *****************************************************************************
!> \brief contract four-center overlap integrals (a,b,c,d) and transfer
!>        to spherical Gaussians
!> \param abcdint ...
!> \param sabcd ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sphi_d ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE abcd_contract(abcdint,sabcd,sphi_a,sphi_b,sphi_c,sphi_d,ncoa,ncob,&
                           ncoc,ncod,nsgfa,nsgfb,nsgfc,nsgfd,error)

    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT)                          :: abcdint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sabcd
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sphi_a, sphi_b, sphi_c, sphi_d
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod, &
                                                nsgfa, nsgfb, nsgfc, nsgfd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'abcd_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, isgfc, isgfd, m1, m2, &
                                                m3, m4, msphia, msphib, &
                                                msphic, msphid, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: temp_cccc, work_cpcc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: temp_cpcc, work_cppc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: cpcc, cppc, cppp

    CALL timeset(routineN,handle)

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)
    msphic = SIZE(sphi_c,1)
    msphid = SIZE(sphi_d,1)

    m1 = SIZE(sabcd,1)
    m2 = SIZE(sabcd,2)
    m3 = SIZE(sabcd,3)
    m4 = SIZE(sabcd,4)

    ALLOCATE(cppp(nsgfa,m2,m3,m4),cppc(nsgfa,m2,m3,nsgfd),&
             cpcc(nsgfa,m2,nsgfc,nsgfd),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(work_cppc(nsgfa,m2,m3),temp_cpcc(nsgfa,m2,nsgfc),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    work_cppc = 0._dp
    temp_cpcc = 0._dp

    ALLOCATE(work_cpcc(nsgfa,m2),temp_cccc(nsgfa,nsgfb),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    work_cpcc = 0._dp
    temp_cccc = 0._dp

    CALL dgemm("T","N",nsgfa,m2*m3*m4,ncoa,1._dp,sphi_a,msphia,sabcd,m1,&
               0.0_dp,cppp,nsgfa)
    CALL dgemm("N","N",nsgfa*m2*m3,nsgfd,ncod,1._dp,cppp,nsgfa*m2*m3,&
               sphi_d,msphid,0.0_dp,cppc,nsgfa*m2*m3)

    DO isgfd=1,nsgfd
      work_cppc(:,:,:) = cppc(:,:,:,isgfd)
      CALL dgemm("N","N",nsgfa*m2,nsgfc,ncoc,1._dp,work_cppc,nsgfa*m2,&
                 sphi_c,msphic,0.0_dp,temp_cpcc,nsgfa*m2)
      cpcc(:,:,:,isgfd) = temp_cpcc(:,:,:)
      DO isgfc=1,nsgfc
        work_cpcc(:,:) = cpcc(:,:,isgfc,isgfd)
        CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,work_cpcc,nsgfa,sphi_b,&
                   msphib,0.0_dp,temp_cccc,nsgfa)
        abcdint(:,:,isgfc,isgfd) = temp_cccc(:,:) 
      END DO
    END DO

    DEALLOCATE(cpcc,cppc,cppp,stat=STAT)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(work_cpcc,work_cppc,temp_cpcc,temp_cccc,stat=STAT)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE abcd_contract

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
!> \param dab ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_calculate_derivative_acoef(lri_env,lrii,iatom,jatom,nba,nbb,&
                                            nfa,nfb,dab,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_int_type), POINTER              :: lrii
    INTEGER, INTENT(IN)                      :: iatom, jatom, nba, nbb, nfa, &
                                                nfb
    REAL(KIND=dp), INTENT(IN)                :: dab
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'lri_calculate_derivative_acoef', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nn, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: abxint_work, Q, R

    failure = .FALSE.
    CALL timeset(routineN,handle)

    nn = nfa + nfb
 
    ALLOCATE(R(nba,nbb,nn),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(Q(nba,nbb,nn),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(abxint_work(nba,nbb,nn),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    R(:,:,:) = 0._dp
    Q(:,:,:) = 0._dp
    abxint_work(:,:,:) = 0._dp

    abxint_work(1:nba,1:nbb,1:nfa)=lrii%abaint(1:nba,1:nbb,1:nfa)    !abaint
    abxint_work(1:nba,1:nbb,nfa+1:nn)=lrii%abbint(1:nba,1:nbb,1:nfb) !abbint
    IF(iatom == jatom.AND.dab < lri_env%delta) THEN
      DO i=1,nba
        DO j=1,nbb
          R(i,j,1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),abxint_work(i,j,1:nfa))
          Q(i,j,1:nfa) = lrii%soo(i,j)/lrii%nsn*lrii%sn(1:nfa) - &
                         DOT_PRODUCT(lrii%sn(1:nfa),abxint_work(i,j,1:nfa))/lrii%nsn&
                         *lrii%sn(1:nfa)
        ENDDO
      ENDDO
      lrii%dacoef(1:nba,1:nbb,1:nfa)= R(1:nba,1:nbb,1:nfa)+Q(1:nba,1:nbb,1:nfa)
    ELSE
      DO i=1,nba
        DO j=1,nbb
          R(i,j,1:nn) =  MATMUL(lrii%sinv(1:nn,1:nn),abxint_work(i,j,1:nn))
          Q(i,j,1:nn) =  lrii%soo(i,j)/lrii%nsn*lrii%sn(1:nn) - &
                         DOT_PRODUCT(lrii%sn(1:nn),abxint_work(i,j,1:nn))/lrii%nsn&
                         *lrii%sn(1:nn)
        ENDDO
      ENDDO
      lrii%dacoef(1:nba,1:nbb,1:nn)= R(1:nba,1:nbb,1:nn)+Q(1:nba,1:nbb,1:nn)
    ENDIF
    
    DEALLOCATE(abxint_work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(R,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(Q,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
  
  END SUBROUTINE lri_calculate_derivative_acoef

! *****************************************************************************
!> \brief debug output 
!> \param lri_env ...
!> \param qs_env ...
!> \param lri_ints ...
!> \param soo_list ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE output_debug_info(lri_env,qs_env,lri_ints,soo_list,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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
                    para_env=para_env,error=error)
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

    logger => cp_error_get_logger(error)
    iunit=cp_print_key_unit_nr(logger,input,"PRINT%PROGRAM_RUN_INFO",&
                               extension=".lridebug",error=error)

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
                    "PRINT%PROGRAM_RUN_INFO", error=error)
    CALL timestop(handle)
  
  END SUBROUTINE output_debug_info

! *****************************************************************************
!> \brief SVD of auxiliary overlap matrix
!> \param lri_env ...
!> \param qs_env ...
!> \param error ...
!> \note not called, routine for TESTING, to be deleted soon
! *****************************************************************************
  SUBROUTINE svd_of_s(lri_env,qs_env,error) 
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'svd_of_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iac, iatom, ikind, ilist, &
                                                info, jatom, jkind, &
                                                jneighbor, lwork, m, n, &
                                                nkind, nlist, nneighbor, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: sig, work
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: s_temp, u, vt
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    NULLIFY(soo_list, lrii, lri_ints, nl_iterator, s_temp, u , vt)
    soo_list => lri_env%soo_list
    lri_ints => lri_env%lri_ints
    CALL neighbor_list_iterator_create(nl_iterator,soo_list)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)

       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
            nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
            iatom=iatom,jatom=jatom,r=rab)
      
       iac = ikind + nkind*(jkind - 1)
       lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
      
       ! re-inverse sinv to get S 
       CALL invmat(lrii%sinv,stat,error)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       m=SIZE(lrii%sinv,1)
       n=m
       lwork = m*(6+4*m)+m 
       !mxm matrix
      
       ALLOCATE(s_temp(m,m),u(m,m),vt(m,m),sig(m),iwork(8*m),work(lwork))
       s_temp(:,:) = lrii%sinv 
      
       ! and change this back 
       CALL invmat(lrii%sinv,stat,error)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! do SVD
       CALL DGESDD('A',m,m,s_temp,m,sig,u,m,vt,m,work,lwork,iwork,info)
       WRITE(*,*) "iatom, jatom", iatom, jatom
       DO i=1,m
        WRITE(*,*) i, sig(i)
       ENDDO
       DEALLOCATE(s_temp, u, vt, sig, iwork, work)
    END DO

    CALL neighbor_list_iterator_release(nl_iterator)

  END SUBROUTINE

END MODULE lri_environment_methods
