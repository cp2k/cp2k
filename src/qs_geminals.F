!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate Geminal integrals
!> \par History
!>      07.2009 created 
!> \author J. Hutter
! *****************************************************************************
MODULE qs_geminals
  
  USE ai_geminals,                     ONLY: g2gemint,&
                                             gemint2,&
                                             gemint2_derivative
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: geminal_basis_set_type,&
                                             get_geminal_basis_set,&
                                             get_gto_basis_set,&
                                             gto_basis_set_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE ri_environment_types,            ONLY: get_ri_env,&
                                             ri_environment_type,&
                                             ri_vector_set,&
                                             ri_vector_type
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: geminal_coulomb, geminal_gto_coulomb
   
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_geminals'

! *****************************************************************************
  
  CONTAINS

! *****************************************************************************
!> \brief Computes two center geminal integrals for all atom pairs
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE geminal_coulomb(rixvec,rirvec,qs_env,energy,calculate_energy,calculate_force,calculate_virial,error)
    TYPE(ri_vector_type), POINTER            :: rixvec, rirvec
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(inout)             :: energy
    LOGICAL, INTENT(IN)                      :: calculate_energy, &
                                                calculate_force, &
                                                calculate_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geminal_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, ia, &
                                                iatom, ikind, istat, ja, &
                                                jatom, jkind, natom, nd, &
                                                nder, nga, ngb, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: ehfx1, fij
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ivec, jvec
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: intab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, rab, ri, rj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ivmat, jvmat
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_a, atomic_kind_b
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind_set,particle_set,local_particles)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    force=force,&
                    virial=virial,&
                    error=error)

    CALL ri_vector_set(rirvec, 0._dp, error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    IF ( calculate_force .OR. calculate_virial ) THEN
      nder = 1
      nd = 4
    ELSE
      nder = 0
      nd = 1
    END IF

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    ehfx1 = 0._dp

    DO ikind=1,nkind
       atomic_kind_a => atomic_kind_set(ikind)
       NULLIFY(geminal_basis_set)
       CALL get_atomic_kind(atomic_kind=atomic_kind_a,geminal_basis_set=geminal_basis_set)

       IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=nga)

       ALLOCATE (ivec(nga),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       ivmat => rixvec%vector(ikind)%vmat

       DO jkind=1,nkind
          atomic_kind_b => atomic_kind_set(jkind)
          NULLIFY(geminal_basis_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind_b,&
                          atom_list=atom_list,&
                          geminal_basis_set=geminal_basis_set)

          IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
          CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=ngb)

          ALLOCATE (intab(nga,ngb,nd),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jvec(ngb),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          jvmat => rixvec%vector(jkind)%vmat

          DO ia=1,local_particles%n_el(ikind)
             iatom=local_particles%list(ikind)%array(ia)
             ri = particle_set(iatom)%r
             DO ja=1,SIZE(atom_list)
                jatom=atom_list(ja)
                rj = particle_set(jatom)%r
                IF (iatom > jatom) THEN
                  IF (MODULO(iatom + jatom,2) == 0) CYCLE
                ELSE
                  IF (MODULO(iatom + jatom,2) /= 0) CYCLE
                END IF
 
                CALL gemgem_integral_atom(intab,atomic_kind_a,ri,atomic_kind_b,rj,nder,error)

                IF ( iatom==jatom) THEN
                  fij=1.0_dp
                ELSE
                  fij=2.0_dp
                END IF

                atom_a = atom_of_kind(iatom)
                atom_b = atom_of_kind(jatom)
                jvec(1:ngb) = fij*MATMUL(ivmat(1:nga,atom_a),intab(1:nga,1:ngb,1))
                ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,1),jvmat(1:ngb,atom_b))

                rirvec%vector(ikind)%vmat(1:nga,atom_a) = rirvec%vector(ikind)%vmat(1:nga,atom_a) + ivec(1:nga)
                rirvec%vector(jkind)%vmat(1:ngb,atom_b) = rirvec%vector(jkind)%vmat(1:ngb,atom_b) + jvec(1:ngb)

                IF ( calculate_energy ) THEN
                  ehfx1 = ehfx1 + DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                END IF
                IF ( calculate_force .OR. calculate_virial ) THEN
                   ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,2),jvmat(1:ngb,atom_b))
                   force_a(1) = DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                   ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,3),jvmat(1:ngb,atom_b))
                   force_a(2) = DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                   ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,4),jvmat(1:ngb,atom_b))
                   force_a(3) = DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                   force_b(1:3) = -force_a(1:3)
                   IF ( calculate_force ) THEN
                      force(ikind)%hfx_ri(:,atom_a)=force(ikind)%hfx_ri(:,atom_a) - force_a(:)
                      force(jkind)%hfx_ri(:,atom_b)=force(jkind)%hfx_ri(:,atom_b) - force_b(:)
                   END IF
                END IF
                IF ( calculate_virial ) THEN
                   rab = ri - rj
                   CALL virial_pair_force ( virial%pv_virial, 1.0_dp, force_a, rab, error)
                END IF

             END DO
          END DO

          DEALLOCATE (intab,jvec,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       END DO

       DEALLOCATE (ivec,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO

    DEALLOCATE (atom_of_kind,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    energy = ehfx1

    CALL timestop(handle)

  END SUBROUTINE geminal_coulomb

! *****************************************************************************
!> \brief Computes Coulomb integrals between geminals and products of 
!>        Gaussian orbitals
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE geminal_gto_coulomb(qs_env,calculate_fock,calculate_energy,calculate_force,calculate_virial,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_fock, &
                                                calculate_energy, &
                                                calculate_force, &
                                                calculate_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geminal_gto_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, atom_c, handle, i, iatom, ikind, istat, jatom, &
      jkind, k, ka, katom, kkind, maxlgto, na, natom, nb, nd, nder, ngem, &
      nkind, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_list, atom_of_kind, &
                                                kind_of
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: ehfx2, f0, fij
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: kdvec
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: intabc
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, force_c, &
                                                rac, rbc, ri, rj, rk
    REAL(KIND=dp), DIMENSION(:), POINTER     :: kcoeff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block1, ks_block2, &
                                                p_block1, p_block2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_a, atomic_kind_b, &
                                                atomic_kind_c
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set_a, &
                                                orb_basis_set_b
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_p
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff, ri_rhs
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind_set,particle_set,local_particles,rho,matrix_ks,ri_env,force,virial)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    rho=rho,&
                    matrix_ks=matrix_ks,&
                    hfx_ri_env=ri_env,&
                    force=force,&
                    virial=virial,&
                    error=error)
    CPPostcondition(ASSOCIATED(ri_env),cp_failure_level,routineP,error,failure)

    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,rhs=ri_rhs)
    CALL ri_vector_set(ri_rhs, 0._dp, error)

    matrix_p => rho%rho_ao
    nspins = SIZE(matrix_p)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    IF ( calculate_force .OR. calculate_virial ) THEN
      nder = 1
      nd   = 7
    ELSE
      nder = 0
      nd   = 1
    END IF

    ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             kind_of=kind_of,&
                             atom_of_kind=atom_of_kind)

    CALL init_orbital_pointers(4*maxlgto+nder)

    ehfx2 = 0._dp

    DO iatom=1,natom
       NULLIFY(block_node)
       block_node => first_block_node(matrix_p(1)%matrix,iatom)
       DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node,block_col=jatom,BLOCK=p_block1)
          IF ( calculate_fock ) THEN
            CALL get_block_node(matrix=matrix_ks(1)%matrix,block_row=iatom,block_col=jatom,&
                                BLOCK=ks_block1)
          END IF
          IF (nspins==2) THEN
            CALL get_block_node(matrix=matrix_p(2)%matrix,block_row=iatom,block_col=jatom,&
                                BLOCK=p_block2)
            IF ( calculate_fock ) THEN
              CALL get_block_node(matrix=matrix_ks(2)%matrix,block_row=iatom,block_col=jatom,&
                                  BLOCK=ks_block2)
            END IF
          END IF
          ikind  = kind_of(iatom)
          jkind  = kind_of(jatom)
          atom_a = atom_of_kind(iatom)
          atom_b = atom_of_kind(jatom)

          fij = 2._dp
          IF ( iatom==jatom ) fij=1._dp

          atomic_kind_a => atomic_kind_set(ikind)
          NULLIFY(orb_basis_set_a)
          CALL get_atomic_kind(atomic_kind=atomic_kind_a,orb_basis_set=orb_basis_set_a)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set_a,nsgf=na)

          atomic_kind_b => atomic_kind_set(jkind)
          NULLIFY(orb_basis_set_b)
          CALL get_atomic_kind(atomic_kind=atomic_kind_b,orb_basis_set=orb_basis_set_b)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set_b,nsgf=nb)

          ri = particle_set(iatom)%r
          rj = particle_set(jatom)%r

          DO kkind=1,nkind
             atomic_kind_c => atomic_kind_set(kkind)
             NULLIFY(geminal_basis_set)
             CALL get_atomic_kind(atomic_kind=atomic_kind_c,&
                                  atom_list=atom_list,&
                                  geminal_basis_set=geminal_basis_set)
             IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
             CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=ngem)

             ALLOCATE (intabc(na,nb,ngem,nd),kdvec(ngem,nd),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

             DO ka=1,SIZE(atom_list)
                katom=atom_list(ka)
                atom_c = atom_of_kind(katom)
                rk = particle_set(katom)%r

                CALL gemgau_integral_atom(intabc,atomic_kind_a,ri,atomic_kind_b,rj,&
                                          atomic_kind_c,rk,nder,error)

                DO i=1,ngem
                   kdvec(i,1) = SUM(intabc(1:na,1:nb,i,1)*p_block1(1:na,1:nb))
                END DO
                ri_rhs%vector(kkind)%vmat(1:ngem,atom_c) = ri_rhs%vector(kkind)%vmat(1:ngem,atom_c) &
                                                           + fij*kdvec(1:ngem,1)

                kcoeff => ri_coeff%vector(kkind)%vmat(:,atom_c)
                IF ( calculate_energy ) THEN
                  ehfx2 = ehfx2 + fij*DOT_PRODUCT(kdvec(1:ngem,1),kcoeff(1:ngem))
                END IF
                IF ( calculate_force .OR. calculate_virial ) THEN
                   DO k=2,7
                      DO i=1,ngem
                         kdvec(i,k) = fij*SUM(intabc(1:na,1:nb,i,k)*p_block1(1:na,1:nb))
                      END DO
                   END DO
                   force_c(1) = -DOT_PRODUCT(kdvec(1:ngem,2),kcoeff(1:ngem))
                   force_c(2) = -DOT_PRODUCT(kdvec(1:ngem,3),kcoeff(1:ngem))
                   force_c(3) = -DOT_PRODUCT(kdvec(1:ngem,4),kcoeff(1:ngem))
                   force_a(1) = -DOT_PRODUCT(kdvec(1:ngem,5),kcoeff(1:ngem))
                   force_a(2) = -DOT_PRODUCT(kdvec(1:ngem,6),kcoeff(1:ngem))
                   force_a(3) = -DOT_PRODUCT(kdvec(1:ngem,7),kcoeff(1:ngem))
                   force_b = -0.5_dp*(force_a + force_c)
                   force_a =  0.5_dp*(force_a - force_c)
                   IF ( calculate_force ) THEN
                      force(ikind)%hfx_ri(:,atom_a)=force(ikind)%hfx_ri(:,atom_a) - force_a(:)
                      force(jkind)%hfx_ri(:,atom_b)=force(jkind)%hfx_ri(:,atom_b) - force_b(:)
                      force(jkind)%hfx_ri(:,atom_c)=force(jkind)%hfx_ri(:,atom_c) - force_c(:)
                   END IF
                END IF
                IF ( calculate_virial ) THEN
                   f0 = 1.0_dp
                   rac = ri - rk
                   rbc = rj - rk
                   CALL virial_pair_force ( virial%pv_virial, f0, force_a, rac, error)
                   CALL virial_pair_force ( virial%pv_virial, f0, force_b, rbc, error)
                END IF
                IF ( calculate_fock ) THEN
                  DO i=1,ngem
                     ks_block1(1:na,1:nb) = ks_block1(1:na,1:nb)-0.5_dp*intabc(1:na,1:nb,i,1)*kcoeff(i)
                  END DO
                END IF

             END DO

             DEALLOCATE (intabc,kdvec,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          END DO
          block_node => next_block_node(block_node)
       END DO
    END DO

    DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ri_env%ehfx2 = -ehfx2

    CALL timestop(handle)

  END SUBROUTINE geminal_gto_coulomb

! *****************************************************************************
!> \brief Computes two center geminal integrals for full basis set
!> \param intab : list of integrals (output)
!> \param atomic_kind_a, atomic_kind_b : info on type of basis
!> \param ra, rb : positions of atoms                         
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE gemgem_integral_atom(intab,atomic_kind_a,ra,atomic_kind_b,rb,nderivative,error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: intab
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rb
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemgem_integral_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii, iira, iirb, iisa, iisb, ij, ipgf, ira, irb, isa, isb, &
      iset, ishell, istat, jj, jpgf, jset, jshell, lrmap, lsap, nder, ngema, &
      ngemb, nn, nna, nnb, nra, nrb, nsa, nsb, nseta, nsetb
    INTEGER, DIMENSION(:), POINTER           :: lrmaxa, lrmaxb, lrmina, &
                                                lrminb, lsa, lsb, npgfa, &
                                                npgfb, nshella, nshellb
    INTEGER, DIMENSION(:, :), POINTER        :: firsta, firstb, la, lb
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: cc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabd
    REAL(KIND=dp), DIMENSION(6)              :: rra, rrb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcca, gccb
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: zeta, zetb, zetbh
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set_a, &
                                                geminal_basis_set_b

    CALL get_atomic_kind(atomic_kind=atomic_kind_a,geminal_basis_set=geminal_basis_set_a)
    CALL get_atomic_kind(atomic_kind=atomic_kind_b,geminal_basis_set=geminal_basis_set_b)

    CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set_a,&
         nset=nseta,npgf=npgfa,lmax=lrmaxa,lmin=lrmina,ls=lsa,nshell=nshella,&
         l=la,first_cgf=firsta,ngeminals=ngema,gcc=gcca,zet=zeta)
    CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set_b,&
         nset=nsetb,npgf=npgfb,lmax=lrmaxb,lmin=lrminb,ls=lsb,nshell=nshellb,&
         l=lb,first_cgf=firstb,ngeminals=ngemb,gcc=gccb,zet=zetb,zeth=zetbh)

    CPPrecondition(nderivative>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nderivative<=1,cp_failure_level,routineP,error,failure)
    nder=ncoset(nderivative)

    CPPrecondition(SIZE(intab,1)>=ngema,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(intab,2)>=ngemb,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(intab,3)>=nder,cp_failure_level,routineP,error,failure)

    rra(1:3)=ra(1:3)
    rra(4:6)=0._dp
    rrb(1:3)=rb(1:3)
    rrb(4:6)=0._dp


    intab = 0._dp
    DO iset=1,nseta
       DO jset=1,nsetb
          lrmap = lrmaxa(iset)+nderivative
          lsap  = lsa(iset)+nderivative
          nn = npgfa(iset)*npgfb(jset)
          nra = ncoset(lrmap)
          nrb = ncoset(lrmaxb(jset))
          nsa = ncoset(lsap)
          nsb = ncoset(lsb(jset))
          ALLOCATE (iab(nn,nra,nsa,nrb,nsb),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          iab = 0._dp

          CALL gemint2(iab,zeta(:,:,:,iset),npgfa(iset),lrmap,lsap,rra,&
                       zetb(:,:,:,jset),npgfb(jset),lrmaxb(jset),lsb(jset),rrb,error)

          CALL gemint2(iab,zeta(:,:,:,iset),npgfa(iset),lrmap,lsap,rra,&
                       zetbh(:,:,:,jset),npgfb(jset),lrmaxb(jset),lsb(jset),rrb,error)

          ! derivatives
          IF ( nderivative > 0 ) THEN
            nra = ncoset(lrmaxa(iset))
            nsa = ncoset(lsa(iset))
            ALLOCATE (iabd(nn,nra,nsa,nrb,nsb,nder),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            CALL gemint2_derivative(iab,iabd,zeta(:,:,:,iset),npgfa(iset),lrmaxa(iset),lsa(iset),&
                 zetb(:,:,:,jset),npgfb(jset),lrmaxb(jset),lsb(jset),nderivative,.TRUE.,error)
          END IF

          ! contraction and storage of requested integrals
          ! we do both contractions and the reshuffeling of the integrals in one sweep
          ! assuming that the contraction depth is not large this should be efficient
          DO ishell=1,nshella(iset)
            nna = nco(lsa(iset))*nco(la(ishell,iset))
            DO jshell=1,nshellb(jset)
              nnb = nco(lsb(jset))*nco(lb(jshell,jset))
              DO ipgf=1,npgfa(iset)
                DO jpgf=1,npgfb(jset)
                  ij = (ipgf-1)*npgfb(jset)+jpgf
                  cc = gcca(ipgf,ishell,iset)*gccb(jpgf,jshell,jset)
                  DO irb=1,nco(lb(jshell,jset))
                    iirb = ncoset(lb(jshell,jset)-1) + irb
                    DO isb=1,nco(lsb(jset))
                      iisb = ncoset(lsb(jset)-1) + isb
                      jj = firstb(jshell,jset) + (irb-1)*nco(lsb(jset)) + isb - 1
                      IF ( nderivative > 0 ) THEN
                        DO ira=1,nco(la(ishell,iset))
                          iira = ncoset(la(ishell,iset)-1) + ira
                          DO isa=1,nco(lsa(iset))
                            iisa = ncoset(lsa(iset)-1) + isa
                            ii = firsta(ishell,iset) + (ira-1)*nco(lsa(iset)) + isa - 1
                            intab(ii,jj,1) = intab(ii,jj,1) + cc*iabd(ij,iira,iisa,iirb,iisb,1)
                            intab(ii,jj,2) = intab(ii,jj,2) + cc*iabd(ij,iira,iisa,iirb,iisb,2)
                            intab(ii,jj,3) = intab(ii,jj,3) + cc*iabd(ij,iira,iisa,iirb,iisb,3)
                            intab(ii,jj,4) = intab(ii,jj,4) + cc*iabd(ij,iira,iisa,iirb,iisb,4)
                          END DO
                        END DO
                      ELSE
                        DO ira=1,nco(la(ishell,iset))
                          iira = ncoset(la(ishell,iset)-1) + ira
                          DO isa=1,nco(lsa(iset))
                            iisa = ncoset(lsa(iset)-1) + isa
                            ii = firsta(ishell,iset) + (ira-1)*nco(lsa(iset)) + isa - 1
                            intab(ii,jj,1) = intab(ii,jj,1) + cc*iab(ij,iira,iisa,iirb,iisb)
                          END DO
                        END DO
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO

          DEALLOCATE (iab,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          IF ( nderivative > 0 ) THEN
            DEALLOCATE (iabd,STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

       END DO
    END DO

  END SUBROUTINE gemgem_integral_atom

! *****************************************************************************
!> \brief Computes three center geminal integrals for full basis set
!> \param intabc : list of integrals (output)
!> \param atomic_kind_a, atomic_kind_b, atomic_kind_c : info on type of basis
!> \param ra, rb, rc : positions of atoms                         
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE gemgau_integral_atom(intabc,atomic_kind_a,ra,atomic_kind_b,rb,&
                                  atomic_kind_c,rc,nderivative,error)
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: intabc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rb
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_c
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rc
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemgau_integral_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: ic, ii, iirc, iisc, il, ipgf, iq, irc, isc, iset, istat, jl, &
      jpgf, jq, jset, k, kpgf, kq, kset, kshell, ncoa, ncob, ncoc, nder, &
      ngemc, nn, nna, nnb, nrc, nsc, nseta, nsetb, nsetc, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, lrmaxc, lrminc, lsc, &
                                                npgfa, npgfb, npgfc, nsgfa, &
                                                nsgfb, nshellc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                firstc, lastc, lc
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: cc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: ab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: abc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabc
    REAL(KIND=dp), DIMENSION(6)              :: rrc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sphia, sphib, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gccc
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: zetc, zetch
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set_c
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set_a, &
                                                orb_basis_set_b

    CPPrecondition(nderivative>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nderivative<=1,cp_failure_level,routineP,error,failure)
    nder = 1
    IF ( nderivative==1 ) nder=7

    CALL get_atomic_kind(atomic_kind=atomic_kind_a,orb_basis_set=orb_basis_set_a)
    CALL get_atomic_kind(atomic_kind=atomic_kind_b,orb_basis_set=orb_basis_set_b)
    CALL get_atomic_kind(atomic_kind=atomic_kind_c,geminal_basis_set=geminal_basis_set_c)

    CALL get_gto_basis_set(gto_basis_set=orb_basis_set_a,nset=nseta,npgf=npgfa,&
         nsgf_set=nsgfa,lmax=la_max,lmin=la_min,sphi=sphia,first_sgf=first_sgfa,zet=zeta)
    CALL get_gto_basis_set(gto_basis_set=orb_basis_set_b,nset=nsetb,npgf=npgfb,&
         nsgf_set=nsgfb,lmax=lb_max,lmin=lb_min,sphi=sphib,first_sgf=first_sgfb,zet=zetb)
    CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set_c,npgf=npgfc,&
         nset=nsetc,lmax=lrmaxc,lmin=lrminc,ls=lsc,nshell=nshellc,&
         l=lc,first_cgf=firstc,last_cgf=lastc,ngeminals=ngemc,gcc=gccc,zet=zetc,zeth=zetch)

    rrc(1:3) = rc(1:3)
    rrc(4:6) = 0._dp

    intabc = 0._dp
    DO iset = 1,nseta
      sgfa = first_sgfa(1,iset)
      DO jset = 1,nsetb
        sgfb = first_sgfb(1,jset)
        DO kset = 1,nsetc
          nn = npgfa(iset)*npgfb(jset)*npgfc(kset)
          nna = ncoset(la_max(iset))
          nnb = ncoset(lb_max(jset))
          nrc = ncoset(lrmaxc(kset))
          nsc = ncoset(lsc(kset))
          ALLOCATE (iabc(nn,nna,nnb,nrc,nsc,nder),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          iabc = 0._dp

          CALL g2gemint(iabc,la_max(iset),npgfa(iset),zeta(:,iset),la_min(iset),ra,&
                        lb_max(jset),npgfb(jset),zetb(:,jset),lb_min(jset),rb,&
                        lrmaxc(kset),lsc(kset),npgfc(kset),zetc(:,:,:,kset),rrc,nderivative,error)
          CALL g2gemint(iabc,la_max(iset),npgfa(iset),zeta(:,iset),la_min(iset),ra,&
                        lb_max(jset),npgfb(jset),zetb(:,jset),lb_min(jset),rb,&
                        lrmaxc(kset),lsc(kset),npgfc(kset),zetch(:,:,:,kset),rrc,nderivative,error)

          ! contraction of geminal and reordering of integrals
          ncoa=nna*npgfa(iset)
          ncob=nnb*npgfb(jset)
          ncoc=lastc(nshellc(kset),kset)-firstc(1,kset)+1
          ALLOCATE (abc(ncoa,ncob,ncoc,nder),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          abc=0._dp
          DO kshell=1,nshellc(kset)
            DO kpgf=1,npgfc(kset)
              cc = gccc(kpgf,kshell,kset)
              DO irc=1,nco(lc(kshell,kset))
                iirc = ncoset(lc(kshell,kset)-1) + irc
                DO isc=1,nco(lsc(kset))
                  iisc = ncoset(lsc(kset)-1) + isc
                  kq=firstc(kshell,kset) + (irc-1)*nco(lsc(kset)) + isc - firstc(1,kset)
                  DO k=1,nder
                    DO il=ncoset(la_min(iset)-1)+1,ncoset(la_max(iset))
                      DO jl=ncoset(lb_min(jset)-1)+1,ncoset(lb_max(jset))
                        DO jpgf=1,npgfb(jset)
                          jq=(jpgf-1)*nnb+jl
                          DO ipgf=1,npgfa(iset)
                            ii=(kpgf-1)*npgfb(jset)*npgfa(iset)+(jpgf-1)*npgfa(iset)+ipgf
                            iq=(ipgf-1)*nna+il
                            abc(iq,jq,kq,k) = abc(iq,jq,kq,k) + cc*iabc(ii,il,jl,iirc,iisc,k)
                          END DO
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
          DEALLOCATE (iabc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (ab(ncoa,ncob),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO k=1,nder
            DO ic=1,ncoc
              ii = firstc(1,kset) + ic - 1

              CALL dgemm("T","N",nsgfa(iset),ncob,ncoa,1.0_dp,sphia(1,sgfa),SIZE(sphia,1),&
                         abc(1,1,ic,k),ncoa,0.0_dp,ab,ncoa)
              CALL dgemm("N","N",nsgfa(iset),nsgfb(iset),ncob,1.0_dp,ab,ncoa,sphib(1,sgfb),&
                         SIZE(sphib,1),0.0_dp,intabc(sgfa,sgfb,ii,k),SIZE(intabc,1))

            END DO
          END DO
          DEALLOCATE (abc,ab,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO
      END DO
    END DO

  END SUBROUTINE gemgau_integral_atom

! *****************************************************************************

END MODULE qs_geminals
