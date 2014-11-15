!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_iterator,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                             hfx_get_mult_cache_elements
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             pi
  USE mathlib,                         ONLY: invmat
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE ri_environment_types,            ONLY: get_ri_env,&
                                             ri_environment_type,&
                                             ri_vector_set,&
                                             ri_vector_sync,&
                                             ri_vector_type
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: geminal_coulomb, geminal_gto_coulomb, geminal_coulomb_diaginv, &
            geminal_charge

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_geminals'

! *****************************************************************************

  CONTAINS

! *****************************************************************************
!> \brief Computes two center geminal integrals for all atom pairs
!> \param rixvec ...
!> \param rirvec ...
!> \param qs_env ...
!> \param energy ...
!> \param calculate_energy ...
!> \param calculate_force ...
!> \param calculate_virial ...
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

    INTEGER :: atom_a, atom_b, handle, ia, iatom, ikind, istat, ja, jatom, &
      jkind, maxlgem(2), ml, natom, nd, nder, nga, ngb, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ehfx1, eps_screening, fij, &
                                                xfraction
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ivec, jvec
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: intab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, rab, ri, rj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ivmat, jvmat
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind_a, qs_kind_b
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(semi_empirical_si_type), POINTER    :: store_int
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind_set,qs_kind_set,particle_set,local_particles,ri_env)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    hfx_ri_env=ri_env,&
                    force=force,&
                    virial=virial,&
                    error=error)
    CPPostcondition(ASSOCIATED(ri_env),cp_failure_level,routineP,error,failure)

    NULLIFY(store_int)
    CALL get_ri_env(ri_env=ri_env,gem_integral_storage=store_int,screening=eps_screening,error=error)
    CALL get_ri_env(ri_env=ri_env,xfraction=xfraction,error=error)

    CALL ri_vector_set(rirvec, 0._dp, error)

    nkind = SIZE(qs_kind_set)
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
    CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
    CALL get_qs_kind_set(qs_kind_set, maxlgem=maxlgem)

    ml = MAXVAL(maxlgem)
    CALL init_orbital_pointers(2*ml+nder)

    ehfx1 = 0._dp

    DO ikind=1,nkind
       qs_kind_a => qs_kind_set(ikind)
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind=qs_kind_a,geminal_basis_set=geminal_basis_set)

       IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=nga)

       ALLOCATE (ivec(nga),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       ivmat => rixvec%vector(ikind)%vmat

       DO jkind=1,nkind
          qs_kind_b => qs_kind_set(jkind)
          NULLIFY(geminal_basis_set)
          CALL get_atomic_kind(atomic_kind_set(jkind), atom_list=atom_list)
          CALL get_qs_kind(qs_kind=qs_kind_b, geminal_basis_set=geminal_basis_set)

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

                CALL gemgem_integral_atom(nga,ngb,nd,intab,qs_kind_a,ri,qs_kind_b,rj,nder,&
                                          store_int,eps_screening,error)

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
                  ehfx1 = ehfx1 + DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a)) * xfraction
                END IF
                IF ( calculate_force .OR. calculate_virial ) THEN
                   ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,2),jvmat(1:ngb,atom_b))
                   force_a(1) = DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                   ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,3),jvmat(1:ngb,atom_b))
                   force_a(2) = DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                   ivec(1:nga) = fij*MATMUL(intab(1:nga,1:ngb,4),jvmat(1:ngb,atom_b))
                   force_a(3) = DOT_PRODUCT(ivec(1:nga),ivmat(1:nga,atom_a))
                   force_a(1:3) = force_a(1:3) * xfraction
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

    CALL get_qs_env(qs_env=qs_env,para_env=para_env,error=error)
    CALL ri_vector_sync(rirvec, para_env, error)
    CALL mp_sum(ehfx1,para_env%group)

    energy = ehfx1

    CALL timestop(handle)

  END SUBROUTINE geminal_coulomb

! *****************************************************************************
!> \brief Computes Coulomb integrals between geminals and products of
!>        Gaussian orbitals
!> \param qs_env ...
!> \param calculate_fock ...
!> \param calculate_energy ...
!> \param calculate_force ...
!> \param calculate_virial ...
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

    INTEGER :: atom_a, atom_b, atom_c, blk, handle, i, iatom, ikind, istat, &
      jatom, jkind, k, ka, katom, kkind, maxlgem(2), maxlgto, ml, na, natom, &
      nb, nd, nder, ngem, nkind, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_list, atom_of_kind, &
                                                kind_of
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: ehfx2, eps_screening, f0, &
                                                fij, xfraction
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
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set_a, &
                                                orb_basis_set_b
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind_a, qs_kind_b, &
                                                qs_kind_c
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff, ri_rhs
    TYPE(semi_empirical_si_type), POINTER    :: store_int
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(qs_kind_set,particle_set,local_particles,rho,matrix_ks,ri_env,force,virial)
    NULLIFY(matrix_ks,matrix_p)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    rho=rho,&
                    matrix_ks=matrix_ks,&
                    hfx_ri_env=ri_env,&
                    force=force,&
                    virial=virial,&
                    error=error)
    CPPostcondition(ASSOCIATED(ri_env),cp_failure_level,routineP,error,failure)

    NULLIFY(store_int)
    CALL get_ri_env(ri_env=ri_env,ggg_integral_storage=store_int,screening=eps_screening,error=error)
    CALL get_ri_env(ri_env=ri_env,xfraction=xfraction,error=error)

    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,rhs=ri_rhs,error=error)
    CALL ri_vector_set(ri_rhs, 0._dp, error)

    matrix_p => rho%rho_ao
    nspins = SIZE(matrix_p)

    nkind = SIZE(qs_kind_set)
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
    CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)
    CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, maxlgem=maxlgem)

    ml = MAXVAL(maxlgem)
    CALL init_orbital_pointers(2*maxlgto+ml+nder)

    ehfx2 = 0._dp

    CALL cp_dbcsr_iterator_start(iter, matrix_p(1)%matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, p_block1, blk)
       IF ( calculate_fock ) THEN
          CALL cp_dbcsr_get_block_p(matrix=matrix_ks(1)%matrix,&
               row=iatom,col=jatom,BLOCK=ks_block1,found=found)
       END IF
       IF (nspins==2) THEN
          CALL cp_dbcsr_get_block_p(matrix=matrix_p(2)%matrix,&
               row=iatom,col=jatom,BLOCK=p_block2,found=found)
          IF ( calculate_fock ) THEN
             CALL cp_dbcsr_get_block_p(matrix=matrix_ks(2)%matrix,&
                  row=iatom,col=jatom,BLOCK=ks_block2,found=found)
          ENDIF
       ENDIF

       ikind  = kind_of(iatom)
       jkind  = kind_of(jatom)
       atom_a = atom_of_kind(iatom)
       atom_b = atom_of_kind(jatom)

       fij = 2._dp
       IF ( iatom==jatom ) fij=1._dp

       qs_kind_a => qs_kind_set(ikind)
       NULLIFY(orb_basis_set_a)
       CALL get_qs_kind(qs_kind_a,orb_basis_set=orb_basis_set_a)
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set_a,nsgf=na)

       qs_kind_b => qs_kind_set(jkind)
       NULLIFY(orb_basis_set_b)
       CALL get_qs_kind(qs_kind_b,orb_basis_set=orb_basis_set_b)
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set_b,nsgf=nb)

       ri = particle_set(iatom)%r
       rj = particle_set(jatom)%r

       DO kkind=1,nkind
          qs_kind_c => qs_kind_set(kkind)
          NULLIFY(geminal_basis_set)
          CALL get_atomic_kind(atomic_kind_set(kkind), atom_list=atom_list)
          CALL get_qs_kind(qs_kind_c, geminal_basis_set=geminal_basis_set)
          IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
          CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=ngem)

          ALLOCATE (intabc(na,nb,ngem,nd),kdvec(ngem,nd),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          DO ka=1,SIZE(atom_list)
             katom=atom_list(ka)
             atom_c = atom_of_kind(katom)
             rk = particle_set(katom)%r

             CALL gemgau_integral_atom(na,nb,ngem,nd,intabc,qs_kind_a,ri,qs_kind_b,rj,&
                  qs_kind_c,rk,nder,store_int,eps_screening,error)

             DO i=1,ngem
                kdvec(i,1) = SUM(intabc(1:na,1:nb,i,1)*p_block1(1:na,1:nb))
             END DO
             ri_rhs%vector(kkind)%vmat(1:ngem,atom_c) = ri_rhs%vector(kkind)%vmat(1:ngem,atom_c) &
                                                           + fij*kdvec(1:ngem,1)

             kcoeff => ri_coeff%vector(kkind)%vmat(:,atom_c)
             IF ( calculate_energy ) THEN
                ehfx2 = ehfx2 + fij*DOT_PRODUCT(kdvec(1:ngem,1),kcoeff(1:ngem)) * xfraction
             END IF
             IF ( calculate_force .OR. calculate_virial ) THEN
                DO k=2,7
                   DO i=1,ngem
                      kdvec(i,k) = xfraction*fij*SUM(intabc(1:na,1:nb,i,k)*p_block1(1:na,1:nb))
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
                   force(kkind)%hfx_ri(:,atom_c)=force(kkind)%hfx_ri(:,atom_c) - force_c(:)
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
                   ks_block1(1:na,1:nb) = ks_block1(1:na,1:nb)-0.5_dp*intabc(1:na,1:nb,i,1)*kcoeff(i)*xfraction
                END DO
             END IF

          END DO

          DEALLOCATE (intabc,kdvec,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       END DO
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,para_env=para_env,error=error)
    CALL ri_vector_sync(ri_rhs, para_env, error)
    CALL mp_sum(ehfx2,para_env%group)

    ri_env%ehfx2 = -ehfx2

    CALL timestop(handle)

  END SUBROUTINE geminal_gto_coulomb

! *****************************************************************************
!> \brief Computes two center geminal integrals for all diagonal atom blocks
!>        and applys its inverse to a vector
!> \param rixvec ...
!> \param rirvec ...
!> \param qs_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE geminal_coulomb_diaginv(rixvec,rirvec,qs_env,error)
    TYPE(ri_vector_type), POINTER            :: rixvec, rirvec
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geminal_coulomb_diaginv', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, cache_size, handle, iat, iatom, ikind, info, istat, &
      mem_compression_counter, memory_usage, natom, new_size, nga, nints, &
      nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    LOGICAL                                  :: buffer_overflow, failure
    REAL(KIND=dp)                            :: eps_screening, eps_storage
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ivec
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: intab
    REAL(KIND=dp), DIMENSION(3)              :: ri
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ivmat
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind_a
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(semi_empirical_si_type), POINTER    :: store_int

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind_set,qs_kind_set,particle_set,local_particles,ri_env)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    hfx_ri_env=ri_env,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    error=error)
    CPPostcondition(ASSOCIATED(ri_env),cp_failure_level,routineP,error,failure)

    CALL get_ri_env(ri_env=ri_env,gemdiag_integral_storage=store_int,screening=eps_screening,error=error)

    CALL ri_vector_set(rirvec, 0._dp, error)

    nkind = SIZE(qs_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)

    DO ikind=1,nkind
       qs_kind_a => qs_kind_set(ikind)
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind=qs_kind_a,geminal_basis_set=geminal_basis_set)

       IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=nga)

       ALLOCATE (ivec(nga),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       ivmat => rixvec%vector(ikind)%vmat

       ALLOCATE (intab(nga,nga,1),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DO iat=1,local_particles%n_el(ikind)
          iatom=local_particles%list(ikind)%array(iat)
          ri = particle_set(iatom)%r

          atom_a = atom_of_kind(iatom)
          ivec(1:nga) = ivmat(1:nga,atom_a)

          IF (store_int%memory_parameter%do_all_on_the_fly) THEN
             CALL gemgem_integral_atom_calc(intab,qs_kind_a,ri,qs_kind_a,ri,0,error)
             CALL invmat(intab(1:nga,1:nga,1),info,error)
             CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
             ivec(1:nga) = MATMUL(intab(1:nga,1:nga,1),ivec(1:nga))
          ELSE
             CPPostcondition(.NOT.store_int%memory_parameter%do_disk_storage,cp_failure_level,routineP,error,failure)
             nints      = nga*nga
             cache_size = store_int%memory_parameter%cache_size
             eps_storage= store_int%memory_parameter%eps_storage_scaling * eps_screening
             eps_storage= MAX(eps_storage,1.e-15_dp)
             IF (store_int%filling_containers) THEN
                mem_compression_counter = store_int%memory_parameter%actual_memory_usage * cache_size
                IF(mem_compression_counter > store_int%memory_parameter%max_compression_counter) THEN
                   buffer_overflow = .TRUE.
                   store_int%memory_parameter%ram_counter = store_int%nbuffer
                ELSE
                   store_int%nbuffer = store_int%nbuffer+1
                   buffer_overflow = .FALSE.
                END IF

                CALL gemgem_integral_atom_calc(intab,qs_kind_a,ri,qs_kind_a,ri,0,error)
                CALL invmat(intab(1:nga,1:nga,1),info,error)
                CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

                ! Store integrals if we did not go overflow
                IF (.NOT.buffer_overflow) THEN

                   memory_usage = store_int%memory_parameter%actual_memory_usage
                   CPPrecondition((nints/1.2_dp)<=HUGE(0)-memory_usage,cp_failure_level,routineP,error,failure)
                   IF (memory_usage+nints>SIZE(store_int%uncompressed_container)) THEN
                      new_size = INT((memory_usage+nints)*1.2_dp)
                      CALL reallocate(store_int%uncompressed_container, 1, new_size)
                   END IF
                   CALL dcopy(nints,intab,1,store_int%uncompressed_container(memory_usage),1)
                   store_int%memory_parameter%actual_memory_usage = memory_usage + nints

                END IF

                ivec(1:nga) = MATMUL(intab(1:nga,1:nga,1),ivec(1:nga))

             ELSE
                ! Get integrals from the containers
                IF(store_int%memory_parameter%ram_counter == store_int%nbuffer) THEN
                   buffer_overflow = .TRUE.
                ELSE
                   store_int%nbuffer = store_int%nbuffer + 1
                   buffer_overflow = .FALSE.
                END IF
                ! Get integrals from cache unless we overflowed
                IF (.NOT.buffer_overflow) THEN

                   memory_usage = store_int%memory_parameter%actual_memory_usage
                   CALL dcopy(nints,store_int%uncompressed_container(memory_usage),1,intab,1)
                   store_int%memory_parameter%actual_memory_usage = memory_usage + nints

                ELSE

                   CALL gemgem_integral_atom_calc(intab,qs_kind_a,ri,qs_kind_a,ri,0,error)
                   CALL invmat(intab(1:nga,1:nga,1),info,error)
                   CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

                END IF

                ivec(1:nga) = MATMUL(intab(1:nga,1:nga,1),ivec(1:nga))

             END IF

          END IF

          rirvec%vector(ikind)%vmat(1:nga,atom_a) = rirvec%vector(ikind)%vmat(1:nga,atom_a) + 2._dp*ivec(1:nga)

       END DO

       DEALLOCATE (intab,ivec,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO

    DEALLOCATE (atom_of_kind,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,para_env=para_env,error=error)
    CALL ri_vector_sync(rirvec, para_env, error)

    CALL timestop(handle)

  END SUBROUTINE geminal_coulomb_diaginv

! *****************************************************************************
!> \brief Computes the total charge from the geminal expansion of the density matrix
!> \param rivec ...
!> \param qs_env ...
!> \param charge ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE geminal_charge(rivec,qs_env,charge,error)
    TYPE(ri_vector_type), POINTER            :: rivec
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(out)               :: charge
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geminal_charge', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ii, ikind, ipgf, ira, &
                                                iset, ishell, istat, ix, iy, &
                                                iz, ll, natom, nga, nkind, &
                                                nseta
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: lrmaxa, lsa, npgfa, nshella
    INTEGER, DIMENSION(:, :), POINTER        :: firsta, la
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: al, cc, ff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ivec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ivmat
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcca
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(ri_environment_type), POINTER       :: ri_env

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(qs_kind_set,particle_set,local_particles,ri_env)

    CALL get_qs_env(qs_env=qs_env,&
                    qs_kind_set=qs_kind_set,&
                    atomic_kind_set=atomic_kind_set,&
                    hfx_ri_env=ri_env,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    error=error)
    CPPostcondition(ASSOCIATED(ri_env),cp_failure_level,routineP,error,failure)

    CALL get_ri_env(ri_env=ri_env,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)

    charge = 0._dp

    DO ikind=1,nkind
       qs_kind => qs_kind_set(ikind)
       NULLIFY(geminal_basis_set)
       CALL get_qs_kind(qs_kind=qs_kind,geminal_basis_set=geminal_basis_set)

       IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,nset=nseta,npgf=npgfa,&
            lmax=lrmaxa,ls=lsa,nshell=nshella,l=la,first_cgf=firsta,ngeminals=nga,gcc=gcca,zet=zeta)

       ALLOCATE (ivec(nga),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ivec = 0._dp

       DO iset=1,nseta
          IF (lsa(iset) /= 0) CYCLE
          DO ishell=1,nshella(iset)
            IF (MOD(la(ishell,iset),2) /= 0) CYCLE
            ll = la(ishell,iset)
            DO ipgf=1,npgfa(iset)
              al = zeta(1,1,ipgf,iset)
              cc = gcca(ipgf,ishell,iset)
              DO ira=1,nco(ll)
                ix = indco(1,ira+ncoset(ll-1))
                iy = indco(2,ira+ncoset(ll-1))
                iz = indco(3,ira+ncoset(ll-1))
                IF ( MOD(ix,2)/=0 .OR. MOD(iy,2)/=0 .OR. MOD(iz,2)/=0 ) CYCLE
                ff = dfac(ix/2)*dfac(iy/2)*dfac(iz/2)
                ii = firsta(ishell,iset) + (ira-1)*nco(lsa(iset))
                ivec(ii) = ivec(ii) + cc*2._dp*pi**1.5_dp/2._dp**(ll/2)/al**(ll/2+1.5_dp)*ff
              END DO
            END DO
          END DO
       END DO

       ivmat => rivec%vector(ikind)%vmat
       charge = charge + SUM(MATMUL(ivec(1:nga),ivmat(1:nga,:)))

       DEALLOCATE (ivec,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO

    DEALLOCATE (atom_of_kind,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE geminal_charge

! *****************************************************************************
!> \brief Computes two center geminal integrals for full basis set
!> \param nga ...
!> \param ngb ...
!> \param nd ...
!> \param intab : list of integrals (output)
!> \param qs_kind_a ...
!> \param ra positions of atoms
!> \param qs_kind_b ...
!> \param rb positions of atoms
!> \param nderivative ...
!> \param store_int ...
!> \param eps_screening ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE gemgem_integral_atom(nga,ngb,nd,intab,qs_kind_a,ra,qs_kind_b,rb,nderivative,&
                                  store_int,eps_screening,error)
    INTEGER, INTENT(IN)                      :: nga, ngb, nd
    REAL(KIND=dp), DIMENSION(nga, ngb, nd), &
      INTENT(OUT)                            :: intab
    TYPE(qs_kind_type), POINTER              :: qs_kind_a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra
    TYPE(qs_kind_type), POINTER              :: qs_kind_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rb
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(semi_empirical_si_type), POINTER    :: store_int
    REAL(KIND=dp), INTENT(IN)                :: eps_screening
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemgem_integral_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: bb, buffer_left, buffer_size, buffer_start, cache_size, ia, &
      ib, mem_compression_counter, nbits, new_size, nints
    LOGICAL                                  :: buffer_overflow, failure
    REAL(KIND=dp)                            :: eps_storage

    failure = .FALSE.
    IF (store_int%memory_parameter%do_all_on_the_fly .OR. nderivative > 0) THEN
       CALL gemgem_integral_atom_calc(intab,qs_kind_a,ra,qs_kind_b,rb,nderivative,error)
    ELSE
       CPPostcondition(store_int%compress,cp_failure_level,routineP,error,failure)
       CPPostcondition(.NOT.store_int%memory_parameter%do_disk_storage,cp_failure_level,routineP,error,failure)
       nints      = nga*ngb
       cache_size = store_int%memory_parameter%cache_size
       eps_storage= store_int%memory_parameter%eps_storage_scaling * eps_screening
       eps_storage= MAX(eps_storage,1.e-15_dp)
       IF (store_int%filling_containers) THEN
          mem_compression_counter = store_int%memory_parameter%actual_memory_usage * cache_size
          IF(mem_compression_counter > store_int%memory_parameter%max_compression_counter) THEN
             buffer_overflow = .TRUE.
             store_int%memory_parameter%ram_counter = store_int%nbuffer
          ELSE
             store_int%nbuffer = store_int%nbuffer+1
             buffer_overflow = .FALSE.
          END IF

          CALL gemgem_integral_atom_calc(intab,qs_kind_a,ra,qs_kind_b,rb,nderivative,error)

          ! Store integrals if we did not go overflow
          IF (.NOT.buffer_overflow) THEN
             ! Store integrals in the containers
             IF (store_int%nbuffer>SIZE(store_int%max_val_buffer)) THEN
                new_size = store_int%nbuffer+1000
                CALL reallocate(store_int%max_val_buffer, 1, new_size)
             END IF
             store_int%max_val_buffer(store_int%nbuffer) = MAXVAL(ABS(intab(:,:,1)))
             nbits        = EXPONENT(store_int%max_val_buffer(store_int%nbuffer)/eps_storage) + 1
             nbits        = MIN(nbits,64)
             nbits        = MAX(nbits,1)
             buffer_left  = nints
             buffer_start = 1
             DO WHILE (buffer_left > 0)
                buffer_size = MIN(buffer_left, cache_size)
                bb=buffer_start-1
                ib=bb/nga+1
                ia=bb-(ib-1)*nga+1
                CALL hfx_add_mult_cache_elements(intab(ia,ib,1),&
                     buffer_size,nbits,store_int%integral_caches(nbits), &
                     store_int%integral_containers(nbits),eps_storage,1.0_dp,&
                     store_int%memory_parameter%actual_memory_usage,.FALSE.,error)
                buffer_left  = buffer_left  - buffer_size
                buffer_start = buffer_start + buffer_size
             END DO
          END IF

       ELSE
          ! Get integrals from the containers
          IF(store_int%memory_parameter%ram_counter == store_int%nbuffer) THEN
             buffer_overflow = .TRUE.
          ELSE
             store_int%nbuffer = store_int%nbuffer + 1
             buffer_overflow = .FALSE.
          END IF
          ! Get integrals from cache unless we overflowed
          IF (.NOT.buffer_overflow) THEN
             ! Get Integrals from containers
             nbits        = EXPONENT(store_int%max_val_buffer(store_int%nbuffer)/eps_storage) + 1
             nbits        = MIN(nbits,64)
             nbits        = MAX(nbits,1)
             buffer_left  = nints
             buffer_start = 1
             DO WHILE (buffer_left > 0)
                buffer_size = MIN(buffer_left, cache_size)
                bb=buffer_start-1
                ib=bb/nga+1
                ia=bb-(ib-1)*nga+1
                CALL hfx_get_mult_cache_elements(intab(ia,ib,1),&
                     buffer_size,nbits,store_int%integral_caches(nbits), &
                     store_int%integral_containers(nbits),eps_storage,1.0_dp, &
                     store_int%memory_parameter%actual_memory_usage,.FALSE.)
                buffer_left  = buffer_left  - buffer_size
                buffer_start = buffer_start + buffer_size
             END DO
          ELSE
             CALL gemgem_integral_atom_calc(intab,qs_kind_a,ra,qs_kind_b,rb,nderivative,error)
          END IF

       END IF
    END IF

  END SUBROUTINE gemgem_integral_atom
! *****************************************************************************
!> \brief ...
!> \param intab ...
!> \param qs_kind_a ...
!> \param ra ...
!> \param qs_kind_b ...
!> \param rb ...
!> \param nderivative ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE gemgem_integral_atom_calc(intab,qs_kind_a,ra,qs_kind_b,rb,nderivative,error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: intab
    TYPE(qs_kind_type), POINTER              :: qs_kind_a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra
    TYPE(qs_kind_type), POINTER              :: qs_kind_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rb
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemgem_integral_atom_calc', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii, iira, iirb, iisa, iisb, ij, ipgf, ira, irb, isa, isb, &
      iset, ishell, istat, jj, jpgf, jset, jshell, lrmap, lsap, ltap, ltb, &
      nder, ngema, ngemb, nn, nna, nnb, nra, nrb, nsa, nsb, nseta, nsetb
    INTEGER, DIMENSION(:), POINTER           :: lrmaxa, lrmaxb, lrmina, &
                                                lrminb, lsa, lsb, npgfa, &
                                                npgfb, nshella, nshellb
    INTEGER, DIMENSION(:, :), POINTER        :: firsta, firstb, la, lb
    LOGICAL                                  :: failure
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

    failure = .FALSE.
    CALL get_qs_kind(qs_kind=qs_kind_a,geminal_basis_set=geminal_basis_set_a)
    CALL get_qs_kind(qs_kind=qs_kind_b,geminal_basis_set=geminal_basis_set_b)

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
          ltap  = lrmap+lsap-nderivative
          nn = npgfa(iset)*npgfb(jset)
          nra = ncoset(lrmap)
          nrb = ncoset(lrmaxb(jset))
          nsa = ncoset(lsap)
          nsb = ncoset(lsb(jset))
          ltb = lrmaxb(jset)+lsb(jset)
          ALLOCATE (iab(nn,nra,nsa,nrb,nsb),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          iab = 0._dp

          CALL gemint2(iab,zeta(:,:,:,iset),npgfa(iset),lrmap,lsap,ltap,rra,&
                       zetb(:,:,:,jset),npgfb(jset),lrmaxb(jset),lsb(jset),ltb,rrb,error)

          CALL gemint2(iab,zeta(:,:,:,iset),npgfa(iset),lrmap,lsap,ltap,rra,&
                       zetbh(:,:,:,jset),npgfb(jset),lrmaxb(jset),lsb(jset),ltb,rrb,error)

          ! derivatives
          IF ( nderivative > 0 ) THEN
            nra = ncoset(lrmaxa(iset))
            nsa = ncoset(lsa(iset))
            ALLOCATE (iabd(nn,nra,nsa,nrb,nsb,nder),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            CALL gemint2_derivative(iab,iabd,zeta(:,:,:,iset),npgfa(iset),lrmaxa(iset),lsa(iset),&
                 npgfb(jset),lrmaxb(jset),lsb(jset),nderivative,.TRUE.,error)
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

  END SUBROUTINE gemgem_integral_atom_calc

! *****************************************************************************
!> \brief Computes three center geminal integrals for full basis set
!> \param na ...
!> \param nb ...
!> \param ngem ...
!> \param nd ...
!> \param intabc : list of integrals (output)
!> \param qs_kind_a ...
!> \param ra positions of atoms
!> \param qs_kind_b ...
!> \param rb positions of atoms
!> \param qs_kind_c ...
!> \param rc positions of atoms
!> \param nderivative ...
!> \param store_int ...
!> \param eps_screening ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!> \author J. Hutter
! *****************************************************************************
  SUBROUTINE gemgau_integral_atom(na,nb,ngem,nd,intabc,qs_kind_a,ra,qs_kind_b,rb,&
             qs_kind_c,rc,nderivative,store_int,eps_screening,error)
    INTEGER, INTENT(IN)                      :: na, nb, ngem, nd
    REAL(KIND=dp), &
      DIMENSION(na, nb, ngem, nd), &
      INTENT(OUT)                            :: intabc
    TYPE(qs_kind_type), POINTER              :: qs_kind_a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra
    TYPE(qs_kind_type), POINTER              :: qs_kind_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rb
    TYPE(qs_kind_type), POINTER              :: qs_kind_c
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rc
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(semi_empirical_si_type), POINTER    :: store_int
    REAL(KIND=dp), INTENT(IN)                :: eps_screening
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemgau_integral_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: bb, buffer_left, buffer_size, buffer_start, cache_size, ia, &
      ib, ig, mem_compression_counter, nbits, new_size, nints
    LOGICAL                                  :: buffer_overflow, failure
    REAL(KIND=dp)                            :: eps_storage

    failure = .FALSE.
    IF (store_int%memory_parameter%do_all_on_the_fly .OR. nderivative > 0) THEN
       CALL gemgau_integral_atom_calc(intabc,qs_kind_a,ra,qs_kind_b,rb,&
                                      qs_kind_c,rc,nderivative,error)
    ELSE
       CPPostcondition(store_int%compress,cp_failure_level,routineP,error,failure)
       CPPostcondition(.NOT.store_int%memory_parameter%do_disk_storage,cp_failure_level,routineP,error,failure)
       nints      = na*nb*ngem
       cache_size = store_int%memory_parameter%cache_size
       eps_storage= store_int%memory_parameter%eps_storage_scaling * eps_screening
       eps_storage= MAX(eps_storage,1.e-15_dp)
       IF (store_int%filling_containers) THEN
          mem_compression_counter = store_int%memory_parameter%actual_memory_usage * cache_size
          IF(mem_compression_counter > store_int%memory_parameter%max_compression_counter) THEN
             buffer_overflow = .TRUE.
             store_int%memory_parameter%ram_counter = store_int%nbuffer
          ELSE
             store_int%nbuffer = store_int%nbuffer+1
             buffer_overflow = .FALSE.
          END IF

          CALL gemgau_integral_atom_calc(intabc,qs_kind_a,ra,qs_kind_b,rb,&
                                         qs_kind_c,rc,nderivative,error)

          ! Store integrals if we did not go overflow
          IF (.NOT.buffer_overflow) THEN
             ! Store integrals in the containers
             IF (store_int%nbuffer>SIZE(store_int%max_val_buffer)) THEN
                new_size = store_int%nbuffer+1000
                CALL reallocate(store_int%max_val_buffer, 1, new_size)
             END IF
             store_int%max_val_buffer(store_int%nbuffer) = MAXVAL(ABS(intabc(:,:,:,1)))
             nbits        = EXPONENT(store_int%max_val_buffer(store_int%nbuffer)/eps_storage) + 1
             nbits        = MIN(nbits,64)
             nbits        = MAX(nbits,1)
             buffer_left  = nints
             buffer_start = 1
             DO WHILE (buffer_left > 0)
                buffer_size = MIN(buffer_left, cache_size)
                bb=buffer_start-1
                ig=bb/(na*nb)+1
                bb=bb-(ig-1)*na*nb
                ib=bb/na+1
                ia=bb-(ib-1)*na+1
                CALL hfx_add_mult_cache_elements(intabc(ia,ib,ig,1),&
                     buffer_size,nbits,store_int%integral_caches(nbits), &
                     store_int%integral_containers(nbits),eps_storage,1.0_dp,&
                     store_int%memory_parameter%actual_memory_usage,.FALSE.,error)
                buffer_left  = buffer_left  - buffer_size
                buffer_start = buffer_start + buffer_size
             END DO
          END IF

       ELSE
          ! Get integrals from the containers
          IF(store_int%memory_parameter%ram_counter == store_int%nbuffer) THEN
             buffer_overflow = .TRUE.
          ELSE
             store_int%nbuffer = store_int%nbuffer + 1
             buffer_overflow = .FALSE.
          END IF
          ! Get integrals from cache unless we overflowed
          IF (.NOT.buffer_overflow) THEN
             ! Get Integrals from containers
             nbits        = EXPONENT(store_int%max_val_buffer(store_int%nbuffer)/eps_storage) + 1
             nbits        = MIN(nbits,64)
             nbits        = MAX(nbits,1)
             buffer_left  = nints
             buffer_start = 1
             DO WHILE (buffer_left > 0)
                buffer_size = MIN(buffer_left, cache_size)
                bb=buffer_start-1
                ig=bb/(na*nb)+1
                bb=bb-(ig-1)*na*nb
                ib=bb/na+1
                ia=bb-(ib-1)*na+1
                CALL hfx_get_mult_cache_elements(intabc(ia,ib,ig,1),&
                     buffer_size,nbits,store_int%integral_caches(nbits), &
                     store_int%integral_containers(nbits),eps_storage,1.0_dp, &
                     store_int%memory_parameter%actual_memory_usage,.FALSE.)
                buffer_left  = buffer_left  - buffer_size
                buffer_start = buffer_start + buffer_size
             END DO
          ELSE
             CALL gemgau_integral_atom_calc(intabc,qs_kind_a,ra,qs_kind_b,rb,&
                                            qs_kind_c,rc,nderivative,error)
          END IF

       END IF
    END IF

  END SUBROUTINE gemgau_integral_atom

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param intabc ...
!> \param qs_kind_a ...
!> \param ra ...
!> \param qs_kind_b ...
!> \param rb ...
!> \param qs_kind_c ...
!> \param rc ...
!> \param nderivative ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE gemgau_integral_atom_calc(intabc,qs_kind_a,ra,qs_kind_b,rb,&
                                  qs_kind_c,rc,nderivative,error)
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: intabc
    TYPE(qs_kind_type), POINTER              :: qs_kind_a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra
    TYPE(qs_kind_type), POINTER              :: qs_kind_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rb
    TYPE(qs_kind_type), POINTER              :: qs_kind_c
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rc
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemgau_integral_atom_calc', &
      routineP = moduleN//':'//routineN

    INTEGER :: ic, ii, iirc, iisc, il, ipgf, iq, irc, isc, iset, istat, jl, &
      jpgf, jq, jset, k, kpgf, kq, kset, kshell, na, nb, ncoa, ncob, ncoc, &
      ncosa, ncosb, nder, ngemc, nn, nna, nnb, nrc, nsc, nseta, nsetb, nsetc, &
      sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, lrmaxc, lrminc, lsc, &
                                                npgfa, npgfb, npgfc, nsgfa, &
                                                nsgfb, nshellc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                firstc, lastc, lc
    LOGICAL                                  :: failure
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

    failure = .FALSE.
    CPPrecondition(nderivative>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nderivative<=1,cp_failure_level,routineP,error,failure)
    nder = 1
    IF ( nderivative==1 ) nder=7

    CALL get_qs_kind(qs_kind=qs_kind_a,orb_basis_set=orb_basis_set_a)
    CALL get_qs_kind(qs_kind=qs_kind_b,orb_basis_set=orb_basis_set_b)
    CALL get_qs_kind(qs_kind=qs_kind_c,geminal_basis_set=geminal_basis_set_c)

    CALL get_gto_basis_set(gto_basis_set=orb_basis_set_a,nset=nseta,npgf=npgfa,&
         nsgf=na,nco_sum=ncosa,&
         nsgf_set=nsgfa,lmax=la_max,lmin=la_min,sphi=sphia,first_sgf=first_sgfa,zet=zeta)
    CALL get_gto_basis_set(gto_basis_set=orb_basis_set_b,nset=nsetb,npgf=npgfb,&
         nsgf=nb,nco_sum=ncosb,&
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

          CALL g2gemint(iabc,la_max(iset),npgfa(iset),zeta(:,iset),ra,&
                        lb_max(jset),npgfb(jset),zetb(:,jset),rb,&
                        lrmaxc(kset),lsc(kset),npgfc(kset),zetc(:,:,:,kset),rrc,nderivative,error)
          CALL g2gemint(iabc,la_max(iset),npgfa(iset),zeta(:,iset),ra,&
                        lb_max(jset),npgfb(jset),zetb(:,jset),rb,&
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
              CALL dgemm("N","N",nsgfa(iset),nsgfb(jset),ncob,1.0_dp,ab,ncoa,sphib(1,sgfb),&
                         SIZE(sphib,1),0.0_dp,intabc(sgfa,sgfb,ii,k),SIZE(intabc,1))

            END DO
          END DO
          DEALLOCATE (abc,ab,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO
      END DO
    END DO

  END SUBROUTINE gemgau_integral_atom_calc

! *****************************************************************************

END MODULE qs_geminals
