!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE atom_grb
   USE ai_onecenter,                    ONLY: sg_conf,&
                                              sg_kinetic,&
                                              sg_nuclear,&
                                              sg_overlap
   USE atom_electronic_structure,       ONLY: calculate_atom
   USE atom_operators,                  ONLY: atom_int_release,&
                                              atom_int_setup,&
                                              atom_ppint_release,&
                                              atom_ppint_setup,&
                                              atom_relint_release,&
                                              atom_relint_setup
   USE atom_types,                      ONLY: &
        CGTO_BASIS, GTO_BASIS, atom_basis_type, atom_integrals, atom_orbitals, atom_p_type, &
        atom_potential_type, atom_state, atom_type, create_atom_orbs, create_atom_type, &
        release_atom_basis, release_atom_type, set_atom
   USE atom_utils,                      ONLY: atom_basis_condnum,&
                                              atom_density
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE input_constants,                 ONLY: barrier_conf,&
                                              do_analytic,&
                                              do_rhf_atom,&
                                              do_rks_atom,&
                                              do_rohf_atom,&
                                              do_uhf_atom,&
                                              do_uks_atom
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE lapack,                          ONLY: lapack_ssygv
   USE mathconstants,                   ONLY: dfac,&
                                              pi
   USE orbital_pointers,                ONLY: deallocate_orbital_pointers,&
                                              init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: deallocate_spherical_harmonics,&
                                              init_spherical_harmonics
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: bohr
   USE powell,                          ONLY: opt_state_type,&
                                              powell_optimize
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom
#include "./base/base_uses.f90"

   IMPLICIT NONE

   TYPE basis_p_type
      TYPE(atom_basis_type), POINTER                :: basis
   END TYPE basis_p_type

   PRIVATE
   PUBLIC  :: atom_grb_construction

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_grb'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************

!> \brief ...
!> \param atom_info ...
!> \param atom_section ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE atom_grb_construction(atom_info, atom_section, iw)

      TYPE(atom_p_type), DIMENSION(:, :), POINTER        :: atom_info
      TYPE(section_vals_type), POINTER                   :: atom_section
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_grb_construction', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=100)                                 :: abas, basname
      CHARACTER(len=100), DIMENSION(1)                   :: basline
      CHARACTER(len=100), DIMENSION(20)                  :: headline
      INTEGER                                            :: i, ider, is, iunit, j, k, l, lhomo, ll, &
                                                            lval, m, maxl, mb, method, mo, n, &
                                                            nder, ngp, nhomo, nr, num_gto, &
                                                            num_pol, quadtype, s1, s2
      INTEGER, DIMENSION(0:3)                            :: next_bas, next_prim
      INTEGER, DIMENSION(0:7)                            :: nbas
      INTEGER, DIMENSION(:), POINTER                     :: num_bas
      REAL(KIND=dp) :: al, amin, aval, cnum, crad, cradx, cval, delta, dene, ear, emax, &
         energy_ex(0:3), energy_ref, energy_vb(0:3), expzet, fhomo, o, prefac, rconf, rk, rmax, &
         scon, zeta, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ale, alp, rho
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ebasis, pbasis, qbasis, rbasis
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: wfn
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: ovlp
      TYPE(atom_basis_type), POINTER                     :: basis, basis_grb, basis_ref, basis_vrb
      TYPE(atom_integrals), POINTER                      :: atint
      TYPE(atom_orbitals), POINTER                       :: orbitals
      TYPE(atom_state), POINTER                          :: state
      TYPE(atom_type), POINTER                           :: atom, atom_ref, atom_test
      TYPE(basis_p_type), DIMENSION(0:10)                :: vbasis
      TYPE(section_vals_type), POINTER                   :: grb_section, powell_section

      IF (iw > 0) WRITE (iw, '(/," ",79("*"),/,T28,A,/," ",79("*"))') "GEOMETRICAL RESPONSE BASIS"

      DO i = 0, 10
         NULLIFY (vbasis(i)%basis)
      END DO
      ! make some basic checks
      is = SIZE(atom_info)
      IF (iw > 0 .AND. is > 1) THEN
         WRITE (iw, '(/,A,/)') " WARNING: Only use first electronic structure/method for basis set generation"
      END IF
      atom_ref => atom_info(1, 1)%atom

      ! check method
      method = atom_ref%method_type
      SELECT CASE (method)
      CASE (do_rks_atom, do_rhf_atom)
         ! restricted methods are okay
      CASE (do_uks_atom, do_uhf_atom, do_rohf_atom)
         CPABORT("Unrestricted methods not allowed for GRB generation")
      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! input for basis optimization
      grb_section => section_vals_get_subs_vals(atom_section, "PRINT%GEOMETRICAL_RESPONSE_BASIS")

      ! generate an atom type
      NULLIFY (atom)
      CALL create_atom_type(atom)
      CALL copy_atom_basics(atom_ref, atom, state=.TRUE., potential=.TRUE., optimization=.TRUE., xc=.TRUE.)
      ! set confinement potential
      atom%potential%confinement = .TRUE.
      atom%potential%conf_type = barrier_conf
      atom%potential%acon = 200._dp
      atom%potential%rcon = 4._dp
      CALL section_vals_val_get(grb_section, "CONFINEMENT", r_val=scon)
      atom%potential%scon = scon
      ! generate main block geometrical exponents
      basis_ref => atom_ref%basis
      ALLOCATE (basis)
      NULLIFY (basis%am, basis%cm, basis%as, basis%ns, basis%bf, basis%dbf, basis%ddbf)
      ! get information on quadrature type and number of grid points
      ! allocate and initialize the atomic grid
      NULLIFY (basis%grid)
      CALL allocate_grid_atom(basis%grid)
      CALL section_vals_val_get(grb_section, "QUADRATURE", i_val=quadtype)
      CALL section_vals_val_get(grb_section, "GRID_POINTS", i_val=ngp)
      IF (ngp <= 0) &
         CPABORT("# point radial grid < 0")
      CALL create_grid_atom(basis%grid, ngp, 1, 1, 0, quadtype)
      basis%grid%nr = ngp
      !
      maxl = atom%state%maxl_occ
      basis%basis_type = GTO_BASIS
      CALL section_vals_val_get(grb_section, "NUM_GTO_CORE", i_val=num_gto)
      basis%nbas = 0
      basis%nbas(0:maxl) = num_gto
      basis%nprim = basis%nbas
      CALL section_vals_val_get(grb_section, "GEOMETRICAL_FACTOR", r_val=cval)
      CALL section_vals_val_get(grb_section, "GEO_START_VALUE", r_val=aval)
      m = MAXVAL(basis%nbas)
      ALLOCATE (basis%am(m, 0:3))
      basis%am = 0._dp
      DO l = 0, 3
         DO i = 1, basis%nbas(l)
            ll = i-1
            basis%am(i, l) = aval*cval**(ll)
         END DO
      END DO

      basis%eps_eig = basis_ref%eps_eig
      basis%geometrical = .TRUE.
      basis%aval = aval
      basis%cval = cval
      basis%start = 0

      ! initialize basis function on a radial grid
      nr = basis%grid%nr
      m = MAXVAL(basis%nbas)
      ALLOCATE (basis%bf(nr, m, 0:3))
      ALLOCATE (basis%dbf(nr, m, 0:3))
      ALLOCATE (basis%ddbf(nr, m, 0:3))
      basis%bf = 0._dp
      basis%dbf = 0._dp
      basis%ddbf = 0._dp
      DO l = 0, 3
         DO i = 1, basis%nbas(l)
            al = basis%am(i, l)
            DO k = 1, nr
               rk = basis%grid%rad(k)
               ear = EXP(-al*basis%grid%rad(k)**2)
               basis%bf(k, i, l) = rk**l*ear
               basis%dbf(k, i, l) = (REAL(l, dp)*rk**(l-1)-2._dp*al*rk**(l+1))*ear
               basis%ddbf(k, i, l) = (REAL(l*(l-1), dp)*rk**(l-2)- &
                                      2._dp*al*REAL(2*l+1, dp)*rk**(l)+4._dp*al*rk**(l+2))*ear
            END DO
         END DO
      END DO

      NULLIFY (orbitals)
      mo = MAXVAL(atom%state%maxn_calc)
      mb = MAXVAL(basis%nbas)
      CALL create_atom_orbs(orbitals, mb, mo)
      CALL set_atom(atom, orbitals=orbitals)

      powell_section => section_vals_get_subs_vals(atom_section, "POWELL")
      CALL atom_fit_grb(atom, basis, iw, powell_section)
      CALL set_atom(atom, basis=basis)

      ! generate response contractions
      CALL section_vals_val_get(grb_section, "DELTA_CHARGE", r_val=delta)
      CALL section_vals_val_get(grb_section, "DERIVATIVES", i_val=nder)
      IF (iw > 0) THEN
         WRITE (iw, '(/,A,T76,I5)') " Generate Response Basis Sets with Order ", nder
      END IF

      state => atom%state
      ! find HOMO
      lhomo = -1
      nhomo = -1
      emax = -HUGE(1._dp)
      DO l = 0, state%maxl_occ
         DO i = 1, state%maxn_occ(l)
            IF (atom%orbitals%ener(i, l) > emax) THEN
               lhomo = l
               nhomo = i
               emax = atom%orbitals%ener(i, l)
               fhomo = state%occupation(l, i)
            END IF
         END DO
      END DO

      s1 = SIZE(atom%orbitals%wfn, 1)
      s2 = SIZE(atom%orbitals%wfn, 2)
      ALLOCATE (wfn(s1, s2, 0:3, -nder:nder))
      s2 = MAXVAL(state%maxn_occ)+nder
      ALLOCATE (rbasis(s1, s2, 0:3), qbasis(s1, s2, 0:3))
      rbasis = 0._dp
      qbasis = 0._dp

      ! calculate integrals
      ALLOCATE (atint)
      CALL atom_int_setup(atint, basis, potential=atom%potential, eri_coulomb=.FALSE., eri_exchange=.FALSE.)
      CALL atom_ppint_setup(atint, basis, potential=atom%potential)
      IF (atom%pp_calc) THEN
         NULLIFY (atint%tzora, atint%hdkh)
      ELSE
         ! relativistic correction terms
         CALL atom_relint_setup(atint, basis, atom%relativistic, zcore=REAL(atom%z, dp))
      END IF
      CALL set_atom(atom, integrals=atint)

      CALL calculate_atom(atom, iw=0)
      DO ider = -nder, nder
         dene = REAL(ider, KIND=dp)*delta
         CPASSERT(fhomo > ABS(dene))
         state%occupation(lhomo, nhomo) = fhomo+dene
         CALL calculate_atom(atom, iw=0, noguess=.TRUE.)
         wfn(:, :, :, ider) = atom%orbitals%wfn
         state%occupation(lhomo, nhomo) = fhomo
      END DO
      IF (iw > 0) THEN
         WRITE (iw, '(A,T76,I5)') " Total number of electronic structure calculations ", 2*nder+1
      END IF

      ovlp => atom%integrals%ovlp

      DO l = 0, state%maxl_occ
         IF (iw > 0) THEN
            WRITE (iw, '(A,T76,I5)') " Response derivatives for l quantum number ", l
         END IF
         ! occupied states
         DO i = 1, MAX(state%maxn_occ(l), 1)
            rbasis(:, i, l) = wfn(:, i, l, 0)
         END DO
         ! differentiation
         DO ider = 1, nder
            i = MAX(state%maxn_occ(l), 1)
            SELECT CASE (ider)
            CASE (1)
               rbasis(:, i+1, l) = 0.5_dp*(wfn(:, i, l, 1)-wfn(:, i, l, -1))/delta
            CASE (2)
               rbasis(:, i+2, l) = 0.25_dp*(wfn(:, i, l, 2)-2._dp*wfn(:, i, l, 0)+wfn(:, i, l, -2))/delta**2
            CASE (3)
               rbasis(:, i+3, l) = 0.125_dp*(wfn(:, i, l, 3)-3._dp*wfn(:, i, l, 1) &
                                             +3._dp*wfn(:, i, l, -1)-wfn(:, i, l, -3))/delta**3
            CASE DEFAULT
               CPABORT("")
            END SELECT
         END DO

         ! orthogonalization, use gram-schmidt in order to keep the natural order (semi-core, valence, response) of the wfn.
         n = state%maxn_occ(l)+nder
         m = atom%basis%nbas(l)
         DO i = 1, n
            DO j = 1, i-1
               o = DOT_PRODUCT(rbasis(1:m, j, l), RESHAPE(MATMUL(ovlp(1:m, 1:m, l), rbasis(1:m, i:i, l)), (/m/)))
               rbasis(1:m, i, l) = rbasis(1:m, i, l)-o*rbasis(1:m, j, l)
            ENDDO
            o = DOT_PRODUCT(rbasis(1:m, i, l), RESHAPE(MATMUL(ovlp(1:m, 1:m, l), rbasis(1:m, i:i, l)), (/m/)))
            rbasis(1:m, i, l) = rbasis(1:m, i, l)/SQRT(o)
         ENDDO

         ! check
         ALLOCATE (amat(n, n))
         amat(1:n, 1:n) = MATMUL(TRANSPOSE(rbasis(1:m, 1:n, l)), MATMUL(ovlp(1:m, 1:m, l), rbasis(1:m, 1:n, l)))
         DO i = 1, n
            amat(i, i) = amat(i, i)-1._dp
         END DO
         IF (MAXVAL(ABS(amat)) > 1.e-12) THEN
            IF (iw > 0) WRITE (iw, '(A,G20.10)') " Orthogonality error  ", MAXVAL(ABS(amat))
         END IF
         DEALLOCATE (amat)

         ! Quickstep normalization
         expzet = 0.25_dp*REAL(2*l+3, dp)
         prefac = SQRT(SQRT(pi)/2._dp**(l+2)*dfac(2*l+1))
         DO i = 1, m
            zeta = (2._dp*atom%basis%am(i, l))**expzet
            qbasis(i, 1:n, l) = rbasis(i, 1:n, l)*prefac/zeta
         END DO

      END DO

      ! check for condition numbers
      IF (iw > 0) WRITE (iw, '(/,A)') " Condition Number of Valence Response Basis Sets"
      CALL init_orbital_pointers(3)
      CALL init_spherical_harmonics(3, 0)
      DO ider = 0, nder
         NULLIFY (basis_vrb)
         ALLOCATE (basis_vrb)
         NULLIFY (basis_vrb%am, basis_vrb%cm, basis_vrb%as, basis_vrb%ns, basis_vrb%bf, &
                  basis_vrb%dbf, basis_vrb%ddbf)
         ! allocate and initialize the atomic grid
         NULLIFY (basis_vrb%grid)
         CALL allocate_grid_atom(basis_vrb%grid)
         CALL create_grid_atom(basis_vrb%grid, ngp, 1, 1, 0, quadtype)
         basis_vrb%grid%nr = ngp
         !
         basis_vrb%eps_eig = basis_ref%eps_eig
         basis_vrb%geometrical = .FALSE.
         basis_vrb%basis_type = CGTO_BASIS
         basis_vrb%nprim = basis%nprim
         basis_vrb%nbas = 0
         DO l = 0, state%maxl_occ
            basis_vrb%nbas(l) = state%maxn_occ(l)+ider
         END DO
         m = MAXVAL(basis_vrb%nprim)
         n = MAXVAL(basis_vrb%nbas)
         ALLOCATE (basis_vrb%am(m, 0:3))
         basis_vrb%am = basis%am
         ! contractions
         ALLOCATE (basis_vrb%cm(m, n, 0:3))
         DO l = 0, state%maxl_occ
            m = basis_vrb%nprim(l)
            n = basis_vrb%nbas(l)
            basis_vrb%cm(1:m, 1:n, l) = rbasis(1:m, 1:n, l)
         END DO

         ! initialize basis function on a radial grid
         nr = basis_vrb%grid%nr
         m = MAXVAL(basis_vrb%nbas)
         ALLOCATE (basis_vrb%bf(nr, m, 0:3))
         ALLOCATE (basis_vrb%dbf(nr, m, 0:3))
         ALLOCATE (basis_vrb%ddbf(nr, m, 0:3))
         basis_vrb%bf = 0._dp
         basis_vrb%dbf = 0._dp
         basis_vrb%ddbf = 0._dp
         DO l = 0, 3
            DO i = 1, basis_vrb%nprim(l)
               al = basis_vrb%am(i, l)
               DO k = 1, nr
                  rk = basis_vrb%grid%rad(k)
                  ear = EXP(-al*basis_vrb%grid%rad(k)**2)
                  DO j = 1, basis_vrb%nbas(l)
                     basis_vrb%bf(k, j, l) = basis_vrb%bf(k, j, l)+rk**l*ear*basis_vrb%cm(i, j, l)
                     basis_vrb%dbf(k, j, l) = basis_vrb%dbf(k, j, l) &
                                              +(REAL(l, dp)*rk**(l-1)-2._dp*al*rk**(l+1))*ear*basis_vrb%cm(i, j, l)
                     basis_vrb%ddbf(k, j, l) = basis_vrb%ddbf(k, j, l)+ &
                                               (REAL(l*(l-1), dp)*rk**(l-2)-2._dp*al*REAL(2*l+1, dp)*rk**(l)+ &
                                                4._dp*al*rk**(l+2))*ear*basis_vrb%cm(i, j, l)
                  END DO
               END DO
            END DO
         END DO

         IF (iw > 0) THEN
            CALL basis_label(abas, basis_vrb%nprim, basis_vrb%nbas)
            WRITE (iw, '(A,A)') " Basis set     ", TRIM(abas)
         END IF
         crad = 2.0_dp*ptable(atom%z)%covalent_radius*bohr
         cradx = crad*1.00_dp
         CALL atom_basis_condnum(basis_vrb, cradx, cnum)
         IF (iw > 0) WRITE (iw, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
         cradx = crad*1.10_dp
         CALL atom_basis_condnum(basis_vrb, cradx, cnum)
         IF (iw > 0) WRITE (iw, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
         cradx = crad*1.20_dp
         CALL atom_basis_condnum(basis_vrb, cradx, cnum)
         IF (iw > 0) WRITE (iw, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
         vbasis(ider)%basis => basis_vrb
      END DO
      CALL deallocate_orbital_pointers
      CALL deallocate_spherical_harmonics

      ! generate polarization sets
      IF (iw > 0) THEN
         WRITE (iw, '(/,A)') " Polarization basis set  "
      END IF
      maxl = atom%state%maxl_occ
      CALL section_vals_val_get(grb_section, "NUM_GTO_POLARIZATION", i_val=num_gto)
      CPASSERT(num_gto > 0)
      num_pol = num_gto
      ALLOCATE (pbasis(num_gto, num_gto, 0:7), alp(num_gto))
      pbasis = 0.0_dp
      ! get density maximum
      ALLOCATE (rho(basis%grid%nr))
      CALL calculate_atom(atom, iw=0, noguess=.TRUE.)
      CALL atom_density(rho(:), atom%orbitals%pmat, atom%basis, maxl, typ="RHO")
      n = SUM(MAXLOC(rho(:)))
      rmax = basis%grid%rad(n)
      DEALLOCATE (rho)
      ! optimize exponents
      lval = maxl+1
      zval = SQRT(REAL(2*lval+2, dp))*REAL(lval+1, dp)/(2._dp*rmax)
      aval = atom%basis%am(1, 0)
      cval = 2.5_dp
      rconf = atom%potential%scon
      CALL atom_fit_pol(zval, rconf, lval, aval, cval, num_gto, iw, powell_section)
      ! calculate contractions
      DO i = 1, num_gto
         alp(i) = aval*cval**(i-1)
      END DO
      ALLOCATE (rho(num_gto))
      DO l = maxl+1, MIN(maxl+num_gto, 7)
         zval = SQRT(REAL(2*l+2, dp))*REAL(l+1, dp)/(2._dp*rmax)
         CALL hydrogenic(zval, rconf, l, alp, num_gto, rho, pbasis(:, :, l))
         IF (iw > 0) WRITE (iw, '(T5,A,i5,T66,A,F10.4)') &
            " Polarization basis set contraction for lval=", l, "zval=", zval
      END DO
      DEALLOCATE (rho)

      ! generate valence expansion sets
      maxl = atom%state%maxl_occ
      CALL section_vals_val_get(grb_section, "NUM_GTO_EXTENDED", i_val=num_gto)
      CALL section_vals_val_get(grb_section, "EXTENSION_BASIS", i_vals=num_bas)
      next_bas(0:3) = 0
      IF (num_bas(1) == -1) THEN
         DO l = 0, maxl
            next_bas(l) = maxl-l+1
         END DO
      ELSE
         n = MIN(SIZE(num_bas, 1), 4)
         next_bas(0:n-1) = num_bas(1:n)
      END IF
      next_prim = 0
      DO l = 0, 3
         IF (next_bas(l) > 0) next_prim(l) = num_gto
      END DO
      IF (iw > 0) THEN
         CALL basis_label(abas, next_prim, next_bas)
         WRITE (iw, '(/,A,A)') " Extension basis set     ", TRIM(abas)
      END IF
      n = MAXVAL(next_prim)
      m = MAXVAL(next_bas)
      ALLOCATE (ebasis(n, n, 0:3), ale(n))
      basis_vrb => vbasis(0)%basis
      amin = atom%basis%aval/atom%basis%cval**1.5_dp
      DO i = 1, n
         ale(i) = amin*atom%basis%cval**(i-1)
      END DO
      ebasis = 0._dp
      ALLOCATE (rho(n))
      rconf = 2.0_dp*atom%potential%scon
      DO l = 0, 3
         IF (next_bas(l) < 1) CYCLE
         zval = SQRT(REAL(2*l+2, dp))*REAL(l+1, dp)/(2._dp*rmax)
         CALL hydrogenic(zval, rconf, l, ale, n, rho, ebasis(:, :, l))
         IF (iw > 0) WRITE (iw, '(T5,A,i5,T66,A,F10.4)') &
            " Extension basis set contraction for lval=", l, "zval=", zval
      END DO
      DEALLOCATE (rho)
      ! check for condition numbers
      IF (iw > 0) WRITE (iw, '(/,A)') " Condition Number of Extended Basis Sets"
      CALL init_orbital_pointers(3)
      CALL init_spherical_harmonics(3, 0)
      DO ider = 0, nder
         NULLIFY (basis_vrb)
         ALLOCATE (basis_vrb)
         NULLIFY (basis_vrb%am, basis_vrb%cm, basis_vrb%as, basis_vrb%ns, basis_vrb%bf, &
                  basis_vrb%dbf, basis_vrb%ddbf)
         ! allocate and initialize the atomic grid
         NULLIFY (basis_vrb%grid)
         CALL allocate_grid_atom(basis_vrb%grid)
         CALL create_grid_atom(basis_vrb%grid, ngp, 1, 1, 0, quadtype)
         basis_vrb%grid%nr = ngp
         !
         basis_vrb%eps_eig = basis_ref%eps_eig
         basis_vrb%geometrical = .FALSE.
         basis_vrb%basis_type = CGTO_BASIS
         basis_vrb%nprim = basis%nprim+next_prim
         basis_vrb%nbas = 0
         DO l = 0, state%maxl_occ
            basis_vrb%nbas(l) = state%maxn_occ(l)+ider+next_bas(l)
         END DO
         m = MAXVAL(basis_vrb%nprim)
         ALLOCATE (basis_vrb%am(m, 0:3))
         ! exponents
         m = SIZE(basis%am, 1)
         basis_vrb%am(1:m, :) = basis%am(1:m, :)
         n = SIZE(ale, 1)
         DO l = 0, state%maxl_occ
            basis_vrb%am(m+1:m+n, l) = ale(1:n)
         END DO
         ! contractions
         m = MAXVAL(basis_vrb%nprim)
         n = MAXVAL(basis_vrb%nbas)
         ALLOCATE (basis_vrb%cm(m, n, 0:3))
         basis_vrb%cm = 0.0_dp
         DO l = 0, state%maxl_occ
            m = basis%nprim(l)
            n = state%maxn_occ(l)+ider
            basis_vrb%cm(1:m, 1:n, l) = rbasis(1:m, 1:n, l)
            basis_vrb%cm(m+1:m+next_prim(l), n+1:n+next_bas(l), l) = ebasis(1:next_prim(l), 1:next_bas(l), l)
         END DO

         ! initialize basis function on a radial grid
         nr = basis_vrb%grid%nr
         m = MAXVAL(basis_vrb%nbas)
         ALLOCATE (basis_vrb%bf(nr, m, 0:3))
         ALLOCATE (basis_vrb%dbf(nr, m, 0:3))
         ALLOCATE (basis_vrb%ddbf(nr, m, 0:3))
         basis_vrb%bf = 0._dp
         basis_vrb%dbf = 0._dp
         basis_vrb%ddbf = 0._dp
         DO l = 0, 3
            DO i = 1, basis_vrb%nprim(l)
               al = basis_vrb%am(i, l)
               DO k = 1, nr
                  rk = basis_vrb%grid%rad(k)
                  ear = EXP(-al*basis_vrb%grid%rad(k)**2)
                  DO j = 1, basis_vrb%nbas(l)
                     basis_vrb%bf(k, j, l) = basis_vrb%bf(k, j, l)+rk**l*ear*basis_vrb%cm(i, j, l)
                     basis_vrb%dbf(k, j, l) = basis_vrb%dbf(k, j, l) &
                                              +(REAL(l, dp)*rk**(l-1)-2._dp*al*rk**(l+1))*ear*basis_vrb%cm(i, j, l)
                     basis_vrb%ddbf(k, j, l) = basis_vrb%ddbf(k, j, l)+ &
                                               (REAL(l*(l-1), dp)*rk**(l-2)-2._dp*al*REAL(2*l+1, dp)*rk**(l)+ &
                                                4._dp*al*rk**(l+2))*ear*basis_vrb%cm(i, j, l)
                  END DO
               END DO
            END DO
         END DO

         IF (iw > 0) THEN
            CALL basis_label(abas, basis_vrb%nprim, basis_vrb%nbas)
            WRITE (iw, '(A,A)') " Basis set     ", TRIM(abas)
         END IF
         crad = 2.0_dp*ptable(atom%z)%covalent_radius*bohr
         cradx = crad*1.00_dp
         CALL atom_basis_condnum(basis_vrb, cradx, cnum)
         IF (iw > 0) WRITE (iw, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
         cradx = crad*1.10_dp
         CALL atom_basis_condnum(basis_vrb, cradx, cnum)
         IF (iw > 0) WRITE (iw, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
         cradx = crad*1.20_dp
         CALL atom_basis_condnum(basis_vrb, cradx, cnum)
         IF (iw > 0) WRITE (iw, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
         vbasis(nder+1+ider)%basis => basis_vrb
      END DO
      CALL deallocate_orbital_pointers
      CALL deallocate_spherical_harmonics

      ! Tests for energy
      energy_ref = atom_ref%energy%etot
      IF (iw > 0) WRITE (iw, '(/,A,A)') " Basis set tests    "
      IF (iw > 0) WRITE (iw, '(T10,A,T59,F22.9)') " Reference Energy [a.u.]  ", energy_ref
      DO ider = 0, 2*nder+1
         ! generate an atom type
         NULLIFY (atom_test)
         CALL create_atom_type(atom_test)
         CALL copy_atom_basics(atom_ref, atom_test, state=.TRUE., potential=.TRUE., &
                               optimization=.TRUE., xc=.TRUE.)
         basis_grb => vbasis(ider)%basis
         NULLIFY (orbitals)
         mo = MAXVAL(atom_test%state%maxn_calc)
         mb = MAXVAL(basis_grb%nbas)
         CALL create_atom_orbs(orbitals, mb, mo)
         CALL set_atom(atom_test, orbitals=orbitals, basis=basis_grb)
         ! calculate integrals
         ALLOCATE (atint)
         CALL atom_int_setup(atint, basis_grb, potential=atom_test%potential, eri_coulomb=.FALSE., eri_exchange=.FALSE.)
         CALL atom_ppint_setup(atint, basis_grb, potential=atom_test%potential)
         IF (atom_test%pp_calc) THEN
            NULLIFY (atint%tzora, atint%hdkh)
         ELSE
            ! relativistic correction terms
            CALL atom_relint_setup(atint, basis_grb, atom_test%relativistic, zcore=REAL(atom_test%z, dp))
         END IF
         CALL set_atom(atom_test, integrals=atint)
         !
         CALL calculate_atom(atom_test, iw=0)
         IF (ider <= nder) THEN
            energy_vb(ider) = atom_test%energy%etot
            IF (iw > 0) WRITE (iw, '(T10,A,i1,A,T40,F13.9,T59,F22.9)') " GRB (VB)", ider, " Energy [a.u.]  ", &
               energy_ref-energy_vb(ider), energy_vb(ider)
         ELSE
            i = ider-nder-1
            energy_ex(i) = atom_test%energy%etot
            IF (iw > 0) WRITE (iw, '(T10,A,i1,A,T40,F13.9,T59,F22.9)') " GRB (EX)", i, " Energy [a.u.]  ", &
               energy_ref-energy_ex(i), energy_ex(i)
         END IF
         CALL atom_int_release(atint)
         CALL atom_ppint_release(atint)
         CALL atom_relint_release(atint)
         DEALLOCATE (atom_test%state, atom_test%potential, atint)
         CALL release_atom_type(atom_test)
      END DO

      ! Quickstep normalization polarization basis
      DO l = 0, 7
         expzet = 0.25_dp*REAL(2*l+3, dp)
         prefac = SQRT(SQRT(pi)/2._dp**(l+2)*dfac(2*l+1))
         DO i = 1, num_pol
            zeta = (2._dp*alp(i))**expzet
            pbasis(i, 1:num_pol, l) = pbasis(i, 1:num_pol, l)*prefac/zeta
         END DO
      END DO
      ! Quickstep normalization extended basis
      DO l = 0, 3
         expzet = 0.25_dp*REAL(2*l+3, dp)
         prefac = SQRT(SQRT(pi)/2._dp**(l+2)*dfac(2*l+1))
         DO i = 1, next_prim(l)
            zeta = (2._dp*ale(i))**expzet
            ebasis(i, 1:next_bas(l), l) = ebasis(i, 1:next_bas(l), l)*prefac/zeta
         END DO
      END DO

      ! Print basis sets
      CALL section_vals_val_get(grb_section, "NAME_BODY", c_val=basname)
      CALL open_file(file_name="GRB_BASIS", file_status="UNKNOWN", file_action="WRITE", unit_number=iunit)
      ! header info
      headline = ""
      headline(1) = "#"
      headline(2) = "# Generated with CP2K Atom Code"
      headline(3) = "#"
      CALL grb_print_basis(header=headline, iunit=iunit)
      ! valence basis
      basline(1) = ""
      WRITE (basline(1), "(T2,A)") ADJUSTL(ptable(atom_ref%z)%symbol)
      DO ider = 0, nder
         basline(1) = ""
         WRITE (basline(1), "(T2,A,T5,A,I1)") ADJUSTL(ptable(atom_ref%z)%symbol), TRIM(ADJUSTL(basname))//"-VAL-", ider
         CALL grb_print_basis(header=basline, nprim=vbasis(ider)%basis%nprim(0), nbas=vbasis(ider)%basis%nbas, &
                              al=vbasis(ider)%basis%am(:, 0), gcc=qbasis, iunit=iunit)
      END DO
      ! polarization basis
      maxl = atom_ref%state%maxl_occ
      DO l = maxl+1, MIN(maxl+num_pol, 7)
         nbas = 0
         DO i = maxl+1, l
            nbas(i) = l-i+1
         END DO
         i = l-maxl
         basline(1) = ""
         WRITE (basline(1), "(T2,A,T5,A,I1)") ADJUSTL(ptable(atom_ref%z)%symbol), TRIM(ADJUSTL(basname))//"-POL-", i
         CALL grb_print_basis(header=basline, nprim=num_pol, nbas=nbas, al=alp, gcc=pbasis, iunit=iunit)
      END DO
      ! extension set
      basline(1) = ""
      WRITE (basline(1), "(T2,A,T5,A)") ADJUSTL(ptable(atom_ref%z)%symbol), TRIM(ADJUSTL(basname))//"-EXT"
      CALL grb_print_basis(header=basline, nprim=next_prim(0), nbas=next_bas, al=ale, gcc=ebasis, iunit=iunit)
      !
      CALL close_file(unit_number=iunit)

      ! clean up
      DEALLOCATE (wfn, rbasis, qbasis, ebasis, pbasis, ale, alp)

      DO ider = 0, 10
         IF (ASSOCIATED(vbasis(ider)%basis)) THEN
            CALL release_atom_basis(vbasis(ider)%basis)
            DEALLOCATE (vbasis(ider)%basis)
         END IF
      END DO

      CALL atom_int_release(atom%integrals)
      CALL atom_ppint_release(atom%integrals)
      CALL atom_relint_release(atom%integrals)
      CALL release_atom_basis(basis)
      DEALLOCATE (atom%potential, atom%state, atom%integrals, basis)
      CALL release_atom_type(atom)

      IF (iw > 0) WRITE (iw, '(" ",79("*"))')

   END SUBROUTINE atom_grb_construction

! **************************************************************************************************
!> \brief ...
!> \param header ...
!> \param nprim ...
!> \param nbas ...
!> \param al ...
!> \param gcc ...
!> \param iunit ...
! **************************************************************************************************
   SUBROUTINE grb_print_basis(header, nprim, nbas, al, gcc, iunit)
      CHARACTER(len=*), DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: header
      INTEGER, INTENT(IN), OPTIONAL                      :: nprim
      INTEGER, DIMENSION(0:), INTENT(IN), OPTIONAL       :: nbas
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: al
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN), &
         OPTIONAL                                        :: gcc
      INTEGER, INTENT(IN)                                :: iunit

      CHARACTER(len=*), PARAMETER :: routineN = 'grb_print_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, l, lmax, lmin, nval

      IF (PRESENT(header)) THEN
         DO i = 1, SIZE(header, 1)
            IF (header(i) .NE. "") THEN
               WRITE (iunit, "(A)") TRIM(header(i))
            END IF
         END DO
      END IF

      IF (PRESENT(nprim)) THEN
         IF (nprim > 0) THEN
            CPASSERT(PRESENT(nbas))
            CPASSERT(PRESENT(al))
            CPASSERT(PRESENT(gcc))
            lmin = 100
            lmax = 0
            DO i = 0, SIZE(nbas, 1)-1
               IF (nbas(i) > 0) THEN
                  lmin = MIN(lmin, i)
                  lmax = MAX(lmax, i)
               END IF
            END DO
            nval = lmax
            WRITE (iunit, *) "  1"
            WRITE (iunit, "(40I3)") nval, lmin, lmax, nprim, (nbas(l), l=lmin, lmax)
            DO i = 1, nprim
               WRITE (iunit, "(G20.12,40F16.10)") al(i), ((gcc(i, j, l), j=1, nbas(l)), l=lmin, lmax)
            END DO
            WRITE (iunit, *)
         END IF
      END IF

   END SUBROUTINE grb_print_basis

! **************************************************************************************************
!> \brief ...
!> \param label ...
!> \param np ...
!> \param nb ...
! **************************************************************************************************
   SUBROUTINE basis_label(label, np, nb)
      CHARACTER(len=*)                                   :: label
      INTEGER, DIMENSION(0:)                             :: np, nb

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_label', routineP = moduleN//':'//routineN
      INTEGER                                            :: i, l, lmax
      CHARACTER(len=1), DIMENSION(0:7), PARAMETER :: &
         lq = (/"s", "p", "d", "f", "g", "h", "i", "k"/)

      label = ""
      lmax = MIN(UBOUND(np, 1), UBOUND(nb, 1), 7)
      i = 1
      label(i:i) = "("
      DO l = 0, lmax
         IF (np(l) > 0) THEN
            i = i+1
            IF (np(l) > 9) THEN
               WRITE (label(i:i+1), "(I2)") np(l)
               i = i+2
            ELSE
               WRITE (label(i:i), "(I1)") np(l)
               i = i+1
            END IF
            label(i:i) = lq(l)
         END IF
      END DO
      i = i+1
      label(i:i+6) = ") --> ["
      i = i+6
      DO l = 0, lmax
         IF (nb(l) > 0) THEN
            i = i+1
            IF (nb(l) > 9) THEN
               WRITE (label(i:i+1), "(I2)") nb(l)
               i = i+2
            ELSE
               WRITE (label(i:i), "(I1)") nb(l)
               i = i+1
            END IF
            label(i:i) = lq(l)
         END IF
      END DO
      i = i+1
      label(i:i) = "]"

   END SUBROUTINE basis_label

! **************************************************************************************************
!> \brief ...
!> \param atom ...
!> \param basis ...
!> \param afun ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE grb_fit(atom, basis, afun, iw)
      TYPE(atom_type), POINTER                           :: atom
      TYPE(atom_basis_type), POINTER                     :: basis
      REAL(dp), INTENT(OUT)                              :: afun
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'grb_fit', routineP = moduleN//':'//routineN

      INTEGER                                            :: do_eric, do_erie, reltyp, zval
      LOGICAL                                            :: eri_c, eri_e
      TYPE(atom_integrals), POINTER                      :: atint
      TYPE(atom_potential_type), POINTER                 :: pot

      ALLOCATE (atint)
      ! calculate integrals
      NULLIFY (pot)
      eri_c = .FALSE.
      eri_e = .FALSE.
      pot => atom%potential
      zval = atom%z
      reltyp = atom%relativistic
      do_eric = atom%coulomb_integral_type
      do_erie = atom%exchange_integral_type
      IF (do_eric == do_analytic) eri_c = .TRUE.
      IF (do_erie == do_analytic) eri_e = .TRUE.
      ! general integrals
      CALL atom_int_setup(atint, basis, potential=pot, eri_coulomb=eri_c, eri_exchange=eri_e)
      ! potential
      CALL atom_ppint_setup(atint, basis, potential=pot)
      IF (atom%pp_calc) THEN
         NULLIFY (atint%tzora, atint%hdkh)
      ELSE
         ! relativistic correction terms
         CALL atom_relint_setup(atint, basis, reltyp, zcore=REAL(zval, dp))
      END IF
      CALL set_atom(atom, basis=basis)
      CALL set_atom(atom, integrals=atint)
      CALL calculate_atom(atom, iw)
      afun = atom%energy%etot
      CALL atom_int_release(atint)
      CALL atom_ppint_release(atint)
      CALL atom_relint_release(atint)
      DEALLOCATE (atint)
   END SUBROUTINE grb_fit

! **************************************************************************************************
!> \brief ...
!> \param atom_ref ...
!> \param atom ...
!> \param state ...
!> \param potential ...
!> \param optimization ...
!> \param xc ...
! **************************************************************************************************
   SUBROUTINE copy_atom_basics(atom_ref, atom, state, potential, optimization, xc)
      TYPE(atom_type), POINTER                           :: atom_ref, atom
      LOGICAL, INTENT(IN), OPTIONAL                      :: state, potential, optimization, xc

      atom%z = atom_ref%z
      atom%zcore = atom_ref%zcore
      atom%pp_calc = atom_ref%pp_calc
      atom%method_type = atom_ref%method_type
      atom%relativistic = atom_ref%relativistic
      atom%coulomb_integral_type = atom_ref%coulomb_integral_type
      atom%exchange_integral_type = atom_ref%exchange_integral_type

      NULLIFY (atom%potential, atom%state, atom%xc_section)
      NULLIFY (atom%basis, atom%integrals, atom%orbitals, atom%fmat)

      IF (PRESENT(state)) THEN
         IF (state) THEN
            ALLOCATE (atom%state)
            atom%state = atom_ref%state
         END IF
      END IF

      IF (PRESENT(potential)) THEN
         IF (potential) THEN
            ALLOCATE (atom%potential)
            atom%potential = atom_ref%potential
         END IF
      END IF

      IF (PRESENT(optimization)) THEN
         IF (optimization) THEN
            atom%optimization = atom_ref%optimization
         END IF
      END IF

      IF (PRESENT(xc)) THEN
         IF (xc) THEN
            atom%xc_section => atom_ref%xc_section
         END IF
      END IF

   END SUBROUTINE copy_atom_basics

! **************************************************************************************************
!> \brief ...
!> \param atom ...
!> \param basis ...
!> \param iunit ...
!> \param powell_section ...
! **************************************************************************************************
   SUBROUTINE atom_fit_grb(atom, basis, iunit, powell_section)
      TYPE(atom_type), POINTER                           :: atom
      TYPE(atom_basis_type), POINTER                     :: basis
      INTEGER, INTENT(IN)                                :: iunit
      TYPE(section_vals_type), POINTER                   :: powell_section

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_grb', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, k, l, ll, n10, nr
      REAL(KIND=dp)                                      :: al, cnum, crad, cradx, ear, fopt, rk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x
      TYPE(opt_state_type)                               :: ostate

      CPASSERT(basis%geometrical)

      CALL section_vals_val_get(powell_section, "ACCURACY", r_val=ostate%rhoend)
      CALL section_vals_val_get(powell_section, "STEP_SIZE", r_val=ostate%rhobeg)
      CALL section_vals_val_get(powell_section, "MAX_FUN", i_val=ostate%maxfun)

      ostate%nvar = 2
      ALLOCATE (x(2))
      x(1) = SQRT(basis%aval)
      x(2) = SQRT(basis%cval)

      ostate%nf = 0
      ostate%iprint = 1
      ostate%unit = iunit

      ostate%state = 0
      IF (iunit > 0) THEN
         WRITE (iunit, '(/," POWELL| Start optimization procedure")')
         WRITE (iunit, '(" POWELL| Total number of parameters in optimization",T71,I10)') ostate%nvar
      END IF
      n10 = MAX(ostate%maxfun/100, 1)

      fopt = HUGE(0._dp)

      DO

         IF (ostate%state == 2) THEN
            basis%am = 0._dp
            DO l = 0, 3
               DO i = 1, basis%nbas(l)
                  ll = i-1+basis%start(l)
                  basis%am(i, l) = x(1)*x(1)*(x(2)*x(2))**(ll)
               END DO
            END DO
            basis%aval = x(1)*x(1)
            basis%cval = x(2)*x(2)
            basis%bf = 0._dp
            basis%dbf = 0._dp
            basis%ddbf = 0._dp
            nr = basis%grid%nr
            DO l = 0, 3
               DO i = 1, basis%nbas(l)
                  al = basis%am(i, l)
                  DO k = 1, nr
                     rk = basis%grid%rad(k)
                     ear = EXP(-al*basis%grid%rad(k)**2)
                     basis%bf(k, i, l) = rk**l*ear
                     basis%dbf(k, i, l) = (REAL(l, dp)*rk**(l-1)-2._dp*al*rk**(l+1))*ear
                     basis%ddbf(k, i, l) = (REAL(l*(l-1), dp)*rk**(l-2)- &
                                            2._dp*al*REAL(2*l+1, dp)*rk**(l)+4._dp*al*rk**(l+2))*ear
                  END DO
               END DO
            END DO
            CALL grb_fit(atom, basis, ostate%f, 0)
            fopt = MIN(fopt, ostate%f)
         END IF

         IF (ostate%state == -1) EXIT

         CALL powell_optimize(ostate%nvar, x, ostate)

         IF (ostate%nf == 2 .AND. iunit > 0) THEN
            WRITE (iunit, '(" POWELL| Inital value of function",T61,F20.10)') ostate%f
         END IF
         IF (MOD(ostate%nf, n10) == 0 .AND. iunit > 0) THEN
            WRITE (iunit, '(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
               INT(REAL(ostate%nf, dp)/REAL(ostate%maxfun, dp)*100._dp), fopt
         END IF

      END DO

      ostate%state = 8
      CALL powell_optimize(ostate%nvar, x, ostate)

      IF (iunit > 0) THEN
         WRITE (iunit, '(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
         WRITE (iunit, '(" POWELL| Final value of function",T61,F20.10)') ostate%fopt
      END IF
      ! x->basis
      basis%am = 0._dp
      DO l = 0, 3
         DO i = 1, basis%nbas(l)
            ll = i-1+basis%start(l)
            basis%am(i, l) = x(1)*x(1)*(x(2)*x(2))**(ll)
         END DO
      END DO
      basis%aval = x(1)*x(1)
      basis%cval = x(2)*x(2)
      basis%bf = 0._dp
      basis%dbf = 0._dp
      basis%ddbf = 0._dp
      nr = basis%grid%nr
      DO l = 0, 3
         DO i = 1, basis%nbas(l)
            al = basis%am(i, l)
            DO k = 1, nr
               rk = basis%grid%rad(k)
               ear = EXP(-al*basis%grid%rad(k)**2)
               basis%bf(k, i, l) = rk**l*ear
               basis%dbf(k, i, l) = (REAL(l, dp)*rk**(l-1)-2._dp*al*rk**(l+1))*ear
               basis%ddbf(k, i, l) = (REAL(l*(l-1), dp)*rk**(l-2)- &
                                      2._dp*al*REAL(2*l+1, dp)*rk**(l)+4._dp*al*rk**(l+2))*ear
            END DO
         END DO
      END DO

      DEALLOCATE (x)

      ! final result
      IF (iunit > 0) THEN
         WRITE (iunit, '(/,A)') " Optimized Geometrical GTO basis set"
         WRITE (iunit, '(A,F15.8,T41,A,F15.8)') " Initial exponent: ", basis%aval, &
            " Proportionality factor: ", basis%cval
         DO l = 0, 3
            WRITE (iunit, '(T41,A,I2,T76,I5)') " Number of exponents for l=", l, basis%nbas(l)
         END DO
      END IF

      IF (iunit > 0) WRITE (iunit, '(/,A)') " Condition number of uncontracted basis set"
      crad = 2.0_dp*ptable(atom%z)%covalent_radius*bohr
      CALL init_orbital_pointers(3)
      CALL init_spherical_harmonics(3, 0)
      cradx = crad*1.00_dp
      CALL atom_basis_condnum(basis, cradx, cnum)
      IF (iunit > 0) WRITE (iunit, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
      cradx = crad*1.10_dp
      CALL atom_basis_condnum(basis, cradx, cnum)
      IF (iunit > 0) WRITE (iunit, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
      cradx = crad*1.20_dp
      CALL atom_basis_condnum(basis, cradx, cnum)
      IF (iunit > 0) WRITE (iunit, '(T5,A,F15.3,T50,A,F14.4)') " Lattice constant:", cradx, "Condition number:", cnum
      CALL deallocate_orbital_pointers
      CALL deallocate_spherical_harmonics

   END SUBROUTINE atom_fit_grb
! **************************************************************************************************
!> \brief ...
!> \param zval ...
!> \param rconf ...
!> \param lval ...
!> \param aval ...
!> \param cval ...
!> \param nbas ...
!> \param iunit ...
!> \param powell_section ...
! **************************************************************************************************
   SUBROUTINE atom_fit_pol(zval, rconf, lval, aval, cval, nbas, iunit, powell_section)
      REAL(KIND=dp), INTENT(IN)                          :: zval, rconf
      INTEGER, INTENT(IN)                                :: lval
      REAL(KIND=dp), INTENT(INOUT)                       :: aval, cval
      INTEGER, INTENT(IN)                                :: nbas, iunit
      TYPE(section_vals_type), POINTER                   :: powell_section

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_pol', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, n10
      REAL(KIND=dp)                                      :: fopt, x(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: am, ener
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: orb
      TYPE(opt_state_type)                               :: ostate

      ALLOCATE (am(nbas), ener(nbas), orb(nbas, nbas))

      CALL section_vals_val_get(powell_section, "ACCURACY", r_val=ostate%rhoend)
      CALL section_vals_val_get(powell_section, "STEP_SIZE", r_val=ostate%rhobeg)
      CALL section_vals_val_get(powell_section, "MAX_FUN", i_val=ostate%maxfun)

      ostate%nvar = 2
      x(1) = SQRT(aval)
      x(2) = SQRT(cval)

      ostate%nf = 0
      ostate%iprint = 1
      ostate%unit = iunit

      ostate%state = 0
      IF (iunit > 0) THEN
         WRITE (iunit, '(/," POWELL| Start optimization procedure")')
         WRITE (iunit, '(" POWELL| Total number of parameters in optimization",T71,I10)') ostate%nvar
      END IF
      n10 = MAX(ostate%maxfun/100, 1)

      fopt = HUGE(0._dp)

      DO

         IF (ostate%state == 2) THEN
            aval = x(1)*x(1)
            cval = x(2)*x(2)
            DO i = 1, nbas
               am(i) = aval*cval**(i-1)
            END DO
            CALL hydrogenic(zval, rconf, lval, am, nbas, ener, orb)
            ostate%f = ener(1)
            fopt = MIN(fopt, ostate%f)
         END IF

         IF (ostate%state == -1) EXIT

         CALL powell_optimize(ostate%nvar, x, ostate)

         IF (ostate%nf == 2 .AND. iunit > 0) THEN
            WRITE (iunit, '(" POWELL| Inital value of function",T61,F20.10)') ostate%f
         END IF
         IF (MOD(ostate%nf, n10) == 0 .AND. iunit > 0) THEN
            WRITE (iunit, '(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
               INT(REAL(ostate%nf, dp)/REAL(ostate%maxfun, dp)*100._dp), fopt
         END IF

      END DO

      ostate%state = 8
      CALL powell_optimize(ostate%nvar, x, ostate)

      IF (iunit > 0) THEN
         WRITE (iunit, '(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
         WRITE (iunit, '(" POWELL| Final value of function",T61,F20.10)') ostate%fopt
      END IF
      ! x->basis
      aval = x(1)*x(1)
      cval = x(2)*x(2)

      ! final result
      IF (iunit > 0) THEN
         WRITE (iunit, '(/,A,T51,A,T76,I5)') " Optimized Polarization basis set", &
            " Number of exponents:", nbas
         WRITE (iunit, '(A,F15.8,T41,A,F15.8)') " Initial exponent: ", aval, &
            " Proportionality factor: ", cval
      END IF

      DEALLOCATE (am, ener, orb)

   END SUBROUTINE atom_fit_pol
! **************************************************************************************************
!> \brief ...
!> \param zval ...
!> \param rconf ...
!> \param lval ...
!> \param am ...
!> \param nbas ...
!> \param ener ...
!> \param orb ...
! **************************************************************************************************
   SUBROUTINE hydrogenic(zval, rconf, lval, am, nbas, ener, orb)
      REAL(KIND=dp), INTENT(IN)                          :: zval, rconf
      INTEGER, INTENT(IN)                                :: lval
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: am
      INTEGER, INTENT(IN)                                :: nbas
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: ener
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: orb

      CHARACTER(len=*), PARAMETER :: routineN = 'hydrogenic', routineP = moduleN//':'//routineN

      INTEGER                                            :: info, k, lwork, n
      REAL(KIND=dp)                                      :: cf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: confmat, hmat, potmat, smat, tmat

      n = nbas
      ALLOCATE (smat(n, n), tmat(n, n), potmat(n, n), confmat(n, n), hmat(n, n))
      ! calclulate overlap matrix
      CALL sg_overlap(smat(1:n, 1:n), lval, am(1:n), am(1:n))
      ! calclulate kinetic energy matrix
      CALL sg_kinetic(tmat(1:n, 1:n), lval, am(1:n), am(1:n))
      ! calclulate core potential matrix
      CALL sg_nuclear(potmat(1:n, 1:n), lval, am(1:n), am(1:n))
      ! calclulate confinement potential matrix
      cf = 0.1_dp
      k = 10
      CALL sg_conf(confmat, rconf, k, lval, am(1:n), am(1:n))
      ! Hamiltionian
      hmat(1:n, 1:n) = tmat(1:n, 1:n)-zval*potmat(1:n, 1:n)+cf*confmat(1:n, 1:n)
      ! solve
      lwork = 100*n
      ALLOCATE (w(n), work(lwork))
      CALL lapack_ssygv(1, "V", "U", n, hmat, n, smat, n, w, work, lwork, info)
      CPASSERT(info == 0)
      orb(1:n, 1:n) = hmat(1:n, 1:n)
      ener(1:n) = w(1:n)
      DEALLOCATE (w, work)
      DEALLOCATE (smat, tmat, potmat, confmat, hmat)

   END SUBROUTINE hydrogenic

END MODULE atom_grb
