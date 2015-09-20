!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Independent helium subroutines shared by other modules
!> \author Lukasz Walewski
!> \date   2009-07-14
!> \note   Avoiding circular deps: do not USE any other helium_* modules here.
! *****************************************************************************
MODULE helium_common

  USE helium_types,                    ONLY: he_mass,&
                                             helium_solvent_type
  USE input_constants,                 ONLY: helium_cell_shape_cube,&
                                             helium_cell_shape_octahedron
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE parallel_rng_types,              ONLY: next_random_number
  USE physcon,                         ONLY: a_mass,&
                                             angstrom,&
                                             boltzmann,&
                                             h_bar,&
                                             kelvin
  USE splines_methods,                 ONLY: spline_value
  USE splines_types,                   ONLY: spline_data_p_type,&
                                             spline_data_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_common'

  PUBLIC :: helium_pbc
  PUBLIC :: helium_boxmean_3d
  PUBLIC :: helium_calc_wnumber
  PUBLIC :: helium_calc_rdf
  PUBLIC :: helium_calc_rho
  PUBLIC :: helium_calc_plength
  PUBLIC :: helium_rotate
  PUBLIC :: helium_eval_expansion
  PUBLIC :: helium_update_transition_matrix
  PUBLIC :: helium_spline
  PUBLIC :: helium_cycle_number
  PUBLIC :: helium_path_length
  PUBLIC :: helium_is_winding

  CONTAINS

! *****************************************************************************
!> \brief  General PBC routine for helium.
!> \param helium ...
!> \param r ...
!> \date   2009-09-30
!> \author Lukasz Walewski
!> \note   Check wheather PBC should be applied, if yes call low level 
!>         routine according to the unit cell shape.
! *****************************************************************************
  SUBROUTINE helium_pbc( helium, r )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(kind=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    IF ( helium%cell_shape .EQ. helium_cell_shape_cube ) THEN
      CALL helium_pbc_cube( helium, r )
    ELSE IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
      CALL helium_pbc_trocta_opt( helium, r )
    END IF

  END SUBROUTINE helium_pbc

! *****************************************************************************
!> \brief  Wrap r back to the helium box (if helium box is periodic)
!> \param    helium - helium environment for which to wrap
!> \param    r      - 3D vector to be wraped back to periodic box
!> \par History
!>      2009-10-02 renamed, originally was helium_box
!>      2009-10-02 redesigned so it is now called as a subroutine [lwalewski]
!>      2009-10-02 redesigned so it now gets/returns a 3D vector [lwalewski]
!> \author hforbert
! *****************************************************************************
  SUBROUTINE helium_pbc_cube( helium, r )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(kind=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    REAL(kind=dp)                            :: s

! x coordinate

      s = helium%cell_size_inv * r(1)
      IF ( s > 0.5_dp ) THEN
        s = s - INT(s+0.5_dp)
      ELSEIF ( s < -0.5_dp ) THEN
        s = s - INT(s-0.5_dp)
      END IF
      r(1) = s * helium%cell_size

      ! y coordinate
      s = helium%cell_size_inv * r(2)
      IF ( s > 0.5_dp ) THEN
        s = s - INT(s+0.5_dp)
      ELSEIF ( s < -0.5_dp ) THEN
        s = s - INT(s-0.5_dp)
      END IF
      r(2) = s * helium%cell_size

      ! z coordinate
      s = helium%cell_size_inv * r(3)
      IF ( s > 0.5_dp ) THEN
        s = s - INT(s+0.5_dp)
      ELSEIF ( s < -0.5_dp ) THEN
        s = s - INT(s-0.5_dp)
      END IF
      r(3) = s * helium%cell_size

    RETURN
  END SUBROUTINE helium_pbc_cube

! *****************************************************************************
!> \brief  Apply PBC within truncated octahedral unit cell.
!> \param helium ...
!> \param r ...
!> \date   2009-10-02
!> \author Lukasz Walewski
!> \note   Original Allen & Tildesley routine adapted for our helium code.
! *****************************************************************************
  SUBROUTINE helium_pbc_trocta( helium, r )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(kind=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    REAL(kind=dp)                            :: corr, r75, rx, ry, rz

! ** PERIODIC BOUNDARY CONDITIONS FOR A TRUNCATED OCTAHEDRON       **
! **                                                               **
! ** THE BOX IS CENTRED AT THE ORIGIN. THE AXES PASS THROUGH THE   **
! ** CENTRES OF THE SIX SQUARE FACES OF THE TRUNCATED OCTAHEDRON   **
! ** (SEE F1G. 1.10(A)). THE CONTAINING CUBE IS OF UNIT LENGTH     **

    PARAMETER ( r75 = 4.0_dp / 3.0_dp )

    rx = r(1) * helium%cell_size_inv
    ry = r(2) * helium%cell_size_inv
    rz = r(3) * helium%cell_size_inv

    rx = rx - ANINT ( rx )
    ry = ry - ANINT ( ry )
    rz = rz - ANINT ( rz )
    corr = 0.5_dp * AINT ( r75 * ( ABS ( rx ) + ABS ( ry ) + ABS ( rz ) ) )
    rx = rx - SIGN ( corr, rx )
    ry = ry - SIGN ( corr, ry )
    rz = rz - SIGN ( corr, rz )

    r(1) = rx * helium%cell_size
    r(2) = ry * helium%cell_size
    r(3) = rz * helium%cell_size

    RETURN
  END SUBROUTINE helium_pbc_trocta

! *****************************************************************************
!> \brief  Apply PBC within truncated octahedral unit cell.
!> \param helium ...
!> \param r ...
!> \date   2009-10-22
!> \author hforbert
!> \note   Version of the original Allen & Tildesley routine optimized for
!>         g95 and intel compilers on x86-64.
! *****************************************************************************
  SUBROUTINE helium_pbc_trocta_opt( helium, r )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(kind=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    REAL(kind=dp)                            :: cell_size, cell_size_inv, &
                                                corr, rx, ry, rz, sx, sy, sz

    cell_size = helium%cell_size
    cell_size_inv = helium%cell_size_inv

    rx = r(1) * cell_size_inv
    IF ( rx > 0.5_dp ) THEN
       rx = rx - INT(rx+0.5_dp)
    ELSEIF ( rx < -0.5_dp ) THEN
       rx = rx - INT(rx-0.5_dp)
    END IF

    ry = r(2) * cell_size_inv
    IF ( ry > 0.5_dp ) THEN
       ry = ry - INT(ry+0.5_dp)
    ELSEIF ( ry < -0.5_dp ) THEN
       ry = ry - INT(ry-0.5_dp)
    END IF

    rz = r(3) * cell_size_inv
    IF ( rz > 0.5_dp ) THEN
       rz = rz - INT(rz+0.5_dp)
    ELSEIF ( rz < -0.5_dp ) THEN
       rz = rz - INT(rz-0.5_dp)
    END IF

    corr = 0.0_dp
    IF ( rx > 0.0_dp ) THEN
       corr = corr + rx
       sx = 0.5_dp
    ELSE
       corr = corr - rx
       sx = -0.5_dp
    END IF
    IF ( ry > 0.0_dp ) THEN
       corr = corr + ry
       sy = 0.5_dp
    ELSE
       corr = corr - ry
       sy = -0.5_dp
    END IF
    IF ( rz > 0.0_dp ) THEN
       corr = corr + rz
       sz = 0.5_dp
    ELSE
       corr = corr - rz
       sz = -0.5_dp
    END IF
    IF ( corr > 0.75_dp ) THEN
       rx = rx - sx
       ry = ry - sy
       rz = rz - sz
    END IF

    r(1) = rx * cell_size
    r(2) = ry * cell_size
    r(3) = rz * cell_size

    RETURN
  END SUBROUTINE helium_pbc_trocta_opt

! *****************************************************************************
!> \brief  Calculate the point equidistant from two other points a and b
!>         within the helium box - 3D version
!> \param    helium - helium environment for which
!> \param a vectors for which to find the mean within the He box
!> \param b vectors for which to find the mean within the He box
!> \param c ...
!> \par History
!>      2009-10-02 renamed, originally was helium_boxmean [lwalewski]
!>      2009-10-02 redesigned so it is now called as a subroutine [lwalewski]
!>      2009-10-02 redesigned so it now gets/returns a 3D vectors [lwalewski]
!> \author hforbert
! *****************************************************************************
  SUBROUTINE helium_boxmean_3d(helium, a, b, c)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: c

    c(:) = b(:) - a(:)
    CALL helium_pbc( helium, c )
    c(:) = a(:) + 0.5_dp * c(:)
    CALL helium_pbc( helium, c )
    RETURN
  END SUBROUTINE helium_boxmean_3d

! ***************************************************************************
!> \brief Calculate the winding number for helium
!> \param helium ...
!> \par   History
!>        2009-10-19 generalized for trunc. octahedron cell shape [lwalewski]
!> \author hforbert
! *****************************************************************************
  SUBROUTINE helium_calc_wnumber( helium )

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_wnumber', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib
    REAL(KIND=dp)                            :: c, mHe, T
    REAL(KIND=dp), DIMENSION(3)              :: r, v, w

!TODO combined with my standalone program - simpify the units
! h_bar - Planck constant over 2 pi in [J*s]
! boltzmann = 1.3806504e-23 - Boltzmann constant in [J/K]
! he_mass - 4He mass in [u]
! a_mass = 1.660538782e-27 - atomic mass unit in [kg]

    mHe  = he_mass * a_mass ! 4He mass in [kg]
    T    = kelvin / helium%tau / helium%beads
    c    = 1e-20 * mHe / h_bar / h_bar * boltzmann * T / 3.0_dp / helium%atoms

    v(:) = 0.0_dp
    ! iterate over all helium atoms
    DO ia = 1, helium%atoms
      ! contribution comming from the last and the first bead
      r(:) = helium%pos(:,ia,helium%beads) - &
        helium%pos(:,helium%permutation(ia),1)
      CALL helium_pbc( helium, r )
      v(:) = v(:) + r(:)
      ! sum of contributions from the rest of bead pairs
      DO ib = 2, helium%beads
        r(:) = helium%pos(:,ia,ib-1) - helium%pos(:,ia,ib)
        CALL helium_pbc( helium, r )
        v(:) = v(:) + r(:)
      END DO
    END DO

    ! calc winding number (cell_m_inv is in Bohr^-1)
    w(1) = &
      helium%cell_m_inv(1,1) * v(1) + &
      helium%cell_m_inv(1,2) * v(2) + &
      helium%cell_m_inv(1,3) * v(3);
    w(2) = &
      helium%cell_m_inv(2,1) * v(1) + &
      helium%cell_m_inv(2,2) * v(2) + &
      helium%cell_m_inv(2,3) * v(3);
    w(3) = &
      helium%cell_m_inv(3,1) * v(1) + &
      helium%cell_m_inv(3,2) * v(2) + &
      helium%cell_m_inv(3,3) * v(3);
    helium%wnumber_inst(:) = w(:)

    ! calc superfluid density
    v(:) = angstrom * v(:)
    helium%sdensity_inst = c * ( v(1)**2 + v(2)**2 + v(3)**2 )

    RETURN
  END SUBROUTINE helium_calc_wnumber

! *****************************************************************************
!> \brief  Calculate helium radial distribution function wrt <r0>.
!> \param helium ...
!> \date   2009-07-22
!> \author Lukasz Walewski
!> \note   Actually calculate the histogram only, the normalization is
!>         postponed to the postprocessing stage.
! *****************************************************************************
  SUBROUTINE helium_calc_rdf( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_rdf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, ia, ib
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ri
    REAL(kind=dp), DIMENSION(3)              :: r, r0

!    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,failure)

    ! calculate the histogram of distances
    helium%rdf_inst(:) = 0.0_dp
    r0(:) = helium%origin(:)
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        r(:) = helium%pos(:,ia,ib) - r0(:)
        CALL helium_pbc( helium, r )
        ri = SQRT( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
        bin = INT(ri/helium%rdf_delr) + 1
        IF (bin .LE. helium%rdf_nbin) THEN
          helium%rdf_inst(bin) = helium%rdf_inst(bin) + 1.0_dp
        END IF
      END DO
    END DO

    ! normalize the histogram to get g(r)
    ! note: helium%density refers to the number of atoms, not the beads
!    norm = 0.0_dp
!    const = 4.0_dp * pi * helium%density / 3.0_dp
!    rpart = REAL(helium%beads)
!    DO bin = 1, helium%rdf_nbin
!      rlower = REAL(bin-1) * helium%rdf_delr
!      rupper = rlower + helium%rdf_delr
!      nideal = const * (rupper**3 - rlower**3)
!      helium%rdf_inst(bin) = helium%rdf_inst(bin) / rpart / nideal
!      norm = norm + helium%rdf_inst(bin)
!    END DO
!    norm = norm * helium%rdf_delr

!    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_calc_rdf


! *****************************************************************************
!> \brief  Calculate helium density distribution function wrt <r0>.
!> \param helium ...
!> \date   2011-06-14
!> \author Lukasz Walewski
!> \note   The calculated density is stored in the helium%rho_inst array. 
!> \note   This version of helium_calc_rho routine calculates only the total
!>         He density but it implements the data structure that is already
!>         suitable for the superfluid density estimator calculation
!>         (full version still in the development branch) [lwalewski]
! *****************************************************************************
  SUBROUTINE helium_calc_rho( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bx, by, bz, handle, ia, ib, &
                                                itmp
    LOGICAL                                  :: failure, ltmp
    REAL(kind=dp)                            :: c, delrinv
    REAL(kind=dp), DIMENSION(3)              :: r, r0

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,failure)

    helium%rho_inst(:,:,:,:) = 0.0_dp
    r0(:) = helium%origin(:)
    delrinv = 1.0_dp / helium%rho_delr
    c = 1.0_dp / helium%beads

    DO ia = 1, helium%atoms

      DO ib = 1, helium%beads

        r(:) = helium%pos(:,ia,ib) - r0(:)
        CALL helium_pbc( helium, r )
        bx = FLOOR( r(1) * delrinv ) - helium%rho_minb + 1
        by = FLOOR( r(2) * delrinv ) - helium%rho_minb + 1
        bz = FLOOR( r(3) * delrinv ) - helium%rho_minb + 1

        ! check that bin numbers are within bounds
        itmp = helium%rho_nbin
        ltmp = (0 .LT. bx) .AND. ( bx .LE. itmp )
        CPPrecondition(ltmp,cp_failure_level,routineP,failure)
        ltmp = (0 .LT. by) .AND. ( by .LE. itmp )
        CPPrecondition(ltmp,cp_failure_level,routineP,failure)
        ltmp = (0 .LT. bz) .AND. ( bz .LE. itmp )
        CPPrecondition(ltmp,cp_failure_level,routineP,failure)

        helium%rho_inst(1,bx,by,bz) = helium%rho_inst(1,bx,by,bz) + c

      END DO
    END DO

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_calc_rho


! ***************************************************************************
!> \brief  Calculate probability distribution of the permutation lengths
!> \param helium ...
!> \date   2010-06-07
!> \author Lukasz Walewski
!> \note  Valid permutation path length is an integer (1, NATOMS), number
!>         of paths of a given length is calculated here and average over
!>         inner loop iterations and helium environments is done in
!>         helium_sample.
! *****************************************************************************
  SUBROUTINE helium_calc_plength( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_plength', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k

    helium%plength_inst(:) = 0.0_dp
    DO i = 1, helium%atoms
      j = helium%permutation(i)
      k = 1
      DO
        IF (j == i) EXIT
        k = k + 1
        j = helium%permutation(j)
      END DO
      helium%plength_inst(k) = helium%plength_inst(k) + 1
    END DO
    helium%plength_inst(:) = helium%plength_inst(:) / helium%atoms

    RETURN
  END SUBROUTINE helium_calc_plength

! ***************************************************************************
!> \brief  Rotate helium particles in imaginary time by nslices
!> \param helium ...
!> \param nslices ...
!> \author hforbert
!> \note  Positions of helium beads in helium%pos array are reorganized such
!>         that the indices are cyclically translated in a permutation-aware
!>         manner. helium%relrot is given a new value that represents the new
!>         'angle' of the beads. This is done modulo helium%beads, so relrot
!>         should be always within 0 (no rotation) and helium%beads-1 (almost
!>         full rotation). [lwalewski]
! *****************************************************************************
  SUBROUTINE helium_rotate(helium, nslices)
    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: nslices

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rotate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, i, j, k, n
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,failure)

    b = helium%beads
    n = helium%atoms
    i = MOD(nslices,b)
    IF (i < 0) i = i + b
    IF ((i>=b).OR.(i<1)) RETURN
    helium%relrot = MOD(helium%relrot+i,b)
    DO k = 1, i
      helium%work(:,:,k) = helium%pos(:,:,k)
    END DO
    DO k = i+1, b
      helium%pos(:,:,k-i) = helium%pos(:,:,k)
    END DO
    DO k = 1, i
      DO j = 1, n
        helium%pos(:,j,b-i+k) = helium%work(:,helium%permutation(j),k)
      END DO
    END DO
    RETURN
  END SUBROUTINE helium_rotate

! *****************************************************************************
!> \brief  Calculate the pair-product action or energy by evaluating the
!>         power series expansion according to Eq. 4.46 in Ceperley 1995.
!> \param helium ...
!> \param r ...
!> \param rp ...
!> \param tab ...
!> \param cut ...
!> \retval res ...
!> \author Harald Forbert
! *****************************************************************************
  FUNCTION helium_eval_expansion(helium,r,rp,tab,cut) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r, rp
    TYPE(spline_data_p_type), &
      DIMENSION(:, :), POINTER               :: tab
    INTEGER, INTENT(IN)                      :: cut
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: ar, arp, q, s, v, z
    REAL(KIND=dp), DIMENSION(3)              :: br, brp

    br(:) = r(:)
    brp(:) = rp(:)
    CALL helium_pbc(helium,br)
    CALL helium_pbc(helium,brp)

    ar = SQRT(br(1)**2+br(2)**2+br(3)**2)
    arp= SQRT(brp(1)**2+brp(2)**2+brp(3)**2)
    q = 0.5_dp*(ar+arp)

    IF (helium%periodic .AND. ((ar > 0.5_dp*helium%cell_size)&
                         .OR.(arp > 0.5_dp*helium%cell_size))) THEN
      v = 0.0_dp
      IF (arp>0.5_dp*helium%cell_size) THEN
        v = v + cut*helium_spline(tab(1,1)%spline_data,0.5_dp*helium%cell_size)
      ELSE
        v = v + helium_spline(tab(1,1)%spline_data,arp)
      END IF
      IF (ar>0.5_dp*helium%cell_size) THEN
        v = v + cut*helium_spline(tab(1,1)%spline_data,0.5_dp*helium%cell_size)
      ELSE
        v = v + helium_spline(tab(1,1)%spline_data,ar)
      END IF
      res = 0.5_dp*v
    ELSE
      ! end-point action (first term):
      v = 0.5_dp*(helium_spline(tab(1,1)%spline_data,ar)+helium_spline(tab(1,1)%spline_data,arp))
      DO i = 1, 3
        br(i) = br(i) - brp(i)
      END DO
      CALL helium_pbc( helium, br )
      s = br(1)**2+br(2)**2+br(3)**2
      z = (ar-arp)**2
      arp = 1.0_dp
      ! j=0 terms
      DO i = 2, SIZE(tab,1)
        arp = arp * s
        v = v + arp*helium_spline(tab(i,1)%spline_data,q)
      END DO
      ar = 1.0_dp
      DO j = 2, SIZE(tab,2)
        ar = ar * z
        arp = ar
        DO i = j, SIZE(tab,1)
          v = v + arp*helium_spline(tab(i,j)%spline_data,q)
          arp = arp * s
        END DO
      END DO
      res = v
    END IF
    RETURN
  END FUNCTION helium_eval_expansion

! *****************************************************************************
!> \brief ...
!> \param helium ...
! *****************************************************************************
  SUBROUTINE helium_update_transition_matrix(helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    INTEGER                                  :: b, c, i, j, k, m, n, nb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: lens, order
    INTEGER, DIMENSION(:), POINTER           :: perm
    INTEGER, DIMENSION(:, :), POINTER        :: nmatrix
    REAL(KIND=dp)                            :: f, q, t, v
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p
    REAL(KIND=dp), DIMENSION(3)              :: r
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ipmatrix, pmatrix, tmatrix
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos

   nb = helium%atoms
   !TODO: check allocation status
   ALLOCATE(p(2*nb),STAT=i)
   ALLOCATE(order(nb),STAT=i)
   ALLOCATE(lens(2*nb),STAT=i)
   b  = helium%beads-helium%bisection+1
   f = -0.5_dp/(helium%hb2m*helium%tau*helium%bisection)
   tmatrix => helium%tmatrix
   pmatrix => helium%pmatrix
   ipmatrix => helium%ipmatrix
   nmatrix => helium%nmatrix
   perm => helium%permutation
   pos => helium%pos
   DO i = 1, nb
      DO j = 1, nb
         v = 0.0_dp
         r(:) = pos(:,i,b) - pos(:,j,1)
         CALL helium_pbc( helium, r )
         v = v + r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
         pmatrix(i,j) = f*v
      END DO
      t = pmatrix(i,perm(i)) ! just some reference
      v = 0.0_dp
      DO j = 1, nb
         tmatrix(i,j) = EXP(pmatrix(i,j)-t)
         v = v + tmatrix(i,j)
      END DO
      ! normalize
      q = t+LOG(v)
      t = 1.0_dp/v
      DO j = 1, nb
         tmatrix(i,j) = tmatrix(i,j)*t
         ipmatrix(i,j) = 1.0_dp/tmatrix(i,j)
      END DO

      ! at this point we have:
      ! tmatrix(i,j) = exp(-f*(r_i^b - r_j^1)**2) normalized such
      !    that sum_j tmatrix(i,j) = 1.
      !    ( tmatrix(k1,k2) = t_{k1,k2} / h_{k1} of ceperly. )
      !    so tmatrix(i,j) is the probability to try to change a permutation
      !    with particle j (assuming particle i is already selected as well)
      ! ipmatrix(i,j) = 1.0/tmatrix(i,j)
      ! pmatrix(i,j) = log(tmatrix(i,j))  + some_offset(i)

      ! generate optimal search tree so we can select which particle j
      ! belongs to a given random_number as fast as possible.
      ! (the traditional approach would be to generate a table
      !  of cumulative probabilities and to search that table)
      ! so for example if we have:
      ! tmatrix(i,:) = ( 0.1 , 0.4 , 0.2 , 0.3 )
      ! traditionally we would build the running sum table:
      !  ( 0.1 , 0.5 , 0.7 , 1.0 ) and for a random number r
      ! would search this table for the lowest index larger than r
      ! (which would then be the particle index chosen by this random number)
      ! we build an optimal binary search tree instead, so here
      ! we would have:
      ! if ( r > 0.6 ) then take index 2,
      ! else if ( r > 0.3 ) then take index 4,
      ! else if ( r > 0.1 ) then take index 3 else index 1.
      ! the search tree is generated in tmatrix and nmatrix.
      ! tmatrix contains the decision values (0.6,0.3,0.1 in this case)
      ! and nmatrix contains the two branches (what to do if lower or higher)
      ! negative numbers in nmatrix mean take minus that index
      ! positive number means go down the tree to that next node, since we
      ! put the root of the tree at the end the arrays in the example would
      ! look like this:
      ! tmatrix(i,:) = ( 0.1 , 0.3 , 0.6 , arbitrary )
      ! namtrix(i,:) = ( -1 , -3 , 1 , -4 , 2 , -2 , arb. , arb. )
      !
      ! the way to generate this tree may not be the best, but the
      ! tree generation itself shouldn't be needed quite that often:
      !
      ! first sort values (with some variant of heap sort)

      DO j = 1, nb
         order(j)=j
         p(j) = tmatrix(i,j)
      END DO
      IF (nb > 1) THEN    ! if nb = 1 it is already sorted.
         k = nb/2+1
         c = nb
         DO
            IF (k > 1) THEN
               ! building up the heap:
               k = k - 1
               n = order(k)
               v = p(k)
            ELSE
               ! removing the top of the heap
               n = order(c)
               v = p(c)
               order(c)=order(1)
               p(c)=p(1)
               c = c - 1
               IF (c == 1) THEN
                  order(1)=n
                  p(1)=v
                  EXIT
               END IF
            END IF
            m = k
            j = 2*k
            ! restoring heap order between k and c
            DO
               IF (j > c) EXIT
               IF (j < c) THEN
                  IF (p(j)<p(j+1)) j = j + 1
               END IF
               IF (v >= p(j)) EXIT
               order(m) = order(j)
               p(m) = p(j)
               m = j
               j = 2*j
            END DO
            order(m) = n
            p(m) = v
         END DO
      END IF

      ! now:
      !    p(1:nb)     : tmatrix(i,1:nb) sorted in ascending order
      !    order(1:nb) : corresponding index: p(j) == tmatrix(i,order(j))
      !                                                       for all j

      ! merge sort with elements as we generate new interior search nodes
      ! by combining older elements/nodes

      ! first fill unused part of array with guard values:
      DO j = nb+1, 2*nb
         p(j)=2.0_dp
      END DO

      ! j   - head of leaf queue
      ! c+1 - head of node queue in p (c in lens)
      ! m+1 - tail of node queue in p (m in lens)
      c = nb+1
      j = 1
      DO m = nb+1, 2*nb-1
         ! get next smallest element
         IF (p(j) < p(c+1)) THEN
            v = p(j)
            lens(j) = m
            j = j + 1
         ELSE
            v = p(c+1)
            lens(c) = m
            c = c + 1
         END IF
         ! get the second next smallest element
         IF (p(j) < p(c+1)) THEN
            p(m+1) = v+p(j)
            lens(j) = m
            j = j + 1
         ELSE
            p(m+1) = v+p(c+1)
            lens(c) = m
            c = c + 1
         END IF
      END DO

      ! lens(:) now has the tree with lens(j) pointing to its parent
      ! the root of the tree is at 2*nb-1
      ! calculate the depth of each node in the tree now: (root = 0)

      lens(2*nb-1) = 0
      DO m = 2*nb-2, 1, -1
         lens(m) = lens(lens(m))+1
      END DO

      ! lens(:) now has the depths of the nodes/leafs

#if 0
      ! calculate average search depth (for information only)
      v = 0.0_dp
      DO j = 1, nb
         v = v + p(j)*lens(j)
      END DO
      PRINT *,"Expected number of comparisons with i=",i,v
#endif

      ! reset the nodes, for the canonical tree we just need the leaf info
      DO j = 1, nb
         lens(j+nb)=0
         p(j+nb) = 0.0_dp
      END DO

      ! build the canonical tree (number of decisions on average are
      ! the same to the tree we build above, but it has better caching behavior

      ! c head of leafs
      ! m head of interior nodes
      c = 1
      m = nb+1
      DO k = 1, 2*nb-2
         j = nb+1+(k-1)/2
         IF (lens(c)>lens(m+1)) THEN
            nmatrix(i,k) = -order(c)
            lens(j+1) = lens(c)-1
            v = p(c)
            c = c + 1
         ELSE
            nmatrix(i,k) = m-nb
            lens(j+1) = lens(m+1)-1
            v = p(m)
            m = m + 1
         END IF
         p(j) = p(j) + v
         IF (MOD(k,2)==1) tmatrix(i,j-nb)=v
      END DO

      ! now:
      !    nmatrix(i,2*j+1) left child of node j
      !    nmatrix(i,2*j+2) right child of node j
      !       children:
      !           negative : leaf with particle index == abs(value)
      !           positive : child node index
      !    p(j) weight of leaf j
      !    p(nb+j) weight of node j
      !    tmatrix(i,j) weight of left child of node j

      ! fix offsets for decision tree:

      p(nb-1)=0.0_dp
      DO m = nb-1,1,-1
         ! if right child is a node, set its offset and
         ! change its decision value
         IF (nmatrix(i,2*m)>0) THEN
            p(nmatrix(i,2*m)) = tmatrix(i,m)
            tmatrix(i,nmatrix(i,2*m)) = tmatrix(i,nmatrix(i,2*m))+tmatrix(i,m)
         END IF
         ! if left child is a node, set its offset and
         ! change its decision value
         IF (nmatrix(i,2*m-1)>0) THEN
            p(nmatrix(i,2*m-1)) = p(m)
            tmatrix(i,nmatrix(i,2*m-1)) = tmatrix(i,nmatrix(i,2*m-1)) + p(m)
         END IF
      END DO

      ! canonical optimal search tree done

#if 0
      !some test code, to check if it gives the right distribution
      DO k = 1, nb
         p(k)=1.0/ipmatrix(i,k)
      END DO
      lens(:)=0
      ! number of random numbers to generate:
      c = 1000000000
      DO j=1, c
         v = next_random_number(helium%rng_stream_uniform)
         ! walk down the search tree:
         k = nb-1
         DO
            IF (tmatrix(i,k) > v) THEN
               k = nmatrix(i,2*k-1)
            ELSE
               k = nmatrix(i,2*k)
            END IF
            IF (k<0) EXIT
         END DO
         k = -k
         ! increment the counter for this particle index
         lens(k) = lens(k)+1
      END DO
      ! search for maximum deviation from expectation value
      ! (relative to the expected variance)
      v = 0.0_dp
      k = -1
      DO j = 1, nb
         q = ABS((lens(j)-c*p(j))/SQRT(c*p(j)))
         !PRINT *,j,lens(j),c*p(j)
         IF (q > v) THEN
            v = q
            k = j
         END IF
         !PRINT *,lens(j),c*p(j),(lens(j)-c*p(j))/sqrt(c*p(j))
      END DO
      PRINT *,"MAXDEV:",k,lens(k),c*p(k),v
      !PRINT *,"TMAT:",tmatrix(i,:)
      !PRINT *,"NMAT:",nmatrix(i,:)
      !STOP
#endif
#if 0
      !additional test code:
      p(:) = -1.0_dp
      p(nb-1) = 0.0_dp
      p(2*nb-1) = 1.0_dp
      DO j = nb-1, 1, -1
         ! right child
         IF (nmatrix(i,2*j) > 0) THEN
            c = nmatrix(i,2*j)
            p(c) = tmatrix(i,j)
            p(c+nb) = p(j+nb)
         ELSE
            c = -nmatrix(i,2*j)
            !PRINT *,c,1.0/ipmatrix(i,c),p(j+nb)-tmatrix(i,j)
            IF (ABS(1.0/ipmatrix(i,c)-(p(j+nb)-tmatrix(i,j))) > &
                  10.0_dp*EPSILON(1.0_dp)) THEN
               PRINT *,"Probability mismatch for particle i->j",i,c
               PRINT *,"Got",p(j+nb)-tmatrix(i,j),"should be",1.0/ipmatrix(i,c)
               STOP
            END IF
         END IF
         ! left child
         IF (nmatrix(i,2*j-1) > 0) THEN
            c = nmatrix(i,2*j-1)
            p(c+nb) = tmatrix(i,j)
            p(c) = p(j)
         ELSE
            c = -nmatrix(i,2*j-1)
            !PRINT *,c,1.0/ipmatrix(i,c),tmatrix(i,j)-p(j)
            IF (ABS(1.0/ipmatrix(i,c)-(tmatrix(i,j)-p(j))) > &
                  10.0_dp*EPSILON(1.0_dp)) THEN
               PRINT *,"Probability mismatch for particle i->j",i,c
               PRINT *,"Got",tmatrix(i,j)-p(j),"should be",1.0/ipmatrix(i,c)
               STOP
            END IF
         END IF
      END DO
      PRINT *,"Probabilities ok"
#endif

   END DO

   ! initialize trial permutation with some identity permutation
   ! (should not be taken, but just in case it does we have something valid)

   helium%pweight = 0.0_dp
   t = next_random_number(helium%rng_stream_uniform)
   helium%ptable(1) = 1+INT(t*nb)
   helium%ptable(2) = -1

   ! recalculate inverse permutation table (just in case)
   DO i = 1, nb
      helium%iperm(perm(i))=i
   END DO

   ! clean up:
   DEALLOCATE(lens,STAT=i)
   DEALLOCATE(order,STAT=i)
   DEALLOCATE(p,STAT=i)

  RETURN
END SUBROUTINE helium_update_transition_matrix

! *****************************************************************************
!> \brief ...
!> \param spl ...
!> \param xx ...
!> \retval res ...
! *****************************************************************************
  FUNCTION helium_spline(spl, xx) RESULT(res)
    TYPE(spline_data_type), POINTER          :: spl
    REAL(KIND=dp), INTENT(IN)                :: xx
    REAL(KIND=dp)                            :: res

    REAL(KIND=dp)                            :: a, b

    IF (xx < spl%x1) THEN
       b = spl%invh*(xx-spl%x1)
       a = 1.0_dp-b
       res = a*spl%y(1)+b*(spl%y(2)-spl%y2(2)*spl%h26)
    ELSE IF (xx > spl%xn) THEN
       b = spl%invh*(xx-spl%xn)+1.0_dp
       a = 1.0_dp-b
       res = b*spl%y(spl%n)+a*(spl%y(spl%n-1)-spl%y2(spl%n-1)*spl%h26)
    ELSE
       res = spline_value(spl,xx)
    END IF
    RETURN
  END FUNCTION helium_spline

! *****************************************************************************
!> \brief ...
!> \param spl ...
!> \param xx ...
!> \retval res ...
! *****************************************************************************
  FUNCTION helium_spline_1(spl, xx) RESULT(res)
    TYPE(spline_data_type), POINTER          :: spl
    REAL(KIND=dp), INTENT(IN)                :: xx
    REAL(KIND=dp)                            :: res

    REAL(KIND=dp)                            :: dummy

    IF (xx < spl%x1) THEN
       res = spl%invh*(spl%y(2)-spl%y(1)-spl%y2(2)*spl%h26)
    ELSE IF (xx > spl%xn) THEN
       res = spl%invh*(spl%y(spl%n)-spl%y(spl%n-1)+spl%y2(spl%n-1)*spl%h26)
    ELSE
       dummy = spline_value(spl,xx,res)
    END IF
    RETURN
  END FUNCTION helium_spline_1

! *****************************************************************************
!> \brief  Return the distance <rij> between bead <ib> of atom <ia>
!>         and bead <jb> of atom <ja>.
!> \param helium ...
!> \param ia ...
!> \param ib ...
!> \param ja ...
!> \param jb ...
!> \retval rij ...
!> \date   2009-07-17
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_bead_rij(helium, ia, ib, ja, jb) RESULT(rij)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: ia, ib, ja, jb
    REAL(kind=dp)                            :: rij

    REAL(kind=dp)                            :: dx, dy, dz

    dx = helium%pos(1,ia,ib) - helium%pos(1,ja,jb)
    dy = helium%pos(2,ia,ib) - helium%pos(2,ja,jb)
    dz = helium%pos(3,ia,ib) - helium%pos(3,ja,jb)
    rij = SQRT(dx*dx+dy*dy+dz*dz)

    RETURN
  END FUNCTION helium_bead_rij

! *****************************************************************************
!> \brief  Given the atom number and permutation state return the cycle
!>         number the atom belongs to.
!> \param helium ...
!> \param atom_number ...
!> \param permutation ...
!> \retval cycle_number ...
!> \date   2009-07-21
!> \author Lukasz Walewski
!> \note   Cycles (or paths) are numbered from 1 to <num_cycles>, where 
!>         <num_cycles> is in the range of (1, <helium%atoms>). 
!>         if (num_cycles .EQ. 1) then all atoms belong to one cycle 
!>         if (num_cycles .EQ. helium%atoms) then there are no cycles of 
!>         length greater than 1 (i.e. no atoms are connected) 
! *****************************************************************************
  FUNCTION helium_cycle_number(helium, atom_number, permutation) RESULT(cycle_number)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: atom_number
    INTEGER, DIMENSION(:), POINTER           :: permutation
    INTEGER                                  :: cycle_number

    INTEGER                                  :: atom_idx, cycle_idx, &
                                                cycle_num, ia, ib, ic, &
                                                num_cycles
    INTEGER, DIMENSION(:), POINTER           :: cycle_index
    LOGICAL                                  :: break, new_cycle

    NULLIFY(cycle_index)
    cycle_index => helium%itmp_atoms_1d
    cycle_index(:) = 0

    num_cycles = 0
    break = .FALSE.
    DO ia = 1, helium%atoms
      ! this loop reaches its maximum iteration count when atom in question
      ! is the last one (i.e. when atom_number .EQ. helium%atoms)

      ! exit if we have found the cycle number for the atom in question
      IF (break) THEN
        EXIT
      END IF

      ! initialize current cycle index with the current atom
      cycle_idx = ia

      atom_idx = ia
      DO ib = 1, helium%atoms * helium%beads
        ! this loop reaches its maximum iteration count when all He atoms
        ! form one cycle (i.e. all beads belong to one path)

        ! proceed along the path
        atom_idx = permutation(atom_idx)

        IF (atom_idx .EQ. ia) THEN
          ! end of cycle detected (looped back to the first atom)

          ! check if this is a new cycle
          new_cycle = .TRUE.
          DO ic = 1, num_cycles
            IF (cycle_index(ic) .EQ. cycle_idx) THEN
              new_cycle = .FALSE.
            END IF
          END DO

          IF (new_cycle) THEN
            ! increase number of cycles and update the current cycle's index
            num_cycles = num_cycles + 1
            cycle_index(num_cycles) = cycle_idx
          END IF

          ! if this was the atom in question
          IF (ia .EQ. atom_number) THEN
            ! save the cycle index it belongs to
            cycle_num = cycle_idx

            ! exit the loop over atoms, we've found what we've been looking for
            break = .TRUE.
          END IF

          ! exit the loop over beads, there are no more (new) beads in this cycle
          EXIT
        END IF

        ! set the cycle index to the lowest atom index in this cycle
        IF (atom_idx .LT. cycle_idx) THEN
          cycle_idx = atom_idx
        END IF

      END DO

    END DO
    ! at this point we know the cycle index for atom <atom_number>
    ! but it is expressed as the atom number of the first atom in that cycle

    ! renumber cycle indices, so that they form a range (1, <num_cycles>)
    ! (don't do it actually - just return the requested <cycle_number>)
    cycle_number = 0
    DO ic = 1, num_cycles
      IF (cycle_index(ic) .EQ. cycle_num) THEN
        cycle_number = ic
        EXIT
      END IF
    END DO

    NULLIFY(cycle_index)

    RETURN
  END FUNCTION helium_cycle_number

! *****************************************************************************
!> \brief  Given the atom number and permutation state return the cycle
!>         number the atom belongs to.
!> \param helium ...
!> \param permutation ...
!> \retval num_cycles ...
!> \date   2009-07-21
!> \author Lukasz Walewski
!> \note   Cycles (or paths) are numbered from 1 to <num_cycles>, where 
!>         <num_cycles> is in the range of (1, <helium%atoms>). 
!>         if (num_cycles .EQ. 1) then all atoms belong to one cycle 
!>         if (num_cycles .EQ. helium%atoms) then there are no cycles of 
!>         length greater than 1 (i.e. no atoms are connected) 
! *****************************************************************************
  FUNCTION helium_ncycles(helium, permutation) RESULT(num_cycles)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, DIMENSION(:), POINTER           :: permutation
    INTEGER                                  :: num_cycles

    INTEGER                                  :: clen, curr_idx, cycle_idx, &
                                                ia, ia2, ib, ic, prev_idx
    INTEGER, DIMENSION(:), POINTER           :: cycle_index, cycle_length
    LOGICAL                                  :: new_cycle
    REAL(KIND=dp), DIMENSION(128, 3)         :: cycle_w
    REAL(KIND=dp), DIMENSION(3)              :: r, v

    NULLIFY(cycle_index)
    cycle_index => helium%itmp_atoms_1d
    cycle_index(:) = 0

    NULLIFY(cycle_length)
    cycle_length => helium%itmp_atoms_np_1d
    cycle_length(:) = 0

    num_cycles = 0
    DO ia = 1, helium%atoms

      ! initialize current cycle index with the current atom
      cycle_idx = ia

      curr_idx = ia
      clen = 0
      v(:) = 0.0_dp
      DO ia2 = 1, helium%atoms

        ! calc winding number contribution from the current atom
        ! contribution comming from the last and the first bead
        r(:) = helium%pos(:,curr_idx,helium%beads) - &
          helium%pos(:,helium%permutation(curr_idx),1)
        CALL helium_pbc( helium, r )
        v(:) = v(:) + r(:)
        ! sum of contributions from the rest of bead pairs
        DO ib = 2, helium%beads
          r(:) = helium%pos(:,curr_idx,ib-1) - helium%pos(:,curr_idx,ib)
          CALL helium_pbc( helium, r )
          v(:) = v(:) + r(:)
        END DO

        clen = clen + 1

        ! proceed along the path
        prev_idx = curr_idx
        curr_idx = permutation(curr_idx)

        IF (curr_idx .EQ. ia) THEN
          ! end of cycle detected (looped back to the first atom)

          ! check if this is a new cycle
          new_cycle = .TRUE.
          DO ic = 1, num_cycles
            IF (cycle_index(ic) .EQ. cycle_idx) THEN
              new_cycle = .FALSE.
            END IF
          END DO

          IF (new_cycle) THEN
            ! increase number of cycles and update the current cycle's index
            num_cycles = num_cycles + 1
            cycle_index(num_cycles) = cycle_idx
            cycle_length(num_cycles) = clen

    ! calc winding number (cell_m_inv is in Bohr^-1)
    cycle_w(num_cycles,1) = &
      helium%cell_m_inv(1,1) * v(1) + &
      helium%cell_m_inv(1,2) * v(2) + &
      helium%cell_m_inv(1,3) * v(3);
    cycle_w(num_cycles,2) = &
      helium%cell_m_inv(2,1) * v(1) + &
      helium%cell_m_inv(2,2) * v(2) + &
      helium%cell_m_inv(2,3) * v(3);
    cycle_w(num_cycles,3) = &
      helium%cell_m_inv(3,1) * v(1) + &
      helium%cell_m_inv(3,2) * v(2) + &
      helium%cell_m_inv(3,3) * v(3);
          END IF

          ! exit the loop over beads, there are no more (new) beads in this cycle
          EXIT
        END IF

        ! set the cycle index to the lowest atom index in this cycle
        IF (curr_idx .LT. cycle_idx) THEN
          cycle_idx = curr_idx
        END IF

      END DO

    END DO

!do ic = 1, num_cycles
!  print *, ic, cycle_index(ic), cycle_length(ic), cycle_w(ic,:)
!end do

    cycle_index(:) = 0
    NULLIFY(cycle_index)

    cycle_length(:) = 0
    NULLIFY(cycle_length)

    RETURN
  END FUNCTION helium_ncycles


! ***************************************************************************
!> \brief  Given the atom number and permutation state return the length of
!>         the path this atom belongs to.
!> \param helium ...
!> \param atom_number ...
!> \param permutation ...
!> \retval path_length ...
!> \date   2009-10-07
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_path_length(helium, atom_number, permutation) RESULT(path_length)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: atom_number
    INTEGER, DIMENSION(:), POINTER           :: permutation
    INTEGER                                  :: path_length

    INTEGER                                  :: atom_idx, ia
    LOGICAL                                  :: path_end_reached

    atom_idx = atom_number
    path_length = 0
    path_end_reached = .FALSE.
    DO ia = 1, helium%atoms
      path_length = path_length + 1
      atom_idx = permutation(atom_idx)
      IF ( atom_idx .EQ. atom_number ) THEN
        path_end_reached = .TRUE.
        EXIT
      END IF
    END DO

    IF ( .NOT. path_end_reached ) THEN
      path_length = -1
    END IF

    RETURN
  END FUNCTION helium_path_length

! ***************************************************************************
!> \brief  Given the atom number and permutation state returns .TRUE. if the
!>         atom belongs to the winding path, .FASLE. otherwise.
!> \param helium ...
!> \param atom_number ...
!> \param pos ...
!> \param permutation ...
!> \retval is_winding ...
!> \date   2010-09-21
!> \author Lukasz Walewski
!> \note   The path winds around the periodic box if it's widing number is
!> \note   greater than 0. This function calculates the winding number of the
!> \note   path the atom in question belongs to by following the permutation
!> \note   cycle; it exits from the path if the starting atom number is
!> \note   reached again.
! *****************************************************************************
  FUNCTION helium_is_winding(helium, atom_number, pos, permutation) RESULT(is_winding)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: atom_number
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos
    INTEGER, DIMENSION(:), POINTER           :: permutation
    LOGICAL                                  :: is_winding

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_is_winding', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: atom_idx, ia, ib, path_length
    LOGICAL                                  :: path_end_reached
    REAL(KIND=dp), DIMENSION(3)              :: r, v, w

    is_winding = .FALSE.
    v(:) = 0.0_dp
    atom_idx = atom_number
    path_length = 0
    path_end_reached = .FALSE.
    ! traverse the path the given atom belongs to
    DO ia = 1, helium%atoms
      ! contribution comming from the last and the first bead
      r(:) = pos(:,atom_idx,helium%beads) - &
             pos(:,permutation(atom_idx),1)
      CALL helium_pbc( helium, r )
      v(:) = v(:) + r(:)
      ! sum of contributions from the rest of bead pairs
      DO ib = 2, helium%beads
        r(:) = pos(:,atom_idx,ib-1) - pos(:,atom_idx,ib)
        CALL helium_pbc( helium, r )
        v(:) = v(:) + r(:)
      END DO
      ! jump to the next atom along the path
      path_length = path_length + 1
      atom_idx = permutation(atom_idx)
      IF ( atom_idx .EQ. atom_number ) THEN
        path_end_reached = .TRUE.
        EXIT
      END IF
    END DO

    stmp = ""
    WRITE(stmp,*) atom_number
    err_str = "Permutation path of atom " // TRIM(ADJUSTL(stmp)) // &
              " is not cyclic."
    CALL cp_assert(path_end_reached,cp_failure_level,&
         cp_assertion_failed, routineP, err_str)

    ! calc winding number (cell_m_inv is in Bohr^-1)
    w(1) = &
      helium%cell_m_inv(1,1) * v(1) + &
      helium%cell_m_inv(1,2) * v(2) + &
      helium%cell_m_inv(1,3) * v(3);
    w(2) = &
      helium%cell_m_inv(2,1) * v(1) + &
      helium%cell_m_inv(2,2) * v(2) + &
      helium%cell_m_inv(2,3) * v(3);
    w(3) = &
      helium%cell_m_inv(3,1) * v(1) + &
      helium%cell_m_inv(3,2) * v(2) + &
      helium%cell_m_inv(3,3) * v(3);

    IF ( w(1)*w(1)+w(2)*w(2)+w(3)*w(3) .GT. 0.5_dp ) THEN
      is_winding = .TRUE.
    END IF

    RETURN
  END FUNCTION helium_is_winding

! ***************************************************************************
!> \brief  Given the atom number and permutation state returns the
!>         contribution of this atom to the total superfluid density W^2/n.
!>         W^2 - squared winding number of the path the atom belongs to
!>         n   - number of atoms that belong to this path
!> \param helium ...
!> \param atom_number ...
!> \retval res ...
!> \date   2011-06-16
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_sdensity_part(helium, atom_number) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: atom_number
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_sdensity_part', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: atom_idx, ia, ib, path_length
    LOGICAL                                  :: path_end_reached
    REAL(KIND=dp)                            :: wsq
    REAL(KIND=dp), DIMENSION(3)              :: r, v, w

    v(:) = 0.0_dp
    atom_idx = atom_number
    path_length = 0
    path_end_reached = .FALSE.
    ! traverse the path given atom belongs to
    DO ia = 1, helium%atoms
      ! contribution comming from the last and the first bead
      r(:) = helium%pos(:,atom_idx,helium%beads) - &
             helium%pos(:,helium%permutation(atom_idx),1)
      CALL helium_pbc( helium, r )
      v(:) = v(:) + r(:)
      ! sum of contributions from the rest of bead pairs
      DO ib = 2, helium%beads
        r(:) = helium%pos(:,atom_idx,ib-1) - helium%pos(:,atom_idx,ib)
        CALL helium_pbc( helium, r )
        v(:) = v(:) + r(:)
      END DO
      ! jump to the next atom along the path
      path_length = path_length + 1
      atom_idx = helium%permutation(atom_idx)
      IF ( atom_idx .EQ. atom_number ) THEN
        path_end_reached = .TRUE.
        EXIT
      END IF
    END DO

    stmp = ""
    WRITE(stmp,*) atom_number
    err_str = "Permutation path of atom " // TRIM(ADJUSTL(stmp)) // &
              " is not cyclic."
    CALL cp_assert(path_end_reached,cp_failure_level,&
         cp_assertion_failed, routineP, err_str)

    ! calc winding number (cell_m_inv is in Bohr^-1)
    w(1) = &
      helium%cell_m_inv(1,1) * v(1) + &
      helium%cell_m_inv(1,2) * v(2) + &
      helium%cell_m_inv(1,3) * v(3);
    w(2) = &
      helium%cell_m_inv(2,1) * v(1) + &
      helium%cell_m_inv(2,2) * v(2) + &
      helium%cell_m_inv(2,3) * v(3);
    w(3) = &
      helium%cell_m_inv(3,1) * v(1) + &
      helium%cell_m_inv(3,2) * v(2) + &
      helium%cell_m_inv(3,3) * v(3);

    wsq = w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
    IF ( ABS(wsq) .LT. EPSILON(0.0_dp) ) wsq = 0.0_dp

    res = wsq / path_length

    RETURN
  END FUNCTION helium_sdensity_part


! ***************************************************************************
!> \brief  Given the permutation state assign cycle lengths to all He atoms.
!> \param helium ...
!> \date   2011-07-06
!> \author Lukasz Walewski
!> \note   The helium%atom_plength array is filled with cycle lengths,
!>         each atom gets the length of the permutation cycle it belongs to.
! *****************************************************************************
  SUBROUTINE helium_calc_atom_cycle_length( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_calc_atom_cycle_length', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str
    INTEGER                                  :: clen, curr_idx, handle, ia, &
                                                start_idx
    INTEGER, DIMENSION(:), POINTER           :: atoms_in_cycle
    LOGICAL                                  :: atoms_left, path_end_reached
    LOGICAL, DIMENSION(:), POINTER           :: atom_was_used

    CALL timeset(routineN,handle)

    NULLIFY(atoms_in_cycle)
    atoms_in_cycle => helium%itmp_atoms_1d
    atoms_in_cycle(:) = 0

    NULLIFY(atom_was_used)
    atom_was_used => helium%ltmp_atoms_1d
    atom_was_used(:) = .FALSE.

    helium%atom_plength(:) = 0

    start_idx = 1
    DO
      clen = 0
      path_end_reached = .FALSE.
      curr_idx = start_idx
      DO ia = 1, helium%atoms
        clen = clen + 1
        atoms_in_cycle(clen) = curr_idx
        atom_was_used(curr_idx) = .TRUE.

        ! follow to the next atom in the cycle
        curr_idx = helium%permutation(curr_idx)
        IF ( curr_idx .EQ. start_idx ) THEN
          path_end_reached = .TRUE.
          EXIT
        END IF
      END DO
      err_str = "Permutation path is not cyclic."
      CALL cp_assert(path_end_reached,cp_failure_level,&
           cp_assertion_failed, routineP, err_str)

      ! assign the cycle length to the collected atoms
      DO ia = 1, clen
        helium%atom_plength(atoms_in_cycle(ia)) = clen
      END DO

      ! go to the next "unused" atom
      atoms_left = .FALSE.
      DO ia = 1, helium%atoms
        IF ( .NOT. atom_was_used(ia) ) THEN
          start_idx = ia
          atoms_left = .TRUE.
          EXIT
        END IF
      END DO

      IF (.NOT. atoms_left) EXIT
    END DO

    atoms_in_cycle(:) = 0
    NULLIFY(atoms_in_cycle)
    atom_was_used(:) = .FALSE.
    NULLIFY(atom_was_used)

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_calc_atom_cycle_length


END MODULE helium_common
