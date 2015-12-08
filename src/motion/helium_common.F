!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Independent helium subroutines shared with other modules
!> \author Lukasz Walewski
!> \date   2009-07-14
!> \note   Avoiding circular deps: do not USE any other helium_* modules here.
! *****************************************************************************
MODULE helium_common

  USE helium_types,                    ONLY: helium_solvent_type,&
                                             int_arr_ptr,&
                                             rho_atom_number,&
                                             rho_moment_of_inertia,&
                                             rho_num,&
                                             rho_projected_area,&
                                             rho_winding_cycle,&
                                             rho_winding_number
  USE input_constants,                 ONLY: denominator_inertia,&
                                             denominator_natoms,&
                                             denominator_rperp2,&
                                             helium_cell_shape_cube,&
                                             helium_cell_shape_octahedron
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE parallel_rng_types,              ONLY: next_random_number
  USE physcon,                         ONLY: angstrom,&
                                             bohr
  USE pint_public,                     ONLY: pint_calc_centroid,&
                                             pint_com_pos
  USE pint_types,                      ONLY: pint_env_type
  USE splines_methods,                 ONLY: spline_value
  USE splines_types,                   ONLY: spline_data_p_type,&
                                             spline_data_type
#include "../base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_common'

  PUBLIC :: helium_pbc
  PUBLIC :: helium_boxmean_3d
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
  PUBLIC :: helium_cycle_of
  PUBLIC :: helium_total_winding_number
  PUBLIC :: helium_total_projected_area
  PUBLIC :: helium_link_projected_area
  PUBLIC :: helium_cycle_projected_area
  PUBLIC :: helium_cycle_projected_area_pbc
  PUBLIC :: helium_total_moment_of_inertia
  PUBLIC :: helium_total_moment_of_inertia_linkwise
  PUBLIC :: helium_com
  PUBLIC :: helium_total_projected_area_linkwise
  PUBLIC :: helium_link_vector
  PUBLIC :: helium_calc_cycles
  PUBLIC :: helium_update_coord_system
  PUBLIC :: helium_set_rdf_coord_system

  CONTAINS


! ***************************************************************************
!> \brief  General PBC routine for helium.
!> \param helium ...
!> \param r ...
!> \param enforce ...
!> \date   2009-09-30
!> \author Lukasz Walewski
!> \note  Check wheather PBC should be applied, if yes call low level 
!>         routine according to the unit cell shape. 
! *****************************************************************************
  SUBROUTINE helium_pbc( helium, r, enforce )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r
    LOGICAL, OPTIONAL                        :: enforce

    IF (.NOT. helium%periodic) THEN
      IF (.NOT. PRESENT(enforce)) THEN
        RETURN
      END IF
    END IF

    SELECT CASE (helium%cell_shape)
      CASE (helium_cell_shape_cube)
        CALL helium_pbc_cube( helium, r )
      CASE (helium_cell_shape_octahedron)
        CALL helium_pbc_trocta_opt( helium, r )
      CASE DEFAULT
        CPABORT("we should never get here")
    END SELECT

    RETURN
  END SUBROUTINE helium_pbc


! ***************************************************************************
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
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    REAL(KIND=dp)                            :: s

! x coordinate

      s = helium%cell_size_inv * r(1)
      IF ( s > 0.5_dp ) THEN
        s = s - REAL(INT(s+0.5_dp),dp)
      ELSEIF ( s < -0.5_dp ) THEN
        s = s - REAL(INT(s-0.5_dp),dp)
      END IF
      r(1) = s * helium%cell_size

      ! y coordinate
      s = helium%cell_size_inv * r(2)
      IF ( s > 0.5_dp ) THEN
        s = s - REAL(INT(s+0.5_dp),dp)
      ELSEIF ( s < -0.5_dp ) THEN
        s = s - REAL(INT(s-0.5_dp),dp)
      END IF
      r(2) = s * helium%cell_size

      ! z coordinate
      s = helium%cell_size_inv * r(3)
      IF ( s > 0.5_dp ) THEN
        s = s - REAL(INT(s+0.5_dp),dp)
      ELSEIF ( s < -0.5_dp ) THEN
        s = s - REAL(INT(s-0.5_dp),dp)
      END IF
      r(3) = s * helium%cell_size

    RETURN
  END SUBROUTINE helium_pbc_cube


#if 0

! ***************************************************************************
!> \brief  Apply PBC within truncated octahedral unit cell.
!> \param helium ...
!> \param r ...
!> \date   2009-10-02
!> \author Lukasz Walewski
!> \note   Original Allen & Tildesley routine adapted for our helium code.
! *****************************************************************************
  SUBROUTINE helium_pbc_trocta( helium, r )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    REAL(KIND=dp)                            :: corr, r75, rx, ry, rz

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

#endif


! ***************************************************************************
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
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: r

    REAL(KIND=dp)                            :: cell_size, cell_size_inv, &
                                                corr, rx, ry, rz, sx, sy, sz

    cell_size = helium%cell_size
    cell_size_inv = helium%cell_size_inv

    rx = r(1) * cell_size_inv
    IF ( rx > 0.5_dp ) THEN
       rx = rx - REAL(INT(rx+0.5_dp),dp)
    ELSEIF ( rx < -0.5_dp ) THEN
       rx = rx - REAL(INT(rx-0.5_dp),dp)
    END IF

    ry = r(2) * cell_size_inv
    IF ( ry > 0.5_dp ) THEN
       ry = ry - REAL(INT(ry+0.5_dp),dp)
    ELSEIF ( ry < -0.5_dp ) THEN
       ry = ry - REAL(INT(ry-0.5_dp),dp)
    END IF

    rz = r(3) * cell_size_inv
    IF ( rz > 0.5_dp ) THEN
       rz = rz - REAL(INT(rz+0.5_dp),dp)
    ELSEIF ( rz < -0.5_dp ) THEN
       rz = rz - REAL(INT(rz-0.5_dp),dp)
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


! ***************************************************************************
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
!> \brief  Given the permutation state assign cycle lengths to all He atoms.
!> \param helium ...
!> \date   2011-07-06
!> \author Lukasz Walewski
!> \note   The helium%atom_plength array is filled with cycle lengths, 
!>         each atom gets the length of the permutation cycle it belongs to. 
! *****************************************************************************
  SUBROUTINE helium_calc_atom_cycle_length( helium )

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
      IF (.NOT.path_end_reached) THEN
        CPABORT(err_str)
      END IF

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


! ***************************************************************************
!> \brief  Decompose a permutation into cycles
!> \param permutation ...
!> \retval cycles ...
!> \date   2013-12-11
!> \author Lukasz Walewski
!> \note   Given a permutation return a pointer to an array of pointers, 
!>         with each element pointing to a cycle of this permutation. 
!> \note   This function allocates memory and returns a pointer to it,
!>         do not forget to deallocate once finished with the results.
! *****************************************************************************
  FUNCTION helium_calc_cycles(permutation) RESULT(cycles)

    INTEGER, DIMENSION(:), POINTER           :: permutation
    TYPE(int_arr_ptr), DIMENSION(:), POINTER :: cycles

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_cycles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: curat, ncycl, nsize, nused
    INTEGER, DIMENSION(:), POINTER           :: cur_cycle, used_indices
    TYPE(int_arr_ptr), DIMENSION(:), POINTER :: my_cycles

    nsize = SIZE(permutation)

    ! most pesimistic scenario: no cycles longer than 1
    ALLOCATE(my_cycles(nsize))

    ! go over the permutation and identify cycles
    curat = 1
    nused = 0
    ncycl = 0
    DO WHILE (curat .LE. nsize)

      ! get the permutation cycle the current atom belongs to
      cur_cycle => helium_cycle_of(curat, permutation)

      ! include the current cycle in the pool of "used" indices
      nused = nused + SIZE(cur_cycle)
      CALL reallocate(used_indices,1,nused)
      used_indices(nused-SIZE(cur_cycle)+1:nused) = cur_cycle(1:SIZE(cur_cycle))

      ! store the pointer to the current cycle
      ncycl = ncycl + 1
      my_cycles(ncycl)%iap => cur_cycle

      ! free the local pointer
      NULLIFY(cur_cycle)

      ! try to increment the current atom index
      DO WHILE (ANY(used_indices .EQ. curat))
        curat = curat + 1
      END DO

    END DO

    DEALLOCATE(used_indices)
    NULLIFY(used_indices)

    ! assign the result
    ALLOCATE(cycles(ncycl))
    cycles(1:ncycl) = my_cycles(1:ncycl)

    DEALLOCATE(my_cycles)
    NULLIFY(my_cycles)

    RETURN
  END FUNCTION helium_calc_cycles


! ***************************************************************************
!> \brief  Calculate helium density distribution functions and store them
!>         in helium%rho_inst
!> \param helium ...
!> \date   2011-06-14
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_calc_rho(helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_rho', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msgstr
    INTEGER                                  :: aa, bx, by, bz, handle, ia, &
                                                ib, ic, id, ii, ir, &
                                                n_out_of_range, nbin
    INTEGER, DIMENSION(3)                    :: nw
    INTEGER, DIMENSION(:), POINTER           :: cycl_len
    LOGICAL                                  :: ltmp1, ltmp2, ltmp3
    REAL(KIND=dp)                            :: invd3r, invdr, invp, rtmp
    REAL(KIND=dp), DIMENSION(3)              :: maxr_half, r, vlink, vtotal, &
                                                wn
    TYPE(int_arr_ptr), DIMENSION(:), POINTER :: cycles

    CALL timeset(routineN,handle)

    maxr_half(:) = helium%rho_maxr / 2.0_dp
    invdr = 1.0_dp / helium%rho_delr
    invd3r = invdr * invdr * invdr
    invp = 1.0_dp / REAL(helium%beads,dp)
    nbin = helium%rho_nbin
    ! assign the cycle lengths to all the atoms
    CALL helium_calc_atom_cycle_length( helium )
    NULLIFY(cycl_len)
    cycl_len => helium%atom_plength
    DO ir = 1, rho_num ! loop over densities ---

      IF (.NOT.helium%rho_property(ir)%is_calculated) THEN
        ! skip densities that are not requested by the user
        CYCLE
      END IF

      SELECT CASE (ir)

        CASE (rho_atom_number)
          ii = helium%rho_property(ir)%component_index(1)
          helium%rho_incr(ii,:,:) = invp

        CASE (rho_projected_area)
          vtotal(:) = helium_total_projected_area(helium)
          DO ia = 1, helium%atoms
            DO ib = 1, helium%beads
              vlink(:) = helium_link_projected_area(helium,ia,ib)
              DO ic = 1, 3
                ii = helium%rho_property(ir)%component_index(ic)
                helium%rho_incr(ii,ia,ib) = vtotal(ic)*vlink(ic)*angstrom*angstrom*angstrom*angstrom
              END DO
            END DO
          END DO

!        CASE (rho_winding_number)
!          vtotal(:) = helium_total_winding_number(helium)
!          DO ia = 1, helium%atoms
!            DO ib = 1, helium%beads
!              vlink(:) = helium_link_winding_number(helium,ia,ib)
!              DO ic = 1, 3
!                ii = helium%rho_property(ir)%component_index(ic)
!                helium%rho_incr(ii,ia,ib) = vtotal(ic)*vlink(ic)*angstrom*angstrom
!              END DO
!            END DO
!          END DO

        CASE (rho_winding_number)
          vtotal(:) = helium_total_winding_number(helium)
          DO id = 1, 3
            ii = helium%rho_property(ir)%component_index(id)
            helium%rho_incr(ii,:,:) = 0.0_dp
          END DO
          NULLIFY(cycles)
          cycles => helium_calc_cycles(helium%permutation)
          DO ic = 1, SIZE(cycles)
            wn = helium_cycle_winding_number(helium,cycles(ic)%iap,helium%pos)
            DO ia = 1, SIZE(cycles(ic)%iap)
              aa = cycles(ic)%iap(ia)
              DO ib = 1, helium%beads
                vlink(:) = helium_link_winding_number(helium,aa,ib)
                DO id = 1, 3
                  IF (ABS(wn(id)) .GT. 100.0_dp*EPSILON(0.0_dp)) THEN
                    ii = helium%rho_property(ir)%component_index(id)
                    helium%rho_incr(ii,aa,ib) = vtotal(id)*vlink(id)*angstrom*angstrom
                  END IF
                END DO
              END DO
            END DO
          END DO
          DEALLOCATE(cycles)

        CASE (rho_winding_cycle)
          vtotal(:) = helium_total_winding_number(helium)
          DO id = 1, 3
            ii = helium%rho_property(ir)%component_index(id)
            helium%rho_incr(ii,:,:) = 0.0_dp
          END DO
          NULLIFY(cycles)
          cycles => helium_calc_cycles(helium%permutation)
          ! compute number of atoms in all winding cycles
          nw(:) = 0
          DO ic = 1, SIZE(cycles)
            wn = helium_cycle_winding_number(helium,cycles(ic)%iap,helium%pos)
            DO id = 1, 3
              IF (ABS(wn(id)) .GT. 100.0_dp*EPSILON(0.0_dp)) THEN
                nw(id) = nw(id) + SIZE(cycles(ic)%iap)
              END IF
            END DO
          END DO
          ! assign contributions to all beads of winding cycles only
          DO ic = 1, SIZE(cycles)
            wn = helium_cycle_winding_number(helium,cycles(ic)%iap,helium%pos)
            DO id = 1, 3
              IF (ABS(wn(id)) .GT. 100.0_dp*EPSILON(0.0_dp)) THEN
                DO ia = 1, SIZE(cycles(ic)%iap)
                  aa = cycles(ic)%iap(ia)
                  DO ib = 1, helium%beads
                    IF (nw(id) .GT. 0) THEN ! this test should always get passed
                      ii = helium%rho_property(ir)%component_index(id)
                      rtmp = invp / REAL(nw(id),dp)
                      helium%rho_incr(ii,aa,ib) = rtmp*vtotal(id)*vtotal(id)*angstrom*angstrom
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO
          DEALLOCATE(cycles)

        CASE (rho_moment_of_inertia)
          vtotal(:) = helium_total_moment_of_inertia(helium)
          DO ia = 1, helium%atoms
            DO ib = 1, helium%beads
              vlink(:) = helium_link_moment_of_inertia(helium,ia,ib)
              DO ic = 1, 3
                ii = helium%rho_property(ir)%component_index(ic)
                helium%rho_incr(ii,ia,ib) = vlink(ic)*angstrom*angstrom
              END DO
            END DO
          END DO

        CASE DEFAULT
          ! do nothing

      END SELECT

    END DO ! loop over densities ---

    n_out_of_range = 0
    helium%rho_inst(:,:,:,:) = 0.0_dp
    DO ia = 1, helium%atoms
    ! bin the bead positions of the current atom using the increments set above
      DO ib = 1, helium%beads
        ! map the current bead position to the corresponding voxel
        r(:) = helium%pos(:,ia,ib) - helium%center(:)
        ! enforce PBC even if this is a non-periodic calc to avoid density leakage
        CALL helium_pbc( helium, r, enforce=.TRUE. )
        ! set up bin indices (translate by L/2 to avoid non-positive array indices)
        bx = INT((r(1) + maxr_half(1)) * invdr) + 1
        by = INT((r(2) + maxr_half(2)) * invdr) + 1
        bz = INT((r(3) + maxr_half(3)) * invdr) + 1
        ! check that the resulting bin numbers are within array bounds
        ltmp1 = (0 .LT. bx) .AND. ( bx .LE. nbin )
        ltmp2 = (0 .LT. by) .AND. ( by .LE. nbin )
        ltmp3 = (0 .LT. bz) .AND. ( bz .LE. nbin )
        IF (ltmp1 .AND. ltmp2 .AND. ltmp3) THEN
          ! increment all the estimators (those that the current atom does not
          ! contribute to have increment incr(ic)==0)
          DO ic = 1, helium%rho_num_act
            helium%rho_inst(ic,bx,by,bz) = helium%rho_inst(ic,bx,by,bz) + helium%rho_incr(ic,ia,ib)
          END DO
        ELSE
          n_out_of_range = n_out_of_range + 1
        END IF
      END DO
    END DO
    ! scale by volume element
    helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) * invd3r

    ! stop if any bead fell out of the range
    ! since enforced PBC should have taken care of such leaks
    WRITE(msgstr,*) n_out_of_range
    msgstr = "Number of bead positions out of range: "//TRIM(ADJUSTL(msgstr))
    IF (n_out_of_range .GT. 0) THEN
      CPABORT(msgstr)
    END IF

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_calc_rho


#if 0
! ***************************************************************************
!> \brief  Normalize superfluid densities according to the input keyword
!>         HELIUM%SUPERFLUID_ESTIMATOR%DENOMINATOR
!> \param helium ...
!> \param rho ...
!> \date   2014-06-24
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_norm_rho(helium,rho)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: rho

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_norm_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ix, iy, iz, ndim
    REAL(KIND=dp)                            :: invatoms, rx, ry, rz
    REAL(KIND=dp), DIMENSION(3)              :: invmoit, invrperp, ro

    SELECT CASE(helium%supest_denominator)

      CASE (denominator_natoms)
        invatoms = 1.0_dp / REAL(helium%atoms,dp)
        rho(2,:,:,:) = rho(2,:,:,:) * invatoms
        rho(3,:,:,:) = rho(3,:,:,:) * invatoms
        rho(4,:,:,:) = rho(4,:,:,:) * invatoms

      CASE (denominator_inertia)
        invmoit(:) = REAL(helium%atoms,dp) / helium%mominer%ravr(:)
        rho(2,:,:,:) = rho(2,:,:,:) * invmoit(1)
        rho(3,:,:,:) = rho(3,:,:,:) * invmoit(2)
        rho(4,:,:,:) = rho(4,:,:,:) * invmoit(3)

      CASE (denominator_rperp2)
        ndim = helium%rho_nbin
        ro(:) = helium%center(:) - 0.5_dp*(helium%rho_maxr-helium%rho_delr)
        DO ix = 1, ndim
          DO iy = 1, ndim
            DO iz = 1, ndim
              rx = ro(1) + REAL(ix-1,dp) * helium%rho_delr
              ry = ro(2) + REAL(iy-1,dp) * helium%rho_delr
              rz = ro(3) + REAL(iz-1,dp) * helium%rho_delr
              invrperp(1) = 1.0_dp / (ry*ry+rz*rz)
              invrperp(2) = 1.0_dp / (rz*rz+rx*rx)
              invrperp(3) = 1.0_dp / (rx*rx+ry*ry)
              rho(2,ix,iy,iz) = rho(2,ix,iy,iz) * invrperp(1)
              rho(3,ix,iy,iz) = rho(3,ix,iy,iz) * invrperp(2)
              rho(4,ix,iy,iz) = rho(4,ix,iy,iz) * invrperp(3)
            END DO
          END DO
        END DO

      CASE DEFAULT
        ! do nothing

    END SELECT

    RETURN
  END SUBROUTINE helium_norm_rho
#endif


! ***************************************************************************
!> \brief  Calculate helium radial distribution functions wrt positions given
!>         by <centers> using the atomic weights given by <weights>.
!> \param helium ...
!> \param centers ...
!> \date   2009-07-22
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_calc_rdf (helium,centers)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(:), POINTER     :: centers

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_rdf', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msgstr
    INTEGER                                  :: bin, handle, ia, ib, ic, id, &
                                                n_out_of_range, nbin
    REAL(KIND=dp)                            :: const, invdr, invp, nideal, &
                                                pi, ri, rlower, rupper, sdens
    REAL(KIND=dp), DIMENSION(3)              :: r, r0
    REAL(KIND=dp), DIMENSION(:), POINTER     :: incr

    CALL timeset(routineN,handle)

    pi = 4.0_dp * ATAN(1.0_dp)
    invdr = 1.0_dp / helium%rdf_delr
    invp = 1.0_dp / helium%beads
    nbin = helium%rdf_nbin
    sdens = helium%wpref * invp * helium%atoms
    ALLOCATE(incr(helium%rdf_num))
    incr(:) = 0.0d0

    ! calculate the histogram of distances
    n_out_of_range = 0
    helium%rdf_inst(:,:,:) = 0.0_dp
    DO ic = 1, SIZE(centers) / 3

      r0(1) = centers(3*(ic-1)+1)
      r0(2) = centers(3*(ic-1)+2)
      r0(3) = centers(3*(ic-1)+3)
      DO ia = 1, helium%atoms

        ! set the increments for this atom
        incr(1) = invp

        DO ib = 1, helium%beads
          r(:) = helium%pos(:,ia,ib) - r0(:)
          CALL helium_pbc( helium, r )
          ri = SQRT( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
          bin = INT(ri * invdr) + 1
          IF ((0 .LT. bin) .AND. ( bin .LE. nbin )) THEN
            ! increment the RDF value for He atoms inside the r_6 sphere
            DO id = 1, helium%rdf_num
              helium%rdf_inst(id,bin,ic) = helium%rdf_inst(id,bin,ic) + incr(id)
            END DO
          ELSE
            !WRITE(msgstr,*) helium%center * angstrom
            !msgstr = "center = "//TRIM(ADJUSTL(msgstr))
            !CALL cp_error_message(cp_warning_level, routineP, msgstr, error)
            !WRITE(msgstr,*) r0 * angstrom
            !msgstr = "r0     = "//TRIM(ADJUSTL(msgstr))
            !CALL cp_error_message(cp_warning_level, routineP, msgstr, error)
            !WRITE(msgstr,*) helium%pos(:,ia,ib) * angstrom
            !msgstr = "pos    = "//TRIM(ADJUSTL(msgstr))
            !CALL cp_error_message(cp_warning_level, routineP, msgstr, error)
            !WRITE(msgstr,*) ri * angstrom
            !msgstr = "ri     = "//TRIM(ADJUSTL(msgstr))
            !CALL cp_error_message(cp_warning_level, routineP, msgstr, error)
            !WRITE(msgstr,*) ri * invdr
            !msgstr = "ri/dr  = "//TRIM(ADJUSTL(msgstr))
            !CALL cp_error_message(cp_warning_level, routineP, msgstr, error)
            !WRITE(msgstr,*) helium%current_step
            !msgstr = "step   = "//TRIM(ADJUSTL(msgstr))
            !CALL cp_error_message(cp_warning_level, routineP, msgstr, error)
            n_out_of_range = n_out_of_range + 1
          END IF
        END DO
      END DO
    END DO

    IF (.NOT. helium%periodic) THEN
      IF (n_out_of_range .GT. 0) THEN
        WRITE(msgstr,*) n_out_of_range
        msgstr = "Number of bead positions out of range: "//TRIM(ADJUSTL(msgstr))
        CPABORT(msgstr)
      END IF
    END IF
    ! for periodic case we intentionally truncate RDFs to spherical volume
    ! so we skip atoms in the corners

    ! normalize the histogram to get g(r)
    ! note: helium%density refers to the number of atoms, not the beads
    const = 4.0_dp * pi * helium%density / 3.0_dp
    DO bin = 1, helium%rdf_nbin
      rlower = REAL(bin-1,dp) * helium%rdf_delr
      rupper = rlower + helium%rdf_delr
      nideal = const * (rupper**3 - rlower**3)
      DO id = 1, helium%rdf_num
        helium%rdf_inst(id,bin,:) = helium%rdf_inst(id,bin,:) / nideal
      END DO
    END DO

    DEALLOCATE(incr)
    NULLIFY(incr)

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_calc_rdf


! ***************************************************************************
!> \brief  Calculate probability distribution of the permutation lengths
!> \param helium ...
!> \date   2010-06-07
!> \author Lukasz Walewski
!> \note   Valid permutation path length is an integer (1, NATOMS), number 
!>         of paths of a given length is calculated here and average over 
!>         inner loop iterations and helium environments is done in 
!>         helium_sample. 
! *****************************************************************************
  SUBROUTINE helium_calc_plength( helium )

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
!> \note   Positions of helium beads in helium%pos array are reorganized such 
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


! ***************************************************************************
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
        v = v + REAL(cut,dp)*helium_spline(tab(1,1)%spline_data,0.5_dp*helium%cell_size)
      ELSE
        v = v + helium_spline(tab(1,1)%spline_data,arp)
      END IF
      IF (ar>0.5_dp*helium%cell_size) THEN
        v = v + REAL(cut,dp)*helium_spline(tab(1,1)%spline_data,0.5_dp*helium%cell_size)
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
!> \author Harald Forbert
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
   ALLOCATE(p(2*nb))
   ALLOCATE(order(nb))
   ALLOCATE(lens(2*nb))
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
   DEALLOCATE(lens)
   DEALLOCATE(order)
   DEALLOCATE(p)

   RETURN
  END SUBROUTINE helium_update_transition_matrix


! *****************************************************************************
!> \brief ...
!> \param spl ...
!> \param xx ...
!> \retval res ...
!> \author Harald Forbert
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


! ***************************************************************************
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
    REAL(KIND=dp)                            :: rij

    REAL(KIND=dp)                            :: dx, dy, dz

    dx = helium%pos(1,ia,ib) - helium%pos(1,ja,jb)
    dy = helium%pos(2,ia,ib) - helium%pos(2,ja,jb)
    dz = helium%pos(3,ia,ib) - helium%pos(3,ja,jb)
    rij = SQRT(dx*dx+dy*dy+dz*dz)

    RETURN
  END FUNCTION helium_bead_rij


! ***************************************************************************
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
    LOGICAL                                  :: found, new_cycle

    NULLIFY(cycle_index)
    cycle_index => helium%itmp_atoms_1d
    cycle_index(:) = 0

    num_cycles = 0
    found = .FALSE.
    cycle_num = -1
    DO ia = 1, helium%atoms
      ! this loop reaches its maximum iteration count when atom in question
      ! is the last one (i.e. when atom_number .EQ. helium%atoms)

      ! exit if we have found the cycle number for the atom in question
      IF (found) THEN
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
            found = .TRUE.
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

!TODO make it cp2k way
    IF (.NOT. found) THEN
        CPWARN("helium_cycle_number: we are going to return -1, problems ahead!")
    END IF

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
!> \brief  Given an element of a permutation return the cycle it belongs to.
!> \param element ...
!> \param permutation ...
!> \retval CYCLE ...
!> \date   2013-12-10
!> \author Lukasz Walewski
!> \note   This function allocates memory and returns a pointer to it,
!>         do not forget to deallocate once finished with the results.
! *****************************************************************************
  FUNCTION helium_cycle_of(element, permutation) RESULT(CYCLE)

    INTEGER, INTENT(IN)                      :: element
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: permutation
    INTEGER, DIMENSION(:), POINTER           :: CYCLE

    INTEGER                                  :: ia, icur, len, nsize
    CHARACTER(len=*), PARAMETER :: routineN = 'helium_cycle_of', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: my_cycle
    LOGICAL                                  :: cycle_end_reached

    nsize = SIZE(permutation)

    ! maximum possible cycle length is the number of atoms
    NULLIFY(my_cycle)
    ALLOCATE(my_cycle(nsize))

    ! traverse the permutation table
    len = 0
    icur = element
    cycle_end_reached = .FALSE.
    DO ia = 1, nsize
      len = len + 1
      my_cycle(len) = icur
      icur = permutation(icur)
      IF ( icur .EQ. element ) THEN
        cycle_end_reached = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT. cycle_end_reached) THEN
      ! return null
      NULLIFY(CYCLE)
    ELSE
      ! assign the result
      ALLOCATE(CYCLE(len))
      CYCLE(1:len) = my_cycle(1:len)
    END IF

    ! clean up
    DEALLOCATE(my_cycle)
    NULLIFY(my_cycle)

    RETURN
  END FUNCTION helium_cycle_of


! ***************************************************************************
!> \brief  Return total winding number
!> \param helium ...
!> \retval wnum ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_total_winding_number(helium) RESULT(wnum)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: wnum

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_winding_number', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: rcur
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ri, rj

    wnum(:) = 0.0_dp
    DO ia = 1, helium%atoms
      ! sum of contributions from the rest of bead pairs
      DO ib = 1, helium%beads-1
        ri => helium%pos(:,ia,ib)
        rj => helium%pos(:,ia,ib+1)
        rcur(:) = ri(:) - rj(:)
        CALL helium_pbc(helium,rcur)
        wnum(:) = wnum(:) + rcur(:)
      END DO
      ! contribution from the last and the first bead
      ri => helium%pos(:,ia,helium%beads)
      rj => helium%pos(:,helium%permutation(ia),1)
      rcur(:) = ri(:) - rj(:)
      CALL helium_pbc(helium,rcur)
      wnum(:) = wnum(:) + rcur(:)
    END DO

  END FUNCTION helium_total_winding_number


! ***************************************************************************
!> \brief  Return link winding number
!> \param helium ...
!> \param ia ...
!> \param ib ...
!> \retval wnum ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_link_winding_number(helium,ia,ib) RESULT(wnum)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: wnum

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_link_winding_number', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ja1, ja2, jb1, jb2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r1, r2

    IF (ib .EQ. helium%beads) THEN
      ja1 = ia
      ja2 = helium%permutation(ia)
      jb1 = ib
      jb2 = 1
    ELSE
      ja1 = ia
      ja2 = ia
      jb1 = ib
      jb2 = ib+1
    END IF
    r1 => helium%pos(:,ja1,jb1)
    r2 => helium%pos(:,ja2,jb2)
    wnum(:) = r1(:) - r2(:)
    CALL helium_pbc(helium,wnum)

    RETURN
  END FUNCTION helium_link_winding_number


! ***************************************************************************
!> \brief  Return total winding number (sum over all links)
!> \param helium ...
!> \retval wnum ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_total_winding_number_linkwise(helium) RESULT(wnum)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: wnum

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_winding_number_linkwise', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib

    wnum(:) = 0.0_dp
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        wnum(:) = wnum(:) + helium_link_winding_number(helium,ia,ib)
      END DO
    END DO

    RETURN
  END FUNCTION helium_total_winding_number_linkwise


! ***************************************************************************
!> \brief  Return cycle winding number
!> \param helium ...
!> \param CYCLE ...
!> \param pos ...
!> \retval wnum ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_cycle_winding_number(helium,CYCLE,pos) RESULT(wnum)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, DIMENSION(:), POINTER           :: CYCLE
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos
    REAL(KIND=dp), DIMENSION(3)              :: wnum

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_cycle_winding_number', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ia, ib, nsize
    REAL(KIND=dp), DIMENSION(3)              :: rcur
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ri, rj

    nsize = SIZE(CYCLE)

    ! traverse the path
    wnum(:) = 0.0_dp
    DO ia = 1, nsize
      ! contributions from all bead pairs of the current atom
      DO ib = 1, helium%beads-1
        ri => pos(:,CYCLE(ia),ib)
        rj => pos(:,CYCLE(ia),ib+1)
        rcur(:) = ri(:) - rj(:)
        CALL helium_pbc( helium, rcur )
        wnum(:) = wnum(:) + rcur(:)
      END DO
      ! contribution from the last bead of the current atom
      ! and the first bead of the next atom
      i1 = CYCLE(ia)
      IF (ia .EQ. nsize) THEN
        i2 = CYCLE(1)
      ELSE
        i2 = CYCLE(ia+1)
      END IF
      ri => pos(:,i1,helium%beads)
      rj => pos(:,i2,1)
      rcur(:) = ri(:) - rj(:)
      CALL helium_pbc( helium, rcur )
      wnum(:) = wnum(:) + rcur(:)
    END DO

    RETURN
  END FUNCTION helium_cycle_winding_number


! ***************************************************************************
!> \brief  Return total winding number (sum over all cycles)
!> \param helium ...
!> \retval wnum ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_total_winding_number_cyclewise(helium) RESULT(wnum)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: wnum

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_winding_number_cyclewise', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic
    REAL(KIND=dp), DIMENSION(3)              :: wn
    TYPE(int_arr_ptr), DIMENSION(:), POINTER :: cycles

! decompose the current permutation state into permutation cycles

    NULLIFY(cycles)
    cycles => helium_calc_cycles(helium%permutation)

    wnum(:) = 0.0_dp
    DO ic = 1, SIZE(cycles)
      wn = helium_cycle_winding_number(helium,cycles(ic)%iap,helium%pos)
      wnum(:) = wnum(:) + wn(:)
    END DO

    DEALLOCATE(cycles)

    RETURN
  END FUNCTION helium_total_winding_number_cyclewise


! ***************************************************************************
!> \brief  Return total projected area
!> \param helium ...
!> \retval area ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_total_projected_area(helium) RESULT(area)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_projected_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: r1, r12, r2, rcur

    area(:) = 0.0_dp
    DO ia = 1, helium%atoms
      ! contributions from all links of the current atom
      DO ib = 1, helium%beads-1
        r1(:) = helium%pos(:,ia,ib)
        r2(:) = helium%pos(:,ia,ib+1)
        ! comment out for non-PBC version -->
        r12(:) = r2(:) - r1(:)
        CALL helium_pbc(helium,r1)
        CALL helium_pbc(helium,r12)
        r2(:) = r1(:) + r12(:)
        ! comment out for non-PBC version <--
        rcur(1) = r1(2)*r2(3)-r1(3)*r2(2)
        rcur(2) = r1(3)*r2(1)-r1(1)*r2(3)
        rcur(3) = r1(1)*r2(2)-r1(2)*r2(1)
        area(:) = area(:) + rcur(:)
      END DO
      ! contribution from the last bead of the current atom
      ! and the first bead of the next atom
      r1(:) = helium%pos(:,ia,helium%beads)
      r2(:) = helium%pos(:,helium%permutation(ia),1)
      ! comment out for non-PBC version -->
      r12(:) = r2(:) - r1(:)
      CALL helium_pbc(helium,r1)
      CALL helium_pbc(helium,r12)
      r2(:) = r1(:) + r12(:)
      ! comment out for non-PBC version <--
      rcur(1) = r1(2)*r2(3)-r1(3)*r2(2)
      rcur(2) = r1(3)*r2(1)-r1(1)*r2(3)
      rcur(3) = r1(1)*r2(2)-r1(2)*r2(1)
      area(:) = area(:) + rcur(:)
    END DO
    area(:) = 0.5_dp * area(:)

    RETURN
  END FUNCTION helium_total_projected_area


! ***************************************************************************
!> \brief  Return link projected area
!> \param helium ...
!> \param ia ...
!> \param ib ...
!> \retval area ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_link_projected_area(helium,ia,ib) RESULT(area)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_link_projected_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ja1, ja2, jb1, jb2
    REAL(KIND=dp), DIMENSION(3)              :: r1, r12, r2

    IF (ib .EQ. helium%beads) THEN
      ja1 = ia
      ja2 = helium%permutation(ia)
      jb1 = ib
      jb2 = 1
    ELSE
      ja1 = ia
      ja2 = ia
      jb1 = ib
      jb2 = ib+1
    END IF
    r1(:) = helium%pos(:,ja1,jb1)
    r2(:) = helium%pos(:,ja2,jb2)
    ! comment out for non-PBC version -->
    r12(:) = r2(:) - r1(:)
    CALL helium_pbc(helium,r1)
    CALL helium_pbc(helium,r12)
    r2(:) = r1(:) + r12(:)
    ! comment out for non-PBC version <--
    area(1) = r1(2)*r2(3)-r1(3)*r2(2)
    area(2) = r1(3)*r2(1)-r1(1)*r2(3)
    area(3) = r1(1)*r2(2)-r1(2)*r2(1)
    area(:) = 0.5_dp * area(:)

    RETURN
  END FUNCTION helium_link_projected_area


! ***************************************************************************
!> \brief  Return total projected area (sum over all links)
!> \param helium ...
!> \retval area ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_total_projected_area_linkwise(helium) RESULT(area)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: area

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_projected_area_linkwise', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib

    area(:) = 0.0_dp
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        area(:) = area(:) + helium_link_projected_area(helium,ia,ib)
      END DO
    END DO

  END FUNCTION helium_total_projected_area_linkwise


! ***************************************************************************
!> \brief  Return cycle projected area
!> \param helium ...
!> \param CYCLE ...
!> \retval area ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_cycle_projected_area(helium,CYCLE) RESULT(area)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, DIMENSION(:), POINTER           :: CYCLE
    REAL(KIND=dp), DIMENSION(3)              :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_cycle_projected_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ia, ib, nsize
    REAL(KIND=dp), DIMENSION(3)              :: rcur, rsum
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ri, rj

    nsize = SIZE(CYCLE)

    ! traverse the path
    rsum(:) = 0.0_dp
    DO ia = 1, nsize
      ! contributions from all bead pairs of the current atom
      DO ib = 1, helium%beads-1
        ri => helium%pos(:,CYCLE(ia),ib)
        rj => helium%pos(:,CYCLE(ia),ib+1)
        rcur(1) = ri(2)*rj(3)-ri(3)*rj(2)
        rcur(2) = ri(3)*rj(1)-ri(1)*rj(3)
        rcur(3) = ri(1)*rj(2)-ri(2)*rj(1)
        rsum(:) = rsum(:) + rcur(:)
      END DO
      ! contribution from the last bead of the current atom
      ! and the first bead of the next atom
      i1 = CYCLE(ia)
      IF (ia .EQ. nsize) THEN
        i2 = CYCLE(1)
      ELSE
        i2 = CYCLE(ia+1)
      END IF
      ri => helium%pos(:,i1,helium%beads)
      rj => helium%pos(:,i2,1)
      rcur(1) = ri(2)*rj(3)-ri(3)*rj(2)
      rcur(2) = ri(3)*rj(1)-ri(1)*rj(3)
      rcur(3) = ri(1)*rj(2)-ri(2)*rj(1)
      rsum(:) = rsum(:) + rcur(:)
    END DO
    area(:) = 0.5_dp * rsum(:)

    RETURN
  END FUNCTION helium_cycle_projected_area


! ***************************************************************************
!> \brief  Return cycle projected area (sum over all links)
!> \param helium ...
!> \param CYCLE ...
!> \retval area ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_cycle_projected_area_pbc(helium,CYCLE) RESULT(area)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, DIMENSION(:), POINTER           :: CYCLE
    REAL(KIND=dp), DIMENSION(3)              :: area

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_cycle_projected_area_pbc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ia, ib, nsize
    REAL(KIND=dp), DIMENSION(3)              :: r1, r12, r2, rcur

    nsize = SIZE(CYCLE)

    ! traverse the path
    area(:) = 0.0_dp
    DO ia = 1, nsize
      ! contributions from all bead pairs of the current atom
      DO ib = 1, helium%beads-1
        r1(:) = helium%pos(:,CYCLE(ia),ib)
        r2(:) = helium%pos(:,CYCLE(ia),ib+1)
        r12(:) = r2(:) - r1(:)
        CALL helium_pbc(helium,r1)
        CALL helium_pbc(helium,r12)
        r2(:) = r1(:) + r12(:)
        rcur(1) = r1(2)*r2(3)-r1(3)*r2(2)
        rcur(2) = r1(3)*r2(1)-r1(1)*r2(3)
        rcur(3) = r1(1)*r2(2)-r1(2)*r2(1)
        area(:) = area(:) + rcur(:)
      END DO
      ! contribution from the last bead of the current atom
      ! and the first bead of the next atom
      i1 = CYCLE(ia)
      IF (ia .EQ. nsize) THEN
        i2 = CYCLE(1)
      ELSE
        i2 = CYCLE(ia+1)
      END IF
      r1(:) = helium%pos(:,i1,helium%beads)
      r2(:) = helium%pos(:,i2,1)
      r12(:) = r2(:) - r1(:)
      CALL helium_pbc(helium,r1)
      CALL helium_pbc(helium,r12)
      r2(:) = r1(:) + r12(:)
      rcur(1) = r1(2)*r2(3)-r1(3)*r2(2)
      rcur(2) = r1(3)*r2(1)-r1(1)*r2(3)
      rcur(3) = r1(1)*r2(2)-r1(2)*r2(1)
      area(:) = area(:) + rcur(:)
    END DO
    area(:) = 0.5_dp * area(:)

    RETURN
  END FUNCTION helium_cycle_projected_area_pbc


! ***************************************************************************
!> \brief  Return total projected area (sum over all cycles)
!> \param helium ...
!> \retval area ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_total_projected_area_cyclewise(helium) RESULT(area)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: area

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_projected_area_cyclewise', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic
    REAL(KIND=dp), DIMENSION(3)              :: pa
    TYPE(int_arr_ptr), DIMENSION(:), POINTER :: cycles

! decompose the current permutation state into permutation cycles

    NULLIFY(cycles)
    cycles => helium_calc_cycles(helium%permutation)

    area(:) = 0.0_dp
    DO ic = 1, SIZE(cycles)
      pa = helium_cycle_projected_area(helium,cycles(ic)%iap)
      area(:) = area(:) + pa(:)
    END DO

    RETURN
  END FUNCTION helium_total_projected_area_cyclewise


! ***************************************************************************
!> \brief  Return total moment of inertia divided by m_He
!> \param helium ...
!> \retval moit ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_total_moment_of_inertia(helium) RESULT(moit)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: moit

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_moment_of_inertia', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: com, r1, r12, r2, rcur

    com(:) = helium%center(:)

    moit(:) = 0.0_dp
    DO ia = 1, helium%atoms
      ! contributions from all the links of the current atom
      DO ib = 1, helium%beads-1
        r1(:) = helium%pos(:,ia,ib) - com(:)
        r2(:) = helium%pos(:,ia,ib+1) - com(:)
        ! comment out for non-PBC version -->
        r12(:) = r2(:) - r1(:)
        CALL helium_pbc(helium,r1)
        CALL helium_pbc(helium,r12)
        r2(:) = r1(:) + r12(:)
        ! comment out for non-PBC version <--
        rcur(1) = r1(2)*r2(2)+r1(3)*r2(3)
        rcur(2) = r1(3)*r2(3)+r1(1)*r2(1)
        rcur(3) = r1(1)*r2(1)+r1(2)*r2(2)
        moit(:) = moit(:) + rcur(:)
      END DO
      ! contribution from the last bead of the current atom
      ! and the first bead of the next atom
      r1(:) = helium%pos(:,ia,helium%beads) - com(:)
      r2(:) = helium%pos(:,helium%permutation(ia),1) - com(:)
      ! comment out for non-PBC version -->
      r12(:) = r2(:) - r1(:)
      CALL helium_pbc(helium,r1)
      CALL helium_pbc(helium,r12)
      r2(:) = r1(:) + r12(:)
      ! comment out for non-PBC version <--
      rcur(1) = r1(2)*r2(2)+r1(3)*r2(3)
      rcur(2) = r1(3)*r2(3)+r1(1)*r2(1)
      rcur(3) = r1(1)*r2(1)+r1(2)*r2(2)
      moit(:) = moit(:) + rcur(:)
    END DO
    moit(:) = moit(:) / helium%beads

    RETURN
  END FUNCTION helium_total_moment_of_inertia


! ***************************************************************************
!> \brief  Return link moment of inertia divided by m_He
!> \param helium ...
!> \param ia ...
!> \param ib ...
!> \retval moit ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_link_moment_of_inertia(helium,ia,ib) RESULT(moit)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: moit

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_link_moment_of_inertia', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ja1, ja2, jb1, jb2
    REAL(KIND=dp), DIMENSION(3)              :: com, r1, r12, r2

    com(:) = helium%center(:)

    IF (ib .EQ. helium%beads) THEN
      ja1 = ia
      ja2 = helium%permutation(ia)
      jb1 = ib
      jb2 = 1
    ELSE
      ja1 = ia
      ja2 = ia
      jb1 = ib
      jb2 = ib+1
    END IF
    r1(:) = helium%pos(:,ja1,jb1) - com(:)
    r2(:) = helium%pos(:,ja2,jb2) - com(:)
    ! comment out for non-PBC version -->
    r12(:) = r2(:) - r1(:)
    CALL helium_pbc(helium,r1)
    CALL helium_pbc(helium,r12)
    r2(:) = r1(:) + r12(:)
    ! comment out for non-PBC version <--
    moit(1) = r1(2)*r2(2)+r1(3)*r2(3)
    moit(2) = r1(3)*r2(3)+r1(1)*r2(1)
    moit(3) = r1(1)*r2(1)+r1(2)*r2(2)
    moit(:) = moit(:) / helium%beads

    RETURN
  END FUNCTION helium_link_moment_of_inertia


! ***************************************************************************
!> \brief  Return total moment of inertia (sum over all links)
!> \param helium ...
!> \retval moit ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_total_moment_of_inertia_linkwise(helium) RESULT(moit)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: moit

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_moment_of_inertia_linkwise', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib

    moit(:) = 0.0_dp
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        moit(:) = moit(:) + helium_link_moment_of_inertia(helium,ia,ib)
      END DO
    END DO

  END FUNCTION helium_total_moment_of_inertia_linkwise


! ***************************************************************************
!> \brief  Return moment of inertia of a cycle divided by m_He
!> \param helium ...
!> \param CYCLE ...
!> \param pos ...
!> \retval moit ...
!> \date   2014-04-24
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_cycle_moment_of_inertia(helium,CYCLE,pos) RESULT(moit)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, DIMENSION(:), POINTER           :: CYCLE
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos
    REAL(KIND=dp), DIMENSION(3)              :: moit

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_cycle_moment_of_inertia', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ia, ib, nsize
    REAL(KIND=dp), DIMENSION(3)              :: com, rcur, ri, rj

    nsize = SIZE(CYCLE)

    ! traverse the path
    moit(:) = 0.0_dp
    com(:) = helium_com(helium)
    DO ia = 1, nsize
      ! contributions from all bead pairs of the current atom
      DO ib = 1, helium%beads-1
        ri = pos(:,CYCLE(ia),ib) - com(:)
        rj = pos(:,CYCLE(ia),ib+1) - com(:)
        rcur(1) = ri(2)*rj(2)+ri(3)*rj(3)
        rcur(2) = ri(3)*rj(3)+ri(1)*rj(1)
        rcur(3) = ri(1)*rj(1)+ri(2)*rj(2)
        moit(:) = moit(:) + rcur(:)
      END DO
      ! contribution from the last bead of the current atom
      ! and the first bead of the next atom
      i1 = CYCLE(ia)
      IF (ia .EQ. nsize) THEN
        i2 = CYCLE(1)
      ELSE
        i2 = CYCLE(ia+1)
      END IF
      ! rotation invariant bead index
      ri = pos(:,i1,helium%beads) - com(:)
      rj = pos(:,i2,1) - com(:)
      rcur(1) = ri(2)*rj(2)+ri(3)*rj(3)
      rcur(2) = ri(3)*rj(3)+ri(1)*rj(1)
      rcur(3) = ri(1)*rj(1)+ri(2)*rj(2)
      moit(:) = moit(:) + rcur(:)
    END DO
    moit(:) = moit(:) / helium%beads

    RETURN
  END FUNCTION helium_cycle_moment_of_inertia


! ***************************************************************************
!> \brief  Return total moment of inertia (sum over all cycles)
!> \param helium ...
!> \retval moit ...
!> \date   2014-04-24
!> \author Lukasz Walewski
!> \note   mostly for sanity checks
! *****************************************************************************
  FUNCTION helium_total_moment_of_inertia_cyclewise(helium) RESULT(moit)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: moit

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'helium_total_moment_of_inertia_cyclewise', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic
    REAL(KIND=dp), DIMENSION(3)              :: pa
    TYPE(int_arr_ptr), DIMENSION(:), POINTER :: cycles

! decompose the current permutation state into permutation cycles

    NULLIFY(cycles)
    cycles => helium_calc_cycles(helium%permutation)

    moit(:) = 0.0_dp
    DO ic = 1, SIZE(cycles)
      pa = helium_cycle_moment_of_inertia(helium,cycles(ic)%iap,helium%pos)
      moit(:) = moit(:) + pa(:)
    END DO

    DEALLOCATE(cycles)

    RETURN
  END FUNCTION helium_total_moment_of_inertia_cyclewise


! ***************************************************************************
!> \brief  Set coordinate system, e.g. for RHO calculations
!> \param helium ...
!> \param pint_env ...
!> \date   2014-04-25
!> \author Lukasz Walewski
!> \note   Sets the origin (center of the coordinate system) wrt which
!>         spatial distribution functions are calculated.
!> \note   It can be extended later to set the axes of the coordinate system
!>         as well, e.g. for dynamic analysis with moving solute.
! *****************************************************************************
  SUBROUTINE helium_update_coord_system(helium,pint_env)

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_update_coord_system', &
      routineP = moduleN//':'//routineN

    IF (helium%solute_present) THEN
      helium%center(:) = pint_com_pos(pint_env)
    ELSE
      IF (helium%periodic) THEN
        helium%center(:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
      ELSE
        helium%center(:) = helium_com(helium)
      END IF
    END IF

    RETURN
  END SUBROUTINE helium_update_coord_system


! ***************************************************************************
!> \brief  Set coordinate system for RDF calculations
!> \param helium ...
!> \param pint_env ...
!> \date   2014-04-25
!> \author Lukasz Walewski
!> \note   Sets the centers wrt which radial distribution functions are
!>         calculated.
! *****************************************************************************
  SUBROUTINE helium_set_rdf_coord_system(helium,pint_env)

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_set_rdf_coord_system', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (helium%solute_present) THEN
      CALL pint_calc_centroid(pint_env)
      i = 3*helium%solute_atoms
      helium%rdf_centers(1:i) = pint_env%centroid(:)
    ELSE
      helium%rdf_centers(:) = helium%center(:)
    END IF

    RETURN
  END SUBROUTINE helium_set_rdf_coord_system


! ***************************************************************************
!> \brief  Calculate center of mass
!> \param helium ...
!> \retval com ...
!> \date   2014-04-09
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_com(helium) RESULT(com)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: com

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_com', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib

    com(:) = 0.0_dp
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        com(:) = com(:) + helium%pos(:,ia,ib)
      END DO
    END DO
    com(:) = com(:) / helium%atoms / helium%beads

  END FUNCTION helium_com


! ***************************************************************************
!> \brief  Return link vector, i.e. displacement vector of two consecutive
!>         beads along the path starting at bead ib of atom ia
!> \param helium ...
!> \param ia ...
!> \param ib ...
!> \retval lvec ...
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_link_vector(helium,ia,ib) RESULT(lvec)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: lvec

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_link_vector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia1, ia2, ib1, ib2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r1, r2

    IF (ib .EQ. helium%beads) THEN
      ia1 = ia
      ia2 = helium%permutation(ia)
      ib1 = ib
      ib2 = 1
    ELSE
      ia1 = ia
      ia2 = ia
      ib1 = ib
      ib2 = ib+1
    END IF
    r1 => helium%pos(:,ia1,ib1)
    r2 => helium%pos(:,ia2,ib2)
    lvec(:) = r2(:) - r1(:)
    CALL helium_pbc(helium,lvec)

  END FUNCTION helium_link_vector


! *****************************************************************************
!> \brief ...
!> \param helium ...
!> \param ia ...
!> \param ib ...
!> \retval rperp ...
! *****************************************************************************
  FUNCTION helium_rperpendicular(helium,ia,ib) RESULT(rperp)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER                                  :: ia, ib
    REAL(KIND=dp), DIMENSION(3)              :: rperp

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rperpendicular', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(:), POINTER     :: ri

    ri => helium%pos(:,ia,ib)
    rperp(1) = SQRT(ri(2)*ri(2)+ri(3)*ri(3))
    rperp(2) = SQRT(ri(3)*ri(3)+ri(1)*ri(1))
    rperp(3) = SQRT(ri(1)*ri(1)+ri(2)*ri(2))

    RETURN
  END FUNCTION helium_rperpendicular


! ***************************************************************************
!> \brief  Convert a winding number vector from real number representation
!>         (in internal units) to integer number representation (in cell
!>         vector units)
!> \param helium ...
!> \param wnum ...
!> \retval inum ...
!> \author Lukasz Walewski
! *****************************************************************************
  FUNCTION helium_wnumber_to_integer(helium,wnum) RESULT(inum)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3)              :: wnum
    INTEGER, DIMENSION(3)                    :: inum

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_wnumber_to_integer', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(3)              :: wcur
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: invcell

    invcell => helium%cell_m_inv
    CALL DGEMV('N',3,3,1.0_dp,invcell,SIZE(invcell,1),wnum,1,0.0_dp,wcur,1)
    inum(:) = NINT(wcur(:))

    RETURN
  END FUNCTION helium_wnumber_to_integer


! ***************************************************************************
!> \brief  Given the atom index and permutation state returns .TRUE. if the
!>         atom belongs to a winding path, .FASLE. otherwise.
!> \param helium ...
!> \param atmidx ...
!> \param pos ...
!> \param permutation ...
!> \retval is_winding ...
!> \date   2010-09-21
!> \author Lukasz Walewski
!> \note   The path winds around the periodic box if any component of its
!>         widing number vector differs from 0.
! *****************************************************************************
  FUNCTION helium_is_winding(helium, atmidx, pos, permutation) RESULT(is_winding)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: atmidx
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos
    INTEGER, DIMENSION(:), POINTER           :: permutation
    LOGICAL                                  :: is_winding

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_is_winding', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic
    INTEGER, DIMENSION(3)                    :: inum
    INTEGER, DIMENSION(:), POINTER           :: CYCLE
    REAL(KIND=dp), DIMENSION(3)              :: wnum

    is_winding = .FALSE.
    NULLIFY(CYCLE)
    CYCLE => helium_cycle_of(atmidx, permutation)
    wnum(:) = bohr * helium_cycle_winding_number(helium,CYCLE,pos)
    inum(:) = helium_wnumber_to_integer(helium,wnum)
    DO ic = 1, 3
      IF (ABS(inum(ic)) .GT. 0) THEN
        is_winding = .TRUE.
        EXIT
      END IF
    END DO
    DEALLOCATE(CYCLE)

    RETURN
  END FUNCTION helium_is_winding


END MODULE helium_common
