!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_nonbond_force [1.0] *
!!
!!   NAME
!!     fist_nonbond_force
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (11 May 2001) : cleaning up of support structures
!!     CJM & HAF (27 July 2001): fixed bug with handling of cutoff larger than
!!                               half the boxsize.
!!
!!   SOURCE
!******************************************************************************

MODULE fist_nonbond_force

  USE atomic_kind_types, ONLY : get_atomic_kind
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY : dbl, wp_size => dp_size
  USE mathconstants, ONLY : pi
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
                             linklist_atoms, linklist_exclusion
  USE pair_potential, ONLY : potential_s, potentialparm_type
  USE particle_types, ONLY : particle_type
  USE qs_neighbor_list_types,   ONLY: first_node,&
                                      get_neighbor_node,&
                                      neighbor_list_type,&
                                      neighbor_node_type,&
                                      next
  USE simulation_cell, ONLY : cell_type, pbc, get_cell_param, &
                              real_to_scaled, scaled_to_real
  USE termination, ONLY : stop_memory
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : include_list

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: force_nonbond, bonded_correct_gaussian

!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE force_nonbond ( ewald_param, part, pnode, box, potparm, &
                           pot_nonbond, f_nonbond, ptens_nonbond )

! Calculates the force and the potential of the minimum image, and
! the pressure tensor
!
  IMPLICIT NONE

! Arguments
  TYPE (ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE (potentialparm_type ), INTENT ( IN ), DIMENSION ( :, : ) :: potparm
  TYPE (particle_type ), INTENT ( IN ), DIMENSION ( : ) :: part
  TYPE (particle_node_type ), INTENT ( IN ), DIMENSION ( : ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: box
  REAL ( dbl ), INTENT ( OUT ) :: pot_nonbond
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_nonbond
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: ptens_nonbond

! Locals
  CHARACTER(LEN=*), PARAMETER :: routine_name = "force_nonbond"
  CHARACTER(LEN=*), PARAMETER :: module_name = "fist_nonbond_force"
  CHARACTER(LEN=*), PARAMETER :: routine =&
    "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"
  INTEGER :: ikind, jkind, nkinds, inode, natoms, nnodes
  INTEGER :: handle, istat, nneighbor, ineighbor, atom_a, atom_b, cell_b ( 3 )
  REAL ( dbl ), DIMENSION (3) :: sab_pbc, rab, rb, sab, ra
  REAL ( dbl ) :: energy, fscalar, rab2, flops
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ) :: rtest
  TYPE ( neighbor_list_type ), POINTER :: neighbor_list
  TYPE ( neighbor_node_type ), POINTER :: neighbor_node
  TYPE ( cell_type ), POINTER :: cell

!------------------------------------------------------------------------------

  CALL timeset ( routine_name, 'E', 'Mflops', handle )
  flops = 0._dbl
! create local pointer to box
  cell => box
! initializing the potential energy, pressure tensor and force
  pot_nonbond = 0._dbl
  f_nonbond = 0._dbl
  ptens_nonbond = 0._dbl

! local copy of cutoffs
  nkinds = SIZE ( potparm, 1 )
  ALLOCATE ( rtest ( nkinds, nkinds ), STAT = istat )
  IF ( istat /= 0 ) CALL stop_memory ( routine, "rtest", nkinds ** 2 )

  DO ikind = 1, nkinds
     DO jkind = 1, nkinds
        rtest ( ikind, jkind ) = potparm ( ikind, jkind ) % rcutsq
     END DO
  END DO

!
! starting the force loop
!
  nnodes = SIZE ( pnode )
  DO inode = 1, nnodes
     atom_a = pnode ( inode ) % p % iatom
     neighbor_list => pnode ( inode ) % nl
     nneighbor = pnode ( inode ) % nneighbor
     ra ( : ) = part ( atom_a ) % r ( : )
     CALL get_atomic_kind ( part ( atom_a ) % atomic_kind,  &
                            kind_number = ikind )
! now do neighbors
     neighbor_node => first_node ( neighbor_list )
     DO ineighbor = 1, nneighbor
        CALL get_neighbor_node ( neighbor_node = neighbor_node, &
                                 neighbor = atom_b, &
                                 cell = cell_b )
        rb ( : ) = part ( atom_b ) % r ( : )
        CALL get_atomic_kind ( part ( atom_b ) % atomic_kind, &
                               kind_number = jkind )
        rab ( : ) = rb ( : ) - ra ( : )
        sab ( : ) = real_to_scaled(rab(:),cell)
        sab_pbc ( : ) = sab ( : ) + REAL ( cell_b ( : ), dbl )
        rab ( : ) = scaled_to_real(sab_pbc(:),cell)
        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

        IF ( rab2 <= rtest ( ikind, jkind ) ) THEN
           CALL potential_s ( rab2, ikind, jkind, energy, fscalar )
!
! summing up the potential energy,the force and pressure tensor
!
           CALL sum_ener_forces ( atom_a, atom_b, pot_nonbond, energy, &
                fscalar, f_nonbond, rab, ptens_nonbond )
           flops = flops + 64.0_dbl
        END IF

       neighbor_node => next ( neighbor_node )

     END DO
  END DO

! computing long range corrections to the potential
!
!when we get the system for handling multiple potentials
!we will add the long range correction
!
!       pot_nonbond=pot_nonbond+lrc*(1./box % deth)
!

  DEALLOCATE ( rtest, STAT = istat )
  IF ( istat /= 0 ) CALL stop_memory ( routine, "rtest" )
!
  flops = flops * 1.E-6_dbl
  CALL timestop ( flops, handle )

END SUBROUTINE force_nonbond

!******************************************************************************

SUBROUTINE sum_ener_forces ( i, j, pot, e, fs, f, rij, pv )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: i, j
  REAL ( dbl ), INTENT ( INOUT ) :: pot
  REAL ( dbl ), INTENT ( IN ) :: e, fs
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: f, pv
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: rij

! Locals
  REAL ( dbl ), DIMENSION(3) :: fr

!------------------------------------------------------------------------------

  pot = pot + e
  f(1,i) = f(1,i) - fs*rij(1)
  f(2,i) = f(2,i) - fs*rij(2)
  f(3,i) = f(3,i) - fs*rij(3)
  f(1,j) = f(1,j) + fs*rij(1)
  f(2,j) = f(2,j) + fs*rij(2)
  f(3,j) = f(3,j) + fs*rij(3)
  fr(1) = rij(1)*fs
  fr(2) = rij(2)*fs
  fr(3) = rij(3)*fs
  pv(1,1) = pv(1,1) + rij(1)*fr(1)
  pv(1,2) = pv(1,2) + rij(1)*fr(2)
  pv(1,3) = pv(1,3) + rij(1)*fr(3)
  pv(2,1) = pv(2,1) + rij(2)*fr(1)
  pv(2,2) = pv(2,2) + rij(2)*fr(2)
  pv(2,3) = pv(2,3) + rij(2)*fr(3)
  pv(3,1) = pv(3,1) + rij(3)*fr(1)
  pv(3,2) = pv(3,2) + rij(3)*fr(2)
  pv(3,3) = pv(3,3) + rij(3)*fr(3)

END SUBROUTINE sum_ener_forces

!******************************************************************************

SUBROUTINE bonded_correct_gaussian ( ewald_param, molecule, v_bonded_corr, &
     pv_bc )

  IMPLICIT NONE

! Arguments
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( molecule_structure_type ), INTENT ( IN ) :: molecule ( : )
  REAL ( dbl ), INTENT ( OUT ) :: v_bonded_corr
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: pv_bc

! Locals
  TYPE ( linklist_exclusion ), POINTER :: llex
  TYPE ( linklist_atoms ), POINTER :: llat
! parameters for the numerical erf in Abramowitz and Stegun Eq. 7.1.26 pg 299
  REAL ( dbl ), PARAMETER :: ac1 = 0.254829592_dbl, &
       ac2 = -0.284496736_dbl, ac3 = 1.421413741_dbl, &
       ac4 = -1.453152027_dbl, ac5 = 1.061405429_dbl, pc = 0.3275911_dbl
  REAL ( dbl ) :: arg, e_arg_arg, tc, errf, alpha, fscalar, flops
  REAL ( dbl ) :: const, ifourpieps0, dij, idij, rijsq, qi, qj
  REAL ( dbl ), DIMENSION (3) :: rij
  INTEGER :: i, j, nmol, imol, iat, handle

!------------------------------------------------------------------------------

  CALL timeset ( 'BOND_CORRECT_GAUSSIAN', 'E', 'Mflops', handle )
  flops = 0.0_dbl

! defining the constants
  pv_bc = 0.0_dbl
  v_bonded_corr = 0.0_dbl

!
!
  alpha = ewald_param % alpha
  const = 2.0_dbl * alpha / SQRT( pi )
  ifourpieps0 = 1.0_dbl / ( 4.0_dbl * pi * ewald_param % eps0 )

  nmol = SIZE ( molecule )

  DO imol = 1, nmol
     llat => molecule ( imol ) % ll_atoms

     DO iat = 1, molecule ( imol ) % natoms_mol
        qi = llat % part % p % prop % charge
        llex => llat % part % ex

        DO j = 1, llat % part % nexcl
!defining the charge
           qj = llex % p % prop % charge
           rij(1) = llat % part % p % r(1) - llex % p % r(1)
           rij(2) = llat % part % p % r(2) - llex % p % r(2)
           rij(3) = llat % part % p % r(3) - llex % p % r(3)
           rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
           dij = SQRT ( rijsq )
           idij = 1.0_dbl / dij
           arg = alpha * dij
           e_arg_arg = EXP ( -arg ** 2 )
           tc = 1.0_dbl / ( 1.0_dbl + pc * arg )

! defining errf=1-erfc
           errf = 1.0_dbl &
                - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg

! getting the potential
           v_bonded_corr = v_bonded_corr - qi*qj*idij*errf

! subtracting the force from the total force
           fscalar = ifourpieps0 * qi * qj * idij ** 2 &
                * ( idij * errf - const * EXP ( -alpha ** 2 * rijsq ) )

           llat % part % p % f(1) = llat % part % p % f(1) - fscalar*rij(1)
           llat % part % p % f(2) = llat % part % p % f(2) - fscalar*rij(2)
           llat % part % p % f(3) = llat % part % p % f(3) - fscalar*rij(3)

! computing the pressure tensor
           pv_bc(1,1) = pv_bc(1,1) - fscalar * rij(1) * llat % part % p % r(1)
           pv_bc(1,2) = pv_bc(1,2) - fscalar * rij(1) * llat % part % p % r(2)
           pv_bc(1,3) = pv_bc(1,3) - fscalar * rij(1) * llat % part % p % r(3)
           pv_bc(2,1) = pv_bc(2,1) - fscalar * rij(2) * llat % part % p % r(1)
           pv_bc(2,2) = pv_bc(2,2) - fscalar * rij(2) * llat % part % p % r(2)
           pv_bc(2,3) = pv_bc(2,3) - fscalar * rij(2) * llat % part % p % r(3)
           pv_bc(3,1) = pv_bc(3,1) - fscalar * rij(3) * llat % part % p % r(1)
           pv_bc(3,2) = pv_bc(3,2) - fscalar * rij(3) * llat % part % p % r(2)
           pv_bc(3,3) = pv_bc(3,3) - fscalar * rij(3) * llat % part % p % r(3)

           flops = flops + 62.0_dbl
           llex => llex % next
        END DO

        llat => llat % next
     END DO
  END DO

! the factor of 1/2 comes from double counting in the exclusion list
  v_bonded_corr = v_bonded_corr * ifourpieps0 * 0.5_dbl

  flops = flops * 1.0E-6_dbl
  CALL timestop ( flops, handle )

END SUBROUTINE bonded_correct_gaussian

!******************************************************************************

END MODULE fist_nonbond_force
