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

  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
       linklist_atoms, linklist_neighbor, linklist_exclusion, linklist_images
  USE pair_potential, ONLY : potential_s, potentialparm_type
  USE simulation_cell, ONLY : cell_type, pbc, get_cell_param
  USE termination, ONLY : stop_memory
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : include_list

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: force_nonbond, bonded_correct_gaussian, find_image

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE force_nonbond ( ewald_param, pnode, box, potparm, pot_nonbond, &
     f_nonbond, ptens_nonbond )

! Calculates the force and the potential of the minimum image, and
! the pressure tensor
!
  IMPLICIT NONE

! Arguments
  TYPE (ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE (potentialparm_type ), INTENT ( IN ), DIMENSION ( :, : ) :: potparm
  TYPE (particle_node_type ), INTENT ( IN ), DIMENSION ( : ) :: pnode
  TYPE (cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( OUT ) :: pot_nonbond
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_nonbond
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: ptens_nonbond

! Locals
  INTEGER :: i, j, ii, jj, iatomtype, jatomtype, imol, nmol, iat
  INTEGER :: id,itype,jtype,isos,natom_types
  INTEGER :: jatom, handle, ios, im
  REAL ( dbl ), DIMENSION (3) :: ri, rij, s, perd, vec
  REAL ( dbl ), DIMENSION (3,3) :: hmat, h_inv
  REAL ( dbl ) :: energy, fscalar, rijsq, flops
  TYPE (linklist_neighbor), POINTER :: current_neighbor
  TYPE (linklist_images), POINTER :: current_image
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ) :: rtest

!------------------------------------------------------------------------------

  CALL timeset ( 'NONBOND', 'E', 'Mflops', handle )
  flops = 0._dbl

! initializing the potential energy, pressure tensor and force
  pot_nonbond = 0._dbl
  f_nonbond = 0._dbl
  ptens_nonbond = 0._dbl

! local copies of the box parameters
  hmat = box % hmat
  h_inv = box % h_inv
  perd = box % perd

! local copy of cutoffs
  natom_types = SIZE ( potparm, 1 )
  ALLOCATE ( rtest ( natom_types, natom_types ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'fist_nonbond_force', &
          'rtest', natom_types ** 2 )

  DO i = 1, natom_types
     DO j = 1, natom_types
        rtest ( i, j ) = potparm ( i, j ) % rcutsq
     END DO
  END DO

!
! starting the force loop
!
  DO iat = 1, SIZE ( pnode )
     iatomtype = pnode ( iat ) % p % prop % ptype
     i = pnode ( iat ) % p % iatom
     ri = pnode ( iat ) % p % r
! now do neighbors
     current_neighbor => pnode ( iat ) % nl
     DO j = 1, pnode ( iat ) % nneighbor
        jatom = current_neighbor % index
        jatomtype = current_neighbor % p % prop % ptype
        rij = current_neighbor % p % r - ri
!
! apply periodic boundary conditions; the following code is the inline version
! of the general function pbc
!
        s ( 1 ) = h_inv ( 1, 1 ) * rij ( 1 ) + h_inv ( 1, 2 ) * rij ( 2 ) &
                + h_inv ( 1, 3 ) * rij ( 3 )
        s ( 2 ) = h_inv ( 2, 1 ) * rij ( 1 ) + h_inv ( 2, 2 ) * rij ( 2 ) &
                + h_inv ( 2, 3 ) * rij ( 3 )
        s ( 3 ) = h_inv ( 3, 1 ) * rij ( 1 ) + h_inv ( 3, 2 ) * rij ( 2 ) & 
                + h_inv ( 3, 3 ) * rij ( 3 )
        IF ( current_neighbor % nimages == 0 ) THEN
!
! if there are no images in the image-list, then check the minimum image
!
          IF ( ABS ( s ( 1 ) ) > 0.5_dbl ) THEN
             s(1) = s(1) - perd(1) * INT(s(1)+SIGN(0.5_dbl,s(1)))
          END IF
          IF ( ABS ( s ( 2 ) ) > 0.5_dbl ) THEN
             s(2) = s(2) - perd(2) * INT(s(2)+SIGN(0.5_dbl,s(2)))
          END IF
          IF ( ABS ( s ( 3 ) ) > 0.5_dbl ) THEN
             s(3) = s(3) - perd(3) * INT(s(3)+SIGN(0.5_dbl,s(3)))
          END IF
          rij(1) = hmat(1,1)*s(1) + hmat(1,2)*s(2) + hmat(1,3)*s(3)
          rij(2) = hmat(2,1)*s(1) + hmat(2,2)*s(2) + hmat(2,3)*s(3)
          rij(3) = hmat(3,1)*s(1) + hmat(3,2)*s(2) + hmat(3,3)*s(3)
          rijsq = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
          flops = flops + 41.0_dbl

          IF ( i /= jatom .AND. rijsq <= rtest ( iatomtype, jatomtype ) ) THEN
             CALL potential_s ( rijsq, iatomtype, jatomtype, &
                                energy, fscalar )
!
! summing up the potential energy,the force and pressure tensor
!
             CALL sum_ener_forces ( i, jatom, pot_nonbond, energy, &
                  fscalar, f_nonbond, rij, ptens_nonbond )
             flops = flops + 64.0_dbl
          END IF
        ELSE
!
! there are images in the image list
!
          flops = flops + 19.0_dbl
!
! now sum over lattice translations of neighbors
!
          current_image => current_neighbor % image
          DO im = 1, current_neighbor % nimages

             CALL find_image ( s, perd, current_image % vec, hmat, rijsq, rij )
             flops = flops + 41.0_dbl
             IF ( rijsq <= rtest ( iatomtype, jatomtype ) ) THEN
                CALL potential_s ( rijsq, iatomtype, jatomtype, &
                                   energy, fscalar )
!
! summing up the potential energy,the force and pressure tensor
!

                CALL sum_ener_forces ( i, jatom, pot_nonbond, energy, &
                     fscalar, f_nonbond, rij, ptens_nonbond )
                flops = flops + 64.0_dbl
             END IF

             current_image => current_image % next
          END DO

        ENDIF

        current_neighbor => current_neighbor % next
     END DO

  END DO

! computing long range corrections to the potential
!
!when we get the system for handling multiple potentials
!we will add the long range correction
!
!       pot_nonbond=pot_nonbond+lrc*(1./box % deth)
!
  DEALLOCATE ( rtest, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'fist_nonbond_force', 'rtest' )
!
  flops = flops * 1.E-6_dbl
  CALL timestop ( flops, handle )

END SUBROUTINE force_nonbond

!******************************************************************************

SUBROUTINE find_image ( s, perd, vec, hmat, rijsq, rij )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: s,perd
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: vec
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: hmat
  REAL ( dbl ), INTENT ( OUT ) :: rijsq
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ) :: rij

! Locals
  REAL ( dbl ) :: strans(3)

!------------------------------------------------------------------------------

  strans ( 1 ) = s ( 1 ) + perd ( 1 ) * vec ( 1 )
  strans ( 2 ) = s ( 2 ) + perd ( 2 ) * vec ( 2 )
  strans ( 3 ) = s ( 3 ) + perd ( 3 ) * vec ( 3 )
  rij(1) = hmat(1,1)*strans(1) + hmat(1,2)*strans(2) + hmat(1,3)*strans(3)
  rij(2) = hmat(2,1)*strans(1) + hmat(2,2)*strans(2) + hmat(2,3)*strans(3)
  rij(3) = hmat(3,1)*strans(1) + hmat(3,2)*strans(2) + hmat(3,3)*strans(3)

  rijsq = rij ( 1 ) * rij ( 1 ) + rij ( 2 ) * rij ( 2 ) + rij ( 3 ) * rij ( 3 )

END SUBROUTINE find_image

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
