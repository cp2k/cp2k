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

  USE atomic_kind_types, ONLY : atomic_kind_type, get_atomic_kind
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY : dbl, dp_size
  USE mathconstants, ONLY : pi
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
                             linklist_atoms, linklist_exclusion
  USE pair_potential_types, ONLY : pair_potential_type
  USE pair_potential, ONLY : potential_s
  USE particle_types, ONLY : particle_type
  USE qs_neighbor_list_types,   ONLY: first_list,&
                                      first_node,&
                                      get_neighbor_list,&
                                      get_neighbor_list_set,&
                                      get_neighbor_node,&
                                      neighbor_list_set_type,&
                                      neighbor_list_type,&
                                      neighbor_node_type,&
                                      next
  USE simulation_cell, ONLY : cell_type, pbc, get_cell_param, &
                              real_to_scaled, scaled_to_real
  USE termination, ONLY : stop_memory
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : include_list

  USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: force_nonbond, bonded_correct_gaussian

!!***
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE force_nonbond ( ewald_param, part, pnode, box, potparm, &
                           pot_nonbond, f_nonbond, ptens_nonbond, &
                           nonbonded, r_last_update )

! Calculates the force and the potential of the minimum image, and
! the pressure tensor
!
  IMPLICIT NONE

! Arguments
  TYPE (ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE (pair_potential_type ), INTENT ( IN ), DIMENSION ( :, : ) :: potparm
  TYPE (particle_type ), INTENT ( IN ), DIMENSION ( : ) :: part
  TYPE (particle_node_type ), INTENT ( IN ), DIMENSION ( : ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: box
  REAL ( dbl ), INTENT ( OUT ) :: pot_nonbond
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_nonbond
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: ptens_nonbond

  TYPE(neighbor_list_set_p_type), DIMENSION(:), INTENT(IN) :: nonbonded
  REAL(dbl), DIMENSION(:,:), INTENT(IN)                    :: r_last_update

! Locals
  CHARACTER(LEN=*), PARAMETER :: routine_name = "force_nonbond"
  CHARACTER(LEN=*), PARAMETER :: module_name = "fist_nonbond_force"
  CHARACTER(LEN=*), PARAMETER :: routine =&
    "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"
  INTEGER :: ikind, jkind, nkinds, inode, natoms, nnodes, ab
  INTEGER :: handle, istat, nneighbor, ineighbor, atom_a, atom_b, cell_b ( 3 )
  REAL ( dbl ), DIMENSION (3) :: sab_pbc, rab, rb, sab, ra
  REAL ( dbl ) :: energy, fscalar, rab2, flops, rab2_max
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ) :: rtest
  TYPE ( neighbor_list_type ), POINTER :: neighbor_list
  TYPE ( neighbor_node_type ), POINTER :: neighbor_node
  TYPE ( cell_type ), POINTER :: cell
!MK
  TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

  INTEGER :: ilist,nlist,nnode

  REAL(dbl), DIMENSION(3) :: dra,drb,fr,rab_last_update
  REAL(dbl), DIMENSION(3) :: sb

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
!MK  ALLOCATE ( rtest ( nkinds, nkinds ), STAT = istat )
!MK  IF ( istat /= 0 ) CALL stop_memory ( routine, "rtest", nkinds ** 2 )
!MK
!MK  DO ikind = 1, nkinds
!MK     DO jkind = 1, nkinds
!MK        rtest ( ikind, jkind ) = potparm ( ikind, jkind ) % rcutsq
!MK     END DO
!MK  END DO

!
! starting the force loop
!
!MK
    DO ikind=1,nkinds
      DO jkind=ikind,nkinds

        ab = ikind + jkind*(jkind - 1)/2

        neighbor_list_set => nonbonded(ab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        rab2_max = potparm(ikind,jkind)%rcutsq

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=atom_a,&
                                 nnode=nnode)

          dra(:) = part(atom_a)%r(:) - scaled_to_real ( r_last_update(:,atom_a), cell )

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=atom_b,&
                                   cell=cell_b,&
                                   r=rab_last_update)

            drb(:) = part(atom_b)%r(:) - scaled_to_real ( r_last_update(:,atom_b), cell )
         !   drb(:) = pbc ( part(atom_b)%r(:), cell ) 
         !   sb ( : ) = real_to_scaled ( drb ( : ), cell )
         !   sb ( : ) = sb ( : ) + REAL ( cell_b, dbl )
         !   drb ( : ) = scaled_to_real ( sb, cell )

            rab ( : ) =  scaled_to_real ( rab_last_update ( : ), cell ) - dra(:) + drb(:)
         !   rab(:) = rab(:) - dra(:) + drb(:)
         !  rab(:) =  drb ( : ) - dra ( : ) 
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

        !    write ( *, * )  rab2, rab2_max, atom_a, atom_b
            IF (rab2 <= rab2_max) THEN
              CALL potential_s(rab2,ikind,jkind,energy,fscalar)
        !      write ( *, * )  rab2, energy
        !      write ( *, * ) part ( atom_a) % r
!MK           CALL sum_ener_forces(atom_a,atom_b,pot_nonbond,energy,&
!MK                                fscalar,f_nonbond,rab,ptens_nonbond)
              pot_nonbond = pot_nonbond + energy
              fr(1) = fscalar*rab(1)
              fr(2) = fscalar*rab(2)
              fr(3) = fscalar*rab(3)
              f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
              f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
              f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)
              f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
              f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
              f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
              ptens_nonbond(1,1) = ptens_nonbond(1,1) + rab(1)*fr(1)
              ptens_nonbond(1,2) = ptens_nonbond(1,2) + rab(1)*fr(2)
              ptens_nonbond(1,3) = ptens_nonbond(1,3) + rab(1)*fr(3)
              ptens_nonbond(2,1) = ptens_nonbond(2,1) + rab(2)*fr(1)
              ptens_nonbond(2,2) = ptens_nonbond(2,2) + rab(2)*fr(2)
              ptens_nonbond(2,3) = ptens_nonbond(2,3) + rab(2)*fr(3)
              ptens_nonbond(3,1) = ptens_nonbond(3,1) + rab(3)*fr(1)
              ptens_nonbond(3,2) = ptens_nonbond(3,2) + rab(3)*fr(2)
              ptens_nonbond(3,3) = ptens_nonbond(3,3) + rab(3)*fr(3)
            END IF

            neighbor_node => next(neighbor_node)

          END DO

          neighbor_list => next(neighbor_list)

        END DO

      END DO

    END DO
!MK  nnodes = SIZE ( pnode )
!MK  DO inode = 1, nnodes
!MK     atom_a = pnode ( inode ) % p % iatom
!MK     neighbor_list => pnode ( inode ) % nl
!MK     nneighbor = pnode ( inode ) % nneighbor
!MK     ra ( : ) = part ( atom_a ) % r ( : )
!MK     CALL get_atomic_kind ( part ( atom_a ) % atomic_kind,  &
!MK                            kind_number = ikind )
!MK! now do neighbors
!MK     neighbor_node => first_node ( neighbor_list )
!MK     DO ineighbor = 1, nneighbor
!MK        CALL get_neighbor_node ( neighbor_node = neighbor_node, &
!MK                                 neighbor = atom_b, &
!MK                                 cell = cell_b )
!MK        rb ( : ) = part ( atom_b ) % r ( : )
!MK        CALL get_atomic_kind ( part ( atom_b ) % atomic_kind, &
!MK                               kind_number = jkind )
!MK        rab ( : ) = rb ( : ) - ra ( : )
!MK        sab ( : ) = real_to_scaled(rab(:),cell)
!MK        sab_pbc ( : ) = sab ( : ) + REAL ( cell_b ( : ), dbl )
!MK        rab ( : ) = scaled_to_real(sab_pbc(:),cell)
!MK        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
!MK
!MK        IF ( rab2 <= rtest ( ikind, jkind ) ) THEN
!MK           CALL potential_s ( rab2, ikind, jkind, energy, fscalar )
!MK!
!MK! summing up the potential energy,the force and pressure tensor
!MK!
!MK           CALL sum_ener_forces ( atom_a, atom_b, pot_nonbond, energy, &
!MK                fscalar, f_nonbond, rab, ptens_nonbond )
!MK           flops = flops + 64.0_dbl
!MK        END IF
!MK
!MK       neighbor_node => next ( neighbor_node )
!MK
!MK     END DO
!MK  END DO

! computing long range corrections to the potential
!
!when we get the system for handling multiple potentials
!we will add the long range correction
!
!       pot_nonbond=pot_nonbond+lrc*(1./box % deth)
!

!MK  DEALLOCATE ( rtest, STAT = istat )
!MK  IF ( istat /= 0 ) CALL stop_memory ( routine, "rtest" )
!MK
!MK  flops = flops * 1.E-6_dbl
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
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind

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
        atomic_kind => llat % part % p % atomic_kind
        CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=qi)
        llex => llat % part % ex

        DO j = 1, llat % part % nexcl
!defining the charge
           atomic_kind => llex % p % atomic_kind
           CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=qj)
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
