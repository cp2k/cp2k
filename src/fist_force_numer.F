!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_force_numer [1.0] *
!!
!!   NAME
!!     fist_force_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Torsions added (DG) 05-Dec-2000
!!     Variable names changed (DG) 05-Dec-2000
!!
!!   SOURCE
!******************************************************************************

MODULE fist_force_numer
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE dg_types,                        ONLY: dg_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE kinds,                           ONLY: dbl,&
                                             dp_size
  USE mol_force,                       ONLY: force_bends,&
                                             force_bonds,&
                                             force_torsions
  USE molecule_types,                  ONLY: linklist_bends,&
                                             linklist_bonds,&
                                             linklist_exclusion,&
                                             linklist_torsions,&
                                             molecule_structure_type,&
                                             particle_node_type
  USE pair_potential,                  ONLY: potential_f
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: particle_type
  USE pme,                             ONLY: pme_evaluate
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_construct,&
                                             pw_grid_destruct,&
                                             pw_grid_setup
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             get_hinv,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory

  TYPE pbc_coord_type
    REAL(dbl), DIMENSION(:,:), POINTER :: r,s
  END TYPE pbc_coord_type

  PRIVATE
  PUBLIC :: force_bond_numer, force_bend_numer, force_torsion_numer, &
            force_nonbond_numer, ewald_force_numer
  PUBLIC :: pv_bond_numer, pv_bend_numer, pv_torsion_numer, ptens_numer
  PUBLIC :: pvg_numer, de_g_numer, potential_g_numer, ewald_energy_numer
  PUBLIC :: pme_energy_numer, pme_force_numer

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fist_force_numer/force_bond_numer [1.0] *
!!
!!   NAME
!!     force_bond_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_bond_numer ( molecule, f_bond_numer, pot_bond, &
     numerical_shift )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: f_bond_numer
    REAL(dbl), INTENT(OUT)                   :: pot_bond
    REAL(dbl), INTENT(IN)                    :: numerical_shift

    INTEGER                                  :: atom1, atom2, ibond, idim, &
                                                imol
    REAL(dbl)                                :: e0, em, ep, fscalar
    REAL(dbl), DIMENSION(3)                  :: delta, rij, rij_minus_delta, &
                                                rij_plus_delta
    TYPE(linklist_bonds), POINTER            :: llbond

!------------------------------------------------------------------------------
! initializing the forces

  f_bond_numer = 0.0_dbl
  pot_bond = 0.0_dbl
  MOL: DO imol = 1, SIZE ( molecule )
     llbond => molecule(imol) % ll_bonds

     BOND: DO ibond = 1, molecule(imol) % nbonds_mol
        atom1 = llbond % index(1)
        atom2 = llbond % index(2)
        rij = llbond % p1 % r - llbond % p2 % r
        delta = 0.0_dbl

! get the energy
        CALL force_bonds ( rij, llbond % bond_param % r0, &
             llbond % bond_param % k, e0, fscalar )
        pot_bond = pot_bond + e0

! get the forces
        DO idim = 1, 3
           delta(idim) = numerical_shift
           rij_minus_delta = rij - delta
           CALL force_bonds ( rij_minus_delta, llbond % bond_param % r0, &
                llbond % bond_param % k, em, fscalar )
           rij_plus_delta = rij + delta
           CALL force_bonds ( rij_plus_delta, llbond % bond_param % r0, &
                llbond % bond_param % k, ep, fscalar )
           f_bond_numer(idim,atom1) = f_bond_numer(idim,atom1) + (em-ep)
           f_bond_numer(idim,atom2) = f_bond_numer(idim,atom2) - (em-ep)
           delta = 0.0_dbl
        END DO

        llbond => llbond % next
     END DO BOND
  END DO MOL

  f_bond_numer = f_bond_numer / ( 2.0_dbl * numerical_shift )

END SUBROUTINE force_bond_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/force_bend_numer [1.0] *
!!
!!   NAME
!!     force_bend_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_bend_numer ( molecule, f_bend_numer, pot_bend, &
     numerical_shift )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: f_bend_numer
    REAL(dbl), INTENT(OUT)                   :: pot_bend
    REAL(dbl), INTENT(IN)                    :: numerical_shift

    INTEGER                                  :: atom1, atom2, atom3, ibend, &
                                                idim, imol
    REAL(dbl)                                :: d12, d12m, d12p, d32, d32m, &
                                                d32p, dist, e0, em, ep, &
                                                fscalar, id12, id12m, id12p, &
                                                id32, id32m, id32p, theta
    REAL(dbl), DIMENSION(3)                  :: b12, b12_minus_delta, &
                                                b12_plus_delta, b32, &
                                                b32_minus_delta, &
                                                b32_plus_delta, delta, g1, &
                                                g2, g3
    TYPE(linklist_bends), POINTER            :: llbend

!------------------------------------------------------------------------------
! initializing the forces

  f_bend_numer = 0.0_dbl
  pot_bend = 0.0_dbl
  MOL: DO imol = 1, SIZE ( molecule )
     llbend => molecule(imol) % ll_bends

     BEND: DO ibend = 1, molecule(imol) % nbends_mol
        atom1 = llbend % index(1)
        atom2 = llbend % index(2)
        atom3 = llbend % index(3)
        b12 = llbend % p1 % r - llbend % p2 % r
        b32 = llbend % p3 % r - llbend % p2 % r
        d12 = sqrt(dot_product(b12,b12))
        id12 = 1.0_dbl/d12
        d32 = sqrt(dot_product(b32,b32))
        id32 = 1.0_dbl/d32
        dist = dot_product(b12,b32)
        theta = acos(dist*id12*id32)

! get the energy
        CALL force_bends ( b12, b32, d12, d32, id12, id32, dist, theta, &
             llbend % bend_param % theta0, llbend % bend_param % k, &
             g1, g2, g3, e0, fscalar )
        pot_bend = pot_bend + e0

! get the force on atom1 and atom2
        delta = 0.0_dbl
        DO idim = 1, 3
           delta(idim) = numerical_shift
           b12_minus_delta = b12 - delta
           d12m = sqrt(dot_product(b12_minus_delta,b12_minus_delta))
           id12m = 1.0_dbl / d12m
           dist = dot_product(b12_minus_delta,b32)
           theta = acos(dist*id12m*id32)
           CALL force_bends ( b12_minus_delta, b32, d12m, d32, id12m, id32, &
                dist, theta, llbend % bend_param % theta0, &
                llbend % bend_param % k, g1, g2, g3, em, fscalar )
           b12_plus_delta = b12 + delta
           d12p = sqrt(dot_product(b12_plus_delta,b12_plus_delta))
           id12p = 1.0_dbl/d12p
           dist = dot_product(b12_plus_delta,b32)
           theta = acos(dist*id12p*id32)
           CALL force_bends ( b12_plus_delta, b32, d12p, d32, id12p, id32, &
                dist, theta, llbend % bend_param % theta0, &
                llbend % bend_param % k, g1, g2, g3, ep, fscalar )
           f_bend_numer(idim,atom1) = f_bend_numer(idim,atom1) + (em-ep)
           f_bend_numer(idim,atom2) = f_bend_numer(idim,atom2) - (em-ep)
           delta = 0.0_dbl
        END DO

! get the force on atom3 and atom2
        delta = 0.0_dbl
        DO idim = 1, 3
           delta(idim) = numerical_shift
           b32_minus_delta = b32 - delta
           d32m = sqrt(dot_product(b32_minus_delta,b32_minus_delta))
           id32m = 1.0_dbl/d32m
           dist = dot_product(b32_minus_delta,b12)
           theta = acos(dist*id32m*id12)
           CALL force_bends(b12,b32_minus_delta,d12,d32m,id12,id32m,dist, &
                theta, llbend % bend_param % theta0, &
                llbend % bend_param % k, g1, g2, g3, em, fscalar )
           b32_plus_delta = b32 + delta
           d32p = sqrt(dot_product(b32_plus_delta,b32_plus_delta))
           id32p = 1.0_dbl/d32p
           dist = dot_product(b32_plus_delta,b12)
           theta = acos(dist*id32p*id12)
           CALL force_bends ( b12, b32_plus_delta, d12, d32p, id12, id32p, &
                dist, theta, llbend % bend_param % theta0, &
                llbend % bend_param % k, g1, g2, g3, ep, fscalar )
           f_bend_numer(idim,atom3) = f_bend_numer(idim,atom3) + (em-ep)
           f_bend_numer(idim,atom2) = f_bend_numer(idim,atom2) - (em-ep)
           delta = 0.0_dbl
        END DO

        llbend => llbend % next
     END DO BEND
  END DO MOL

  f_bend_numer = f_bend_numer / ( 2.0_dbl * numerical_shift )

END SUBROUTINE force_bend_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/force_torsion_numer [1.0] *
!!
!!   NAME
!!     force_torsion_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_torsion_numer ( molecule, f_torsion_numer, pot_torsion, &
     numerical_shift )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: f_torsion_numer
    REAL(dbl), INTENT(OUT)                   :: pot_torsion
    REAL(dbl), INTENT(IN)                    :: numerical_shift

    INTEGER                                  :: atom1, atom2, atom3, atom4, &
                                                idim, imol, itorsion
    REAL(dbl) :: dist1, dist1m, dist1p, dist2, dist2m, dist2p, e0, em, ep, &
      fscalar, is32, is32m, is32p, iw1, iw1m, iw1p, iw2, iw2m, iw2p, s32, &
      s32m, s32p, w1, w1m, w1p, w2, w2m, w2p
    REAL(dbl), DIMENSION(3) :: delta, gt1, gt2, gt3, gt4, t12, &
      t12_minus_delta, t12_plus_delta, t32, t32_minus_delta, t32_plus_delta, &
      t43, t43_minus_delta, t43_plus_delta, u1, u1m, u1p, u2, u2m, u2p
    TYPE(linklist_torsions), POINTER         :: lltorsion

!------------------------------------------------------------------------------
! initializing the forces

  f_torsion_numer = 0.0_dbl
  pot_torsion = 0.0_dbl
  MOL: DO imol = 1, SIZE ( molecule )
     lltorsion => molecule(imol) % ll_torsions

     torsion: DO itorsion = 1, molecule(imol) % ntorsions_mol
        atom1 = lltorsion % index(1)
        atom2 = lltorsion % index(2)
        atom3 = lltorsion % index(3)
        atom4 = lltorsion % index(4)
        t12 = lltorsion%p1%r - lltorsion%p2%r
        t32 = lltorsion%p3%r - lltorsion%p2%r
        t43 = lltorsion%p4%r - lltorsion%p3%r
        s32 = sqrt(dot_product(t32,t32))
        is32 = 1._dbl/s32
        dist1 = dot_product(t12,t32)
        dist2 = dot_product(t43,t32)
        u1 = t12 - dist1*t32*is32**2
        u2 = t43 - dist2*t32*is32**2
        w1 = sqrt(dot_product(u1,u1))
        w2 = sqrt(dot_product(u2,u2))
        iw1 = 1._dbl/w1
        iw2 = 1._dbl/w2

! get the energy
        CALL force_torsions(is32,dist1,dist2,u1,u2,iw1,iw2, &
             lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,e0,fscalar)
        pot_torsion = pot_torsion + e0

! get the force on atom1 and atom2
        delta = 0.0_dbl
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t12_minus_delta = t12 - delta
           s32 = sqrt(dot_product(t32,t32))
           is32 = 1._dbl/s32
           dist1m = dot_product(t12_minus_delta,t32)
           dist2 = dot_product(t43,t32)
           u1m = t12_minus_delta - dist1m*t32*is32**2
           u2 = t43 - dist2*t32*is32**2
           w1m = sqrt(dot_product(u1m,u1m))
           w2 = sqrt(dot_product(u2,u2))
           iw1m = 1._dbl/w1m
           iw2 = 1._dbl/w2
           CALL force_torsions(is32,dist1m,dist2,u1m,u2,iw1m,iw2, &
                lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,em,fscalar)
           t12_plus_delta = t12 + delta
           s32 = sqrt(dot_product(t32,t32))
           is32 = 1._dbl/s32
           dist1p = dot_product(t12_plus_delta,t32)
           dist2 = dot_product(t43,t32)
           u1p = t12_plus_delta - dist1p*t32*is32**2
           u2 = t43 - dist2*t32*is32**2
           w1p = sqrt(dot_product(u1p,u1p))
           w2 = sqrt(dot_product(u2,u2))
           iw1p = 1._dbl/w1p
           iw2 = 1._dbl/w2
           CALL force_torsions(is32,dist1p,dist2,u1p,u2,iw1p,iw2, &
                lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,ep,fscalar)
           f_torsion_numer(idim,atom1) = f_torsion_numer(idim,atom1) + (em-ep)
           f_torsion_numer(idim,atom2) = f_torsion_numer(idim,atom2) - (em-ep)
           delta = 0.0_dbl
        END DO

! get the force on atom3 and atom2
        delta = 0.0_dbl
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t32_minus_delta = t32 - delta
           s32m = sqrt(dot_product(t32_minus_delta,t32_minus_delta))
           is32m = 1.0_dbl/s32m
           dist1m = dot_product(t12,t32_minus_delta)
           dist2m = dot_product(t43,t32_minus_delta)
           u1m = t12 - dist1m*t32_minus_delta*is32m**2
           u2m = t43 - dist2m*t32_minus_delta*is32m**2
           w1m = sqrt(dot_product(u1m,u1m))
           w2m = sqrt(dot_product(u2m,u2m))
           iw1m = 1._dbl/w1m
           iw2m = 1._dbl/w2m
           CALL force_torsions(is32m,dist1m,dist2m,u1m,u2m,iw1m,iw2m, &
                lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,em,fscalar)
           t32_plus_delta = t32 + delta
           s32p = sqrt(dot_product(t32_plus_delta,t32_plus_delta))
           is32p = 1.0_dbl/s32p
           dist1p = dot_product(t12,t32_plus_delta)
           dist2p = dot_product(t43,t32_plus_delta)
           u1p = t12 - dist1p*t32_plus_delta*is32p**2
           u2p = t43 - dist2p*t32_plus_delta*is32p**2
           w1p = sqrt(dot_product(u1p,u1p))
           w2p = sqrt(dot_product(u2p,u2p))
           iw1p = 1._dbl/w1p
           iw2p = 1._dbl/w2p
           CALL force_torsions(is32p,dist1p,dist2p,u1p,u2p,iw1p,iw2p, &
                lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,ep,fscalar)
           f_torsion_numer(idim,atom3) = f_torsion_numer(idim,atom3) + (em-ep)
           f_torsion_numer(idim,atom2) = f_torsion_numer(idim,atom2) - (em-ep)
           delta = 0.0_dbl
        END DO

! get the force on atom4 and atom3
        delta = 0.0_dbl
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t43_minus_delta = t43 - delta
           s32 = sqrt(dot_product(t32,t32))
           is32 = 1._dbl/s32
           dist1 = dot_product(t12,t32)
           dist2m = dot_product(t43_minus_delta,t32)
           u1 = t12 - dist1*t32*is32**2
           u2m = t43_minus_delta - dist2m*t32*is32**2
           w1 = sqrt(dot_product(u1,u1))
           w2m = sqrt(dot_product(u2m,u2m))
           iw1 = 1._dbl/w1
           iw2m = 1._dbl/w2m
           CALL force_torsions(is32,dist1,dist2m,u1,u2m,iw1,iw2m, &
                lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,em,fscalar)
           t43_plus_delta = t43 + delta
           s32 = sqrt(dot_product(t32,t32))
           is32 = 1._dbl/s32
           dist1 = dot_product(t12,t32)
           dist2p = dot_product(t43_plus_delta,t32)
           u1 = t12 - dist1*t32*is32**2
           u2p = t43_plus_delta - dist2p*t32*is32**2
           w1 = sqrt(dot_product(u1,u1))
           w2p = sqrt(dot_product(u2p,u2p))
           iw1 = 1._dbl/w1
           iw2p = 1._dbl/w2p
           CALL force_torsions(is32,dist1,dist2p,u1,u2p,iw1,iw2p, &
                lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,ep,fscalar)
           f_torsion_numer(idim,atom4) = f_torsion_numer(idim,atom4) + (em-ep)
           f_torsion_numer(idim,atom3) = f_torsion_numer(idim,atom3) - (em-ep)
           delta = 0.0_dbl
        END DO

        lltorsion => lltorsion % next
     END DO torsion
  END DO MOL

  f_torsion_numer = f_torsion_numer / ( 2.0_dbl * numerical_shift )

END SUBROUTINE force_torsion_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/force_nonbond_numer [1.0] *
!!
!!   NAME
!!     force_nonbond_numer
!!
!!   FUNCTION
!!     Computes the force and the potential of the minimum image,
!!     and the pressure tensor
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_nonbond_numer ( box, ewald_param, part, pnode, potparm,  &
                                 numerical_shift, f_numer, energy_numer )


    TYPE(cell_type), POINTER                 :: box
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(pair_potential_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    REAL(dbl), INTENT(IN)                    :: numerical_shift
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: f_numer
    REAL(dbl), INTENT(OUT)                   :: energy_numer

    CHARACTER(LEN=*), PARAMETER :: module_name = "fist_force_numer", &
      routine_name = "force_nonbond_numer", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER :: atom_a, atom_b, cell_b( 3 ), icell, id, iexclude, ikind, &
      inode, istat, jcell, jkind, jnode, kcell, natoms, ncell( 3 ), nexcl, &
      nnode, periodic( 3 )
    LOGICAL                                  :: match
    REAL(dbl) :: energy, energy_minus, energy_plus, qa, qb, rab2, rab2_max, &
      rab_max, rij_minus_delta_sq, rij_plus_delta_sq
    REAL(dbl), DIMENSION(3)                  :: delta, r, ra_pbc, rab, rb, &
                                                sab_max, sb, sb_pbc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(linklist_exclusion), POINTER        :: ll_exclude
    TYPE(pbc_coord_type)                     :: pbc_coord

!------------------------------------------------------------------------------

 energy_numer = 0.0_dbl
 f_numer = 0.0_dbl
 delta = 0.0_dbl

 CALL get_cell ( cell = box, periodic = periodic )

 natoms = SIZE ( part )
 ALLOCATE ( pbc_coord % s ( 3, natoms ), STAT = istat )
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s",&
                                  3*natoms*dp_size)
 ALLOCATE ( pbc_coord % r ( 3, natoms ), STAT = istat )
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%r" )

 nnode = SIZE ( pnode )

 DO inode=1,nnode
   atom_a = pnode ( inode ) % p % iatom
   ra_pbc ( : ) = pbc ( pnode ( inode ) % p % r ( : ), box )
   pbc_coord % r ( :, atom_a ) = ra_pbc ( : )
   pbc_coord % s ( :, atom_a ) = real_to_scaled ( ra_pbc ( : ), box )
 END DO

!
! starting the force loop
!
 ILOOP:   DO inode = 1, nnode
            atomic_kind => pnode ( inode ) % p % atomic_kind
            CALL get_atomic_kind ( atomic_kind = atomic_kind, &
                             kind_number=ikind, qeff=qa )
            atom_a = pnode ( inode ) % p % iatom
            nexcl = pnode ( inode ) % nexcl
            ll_exclude => pnode ( inode ) % ex
 JLOOP:     DO jnode = 1, nnode
              atom_b = pnode ( jnode ) % p % iatom
              IF ( atom_b < atom_a ) CYCLE
              sb_pbc(:) = pbc_coord % s ( :, atom_b )
              atomic_kind => pnode ( jnode ) % p % atomic_kind
              CALL get_atomic_kind ( atomic_kind = atomic_kind, &
                             kind_number=jkind, qeff=qb )
              rab2_max = potparm ( ikind, jkind ) % rcutsq 
              rab_max = SQRT ( potparm ( ikind, jkind ) % rcutsq )
              r(:) = rab_max
              sab_max ( : ) = real_to_scaled ( r ( : ), box )
              ncell ( : ) = ( INT ( sab_max ( : ) ) + 1 ) * periodic ( : )

              DO icell=-ncell(1),ncell(1)
                cell_b(1) = icell
               DO jcell=-ncell(2),ncell(2)
                 cell_b(2) = jcell
                 DO kcell=-ncell(3),ncell(3)
                   cell_b(3) = kcell

                   sb(:) = sb_pbc(:) + REAL(cell_b(:),dbl)
                   rb(:) = scaled_to_real(sb(:),box)

                   rab ( : ) = rb( : ) - pbc_coord % r ( :, atom_a )
                   rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

                   IF ( rab2 < rab2_max .AND. rab2 > 1.e-6_dbl ) THEN
                     ll_exclude => pnode ( inode ) % ex
                     match = .FALSE.
                     EXCL: DO iexclude = 1, nexcl
                       IF ( ll_exclude % p % iatom == atom_b ) THEN
                         match = .TRUE.
                         EXIT EXCL
                       END IF
                       ll_exclude => ll_exclude % next
                     END DO EXCL
                     NOTMATCH: IF ( .NOT. match ) THEN
                       delta = 0.0_dbl
                       DIM_LOOP: DO id = 1, 3
                         delta(id) = numerical_shift
                         rij_minus_delta_sq = DOT_PRODUCT(rab-delta,rab-delta)
                         rij_plus_delta_sq = DOT_PRODUCT(rab+delta,rab+delta)
                         IF (qa==0.AND.qb==0) THEN
                           CALL potential_f(rij_minus_delta_sq,potparm,qa,qb, &
                                ikind, jkind, energy_minus )
                           CALL potential_f(rij_plus_delta_sq,potparm,qa,qb, &
                              ikind, jkind, energy_plus )
                         ELSE
                           CALL potential_f(rij_minus_delta_sq,potparm,qa,qb, &
                              ikind, jkind, energy_minus, ewald_param )
                           CALL potential_f(rij_plus_delta_sq,potparm,qa,qb, &
                              ikind, jkind, energy_plus, ewald_param )
                         ENDIF
                           f_numer(id,atom_a) = f_numer(id,atom_a) + (energy_plus-energy_minus)
                           f_numer(id,atom_b) = f_numer(id,atom_b) - (energy_plus-energy_minus)
                           delta = 0.0_dbl
                       END DO DIM_LOOP
                       IF ( qa == 0 .AND. qb == 0 ) THEN
                        CALL potential_f ( rab2, potparm, qa, qb, &
                           ikind, jkind, energy )
                       ELSE
                         CALL potential_f ( rab2, potparm, qa, qb, &
                           ikind, jkind, energy, ewald_param )
                       ENDIF
                       energy_numer = energy_numer + energy
                     END IF NOTMATCH
                   ENDIF
                   END DO
                 END DO
               END DO
            END DO JLOOP
          END DO ILOOP

! computing long range corrections to the potential
! when we get the system for handling multiple potentials
! we will add the long range correction
! energy_numer=energy_numer+lrc*(1./box % deth)

 f_numer = f_numer/(2.0_dbl*numerical_shift)

 DEALLOCATE (pbc_coord%r,STAT=istat)
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%r")
 DEALLOCATE (pbc_coord%s,STAT=istat)
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s")

END SUBROUTINE force_nonbond_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/ewald_force_numer [1.0] *
!!
!!   NAME
!!     ewald_force_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_force_numer ( ewald_param, pnode, box, pw_grid,  &
                               numerical_shift, f_numer )


    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    REAL(dbl), INTENT(IN)                    :: numerical_shift
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_numer

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: sum_igr
    INTEGER                                  :: i, idim, isos, natoms, ngtot
    REAL(dbl)                                :: alpha, charge, em, ep, &
                                                epsilon0
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: gauss, r_delta
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!------------------------------------------------------------------------------
! allocating

  ngtot = pw_grid % ngpts_cut
  ALLOCATE (sum_igr(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'sum_igr', ngtot )
  ALLOCATE (gauss(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'gauss', ngtot )
  ALLOCATE (r_delta(3),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'r_delta', 3 )
! defining alpha and epsilon0
  alpha = ewald_param % alpha
  epsilon0 = ewald_param % eps0

! first initialize the arrays gauss and sum_igr
  CALL potential_g_numer ( ep, pnode, sum_igr, gauss, alpha,  &
                           pw_grid % g, ngtot )

! initializing numerical force
  f_numer( :, : ) = 0.0_dbl

! computing the numerical force on each atom
  natoms = size(pnode)
  DO i = 1, natoms
     atomic_kind => pnode(i) % p % atomic_kind
     CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=charge)
     r_delta ( : ) = pnode(i) % p % r ( : )
     DO idim = 1, 3
        r_delta(idim) = r_delta(idim) + numerical_shift
        CALL de_g_numer(ep, sum_igr, gauss, pnode(i) % p % r ,  &
                        r_delta, charge, pw_grid % g, ngtot)
        r_delta(idim) = r_delta(idim) - 2.0_dbl*numerical_shift
        CALL de_g_numer(em, sum_igr, gauss, pnode(i) % p % r ,  &
                        r_delta, charge, pw_grid % g, ngtot)
        f_numer(idim,i) = ( em - ep ) / epsilon0 / box % deth
        r_delta(idim) = r_delta(idim) + numerical_shift
     END DO
  END DO
  f_numer = f_numer /  2.0_dbl / numerical_shift

  DEALLOCATE (r_delta,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer','r_delta')
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer','sum_igr')
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'gauss' )

END SUBROUTINE ewald_force_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_bond_numer [1.0] *
!!
!!   NAME
!!     pv_bond_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_bond_numer ( molecule, pv_test, box, pnode, delta )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pv_test
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    REAL(dbl), INTENT(IN)                    :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(dbl)                                :: idelta, vm, vp
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: s, x
    REAL(dbl), DIMENSION(3, 3)               :: dvdh
    TYPE(cell_type)                          :: box_local

!------------------------------------------------------------------------------
! assigning the local variables

  natoms = size(pnode)
  ALLOCATE (x(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 'x', 3 * natoms )
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 's', 3 * natoms )
  DO i = 1, natoms
     x(:,i) = pnode(i) % p % r ( : )
     s(:,i) = matmul(box % h_inv,pnode(i) % p % r)
  END DO
  box_local = box

! Initializing pv_test
  pv_test = 0.0_dbl
! Defining the increments
  idelta = 1.0_dbl/(2.0_dbl*delta)
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv(box_local)

! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL potential_bond(molecule,x,vp)

! tweak again
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl*delta

! get inverse
        CALL get_hinv ( box_local )

! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL potential_bond(molecule,x,vm)

! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
        CALL get_hinv(box_local)
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*box_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (x,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 'x' )
  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 's' )

END SUBROUTINE pv_bond_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_bend_numer [1.0] *
!!
!!   NAME
!!     pv_bend_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_bend_numer(molecule,pv_test,box,pnode,delta)


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pv_test
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    REAL(dbl), INTENT(IN)                    :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(dbl)                                :: idelta, vm, vp
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: s, x
    REAL(dbl), DIMENSION(3, 3)               :: dvdh
    TYPE(cell_type)                          :: box_local

!------------------------------------------------------------------------------
! assigning the local variables

  natoms = size(pnode)
  ALLOCATE (x(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bend_numer', 'x', 3 * natoms )
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bend_numer', 's', 3 * natoms )
  DO i = 1, natoms
     x(:,i) = pnode(i) % p % r ( : )
     s(:,i) = matmul(box % h_inv,pnode(i) % p % r)
  END DO
  box_local = box

! Initializing pv_test
  pv_test = 0.0_dbl

! Defining the increments
  idelta = 1.0_dbl/(2.0_dbl*delta)
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
! get inverse
        CALL get_hinv(box_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL potential_bend(molecule,x,vp)
! tweak again
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl*delta
! get inverse
        CALL get_hinv(box_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL potential_bend(molecule,x,vm)
! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
        CALL get_hinv(box_local)
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*box_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (x,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bend_numer', 'x' )
  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bend_numer', 's' )

END SUBROUTINE pv_bend_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_torsion_numer [1.0] *
!!
!!   NAME
!!     pv_torsion_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_torsion_numer(molecule,pv_test,box,pnode,delta)


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pv_test
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    REAL(dbl), INTENT(IN)                    :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(dbl)                                :: idelta, vm, vp
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: s, x
    REAL(dbl), DIMENSION(3, 3)               :: dvdh
    TYPE(cell_type)                          :: box_local

!------------------------------------------------------------------------------
! assigning the local variables

  natoms = size(pnode)
  ALLOCATE (x(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_torsion_numer', 'x', 3 * natoms )
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_torsion_numer', 's', 3 * natoms )
  DO i = 1, natoms
     x(:,i) = pnode(i) % p % r ( : )
     s(:,i) = matmul(box % h_inv,pnode(i) % p % r)
  END DO
  box_local = box

! Initializing pv_test
  pv_test = 0.0_dbl

! Defining the increments
  idelta = 1.0_dbl/(2.0_dbl*delta)
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
! get inverse
        CALL get_hinv(box_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL potential_torsion(molecule,x,vp)
! tweak again
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl*delta
! get inverse
        CALL get_hinv(box_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL potential_torsion(molecule,x,vm)
! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
        CALL get_hinv(box_local)
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*box_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (x,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_torsion_numer', 'x' )
  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_torsion_numer', 's' )

END SUBROUTINE pv_torsion_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/ptens_numer [1.0] *
!!
!!   NAME
!!     ptens_numer
!!
!!   FUNCTION
!!     Computes the numerical pressure tensor
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ptens_numer(box,ewald_param,part,pnode,potparm,pv_test,delta)


    TYPE(cell_type), POINTER                 :: box
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(pair_potential_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    REAL(dbl), DIMENSION(3, 3), INTENT(OUT)  :: pv_test
    REAL(dbl), INTENT(IN)                    :: delta

    CHARACTER(LEN=*), PARAMETER :: module_name = "fist_force_numer", &
      routine_name = "ptens_numer", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: atom_a, i, ii, inode, istat, &
                                                j, jj, k, natoms, nnode
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(dbl)                                :: idelta, vm, vnb, vp
    REAL(dbl), DIMENSION(3)                  :: ra_pbc, sa
    REAL(dbl), DIMENSION(3, 3)               :: dvdh
    TYPE(pbc_coord_type)                     :: pbc_coord

!------------------------------------------------------------------------------

 CALL get_cell ( cell = box, periodic = periodic )
 box%orthorhombic = .FALSE.

 natoms = SIZE ( part )

 ALLOCATE ( pbc_coord % s ( 3, natoms ), STAT = istat )
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s",&
                                  3*natoms*dp_size)
 ALLOCATE ( pbc_coord % r ( 3, natoms ), STAT = istat )
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%r",&
                                  3*natoms*dp_size)

 nnode = SIZE ( pnode )

 DO inode=1,nnode
   atom_a = pnode ( inode ) % p % iatom
   ra_pbc ( : ) = pbc ( pnode ( inode ) % p % r ( : ), box )
   pbc_coord % r ( :, atom_a ) = ra_pbc ( : )
   pbc_coord % s ( :, atom_a ) = real_to_scaled ( ra_pbc ( : ), box )
 END DO

! Initializing pv_test
  pv_test = 0.0_dbl

! Defining the increments
  idelta = 1.0_dbl / (2.0_dbl * delta )
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        box % hmat(ii,jj) = box % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv ( box )
! get coords from perturbed hmat
        DO inode = 1, nnode
           atom_a = pnode ( inode ) % p % iatom
           sa ( : ) = pbc_coord % s ( :, atom_a )
           pbc_coord % r ( :, atom_a ) = scaled_to_real( sa ( : ), box )
        END DO
        CALL getv_ptens(ewald_param,pnode,pbc_coord,box,periodic,potparm,vnb)
        vp = vnb

! tweak again
        box % hmat(ii,jj) = box % hmat(ii,jj) - 2.0_dbl*delta
        CALL get_hinv(box)
        DO inode = 1, nnode
           atom_a = pnode ( inode ) % p % iatom
           sa ( : ) = pbc_coord % s ( :, atom_a )
           pbc_coord % r ( :, atom_a ) = scaled_to_real( sa ( : ), box )
        END DO
        CALL getv_ptens(ewald_param,pnode,pbc_coord,box,periodic,potparm,vnb)
        vm = vnb

! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta

! resetting the hmat
        box % hmat(ii,jj) = box % hmat(ii,jj) + delta
        CALL get_hinv(box)
        DO inode = 1, nnode
           atom_a = pnode ( inode ) % p % iatom
           sa ( : ) = pbc_coord % s ( :, atom_a )
           pbc_coord % r ( :, atom_a ) = scaled_to_real ( sa ( : ), box )
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*box % hmat(j,k)
        END DO
     END DO
  END DO

 DEALLOCATE (pbc_coord%r,STAT=istat)
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%r")
 DEALLOCATE (pbc_coord%s,STAT=istat)
 IF (istat /= 0) CALL stop_memory(routine,"pbc_coord%s")

END SUBROUTINE ptens_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pvg_numer [1.0] *
!!
!!   NAME
!!     pvg_numer
!!
!!   FUNCTION
!!     Computes the numerical pressure tensor
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pvg_numer(ewald_param,pnode,box,pw_grid,pv_test,delta)


    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pv_test
    REAL(dbl), INTENT(IN)                    :: delta

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: sum_igr
    INTEGER                                  :: gmax, i, ii, isos, j, jj, k, &
                                                natoms, ngtot, npts_s(3)
    REAL(dbl)                                :: alpha, cutoff, ecut, &
                                                epsilon0, idelta, vm, vp
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: gauss
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: s
    REAL(dbl), DIMENSION(3, 3)               :: dvdh
    TYPE(cell_type)                          :: box_local
    TYPE(particle_node_type), ALLOCATABLE, &
      DIMENSION(:)                           :: pnode_local
    TYPE(pw_grid_type)                       :: pw_grid_local

!------------------------------------------------------------------------------
! allocating

  natoms = SIZE ( pnode )
  ALLOCATE ( s ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 's', 3 * natoms )
  ngtot = pw_grid % ngpts_cut
  ALLOCATE (sum_igr(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 'sum_igr', ngtot )
  ALLOCATE (gauss(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 'gauss', ngtot )
  ALLOCATE (pnode_local(natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 'pnode_local', natoms )

! assigning the local variables
  box_local = box
  pnode_local = pnode
  pw_grid_local = pw_grid
  DO i = 1, natoms
     s(:,i) = matmul(box % h_inv,pnode(i) % p % r)
  END DO

! defining alpha and epsilon0 (see ewald.f for details)
  alpha = ewald_param % alpha
  epsilon0 = ewald_param % eps0

! Initializing pv_test
  pv_test = 0.0_dbl

  gmax = ewald_param % gmax
  npts_s = (/ gmax, gmax, gmax /)
  CALL pw_find_cutoff ( npts_s, box, cutoff )
  ecut = 0.5_dbl * cutoff * cutoff

! Defining the increments
  idelta = 1.0_dbl/(2.0_dbl*delta)
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv(box_local)

! recompute new partical postions from old scaled coords.
        DO i = 1, natoms
           pnode_local(i) % p % r = matmul(box_local % hmat,s(:,i))
        END DO

! compute the g-vectors
        CALL pw_grid_construct( pw_grid_local )
        CALL pw_grid_setup( box_local, pw_grid_local, ecut)

        CALL potential_g_numer ( vp, pnode_local, sum_igr, gauss, alpha, &
             pw_grid_local % g, ngtot )
        vp = vp / ( epsilon0 * box_local % deth )

! tweak again
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl * delta

! get inverse
        CALL get_hinv ( box_local )

! recompute new partical postions from old scaled coords.
        DO i = 1, natoms
           pnode_local(i) % p % r = matmul(box_local % hmat,s(:,i))
        END DO

! compute the g-vectors
        CALL pw_grid_setup( box_local, pw_grid_local, ecut)

        CALL potential_g_numer ( vm, pnode_local, sum_igr, gauss, alpha, &
             pw_grid_local % g, ngtot )
        vm = vm / epsilon0 / box_local % deth

        CALL pw_grid_destruct( pw_grid_local )

! calculate the derivative
        dvdh ( ii, jj ) = ( vm - vp ) * idelta

! resetting the hmat
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*box_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'pvg_numer', 's' )
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'pvg_numer', 'sum_igr' )
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'pvg_numer', 'gauss' )

END SUBROUTINE pvg_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_g_numer [1.0] *
!!
!!   NAME
!!     potential_g_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_g_numer(energy,pnode,sum_igr,gauss,alpha, &
     glocal,ngtot)


    REAL(dbl), INTENT(OUT)                   :: energy
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    COMPLEX(dbl), DIMENSION(:), INTENT(OUT)  :: sum_igr
    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: gauss
    REAL(dbl), INTENT(IN)                    :: alpha
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: glocal
    INTEGER, INTENT(IN)                      :: ngtot

    INTEGER                                  :: iatom, ig, natoms
    REAL(dbl)                                :: charge, gdotr, gsq
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!------------------------------------------------------------------------------
! numerical energy

  natoms = SIZE ( pnode )
  energy = 0.0_dbl
  sum_igr = 0.0_dbl
  DO ig = 1, ngtot
     gsq = DOT_PRODUCT(glocal(:,ig),glocal(:,ig))
     IF ( gsq <= 1.0E-10_dbl ) CYCLE
     gauss(ig) = exp(-gsq*0.25_dbl/alpha/alpha)/gsq
     DO iatom = 1, natoms
        atomic_kind => pnode(iatom) % p % atomic_kind
        CALL get_atomic_kind (atomic_kind=atomic_kind, qeff=charge)
        gdotr = DOT_PRODUCT(pnode(iatom) % p % r ( : ), glocal(:,ig))
        sum_igr(ig) = sum_igr(ig) + charge*CMPLX(COS(gdotr),SIN(gdotr),dbl)
     END DO

! computing the potential energy
     energy = energy + gauss(ig)* REAL ( sum_igr ( ig ) *  &
                                        CONJG( sum_igr ( ig ) ), dbl )
  END DO

END SUBROUTINE potential_g_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/de_g_numer [1.0] *
!!
!!   NAME
!!     de_g_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE de_g_numer(energy,sum_igr,gauss,r,r_delta,charge,glocal, &
     igtot)


    REAL(dbl), INTENT(OUT)                   :: energy
    COMPLEX(dbl), DIMENSION(:), INTENT(IN)   :: sum_igr
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: gauss, r, r_delta
    REAL(dbl), INTENT(IN)                    :: charge
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: glocal
    INTEGER, INTENT(IN)                      :: igtot

    COMPLEX(dbl)                             :: sum
    INTEGER                                  :: ig
    REAL(dbl)                                :: gdotr, gdotr_delta, gsq

!------------------------------------------------------------------------------
! initialize energy

  energy =  0._dbl
  sum = CMPLX( 0._dbl, 0._dbl, dbl )
  DO ig = 1, igtot
     gsq = DOT_PRODUCT( glocal(:,ig), glocal(:,ig))
     IF ( gsq <= 1.0E-10_dbl ) CYCLE
! compute g.r and g.(r+delta)
     gdotr = DOT_PRODUCT ( r ( : ), glocal ( :, ig ) )
     gdotr_delta = DOT_PRODUCT ( r_delta ( : ), glocal ( :, ig ) )
! subtract off exp(ig.r) and add exp(ig.(r+delta))
     sum = sum_igr ( ig ) - charge * &
           CMPLX ( COS( gdotr ), SIN( gdotr ), dbl ) + &
           charge * CMPLX( COS( gdotr_delta ), SIN( gdotr_delta ), dbl )
! recompute energy
     energy = energy + gauss(ig) * REAL ( sum * CONJG ( sum ), dbl )
  END DO

END SUBROUTINE de_g_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pme_force_numer [1.0] *
!!
!!   NAME
!!     pme_force_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pme_force_numer ( delta, f_part, dg, part, box, grid_s, grid_b,  &
                             ewald_param )
    REAL(dbl), INTENT(IN)                    :: delta
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: f_part
    TYPE(dg_type), INTENT(IN)                :: dg
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param

    INTEGER                                  :: i, id, natoms
    REAL(dbl)                                :: eminus, eplus

  natoms = SIZE ( part ) 
  DO i = 1, natoms 
    DO id = 1, 3
      part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta 
      CALL pme_energy_numer ( dg, part, box, grid_s, grid_b, ewald_param, eplus )
      part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2._dbl * delta 
      CALL pme_energy_numer ( dg, part, box, grid_s, grid_b, ewald_param, eminus )
      f_part ( id, i ) = eminus - eplus
      part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta 
    ENDDO
  ENDDO
  f_part = f_part / 2._dbl / delta
   
END SUBROUTINE pme_force_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pme_energy_numer [1.0] *
!!
!!   NAME
!!     pme_energy_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pme_energy_numer ( dg, part, box, grid_s, grid_b, ewald_param, energy )
    TYPE(dg_type), INTENT(IN)                :: dg
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(dbl), INTENT(OUT)                   :: energy

  CALL pme_evaluate ( dg, part, box, grid_s, grid_b, ewald_param, energy )
   
END SUBROUTINE pme_energy_numer
!!*****
!******************************************************************************
!!****** fist_force_numer/ewald_energy_numer [1.0] *
!!
!!   NAME
!!     ewald_energy_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_energy_numer ( ewald_param, pnode, box, pw_grid, &
     energy_numer, ngtot )


    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    REAL(dbl), INTENT(OUT)                   :: energy_numer
    INTEGER, INTENT(IN)                      :: ngtot

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: sum_igr
    INTEGER                                  :: isos
    REAL(dbl)                                :: alpha, epsilon0
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: gauss

!------------------------------------------------------------------------------
! allocating

  ALLOCATE ( sum_igr ( ngtot ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'sum_igr', ngtot )
  ALLOCATE ( gauss ( ngtot ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'gauss', ngtot )

! defining alpha and epsilon0
  alpha = ewald_param % alpha
  epsilon0 = ewald_param % eps0

! initializing numerical energy
  energy_numer = 0.0_dbl

! computing the energy
  CALL potential_g_numer(energy_numer,pnode,sum_igr,gauss,alpha, &
       pw_grid % g, ngtot)
  energy_numer = energy_numer/epsilon0/box % deth
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'sum_igr' )
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'gauss' )

END SUBROUTINE ewald_energy_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_bond [1.0] *
!!
!!   NAME
!!     potential_bond
!!
!!   FUNCTION
!!     Computes analytically the potential from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_bond ( molecule, x, energy )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: x
    REAL(dbl), INTENT(OUT)                   :: energy

    INTEGER                                  :: ibond, imol
    REAL(dbl)                                :: e0, fscalar
    REAL(dbl), DIMENSION(3)                  :: rij
    TYPE(linklist_bonds), POINTER            :: llbond

!------------------------------------------------------------------------------
! initializing

  energy = 0.0_dbl
  MOL: DO imol = 1, SIZE ( molecule )
     llbond => molecule ( imol ) % ll_bonds
     BOND: DO ibond = 1, molecule ( imol ) % nbonds_mol
        rij = x ( :, llbond % index ( 1 ) ) - x ( :, llbond % index ( 2 ) )
        CALL force_bonds ( rij, llbond % bond_param % r0, &
             llbond % bond_param % k, e0, fscalar )
        energy = energy + e0
        llbond => llbond % next
     END DO BOND
  END DO MOL

END SUBROUTINE potential_bond

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_bend [1.0] *
!!
!!   NAME
!!     potential_bend
!!
!!   FUNCTION
!!     Computes analytically the potential from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_bend ( molecule, x, energy )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: x
    REAL(dbl), INTENT(OUT)                   :: energy

    INTEGER                                  :: ibend, imol
    REAL(dbl)                                :: d12, d32, dist, e0, fscalar, &
                                                id12, id32, theta
    REAL(dbl), DIMENSION(3)                  :: b12, b32, g1, g2, g3
    TYPE(linklist_bends), POINTER            :: llbend

!------------------------------------------------------------------------------
! initializing

  energy = 0.0_dbl
  MOL: DO imol = 1, SIZE ( molecule )
     llbend => molecule ( imol ) % ll_bends
     BEND: DO ibend = 1, molecule ( imol ) % nbends_mol
        b12 = x(:,llbend % index(1)) - x(:,llbend % index(2))
        b32 = x(:,llbend % index(3)) - x(:,llbend % index(2))
        d12 = SQRT ( DOT_PRODUCT ( b12, b12 ) )
        d32 = SQRT ( DOT_PRODUCT ( b32, b32 ) )
        id12 = 1.0_dbl / d12
        id32 = 1.0_dbl / d32
        dist = DOT_PRODUCT ( b12, b32 )
        theta = ACOS ( dist * id12 * id32 )
        CALL force_bends ( b12, b32, d12, d32, id12, id32, dist, theta, &
             llbend % bend_param % theta0, llbend % bend_param % k, &
             g1, g2, g3, e0, fscalar )
        energy = energy + e0
        llbend => llbend % next
     END DO BEND
  END DO MOL

END SUBROUTINE potential_bend

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_torsion [1.0] *
!!
!!   NAME
!!     potential_torsion
!!
!!   FUNCTION
!!     Computes analytically the potential from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_torsion ( molecule, x, energy )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: x
    REAL(dbl), INTENT(OUT)                   :: energy

    INTEGER                                  :: imol, itorsion
    REAL(dbl)                                :: dist1, dist2, e0, fscalar, &
                                                is32, iw1, iw2, s32, w1, w2
    REAL(dbl), DIMENSION(3)                  :: gt1, gt2, gt3, gt4, t12, t32, &
                                                t43, u1, u2
    TYPE(linklist_torsions), POINTER         :: lltorsion

!------------------------------------------------------------------------------
! initializing

  energy = 0.0_dbl
  MOL: DO imol = 1, SIZE ( molecule )
     lltorsion => molecule ( imol ) % ll_torsions
     torsion: DO itorsion = 1, molecule ( imol ) % ntorsions_mol
        t12 = x(:,lltorsion % index(1)) - x(:,lltorsion % index(2))
        t32 = x(:,lltorsion % index(3)) - x(:,lltorsion % index(2))
        t43 = x(:,lltorsion % index(4)) - x(:,lltorsion % index(3))
        s32 = SQRT ( DOT_PRODUCT ( t32, t32 ) )
        is32 = 1.0_dbl / s32
        dist1 = DOT_PRODUCT ( t12, t32 )
        dist2 = DOT_PRODUCT ( t43, t32 )
        u1 = t12 - dist1 * t32 * is32**2
        u2 = t43 - dist2 * t32 * is32**2
        w1 = SQRT ( DOT_PRODUCT ( u1, u1 ) )
        w2 = SQRT ( DOT_PRODUCT ( u2, u2 ) )
        iw1 = 1.0_dbl / w1
        iw2 = 1.0_dbl / w2
        CALL force_torsions ( is32, dist1, dist2, u1, u2, iw1, iw2, &
             lltorsion % torsion_param % a, gt1, gt2, gt3, gt4, e0, fscalar)
        energy = energy + e0
        lltorsion => lltorsion % next
     END DO torsion
  END DO MOL

END SUBROUTINE potential_torsion

!!*****
!******************************************************************************
!!****** fist_force_numer/getv_ptens [1.0] *
!!
!!   NAME
!!     getv_ptens
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE getv_ptens ( ewald_param, pnode, pbc_coord, box, periodic, potparm, vnb )


    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(pbc_coord_type), INTENT(IN)         :: pbc_coord
    TYPE(cell_type), POINTER                 :: box
    INTEGER, DIMENSION(:), INTENT(IN)        :: periodic
    TYPE(pair_potential_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    REAL(dbl), INTENT(OUT)                   :: vnb

    INTEGER :: atom_a, atom_b, cell_b( 3 ), icell, iexclude, ikind, inode, &
      jcell, jkind, jnode, kcell, ncell( 3 ), nexcl, nnode
    LOGICAL                                  :: match
    REAL(dbl)                                :: energy, qa, qb, rab2, &
                                                rab2_max, rab_max
    REAL(dbl), DIMENSION(3)                  :: r, rab, rb, sab_max, sb, &
                                                sb_pbc
    TYPE(linklist_exclusion), POINTER        :: ll_exclude

!------------------------------------------------------------------------------
!
! initializing the potential energy
!

  vnb = 0.0_dbl
!
! starting the loop
!
 nnode = SIZE ( pnode )
 ILOOP:   DO inode = 1, nnode
            CALL get_atomic_kind ( pnode ( inode ) % p % atomic_kind, &
                              kind_number = ikind, qeff=qa )
            atom_a = pnode ( inode ) % p % iatom
            nexcl = pnode ( inode ) % nexcl
            ll_exclude => pnode ( inode ) % ex
 JLOOP:     DO jnode = 1, nnode
              atom_b = pnode ( jnode ) % p % iatom
              IF ( atom_b < atom_a ) CYCLE
              sb_pbc(:) = pbc_coord % s ( :, atom_b )
              CALL get_atomic_kind ( pnode ( jnode ) % p % atomic_kind, &
                                  kind_number = jkind, qeff=qb )
              rab2_max = potparm ( ikind, jkind ) % rcutsq 
              rab_max = SQRT ( potparm ( ikind, jkind ) % rcutsq )
              r(:) = rab_max
              sab_max ( : ) = real_to_scaled ( r ( : ), box )
              ncell ( : ) = ( INT ( sab_max ( : ) ) + 1 ) * periodic ( : )

              DO icell=-ncell(1),ncell(1)
                cell_b(1) = icell
               DO jcell=-ncell(2),ncell(2)
                 cell_b(2) = jcell
                 DO kcell=-ncell(3),ncell(3)
                   cell_b(3) = kcell

                   sb(:) = sb_pbc(:) + REAL(cell_b(:),dbl)
                   rb(:) = scaled_to_real(sb(:),box)

                   rab ( : ) = rb( : ) - pbc_coord % r ( :, atom_a )
                   rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

                   IF ( rab2 < rab2_max .AND. rab2 > 1.e-6_dbl ) THEN
                     ll_exclude => pnode ( inode ) % ex
                     match = .FALSE.
                     EXCL: DO iexclude = 1, nexcl
                       IF ( ll_exclude % p % iatom == atom_b ) THEN
                         match = .TRUE.
                         EXIT EXCL
                       END IF
                       ll_exclude => ll_exclude % next
                     END DO EXCL
                     NOTMATCH: IF ( .NOT. match ) THEN
                         IF ( qa == 0 .AND. qb == 0 ) THEN
                          CALL potential_f ( rab2, potparm, qa, qb, &
                             ikind, jkind, energy )
                         ELSE
                           CALL potential_f ( rab2, potparm, qa, qb, &
                             ikind, jkind, energy, ewald_param )
                         ENDIF
                         vnb = vnb + energy
                     END IF NOTMATCH
                   ENDIF
                 END DO
               END DO
             END DO
          END DO JLOOP
        END DO ILOOP

END SUBROUTINE getv_ptens

!!*****
!******************************************************************************

END MODULE fist_force_numer

!******************************************************************************
