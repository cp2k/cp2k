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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fist_force_numer
  
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_nonbond_force, ONLY : find_image
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi
  USE molecule_types, ONLY : molecule_structure_type, linklist_bonds, &
       linklist_bends, linklist_exclusion, particle_node_type
  USE mol_force, ONLY : force_bonds, force_bends
  USE pair_potential, ONLY : potential_f, potentialparm_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_find_cutoff, pw_grid_setup
  USE simulation_cell, ONLY : cell_type, pbc, get_hinv, get_cell_param
  USE stop_program, ONLY : stop_memory
  
  PRIVATE
  PUBLIC :: force_bend_numer, force_bond_numer, force_nonbond_numer
  PUBLIC :: force_recip_numer, pvbond_numer, pvbend_numer, ptens_numer
  PUBLIC :: pvg_numer, de_g_numer, potential_g_numer, energy_recip_numer
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE force_bond_numer ( molecule, f_bond_numer, v_bond_numer, &
     numerical_shift )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: f_bond_numer
  REAL ( dbl ), INTENT ( OUT ) :: v_bond_numer
  REAL ( dbl ), INTENT ( IN ) :: numerical_shift
  
! Locals
  TYPE ( linklist_bonds ), POINTER :: llbond
  INTEGER :: imol, ibond, idim, atom1, atom2
  REAL ( dbl ), DIMENSION (3) :: rij, rij_minus_delta, rij_plus_delta
  REAL ( dbl ) :: e0, ep, em, fscalar
  REAL ( dbl ), DIMENSION (3) :: delta
  
!------------------------------------------------------------------------------
  
! initializing the forces
  f_bond_numer = 0.0_dbl
  v_bond_numer = 0.0_dbl
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
        v_bond_numer = v_bond_numer + e0
        
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

!******************************************************************************

SUBROUTINE force_bend_numer ( molecule, f_bend_numer, v_bend_numer, &
     numerical_shift )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: f_bend_numer
  REAL ( dbl ), INTENT ( OUT ) :: v_bend_numer
  REAL ( dbl ), INTENT ( IN ) :: numerical_shift
  
! Locals
  TYPE ( linklist_bends ), POINTER :: llbend
  INTEGER :: idim, imol, ibend, atom1, atom2, atom3
  REAL ( dbl ) :: d12, d12p, d12m, id12, id12p, id12m
  REAL ( dbl ) :: d32, d32p, d32m, id32, id32p, id32m
  REAL ( dbl ) :: dist, theta
  REAL ( dbl ) :: em, ep, e0, fscalar
  REAL ( dbl ), DIMENSION (3) :: g1, g2, g3
  REAL ( dbl ), DIMENSION (3) :: b12, b12_plus_delta, b12_minus_delta
  REAL ( dbl ), DIMENSION (3) :: b32, b32_plus_delta, b32_minus_delta
  REAL ( dbl ), DIMENSION (3) :: delta
  
!------------------------------------------------------------------------------
  
! initializing the forces
  f_bend_numer = 0.0_dbl
  v_bend_numer = 0.0_dbl
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
        v_bend_numer = v_bend_numer + e0
        
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

!******************************************************************************

SUBROUTINE force_nonbond_numer ( ewald_param, pnode, box, potparm,  &
     numerical_shift, f_numer, energy_numer )
  
!
! Calculates the force and the potential of the minimum image, and
! the pressure tensor
!
  IMPLICIT NONE
  
! Arguments
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( cell_type ), INTENT ( INOUT ) :: box
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( IN ) :: potparm
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_numer
  REAL ( dbl ), INTENT ( OUT ) :: energy_numer
  REAL ( dbl ), INTENT ( IN ) :: numerical_shift
  
! Locals
  INTEGER :: i, j, ii, jj, iatomtype, jatomtype
  INTEGER :: icell,jcell,kcell,itype,jtype,id
  INTEGER :: iexclude, natoms, n, isos
  INTEGER :: xx, yy, zz
  REAL ( dbl ), DIMENSION (3) ::vec, ri, rij, cell, ang, perd, dist
  REAL ( dbl ) :: energy, fscalar, rijsq, rcut
  REAL ( dbl ) :: rij_minus_delta_sq, rij_plus_delta_sq
  REAL ( dbl ), DIMENSION (3) :: delta, s, quotient
  REAL ( dbl ) :: energy_plus, energy_minus, qi, qj
  INTEGER, DIMENSION ( :, :, : ), ALLOCATABLE :: n_images
  TYPE ( linklist_exclusion ), POINTER :: ll_exclude
  LOGICAL :: match
  
!------------------------------------------------------------------------------
  
  n = SIZE ( potparm ( :, 1 ) )
  IF ( .NOT. ALLOCATED ( n_images ) ) &
       ALLOCATE ( n_images ( n, n, 3 ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_nonbond_numer', 'n_images', n * n * 3 )
  
  CALL get_cell_param ( box, cell )
  
  DO itype=1, n
     DO jtype=1, n
        rcut = SQRT ( potparm ( itype, jtype ) % rcutsq )
        quotient ( : ) = rcut / cell ( : )
        DO id=1,3
           IF ( quotient ( id ) <= 0.5_dbl ) THEN
              n_images ( itype, jtype, id )=0
           ELSE
              n_images ( itype, jtype, id ) = CEILING ( rcut / cell ( id ) )
           END IF
        END DO
     END DO
  END DO
  
  energy_numer = 0.0_dbl
  f_numer = 0.0_dbl
  
  delta = 0.0_dbl
!
! starting the force loop
!
  CALL get_hinv ( box )
  CALL get_cell_param(box,cell,ang)
  perd=box % perd
  natoms = size(pnode)
  ILOOP:  DO i = 1, natoms
     iatomtype = pnode(i) % p % prop % ptype
     qi = pnode(i) % p % prop % charge
     ri = pnode(i) % p % r
     JLOOP: DO j = 1, natoms
        IF (j>i) THEN
           jatomtype = pnode(j) % p % prop % ptype
           qj = pnode(j) % p % prop % charge
           rij = pnode(j) % p % r - ri
           
!
! imaging,transform back to unscaled distances and get potential and force
!
           IF (MAXVAL(n_images(iatomtype,jatomtype,:))==0) THEN
              rij = pbc(rij,box)
           END IF
           rijsq = dot_product(rij,rij)
           
           IF ( rijsq <= potparm ( iatomtype, jatomtype ) % rcutsq ) THEN
! pointing to the start of the exclusion list of particle i
              ll_exclude => pnode(i) % ex
              match = .FALSE.
              EXCL: DO iexclude = 1, pnode(i) % nexcl
                 IF (ll_exclude % p % iatom==j) THEN
                    match = .TRUE.
                    EXIT EXCL
                 END IF
                 ll_exclude => ll_exclude % next
              END DO EXCL
              NOTMATCH: IF ( .NOT. match) THEN
                 delta = 0.0_dbl
                 DIM_LOOP: DO id = 1, 3
                    delta(id) = numerical_shift
                    rij_minus_delta_sq = dot_product(rij-delta,rij-delta)
                    rij_plus_delta_sq = dot_product(rij+delta,rij+delta)
                    IF (qi==0.AND.qj==0) THEN
                      CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
                         iatomtype, jatomtype, energy_minus )
                      CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
                         iatomtype, jatomtype, energy_plus )
                    ELSE
                      CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
                         iatomtype, jatomtype, energy_minus, ewald_param )
                      CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
                         iatomtype, jatomtype, energy_plus, ewald_param )
                    ENDIF
                    f_numer(id,i) = f_numer(id,i) + (energy_plus-energy_minus)
                    f_numer(id,j) = f_numer(id,j) - (energy_plus-energy_minus)
                    delta = 0.0_dbl
                 END DO DIM_LOOP
                 IF (qi==0.AND.qj==0) THEN
                   CALL potential_f ( rijsq, potparm, qi, qj, &
                      iatomtype, jatomtype, energy )
                 ELSE
                   CALL potential_f ( rijsq, potparm, qi, qj, &
                      iatomtype, jatomtype, energy, ewald_param )
                 ENDIF
                 
                 energy_numer = energy_numer + energy
              END IF NOTMATCH
           END IF
           
           zz = n_images ( iatomtype, jatomtype, 3 )
           yy = n_images ( iatomtype, jatomtype, 2 )
           xx = n_images ( iatomtype, jatomtype, 1 )
           s(1) = box % h_inv(1,1)*rij(1) &
                + box % h_inv(1,2)*rij(2) &
                + box % h_inv(1,3)*rij(3)
           s(2) = box % h_inv(2,1)*rij(1) &
                + box % h_inv(2,2)*rij(2) &
                + box % h_inv(2,3)*rij(3)
           s(3) = box % h_inv(3,1)*rij(1) &
                + box % h_inv(3,2)*rij(2) &
                + box % h_inv(3,3)*rij(3)
           
           DO kcell = -zz, zz
              DO jcell = -yy, yy
                 DO icell = -xx, xx
                    IF ( icell == 0 .AND. jcell == 0 .AND. kcell == 0 ) CYCLE
                    
                    vec ( : ) = (/ icell, jcell, kcell /)
                    CALL find_image(s,perd,vec,box % hmat,rijsq,rij)
                    IF ( rijsq <= potparm(iatomtype,jatomtype ) % rcutsq ) THEN
                       delta = 0.0_dbl
                       DIM_LOOP2: DO id = 1, 3
                          delta(id) = numerical_shift
                          rij_minus_delta_sq = dot_product(rij-delta,rij-delta)
                          rij_plus_delta_sq = dot_product(rij+delta,rij+delta)
                          IF (qi==0.AND.qj==0) THEN
                            CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_minus )
                            CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_plus )
                          ELSE
                            CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_minus, ewald_param )
                            CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_plus, ewald_param )
                          ENDIF
                          f_numer(id,i) = f_numer(id,i) &
                               + (energy_plus-energy_minus)
                          f_numer(id,j) = f_numer(id,j) &
                               - (energy_plus-energy_minus)
                          delta = 0.0_dbl
                       END DO DIM_LOOP2
                       
                       IF (qi==0.AND.qj==0) THEN
                         CALL potential_f ( rijsq, potparm, qi, qj, &
                            iatomtype, jatomtype, energy )
                       ELSE
                         CALL potential_f ( rijsq, potparm, qi, qj, &
                            iatomtype, jatomtype, energy, ewald_param )
                       ENDIF
                       energy_numer = energy_numer + energy
                    END IF
                 END DO
              END DO
           END DO
           
        ELSE IF (i==j) THEN
! search images if needed
           zz=n_images(iatomtype,iatomtype,3)
           yy=n_images(iatomtype,iatomtype,2)
           xx=n_images(iatomtype,iatomtype,1)
           s = 0.0_dbl
           DO kcell=-zz,zz
              DO jcell=-yy,yy
                 DO icell=-xx,xx
                    IF (icell==0.AND.jcell==0.AND.kcell==0) CYCLE
                    vec ( : )=(/icell,jcell,kcell/)
                    CALL find_image(s,perd,vec,box % hmat,rijsq,rij)
                    IF ( rijsq <= potparm ( iatomtype, iatomtype ) % rcutsq ) &
                         THEN
                       delta = 0.0_dbl
                       DIM_LOOP3: DO id = 1, 3
                          delta(id) = numerical_shift
                          rij_minus_delta_sq = dot_product(rij-delta,rij-delta)
                          rij_plus_delta_sq = dot_product(rij+delta,rij+delta)
                          IF (qi==0.AND.qj==0) THEN
                            CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_minus )
                            CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_plus )
                          ELSE
                            CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_minus, ewald_param )
                            CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
                               iatomtype, jatomtype, energy_plus, ewald_param )
                          ENDIF
                          f_numer(id,i) = f_numer(id,i) &
                               + (energy_plus-energy_minus)
                          f_numer(id,i) = f_numer(id,i) &
                               - (energy_plus-energy_minus)
                          delta = 0.0_dbl
                       END DO DIM_LOOP3
                       
                       IF (qi==0.AND.qj==0) THEN
                         CALL potential_f ( rijsq, potparm, qi, qj, &
                            iatomtype, jatomtype, energy )
                       ELSE
                         CALL potential_f ( rijsq, potparm, qi, qj, &
                            iatomtype, jatomtype, energy, ewald_param )
                       ENDIF
                       energy_numer = energy_numer + energy
                    END IF
                 END DO
              END DO
           END DO
        END IF
        
     END DO JLOOP
  END DO ILOOP
  
! computing long range corrections to the potential
!when we get the system for handling multiple potentials
!we will add the long range correction
!        energy_numer=energy_numer+lrc*(1./box % deth)
  
  f_numer = f_numer/(2.0_dbl*numerical_shift)
  
END SUBROUTINE force_nonbond_numer

!******************************************************************************

SUBROUTINE force_recip_numer ( ewald_param, pnode, box, pw_grid,  &
                               numerical_shift, f_numer )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( pw_grid_type ), INTENT ( IN ) :: pw_grid
  REAL ( dbl ), INTENT ( IN ) :: numerical_shift
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: f_numer
  
! Locals
  COMPLEX ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: sum_igr
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: gauss
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: r_delta
  REAL ( dbl ) :: alpha, epsilon0, ep, em, charge
  INTEGER :: gpt, i, idim, natoms, ngtot, isos, lp, mp, np
  
!------------------------------------------------------------------------------

! allocating
  ngtot = pw_grid % ngpts_cut
  ALLOCATE (sum_igr(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'force_recip_numer', 'sum_igr', ngtot )
  ALLOCATE (gauss(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'force_recip_numer', 'gauss', ngtot )
  ALLOCATE (r_delta(3),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'force_recip_numer', 'r_delta', 3 )
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
     charge = pnode(i) % p % prop % charge
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
  IF ( isos /= 0 ) CALL stop_memory ( 'force_recip_numer','r_delta')
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'force_recip_numer','sum_igr')
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'force_recip_numer', 'gauss' )
  
END SUBROUTINE force_recip_numer

!******************************************************************************

SUBROUTINE pvbond_numer ( molecule, pv_test, box, pnode, delta )
  
! computes the numerical pressure tensor
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: pv_test
  REAL ( dbl ), INTENT ( IN ) :: delta
  
! Locals
  TYPE ( cell_type ) :: box_local
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: x, s
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: dvdh
  REAL ( dbl ) :: idelta, vm, vp, vnb
  INTEGER :: i, j, k, ii, jj, natoms, isos
  
!------------------------------------------------------------------------------
  
! assigning the local variables
  natoms = size(pnode)
  ALLOCATE (x(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbond_numer', 'x', 3 * natoms )
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbond_numer', 's', 3 * natoms )
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
        CALL vnumer_bonds(molecule,x,vp)
        
! tweak again 
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl*delta
        
! get inverse
        CALL get_hinv ( box_local )
        
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL vnumer_bonds(molecule,x,vm)
        
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
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbond_numer', 'x' )
  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbond_numer', 's' )
  
END SUBROUTINE pvbond_numer

!******************************************************************************

SUBROUTINE pvbend_numer(molecule,pv_test,box,pnode,delta)
  
! computes the numerical pressure tensor
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), INTENT ( IN ), DIMENSION ( : ) :: molecule
  TYPE ( particle_node_type ), INTENT ( IN ), DIMENSION ( : ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: pv_test
  REAL ( dbl ), INTENT ( IN ) :: delta
  
! Locals
  TYPE (cell_type) :: box_local
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: x, s
  REAL ( dbl ), DIMENSION (3,3) :: dvdh
  REAL ( dbl ) :: idelta, vm, vp, vnb
  INTEGER :: i, j, k, ii, jj, natoms, isos
  
!------------------------------------------------------------------------------
  
! assigning the local variables
  natoms = size(pnode)
  ALLOCATE (x(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbend_numer', 'x', 3 * natoms )
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbend_numer', 's', 3 * natoms )
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
        CALL vnumer_bends(molecule,x,vp)
! tweak again 
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl*delta
! get inverse
        CALL get_hinv(box_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL vnumer_bends(molecule,x,vm)
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
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbend_numer', 'x' )
  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvbend_numer', 's' )
  
END SUBROUTINE pvbend_numer

!******************************************************************************

SUBROUTINE ptens_numer(ewald_param,pnode,box,potparm,pv_test,delta)
  
! computes the numerical pressure tensor
  IMPLICIT NONE
  
! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( IN ) :: delta
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ) :: pv_test
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( IN ) :: potparm
  
! Locals
  TYPE ( cell_type ) :: box_local
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: x, s
  REAL ( dbl ), DIMENSION (3,3) :: dvdh
  REAL ( dbl ) :: idelta, vm, vp, vnb
  INTEGER :: i, j, k, ii, jj, natoms, isos
  
!------------------------------------------------------------------------------
  
! assigning the local variables
  natoms = size(pnode)
  ALLOCATE (x(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ptens_numer', 'x', 3 * natoms )
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ptens_numer', 's', 3 * natoms )
  
  DO i = 1, natoms
     x(:,i) = pnode(i) % p % r ( : )
     s(:,i) = matmul(box % h_inv,pnode(i) % p % r)
  END DO
  box_local = box
  
! Initializing pv_test
  pv_test = 0.0_dbl
  
! Defining the increments
  idelta = 1.0_dbl / (2.0_dbl * delta )
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) + delta
        
! get inverse
        CALL get_hinv(box_local)
        
! get coords from perturbed hmat
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL getv_ptens(ewald_param,pnode,x,vnb,box_local,potparm)
        vp = vnb
        
! tweak again 
        box_local % hmat(ii,jj) = box_local % hmat(ii,jj) - 2.0_dbl*delta
        CALL get_hinv(box_local)
        DO i = 1, natoms
           x(:,i) = matmul(box_local % hmat,s(:,i))
        END DO
        CALL getv_ptens(ewald_param,pnode,x,vnb,box_local,potparm)
        vm = vnb
        
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
  IF ( isos /= 0 ) CALL stop_memory( 'ptens_numer', 'x' )
  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'ptens_numer', 's' )
  
END SUBROUTINE ptens_numer

!******************************************************************************

SUBROUTINE pvg_numer(ewald_param,pnode,box,pw_grid,pv_test,delta)
  
! computes the numerical pressure tensor
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( pw_grid_type ), INTENT ( IN ) :: pw_grid
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( IN ) :: delta
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: pv_test
  
! Locals
  TYPE ( pw_grid_type ) :: pw_grid_local
  TYPE ( cell_type ) :: box_local
  TYPE ( particle_node_type ), DIMENSION ( : ), ALLOCATABLE :: pnode_local
  REAL ( dbl ), DIMENSION (3,3) :: dvdh
  REAL ( dbl ) :: idelta, vm, vp
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: s
  REAL ( dbl ) :: alpha, epsilon0, cutoff
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: gauss
  COMPLEX ( dbl ), DIMENSION ( : ), ALLOCATABLE :: sum_igr
  INTEGER :: i, j, k, ii, jj, natoms, isos, gmax, npts_s(3), ngtot
  
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
        CALL pw_grid_setup( box_local, pw_grid_local, cutoff)
        
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
        CALL pw_grid_setup( box_local, pw_grid_local, cutoff)

        CALL potential_g_numer ( vm, pnode_local, sum_igr, gauss, alpha, &
             pw_grid_local % g, ngtot )
        vm = vm / epsilon0 / box_local % deth
        
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

!******************************************************************************

SUBROUTINE potential_g_numer(energy,pnode,sum_igr,gauss,alpha, &
     glocal,ngtot)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  REAL ( dbl ), INTENT ( OUT ) :: energy
  REAL ( dbl ), INTENT ( IN ) :: alpha
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: glocal
  COMPLEX ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: sum_igr
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: gauss
  INTEGER, INTENT ( IN ) :: ngtot
  
! Locals
  REAL ( dbl ) :: gdotr, gsq, charge
  INTEGER :: iatom, ig, natoms
  
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
        gdotr = DOT_PRODUCT(pnode(iatom) % p % r ( : ), glocal(:,ig))
        charge = pnode(iatom) % p % prop % charge
        sum_igr(ig) = sum_igr(ig) + charge*CMPLX(COS(gdotr),SIN(gdotr))
     END DO
     
! computing the potential energy
     energy = energy + gauss(ig)* REAL ( sum_igr ( ig ) *  &
                                        CONJG( sum_igr ( ig ) ) )
  END DO
  
END SUBROUTINE potential_g_numer

!******************************************************************************

SUBROUTINE de_g_numer(energy,sum_igr,gauss,r,r_delta,charge,glocal, &
     igtot)
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ) :: charge
  COMPLEX ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: sum_igr
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: gauss
  REAL ( dbl ), INTENT ( OUT ) :: energy
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: glocal
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: r, r_delta
  INTEGER, INTENT ( IN ) :: igtot
  
! Locals
  INTEGER :: ig
  REAL ( dbl ) :: gdotr, gdotr_delta, gsq
  COMPLEX ( dbl ) :: sum
  
!------------------------------------------------------------------------------
  
! initialize energy
  energy =  0._dbl
  sum = ( 0._dbl, 0._dbl )
  DO ig = 1, igtot
     gsq = DOT_PRODUCT( glocal(:,ig), glocal(:,ig))
     IF ( gsq <= 1.0E-10_dbl ) CYCLE
! compute g.r and g.(r+delta)
     gdotr = DOT_PRODUCT ( r ( : ), glocal ( :, ig ) )
     gdotr_delta = DOT_PRODUCT ( r_delta ( : ), glocal ( :, ig ) )
! subtract off exp(ig.r) and add exp(ig.(r+delta))
     sum = sum_igr ( ig ) - charge * CMPLX ( COS( gdotr ), SIN( gdotr ) ) + &
          charge * CMPLX( COS( gdotr_delta ), SIN( gdotr_delta ) )
! recompute energy
     energy = energy + gauss(ig) * REAL ( sum * CONJG ( sum ) )
  END DO
  
END SUBROUTINE de_g_numer

!******************************************************************************

SUBROUTINE energy_recip_numer ( ewald_param, pnode, box, pw_grid, &
     energy_numer, ngtot )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( pw_grid_type ), INTENT ( IN ):: pw_grid
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( OUT ) :: energy_numer
  INTEGER, INTENT ( IN ) :: ngtot
  
! Locals
  COMPLEX ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: sum_igr
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: gauss
  REAL ( dbl ) :: alpha, epsilon0
  INTEGER :: ig, i, idim, isos, lp, mp, np
  
!------------------------------------------------------------------------------

! allocating

  ALLOCATE ( sum_igr ( ngtot ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'energy_recip_numer', 'sum_igr', ngtot )
  ALLOCATE ( gauss ( ngtot ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'energy_recip_numer', 'gauss', ngtot )

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
  IF ( isos /= 0 ) CALL stop_memory ( 'energy_recip_numer', 'sum_igr' )
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'energy_recip_numer', 'gauss' )
  
END SUBROUTINE energy_recip_numer

!******************************************************************************

SUBROUTINE vnumer_bonds ( molecule, x, energy )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( OUT ) :: energy
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: x
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  
! Locals
  REAL ( dbl ), DIMENSION ( 3 ) :: rij
  REAL ( dbl ) :: e0, fscalar
  INTEGER :: imol, ibond
  TYPE ( linklist_bonds ), POINTER :: llbond
  
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
  
END SUBROUTINE vnumer_bonds

!******************************************************************************

SUBROUTINE vnumer_bends ( molecule, x, energy )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( OUT ) :: energy
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: x
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  
! Locals
  REAL ( dbl ), DIMENSION (3) :: b12, b32, g1, g2, g3
  REAL ( dbl ) :: d12, d32, id12, id32, dist, theta, e0, fscalar
  INTEGER :: imol, ibend
  TYPE ( linklist_bends ), POINTER :: llbend
  
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
  
END SUBROUTINE vnumer_bends

!******************************************************************************

SUBROUTINE getv_ptens ( ewald_param, pnode, x, vnb, box, potparm )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: x
  REAL ( dbl ), INTENT ( INOUT ) :: vnb
  TYPE (cell_type), INTENT ( INOUT ) :: box
  TYPE (particle_node_type ), INTENT ( IN ), DIMENSION ( : ) :: pnode
  TYPE (potentialparm_type), INTENT ( IN ), DIMENSION ( :, : ) :: potparm
  TYPE (ewald_parameters_type), INTENT (IN) :: ewald_param
  
! Locals
  INTEGER :: i, j, ii, jj, iatomtype, jatomtype, iexclude, natoms, id
  INTEGER :: zz,yy,xx
  INTEGER :: icell, jcell, kcell, n, itype, jtype, isos
  REAL ( dbl ), DIMENSION (3) :: vec,ri,rij,cell,s,perd,quotient
  REAL ( dbl ) :: energy, fscalar, rijsq, qi, qj, rcut
  TYPE (linklist_exclusion), POINTER :: ll_exclude
  LOGICAL :: match
  INTEGER, DIMENSION ( :, :, : ), ALLOCATABLE :: ncell
  
!------------------------------------------------------------------------------
  
  n = SIZE ( potparm ( :, 1 ) )
  
  isos = 0
  IF ( .NOT. ALLOCATED ( ncell ) ) ALLOCATE ( ncell ( n, n, 3 ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_nonbond_numer', 'ncell', n * n * 3 )
  
  CALL get_cell_param ( box, cell )
  DO itype = 1, n
     DO jtype = 1, n
        rcut = SQRT ( potparm ( itype, jtype ) % rcutsq )
        quotient ( : ) = rcut / cell ( : )
        DO id = 1, 3
           IF ( quotient ( id ) <= 0.5_dbl ) THEN
              ncell ( itype, jtype, id ) = 0
           ELSE
              ncell ( itype, jtype, id ) = CEILING ( rcut/ cell ( id ) )
           END IF
        END DO
     END DO
  END DO
  
!
! initializing the potential energy
!
  vnb = 0.0_dbl
  
!
! starting the loop
!
  CALL get_hinv ( box )
  CALL get_cell_param ( box, cell )
  perd = box % perd
  natoms = size(pnode)
  ILOOP: DO i = 1, natoms 
     iatomtype = pnode(i) % p % prop % ptype
     qi = pnode(i) % p % prop % charge
     ri = x(:,i)
     JLOOP: DO j = 1, natoms
        IF (j>i) THEN
           jatomtype = pnode(j) % p % prop % ptype
           qj = pnode(j) % p % prop % charge
           rij = x(:,j) - ri
           
!
! imaging,transform back to unscaled distances and get potential
!
           IF ( MAXVAL ( ncell ( iatomtype, jatomtype, : ) ) == 0 ) THEN
              rij = pbc ( rij, box )
           END IF
           
           rijsq = DOT_PRODUCT ( rij, rij )
           
           IF ( rijsq <= potparm ( iatomtype, jatomtype ) % rcutsq ) THEN
              
! pointing to the start of the exclusion list of particle i
              ll_exclude => pnode(i) % ex
              match = .FALSE.
              
              EXCL: DO iexclude = 1, pnode(i) % nexcl
                 IF (ll_exclude % p % iatom==j) THEN
                    match = .TRUE.
                    EXIT EXCL
                 END IF
                 ll_exclude => ll_exclude % next
              END DO EXCL
              
              NOTMATCH: IF ( .NOT. match) THEN
                 IF (qi==0.AND.qj==0) THEN
                   CALL potential_f ( rijsq, potparm, qi, qj, &
                      iatomtype, jatomtype, energy )
                 ELSE
                   CALL potential_f ( rijsq, potparm, qi, qj, &
                      iatomtype, jatomtype, energy, ewald_param )
                 ENDIF
!     
! summing up the potential energy
!
                 vnb = vnb + energy
              END IF NOTMATCH
           END IF
           
! searching images if needed
           zz = ncell(iatomtype,jatomtype,3)
           yy = ncell(iatomtype,jatomtype,2)
           xx = ncell(iatomtype,jatomtype,1)
           s(1) = box % h_inv(1,1)*rij(1) &
                + box % h_inv(1,2)*rij(2) &
                + box % h_inv(1,3)*rij(3)
           s(2) = box % h_inv(2,1)*rij(1) &
                + box % h_inv(2,2)*rij(2) &
                + box % h_inv(2,3)*rij(3)
           s(3) = box % h_inv(3,1)*rij(1) &
                + box % h_inv(3,2)*rij(2) &
                + box % h_inv(3,3)*rij(3)
           DO kcell = -zz,zz
              DO jcell = -yy,yy
                 DO icell = -xx,xx
                    IF (icell==0.AND.jcell==0.AND.kcell==0) CYCLE
                    vec ( : ) = (/icell,jcell,kcell/)
                    CALL find_image(s,perd,vec,box % hmat,rijsq,rij)
                    IF ( rijsq <= potparm ( iatomtype, jatomtype ) % rcutsq ) &
                         THEN
                      IF (qi==0.AND.qj==0) THEN
                        CALL potential_f ( rijsq, potparm, qi, qj, &
                           iatomtype, jatomtype, energy )
                      ELSE
                        CALL potential_f ( rijsq, potparm, qi, qj, &
                           iatomtype, jatomtype, energy, ewald_param )
                      ENDIF
                      vnb = vnb + energy
                    END IF
                 END DO
              END DO
           END DO
           
        ELSE IF ( i == j ) THEN
! search images if needed
           zz = ncell(iatomtype,iatomtype,3)
           yy = ncell(iatomtype,iatomtype,2)
           xx = ncell(iatomtype,iatomtype,1)
           s = 0.0_dbl
           DO kcell = -zz,zz
              DO jcell = -yy,yy
                 DO icell = -xx,xx
                    IF ( icell == 0 .AND. jcell == 0 .AND. kcell == 0 ) CYCLE
                    
                    vec ( : ) = (/ icell, jcell, kcell /)
                    CALL find_image ( s, perd, vec, box % hmat, rijsq, rij )
                    IF ( rijsq <= potparm ( iatomtype, iatomtype ) % rcutsq ) &
                         THEN
                      IF (qi==0.AND.qj==0) THEN
                        CALL potential_f ( rijsq, potparm, qi, qj, &
                           iatomtype, jatomtype, energy )
                      ELSE
                        CALL potential_f ( rijsq, potparm, qi, qj, &
                           iatomtype, jatomtype, energy, ewald_param )
                      ENDIF
                      vnb = vnb + energy
                    END IF
                    
                 END DO
              END DO
           END DO
        END IF
        
     END DO JLOOP
  END DO ILOOP
  
END SUBROUTINE getv_ptens

!******************************************************************************

END MODULE fist_force_numer

!******************************************************************************
