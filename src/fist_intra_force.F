!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE fist_intra_force
  
  USE kinds, ONLY : dbl
  USE mol_force, ONLY : force_bonds, force_bends, force_torsions, &
       get_pv_bond, get_pv_bend, get_pv_torsion
  USE molecule_types, ONLY : molecule_structure_type, linklist_bonds, &
       linklist_bends
  USE timings, ONLY : timeset, timestop
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: force_intra_control
  
CONTAINS

!******************************************************************************

SUBROUTINE force_intra_control ( molecule, v_bond, v_bend, pvbond, pvbend, &
          f_bond, f_bend )
  
! parses up intramolecular properties
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  REAL ( dbl ), INTENT ( INOUT ) :: v_bond, v_bend
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: pvbond, pvbend
  REAL ( dbl ), DIMENSION ( :, : ), OPTIONAL, INTENT ( OUT ) :: f_bond, f_bend
  
! Locals
  REAL ( dbl ) :: d12, d32, id12, id32, dist, theta
  REAL ( dbl ) :: energy, fscalar
  REAL ( dbl ), DIMENSION (3) :: rij, b12, b32, g1, g2, g3
  TYPE (linklist_bonds), POINTER :: llbond
  TYPE (linklist_bends), POINTER :: llbend
  INTEGER :: ibond, ibend, imol, handle
  
!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE_INTRA_CONTROL','I',' ',handle )
  
  IF ( PRESENT ( f_bond)) f_bond = 0._dbl
  IF ( PRESENT ( f_bend)) f_bend = 0._dbl
  v_bond =  0._dbl
  v_bend =  0._dbl
  
  MOL: DO imol = 1, size(molecule)
     llbond => molecule(imol) %ll_bonds
     llbend => molecule(imol) %ll_bends
     
     BOND: DO ibond = 1, molecule(imol) %nbonds_mol
        rij = llbond%p1%r - llbond%p2%r
        CALL force_bonds(rij,llbond%bond_param%r0,llbond%bond_param%k, &
             energy,fscalar)
        v_bond = v_bond + energy
        llbond%p1%f = llbond%p1%f - rij*fscalar
        llbond%p2%f = llbond%p2%f + rij*fscalar
        
! computing the pressure tensor
        CALL get_pv_bond ( -rij * fscalar, rij * fscalar, &
             llbond%p1%r, llbond%p2%r, pvbond )
        
! the contribution from the bonds. ONLY FOR DEBUG
        IF ( PRESENT ( f_bond)) THEN
           f_bond(:,llbond%index(1)) = f_bond(:,llbond%index(1)) - &
                rij*fscalar
           f_bond(:,llbond%index(2)) = f_bond(:,llbond%index(2)) + &
                rij*fscalar
        END IF
        
        llbond => llbond%next
     END DO BOND
     
     BEND: DO ibend = 1, molecule(imol) %nbends_mol
        b12 = llbend%p1%r - llbend%p2%r
        b32 = llbend%p3%r - llbend%p2%r
        d12 = sqrt(dot_product(b12,b12))
        id12 = 1._dbl/d12
        d32 = sqrt(dot_product(b32,b32))
        id32 = 1._dbl/d32
        dist = dot_product(b12,b32)
        theta = acos(dist*id12*id32)
        CALL force_bends(b12,b32,d12,d32,id12,id32,dist,theta, &
             llbend%bend_param%theta0,llbend%bend_param%k,g1,g2,g3,energy, &
             fscalar)
        v_bend = v_bend + energy
        llbend%p1%f = llbend%p1%f + fscalar*g1
        llbend%p2%f = llbend%p2%f + fscalar*g2
        llbend%p3%f = llbend%p3%f + fscalar*g3
        
! computing the pressure tensor
        CALL get_pv_bend(fscalar*g1,fscalar*g2,fscalar*g3,llbend%p1%r, &
             llbend%p2%r,llbend%p3%r,pvbend)
        
! the contribution from the bends. ONLY FOR DEBUG
        IF ( PRESENT ( f_bend)) THEN
           f_bend(:,llbend%index(1)) = f_bend(:,llbend%index(1)) + &
                fscalar*g1
           f_bend(:,llbend%index(2)) = f_bend(:,llbend%index(2)) + &
                fscalar*g2
           f_bend(:,llbend%index(3)) = f_bend(:,llbend%index(3)) + &
                fscalar*g3
        END IF
        
        llbend => llbend%next
     END DO BEND
     
  END DO MOL

  CALL timestop ( 0._dbl, handle )
  
END SUBROUTINE force_intra_control

!******************************************************************************

END MODULE fist_intra_force
