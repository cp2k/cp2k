!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/constraint [1.0] *
!!
!!   NAME
!!     constraint
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE constraint
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE molecule_types, ONLY : molecule_structure_type, &
       linklist_dist_constraints, linklist_g3x3_constraints
  USE message_passing, ONLY : mp_sum
  USE stop_program, ONLY : stop_prg
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : getinv3x3, dotprod_3d, matvec_3x3, matmul_3x3, transpose_3d
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: shake_control, rattle_control, pv_constraint, getold
  PUBLIC :: roll, initialize_roll, shake_roll_control, rattle_roll_control
  
  REAL ( dbl ) :: scale
  REAL ( dbl ), DIMENSION (3,3) :: r_rattle

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************

SUBROUTINE shake_control(molecule,simpar)
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  
! Locals
  REAL ( dbl ) :: max_sigma
  INTEGER :: imol, iconst, ishake, nmol, handle
  
!------------------------------------------------------------------------------
  
! initialize locals
  CALL timeset ( 'SHAKE_CONTROL', 'I', ' ', handle )
  
  nmol = SIZE ( molecule )
  MOL: DO imol = 1, nmol
     
! 3x3
     IF (molecule(imol)%ng3x3const_mol /= 0 ) THEN
        CALL shake_3x3(molecule(imol),simpar,ishake)
     END IF
     
     ishake = 0
     
! dist
     IF (molecule(imol)%ndistconst_mol /= 0 ) THEN
        max_sigma = -1.0E+10_dbl
        CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
           max_sigma = 0._dbl
           ishake = ishake + 1
           CALL shake(molecule(imol),simpar,max_sigma,ishake)
        END DO CONS
     END IF
  END DO MOL
  
  CALL timestop ( 0.0_dbl, handle )
  
END SUBROUTINE shake_control

!******************************************************************************

SUBROUTINE rattle_control ( molecule, simpar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  
! Locals
  REAL ( dbl ) :: max_sigma
  INTEGER :: imol, iconst, nmol, irattle, handle
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'RATTLE_CONTROL', 'I', ' ', handle )
! initialize locals
  nmol = SIZE ( molecule )
  
  MOL: DO imol = 1, nmol
! 3x3
     IF ( molecule ( imol ) % ng3x3const_mol /= 0 ) THEN
        CALL rattle_3x3 ( molecule ( imol ), simpar )
     END IF
     irattle = 0
     
! dist
     IF ( molecule ( imol ) % ndistconst_mol /= 0 ) THEN
        max_sigma = -1.0E+10_dbl
        CONS: DO WHILE ( ABS ( max_sigma ) >= simpar % shake_tol )
           max_sigma = 0.0_dbl
           irattle = irattle + 1
           CALL rattle ( molecule ( imol ), simpar, max_sigma, irattle )
        END DO CONS
     END IF
  END DO MOL
  
  CALL timestop ( 0.0_dbl, handle )
  
END SUBROUTINE rattle_control

!******************************************************************************

SUBROUTINE shake_roll_control ( molecule, simpar, thermo, roll_tol, iroll, &
     vector, group, u )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( INOUT ), OPTIONAL :: thermo
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: vector
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  INTEGER, INTENT ( INOUT ) :: iroll
  INTEGER, INTENT ( IN ) :: group
  REAL ( dbl ), INTENT ( IN ), OPTIONAL, DIMENSION ( :, : ) :: u
  
! Locals
  REAL ( dbl ) :: max_sigma
  INTEGER :: imol, iconst, ishake, nmol, handle
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: r_shake
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'SHAKE_ROLL', 'I', ' ', handle )
  
! setting up for roll
  IF ( simpar % ensemble == 'NPT_I' ) THEN
     CALL get_r_shake ( r_shake, vector )
  ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
     CALL get_r_shake ( r_shake, vector, u )
  END IF
  
! initialize locals
  nmol = SIZE ( molecule )
  MOL: DO imol = 1, nmol
     
! 3x3
     IF (molecule(imol)%ng3x3const_mol /= 0 ) THEN
        CALL shake_roll_3x3(molecule(imol),simpar,r_shake,ishake)
     END IF
     ishake = 0
     
! dist
     IF (molecule(imol)%ndistconst_mol /= 0 ) THEN
        max_sigma = -1.0E+10_dbl
        CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
           ishake = ishake + 1
           max_sigma = 0._dbl
           CALL shake_roll(molecule(imol),simpar,r_shake,max_sigma,ishake)
        END DO CONS
     END IF
  END DO MOL
  
  IF ( simpar % ensemble == 'NPT_I' ) &
       CALL check_tol ( roll_tol, iroll, 'SHAKE', vector )
  CALL pv_constraint ( molecule, thermo, 'FULL', group )
  
  CALL timestop ( 0.0_dbl, handle )
  
END SUBROUTINE shake_roll_control

!******************************************************************************

SUBROUTINE rattle_roll_control ( molecule, simpar, thermo, veps, &
     roll_tol, iroll, group )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: veps
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  INTEGER, INTENT ( INOUT ) :: iroll
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  REAL ( dbl ) :: max_sigma, feps
  INTEGER :: imol, iconst, irattle, nmol, handle
  
!------------------------------------------------------------------------------
  
  CALL timeset( 'RATTLE_ROLL','I',' ',handle)
  
! initialize locals
  irattle = 0
  nmol = SIZE ( molecule )
  max_sigma = -1.0E+10_dbl
  
! Initialize scale and r_rattle
  CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
     irattle = irattle + 1
     max_sigma = 0._dbl
     MOL: DO imol = 1, nmol
        IF (molecule(imol)%ndistconst_mol /= 0 ) THEN
           
! dist
           CALL rattle_roll(molecule(imol),simpar,veps,max_sigma,irattle)
        ELSE IF (molecule(imol)%ng3x3const_mol /= 0 ) THEN
           
! 3x3
           CALL rattle_roll_3x3(molecule(imol),simpar,veps,max_sigma,irattle)
        END IF
     END DO MOL
     
! update veps
     IF ( irattle /= 1 ) CALL update_veps ( molecule, &
          veps, thermo, simpar % ensemble, group )
  END DO CONS
  
  CALL check_tol ( roll_tol, iroll, 'RATTLE', molecule = molecule )
  CALL pv_constraint ( molecule, thermo, 'FULL', group )
  
  CALL timestop ( 0.0_dbl, handle )
  
END SUBROUTINE rattle_roll_control

!******************************************************************************

SUBROUTINE shake_3x3 ( molecule, simpar, ishake )
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  INTEGER, INTENT ( OUT ) :: ishake
  
! Locals
  REAL ( dbl ) :: idtsq, dtsqby2, dtby2, imass1, imass2, imass3
  REAL ( dbl ) :: sigma, max_sigma
  REAL ( dbl ), DIMENSION(3) :: r12, r23, r13, vec, bvec
  REAL ( dbl ), DIMENSION(3) :: r0_12, r0_23, r0_13
  REAL ( dbl ), DIMENSION(3) :: r1, r2, r3, v1, v2, v3, fc1, fc2, fc3
  REAL ( dbl ), DIMENSION(3,3) :: amat, ainv
  TYPE (linklist_g3x3_constraints), POINTER :: llconst
  INTEGER :: info, iconst
  
!------------------------------------------------------------------------------
  
! define local variables
  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
  idtsq = 1._dbl/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  llconst => molecule%ll_g3x3_const
  DO iconst = 1, molecule%ng3x3const_mol
     ishake=0
     r0_12=llconst%p1%r-llconst%p2%r
     r0_13=llconst%p1%r-llconst%p3%r
     r0_23=llconst%p2%r-llconst%p3%r
     imass1 = 1._dbl/llconst%p1%prop%mass
     imass2 = 1._dbl/llconst%p2%prop%mass
     imass3 = 1._dbl/llconst%p3%prop%mass
     llconst%fc1 = -2._dbl*(llconst%rold1-llconst%rold2)
     llconst%fc2 = -2._dbl*(llconst%rold1-llconst%rold3)
     llconst%fc3 = -2._dbl*(llconst%rold2-llconst%rold3)
     
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,llconst%fc1)
     amat(1,2)=imass1*DOTPROD_3D(r0_12,llconst%fc2)
     amat(1,3)=-imass2*DOTPROD_3D(r0_12,llconst%fc3)
     amat(2,1)=imass1*DOTPROD_3D(r0_13,llconst%fc1)
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,llconst%fc2)
     amat(2,3)=imass3*DOTPROD_3D(r0_13,llconst%fc3)
     amat(3,1)=-imass2*DOTPROD_3D(r0_23,llconst%fc1)
     amat(3,2)=imass3*DOTPROD_3D(r0_23,llconst%fc2)
     amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,llconst%fc3)
     
! Iterate until convergence:
     max_sigma = -1.0E+10_dbl
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0._dbl
        ishake=ishake+1
        vec= llconst%lambda(1)*llconst%fc1*(imass1+imass2) + &
             llconst%lambda(2)*imass1*llconst%fc2 - &
             llconst%lambda(3)*imass2*llconst%fc3
        bvec(1)=llconst%d12*llconst%d12-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
             -DOTPROD_3D(r0_12,r0_12)
        vec= llconst%lambda(1)*llconst%fc1*imass1 + &
             llconst%lambda(2)*(imass1+imass3)*llconst%fc2 + &
             llconst%lambda(3)*imass3*llconst%fc3
        bvec(2)=llconst%d13*llconst%d13-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
             -DOTPROD_3D(r0_13,r0_13)
        vec= -llconst%lambda(1)*llconst%fc1*imass2 + &
             llconst%lambda(2)*imass3*llconst%fc2 + &
             llconst%lambda(3)*(imass2+imass3)*llconst%fc3
        bvec(3)=llconst%d23*llconst%d23-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
             -DOTPROD_3D(r0_23,r0_23)
        bvec = bvec * idtsq
        
! find inverse of amat
        CALL getinv3x3 ( amat, ainv )
        
! get lambda
        llconst%lambda ( : )=MATVEC_3x3(ainv,bvec)
        fc1=llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(2)*llconst%fc2
        fc2=-llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(3)*llconst%fc3
        fc3=-llconst%lambda(2)*llconst%fc2- &
             llconst%lambda(3)*llconst%fc3
        r1 = llconst%p1%r + imass1*dtsqby2*fc1
        r2 = llconst%p2%r + imass2*dtsqby2*fc2
        r3 = llconst%p3%r + imass3*dtsqby2*fc3
        v1 = llconst%p1%v + imass1*dtby2*fc1
        v2 = llconst%p2%v + imass2*dtby2*fc2
        v3 = llconst%p3%v + imass3*dtby2*fc3
        r12=r1-r2
        r13=r1-r3
        r23=r2-r3
        
! compute the tolerance:
        sigma = dot_product(r12,r12) - llconst%d12*llconst%d12 
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r13,r13) - llconst%d13*llconst%d13 
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r23,r23) - llconst%d23*llconst%d23 
        max_sigma=MAX(max_sigma,ABS(sigma))
     END DO CONS
     
! update positions with full multiplier
     llconst%p1%r = r1
     llconst%p2%r = r2
     llconst%p3%r = r3
     
! update velocites with full multiplier
     llconst%p1%v = v1
     llconst%p2%v = v2
     llconst%p3%v = v3
     llconst => llconst%next
     
  END DO
  
END SUBROUTINE shake_3x3

!******************************************************************************

SUBROUTINE shake_roll_3x3(molecule,simpar,r_shake,ishake)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: r_shake
  INTEGER, INTENT ( OUT ) :: ishake
  
! Locals
  REAL ( dbl ) :: idtsq, dtsqby2, dtby2, imass1, imass2, imass3
  REAL ( dbl ) :: sigma, max_sigma
  REAL ( dbl ), DIMENSION(3) :: r12, r23, r13, vec, bvec
  REAL ( dbl ), DIMENSION(3) :: r0_12, r0_23, r0_13
  REAL ( dbl ), DIMENSION(3) :: r1, r2, r3, v1, v2, v3, fc1, fc2, fc3, f_roll1, &
       f_roll2, f_roll3
  REAL ( dbl ), DIMENSION(3,3) :: amat, ainv
  TYPE (linklist_g3x3_constraints), POINTER :: llconst
  INTEGER :: info, iconst
  
!------------------------------------------------------------------------------
  
! define local variables
  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
  idtsq = 1._dbl/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  llconst => molecule%ll_g3x3_const
  
  DO iconst = 1, molecule%ng3x3const_mol
     ishake=0
     r0_12=llconst%p1%r-llconst%p2%r
     r0_13=llconst%p1%r-llconst%p3%r
     r0_23=llconst%p2%r-llconst%p3%r
     imass1 = 1._dbl/llconst%p1%prop%mass
     imass2 = 1._dbl/llconst%p2%prop%mass
     imass3 = 1._dbl/llconst%p3%prop%mass
     llconst%fc1 = -2._dbl*(llconst%rold1-llconst%rold2)
     llconst%fc2 = -2._dbl*(llconst%rold1-llconst%rold3)
     llconst%fc3 = -2._dbl*(llconst%rold2-llconst%rold3)
     
! rotate fconst:  
     f_roll1 = matvec_3x3(r_shake,llconst%fc1)
     f_roll2 = matvec_3x3(r_shake,llconst%fc2)
     f_roll3 = matvec_3x3(r_shake,llconst%fc3)
     
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,f_roll1)
     amat(1,2)=imass1*DOTPROD_3D(r0_12,f_roll2)
     amat(1,3)=-imass2*DOTPROD_3D(r0_12,f_roll3)
     amat(2,1)=imass1*DOTPROD_3D(r0_13,f_roll1)
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,f_roll2)
     amat(2,3)=imass3*DOTPROD_3D(r0_13,f_roll3)
     amat(3,1)=-imass2*DOTPROD_3D(r0_23,f_roll1)
     amat(3,2)=imass3*DOTPROD_3D(r0_23,f_roll2)
     amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,f_roll3)
     
! Iterate until convergence:
     max_sigma = -1.0E+10_dbl
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0._dbl
        ishake=ishake+1
        vec= llconst%lambda(1)*f_roll1*(imass1+imass2) + &
             llconst%lambda(2)*imass1*f_roll2 - &
             llconst%lambda(3)*imass2*f_roll3
        bvec(1)=llconst%d12*llconst%d12-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
             -DOTPROD_3D(r0_12,r0_12)
        vec= llconst%lambda(1)*f_roll1*imass1 + &
             llconst%lambda(2)*(imass1+imass3)*f_roll2 + &
             llconst%lambda(3)*imass3*f_roll3
        bvec(2)=llconst%d13*llconst%d13-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
             -DOTPROD_3D(r0_13,r0_13)
        vec= -llconst%lambda(1)*f_roll1*imass2 + &
             llconst%lambda(2)*imass3*f_roll2 + &
             llconst%lambda(3)*(imass2+imass3)*f_roll3
        bvec(3)=llconst%d23*llconst%d23-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
             -DOTPROD_3D(r0_23,r0_23)
        bvec=bvec*idtsq
        
! find inverse of amat
        CALL getinv3x3 ( amat, ainv )
        
! get lambda
        llconst%lambda ( : )=MATVEC_3x3(ainv,bvec)
        fc1=llconst%lambda(1)*f_roll1+ &
             llconst%lambda(2)*f_roll2
        fc2=-llconst%lambda(1)*f_roll1+ &
             llconst%lambda(3)*f_roll3
        fc3=-llconst%lambda(2)*f_roll2- &
             llconst%lambda(3)*f_roll3
        r1 = llconst%p1%r + imass1*dtsqby2*fc1
        r2 = llconst%p2%r + imass2*dtsqby2*fc2
        r3 = llconst%p3%r + imass3*dtsqby2*fc3
        v1 = llconst%p1%v + imass1*dtby2*fc1
        v2 = llconst%p2%v + imass2*dtby2*fc2
        v3 = llconst%p3%v + imass3*dtby2*fc3
        r12=r1-r2
        r13=r1-r3
        r23=r2-r3
        
! compute the tolerance:
        sigma = dot_product(r12,r12) - llconst%d12*llconst%d12 
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r13,r13) - llconst%d13*llconst%d13 
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r23,r23) - llconst%d23*llconst%d23 
        max_sigma=MAX(max_sigma,ABS(sigma))
     END DO CONS
     
! update positions with full multiplier
     llconst%p1%r = r1
     llconst%p2%r = r2
     llconst%p3%r = r3
     
! update velocites with full multiplier
     llconst%p1%v = v1
     llconst%p2%v = v2
     llconst%p3%v = v3
     llconst => llconst%next
     
  END DO
  
END SUBROUTINE shake_roll_3x3

!******************************************************************************

SUBROUTINE rattle_roll_3x3(molecule,simpar,veps,max_sigma,irattle)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  INTEGER, INTENT ( IN ) :: irattle
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: veps
  
! Locals
  INTEGER :: iconst, info
  REAL ( dbl ) :: idt, dtby2, imass1, imass2, imass3, sigma_dot
  REAL ( dbl ), DIMENSION(3) :: r12, r23, r13, vec, bvec
  REAL ( dbl ), DIMENSION(3) :: f_roll1, f_roll2, f_roll3
  REAL ( dbl ), DIMENSION(3) :: v12_0, v23_0, v13_0
  REAL ( dbl ), DIMENSION(3) :: v12, v23, v13
  REAL ( dbl ), DIMENSION(3) :: v1, v2, v3
  REAL ( dbl ), DIMENSION(3) :: v1_0, v2_0, v3_0
  REAL ( dbl ), DIMENSION(3) :: fc1, fc2, fc3, lambda
  REAL ( dbl ), DIMENSION(3,3) :: amat, ainv
  TYPE (linklist_g3x3_constraints), POINTER :: llconst
  
!------------------------------------------------------------------------------
  
! define local variables
  idt = 1._dbl/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  IF (irattle==1) THEN
     llconst => molecule%ll_g3x3_const
     DO iconst = 1, molecule%ng3x3const_mol
        llconst%vold1=llconst%p1%v
        llconst%vold2=llconst%p2%v
        llconst%vold3=llconst%p3%v
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        imass3 = 1._dbl/llconst%p3%prop%mass
        fc1=llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(2)*llconst%fc2
        fc2=-llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(3)*llconst%fc3
        fc3=-llconst%lambda(2)*llconst%fc2- &
             llconst%lambda(3)*llconst%fc3
        llconst%p1%v=llconst%p1%v+imass1*dtby2*fc1
        llconst%p2%v=llconst%p2%v+imass2*dtby2*fc2
        llconst%p3%v=llconst%p3%v+imass3*dtby2*fc3
        llconst => llconst%next
     END DO
     max_sigma = -1.0E+10_dbl
     
  ELSE
     
     llconst => molecule%ll_g3x3_const
     DO iconst = 1, molecule%ng3x3const_mol
        v12_0=llconst%vold1-llconst%vold2
        v13_0=llconst%vold1-llconst%vold3
        v23_0=llconst%vold2-llconst%vold3
        r12=llconst%p1%r-llconst%p2%r
        r13=llconst%p1%r-llconst%p3%r
        r23=llconst%p2%r-llconst%p3%r
        
! roll the local velocities
        v12_0=matvec_3x3(r_rattle,v12_0)*llconst%scale
        v13_0=matvec_3x3(r_rattle,v13_0)*llconst%scale
        v23_0=matvec_3x3(r_rattle,v23_0)*llconst%scale
        v1_0=matvec_3x3(r_rattle,llconst%vold1)*llconst%scale
        v2_0=matvec_3x3(r_rattle,llconst%vold2)*llconst%scale
        v3_0=matvec_3x3(r_rattle,llconst%vold3)*llconst%scale
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        imass3 = 1._dbl/llconst%p3%prop%mass
        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
        llconst%fc3 = -2._dbl*(llconst%p2%r-llconst%p3%r)
        
! roll the fc
        f_roll1=matvec_3x3(r_rattle,llconst%fc1)*llconst%scale
        f_roll2=matvec_3x3(r_rattle,llconst%fc2)*llconst%scale
        f_roll3=matvec_3x3(r_rattle,llconst%fc3)*llconst%scale
        
! construct matrix
        amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,f_roll1)
        amat(1,2)=imass1*DOTPROD_3D(r12,f_roll2)
        amat(1,3)=-imass2*DOTPROD_3D(r12,f_roll3)
        amat(2,1)=imass1*DOTPROD_3D(r13,f_roll1)
        amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,f_roll2)
        amat(2,3)=imass3*DOTPROD_3D(r13,f_roll3)
        amat(3,1)=-imass2*DOTPROD_3D(r23,f_roll1)
        amat(3,2)=imass3*DOTPROD_3D(r23,f_roll2)
        amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23,f_roll3)
        
! construct solution vector
        bvec(1)=DOTPROD_3D(r12,v12_0+MATVEC_3X3(veps,r12))
        bvec(2)=DOTPROD_3D(r13,v13_0+MATVEC_3X3(veps,r13))
        bvec(3)=DOTPROD_3D(r23,v23_0+MATVEC_3X3(veps,r23))
        bvec = -bvec * 2.0_dbl * idt
        
! find inverse of amat
        CALL getinv3x3 ( amat, ainv )
        
! get lambda
        lambda ( : )=MATVEC_3X3(ainv,bvec)
        llconst%del_lam=lambda-llconst%lambda
        
! update local velocities for calculation of sigma_dot
! NOTE: use full multiplier because always starting
! from reference velocity, v1_0, v2_0, v3_0.
        fc1=llconst%lambda(1)*f_roll1+ &
             llconst%lambda(2)*f_roll2
        fc2=-llconst%lambda(1)*f_roll1+ &
             llconst%lambda(3)*f_roll3
        fc3=-llconst%lambda(2)*f_roll2- &
             llconst%lambda(3)*f_roll3
        v1=v1_0+imass1*dtby2*fc1
        v2=v2_0+imass2*dtby2*fc2
        v3=v3_0+imass3*dtby2*fc3
        
! update true velocities
        fc1=llconst%del_lam(1)*llconst%fc1+ &
             llconst%del_lam(2)*llconst%fc2
        fc2=-llconst%del_lam(1)*llconst%fc1+ &
             llconst%del_lam(3)*llconst%fc3
        fc3=-llconst%del_lam(2)*llconst%fc2- &
             llconst%del_lam(3)*llconst%fc3
        llconst%p1%v=llconst%p1%v+imass1*dtby2*fc1
        llconst%p2%v=llconst%p2%v+imass2*dtby2*fc2
        llconst%p3%v=llconst%p3%v+imass3*dtby2*fc3
        
! computing the convergence of the sigma_dot
        v12=v1-v2
        sigma_dot=2._dbl*DOTPROD_3D(r12,v12+MATVEC_3X3(veps,r12))
        max_sigma = MAX(abs(sigma_dot),max_sigma)
        v13=v1-v3
        sigma_dot=2._dbl*DOTPROD_3D(r13,v13+MATVEC_3X3(veps,r13))
        max_sigma = MAX(abs(sigma_dot),max_sigma)
        v23=v2-v3
        sigma_dot=2._dbl*DOTPROD_3D(r23,v23+MATVEC_3X3(veps,r23))
        max_sigma = MAX(abs(sigma_dot),max_sigma)
        llconst%lambda=lambda
        llconst => llconst%next
     END DO
     
  END IF
  
END SUBROUTINE rattle_roll_3x3

!******************************************************************************

SUBROUTINE rattle_3x3(molecule,simpar)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  
! Locals
  INTEGER :: iconst, info
  REAL ( dbl ) :: idt, dtby2, imass1, imass2, imass3
  REAL ( dbl ), DIMENSION(3) :: r12, r23, r13, vec, bvec
  REAL ( dbl ), DIMENSION(3) :: v12, v23, v13
  REAL ( dbl ), DIMENSION(3) :: fc1, fc2, fc3
  REAL ( dbl ), DIMENSION(3,3) :: amat, ainv
  TYPE (linklist_g3x3_constraints), POINTER :: llconst
  
!------------------------------------------------------------------------------
  
! define local variables
  idt = 1._dbl/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  llconst => molecule%ll_g3x3_const
  DO iconst = 1, molecule%ng3x3const_mol
     v12=llconst%p1%v-llconst%p2%v
     v13=llconst%p1%v-llconst%p3%v
     v23=llconst%p2%v-llconst%p3%v
     r12=llconst%p1%r-llconst%p2%r
     r13=llconst%p1%r-llconst%p3%r
     r23=llconst%p2%r-llconst%p3%r
     imass1 = 1._dbl/llconst%p1%prop%mass
     imass2 = 1._dbl/llconst%p2%prop%mass
     imass3 = 1._dbl/llconst%p3%prop%mass
     llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
     llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
     llconst%fc3 = -2._dbl*(llconst%p2%r-llconst%p3%r)
     
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,llconst%fc1)
     amat(1,2)=imass1*DOTPROD_3D(r12,llconst%fc2)
     amat(1,3)=-imass2*DOTPROD_3D(r12,llconst%fc3)
     amat(2,1)=imass1*DOTPROD_3D(r13,llconst%fc1)
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,llconst%fc2)
     amat(2,3)=imass3*DOTPROD_3D(r13,llconst%fc3)
     amat(3,1)=-imass2*DOTPROD_3D(r23,llconst%fc1)
     amat(3,2)=imass3*DOTPROD_3D(r23,llconst%fc2)
     amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23,llconst%fc3)
     
! construct solution vector
     bvec(1)=DOTPROD_3D(r12,v12)
     bvec(2)=DOTPROD_3D(r13,v13)
     bvec(3)=DOTPROD_3D(r23,v23)
     bvec=-bvec*2._dbl*idt
     
! find inverse of amat
     CALL getinv3x3(amat,ainv)
     
! get lambda
     llconst%lambda ( : )=MATVEC_3X3(ainv,bvec)
     fc1=llconst%lambda(1)*llconst%fc1+ &
          llconst%lambda(2)*llconst%fc2
     fc2=-llconst%lambda(1)*llconst%fc1+ &
          llconst%lambda(3)*llconst%fc3
     fc3=-llconst%lambda(2)*llconst%fc2- &
          llconst%lambda(3)*llconst%fc3
     llconst%p1%v=llconst%p1%v+imass1*dtby2*fc1
     llconst%p2%v=llconst%p2%v+imass2*dtby2*fc2
     llconst%p3%v=llconst%p3%v+imass3*dtby2*fc3
     llconst => llconst%next
  END DO
  
END SUBROUTINE rattle_3x3

!******************************************************************************

SUBROUTINE shake(molecule,simpar,max_sigma,ishake)
  
! updates the multiplier one molecule type at a time
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  INTEGER, INTENT ( IN ) :: ishake
  REAL ( dbl ) :: dtsqby2, dtby2, imass1, imass2, fdotf, fdotf_sum, del_lam
  REAL ( dbl ), DIMENSION (3) :: fc_t, r12
  TYPE (linklist_dist_constraints), POINTER :: llconst
  INTEGER :: imol, iconst
  
!------------------------------------------------------------------------------
  
! define local variables
  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
  dtby2 = simpar%dt*.5_dbl
  IF (ishake==1) THEN
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        llconst%fc = -2._dbl*(llconst%rold1-llconst%rold2)
        
! update positions with full multiplier
        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%lambda* &
             llconst%fc
        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%lambda* &
             llconst%fc
        
! update velocites with full multiplier
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%lambda*llconst% &
             fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%lambda*llconst% &
             fc
        llconst => llconst%next
     END DO
     
  ELSE
     
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        
! computing the constraint
        r12 = llconst%p1%r - llconst%p2%r
        llconst%sigma = dot_product(r12,r12) - llconst%distance*llconst% &
             distance
        
! computing fc(dt)
        fc_t = -2._dbl*(llconst%p1%r-llconst%p2%r)
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        fdotf = dot_product(llconst%fc,fc_t)
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2._dbl*llconst%sigma/(simpar%dt*simpar%dt*fdotf_sum)
        llconst%lambda = llconst%lambda + del_lam
        
! update positions by delta_lam
        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*del_lam*llconst%fc
        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*del_lam*llconst%fc
        
! update velocites by delta_lam
        llconst%p1%v = llconst%p1%v + imass1*dtby2*del_lam*llconst%fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*del_lam*llconst%fc
        llconst => llconst%next
     END DO
  END IF
  
! computing the constraint and value of tolerance 
  llconst => molecule%ll_dist_const
  DO iconst = 1, molecule%ndistconst_mol
     r12 = llconst%p1%r - llconst%p2%r
     llconst%sigma = dot_product(r12,r12) - llconst%distance*llconst% &
          distance
     max_sigma = max(abs(llconst%sigma),max_sigma)
     llconst => llconst%next
  END DO
  
END SUBROUTINE shake

!******************************************************************************

SUBROUTINE rattle ( molecule, simpar, max_sigma, irattle )
  
! updates the multiplier one molecule type at a time
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  INTEGER, INTENT ( IN ) :: irattle
  
! Locals
  REAL ( dbl ) :: dtby2, imass1, imass2, fdotf, fdotf_sum, del_lam
  REAL ( dbl ), DIMENSION (3) :: r12, v12
  TYPE (linklist_dist_constraints), POINTER :: llconst
  INTEGER :: imol, iconst
  
!------------------------------------------------------------------------------
  
! define local variables
  dtby2 = simpar%dt*.5_dbl
  IF (irattle==1) THEN
     llconst => molecule%ll_dist_const
! if statement to catch molecule types with no constraints.:
     DO iconst = 1, molecule%ndistconst_mol
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        llconst%fc = -2._dbl*(llconst%p1%r-llconst%p2%r)
        
! update velocites with full multiplier
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%lambda*llconst% &
             fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%lambda*llconst% &
             fc
        llconst => llconst%next
     END DO
     
  ELSE
     
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        
! computing sigma
        r12 = llconst%p1%r - llconst%p2%r
        v12 = llconst%p1%v - llconst%p2%v
        llconst%sigma = 2._dbl*dot_product(r12,v12)
        
! computing f(dt)
        llconst%fc = -2._dbl*(llconst%p1%r-llconst%p2%r)
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        fdotf = dot_product(llconst%fc,llconst%fc)
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2._dbl*llconst%sigma/(simpar%dt*fdotf_sum)
        llconst%lambda = llconst%lambda + del_lam
        
!update velocites with del_lam
        llconst%p1%v = llconst%p1%v + imass1*dtby2*del_lam*llconst%fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*del_lam*llconst%fc
        llconst => llconst%next
     END DO
  END IF
  
! computing the constraint and value of the tolerance
  llconst => molecule%ll_dist_const
  DO iconst = 1, molecule%ndistconst_mol
     
! computing the updated constraint
     r12 = llconst%p1%r - llconst%p2%r
     v12 = llconst%p1%v - llconst%p2%v
     llconst%sigma = 2._dbl*dot_product(r12,v12)
     max_sigma = max(abs(llconst%sigma),max_sigma)
     llconst => llconst%next
  END DO
  
END SUBROUTINE rattle

!******************************************************************************

SUBROUTINE rattle_roll(molecule,simpar,veps,max_sigma,irattle)

! updates the multiplier one molecule type at a time
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  INTEGER, INTENT ( IN ) :: irattle
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: veps
  
! Locals
  REAL ( dbl ) :: dtby2, imass1, imass2, fdotf, fdotf_sum, del_lam
  REAL ( dbl ), DIMENSION (3) :: r12, v12, f_roll
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: v1, v2
  TYPE (linklist_dist_constraints), POINTER :: llconst
  INTEGER :: imol, iconst, isos
  
!------------------------------------------------------------------------------
  
! define local variables
  dtby2 = simpar%dt*.5_dbl
  IF (irattle==1) THEN
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        llconst%fc = -2._dbl*(llconst%p1%r-llconst%p2%r)
        
! update velocites with full multiplier
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%lambda*llconst% &
             fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%lambda*llconst% &
             fc
        llconst => llconst%next
     END DO
     max_sigma = -1.E+10_dbl
     
  ELSE
     llconst => molecule%ll_dist_const
     
! allocating local vs
     ALLOCATE (v1(3,molecule%ndistconst_mol),STAT=isos)
     IF ( isos /= 0 ) CALL stop_prg( 'rattle_roll','failed to allocate v1')
     ALLOCATE (v2(3,molecule%ndistconst_mol),STAT=isos)
     IF ( isos /= 0 ) CALL stop_prg( 'rattle_roll','failed to allocate v2')
     
     DO iconst = 1, molecule%ndistconst_mol
        
! assign local variable
        v1(:,iconst) = llconst%p1%v
        v2(:,iconst) = llconst%p2%v
        r12 = llconst%p1%r - llconst%p2%r
        
! Roll the local vs and fc
        v1(:,iconst) = matvec_3x3(r_rattle ( :, : ),v1(:,iconst))
        v2(:,iconst) = matvec_3x3(r_rattle ( :, : ),v2(:,iconst))
        v1(:,iconst) = v1(:,iconst)*llconst%scale
        v2(:,iconst) = v2(:,iconst)*llconst%scale
        llconst%fc = -2._dbl*r12
        f_roll = matvec_3x3(r_rattle,llconst%fc)
        f_roll = f_roll*llconst%scale
        
! computing sigma
        v12 ( : ) = v1(:,iconst) + matvec_3x3(veps,llconst%p1%r) &
             - v2(:,iconst) - matvec_3x3(veps,llconst%p2%r)
        llconst%sigma = 2._dbl*dot_product(r12,v12)
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        fdotf = dot_product(f_roll,llconst%fc)
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        llconst%del_lam = 2._dbl*llconst%sigma/(simpar%dt*fdotf_sum)
        llconst%lambda = llconst%lambda + llconst%del_lam
        
!update local velocities with del_lam
        v1(:,iconst) = v1(:,iconst) &
             + imass1*dtby2*llconst%del_lam*f_roll ( : )
        v2(:,iconst) = v2(:,iconst) &
             - imass2*dtby2*llconst%del_lam*f_roll ( : )
        
!update true velocites with del_lam
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%del_lam*llconst &
             %fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%del_lam*llconst &
             %fc
        llconst => llconst%next
     END DO
     
! computing the constraint and value of the tolerance with local vs
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        
! computing the updated constraint
        r12 = llconst%p1%r - llconst%p2%r
        v12 ( : ) = v1(:,iconst)-v2(:,iconst) 
        llconst%sigma = 2._dbl*DOTPROD_3D(r12,v12+MATVEC_3X3(veps,r12))
        max_sigma = max(abs(llconst%sigma),max_sigma)
        llconst => llconst%next
     END DO
     
     DEALLOCATE (v1,STAT=isos)
     IF ( isos /= 0 ) CALL stop_prg( 'rattle_roll','failed to deallocate v1')
     DEALLOCATE (v2,STAT=isos)
     IF ( isos /= 0 ) CALL stop_prg( 'rattle_roll','failed to deallocate v2')
     
  END IF
  
END SUBROUTINE rattle_roll

!******************************************************************************

SUBROUTINE getold(molecule)
  
  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: molecule
  
! Locals
  INTEGER :: imol, iconst
  TYPE (linklist_dist_constraints), POINTER :: ll_dist
  TYPE (linklist_g3x3_constraints), POINTER :: ll_g3x3
  
!------------------------------------------------------------------------------
  
! saves all of the old variables
  DO imol = 1, size(molecule)
     ll_dist => molecule(imol) %ll_dist_const
     DO iconst = 1, molecule(imol) %ndistconst_mol
        ll_dist%rold1(1) = ll_dist%p1%r(1)
        ll_dist%rold1(2) = ll_dist%p1%r(2)
        ll_dist%rold1(3) = ll_dist%p1%r(3)
        ll_dist%rold2(1) = ll_dist%p2%r(1)
        ll_dist%rold2(2) = ll_dist%p2%r(2)
        ll_dist%rold2(3) = ll_dist%p2%r(3)
        ll_dist => ll_dist%next
     END DO
     
     ll_g3x3 => molecule(imol) %ll_g3x3_const
     DO iconst = 1, molecule(imol) %ng3x3const_mol
        ll_g3x3%rold1(1) = ll_g3x3%p1%r(1)
        ll_g3x3%rold1(2) = ll_g3x3%p1%r(2)
        ll_g3x3%rold1(3) = ll_g3x3%p1%r(3)
        ll_g3x3%rold2(1) = ll_g3x3%p2%r(1)
        ll_g3x3%rold2(2) = ll_g3x3%p2%r(2)
        ll_g3x3%rold2(3) = ll_g3x3%p2%r(3)
        ll_g3x3%rold3(1) = ll_g3x3%p3%r(1)
        ll_g3x3%rold3(2) = ll_g3x3%p3%r(2)
        ll_g3x3%rold3(3) = ll_g3x3%p3%r(3)
        ll_g3x3 => ll_g3x3%next
     END DO
     
  END DO
  
END SUBROUTINE getold

!******************************************************************************

SUBROUTINE pv_constraint ( molecule, thermo, char, group )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  REAL ( dbl ) :: lambda, pv ( 3, 3 ), lambda_3x3 ( 3 )
  REAL ( dbl ) :: fc1 ( 3 ), fc2 ( 3 ), fc3 ( 3 )
  TYPE ( linklist_dist_constraints ), POINTER :: llconst
  TYPE ( linklist_g3x3_constraints ), POINTER :: ll3x3
  INTEGER :: imol, iconst
  
!------------------------------------------------------------------------------
  
  pv = 0.0_dbl
  DO imol = 1, size(molecule)
     
! dist
     llconst => molecule(imol) %ll_dist_const
     DO iconst = 1, molecule(imol) %ndistconst_mol
        
        SELECT CASE (char)
        CASE ( 'FULL')
!  pv gets updated with FULL multiplier
           lambda = llconst%lambda
           
        CASE ( 'DEL')
!  pv gets updated with del_lam
           lambda = llconst%del_lam
        END SELECT
        
!pv(1,1)
        pv(1,1) = pv(1,1) + lambda*llconst%fc(1)*llconst%p1%r(1)
        pv(1,1) = pv(1,1) - lambda*llconst%fc(1)*llconst%p2%r(1)
!pv(1,2)
        pv(1,2) = pv(1,2) + lambda*llconst%fc(1)*llconst%p1%r(2)
        pv(1,2) = pv(1,2) - lambda*llconst%fc(1)*llconst%p2%r(2)
!pv(1,3)
        pv(1,3) = pv(1,3) + lambda*llconst%fc(1)*llconst%p1%r(3)
        pv(1,3) = pv(1,3) - lambda*llconst%fc(1)*llconst%p2%r(3)
!pv(2,1) 
        pv(2,1) = pv(2,1) + lambda*llconst%fc(2)*llconst%p1%r(1)
        pv(2,1) = pv(2,1) - lambda*llconst%fc(2)*llconst%p2%r(1)
!pv(2,2)
        pv(2,2) = pv(2,2) + lambda*llconst%fc(2)*llconst%p1%r(2)
        pv(2,2) = pv(2,2) - lambda*llconst%fc(2)*llconst%p2%r(2)
!pv(2,3)
        pv(2,3) = pv(2,3) + lambda*llconst%fc(2)*llconst%p1%r(3)
        pv(2,3) = pv(2,3) - lambda*llconst%fc(2)*llconst%p2%r(3)
!pv(3,1)
        pv(3,1) = pv(3,1) + lambda*llconst%fc(3)*llconst%p1%r(1)
        pv(3,1) = pv(3,1) - lambda*llconst%fc(3)*llconst%p2%r(1)
!pv(3,2)
        pv(3,2) = pv(3,2) + lambda*llconst%fc(3)*llconst%p1%r(2)
        pv(3,2) = pv(3,2) - lambda*llconst%fc(3)*llconst%p2%r(2)
!pv(3,3)
        pv(3,3) = pv(3,3) + lambda*llconst%fc(3)*llconst%p1%r(3)
        pv(3,3) = pv(3,3) - lambda*llconst%fc(3)*llconst%p2%r(3)
        llconst => llconst%next
     END DO
     
! 3x3
     ll3x3 => molecule(imol)%ll_g3x3_const
     DO iconst = 1, molecule(imol)%ng3x3const_mol
        
        SELECT CASE (char)
        CASE ( 'FULL')
!  pv gets updated with FULL multiplier
           lambda_3x3 = ll3x3%lambda
           
        CASE ( 'DEL')
!  pv gets updated with del_lam
           lambda_3x3 = ll3x3%del_lam
        END SELECT
        
        fc1=lambda_3x3(1)*ll3x3%fc1+ &
             lambda_3x3(2)*ll3x3%fc2
        fc2=-lambda_3x3(1)*ll3x3%fc1+ &
             lambda_3x3(3)*ll3x3%fc3
        fc3=-lambda_3x3(2)*ll3x3%fc2- &
             lambda_3x3(3)*ll3x3%fc3
        
!pv(1,1)
        pv(1,1) = pv(1,1) + fc1(1)*ll3x3%p1%r(1)
        pv(1,1) = pv(1,1) + fc2(1)*ll3x3%p2%r(1)
        pv(1,1) = pv(1,1) + fc3(1)*ll3x3%p3%r(1)
!pv(1,2)
        pv(1,2) = pv(1,2) + fc1(1)*ll3x3%p1%r(2)
        pv(1,2) = pv(1,2) + fc2(1)*ll3x3%p2%r(2)
        pv(1,2) = pv(1,2) + fc3(1)*ll3x3%p3%r(2)
!pv(1,3)
        pv(1,3) = pv(1,3) + fc1(1)*ll3x3%p1%r(3)
        pv(1,3) = pv(1,3) + fc2(1)*ll3x3%p2%r(3)
        pv(1,3) = pv(1,3) + fc3(1)*ll3x3%p3%r(3)
!pv(2,1) 
        pv(2,1) = pv(2,1) + fc1(2)*ll3x3%p1%r(1)
        pv(2,1) = pv(2,1) + fc2(2)*ll3x3%p2%r(1)
        pv(2,1) = pv(2,1) + fc3(2)*ll3x3%p3%r(1)
!pv(2,2)
        pv(2,2) = pv(2,2) + fc1(2)*ll3x3%p1%r(2)
        pv(2,2) = pv(2,2) + fc2(2)*ll3x3%p2%r(2)
        pv(2,2) = pv(2,2) + fc3(2)*ll3x3%p3%r(2)
!pv(2,3)
        pv(2,3) = pv(2,3) + fc1(2)*ll3x3%p1%r(3)
        pv(2,3) = pv(2,3) + fc2(2)*ll3x3%p2%r(3)
        pv(2,3) = pv(2,3) + fc3(2)*ll3x3%p3%r(3)
!pv(3,1)
        pv(3,1) = pv(3,1) + fc1(3)*ll3x3%p1%r(1)
        pv(3,1) = pv(3,1) + fc2(3)*ll3x3%p2%r(1)
        pv(3,1) = pv(3,1) + fc3(3)*ll3x3%p3%r(1)
!pv(3,2)
        pv(3,2) = pv(3,2) + fc1(3)*ll3x3%p1%r(2)
        pv(3,2) = pv(3,2) + fc2(3)*ll3x3%p2%r(2)
        pv(3,2) = pv(3,2) + fc3(3)*ll3x3%p3%r(2)
!pv(3,3)
        pv(3,3) = pv(3,3) + fc1(3)*ll3x3%p1%r(3)
        pv(3,3) = pv(3,3) + fc2(3)*ll3x3%p2%r(3)
        pv(3,3) = pv(3,3) + fc3(3)*ll3x3%p3%r(3)
        
        ll3x3 => ll3x3%next
     END DO
     
  END DO
  
#if defined(__parallel)
  CALL mp_sum ( pv, group )
#endif
  
  thermo % pv_const = pv
  
END SUBROUTINE pv_constraint

!******************************************************************************

SUBROUTINE shake_roll(molecule,simpar,r_shake,max_sigma,ishake)
  
! updates the multiplier one molecule type at a time
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: r_shake
  INTEGER, INTENT ( IN ) :: ishake
  
! Locals
  REAL ( dbl ) :: dtsqby2, dtby2, imass1, imass2, fdotf, fdotf_sum
  REAL ( dbl ), DIMENSION (3) :: fc_t, r12, f_roll
  TYPE (linklist_dist_constraints), POINTER :: llconst
  INTEGER :: imol, iconst
  
!------------------------------------------------------------------------------
  
! define local variables
  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
  dtby2 = simpar%dt*.5_dbl
  IF (ishake==1) THEN
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        llconst%fc = -2._dbl*(llconst%rold1-llconst%rold2)
        
! rotate fconst:  
        f_roll = matvec_3x3(r_shake,llconst%fc)
        
! update positions with full multiplier
        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%lambda*f_roll
        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%lambda*f_roll
        
! update velocites with full multiplier
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%lambda*llconst% &
             fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%lambda*llconst% &
             fc
        llconst => llconst%next
     END DO
     
  ELSE
     
     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol
        
! computing the constraint
        r12 = llconst%p1%r - llconst%p2%r
        llconst%sigma = dot_product(r12,r12) - llconst%distance*llconst% &
             distance
        
! computing fc(dt)
        fc_t = -2._dbl*(llconst%p1%r-llconst%p2%r)
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        
! rotate fc: 
        f_roll = matvec_3x3(r_shake,llconst%fc)
        fdotf = dot_product(f_roll,fc_t)
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        llconst%del_lam = 2._dbl*llconst%sigma/(simpar%dt*simpar%dt* &
             fdotf_sum)
        llconst%lambda = llconst%lambda + llconst%del_lam
        
! update positions by delta_lam
        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%del_lam* &
             f_roll
        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%del_lam* &
             f_roll
        
! update velocites by delta_lam
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%del_lam*llconst &
             %fc
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%del_lam*llconst &
             %fc
        
        llconst => llconst%next
     END DO
     
  END IF
  
! computing the constraint and value of tolerance 
  llconst => molecule%ll_dist_const
  DO iconst = 1, molecule%ndistconst_mol
     r12 = llconst%p1%r - llconst%p2%r
     llconst%sigma = dot_product(r12,r12) - llconst%distance*llconst% &
          distance
     max_sigma = max(abs(llconst%sigma),max_sigma)
     llconst => llconst%next
  END DO
END SUBROUTINE shake_roll

!******************************************************************************

SUBROUTINE roll(molecule,veps,veta,wb,dt,nnhc,nfree,eval,u)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: molecule
  REAL ( dbl ), INTENT ( INOUT ), DIMENSION ( :, : ) :: veta, veps, wb
  REAL ( dbl ), INTENT ( IN ) :: dt
  INTEGER, INTENT ( IN ) :: nfree, nnhc
  REAL ( dbl ), INTENT ( IN ), OPTIONAL, DIMENSION ( :, : ) :: u
  REAL ( dbl ), INTENT ( IN ), OPTIONAL, DIMENSION ( : ) :: eval
  
! Locals
  REAL ( dbl ) :: aa, bb, infree, trvg
  REAL ( dbl ), DIMENSION (3,3) :: roll_local, diag
  TYPE (linklist_dist_constraints), POINTER :: llconst
  TYPE (linklist_g3x3_constraints), POINTER :: ll3x3
  INTEGER :: imol, iconst
  
!------------------------------------------------------------------------------
  
  infree = 1.0_dbl / REAL ( nfree, dbl )
  
! updating the scale factor for barostat velocity
  aa = exp(-.125_dbl*dt*veta(1,nnhc))
  bb = exp(-3._dbl*.125_dbl*dt*veta(1,nnhc))
  scale = scale*aa*aa*aa*aa + (aa+bb)*.25_dbl*dt/wb(1,1)
  
! updating the roll matrix for the velocities
  IF ( PRESENT ( u)) THEN
     diag = 0._dbl
     diag(1,1) = eval(1)
     diag(2,2) = eval(2)
     diag(3,3) = eval(3)
     trvg = infree*(veps(1,1)+veps(2,2)+veps(3,3))
     roll_local = matmul_3x3(matmul_3x3(u,diag),transpose_3d(u))
     r_rattle = matmul_3x3(roll_local,r_rattle)
     
     DO imol = 1, size(molecule)
        
! dist
        llconst => molecule(imol) %ll_dist_const
        DO iconst = 1, molecule(imol) %ndistconst_mol
           llconst%scale = llconst%scale*exp(-dt*(veta(1, &
                molecule(imol)%nvt_num)+trvg))
           llconst => llconst%next
        END DO
        
! 3x3
        ll3x3 => molecule(imol) %ll_g3x3_const
        DO iconst = 1, molecule(imol) %ng3x3const_mol
           ll3x3%scale = ll3x3%scale*exp(-dt*(veta(1, &
                molecule(imol)%nvt_num)+trvg))
           ll3x3 => ll3x3%next
        END DO
     END DO
     
  ELSE
     
     DO imol = 1, size(molecule)
! dist
        llconst => molecule(imol) %ll_dist_const
        DO iconst = 1, molecule(imol) %ndistconst_mol
           llconst%scale = llconst%scale*exp(-dt*(veta(1, &
                molecule(imol)%nvt_num)+(1._dbl+3._dbl*infree)*veps(1,1)))
           llconst => llconst%next
        END DO
        
! 3x3
        ll3x3 => molecule(imol) %ll_g3x3_const
        DO iconst = 1, molecule(imol) %ng3x3const_mol
           ll3x3%scale = ll3x3%scale*exp(-dt*(veta(1, &
                molecule(imol)%nvt_num)+(1._dbl+3._dbl*infree)*veps(1,1)))
           ll3x3 => ll3x3%next
        END DO
     END DO
     
  END IF
  
END SUBROUTINE roll

!******************************************************************************

SUBROUTINE check_tol(roll_tol,iroll,char,vector,molecule)
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  INTEGER, INTENT ( INOUT ) :: iroll
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ), OPTIONAL :: vector
  TYPE ( molecule_structure_type ), INTENT ( IN ), DIMENSION ( : ), &
       OPTIONAL :: molecule
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char
  
! Locals
  REAL ( dbl ), SAVE :: scale_old
  REAL ( dbl ), DIMENSION (3), SAVE :: vector_old
  REAL ( dbl ), DIMENSION (3,3), SAVE :: r_rattle_old
  REAL ( dbl ), DIMENSION (3,3) :: diff_rattle
  REAL ( dbl ), DIMENSION (3) :: diff_shake
  REAL ( dbl ) :: local_tol
  TYPE (linklist_dist_constraints), POINTER :: llconst
  TYPE (linklist_g3x3_constraints), POINTER :: ll3x3
  INTEGER :: imol, iconst
  
!------------------------------------------------------------------------------
  
  SELECT CASE (char)
  CASE ( 'SHAKE')
     IF (iroll==1) THEN
        vector_old = vector
        roll_tol = -1.E10_dbl
     ELSE
        roll_tol = 0._dbl
        diff_shake = abs(vector_old-vector)
        local_tol = maxval(diff_shake)
        roll_tol = max(roll_tol,local_tol)
        vector_old = vector
     END IF
     iroll = iroll + 1
     
  CASE ( 'RATTLE')
     IF (iroll==1) THEN
        scale_old = scale
        r_rattle_old = r_rattle
        DO imol = 1, size(molecule)
           
! dist
           llconst => molecule(imol) %ll_dist_const
           DO iconst = 1, molecule(imol) %ndistconst_mol
              llconst%scale_old = llconst%scale
              llconst => llconst%next
           END DO
           
! 3x3
           ll3x3 => molecule(imol) %ll_g3x3_const
           DO iconst = 1, molecule(imol) %ng3x3const_mol
              ll3x3%scale_old = ll3x3%scale
              ll3x3 => ll3x3%next
           END DO
        END DO
        roll_tol = -1.E+10_dbl
        
     ELSE
        
        roll_tol = 0._dbl
        
! compute tolerance on 
        DO imol = 1, size(molecule)
           
! dist
           llconst => molecule(imol) %ll_dist_const
           DO iconst = 1, molecule(imol) %ndistconst_mol
              local_tol = abs(llconst%scale-llconst%scale_old)
              roll_tol = max(roll_tol,local_tol)
              llconst%scale_old = llconst%scale
              llconst => llconst%next
           END DO
           
! 3x3
           ll3x3 => molecule(imol) %ll_g3x3_const
           DO iconst = 1, molecule(imol) %ng3x3const_mol
              local_tol = abs(ll3x3%scale-ll3x3%scale_old)
              roll_tol = max(roll_tol,local_tol)
              ll3x3%scale_old = ll3x3%scale
              ll3x3 => ll3x3%next
           END DO
           
        END DO
        
        diff_rattle = abs(r_rattle-r_rattle_old)
        local_tol = maxval(diff_rattle)
        roll_tol = max(roll_tol,local_tol)
        local_tol = abs(scale-scale_old)
        roll_tol = max(roll_tol,local_tol)
     END IF
     
     iroll = iroll + 1
  END SELECT
  
END SUBROUTINE check_tol

!******************************************************************************

SUBROUTINE initialize_roll(molecule)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: molecule
  
! Locals
  INTEGER :: imol, iconst
  REAL ( dbl ), DIMENSION (3,3) :: unit
  TYPE (linklist_dist_constraints), POINTER :: llconst
  TYPE (linklist_g3x3_constraints), POINTER :: ll3x3
  
!------------------------------------------------------------------------------
  
! initialize unit matrix
  unit = 0._dbl
  unit(1,1) = 1._dbl
  unit(2,2) = 1._dbl
  unit(3,3) = 1._dbl
  
! initialize scale variable for veps
  scale = 0._dbl
  r_rattle = unit
  
! initialize r_rattle
  DO imol = 1, SIZE ( molecule )
     
! dist
     llconst => molecule(imol) %ll_dist_const
     DO iconst = 1, molecule(imol) %ndistconst_mol
        llconst%scale = 1._dbl
        llconst => llconst%next
     END DO
     
! 3x3
     ll3x3 => molecule(imol) %ll_g3x3_const
     DO iconst = 1, molecule(imol) %ng3x3const_mol
        ll3x3%scale = 1._dbl
        ll3x3 => ll3x3%next
     END DO
  END DO
  
END SUBROUTINE initialize_roll

!******************************************************************************

SUBROUTINE update_veps ( molecule, veps, thermo, ensemble, group )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: veps
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  CHARACTER ( LEN = * ), INTENT ( IN ) :: ensemble
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  REAL ( dbl ) :: feps
  
!------------------------------------------------------------------------------
  
! increment pressure pv_const
  CALL pv_constraint ( molecule, thermo, 'DEL', group )
  IF ( ensemble == 'NPT_I' ) THEN
     feps = thermo%pv_const(1,1) &
          + thermo%pv_const(2,2) &
          + thermo%pv_const(3,3)
     veps(1,1) = veps(1,1) + scale * feps
     veps(2,2) = veps(1,1)
     veps(3,3) = veps(1,1)
  ELSE IF ( ensemble == 'NPT_F' ) THEN
     veps ( :, : ) = veps ( :, : ) + scale * thermo % pv_const ( :, : )
  END IF
  
END SUBROUTINE update_veps

!******************************************************************************

SUBROUTINE get_r_shake(r_shake,vector,u)
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: r_shake
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: vector
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ), OPTIONAL :: u
  
! Locals
  INTEGER :: i
  REAL ( dbl ), DIMENSION (3,3) :: diag
  
!------------------------------------------------------------------------------
  
  r_shake = 0.0_dbl
  diag = 0.0_dbl
  
  IF ( PRESENT ( u ) ) THEN
     diag(1,1) = vector(1)
     diag(2,2) = vector(2)
     diag(3,3) = vector(3)
     r_shake = matmul_3x3(matmul_3x3(u,diag),transpose_3d(u))
  ELSE
     DO i = 1, 3
        r_shake(i,i) = vector(i)
     END DO
  END IF
  
END SUBROUTINE get_r_shake

!******************************************************************************

END MODULE constraint
