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
  USE md, ONLY : simulation_parameters_type, virial_type
  USE molecule_types, ONLY : molecule_structure_type, &
       linklist_dist_constraints, linklist_g3x3_constraints
  USE message_passing, ONLY : mp_sum
  USE termination, ONLY : stop_program
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : getinv3x3, dotprod_3d, matvec_3x3, matmul_3x3, transpose_3d

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_control, rattle_control, pv_constraint, getold
  PUBLIC :: shake_roll_control, rattle_roll_control

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
     vector_r, vector_v, group, u )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( virial_type ), INTENT ( INOUT ), OPTIONAL :: thermo
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: vector_r, vector_v
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  INTEGER, INTENT ( INOUT ) :: iroll
  INTEGER, INTENT ( IN ) :: group
  REAL ( dbl ), INTENT ( IN ), OPTIONAL, DIMENSION ( :, : ) :: u

! Locals
  REAL ( dbl ) :: max_sigma
  INTEGER :: imol, iconst, ishake, nmol, handle
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: r_shake, v_shake

!------------------------------------------------------------------------------

  CALL timeset ( 'SHAKE_ROLL', 'I', ' ', handle )

! setting up for roll
  IF ( simpar % ensemble == 'NPT_I' ) THEN
     CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v )
  ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
     CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v, u )
  END IF

! initialize locals
  nmol = SIZE ( molecule )
  MOL: DO imol = 1, nmol

! 3x3
     IF (molecule(imol)%ng3x3const_mol /= 0 ) THEN
        CALL shake_roll_3x3(molecule(imol),simpar,r_shake,v_shake,ishake)
     END IF
     ishake = 0

! dist
     IF (molecule(imol)%ndistconst_mol /= 0 ) THEN
        max_sigma = -1.0E+10_dbl
        CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
           ishake = ishake + 1
           max_sigma = 0._dbl
           CALL shake_roll(molecule(imol),simpar,r_shake,v_shake,max_sigma,ishake)
        END DO CONS
     END IF
  END DO MOL

  CALL check_tol ( roll_tol, iroll, 'SHAKE', r_shake )

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE shake_roll_control

!******************************************************************************

SUBROUTINE rattle_roll_control ( molecule, simpar, thermo, veps, mass, &
     roll_tol, iroll, vector, group, u )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( virial_type ), INTENT ( INOUT ) :: thermo
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: veps
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: mass
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: vector
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ), OPTIONAL :: u
  INTEGER, INTENT ( INOUT ) :: iroll
  INTEGER, INTENT ( IN ) :: group

! Locals
  REAL ( dbl ) :: max_sigma, feps
  INTEGER :: imol, iconst, irattle, nmol, handle
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: r_rattle

!------------------------------------------------------------------------------

  CALL timeset( 'RATTLE_ROLL','I',' ',handle)

! setting up for roll
  IF ( simpar % ensemble == 'NPT_I' ) THEN
     CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector )
  ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
     CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector, u=u )
  END IF

! initialize locals
  irattle = 0
  nmol = SIZE ( molecule )
  max_sigma = -1.0E+10_dbl


! Initialize scale and r_rattle
  CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
     irattle = irattle + 1
     max_sigma = 0._dbl
     MOL: DO imol = 1, nmol
        IF ( molecule ( imol ) % ndistconst_mol /= 0 ) THEN

! dist
           CALL rattle_roll ( molecule ( imol ), simpar, r_rattle, veps, &
                             max_sigma, irattle )
        ELSE IF ( molecule ( imol ) % ng3x3const_mol /= 0 ) THEN

! 3x3
           CALL rattle_roll_3x3 ( molecule ( imol ), simpar, r_rattle, veps, &
                                  max_sigma, irattle )
        END IF
     END DO MOL

! update veps
     IF ( irattle /= 1 ) CALL update_veps ( molecule, &
          veps, mass, thermo, simpar % ensemble, simpar % dt, group )
  END DO CONS

  CALL check_tol ( roll_tol, iroll, 'RATTLE', veps = veps )

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

SUBROUTINE shake_roll_3x3(molecule,simpar,r_shake,v_shake,ishake)

  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: r_shake, v_shake
  INTEGER, INTENT ( OUT ) :: ishake

! Locals
  REAL ( dbl ) :: idtsq, dtsqby2, dtby2, imass1, imass2, imass3
  REAL ( dbl ) :: sigma, max_sigma
  REAL ( dbl ), DIMENSION(3) :: r12, r23, r13, vec, bvec
  REAL ( dbl ), DIMENSION(3) :: r0_12, r0_23, r0_13
  REAL ( dbl ), DIMENSION(3) :: r1, r2, r3, v1, v2, v3, fc1, fc2, fc3,  &
  f_roll1, f_roll2, f_roll3
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
     f_roll1 = matvec_3x3 ( r_shake, llconst % fc1 )
     f_roll2 = matvec_3x3 ( r_shake, llconst % fc2 )
     f_roll3 = matvec_3x3 ( r_shake, llconst % fc3 )

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
        llconst%lambda ( : ) = MATVEC_3x3 ( ainv, bvec )
        fc1=llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(2)*llconst%fc2
        fc2=-llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(3)*llconst%fc3
        fc3=-llconst%lambda(2)*llconst%fc2- &
             llconst%lambda(3)*llconst%fc3
        r1 = llconst%p1%r + imass1*dtsqby2*MATVEC_3X3 ( r_shake, fc1 )
        r2 = llconst%p2%r + imass2*dtsqby2*MATVEC_3X3 ( r_shake, fc2 )
        r3 = llconst%p3%r + imass3*dtsqby2*MATVEC_3X3 ( r_shake, fc3 )
        v1 = llconst%p1%v + imass1*dtby2*MATVEC_3X3 ( v_shake, fc1 )
        v2 = llconst%p2%v + imass2*dtby2*MATVEC_3X3 ( v_shake, fc2 )
        v3 = llconst%p3%v + imass3*dtby2*MATVEC_3X3 ( v_shake, fc3 )
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

SUBROUTINE rattle_roll_3x3(molecule,simpar,r_rattle,veps,max_sigma,irattle)

  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  INTEGER, INTENT ( IN ) :: irattle
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: veps
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: r_rattle

! Locals
  INTEGER :: iconst, info
  REAL ( dbl ) :: idt, dtby2, imass1, imass2, imass3, sigma_dot
  REAL ( dbl ), DIMENSION(3) :: r12, r23, r13, vec, bvec
  REAL ( dbl ), DIMENSION(3) :: f_roll1, f_roll2, f_roll3
  REAL ( dbl ), DIMENSION(3) :: v12_0, v23_0, v13_0
  REAL ( dbl ), DIMENSION(3) :: v12, v23, v13
  REAL ( dbl ), DIMENSION(3) :: fc1, fc2, fc3, lambda
  REAL ( dbl ), DIMENSION(3,3) :: amat, ainv
  TYPE (linklist_g3x3_constraints), POINTER :: llconst

!------------------------------------------------------------------------------

! define local variables
  idt = 1._dbl/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  IF (irattle==1) THEN
     llconst => molecule % ll_g3x3_const
     DO iconst = 1, molecule % ng3x3const_mol
        llconst % vold1 = llconst % p1 % v
        llconst % vold2 = llconst % p2 % v
        llconst % vold3 = llconst % p3 % v
        imass1 = 1._dbl / llconst % p1 % prop % mass
        imass2 = 1._dbl / llconst % p2 % prop % mass
        imass3 = 1._dbl / llconst % p3 % prop % mass
        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
        llconst%fc3 = -2._dbl*(llconst%p2%r-llconst%p3%r)
        fc1=llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(2)*llconst%fc2
        fc2=-llconst%lambda(1)*llconst%fc1+ &
             llconst%lambda(3)*llconst%fc3
        fc3=-llconst%lambda(2)*llconst%fc2- &
             llconst%lambda(3)*llconst%fc3
        llconst%p1%v=llconst%p1%v+imass1*dtby2*MATVEC_3X3 ( r_rattle, fc1 )
        llconst%p2%v=llconst%p2%v+imass2*dtby2*MATVEC_3X3 ( r_rattle, fc2 )
        llconst%p3%v=llconst%p3%v+imass3*dtby2*MATVEC_3X3 ( r_rattle, fc3 )
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
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass
        imass3 = 1._dbl/llconst%p3%prop%mass
        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
        llconst%fc3 = -2._dbl*(llconst%p2%r-llconst%p3%r)

! roll the fc
        f_roll1 = MATVEC_3x3 ( r_rattle, llconst % fc1 )
        f_roll2 = MATVEC_3x3 ( r_rattle, llconst % fc2 )
        f_roll3 = MATVEC_3x3 ( r_rattle, llconst % fc3 )

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
        bvec ( 1 ) = DOTPROD_3D ( r12, v12_0 + MATVEC_3X3 ( veps, r12 ) )
        bvec ( 2 ) = DOTPROD_3D ( r13, v13_0 + MATVEC_3X3 ( veps, r13 ) )
        bvec ( 3 ) = DOTPROD_3D ( r23, v23_0 + MATVEC_3X3 ( veps, r23 ) )
        bvec = -bvec * 2.0_dbl * idt

! find inverse of amat
        CALL getinv3x3 ( amat, ainv )

! get lambda
        lambda ( : ) = MATVEC_3X3 ( ainv, bvec )
        llconst % del_lam = lambda - llconst % lambda

! update velocities
        fc1 = llconst % del_lam ( 1 ) * f_roll1 + &
              llconst % del_lam ( 2 ) * f_roll2
        fc2 = -llconst % del_lam ( 1 ) * f_roll1 + &
               llconst % del_lam ( 3 ) * f_roll3
        fc3 = -llconst % del_lam ( 2 ) * f_roll2 - &
               llconst % del_lam ( 3 ) * f_roll3

        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 * fc1
        llconst % p2 % v = llconst % p2 % v + imass2 * dtby2 * fc2
        llconst % p3 % v = llconst % p3 % v + imass3 * dtby2 * fc3

! computing the convergence of the sigma_dot
        v12 = llconst % p1 % v - llconst % p2 % v
        sigma_dot=2._dbl*DOTPROD_3D(r12,v12+MATVEC_3X3(veps,r12))
        max_sigma = MAX(ABS(sigma_dot),max_sigma)

        v13=llconst%p1%v-llconst%p3%v
        sigma_dot=2._dbl*DOTPROD_3D(r13,v13+MATVEC_3X3(veps,r13))
        max_sigma = MAX(ABS(sigma_dot),max_sigma)

        v23=llconst%p2%v-llconst%p3%v
        sigma_dot=2._dbl*DOTPROD_3D(r23,v23+MATVEC_3X3(veps,r23))
        max_sigma = MAX(ABS(sigma_dot),max_sigma)

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

SUBROUTINE rattle_roll(molecule,simpar,r_rattle,veps,max_sigma,irattle)

! updates the multiplier one molecule type at a time

  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: r_rattle
  INTEGER, INTENT ( IN ) :: irattle
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: veps

! Locals
  REAL ( dbl ) :: dtby2, imass1, imass2, fdotf, fdotf_sum, del_lam
  REAL ( dbl ), DIMENSION (3) :: r12, v12, f_roll
  TYPE (linklist_dist_constraints), POINTER :: llconst
  INTEGER :: imol, iconst, isos

!------------------------------------------------------------------------------

! define local variables
  dtby2 = simpar % dt * 0.5_dbl
  IF ( irattle == 1 ) THEN
     llconst => molecule % ll_dist_const
     DO iconst = 1, molecule % ndistconst_mol
        imass1 = 1._dbl / llconst % p1 % prop % mass
        imass2 = 1._dbl / llconst % p2 % prop % mass
        llconst % fc = -2._dbl * ( llconst % p1 % r - llconst % p2 % r )

! update velocites with full multiplier
        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 *  &
                  llconst % lambda * MATVEC_3X3 ( r_rattle, llconst% fc )
        llconst % p2 % v = llconst % p2 % v - imass2 * dtby2 *  &
                  llconst % lambda * MATVEC_3X3 ( r_rattle, llconst% fc )
        llconst => llconst % next
     END DO
     max_sigma = -1.E+10_dbl

  ELSE
     llconst => molecule % ll_dist_const

     DO iconst = 1, molecule % ndistconst_mol

! assign local variable

        r12 = llconst % p1 % r - llconst % p2 % r

! Roll fc

        llconst % fc = -2._dbl * r12
        f_roll = MATVEC_3x3 ( r_rattle, llconst % fc )

! computing sigma
        v12 ( : ) = llconst % p1 % v + MATVEC_3X3 ( veps, llconst % p1 % r ) &
             - llconst % p2 % v - MATVEC_3X3 ( veps, llconst % p2 % r )
        llconst % sigma = 2._dbl * DOTPROD_3D ( r12, v12 )
        imass1 = 1._dbl / llconst % p1 % prop % mass
        imass2 = 1._dbl / llconst % p2 % prop % mass
        fdotf = DOTPROD_3D (f_roll, llconst % fc )
        fdotf_sum = fdotf * imass1 + fdotf * imass2
        llconst % del_lam = 2._dbl * llconst % sigma / ( simpar % dt * fdotf_sum )
        llconst % lambda = llconst % lambda + llconst % del_lam

!update true velocites with del_lam
        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 *  &
                           llconst % del_lam * f_roll
        llconst % p2 % v = llconst % p2 % v - imass2 * dtby2 *  &
                           llconst % del_lam * f_roll
        llconst => llconst % next
     END DO

! computing the constraint and value of the tolerance
     llconst => molecule % ll_dist_const
     DO iconst = 1, molecule % ndistconst_mol

! computing the updated constraint
        r12 = llconst % p1 % r - llconst % p2 % r
        v12 ( : ) = llconst % p1 % v + MATVEC_3X3 ( veps, llconst % p1 % r ) &
             - llconst % p2 % v - MATVEC_3X3 ( veps, llconst % p2 % r )
        llconst % sigma = 2._dbl * DOTPROD_3D ( r12, v12 )
        max_sigma = MAX ( ABS ( llconst % sigma ), max_sigma )
        llconst => llconst % next
     END DO
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
     ll_dist => molecule(imol) % ll_dist_const
     DO iconst = 1, molecule(imol) % ndistconst_mol
        ll_dist % rold1 ( 1 ) = ll_dist % p1 % r ( 1 )
        ll_dist % rold1 ( 2 ) = ll_dist % p1 % r ( 2 )
        ll_dist % rold1 ( 3 ) = ll_dist % p1 % r ( 3 )
        ll_dist % rold2 ( 1 ) = ll_dist % p2 % r ( 1 )
        ll_dist % rold2 ( 2 ) = ll_dist % p2 % r ( 2 )
        ll_dist % rold2 ( 3 ) = ll_dist % p2 % r ( 3 )
        ll_dist => ll_dist%next
     END DO

     ll_g3x3 => molecule(imol) % ll_g3x3_const
     DO iconst = 1, molecule(imol) % ng3x3const_mol
        ll_g3x3 % rold1 ( 1 ) = ll_g3x3 % p1 % r ( 1 )
        ll_g3x3 % rold1 ( 2 ) = ll_g3x3 % p1 % r ( 2 )
        ll_g3x3 % rold1 ( 3 ) = ll_g3x3 % p1 % r ( 3 )
        ll_g3x3 % rold2 ( 1 ) = ll_g3x3 % p2 % r ( 1 )
        ll_g3x3 % rold2 ( 2 ) = ll_g3x3 % p2 % r ( 2 )
        ll_g3x3 % rold2 ( 3 ) = ll_g3x3 % p2 % r ( 3 )
        ll_g3x3 % rold3 ( 1 ) = ll_g3x3 % p3 % r ( 1 )
        ll_g3x3 % rold3 ( 2 ) = ll_g3x3 % p3 % r ( 2 )
        ll_g3x3 % rold3 ( 3 ) = ll_g3x3 % p3 % r ( 3 )
        ll_g3x3 => ll_g3x3%next
     END DO

  END DO

END SUBROUTINE getold

!******************************************************************************

SUBROUTINE pv_constraint ( molecule, thermo, char, group )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  TYPE ( virial_type ), INTENT ( INOUT ) :: thermo
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

  thermo % pv_constraint = pv

END SUBROUTINE pv_constraint

!******************************************************************************

SUBROUTINE shake_roll(molecule,simpar,r_shake,v_shake,max_sigma,ishake)

! updates the multiplier one molecule type at a time
  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type ), INTENT ( INOUT ) :: molecule
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( INOUT ) :: max_sigma
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: r_shake, v_shake
  INTEGER, INTENT ( IN ) :: ishake

! Locals
  REAL ( dbl ) :: dtsqby2, dtby2, imass1, imass2, fdotf, fdotf_sum
  REAL ( dbl ), DIMENSION (3) :: fc_t, r12, f_roll_r, f_roll_v
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
        f_roll_r = MATVEC_3X3 ( r_shake, llconst % fc )
        f_roll_v = MATVEC_3X3 ( v_shake, llconst % fc )

! update positions with full multiplier
        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%lambda*f_roll_r
        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%lambda*f_roll_r

! update velocites with full multiplier
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%lambda*f_roll_v
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%lambda*f_roll_v
        llconst => llconst%next
     END DO

  ELSE

     llconst => molecule%ll_dist_const
     DO iconst = 1, molecule%ndistconst_mol

! computing the constraint
        r12 = llconst % p1 % r - llconst % p2 % r
        llconst%sigma = DOTPROD_3D ( r12, r12 ) - &
                       llconst % distance * llconst % distance

! computing fc(dt)
        fc_t = -2._dbl*(llconst%p1%r-llconst%p2%r)
        imass1 = 1._dbl/llconst%p1%prop%mass
        imass2 = 1._dbl/llconst%p2%prop%mass

! rotate fc:
        f_roll_r = MATVEC_3X3 ( r_shake, llconst % fc )
        f_roll_v = MATVEC_3X3 ( v_shake, llconst % fc )
        fdotf = DOTPROD_3D ( f_roll_r, fc_t )
        fdotf_sum = fdotf * imass1 + fdotf * imass2
        llconst%del_lam = 2._dbl*llconst%sigma/(simpar%dt*simpar%dt* &
             fdotf_sum)
        llconst % lambda = llconst % lambda + llconst % del_lam

! update positions by delta_lam
        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%del_lam* &
             f_roll_r
        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%del_lam* &
             f_roll_r

! update velocites by delta_lam
        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%del_lam*f_roll_v
        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%del_lam*f_roll_v

        llconst => llconst%next
     END DO

  END IF

! computing the constraint and value of tolerance
  llconst => molecule % ll_dist_const
  DO iconst = 1, molecule % ndistconst_mol
     r12 = llconst % p1 % r - llconst % p2 % r
     llconst%sigma = DOTPROD_3D ( r12, r12 ) -  &
                     llconst % distance * llconst % distance
     max_sigma = MAX ( ABS ( llconst % sigma ), max_sigma )
     llconst => llconst % next
  END DO
END SUBROUTINE shake_roll

!******************************************************************************

SUBROUTINE check_tol(roll_tol,iroll,char,matrix,veps)

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  INTEGER, INTENT ( INOUT ) :: iroll
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ), OPTIONAL :: matrix
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ), OPTIONAL :: veps
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char

! Locals
  REAL ( dbl ), DIMENSION (3,3), SAVE :: matrix_old
  REAL ( dbl ), DIMENSION (3,3), SAVE :: veps_old
  REAL ( dbl ), DIMENSION (3,3) :: diff_rattle
  REAL ( dbl ), DIMENSION (3,3) :: diff_shake
  REAL ( dbl ) :: local_tol
  INTEGER :: imol, iconst

!------------------------------------------------------------------------------

  SELECT CASE ( char )
  CASE ( 'SHAKE' )
     IF (iroll==1) THEN
        matrix_old = matrix
        roll_tol = -1.E10_dbl
     ELSE
        roll_tol = 0._dbl
        diff_shake = ABS ( matrix_old - matrix )
        local_tol = MAXVAL ( diff_shake )
        roll_tol = MAX ( roll_tol, local_tol )
        matrix_old = matrix
     END IF
     iroll = iroll + 1

  CASE ( 'RATTLE' )
     IF (iroll==1) THEN
        veps_old = veps
        roll_tol = -1.E+10_dbl
     ELSE
        roll_tol = 0._dbl

! compute tolerance on veps

        diff_rattle = ABS ( veps - veps_old )
        local_tol = MAXVAL ( diff_rattle )
        roll_tol = MAX ( roll_tol, local_tol )
        veps_old = veps

     END IF
     iroll = iroll + 1

  END SELECT

END SUBROUTINE check_tol

!******************************************************************************

SUBROUTINE get_roll_matrix(char,r_shake,v_shake,vector_r,vector_v,u)

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ), OPTIONAL :: r_shake,v_shake
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ), OPTIONAL :: vector_r, vector_v
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ), OPTIONAL :: u
  CHARACTER ( len = * ), INTENT ( IN ) :: char

! Locals
  INTEGER :: i
  REAL ( dbl ), DIMENSION (3,3) :: diag

!------------------------------------------------------------------------------

  IF ( PRESENT ( r_shake ) ) r_shake = 0.0_dbl
  IF ( PRESENT ( v_shake ) ) v_shake = 0.0_dbl
  diag = 0.0_dbl

  SELECT CASE ( char )
  CASE ( 'SHAKE' )
    IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
                             PRESENT ( vector_r ) ) THEN
       diag ( 1, 1 ) = vector_r ( 1 )
       diag ( 2, 2 ) = vector_r ( 2 )
       diag ( 3, 3 ) = vector_r ( 3 )
       r_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ),TRANSPOSE_3D ( u ) )
       diag ( 1, 1 ) = vector_v ( 1 )
       diag ( 2, 2 ) = vector_v ( 2 )
       diag ( 3, 3 ) = vector_v ( 3 )
       v_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
       diag = matmul_3x3 ( r_shake, v_shake )
       r_shake = diag
    ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
                             PRESENT ( vector_r ) ) THEN
       DO i = 1, 3
          r_shake ( i, i ) = vector_r ( i ) * vector_v ( i )
          v_shake ( i, i ) = vector_v ( i )
       END DO
    ELSE
       CALL stop_program ( "get_roll_matrix", "not sufficient parameters" )
    END IF
  CASE ( 'RATTLE' )
    IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
       diag ( 1, 1 ) = vector_v ( 1 )
       diag ( 2, 2 ) = vector_v ( 2 )
       diag ( 3, 3 ) = vector_v ( 3 )
       v_shake = MATMUL_3x3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
    ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
       DO i = 1, 3
          v_shake ( i, i ) = vector_v ( i )
       END DO
    ELSE
       CALL stop_program ( "get_roll_matrix", "not sufficient parameters" )
    END IF
  END SELECT

END SUBROUTINE get_roll_matrix

!******************************************************************************

SUBROUTINE update_veps ( molecule, veps, mass, thermo, ensemble, dt, group )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: veps
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: mass
  REAL ( dbl ), INTENT ( IN ) :: dt
  TYPE ( virial_type ), INTENT ( INOUT ) :: thermo
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  CHARACTER ( LEN = * ), INTENT ( IN ) :: ensemble
  INTEGER, INTENT ( IN ) :: group

! Locals
  REAL ( dbl ) :: feps

!------------------------------------------------------------------------------

! increment pressure pv_const
  CALL pv_constraint ( molecule, thermo, 'DEL', group )

  IF ( ensemble == 'NPT_I' ) THEN
     feps = thermo%pv_constraint(1,1) &
          + thermo%pv_constraint(2,2) &
          + thermo%pv_constraint(3,3)
     veps(1,1) = veps(1,1) + feps * 0.5_dbl * dt / mass ( 1, 1 )
     veps(2,2) = veps(1,1)
     veps(3,3) = veps(1,1)
  ELSE IF ( ensemble == 'NPT_F' ) THEN
     veps ( :, : ) = veps ( :, : ) +  0.5_dbl *  dt *  &
                     thermo % pv_constraint ( :, : ) / mass ( :, : )
  END IF

END SUBROUTINE update_veps

!******************************************************************************

END MODULE constraint
