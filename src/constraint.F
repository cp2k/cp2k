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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dbl
  USE linear_systems,                  ONLY: solve_system
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             g3x3_constraint_type, &
                                             dist_constraint_type
  USE molecule_types_new,              ONLY: molecule_type, & 
                                             get_molecule, &
                                             local_dist_constraint_type, &
                                             local_g3x3_constraint_type 
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: dotprod_3d,&
                                             matmul_3x3,&
                                             matvec_3x3,&
                                             transpose_3d
  USE virial_types,                    ONLY: virial_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_control, rattle_control,  getold !CJM, pv_constraint
!CJM  PUBLIC :: shake_roll_control, rattle_roll_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE shake_control( local_molecules, molecule_set, molecule_kind_set, &
                          particle_set, simpar )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set ( : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
! locals
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, n3x3con, ndcon
    INTEGER                                  :: handle, ishake, nmol_per_kind
    REAL(dbl)                                :: max_sigma
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind

!------------------------------------------------------------------------------
     CALL timeset ( 'SHAKE_CONTROL', 'I', ' ', handle )
     nkind = SIZE ( molecule_kind_set )
MOL:  DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndcon, ng3x3 = n3x3con )

! 3x3
           IF ( n3x3con /= 0 ) CALL shake_3x3( molecule, particle_set,  &
                                               simpar, ishake )

!! 4x6
!           IF ( n4x6con /= 0 ) CALL shake_4x6( molecule, particle_set, &
!                                               simpar, ishake )
!
           ishake = 0

! dist
          IF ( ndcon /= 0 ) THEN
            max_sigma = -1.0E+10_dbl
            CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
               max_sigma = 0._dbl
               ishake = ishake + 1
               CALL shake( molecule, particle_set, simpar, &
                           max_sigma,ishake )
            END DO CONS
         END IF
       END DO
     END DO MOL

     CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE shake_control

!******************************************************************************

SUBROUTINE rattle_control( local_molecules, molecule_set, molecule_kind_set, &
                          particle_set, simpar )


    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set ( : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
! locals
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, irattle, n3x3con, ndcon
    INTEGER                                  :: handle, ishake, nmol_per_kind
    REAL(dbl)                                :: max_sigma
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind

!------------------------------------------------------------------------------
     CALL timeset ( 'RATTLE_CONTROL', 'I', ' ', handle )
     nkind = SIZE ( molecule_kind_set )
MOL:  DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndcon, ng3x3 = n3x3con )

! 3x3
         IF ( n3x3con /= 0 ) CALL rattle_3x3 ( molecule, particle_set, &
                                               simpar )

!! 4x6
!         IF ( n4x6con /= 0 ) CALL rattle_4x6 ( molecule, particle_set, &
!                                                simpar )

         irattle = 0
! dist
         IF ( ndcon /= 0 ) THEN
            max_sigma = -1.0E+10_dbl
            CONS: DO WHILE ( ABS ( max_sigma ) >= simpar % shake_tol )
               max_sigma = 0.0_dbl
               irattle = irattle + 1
               CALL rattle ( molecule, particle_set, simpar, max_sigma, irattle )
            END DO CONS
         END IF
       END DO
     END DO MOL
     CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE rattle_control

!******************************************************************************

!SUBROUTINE shake_roll_control ( molecule, simpar, thermo, roll_tol, iroll, &
!     vector_r, vector_v, group, u )
!
!
!    TYPE(molecule_structure_type), &
!      DIMENSION(:), INTENT(INOUT)            :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    TYPE(virial_type), INTENT(INOUT), &
!      OPTIONAL                               :: thermo
!    REAL(dbl), INTENT(OUT)                   :: roll_tol
!    INTEGER, INTENT(INOUT)                   :: iroll
!    REAL(dbl), DIMENSION(:), INTENT(IN)      :: vector_r, vector_v
!    INTEGER, INTENT(IN)                      :: group
!    REAL(dbl), DIMENSION(:, :), INTENT(IN), &
!      OPTIONAL                               :: u
!
!    INTEGER                                  :: handle, imol, ishake, nmol
!    REAL(dbl)                                :: max_sigma
!    REAL(dbl), DIMENSION(3, 3)               :: r_shake, v_shake
!
!!------------------------------------------------------------------------------
!
!  CALL timeset ( 'SHAKE_ROLL', 'I', ' ', handle )
!
!! setting up for roll
!  IF ( simpar % ensemble == 'NPT_I' ) THEN
!     CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v )
!  ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
!     CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v, u )
!  END IF
!
!! initialize locals
!  nmol = SIZE ( molecule )
!  MOL: DO imol = 1, nmol
!
!! 3x3
!     IF (molecule(imol)%ng3x3const_mol /= 0 ) THEN
!        CALL shake_roll_3x3(molecule(imol),simpar,r_shake,v_shake,ishake)
!     END IF
!
!! 4x6
!     IF (molecule(imol)%ng4x6const_mol /= 0 ) THEN
!        CALL shake_roll_4x6(molecule(imol),simpar,r_shake,v_shake,ishake)
!     END IF
!     ishake = 0
!
!! dist
!     IF (molecule(imol)%ndistconst_mol /= 0 ) THEN
!        max_sigma = -1.0E+10_dbl
!        CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
!           ishake = ishake + 1
!           max_sigma = 0._dbl
!           CALL shake_roll(molecule(imol),simpar,r_shake,v_shake,max_sigma,ishake)
!        END DO CONS
!     END IF
!  END DO MOL
!
!  CALL check_tol ( roll_tol, iroll, 'SHAKE', r_shake )
!
!  CALL timestop ( 0.0_dbl, handle )
!
!END SUBROUTINE shake_roll_control
!
!!******************************************************************************

!SUBROUTINE rattle_roll_control ( molecule, simpar, thermo, veps, mass, &
!     roll_tol, iroll, vector, group, u )
!
!
!    TYPE(molecule_structure_type), &
!      DIMENSION(:), INTENT(INOUT)            :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    TYPE(virial_type), INTENT(INOUT)         :: thermo
!    REAL(dbl), DIMENSION(:, :), &
!      INTENT(INOUT)                          :: veps
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: mass
!    REAL(dbl), INTENT(OUT)                   :: roll_tol
!    INTEGER, INTENT(INOUT)                   :: iroll
!    REAL(dbl), DIMENSION(:), INTENT(IN)      :: vector
!    INTEGER, INTENT(IN)                      :: group
!    REAL(dbl), DIMENSION(:, :), INTENT(IN), &
!      OPTIONAL                               :: u
!
!    INTEGER                                  :: handle, imol, irattle, nmol
!    REAL(dbl)                                :: max_sigma
!    REAL(dbl), DIMENSION(3, 3)               :: r_rattle
!
!!------------------------------------------------------------------------------
!
!  CALL timeset( 'RATTLE_ROLL','I',' ',handle)
!
!! setting up for roll
!  IF ( simpar % ensemble == 'NPT_I' ) THEN
!     CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector )
!  ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
!     CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector, u=u )
!  END IF
!
!! initialize locals
!  irattle = 0
!  nmol = SIZE ( molecule )
!  max_sigma = -1.0E+10_dbl
!
!
!! Initialize scale and r_rattle
!  CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
!     irattle = irattle + 1
!     max_sigma = 0._dbl
!     MOL: DO imol = 1, nmol
!        IF ( molecule ( imol ) % ndistconst_mol /= 0 ) THEN
!
!! dist
!           CALL rattle_roll ( molecule ( imol ), simpar, r_rattle, veps, &
!                             max_sigma, irattle )
!        ELSE IF ( molecule ( imol ) % ng3x3const_mol /= 0 ) THEN
!
!! 3x3
!           CALL rattle_roll_3x3 ( molecule ( imol ), simpar, r_rattle, veps, &
!                                  max_sigma, irattle )
!        ELSE IF ( molecule ( imol ) % ng4x6const_mol /= 0 ) THEN
!
!! 4x6
!           CALL rattle_roll_4x6 ( molecule ( imol ), simpar, r_rattle, veps, &
!                                  max_sigma, irattle )
!        END IF
!     END DO MOL
!
!! update veps
!    IF ( irattle /= 1 ) CALL update_veps ( molecule, &
!          veps, mass, thermo, simpar % ensemble, simpar % dt, group )
! END DO CONS

! CALL check_tol ( roll_tol, iroll, 'RATTLE', veps = veps )

! CALL timestop ( 0.0_dbl, handle )

!END SUBROUTINE rattle_roll_control

!******************************************************************************

SUBROUTINE shake_3x3( molecule, particle_set, simpar, ishake )

    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE(simulation_parameters_type), &
      INTENT(IN)                                      :: simpar
    INTEGER, INTENT(OUT)                              :: ishake

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ng3x3
    REAL(dbl)                                :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, mass, &
                                                max_sigma, sigma
    REAL(dbl), DIMENSION(3)                  :: fc1, fc2, fc3, r0_12, r0_13, &
                                                r0_23, r1, r12, r13, r2, r23, &
                                                r3, v1, v2, v3, vec
    REAL(dbl), DIMENSION(3, 1)               :: bvec
    REAL(dbl), DIMENSION(3, 3)               :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
  idtsq = 1._dbl/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
  DO iconst = 1, ng3x3
     ishake=0
     index_a = g3x3_list ( iconst ) % a + first_atom -1
     index_b = g3x3_list ( iconst ) % b + first_atom -1
     index_c = g3x3_list ( iconst ) % c + first_atom -1
     r0_12=particle_set ( index_a ) % r - particle_set ( index_b ) % r
     r0_13=particle_set ( index_a ) % r - particle_set ( index_c ) % r
     r0_23=particle_set ( index_b ) % r - particle_set ( index_c ) % r
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1._dbl/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1._dbl/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1._dbl/mass
     lg3x3 ( iconst ) % fa = -2._dbl*( lg3x3 ( iconst ) % ra_old - &
                                       lg3x3 ( iconst ) % rb_old )
     lg3x3 ( iconst ) % fb = -2._dbl*( lg3x3 ( iconst ) % ra_old - &
                                       lg3x3 ( iconst ) % rc_old )
     lg3x3 ( iconst ) % fc = -2._dbl*( lg3x3 ( iconst ) % rb_old - &
                                       lg3x3 ( iconst ) % rc_old )
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,lg3x3 ( iconst ) % fa )
     amat(1,2)=imass1*DOTPROD_3D(r0_12, lg3x3 ( iconst) % fb )
     amat(1,3)=-imass2*DOTPROD_3D(r0_12, lg3x3 ( iconst ) % fc )
     amat(2,1)=imass1*DOTPROD_3D(r0_13, lg3x3 ( iconst ) % fa )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,lg3x3 ( iconst ) % fb )
     amat(2,3)=imass3*DOTPROD_3D(r0_13,lg3x3 ( iconst ) % fc )
     amat(3,1)=-imass2*DOTPROD_3D(r0_23, lg3x3 ( iconst ) % fa )
     amat(3,2)=imass3*DOTPROD_3D(r0_23, lg3x3 ( iconst ) % fb )
     amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,lg3x3 ( iconst ) % fc )

! Iterate until convergence:
     max_sigma = -1.0E+10_dbl
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0._dbl
        ishake=ishake+1
        vec= lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa *(imass1+imass2) + &
             lg3x3 ( iconst ) % lambda ( 2 ) * imass1*lg3x3 ( iconst ) % fb - &
             lg3x3 ( iconst ) % lambda ( 3 ) * imass2*lg3x3 ( iconst ) % fc
        bvec(1,1)=g3x3_list ( iconst ) % dab * g3x3_list ( iconst) % dab &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
        vec= lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa * imass1 + &
             lg3x3 ( iconst ) % lambda ( 2 ) * ( imass1 + imass3 ) * lg3x3 ( iconst ) % fb + &
             lg3x3 ( iconst ) % lambda ( 3 ) * imass3 * lg3x3 ( iconst ) % fc
        bvec(2,1)=g3x3_list ( iconst ) % dac * g3x3_list ( iconst) % dac &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
        vec= -lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa * imass2 + &
             lg3x3 ( iconst ) % lambda ( 2 ) * imass3 * lg3x3 ( iconst ) % fb + &
             lg3x3 ( iconst ) % lambda ( 3 ) * ( imass2 + imass3 ) * lg3x3 ( iconst ) % fc
        bvec(3,1)=g3x3_list ( iconst ) % dbc * g3x3_list ( iconst) % dbc &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
        bvec = bvec * idtsq

! get lambda
        atemp = amat
        CALL solve_system ( atemp, 3, bvec )
        lg3x3 ( iconst ) % lambda ( : ) = bvec ( :, 1 )

        fc1=lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb
        fc2=-lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
        fc3=-lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb - &
            lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
        r1 = particle_set ( index_a ) % r + imass1*dtsqby2*fc1
        r2 = particle_set ( index_b ) % r + imass2*dtsqby2*fc2
        r3 = particle_set ( index_c ) % r + imass3*dtsqby2*fc3
        v1 = particle_set ( index_a ) % v + imass1*dtby2*fc1
        v2 = particle_set ( index_b ) % v + imass2*dtby2*fc2
        v3 = particle_set ( index_c ) % v + imass3*dtby2*fc3
        r12=r1-r2
        r13=r1-r3
        r23=r2-r3

! compute the tolerance:
        sigma = dot_product(r12,r12) - g3x3_list ( iconst ) % dab *  &
                                       g3x3_list ( iconst ) % dab
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r13,r13) -  g3x3_list ( iconst ) % dac * & 
                                        g3x3_list ( iconst ) % dac
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r23,r23) - g3x3_list ( iconst ) % dbc *  &
                                       g3x3_list ( iconst ) % dbc
        max_sigma=MAX(max_sigma,ABS(sigma))
     END DO CONS

! update positions with full multiplier
     particle_set ( index_a ) % r = r1
     particle_set ( index_b ) % r = r2
     particle_set ( index_c ) % r = r3

! update velocites with full multiplier
     particle_set ( index_a ) % v = v1
     particle_set ( index_b ) % v = v2
     particle_set ( index_c ) % v = v3
  END DO

END SUBROUTINE shake_3x3

!******************************************************************************

!SUBROUTINE shake_roll_3x3(molecule,simpar,r_shake,v_shake,ishake)


!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: r_shake, v_shake
!    INTEGER, INTENT(OUT)                     :: ishake
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, dtsqby2, idtsq, &
!                                                imass1, imass2, imass3, mass, &
!                                                max_sigma, sigma
!    REAL(dbl), DIMENSION(3)                  :: f_roll1, f_roll2, f_roll3, &
!                                                fc1, fc2, fc3, r0_12, r0_13, &
!                                                r0_23, r1, r12, r13, r2, r23, &
!                                                r3, v1, v2, v3, vec
!    REAL(dbl), DIMENSION(3, 1)               :: bvec
!    REAL(dbl), DIMENSION(3, 3)               :: amat, atemp
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_g3x3_constraints), POINTER :: llconst
!
!------------------------------------------------------------------------------
! define local variables

!  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
!  idtsq = 1._dbl/simpar%dt/simpar%dt
!  dtby2 = simpar%dt*.5_dbl
!  llconst => molecule%ll_g3x3_const
!
!  DO iconst = 1, molecule%ng3x3const_mol
!     ishake=0
!     r0_12=llconst%p1%r-llconst%p2%r
!     r0_13=llconst%p1%r-llconst%p3%r
!     r0_23=llconst%p2%r-llconst%p3%r
!     atomic_kind=>llconst%p1%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass1 = 1._dbl/mass
!     atomic_kind=>llconst%p2%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass2 = 1._dbl/mass
!     atomic_kind=>llconst%p3%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass3 = 1._dbl/mass
!     llconst%fc1 = -2._dbl*(llconst%rold1-llconst%rold2)
!     llconst%fc2 = -2._dbl*(llconst%rold1-llconst%rold3)
!     llconst%fc3 = -2._dbl*(llconst%rold2-llconst%rold3)
!
!! rotate fconst:
!     f_roll1 = matvec_3x3 ( r_shake, llconst % fc1 )
!     f_roll2 = matvec_3x3 ( r_shake, llconst % fc2 )
!     f_roll3 = matvec_3x3 ( r_shake, llconst % fc3 )
!
!! construct matrix
!     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,f_roll1)
!     amat(1,2)=imass1*DOTPROD_3D(r0_12,f_roll2)
!     amat(1,3)=-imass2*DOTPROD_3D(r0_12,f_roll3)
!     amat(2,1)=imass1*DOTPROD_3D(r0_13,f_roll1)
!     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,f_roll2)
!     amat(2,3)=imass3*DOTPROD_3D(r0_13,f_roll3)
!     amat(3,1)=-imass2*DOTPROD_3D(r0_23,f_roll1)
!     amat(3,2)=imass3*DOTPROD_3D(r0_23,f_roll2)
!     amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,f_roll3)
!
!! Iterate until convergence:
!     max_sigma = -1.0E+10_dbl
!     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
!        max_sigma=0._dbl
!        ishake=ishake+1
!        vec= llconst%lambda(1)*f_roll1*(imass1+imass2) + &
!             llconst%lambda(2)*imass1*f_roll2 - &
!             llconst%lambda(3)*imass2*f_roll3
!        bvec(1,1)=llconst%d12*llconst%d12-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_12,r0_12)
!        vec= llconst%lambda(1)*f_roll1*imass1 + &
!             llconst%lambda(2)*(imass1+imass3)*f_roll2 + &
!             llconst%lambda(3)*imass3*f_roll3
!        bvec(2,1)=llconst%d13*llconst%d13-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_13,r0_13)
!        vec= -llconst%lambda(1)*f_roll1*imass2 + &
!             llconst%lambda(2)*imass3*f_roll2 + &
!             llconst%lambda(3)*(imass2+imass3)*f_roll3
!        bvec(3,1)=llconst%d23*llconst%d23-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_23,r0_23)
!        bvec=bvec*idtsq
!
!! get lambda
!        atemp = amat
!        CALL solve_system ( atemp, 3, bvec )
!        llconst%lambda ( : ) = bvec ( :, 1 )
!
!        fc1=llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(2)*llconst%fc2
!        fc2=-llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(3)*llconst%fc3
!        fc3=-llconst%lambda(2)*llconst%fc2- &
!             llconst%lambda(3)*llconst%fc3
!        r1 = llconst%p1%r + imass1*dtsqby2*MATVEC_3X3 ( r_shake, fc1 )
!        r2 = llconst%p2%r + imass2*dtsqby2*MATVEC_3X3 ( r_shake, fc2 )
!        r3 = llconst%p3%r + imass3*dtsqby2*MATVEC_3X3 ( r_shake, fc3 )
!        v1 = llconst%p1%v + imass1*dtby2*MATVEC_3X3 ( v_shake, fc1 )
!        v2 = llconst%p2%v + imass2*dtby2*MATVEC_3X3 ( v_shake, fc2 )
!        v3 = llconst%p3%v + imass3*dtby2*MATVEC_3X3 ( v_shake, fc3 )
!        r12=r1-r2
!        r13=r1-r3
!        r23=r2-r3
!
!! compute the tolerance:
!        sigma = dot_product(r12,r12) - llconst%d12*llconst%d12
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r13,r13) - llconst%d13*llconst%d13
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r23,r23) - llconst%d23*llconst%d23
!        max_sigma=MAX(max_sigma,ABS(sigma))
!     END DO CONS
!
!! update positions with full multiplier
!     llconst%p1%r = r1
!     llconst%p2%r = r2
!     llconst%p3%r = r3
!
!! update velocites with full multiplier
!     llconst%p1%v = v1
!     llconst%p2%v = v2
!     llconst%p3%v = v3
!     llconst => llconst%next
!
!  END DO
!
!END SUBROUTINE shake_roll_3x3
!
!!******************************************************************************

!SUBROUTINE shake_4x6 ( molecule, simpar, ishake )
!
!
!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    INTEGER, INTENT(OUT)                     :: ishake
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, dtsqby2, idtsq, &
!                                                imass1, imass2, imass3, &
!                                                imass4, mass, max_sigma, sigma
!    REAL(dbl), DIMENSION(3) :: fc1, fc2, fc3, fc4, r0_12, r0_13, r0_14, &
!      r0_23, r0_24, r0_34, r1, r12, r13, r14, r2, r23, r24, r3, r34, r4, v1, &
!      v2, v3, v4, vec
!    REAL(dbl), DIMENSION(6, 1)               :: bvec
!    REAL(dbl), DIMENSION(6, 6)               :: amat, atemp
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_g4x6_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
!  idtsq = 1._dbl/simpar%dt/simpar%dt
!  dtby2 = simpar%dt*.5_dbl
!  llconst => molecule%ll_g4x6_const
!  DO iconst = 1, molecule%ng4x6const_mol
!     ishake=0
!     r0_12=llconst%p1%r-llconst%p2%r
!     r0_13=llconst%p1%r-llconst%p3%r
!     r0_23=llconst%p2%r-llconst%p3%r
!     r0_14=llconst%p1%r-llconst%p4%r
!     r0_24=llconst%p2%r-llconst%p4%r
!     r0_34=llconst%p3%r-llconst%p4%r
!
!     atomic_kind=>llconst%p1%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass1 = 1._dbl/mass
!     atomic_kind=>llconst%p2%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass2 = 1._dbl/mass
!     atomic_kind=>llconst%p3%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass3 = 1._dbl/mass
!     atomic_kind=>llconst%p4%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass4 = 1._dbl/mass
!
!     llconst%fc1 = -2._dbl*(llconst%rold1-llconst%rold2)
!     llconst%fc2 = -2._dbl*(llconst%rold1-llconst%rold3)
!     llconst%fc3 = -2._dbl*(llconst%rold1-llconst%rold4)
!     llconst%fc4 = -2._dbl*(llconst%rold2-llconst%rold3)
!     llconst%fc5 = -2._dbl*(llconst%rold2-llconst%rold4)
!     llconst%fc6 = -2._dbl*(llconst%rold3-llconst%rold4)
!
!! construct matrix
!     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,llconst%fc1)
!     amat(1,2)=imass1*DOTPROD_3D(r0_12,llconst%fc2)
!     amat(1,3)=imass1*DOTPROD_3D(r0_12,llconst%fc3)
!     amat(1,4)=-imass2*DOTPROD_3D(r0_12,llconst%fc4)
!     amat(1,5)=-imass2*DOTPROD_3D(r0_12,llconst%fc5)
!     amat(1,6)=0._dbl
!
!     amat(2,1)=imass1*DOTPROD_3D(r0_13,llconst%fc1)
!     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,llconst%fc2)
!     amat(2,3)=imass1*DOTPROD_3D(r0_13,llconst%fc3)
!     amat(2,4)=imass3*DOTPROD_3D(r0_13,llconst%fc4)
!     amat(2,5)=0._dbl
!     amat(2,6)=-imass3*DOTPROD_3D(r0_13,llconst%fc6)
!
!     amat(3,1)=-imass1*DOTPROD_3D(r0_14,llconst%fc1)
!     amat(3,2)=imass1*DOTPROD_3D(r0_14,llconst%fc2)
!     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r0_14,llconst%fc3)
!     amat(3,4)=0._dbl
!     amat(3,5)=imass4*DOTPROD_3D(r0_14,llconst%fc5)
!     amat(3,6)=imass4*DOTPROD_3D(r0_14,llconst%fc6)
!
!     amat(4,1)=-imass2*DOTPROD_3D(r0_23,llconst%fc1)
!     amat(4,2)=imass3*DOTPROD_3D(r0_23,llconst%fc2)
!     amat(4,3)=0._dbl
!     amat(4,4)=(imass3+imass2)*DOTPROD_3D(r0_23,llconst%fc4)
!     amat(4,5)=imass2*DOTPROD_3D(r0_23,llconst%fc5)
!     amat(4,6)=-imass3*DOTPROD_3D(r0_23,llconst%fc6)
!
!     amat(5,1)=-imass2*DOTPROD_3D(r0_24,llconst%fc1)
!     amat(5,2)=0._dbl
!     amat(5,3)=imass4*DOTPROD_3D(r0_24,llconst%fc3)
!     amat(5,4)=imass2*DOTPROD_3D(r0_24,llconst%fc4)
!     amat(5,5)=(imass4+imass2)*DOTPROD_3D(r0_24,llconst%fc5)
!     amat(5,6)=imass4*DOTPROD_3D(r0_24,llconst%fc6)
!
!     amat(6,1)=0._dbl
!     amat(6,2)=-imass3*DOTPROD_3D(r0_34,llconst%fc2)
!     amat(6,3)=imass4*DOTPROD_3D(r0_34,llconst%fc3)
!     amat(6,4)=-imass3*DOTPROD_3D(r0_34,llconst%fc4)
!     amat(6,5)=imass4*DOTPROD_3D(r0_34,llconst%fc5)
!     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r0_34,llconst%fc6)
!
!
!! Iterate until convergence:
!     max_sigma = -1.0E+10_dbl
!     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
!        max_sigma=0._dbl
!        ishake=ishake+1
!        vec= llconst%lambda(1)*llconst%fc1*(imass1+imass2) + &
!             llconst%lambda(2)*imass1*llconst%fc2 + &
!             llconst%lambda(3)*imass1*llconst%fc3 - &
!             llconst%lambda(4)*imass2*llconst%fc4 - &
!             llconst%lambda(5)*imass2*llconst%fc5 
!        bvec(1,1)=llconst%d12*llconst%d12-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_12,r0_12)
!        vec= llconst%lambda(2)*llconst%fc2*(imass1+imass3) + &
!             llconst%lambda(1)*imass1*llconst%fc1 + &
!             llconst%lambda(3)*imass1*llconst%fc3 + &
!             llconst%lambda(4)*imass3*llconst%fc4 - &
!             llconst%lambda(6)*imass3*llconst%fc6 
!        bvec(2,1)=llconst%d13*llconst%d13-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_13,r0_13)
!        vec= llconst%lambda(3)*llconst%fc3*(imass1+imass4) + &
!             llconst%lambda(1)*imass1*llconst%fc1 + &
!             llconst%lambda(2)*imass1*llconst%fc2 + &
!             llconst%lambda(5)*imass4*llconst%fc5 + &
!             llconst%lambda(6)*imass4*llconst%fc6 
!        bvec(3,1)=llconst%d14*llconst%d14-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_14,r0_14)
!        vec= llconst%lambda(4)*llconst%fc4*(imass2+imass3) - &
!             llconst%lambda(1)*imass2*llconst%fc1 + &
!             llconst%lambda(2)*imass3*llconst%fc2 + &
!             llconst%lambda(5)*imass2*llconst%fc5 - &
!             llconst%lambda(6)*imass3*llconst%fc6 
!        bvec(4,1)=llconst%d23*llconst%d23-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_23,r0_23)
!        vec= llconst%lambda(5)*llconst%fc5*(imass2+imass4) - &
!             llconst%lambda(1)*imass2*llconst%fc1 + &
!             llconst%lambda(3)*imass4*llconst%fc3 + &
!             llconst%lambda(4)*imass2*llconst%fc4 + &
!             llconst%lambda(6)*imass4*llconst%fc6 
!        bvec(5,1)=llconst%d24*llconst%d24-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_24,r0_24)
!        vec= llconst%lambda(6)*llconst%fc6*(imass3+imass4) - &
!             llconst%lambda(2)*imass3*llconst%fc2 + &
!             llconst%lambda(3)*imass4*llconst%fc3 - &
!             llconst%lambda(4)*imass3*llconst%fc4 + &
!             llconst%lambda(5)*imass4*llconst%fc5 
!        bvec(6,1)=llconst%d34*llconst%d34-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_34,r0_34)
!
!        bvec = bvec * idtsq
!
!! get lambda
!        atemp = amat
!        CALL solve_system ( atemp, 6, bvec )
!        llconst%lambda ( : ) = bvec ( :, 1 )
!
!        fc1=llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(2)*llconst%fc2+ &
!             llconst%lambda(3)*llconst%fc3
!        fc2=-llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(4)*llconst%fc4+ &
!             llconst%lambda(5)*llconst%fc5
!        fc3=-llconst%lambda(2)*llconst%fc2- &
!             llconst%lambda(4)*llconst%fc4+ &
!             llconst%lambda(6)*llconst%fc6
!        fc4=-llconst%lambda(3)*llconst%fc3- &
!             llconst%lambda(5)*llconst%fc5- &
!             llconst%lambda(6)*llconst%fc6
!        r1 = llconst%p1%r + imass1*dtsqby2*fc1
!        r2 = llconst%p2%r + imass2*dtsqby2*fc2
!        r3 = llconst%p3%r + imass3*dtsqby2*fc3
!        r4 = llconst%p4%r + imass4*dtsqby2*fc4
!        v1 = llconst%p1%v + imass1*dtby2*fc1
!        v2 = llconst%p2%v + imass2*dtby2*fc2
!        v3 = llconst%p3%v + imass3*dtby2*fc3
!        v4 = llconst%p4%v + imass3*dtby2*fc4
!        r12=r1-r2
!        r13=r1-r3
!        r14=r1-r4
!        r23=r2-r3
!        r24=r2-r4
!        r34=r3-r4
!
!! compute the tolerance:
!        sigma = dot_product(r12,r12) - llconst%d12*llconst%d12
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r13,r13) - llconst%d13*llconst%d13
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r14,r14) - llconst%d14*llconst%d14
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r23,r23) - llconst%d23*llconst%d23
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r24,r24) - llconst%d24*llconst%d24
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r34,r34) - llconst%d34*llconst%d34
!        max_sigma=MAX(max_sigma,ABS(sigma))
!     END DO CONS
!
!! update positions with full multiplier
!     llconst%p1%r = r1
!     llconst%p2%r = r2
!     llconst%p3%r = r3
!     llconst%p4%r = r4
!
!! update velocites with full multiplier
!     llconst%p1%v = v1
!     llconst%p2%v = v2
!     llconst%p3%v = v3
!     llconst%p4%v = v4
!     llconst => llconst%next
!
!  END DO
!
!END SUBROUTINE shake_4x6
!******************************************************************************

!SUBROUTINE shake_roll_4x6(molecule,simpar,r_shake,v_shake,ishake)
!
!
!!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: r_shake, v_shake
!    INTEGER, INTENT(OUT)                     :: ishake
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, dtsqby2, idtsq, &
!                                                imass1, imass2, imass3, &
!                                                imass4, mass, max_sigma, sigma
!    REAL(dbl), DIMENSION(3) :: f_roll1, f_roll2, f_roll3, f_roll4, f_roll5, &
!      f_roll6, fc1, fc2, fc3, fc4, r0_12, r0_13, r0_14, r0_23, r0_24, r0_34, &
!      r1, r12, r13, r14, r2, r23, r24, r3, r34, r4, v1, v2, v3, v4, vec
!    REAL(dbl), DIMENSION(6, 1)               :: bvec
!    REAL(dbl), DIMENSION(6, 6)               :: amat, atemp
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_g4x6_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
!  idtsq = 1._dbl/simpar%dt/simpar%dt
!  dtby2 = simpar%dt*.5_dbl
!  llconst => molecule%ll_g4x6_const
!
!  DO iconst = 1, molecule%ng4x6const_mol
!     ishake=0
!     r0_12=llconst%p1%r-llconst%p2%r
!     r0_13=llconst%p1%r-llconst%p3%r
!     r0_23=llconst%p2%r-llconst%p3%r
!     r0_14=llconst%p1%r-llconst%p4%r
!     r0_24=llconst%p2%r-llconst%p4%r
!     r0_34=llconst%p3%r-llconst%p4%r
!     atomic_kind=>llconst%p1%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass1 = 1._dbl/mass
!     atomic_kind=>llconst%p2%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass2 = 1._dbl/mass
!     atomic_kind=>llconst%p3%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass3 = 1._dbl/mass
!     atomic_kind=>llconst%p4%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass4 = 1._dbl/mass
!
!     llconst%fc1 = -2._dbl*(llconst%rold1-llconst%rold2)
!     llconst%fc2 = -2._dbl*(llconst%rold1-llconst%rold3)
!     llconst%fc3 = -2._dbl*(llconst%rold1-llconst%rold4)
!     llconst%fc4 = -2._dbl*(llconst%rold2-llconst%rold3)
!     llconst%fc5 = -2._dbl*(llconst%rold2-llconst%rold4)
!     llconst%fc6 = -2._dbl*(llconst%rold3-llconst%rold4)
!
!! rotate fconst:
!     f_roll1 = matvec_3x3 ( r_shake, llconst % fc1 )
!     f_roll2 = matvec_3x3 ( r_shake, llconst % fc2 )
!     f_roll3 = matvec_3x3 ( r_shake, llconst % fc3 )
!     f_roll4 = matvec_3x3 ( r_shake, llconst % fc4 )
!     f_roll5 = matvec_3x3 ( r_shake, llconst % fc5 )
!     f_roll6 = matvec_3x3 ( r_shake, llconst % fc6 )
!
!! construct matrix
!     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,f_roll1)
!     amat(1,2)=imass1*DOTPROD_3D(r0_12,f_roll2)
!     amat(1,3)=imass1*DOTPROD_3D(r0_12,f_roll3)
!     amat(1,4)=-imass2*DOTPROD_3D(r0_12,f_roll4)
!     amat(1,5)=-imass2*DOTPROD_3D(r0_12,f_roll5)
!     amat(1,6)=0._dbl
!
!     amat(2,1)=imass1*DOTPROD_3D(r0_13,f_roll1)
!     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,f_roll2)
!     amat(2,3)=imass1*DOTPROD_3D(r0_13,f_roll3)
!     amat(2,4)=imass3*DOTPROD_3D(r0_13,f_roll4)
!     amat(2,5)=0._dbl
!     amat(2,6)=-imass3*DOTPROD_3D(r0_13,f_roll6)
!
!     amat(3,1)=-imass1*DOTPROD_3D(r0_14,f_roll1)
!     amat(3,2)=imass1*DOTPROD_3D(r0_14,f_roll2)
!     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r0_14,f_roll3)
!     amat(3,4)=0._dbl
!     amat(3,5)=imass4*DOTPROD_3D(r0_14,f_roll5)
!     amat(3,6)=imass4*DOTPROD_3D(r0_14,f_roll6)
!
!     amat(4,1)=-imass2*DOTPROD_3D(r0_23,f_roll1)
!     amat(4,2)=imass3*DOTPROD_3D(r0_23,f_roll2)
!     amat(4,3)=0._dbl
!     amat(4,4)=(imass3+imass2)*DOTPROD_3D(r0_23,f_roll4)
!     amat(4,5)=imass2*DOTPROD_3D(r0_23,f_roll5)
!     amat(4,6)=-imass3*DOTPROD_3D(r0_23,f_roll6)
!
!     amat(5,1)=-imass2*DOTPROD_3D(r0_24,f_roll1)
!     amat(5,2)=0._dbl
!     amat(5,3)=imass4*DOTPROD_3D(r0_24,f_roll3)
!     amat(5,4)=imass2*DOTPROD_3D(r0_24,f_roll4)
!     amat(5,5)=(imass4+imass2)*DOTPROD_3D(r0_24,f_roll5)
!     amat(5,6)=imass4*DOTPROD_3D(r0_24,f_roll6)
!
!     amat(6,1)=0._dbl
!     amat(6,2)=-imass3*DOTPROD_3D(r0_34,f_roll2)
!     amat(6,3)=imass4*DOTPROD_3D(r0_34,f_roll3)
!     amat(6,4)=-imass3*DOTPROD_3D(r0_34,f_roll4)
!     amat(6,5)=imass4*DOTPROD_3D(r0_34,f_roll5)
!     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r0_34,f_roll6)
!
!! Iterate until convergence:
!     max_sigma = -1.0E+10_dbl
!     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
!        max_sigma=0._dbl
!        ishake=ishake+1
!        vec= llconst%lambda(1)*f_roll1*(imass1+imass2) + &
!             llconst%lambda(2)*imass1*f_roll2 + &
!             llconst%lambda(3)*imass1*f_roll3 - &
!             llconst%lambda(4)*imass2*f_roll4 - &
!             llconst%lambda(5)*imass2*f_roll5 
!        bvec(1,1)=llconst%d12*llconst%d12-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_12,r0_12)
!        vec= llconst%lambda(2)*f_roll2*(imass1+imass3) + &
!             llconst%lambda(1)*imass1*f_roll1 + &
!             llconst%lambda(3)*imass1*f_roll3 + &
!             llconst%lambda(4)*imass3*f_roll4 - &
!             llconst%lambda(6)*imass3*f_roll6 
!        bvec(2,1)=llconst%d13*llconst%d13-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_13,r0_13)
!        vec= llconst%lambda(3)*f_roll3*(imass1+imass4) + &
!             llconst%lambda(1)*imass1*f_roll1 + &
!             llconst%lambda(2)*imass1*f_roll2 + &
!             llconst%lambda(5)*imass4*f_roll5 + &
!             llconst%lambda(6)*imass4*f_roll6 
!        bvec(3,1)=llconst%d14*llconst%d14-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_14,r0_14)
!        vec= llconst%lambda(4)*f_roll4*(imass2+imass3) - &
!             llconst%lambda(1)*imass2*f_roll1 + &
!             llconst%lambda(2)*imass3*f_roll2 + &
!             llconst%lambda(5)*imass2*f_roll5 - &
!             llconst%lambda(6)*imass3*f_roll6 
!        bvec(4,1)=llconst%d23*llconst%d23-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_23,r0_23)
!        vec= llconst%lambda(5)*f_roll5*(imass2+imass4) - &
!             llconst%lambda(1)*imass2*f_roll1 + &
!             llconst%lambda(3)*imass4*f_roll3 + &
!             llconst%lambda(4)*imass2*f_roll4 + &
!             llconst%lambda(6)*imass4*f_roll6 
!        bvec(5,1)=llconst%d24*llconst%d24-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_24,r0_24)
!        vec= llconst%lambda(6)*f_roll6*(imass3+imass4) - &
!             llconst%lambda(2)*imass3*f_roll2 + &
!             llconst%lambda(3)*imass4*f_roll3 - &
!             llconst%lambda(4)*imass3*f_roll4 + &
!             llconst%lambda(5)*imass4*f_roll5
!        bvec(6,1)=llconst%d34*llconst%d34-dtsqby2*dtsqby2*DOTPROD_3D(vec,vec) &
!             -DOTPROD_3D(r0_34,r0_34)
!
!        bvec=bvec*idtsq
!
!! get lambda
!        atemp = amat
!        CALL solve_system ( atemp, 6, bvec )
!        llconst%lambda ( : ) = bvec ( :, 1 )
!
!        fc1=llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(2)*llconst%fc2+ &
!             llconst%lambda(3)*llconst%fc3
!        fc2=-llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(4)*llconst%fc4+ &
!             llconst%lambda(5)*llconst%fc5
!        fc3=-llconst%lambda(2)*llconst%fc2- &
!             llconst%lambda(4)*llconst%fc4+ &
!             llconst%lambda(6)*llconst%fc6
!        fc4=-llconst%lambda(3)*llconst%fc3- &
!             llconst%lambda(5)*llconst%fc5- &
!             llconst%lambda(6)*llconst%fc6
!        r1 = llconst%p1%r + imass1*dtsqby2*MATVEC_3X3 ( r_shake, fc1 )
!        r2 = llconst%p2%r + imass2*dtsqby2*MATVEC_3X3 ( r_shake, fc2 )
!        r3 = llconst%p3%r + imass3*dtsqby2*MATVEC_3X3 ( r_shake, fc3 )
!        r4 = llconst%p4%r + imass3*dtsqby2*MATVEC_3X3 ( r_shake, fc4 )
!        v1 = llconst%p1%v + imass1*dtby2*MATVEC_3X3 ( v_shake, fc1 )
!        v2 = llconst%p2%v + imass2*dtby2*MATVEC_3X3 ( v_shake, fc2 )
!        v3 = llconst%p3%v + imass3*dtby2*MATVEC_3X3 ( v_shake, fc3 )
!        v4 = llconst%p4%v + imass3*dtby2*MATVEC_3X3 ( v_shake, fc4 )
!        r12=r1-r2
!        r13=r1-r3
!        r14=r1-r4
!        r23=r2-r3
!        r24=r2-r4
!        r34=r3-r4
!
!! compute the tolerance:
!        sigma = dot_product(r12,r12) - llconst%d12*llconst%d12
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r13,r13) - llconst%d13*llconst%d13
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r14,r14) - llconst%d14*llconst%d14
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r23,r23) - llconst%d23*llconst%d23
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r24,r24) - llconst%d24*llconst%d24
!        max_sigma=MAX(max_sigma,ABS(sigma))
!        sigma = dot_product(r34,r34) - llconst%d34*llconst%d34
!        max_sigma=MAX(max_sigma,ABS(sigma))
!     END DO CONS
!
!! update positions with full multiplier
!     llconst%p1%r = r1
!     llconst%p2%r = r2
!     llconst%p3%r = r3
!     llconst%p4%r = r4
!
!! update velocites with full multiplier
!     llconst%p1%v = v1
!     llconst%p2%v = v2
!     llconst%p3%v = v3
!     llconst%p4%v = v4
!     llconst => llconst%next
!
!  END DO
!
!END SUBROUTINE shake_roll_4x6
!
!!******************************************************************************

!SUBROUTINE rattle_roll_3x3(molecule,simpar,r_rattle,veps,max_sigma,irattle)
!
!
!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: r_rattle, veps
!    REAL(dbl), INTENT(INOUT)                 :: max_sigma
!    INTEGER, INTENT(IN)                      :: irattle
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, idt, imass1, imass2, &
!                                                imass3, mass, sigma_dot
!    REAL(dbl), DIMENSION(3)                  :: f_roll1, f_roll2, f_roll3, &
!                                                fc1, fc2, fc3, lambda, r12, &
!                                                r13, r23, v12, v12_0, v13, &
!                                                v13_0, v23, v23_0
!    REAL(dbl), DIMENSION(3, 1)               :: bvec
!    REAL(dbl), DIMENSION(3, 3)               :: amat
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_g3x3_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  idt = 1._dbl/simpar%dt
!  dtby2 = simpar%dt*.5_dbl
!  IF (irattle==1) THEN
!     llconst => molecule % ll_g3x3_const
!     DO iconst = 1, molecule % ng3x3const_mol
!        llconst % vold1 = llconst % p1 % v
!        llconst % vold2 = llconst % p2 % v
!        llconst % vold3 = llconst % p3 % v
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!        atomic_kind=>llconst%p3%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass3 = 1._dbl/mass
!        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
!        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
!        llconst%fc3 = -2._dbl*(llconst%p2%r-llconst%p3%r)
!        fc1=llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(2)*llconst%fc2
!        fc2=-llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(3)*llconst%fc3
!        fc3=-llconst%lambda(2)*llconst%fc2- &
!             llconst%lambda(3)*llconst%fc3
!        llconst%p1%v=llconst%p1%v+imass1*dtby2*MATVEC_3X3 ( r_rattle, fc1 )
!        llconst%p2%v=llconst%p2%v+imass2*dtby2*MATVEC_3X3 ( r_rattle, fc2 )
!        llconst%p3%v=llconst%p3%v+imass3*dtby2*MATVEC_3X3 ( r_rattle, fc3 )
!        llconst => llconst%next
!     END DO
!     max_sigma = -1.0E+10_dbl
!
!  ELSE
!
!     llconst => molecule%ll_g3x3_const
!     DO iconst = 1, molecule%ng3x3const_mol
!        v12_0=llconst%vold1-llconst%vold2
!        v13_0=llconst%vold1-llconst%vold3
!        v23_0=llconst%vold2-llconst%vold3
!        r12=llconst%p1%r-llconst%p2%r
!        r13=llconst%p1%r-llconst%p3%r
!        r23=llconst%p2%r-llconst%p3%r
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!        atomic_kind=>llconst%p3%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass3 = 1._dbl/mass
!        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
!        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
!        llconst%fc3 = -2._dbl*(llconst%p2%r-llconst%p3%r)
!
!! roll the fc
!        f_roll1 = MATVEC_3x3 ( r_rattle, llconst % fc1 )
!        f_roll2 = MATVEC_3x3 ( r_rattle, llconst % fc2 )
!        f_roll3 = MATVEC_3x3 ( r_rattle, llconst % fc3 )
!
!! construct matrix
!        amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,f_roll1)
!        amat(1,2)=imass1*DOTPROD_3D(r12,f_roll2)
!        amat(1,3)=-imass2*DOTPROD_3D(r12,f_roll3)
!        amat(2,1)=imass1*DOTPROD_3D(r13,f_roll1)
!        amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,f_roll2)
!        amat(2,3)=imass3*DOTPROD_3D(r13,f_roll3)
!        amat(3,1)=-imass2*DOTPROD_3D(r23,f_roll1)
!        amat(3,2)=imass3*DOTPROD_3D(r23,f_roll2)
!        amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23,f_roll3)
!
!! construct solution vector
!        bvec ( 1, 1 ) = DOTPROD_3D ( r12, v12_0 + MATVEC_3X3 ( veps, r12 ) )
!        bvec ( 2, 1 ) = DOTPROD_3D ( r13, v13_0 + MATVEC_3X3 ( veps, r13 ) )
!        bvec ( 3, 1 ) = DOTPROD_3D ( r23, v23_0 + MATVEC_3X3 ( veps, r23 ) )
!        bvec = -bvec * 2.0_dbl * idt
!
!! get lambda
!        CALL solve_system ( amat, 3, bvec )
!        lambda ( : ) = bvec ( :, 1 )
!        llconst % del_lam = lambda - llconst % lambda
!
!! update velocities
!        fc1 = llconst % del_lam ( 1 ) * f_roll1 + &
!              llconst % del_lam ( 2 ) * f_roll2
!        fc2 = -llconst % del_lam ( 1 ) * f_roll1 + &
!               llconst % del_lam ( 3 ) * f_roll3
!        fc3 = -llconst % del_lam ( 2 ) * f_roll2 - &
!               llconst % del_lam ( 3 ) * f_roll3
!
!        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 * fc1
!        llconst % p2 % v = llconst % p2 % v + imass2 * dtby2 * fc2
!        llconst % p3 % v = llconst % p3 % v + imass3 * dtby2 * fc3
!
!! computing the convergence of the sigma_dot
!        v12 = llconst % p1 % v - llconst % p2 % v
!        sigma_dot=2._dbl*DOTPROD_3D(r12,v12+MATVEC_3X3(veps,r12))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v13=llconst%p1%v-llconst%p3%v
!        sigma_dot=2._dbl*DOTPROD_3D(r13,v13+MATVEC_3X3(veps,r13))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v23=llconst%p2%v-llconst%p3%v
!        sigma_dot=2._dbl*DOTPROD_3D(r23,v23+MATVEC_3X3(veps,r23))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        llconst%lambda=lambda
!        llconst => llconst%next
!     END DO
!
!  END IF
!
!END SUBROUTINE rattle_roll_3x3
!
!******************************************************************************

SUBROUTINE rattle_3x3( molecule, particle_set, simpar)


    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
    INTENT(IN)                             :: simpar

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ng3x3
    REAL(dbl)                                :: dtby2, idt, imass1, imass2, &
                                                imass3, mass
    REAL(dbl), DIMENSION(3)                  :: fc1, fc2, fc3, r12, r13, r23, &
                                                v12, v13, v23
    REAL(dbl), DIMENSION(3, 1)               :: bvec
    REAL(dbl), DIMENSION(3, 3)               :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

  idt = 1._dbl/simpar%dt
  dtby2 = simpar%dt*.5_dbl
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
  DO iconst = 1, SIZE ( g3x3_list )
     index_a = g3x3_list ( iconst ) % a + first_atom - 1
     index_b = g3x3_list ( iconst ) % b + first_atom - 1
     index_c = g3x3_list ( iconst ) % c + first_atom - 1
     v12=particle_set ( index_a ) % v - particle_set ( index_b ) % v
     v13=particle_set ( index_a ) % v - particle_set ( index_c ) % v
     v23=particle_set ( index_b ) % v - particle_set ( index_c ) % v
     r12=particle_set ( index_a ) % r - particle_set ( index_b ) % r
     r13=particle_set ( index_a ) % r - particle_set ( index_c ) % r
     r23=particle_set ( index_b ) % r - particle_set ( index_c ) % r
     atomic_kind=>particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1._dbl/mass
     atomic_kind=>particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1._dbl/mass
     atomic_kind=>particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1._dbl/mass
     lg3x3 ( iconst ) % fa = -2._dbl* r12
     lg3x3 ( iconst ) % fb = -2._dbl* r13
     lg3x3 ( iconst ) % fc = -2._dbl* r23

! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,lg3x3 ( iconst ) % fa )
     amat(1,2)=imass1*DOTPROD_3D(r12,lg3x3 ( iconst ) % fb )
     amat(1,3)=-imass2*DOTPROD_3D(r12,lg3x3 ( iconst ) % fc )
     amat(2,1)=imass1*DOTPROD_3D(r13,lg3x3 ( iconst ) % fa )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,lg3x3 ( iconst ) % fb )
     amat(2,3)=imass3*DOTPROD_3D(r13,lg3x3 ( iconst ) % fc )
     amat(3,1)=-imass2*DOTPROD_3D(r23,lg3x3 ( iconst ) % fa )
     amat(3,2)=imass3*DOTPROD_3D(r23,lg3x3 ( iconst ) % fb )
     amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23,lg3x3 ( iconst ) % fc )

! construct solution vector
     bvec(1,1)=DOTPROD_3D(r12,v12)
     bvec(2,1)=DOTPROD_3D(r13,v13)
     bvec(3,1)=DOTPROD_3D(r23,v23)
     bvec=-bvec*2._dbl*idt


! get lambda
     CALL solve_system ( amat, 3, bvec )
     lg3x3 ( iconst ) % lambda ( : )=bvec ( :, 1 )

     fc1=lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
         lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb
     fc2=-lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
         lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
     fc3=-lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb - &
         lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
     particle_set ( index_a ) % v = particle_set ( index_a ) % v + imass1 * dtby2 * fc1
     particle_set ( index_b ) % v = particle_set ( index_b ) % v + imass2 * dtby2 * fc2
     particle_set ( index_c ) % v = particle_set ( index_c ) % v + imass3 * dtby2 * fc3
  END DO

END SUBROUTINE rattle_3x3

!******************************************************************************
!SUBROUTINE rattle_4x6(molecule,simpar)
!
!
!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, idt, imass1, imass2, &
!                                                imass3, imass4, mass
!    REAL(dbl), DIMENSION(3)                  :: fc1, fc2, fc3, fc4, r12, r13, &
!                                                r14, r23, r24, r34, v12, v13, &
!                                                v14, v23, v24, v34
!    REAL(dbl), DIMENSION(6, 1)               :: bvec
!    REAL(dbl), DIMENSION(6, 6)               :: amat
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_g4x6_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  idt = 1._dbl/simpar%dt
!  dtby2 = simpar%dt*.5_dbl
!  llconst => molecule%ll_g4x6_const
!  DO iconst = 1, molecule%ng4x6const_mol
!     v12=llconst%p1%v-llconst%p2%v
!     v13=llconst%p1%v-llconst%p3%v
!     v23=llconst%p2%v-llconst%p3%v
!     v14=llconst%p1%v-llconst%p4%v
!     v24=llconst%p2%v-llconst%p4%v
!     v34=llconst%p3%v-llconst%p4%v
!
!     r12=llconst%p1%r-llconst%p2%r
!     r13=llconst%p1%r-llconst%p3%r
!     r23=llconst%p2%r-llconst%p3%r
!     r14=llconst%p1%r-llconst%p4%r
!     r24=llconst%p2%r-llconst%p4%r
!     r34=llconst%p3%r-llconst%p4%r
!
!     atomic_kind=>llconst%p1%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass1 = 1._dbl/mass
!     atomic_kind=>llconst%p2%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass2 = 1._dbl/mass
!     atomic_kind=>llconst%p3%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass3 = 1._dbl/mass
!     atomic_kind=>llconst%p4%atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!     imass4 = 1._dbl/mass
!
!     llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
!     llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
!     llconst%fc3 = -2._dbl*(llconst%p1%r-llconst%p4%r)
!     llconst%fc4 = -2._dbl*(llconst%p2%r-llconst%p3%r)
!     llconst%fc5 = -2._dbl*(llconst%p2%r-llconst%p4%r)
!     llconst%fc6 = -2._dbl*(llconst%p3%r-llconst%p4%r)
!
!! construct matrix
!     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,llconst%fc1)
!     amat(1,2)=imass1*DOTPROD_3D(r12,llconst%fc2)
!     amat(1,3)=imass1*DOTPROD_3D(r12,llconst%fc3)
!     amat(1,4)=-imass2*DOTPROD_3D(r12,llconst%fc4)
!     amat(1,5)=-imass2*DOTPROD_3D(r12,llconst%fc5)
!     amat(1,6)=0._dbl
!
!     amat(2,1)=imass1*DOTPROD_3D(r13,llconst%fc1)
!     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,llconst%fc2)
!     amat(2,3)=imass1*DOTPROD_3D(r13,llconst%fc3)
!     amat(2,4)=imass3*DOTPROD_3D(r13,llconst%fc4)
!     amat(2,5)=0._dbl
!     amat(2,6)=imass3*DOTPROD_3D(r13,llconst%fc6)
!
!     amat(3,1)=imass1*DOTPROD_3D(r14,llconst%fc1)
!     amat(3,2)=imass1*DOTPROD_3D(r14,llconst%fc2)
!     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r14,llconst%fc3)
!     amat(3,4)=0._dbl
!     amat(3,5)=imass4*DOTPROD_3D(r14,llconst%fc5)
!     amat(3,6)=imass4*DOTPROD_3D(r14,llconst%fc6)
!
!     amat(4,1)=-imass2*DOTPROD_3D(r23,llconst%fc1)
!     amat(4,2)=imass3*DOTPROD_3D(r23,llconst%fc2)
!     amat(4,3)=0._dbl
!     amat(4,4)=(imass2+imass3)*DOTPROD_3D(r23,llconst%fc4)
!     amat(4,5)=imass2*DOTPROD_3D(r23,llconst%fc5)
!     amat(4,6)=-imass3*DOTPROD_3D(r23,llconst%fc6)
!
!     amat(5,1)=-imass2*DOTPROD_3D(r24,llconst%fc1)
!     amat(5,2)=0._dbl
!     amat(5,3)=imass4*DOTPROD_3D(r24,llconst%fc3)
!     amat(5,4)=imass2*DOTPROD_3D(r24,llconst%fc4)
!     amat(5,5)=(imass2+imass4)*DOTPROD_3D(r24,llconst%fc5)
!     amat(5,6)=imass4*DOTPROD_3D(r24,llconst%fc6)
!
!     amat(6,1)=0._dbl
!     amat(6,2)=-imass3*DOTPROD_3D(r34,llconst%fc2)
!     amat(6,3)=imass4*DOTPROD_3D(r34,llconst%fc3)
!     amat(6,4)=-imass3*DOTPROD_3D(r34,llconst%fc4)
!     amat(6,5)=imass4*DOTPROD_3D(r34,llconst%fc5)
!     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r34,llconst%fc6)
!
!! construct solution vector
!     bvec(1,1)=DOTPROD_3D(r12,v12)
!     bvec(2,1)=DOTPROD_3D(r13,v13)
!     bvec(3,1)=DOTPROD_3D(r14,v14)
!     bvec(4,1)=DOTPROD_3D(r23,v23)
!     bvec(5,1)=DOTPROD_3D(r24,v24)
!     bvec(6,1)=DOTPROD_3D(r34,v34)
!     bvec=-bvec*2._dbl*idt
!
!! get lambda
!     CALL solve_system ( amat, 6, bvec )
!     llconst%lambda ( : )=bvec(:,1)
!
!     fc1=llconst%lambda(1)*llconst%fc1+ &
!          llconst%lambda(2)*llconst%fc2+ &
!          llconst%lambda(3)*llconst%fc3
!     fc2=-llconst%lambda(1)*llconst%fc1+ &
!          llconst%lambda(4)*llconst%fc4+ &
!          llconst%lambda(5)*llconst%fc5
!     fc3=-llconst%lambda(2)*llconst%fc2- &
!          llconst%lambda(4)*llconst%fc4+ &
!          llconst%lambda(6)*llconst%fc6
!     fc4=-llconst%lambda(3)*llconst%fc3- &
!          llconst%lambda(5)*llconst%fc5- &
!          llconst%lambda(6)*llconst%fc6
!     llconst%p1%v=llconst%p1%v+imass1*dtby2*fc1
!     llconst%p2%v=llconst%p2%v+imass2*dtby2*fc2
!     llconst%p3%v=llconst%p3%v+imass3*dtby2*fc3
!     llconst%p4%v=llconst%p4%v+imass4*dtby2*fc4
!     llconst => llconst%next
!  END DO
!
!END SUBROUTINE rattle_4x6
!!******************************************************************************
!
!!SUBROUTINE rattle_roll_4x6(molecule,simpar,r_rattle,veps,max_sigma,irattle)
!
!
!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: r_rattle, veps
!    REAL(dbl), INTENT(INOUT)                 :: max_sigma
!    INTEGER, INTENT(IN)                      :: irattle
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, idt, imass1, imass2, &
!                                                imass3, imass4, mass, &
!                                                sigma_dot
!    REAL(dbl), DIMENSION(3) :: f_roll1, f_roll2, f_roll3, f_roll4, f_roll5, &
!      f_roll6, fc1, fc2, fc3, fc4, r12, r13, r14, r23, r24, r34, v12, v12_0, &
!      v13, v13_0, v14, v14_0, v23, v23_0, v24, v24_0, v34, v34_0
!    REAL(dbl), DIMENSION(6)                  :: lambda
!    REAL(dbl), DIMENSION(6, 1)               :: bvec
!    REAL(dbl), DIMENSION(6, 6)               :: amat
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_g4x6_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  idt = 1._dbl/simpar%dt
!  dtby2 = simpar%dt*.5_dbl
!  IF (irattle==1) THEN
!     llconst => molecule % ll_g4x6_const
!     DO iconst = 1, molecule % ng4x6const_mol
!        llconst % vold1 = llconst % p1 % v
!        llconst % vold2 = llconst % p2 % v
!        llconst % vold3 = llconst % p3 % v
!        llconst % vold4 = llconst % p4 % v
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!        atomic_kind=>llconst%p3%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass3 = 1._dbl/mass
!        atomic_kind=>llconst%p4%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass4 = 1._dbl/mass
!
!        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
!        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
!        llconst%fc3 = -2._dbl*(llconst%p1%r-llconst%p4%r)
!        llconst%fc4 = -2._dbl*(llconst%p2%r-llconst%p3%r)
!        llconst%fc5 = -2._dbl*(llconst%p2%r-llconst%p4%r)
!        llconst%fc6 = -2._dbl*(llconst%p3%r-llconst%p4%r)
!
!        fc1=llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(2)*llconst%fc2+ &
!             llconst%lambda(3)*llconst%fc3
!        fc2=-llconst%lambda(1)*llconst%fc1+ &
!             llconst%lambda(4)*llconst%fc4+ &
!             llconst%lambda(5)*llconst%fc5
!        fc3=-llconst%lambda(2)*llconst%fc2- &
!             llconst%lambda(4)*llconst%fc4+ &
!             llconst%lambda(6)*llconst%fc6
!        fc4=-llconst%lambda(3)*llconst%fc3- &
!             llconst%lambda(5)*llconst%fc5- &
!             llconst%lambda(6)*llconst%fc6
!
!        llconst%p1%v=llconst%p1%v+imass1*dtby2*MATVEC_3X3 ( r_rattle, fc1 )
!        llconst%p2%v=llconst%p2%v+imass2*dtby2*MATVEC_3X3 ( r_rattle, fc2 )
!        llconst%p3%v=llconst%p3%v+imass3*dtby2*MATVEC_3X3 ( r_rattle, fc3 )
!        llconst%p4%v=llconst%p4%v+imass4*dtby2*MATVEC_3X3 ( r_rattle, fc4 )
!        llconst => llconst%next
!     END DO
!     max_sigma = -1.0E+10_dbl
!
!  ELSE
!
!     llconst => molecule%ll_g4x6_const
!     DO iconst = 1, molecule%ng4x6const_mol
!        v12_0=llconst%vold1-llconst%vold2
!        v13_0=llconst%vold1-llconst%vold3
!        v23_0=llconst%vold2-llconst%vold3
!        v14_0=llconst%vold1-llconst%vold4
!        v24_0=llconst%vold2-llconst%vold4
!        v34_0=llconst%vold3-llconst%vold4
!
!        r12=llconst%p1%r-llconst%p2%r
!        r13=llconst%p1%r-llconst%p3%r
!        r23=llconst%p2%r-llconst%p3%r
!        r14=llconst%p1%r-llconst%p4%r
!        r24=llconst%p2%r-llconst%p4%r
!        r34=llconst%p3%r-llconst%p4%r
!
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!        atomic_kind=>llconst%p3%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass3 = 1._dbl/mass
!        atomic_kind=>llconst%p4%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass4 = 1._dbl/mass
!
!        llconst%fc1 = -2._dbl*(llconst%p1%r-llconst%p2%r)
!        llconst%fc2 = -2._dbl*(llconst%p1%r-llconst%p3%r)
!        llconst%fc3 = -2._dbl*(llconst%p1%r-llconst%p4%r)
!        llconst%fc4 = -2._dbl*(llconst%p2%r-llconst%p3%r)
!        llconst%fc5 = -2._dbl*(llconst%p2%r-llconst%p4%r)
!        llconst%fc6 = -2._dbl*(llconst%p3%r-llconst%p4%r)
!
!! roll the fc
!        f_roll1 = MATVEC_3x3 ( r_rattle, llconst % fc1 )
!        f_roll2 = MATVEC_3x3 ( r_rattle, llconst % fc2 )
!        f_roll3 = MATVEC_3x3 ( r_rattle, llconst % fc3 )
!        f_roll4 = MATVEC_3x3 ( r_rattle, llconst % fc4 )
!        f_roll5 = MATVEC_3x3 ( r_rattle, llconst % fc5 )
!        f_roll6 = MATVEC_3x3 ( r_rattle, llconst % fc6 )
!
!! construct matrix
!        amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,f_roll1)
!        amat(1,2)=imass1*DOTPROD_3D(r12,f_roll2)
!        amat(1,3)=imass1*DOTPROD_3D(r12,f_roll3)
!        amat(1,4)=-imass2*DOTPROD_3D(r12,f_roll4)
!        amat(1,5)=-imass2*DOTPROD_3D(r12,f_roll5)
!        amat(1,6)=0._dbl
!
!        amat(2,1)=imass1*DOTPROD_3D(r13,f_roll1)
!        amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,f_roll2)
!        amat(2,3)=imass1*DOTPROD_3D(r13,f_roll3)
!        amat(2,4)=imass3*DOTPROD_3D(r13,f_roll4)
!        amat(2,5)=0._dbl
!        amat(2,6)=imass3*DOTPROD_3D(r13,f_roll5)
!
!        amat(3,1)=imass1*DOTPROD_3D(r14,f_roll1)
!        amat(3,2)=imass1*DOTPROD_3D(r14,f_roll2)
!        amat(3,3)=(imass1+imass4)*DOTPROD_3D(r14,f_roll3)
!        amat(3,4)=0._dbl
!        amat(3,5)=imass4*DOTPROD_3D(r14,f_roll5)
!        amat(3,6)=imass4*DOTPROD_3D(r14,f_roll6)
!
!        amat(4,1)=-imass2*DOTPROD_3D(r23,f_roll1)
!        amat(4,2)=imass3*DOTPROD_3D(r23,f_roll2)
!        amat(4,3)=0._dbl
!        amat(4,4)=(imass2+imass3)*DOTPROD_3D(r23,f_roll4)
!        amat(4,5)=imass2*DOTPROD_3D(r23,f_roll5)
!        amat(4,6)=-imass3*DOTPROD_3D(r23,f_roll6)
!
!        amat(5,1)=-imass2*DOTPROD_3D(r24,f_roll1)
!        amat(5,2)=0._dbl
!        amat(5,3)=imass4*DOTPROD_3D(r24,f_roll3)
!        amat(5,4)=imass2*DOTPROD_3D(r24,f_roll4)
!        amat(5,5)=(imass2+imass4)*DOTPROD_3D(r24,f_roll5)
!        amat(5,6)=imass4*DOTPROD_3D(r24,f_roll6)
!
!        amat(6,1)=0._dbl
!        amat(6,2)=-imass3*DOTPROD_3D(r34,f_roll2)
!        amat(6,3)=imass4*DOTPROD_3D(r34,f_roll3)
!        amat(6,4)=-imass3*DOTPROD_3D(r34,f_roll4)
!        amat(6,5)=imass4*DOTPROD_3D(r34,f_roll5)
!        amat(6,6)=(imass3+imass4)*DOTPROD_3D(r34,f_roll6)
!
!! construct solution vector
!        bvec ( 1, 1 ) = DOTPROD_3D ( r12, v12_0 + MATVEC_3X3 ( veps, r12 ) )
!        bvec ( 2, 1 ) = DOTPROD_3D ( r13, v13_0 + MATVEC_3X3 ( veps, r13 ) )
!        bvec ( 3, 1 ) = DOTPROD_3D ( r14, v14_0 + MATVEC_3X3 ( veps, r14 ) )
!        bvec ( 4, 1 ) = DOTPROD_3D ( r23, v23_0 + MATVEC_3X3 ( veps, r23 ) )
!        bvec ( 5, 1 ) = DOTPROD_3D ( r24, v24_0 + MATVEC_3X3 ( veps, r24 ) )
!        bvec ( 6, 1 ) = DOTPROD_3D ( r34, v34_0 + MATVEC_3X3 ( veps, r34 ) )
!        bvec = -bvec * 2.0_dbl * idt
!
!
!! get lambda
!        CALL solve_system ( amat, 6, bvec )
!        lambda ( : ) = bvec ( :, 1 )
!        llconst % del_lam = lambda - llconst % lambda
!
!! update velocities
!        fc1=llconst%del_lam(1)*f_roll1+ &
!             llconst%del_lam(2)*f_roll2+ &
!             llconst%del_lam(3)*f_roll3
!        fc2=-llconst%del_lam(1)*f_roll1+ &
!             llconst%del_lam(4)*f_roll4+ &
!             llconst%del_lam(5)*f_roll5
!        fc3=-llconst%del_lam(2)*f_roll2- &
!             llconst%del_lam(4)*f_roll4+ &
!             llconst%del_lam(6)*f_roll6
!        fc4=-llconst%del_lam(3)*f_roll3- &
!             llconst%del_lam(5)*f_roll5- &
!             llconst%del_lam(6)*f_roll6
!
!        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 * fc1
!        llconst % p2 % v = llconst % p2 % v + imass2 * dtby2 * fc2
!        llconst % p3 % v = llconst % p3 % v + imass3 * dtby2 * fc3
!        llconst % p4 % v = llconst % p4 % v + imass4 * dtby2 * fc4
!
!! computing the convergence of the sigma_dot
!        v12 = llconst % p1 % v - llconst % p2 % v
!        sigma_dot=2._dbl*DOTPROD_3D(r12,v12+MATVEC_3X3(veps,r12))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v13=llconst%p1%v-llconst%p3%v
!        sigma_dot=2._dbl*DOTPROD_3D(r13,v13+MATVEC_3X3(veps,r13))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v14=llconst%p1%v-llconst%p4%v
!        sigma_dot=2._dbl*DOTPROD_3D(r14,v14+MATVEC_3X3(veps,r14))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v23=llconst%p2%v-llconst%p3%v
!        sigma_dot=2._dbl*DOTPROD_3D(r23,v23+MATVEC_3X3(veps,r23))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v24=llconst%p2%v-llconst%p4%v
!        sigma_dot=2._dbl*DOTPROD_3D(r24,v24+MATVEC_3X3(veps,r24))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        v34=llconst%p3%v-llconst%p4%v
!        sigma_dot=2._dbl*DOTPROD_3D(r34,v34+MATVEC_3X3(veps,r34))
!        max_sigma = MAX(ABS(sigma_dot),max_sigma)
!
!        llconst%lambda=lambda
!        llconst => llconst%next
!     END DO
!
!  END IF
!
!END SUBROUTINE rattle_roll_4x6

!******************************************************************************
SUBROUTINE shake( molecule, particle_set, simpar, max_sigma, ishake )
! updates the multiplier one molecule type at a time
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    INTEGER, INTENT(INOUT)                     :: ishake
    REAL(dbl), INTENT(INOUT)                 :: max_sigma

! locals

    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ndist
    REAL(dbl)                                :: del_lam, dtby2, dtsqby2, &
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(dbl), DIMENSION(3)                  :: fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )
!
!------------------------------------------------------------------------------
! define local variables

  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
  dtby2 = simpar%dt*.5_dbl
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
  CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
  IF (ishake==1) THEN
     DO iconst = 1, ndist
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1._dbl/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1._dbl/mass
        ldist ( iconst ) % f = -2._dbl * ( ldist ( iconst ) % ra_old - &
                                           ldist ( iconst ) % rb_old )

! update positions with full multiplier
         particle_set ( index_a ) % r = particle_set ( index_a ) % r + &
             imass1*dtsqby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f
         particle_set ( index_b ) % r = particle_set ( index_b ) % r - &
             imass2*dtsqby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f

! update velocites with full multiplier
         particle_set ( index_a ) % v = particle_set ( index_a ) % v + &
             imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f
         particle_set ( index_b ) % v = particle_set ( index_b ) % v - &
             imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f
     END DO

  ELSE

     DO iconst = 1, ndist
! computing the constraint
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        r12 = particle_set ( index_a ) % r -  particle_set ( index_b ) % r
        ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
                                                          dist_list ( iconst ) % dab
! computing fc(dt)
        fc_t = -2._dbl*r12
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1._dbl/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1._dbl/mass
        fdotf = dot_product( ldist ( iconst ) % f , fc_t )
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2._dbl*ldist ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
        ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

! update positions by delta_lam
         particle_set ( index_a ) % r = particle_set ( index_a ) % r + &
             imass1*dtsqby2*del_lam * ldist ( iconst ) % f
         particle_set ( index_b ) % r = particle_set ( index_b ) % r - &
             imass2*dtsqby2*del_lam * ldist ( iconst ) % f

! update velocites by delta_lam
         particle_set ( index_a ) % v = particle_set ( index_a ) % v + &
             imass1*dtby2*del_lam * ldist ( iconst ) % f
         particle_set ( index_b ) % v = particle_set ( index_b ) % v - &
             imass2*dtby2*del_lam * ldist ( iconst ) % f
     END DO
  END IF

! computing the constraint and value of tolerance
  DO iconst = 1, ndist
     index_a = dist_list ( iconst ) % a + first_atom -1
     index_b = dist_list ( iconst ) % b + first_atom -1
     r12 = particle_set ( index_a ) % r -  particle_set ( index_b ) % r
     ldist ( iconst ) % sigma = dot_product(r12,r12) - dist_list ( iconst ) % dab * &
                                                       dist_list ( iconst ) % dab
     max_sigma = MAX(ABS(ldist ( iconst ) % sigma),max_sigma)
  END DO

END SUBROUTINE shake

!!******************************************************************************
SUBROUTINE rattle ( molecule, particle_set, simpar, max_sigma, irattle )
! updates the multiplier one molecule type at a time
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    INTEGER, INTENT(INOUT)                     :: irattle
    REAL(dbl), INTENT(INOUT)                 :: max_sigma

    INTEGER                                  :: iconst, index_a, index_b
    INTEGER                                  :: first_atom, ndist
    REAL(dbl)                                :: del_lam, dtby2, fdotf, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(dbl), DIMENSION(3)                  :: r12, v12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

   dtby2 = simpar%dt*.5_dbl
   molecule_kind => molecule % molecule_kind
   CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
   CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
  IF (irattle==1) THEN
     DO iconst = 1, ndist
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1._dbl/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1._dbl/mass
        ldist ( iconst ) % f = -2._dbl * ( particle_set ( index_a ) % r - &
                                           particle_set ( index_b ) % r )
! update velocites with full multiplier
         particle_set ( index_a ) % v = particle_set ( index_a ) % v + &
             imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f
         particle_set ( index_b ) % v = particle_set ( index_b ) % v - &
             imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f
     END DO

  ELSE

     DO iconst = 1, ndist

! computing sigma
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        r12 = particle_set ( index_a ) % r -  particle_set ( index_b ) % r
        v12 = particle_set ( index_a ) % v -  particle_set ( index_b ) % v
        ldist ( iconst ) % sigma = 2._dbl*dot_product(r12,v12)

! computing f(dt)
        ldist ( iconst ) % f = -2._dbl * r12
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1._dbl/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1._dbl/mass
        fdotf = DOT_PRODUCT( ldist ( iconst ) % f , ldist ( iconst ) % f )
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2._dbl*ldist ( iconst ) % sigma/(simpar%dt*fdotf_sum)
        ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

!update velocites with del_lam
         particle_set ( index_a ) % v = particle_set ( index_a ) % v + &
             imass1*dtby2*del_lam * ldist ( iconst ) % f
         particle_set ( index_b ) % v = particle_set ( index_b ) % v - &
             imass2*dtby2*del_lam * ldist ( iconst ) % f
     END DO
  END IF

! computing the constraint and value of the tolerance
  DO iconst = 1, ndist
     index_a = dist_list ( iconst ) % a + first_atom -1
     index_b = dist_list ( iconst ) % b + first_atom -1
     r12 = particle_set ( index_a ) % r -  particle_set ( index_b ) % r
     v12 = particle_set ( index_a ) % v -  particle_set ( index_b ) % v
     ldist ( iconst ) % sigma = dot_product(r12,v12)
     max_sigma = MAX(ABS(ldist ( iconst ) % sigma),max_sigma)
  END DO

END SUBROUTINE rattle

!******************************************************************************

!SUBROUTINE rattle_roll(molecule,simpar,r_rattle,veps,max_sigma,irattle)
!
!! updates the multiplier one molecule type at a time
!
!
!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: r_rattle, veps
!    REAL(dbl), INTENT(INOUT)                 :: max_sigma
!    INTEGER, INTENT(IN)                      :: irattle
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: del_lam, dtby2, fdotf, &
!                                                fdotf_sum, imass1, imass2, &
!                                                mass
!    REAL(dbl), DIMENSION(3)                  :: f_roll, r12, v12
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_dist_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  dtby2 = simpar % dt * 0.5_dbl
!  IF ( irattle == 1 ) THEN
!     llconst => molecule % ll_dist_const
!     DO iconst = 1, molecule % ndistconst_mol
!
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!
!        llconst % fc = -2._dbl * ( llconst % p1 % r - llconst % p2 % r )
!
!! update velocites with full multiplier
!        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 *  &
!                  llconst % lambda * MATVEC_3X3 ( r_rattle, llconst% fc )
!        llconst % p2 % v = llconst % p2 % v - imass2 * dtby2 *  &
!                  llconst % lambda * MATVEC_3X3 ( r_rattle, llconst% fc )
!        llconst => llconst % next
!     END DO
!     max_sigma = -1.E+10_dbl
!
!  ELSE
!     llconst => molecule % ll_dist_const
!
!     DO iconst = 1, molecule % ndistconst_mol
!
!! assign local variable
!
!        r12 = llconst % p1 % r - llconst % p2 % r
!
!! Roll fc
!
!        llconst % fc = -2._dbl * r12
!        f_roll = MATVEC_3x3 ( r_rattle, llconst % fc )
!
!! computing sigma
!        v12 ( : ) = llconst % p1 % v + MATVEC_3X3 ( veps, llconst % p1 % r ) &
!             - llconst % p2 % v - MATVEC_3X3 ( veps, llconst % p2 % r )
!        llconst % sigma = 2._dbl * DOTPROD_3D ( r12, v12 )
!
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!
!        fdotf = DOTPROD_3D (f_roll, llconst % fc )
!        fdotf_sum = fdotf * imass1 + fdotf * imass2
!        llconst % del_lam = 2._dbl * llconst % sigma / ( simpar % dt * fdotf_sum )
!        llconst % lambda = llconst % lambda + llconst % del_lam
!
!!update true velocites with del_lam
!        llconst % p1 % v = llconst % p1 % v + imass1 * dtby2 *  &
!                           llconst % del_lam * f_roll
!        llconst % p2 % v = llconst % p2 % v - imass2 * dtby2 *  &
!                           llconst % del_lam * f_roll
!        llconst => llconst % next
!     END DO
!
!! computing the constraint and value of the tolerance
!     llconst => molecule % ll_dist_const
!     DO iconst = 1, molecule % ndistconst_mol
!
!! computing the updated constraint
!        r12 = llconst % p1 % r - llconst % p2 % r
!        v12 ( : ) = llconst % p1 % v + MATVEC_3X3 ( veps, llconst % p1 % r ) &
!             - llconst % p2 % v - MATVEC_3X3 ( veps, llconst % p2 % r )
!        llconst % sigma = 2._dbl * DOTPROD_3D ( r12, v12 )
!        max_sigma = MAX ( ABS ( llconst % sigma ), max_sigma )
!        llconst => llconst % next
!     END DO
!  END IF
!
!END SUBROUTINE rattle_roll
!
!!******************************************************************************

SUBROUTINE getold( local_molecules, molecule_set, molecule_kind_set, particle_set )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set ( : )
! locals
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, n3x3con, ndcon, index
    INTEGER                                  :: nmol_per_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( molecule_type ), POINTER     :: molecule

!------------------------------------------------------------------------------
! saves all of the old variables
     nkind = SIZE ( molecule_kind_set )
     DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i )
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist=ndcon, ng3x3=n3x3con, & 
                                    dist_list=dist_list, g3x3_list=g3x3_list )
           CALL get_molecule ( molecule, first_atom=first_atom, &
                               ldist=ldist, lg3x3=lg3x3 )
! distance constraints
           DO iconst = 1, ndcon
             index = dist_list ( iconst ) % a + first_atom - 1
             ldist ( iconst ) % ra_old = particle_set ( index ) % r 
             index = dist_list ( iconst ) % b + first_atom - 1
             ldist ( iconst ) % rb_old = particle_set ( index ) % r 
           ENDDO
! 3x3 constraints
           DO iconst = 1, n3x3con
             index = g3x3_list ( iconst ) % a + first_atom - 1
             lg3x3 ( iconst ) % ra_old = particle_set ( index ) % r 
             index = g3x3_list ( iconst ) % b + first_atom - 1
             lg3x3 ( iconst ) % rb_old = particle_set ( index ) % r 
             index = g3x3_list ( iconst ) % c + first_atom - 1
             lg3x3 ( iconst ) % rc_old = particle_set ( index ) % r 
           ENDDO
! 4x6 constraints
!           DO iconst = 1, n4x6con
!             index = g4x6_list ( iconst ) % a + first_atom - 1
!             lg4x6 ( iconst ) % ra_old = particle_set ( index ) % r 
!             index = g4x6_list ( iconst ) % b + first_atom - 1
!             lg4x6 ( iconst ) % rb_old = particle_set ( index ) % r 
!             index = g4x6_list ( iconst ) % c + first_atom - 1
!             lg4x6 ( iconst ) % rc_old = particle_set ( index ) % r 
!             index = g4x6_list ( iconst ) % d + first_atom - 1
!             lg4x6 ( iconst ) % rd_old = particle_set ( index ) % r 
!           ENDDO
         END DO
       END DO

END SUBROUTINE getold

!******************************************************************************

!SUBROUTINE pv_constraint ( molecule, thermo, char, group )
!
!
!    TYPE(molecule_structure_type), &
!      DIMENSION(:), INTENT(IN)               :: molecule
!    TYPE(virial_type), INTENT(INOUT)         :: thermo
!    CHARACTER(LEN=*), INTENT(IN)             :: char
!    INTEGER, INTENT(IN)                      :: group
!
!    INTEGER                                  :: iconst, imol
!    REAL(dbl) :: fc1( 3 ), fc2( 3 ), fc3( 3 ), fc4( 3 ), lambda, &
!      lambda_3x3( 3 ), lambda_4x6( 6 ), pv( 3, 3 )
!    TYPE(linklist_dist_constraints), POINTER :: llconst
!    TYPE(linklist_g3x3_constraints), POINTER :: ll3x3
!    TYPE(linklist_g4x6_constraints), POINTER :: ll4x6
!
!!------------------------------------------------------------------------------
!
!  pv = 0.0_dbl
!  DO imol = 1, size(molecule)
!
!! dist
!     llconst => molecule(imol) %ll_dist_const
!     DO iconst = 1, molecule(imol) %ndistconst_mol
!
!        SELECT CASE (char)
!        CASE ( 'FULL')
!!  pv gets updated with FULL multiplier
!           lambda = llconst%lambda
!
!        CASE ( 'DEL')
!!  pv gets updated with del_lam
!           lambda = llconst%del_lam
!        END SELECT
!
!!pv(1,1)
!        pv(1,1) = pv(1,1) + lambda*llconst%fc(1)*llconst%p1%r(1)
!        pv(1,1) = pv(1,1) - lambda*llconst%fc(1)*llconst%p2%r(1)
!!pv(1,2)
!        pv(1,2) = pv(1,2) + lambda*llconst%fc(1)*llconst%p1%r(2)
!        pv(1,2) = pv(1,2) - lambda*llconst%fc(1)*llconst%p2%r(2)
!!pv(1,3)
!        pv(1,3) = pv(1,3) + lambda*llconst%fc(1)*llconst%p1%r(3)
!        pv(1,3) = pv(1,3) - lambda*llconst%fc(1)*llconst%p2%r(3)
!!pv(2,1)
!        pv(2,1) = pv(2,1) + lambda*llconst%fc(2)*llconst%p1%r(1)
!        pv(2,1) = pv(2,1) - lambda*llconst%fc(2)*llconst%p2%r(1)
!!pv(2,2)
!        pv(2,2) = pv(2,2) + lambda*llconst%fc(2)*llconst%p1%r(2)
!        pv(2,2) = pv(2,2) - lambda*llconst%fc(2)*llconst%p2%r(2)
!!pv(2,3)
!        pv(2,3) = pv(2,3) + lambda*llconst%fc(2)*llconst%p1%r(3)
!        pv(2,3) = pv(2,3) - lambda*llconst%fc(2)*llconst%p2%r(3)
!!pv(3,1)
!        pv(3,1) = pv(3,1) + lambda*llconst%fc(3)*llconst%p1%r(1)
!        pv(3,1) = pv(3,1) - lambda*llconst%fc(3)*llconst%p2%r(1)
!!pv(3,2)
!        pv(3,2) = pv(3,2) + lambda*llconst%fc(3)*llconst%p1%r(2)
!        pv(3,2) = pv(3,2) - lambda*llconst%fc(3)*llconst%p2%r(2)
!!pv(3,3)
!        pv(3,3) = pv(3,3) + lambda*llconst%fc(3)*llconst%p1%r(3)
!        pv(3,3) = pv(3,3) - lambda*llconst%fc(3)*llconst%p2%r(3)
!        llconst => llconst%next
!     END DO
!
!! 3x3
!     ll3x3 => molecule(imol)%ll_g3x3_const
!     DO iconst = 1, molecule(imol)%ng3x3const_mol
!
!        SELECT CASE (char)
!        CASE ( 'FULL')
!!  pv gets updated with FULL multiplier
!           lambda_3x3 = ll3x3%lambda
!
!        CASE ( 'DEL')
!!  pv gets updated with del_lam
!           lambda_3x3 = ll3x3%del_lam
!        END SELECT
!
!        fc1=lambda_3x3(1)*ll3x3%fc1+ &
!             lambda_3x3(2)*ll3x3%fc2
!        fc2=-lambda_3x3(1)*ll3x3%fc1+ &
!             lambda_3x3(3)*ll3x3%fc3
!        fc3=-lambda_3x3(2)*ll3x3%fc2- &
!             lambda_3x3(3)*ll3x3%fc3
!
!!pv(1,1)
!        pv(1,1) = pv(1,1) + fc1(1)*ll3x3%p1%r(1)
!        pv(1,1) = pv(1,1) + fc2(1)*ll3x3%p2%r(1)
!        pv(1,1) = pv(1,1) + fc3(1)*ll3x3%p3%r(1)
!!pv(1,2)
!        pv(1,2) = pv(1,2) + fc1(1)*ll3x3%p1%r(2)
!        pv(1,2) = pv(1,2) + fc2(1)*ll3x3%p2%r(2)
!        pv(1,2) = pv(1,2) + fc3(1)*ll3x3%p3%r(2)
!!pv(1,3)
!        pv(1,3) = pv(1,3) + fc1(1)*ll3x3%p1%r(3)
!        pv(1,3) = pv(1,3) + fc2(1)*ll3x3%p2%r(3)
!        pv(1,3) = pv(1,3) + fc3(1)*ll3x3%p3%r(3)
!!pv(2,1)
!        pv(2,1) = pv(2,1) + fc1(2)*ll3x3%p1%r(1)
!        pv(2,1) = pv(2,1) + fc2(2)*ll3x3%p2%r(1)
!        pv(2,1) = pv(2,1) + fc3(2)*ll3x3%p3%r(1)
!!pv(2,2)
!        pv(2,2) = pv(2,2) + fc1(2)*ll3x3%p1%r(2)
!        pv(2,2) = pv(2,2) + fc2(2)*ll3x3%p2%r(2)
!        pv(2,2) = pv(2,2) + fc3(2)*ll3x3%p3%r(2)
!!pv(2,3)
!        pv(2,3) = pv(2,3) + fc1(2)*ll3x3%p1%r(3)
!        pv(2,3) = pv(2,3) + fc2(2)*ll3x3%p2%r(3)
!        pv(2,3) = pv(2,3) + fc3(2)*ll3x3%p3%r(3)
!!pv(3,1)
!        pv(3,1) = pv(3,1) + fc1(3)*ll3x3%p1%r(1)
!        pv(3,1) = pv(3,1) + fc2(3)*ll3x3%p2%r(1)
!        pv(3,1) = pv(3,1) + fc3(3)*ll3x3%p3%r(1)
!!pv(3,2)
!        pv(3,2) = pv(3,2) + fc1(3)*ll3x3%p1%r(2)
!        pv(3,2) = pv(3,2) + fc2(3)*ll3x3%p2%r(2)
!        pv(3,2) = pv(3,2) + fc3(3)*ll3x3%p3%r(2)
!!pv(3,3)
!        pv(3,3) = pv(3,3) + fc1(3)*ll3x3%p1%r(3)
!        pv(3,3) = pv(3,3) + fc2(3)*ll3x3%p2%r(3)
!        pv(3,3) = pv(3,3) + fc3(3)*ll3x3%p3%r(3)
!
!        ll3x3 => ll3x3%next
!     END DO
!
!! 4x6
!     ll4x6 => molecule(imol)%ll_g4x6_const
!     DO iconst = 1, molecule(imol)%ng4x6const_mol
!
!        SELECT CASE (char)
!        CASE ( 'FULL')
!!  pv gets updated with FULL multiplier
!           lambda_4x6 = ll4x6%lambda
!
!        CASE ( 'DEL')
!!  pv gets updated with del_lam
!           lambda_4x6 = ll4x6%del_lam
!        END SELECT
!        fc1=lambda_4x6(1)*ll4x6%fc1+ &
!             lambda_4x6(2)*ll4x6%fc2+ &
!             lambda_4x6(3)*ll4x6%fc3
!        fc2=-lambda_4x6(1)*ll4x6%fc1+ &
!             lambda_4x6(4)*ll4x6%fc4+ &
!             lambda_4x6(5)*ll4x6%fc5
!        fc3=-lambda_4x6(2)*ll4x6%fc2- &
!             lambda_4x6(4)*ll4x6%fc4+ &
!             lambda_4x6(6)*ll4x6%fc6
!        fc4=-lambda_4x6(3)*ll4x6%fc3- &
!             lambda_4x6(5)*ll4x6%fc5- &
!             lambda_4x6(6)*ll4x6%fc6
!
!!pv(1,1)
!        pv(1,1) = pv(1,1) + fc1(1)*ll4x6%p1%r(1)
!        pv(1,1) = pv(1,1) + fc2(1)*ll4x6%p2%r(1)
!        pv(1,1) = pv(1,1) + fc3(1)*ll4x6%p3%r(1)
!        pv(1,1) = pv(1,1) + fc4(1)*ll4x6%p4%r(1)
!!pv(1,2)
!        pv(1,2) = pv(1,2) + fc1(1)*ll4x6%p1%r(2)
!        pv(1,2) = pv(1,2) + fc2(1)*ll4x6%p2%r(2)
!        pv(1,2) = pv(1,2) + fc3(1)*ll4x6%p3%r(2)
!        pv(1,2) = pv(1,2) + fc4(1)*ll4x6%p4%r(2)
!!pv(1,3)
!        pv(1,3) = pv(1,3) + fc1(1)*ll4x6%p1%r(3)
!        pv(1,3) = pv(1,3) + fc2(1)*ll4x6%p2%r(3)
!        pv(1,3) = pv(1,3) + fc3(1)*ll4x6%p3%r(3)
!        pv(1,3) = pv(1,3) + fc4(1)*ll4x6%p4%r(3)
!!pv(2,1)
!        pv(2,1) = pv(2,1) + fc1(2)*ll4x6%p1%r(1)
!        pv(2,1) = pv(2,1) + fc2(2)*ll4x6%p2%r(1)
!        pv(2,1) = pv(2,1) + fc3(2)*ll4x6%p3%r(1)
!        pv(2,1) = pv(2,1) + fc4(2)*ll4x6%p4%r(1)
!!pv(2,2)
!        pv(2,2) = pv(2,2) + fc1(2)*ll4x6%p1%r(2)
!        pv(2,2) = pv(2,2) + fc2(2)*ll4x6%p2%r(2)
!        pv(2,2) = pv(2,2) + fc3(2)*ll4x6%p3%r(2)
!        pv(2,2) = pv(2,2) + fc4(2)*ll4x6%p4%r(2)
!!pv(2,3)
!        pv(2,3) = pv(2,3) + fc1(2)*ll4x6%p1%r(3)
!        pv(2,3) = pv(2,3) + fc2(2)*ll4x6%p2%r(3)
!        pv(2,3) = pv(2,3) + fc3(2)*ll4x6%p3%r(3)
!        pv(2,3) = pv(2,3) + fc4(2)*ll4x6%p4%r(3)
!!pv(3,1)
!        pv(3,1) = pv(3,1) + fc1(3)*ll4x6%p1%r(1)
!        pv(3,1) = pv(3,1) + fc2(3)*ll4x6%p2%r(1)
!        pv(3,1) = pv(3,1) + fc3(3)*ll4x6%p3%r(1)
!        pv(3,1) = pv(3,1) + fc4(3)*ll4x6%p4%r(1)
!!pv(3,2)
!        pv(3,2) = pv(3,2) + fc1(3)*ll4x6%p1%r(2)
!        pv(3,2) = pv(3,2) + fc2(3)*ll4x6%p2%r(2)
!        pv(3,2) = pv(3,2) + fc3(3)*ll4x6%p3%r(2)
!        pv(3,2) = pv(3,2) + fc4(3)*ll4x6%p4%r(2)
!!pv(3,3)
!        pv(3,3) = pv(3,3) + fc1(3)*ll4x6%p1%r(3)
!        pv(3,3) = pv(3,3) + fc2(3)*ll4x6%p2%r(3)
!        pv(3,3) = pv(3,3) + fc3(3)*ll4x6%p3%r(3)
!        pv(3,3) = pv(3,3) + fc4(3)*ll4x6%p4%r(3)
!
!        ll4x6 => ll4x6%next
!     END DO
!
!  END DO
!
!#if defined(__parallel)
!  CALL mp_sum ( pv, group )
!#endif
!
!  thermo % pv_constraint = pv
!
!END SUBROUTINE pv_constraint
!
!!******************************************************************************

!SUBROUTINE shake_roll(molecule,simpar,r_shake,v_shake,max_sigma,ishake)
!
!! updates the multiplier one molecule type at a time
!
!    TYPE(molecule_structure_type), &
!      INTENT(INOUT)                          :: molecule
!    TYPE(simulation_parameters_type), &
!      INTENT(IN)                             :: simpar
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: r_shake, v_shake
!    REAL(dbl), INTENT(INOUT)                 :: max_sigma
!    INTEGER, INTENT(IN)                      :: ishake
!
!    INTEGER                                  :: iconst
!    REAL(dbl)                                :: dtby2, dtsqby2, fdotf, &
!                                                fdotf_sum, imass1, imass2, &
!                                                mass
!    REAL(dbl), DIMENSION(3)                  :: f_roll_r, f_roll_v, fc_t, r12
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(linklist_dist_constraints), POINTER :: llconst
!
!!------------------------------------------------------------------------------
!! define local variables
!
!  dtsqby2 = simpar%dt*simpar%dt*.5_dbl
!  dtby2 = simpar%dt*.5_dbl
!  IF (ishake==1) THEN
!     llconst => molecule%ll_dist_const
!     DO iconst = 1, molecule%ndistconst_mol
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!
!        llconst%fc = -2._dbl*(llconst%rold1-llconst%rold2)
!
!! rotate fconst:
!        f_roll_r = MATVEC_3X3 ( r_shake, llconst % fc )
!        f_roll_v = MATVEC_3X3 ( v_shake, llconst % fc )
!
!! update positions with full multiplier
!        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%lambda*f_roll_r
!        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%lambda*f_roll_r
!
!! update velocites with full multiplier
!        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%lambda*f_roll_v
!        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%lambda*f_roll_v
!        llconst => llconst%next
!     END DO
!
!  ELSE
!
!     llconst => molecule%ll_dist_const
!     DO iconst = 1, molecule%ndistconst_mol
!
!! computing the constraint
!        r12 = llconst % p1 % r - llconst % p2 % r
!        llconst%sigma = DOTPROD_3D ( r12, r12 ) - &
!                       llconst % distance * llconst % distance
!
!! computing fc(dt)
!        fc_t = -2._dbl*(llconst%p1%r-llconst%p2%r)
!
!        atomic_kind=>llconst%p1%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass1 = 1._dbl/mass
!        atomic_kind=>llconst%p2%atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
!        imass2 = 1._dbl/mass
!
!! rotate fc:
!        f_roll_r = MATVEC_3X3 ( r_shake, llconst % fc )
!        f_roll_v = MATVEC_3X3 ( v_shake, llconst % fc )
!        fdotf = DOTPROD_3D ( f_roll_r, fc_t )
!        fdotf_sum = fdotf * imass1 + fdotf * imass2
!        llconst%del_lam = 2._dbl*llconst%sigma/(simpar%dt*simpar%dt* &
!             fdotf_sum)
!        llconst % lambda = llconst % lambda + llconst % del_lam
!
!! update positions by delta_lam
!        llconst%p1%r = llconst%p1%r + imass1*dtsqby2*llconst%del_lam* &
!             f_roll_r
!        llconst%p2%r = llconst%p2%r - imass2*dtsqby2*llconst%del_lam* &
!             f_roll_r
!
!! update velocites by delta_lam
!        llconst%p1%v = llconst%p1%v + imass1*dtby2*llconst%del_lam*f_roll_v
!        llconst%p2%v = llconst%p2%v - imass2*dtby2*llconst%del_lam*f_roll_v
!
!        llconst => llconst%next
!     END DO
!
!  END IF
!
!! computing the constraint and value of tolerance
!  llconst => molecule % ll_dist_const
!  DO iconst = 1, molecule % ndistconst_mol
!     r12 = llconst % p1 % r - llconst % p2 % r
!     llconst%sigma = DOTPROD_3D ( r12, r12 ) -  &
!                     llconst % distance * llconst % distance
!     max_sigma = MAX ( ABS ( llconst % sigma ), max_sigma )
!     llconst => llconst % next
!  END DO
!END SUBROUTINE shake_roll
!
!!******************************************************************************
!
!SUBROUTINE check_tol(roll_tol,iroll,char,matrix,veps)


!    REAL(dbl), INTENT(OUT)                   :: roll_tol
!    INTEGER, INTENT(INOUT)                   :: iroll
!    CHARACTER(LEN=*), INTENT(IN)             :: char
!    REAL(dbl), DIMENSION(:, :), INTENT(IN), &
!      OPTIONAL                               :: matrix, veps
!
!    REAL(dbl)                                :: local_tol
!    REAL(dbl), DIMENSION(3, 3)               :: diff_rattle, diff_shake
!    REAL(dbl), DIMENSION(3, 3), SAVE         :: matrix_old, veps_old
!
!!------------------------------------------------------------------------------
!
!  SELECT CASE ( char )
!  CASE ( 'SHAKE' )
!     IF (iroll==1) THEN
!        matrix_old = matrix
!        roll_tol = -1.E10_dbl
!     ELSE
!        roll_tol = 0._dbl
!        diff_shake = ABS ( matrix_old - matrix )
!        local_tol = MAXVAL ( diff_shake )
!        roll_tol = MAX ( roll_tol, local_tol )
!        matrix_old = matrix
!     END IF
!     iroll = iroll + 1
!
!  CASE ( 'RATTLE' )
!     IF (iroll==1) THEN
!        veps_old = veps
!        roll_tol = -1.E+10_dbl
!     ELSE
!        roll_tol = 0._dbl
!
!! compute tolerance on veps
!
!        diff_rattle = ABS ( veps - veps_old )
!        local_tol = MAXVAL ( diff_rattle )
!        roll_tol = MAX ( roll_tol, local_tol )
!        veps_old = veps
!
!     END IF
!     iroll = iroll + 1
!
!  END SELECT
!
!END SUBROUTINE check_tol
!
!!******************************************************************************

!SUBROUTINE get_roll_matrix(char,r_shake,v_shake,vector_r,vector_v,u)
!
!
!    CHARACTER(len=*), INTENT(IN)             :: char
!    REAL(dbl), DIMENSION(:, :), &
!      INTENT(OUT), OPTIONAL                  :: r_shake, v_shake
!    REAL(dbl), DIMENSION(:), INTENT(IN), &
!      OPTIONAL                               :: vector_r, vector_v
!    REAL(dbl), DIMENSION(:, :), INTENT(IN), &
!      OPTIONAL                               :: u
!
!    INTEGER                                  :: i
!    REAL(dbl), DIMENSION(3, 3)               :: diag
!
!!------------------------------------------------------------------------------
!
!  IF ( PRESENT ( r_shake ) ) r_shake = 0.0_dbl
!  IF ( PRESENT ( v_shake ) ) v_shake = 0.0_dbl
!  diag = 0.0_dbl
!
!  SELECT CASE ( char )
!  CASE ( 'SHAKE' )
!    IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
!                             PRESENT ( vector_r ) ) THEN
!       diag ( 1, 1 ) = vector_r ( 1 )
!       diag ( 2, 2 ) = vector_r ( 2 )
!       diag ( 3, 3 ) = vector_r ( 3 )
!       r_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ),TRANSPOSE_3D ( u ) )
!       diag ( 1, 1 ) = vector_v ( 1 )
!       diag ( 2, 2 ) = vector_v ( 2 )
!       diag ( 3, 3 ) = vector_v ( 3 )
!       v_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
!       diag = matmul_3x3 ( r_shake, v_shake )
!       r_shake = diag
!    ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
!                             PRESENT ( vector_r ) ) THEN
!       DO i = 1, 3
!          r_shake ( i, i ) = vector_r ( i ) * vector_v ( i )
!          v_shake ( i, i ) = vector_v ( i )
!       END DO
!    ELSE
!       CALL stop_program ( "get_roll_matrix", "not sufficient parameters" )
!    END IF
!  CASE ( 'RATTLE' )
!    IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
!       diag ( 1, 1 ) = vector_v ( 1 )
!       diag ( 2, 2 ) = vector_v ( 2 )
!       diag ( 3, 3 ) = vector_v ( 3 )
!       v_shake = MATMUL_3x3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
!    ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
!       DO i = 1, 3
!          v_shake ( i, i ) = vector_v ( i )
!       END DO
!    ELSE
!       CALL stop_program ( "get_roll_matrix", "not sufficient parameters" )
!    END IF
!  END SELECT
!
!END SUBROUTINE get_roll_matrix
!
!******************************************************************************

!SUBROUTINE update_veps ( molecule, veps, mass, thermo, ensemble, dt, group )
!
!
!    TYPE(molecule_structure_type), &
!      DIMENSION(:), INTENT(IN)               :: molecule
!    REAL(dbl), DIMENSION(:, :), &
!      INTENT(INOUT)                          :: veps
!    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: mass
!    TYPE(virial_type), INTENT(INOUT)         :: thermo
!    CHARACTER(LEN=*), INTENT(IN)             :: ensemble
!    REAL(dbl), INTENT(IN)                    :: dt
!    INTEGER, INTENT(IN)                      :: group
!
!    REAL(dbl)                                :: feps
!
!!------------------------------------------------------------------------------
!! increment pressure pv_const
!
!  CALL pv_constraint ( molecule, thermo, 'DEL', group )
!
!  IF ( ensemble == 'NPT_I' ) THEN
!     feps = thermo%pv_constraint(1,1) &
!          + thermo%pv_constraint(2,2) &
!          + thermo%pv_constraint(3,3)
!     veps(1,1) = veps(1,1) + feps * 0.5_dbl * dt / mass ( 1, 1 )
!     veps(2,2) = veps(1,1)
!     veps(3,3) = veps(1,1)
!  ELSE IF ( ensemble == 'NPT_F' ) THEN
!     veps ( :, : ) = veps ( :, : ) +  0.5_dbl *  dt *  &
!                     thermo % pv_constraint ( :, : ) / mass ( :, : )
!  END IF
!
!END SUBROUTINE update_veps
!
!******************************************************************************

END MODULE constraint
