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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp, &
                                             dp_size
  USE linear_systems,                  ONLY: solve_system
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_bcast, mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             g3x3_constraint_type, &
                                             g4x6_constraint_type, &
                                             dist_constraint_type,&
                                             fixd_constraint_type
  USE molecule_types_new,              ONLY: molecule_type, & 
                                             get_molecule, &
                                             local_dist_constraint_type, &
                                             local_g3x3_constraint_type, & 
                                             local_g4x6_constraint_type 
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory, stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: dotprod_3d,&
                                             matmul_3x3,&
                                             matvec_3x3,&
                                             transpose_3d
  USE virial_types,                    ONLY: virial_type
  USE force_env_types,                 ONLY: force_env_type,&
                                             force_env_get
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsys_get
  USE particle_list_types,             ONLY: particle_list_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_control, rattle_control,  getold,  pv_constraint
  PUBLIC :: shake_roll_control, rattle_roll_control, fix_atom_control
  CHARACTER(LEN=*), PARAMETER :: moduleN = "constraint"


!!*****
!******************************************************************************

CONTAINS


  SUBROUTINE fix_atom_control( force_env, error)
    IMPLICIT NONE    
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'fix_atom_control', &
         routineP = moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: handle, isubsys
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    INTEGER :: nfixed_atoms, nkind, ikind, i,&
         nparticle_local, iparticle_local, iparticle
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER        :: fixd_list
    TYPE(distribution_1d_type), POINTER      :: local_particles

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    IF (.NOT. failure) THEN
!       WRITE(*,*)"FORCE_ENV%IN_USE",force_env%in_use
       CALL force_env_get(force_env=force_env, subsys=subsys)

       DO isubsys = 1, SIZE(subsys)
          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                             particles=particles, &
                             local_particles=local_particles,&
                             molecule_kinds_new=molecule_kinds )
          
          nkind             =  molecule_kinds % n_els
          molecule_kind_set => molecule_kinds % els
          particle_set      => particles%els
          
          Mol_kind:  DO ikind = 1, nkind
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list )
          
             IF (nfixed_atoms ==0) CYCLE Mol_Kind
             DO i = 1, SIZE(local_particles%n_el)
                nparticle_local   = local_particles%n_el(i)
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(i)%array(iparticle_local)
                   IF (ANY(fixd_list(:)%fixd == iparticle)) particle_set(iparticle)%f(:) = 0.0_dp
                END DO
             END DO
          END DO Mol_kind

       END DO
    END IF
!    WRITE(*,'(3f15.6)')(particle_set(i)%f(:),i=1,SIZE(particle_set))
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE fix_atom_control

!******************************************************************************

SUBROUTINE shake_control( local_molecules, molecule_set, molecule_kind_set, &
                          particle_set, pos, vel, simpar )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE ( particle_type ), POINTER          :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )           :: pos ( :, : ), vel ( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
! locals
    CHARACTER(LEN=*), PARAMETER              :: routineN = "shake_control"
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, n3x3con, n4x6con, ndcon
    INTEGER                                  :: handle, ishake, nmol_per_kind
    INTEGER                                  :: isos
    REAL(KIND=dp)                                :: max_sigma
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( atomic_kind_type ), POINTER       :: atomic_kind

!------------------------------------------------------------------------------
     CALL timeset ( 'SHAKE_CONTROL', 'I', ' ', handle )
     nkind = SIZE ( molecule_kind_set )

MOL:  DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndcon,  &
                                    ng3x3 = n3x3con, ng4x6=n4x6con )

! 3x3
           IF ( n3x3con /= 0 ) CALL shake_3x3( molecule, particle_set,  &
                                               pos, vel, simpar, ishake )

! 4x6
           IF ( n4x6con /= 0 ) CALL shake_4x6( molecule, particle_set, &
                                               pos, vel, simpar, ishake )

           ishake = 0

! dist
          IF ( ndcon /= 0 ) THEN
            max_sigma = -1.0E+10_dp
            CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
               max_sigma = 0.0_dp
               ishake = ishake + 1
               CALL shake( molecule, particle_set, pos, vel, simpar, &
                           max_sigma,ishake )
            END DO CONS
         END IF
       END DO
     END DO MOL

     CALL timestop ( 0.0_dp, handle )

END SUBROUTINE shake_control

!******************************************************************************

SUBROUTINE rattle_control( local_molecules, molecule_set, molecule_kind_set,  &
                           particle_set, vel, simpar )


    TYPE ( particle_type ), POINTER          :: particle_set ( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    REAL (KIND=dp), INTENT ( INOUT )           :: vel ( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
! locals
    CHARACTER(LEN=*), PARAMETER              :: routineN = "rattle_control"
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, irattle, n3x3con, n4x6con, ndcon
    INTEGER                                  :: handle, ishake, nmol_per_kind
    INTEGER                                  :: isos 
    REAL(KIND=dp)                                :: max_sigma
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( atomic_kind_type ), POINTER       :: atomic_kind

!------------------------------------------------------------------------------
     CALL timeset ( 'RATTLE_CONTROL', 'I', ' ', handle )
     nkind = SIZE ( molecule_kind_set )
MOL: DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndcon, ng3x3 = n3x3con, ng4x6 = n4x6con )

! 3x3
         IF ( n3x3con /= 0 ) CALL rattle_3x3 ( molecule, particle_set, vel, simpar )

! 4x6
         IF ( n4x6con /= 0 ) CALL rattle_4x6 ( molecule, particle_set, vel, simpar )

         irattle = 0
! dist
         IF ( ndcon /= 0 ) THEN
            max_sigma = -1.0E+10_dp
            CONS: DO WHILE ( ABS ( max_sigma ) >= simpar % shake_tol )
               max_sigma = 0.0_dp
               irattle = irattle + 1
               CALL rattle ( molecule, particle_set, vel, simpar, max_sigma, irattle )
            END DO CONS
         END IF
       END DO
     END DO MOL
     CALL timestop ( 0.0_dp, handle )

END SUBROUTINE rattle_control

!******************************************************************************

SUBROUTINE shake_roll_control( local_molecules, molecule_set, molecule_kind_set, &
                               particle_set, pos, vel, simpar, roll_tol, iroll,  &
                               vector_r, vector_v, group, u )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE ( particle_type ), POINTER          :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )           :: pos ( :, : ), vel ( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(OUT)                   :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: vector_r, vector_v
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: u

! locals
    CHARACTER(LEN=*), PARAMETER              :: routineN = "shake_roll_control"
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, n3x3con, n4x6con, ndcon
    INTEGER                                  :: handle, ishake, nmol_per_kind
    INTEGER                                  :: isos
    REAL(KIND=dp)                                :: max_sigma
    REAL (KIND=dp), DIMENSION ( 3, 3 )         :: r_shake, v_shake
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( atomic_kind_type ), POINTER       :: atomic_kind

!------------------------------------------------------------------------------

  CALL timeset ( 'SHAKE_ROLL_CONTROL', 'I', ' ', handle )
  nkind = SIZE ( molecule_kind_set )
! setting up for roll
  IF ( simpar % ensemble == 'NPT_I' ) THEN
     CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v )
  ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
     CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v, u )
  END IF

MOL:  DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndcon,  &
                                    ng3x3 = n3x3con, ng4x6=n4x6con )

! 3x3
           IF ( n3x3con /= 0 ) CALL shake_roll_3x3( molecule, particle_set,  &
                                               pos, vel, r_shake, v_shake, &
                                               simpar, ishake )

! 4x6
           IF ( n4x6con /= 0 ) CALL shake_roll_4x6( molecule, particle_set, &
                                               pos, vel, r_shake, v_shake, &
                                               simpar, ishake )

           ishake = 0

! dist
          IF ( ndcon /= 0 ) THEN
            max_sigma = -1.0E+10_dp
            CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
               max_sigma = 0.0_dp
               ishake = ishake + 1
               CALL shake_roll( molecule, particle_set, pos, vel, &
                                r_shake, v_shake, simpar, max_sigma,ishake )
            END DO CONS
         END IF
       END DO
     END DO MOL

     CALL check_tol ( roll_tol, iroll, 'SHAKE', r_shake )

     CALL timestop ( 0.0_dp, handle )

END SUBROUTINE shake_roll_control

!!******************************************************************************

SUBROUTINE rattle_roll_control( local_molecules, molecule_set, molecule_kind_set, &
                                particle_set, vel, simpar, virial, vector,   &
                                veps, mass, roll_tol, iroll, para_env, u )


    TYPE ( particle_type ), POINTER          :: particle_set ( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE ( virial_type ), POINTER            :: virial
    TYPE ( cp_para_env_type ), INTENT ( IN ) :: para_env
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL (KIND=dp), INTENT ( INOUT )           :: vel ( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: veps
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: mass
    REAL(KIND=dp), INTENT(OUT)                   :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: vector
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: u

! locals
    CHARACTER(LEN=*), PARAMETER              :: routineN = "rattle_roll_control"
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, irattle, n3x3con, n4x6con, ndcon
    INTEGER                                  :: handle, ishake, nmol_per_kind
    INTEGER                                  :: isos
    REAL(KIND=dp)                                :: max_sigma
    REAL(KIND=dp), DIMENSION(3, 3)               :: r_rattle
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( atomic_kind_type ), POINTER       :: atomic_kind

!------------------------------------------------------------------------------
   CALL timeset ( 'RATTLE_ROLL_CONTROL', 'I', ' ', handle )
! initialize locals
   nkind = SIZE ( molecule_kind_set )
! setting up for roll
   IF ( simpar % ensemble == 'NPT_I' ) THEN
      CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector )
   ELSE IF ( simpar % ensemble == 'NPT_F' ) THEN
      CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector, u=u )
   END IF
   MOL: DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndcon,  &
                                    ng3x3 = n3x3con, ng4x6 = n4x6con )
! 3x3
           IF ( n3x3con /= 0 ) CALL rattle_roll_3x3 ( molecule, particle_set, vel, r_rattle, &
                                                      simpar, veps )

! 4x6
           IF ( n4x6con /= 0 ) CALL rattle_roll_4x6 ( molecule, particle_set, vel, r_rattle, &
                                                      simpar, veps )
           irattle = 0 
! dist
           IF ( ndcon /= 0 )  THEN
             max_sigma = -1.0E+10_dp
             CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
                max_sigma = 0.0_dp
                irattle = irattle + 1
                CALL rattle_roll ( molecule, particle_set, vel, r_rattle, &
                                   simpar, veps, max_sigma, irattle )
             END DO CONS
           ENDIF
        END DO
     END DO MOL


 CALL check_tol ( roll_tol, iroll, 'RATTLE', veps = veps )

 CALL timestop ( 0.0_dp, handle )

END SUBROUTINE rattle_roll_control

!******************************************************************************

SUBROUTINE shake_3x3( molecule, particle_set, pos, vel, simpar, ishake )
    
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: pos ( :, : ), vel ( :, : )
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE(simulation_parameters_type), &
      INTENT(IN)                                      :: simpar
    INTEGER, INTENT(OUT)                              :: ishake

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ng3x3
    REAL(KIND=dp)                                :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, mass, &
                                                max_sigma, sigma
    REAL(KIND=dp), DIMENSION(3)                  :: fc1, fc2, fc3, r0_12, r0_13, &
                                                r0_23, r1, r12, r13, r2, r23, &
                                                r3, v1, v2, v3, vec
    REAL(KIND=dp), DIMENSION(3, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)               :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

  dtsqby2 = simpar%dt*simpar%dt*.5_dp
  idtsq = 1.0_dp/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
  DO iconst = 1, ng3x3
     ishake=0
     index_a = g3x3_list ( iconst ) % a + first_atom -1
     index_b = g3x3_list ( iconst ) % b + first_atom -1
     index_c = g3x3_list ( iconst ) % c + first_atom -1
     r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b ) 
     r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c ) 
     r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c ) 
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     lg3x3 ( iconst ) % fa = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
                                       lg3x3 ( iconst ) % rb_old )
     lg3x3 ( iconst ) % fb = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
                                       lg3x3 ( iconst ) % rc_old )
     lg3x3 ( iconst ) % fc = -2.0_dp*( lg3x3 ( iconst ) % rb_old - &
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
     max_sigma = -1.0E+10_dp
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0.0_dp
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
        r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*fc1 ( : )
        r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*fc2 ( : )
        r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*fc3 ( : )
        v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*fc1 ( : )
        v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*fc2 ( : )
        v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*fc3 ( : )
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
     pos ( :, index_a ) = r1 ( : )
     pos ( :, index_b ) = r2 ( : )
     pos ( :, index_c ) = r3 ( : )

! update velocites with full multiplier
     vel ( :, index_a ) = v1 ( : )
     vel ( :, index_b ) = v2 ( : )
     vel ( :, index_c ) = v3 ( : )
  END DO

END SUBROUTINE shake_3x3

!******************************************************************************

SUBROUTINE shake_roll_3x3( molecule, particle_set, pos, vel, r_shake,  &
                          v_shake, simpar, ishake )


    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: pos ( :, : ), vel ( :, : )
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: r_shake, v_shake
    INTEGER, INTENT(OUT)                     :: ishake

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ng3x3
    REAL(KIND=dp)                                :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, mass, &
                                                max_sigma, sigma
    REAL(KIND=dp), DIMENSION(3)                  :: fc1, fc2, fc3, r0_12, r0_13, &
                                                r0_23, r1, r12, r13, r2, r23, &
                                                r3, v1, v2, v3, vec
    REAL(KIND=dp), DIMENSION(3)                  :: f_roll1, f_roll2, f_roll3
    REAL(KIND=dp), DIMENSION(3, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)               :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

  dtsqby2 = simpar%dt*simpar%dt*.5_dp
  idtsq = 1.0_dp/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
  DO iconst = 1, ng3x3
     ishake=0
     index_a = g3x3_list ( iconst ) % a + first_atom -1
     index_b = g3x3_list ( iconst ) % b + first_atom -1
     index_c = g3x3_list ( iconst ) % c + first_atom -1
     r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b ) 
     r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c ) 
     r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c ) 
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     lg3x3 ( iconst ) % fa = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
                                       lg3x3 ( iconst ) % rb_old )
     lg3x3 ( iconst ) % fb = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
                                       lg3x3 ( iconst ) % rc_old )
     lg3x3 ( iconst ) % fc = -2.0_dp*( lg3x3 ( iconst ) % rb_old - &
                                       lg3x3 ( iconst ) % rc_old )
! rotate fconst:
     f_roll1 = matvec_3x3 ( r_shake, lg3x3 ( iconst ) % fa )
     f_roll2 = matvec_3x3 ( r_shake, lg3x3 ( iconst ) % fb )
     f_roll3 = matvec_3x3 ( r_shake, lg3x3 ( iconst ) % fc )

! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,f_roll1 )
     amat(1,2)=imass1*DOTPROD_3D(r0_12, f_roll2 )
     amat(1,3)=-imass2*DOTPROD_3D(r0_12, f_roll3 )
     amat(2,1)=imass1*DOTPROD_3D(r0_13, f_roll1 )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,f_roll2 )
     amat(2,3)=imass3*DOTPROD_3D(r0_13,f_roll3)
     amat(3,1)=-imass2*DOTPROD_3D(r0_23, f_roll1 )
     amat(3,2)=imass3*DOTPROD_3D(r0_23, f_roll2 )
     amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,f_roll3)

! Iterate until convergence:
     max_sigma = -1.0E+10_dp
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0.0_dp
        ishake=ishake+1
        vec= lg3x3 ( iconst ) % lambda ( 1 ) * f_roll1 *(imass1+imass2) + &
             lg3x3 ( iconst ) % lambda ( 2 ) * imass1* f_roll2 - &
             lg3x3 ( iconst ) % lambda ( 3 ) * imass2* f_roll3
        bvec(1,1)=g3x3_list ( iconst ) % dab * g3x3_list ( iconst) % dab &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
        vec= lg3x3 ( iconst ) % lambda ( 1 ) * f_roll1 * imass1 + &
             lg3x3 ( iconst ) % lambda ( 2 ) * ( imass1 + imass3 ) * f_roll2 + &
             lg3x3 ( iconst ) % lambda ( 3 ) * imass3 * f_roll3
        bvec(2,1)=g3x3_list ( iconst ) % dac * g3x3_list ( iconst) % dac &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
        vec= -lg3x3 ( iconst ) % lambda ( 1 ) * f_roll1 * imass2 + &
             lg3x3 ( iconst ) % lambda ( 2 ) * imass3 * f_roll2 + &
             lg3x3 ( iconst ) % lambda ( 3 ) * ( imass2 + imass3 ) * f_roll3
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
        r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*MATVEC_3x3 ( r_shake, fc1 )
        r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*MATVEC_3x3 ( r_shake, fc2 )
        r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*MATVEC_3x3 ( r_shake, fc3 )
        v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*MATVEC_3x3 ( v_shake, fc1 )
        v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*MATVEC_3x3 ( v_shake, fc2 )
        v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*MATVEC_3x3 ( v_shake, fc3 )
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
     pos ( :, index_a ) = r1 ( : )
     pos ( :, index_b ) = r2 ( : )
     pos ( :, index_c ) = r3 ( : )

! update velocites with full multiplier
     vel ( :, index_a ) = v1 ( : )
     vel ( :, index_b ) = v2 ( : )
     vel ( :, index_c ) = v3 ( : )
  END DO

END SUBROUTINE shake_roll_3x3
!
!!******************************************************************************

SUBROUTINE shake_4x6( molecule, particle_set, pos, vel, simpar, ishake )
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: pos ( :, : ), vel ( :, : )
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE(simulation_parameters_type), &
      INTENT(IN)                                      :: simpar
    INTEGER, INTENT(OUT)                              :: ishake

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c,  index_d
    INTEGER                                  :: first_atom, ng4x6
    REAL(KIND=dp)                                :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, &
                                                imass4, mass, max_sigma, sigma
    REAL(KIND=dp), DIMENSION(3) :: fc1, fc2, fc3, fc4, r0_12, r0_13, r0_14, &
      r0_23, r0_24, r0_34, r1, r12, r13, r14, r2, r23, r24, r3, r34, r4, v1, &
      v2, v3, v4, vec
    REAL(KIND=dp), DIMENSION(6, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)               :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g4x6_constraint_type ), POINTER   :: lg4x6 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g4x6_constraint_type ), POINTER   :: g4x6_list ( : )
!!------------------------------------------------------------------------------

  dtsqby2 = simpar%dt*simpar%dt*.5_dp
  idtsq = 1.0_dp/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
  DO iconst = 1, ng4x6
     ishake=0
     index_a = g4x6_list ( iconst ) % a + first_atom -1
     index_b = g4x6_list ( iconst ) % b + first_atom -1
     index_c = g4x6_list ( iconst ) % c + first_atom -1
     index_d = g4x6_list ( iconst ) % d + first_atom -1
     r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b ) 
     r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c ) 
     r0_14 ( : ) = pos ( :, index_a ) - pos ( :,index_d ) 
     r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c ) 
     r0_24 ( : ) = pos ( :, index_b ) - pos ( :,index_d ) 
     r0_34 ( : ) = pos ( :, index_c ) - pos ( :,index_d ) 
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     atomic_kind => particle_set ( index_d ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass4 = 1.0_dp/mass
     lg4x6 ( iconst ) % fa = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
                                       lg4x6 ( iconst ) % rb_old )
     lg4x6 ( iconst ) % fb = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
                                       lg4x6 ( iconst ) % rc_old )
     lg4x6 ( iconst ) % fc = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
                                       lg4x6 ( iconst ) % rd_old )
     lg4x6 ( iconst ) % fd = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
                                       lg4x6 ( iconst ) % rc_old )
     lg4x6 ( iconst ) % fe = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
                                       lg4x6 ( iconst ) % rd_old )
     lg4x6 ( iconst ) % ff = -2.0_dp*( lg4x6 ( iconst ) % rc_old - &
                                       lg4x6 ( iconst ) % rd_old )
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,lg4x6 ( iconst ) % fa )
     amat(1,2)=imass1*DOTPROD_3D(r0_12, lg4x6 ( iconst) % fb )
     amat(1,3)=imass1*DOTPROD_3D(r0_12, lg4x6 ( iconst ) % fc )
     amat(1,4)=-imass2*DOTPROD_3D(r0_12, lg4x6 ( iconst ) % fd )
     amat(1,5)=-imass2*DOTPROD_3D(r0_12, lg4x6 ( iconst ) % fe )
     amat(1,6)=0.0_dp

     amat(2,1)=imass1*DOTPROD_3D(r0_13,lg4x6 ( iconst ) % fa )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13, lg4x6 ( iconst) % fb )
     amat(2,3)=imass1*DOTPROD_3D(r0_13, lg4x6 ( iconst ) % fc )
     amat(2,4)=imass3*DOTPROD_3D(r0_13, lg4x6 ( iconst ) % fd )
     amat(2,5)=0.0_dp
     amat(2,6)=-imass3*DOTPROD_3D(r0_13, lg4x6 ( iconst ) % ff )

     amat(3,1)=imass1*DOTPROD_3D(r0_14,lg4x6 ( iconst ) % fa )
     amat(3,2)=imass1*DOTPROD_3D(r0_14, lg4x6 ( iconst) % fb )
     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r0_14, lg4x6 ( iconst ) % fc )
     amat(3,4)=0.0_dp
     amat(3,5)=imass4*DOTPROD_3D(r0_14, lg4x6 ( iconst ) % fe )
     amat(3,6)=imass4*DOTPROD_3D(r0_14, lg4x6 ( iconst ) % ff )

     amat(4,1)=-imass2*DOTPROD_3D(r0_23,lg4x6 ( iconst ) % fa )
     amat(4,2)=imass3*DOTPROD_3D(r0_23, lg4x6 ( iconst) % fb )
     amat(4,3)=0.0_dp
     amat(4,4)=(imass3+imass2)*DOTPROD_3D(r0_23, lg4x6 ( iconst ) % fd )
     amat(4,5)=imass2*DOTPROD_3D(r0_23, lg4x6 ( iconst ) % fe )
     amat(4,6)=-imass3*DOTPROD_3D(r0_23, lg4x6 ( iconst ) % ff )

     amat(5,1)=-imass2*DOTPROD_3D(r0_24,lg4x6 ( iconst ) % fa )
     amat(5,2)=0.0_dp
     amat(5,3)=imass4*DOTPROD_3D(r0_24, lg4x6 ( iconst) % fc )
     amat(5,4)=imass2*DOTPROD_3D(r0_24, lg4x6 ( iconst ) % fd )
     amat(5,5)=(imass4+imass2)*DOTPROD_3D(r0_24, lg4x6 ( iconst ) % fe )
     amat(5,6)=imass4*DOTPROD_3D(r0_24, lg4x6 ( iconst ) % ff )

     amat(6,1)=0.0_dp
     amat(6,2)=-imass3*DOTPROD_3D(r0_34,lg4x6 ( iconst ) % fb )
     amat(6,3)=imass4*DOTPROD_3D(r0_34, lg4x6 ( iconst) % fc )
     amat(6,4)=-imass3*DOTPROD_3D(r0_34, lg4x6 ( iconst ) % fd )
     amat(6,5)=imass4*DOTPROD_3D(r0_34, lg4x6 ( iconst ) % fe )
     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r0_34, lg4x6 ( iconst ) % ff )

! Iterate until convergence:
     max_sigma = -1.0E+10_dp
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0.0_dp
        ishake=ishake+1
        vec= lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa *(imass1+imass2) + &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass1*lg4x6 ( iconst ) % fb + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass1*lg4x6 ( iconst ) % fc - &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass2*lg4x6 ( iconst ) % fd - &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass2*lg4x6 ( iconst ) % fe
        bvec(1,1)=g4x6_list ( iconst ) % dab * g4x6_list ( iconst) % dab &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
        vec= lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb *(imass1+imass3) + &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass1*lg4x6 ( iconst ) % fa + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass1*lg4x6 ( iconst ) % fc + &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass3*lg4x6 ( iconst ) % fd - &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass3*lg4x6 ( iconst ) % ff
        bvec(2,1)=g4x6_list ( iconst ) % dac * g4x6_list ( iconst) % dac &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
        vec= lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc *(imass1+imass4) + &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass1*lg4x6 ( iconst ) % fa + &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass1*lg4x6 ( iconst ) % fb + &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass4*lg4x6 ( iconst ) % fe + &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass4*lg4x6 ( iconst ) % ff
        bvec(3,1)=g4x6_list ( iconst ) % dad * g4x6_list ( iconst) % dad &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_14,r0_14)
        vec= lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd *(imass2+imass3) - &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass2*lg4x6 ( iconst ) % fa + &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass3*lg4x6 ( iconst ) % fb + &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass2*lg4x6 ( iconst ) % fe - &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass3*lg4x6 ( iconst ) % ff
        bvec(4,1)=g4x6_list ( iconst ) % dbc * g4x6_list ( iconst) % dbc &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
        vec= lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe *(imass2+imass4) - &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass2*lg4x6 ( iconst ) % fa + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass4*lg4x6 ( iconst ) % fc + &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass2*lg4x6 ( iconst ) % fd + &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass4*lg4x6 ( iconst ) % ff
        bvec(5,1)=g4x6_list ( iconst ) % dbd * g4x6_list ( iconst) % dbd &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_24,r0_24)
        vec= lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff *(imass3+imass4) - &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass3*lg4x6 ( iconst ) % fb + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass4*lg4x6 ( iconst ) % fc - &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass3*lg4x6 ( iconst ) % fd + &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass4*lg4x6 ( iconst ) % fe
        bvec(6,1)=g4x6_list ( iconst ) % dcd * g4x6_list ( iconst) % dcd &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_34,r0_34)
        bvec = bvec * idtsq

! get lambda
        atemp = amat
        CALL solve_system ( atemp, 6, bvec )
        lg4x6 ( iconst ) % lambda ( : ) = bvec ( :, 1 )

        fc1=lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
            lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb + &
            lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc 
        fc2=-lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
            lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
            lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe 
        fc3=-lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb - &
            lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
            lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
        fc4=-lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc - &
            lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe - &
            lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
        r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*fc1 ( : )
        r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*fc2 ( : )
        r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*fc3 ( : )
        r4 ( : ) = pos ( :, index_d ) + imass4*dtsqby2*fc4 ( : )
        v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*fc1 ( : )
        v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*fc2 ( : )
        v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*fc3 ( : )
        v4 ( : ) = vel ( :, index_d ) + imass4*dtby2*fc4 ( : )
        r12=r1-r2
        r13=r1-r3
        r14=r1-r4
        r23=r2-r3
        r24=r2-r4
        r34=r3-r4
! compute the tolerance:
        sigma = dot_product(r12,r12) - g4x6_list ( iconst ) % dab *  &
                                       g4x6_list ( iconst ) % dab
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r13,r13) - g4x6_list ( iconst ) % dac *  &
                                       g4x6_list ( iconst ) % dac
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r14,r14) - g4x6_list ( iconst ) % dad *  &
                                       g4x6_list ( iconst ) % dad
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r23,r23) - g4x6_list ( iconst ) % dbc *  &
                                       g4x6_list ( iconst ) % dbc
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r24,r24) - g4x6_list ( iconst ) % dbd *  &
                                       g4x6_list ( iconst ) % dbd
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r34,r34) - g4x6_list ( iconst ) % dcd *  &
                                       g4x6_list ( iconst ) % dcd
        max_sigma=MAX(max_sigma,ABS(sigma))
     END DO CONS

! update positions with full multiplier
     pos ( :, index_a ) = r1 ( : )
     pos ( :, index_b ) = r2 ( : )
     pos ( :, index_c ) = r3 ( : )
     pos ( :, index_d ) = r4 ( : )

! update velocites with full multiplier
     vel ( :, index_a ) = v1 ( : )
     vel ( :, index_b ) = v2 ( : )
     vel ( :, index_c ) = v3 ( : )
     vel ( :, index_d ) = v4 ( : )
  END DO

END SUBROUTINE shake_4x6

!******************************************************************************
SUBROUTINE shake_roll_4x6( molecule, particle_set, pos, vel, r_shake,  &
                           v_shake, simpar, ishake )
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: pos ( :, : ), vel ( :, : )
    TYPE ( molecule_type ), POINTER                   :: molecule
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)            :: r_shake, v_shake
    TYPE(simulation_parameters_type), &
      INTENT(IN)                                      :: simpar
    INTEGER, INTENT(OUT)                              :: ishake

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c,  index_d
    INTEGER                                  :: first_atom, ng4x6
    REAL(KIND=dp)                            :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, &
                                                imass4, mass, max_sigma, sigma
    REAL(KIND=dp), DIMENSION(3) :: fc1, fc2, fc3, fc4, r0_12, r0_13, r0_14, &
      r0_23, r0_24, r0_34, r1, r12, r13, r14, r2, r23, r24, r3, r34, r4, v1, &
      v2, v3, v4, vec
    REAL(KIND=dp), DIMENSION(3) :: f_roll1, f_roll2, f_roll3, f_roll4, f_roll5, &
                               f_roll6
    REAL(KIND=dp), DIMENSION(6, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)               :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g4x6_constraint_type ), POINTER   :: lg4x6 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g4x6_constraint_type ), POINTER   :: g4x6_list ( : )
!!------------------------------------------------------------------------------

  dtsqby2 = simpar%dt*simpar%dt*.5_dp
  idtsq = 1.0_dp/simpar%dt/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
  DO iconst = 1, ng4x6
     ishake=0
     index_a = g4x6_list ( iconst ) % a + first_atom -1
     index_b = g4x6_list ( iconst ) % b + first_atom -1
     index_c = g4x6_list ( iconst ) % c + first_atom -1
     index_d = g4x6_list ( iconst ) % d + first_atom -1
     r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b ) 
     r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c ) 
     r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c ) 
     r0_14 ( : ) = pos ( :, index_a ) - pos ( :,index_d ) 
     r0_24 ( : ) = pos ( :, index_b ) - pos ( :,index_d ) 
     r0_34 ( : ) = pos ( :, index_c ) - pos ( :,index_d ) 
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     atomic_kind => particle_set ( index_d ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass4 = 1.0_dp/mass
     lg4x6 ( iconst ) % fa = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
                                       lg4x6 ( iconst ) % rb_old )
     lg4x6 ( iconst ) % fb = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
                                       lg4x6 ( iconst ) % rc_old )
     lg4x6 ( iconst ) % fc = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
                                       lg4x6 ( iconst ) % rd_old )
     lg4x6 ( iconst ) % fd = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
                                       lg4x6 ( iconst ) % rc_old )
     lg4x6 ( iconst ) % fe = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
                                       lg4x6 ( iconst ) % rd_old )
     lg4x6 ( iconst ) % ff = -2.0_dp*( lg4x6 ( iconst ) % rc_old - &
                                       lg4x6 ( iconst ) % rd_old )
! rotate fconst:
     f_roll1 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fa )
     f_roll2 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fb )
     f_roll3 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fc )
     f_roll4 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fd )
     f_roll5 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fe )
     f_roll6 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % ff )

! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12, f_roll1)
     amat(1,2)=imass1*DOTPROD_3D(r0_12, f_roll2 )
     amat(1,3)=imass1*DOTPROD_3D(r0_12, f_roll3 )
     amat(1,4)=-imass2*DOTPROD_3D(r0_12, f_roll4 )
     amat(1,5)=-imass2*DOTPROD_3D(r0_12, f_roll5 )
     amat(1,6)=0.0_dp

     amat(2,1)=imass1*DOTPROD_3D(r0_13, f_roll1 )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13, f_roll2 )
     amat(2,3)=imass1*DOTPROD_3D(r0_13, f_roll3 )
     amat(2,4)=imass3*DOTPROD_3D(r0_13, f_roll4 )
     amat(2,5)=0.0_dp
     amat(2,6)=-imass3*DOTPROD_3D(r0_13, f_roll6 )

     amat(3,1)=imass1*DOTPROD_3D(r0_14, f_roll1 )
     amat(3,2)=imass1*DOTPROD_3D(r0_14, f_roll2 )
     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r0_14, f_roll3 )
     amat(3,4)=0.0_dp
     amat(3,5)=imass4*DOTPROD_3D(r0_14, f_roll5 )
     amat(3,6)=imass4*DOTPROD_3D(r0_14, f_roll6 )

     amat(4,1)=-imass2*DOTPROD_3D(r0_23, f_roll1 )
     amat(4,2)=imass3*DOTPROD_3D(r0_23, f_roll2)
     amat(4,3)=0.0_dp
     amat(4,4)=(imass3+imass2)*DOTPROD_3D(r0_23, f_roll4 )
     amat(4,5)=imass2*DOTPROD_3D(r0_23, f_roll5 )
     amat(4,6)=-imass3*DOTPROD_3D(r0_23, f_roll6 )

     amat(5,1)=-imass2*DOTPROD_3D(r0_24, f_roll1 )
     amat(5,2)=0.0_dp
     amat(5,3)=imass4*DOTPROD_3D(r0_24, f_roll3 )
     amat(5,4)=imass2*DOTPROD_3D(r0_24, f_roll4 )
     amat(5,5)=(imass4+imass2)*DOTPROD_3D(r0_24, f_roll5 )
     amat(5,6)=imass4*DOTPROD_3D(r0_24, f_roll6 )

     amat(6,1)=0.0_dp
     amat(6,2)=-imass3*DOTPROD_3D(r0_34, f_roll2 )
     amat(6,3)=imass4*DOTPROD_3D(r0_34, f_roll3 )
     amat(6,4)=-imass3*DOTPROD_3D(r0_34, f_roll4 )
     amat(6,5)=imass4*DOTPROD_3D(r0_34, f_roll5 )
     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r0_34, f_roll6 )

! Iterate until convergence:
     max_sigma = -1.0E+10_dp
     CONS: DO WHILE (abs(max_sigma)>=simpar%shake_tol)
        max_sigma=0.0_dp
        ishake=ishake+1
        vec= lg4x6 ( iconst ) % lambda ( 1 ) * f_roll1*(imass1+imass2) + &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass1*f_roll2 + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass1*f_roll3 - &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass2*f_roll4 - &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass2*f_roll5
        bvec(1,1)=g4x6_list ( iconst ) % dab * g4x6_list ( iconst) % dab &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
        vec= lg4x6 ( iconst ) % lambda ( 2 ) * f_roll2 *(imass1+imass3) + &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass1*f_roll1 + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass1*f_roll3 + &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass3*f_roll4 - &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass3*f_roll6
        bvec(2,1)=g4x6_list ( iconst ) % dac * g4x6_list ( iconst) % dac &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
        vec= lg4x6 ( iconst ) % lambda ( 3 ) * f_roll3 *(imass1+imass4) + &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass1*f_roll1 + &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass1*f_roll2 + &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass4*f_roll5 + &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass4*f_roll6
        bvec(3,1)=g4x6_list ( iconst ) % dad * g4x6_list ( iconst) % dad &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_14,r0_14)
        vec= lg4x6 ( iconst ) % lambda ( 4 ) * f_roll4 *(imass2+imass3) - &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass2*f_roll1 + &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass3*f_roll2 + &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass2*f_roll5 - &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass3*f_roll6
        bvec(4,1)=g4x6_list ( iconst ) % dbc * g4x6_list ( iconst) % dbc &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
        vec= lg4x6 ( iconst ) % lambda ( 5 ) * f_roll5 *(imass2+imass4) - &
             lg4x6 ( iconst ) % lambda ( 1 ) * imass2*f_roll1 + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass4*f_roll3 + &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass2*f_roll4 + &
             lg4x6 ( iconst ) % lambda ( 6 ) * imass4*f_roll6
        bvec(5,1)=g4x6_list ( iconst ) % dbd * g4x6_list ( iconst) % dbd &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_24,r0_24)
        vec= lg4x6 ( iconst ) % lambda ( 6 ) * f_roll6 *(imass3+imass4) - &
             lg4x6 ( iconst ) % lambda ( 2 ) * imass3*f_roll2 + &
             lg4x6 ( iconst ) % lambda ( 3 ) * imass4*f_roll3 - &
             lg4x6 ( iconst ) % lambda ( 4 ) * imass3*f_roll4 + &
             lg4x6 ( iconst ) % lambda ( 5 ) * imass4*f_roll5
        bvec(6,1)=g4x6_list ( iconst ) % dcd * g4x6_list ( iconst) % dcd &
                 -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_34,r0_34)
        bvec = bvec * idtsq

! get lambda
        atemp = amat
        CALL solve_system ( atemp, 6, bvec )
        lg4x6 ( iconst ) % lambda ( : ) = bvec ( :, 1 )

        fc1=lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
            lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb + &
            lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc 
        fc2=-lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
            lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
            lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe 
        fc3=-lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb - &
            lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
            lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
        fc4=-lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc - &
            lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe - &
            lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
        r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*MATVEC_3x3(r_shake,fc1)
        r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*MATVEC_3x3(r_shake,fc2)
        r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*MATVEC_3x3(r_shake,fc3)
        r4 ( : ) = pos ( :, index_d ) + imass4*dtsqby2*MATVEC_3x3(r_shake,fc4)
        v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*MATVEC_3x3(r_shake,fc1)
        v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*MATVEC_3x3(r_shake,fc2)
        v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*MATVEC_3x3(r_shake,fc3)
        v4 ( : ) = vel ( :, index_d ) + imass4*dtby2*MATVEC_3x3(r_shake,fc4)
        r12=r1-r2
        r13=r1-r3
        r23=r2-r3
        r14=r1-r4
        r24=r2-r4
        r34=r3-r4
! compute the tolerance:
        sigma = dot_product(r12,r12) - g4x6_list ( iconst ) % dab *  &
                                       g4x6_list ( iconst ) % dab
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r13,r13) - g4x6_list ( iconst ) % dac *  &
                                       g4x6_list ( iconst ) % dac
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r14,r14) - g4x6_list ( iconst ) % dad *  &
                                       g4x6_list ( iconst ) % dad
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r23,r23) - g4x6_list ( iconst ) % dbc *  &
                                       g4x6_list ( iconst ) % dbc
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r24,r24) - g4x6_list ( iconst ) % dbd *  &
                                       g4x6_list ( iconst ) % dbd
        max_sigma=MAX(max_sigma,ABS(sigma))
        sigma = dot_product(r34,r34) - g4x6_list ( iconst ) % dcd *  &
                                       g4x6_list ( iconst ) % dcd
        max_sigma=MAX(max_sigma,ABS(sigma))
     END DO CONS

! update positions with full multiplier
     pos ( :, index_a ) = r1 ( : )
     pos ( :, index_b ) = r2 ( : )
     pos ( :, index_c ) = r3 ( : )
     pos ( :, index_d ) = r4 ( : )

! update velocites with full multiplier
     vel ( :, index_a ) = v1 ( : )
     vel ( :, index_b ) = v2 ( : )
     vel ( :, index_c ) = v3 ( : )
     vel ( :, index_d ) = v4 ( : )
  END DO

END SUBROUTINE shake_roll_4x6
!!******************************************************************************

SUBROUTINE rattle_roll_3x3( molecule, particle_set, vel, r_rattle, simpar, &
                            veps )
    IMPLICIT NONE

    TYPE ( molecule_type ), POINTER                   :: molecule
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : )
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
    INTENT(IN)                             :: simpar
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: r_rattle, veps

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ng3x3
    REAL(KIND=dp)                                :: dtby2, idt, imass1, imass2, &
                                                imass3, mass, sigma_dot
    REAL(KIND=dp), DIMENSION(3)                  :: f_roll1, f_roll2, f_roll3
    REAL(KIND=dp), DIMENSION(3)                  :: fc1, fc2, fc3, r12, r13, r23, &
                                                    v12, v13, v23, lambda
    REAL(KIND=dp), DIMENSION(3, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)               :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

  idt = 1.0_dp/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
  DO iconst = 1, SIZE ( g3x3_list )
     index_a = g3x3_list ( iconst ) % a + first_atom - 1
     index_b = g3x3_list ( iconst ) % b + first_atom - 1
     index_c = g3x3_list ( iconst ) % c + first_atom - 1
     v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b ) 
     v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c ) 
     v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c ) 
     r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
     r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
     r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
     atomic_kind=>particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind=>particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind=>particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     lg3x3 ( iconst ) % fa = -2.0_dp* r12
     lg3x3 ( iconst ) % fb = -2.0_dp* r13
     lg3x3 ( iconst ) % fc = -2.0_dp* r23

! roll the fc
     f_roll1 = MATVEC_3x3 ( r_rattle, lg3x3 ( iconst ) % fa )
     f_roll2 = MATVEC_3x3 ( r_rattle, lg3x3 ( iconst ) % fb )
     f_roll3 = MATVEC_3x3 ( r_rattle, lg3x3 ( iconst ) % fc )

! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12, f_roll1 )
     amat(1,2)=imass1*DOTPROD_3D(r12, f_roll2 )
     amat(1,3)=-imass2*DOTPROD_3D(r12,f_roll3 )
     amat(2,1)=imass1*DOTPROD_3D(r13, f_roll1 )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13, f_roll2 )
     amat(2,3)=imass3*DOTPROD_3D(r13, f_roll3 )
     amat(3,1)=-imass2*DOTPROD_3D(r23,f_roll1 )
     amat(3,2)=imass3*DOTPROD_3D(r23, f_roll2 )
     amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23, f_roll3 )

! construct solution vector
     bvec ( 1, 1 ) = DOTPROD_3D ( r12, v12 + MATVEC_3X3 ( veps, r12 ) )
     bvec ( 2, 1 ) = DOTPROD_3D ( r13, v13 + MATVEC_3X3 ( veps, r13 ) )
     bvec ( 3, 1 ) = DOTPROD_3D ( r23, v23 + MATVEC_3X3 ( veps, r23 ) )
     bvec = -bvec * 2.0_dp * idt

! get lambda
     CALL solve_system ( amat, 3, bvec )
     lambda ( : ) = bvec ( :, 1 )
     lg3x3 ( iconst ) % lambda ( : )= lambda 

     fc1= lambda ( 1 ) * f_roll1 + &
          lambda ( 2 ) * f_roll2
     fc2=-lambda ( 1 ) * f_roll1 + &
          lambda ( 3 ) * f_roll3
     fc3=-lambda ( 2 ) * f_roll2 - &
          lambda ( 3 ) * f_roll3
     vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
     vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
     vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
  END DO

END SUBROUTINE rattle_roll_3x3
!******************************************************************************

SUBROUTINE rattle_3x3( molecule, particle_set, vel, simpar)


    TYPE ( molecule_type ), POINTER                   :: molecule
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : )
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
    INTENT(IN)                             :: simpar

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ng3x3
    REAL(KIND=dp)                                :: dtby2, idt, imass1, imass2, &
                                                imass3, mass
    REAL(KIND=dp), DIMENSION(3)                  :: fc1, fc2, fc3, r12, r13, r23, &
                                                v12, v13, v23
    REAL(KIND=dp), DIMENSION(3, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)               :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )

!------------------------------------------------------------------------------
! define local variables

  idt = 1.0_dp/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
  DO iconst = 1, SIZE ( g3x3_list )
     index_a = g3x3_list ( iconst ) % a + first_atom - 1
     index_b = g3x3_list ( iconst ) % b + first_atom - 1
     index_c = g3x3_list ( iconst ) % c + first_atom - 1
     v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b ) 
     v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c ) 
     v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c ) 
     r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
     r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
     r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
     atomic_kind=>particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind=>particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind=>particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     lg3x3 ( iconst ) % fa = -2.0_dp* r12
     lg3x3 ( iconst ) % fb = -2.0_dp* r13
     lg3x3 ( iconst ) % fc = -2.0_dp* r23

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
     bvec=-bvec*2.0_dp*idt


! get lambda
     CALL solve_system ( amat, 3, bvec )
     lg3x3 ( iconst ) % lambda ( : )=bvec ( :, 1 )

     fc1=lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
         lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb
     fc2=-lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
         lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
     fc3=-lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb - &
         lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
     vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
     vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
     vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
  END DO

END SUBROUTINE rattle_3x3

!******************************************************************************
SUBROUTINE rattle_4x6( molecule, particle_set, vel, simpar)
    TYPE ( molecule_type ), POINTER                   :: molecule
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : )
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
    INTENT(IN)                             :: simpar

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c, index_d
    INTEGER                                  :: first_atom, ng4x6
    REAL(KIND=dp)                                :: dtby2, idt, imass1, imass2, &
                                                imass3, imass4, mass
    REAL(KIND=dp), DIMENSION(3)                  :: fc1, fc2, fc3, fc4, r12, r13, &
                                                r14, r23, r24, r34, v12, v13, &
                                                v14, v23, v24, v34
    REAL(KIND=dp), DIMENSION(6, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)               :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_g4x6_constraint_type ), POINTER   :: lg4x6 ( : )
    TYPE ( g4x6_constraint_type ), POINTER   :: g4x6_list ( : )

!------------------------------------------------------------------------------
! define local variables

  idt = 1.0_dp/simpar%dt
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
  CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
  DO iconst = 1, SIZE ( g4x6_list )
     index_a = g4x6_list ( iconst ) % a + first_atom -1
     index_b = g4x6_list ( iconst ) % b + first_atom -1
     index_c = g4x6_list ( iconst ) % c + first_atom -1
     index_d = g4x6_list ( iconst ) % d + first_atom -1
     v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b ) 
     v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c ) 
     v14 ( : ) = vel ( :, index_a ) - vel ( :, index_d ) 
     v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c ) 
     v24 ( : ) = vel ( :, index_b ) - vel ( :, index_d ) 
     v34 ( : ) = vel ( :, index_c ) - vel ( :, index_d ) 

     r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
     r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
     r14 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_d ) % r ( : )
     r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
     r24 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_d ) % r ( : )
     r34 ( : ) = particle_set ( index_c ) % r ( : ) - particle_set ( index_d ) % r ( : )
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     atomic_kind => particle_set ( index_d ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass4 = 1.0_dp/mass
     lg4x6 ( iconst ) % fa = -2.0_dp* r12
     lg4x6 ( iconst ) % fb = -2.0_dp* r13
     lg4x6 ( iconst ) % fc = -2.0_dp* r14
     lg4x6 ( iconst ) % fd = -2.0_dp* r23
     lg4x6 ( iconst ) % fe = -2.0_dp* r24
     lg4x6 ( iconst ) % ff = -2.0_dp* r34
!
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,lg4x6 ( iconst ) % fa )
     amat(1,2)=imass1*DOTPROD_3D(r12, lg4x6 ( iconst) % fb )
     amat(1,3)=imass1*DOTPROD_3D(r12, lg4x6 ( iconst ) % fc )
     amat(1,4)=-imass2*DOTPROD_3D(r12, lg4x6 ( iconst ) % fd )
     amat(1,5)=-imass2*DOTPROD_3D(r12, lg4x6 ( iconst ) % fe )
     amat(1,6)=0.0_dp

     amat(2,1)=imass1*DOTPROD_3D(r13,lg4x6 ( iconst ) % fa )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13, lg4x6 ( iconst) % fb )
     amat(2,3)=imass1*DOTPROD_3D(r13, lg4x6 ( iconst ) % fc )
     amat(2,4)=imass3*DOTPROD_3D(r13, lg4x6 ( iconst ) % fd )
     amat(2,5)=0.0_dp
     amat(2,6)=-imass3*DOTPROD_3D(r13, lg4x6 ( iconst ) % ff )

     amat(3,1)=imass1*DOTPROD_3D(r14,lg4x6 ( iconst ) % fa )
     amat(3,2)=imass1*DOTPROD_3D(r14, lg4x6 ( iconst) % fb )
     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r14, lg4x6 ( iconst ) % fc )
     amat(3,4)=0.0_dp
     amat(3,5)=imass4*DOTPROD_3D(r14, lg4x6 ( iconst ) % fe )
     amat(3,6)=imass4*DOTPROD_3D(r14, lg4x6 ( iconst ) % ff )

     amat(4,1)=-imass2*DOTPROD_3D(r23,lg4x6 ( iconst ) % fa )
     amat(4,2)=imass3*DOTPROD_3D(r23, lg4x6 ( iconst) % fb )
     amat(4,3)=0.0_dp
     amat(4,4)=(imass3+imass2)*DOTPROD_3D(r23, lg4x6 ( iconst ) % fd )
     amat(4,5)=imass2*DOTPROD_3D(r23, lg4x6 ( iconst ) % fe )
     amat(4,6)=-imass3*DOTPROD_3D(r23, lg4x6 ( iconst ) % ff )

     amat(5,1)=-imass2*DOTPROD_3D(r24,lg4x6 ( iconst ) % fa )
     amat(5,2)=0.0_dp
     amat(5,3)=imass4*DOTPROD_3D(r24, lg4x6 ( iconst) % fc )
     amat(5,4)=imass2*DOTPROD_3D(r24, lg4x6 ( iconst ) % fd )
     amat(5,5)=(imass4+imass2)*DOTPROD_3D(r24, lg4x6 ( iconst ) % fe )
     amat(5,6)=imass4*DOTPROD_3D(r24, lg4x6 ( iconst ) % ff )

     amat(6,1)=0.0_dp
     amat(6,2)=-imass3*DOTPROD_3D(r34,lg4x6 ( iconst ) % fb )
     amat(6,3)=imass4*DOTPROD_3D(r34, lg4x6 ( iconst) % fc )
     amat(6,4)=-imass3*DOTPROD_3D(r34, lg4x6 ( iconst ) % fd )
     amat(6,5)=imass4*DOTPROD_3D(r34, lg4x6 ( iconst ) % fe )
     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r34, lg4x6 ( iconst ) % ff )

! construct solution vector
     bvec(1,1)=DOTPROD_3D(r12,v12)
     bvec(2,1)=DOTPROD_3D(r13,v13)
     bvec(3,1)=DOTPROD_3D(r14,v14)
     bvec(4,1)=DOTPROD_3D(r23,v23)
     bvec(5,1)=DOTPROD_3D(r24,v24)
     bvec(6,1)=DOTPROD_3D(r34,v34)
     bvec=-bvec*2.0_dp*idt

! get lambda
     CALL solve_system ( amat, 6, bvec )
     lg4x6 ( iconst ) % lambda ( : )=bvec ( :, 1 )

     fc1=lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
         lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb + &
         lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc 
     fc2=-lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
         lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
         lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe 
     fc3=-lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb - &
         lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
         lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
     fc4=-lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc - &
         lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe - &
         lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
     vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
     vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
     vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
     vel ( :, index_d ) = vel ( :, index_d ) + imass4 * dtby2 * fc4 ( : )
  END DO

END SUBROUTINE rattle_4x6
!!******************************************************************************
!
SUBROUTINE rattle_roll_4x6( molecule, particle_set, vel, r_rattle, simpar, &
                            veps )
    IMPLICIT NONE
    TYPE ( molecule_type ), POINTER                   :: molecule
    REAL (KIND=dp), INTENT ( INOUT )                  :: vel ( :, : )
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    TYPE(simulation_parameters_type), &
    INTENT(IN)                             :: simpar
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: r_rattle, veps

! locals
    INTEGER                                  :: iconst, index_a, index_b, index_c, index_d
    INTEGER                                  :: first_atom, ng4x6
    REAL(KIND=dp)                                :: dtby2, idt, imass1, imass2, &
                                                imass3, imass4, mass, sigma_dot
    REAL(KIND=dp), DIMENSION(3)                  :: fc1, fc2, fc3, fc4, r12, r13, &
                                                r14, r23, r24, r34, v12, v13, &
                                                v14, v23, v24, v34
    REAL(KIND=dp), DIMENSION(6, 1)               :: bvec
    REAL(KIND=dp), DIMENSION(3) :: f_roll1, f_roll2, f_roll3, f_roll4, f_roll5, &
                                   f_roll6
    REAL(KIND=dp), DIMENSION(6)                  :: lambda
    REAL(KIND=dp), DIMENSION(6, 6)               :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_g4x6_constraint_type ), POINTER   :: lg4x6 ( : )
    TYPE ( g4x6_constraint_type ), POINTER   :: g4x6_list ( : )
!------------------------------------------------------------------------------
! define local variables
    idt = 1.0_dp/simpar%dt
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
    DO iconst = 1, SIZE ( g4x6_list )
     index_a = g4x6_list ( iconst ) % a + first_atom -1
     index_b = g4x6_list ( iconst ) % b + first_atom -1
     index_c = g4x6_list ( iconst ) % c + first_atom -1
     index_d = g4x6_list ( iconst ) % d + first_atom -1
     v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b ) 
     v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c ) 
     v14 ( : ) = vel ( :, index_a ) - vel ( :, index_d ) 
     v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c ) 
     v24 ( : ) = vel ( :, index_b ) - vel ( :, index_d ) 
     v34 ( : ) = vel ( :, index_c ) - vel ( :, index_d ) 

     r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
     r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
     r14 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_d ) % r ( : )
     r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
     r24 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_d ) % r ( : )
     r34 ( : ) = particle_set ( index_c ) % r ( : ) - particle_set ( index_d ) % r ( : )
     atomic_kind => particle_set ( index_a ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass1 = 1.0_dp/mass
     atomic_kind => particle_set ( index_b ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass2 = 1.0_dp/mass
     atomic_kind => particle_set ( index_c ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass3 = 1.0_dp/mass
     atomic_kind => particle_set ( index_d ) % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
     imass4 = 1.0_dp/mass
     lg4x6 ( iconst ) % fa = -2.0_dp* r12
     lg4x6 ( iconst ) % fb = -2.0_dp* r13
     lg4x6 ( iconst ) % fc = -2.0_dp* r14
     lg4x6 ( iconst ) % fd = -2.0_dp* r23
     lg4x6 ( iconst ) % fe = -2.0_dp* r24
     lg4x6 ( iconst ) % ff = -2.0_dp* r34
! roll the fc
     f_roll1 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fa )
     f_roll2 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fb )
     f_roll3 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fc )
     f_roll4 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fd )
     f_roll5 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fe )
     f_roll6 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % ff )
! construct matrix
     amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,f_roll1 )
     amat(1,2)=imass1*DOTPROD_3D(r12, f_roll2 )
     amat(1,3)=imass1*DOTPROD_3D(r12, f_roll3 )
     amat(1,4)=-imass2*DOTPROD_3D(r12, f_roll4 )
     amat(1,5)=-imass2*DOTPROD_3D(r12, f_roll5 )
     amat(1,6)=0.0_dp

     amat(2,1)=imass1*DOTPROD_3D(r13,f_roll1 )
     amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13, f_roll2)
     amat(2,3)=imass1*DOTPROD_3D(r13, f_roll3 )
     amat(2,4)=imass3*DOTPROD_3D(r13, f_roll4 )
     amat(2,5)=0.0_dp
     amat(2,6)=-imass3*DOTPROD_3D(r13, f_roll6 )

     amat(3,1)=imass1*DOTPROD_3D(r14,f_roll1 )
     amat(3,2)=imass1*DOTPROD_3D(r14, f_roll2 )
     amat(3,3)=(imass1+imass4)*DOTPROD_3D(r14, f_roll3 )
     amat(3,4)=0.0_dp
     amat(3,5)=imass4*DOTPROD_3D(r14, f_roll5 )
     amat(3,6)=imass4*DOTPROD_3D(r14, f_roll6 )

     amat(4,1)=-imass2*DOTPROD_3D(r23,f_roll1 )
     amat(4,2)=imass3*DOTPROD_3D(r23, f_roll2 )
     amat(4,3)=0.0_dp
     amat(4,4)=(imass3+imass2)*DOTPROD_3D(r23, f_roll4 )
     amat(4,5)=imass2*DOTPROD_3D(r23, f_roll5 )
     amat(4,6)=-imass3*DOTPROD_3D(r23, f_roll6 )

     amat(5,1)=-imass2*DOTPROD_3D(r24,f_roll1 )
     amat(5,2)=0.0_dp
     amat(5,3)=imass4*DOTPROD_3D(r24, f_roll3)
     amat(5,4)=imass2*DOTPROD_3D(r24, f_roll4 )
     amat(5,5)=(imass4+imass2)*DOTPROD_3D(r24, f_roll5 )
     amat(5,6)=imass4*DOTPROD_3D(r24, f_roll6 )

     amat(6,1)=0.0_dp
     amat(6,2)=-imass3*DOTPROD_3D(r34,f_roll2 )
     amat(6,3)=imass4*DOTPROD_3D(r34, f_roll3 )
     amat(6,4)=-imass3*DOTPROD_3D(r34, f_roll4 )
     amat(6,5)=imass4*DOTPROD_3D(r34, f_roll5 )
     amat(6,6)=(imass3+imass4)*DOTPROD_3D(r34, f_roll6 )

! construct solution vector
     bvec ( 1, 1 ) = DOTPROD_3D ( r12, v12 + MATVEC_3X3 ( veps, r12 ) )
     bvec ( 2, 1 ) = DOTPROD_3D ( r13, v13 + MATVEC_3X3 ( veps, r13 ) )
     bvec ( 3, 1 ) = DOTPROD_3D ( r14, v14 + MATVEC_3X3 ( veps, r14 ) )
     bvec ( 4, 1 ) = DOTPROD_3D ( r23, v23 + MATVEC_3X3 ( veps, r23 ) )
     bvec ( 5, 1 ) = DOTPROD_3D ( r24, v24 + MATVEC_3X3 ( veps, r24 ) )
     bvec ( 6, 1 ) = DOTPROD_3D ( r34, v34 + MATVEC_3X3 ( veps, r34 ) )
     bvec=-bvec*2.0_dp*idt


! get lambda
     CALL solve_system ( amat, 6, bvec )
     lambda ( : ) = bvec ( :, 1 )
     lg4x6 ( iconst ) % lambda ( : )= lambda 

     fc1=lambda ( 1 ) * f_roll1 + &
         lambda ( 2 ) * f_roll2 + &
         lambda ( 3 ) * f_roll3 
     fc2=-lambda ( 1 ) * f_roll1 + &
          lambda ( 4 ) * f_roll4 + &
          lambda ( 5 ) * f_roll5 
     fc3=-lambda ( 2 ) * f_roll2 - &
          lambda ( 4 ) * f_roll4 + &
          lambda ( 6 ) * f_roll6 
     fc4=-lambda ( 3 ) * f_roll3 - &
          lambda ( 5 ) * f_roll5 - &
          lambda ( 6 ) * f_roll6 
     vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
     vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
     vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
     vel ( :, index_d ) = vel ( :, index_d ) + imass4 * dtby2 * fc4 ( : )
  END DO

END SUBROUTINE rattle_roll_4x6

!******************************************************************************
SUBROUTINE shake( molecule, particle_set, pos, vel, simpar, max_sigma, ishake )
! updates the multiplier one molecule type at a time
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : ), pos ( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    INTEGER, INTENT(INOUT)                     :: ishake
    REAL(KIND=dp), INTENT(INOUT)                 :: max_sigma

! locals

    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ndist
    REAL(KIND=dp)                                :: del_lam, dtby2, dtsqby2, &
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(KIND=dp), DIMENSION(3)                  :: fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )
!
!------------------------------------------------------------------------------
! define local variables

  dtsqby2 = simpar%dt*simpar%dt*.5_dp
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
  CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
  IF (ishake==1) THEN
     DO iconst = 1, ndist
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
        ldist ( iconst ) % f = -2.0_dp * ( ldist ( iconst ) % ra_old - &
                                           ldist ( iconst ) % rb_old )

! update positions with full multiplier
         pos ( :, index_a ) = pos ( :, index_a ) + &
             imass1*dtsqby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
         pos ( :, index_b ) = pos ( :, index_b ) - &
             imass2*dtsqby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )

! update velocites with full multiplier
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : ) 
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : ) 
     END DO

  ELSE

     DO iconst = 1, ndist
! computing the constraint
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        r12 = pos ( :, index_a )  -  pos ( :, index_b )
        ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
                                                          dist_list ( iconst ) % dab
! computing fc(dt)
        fc_t = -2.0_dp*r12
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
        fdotf = dot_product( ldist ( iconst ) % f , fc_t )
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
        ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

! update positions by delta_lam
         pos ( :, index_a ) = pos ( :, index_a ) + &
             imass1*dtsqby2*del_lam * ldist ( iconst ) % f ( : )
         pos ( :, index_b ) = pos ( :, index_b ) - &
             imass2*dtsqby2*del_lam * ldist ( iconst ) % f ( : )

! update velocites by delta_lam
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*del_lam * ldist ( iconst ) % f ( : ) 
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*del_lam * ldist ( iconst ) % f ( : ) 
     END DO
  END IF

! computing the constraint and value of tolerance
  DO iconst = 1, ndist
     index_a = dist_list ( iconst ) % a + first_atom -1
     index_b = dist_list ( iconst ) % b + first_atom -1
     r12 ( : ) = pos ( :, index_a ) -  pos ( :, index_b )
     ldist ( iconst ) % sigma = dot_product(r12,r12) - dist_list ( iconst ) % dab * &
                                                       dist_list ( iconst ) % dab
     max_sigma = MAX(ABS(ldist ( iconst ) % sigma),max_sigma)
  END DO

END SUBROUTINE shake

!******************************************************************************
SUBROUTINE rattle ( molecule, particle_set, vel, simpar, max_sigma, irattle )
! updates the multiplier one molecule type at a time
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    INTEGER, INTENT(INOUT)                   :: irattle
    REAL(KIND=dp), INTENT(INOUT)                 :: max_sigma

    INTEGER                                  :: iconst, index_a, index_b
    INTEGER                                  :: first_atom, ndist
    REAL(KIND=dp)                                :: del_lam, dtby2, fdotf, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(KIND=dp), DIMENSION(3)                  :: r12, v12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )

!------------------------------------------------------------------------------
! define local variables

   dtby2 = simpar%dt*.5_dp
   molecule_kind => molecule % molecule_kind
   CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
   CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
  IF (irattle==1) THEN
     DO iconst = 1, ndist
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
        ldist ( iconst ) % f = -2.0_dp * ( particle_set ( index_a ) % r - &
                                           particle_set ( index_b ) % r )
! update velocites with full multiplier
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
     END DO

  ELSE

     DO iconst = 1, ndist

! computing sigma
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        r12 ( : ) = particle_set ( index_a ) % r ( : ) -  particle_set ( index_b ) % r ( : )
        v12 ( : ) = vel ( :, index_a ) -  vel ( :, index_b ) 
        ldist ( iconst ) % sigma = 2.0_dp*dot_product(r12,v12)

! computing f(dt)
        ldist ( iconst ) % f = -2.0_dp * r12
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
        fdotf = DOT_PRODUCT( ldist ( iconst ) % f , ldist ( iconst ) % f )
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*fdotf_sum)
        ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

!update velocites with del_lam
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*del_lam * ldist ( iconst ) % f ( : )
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*del_lam * ldist ( iconst ) % f ( : )
     END DO
  END IF

! computing the constraint and value of the tolerance
  DO iconst = 1, ndist
     index_a = dist_list ( iconst ) % a + first_atom -1
     index_b = dist_list ( iconst ) % b + first_atom -1
     r12 = particle_set ( index_a ) % r -  particle_set ( index_b ) % r
     v12 = vel ( :, index_a ) -  vel ( :, index_b )
     ldist ( iconst ) % sigma = dot_product(r12,v12)
     max_sigma = MAX ( ABS ( ldist ( iconst ) % sigma ), max_sigma )
  END DO

END SUBROUTINE rattle

!******************************************************************************
SUBROUTINE rattle_roll ( molecule, particle_set, vel, r_rattle,  &
                         simpar, veps, max_sigma, irattle )
! updates the multiplier one molecule type at a time
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : )
    REAL(KIND=dp), INTENT ( IN )                          :: veps ( :, : )
    REAL(KIND=dp), INTENT ( IN )                          :: r_rattle ( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    INTEGER, INTENT(INOUT)                   :: irattle
    REAL(KIND=dp), INTENT(INOUT)                 :: max_sigma

    INTEGER                                  :: iconst, index_a, index_b
    INTEGER                                  :: first_atom, ndist
    REAL(KIND=dp)                            :: dtby2, fdotf, del_lam, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(KIND=dp), DIMENSION(3)                  :: r12, v12, r1, r2, f_roll
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )

!------------------------------------------------------------------------------
! define local variables

   dtby2 = simpar%dt*.5_dp
   molecule_kind => molecule % molecule_kind
   CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
   CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
   IF (irattle==1) THEN
     DO iconst = 1, SIZE ( dist_list )
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
        ldist ( iconst ) % f = -2.0_dp * ( particle_set ( index_a ) % r - &
                                           particle_set ( index_b ) % r )
! update velocites with full multiplier
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
     END DO
     max_sigma = -1.E+10_dp
  ELSE
     DO iconst = 1, SIZE ( dist_list )
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        r1 ( : ) = particle_set ( index_a ) % r ( : ) 
        r2 ( : ) = particle_set ( index_b ) % r ( : ) 
        r12 ( : ) = r1 ( : ) - r2 ( : )
        v12 ( : ) = vel ( :, index_a ) -  vel ( :, index_b ) 

!  rolling f(dt)
        ldist ( iconst ) % f = -2.0_dp * r12
        f_roll = MATVEC_3x3 ( r_rattle, ldist ( iconst ) % f )

! computing sigma
        v12 ( : ) = vel ( :, index_a ) + MATVEC_3X3 ( veps, r1 ( : ) ) &
             - vel ( :, index_b ) - MATVEC_3X3 ( veps, r2 ( : ) )
        ldist ( iconst ) % sigma = 2.0_dp * DOTPROD_3D ( r12, v12 )

        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass

        fdotf = DOT_PRODUCT( f_roll , ldist ( iconst ) % f )
        fdotf_sum = fdotf*imass1 + fdotf*imass2
        del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*fdotf_sum)
        ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

!update velocites with del_lam
         vel ( :, index_a ) = vel ( :, index_a ) + &
             Imass1*dtby2*del_lam * f_roll ( : )
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*del_lam * f_roll ( : )
     END DO
! computing the constraint and value of the tolerance
     DO iconst = 1, SIZE ( dist_list )
       index_a = dist_list ( iconst ) % a + first_atom -1
       index_b = dist_list ( iconst ) % b + first_atom -1
       r1 ( : ) = particle_set ( index_a ) % r ( : ) 
       r2 ( : ) = particle_set ( index_b ) % r ( : ) 
       r12 ( : ) = r1 ( : ) - r2 ( : )
       v12 ( : ) = vel ( :, index_a ) + MATVEC_3X3 ( veps, r1 ( : ) ) &
                 - vel ( :, index_b ) - MATVEC_3X3 ( veps, r2 ( : ) )
       ldist ( iconst ) % sigma = 2.0_dp * DOTPROD_3D ( r12, v12 )
       max_sigma = MAX ( ABS ( ldist ( iconst ) % sigma ), max_sigma )
     END DO
  END IF
END SUBROUTINE rattle_roll
!!******************************************************************************

SUBROUTINE getold( local_molecules, molecule_set, molecule_kind_set, particle_set )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set ( : )
! locals
    INTEGER                                  :: iconst, ikind, nkind, i, first_atom
    INTEGER                                  :: imol, n3x3con, n4x6con, ndcon, index
    INTEGER                                  :: nmol_per_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( local_g4x6_constraint_type ), POINTER   :: lg4x6 ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
    TYPE ( g3x3_constraint_type ), POINTER   :: g3x3_list ( : )
    TYPE ( g4x6_constraint_type ), POINTER   :: g4x6_list ( : )
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
           CALL get_molecule_kind ( molecule_kind, ndist=ndcon, ng3x3=n3x3con, ng4x6=n4x6con, & 
                                    dist_list=dist_list, g3x3_list=g3x3_list, g4x6_list=g4x6_list )
           CALL get_molecule ( molecule, first_atom=first_atom, &
                               ldist=ldist, lg3x3=lg3x3, lg4x6=lg4x6 )
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
           DO iconst = 1, n4x6con
             index = g4x6_list ( iconst ) % a + first_atom - 1
             lg4x6 ( iconst ) % ra_old = particle_set ( index ) % r 
             index = g4x6_list ( iconst ) % b + first_atom - 1
             lg4x6 ( iconst ) % rb_old = particle_set ( index ) % r 
             index = g4x6_list ( iconst ) % c + first_atom - 1
             lg4x6 ( iconst ) % rc_old = particle_set ( index ) % r 
             index = g4x6_list ( iconst ) % d + first_atom - 1
             lg4x6 ( iconst ) % rd_old = particle_set ( index ) % r 
           ENDDO
         END DO
       END DO

END SUBROUTINE getold

!******************************************************************************

SUBROUTINE pv_constraint ( local_molecules, molecule_set, molecule_kind_set, &
                           particle_set, virial, group )

    IMPLICIT NONE
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE ( particle_type ), POINTER          :: particle_set ( : )
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: nkind, iconst, imol, ikind, nmol_per_kind, i
    INTEGER                                  :: ndist, ng3x3, ng4x6
    INTEGER                                  :: index_a, index_b, index_c, index_d
    INTEGER                                  :: first_atom
    REAL(KIND=dp) :: fc1( 3 ), fc2( 3 ), fc3( 3 ), fc4( 3 ), lambda, &
      lambda_3x3( 3 ), lambda_4x6( 6 ), pv( 3, 3 )
    TYPE ( local_g3x3_constraint_type ), POINTER   :: lg3x3 ( : )
    TYPE ( local_g4x6_constraint_type ), POINTER   :: lg4x6 ( : )
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( g3x3_constraint_type )      , POINTER   :: g3x3_list ( : )
    TYPE ( g4x6_constraint_type )      , POINTER   :: g4x6_list ( : )
    TYPE ( dist_constraint_type )      , POINTER   :: dist_list ( : )
    TYPE ( molecule_type ), POINTER          :: molecule
    TYPE ( molecule_kind_type ), POINTER     :: molecule_kind
!
!------------------------------------------------------------------------------
!
   pv = 0.0_dp
   nkind = SIZE ( molecule_kind_set )
   MOL: DO ikind = 1, nkind
         nmol_per_kind = local_molecules % n_el ( ikind )
         DO imol = 1, nmol_per_kind
           i = local_molecules % list ( ikind ) % array ( imol )
           molecule => molecule_set ( i ) 
           molecule_kind => molecule % molecule_kind
           CALL get_molecule_kind ( molecule_kind, ndist = ndist, ng3x3 = ng3x3, &
                                    ng4x6 = ng4x6, g3x3_list=g3x3_list, g4x6_list=g4x6_list, &
                                    dist_list = dist_list )
           CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3, &
                               lg4x6 = lg4x6, ldist = ldist )
! dist
           DO iconst = 1, ndist
!  pv gets updated with FULL multiplier
            lambda = ldist ( iconst ) % lambda

            index_a = dist_list ( iconst ) % a + first_atom -1
            index_b = dist_list ( iconst ) % b + first_atom -1
!pv(1,1)
            pv(1,1) = pv(1,1) + lambda*ldist ( iconst ) %f(1)*particle_set ( index_a ) % r ( 1 )
            pv(1,1) = pv(1,1) - lambda*ldist ( iconst ) %f(1)*particle_set ( index_b ) % r ( 1 )
!pv(1,2)
            pv(1,2) = pv(1,2) + lambda*ldist ( iconst ) %f(1)*particle_set ( index_a ) % r ( 2 )
            pv(1,2) = pv(1,2) - lambda*ldist ( iconst ) %f(1)*particle_set ( index_b ) % r ( 2 )
!pv(1,3)
            pv(1,3) = pv(1,3) + lambda*ldist ( iconst ) %f(1)*particle_set ( index_a ) % r ( 3 )
            pv(1,3) = pv(1,3) - lambda*ldist ( iconst ) %f(1)*particle_set ( index_b ) % r ( 3 )
!pv(2,1)
            pv(2,1) = pv(2,1) + lambda*ldist ( iconst ) %f(2)*particle_set ( index_a ) % r ( 1 )
            pv(2,1) = pv(2,1) - lambda*ldist ( iconst ) %f(2)*particle_set ( index_b ) % r ( 1 )
!pv(2,2)
            pv(2,2) = pv(2,2) + lambda*ldist ( iconst ) %f(2)*particle_set ( index_a ) % r ( 2 )
            pv(2,2) = pv(2,2) - lambda*ldist ( iconst ) %f(2)*particle_set ( index_b ) % r ( 2 )
!pv(2,3)
            pv(2,3) = pv(2,3) + lambda*ldist ( iconst ) %f(2)*particle_set ( index_a ) % r ( 3 )
            pv(2,3) = pv(2,3) - lambda*ldist ( iconst ) %f(2)*particle_set ( index_b ) % r ( 3 )
!pv(3,1)
            pv(3,1) = pv(3,1) + lambda*ldist ( iconst ) %f(3)*particle_set ( index_a ) % r ( 1 )
            pv(3,1) = pv(3,1) - lambda*ldist ( iconst ) %f(3)*particle_set ( index_b ) % r ( 1 )
!pv(3,2)
            pv(3,2) = pv(3,2) + lambda*ldist ( iconst ) %f(3)*particle_set ( index_a ) % r ( 2 )
            pv(3,2) = pv(3,2) - lambda*ldist ( iconst ) %f(3)*particle_set ( index_b ) % r ( 2 )
!pv(3,3)
            pv(3,3) = pv(3,3) + lambda*ldist ( iconst ) %f(3)*particle_set ( index_a ) % r ( 3 )
            pv(3,3) = pv(3,3) - lambda*ldist ( iconst ) %f(3)*particle_set ( index_b ) % r ( 3 )
          END DO
! 3x3
          DO iconst = 1, ng3x3

!  pv gets updated with FULL multiplier
            lambda_3x3 = lg3x3 ( iconst ) % lambda


             fc1=lambda_3x3(1)*lg3x3 ( iconst ) %fa+ &
                  lambda_3x3(2)*lg3x3 ( iconst ) %fb
             fc2=-lambda_3x3(1)*lg3x3 ( iconst ) %fa+ &
                  lambda_3x3(3)*lg3x3 ( iconst ) %fc
             fc3=-lambda_3x3(2)*lg3x3 ( iconst ) %fb- &
                  lambda_3x3(3)*lg3x3 ( iconst ) %fc
             index_a = g3x3_list ( iconst ) % a + first_atom -1
             index_b = g3x3_list ( iconst ) % b + first_atom -1
             index_c = g3x3_list ( iconst ) % c + first_atom -1

!pv(1,1)
             pv(1,1) = pv(1,1) + fc1(1)*particle_set ( index_a ) % r ( 1 )
             pv(1,1) = pv(1,1) + fc2(1)*particle_set ( index_b ) % r ( 1 )
             pv(1,1) = pv(1,1) + fc3(1)*particle_set ( index_c ) % r ( 1 )
!pv(1,2)
             pv(1,2) = pv(1,2) + fc1(1)*particle_set ( index_a ) % r ( 2 )
             pv(1,2) = pv(1,2) + fc2(1)*particle_set ( index_b ) % r ( 2 )
             pv(1,2) = pv(1,2) + fc3(1)*particle_set ( index_c ) % r ( 2 )
!pv(1,3)
             pv(1,3) = pv(1,3) + fc1(1)*particle_set ( index_a ) % r ( 3 )
             pv(1,3) = pv(1,3) + fc2(1)*particle_set ( index_b ) % r ( 3 )
             pv(1,3) = pv(1,3) + fc3(1)*particle_set ( index_c ) % r ( 3 )
!pv(2,1)
             pv(2,1) = pv(2,1) + fc1(2)*particle_set ( index_a ) % r ( 1 )
             pv(2,1) = pv(2,1) + fc2(2)*particle_set ( index_b ) % r ( 1 )
             pv(2,1) = pv(2,1) + fc3(2)*particle_set ( index_c ) % r ( 1 )
!pv(2,2)
             pv(2,2) = pv(2,2) + fc1(2)*particle_set ( index_a ) % r ( 2 )
             pv(2,2) = pv(2,2) + fc2(2)*particle_set ( index_b ) % r ( 2 )
             pv(2,2) = pv(2,2) + fc3(2)*particle_set ( index_c ) % r ( 2 )
!pv(2,3)
             pv(2,3) = pv(2,3) + fc1(2)*particle_set ( index_a ) % r ( 3 )
             pv(2,3) = pv(2,3) + fc2(2)*particle_set ( index_b ) % r ( 3 )
             pv(2,3) = pv(2,3) + fc3(2)*particle_set ( index_c ) % r ( 3 )
!pv(3,1)
             pv(3,1) = pv(3,1) + fc1(3)*particle_set ( index_a ) % r ( 1 )
             pv(3,1) = pv(3,1) + fc2(3)*particle_set ( index_b ) % r ( 1 )
             pv(3,1) = pv(3,1) + fc3(3)*particle_set ( index_c ) % r ( 1 )
!pv(3,2)
             pv(3,2) = pv(3,2) + fc1(3)*particle_set ( index_a ) % r ( 2 )
             pv(3,2) = pv(3,2) + fc2(3)*particle_set ( index_b ) % r ( 2 )
             pv(3,2) = pv(3,2) + fc3(3)*particle_set ( index_c ) % r ( 2 )
!pv(3,3)
             pv(3,3) = pv(3,3) + fc1(3)*particle_set ( index_a ) % r ( 3 )
             pv(3,3) = pv(3,3) + fc2(3)*particle_set ( index_b ) % r ( 3 )
             pv(3,3) = pv(3,3) + fc3(3)*particle_set ( index_c ) % r ( 3 )
           END DO

! 4x6
           DO iconst = 1, ng4x6

!  pv gets updated with FULL multiplier
            lambda_4x6 = lg4x6 ( iconst ) % lambda

            fc1=lambda_4x6(1)*lg4x6 ( iconst ) %fa+ &
                 lambda_4x6(2)*lg4x6 ( iconst ) %fb+ &
                 lambda_4x6(3)*lg4x6 ( iconst ) %fc
            fc2=-lambda_4x6(1)*lg4x6 ( iconst ) %fa+ &
                 lambda_4x6(4)*lg4x6 ( iconst ) %fd+ &
                 lambda_4x6(5)*lg4x6 ( iconst ) %fe
            fc3=-lambda_4x6(2)*lg4x6 ( iconst ) %fb- &
                 lambda_4x6(4)*lg4x6 ( iconst ) %fd+ &
                 lambda_4x6(6)*lg4x6 ( iconst ) %ff
            fc4=-lambda_4x6(3)*lg4x6 ( iconst ) %fc- &
                 lambda_4x6(5)*lg4x6 ( iconst ) %fe- &
                 lambda_4x6(6)*lg4x6 ( iconst ) %ff
            index_a = g4x6_list ( iconst ) % a + first_atom -1
            index_b = g4x6_list ( iconst ) % b + first_atom -1
            index_c = g4x6_list ( iconst ) % c + first_atom -1
            index_d = g4x6_list ( iconst ) % d + first_atom -1

!pv(1,1)
            pv(1,1) = pv(1,1) + fc1(1)*particle_set(index_a)%r(1)
            pv(1,1) = pv(1,1) + fc2(1)*particle_set(index_b)%r(1)
            pv(1,1) = pv(1,1) + fc3(1)*particle_set(index_c)%r(1)
            pv(1,1) = pv(1,1) + fc4(1)*particle_set(index_d)%r(1)
!pv(1,2)
            pv(1,2) = pv(1,2) + fc1(1)*particle_set(index_a)%r(2)
            pv(1,2) = pv(1,2) + fc2(1)*particle_set(index_b)%r(2)
            pv(1,2) = pv(1,2) + fc3(1)*particle_set(index_c)%r(2)
            pv(1,2) = pv(1,2) + fc4(1)*particle_set(index_d)%r(2)
!pv(1,3)
            pv(1,3) = pv(1,3) + fc1(1)*particle_set(index_a)%r(3)
            pv(1,3) = pv(1,3) + fc2(1)*particle_set(index_b)%r(3)
            pv(1,3) = pv(1,3) + fc3(1)*particle_set(index_c)%r(3)
            pv(1,3) = pv(1,3) + fc4(1)*particle_set(index_d)%r(3)
!pv(2,1)
            pv(2,1) = pv(2,1) + fc1(2)*particle_set(index_a)%r(1)
            pv(2,1) = pv(2,1) + fc2(2)*particle_set(index_b)%r(1)
            pv(2,1) = pv(2,1) + fc3(2)*particle_set(index_c)%r(1)
            pv(2,1) = pv(2,1) + fc4(2)*particle_set(index_d)%r(1)
!pv(2,2)
            pv(2,2) = pv(2,2) + fc1(2)*particle_set(index_a)%r(2)
            pv(2,2) = pv(2,2) + fc2(2)*particle_set(index_b)%r(2)
            pv(2,2) = pv(2,2) + fc3(2)*particle_set(index_c)%r(2)
            pv(2,2) = pv(2,2) + fc4(2)*particle_set(index_d)%r(2)
!pv(2,3)
            pv(2,3) = pv(2,3) + fc1(2)*particle_set(index_a)%r(3)
            pv(2,3) = pv(2,3) + fc2(2)*particle_set(index_b)%r(3)
            pv(2,3) = pv(2,3) + fc3(2)*particle_set(index_c)%r(3)
            pv(2,3) = pv(2,3) + fc4(2)*particle_set(index_d)%r(3)
!pv(3,1)
            pv(3,1) = pv(3,1) + fc1(3)*particle_set(index_a)%r(1)
            pv(3,1) = pv(3,1) + fc2(3)*particle_set(index_b)%r(1)
            pv(3,1) = pv(3,1) + fc3(3)*particle_set(index_c)%r(1)
            pv(3,1) = pv(3,1) + fc4(3)*particle_set(index_d)%r(1)
!pv(3,2)
            pv(3,2) = pv(3,2) + fc1(3)*particle_set(index_a)%r(2)
            pv(3,2) = pv(3,2) + fc2(3)*particle_set(index_b)%r(2)
            pv(3,2) = pv(3,2) + fc3(3)*particle_set(index_c)%r(2)
            pv(3,2) = pv(3,2) + fc4(3)*particle_set(index_d)%r(2)
!pv(3,3)
            pv(3,3) = pv(3,3) + fc1(3)*particle_set(index_a)%r(3)
            pv(3,3) = pv(3,3) + fc2(3)*particle_set(index_b)%r(3)
            pv(3,3) = pv(3,3) + fc3(3)*particle_set(index_c)%r(3)
            pv(3,3) = pv(3,3) + fc4(3)*particle_set(index_d)%r(3)

          END DO
        END DO
      END DO MOL
      CALL mp_sum ( pv, group )
      virial % pv_constraint = pv

END SUBROUTINE pv_constraint
!******************************************************************************

SUBROUTINE shake_roll( molecule, particle_set, pos, vel, r_shake, v_shake, &
                       simpar, max_sigma, ishake )
! updates the multiplier one molecule type at a time
    TYPE ( molecule_type ), POINTER                   :: molecule
    TYPE ( particle_type ), POINTER                   :: particle_set ( : )
    REAL (KIND=dp), INTENT ( INOUT )                    :: vel ( :, : ), pos ( :, : )
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)            :: r_shake, v_shake
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    INTEGER, INTENT(INOUT)                     :: ishake
    REAL(KIND=dp), INTENT(INOUT)                 :: max_sigma

! locals

    INTEGER                                  :: iconst, index_a, index_b, index_c
    INTEGER                                  :: first_atom, ndist
    REAL(KIND=dp)                                :: del_lam, dtby2, dtsqby2, &
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(KIND=dp), DIMENSION(3)                  :: f_roll_r, f_roll_v
    REAL(KIND=dp), DIMENSION(3)                  :: fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type ), POINTER       :: molecule_kind
    TYPE ( local_dist_constraint_type ), POINTER   :: ldist ( : )
    TYPE ( dist_constraint_type ), POINTER   :: dist_list ( : )
!
!------------------------------------------------------------------------------
! define local variables

  dtsqby2 = simpar%dt*simpar%dt*.5_dp
  dtby2 = simpar%dt*.5_dp
  molecule_kind => molecule % molecule_kind
  CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
  CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
  IF (ishake==1) THEN
     DO iconst = 1, ndist
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
        ldist ( iconst ) % f = -2.0_dp * ( ldist ( iconst ) % ra_old - &
                                           ldist ( iconst ) % rb_old )

! rotate fconst:
        f_roll_r = MATVEC_3X3 ( r_shake, ldist ( iconst ) % f )
        f_roll_v = MATVEC_3X3 ( v_shake, ldist ( iconst ) % f )
! update positions with full multiplier
         pos ( :, index_a ) = pos ( :, index_a ) + &
             imass1*dtsqby2*ldist ( iconst ) % lambda * f_roll_r
         pos ( :, index_b ) = pos ( :, index_b ) - &
             imass2*dtsqby2*ldist ( iconst ) % lambda * f_roll_r

! update velocites with full multiplier
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*ldist ( iconst ) % lambda * f_roll_v
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*ldist ( iconst ) % lambda * f_roll_v
     END DO

  ELSE

     DO iconst = 1, ndist
! computing the constraint
        index_a = dist_list ( iconst ) % a + first_atom -1
        index_b = dist_list ( iconst ) % b + first_atom -1
        r12 = pos ( :, index_a ) -  pos ( :, index_b )
        ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
                                                          dist_list ( iconst ) % dab
! computing fc(dt)
        fc_t = -2.0_dp*r12
        atomic_kind=>particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass1 = 1.0_dp/mass
        atomic_kind=>particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        imass2 = 1.0_dp/mass
! rotate fc:
        f_roll_r = MATVEC_3X3 ( r_shake, ldist ( iconst ) % f )
        f_roll_v = MATVEC_3X3 ( v_shake, ldist ( iconst ) % f )
        fdotf = DOTPROD_3D ( f_roll_r, fc_t )
        fdotf_sum = fdotf * imass1 + fdotf * imass2
        del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
        ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

! update positions by delta_lam
         pos ( :, index_a ) = pos ( :, index_a ) + &
             imass1*dtsqby2*del_lam * f_roll_r
         pos ( :, index_b ) = pos ( :, index_b ) - &
             imass2*dtsqby2*del_lam * f_roll_r

! update velocites by delta_lam
         vel ( :, index_a ) = vel ( :, index_a ) + &
             imass1*dtby2*del_lam * f_roll_v
         vel ( :, index_b ) = vel ( :, index_b ) - &
             imass2*dtby2*del_lam * f_roll_v
     END DO
  END IF

! computing the constraint and value of tolerance
  DO iconst = 1, ndist
     index_a = dist_list ( iconst ) % a + first_atom -1
     index_b = dist_list ( iconst ) % b + first_atom -1
     r12 ( : ) = pos ( :, index_a ) -  pos ( :, index_b )
     ldist ( iconst ) % sigma = dot_product(r12,r12) - dist_list ( iconst ) % dab * &
                                                       dist_list ( iconst ) % dab
     max_sigma = MAX(ABS(ldist ( iconst ) % sigma),max_sigma)
  END DO

END SUBROUTINE shake_roll
!
!!******************************************************************************
!
SUBROUTINE check_tol(roll_tol,iroll,char,matrix,veps)


    REAL(KIND=dp), INTENT(OUT)                   :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    CHARACTER(LEN=*), INTENT(IN)             :: char
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: matrix, veps

    REAL(KIND=dp)                                :: local_tol
    REAL(KIND=dp), DIMENSION(3, 3)               :: diff_rattle, diff_shake
    REAL(KIND=dp), DIMENSION(3, 3), SAVE         :: matrix_old, veps_old

!------------------------------------------------------------------------------

  SELECT CASE ( char )
  CASE ( 'SHAKE' )
     IF (iroll==1) THEN
        matrix_old = matrix
        roll_tol = -1.E10_dp
     ELSE
        roll_tol = 0.0_dp
        diff_shake = ABS ( matrix_old - matrix )
        local_tol = MAXVAL ( diff_shake )
        roll_tol = MAX ( roll_tol, local_tol )
        matrix_old = matrix
     END IF
     iroll = iroll + 1

  CASE ( 'RATTLE' )
     IF (iroll==1) THEN
        veps_old = veps
        roll_tol = -1.E+10_dp
     ELSE
        roll_tol = 0.0_dp

! compute tolerance on veps

        diff_rattle = ABS ( veps - veps_old )
        local_tol = MAXVAL ( diff_rattle )
        roll_tol = MAX ( roll_tol, local_tol )
        veps_old = veps

     END IF
     iroll = iroll + 1

  END SELECT

END SUBROUTINE check_tol

!!******************************************************************************

SUBROUTINE get_roll_matrix(char,r_shake,v_shake,vector_r,vector_v,u)


    CHARACTER(len=*), INTENT(IN)             :: char
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: r_shake, v_shake
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: u

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3, 3)               :: diag

!------------------------------------------------------------------------------

  IF ( PRESENT ( r_shake ) ) r_shake = 0.0_dp
  IF ( PRESENT ( v_shake ) ) v_shake = 0.0_dp
  diag = 0.0_dp

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
       diag = MATMUL_3x3 ( r_shake, v_shake )
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

END MODULE constraint
