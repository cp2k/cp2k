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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_get
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type, &
                                             fist_nonbond_env_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing, ONLY : mp_sum
  USE pair_potential,                  ONLY: potential_s
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE simulation_cell,                 ONLY: cell_type,&
                                             scaled_to_real
  USE splines,                         ONLY: spline_environment_type, &
                                             spline_data
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: force_nonbond, bonded_correct_gaussian

!!***
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE force_nonbond ( fist_nonbond_env, particle_set, cell,  &
                           pot_nonbond, f_nonbond, ptens_nonbond )

! Calculates the force and the potential of the minimum image, and
! the pressure tensor
!

    TYPE ( fist_nonbond_env_type ), POINTER :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: cell 
    REAL(KIND=dp), INTENT(OUT)                   :: pot_nonbond
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: f_nonbond, ptens_nonbond

    CHARACTER(LEN=*), PARAMETER :: module_name = "fist_nonbond_force", &
      routine_name = "force_nonbond", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: ab, atom_a, atom_b, &
                                                cell_b( 3 ), handle, ikind, &
                                                ilist, inode, jkind, nkinds, &
                                                nlist, nnode, nexl
    REAL(KIND=dp)                                :: energy, flops, fscalar, rab2, &
                                                rab2_max
    REAL(KIND=dp), DIMENSION(3)                  :: dra, drb, fr, rab
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    INTEGER, DIMENSION ( :, : ), POINTER     :: spltab
    INTEGER, DIMENSION ( : ), POINTER        :: list
    TYPE ( spline_data ),  DIMENSION ( : ), POINTER :: spline_table
    TYPE ( neighbor_list_set_p_type ), &
      DIMENSION ( : ),  POINTER                  :: nonbonded
    REAL(KIND=dp), DIMENSION(:, :), POINTER   :: r_last_update
    REAL(KIND=dp), DIMENSION(:), POINTER   :: rab_last_update
    TYPE ( spline_environment_type ), POINTER     :: spline_env
    TYPE ( pair_potential_type ), POINTER     :: potparm ( :, : )

!MK
!------------------------------------------------------------------------------

  CALL timeset ( routine_name, 'E', 'Mflops', handle )
  flops = 0.0_dp
  CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, &
                             r_last_update = r_last_update,            &
                             spline_env = spline_env, natom_types = nkinds )
! initializing the potential energy, pressure tensor and force and spline pointers
  pot_nonbond = 0.0_dp
  f_nonbond = 0.0_dp
  ptens_nonbond = 0.0_dp
  spline_table => spline_env % spline_table
  spltab => spline_env % spltab
  potparm => spline_env % potparm

!
! starting the force loop
!
    DO ikind=1,nkinds
      DO jkind=1,nkinds

        ab = ikind + nkinds*(jkind - 1)

        neighbor_list_set => nonbonded(ab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        rab2_max = potparm ( ikind, jkind ) % rcutsq 

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=atom_a,&
                                 nnode=nnode)

          dra(:) = particle_set(atom_a)%r(:) - &
                   scaled_to_real ( r_last_update(:,atom_a), cell )

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=atom_b,&
                                   cell=cell_b,&
                                   s=rab)

            drb(:) = particle_set(atom_b)%r(:) - &
                     scaled_to_real ( r_last_update(:,atom_b), cell )
            rab ( : ) =  scaled_to_real ( rab ( : ), cell ) - dra(:) + drb(:)
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

            IF (rab2 <= rab2_max) THEN
              CALL potential_s ( spline_table, spltab, rab2,  &
                                ikind, jkind, energy, fscalar )
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

! computing long range corrections to the potential
!
!when we get the system for handling multiple potentials
!we will add the long range correction
!
!       pot_nonbond=pot_nonbond+lrc*(1./box % deth)
!

  CALL timestop ( flops, handle )

END SUBROUTINE force_nonbond

!******************************************************************************

SUBROUTINE bonded_correct_gaussian ( atomic_kind_set, local_particles, &
                                     particle_set, ewald_env, &
                                     v_bonded_corr, pv_bc )


    TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
    TYPE ( distribution_1d_type ), POINTER :: local_particles
    TYPE ( particle_type ), POINTER :: particle_set ( : )   
    TYPE(ewald_environment_type), POINTER  :: ewald_env
    REAL(KIND=dp), INTENT(OUT)                   :: v_bonded_corr
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: pv_bc

    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    INTEGER                                  :: handle, iat, imol, i, j, nmol
    REAL(KIND=dp)                                :: alpha, arg, const, dij, &
                                                e_arg_arg, errf, flops, &
                                                fscalar, idij, ifourpieps0, &
                                                qi, qj, rijsq, tc, eps0
    REAL(KIND=dp), DIMENSION(3)                  :: ri, rj, rij
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    INTEGER                                  :: natoms_per_kind, ikind, nkind, &
                                                jkind, iexl, nexl, iatom, group
    INTEGER, DIMENSION ( : ), POINTER        :: list

! parameters for the numerical erf in Abramowitz and Stegun Eq. 7.1.26 pg 299
!------------------------------------------------------------------------------

  CALL timeset ( 'BOND_CORRECT_GAUSSIAN', 'E', 'Mflops', handle )
  flops = 0.0_dp

! defining the constants
  pv_bc = 0.0_dp
  v_bonded_corr = 0.0_dp

!
!
  CALL ewald_env_get ( ewald_env, eps0 = eps0, alpha = alpha, group = group )
  const = 2.0_dp * alpha / SQRT( pi )
  ifourpieps0 = 1.0_dp / ( 4.0_dp * pi * eps0 )

  nkind = SIZE ( atomic_kind_set )
  DO ikind = 1, nkind
     natoms_per_kind = local_particles % n_el ( ikind )
     DO iatom = 1, natoms_per_kind
       i = local_particles % list ( ikind ) % array ( iatom )
       atomic_kind => particle_set ( i ) % atomic_kind
       ri ( : ) = particle_set ( i ) % r ( : )
       CALL get_atomic_kind ( atomic_kind,qeff = qi )
       DO jkind = 1, nkind 
         list => particle_set ( i ) % lpi % exclusion ( jkind ) % list
         IF(.NOT.ASSOCIATED(list)) CYCLE
         nexl = SIZE ( list )
         DO iexl = 1, nexl
           j =  list ( iexl )
           IF ( j==i ) CYCLE
           atomic_kind => particle_set ( j ) % atomic_kind
           rj ( : ) = particle_set ( j ) % r ( : )
           CALL get_atomic_kind ( atomic_kind,qeff = qj )
           rij(1) = ri(1) - rj(1)
           rij(2) = ri(2) - rj(2)
           rij(3) = ri(3) - rj(3)
           rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
           dij = SQRT ( rijsq )
           idij = 1.0_dp / dij
           arg = alpha * dij
           e_arg_arg = EXP ( -arg ** 2 )
           tc = 1.0_dp / ( 1.0_dp + pc * arg )

! defining errf=1-erfc
           errf = 1.0_dp &
                - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg

! getting the potential
           v_bonded_corr = v_bonded_corr - qi*qj*idij*errf

! subtracting the force from the total force
           fscalar = ifourpieps0 * qi * qj * idij ** 2 &
                * ( idij * errf - const * EXP ( -alpha ** 2 * rijsq ) )

           particle_set ( i ) % f(1) = particle_set ( i ) % f(1) - fscalar*rij(1)
           particle_set ( i ) % f(2) = particle_set ( i ) % f(2) - fscalar*rij(2)
           particle_set ( i ) % f(3) = particle_set ( i ) % f(3) - fscalar*rij(3)

           pv_bc(1,1) = pv_bc(1,1) - fscalar * rij(1) * particle_set ( i ) % r(1)
           pv_bc(1,2) = pv_bc(1,2) - fscalar * rij(1) * particle_set ( i ) % r(2)
           pv_bc(1,3) = pv_bc(1,3) - fscalar * rij(1) * particle_set ( i ) % r(3)
           pv_bc(2,1) = pv_bc(2,1) - fscalar * rij(2) * particle_set ( i ) % r(1)
           pv_bc(2,2) = pv_bc(2,2) - fscalar * rij(2) * particle_set ( i ) % r(2)
           pv_bc(2,3) = pv_bc(2,3) - fscalar * rij(2) * particle_set ( i ) % r(3)
           pv_bc(3,1) = pv_bc(3,1) - fscalar * rij(3) * particle_set ( i ) % r(1)
           pv_bc(3,2) = pv_bc(3,2) - fscalar * rij(3) * particle_set ( i ) % r(2)
           pv_bc(3,3) = pv_bc(3,3) - fscalar * rij(3) * particle_set ( i ) % r(3)

           flops = flops + 62.0_dp
         ENDDO
       ENDDO
     ENDDO
   ENDDO

! the factor of 1/2 comes from double counting in the exclusion list
  v_bonded_corr = v_bonded_corr * ifourpieps0 * 0.5_dp
  CALL mp_sum ( v_bonded_corr, group )

  flops = flops * 1.0E-6_dp
  CALL timestop ( flops, handle )

END SUBROUTINE bonded_correct_gaussian

!******************************************************************************

END MODULE fist_nonbond_force
