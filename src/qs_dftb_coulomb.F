!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_dftb_coulomb [1.1] *
!!
!!   FUNCTION
!!     Calculation of Coulomb contributions in DFTB
!!
!!   AUTHOR
!!     JGH
!!
!!***
!******************************************************************************

MODULE qs_dftb_coulomb

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc,&
                                             real_to_scaled
  USE cp_control_types,                ONLY: dft_control_type,&
                                             dftb_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dgs,                             ONLY: dg_sum_patch,&
                                             dg_sum_patch_force_3d
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE input_constants,                 ONLY: do_rsgrid_distributed
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_charges
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE pme_tools,                       ONLY: get_center,&
                                             set_list
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_poisson_methods,              ONLY: pw_poisson_rebuild,&
                                             pw_poisson_solve
  USE pw_poisson_types,                ONLY: greens_fn_type,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_integral_a2b,&
                                             pw_p_type,&
                                             pw_transfer,&
                                             pw_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                             qs_dftb_pairpot_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             set_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_p_type, realspace_grid_type, rs2pw, &
       rs_grid_create, rs_grid_p_create, rs_grid_p_release, rs_grid_release, &
       rs_pw_transfer
  USE sparse_matrix_output,            ONLY: write_sparse_matrix
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix_set,&
                                             allocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             sparse_times_local,&
                                             cp_sm_scale_and_add,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_coulomb'

  PUBLIC :: build_dftb_coulomb

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_dftb_coulomb(ks_env,qs_env,ks_matrix,rho,mcharge,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(dp), DIMENSION(:)                   :: mcharge
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'build_dftb_coulomb', routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, do_ewald
    INTEGER                                  :: handle, istat, natom, &
      nbrow, nbcol, irow, icol, iatom, jatom, nmat, i, ikind, atom_i, &
      jkind, atom_j
    INTEGER, DIMENSION(:), POINTER           :: kind_of, atom_of_kind
    REAL(KIND=dp)                            :: dr, gmij, fi
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: mcin, mcout, gmcharge
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: gblock, ksblock, sblock, pblock,&
                                                dsblock
    REAL(KIND=dp), DIMENSION(3)              :: rij
    INTEGER, DIMENSION(:), POINTER           :: first_row, last_row, &
                                                first_col, last_col
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_p, gamma_matrix

!   ----------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.

    IF ( calculate_forces ) THEN
      nmat = 4
    ELSE
      nmat = 1
    END IF

    natom = SIZE(mcharge)
    ALLOCATE(gmcharge(natom,nmat),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    gmcharge = 0._dp

    gamma_matrix => qs_env%dft_control%qs_control%dftb_control%gamma_matrix

    ALLOCATE(mcin(1,natom),mcout(1,natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mcin(1,:) = mcharge(:)
    DO i=1,nmat
      mcout = 0._dp
      CALL sparse_times_local(gamma_matrix(i)%matrix,mcin,mcout,1,1)
      gmcharge(:,i) = mcout(1,:)
    END DO
    DEALLOCATE(mcin,mcout,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! 1/R contribution
    do_ewald = qs_env%dft_control%qs_control%dftb_control%do_ewald
    IF (do_ewald) THEN
      ! Ewald sum
      ewald_env => qs_env%dft_control%qs_control%dftb_control%ewald_env
!deb
      write(6,*) "DO EWALD"
!deb
    ELSE
      ! direct sum
      CALL get_qs_env(qs_env=qs_env,&
                      cell=cell,&
                      particle_set=particle_set,error=error)
      DO iatom=1,SIZE(particle_set)
        DO jatom=1,iatom-1
          rij = particle_set(iatom)%r - particle_set(jatom)%r
          rij = pbc(rij,cell)
          dr = SQRT(SUM(rij(:)**2))
          gmcharge(iatom,1)=gmcharge(iatom,1)+mcharge(jatom)/dr
          gmcharge(jatom,1)=gmcharge(jatom,1)+mcharge(iatom)/dr
          DO i=2,nmat
            gmcharge(iatom,i)=gmcharge(iatom,i)-rij(i-1)*mcharge(jatom)/dr**3
            gmcharge(jatom,i)=gmcharge(jatom,i)+rij(i-1)*mcharge(iatom)/dr**3
          END DO
        END DO
      END DO
      energy%core_overlap = 0._dp
      energy%core_self = 0._dp
    END IF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    force=force,para_env=para_env,error=error)
    CALL mp_sum(gmcharge,para_env%group)

    energy%hartree = energy%hartree + SUM(mcharge(:)*gmcharge(:,1))
    IF ( calculate_forces ) THEN
      ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               kind_of=kind_of,&
                               atom_of_kind=atom_of_kind)

      gmcharge(:,2)=gmcharge(:,2)*mcharge(:)
      gmcharge(:,3)=gmcharge(:,3)*mcharge(:)
      gmcharge(:,4)=gmcharge(:,4)*mcharge(:)
      DO iatom=1,natom
        ikind  = kind_of(iatom)
        atom_i = atom_of_kind(iatom)
        force(ikind)%rho_elec(1,atom_i) =&
            force(ikind)%rho_elec(1,atom_i) + gmcharge(iatom,2)
        force(ikind)%rho_elec(2,atom_i) =&
            force(ikind)%rho_elec(2,atom_i) + gmcharge(iatom,3)
        force(ikind)%rho_elec(3,atom_i) =&
            force(ikind)%rho_elec(3,atom_i) + gmcharge(iatom,4)
      END DO
    END IF

    IF ( .NOT. just_energy ) THEN
      CALL get_qs_env(qs_env=qs_env,&
                      matrix_s=matrix_s,error=error)
      matrix_p => rho%rho_ao
      IF ( calculate_forces .AND. SIZE(matrix_p) == 2) THEN
        CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,&
                matrix_p(2)%matrix,1.0_dp,error=error)
      END IF

      DO iatom=1,natom
        NULLIFY(block_node)
        block_node => first_block_node(ks_matrix(1)%matrix,iatom)
        DO WHILE (ASSOCIATED(block_node))
          NULLIFY(ksblock,sblock)
          CALL get_block_node(block_node=block_node,block_col=jatom,block=ksblock)
          CALL get_block_node(matrix=matrix_s(1)%matrix,block_row=iatom,&
                 block_col=jatom,block=sblock)
          gmij = 0.5_dp*(gmcharge(iatom,1)+gmcharge(jatom,1))
          ksblock = ksblock + gmij*sblock
          IF ( calculate_forces ) THEN
            ikind  = kind_of(iatom)
            atom_i = atom_of_kind(iatom)
            jkind  = kind_of(jatom)
            atom_j = atom_of_kind(jatom)
            NULLIFY(pblock)
            CALL get_block_node(matrix=matrix_p(1)%matrix,block_row=iatom,&
                 block_col=jatom,block=pblock)
            DO i=1,3
              NULLIFY(dsblock)
              CALL get_block_node(matrix=matrix_s(1+i)%matrix,block_row=iatom,&
                   block_col=jatom,block=dsblock)
              fi = gmij*SUM(pblock*dsblock)
              force(ikind)%rho_elec(i,atom_i) =&
                  force(ikind)%rho_elec(i,atom_i) + fi
              force(jkind)%rho_elec(i,atom_j) =&
                  force(jkind)%rho_elec(i,atom_j) - fi
            END DO
          END IF
          block_node => next_block_node(block_node)
        END DO
      END DO
      IF ( calculate_forces .AND. SIZE(matrix_p) == 2) THEN
        CALL cp_sm_scale_and_add(matrix_p(1)%matrix,1.0_dp,&
               matrix_p(2)%matrix,-1.0_dp,error=error)
      END IF
    END IF

    IF ( calculate_forces ) THEN
      DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(gmcharge,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(gamma_matrix)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_dftb_coulomb

! *****************************************************************************

  SUBROUTINE dftb_spme_evaluate (ewald_env,ewald_pw,particle_set,box,&
        fg_coulomb,vg_coulomb,pv_g,error)

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: box
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fg_coulomb
    REAL(KIND=dp), INTENT(OUT)               :: vg_coulomb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_g
!   TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dftb_spme_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, o_spline, rs_distributed_id, group, ierr, n, npart, &
               nsmax, ipart, jpart, p1, i, j
    INTEGER, DIMENSION(3)                    :: ndim
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, eps0, inv_fourpi_eps0,&
                                                ffa, dvols
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: center
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rhos
    REAL(KIND=dp), DIMENSION(3)              :: fat
    REAL(KIND=dp), DIMENSION(3, 3)           :: f_stress, h_stress
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(pw_grid_type), POINTER              :: grid_spme
    TYPE(pw_p_type), DIMENSION(3)            :: dphi_g
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: phi_g, rhob_g, rhob_r
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: drpot
    TYPE(realspace_grid_type), POINTER       :: rden

!------------------------------------------------------------------------------
!-------------- INITIALISATION ---------------------

  CALL timeset ( 'DFTB_SPME', 'I', ' ', handle )

  failure = .FALSE.

  CALL ewald_env_get ( ewald_env, alpha=alpha, eps0 = eps0, &
                       o_spline = o_spline, group = group, &
                       rs_distributed_id = rs_distributed_id, error=error)
  NULLIFY(green, poisson_env, pw_pool)
  CALL ewald_pw_get ( ewald_pw, pw_big_pool=pw_pool, poisson_env=poisson_env)
  CALL pw_poisson_rebuild(poisson_env,error=error)
  green => poisson_env%green_fft
  grid_spme => pw_pool % pw_grid

  inv_fourpi_eps0 = 1.0_dp / ( fourpi * eps0 )
  npart = SIZE ( particle_set )

  n = o_spline
  ALLOCATE ( rhos ( n, n, n ), STAT = ierr )
  CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
  ndim = n
  IF ( rs_distributed_id == do_rsgrid_distributed ) THEN
    nsmax = MAXVAL ( ndim(1:3) )
  ELSE
    nsmax = -1
  END IF
  CALL rs_grid_create ( rden, grid_spme, nsmax, error=error)
  rden % r = 0.0_dp

  ALLOCATE ( center ( 3, npart ), STAT = ierr )
  CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
  CALL get_center ( particle_set, box, center, grid_spme % npts, n )

!-------------- DENSITY CALCULATION ----------------

  ipart = 0
  jpart = 0

  DO

     CALL set_list ( particle_set, npart, center, p1, rden, ipart, jpart )
     IF ( p1 == 0 ) EXIT

     ! calculate function on small boxes
     CALL get_patch ( particle_set, box, green, grid_spme % npts, p1, rhos)

     ! add boxes to real space grid (big box)
     CALL dg_sum_patch ( rden, rhos, center(:,p1) )
     
  END DO     

  CALL pw_pool_create_pw ( pw_pool, rhob_r,  &
                           use_data = REALDATA3D, &
                           in_space = REALSPACE ,error=error)

  CALL rs_pw_transfer ( rden, rhob_r, rs2pw )

  ! transform density to G space and add charge function
  CALL pw_pool_create_pw ( pw_pool, rhob_g,  &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE ,error=error)
  CALL pw_transfer ( rhob_r, rhob_g )
  ! update charge function
  rhob_g % cc = rhob_g % cc * green % p3m_charge % cr

!-------------- ELECTROSTATIC CALCULATION -----------

  ! allocate intermediate arrays
  DO i = 1, 3
     NULLIFY(dphi_g(i)%pw)
     CALL pw_pool_create_pw ( pw_pool, dphi_g ( i )%pw,  &
                              use_data = COMPLEXDATA1D,   &
                              in_space = RECIPROCALSPACE ,error=error)
  END DO
  CALL pw_pool_create_pw ( pw_pool, phi_g,  &
                           use_data = COMPLEXDATA1D,   &
                           in_space = RECIPROCALSPACE ,error=error)
  CALL pw_poisson_solve ( poisson_env, rhob_g, vg_coulomb, phi_g, dphi_g, &
                          h_stress ,error=error)

  CALL pw_pool_give_back_pw ( pw_pool, rhob_g ,error=error)
  vg_coulomb = vg_coulomb * inv_fourpi_eps0

!---------- END OF ELECTROSTATIC CALCULATION --------

!------------- STRESS TENSOR CALCULATION ------------

  DO i = 1, 3
    DO j = i, 3
       f_stress(i,j) = pw_integral_a2b(dphi_g(i)%pw,dphi_g(j)%pw)
       f_stress(j,i) = f_stress(i,j)
    END DO
  END DO
  ffa = ( 1.0_dp / fourpi ) * ( 0.5_dp / alpha ) ** 2
  f_stress = -ffa * f_stress
  pv_g = inv_fourpi_eps0 * ( h_stress + f_stress )

!--------END OF STRESS TENSOR CALCULATION -----------

  ! move derivative of potential to real space grid and
  ! multiply by charge function in g-space
  CALL rs_grid_p_create ( drpot, 3, grid_spme, nsmax, error=error)
  DO i = 1, 3
     dphi_g ( i ) % pw % cc = dphi_g ( i ) % pw % cc * &
                              green % p3m_charge % cr
     CALL pw_transfer ( dphi_g ( i )%pw, rhob_r )
     CALL pw_pool_give_back_pw ( pw_pool, dphi_g ( i )%pw ,error=error)
     CALL rs_pw_transfer ( drpot ( i ) % rs_grid, rhob_r, pw2rs )
  END DO

  CALL pw_pool_give_back_pw ( pw_pool, rhob_r ,error=error)
  CALL pw_pool_give_back_pw ( pw_pool, phi_g ,error=error)

!----------------- FORCE CALCULATION ----------------

! initialize the forces
  fg_coulomb = 0.0_dp
  dvols = grid_spme % dvol * inv_fourpi_eps0

  ipart = 0
  jpart = 0
  DO

     CALL set_list ( particle_set, npart, center, p1, rden, ipart, jpart )
     IF ( p1 == 0 ) EXIT

     ! calculate function on small boxes
     CALL get_patch ( particle_set, box, green, grid_spme % npts, p1, rhos )

     ! add boxes to real space grid (big box)
     CALL dg_sum_patch_force_3d ( drpot, rhos, center(:,p1), fat )
     fg_coulomb ( 1, p1 ) = fg_coulomb ( 1, p1 ) - fat ( 1 ) * dvols
     fg_coulomb ( 2, p1 ) = fg_coulomb ( 2, p1 ) - fat ( 2 ) * dvols
     fg_coulomb ( 3, p1 ) = fg_coulomb ( 3, p1 ) - fat ( 3 ) * dvols

  END DO

!--------------END OF FORCE CALCULATION -------------

!------------------CLEANING UP ----------------------

  CALL rs_grid_release ( rden ,error=error)
  CALL rs_grid_p_release ( drpot ,error=error)
  DEALLOCATE ( rhos, STAT = ierr )
  CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE ( center, STAT = ierr )
  CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  CALL timestop ( 0.0_dp, handle )

END SUBROUTINE dftb_spme_evaluate
!******************************************************************************

SUBROUTINE get_patch ( part, box, green, npts, p, rhos )

    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), POINTER                 :: box
    TYPE(greens_fn_type), POINTER            :: green
    INTEGER, DIMENSION(3), INTENT(IN)        :: npts
    INTEGER, INTENT(IN)                      :: p
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: rhos

    INTEGER                                  :: nbox
    REAL(KIND=dp)                            :: q
    REAL(KIND=dp), DIMENSION(3)              :: delta, r

!------------------------------------------------------------------------------

  nbox = SIZE ( rhos, 1 )
  r = part(p)% r
  CALL get_atomic_kind ( atomic_kind=part(p)%atomic_kind, qeff=q)
 
  CALL get_delta ( box, r, npts, delta, nbox )

  CALL spme_get_patch ( rhos, nbox, delta, q, green % p3m_coeff )

END SUBROUTINE get_patch

SUBROUTINE spme_get_patch ( rhos, n, delta, q, coeff )

    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: rhos
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: delta
    REAL(KIND=dp), INTENT(IN)                :: q
    REAL(KIND=dp), &
      DIMENSION(-(n-1):n-1, 0:n-1), &
      INTENT(IN)                             :: coeff

    INTEGER, PARAMETER                       :: nmax = 12

    INTEGER                                  :: i, i1, i2, i3, j, l
    REAL(KIND=dp)                            :: r2, r3
    REAL(KIND=dp), DIMENSION(3, -nmax:nmax)  :: w_assign
    REAL(KIND=dp), DIMENSION(3, 0:nmax-1)    :: deltal
    REAL(KIND=dp), DIMENSION(3, 1:nmax)      :: f_assign

!------------------------------------------------------------------------------

!  CPPostcondition(n <= nmax,cp_failure_level,routineP,error,failure)

! calculate the assignment function values and
! the charges on the grid (small box)

  deltal ( 1, 0 ) = 1.0_dp
  deltal ( 2, 0 ) = 1.0_dp
  deltal ( 3, 0 ) = 1.0_dp
  DO l = 1, n-1
     deltal ( 1, l ) = deltal ( 1, l-1 ) * delta ( 1 )
     deltal ( 2, l ) = deltal ( 2, l-1 ) * delta ( 2 )
     deltal ( 3, l ) = deltal ( 3, l-1 ) * delta ( 3 )
  END DO

  w_assign = 0.0_dp
  DO j = -(n-1), n-1, 2
    DO l = 0, n-1
      w_assign ( 1, j ) =  w_assign ( 1, j ) + coeff ( j, l ) * deltal ( 1, l )
      w_assign ( 2, j ) =  w_assign ( 2, j ) + coeff ( j, l ) * deltal ( 2, l )
      w_assign ( 3, j ) =  w_assign ( 3, j ) + coeff ( j, l ) * deltal ( 3, l )
    END DO
  END DO
  DO i = 1, n
    j = n + 1 - 2 * i
    f_assign (1, i ) = w_assign ( 1, j )
    f_assign (2, i ) = w_assign ( 2, j )
    f_assign (3, i ) = w_assign ( 3, j )
  END DO

  DO i3 = 1, n
    r3 = q * f_assign ( 3, i3 )
    DO i2 = 1, n
      r2 = r3 * f_assign ( 2, i2 )
      DO i1 = 1, n
        rhos ( i1, i2, i3 ) = r2 * f_assign ( 1, i1 )
      END DO
    END DO
  END DO

END SUBROUTINE spme_get_patch

SUBROUTINE get_delta ( box, r, npts, delta, n )

    TYPE(cell_type), POINTER                 :: box
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    INTEGER, DIMENSION(3), INTENT(IN)        :: npts
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: delta
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: mp
    INTEGER, DIMENSION(3)                    :: center
    REAL(KIND=dp)                            :: rmp
    REAL(KIND=dp), DIMENSION(3)              :: ca, grid_i, s

!------------------------------------------------------------------------------

  mp = MAXVAL ( npts(:) )
  rmp = REAL ( mp,KIND=dp)
! compute the scaled coordinate of atomi
  CALL real_to_scaled(s,r,box)
  s = s - REAL ( NINT ( s ),KIND=dp)

! find the continuous ``grid'' point
  grid_i ( 1:3 ) = REAL ( npts ( 1:3 ),KIND=dp) * s ( 1:3 )

! find the closest grid point

  IF ( MOD ( n, 2 ) == 0 ) THEN
    center ( : ) = INT ( grid_i ( : ) + rmp ) - mp
    ca ( : ) = REAL ( center ( : ) ) + 0.5_dp
  ELSE
    center ( : ) = NINT ( grid_i ( : ) )
    ca ( : ) = REAL ( center ( : ) )
  END IF

! find the distance vector
  delta ( : ) = grid_i ( : ) - ca ( : )

END SUBROUTINE get_delta

! *****************************************************************************

END MODULE qs_dftb_coulomb

! *****************************************************************************

