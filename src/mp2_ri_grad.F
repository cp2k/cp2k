!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate gradients of RI-GPW-MP2 energy using pw
!> \par History
!>      10.2013 created [Mauro Del Ben]
! *****************************************************************************
MODULE mp2_ri_grad
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_copy_into_existing, &
       cp_dbcsr_create, cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_set, cp_dbcsr_transposed, cp_dbcsr_type, &
       dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_indxg2l,&
                                             cp_fm_indxg2p,&
                                             cp_fm_indxl2g,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE message_passing,                 ONLY: &
       mp_alltoall, mp_comm_split_direct, mp_irecv, mp_isend, &
       mp_request_null, mp_sendrecv, mp_sum, mp_wait, mp_waitall
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                             integ_mat_buffer_type_2D,&
                                             mp2_type
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_derive,&
                                             pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                             calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace,&
                                             integrate_v_rspace
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_p_type,&
                                             rs_grid_release,&
                                             rs_grid_retain
  USE rs_pw_interface,                 ONLY: potential_pw2rs
  USE task_list_types,                 ONLY: task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_grad'

  PUBLIC :: calc_ri_mp2_nonsep


  CONTAINS

! *****************************************************************************
!> \brief Calcualte the non-separable part of the gradients and update the 
!>        Lagrangian
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param dft_control ...
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param Eigenval ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param mat_munu ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param task_list_sub ...
!> \param blacs_env_sub ...
!> \param error ...
!> \param error_sub ...
!> \author Mauro Del Ben 
! *****************************************************************************
  SUBROUTINE calc_ri_mp2_nonsep(qs_env,mp2_env,para_env,para_env_sub,dft_control,cell,particle_set,&
                                atomic_kind_set,qs_kind_set,mo_coeff,nmo,homo,dimen_RI,Eigenval,&
                                my_group_L_start,my_group_L_end,my_group_L_size,rho_r,rho_g,pot_g,&
                                mat_munu,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                                blacs_env_sub,error,error_sub)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER                                  :: nmo, homo, dimen_RI
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: my_group_L_start, &
                                                my_group_L_end, &
                                                my_group_L_size
    TYPE(pw_p_type)                          :: rho_r, rho_g, pot_g
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(pw_env_type), POINTER               :: pw_env_sub
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(task_list_type), POINTER            :: task_list_sub
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ri_mp2_nonsep', &
      routineP = moduleN//':'//routineN

    INTEGER :: alpha, atom_a, beta, dimen, dir, handle, handle2, handle3, i, &
      iatom, igrid_level, ikind, iorb, ipgf, iset, istat, itmp(2), L_counter, &
      lb(3), LLL, location(3), my_P_end, my_P_size, my_P_start, na1, na2, &
      natom, ncoa, nkind, nseta, offset, sgfa, stat, tp(3), ub(3), virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    INTEGER, DIMENSION(3)                    :: comp
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, map_it_here, &
                                                skip_shell, use_virial
    REAL(KIND=dp)                            :: e_hartree, eps_filter, &
                                                pair_energy, rab2, total_rho
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: wf_vector
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: G_PQ_local, I_ab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress, my_virial_a, &
                                                my_virial_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: I_tmp2, pab, rpgfa, sphi_a, &
                                                zeta
    TYPE(cp_dbcsr_p_type)                    :: Lag_mu_i_1, Lag_nu_a_2, &
                                                matrix_P_inu, matrix_P_munu, &
                                                matrix_P_munu_nosym
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: G_P_ia
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_v
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: L1_mu_i, L2_nu_a
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(pw_p_type)                          :: dvg(3), psi_L, temp_pw_g
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)
    failure=.FALSE.

    dimen=nmo
    virtual=dimen-homo
    eps_filter=mp2_env%mp2_gpw%eps_filter
    NULLIFY(mo_coeff_o, mo_coeff_v, G_P_ia, ks_env)
    mo_coeff_o => mp2_env%ri_grad%mo_coeff_o
    mo_coeff_v => mp2_env%ri_grad%mo_coeff_v
    G_P_ia => mp2_env%ri_grad%G_P_ia

    CALL get_qs_env(qs_env, ks_env=ks_env, molecule_set=molecule_set, &
                    molecule_kind_set=molecule_kind_set, error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (kind_of(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

    itmp=get_limit(dimen_RI,para_env_sub%num_pe,para_env_sub%mepos)
    my_P_start=itmp(1)
    my_P_end=itmp(2)
    my_P_size=itmp(2)-itmp(1)+1
    ALLOCATE(G_PQ_local(dimen_RI,my_group_L_size),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    G_PQ_local=0.0_dp
    G_PQ_local(my_P_start:my_P_end,1:my_group_L_size)=mp2_env%ri_grad%Gamma_PQ
    DEALLOCATE(mp2_env%ri_grad%Gamma_PQ,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL mp_sum(G_PQ_local,para_env_sub%group)

    ! deallocate here PQ_half, maybe usefull in the future
    DEALLOCATE(mp2_env%ri_grad%PQ_half,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! create matrix holding the back transformation (G_P_inu)
    ALLOCATE(matrix_P_inu%matrix)
    CALL cp_dbcsr_init(matrix_P_inu%matrix,error=error_sub)
    CALL cp_dbcsr_create(matrix_P_inu%matrix,template=mo_coeff_o,error=error_sub)

    ! non symmetric matrix
    ALLOCATE(matrix_P_munu_nosym%matrix)
    CALL cp_dbcsr_init(matrix_P_munu_nosym%matrix,error=error_sub)
    CALL cp_dbcsr_create(matrix_P_munu_nosym%matrix,template=mat_munu%matrix,&
                         matrix_type=dbcsr_type_no_symmetry,error=error_sub)
    ! CALL cp_dbcsr_create(matrix_P_munu_nosym%matrix,template=mat_munu%matrix,error=error_sub)
    ! CALL cp_dbcsr_copy(matrix_P_munu_nosym%matrix,mat_munu%matrix,&
    !                    name="matrix_P_munu_nosym",error=error_sub)
    !XXXXXX
    ! CALL cp_dbcsr_desymmetrize(mat_munu%matrix,matrix_P_munu_nosym%matrix,error=error_sub)

    ! create Lagrangian matrices in mixed AO/MO formalism
    ALLOCATE(Lag_mu_i_1%matrix)
    CALL cp_dbcsr_init(Lag_mu_i_1%matrix,error=error_sub)
    CALL cp_dbcsr_create(Lag_mu_i_1%matrix,template=mo_coeff_o,error=error_sub)
    CALL cp_dbcsr_set(Lag_mu_i_1%matrix,0.0_dp,error=error_sub)

    ALLOCATE(Lag_nu_a_2%matrix)
    CALL cp_dbcsr_init(Lag_nu_a_2%matrix,error=error_sub)
    CALL cp_dbcsr_create(Lag_nu_a_2%matrix,template=mo_coeff_v,error=error_sub)
    CALL cp_dbcsr_set(Lag_nu_a_2%matrix,0.0_dp,error=error_sub)

    ALLOCATE(matrix_P_munu%matrix)
    CALL cp_dbcsr_init(matrix_P_munu%matrix,error=error_sub)
    CALL cp_dbcsr_create(matrix_P_munu%matrix,template=mat_munu%matrix,&
                         matrix_type=dbcsr_type_symmetric,error=error_sub)
    ! CALL cp_dbcsr_set(mat_munu%matrix,0.0_dp,error=error_sub)
    ! CALL cp_dbcsr_copy(matrix_P_munu%matrix,mat_munu%matrix,&
    !                    name="matrix_P_munu",error=error_sub)   

    ! wave function vector and supporting stuff
    ALLOCATE(wf_vector(dimen_RI),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(psi_L%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,psi_L%pw,&
                           use_data=REALDATA3D,&
                           in_space=REALSPACE,error=error_sub)

    ! get forces
    NULLIFY (force, virial)
    CALL get_qs_env(qs_env=qs_env,force=force,virial=virial,error=error)

    ! check if we want to calculate the virial
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
 
    ! in case virial is required we need auxilliary pw 
    ! for calculate the MP2-volume contribution to the virial
    ! (hartree potential derivatives)   
    IF(use_virial) THEN
      NULLIFY(temp_pw_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,temp_pw_g%pw,&
                             use_data=COMPLEXDATA1D,&
                             in_space=RECIPROCALSPACE,error=error_sub)
      DO i=1, 3
        NULLIFY(dvg(i)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool,dvg(i)%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,error=error_sub)
      END DO
    END IF

    ! start main loop over auxiliary basis functions
    CALL timeset(routineN//"_loop",handle2)
    L_counter=0
    DO LLL=my_group_L_start, my_group_L_end
      L_counter=L_counter+1

      ! first back-transformation a->nu
      CALL timeset(routineN//"_back_v",handle3)
      CALL cp_dbcsr_set(matrix_P_inu%matrix,0.0_dp,error=error_sub)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, mo_coeff_v, G_P_ia(L_counter)%matrix, &
                              0.0_dp, matrix_P_inu%matrix, filter_eps=eps_filter, error=error_sub)
      CALL timestop(handle3)

      ! second back-transformation i->mu
      CALL timeset(routineN//"_back_o",handle3)
      CALL cp_dbcsr_set(matrix_P_munu_nosym%matrix,0.0_dp,error=error_sub)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_P_inu%matrix, mo_coeff_o, &
                              0.0_dp, matrix_P_munu_nosym%matrix, filter_eps=eps_filter, error=error_sub)

      ! ! first back-transformation a->nu
      ! CALL timeset(routineN//"_back_v",handle3)
      ! CALL cp_dbcsr_set(matrix_P_inu%matrix,0.0_dp,error=error_sub)
      ! CALL cp_dbcsr_multiply("N", "T", 1.0_dp, G_P_ia(L_counter)%matrix, mo_coeff_v, &
      !                         0.0_dp, matrix_P_inu%matrix, filter_eps=eps_filter, error=error_sub)
      ! CALL timestop(handle3)

      ! ! second back-transformation i->mu
      ! CALL timeset(routineN//"_back_o",handle3)
      ! CALL cp_dbcsr_set(matrix_P_munu_nosym%matrix,0.0_dp,error=error_sub)
      ! CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mo_coeff_o, matrix_P_inu%matrix, &
      !                         0.0_dp, matrix_P_munu_nosym%matrix, filter_eps=eps_filter, error=error_sub)
      ! symmetrize
      CALL cp_dbcsr_set(matrix_P_munu%matrix,0.0_dp,error=error_sub)
      CALL cp_dbcsr_transposed(matrix_P_munu%matrix,matrix_P_munu_nosym%matrix,&
                               error=error_sub)
      CALL cp_dbcsr_add(matrix_P_munu%matrix,matrix_P_munu_nosym%matrix,&
                        alpha_scalar=2.0_dp,beta_scalar=2.0_dp,error=error)
      ! this is a trick to avoid that integrate_v_rspace starts to cry
      CALL cp_dbcsr_copy_into_existing(mat_munu%matrix, matrix_P_munu%matrix, error)
      CALL cp_dbcsr_copy(matrix_P_munu%matrix,mat_munu%matrix,error=error_sub)

      CALL timestop(handle3)

      ! calculate potential associted to the single aux function
      CALL timeset(routineN//"_wf_pot",handle3)
      wf_vector=0.0_dp
      wf_vector(LLL)=1.0_dp
      ! pseudo psi_L
      CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set,&
                                  qs_kind_set,cell,dft_control,particle_set,pw_env_sub,&
                                  basis_type="RI_AUX",&
                                  external_vector=wf_vector,&
                                  error=error_sub)

      rho_r%pw%cr3d = psi_L%pw%cr3d
      CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
      CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
      CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
      CALL timestop(handle3)

      IF(use_virial) THEN
        ! make a copy of the density in G space
        ! calculate the potential derivatives in G space 
        CALL timeset(routineN//"_Virial",handle3)
        CALL pw_copy(rho_g%pw, temp_pw_g%pw, error=error)
        DO i=1, 3
          comp=0
          comp(i)=1 
          CALL pw_copy(pot_g%pw, dvg(i)%pw, error=error_sub)
          CALL pw_derive(dvg(i)%pw, comp, error=error_sub)
        END DO
        CALL timestop(handle3)
      END IF
 
      ! integrate the potential of the single gaussian and update 
      ! 2-center forces with Gamma_PQ
      CALL timeset(routineN//"_int_PQ",handle3)
      NULLIFY(rs_v)
      NULLIFY(rs_descs)
      CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v, error=error_sub)
      DO i=1,SIZE(rs_v)
        CALL rs_grid_retain(rs_v(i)%rs_grid,error=error_sub)
      END DO
      CALL potential_pw2rs(rs_v,rho_r,pw_env_sub,error_sub)

      offset=0
      DO iatom=1, natom
        ikind=kind_of(iatom)
        atom_a = atom_of_kind(iatom)
        CALL get_qs_kind(qs_kind=qs_kind_set(ikind),basis_set=basis_set_a,&
                         basis_type="RI_AUX",error=error)

        first_sgfa   =>  basis_set_a%first_sgf
        la_max       =>  basis_set_a%lmax
        la_min       =>  basis_set_a%lmin
        npgfa        =>  basis_set_a%npgf
        nseta        =   basis_set_a%nset
        nsgfa        =>  basis_set_a%nsgf_set
        rpgfa        =>  basis_set_a%pgf_radius
        set_radius_a =>  basis_set_a%set_radius
        sphi_a       =>  basis_set_a%sphi
        zeta         =>  basis_set_a%zet

        ra(:) = pbc(particle_set(iatom)%r,cell)
        rab=0.0_dp
        rab2=0.0_dp

        force_a(:) = 0.0_dp
        force_b(:) = 0.0_dp
        IF (use_virial) THEN
          my_virial_a = 0.0_dp
          my_virial_b = 0.0_dp
        END IF

        DO iset=1, nseta
         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         ALLOCATE(I_tmp2(ncoa,1))
         I_tmp2=0.0_dp
         ALLOCATE(I_ab(nsgfa(iset),1))
         I_ab=0.0_dp
         ALLOCATE(pab(ncoa,1))
         pab=0.0_dp

         I_ab(1:nsgfa(iset),1)=-4.0_dp*G_PQ_local(offset+1:offset+nsgfa(iset),L_counter)

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                     1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                     I_ab(1,1),SIZE(I_ab,1),&
                     0.0_dp,pab(1,1),SIZE(pab,1))

         I_ab=0.0_dp

         igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,MINVAL(zeta(:,iset)))
         map_it_here=.FALSE.
         IF (.NOT. ALL (rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
            DO dir = 1,3
                  ! bounds of local grid (i.e. removing the 'wings'), if periodic
                  tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                  tp(dir) = MODULO ( tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir) )
                  IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                     lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir ) + rs_v(igrid_level)%rs_grid%desc%border
                     ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir ) - rs_v(igrid_level)%rs_grid%desc%border
                  ELSE
                     lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir )
                     ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir )
                  ENDIF
                  ! distributed grid, only map if it is local to the grid
                  location(dir)=tp(dir)+rs_v(igrid_level)%rs_grid%desc%lb(dir)
            ENDDO
            IF  (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                 lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                 lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
               map_it_here=.TRUE.
            ENDIF
         ELSE
            ! not distributed, just a round-robin distribution over the full set of CPUs
            IF (MODULO(offset,para_env_sub%num_pe)==para_env_sub%mepos) map_it_here=.TRUE.
         ENDIF

         offset=offset+nsgfa(iset)

         IF (map_it_here) THEN
           DO ipgf=1, npgfa(iset)
             na1=(ipgf - 1)*ncoset(la_max(iset)) + 1
             na2=ipgf*ncoset(la_max(iset))
             igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,zeta(ipgf,iset))

             CALL integrate_pgf_product_rspace(la_max=la_max(iset),zeta=zeta(ipgf,iset)/2.0_dp,la_min=la_min(iset),&
                                               lb_max=0,zetb=zeta(ipgf,iset)/2.0_dp,lb_min=0,&
                                               ra=ra,rab=rab,rab2=rab2,&
                                               rsgrid=rs_v(igrid_level)%rs_grid,&
                                               cell=cell,&
                                               cube_info=pw_env_sub%cube_info(igrid_level),&
                                               hab=I_tmp2,&
                                               pab=pab,&
                                               o1=na1-1,&
                                               o2=0,&
                                               map_consistent=.TRUE.,&
                                               eps_gvg_rspace=dft_control%qs_control%eps_gvg_rspace,&
                                               calculate_forces=.TRUE.,&
                                               force_a=force_a, force_b=force_b, &
                                               use_virial=use_virial,my_virial_a=my_virial_a,my_virial_b=my_virial_b,&
                                               error=error_sub)

           END DO

           ! CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
           !             1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
           !             I_tmp2(1,1),SIZE(I_tmp2,1),&
           !             1.0_dp,I_ab(1,1),SIZE(I_ab,1))
           ! L_local_col(offset-nsgfa(iset)+1:offset,i_counter)=I_ab(1:nsgfa(iset),1)
         END IF

         DEALLOCATE(I_tmp2)
         DEALLOCATE(I_ab)
         DEALLOCATE(pab)

        END DO

        force(ikind)%rho_elec(:,atom_a) =&
                    force(ikind)%rho_elec(:,atom_a) + force_a(:) + force_b
        IF(use_virial) THEN
          virial%pv_virial = virial%pv_virial + my_virial_a + my_virial_b
        END IF
      END DO

      DO i=1,SIZE(rs_v)
        CALL rs_grid_release(rs_v(i)%rs_grid, error=error_sub)
      END DO
      CALL timestop(handle3)
      ! here we are done with the 2 centers
     
      ! CALL cp_dbcsr_write_sparse_matrix(matrix_P_munu%matrix,4,12,qs_env,para_env_sub,&
      !                                   output_unit=unit_nr,error=error_sub)

      ! integrate the potential of the single gaussian and update 
      ! 3-center forces
      CALL timeset(routineN//"_int",handle3)
      CALL cp_dbcsr_set(mat_munu%matrix,0.0_dp,error=error_sub)
      CALL integrate_v_rspace(rho_r,hmat=mat_munu,pmat=matrix_P_munu,&
                              qs_env=qs_env,calculate_forces=.TRUE.,compute_tau=.FALSE.,gapw=.FALSE.,&
                              pw_env_external=pw_env_sub,&
                              task_list_external=task_list_sub, error=error_sub)
      CALL timestop(handle3) 

      ! update lagrangian
      CALL timeset(routineN//"_Lag",handle3)
      ! first contract mat_munu with the half back transformed Gamma_i_nu
      ! in order to update Lag_mu_i_1
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, matrix_P_inu%matrix, &
                             1.0_dp, Lag_mu_i_1%matrix, filter_eps=eps_filter, error=error_sub)

      ! transform first index of mat_munu and store the result into matrix_P_inu
      CALL cp_dbcsr_set(matrix_P_inu%matrix,0.0_dp,error=error_sub)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                             0.0_dp, matrix_P_inu%matrix, filter_eps=eps_filter, error=error_sub)

      ! contract the transformend matrix_P_inu with the untransformend Gamma_i_a
      ! in order to update Lag_nu_a_2
      CALL cp_dbcsr_multiply("N", "N", -1.0_dp, matrix_P_inu%matrix, G_P_ia(L_counter)%matrix, &
                             1.0_dp, Lag_nu_a_2%matrix, filter_eps=eps_filter, error=error_sub)

      ! release the actual gamma_P_ia
      CALL cp_dbcsr_release(G_P_ia(L_counter)%matrix,error=error_sub)
      DEALLOCATE(G_P_ia(L_counter)%matrix,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL timestop(handle3)

      IF(use_virial) THEN
        ! add the volume contribution to the virial due to
        ! the (P|Q) integrals, first we put the full gamme_PQ
        ! pseudo wave-function into grid in order to calculate the 
        ! hartree potential derivatives
        CALL timeset(routineN//"_Virial",handle3)
        wf_vector=0.0_dp
        wf_vector(:)=-2.0_dp*G_PQ_local(:,L_counter)
        CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set, &
                                    qs_kind_set,cell,dft_control,particle_set,pw_env_sub,&
                                    basis_type="RI_AUX",&
                                    external_vector=wf_vector,&
                                    error=error_sub)
        ! transfer to reciprocal space and calculate potential
        rho_r%pw%cr3d = psi_L%pw%cr3d
        CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
        CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
        ! update virial with volume term (first calculate hartree like energy (diagonal part of the virial))
        e_hartree=0.0_dp
        h_stress=0.0_dp
        e_hartree=pw_integral_ab(temp_pw_g%pw, pot_g%pw, error=error_sub)
        DO alpha=1, 3
          comp=0
          comp(alpha)=1
          CALL pw_copy(pot_g%pw, rho_g%pw, error=error_sub)
          CALL pw_derive(rho_g%pw, comp, error=error_sub)   
          h_stress(alpha,alpha)=-e_hartree
          DO beta=alpha, 3
            h_stress(alpha,beta)=h_stress(alpha,beta) &
                       -2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw, error=error_sub)/fourpi
            h_stress (beta,alpha)=h_stress(alpha,beta)
          END DO
        END DO
        virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env_sub%num_pe,dp)
        CALL timestop(handle3)
      END IF

      ! put the gamma density on grid
      CALL timeset(routineN//"_Gpot",handle3)
      CALL calculate_rho_elec(matrix_p=matrix_P_munu%matrix,&
                              rho=rho_r,&
                              rho_gspace=rho_g,&
                              total_rho=total_rho,&
                              task_list_external=task_list_sub,&
                              pw_env_external=pw_env_sub,&
                              ks_env=ks_env,error=error_sub)
      ! calculate associated hartree potential
      ! CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
      CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
      CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
      CALL timestop(handle3)

      IF(use_virial) THEN
        ! add the volume contribution to the virial coming from 
        ! the 3-center integrals (mu nu|P) 
        CALL timeset(routineN//"_Virial",handle3)
        e_hartree=0.0_dp
        h_stress=0.0_dp
        e_hartree=pw_integral_ab(temp_pw_g%pw, pot_g%pw, error=error_sub)
        DO alpha=1, 3
          comp=0
          comp(alpha)=1
          CALL pw_copy(pot_g%pw, rho_g%pw, error=error_sub)
          CALL pw_derive(rho_g%pw, comp, error=error_sub)
          h_stress(alpha,alpha)=-e_hartree
          DO beta=alpha, 3
            h_stress(alpha,beta)=h_stress(alpha,beta) &
                       -2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw, error=error_sub)/fourpi
            h_stress (beta,alpha)=h_stress(alpha,beta)
          END DO
        END DO
        virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env_sub%num_pe,dp)
        CALL timestop(handle3)
      END IF

      ! integrate potential with auxiliary basis function derivatives
      NULLIFY(rs_v)
      NULLIFY(rs_descs)
      CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v, error=error_sub)
      DO i=1,SIZE(rs_v)
        CALL rs_grid_retain(rs_v(i)%rs_grid,error=error_sub)
      END DO
      CALL potential_pw2rs(rs_v,rho_r,pw_env_sub,error_sub)

      offset=0
      DO iatom=1, natom
        ikind=kind_of(iatom)
        atom_a = atom_of_kind(iatom)
        CALL get_qs_kind(qs_kind=qs_kind_set(ikind),basis_set=basis_set_a,&
                         basis_type="RI_AUX",error=error)

        first_sgfa   =>  basis_set_a%first_sgf
        la_max       =>  basis_set_a%lmax
        la_min       =>  basis_set_a%lmin
        npgfa        =>  basis_set_a%npgf
        nseta        =   basis_set_a%nset
        nsgfa        =>  basis_set_a%nsgf_set
        rpgfa        =>  basis_set_a%pgf_radius
        set_radius_a =>  basis_set_a%set_radius
        sphi_a       =>  basis_set_a%sphi
        zeta         =>  basis_set_a%zet

        ra(:) = pbc(particle_set(iatom)%r,cell)
        rab=0.0_dp
        rab2=0.0_dp

        force_a(:) = 0.0_dp
        force_b(:) = 0.0_dp
        IF (use_virial) THEN
          my_virial_a = 0.0_dp
          my_virial_b = 0.0_dp
        END IF

        DO iset=1, nseta
         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         ALLOCATE(I_tmp2(ncoa,1))
         I_tmp2=0.0_dp
         ALLOCATE(I_ab(nsgfa(iset),1))
         I_ab=0.0_dp
         ALLOCATE(pab(ncoa,1))
         pab=0.0_dp

         skip_shell=.TRUE.
         DO iorb=1, nsgfa(iset)
           IF(iorb+offset==LLL) THEN
             I_ab(iorb,1)=1.0_dp
             skip_shell=.FALSE.
           END IF
         END DO  

         IF(skip_shell) THEN
           offset=offset+nsgfa(iset)
           DEALLOCATE(I_tmp2)
           DEALLOCATE(I_ab)
           DEALLOCATE(pab)
           CYCLE
         END IF

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                     1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                     I_ab(1,1),SIZE(I_ab,1),&
                     0.0_dp,pab(1,1),SIZE(pab,1))
         I_ab=0.0_dp

         igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,MINVAL(zeta(:,iset)))
         map_it_here=.FALSE.
         IF (.NOT. ALL (rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
            DO dir = 1,3
                  ! bounds of local grid (i.e. removing the 'wings'), if periodic
                  tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                  tp(dir) = MODULO ( tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir) )
                  IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                     lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir ) + rs_v(igrid_level)%rs_grid%desc%border
                     ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir ) - rs_v(igrid_level)%rs_grid%desc%border
                  ELSE
                     lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir )
                     ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir )
                  ENDIF
                  ! distributed grid, only map if it is local to the grid
                  location(dir)=tp(dir)+rs_v(igrid_level)%rs_grid%desc%lb(dir)
            ENDDO
            IF  (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                 lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                 lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
               map_it_here=.TRUE.
            ENDIF
         ELSE
            ! not distributed, just a round-robin distribution over the full set of CPUs
            IF (MODULO(offset,para_env_sub%num_pe)==para_env_sub%mepos) map_it_here=.TRUE.
         ENDIF

         offset=offset+nsgfa(iset)

         IF (map_it_here) THEN
           DO ipgf=1, npgfa(iset)
             na1=(ipgf - 1)*ncoset(la_max(iset)) + 1
             na2=ipgf*ncoset(la_max(iset))
             igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,zeta(ipgf,iset))

             CALL integrate_pgf_product_rspace(la_max=la_max(iset),zeta=zeta(ipgf,iset)/2.0_dp,la_min=la_min(iset),&
                                               lb_max=0,zetb=zeta(ipgf,iset)/2.0_dp,lb_min=0,&
                                               ra=ra,rab=rab,rab2=rab2,&
                                               rsgrid=rs_v(igrid_level)%rs_grid,&
                                               cell=cell,&
                                               cube_info=pw_env_sub%cube_info(igrid_level),&
                                               hab=I_tmp2,&
                                               pab=pab,&
                                               o1=na1-1,&
                                               o2=0,&
                                               map_consistent=.TRUE.,&
                                               eps_gvg_rspace=dft_control%qs_control%eps_gvg_rspace,&
                                               calculate_forces=.TRUE.,&
                                               force_a=force_a, force_b=force_b, &
                                               use_virial=use_virial,my_virial_a=my_virial_a,my_virial_b=my_virial_b,&
                                               error=error_sub)

           END DO

           ! CALL dgemm("T","N".0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
           !             I_tmp2(1,1),SIZE(I_tmp2,1),&
           !             1.0_dp,I_ab(1,1),SIZE(I_ab,1))
           ! L_local_col(offset-nsgfa(iset)+1:offset,i_counter)=I_ab(1:nsgfa(iset),1)
         END IF

         DEALLOCATE(I_tmp2)
         DEALLOCATE(I_ab)
         DEALLOCATE(pab)

        END DO

        force(ikind)%rho_elec(:,atom_a) =&
                      force(ikind)%rho_elec(:,atom_a) + force_a(:) + force_b(:)
        IF(use_virial) THEN
          virial%pv_virial = virial%pv_virial + my_virial_a + my_virial_b
        END IF
      END DO

      DO i=1,SIZE(rs_v)
        CALL rs_grid_release(rs_v(i)%rs_grid, error=error_sub)
      END DO

    END DO
    CALL timestop(handle2)

    DEALLOCATE(G_PQ_local,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(wf_vector,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_L%pw,error=error_sub)

    IF(use_virial) THEN
      CALL pw_pool_give_back_pw(auxbas_pw_pool,temp_pw_g%pw,error=error_sub)
      DO i=1, 3
        CALL pw_pool_give_back_pw(auxbas_pw_pool,dvg(i)%pw,error=error_sub)
      END DO
    END IF

    CALL cp_dbcsr_release(matrix_P_inu%matrix,error=error_sub)
    DEALLOCATE(matrix_P_inu%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_release(matrix_P_munu_nosym%matrix,error=error_sub)
    DEALLOCATE(matrix_P_munu_nosym%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_dbcsr_release(matrix_P_munu%matrix,error=error_sub)
    DEALLOCATE(matrix_P_munu%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! release the full gamma_P_ia structure
    DEALLOCATE(G_P_ia,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(kind_of)
    DEALLOCATE(atom_of_kind)

    ! move the forces in the correct place
    DO ikind=1,SIZE(force)
      force(ikind)%mp2_non_sep(:,:) = force(ikind)%rho_elec(:,:) 
      force(ikind)%rho_elec(:,:) = 0.0_dp
    END DO

    ! Now we move from the local matrices to the global ones
    ! defined over all MPI tasks
    ! Start with moving from the DBCSR to FM for the lagrangians
    NULLIFY(L1_mu_i, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env_sub,context=blacs_env_sub, &
                             nrow_global=dimen,ncol_global=homo,error=error_sub)
    CALL cp_fm_create(L1_mu_i, fm_struct_tmp,name="Lag_mu_i",error=error_sub)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error_sub)
    CALL cp_fm_set_all(L1_mu_i,0.0_dp,error=error_sub)
    CALL copy_dbcsr_to_fm(matrix=Lag_mu_i_1%matrix, fm=L1_mu_i, error=error_sub)

    ! release Lag_mu_i_1
    CALL cp_dbcsr_release(Lag_mu_i_1%matrix,error=error_sub)
    DEALLOCATE(Lag_mu_i_1%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(L2_nu_a, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env_sub,context=blacs_env_sub, &
                             nrow_global=dimen,ncol_global=virtual,error=error_sub)
    CALL cp_fm_create(L2_nu_a, fm_struct_tmp,name="Lag_nu_a",error=error_sub)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error_sub)
    CALL cp_fm_set_all(L2_nu_a,0.0_dp,error=error_sub)
    CALL copy_dbcsr_to_fm(matrix=Lag_nu_a_2%matrix, fm=L2_nu_a, error=error_sub)

    ! release Lag_nu_a_2
    CALL cp_dbcsr_release(Lag_nu_a_2%matrix,error=error_sub)
    DEALLOCATE(Lag_nu_a_2%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL create_W_P(qs_env,mp2_env,mo_coeff,homo,virtual,dimen,para_env,para_env_sub,&
                    Eigenval,L1_mu_i,L2_nu_a,error,error_sub)

    CALL timestop(handle)

  END SUBROUTINE calc_ri_mp2_nonsep

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param mo_coeff ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param Eigenval ...
!> \param L1_mu_i ...
!> \param L2_nu_a ...
!> \param error ...
!> \param error_sub ...
! *****************************************************************************
  SUBROUTINE create_W_P(qs_env,mp2_env,mo_coeff,homo,virtual,dimen,para_env,para_env_sub,&
                        Eigenval,L1_mu_i,L2_nu_a,error,error_sub)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER                                  :: homo, virtual, dimen
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    TYPE(cp_fm_type), POINTER                :: L1_mu_i, L2_nu_a
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_W_P', &
      routineP = moduleN//':'//routineN

    INTEGER :: color_exchange, comm_exchange, dummy_proc, handle, handle2, &
      handle3, i, i_global, i_local, iiB, iii, iproc, itmp(2), j_global, &
      j_local, jjB, max_col_size, max_row_size, my_B_virtual_end, &
      my_B_virtual_start, mypcol, mypcol_1i, mypcol_2a, myprow, myprow_1i, &
      myprow_2a, ncol_block, ncol_block_1i, ncol_block_2a, ncol_local, &
      ncol_local_1i, ncol_local_2a, npcol, npcol_1i, npcol_2a, nprow, &
      nprow_1i, nprow_2a, nrow_block, nrow_block_1i, nrow_block_2a, &
      nrow_local, nrow_local_1i, nrow_local_2a, number_of_rec, &
      number_of_send, proc_receive, proc_receive_static, proc_send, &
      proc_send_ex
    INTEGER :: proc_send_static, proc_send_sub, proc_shift, rec_col_size, &
      rec_counter, rec_row_size, send_col_size, send_counter, send_pcol, &
      send_prow, send_row_size, size_rec_buffer, size_send_buffer, stat
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iii_vet, map_rec_size, &
      map_send_size, pos_info, pos_info_ex, proc_2_send_pos, proc_map_ex, &
      req_send, sub_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: grid_2_mepos, mepos_2_grid, &
                                                mepos_2_grid_1i, &
                                                mepos_2_grid_2a, sizes, &
                                                sizes_1i, sizes_2a
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: col_indeces_info_1i, &
      col_indeces_info_2a, row_indeces_info_1i, row_indeces_info_2a
    INTEGER, DIMENSION(:), POINTER :: col_indices, col_indices_1i, &
      col_indices_2a, row_indices, row_indices_1i, row_indices_2a
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: ab_rec, ab_send, mat_deb, &
                                                mat_rec, mat_send
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_P_ij, L_mu_q
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send
    TYPE(integ_mat_buffer_type_2D), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_cyclic

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! create the globally distributed mixed lagrangian
    NULLIFY(blacs_env)
    CALL get_qs_env(qs_env, blacs_env=blacs_env, error=error)

    NULLIFY(L_mu_q, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=dimen,error=error)
    CALL cp_fm_create(L_mu_q, fm_struct_tmp, name="Lag_mu_q",error=error_sub)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(L_mu_q, 0.0_dp,error=error_sub)

    ! create all information array
    ALLOCATE(pos_info(0:para_env%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    pos_info=0
    pos_info(para_env%mepos)=para_env_sub%mepos
    CALL mp_sum(pos_info,para_env%group)

    ALLOCATE(sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    sub_proc_map=0
    DO i=0,para_env_sub%num_pe-1
      sub_proc_map(i)=i
      sub_proc_map(-i-1)=para_env_sub%num_pe-i-1
      sub_proc_map(para_env_sub%num_pe+i)=i
    END DO
 
    ! get matrix information for the global
    CALL cp_fm_get_info(matrix=L_mu_q,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block,&
                        error=error)
    myprow=L_mu_q%matrix_struct%context%mepos(1)
    mypcol=L_mu_q%matrix_struct%context%mepos(2)
    nprow =L_mu_q%matrix_struct%context%num_pe(1)
    npcol =L_mu_q%matrix_struct%context%num_pe(2)

    ALLOCATE(grid_2_mepos(0:nprow-1,0:npcol-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    grid_2_mepos=0
    grid_2_mepos(myprow,mypcol)=para_env%mepos
    CALL mp_sum(grid_2_mepos,para_env%group)

    ! get matrix information for L1_mu_i
    CALL cp_fm_get_info(matrix=L1_mu_i,&
                        nrow_local=nrow_local_1i,&
                        ncol_local=ncol_local_1i,&
                        row_indices=row_indices_1i,&
                        col_indices=col_indices_1i,&
                        nrow_block=nrow_block_1i,&
                        ncol_block=ncol_block_1i,&
                        error=error_sub)
    myprow_1i=L1_mu_i%matrix_struct%context%mepos(1)
    mypcol_1i=L1_mu_i%matrix_struct%context%mepos(2)
    nprow_1i =L1_mu_i%matrix_struct%context%num_pe(1)
    npcol_1i =L1_mu_i%matrix_struct%context%num_pe(2)
    
    ALLOCATE(mepos_2_grid_1i(0:para_env_sub%num_pe-1,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    mepos_2_grid_1i=0
    mepos_2_grid_1i(para_env_sub%mepos,1)=myprow_1i
    mepos_2_grid_1i(para_env_sub%mepos,2)=mypcol_1i
    CALL mp_sum(mepos_2_grid_1i,para_env_sub%group)

    ALLOCATE(sizes_1i(2,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sizes_1i=0
    sizes_1i(1,para_env_sub%mepos)=nrow_local_1i
    sizes_1i(2,para_env_sub%mepos)=ncol_local_1i
    CALL mp_sum(sizes_1i,para_env_sub%group)

    ! get matrix information for L2_nu_a
    CALL cp_fm_get_info(matrix=L2_nu_a,&
                        nrow_local=nrow_local_2a,&
                        ncol_local=ncol_local_2a,&
                        row_indices=row_indices_2a,&
                        col_indices=col_indices_2a,&
                        nrow_block=nrow_block_2a,&
                        ncol_block=ncol_block_2a,&
                        error=error_sub)
    myprow_2a=L2_nu_a%matrix_struct%context%mepos(1)
    mypcol_2a=L2_nu_a%matrix_struct%context%mepos(2)
    nprow_2a =L2_nu_a%matrix_struct%context%num_pe(1)
    npcol_2a =L2_nu_a%matrix_struct%context%num_pe(2)

    ALLOCATE(mepos_2_grid_2a(0:para_env_sub%num_pe-1,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    mepos_2_grid_2a=0
    mepos_2_grid_2a(para_env_sub%mepos,1)=myprow_2a
    mepos_2_grid_2a(para_env_sub%mepos,2)=mypcol_2a
    CALL mp_sum(mepos_2_grid_2a,para_env_sub%group)

    ALLOCATE(sizes_2a(2,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sizes_2a=0
    sizes_2a(1,para_env_sub%mepos)=nrow_local_2a
    sizes_2a(2,para_env_sub%mepos)=ncol_local_2a
    CALL mp_sum(sizes_2a,para_env_sub%group)

    ! Here we perform a ring communication scheme taking into account
    ! for the sub-group distribution of the source matrices.
    ! as a first step we need to redistribute the data within
    ! the subgroup.
    ! In order to do so we have to allocate the structure
    ! that will hold the local data involved in communication, this
    ! structure will be the same for processes in different subgroups
    ! sharing the same position in the subgroup.
    ! -1) create the exchange para_env
    color_exchange=para_env_sub%mepos
    CALL mp_comm_split_direct(para_env%group,comm_exchange,color_exchange)
    NULLIFY(para_env_exchange)
    CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error)
    ! crate the proc maps exchange and info
    ALLOCATE(proc_map_ex(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
    DO i=0,para_env_exchange%num_pe-1
      proc_map_ex(i)=i
      proc_map_ex(-i-1)=para_env_exchange%num_pe-i-1
      proc_map_ex(para_env_exchange%num_pe+i)=i
    END DO
    ALLOCATE(pos_info_ex(0:para_env%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    pos_info_ex=0
    pos_info_ex(para_env%mepos)=para_env_exchange%mepos
    CALL mp_sum(pos_info_ex,para_env%group)
    ALLOCATE(sizes(2,0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    sizes=0
    sizes(1,para_env_exchange%mepos)=nrow_local
    sizes(2,para_env_exchange%mepos)=ncol_local
    CALL mp_sum(sizes,para_env_exchange%group)

    ! 0) store some info about indeces of the fm matrices (subgroup)
    CALL timeset(routineN//"_inx",handle2)
    ! matrix L1_mu_i
    max_row_size=MAXVAL(sizes_1i(1,:))
    max_col_size=MAXVAL(sizes_1i(2,:))
    ALLOCATE(row_indeces_info_1i(2,max_row_size,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    ALLOCATE(col_indeces_info_1i(2,max_col_size,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    row_indeces_info_1i=0
    col_indeces_info_1i=0
    ! row
    DO iiB=1, nrow_local_1i
      i_global=row_indices_1i(iiB)
      send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
                              L_mu_q%matrix_struct%first_p_pos(1),nprow)
      i_local=cp_fm_indxg2l(i_global,nrow_block,dummy_proc,&
                            L_mu_q%matrix_struct%first_p_pos(1),nprow)
      row_indeces_info_1i(1,iiB,para_env_sub%mepos)=send_prow
      row_indeces_info_1i(2,iiB,para_env_sub%mepos)=i_local
    END DO
    ! col
    DO jjB=1, ncol_local_1i
      j_global=col_indices_1i(jjB)
      send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
                              L_mu_q%matrix_struct%first_p_pos(2),npcol)
      j_local=cp_fm_indxg2l(j_global,ncol_block,dummy_proc,&
                            L_mu_q%matrix_struct%first_p_pos(2),npcol)
      col_indeces_info_1i(1,jjB,para_env_sub%mepos)=send_pcol
      col_indeces_info_1i(2,jjB,para_env_sub%mepos)=j_local
    END DO
    CALL mp_sum(row_indeces_info_1i,para_env_sub%group)
    CALL mp_sum(col_indeces_info_1i,para_env_sub%group)

    ! matrix L2_nu_a
    max_row_size=MAXVAL(sizes_2a(1,:))
    max_col_size=MAXVAL(sizes_2a(2,:))
    ALLOCATE(row_indeces_info_2a(2,max_row_size,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    ALLOCATE(col_indeces_info_2a(2,max_col_size,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    row_indeces_info_2a=0
    col_indeces_info_2a=0
    ! row
    DO iiB=1, nrow_local_2a
      i_global=row_indices_2a(iiB)
      send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
                              L_mu_q%matrix_struct%first_p_pos(1),nprow)
      i_local=cp_fm_indxg2l(i_global,nrow_block,dummy_proc,&
                            L_mu_q%matrix_struct%first_p_pos(1),nprow)
      row_indeces_info_2a(1,iiB,para_env_sub%mepos)=send_prow
      row_indeces_info_2a(2,iiB,para_env_sub%mepos)=i_local
    END DO
    ! col
    DO jjB=1, ncol_local_2a
      j_global=col_indices_2a(jjB)+homo
      send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
                              L_mu_q%matrix_struct%first_p_pos(2),npcol)
      j_local=cp_fm_indxg2l(j_global,ncol_block,dummy_proc,&
                            L_mu_q%matrix_struct%first_p_pos(2),npcol)
      col_indeces_info_2a(1,jjB,para_env_sub%mepos)=send_pcol
      col_indeces_info_2a(2,jjB,para_env_sub%mepos)=j_local
    END DO
    CALL mp_sum(row_indeces_info_2a,para_env_sub%group)
    CALL mp_sum(col_indeces_info_2a,para_env_sub%group)
    CALL timestop(handle2)

    ! 1) define the map for sending data in the subgroup starting with L1_mu_i
    CALL timeset(routineN//"_subinfo",handle2)
    ALLOCATE(map_send_size(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_send_size=0
    DO jjB=1, ncol_local_1i
      ! j_global=col_indices_1i(jjB)
      ! send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
      !                         L_mu_q%matrix_struct%first_p_pos(2),npcol)
      send_pcol=col_indeces_info_1i(1,jjB,para_env_sub%mepos)
      DO iiB=1, nrow_local_1i
        ! i_global=row_indices_1i(iiB)
        ! send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
        !                         L_mu_q%matrix_struct%first_p_pos(1),nprow)
        send_prow=row_indeces_info_1i(1,iiB,para_env_sub%mepos)
        proc_send=grid_2_mepos(send_prow,send_pcol)
        proc_send_sub=pos_info(proc_send)
        map_send_size(proc_send_sub)=map_send_size(proc_send_sub)+1
      END DO
    END DO
    ! and the same for L2_nu_a
    DO jjB=1, ncol_local_2a
      ! j_global=col_indices_2a(jjB)+homo
      ! send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
      !                         L_mu_q%matrix_struct%first_p_pos(2),npcol)
      send_pcol=col_indeces_info_2a(1,jjB,para_env_sub%mepos)
      DO iiB=1, nrow_local_2a
        ! i_global=row_indices_2a(iiB)
        ! send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
        !                         L_mu_q%matrix_struct%first_p_pos(1),nprow)
        send_prow=row_indeces_info_2a(1,iiB,para_env_sub%mepos)
        proc_send=grid_2_mepos(send_prow,send_pcol)
        proc_send_sub=pos_info(proc_send)
        map_send_size(proc_send_sub)=map_send_size(proc_send_sub)+1
      END DO
    END DO
    ! and exchange data in order to create map_rec_size
    ALLOCATE(map_rec_size(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_rec_size=0
    CALL mp_alltoall(map_send_size,map_rec_size,1,para_env_sub%group)
    CALL timestop(handle2)
 
    ! 2) reorder data in sending buffer
    CALL timeset(routineN//"_sub_Bsend",handle2)
    ! count the number of messages (include myself)
    number_of_send=0
    DO proc_shift=0, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        number_of_send=number_of_send+1
      END IF
    END DO
    ! allocate the structure that will hold the messages to be sent
    ALLOCATE(buffer_send(number_of_send),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    send_counter=0
    ALLOCATE(proc_2_send_pos(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    proc_2_send_pos=0
    DO proc_shift=0, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      size_send_buffer=map_send_size(proc_send)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        ! allocate the sending buffer (msg)
        ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        buffer_send(send_counter)%msg=0.0_dp
        buffer_send(send_counter)%proc=proc_send
        proc_2_send_pos(proc_send)=send_counter
      END IF
    END DO
    ! loop over the locally held data and fill the buffer_send
    ! for doing that we need an array that keep track if the 
    ! sequential increase of the index for each message
    ALLOCATE(iii_vet(number_of_send),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    iii_vet=0
    DO jjB=1, ncol_local_1i
      ! j_global=col_indices_1i(jjB)
      ! send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
      !                         L_mu_q%matrix_struct%first_p_pos(2),npcol)
      send_pcol=col_indeces_info_1i(1,jjB,para_env_sub%mepos)
      DO iiB=1, nrow_local_1i
        ! i_global=row_indices_1i(iiB)
        ! send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
        !                         L_mu_q%matrix_struct%first_p_pos(1),nprow)
        send_prow=row_indeces_info_1i(1,iiB,para_env_sub%mepos)
        proc_send=grid_2_mepos(send_prow,send_pcol)
        proc_send_sub=pos_info(proc_send)
        send_counter=proc_2_send_pos(proc_send_sub)
        iii_vet(send_counter)=iii_vet(send_counter)+1
        iii=iii_vet(send_counter)
        buffer_send(send_counter)%msg(iii)=L1_mu_i%local_data(iiB,jjB)
      END DO
    END DO
    ! release the local data of L1_mu_i
    DEALLOCATE(L1_mu_i%local_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! and the same for L2_nu_a
    DO jjB=1, ncol_local_2a
      ! j_global=col_indices_2a(jjB)+homo
      ! send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
      !                         L_mu_q%matrix_struct%first_p_pos(2),npcol)
      send_pcol=col_indeces_info_2a(1,jjB,para_env_sub%mepos)
      DO iiB=1, nrow_local_2a
        ! i_global=row_indices_2a(iiB)
        ! send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
        !                         L_mu_q%matrix_struct%first_p_pos(1),nprow)
        send_prow=row_indeces_info_2a(1,iiB,para_env_sub%mepos)
        proc_send=grid_2_mepos(send_prow,send_pcol)
        proc_send_sub=pos_info(proc_send)
        send_counter=proc_2_send_pos(proc_send_sub)
        iii_vet(send_counter)=iii_vet(send_counter)+1
        iii=iii_vet(send_counter)
        buffer_send(send_counter)%msg(iii)=L2_nu_a%local_data(iiB,jjB)
      END DO
    END DO
    DEALLOCATE(L2_nu_a%local_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(proc_2_send_pos)
    DEALLOCATE(iii_vet)
    CALL timestop(handle2)

    ! 3) create the buffer for receive, post the message with irecv
    !    and send the messages with mp_isend
    CALL timeset(routineN//"_sub_isendrecv",handle2)
    ! count the number of messages to be received
    number_of_rec=0
    DO proc_shift=0, para_env_sub%num_pe-1
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
      IF(map_rec_size(proc_receive)>0) THEN
        number_of_rec=number_of_rec+1
      END IF
    END DO
    ALLOCATE(buffer_rec(number_of_rec),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
    rec_counter=0
    DO proc_shift=0, para_env_sub%num_pe-1
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)
      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1
        ! prepare the buffer for receive
        ALLOCATE(buffer_rec(rec_counter)%msg(size_rec_buffer),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        buffer_rec(rec_counter)%msg=0.0_dp
        buffer_rec(rec_counter)%proc=proc_receive
        ! post the message to be received (not need to send to myself)
        IF(proc_receive/=para_env_sub%mepos) THEN
          CALL mp_irecv(buffer_rec(rec_counter)%msg,proc_receive,para_env_sub%group,&
                        buffer_rec(rec_counter)%msg_req)
        END IF
      END IF
    END DO
    ! send messages
    ALLOCATE(req_send(number_of_send),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    req_send=mp_request_null
    send_counter=0
    DO proc_shift=0, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        IF(proc_send==para_env_sub%mepos) THEN
          buffer_rec(send_counter)%msg=buffer_send(send_counter)%msg
        ELSE
          CALL mp_isend(buffer_send(send_counter)%msg,proc_send,para_env_sub%group,&
                        buffer_send(send_counter)%msg_req)
          req_send(send_counter)=buffer_send(send_counter)%msg_req
        END IF
      END IF
    END DO
    DEALLOCATE(map_send_size)
    CALL timestop(handle2)

    ! 4) (if memory is a problem we should move this part after point 5)
    !    Here we create the new buffer for cyclic(ring) communication and
    !    we fill it with the data received from the other member of the 
    !    subgroup
    CALL timeset(routineN//"_Bcyclic",handle2)
    ! first allocata new structure
    ALLOCATE(buffer_cyclic(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iproc=0, para_env_exchange%num_pe-1
      rec_row_size=sizes(1,iproc)
      rec_col_size=sizes(2,iproc)
      ALLOCATE(buffer_cyclic(iproc)%msg(rec_row_size,rec_col_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      buffer_cyclic(iproc)%msg=0.0_dp
    END DO
    ! now collect data from other member of the subgroup and fill 
    ! buffer_cyclic
    rec_counter=0
    DO proc_shift=0, para_env_sub%num_pe-1
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)
      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1

        ! wait for the message       
        IF(proc_receive/=para_env_sub%mepos) CALL mp_wait(buffer_rec(rec_counter)%msg_req)

        CALL timeset(routineN//"_fill",handle3)
        iii=0
        DO jjB=1, sizes_1i(2,proc_receive)
          ! j_global=cp_fm_indxl2g(jjB,ncol_block_1i,mepos_2_grid_1i(proc_receive,2),&
          !                        L1_mu_i%matrix_struct%first_p_pos(2),npcol_1i)
          ! send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
          !                         L_mu_q%matrix_struct%first_p_pos(2),npcol)
          ! j_local=cp_fm_indxg2l(j_global,ncol_block,dummy_proc,&
          !                       L_mu_q%matrix_struct%first_p_pos(2),npcol)
          send_pcol = col_indeces_info_1i(1,jjB,proc_receive)
          j_local   = col_indeces_info_1i(2,jjB,proc_receive)
          DO iiB=1, sizes_1i(1,proc_receive)
            ! i_global=cp_fm_indxl2g(iiB,nrow_block_1i,mepos_2_grid_1i(proc_receive,1),&
            !                        L1_mu_i%matrix_struct%first_p_pos(1),nprow_1i)
            ! send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
            !                         L_mu_q%matrix_struct%first_p_pos(1),nprow)
            send_prow=row_indeces_info_1i(1,iiB,proc_receive)
            proc_send=grid_2_mepos(send_prow,send_pcol)
            proc_send_sub=pos_info(proc_send)
            IF(proc_send_sub/=para_env_sub%mepos) CYCLE
            iii=iii+1
            ! i_local=cp_fm_indxg2l(i_global,nrow_block,dummy_proc,&
            !                       L_mu_q%matrix_struct%first_p_pos(1),nprow)
            i_local=row_indeces_info_1i(2,iiB,proc_receive)
            proc_send_ex=pos_info_ex(proc_send)
            buffer_cyclic(proc_send_ex)%msg(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
          END DO
        END DO
        ! and the same for L2_nu_a
        DO jjB=1, sizes_2a(2,proc_receive)
          ! j_global=cp_fm_indxl2g(jjB,ncol_block_2a,mepos_2_grid_2a(proc_receive,2),&
          !                        L2_nu_a%matrix_struct%first_p_pos(2),npcol_2a)+homo
          ! send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
          !                         L_mu_q%matrix_struct%first_p_pos(2),npcol)
          ! j_local=cp_fm_indxg2l(j_global,ncol_block,dummy_proc,&
          !                       L_mu_q%matrix_struct%first_p_pos(2),npcol)
          send_pcol = col_indeces_info_2a(1,jjB,proc_receive)
          j_local   = col_indeces_info_2a(2,jjB,proc_receive)
          DO iiB=1, sizes_2a(1,proc_receive)
            ! i_global=cp_fm_indxl2g(iiB,nrow_block_2a,mepos_2_grid_2a(proc_receive,1),&
            !                        L2_nu_a%matrix_struct%first_p_pos(1),nprow_2a)
            ! send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
            !                         L_mu_q%matrix_struct%first_p_pos(1),nprow)
            send_prow=row_indeces_info_2a(1,iiB,proc_receive)
            proc_send=grid_2_mepos(send_prow,send_pcol)
            proc_send_sub=pos_info(proc_send)
            IF(proc_send_sub/=para_env_sub%mepos) CYCLE
            iii=iii+1
            ! i_local=cp_fm_indxg2l(i_global,nrow_block,dummy_proc,&
            !                       L_mu_q%matrix_struct%first_p_pos(1),nprow)
            i_local=row_indeces_info_2a(2,iiB,proc_receive)
            proc_send_ex=pos_info_ex(proc_send)
            buffer_cyclic(proc_send_ex)%msg(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
          END DO
        END DO
        CALL timestop(handle3)

        ! deallocate the received message
        DEALLOCATE(buffer_rec(rec_counter)%msg)
      END IF
    END DO
    DEALLOCATE(row_indeces_info_1i)
    DEALLOCATE(col_indeces_info_1i)
    DEALLOCATE(row_indeces_info_2a)
    DEALLOCATE(col_indeces_info_2a)
    DEALLOCATE(buffer_rec)
    DEALLOCATE(map_rec_size)
    CALL timestop(handle2)

    ! 5)  Wait for all messeges to be sent in the subgroup
    CALL timeset(routineN//"_sub_waitall",handle2)
    CALL mp_waitall(req_send(:))
    DO send_counter=1, number_of_send
      DEALLOCATE(buffer_send(send_counter)%msg)
    END DO
    DEALLOCATE(buffer_send)
    DEALLOCATE(req_send)
    CALL timestop(handle2)

    ! 6) Start with ring communication
    CALL timeset(routineN//"_ring",handle2)
    proc_send_static=proc_map_ex(para_env_exchange%mepos+1)
    proc_receive_static=proc_map_ex(para_env_exchange%mepos-1)
    max_row_size=MAXVAL(sizes(1,:))
    max_col_size=MAXVAL(sizes(2,:))
    ALLOCATE(mat_send(max_row_size,max_col_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(mat_rec(max_row_size,max_col_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    mat_send=0.0_dp
    mat_send(1:nrow_local,1:ncol_local)=buffer_cyclic(para_env_exchange%mepos)%msg(:,:)
    DEALLOCATE(buffer_cyclic(para_env_exchange%mepos)%msg)
    DO proc_shift=1, para_env_exchange%num_pe-1
      proc_send=proc_map_ex(para_env_exchange%mepos+proc_shift)
      proc_receive=proc_map_ex(para_env_exchange%mepos-proc_shift)

      rec_row_size=sizes(1,proc_receive)
      rec_col_size=sizes(2,proc_receive)

      mat_rec=0.0_dp
      CALL  mp_sendrecv(mat_send,proc_send_static,&
                        mat_rec,proc_receive_static,&
                        para_env_exchange%group)

      mat_send=0.0_dp
      mat_send(1:rec_row_size,1:rec_col_size)=mat_rec(1:rec_row_size,1:rec_col_size)+&
                                              buffer_cyclic(proc_receive)%msg(:,:)

      DEALLOCATE(buffer_cyclic(proc_receive)%msg)
    END DO
    ! and finaly
    CALL  mp_sendrecv(mat_send,proc_send_static,&
                      mat_rec,proc_receive_static,&
                      para_env_exchange%group)
    L_mu_q%local_data(1:nrow_local,1:ncol_local)=mat_rec(1:nrow_local,1:ncol_local)
    DEALLOCATE(buffer_cyclic)
    DEALLOCATE(mat_send)
    DEALLOCATE(mat_rec)
    CALL timestop(handle2)

    DEALLOCATE(proc_map_ex)
    ! release para_env_exchange
    CALL cp_para_env_release(para_env_exchange,error=error_sub)

    CALL cp_fm_release(L1_mu_i,error=error_sub)
    CALL cp_fm_release(L2_nu_a,error=error_sub)
    DEALLOCATE(pos_info_ex)
    DEALLOCATE(grid_2_mepos)
    DEALLOCATE(mepos_2_grid_1i)
    DEALLOCATE(mepos_2_grid_2a)
    DEALLOCATE(sizes)
    DEALLOCATE(sizes_1i)
    DEALLOCATE(sizes_2a)

    ! update the P_ij block of P_MP2 with the 
    ! non-singular ij pairs
    CALL timeset(routineN//"_Pij",handle2)
    NULLIFY(fm_P_ij, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=homo,ncol_global=homo,error=error)
    CALL cp_fm_create(fm_P_ij, fm_struct_tmp, name="fm_P_ij",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(fm_P_ij, 0.0_dp,error=error)

    CALL cp_gemm('T','N',homo,homo,dimen,1.0_dp,&
                    mo_coeff,L_mu_q,0.0_dp,fm_P_ij,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    ! don't know if it is better to transpose (for communication reasons)
    ! or just recompute the transposed matrix
    CALL cp_gemm('T','N',homo,homo,dimen,-2.0_dp,&
                    L_mu_q,mo_coeff,2.0_dp,fm_P_ij,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    ! we have it, update P_ij local
    CALL cp_fm_get_info(matrix=fm_P_ij,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error)   
    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      DO iiB=1, nrow_local
        i_global=row_indices(iiB)
        ! diagonal elements alreasy updated
        IF(j_global==i_global) CYCLE
        ! check if the given element is above the threshold
        IF(ABS(Eigenval(j_global)-Eigenval(i_global))<mp2_env%ri_mp2%eps_canonical) CYCLE
        mp2_env%ri_grad%P_ij(i_global,j_global)=&
                      fm_P_ij%local_data(iiB,jjB)/(Eigenval(j_global)-Eigenval(i_global))
      END DO
    END DO
    ! release fm_P_ij
    CALL cp_fm_release(fm_P_ij,error=error)
    ! mp_sum it (we can avoid mp_sum, but for now let's keep it easy)
    CALL mp_sum(mp2_env%ri_grad%P_ij,para_env%group)
    CALL timestop(handle2)

    ! Now create and fill the P matrix (MO)
    ! FOR NOW WE ASSUME P_ab AND P_ij ARE REPLICATED OVER EACH MPI RANK
    CALL timeset(routineN//"_PMO",handle2)
    NULLIFY(mp2_env%ri_grad%P_mo, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=dimen,error=error)
    CALL cp_fm_create(mp2_env%ri_grad%P_mo, fm_struct_tmp, name="P_MP2_MO",error=error)
    CALL cp_fm_set_all(mp2_env%ri_grad%P_mo, 0.0_dp,error=error)

    ! start with the (easy) occ-occ block and locally held P_ab elements
    itmp=get_limit(virtual,para_env_sub%num_pe,para_env_sub%mepos)
    my_B_virtual_start=itmp(1)
    my_B_virtual_end=itmp(2)

    CALL cp_fm_get_info(matrix=mp2_env%ri_grad%P_mo,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block,&
                        error=error)
    myprow=mp2_env%ri_grad%P_mo%matrix_struct%context%mepos(1)
    mypcol=mp2_env%ri_grad%P_mo%matrix_struct%context%mepos(2)
    nprow =mp2_env%ri_grad%P_mo%matrix_struct%context%num_pe(1)
    npcol =mp2_env%ri_grad%P_mo%matrix_struct%context%num_pe(2)

    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      DO iiB=1, nrow_local
        i_global=row_indices(iiB)
        IF(i_global<=homo.AND.j_global<=homo) THEN
          mp2_env%ri_grad%P_mo%local_data(iiB,jjB)=mp2_env%ri_grad%P_ij(i_global,j_global)
        END IF
        IF((my_B_virtual_start<=i_global-homo.AND.i_global-homo<=my_B_virtual_end).AND.(j_global>homo)) THEN
          mp2_env%ri_grad%P_mo%local_data(iiB,jjB)=&
                             mp2_env%ri_grad%P_ab(i_global-homo-my_B_virtual_start+1,j_global-homo)
        END IF
      END DO
    END DO
    ! deallocate the local P_ij
    DEALLOCATE(mp2_env%ri_grad%P_ij,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
    ! send around the sub_group the local data and check if we 
    ! have to update our block with external elements
    ALLOCATE(mepos_2_grid(0:para_env_sub%num_pe-1,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    mepos_2_grid=0
    mepos_2_grid(para_env_sub%mepos,1)=myprow
    mepos_2_grid(para_env_sub%mepos,2)=mypcol
    CALL mp_sum(mepos_2_grid,para_env_sub%group)

    ALLOCATE(sizes(2,0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sizes=0
    sizes(1,para_env_sub%mepos)=nrow_local
    sizes(2,para_env_sub%mepos)=ncol_local
    CALL mp_sum(sizes,para_env_sub%group)

    ALLOCATE(ab_rec(nrow_local,ncol_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
    
      send_prow=mepos_2_grid(proc_send,1)
      send_pcol=mepos_2_grid(proc_send,2)

      send_row_size=sizes(1,proc_send)
      send_col_size=sizes(2,proc_send)

      ALLOCATE(ab_send(send_row_size,send_col_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ab_send=0.0_dp

      ! first loop over row since in this way we can cycle
      DO iiB=1, send_row_size
        i_global=cp_fm_indxl2g(iiB,nrow_block,send_prow,&
                               mp2_env%ri_grad%P_mo%matrix_struct%first_p_pos(1),nprow)
        IF(i_global<=homo) CYCLE
        i_global=i_global-homo
        IF(.NOT.(my_B_virtual_start<=i_global.AND.i_global<=my_B_virtual_end)) CYCLE
        DO jjB=1, send_col_size
          j_global=cp_fm_indxl2g(jjB,ncol_block,send_pcol,&
                                 mp2_env%ri_grad%P_mo%matrix_struct%first_p_pos(2),npcol)
          IF(j_global<=homo) CYCLE
          j_global=j_global-homo
          ab_send(iiB,jjB)=mp2_env%ri_grad%P_ab(i_global-my_B_virtual_start+1,j_global)
        END DO
      END DO

      ab_rec=0.0_dp
      CALL  mp_sendrecv(ab_send,proc_send,&
                        ab_rec,proc_receive,&
                        para_env_sub%group)

      mp2_env%ri_grad%P_mo%local_data(1:nrow_local,1:ncol_local)=&
      mp2_env%ri_grad%P_mo%local_data(1:nrow_local,1:ncol_local)+ab_rec(1:nrow_local,1:ncol_local)

      DEALLOCATE(ab_send)
    END DO
    DEALLOCATE(ab_rec)
    DEALLOCATE(mepos_2_grid)
    DEALLOCATE(sizes)

    ! deallocate the local P_ab
    DEALLOCATE(mp2_env%ri_grad%P_ab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)

    ! create also W_MP2_MO
    CALL timeset(routineN//"_WMO",handle2)
    NULLIFY(mp2_env%ri_grad%W_mo)
    CALL cp_fm_create(mp2_env%ri_grad%W_mo, fm_struct_tmp, name="W_MP2_MO",error=error)
    CALL cp_fm_set_all(mp2_env%ri_grad%W_mo, 0.0_dp,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! all block
    CALL cp_gemm('T','N',dimen,dimen,dimen,2.0_dp,&
                    L_mu_q,mo_coeff,0.0_dp,mp2_env%ri_grad%W_mo,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)

    ! occ-occ block
    CALL cp_gemm('T','N',homo,homo,dimen,-2.0_dp,&
                    L_mu_q,mo_coeff,0.0_dp,mp2_env%ri_grad%W_mo,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)

    ! occ-virt block
    CALL cp_gemm('T','N',homo,virtual,dimen,2.0_dp,&
                    mo_coeff,L_mu_q,0.0_dp,mp2_env%ri_grad%W_mo,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=homo+1,&
                    b_first_row=1,&
                    c_first_col=homo+1,&
                    c_first_row=1)
    CALL timestop(handle2)

    ! IF(para_env%mepos==0) CALL write_array(mp2_env%ri_grad%P_ij(1:homo,1:homo))

    ! Calculate occ-virt block of the lagrangian in MO
    CALL timeset(routineN//"_Ljb",handle2)
    NULLIFY(mp2_env%ri_grad%L_jb, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=homo,ncol_global=virtual,error=error)
    CALL cp_fm_create(mp2_env%ri_grad%L_jb, fm_struct_tmp, name="fm_L_jb",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(mp2_env%ri_grad%L_jb, 0.0_dp,error=error)
 
    ! first Virtual
    CALL cp_gemm('T','N',homo,virtual,dimen,1.0_dp,&
                    L_mu_q,mo_coeff,0.0_dp,mp2_env%ri_grad%L_jb,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=homo+1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    ! then occupied
    CALL cp_gemm('T','N',homo,virtual,dimen,2.0_dp,&
                    mo_coeff,L_mu_q,2.0_dp,mp2_env%ri_grad%L_jb,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=homo+1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    ! finally release L_mu_q
    CALL cp_fm_release(L_mu_q,error=error)
    CALL timestop(handle2)

    ! here we should be done next CPHF

    IF(.FALSE.) THEN   
      ALLOCATE(ab_send(homo,virtual))
      ab_send=0.0_dp
      CALL cp_fm_get_info(matrix=mp2_env%ri_grad%L_jb,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          ab_send(i_global,j_global)=mp2_env%ri_grad%L_jb%local_data(iiB,jjB)
        END DO
      END DO
      CALL mp_sum(ab_send,para_env%group)
      IF(para_env%mepos==0) THEN
        WRITE(*,*) 'Ljb'
        CALL write_array(ab_send(1:homo,1:virtual))
      END IF
      DEALLOCATE(ab_send)
      ALLOCATE(ab_send(homo,homo))
      ALLOCATE(ab_rec(virtual,virtual))
      ALLOCATE(mat_deb(homo,virtual))
      ab_send=0.0_dp
      ab_rec=0.0_dp
      mat_deb=0.0_dp
      CALL cp_fm_get_info(matrix=mp2_env%ri_grad%P_mo,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          IF(i_global<=homo.AND.j_global<=homo) THEN
            ab_send(i_global,j_global)=mp2_env%ri_grad%P_mo%local_data(iiB,jjB)
          ELSE
            IF(i_global>homo.AND.j_global>homo) THEN
              ab_rec(i_global-homo,j_global-homo)=mp2_env%ri_grad%P_mo%local_data(iiB,jjB)
            ELSE
              mat_deb(MIN(i_global,j_global),MAX(i_global,j_global)-homo)=&
              mat_deb(MIN(i_global,j_global),MAX(i_global,j_global)-homo)+mp2_env%ri_grad%P_mo%local_data(iiB,jjB)
            END IF
          END IF
        END DO
      END DO
      CALL mp_sum(ab_send,para_env%group)
      CALL mp_sum(ab_rec,para_env%group)
      CALL mp_sum(mat_deb,para_env%group)
      IF(para_env%mepos==0) THEN
        WRITE(*,*) 'Pij'
        CALL write_array(ab_send(1:homo,1:homo))
        WRITE(*,*) 'Pab'
        CALL write_array(ab_rec(1:virtual,1:virtual))
        WRITE(*,*) 'Pia' 
        CALL write_array(mat_deb(1:homo,1:virtual))
      END IF
      DEALLOCATE(ab_send, ab_rec, mat_deb)
      ALLOCATE(ab_send(dimen,dimen))
      ab_send=0.0_dp
      CALL cp_fm_get_info(matrix=mp2_env%ri_grad%W_mo,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          ab_send(i_global,j_global)=mp2_env%ri_grad%W_mo%local_data(iiB,jjB)
        END DO
      END DO
      CALL mp_sum(ab_send,para_env%group)
      IF(para_env%mepos==0) THEN
        WRITE(*,*) 'Wij'
        CALL write_array(ab_send(1:dimen,1:dimen))
      END IF
      DEALLOCATE(ab_send)
    END IF

    DEALLOCATE(sub_proc_map)
    DEALLOCATE(pos_info)

    CALL timestop(handle)

  END SUBROUTINE create_W_P


! *****************************************************************************
!> \brief ...
!> \param mat ...
!> \param unitout ...
! *****************************************************************************
   SUBROUTINE write_array(mat,unitout)
    REAL(KIND=dp), DIMENSION(:, :)           :: mat
    INTEGER, OPTIONAL                        :: unitout

    INTEGER                                  :: iii, jjj

     WRITE(*,*)
     DO iii=1, SIZE(mat,1)
       WRITE(*,*) iii
       DO jjj=1, SIZE(mat,2), 20
         IF(PRESENT(unitout)) THEN
           WRITE(1000+unitout,'(1000F10.5)') mat(iii,jjj:MIN(SIZE(mat,2),jjj+19))
         ELSE
           WRITE(*,'(1000F10.5)') mat(iii,jjj:MIN(SIZE(mat,2),jjj+19))
         END IF
       END DO
       WRITE(*,*)
     END DO
     WRITE(*,*)
  END SUBROUTINE

END MODULE mp2_ri_grad
