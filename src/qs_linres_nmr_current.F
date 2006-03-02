!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_nmr_current *
!!
!!   NAME 
!!     qs_linres_nmr_current
!!
!!   FUNCTION
!!     given the response wavefunctions obtained by the application
!!     of the (rxp), p, and ((dk-dl)xp) operators,
!!     here the current density vector (jx, jy, jz)
!!     is computed for the 3 directions of the magnetic field (Bx, By, Bz)   
!!     
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     created 02-2006 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_current


  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_scale,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: l_info_type,&
                                             return_l_info
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_collocate_density,            ONLY: calculate_total_rho,&
                                             density_on_full_grid,&
                                             lgrid_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_linres_nmr_op,                ONLY: fac_vecp,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_nmr_shift,             ONLY: interpolate_shift_pwgrid
  USE qs_linres_types,                 ONLY: nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_get_loop_vars,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_to_cube
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_type,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE


  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: nmr_response_current

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_current'

!!***
! *****************************************************************************

CONTAINS

!!****f*  qs_linres_nmr_current%nmr_response_current
!!
!!  NAME 
!!      nmr_response_current
!!
!!  FUNCTION
!!      First calculate the density matrixes, for each component of the current
!!      they are 3 because of the r dependent terms
!!      Next it collocates on the grid to have J(r)
!!      In the GAPW case one need to collocate on the PW grid only the soft part
!!      while the rest goes on Lebedev grids
!!      The contributions to the shift and to the susceptibility will be
!!      calulated separatedly and added only at the end
!!      The calculation of the shift tensor is performed on the position of the atoms
!!      and on other specified points in real space summing up the contributions
!!      from the PW grid current density and the local densities
!!     
!!  ARGUMENTS
!!      nmr_env
!!      qs_env
!!      psi1, p_psi1 : scratch MOS coefficients 
!!      error
!!
!!  NOTES
!!      To calculate the susceptibility on the PW gris it is necessary to apply
!!      the position operator yet another time.
!!      This cannot be done on directly on J(r) because it is not localized and 
!!      the usual trick of the wfn center cannot be used.
!!      The idea would be to construct additional density matrices that somehow
!!      contain the information on the wannier centers. 
!!      This would mean that the entire produre used to compute on the pw grid 
!!      the terms of the current  density of the type (r'-dk) x psi1_p(r'),
!!      should be repeated to compute the function r x J(r)
!!
!!  AUTHOR
!!    MI
!!
!!*** **********************************************************************

  SUBROUTINE nmr_response_current(nmr_env,qs_env,psi1,p_psi1,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi1, p_psi1
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_response_current', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: ext, filename
    INTEGER :: handle, homo, i_B, iao, iatom, idir, idir2, idir3, ii_B, &
      iii_B, ispin, istat, istate, natom, nmo, nspins, unit_nr
    LOGICAL                                  :: failure, gapw, ionode, uni_occ
    REAL(dp)                                 :: alpha, dk(3), jrho_tot, &
                                                maxocc, rmu(3), rmu_dk(3), &
                                                scale_fac
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_psi1
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: pw_gspace_work, &
                                                shift_pw_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: shift_pw_gspace
    TYPE(coeff_type), POINTER                :: rho_gspace, rho_rspace, wf_r
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: vecbuf
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(real_matrix_type), POINTER          :: density_matrix, &
                                                density_matrix_2_ii, &
                                                density_matrix_2_iii
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(section_vals_type), POINTER         :: nmr_section

!   ---------------------------------------------------------------------------

    failure = .FALSE.
    gapw=.FALSE.
    NULLIFY(rho_rspace,rho_gspace,logger,nmr_section)
    NULLIFY(density_matrix,density_matrix_2_ii,density_matrix_2_iii)
    NULLIFY(cell,dft_control,mos,rho_struct,particle_set,pw_env)
    NULLIFY(auxbas_rs_pool,auxbas_pw_pool,pw_pools,wf_r)
    NULLIFY(shift_pw_gspace)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
!    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

    nspins = SIZE(psi1,1)
   
   
    CALL get_qs_env(qs_env=qs_env,rho=rho_struct,&
         cell=cell, dft_control=dft_control,mos=mos,&
         particle_set=particle_set)
    gapw = dft_control%qs_control%gapw
    natom = SIZE(particle_set,1)
    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    nmr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES%NMR",&
                   error=error)

    IF (.NOT. failure) THEN
      CALL timeset("nmr_response_current","I"," ",handle)

      IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
                  "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN

          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_rs_pool=auxbas_rs_pool,&
               auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools)
          CALL rs_pool_create_rs(auxbas_rs_pool,rs,force_env_section=qs_env%input,&
               error=error)
          
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
                  use_data = REALDATA3D,&
                  in_space = REALSPACE, error=error)
      END IF

!-----------------------------------------------------------------------!
!     Allocate grids for the calculation of jrho and the shift 
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,pw_gspace_work,&
               error=error)
      ALLOCATE(shift_pw_gspace(3),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO idir = 1,3
        CALL pw_pool_init_coeff(auxbas_pw_pool,shift_pw_gspace(idir), &
                                use_data = COMPLEXDATA1D,&
                                in_space = RECIPROCALSPACE, error=error)
      END DO
      CALL pw_pool_init_coeff(auxbas_pw_pool,shift_pw_rspace,&
                 use_data=REALDATA3D, in_space=REALSPACE)
!
!-----------------------------------------------------------------------!


!-----------------------------------------------------------------------!
!     Initialize shift and chi
      nmr_env%chemical_shift = 0.0_dp
      nmr_env%chemical_shift_loc = 0.0_dp
      IF(nmr_env%do_nics) THEN
        nmr_env%chemical_shift_nics = 0.0_dp
        nmr_env%chemical_shift_loc_nics = 0.0_dp
      END IF
      nmr_env%chi_tensor = 0.0_dp
      nmr_env%chi_tensor_loc = 0.0_dp
!-----------------------------------------------------------------------!

! Loop on the field direction 
      DO i_B = 1,3
       DO ispin = 1,dft_control%nspins
         CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
         CALL cp_fm_set_all(p_psi1(ispin)%matrix,0.0_dp,error=error)
       END DO

       CALL set_vecp(i_B,ii_B,iii_B)
       DO ispin = 1,nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
              occupation_numbers=occupation, homo=homo, nmo=nmo,&
              uniform_occupation=uni_occ,maxocc=maxocc)

         !First contribution from the response psi1_p
         DO istate = 1,nmr_env%nstates(ispin)
           CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,ii_B)%matrix,&
                vecbuf(ii_B)%array,1,istate,nmr_env%nao,1,transpose=.TRUE.,&
                error=error)
           CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,iii_B)%matrix,&
                vecbuf(iii_B)%array,1,istate,nmr_env%nao,1,transpose=.TRUE.,&
                error=error)
           dk(1:3) = nmr_env%centers_set(ispin)%array(1:3,istate)
           DO iao = 1,nmr_env%nao
             Rmu(1:3) =  nmr_env%basisfun_center(1:3,iao)
             Rmu_dk = pbc(dk,Rmu,cell)
             vecbuf_psi1(1,iao) = Rmu_dk(ii_B)*vecbuf(iii_B)%array(1,iao)-&
                                  Rmu_dk(iii_B)*vecbuf(ii_B)%array(1,iao)
           END DO  ! iao
           ! Copy the vector in the full matrix psi1
           CALL cp_fm_set_submatrix(psi1(ispin)%matrix,vecbuf_psi1,&
                1,istate,nmr_env%nao,1,transpose=.TRUE.,error=error)
         END DO  ! istate

         !Add to psi the psi1_rxp contribution
         CALL  cp_fm_scale_and_add(-1.0_dp,psi1(ispin)%matrix,1.0_dp,nmr_env%psi1_rxp(ispin,i_B)%matrix,error=error)
         IF(nmr_env%full_nmr) THEN
            CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,-1.0_dp,nmr_env%psi1_D(ispin,i_B)%matrix,error=error)
         END IF

         DO idir = 1,3
         ! Apply the p operator
           CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
                p_psi1(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
           density_matrix => nmr_env%jp1_ao
           CALL set_matrix(density_matrix,0.0_dp)

         ! Contribution to the current density matrix that does not depend 
         !   on the position r where J(r) is computed
           IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
             CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
           ELSE
             alpha = maxocc
           END IF
           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=nmr_env%p_psi0(ispin,idir)%matrix,&
                            matrix_g=psi1(ispin)%matrix,ncol=homo,&
                            alpha=alpha)

           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=mo_coeff,&
                            matrix_g=p_psi1(ispin)%matrix,ncol=homo,&
                            alpha=-alpha)
          
           ! For the position dependent term two more density matrix are needed
           density_matrix_2_ii => nmr_env%jp2_ao(1)%matrix
           CALL set_matrix(density_matrix_2_ii,0.0_dp)
           density_matrix_2_iii => nmr_env%jp2_ao(2)%matrix
           CALL set_matrix(density_matrix_2_iii,0.0_dp)

           CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,iii_B)%matrix,psi1(ispin)%matrix)
           CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
                p_psi1(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
           IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
             CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
           ELSE
             alpha = maxocc
           END IF
           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_ii,&
                            matrix_v=nmr_env%p_psi0(ispin,idir)%matrix,&
                            matrix_g=psi1(ispin)%matrix,ncol=homo,&
                            alpha=alpha)

           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_ii,&
                            matrix_v=mo_coeff,&
                            matrix_g=p_psi1(ispin)%matrix,ncol=homo,&
                            alpha=-alpha)

           CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,ii_B)%matrix,psi1(ispin)%matrix)
           CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
                p_psi1(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
           IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
             CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
           ELSE
             alpha = maxocc
           END IF
           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_iii,&
                            matrix_v=nmr_env%p_psi0(ispin,idir)%matrix,&
                            matrix_g=psi1(ispin)%matrix,ncol=homo,&
                            alpha=alpha)

           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_iii,&
                            matrix_v=mo_coeff,&
                            matrix_g=p_psi1(ispin)%matrix,ncol=homo,&
                            alpha=-alpha)

        ! Calculate the current density on the pw grid (only soft if GAPW)
        ! This is vector component idir of the response current density 
        ! generated by the magnetic field in cartesian direction i_B
        ! Use the qs_rho_type already  used for rho during the scf

            IF(nmr_env%store_current) THEN
              rho_rspace => nmr_env%jrho1_set(idir,i_B)%rho%rho_r(ispin)
              rho_gspace => nmr_env%jrho1_set(idir,i_B)%rho%rho_g(ispin)
            ELSE
              rho_rspace => rho_struct%rho_r(ispin)
              rho_gspace => rho_struct%rho_g(ispin)
            END IF
            CALL calculate_jrho_resp(density_matrix, density_matrix_2_ii,&
                 density_matrix_2_iii, ii_B, iii_B, rho_rspace, rho_gspace,qs_env,&
                 gapw, error=error)

            jrho_tot = calculate_total_rho(rho_gspace)

            ! loop over the Gvec  components: x,y,z
            DO idir2 = 1,3
              IF(idir /= idir2) THEN
                ! in reciprocal space multiply (G_idir2(i)/G(i)^2)J_(idir)(G(i))
!                CALL mult_G_ov_G2_grid(cell,rho_gspace,pw_gspace_work,idir2,error=error)

                ! scale and add to the correct component of the shift column
                 CALL set_vecp_rev(idir,idir2,idir3)
                 scale_fac=fac_vecp(idir3,idir2,idir)
                 CALL coeff_scale(pw_gspace_work,scale_fac)
                 CALL coeff_sumup(pw_gspace_work,shift_pw_gspace(idir3))
              END IF
            END DO

            IF(gapw) THEN
            ! compute the atomic response current densities on the spherical grids
            END IF

         END DO  ! idir
         
       END DO  !  ispin

       ! Tranform the column i_B of the shift tensor from reciprocal to real space
       ! get the values on the atomic positions and on the other required points (if any) 
       ! this can be done by interpolation of the values of the grid 
       ! on the required positions in real space

       DO idir = 1,3
         CALL coeff_transform_space(shift_pw_gspace(idir),shift_pw_rspace)
         CALL interpolate_shift_pwgrid(nmr_env,particle_set,cell,shift_pw_rspace,&
              i_B,idir,error=error)
       END DO 
       IF(gapw) THEN
         ! Add the local contributions by numerical integration over the spherical grids
       END IF

       IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
                  "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN

          DO idir = 1,3
            CALL coeff_zero(wf_r)
            DO ispin =1 ,nspins
              CALL coeff_copy(nmr_env%jrho1_set(idir,i_B)%rho%rho_r(ispin),wf_r)
            END DO 
            IF(gapw) THEN
            ! Add the local hard and soft contributions
            ! This can be done atom by atom by an spline extrapolation of the  values
            ! on the spherical grid, on the classical grid points. It is not very
            ! accurate, but it is used only for the printing of cube files.
              DO iatom = 1,natom
              END DO
            END IF
            filename="jresp"
            WRITE(ext,'(a2,I1,a2,I1,a5)')  "iB",i_B,"_d",idir,".cube"
            unit_nr=cp_print_key_unit_nr(logger,nmr_section,"PRINT%CURRENT_CUBES",&
                    extension=TRIM(ext),middle_name=TRIM(filename),&
                    log_filename=.FALSE.,error=error)
            CALL rs_pw_to_cube(wf_r%pw,unit_nr,ionode,"RESPONSE CURRENT DENSITY ",&
                    rs,stride=section_get_ival(nmr_section,"PRINT%CURRENT_CUBE%STRIDE"),&
                    force_env_section=qs_env%input, error=error)
            CALL cp_print_key_finished_output(unit_nr,logger,nmr_section,&
                  "PRINT%CURRENT_CUBES",error=error)
          END DO 

       END IF
      END DO  ! i_B

      IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
                  "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r,&
               error=error)
        CALL rs_pool_give_back_rs(auxbas_rs_pool,rs, error=error)
      END IF

!-----------------------------------------------------------------------!
!     Dellocate grids for the calculation of jrho and the shift 
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,pw_gspace_work,&
               error=error)
      DO idir = 1,3
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,shift_pw_gspace(idir), &
                                 error=error)
      END DO
      DEALLOCATE(shift_pw_gspace,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,shift_pw_rspace,&
               error=error)
!
!-----------------------------------------------------------------------!

      CALL timestop(0.0_dp,handle)
    END IF  !  failure

  END SUBROUTINE nmr_response_current

 ! *****************************************************************************
 
  SUBROUTINE calculate_jrho_resp(mat_jp,mat_jp_rii,mat_jp_riii,ii_B,iii_B,rho_rs, rho_gs, qs_env, soft_valid, error)

    TYPE(real_matrix_type), POINTER          :: mat_jp, mat_jp_rii, &
                                                mat_jp_riii
    INTEGER, INTENT(IN)                      :: ii_B, iii_B
    TYPE(coeff_type), INTENT(INOUT)          :: rho_rs, rho_gs
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_resp', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, bcol, brow, curr_tasks, dir, handle, i, iatom, &
      igrid_level, ijatoms, ijsets, ikind, ilist, inode, ipgf, iset, istat, &
      itask, ithread, j, jatom, jkind, jpgf, jset, k, kind_i, maxco, maxsgf, &
      maxsgf_set, n, na1, na2, natom_pairs, nb1, nb2, ncoi, ncoj, nkind, &
      nlist, nnode, npme, nseta, nsetb, nthread, omp_get_max_threads, &
      omp_get_thread_num, sgfb, sgfi, sgfj, stat, tp
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb, lb_max, lb_min, &
      li_max, li_min, lj_max, lj_min, npgfa, npgfb, nsgfa, nsgfb, nsgfi, &
      nsgfj, ntasks, ub
    INTEGER, DIMENSION(:, :), POINTER        :: asets, atasks, first_sgfa, &
                                                first_sgfb, ival, latom, &
                                                tasks_local
    INTEGER, DIMENSION(:, :, :), POINTER     :: tasks
    LOGICAL                                  :: distributed_rs_grids, &
                                                failure, map_consistent, &
                                                my_soft
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, &
                                                kind_radius_a, kind_radius_b, &
                                                rab2, scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb, rp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: dab_local, jp_block, &
      jp_block_rii, jp_block_riii, jpab, jpab_ii, jpab_iii, jpblock, &
      jpblock_rii, jpblock_riii, rpgfa, rpgfb, sphi_a, sphi_b, sphi_i, &
      sphi_j, work, zeta, zetb, zeti, zetj
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dist_ab, jpabt, jpabt_ii, &
                                                jpabt_iii, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind, atomic_kind_i
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(l_info_type)                        :: l_info
    TYPE(lgrid_type)                         :: lgrid
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(real_matrix_type), POINTER          :: deltajp, deltajp_rii, &
                                                deltajp_riii
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(section_vals_type), POINTER         :: input, interp_section

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kind,cell,dft_control,orb_basis_set,sab_orb_neighbor_list,&
         sab_orb_neighbor_node,atomic_kind_set,sab_orb,particle_set,&
         rs_rho,pw_env,rs_pools,para_env,dist_ab,dab_local,&
         set_radius_a,set_radius_b,la_max,la_min,&
         lb_max,lb_min,npgfa,npgfb,nsgfa,nsgfb,&
         rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb,first_sgfa,first_sgfb,&
         tasks,tasks_local,ival,latom,ntasks,asets,atasks,workt)
    NULLIFY(li_max,li_min, lj_max, lj_min,nsgfi,nsgfj,zeti,zetj)
    NULLIFY(deltajp,deltajp_rii,deltajp_riii)
    NULLIFY(jp_block,jp_block_rii,jp_block_riii)
    NULLIFY(jpblock,jpblock_rii,jpblock_riii)
    NULLIFY(jpabt,jpabt_ii,jpabt_iii)
    NULLIFY(lgrid%r,atomic_kind_i)

!    debug_count=debug_count+1

    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    CALL timeset("calculate__jrho_resp","I"," ",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    para_env=para_env,&
                    input=input,&
                    pw_env=pw_env)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info
    l_info=pw_env%l_info

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
!    CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)

    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, error=error)


    ! *** set up the rs multi-grids
    distributed_rs_grids=.FALSE.
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, force_env_section=input,&
         error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       IF ( rs_rho(igrid_level)%rs_grid%direction /= 0 ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf,&
                             maxsgf_set=maxsgf_set)

    IF ( nthread > 1 ) THEN
      n=0
      DO igrid_level = 1,gridlevel_info%ngrid_levels
        n = MAX(n,rs_rho(igrid_level)%rs_grid%ngpts_local)
      END DO
      n = n*nthread
      CALL reallocate(lgrid%r,1,n)
    END IF

    nkind = SIZE(atomic_kind_set)

    CALL reallocate(jpabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_ii,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_iii,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(ntasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks,1,8,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(dist_ab,1,3,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks_local,1,2,1,max_tasks)
    CALL reallocate(ival,1,6,1,max_tasks)
    CALL reallocate(latom,1,2,1,max_tasks)
    CALL reallocate(dab_local,1,3,1,max_tasks)
    CALL reallocate(atasks,1,2,1,max_tasks)
    CALL reallocate(asets,1,2,1,max_tasks)
    curr_tasks = max_tasks

!   *** Initialize working density matrix ***

    ! distributed rs grids require a matrix that will be changed (rs_get_my_tasks)
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
        CALL allocate_matrix(matrix=deltajp,&
                         nrow=mat_jp%nrow,&
                         ncol=mat_jp%ncol,&
                         nblock_row=mat_jp%nblock_row,&
                         nblock_col=mat_jp%nblock_col,&
                         first_row=mat_jp%first_row(:),&
                         last_row=mat_jp%last_row(:),&
                         first_col=mat_jp%first_col(:),&
                         last_col=mat_jp%last_col(:),&
                         matrix_name="DeltaP",&
                         sparsity_id=-1, &   ! basically unknown sparsity in parallel
                         matrix_symmetry=mat_jp%symmetry)
        CALL allocate_matrix(matrix=deltajp_rii,&
                         nrow=mat_jp%nrow,&
                         ncol=mat_jp%ncol,&
                         nblock_row=mat_jp%nblock_row,&
                         nblock_col=mat_jp%nblock_col,&
                         first_row=mat_jp%first_row(:),&
                         last_row=mat_jp%last_row(:),&
                         first_col=mat_jp%first_col(:),&
                         last_col=mat_jp%last_col(:),&
                         matrix_name="DeltaP",&
                         sparsity_id=-1, &   ! basically unknown sparsity in parallel
                         matrix_symmetry=mat_jp%symmetry)
        CALL allocate_matrix(matrix=deltajp_riii,&
                         nrow=mat_jp%nrow,&
                         ncol=mat_jp%ncol,&
                         nblock_row=mat_jp%nblock_row,&
                         nblock_col=mat_jp%nblock_col,&
                         first_row=mat_jp%first_row(:),&
                         last_row=mat_jp%last_row(:),&
                         first_col=mat_jp%first_col(:),&
                         last_col=mat_jp%last_col(:),&
                         matrix_name="DeltaP",&
                         sparsity_id=-1, &   ! basically unknown sparsity in parallel
                         matrix_symmetry=mat_jp%symmetry)
    ELSE
        deltajp=>mat_jp
        deltajp_rii=>mat_jp
        deltajp_riii=>mat_jp
    ENDIF


    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           softb = my_soft, &
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             kind_radius=kind_radius_a,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             set_radius=set_radius_a,&
                             sphi=sphi_a,&
                             zet=zeta)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           softb = my_soft, &
                           orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgfb,&
                               kind_radius=kind_radius_b,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        ab = ikind + nkind*(jkind - 1)

        IF (ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN

           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_orb(ab)%neighbor_list_set,&
                                      nlist=nlist)
           sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
        ELSE
           nlist=0
        END IF

        ntasks = 0
        tasks = 0

       DO ilist = 1, nlist

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom,nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO inode = 1, nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

            brow = iatom
            bcol = jatom

             CALL get_block_node(matrix=mat_jp,&
                                 block_row=brow,&
                                 block_col=bcol,&
                                 BLOCK=jp_block)
             CALL get_block_node(matrix=mat_jp_rii,&
                                 block_row=brow,&
                                 block_col=bcol,&
                                 BLOCK=jp_block_rii)
             CALL get_block_node(matrix=mat_jp_riii,&
                                 block_row=brow,&
                                 block_col=bcol,&
                                 BLOCK=jp_block_riii)

             IF (.NOT.ASSOCIATED(jp_block)) THEN
               sab_orb_neighbor_node => next(sab_orb_neighbor_node)
               CYCLE
             END IF

             IF (distributed_rs_grids) THEN
                 NULLIFY (jpblock,jpblock_rii,jp_block_riii )
                 CALL add_block_node ( deltajp, brow, bcol, jpblock )
                 jpblock = jp_block
                 CALL add_block_node ( deltajp_rii, brow, bcol, jpblock_rii )
                 jpblock_rii = jp_block_rii
                 CALL add_block_node ( deltajp_riii, brow, bcol, jpblock_riii )
                 jpblock_riii = jp_block_riii
             ELSE
                 jpblock => jp_block
                 jpblock_rii => jp_block_rii
                 jpblock_riii => jp_block_riii
             ENDIF

             IF (.NOT. map_consistent) THEN
                IF ( ALL ( 100.0_dp*ABS(jpblock) < eps_rho_rspace) ) THEN
                  sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                  CYCLE
                END IF
             END IF

             rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
             dab = SQRT(rab2)

             DO iset=1,nseta
               IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE

               DO jset = 1,nsetb
                 IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                 DO ipgf=1,npgfa(iset)
                   IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

                   DO jpgf=1,npgfb(jset)
                     IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                     zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

                     IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                         igrid_level = 1
                     ELSE
                         igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                     ENDIF

                     ntasks(igrid_level) = ntasks(igrid_level) + 1
                     n = ntasks(igrid_level)
                     IF ( n > curr_tasks ) THEN
                       curr_tasks = curr_tasks*mult_tasks
                       CALL reallocate(tasks,1,8,1,curr_tasks,&
                                       1,gridlevel_info%ngrid_levels)
                       CALL reallocate(dist_ab,1,3,1,curr_tasks,&
                                       1,gridlevel_info%ngrid_levels)
                     END IF
                     dir = rs_rho(igrid_level)%rs_grid%direction
                     tasks (1,n,igrid_level) = n
                     IF ( dir /= 0) THEN
                       rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                       rp(:) = pbc(rp,cell)
                       tp = FLOOR(rp(dir)/rs_rho(igrid_level)%rs_grid%dr(dir))
                       tp = MODULO ( tp, rs_rho(igrid_level)%rs_grid%npts(dir) )
                       tasks (2,n,igrid_level) = tp + rs_rho(igrid_level)%rs_grid%lb(dir)
                     END IF
                     tasks (3,n,igrid_level) = iatom
                     tasks (4,n,igrid_level) = jatom
                     tasks (5,n,igrid_level) = iset
                     tasks (6,n,igrid_level) = jset
                     tasks (7,n,igrid_level) = ipgf
                     tasks (8,n,igrid_level) = jpgf
                     dist_ab (1,n,igrid_level) = rab(1)
                     dist_ab (2,n,igrid_level) = rab(2)
                     dist_ab (3,n,igrid_level) = rab(3)

                   END DO  ! jpgf

                 END DO  ! ipgf

               END DO  ! jset

             END DO  ! iset
!!!!!!!!!!
!!!!!!!!!!  The response current density is not a symmetric function
!!!!!!!!!!
            IF(iatom/=jatom) THEN
              brow = jatom
              bcol = iatom
  
               CALL get_block_node(matrix=mat_jp,&
                                   block_row=brow,&
                                   block_col=bcol,&
                                   BLOCK=jp_block)
               CALL get_block_node(matrix=mat_jp_rii,&
                                   block_row=brow,&
                                   block_col=bcol,&
                                   BLOCK=jp_block_rii)
               CALL get_block_node(matrix=mat_jp_riii,&
                                   block_row=brow,&
                                   block_col=bcol,&
                                   BLOCK=jp_block_riii)

               IF (.NOT.ASSOCIATED(jp_block)) THEN
                 sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                 CYCLE
               END IF

               IF (distributed_rs_grids) THEN
                 NULLIFY (jpblock,jpblock_rii,jp_block_riii )
                 CALL add_block_node ( deltajp, brow, bcol, jpblock )
                 jpblock = jp_block
                 CALL add_block_node ( deltajp_rii, brow, bcol, jpblock_rii )
                 jpblock_rii = jp_block_rii
                 CALL add_block_node ( deltajp_riii, brow, bcol, jpblock_riii )
                 jpblock_riii = jp_block_riii
               ELSE
                 jpblock => jp_block
                 jpblock_rii => jp_block_rii
                 jpblock_riii => jp_block_riii
               ENDIF

               IF (.NOT. map_consistent) THEN
                 IF ( ALL ( 100.0_dp*ABS(jpblock) < eps_rho_rspace) ) THEN
                   sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                   CYCLE
                 END IF
               END IF

               rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
               dab = SQRT(rab2)

               DO jset=1,nsetb
                 IF (set_radius_b(jset) + kind_radius_a < dab) CYCLE
 
                 DO iset = 1,nseta
                   IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                   DO jpgf=1,npgfb(jset)
                     IF (rpgfb(jpgf,jset) + set_radius_a(iset) < dab) CYCLE

                     DO ipgf=1,npgfa(iset)
                       IF (rpgfb(jpgf,jset) + rpgfa(ipgf,iset) < dab) CYCLE

                       zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

                       IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                           igrid_level = 1
                       ELSE
                           igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                       ENDIF

                       ntasks(igrid_level) = ntasks(igrid_level) + 1
                       n = ntasks(igrid_level)
                       IF ( n > curr_tasks ) THEN
                         curr_tasks = curr_tasks*mult_tasks
                         CALL reallocate(tasks,1,8,1,curr_tasks,&
                                         1,gridlevel_info%ngrid_levels)
                         CALL reallocate(dist_ab,1,3,1,curr_tasks,&
                                         1,gridlevel_info%ngrid_levels)
                       END IF
                       dir = rs_rho(igrid_level)%rs_grid%direction
                       tasks (1,n,igrid_level) = n
                       IF ( dir /= 0) THEN
                         rb(:) = ra(:) + rab(:)
                         rp(:) = rb(:) + zeta(ipgf,iset)/zetp*rab(:)
                         rp(:) = pbc(rp,cell)
                         tp = FLOOR(rp(dir)/rs_rho(igrid_level)%rs_grid%dr(dir))
                         tp = MODULO ( tp, rs_rho(igrid_level)%rs_grid%npts(dir) )
                         tasks (2,n,igrid_level) = tp + rs_rho(igrid_level)%rs_grid%lb(dir)
                       END IF
                       tasks (3,n,igrid_level) = jatom
                       tasks (4,n,igrid_level) = iatom
                       tasks (5,n,igrid_level) = jset
                       tasks (6,n,igrid_level) = iset
                       tasks (7,n,igrid_level) = jpgf
                       tasks (8,n,igrid_level) = ipgf
                       dist_ab (1,n,igrid_level) = -rab(1)
                       dist_ab (2,n,igrid_level) = -rab(2)
                       dist_ab (3,n,igrid_level) = -rab(3)

                     END DO  ! jpgf

                   END DO  ! ipgf

                 END DO  ! jset

               END DO  ! iset
             END IF  ! iatom/=jatom
!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
             sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

        DO igrid_level = 1, gridlevel_info%ngrid_levels
          n = ntasks ( igrid_level )
          IF ( n > SIZE ( tasks_local, 2 ) ) &
            CALL reallocate(tasks_local,1,2,1,n)
          IF ( n > SIZE ( ival, 2 ) ) &
            CALL reallocate(ival,1,6,1,n)
          IF ( n > SIZE ( latom, 2 ) ) &
            CALL reallocate(latom,1,2,1,n)
          IF ( n > SIZE ( dab_local, 2 ) ) &
            CALL reallocate(dab_local,1,3,1,n)

!$OMP parallel do private(i)

          DO i=1,n
            tasks_local(1,i) = tasks(1,i,igrid_level)
            tasks_local(2,i) = tasks(2,i,igrid_level)
            latom(1,i) = tasks(3,i,igrid_level)
            latom(2,i) = tasks(4,i,igrid_level)
            ival(1,i) = tasks(3,i,igrid_level)
            ival(2,i) = tasks(4,i,igrid_level)
            ival(3,i) = tasks(5,i,igrid_level)
            ival(4,i) = tasks(6,i,igrid_level)
            ival(5,i) = tasks(7,i,igrid_level)
            ival(6,i) = tasks(8,i,igrid_level)
            dab_local(1,i) = dist_ab(1,i,igrid_level)
            dab_local(2,i) = dist_ab(2,i,igrid_level)
            dab_local(3,i) = dist_ab(3,i,igrid_level)
          END DO
!$OMP parallel do private(i)
          DO i=n+1,SIZE(tasks_local,2)
            tasks_local(1,i) = 0
            tasks_local(2,i) = 0
          END DO
          npme = 0
          ! modifies deltap if distributed_rs_grids
          CALL rs_get_my_tasks ( rs_rho(igrid_level)%rs_grid, tasks_local, &
               npme, ival=ival, rval=dab_local, pmat=deltajp, pmat2=deltajp_rii,&
               pmat3=deltajp_riii, pcor=latom, symmetric=.FALSE. )
          CALL rs_get_loop_vars ( npme, ival, natom_pairs, asets, atasks )

          IF ( nthread > 1 .AND. natom_pairs > 0 ) THEN
            lb => rs_rho(igrid_level)%rs_grid%lb_local
            ub => rs_rho(igrid_level)%rs_grid%ub_local
            lgrid%ldim = rs_rho(igrid_level)%rs_grid%ngpts_local
!$OMP parallel private(ithread,n)
!$          ithread = omp_get_thread_num()
            n = ithread*lgrid%ldim + 1
            CALL dcopy(lgrid%ldim,0._dp,0,lgrid%r(n),1)
!$OMP end parallel
          END IF
!$OMP parallel &
!$OMP default(none) &
!$OMP private(ijatoms,ithread,itask,iatom,jatom,ra,brow,bcol,p_block) &
!$OMP private(ijsets,iset,jset,ncoi,ncoj,sgfi,sgfj) &
!$OMP private(work,jpab,jpab_ii,jpab_iii,istat) &
!$OMP private(li_max,lj_max,li_min,lj_min,sphi_i,sphi_j,zeti,zetj)
!$OMP private(rb,rab,rab2,ipgf,jpgf,na1,na2,nb1,nb2,scale) &
!$OMP shared(maxco,maxsgf_set,natom_pairs,atasks,asets,ival,particle_set,cell) &
!$OMP shared(deltajp,deltajp_rii,deltajp_riii,npgfa,npgfb)&
!$OMP shared(ncoset,la_max,lb_max,first_sgfa,first_sgfb) &
!$OMP shared(nsgfa,nsgfb,sphi_a,dab_local)&
!$OMP shared(la_min,lb_min,zeta,zetb) &
!$OMP shared(rs_rho,igrid_level,cube_info,l_info,eps_rho_rspace,lgrid,nthread) &
!$OMP shared(workt,jpabt,jpabt_ii,jpabt_iii,map_consistent)
          ithread = 0
!$        ithread = omp_get_thread_num()
          jpab => jpabt(:,:,ithread)
          jpab_ii => jpabt_ii(:,:,ithread)
          jpab_iii => jpabt_iii(:,:,ithread)
          work => workt(:,:,ithread)
!$OMP do
          DO ijatoms = 1,natom_pairs
            itask = atasks(1,asets(1,ijatoms))
            iatom  = ival (1,itask)
            jatom  = ival (2,itask)
            ra(:) = pbc(particle_set(iatom)%r,cell)
!!!!!!!!All the pairs should be included at this point
            brow = iatom
            bcol = jatom
            CALL get_block_node(matrix=deltajp,&
                                block_row=brow,&
                                block_col=bcol,&
                                BLOCK=jp_block)
            CALL get_block_node(matrix=deltajp_rii,&
                                block_row=brow,&
                                block_col=bcol,&
                                BLOCK=jp_block_rii)
            CALL get_block_node(matrix=deltajp_riii,&
                                block_row=brow,&
                                block_col=bcol,&
                                BLOCK=jp_block_riii)
            IF (.NOT.ASSOCIATED(jp_block)) &
               CALL stop_program(routineP,"p_block not associated in deltap")
            DO ijsets = asets(1,ijatoms), asets(2,ijatoms)
              itask = atasks(1,ijsets)
              iset   = ival (3,itask)
              jset   = ival (4,itask)
              atomic_kind_i = particle_set(iatom)%atomic_kind
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                kind_number=kind_i)
              IF(kind_i == ikind) THEN
                ncoi = npgfa(iset)*ncoset(la_max(iset))
                sgfi = first_sgfa(1,iset)
                ncoj = npgfb(jset)*ncoset(lb_max(jset))
                sgfj = first_sgfb(1,jset) 
                sphi_i => sphi_a
                sphi_j => sphi_b
                nsgfi  => nsgfa 
                nsgfj  =>  nsgfb
                li_max => la_max
                li_min => la_min
                lj_max => lb_max
                lj_min => lb_min
                zeti   => zeta
                zetj   => zetb
              ELSEIF(kind_i == jkind) THEN
                ncoi = npgfb(iset)*ncoset(lb_max(iset))
                sgfi = first_sgfb(1,iset)
                ncoj = npgfa(jset)*ncoset(la_max(jset))
                sgfj = first_sgfa(1,jset)
                sphi_i => sphi_b
                sphi_j => sphi_a
                nsgfi  => nsgfb 
                nsgfj  =>  nsgfa
                li_max => lb_max
                li_min => lb_min
                lj_max => la_max
                lj_min => la_min
                zeti   => zetb
                zetj   => zeta
              END IF! ikind or jkind
!!
              CALL dgemm("N","N",ncoi,nsgfj(jset),nsgfi(iset),&
                         1.0_dp,sphi_i(1,sgfi),SIZE(sphi_i,1),&
                         jp_block(sgfi,sgfj),SIZE(jp_block,1),&
                         0.0_dp,work(1,1),maxco)
              CALL dgemm("N","T",ncoi,ncoj,nsgfj(jset),&
                         1.0_dp,work(1,1),maxco,&
                         sphi_j(1,sgfb),SIZE(sphi_j,1),&
                         0.0_dp,jpab(1,1),maxco)
!!
              CALL dgemm("N","N",ncoi,nsgfj(jset),nsgfi(iset),&
                         1.0_dp,sphi_i(1,sgfi),SIZE(sphi_i,1),&
                         jp_block_rii(sgfi,sgfj),SIZE(jp_block_rii,1),&
                         0.0_dp,work(1,1),maxco)
              CALL dgemm("N","T",ncoi,ncoj,nsgfj(jset),&
                         1.0_dp,work(1,1),maxco,&
                         sphi_j(1,sgfb),SIZE(sphi_j,1),&
                         0.0_dp,jpab_ii(1,1),maxco)
!!
              CALL dgemm("N","N",ncoi,nsgfj(jset),nsgfi(iset),&
                         1.0_dp,sphi_i(1,sgfi),SIZE(sphi_i,1),&
                         jp_block_riii(sgfi,sgfj),SIZE(jp_block_riii,1),&
                         0.0_dp,work(1,1),maxco)
              CALL dgemm("N","T",ncoi,ncoj,nsgfj(jset),&
                         1.0_dp,work(1,1),maxco,&
                         sphi_j(1,sgfb),SIZE(sphi_j,1),&
                         0.0_dp,jpab_iii(1,1),maxco)
!!
              DO itask = atasks(1,ijsets),atasks(2,ijsets)
                rab(:) = dab_local (:,itask)
                rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                rb(:) = ra(:) + rab(:)
                ipgf   = ival (5,itask)
                jpgf   = ival (6,itask)
                na1 = (ipgf - 1)*ncoset(li_max(iset)) + 1
                na2 = ipgf*ncoset(li_max(iset))
                nb1 = (jpgf - 1)*ncoset(lj_max(jset)) + 1
                nb2 = jpgf*ncoset(lj_max(jset))

                scale = 1.0_dp

                IF ( nthread > 1 ) THEN
                  CALL collocate_pgf_response_current(li_max(iset),zeti(ipgf,iset),&
                       li_min(iset),lj_max(jset),zetj(jpgf,jset),lj_min(jset),&
                       ra,rab,rab2,scale,ii_B,iii_B,jpab,&
                       jpab_ii,jpab_iii,na1-1,nb1-1,rs_rho(igrid_level)%rs_grid,&
                       cube_info(igrid_level),l_info,eps_rho_rspace,lgrid=lgrid,&
                       ithread=ithread,map_consistent=map_consistent)
                ELSE
                  CALL collocate_pgf_response_current(li_max(iset),zeti(ipgf,iset),&
                       li_min(iset),lj_max(jset),zetj(jpgf,jset),lj_min(jset),&
                       ra,rab,rab2,scale,ii_B,iii_B,jpab,&
                       jpab_ii,jpab_iii,na1-1,nb1-1,rs_rho(igrid_level)%rs_grid,&
                       cube_info(igrid_level),l_info,eps_rho_rspace,&
                       map_consistent=map_consistent)
                END IF  ! nthread
             END DO

            END DO

          END DO
!$OMP end parallel
          IF ( nthread > 1 .AND. natom_pairs > 0 ) THEN
            n = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)
            DO i=1,nthread
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,k) &
!$OMP shared(i,lb,ub,lgrid,rs_rho,n,igrid_level)
              DO j=lb(3),ub(3)
                k = lgrid%ldim*(i-1) + n*(j-lb(3)) + 1
                CALL daxpy (n,1._dp,lgrid%r(k),1,&
                  rs_rho(igrid_level)%rs_grid%r(lb(1),lb(2),j),1)
              END DO
            END DO
          END IF

        END DO

      END DO

    END DO

!   *** Release work storage ***

    IF (distributed_rs_grids) THEN
        CALL deallocate_matrix ( deltajp )
        CALL deallocate_matrix ( deltajp_rii )
        CALL deallocate_matrix ( deltajp_riii )
     END IF

    IF ( nthread > 1 ) THEN
      DEALLOCATE (lgrid%r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"lgrid%r")
    END IF

    DEALLOCATE (jpabt,jpabt_ii,jpabt_iii,workt,ntasks,tasks,tasks_local,ival,latom,&
        dist_ab,dab_local,asets,atasks,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"jpabt,workt,ntasks,"//&
        "tasks,tasks_local,ival,latom,dist_ab,dab_local,asets,atasks")

    CALL density_on_full_grid(pw_env,rs_rho,rho_rs,rho_gs,interp_section=interp_section,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_jrho_resp


! *****************************************************************************

  SUBROUTINE collocate_pgf_response_current(la_max,zeta,la_min,&
                                            lb_max,zetb,lb_min,&
                                            ra,rab,rab2,scale,ii,iii,&
                                            pab,pab_ii,pab_iii,&
                                            o1,o2,rsgrid,cube_info,l_info,&
                                            eps_rho_rspace,lgrid,ithread,map_consistent)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2, scale
    INTEGER, INTENT(IN)                      :: ii, iii
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, pab_ii, pab_iii
    INTEGER                                  :: o1, o2
    TYPE(realspace_grid_type), POINTER       :: rsgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(KIND=dp), INTENT(IN)                :: eps_rho_rspace
    TYPE(lgrid_type), OPTIONAL               :: lgrid
    INTEGER, INTENT(IN), OPTIONAL            :: ithread
    LOGICAL, INTENT(IN), OPTIONAL            :: map_consistent

    INTEGER :: cmax, coef_max, gridbounds(2,3), i, ico, icoef, ig, ithread_l, &
      jco, k, l, length, lx, lx_max, lxa, lxb, lxy, lxy_max, lxyz, lxyz_max, &
      lya, lyb, lza, lzb, offset, start
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    INTEGER, DIMENSION(:, :), POINTER        :: map
    INTEGER, POINTER                         :: ipzyx(:,:,:,:,:,:)
    LOGICAL                                  :: my_map_consistent
    REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, cutoff, f, pg, &
      prefactor, radius, rpg, xrrb, ya, yap, yb, ybp, yrrb, za, zap, zb, zbp, &
      zetp, zrrb
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rboffset, rbp, &
                                                roffset, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    REAL(KIND=dp), POINTER :: alpha(:,:), dpy(:,:), dpz(:,:), polx(:,:), &
      polx_ii(:,:), polx_iii(:,:), poly(:,:), poly_ii(:,:), poly_iii(:,:), &
      polz(:,:), polz_ii(:,:), polz_iii(:,:), pzyx(:), pzyx_ii(:), pzyx_iii(:)

!   ---------------------------------------------------------------------------

    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF

    ! use identical radii for integrate and collocate ?
    IF (PRESENT(map_consistent)) THEN
       my_map_consistent=map_consistent
    ELSE
       my_map_consistent=.FALSE.
    ENDIF

    zetp      = zeta + zetb
    f         = zetb/zetp
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)
    rb(:)     = ra(:)+rab(:)

   IF (my_map_consistent) THEN
       cutoff    = 1.0_dp
       prefactor = EXP(-zeta*f*rab2)
       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
                     zetp=zetp,eps=eps_rho_rspace,prefactor=prefactor,cutoff=cutoff)
       prefactor = scale*EXP(-zeta*f*rab2)
    ELSE
       cutoff    = 0.0_dp
       prefactor = scale*EXP(-zeta*f*rab2)
       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                                       zetp,eps_rho_rspace,prefactor,cutoff)
    ENDIF

    IF (radius .EQ. 0.0_dp ) THEN
      RETURN
    END IF

    coef_max=la_max+lb_max+1
!   *** properties of the grid ***
    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)

    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min,la_max,lb_min,lb_max,ithread_l,lx_max, &
                       lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),dp)*dr(:)

    rboffset(:) = rb(:) - REAL(cubecenter(:),dp)*dr(:)

!   *** a mapping so that the ig corresponds to the right grid point
    DO i=1,3
      IF ( rsgrid % perd ( i ) == 1 ) THEN
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         IF (start+length.GE.ub_cube(i)) EXIT
         start=start+length+1
        END DO
      ELSE
        ! this takes partial grid + border regions into account
        offset=MODULO(cubecenter(i),ng(i))+rsgrid%lb(i)
        offset=offset-rsgrid%lb_local(i)+1
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO

!   *** initialise the p terms and loop logic
    lxyz=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             pzyx(lxyz)=prefactor*pab(o1+ico,o2+jco)
             pzyx_ii(lxyz)=-prefactor*pab_ii(o1+ico,o2+jco)
             pzyx_iii(lxyz)=prefactor*pab_iii(o1+ico,o2+jco)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,dp)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      zrrb = REAL(ig,dp)*dr(3) - rboffset(3)
      DO lza=0,la_max
       zbp=1.0_dp
       DO lzb=0,lb_max
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
             IF(ii==3) THEN
               polz_ii(lxyz,ig)=dpz(lzb,lza)*zrrb
             ELSE
               polz_ii(lxyz,ig)=dpz(lzb,lza)
             END IF
             IF(iii==3) THEN
               polz_iii(lxyz,ig)=dpz(lzb,lza)*zrrb
             ELSE
               polz_iii(lxyz,ig)=dpz(lzb,lza)
             END IF
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

    DO ig=lb_cube(2),ub_cube(2)
      rpg = REAL(ig,dp)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      yrrb =  REAL(ig,dp)*dr(2) - rboffset(2)
      DO lya=0,la_max
       ybp=1.0_dp
       DO lyb=0,lb_max
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          lxy=lxy+1
          poly(lxy,ig)=dpy(lyb,lya)
          IF(ii==2) THEN
            poly_ii(lxy,ig)=dpy(lyb,lya)*yrrb
          ELSE
            poly_ii(lxy,ig)=dpy(lyb,lya)
          END IF
          IF(iii==2) THEN
            poly_iii(lxy,ig)=dpy(lyb,lya)*yrrb
          ELSE
            poly_iii(lxy,ig)=dpy(lyb,lya)
          END IF
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

!   *** make the alpha matrix ***
    alpha(:,:)=0.0_dp
    lx=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       lx=lx+1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO

    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,dp)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      xrrb = REAL(ig,dp)*dr(1) - rboffset(1)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         IF(ii==1) THEN
           polx_ii(icoef,ig)=pg*xrrb
         ELSE
           polx_ii(icoef,ig)=pg
         END IF
         IF(iii==1) THEN
           polx_iii(icoef,ig)=pg*xrrb
         ELSE
           polx_iii(icoef,ig)=pg
         END IF
         pg=pg*(rpg)
      ENDDO
    ENDDO

    IF ( PRESENT ( lgrid ) ) THEN
      ig = lgrid%ldim * ithread_l + 1
      CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax),&
              polz(1,-cmax),lgrid%r(ig),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
      CALL collocate_core(pzyx_ii(1),polx_ii(1,-cmax),poly_ii(1,-cmax),&
              polz_ii(1,-cmax),lgrid%r(ig),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
      CALL collocate_core(pzyx_iii(1),polx_iii(1,-cmax),poly_iii(1,-cmax),&
              polz_iii(1,-cmax),lgrid%r(ig),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
    ELSE
!   *** do the loop over the grid
!   notice this is not the same as critical or so, since we may have
!   several different rsgrids that all use the same function
!   I guess we need a flush (of the grid) to guarantee that we are working
!   on an uptodate copy of the grid
      CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax),&
              polz(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
      CALL collocate_core(pzyx_ii(1),polx_ii(1,-cmax),poly_ii(1,-cmax),&
              polz_ii(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
      CALL collocate_core(pzyx_iii(1),polx_iii(1,-cmax),poly_iii(1,-cmax),&
              polz_iii(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
    END IF

  END SUBROUTINE  collocate_pgf_response_current

END MODULE qs_linres_nmr_current
