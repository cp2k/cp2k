!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate the 3 and 2 center ERI's needed in the RI 
!>        approximation using libint
!> \par History
!>      08.2013 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
MODULE mp2_ri_libint

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type,&
                                             init_aux_basis_set
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_files,                        ONLY: get_unit_number
  USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_release,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_energy_potential,            ONLY: coulomb4
  USE hfx_libint_wrapper_types,        ONLY: lib_int
  USE hfx_pair_list_methods,           ONLY: build_pair_list_mp2
  USE hfx_screening_methods,           ONLY: calc_pair_dist_radii,&
                                             calc_screening_functions
  USE hfx_types,                       ONLY: &
       hfx_basis_info_type, hfx_basis_type, hfx_create_neighbor_cells, &
       hfx_general_type, hfx_load_balance_type, hfx_memory_type, &
       hfx_pgf_list, hfx_pgf_product_list, hfx_potential_type, &
       hfx_screen_coeff_type, hfx_screening_type, hfx_type, log_zero, &
       pair_set_list_type
  USE input_constants,                 ONLY: do_mp2_potential_TShPSC,&
                                             use_ri_aux_basis_set
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_sum
  USE mp2_types,                       ONLY: init_TShPSC_lmax,&
                                             mp2_biel_type,&
                                             mp2_type,&
                                             pair_list_type_mp2
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset,&
                                             nso
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE t_sh_p_s_c,                      ONLY: free_C0,&
                                             init
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC ::  libint_ri_mp2, read_RI_basis_set, release_RI_basis_set

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_libint'

!***
  
  CONTAINS

  SUBROUTINE libint_ri_mp2(dimen,RI_dimen,occupied,natom,mp2_biel,mp2_env,C,&
                           kind_of,basis_parameter,particle_set,&     
                           RI_basis_parameter,RI_basis_info,basis_S0,RI_index_table,&                                       
                           qs_env,hfx_section,para_env, &                        
                           unit_nr,Lai,error)                                            
    INTEGER                                  :: dimen, RI_dimen, occupied, &
                                                natom
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C
    INTEGER, DIMENSION(:)                    :: kind_of
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: RI_index_table
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: unit_nr
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Lai
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'libint_ri_mp2', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: due = 2.D+00, zero = 0.D+00

    INTEGER                                  :: handle, nkind
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.              

    ! Get the RI basis set and store in to a nice form 
    IF(.NOT.(ASSOCIATED(RI_basis_parameter))) THEN
      IF(ALLOCATED(RI_index_table)) DEALLOCATE(RI_index_table)
      IF(ASSOCIATED(basis_S0)) DEALLOCATE(basis_S0)
      CALL read_RI_basis_set(qs_env,RI_basis_parameter,RI_basis_info,&
                             natom,nkind,kind_of,RI_index_table,RI_dimen,&
                             basis_S0,error)
    END IF
 
    CALL  calc_lai_libint(mp2_env,qs_env,hfx_section,para_env,&
                          mp2_biel,dimen,C,occupied,&
                          RI_basis_parameter,RI_basis_info,RI_index_table,RI_dimen,basis_S0,&
                          Lai,error)   

    CALL timestop(handle)

  END SUBROUTINE  libint_ri_mp2

! *****************************************************************************
!> \brief Read the auxiliary basis set for RI approxiamtion
!> \par History
!>      08.2013 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
  SUBROUTINE read_RI_basis_set(qs_env,RI_basis_parameter,RI_basis_info,&
                               natom,nkind,kind_of,RI_index_table,RI_dimen,&
                               basis_S0,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    INTEGER                                  :: natom, nkind
    INTEGER, DIMENSION(:)                    :: kind_of
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: RI_index_table
    INTEGER                                  :: RI_dimen
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_RI_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER :: co_counter, counter, i, iatom, ikind, ipgf, iset, j, k, la, &
      max_am_kind, max_coeff, max_nsgfl, max_pgf, max_pgf_kind, max_set, &
      nl_count, nset, nseta, offset_a, offset_a1, s_offset_nl_a, sgfa, &
      so_counter, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, nl_a
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_a
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a

    failure = .FALSE.

    NULLIFY(RI_basis_parameter)

    NULLIFY(atomic_kind_set)
    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, error=error)

    nkind = SIZE(atomic_kind_set,1)
    ALLOCATE(RI_basis_parameter(nkind), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(basis_S0(nkind), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    max_set = 0
    DO ikind = 1,nkind
      NULLIFY(atom_kind)
      atom_kind => atomic_kind_set(ikind)
      IF(.NOT.(ASSOCIATED(atomic_kind_set(ikind)%ri_aux_basis_set))) THEN
        CALL stop_program(routineN,moduleN,__LINE__,"Initial RI auxiliary basis not specified.")
      END IF
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           ri_aux_basis_set=orb_basis_a)
      ! here we reset the initial RI basis such that we can
      ! work with non-normalized auxiliary basis functions
      orb_basis_a%gcc=1.0_dp
      orb_basis_a%norm_type = 1
      CALL init_aux_basis_set(orb_basis_a,error)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               maxsgf=RI_basis_info%max_sgf,&
                               maxnset=RI_basis_info%max_set,&
                               maxlgto=RI_basis_info%max_am,&
                               basis_set_id=use_ri_aux_basis_set)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, &
                             lmax=RI_basis_parameter(ikind)%lmax, &
                             lmin=RI_basis_parameter(ikind)%lmin, &
                             npgf=RI_basis_parameter(ikind)%npgf, &
                             nset=RI_basis_parameter(ikind)%nset, &
                             zet=RI_basis_parameter(ikind)%zet, &
                             nsgf_set=RI_basis_parameter(ikind)%nsgf, &
                             first_sgf=RI_basis_parameter(ikind)%first_sgf, &
                             sphi=RI_basis_parameter(ikind)%sphi, &
                             nsgf=RI_basis_parameter(ikind)%nsgf_total,&
                             l=RI_basis_parameter(ikind)%nl,&
                             nshell=RI_basis_parameter(ikind)%nshell,&
                             set_radius=RI_basis_parameter(ikind)%set_radius,&
                             pgf_radius=RI_basis_parameter(ikind)%pgf_radius,&
                             kind_radius=RI_basis_parameter(ikind)%kind_radius)

      max_set = MAX(max_set,RI_basis_parameter(ikind)%nset)
       
       basis_S0(ikind)%kind_radius=RI_basis_parameter(ikind)%kind_radius
       basis_S0(ikind)%nset=1
       basis_S0(ikind)%nsgf_total=1
       ALLOCATE(basis_S0(ikind)%set_radius(1))
       basis_S0(ikind)%set_radius(1)=RI_basis_parameter(ikind)%kind_radius
       ALLOCATE(basis_S0(ikind)%lmax(1))
       basis_S0(ikind)%lmax(1)=0
       ALLOCATE(basis_S0(ikind)%lmin(1))
       basis_S0(ikind)%lmin(1)=0
       ALLOCATE(basis_S0(ikind)%npgf(1))
       basis_S0(ikind)%npgf(1)=1
       ALLOCATE(basis_S0(ikind)%nsgf(1))
       basis_S0(ikind)%nsgf(1)=1
       ALLOCATE(basis_S0(ikind)%nshell(1))
       basis_S0(ikind)%nshell(1)=1
       ALLOCATE(basis_S0(ikind)%pgf_radius(1,1))
       basis_S0(ikind)%pgf_radius(1,1)=RI_basis_parameter(ikind)%kind_radius
       ALLOCATE(basis_S0(ikind)%sphi(1,1))
       basis_S0(ikind)%sphi(1,1)=1.0_dp
       ALLOCATE(basis_S0(ikind)%zet(1,1))
       basis_S0(ikind)%zet(1,1)=0.0_dp
       ALLOCATE(basis_S0(ikind)%first_sgf(1,1))
       basis_S0(ikind)%first_sgf(1,1)=1
       ALLOCATE(basis_S0(ikind)%nl(1,1))
       basis_S0(ikind)%nl(1,1)=0

       ALLOCATE(basis_S0(ikind)%nsgfl(0:0,1))
       basis_S0(ikind)%nsgfl=1
       ALLOCATE(basis_S0(ikind)%sphi_ext(1,0:0,1,1))
       basis_S0(ikind)%sphi_ext(1,0,1,1)=1.0_dp

    END DO
    RI_basis_info%max_set=max_set

    DO ikind = 1,nkind
      ALLOCATE(RI_basis_parameter(ikind)%nsgfl(0:RI_basis_info%max_am,max_set))
      RI_basis_parameter(ikind)%nsgfl = 0
      nset = RI_basis_parameter(ikind)%nset
      nshell => RI_basis_parameter(ikind)%nshell
      DO iset = 1,nset
        DO i=0,RI_basis_info%max_am
          nl_count = 0
          DO j=1,nshell(iset)
            IF( RI_basis_parameter(ikind)%nl(j,iset) == i) nl_count = nl_count + 1
          END DO
          RI_basis_parameter(ikind)%nsgfl(i,iset) = nl_count
        END DO
      END DO
    END DO

    max_nsgfl = 0
    max_pgf = 0
    DO ikind = 1,nkind
      max_coeff = 0
      max_am_kind = 0
      max_pgf_kind = 0
      npgfa => RI_basis_parameter(ikind)%npgf
      nseta = RI_basis_parameter(ikind)%nset
      nl_a => RI_basis_parameter(ikind)%nsgfl
      la_max => RI_basis_parameter(ikind)%lmax
      la_min => RI_basis_parameter(ikind)%lmin
      DO iset = 1,nseta
        max_pgf_kind = MAX(max_pgf_kind,npgfa(iset))
        max_pgf = MAX(max_pgf, npgfa(iset))
        DO la = la_min(iset), la_max(iset)
          max_nsgfl = MAX(max_nsgfl, nl_a(la,iset))
          max_coeff = MAX(max_coeff, nso(la)*nl_a(la,iset)*nco(la))
          max_am_kind = MAX(max_am_kind, la)
        END DO
      END DO
      ALLOCATE(RI_basis_parameter(ikind)%sphi_ext(max_coeff,0:max_am_kind, max_pgf_kind, nseta), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      RI_basis_parameter(ikind)%sphi_ext = 0.0_dp
    END DO 

    DO ikind = 1,nkind
      sphi_a => RI_basis_parameter(ikind)%sphi
      nseta = RI_basis_parameter(ikind)%nset
      la_max => RI_basis_parameter(ikind)%lmax
      la_min => RI_basis_parameter(ikind)%lmin
      npgfa => RI_basis_parameter(ikind)%npgf
      first_sgfa => RI_basis_parameter(ikind)%first_sgf
      nl_a => RI_basis_parameter(ikind)%nsgfl
      DO iset = 1,nseta
        sgfa = first_sgfa(1,iset)
        DO ipgf = 1,npgfa(iset)
          offset_a1 = (ipgf-1)*ncoset(la_max(iset))
          s_offset_nl_a = 0
          DO la = la_min(iset),la_max(iset)
            offset_a = offset_a1 + ncoset(la-1)
            co_counter = 0
            co_counter = co_counter + 1
            so_counter = 0
            DO k = sgfa+s_offset_nl_a,sgfa+s_offset_nl_a+nso(la)*nl_a(la,iset)-1
              DO i = offset_a+1,offset_a+nco(la)
                so_counter = so_counter + 1
                RI_basis_parameter(ikind)%sphi_ext(so_counter,la,ipgf,iset) = sphi_a(i,k)
              END DO
            END DO
            s_offset_nl_a = s_offset_nl_a + nso(la) * (nl_a(la,iset))
          END DO
        END DO
      END DO
    END DO

    ALLOCATE(RI_index_table(natom,max_set))
    RI_index_table=-HUGE(0)
    counter=0
    RI_dimen=0
    DO iatom=1, natom
      ikind=kind_of(iatom)
      nset = RI_basis_parameter(ikind)%nset
      DO iset=1, nset
        RI_index_table(iatom,iset)=counter+1
        counter=counter+RI_basis_parameter(ikind)%nsgf(iset)
        RI_dimen=RI_dimen+RI_basis_parameter(ikind)%nsgf(iset)
      END DO
    END DO

  END SUBROUTINE read_RI_basis_set

! *****************************************************************************
!> \brief Release the auxiliary basis set for RI approxiamtion (to be used
!>        only in the case of basis optimization)
!> \par History
!>      08.2013 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
  SUBROUTINE release_RI_basis_set(RI_basis_parameter,basis_S0,error)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter, basis_S0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_RI_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ! RI basis 
    DO i=1,SIZE(RI_basis_parameter)
       DEALLOCATE(RI_basis_parameter(i)%nsgfl,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(RI_basis_parameter(i)%sphi_ext,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(RI_basis_parameter,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! S0 basis
    DO i=1, SIZE(basis_S0)
       DEALLOCATE(basis_S0(i)%set_radius,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%lmax,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%lmin,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%npgf,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%nsgf,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%nshell,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%pgf_radius,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%sphi,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%zet,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%first_sgf,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%nl,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%nsgfl,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_S0(i)%sphi_ext,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(basis_S0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE release_RI_basis_set 

! *****************************************************************************
!> \brief 
!> \param qs_env  
!> \param energy
!> \param hfx_section input_section HFX
!> \param para_env
!> \param geometry_did_change flag that indicates we have to recalc integrals
!> \param irep Index for the HFX replica
!> \param distribute_fock_matrix Flag that indicates whether to communicate the
!>        new fock matrix or not
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2013 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
  SUBROUTINE calc_lai_libint(mp2_env,qs_env,hfx_section,para_env,& 
                             mp2_biel,dimen,C,occupied,&
                             RI_basis_parameter,RI_basis_info,RI_index_table,RI_dimen,basis_S0,&
                             Lai,error)

    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mp2_biel_type)                      :: mp2_biel
    INTEGER                                  :: dimen
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C
    INTEGER                                  :: occupied
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: RI_index_table
    INTEGER                                  :: RI_dimen
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Lai
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_lai_libint', &
      routineP = moduleN//':'//routineN

    INTEGER :: bits_max_val, cache_size, counter_L_blocks, handle, i, &
      i_list_kl, i_set_list_kl, i_set_list_kl_start, i_set_list_kl_stop, &
      i_thread, iatom, iatom_end, iatom_start, iiB, ikind, info_chol, irep, &
      iset, jatom, jatom_end, jatom_start, jjB, jkind, jset, katom, &
      katom_end, katom_start, kkB, kkind, kset, kset_start, L_B_i_end, &
      L_B_i_start, L_B_k_end, L_B_k_start, l_max, latom, latom_end, &
      latom_start, lkind, llB, lset, max_am, max_pgf, max_set, n_threads, &
      natom, ncob, ncos_max, nints, nkind, nneighbors, nseta, nsetb, nsetc, &
      nsetd, nsgf_max, nspins, orb_k_end, orb_k_start, orb_l_end, &
      orb_l_start, primitive_counter
    INTEGER :: sgfb, sphi_a_u1, sphi_a_u2, sphi_a_u3, sphi_b_u1, sphi_b_u2, &
      sphi_b_u3, sphi_c_u1, sphi_c_u2, sphi_c_u3, sphi_d_u1, sphi_d_u2, &
      sphi_d_u3, stat, unit_id, virtual
    INTEGER(int_8) :: estimate_to_store_int, max_val_memory, mem_eris, &
      mem_eris_disk, mem_max_val, n_processes, ncpu, neris_disk, &
      neris_incore, neris_onthefly, neris_tmp, neris_total, nprim_ints, &
      shm_neris_disk, shm_neris_incore, shm_neris_onthefly, shm_neris_total, &
      shm_nprim_ints, shm_stor_count_int_disk, shm_stor_count_max_val, &
      shm_storage_counter_integrals
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, last_sgf_global, &
                                                nimages
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd, shm_block_offset
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, nsgfl_a, nsgfl_b, &
                                                nsgfl_c, nsgfl_d, &
                                                shm_atomic_block_offset
    INTEGER, DIMENSION(:, :), POINTER, SAVE  :: shm_is_assoc_atomic_block
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: shm_set_offset
    INTEGER, SAVE                            :: shm_number_of_p_entries
    LOGICAL                                  :: do_periodic, failure, &
                                                treat_lsd_in_core
    LOGICAL, DIMENSION(:, :), POINTER        :: shm_atomic_pair_list
    REAL(KIND=dp) :: cartesian_estimate, coeffs_kind_max0, &
      compression_factor, compression_factor_disk, eps_schwarz, eps_storage, &
      hf_fraction, ln_10, log10_eps_schwarz, log10_pmax, max_contraction_val, &
      pmax_atom, pmax_entry, ra(3), rab2, rb(3), rc(3), rcd2, rd(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ee_buffer1, ee_buffer2, &
                                                ee_primitives_tmp, ee_work, &
                                                ee_work2, primitive_integrals
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: L_block, L_full_matrix
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BI1, MNRS
    REAL(KIND=dp), DIMENSION(:), POINTER     :: p_work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: max_contraction, &
                                                shm_pmax_block, sphi_b, zeta, &
                                                zetb, zetc, zetd
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: sphi_a_ext_set, &
                                                sphi_b_ext_set, &
                                                sphi_c_ext_set, sphi_d_ext_set
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: sphi_a_ext, sphi_b_ext, &
                                                sphi_c_ext, sphi_d_ext
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_L
    TYPE(cp_fm_type), POINTER                :: fm_matrix_L
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hfx_basis_info_type), POINTER       :: basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_general_type)                   :: general_parameter
    TYPE(hfx_load_balance_type), POINTER     :: load_balance_parameter
    TYPE(hfx_memory_type), POINTER           :: memory_parameter
    TYPE(hfx_pgf_list), ALLOCATABLE, &
      DIMENSION(:)                           :: pgf_list_ij, pgf_list_kl
    TYPE(hfx_pgf_product_list), &
      ALLOCATABLE, DIMENSION(:)              :: pgf_product_list
    TYPE(hfx_potential_type)                 :: mp2_potential_parameter, &
                                                potential_parameter
    TYPE(hfx_screen_coeff_type), &
      ALLOCATABLE, DIMENSION(:, :), TARGET   :: radii_pgf_large, &
                                                screen_coeffs_pgf_large
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: screen_coeffs_kind, tmp_R_1, &
                                                tmp_R_2, tmp_screen_pgf1, &
                                                tmp_screen_pgf2
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: screen_coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: radii_pgf, screen_coeffs_pgf
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data, &
                                                shm_master_x_data
    TYPE(lib_int)                            :: private_lib
    TYPE(pair_list_type_mp2)                 :: list_kl
    TYPE(pair_set_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: set_list_kl
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)

    log10_eps_schwarz=log_zero

    irep=1

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, error=error)

    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(qs_env%x_data(1,1)%basis_parameter(ikind)%lmax))
      l_max = MAX(l_max,MAXVAL(RI_basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max
    CALL init_md_ftable(l_max)

    IF(mp2_env%potential_parameter%potential_type == do_mp2_potential_TShPSC) THEN
      IF(l_max>init_TShPSC_lmax) THEN
        IF( para_env%mepos == 0 ) THEN
          unit_id = get_unit_number()
          OPEN(unit_id,FILE=mp2_env%potential_parameter%filename)
        END IF
        CALL init(l_max, unit_id, para_env%mepos, para_env%group)
        IF(para_env%mepos == 0 ) THEN
          CLOSE(unit_id)
        END IF
        init_TShPSC_lmax = l_max
      END IF
      mp2_potential_parameter%cutoff_radius=mp2_env%potential_parameter%truncation_radius/2.0_dp
    END IF

    mp2_potential_parameter%potential_type=mp2_env%potential_parameter%potential_type

    n_threads = 1

    shm_neris_total = 0
    shm_nprim_ints = 0
    shm_neris_onthefly = 0
    shm_storage_counter_integrals = 0
    shm_stor_count_int_disk = 0
    shm_neris_incore = 0
    shm_neris_disk = 0
    shm_stor_count_max_val = 0

    ln_10 = LOG(10.0_dp)
    i_thread = 0

    actual_x_data => qs_env%x_data(irep, i_thread + 1)

    shm_master_x_data => qs_env%x_data(irep,1)

    do_periodic = actual_x_data%periodic_parameter%do_periodic

    IF( do_periodic ) THEN
      ! ** Rebuild neighbor lists in case the cell has changed (i.e. NPT MD)
      actual_x_data%periodic_parameter%number_of_shells = actual_x_data%periodic_parameter%mode
      CALL hfx_create_neighbor_cells(actual_x_data, actual_x_data%periodic_parameter%number_of_shells_from_input,&
                                     cell, i_thread, &
                                     error)
    END IF

    screening_parameter    = actual_x_data%screening_parameter
    potential_parameter    = actual_x_data%potential_parameter

    general_parameter      = actual_x_data%general_parameter
    load_balance_parameter => actual_x_data%load_balance_parameter
    memory_parameter       => actual_x_data%memory_parameter

    cache_size = memory_parameter%cache_size
    bits_max_val = memory_parameter%bits_max_val

    basis_parameter   => actual_x_data%basis_parameter
    basis_info        => actual_x_data%basis_info

    treat_lsd_in_core = general_parameter%treat_lsd_in_core

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads

    !! initalize some counters
    neris_total = 0_int_8
    neris_incore = 0_int_8
    neris_disk = 0_int_8
    neris_onthefly = 0_int_8
    mem_eris = 0_int_8
    mem_eris_disk = 0_int_8
    mem_max_val = 0_int_8
    compression_factor = 0.0_dp
    compression_factor_disk = 0.0_dp
    nprim_ints = 0_int_8
    neris_tmp = 0_int_8
    max_val_memory = 1_int_8

    max_am = MAX(basis_info%max_am,RI_basis_info%max_am)

    max_set = MAX(basis_info%max_set,RI_basis_info%max_set)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)

    natom = SIZE(particle_set,1)

    ALLOCATE(kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of)

    !! precompute maximum nco and allocate scratch
    ncos_max=0
    nsgf_max=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      nseta = basis_parameter(ikind)%nset
      npgfa => basis_parameter(ikind)%npgf
      la_max => basis_parameter(ikind)%lmax
      nsgfa => basis_parameter(ikind)%nsgf
      DO iset = 1, nseta
         ncos_max = MAX(ncos_max,ncoset(la_max(iset)))
         nsgf_max = MAX(nsgf_max,nsgfa(iset))
      ENDDO
      nseta = RI_basis_parameter(ikind)%nset
      npgfa => RI_basis_parameter(ikind)%npgf
      la_max => RI_basis_parameter(ikind)%lmax
      nsgfa => RI_basis_parameter(ikind)%nsgf
      DO iset = 1, nseta
         ncos_max = MAX(ncos_max,ncoset(la_max(iset)))
         nsgf_max = MAX(nsgf_max,nsgfa(iset))
      ENDDO
    ENDDO
    !! Allocate the arrays for the integrals.
    ALLOCATE(primitive_integrals(nsgf_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_integrals = 0.0_dp

    ALLOCATE(ee_work(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_work2(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_buffer1(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_buffer2(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_primitives_tmp(nsgf_max**4),STAT=stat) ! XXXXX could be wrong
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    nspins = qs_env%dft_control%nspins

    ALLOCATE(max_contraction(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    max_contraction=0.0_dp
    max_pgf = 0
    DO jatom=1,natom
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      nsetb = basis_parameter(jkind)%nset
      npgfb => basis_parameter(jkind)%npgf
      first_sgfb => basis_parameter(jkind)%first_sgf
      sphi_b => basis_parameter(jkind)%sphi
      nsgfb => basis_parameter(jkind)%nsgf
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        max_contraction(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
        max_pgf = MAX(max_pgf,npgfb(jset))
      ENDDO
    ENDDO

    ! ** Allocate buffers for pgf_lists
    nneighbors = SIZE(actual_x_data%neighbor_cells)
    ALLOCATE(pgf_list_ij(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pgf_list_kl(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pgf_product_list(nneighbors**3), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nimages(max_pgf**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,max_pgf**2
      ALLOCATE(pgf_list_ij(i)%image_list(nneighbors), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pgf_list_kl(i)%image_list(nneighbors), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO

    ! ** Set pointers
    shm_number_of_p_entries = shm_master_x_data%number_of_p_entries
    shm_is_assoc_atomic_block => shm_master_x_data%is_assoc_atomic_block
    shm_atomic_block_offset => shm_master_x_data%atomic_block_offset
    shm_set_offset => shm_master_x_data%set_offset
    shm_block_offset => shm_master_x_data%block_offset

    !!! Skipped part

    !! Get screening parameter
    eps_schwarz = screening_parameter%eps_schwarz
    IF( eps_schwarz <= 0.0_dp) THEN
      log10_eps_schwarz = log_zero
    ELSE
      log10_eps_schwarz = LOG10(eps_schwarz)
    END IF
    !! get storage epsilon
    eps_storage = eps_schwarz*memory_parameter%eps_storage_scaling

    !! If we have a hybrid functional, we may need only a fraction of exact exchange
    hf_fraction = general_parameter%fraction

    !! The number of integrals that fit into the given MAX_MEMORY

    !! Parameters related to the potential 1/r, erf(wr)/r, erfc(wr/r)
    potential_parameter = actual_x_data%potential_parameter

    logger => cp_error_get_logger(error)

    private_lib = actual_x_data%lib

    !!! Helper array to map local basis function indeces to global ones
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+ basis_parameter(ikind)%nsgf_total
    END DO

    IF( .NOT. shm_master_x_data%screen_funct_is_initialized ) THEN
      CALL calc_pair_dist_radii(qs_env, basis_parameter,&
                                shm_master_x_data%pair_dist_radii_pgf, max_set, max_pgf, eps_schwarz,&
                                n_threads, i_thread, error)
      CALL calc_screening_functions(qs_env, basis_parameter, private_lib, shm_master_x_data%potential_parameter,&
                                    shm_master_x_data%screen_funct_coeffs_set,&
                                    shm_master_x_data%screen_funct_coeffs_kind, &
                                    shm_master_x_data%screen_funct_coeffs_pgf, &
                                    shm_master_x_data%pair_dist_radii_pgf,&
                                    max_set, max_pgf, n_threads, i_thread, p_work, error)

      shm_master_x_data%screen_funct_is_initialized = .TRUE.
   END IF

    screen_coeffs_set   => shm_master_x_data%screen_funct_coeffs_set
    screen_coeffs_kind  => shm_master_x_data%screen_funct_coeffs_kind
    screen_coeffs_pgf   => shm_master_x_data%screen_funct_coeffs_pgf
    radii_pgf           => shm_master_x_data%pair_dist_radii_pgf

    ALLOCATE(radii_pgf_large(SIZE(radii_pgf,1),SIZE(radii_pgf,2)))
    ALLOCATE(screen_coeffs_pgf_large(SIZE(screen_coeffs_pgf,1),SIZE(screen_coeffs_pgf,2)))
    DO iiB=1, SIZE(radii_pgf,1)
      DO jjB=1, SIZE(radii_pgf,2)
        radii_pgf_large(iiB,jjB)%x=100_dp
      END DO
    END DO
    DO iiB=1, SIZE(screen_coeffs_pgf,1)
      DO jjB=1, SIZE(screen_coeffs_pgf,2)
        screen_coeffs_pgf_large(iiB,jjB)%x=5000_dp
      END DO
    END DO
    tmp_R_1 => radii_pgf_large(:,:)
    tmp_R_2 => radii_pgf_large(:,:)
    tmp_screen_pgf1 => screen_coeffs_pgf_large(:,:)
    tmp_screen_pgf2 => screen_coeffs_pgf_large(:,:)


    ! start computing the L matrix
    ALLOCATE(L_full_matrix(RI_dimen,RI_dimen), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    L_full_matrix=0.0_dp

    counter_L_blocks=0
    DO iatom=1, natom
      jatom=iatom

      ikind=kind_of(iatom)
      jkind=kind_of(jatom)
      ra=particle_set(iatom)%r(:)
      rb=particle_set(jatom)%r(:)
      !rab2=(ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
      rab2=0.0_dp

      la_max => RI_basis_parameter(ikind)%lmax
      la_min => RI_basis_parameter(ikind)%lmin
      npgfa => RI_basis_parameter(ikind)%npgf
      nseta = RI_basis_parameter(ikind)%nset
      zeta => RI_basis_parameter(ikind)%zet
      nsgfa => RI_basis_parameter(ikind)%nsgf
      sphi_a_ext => RI_basis_parameter(ikind)%sphi_ext(:,:,:,:)
      nsgfl_a => RI_basis_parameter(ikind)%nsgfl
      sphi_a_u1=UBOUND(sphi_a_ext,1)
      sphi_a_u2=UBOUND(sphi_a_ext,2)
      sphi_a_u3=UBOUND(sphi_a_ext,3)

      lb_max => basis_S0(jkind)%lmax
      lb_min => basis_S0(jkind)%lmin
      npgfb => basis_S0(jkind)%npgf
      nsetb = basis_S0(jkind)%nset
      zetb => basis_S0(jkind)%zet
      nsgfb => basis_S0(jkind)%nsgf
      sphi_b_ext => basis_S0(jkind)%sphi_ext(:,:,:,:)
      nsgfl_b => basis_S0(jkind)%nsgfl
      sphi_b_u1=UBOUND(sphi_b_ext,1)
      sphi_b_u2=UBOUND(sphi_b_ext,2)
      sphi_b_u3=UBOUND(sphi_b_ext,3)

      DO katom=iatom, natom
        latom=katom

        kkind=kind_of(katom)
        lkind=kind_of(latom)
        rc=particle_set(katom)%r(:)
        rd=particle_set(latom)%r(:)
        !rcd2=(rc(1)-rd(1))**2+(rc(2)-rd(2))**2+(rc(3)-rd(3))**2
        rcd2=0.0_dp

        !pmax_atom = 0.0_dp
        !screen_kind_ij = screen_coeffs_kind(jkind,ikind)%x(1)*rab2+&
        !                 screen_coeffs_kind(jkind,ikind)%x(2)
        !screen_kind_kl = screen_coeffs_kind(lkind,kkind)%x(1)*rcd2+&
        !                 screen_coeffs_kind(lkind,kkind)%x(2)
        !IF( max_val1 + screen_kind_kl + pmax_atom < log10_eps_schwarz) CYCLE
        !IF( screen_kind_ij + screen_kind_kl + pmax_atom < log10_eps_schwarz ) CYCLE

        lc_max => RI_basis_parameter(kkind)%lmax
        lc_min => RI_basis_parameter(kkind)%lmin
        npgfc => RI_basis_parameter(kkind)%npgf
        nsetc = RI_basis_parameter(kkind)%nset
        zetc => RI_basis_parameter(kkind)%zet
        nsgfc => RI_basis_parameter(kkind)%nsgf
        sphi_c_ext => RI_basis_parameter(kkind)%sphi_ext(:,:,:,:)
        nsgfl_c => RI_basis_parameter(kkind)%nsgfl
        sphi_c_u1=UBOUND(sphi_c_ext,1)
        sphi_c_u2=UBOUND(sphi_c_ext,2)
        sphi_c_u3=UBOUND(sphi_c_ext,3)

        ld_max => basis_S0(lkind)%lmax
        ld_min => basis_S0(lkind)%lmin
        npgfd => basis_S0(lkind)%npgf
        nsetd = RI_basis_parameter(lkind)%nset
        zetd => basis_S0(lkind)%zet
        nsgfd => basis_S0(lkind)%nsgf
        sphi_d_ext => basis_S0(lkind)%sphi_ext(:,:,:,:)
        nsgfl_d => basis_S0(lkind)%nsgfl
        sphi_d_u1=UBOUND(sphi_d_ext,1)
        sphi_d_u2=UBOUND(sphi_d_ext,2)
        sphi_d_u3=UBOUND(sphi_d_ext,3)

        jset=1
        lset=1
        DO iset=1, nseta

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          !max_val1 = screen_coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
          !           screen_coeffs_set(jset,iset,jkind,ikind)%x(2)
          sphi_a_ext_set => sphi_a_ext(:,:,:,iset)
          sphi_b_ext_set => sphi_b_ext(:,:,:,jset)

          L_B_i_start=RI_index_table(iatom,iset)
          L_B_i_end=RI_index_table(iatom,iset)+nsgfa(iset)-1

          kset_start=1
          IF(iatom==katom) kset_start=iset
          DO kset=kset_start, nsetc

            counter_L_blocks=counter_L_blocks+1
            IF(MOD(counter_L_blocks,para_env%num_pe)/=para_env%mepos) CYCLE
            !max_val2_set = (screen_coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
            !                screen_coeffs_set(lset,kset,lkind,kkind)%x(2) )
            !max_val2 = max_val1 + max_val2_set
            !!! Near field screening
            !IF(max_val2 + pmax_atom <log10_eps_schwarz) CYCLE

            sphi_c_ext_set => sphi_c_ext(:,:,:,kset)
            sphi_d_ext_set => sphi_d_ext(:,:,:,lset)

            L_B_k_start=RI_index_table(katom,kset)
            L_B_k_end=RI_index_table(katom,kset)+nsgfc(kset)-1

            pmax_entry = 0.0_dp
            log10_pmax = pmax_entry
            log10_eps_schwarz=log_zero
            !max_val2 = max_val2 + log10_pmax
            !IF(max_val2<log10_eps_schwarz) CYCLE
            !pmax_entry = EXP(log10_pmax*ln_10)  

            max_contraction_val = 1.0000_dp
            ! tmp_R_1 => radii_pgf(:,:,1,1,1,1)
            ! tmp_R_2 => radii_pgf(:,:,1,1,1,1)
            ! tmp_screen_pgf1 => screen_coeffs_pgf(:,:,1,1,1,1)
            ! tmp_screen_pgf2 => screen_coeffs_pgf(:,:,1,1,1,1)

            ALLOCATE(L_block(nsgfa(iset),nsgfc(kset)))

            CALL coulomb4(private_lib, ra, rb, rc, rd, npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset), &
                          la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                          lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                          nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset),&
                          sphi_a_u1,sphi_a_u2,sphi_a_u3,&
                          sphi_b_u1,sphi_b_u2,sphi_b_u3,&
                          sphi_c_u1,sphi_c_u2,sphi_c_u3,&
                          sphi_d_u1,sphi_d_u2,sphi_d_u3,&
                          zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset),&
                          zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                          primitive_integrals,&
                          mp2_potential_parameter, &
                          actual_x_data%neighbor_cells, screen_coeffs_set(1,1,1,1)%x,&
                          screen_coeffs_set(1,1,1,1)%x, eps_schwarz, &
                          max_contraction_val, cartesian_estimate,  cell, neris_tmp,&
                          log10_pmax, log10_eps_schwarz, &
                          tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2,&
                          pgf_list_ij, pgf_list_kl, pgf_product_list, &
                          nsgfl_a(:,iset), nsgfl_b(:,jset),&
                          nsgfl_c(:,kset), nsgfl_d(:,lset),&
                          sphi_a_ext_set,&
                          sphi_b_ext_set,&
                          sphi_c_ext_set,&
                          sphi_d_ext_set,&
                          ee_work,ee_work2,ee_buffer1,ee_buffer2, ee_primitives_tmp,&
                          nimages,do_periodic, p_work)

            primitive_counter=0
            DO llB=1, nsgfd(lset)
              DO kkB=1, nsgfc(kset)
                DO jjB=1, nsgfb(jset)
                  DO iiB=1, nsgfa(iset)
                    primitive_counter=primitive_counter+1
                    L_block(iiB,kkB)=primitive_integrals(primitive_counter)
                    ! WRITE(1001,*) L_B_i_start+iiB-1, L_B_k_start+kkB-1, L_block(iiB,kkB)
                  END DO
                END DO
              END DO
            END DO

            L_full_matrix(L_B_i_start:L_B_i_end,L_B_k_start:L_B_k_end)=L_block
            L_full_matrix(L_B_k_start:L_B_k_end,L_B_i_start:L_B_i_end)=TRANSPOSE(L_block)

            DEALLOCATE(L_block)

          END DO ! kset
        END DO ! iset

      END DO !katom
    END DO !iatom

    ! now create a fm_matrix for the L matrix
    CALL mp_sum(L_full_matrix,para_env%group)

    ! create a sub blacs_env
    NULLIFY(blacs_env)
    CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env, error=error)  
    
    NULLIFY(fm_matrix_L)
    NULLIFY(fm_struct_L)
    CALL cp_fm_struct_create(fm_struct_L,context=blacs_env,nrow_global=RI_dimen,&
                             ncol_global=RI_dimen,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix_L,fm_struct_L,name="fm_matrix_L",error=error)
    CALL cp_fm_struct_release(fm_struct_L,error=error)
    CALL cp_blacs_env_release(blacs_env, error)

    CALL cp_fm_set_submatrix(fm=fm_matrix_L,new_values=L_full_matrix,start_row=1,start_col=1,&
                               n_rows=RI_dimen, n_cols=RI_dimen, error=error)

    info_chol=0
    CALL cp_fm_cholesky_decompose(matrix=fm_matrix_L, n=RI_dimen, info_out=info_chol, error=error)
    CPPostcondition(info_chol==0,cp_failure_level,routineP,error,failure)

    ! triangual invert
    CALL cp_fm_triangular_invert(matrix_a=fm_matrix_L,uplo_tr='U',error=error)

    ! replicate L matrix to each proc
    L_full_matrix=0.0_dp
    CALL cp_fm_get_submatrix(fm_matrix_L, L_full_matrix, 1, 1, RI_dimen, RI_dimen, .FALSE., error)
    CALL cp_fm_release(fm_matrix_L,error=error)

    ! clean lower part
    DO iiB=1, RI_dimen
      L_full_matrix(iiB+1:RI_dimen,iiB)=0.0_dp
    END DO

    ALLOCATE(list_kl%elements(natom**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeffs_kind_max0=MAXVAL(screen_coeffs_kind(:,:)%x(2))
    ALLOCATE(set_list_kl((max_set*natom)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !! precalculate maximum density matrix elements in blocks
    actual_x_data%pmax_block = 0.0_dp
    shm_pmax_block => actual_x_data%pmax_block

    shm_atomic_pair_list => actual_x_data%atomic_pair_list
   
    iatom_start=1
    iatom_end=natom 
    jatom_start=1
    jatom_end=natom 
    katom_start=1
    katom_end=natom 
    latom_start=1
    latom_end=natom 

    CALL build_pair_list_mp2(natom, list_kl,set_list_kl, katom_start, katom_end, &
                     latom_start, latom_end, &
                     kind_of,basis_parameter,particle_set, &
                     do_periodic,screen_coeffs_set,screen_coeffs_kind, &
                     coeffs_kind_max0,log10_eps_schwarz, cell, 0.D+00, &
                     shm_atomic_pair_list)

     virtual=dimen-occupied

     ALLOCATE(Lai(RI_dimen,virtual,occupied),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     Lai=0.0_dp

     DO iatom=1, natom
       jatom=iatom

       ikind=kind_of(iatom)
       jkind=kind_of(jatom)
       ra=particle_set(iatom)%r(:)
       rb=particle_set(jatom)%r(:)
       !rab2=(ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
       rab2=0.0_dp

       la_max => RI_basis_parameter(ikind)%lmax
       la_min => RI_basis_parameter(ikind)%lmin
       npgfa => RI_basis_parameter(ikind)%npgf
       nseta = RI_basis_parameter(ikind)%nset
       zeta => RI_basis_parameter(ikind)%zet
       nsgfa => RI_basis_parameter(ikind)%nsgf
       sphi_a_ext => RI_basis_parameter(ikind)%sphi_ext(:,:,:,:)
       nsgfl_a => RI_basis_parameter(ikind)%nsgfl
       sphi_a_u1=UBOUND(sphi_a_ext,1)
       sphi_a_u2=UBOUND(sphi_a_ext,2)
       sphi_a_u3=UBOUND(sphi_a_ext,3)

       lb_max => basis_S0(jkind)%lmax
       lb_min => basis_S0(jkind)%lmin
       npgfb => basis_S0(jkind)%npgf
       nsetb = basis_S0(jkind)%nset
       zetb => basis_S0(jkind)%zet
       nsgfb => basis_S0(jkind)%nsgf
       sphi_b_ext => basis_S0(jkind)%sphi_ext(:,:,:,:)
       nsgfl_b => basis_S0(jkind)%nsgfl
       sphi_b_u1=UBOUND(sphi_b_ext,1)
       sphi_b_u2=UBOUND(sphi_b_ext,2)
       sphi_b_u3=UBOUND(sphi_b_ext,3)

       jset=1
       DO iset=1, nseta

         counter_L_blocks=counter_L_blocks+1
         IF(MOD(counter_L_blocks,para_env%num_pe)/=para_env%mepos) CYCLE

         ncob = npgfb(jset)*ncoset(lb_max(jset))
         !max_val1 = screen_coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
         !           screen_coeffs_set(jset,iset,jkind,ikind)%x(2)
         sphi_a_ext_set => sphi_a_ext(:,:,:,iset)
         sphi_b_ext_set => sphi_b_ext(:,:,:,jset)

         L_B_i_start=RI_index_table(iatom,iset)
         L_B_i_end=RI_index_table(iatom,iset)+nsgfa(iset)-1

         ALLOCATE(BI1(dimen,dimen,nsgfa(iset)),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         BI1=0.0_dp

         DO i_list_kl=1,list_kl%n_element

           katom=list_kl%elements(i_list_kl)%pair(1)
           latom=list_kl%elements(i_list_kl)%pair(2)

           i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
           i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
           kkind=list_kl%elements(i_list_kl)%kind_pair(1)
           lkind=list_kl%elements(i_list_kl)%kind_pair(2)
           rc=list_kl%elements(i_list_kl)%r1
           rd=list_kl%elements(i_list_kl)%r2
           rcd2=list_kl%elements(i_list_kl)%dist2

           pmax_atom = 0.0_dp

           !screen_kind_ij = screen_coeffs_kind(jkind,ikind)%x(1)*rab2+&
           !                 screen_coeffs_kind(jkind,ikind)%x(2)
           !screen_kind_kl = screen_coeffs_kind(lkind,kkind)%x(1)*rcd2+&
           !                 screen_coeffs_kind(lkind,kkind)%x(2)

           !!!!!! Change the loop order
           !IF( max_val1 + screen_kind_kl + pmax_atom < log10_eps_schwarz) CYCLE
           !!!!!!
           !IF( screen_kind_ij + screen_kind_kl + pmax_atom < log10_eps_schwarz ) CYCLE

           lc_max => basis_parameter(kkind)%lmax
           lc_min => basis_parameter(kkind)%lmin
           npgfc => basis_parameter(kkind)%npgf
           zetc => basis_parameter(kkind)%zet
           nsgfc => basis_parameter(kkind)%nsgf
           sphi_c_ext => basis_parameter(kkind)%sphi_ext(:,:,:,:)
           nsgfl_c => basis_parameter(kkind)%nsgfl
           sphi_c_u1=UBOUND(sphi_c_ext,1)
           sphi_c_u2=UBOUND(sphi_c_ext,2)
           sphi_c_u3=UBOUND(sphi_c_ext,3)

           ld_max => basis_parameter(lkind)%lmax
           ld_min => basis_parameter(lkind)%lmin
           npgfd => basis_parameter(lkind)%npgf
           zetd => basis_parameter(lkind)%zet
           nsgfd => basis_parameter(lkind)%nsgf
           sphi_d_ext => basis_parameter(lkind)%sphi_ext(:,:,:,:)
           nsgfl_d => basis_parameter(lkind)%nsgfl
           sphi_d_u1=UBOUND(sphi_d_ext,1)
           sphi_d_u2=UBOUND(sphi_d_ext,2)
           sphi_d_u3=UBOUND(sphi_d_ext,3)

           DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
             kset=set_list_kl(i_set_list_kl)%pair(1)
             lset=set_list_kl(i_set_list_kl)%pair(2)

             IF(katom==latom .AND. lset<kset) CYCLE

             orb_k_start=mp2_biel%index_table(katom,kset)
             orb_k_end=orb_k_start+nsgfc(kset)-1
             orb_l_start=mp2_biel%index_table(latom,lset)
             orb_l_end=orb_l_start+nsgfd(lset)-1

             !max_val2_set = (screen_coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
             !                screen_coeffs_set(lset,kset,lkind,kkind)%x(2) )
             !max_val2 = max_val1 + max_val2_set

             !!! Near field screening
             !IF(max_val2 + pmax_atom <log10_eps_schwarz) CYCLE
             !! get max_vals if we screen on initial density
             pmax_entry = 0.0_dp

             sphi_c_ext_set => sphi_c_ext(:,:,:,kset)
             sphi_d_ext_set => sphi_d_ext(:,:,:,lset)

             log10_pmax = pmax_entry
             log10_eps_schwarz=log_zero
             !max_val2 = max_val2 + log10_pmax
             !!IF(max_val2<log10_eps_schwarz) CYCLE
             !pmax_entry = EXP(log10_pmax*ln_10)

             IF(ALLOCATED(MNRS)) DEALLOCATE(MNRS)
             ALLOCATE(MNRS(nsgfd(lset),nsgfc(kset),nsgfa(iset)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

             MNRS=0.0_dp
 
               !max_contraction_val =  max_contraction(iset,iatom) * &
               !                       max_contraction(jset,jatom) * &
               !                       max_contraction(kset,katom) * &
               !                       max_contraction(lset,latom) * pmax_entry
               !tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
               !tmp_R_2 => radii_pgf(:,:,lset,kset,lkind,kkind)
               !tmp_screen_pgf1 => screen_coeffs_pgf(:,:,jset,iset,jkind,ikind)
               !tmp_screen_pgf2 => screen_coeffs_pgf(:,:,lset,kset,lkind,kkind)

               max_contraction_val = max_contraction(kset,katom) * max_contraction(lset,latom)
               ! tmp_R_1 => radii_pgf(:,:,kset,kset,kkind,kkind)
               ! tmp_R_2 => radii_pgf(:,:,kset,kset,kkind,kkind)
               ! tmp_screen_pgf1 => screen_coeffs_pgf(:,:,kset,kset,kkind,kkind)
               ! tmp_screen_pgf2 => screen_coeffs_pgf(:,:,kset,kset,kkind,kkind)

               CALL coulomb4(private_lib, ra, rb, rc, rd, npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset), &
                             la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                             lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                             nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset),&
                             sphi_a_u1,sphi_a_u2,sphi_a_u3,&
                             sphi_b_u1,sphi_b_u2,sphi_b_u3,&
                             sphi_c_u1,sphi_c_u2,sphi_c_u3,&
                             sphi_d_u1,sphi_d_u2,sphi_d_u3,&
                             zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset),&
                             zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                             primitive_integrals,&
                             mp2_potential_parameter, &
                             actual_x_data%neighbor_cells, screen_coeffs_set(kset,kset,kkind,kkind)%x,&
                             screen_coeffs_set(kset,kset,kkind,kkind)%x, eps_schwarz, &
                             max_contraction_val, cartesian_estimate,  cell, neris_tmp,&
                             log10_pmax, log10_eps_schwarz, &
                             tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2,&
                             pgf_list_ij, pgf_list_kl, pgf_product_list, &
                             nsgfl_a(:,iset), nsgfl_b(:,jset),&
                             nsgfl_c(:,kset), nsgfl_d(:,lset),&
                             sphi_a_ext_set,&
                             sphi_b_ext_set,&
                             sphi_c_ext_set,&
                             sphi_d_ext_set,&
                             ee_work,ee_work2,ee_buffer1,ee_buffer2, ee_primitives_tmp,&
                             nimages,do_periodic, p_work)

             nints = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
             neris_total = neris_total + nints
             nprim_ints = nprim_ints + neris_tmp
             IF(cartesian_estimate == 0.0_dp) cartesian_estimate = TINY(cartesian_estimate)
             estimate_to_store_int = EXPONENT(cartesian_estimate)
             estimate_to_store_int = MAX(estimate_to_store_int,-15_int_8)
             cartesian_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)

 
             !IF(cartesian_estimate<eps_schwarz) CYCLE
             
             ! WRITE(1111,*) nsgfd(lset),nsgfc(kset),nsgfa(iset)  
 
             primitive_counter=0
             DO llB=1, nsgfd(lset)
               DO kkB=1, nsgfc(kset)
                 DO jjB=1, nsgfb(jset)
                   DO iiB=1, nsgfa(iset)
                     primitive_counter=primitive_counter+1
                     MNRS(llB,kkB,iiB)=primitive_integrals(primitive_counter)
                     ! WRITE(1111,*) orb_l_start+llB-1,orb_k_start+kkB-1,L_B_i_start+iiB-1,MNRS(llB,kkB,iiB)
                   END DO
                 END DO
               END DO
             END DO

             DO iiB=1, nsgfa(iset)
               BI1(orb_l_start:orb_l_end,orb_k_start:orb_k_end,iiB)=MNRS(:,:,iiB)
               BI1(orb_k_start:orb_k_end,orb_l_start:orb_l_end,iiB)=TRANSPOSE(MNRS(:,:,iiB))
             END DO

           END DO ! i_set_list_kl
         END DO ! i_list_kl

         DO iiB=1, nsgfa(iset)
           BI1(1:virtual,1:occupied,iiB)=MATMUL(TRANSPOSE(C(1:dimen,occupied+1:dimen)),&
                                              MATMUL(BI1(1:dimen,1:dimen,iiB),C(1:dimen,1:occupied)))
           Lai(L_B_i_start+iiB-1,1:virtual,1:occupied)=BI1(1:virtual,1:occupied,iiB)
         END DO

         DEALLOCATE(BI1)

      END DO
    END DO

    CALL mp_sum(Lai,para_env%group)

    DO iiB=1, occupied
      IF (MOD(iiB,para_env%num_pe)==para_env%mepos) THEN
        Lai(1:RI_dimen,1:virtual,iiB)=MATMUL(TRANSPOSE(L_full_matrix),Lai(1:RI_dimen,1:virtual,iiB))
      ELSE
        Lai(:,:,iiB)=0.0_dp
      END IF
    END DO

    CALL mp_sum(Lai,para_env%group)

    ! DO iiB=1, occupied
    !   DO jjB=1, virtual
    !     DO kkB=1, RI_dimen
    !       WRITE(3011,*) kkB,jjB,iiB,Lai(kkB,jjB,iiB)
    !     END DO
    !   END DO
    ! END DO

    DEALLOCATE(set_list_kl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,max_pgf**2
      DEALLOCATE(pgf_list_ij(i)%image_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(pgf_list_kl(i)%image_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO

    DEALLOCATE(pgf_list_ij, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pgf_list_kl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pgf_product_list, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(max_contraction, kind_of, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(ee_work, ee_work2, ee_buffer1, ee_buffer2, ee_primitives_tmp, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(nimages, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(mp2_env%potential_parameter%potential_type == do_mp2_potential_TShPSC) THEN
      init_TShPSC_lmax = -1
      CALL free_C0()
    END IF

    CALL timestop(handle)

  END SUBROUTINE calc_lai_libint

END MODULE mp2_ri_libint
