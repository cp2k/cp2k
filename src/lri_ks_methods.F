!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the Kohn-Sham matrix for the LRIGPW 
!>      and xc parts
!> \par History
!>      09.2013 created [Dorothea Golze]
!> \author Dorothea Golze
! *****************************************************************************
MODULE lri_ks_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add_block_node,&
                                             cp_dbcsr_finalize,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE lri_environment_types,           ONLY: lri_environment_type,&
                                             lri_int_type,&
                                             lri_kind_type
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_scale
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_ks_methods'

  PUBLIC :: calculate_lri_ks_matrix

CONTAINS

!*****************************************************************************
!> \brief update of LRIPGW KS matrix
!> \param h_matrix KS matrix, on entry containing the core hamiltonian
!> \param lri_v_int integrals of potential * ri basis set
!> \note including this in lri_environment_methods?
!*****************************************************************************
  SUBROUTINE calculate_lri_ks_matrix(lri_env, lri_v_int, h_matrix, &
                                     atomic_kind_set, para_env,  error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: lri_v_int
    TYPE(cp_dbcsr_type), POINTER             :: h_matrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'calculate_lri_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, col, handle, i, iac, iatom, ikind, ilist, &
      jatom, jkind, jneighbor, natom, nba, nbb, nfa, nfb, nkind, nlist, nn, &
      nneighbor, row, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: h_work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(atom_of_kind, h_block, lrii, nl_iterator, soo_list)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list
 
      nkind = lri_env%lri_ints%nkind

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=natom)
      ALLOCATE(atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
      CALL neighbor_list_iterator_create(nl_iterator,soo_list)

      DO WHILE (neighbor_list_iterate(nl_iterator)==0)
         CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,&
              jatom=jatom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor)

         iac = ikind + nkind*(jkind - 1)
 
         IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

         lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
         
         nfa=lrii%nfa
         nfb=lrii%nfb
         nba = lrii%nba
         nbb = lrii%nbb
         nn=nfa+nfb

         ALLOCATE(h_work(nba,nbb),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         h_work(:,:)=0.0_dp

         IF (iatom <= jatom) THEN
            row = iatom
            col = jatom
            trans = .FALSE.
         ELSE
            row = jatom
            col = iatom
            trans = .TRUE.
         END IF

         NULLIFY(h_block)
         CALL cp_dbcsr_get_block_p(h_matrix,row,col,h_block,found)
         IF (.NOT.ASSOCIATED(h_block)) THEN
              CALL cp_dbcsr_add_block_node (h_matrix, row, col, h_block ,error=error)
         END IF

         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)

         ! sum: Hab(a)+Hab(b)
         ! Hab(a)
         DO i=1,nfa
            h_work(1:nba, 1:nbb) = h_work(1:nba,1:nbb) + &
                                   lri_v_int(ikind)%v_int(atom_a,i)*lrii%dacoef(1:nba,1:nbb,i)
         ENDDO
         ! Hab(b)
         IF(iatom /= jatom) THEN
          DO i=1,nfb
              h_work(1:nba, 1:nbb) = h_work(1:nba,1:nbb) + &
                                     lri_v_int(jkind)%v_int(atom_b,i)*lrii%dacoef(1:nba,1:nbb,nfa+i) 
          ENDDO
         ENDIF

         ! add h_work to core hamiltonian
         IF(trans) THEN
           h_work(1:nba,1:nbb)= h_work(1:nba,1:nbb) + TRANSPOSE(h_block(1:nbb,1:nba))
           h_block(1:nbb,1:nba)=TRANSPOSE(h_work(1:nba,1:nbb))
         ELSE
           h_block(1:nba,1:nbb)=h_block(1:nba,1:nbb)+h_work(1:nba,1:nbb)
         ENDIF
 
         DEALLOCATE(h_work,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO

      CALL cp_dbcsr_finalize(h_matrix, error=error)
 
      CALL neighbor_list_iterator_release(nl_iterator)
      DEALLOCATE(atom_of_kind,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_ks_matrix

! *****************************************************************************
!> routine for testing, to be deleted later
!******************************************************************************
  SUBROUTINE calc_lri_rho_tot_gspace(rho_tot_gspace, qs_env, lri_rho_g,&
                                     skip_nuclear_density, error)
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_tot_gspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: lri_rho_g
    LOGICAL, INTENT(IN), OPTIONAL            :: skip_nuclear_density
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'calc_lri_rho_tot_gspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: my_skip
    TYPE(pw_p_type), POINTER                 :: rho_core

    NULLIFY (rho_core)
    my_skip = .FALSE.
    IF (PRESENT(skip_nuclear_density)) my_skip=skip_nuclear_density
    CALL get_qs_env(qs_env=qs_env,&
                    rho_core=rho_core,&
                    error=error)

    IF (.NOT.my_skip) THEN
       CALL pw_copy(rho_core%pw,rho_tot_gspace%pw, error=error)
       CALL pw_axpy(lri_rho_g(1)%pw,rho_tot_gspace%pw, error=error)
    ELSE
       CALL pw_axpy(lri_rho_g(1)%pw, rho_tot_gspace%pw, error=error)
    END IF
    DO ispin=2, qs_env%dft_control%nspins
       CALL pw_axpy(lri_rho_g(ispin)%pw, rho_tot_gspace%pw, error=error)
    END DO
 
  END SUBROUTINE calc_lri_rho_tot_gspace
! *****************************************************************************
!> only NSPIN=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> routine for testing, to be deleted later
!******************************************************************************
  SUBROUTINE print_lri_density_cube(qs_env, rho_r, lri, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    LOGICAL, INTENT(IN)                      :: lri
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'print_lri_density_cube', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, my_pos_cube
    INTEGER                                  :: output_unit, unit_nr
    LOGICAL                                  :: append_cube
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: dft_section, input

    NULLIFY(dft_section, logger, particles, subsys)

    CALL get_qs_env(qs_env=qs_env, input=input, subsys=subsys,  error=error) 
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    CALL cp_subsys_get(subsys,particles=particles,error=error)
 
    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
       "DFT%PRINT%E_DENSITY_CUBE", error=error),cp_p_file)) THEN

       append_cube = section_get_lval(input,"DFT%PRINT%E_DENSITY_CUBE%APPEND",error=error)
       my_pos_cube="REWIND"
       IF(append_cube) THEN
         my_pos_cube="APPEND"
       END IF 
  
       IF(lri) THEN
        filename = "E_DENSITY_LRI"
       ELSE
        filename = "E_DENSITY_NORMAL"
       ENDIF
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
               extension=".cube",middle_name=TRIM(filename),file_position=my_pos_cube,&
               log_filename=.FALSE.,error=error)
       CALL cp_pw_to_cube(rho_r(1)%pw,unit_nr,"ELECTRONIC DENSITY",&
            particles=particles,&
            stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
            error=error)
       CALL cp_print_key_finished_output(unit_nr,logger,input,&
            "DFT%PRINT%E_DENSITY_CUBE",error=error)
    END IF

  END SUBROUTINE print_lri_density_cube
! *****************************************************************************
!> only NSPIN=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> routine for testing, to be deleted later
!******************************************************************************
  SUBROUTINE print_v_xc_cube(qs_env, v_rspace_new, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: v_rspace_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'print_v_xc_cube', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: my_pos_cube
    INTEGER                                  :: output_unit, unit_nr
    LOGICAL                                  :: append_cube
    REAL(KIND=dp)                            :: udvol
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: aux_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: dft_section, input

    NULLIFY(dft_section, logger, particles, subsys)
    NULLIFY(auxbas_pw_pool, pw_env)

    CALL get_qs_env(qs_env=qs_env, input=input, subsys=subsys,  error=error) 
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    CALL cp_subsys_get(subsys,particles=particles,error=error)
 
    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
         "DFT%PRINT%V_HARTREE_CUBE",error=error),cp_p_file)) THEN

       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,aux_r%pw,&
                               use_data = REALDATA3D,&
                               in_space = REALSPACE, error=error)

       append_cube = section_get_lval(input,"DFT%PRINT%V_HARTREE_CUBE%APPEND",error=error)
       my_pos_cube="REWIND"
       IF(append_cube) THEN
         my_pos_cube="APPEND"
       END IF
       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
       CALL pw_env_get(pw_env, error=error)
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%V_HARTREE_CUBE",&
            extension=".cube",middle_name="v_xc",file_position=my_pos_cube,error=error)
       udvol = 1.0_dp/v_rspace_new%pw%pw_grid%dvol

       CALL pw_copy(v_rspace_new%pw,aux_r%pw, error=error)
       CALL pw_scale(aux_r%pw,udvol,error=error)

       CALL cp_pw_to_cube(aux_r%pw,unit_nr,"Exchange POTENTIAL",particles=particles,&
            stride=section_get_ivals(dft_section,&
            "PRINT%V_HARTREE_CUBE%STRIDE",error=error),&
            error=error)
       CALL cp_print_key_finished_output(unit_nr,logger,input,&
            "DFT%PRINT%V_HARTREE_CUBE",error=error)

       CALL pw_pool_give_back_pw(auxbas_pw_pool,aux_r%pw, error=error)
    ENDIF
  END SUBROUTINE print_v_xc_cube
!
END MODULE lri_ks_methods
