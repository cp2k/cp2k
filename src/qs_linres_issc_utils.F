!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Chemical shift calculation by dfpt
!>      Initialization of the issc_env, creation of the special neighbor lists
!>      Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!>      Write output
!>      Deallocate everything
!> \note
!>      The psi0 should be localized
!>      the Sebastiani method works within the assumption that the orbitals are
!>      completely contained in the simulation box
! *****************************************************************************
MODULE qs_linres_issc_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, &
       convert_offsets_to_sizes, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_set, dbcsr_distribution_obj, &
       dbcsr_type_antisymmetric, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                             cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: diamat_all
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE physcon,                         ONLY: a_fine,&
                                             e_mass,&
                                             hertz,&
                                             p_mass
  USE qs_elec_field,                   ONLY: build_efg_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_fermi_contact,                ONLY: build_fermi_contact_matrix
  USE qs_linres_methods,               ONLY: linres_solver
  USE qs_linres_types,                 ONLY: get_issc_env,&
                                             issc_env_create,&
                                             issc_env_type,&
                                             linres_control_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_spin_orbit,                   ONLY: build_pso_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: issc_env_cleanup, issc_env_init, issc_response, issc_issc, issc_print

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_issc_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the issc environment
! *****************************************************************************
  SUBROUTINE issc_response(issc_env,p_env,qs_env,error)
    !
    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_response', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ijdir, ispin, &
                                                istat, jdir, nao, nmo, &
                                                nspins, output_unit
    LOGICAL                                  :: do_dso, do_fc, do_pso, do_sd, &
                                                failure, should_stop
    REAL(dp)                                 :: chk, fro
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fc_psi0, h1_psi0, psi0_order, &
                                                psi1, psi1_fc
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dso_psi0, efg_psi0, psi1_dso, &
                                                psi1_efg, psi1_pso, pso_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: issc_section, lr_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    !
    NULLIFY(dft_control, linres_control, lr_section, issc_section)
    NULLIFY(logger, mpools, psi1,h1_psi0, mo_coeff,para_env)
    NULLIFY(tmp_fm_struct, psi1_fc, psi1_efg, psi1_pso, pso_psi0, fc_psi0, efg_psi0, psi0_order)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    issc_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%SPINSPIN",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)
    IF(output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "*** Self consistent optimization of the response wavefunctions ***"
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    mpools=mpools,&
                    linres_control=linres_control,&
                    mos=mos,&
                    para_env=para_env,&
                    error=error)

    nspins = dft_control%nspins

    CALL get_issc_env(issc_env=issc_env, &
                      !list_cubes=list_cubes, &
                      psi1_efg=psi1_efg, &
                      psi1_pso=psi1_pso, &
                      psi1_dso=psi1_dso, &
                      psi1_fc=psi1_fc, &
                      efg_psi0=efg_psi0, &
                      pso_psi0=pso_psi0, &
                      dso_psi0=dso_psi0, &
                      fc_psi0=fc_psi0, &
                      do_fc=do_fc, &
                      do_sd=do_sd, &
                      do_pso=do_pso, &
                      do_dso=do_dso, &
                      error=error)
    !
    ! allocate the vectors
    ALLOCATE(psi0_order(nspins))
    ALLOCATE(psi1(nspins),h1_psi0(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       psi0_order(ispin)%matrix => mo_coeff
       CALL cp_fm_get_info(mo_coeff,ncol_global=nmo,nrow_global=nao,error=error)
       NULLIFY(tmp_fm_struct,psi1(ispin)%matrix,h1_psi0(ispin)%matrix)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
            &                   ncol_global=nmo,&
            &                   context=mo_coeff%matrix_struct%context,&
            &                   error=error)
       CALL cp_fm_create(psi1(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(h1_psi0(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    ENDDO
    chk = 0.0_dp
    should_stop =.FALSE.
    !
    ! operator efg
    IF(do_sd) THEN
       ijdir = 0
       DO idir = 1,3
       DO jdir = idir,3
          ijdir = ijdir+1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1_efg(ispin,ijdir)%matrix,0.0_dp,error=error)
          ENDDO
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator efg_"//ACHAR(idir+119)//ACHAR(jdir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(p_psi0(ispin,ijdir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(efg_psi0(ispin,ijdir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          !
          linres_control%lr_triplet = .FALSE.
          linres_control%do_kernel  = .FALSE.
          linres_control%converged  = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,should_stop,error=error)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_efg(ispin,ijdir)%matrix,error=error)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
             chk = chk + fro
          ENDDO
          !
          ! print response functions
          !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
          !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          !   ncubes = SIZE(list_cubes,1)
          !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
          !   DO ispin = 1,nspins
          !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
          !            centers_set(ispin)%array,print_key,'psi1_efg',&
          !            idir=ijdir,ispin=ispin,error=error)
          !   ENDDO ! ispin
          !ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          !
          ! Write the result in the restart file
       ENDDO ! jdir
       ENDDO ! idir
    ENDIF
    !
    ! operator pso
    IF(do_pso) THEN
       DO idir = 1,3
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1_pso(ispin,idir)%matrix,0.0_dp,error=error)
          ENDDO
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator pso_"//ACHAR(idir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(pso_psi0(ispin,idir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          !
          linres_control%lr_triplet = .FALSE. ! we do singlet response
          linres_control%do_kernel  = .FALSE. ! we do uncoupled response
          linres_control%converged  = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,should_stop,error=error)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_pso(ispin,idir)%matrix,error=error)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
             chk = chk + fro
          ENDDO
          !
          ! print response functions
          !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
          !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          !   ncubes = SIZE(list_cubes,1)
          !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
          !   DO ispin = 1,nspins
          !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
          !           centers_set(ispin)%array,print_key,'psi1_pso',&
          !           idir=idir,ispin=ispin,error=error)
          !   ENDDO ! ispin
          !ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          !
          ! Write the result in the restart file
       ENDDO ! idir
    ENDIF
    !
    ! operator fc
    IF(do_fc) THEN
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1_fc(ispin)%matrix,0.0_dp,error=error)
       ENDDO
       IF(output_unit>0) THEN
          WRITE(output_unit,"(T10,A)") "Response to the perturbation operator fc"
       ENDIF
       !
       !Initial guess for psi1
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
          !CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
          !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
       ENDDO
       !
       !DO scf cycle to optimize psi1
       DO ispin = 1,nspins
          CALL cp_fm_to_fm(fc_psi0(ispin)%matrix,h1_psi0(ispin)%matrix,error=error)
       ENDDO
       !
       !
       linres_control%lr_triplet = .TRUE. ! we do triplet response
       linres_control%do_kernel  = .TRUE. ! we do coupled response
       linres_control%converged  = .FALSE.
       CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,should_stop,error=error)
       !
       !
       ! copy the response
       DO ispin=1,nspins
          CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_fc(ispin)%matrix,error=error)
          CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
          chk = chk + fro
       ENDDO
       !
       ! print response functions
       !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
       !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
       !   ncubes = SIZE(list_cubes,1)
       !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
       !   DO ispin = 1,nspins
       !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
       !           centers_set(ispin)%array,print_key,'psi1_pso',&
       !           idir=idir,ispin=ispin,error=error)
       !   ENDDO ! ispin
       !ENDIF ! print response functions
       !
       !
       IF(output_unit>0) THEN
          WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
       ENDIF
       !
       ! Write the result in the restart file
    ENDIF

    !>>>> debugging only
    !
    ! here we have the operator r and compute the polarizability for debugging the kernel only
    IF(do_dso) THEN
       DO idir = 1,3
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1_dso(ispin,idir)%matrix,0.0_dp,error=error)
          ENDDO
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Response to the perturbation operator r_"//ACHAR(idir+119)
          ENDIF
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
             !CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix, psi1(ispin)%matrix,error=error)
             !CALL cp_fm_scale(-1.0_dp,psi1(ispin)%matrix,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(dso_psi0(ispin,idir)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          !
          linres_control%lr_triplet = .FALSE. ! we do singlet response
          linres_control%do_kernel  = .TRUE.  ! we do uncoupled response
          linres_control%converged  = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order,should_stop,error=error)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_dso(ispin,idir)%matrix,error=error)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
             chk = chk + fro
          ENDDO
          !
          ! print response functions
          !IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
          !     &   "PRINT%RESPONSE_FUNCTION_CUBES",error=error),cp_p_file)) THEN
          !   ncubes = SIZE(list_cubes,1)
          !   print_key => section_vals_get_subs_vals(issc_section,"PRINT%RESPONSE_FUNCTION_CUBES",error=error)
          !   DO ispin = 1,nspins
          !      CALL qs_print_cubes(qs_env,psi1(ispin)%matrix,ncubes,list_cubes,&
          !           centers_set(ispin)%array,print_key,'psi1_pso',&
          !           idir=idir,ispin=ispin,error=error)
          !   ENDDO ! ispin
          !ENDIF ! print response functions
          !
          !
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
          ENDIF
          !
          ! Write the result in the restart file
       ENDDO ! idir
    ENDIF
    !<<<< debugging only

    !
    !
    ! print the checksum
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| response: CheckSum =',chk
    ENDIF
    !
    !
    ! clean up
    DO ispin = 1,nspins
       CALL cp_fm_release(psi1(ispin)%matrix,error=error)
       CALL cp_fm_release(h1_psi0(ispin)%matrix,error=error)
    ENDDO
    DEALLOCATE(psi1,h1_psi0,psi0_order,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE issc_response


  SUBROUTINE issc_issc(issc_env,qs_env,iatom,error)

    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_issc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, ixyz, jatom, &
                                                jxyz, natom, nmo, nspins
    LOGICAL                                  :: do_dso, do_fc, do_pso, do_sd, &
                                                failure, gapw
    REAL(dp)                                 :: buf, facdso, facfc, facpso, &
                                                facsd, g, issc_dso, issc_fc, &
                                                issc_pso, issc_sd, maxocc
    REAL(dp), DIMENSION(3)                   :: r_i, r_j
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: issc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dso, matrix_efg, &
                                                matrix_fc, matrix_pso
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fc_psi0, psi1_fc
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_dso, psi1_efg, psi1_pso
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: issc_section

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(cell,dft_control,particle_set,issc,psi1_fc,psi1_efg,psi1_pso)
    NULLIFY(matrix_efg, matrix_fc, matrix_pso, mos, mo_coeff, fc_psi0)

    CALL get_qs_env(qs_env=qs_env,&
         &          cell=cell,&
         &          dft_control=dft_control,&
         &          particle_set=particle_set,&
         &          mos=mos,&
         &          error=error)

    gapw = dft_control%qs_control%gapw
    natom = SIZE(particle_set,1)
    nspins = dft_control%nspins

    CALL get_issc_env(issc_env=issc_env,&
                      matrix_efg=matrix_efg, &
                      matrix_pso=matrix_pso, &
                      matrix_fc=matrix_fc, &
                      matrix_dso=matrix_dso,&
                      psi1_fc=psi1_fc, &
                      psi1_efg=psi1_efg, &
                      psi1_pso=psi1_pso, &
                      psi1_dso=psi1_dso, &
                      fc_psi0=fc_psi0, &
                      issc=issc, &
                      do_fc=do_fc, &
                      do_sd=do_sd, &
                      do_pso=do_pso, &
                      do_dso=do_dso, &
                      error=error)

    g = e_mass / ( 2.0_dp * p_mass )
    facfc  = hertz * g**2 * a_fine**4
    facpso = hertz * g**2 * a_fine**4
    facsd  = hertz * g**2 * a_fine**4
    facdso = hertz * g**2 * a_fine**4

    !
    !
    issc_section => section_vals_get_subs_vals(qs_env%input, &
         & "PROPERTIES%LINRES%SPINSPIN",error=error)
    !
    ! Initialize
    DO ispin = 1,nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,maxocc=maxocc)
       CALL cp_fm_get_info(mo_coeff,ncol_global=nmo,error=error)

       DO jatom = 1,natom
          r_i = particle_set(iatom)%r
          r_j = particle_set(jatom)%r
          r_j = pbc(r_i,r_j,cell)+r_i
          !
          !
          !
          !write(*,*) 'iatom =',iatom,' r_i=',r_i
          !write(*,*) 'jatom =',jatom,' r_j=',r_j
          !
          ! FC term
          !
          IF(do_fc.AND.iatom.NE.jatom) THEN
             !
             ! build the integral for the jatom
             CALL cp_dbcsr_set(matrix_fc(1)%matrix,0.0_dp,error=error)
             CALL build_fermi_contact_matrix(qs_env,matrix_fc,r_j,error)
             CALL cp_dbcsr_sm_fm_multiply(matrix_fc(1)%matrix,mo_coeff,&
                  &                 fc_psi0(ispin)%matrix,ncol=nmo,& ! fc_psi0 a buffer
                  &                 alpha=1.0_dp,error=error)

             CALL cp_fm_trace(fc_psi0(ispin)%matrix,mo_coeff,buf,error=error)
             WRITE(*,*) ' jatom',jatom,'tr(P*fc)=',buf

             CALL cp_fm_trace(fc_psi0(ispin)%matrix,psi1_fc(ispin)%matrix,buf,error=error)
             issc_fc = 2.0_dp * 2.0_dp * maxocc * facfc * buf
             issc(1,1,iatom,jatom,1) = issc(1,1,iatom,jatom,1) + issc_fc
             issc(2,2,iatom,jatom,1) = issc(2,2,iatom,jatom,1) + issc_fc
             issc(3,3,iatom,jatom,1) = issc(3,3,iatom,jatom,1) + issc_fc
          ENDIF
          !
          ! SD term
          !
          IF(do_sd.AND.iatom.NE.jatom) THEN
             !
             ! build the integral for the jatom
             CALL cp_dbcsr_set(matrix_efg(1)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_efg(2)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_efg(3)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_efg(4)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_efg(5)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_efg(6)%matrix,0.0_dp,error=error)
             CALL build_efg_matrix(qs_env,matrix_efg,r_j,error)
             DO ixyz = 1,6
                CALL cp_dbcsr_sm_fm_multiply(matrix_efg(ixyz)%matrix,mo_coeff,&
                     &                 fc_psi0(ispin)%matrix,ncol=nmo,& ! fc_psi0 a buffer
                     &                 alpha=1.0_dp,beta=0.0_dp,error=error)
                CALL cp_fm_trace(fc_psi0(ispin)%matrix,mo_coeff,buf,error=error)
                WRITE(*,*) ' jatom',jatom,ixyz,'tr(P*efg)=',buf
             DO jxyz = 1,6
                CALL cp_fm_trace(fc_psi0(ispin)%matrix,psi1_efg(ispin,jxyz)%matrix,buf,error=error)
                issc_sd = 2.0_dp * maxocc * facsd * buf
                !issc(ixyz,jxyz,iatom,jatom) = issc_sd
                !write(*,*) 'pso_',ixyz,jxyz,' iatom',iatom,'jatom',jatom,issc_pso
             ENDDO
             ENDDO
          ENDIF
          !
          ! PSO term
          !
          IF(do_pso.AND.iatom.NE.jatom) THEN
             !
             ! build the integral for the jatom
             CALL cp_dbcsr_set(matrix_pso(1)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_pso(2)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(matrix_pso(3)%matrix,0.0_dp,error=error)
             CALL build_pso_matrix(qs_env,matrix_pso,r_j,error)
             DO ixyz = 1,3
                CALL cp_dbcsr_sm_fm_multiply(matrix_pso(ixyz)%matrix,mo_coeff,&
                     &                 fc_psi0(ispin)%matrix,ncol=nmo,& ! fc_psi0 a buffer
                     &                 alpha=1.0_dp,beta=0.0_dp,error=error)
             DO jxyz = 1,3
                CALL cp_fm_trace(fc_psi0(ispin)%matrix,psi1_pso(ispin,jxyz)%matrix,buf,error=error)
                issc_pso = -2.0_dp * maxocc * facpso * buf
                issc(ixyz,jxyz,iatom,jatom,3) = issc(ixyz,jxyz,iatom,jatom,3) + issc_pso
             ENDDO
             ENDDO
          ENDIF
          !
          ! DSO term
          !
          !>>>>> for debugging we compute here the polarizability and NOT the DSO term!
          IF(do_dso.AND.iatom.EQ.natom.AND.jatom.EQ.natom) THEN
             !
             ! build the integral for the jatom
             !CALL cp_dbcsr_set(matrix_dso(1)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(2)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(3)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(4)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(5)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(6)%matrix,0.0_dp,error=error)
             !CALL build_dso_matrix(qs_env,matrix_dso,r_i,r_j,error)
             !DO ixyz = 1,6
             !   CALL cp_sm_fm_multiply(matrix_dso(ixyz)%matrix,mo_coeff,&
             !        &                 fc_psi0(ispin)%matrix,ncol=nmo,& ! fc_psi0 a buffer
             !        &                 alpha=1.0_dp,beta=0.0_dp,error=error)
             !   CALL cp_fm_trace(fc_psi0(ispin)%matrix,mo_coeff,buf,error=error)
             !   issc_dso = 2.0_dp * maxocc * facdso * buf
             !   issc(ixyz,jxyz,iatom,jatom,4) = issc_dso
             !ENDDO
             !CALL cp_dbcsr_set(matrix_dso(1)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(2)%matrix,0.0_dp,error=error)
             !CALL cp_dbcsr_set(matrix_dso(3)%matrix,0.0_dp,error=error)
             !CALL rRc_xyz_ao(matrix_dso,qs_env,(/0.0_dp,0.0_dp,0.0_dp/),1,error=error)
             DO ixyz = 1,3
                CALL cp_dbcsr_sm_fm_multiply(matrix_dso(ixyz)%matrix,mo_coeff,&
                     &                 fc_psi0(ispin)%matrix,ncol=nmo,& ! fc_psi0 a buffer
                     &                 alpha=1.0_dp,beta=0.0_dp,error=error)
             DO jxyz = 1,3
                CALL cp_fm_trace(psi1_dso(ispin,jxyz)%matrix,fc_psi0(ispin)%matrix,buf,error=error)
                ! we save the polarizability for a checksum later on !
                issc_dso = 2.0_dp * maxocc * buf
                !WRITE(*,*) ixyz,jxyz,'tr(P_r*r)=',2.0_dp * maxocc * buf
                issc(ixyz,jxyz,iatom,jatom,4) = issc(ixyz,jxyz,iatom,jatom,4) + issc_dso
             ENDDO
             ENDDO

          ENDIF
          !
       ENDDO ! jatom
    ENDDO ! ispin
    !
    !
    ! Finalize
    CALL timestop(handle)
    !
  END SUBROUTINE issc_issc

  SUBROUTINE issc_print(issc_env,qs_env,error)
    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'issc_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol_i, &
                                                element_symbol_j
    CHARACTER(LEN=default_string_length)     :: name_i, name_j, title
    INTEGER                                  :: iatom, jatom, natom, &
                                                output_unit, unit_atoms
    LOGICAL                                  :: do_dso, do_fc, do_pso, do_sd, &
                                                failure, gapw
    REAL(dp) :: eig(3), issc_iso_dso, issc_iso_fc, issc_iso_pso, issc_iso_sd, &
      issc_iso_tot, issc_tmp(3,3)
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: issc
    REAL(dp), EXTERNAL                       :: DDOT
    TYPE(atomic_kind_type), POINTER          :: atom_kind_i, atom_kind_j
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: issc_section

    failure = .FALSE.

    NULLIFY(logger,particle_set,atom_kind_i,atom_kind_j,dft_control)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    issc_section => section_vals_get_subs_vals(qs_env%input, &
         "PROPERTIES%LINRES%SPINSPIN",error=error)

    CALL get_issc_env(issc_env=issc_env, &
         &            issc=issc, &
                      do_fc=do_fc, &
                      do_sd=do_sd, &
                      do_pso=do_pso, &
                      do_dso=do_dso, &
         &            error=error)
    !
    CALL get_qs_env(qs_env=qs_env, &
         &          dft_control=dft_control, &
         &          particle_set=particle_set, &
         &          error=error)

    natom = SIZE(particle_set,1)
    gapw = dft_control%qs_control%gapw

    !
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E14.6)') 'ISSC| CheckSum K =',&
            & SQRT(DDOT(SIZE(issc),issc,1,issc,1))
    ENDIF
    !
    IF(BTEST(cp_print_key_should_output(logger%iter_info,issc_section,&
             "PRINT%K_MATRIX",error=error),cp_p_file)) THEN

       unit_atoms=cp_print_key_unit_nr(logger,issc_section,"PRINT%K_MATRIX",&
                  extension=".data",middle_name="K",log_filename=.FALSE.,&
                  error=error)

       IF(unit_atoms > 0) THEN
          WRITE(unit_atoms,*)
          WRITE(unit_atoms,*)
          WRITE(title,'(A)') "Indirect spin-spin coupling matrix"
          WRITE(unit_atoms,'(T2,A)') title
          DO iatom = 1,natom
             atom_kind_i => particle_set(iatom)%atomic_kind
             CALL get_atomic_kind(atom_kind_i,name=name_i,element_symbol=element_symbol_i)
             DO jatom = 1,natom
                atom_kind_j => particle_set(jatom)%atomic_kind
                CALL get_atomic_kind(atom_kind_j,name=name_j,element_symbol=element_symbol_j)
                !
                IF(iatom.EQ.jatom.AND..NOT.do_dso) CYCLE
                !
                !
                ! FC
                issc_tmp(:,:) = issc(:,:,iatom,jatom,1)
                issc_tmp(:,:) = 0.5_dp * ( issc_tmp(:,:) + TRANSPOSE( issc_tmp(:,:) ) )
                CALL diamat_all(issc_tmp,eig,error=error)
                issc_iso_fc = (eig(1)+ eig(2)+eig(3))/3.0_dp
                !
                ! SD
                issc_tmp(:,:) = issc(:,:,iatom,jatom,2)
                issc_tmp(:,:) = 0.5_dp * ( issc_tmp(:,:) + TRANSPOSE( issc_tmp(:,:) ) )
                CALL diamat_all(issc_tmp,eig,error=error)
                issc_iso_sd = (eig(1)+ eig(2)+eig(3))/3.0_dp
                !
                ! PSO
                issc_tmp(:,:) = issc(:,:,iatom,jatom,3)
                issc_tmp(:,:) = 0.5_dp * ( issc_tmp(:,:) + TRANSPOSE( issc_tmp(:,:) ) )
                CALL diamat_all(issc_tmp,eig,error=error)
                issc_iso_pso = (eig(1)+ eig(2)+eig(3))/3.0_dp
                !
                ! DSO
                issc_tmp(:,:) = issc(:,:,iatom,jatom,4)
                issc_tmp(:,:) = 0.5_dp * ( issc_tmp(:,:) + TRANSPOSE( issc_tmp(:,:) ) )
                CALL diamat_all(issc_tmp,eig,error=error)
                issc_iso_dso = (eig(1)+ eig(2)+eig(3))/3.0_dp
                !
                ! TOT
                issc_iso_tot = issc_iso_fc + issc_iso_sd + issc_iso_dso + issc_iso_pso
                !
                !
                WRITE(unit_atoms,*)
                WRITE(unit_atoms,'(T2,2(A,I5,A,2X,A2))') 'Indirect spin-spin coupling between ',&
                     iatom,TRIM(name_i),element_symbol_i,' and ',&
                     jatom,TRIM(name_j),element_symbol_j
                !
                IF(do_fc ) WRITE(unit_atoms,'(T1,A,f12.4,A)') ' Isotropic FC contribution  = ',issc_iso_fc,' Hz'
                IF(do_sd ) WRITE(unit_atoms,'(T1,A,f12.4,A)') ' Isotropic SD contribution  = ',issc_iso_sd,' Hz'
                IF(do_pso) WRITE(unit_atoms,'(T1,A,f12.4,A)') ' Isotropic PSO contribution = ',issc_iso_pso,' Hz'
                !IF(do_dso) WRITE(unit_atoms,'(T1,A,f12.4,A)') ' Isotropic DSO contribution = ',issc_iso_dso,' Hz'
                IF(do_dso) WRITE(unit_atoms,'(T1,A,f12.4,A)') ' !!! POLARIZABILITY (for the moment) = ',issc_iso_dso,' Hz'
                IF(.NOT.do_dso)WRITE(unit_atoms,'(T1,A,f12.4,A)') ' Isotropic coupling         = ',issc_iso_tot,' Hz'
             ENDDO
          ENDDO
       ENDIF
       CALL cp_print_key_finished_output(unit_atoms,logger,issc_section,&
            &                            "PRINT%K_MATRIX",error=error)
    ENDIF
    !
    !
  END SUBROUTINE issc_print


! *****************************************************************************
!> \brief Initialize the issc environment
! *****************************************************************************
  SUBROUTINE issc_env_init(issc_env,qs_env,error)
    !
    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, idir, ini, ir, &
                                                ispin, istat, m, n, n_rep, &
                                                nao, natom, nspins, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: list, rbs
    LOGICAL                                  :: failure, gapw
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: issc_section, lr_section

!

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(linres_control)
    NULLIFY(logger,issc_section)
    NULLIFY(tmp_fm_struct)
    NULLIFY(particle_set)
    NULLIFY(sab_orb)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         &                             extension=".linresLog",error=error)

    IF(issc_env%ref_count /= 0) THEN
       CALL issc_env_cleanup(issc_env,error=error)
    ENDIF

    IF(output_unit>0) THEN
       WRITE(output_unit,"(/,T20,A,/)") "*** Start indirect spin-spin coupling Calculation ***"
       WRITE(output_unit,"(T10,A,/)") "Inizialization of the ISSC environment"
    ENDIF

    CALL issc_env_create(issc_env,error=error)
    !
    issc_section => section_vals_get_subs_vals(qs_env%input, &
         &          "PROPERTIES%LINRES%SPINSPIN",error=error)
    !CALL section_vals_val_get(nmr_section,"INTERPOLATE_SHIFT",l_val=nmr_env%interpolate_shift,error=error)
    !CALL section_vals_val_get(nmr_section,"SHIFT_GAPW_RADIUS",r_val=nmr_env%shift_gapw_radius,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         &          dft_control=dft_control,&
         &          linres_control=linres_control,&
         &          mos=mos,&
         &          sab_orb=sab_orb,&
         &          particle_set=particle_set,&
         &          dbcsr_dist=dbcsr_dist,&
         &          error=error)
    !
    !
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure.AND.(output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get indirect spin-spin coupling parameters within PBC you need to localize zero order orbitals '
    ENDIF
    !
    !
    ! read terms need to be calculated
    ! FC
    CALL section_vals_val_get(issc_section,"DO_FC",l_val=issc_env%do_fc,error=error)
    ! SD
    CALL section_vals_val_get(issc_section,"DO_SD",l_val=issc_env%do_sd,error=error)
    ! PSO
    CALL section_vals_val_get(issc_section,"DO_PSO",l_val=issc_env%do_pso,error=error)
    ! DSO
    CALL section_vals_val_get(issc_section,"DO_DSO",l_val=issc_env%do_dso,error=error)
    !
    !
    ! read the list of atoms on which the issc need to be calculated
    CALL section_vals_val_get(issc_section,"ISSC_ON_ATOM_LIST",n_rep_val=n_rep,error=error)
    !
    !
    NULLIFY(issc_env%issc_on_atom_list)
    n = 0
    DO ir = 1,n_rep
       NULLIFY(list)
       CALL section_vals_val_get(issc_section,"ISSC_ON_ATOM_LIST",i_rep_val=ir,i_vals=list,error=error)
       IF(ASSOCIATED(list)) THEN
          CALL reallocate(issc_env%issc_on_atom_list,1,n+SIZE(list))
          DO ini = 1,SIZE(list)
             issc_env%issc_on_atom_list(ini+n) = list(ini)
          ENDDO
          n = n + SIZE(list)
       ENDIF
    ENDDO
    !
    IF(.NOT.ASSOCIATED(issc_env%issc_on_atom_list)) THEN
       ALLOCATE(issc_env%issc_on_atom_list(natom),STAT=istat)
       CPPrecondition(istat.EQ.0,cp_warning_level,routineP,error,failure)
       DO iatom = 1,natom
          issc_env%issc_on_atom_list(iatom) = iatom
       ENDDO
    ENDIF
    issc_env%issc_natms = SIZE(issc_env%issc_on_atom_list)
    !
    !
    ! Initialize the issc tensor
    ALLOCATE(issc_env%issc(3,3,issc_env%issc_natms,issc_env%issc_natms,4),&
         &   issc_env%issc_loc(3,3,issc_env%issc_natms,issc_env%issc_natms,4),&
         &   STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    issc_env%issc(:,:,:,:,:) = 0.0_dp
    issc_env%issc_loc(:,:,:,:,:) = 0.0_dp
    !
    ! allocation
    ALLOCATE(issc_env%efg_psi0(nspins,6),issc_env%pso_psi0(nspins,3),issc_env%fc_psi0(nspins),&
         &   issc_env%psi1_efg(nspins,6),issc_env%psi1_pso(nspins,3),issc_env%psi1_fc(nspins),&
         &   issc_env%dso_psi0(nspins,3),issc_env%psi1_dso(nspins,3),&
         &   STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       !mo_coeff => current_env%psi0_order(ispin)%matrix
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,ncol_global=m,nrow_global=nao,error=error)

       NULLIFY(tmp_fm_struct)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
            &                   ncol_global=m,&
            &                   context=mo_coeff%matrix_struct%context,&
            &                   error=error)
       DO idir = 1,6
          NULLIFY(issc_env%psi1_efg(ispin,idir)%matrix,issc_env%efg_psi0(ispin,idir)%matrix)
          CALL cp_fm_create(issc_env%psi1_efg(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(issc_env%efg_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
       ENDDO
       DO idir = 1,3
          NULLIFY(issc_env%psi1_pso(ispin,idir)%matrix,issc_env%pso_psi0(ispin,idir)%matrix,&
               &  issc_env%dso_psi0(ispin,idir)%matrix)
          CALL cp_fm_create(issc_env%psi1_pso(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(issc_env%pso_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(issc_env%psi1_dso(ispin,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(issc_env%dso_psi0(ispin,idir)%matrix,tmp_fm_struct,error=error)
       ENDDO
       NULLIFY(issc_env%psi1_fc(ispin)%matrix,issc_env%fc_psi0(ispin)%matrix)
       CALL cp_fm_create(issc_env%psi1_fc(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(issc_env%fc_psi0(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    ENDDO
    !
    ! prepare for allocation
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf,error=error)
    ALLOCATE (rbs(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    !
    ! efg, pso and fc operators
    CALL cp_dbcsr_allocate_matrix_set(issc_env%matrix_efg,6,error=error)
    ALLOCATE(issc_env%matrix_efg(1)%matrix)
    CALL cp_dbcsr_init(issc_env%matrix_efg(1)%matrix,error=error)
    CALL cp_dbcsr_create(matrix=issc_env%matrix_efg(1)%matrix, &
         name="efg (3xx-rr)/3", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(issc_env%matrix_efg(1)%matrix,sab_orb,error=error)

    ALLOCATE(issc_env%matrix_efg(2)%matrix,&
             issc_env%matrix_efg(3)%matrix,issc_env%matrix_efg(4)%matrix,&
             issc_env%matrix_efg(5)%matrix,issc_env%matrix_efg(6)%matrix)
    CALL cp_dbcsr_init (issc_env%matrix_efg(2)%matrix,error=error)
    CALL cp_dbcsr_init (issc_env%matrix_efg(3)%matrix,error=error)
    CALL cp_dbcsr_init (issc_env%matrix_efg(4)%matrix,error=error)
    CALL cp_dbcsr_init (issc_env%matrix_efg(5)%matrix,error=error)
    CALL cp_dbcsr_init (issc_env%matrix_efg(6)%matrix,error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_efg(2)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'efg xy',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_efg(3)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'efg xz',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_efg(4)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'efg (3yy-rr)/3',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_efg(5)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'efg yz',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_efg(6)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'efg (3zz-rr)/3',error=error)

    CALL cp_dbcsr_set(issc_env%matrix_efg(1)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_efg(2)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_efg(3)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_efg(4)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_efg(5)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_efg(6)%matrix,0.0_dp,error=error)
    !
    ! PSO
    CALL cp_dbcsr_allocate_matrix_set(issc_env%matrix_pso,3,error=error)
    ALLOCATE(issc_env%matrix_pso(1)%matrix)
    CALL cp_dbcsr_init(issc_env%matrix_pso(1)%matrix,error=error)
    CALL cp_dbcsr_create(matrix=issc_env%matrix_pso(1)%matrix, &
         name="pso x", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(issc_env%matrix_pso(1)%matrix,sab_orb,error=error)

    ALLOCATE(issc_env%matrix_pso(2)%matrix,issc_env%matrix_pso(3)%matrix)
    CALL cp_dbcsr_init (issc_env%matrix_pso(2)%matrix,error=error)
    CALL cp_dbcsr_init (issc_env%matrix_pso(3)%matrix,error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_pso(2)%matrix,issc_env%matrix_pso(1)%matrix,&
                   'pso y',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_pso(3)%matrix,issc_env%matrix_pso(1)%matrix,&
                   'pso z',error=error)
    CALL cp_dbcsr_set(issc_env%matrix_pso(1)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_pso(2)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_pso(3)%matrix,0.0_dp,error=error)
    !
    ! DSO
    CALL cp_dbcsr_allocate_matrix_set(issc_env%matrix_dso,3,error=error)
    ALLOCATE(issc_env%matrix_dso(1)%matrix,issc_env%matrix_dso(2)%matrix,issc_env%matrix_dso(3)%matrix)
    CALL cp_dbcsr_init(issc_env%matrix_dso(1)%matrix,error=error)
    CALL cp_dbcsr_init(issc_env%matrix_dso(2)%matrix,error=error)
    CALL cp_dbcsr_init(issc_env%matrix_dso(3)%matrix,error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_dso(1)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'dso x',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_dso(2)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'dso y',error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_dso(3)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'dso z',error=error)
    CALL cp_dbcsr_set(issc_env%matrix_dso(1)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_dso(2)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_set(issc_env%matrix_dso(3)%matrix,0.0_dp,error=error)
    !
    ! FC
    CALL cp_dbcsr_allocate_matrix_set(issc_env%matrix_fc,1,error=error)
    ALLOCATE(issc_env%matrix_fc(1)%matrix)
    CALL cp_dbcsr_init(issc_env%matrix_fc(1)%matrix,error=error)
    CALL cp_dbcsr_copy(issc_env%matrix_fc(1)%matrix,issc_env%matrix_efg(1)%matrix,&
                   'fc',error=error)
    CALL cp_dbcsr_set(issc_env%matrix_fc(1)%matrix,0.0_dp,error=error)

    CALL array_release (row_blk_sizes)
    !
    ! Conversion factors
    IF(output_unit>0) THEN
       WRITE(output_unit,"(T2,A,T60,I4,A)")&
            & "ISSC| spin-spin coupling computed for ",issc_env%issc_natms,' atoms'
    ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE issc_env_init

! *****************************************************************************
!> \brief Deallocate the issc environment
!> \par History
! *****************************************************************************
  SUBROUTINE issc_env_cleanup(issc_env,error)

    TYPE(issc_env_type)                      :: issc_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir, ispin, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    issc_env%ref_count = issc_env%ref_count - 1
    IF(issc_env%ref_count == 0 ) THEN
       IF(ASSOCIATED(issc_env%issc_on_atom_list)) THEN
          DEALLOCATE(issc_env%issc_on_atom_list,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(issc_env%issc)) THEN
          DEALLOCATE(issc_env%issc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(issc_env%issc_loc)) THEN
          DEALLOCATE(issc_env%issc_loc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !efg_psi0
       IF(ASSOCIATED(issc_env%efg_psi0)) THEN
          DO idir = 1,SIZE(issc_env%efg_psi0,2)
             DO ispin = 1,SIZE(issc_env%efg_psi0,1)
                CALL cp_fm_release(issc_env%efg_psi0(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%efg_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !pso_psi0
       IF(ASSOCIATED(issc_env%pso_psi0)) THEN
          DO idir = 1,SIZE(issc_env%pso_psi0,2)
             DO ispin = 1,SIZE(issc_env%pso_psi0,1)
                CALL cp_fm_release(issc_env%pso_psi0(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%pso_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !dso_psi0
       IF(ASSOCIATED(issc_env%dso_psi0)) THEN
          DO idir = 1,SIZE(issc_env%dso_psi0,2)
             DO ispin = 1,SIZE(issc_env%dso_psi0,1)
                CALL cp_fm_release(issc_env%dso_psi0(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%dso_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !fc_psi0
       IF(ASSOCIATED(issc_env%fc_psi0)) THEN
          DO ispin = 1,SIZE(issc_env%fc_psi0,1)
             CALL cp_fm_release(issc_env%fc_psi0(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(issc_env%fc_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_efg
       IF(ASSOCIATED(issc_env%psi1_efg)) THEN
          DO idir = 1,SIZE(issc_env%psi1_efg,2)
             DO ispin = 1,SIZE(issc_env%psi1_efg,1)
                CALL cp_fm_release(issc_env%psi1_efg(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%psi1_efg,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_pso
       IF(ASSOCIATED(issc_env%psi1_pso)) THEN
          DO idir = 1,SIZE(issc_env%psi1_pso,2)
             DO ispin = 1,SIZE(issc_env%psi1_pso,1)
                CALL cp_fm_release(issc_env%psi1_pso(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%psi1_pso,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_dso
       IF(ASSOCIATED(issc_env%psi1_dso)) THEN
          DO idir = 1,SIZE(issc_env%psi1_dso,2)
             DO ispin = 1,SIZE(issc_env%psi1_dso,1)
                CALL cp_fm_release(issc_env%psi1_dso(ispin,idir)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(issc_env%psi1_dso,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       !psi1_fc
       IF(ASSOCIATED(issc_env%psi1_fc)) THEN
          DO ispin = 1,SIZE(issc_env%psi1_fc,1)
             CALL cp_fm_release(issc_env%psi1_fc(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(issc_env%psi1_fc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !
       ! cubes
       !IF(ASSOCIATED(issc_env%list_cubes)) THEN
       !   DEALLOCATE(issc_env%list_cubes,STAT=istat)
       !   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !ENDIF
       !
       !matrix_efg
       IF(ASSOCIATED(issc_env%matrix_efg)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(issc_env%matrix_efg,error=error)
       ENDIF
       !
       !matrix_pso
       IF(ASSOCIATED(issc_env%matrix_pso)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(issc_env%matrix_pso,error=error)
       ENDIF
       !
       !matrix_dso
       IF(ASSOCIATED(issc_env%matrix_dso)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(issc_env%matrix_dso,error=error)
       ENDIF
       !
       !matrix_fc
       IF(ASSOCIATED(issc_env%matrix_fc)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(issc_env%matrix_fc,error=error)
       ENDIF

    ENDIF ! ref count

  END SUBROUTINE issc_env_cleanup

END MODULE qs_linres_issc_utils
