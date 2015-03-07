!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of STM image as post processing of an electronic
!>     structure calculation,
!> \par History
!>      Started as a copy from the code in qs_scf_post
!> \author Joost VandeVondele 7.2008, MI 02.2009
! *****************************************************************************
MODULE stm_images
  USE cp_array_utils_r,                ONLY: cp_1d_r_p_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'stm_images'
  PUBLIC :: th_stm_image


CONTAINS
! *****************************************************************************
!> \brief Driver for the calculation of STM image, as post processing of a
!>        ground-state electronic structure calculation.
!> \param qs_env ...
!> \param stm_section ...
!> \param particles ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param error ...
!> \param
!> \par History
!>      02.2009 Created [MI]
!> \author MI
!> \note
!>   The Tersoff-Hamman
!>        approximation is applied, occupied and a sufficient number of
!>        unoccupied eigenstates are needed (depending on the given Bias potential)
!>        and should be computed in advance. Unoccupied states are calculated
!>        before enetering this module when NLUMO =/ 0
! *****************************************************************************

  SUBROUTINE th_stm_image(qs_env, stm_section, particles, unoccupied_orbs, &
             unoccupied_evals, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: stm_section
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_evals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'th_stm_image', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, irep, ispin, istat, &
                                                n_rep, ndim, nmo, nspin, &
                                                output_unit
    INTEGER, DIMENSION(:), POINTER           :: nadd_unocc, stm_th_torb
    LOGICAL                                  :: append_cube, failure, &
                                                use_ref_energy
    REAL(KIND=dp)                            :: efermi, ref_energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, mo_occ, &
                                                stm_biases
    TYPE(cp_1d_r_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: evals, occupation
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(cp_dbcsr_type), POINTER             :: stm_density_ao
    TYPE(cp_fm_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: mo_arrays
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

     NULLIFY (ks_env, mos, rho, rho_ao, pw_env, stm_th_torb, fm_struct_tmp)
     NULLIFY (auxbas_pw_pool, pw_pools, stm_density_ao, mo_coeff)

     CALL get_qs_env(qs_env,&
                     ks_env=ks_env,&
                     mos=mos,&
                     rho=rho,&
                     pw_env=pw_env,&
                     error=error)

     CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)

     CALL section_vals_val_get(stm_section,"APPEND", l_val=append_cube,error=error)
     CALL section_vals_val_get(stm_section,"BIAS", r_vals=stm_biases,error=error)
     CALL section_vals_val_get(stm_section,"REF_ENERGY",r_val=ref_energy,explicit=use_ref_energy,error=error)
     CALL section_vals_val_get(stm_section,"TH_TORB", n_rep_val=n_rep, error=error)
     IF(n_rep==0) THEN
        ALLOCATE (stm_th_torb(1), STAT=istat)
        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
        stm_th_torb(1) = 0
     ELSE
        ALLOCATE (stm_th_torb(n_rep), STAT=istat)
        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       DO irep = 1,n_rep
          CALL section_vals_val_get(stm_section,"TH_TORB",&
               i_rep_val=irep,i_val=stm_th_torb(irep), error=error)
       END DO
     END IF

     ALLOCATE(stm_density_ao)
     CALL cp_dbcsr_init(stm_density_ao, error=error)
     CALL cp_dbcsr_copy(stm_density_ao,rho_ao(1)%matrix,&
          name="stm_density_ao",error=error)

     CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                     pw_pools=pw_pools,error=error)
     CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
     CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)

     nspin = SIZE(mos,1)
     ALLOCATE(nadd_unocc(nspin),STAT=istat)
     CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
     nadd_unocc=0
     IF(ASSOCIATED(unoccupied_orbs)) THEN
       DO ispin =1,nspin
         nadd_unocc(ispin)=SIZE(unoccupied_evals(ispin)%array)
       END DO
     END IF

     ALLOCATE(mo_arrays(nspin), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(evals(nspin),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(occupation(nspin),STAT=istat)
     DO ispin=1,nspin
          IF(nadd_unocc(ispin)==0) THEN
            CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                 eigenvalues=mo_eigenvalues, nmo=nmo, mu=efermi, occupation_numbers=mo_occ)
            mo_arrays(ispin)%matrix=>mo_coeff
            evals(ispin)%array=>mo_eigenvalues
            occupation(ispin)%array => mo_occ
          ELSE
            CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                 eigenvalues=mo_eigenvalues, nmo=nmo, mu=efermi, occupation_numbers=mo_occ)
             ndim = nmo+nadd_unocc(ispin)
             ALLOCATE(evals(ispin)%array(ndim),STAT=istat)
             evals(ispin)%array(1:nmo)=mo_eigenvalues(1:nmo)
             evals(ispin)%array(1+nmo:ndim)=unoccupied_evals(ispin)%array(1:nadd_unocc(ispin))
             ALLOCATE(occupation(ispin)%array(ndim),STAT=istat)
             occupation(ispin)%array(1:nmo)=mo_occ(1:nmo)
             occupation(ispin)%array(1+nmo:ndim)=0.0_dp
             CALL cp_fm_struct_create(fm_struct_tmp, ncol_global=ndim, &
                          template_fmstruct=mo_coeff%matrix_struct, error=error)
             CALL cp_fm_create(mo_arrays(ispin)%matrix,fm_struct_tmp, name="mo_arrays",error=error)
             CALL cp_fm_struct_release(fm_struct_tmp,error=error)
             CALL cp_fm_to_fm(mo_coeff, mo_arrays(ispin)%matrix, nmo, 1, 1)
             CALL cp_fm_to_fm(unoccupied_orbs(ispin)%matrix, mo_arrays(ispin)%matrix, &
                    nadd_unocc(ispin), 1, nmo+1)
          END IF
     ENDDO
     IF(use_ref_energy) efermi = ref_energy
     
     CALL stm_cubes(ks_env, stm_section, stm_density_ao, wf_r, wf_g, mo_arrays, evals, &
                    occupation, efermi, stm_biases, stm_th_torb,  particles, &
                    output_unit, append_cube,  error=error)
     DO ispin = 1,nspin
       IF(nadd_unocc(ispin)>0) THEN
            DEALLOCATE(evals(ispin)%array,STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
            DEALLOCATE(occupation(ispin)%array,STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
            CALL cp_fm_release( mo_arrays(ispin)%matrix, error=error)
       END IF
     END DO
     DEALLOCATE(mo_arrays,STAT=istat)
        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE(evals,STAT=istat)
        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE(occupation,STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL cp_dbcsr_deallocate_matrix(stm_density_ao,error=error)
     CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
     CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw, error=error)

     DEALLOCATE (stm_th_torb, STAT=istat)
     CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE (nadd_unocc, STAT=istat)
     CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE th_stm_image

! *****************************************************************************
!> \brief computes a simple approximation to the tunneling current for STM
!> \param ks_env ...
!> \param stm_section ...
!> \param stm_density_ao ...
!> \param wf_r ...
!> \param wf_g ...
!> \param mo_arrays ...
!> \param evals ...
!> \param occupation ...
!> \param efermi ...
!> \param stm_biases ...
!> \param stm_th_torb ...
!> \param particles ...
!> \param output_unit ...
!> \param append_cube ...
!> \param error ...
!> \param
!> \par History
!>      7.2008 Created [Joost VandeVondele]
!>       07.2009 modified MI
!> \author Joost VandeVondele
!> \note
!>      requires the MOs that are passed to be eigenstates, and energy ordered
! *****************************************************************************
  SUBROUTINE stm_cubes(ks_env, stm_section, stm_density_ao, wf_r, wf_g, mo_arrays, evals, &
                       occupation, efermi, stm_biases,  stm_th_torb, particles, &
                       output_unit, append_cube,  error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(section_vals_type), POINTER         :: stm_section
    TYPE(cp_dbcsr_type), POINTER             :: stm_density_ao
    TYPE(pw_p_type)                          :: wf_r, wf_g
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(IN)                             :: mo_arrays
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      INTENT(IN)                             :: evals, occupation
    REAL(KIND=dp)                            :: efermi
    REAL(KIND=dp), DIMENSION(:), POINTER     :: stm_biases
    INTEGER, DIMENSION(:), POINTER           :: stm_th_torb
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: append_cube
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), DIMENSION(0:9), PARAMETER :: torb_string = (/"  s"," px"&
      ," py", " pz", "dxy", "dyz", "dzx", "dx2", "dy2", "dz2"/)
    CHARACTER(len=*), PARAMETER :: routineN = 'stm_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: my_pos, oname, title
    INTEGER                                  :: handle, i, ibias, imo, iorb, &
                                                ispin, istat, istates, nmo, &
                                                nspin, nstates(2), &
                                                state_start(2), unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, total_rho
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: occ_tot
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_v, matrix_vf
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(fm_struct_tmp)

    nspin = SIZE(mo_arrays)

    IF (output_unit>0) WRITE(output_unit,'(T2,A)')       ""
    IF (output_unit>0) WRITE(output_unit,'(T2,A,F12.6, A)') "STM : Reference energy ", efermi , " a.u. "
    DO ibias=1,SIZE(stm_biases)

       IF (output_unit>0) WRITE(output_unit,'(T2,A)')       ""
       IF (output_unit>0) WRITE(output_unit,'(T2,A,F16.6)') &
                            "Preparing for STM image at bias [a.u.] ",stm_biases(ibias)

       istates=0
       nstates = 0
       state_start=0
       DO ispin=1,nspin
          IF (stm_biases(ibias)<0.0_dp) THEN
             nmo = SIZE(evals(ispin)%array)
             DO imo =  1,nmo
               IF(evals(ispin)%array(imo) >(efermi+stm_biases(ibias)) .AND. &
                  evals(ispin)%array(imo)<= efermi) THEN
                  IF(nstates(ispin)==0) state_start(ispin) = imo
                  nstates(ispin) = nstates(ispin)+1
               END IF
             END DO
             IF ((output_unit>0).AND.evals(ispin)%array(1)> efermi+stm_biases(ibias)) &
                 WRITE(output_unit,'(T4,A)') "Warning: EFermi+bias below lowest computed occupied MO"
          ELSE
             nmo = SIZE(evals(ispin)%array)
             DO imo =  1,nmo
               IF(evals(ispin)%array(imo) <=(efermi+stm_biases(ibias)) .AND. &
                  evals(ispin)%array(imo)> efermi) THEN
                  IF(nstates(ispin)==0) state_start(ispin) = imo
                  nstates(ispin) = nstates(ispin)+1
               END IF
             END DO
             IF ((output_unit>0) .AND.  evals(ispin)%array(nmo) < efermi+stm_biases(ibias)) &
                 WRITE(output_unit,'(T4,A)') "Warning: E-Fermi+bias above highest computed unoccupied MO"
          ENDIF
          istates=istates+nstates(ispin)
       ENDDO
       IF ((output_unit>0)) WRITE(output_unit,'(T4,A,I0,A)') "Using a total of ",istates," states"
       IF(istates==0) CYCLE

       CALL cp_fm_struct_create(fm_struct_tmp, ncol_global=istates, &
                       template_fmstruct=mo_arrays(1)%matrix%matrix_struct, error=error)
       CALL cp_fm_create(matrix_v,fm_struct_tmp, name="matrix_v",error=error)
       CALL cp_fm_create(matrix_vf,fm_struct_tmp, name="matrix_vf",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       ALLOCATE(occ_tot(istates),STAT=istat)

       ! we sum both alpha and beta electrons together for this density of states
       istates=0
       alpha=1.0_dp
       IF(nspin==1) alpha=2.0_dp
       DO ispin=1,nspin
          CALL cp_fm_to_fm(mo_arrays(ispin)%matrix,matrix_v,nstates(ispin),state_start(ispin),istates+1)
          CALL cp_fm_to_fm(mo_arrays(ispin)%matrix,matrix_vf,nstates(ispin),state_start(ispin),istates+1)
          IF (stm_biases(ibias)<0.0_dp) THEN
            occ_tot(istates+1:istates+nstates(ispin)) = &
                 occupation(ispin)%array(state_start(ispin):state_start(ispin)-1+nstates(ispin))
          ELSE
            occ_tot(istates+1:istates+nstates(ispin)) = alpha - &
                 occupation(ispin)%array(state_start(ispin):state_start(ispin)-1+nstates(ispin))
          END IF
          istates=istates+nstates(ispin)
       ENDDO

       CALL cp_fm_column_scale(matrix_vf,occ_tot(1:istates))
       alpha=1.0_dp

       CALL cp_dbcsr_set(stm_density_ao,0.0_dp,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(stm_density_ao,matrix_v=matrix_v, matrix_g=matrix_vf,ncol=istates,&
                               alpha=alpha,error=error)

       DO i = 1,SIZE(stm_th_torb)
         iorb=stm_th_torb(i)
         CALL calculate_rho_elec(matrix_p=stm_density_ao,&
                 rho=wf_r,rho_gspace=wf_g, total_rho=total_rho,&
                 ks_env=ks_env,  der_type=iorb, error=error)

         oname = torb_string(iorb)
!         fname = "STM_"//TRIM(torb_string(iorb))
         WRITE(filename,'(a4,I2.2,a1,I5.5)')"STM_d",iorb,"_",ibias
         my_pos = "REWIND"
         IF (append_cube) THEN
           my_pos = "APPEND"
         END IF

         unit_nr=cp_print_key_unit_nr(logger,stm_section,extension=".cube",&
              middle_name=TRIM(filename),file_position=my_pos, file_action="WRITE",&
              log_filename=.FALSE.,error=error)
         WRITE(title,'(A,I0,A,I0,A,F16.8)') "STM cube ",ibias," wfn deriv. ",iorb," at bias ",stm_biases(ibias)
         CALL cp_pw_to_cube(wf_r%pw,unit_nr,title,particles=particles,&
             stride=section_get_ivals(stm_section,"STRIDE",error=error),zero_tails=.TRUE., error=error)

         CALL cp_print_key_finished_output(unit_nr,logger,stm_section,error=error)
       END DO

       CALL cp_fm_release(matrix_v, error=error)
       CALL cp_fm_release(matrix_vf, error=error)
       DEALLOCATE(occ_tot,STAT=istat)

    ENDDO

    CALL timestop(handle)

  END SUBROUTINE stm_cubes

END MODULE stm_images
