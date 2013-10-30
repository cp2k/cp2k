!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Polarizability calculation by dfpt
!>      Initialization of the polar_env, 
!>      Perturbation Hamiltonian by application of the Berry phase operator to psi0
!>      Write output
!>      Deallocate everything
!> periodic Raman SL February 2013
!> \note
! *****************************************************************************
MODULE qs_linres_polar_utils
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_create,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_alloc_block_from_nbl,&
                                             cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type, &
                                             cp_dbcsr_type
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                             cp_fm_transpose, & 
                                             cp_fm_gemm, &
                                             cp_fm_trace, & 
                                             cp_fm_scale_and_add
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
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             put_results
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_type_antisymmetric,&
                                             dbcsr_type_symmetric,&
                                             dbcsr_type_no_symmetry
  USE dbcsr_util,                      ONLY: convert_offsets_to_sizes
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE mathlib,                         ONLY: diamat_all
  USE memory_utilities,                ONLY: reallocate
  USE physcon,                         ONLY: angstrom
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_methods,               ONLY: linres_solver
  USE qs_linres_types,                 ONLY: get_polar_env,&
                                             polar_env_create,&
                                             polar_env_type,&
                                             linres_control_type
 USE qs_loc_methods,                  ONLY: qs_print_cubes
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: polar_env_cleanup, polar_env_init, polar_response, polar_polar, polar_print

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_polar_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the polar environment
! *****************************************************************************
  SUBROUTINE polar_response(polar_env,p_env,qs_env,error)
   
    INTEGER, DIMENSION(:), POINTER           :: list_cubes
    TYPE(polar_env_type)                      :: polar_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_response', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ijdir, ispin, &
                                                istat, jdir, nao, nmo, &
                                                nspins, output_unit,ncubes,i
    LOGICAL                                  :: do_raman, failure, should_stop
    REAL(dp)                                 :: chk, fro
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_dBerry, dBerry_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: h1_psi0, psi0_order, &
                                                psi1
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(section_vals_type), POINTER         :: polar_section, lr_section,print_key
    REAL(dp), DIMENSION(:, :), POINTER       :: centers


    failure = .FALSE.
    CALL timeset(routineN,handle)
    !
    NULLIFY(dft_control, linres_control, lr_section, polar_section)
    NULLIFY(logger, mpools, psi1,h1_psi0, mo_coeff,para_env,list_cubes)
    NULLIFY(tmp_fm_struct, psi1_dBerry,dBerry_psi0)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    polar_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%POLAR",error=error)

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

    CALL get_polar_env(polar_env=polar_env, &
                      do_raman=do_raman, &
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
                               ncol_global=nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               error=error)
       CALL cp_fm_create(psi1(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(h1_psi0(ispin)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    ENDDO
    chk = 0.0_dp
    !
    IF(do_raman) THEN
       CALL get_polar_env(polar_env=polar_env, &
                      psi1_dBerry=psi1_dBerry, &
                      dBerry_psi0=dBerry_psi0, &
                      error=error)
       DO idir = 1,3
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T10,A)") &
               "Response to the perturbation operator Berry phase_"//ACHAR(idir+119)
          ENDIF
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1_dBerry(idir,ispin)%matrix,0.0_dp,error=error)
          ENDDO
          !
          !Initial guess for psi1
          DO ispin = 1,nspins
             CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
          ENDDO
          !
          !DO scf cycle to optimize psi1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(dBerry_psi0(idir,ispin)%matrix,h1_psi0(ispin)%matrix,error=error)
          ENDDO
          !
          !
          linres_control%lr_triplet = .FALSE. ! we do singlet response
          linres_control%do_kernel  = .TRUE.  ! we do coupled response
          linres_control%converged  = .FALSE.
          CALL linres_solver(p_env,qs_env,psi1,h1_psi0,psi0_order, should_stop,error=error)
          !
          !
          ! copy the response
          DO ispin=1,nspins
             CALL cp_fm_to_fm(psi1(ispin)%matrix,psi1_dBerry(idir,ispin)%matrix,error=error)
             CALL cp_fm_frobenius_norm(psi1(ispin)%matrix,fro,error=error)
             chk = chk + fro
          ENDDO
          !
         IF(output_unit>0) THEN
           WRITE(output_unit,"(T10,A)") "Write the resulting psi1 in restart file... not implemented yet"
         ENDIF
          !
          ! Write the result in the restart file
       ENDDO ! idir
    ENDIF  !raman

    !
    ! print the checksum
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'polar| response: CheckSum =',chk
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
                                     "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE polar_response


  SUBROUTINE polar_polar(polar_env,qs_env,error)

    TYPE(polar_env_type)                      :: polar_env
    INTEGER                                  :: unit_p, output_unit
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_polar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao, j,m, &
                                                istat, natom, nmo, nspins, x,y,z,i,idir
    LOGICAL                                  :: do_raman, failure
    REAL(dp)                                 :: buf, maxocc
    TYPE(cp_logger_type), POINTER            :: logger

    REAL(dp), DIMENSION(:, :), &
      POINTER                                :: polar,polar_tmp
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_dBerry, dBerry_psi0,mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(section_vals_type), POINTER         :: polar_section
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct


    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(cell,dft_control,polar,psi1_dBerry,logger)
    NULLIFY( mos, mo_coeff, dBerry_psi0,mo_derivs)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL get_qs_env(qs_env=qs_env,&
                   cell=cell,&
                   dft_control=dft_control,&
                   mos=mos,&
                   error=error)

    nspins = dft_control%nspins

    CALL get_polar_env(polar_env=polar_env,&
                      do_raman=do_raman, &
                      error=error)
!
    IF(do_raman) THEN
    
        CALL get_polar_env(polar_env=polar_env,&
                      psi1_dBerry=psi1_dBerry, &
                      dBerry_psi0=dBerry_psi0, &
                      polar=polar, &
                      mo_derivs=mo_derivs, &
                      error=error)
    !
    !
!    polar_section => section_vals_get_subs_vals(qs_env%input, &
!         & "PROPERTIES%LINRES%POLAR",error=error)
    !
    ! Initialize
    ALLOCATE ( polar_tmp( 3, 3), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    polar_tmp(:,:)=0.0_dp

    DO i = 1,3 ! directions of electric field
     DO z=1,3 !dipole directions
       DO ispin=1,dft_control%nspins
!SL compute trace 
            CALL cp_fm_trace(psi1_dBerry(i,ispin)%matrix,mo_derivs(z,ispin)%matrix,polar_tmp(i,z),error=error)
         polar_tmp(i,z)=polar_tmp(i,z) + polar_tmp(i,z)
       END DO
          polar_tmp(i,z) = polar_tmp(i,z)/(twopi * twopi)
      END DO
     END DO !spin
   
     polar(:,:) = MATMUL (MATMUL (cell%hmat,polar_tmp),TRANSPOSE(cell%hmat)) 
!SL evtl maxocc instead?
      IF (dft_control%nspins .eq. 1) THEN
        polar(:,:) = 2.0_dp * polar(:,:)
      END IF

    !
       IF(ASSOCIATED(polar_tmp)) THEN
         DEALLOCATE(polar_tmp, STAT = istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
       
    ENDIF  !raman
    !
    ! Finalize
    CALL timestop(handle)
    !
  END SUBROUTINE polar_polar

  SUBROUTINE polar_print(polar_env,qs_env,error)
    TYPE(polar_env_type)                      :: polar_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'polar_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    INTEGER                                  :: unit_p, output_unit
    LOGICAL                                  :: do_raman, failure
    REAL(dp), DIMENSION(:, : ), &
      POINTER                                :: polar
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: polar_section,input
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: i

    failure = .FALSE.

    NULLIFY(logger,dft_control,para_env)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    polar_section => section_vals_get_subs_vals(qs_env%input, &
         "PROPERTIES%LINRES%POLAR",error=error)

    !
   
    CALL get_polar_env(polar_env=polar_env, &
                     polar=polar, do_raman=do_raman, error=error)
               description="[POLAR]"
       CALL cp_results_erase(results=qs_env%results,description=description,error=error)
       CALL put_results(results=qs_env%results,description=description,&
                        values=polar(:,:),error=error)                     
    CALL get_qs_env(qs_env=qs_env, &
                   dft_control=dft_control, para_env=para_env, error=error)


    !
    IF(BTEST(cp_print_key_should_output(logger%iter_info,qs_env%input,&
             "DFT%PRINT%LINRES",error=error),cp_p_file)) THEN

       unit_p=cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%LINRES",&
                  extension=".data",middle_name="raman",log_filename=.FALSE.,&
                  error=error)
      IF(unit_p>1) THEN 
         write(unit_p,*)
         write (unit_p,*) 'POLARIZABILITY TENSOR (atomic units):'
         write (unit_p,'(A,3F15.5)')"xx,yy,zz", polar(1,1),polar(2,2),polar(3,3)
         write (unit_p,'(A,3F15.5)')"xy,xz,yz",polar(1,2),polar(1,3),polar(2,3)
         write (unit_p,'(A,3F15.5)')"yx,zx,zy",polar(2,1),polar(3,1),polar(3,2)
         write (unit_p,*) 'POLARIZABILITY TENSOR (Angstrom^3):'
         write (unit_p,'(A,3F15.5)')"xx,yy,zz", polar(1,1)*angstrom**3,polar(2,2)*angstrom**3,&
                        polar(3,3)*angstrom**3
         write (unit_p,'(A,3F15.5)')"xy,xz,yz",polar(1,2)*angstrom**3,polar(1,3)*angstrom**3,&
                        polar(2,3)*angstrom**3
         write (unit_p,'(A,3F15.5)')"yx,zx,zy",polar(2,1)*angstrom**3,polar(3,1)*angstrom**3,&
                        polar(3,2)*angstrom**3
       CALL cp_print_key_finished_output(unit_p,logger,qs_env%input,"DFT%PRINT%LINRES",error=error)
      ENDIF
    ENDIF
    IF(do_raman)  THEN
         write (*,*)
         write (*,*) 'POLARIZABILITY TENSOR (atomic units):'
         write (*,*)"xx,yy,zz", polar(1,1),polar(2,2),polar(3,3)
         write (*,*)"xy,xz,yz",polar(1,2),polar(1,3),polar(2,3)
         write (*,'(A,3F15.5)')"xx,yy,zz", polar(1,1),polar(2,2),polar(3,3)
         write (*,'(A,3F15.5)')"xy,xz,yz",polar(1,2),polar(1,3),polar(2,3)
         write (*,'(A,3F15.5)')"yx,zx,zy",polar(2,1),polar(3,1),polar(3,2)
         write (*,*) 'POLARIZABILITY TENSOR (Angstrom^3):'
         write (*,'(A,3F15.5)')"xx,yy,zz", polar(1,1)*angstrom**3 ,&
                 polar(2,2)*angstrom**3 ,polar(3,3)*angstrom**3 
         write (*,'(A,3F15.5)')"xy,xz,yz",polar(1,2)*angstrom**3 ,&
                 polar(1,3)*angstrom**3 ,polar(2,3 )*angstrom**3
         write (*,'(A,3F15.5)')"yx,zx,zy",polar(2,1)*angstrom**3 ,&
                 polar(3,1)*angstrom**3 ,polar(3,2 )*angstrom**3
     END IF
    !
    !
  END SUBROUTINE polar_print


! *****************************************************************************
!> \brief Initialize the polar environment
! *****************************************************************************
  SUBROUTINE polar_env_init(polar_env,qs_env,error)
    !
    TYPE(polar_env_type)                      :: polar_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_env_polar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, idir, ini, ir, &
                                                ispin, istat, m, n, n_rep, &
                                                nao, natom, nspins, &
                                                output_unit,max_states,nmo
    LOGICAL                                  :: failure
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(section_vals_type), POINTER         :: polar_section, lr_section

!

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(linres_control)
    NULLIFY(logger,polar_section)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)


    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                     extension=".linresLog",error=error)

   IF(polar_env%ref_count /= 0) THEN
      CALL polar_env_cleanup(polar_env,error=error)
   ENDIF

   IF(output_unit>0) THEN
      WRITE(output_unit,"(/,T20,A,/)") "*** Start polarizability calculation ***"
      WRITE(output_unit,"(T10,A,/)") "Initialization of the polar environment"
   ENDIF

   CALL polar_env_create(polar_env,error=error)
  !
   polar_section => section_vals_get_subs_vals(qs_env%input, &
                  "PROPERTIES%LINRES%POLAR",error=error)

   CALL get_qs_env(qs_env=qs_env,&
                  dft_control=dft_control,&
                  matrix_s=matrix_s, &
                  linres_control=linres_control,&
                  mos=mos,&
                  error=error)
   !
   !
    nspins = dft_control%nspins
   !
!    CALL section_vals_val_get(polar_section,"DO_RAMAN",l_val=polar_env%do_raman,error=error)
   !
   !
   !
   ALLOCATE(polar_env%polar(3,3), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    polar_env%polar(:,:) = 0.0_dp
  !
      ! allocation
   ALLOCATE(polar_env%dBerry_psi0(3,nspins), &
           polar_env%psi1_dBerry(3,nspins), polar_env%mo_derivs(3,nspins), STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
   DO ispin = 1,nspins
      CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
      CALL cp_fm_get_info(mo_coeff,ncol_global=m,nrow_global=nao,error=error)

      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                              ncol_global=m,&
                              context=mo_coeff%matrix_struct%context,&
                              error=error)
      DO idir = 1,3
         NULLIFY(polar_env%psi1_dBerry(idir,ispin)%matrix,polar_env%dBerry_psi0(idir,ispin)%matrix,&
                 polar_env%mo_derivs(idir,ispin)%matrix)
         CALL cp_fm_create(polar_env%psi1_dBerry(idir,ispin)%matrix,tmp_fm_struct,error=error)
         CALL cp_fm_create(polar_env%dBerry_psi0(idir,ispin)%matrix,tmp_fm_struct,error=error)
         CALL cp_fm_create(polar_env%mo_derivs(idir,ispin)%matrix,tmp_fm_struct,error=error)
      ENDDO
      CALL cp_fm_struct_release(tmp_fm_struct,error=error)

   !
   END DO

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
                                     "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE polar_env_init

! *****************************************************************************
!> \brief Deallocate the polar environment
!> \par History
! *****************************************************************************
  SUBROUTINE polar_env_cleanup(polar_env,error)

    TYPE(polar_env_type)                      :: polar_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir, ispin, istat,i,j
    LOGICAL                                  :: failure

    failure = .FALSE.

    polar_env%ref_count = polar_env%ref_count - 1
    IF(polar_env%ref_count == 0 ) THEN
       IF(ASSOCIATED(polar_env%polar)) THEN
          DEALLOCATE(polar_env%polar,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(polar_env%dBerry_psi0)) THEN
          DO idir = 1,SIZE(polar_env%dBerry_psi0,1)
             DO ispin = 1,SIZE(polar_env%dBerry_psi0,2)
                CALL cp_fm_release(polar_env%dBerry_psi0(idir,ispin)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(polar_env%dBerry_psi0,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(polar_env%mo_derivs)) THEN
          DO idir = 1,SIZE(polar_env%mo_derivs,1)
             DO ispin = 1,SIZE(polar_env%mo_derivs,2)
                CALL cp_fm_release(polar_env%mo_derivs(idir,ispin)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(polar_env%mo_derivs,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(polar_env%psi1_dBerry)) THEN
          DO idir = 1,SIZE(polar_env%psi1_dBerry,1)
             DO ispin = 1,SIZE(polar_env%psi1_dBerry,2)
            write(*,*)'release it'
                CALL cp_fm_release(polar_env%psi1_dBerry(idir,ispin)%matrix,error=error)
             ENDDO
          ENDDO
          DEALLOCATE(polar_env%psi1_dBerry,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       !



    ENDIF ! ref count

  END SUBROUTINE polar_env_cleanup

END MODULE qs_linres_polar_utils
