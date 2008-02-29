!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation and writing of projected density of  states
!>         The DOS is computed per angular momentum and per kind 
!> \par History
!>      - 
!> \author Marcella (29.02.2008,MK)
! *****************************************************************************
MODULE qs_pdos

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE orbital_pointers,                ONLY: nso
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_pdos'

! *****************************************************************************
 ! *** Public subroutines ***

  PUBLIC :: calculate_projected_dos


CONTAINS

! *****************************************************************************
!> \brief   Compute and write projected density of states
!> \author  MI
!> \date    26.02.2008
!> \par History: 
!>       - 
!> \par Variables
!>       - 
!>       - 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_projected_dos(mo_set,atomic_kind_set,particle_set,qs_env,&
                                 dft_section,ispin,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: dft_section
    INTEGER, INTENT(IN), OPTIONAL            :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_projected_dos', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=26)                        :: fmtstr
    CHARACTER(LEN=default_string_length)     :: my_act, my_mittle, my_pos, &
                                                spin(2)
    INTEGER :: handle, iatom, imo, irow, iset, isgf, ishell, iso, istat, &
      iterstep, iw, kind_number, lshell, maxl, maxlgto, mepos, n_dependent, &
      nao, natom, ncol_global, nkind, nrow_global, nset, nsgf
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    LOGICAL                                  :: failure, ionode, should_output
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pdos_array
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_shalf, matrix_shalfc, &
                                                matrix_work, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s_matrix

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    failure = .FALSE.
    should_output = BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
                          "PRINT%PDOS",error=error),cp_p_file)
    spin(1)="ALPHA"
    spin(2)="BETA"
    IF ((.NOT.should_output)) RETURN

    NULLIFY (atomic_kind, context, s_matrix, orb_basis_set, para_env, pdos_array)
    NULLIFY (fm_struct_tmp,matrix_work,matrix_shalf,mo_coeff,vecbuffer)

    CALL timeset(routineN,"I","",handle)
    iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel) 


    CALL get_qs_env(qs_env=qs_env,&
         matrix_s=s_matrix,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
           natom=natom, nsgf=nsgf, maxlgto=maxlgto)
    nkind = SIZE(atomic_kind_set)
    
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff )
    CALL cp_fm_get_info(mo_coeff,&
         context=context, para_env=para_env,&
         nrow_global=nrow_global,&
         ncol_global=ncol_global,error=error)

    ! Create S^1/2 : from sparse to full matrix 
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                  nrow_global=nrow_global,ncol_global=nrow_global,error=error)
    CALL cp_fm_create(matrix_shalf, fm_struct_tmp,name="matrix_shalf",error=error)
    CALL cp_fm_create(matrix_work, fm_struct_tmp,name="matrix_work",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL copy_sm_to_fm(s_matrix(1)%matrix,matrix_shalf,error=error)
    CALL cp_fm_power(matrix_shalf,matrix_work,0.5_dp,EPSILON(0.0_dp),n_dependent,error=error)

    ! Multiply S^(1/2) time the mOS coefficients to get orthonormalized MOS
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
         nrow_global=nrow_global,ncol_global=ncol_global,error=error)
    CALL cp_fm_create(matrix_shalfc,fm_struct_tmp,name="matrix_shalfc",error=error)
    CALL cp_fm_gemm("N","N",nrow_global,ncol_global,nrow_global, &
                   1.0_dp,matrix_shalf,mo_coeff,0.0_dp,matrix_shalfc,error=error)

    CALL cp_fm_release(matrix_work,error=error)
    CALL cp_fm_release(matrix_shalf,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ! Array to store the PDOS per kind and angular momentum
    ALLOCATE(pdos_array(0:maxlgto,nkind,mo_set%nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
    pdos_array = 0.0_dp
    nao = mo_set%nao
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp

    DO imo =  1,mo_set%nmo
      
      ! Extract the eigenvector from the distributed full matrix
      CALL cp_fm_get_submatrix(matrix_shalfc,vecbuffer,1,imo,&
                     nao,1,transpose=.TRUE.,error=error)
      irow = 1
      DO iatom = 1,natom

          NULLIFY(orb_basis_set)
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
               kind_number=kind_number,&
               orb_basis_set=orb_basis_set)

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                          nset=nset,&
                          nshell=nshell,&
                          l=l)

          isgf = 1
          DO iset=1,nset
            DO ishell=1,nshell(iset)
               lshell = l(ishell,iset)
               DO iso=1,nso(lshell)
                 pdos_array(lshell,kind_number,imo) = &
                      pdos_array(lshell,kind_number,imo) +&
                      vecbuffer(1,irow) * vecbuffer(1,irow) 
                 irow = irow + 1
               END DO  ! iso
            END DO  ! ishell
          END DO  ! iset

      END DO  ! iatom

    END DO  ! imo

    CALL cp_fm_release(matrix_shalfc,error=error)
    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    my_pos = "REWIND"
    my_act = "WRITE"
    DO kind_number = 1,nkind
    
      NULLIFY(atomic_kind, orb_basis_set)
      atomic_kind => atomic_kind_set(kind_number)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
           orb_basis_set=orb_basis_set)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
           maxl=maxl)

      IF(PRESENT(ispin)) THEN
        my_mittle =TRIM(spin(ispin))//"_k"//TRIM(ADJUSTL(cp_to_string(kind_number)))
      ELSE
        my_mittle ="k"//TRIM(ADJUSTL(cp_to_string(kind_number)))
      END IF
      
      iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PDOS",&
           extension=".pdos", file_position=my_pos, file_action=my_act,&
           file_form="FORMATTED", middle_name=TRIM(my_mittle), &
           error=error)
      IF(iw>0) THEN

        fmtstr = "(I8,2X,F16.6,  (2X,F16.8))"
        WRITE (UNIT=fmtstr(14:15),FMT="(I2)") maxl+1
 
        WRITE(UNIT=iw,FMT="(A,I5,A,I10)") "#  Projected DOS for kind number ", &
              kind_number, " at iteration step i = ", iterstep
        DO imo = 1,mo_set%nmo
          WRITE(UNIT=iw,FMT=fmtstr) imo, mo_set%eigenvalues(imo),&
                (pdos_array(lshell,kind_number,imo), lshell=0,maxl)
        END DO
      END IF
      CALL cp_print_key_finished_output(iw,logger,dft_section,&
                "PRINT%PDOS", error=error)

    END DO  ! ikind

    DEALLOCATE(pdos_array, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE  calculate_projected_dos

END MODULE qs_pdos

