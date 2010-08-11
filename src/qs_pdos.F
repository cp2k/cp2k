!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: nso,&
                                             nsoset
  USE orbital_symbols,                 ONLY: l_sym,&
                                             sgf_symbol
  USE particle_types,                  ONLY: particle_type
  USE preconditioner_types,            ONLY: preconditioner_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_type
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
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
                                 scf_env,dft_section,ispin,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(section_vals_type), POINTER         :: dft_section
    INTEGER, INTENT(IN), OPTIONAL            :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_projected_dos', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=16)                        :: fmtstr2
    CHARACTER(LEN=27)                        :: fmtstr1
    CHARACTER(LEN=6), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: tmp_str
    CHARACTER(LEN=default_string_length)     :: kind_name, my_act, my_mittle, &
                                                my_pos, spin(2)
    INTEGER :: handle, homo, i, iatom, il, im, imo, ir, irow, iset, isgf, &
      ishell, iso, iterstep, iw, j, kind_number, lcomponent, lshell, maxl, &
      maxl_list, maxlgto, my_spin, n_dependent, n_rep, nao, natom, &
      ncol_global, nkind, nlist, nlumo, nmo, nrow_global, nset, nsgf, nvirt, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: index_list, list, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    LOGICAL :: add_to_list, append, do_list, do_virt, failure, ionode, &
      separate_components, should_output
    REAL(KIND=dp)                            :: e_fermi
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, evals_virt, &
                                                occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pdos_array_list, vecbuffer
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pdos_array
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: s_matrix
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_shalf, matrix_shalfc, &
                                                matrix_work, mo_coeff, mo_virt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

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
    NULLIFY (eigenvalues,fm_struct_tmp,matrix_work,matrix_shalf,mo_coeff,vecbuffer)
    NULLIFY (index_list, list, pdos_array_list)
    NULLIFY (occupation_numbers)

    CALL timeset(routineN,handle)
    iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel) 


    CALL get_qs_env(qs_env=qs_env,&
         matrix_s=s_matrix,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
           natom=natom, nsgf=nsgf, maxlgto=maxlgto)
    nkind = SIZE(atomic_kind_set)
    
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, homo=homo, nao=nao, nmo=nmo,&
                    mu=e_fermi)
    CALL cp_fm_get_info(mo_coeff,&
         context=context, para_env=para_env,&
         nrow_global=nrow_global,&
         ncol_global=ncol_global,error=error)

    
    CALL section_vals_val_get(dft_section,"PRINT%PDOS%NLUMO",i_val=nlumo,error=error)
    IF (nlumo == -1) nlumo = nao - homo
    do_virt = (nlumo>(nmo-homo))
    nvirt = nlumo - (nmo-homo)
    ! Generate virtual orbitals
    IF(do_virt) THEN 
       IF(PRESENT(ispin)) THEN
         my_spin = ispin
       ELSE
         my_spin = 1
       END IF

       CALL generate_virtual_mo(qs_env,scf_env,mo_set, evals_virt, mo_virt, nvirt, ispin=my_spin, error=error )
    ELSE
       NULLIFY (evals_virt,mo_virt)
       nvirt = 0
    END IF

    ! Create S^1/2 : from sparse to full matrix 
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                  nrow_global=nrow_global,ncol_global=nrow_global,error=error)
    CALL cp_fm_create(matrix_shalf, fm_struct_tmp,name="matrix_shalf",error=error)
    CALL cp_fm_create(matrix_work, fm_struct_tmp,name="matrix_work",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL copy_dbcsr_to_fm(s_matrix(1)%matrix,matrix_shalf,error=error)
    CALL cp_fm_power(matrix_shalf,matrix_work,0.5_dp,EPSILON(0.0_dp),n_dependent,error=error)
    CALL cp_fm_release(matrix_work,error=error)

    ! Multiply S^(1/2) time the mOS coefficients to get orthonormalized MOS
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
         nrow_global=nrow_global,ncol_global=ncol_global,error=error)
    CALL cp_fm_create(matrix_shalfc,fm_struct_tmp,name="matrix_shalfc",error=error)
    CALL cp_fm_gemm("N","N",nrow_global,ncol_global,nrow_global, &
                   1.0_dp,matrix_shalf,mo_coeff,0.0_dp,matrix_shalfc,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    IF (do_virt) THEN
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
            nrow_global=nrow_global,ncol_global=nvirt,error=error)
      CALL cp_fm_create(matrix_work,fm_struct_tmp,name="matrix_shalfc",error=error)
      CALL cp_fm_gemm("N","N",nrow_global,nvirt,nrow_global, &
                     1.0_dp,matrix_shalf,mo_virt,0.0_dp,matrix_work,error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    END IF

    CALL cp_fm_release(matrix_shalf,error=error)
    ! Array to store the PDOS per kind and angular momentum
    do_list = .FALSE.
    CALL section_vals_val_get(dft_section,"PRINT%PDOS%LIST", n_rep_val=n_rep, error=error)  
    IF(n_rep>0) THEN
       do_list = .TRUE.
       nlist = 0
       DO ir = 1,n_rep
          NULLIFY (list)
          CALL section_vals_val_get(dft_section,"PRINT%PDOS%LIST",i_rep_val=ir,i_vals=list, error=error)
          IF(ASSOCIATED(list)) THEN
             CALL reallocate(index_list,1,nlist+SIZE(list))
             DO i = 1,SIZE(list)
                index_list(i+nlist) = list(i)
             END DO
             nlist = nlist + SIZE(list)
          END IF
       END DO
    END IF
    CALL section_vals_val_get(dft_section,"PRINT%PDOS%COMPONENTS",l_val=separate_components,error=error)
    IF (separate_components) THEN
      ALLOCATE(pdos_array(nsoset(maxlgto),nkind,nmo+nvirt),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
      IF(do_list) THEN
         ALLOCATE(pdos_array_list(nsoset(maxlgto),nmo+nvirt),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
    ELSE
      ALLOCATE(pdos_array(0:maxlgto,nkind,nmo+nvirt),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
      IF(do_list) THEN
         ALLOCATE(pdos_array_list(0:maxlgto,nmo+nvirt),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF
    IF(do_virt) THEN
       ALLOCATE(eigenvalues(nmo+nvirt),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       eigenvalues(1:nmo) = mo_set%eigenvalues(1:nmo)
       eigenvalues(nmo+1:nmo+nvirt) = evals_virt(1:nvirt)
       ALLOCATE(occupation_numbers(nmo+nvirt),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       occupation_numbers(:) = 0.0_dp
       occupation_numbers(1:nmo) = mo_set%occupation_numbers(1:nmo)
    ELSE
       eigenvalues => mo_set%eigenvalues
       occupation_numbers => mo_set%occupation_numbers
    END IF

    pdos_array = 0.0_dp
    IF(do_list) pdos_array_list = 0.0_dp
    maxl_list = -1
    nao = mo_set%nao
    ALLOCATE(vecbuffer(1,nao),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp

    DO imo =  1,nmo+nvirt
      
      ! Extract the eigenvector from the distributed full matrix
      IF(imo>nmo) THEN
        CALL cp_fm_get_submatrix(matrix_work,vecbuffer,1,imo-nmo,&
                       nao,1,transpose=.TRUE.,error=error)
      ELSE
        CALL cp_fm_get_submatrix(matrix_shalfc,vecbuffer,1,imo,&
                       nao,1,transpose=.TRUE.,error=error)
      END IF
      irow = 1
      DO iatom = 1,natom

          NULLIFY(orb_basis_set)
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
               kind_number=kind_number,&
               orb_basis_set=orb_basis_set)

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                          nset=nset,&
                          nshell=nshell,&
                          l=l, maxl=maxl)
           add_to_list = .FALSE.                
          IF(do_list) THEN
             DO il = 1,nlist
                IF(iatom==index_list(il))  THEN
                   add_to_list = .TRUE.
                   maxl_list = MAX(maxl,maxl_list)
                   EXIT
                END IF
             ENDDO
          END IF
          IF (separate_components) THEN
            isgf = 1
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                 lshell = l(ishell,iset)
                 DO iso=1,nso(lshell)
                   lcomponent = nsoset(lshell-1)+iso
                   pdos_array(lcomponent,kind_number,imo) = &
                        pdos_array(lcomponent,kind_number,imo) +&
                        vecbuffer(1,irow) * vecbuffer(1,irow) 
                   IF(add_to_list) pdos_array_list(lcomponent,imo) = &
                             pdos_array_list(lcomponent,imo) +&
                             vecbuffer(1,irow) * vecbuffer(1,irow) 
                   irow = irow + 1
                 END DO  ! iso
              END DO  ! ishell
            END DO  ! iset
          ELSE
            isgf = 1
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                 lshell = l(ishell,iset)
                 DO iso=1,nso(lshell)
                   pdos_array(lshell,kind_number,imo) = &
                        pdos_array(lshell,kind_number,imo) +&
                        vecbuffer(1,irow) * vecbuffer(1,irow) 
                   IF(add_to_list) pdos_array_list(lshell,imo) = &
                             pdos_array_list(lshell,imo) +&
                             vecbuffer(1,irow) * vecbuffer(1,irow) 
                   irow = irow + 1
                 END DO  ! iso
              END DO  ! ishell
            END DO  ! iset
          END IF

      END DO  ! iatom

    END DO  ! imo

    CALL cp_fm_release(matrix_shalfc,error=error)
    DEALLOCATE(vecbuffer,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    CALL section_vals_val_get(dft_section,"PRINT%PDOS%APPEND",l_val=append,error=error)
    IF(append .AND. iterstep > 1 ) THEN
      my_pos = "APPEND"
    ELSE
      my_pos = "REWIND"
    END IF
    my_act = "WRITE"
    DO kind_number = 1,nkind
    
      NULLIFY(atomic_kind, orb_basis_set)
      atomic_kind => atomic_kind_set(kind_number)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
           orb_basis_set=orb_basis_set,name=kind_name)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
           maxl=maxl) 

      ! basis none has no associated maxl, and no pdos
      IF (maxl<0) CYCLE

      IF(PRESENT(ispin)) THEN
        my_mittle =TRIM(spin(ispin))//"_k"//TRIM(ADJUSTL(cp_to_string(kind_number)))
        my_spin = ispin
      ELSE
        my_mittle ="k"//TRIM(ADJUSTL(cp_to_string(kind_number)))
        my_spin = 1
      END IF
      
      iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PDOS",&
           extension=".pdos", file_position=my_pos, file_action=my_act,&
           file_form="FORMATTED", middle_name=TRIM(my_mittle), &
           error=error)
      IF(iw>0) THEN

        fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
        fmtstr2 = "(A42,  (10X,A8))"
        IF(separate_components) THEN
          WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") nsoset(maxl)
          WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") nsoset(maxl)
        ELSE
          WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") maxl+1
          WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") maxl+1
        END IF
 
        WRITE (UNIT=iw,FMT="(A,I0,A,F12.6,A)")&
          "# Projected DOS for atomic kind "//TRIM(kind_name)//" at iteration step i = ",&
          iterstep,", E(Fermi) = ",e_fermi," a.u."
        IF (separate_components) THEN
          ALLOCATE(tmp_str(0:0,0:maxl,-maxl:maxl),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          tmp_str = ""
          DO j = 0, maxl
             DO i = -j,j
                tmp_str(0,j,i) = sgf_symbol(0,j,i)
             END DO
          END DO

          WRITE (UNIT=iw,FMT=fmtstr2)&
            "#     MO Eigenvalue [a.u.]      Occupation",&
            ((TRIM(tmp_str(0,il,im)),im=-il,il),il=0,maxl)
          DO imo = 1,nmo+nvirt
            WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
              (pdos_array(lshell,kind_number,imo),lshell=1,nsoset(maxl))
          END DO
          DEALLOCATE(tmp_str,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ELSE
          WRITE (UNIT=iw,FMT=fmtstr2)&
            "#     MO Eigenvalue [a.u.]      Occupation",&
            (TRIM(l_sym(il)),il=0,maxl)
          DO imo = 1,nmo+nvirt
            WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
              (pdos_array(lshell,kind_number,imo),lshell=0,maxl)
          END DO
        END IF
      END IF
      CALL cp_print_key_finished_output(iw,logger,dft_section,&
                "PRINT%PDOS", error=error)

    END DO  ! ikind


   ! basis none has no associated maxl, and no pdos
    IF (maxl_list>0) THEN

      IF(PRESENT(ispin)) THEN
        my_mittle =TRIM(spin(ispin))//"_list"
        my_spin = ispin
      ELSE
        my_mittle ="list"
        my_spin = 1
      END IF
      
      iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PDOS",&
           extension=".pdos", file_position=my_pos, file_action=my_act,&
           file_form="FORMATTED", middle_name=TRIM(my_mittle), &
           error=error)
      IF(iw>0) THEN

        fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
        fmtstr2 = "(A42,  (10X,A8))"
        IF(separate_components) THEN
          WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") nsoset(maxl_list)
          WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") nsoset(maxl_list)
        ELSE
          WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") maxl_list+1
          WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") maxl_list+1
        END IF
 
        WRITE (UNIT=iw,FMT="(A,I0,A,I0,A,F12.6,A)")&
          "# Projected DOS for a list of ",nlist," atoms, at iteration step i = ",iterstep,&
          ", E(Fermi) = ",e_fermi," a.u."
        IF(separate_components) THEN
          ALLOCATE(tmp_str(0:0,0:maxl_list,-maxl_list:maxl_list),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          tmp_str = ""
          DO j = 0, maxl_list
             DO i = -j,j
                tmp_str(0,j,i) = sgf_symbol(0,j,i)
             END DO
          END DO

          WRITE (UNIT=iw,FMT=fmtstr2)&
            "#     MO Eigenvalue [a.u.]      Occupation",&
            ((TRIM(tmp_str(0,il,im)),im=-il,il),il=0,maxl_list)
          DO imo = 1,nmo+nvirt
            WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
              (pdos_array_list(lshell,imo),lshell=1,nsoset(maxl_list))
          END DO
          DEALLOCATE(tmp_str,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ELSE
          WRITE (UNIT=iw,FMT=fmtstr2)&
            "#     MO Eigenvalue [a.u.]      Occupation",&
            (TRIM(l_sym(il)),il=0,maxl_list)
          DO imo = 1,nmo+nvirt
            WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
              (pdos_array_list(lshell,imo),lshell=0,maxl_list)
          END DO
        END IF
      END IF
      CALL cp_print_key_finished_output(iw,logger,dft_section,&
                "PRINT%PDOS", error=error)

    END IF  ! do_list

    DEALLOCATE(pdos_array, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
    IF(do_list) THEN
       DEALLOCATE(pdos_array_list, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(index_list, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(do_virt) THEN
      DEALLOCATE(evals_virt,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
      CALL cp_fm_release(mo_virt,error=error)
      CALL cp_fm_release(matrix_work,error=error)
      DEALLOCATE(eigenvalues,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
      DEALLOCATE(occupation_numbers,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE  calculate_projected_dos

! *****************************************************************************
!> \brief   Compute additional virtual states  starting from the available MOS
!> \author  MI
!> \date    08.03.2008
!> \par History: 
!>       - 
!> \par Variables
!>       - 
!>       - 
!> \version 1.0
! *****************************************************************************
  
  SUBROUTINE generate_virtual_mo (qs_env, scf_env, mo_set, evals_virt, mo_virt, &
             nvirt, ispin, error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals_virt
    TYPE(cp_fm_type), POINTER                :: mo_virt
    INTEGER, INTENT(IN)                      :: nvirt, ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_virtual_mo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nmo, nrow_global, stat
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, s_matrix
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(preconditioner_type), POINTER       :: local_preconditioner
    TYPE(scf_control_type), POINTER          :: scf_control

    failure = .FALSE.
    NULLIFY(evals_virt,mo_virt)
    ALLOCATE(evals_virt(nvirt),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 

    CALL get_qs_env(qs_env,matrix_ks=ks_matrix, matrix_s=s_matrix, &
         scf_control=scf_control, error=error)
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, nmo=nmo )
    CALL cp_fm_get_info(mo_coeff, context=context, para_env=para_env,&
         nrow_global=nrow_global, error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
         nrow_global=nrow_global,ncol_global=nvirt,error=error)
    CALL cp_fm_create(mo_virt, fm_struct_tmp,name="virtual",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_init_random(mo_virt,nvirt,error=error)

    NULLIFY (local_preconditioner)

    CALL ot_eigensolver(matrix_h=ks_matrix(ispin)%matrix,matrix_s=s_matrix(1)%matrix, &
         matrix_c_fm=mo_virt, matrix_orthogonal_space_fm=mo_coeff, &
         eps_gradient=scf_control%eps_lumos, &
         preconditioner=local_preconditioner, &
         iter_max=scf_control%max_iter_lumos,&
         size_ortho_space=nmo,error=error)

    CALL calculate_subspace_eigenvalues(mo_virt,ks_matrix(ispin)%matrix,&
         evals_virt, error=error)


  END SUBROUTINE generate_virtual_mo

END MODULE qs_pdos

