!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/molecular_states [1.0] *
!!
!!   NAME
!!     molecular_states
!!
!!   FUNCTION
!!     Routines for the calculation of molecular states
!!
!!   AUTHOR
!!     Manu Manu (transcribed by CJM) 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE molecular_states
! *****************************************************************************
  USE coefficient_types,               ONLY: coeff_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_type,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm, &
                                             cp_fm_write
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp
  USE molecule_types,                  ONLY: molecule_structure_type
  USE pw_env_types,                    ONLY: pw_env_type 
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: qs_environment_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop,&
                                             trace_debug
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "molecular_states"

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Public subroutines ***

  PUBLIC :: construct_molecular_states

! *****************************************************************************

CONTAINS

! *****************************************************************************

        SUBROUTINE construct_molecular_states(molecule,mo_eigenvectors, &
                                   Hks,globenv,pw_env,qs_env,wf_r,wf_g, &
                                   mgrid_rspace,mgrid_gspace)

        IMPLICIT NONE

        TYPE ( molecule_structure_type ), DIMENSION ( : ), POINTER      :: molecule
        TYPE ( cp_fm_type ), POINTER                          :: mo_eigenvectors
        TYPE ( real_matrix_type ), POINTER                    :: Hks
        TYPE ( global_environment_type ), INTENT ( IN )       :: globenv
        TYPE ( pw_env_type ), POINTER                         :: pw_env
        TYPE ( qs_environment_type ), INTENT ( IN )           :: qs_env
        TYPE ( coeff_type ), INTENT ( INOUT )                 :: wf_r, wf_g
        TYPE ( coeff_type ), DIMENSION ( : ), POINTER         :: mgrid_gspace,&
                                                                 mgrid_rspace
! locals
        CHARACTER (LEN=20)                              :: filename
        INTEGER                                         :: i, j, k, imol, isos, &
                                                           nrow_global, ncol_global, ns
        INTEGER, DIMENSION ( : ), POINTER               :: ind
        REAL(wp), ALLOCATABLE                           :: evals( : )
        TYPE(cp_fm_type), POINTER                       :: b,c,d
        TYPE(cp_fm_type), POINTER                       :: e_vectors,&
                                                           rot_e_vectors
        TYPE(cp_fm_struct_type), POINTER                :: fm_struct_tmp
        TYPE(cp_error_type)                             :: error

        DO imol = 1, SIZE(molecule)
          ind => molecule ( imol ) % states
          ns = molecule ( imol ) % nstates
          ALLOCATE ( evals ( ns ), STAT = isos )
          IF ( isos /= 0 ) CALL stop_memory ( "construct_molecular_states", &
                                              "eval", ns )

          NULLIFY(b,c,d,fm_struct_tmp,e_vectors,rot_e_vectors)

          CALL cp_fm_get_info(mo_eigenvectors, &
                           ncol_global=ncol_global, &
                           nrow_global=nrow_global )
          
          CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_global,  ncol_global=ns, &
                                  para_env=mo_eigenvectors%matrix_struct%para_env, &
                                  context=mo_eigenvectors%matrix_struct%context)

          CALL cp_fm_create(b,fm_struct_tmp, name="b")
          CALL cp_fm_create(c,fm_struct_tmp, name="c")
          CALL cp_fm_create(rot_e_vectors,fm_struct_tmp, name="rot_e_vectors")

          CALL cp_fm_struct_release(fm_struct_tmp)

          CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns,  ncol_global=ns, &
                                  para_env=mo_eigenvectors%matrix_struct%para_env, &
                                  context=mo_eigenvectors%matrix_struct%context)

          CALL cp_fm_create(d,fm_struct_tmp, name="d")
          CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e_vectors")
          CALL cp_fm_struct_release(fm_struct_tmp)
                
          DO i=1,ns
            CALL cp_fm_to_fm ( mo_eigenvectors, b, 1, ind ( i ), i)
          END DO

          CALL cp_sm_fm_multiply(Hks,b,c,ns,para_env=qs_env%para_env)     

          CALL cp_fm_gemm('T','N',ns,ns,nrow_global,1.0_wp, &             
                         b,c,0.0_wp,d)    

          CALL cp_fm_syevd ( d, e_vectors, evals )

          IF (globenv%ionode) WRITE(globenv%scr,*)"" 
          IF (globenv%ionode) WRITE(globenv%scr,*)"MOLECULE ",imol
          IF (globenv%ionode) WRITE(globenv%scr,*)"NUMBER OF STATES  ", ns
          IF (globenv%ionode) WRITE(globenv%scr,*)"EIGENVALUES"
          IF (globenv%ionode) WRITE(globenv%scr,*)"" 

          DO k=1,ns
             IF (globenv%ionode) WRITE(globenv%scr,*) evals(k)
          END DO
                  
          CALL cp_fm_gemm('N','N',nrow_global,ns,ns,1.0_wp, &                             
                          b,e_vectors,0.0_wp,rot_e_vectors)

          DEALLOCATE ( evals,  STAT = isos )
          IF ( isos /= 0 ) CALL stop_memory ( "construct_molecular_states", &
                                              "eval", ns )
          
          IF ( globenv % print % cubes ) THEN
            DO i=1,ns
              j=ind(i)
              CALL calculate_wavefunction(rot_e_vectors,i,wf_r,              &
                                     wf_g,mgrid_rspace, mgrid_gspace,        &
                                     pw_env%gridlevel_info,pw_env%cube_info, &
                                     pw_env%l_info, qs_env, error = error )


              WRITE(filename,'(a9,I2.2,a1,I3.3,a5)')  &
              "MOLECULE_",imol,"_",j,".cube"
              IF (globenv%ionode) OPEN(UNIT=22,FILE=filename)
              CALL rs_pw_to_cube(wf_r%pw,22,globenv%ionode,filename)
              IF (globenv%ionode) CLOSE ( 22 )
              IF (globenv%ionode) WRITE(globenv%scr,*) filename
            END DO
          ENDIF  

          CALL cp_fm_release ( b )
          CALL cp_fm_release ( c )
          CALL cp_fm_release ( d )
          CALL cp_fm_release ( e_vectors )
          CALL cp_fm_release ( rot_e_vectors )

        END DO 
!--------------------------------------------------------------------------------------------------
        END SUBROUTINE construct_molecular_states
!--------------------------------------------------------------------------------------------------
END MODULE molecular_states
