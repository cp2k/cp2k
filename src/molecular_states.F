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
                                             cp_fm_get_element,&
                                             cp_fm_release,&
                                             cp_fm_to_fm, &
                                             cp_fm_write
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_max,mp_maxloc,mp_sum,&
                                             mp_allgather,mp_bcast
  USE molecule_types_new,              ONLY: molecule_type
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

SUBROUTINE construct_molecular_states(molecule_set, mo_localized, &
                                   mo_eigenvectors, mo_eigenvalues, &
                                   Hks,matrix_S,globenv,pw_env,qs_env,wf_r,wf_g, &
                                   mgrid_rspace,mgrid_gspace)

  IMPLICIT NONE

  TYPE ( molecule_type ), DIMENSION ( : ), POINTER      :: molecule_set
  TYPE ( cp_fm_type ), POINTER                          :: mo_eigenvectors, mo_localized
  REAL (KIND = dp), POINTER, DIMENSION(:)                    :: mo_eigenvalues
  TYPE ( real_matrix_type ), POINTER                    :: Hks,matrix_S
  TYPE ( global_environment_type ), INTENT ( IN )       :: globenv
  TYPE ( pw_env_type ), POINTER                         :: pw_env
  TYPE ( qs_environment_type ), INTENT ( IN )           :: qs_env
  TYPE ( coeff_type ), INTENT ( INOUT )                 :: wf_r, wf_g
  TYPE ( coeff_type ), DIMENSION ( : ), POINTER         :: mgrid_gspace,&
                                                           mgrid_rspace
  
! locals
  TYPE ( cp_para_env_type ), POINTER :: para_env
  CHARACTER (LEN=default_string_length)                            :: filename
  INTEGER                                       :: i, j, k, imol, isos, &
                                                   nrow_global, ncol_global, ns
  REAL(KIND = dp), ALLOCATABLE                           :: evals( : )
  TYPE(cp_fm_type), POINTER                       :: b,c,d
  TYPE(cp_fm_type), POINTER                       :: e_vectors,&
                                                     rot_e_vectors
  TYPE(cp_fm_type), POINTER                       :: D_igk
  TYPE ( cp_fm_type ), POINTER                    :: smo
  TYPE(cp_fm_struct_type), POINTER                :: fm_struct_tmp
  TYPE(cp_error_type)                             :: error

  INTEGER  :: nproc, iproc,nmol,istates
  INTEGER, POINTER :: nstates(:)
  INTEGER, POINTER :: states(:)
  INTEGER, DIMENSION ( : ), POINTER               :: ind
  INTEGER :: unit_report
  REAL(KIND = dp) :: tmp
  INTEGER :: ikind,nkind
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. 
  !-----------------------------------------------------------------------------
  para_env => qs_env % para_env
  nproc = para_env%num_pe
 

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. 
  !-----------------------------------------------------------------------------
  IF (globenv%ionode) THEN 
     unit_report=171
     OPEN(unit_report,FILE="Molecular_DOS")
     WRITE(unit_report,*) SIZE(mo_eigenvalues)," number of states " 
     DO i=1,SIZE(mo_eigenvalues)
        WRITE(unit_report,*) mo_eigenvalues(i)
     ENDDO
  ENDIF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. 
  !-----------------------------------------------------------------------------
  CALL cp_fm_get_info(mo_localized, &
                      ncol_global=ncol_global, &
                      nrow_global=nrow_global )
  NULLIFY(smo)
  CALL cp_fm_create(smo,mo_eigenvectors%matrix_struct)
  CALL cp_sm_fm_multiply(matrix_S,mo_eigenvectors,smo,ncol_global,qs_env%para_env)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. 
  !-----------------------------------------------------------------------------
  ALLOCATE(nstates(2),STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nstates')

  DO imol = 1, SIZE(molecule_set)
    nstates(1) = molecule_set(imol)%lmi%nstates
    nstates(2) = para_env%mepos

    CALL mp_maxloc(nstates,para_env%group)

    NULLIFY(states)
    ALLOCATE(states(nstates(1)),STAT=isos)
    IF(isos/=0) CALL stop_memory('construct_molecular_states','states')
    states(:) = 0

    iproc=nstates(2)
    IF(iproc == para_env%mepos) THEN
      states(:) = molecule_set(imol)%lmi%states(:)
    END IF
    !!BCAST from here root = iproc
    CALL mp_bcast(states,iproc,para_env%group)

    ns = nstates(1)
    ind => states( : )

    ALLOCATE ( evals ( ns ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "construct_molecular_states", &
                                        "eval", ns )

    NULLIFY(b,c,d,fm_struct_tmp,e_vectors,rot_e_vectors)

          
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_global, &
                             ncol_global=ns, &
                             para_env=mo_localized%matrix_struct%para_env, &
                             context=mo_localized%matrix_struct%context)

    CALL cp_fm_create(b,fm_struct_tmp, name="b")
    CALL cp_fm_create(c,fm_struct_tmp, name="c")
    CALL cp_fm_create(rot_e_vectors,fm_struct_tmp, name="rot_e_vectors")

    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns,  ncol_global=ns, &
                             para_env=mo_localized%matrix_struct%para_env, &
                             context=mo_localized%matrix_struct%context)

    CALL cp_fm_create(d,fm_struct_tmp, name="d")
    CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e_vectors")
    CALL cp_fm_struct_release(fm_struct_tmp)
                
    DO i=1,ns
      CALL cp_fm_to_fm ( mo_localized, b, 1, ind ( i ), i)
    END DO

    CALL cp_sm_fm_multiply(Hks,b,c,ns,para_env=qs_env%para_env)     

    CALL cp_fm_gemm('T','N',ns,ns,nrow_global,1.0_dp, &             
                         b,c,0.0_dp,d)    

    CALL cp_fm_syevd ( d, e_vectors, evals )

    IF (globenv%ionode) WRITE(globenv%scr,*)"" 
    IF (globenv%ionode) WRITE(globenv%scr,*)"MOLECULE ",imol
    IF (globenv%ionode) WRITE(globenv%scr,*)"NUMBER OF STATES  ", ns
    IF (globenv%ionode) WRITE(globenv%scr,*)"EIGENVALUES"
    IF (globenv%ionode) WRITE(globenv%scr,*)"" 

    DO k=1,ns
      IF (globenv%ionode) WRITE(globenv%scr,*) evals(k)
    END DO
    IF (globenv%ionode) THEN
       WRITE(unit_report,*) imol,ns, " imol, number of states"
       DO k=1,ns
          WRITE(unit_report,*) evals(k)
       ENDDO
    ENDIF
                  
    CALL cp_fm_gemm('N','N',nrow_global,ns,ns,1.0_dp, &
                    b,e_vectors,0.0_dp,rot_e_vectors)

    DEALLOCATE ( evals,  STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "construct_molecular_states", &
                                        "eval", ns )

    ! compute Eq. 6 from P. Hunt et al. (CPL 376, p. 68-74)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns, &
                             ncol_global=ncol_global, &
                             para_env=mo_localized%matrix_struct%para_env, &
                             context=mo_localized%matrix_struct%context)
    CALL cp_fm_create(D_igk,fm_struct_tmp)
    CALL cp_fm_struct_release(fm_struct_tmp) 
    CALL cp_fm_gemm('T','N',ns,ncol_global,nrow_global,1.0_dp, &             
                         rot_e_vectors,smo,0.0_dp,D_igk)    
    DO i=1,ns
       DO k=1,ncol_global
          CALL cp_fm_get_element(D_igk,i,k,tmp)
          IF (globenv%ionode) THEN
             write(unit_report,*) tmp**2
          ENDIF
       ENDDO
    ENDDO
    CALL cp_fm_release(D_igk)
          
    IF ( globenv % print % cubes ) THEN
      DO i=1,ns
        j=ind(i)
        CALL calculate_wavefunction(rot_e_vectors,i,wf_r, &
                                    wf_g,mgrid_rspace, mgrid_gspace, &
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

    DEALLOCATE(states,STAT=isos)
    IF(isos/=0) CALL stop_memory('construct_molecular_states','states')

  END DO 
  DEALLOCATE(nstates,STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nstates')
  IF (globenv%ionode) CLOSE(unit_report)
  !CALL stop_program("molecular_states","checkpoint")
!------------------------------------------------------------------------------
END SUBROUTINE construct_molecular_states

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END MODULE molecular_states
