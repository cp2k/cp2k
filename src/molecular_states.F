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
                                             dp
  USE message_passing,                 ONLY: mp_max,mp_minloc,mp_sum,&
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
  CHARACTER (LEN=20)                            :: filename
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

  INTEGER  :: my_nmol, nproc, counter,iproc,nstates_max,nmol,istates
  INTEGER, POINTER  :: tmp_index(:)
  INTEGER, POINTER  :: nstates_array(:),nmol_array(:)
  INTEGER, DIMENSION ( : ), POINTER               :: ind
  INTEGER, POINTER  :: states_index(:,:)
  INTEGER :: unit_report
  REAL(KIND = dp) :: tmp
  TYPE(distribution_1d_type), POINTER :: local_molecules
  INTEGER :: ikind,nkind

  IF (globenv%ionode) THEN 
     unit_report=171
     OPEN(unit_report,FILE="Molecular_DOS")
     WRITE(unit_report,*) SIZE(mo_eigenvalues)," number of states " 
     DO i=1,SIZE(mo_eigenvalues)
        WRITE(unit_report,*) mo_eigenvalues(i)
     ENDDO
  ENDIF

#if defined ( __parallel )

  para_env => qs_env % para_env
  local_molecules => qs_env % subsys % local_molecules_new

  nproc = para_env%num_pe

  !WRITE(*,*) "nproc = ",nproc

  my_nmol = 0
  nkind = SIZE ( local_molecules%list )
  DO ikind=1,nkind
    nmol = SIZE ( local_molecules%list(ikind)%array )
    my_nmol = my_nmol + nmol
  END DO
  nmol = my_nmol

  !WRITE(*,*) "my_nmol = ",my_nmol

  CALL mp_sum(nmol,para_env%group)

  !WRITE(*,*) "nmol = ",nmol

  ALLOCATE(nstates_array(nmol),STAT=isos) !hold the nstates each imol has
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nstates_array')
  ALLOCATE(nmol_array(nproc),STAT=isos)   !hold the nmol in each procs
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nmol_array')

  nstates_array(:) = 0

  counter = 0
  DO ikind=1,nkind
    DO imol=1, SIZE ( local_molecules%list(ikind)%array )
      counter = counter + 1
      i = local_molecules%list(ikind)%array(imol)
      nstates_array(counter) = molecule_set(i)%lmi%nstates
    END DO
  END DO
  !WRITE(*,*) "counter = ",counter

  nstates_max = MAXVAL(nstates_array)
  CALL mp_max(nstates_max,para_env%group)
  !WRITE(*,*) "nstates_max = ",nstates_max

  ALLOCATE(states_index(1:nmol,1:nstates_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','states_index')
  states_index(:,:) = 0

  !Set up nmol_array on all procs
  CALL mp_allgather(my_nmol,nmol_array,para_env%group)
  !WRITE(*,*) "nmol_array ",nmol_array

  !WRITE(*,*) "AAAAAAAAAAAA"
  ALLOCATE(tmp_index(nstates_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','tmp_index')
  counter = 0
  !WRITE(*,*) "BAAAAAAAAAAA"
  DO iproc=0,nproc-1
    DO imol = 1,nmol_array(iproc+1)
      counter = counter + 1
      tmp_index(:) = 0
      IF(iproc == para_env%mepos) THEN
        i = local_molecules%list(1)%array(imol)
        tmp_index(:) = molecule_set(i)%lmi%states
      END IF
      !BCAST from here root = iproc
      CALL mp_bcast(tmp_index,iproc,para_env%group)
      !WRITE(*,*) "oi ",counter,tmp_index(:)
      states_index(counter,1:nstates_max) = tmp_index(1:nstates_max)
      !WRITE(*,*) "si ",counter,states_index(counter,:)
    END DO
  END DO
  !WRITE(*,*) "CAAAAAAAAAAA",counter
  DEALLOCATE(tmp_index,STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','tmp_index')

  !WRITE(*,*) "DAAAAAAAAAAA"
  DO imol=1,nmol
    !Count number of states for each molecule using the states_index(:,:) info
    counter = 0
    DO istates=1,nstates_max
      IF(states_index(imol,istates) /= 0) counter = counter + 1
    END DO
    nstates_array(imol) = counter
  END DO
  !WRITE(*,*) "EAAAAAAAAAAA"

#else

  !WRITE(*,*) "Enter serial section"
  nmol = SIZE(molecule_set)
  ALLOCATE(nstates_array(nmol),STAT=isos) !hold the nstates each imol has
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nstates_array')

  DO imol = 1,nmol
    nstates_array(imol) = molecule_set(imol)%lmi%nstates
  END DO

  nstates_max = MAXVAL(nstates_array)

  ALLOCATE(states_index(nmol,nstates_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','states_index')

  DO imol = 1,nmol
    states_index(imol,:) = molecule_set(imol)%lmi%states
  END DO

#endif

  IF (globenv%ionode) THEN 
    DO imol = 1,nmol
      WRITE(*,*) imol,states_index(imol,:)
    END DO
  END IF

  CALL cp_fm_get_info(mo_localized, &
                      ncol_global=ncol_global, &
                      nrow_global=nrow_global )
  NULLIFY(smo)
  CALL cp_fm_create(smo,mo_eigenvectors%matrix_struct)
  CALL cp_sm_fm_multiply(matrix_S,mo_eigenvectors,smo,ncol_global,qs_env%para_env)

  IF (globenv%ionode) THEN
     WRITE(unit_report,*) nmol,nstates_max, " number of molecules, maximum number of states "
  ENDIF

  DO imol = 1, nmol
    ns = nstates_array( imol )
    ind => states_index ( imol, : )
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

  END DO 
  IF (globenv%ionode) CLOSE(unit_report)
!------------------------------------------------------------------------------
END SUBROUTINE construct_molecular_states

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END MODULE molecular_states
