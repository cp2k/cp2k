!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_nmr_utils *
!!
!!   NAME
!!     qs_linres_nmr_utils
!!
!!   FUNCTION
!!     Chemical shift calculation by dfpt
!!     Initialization of the nmr_env, creation of the special neighbor lists
!!     Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!!     Write output
!!     Deallocate everything 
!!
!!   NOTE
!!     The psi0 should be localized
!!     the Sebastiani method works within the assumption that the orbitals are 
!!     completely contained in the simulation box 
!!  
!!   AUTHOR
!!      MI 
!!
!!   MODIFICATION HISTORY
!!      created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: restart_guess
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_methods,               ONLY: linres_localize, linres_scf
  USE qs_linres_nmr_op,                ONLY: set_vecp
  USE qs_linres_types,                 ONLY: linres_control_type,&
                                             nmr_env_create, nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_operators_ao,                 ONLY: set_up_op_rmd
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type, pbc
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,cp_sm_scale_and_add,&
                                             deallocate_matrix_set,&
                                             deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             transfer_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmr_do_nmr, nmr_env_cleanup, nmr_env_init

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_nmr_utils'

CONTAINS

  SUBROUTINE nmr_do_nmr(nmr_env,p_env,qs_env,globenv,error)

    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(global_environment_type), &
      POINTER                                   :: globenv
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_do_nmr', &
                                   routineP = moduleN//'/'//routineN

    INTEGER ::  handle, iao, ii, iii, idir, idir2, ispin, istat, ist_true,  istate, istate2, output_unit
    REAL(dp) :: dk(3), dkl(3), dl(3)
    REAL(dp), DIMENSION(:,:), POINTER           :: vecbuf_dklxp0
    TYPE(cp_2d_r_p_type), DIMENSION(3)          :: vecbuf_p0
    LOGICAL                                     :: ionode, failure
    TYPE(cell_type), POINTER                    :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: psi1,h1_psi0
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                   :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER               :: logger 
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(qs_matrix_pools_type), POINTER         :: mpools


    failure   = .FALSE.
    CALL timeset(routineN,"I"," ",handle)

    NULLIFY(ao_mo_fm_pools, cell, dft_control, linres_control)
    NULLIFY(logger, mpools, psi1,h1_psi0, vecbuf_dklxp0)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    IF(ionode .AND. logger%print_keys%linres) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "*** Self consistent optimization of the response wavefunctions ***"
    END IF

    CALL get_qs_env(qs_env=qs_env,&
          dft_control= dft_control,&
          mpools=mpools,cell=cell,&
          linres_control=linres_control)

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    CALL mpools_get(mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, psi1, name=routineP//":psi1")
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, h1_psi0, name=routineP//":h1_psi0")

    !Check if restart and from where
    IF(nmr_env%restart_nmr) THEN
!       CALL nmr_restart_nmr(nmr_env,qs_env,action="read",error=error)
    END IF

    ! operator p
    DO idir = 1,3
      IF(.NOT. nmr_env%simplenmr_done(idir)) THEN

         IF(ionode .AND. logger%print_keys%linres) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,I4/)")&
           "Response to the perturbation operator p in dir. ", idir
         END IF

        !Initial guess for psi1
        !Maybe it could be better than this
         DO ispin = 1,dft_control%nspins
           CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
         END DO 

         !DO scf cycle to optimize psi1
         DO ispin = 1,dft_control%nspins
           CALL cp_fm_to_fm(nmr_env%p_psi0(ispin,idir)%matrix, h1_psi0(ispin)%matrix)
         END DO
         linres_control%converged = .FALSE.
         CALL linres_scf(p_env, qs_env, psi1, h1_psi0,error=error)
 stop 'wait'

         IF(linres_control%converged) THEN
         !Calculate the contribution to the chemical shift
!         IF(ionode .AND. logger%print_keys%linres) THEN
!           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
!           "Calculate the contribution to the Chemical Shift tensor, Chi, and J matrixes "
!         END IF
!        CALL chemshift_p(nmr_env%chemical_shift,psi1,idir,error=error)

           IF(nmr_env%current_density) THEN
            ! the optimized wfns are copied in the  fm
             IF(ionode .AND. logger%print_keys%linres) THEN
               WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
               "Store the psi1 for the calculation of the response current density "
             END IF
             DO ispin=1, dft_control%nspins
               CALL cp_fm_to_fm(psi1(ispin)%matrix, nmr_env%psi1_p(ispin,idir)%matrix)
             END DO
           END IF
         END IF

         IF(ionode .AND. logger%print_keys%linres) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
            "Write the resulting psi1 in restart file "
         END IF
         ! Write the result inthe restart file
         nmr_env%simplenmr_done(idir) = .TRUE.
         nmr_env%simplenmr_converged(idir) = linres_control%converged
!        CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir,error=error)

      END IF  
    END DO ! idir

    ! operator rxp
    DO idir = 1,3 

      IF(.NOT. nmr_env%simplenmr_done(idir+3)) THEN

         IF(ionode .AND. logger%print_keys%linres) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,I4/)")&
           "Response to the perturbation operator rxp in dir. ", idir
         END IF

         !Initial guess for psi1
         !Maybe it could be beter than this
         DO ispin = 1,dft_control%nspins
           CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
         END DO 

         !DO scf cycle to optimize psi1
         DO ispin = 1,dft_control%nspins
            CALL cp_fm_to_fm(nmr_env%rxp_psi0(ispin,idir)%matrix, h1_psi0(ispin)%matrix)
         END DO

         ! Optimize the response wavefunctions 
         linres_control%converged = .FALSE.
         CALL linres_scf(p_env, qs_env, psi1, h1_psi0,error=error)

         IF(linres_control%converged) THEN
         !Calculate the contribution to the chemical shift
!         IF(ionode .AND. logger%print_keys%linres) THEN
!           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
!           "Calculate the contribution to the Chemical Shift tensor, Chi, and J matrixes "
!         END IF
!         CALL chemshift_p(nmr_env%chemical_shift,psi1,idir,error=error)

           IF(nmr_env%current_density) THEN
           !   the optimized wfns are copied in the  fm
             IF(ionode .AND. logger%print_keys%linres) THEN
               WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
               "Store the psi1 for the calculation of the response current density "
             END IF
             DO ispin=1, dft_control%nspins
               CALL cp_fm_to_fm(psi1(ispin)%matrix, nmr_env%psi1_rxp(ispin,idir)%matrix)
             END DO
           END IF
         END IF

         IF(ionode .AND. logger%print_keys%linres) THEN
           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
            "Write the resulting psi1 in restart file "
         END IF
         ! Write the result inthe restart file
         nmr_env%simplenmr_done(idir+3) = .TRUE.
!        CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir+3,error=error)
        END IF
     END DO ! idir

     IF(nmr_env%full_nmr) THEN
       DO ispin = 1,dft_control%nspins
         CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
       END DO 

!      The correction is state depedent a loop over the states is necessary
       DO idir = 1,3 
         NULLIFY(vecbuf_p0(idir)%array,vecbuf_dklxp0)
         ALLOCATE(vecbuf_p0(idir)%array(1,nmr_env%nao),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         vecbuf_p0(idir)%array(1,nmr_env%nao) = 0.0_dp
       END DO
       ALLOCATE(vecbuf_dklxp0(1,nmr_env%nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       vecbuf_dklxp0(1,nmr_env%nao) = 0.0_dp

       DO ispin=1, dft_control%nspins
         DO idir = 1,3
           DO istate = 1,nmr_env%nstates(ispin)
              ist_true = nmr_env%statetrueindex(idir,istate)
              !the initial guess is the previous set of psi1, just optimized
              CALL set_vecp(idir,ii,iii)
              dk(1:3) = nmr_env%centers_set(ispin)%array(1:3,ist_true)
              
              DO istate2 = 1,nmr_env%nstates(ispin)
                dl(1:3)  = nmr_env%centers_set(ispin)%array(1:3,istate2)
                dkl = pbc(dl,dk,cell)

                DO idir2 = 1,3
                  CALL cp_fm_get_submatrix(nmr_env%p_psi0(ispin,idir2)%matrix,&
                       vecbuf_p0(idir2)%array,&
                       1,istate2,nmr_env%nao,1,transpose=.TRUE.,error=error)
                END DO

                DO iao = 1,nmr_env%nao
                  vecbuf_dklxp0(1,iao) = dkl(ii)* vecbuf_p0(iii)%array(1,iao)-&
                                         dkl(iii)* vecbuf_p0(ii)%array(1,iao)
                END DO  ! iao
                CALL cp_fm_set_submatrix(h1_psi0(ispin)%matrix, vecbuf_dklxp0,&
                     1,istate2,nmr_env%nao,1,transpose=.TRUE.,error=error)

              END DO ! istate2 

              ! Optimize the response wavefunctions 
              CALL linres_scf(p_env, qs_env, psi1, h1_psi0,error=error)

         !Calculate the contribution to the chemical shift
!         IF(ionode .AND. logger%print_keys%linres) THEN
!           WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
!           "Calculate the contribution to the Chemical Shift tensor, Chi, and J matrixes "
!         END IF
!         CALL chemshift_p(nmr_env%chemical_shift,psi1,idir,error=error)


             IF(nmr_env%current_density) THEN
               ! the optimized wfns are copied in the  fm
               CALL cp_fm_get_submatrix(psi1(ispin)%matrix,vecbuf_dklxp0 ,&
                    1,ist_true,nmr_env%nao,1,transpose=.TRUE.,error=error)
               CALL cp_fm_set_submatrix(nmr_env%psi1_D(ispin,idir)%matrix, vecbuf_dklxp0,&
                    1,ist_true,nmr_env%nao,1,transpose=.TRUE.,error=error)
             END IF
             nmr_env%fullnmr_done(idir*ispin,istate) = .TRUE.

          END DO  ! istate

!        CALL nmr_restart_nmr(nmr_env,qs_env,action="write",ivec=idir+6,error=error)
          
        END DO  ! idir
       END DO  ! ispin
     END IF

    !Calculate current density and induced field
     IF(nmr_env%current_density) THEN
!      CALL current_density(nmr_env,error=error)
     END IF

    !----------!
    ! clean up !
    !----------!
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, psi1)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, h1_psi0)
     NULLIFY(psi1,h1_psi0)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE nmr_do_nmr

! *****************************************************************************
!!****f* cp2k/qs_linres_nmr_utils/nmr_env_init
!!
!!   NAME
!!     nmr_env_init
!!
!!   FUNCTION
!!     Initialize the nmr environment
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE nmr_env_init(nmr_env,qs_env,globenv,error)

    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(global_environment_type), &
      POINTER                                   :: globenv
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_env_init', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER                                     :: homo, icount,idir,ini, ispin,istat,&
                                                   istate,istate2,istate_next,j,k,&
                                                   n_mo(2),nao,natom,nmoloc,&
                                                   nstates, output_unit
    LOGICAL                                     :: ionode, failure
    LOGICAL, DIMENSION(:,:), ALLOCATABLE          :: state_done
    REAL(dp) :: center(3), center2(3), dist, mdist, rab(3)
    REAL(dp), DIMENSION(:,:), POINTER           :: centers
    TYPE(atomic_kind_type), DIMENSION(:),&
      POINTER                                   :: atomic_kind_set
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                   :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:),           POINTER           :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER  :: particle_set
    TYPE(qs_matrix_pools_type), POINTER         :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s
    TYPE(scf_control_type), POINTER             :: scf_control
    TYPE(section_vals_type), POINTER            :: nmr_section
    TYPE(cell_type), POINTER                    :: cell

    failure = .FALSE.

    NULLIFY(atomic_kind_set,cell,centers,dft_control,linres_control,scf_control)
    NULLIFY(logger,matrix_s,mos,mpools,nmr_section,particle_set,sab_orb)
    n_mo(1:2) = 0
    nao = 0
    nmoloc = 0

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    IF(nmr_env%ref_count /= 0) THEN
      CALL nmr_env_cleanup(nmr_env,error=error)
    END IF


    IF(ionode .AND. logger%print_keys%linres) THEN

      WRITE (UNIT=output_unit,FMT="(/,T20,A,/)")&
         "*** Start NMR Chemical Shift Calculation ***"
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "Inizialization of the NMR environment"

    END IF

    CALL nmr_env_create(nmr_env,error=error)

    ! If surrent_density or full_nmr different allocations are required
    nmr_section => section_vals_get_subs_vals(globenv%input_file, &
                                                 "FORCE_EVAL%PROPERTIES%LINRES%NMR")
    CALL section_vals_val_get(nmr_section,"CURRENT_DENSITY", &
         l_val=nmr_env%current_density,error=error)
    CALL section_vals_val_get(nmr_section,"FULL",l_val=nmr_env%full_nmr,error=error)
    CALL section_vals_val_get(nmr_section,"RESTART_NMR",l_val=nmr_env%restart_nmr,error=error)

    ! Check if restat also psi0 shoul be restarted
    IF(nmr_env%restart_nmr .AND. scf_control%density_guess/=restart_guess)THEN
      CALL stop_program("nmr_env_init","restart_nmr requires density_guess=restart")
    END IF


    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    matrix_s=matrix_s,&
                    mos=mos,&
                    mpools=mpools,&
                    particle_set=particle_set,&
                    scf_control=scf_control)

    ! check that the psi0 are localized and you have all the centers
     CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
     IF(failure .AND. ionode) THEN
       WRITE (UNIT=output_unit,FMT='(A)') &
          ' To get NMR parameters within PBC you need localized zero order orbitals '
     END IF
     ALLOCATE(nmr_env%centers_set(dft_control%nspins),STAT=istat)
     nstates = 0
     nmr_env%nstates = 0
     DO ispin = 1,dft_control%nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=n_mo(ispin),homo=homo)
!       centers => linres_control%localized_wfn_control%centers_set(ispin)%array
       nmoloc = size(linres_control%localized_wfn_control%centers_set(ispin)%array,2)
       CPPrecondition(nmoloc==homo,cp_warning_level,routineP,error,failure)
       ALLOCATE(nmr_env%centers_set(ispin)%array(3,nmoloc),STAT=istat)
       CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
       ! point to the psi0 centers
       DO idir = 1,3
         nmr_env%centers_set(ispin)%array(idir,:) = &
           linres_control%localized_wfn_control%centers_set(ispin)%array(idir,:)
       end do
       nmr_env%nstates(ispin) = homo
       nstates = MAX(nstates,homo)
     END DO  
     nmr_env%nao =nao
     ALLOCATE(nmr_env%statetrueindex(3,nstates),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     nmr_env%statetrueindex = 0
     ALLOCATE(state_done(3,nstates),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     state_done = .false.

     nmr_env%statetrueindex(1,1) = 1
     center(1)  = nmr_env%centers_set(1)%array(1,1)
     center(2)  = nmr_env%centers_set(1)%array(2,1)
     center(3)  = nmr_env%centers_set(1)%array(3,1)
     state_done(1,1) = .true. 
     icount = 1

     DO idir = 1,3
       ini = 1
       if(idir == 1) ini = 2
       DO istate = ini, nstates
         mdist = 100.0_dp

         DO istate2 = 1,nstates
           IF(.NOT.state_done(idir,istate2)) THEN
             center2(1)  = nmr_env%centers_set(1)%array(1,istate2)
             center2(2)  = nmr_env%centers_set(1)%array(2,istate2)
             center2(3)  = nmr_env%centers_set(1)%array(3,istate2)

             rab = pbc(center,center2,cell)
             CALL set_vecp(idir,j,k)
             dist = sqrt(rab(j)*rab(j)+rab(k)*rab(k))

             IF(dist .LT. mdist) THEN
                mdist =dist
                istate_next = istate2
             END IF
           END IF
         END DO  ! istate2

         icount = icount + 1
         state_done(idir,istate_next) = .true. 
         nmr_env%statetrueindex(idir,icount) = istate_next

         center(1)  = nmr_env%centers_set(1)%array(1,istate_next)
         center(2)  = nmr_env%centers_set(1)%array(2,istate_next)
         center(3)  = nmr_env%centers_set(1)%array(3,istate_next)
       END DO  ! istate
       icount = 0
     END DO  ! idir
     DEALLOCATE(state_done,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=nmr_env%ao_mo_fm_pools, error=error)
    ! for the chemical shift we need 6 psi1, i.e. 6 optimization procedures
    ! They become 9 if full nmr is calculated, i.e. with the correction term too
    ! All of them are required at the end of the optimization procedure
    ! if the current density and the induced fiels have to be calculated 
    ! If instead only the shift is needed, only one psi1 should be enough, providing
    ! that after every optimization the corresponding shift contribution is calculated
    ! prepare the psi1
     IF(nmr_env%current_density) THEN

       ALLOCATE(nmr_env%psi1_p(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(nmr_env%psi1_rxp(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,dft_control%nspins
         DO idir = 1,3
           NULLIFY(nmr_env%psi1_p(ispin,idir)%matrix)
           CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%psi1_p(ispin,idir)%matrix)
           NULLIFY(nmr_env%psi1_rxp(ispin,idir)%matrix)
           CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%psi1_rxp(ispin,idir)%matrix)
         END DO 
       END DO 

       IF(nmr_env%full_nmr) THEN
         ALLOCATE(nmr_env%psi1_D(dft_control%nspins,3),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO ispin = 1,dft_control%nspins
           DO idir = 1,3
             NULLIFY(nmr_env%psi1_D(ispin,idir)%matrix)
             CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%psi1_D(ispin,idir)%matrix)
           END DO
         END DO
       END IF
     END IF

     ! px py pz operator acting on the atomic orbitals
     CALL allocate_matrix_set(nmr_env%op_p_ao,3)
     DO idir = 1,3 
       CALL replicate_matrix_structure(matrix_s(1)%matrix, &
            nmr_env%op_p_ao(idir)%matrix,"nmr_env%op_p_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),&
            target_symmetry="symmetric")
       CALL set_matrix(nmr_env%op_p_ao(idir)%matrix,0.0_dp)
     END DO

     ! for the rxp we cannot calculate it a priori because it is in facts (r-dk)xp 
     ! where dk is the center of the orbital it is applied to. We would need nstate operators
     ! What we can store is (r-dk)xp|psi0k> for each k, which is a full matrix only
     ! Therefore we prepare here the full matrix p_psi0 and rxp_psi0
     ! We also need a temporary sparse matrix where to store the integrals during the calculation
     ALLOCATE(nmr_env%p_psi0(dft_control%nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(nmr_env%rxp_psi0(dft_control%nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(nmr_env%p_psi1(dft_control%nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO ispin = 1,dft_control%nspins
       DO idir = 1,3
         NULLIFY(nmr_env%p_psi0(ispin,idir)%matrix)
         CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%p_psi0(ispin,idir)%matrix)
         NULLIFY(nmr_env%rxp_psi0(ispin,idir)%matrix)
         CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%rxp_psi0(ispin,idir)%matrix)
         NULLIFY(nmr_env%p_psi1(ispin,idir)%matrix)
         CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%p_psi1(ispin,idir)%matrix)
       END DO 
     END DO 

     CALL allocate_matrix_set(nmr_env%op_rmd_ao,3)

     CALL set_up_op_rmd(nmr_env%op_rmd_ao(1)%matrix,qs_env,error=error)

     DO idir = 2,3 
       CALL replicate_matrix_structure( nmr_env%op_rmd_ao(1)%matrix, &
            nmr_env%op_rmd_ao(idir)%matrix,"nmr_env%op_rmd_ao"//&
            "-"//TRIM(ADJUSTL(cp_to_string(idir))),target_symmetry="none")
       CALL set_matrix(nmr_env%op_rmd_ao(idir)%matrix,0.0_dp)
     END DO


    IF(nmr_env%current_density) THEN
    ! Current density matrix in x y and z:
       ALLOCATE( nmr_env%jp1_Bx_ao(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE( nmr_env%jp1_By_ao(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE( nmr_env%jp1_Bz_ao(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DO idir = 1,3 
         DO ispin = 1,dft_control%nspins
           NULLIFY(nmr_env%jp1_Bx_ao(ispin,idir)%matrix)
           NULLIFY(nmr_env%jp1_By_ao(ispin,idir)%matrix)
           NULLIFY(nmr_env%jp1_Bz_ao(ispin,idir)%matrix)
           CALL replicate_matrix_structure( nmr_env%op_rmd_ao(1)%matrix, &
              nmr_env%jp1_Bx_ao(ispin,idir)%matrix,"nmr_env%jp1_Bx_ao"//&
              "-"//TRIM(ADJUSTL(cp_to_string(idir))),target_symmetry="none")
           CALL set_matrix(nmr_env%jp1_Bx_ao(ispin,idir)%matrix,0.0_dp)
           CALL replicate_matrix_structure( nmr_env%op_rmd_ao(1)%matrix, &
              nmr_env%jp1_By_ao(ispin,idir)%matrix,"nmr_env%jp1_Bx_ao"//&
              "-"//TRIM(ADJUSTL(cp_to_string(idir))),target_symmetry="none")
           CALL set_matrix(nmr_env%jp1_By_ao(ispin,idir)%matrix,0.0_dp)
           CALL replicate_matrix_structure( nmr_env%op_rmd_ao(1)%matrix, &
              nmr_env%jp1_Bz_ao(ispin,idir)%matrix,"nmr_env%jp1_Bz_ao"//&
              "-"//TRIM(ADJUSTL(cp_to_string(idir))),target_symmetry="none")
           CALL set_matrix(nmr_env%jp1_Bz_ao(ispin,idir)%matrix,0.0_dp)
         END DO
       END DO
     END IF

 
    ! Initialize the chemical shift tensor
    nmr_env%chi_tensor = 0.0_dp
    natom = SIZE(particle_set,1)
    ALLOCATE(nmr_env%chemical_shift(3,3,natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    nmr_env%chemical_shift = 0.0_dp
    nmr_env%simplenmr_done(1:6) = .FALSE.
    ALLOCATE(nmr_env%fullnmr_done(3*dft_control%nspins,nstates),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    nmr_env%fullnmr_done = .FALSE.

  END SUBROUTINE nmr_env_init

! *****************************************************************************
!!****f* cp2k/qs_linres_nmr_utils/nmr_env_cleanup
!!
!!   NAME
!!     nmr_env_cleanup
!!
!!   FUNCTION
!!     Deallocate the nmr environment
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE nmr_env_cleanup(nmr_env,error)

    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_env_cleanup', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER                                     :: idir,ispin,istat
    LOGICAL                                     :: failure

    failure=.FALSE.
    IF(.NOT. failure) THEN
      nmr_env%ref_count = nmr_env%ref_count - 1
      IF(nmr_env%ref_count == 0 ) THEN
        !psi1_p
        IF(ASSOCIATED(nmr_env%psi1_p)) THEN
          DO idir = 1,SIZE(nmr_env%psi1_p,2)
            DO ispin = 1,SIZE(nmr_env%psi1_p,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%psi1_p(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%psi1_p, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !psi1_rxp
        IF(ASSOCIATED(nmr_env%psi1_rxp)) THEN
          DO idir = 1,SIZE(nmr_env%psi1_rxp,2)
            DO ispin = 1,SIZE(nmr_env%psi1_rxp,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%psi1_rxp(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%psi1_rxp, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !psi1_D
        IF(ASSOCIATED(nmr_env%psi1_D)) THEN
          DO idir = 1,SIZE(nmr_env%psi1_D,2)
            DO ispin = 1,SIZE(nmr_env%psi1_D,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%psi1_D(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%psi1_D, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !p_psi0
        IF(ASSOCIATED(nmr_env%p_psi0)) THEN
          DO idir = 1,SIZE(nmr_env%p_psi0,2)
            DO ispin = 1,SIZE(nmr_env%p_psi0,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%p_psi0(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%p_psi0, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !rxp_psi0
        IF(ASSOCIATED(nmr_env%rxp_psi0)) THEN
          DO idir = 1,SIZE(nmr_env%rxp_psi0,2)
            DO ispin = 1,SIZE(nmr_env%rxp_psi0,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%rxp_psi0(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%rxp_psi0, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
       !p_psi1
        IF(ASSOCIATED(nmr_env%p_psi1)) THEN
          DO idir = 1,SIZE(nmr_env%p_psi1,2)
            DO ispin = 1,SIZE(nmr_env%p_psi1,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%p_psi1(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%p_psi1, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        NULLIFY(nmr_env%ao_mo_fm_pools)

        DO ispin = 1,SIZE(nmr_env%centers_set,1)
          DEALLOCATE(nmr_env%centers_set(ispin)%array,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO 
        DEALLOCATE(nmr_env%centers_set,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        !op_p_ao
        IF(ASSOCIATED(nmr_env%op_p_ao)) THEN
          CALL deallocate_matrix_set(nmr_env%op_p_ao)
        END IF
        !op_rmd_ao
        IF(ASSOCIATED(nmr_env%op_rmd_ao)) THEN
          CALL deallocate_matrix_set(nmr_env%op_rmd_ao)
        END IF
        !jp1_Bx_ao jp1_By_ao jp1_Bz_ao
        IF(ASSOCIATED(nmr_env%jp1_Bx_ao)) THEN
          DO idir = 1,3
            DO ispin = 1,SIZE(nmr_env%jp1_Bx_ao,1)
              CALL deallocate_matrix(nmr_env%jp1_Bx_ao(ispin,idir)%matrix)
              CALL deallocate_matrix(nmr_env%jp1_By_ao(ispin,idir)%matrix)
              CALL deallocate_matrix(nmr_env%jp1_Bz_ao(ispin,idir)%matrix)
            END DO
          END DO
          DEALLOCATE(nmr_env%jp1_Bx_ao,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(nmr_env%jp1_By_ao,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(nmr_env%jp1_Bz_ao,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !fullnmr_done
        IF(ASSOCIATED(nmr_env%fullnmr_done)) THEN
          DEALLOCATE(nmr_env%fullnmr_done,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !chemical_shift
        IF(ASSOCIATED(nmr_env%chemical_shift)) THEN
          DEALLOCATE(nmr_env%chemical_shift,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(nmr_env%statetrueindex)) THEN
          DEALLOCATE(nmr_env%statetrueindex,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
      END IF
    END IF
    
  END SUBROUTINE nmr_env_cleanup


END MODULE  qs_linres_nmr_utils
