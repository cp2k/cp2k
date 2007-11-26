!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/vibrational_analysis [1.0] *
!!
!!   NAME
!!     vibrational_analysis
!!
!!   FUNCTION
!!     Module performing a vibrational analysis
!!
!!   NOTES
!!     Numerical accuracy for parallel runs:
!!      Each replica starts the SCF run from the one optimized
!!      in a previous run. It may happen then energies and derivatives
!!      of a serial run and a parallel run could be slightly different
!!      'cause of a different starting density matrix.
!!      Exact results are obtained using:
!!         EXTRAPOLATION USE_GUESS in QS section (Teo 08.2006)
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
MODULE vibrational_analysis
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_matrix
  USE physcon,                         ONLY: angstrom,&
                                             massunit,&
                                             vibfac
  USE replica_types,                   ONLY: rep_env_calc_e_f,&
                                             rep_env_create,&
                                             rep_env_release,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'vibrational_analysis'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: thrs = 100.0_dp*EPSILON(0.0_dp)
  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.

  PUBLIC :: vb_anal

CONTAINS

  !!****f* cp2k/vibrational_analysis [1.0] *
  !!
  !!   NAME
  !!     vibrational_analysis
  !!
  !!   FUNCTION
  !!     Module performing a vibrational analysis
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Teodoro Laino 08.2006
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE vb_anal(input, para_env, globenv, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vb_anal', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1), DIMENSION(3)           :: lab = (/"X","Y","Z"/)
    INTEGER :: handle, i, icoord, icoordm, icoordp, ierr, imap, ip1, ip2, &
      iparticle1, iparticle2, iseq, iw, j, mat_size, n_rep, natoms, ncoord, &
      nrep, nRotTrM, nvib, output_molden, output_unit, precondit, prep, stat
    INTEGER, DIMENSION(:), POINTER           :: Clist, Mlist
    LOGICAL                                  :: converged, do_mode_tracking, &
                                                failure, keep_rotations, &
                                                something_frozen
    REAL(KIND=dp)                            :: dx, norm, sel_freq, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: H_eigval1, H_eigval2, konst, &
                                                mass, pos0, rmass, step_b, &
                                                step_r
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: b_vec, delta_vec, Hessian, &
                                                Hint1, Hint2, ms_force
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eig_bfgs
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: b_mat, D, hes_bfgs, RotTrM, &
                                                s_mat
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                mode_tracking_section

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(D, RotTrM, logger, subsys, f_env, particles, rep_env, eig_bfgs, hes_bfgs, b_mat, s_mat)
       logger => cp_error_get_logger(error)
       output_unit=cp_print_key_unit_nr(logger,input,"VIBRATIONAL_ANALYSIS%PROGRAM_RUN_INFO",&
            extension=".vibLog",error=error)
       output_molden=cp_print_key_unit_nr(logger,input,"VIBRATIONAL_ANALYSIS%MOLDEN_VIB",&
            extension=".mol",file_status='REPLACE',error=error)
       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%DX",&
            r_val=dx, error=error)
       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%NPROC_REP",&
            i_val=prep, error=error)
       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%FULLY_PERIODIC",&
            l_val=keep_rotations, error=error)
       mode_tracking_section=>section_vals_get_subs_vals(input,"VIBRATIONAL_ANALYSIS%MODE_SELECTIVE",&
                              error=error)
       CALL section_vals_get(mode_tracking_section,explicit=do_mode_tracking, error=error)
       nrep = MAX(1,para_env%num_pe/prep)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"VIB|",3(A,I5))')"NUMBER OF REPLICAS:",nrep,&
               "NUMBER OF PROCS:",para_env%num_pe,"NUMBER OF PROCS PER REPLICA:",prep
       END IF
       ! Just one force_env allowed
       force_env_section => section_vals_get_subs_vals(input,"FORCE_EVAL",error=error)
       CALL section_vals_get(force_env_section,n_repetition=n_rep,error=error)
       CPPostcondition(n_rep==1,cp_failure_level,routineP,error,failure)       
       ! Create Replica Environments
       CALL rep_env_create(rep_env, para_env=para_env, input=input,&
            nrep=nrep,prep=prep,error=error)
       IF (ASSOCIATED(rep_env)) THEN
          CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
               new_error=new_error, failure=failure)
          CALL force_env_get(f_env%force_env,subsys=subsys,error=error)
          particles => subsys%particles%els
          ! Decide which kind of Vibrational Analysis to perform
          IF (do_mode_tracking)THEN
             converged=.FALSE.
             natoms=SIZE(particles)
             ncoord=3*natoms
             ALLOCATE(mass(3*natoms), stat=stat)            
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO i = 1, natoms
                DO j=1,3
                   mass((i-1)*3+j) = particles(i)%atomic_kind%mass
                   mass((i-1)*3+j)= SQRT(mass((i-1)*3+j))
                END DO
             END DO
             ! Allocate working arrays
             ALLOCATE(delta_vec(ncoord,nrep),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(b_vec(ncoord,nrep),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(step_r(nrep),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(step_b(nrep),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

             CALL MS_initial_moves( para_env, nrep , input, globenv,&
                                    particles,delta_vec,b_vec, step_r,&
                                    step_b,sel_freq ,mass,hes_bfgs,eig_bfgs,&
                                    precondit,dx,output_unit,&
                                    b_mat,s_mat,mat_size,error)
             ncoord=3*natoms
             ALLOCATE(pos0(ncoord), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(ms_force(ncoord,nrep),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO i = 1, natoms
                DO j=1,3
                   pos0((i-1)*3+j) = particles((i))%r(j)
                END DO
             END DO
             ncoord=3*natoms
             DO
                ms_force = HUGE(0.0_dp)
                DO i=1,nrep
                   DO j=1,ncoord
                      rep_env%r(j,i)=pos0(j)+step_r(i)*delta_vec(j,i)
                   END DO
                END DO
                CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)
                DO i = 1, nrep
                   DO j = 1, ncoord
                      ms_force(j,i) = rep_env%f(j,i)
                   END DO
                END DO
                DO i=1,nrep
                   DO j=1,ncoord
                      rep_env%r(j,i)=pos0(j)-step_r(i)*delta_vec(j,i)
                   END DO
                END DO
                CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)
                CALL evaluate_H_update_b(mat_size,ms_force,rep_env,b_vec,nrep,&
                                         particles,b_mat,s_mat,step_b,step_r,&
                                         delta_vec,mass,sel_freq,converged,hes_bfgs,&
                                         eig_bfgs,precondit,dx,output_unit,output_molden, error)
                IF(converged)EXIT
             END DO
             CALL f_env_rm_defaults(f_env,new_error,ierr)

             DEALLOCATE(ms_force, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(pos0, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(step_r, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(step_b, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(b_vec, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(delta_vec, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(mass, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(b_mat, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(s_mat, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF (ASSOCIATED(eig_bfgs)) THEN 
                DEALLOCATE(eig_bfgs, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(hes_bfgs)) THEN 
                DEALLOCATE(hes_bfgs, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          ELSE
             CALL get_moving_atoms(force_env=f_env%force_env,Ilist=Mlist,error=error)
             something_frozen = SIZE(particles).NE.SIZE(Mlist)
             natoms=SIZE(Mlist)
             ncoord=natoms*3
             ALLOCATE(Clist(ncoord), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
             ALLOCATE(mass(natoms), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(pos0(ncoord), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(Hessian(ncoord,ncoord),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
             Clist = 0
             DO i = 1, natoms
                imap             = Mlist(i)
                Clist((i-1)*3+1) = (imap-1)*3+1
                Clist((i-1)*3+2) = (imap-1)*3+2
                Clist((i-1)*3+3) = (imap-1)*3+3
                mass(i) = particles(imap)%atomic_kind%mass
                CPPostcondition(mass(i)>0.0_dp,cp_failure_level,routineP,error,failure)
                mass(i)= SQRT(mass(i))
                pos0((i-1)*3+1) = particles(imap)%r(1)
                pos0((i-1)*3+2) = particles(imap)%r(2)
                pos0((i-1)*3+3) = particles(imap)%r(3)
             END DO
             ! 
             ! Determine the principal axes of inertia.
             ! Generation of coordinates in the rotating and translating frame
             !
             IF (something_frozen) THEN
                nRotTrM = 0
                ALLOCATE(RotTrM(natoms*3,nRotTrM),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ELSE
                CALL rot_ana(particles,RotTrM,nRotTrM,output_unit,&
                     keep_rotations,natoms,error)
             END IF
             ! Generate the suitable rototranslating basis set
             CALL build_D_matrix(particles,RotTrM,nRotTrM,D,full=.FALSE.,&
                  natoms=natoms,error=error)
             !
             ! Loop on atoms and coordinates
             !
             Hessian = HUGE(0.0_dp)
             IF (output_unit>0) WRITE(output_unit,'(/,T2,A)')"VIB| Vibrational Analysis Info"
             DO icoordp = 1, ncoord, nrep
                icoord  = icoordp-1
                DO j = 1, nrep
                   DO i = 1, ncoord
                      imap = Clist(i)
                      rep_env%r(imap,j)=pos0(i)
                   END DO
                   IF (icoord+j <= ncoord) THEN
                      imap = Clist(icoord+j)
                      rep_env%r(imap,j) = rep_env%r(imap,j) + Dx
                   END IF
                END DO
                CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)
                DO j = 1, nrep
                   IF (icoord+j <= ncoord) THEN
                      DO i = 1, ncoord
                         imap = Clist(i)
                         Hessian(i,icoord+j) = rep_env%f(imap,j)
                      END DO
                      imap = Clist(icoord+j)
                      ! Dump Info
                      IF (output_unit>0) THEN
                         iparticle1 = imap/3
                         IF (MOD(imap,3)/=0) iparticle1=iparticle1+1
                         WRITE(output_unit,'(T2,A,I5,A,I5,3A)')&
                              "VIB| REPLICA Nr.",j,"- Energy and Forces for particle:",&
                              iparticle1,"  coordinate: ",lab(imap-(iparticle1-1)*3),&
                              " + D"//TRIM(lab(imap-(iparticle1-1)*3))
                         !
                         WRITE(output_unit,'(T2,A,T42,A,9X,F15.9)')&
                              "VIB|"," Total Energy: ",rep_env%f(rep_env%ndim+1,j)
                         WRITE(output_unit,'(T2,"VIB|",T10,"ATOM",T33,3(9X,A,7X))')lab(1),lab(2),lab(3)
                         DO i = 1, natoms
                            imap = Mlist(i)
                            WRITE(output_unit,'(T2,"VIB|",T12,A,T30,3(2X,F15.9))')&
                                 particles(imap)%atomic_kind%name,&
                                 rep_env%f((imap-1)*3+1:(imap-1)*3+3,j)
                         END DO
                      END IF
                   END IF
                END DO
             END DO
             DO icoordm = 1, ncoord, nrep
                icoord  = icoordm-1
                DO j = 1, nrep
                   DO i = 1, ncoord
                      imap = Clist(i)
                      rep_env%r(imap,j)=pos0(i)
                   END DO
                   IF (icoord+j <= ncoord) THEN
                      imap = Clist(icoord+j)
                      rep_env%r(imap,j) = rep_env%r(imap,j) - Dx
                   END IF
                END DO
                CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)
                DO j = 1, nrep
                   IF (icoord+j <= ncoord) THEN
                      imap = Clist(icoord+j)
                      iparticle1 = imap/3
                      IF (MOD(imap,3)/=0) iparticle1=iparticle1+1
                      ip1  = (icoord+j)/3
                      IF (MOD(icoord+j,3)/=0) ip1=ip1+1
                      ! Dump Info
                      IF (output_unit>0) THEN
                         WRITE(output_unit,'(T2,A,I5,A,I5,3A)')&
                              "VIB| REPLICA Nr.",j,"- Energy and Forces for particle:",&
                              iparticle1,"  coordinate: ",lab(imap-(iparticle1-1)*3),&
                              " - D"//TRIM(lab(imap-(iparticle1-1)*3))
                         !
                         WRITE(output_unit,'(T2,A,T42,A,9X,F15.9)')&
                              "VIB|"," Total Energy: ",rep_env%f(rep_env%ndim+1,j)
                         WRITE(output_unit,'(T2,"VIB|",T10,"ATOM",T33,3(9X,A,7X))')lab(1),lab(2),lab(3)
                         DO i = 1, natoms
                            imap = Mlist(i)
                            WRITE(output_unit,'(T2,"VIB|",T12,A,T30,3(2X,F15.9))')&
                                 particles(imap)%atomic_kind%name,&
                                 rep_env%f((imap-1)*3+1:(imap-1)*3+3,j)
                         END DO
                      END IF
                      DO iseq = 1, ncoord
                         imap = Clist(iseq)
                         iparticle2 = imap/3
                         IF (MOD(imap,3)/=0) iparticle2=iparticle2+1
                         ip2  = iseq/3
                         IF (MOD(iseq,3)/=0) ip2=ip2+1
                         tmp  = Hessian(iseq,icoord+j) - rep_env%f(imap,j)
                         tmp  = - tmp/(2.0_dp*Dx*mass(ip1)*mass(ip2))*1E6_dp
                         ! Mass weighted Hessian
                         Hessian(iseq,icoord+j) = tmp
                      END DO
                   END IF
                END DO
             END DO
             ! Dump Info
             IF (output_unit>0) THEN
                WRITE(output_unit,'(T2,A)')"VIB| Hessian in cartesian coordinates"
                CALL write_particle_matrix(Hessian,particles,output_unit,el_per_part=3,&
                     Ilist=Mlist,error=error)
             END IF
             ! Enforce symmetry in the Hessian
             DO i =1, ncoord
                DO j = i, ncoord
                   ! Take the upper diagonal part
                   Hessian(j,i) = Hessian(i,j)
                END DO
             END DO
             nvib = ncoord-nRotTrM
             ALLOCATE(H_eigval1(ncoord),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(H_eigval2(SIZE(D,2)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(Hint1(ncoord,ncoord),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(Hint2(SIZE(D,2),SIZE(D,2)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(rmass(SIZE(D,2)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(konst(SIZE(D,2)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             Hint1 = Hessian
             CALL diamat_all(Hint1,H_eigval1)
             IF (output_unit>0) THEN
                WRITE(output_unit,'(T2,"VIB| Cartesian Low frequencies ---",4G12.5)')&
                     (H_eigval1(i),i=1,MIN(9,ncoord))
             END IF
             IF (nvib/=0) THEN
                Hint2=MATMUL(TRANSPOSE(D),MATMUL(Hessian,D))
                CALL diamat_all(Hint2,H_eigval2)
                IF (output_unit>0) THEN
                   WRITE(output_unit,'(T2,"VIB| Frequencies after removal of the rotations and translations")')
                   ! Frequency at the moment are in a.u.
                   WRITE(output_unit,'(T2,"VIB| Internal  Low frequencies ---",4G12.5)')H_eigval2
                END IF
                Hessian = 0.0_dp
                DO i = 1, natoms
                   DO j = 1, 3
                      Hessian((i-1)*3+j,(i-1)*3+j) = 1.0_dp/mass(i)
                   END DO
                END DO
                ! Cartesian displacements of the normal modes
                D = MATMUL(Hessian,MATMUL(D,Hint2))
                DO i = 1, nvib
                   norm     = 1.0_dp/SUM(D(:,i)*D(:,i))
                   ! Reduced Masess
                   rmass(i) = norm/massunit
                   ! Renormalize displacements and convert in Angstrom
                   D(:,i)   = SQRT(norm)*D(:,i)             
                   ! Force constants
                   konst(i) = SIGN(1.0_dp,H_eigval2(i))*2.0_dp*pi**2*(ABS(H_eigval2(i))/massunit)**2*rmass(i)
                   ! Convert frequencies to cm^-1
                   H_eigval2(i) =  SIGN(1.0_dp,H_eigval2(i))*SQRT(ABS(H_eigval2(i))*massunit)*vibfac/1000.0_dp
                END DO
                ! Dump Info
                IF (para_env%mepos==para_env%source) THEN             
                   iw = cp_logger_get_default_unit_nr(logger)
                   CALL vib_out(iw,nvib,D,konst,rmass,H_eigval2,particles,Mlist)
                END IF
                IF(.NOT.something_frozen) CALL molden_out(output_molden,particles,H_eigval2,D,error)
             ELSE
                IF (output_unit>0) THEN
                   WRITE(output_unit,'(T2,"VIB| No further vibrational info. Detected a single atom")')
                END IF
             END IF
             ! Deallocate working arrays
             DEALLOCATE(Clist,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(Mlist,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(H_eigval1,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(H_eigval2,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(Hint1,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(Hint2,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(rmass,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(konst,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(mass, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(pos0, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(D,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)                    
             DEALLOCATE(Hessian,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
             CALL f_env_rm_defaults(f_env,new_error,ierr)
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,input,"VIBRATIONAL_ANALYSIS%PROGRAM_RUN_INFO",error=error)
       CALL cp_print_key_finished_output(output_molden,logger,input,"VIBRATIONAL_ANALYSIS%MOLDEN_VIB",error=error)
       CALL rep_env_release(rep_env,error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE vb_anal

  !!****f* cp2k/get_moving_atoms [1.0] *
  !!
  !!   NAME
  !!     get_moving_atoms
  !!
  !!   FUNCTION
  !!     give back a list of moving atoms
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Teodoro Laino 08.2006
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE get_moving_atoms(force_env, Ilist, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, DIMENSION(:), POINTER           :: Ilist
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_moving_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ii, ikind, iparticle, iparticle_local, ndim, &
      nfixed_atoms, nkind, nparticle_local, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: fixme
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    IF (.NOT. failure) THEN
       CALL force_env_get(force_env=force_env, subsys=subsys, error=error)

       CALL cp_subsys_get(subsys=subsys,particles=particles,local_particles=local_particles,&
            molecule_kinds_new=molecule_kinds, error=error)

       nkind             =  molecule_kinds % n_els
       molecule_kind_set => molecule_kinds % els
       particle_set      => particles%els
       ALLOCATE(fixme(SIZE(particle_set)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
       fixme=0

       Mol_kind:  DO ikind = 1, nkind
          molecule_kind => molecule_kind_set(ikind)
          CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list )

          IF (nfixed_atoms ==0) CYCLE Mol_Kind
          DO i = 1, SIZE(local_particles%n_el)
             nparticle_local   = local_particles%n_el(i)
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(i)%array(iparticle_local)
                DO ii = 1, SIZE(fixd_list)
                   IF (fixd_list(ii)%fixd == iparticle) THEN
                      IF (.NOT.fixd_list(ii)%restraint%active) THEN
                         ! Identify fixed atoms..
                         fixme(iparticle)=1
                      END IF
                      EXIT
                   END IF
                END DO
             END DO
          END DO
       END DO Mol_kind

       CALL mp_sum(fixme,force_env%para_env%group)
       ndim = COUNT(fixme==0)
       ALLOCATE(Ilist(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
       i = 0
       DO iparticle=1,SIZE(particle_set)
          IF (fixme(iparticle)==0) THEN
             i = i + 1
             Ilist(i) = iparticle
          END IF
       ENDDO
       DEALLOCATE(fixme,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE get_moving_atoms

  !!****f* cp2k/vib_out [1.0] *
  !!
  !!   NAME
  !!     vib_out
  !!
  !!   FUNCTION
  !!     Dumps results of the vibrational analysis
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Teodoro Laino 08.2006
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE vib_out(iw,nvib,D,k,m,freq,particles,Mlist)
    INTEGER, INTENT(IN)                      :: iw, nvib
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: D
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: k, m, freq
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER, DIMENSION(:), POINTER           :: Mlist

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: from, iatom, icol, j, jatom, &
                                                katom, natom, to

    natom = SIZE(D,1)
    WRITE(UNIT=iw,FMT="(/,T2,'VIB|',T30,'NORMAL MODES - CARTESIAN DISPLACEMENTS')")
    WRITE(UNIT=iw,FMT="(T2,'VIB|')")
    DO jatom=1,nvib,3
       from = jatom
       to = MIN(from+2,nvib)
       WRITE(UNIT=iw,FMT="(T2,'VIB|',13X,3(8X,I5,8X))")&
            (icol,icol=from,to)
       WRITE(UNIT=iw,FMT="(T2,'VIB|Frequency (cm^-1)',3(1X,F12.6,8X))")&
            (freq(icol),icol=from,to)
       WRITE(UNIT=iw,FMT="(T2,'VIB|Red.Masses (a.u.)',3(1X,F12.6,8X))")&
            (m(icol),icol=from,to)
       WRITE(UNIT=iw,FMT="(T2,'VIB|Frc consts (a.u.)',3(1X,F12.6,8X))")&
            (k(icol),icol=from,to)
       WRITE(UNIT=iw,FMT="(T2,' ATOM',2X,'EL',7X,3(4X,'  X  ',1X,'  Y  ',1X,'  Z  '))")
       DO iatom=1,natom,3
          katom = iatom/3
          IF (MOD(iatom,3)/=0) katom = katom+1
          CALL get_atomic_kind(atomic_kind=particles(Mlist(katom))%atomic_kind,&
               element_symbol=element_symbol)
          WRITE (UNIT=iw,FMT="(T2,I5,2X,A2,7X,3(4X,2(F5.2,1X),F5.2))")&
               Mlist(katom),element_symbol,&
               ((D(iatom+j,icol),j=0,2),icol=from,to)
       END DO
       WRITE(UNIT=iw,FMT="(/)")
    END DO

  END SUBROUTINE vib_out

  !!****f* vibrational_analysis/rot_ana [1.0] *
  !!
  !!   NAME
  !!     rot_ana
  !!
  !!   FUNCTION
  !!     Performs an analysis of the principal inertia axis
  !!     Getting back the generators of the translating and
  !!     rotating frame
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Teodoro Laino 08.2006
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE rot_ana(particles, mat, dof, iw, keep_rotations, natoms,&
       error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat
    INTEGER, INTENT(OUT)                     :: dof
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN)                      :: keep_rotations
    INTEGER, INTENT(IN)                      :: natoms
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rot_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iparticle, iseq, &
                                                iw2, j, k, lrot(3), stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cp(3), Ip(3,3), Ip_eigval(3), &
                                                mass, masst, norm, rcom(3), &
                                                rm(3)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Rot, Tr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I","",handle)
    logger => cp_error_get_logger(error)
    iw2 = cp_logger_get_default_unit_nr(logger)
    failure = .FALSE.
    CPPostcondition(ASSOCIATED(particles),cp_failure_level,routineP,error,failure)    
    IF (.NOT.keep_rotations) THEN
       rcom  = 0.0_dp
       masst = 0.0_dp
       ! Center of mass
       DO iparticle = 1, natoms
          mass = particles(iparticle)%atomic_kind%mass
          CPPostcondition(mass>0.0_dp,cp_failure_level,routineP,error,failure)
          masst = masst + mass
          rcom =  particles(iparticle)%r * mass + rcom
       END DO
       rcom = rcom / masst
       ! Intertia Tensor
       Ip   = 0.0_dp
       DO iparticle = 1, natoms
          mass = particles(iparticle)%atomic_kind%mass
          rm = particles(iparticle)%r - rcom
          Ip(1,1) = Ip(1,1) + mass * ( rm(2)**2 + rm(3)**2 )
          Ip(2,2) = Ip(2,2) + mass * ( rm(1)**2 + rm(3)**2 )
          Ip(3,3) = Ip(3,3) + mass * ( rm(1)**2 + rm(2)**2 )
          Ip(1,2) = Ip(1,2) - mass * ( rm(1)*rm(2) )
          Ip(1,3) = Ip(1,3) - mass * ( rm(1)*rm(3) )
          Ip(2,3) = Ip(2,3) - mass * ( rm(2)*rm(3) )
       END DO
       ! Diagonalize the Inertia Tensor
       CALL diamat_all(Ip,Ip_eigval)
       IF (iw>0) THEN
          WRITE(iw,'(T2,A,3F12.6)')"ROT| Rotational Analysis Info "
          WRITE(iw,'(T2,A)')"ROT| Principal axes and moments of inertia in atomic units:"
          WRITE(iw,'(T2,A,T32,3(3X,I6,6X))')"ROT|",1,2,3
          WRITE(iw,'(T2,A,T30,3(3X,F12.6))')"ROT| EIGENVALUES",Ip_eigval
          WRITE(iw,'(T2,A,T30,3(3X,F12.6))')"ROT|      X",Ip(1,1),Ip(1,2),Ip(1,3)
          WRITE(iw,'(T2,A,T30,3(3X,F12.6))')"ROT|      Y",Ip(2,1),Ip(2,2),Ip(2,3)
          WRITE(iw,'(T2,A,T30,3(3X,F12.6))')"ROT|      Z",Ip(3,1),Ip(3,2),Ip(3,3)
       END IF
       IF (iw>0) THEN
          WRITE(iw,'(/,T2,A)')"ROT| Standard Molecule Orientation - ANGSTROM "
          DO iparticle = 1, natoms
             WRITE(iw,'(T2,"ROT|",T20,A,T30,3(3X,F12.6))')&
                  TRIM(particles(iparticle)%atomic_kind%name),&
                  MATMUL(particles(iparticle)%r,Ip)*angstrom
          END DO
       END IF
    END IF
    ! Build up the Translational vectors
    ALLOCATE(Tr(natoms*3,3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    Tr   = 0.0_dp
    DO k = 1, 3
       iseq = 0
       DO iparticle = 1, natoms
          mass = SQRT(particles(iparticle)%atomic_kind%mass)
          DO j = 1, 3
             iseq = iseq + 1
             IF (j==k) Tr(iseq,k) = mass
          END DO
       END DO
    END DO
    ! Normalize Translations
    DO i = 1, 3
       norm = SQRT(DOT_PRODUCT(Tr(:,i),Tr(:,i)))
       Tr(:,i) = Tr(:,i)/norm
    END DO
    dof = 3
    ! Build up the Rotational vectors
    ALLOCATE(Rot(natoms*3,3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    lrot = 0
    IF (.NOT.keep_rotations) THEN
       DO iparticle = 1, natoms
          mass = SQRT(particles(iparticle)%atomic_kind%mass)
          rm = particles(iparticle)%r - rcom
          cp(1) = rm(1)*Ip(1,1) + rm(2) * Ip(2,1) + rm(3) * Ip(3,1)
          cp(2) = rm(1)*Ip(1,2) + rm(2) * Ip(2,2) + rm(3) * Ip(3,2)
          cp(3) = rm(1)*Ip(1,3) + rm(2) * Ip(2,3) + rm(3) * Ip(3,3)
          ! X Rot
          Rot((iparticle-1)*3+1,1) = (cp(2)*Ip(1,3)-Ip(1,2)*cp(3))*mass
          Rot((iparticle-1)*3+2,1) = (cp(2)*Ip(2,3)-Ip(2,2)*cp(3))*mass
          Rot((iparticle-1)*3+3,1) = (cp(2)*Ip(3,3)-Ip(3,2)*cp(3))*mass
          ! Y Rot
          Rot((iparticle-1)*3+1,2) = (cp(3)*Ip(1,1)-Ip(1,3)*cp(1))*mass
          Rot((iparticle-1)*3+2,2) = (cp(3)*Ip(2,1)-Ip(2,3)*cp(1))*mass
          Rot((iparticle-1)*3+3,2) = (cp(3)*Ip(3,1)-Ip(3,3)*cp(1))*mass
          ! Z Rot
          Rot((iparticle-1)*3+1,3) = (cp(1)*Ip(1,2)-Ip(1,1)*cp(2))*mass
          Rot((iparticle-1)*3+2,3) = (cp(1)*Ip(2,2)-Ip(2,1)*cp(2))*mass
          Rot((iparticle-1)*3+3,3) = (cp(1)*Ip(3,2)-Ip(3,1)*cp(2))*mass
       END DO

       ! Normalize Rotations and count the number of degree of freedom
       lrot = 1
       DO i = 1, 3
          norm = SQRT(DOT_PRODUCT(Rot(:,i),Rot(:,i)))
          IF (norm <= thrs) THEN
             lrot(i) = 0
             CYCLE
          END IF
          Rot(:,i) = Rot(:,i)/norm
          ! Clean Rotational modes for spurious/numerical contamination
          IF (i<3) THEN
             DO j = 1,i
                Rot(:,i+1) = Rot(:,i+1) - DOT_PRODUCT(Rot(:,i+1),Rot(:,j)) * Rot(:,j)
             END DO
          END IF
       END DO
    END IF
    dof = dof + COUNT(lrot==1)
    IF (iw>0) THEN
       WRITE(iw,'(T2,A,I6)')"ROT| Numer of Rotovibrational vectors:",dof
       IF (dof==5) WRITE(iw,'(T2,A)')"ROT| Linear Molecule detected.."
       IF (dof==3.AND.(.NOT.keep_rotations)) WRITE(iw,'(T2,A)')"ROT| Single Atom detected.." 
    END IF
    ! Give back the vectors generating the rototranslating Frame
    ALLOCATE(mat(natoms*3,dof),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    iseq = 0
    DO i = 1, 3
       mat(:,i)=Tr(:,i)
       IF (lrot(i)==1) THEN 
          iseq = iseq + 1
          mat(:,3+iseq) = Rot(:,i) 
       END IF
    END DO
    DEALLOCATE(Tr,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)   
    DEALLOCATE(Rot,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE rot_ana

  !!****f* vibrational_analysis/build_D_matrix [1.0] *
  !!
  !!   NAME
  !!     build_D_matrix
  !!
  !!   FUNCTION
  !!     Generates the transformation matrix from hessian in cartesian into
  !!     internal coordinates (based on Gram-Schmidt orthogonalization)
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Teodoro Laino 08.2006
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE build_D_matrix(particles,mat,dof,Dout,full,natoms,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat
    INTEGER, INTENT(IN)                      :: dof
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Dout
    LOGICAL, OPTIONAL                        :: full
    INTEGER, INTENT(IN)                      :: natoms
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_D_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ifound, iseq, j, &
                                                nvib, stat
    LOGICAL                                  :: failure, my_full
    REAL(KIND=dp)                            :: norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: D

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    my_full = .TRUE.
    IF (PRESENT(full)) my_full = full
    IF (.NOT.failure) THEN
       ! Generate the missing vectors of the orthogonal basis set
       nvib = 3*natoms-dof
       ALLOCATE(work(3*natoms),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(D(3*natoms,3*natoms),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Check First orthogonality in the first element of the basis set
       DO i = 1, dof
          D(:,i) = mat(:,i)
          DO j = i+1, dof
             norm = DOT_PRODUCT(mat(:,i),mat(:,j))
             CPPostcondition(ABS(norm)<thrs,cp_failure_level,routineP,error,failure)
          END DO
       END DO
       ! Generate the nvib orthogonal vectors
       iseq   = 0
       ifound = 0
       DO WHILE (ifound /= nvib)
          iseq = iseq + 1
          CPPostcondition(iseq<=3*natoms,cp_failure_level,routineP,error,failure)
          work       = 0.0_dp
          work(iseq) = 1.0_dp
          ! Gram Schmidt orthogonalization
          DO i = 1, dof+ifound
             norm = DOT_PRODUCT(work,D(:,i))
             work = work - norm * D(:,i)
          END DO
          ! Check norm of the new generated vector
          norm = SQRT(DOT_PRODUCT(work,work))
          IF (norm>=10E4_dp*thrs) THEN
             ! Accept new vector
             ifound = ifound + 1
             D(:,dof+ifound) = work / norm
          END IF
       END DO
       CPPostcondition(dof+ifound==3*natoms,cp_failure_level,routineP,error,failure)
       IF (my_full) THEN
          ALLOCATE(Dout(3*natoms,3*natoms),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
          Dout = D
       ELSE
          ALLOCATE(Dout(3*natoms,nvib),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
          Dout = D(:,dof+1:)
       END IF
       DEALLOCATE(work,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(D,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE build_D_matrix

  !!****f* vibrational_analysis/MS_initial_moves [1.0] *
  !!
  !!   NAME
  !!     MS_initial_moves
  !!
  !!   FUNCTION
  !!     Generates the first displacement vector for a mode selctive vibrational
  !!     analysis. At the moment this is a random number for selected atoms 
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Florian Schiffmann 11.2007
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE MS_initial_moves(para_env, nrep , input, globenv ,particles ,&
                              delta_vec , b_vec, step_r,step_b, sel_freq,&
                              mass,hes_bfgs,eig_bfgs,precondit, dx,output_unit_ms,&
                              b_mat,s_mat,mat_size, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: nrep
    TYPE(section_vals_type), POINTER         :: input
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:, :)           :: delta_vec, b_vec
    REAL(Kind=dp), DIMENSION(:)              :: step_r, step_b
    REAL(KIND=dp)                            :: sel_freq
    REAL(Kind=dp), DIMENSION(:)              :: mass
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hes_bfgs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eig_bfgs
    INTEGER                                  :: precondit
    REAL(KIND=dp)                            :: dx
    INTEGER                                  :: output_unit_ms
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: b_mat, s_mat
    INTEGER                                  :: mat_size
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'MS_initial_moves', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: guess, handle, hesunit, i, j, &
                                                jj, k, m, n_rep_val, natoms, &
                                                ncoord, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map_atoms
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure, ionode
    REAL(KIND=dp)                            :: my_val, norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: hes, tmp
    TYPE(section_vals_type), POINTER         :: ms_vib_section

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    IF (.NOT. failure) THEN
       ms_vib_section => section_vals_get_subs_vals(input,"VIBRATIONAL_ANALYSIS%MODE_SELECTIVE", error=error)
       CALL section_vals_val_get(ms_vib_section,"INITIAL_GUESS",i_val=guess,error=error)
       CALL section_vals_val_get(ms_vib_section,"PRECONDITIONER",i_val=precondit,error=error)      
       CALL section_vals_val_get(ms_vib_section,"FREQUENCY",r_val=sel_freq, error=error)
       
       ionode=para_env%mepos==para_env%source
       SELECT CASE(guess)
       CASE(ms_guess_atomic)
          CALL section_vals_val_get(ms_vib_section,"ATOMS",n_rep_val=n_rep_val,error=error)
          jj = 0
          DO k = 1,n_rep_val
             CALL section_vals_val_get(ms_vib_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
             DO  j = 1,SIZE(tmplist)
                jj = jj+1
             END DO
          END DO
          IF (jj<1)THEN
             natoms=SIZE(particles)
             ALLOCATE(map_atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO j=1,natoms
                map_atoms(j)=j
             END DO
          ELSE
             natoms=jj
             ALLOCATE(map_atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             jj = 0
             DO m=1,n_rep_val
                CALL section_vals_val_get(ms_vib_section,"ATOMS",i_rep_val=m,i_vals=tmplist, error=error)
                DO  j = 1,SIZE(tmplist)
                   map_atoms (j) = tmplist(j)
                END DO
             END DO
          END IF
          
          ! apply random displacement along the mass weighted nuclear cartesian coordinates
          b_vec=0._dp
          delta_vec=0._dp
          jj=0
          
          DO i=1,nrep
             DO j=1,natoms
                DO k=1,3
                   jj=(map_atoms(j)-1)*3+k
                   b_vec(jj,i)=ABS(next_random_number(globenv%gaussian_rng_stream,error=error))
                END DO
             END DO
             norm=SQRT(DOT_PRODUCT(b_vec(:,i),b_vec(:,i)))
             b_vec(:,i)=b_vec(:,i)/norm
          END DO
          
          IF(nrep.GT.1)THEN
             DO k=1,10
                DO j=1,nrep
                   DO i=1,nrep
                      IF(i.NE.j)THEN
                         b_vec(:,j)=b_vec(:,j)-DOT_PRODUCT(b_vec(:,j),b_vec(:,i))* b_vec(:,i)
                         b_vec(:,j)=b_vec(:,j)/SQRT(DOT_PRODUCT(b_vec(:,j),b_vec(:,j)))
                      END IF
                   END DO
                END DO
             END DO
          END IF
          DO i=1,SIZE(b_vec,1)
             delta_vec(i,:)=b_vec(i,:)/mass(i)
          END DO
          mat_size=0
       CASE(ms_guess_bfgs)
          natoms=SIZE(particles)
          ncoord=3*natoms

          ALLOCATE(hes_bfgs(ncoord,ncoord),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(eig_bfgs(ncoord),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          IF(ionode)THEN
             CALL open_file(file_name="HESSIAN",&
                  file_status="UNKNOWN",&
                  file_form="UNFORMATTED",&
                  file_action="READ",&
                  unit_number=hesunit)
             ALLOCATE (hes(ncoord*(ncoord+1)/2),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmp(ncoord),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmplist(ncoord),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

             READ(UNIT=hesunit,IOSTAT=stat) hes
             CALL close_file(hesunit)
             IF (stat /= 0) THEN
                WRITE(6,FMT="(/,T2,A)")"**  Error while reading HESSIAN **"
             ELSE
                WRITE(6,FMT="(/,T2,A)")&
                     "*** Initial Hessian has been read successfully ***"
             END IF
             DO i = 1, ncoord
                DO j = 1, i
                   jj = i*(i-1)/2 + j
                   hes_bfgs(i,j) = hes(jj)
                   hes_bfgs(j,i) = hes(jj)
                END DO
             END DO
             DO i = 1, ncoord
                DO j = 1, ncoord
                   hes_bfgs(i,j) = hes_bfgs(i,j)/(mass(i)*mass(j))
                END DO
             END DO
             CALL diamat_all(hes_bfgs,eig_bfgs)
             my_val=(sel_freq/(vibfac))**2/massunit
             DO i=1,ncoord
                tmp(i)=ABS(my_val-eig_bfgs(i))
             END DO
             CALL sort(tmp,ncoord,tmplist)
             DO i=1,nrep
                b_vec(:,i)=hes_bfgs(:,tmplist(i))
                norm=SQRT(DOT_PRODUCT(b_vec(:,i),b_vec(:,i)))
                b_vec(:,i)=b_vec(:,i)/norm
             END DO
             DO i=1,SIZE(b_vec,1)
                delta_vec(i,:)=b_vec(i,:)/mass(i)
             END DO
             DEALLOCATE(hes,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(tmp,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(tmplist,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(precondit.NE.bfgs_ms_precond)THEN
             DEALLOCATE(hes_bfgs,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(eig_bfgs,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          mat_size=0
       CASE(ms_guess_restart)
          natoms=SIZE(particles)
          ncoord=3*natoms
          IF(ionode)THEN
             CALL open_file(file_name="MS_RESTART",&
                  file_status="UNKNOWN",&
                  file_form="UNFORMATTED",&
                  file_action="READ",&
                  unit_number=hesunit) 
             READ(UNIT=hesunit,IOSTAT=stat)mat_size
          END IF
          CALL mp_bcast(mat_size,para_env%source,para_env%group)
          ALLOCATE(b_mat(ncoord,mat_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(s_mat(ncoord,mat_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF(ionode)THEN
             READ(UNIT=hesunit,IOSTAT=stat)b_mat
             READ(UNIT=hesunit,IOSTAT=stat)s_mat
             READ(UNIT=hesunit,IOSTAT=stat)b_vec
             READ(UNIT=hesunit,IOSTAT=stat)delta_vec 
             CALL close_file(hesunit)
             IF (stat /= 0) THEN
                WRITE(6,FMT="(/,T2,A)")"**  Error while reading MS_RESTART **"
             ELSE
                WRITE(6,FMT="(/,T2,A)")"*** RESTART has been read successfully ***"
             END IF
          END IF
          CALL mp_bcast(b_mat,para_env%source,para_env%group)
          CALL mp_bcast(s_mat,para_env%source,para_env%group)
       END SELECT
       CALL mp_bcast(b_vec,para_env%source,para_env%group)
       CALL mp_bcast(delta_vec,para_env%source,para_env%group)
    END IF
    DO i=1,nrep
       step_r(i)=dx/SQRT(DOT_PRODUCT( delta_vec(:,i), delta_vec(:,i)))
       step_b(i)=SQRT(DOT_PRODUCT(step_r(i)*b_vec(:,i),step_r(i)*b_vec(:,i)))
    END DO
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE MS_initial_moves


  !!****f* vibrational_analysis/evaluate_H_update_b [1.0] *
  !!
  !!   NAME
  !!     evaluate_H_update_b
  !!
  !!   FUNCTION
  !!     Davidson algorithm for to generate a approximate Hessian for mode
  !!     selective vibrational analysis
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Florian Schiffmann 11.2007
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************


  SUBROUTINE  evaluate_H_update_b(mat_size,ms_force,rep_env,b_vec,nrep,&
                                  particles,b_mat,s_mat,step_b,&
                                  step_r,delta_vec, mass,sel_freq, &
                                  converged,hes_bfgs,eig_bfgs,precondit,dx,output_unit_ms,&
                                  output_molden,error)     
    INTEGER                                  :: mat_size
    REAL(KIND=dp), DIMENSION(:, :)           :: ms_force
    TYPE(replica_env_type), POINTER          :: rep_env
    REAL(KIND=dp), DIMENSION(:, :)           :: b_vec
    INTEGER                                  :: nrep
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: b_mat, s_mat
    REAL(Kind=dp), DIMENSION(:)              :: step_b, step_r
    REAL(KIND=dp), DIMENSION(:, :)           :: delta_vec
    REAL(Kind=dp), DIMENSION(:)              :: mass
    REAL(KIND=dp)                            :: sel_freq
    LOGICAL                                  :: converged
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hes_bfgs
    REAL(Kind=dp), DIMENSION(:), POINTER     :: eig_bfgs
    INTEGER                                  :: precondit
    REAL(KIND=dp)                            :: dx
    INTEGER                                  :: output_unit_ms, output_molden
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_H_update_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, jj, k, natoms, ncoord, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_val, norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenval, freq, tmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: approx_H, H_save, residuum, &
                                                selected, tmp_b, tmp_s, Xpre
    REAL(KIND=dp), DIMENSION(2, nrep)        :: criteria

    natoms=SIZE(particles)
    ncoord=3*natoms
    nrep=SIZE(rep_env%f,2)
    failure=.FALSE.

    !!!!!!!!   reallocate and update the davidson matrices   !!!!!!!!!!
     IF(mat_size.NE.0)THEN

       ALLOCATE(tmp_b(3*natoms,mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_s(3*natoms,mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       tmp_b=b_mat
       tmp_s=s_mat

       DEALLOCATE(b_mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(s_mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    
    ALLOCATE(b_mat(3*natoms,mat_size+nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(s_mat(3*natoms,mat_size+nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    s_mat=0.0_dp

    DO i=1,3*natoms
       IF(mat_size.NE.0)THEN
          DO j=1,mat_size
             b_mat(i,j)=tmp_b(i,j)
             s_mat(i,j)=tmp_s(i,j) 
          END DO
       END IF 
       DO j=1,nrep
          b_mat(i,mat_size+j)=b_vec(i,j) 
       END DO
    END DO

    IF(mat_size.NE.0)THEN
       DEALLOCATE(tmp_s,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_b,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    mat_size=mat_size+nrep

    ALLOCATE(approx_H(mat_size,mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(H_save(mat_size,mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eigenval(mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    !!!!!!!!!!!!  calculate the new derivativand the approximate hessian

    DO i=1,nrep
       DO j=1,3*natoms
          s_mat(j,mat_size-nrep+i)= -(ms_force(j,i)-rep_env%f(j,i))/(2*step_b(i)*mass(j))
       END DO
    END DO

    CALL dgemm('T','N',mat_size,mat_size,SIZE(s_mat,1),1._dp,b_mat,SIZE(b_mat,1),&
         s_mat,SIZE(s_mat,1),0._dp,approx_H,mat_size)
    H_save=approx_H
    CALL diamat_all(approx_H,eigenval)

    !!!!!!!!!!!! select eigenvalue(s) and vector(s) and calculate the new displacement vector

    ALLOCATE(selected(mat_size,nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp(mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ind(mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(residuum(SIZE(s_mat,1),nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    my_val=(sel_freq/(vibfac))**2/massunit
    DO i=1,mat_size
       tmp(i)=ABS(my_val-eigenval(i))
    END DO
    CALL sort(tmp,(mat_size),ind)
    residuum=0._dp

    DO j=1,nrep
       DO i=1,mat_size
          residuum(:,j)=residuum(:,j)+approx_H(i,ind(j))*(s_mat(:,i)-eigenval(ind(j))*b_mat(:,i))
       END DO
    END DO   
    


    DO j=1,nrep
       DO i=1,mat_size
          residuum(:,j)=residuum(:,j)-DOT_PRODUCT(residuum(:,j),b_mat(:,i))* b_mat(:,i)
       END DO
    END DO
    
    DO i=1,nrep
       criteria(1,i)=MAXVAL((residuum(:,i)))
       criteria(2,i)=SQRT(DOT_PRODUCT(residuum(:,i),residuum(:,i)))
    END DO
    DO k=1,10
       DO j=1,nrep
             DO i=1,mat_size
                residuum(:,j)=residuum(:,j)-DOT_PRODUCT(residuum(:,j),b_mat(:,i))* b_mat(:,i)
                residuum(:,j)=residuum(:,j)/SQRT(DOT_PRODUCT(residuum(:,j),residuum(:,j)))
             END DO
             IF(nrep.GT.1)THEN
                DO i=1,nrep
                   IF(i.NE.j)THEN
                      residuum(:,j)=residuum(:,j)-DOT_PRODUCT(residuum(:,j),residuum(:,i))* residuum(:,i)
                      residuum(:,j)=residuum(:,j)/SQRT(DOT_PRODUCT(residuum(:,j),residuum(:,j)))
                   END IF
                END DO
             END IF
       END DO
    END DO

    SELECT CASE(precondit)
    CASE(no_ms_precond)
       DO i=1,nrep
          norm=SQRT(DOT_PRODUCT(residuum(:,i),residuum(:,i)))     
          b_vec(:,i)=residuum(:,i)/norm
       END DO
    CASE(bfgs_ms_precond)
       CPPostcondition(SIZE(hes_bfgs,1)==ncoord,cp_failure_level,routineP,error,failure)
       CPPostcondition(SIZE(eig_bfgs,1)==ncoord,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_s(ncoord,ncoord),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_b(ncoord,ncoord),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Xpre(ncoord,ncoord),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO i=1,ncoord
          tmp_s(i,i)=eig_bfgs(i)-eigenval(ind(1))
          tmp_s(i,i)=1/(MAX(0.0001_dp,ABS(tmp_s(i,i))))
       END DO
       
       CALL dgemm('N','T',ncoord,ncoord,ncoord,1._dp,tmp_s,ncoord,&
            hes_bfgs,ncoord,0._dp,tmp_b,ncoord)
       CALL dgemm('N','N',ncoord,ncoord,ncoord,1._dp,hes_bfgs,ncoord,&
            tmp_s,ncoord,0._dp,Xpre,ncoord)
              
       DO i=1,nrep
          b_vec(:,i)=MATMUL(Xpre,residuum(:,i))
          norm=SQRT(DOT_PRODUCT(b_vec(:,i),b_vec(:,i)))
          b_vec(:,i)=b_vec(:,i)/norm
       END DO
       DEALLOCATE(tmp_s,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_b,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(Xpre,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       

    END SELECT
    
    DO i=1,nrep
       DO j=1,natoms
          DO k=1,3
             jj=(j-1)*3+k
             delta_vec(jj,i)= b_vec(jj,i)/mass(jj)
          END DO
       END DO
    END DO

    DO i=1,nrep
       step_r(i)=dx/SQRT(DOT_PRODUCT( delta_vec(:,i), delta_vec(:,i)))
       step_b(i)=SQRT(DOT_PRODUCT(step_r(i)*b_vec(:,i),step_r(i)*b_vec(:,i)))
    END DO
    converged=.FALSE.
    IF(MAXVAL(criteria(1,:)).LE.1E-7.AND.MAXVAL(criteria(2,:)).LE.1E-6.OR.mat_size==ncoord)converged=.TRUE.
    ALLOCATE(freq(nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
    DO i=1,nrep
       freq(i)=SQRT(ABS(eigenval(ind(i)))*massunit)*vibfac
    END DO

    !!!   write information and output   !!!
    CALL ms_out(output_unit_ms,converged,freq,criteria ,b_mat,s_mat,b_vec,delta_vec,mat_size,nrep,error)
    IF(converged)THEN
       eigenval(:)=SQRT(ABS(eigenval(:))*massunit)*vibfac
       ALLOCATE(tmp_b(ncoord,mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_b=0._dp
       DO i=1,mat_size
          DO j=1,mat_size
             tmp_b(:,i)=tmp_b(:,i)+approx_H(j,i)*b_mat(:,j)
          END DO
       END DO
       CALL molden_out(output_molden,particles,eigenval,tmp_b,error)
       DEALLOCATE(tmp_b,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    DEALLOCATE(freq,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(approx_H,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eigenval,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(selected,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ind,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(residuum,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE evaluate_H_update_b

  !!****f* vibrational_analysis/ms_out [1.0] *
  !!
  !!   NAME
  !!     ms_out
  !!
  !!   FUNCTION
  !!     writes the output for a mode tracking calculation
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Florian Schiffmann 11.2007
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE ms_out(iw,converged,freq,criter ,b_mat,s_mat,b_vec,delta_vec,mat_size,nrep,error) 
    INTEGER                                  :: iw
    LOGICAL                                  :: converged
    REAL(KIND=dp), DIMENSION(:)              :: freq
    REAL(KIND=dp), DIMENSION(:, :)           :: criter
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: b_mat, s_mat
    REAL(KIND=dp), DIMENSION(:, :)           :: b_vec, delta_vec
    INTEGER                                  :: mat_size, nrep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ms_out', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, msunit, stat

    IF(converged)THEN
       IF(iw.GT.0)THEN
          WRITE(iw,'(T2,A)')"MS| DAVIDSON ALGORITHM CONVERGED"   
          DO i=1,nrep
             WRITE(iw,'(T2,"MS| TRACKED FREQUENCY (",I0,") IS:",F12.6,3X,A)') i,freq(i),'cm-1'
          END DO
       END IF
    ELSE
       IF(iw.GT.0)THEN
          CALL open_file(file_name="MS_RESTART",&
               file_status="REPLACE",&
               file_form="UNFORMATTED",&
               file_action="WRITE",&
               unit_number=msunit)
          
          WRITE(UNIT=msunit,IOSTAT=stat)mat_size
          WRITE(UNIT=msunit,IOSTAT=stat)b_mat
          WRITE(UNIT=msunit,IOSTAT=stat)s_mat
          WRITE(UNIT=msunit,IOSTAT=stat)b_vec
          WRITE(UNIT=msunit,IOSTAT=stat)delta_vec 
          
          CALL close_file(unit_number=msunit)
       
          WRITE(iw,'(T2,A,3X,I6)')"MS| ITERATION STEP", mat_size/nrep
          DO i=1,nrep
             IF(criter(1,i).LE.1E-7.AND.(criter(2,i)).LE.1E-6)THEN
                WRITE(iw,'(T2,A,3X,F12.6,A)')"MS| TRACKED MODE ", freq(i) , "cm-1  IS  CONVERGED"
             ELSE
                WRITE(iw,'(T2,A,3X,F12.6,A)')"MS| TRACKED MODE ", freq(i) , "cm-1  NOT  CONVERGED"
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE ms_out

  !!****f* vibrational_analysis/ms_out [1.0] *
  !!
  !!   NAME
  !!     ms_out
  !!
  !!   FUNCTION
  !!     writes the output for a mode tracking calculation
  !!
  !!   NOTES
  !!     
  !!
  !!   AUTHOR
  !!     Florian Schiffmann 11.2007
  !!
  !!   MODIFICATION HISTORY
  !!             
  !!
  !!   SOURCE
  !****************************************************************************
  SUBROUTINE molden_out(iw,particles,freq,eigen_vec,error)
    INTEGER                                  :: iw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:)              :: freq
    REAL(KIND=dp), DIMENSION(:, :)           :: eigen_vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'molden_out', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: i, j, k
    LOGICAL                                  :: failure

    IF(iw.GT.0)THEN
       CPPostcondition(MOD(SIZE(eigen_vec,1),3)==0,cp_failure_level,routineP,error,failure)    
       CPPostcondition(SIZE(particles)==SIZE(eigen_vec,1)/3,cp_failure_level,routineP,error,failure)
       CPPostcondition(SIZE(freq,1)==SIZE(eigen_vec,2),cp_failure_level,routineP,error,failure)       
       WRITE(iw,'(T2,A)')"[Molden Format]"   
       WRITE(iw,'(T2,A)')"[FREQ]"
       DO i=1,SIZE(freq,1)
          WRITE(iw,'(T5,F12.6)') freq(i)
       END DO
       WRITE(iw,'(T2,A)')"[FR-COORD]"       
       DO i=1,SIZE(particles)
          CALL get_atomic_kind(atomic_kind=particles(i)%atomic_kind,&
               element_symbol=element_symbol) 
               WRITE(iw,'(T2,A2,3X,3(F12.6,3X))')&
                    element_symbol, particles((i))%r(:)
       END DO
       WRITE(iw,'(T2,A)')"[FR-NORM-COORD]"
       DO i=1,SIZE(eigen_vec,2)
          WRITE(iw,'(T2,A,1X,I6)')"vibration",i
          DO j=1,SIZE(eigen_vec,1)/3
             k=(j-1)*3
             WRITE(iw,'(T2,3(F12.6,3X))')eigen_vec(k+1,i),eigen_vec(k+2,i),eigen_vec(k+3,i)
          END DO
       END DO
    END IF
  END SUBROUTINE molden_out

END MODULE vibrational_analysis

