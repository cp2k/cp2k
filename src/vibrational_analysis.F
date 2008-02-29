!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module performing a vibrational analysis
!> \note
!>      Numerical accuracy for parallel runs:
!>       Each replica starts the SCF run from the one optimized
!>       in a previous run. It may happen then energies and derivatives
!>       of a serial run and a parallel run could be slightly different
!>       'cause of a different starting density matrix.
!>       Exact results are obtained using:
!>          EXTRAPOLATION USE_GUESS in QS section (Teo 08.2006)
!> \author Teodoro Laino 08.2006
! *****************************************************************************
MODULE vibrational_analysis
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: description_length,&
                                             get_results
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
  USE message_passing,                 ONLY: mp_sum
  USE mode_selective,                  ONLY: molden_out,&
                                             ms_vb_anal
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
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
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'vibrational_analysis'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: thrs = 100.0_dp*EPSILON(0.0_dp)
  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.

  PUBLIC :: vb_anal

CONTAINS

! *****************************************************************************
!> \brief Module performing a vibrational analysis
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE vb_anal(input, para_env, globenv, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vb_anal', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1), DIMENSION(3)           :: lab = (/"X","Y","Z"/)
    CHARACTER(LEN=description_length)        :: description
    INTEGER :: handle, i, icoord, icoordm, icoordp, ierr, imap, ip1, ip2, &
      iparticle1, iparticle2, iseq, iw, j, k, n_rep, natoms, ncoord, nrep, &
      nres, nRotTrM, nvib, output_molden, output_unit, prep, stat
    INTEGER, DIMENSION(:), POINTER           :: Clist, Mlist
    LOGICAL                                  :: calc_intens, &
                                                do_mode_tracking, failure, &
                                                keep_rotations, &
                                                something_frozen
    REAL(KIND=dp)                            :: dx, norm, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: H_eigval1, H_eigval2, konst, &
                                                mass, pos0, rmass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Hessian, Hint1, Hint2
    REAL(KIND=dp), DIMENSION(3)              :: D_deriv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: intensities
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: D, dip_deriv, RotTrM
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: tmp_dip
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
       NULLIFY(D, RotTrM, logger, subsys, f_env, particles, rep_env, &
            intensities)
       logger => cp_error_get_logger(error)
       output_unit=cp_print_key_unit_nr(logger,input,"VIBRATIONAL_ANALYSIS%PROGRAM_RUN_INFO",&
            extension=".vibLog",error=error)

       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%DX",&
            r_val=dx, error=error)
       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%NPROC_REP",&
            i_val=prep, error=error)
       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%FULLY_PERIODIC",&
            l_val=keep_rotations, error=error)
       CALL section_vals_val_get(input,"VIBRATIONAL_ANALYSIS%INTENSITIES",&
            l_val=calc_intens, error=error)
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
             CALL ms_vb_anal(input,rep_env, para_env, globenv, particles,&
                  nrep,calc_intens,dx,output_unit,output_molden, error)
             CALL f_env_rm_defaults(f_env,new_error,ierr)
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
             IF(calc_intens)THEN
                description='[DIPOLE]'
                ALLOCATE(tmp_dip(ncoord,3,2),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
                tmp_dip=0._dp
             END IF
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
                   IF(calc_intens)THEN
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           n_rep=nres,error=error)
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           values=tmp_dip(icoord+j,:,1),&
                           nval=nres,error=error)
                   END IF
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
                   IF(calc_intens)THEN
                      k=(icoord+j+2)/3
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           n_rep=nres,error=error)
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           values=tmp_dip(icoord+j,:,2),&
                           nval=nres,error=error)
                      tmp_dip(icoord+j,:,1)=(tmp_dip(icoord+j,:,1)-tmp_dip(icoord+j,:,2))/(2.0_dp*Dx*mass(k))
                   END IF
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
             IF(calc_intens)THEN
                ALLOCATE(dip_deriv(3,SIZE(D,2)),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                dip_deriv=0.0_dp
             END IF
             ALLOCATE(intensities(SIZE(D,2)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             intensities=  0._dp
             Hint1      = Hessian
             CALL diamat_all(Hint1,H_eigval1)
             IF (output_unit>0) THEN
                WRITE(output_unit,'(T2,"VIB| Cartesian Low frequencies ---",4G12.5)')&
                     (H_eigval1(i),i=1,MIN(9,ncoord))
             END IF
             IF (nvib/=0) THEN
                Hint2=MATMUL(TRANSPOSE(D),MATMUL(Hessian,D))
                IF(calc_intens)THEN
                   DO i=1,3
                      dip_deriv(i,:)=MATMUL(tmp_dip(:,i,1),D)
                   END DO
                END IF
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

                   IF(calc_intens)THEN
                      D_deriv=0._dp
                      DO j=1,nvib
                         D_deriv (:)=D_deriv(:)+dip_deriv(:,j)*Hint2(j,i)
                      END DO
                      intensities(i)=SQRT(DOT_PRODUCT(D_deriv,D_deriv))
                   END IF
                   ! Convert frequencies to cm^-1
                   H_eigval2(i) =  SIGN(1.0_dp,H_eigval2(i))*SQRT(ABS(H_eigval2(i))*massunit)*vibfac/1000.0_dp
                END DO
                ! Dump Info
                IF (para_env%mepos==para_env%source) THEN             
                   iw = cp_logger_get_default_unit_nr(logger)
                   CALL vib_out(iw,nvib,D,konst,rmass,H_eigval2,particles,Mlist, intensities)
                END IF
                IF(.NOT.something_frozen) CALL molden_out(input,particles,H_eigval2,D,intensities,calc_intens,error)
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
             IF(calc_intens)THEN
                DEALLOCATE(dip_deriv,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(tmp_dip,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             DEALLOCATE(intensities,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL f_env_rm_defaults(f_env,new_error,ierr)
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,input,"VIBRATIONAL_ANALYSIS%PROGRAM_RUN_INFO",error=error)
       CALL rep_env_release(rep_env,error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE vb_anal

! *****************************************************************************
!> \brief give back a list of moving atoms
!> \author Teodoro Laino 08.2006
! *****************************************************************************
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

! *****************************************************************************
!> \brief Dumps results of the vibrational analysis
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE vib_out(iw,nvib,D,k,m,freq,particles,Mlist,intensities)
    INTEGER, INTENT(IN)                      :: iw, nvib
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: D
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: k, m, freq
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER, DIMENSION(:), POINTER           :: Mlist
    REAL(KIND=dp), DIMENSION(:), POINTER     :: intensities

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
       IF (ASSOCIATED(intensities)) THEN
          WRITE(UNIT=iw,FMT="(T2,'VIB|Intensities      ',3(1X,F12.6,8X))")&
               (intensities(icol),icol=from,to)
       END IF
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

! *****************************************************************************
!> \brief Performs an analysis of the principal inertia axis
!>      Getting back the generators of the translating and
!>      rotating frame
!> \author Teodoro Laino 08.2006
! *****************************************************************************
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

! *****************************************************************************
!> \brief Generates the transformation matrix from hessian in cartesian into
!>      internal coordinates (based on Gram-Schmidt orthogonalization)
!> \author Teodoro Laino 08.2006
! *****************************************************************************
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


END MODULE vibrational_analysis

