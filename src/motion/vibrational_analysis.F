!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE cp_result_methods,               ONLY: get_results
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: vib_header
  USE input_constants,                 ONLY: do_rep_blocked
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: diamat_all
  USE mode_selective,                  ONLY: ms_vb_anal
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molden_utils,                    ONLY: molden_out
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE motion_utils,                    ONLY: rot_ana,&
                                             thrs_motion
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_methods,                ONLY: write_particle_matrix
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: massunit,&
                                             vibfac
  USE replica_methods,                 ONLY: rep_env_calc_e_f,&
                                             rep_env_create
  USE replica_types,                   ONLY: rep_env_release,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'vibrational_analysis'
  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.

  PUBLIC :: vb_anal

CONTAINS

! *****************************************************************************
!> \brief Module performing a vibrational analysis
!> \param input ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE vb_anal(input, input_declaration, para_env, globenv)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'vb_anal', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), DIMENSION(3), &
      PARAMETER                              :: lab = (/"X","Y","Z"/)

    CHARACTER(LEN=default_string_length)     :: description
    INTEGER :: handle, i, icoord, icoordm, icoordp, ierr, imap, ip1, ip2, &
      iparticle1, iparticle2, iseq, iw, j, k, natoms, ncoord, nrep, nres, &
      nRotTrM, nvib, output_unit, prep, proc_dist_type
    INTEGER, DIMENSION(:), POINTER           :: Clist, Mlist
    LOGICAL                                  :: calc_intens, &
                                                do_mode_tracking, failure, &
                                                keep_rotations, row_force, &
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
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                mode_tracking_section, &
                                                print_section, vib_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(D, RotTrM, logger, subsys, f_env, particles, rep_env, intensities,&
            vib_section, print_section)
    logger => cp_get_default_logger()
    vib_section   => section_vals_get_subs_vals(input,"VIBRATIONAL_ANALYSIS")
    print_section => section_vals_get_subs_vals(vib_section,"PRINT")
    output_unit=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",extension=".vibLog")

    CALL section_vals_val_get(vib_section,"DX",r_val=dx)
    CALL section_vals_val_get(vib_section,"NPROC_REP",i_val=prep)
    CALL section_vals_val_get(vib_section,"PROC_DIST_TYPE",i_val=proc_dist_type)
    row_force = (proc_dist_type==do_rep_blocked)
    CALL section_vals_val_get(vib_section,"FULLY_PERIODIC",l_val=keep_rotations)
    CALL section_vals_val_get(vib_section,"INTENSITIES",l_val=calc_intens)
    mode_tracking_section=>section_vals_get_subs_vals(vib_section,"MODE_SELECTIVE")
    CALL section_vals_get(mode_tracking_section,explicit=do_mode_tracking)
    nrep = MAX(1,para_env%num_pe/prep)
    prep = para_env%num_pe/nrep
    iw=cp_print_key_unit_nr(logger,print_section,"BANNER",extension=".vibLog")
    CALL vib_header(iw, nrep, prep)
    CALL cp_print_key_finished_output(iw,logger,print_section,"BANNER")
    ! Just one force_env allowed
    force_env_section => section_vals_get_subs_vals(input,"FORCE_EVAL")
    ! Create Replica Environments
    CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          input_declaration=input_declaration,nrep=nrep,prep=prep,row_force=row_force)
    IF (ASSOCIATED(rep_env)) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            failure=failure)
       CALL force_env_get(f_env%force_env,subsys=subsys)
       particles => subsys%particles%els
       ! Decide which kind of Vibrational Analysis to perform
       IF (do_mode_tracking)THEN
          CALL ms_vb_anal(input,rep_env, para_env, globenv, particles,&
                          nrep,calc_intens,dx,output_unit, logger)
          CALL f_env_rm_defaults(f_env,ierr)
       ELSE
          CALL get_moving_atoms(force_env=f_env%force_env,Ilist=Mlist)
          something_frozen = SIZE(particles).NE.SIZE(Mlist)
          natoms=SIZE(Mlist)
          ncoord=natoms*3
          ALLOCATE(Clist(ncoord))
          ALLOCATE(mass(natoms))
          ALLOCATE(pos0(ncoord))
          ALLOCATE(Hessian(ncoord,ncoord))
          IF(calc_intens)THEN
             description='[DIPOLE]'
             ALLOCATE(tmp_dip(ncoord,3,2))
             tmp_dip=0._dp
          END IF
          Clist = 0
          DO i = 1, natoms
             imap             = Mlist(i)
             Clist((i-1)*3+1) = (imap-1)*3+1
             Clist((i-1)*3+2) = (imap-1)*3+2
             Clist((i-1)*3+3) = (imap-1)*3+3
             mass(i) = particles(imap)%atomic_kind%mass
             CPASSERT(mass(i)>0.0_dp)
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
             ALLOCATE(RotTrM(natoms*3,nRotTrM))
          ELSE
             CALL rot_ana(particles,RotTrM,nRotTrM,print_section,&
                  keep_rotations,mass_weighted=.TRUE.,natoms=natoms)
          END IF
          ! Generate the suitable rototranslating basis set
          CALL build_D_matrix(RotTrM,nRotTrM,D,full=.FALSE.,&
               natoms=natoms)
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
             CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.)

             DO j = 1, nrep
                IF(calc_intens)THEN
                   IF (icoord+j <= ncoord) THEN
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           n_rep=nres)
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           values=tmp_dip(icoord+j,:,1),&
                           nval=nres)
                   END IF
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
             CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.)

             DO j = 1, nrep
                IF (calc_intens) THEN
                   IF (icoord+j <= ncoord) THEN
                      k=(icoord+j+2)/3
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           n_rep=nres)
                      CALL get_results(results=rep_env%results(j)%results,&
                           description=description,&
                           values=tmp_dip(icoord+j,:,2),&
                           nval=nres)
                      tmp_dip(icoord+j,:,1)=(tmp_dip(icoord+j,:,1)-tmp_dip(icoord+j,:,2))/(2.0_dp*Dx*mass(k))
                   END IF
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

          ! restore original particle positions for output
          DO i = 1, natoms
             imap             = Mlist(i)
             particles(imap)%r(1:3)=pos0((i-1)*3+1:(i-1)*3+3)
          ENDDO
          ! Dump Info
          IF (output_unit>0) THEN
             WRITE(output_unit,'(T2,A)')"VIB| Hessian in cartesian coordinates"
             CALL write_particle_matrix(Hessian,particles,output_unit,el_per_part=3,&
                  Ilist=Mlist)
          END IF
          ! Enforce symmetry in the Hessian
          DO i =1, ncoord
             DO j = i, ncoord
                ! Take the upper diagonal part
                Hessian(j,i) = Hessian(i,j)
             END DO
          END DO
          nvib = ncoord-nRotTrM
          ALLOCATE(H_eigval1(ncoord))
          ALLOCATE(H_eigval2(SIZE(D,2)))
          ALLOCATE(Hint1(ncoord,ncoord))
          ALLOCATE(Hint2(SIZE(D,2),SIZE(D,2)))
          ALLOCATE(rmass(SIZE(D,2)))
          ALLOCATE(konst(SIZE(D,2)))
          IF(calc_intens)THEN
             ALLOCATE(dip_deriv(3,SIZE(D,2)))
             dip_deriv=0.0_dp
          END IF
          ALLOCATE(intensities(SIZE(D,2)))
          intensities=  0._dp
          Hint1(:,:) = Hessian
          CALL diamat_all(Hint1,H_eigval1)
          IF (output_unit>0) THEN
             WRITE(output_unit,'(T2,"VIB| Cartesian Low frequencies ---",4G12.5)')&
                  (H_eigval1(i),i=1,MIN(9,ncoord))
          END IF
          IF (nvib/=0) THEN
             Hint2(:,:)=MATMUL(TRANSPOSE(D),MATMUL(Hessian,D))
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
             iw = cp_logger_get_default_io_unit(logger)
             IF (iw>0) THEN
                CALL vib_out(iw,nvib,D,konst,rmass,H_eigval2,particles,Mlist, intensities)
             END IF
             IF(.NOT.something_frozen) CALL molden_out(input,particles,H_eigval2,D,intensities,calc_intens,&
                                            dump_only_positive=.FALSE.,logger=logger)
          ELSE
             IF (output_unit>0) THEN
                WRITE(output_unit,'(T2,"VIB| No further vibrational info. Detected a single atom")')
             END IF
          END IF
          ! Deallocate working arrays
          DEALLOCATE(Clist)
          DEALLOCATE(Mlist)
          DEALLOCATE(H_eigval1)
          DEALLOCATE(H_eigval2)
          DEALLOCATE(Hint1)
          DEALLOCATE(Hint2)
          DEALLOCATE(rmass)
          DEALLOCATE(konst)
          DEALLOCATE(mass)
          DEALLOCATE(pos0)
          DEALLOCATE(D)
          DEALLOCATE(Hessian)
          IF(calc_intens)THEN
             DEALLOCATE(dip_deriv)
             DEALLOCATE(tmp_dip)
          END IF
          DEALLOCATE(intensities)
          CALL f_env_rm_defaults(f_env,ierr)
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,print_section,"PROGRAM_RUN_INFO")
    CALL rep_env_release(rep_env)
    CALL timestop(handle)
  END SUBROUTINE vb_anal

! *****************************************************************************
!> \brief give back a list of moving atoms
!> \param force_env ...
!> \param Ilist ...
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE get_moving_atoms(force_env, Ilist)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, DIMENSION(:), POINTER           :: Ilist

    CHARACTER(len=*), PARAMETER :: routineN = 'get_moving_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, ikind, j, &
                                                ndim, nfixed_atoms, &
                                                nfixed_atoms_total, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ifixd_list, work
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
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
    CALL timeset(routineN,handle)
    CALL force_env_get(force_env=force_env, subsys=subsys)

    CALL cp_subsys_get(subsys=subsys,particles=particles, &
                       molecule_kinds_new=molecule_kinds)

    nkind             =  molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els
    particle_set      => particles%els

    ! Count the number of fixed atoms
    nfixed_atoms_total = 0
    DO ikind = 1, nkind
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind, nfixd=nfixed_atoms)
       nfixed_atoms_total = nfixed_atoms_total + nfixed_atoms
    END DO
    ndim = SIZE(particle_set)-nfixed_atoms_total
    CPASSERT(ndim>=0)
    ALLOCATE(Ilist(ndim))

    IF (nfixed_atoms_total/=0) THEN
       ALLOCATE(ifixd_list(nfixed_atoms_total))
       ALLOCATE(work(nfixed_atoms_total))
       nfixed_atoms_total = 0
       DO ikind = 1, nkind
          molecule_kind => molecule_kind_set(ikind)
          CALL get_molecule_kind(molecule_kind, fixd_list=fixd_list)
          IF (ASSOCIATED(fixd_list)) THEN
             DO ii = 1, SIZE(fixd_list)
                IF (.NOT.fixd_list(ii)%restraint%active) THEN
                   nfixed_atoms_total = nfixed_atoms_total + 1
                   ifixd_list(nfixed_atoms_total) = fixd_list(ii)%fixd
                END IF
             END DO
          END IF
       END DO
       CALL sort(ifixd_list,  nfixed_atoms_total, work)

       ndim = 0
       j    = 1
       Loop_count: DO i = 1, SIZE(particle_set)
          DO WHILE (i>ifixd_list(j))
             j = j + 1
             IF (j>nfixed_atoms_total) EXIT Loop_count
          END DO
          IF (i/=ifixd_list(j)) THEN
             ndim = ndim + 1
             Ilist(ndim) = i
          END IF
       END DO Loop_count
       DEALLOCATE(ifixd_list)
       DEALLOCATE(work)
    ELSE
       i    = 1
       ndim = 0
    END IF
    DO j = i, SIZE(particle_set)
       ndim = ndim + 1
       Ilist(ndim) = j
    END DO
    CALL timestop(handle)

  END SUBROUTINE get_moving_atoms

! *****************************************************************************
!> \brief Dumps results of the vibrational analysis
!> \param iw ...
!> \param nvib ...
!> \param D ...
!> \param k ...
!> \param m ...
!> \param freq ...
!> \param particles ...
!> \param Mlist ...
!> \param intensities ...
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
!> \brief Generates the transformation matrix from hessian in cartesian into
!>      internal coordinates (based on Gram-Schmidt orthogonalization)
!> \param mat ...
!> \param dof ...
!> \param Dout ...
!> \param full ...
!> \param natoms ...
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE build_D_matrix(mat,dof,Dout,full,natoms)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat
    INTEGER, INTENT(IN)                      :: dof
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Dout
    LOGICAL, OPTIONAL                        :: full
    INTEGER, INTENT(IN)                      :: natoms

    CHARACTER(len=*), PARAMETER :: routineN = 'build_D_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ifound, iseq, j, &
                                                nvib
    LOGICAL                                  :: failure, my_full
    REAL(KIND=dp)                            :: norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: D

    CALL timeset(routineN,handle)
    failure = .FALSE.
    my_full = .TRUE.
    IF (PRESENT(full)) my_full = full
    ! Generate the missing vectors of the orthogonal basis set
    nvib = 3*natoms-dof
    ALLOCATE(work(3*natoms))
    ALLOCATE(D(3*natoms,3*natoms))
    ! Check First orthogonality in the first element of the basis set
    DO i = 1, dof
       D(:,i) = mat(:,i)
       DO j = i+1, dof
          norm = DOT_PRODUCT(mat(:,i),mat(:,j))
          CPASSERT(ABS(norm)<thrs_motion)
       END DO
    END DO
    ! Generate the nvib orthogonal vectors
    iseq   = 0
    ifound = 0
    DO WHILE (ifound /= nvib)
       iseq = iseq + 1
       CPASSERT(iseq<=3*natoms)
       work       = 0.0_dp
       work(iseq) = 1.0_dp
       ! Gram Schmidt orthogonalization
       DO i = 1, dof+ifound
          norm = DOT_PRODUCT(work,D(:,i))
          work(:) = work - norm * D(:,i)
       END DO
       ! Check norm of the new generated vector
       norm = SQRT(DOT_PRODUCT(work,work))
       IF (norm>=10E4_dp*thrs_motion) THEN
          ! Accept new vector
          ifound = ifound + 1
          D(:,dof+ifound) = work / norm
       END IF
    END DO
    CPASSERT(dof+ifound==3*natoms)
    IF (my_full) THEN
       ALLOCATE(Dout(3*natoms,3*natoms))
       Dout = D
    ELSE
       ALLOCATE(Dout(3*natoms,nvib))
       Dout = D(:,dof+1:)
    END IF
    DEALLOCATE(work)
    DEALLOCATE(D)
    DEALLOCATE(mat)
    CALL timestop(handle)
  END SUBROUTINE build_D_matrix

END MODULE vibrational_analysis

