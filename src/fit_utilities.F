!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart               !
!------------------------------------------------------------------------------!
    MODULE fit_utilities
!-----------------------------------------------------------------------------!
!!****** cp2k/fit_utilities [1.0] *
!!
!!   NAME
!!     fit_utilities
!!
!!   FUNCTION
!!
!!   enables to do different types of fit (rigid_ion, polarizable ..)
!!   using the amoeba procedire
!!
!*****************************************************************
!------------------------------------------------------------------------------!
!
  USE fit_parameters_types,            ONLY: fit_averages_type,&
                                             fit_info_type,&
                                             fit_parameter_type,&
                                             linklist_parameters_type,&
                                             linklist_real_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE md,                              ONLY: simulation_parameters_type
  USE pair_potential,                  ONLY: potentialparm_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv
  USE structure_types,                 ONLY: structure_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE
!
  PUBLIC :: force_ref, fit_print_info, get_will_params, get_lj_params, &
            fit_read_restart, fit_write_restart
  PRIVATE
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE force_ref(struc,potparm,simpar,fit_param,finfo, box_ref)
    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref

    INTEGER                                  :: ii, ipart, isos, it, itime, &
                                                jj, kk, npart
    LOGICAL                                  :: copy
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r, v
    REAL(dbl), DIMENSION(3)                  :: fref
    REAL(dbl), DIMENSION(3, 3)               :: hmat

! From the first principles trajectory file obtain the total force
! on the particle by inverting the velcocity Verlet equations

    npart = size(struc%pnode)
    itime = fit_param%nsteps/fit_param%iskip

    ALLOCATE (box_ref(itime), stat = isos) 
    IF ( isos /= 0 ) CALL stop_memory ( 'force_ref', 'box_ref')

    IF ( fit_param % datatype == 'POS_FORCE' ) THEN
       ALLOCATE ( finfo ( itime, npart ), STAT = isos )
       IF ( isos /= 0 ) CALL stop_memory ( 'force_ref', 'finfo')

    ELSE IF ( fit_param % datatype == 'POS_VEL' ) THEN

! Can only compute force up to nsteps-1 step: thus allocated so
       IF ( MOD ( fit_param % nsteps, itime ) == 0 ) THEN
          ALLOCATE ( finfo ( itime-1, npart ), STAT = isos )
          IF ( isos /= 0 ) CALL stop_memory ( 'force_ref',  'finfo')
       ELSE
          ALLOCATE ( finfo ( itime, npart ), STAT = isos )
          IF ( isos /= 0 ) CALL stop_memory ( 'force_ref',  'finfo')
       END IF
    END IF

    SELECT CASE (fit_param%datatype)
! FOR DATA TYPE FORCES
    CASE ('POS_FORCE')

      ALLOCATE (r(1,1,3),STAT=isos)
      ALLOCATE (v(1,1,3),STAT=isos) !  reads force data
      IF (isos/=0) CALL stop_memory('force_ref','failed to allocate r and v')
      copy = .FALSE.
      it = 0
      DO itime = 1, fit_param%nsteps
        IF (mod(itime,fit_param%iskip)==0) THEN
          copy = .TRUE.
          it = it + 1
        END IF
!  particle position and forces from ab-initio trajectory    
        DO ipart = 1, npart
          READ (250,*) ii, r(1,1,:), v(1,1,:)
          IF (copy) THEN
            finfo(it,ipart) %r(:) = r(1,1,:)
            finfo(it,ipart) %f_ref(:) = v(1,1,:)
          END IF
        END DO
! read the box    
        DO ii=1,3
          READ (250,*) kk, (hmat(ii,jj),jj=1,3)
        ENDDO
        IF (copy) THEN
          box_ref(it)%hmat=hmat
          box_ref(it)%perd=1
          CALL get_hinv(box_ref(it))
        ENDIF
        copy = .FALSE.
      END DO
! Reassign fit_param%nsteps=it
      fit_param%nsteps = it
      DEALLOCATE (r,STAT=isos)
      DEALLOCATE (v,STAT=isos)
      IF (isos/=0) CALL stop_memory('force_ref', 'r and v')
      
! FOR DATA TYPE TRAJECTORY
    CASE ('POS_VEL')
      ALLOCATE (r(2,npart,3),STAT=isos)
      ALLOCATE (v(2,npart,3),STAT=isos) ! reads velocity data
      IF (isos/=0) CALL stop_memory('force_ref','r and v')
! Read in particles from t=1
      DO ipart = 1, npart
        READ (200,*) ii, r(1,ipart,:), v(1,ipart,:)
      END DO
! initialize counter and logical
      copy = .FALSE.
      it = 1
      DO itime = 2, fit_param%nsteps
! set flag
        IF (mod(itime-1,fit_param%iskip)==0) copy = .TRUE.
        DO ipart = 1, npart
          READ (200,*) ii, r(2,ipart,:), v(2,ipart,:)
! Invert velocity verlet equations
          fref(:) = 2._dbl*struc%pnode(ipart) %p%prop%mass* &
            (r(2,ipart,:)-r(1,ipart,:)-simpar%dt*v(1,ipart,:))/ &
            simpar%dt/simpar%dt
! reassigning r(1) to r(2) and v(1) to v(2) for next force calculation
          IF (copy) THEN
            finfo(it,ipart) %r(:) = r(1,ipart,:)
            finfo(it,ipart) %f_ref(:) = fref(:)
          END IF
        END DO
! increment the counter
        IF ((itime/=fit_param%nsteps) .AND. copy) it = it + 1
! reset logical
        copy = .FALSE.
        r(1,:,:) = r(2,:,:)
        v(1,:,:) = v(2,:,:)
      END DO
! Reassign fit_param%nsteps=it-1 
      fit_param%nsteps = it - 1
      DEALLOCATE (r,STAT=isos)
      DEALLOCATE (v,STAT=isos)
      IF (isos/=0) CALL stop_program('force_ref', &
        'failed to deallocate r and v')
    END SELECT
  END SUBROUTINE force_ref
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE fit_print_info(tol,potparm,fit_param,y,atom_names,&
                            fit_ave,globenv,iter)
    REAL(dbl), INTENT(IN)                    :: tol
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    TYPE(fit_parameter_type), INTENT(IN)     :: fit_param
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: y
    CHARACTER(LEN=*), DIMENSION(:), POINTER  :: atom_names
    TYPE(fit_averages_type)                  :: fit_ave
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER                                  :: iter

    INTEGER                                  :: iatom, iw, jatom

    iw = globenv % scr 

    WRITE (iw,*)
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
    WRITE (iw,'(A,T71,A)') ' FIT TYPE                = ', &
      adjustr(fit_param%fit_type)
    WRITE (iw,'(A,T71,A)') ' OUT UNIT                = ', &
      adjustr(fit_param%outunit)
    WRITE (iw,'(A,T71,I10)') ' ITERATION NUMBER        = ', iter
    WRITE (iw,'(A,T66,E15.7)') ' |(F_REF-F_FIT)|^2       = ', &
      sum(y)/size(y)
    WRITE (iw,'(A,T66,E15.7)') ' TOLERANCE               = ', tol
    WRITE (iw,*)
    IF (fit_param%pot_type=='LJ') THEN
      WRITE (iw,'(T5,A,T15,A,T50,A,T73,A)') 'ATOMI', 'ATOMJ', 'EPSILON', &
        'SIGMA'
      WRITE (iw,*)
      DO iatom = 1, size(potparm,1)
        DO jatom = iatom, size(potparm,2)
          WRITE (iw,'(T5,A,T15,A,T45,E15.8,T66,E15.8)') atom_names(iatom), &
            atom_names(jatom), fit_ave% ave_eps(iatom,jatom), &
             fit_ave% ave_sig(iatom,jatom)
        END DO
      END DO
    ELSE IF (fit_param%pot_type=='WILL') THEN
      WRITE (iw,'(T5,A,T15,A,T44,A,T60,A,T75,A)') 'ATOMI', 'ATOMJ', 'A', &
        'B', 'C'
      WRITE (iw,*)
      DO iatom = 1, size(potparm,1)
        DO jatom = iatom, size(potparm,2)
          WRITE (iw,'(T5,A,T15,A,T39,E12.6,T54,E12.6,T69,E12.6)') &
            atom_names(iatom), atom_names(jatom), fit_ave%ave_a(iatom,jatom), &
            fit_ave% ave_b(iatom,jatom), fit_ave %ave_c(iatom,jatom)
        END DO
      END DO
    END IF
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
  END SUBROUTINE fit_print_info
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE get_lj_params(potparm, y, param, ave_sig, ave_eps)
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: y
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(IN)            :: param
    REAL(dbl), DIMENSION(:, :), POINTER      :: ave_sig, ave_eps

    INTEGER                                  :: i, iatom, jatom
    TYPE(linklist_real_type), POINTER        :: lleps, llsig

    ave_sig =  0_dbl
    ave_eps =  0_dbl
    DO iatom = 1, size(potparm,1)
      DO jatom = iatom, size(potparm,2)
        lleps => param(iatom,jatom) %lleps
        llsig => param(iatom,jatom) %llsig
        DO i = 1, size(y)
          ave_sig(iatom,jatom) = ave_sig(iatom,jatom) + &
                                 llsig%real_value*llsig%real_value
          ave_eps(iatom,jatom) = ave_eps(iatom,jatom) + &
                                 lleps%real_value*lleps%real_value
          llsig => llsig%next
          lleps => lleps%next
        END DO
        ave_sig(iatom,jatom) = ave_sig(iatom,jatom)/size(y)
        ave_eps(iatom,jatom) = ave_eps(iatom,jatom)/size(y)
      END DO
    END DO
    ave_sig = ave_sig**(1._dbl/6._dbl)
  END SUBROUTINE get_lj_params
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE get_will_params(potparm, y, param, ave_a, ave_b, ave_c)
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: y
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(IN)            :: param
    REAL(dbl), DIMENSION(:, :), POINTER      :: ave_a, ave_b, ave_c

    INTEGER                                  :: i, iatom, jatom
    TYPE(linklist_real_type), POINTER        :: lla, llb, llc

    ave_a =  0_dbl
    ave_b =  0_dbl
    ave_c =  0_dbl
    DO iatom = 1, size(potparm,1)
      DO jatom = iatom, size(potparm,2)
        lla => param(iatom,jatom) %lla
        llb => param(iatom,jatom) %llb
        llc => param(iatom,jatom) %llc
        DO i = 1, size(y)
          ave_a(iatom,jatom) = ave_a(iatom,jatom) + lla%real_value*lla%real_value 
          ave_b(iatom,jatom) = ave_b(iatom,jatom) + llb%real_value*llb%real_value
          ave_c(iatom,jatom) = ave_c(iatom,jatom) + llc%real_value*llc%real_value
          lla => lla%next
          llb => llb%next
          llc => llc%next
        END DO
        ave_a(iatom,jatom) = ave_a(iatom,jatom)/size(y)
        ave_b(iatom,jatom) = ave_b(iatom,jatom)/size(y)
        ave_c(iatom,jatom) = ave_c(iatom,jatom)/size(y)
      END DO
    END DO
  END SUBROUTINE get_will_params
!******************************************************************************

SUBROUTINE fit_write_restart ( y, p )


    REAL(dbl), DIMENSION(:), INTENT(IN)      :: y
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: p

    INTEGER                                  :: i, iochannel = 99, j

!------------------------------------------------------------------------------
! write to the restart file

  OPEN ( UNIT = iochannel, FILE = 'amoeba.rst', STATUS = 'UNKNOWN' )

  DO i = 1, SIZE ( y )
     WRITE ( iochannel, FMT = * ) y ( i )
  END DO

  DO i = 1, SIZE ( y )
     WRITE ( iochannel, FMT = * ) ( p ( i, j ), j = 1, SIZE ( p, 2 )  )
  END DO

  CLOSE ( UNIT = iochannel )

END SUBROUTINE fit_write_restart

!******************************************************************************

SUBROUTINE fit_read_restart ( y, p )


    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: y
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: p

    INTEGER                                  :: i, iochannel = 99, j

!------------------------------------------------------------------------------
! read from the restart file

  OPEN ( UNIT = iochannel, FILE = 'amoeba.rst', STATUS = 'UNKNOWN' )

  DO i = 1, SIZE ( y )
     READ ( iochannel, FMT = * ) y ( i )
  END DO

  DO i = 1, SIZE ( y )
     READ ( iochannel, FMT = * ) ( p ( i, j ), j = 1, SIZE ( p, 2 )  )
  END DO

  CLOSE ( UNIT = iochannel )

END SUBROUTINE fit_read_restart

!******************************************************************************
!
    END MODULE fit_utilities
!------------------------------------------------------------------------------!
