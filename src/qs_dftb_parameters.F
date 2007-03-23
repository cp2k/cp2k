!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_dftb_parameters [1.0] *
!!
!!   NAME
!!     qs_dftb_parameters
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (27.02.2007)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_dftb_parameters

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_control_types,                ONLY: dftb_control_type
  USE cp_files,                        ONLY: close_file, &
                                             open_file, &
                                             get_unit_number
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE external_potential_types,        ONLY: set_potential
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE message_passing,                 ONLY: mp_bcast
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                             qs_dftb_pairpot_type,&
                                             qs_dftb_pairpot_init,&
                                             qs_dftb_pairpot_create
  USE qs_dftb_utils,                   ONLY: allocate_dftb_atom_param,&
                                             get_dftb_atom_param,&
                                             set_dftb_atom_param


#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_parameters'

  REAL(dp), PARAMETER                  :: slako_d0 = 1._dp

! *** Public subroutines ***

  PUBLIC :: qs_dftb_param_init

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE qs_dftb_param_init (atomic_kind_set,dftb_control,&
                                 force_env_section,para_env,error)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dftb_control_type), INTENT(inout)   :: dftb_control
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_dftb_param_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, found, search
    INTEGER                                  :: k, j, ikind, jkind, nkind, istat,&
       output_unit, runit, ngrd, n_urpoly, z, spdim, isp, l, lp, lmax, &
       lmax_a, lmax_b, l1, l2, m, llm, nskcut
    REAL(dp)                                 :: dgrd, energy, s_cut, radmax, &
                                                zeff, ra, rb
    REAL(dp),DIMENSION(10)                   :: uwork, fwork, swork
    REAL(dp),DIMENSION(0:3)                  :: skself, eta, occupation
    REAL(dp),DIMENSION(1:3)                  :: srep
    REAL(dp),DIMENSION(1:2)                  :: surr
    CHARACTER(LEN=default_string_length)     :: iname, jname, name_a, name_b
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_a, atomic_kind_b
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_atom_a, dftb_atom_b
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:,:), POINTER                :: dftb_potential
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(LEN=6)                         :: cspline
    REAL(dp),DIMENSION(:,:),ALLOCATABLE      :: fmat, smat, spxr, scoeff

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:,:), ALLOCATABLE            :: sk_files

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    output_unit = -1
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
         "SUBSYS%PRINT%KINDS/BASIS_SET",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,force_env_section,&
            "SUBSYS%PRINT%KINDS",extension=".Log",error=error)

    END IF

    nkind  = SIZE(atomic_kind_set)
    ALLOCATE(sk_files(nkind,nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ! allocate potential structures
    ALLOCATE(dftb_control%dftb_potential(nkind,nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    dftb_potential => dftb_control%dftb_potential
    CALL qs_dftb_pairpot_init(dftb_potential,error)

    DO ikind = 1, nkind
      atomic_kind_a => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind_a,name=iname)
      DO jkind = 1, ikind
        atomic_kind_b => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind_b,name=jname)
        found = .FALSE.
        DO k=1,SIZE(dftb_control%sk_pair_list,2)
          name_a=TRIM(dftb_control%sk_pair_list(1,k))
          name_b=TRIM(dftb_control%sk_pair_list(2,k))
          IF ( (iname==name_a .AND. jname==name_b) .OR. &
               (iname==name_b .AND. jname==name_a) ) THEN
            sk_files(ikind,jkind) = TRIM(dftb_control%sk_file_path)//"/"//&
                                    TRIM(dftb_control%sk_pair_list(3,k))
            sk_files(jkind,ikind) = TRIM(dftb_control%sk_file_path)//"/"//&
                                    TRIM(dftb_control%sk_pair_list(3,k))
            found = .TRUE.
            EXIT
          END IF
        END DO
        CPPostcondition(found,cp_failure_level,routineP,error,failure)
      END DO
    END DO

    ! reading the files
    IF (output_unit > 0) THEN
      WRITE(output_unit,"(/,A)") " DFTB| Reading parameter files "
    END IF
    DO ikind = 1, nkind
      atomic_kind_a => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind_a,dftb_parameter=dftb_atom_a,&
                           z=z,name=iname)

      IF (.NOT.ASSOCIATED(dftb_atom_a)) THEN
        CALL allocate_dftb_atom_param(atomic_kind_a%dftb_parameter,error=error)
        CALL get_atomic_kind(atomic_kind_a,dftb_parameter=dftb_atom_a)
      END IF

      DO jkind = ikind,1,-1
        atomic_kind_b => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind_b,dftb_parameter=dftb_atom_b,&
                             name=jname)

        IF (output_unit > 0) THEN
          WRITE(output_unit,"(A,A)") " DFTB| Reading parameter file ",sk_files(jkind,ikind)
        END IF
        skself = 0._dp
        eta    = 0._dp
        occupation = 0._dp
        IF ( para_env%ionode ) THEN
          runit = get_unit_number()
          CALL open_file(file_name=sk_files(jkind,ikind),unit_number=runit)
          ! grid density and number of grid poin ts
          READ(runit,fmt=*,end=1,err=1) dgrd,ngrd
!
! ngrd -1 ?
! In Slako tables, the grid starts at 0.0, in deMon it starts with dgrd
!
          ngrd = ngrd - 1
!
          IF ( ikind == jkind) THEN
            ! orbital energy, total energy, hardness, occupation
            READ (runit,fmt=*,end=1,err=1) skself(2:0:-1),energy, &
                                           eta(2:0:-1),occupation(2:0:-1)
          END IF
          ! repulsive potential as polynomial
          READ (runit,fmt=*,end=1,err=1) uwork(1:10)
          n_urpoly = 0
          IF ( DOT_PRODUCT(uwork(2:10),uwork(2:10)) >= TINY(1._dp)) THEN
            n_urpoly = 1
            DO k=2,10
              IF ( ABS(uwork(k)) >= TINY(1._dp) ) n_urpoly = k
            END DO
          END IF
        END IF

        CALL mp_bcast(n_urpoly,para_env%source,para_env%group)
        CALL mp_bcast(uwork,para_env%source,para_env%group)
        CALL mp_bcast(ngrd,para_env%source,para_env%group)
        CALL mp_bcast(dgrd,para_env%source,para_env%group)

        IF ( ikind == jkind) THEN
          CALL mp_bcast(skself,para_env%source,para_env%group)
          CALL mp_bcast(energy,para_env%source,para_env%group)
          CALL mp_bcast(eta,para_env%source,para_env%group)
          CALL mp_bcast(occupation,para_env%source,para_env%group)
  
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,&
               z=z,zeff=SUM(occupation),defined=.TRUE.,&
               skself=skself,energy=energy,eta=eta,occupation=occupation,&
               error=error)
        END IF

        ! Slater-Koster table
        ALLOCATE(fmat(ngrd,10),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(smat(ngrd,10),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        IF ( para_env%ionode ) THEN
          DO k=1,ngrd
            READ (runit,fmt=*,end=1,err=1) fwork(1:10),swork(1:10)
            fmat(k,1:10) = fwork(1:10)
            smat(k,1:10) = swork(1:10)
          END DO
        END IF
        CALL mp_bcast(fmat,para_env%source,para_env%group)
        CALL mp_bcast(smat,para_env%source,para_env%group)

        IF ( ikind == jkind ) THEN
          !
          ! Determine lmax for atom type.
          ! An atomic orbital is 'active' if either its onsite energy is different from zero, 
          ! or
          ! if this matrix element contains non-zero elements. 
          ! The sigma interactions are sufficient for that.
          ! In the DFTB-Slako convention they are on orbital 10 (s-s-sigma), 
          ! 7 (p-p-sigma) and 3 (d-d-sigma).
          !
          lmax=0
          DO l=0,3
            SELECT CASE (l)
              CASE DEFAULT
                stop
              CASE (0)
                lp = 10
              CASE (1)
                lp = 7
              CASE (2)
                lp = 3
              CASE (3)
                lp = 3 ! this is wrong but we don't allow f anyway 
            END SELECT
            ! Technical note: In some slako files dummies are included in the 
            ! first matrix elements, so remove them.
!? in the comment above it says .OR.
            IF ( (ABS(skself(l)) > 0._dp) .AND. &
                 (SUM(ABS(fmat(ngrd/10:ngrd,lp))) > 0._dp) ) lmax=l
          END DO
          ! l=2 (d) is maximum
          lmax = MIN ( 2, lmax )
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,&
             lmax=lmax, natorb=(lmax+1)**2, error=error)

        END IF

        spdim = 0
        IF ( n_urpoly == 0 ) THEN
          IF ( para_env%ionode ) THEN
            ! Look for spline representation of repulsive potential
            search = .true.
            DO WHILE (search)
              READ (runit,fmt='(A6)',end=1,err=1) cspline
              IF (cspline == 'Spline') THEN
                search = .false.
                ! spline dimension and left-hand cutoff
                READ (runit,fmt=*,end=1,err=1) spdim,s_cut
                ALLOCATE(spxr(spdim,2),STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(scoeff(spdim,4),STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                ! e-functions describing left-hand extrapolation
                READ (runit,fmt=*,end=1,err=1) srep(1:3)
                DO isp = 1,spdim-1
                   ! location and coefficients of 'normal' spline range
                   READ (runit,fmt=*,end=1,err=1) spxr(isp,1:2),scoeff(isp,1:4)
                END DO
                ! last point has 2 more coefficients
                READ (runit,fmt=*,end=1,err=1) spxr(spdim,1:2),scoeff(spdim,1:4),surr(1:2)
              END IF
            END DO
          END IF
        END IF

        IF ( para_env%ionode ) THEN
          CALL close_file(unit_number=runit)
        END IF

        CALL mp_bcast(spdim,para_env%source,para_env%group)
        IF ( spdim > 0 .AND. (.NOT. para_env%ionode)) THEN
          ALLOCATE(spxr(spdim,2),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(scoeff(spdim,4),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(spdim > 0) THEN
          CALL mp_bcast(spxr,para_env%source,para_env%group)
          CALL mp_bcast(scoeff,para_env%source,para_env%group)
          CALL mp_bcast(surr,para_env%source,para_env%group)
          CALL mp_bcast(srep,para_env%source,para_env%group)
          CALL mp_bcast(s_cut,para_env%source,para_env%group)
        END IF

        ! store potential data
        ! allocate data
        CALL get_dftb_atom_param(dftb_parameter=dftb_atom_a,lmax=lmax_a)
        CALL get_dftb_atom_param(dftb_parameter=dftb_atom_b,lmax=lmax_b)
        llm = 0
        do l1 = 0,MAX(lmax_a,lmax_b)
          do l2 = 0,MIN(l1,lmax_a,lmax_b)
            do m = 0,l2
              llm = llm + 1
            end do
          end do
        end do
        CALL qs_dftb_pairpot_create(dftb_potential(ikind,jkind),&
                                    ngrd,llm,spdim,error)

        ! repulsive potential
        dftb_potential(ikind,jkind)%n_urpoly = n_urpoly
        dftb_potential(ikind,jkind)%urep_cut = uwork(10)
        dftb_potential(ikind,jkind)%urep(1:n_urpoly-1) = uwork(2:n_urpoly)

        ! Slater-Koster tables
        dftb_potential(ikind,jkind)%dgrd = dgrd
        CALL skreorder(fmat,lmax_a,lmax_b)
        dftb_potential(ikind,jkind)%fmat(:,1:llm) = fmat(:,1:llm)
        CALL skreorder(smat,lmax_a,lmax_b)
        dftb_potential(ikind,jkind)%smat(:,1:llm) = smat(:,1:llm)
! how to do this correctly???????????????
        dftb_potential(ikind,jkind)%ngrdcut = ngrd + slako_d0/dgrd
        ! Splines
        IF(spdim > 0) THEN
          dftb_potential(ikind,jkind)%s_cut = s_cut
          dftb_potential(ikind,jkind)%srep = srep
          dftb_potential(ikind,jkind)%spxr = spxr
          dftb_potential(ikind,jkind)%scoeff = scoeff
          dftb_potential(ikind,jkind)%surr = surr
        END IF

        DEALLOCATE(fmat,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(smat,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        IF(spdim > 0) THEN
          DEALLOCATE(spxr,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(scoeff,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

      END DO
    END DO

    DEALLOCATE(sk_files,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! extract simple atom interaction radii
    DO ikind = 1, nkind
      atomic_kind_a => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind_a,dftb_parameter=dftb_atom_a)
      radmax = (dftb_potential(ikind,ikind)%ngrdcut + 1) * &
               dftb_potential(ikind,ikind)%dgrd*0.5_dp
      CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,cutoff=radmax,error=error)
    END DO
    DO ikind = 1, nkind
      atomic_kind_a => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind_a,dftb_parameter=dftb_atom_a)
      CALL get_dftb_atom_param(dftb_parameter=dftb_atom_a,cutoff=ra,error=error)
      DO jkind = 1, ikind
        atomic_kind_b => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind_b,dftb_parameter=dftb_atom_b)
        CALL get_dftb_atom_param(dftb_parameter=dftb_atom_b,cutoff=rb,error=error)
        radmax = (dftb_potential(ikind,jkind)%ngrdcut + 1) * &
                 dftb_potential(ikind,jkind)%dgrd
        IF ( ra+rb < radmax ) THEN
          ra = ra + (radmax-ra-rb)*0.5_dp
          rb = rb + (radmax-ra-rb)*0.5_dp
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,cutoff=ra,error=error)
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_b,cutoff=rb,error=error)
        END IF
      END DO
    END DO

    ! set correct core charge in potential
    DO ikind = 1, nkind
      atomic_kind_a => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind_a,dftb_parameter=dftb_atom_a)
      CALL get_dftb_atom_param(dftb_parameter=dftb_atom_a,zeff=zeff,error=error)
      CALL set_potential(potential=atomic_kind_a%all_potential,&
                  zeff=zeff,zeff_correction=0.0_dp)
    END DO
    
    RETURN

1   CONTINUE
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)

  END SUBROUTINE qs_dftb_param_init

! *****************************************************************************
  SUBROUTINE skreorder(xmat,la,lb)
      REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: xmat
      INTEGER, INTENT(IN)                     :: la,lb

      REAL(dp)                :: skllm(0:3,0:3,0:3)
      INTEGER                 :: i,l1,l2,llm,m
      !
      ! Purpose: Transform Slako format in l1/l2/m format
      !
      ! Slako tables from Dresden/Paderborn/Heidelberg groups are
      ! stored in the following native format:
      !
      ! Convention: Higher angular momenta are always on the right-hand side
      !
      ! 1: d - d - delta
      ! 2: d - d - pi
      ! 3: d - d - sigma
      ! 4: p - d - pi
      ! 5: p - d - sigma
      ! 6: p - p - pi
      ! 7: p - p - sigma
      ! 8: d - s - sigma
      ! 9: p - s - sigma
      !10: s - s - sigma
      !
      !
      DO i=1,SIZE(xmat,1)
        skllm = 0._dp
        skllm(0,0,0) = xmat(i,10)
        skllm(1,0,0) = xmat(i,9)
        skllm(2,0,0) = xmat(i,8)
        skllm(1,1,1) = xmat(i,7)
        skllm(1,1,0) = xmat(i,6)
        skllm(2,1,1) = xmat(i,5)
        skllm(2,1,0) = xmat(i,4)
        skllm(2,2,2) = xmat(i,3)
        skllm(2,2,1) = xmat(i,2)
        skllm(2,2,0) = xmat(i,1)
        llm=0
        DO l1=0,MAX(la,lb)
          DO l2=0,MIN(l1,la,lb)
            DO m=0,l2
              llm=llm+1
              xmat(i,llm)=skllm(l1,l2,m)
            END DO
          END DO
        END DO
      END DO
      !
   END SUBROUTINE skreorder

! *****************************************************************************
  
END MODULE qs_dftb_parameters

! *****************************************************************************
