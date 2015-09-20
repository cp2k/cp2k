!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author JGH (27.02.2007)
! *****************************************************************************
MODULE qs_dftb_parameters

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_control_types,                ONLY: dftb_control_type
  USE cp_files,                        ONLY: close_file,&
                                             get_unit_number,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE external_potential_types,        ONLY: set_potential
  USE input_constants,                 ONLY: dispersion_uff
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_bcast
  USE physcon,                         ONLY: angstrom,&
                                             kcalmol
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                             qs_dftb_pairpot_create,&
                                             qs_dftb_pairpot_init,&
                                             qs_dftb_pairpot_type
  USE qs_dftb_utils,                   ONLY: allocate_dftb_atom_param,&
                                             get_dftb_atom_param,&
                                             set_dftb_atom_param
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type,&
                                             set_qs_kind
  USE string_utilities,                ONLY: uppercase
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_parameters'

  REAL(dp), PARAMETER                  :: slako_d0 = 1._dp

! *** Public subroutines ***

  PUBLIC :: qs_dftb_param_init

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param dftb_control ...
!> \param dftb_potential ...
!> \param subsys_section ...
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE qs_dftb_param_init (atomic_kind_set,qs_kind_set,dftb_control,dftb_potential,&
                                 subsys_section,para_env)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(dftb_control_type), INTENT(inout)   :: dftb_control
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), POINTER               :: dftb_potential
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_dftb_param_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: iel, jel
    CHARACTER(LEN=6)                         :: cspline
    CHARACTER(LEN=default_path_length)       :: file_name
    CHARACTER(LEN=default_path_length), &
      ALLOCATABLE, DIMENSION(:, :)           :: sk_files
    CHARACTER(LEN=default_string_length)     :: iname, jname, name_a, name_b, &
                                                skfn
    INTEGER :: ikind, isp, istat, jkind, k, l, l1, l2, llm, lmax, lmax_a, &
      lmax_b, lp, m, n_urpoly, ngrd, nkind, output_unit, runit, spdim, z
    LOGICAL                                  :: at_end, failure, found, ldum, &
                                                search, sklist
    REAL(dp)                                 :: da, db, dgrd, dij, energy, &
                                                eps_disp, ra, radmax, rb, &
                                                rcdisp, rmax6, s_cut, xij, &
                                                zeff
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: fmat, scoeff, smat, spxr
    REAL(dp), DIMENSION(0:3)                 :: eta, occupation, skself
    REAL(dp), DIMENSION(10)                  :: fwork, swork, uwork
    REAL(dp), DIMENSION(1:2)                 :: surr
    REAL(dp), DIMENSION(1:3)                 :: srep
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_atom_a, dftb_atom_b

    failure = .FALSE.
    output_unit = -1
    NULLIFY(logger)
    logger => cp_get_default_logger()
    IF (BTEST(cp_print_key_should_output(logger%iter_info,subsys_section,&
         "PRINT%KINDS/BASIS_SET"),cp_p_file)) THEN
       output_unit = cp_print_key_unit_nr(logger,subsys_section,&
            "PRINT%KINDS",extension=".Log")
       IF ( output_unit > 0 ) THEN
         WRITE(output_unit,"(/,A)") " DFTB| A set of relativistic DFTB "//&
              "parameters for material sciences."
         WRITE(output_unit,"(A)") " DFTB| J. Frenzel, N. Jardillier, A.F. Oliveira,"//&
              " T. Heine, G. Seifert"
         WRITE(output_unit,"(A)") " DFTB| TU Dresden, 2002-2007"
         WRITE(output_unit,"(/,A)") " DFTB| Non-SCC parameters "
         WRITE(output_unit,"(A,T25,A)") " DFTB| C,H         :",&
            " D. Porezag et al, PRB 51 12947 (1995)"
         WRITE(output_unit,"(A,T25,A)") " DFTB| B,N         :",&
            " J. Widany et al, PRB 53 4443 (1996)"
         WRITE(output_unit,"(A,T25,A)") " DFTB| Li,Na,K,Cl  :",&
            " S. Hazebroucq et al, JCP 123 134510 (2005)"
         WRITE(output_unit,"(A,T25,A)") " DFTB| F           :",&
            " T. Heine et al, JCSoc-Perkins Trans 2 707 (1999)"
         WRITE(output_unit,"(A,T25,A)") " DFTB| Mo,S        :",&
            " G. Seifert et al, PRL 85 146 (2000)"
         WRITE(output_unit,"(A,T25,A)") " DFTB| P           :",&
            " G. Seifert et al, EPS 16 341 (2001)"
         WRITE(output_unit,"(A,T25,A)") " DFTB| Sc,N,C      :",&
            " M. Krause et al, JCP 115 6596 (2001)"
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
            "PRINT%KINDS")
    END IF

    sklist = (dftb_control%sk_file_list /= "")

    nkind  = SIZE(atomic_kind_set)
    ALLOCATE(sk_files(nkind,nkind),STAT=istat)
    CPASSERT(istat==0)
    ! allocate potential structures
    ALLOCATE(dftb_potential(nkind,nkind),STAT=istat)
    CPASSERT(istat==0)
    CALL qs_dftb_pairpot_init(dftb_potential)

    DO ikind = 1, nkind
      CALL get_atomic_kind(atomic_kind_set(ikind),name=iname,element_symbol=iel)
      CALL uppercase(iname)
      CALL uppercase(iel)
      ldum = qmmm_ff_precond_only_qm(iname)
      DO jkind = 1, nkind
        CALL get_atomic_kind(atomic_kind_set(jkind),name=jname,element_symbol=jel)
        CALL uppercase(jname)
        CALL uppercase(jel)
        ldum = qmmm_ff_precond_only_qm(jname)
        found = .FALSE.
        DO k=1,SIZE(dftb_control%sk_pair_list,2)
          name_a=TRIM(dftb_control%sk_pair_list(1,k))
          name_b=TRIM(dftb_control%sk_pair_list(2,k))
          CALL uppercase(name_a)
          CALL uppercase(name_b)
          IF ( (iname==name_a .AND. jname==name_b) ) THEN
            sk_files(ikind,jkind) = TRIM(dftb_control%sk_file_path)//"/"//&
                                    TRIM(dftb_control%sk_pair_list(3,k))
            found = .TRUE.
            EXIT
          END IF
        END DO
        IF ( .NOT. found .AND. sklist ) THEN
          file_name = TRIM(dftb_control%sk_file_path)//"/"//&
                      TRIM(dftb_control%sk_file_list)
          NULLIFY(parser)
          CALL parser_create(parser,file_name,para_env=para_env)
          DO
            at_end = .FALSE.
            CALL parser_get_next_line(parser,1,at_end)
            IF ( at_end ) EXIT
            CALL parser_get_object(parser,name_a,lower_to_upper=.TRUE.)
            CALL parser_get_object(parser,name_b,lower_to_upper=.TRUE.)
            !Checking Names
            IF ( (iname==name_a .AND. jname==name_b) ) THEN
              CALL parser_get_object(parser,skfn,string_length=8)
              sk_files(ikind,jkind) = TRIM(dftb_control%sk_file_path)//"/"//&
                                      TRIM(skfn)
              found = .TRUE.
              EXIT
            END IF
            !Checking Element
            IF ( (iel==name_a .AND. jel==name_b) ) THEN
              CALL parser_get_object(parser,skfn,string_length=8)
              sk_files(ikind,jkind) = TRIM(dftb_control%sk_file_path)//"/"//&
                                      TRIM(skfn)
              found = .TRUE.
              EXIT
            END IF
          END DO
          CALL parser_release(parser)
        END IF
        CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
             "Failure in assigning KINDS <"//TRIM(iname)//"> and <"//TRIM(jname)//&
             "> to a DFTB interaction pair!"//&
CPSourceFileRef,&
             only_ionode=.TRUE.)
      END DO
    END DO
    ! reading the files
    ! read all pairs, equal kind first
    DO ikind = 1, nkind
      CALL get_atomic_kind(atomic_kind_set(ikind), z=z,name=iname)

      CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)
      IF (.NOT.ASSOCIATED(dftb_atom_a)) THEN
        CALL allocate_dftb_atom_param(dftb_atom_a)
        CALL set_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)
      END IF

      ! read all pairs, equal kind first
      jkind = ikind

      CALL get_atomic_kind(atomic_kind_set(jkind), name=jname)
      CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_atom_b)

      IF (output_unit > 0) THEN
        WRITE(output_unit,"(A,T30,A50)") " DFTB| Reading parameter file ",&
           ADJUSTR(TRIM(sk_files(jkind,ikind)))
      END IF
      skself = 0._dp
      eta    = 0._dp
      occupation = 0._dp
      IF ( para_env%ionode ) THEN
        runit = get_unit_number()
        CALL open_file(file_name=sk_files(jkind,ikind),unit_number=runit)
        ! grid density and number of grid poin ts
        READ(runit,fmt=*,END=1,err=1) dgrd,ngrd
!
! ngrd -1 ?
! In Slako tables, the grid starts at 0.0, in deMon it starts with dgrd
!
        ngrd = ngrd - 1
!
        ! orbital energy, total energy, hardness, occupation
        READ (runit,fmt=*,END=1,err=1) skself(2:0:-1),energy, &
                                       eta(2:0:-1),occupation(2:0:-1)
        ! repulsive potential as polynomial
        READ (runit,fmt=*,END=1,err=1) uwork(1:10)
        n_urpoly = 0
        IF ( DOT_PRODUCT(uwork(2:10),uwork(2:10)) >= 1.e-12_dp) THEN
          n_urpoly = 1
          DO k=2,9
            IF ( ABS(uwork(k)) >= 1.e-12_dp ) n_urpoly = k
          END DO
        END IF
      END IF

      CALL mp_bcast(n_urpoly,para_env%source,para_env%group)
      CALL mp_bcast(uwork,para_env%source,para_env%group)
      CALL mp_bcast(ngrd,para_env%source,para_env%group)
      CALL mp_bcast(dgrd,para_env%source,para_env%group)

      CALL mp_bcast(skself,para_env%source,para_env%group)
      CALL mp_bcast(energy,para_env%source,para_env%group)
      CALL mp_bcast(eta,para_env%source,para_env%group)
      CALL mp_bcast(occupation,para_env%source,para_env%group)

      CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,&
           z=z,zeff=SUM(occupation),defined=.TRUE.,&
           skself=skself,energy=energy,eta=eta,occupation=occupation)

      ! Slater-Koster table
      ALLOCATE(fmat(ngrd,10),STAT=istat)
      CPASSERT(istat==0)
      ALLOCATE(smat(ngrd,10),STAT=istat)
      CPASSERT(istat==0)
      IF ( para_env%ionode ) THEN
        DO k=1,ngrd
          READ (runit,fmt=*,END=1,err=1) fwork(1:10),swork(1:10)
          fmat(k,1:10) = fwork(1:10)
          smat(k,1:10) = swork(1:10)
        END DO
      END IF
      CALL mp_bcast(fmat,para_env%source,para_env%group)
      CALL mp_bcast(smat,para_env%source,para_env%group)

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
            STOP
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
        IF ( (ABS(skself(l)) > 0._dp) .OR. &
             (SUM(ABS(fmat(ngrd/10:ngrd,lp))) > 0._dp) ) lmax=l
      END DO
      ! l=2 (d) is maximum
      lmax = MIN ( 2, lmax )
      CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,&
         lmax=lmax, natorb=(lmax+1)**2)

      spdim = 0
      IF ( n_urpoly == 0 ) THEN
        IF ( para_env%ionode ) THEN
          ! Look for spline representation of repulsive potential
          search = .TRUE.
          DO WHILE (search)
            READ (runit,fmt='(A6)',END=1,err=1) cspline
            IF (cspline == 'Spline') THEN
              search = .FALSE.
              ! spline dimension and left-hand cutoff
              READ (runit,fmt=*,END=1,err=1) spdim,s_cut
              ALLOCATE(spxr(spdim,2),STAT=istat)
              CPASSERT(istat==0)
              ALLOCATE(scoeff(spdim,4),STAT=istat)
              CPASSERT(istat==0)
              ! e-functions describing left-hand extrapolation
              READ (runit,fmt=*,END=1,err=1) srep(1:3)
              DO isp = 1,spdim-1
                 ! location and coefficients of 'normal' spline range
                 READ (runit,fmt=*,END=1,err=1) spxr(isp,1:2),scoeff(isp,1:4)
              END DO
              ! last point has 2 more coefficients
              READ (runit,fmt=*,END=1,err=1) spxr(spdim,1:2),scoeff(spdim,1:4),surr(1:2)
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
        CPASSERT(istat==0)
        ALLOCATE(scoeff(spdim,4),STAT=istat)
        CPASSERT(istat==0)
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
      DO l1 = 0,MAX(lmax_a,lmax_b)
        DO l2 = 0,MIN(l1,lmax_a,lmax_b)
          DO m = 0,l2
            llm = llm + 1
          END DO
        END DO
      END DO
      CALL qs_dftb_pairpot_create(dftb_potential(ikind,jkind),&
                                  ngrd,llm,spdim)

      ! repulsive potential
      dftb_potential(ikind,jkind)%n_urpoly = n_urpoly
      dftb_potential(ikind,jkind)%urep_cut = uwork(10)
      dftb_potential(ikind,jkind)%urep(:) = 0._dp
      dftb_potential(ikind,jkind)%urep(1) = uwork(10)
      dftb_potential(ikind,jkind)%urep(2:n_urpoly) = uwork(2:n_urpoly)

      ! Slater-Koster tables
      dftb_potential(ikind,jkind)%dgrd = dgrd
      CALL skreorder(fmat,lmax_a,lmax_b)
      dftb_potential(ikind,jkind)%fmat(:,1:llm) = fmat(:,1:llm)
      CALL skreorder(smat,lmax_a,lmax_b)
      dftb_potential(ikind,jkind)%smat(:,1:llm) = smat(:,1:llm)
      dftb_potential(ikind,jkind)%ngrdcut = ngrd + INT(slako_d0/dgrd)
      ! Splines
      IF(spdim > 0) THEN
        dftb_potential(ikind,jkind)%s_cut = s_cut
        dftb_potential(ikind,jkind)%srep = srep
        dftb_potential(ikind,jkind)%spxr = spxr
        dftb_potential(ikind,jkind)%scoeff = scoeff
        dftb_potential(ikind,jkind)%surr = surr
      END IF

      DEALLOCATE(fmat,STAT=istat)
      CPASSERT(istat==0)
      DEALLOCATE(smat,STAT=istat)
      CPASSERT(istat==0)
      IF(spdim > 0) THEN
        DEALLOCATE(spxr,STAT=istat)
        CPASSERT(istat==0)
        DEALLOCATE(scoeff,STAT=istat)
        CPASSERT(istat==0)
      END IF

    END DO

    ! no all other pairs
    DO ikind = 1, nkind
      CALL get_atomic_kind(atomic_kind_set(ikind), z=z,name=iname)
      CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)

      IF (.NOT.ASSOCIATED(dftb_atom_a)) THEN
        CALL allocate_dftb_atom_param(dftb_atom_a)
        CALL set_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)
      END IF

      DO jkind = 1, nkind

        IF ( ikind == jkind ) CYCLE
        CALL get_atomic_kind(atomic_kind_set(jkind), name=jname)
        CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_atom_b)

        IF (output_unit > 0) THEN
          WRITE(output_unit,"(A,T30,A50)") " DFTB| Reading parameter file ",&
           ADJUSTR(TRIM(sk_files(ikind,jkind)))
        END IF
        skself = 0._dp
        eta    = 0._dp
        occupation = 0._dp
        IF ( para_env%ionode ) THEN
          runit = get_unit_number()
          CALL open_file(file_name=sk_files(ikind,jkind),unit_number=runit)
          ! grid density and number of grid poin ts
          READ(runit,fmt=*,END=1,err=1) dgrd,ngrd
!
! ngrd -1 ?
! In Slako tables, the grid starts at 0.0, in deMon it starts with dgrd
!
          ngrd = ngrd - 1
!
          IF ( ikind == jkind) THEN
            ! orbital energy, total energy, hardness, occupation
            READ (runit,fmt=*,END=1,err=1) skself(2:0:-1),energy, &
                                           eta(2:0:-1),occupation(2:0:-1)
          END IF
          ! repulsive potential as polynomial
          READ (runit,fmt=*,END=1,err=1) uwork(1:10)
          n_urpoly = 0
          IF ( DOT_PRODUCT(uwork(2:10),uwork(2:10)) >= 1.e-12_dp) THEN
            n_urpoly = 1
            DO k=2,9
              IF ( ABS(uwork(k)) >= 1.e-12_dp ) n_urpoly = k
            END DO
          END IF
        END IF

        CALL mp_bcast(n_urpoly,para_env%source,para_env%group)
        CALL mp_bcast(uwork,para_env%source,para_env%group)
        CALL mp_bcast(ngrd,para_env%source,para_env%group)
        CALL mp_bcast(dgrd,para_env%source,para_env%group)

        ! Slater-Koster table
        ALLOCATE(fmat(ngrd,10),STAT=istat)
        CPASSERT(istat==0)
        ALLOCATE(smat(ngrd,10),STAT=istat)
        CPASSERT(istat==0)
        IF ( para_env%ionode ) THEN
          DO k=1,ngrd
            READ (runit,fmt=*,END=1,err=1) fwork(1:10),swork(1:10)
            fmat(k,1:10) = fwork(1:10)
            smat(k,1:10) = swork(1:10)
          END DO
        END IF
        CALL mp_bcast(fmat,para_env%source,para_env%group)
        CALL mp_bcast(smat,para_env%source,para_env%group)

        spdim = 0
        IF ( n_urpoly == 0 ) THEN
          IF ( para_env%ionode ) THEN
            ! Look for spline representation of repulsive potential
            search = .TRUE.
            DO WHILE (search)
              READ (runit,fmt='(A6)',END=1,err=1) cspline
              IF (cspline == 'Spline') THEN
                search = .FALSE.
                ! spline dimension and left-hand cutoff
                READ (runit,fmt=*,END=1,err=1) spdim,s_cut
                ALLOCATE(spxr(spdim,2),STAT=istat)
                CPASSERT(istat==0)
                ALLOCATE(scoeff(spdim,4),STAT=istat)
                CPASSERT(istat==0)
                ! e-functions describing left-hand extrapolation
                READ (runit,fmt=*,END=1,err=1) srep(1:3)
                DO isp = 1,spdim-1
                   ! location and coefficients of 'normal' spline range
                   READ (runit,fmt=*,END=1,err=1) spxr(isp,1:2),scoeff(isp,1:4)
                END DO
                ! last point has 2 more coefficients
                READ (runit,fmt=*,END=1,err=1) spxr(spdim,1:2),scoeff(spdim,1:4),surr(1:2)
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
          CPASSERT(istat==0)
          ALLOCATE(scoeff(spdim,4),STAT=istat)
          CPASSERT(istat==0)
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
        DO l1 = 0,MAX(lmax_a,lmax_b)
          DO l2 = 0,MIN(l1,lmax_a,lmax_b)
            DO m = 0,l2
              llm = llm + 1
            END DO
          END DO
        END DO
        CALL qs_dftb_pairpot_create(dftb_potential(ikind,jkind),&
                                    ngrd,llm,spdim)

        ! repulsive potential
        dftb_potential(ikind,jkind)%n_urpoly = n_urpoly
        dftb_potential(ikind,jkind)%urep_cut = uwork(10)
        dftb_potential(ikind,jkind)%urep(:) = 0._dp
        dftb_potential(ikind,jkind)%urep(1) = uwork(10)
        dftb_potential(ikind,jkind)%urep(2:n_urpoly) = uwork(2:n_urpoly)

        ! Slater-Koster tables
        dftb_potential(ikind,jkind)%dgrd = dgrd
        CALL skreorder(fmat,lmax_a,lmax_b)
        dftb_potential(ikind,jkind)%fmat(:,1:llm) = fmat(:,1:llm)
        CALL skreorder(smat,lmax_a,lmax_b)
        dftb_potential(ikind,jkind)%smat(:,1:llm) = smat(:,1:llm)
        dftb_potential(ikind,jkind)%ngrdcut = ngrd +INT(slako_d0/dgrd)
        ! Splines
        IF(spdim > 0) THEN
          dftb_potential(ikind,jkind)%s_cut = s_cut
          dftb_potential(ikind,jkind)%srep = srep
          dftb_potential(ikind,jkind)%spxr = spxr
          dftb_potential(ikind,jkind)%scoeff = scoeff
          dftb_potential(ikind,jkind)%surr = surr
        END IF

        DEALLOCATE(fmat,STAT=istat)
        CPASSERT(istat==0)
        DEALLOCATE(smat,STAT=istat)
        CPASSERT(istat==0)
        IF(spdim > 0) THEN
          DEALLOCATE(spxr,STAT=istat)
          CPASSERT(istat==0)
          DEALLOCATE(scoeff,STAT=istat)
          CPASSERT(istat==0)
        END IF

      END DO
    END DO

    DEALLOCATE(sk_files,STAT=istat)
    CPASSERT(istat==0)

    ! read dispersion parameters (UFF type)
    IF ( dftb_control%dispersion ) THEN

      IF ( dftb_control%dispersion_type == dispersion_uff ) THEN
        file_name = TRIM(dftb_control%sk_file_path)//"/"//&
                    TRIM(dftb_control%uff_force_field)
        DO ikind = 1, nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), name=iname)
          CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)

          m = LEN_TRIM(iname)
          NULLIFY(parser)
          CALL parser_create(parser,file_name,para_env=para_env)
          found = .FALSE.
          DO
            at_end = .FALSE.
            CALL parser_get_next_line(parser,1,at_end)
            IF ( at_end ) EXIT
            CALL parser_get_object(parser,name_a)
            ! parser is no longer removing leading quotes
            IF(name_a(1:1) == '"') name_a(1:m) = name_a(2:m+1)
            IF ( name_a(1:m) == TRIM(iname) ) THEN
              CALL parser_get_object(parser,rb)
              CALL parser_get_object(parser,rb)
              CALL parser_get_object(parser,ra)
              CALL parser_get_object(parser,da)
              found = .TRUE.
              ra = ra/angstrom
              da = da/kcalmol
              CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,name=iname,xi=ra,di=da)
              EXIT
            END IF
          END DO
          CALL parser_release(parser)
        END DO
      END IF

    END IF

    ! extract simple atom interaction radii
    DO ikind = 1, nkind
      CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)
      radmax = (dftb_potential(ikind,ikind)%ngrdcut + 1) * &
               dftb_potential(ikind,ikind)%dgrd*0.5_dp
      CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,cutoff=radmax)
    END DO
    DO ikind = 1, nkind
      CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)
      CALL get_dftb_atom_param(dftb_parameter=dftb_atom_a,cutoff=ra)
      DO jkind = 1, nkind
        CALL get_qs_kind(qs_kind_set(jkind),dftb_parameter=dftb_atom_b)
        CALL get_dftb_atom_param(dftb_parameter=dftb_atom_b,cutoff=rb)
        radmax = (dftb_potential(ikind,jkind)%ngrdcut + 1) * &
                 dftb_potential(ikind,jkind)%dgrd
        IF ( ra+rb < radmax ) THEN
          ra = ra + (radmax-ra-rb)*0.5_dp
          rb = rb + (radmax-ra-rb)*0.5_dp
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,cutoff=ra)
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_b,cutoff=rb)
        END IF
      END DO
    END DO

    ! set correct core charge in potential
    DO ikind = 1, nkind
      CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom_a)
      CALL get_dftb_atom_param(dftb_parameter=dftb_atom_a,zeff=zeff)
      CALL set_potential(potential=qs_kind_set(ikind)%all_potential,&
                  zeff=zeff,zeff_correction=0.0_dp)
    END DO

    ! setup DFTB3 parameters
    IF ( dftb_control%dftb3_diagonal ) THEN
      DO ikind = 1, nkind
        CALL get_qs_kind(qs_kind_set(ikind), dftb3_param=db)
        CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_atom_a)
        CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,dudq=db)
      END DO
    END IF

    ! setup dispersion parameters (UFF type)
    IF ( dftb_control%dispersion ) THEN
      IF ( dftb_control%dispersion_type == dispersion_uff ) THEN
        eps_disp = dftb_control%eps_disp
        DO ikind = 1, nkind
          CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_atom_a)
          CALL get_dftb_atom_param(dftb_parameter=dftb_atom_a,xi=ra,di=da)
          rcdisp = 0._dp
          DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind),dftb_parameter=dftb_atom_b)
            CALL get_dftb_atom_param(dftb_parameter=dftb_atom_b,xi=rb,di=db)
            xij = SQRT(ra*rb)
            dij = SQRT(da*db)
            dftb_potential(ikind,jkind)%xij  = xij
            dftb_potential(ikind,jkind)%dij  = dij
            dftb_potential(ikind,jkind)%x0ij = xij*(0.5_dp**(1.0_dp/6.0_dp))
            dftb_potential(ikind,jkind)%a    = dij*396.0_dp/25.0_dp
            dftb_potential(ikind,jkind)%b    = &
                  dij/(xij**5)*672.0_dp*2.0_dp**(5.0_dp/6.0_dp)/25.0_dp
            dftb_potential(ikind,jkind)%c    = &
                  -dij/(xij**10)*2.0_dp**(2.0_dp/3.0_dp)*552.0_dp/25.0_dp
            rmax6 = ((8._dp*pi*dij/eps_disp)*xij**6)**0.25_dp
            rcdisp = MAX(rcdisp,rmax6*0.5_dp)
          END DO
          CALL set_dftb_atom_param(dftb_parameter=dftb_atom_a,rcdisp=rcdisp)
        END DO
      END IF
    END IF

    RETURN

1   CONTINUE
    CPASSERT(.FALSE.)

  END SUBROUTINE qs_dftb_param_init

! *****************************************************************************
!> \brief   Transform Slako format in l1/l2/m format
!> \param xmat ...
!> \param la ...
!> \param lb ...
!> \par Notes
!>         Slako tables from Dresden/Paderborn/Heidelberg groups are
!>         stored in the following native format:
!>
!>         Convention: Higher angular momenta are always on the right-hand side
!>
!>         1: d - d - delta
!>         2: d - d - pi
!>         3: d - d - sigma
!>         4: p - d - pi
!>         5: p - d - sigma
!>         6: p - p - pi
!>         7: p - p - sigma
!>         8: d - s - sigma
!>         9: p - s - sigma
!>        10: s - s - sigma
!> \version 1.0
! *****************************************************************************
  SUBROUTINE skreorder(xmat,la,lb)
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: xmat
    INTEGER, INTENT(IN)                      :: la, lb

    INTEGER                                  :: i, l1, l2, llm, m
    REAL(dp)                                 :: skllm(0:3,0:3,0:3)

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

END MODULE qs_dftb_parameters

