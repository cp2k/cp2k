!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/external_potential_types [1.0] *
!!
!!   NAME
!!     external_potential_types
!!
!!   FUNCTION
!!     Definition of the atomic potential types.
!!
!!   AUTHOR
!!     Matthias Krack (04.07.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE external_potential_types

! *****************************************************************************

! Index:

! SUBROUTINE allocate_all_potential(all_potential)
! SUBROUTINE allocate_gth_potential(gth_potential)
! SUBROUTINE deallocate_all_potential(all_potential)
! SUBROUTINE deallocate_gth_potential(gth_potential)
! SUBROUTINE init_cprj_ppnl(gth_potential)
! SUBROUTINE init_gth_potential(gth_potential)
! SUBROUTINE init_vprj_ppnl(gth_potential)
! SUBROUTINE read_all_potential(element_symbol,potential_name,alpha_erf,&
!                               elec_conf,all_potential,globenv)
! SUBROUTINE read_gth_potential(element_symbol,potential_name,alpha_erf,&
!                               elec_conf,gth_potential,globenv)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  IMPLICIT NONE

  PRIVATE

! *** Define the all-electron potential type ***

! *** Literature: M. Krack and M. Parrinello,
!                 Phys. Chem. Chem. Phys. 2, 2105 (2000)

  TYPE all_potential_type
    CHARACTER(LEN=60), DIMENSION(2) :: description
    REAL(wp)                        :: erfc_radius
  END TYPE all_potential_type

! *** Define the GTH potential type ***

! *** Literature: - S. Goedecker, M. Teter and J. Hutter,
!                   Phys. Rev. B 54, 1703 (1996)
!                 - C. Hartwigsen, S. Goedecker and J. Hutter,
!                   Phys. Rev. B 58, 3641 (1998)

  TYPE gth_potential_type
    CHARACTER(LEN=60), DIMENSION(3)     :: description
    REAL(wp)                            :: alpha_exp_ppl,erf_ppl_radius,&
                                           exp_ppl_radius,prj_ppnl_radius
    INTEGER                             :: lppnl,lprj_ppnl_max,nexp_ppl,&
                                           nppnl,nprj_ppnl_max
    REAL(wp), DIMENSION(:), POINTER     :: alpha_prj_ppnl,cexp_ppl
    INTEGER, DIMENSION(:), POINTER      :: elec_conf,nprj_ppnl
    REAL(wp), DIMENSION(:,:), POINTER   :: cprj,cprj_ppnl,vprj_ppnl
    REAL(wp), DIMENSION(:,:,:), POINTER :: hprj_ppnl
  END TYPE gth_potential_type

  INTEGER :: maxlppl = 0,&
             maxlppnl = 0

! *** Public subroutines ***

  PUBLIC :: allocate_potential,&
            deallocate_potential,&
            init_potential,&
            read_potential

! *** Public data types ***

  PUBLIC :: all_potential_type,&
            gth_potential_type

! *** Public variables ***

  PUBLIC :: maxlppl,&
            maxlppnl

! *****************************************************************************

  INTERFACE allocate_potential
    MODULE PROCEDURE allocate_all_potential,allocate_gth_potential
  END INTERFACE

  INTERFACE deallocate_potential
    MODULE PROCEDURE deallocate_all_potential,deallocate_gth_potential
  END INTERFACE

  INTERFACE init_potential
    MODULE PROCEDURE init_gth_potential
  END INTERFACE

  INTERFACE read_potential
    MODULE PROCEDURE read_all_potential,read_gth_potential
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_all_potential(all_potential)

!   Purpose: Allocate an atomic all-electron potential data structure.

!   History: Creation (25.07.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(all_potential_type), POINTER :: all_potential

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_all_potential (MODULE external_potential_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(all_potential)) CALL deallocate_all_potential(all_potential)

    ALLOCATE (all_potential,STAT=istat)

    IF (istat /= 0) CALL stop_memory(routine,"all_potential",0)

    all_potential%description(1) = "All-Electron potential"
    all_potential%description(2) =&
      "- Krack, Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)"

  END SUBROUTINE allocate_all_potential

! *****************************************************************************

  SUBROUTINE allocate_gth_potential(gth_potential)

!   Purpose: Allocate an atomic GTH potential data structure.

!   History: Creation (25.07.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(gth_potential_type), POINTER :: gth_potential

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_gth_potential (MODULE external_potential_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(gth_potential)) CALL deallocate_gth_potential(gth_potential)

    ALLOCATE (gth_potential,STAT=istat)

    IF (istat /= 0) CALL stop_memory(routine,"gth_potential",0)

    NULLIFY (gth_potential%alpha_prj_ppnl)
    NULLIFY (gth_potential%cexp_ppl)
    NULLIFY (gth_potential%elec_conf)
    NULLIFY (gth_potential%cprj)
    NULLIFY (gth_potential%cprj_ppnl)
    NULLIFY (gth_potential%hprj_ppnl)
    NULLIFY (gth_potential%nprj_ppnl)
    NULLIFY (gth_potential%vprj_ppnl)

    gth_potential%description(1) = "Goedecker-Teter-Hutter pseudopotential"
    gth_potential%description(2) =&
      "- Goedecker, Teter, Hutter, PRB 54, 1703 (1996)"
    gth_potential%description(3) =&
      "- Hartwigsen, Goedecker, Hutter, PRB 58, 3641 (1998)"

  END SUBROUTINE allocate_gth_potential

! *****************************************************************************

  SUBROUTINE deallocate_all_potential(all_potential)

!   Purpose: Deallocate an atomic all-electron potential data structure.

!   History: Creation (03.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(all_potential_type), POINTER :: all_potential

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_all_potential (MODULE external_potential_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (all_potential,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"all_potential")

  END SUBROUTINE deallocate_all_potential

! *****************************************************************************

  SUBROUTINE deallocate_gth_potential(gth_potential)

!   Purpose: Deallocate an atomic GTH potential data structure.

!   History: Creation (03.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(gth_potential_type), POINTER :: gth_potential

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_gth_potential (MODULE external_potential_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (gth_potential%elec_conf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gth_potential%elec_conf")

!   *** Deallocate the parameters of the local part ***

    IF (gth_potential%nexp_ppl > 0) THEN
      DEALLOCATE (gth_potential%cexp_ppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%cexp_ppl")
    END IF

!   *** Deallocate the parameters of the non-local part ***

    IF (gth_potential%nppnl > 0) THEN
      DEALLOCATE (gth_potential%alpha_prj_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%alpha_prj_ppnl")
      DEALLOCATE (gth_potential%cprj,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%cprj")
      DEALLOCATE (gth_potential%cprj_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%cprj_ppnl")
      DEALLOCATE (gth_potential%hprj_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%hprj_ppnl")
      DEALLOCATE (gth_potential%nprj_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%nprj_ppnl")
      DEALLOCATE (gth_potential%vprj_ppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"gth_potential%vprj_ppnl")
    END IF

    DEALLOCATE (gth_potential,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gth_potential")

  END SUBROUTINE deallocate_gth_potential

! *****************************************************************************

  SUBROUTINE init_cprj_ppnl(gth_potential)

!   Purpose: Initialize the coefficients of the projectors of the non-local
!            part of the GTH pseudopotential and the transformation matrices
!            for Cartesian overlap integrals between the orbital basis
!            functions and the projector functions.

!   History: - Creation (16.10.2000, Matthias Krack)

!   ***************************************************************************

    USE mathconstants,    ONLY: dfac,fac,pi
    USE orbital_pointers, ONLY: co,coset,nco

    TYPE(gth_potential_type), POINTER :: gth_potential

!   *** Local variables ***

    REAL(wp) :: alpha_prj_ppnl,cp
    INTEGER  :: cpx,cpy,cpz,cx,cy,cz,ico,iprj,iprj_ppnl,l,lp,lprj_ppnl,nprj,&
                px,py,pz

!   ---------------------------------------------------------------------------

    nprj = 0

    DO l=0,gth_potential%lppnl
      alpha_prj_ppnl = gth_potential%alpha_prj_ppnl(l)
      DO iprj_ppnl=1,gth_potential%nprj_ppnl(l)
        lp = iprj_ppnl - 1
        lprj_ppnl = l + 2*lp
        cp = SQRT(2.0_wp**(2.0_wp*REAL(lprj_ppnl,wp) + 3.5_wp)*&
                  alpha_prj_ppnl**(REAL(lprj_ppnl,wp) + 1.5_wp)/&
                  (SQRT(pi)*dfac(2*lprj_ppnl + 1)))
        gth_potential%cprj_ppnl(iprj_ppnl,l) = cp
        DO cx=0,l
          DO cy=0,l-cx
            cz = l - cx - cy
            iprj = nprj + co(cx,cy,cz)
            DO px=0,lp
              DO py=0,lp-px
                pz = lp - px - py
                cpx = cx + 2*px
                cpy = cy + 2*py
                cpz = cz + 2*pz
                ico = coset(cpx,cpy,cpz)
                gth_potential%cprj(ico,iprj) =&
                  cp*fac(lp)/(fac(px)*fac(py)*fac(pz))
              END DO
            END DO
          END DO
        END DO
        nprj = nprj + nco(l)
      END DO
    END DO

  END SUBROUTINE init_cprj_ppnl

! *****************************************************************************

  SUBROUTINE init_gth_potential(gth_potential)

!   Purpose: Initialize a GTH potential data structure.

!   History: Creation (27.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER :: gth_potential

!   ---------------------------------------------------------------------------

!   *** Initialize the projector coefficients of the    ***
!   *** non-local part of the GTH pseudopotential and   ***
!   *** the transformation matrices "pgf" -> "prj_ppnl" ***

    CALL init_cprj_ppnl(gth_potential)

!   *** Initialize the h(i,j) projector coefficients of ***
!   *** the non-local part of the GTH pseudopotential   ***

    CALL init_vprj_ppnl(gth_potential)

  END SUBROUTINE init_gth_potential

! *****************************************************************************

  SUBROUTINE init_vprj_ppnl(gth_potential)

!   Purpose: Initialize the h(i,j) projector coefficients of the non-local part
!            of the GTH pseudopotential.

!   History: - Creation (24.10.2000, Matthias Krack)

!   ***************************************************************************

    USE orbital_pointers,                ONLY: nco,nso
    USE orbital_transformation_matrices, ONLY: orbtramat

    TYPE(gth_potential_type), POINTER :: gth_potential

!   *** Local variables ***

    INTEGER :: i,ico,iprj,iprj_ppnl,iso,j,jco,jprj,jprj_ppnl,l,nprj

!   ---------------------------------------------------------------------------

    nprj = 0

    DO l=0,gth_potential%lppnl
      DO iprj_ppnl=1,gth_potential%nprj_ppnl(l)
        iprj = nprj + (iprj_ppnl - 1)*nco(l)
        DO jprj_ppnl=1,gth_potential%nprj_ppnl(l)
          jprj = nprj + (jprj_ppnl - 1)*nco(l)
          DO ico=1,nco(l)
            i = iprj + ico
            DO jco=1,nco(l)
              j = jprj + jco
              DO iso=1,nso(l)
                gth_potential%vprj_ppnl(i,j) =&
                  gth_potential%vprj_ppnl(i,j) +&
                  orbtramat(l)%slm(iso,ico)*&
                  gth_potential%hprj_ppnl(iprj_ppnl,jprj_ppnl,l)*&
                  orbtramat(l)%slm(iso,jco)
              END DO
            END DO
          END DO
        END DO
      END DO
      nprj = nprj + gth_potential%nprj_ppnl(l)*nco(l)
    END DO

  END SUBROUTINE init_vprj_ppnl

! *****************************************************************************

  SUBROUTINE read_all_potential(element_symbol,potential_name,alpha_erf,&
                                elec_conf,all_potential,globenv)

!   Purpose: Read an atomic all-electron potential data structure from the
!            database file.

!   History: - Creation (14.05.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,     ONLY: global_environment_type
    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session,&
                                test_object
    USE memory_utilities, ONLY: reallocate
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(all_potential_type), POINTER         :: all_potential
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: element_symbol,potential_name
    REAL(wp), INTENT(OUT)                     :: alpha_erf
    INTEGER, DIMENSION(:), POINTER            :: elec_conf

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_all_potential (MODULE external_potential_types)"

!   *** Local variables ***

    CHARACTER(LEN=LEN(potential_name)) :: apname
    CHARACTER(LEN=240)                 :: line
    CHARACTER(LEN=LEN(element_symbol)) :: symbol
    CHARACTER(LEN=11)                  :: test_result
    REAL(wp)                           :: r
    INTEGER                            :: l
    LOGICAL                            :: found

!   ---------------------------------------------------------------------------

    CALL start_input_session(globenv%potential_file_name,globenv)

!   *** Search for the requested potential in the potential file   ***
!   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol

    search_loop: DO

      CALL search(apname,.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(apname)

        IF ((INDEX(" "//line//" "," "//TRIM(apname)//" ") > 0).AND.&
            (INDEX(" "//line//" "," "//TRIM(symbol)//" ") > 0)) THEN

!         *** Read the electronic configuration ***

          NULLIFY (elec_conf)

          l = 0
          elec_conf => reallocate(elec_conf,0,l)
          CALL read_object(elec_conf(l),skip_lines=1)
          CALL test_object(test_result)

          DO WHILE (test_result == "integer")
            l = l + 1
            elec_conf => reallocate(elec_conf,0,l)
            CALL read_object(elec_conf(l))
            CALL test_object(test_result)
          END DO

!         *** Read r(loc) to define the exponent ***
!         *** of the nuclear charge distribution ***

          CALL read_object(r,skip_lines=1)

          alpha_erf = 1.0_wp/(SQRT(2.0_wp)*r)

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine,&
                          "The requested atomic potential <"//&
                          TRIM(potential_name)//&
                          "> was not found in the potential file <"//&
                          TRIM(globenv%potential_file_name)//">")

      END IF

    END DO search_loop

    CALL finish_input_session

  END SUBROUTINE read_all_potential

! *****************************************************************************

  SUBROUTINE read_gth_potential(element_symbol,potential_name,alpha_erf,&
                                elec_conf,gth_potential,globenv)

!   Purpose: Read an atomic GTH potential data structure from the database
!            file.

!   Literature: - S. Goedecker, M. Teter and J. Hutter,
!                 Phys. Rev. B 54, 1703 (1996)
!               - C. Hartwigsen, S. Goedecker and J. Hutter,
!                 Phys. Rev. B 58, 3641 (1998)

!   History: - Creation (14.05.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,     ONLY: global_environment_type
    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session,&
                                test_object
    USE mathlib,          ONLY: symmetrize_matrix
    USE memory_utilities, ONLY: reallocate
    USE orbital_pointers, ONLY: init_orbital_pointers,nco,ncoset
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_memory,stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(gth_potential_type), POINTER         :: gth_potential
    CHARACTER(LEN=*), INTENT(IN)              :: element_symbol,potential_name
    REAL(wp), INTENT(OUT)                     :: alpha_erf
    INTEGER, DIMENSION(:), POINTER            :: elec_conf

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_gth_potential (MODULE external_potential_types)"

!   *** Local variables ***

    CHARACTER(LEN=LEN(potential_name)) :: apname
    CHARACTER(LEN=240)                 :: line
    CHARACTER(LEN=LEN(element_symbol)) :: symbol
    CHARACTER(LEN=11)                  :: test_result
    REAL(wp)                           :: alpha,r
    INTEGER                            :: i,istat,j,l,lppnl,lprj_ppnl,&
                                          lprj_ppnl_max,n,nppnl,nprj_ppnl,&
                                          nprj_ppnl_max
    LOGICAL                            :: found

    INTEGER, DIMENSION(:), POINTER      :: all_elec_conf,gth_elec_conf
    REAL(wp), DIMENSION(:,:,:), POINTER :: hprj_ppnl

!   ---------------------------------------------------------------------------

    CALL start_input_session(globenv%potential_file_name,globenv)

!   *** Search for the requested potential in the potential file   ***
!   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol

    search_loop: DO

      CALL search(apname,.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(apname)

        IF ((INDEX(" "//line//" "," "//TRIM(apname)//" ") > 0).AND.&
            (INDEX(" "//line//" "," "//TRIM(symbol)//" ") > 0)) THEN

!         *** Read the electronic configuration (core and valence part) ***

          NULLIFY (all_elec_conf)

          l = 0
          all_elec_conf => reallocate(all_elec_conf,0,l)
          CALL read_object(all_elec_conf(l),skip_lines=1)
          CALL test_object(test_result)

          DO WHILE (test_result == "integer")
            l = l + 1
            all_elec_conf => reallocate(all_elec_conf,0,l)
            CALL read_object(all_elec_conf(l))
            CALL test_object(test_result)
          END DO

          elec_conf => reallocate(elec_conf,0,l)
          elec_conf(:) = all_elec_conf(:)

          DEALLOCATE (all_elec_conf)

!         *** Read the electronic configuration of    ***
!         *** the GTH-pseudo atom (valence part only) ***

          NULLIFY (gth_elec_conf)

          l = 0
          gth_elec_conf => reallocate(gth_elec_conf,0,l)
          CALL read_object(gth_elec_conf(l),skip_lines=1)
          CALL test_object(test_result)

          DO WHILE (test_result == "integer")
            l = l + 1
            gth_elec_conf => reallocate(gth_elec_conf,0,l)
            CALL read_object(gth_elec_conf(l))
            CALL test_object(test_result)
          END DO

          gth_potential%elec_conf => reallocate(gth_potential%elec_conf,0,l)
          gth_potential%elec_conf(:) = gth_elec_conf(:)

          DEALLOCATE (gth_elec_conf)

!         *** Read r(loc) to define the exponent ***
!         *** of the nuclear charge distribution ***

          CALL read_object(r,skip_lines=1)

          alpha_erf = 1.0_wp/(SQRT(2.0_wp)*r)

          gth_potential%alpha_exp_ppl = alpha_erf**2

!         *** Read the parameters for the local part ***
!         *** of the GTH pseudopotential (ppl)       ***

          CALL read_object(n)
          gth_potential%nexp_ppl = n
          gth_potential%cexp_ppl => reallocate(gth_potential%cexp_ppl,1,n)

          DO i=1,n
            CALL read_object(gth_potential%cexp_ppl(i))
            gth_potential%cexp_ppl(i) =&
              (2.0_wp*gth_potential%alpha_exp_ppl)**(i - 1)*&
              gth_potential%cexp_ppl(i)
          END DO

          maxlppl = MAX(maxlppl,2*(n - 1))

          CALL init_orbital_pointers(maxlppl)

!         *** Read the parameters for the non-local  ***
!         *** part of the GTH pseudopotential (ppnl) ***

          CALL read_object(n,skip_lines=1)

          gth_potential%lppnl = n - 1
          gth_potential%nppnl = 0

          gth_potential%lprj_ppnl_max = n - 1
          gth_potential%nprj_ppnl_max = 0

          IF (n > 0) THEN

            lppnl = gth_potential%lppnl
            nppnl = gth_potential%nppnl

            CALL init_orbital_pointers(lppnl)

            NULLIFY (hprj_ppnl)

!           *** Load the parameter for n non-local projectors ***

            gth_potential%alpha_prj_ppnl =>&
              reallocate(gth_potential%alpha_prj_ppnl,0,lppnl)
            gth_potential%nprj_ppnl =>&
              reallocate(gth_potential%nprj_ppnl,0,lppnl)

            lprj_ppnl_max = -1
            nprj_ppnl_max = 0

            DO l=0,lppnl
              CALL read_object(r,skip_lines=1)
              gth_potential%alpha_prj_ppnl(l) = 1.0_wp/(2.0_wp*r**2)
              CALL read_object(nprj_ppnl)
              gth_potential%nprj_ppnl(l) = nprj_ppnl
              nppnl = nppnl + nprj_ppnl*nco(l)
              IF (nprj_ppnl > nprj_ppnl_max) THEN
                nprj_ppnl_max = nprj_ppnl
                hprj_ppnl => reallocate(hprj_ppnl,1,nprj_ppnl_max,&
                                        1,nprj_ppnl_max,0,lppnl)
              END IF
              DO i=1,nprj_ppnl
                IF (i == 1) THEN
                  CALL read_object(hprj_ppnl(i,i,l))
                ELSE
                  CALL read_object(hprj_ppnl(i,i,l),skip_lines=1)
                END IF
                DO j=i+1,nprj_ppnl
                  CALL read_object(hprj_ppnl(i,j,l))
                END DO
              END DO
              IF (nprj_ppnl > 1) THEN
                CALL symmetrize_matrix(hprj_ppnl(:,:,l),"upper_to_lower")
              END IF
              lprj_ppnl_max = MAX(lprj_ppnl_max,l + 2*(nprj_ppnl - 1))
            END DO

            gth_potential%nppnl = nppnl

            maxlppnl = MAX(maxlppnl,lprj_ppnl_max)

            CALL init_orbital_pointers(maxlppnl)

            gth_potential%lprj_ppnl_max = lprj_ppnl_max
            gth_potential%nprj_ppnl_max = nprj_ppnl_max

            gth_potential%hprj_ppnl =>&
              reallocate(gth_potential%hprj_ppnl,1,nprj_ppnl_max,&
                         1,nprj_ppnl_max,0,lppnl)
            gth_potential%hprj_ppnl(:,:,:) = hprj_ppnl(:,:,:)

            gth_potential%cprj =>&
              reallocate(gth_potential%cprj,1,ncoset(lprj_ppnl_max),1,nppnl)
            gth_potential%cprj_ppnl =>&
              reallocate(gth_potential%cprj_ppnl,1,nprj_ppnl_max,0,lppnl)
            gth_potential%vprj_ppnl =>&
              reallocate(gth_potential%vprj_ppnl,1,nppnl,1,nppnl)

            DEALLOCATE (hprj_ppnl)

          END IF

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine,&
                          "The requested atomic potential <"//&
                          TRIM(potential_name)//&
                          "> was not found in the potential file <"//&
                          TRIM(globenv%potential_file_name)//">")

      END IF

    END DO search_loop

    CALL finish_input_session

  END SUBROUTINE read_gth_potential

! *****************************************************************************

END MODULE external_potential_types
