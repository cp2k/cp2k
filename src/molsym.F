!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Molecular symmetry routines
!> \par History
!>      2008 adapted from older routines by Matthias Krack
!> \author jgh
! *****************************************************************************
MODULE molsym
  
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: angle,&
                                             build_rotmat,&
                                             jacobi,&
                                             reflect_vector,&
                                             rotate_vector,&
                                             unit_matrix,&
                                             vector_product
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molsym'

  INTEGER, PARAMETER :: maxcn = 20,&
                        maxsec = maxcn + 1,&
                        maxses = 2*maxcn + 1,&
                        maxsig = maxcn + 1,&
                        maxsn  = 2*maxcn

  PUBLIC :: molsym_type
  PUBLIC :: release_molsym, molecular_symmetry, print_symmetry

!***

! *****************************************************************************
!> \brief Container for information about molecular symmetry
!> \param ain               : Atomic identification numbers (symmetry code).
!> \param aw                : Atomic weights for the symmetry analysis.
!> \param eps_geo           : Accuracy requested for analysis
!> \param inptostd          : Transformation matrix for input to standard orientation.
!> \param point_group_symbol: Point group symbol.
!> \param rotmat            : Rotation matrix.
!> \param sec               : List of C axes
!> \param                     (sec(:,i,j) => x,y,z of the ith j-fold C axis).
!> \param ses               : List of S axes
!> \param                     (ses(:,i,j) => x,y,z of the ith j-fold S axis).
!> \param sig               : List of mirror planes
!> \param                     (sig(:,i) => x,y,z of the ith mirror plane).
!> \param center_of_mass    : Shift vector for the center of mass.
!> \param tenmat            : Molecular tensor of inertia.
!> \param tenval            : Eigenvalues of the molecular tensor of inertia.
!> \param tenvec            : Eigenvectors of the molecular tensor of inertia.
!> \param group_of          : Group of equivalent atom.
!> \param llequatom         : Lower limit of a group in symequ_list.
!> \param ncn               : Degree of the C axis with highest degree.
!> \param ndod              : Number of substituted dodecahedral angles.
!> \param nequatom          : Number of equivalent atoms.
!> \param ngroup            : Number of groups of equivalent atoms.
!> \param nico              : Number of substituted icosahedral angles.
!> \param nlin              : Number of substituted angles of 180 degrees.
!> \param nsec              : Number of C axes.
!> \param nses              : Number of S axes.
!> \param nsig              : Number of mirror planes.
!> \param nsn               : Degree of the S axis with highest degree.
!> \param ntet              : Number of substituted tetrahedral angles.
!> \param point_group_order : Group order.
!> \param symequ_list       : List of all atoms ordered in groups of equivalent atoms.
!> \param ulequatom         : Upper limit of a group in symequ_list.
!> \param cubic : .TRUE., if a cubic point group was found (T,Th,Td,O,Oh,I,Ih).
!> \param dgroup: .TRUE., if a point group of D symmetry was found (Dn,Dnh,Dnd)
!> \param igroup: .TRUE., if a point group of icosahedral symmetry was found (I,Ih).
!> \param invers: .TRUE., if the molecule has a center of inversion.
!> \param linear: .TRUE., if the molecule is linear.
!> \param maxis : .TRUE., if the molecule has a main axis.
!> \param ogroup: .TRUE., if a point group of octahedral symmetry was found (O,Oh)
!> \param sgroup: .TRUE., if a point group of S symmetry was found.
!> \param sigmad: .TRUE., if there is a sigma_d mirror plane.
!> \param sigmah: .TRUE., if there is a sigma_h mirror plane.
!> \param sigmav: .TRUE., if there is a sigma_v mirror plane.
!> \param tgroup: .TRUE., if a point group of tetrahedral symmetry was found (T,Th,Td).
!> \par History
!>      05.2008 created
!> \author jgh
! *****************************************************************************
  TYPE molsym_type
    CHARACTER(4)                               :: point_group_symbol
    INTEGER                                    :: point_group_order
    INTEGER                                    :: ncn,ndod,ngroup,nico,nlin,nsig,nsn,ntet
    LOGICAL                                    :: cubic,dgroup,igroup,invers,linear,maxis,&
                                                  ogroup,sgroup,sigmad,sigmah,sigmav,tgroup
    REAL(KIND=dp)                              :: eps_geo
    REAL(KIND=dp), DIMENSION(3)                :: center_of_mass,tenval
    REAL(KIND=dp), DIMENSION(3)                :: x_axis,y_axis,z_axis
    REAL(KIND=dp), DIMENSION(:), POINTER       :: ain,aw
    REAL(KIND=dp), DIMENSION(3,3)              :: inptostd,rotmat,tenmat,tenvec
    REAL(KIND=dp), DIMENSION(3,maxsig)         :: sig
    REAL(KIND=dp), DIMENSION(3,maxsec,2:maxcn) :: sec
    REAL(KIND=dp), DIMENSION(3,maxses,2:maxsn) :: ses
    INTEGER, DIMENSION(maxcn)                  :: nsec
    INTEGER, DIMENSION(maxsn)                  :: nses
    INTEGER, DIMENSION(:), POINTER             :: group_of,llequatom,nequatom,&
                                                  symequ_list,ulequatom
  END TYPE molsym_type

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief  Create an object of molsym type
!> \param sym ...
!> \param natoms ...
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_molsym(sym,natoms)
    TYPE(molsym_type), POINTER               :: sym
    INTEGER, INTENT(IN)                      :: natoms

    CHARACTER(len=*), PARAMETER :: routineN = 'create_molsym', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( ASSOCIATED(sym) ) CALL release_molsym(sym)

    ALLOCATE(sym,STAT=ierr)
    CPASSERT(ierr==0)

    ALLOCATE (sym%ain(natoms),sym%aw(natoms),sym%group_of(natoms),sym%llequatom(natoms),&
              sym%nequatom(natoms),sym%symequ_list(natoms),sym%ulequatom(natoms),STAT=ierr)
    CPASSERT(ierr==0)

  END SUBROUTINE create_molsym

! *****************************************************************************
!> \brief  release an object of molsym type
!> \param sym ...
!> \author jgh
! *****************************************************************************
  SUBROUTINE release_molsym(sym)
    TYPE(molsym_type), POINTER               :: sym

    CHARACTER(len=*), PARAMETER :: routineN = 'release_molsym', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(sym))

    IF ( ASSOCIATED(sym%aw) ) THEN
      DEALLOCATE(sym%aw,STAT=ierr)
      CPASSERT(ierr==0)
    END IF
    IF ( ASSOCIATED(sym%ain) ) THEN
      DEALLOCATE(sym%ain,STAT=ierr)
      CPASSERT(ierr==0)
    END IF
    IF ( ASSOCIATED(sym%group_of) ) THEN
      DEALLOCATE(sym%group_of,STAT=ierr)
      CPASSERT(ierr==0)
    END IF
    IF ( ASSOCIATED(sym%llequatom) ) THEN
      DEALLOCATE(sym%llequatom,STAT=ierr)
      CPASSERT(ierr==0)
    END IF
    IF ( ASSOCIATED(sym%nequatom) ) THEN
      DEALLOCATE(sym%nequatom,STAT=ierr)
      CPASSERT(ierr==0)
    END IF
    IF ( ASSOCIATED(sym%symequ_list) ) THEN
      DEALLOCATE(sym%symequ_list,STAT=ierr)
      CPASSERT(ierr==0)
    END IF
    IF ( ASSOCIATED(sym%ulequatom) ) THEN
      DEALLOCATE(sym%ulequatom,STAT=ierr)
      CPASSERT(ierr==0)
    END IF

    DEALLOCATE(sym,STAT=ierr)
    CPASSERT(ierr==0)

  END SUBROUTINE release_molsym

! *****************************************************************************

! *****************************************************************************
!> \brief Add a new C_n axis to the list sec, but first check, if the
!>         Cn axis is already in the list.
!> \param n ...
!> \param a ...
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE addsec(n,a,sym)
    INTEGER, INTENT(IN)                      :: n
    REAL(dp), DIMENSION(3), INTENT(IN)       :: a
    TYPE(molsym_type), INTENT(inout)         :: sym

    CHARACTER(len=*), PARAMETER :: routineN = 'addsec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isec
    LOGICAL                                  :: failure
    REAL(dp)                                 :: length_of_a, scapro
    REAL(dp), DIMENSION(3)                   :: d

    failure = .FALSE.
    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
    d(:) = a(:)/length_of_a

    ! Check, if the current Cn axis is already in the list
    DO isec=1,sym%nsec(n)
      scapro = sym%sec(1,isec,n)*d(1) + sym%sec(2,isec,n)*d(2) + sym%sec(3,isec,n)*d(3)
      IF (ABS(ABS(scapro) - 1.0_dp) < sym%eps_geo) RETURN
    END DO
    sym%ncn = MAX(sym%ncn,n)

    ! Add the new Cn axis to the list sec
    CPASSERT(sym%nsec(n) < maxsec)
    sym%nsec(1) = sym%nsec(1) + 1
    sym%nsec(n) = sym%nsec(n) + 1
    sym%sec(:,sym%nsec(n),n) = d(:)

  END SUBROUTINE addsec

! *****************************************************************************
!> \brief  Add a new Sn axis to the list ses, but first check, if the
!>         Sn axis is already in the list.
!> \param n ...
!> \param a ...
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE addses(n,a,sym)
    INTEGER, INTENT(IN)                      :: n
    REAL(dp), DIMENSION(3), INTENT(IN)       :: a
    TYPE(molsym_type), INTENT(inout)         :: sym

    CHARACTER(len=*), PARAMETER :: routineN = 'addses', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ises
    LOGICAL                                  :: failure
    REAL(dp)                                 :: length_of_a, scapro
    REAL(dp), DIMENSION(3)                   :: d

    failure = .FALSE.
    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
    d(:) = a(:)/length_of_a

    ! Check, if the current Sn axis is already in the list
    DO ises=1,sym%nses(n)
      scapro = sym%ses(1,ises,n)*d(1) + sym%ses(2,ises,n)*d(2) + sym%ses(3,ises,n)*d(3)
      IF (ABS(ABS(scapro) - 1.0_dp) < sym%eps_geo) RETURN
    END DO
    sym%nsn = MAX(sym%nsn,n)

    ! Add the new Sn axis to the list ses
    CPASSERT(sym%nses(n) < maxses)
    sym%nses(1) = sym%nses(1) + 1
    sym%nses(n) = sym%nses(n) + 1
    sym%ses(:,sym%nses(n),n) = d(:)

  END SUBROUTINE addses

! *****************************************************************************
!> \brief  Add a new mirror plane to the list sig, but first check, if the
!>         normal vector of the mirror plane is already in the list.
!> \param a ...
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE addsig(a,sym)
    REAL(dp), DIMENSION(3), INTENT(IN)       :: a
    TYPE(molsym_type), INTENT(inout)         :: sym

    CHARACTER(len=*), PARAMETER :: routineN = 'addsig', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isig
    LOGICAL                                  :: failure
    REAL(dp)                                 :: length_of_a, scapro
    REAL(dp), DIMENSION(3)                   :: d

    failure = .FALSE.
    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
    d(:) = a(:)/length_of_a

    ! Check, if the normal vector of the current mirror plane is already in the list
    DO isig=1,sym%nsig
      scapro = sym%sig(1,isig)*d(1) + sym%sig(2,isig)*d(2) + sym%sig(3,isig)*d(3)
      IF (ABS(ABS(scapro) - 1.0_dp) < sym%eps_geo) RETURN
    END DO

    ! Add the normal vector of the new mirror plane to the list sig
    CPASSERT(sym%nsig < maxsig)
    sym%nsig = sym%nsig + 1
    sym%sig(:,sym%nsig) = d(:)

  END SUBROUTINE addsig

! *****************************************************************************
!> \brief  Symmetry handling for only one atom.
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE atomsym(sym)
    TYPE(molsym_type), INTENT(inout)         :: sym

    CHARACTER(len=*), PARAMETER :: routineN = 'atomsym', &
      routineP = moduleN//':'//routineN

! Set point group symbol

    sym%point_group_symbol = "R(3)"

    ! Set variables like D*h
    sym%linear = .TRUE.
    sym%invers = .TRUE.
    sym%dgroup = .TRUE.
    sym%sigmah = .TRUE.

  END SUBROUTINE atomsym

! *****************************************************************************
!> \brief Search for point groups with AXial SYMmetry (Cn,Cnh,Cnv,Dn,Dnh,Dnd,S2n).
!> \param coord ...
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE axsym(coord,sym)
    REAL(Kind=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    TYPE(molsym_type), INTENT(inout)         :: sym

    CHARACTER(len=*), PARAMETER :: routineN = 'axsym', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, isig, jatom, m, n, &
                                                natoms, nx, ny, nz
    REAL(dp)                                 :: phi
    REAL(dp), DIMENSION(3)                   :: a, b, d

! Find the correct main axis and rotate main axis to z axis

    IF ((sym%ncn == 2).AND.(sym%nses(2) == 3)) THEN
      IF (sym%nsn == 4) THEN
        ! Special case: D2d
        phi = angle(sym%ses(:,1,sym%nsn),sym%z_axis(:))
        d(:) = vector_product(sym%ses(:,1,sym%nsn),sym%z_axis(:))
        CALL rotate_molecule(phi,d(:),sym,coord)
      ELSE
        ! Special cases: D2 and D2h
        phi = 0.5_dp*pi
        nx = naxis(sym%x_axis(:),coord,sym)
        ny = naxis(sym%y_axis(:),coord,sym)
        nz = naxis(sym%z_axis(:),coord,sym)
        IF ((nx > ny).AND.(nx > nz)) THEN
          CALL rotate_molecule(-phi,sym%y_axis(:),sym,coord)
        ELSE IF ((ny > nz).AND.(ny > nx)) THEN
          CALL rotate_molecule(phi,sym%x_axis(:),sym,coord)
        END IF
      END IF
    ELSE
      phi = angle(sym%sec(:,1,sym%ncn),sym%z_axis(:))
      d(:) = vector_product(sym%sec(:,1,sym%ncn),sym%z_axis(:))
      CALL rotate_molecule(phi,d(:),sym,coord)
    END IF

    ! Search for C2 axes perpendicular to the main axis
    ! Loop over all pairs of atoms (except on the z axis)
    natoms = SIZE(coord,2)
    DO iatom=1,natoms
      a(:) = coord(:,iatom)
      IF ((ABS(a(1)) > sym%eps_geo).OR.(ABS(a(2)) > sym%eps_geo)) THEN
        a(3) = 0.0_dp
        IF (caxis(2,a(:),sym,coord)) CALL addsec(2,a(:),sym)
        d(:) = vector_product(a(:),sym%z_axis(:))
        IF (sigma(d(:),sym,coord)) CALL addsig(d(:),sym)
        DO jatom=iatom+1,natoms
          b(:) = coord(:,jatom)
          IF ((ABS(b(1)) > sym%eps_geo).OR.(ABS(b(2)) > sym%eps_geo)) THEN
            b(3) = 0.0_dp
            phi = 0.5_dp*angle(a(:),b(:))
            d(:) = vector_product(a(:),b(:))
            b(:) = rotate_vector(a(:),phi,d(:))
            IF (caxis(2,b(:),sym,coord)) CALL addsec(2,b(:),sym)
            d(:) = vector_product(b(:),sym%z_axis)
            IF (sigma(d(:),sym,coord)) CALL addsig(d(:),sym)
          END IF
        END DO
      END IF
    END DO

    ! Check the xy plane for mirror plane
    IF (sigma(sym%z_axis(:),sym,coord)) THEN
      CALL addsig(sym%z_axis(:),sym)
      sym%sigmah = .TRUE.
    END IF

    ! Set logical variables for point group determination ***
    IF (sym%nsec(2) > 1) THEN
      sym%dgroup = .TRUE.
      IF ((.NOT.sym%sigmah).AND.(sym%nsig > 0)) sym%sigmad = .TRUE.
    ELSE
      IF ((.NOT.sym%sigmah).AND.(sym%nsig > 0)) sym%sigmav = .TRUE.
      ! Cnh groups with n>3 were wrong
      ! IF (sym%nsn > 3) sym%sgroup = .TRUE.
      IF ((.NOT.sym%sigmah).AND.(sym%nsn > 3)) sym%sgroup = .TRUE.
    END IF

    ! Rotate to standard orientation
    n = 0
    m = 0
    IF ((sym%dgroup.AND.sym%sigmah).OR.(sym%dgroup.AND.sym%sigmad).OR.sym%sigmav) THEN
    ! Dnh, Dnd or Cnv
      DO isig=1,sym%nsig
        IF (ABS(sym%sig(3,isig)) < sym%eps_geo) THEN
          n = nsigma(sym%sig(:,isig),sym,coord)
          IF (n > m) THEN
            m = n
            a(:) = sym%sig(:,isig)
          END IF
        END IF
      END DO
      IF (m > 0) THEN
        ! Check for special case: C2v with all atoms in a plane
        IF (sym%sigmav.AND.(m == natoms)) THEN
          phi = angle(a(:),sym%x_axis(:))
          d(:) = vector_product(a(:),sym%x_axis(:))
        ELSE
          phi = angle(a(:),sym%y_axis(:))
          d(:) = vector_product(a(:),sym%y_axis(:))
        END IF
        CALL rotate_molecule(phi,d(:),sym,coord)
      END IF
    ELSE IF (sym%sigmah) THEN
      ! Cnh
      DO iatom=1,natoms
        IF (ABS(coord(3,iatom)) < sym%eps_geo) THEN
          n = naxis(coord(:,iatom),coord,sym)
          IF (n > m) THEN
            m = n
            a(:) = coord(:,iatom)
          END IF
        END IF
      END DO
      IF (m > 0) THEN
        phi = angle(a(:),sym%x_axis(:))
        d(:) = vector_product(a(:),sym%x_axis(:))
        CALL rotate_molecule(phi,d(:),sym,coord)
      END IF
    ! No action for Dn, Cn or S2n ***
    END IF

  END SUBROUTINE axsym

! *****************************************************************************
!> \brief Generate a symmetry list to identify equivalent atoms.
!> \param sym ...
!> \param coord ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE build_symequ_list(sym,coord)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(Kind=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord

    CHARACTER(len=*), PARAMETER :: routineN = 'build_symequ_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, icn, iequatom, &
                                                incr, isec, ises, isig, isn, &
                                                jatom, jcn, jsn, natoms
    REAL(KIND=dp)                            :: phi
    REAL(KIND=dp), DIMENSION(3)              :: a

    natoms = SIZE(coord,2)

    IF (sym%sigmah) THEN
      incr = 1
    ELSE
      incr = 2
    END IF

    ! Initialize the atom and the group counter
    iatom = 1
    sym%ngroup = 1

    loop: DO

      ! Loop over all mirror planes
      DO isig=1,sym%nsig
        a(:) = reflect_vector(coord(:,iatom),sym%sig(:,isig))
        iequatom = equatom(iatom,a(:),sym,coord)
        IF ((iequatom > 0).AND.(.NOT.in_symequ_list(iequatom,sym))) THEN
          sym%ulequatom(sym%ngroup) = sym%ulequatom(sym%ngroup) + 1
          sym%symequ_list(sym%ulequatom(sym%ngroup)) = iequatom
          sym%nequatom(sym%ngroup) = sym%nequatom(sym%ngroup) + 1
        END IF
      END DO

      ! Loop over all Cn axes
      DO icn=2,sym%ncn
        DO isec=1,sym%nsec(icn)
          DO jcn=1,icn-1
            IF (newse(jcn,icn)) THEN
              phi = 2.0_dp*pi*REAL(jcn,KIND=dp)/REAL(icn,KIND=dp)
              a(:) = rotate_vector(coord(:,iatom),phi,sym%sec(:,isec,icn))
              iequatom = equatom(iatom,a(:),sym,coord)
              IF ((iequatom > 0).AND.(.NOT.in_symequ_list(iequatom,sym))) THEN
                sym%ulequatom(sym%ngroup) = sym%ulequatom(sym%ngroup) + 1
                sym%symequ_list(sym%ulequatom(sym%ngroup)) = iequatom
                sym%nequatom(sym%ngroup) = sym%nequatom(sym%ngroup) + 1
              END IF
            END IF
          END DO
        END DO
      END DO

      ! Loop over all Sn axes
      DO isn=2,sym%nsn
        DO ises=1,sym%nses(isn)
          DO jsn=1,isn-1,incr
            IF (newse(jsn,isn)) THEN
              phi = 2.0_dp*pi*REAL(jsn,KIND=dp)/REAL(isn,KIND=dp)
              a(:) = rotate_vector(coord(:,iatom),phi,sym%ses(:,ises,isn))
              a(:) = reflect_vector(a(:),sym%ses(:,ises,isn))
              iequatom = equatom(iatom,a(:),sym,coord)
              IF ((iequatom > 0).AND.(.NOT.in_symequ_list(iequatom,sym))) THEN
                sym%ulequatom(sym%ngroup) = sym%ulequatom(sym%ngroup) + 1
                sym%symequ_list(sym%ulequatom(sym%ngroup)) = iequatom
                sym%nequatom(sym%ngroup) = sym%nequatom(sym%ngroup) + 1
              END IF
            END IF
          END DO
        END DO
      END DO

      ! Exit loop, if all atoms are in the list
      IF (sym%ulequatom(sym%ngroup) == natoms) EXIT

      ! Search for the next atom which is not in the list yet
      DO jatom=2,natoms
        IF (.NOT.in_symequ_list(jatom,sym)) THEN
          iatom = jatom
          sym%ngroup = sym%ngroup + 1
          sym%llequatom(sym%ngroup) = sym%ulequatom(sym%ngroup-1) + 1
          sym%ulequatom(sym%ngroup) = sym%llequatom(sym%ngroup)
          sym%symequ_list(sym%llequatom(sym%ngroup)) = iatom
          CYCLE loop
        END IF
      END DO

    END DO loop

    ! Generate list vector group_of
    DO i=1,sym%ngroup
      DO iequatom=sym%llequatom(i),sym%ulequatom(i)
        sym%group_of(sym%symequ_list(iequatom)) = i
      END DO
    END DO

  END SUBROUTINE build_symequ_list

! *****************************************************************************
!> \brief Rotate the molecule about a n-fold axis defined by the vector a
!>        using a rotation angle of 2*pi/n. Check, if the axis is a Cn axis.
!> \param n ...
!> \param a ...
!> \param sym ...
!> \param coord ...
!> \retval caxis ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION caxis(n,a,sym,coord)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(Kind=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    LOGICAL                                  :: caxis

    INTEGER                                  :: iatom, natoms
    REAL(KIND=dp)                            :: length_of_a, phi
    REAL(KIND=dp), DIMENSION(3)              :: b

    caxis = .FALSE.
    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

    ! Check the length of the axis vector a
    natoms = SIZE(coord,2)
    IF (length_of_a > sym%eps_geo) THEN
      ! Calculate the rotation angle phi and build up the rotation matrix rotmat
      phi = 2.0_dp*pi/REAL(n,dp)
      CALL build_rotmat(phi,a(:),sym%rotmat(:,:))
      ! Check all atoms
      DO iatom=1,natoms
        ! Rotate the actual atom by 2*pi/n about a
        b(:) = MATMUL(sym%rotmat(:,:),coord(:,iatom))
        ! Check if the coordinates of the rotated atom
        ! are in the coordinate set of the molecule
        IF (equatom(iatom,b(:),sym,coord) == 0) RETURN
      END DO
      ! The axis defined by vector a is a Cn axis, thus return with caxis = .TRUE.
      caxis = .TRUE.
    END IF

  END FUNCTION caxis

! *****************************************************************************
!> \brief Check for a point group of cubic symmetry (I,Ih,O,Oh,T,Th,Td).
!> \param sym ...
!> \param coord ...
!> \param failed ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE cubsym(sym,coord,failed)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    LOGICAL, INTENT(INOUT)                   :: failed

    INTEGER                                  :: i, iatom, iax, ic3, isec, &
                                                jatom, jc3, jsec, katom, &
                                                natoms, nc3
    REAL(KIND=dp)                            :: phi1, phi2, phidd, phidi, &
                                                phiii
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d

! Angle between two adjacent atoms of the dodecahedron => <(C3,C3)

    phidd = ATAN(0.4_dp*SQRT(5.0_dp))
    ! Angle between two adjacent atoms of the icosahedron and the dodecahedron => <(C5,C3)
    phidi = ATAN(3.0_dp - SQRT(5.0_dp))
    ! Angle between two adjacent atoms of the icosahedron <(C5,C5)
    phiii = ATAN(2.0_dp)

    ! Find a pair of C3 axes
    natoms = SIZE(coord,2)
    DO iatom=1,natoms
      ! Check all atomic vectors for C3 axis
      IF (caxis(3,coord(:,iatom),sym,coord)) THEN
        CALL addsec(3,coord(:,iatom),sym)
        IF (sym%nsec(3) > 1) EXIT
      END IF
    END DO

    ! Check the center of mass vector of a triangle
    ! defined by three atomic vectors for C3 axis
    IF (sym%nsec(3) < 2) THEN

      loop: DO iatom=1,natoms
        a(:) = coord(:,iatom)
        DO jatom=iatom+1,natoms
          DO katom=jatom+1,natoms
            IF ((ABS(dist(coord(:,iatom),coord(:,jatom)) &
              - dist(coord(:,iatom),coord(:,katom))) < sym%eps_geo).AND.&
                (ABS(dist(coord(:,iatom),coord(:,jatom)) &
              - dist(coord(:,jatom),coord(:,katom))) < sym%eps_geo)) THEN
              b(:) = a(:) + coord(:,jatom) + coord(:,katom)
              IF (caxis(3,b(:),sym,coord)) THEN
                CALL addsec(3,b(:),sym)
                IF (sym%nsec(3) > 1) EXIT loop
              END IF
            END IF
          END DO
        END DO
      END DO loop

    END IF

    ! Return after unsuccessful search for a pair of C3 axes
    IF (sym%nsec(3) < 2) THEN
      failed = .TRUE.
      RETURN
    END IF

    ! Generate the remaining C3 axes
    nc3 = 0
    DO
      nc3 = sym%nsec(3)
      DO ic3=1,nc3
        a(:) = sym%sec(:,ic3,3)
        DO jc3=1,nc3
          IF (ic3 /= jc3) THEN
            d(:) = sym%sec(:,jc3,3)
            DO i=1,2
              phi1 = 2.0_dp*REAL(i,KIND=dp)*pi/3.0_dp
              b(:) = rotate_vector(a(:),phi1,d(:))
              CALL addsec(3,b(:),sym)
            END DO
          END IF
        END DO
      END DO

      IF (sym%nsec(3) > nc3) THEN
        ! New C3 axes were found
        CYCLE
      ELSE
        ! In the case of I or Ih there have to be more C3 axes
        IF (sym%nsec(3) == 4) THEN
          a(:) = sym%sec(:,1,3)
          b(:) = sym%sec(:,2,3)
          phi1 = 0.5_dp*angle(a(:),b(:))
          IF (phi1 < 0.5_dp*pi) phi1 = phi1 - 0.5_dp*pi
          d(:) = vector_product(a(:),b(:))
          b(:) = rotate_vector(a(:),phi1,d(:))
          c(:) = sym%sec(:,3,3)
          phi1 = 0.5_dp*angle(a(:),c(:))
          IF (phi1 < 0.5_dp*pi) phi1 = phi1 - 0.5_dp*pi
          d(:) = vector_product(a(:),c(:))
          c(:) = rotate_vector(a(:),phi1,d(:))
          d(:) = vector_product(b(:),c(:))
          phi1 = 0.5_dp*phidd
          a(:) = rotate_vector(b(:),phi1,d(:))
          IF (caxis(3,a(:),sym,coord)) THEN
            CALL addsec(3,a(:),sym)
          ELSE
            phi2 = 0.5_dp*pi - phi1
            a(:) = rotate_vector(b(:),phi2,d(:))
            IF (caxis(3,a(:),sym,coord)) CALL addsec(3,a(:),sym)
          END IF

          IF (sym%nsec(3) > 4) CYCLE

        ELSE IF (sym%nsec(3) /= 10) THEN

          !  C3 axes found, but only 4 or 10 are possible
          failed = .TRUE.
          RETURN

        END IF

        ! Exit loop, if 4 or 10 C3 axes were found
        EXIT

      END IF

    END DO

    ! Loop over all pairs of C3 axes to find all C5, C4 and C2 axes
    DO isec=1,sym%nsec(3)

      a(:) = sym%sec(:,isec,3)
      DO jsec=isec+1,sym%nsec(3)
        phi1 = 0.5_dp*angle(a(:),sym%sec(:,jsec,3))
        d(:) = vector_product(a(:),sym%sec(:,jsec,3))

        ! Check for C5 axis (I,Ih)
        IF (sym%nsec(3) == 10) THEN
          b(:) = rotate_vector(a(:),phidi,d(:))
          IF (caxis(5,b(:),sym,coord)) THEN
            CALL addsec(5,b(:),sym)
            phi1 = phidi + phiii
            b(:) = rotate_vector(a(:),phi1,d(:))
            IF (caxis(5,b(:),sym,coord)) CALL addsec(5,b(:),sym)
          END IF
        END IF

        ! Check for C4 (O,Oh), C2 and S2 (center of inversion) axis
        DO i=0,1
          phi2 = phi1 - 0.5_dp*REAL(i,KIND=dp)*pi
          b(:) = rotate_vector(a(:),phi2,d(:))
          IF (caxis(2,b(:),sym,coord)) THEN
            CALL addsec(2,b(:),sym)
            IF (sym%nsec(3) == 4) THEN
              IF (caxis(4,b(:),sym,coord)) CALL addsec(4,b(:),sym)
            END IF
            IF (saxis(2,b(:),sym,coord)) THEN
              CALL addses(2,b(:),sym)
              sym%invers = .TRUE.
            END IF
          END IF
          IF (sigma(b(:),sym,coord)) CALL addsig(b(:),sym)
        END DO

        ! Check for mirror plane
        IF (sigma(d(:),sym,coord)) CALL addsig(d(:),sym)

      END DO

    END DO

    ! Set the logical variables for point group determination
    iax = sym%ncn

    IF ((sym%ncn == 5).AND.(sym%nsec(5) > 1)) THEN
      sym%igroup = .TRUE.
    ELSE IF ((sym%ncn == 4).AND.(sym%nsec(4) > 1)) THEN
      sym%ogroup = .TRUE.
    ELSE IF ((sym%ncn == 3).AND.(sym%nsec(3) > 1)) THEN
      sym%tgroup = .TRUE.
      IF ((.NOT.sym%invers).AND.(sym%nsig > 0)) sym%sigmad = .TRUE.
      iax = 2
    ELSE
      ! No main axis found
      failed = .TRUE.
      RETURN
    END IF

    ! Rotate molecule to standard orientation
    phi1 = angle(sym%sec(:,1,iax),sym%z_axis(:))
    d(:) = vector_product(sym%sec(:,1,iax),sym%z_axis(:))
    CALL rotate_molecule(phi1,d(:),sym,coord)

    a(:) = sym%sec(:,2,iax)
    a(3) = 0.0_dp
    phi2 = angle(a(:),sym%x_axis(:))
    d(:) = vector_product(a(:),sym%x_axis(:))
    CALL rotate_molecule(phi2,d(:),sym,coord)

  END SUBROUTINE cubsym

! *****************************************************************************
!> \brief The number of a equivalent atom to atoma is returned. If there
!>        is no equivalent atom, then zero is returned. The cartesian
!>        coordinates of the equivalent atom are defined by the vector a.
!> \param atoma ...
!> \param a ...
!> \param sym ...
!> \param coord ...
!> \retval equatom ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION equatom(atoma,a,sym,coord)
    INTEGER, INTENT(IN)                      :: atoma
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    INTEGER                                  :: equatom

    INTEGER                                  :: iatom, natoms

    equatom = 0
    natoms = SIZE(coord,2)
    DO iatom=1,natoms
      IF ((ABS(sym%ain(iatom) - sym%ain(atoma)) < TINY(0.0_dp)).AND.&
          (ABS(a(1) - coord(1,iatom)) < sym%eps_geo).AND.&
          (ABS(a(2) - coord(2,iatom)) < sym%eps_geo).AND.&
          (ABS(a(3) - coord(3,iatom)) < sym%eps_geo)) THEN
        equatom = iatom
        RETURN
      END IF
    END DO

  END FUNCTION equatom

! *****************************************************************************
!> \brief Calculate the order of the point group.
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE get_point_group_order(sym)
    TYPE(molsym_type), INTENT(inout)         :: sym

    INTEGER                                  :: icn, incr, isec, ises, isn, &
                                                jcn, jsn

! Count all symmetry elements of the symmetry group
! First E and all mirror planes

    sym%point_group_order = 1 + sym%nsig

    ! Loop over all C axes
    DO icn=2,sym%ncn
      DO isec=1,sym%nsec(icn)
        DO jcn=1,icn-1
          IF (newse(jcn,icn)) sym%point_group_order = sym%point_group_order + 1
        END DO
      END DO
    END DO

    ! Loop over all S axes
    IF (sym%sigmah) THEN
      incr = 1
    ELSE
      incr = 2
    END IF

    DO isn=2,sym%nsn
      DO ises=1,sym%nses(isn)
        DO jsn=1,isn-1,incr
          IF (newse(jsn,isn)) sym%point_group_order = sym%point_group_order + 1
        END DO
      END DO
    END DO

  END SUBROUTINE get_point_group_order

! *****************************************************************************
!> \brief Get the point group symbol.
!> \param sym ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE get_point_group_symbol(sym)
    TYPE(molsym_type), INTENT(inout)         :: sym

    CHARACTER(4)                             :: degree

    IF (sym%linear) THEN
      IF (sym%invers) THEN
        sym%point_group_symbol = "D*h"
      ELSE
        sym%point_group_symbol = "C*v"
      END IF
    ELSE IF (sym%cubic) THEN
      IF (sym%igroup) THEN
        sym%point_group_symbol = "I"
      ELSE IF (sym%ogroup) THEN
        sym%point_group_symbol = "O"
      ELSE IF (sym%tgroup) THEN
        sym%point_group_symbol = "T"
      END IF
      IF (sym%invers) THEN
        sym%point_group_symbol = TRIM(sym%point_group_symbol)//"h"
      ELSE IF (sym%sigmad) THEN
        sym%point_group_symbol = TRIM(sym%point_group_symbol)//"d"
      END IF
    ELSE IF (sym%maxis) THEN
      IF (sym%dgroup) THEN
        WRITE (degree,"(I3)") sym%ncn
        sym%point_group_symbol = "D"//TRIM(ADJUSTL(degree))
        IF (sym%sigmah) THEN
          sym%point_group_symbol = TRIM(sym%point_group_symbol)//"h"
        ELSE IF (sym%sigmad) THEN
          sym%point_group_symbol = TRIM(sym%point_group_symbol)//"d"
        END IF
      ELSE IF (sym%sgroup) THEN
        WRITE (degree,"(I3)") sym%nsn
        sym%point_group_symbol = "S"//TRIM(ADJUSTL(degree))
      ELSE
        WRITE (degree,"(I3)") sym%ncn
        sym%point_group_symbol = "C"//TRIM(ADJUSTL(degree))
        IF (sym%sigmah) THEN
          sym%point_group_symbol = TRIM(sym%point_group_symbol)//"h"
        ELSE IF (sym%sigmav) THEN
          sym%point_group_symbol = TRIM(sym%point_group_symbol)//"v"
        END IF
      END IF
    ELSE
      IF (sym%invers) THEN
        sym%point_group_symbol = "Ci"
      ELSE IF (sym%sigmah) THEN
        sym%point_group_symbol = "Cs"
      ELSE
        sym%point_group_symbol = "C1"
      END IF
    ENDIF

  END SUBROUTINE get_point_group_symbol

! *****************************************************************************
!> \brief Initialization of the global variables of module symmetry.
!> \param sym ...
!> \param atype ...
!> \param weight ...
!> \par History
!>        Creation (19.10.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE init_symmetry(sym,atype,weight)
    TYPE(molsym_type), INTENT(inout)         :: sym
    INTEGER, DIMENSION(:), INTENT(IN)        :: atype
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: weight

    CHARACTER(len=*), PARAMETER :: routineN = 'init_symmetry', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, natoms

! Define the Cartesian axis vectors

    sym%x_axis = (/1.0_dp,0.0_dp,0.0_dp/)
    sym%y_axis = (/0.0_dp,1.0_dp,0.0_dp/)
    sym%z_axis = (/0.0_dp,0.0_dp,1.0_dp/)

    ! Initialize lists for symmetry elements
    sym%sec(:,:,:) = 0.0_dp
    sym%ses(:,:,:) = 0.0_dp
    sym%sig(:,:) = 0.0_dp

    sym%center_of_mass(:) = 0.0_dp

    ! Initialize symmetry element counters
    sym%ncn = 1
    sym%nsn = 1
    sym%nsec(:) = 0
    sym%nses(:) = 0
    sym%nsig = 0

    ! Initialize logical variables
    sym%cubic = .FALSE.
    sym%dgroup = .FALSE.
    sym%igroup = .FALSE.
    sym%invers = .FALSE.
    sym%linear = .FALSE.
    sym%maxis = .FALSE.
    sym%ogroup = .FALSE.
    sym%sgroup = .FALSE.
    sym%sigmad = .FALSE.
    sym%sigmah = .FALSE.
    sym%sigmav = .FALSE.
    sym%tgroup = .FALSE.

    ! Initialize list of equivalent atoms (C1 symmetry)
    natoms = SIZE(sym%nequatom)
    sym%ngroup = natoms
    sym%nequatom(:) = 1
    DO i=1,sym%ngroup
      sym%group_of(i) = i
      sym%llequatom(i) = i
      sym%symequ_list(i) = i
      sym%ulequatom(i) = i
    END DO

    sym%point_group_order = -1

    DO iatom=1,natoms
      sym%aw(iatom) = weight(iatom)
    END DO

    ! Generate atomic identification numbers (symmetry code) ***
    sym%ain(:) = REAL(atype(:),KIND=dp) + sym%aw(:)

    ! Initialize the transformation matrix for input orientation -> standard orientation
    CALL unit_matrix(sym%inptostd(:,:))

  END SUBROUTINE init_symmetry

! *****************************************************************************
!> \brief  Return .TRUE., if the atom atoma is in the list symequ_list.
!> \param atoma ...
!> \param sym ...
!> \retval in_symequ_list ...
!> \par History
!>        Creation (21.04.95, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION in_symequ_list(atoma,sym)
    INTEGER, INTENT(IN)                      :: atoma
    TYPE(molsym_type), INTENT(inout)         :: sym
    LOGICAL                                  :: in_symequ_list

    INTEGER                                  :: iatom

    in_symequ_list = .FALSE.

    DO iatom=1,sym%ulequatom(sym%ngroup)
      IF (sym%symequ_list(iatom) == atoma) THEN
        in_symequ_list = .TRUE.
        RETURN
      END IF
    END DO

  END FUNCTION in_symequ_list

! *****************************************************************************
!> \brief Search for point groups of LOW SYMmetry (Ci,Cs).
!> \param sym ...
!> \param coord ...
!> \par History
!>        Creation (21.04.95, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE lowsym(sym,coord)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord

    CHARACTER(len=*), PARAMETER :: routineN = 'lowsym', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: phi
    REAL(KIND=dp), DIMENSION(3)              :: d

    IF (sym%nsn == 2) THEN

      ! Ci
      sym%invers = .TRUE.
      phi = angle(sym%ses(:,1,2),sym%z_axis(:))
      d(:) = vector_product(sym%ses(:,1,2),sym%z_axis(:))
      CALL rotate_molecule(phi,d(:),sym,coord)

    ELSE IF (sym%nsig == 1) THEN

      ! Cs
      sym%sigmah = .TRUE.
      phi = angle(sym%sig(:,1),sym%z_axis(:))
      d(:) = vector_product(sym%sig(:,1),sym%z_axis(:))
      CALL rotate_molecule(phi,d(:),sym,coord)

    END IF

  END SUBROUTINE lowsym

! *****************************************************************************
!> \brief Main program for symmetry analysis.
!> \param sym ...
!> \param eps_geo ...
!> \param coord ...
!> \param atype ...
!> \param weight ...
!> \par History
!>        Creation (14.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE molecular_symmetry(sym,eps_geo,coord,atype,weight)
    TYPE(molsym_type), POINTER               :: sym
    REAL(KIND=dp), INTENT(IN)                :: eps_geo
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    INTEGER, DIMENSION(:), INTENT(IN)        :: atype
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: weight

    CHARACTER(len=*), PARAMETER :: routineN = 'molecular_symmetry', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, natoms

    CALL timeset(routineN,handle)

    ! Perform memory allocation for the symmetry analysis
    natoms = SIZE(coord,2)
    CALL create_molsym(sym,natoms)
    sym%eps_geo=eps_geo

    ! Initialization of symmetry analysis
    CALL init_symmetry(sym,atype,weight)

    IF (natoms == 1) THEN
      ! Special case: only one atom
      CALL atomsym(sym)
    ELSE
      ! Find molecular symmetry
      CALL moleculesym(sym,coord)
      ! Get point group and load point group table
      CALL get_point_group_symbol(sym)
    END IF

    ! Calculate group order
    IF (.NOT.sym%linear) CALL get_point_group_order(sym)

    ! Generate a list of equivalent atoms
    CALL build_symequ_list(sym,coord)

    CALL timestop(handle)

  END SUBROUTINE molecular_symmetry

! *****************************************************************************
!> \brief Find the molecular symmetry.
!> \param sym ...
!> \param coord ...
!> \par History
!>        Creation (14.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE moleculesym(sym,coord)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord

    CHARACTER(len=*), PARAMETER :: routineN = 'moleculesym', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icn, isec, isn
    LOGICAL                                  :: failed
    REAL(KIND=dp)                            :: eps_tenval

! Tolerance of degenerate eigenvalues for the molecular tensor of inertia

    eps_tenval = 0.01_dp*sym%eps_geo

    ! Calculate the molecular tensor of inertia
    CALL tensor(sym,coord)
    ! Use symmetry information from the eigenvalues of the molecular tensor of inertia
    IF ((sym%tenval(3) - sym%tenval(1)) < eps_tenval) THEN
      ! 0 < tenval(1) = tenval(2) = tenval(3)
      sym%cubic = .TRUE.
    ELSE IF ((sym%tenval(3) - sym%tenval(2)) < eps_tenval) THEN
      ! 0 < tenval(1) < tenval(2) = tenval(3)
      ! special case: 0 = tenval(1) < tenval(2) = tenval(3)
      IF (sym%tenval(1) < eps_tenval) sym%linear = .TRUE.
      CALL tracar(sym,coord,(/3,1,2/))
    ELSE
      ! 0 < tenval(1) = tenval(2) < tenval(3) or 0 < tenval(1) < tenval(2) < tenval(3)
      CALL tracar(sym,coord,(/1,2,3/))
    END IF

    ! Continue with the new coordinate axes
    DO
      failed = .FALSE.
      IF (sym%cubic) THEN
        CALL cubsym(sym,coord,failed)
        IF (failed) THEN
          sym%cubic = .FALSE.
          CYCLE
        END IF

      ELSE IF (sym%linear) THEN

        ! Linear molecule
        IF (saxis(2,sym%z_axis(:),sym,coord)) THEN
          sym%invers = .TRUE.
          sym%dgroup = .TRUE.
          sym%sigmah = .TRUE.
        END IF

      ELSE

        ! Check the new coordinate axes for Cn axes
        DO icn=2,maxcn
          IF (caxis(icn,sym%z_axis(:),sym,coord)) CALL addsec(icn,sym%z_axis(:),sym)
          IF (caxis(icn,sym%x_axis(:),sym,coord)) CALL addsec(icn,sym%x_axis(:),sym)
          IF (caxis(icn,sym%y_axis(:),sym,coord)) CALL addsec(icn,sym%y_axis(:),sym)
        END DO

        ! Check the new coordinate axes for Sn axes
        DO isn=2,maxsn
          IF (saxis(isn,sym%z_axis(:),sym,coord)) CALL addses(isn,sym%z_axis(:),sym)
          IF (saxis(isn,sym%x_axis(:),sym,coord)) CALL addses(isn,sym%x_axis(:),sym)
          IF (saxis(isn,sym%y_axis(:),sym,coord)) CALL addses(isn,sym%y_axis(:),sym)
        END DO

        ! Check the new coordinate planes for mirror planes
        IF (sigma(sym%z_axis(:),sym,coord)) CALL addsig(sym%z_axis(:),sym)
        IF (sigma(sym%x_axis(:),sym,coord)) CALL addsig(sym%x_axis(:),sym)
        IF (sigma(sym%y_axis(:),sym,coord)) CALL addsig(sym%y_axis(:),sym)

        ! There is a main axis (MAXIS = .TRUE.)
        IF ((sym%ncn > 1).OR.(sym%nsn > 3)) THEN
          sym%maxis = .TRUE.
          CALL axsym(coord,sym)
        ELSE
          ! Only low or no symmetry (Ci, Cs or C1)
          CALL lowsym(sym,coord)
        END IF

      END IF

      IF (.NOT.failed) EXIT

    END DO

    ! Find the remaining S axes
    IF (.NOT.sym%linear) THEN
      DO icn=2,sym%ncn
        DO isec=1,sym%nsec(icn)
          IF (saxis(icn,sym%sec(:,isec,icn),sym,coord)) &
                    CALL addses(icn,sym%sec(:,isec,icn),sym)
          IF (saxis(2*icn,sym%sec(:,isec,icn),sym,coord)) &
                    CALL addses(2*icn,sym%sec(:,isec,icn),sym)
        END DO
      END DO
    END IF

    ! Remove redundant S2 axes
    IF (sym%nses(2) > 0) THEN
      sym%nses(1) = sym%nses(1) - sym%nses(2)
      sym%nses(2) = 0
      CALL addses(2,sym%z_axis(:),sym)
    END IF

  END SUBROUTINE moleculesym

! *****************************************************************************
!> \brief The number of atoms which are placed on an axis defined by the vector a is returned.
!> \param a ...
!> \param coord ...
!> \param sym ...
!> \retval naxis ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION naxis(a,coord,sym)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    TYPE(molsym_type), INTENT(inout)         :: sym
    INTEGER                                  :: naxis

    INTEGER                                  :: iatom, natoms
    REAL(KIND=dp)                            :: length_of_a, length_of_b, &
                                                scapro
    REAL(KIND=dp), DIMENSION(3)              :: a_norm, b, b_norm

    naxis = 0
    natoms = SIZE(coord,2)

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

    ! Check the length of vector a
    IF (length_of_a > sym%eps_geo) THEN

      DO iatom=1,natoms
        b(:) = coord(:,iatom)
        length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
        ! An atom in the origin counts for each axis
        IF (length_of_b < sym%eps_geo) THEN
          naxis = naxis + 1
        ELSE
         a_norm = a(:)/length_of_a
         b_norm = b(:)/length_of_b
         scapro = a_norm(1)*b_norm(1) + a_norm(2)*b_norm(2) + a_norm(3)*b_norm(3)
         IF (ABS(ABS(scapro) - 1.0_dp) < sym%eps_geo) naxis = naxis + 1
        END IF
      END DO

    END IF

  END FUNCTION naxis

! *****************************************************************************
!> \brief Return .FALSE., if the symmetry elements se1 and se2 are a pair of
!>         redundant symmetry elements.
!> \param se1 ...
!> \param se2 ...
!> \retval newse ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION newse(se1,se2)
    INTEGER, INTENT(IN)                      :: se1, se2
    LOGICAL                                  :: newse

    INTEGER                                  :: ise

    newse = .TRUE.

    DO ise=2,MIN(se1,se2)
      IF ((MODULO(se1,ise) == 0).AND.(MODULO(se2,ise) == 0)) THEN
        newse = .FALSE.
        RETURN
      END IF
    END DO

  END FUNCTION newse

! *****************************************************************************
!> \brief The number of atoms which are placed in a plane defined by the
!>         the normal vector a is returned.
!> \param a ...
!> \param sym ...
!> \param coord ...
!> \retval nsigma ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION nsigma(a,sym,coord)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    INTEGER                                  :: nsigma

    INTEGER                                  :: iatom, natoms
    REAL(KIND=dp)                            :: length_of_a, length_of_b, &
                                                scapro
    REAL(KIND=dp), DIMENSION(3)              :: a_norm, b, b_norm

    nsigma = 0

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

    ! Check the length of vector a
    IF (length_of_a > sym%eps_geo) THEN
      natoms = SIZE(coord,2)
      DO iatom=1,natoms
        b(:) = coord(:,iatom)
        length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
        ! An atom in the origin counts for each mirror plane
        IF (length_of_b < sym%eps_geo) THEN
          nsigma = nsigma + 1
        ELSE
         a_norm = a(:)/length_of_a
         b_norm = b(:)/length_of_b
         scapro = a_norm(1)*b_norm(1) + a_norm(2)*b_norm(2) + a_norm(3)*b_norm(3)
         IF (ABS(scapro) < sym%eps_geo) nsigma = nsigma + 1
        END IF
      END DO
    END IF

  END FUNCTION nsigma

! *****************************************************************************
!> \brief Style the output of the symmetry elements.
!> \param se ...
!> \param eps ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE outse(se,eps)
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: se
    REAL(KIND=dp), INTENT(IN)                :: eps

! Positive z component for all vectors

    IF (ABS(se(3)) < eps) THEN
      IF (ABS(se(1)) < eps) THEN
        se(2) = -se(2)
      ELSE
        IF (se(1) < 0.0_dp) se(1:2) = -se(1:2)
      END IF
    ELSE
      IF (se(3) < 0.0_dp) se(:) = -se(:)
    END IF

  END SUBROUTINE outse

! *****************************************************************************
!> \brief Print the output of the symmetry analysis.
!> \param sym ...
!> \param coord ...
!> \param atype ...
!> \param element ...
!> \param z ...
!> \param weight ...
!> \param iw ...
!> \param plevel ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE print_symmetry(sym,coord,atype,element,z,weight,iw,plevel)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: coord
    INTEGER, DIMENSION(:), INTENT(in)        :: atype
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: element
    INTEGER, DIMENSION(:), INTENT(in)        :: z
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: weight
    INTEGER, INTENT(IN)                      :: iw, plevel

    CHARACTER(len=*), PARAMETER :: routineN = 'print_symmetry', &
      routineP = moduleN//':'//routineN

    CHARACTER(3)                             :: string
    INTEGER                                  :: i, icn, iequatom, isec, ises, &
                                                isig, isn, j, secount
    REAL(KIND=dp)                            :: null

! Print point group symbol and point group order

    IF (sym%point_group_order == -1) THEN
      WRITE (iw,"(/,T2,A,T40,A,T77,A4)") "Molecular Symmetry Information", &
        "Point Group: ",ADJUSTR(sym%point_group_symbol)
    ELSE
      WRITE (iw,"(/,T2,A,T40,A,T77,A4,/,T40,A,T77,I4)")&
        "Molecular Symmetry Information",&
        "Point Group: ",ADJUSTR(sym%point_group_symbol),&
        "Group Order: ",sym%point_group_order
    END IF

    IF ( MOD(plevel,10) == 1 ) THEN
      ! Print the Cartesian coordinates of the standard orientation
      CALL write_coordinates(coord,atype,element,z,weight,iw)

      WRITE (iw,"(/,T3,A,(T41,10I4))") "Group Number:      1   Group Members:",&
        (sym%symequ_list(iequatom),iequatom=sym%llequatom(1),sym%ulequatom(1))
      DO i=2,sym%ngroup
        WRITE (iw,"(T15,I8,(T41,10I4))")&
          i,(sym%symequ_list(iequatom),iequatom=sym%llequatom(i),sym%ulequatom(i))
      END DO
    END IF

    IF ( MOD(plevel/100,10) == 1 ) THEN
      !   *** Print all symmetry elements ***
      WRITE (iw,"(/,T24,A,9X,A,12X,A,12X,A)") "Symmetry Elements","X","Y","Z"

      secount = 1
      null=0._dp
      WRITE (iw,"(T24,2I5,2X,A1,5x,3F13.6)") secount,secount,"E",null,null,null

      string = "@  "
      DO isig=1,sym%nsig
        secount = secount + 1
        CALL outse(sym%sig(:,isig),sym%eps_geo)
        WRITE (iw,"(T24,2I5,2X,A3,3X,3F13.6)")&
          secount,isig,string,(sym%sig(i,isig),i=1,3)
      END DO

      DO icn=2,sym%ncn
        IF (icn < 10) THEN
          WRITE (string,"(A1,I1)") "C",icn
        ELSE
          WRITE (string,"(A1,I2)") "C",icn
        END IF
        DO isec=1,sym%nsec(icn)
          secount = secount + 1
          CALL outse(sym%sec(:,isec,icn),sym%eps_geo)
          WRITE (iw,"(T24,2I5,2X,A3,3X,3F13.6)")&
            secount,isec,string,(sym%sec(i,isec,icn),i=1,3)
        END DO
      END DO

      DO isn=2,sym%nsn
        IF (isn == 2) THEN
          WRITE (string,"(A3)") "i  "
        ELSE IF (isn < 10) THEN
          WRITE (string,"(A1,I1)") "S",isn
        ELSE
          WRITE (string,"(A1,I2)") "S",isn
        END IF
        DO ises=1,sym%nses(isn)
          secount = secount + 1
          CALL outse(sym%ses(:,ises,isn),sym%eps_geo)
          WRITE (iw,"(T24,2I5,2X,A3,3X,3F13.6)")&
            secount,ises,string,(sym%ses(i,ises,icn),i=1,3)
        END DO
      END DO
    END IF

    IF ( MOD(plevel/10,10) == 1 .AND. SIZE(coord,2)>1 ) THEN
      ! Print the moments of the molecular inertia tensor
      WRITE (iw,"(/,T24,A,/,T40,A,2(16X,A),3(/,T25,A,3F18.5))")&
        "Moments of the Molecular Inertia Tensor","Ix","Iy","Iz",&
        "Ix",(sym%tenmat(1,i),i=1,3),&
        "Iy",(sym%tenmat(2,i),i=1,3),&
        "Iz",(sym%tenmat(3,i),i=1,3)
      WRITE (iw,"(/,T24,A,/,T27,3F18.5,/,(T27,3F18.5))")&
        "Principal Moments and Axes of Inertia",(sym%tenval(i),i=1,3),&
        ((sym%tenvec(i,j),j=1,3),i=1,3)
    END IF

  END SUBROUTINE print_symmetry

! *****************************************************************************
!> \brief Rotate the molecule about an axis defined by the vector a. The
!>        rotation angle is phi (radians).
!> \param phi ...
!> \param a ...
!> \param sym ...
!> \param coord ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE rotate_molecule(phi,a,sym,coord)
    REAL(KIND=dp), INTENT(IN)                :: phi
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord

    CHARACTER(len=*), PARAMETER :: routineN = 'rotate_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: length_of_a

! Check the length of vector a

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
    IF (length_of_a > sym%eps_geo) THEN

      ! Build up the rotation matrix
      CALL build_rotmat(phi,a(:),sym%rotmat(:,:))

      ! Rotate the molecule by phi around a
      coord(:,:) = MATMUL(sym%rotmat(:,:),coord(:,:))

      ! Rotate the normal vectors of the mirror planes which are already found
      sym%sig(:,1:sym%nsig) = MATMUL(sym%rotmat(:,:),sym%sig(:,1:sym%nsig))

      ! Rotate the Cn axes which are already found
      DO i=2,sym%ncn
        sym%sec(:,1:sym%nsec(i),i) = MATMUL(sym%rotmat(:,:),sym%sec(:,1:sym%nsec(i),i))
      END DO

      ! Rotate the Sn axes which are already found
      DO i=2,sym%nsn
        sym%ses(:,1:sym%nses(i),i) = MATMUL(sym%rotmat(:,:),sym%ses(:,1:sym%nses(i),i))
      END DO

      ! Store current rotation
      sym%inptostd(:,:) = MATMUL(sym%rotmat(:,:),sym%inptostd(:,:))

    END IF

  END SUBROUTINE rotate_molecule

! *****************************************************************************
!> \brief Rotate the molecule about a n-fold axis defined by the vector a
!>        using a rotation angle of 2*pi/n. Check, if the axis is a Cn axis.
!> \param n ...
!> \param a ...
!> \param sym ...
!> \param coord ...
!> \retval saxis ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION saxis(n,a,sym,coord)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    LOGICAL                                  :: saxis

    INTEGER                                  :: iatom, natoms
    REAL(KIND=dp)                            :: length_of_a, phi
    REAL(KIND=dp), DIMENSION(3)              :: b

    saxis = .FALSE.

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

    natoms = SIZE(coord,2)

    ! Check the length of the axis vector a
    IF (length_of_a > sym%eps_geo) THEN
      ! Calculate the rotation angle phi and build up the rotation matrix rotmat
      phi = 2.0_dp*pi/REAL(n,KIND=dp)
      CALL build_rotmat(phi,a(:),sym%rotmat(:,:))
      ! Check all atoms
      DO iatom=1,natoms
        ! Rotate the actual atom by 2*pi/n about a
        b(:) = MATMUL(sym%rotmat(:,:),coord(:,iatom))
        ! Reflect the coordinates of the rotated atom
        b(:) = reflect_vector(b(:),a(:))
        ! Check if the coordinates of the rotated atom are in the coordinate set of the molecule
        IF (equatom(iatom,b(:),sym,coord) == 0) RETURN
      END DO
      ! The axis defined by a is a Sn axis, thus return with saxis = .TRUE.
      saxis = .TRUE.
    END IF

  END FUNCTION saxis

! *****************************************************************************
!> \brief Reflect all atoms of the molecule through a mirror plane defined
!>        by the normal vector a.
!> \param a ...
!> \param sym ...
!> \param coord ...
!> \retval sigma ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  FUNCTION sigma(a,sym,coord)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    LOGICAL                                  :: sigma

    INTEGER                                  :: iatom, natoms
    REAL(KIND=dp)                            :: length_of_a
    REAL(KIND=dp), DIMENSION(3)              :: b

    sigma = .FALSE.

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

    ! Check the length of vector a
    IF (length_of_a > sym%eps_geo) THEN
      natoms = SIZE(coord,2)
      DO iatom=1,natoms
        ! Reflect the actual atom
        b(:) = reflect_vector(coord(:,iatom),a(:))
        ! Check if the coordinates of the reflected atom are in the coordinate set of the molecule
        IF (equatom(iatom,b(:),sym,coord) == 0) RETURN
      END DO
      ! The plane defined by a is a mirror plane, thus return with sigma = .TRUE.
      sigma = .TRUE.
    END IF

  END FUNCTION sigma

! *****************************************************************************
!> \brief Calculate the molecular tensor of inertia and the vector to the
!>        center of mass of the molecule.
!> \param sym ...
!> \param coord ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE tensor(sym,coord)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord

    CHARACTER(len=*), PARAMETER :: routineN = 'tensor', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natoms
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tt

    failure = .FALSE.
! Find the vector center_of_mass to molecular center of mass

    natoms = SIZE(coord,2)
    sym%center_of_mass(:) = MATMUL(coord(:,1:natoms),sym%aw(1:natoms))/SUM(sym%aw(1:natoms))

    ! Translate the center of mass of the molecule to the origin
    coord(:,1:natoms) = coord(:,1:natoms) - SPREAD(sym%center_of_mass(:),2,natoms)

    ! Build up the molecular tensor of inertia

    sym%tenmat(1,1) = DOT_PRODUCT(sym%aw(1:natoms),(coord(2,1:natoms)**2 + coord(3,1:natoms)**2))
    sym%tenmat(2,2) = DOT_PRODUCT(sym%aw(1:natoms),(coord(3,1:natoms)**2 + coord(1,1:natoms)**2))
    sym%tenmat(3,3) = DOT_PRODUCT(sym%aw(1:natoms),(coord(1,1:natoms)**2 + coord(2,1:natoms)**2))

    sym%tenmat(1,2) = -DOT_PRODUCT(sym%aw(1:natoms),(coord(1,1:natoms)* coord(2,1:natoms)))
    sym%tenmat(1,3) = -DOT_PRODUCT(sym%aw(1:natoms),(coord(1,1:natoms)* coord(3,1:natoms)))
    sym%tenmat(2,3) = -DOT_PRODUCT(sym%aw(1:natoms),(coord(2,1:natoms)* coord(3,1:natoms)))

    ! Symmetrize tensor matrix
    sym%tenmat(2,1) = sym%tenmat(1,2)
    sym%tenmat(3,1) = sym%tenmat(1,3)
    sym%tenmat(3,2) = sym%tenmat(2,3)

    ! Diagonalize the molecular tensor of inertia
    CALL jacobi(sym%tenmat(:,:),sym%tenval(:),sym%tenvec(:,:))

    ! Secure that the principal axes are right-handed
    sym%tenvec(:,3) = vector_product(sym%tenvec(:,1),sym%tenvec(:,2))

    tt = SQRT(sym%tenval(1)**2 + sym%tenval(2)**2 + sym%tenval(3)**2)
    CPASSERT(tt /= 0)

  END SUBROUTINE tensor

! *****************************************************************************
!> \brief Transformation the Cartesian coodinates with the matrix tenvec.
!>        The coordinate axes can be switched according to the index
!>        vector idx. If idx(1) is equal to 3 instead to 1, then the first
!>        eigenvector (smallest eigenvalue) of the molecular tensor of
!>        inertia is connected to the z axis instead of the x axis. In
!>        addition the local atomic coordinate systems are initialized,
!>        if the symmetry is turned off.
!> \param sym ...
!> \param coord ...
!> \param idx ...
!> \par History
!>        Creation (16.11.98, Matthias Krack)
!> \author jgh
! *****************************************************************************
  SUBROUTINE tracar(sym,coord,idx)
    TYPE(molsym_type), INTENT(inout)         :: sym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: coord
    INTEGER, DIMENSION(3), INTENT(IN)        :: idx

! Transformation of the atomic coordinates  ***

    coord(idx,:) = MATMUL(TRANSPOSE(sym%tenvec(:,:)),coord(:,:))

    ! Modify the transformation matrix for input orientation -> standard orientation
    sym%inptostd(idx,:) = TRANSPOSE(sym%tenvec(:,:))

  END SUBROUTINE tracar

! *****************************************************************************
!> \brief Distance between two points
!> \param a ...
!> \param b ...
!> \retval d ...
!> \author jgh
! *****************************************************************************
  FUNCTION dist(a,b) RESULT(d)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp)                            :: d

    d = SQRT(SUM((a-b)**2))

  END FUNCTION
! *****************************************************************************
!> \brief   Write atomic coordinates to output unit iw.
!> \param coord ...
!> \param atype ...
!> \param element ...
!> \param z ...
!> \param weight ...
!> \param iw ...
!> \date    08.05.2008
!> \author  JGH
! *****************************************************************************
  SUBROUTINE write_coordinates(coord,atype,element,z,weight,iw)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: coord
    INTEGER, DIMENSION(:), INTENT(in)        :: atype
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: element
    INTEGER, DIMENSION(:), INTENT(in)        :: z
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: weight
    INTEGER, INTENT(IN)                      :: iw

    INTEGER                                  :: iatom, natom

    IF (iw>0) THEN

      WRITE (UNIT=iw,FMT="(/,T2,A)")&
        "Atomic Coordinates [Bohr] in Standard Orientation"
      WRITE (UNIT=iw,&
           FMT="(/,T3,A,7X,2(A1,11X),A1,6X,A10,4X,A4,/)")&
        "Atom  Kind  Element","X","Y","Z","Atomic-No.","Mass"

      natom = SIZE(coord,2)
      DO iatom=1,natom
        WRITE (UNIT=iw,FMT="(T2,I5,1X,I4,3X,A2,5X,3F12.6,4X,I4,2X,F11.4)")&
          iatom,atype(iatom),element(iatom),coord(1:3,iatom),z(iatom),weight(iatom)
      END DO

    END IF

  END SUBROUTINE write_coordinates
! *****************************************************************************

END MODULE molsym
