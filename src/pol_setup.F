!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/pol_setup [1.0] *
!!
!!   NAME
!!     pol_setup
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria 080501
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!
!!   SOURCE
!******************************************************************************

MODULE pol_setup
  
  USE atomic_kinds,                    ONLY: kind_info_type
  USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                             deallocate_gto_basis_set,&
                                             gto_basis_set_type,&
                                             read_gto_basis_set
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE molecule_types,                  ONLY: molecule_type,&
                                             read_basis_type
  USE orbital_pointers,                ONLY: coset,&
                                             init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_interactions,                 ONLY: exp_radius
  USE string_utilities,                ONLY: str_search
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PUBLIC :: read_pol_basis_sets, get_rcutsq_cgf

  PRIVATE


  CONTAINS

!******************************************************************************
 SUBROUTINE read_pol_basis_sets(rho0_basis_info, drho_basis_info, &
                 mol_setup,part,atom_names,ncoef,globenv)

    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: rho0_basis_info, &
                                                drho_basis_info
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: mol_setup
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    CHARACTER(LEN=*), DIMENSION(:), POINTER  :: atom_names
    integer, intent(out)                     :: ncoef
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, icount, imol, isos, ndrho, &
                                                nrho0
    LOGICAL                                  :: flag
    TYPE(read_basis_type), DIMENSION(:), &
      POINTER                                :: pol_basis, rho0_basis

!    TYPE (ao_type) :: ao
!   count the total number of basis functions for 
!   rho0 and polarization density

    nrho0 = 0
    ndrho  = 0
    icount = 0
    
    DO imol=1, SIZE(mol_setup)
      nrho0 = nrho0 + mol_setup(imol)%molpar%nbasis_rho0
      ndrho = ndrho + mol_setup(imol)%molpar%nbasis_pol
    END DO

    ALLOCATE (pol_basis(ndrho), STAT = isos)

    DO imol=1, SIZE(mol_setup)
       DO i = 1, mol_setup(imol)%molpar%nbasis_pol
          icount = icount + 1
          pol_basis(icount)%aname = mol_setup(imol)%molpar%pol_basis(i)%aname
          pol_basis(icount)%libname = mol_setup(imol)%molpar%pol_basis(i)%libname
          pol_basis(icount)%eps = mol_setup(imol)%molpar%pol_basis(i)%eps
          pol_basis(icount)%nmax = mol_setup(imol)%molpar%pol_basis(i)%nmax
       END DO
    END DO
    
    ALLOCATE (rho0_basis(nrho0), STAT = isos)
    
    icount = 0
    DO imol=1, SIZE(mol_setup)
       DO i = 1, mol_setup(imol)%molpar%nbasis_rho0
          icount = icount + 1
          rho0_basis(icount)%aname = mol_setup(imol)%molpar%rho0_basis(i)%aname
          rho0_basis(icount)%libname = mol_setup(imol)%molpar%rho0_basis(i)%libname
          rho0_basis(icount)%eps = mol_setup(imol)%molpar%rho0_basis(i)%eps
          rho0_basis(icount)%nmax = mol_setup(imol)%molpar%rho0_basis(i)%nmax
       END DO
    END DO
     
!   allocate the pointer drho_basis_info
    NULLIFY (drho_basis_info)
    CALL  allocate_pol_basis_info(drho_basis_info,ndrho)

!   allocate the pointer rho0_basis_info
    NULLIFY (rho0_basis_info)
    CALL  allocate_pol_basis_info(rho0_basis_info,nrho0)

!   fill the pointer drho_basis_info
    flag = .TRUE.
    CALL  initialize_pol_basis_info(drho_basis_info,ndrho,pol_basis,part,atom_names, &
                                            flag, globenv)

!   fill the pointer rho0_basis_info
    flag = .FALSE.
    CALL  initialize_pol_basis_info(rho0_basis_info,nrho0,rho0_basis,part,atom_names, &
                                            flag, globenv)

!   count the total number of basis functions (ncoef) needed to allocate
!   the coefficents vectors
    CALL get_number_of_coefs(drho_basis_info, part, ncoef)

  END SUBROUTINE read_pol_basis_sets

!   ***************************************************************************

  SUBROUTINE allocate_pol_basis_info(ki,number_of_kinds)

    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: ki
    INTEGER, INTENT(IN)                      :: number_of_kinds

    INTEGER                                  :: ikind, istat, nkind

!   ---------------------------------------------------------------------------
!   *** Deallocate the old data structure for the atomic kGind information ***

    IF (ASSOCIATED(ki)) THEN
      nkind = size(ki)
      DO ikind=1,nkind
        IF (ASSOCIATED(ki(ikind)%aux_basis_set)) THEN
          CALL deallocate_gto_basis_set(ki(ikind)%aux_basis_set)
        END IF
        IF (ASSOCIATED(ki(ikind)%orb_basis_set)) THEN
          CALL deallocate_gto_basis_set(ki(ikind)%orb_basis_set)
        END IF
        DEALLOCATE (ki(ikind)%atom_list,STAT=istat)
        IF (istat /= 0) CALL stop_memory("allocate_pol_basis_info", &
                    "ki(ikind)%atom_list")
        DEALLOCATE (ki(ikind)%elec_conf,STAT=istat)
        IF (istat /= 0) CALL stop_memory("allocate_pol_basis_info","ki(ikind)%elec_conf")
      END DO
      DEALLOCATE (ki,STAT=istat)
      IF (istat /= 0) CALL stop_memory("allocate_pol_basis_info","ki")
    END IF

    nkind = number_of_kinds

    IF (nkind >= 0) THEN
      ALLOCATE (ki(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory("allocate_pol_basis_info","ki(nkind)",0)
      DO ikind=1,nkind
        NULLIFY (ki ( ikind ) % aux_basis_set)
        NULLIFY (ki ( ikind ) % orb_basis_set)
        NULLIFY (ki ( ikind ) % atom_list)
        NULLIFY (ki ( ikind) % elec_conf)
      END DO
    ELSE
      CALL stop_program("allocate_pol_basis_info",&
                        "invalid number of atomic kinds")
    END IF

  END SUBROUTINE allocate_pol_basis_info

! *****************************************************************************

  SUBROUTINE initialize_pol_basis_info(ki,basis_size,molinfo,part,atom_names,&
                                       flag,globenv)
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: ki
    INTEGER                                  :: basis_size
    TYPE(read_basis_type), DIMENSION(:), &
      POINTER                                :: molinfo
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    CHARACTER(LEN=*), DIMENSION(:), POINTER  :: atom_names
    LOGICAL, INTENT(IN)                      :: flag
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=20)                        :: string
    INTEGER                                  :: ikind, itype, nt

    nt = size (atom_names)

    DO ikind= 1, size(ki)

      ki(ikind)%element_symbol = molinfo(ikind)%aname
      ki(ikind)%orb_basis_set_name = molinfo(ikind)%libname
      ki(ikind)%number_of_grid_points = molinfo(ikind)%nmax
      
      string=molinfo(ikind)%aname
      itype = str_search(atom_names,nt,string)  
 
      CALL get_atom_list(ki(ikind),part,itype,ikind,flag)

      CALL allocate_gto_basis_set(ki(ikind)%orb_basis_set)
 
      CALL read_gto_basis_set(ki(ikind)%element_symbol,ki(ikind)%orb_basis_set_name, &
                        ki(ikind)%orb_basis_set, globenv)

      CALL init_cphi(ki(ikind)% orb_basis_set)

      CALL get_radii(ki(ikind), moliNFo(ikind)%eps)

    ENDDO

  END SUBROUTINE initialize_pol_basis_info                                           
  
!------------------------------------------------------------------------------!
  SUBROUTINE get_atom_list(ki,part,itype,ikind,flag)
  

    TYPE(kind_info_type), INTENT(inout)      :: ki
    TYPE(particle_type), DIMENSION(:), &
      INTENT(inout)                          :: part
    INTEGER, INTENT(in)                      :: itype, ikind
    LOGICAL, INTENT(IN)                      :: flag

    INTEGER                                  :: i, icount, isos

    icount=0
    DO i=1, SIZE(part)
      IF (part(i)%prop%ptype==itype) icount=icount+1
    ENDDO
    ki%natom = icount
    ALLOCATE(ki%atom_list(icount), STAT=isos)
    IF (isos/=0) CALL stop_memory('get_atom_list', 'atom_list')
    icount=0
    DO i=1, SIZE(part)
      IF (part(i)%prop%ptype==itype) THEN
        icount=icount+1
        ki%atom_list(icount)=i
        IF (flag) THEN
          part(i) % kind = ikind
        END IF
      ENDIF
    ENDDO
  END SUBROUTINE get_atom_list                                
!------------------------------------------------------------------------------!
  SUBROUTINE get_radii(ki,eps)

   
    TYPE(kind_info_type), INTENT(inout)      :: ki
    REAL(dbl), INTENT(inout)                 :: eps

    INTEGER                                  :: ipgf, iset, ishell, l
    REAL(dbl)                                :: gcca, kind_radius, &
                                                pgf_radius, set_radius, zeta

!---------------------------------------------------------------------------

  kind_radius = 0.0_dbl
  DO iset=1,ki%orb_basis_set%nset
    set_radius = 0.0_dbl
    DO ipgf=1,ki%orb_basis_set%npgf(iset)
      pgf_radius = 0.0_dbl
      DO ishell=1,ki%orb_basis_set%nshell(iset)
        l = ki%orb_basis_set%l(ishell,iset)
        gcca = ki%orb_basis_set%gcc(ipgf,ishell,iset)
        zeta = ki%orb_basis_set%zet(ipgf,iset)
        pgf_radius = MAX(pgf_radius,exp_radius(l,zeta,eps,gcca))
      END DO
      ki%orb_basis_set%pgf_radius(ipgf,iset) = pgf_radius
      set_radius = MAX(set_radius,pgf_radius)
    END DO
    ki%orb_basis_set%set_radius(iset) = set_radius
    kind_radius = MAX(kind_radius,set_radius)
  END DO
  ki%orb_basis_set%kind_radius = kind_radius
 
  END SUBROUTINE get_radii
!------------------------------------------------------------------------------  

  SUBROUTINE get_number_of_coefs(ki, part, ncoefs)
  
  
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(in)                             :: ki
    TYPE(particle_type), DIMENSION(:), &
      INTENT(inout)                          :: part
    INTEGER, INTENT(OUT)                     :: ncoefs

    INTEGER                                  :: i, icgf, icount, ii, ikind, &
                                                isos, ncgf

    ncoefs = 0
    icount = 0
    
    DO ikind = 1, size(ki)

      DO ii = 1, ki(ikind)%natom
      
        i = ki(ikind)%atom_list(ii)

        ncgf =  ki(ikind)% orb_basis_set % ncgf
        ncoefs = ncoefs + ki(ikind)% orb_basis_set % ncgf

        IF (.NOT.ASSOCIATED (part (i) % coef_list)) THEN
          ALLOCATE(part (i) %coef_list(ncgf), STAT=isos)
            IF (isos/=0) CALL stop_memory('get_number_of_coef', 'coeflist')
         END IF

        DO icgf = 1, ncgf
          icount = icount + 1
          part (i) % coef_list(icgf)=icount
       ENDDO                          
       
    END DO
    
  END DO
  
  END SUBROUTINE get_number_of_coefs
!------------------------------------------------------------------------------ 
  SUBROUTINE get_rcutsq_cgf(ki,rsq)

   
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(in)                             :: kI
    REAL(dbl), DIMENSION(:, :), INTENT(out)  :: rsq

    INTEGER                                  :: i, j
    REAL(dbl)                                :: rcut

! calculate the cutoff radius for the Gaussians needed for the overlap list

    DO i = 1, size(ki)
    
       DO j = 1, size (ki)

          rcut = ki(i)%orb_basis_set%kind_radius + &
                  ki(j)%orb_basis_set%kind_radius 
          rsq (i,j) = rcut * rcut
           
       END DO

    END DO

  END SUBROUTINE get_rcutsq_cgf

! *****************************************************************************

  SUBROUTINE init_cphi(orb_basis_set)

!   Purpose: Initialize the matrices for the transformation of primitive
!            Cartesian Gaussian-type functions to contracted Cartesian
!            (cphi)  Gaussian-type functions.

!   ***************************************************************************
 
 
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    INTEGER                                  :: icgf, ico, ipgf, iset, &
                                                ishell, n

!   ---------------------------------------------------------------------------
!   *** Build the Cartesian transformation matrix "cphi" ***

    DO iset=1,orb_basis_set%nset
      CALL init_orbital_pointers(orb_basis_set%lmax(iset)+1)
      n = ncoset(orb_basis_set%lmax(iset))
      DO ishell=1,orb_basis_set%nshell(iset)
        DO icgf=orb_basis_set%first_cgf(ishell,iset),&
                orb_basis_set%last_cgf(ishell,iset)
          ico = coset(orb_basis_set%lx(icgf),&
                      orb_basis_set%ly(icgf),&
                      orb_basis_set%lz(icgf))
          DO ipgf=1,orb_basis_set%npgf(iset)
            orb_basis_set%cphi(ico,icgf) = orb_basis_set%gcc(ipgf,ishell,iset)
            ico = ico + n
          END DO
        END DO
      END DO
    END DO
 
  END SUBROUTINE init_cphi

!------------------------------------------------------------------------------!

   END MODULE pol_setup
!------------------------------------------------------------------------------!
