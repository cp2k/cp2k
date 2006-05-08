!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****m* cp2k/qs_electric_field_gradient *
!!
!!   NAME
!!     qs_electric_field_gradient
!!
!!   FUNCTION
!!     Calculates electric field gradients
!!     H.M. Petrili, P.E. Blochl, P. Blaha, K. Schwarz, PRB 57, 14690 (1998)
!!
!!   NOTE
!!  
!!   AUTHOR
!!     JGH (03-05-2006)
!!
!!   MODIFICATION HISTORY
!!     none
!!  
!!****

MODULE qs_electric_field_gradient
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
! USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE eigenvalueproblems,              ONLY: diagonalise
! USE input_section_types,             ONLY: section_vals_get_subs_vals,&
!                                            section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac,&
                                             fourpi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: a_bohr,&
                                             e_charge,&
                                             joule
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_dr2,&
                                             pw_integral_ab,&
                                             pw_structure_factor
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE util,                            ONLY: get_limit,&
                                             sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qs_efg_calc

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electric_field_gradient'

!*******************************************************************************

CONTAINS

!*******************************************************************************
  
  SUBROUTINE qs_efg_calc(qs_env,unit_nr,interpolation,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: unit_nr
    LOGICAL, INTENT(IN)                      :: interpolation
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efg_calc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: i, iat, iatom, ij, ikind, &
                                                istat, j, m, natom, &
                                                natomkind, nkind, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, gapw, paw_atom
    REAL(KIND=dp)                            :: efg_units, ehartree, f1, f2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: efg_diagval, vh0
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: efg_tensor
    REAL(KIND=dp), DIMENSION(3)              :: eigenvalues, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: eigenvectors
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(coeff_type)                         :: rho_tot_gspace, &
                                                structure_factor, &
                                                v_hartree_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: dvr2
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

    NULLIFY(atomic_kind_set, dft_control, para_env, particle_set, rho, &
         rho_atom_set)

    IF (unit_nr > 0) THEN
      WRITE(unit_nr,"(/,A,/)") " ELECTRIC FIELD GRADIENTS [V/m^2]"
    END IF

    efg_units = Joule/a_bohr**2/e_charge * 1.e-21_dp
    f1 = SQRT(15._dp/fourpi)
    f2 = SQRT(5._dp/fourpi)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    rho=rho,atomic_kind_set=atomic_kind_set,&
                    rho_atom_set=rho_atom_set,pw_env=pw_env,&
                    particle_set=particle_set,para_env=para_env,&
                    error=error)

    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins

    natom = SIZE(particle_set,1)
    ALLOCATE (efg_tensor(3,3,natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    efg_tensor = 0._dp
    ALLOCATE (efg_diagval(3,natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    efg_diagval = 0._dp

    ALLOCATE (vh0(natom,-2:2),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (dvr2(6),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    !prepare calculation
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                    poisson_env=poisson_env,error=error)
    IF (gapw) CALL prepare_gapw_den(qs_env,do_rho0=.TRUE.,error=error)

    !calculate electrostatic potential
    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_hartree_gspace, &
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=rho_tot_gspace,&
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

    CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,ehartree,&
                          v_hartree_gspace%pw)

    CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_gspace,error=error)

    DO i=1,3
      DO j=1,i
        ij=(i*(i-1))/2+j
        CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=dvr2(ij), &
                                use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                error=error)
        CALL pw_dr2 ( v_hartree_gspace%pw, dvr2(ij)%pw, i, j )
      END DO
    END DO
    IF ( .NOT. interpolation ) THEN
      CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=structure_factor, &
                              use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                              error=error)
    END IF

    nkind = SIZE(atomic_kind_set) 

    DO ikind = 1,nkind 
      atomic_kind => atomic_kind_set(ikind)
      NULLIFY(grid_atom,atom_list)
      CALL get_atomic_kind(atomic_kind=atomic_kind,grid_atom=grid_atom,&
                           atom_list=atom_list,natom=natomkind,paw_atom=paw_atom)
      IF( paw_atom) THEN
        CALL vlimit_atom(para_env,vh0,rho_atom_set,atomic_kind,&
                         atom_list,grid_atom,natomkind,nspins,error)
      END IF

      DO iat = 1,natomkind 
        iatom = atom_list(iat)
        IF ( interpolation ) THEN
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        ELSE
          CALL pw_structure_factor(structure_factor%pw,ra)
          DO i=1,3
            DO j=1,i
              ij=(i*(i-1))/2+j
              efg_tensor(i,j,iatom)=-pw_integral_ab(dvr2(ij)%pw,structure_factor%pw)
              efg_tensor(j,i,iatom)=efg_tensor(i,j,iatom)
            END DO
          END DO
          efg_tensor(:,:,iatom)=efg_tensor(:,:,iatom)/structure_factor%pw%pw_grid%vol
        END IF
        IF(paw_atom) THEN
          efg_tensor(1,1,iatom)=efg_tensor(1,1,iatom)+f1*(vh0(iatom,2))-f2*(vh0(iatom,0))
          efg_tensor(2,2,iatom)=efg_tensor(2,2,iatom)-f1*(vh0(iatom,2))-f2*(vh0(iatom,0))
          efg_tensor(3,3,iatom)=efg_tensor(3,3,iatom)+2._dp*f2*(vh0(iatom,0))
          efg_tensor(1,2,iatom)=efg_tensor(1,2,iatom)+f1*(vh0(iatom,-2))
          efg_tensor(2,1,iatom)=efg_tensor(2,1,iatom)+f1*(vh0(iatom,-2))
          efg_tensor(1,3,iatom)=efg_tensor(1,3,iatom)+f1*(vh0(iatom,1))
          efg_tensor(3,1,iatom)=efg_tensor(3,1,iatom)+f1*(vh0(iatom,1))
          efg_tensor(2,3,iatom)=efg_tensor(2,3,iatom)+f1*(vh0(iatom,-1))
          efg_tensor(3,2,iatom)=efg_tensor(3,2,iatom)+f1*(vh0(iatom,-1))
        END IF
        CALL diagonalise ( efg_tensor(:,:,iatom), 3, "UPPER", eigenvalues, eigenvectors )
        CALL efgsort(eigenvalues,efg_diagval(:,iatom))
      END DO ! iat
    END DO ! ikind

    efg_tensor(:,:,:) = efg_tensor(:,:,:)*efg_units
    efg_diagval(:,:) = efg_diagval(:,:)*efg_units

    DO iatom=1,natom
      IF (unit_nr > 0) THEN
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                            element_symbol=element_symbol)
          WRITE (UNIT=unit_nr,FMT="(T2,I5,T8,A,T12,A,3(T39,3F14.7,/))")&
           iatom, element_symbol,"EFG Tensor",(efg_tensor(i,:,iatom),i=1,3)
          WRITE (UNIT=unit_nr,FMT="(T12,A,T39,3F14.7)")&
           "EFG Tensor eigenvalues",efg_diagval(:,iatom)
          WRITE (UNIT=unit_nr,FMT="(T12,A,T67,F14.7)") "EFG Tensor anisotropy",&
             (efg_diagval(1,iatom)-efg_diagval(2,iatom))/efg_diagval(3,iatom)
          WRITE (UNIT=unit_nr,FMT=*) 
      END IF
    ENDDO

    CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,error=error)
    IF ( .NOT. interpolation ) THEN
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,structure_factor,error=error)
    END IF

    DO i=1,6
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,dvr2(i),error=error)
    END DO

    DEALLOCATE (efg_tensor,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (vh0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (dvr2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!deb
    WRITE(6,*) "END EFG CALCULATION"
!deb

  END SUBROUTINE qs_efg_calc

!*******************************************************************************
  SUBROUTINE vlimit_atom(para_env,vlimit,rho_atom_set,atom_kind,&
                         atom_list,grid_atom,natom,nspins,error)

  ! calculate : Limit(r->0) V_hartree(r)/r^2

    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:, -2:), &
      INTENT(inout)                          :: vlimit
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    INTEGER, DIMENSION(:), INTENT(IN)        :: atom_list
    TYPE(grid_atom_type), INTENT(IN)         :: grid_atom
    INTEGER, INTENT(IN)                      :: natom, nspins
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vlimit_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, i1, i2, iat, iatom, icg, ipgf1, ipgf2, iset1, iset2, iso, &
      iso1, iso1_first, iso1_last, iso2, iso2_first, iso2_last, istat, l, &
      l_iso, m1s, m2s, m_iso, max_iso_not0, maxl, maxso, mepos, n1s, n2s, &
      nset, num_pe, size1, size2
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf
    LOGICAL                                  :: failure
    REAL(dp)                                 :: zet12
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: CPC_sphere
    REAL(dp), DIMENSION(20)                  :: vgg
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics

! *****************************************************************************

    failure = .FALSE.
    NULLIFY(orb_basis)
    NULLIFY(harmonics)
    NULLIFY(lmin,lmax,npgf,zet,my_CG,coeff)

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         harmonics=harmonics)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxl=maxl,npgf=npgf,nset=nset,zet=zet,&
                           maxso=maxso)

    max_iso_not0 = harmonics%max_iso_not0

!   Distribute the atoms of this kind
    num_pe = para_env%num_pe
    mepos  = para_env%mepos
    bo = get_limit( natom, num_pe, mepos )

    my_CG        => harmonics%my_CG

    ALLOCATE(CPC_sphere(nsoset(maxl),nsoset(maxl)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!******************************************************************************

    m1s = 0

    DO iset1 = 1,nset
      n1s = nsoset(lmax(iset1))
      DO ipgf1  = 1,npgf(iset1)

        iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
        iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
        size1 = iso1_last - iso1_first + 1
        i1 = nsoset(lmin(iset1)-1)+1

        m2s = 0
        DO iset2 = 1,nset
          n2s=nsoset(lmax(iset2))
          DO ipgf2 = 1,npgf(iset2)
            iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
            iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
            size2 = iso2_last - iso2_first + 1
            i2 = nsoset(lmin(iset2)-1)+1

            zet12 = zet(ipgf1,iset1)+zet(ipgf2,iset2)

            vgg = 0.0_dp
           
            DO iso = 1,max_iso_not0
               l_iso = indso(1,iso)
               IF ( l_iso /= 2 ) CYCLE
               DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                 iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                 iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                 l = indso(1,iso1) + indso(1,iso2)
                 IF ( MOD(l,2) == 0 .AND. l > 0 ) THEN
                   vgg(l/2) = fourpi/10._dp * fac(l-2) / zet12**(l/2)
                 END IF
               END DO
            END DO

            DO iat = bo(1),bo(2)
              iatom = atom_list(iat)

              CPC_sphere = 0.0_dp
              DO i=1,nspins
                coeff => rho_atom_set(iatom)%cpc_h(i)%r_coef
                CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) = &
                  CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) + &
                  coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                coeff => rho_atom_set(iatom)%cpc_s(i)%r_coef
                CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) = &
                  CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) - &
                  coeff(iso1_first:iso1_last,iso2_first:iso2_last)
               ENDDO ! i

               DO iso = 1,max_iso_not0
                 l_iso = indso(1,iso) 
                 m_iso = indso(2,iso) 
                 IF ( l_iso /= 2 ) CYCLE
                 DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                   iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                   iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)

                   l = indso(1,iso1) + indso(1,iso2)
                   IF ( MOD(l,2) == 0 .AND. l > 0 ) THEN
                     vlimit(iatom,m_iso) = vlimit(iatom,m_iso) +&
                       vgg(l/2)*CPC_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)
                   END IF
                 ENDDO  ! icg
               ENDDO  ! iso 

            ENDDO  ! iat
        
          ENDDO ! ipgf2
          m2s = m2s+maxso

        ENDDO  ! iset2
      ENDDO   ! ipgf1
      m1s = m1s+maxso
    ENDDO    ! iset1

    CALL mp_sum(vlimit,para_env%group)

    DEALLOCATE(CPC_sphere,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE vlimit_atom

!*******************************************************************************
  SUBROUTINE efgsort(ein,eout)
    REAL(dp), DIMENSION(3), INTENT(in)       :: ein
    REAL(dp), DIMENSION(3), INTENT(inout)    :: eout

    INTEGER                                  :: i
    INTEGER, DIMENSION(3)                    :: ind
    REAL(dp), DIMENSION(3)                   :: eab

    eab = ABS(ein)
    CALL sort(eab,3,ind)
    DO i=1,3
      eout(i) = ein(ind(i))
    END DO
    
  END SUBROUTINE efgsort
!*******************************************************************************

END MODULE  qs_electric_field_gradient

!*******************************************************************************
