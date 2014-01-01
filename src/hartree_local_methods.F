!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE  hartree_local_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_array_init
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hartree_local_types,             ONLY: allocate_ecoul_1center,&
                                             ecoul_1center_type,&
                                             hartree_local_type,&
                                             set_ecoul_1c
  USE input_constants,                 ONLY: tddfpt_singlet,&
                                             use_periodic
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             pi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                             harmonics_atom_type
  USE qs_local_rho_types,              ONLY: rhoz_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_atom_type,&
                                             rho0_mpole_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hartree_local_methods'

  ! Public Subroutine

  PUBLIC :: init_coulomb_local, calculate_Vh_1center, Vh_1c_gg_integrals

CONTAINS

! *****************************************************************************
  SUBROUTINE init_coulomb_local(hartree_local,natom,error)

    TYPE(hartree_local_type), POINTER        :: hartree_local
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_coulomb_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c

    CALL timeset(routineN,handle)

    NULLIFY(ecoul_1c)
    !   Allocate and Initialize 1-center Potentials and Integrals
    CALL allocate_ecoul_1center(ecoul_1c,natom, error)
    hartree_local%ecoul_1c => ecoul_1c

    CALL timestop(handle)

  END SUBROUTINE init_coulomb_local

! *****************************************************************************
!> \brief Calculates Hartree potential for hard and soft densities (including
!>        nuclear charge and compensation charges) using numerical integration
!> \par History
!>      05.2012 JGH refactoring
!> \author ??
! *****************************************************************************
  SUBROUTINE calculate_Vh_1center(vrad_h,vrad_s,rrad_h,rrad_s,rrad_0,rrad_z,grid_atom)

    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: vrad_h, vrad_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      INTENT(IN)                             :: rrad_h, rrad_s
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: rrad_0
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rrad_z
    TYPE(grid_atom_type), POINTER            :: grid_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_Vh_1center', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ir, iso, ispin, &
                                                l_ang, max_s_harm, nchannels, &
                                                nr, nspins
    REAL(dp)                                 :: I1_down, I1_up, I2_down, &
                                                I2_up, prefactor
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rho_1, rho_2
    REAL(dp), DIMENSION(:), POINTER          :: wr
    REAL(dp), DIMENSION(:, :), POINTER       :: oor2l, r2l

    CALL timeset(routineN,handle)

    nr = grid_atom%nr
    max_s_harm = SIZE(vrad_h,2)
    nspins = SIZE(rrad_h,1)
    nchannels = SIZE(rrad_0,2)

    r2l => grid_atom%rad2l
    oor2l => grid_atom%oorad2l
    wr => grid_atom%wr

    ALLOCATE(rho_1(nr,max_s_harm),rho_2(nr,max_s_harm))
    rho_1 = 0.0_dp
    rho_2 = 0.0_dp

    !   Case lm = 0
    rho_1(:,1) = rrad_z(:)
    rho_2(:,1) = rrad_0(:,1)

    DO iso = 2,nchannels
       rho_2(:,iso) = rrad_0(:,iso)
    END DO

    DO iso = 1,max_s_harm
       DO ispin = 1,nspins
          rho_1(:,iso) = rho_1(:,iso) + rrad_h(ispin)%r_coef(:,iso)
          rho_2(:,iso) = rho_2(:,iso) + rrad_s(ispin)%r_coef(:,iso)
       END DO

       l_ang = indso(1,iso)
       prefactor = fourpi/(2._dp*l_ang+1._dp)

       rho_1(:,iso) = rho_1(:,iso)*wr(:)
       rho_2(:,iso) = rho_2(:,iso)*wr(:)

       I1_up = 0.0_dp
       I1_down = 0.0_dp
       I2_up = 0.0_dp
       I2_down = 0.0_dp

       I1_up = r2l(nr,l_ang)*rho_1(nr,iso)
       I2_up = r2l(nr,l_ang)*rho_2(nr,iso)

       DO ir = nr-1,1,-1
          I1_down = I1_down + oor2l(ir,l_ang+1)*rho_1(ir,iso)
          I2_down = I2_down + oor2l(ir,l_ang+1)*rho_2(ir,iso)
       END DO

       vrad_h(nr,iso) = vrad_h(nr,iso) + prefactor*&
            (oor2l(nr,l_ang+1)*I1_up + r2l(nr,l_ang)*I1_down)
       vrad_s(nr,iso) = vrad_s(nr,iso) + prefactor*&
            (oor2l(nr,l_ang+1)*I2_up + r2l(nr,l_ang)*I2_down)

       DO ir = nr-1,1,-1
          I1_up = I1_up + r2l(ir,l_ang)*rho_1(ir,iso)
          I1_down = I1_down -oor2l(ir,l_ang+1)*rho_1(ir,iso)
          I2_up = I2_up + r2l(ir,l_ang)*rho_2(ir,iso)
          I2_down = I2_down -oor2l(ir,l_ang+1)*rho_2(ir,iso)

          vrad_h(ir,iso) = vrad_h(ir,iso) + prefactor*&
               (oor2l(ir,l_ang+1)*I1_up + r2l(ir,l_ang)*I1_down)
          vrad_s(ir,iso) = vrad_s(ir,iso) + prefactor*&
               (oor2l(ir,l_ang+1)*I2_up + r2l(ir,l_ang)*I2_down)

       END DO

    END DO

    DEALLOCATE(rho_1,rho_2)

    CALL timestop(handle)

  END SUBROUTINE calculate_Vh_1center

! *****************************************************************************
!> \brief Calculates one center GAPW Hartree energies and matrix elements
!>        Hartree potentials are input
!>        Takes possible background charge into account
!>        Special case for densities without core charge
!> \par History
!>      05.2012 JGH refactoring
!> \author ??
! *****************************************************************************
  SUBROUTINE Vh_1c_gg_integrals(qs_env,energy_hartree_1c,tddft,do_triplet,p_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(kind=dp), INTENT(out)               :: energy_hartree_1c
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft, do_triplet
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: p_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Vh_1c_gg_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2), handle, iat, iatom, ikind, ipgf1, is1, iset1, iso, &
      istat, l_ang, llmax, lmax0, lmax_0, m1, max_iso, max_iso_not0, &
      max_s_harm, maxl, maxso, mepos, n1, nat, natom, nchan_0, nkind, nr, &
      nset, nsotot, nspins, num_pe
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cg_n_list
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cg_list
    INTEGER, DIMENSION(:), POINTER           :: atom_list, lmax, lmin, npgf
    LOGICAL                                  :: atenergy, core_charge, &
                                                failure, my_periodic, &
                                                my_tddft, paw_atom
    REAL(dp)                                 :: back_ch, factor
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: gexp, sqrtwr
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: aVh1b_00, aVh1b_hh, aVh1b_ss, &
                                                g0_h_w
    REAL(dp), DIMENSION(:), POINTER          :: rrad_z, vrrad_z
    REAL(dp), DIMENSION(:, :), POINTER       :: g0_h, gsph, rrad_0, Vh1_h, &
                                                Vh1_s, vrrad_0, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG, Qlm_gg
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(rho0_atom_type), DIMENSION(:), &
      POINTER                                :: rho0_atom_set
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(rhoz_type), DIMENSION(:), POINTER   :: rhoz_set

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(cell,dft_control,para_env,poisson_env,pw_env,qs_charges)
    NULLIFY(atomic_kind_set, rho_atom_set, rho0_atom_set)
    NULLIFY(rho0_mpole,rhoz_set,ecoul_1c)
    NULLIFY(atom_kind,atom_list,grid_atom,harmonics)
    NULLIFY(orb_basis,lmin,lmax,npgf,zet)
    NULLIFY(gsph)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft
    core_charge = .NOT.my_tddft

    CALL get_qs_env(qs_env=qs_env,&
         cell=cell,dft_control=dft_control,&
         para_env=para_env,&
         atomic_kind_set=atomic_kind_set,&
         pw_env=pw_env,qs_charges=qs_charges,error=error)

    CALL pw_env_get(pw_env,poisson_env=poisson_env,error=error)
    my_periodic= (poisson_env%method==use_periodic)

    back_ch = qs_charges%background*cell%deth

    IF (my_tddft) THEN
       rho_atom_set => p_env%local_rho_set%rho_atom_set
       rho0_atom_set=> p_env%local_rho_set%rho0_atom_set
       rho0_mpole   => p_env%local_rho_set%rho0_mpole
       ecoul_1c     => p_env%hartree_local%ecoul_1c
    ELSE
       CALL get_qs_env(qs_env=qs_env, &
            rho_atom_set= rho_atom_set,&
            rho0_atom_set=rho0_atom_set, &
            rho0_mpole=rho0_mpole,&
            rhoz_set=rhoz_set,&
            ecoul_1c=ecoul_1c,error=error)
    END IF

    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    IF (my_tddft) THEN
       IF(PRESENT(do_triplet)) THEN
          IF (nspins==1.AND.do_triplet) RETURN
       ELSE
          IF (nspins==1.AND.dft_control%tddfpt_control%res_etype/=tddfpt_singlet) RETURN
       ENDIF
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxg_iso_not0=max_iso)
    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,lmax_0=lmax_0)

    atenergy=.FALSE.
    IF (ASSOCIATED(qs_env%atprop)) THEN
      atenergy = qs_env%atprop%energy
      IF (atenergy) THEN
        CALL get_qs_env(qs_env=qs_env,natom=natom,error=error)
        CALL atprop_array_init(qs_env%atprop%ate1c,natom,error)
      END IF
    END IF

    !   Put to 0 the local hartree energy contribution from 1 center integrals
    energy_hartree_1c = 0.0_dp

    !   Here starts the loop over all the atoms
    DO ikind = 1,nkind

      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
            orb_basis_set=orb_basis,&
            natom=nat,grid_atom=grid_atom,&
            harmonics=harmonics,ngrid_rad=nr,&
            max_iso_not0=max_iso_not0,paw_atom=paw_atom)
      IF(paw_atom) THEN
!===========    PAW   ===============
         CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
              maxso=maxso,npgf=npgf,maxl=maxl,&
              nset=nset,zet=zet)

         max_s_harm = harmonics%max_s_harm
         llmax = harmonics%llmax

         nsotot = maxso*nset
         ALLOCATE(gsph(nr,nsotot),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(gexp(nr),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(sqrtwr(nr),g0_h_w(nr,0:lmax_0),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

         NULLIFY(Vh1_h,Vh1_s)
         ALLOCATE(Vh1_h(nr,max_iso_not0),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(Vh1_s(nr,max_iso_not0),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

         ALLOCATE(aVh1b_hh(nsotot,nsotot),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(aVh1b_ss(nsotot,nsotot),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(aVh1b_00(nsotot,nsotot),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(cg_list(2,nsoset(maxl)**2,max_s_harm),cg_n_list(max_s_harm),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

         NULLIFY(Qlm_gg,g0_h)
         CALL get_rho0_mpole(rho0_mpole=rho0_mpole, ikind=ikind, &
              l0_ikind=lmax0, &
              Qlm_gg=Qlm_gg, g0_h=g0_h)

         nchan_0 = nsoset(lmax0)

         IF(nchan_0 > max_iso_not0) CALL stop_program(routineN,moduleN,__LINE__,&
            "channels for rho0 > # max of spherical harmonics",para_env)

         NULLIFY(rrad_z,my_CG)
         my_CG  => harmonics%my_CG

         !     set to zero temporary arrays
         sqrtwr=0.0_dp
         g0_h_w=0.0_dp
         gexp=0.0_dp
         gsph=0.0_dp

         sqrtwr(1:nr) = SQRT(grid_atom%wr(1:nr))
         DO l_ang = 0,lmax0
            g0_h_w(1:nr,l_ang) = g0_h(1:nr,l_ang)*grid_atom%wr(1:nr)
         END DO

         m1 = 0
         DO iset1 = 1,nset
            n1 = nsoset(lmax(iset1))
            DO ipgf1  = 1,npgf(iset1)
               gexp(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))*sqrtwr(1:nr)
               DO is1 = nsoset(lmin(iset1)-1)+1,nsoset(lmax(iset1))
                  iso = is1 + (ipgf1-1)*n1 + m1
                  l_ang = indso(1,is1)
                  gsph(1:nr,iso) = grid_atom%rad2l(1:nr,l_ang)*gexp(1:nr)
               END DO  ! is1
            END DO  ! ipgf1
            m1 = m1 + maxso
         END DO  ! iset1

         !     Distribute the atoms of this kind
         num_pe = para_env%num_pe
         mepos  = para_env%mepos
         bo = get_limit( nat, num_pe, mepos )

         DO iat = bo(1), bo(2) !1,nat
            iatom = atom_list(iat)
            rho_atom => rho_atom_set(iatom)

            NULLIFY(rrad_z,vrrad_z,rrad_0,vrrad_0)
            IF(core_charge) THEN
               rrad_z  => rhoz_set(ikind)%r_coef
               vrrad_z => rhoz_set(ikind)%vr_coef
            END IF
            rrad_0 => rho0_atom_set(iatom)%rho0_rad_h%r_coef
            vrrad_0 => rho0_atom_set(iatom)%vrho0_rad_h%r_coef
            IF(my_periodic .AND. back_ch .GT. 1.E-3_dp) THEN
               factor = -2.0_dp*pi/3.0_dp*SQRT(fourpi)*qs_charges%background
            ELSE
               factor = 0._dp
            END IF

            CALL Vh_1c_atom_potential(rho_atom,vrrad_0,&
                 grid_atom,iatom,core_charge,vrrad_z,Vh1_h,Vh1_s,&
                 nchan_0,nspins,max_iso_not0,factor,error)

            CALL Vh_1c_atom_energy(energy_hartree_1c,ecoul_1c,rho_atom,rrad_0,&
                 grid_atom,iatom,core_charge,rrad_z,Vh1_h,Vh1_s,&
                 nchan_0,nspins,max_iso_not0,atenergy,qs_env%atprop%ate1c,error)

            CALL Vh_1c_atom_integrals(rho_atom,&
                 aVh1b_hh,aVh1b_ss,aVh1b_00,Vh1_h,Vh1_s,max_iso_not0,&
                 max_s_harm,llmax,cg_list,cg_n_list,&
                 nset,npgf,lmin,lmax,nsotot,maxso,nspins,nchan_0,gsph,g0_h_w,my_CG,Qlm_gg,error)

         END DO  ! iat

         DEALLOCATE(aVh1b_hh,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(aVh1b_ss,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(aVh1b_00,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(Vh1_h,Vh1_s, STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(cg_list,cg_n_list,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(gsph,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(gexp,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(sqrtwr,g0_h_w,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      ELSE
!===========   NO  PAW   ===============
!  This term is taken care of using the core density as in GPW
           CYCLE
      END IF  ! paw
    END DO  ! ikind

    CALL mp_sum(energy_hartree_1c,para_env%group)

    CALL timestop(handle)

  END SUBROUTINE Vh_1c_gg_integrals

! *****************************************************************************

  SUBROUTINE Vh_1c_atom_potential(rho_atom,vrrad_0,&
             grid_atom,iatom,core_charge,vrrad_z,Vh1_h,Vh1_s,&
             nchan_0,nspins,max_iso_not0,bfactor,error)

    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(:, :), POINTER       :: vrrad_0
    TYPE(grid_atom_type), POINTER            :: grid_atom
    INTEGER, INTENT(IN)                      :: iatom
    LOGICAL, INTENT(IN)                      :: core_charge
    REAL(dp), DIMENSION(:), POINTER          :: vrrad_z
    REAL(dp), DIMENSION(:, :), POINTER       :: Vh1_h, Vh1_s
    INTEGER, INTENT(IN)                      :: nchan_0, nspins, max_iso_not0
    REAL(dp), INTENT(IN)                     :: bfactor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Vh_1c_atom_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ir, iso, ispin, nr
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: vr_h, vr_s

    nr = grid_atom%nr

    NULLIFY(vr_h,vr_s)
    CALL get_rho_atom(rho_atom=rho_atom,vrho_rad_h=vr_h,vrho_rad_s=vr_s)

    Vh1_h = 0.0_dp
    Vh1_s = 0.0_dp

    IF (core_charge) Vh1_h(:,1) = vrrad_z(:)

    DO iso = 1,nchan_0
       Vh1_s(:,iso) = vrrad_0(:,iso)
    END DO

    DO ispin = 1,nspins
       DO iso =1,max_iso_not0
          Vh1_h(:,iso) = Vh1_h(:,iso) + vr_h(ispin)%r_coef(:,iso)
          Vh1_s(:,iso) = Vh1_s(:,iso) + vr_s(ispin)%r_coef(:,iso)
       END DO
    END DO

    IF(bfactor /= 0._dp) THEN
       DO ir = 1,nr
          Vh1_h(ir,1) =  Vh1_h(ir,1)  + bfactor * grid_atom%rad2(ir)*grid_atom%wr(ir)
          Vh1_s(ir,1) =  Vh1_s(ir,1)  + bfactor * grid_atom%rad2(ir)*grid_atom%wr(ir)
       END DO
    END IF

  END SUBROUTINE Vh_1c_atom_potential

! *****************************************************************************

  SUBROUTINE Vh_1c_atom_energy(energy_hartree_1c,ecoul_1c,rho_atom,rrad_0,&
             grid_atom,iatom,core_charge,rrad_z,Vh1_h,Vh1_s,&
             nchan_0,nspins,max_iso_not0,atenergy,ate1c,error)

    REAL(dp), INTENT(INOUT)                  :: energy_hartree_1c
    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(:, :), POINTER       :: rrad_0
    TYPE(grid_atom_type), POINTER            :: grid_atom
    INTEGER, INTENT(IN)                      :: iatom
    LOGICAL, INTENT(IN)                      :: core_charge
    REAL(dp), DIMENSION(:), POINTER          :: rrad_z
    REAL(dp), DIMENSION(:, :), POINTER       :: Vh1_h, Vh1_s
    INTEGER, INTENT(IN)                      :: nchan_0, nspins, max_iso_not0
    LOGICAL, INTENT(IN)                      :: atenergy
    REAL(dp), DIMENSION(:), POINTER          :: ate1c
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Vh_1c_atom_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iso, ispin, nr
    REAL(dp)                                 :: ecoul_1_0, ecoul_1_h, &
                                                ecoul_1_s, ecoul_1_z
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: r_h, r_s

    nr = grid_atom%nr

    NULLIFY(r_h,r_s)
    CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,rho_rad_s=r_s)

    !       Calculate the contributions to Ecoul coming from Vh1_h*rhoz
    ecoul_1_z = 0.0_dp
    IF(core_charge) THEN
       ecoul_1_z = 0.5_dp*SUM(Vh1_h(:,1)*rrad_z(:)*grid_atom%wr(:))
    END IF

    !       Calculate the contributions to Ecoul coming from  Vh1_s*rho0
    ecoul_1_0 = 0.0_dp
    DO iso = 1,nchan_0
       ecoul_1_0 = ecoul_1_0 + 0.5_dp*SUM(Vh1_s(:,iso)*rrad_0(:,iso)*grid_atom%wr(:))
    END DO

    !       Calculate the contributions to Ecoul coming from Vh1_h*rho1_h and Vh1_s*rho1_s
    ecoul_1_s = 0.0_dp
    ecoul_1_h = 0.0_dp
    DO ispin = 1,nspins
       DO iso = 1,max_iso_not0
          ecoul_1_s = ecoul_1_s + 0.5_dp*SUM(Vh1_s(:,iso)*r_s(ispin)%r_coef(:,iso)*grid_atom%wr(:))
          ecoul_1_h = ecoul_1_h + 0.5_dp*SUM(Vh1_h(:,iso)*r_h(ispin)%r_coef(:,iso)*grid_atom%wr(:))
       END DO 
    END DO 

    CALL set_ecoul_1c(ecoul_1c,iatom,ecoul_1_z=ecoul_1_z,ecoul_1_0=ecoul_1_0)
    CALL set_ecoul_1c(ecoul_1c=ecoul_1c,iatom=iatom,ecoul_1_h=ecoul_1_h,ecoul_1_s=ecoul_1_s)

    energy_hartree_1c = energy_hartree_1c + ecoul_1_z - ecoul_1_0
    energy_hartree_1c = energy_hartree_1c + ecoul_1_h - ecoul_1_s

    ! atomic energy contribution
    IF (atenergy) THEN
       ate1c(iatom) = ate1c(iatom) + ecoul_1_z - ecoul_1_0
    END IF

  END SUBROUTINE Vh_1c_atom_energy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE Vh_1c_atom_integrals(rho_atom,&
             aVh1b_hh,aVh1b_ss,aVh1b_00,Vh1_h,Vh1_s,max_iso_not0,&
             max_s_harm,llmax,cg_list,cg_n_list,&
             nset,npgf,lmin,lmax,nsotot,maxso,nspins,nchan_0,gsph,g0_h_w,my_CG,Qlm_gg,error)

    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(:, :)                :: aVh1b_hh, aVh1b_ss, aVh1b_00
    REAL(dp), DIMENSION(:, :), POINTER       :: Vh1_h, Vh1_s
    INTEGER, INTENT(IN)                      :: max_iso_not0, max_s_harm, &
                                                llmax
    INTEGER, DIMENSION(:, :, :)              :: cg_list
    INTEGER, DIMENSION(:)                    :: cg_n_list
    INTEGER, INTENT(IN)                      :: nset
    INTEGER, DIMENSION(:), POINTER           :: npgf, lmin, lmax
    INTEGER, INTENT(IN)                      :: nsotot, maxso, nspins, nchan_0
    REAL(dp), DIMENSION(:, :), POINTER       :: gsph
    REAL(dp), DIMENSION(:, 0:)               :: g0_h_w
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG, Qlm_gg
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Vh_1c_atom_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER :: icg, ipgf1, ipgf2, ir, is1, is2, iset1, iset2, iso, iso1, &
      iso2, ispin, l_ang, m1, m2, max_iso_not0_local, n1, n2, nr
    REAL(dp)                                 :: gVg_0, gVg_h, gVg_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h, int_local_s

    NULLIFY(int_local_h,int_local_s)
    CALL get_rho_atom(rho_atom=rho_atom,&
         ga_Vlocal_gb_h=int_local_h,&
         ga_Vlocal_gb_s=int_local_s)

    !       Calculate the integrals of the potential with 2 primitives
    aVh1b_hh = 0.0_dp
    aVh1b_ss = 0.0_dp
    aVh1b_00 = 0.0_dp

    nr = SIZE(gsph,1)

    m1 = 0
    DO iset1 = 1,nset
       m2 = 0
       DO iset2 = 1,nset
          CALL get_none0_cg_list(my_CG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                 max_s_harm,llmax,cg_list,cg_n_list,max_iso_not0_local,error)

          n1 = nsoset(lmax(iset1))
          DO ipgf1  = 1,npgf(iset1)
             n2 = nsoset(lmax(iset2))
             DO ipgf2  = 1,npgf(iset2)
                !               with contributions to  V1_s*rho0
                DO iso = 1,nchan_0
                   l_ang = indso(1,iso)
                   gVg_0 = SUM(Vh1_s(:,iso)*g0_h_w(:,l_ang))
                   DO icg = 1,cg_n_list(iso)
                      is1 = cg_list(1,icg,iso)
                      is2 = cg_list(2,icg,iso)

                      iso1 = is1 + n1*(ipgf1-1) + m1
                      iso2 = is2 + n2*(ipgf2-1) + m2
                      gVg_h = 0.0_dp
                      gVg_s = 0.0_dp

                      DO ir = 1,nr
                         gVg_h = gVg_h + gsph(ir,iso1)*gsph(ir,iso2)* Vh1_h(ir,iso)
                         gVg_s = gVg_s + gsph(ir,iso1)*gsph(ir,iso2)* Vh1_s(ir,iso)
                      END DO  ! ir

                      aVh1b_hh(iso1,iso2) = aVh1b_hh(iso1,iso2) + gVg_h*my_CG(is1,is2,iso)
                      aVh1b_ss(iso1,iso2) = aVh1b_ss(iso1,iso2) + gVg_s*my_CG(is1,is2,iso)
                      aVh1b_00(iso1,iso2) = aVh1b_00(iso1,iso2) + gVg_0*Qlm_gg(iso1,iso2,iso)

                   END DO  !icg
                END DO  ! iso
                !               without contributions to  V1_s*rho0
                DO iso = nchan_0+1,max_iso_not0
                   DO icg = 1,cg_n_list(iso)
                      is1 = cg_list(1,icg,iso)
                      is2 = cg_list(2,icg,iso)

                      iso1 = is1 + n1*(ipgf1-1) + m1
                      iso2 = is2 + n2*(ipgf2-1) + m2
                      gVg_h = 0.0_dp
                      gVg_s = 0.0_dp

                      DO ir = 1,nr
                         gVg_h = gVg_h + gsph(ir,iso1)*gsph(ir,iso2)* Vh1_h(ir,iso)
                         gVg_s = gVg_s + gsph(ir,iso1)*gsph(ir,iso2)* Vh1_s(ir,iso)
                      END DO  ! ir

                      aVh1b_hh(iso1,iso2) = aVh1b_hh(iso1,iso2) + gVg_h*my_CG(is1,is2,iso)
                      aVh1b_ss(iso1,iso2) = aVh1b_ss(iso1,iso2) + gVg_s*my_CG(is1,is2,iso)

                   END DO  !icg
                END DO  ! iso
             END DO  ! ipgf2
          END DO  ! ipgf1
          m2 = m2 + maxso
       END DO  ! iset2
       m1 = m1 + maxso
    END DO  !iset1

    DO ispin = 1,nspins
       CALL daxpy(nsotot*nsotot,1.0_dp,aVh1b_hh,1,int_local_h(ispin)%r_coef,1)
       CALL daxpy(nsotot*nsotot,1.0_dp,aVh1b_ss,1,int_local_s(ispin)%r_coef,1)
       CALL daxpy(nsotot*nsotot,-1.0_dp,aVh1b_00,1,int_local_h(ispin)%r_coef,1)
       CALL daxpy(nsotot*nsotot,-1.0_dp,aVh1b_00,1,int_local_s(ispin)%r_coef,1)
    END DO  ! ispin

  END SUBROUTINE Vh_1c_atom_integrals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE hartree_local_methods

