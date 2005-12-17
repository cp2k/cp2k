!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/xray_diffraction [1.0] *
!!
!!   NAME
!!     xray_diffraction
!!
!!   FUNCTION
!!
!!   LITERATURE
!!
!!   AUTHOR
!!     Matthias Krack (30.11.2005)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE xray_diffraction

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_transform_space,&
                                             coeff_scale,&
                                             coeff_sumup,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE
  USE qs_collocate_density,            ONLY: calculate_total_rho,&
                                             collocate_pgf_product_gspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters in this module

  CHARACTER(LEN=*), PARAMETER :: moduleN = "xray_diffraction"

  ! Public subroutines

  PUBLIC :: xray_diffraction_spectrum

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE xray_diffraction_spectrum(qs_env,file_name,q_max,error)

    ! Calculate the coherent X-ray diffraction spectrum.
    ! The total electronic density in reciprocal space (g-space) is calculated.

    TYPE(qs_environment_type), POINTER           :: qs_env
    CHARACTER(LEN=*), INTENT(IN)                 :: file_name
    REAL(KIND=dp), INTENT(IN)                    :: q_max
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    ! Local parameter

    CHARACTER(LEN=*), PARAMETER :: routineN = "xray_diffraction_spectrum",&
                                   routineP = moduleN//":"//routineN
    INTEGER, PARAMETER          :: nblock = 100

    ! Local variables

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho_elec_gspace,&
                                                rhotot_elec_gspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(KIND=dp)                            :: eps_rho_gspace,f,f2,f2sum,&
      f4sum,fsum,fmax,fmin,ng_shell,q,rho_hard,rho_soft,rho_total,scale
    INTEGER                                  :: atom,group,handle,i,iatom,&
      ico,ico1_pgf,ico1_set,ig,ig_shell,ikind,ipgf,iset,ishell,iso,iso1_pgf,&
      iso1_set,ispin,istat,jco,jco1_pgf,jco1_set,jg,jpgf,jset,jso,jso1_pgf,&
      jso1_set,la,lb,maxco,maxso,na,nb,natom,ncoa,ncob,ncotot,ng,nkind,&
      nset,nshell,nshell_printed,nsoa,nsob,nsotot,nspin,unit_number
    LOGICAL                                  :: failure,ionode,paw_atom

    REAL(KIND=dp), DIMENSION(3)              :: dg,ra

    TYPE(atomic_kind_type), DIMENSION(:),&
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:),&
      POINTER                                :: particle_set
    TYPE(rho_atom_coeff), DIMENSION(:),&
      POINTER                                :: cpc_h,cpc_s
    TYPE(rho_atom_type), DIMENSION(:),&
      POINTER                                :: rho_atom_set
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: delta_cpc,pab,work,zet
    INTEGER, DIMENSION(:), POINTER           :: atom_list,idx,iwork,lmax,lmin,&
                                                npgf

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      CALL timeset(routineN,"I","",handle)

      NULLIFY (atom_list,&
               atomic_kind,&
               atomic_kind_set,&
               auxbas_pw_pool,&
               cpc_h,&
               cpc_s,&
               delta_cpc,&
               idx,&
               iwork,&
               lmax,&
               lmin,&
               npgf,&
               orb_basis_set,&
               pab,&
               particle_set,&
               rho,&
               rho_atom,&
               rho_atom_set,&
               rho_elec_gspace,&
               rhotot_elec_gspace,&
               work,&
               zet)

      CALL get_qs_env(qs_env=qs_env,&
                      atomic_kind_set=atomic_kind_set,&
                      cell=cell,&
                      dft_control=dft_control,&
                      para_env=para_env,&
                      particle_set=particle_set,&
                      pw_env=pw_env,&
                      rho=rho,&
                      rho_atom_set=rho_atom_set)

      CALL pw_env_get(pw_env=pw_env,&
                      auxbas_pw_pool=auxbas_pw_pool,&
                      error=error)

      eps_rho_gspace = dft_control%qs_control%eps_rho_gspace
      group = para_env%group
      ionode = (para_env%mepos == para_env%source)
      nkind = SIZE(atomic_kind_set)
      nspin = dft_control%nspins

      ! Plane waves grid to assemble the total electronic density

      ALLOCATE (rhotot_elec_gspace,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL pw_pool_init_coeff(pool=auxbas_pw_pool,&
                              coeff=rhotot_elec_gspace,&
                              use_data=COMPLEXDATA1D,&
                              in_space=RECIPROCALSPACE,&
                              error=error)
      CALL coeff_zero(rhotot_elec_gspace)

      pw_grid => rhotot_elec_gspace%pw%pw_grid

      ! Load the soft contribution of the electronic density

      ALLOCATE (rho_elec_gspace,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL pw_pool_init_coeff(pool=auxbas_pw_pool,&
                              coeff=rho_elec_gspace,&
                              use_data=COMPLEXDATA1D,&
                              in_space=RECIPROCALSPACE,&
                              error=error)

      DO ispin=1,nspin
        CALL coeff_zero(rho_elec_gspace)
        CALL coeff_transform_space(rho%rho_r(ispin),rho_elec_gspace)
        CALL coeff_sumup(rho_elec_gspace,rhotot_elec_gspace)
      END DO

      CALL pw_pool_give_back_coeff(pool=auxbas_pw_pool,&
                                   coeff=rho_elec_gspace,&
                                   error=error)
      DEALLOCATE (rho_elec_gspace,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      rho_soft = calculate_total_rho(rhotot_elec_gspace)
      CALL coeff_scale(rhotot_elec_gspace,pw_grid%vol)

      ! Add the hard contribution of the electronic density

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             atom_list=atom_list,&
                             natom=natom,&
                             orb_basis_set=orb_basis_set,&
                             paw_atom=paw_atom)

        IF (.NOT.paw_atom) CYCLE ! no PAW atom: nothing to do

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               lmax=lmax,&
                               lmin=lmin,&
                               maxco=maxco,&
                               maxso=maxso,&
                               npgf=npgf,&
                               nset=nset,&
                               zet=zet)

        ncotot = maxco*nset
        nsotot = maxso*nset
        CALL reallocate(delta_cpc,1,nsotot,1,nsotot)
        CALL reallocate(pab,1,ncotot,1,ncotot)
        CALL reallocate(work,1,maxso,1,maxco)

        DO iatom=1,natom

          atom = atom_list(iatom)
          rho_atom => rho_atom_set(atom)

          CALL get_rho_atom(rho_atom=rho_atom,&
                            cpc_h=cpc_h,&
                            cpc_s=cpc_s)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          delta_cpc = 0.0_dp

          DO ispin=1,nspin
            delta_cpc = delta_cpc + cpc_h(ispin)%r_coef - cpc_s(ispin)%r_coef
          END DO

          scale = 1.0_dp

          DO iset=1,nset
            ico1_set = (iset - 1)*maxco + 1
            iso1_set = (iset - 1)*maxso + 1
            ncoa = ncoset(lmax(iset))
            nsoa = nsoset(lmax(iset))
            DO jset=1,nset
              jco1_set = (jset - 1)*maxco + 1
              jso1_set = (jset - 1)*maxso + 1
              ncob = ncoset(lmax(jset))
              nsob = nsoset(lmax(jset))
              DO ipgf=1,npgf(iset)
                ico1_pgf = ico1_set + (ipgf - 1)*ncoa
                iso1_pgf = iso1_set + (ipgf - 1)*nsoa
                DO jpgf=1,npgf(jset)
                  jco1_pgf = jco1_set + (jpgf - 1)*ncob
                  jso1_pgf = jso1_set + (jpgf - 1)*nsob
                  ico = ico1_pgf + ncoset(lmin(iset)-1)
                  iso = iso1_pgf + nsoset(lmin(iset)-1)

                  ! Transformation spherical to Cartesian

                  DO la=lmin(iset),lmax(iset)
                    jco = jco1_pgf + ncoset(lmin(jset)-1)
                    jso = jso1_pgf + nsoset(lmin(jset)-1)
                    DO lb=lmin(jset),lmax(jset)
                      CALL dgemm("N","N",nso(la),nco(lb),nso(lb),1.0_dp,&
                                 delta_cpc(iso,jso),SIZE(delta_cpc,1),&
                                 orbtramat(lb)%slm,nso(lb),0.0_dp,work,&
                                 maxso)
                      CALL dgemm("T","N",nco(la),nco(lb),nso(la),1.0_dp,&
                                 orbtramat(la)%slm,nso(la),work,maxso,&
                                 0.0_dp,pab(ico,jco),SIZE(pab,1))
                      jco = jco + nco(lb)
                      jso = jso + nso(lb)
                    END DO ! next lb
                    ico = ico + nco(la)
                    iso = iso + nso(la)
                  END DO ! next la

                  ! Collocate current product of primitive Cartesian functions

                  na = ico1_pgf - 1
                  nb = jco1_pgf - 1

                  CALL collocate_pgf_product_gspace(&
                    la_max=lmax(iset),&
                    zeta=zet(ipgf,iset),&
                    la_min=lmin(iset),&
                    lb_max=lmax(jset),&
                    zetb=zet(jpgf,jset),&
                    lb_min=lmin(jset),&
                    ra=ra,&
                    rab=(/0.0_dp,0.0_dp,0.0_dp/),&
                    rab2=0.0_dp,&
                    scale=scale,&
                    pab=pab,&
                    na=na,&
                    nb=nb,&
                    eps_rho_gspace=eps_rho_gspace,&
                    pw=rhotot_elec_gspace%pw)

                END DO ! next primitive Gaussian function "jpgf"
              END DO ! next primitive Gaussian function "ipgf"
            END DO ! next shell set "jset"
          END DO ! next shell set "iset"
        END DO ! next atom "iatom" of atomic kind "ikind"
      END DO ! next atomic kind "ikind"

      ! Release work storage

      IF (ASSOCIATED(delta_cpc)) THEN
        DEALLOCATE (delta_cpc,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(work)) THEN
        DEALLOCATE (work,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(pab)) THEN
        DEALLOCATE (pab,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      CALL mp_sum(rhotot_elec_gspace%pw%cc,group)

      rho_total = -rhotot_elec_gspace%pw%cc(1)
      rho_hard = rho_total - rho_soft

      ! Calculate the coherent x-ray spectrum

      IF (ionode) THEN

        CALL open_file(file_name=TRIM(file_name),&
                       file_action="WRITE",&
                       file_form="FORMATTED",&
                       file_status="REPLACE",&
                       unit_number=unit_number)

        ng = pw_grid%ngpts
        dg(:) = twopi/(pw_grid%npts(:)*pw_grid%dr(:))

        ! Scan the g-vector shells

        CALL reallocate(idx,1,nblock)
        CALL reallocate(iwork,1,nblock)

        idx(1) = 1
        iwork(1) = 1
        jg = 1
        nshell = 1
        nshell_printed = 1

        DO ig=2,ng
          IF (ABS(pw_grid%gsq(ig) - pw_grid%gsq(jg)) > 1.0E-12_dp) THEN
            nshell = nshell + 1
            q = SQRT(pw_grid%gsq(ig))
            IF (q <= q_max) nshell_printed = nshell
            IF (nshell > SIZE(idx)) THEN
              CALL reallocate(idx,1,SIZE(idx)+nblock)
              CALL reallocate(iwork,1,SIZE(idx)+nblock)
            END IF
            idx(nshell) = ig
            iwork(nshell) = 1
            jg = ig
          ELSE
            iwork(nshell) = iwork(nshell) + 1
          END IF
        END DO

        CALL reallocate(idx,1,nshell)
        CALL reallocate(iwork,1,nshell)

        ! Print header of the X-ray diffraction spectrum

        WRITE (UNIT=unit_number,FMT="(A)")&
          "#",&
          "# Coherent X-ray diffraction spectrum",&
          "#"
        WRITE (UNIT=unit_number,FMT="(A,1X,F20.10)")&
          "# Soft electronic charge (G-space) :",rho_soft,&
          "# Hard electronic charge (G-space) :",rho_hard,&
          "# Total electronic charge (G-space):",rho_total,&
          "# Density cutoff [Rydberg]         :",2.0_dp*pw_grid%cutoff,&
          "# q(min) [1/Angstrom]              :",SQRT(pw_grid%gsq(2))/angstrom,&
          "# q(max) [1/Angstrom]              :",SQRT(pw_grid%gsq(ng))/angstrom,&
          "# q(max) [1/Angstrom] (printed)    :",q_max/angstrom
        WRITE (UNIT=unit_number,FMT="(A,2X,I8)")&
          "# Number of g-vectors (grid points):",ng,&
          "# Number of g-vector shells        :",nshell,&
          "# Number of printed g-vector shells:",nshell_printed
        WRITE (UNIT=unit_number,FMT="(A)")&
          "#",&
          "# shell  points         q [1/A]      <|F(q)|^2>     Min(|F(q)|) "//&
          "    Max(|F(q)|)      <|F(q)|>^2      <|F(q)|^4>"

        ! Print X-ray diffraction spectrum

        DO ishell=1,nshell
          q = SQRT(pw_grid%gsq(idx(ishell)))
          IF (q > q_max) EXIT
          ng_shell = REAL(iwork(ishell),dp)
          fmin = HUGE(0.0_dp)
          fmax = 0.0_dp
          fsum = 0.0_dp
          f2sum = 0.0_dp
          f4sum = 0.0_dp
          DO ig_shell=0,iwork(ishell)-1
            i = idx(ishell) + ig_shell
            f = ABS(rhotot_elec_gspace%pw%cc(i))
            fmin = MIN(fmin,f)
            fmax = MAX(fmax,f)
            fsum = fsum + f
            f2 = f*f
            f2sum = f2sum + f2
            f4sum = f4sum + f2*f2
          END DO
          WRITE (UNIT=unit_number,FMT="(T2,I6,2X,I6,5(1X,F15.6),1X,E15.6)")&
            ishell,iwork(ishell),q/angstrom,f2sum/ng_shell,fmin,fmax,&
            (fsum/ng_shell)**2,f4sum/ng_shell
        END DO

        ! Release work storage

        IF (ASSOCIATED(iwork)) THEN
          DEALLOCATE (iwork,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(idx)) THEN
          DEALLOCATE (idx,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        CALL close_file(unit_number)

      END IF

      CALL pw_pool_give_back_coeff(pool=auxbas_pw_pool,&
                                   coeff=rhotot_elec_gspace,&
                                   error=error)
      DEALLOCATE (rhotot_elec_gspace,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL timestop(0.0_dp,handle)

    END IF

  END SUBROUTINE xray_diffraction_spectrum

  ! ***************************************************************************

END MODULE xray_diffraction
