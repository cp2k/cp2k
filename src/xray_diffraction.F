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
!!     M. Krack, A. Gambirasio, and M.\ Parrinello,
!!     "Ab-initio x-ray scattering of liquid water",
!!     J. Chem. Phys. 117, 9409 (2002)
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
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE coefficient_types,               ONLY: coeff_scale,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_gather,&
                                             mp_gatherv,&
                                             mp_sum
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
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters in this module

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xray_diffraction'

  ! Public subroutines

  PUBLIC :: xray_diffraction_spectrum

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE xray_diffraction_spectrum(qs_env,unit_number,q_max,error)

    ! Calculate the coherent X-ray diffraction spectrum.
    ! The total electronic density in reciprocal space (g-space) is calculated.

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: unit_number
    REAL(KIND=dp), INTENT(IN)                :: q_max
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xray_diffraction_spectrum', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 10

    INTEGER :: atom, group, handle, i, iatom, ico, ico1_pgf, ico1_set, ig, &
      ig_shell, ikind, ipe, ipgf, iset, ishell, iso, iso1_pgf, iso1_set, &
      ispin, istat, jco, jco1_pgf, jco1_set, jg, jpgf, jset, jso, jso1_pgf, &
      jso1_set, la, lb, maxco, maxso, mype, na, natom, nb, ncoa, ncob, &
      ncotot, ng, nkind, npe, nset, nshell, nshell_gather, nsoa, nsob, &
      nsotot, nspin, source
    INTEGER, DIMENSION(:), POINTER           :: atom_list, index, lmax, lmin, &
                                                ng_shell, ng_shell_gather, &
                                                npgf, nshell_pe, offset_pe
    LOGICAL                                  :: failure, paw_atom
    REAL(KIND=dp)                            :: eps_rho_gspace, f, f2, q, &
                                                rho_hard, rho_soft, &
                                                rho_total, scale
    REAL(KIND=dp), DIMENSION(3)              :: dg, ra
    REAL(KIND=dp), DIMENSION(:), POINTER :: f2sum, f2sum_gather, f4sum, &
      f4sum_gather, fmax, fmax_gather, fmin, fmin_gather, fsum, fsum_gather, &
      q_shell, q_shell_gather
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: delta_cpc, pab, work, zet
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho_elec_gspace, &
                                                rhotot_elec_gspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: cpc_h, cpc_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom

! Local parameter
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
               f2sum,&
               f2sum_gather,&
               f4sum,&
               f4sum_gather,&
               fmax,&
               fmax_gather,&
               fmin,&
               fmin_gather,&
               fsum,&
               fsum_gather,&
               index,&
               lmax,&
               lmin)

      NULLIFY (ng_shell,&
               ng_shell_gather,&
               npgf,&
               nshell_pe,&
               offset_pe,&
               orb_basis_set,&
               pab,&
               particle_set,&
               q_shell,&
               q_shell_gather,&
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
      mype = para_env%mepos
      npe = para_env%num_pe
      source = para_env%source
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

      ! Each process has to loop over all PAW atoms, since the g-space grid
      ! is already distributed over all processes

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
                    gsq_max=q_max*q_max,&
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

      IF (pw_grid%have_g0) THEN
        rho_total = -rhotot_elec_gspace%pw%cc(1)
      ELSE
        rho_total = 0.0_dp
      END IF
      CALL mp_sum(rho_total,group)

      rho_hard = rho_total - rho_soft

      ! Calculate the coherent x-ray spectrum

      ! Now we have to gather the data from all processes, since each
      ! process has only worked his sub-grid

      dg(:) = twopi/(pw_grid%npts(:)*pw_grid%dr(:))

      ! Scan the g-vector shells

      CALL reallocate(q_shell,1,nblock)
      CALL reallocate(ng_shell,1,nblock)

      ng = SIZE(pw_grid%gsq)

      jg = 1
      nshell = 1
      q_shell(1) = SQRT(pw_grid%gsq(1))
      ng_shell(1) = 1

      DO ig=2,ng
        IF (ABS(pw_grid%gsq(ig) - pw_grid%gsq(jg)) > 1.0E-12_dp) THEN
          nshell = nshell + 1
          IF (nshell > SIZE(q_shell)) THEN
            CALL reallocate(q_shell,1,SIZE(q_shell)+nblock)
            CALL reallocate(ng_shell,1,SIZE(ng_shell)+nblock)
          END IF
          q = SQRT(pw_grid%gsq(ig))
          IF (q > q_max) THEN
            nshell = nshell - 1
            EXIT
          END IF
          q_shell(nshell) = q
          ng_shell(nshell) = 1
          jg = ig
        ELSE
          ng_shell(nshell) = ng_shell(nshell) + 1
        END IF
      END DO

      CALL reallocate(q_shell,1,nshell)
      CALL reallocate(ng_shell,1,nshell)
      CALL reallocate(fmin,1,nshell)
      CALL reallocate(fmax,1,nshell)
      CALL reallocate(fsum,1,nshell)
      CALL reallocate(f2sum,1,nshell)
      CALL reallocate(f4sum,1,nshell)

      ig = 0
      DO ishell=1,nshell
        fmin(ishell) = HUGE(0.0_dp)
        fmax(ishell) = 0.0_dp
        fsum(ishell) = 0.0_dp
        f2sum(ishell) = 0.0_dp
        f4sum(ishell) = 0.0_dp
        DO ig_shell=1,ng_shell(ishell)
          f = ABS(rhotot_elec_gspace%pw%cc(ig+ig_shell))
          fmin(ishell) = MIN(fmin(ishell),f)
          fmax(ishell) = MAX(fmax(ishell),f)
          fsum(ishell) = fsum(ishell) + f
          f2 = f*f
          f2sum(ishell) = f2sum(ishell) + f2
          f4sum(ishell) = f4sum(ishell) + f2*f2
        END DO
        ig = ig + ng_shell(ishell)
      END DO

      CALL reallocate(nshell_pe,0,npe-1)
      CALL reallocate(offset_pe,0,npe-1)

      ! Root (source) process gathers the number of shell of each process

      CALL mp_gather(nshell,nshell_pe,source,group)

      ! Only the root process which has to print the full spectrum has to
      ! allocate here the receive buffers with their real sizes

      IF (unit_number > 0) THEN
        nshell_gather = SUM(nshell_pe)
        offset_pe(0) = 0
        DO ipe=1,npe-1
          offset_pe(ipe) = offset_pe(ipe-1) + nshell_pe(ipe-1)
        END DO
      ELSE
        nshell_gather = 1 ! dummy value for the non-root processes
      END IF

      CALL reallocate(q_shell_gather,1,nshell_gather)
      CALL reallocate(ng_shell_gather,1,nshell_gather)
      CALL reallocate(fmin_gather,1,nshell_gather)
      CALL reallocate(fmax_gather,1,nshell_gather)
      CALL reallocate(fsum_gather,1,nshell_gather)
      CALL reallocate(f2sum_gather,1,nshell_gather)
      CALL reallocate(f4sum_gather,1,nshell_gather)

      CALL mp_gatherv(q_shell,q_shell_gather,nshell_pe,offset_pe,source,group)
      CALL mp_gatherv(ng_shell,ng_shell_gather,nshell_pe,offset_pe,source,group)
      CALL mp_gatherv(fmax,fmax_gather,nshell_pe,offset_pe,source,group)
      CALL mp_gatherv(fmin,fmin_gather,nshell_pe,offset_pe,source,group)
      CALL mp_gatherv(fsum,fsum_gather,nshell_pe,offset_pe,source,group)
      CALL mp_gatherv(f2sum,f2sum_gather,nshell_pe,offset_pe,source,group)
      CALL mp_gatherv(f4sum,f4sum_gather,nshell_pe,offset_pe,source,group)

      IF (ASSOCIATED(offset_pe)) THEN
        DEALLOCATE (offset_pe,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(nshell_pe)) THEN
        DEALLOCATE (nshell_pe,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      ! Print X-ray diffraction spectrum (I/O node only)

      IF (unit_number > 0) THEN

        CALL reallocate(index,1,nshell_gather)

        ! Sort the gathered shells

        CALL sort(q_shell_gather,nshell_gather,index)

        ! Allocate final arrays of sufficient size, i.e. nshell_gather
        ! is always greater or equal the final nshell value

        CALL reallocate(q_shell,1,nshell_gather)
        CALL reallocate(ng_shell,1,nshell_gather)
        CALL reallocate(fmin,1,nshell_gather)
        CALL reallocate(fmax,1,nshell_gather)
        CALL reallocate(fsum,1,nshell_gather)
        CALL reallocate(f2sum,1,nshell_gather)
        CALL reallocate(f4sum,1,nshell_gather)

        jg = 1
        nshell = 1
        q_shell(1) = q_shell_gather(1)
        i = INDEX(1)
        ng_shell(1) = ng_shell_gather(i)
        fmin(1) = fmin_gather(i)
        fmax(1) = fmax_gather(i)
        fsum(1) = fsum_gather(i)
        f2sum(1) = f2sum_gather(i)
        f4sum(1) = f4sum_gather(i)

        DO ig=2,nshell_gather
          i = INDEX(ig)
          IF (ABS(q_shell_gather(ig) - q_shell_gather(jg)) > 1.0E-12_dp) THEN
            nshell = nshell + 1
            q_shell(nshell) = q_shell_gather(ig)
            ng_shell(nshell) = ng_shell_gather(i)
            fmin(nshell) = fmin_gather(i)
            fmax(nshell) = fmax_gather(i)
            fsum(nshell) = fsum_gather(i)
            f2sum(nshell) = f2sum_gather(i)
            f4sum(nshell) = f4sum_gather(i)
            jg = ig
          ELSE
            ng_shell(nshell) = ng_shell(nshell) + ng_shell_gather(i)
            fmin(nshell) = MIN(fmin(nshell),fmin_gather(i))
            fmax(nshell) = MAX(fmax(nshell),fmax_gather(i))
            fsum(nshell) = fsum(nshell) + fsum_gather(i)
            f2sum(nshell) = f2sum(nshell) + f2sum_gather(i)
            f4sum(nshell) = f4sum(nshell) + f4sum_gather(i)
          END IF
        END DO

        ! The auxiliary index array is no longer needed now

        IF (ASSOCIATED(index)) THEN
          DEALLOCATE (index,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        ! Allocate the final arrays for printing with their real size

        CALL reallocate(q_shell,1,nshell)
        CALL reallocate(ng_shell,1,nshell)
        CALL reallocate(fmin,1,nshell)
        CALL reallocate(fmax,1,nshell)
        CALL reallocate(fsum,1,nshell)
        CALL reallocate(f2sum,1,nshell)
        CALL reallocate(f4sum,1,nshell)

        ! Write the X-ray diffraction spectrum to the specified file

        WRITE (UNIT=unit_number,FMT="(A)")&
          "#",&
          "# Coherent X-ray diffraction spectrum",&
          "#"
        WRITE (UNIT=unit_number,FMT="(A,1X,F20.10)")&
          "# Soft electronic charge (G-space) :",rho_soft,&
          "# Hard electronic charge (G-space) :",rho_hard,&
          "# Total electronic charge (G-space):",rho_total,&
          "# Density cutoff [Rydberg]         :",2.0_dp*pw_grid%cutoff,&
          "# q(min) [1/Angstrom]              :",q_shell(2)/angstrom,&
          "# q(max) [1/Angstrom]              :",q_shell(nshell)/angstrom
        WRITE (UNIT=unit_number,FMT="(A,2X,I8)")&
          "# Number of g-vectors (grid points):",pw_grid%ngpts,&
          "# Number of g-vector shells        :",nshell
        WRITE (UNIT=unit_number,FMT="(A,3X,3I6)")&
          "# Grid size (a,b,c)                :",pw_grid%npts(1:3)
        WRITE (UNIT=unit_number,FMT="(A)")&
          "#",&
          "# shell  points         q [1/A]      <|F(q)|^2>     Min(|F(q)|) "//&
          "    Max(|F(q)|)      <|F(q)|>^2      <|F(q)|^4>"

        DO ishell=1,nshell
          WRITE (UNIT=unit_number,FMT="(T2,I6,2X,I6,5(1X,F15.6),1X,ES15.6)")&
            ishell,&
            ng_shell(ishell),&
            q_shell(ishell)/angstrom,&
            f2sum(ishell)/REAL(ng_shell(ishell),KIND=dp),&
            fmin(ishell),&
            fmax(ishell),&
            (fsum(ishell)/REAL(ng_shell(ishell),KIND=dp))**2,&
            f4sum(ishell)/REAL(ng_shell(ishell),KIND=dp)
        END DO

      END IF

      ! Release work storage

      IF (ASSOCIATED(fmin)) THEN
        DEALLOCATE (fmin,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(fmax)) THEN
        DEALLOCATE (fmax,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(fsum)) THEN
        DEALLOCATE (fsum,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(f2sum)) THEN
        DEALLOCATE (f2sum,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(f4sum)) THEN
        DEALLOCATE (f4sum,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(ng_shell)) THEN
        DEALLOCATE (ng_shell,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(q_shell)) THEN
        DEALLOCATE (q_shell,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(fmin_gather)) THEN
        DEALLOCATE (fmin_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(fmax_gather)) THEN
        DEALLOCATE (fmax_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(fsum_gather)) THEN
        DEALLOCATE (fsum_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(f2sum_gather)) THEN
        DEALLOCATE (f2sum_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(f4sum_gather)) THEN
        DEALLOCATE (f4sum_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(ng_shell_gather)) THEN
        DEALLOCATE (ng_shell_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(q_shell_gather)) THEN
        DEALLOCATE (q_shell_gather,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
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
