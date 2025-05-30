!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par Literature
!>      M. Krack, A. Gambirasio, and M. Parrinello,
!>      "Ab-initio x-ray scattering of liquid water",
!>      J. Chem. Phys. 117, 9409 (2002)
!> \author Matthias Krack
!> \date   30.11.2005
! **************************************************************************************************
MODULE xray_diffraction
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Krack2002,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE physcon,                         ONLY: angstrom
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: get_pw_grid_info
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xray_diffraction'

   PUBLIC :: calculate_rhotot_elec_gspace, &
             xray_diffraction_spectrum

CONTAINS

! **************************************************************************************************
!> \brief Calculate the coherent X-ray diffraction spectrum using the total
!>        electronic density in reciprocal space (g-space).
!> \param qs_env ...
!> \param unit_number ...
!> \param q_max ...
!> \date   30.11.2005
!> \author Matthias Krack
! **************************************************************************************************
   SUBROUTINE xray_diffraction_spectrum(qs_env, unit_number, q_max)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_number
      REAL(KIND=dp), INTENT(IN)                          :: q_max

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xray_diffraction_spectrum'
      INTEGER, PARAMETER                                 :: nblock = 100

      INTEGER                                            :: handle, i, ig, ig_shell, ipe, ishell, &
                                                            jg, ng, npe, nshell, nshell_gather
      INTEGER(KIND=int_8)                                :: ngpts
      INTEGER, DIMENSION(3)                              :: npts
      INTEGER, DIMENSION(:), POINTER                     :: aux_index, ng_shell, ng_shell_gather, &
                                                            nshell_pe, offset_pe
      REAL(KIND=dp)                                      :: cutoff, f, f2, q, rho_hard, rho_soft, &
                                                            rho_total
      REAL(KIND=dp), DIMENSION(3)                        :: dg, dr
      REAL(KIND=dp), DIMENSION(:), POINTER :: f2sum, f2sum_gather, f4sum, f4sum_gather, fmax, &
         fmax_gather, fmin, fmin_gather, fsum, fsum_gather, gsq, q_shell, q_shell_gather
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: rhotot_elec_gspace
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set

      CPASSERT(ASSOCIATED(qs_env))

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)
      NULLIFY (aux_index)
      NULLIFY (auxbas_pw_pool)
      NULLIFY (dft_control)
      NULLIFY (f2sum)
      NULLIFY (f2sum_gather)
      NULLIFY (f4sum)
      NULLIFY (f4sum_gather)
      NULLIFY (fmax)
      NULLIFY (fmax_gather)
      NULLIFY (fmin)
      NULLIFY (fmin_gather)
      NULLIFY (fsum)
      NULLIFY (fsum_gather)
      NULLIFY (gsq)
      NULLIFY (ng_shell)
      NULLIFY (ng_shell_gather)
      NULLIFY (nshell_pe)
      NULLIFY (offset_pe)
      NULLIFY (para_env)
      NULLIFY (particle_set)
      NULLIFY (pw_env)
      NULLIFY (q_shell)
      NULLIFY (q_shell_gather)
      NULLIFY (rho)
      NULLIFY (rho_atom_set)

      CALL cite_reference(Krack2002)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      rho=rho, &
                      rho_atom_set=rho_atom_set)

      CALL pw_env_get(pw_env=pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool)

      npe = para_env%num_pe

      ! Plane waves grid to assemble the total electronic density

      CALL auxbas_pw_pool%create_pw(pw=rhotot_elec_gspace)
      CALL pw_zero(rhotot_elec_gspace)

      CALL get_pw_grid_info(pw_grid=rhotot_elec_gspace%pw_grid, &
                            dr=dr, &
                            npts=npts, &
                            cutoff=cutoff, &
                            ngpts=ngpts, &
                            gsquare=gsq)

      dg(:) = twopi/(npts(:)*dr(:))

      ! Build the total electronic density in reciprocal space

      CALL calculate_rhotot_elec_gspace(qs_env=qs_env, &
                                        auxbas_pw_pool=auxbas_pw_pool, &
                                        rhotot_elec_gspace=rhotot_elec_gspace, &
                                        q_max=q_max, &
                                        rho_hard=rho_hard, &
                                        rho_soft=rho_soft)

      rho_total = rho_hard + rho_soft

      ! Calculate the coherent X-ray spectrum

      ! Now we have to gather the data from all processes, since each
      ! process has only worked his sub-grid

      ! Scan the g-vector shells

      CALL reallocate(q_shell, 1, nblock)
      CALL reallocate(ng_shell, 1, nblock)

      ng = SIZE(gsq)

      jg = 1
      nshell = 1
      q_shell(1) = SQRT(gsq(1))
      ng_shell(1) = 1

      DO ig = 2, ng
         CPASSERT(gsq(ig) >= gsq(jg))
         IF (ABS(gsq(ig) - gsq(jg)) > 1.0E-12_dp) THEN
            nshell = nshell + 1
            IF (nshell > SIZE(q_shell)) THEN
               CALL reallocate(q_shell, 1, SIZE(q_shell) + nblock)
               CALL reallocate(ng_shell, 1, SIZE(ng_shell) + nblock)
            END IF
            q = SQRT(gsq(ig))
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

      CALL reallocate(q_shell, 1, nshell)
      CALL reallocate(ng_shell, 1, nshell)
      CALL reallocate(fmin, 1, nshell)
      CALL reallocate(fmax, 1, nshell)
      CALL reallocate(fsum, 1, nshell)
      CALL reallocate(f2sum, 1, nshell)
      CALL reallocate(f4sum, 1, nshell)

      ig = 0
      DO ishell = 1, nshell
         fmin(ishell) = HUGE(0.0_dp)
         fmax(ishell) = 0.0_dp
         fsum(ishell) = 0.0_dp
         f2sum(ishell) = 0.0_dp
         f4sum(ishell) = 0.0_dp
         DO ig_shell = 1, ng_shell(ishell)
            f = ABS(rhotot_elec_gspace%array(ig + ig_shell))
            fmin(ishell) = MIN(fmin(ishell), f)
            fmax(ishell) = MAX(fmax(ishell), f)
            fsum(ishell) = fsum(ishell) + f
            f2 = f*f
            f2sum(ishell) = f2sum(ishell) + f2
            f4sum(ishell) = f4sum(ishell) + f2*f2
         END DO
         ig = ig + ng_shell(ishell)
      END DO

      CALL reallocate(nshell_pe, 0, npe - 1)
      CALL reallocate(offset_pe, 0, npe - 1)

      ! Root (source) process gathers the number of shell of each process

      CALL para_env%gather(nshell, nshell_pe)

      ! Only the root process which has to print the full spectrum has to
      ! allocate here the receive buffers with their real sizes

      IF (unit_number > 0) THEN
         nshell_gather = SUM(nshell_pe)
         offset_pe(0) = 0
         DO ipe = 1, npe - 1
            offset_pe(ipe) = offset_pe(ipe - 1) + nshell_pe(ipe - 1)
         END DO
      ELSE
         nshell_gather = 1 ! dummy value for the non-root processes
      END IF

      CALL reallocate(q_shell_gather, 1, nshell_gather)
      CALL reallocate(ng_shell_gather, 1, nshell_gather)
      CALL reallocate(fmin_gather, 1, nshell_gather)
      CALL reallocate(fmax_gather, 1, nshell_gather)
      CALL reallocate(fsum_gather, 1, nshell_gather)
      CALL reallocate(f2sum_gather, 1, nshell_gather)
      CALL reallocate(f4sum_gather, 1, nshell_gather)

      CALL para_env%gatherv(q_shell, q_shell_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(ng_shell, ng_shell_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(fmax, fmax_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(fmin, fmin_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(fsum, fsum_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(f2sum, f2sum_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(f4sum, f4sum_gather, nshell_pe, offset_pe)

      IF (ASSOCIATED(offset_pe)) THEN
         DEALLOCATE (offset_pe)
      END IF

      IF (ASSOCIATED(nshell_pe)) THEN
         DEALLOCATE (nshell_pe)
      END IF

      ! Print X-ray diffraction spectrum (I/O node only)

      IF (unit_number > 0) THEN

         CALL reallocate(aux_index, 1, nshell_gather)

         ! Sort the gathered shells

         CALL sort(q_shell_gather, nshell_gather, aux_index)

         ! Allocate final arrays of sufficient size, i.e. nshell_gather
         ! is always greater or equal the final nshell value

         CALL reallocate(q_shell, 1, nshell_gather)
         CALL reallocate(ng_shell, 1, nshell_gather)
         CALL reallocate(fmin, 1, nshell_gather)
         CALL reallocate(fmax, 1, nshell_gather)
         CALL reallocate(fsum, 1, nshell_gather)
         CALL reallocate(f2sum, 1, nshell_gather)
         CALL reallocate(f4sum, 1, nshell_gather)

         jg = 1
         nshell = 1
         q_shell(1) = q_shell_gather(1)
         i = aux_index(1)
         ng_shell(1) = ng_shell_gather(i)
         fmin(1) = fmin_gather(i)
         fmax(1) = fmax_gather(i)
         fsum(1) = fsum_gather(i)
         f2sum(1) = f2sum_gather(i)
         f4sum(1) = f4sum_gather(i)

         DO ig = 2, nshell_gather
            i = aux_index(ig)
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
               fmin(nshell) = MIN(fmin(nshell), fmin_gather(i))
               fmax(nshell) = MAX(fmax(nshell), fmax_gather(i))
               fsum(nshell) = fsum(nshell) + fsum_gather(i)
               f2sum(nshell) = f2sum(nshell) + f2sum_gather(i)
               f4sum(nshell) = f4sum(nshell) + f4sum_gather(i)
            END IF
         END DO

         ! The auxiliary index array is no longer needed now

         IF (ASSOCIATED(aux_index)) THEN
            DEALLOCATE (aux_index)
         END IF

         ! Allocate the final arrays for printing with their real size

         CALL reallocate(q_shell, 1, nshell)
         CALL reallocate(ng_shell, 1, nshell)
         CALL reallocate(fmin, 1, nshell)
         CALL reallocate(fmax, 1, nshell)
         CALL reallocate(fsum, 1, nshell)
         CALL reallocate(f2sum, 1, nshell)
         CALL reallocate(f4sum, 1, nshell)

         ! Write the X-ray diffraction spectrum to the specified file

         WRITE (UNIT=unit_number, FMT="(A)") &
            "#", &
            "# Coherent X-ray diffraction spectrum", &
            "#"
         WRITE (UNIT=unit_number, FMT="(A,1X,F20.10)") &
            "# Soft electronic charge (G-space) :", rho_soft, &
            "# Hard electronic charge (G-space) :", rho_hard, &
            "# Total electronic charge (G-space):", rho_total, &
            "# Density cutoff [Rydberg]         :", 2.0_dp*cutoff, &
            "# q(min) [1/Angstrom]              :", q_shell(2)/angstrom, &
            "# q(max) [1/Angstrom]              :", q_shell(nshell)/angstrom, &
            "# q(max) [1/Angstrom] (requested)  :", q_max/angstrom
         WRITE (UNIT=unit_number, FMT="(A,2X,I8)") &
            "# Number of g-vectors (grid points):", ngpts, &
            "# Number of g-vector shells        :", nshell
         WRITE (UNIT=unit_number, FMT="(A,3(1X,I6))") &
            "# Grid size (a,b,c)                :", npts(1:3)
         WRITE (UNIT=unit_number, FMT="(A,3F7.3)") &
            "# dg [1/Angstrom]                  :", dg(1:3)/angstrom, &
            "# dr [Angstrom]                    :", dr(1:3)*angstrom
         WRITE (UNIT=unit_number, FMT="(A)") &
            "#", &
            "# shell  points         q [1/A]      <|F(q)|^2>     Min(|F(q)|)"// &
            "     Max(|F(q)|)      <|F(q)|>^2      <|F(q)|^4>"

         DO ishell = 1, nshell
            WRITE (UNIT=unit_number, FMT="(T2,I6,2X,I6,5(1X,F15.6),1X,ES15.6)") &
               ishell, &
               ng_shell(ishell), &
               q_shell(ishell)/angstrom, &
               f2sum(ishell)/REAL(ng_shell(ishell), KIND=dp), &
               fmin(ishell), &
               fmax(ishell), &
               (fsum(ishell)/REAL(ng_shell(ishell), KIND=dp))**2, &
               f4sum(ishell)/REAL(ng_shell(ishell), KIND=dp)
         END DO

      END IF

      ! Release work storage

      IF (ASSOCIATED(fmin)) THEN
         DEALLOCATE (fmin)
      END IF

      IF (ASSOCIATED(fmax)) THEN
         DEALLOCATE (fmax)
      END IF

      IF (ASSOCIATED(fsum)) THEN
         DEALLOCATE (fsum)
      END IF

      IF (ASSOCIATED(f2sum)) THEN
         DEALLOCATE (f2sum)
      END IF

      IF (ASSOCIATED(f4sum)) THEN
         DEALLOCATE (f4sum)
      END IF

      IF (ASSOCIATED(ng_shell)) THEN
         DEALLOCATE (ng_shell)
      END IF

      IF (ASSOCIATED(q_shell)) THEN
         DEALLOCATE (q_shell)
      END IF

      IF (ASSOCIATED(fmin_gather)) THEN
         DEALLOCATE (fmin_gather)
      END IF

      IF (ASSOCIATED(fmax_gather)) THEN
         DEALLOCATE (fmax_gather)
      END IF

      IF (ASSOCIATED(fsum_gather)) THEN
         DEALLOCATE (fsum_gather)
      END IF

      IF (ASSOCIATED(f2sum_gather)) THEN
         DEALLOCATE (f2sum_gather)
      END IF

      IF (ASSOCIATED(f4sum_gather)) THEN
         DEALLOCATE (f4sum_gather)
      END IF

      IF (ASSOCIATED(ng_shell_gather)) THEN
         DEALLOCATE (ng_shell_gather)
      END IF

      IF (ASSOCIATED(q_shell_gather)) THEN
         DEALLOCATE (q_shell_gather)
      END IF

      CALL auxbas_pw_pool%give_back_pw(rhotot_elec_gspace)

      CALL timestop(handle)

   END SUBROUTINE xray_diffraction_spectrum

! **************************************************************************************************
!> \brief  The total electronic density in reciprocal space (g-space) is
!>         calculated.
!> \param qs_env ...
!> \param auxbas_pw_pool ...
!> \param rhotot_elec_gspace ...
!> \param q_max ...
!> \param rho_hard ...
!> \param rho_soft ...
!> \param fsign ...
!> \date   14.03.2008 (splitted from the routine xray_diffraction_spectrum)
!> \author Matthias Krack
!> \note   This code assumes that the g-vectors are ordered (in gsq and %cc)
! **************************************************************************************************
   SUBROUTINE calculate_rhotot_elec_gspace(qs_env, auxbas_pw_pool, &
                                           rhotot_elec_gspace, q_max, rho_hard, &
                                           rho_soft, fsign)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rhotot_elec_gspace
      REAL(KIND=dp), INTENT(IN)                          :: q_max
      REAL(KIND=dp), INTENT(OUT)                         :: rho_hard, rho_soft
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: fsign

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_rhotot_elec_gspace'

      INTEGER :: atom, handle, iatom, ico, ico1_pgf, ico1_set, ikind, ipgf, iset, iso, iso1_pgf, &
         iso1_set, ison, ispin, jco, jco1_pgf, jco1_set, jpgf, jset, jso, jso1_pgf, jso1_set, &
         json, la, lb, maxco, maxso, na, natom, nb, ncoa, ncob, ncotot, nkind, nsatbas, nset, &
         nsoa, nsob, nsotot, nspin
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, lmax, lmin, npgf, o2nindex
      LOGICAL                                            :: orthorhombic, paw_atom
      REAL(KIND=dp)                                      :: alpha, eps_rho_gspace, rho_total, scale, &
                                                            volume
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: delta_cpc, pab, work, zet
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: rho_elec_gspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: cpc_h, cpc_s
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(auxbas_pw_pool))

      CALL timeset(routineN, handle)

      NULLIFY (atom_list)
      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (cell)
      NULLIFY (cpc_h)
      NULLIFY (cpc_s)
      NULLIFY (delta_cpc)
      NULLIFY (dft_control)
      NULLIFY (lmax)
      NULLIFY (lmin)
      NULLIFY (npgf)
      NULLIFY (basis_1c_set)
      NULLIFY (pab)
      NULLIFY (particle_set)
      NULLIFY (rho, rho_r)
      NULLIFY (rho_atom)
      NULLIFY (rho_atom_set)
      NULLIFY (work)
      NULLIFY (zet)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      rho=rho, &
                      rho_atom_set=rho_atom_set)

      CALL qs_rho_get(rho, rho_r=rho_r)
      eps_rho_gspace = dft_control%qs_control%eps_rho_gspace
      nkind = SIZE(atomic_kind_set)
      nspin = dft_control%nspins

      ! Load the soft contribution of the electronic density

      CALL auxbas_pw_pool%create_pw(pw=rho_elec_gspace)

      CALL pw_zero(rhotot_elec_gspace)

      DO ispin = 1, nspin
         CALL pw_zero(rho_elec_gspace)
         CALL pw_transfer(rho_r(ispin), rho_elec_gspace)
         IF (PRESENT(fsign) .AND. (ispin == 2)) THEN
            alpha = fsign
         ELSE
            alpha = 1.0_dp
         END IF
         CALL pw_axpy(rho_elec_gspace, rhotot_elec_gspace, alpha=alpha)
      END DO

      ! Release the auxiliary PW grid for the calculation of the soft
      ! contribution

      CALL auxbas_pw_pool%give_back_pw(rho_elec_gspace)

      rho_soft = pw_integrate_function(rhotot_elec_gspace, isign=-1)

      CALL get_pw_grid_info(pw_grid=rhotot_elec_gspace%pw_grid, &
                            vol=volume, &
                            orthorhombic=orthorhombic)
      IF (.NOT. orthorhombic) THEN
         CALL cp_abort(__LOCATION__, &
                       "The calculation of XRD spectra for non-orthorhombic cells is not implemented")
      END IF

      CALL pw_scale(rhotot_elec_gspace, volume)

      ! Add the hard contribution of the electronic density

      ! Each process has to loop over all PAW atoms, since the g-space grid
      ! is already distributed over all processes

      DO ikind = 1, nkind

         CALL get_atomic_kind(atomic_kind_set(ikind), &
                              atom_list=atom_list, &
                              natom=natom)

         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=basis_1c_set, &
                          basis_type="GAPW_1C", &
                          paw_atom=paw_atom)

         IF (.NOT. paw_atom) CYCLE ! no PAW atom: nothing to do

         CALL get_paw_basis_info(basis_1c_set, o2nindex=o2nindex, nsatbas=nsatbas)

         CALL get_gto_basis_set(gto_basis_set=basis_1c_set, &
                                lmax=lmax, &
                                lmin=lmin, &
                                maxco=maxco, &
                                maxso=maxso, &
                                npgf=npgf, &
                                nset=nset, &
                                zet=zet)

         ncotot = maxco*nset
         nsotot = maxso*nset
         CALL reallocate(delta_cpc, 1, nsatbas, 1, nsatbas)
         CALL reallocate(pab, 1, ncotot, 1, ncotot)
         CALL reallocate(work, 1, maxso, 1, maxco)

         DO iatom = 1, natom

            atom = atom_list(iatom)
            rho_atom => rho_atom_set(atom)

            CALL get_rho_atom(rho_atom=rho_atom, &
                              cpc_h=cpc_h, &
                              cpc_s=cpc_s)

            ra(:) = pbc(particle_set(iatom)%r, cell)

            delta_cpc = 0.0_dp

            DO ispin = 1, nspin
               IF (PRESENT(fsign) .AND. (ispin == 2)) THEN
                  alpha = fsign
               ELSE
                  alpha = 1.0_dp
               END IF
               delta_cpc = delta_cpc + alpha*(cpc_h(ispin)%r_coef - cpc_s(ispin)%r_coef)
            END DO

            scale = 1.0_dp

            DO iset = 1, nset
               ico1_set = (iset - 1)*maxco + 1
               iso1_set = (iset - 1)*maxso + 1
               ncoa = ncoset(lmax(iset))
               nsoa = nsoset(lmax(iset))
               DO jset = 1, nset
                  jco1_set = (jset - 1)*maxco + 1
                  jso1_set = (jset - 1)*maxso + 1
                  ncob = ncoset(lmax(jset))
                  nsob = nsoset(lmax(jset))
                  DO ipgf = 1, npgf(iset)
                     ico1_pgf = ico1_set + (ipgf - 1)*ncoa
                     iso1_pgf = iso1_set + (ipgf - 1)*nsoa
                     DO jpgf = 1, npgf(jset)
                        jco1_pgf = jco1_set + (jpgf - 1)*ncob
                        jso1_pgf = jso1_set + (jpgf - 1)*nsob
                        ico = ico1_pgf + ncoset(lmin(iset) - 1)
                        iso = iso1_pgf + nsoset(lmin(iset) - 1)

                        ! Transformation spherical to Cartesian

                        DO la = lmin(iset), lmax(iset)
                           jco = jco1_pgf + ncoset(lmin(jset) - 1)
                           jso = jso1_pgf + nsoset(lmin(jset) - 1)
                           DO lb = lmin(jset), lmax(jset)
                              ison = o2nindex(iso)
                              json = o2nindex(jso)
                              CALL dgemm("N", "N", nso(la), nco(lb), nso(lb), 1.0_dp, &
                                         delta_cpc(ison:ison + nso(la) - 1, json), SIZE(delta_cpc, 1), &
                                         orbtramat(lb)%slm, nso(lb), 0.0_dp, work, &
                                         maxso)
                              CALL dgemm("T", "N", nco(la), nco(lb), nso(la), 1.0_dp, &
                                         orbtramat(la)%slm, nso(la), work, maxso, &
                                         0.0_dp, pab(ico:ico + nco(la) - 1, jco), SIZE(pab, 1))
                              jco = jco + nco(lb)
                              jso = jso + nso(lb)
                           END DO ! next lb
                           ico = ico + nco(la)
                           iso = iso + nso(la)
                        END DO ! next la

                        ! Collocate current product of primitive Cartesian functions

                        na = ico1_pgf - 1
                        nb = jco1_pgf - 1

                        CALL collocate_pgf_product_gspace( &
                           la_max=lmax(iset), &
                           zeta=zet(ipgf, iset), &
                           la_min=lmin(iset), &
                           lb_max=lmax(jset), &
                           zetb=zet(jpgf, jset), &
                           lb_min=lmin(jset), &
                           ra=ra, &
                           rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                           rab2=0.0_dp, &
                           scale=scale, &
                           pab=pab, &
                           na=na, &
                           nb=nb, &
                           eps_rho_gspace=eps_rho_gspace, &
                           gsq_max=q_max*q_max, &
                           pw=rhotot_elec_gspace)

                     END DO ! next primitive Gaussian function "jpgf"
                  END DO ! next primitive Gaussian function "ipgf"
               END DO ! next shell set "jset"
            END DO ! next shell set "iset"
         END DO ! next atom "iatom" of atomic kind "ikind"
         DEALLOCATE (o2nindex)
      END DO ! next atomic kind "ikind"

      rho_total = pw_integrate_function(rhotot_elec_gspace, isign=-1)/volume

      rho_hard = rho_total - rho_soft

      ! Release work storage

      IF (ASSOCIATED(delta_cpc)) THEN
         DEALLOCATE (delta_cpc)
      END IF

      IF (ASSOCIATED(work)) THEN
         DEALLOCATE (work)
      END IF

      IF (ASSOCIATED(pab)) THEN
         DEALLOCATE (pab)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_rhotot_elec_gspace

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions in g-space
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rab2 ...
!> \param scale ...
!> \param pab ...
!> \param na ...
!> \param nb ...
!> \param eps_rho_gspace ...
!> \param gsq_max ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_gspace(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, rab2, scale, pab, na, nb, &
                                           eps_rho_gspace, gsq_max, pw)

      ! NOTE: this routine is much slower than the real-space version of collocate_pgf_product

      INTEGER, INTENT(IN)                                :: la_max
      REAL(dp), INTENT(IN)                               :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(dp), INTENT(IN)                               :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rab
      REAL(dp), INTENT(IN)                               :: rab2, scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: pab
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(dp), INTENT(IN)                               :: eps_rho_gspace, gsq_max
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collocate_pgf_product_gspace'

      COMPLEX(dp), DIMENSION(3)                          :: phasefactor
      COMPLEX(dp), DIMENSION(:), POINTER                 :: rag, rbg
      COMPLEX(dp), DIMENSION(:, :, :, :), POINTER        :: cubeaxis
      INTEGER                                            :: ax, ay, az, bx, by, bz, handle, i, ico, &
                                                            ig, ig2, jco, jg, kg, la, lb, &
                                                            lb_cube_min, lb_grid, ub_cube_max, &
                                                            ub_grid
      INTEGER, DIMENSION(3)                              :: lb_cube, ub_cube
      REAL(dp)                                           :: f, fa, fb, pij, prefactor, rzetp, &
                                                            twozetp, zetp
      REAL(dp), DIMENSION(3)                             :: dg, expfactor, fap, fbp, rap, rbp, rp
      REAL(dp), DIMENSION(:), POINTER                    :: g

      CALL timeset(routineN, handle)

      dg(:) = twopi/(pw%pw_grid%npts(:)*pw%pw_grid%dr(:))

      zetp = zeta + zetb
      rzetp = 1.0_dp/zetp
      f = zetb*rzetp
      rap(:) = f*rab(:)
      rbp(:) = rap(:) - rab(:)
      rp(:) = ra(:) + rap(:)
      twozetp = 2.0_dp*zetp
      fap(:) = twozetp*rap(:)
      fbp(:) = twozetp*rbp(:)

      prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)
      phasefactor(:) = EXP(CMPLX(0.0_dp, -rp(:)*dg(:), KIND=dp))
      expfactor(:) = EXP(-0.25*rzetp*dg(:)*dg(:))

      lb_cube(:) = pw%pw_grid%bounds(1, :)
      ub_cube(:) = pw%pw_grid%bounds(2, :)

      lb_cube_min = MINVAL(lb_cube(:))
      ub_cube_max = MAXVAL(ub_cube(:))

      NULLIFY (cubeaxis, g, rag, rbg)

      CALL reallocate(cubeaxis, lb_cube_min, ub_cube_max, 1, 3, 0, la_max, 0, lb_max)
      CALL reallocate(g, lb_cube_min, ub_cube_max)
      CALL reallocate(rag, lb_cube_min, ub_cube_max)
      CALL reallocate(rbg, lb_cube_min, ub_cube_max)

      lb_grid = LBOUND(pw%array, 1)
      ub_grid = UBOUND(pw%array, 1)

      DO i = 1, 3

         DO ig = lb_cube(i), ub_cube(i)
            ig2 = ig*ig
            cubeaxis(ig, i, 0, 0) = expfactor(i)**ig2*phasefactor(i)**ig
         END DO

         IF (la_max > 0) THEN
            DO ig = lb_cube(i), ub_cube(i)
               g(ig) = REAL(ig, dp)*dg(i)
               rag(ig) = CMPLX(fap(i), -g(ig), KIND=dp)
               cubeaxis(ig, i, 1, 0) = rag(ig)*cubeaxis(ig, i, 0, 0)
            END DO
            DO la = 2, la_max
               fa = REAL(la - 1, dp)*twozetp
               DO ig = lb_cube(i), ub_cube(i)
                  cubeaxis(ig, i, la, 0) = rag(ig)*cubeaxis(ig, i, la - 1, 0) + &
                                           fa*cubeaxis(ig, i, la - 2, 0)
               END DO
            END DO
            IF (lb_max > 0) THEN
               fa = twozetp
               DO ig = lb_cube(i), ub_cube(i)
                  rbg(ig) = CMPLX(fbp(i), -g(ig), KIND=dp)
                  cubeaxis(ig, i, 0, 1) = rbg(ig)*cubeaxis(ig, i, 0, 0)
                  cubeaxis(ig, i, 1, 1) = rbg(ig)*cubeaxis(ig, i, 1, 0) + &
                                          fa*cubeaxis(ig, i, 0, 0)
               END DO
               DO lb = 2, lb_max
                  fb = REAL(lb - 1, dp)*twozetp
                  DO ig = lb_cube(i), ub_cube(i)
                     cubeaxis(ig, i, 0, lb) = rbg(ig)*cubeaxis(ig, i, 0, lb - 1) + &
                                              fb*cubeaxis(ig, i, 0, lb - 2)
                     cubeaxis(ig, i, 1, lb) = rbg(ig)*cubeaxis(ig, i, 1, lb - 1) + &
                                              fb*cubeaxis(ig, i, 1, lb - 2) + &
                                              fa*cubeaxis(ig, i, 0, lb - 1)
                  END DO
               END DO
               DO la = 2, la_max
                  fa = REAL(la, dp)*twozetp
                  DO ig = lb_cube(i), ub_cube(i)
                     cubeaxis(ig, i, la, 1) = rbg(ig)*cubeaxis(ig, i, la, 0) + &
                                              fa*cubeaxis(ig, i, la - 1, 0)
                  END DO
                  DO lb = 2, lb_max
                     fb = REAL(lb - 1, dp)*twozetp
                     DO ig = lb_cube(i), ub_cube(i)
                        cubeaxis(ig, i, la, lb) = rbg(ig)*cubeaxis(ig, i, la, lb - 1) + &
                                                  fb*cubeaxis(ig, i, la, lb - 2) + &
                                                  fa*cubeaxis(ig, i, la - 1, lb - 1)
                     END DO
                  END DO
               END DO
            END IF
         ELSE
            IF (lb_max > 0) THEN
               DO ig = lb_cube(i), ub_cube(i)
                  g(ig) = REAL(ig, dp)*dg(i)
                  rbg(ig) = CMPLX(fbp(i), -g(ig), KIND=dp)
                  cubeaxis(ig, i, 0, 1) = rbg(ig)*cubeaxis(ig, i, 0, 0)
               END DO
               DO lb = 2, lb_max
                  fb = REAL(lb - 1, dp)*twozetp
                  DO ig = lb_cube(i), ub_cube(i)
                     cubeaxis(ig, i, 0, lb) = rbg(ig)*cubeaxis(ig, i, 0, lb - 1) + &
                                              fb*cubeaxis(ig, i, 0, lb - 2)
                  END DO
               END DO
            END IF
         END IF

      END DO

      DO la = 0, la_max
         DO lb = 0, lb_max
            IF (la + lb == 0) CYCLE
            fa = (1.0_dp/twozetp)**(la + lb)
            DO i = 1, 3
               DO ig = lb_cube(i), ub_cube(i)
                  cubeaxis(ig, i, la, lb) = fa*cubeaxis(ig, i, la, lb)
               END DO
            END DO
         END DO
      END DO

      ! Add the current primitive Gaussian function product to grid

      DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)

         ax = indco(1, ico)
         ay = indco(2, ico)
         az = indco(3, ico)

         DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)

            pij = prefactor*pab(na + ico, nb + jco)

            IF (ABS(pij) < eps_rho_gspace) CYCLE

            bx = indco(1, jco)
            by = indco(2, jco)
            bz = indco(3, jco)

            DO i = lb_grid, ub_grid
               IF (pw%pw_grid%gsq(i) > gsq_max) CYCLE
               ig = pw%pw_grid%g_hat(1, i)
               jg = pw%pw_grid%g_hat(2, i)
               kg = pw%pw_grid%g_hat(3, i)
               pw%array(i) = pw%array(i) + pij*cubeaxis(ig, 1, ax, bx)* &
                             cubeaxis(jg, 2, ay, by)* &
                             cubeaxis(kg, 3, az, bz)
            END DO

         END DO

      END DO

      DEALLOCATE (cubeaxis)
      DEALLOCATE (g)
      DEALLOCATE (rag)
      DEALLOCATE (rbg)

      CALL timestop(handle)

   END SUBROUTINE collocate_pgf_product_gspace

END MODULE xray_diffraction
