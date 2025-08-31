!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types used by CNEO-DFT
!>      (see J. Chem. Theory Comput. 2025, 21, 16, 7865–7877)
!> \par History
!>      08.2025 created [zc62]
!> \author Zehua Chen
! **************************************************************************************************
MODULE qs_cneo_types
   USE kinds,                           ONLY: dp
   USE periodic_table,                  ONLY: ptable
   USE qs_harmonics_atom,               ONLY: deallocate_harmonics_atom,&
                                              harmonics_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cneo_types'

   ! Essential matrices, density and potential for each quantum nucleus
   TYPE rhoz_cneo_type
      LOGICAL    :: ready = .FALSE.                ! if pmat is ready. Useful in first iter
      REAL(dp), DIMENSION(:, :), &
         POINTER :: pmat => Null(), &              ! nuclear density matrix
                    core => Null(), &              ! nuclear core Hamiltonian
                    vmat => Null(), &              ! nuclear Hartree from soft basis
                    fmat => Null(), &              ! Fock = core + Hartree
                    wfn => Null()                  ! nuclear orbital coefficients
      REAL(dp), DIMENSION(3) &
         :: f = (/0.0_dp, 0.0_dp, 0.0_dp/) ! Lagrange multiplier for CNEO
      REAL(dp)   :: e_core = 0.0_dp                ! nuclear core energy
      REAL(dp), DIMENSION(:, :), &
         POINTER :: cpc_h => Null(), &             ! decontracted density matrix
                    cpc_s => Null(), &             ! decontracted density matrix, soft tail
                    rho_rad_h => Null(), &         ! density on radial grid
                    rho_rad_s => Null(), &         ! density on radial grid, soft tail
                    vrho_rad_h => Null(), &        ! potential on radial grid
                    vrho_rad_s => Null(), &        ! potential on radial grid, soft tail
                    ga_Vlocal_gb_h => Null(), &    ! local Hartree integral
                    ga_Vlocal_gb_s => Null()       ! local Hartree integral, soft tail
   END TYPE rhoz_cneo_type

   ! S, T, transformation matrices and distance are shared by the same kind,
   ! since they are not affected by the position of basis center
   TYPE cneo_potential_type
      INTEGER    :: z = 0                          ! atomic number
      REAL(dp)   :: zeff = 0.0_dp, &               ! zeff = REAL(z)
                    mass = 0.0_dp                  ! atomic mass in dalton
      INTEGER, DIMENSION(:), &
         POINTER :: elec_conf => Null()
      INTEGER    :: nsgf = 0, &                    ! nuclear basis set is usually uncontracted
                    nne = 0, &                     ! number of linear-independent basis functions
                    npsgf = 0, &                   ! number of primitive SGFs
                    nsotot = 0                     ! maxso * nset
      REAL(dp), DIMENSION(:, :), &
         POINTER :: my_gcc_h => Null(), &          ! 3D-normalized contraction coefficients
                    my_gcc_s => Null(), &          ! contraction coefficients for the soft tail
                    ovlp => Null(), &              ! nuclear basis overlap matrix (unused)
                    kin => Null(), &               ! nuclear kinetic energy matrix
                    utrans => Null()               ! nuclear basis transformation matrix
      REAL(dp), DIMENSION(:, :, :), &
         POINTER :: distance => Null()             ! distance from nuclear basis center
      TYPE(harmonics_atom_type), &
         POINTER :: harmonics => Null()            ! most of the data will be missing
      REAL(dp), DIMENSION(:, :, :), &
         POINTER :: Qlm_gg => Null(), &            ! multipole expansion of nuclear gg
                    gg => Null()                   ! precompute and store gg on radial grid
      REAL(dp), DIMENSION(:, :, :, :), &
         POINTER :: vgg => Null()                  ! precompute and store vgg on grid
      INTEGER, DIMENSION(:), &
         POINTER :: n2oindex => Null(), &          ! new to old index
                    o2nindex => Null()             ! old to new index
      REAL(dp), DIMENSION(:, :), &
         POINTER :: rad2l => Null(), &             ! store my own rad2l
                    oorad2l => Null()              ! store my own oorad2l
   END TYPE cneo_potential_type

! Public Types

   PUBLIC :: cneo_potential_type, rhoz_cneo_type

! Public Subroutine

   PUBLIC :: allocate_cneo_potential, allocate_rhoz_cneo_set, deallocate_cneo_potential, &
             deallocate_rhoz_cneo_set, get_cneo_potential, set_cneo_potential, write_cneo_potential

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rhoz_cneo_set ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE allocate_rhoz_cneo_set(rhoz_cneo_set, natom)

      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set
      INTEGER, INTENT(IN)                                :: natom

      IF (ASSOCIATED(rhoz_cneo_set)) THEN
         CALL deallocate_rhoz_cneo_set(rhoz_cneo_set)
      END IF

      ALLOCATE (rhoz_cneo_set(natom))

   END SUBROUTINE allocate_rhoz_cneo_set

! **************************************************************************************************
!> \brief ...
!> \param rhoz_cneo ...
! **************************************************************************************************
   SUBROUTINE deallocate_rhoz_cneo(rhoz_cneo)

      TYPE(rhoz_cneo_type), POINTER                      :: rhoz_cneo

      IF (ASSOCIATED(rhoz_cneo)) THEN
         IF (ASSOCIATED(rhoz_cneo%pmat)) &
            DEALLOCATE (rhoz_cneo%pmat)
         IF (ASSOCIATED(rhoz_cneo%core)) &
            DEALLOCATE (rhoz_cneo%core)
         IF (ASSOCIATED(rhoz_cneo%vmat)) &
            DEALLOCATE (rhoz_cneo%vmat)
         IF (ASSOCIATED(rhoz_cneo%fmat)) &
            DEALLOCATE (rhoz_cneo%fmat)
         IF (ASSOCIATED(rhoz_cneo%wfn)) &
            DEALLOCATE (rhoz_cneo%wfn)
         IF (ASSOCIATED(rhoz_cneo%cpc_h)) &
            DEALLOCATE (rhoz_cneo%cpc_h)
         IF (ASSOCIATED(rhoz_cneo%cpc_s)) &
            DEALLOCATE (rhoz_cneo%cpc_s)
         IF (ASSOCIATED(rhoz_cneo%rho_rad_h)) &
            DEALLOCATE (rhoz_cneo%rho_rad_h)
         IF (ASSOCIATED(rhoz_cneo%rho_rad_s)) &
            DEALLOCATE (rhoz_cneo%rho_rad_s)
         IF (ASSOCIATED(rhoz_cneo%vrho_rad_h)) &
            DEALLOCATE (rhoz_cneo%vrho_rad_h)
         IF (ASSOCIATED(rhoz_cneo%vrho_rad_s)) &
            DEALLOCATE (rhoz_cneo%vrho_rad_s)
         IF (ASSOCIATED(rhoz_cneo%ga_Vlocal_gb_h)) &
            DEALLOCATE (rhoz_cneo%ga_Vlocal_gb_h)
         IF (ASSOCIATED(rhoz_cneo%ga_Vlocal_gb_s)) &
            DEALLOCATE (rhoz_cneo%ga_Vlocal_gb_s)
      END IF

   END SUBROUTINE deallocate_rhoz_cneo

! **************************************************************************************************
!> \brief ...
!> \param rhoz_cneo_set ...
! **************************************************************************************************
   SUBROUTINE deallocate_rhoz_cneo_set(rhoz_cneo_set)

      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set

      INTEGER                                            :: iat, natom
      TYPE(rhoz_cneo_type), POINTER                      :: rhoz_cneo

      IF (ASSOCIATED(rhoz_cneo_set)) THEN
         natom = SIZE(rhoz_cneo_set)
         DO iat = 1, natom
            rhoz_cneo => rhoz_cneo_set(iat)
            CALL deallocate_rhoz_cneo(rhoz_cneo)
         END DO
         DEALLOCATE (rhoz_cneo_set)
      END IF

   END SUBROUTINE deallocate_rhoz_cneo_set

! **************************************************************************************************
!> \brief ...
!> \param potential ...
! **************************************************************************************************
   SUBROUTINE allocate_cneo_potential(potential)

      TYPE(cneo_potential_type), POINTER                 :: potential

      IF (ASSOCIATED(potential)) &
         CALL deallocate_cneo_potential(potential)

      ALLOCATE (potential)

   END SUBROUTINE allocate_cneo_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
! **************************************************************************************************
   SUBROUTINE deallocate_cneo_potential(potential)

      TYPE(cneo_potential_type), POINTER                 :: potential

      IF (ASSOCIATED(potential)) THEN
         IF (ASSOCIATED(potential%elec_conf)) &
            DEALLOCATE (potential%elec_conf)
         IF (ASSOCIATED(potential%my_gcc_h)) &
            DEALLOCATE (potential%my_gcc_h)
         IF (ASSOCIATED(potential%my_gcc_s)) &
            DEALLOCATE (potential%my_gcc_s)
         IF (ASSOCIATED(potential%ovlp)) &
            DEALLOCATE (potential%ovlp)
         IF (ASSOCIATED(potential%kin)) &
            DEALLOCATE (potential%kin)
         IF (ASSOCIATED(potential%utrans)) &
            DEALLOCATE (potential%utrans)
         IF (ASSOCIATED(potential%distance)) &
            DEALLOCATE (potential%distance)
         IF (ASSOCIATED(potential%harmonics)) &
            CALL deallocate_harmonics_atom(potential%harmonics)
         IF (ASSOCIATED(potential%Qlm_gg)) &
            DEALLOCATE (potential%Qlm_gg)
         IF (ASSOCIATED(potential%gg)) &
            DEALLOCATE (potential%gg)
         IF (ASSOCIATED(potential%vgg)) &
            DEALLOCATE (potential%vgg)
         IF (ASSOCIATED(potential%n2oindex)) &
            DEALLOCATE (potential%n2oindex)
         IF (ASSOCIATED(potential%o2nindex)) &
            DEALLOCATE (potential%o2nindex)
         IF (ASSOCIATED(potential%rad2l)) &
            DEALLOCATE (potential%rad2l)
         IF (ASSOCIATED(potential%oorad2l)) &
            DEALLOCATE (potential%oorad2l)
         DEALLOCATE (potential)
      END IF

   END SUBROUTINE deallocate_cneo_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param z ...
!> \param zeff ...
!> \param mass ...
!> \param elec_conf ...
!> \param nsgf ...
!> \param nne ...
!> \param npsgf ...
!> \param nsotot ...
!> \param my_gcc_h ...
!> \param my_gcc_s ...
!> \param ovlp ...
!> \param kin ...
!> \param utrans ...
!> \param distance ...
!> \param harmonics ...
!> \param Qlm_gg ...
!> \param gg ...
!> \param vgg ...
!> \param n2oindex ...
!> \param o2nindex ...
!> \param rad2l ...
!> \param oorad2l ...
! **************************************************************************************************
   SUBROUTINE get_cneo_potential(potential, z, zeff, mass, elec_conf, nsgf, nne, npsgf, &
                                 nsotot, my_gcc_h, my_gcc_s, ovlp, kin, utrans, distance, &
                                 harmonics, Qlm_gg, gg, vgg, n2oindex, o2nindex, rad2l, oorad2l)

      TYPE(cneo_potential_type), POINTER                 :: potential
      INTEGER, INTENT(OUT), OPTIONAL                     :: z
      REAL(dp), INTENT(OUT), OPTIONAL                    :: zeff, mass
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf
      INTEGER, INTENT(OUT), OPTIONAL                     :: nsgf, nne, npsgf, nsotot
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: my_gcc_h, my_gcc_s, ovlp, kin, utrans
      REAL(dp), DIMENSION(:, :, :), OPTIONAL, POINTER    :: distance
      TYPE(harmonics_atom_type), OPTIONAL, POINTER       :: harmonics
      REAL(dp), DIMENSION(:, :, :), OPTIONAL, POINTER    :: Qlm_gg, gg
      REAL(dp), DIMENSION(:, :, :, :), OPTIONAL, POINTER :: vgg
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: n2oindex, o2nindex
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: rad2l, oorad2l

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(z)) z = potential%z
         IF (PRESENT(zeff)) zeff = potential%zeff
         IF (PRESENT(mass)) mass = potential%mass
         IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf
         IF (PRESENT(nsgf)) nsgf = potential%nsgf
         IF (PRESENT(nne)) nne = potential%nne
         IF (PRESENT(npsgf)) npsgf = potential%npsgf
         IF (PRESENT(nsotot)) nsotot = potential%nsotot
         IF (PRESENT(my_gcc_h)) my_gcc_h => potential%my_gcc_h
         IF (PRESENT(my_gcc_s)) my_gcc_s => potential%my_gcc_s
         IF (PRESENT(ovlp)) ovlp => potential%ovlp
         IF (PRESENT(kin)) kin => potential%kin
         IF (PRESENT(ovlp)) ovlp => potential%ovlp
         IF (PRESENT(utrans)) utrans => potential%utrans
         IF (PRESENT(distance)) distance => potential%distance
         IF (PRESENT(harmonics)) harmonics => potential%harmonics
         IF (PRESENT(Qlm_gg)) Qlm_gg => potential%Qlm_gg
         IF (PRESENT(gg)) gg => potential%gg
         IF (PRESENT(vgg)) vgg => potential%vgg
         IF (PRESENT(n2oindex)) n2oindex => potential%n2oindex
         IF (PRESENT(o2nindex)) o2nindex => potential%o2nindex
         IF (PRESENT(rad2l)) rad2l => potential%rad2l
         IF (PRESENT(oorad2l)) oorad2l => potential%oorad2l

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE get_cneo_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param z ...
!> \param mass ...
!> \param elec_conf ...
!> \param nsgf ...
!> \param nne ...
!> \param npsgf ...
!> \param nsotot ...
!> \param my_gcc_h ...
!> \param my_gcc_s ...
!> \param ovlp ...
!> \param kin ...
!> \param utrans ...
!> \param distance ...
!> \param harmonics ...
!> \param Qlm_gg ...
!> \param gg ...
!> \param vgg ...
!> \param n2oindex ...
!> \param o2nindex ...
!> \param rad2l ...
!> \param oorad2l ...
! **************************************************************************************************
   SUBROUTINE set_cneo_potential(potential, z, mass, elec_conf, nsgf, nne, npsgf, &
                                 nsotot, my_gcc_h, my_gcc_s, ovlp, kin, utrans, distance, &
                                 harmonics, Qlm_gg, gg, vgg, n2oindex, o2nindex, rad2l, oorad2l)

      TYPE(cneo_potential_type), POINTER                 :: potential
      INTEGER, INTENT(IN), OPTIONAL                      :: z
      REAL(dp), INTENT(IN), OPTIONAL                     :: mass
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf
      INTEGER, INTENT(IN), OPTIONAL                      :: nsgf, nne, npsgf, nsotot
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: my_gcc_h, my_gcc_s, ovlp, kin, utrans
      REAL(dp), DIMENSION(:, :, :), OPTIONAL, POINTER    :: distance
      TYPE(harmonics_atom_type), OPTIONAL, POINTER       :: harmonics
      REAL(dp), DIMENSION(:, :, :), OPTIONAL, POINTER    :: Qlm_gg, gg
      REAL(dp), DIMENSION(:, :, :, :), OPTIONAL, POINTER :: vgg
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: n2oindex, o2nindex
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: rad2l, oorad2l

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(z)) THEN
            potential%z = z
            potential%zeff = REAL(z, dp)
            IF (ASSOCIATED(potential%elec_conf)) &
               CPABORT("elec_conf is already associated")
            ALLOCATE (potential%elec_conf(0:3))
            potential%elec_conf(0:3) = ptable(z)%e_conv(0:3)
            CPASSERT(potential%mass == 0.0_dp)
            IF (z == 1) THEN
               ! Hydrogen is 1.007825, not 1.00794
               ! subtract the electron mass to get the proton mass
               potential%mass = 1.007825_dp - 0.000548579909_dp
            ELSE
               ! In principle, the most abundant pure isotope mass
               ! should be used, but no such data is available in ptable
               potential%mass = ptable(z)%amass - 0.000548579909_dp*REAL(z, dp)
            END IF
         END IF
         IF (PRESENT(mass)) THEN
            potential%mass = mass
         END IF
         IF (PRESENT(elec_conf)) THEN
            IF (ASSOCIATED(potential%elec_conf)) THEN
               DEALLOCATE (potential%elec_conf)
            END IF
            ALLOCATE (potential%elec_conf(0:SIZE(elec_conf) - 1))
            potential%elec_conf(:) = elec_conf(:)
         END IF
         IF (PRESENT(nsgf)) potential%nsgf = nsgf
         IF (PRESENT(nne)) potential%nne = nne
         IF (PRESENT(npsgf)) potential%npsgf = npsgf
         IF (PRESENT(nsotot)) potential%nsotot = nsotot
         IF (PRESENT(my_gcc_h)) potential%my_gcc_h => my_gcc_h
         IF (PRESENT(my_gcc_s)) potential%my_gcc_s => my_gcc_s
         IF (PRESENT(ovlp)) potential%ovlp => ovlp
         IF (PRESENT(kin)) potential%kin => kin
         IF (PRESENT(utrans)) potential%utrans => utrans
         IF (PRESENT(distance)) potential%distance => distance
         IF (PRESENT(harmonics)) potential%harmonics => harmonics
         IF (PRESENT(Qlm_gg)) potential%Qlm_gg => Qlm_gg
         IF (PRESENT(gg)) potential%gg => gg
         IF (PRESENT(vgg)) potential%vgg => vgg
         IF (PRESENT(n2oindex)) potential%n2oindex => n2oindex
         IF (PRESENT(o2nindex)) potential%o2nindex => o2nindex
         IF (PRESENT(rad2l)) potential%rad2l => rad2l
         IF (PRESENT(oorad2l)) potential%oorad2l => oorad2l

      ELSE

         CPABORT("The pointer potential is not associated")

      END IF

   END SUBROUTINE set_cneo_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE write_cneo_potential(potential, output_unit)

      TYPE(cneo_potential_type), POINTER                 :: potential
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=20)                                  :: string

      IF (output_unit > 0 .AND. ASSOCIATED(potential)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T6,A,/)") &
            "CNEO Potential information"
         WRITE (UNIT=output_unit, FMT="(T8,A,T41,A,I4,A,F11.6)") &
            "Description: ", "Z =", potential%z, &
            ", nuclear mass =", potential%mass
         WRITE (UNIT=string, FMT="(5I4)") potential%elec_conf
         WRITE (UNIT=output_unit, FMT="(T8,A,T61,A20)") &
            "Electronic configuration (s p d ...):", &
            ADJUSTR(TRIM(string))
      END IF

   END SUBROUTINE write_cneo_potential

END MODULE qs_cneo_types
