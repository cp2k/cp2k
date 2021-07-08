!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Implementation of the GAL21 potential
!>
!> \author Clabaut Paul
! **************************************************************************************************
MODULE manybody_gal21

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: pos_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: matvec_3x3
   USE message_passing,                 ONLY: mp_sum
   USE pair_potential_types,            ONLY: gal21_pot_type,&
                                              gal21_type,&
                                              pair_potential_pp_type,&
                                              pair_potential_single_type
   USE particle_types,                  ONLY: particle_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: setup_gal21_arrays, destroy_gal21_arrays, &
             gal21_energy, gal21_forces, &
             print_nr_ions_gal21
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_gal21'

CONTAINS

! **************************************************************************************************
!> \brief  Main part of the energy evaluation of GAL2119
!> \param pot_loc value of total potential energy
!> \param gal21 all parameters of GAL2119
!> \param r_last_update_pbc position of every atoms on previous frame
!> \param iparticle first index of the atom of the evaluated pair
!> \param jparticle second index of the atom of the evaluated pair
!> \param cell dimension of the pbc cell
!> \param particle_set full list of atoms of the system
!> \param mm_section ...
!> \author Clabaut Paul - 2019 - ENS de Lyon
! **************************************************************************************************
   SUBROUTINE gal21_energy(pot_loc, gal21, r_last_update_pbc, iparticle, jparticle, &
                         cell, particle_set, mm_section)

      REAL(KIND=dp), INTENT(OUT)                         :: pot_loc
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: iparticle, jparticle
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: mm_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gal21_energy', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: index_outfile
      REAL(KIND=dp)                                      :: anglepart, cosalpha, drji2, &
                                                            nvec(3), rji(3), &
                                                            sinalpha, sum_weight, Vang, Vgaussian, &
                                                            VTT, VH, weight, eps, bxy, bz, AO, BO
      TYPE(cp_logger_type), POINTER                      :: logger

      pot_loc = 0.0_dp
      CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind, &
                           element_symbol=element_symbol) !Read the atom type of i

      IF (element_symbol == "O") THEN !To avoid counting two times each pair

         rji(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(iparticle)%r(:), cell) !Vector in pbc from j to i

         IF (.NOT. ALLOCATED(gal21%n_vectors)) THEN !First calling of the forcefield only
            ALLOCATE (gal21%n_vectors(3, SIZE(particle_set)))
            gal21%n_vectors(:, :) = 0.0_dp
         END IF


         IF (gal21%express) THEN
                  logger => cp_get_default_logger()
                  index_outfile = cp_print_key_unit_nr(logger, mm_section, &
                                                       "PRINT%PROGRAM_RUN_INFO", extension=".mmLog")
                  IF (index_outfile > 0) WRITE (index_outfile, *) "GCN", gal21%gcn(jparticle)
                  CALL cp_print_key_finished_output(index_outfile, logger, mm_section, &
                                                    "PRINT%PROGRAM_RUN_INFO")
         END IF

         !Build epsilon attraction and the parameters of the gaussian attraction as a function of gcn
         eps = gal21%epsilon1*gal21%gcn(jparticle)*gal21%gcn(jparticle) + gal21%epsilon2*gal21%gcn(jparticle) + gal21%epsilon3
         bxy = gal21%bxy1 + gal21%bxy2*gal21%gcn(jparticle)
         bz = gal21%bz1 + gal21%bz2*gal21%gcn(jparticle)

         !Angular dependance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Vang = 0.0_dp

         !Calculation of the normal vector centered on the Me atom of the pair, only the first time that an interaction with the metal atom of the pair is evaluated
         IF (gal21%n_vectors(1, jparticle) == 0.0_dp .AND. &
             gal21%n_vectors(2, jparticle) == 0.0_dp .AND. &
             gal21%n_vectors(3, jparticle) == 0.0_dp) THEN
            gal21%n_vectors(:, jparticle) = normale(gal21, r_last_update_pbc, jparticle, &
                                                  particle_set, cell)
         END IF

         nvec(:) = gal21%n_vectors(:, jparticle) !Else, retrive it, should not have moved sinc metal is supposed to be frozen

         !Calculation of the sum of the expontial weights of each Me surrounding the principal one
         sum_weight = somme(gal21, r_last_update_pbc, iparticle, particle_set, cell)

         !Exponential damping weight for angular dependance
         weight = EXP(-SQRT(DOT_PRODUCT(rji, rji))/gal21%r1)

         !Calculation of the truncated fourier series of the water-dipole/surface-normal angle
         anglepart = 0.0_dp
         VH = 0.0_dp
         CALL angular(anglepart, VH, gal21, r_last_update_pbc, iparticle, jparticle, cell, particle_set, nvec, &
                       .TRUE., mm_section)

         !Build the complete angular potential while avoiding division by 0
         IF (weight /= 0) THEN
            Vang = weight*weight*anglepart/sum_weight
            IF (gal21%express) THEN
               logger => cp_get_default_logger()
               index_outfile = cp_print_key_unit_nr(logger, mm_section, &
                                                    "PRINT%PROGRAM_RUN_INFO", extension=".mmLog")
               IF (index_outfile > 0) WRITE (index_outfile, *) "Fermi", weight*weight/sum_weight
               CALL cp_print_key_finished_output(index_outfile, logger, mm_section, &
                                                 "PRINT%PROGRAM_RUN_INFO")
            END IF
         END IF
         !END Angular%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !Attractive Gaussian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Vgaussian = 0.0_dp
         drji2 = DOT_PRODUCT(rji, rji)
         !Alpha is the angle of the Me-O vector with the normale vector. Used for gaussian attaction

         cosalpha = DOT_PRODUCT(rji, nvec)/SQRT(drji2)
         IF (cosalpha < -1.0_dp) cosalpha = -1.0_dp
         IF (cosalpha > +1.0_dp) cosalpha = +1.0_dp
         sinalpha = SIN(ACOS(cosalpha))

         !Gaussian component of the energy
         Vgaussian = -1.0_dp*eps*EXP(-bz*drji2*cosalpha*cosalpha &
                                                         - bxy*drji2*sinalpha*sinalpha)
         !END Gaussian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         AO = gal21%AO1 + gal21%AO2*gal21%gcn(jparticle)
         BO = gal21%BO1 + gal21%BO2*gal21%gcn(jparticle)

         !Tang and toennies potential for physisorption
         VTT = AO*EXP(-BO*SQRT(drji2)) - (1.0 - EXP(-BO*SQRT(drji2)) &
                                                - BO*SQRT(drji2)*EXP(-BO*SQRT(drji2)) &
                                                - (((BO*SQRT(drji2))**2)/2)*EXP(-BO*SQRT(drji2)) &
                                                - (((BO*SQRT(drji2))**3)/6)*EXP(-BO*SQRT(drji2)) &
                                                - (((BO*SQRT(drji2))**4)/24)*EXP(-BO*SQRT(drji2)) &
                                                - (((BO*SQRT(drji2))**5)/120)*EXP(-BO*SQRT(drji2)) &
                                                - (((BO*SQRT(drji2))**6)/720)*EXP(-BO*SQRT(drji2))) &
               *gal21%c/(SQRT(drji2)**6)

         !For fit purpose only
         IF (gal21%express) THEN
            logger => cp_get_default_logger()
            index_outfile = cp_print_key_unit_nr(logger, mm_section, &
                                                 "PRINT%PROGRAM_RUN_INFO", extension=".mmLog")
            IF (index_outfile > 0) WRITE (index_outfile, *) "Gau", -1.0_dp*EXP(-bz*drji2*cosalpha*cosalpha &
                                                                                           - bxy*drji2*sinalpha*sinalpha)
            IF (weight == 0 .AND. index_outfile > 0) WRITE (index_outfile, *) "Fermi  0"
            IF (index_outfile > 0) WRITE (index_outfile, *) "expO", EXP(-BO*SQRT(drji2))
            IF (index_outfile > 0) WRITE (index_outfile, *) "cstpart", -(1.0 - EXP(-BO*SQRT(drji2)) &
                                                                         - BO*SQRT(drji2)*EXP(-BO*SQRT(drji2)) &
                                                                         - (((BO*SQRT(drji2))**2)/2)*EXP(-BO*SQRT(drji2)) &
                                                                         - (((BO*SQRT(drji2))**3)/6)*EXP(-BO*SQRT(drji2)) &
                                                                         - (((BO*SQRT(drji2))**4)/24)*EXP(-BO*SQRT(drji2)) &
                                                                         - (((BO*SQRT(drji2))**5)/120)*EXP(-BO*SQRT(drji2)) &
                                                                         - (((BO*SQRT(drji2))**6)/720)*EXP(-BO*SQRT(drji2))) &
               *gal21%c/(SQRT(drji2)**6)
            IF (index_outfile > 0) WRITE (index_outfile, *) "params_lin_eps", gal21%epsilon1, gal21%epsilon2, gal21%epsilon3
            IF (index_outfile > 0) WRITE (index_outfile, *) "params_lin_A0", AO 
            CALL cp_print_key_finished_output(index_outfile, logger, mm_section, &
                                              "PRINT%PROGRAM_RUN_INFO")
         END IF
         !Compute the total energy
         pot_loc = Vgaussian + Vang + VTT + VH

      END IF

   END SUBROUTINE gal21_energy

! **************************************************************************************************
! The idea is to build a vector normal to the local surface by using the symetry of the surface that
! make the opposite vectors compensate themself. The vector is therefore in the direction of the
! missing atoms of a large coordination sphere
! **************************************************************************************************
!> \brief ...
!> \param gal21 ...
!> \param r_last_update_pbc ...
!> \param jparticle ...
!> \param particle_set ...
!> \param cell ...
!> \return ...
!> \retval normale ...
! **************************************************************************************************
   FUNCTION normale(gal21, r_last_update_pbc, jparticle, particle_set, cell)
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: jparticle
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp)                                      :: normale(3)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'normale', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol_k
      INTEGER                                            :: kparticle, natom
      REAL(KIND=dp)                                      :: drjk, rjk(3)

      natom = SIZE(particle_set)
      normale(:) = 0.0_dp

      DO kparticle = 1, natom !Loop on every atom of the system
         IF (kparticle == jparticle) CYCLE !Avoid the principal Me atom (j) in the counting
         CALL get_atomic_kind(atomic_kind=particle_set(kparticle)%atomic_kind, &
                              element_symbol=element_symbol_k)
         IF (element_symbol_k /= gal21%met1 .AND. element_symbol_k /= gal21%met2) CYCLE !Keep only metals
         rjk(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(kparticle)%r(:), cell)
         drjk = SQRT(DOT_PRODUCT(rjk, rjk))
         IF (drjk > gal21%rcutsq) CYCLE !Keep only those within square root of the force-field cutoff distance of the metallic atom of the evaluated pair
         normale(:) = normale(:)-rjk(:)/(drjk*drjk*drjk*drjk*drjk) !Build the normal, vector by vector
      END DO

      ! Normalisation of the vector
      normale(:) = normale(:)/SQRT(DOT_PRODUCT(normale, normale))

   END FUNCTION normale

! **************************************************************************************************
! Scan all the Me atoms that have been counted in the O-Me paires and sum their exponential weights
! **************************************************************************************************
!> \brief ...
!> \param gal21 ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param particle_set ...
!> \param cell ...
!> \return ...
!> \retval somme ...
! **************************************************************************************************
   FUNCTION somme(gal21, r_last_update_pbc, iparticle, particle_set, cell)
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: iparticle
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp)                                      :: somme

      CHARACTER(LEN=*), PARAMETER :: routineN = 'somme', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol_k
      INTEGER                                            :: kparticle, natom
      REAL(KIND=dp)                                      :: rki(3)

      natom = SIZE(particle_set)
      somme = 0.0_dp

      DO kparticle = 1, natom !Loop on every atom of the system
         CALL get_atomic_kind(atomic_kind=particle_set(kparticle)%atomic_kind, &
                              element_symbol=element_symbol_k)
         IF (element_symbol_k /= gal21%met1 .AND. element_symbol_k /= gal21%met2) CYCLE !Keep only metals
         rki(:) = pbc(r_last_update_pbc(kparticle)%r(:), r_last_update_pbc(iparticle)%r(:), cell)
         IF (SQRT(DOT_PRODUCT(rki, rki)) > gal21%rcutsq) CYCLE !Keep only those within cutoff distance of the oxygen atom of the evaluated pair (the omega ensemble)
         IF (element_symbol_k == gal21%met1) somme = somme + EXP(-SQRT(DOT_PRODUCT(rki, rki))/gal21%r1) !Build the sum of the exponential weights
         IF (element_symbol_k == gal21%met2) somme = somme + EXP(-SQRT(DOT_PRODUCT(rki, rki))/gal21%r2) !Build the sum of the exponential weights
      END DO

   END FUNCTION somme

! **************************************************************************************************

! **************************************************************************************************
! Compute the angular dependance (on theta) of the forcefield
! **************************************************************************************************
!> \brief ...
!> \param anglepart ...
!> \param VH ...
!> \param gal21 ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param cell ...
!> \param particle_set ...
!> \param nvec ...
!> \param energy ...
!> \param mm_section ...
!> \return ...
!> \retval angular ...
! **************************************************************************************************
   SUBROUTINE angular(anglepart, VH, gal21, r_last_update_pbc, iparticle, jparticle, cell, &
                                                        particle_set, nvec, energy, mm_section) 
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: iparticle, jparticle
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(3)                        :: nvec
      LOGICAL                                            :: energy
      TYPE(section_vals_type), POINTER                   :: mm_section
      REAL(KIND=dp)                                      :: anglepart, VH

      CHARACTER(LEN=*), PARAMETER :: routineN = 'angular', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: count_h, iatom, index_h1, index_h2, &
                                                            index_outfile, natom
      REAL(KIND=dp)                                      :: costheta, h_max_dist, rih(3), rih1(3), &
                                                            rih2(3), rix(3), theta, rjh1(3), rjh2(3),&
                                                            a1, a2, a3, a4, BH
      TYPE(cp_logger_type), POINTER                      :: logger

      count_h = 0
      index_h1 = 0
      index_h2 = 0
      h_max_dist = 2.1_dp ! 1.1 angstrom
      natom = SIZE(particle_set)

      DO iatom = 1, natom !Loop on every atom of the system
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              element_symbol=element_symbol)
         IF (element_symbol /= "H") CYCLE !Kepp only hydrogen
         rih(:) = pbc(r_last_update_pbc(iparticle)%r(:), r_last_update_pbc(iatom)%r(:), cell)
         IF (SQRT(DOT_PRODUCT(rih, rih)) >= h_max_dist) CYCLE !Keep only hydrogen that are bounded to the considered O
         count_h = count_h + 1
         IF (count_h == 1) THEN
            index_h1 = iatom
         ELSEIF (count_h == 2) THEN
            index_h2 = iatom
         ENDIF
      ENDDO

      ! Abort if the oxygen is not part of a water molecule (2 H)
      IF (count_h /= 2) THEN
         CALL cp_abort(__LOCATION__, &
                       " Error: Found "//cp_to_string(count_h)//" H atoms for O atom "//cp_to_string(iparticle))
      ENDIF

      a1 = gal21%a11 + gal21%a12*gal21%gcn(jparticle) + gal21%a13*gal21%gcn(jparticle)*gal21%gcn(jparticle)
      a2 = gal21%a21 + gal21%a22*gal21%gcn(jparticle) + gal21%a23*gal21%gcn(jparticle)*gal21%gcn(jparticle)
      a3 = gal21%a31 + gal21%a32*gal21%gcn(jparticle) + gal21%a33*gal21%gcn(jparticle)*gal21%gcn(jparticle)
      a4 = gal21%a41 + gal21%a42*gal21%gcn(jparticle) + gal21%a43*gal21%gcn(jparticle)*gal21%gcn(jparticle)

      rih1(:) = pbc(r_last_update_pbc(iparticle)%r(:), r_last_update_pbc(index_h1)%r(:), cell)
      rih2(:) = pbc(r_last_update_pbc(iparticle)%r(:), r_last_update_pbc(index_h2)%r(:), cell)
      rix(:) = rih1(:) + rih2(:) ! build the dipole vector rix of the H2O molecule
      costheta = DOT_PRODUCT(rix, nvec)/SQRT(DOT_PRODUCT(rix, rix))
      IF (costheta < -1.0_dp) costheta = -1.0_dp
      IF (costheta > +1.0_dp) costheta = +1.0_dp
      theta = ACOS(costheta) ! Theta is the angle between the normal to the surface and the dipole
      anglepart = a1*costheta + a2*COS(2.0_dp*theta) + a3*COS(3.0_dp*theta) &
                + a4*COS(4.0_dp*theta) ! build the fourier series

      BH = gal21%BH1 + gal21%gcn(jparticle)*gal21%BH2

      rjh1(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(index_h1)%r(:), cell)
      rjh2(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(index_h2)%r(:), cell)
      VH = (gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)*(exp(-BH*SQRT(DOT_PRODUCT(rjh1,rjh1))) + &
                                                    exp(-BH*SQRT(DOT_PRODUCT(rjh2,rjh2))) ) 
           

      ! For fit purpose
      IF (gal21%express .AND. energy) THEN
         logger => cp_get_default_logger()
         index_outfile = cp_print_key_unit_nr(logger, mm_section, &
                                              "PRINT%PROGRAM_RUN_INFO", extension=".mmLog")

         IF (index_outfile > 0) WRITE (index_outfile, *) "Fourier", costheta, COS(2.0_dp*theta), COS(3.0_dp*theta), &
            COS(4.0_dp*theta) !, theta
         IF (index_outfile > 0) WRITE (index_outfile, *) "H_rep", exp(-BH*SQRT(DOT_PRODUCT(rjh1, rjh1))) + &
                                  exp(-BH*SQRT(DOT_PRODUCT(rjh2, rjh2)))
         !IF (index_outfile > 0) WRITE (index_outfile, *) "H_r6", -1/DOT_PRODUCT(rjh1,rjh1)**3 -1/DOT_PRODUCT(rjh2,rjh2)**3

         CALL cp_print_key_finished_output(index_outfile, logger, mm_section, &
                                           "PRINT%PROGRAM_RUN_INFO")
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief forces generated by the GAL2119 potential
!> \param gal21 all parameters of GAL2119
!> \param r_last_update_pbc position of every atoms on previous frame
!> \param iparticle first index of the atom of the evaluated pair
!> \param jparticle second index of the atom of the evaluated pair
!> \param f_nonbond all the forces applying on the system
!> \param pv_nonbond ...
!> \param use_virial request of usage of virial (for barostat)
!> \param cell dimension of the pbc cell
!> \param particle_set full list of atoms of the system
!> \author Clabaut Paul - 2019 - ENS de Lyon
! **************************************************************************************************
   SUBROUTINE gal21_forces(gal21, r_last_update_pbc, iparticle, jparticle, f_nonbond, pv_nonbond, &
                          use_virial, cell, particle_set)
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: iparticle, jparticle
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gal21_forces', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      REAL(KIND=dp) :: anglepart, VH, cosalpha, dGauss(3), drji, drjicosalpha(3), drjisinalpha(3), &
         dTT(3), dweight(3), nvec(3), prefactor, rji(3), rji_hat(3), &
         sinalpha, sum_weight, Vgaussian, weight, eps, bxy, bz, AO, BO
      TYPE(section_vals_type), POINTER                   :: mm_section

      CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind, &
                           element_symbol=element_symbol)

      IF (element_symbol == "O") THEN !To avoid counting two times each pair

         rji(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(iparticle)%r(:), cell)
         drji = SQRT(DOT_PRODUCT(rji, rji))
         rji_hat(:) = rji(:)/drji ! hat = pure directional component of a given vector

         IF (.NOT. ALLOCATED(gal21%n_vectors)) THEN !First calling of the forcefield only
            ALLOCATE (gal21%n_vectors(3, SIZE(particle_set)))
            gal21%n_vectors(:, :) = 0.0_dp
         END IF


         !Build epsilon attraction and the a parameters of the Fourier serie as quadratic fucntion of gcn
         eps = gal21%epsilon1*gal21%gcn(jparticle)*gal21%gcn(jparticle) + gal21%epsilon2*gal21%gcn(jparticle) + gal21%epsilon3
         bxy = gal21%bxy1 + gal21%bxy2*gal21%gcn(jparticle)
         bz = gal21%bz1 + gal21%bz2*gal21%gcn(jparticle)

         !Angular dependance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !Calculation of the normal vector centered on the Me atom of the pair, only the first time that an interaction with the metal atom of the pair is evaluated
         IF (gal21%n_vectors(1, jparticle) == 0.0_dp .AND. &
             gal21%n_vectors(2, jparticle) == 0.0_dp .AND. &
             gal21%n_vectors(3, jparticle) == 0.0_dp) THEN
            gal21%n_vectors(:, jparticle) = normale(gal21, r_last_update_pbc, jparticle, &
                                                  particle_set, cell)
         END IF

         nvec(:) = gal21%n_vectors(:, jparticle) !Else, retrive it, should not have moved sinc metal is supposed to be frozen

         !Calculation of the sum of the expontial weights of each Me surrounding the principal one
         sum_weight = somme(gal21, r_last_update_pbc, iparticle, particle_set, cell)

         !Exponential damping weight for angular dependance
         weight = EXP(-drji/gal21%r1)
         dweight(:) = 1.0_dp/gal21%r1*weight*rji_hat(:) !Derivativ of it

         !Calculation of the truncated fourier series of the water-dipole/surface-normal angle
         NULLIFY (mm_section)
         anglepart = 0.0_dp
         VH = 0.0_dp
         CALL angular(anglepart, VH, gal21, r_last_update_pbc, iparticle, jparticle, cell, particle_set, nvec, &
                             .FALSE., mm_section)

         !Build the average of the exponential weight while avoiding division by 0
         IF (weight /= 0) THEN
            ! Calculate the first component of the derivativ of the angular term
            f_nonbond(1:3, iparticle) = f_nonbond(1:3, iparticle) + 2.0_dp*dweight(1:3)*weight* &
                                        anglepart/sum_weight

            IF (use_virial) THEN
              pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) + rji(1)*2.0_dp*dweight(1:3)*weight* &
                                        anglepart/sum_weight
              pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) + rji(2)*2.0_dp*dweight(1:3)*weight* &
                                        anglepart/sum_weight
              pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) + rji(3)*2.0_dp*dweight(1:3)*weight* &
                                        anglepart/sum_weight
            END IF

            ! Calculate the second component of the derivativ of the angular term
            CALL somme_d(gal21, r_last_update_pbc, iparticle, jparticle, &
                 f_nonbond, pv_nonbond, use_virial, particle_set, cell, anglepart, sum_weight)

            prefactor = (-1.0_dp)*weight*weight/sum_weight ! Avoiding division by 0

            ! Calculate the third component of the derivativ of the angular term
            CALL angular_d(gal21, r_last_update_pbc, iparticle, jparticle, &
                 f_nonbond, pv_nonbond, use_virial, prefactor, cell, particle_set, nvec)
         END IF
         !END Angular%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !Attractive Gaussian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !Alpha is the angle of the Me-O vector with the normale vector. Used for gaussian attaction
         cosalpha = DOT_PRODUCT(rji, nvec)/drji
         IF (cosalpha < -1.0_dp) cosalpha = -1.0_dp
         IF (cosalpha > +1.0_dp) cosalpha = +1.0_dp
         sinalpha = SIN(ACOS(cosalpha))

         !Gaussian component of the energy
         Vgaussian = -1.0_dp*eps*EXP(-bz*DOT_PRODUCT(rji, rji)*cosalpha*cosalpha &
                                                         - bxy*DOT_PRODUCT(rji, rji)*sinalpha*sinalpha)

         ! Calculation of partial derivativ of the gaussian components
         drjicosalpha(:) = rji_hat(:)*cosalpha + nvec(:) - cosalpha*rji_hat(:)
         drjisinalpha(:) = rji_hat(:)*sinalpha - (cosalpha/sinalpha)*(nvec(:) - cosalpha*rji_hat(:))
         dGauss(:) = (-1.0_dp*bz*2*drji*cosalpha*drjicosalpha - &
                      1.0_dp*bxy*2*drji*sinalpha*drjisinalpha)*Vgaussian*(-1.0_dp)

         ! Force due to gaussian term
         f_nonbond(1:3, iparticle) = f_nonbond(1:3, iparticle) + dGauss(1:3)

         IF (use_virial) THEN
           pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) + rji(1)*dGauss(1:3)
           pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) + rji(2)*dGauss(1:3)
           pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) + rji(3)*dGauss(1:3)
         END IF
         !END Gaussian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         AO = gal21%AO1 + gal21%AO2*gal21%gcn(jparticle)
         BO = gal21%BO1 + gal21%BO2*gal21%gcn(jparticle)

         !Derivativ of the Tang and Toennies term
         dTT(:) = (-(AO*BO + (BO**7)*gal21%c/720)*EXP(-BO*drji) + 6*(gal21%c/drji**7)* &
                   (1.0 - EXP(-BO*drji) &
                    - BO*drji*EXP(-BO*drji) &
                    - (((BO*drji)**2)/2)*EXP(-BO*drji) &
                    - (((BO*drji)**3)/6)*EXP(-BO*drji) &
                    - (((BO*drji)**4)/24)*EXP(-BO*drji) &
                    - (((BO*drji)**5)/120)*EXP(-BO*drji) &
                    - (((BO*drji)**6)/720)*EXP(-BO*drji)) &
                   )*rji_hat(:)

         ! Force of Tang & Toennies
         f_nonbond(1:3, iparticle) = f_nonbond(1:3, iparticle) - dTT(1:3)

         IF (use_virial) THEN
           pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) - rji(1)*dTT(1:3)
           pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) - rji(2)*dTT(1:3)
           pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) - rji(3)*dTT(1:3)
         END IF


      END IF

   END SUBROUTINE gal21_forces
! **************************************************************************************************
! Derivativ of the second component of angular dependance
! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param gal21 ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
!> \param particle_set ...
!> \param cell ...
!> \param anglepart ...
!> \param sum_weight ...
! **************************************************************************************************
   SUBROUTINE somme_d(gal21, r_last_update_pbc, iparticle, jparticle, &
                 f_nonbond, pv_nonbond, use_virial, particle_set, cell, anglepart, sum_weight)
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: iparticle, jparticle
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: anglepart, sum_weight

      CHARACTER(LEN=*), PARAMETER :: routineN = 'somme_d', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol_k
      INTEGER                                            :: kparticle, natom
      REAL(KIND=dp)                                      :: drki, dwdr(3), rji(3), rki(3), &
                                                            rki_hat(3), weight_rji

      rji(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(iparticle)%r(:), cell)
      weight_rji = EXP(-SQRT(DOT_PRODUCT(rji, rji))/gal21%r1)

      natom = SIZE(particle_set)
      DO kparticle = 1, natom !Loop on every atom of the system
         CALL get_atomic_kind(atomic_kind=particle_set(kparticle)%atomic_kind, &
                              element_symbol=element_symbol_k)
         IF (element_symbol_k /= gal21%met1 .AND. element_symbol_k /= gal21%met2) CYCLE !Keep only metals
         rki(:) = pbc(r_last_update_pbc(kparticle)%r(:), r_last_update_pbc(iparticle)%r(:), cell)
         IF (SQRT(DOT_PRODUCT(rki, rki)) > gal21%rcutsq) CYCLE !Keep only those within cutoff distance of the oxygen atom of the evaluated pair (the omega ensemble)
         drki = SQRT(DOT_PRODUCT(rki, rki))
         rki_hat(:) = rki(:)/drki

         IF (element_symbol_k == gal21%met1) dwdr(:) = (-1.0_dp)*(1.0_dp/gal21%r1)*EXP(-drki/gal21%r1)*rki_hat(:) !Build the sum of derivativs
         IF (element_symbol_k == gal21%met2) dwdr(:) = (-1.0_dp)*(1.0_dp/gal21%r2)*EXP(-drki/gal21%r2)*rki_hat(:)

         f_nonbond(1:3, iparticle) = f_nonbond(1:3, iparticle) + dwdr(1:3)*weight_rji &
                                     *weight_rji*anglepart/(sum_weight**2)

         IF (use_virial) THEN
           pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) + rki(1)*dwdr(1:3)*weight_rji &
                                       *weight_rji*anglepart/(sum_weight**2)
           pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) + rki(2)*dwdr(1:3)*weight_rji &
                                       *weight_rji*anglepart/(sum_weight**2)
           pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) + rki(3)*dwdr(1:3)*weight_rji &
                                       *weight_rji*anglepart/(sum_weight**2)
         END IF

      END DO

   END SUBROUTINE somme_d

! **************************************************************************************************
! Derivativ of the third component of angular term
! **************************************************************************************************
!> \brief ...
!> \param gal21 ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
!> \param prefactor ...
!> \param cell ...
!> \param particle_set ...
!> \param nvec ...
! **************************************************************************************************
   SUBROUTINE angular_d(gal21, r_last_update_pbc, iparticle, jparticle, f_nonbond, &
                 pv_nonbond, use_virial, prefactor, cell, particle_set, nvec)
      TYPE(gal21_pot_type), POINTER                        :: gal21
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, INTENT(IN)                                :: iparticle, jparticle
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(3)                        :: nvec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'angular_d', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: count_h, iatom, index_h1, index_h2, natom
      REAL(KIND=dp)                                      :: costheta, dsumdtheta, h_max_dist, theta,&
                                                            a1, a2, a3, a4, BH
      REAL(KIND=dp), DIMENSION(3)                        :: dangular, dcostheta, rih, rih1, rih2, &
                                                            rix, rix_hat, rji, rji_hat, rjh1, rjh2

      count_h = 0
      index_h1 = 0
      index_h2 = 0
      h_max_dist = 2.1_dp ! 1.1 angstrom
      natom = SIZE(particle_set)

      DO iatom = 1, natom !Loop on every atom of the system
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              element_symbol=element_symbol)
         IF (element_symbol /= "H") CYCLE !Kepp only hydrogen
         rih(:) = pbc(r_last_update_pbc(iparticle)%r(:), r_last_update_pbc(iatom)%r(:), cell)
         IF (SQRT(DOT_PRODUCT(rih, rih)) >= h_max_dist) CYCLE !Keep only hydrogen that are bounded to the considered O
         count_h = count_h + 1
         IF (count_h == 1) THEN
            index_h1 = iatom
         ELSEIF (count_h == 2) THEN
            index_h2 = iatom
         ENDIF
      ENDDO

      ! Abort if the oxygen is not part of a water molecule (2 H)
      IF (count_h /= 2) THEN
         CALL cp_abort(__LOCATION__, &
                       " Error: Found "//cp_to_string(count_h)//" H atoms for O atom "//cp_to_string(iparticle))
      ENDIF

      a1 = gal21%a11 + gal21%a12*gal21%gcn(jparticle) + gal21%a13*gal21%gcn(jparticle)*gal21%gcn(jparticle)
      a2 = gal21%a21 + gal21%a22*gal21%gcn(jparticle) + gal21%a23*gal21%gcn(jparticle)*gal21%gcn(jparticle)
      a3 = gal21%a31 + gal21%a32*gal21%gcn(jparticle) + gal21%a33*gal21%gcn(jparticle)*gal21%gcn(jparticle)
      a4 = gal21%a41 + gal21%a42*gal21%gcn(jparticle) + gal21%a43*gal21%gcn(jparticle)*gal21%gcn(jparticle)

      rji(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(iparticle)%r(:), cell)
      rji_hat(:) = rji(:)/SQRT(DOT_PRODUCT(rji, rji)) ! hat = pure directional component of a given vector

      !dipole vector rix of the H2O molecule
      rih1(:) = pbc(r_last_update_pbc(iparticle)%r(:), r_last_update_pbc(index_h1)%r(:), cell)
      rih2(:) = pbc(r_last_update_pbc(iparticle)%r(:), r_last_update_pbc(index_h2)%r(:), cell)
      rix(:) = rih1(:) + rih2(:) ! build the dipole vector rix of the H2O molecule
      rix_hat(:) = rix(:)/SQRT(DOT_PRODUCT(rix, rix)) ! hat = pure directional component of a given vector
      costheta = DOT_PRODUCT(rix, nvec)/SQRT(DOT_PRODUCT(rix, rix)) ! Theta is the angle between the normal to the surface and the dipole
      IF (costheta < -1.0_dp) costheta = -1.0_dp
      IF (costheta > +1.0_dp) costheta = +1.0_dp
      theta = ACOS(costheta) ! Theta is the angle between the normal to the surface and the dipole

      ! Calculation of partial derivativ of the angular components
      dsumdtheta = -1.0_dp*a1*SIN(theta) - a2*2.0_dp*SIN(2.0_dp*theta) - &
                   a3*3.0_dp*SIN(3.0_dp*theta) - a4*4.0_dp*SIN(4.0_dp*theta)
      dcostheta(:) = (1.0_dp/SQRT(DOT_PRODUCT(rix, rix)))*(nvec(:) - costheta*rix_hat(:))
      dangular(:) = prefactor*dsumdtheta*(-1.0_dp/SIN(theta))*dcostheta(:)

      !Force due to the third component of the derivativ of the angular term
      f_nonbond(1:3, iparticle) = f_nonbond(1:3, iparticle) - dangular(1:3)*2.0_dp !(one per H)
      f_nonbond(1:3, index_h1) = f_nonbond(1:3, index_h1) + dangular(1:3)
      f_nonbond(1:3, index_h2) = f_nonbond(1:3, index_h2) + dangular(1:3)

      IF (use_virial) THEN
        pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) + rix(1)*dangular(1:3)
        pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) + rix(2)*dangular(1:3)
        pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) + rix(3)*dangular(1:3)
      END IF

      BH = gal21%BH1 + gal21%gcn(jparticle)*gal21%BH2

      rjh1(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(index_h1)%r(:), cell)
      f_nonbond(1:3, index_h1) = f_nonbond(1:3, index_h1) + (gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                  BH*exp(-BH*SQRT(DOT_PRODUCT(rjh1,rjh1)))*rjh1(:)/SQRT(DOT_PRODUCT(rjh1, rjh1))

      IF (use_virial) THEN
        pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) + rjh1(1)*((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                    BH*exp(-BH*SQRT(DOT_PRODUCT(rjh1,rjh1)))) &
                    *rjh1(:)/SQRT(DOT_PRODUCT(rjh1, rjh1))
        pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) + rjh1(2)*((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                    BH*exp(-BH*SQRT(DOT_PRODUCT(rjh1,rjh1))))&
                    *rjh1(:)/SQRT(DOT_PRODUCT(rjh1, rjh1))
        pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) + rjh2(3)*((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                    BH*exp(-BH*SQRT(DOT_PRODUCT(rjh1,rjh1))))&
                    *rjh1(:)/SQRT(DOT_PRODUCT(rjh1, rjh1))
      END IF
            
      rjh2(:) = pbc(r_last_update_pbc(jparticle)%r(:), r_last_update_pbc(index_h2)%r(:), cell)
      f_nonbond(1:3, index_h2) = f_nonbond(1:3, index_h2) + ((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                  BH*exp(-BH*SQRT(DOT_PRODUCT(rjh2,rjh2))))&
                  *rjh2(:)/SQRT(DOT_PRODUCT(rjh2, rjh2))

      IF (use_virial) THEN
        pv_nonbond(1,1:3)  = pv_nonbond(1,1:3) + rjh2(1)*((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                    BH*exp(-BH*SQRT(DOT_PRODUCT(rjh2,rjh2))))&
                    *rjh2(:)/SQRT(DOT_PRODUCT(rjh2, rjh2))
        pv_nonbond(2,1:3)  = pv_nonbond(2,1:3) + rjh2(2)*((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                    BH*exp(-BH*SQRT(DOT_PRODUCT(rjh2,rjh2))))&
                    *rjh2(:)/SQRT(DOT_PRODUCT(rjh2, rjh2))
        pv_nonbond(3,1:3)  = pv_nonbond(3,1:3) + rjh2(3)*((gal21%AH2*gal21%gcn(jparticle) + gal21%AH1)* &
                    BH*exp(-BH*SQRT(DOT_PRODUCT(rjh2,rjh2))))&
                    *rjh2(:)/SQRT(DOT_PRODUCT(rjh2, rjh2))
      END IF


   END SUBROUTINE angular_d

! **************************************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param potparm ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \param cell ...
!> \par History
! **************************************************************************************************
   SUBROUTINE setup_gal21_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, &
                               glob_loc_list_a, cell)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_gal21_arrays', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iend, igrp, ikind, ilist, &
                                                            ipair, istart, jkind, nkinds, npairs, &
                                                            npairs_tot
      INTEGER, DIMENSION(:), POINTER                     :: work_list, work_list2
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rwork_list
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      CPASSERT(.NOT. ASSOCIATED(glob_loc_list))
      CPASSERT(.NOT. ASSOCIATED(glob_loc_list_a))
      CPASSERT(.NOT. ASSOCIATED(glob_cell_v))
      CALL timeset(routineN, handle)
      npairs_tot = 0
      nkinds = SIZE(potparm%pot, 1)
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Group_Loop1: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            pot => potparm%pot(ikind, jkind)%pot
            npairs = iend - istart + 1
            IF (pot%no_mb) CYCLE
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == gal21_type) npairs_tot = npairs_tot + npairs
            END DO
         END DO Kind_Group_Loop1
      END DO
      ALLOCATE (work_list(npairs_tot))
      ALLOCATE (work_list2(npairs_tot))
      ALLOCATE (glob_loc_list(2, npairs_tot))
      ALLOCATE (glob_cell_v(3, npairs_tot))
      ! Fill arrays with data
      npairs_tot = 0
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Group_Loop2: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            list => neighbor_kind_pair%list
            cvi = neighbor_kind_pair%cell_vector
            pot => potparm%pot(ikind, jkind)%pot
            npairs = iend - istart + 1
            IF (pot%no_mb) CYCLE
            CALL matvec_3x3(cell_v, cell%hmat, cvi)
            DO i = 1, SIZE(pot%type)
               ! gal21
               IF (pot%type(i) == gal21_type) THEN
                  DO ipair = 1, npairs
                     glob_loc_list(:, npairs_tot + ipair) = list(:, istart - 1 + ipair)
                     glob_cell_v(1:3, npairs_tot + ipair) = cell_v(1:3)
                  END DO
                  npairs_tot = npairs_tot + npairs
               END IF
            END DO
         END DO Kind_Group_Loop2
      END DO
      ! Order the arrays w.r.t. the first index of glob_loc_list
      CALL sort(glob_loc_list(1, :), npairs_tot, work_list)
      DO ipair = 1, npairs_tot
         work_list2(ipair) = glob_loc_list(2, work_list(ipair))
      END DO
      glob_loc_list(2, :) = work_list2
      DEALLOCATE (work_list2)
      ALLOCATE (rwork_list(3, npairs_tot))
      DO ipair = 1, npairs_tot
         rwork_list(:, ipair) = glob_cell_v(:, work_list(ipair))
      END DO
      glob_cell_v = rwork_list
      DEALLOCATE (rwork_list)
      DEALLOCATE (work_list)
      ALLOCATE (glob_loc_list_a(npairs_tot))
      glob_loc_list_a = glob_loc_list(1, :)
      CALL timestop(handle)
   END SUBROUTINE setup_gal21_arrays

! **************************************************************************************************
!> \brief ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
! **************************************************************************************************
   SUBROUTINE destroy_gal21_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a)
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a

      CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_gal21_arrays', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(glob_loc_list)) THEN
         DEALLOCATE (glob_loc_list)
      END IF
      IF (ASSOCIATED(glob_loc_list_a)) THEN
         DEALLOCATE (glob_loc_list_a)
      END IF
      IF (ASSOCIATED(glob_cell_v)) THEN
         DEALLOCATE (glob_cell_v)
      END IF

   END SUBROUTINE destroy_gal21_arrays

! **************************************************************************************************
!> \brief prints the number of OH- ions or H3O+ ions near surface
!> \param nr_ions number of ions
!> \param mm_section ...
!> \param para_env ...
!> \param print_oh flag indicating if number OH- is printed
!> \param print_h3o flag indicating if number H3O+ is printed
!> \param print_o flag indicating if number O^(2-) is printed
! **************************************************************************************************
   SUBROUTINE print_nr_ions_gal21(nr_ions, mm_section, para_env, print_oh, &
                                print_h3o, print_o)
      INTEGER, INTENT(INOUT)                             :: nr_ions
      TYPE(section_vals_type), POINTER                   :: mm_section
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env
      LOGICAL, INTENT(IN)                                :: print_oh, print_h3o, print_o

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_nr_ions_gal21', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iw
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)

      CALL mp_sum(nr_ions, para_env%group)
      logger => cp_get_default_logger()

      iw = cp_print_key_unit_nr(logger, mm_section, "PRINT%PROGRAM_RUN_INFO", &
                                extension=".mmLog")

      IF (iw > 0 .AND. nr_ions > 0 .AND. print_oh) THEN
         WRITE (iw, '(/,A,T71,I10,/)') " gal21: number of OH- ions at surface", nr_ions
      ENDIF
      IF (iw > 0 .AND. nr_ions > 0 .AND. print_h3o) THEN
         WRITE (iw, '(/,A,T71,I10,/)') " gal21: number of H3O+ ions at surface", nr_ions
      ENDIF
      IF (iw > 0 .AND. nr_ions > 0 .AND. print_o) THEN
         WRITE (iw, '(/,A,T71,I10,/)') " gal21: number of O^2- ions at surface", nr_ions
      ENDIF

      CALL cp_print_key_finished_output(iw, logger, mm_section, "PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE print_nr_ions_gal21

END MODULE manybody_gal21
