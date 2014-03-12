!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods that handle helium-solvent and helium-helium interactions
!> \author Lukasz Walewski
!> \date   2009-06-10
! *****************************************************************************
MODULE helium_interactions

  USE helium_common,                   ONLY: helium_eval_expansion,&
                                             helium_pbc
  USE helium_types,                    ONLY: e_id_interact,&
                                             e_id_kinetic,&
                                             e_id_potential,&
                                             e_id_thermo,&
                                             e_id_total,&
                                             e_id_virial,&
                                             helium_solvent_type
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE pint_types,                      ONLY: pint_env_type
  USE splines_types,                   ONLY: spline_data_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_interactions'

  PUBLIC :: helium_calc_energy
  PUBLIC :: helium_solute_e_f
  PUBLIC :: helium_bead_solute_e
  PUBLIC :: helium_solute_e

  CONTAINS

  ! ***************************************************************************
  !> \brief  Calculate the helium energy (including helium-solute interaction)
  !> \param    pint_env - path integral environment
  !> \param    helium   - helium environment
  !> \author hforbert
  !> \par History
  !>         2009-06 moved I/O out from here [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_calc_energy(helium,pint_env)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, bead, c, i, j, n
    INTEGER, DIMENSION(:), POINTER           :: perm
    LOGICAL                                  :: nperiodic
    REAL(KIND=dp)                            :: a, cell_size, en, interac, &
                                                kin, pot, rmax, rmin, vkin
    REAL(KIND=dp), DIMENSION(3)              :: r, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos
    TYPE(spline_data_p_type), &
      DIMENSION(:, :), POINTER               :: eij

   pos => helium%pos
   perm => helium%permutation
   eij => helium%eij
   cell_size = 0.5_dp*helium%cell_size
   nperiodic = .NOT.helium%periodic
   n = helium%atoms
   b = helium%beads
   en = 0.0_dp
   pot = 0.0_dp
   rmin = 1.0e20_dp
   rmax = 0.0_dp
   DO i = 1, n-1
      DO j = i+1, n
         rp(:) = pos(:,i,1) - pos(:,j,1)
         CALL helium_pbc( helium, rp )
         DO bead = 2, b
            a = 0.0_dp
            DO c = 1, 3
               r(c) = rp(c)
               a = a + r(c)**2
               rp(c) = pos(c,i,bead) - pos(c,j,bead)
            END DO
            CALL helium_pbc( helium, rp )
            en = en+helium_eval_expansion(helium,r,rp,eij,0)
            a = SQRT(a)
            IF (a < rmin) rmin = a
            IF (a > rmax) rmax = a
            IF ((a < cell_size).OR.nperiodic) THEN
               pot = pot + helium_vij(a)
            END IF
         END DO
         a = 0.0_dp
         DO c = 1, 3
            r(c) = rp(c)
            a = a + r(c)**2
            rp(c) = pos(c,perm(i),1) - pos(c,perm(j),1)
         END DO
         CALL helium_pbc( helium, rp )
         en = en + helium_eval_expansion(helium,r,rp,eij,0)
         a = SQRT(a)
         IF (a < rmin) rmin = a
         IF (a > rmax) rmax = a
         IF ((a < cell_size).OR.nperiodic) THEN
            pot = pot + helium_vij(a)
         END IF
      END DO
   END DO
   pot = pot / b
   en = en / b

    ! helium-solute interaction energy (all beads of all particles)
    interac = 0.0_dp
    IF (helium%solute_present) THEN
      CALL helium_solute_e(pint_env, helium, interac)
    END IF
    interac = interac / b

!TODO:
vkin = 0.0_dp
!   vkin = helium_virial_energy(helium)

   kin = 0.0_dp
   DO i = 1, n
      r(:) = pos(:,i,b) - pos(:,perm(i),1)
      CALL helium_pbc( helium, r )
      kin = kin + r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      DO bead = 2, b
         r(:) = pos(:,i,bead-1) - pos(:,i,bead)
         CALL helium_pbc( helium, r )
         kin = kin + r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      END DO
   END DO
   kin = 1.5_dp*n/helium%tau - 0.5*kin/(b*helium%tau**2*helium%hb2m)

! TODO: move printing somwhere else ?
!   print *,"POT = ",(pot/n+helium%e_corr)*kelvin,"K"
!   print *,"INTERAC = ",interac*kelvin,"K"
!   print *,"RMIN= ",rmin*angstrom,"A"
!   print *,"RMAX= ",rmax*angstrom,"A"
!   print *,"EVIRIAL not valid!"
!   print *,"ETHERMO= ",((en+kin)/n+helium%e_corr)*kelvin,"K"
!   print *,"ECORR= ",helium%e_corr*kelvin,"K"
!!   kin = helium_total_action(helium)
!!   print *,"ACTION= ",kin
!   print *,"WINDING#= ",helium_calc_winding(helium)

   helium%energy_inst(e_id_potential) = pot/n+helium%e_corr
   helium%energy_inst(e_id_kinetic) = (en-pot+kin)/n
   helium%energy_inst(e_id_interact) = interac
   helium%energy_inst(e_id_thermo) = (en+kin)/n+helium%e_corr
   helium%energy_inst(e_id_virial) = 0.0_dp !(en+vkin)/n+helium%e_corr
   helium%energy_inst(e_id_total) = (en+vkin)/n+helium%e_corr

   RETURN
  END SUBROUTINE helium_calc_energy

  ! ***************************************************************************
  !> \descr  Calculate helium-solute interaction energy and forces contribution
  !>         from one helium bead of one helium particle and the corresponding
  !>         solute time slice.
  !> \note   This is a stub, Lennard-Jones potential with no warranty. It should
  !>         be replaced with realistic potential for real-world simulations.
  !> \author Lukasz Walewski
  !> \date   2009-12-11
  ! ***************************************************************************
  SUBROUTINE helium_bead_solute_e_f(pint_env, helium, helium_part_index, &
    helium_slice_index, helium_r_opt, energy, force)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: helium_part_index, &
                                                helium_slice_index
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(IN), OPTIONAL                   :: helium_r_opt
    REAL(KIND=dp), INTENT(OUT)               :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: force

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_bead_solute_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: hbeads, hi, i, ig, &
                                                num_chlorine, num_hydrogen, &
                                                num_oxygen, qi
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d, d2, dd, ep, eps, s1, s2, &
                                                sig
    REAL(KIND=dp), DIMENSION(3)              :: dr, helium_r, solute_r
    TYPE(cp_error_type)                      :: error

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    num_chlorine = helium%solute_number(1)
    num_oxygen   = helium%solute_number(2)
    num_hydrogen = helium%solute_number(3)

    hbeads = helium%beads
    ! helium bead index that is invariant wrt the rotations
    hi = MOD(helium_slice_index-1+hbeads+helium%relrot,hbeads) + 1
    ! solute bead index that belongs to hi helium index
    qi = ((hi-1)*pint_env%p)/hbeads+1

    ! coordinates of the helium bead
    IF (PRESENT(helium_r_opt)) THEN
      helium_r(:) = helium_r_opt(:)
    ELSE
      helium_r(:) = helium%pos(:,helium_part_index,helium_slice_index)
    END IF

    energy = 0.0_dp
    force(:,:) = 0.0_dp

    sig = 2.69_dp     ! 1.4 Angstrom
    eps = 60.61e-6_dp ! 19 K
    s1 = 0.0_dp
    DO i = 1, num_hydrogen
      ig = helium%solute_index(3,i)-1 ! global hydrogen index (3 == H)
      solute_r(1) = pint_env%x(qi,3*ig+1)
      solute_r(2) = pint_env%x(qi,3*ig+2)
      solute_r(3) = pint_env%x(qi,3*ig+3)
      dr(:) = solute_r(:) - helium_r(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
      s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
      force(qi,3*ig+1) = force(qi,3*ig+1) + s2*dr(1)
      force(qi,3*ig+2) = force(qi,3*ig+2) + s2*dr(2)
      force(qi,3*ig+3) = force(qi,3*ig+3) + s2*dr(3)
    END DO ! i = 1, num_hydrogen
    energy = energy + s1

    sig = 5.01_dp     ! 2.6 Angstrom
    eps = 79.17e-6_dp ! 25 K
    s1 = 0.0_dp
    DO i = 1, num_chlorine
      ig = helium%solute_index(1,i)-1 ! global chlorine index (1 == Cl)
      solute_r(1) = pint_env%x(qi,3*ig+1)
      solute_r(2) = pint_env%x(qi,3*ig+2)
      solute_r(3) = pint_env%x(qi,3*ig+3)
      dr(:) = solute_r(:) - helium_r(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
      s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
      force(qi,3*ig+1) = force(qi,3*ig+1) + s2*dr(1)
      force(qi,3*ig+2) = force(qi,3*ig+2) + s2*dr(2)
      force(qi,3*ig+3) = force(qi,3*ig+3) + s2*dr(3)
    END DO ! i = 1, num_chlorine
    energy = energy + s1

    sig = 5.01_dp     ! 2.6 Angstrom
    eps = 104.5e-6_dp ! 33 K
    s1 = 0.0_dp
    DO i = 1, num_oxygen
      ig = helium%solute_index(2,i)-1 ! global oxygen index (2 == O)
      solute_r(1) = pint_env%x(qi,3*ig+1)
      solute_r(2) = pint_env%x(qi,3*ig+2)
      solute_r(3) = pint_env%x(qi,3*ig+3)
      dr(:) = solute_r(:) - helium_r(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
      s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
      force(qi,3*ig+1) = force(qi,3*ig+1) + s2*dr(1)
      force(qi,3*ig+2) = force(qi,3*ig+2) + s2*dr(2)
      force(qi,3*ig+3) = force(qi,3*ig+3) + s2*dr(3)
    END DO ! i = 1, num_chlorine
    energy = energy + s1

    RETURN
  END SUBROUTINE helium_bead_solute_e_f

  ! ***************************************************************************
  !> \brief Calculate total helium-solute interaction energy and forces.
  !> \param   pint_env - path integral environment
  !> \param   energy   - calculated interaction energy
  !> \param   force    - calculated forces (P,NDIM)
  !> \author Lukasz Walewski
  ! ***************************************************************************
  SUBROUTINE helium_solute_e_f(pint_env, helium, energy)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), INTENT(OUT)               :: energy

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_solute_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibead, ipart, j
    REAL(KIND=dp)                            :: my_energy
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: force

    NULLIFY(force)
    force => helium%force_inst

    energy = 0.0_dp
    force(:,:) = 0.0_dp

    ! calculate the total interaction energy and gradients between the
    ! solute and the helium, sum over all beads of all He particles
    DO ipart = 1, helium%atoms
      DO ibead = 1, helium%beads
        CALL helium_bead_solute_e_f(pint_env, helium, ipart, &
          ibead, energy=my_energy, force=helium%rtmp_p_ndim_2d)
        energy = energy + my_energy
        DO i = 1, pint_env%p
          DO j = 1, pint_env%ndim
            force(i,j) = force(i,j) + helium%rtmp_p_ndim_2d(i,j)
          END DO
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE helium_solute_e_f

  ! ***************************************************************************
  !> \descr  Calculate helium-solute interaction energy contribution from one
  !>         helium bead of one helium particle and the corresponding solute
  !>         time slice.
  !> \note   This is a stub, Lennard-Jones potential with no warranty. It should
  !>         be replaced with realistic potential for real-world simulations.
  !> \author Lukasz Walewski
  !> \date   2009-12-11
  ! ***************************************************************************
  SUBROUTINE helium_bead_solute_e(pint_env, helium, helium_part_index, &
    helium_slice_index, helium_r_opt, energy)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: helium_part_index, &
                                                helium_slice_index
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(IN), OPTIONAL                   :: helium_r_opt
    REAL(KIND=dp), INTENT(OUT)               :: energy

    INTEGER                                  :: hbeads, hi, i, ig, &
                                                num_chlorine, num_hydrogen, &
                                                num_oxygen, qi
    REAL(KIND=dp)                            :: d, d2, dd, ep, eps, s1, sig
    REAL(KIND=dp), DIMENSION(3)              :: dr, helium_r, solute_r

    num_chlorine = helium%solute_number(1)
    num_oxygen   = helium%solute_number(2)
    num_hydrogen = helium%solute_number(3)

    hbeads = helium%beads
    ! helium bead index that is invariant wrt the rotations
    hi = MOD(helium_slice_index-1+hbeads+helium%relrot,hbeads) + 1
    ! solute bead index that belongs to hi helium index
    qi = ((hi-1)*pint_env%p)/hbeads+1

    ! coordinates of the helium bead
    IF (PRESENT(helium_r_opt)) THEN
      helium_r(:) = helium_r_opt(:)
    ELSE
      helium_r(:) = helium%pos(:,helium_part_index,helium_slice_index)
    END IF

    energy = 0.0_dp

    sig = 2.69_dp     ! 1.4 Angstrom
    eps = 60.61e-6_dp ! 19 K
    s1 = 0.0_dp
    DO i = 1, num_hydrogen
      ig = helium%solute_index(3,i)-1 ! global hydrogen index (3 == H)
      solute_r(1) = pint_env%x(qi,3*ig+1)
      solute_r(2) = pint_env%x(qi,3*ig+2)
      solute_r(3) = pint_env%x(qi,3*ig+3)
      dr(:) = solute_r(:) - helium_r(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
    END DO ! i = 1, num_hydrogen
    energy = energy + s1

    sig = 5.01_dp     ! 2.6 Angstrom
    eps = 79.17e-6_dp ! 25 K
    s1 = 0.0_dp
    DO i = 1, num_chlorine
      ig = helium%solute_index(1,i)-1 ! global chlorine index (1 == Cl)
      solute_r(1) = pint_env%x(qi,3*ig+1)
      solute_r(2) = pint_env%x(qi,3*ig+2)
      solute_r(3) = pint_env%x(qi,3*ig+3)
      dr(:) = solute_r(:) - helium_r(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
    END DO ! i = 1, num_chlorine
    energy = energy + s1

    sig = 5.01_dp     ! 2.6 Angstrom
    eps = 104.5e-6_dp ! 33 K
    s1 = 0.0_dp
    DO i = 1, num_oxygen
      ig = helium%solute_index(2,i)-1 ! global oxygen index (2 == O)
      solute_r(1) = pint_env%x(qi,3*ig+1)
      solute_r(2) = pint_env%x(qi,3*ig+2)
      solute_r(3) = pint_env%x(qi,3*ig+3)
      dr(:) = solute_r(:) - helium_r(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
    END DO ! i = 1, num_oxygen
    energy = energy + s1

    RETURN
  END SUBROUTINE helium_bead_solute_e

  ! ***************************************************************************
  !> \brief Calculate total helium-solute interaction energy.
  !> \param   pint_env - path integral environment
  !> \param   energy   - calculated interaction energy
  !> \author Lukasz Walewski
  ! ***************************************************************************
  SUBROUTINE helium_solute_e(pint_env, helium, energy)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), INTENT(OUT)               :: energy

    INTEGER                                  :: ibead, ipart
    REAL(KIND=dp)                            :: my_energy

    energy = 0.0_dp

    DO ipart = 1, helium%atoms
      DO ibead = 1, helium%beads
        CALL helium_bead_solute_e(pint_env, helium, ipart, &
          ibead, energy=my_energy)
        energy = energy + my_energy
      END DO
    END DO

    RETURN
  END SUBROUTINE helium_solute_e

  ! ***************************************************************************
  !> \brief Helium-helium pair interaction potential.
  ! ***************************************************************************
  FUNCTION helium_vij(r) RESULT(vij)

    REAL(kind=dp), INTENT(IN)                :: r
    REAL(kind=dp)                            :: vij

    REAL(kind=dp)                            :: f, x, x2

    x = angstrom*r/2.9673_dp
    IF (x < 1.241314_dp) THEN
      x2 = 1.241314_dp/x-1.0_dp
      f = EXP(-x2*x2)
    ELSE
      f = 1.0_dp
    END IF
    x2 = 1.0_dp/(x*x)
    vij = 10.8_dp/kelvin*(544850.4_dp*EXP(-13.353384_dp*x)-f* &
          ((0.1781_dp*x2+0.4253785_dp)*x2+1.3732412_dp)*x2*x2*x2)
    RETURN
  END FUNCTION helium_vij

  ! ***************************************************************************
  !> \brief Helium-helium pair interaction potential's derivative.
  ! ***************************************************************************
  FUNCTION helium_d_vij(r) RESULT(dvij)

    REAL(kind=dp), INTENT(IN)                :: r
    REAL(kind=dp)                            :: dvij

    REAL(kind=dp)                            :: f, fp, x, x2, y

    x = angstrom*r/2.9673_dp
    x = r/2.9673_dp
    x2 = 1.0_dp/(x*x)
    IF (x < 1.241314_dp) THEN
      y = 1.241314_dp/x-1.0_dp
      f = EXP(-y*y)
      fp= 2.0_dp*1.241314_dp*f*y* &
          ((0.1781_dp*x2+0.4253785_dp)*x2+1.3732412_dp)*x2*x2*x2*x2
    ELSE
      f = 1.0_dp
      fp = 0.0_dp
    END IF

    dvij = angstrom*(10.8_dp/2.9673_dp)*((-13.353384_dp*544850.4_dp)* &
           EXP(-13.353384_dp*x)-fp+f*(((10.0_dp*0.1781_dp)*x2+ &
           (8.0_dp*0.4253785_dp))*x2+(6.0_dp*1.3732412_dp))*x2*x2*x2/x)/( &
           r*kelvin )
    RETURN
  END FUNCTION helium_d_vij

END MODULE helium_interactions
