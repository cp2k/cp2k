!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods that handle helium-solvent and helium-helium interactions
!> \author Lukasz Walewski
!> \date   2009-06-10
! *****************************************************************************
MODULE helium_interactions

  USE helium_common,                   ONLY: helium_eval_expansion,&
                                             helium_pbc
  USE helium_types,                    ONLY: &
       e_id_interact, e_id_kinetic, e_id_potential, e_id_thermo, e_id_total, &
       e_id_virial, helium_solvent_type, hid_chlorine, hid_hydrogen, &
       hid_oxygen, int_arr_ptr
  USE input_constants,                 ONLY: helium_solute_intpot_mwater,&
                                             helium_solute_intpot_none
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE pint_types,                      ONLY: pint_env_type
  USE splines_types,                   ONLY: spline_data_p_type
#include "../base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_interactions'

  PUBLIC :: helium_calc_energy
  PUBLIC :: helium_solute_e_f
  PUBLIC :: helium_bead_solute_e_f
  PUBLIC :: helium_intpot_scan

  CONTAINS

! ***************************************************************************
!> \brief  Calculate the helium energy (including helium-solute interaction)
!> \param    helium   - helium environment
!> \param    pint_env - path integral environment
!> \par History
!>         2009-06 moved I/O out from here [lwalewski]
!> \author hforbert
! *****************************************************************************
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
!> \brief Calculate general helium-solute interaction energy (and forces)
!>        between one helium bead and the corresponding solute time slice.
!> \param pint_env           path integral environment
!> \param helium ...
!> \param helium_part_index  helium particle index
!> \param helium_slice_index helium time slice index
!> \param helium_r_opt       explicit helium bead coordinates (optional)
!> \param energy             calculated energy
!> \param force              calculated force (if requested)
!> \author Lukasz Walewski
! *****************************************************************************
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
      INTENT(OUT), OPTIONAL, POINTER         :: force

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_bead_solute_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: hbeads, hi, qi
    REAL(KIND=dp), DIMENSION(3)              :: helium_r
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_force

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

    SELECT CASE (helium%solute_interaction)


      CASE (helium_solute_intpot_mwater)
        IF (PRESENT(force)) THEN
          force(:,:) = 0.0_dp
          my_force => force(qi,:)
          CALL helium_intpot_model_water( &
            pint_env%x(qi,:), &
            helium%solute_i, &
            helium, &
            helium_r,&
            energy, &
            my_force &
          )
        ELSE
          CALL helium_intpot_model_water( &
          pint_env%x(qi,:), &
          helium%solute_i, &
          helium, &
          helium_r,&
          energy &
          )
        END IF

      CASE (helium_solute_intpot_none)
        energy = 0.0_dp
        IF (PRESENT(force)) THEN
          force(:,:) = 0.0_dp
        END IF



      CASE DEFAULT

    END SELECT

    RETURN
  END SUBROUTINE helium_bead_solute_e_f


! ***************************************************************************
!> \brief Calculate total helium-solute interaction energy and forces.
!> \param   pint_env - path integral environment
!> \param helium ...
!> \param   energy   - calculated interaction energy
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_solute_e_f(pint_env, helium, energy)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), INTENT(OUT)               :: energy

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_solute_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib, jb, jc
    REAL(KIND=dp)                            :: my_energy
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: force

    NULLIFY(force)
    force => helium%force_inst

    energy = 0.0_dp
    force(:,:) = 0.0_dp

    ! calculate the total interaction energy and gradients between the
    ! solute and the helium, sum over all beads of all He particles
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        CALL helium_bead_solute_e_f(pint_env, helium, ia, ib, &
          energy=my_energy, force=helium%rtmp_p_ndim_2d)
          energy = energy + my_energy
        DO jb = 1, pint_env%p
          DO jc = 1, pint_env%ndim
            force(jb,jc) = force(jb,jc) + helium%rtmp_p_ndim_2d(jb,jc)
          END DO
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE helium_solute_e_f


! ***************************************************************************
!> \brief Calculate total helium-solute interaction energy.
!> \param   pint_env - path integral environment
!> \param helium ...
!> \param   energy   - calculated interaction energy
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_solute_e(pint_env, helium, energy)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), INTENT(OUT)               :: energy

    INTEGER                                  :: ia, ib
    REAL(KIND=dp)                            :: my_energy

    energy = 0.0_dp

    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        CALL helium_bead_solute_e_f(pint_env, helium, ia, ib, energy=my_energy)
        energy = energy + my_energy
      END DO
    END DO

    RETURN
  END SUBROUTINE helium_solute_e


! ***************************************************************************
!> \brief Calculate l-th Legendre polynomial P_{l}(x) at point x
!> \param x ...
!> \param n ...
!> \retval Pl ...
! *****************************************************************************
  FUNCTION Pl(x,n)
    REAL(KIND=dp), INTENT(IN)                :: x
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp)                            :: Pl

    INTEGER                                  :: k
    REAL(KIND=dp)                            :: pln(0:n)

    pln(0) = 1.0_dp
    pln(1) = x

    IF (n <= 1) THEN
      Pl = pln(n)
    ELSE
      DO k=1,n-1
        pln(k+1) = ((2.0*k+1.0)*x*pln(k) - REAL(k,dp)*pln(k-1))/(REAL(k+1,dp))
      END DO
      Pl = pln(n)
    END IF
    RETURN
  END FUNCTION Pl


! ***************************************************************************
!> \brief  Scan the helium-solute interaction energy within the periodic cell
!> \param pint_env ...
!> \param helium ...
!> \date   2014-01-22
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_intpot_scan(pint_env, helium)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_intpot_scan', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ic, ix, iy, iz, nbin
    LOGICAL                                  :: wrapped
    REAL(KIND=dp)                            :: delr, my_en, ox, oy, oz
    REAL(kind=dp), DIMENSION(3)              :: pbc1, pbc2, pos

    CALL timeset(routineN,handle)

    nbin = helium%rho_nbin
    delr = helium%rho_delr
    ox = helium%center(1) - helium%rho_maxr / 2.0_dp
    oy = helium%center(2) - helium%rho_maxr / 2.0_dp
    oz = helium%center(3) - helium%rho_maxr / 2.0_dp

    DO ix = 1, nbin
      DO iy = 1, nbin
        DO iz = 1, nbin

          ! put the probe in the center of the current voxel
          pos(:) = (/ox+(ix-0.5_dp)*delr, oy+(iy-0.5_dp)*delr, oz+(iz-0.5_dp)*delr/)

          ! calc interaction energy for the current probe position
          helium%pos(:,1,1) = pos(:)
          CALL helium_bead_solute_e_f(pint_env,helium,1,1,energy=my_en)

          ! check if the probe fits within the unit cell
          pbc1(:) = pos(:) - helium%center
          pbc2(:) = pbc1(:)
          CALL helium_pbc( helium, pbc2 )
          wrapped = .FALSE.
          DO ic = 1, 3
            IF (ABS(pbc1(ic)-pbc2(ic)) .GT. 10.0_dp * EPSILON(0.0_dp)) THEN
              wrapped = .TRUE.
            END IF
          END DO

          ! set the interaction energy value
          IF (wrapped) THEN
            helium%rho_inst(1,ix,iy,iz) = 0.0_dp
          ELSE
            helium%rho_inst(1,ix,iy,iz) = my_en
          END IF

        END DO
      END DO
    END DO

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_intpot_scan


! ***************************************************************************
!> \brief Calculate model helium-solute interaction energy and forces 
!>        between one helium bead and the corresponding solute time 
!>        slice asuming water solute.
!> \param solute_x  solute positions ARR(3*NATOMS)
!> \param solute_i  solute indices, mapping from instances of a given element
!>        to global atom indices
!> \param helium    only needed for helium_pbc call at the moment
!> \param helium_x  helium bead position ARR(3)
!> \param energy    calculated interaction energy
!> \param force ...
!> \author Felix Uhl
! *****************************************************************************
  SUBROUTINE helium_intpot_model_water(solute_x, solute_i, helium, helium_x, energy, force)


    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: solute_x
    TYPE(int_arr_ptr), DIMENSION(:), &
      INTENT(IN)                             :: solute_i
    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: helium_x
    REAL(KIND=dp), INTENT(OUT)               :: energy
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(OUT), OPTIONAL, POINTER         :: force

    CHARACTER(LEN=*), PARAMETER :: routineN = 'helium_intpot_model_water', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ig, num_chlorine, &
                                                num_hydrogen, num_oxygen
    REAL(KIND=dp)                            :: d, d2, dd, ep, eps, s1, s2, &
                                                sig
    REAL(KIND=dp), DIMENSION(3)              :: dr, solute_r

      IF (ASSOCIATED(solute_i(hid_chlorine)%iap)) THEN
        num_chlorine = SIZE(solute_i(hid_chlorine)%iap)
      ELSE
        num_chlorine = 0
      END IF

      IF (ASSOCIATED(solute_i(hid_oxygen)%iap)) THEN
        num_oxygen = SIZE(solute_i(hid_oxygen)%iap)
      ELSE
        num_oxygen = 0
      END IF

      IF (ASSOCIATED(solute_i(hid_hydrogen)%iap)) THEN
        num_hydrogen = SIZE(solute_i(hid_hydrogen)%iap)
      ELSE
        num_hydrogen = 0
      END IF
       
      energy = 0.0_dp
      IF (PRESENT(force)) THEN
          force(:) = 0.0_dp
      END IF

      sig = 2.69_dp     ! 1.4 Angstrom
      eps = 60.61e-6_dp ! 19 K
      s1 = 0.0_dp
      DO i = 1, num_hydrogen
      ig = solute_i(hid_hydrogen)%iap(i)-1 ! global hydrogen index (3 == H)
      solute_r(1) = solute_x(3*ig+1)
      solute_r(2) = solute_x(3*ig+2)
      solute_r(3) = solute_x(3*ig+3)
      dr(:) = solute_r(:) - helium_x(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
      s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
      IF (PRESENT(force)) THEN
          force(3*ig+1) = force(3*ig+1) + s2*dr(1)
          force(3*ig+2) = force(3*ig+2) + s2*dr(2)
          force(3*ig+3) = force(3*ig+3) + s2*dr(3)
      END IF
      END DO ! i = 1, num_hydrogen
      energy = energy + s1

      sig = 5.01_dp     ! 2.6 Angstrom
      eps = 104.5e-6_dp ! 33 K
      s1 = 0.0_dp
      DO i = 1, num_oxygen
      ig = solute_i(hid_oxygen)%iap(i)-1 ! global oxygen index (2 == O)
      solute_r(1) = solute_x(3*ig+1)
      solute_r(2) = solute_x(3*ig+2)
      solute_r(3) = solute_x(3*ig+3)
      dr(:) = solute_r(:) - helium_x(:)
      CALL helium_pbc( helium, dr )
      d2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
      d = SQRT(d2)
      dd = (sig/d)**6
      ep = 4.0_dp*eps*dd*(dd-1.0_dp)
      s1 = s1 + ep
      s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
      IF (PRESENT(force)) THEN
          force(3*ig+1) = force(3*ig+1) + s2*dr(1)
          force(3*ig+2) = force(3*ig+2) + s2*dr(2)
          force(3*ig+3) = force(3*ig+3) + s2*dr(3)
      END IF
      END DO ! i = 1, num_chlorine
      energy = energy + s1

      RETURN

END SUBROUTINE helium_intpot_model_water



! ***************************************************************************
!> \brief Helium-helium pair interaction potential.
!> \param r ...
!> \retval vij ...
! *****************************************************************************
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


#if 0

  ! this block is currently turned off

! ***************************************************************************
!> \brief Helium-helium pair interaction potential's derivative.
!> \param r ...
!> \retval dvij ...
! *****************************************************************************
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


! *****************************************************************************
!> \brief ...
!> \param helium ...
!> \param n ...
!> \param i ...
!> \retval res ...
! *****************************************************************************
  FUNCTION helium_atom_action(helium,n,i) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: n, i
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: c, j
    REAL(KIND=dp)                            :: r(3), rp(3), s, t

    s = 0.0_dp
    t = 0.0_dp
    IF (n < helium%beads) THEN
      DO c = 1, 3
        r(c) = helium%pos(c,i,n) - helium%pos(c,i,n+1)
      END DO
      CALL helium_pbc( helium, r )
      t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      DO j = 1, i-1
        DO c = 1, 3
          r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
          rp(c) = helium%pos(c,i,n+1) - helium%pos(c,j,n+1)
        END DO
        CALL helium_pbc( helium, r )
        CALL helium_pbc( helium, rp )
        s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
      DO j = i+1, helium%atoms
        DO c = 1, 3
          r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
          rp(c) = helium%pos(c,i,n+1) - helium%pos(c,j,n+1)
        END DO
        CALL helium_pbc( helium, r )
        CALL helium_pbc( helium, rp )
        s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
    ELSE
      DO c = 1, 3
        r(c) = helium%pos(c,i,n) - helium%pos(c,helium%permutation(i),1)
      END DO
      CALL helium_pbc( helium, r )
      t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      DO j = 1, i-1
        DO c = 1, 3
          r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
          rp(c) = helium%pos(c,helium%permutation(i),1) - helium%pos(c,helium%permutation(j),1)
        END DO
        CALL helium_pbc( helium, r )
        CALL helium_pbc( helium, rp )
        s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
      DO j = i+1, helium%atoms
        DO c = 1, 3
          r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
          rp(c) = helium%pos(c,helium%permutation(i),1) - helium%pos(c,helium%permutation(j),1)
        END DO
        CALL helium_pbc( helium, r )
        CALL helium_pbc( helium, rp )
        s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
    END IF
    t = t / (2.0_dp*helium%tau*helium%hb2m)
    s = s * 0.5_dp
    res = s+t
    RETURN

  END FUNCTION helium_atom_action

! *****************************************************************************
!> \brief ...
!> \param helium ...
!> \param n ...
!> \retval res ...
! *****************************************************************************
  FUNCTION helium_link_action(helium,n) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: c, i, j
    REAL(KIND=dp)                            :: r(3), rp(3), s, t

    s = 0.0_dp
    t = 0.0_dp
    IF (n < helium%beads) THEN
      DO i = 1, helium%atoms
        DO c = 1, 3
          r(c) = helium%pos(c,i,n) - helium%pos(c,i,n+1)
        END DO
        CALL helium_pbc( helium, r )
        t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
        DO j = 1, i-1
          DO c = 1, 3
            r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
            rp(c) = helium%pos(c,i,n+1) - helium%pos(c,j,n+1)
          END DO
          CALL helium_pbc( helium, r )
          CALL helium_pbc( helium, rp )
          s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
        END DO
      END DO
    ELSE
      DO i = 1, helium%atoms
        DO c = 1, 3
          r(c) = helium%pos(c,i,n) - helium%pos(c,helium%permutation(i),1)
        END DO
        CALL helium_pbc( helium, r )
        t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
        DO j = 1, i-1
          DO c = 1, 3
            r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
            rp(c) = helium%pos(c,helium%permutation(i),1) - helium%pos(c,helium%permutation(j),1)
          END DO
          CALL helium_pbc( helium, r )
          CALL helium_pbc( helium, rp )
          s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
        END DO
      END DO
    END IF
    t = t / (2.0_dp*helium%tau*helium%hb2m)
    res = s+t
    RETURN

  END FUNCTION helium_link_action

! *****************************************************************************
!> \brief ...
!> \param helium ...
!> \retval res ...
! *****************************************************************************
  FUNCTION helium_total_action(helium) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: s

    s = 0.0_dp
    DO i = 1, helium%beads
      s = s + helium_link_action(helium,i)
    END DO
    res = s
    RETURN

  END FUNCTION helium_total_action

! *****************************************************************************
!> \brief ...
!> \param helium ...
!> \param part ...
!> \param ref_bead ...
!> \param delta_bead ...
!> \param d ...
! *****************************************************************************
  SUBROUTINE helium_delta_pos(helium,part,ref_bead,delta_bead,d)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: part, ref_bead, delta_bead
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: d

    INTEGER                                  :: b, bead, db, nbead, np, p
    REAL(KIND=dp), DIMENSION(3)              :: r

    b = helium%beads

    d(:) = 0.0_dp
    IF (delta_bead > 0) THEN
      bead = ref_bead
      p = part
      db = delta_bead
      DO
        IF (db < 1) EXIT
        nbead = bead + 1
        np = p
        IF (nbead > b) THEN
          nbead = nbead - b
          np = helium%permutation(np)
        END IF
        r(:) = helium%pos(:,p,bead) - helium%pos(:,np,nbead)
        CALL helium_pbc( helium, r )
        d(:) = d(:) + r(:)
        bead = nbead
        p = np
        db = db-1
      END DO
    ELSEIF ( delta_bead < 0) THEN
      bead = ref_bead
      p = part
      db = delta_bead
      DO
        IF (db >= 0) EXIT
        nbead = bead - 1
        np = p
        IF (nbead < 1) THEN
          nbead = nbead + b
          np = helium%iperm(np)
        END IF
        r(:) = helium%pos(:,p,bead) - helium%pos(:,np,nbead)
        CALL helium_pbc( helium, r )
        d(:) = d(:) + r(:)
        bead = nbead
        p = np
        db = db + 1
      END DO
    END IF
    RETURN
  END SUBROUTINE helium_delta_pos

#endif


END MODULE helium_interactions
