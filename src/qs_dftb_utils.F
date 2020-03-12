!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Working with the DFTB parameter types.
!> \author JGH (24.02.2007)
! **************************************************************************************************
MODULE qs_dftb_utils

   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_utils'

   ! Maximum number of points used for interpolation
   INTEGER, PARAMETER                     :: max_inter = 5
   ! Maximum number of points used for extrapolation
   INTEGER, PARAMETER                     :: max_extra = 9
   ! see also qs_dftb_parameters
   REAL(dp), PARAMETER                    :: slako_d0 = 1._dp
   ! pointer to skab
   INTEGER, DIMENSION(0:3, 0:3, 0:3, 0:3, 0:3):: iptr
   ! small real number
   REAL(dp), PARAMETER                    :: rtiny = 1.e-10_dp
   ! eta(0) for mm atoms and non-scc qm atoms
   REAL(dp), PARAMETER                    :: eta_mm = 0.47_dp
   ! step size for qmmm finite difference
   REAL(dp), PARAMETER                    :: ddrmm = 0.0001_dp

   PUBLIC :: allocate_dftb_atom_param, &
             deallocate_dftb_atom_param, &
             get_dftb_atom_param, &
             set_dftb_atom_param, &
             write_dftb_atom_param
   PUBLIC :: compute_block_sk, &
             urep_egr, iptr

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dftb_parameter ...
! **************************************************************************************************
   SUBROUTINE allocate_dftb_atom_param(dftb_parameter)

      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dftb_atom_param', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(dftb_parameter)) &
         CALL deallocate_dftb_atom_param(dftb_parameter)

      ALLOCATE (dftb_parameter)

      dftb_parameter%defined = .FALSE.
      dftb_parameter%name = ""
      dftb_parameter%typ = "NONE"
      dftb_parameter%z = -1
      dftb_parameter%zeff = -1.0_dp
      dftb_parameter%natorb = 0
      dftb_parameter%lmax = -1
      dftb_parameter%skself = 0.0_dp
      dftb_parameter%occupation = 0.0_dp
      dftb_parameter%eta = 0.0_dp
      dftb_parameter%energy = 0.0_dp
      dftb_parameter%xi = 0.0_dp
      dftb_parameter%di = 0.0_dp
      dftb_parameter%rcdisp = 0.0_dp
      dftb_parameter%dudq = 0.0_dp

   END SUBROUTINE allocate_dftb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param dftb_parameter ...
! **************************************************************************************************
   SUBROUTINE deallocate_dftb_atom_param(dftb_parameter)

      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dftb_atom_param', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(dftb_parameter))
      DEALLOCATE (dftb_parameter)

   END SUBROUTINE deallocate_dftb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param dftb_parameter ...
!> \param name ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param lmax ...
!> \param skself ...
!> \param occupation ...
!> \param eta ...
!> \param energy ...
!> \param cutoff ...
!> \param xi ...
!> \param di ...
!> \param rcdisp ...
!> \param dudq ...
! **************************************************************************************************
   SUBROUTINE get_dftb_atom_param(dftb_parameter, name, typ, defined, z, zeff, natorb, &
                                  lmax, skself, occupation, eta, energy, cutoff, xi, di, rcdisp, dudq)

      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: name, typ
      LOGICAL, INTENT(OUT), OPTIONAL                     :: defined
      INTEGER, INTENT(OUT), OPTIONAL                     :: z
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: zeff
      INTEGER, INTENT(OUT), OPTIONAL                     :: natorb, lmax
      REAL(KIND=dp), DIMENSION(0:3), OPTIONAL            :: skself, occupation, eta
      REAL(KIND=dp), OPTIONAL                            :: energy, cutoff, xi, di, rcdisp, dudq

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dftb_atom_param', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(dftb_parameter))

      IF (PRESENT(name)) name = dftb_parameter%name
      IF (PRESENT(typ)) typ = dftb_parameter%typ
      IF (PRESENT(defined)) defined = dftb_parameter%defined
      IF (PRESENT(z)) z = dftb_parameter%z
      IF (PRESENT(zeff)) zeff = dftb_parameter%zeff
      IF (PRESENT(natorb)) natorb = dftb_parameter%natorb
      IF (PRESENT(lmax)) lmax = dftb_parameter%lmax
      IF (PRESENT(skself)) skself = dftb_parameter%skself
      IF (PRESENT(eta)) eta = dftb_parameter%eta
      IF (PRESENT(energy)) energy = dftb_parameter%energy
      IF (PRESENT(cutoff)) cutoff = dftb_parameter%cutoff
      IF (PRESENT(occupation)) occupation = dftb_parameter%occupation
      IF (PRESENT(xi)) xi = dftb_parameter%xi
      IF (PRESENT(di)) di = dftb_parameter%di
      IF (PRESENT(rcdisp)) rcdisp = dftb_parameter%rcdisp
      IF (PRESENT(dudq)) dudq = dftb_parameter%dudq

   END SUBROUTINE get_dftb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param dftb_parameter ...
!> \param name ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param lmax ...
!> \param skself ...
!> \param occupation ...
!> \param eta ...
!> \param energy ...
!> \param cutoff ...
!> \param xi ...
!> \param di ...
!> \param rcdisp ...
!> \param dudq ...
! **************************************************************************************************
   SUBROUTINE set_dftb_atom_param(dftb_parameter, name, typ, defined, z, zeff, natorb, &
                                  lmax, skself, occupation, eta, energy, cutoff, xi, di, rcdisp, dudq)

      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: name, typ
      LOGICAL, INTENT(IN), OPTIONAL                      :: defined
      INTEGER, INTENT(IN), OPTIONAL                      :: z
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: zeff
      INTEGER, INTENT(IN), OPTIONAL                      :: natorb, lmax
      REAL(KIND=dp), DIMENSION(0:3), OPTIONAL            :: skself, occupation, eta
      REAL(KIND=dp), OPTIONAL                            :: energy, cutoff, xi, di, rcdisp, dudq

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_dftb_atom_param', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(dftb_parameter))

      IF (PRESENT(name)) dftb_parameter%name = name
      IF (PRESENT(typ)) dftb_parameter%typ = typ
      IF (PRESENT(defined)) dftb_parameter%defined = defined
      IF (PRESENT(z)) dftb_parameter%z = z
      IF (PRESENT(zeff)) dftb_parameter%zeff = zeff
      IF (PRESENT(natorb)) dftb_parameter%natorb = natorb
      IF (PRESENT(lmax)) dftb_parameter%lmax = lmax
      IF (PRESENT(skself)) dftb_parameter%skself = skself
      IF (PRESENT(eta)) dftb_parameter%eta = eta
      IF (PRESENT(occupation)) dftb_parameter%occupation = occupation
      IF (PRESENT(energy)) dftb_parameter%energy = energy
      IF (PRESENT(cutoff)) dftb_parameter%cutoff = cutoff
      IF (PRESENT(xi)) dftb_parameter%xi = xi
      IF (PRESENT(di)) dftb_parameter%di = di
      IF (PRESENT(rcdisp)) dftb_parameter%rcdisp = rcdisp
      IF (PRESENT(dudq)) dftb_parameter%dudq = dudq

   END SUBROUTINE set_dftb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param dftb_parameter ...
!> \param subsys_section ...
! **************************************************************************************************
   SUBROUTINE write_dftb_atom_param(dftb_parameter, subsys_section)

      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter
      TYPE(section_vals_type), POINTER                   :: subsys_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_dftb_atom_param', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: name, typ
      INTEGER                                            :: lmax, natorb, output_unit, z
      LOGICAL                                            :: defined
      REAL(dp)                                           :: zeff
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (ASSOCIATED(dftb_parameter) .AND. &
          BTEST(cp_print_key_should_output(logger%iter_info, subsys_section, &
                                           "PRINT%KINDS/POTENTIAL"), cp_p_file)) THEN

         output_unit = cp_print_key_unit_nr(logger, subsys_section, "PRINT%KINDS", &
                                            extension=".Log")

         IF (output_unit > 0) THEN
            CALL get_dftb_atom_param(dftb_parameter, name=name, typ=typ, defined=defined, &
                                     z=z, zeff=zeff, natorb=natorb, lmax=lmax)

            WRITE (UNIT=output_unit, FMT="(/,A,T67,A14)") &
               " DFTB  parameters: ", TRIM(name)
            IF (defined) THEN
               WRITE (UNIT=output_unit, FMT="(T16,A,T71,F10.2)") &
                  "Effective core charge:", zeff
               WRITE (UNIT=output_unit, FMT="(T16,A,T71,I10)") &
                  "Number of orbitals:", natorb
            ELSE
               WRITE (UNIT=output_unit, FMT="(T55,A)") &
                  "Parameters are not defined"
            END IF
         END IF
         CALL cp_print_key_finished_output(output_unit, logger, subsys_section, &
                                           "PRINT%KINDS")
      END IF

   END SUBROUTINE write_dftb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param block ...
!> \param smatij ...
!> \param smatji ...
!> \param rij ...
!> \param ngrd ...
!> \param ngrdcut ...
!> \param dgrd ...
!> \param llm ...
!> \param lmaxi ...
!> \param lmaxj ...
!> \param irow ...
!> \param iatom ...
! **************************************************************************************************
   SUBROUTINE compute_block_sk(block, smatij, smatji, rij, ngrd, ngrdcut, dgrd, &
                               llm, lmaxi, lmaxj, irow, iatom)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, smatij, smatji
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      INTEGER                                            :: ngrd, ngrdcut
      REAL(KIND=dp)                                      :: dgrd
      INTEGER                                            :: llm, lmaxi, lmaxj, irow, iatom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_block_sk', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: dr
      REAL(KIND=dp), DIMENSION(20)                       :: skabij, skabji

      dr = SQRT(SUM(rij(:)**2))
      CALL getskz(smatij, skabij, dr, ngrd, ngrdcut, dgrd, llm)
      CALL getskz(smatji, skabji, dr, ngrd, ngrdcut, dgrd, llm)
      IF (irow == iatom) THEN
         CALL turnsk(block, skabji, skabij, rij, dr, lmaxi, lmaxj)
      ELSE
         CALL turnsk(block, skabij, skabji, -rij, dr, lmaxj, lmaxi)
      END IF

   END SUBROUTINE compute_block_sk

! **************************************************************************************************
!> \brief Gets matrix elements on z axis, as they are stored in the tables
!> \param slakotab ...
!> \param skpar ...
!> \param dx ...
!> \param ngrd ...
!> \param ngrdcut ...
!> \param dgrd ...
!> \param llm ...
!> \author 07. Feb. 2004, TH
! **************************************************************************************************
   SUBROUTINE getskz(slakotab, skpar, dx, ngrd, ngrdcut, dgrd, llm)
      REAL(dp), INTENT(in)                               :: slakotab(:, :), dx
      INTEGER, INTENT(in)                                :: ngrd, ngrdcut
      REAL(dp), INTENT(in)                               :: dgrd
      INTEGER, INTENT(in)                                :: llm
      REAL(dp), INTENT(out)                              :: skpar(llm)

      INTEGER                                            :: clgp

      skpar = 0._dp
      !
      ! Determine closest grid point
      !
      clgp = NINT(dx/dgrd)
      !
      ! Screen elements which are too far away
      !
      IF (clgp > ngrdcut) RETURN
      !
      ! The grid point is either contained in the table --> matrix element
      ! can be interpolated, or it is outside the table --> matrix element
      ! needs to be extrapolated.
      !
      IF (clgp > ngrd) THEN
         !
         ! Extrapolate external matrix elements if table does not finish with zero
         !
         CALL extrapol(slakotab, skpar, dx, ngrd, dgrd, llm)
      ELSE
         !
         ! Interpolate tabulated matrix elements
         !
         CALL interpol(slakotab, skpar, dx, ngrd, dgrd, llm, clgp)
      END IF
   END SUBROUTINE getskz

! **************************************************************************************************
!> \brief ...
!> \param slakotab ...
!> \param skpar ...
!> \param dx ...
!> \param ngrd ...
!> \param dgrd ...
!> \param llm ...
!> \param clgp ...
! **************************************************************************************************
   SUBROUTINE interpol(slakotab, skpar, dx, ngrd, dgrd, llm, clgp)
      REAL(dp), INTENT(in)                               :: slakotab(:, :), dx
      INTEGER, INTENT(in)                                :: ngrd
      REAL(dp), INTENT(in)                               :: dgrd
      INTEGER, INTENT(in)                                :: llm
      REAL(dp), INTENT(out)                              :: skpar(llm)
      INTEGER, INTENT(in)                                :: clgp

      INTEGER                                            :: fgpm, k, l, lgpm
      REAL(dp)                                           :: error, xa(max_inter), ya(max_inter)

      lgpm = MIN(clgp + INT(max_inter/2.0), ngrd)
      fgpm = lgpm - max_inter + 1
      DO k = 0, max_inter - 1
         xa(k + 1) = (fgpm + k)*dgrd
      END DO
      !
      ! Interpolate matrix elements for all orbitals
      !
      DO l = 1, llm
         !
         ! Read SK parameters from table
         !
         ya(1:max_inter) = slakotab(fgpm:lgpm, l)
         CALL polint(xa, ya, max_inter, dx, skpar(l), error)
      END DO
   END SUBROUTINE interpol

! **************************************************************************************************
!> \brief ...
!> \param slakotab ...
!> \param skpar ...
!> \param dx ...
!> \param ngrd ...
!> \param dgrd ...
!> \param llm ...
! **************************************************************************************************
   SUBROUTINE extrapol(slakotab, skpar, dx, ngrd, dgrd, llm)
      REAL(dp), INTENT(in)                               :: slakotab(:, :), dx
      INTEGER, INTENT(in)                                :: ngrd
      REAL(dp), INTENT(in)                               :: dgrd
      INTEGER, INTENT(in)                                :: llm
      REAL(dp), INTENT(out)                              :: skpar(llm)

      INTEGER                                            :: fgp, k, l, lgp, ntable, nzero
      REAL(dp)                                           :: error, xa(max_extra), ya(max_extra)

      nzero = max_extra/3
      ntable = max_extra - nzero
      !
      ! Get the three last distances from the table
      !
      DO k = 1, ntable
         xa(k) = (ngrd - (max_extra - 3) + k)*dgrd
      END DO
      DO k = 1, nzero
         xa(ntable + k) = (ngrd + k - 1)*dgrd + slako_d0
         ya(ntable + k) = 0.0
      END DO
      !
      ! Extrapolate matrix elements for all orbitals
      !
      DO l = 1, llm
         !
         ! Read SK parameters from table
         !
         fgp = ngrd + 1 - (max_extra - 3)
         lgp = ngrd
         ya(1:max_extra - 3) = slakotab(fgp:lgp, l)
         CALL polint(xa, ya, max_extra, dx, skpar(l), error)
      END DO
   END SUBROUTINE extrapol

! **************************************************************************************************
!> \brief   Turn matrix element from z-axis to orientation of dxv
!> \param mat ...
!> \param skab1 ...
!> \param skab2 ...
!> \param dxv ...
!> \param dx ...
!> \param lmaxa ...
!> \param lmaxb ...
!> \date    13. Jan 2004
!> \par Notes
!>          These routines are taken from an old TB code (unknown to TH).
!>          They are highly optimised and taken because they are time critical.
!>          They are explicit, so not recursive, and work up to d functions.
!>
!>          Set variables necessary for rotation of matrix elements
!>
!>          r_i^2/r, replicated in rr2(4:6) for index convenience later
!>          r_i/r, direction vector, rr(4:6) are replicated from 1:3
!>          lmax of A and B
!> \author  TH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE turnsk(mat, skab1, skab2, dxv, dx, lmaxa, lmaxb)
      REAL(dp), INTENT(inout)                  :: mat(:, :)
      REAL(dp), INTENT(in)                     :: skab1(:), skab2(:), dxv(3), dx
      INTEGER, INTENT(in)                      :: lmaxa, lmaxb

      CHARACTER(len=*), PARAMETER :: routineN = 'turnsk', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: lmaxab, minlmaxab
      REAL(dp)                                 :: rinv, rr(6), rr2(6)

      lmaxab = MAX(lmaxa, lmaxb)
      ! Determine l quantum limits.
      IF (lmaxab .GT. 2) CPABORT('lmax=2')
      minlmaxab = MIN(lmaxa, lmaxb)
      !
      ! s-s interaction
      !
      CALL skss(skab1, mat)
      !
      IF (lmaxab .LE. 0) RETURN
      !
      rr2(1:3) = dxv(1:3)**2
      rr(1:3) = dxv(1:3)
      rinv = 1.0_dp/dx
      !
      rr(1:3) = rr(1:3)*rinv
      rr(4:6) = rr(1:3)
      rr2(1:3) = rr2(1:3)*rinv**2
      rr2(4:6) = rr2(1:3)
      !
      ! s-p, p-s and p-p interaction
      !
      IF (minlmaxab .GE. 1) THEN
         CALL skpp(skab1, mat, iptr(:, :, :, lmaxa, lmaxb))
         CALL sksp(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
         CALL sksp(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
      ELSE
         IF (lmaxb .GE. 1) THEN
            CALL sksp(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
         ELSE
            CALL sksp(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
         END IF
      END IF
      !
      ! If there is only s-p interaction we have finished
      !
      IF (lmaxab .LE. 1) RETURN
      !
      ! at least one atom has d functions
      !
      IF (minlmaxab .EQ. 2) THEN
         !
         ! in case both atoms have d functions
         !
         CALL skdd(skab1, mat, iptr(:, :, :, lmaxa, lmaxb))
         CALL sksd(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
         CALL sksd(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
         CALL skpd(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
         CALL skpd(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
      ELSE
         !
         ! One atom has d functions, the other has s or s and p functions
         !
         IF (lmaxa .EQ. 0) THEN
            !
            ! atom b has d, the atom a only s functions
            !
            CALL sksd(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
         ELSE IF (lmaxa .EQ. 1) THEN
            !
            ! atom b has d, the atom a s and p functions
            !
            CALL sksd(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
            CALL skpd(skab2, mat, iptr(:, :, :, lmaxa, lmaxb), .TRUE.)
         ELSE
            !
            ! atom a has d functions
            !
            IF (lmaxb .EQ. 0) THEN
               !
               ! atom a has d, atom b has only s functions
               !
               CALL sksd(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
            ELSE
               !
               ! atom a has d, atom b has s and p functions
               !
               CALL sksd(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
               CALL skpd(skab1, mat, iptr(:, :, :, lmaxa, lmaxb), .FALSE.)
            END IF
         END IF
      END IF
      !
   CONTAINS
      !
      ! The subroutines to turn the matrix elements are taken as internal subroutines
      ! as it is beneficial to inline them.
      !
      ! They are both turning the matrix elements and placing them appropriately
      ! into the matrix block
      !
! **************************************************************************************************
!> \brief   s-s interaction (no rotation necessary)
!> \param skpar ...
!> \param mat ...
!> \version 1.0
! **************************************************************************************************
      SUBROUTINE skss(skpar, mat)
      REAL(dp), INTENT(in)                               :: skpar(:)
      REAL(dp), INTENT(inout)                            :: mat(:, :)

         mat(1, 1) = mat(1, 1) + skpar(1)
         !
      END SUBROUTINE skss

! **************************************************************************************************
!> \brief  s-p interaction (simple rotation)
!> \param skpar ...
!> \param mat ...
!> \param ind ...
!> \param transposed ...
!> \version 1.0
! **************************************************************************************************
      SUBROUTINE sksp(skpar, mat, ind, transposed)
      REAL(dp), INTENT(in)                               :: skpar(:)
      REAL(dp), INTENT(inout)                            :: mat(:, :)
      INTEGER, INTENT(in)                                :: ind(0:, 0:, 0:)
      LOGICAL, INTENT(in)                                :: transposed

      INTEGER                                            :: l
      REAL(dp)                                           :: skp

         skp = skpar(ind(1, 0, 0))
         IF (transposed) THEN
            DO l = 1, 3
               mat(1, l + 1) = mat(1, l + 1) + rr(l)*skp
            END DO
         ELSE
            DO l = 1, 3
               mat(l + 1, 1) = mat(l + 1, 1) - rr(l)*skp
            END DO
         END IF
         !
      END SUBROUTINE sksp

! **************************************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
! **************************************************************************************************
      SUBROUTINE skpp(skpar, mat, ind)
      REAL(dp), INTENT(in)                               :: skpar(:)
      REAL(dp), INTENT(inout)                            :: mat(:, :)
      INTEGER, INTENT(in)                                :: ind(0:, 0:, 0:)

      INTEGER                                            :: ii, ir, is, k, l
      REAL(dp)                                           :: epp(6), matel(6), skppp, skpps

         epp(1:3) = rr2(1:3)
         DO l = 1, 3
            epp(l + 3) = rr(l)*rr(l + 1)
         END DO
         skppp = skpar(ind(1, 1, 1))
         skpps = skpar(ind(1, 1, 0))
         !
         DO l = 1, 3
            matel(l) = epp(l)*skpps + (1._dp - epp(l))*skppp
         END DO
         DO l = 4, 6
            matel(l) = epp(l)*(skpps - skppp)
         END DO
         !
         DO ir = 1, 3
            DO is = 1, ir - 1
               ii = ir - is
               k = 3*ii - (ii*(ii - 1))/2 + is
               mat(is + 1, ir + 1) = mat(is + 1, ir + 1) + matel(k)
               mat(ir + 1, is + 1) = mat(ir + 1, is + 1) + matel(k)
            END DO
            mat(ir + 1, ir + 1) = mat(ir + 1, ir + 1) + matel(ir)
         END DO
      END SUBROUTINE skpp

! **************************************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
!> \param transposed ...
! **************************************************************************************************
      SUBROUTINE sksd(skpar, mat, ind, transposed)
      REAL(dp), INTENT(in)                               :: skpar(:)
      REAL(dp), INTENT(inout)                            :: mat(:, :)
      INTEGER, INTENT(in)                                :: ind(0:, 0:, 0:)
      LOGICAL, INTENT(in)                                :: transposed

      INTEGER                                            :: l
      REAL(dp)                                           :: d4, d5, es(5), r3, sksds

         sksds = skpar(ind(2, 0, 0))
         r3 = SQRT(3._dp)
         d4 = rr2(3) - 0.5_dp*(rr2(1) + rr2(2))
         d5 = rr2(1) - rr2(2)
         !
         DO l = 1, 3
            es(l) = r3*rr(l)*rr(l + 1)
         END DO
         es(4) = 0.5_dp*r3*d5
         es(5) = d4
         !
         IF (transposed) THEN
            DO l = 1, 5
               mat(1, l + 4) = mat(1, l + 4) + es(l)*sksds
            END DO
         ELSE
            DO l = 1, 5
               mat(l + 4, 1) = mat(l + 4, 1) + es(l)*sksds
            END DO
         END IF
      END SUBROUTINE sksd

! **************************************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
!> \param transposed ...
! **************************************************************************************************
      SUBROUTINE skpd(skpar, mat, ind, transposed)
      REAL(dp), INTENT(in)                               :: skpar(:)
      REAL(dp), INTENT(inout)                            :: mat(:, :)
      INTEGER, INTENT(in)                                :: ind(0:, 0:, 0:)
      LOGICAL, INTENT(in)                                :: transposed

      INTEGER                                            :: ir, is, k, l, m
      REAL(dp)                                           :: d3, d4, d5, d6, dm(15), epd(13, 2), r3, &
                                                            sktmp

         r3 = SQRT(3.0_dp)
         d3 = rr2(1) + rr2(2)
         d4 = rr2(3) - 0.5_dp*d3
         d5 = rr2(1) - rr2(2)
         d6 = rr(1)*rr(2)*rr(3)
         DO l = 1, 3
            epd(l, 1) = r3*rr2(l)*rr(l + 1)
            epd(l, 2) = rr(l + 1)*(1.0_dp - 2._dp*rr2(l))
            epd(l + 4, 1) = r3*rr2(l)*rr(l + 2)
            epd(l + 4, 2) = rr(l + 2)*(1.0_dp - 2*rr2(l))
            epd(l + 7, 1) = 0.5_dp*r3*rr(l)*d5
            epd(l + 10, 1) = rr(l)*d4
         END DO
         !
         epd(4, 1) = r3*d6
         epd(4, 2) = -2._dp*d6
         epd(8, 2) = rr(1)*(1.0_dp - d5)
         epd(9, 2) = -rr(2)*(1.0_dp + d5)
         epd(10, 2) = -rr(3)*d5
         epd(11, 2) = -r3*rr(1)*rr2(3)
         epd(12, 2) = -r3*rr(2)*rr2(3)
         epd(13, 2) = r3*rr(3)*d3
         !
         dm(1:15) = 0.0_dp
         !
         DO m = 1, 2
            sktmp = skpar(ind(2, 1, m - 1))
            dm(1) = dm(1) + epd(1, m)*sktmp
            dm(2) = dm(2) + epd(6, m)*sktmp
            dm(3) = dm(3) + epd(4, m)*sktmp
            dm(5) = dm(5) + epd(2, m)*sktmp
            dm(6) = dm(6) + epd(7, m)*sktmp
            dm(7) = dm(7) + epd(5, m)*sktmp
            dm(9) = dm(9) + epd(3, m)*sktmp
            DO l = 8, 13
               dm(l + 2) = dm(l + 2) + epd(l, m)*sktmp
            END DO
         END DO
         !
         dm(4) = dm(3)
         dm(8) = dm(3)
         !
         IF (transposed) THEN
            DO ir = 1, 5
               DO is = 1, 3
                  k = 3*(ir - 1) + is
                  mat(is + 1, ir + 4) = mat(is + 1, ir + 4) + dm(k)
               END DO
            END DO
         ELSE
            DO ir = 1, 5
               DO is = 1, 3
                  k = 3*(ir - 1) + is
                  mat(ir + 4, is + 1) = mat(ir + 4, is + 1) - dm(k)
               END DO
            END DO
         END IF
         !
      END SUBROUTINE skpd

! **************************************************************************************************
!> \brief ...
!> \param skpar ...
!> \param mat ...
!> \param ind ...
! **************************************************************************************************
      SUBROUTINE skdd(skpar, mat, ind)
      REAL(dp), INTENT(in)                               :: skpar(:)
      REAL(dp), INTENT(inout)                            :: mat(:, :)
      INTEGER, INTENT(in)                                :: ind(0:, 0:, 0:)

      INTEGER                                            :: ii, ir, is, k, l, m
      REAL(dp)                                           :: d3, d4, d5, dd(3), dm(15), e(15, 3), r3

         r3 = SQRT(3._dp)
         d3 = rr2(1) + rr2(2)
         d4 = rr2(3) - 0.5_dp*d3
         d5 = rr2(1) - rr2(2)
         DO l = 1, 3
            e(l, 1) = rr2(l)*rr2(l + 1)
            e(l, 2) = rr2(l) + rr2(l + 1) - 4._dp*e(l, 1)
            e(l, 3) = rr2(l + 2) + e(l, 1)
            e(l, 1) = 3._dp*e(l, 1)
         END DO
         e(4, 1) = d5**2
         e(4, 2) = d3 - e(4, 1)
         e(4, 3) = rr2(3) + 0.25_dp*e(4, 1)
         e(4, 1) = 0.75_dp*e(4, 1)
         e(5, 1) = d4**2
         e(5, 2) = 3._dp*rr2(3)*d3
         e(5, 3) = 0.75_dp*d3**2
         dd(1) = rr(1)*rr(3)
         dd(2) = rr(2)*rr(1)
         dd(3) = rr(3)*rr(2)
         DO l = 1, 2
            e(l + 5, 1) = 3._dp*rr2(l + 1)*dd(l)
            e(l + 5, 2) = dd(l)*(1._dp - 4._dp*rr2(l + 1))
            e(l + 5, 3) = dd(l)*(rr2(l + 1) - 1._dp)
         END DO
         e(8, 1) = dd(1)*d5*1.5_dp
         e(8, 2) = dd(1)*(1.0_dp - 2.0_dp*d5)
         e(8, 3) = dd(1)*(0.5_dp*d5 - 1.0_dp)
         e(9, 1) = d5*0.5_dp*d4*r3
         e(9, 2) = -d5*rr2(3)*r3
         e(9, 3) = d5*0.25_dp*(1.0_dp + rr2(3))*r3
         e(10, 1) = rr2(1)*dd(3)*3.0_dp
         e(10, 2) = (0.25_dp - rr2(1))*dd(3)*4.0_dp
         e(10, 3) = dd(3)*(rr2(1) - 1.0_dp)
         e(11, 1) = 1.5_dp*dd(3)*d5
         e(11, 2) = -dd(3)*(1.0_dp + 2.0_dp*d5)
         e(11, 3) = dd(3)*(1.0_dp + 0.5_dp*d5)
         e(13, 3) = 0.5_dp*d5*dd(2)
         e(13, 2) = -2.0_dp*dd(2)*d5
         e(13, 1) = e(13, 3)*3.0_dp
         e(12, 1) = d4*dd(1)*r3
         e(14, 1) = d4*dd(3)*r3
         e(15, 1) = d4*dd(2)*r3
         e(15, 2) = -2.0_dp*r3*dd(2)*rr2(3)
         e(15, 3) = 0.5_dp*r3*(1.0_dp + rr2(3))*dd(2)
         e(14, 2) = r3*dd(3)*(d3 - rr2(3))
         e(14, 3) = -r3*0.5_dp*dd(3)*d3
         e(12, 2) = r3*dd(1)*(d3 - rr2(3))
         e(12, 3) = -r3*0.5_dp*dd(1)*d3
         !
         dm(1:15) = 0._dp
         DO l = 1, 15
            DO m = 1, 3
               dm(l) = dm(l) + e(l, m)*skpar(ind(2, 2, m - 1))
            END DO
         END DO
         !
         DO ir = 1, 5
            DO is = 1, ir - 1
               ii = ir - is
               k = 5*ii - (ii*(ii - 1))/2 + is
               mat(ir + 4, is + 4) = mat(ir + 4, is + 4) + dm(k)
               mat(is + 4, ir + 4) = mat(is + 4, ir + 4) + dm(k)
            END DO
            mat(ir + 4, ir + 4) = mat(ir + 4, ir + 4) + dm(ir)
         END DO
      END SUBROUTINE skdd
      !
   END SUBROUTINE turnsk

! **************************************************************************************************
!> \brief ...
!> \param xa ...
!> \param ya ...
!> \param n ...
!> \param x ...
!> \param y ...
!> \param dy ...
! **************************************************************************************************
   SUBROUTINE polint(xa, ya, n, x, y, dy)
      INTEGER, INTENT(in)                                :: n
      REAL(dp), INTENT(in)                               :: ya(n), xa(n), x
      REAL(dp), INTENT(out)                              :: y, dy

      CHARACTER(len=*), PARAMETER :: routineN = 'polint', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, m, ns
      REAL(dp)                                           :: c(n), d(n), den, dif, dift, ho, hp, w

!
!

      ns = 1

      dif = ABS(x - xa(1))
      DO i = 1, n
         dift = ABS(x - xa(i))
         IF (dift .LT. dif) THEN
            ns = i
            dif = dift
         ENDIF
         c(i) = ya(i)
         d(i) = ya(i)
      END DO
      !
      y = ya(ns)
      ns = ns - 1
      DO m = 1, n - 1
         DO i = 1, n - m
            ho = xa(i) - x
            hp = xa(i + m) - x
            w = c(i + 1) - d(i)
            den = ho - hp
            CPASSERT(den /= 0.0_dp)
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         END DO
         IF (2*ns .LT. n - m) THEN
            dy = c(ns + 1)
         ELSE
            dy = d(ns)
            ns = ns - 1
         ENDIF
         y = y + dy
      END DO
!
      RETURN
   END SUBROUTINE polint

! **************************************************************************************************
!> \brief ...
!> \param rv ...
!> \param r ...
!> \param erep ...
!> \param derep ...
!> \param n_urpoly ...
!> \param urep ...
!> \param spdim ...
!> \param s_cut ...
!> \param srep ...
!> \param spxr ...
!> \param scoeff ...
!> \param surr ...
!> \param dograd ...
! **************************************************************************************************
   SUBROUTINE urep_egr(rv, r, erep, derep, &
                       n_urpoly, urep, spdim, s_cut, srep, spxr, scoeff, surr, dograd)

      REAL(dp), INTENT(in)                               :: rv(3), r
      REAL(dp), INTENT(inout)                            :: erep, derep(3)
      INTEGER, INTENT(in)                                :: n_urpoly
      REAL(dp), INTENT(in)                               :: urep(:)
      INTEGER, INTENT(in)                                :: spdim
      REAL(dp), INTENT(in)                               :: s_cut, srep(3)
      REAL(dp), POINTER                                  :: spxr(:, :), scoeff(:, :)
      REAL(dp), INTENT(in)                               :: surr(2)
      LOGICAL, INTENT(in)                                :: dograd

      INTEGER                                            :: ic, isp, jsp, nsp
      REAL(dp)                                           :: de_z, rz

      derep = 0._dp
      de_z = 0._dp
      IF (n_urpoly > 0) THEN
         !
         ! polynomial part
         !
         rz = urep(1) - r
         IF (rz <= rtiny) RETURN
         DO ic = 2, n_urpoly
            erep = erep + urep(ic)*rz**(ic)
         END DO
         IF (dograd) THEN
            DO ic = 2, n_urpoly
               de_z = de_z - ic*urep(ic)*rz**(ic - 1)
            END DO
         END IF
      ELSE IF (spdim > 0) THEN
         !
         ! spline part
         !
         ! This part is kind of proprietary Paderborn code and I won't reverse-engineer
         ! everything in detail. What is obvious is documented.
         !
         ! This part has 4 regions:
         ! a) very long range is screened
         ! b) short-range is extrapolated with e-functions
         ! ca) normal range is approximated with a spline
         ! cb) longer range is extrapolated with an higher degree spline
         !
         IF (r > s_cut) RETURN ! screening (condition a)
         !
         IF (r < spxr(1, 1)) THEN
            ! a) short range
            erep = erep + EXP(-srep(1)*r + srep(2)) + srep(3)
            IF (dograd) de_z = de_z - srep(1)*EXP(-srep(1)*r + srep(2))
         ELSE
            !
            ! condition c). First determine between which places the spline is located:
            !
            ispg: DO isp = 1, spdim ! condition ca)
               IF (r < spxr(isp, 1)) CYCLE ispg ! distance is smaller than this spline range
               IF (r >= spxr(isp, 2)) CYCLE ispg ! distance is larger than this spline range
               ! at this point we have found the correct spline interval
               rz = r - spxr(isp, 1)
               IF (isp /= spdim) THEN
                  nsp = 3 ! condition ca
                  DO jsp = 0, nsp
                     erep = erep + scoeff(isp, jsp + 1)*rz**(jsp)
                  END DO
                  IF (dograd) THEN
                     DO jsp = 1, nsp
                        de_z = de_z + jsp*scoeff(isp, jsp + 1)*rz**(jsp - 1)
                     END DO
                  END IF
               ELSE
                  nsp = 5 ! condition cb
                  DO jsp = 0, nsp
                     IF (jsp <= 3) THEN
                        erep = erep + scoeff(isp, jsp + 1)*rz**(jsp)
                     ELSE
                        erep = erep + surr(jsp - 3)*rz**(jsp)
                     ENDIF
                  END DO
                  IF (dograd) THEN
                     DO jsp = 1, nsp
                        IF (jsp <= 3) THEN
                           de_z = de_z + jsp*scoeff(isp, jsp + 1)*rz**(jsp - 1)
                        ELSE
                           de_z = de_z + jsp*surr(jsp - 3)*rz**(jsp - 1)
                        ENDIF
                     END DO
                  END IF
               END IF
               EXIT ispg
            END DO ispg
         END IF
      END IF
      !
      IF (dograd) THEN
         IF (r > 1.e-12_dp) derep(1:3) = (de_z/r)*rv(1:3)
      END IF

   END SUBROUTINE urep_egr

END MODULE qs_dftb_utils

