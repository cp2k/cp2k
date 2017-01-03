!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Some basic routines for atomic calculations
!> \author  jgh
!> \date    01.04.2008
!> \version 1.0
!>
! **************************************************************************************************
MODULE atom_utils
   USE ai_onecenter,                    ONLY: sg_overlap,&
                                              sto_overlap
   USE ai_overlap,                      ONLY: overlap_ab_s,&
                                              overlap_ab_sp
   USE ao_util,                         ONLY: exp_radius
   USE atom_types,                      ONLY: &
        CGTO_BASIS, GTO_BASIS, NUM_BASIS, STO_BASIS, atom_basis_type, atom_gthpot_type, &
        atom_potential_type, atom_state, atom_type, ecp_pseudo, eri, gth_pseudo, no_pseudo, &
        upf_pseudo
   USE basis_set_types,                 ONLY: srules
   USE cp_files,                        ONLY: close_file,&
                                              get_unit_number,&
                                              open_file
   USE erf_fn,                          ONLY: erf
   USE input_constants,                 ONLY: do_rhf_atom,&
                                              do_rks_atom,&
                                              do_rohf_atom,&
                                              do_uhf_atom,&
                                              do_uks_atom
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lapack,                          ONLY: lapack_ssyev
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              fourpi,&
                                              maxfac,&
                                              rootpi
   USE mathlib,                         ONLY: invmat_symm
   USE orbital_pointers,                ONLY: deallocate_orbital_pointers,&
                                              init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: deallocate_spherical_harmonics,&
                                              init_spherical_harmonics
   USE periodic_table,                  ONLY: nelem,&
                                              ptable
   USE physcon,                         ONLY: bohr
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE splines,                         ONLY: spline3ders
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_utils'

   PUBLIC :: atom_condnumber, atom_completeness, atom_basis_condnum
   PUBLIC :: atom_set_occupation, get_maxl_occ, get_maxn_occ
   PUBLIC :: atom_denmat, atom_density, atom_core_density
   PUBLIC :: integrate_grid, atom_trace, atom_solve
   PUBLIC :: coulomb_potential_numeric, coulomb_potential_analytic
   PUBLIC :: exchange_numeric, exchange_semi_analytic
   PUBLIC :: numpot_matrix, ceri_contract, eeri_contract, err_matrix
   PUBLIC :: slater_density, wigner_slater_functional, atom_local_potential
   PUBLIC :: atom_orbital_charge, atom_orbital_nodes, atom_consistent_method
   PUBLIC :: atom_orbital_max, atom_wfnr0, get_rho0
   PUBLIC :: contract2, contract4, contract2add
! ZMP added public subroutines
   PUBLIC :: atom_read_external_density
   PUBLIC :: atom_read_external_vxc
   PUBLIC :: atom_read_zmp_restart
   PUBLIC :: atom_write_zmp_restart

!-----------------------------------------------------------------------------!

   INTERFACE integrate_grid
      MODULE PROCEDURE integrate_grid_function1, &
         integrate_grid_function2, &
         integrate_grid_function3
   END INTERFACE

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ostring ...
!> \param occupation ...
!> \param wfnocc ...
!> \param multiplicity ...
! **************************************************************************************************
   SUBROUTINE atom_set_occupation(ostring, occupation, wfnocc, multiplicity)
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: ostring
      REAL(Kind=dp), DIMENSION(0:3, 10)                  :: occupation, wfnocc
      INTEGER, INTENT(OUT), OPTIONAL                     :: multiplicity

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_set_occupation', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: elem
      CHARACTER(LEN=default_string_length)               :: pstring
      INTEGER                                            :: i, i1, i2, ielem, is, jd, jf, jp, js, k, &
                                                            l, mult, n, no
      REAL(Kind=dp)                                      :: e0, el, oo

      occupation = 0._dp

      CPASSERT(ASSOCIATED(ostring))
      CPASSERT(SIZE(ostring) > 0)

      no = SIZE(ostring)

      is = 1
      ! look for multiplicity
      mult = -1 !not specified
      IF (is <= no) THEN
         IF (INDEX(ostring(is), "(") /= 0) THEN
            i1 = INDEX(ostring(is), "(")
            i2 = INDEX(ostring(is), ")")
            CPASSERT((i2-i1-1 > 0) .AND. (i2-i1-1 < 3))
            elem = ostring(is) (i1+1:i2-1)
            IF (INDEX(elem, "HS") /= 0) THEN
               mult = -2 !High spin
            ELSE IF (INDEX(elem, "LS") /= 0) THEN
               mult = -3 !Low spin
            ELSE
               READ (elem, *) mult
            END IF
            is = is+1
         END IF
      END IF

      IF (is <= no) THEN
         IF (INDEX(ostring(is), "CORE") /= 0) is = is+1 !Pseudopotential detected
      END IF
      IF (is <= no) THEN
         IF (INDEX(ostring(is), "none") /= 0) is = is+1 !no electrons, used with CORE
      END IF
      IF (is <= no) THEN
         IF (INDEX(ostring(is), "[") /= 0) THEN
            ! core occupation from element [XX]
            i1 = INDEX(ostring(is), "[")
            i2 = INDEX(ostring(is), "]")
            CPASSERT((i2-i1-1 > 0) .AND. (i2-i1-1 < 3))
            elem = ostring(is) (i1+1:i2-1)
            ielem = 0
            DO k = 1, nelem
               IF (elem == ptable(k)%symbol) THEN
                  ielem = k
                  EXIT
               END IF
            END DO
            CPASSERT(ielem /= 0)
            DO l = 0, 3
               el = 2._dp*(2._dp*REAL(l, dp)+1._dp)
               e0 = ptable(ielem)%e_conv(l)
               DO k = 1, 10
                  occupation(l, k) = MIN(el, e0)
                  e0 = e0-el
                  IF (e0 <= 0._dp) EXIT
               END DO
            END DO
            is = is+1
         END IF

      END IF

      DO i = is, no
         pstring = ostring(i)
         CALL uppercase(pstring)
         js = INDEX(pstring, "S")
         jp = INDEX(pstring, "P")
         jd = INDEX(pstring, "D")
         jf = INDEX(pstring, "F")
         CPASSERT(js+jp+jd+jf > 0)
         IF (js > 0) THEN
            CPASSERT(jp+jd+jf == 0)
            READ (pstring(1:js-1), *) n
            READ (pstring(js+1:), *) oo
            CPASSERT(n > 0)
            CPASSERT(oo >= 0._dp)
            CPASSERT(occupation(0, n) == 0)
            occupation(0, n) = oo
         END IF
         IF (jp > 0) THEN
            CPASSERT(js+jd+jf == 0)
            READ (pstring(1:jp-1), *) n
            READ (pstring(jp+1:), *) oo
            CPASSERT(n > 1)
            CPASSERT(oo >= 0._dp)
            CPASSERT(occupation(1, n-1) == 0)
            occupation(1, n-1) = oo
         END IF
         IF (jd > 0) THEN
            CPASSERT(js+jp+jf == 0)
            READ (pstring(1:jd-1), *) n
            READ (pstring(jd+1:), *) oo
            CPASSERT(n > 2)
            CPASSERT(oo >= 0._dp)
            CPASSERT(occupation(2, n-2) == 0)
            occupation(2, n-2) = oo
         END IF
         IF (jf > 0) THEN
            CPASSERT(js+jp+jd == 0)
            READ (pstring(1:jf-1), *) n
            READ (pstring(jf+1:), *) oo
            CPASSERT(n > 3)
            CPASSERT(oo >= 0._dp)
            CPASSERT(occupation(3, n-3) == 0)
            occupation(3, n-3) = oo
         END IF

      END DO

      wfnocc = 0._dp
      DO l = 0, 3
         k = 0
         DO i = 1, 10
            IF (occupation(l, i) /= 0._dp) THEN
               k = k+1
               wfnocc(l, k) = occupation(l, i)
            END IF
         END DO
      END DO

      !Check for consistency with multiplicity
      IF (mult /= -1) THEN
         ! count open shells
         js = 0
         DO l = 0, 3
            k = 2*(2*l+1)
            DO i = 1, 10
               IF (wfnocc(l, i) /= 0._dp .AND. wfnocc(l, i) /= REAL(k, dp)) THEN
                  js = js+1
                  i1 = l
                  i2 = i
               END IF
            END DO
         END DO

         IF (js == 0 .AND. mult == -2) mult = 1
         IF (js == 0 .AND. mult == -3) mult = 1
         IF (js == 0) THEN
            CPASSERT(mult == 1)
         END IF
         IF (js == 1) THEN
            l = i1
            i = i2
            k = NINT(wfnocc(l, i))
            IF (k > (2*l+1)) k = 2*(2*l+1)-k
            IF (mult == -2) mult = k+1
            IF (mult == -3) mult = MOD(k, 2)+1
            CPASSERT(MOD(k+1-mult, 2) == 0)
         END IF
         IF (js > 1 .AND. mult /= -2) THEN
            CPASSERT(mult == -2)
         END IF

      END IF

      IF (PRESENT(multiplicity)) multiplicity = mult

   END SUBROUTINE atom_set_occupation

! **************************************************************************************************
!> \brief ...
!> \param occupation ...
!> \retval maxl ...
! **************************************************************************************************
   FUNCTION get_maxl_occ(occupation) RESULT(maxl)
      REAL(Kind=dp), DIMENSION(0:3, 10)                  :: occupation
      INTEGER                                            :: maxl

      INTEGER                                            :: l

      maxl = 0
      DO l = 0, 3
         IF (SUM(occupation(l, :)) /= 0._dp) maxl = l
      END DO

   END FUNCTION get_maxl_occ

! **************************************************************************************************
!> \brief ...
!> \param occupation ...
!> \retval maxn ...
! **************************************************************************************************
   FUNCTION get_maxn_occ(occupation) RESULT(maxn)
      REAL(Kind=dp), DIMENSION(0:3, 10)                  :: occupation
      INTEGER, DIMENSION(0:3)                            :: maxn

      INTEGER                                            :: k, l

      maxn = 0
      DO l = 0, 3
         DO k = 1, 10
            IF (occupation(l, k) /= 0._dp) maxn(l) = maxn(l)+1
         END DO
      END DO

   END FUNCTION get_maxn_occ

! **************************************************************************************************
!> \brief ...
!> \param pmat ...
!> \param wfn ...
!> \param nbas ...
!> \param occ ...
!> \param maxl ...
!> \param maxn ...
! **************************************************************************************************
   SUBROUTINE atom_denmat(pmat, wfn, nbas, occ, maxl, maxn)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: pmat
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: wfn
      INTEGER, DIMENSION(0:3), INTENT(IN)                :: nbas
      REAL(KIND=dp), DIMENSION(0:, :), INTENT(IN)        :: occ
      INTEGER, INTENT(IN)                                :: maxl
      INTEGER, DIMENSION(0:3), INTENT(IN)                :: maxn

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_denmat', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, l, n

      pmat = 0._dp
      n = SIZE(wfn, 2)
      DO l = 0, maxl
         DO i = 1, MIN(n, maxn(l))
            DO k = 1, nbas(l)
               DO j = 1, nbas(l)
                  pmat(j, k, l) = pmat(j, k, l)+occ(l, i)*wfn(j, i, l)*wfn(k, i, l)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE atom_denmat
! **************************************************************************************************
!> \brief ...
!> \param density ...
!> \param pmat ...
!> \param basis ...
!> \param maxl ...
!> \param typ ...
!> \param rr ...
! **************************************************************************************************
   SUBROUTINE atom_density(density, pmat, basis, maxl, typ, rr)
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: density
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: pmat
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      INTEGER, INTENT(IN)                                :: maxl
      CHARACTER(LEN=*), OPTIONAL                         :: typ
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: rr

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_density', routineP = moduleN//':'//routineN

      CHARACTER(LEN=3)                                   :: my_typ
      INTEGER                                            :: i, j, l, n
      REAL(KIND=dp)                                      :: ff

      my_typ = "RHO"
      IF (PRESENT(typ)) my_typ = typ(1:3)
      IF (my_typ == "KIN") THEN
         CPASSERT(PRESENT(rr))
      END IF

      density = 0._dp
      DO l = 0, maxl
         n = basis%nbas(l)
         DO i = 1, n
            DO j = i, n
               ff = pmat(i, j, l)
               IF (i /= j) ff = 2._dp*pmat(i, j, l)
               IF (my_typ == "RHO") THEN
                  density(:) = density(:)+ff*basis%bf(:, i, l)*basis%bf(:, j, l)
               ELSE IF (my_typ == "DER") THEN
                  density(:) = density(:)+ff*basis%dbf(:, i, l)*basis%bf(:, j, l) &
                               +ff*basis%bf(:, i, l)*basis%dbf(:, j, l)
               ELSE IF (my_typ == "KIN") THEN
                  density(:) = density(:)+0.5_dp*ff*( &
                               basis%dbf(:, i, l)*basis%dbf(:, j, l)+ &
                               REAL(l*(l+1), dp)*basis%bf(:, i, l)*basis%bf(:, j, l)/rr(:))
               ELSE IF (my_typ == "LAP") THEN
                  density(:) = density(:)+ff*basis%ddbf(:, i, l)*basis%bf(:, j, l) &
                               +ff*basis%bf(:, i, l)*basis%ddbf(:, j, l) &
                               +2._dp*ff*basis%dbf(:, i, l)*basis%bf(:, j, l)/rr(:) &
                               +2._dp*ff*basis%bf(:, i, l)*basis%dbf(:, j, l)/rr(:)
               ELSE
                  CPABORT("")
               END IF
            END DO
         END DO
      END DO
      ! this factor from the product of two spherical harmonics
      density = density/fourpi

   END SUBROUTINE atom_density
! **************************************************************************************************
!> \brief   ZMP subroutine to write external restart file
!>
!> \param atom ...
!> \date    07.10.2013
!> \author  D. Varsano [daniele.varsano@nano.cnr.it]
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE atom_write_zmp_restart(atom)
      TYPE(atom_type), INTENT(IN)                        :: atom

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_write_zmp_restart', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: extunit, i, k, l, n

      extunit = get_unit_number()
      CALL open_file(file_name=atom%zmp_restart_file, file_status="UNKNOWN", &
                     file_form="FORMATTED", file_action="WRITE", &
                     unit_number=extunit)

      n = SIZE(atom%orbitals%wfn, 2)
      WRITE (extunit, *) atom%basis%nbas
      DO l = 0, atom%state%maxl_occ
         DO i = 1, MIN(n, atom%state%maxn_occ(l))
            DO k = 1, atom%basis%nbas(l)
               WRITE (extunit, *) atom%orbitals%wfn(k, i, l)
            END DO
         END DO
      END DO

      CALL close_file(unit_number=extunit)

   END SUBROUTINE atom_write_zmp_restart

! **************************************************************************************************
!> \brief   ZMP subroutine to read external restart file
!>
!> \param atom ...
!> \param doguess ...
!> \param iw ...
!> \date    07.10.2013
!> \author  D. Varsano [daniele.varsano@nano.cnr.it]
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE atom_read_zmp_restart(atom, doguess, iw)
      TYPE(atom_type), INTENT(INOUT)                     :: atom
      LOGICAL, INTENT(INOUT)                             :: doguess
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_zmp_restart', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: er, extunit, i, k, l, maxl, n
      INTEGER, DIMENSION(0:3)                            :: maxn, nbas

      INQUIRE (file=atom%zmp_restart_file, exist=atom%doread)

      IF (atom%doread) THEN
         WRITE (iw, fmt="(' ZMP       | Restart file found')")
         extunit = get_unit_number()
         CALL open_file(file_name=atom%zmp_restart_file, file_status="OLD", &
                        file_form="FORMATTED", file_action="READ", &
                        unit_number=extunit)

         READ (extunit, *, IOSTAT=er) nbas

         IF (er .NE. 0) THEN
            WRITE (iw, fmt="(' ZMP       | ERROR! Restart file unreadable')")
            WRITE (iw, fmt="(' ZMP       | ERROR! Starting ZMP calculation form initial atomic guess')")
            doguess = .TRUE.
            atom%doread = .FALSE.
         ELSE IF (nbas(1) .NE. atom%basis%nbas(1)) THEN
            WRITE (iw, fmt="(' ZMP       | ERROR! Restart file contains a different basis set')")
            WRITE (iw, fmt="(' ZMP       | ERROR! Starting ZMP calculation form initial atomic guess')")
            doguess = .TRUE.
            atom%doread = .FALSE.
         ELSE
            nbas = atom%basis%nbas
            maxl = atom%state%maxl_occ
            maxn = atom%state%maxn_occ
            n = SIZE(atom%orbitals%wfn, 2)
            DO l = 0, atom%state%maxl_occ
               DO i = 1, MIN(n, atom%state%maxn_occ(l))
                  DO k = 1, atom%basis%nbas(l)
                     READ (extunit, *) atom%orbitals%wfn(k, i, l)
                  END DO
               END DO
            END DO
            doguess = .FALSE.
         END IF
         CALL close_file(unit_number=extunit)
      ELSE
         WRITE (iw, fmt="(' ZMP       | WARNING! Restart file not found')")
         WRITE (iw, fmt="(' ZMP       | WARNING! Starting ZMP calculation form initial atomic guess')")
      ENDIF
   END SUBROUTINE atom_read_zmp_restart
! **************************************************************************************************
!> \brief   ZMP subroutine to read external density from linear grid of
!>          density matrix.
!> \param density ...
!> \param atom ...
!> \param iw ...
!> \date    07.10.2013
!> \author  D. Varsano [daniele.varsano@nano.cnr.it]
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE atom_read_external_density(density, atom, iw)
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: density
      TYPE(atom_type), INTENT(INOUT)                     :: atom
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_external_density', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: filename
      INTEGER                                            :: extunit, ir, j, k, l, maxl_occ, maxnbas, &
                                                            nbas, nr
      LOGICAL                                            :: ldm
      REAL(KIND=dp)                                      :: rr
      REAL(KIND=dp), ALLOCATABLE                         :: pmatread(:, :, :)

      filename = atom%ext_file
      ldm = atom%dm
      extunit = get_unit_number()

      CALL open_file(file_name=filename, file_status="OLD", &
                     file_form="FORMATTED", file_action="READ", &
                     unit_number=extunit)

      IF (.NOT. ldm) THEN
         READ (extunit, *) nr

         IF (nr .NE. atom%basis%grid%nr) THEN
            IF (iw > 0) WRITE (iw, fmt="(' ZMP       | ERROR! External grid dimension ',I4,' differs from internal grid ',I4)") &
               nr, atom%basis%grid%nr
            IF (iw > 0) WRITE (iw, fmt="(' ZMP       | ERROR! Stopping ZMP calculation')")
            CPABORT("")
         ENDIF

         DO ir = 1, nr
            READ (extunit, *) rr, density(ir)
            IF (ABS(rr-atom%basis%grid%rad(ir)) .GT. atom%zmpgrid_tol) THEN
               IF (iw > 0) WRITE (iw, fmt="(' ZMP       | ERROR! Grid points do not coincide: ')")
               IF (iw > 0) WRITE (iw, fmt='(" ZMP       |",T20,"R_out[bohr]",T36,"R_in[bohr]",T61,"R_diff[bohr]")')
               IF (iw > 0) WRITE (iw, fmt='(" ZMP       |",T14,E24.15,T39,E24.15,T64,E24.15)') &
                  rr, atom%basis%grid%rad(ir), ABS(rr-atom%basis%grid%rad(ir))
               CPABORT("")
            ENDIF
         ENDDO
         CALL close_file(unit_number=extunit)
      ELSE
         READ (extunit, *) maxl_occ
         maxnbas = MAXVAL(atom%basis%nbas)
         ALLOCATE (pmatread(maxnbas, maxnbas, 0:maxl_occ))
         pmatread = 0.0
         DO l = 0, maxl_occ
            nbas = atom%basis%nbas(l)
            READ (extunit, *) ! Read empty line
            DO k = 1, nbas
               READ (extunit, *) (pmatread(j, k, l), j=1, k)
               DO j = 1, k
                  pmatread(k, j, l) = pmatread(j, k, l)
               ENDDO
            ENDDO
         ENDDO

         CALL close_file(unit_number=extunit)

         CALL atom_density(density, pmatread, atom%basis, maxl_occ, typ="RHO")

         extunit = get_unit_number()
         CALL open_file(file_name="rho_target.dat", file_status="UNKNOWN", &
                        file_form="FORMATTED", file_action="WRITE", unit_number=extunit)

         IF (iw > 0) WRITE (iw, fmt="(' ZMP       | Writing target density from density matrix')")

         WRITE (extunit, fmt='("# Target density built from density matrix : ",A20)') filename
         WRITE (extunit, fmt='("#",T10,"R[bohr]",T36,"Rho[au]")')

         nr = atom%basis%grid%nr

         DO ir = 1, nr
            WRITE (extunit, fmt='(T1,E24.15,T26,E24.15)') &
               atom%basis%grid%rad(ir), density(ir)
         ENDDO
         DEALLOCATE (pmatread)

         CALL close_file(unit_number=extunit)

      ENDIF

   END SUBROUTINE atom_read_external_density

! **************************************************************************************************
!> \brief ZMP subroutine to read external v_xc in the atomic code
!>
!> \param vxc ...
!> \param atom ...
!> \param iw ...
!> \author D. Varsano [daniele.varsano@nano.cnr.it]
! **************************************************************************************************
   SUBROUTINE atom_read_external_vxc(vxc, atom, iw)
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: vxc
      TYPE(atom_type), INTENT(INOUT)                     :: atom
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_external_vxc', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: adum, filename
      INTEGER                                            :: extunit, ir, nr
      REAL(KIND=dp)                                      :: rr

      filename = atom%ext_vxc_file
      extunit = get_unit_number()

      CALL open_file(file_name=filename, file_status="OLD", &
                     file_form="FORMATTED", file_action="READ", &
                     unit_number=extunit)

      READ (extunit, *)
      READ (extunit, *) adum, nr

      IF (nr .NE. atom%basis%grid%nr) THEN
         IF (iw > 0) WRITE (iw, fmt="(' ZMP       | ERROR! External grid dimension ',I4,' differs from internal grid ',I4)") &
            nr, atom%basis%grid%nr
         IF (iw > 0) WRITE (iw, fmt="(' ZMP       | ERROR! Stopping ZMP calculation')")
         CPABORT("")
      ENDIF
      DO ir = 1, nr
         READ (extunit, *) rr, vxc(ir)
         IF (ABS(rr-atom%basis%grid%rad(ir)) .GT. atom%zmpvxcgrid_tol) THEN
            IF (iw > 0) WRITE (iw, fmt="(' ZMP       | ERROR! Grid points do not coincide: ')")
            IF (iw > 0) WRITE (iw, fmt='(" ZMP       |",T20,"R_out[bohr]",T36,"R_in[bohr]",T61,"R_diff[bohr]")')
            IF (iw > 0) WRITE (iw, fmt='(" ZMP       |",T14,E24.15,T39,E24.15,T64,E24.15)') &
               rr, atom%basis%grid%rad(ir), ABS(rr-atom%basis%grid%rad(ir))
            CPABORT("")
         ENDIF
      ENDDO

   END SUBROUTINE atom_read_external_vxc

! **************************************************************************************************
!> \brief ...
!> \param charge ...
!> \param wfn ...
!> \param rcov ...
!> \param l ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE atom_orbital_charge(charge, wfn, rcov, l, basis)
      REAL(KIND=dp), INTENT(OUT)                         :: charge
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wfn
      REAL(KIND=dp), INTENT(IN)                          :: rcov
      INTEGER, INTENT(IN)                                :: l
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_orbital_charge', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, m, n
      REAL(KIND=dp)                                      :: ff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: den

      charge = 0._dp
      m = SIZE(basis%bf, 1)
      ALLOCATE (den(m))
      n = basis%nbas(l)
      den = 0._dp
      DO i = 1, n
         DO j = 1, n
            ff = wfn(i)*wfn(j)
            den(1:m) = den(1:m)+ff*basis%bf(1:m, i, l)*basis%bf(1:m, j, l)
         END DO
      END DO
      DO i = 1, m
         IF (basis%grid%rad(i) > rcov) den(i) = 0._dp
      END DO
      charge = SUM(den(1:m)*basis%grid%wr(1:m))
      DEALLOCATE (den)

   END SUBROUTINE atom_orbital_charge
! **************************************************************************************************
!> \brief ...
!> \param corden ...
!> \param potential ...
!> \param typ ...
!> \param rr ...
! **************************************************************************************************
   SUBROUTINE atom_core_density(corden, potential, typ, rr)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: corden
      TYPE(atom_potential_type), INTENT(IN)              :: potential
      CHARACTER(LEN=*), OPTIONAL                         :: typ
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rr

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_core_density', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=3)                                   :: my_typ
      INTEGER                                            :: i, j, m, n
      LOGICAL                                            :: reverse
      REAL(KIND=dp)                                      :: a, a2, cval, fb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fe, rc, rhoc, rval

      my_typ = "RHO"
      IF (PRESENT(typ)) my_typ = typ(1:3)

      SELECT CASE (potential%ppot_type)
      CASE (no_pseudo, ecp_pseudo)
         ! we do nothing
      CASE (gth_pseudo)
         IF (potential%gth_pot%nlcc) THEN
            m = SIZE(corden)
            ALLOCATE (fe(m), rc(m))
            n = potential%gth_pot%nexp_nlcc
            DO i = 1, n
               a = potential%gth_pot%alpha_nlcc(i)
               a2 = a*a
               ! note that all terms are computed with rc, not rr
               rc(:) = rr(:)/a
               fe(:) = EXP(-0.5_dp*rc(:)*rc(:))
               DO j = 1, potential%gth_pot%nct_nlcc(i)
                  cval = potential%gth_pot%cval_nlcc(j, i)
                  IF (my_typ == "RHO") THEN
                     corden(:) = corden(:)+fe(:)*rc**(2*j-2)*cval
                  ELSE IF (my_typ == "DER") THEN
                     corden(:) = corden(:)-fe(:)*rc**(2*j-1)*cval/a
                     IF (j > 1) THEN
                        corden(:) = corden(:)+REAL(2*j-2, dp)*fe(:)*rc**(2*j-3)*cval/a
                     END IF
                  ELSE IF (my_typ == "LAP") THEN
                     fb = 2._dp*cval/a
                     corden(:) = corden(:)-fb*fe(:)/rr(:)*rc**(2*j-1)
                     corden(:) = corden(:)+fe(:)*rc**(2*j)*cval/a2
                     IF (j > 1) THEN
                        corden(:) = corden(:)+fb*REAL(2*j-2, dp)*fe(:)/rr(:)*rc**(2*j-3)
                        corden(:) = corden(:)+REAL((2*j-2)*(2*j-3), dp)*fe(:)*rc**(2*j-4)*cval/a2
                        corden(:) = corden(:)-REAL(2*j-2, dp)*fe(:)*rc**(2*j-2)*cval/a2
                     END IF
                  ELSE
                     CPABORT("")
                  END IF
               END DO
            END DO
            DEALLOCATE (fe, rc)
         END IF
      CASE (upf_pseudo)
         IF (potential%upf_pot%core_correction) THEN
            m = SIZE(corden)
            n = potential%upf_pot%mesh_size
            reverse = .FALSE.
            IF (rr(1) > rr(m)) reverse = .TRUE.
            ALLOCATE (rhoc(m), rval(m))
            IF (reverse) THEN
               DO i = 1, m
                  rval(i) = rr(m-i+1)
               END DO
            ELSE
               rval(1:m) = rr(1:m)
            END IF
            IF (my_typ == "RHO") THEN
               CALL spline3ders(potential%upf_pot%r(1:n), potential%upf_pot%rho_nlcc(1:n), rval(1:m), &
                                ynew=rhoc(1:m))
            ELSE IF (my_typ == "DER") THEN
               CALL spline3ders(potential%upf_pot%r(1:n), potential%upf_pot%rho_nlcc(1:n), rval(1:m), &
                                dynew=rhoc(1:m))
            ELSE IF (my_typ == "LAP") THEN
               CALL spline3ders(potential%upf_pot%r(1:n), potential%upf_pot%rho_nlcc(1:n), rval(1:m), &
                                d2ynew=rhoc(1:m))
            ELSE
               CPABORT("")
            END IF
            IF (reverse) THEN
               DO i = 1, m
                  rval(i) = rr(m-i+1)
                  corden(i) = corden(i)+rhoc(m-i+1)
               END DO
            ELSE
               corden(1:m) = corden(1:m)+rhoc(1:m)
            END IF
            DEALLOCATE (rhoc, rval)
         END IF
      CASE DEFAULT
         CPABORT("Unknown PP type")
      END SELECT

   END SUBROUTINE atom_core_density
! **************************************************************************************************
!> \brief ...
!> \param locpot ...
!> \param gthpot ...
!> \param rr ...
! **************************************************************************************************
   SUBROUTINE atom_local_potential(locpot, gthpot, rr)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: locpot
      TYPE(atom_gthpot_type), INTENT(IN)                 :: gthpot
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rr

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_local_potential', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, m, n
      REAL(KIND=dp)                                      :: a
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fe, rc

      m = SIZE(locpot)
      ALLOCATE (fe(m), rc(m))
      rc(:) = rr(:)/gthpot%rc
      DO i = 1, m
         locpot(i) = -gthpot%zion*erf(rc(i)/SQRT(2._dp))/rr(i)
      END DO
      n = gthpot%ncl
      fe(:) = EXP(-0.5_dp*rc(:)*rc(:))
      DO i = 1, n
         locpot(:) = locpot(:)+fe(:)*rc**(2*i-2)*gthpot%cl(i)
      END DO
      IF (gthpot%lpotextended) THEN
         DO j = 1, gthpot%nexp_lpot
            a = gthpot%alpha_lpot(j)
            rc(:) = rr(:)/a
            fe(:) = EXP(-0.5_dp*rc(:)*rc(:))
            n = gthpot%nct_lpot(j)
            DO i = 1, n
               locpot(:) = locpot(:)+fe(:)*rc**(2*i-2)*gthpot%cval_lpot(i, j)
            END DO
         END DO
      END IF
      DEALLOCATE (fe, rc)

   END SUBROUTINE atom_local_potential
! **************************************************************************************************
!> \brief ...
!> \param rmax ...
!> \param wfn ...
!> \param rcov ...
!> \param l ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE atom_orbital_max(rmax, wfn, rcov, l, basis)
      REAL(KIND=dp), INTENT(OUT)                         :: rmax
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wfn
      REAL(KIND=dp), INTENT(IN)                          :: rcov
      INTEGER, INTENT(IN)                                :: l
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_orbital_max', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, m, n
      REAL(KIND=dp)                                      :: ff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dorb

      m = SIZE(basis%bf, 1)
      ALLOCATE (dorb(m))
      n = basis%nbas(l)
      dorb = 0._dp
      DO i = 1, n
         ff = wfn(i)
         dorb(1:m) = dorb(1:m)+ff*basis%dbf(1:m, i, l)
      END DO
      rmax = -1._dp
      DO i = 1, m-1
         IF (basis%grid%rad(i) < 2*rcov) THEN
            IF (dorb(i)*dorb(i+1) < 0._dp) THEN
               rmax = MAX(rmax, basis%grid%rad(i))
            END IF
         END IF
      END DO
      DEALLOCATE (dorb)

   END SUBROUTINE atom_orbital_max
! **************************************************************************************************
!> \brief ...
!> \param node ...
!> \param wfn ...
!> \param rcov ...
!> \param l ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE atom_orbital_nodes(node, wfn, rcov, l, basis)
      INTEGER, INTENT(OUT)                               :: node
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wfn
      REAL(KIND=dp), INTENT(IN)                          :: rcov
      INTEGER, INTENT(IN)                                :: l
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_orbital_nodes', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, m, n
      REAL(KIND=dp)                                      :: ff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: orb

      node = 0
      m = SIZE(basis%bf, 1)
      ALLOCATE (orb(m))
      n = basis%nbas(l)
      orb = 0._dp
      DO i = 1, n
         ff = wfn(i)
         orb(1:m) = orb(1:m)+ff*basis%bf(1:m, i, l)
      END DO
      DO i = 1, m-1
         IF (basis%grid%rad(i) < rcov) THEN
            IF (orb(i)*orb(i+1) < 0._dp) node = node+1
         END IF
      END DO
      DEALLOCATE (orb)

   END SUBROUTINE atom_orbital_nodes
! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param wfn ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE atom_wfnr0(value, wfn, basis)
      REAL(KIND=dp), INTENT(OUT)                         :: value
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wfn
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_wfnr0', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, m, n

      value = 0._dp
      m = MAXVAL(MINLOC(basis%grid%rad))
      n = basis%nbas(0)
      DO i = 1, n
         value = value+wfn(i)*basis%bf(m, i, 0)
      END DO
   END SUBROUTINE atom_wfnr0
! **************************************************************************************************
!> \brief ...
!> \param hmat ...
!> \param umat ...
!> \param orb ...
!> \param ener ...
!> \param nb ...
!> \param nv ...
!> \param maxl ...
! **************************************************************************************************
   SUBROUTINE atom_solve(hmat, umat, orb, ener, nb, nv, maxl)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: hmat, umat
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: orb
      REAL(KIND=dp), DIMENSION(:, 0:), INTENT(INOUT)     :: ener
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: nb, nv
      INTEGER, INTENT(IN)                                :: maxl

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_solve', routineP = moduleN//':'//routineN

      INTEGER                                            :: info, l, lwork, m, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a

      CPASSERT(ALL(nb >= nv))

      orb = 0._dp
      DO l = 0, maxl
         n = nb(l)
         m = nv(l)
         IF (n > 0 .AND. m > 0) THEN
            lwork = 10*m
            ALLOCATE (a(n, n), w(n), work(lwork))
            a(1:m, 1:m) = MATMUL(TRANSPOSE(umat(1:n, 1:m, l)), MATMUL(hmat(1:n, 1:n, l), umat(1:n, 1:m, l)))
            CALL lapack_ssyev("V", "U", m, a(1:m, 1:m), m, w(1:m), work, lwork, info)
            a(1:n, 1:m) = MATMUL(umat(1:n, 1:m, l), a(1:m, 1:m))

            m = MIN(m, SIZE(orb, 2))
            orb(1:n, 1:m, l) = a(1:n, 1:m)
            ener(1:m, l) = w(1:m)

            DEALLOCATE (a, w, work)
         END IF
      END DO

   END SUBROUTINE atom_solve
! **************************************************************************************************
!> \brief ...
!> \param fun ...
!> \param deps ...
!> \retval nc ...
! **************************************************************************************************
   FUNCTION prune_grid(fun, deps) RESULT(nc)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fun
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: deps
      INTEGER                                            :: nc

      INTEGER                                            :: i, nr
      REAL(KIND=dp)                                      :: meps

      meps = 1.e-14_dp
      IF (PRESENT(deps)) meps = deps

      nr = SIZE(fun)
      nc = 0
      DO i = nr, 1, -1
         IF (ABS(fun(i)) > meps) THEN
            nc = i
            EXIT
         END IF
      END DO

   END FUNCTION prune_grid
! **************************************************************************************************
!> \brief ...
!> \param fun ...
!> \param grid ...
!> \retval integral ...
! **************************************************************************************************
   FUNCTION integrate_grid_function1(fun, grid) RESULT(integral)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fun
      TYPE(grid_atom_type)                               :: grid
      REAL(KIND=dp)                                      :: integral

      INTEGER                                            :: nc

      nc = SIZE(fun)
      integral = SUM(fun(1:nc)*grid%wr(1:nc))

   END FUNCTION integrate_grid_function1

! **************************************************************************************************
!> \brief ...
!> \param fun1 ...
!> \param fun2 ...
!> \param grid ...
!> \retval integral ...
! **************************************************************************************************
   FUNCTION integrate_grid_function2(fun1, fun2, grid) RESULT(integral)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fun1, fun2
      TYPE(grid_atom_type)                               :: grid
      REAL(KIND=dp)                                      :: integral

      INTEGER                                            :: nc

      nc = MIN(SIZE(fun1), SIZE(fun2))
      integral = SUM(fun1(1:nc)*fun2(1:nc)*grid%wr(1:nc))

   END FUNCTION integrate_grid_function2

! **************************************************************************************************
!> \brief ...
!> \param fun1 ...
!> \param fun2 ...
!> \param fun3 ...
!> \param grid ...
!> \retval integral ...
! **************************************************************************************************
   FUNCTION integrate_grid_function3(fun1, fun2, fun3, grid) RESULT(integral)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fun1, fun2, fun3
      TYPE(grid_atom_type)                               :: grid
      REAL(KIND=dp)                                      :: integral

      INTEGER                                            :: nc

      nc = MIN(SIZE(fun1), SIZE(fun2), SIZE(fun3))
      integral = SUM(fun1(1:nc)*fun2(1:nc)*fun3(1:nc)*grid%wr(1:nc))

   END FUNCTION integrate_grid_function3
! **************************************************************************************************
!> \brief ...
!> \param cpot ...
!> \param density ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE coulomb_potential_numeric(cpot, density, grid)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: cpot
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: density
      TYPE(grid_atom_type), INTENT(IN)                   :: grid

      CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_potential_numeric', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, nc
      REAL(dp)                                           :: int1, int2
      REAL(dp), DIMENSION(:), POINTER                    :: r, wr

      nc = MIN(SIZE(cpot), SIZE(density))
      r => grid%rad
      wr => grid%wr

      int1 = fourpi*integrate_grid(density, grid)
      int2 = 0._dp
      cpot(nc:) = int1/r(nc:)

      ! test that grid is decreasing
      CPASSERT(r(1) > r(nc))
      DO i = 1, nc
         cpot(i) = int1/r(i)+int2
         int1 = int1-fourpi*density(i)*wr(i)
         int2 = int2+fourpi*density(i)*wr(i)/r(i)
      END DO

   END SUBROUTINE coulomb_potential_numeric

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param cpot ...
!> \param pmat ...
!> \param basis ...
!> \param grid ...
!> \param maxl ...
! **************************************************************************************************
   SUBROUTINE coulomb_potential_analytic(cpot, pmat, basis, grid, maxl)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: cpot
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: pmat
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      TYPE(grid_atom_type)                               :: grid
      INTEGER, INTENT(IN)                                :: maxl

      CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_potential_analytic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, l, m, n
      REAL(KIND=dp)                                      :: a, b, ff, oab, sab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: erfa, expa, z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: unp

      m = SIZE(cpot)
      ALLOCATE (erfa(1:m), expa(1:m), z(1:m))

      cpot = 0._dp

      DO l = 0, maxl
         IF (MAXVAL(ABS(pmat(:, :, l))) < 1.e-14_dp) CYCLE
         SELECT CASE (basis%basis_type)
         CASE DEFAULT
            CPABORT("")
         CASE (GTO_BASIS)
            DO i = 1, basis%nbas(l)
               DO j = i, basis%nbas(l)
                  IF (ABS(pmat(i, j, l)) < 1.e-14_dp) CYCLE
                  ff = pmat(i, j, l)
                  IF (i /= j) ff = 2._dp*ff
                  a = basis%am(i, l)
                  b = basis%am(j, l)
                  sab = SQRT(a+b)
                  oab = rootpi/(a+b)**(l+1.5_dp)*ff
                  z(:) = sab*grid%rad(:)
                  DO k = 1, SIZE(erfa)
                     erfa(k) = oab*erf(z(k))/grid%rad(k)
                  END DO
                  expa(:) = EXP(-z(:)**2)*ff/(a+b)**(l+1)
                  SELECT CASE (l)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     cpot(:) = cpot(:)+0.25_dp*erfa(:)
                  CASE (1)
                     cpot(:) = cpot(:)+0.375_dp*erfa(:)-0.25_dp*expa(:)
                  CASE (2)
                     cpot(:) = cpot(:)+0.9375_dp*erfa(:)-expa(:)*(0.875_dp+0.25_dp*z(:)**2)
                  CASE (3)
                     cpot(:) = cpot(:)+3.28125_dp*erfa(:)-expa(:)*(3.5625_dp+1.375_dp*z(:)**2+0.25*z(:)**4)
                  END SELECT
               END DO
            END DO
         CASE (CGTO_BASIS)
            n = basis%nprim(l)
            m = basis%nbas(l)
            ALLOCATE (unp(n, n))

            unp(1:n, 1:n) = MATMUL(MATMUL(basis%cm(1:n, 1:m, l), pmat(1:m, 1:m, l)), &
                                   TRANSPOSE(basis%cm(1:n, 1:m, l)))
            DO i = 1, basis%nprim(l)
               DO j = i, basis%nprim(l)
                  IF (ABS(unp(i, j)) < 1.e-14_dp) CYCLE
                  ff = unp(i, j)
                  IF (i /= j) ff = 2._dp*ff
                  a = basis%am(i, l)
                  b = basis%am(j, l)
                  sab = SQRT(a+b)
                  oab = rootpi/(a+b)**(l+1.5_dp)*ff
                  z(:) = sab*grid%rad(:)
                  DO k = 1, SIZE(erfa)
                     erfa(k) = oab*erf(z(k))/grid%rad(k)
                  END DO
                  expa(:) = EXP(-z(:)**2)*ff/(a+b)**(l+1)
                  SELECT CASE (l)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     cpot(:) = cpot(:)+0.25_dp*erfa(:)
                  CASE (1)
                     cpot(:) = cpot(:)+0.375_dp*erfa(:)-0.25_dp*expa(:)
                  CASE (2)
                     cpot(:) = cpot(:)+0.9375_dp*erfa(:)-expa(:)*(0.875_dp+0.25_dp*z(:)**2)
                  CASE (3)
                     cpot(:) = cpot(:)+3.28125_dp*erfa(:)-expa(:)*(3.5625_dp+1.375_dp*z(:)**2+0.25*z(:)**4)
                  END SELECT
               END DO
            END DO

            DEALLOCATE (unp)
         END SELECT
      END DO
      DEALLOCATE (erfa, expa, z)

   END SUBROUTINE coulomb_potential_analytic

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param kmat ...
!> \param state ...
!> \param occ ...
!> \param wfn ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE exchange_numeric(kmat, state, occ, wfn, basis)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: kmat
      TYPE(atom_state), INTENT(IN)                       :: state
      REAL(KIND=dp), DIMENSION(0:, :), INTENT(IN)        :: occ
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: wfn
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'exchange_numeric', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ia, ib, k, lad, lbc, lh, ll, nbas, &
                                                            norb, nr, nu
      REAL(KIND=dp)                                      :: almn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpot, nai, nbi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: orb
      REAL(KIND=dp), DIMENSION(0:maxfac)                 :: arho

      arho = 0._dp
      DO ll = 0, maxfac, 2
         lh = ll/2
         arho(ll) = fac(ll)/fac(lh)**2
      END DO

      kmat = 0._dp

      nr = basis%grid%nr
      ALLOCATE (nai(nr), nbi(nr), cpot(nr))

      DO lad = 0, state%maxl_calc
         DO lbc = 0, state%maxl_occ
            norb = state%maxn_occ(lbc)
            nbas = basis%nbas(lbc)
            ! calculate orbitals for angmom lbc
            ALLOCATE (orb(nr, norb))
            orb = 0._dp
            DO i = 1, norb
               DO k = 1, nbas
                  orb(:, i) = orb(:, i)+wfn(k, i, lbc)*basis%bf(:, k, lbc)
               END DO
            END DO
            DO nu = ABS(lad-lbc), lad+lbc, 2
               almn = arho(-lad+lbc+nu)*arho(lad-lbc+nu)*arho(lad+lbc-nu)/(REAL(lad+lbc+nu+1, dp)*arho(lad+lbc+nu))
               almn = -0.5_dp*almn

               DO ia = 1, basis%nbas(lad)
                  DO i = 1, norb
                     nai(:) = orb(:, i)*basis%bf(:, ia, lad)
                     CALL potential_numeric(cpot, nai, nu, basis%grid)
                     DO ib = 1, basis%nbas(lad)
                        kmat(ia, ib, lad) = kmat(ia, ib, lad)+almn*occ(lbc, i)* &
                                            integrate_grid(cpot, orb(:, i), basis%bf(:, ib, lad), basis%grid)
                     END DO
                  END DO
               END DO

            END DO
            DEALLOCATE (orb)
         END DO
      END DO

      DEALLOCATE (nai, nbi, cpot)

   END SUBROUTINE exchange_numeric

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param kmat ...
!> \param state ...
!> \param occ ...
!> \param wfn ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE exchange_semi_analytic(kmat, state, occ, wfn, basis)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: kmat
      TYPE(atom_state), INTENT(IN)                       :: state
      REAL(KIND=dp), DIMENSION(0:, :), INTENT(IN)        :: occ
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: wfn
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'exchange_semi_analytic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ia, ib, k, lad, lbc, lh, ll, nbas, &
                                                            norb, nr, nu
      REAL(KIND=dp)                                      :: almn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpot, nai, nbi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: orb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: pot
      REAL(KIND=dp), DIMENSION(0:maxfac)                 :: arho

      arho = 0._dp
      DO ll = 0, maxfac, 2
         lh = ll/2
         arho(ll) = fac(ll)/fac(lh)**2
      END DO

      kmat = 0._dp

      nr = basis%grid%nr
      nbas = MAXVAL(basis%nbas)
      ALLOCATE (pot(nr, nbas, nbas))
      ALLOCATE (nai(nr), nbi(nr), cpot(nr))

      DO lad = 0, state%maxl_calc
         DO lbc = 0, state%maxl_occ
            norb = state%maxn_occ(lbc)
            nbas = basis%nbas(lbc)
            ! calculate orbitals for angmom lbc
            ALLOCATE (orb(nr, norb))
            orb = 0._dp
            DO i = 1, norb
               DO k = 1, nbas
                  orb(:, i) = orb(:, i)+wfn(k, i, lbc)*basis%bf(:, k, lbc)
               END DO
            END DO
            DO nu = ABS(lad-lbc), lad+lbc, 2
               almn = arho(-lad+lbc+nu)*arho(lad-lbc+nu)*arho(lad+lbc-nu)/(REAL(lad+lbc+nu+1, dp)*arho(lad+lbc+nu))
               almn = -0.5_dp*almn
               ! calculate potential for basis function pair (lad,lbc)
               pot = 0._dp
               CALL potential_analytic(pot, lad, lbc, nu, basis)
               DO ia = 1, basis%nbas(lad)
                  DO i = 1, norb
                     cpot = 0._dp
                     DO k = 1, nbas
                        cpot(:) = cpot(:)+pot(:, ia, k)*wfn(k, i, lbc)
                     END DO
                     DO ib = 1, basis%nbas(lad)
                        kmat(ia, ib, lad) = kmat(ia, ib, lad)+almn*occ(lbc, i)* &
                                            integrate_grid(cpot, orb(:, i), basis%bf(:, ib, lad), basis%grid)
                     END DO
                  END DO
               END DO
            END DO
            DEALLOCATE (orb)
         END DO
      END DO

      DEALLOCATE (pot)
      DEALLOCATE (nai, nbi, cpot)

   END SUBROUTINE exchange_semi_analytic

! **************************************************************************************************
!> \brief ...
!> \param cpot ...
!> \param density ...
!> \param nu ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE potential_numeric(cpot, density, nu, grid)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: cpot
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: density
      INTEGER, INTENT(IN)                                :: nu
      TYPE(grid_atom_type), INTENT(IN)                   :: grid

      CHARACTER(len=*), PARAMETER :: routineN = 'potential_numeric', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, nc
      REAL(dp)                                           :: int1, int2
      REAL(dp), DIMENSION(:), POINTER                    :: r, wr

      nc = MIN(SIZE(cpot), SIZE(density))
      r => grid%rad
      wr => grid%wr

      int1 = integrate_grid(density, r**nu, grid)
      int2 = 0._dp
      cpot(nc:) = int1/r(nc:)**(nu+1)

      ! test that grid is decreasing
      CPASSERT(r(1) > r(nc))
      DO i = 1, nc
         cpot(i) = int1/r(i)**(nu+1)+int2*r(i)**nu
         int1 = int1-r(i)**(nu)*density(i)*wr(i)
         int2 = int2+density(i)*wr(i)/r(i)**(nu+1)
      END DO

   END SUBROUTINE potential_numeric
! **************************************************************************************************
!> \brief ...
!> \param cpot ...
!> \param la ...
!> \param lb ...
!> \param nu ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE potential_analytic(cpot, la, lb, nu, basis)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: cpot
      INTEGER, INTENT(IN)                                :: la, lb, nu
      TYPE(atom_basis_type), INTENT(IN)                  :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'potential_analytic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, l, ll, m
      REAL(KIND=dp)                                      :: a, b, oab, sab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: erfa, expa, z

      m = SIZE(cpot, 1)
      ALLOCATE (erfa(1:m), expa(1:m), z(1:m))

      ll = la+lb

      cpot = 0._dp

      SELECT CASE (basis%basis_type)
      CASE DEFAULT
         CPABORT("")
      CASE (GTO_BASIS)
         DO i = 1, basis%nbas(la)
            DO j = 1, basis%nbas(lb)
               a = basis%am(i, la)
               b = basis%am(j, lb)
               sab = SQRT(a+b)
               oab = dfac(ll+nu+1)*rootpi/(sab**(ll+2))/(2._dp**((ll+nu)/2+2))
               z(:) = sab*basis%grid%rad(:)
               DO k = 1, SIZE(erfa)
                  erfa(k) = oab*erf(z(k))/z(k)**(nu+1)
               END DO
               cpot(:, i, j) = erfa(:)
               expa(:) = EXP(-z(:)**2)/(sab**(ll+2))/(2._dp**((ll+nu)/2+2))
               SELECT CASE (ll)
               CASE DEFAULT
                  CPABORT("")
               CASE (0)
                  CPASSERT(nu == 0)
               CASE (1)
                  CPASSERT(nu == 1)
                  cpot(:, i, j) = cpot(:, i, j)-6._dp*expa(:)/z(:)
               CASE (2)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     cpot(:, i, j) = cpot(:, i, j)-2._dp*expa(:)
                  CASE (2)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(20._dp+30._dp/z(:)**2)
                  END SELECT
               CASE (3)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (1)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(12._dp*z(:)+30._dp/z(:))
                  CASE (3)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(56._dp*z(:)+140._dp/z(:)+210._dp*z(:)**3)
                  END SELECT
               CASE (4)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(4._dp*z(:)**2+14._dp)
                  CASE (2)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(40._dp*z(:)**2+140._dp+210._dp/z(:)**2)
                  CASE (4)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(144._dp*z(:)**2+504._dp+1260._dp/z(:)**2+1890._dp/z(:)**4)
                  END SELECT
               CASE (5)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (1)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(24._dp*z(:)**3+108._dp*z(:)+210._dp/z(:))
                  CASE (3)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(112._dp*z(:)**3+504._dp*z(:)+1260._dp/z(:)+1890._dp/z(:)**3)
                  CASE (5)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(352._dp*z(:)**3+1584._dp*z(:)+5544._dp/z(:)+ &
                                                            13860._dp/z(:)**3+20790._dp/z(:)**5)
                  END SELECT
               CASE (6)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(8._dp*z(:)**4+44._dp*z(:)**2+114._dp)
                  CASE (2)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(80._dp*z(:)**4+440._dp*z(:)**2+1260._dp+1896._dp/z(:)**2)
                  CASE (4)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(288._dp*z(:)**4+1584._dp*z(:)**2+5544._dp+ &
                                                            13860._dp/z(:)**2+20790._dp/z(:)**4)
                  CASE (6)
                     cpot(:, i, j) = cpot(:, i, j)-expa(:)*(832._dp*z(:)**4+4576._dp*z(:)**2+20592._dp+ &
                                                            72072._dp/z(:)**2+180180._dp/z(:)**4+270270._dp/z(:)**6)
                  END SELECT
               END SELECT
            END DO
         END DO
      CASE (CGTO_BASIS)
         DO i = 1, basis%nprim(la)
            DO j = 1, basis%nprim(lb)
               a = basis%am(i, la)
               b = basis%am(j, lb)
               sab = SQRT(a+b)
               oab = dfac(ll+nu+1)*rootpi/sab**(ll+2)/2._dp**((ll+nu)/2+2)
               z(:) = sab*basis%grid%rad(:)
               DO k = 1, SIZE(erfa)
                  erfa(k) = oab*erf(z(k))/z(k)**(nu+1)
               END DO
               expa(:) = EXP(-z(:)**2)/sab**(ll+2)/2._dp**((ll+nu)/2+2)
               SELECT CASE (ll)
               CASE DEFAULT
                  CPABORT("")
               CASE (0)
                  CPASSERT(nu == 0)
               CASE (1)
                  CPASSERT(nu == 1)
                  erfa(:) = erfa(:)-6._dp*expa(:)/z(:)
               CASE (2)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     erfa(:) = erfa(:)-2._dp*expa(:)
                  CASE (2)
                     erfa(:) = erfa(:)-expa(:)*(20._dp+30._dp/z(:)**2)
                  END SELECT
               CASE (3)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (1)
                     erfa(:) = erfa(:)-expa(:)*(12._dp*z(:)+30._dp/z(:))
                  CASE (3)
                     erfa(:) = erfa(:)-expa(:)*(56._dp*z(:)+140._dp/z(:)+210._dp*z(:)**3)
                  END SELECT
               CASE (4)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     erfa(:) = erfa(:)-expa(:)*(4._dp*z(:)**2+14._dp)
                  CASE (2)
                     erfa(:) = erfa(:)-expa(:)*(40._dp*z(:)**2+140._dp+210._dp/z(:)**2)
                  CASE (4)
                     erfa(:) = erfa(:)-expa(:)*(144._dp*z(:)**2+504._dp+1260._dp/z(:)**2+1890._dp/z(:)**4)
                  END SELECT
               CASE (5)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (1)
                     erfa(:) = erfa(:)-expa(:)*(24._dp*z(:)**3+108._dp*z(:)+210._dp/z(:))
                  CASE (3)
                     erfa(:) = erfa(:)-expa(:)*(112._dp*z(:)**3+504._dp*z(:)+1260._dp/z(:)+1890._dp/z(:)**3)
                  CASE (5)
                     erfa(:) = erfa(:)-expa(:)*(352._dp*z(:)**3+1584._dp*z(:)+5544._dp/z(:)+ &
                                                13860._dp/z(:)**3+20790._dp/z(:)**5)
                  END SELECT
               CASE (6)
                  SELECT CASE (nu)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (0)
                     erfa(:) = erfa(:)-expa(:)*(8._dp*z(:)**4+44._dp*z(:)**2+114._dp)
                  CASE (2)
                     erfa(:) = erfa(:)-expa(:)*(80._dp*z(:)**4+440._dp*z(:)**2+1260._dp+1896._dp/z(:)**2)
                  CASE (4)
                     erfa(:) = erfa(:)-expa(:)*(288._dp*z(:)**4+1584._dp*z(:)**2+5544._dp+ &
                                                13860._dp/z(:)**2+20790._dp/z(:)**4)
                  CASE (6)
                     erfa(:) = erfa(:)-expa(:)*(832._dp*z(:)**4+4576._dp*z(:)**2+20592._dp+ &
                                                72072._dp/z(:)**2+180180._dp/z(:)**4+270270._dp/z(:)**6)
                  END SELECT
               END SELECT
               DO k = 1, basis%nbas(la)
                  DO l = 1, basis%nbas(lb)
                     cpot(:, k, l) = cpot(:, k, l)+erfa(:)*basis%cm(i, k, la)*basis%cm(j, l, lb)
                  END DO
               END DO
            END DO
         END DO

      END SELECT

      DEALLOCATE (erfa, expa, z)

   END SUBROUTINE potential_analytic
! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param opmat ...
!> \param pmat ...
!> \retval trace ...
! **************************************************************************************************
   FUNCTION atom_trace(opmat, pmat) RESULT(trace)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: opmat, pmat
      REAL(KIND=dp)                                      :: trace

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_trace', routineP = moduleN//':'//routineN

      trace = SUM(opmat(:, :, :)*pmat(:, :, :))

   END FUNCTION atom_trace

! **************************************************************************************************
!> \brief ...
!> \param imat ...
!> \param cpot ...
!> \param basis ...
!> \param derivatives ...
! **************************************************************************************************
   SUBROUTINE numpot_matrix(imat, cpot, basis, derivatives)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: imat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cpot
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis
      INTEGER, INTENT(IN)                                :: derivatives

      CHARACTER(len=*), PARAMETER :: routineN = 'numpot_matrix', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, l, n

      SELECT CASE (derivatives)
      CASE (0)
         DO l = 0, 3
            n = basis%nbas(l)
            DO i = 1, n
               DO j = i, n
                  imat(i, j, l) = imat(i, j, l)+ &
                                  integrate_grid(cpot, basis%bf(:, i, l), basis%bf(:, j, l), basis%grid)
                  imat(j, i, l) = imat(i, j, l)
               END DO
            END DO
         END DO
      CASE (1)
         DO l = 0, 3
            n = basis%nbas(l)
            DO i = 1, n
               DO j = i, n
                  imat(i, j, l) = imat(i, j, l)+ &
                                  integrate_grid(cpot, basis%dbf(:, i, l), basis%bf(:, j, l), basis%grid)
                  imat(i, j, l) = imat(i, j, l)+ &
                                  integrate_grid(cpot, basis%bf(:, i, l), basis%dbf(:, j, l), basis%grid)
                  imat(j, i, l) = imat(i, j, l)
               END DO
            END DO
         END DO
      CASE (2)
         DO l = 0, 3
            n = basis%nbas(l)
            DO i = 1, n
               DO j = i, n
                  imat(i, j, l) = imat(i, j, l)+ &
                                  integrate_grid(cpot, basis%dbf(:, i, l), basis%dbf(:, j, l), basis%grid)
                  imat(j, i, l) = imat(i, j, l)
               END DO
            END DO
         END DO
      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE numpot_matrix
! **************************************************************************************************
!> \brief ...
!> \param jmat ...
!> \param erint ...
!> \param pmat ...
!> \param nsize ...
!> \param all_nu ...
! **************************************************************************************************
   SUBROUTINE ceri_contract(jmat, erint, pmat, nsize, all_nu)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: jmat
      TYPE(eri), DIMENSION(:), INTENT(IN)                :: erint
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: pmat
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: nsize
      LOGICAL, INTENT(IN), OPTIONAL                      :: all_nu

      CHARACTER(len=*), PARAMETER :: routineN = 'ceri_contract', routineP = moduleN//':'//routineN

      INTEGER                                            :: i1, i2, ij1, ij2, j1, j2, l1, l2, ll, &
                                                            n1, n2
      LOGICAL                                            :: have_all_nu
      REAL(KIND=dp)                                      :: eint, f1, f2

      IF (PRESENT(all_nu)) THEN
         have_all_nu = all_nu
      ELSE
         have_all_nu = .FALSE.
      END IF

      jmat(:, :, :) = 0._dp
      ll = 0
      DO l1 = 0, 3
         n1 = nsize(l1)
         DO l2 = 0, l1
            n2 = nsize(l2)
            ll = ll+1
            ij1 = 0
            DO i1 = 1, n1
               DO j1 = i1, n1
                  ij1 = ij1+1
                  f1 = 1._dp
                  IF (i1 /= j1) f1 = 2._dp
                  ij2 = 0
                  DO i2 = 1, n2
                     DO j2 = i2, n2
                        ij2 = ij2+1
                        f2 = 1._dp
                        IF (i2 /= j2) f2 = 2._dp
                        eint = erint(ll)%int(ij1, ij2)
                        IF (l1 == l2) THEN
                           jmat(i1, j1, l1) = jmat(i1, j1, l1)+f2*pmat(i2, j2, l2)*eint
                        ELSE
                           jmat(i1, j1, l1) = jmat(i1, j1, l1)+f2*pmat(i2, j2, l2)*eint
                           jmat(i2, j2, l2) = jmat(i2, j2, l2)+f1*pmat(i1, j1, l1)*eint
                        END IF
                     END DO
                  END DO
               END DO
            END DO
            IF (have_all_nu) THEN
               ! skip integral blocks with nu/=0
               ll = ll+l2
            END IF
         END DO
      END DO
      DO l1 = 0, 3
         n1 = nsize(l1)
         DO i1 = 1, n1
            DO j1 = i1, n1
               jmat(j1, i1, l1) = jmat(i1, j1, l1)
            END DO
         END DO
      END DO

   END SUBROUTINE ceri_contract
! **************************************************************************************************
!> \brief ...
!> \param kmat ...
!> \param erint ...
!> \param pmat ...
!> \param nsize ...
! **************************************************************************************************
   SUBROUTINE eeri_contract(kmat, erint, pmat, nsize)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(INOUT)  :: kmat
      TYPE(eri), DIMENSION(:), INTENT(IN)                :: erint
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: pmat
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: nsize

      CHARACTER(len=*), PARAMETER :: routineN = 'eeri_contract', routineP = moduleN//':'//routineN

      INTEGER                                            :: i1, i2, ij1, ij2, j1, j2, l1, l2, lh, &
                                                            ll, n1, n2, nu
      REAL(KIND=dp)                                      :: almn, eint, f1, f2
      REAL(KIND=dp), DIMENSION(0:maxfac)                 :: arho

      arho = 0._dp
      DO ll = 0, maxfac, 2
         lh = ll/2
         arho(ll) = fac(ll)/fac(lh)**2
      END DO

      kmat(:, :, :) = 0._dp
      ll = 0
      DO l1 = 0, 3
         n1 = nsize(l1)
         DO l2 = 0, l1
            n2 = nsize(l2)
            DO nu = ABS(l1-l2), l1+l2, 2
               almn = arho(-l1+l2+nu)*arho(l1-l2+nu)*arho(l1+l2-nu)/(REAL(l1+l2+nu+1, dp)*arho(l1+l2+nu))
               almn = -0.5_dp*almn
               ll = ll+1
               ij1 = 0
               DO i1 = 1, n1
                  DO j1 = i1, n1
                     ij1 = ij1+1
                     f1 = 1._dp
                     IF (i1 /= j1) f1 = 2._dp
                     ij2 = 0
                     DO i2 = 1, n2
                        DO j2 = i2, n2
                           ij2 = ij2+1
                           f2 = 1._dp
                           IF (i2 /= j2) f2 = 2._dp
                           eint = erint(ll)%int(ij1, ij2)
                           IF (l1 == l2) THEN
                              kmat(i1, j1, l1) = kmat(i1, j1, l1)+f2*almn*pmat(i2, j2, l2)*eint
                           ELSE
                              kmat(i1, j1, l1) = kmat(i1, j1, l1)+f2*almn*pmat(i2, j2, l2)*eint
                              kmat(i2, j2, l2) = kmat(i2, j2, l2)+f1*almn*pmat(i1, j1, l1)*eint
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      DO l1 = 0, 3
         n1 = nsize(l1)
         DO i1 = 1, n1
            DO j1 = i1, n1
               kmat(j1, i1, l1) = kmat(i1, j1, l1)
            END DO
         END DO
      END DO

   END SUBROUTINE eeri_contract
! **************************************************************************************************
!> \brief ...
!> \param emat ...
!> \param demax ...
!> \param kmat ...
!> \param pmat ...
!> \param umat ...
!> \param upmat ...
!> \param nval ...
!> \param nbs ...
! **************************************************************************************************
   SUBROUTINE err_matrix(emat, demax, kmat, pmat, umat, upmat, nval, nbs)
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(OUT)    :: emat
      REAL(KIND=dp), INTENT(OUT)                         :: demax
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)     :: kmat, pmat, umat, upmat
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: nval, nbs

      CHARACTER(len=*), PARAMETER :: routineN = 'err_matrix', routineP = moduleN//':'//routineN

      INTEGER                                            :: l, m, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tkmat, tpmat

      emat = 0._dp
      DO l = 0, 3
         n = nval(l)
         m = nbs(l)
         IF (m > 0) THEN
            ALLOCATE (tkmat(1:m, 1:m), tpmat(1:m, 1:m))
            tkmat = 0._dp
            tpmat = 0._dp
            tkmat(1:m, 1:m) = MATMUL(TRANSPOSE(umat(1:n, 1:m, l)), MATMUL(kmat(1:n, 1:n, l), umat(1:n, 1:m, l)))
            tpmat(1:m, 1:m) = MATMUL(TRANSPOSE(umat(1:n, 1:m, l)), MATMUL(pmat(1:n, 1:n, l), umat(1:n, 1:m, l)))
            tpmat(1:m, 1:m) = MATMUL(upmat(1:m, 1:m, l), MATMUL(tpmat(1:m, 1:m), upmat(1:m, 1:m, l)))

            emat(1:m, 1:m, l) = MATMUL(tkmat(1:m, 1:m), tpmat(1:m, 1:m))-MATMUL(tpmat(1:m, 1:m), tkmat(1:m, 1:m))

            DEALLOCATE (tkmat, tpmat)
         END IF
      END DO
      demax = MAXVAL(ABS(emat))

   END SUBROUTINE err_matrix
! **************************************************************************************************
!> \brief ...
!> \param density1 ...
!> \param density2 ...
!> \param zcore ...
!> \param state ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE slater_density(density1, density2, zcore, state, grid)
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: density1, density2
      INTEGER, INTENT(IN)                                :: zcore
      TYPE(atom_state), INTENT(IN)                       :: state
      TYPE(grid_atom_type), INTENT(IN)                   :: grid

      CHARACTER(len=*), PARAMETER :: routineN = 'slater_density', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, l, n
      INTEGER, DIMENSION(4, 7)                           :: ne
      REAL(KIND=dp)                                      :: a, pf

      density1 = 0._dp
      density2 = 0._dp
      ne(1:4, 1:7) = NINT(state%core(0:3, 1:7)+state%occ(0:3, 1:7))
      DO l = 0, state%maxl_occ
         DO i = 1, SIZE(state%occ, 2)
            IF (state%occ(l, i) > 0._dp) THEN
               n = i+l
               a = srules(zcore, ne, n, l)
               pf = 1._dp/SQRT(fac(2*n))*(2._dp*a)**(n+0.5_dp)
               IF (state%multiplicity == -1) THEN
                  density1(:) = density1(:)+state%occ(l, i)/fourpi*(grid%rad(:)**(n-1)*EXP(-a*grid%rad(:))*pf)**2
               ELSE
                  density1(:) = density1(:)+state%occa(l, i)/fourpi*(grid%rad(:)**(n-1)*EXP(-a*grid%rad(:))*pf)**2
                  density2(:) = density2(:)+state%occb(l, i)/fourpi*(grid%rad(:)**(n-1)*EXP(-a*grid%rad(:))*pf)**2
               END IF
            END IF
         END DO
      END DO

   END SUBROUTINE slater_density
! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param vxc ...
! **************************************************************************************************
   SUBROUTINE wigner_slater_functional(rho, vxc)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rho
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: vxc

      CHARACTER(len=*), PARAMETER :: routineN = 'wigner_slater_functional', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: ec, ex, rs, vc, vx

      vxc = 0._dp
      DO i = 1, SIZE(rho)
         IF (rho(i) > 1.e-20_dp) THEN
            ex = -0.7385588_dp*rho(i)**0.333333333_dp
            vx = 1.333333333_dp*ex
            rs = (3._dp/fourpi/rho(i))**0.333333333_dp
            ec = -0.88_dp/(rs+7.8_dp)
            vc = ec*(1._dp+rs/(3._dp*(rs+7.8_dp)))
            vxc(i) = vx+vc
         END IF
      END DO

   END SUBROUTINE wigner_slater_functional
! **************************************************************************************************
!> \brief ...
!> \param method ...
!> \param multiplicity ...
!> \retval consistent ...
! **************************************************************************************************
   FUNCTION atom_consistent_method(method, multiplicity) RESULT(consistent)
      INTEGER, INTENT(IN)                                :: method, multiplicity
      LOGICAL                                            :: consistent

      SELECT CASE (method)
      CASE DEFAULT
         consistent = .FALSE.
      CASE (do_rks_atom)
         consistent = (multiplicity == -1)
      CASE (do_uks_atom)
         consistent = (multiplicity /= -1)
      CASE (do_rhf_atom)
         consistent = (multiplicity == -1)
      CASE (do_uhf_atom)
         consistent = (multiplicity /= -1)
      CASE (do_rohf_atom)
         consistent = .FALSE.
      END SELECT

   END FUNCTION atom_consistent_method
! **************************************************************************************************
!> \brief ...
!> \param atom ...
!> \param rho0 ...
!> \retval nc ...
! **************************************************************************************************
   SUBROUTINE get_rho0(atom, rho0)
      TYPE(atom_type), INTENT(IN)                        :: atom
      REAL(KIND=dp), INTENT(OUT)                         :: rho0

      INTEGER                                            :: m0, m1, m2, method, nr
      LOGICAL                                            :: nlcc, spinpol
      REAL(KIND=dp)                                      :: d0, d1, d2, r0, r1, r2, w0, w1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xfun
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rho

      method = atom%method_type
      SELECT CASE (method)
      CASE (do_rks_atom, do_rhf_atom)
         spinpol = .FALSE.
      CASE (do_uks_atom, do_uhf_atom)
         spinpol = .TRUE.
      CASE (do_rohf_atom)
         CPABORT("")
      CASE DEFAULT
         CPABORT("")
      END SELECT

      nr = atom%basis%grid%nr
      nlcc = .FALSE.
      IF (atom%potential%ppot_type == gth_pseudo) THEN
         nlcc = atom%potential%gth_pot%nlcc
      ELSE IF (atom%potential%ppot_type == upf_pseudo) THEN
         nlcc = atom%potential%upf_pot%core_correction
      END IF
      IF (nlcc) THEN
         ALLOCATE (xfun(nr))
      END IF

      m0 = MAXVAL(MINLOC(atom%basis%grid%rad))
      IF (m0 == nr) THEN
         m1 = m0-1
         m2 = m0-2
      ELSE IF (m0 == 1) THEN
         m1 = 2
         m2 = 3
      ELSE
         CPABORT("GRID Definition incompatible")
      END IF
      r0 = atom%basis%grid%rad(m0)
      r1 = atom%basis%grid%rad(m1)
      r2 = atom%basis%grid%rad(m2)
      w0 = r1/(r1-r0)
      w1 = 1-w0

      IF (spinpol) THEN
         ALLOCATE (rho(nr, 2))
         CALL atom_density(rho(:, 1), atom%orbitals%pmata, atom%basis, atom%state%maxl_occ, typ="RHO")
         CALL atom_density(rho(:, 2), atom%orbitals%pmatb, atom%basis, atom%state%maxl_occ, typ="RHO")
         IF (nlcc) THEN
            xfun = 0.0_dp
            CALL atom_core_density(xfun(:), atom%potential, typ="RHO", rr=atom%basis%grid%rad)
            rho(:, 1) = rho(:, 1)+0.5_dp*xfun(:)
            rho(:, 2) = rho(:, 2)+0.5_dp*xfun(:)
         END IF
         rho(:, 1) = rho(:, 1)+rho(:, 2)
      ELSE
         ALLOCATE (rho(nr, 1))
         CALL atom_density(rho(:, 1), atom%orbitals%pmat, atom%basis, atom%state%maxl_occ, typ="RHO")
         IF (nlcc) THEN
            CALL atom_core_density(rho(:, 1), atom%potential, typ="RHO", rr=atom%basis%grid%rad)
         END IF
      END IF
      d0 = rho(m0, 1)
      d1 = rho(m1, 1)
      d2 = rho(m2, 1)

      rho0 = w0*d0+w1*d1
      rho0 = MAX(rho0, 0.0_dp)

      DEALLOCATE (rho)
      IF (nlcc) THEN
         DEALLOCATE (xfun)
      END IF

   END SUBROUTINE get_rho0

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param crad ...
!> \param iw    ...
! **************************************************************************************************
   SUBROUTINE atom_condnumber(basis, crad, iw)
      TYPE(atom_basis_type), POINTER                     :: basis
      REAL(KIND=dp)                                      :: crad
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: ci
      REAL(KIND=dp), DIMENSION(10)                       :: cnum, rad

      WRITE (iw, '(/,A,F8.4)') " Basis Set Condition Numbers: 2*covalent radius=", 2*crad
      CALL init_orbital_pointers(3)
      CALL init_spherical_harmonics(3, 0)
      cnum = 0.0_dp
      DO i = 1, 9
         ci = 2.0_dp*(0.85_dp+i*0.05_dp)
         rad(i) = crad*ci
         CALL atom_basis_condnum(basis, rad(i), cnum(i))
         WRITE (iw, '(A,F15.3,T50,A,F14.4)') " Lattice constant:", &
            rad(i), "Condition number:", cnum(i)
      END DO
      rad(10) = 0.01_dp
      CALL atom_basis_condnum(basis, rad(10), cnum(10))
      WRITE (iw, '(A,A,T50,A,F14.4)') " Lattice constant:", &
         "            Inf", "Condition number:", cnum(i)
      CALL deallocate_orbital_pointers
      CALL deallocate_spherical_harmonics

   END SUBROUTINE atom_condnumber
! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param zv ...
!> \param iw    ...
! **************************************************************************************************
   SUBROUTINE atom_completeness(basis, zv, iw)
      TYPE(atom_basis_type), POINTER                     :: basis
      INTEGER, INTENT(IN)                                :: zv, iw

      INTEGER                                            :: i, j, l, ll, m, n, nbas, nl, nr
      INTEGER, DIMENSION(0:3)                            :: nelem, nlmax, nlmin
      INTEGER, DIMENSION(4, 7)                           :: ne
      REAL(KIND=dp)                                      :: al, c1, c2, pf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sfun
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: bmat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: omat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: sint
      REAL(KIND=dp), DIMENSION(0:3, 10)                  :: snl
      REAL(KIND=dp), DIMENSION(2)                        :: sse
      REAL(KIND=dp), DIMENSION(4, 7)                     :: sexp

      ne = 0
      nelem(0:3) = ptable(zv)%e_conv(0:3)
      DO l = 0, 3
         ll = 2*(2*l+1)
         DO i = 1, 7
            IF (nelem(l) >= ll) THEN
               ne(l+1, i) = ll
               nelem(l) = nelem(l)-ll
            ELSE IF (nelem(l) > 0) THEN
               ne(l+1, i) = nelem(l)
               nelem(l) = 0
            ELSE
               EXIT
            END IF
         END DO
      END DO

      nlmin = 1
      nlmax = 1
      DO l = 0, 3
         nlmin(l) = l+1
         DO i = 1, 7
            IF (ne(l+1, i) > 0) THEN
               nlmax(l) = i+l
            END IF
         END DO
         nlmax(l) = MAX(nlmax(l), nlmin(l)+1)
      END DO

      ! Slater exponents
      sexp = 0.0_dp
      DO l = 0, 3
         sse(1) = 0.05_dp
         sse(2) = 10.0_dp
         DO i = l+1, 7
            sexp(l+1, i) = srules(zv, ne, i, l)
            IF (ne(l+1, i-l) > 0) THEN
               sse(1) = MAX(sse(1), sexp(l+1, i))
               sse(2) = MIN(sse(2), sexp(l+1, i))
            END IF
         END DO
         DO i = 1, 10
            snl(l, i) = ABS(2._dp*sse(1)-0.5_dp*sse(2))/9._dp*REAL(i-1, KIND=dp)+0.5_dp*MIN(sse(1), sse(2))
         END DO
      END DO

      nbas = MAXVAL(basis%nbas)
      ALLOCATE (omat(nbas, nbas, 0:3))
      nr = SIZE(basis%bf, 1)
      ALLOCATE (sfun(nr), sint(10, 2, nbas, 0:3))
      sint = 0._dp

      ! calculate overlaps between test functions and basis
      DO l = 0, 3
         DO i = 1, 10
            al = snl(l, i)
            nl = nlmin(l)
            pf = (2._dp*al)**nl*SQRT(2._dp*al/fac(2*nl))
            sfun(1:nr) = pf*basis%grid%rad(1:nr)**(nl-1)*EXP(-al*basis%grid%rad(1:nr))
            DO j = 1, basis%nbas(l)
               sint(i, 1, j, l) = SUM(sfun(1:nr)*basis%bf(1:nr, j, l)*basis%grid%wr(1:nr))
            END DO
            nl = nlmax(l)
            pf = (2._dp*al)**nl*SQRT(2._dp*al/fac(2*nl))
            sfun(1:nr) = pf*basis%grid%rad(1:nr)**(nl-1)*EXP(-al*basis%grid%rad(1:nr))
            DO j = 1, basis%nbas(l)
               sint(i, 2, j, l) = SUM(sfun(1:nr)*basis%bf(1:nr, j, l)*basis%grid%wr(1:nr))
            END DO
         END DO
      END DO

      DO l = 0, 3
         n = basis%nbas(l)
         IF (n <= 0) CYCLE
         m = basis%nprim(l)
         SELECT CASE (basis%basis_type)
         CASE DEFAULT
            CPABORT("")
         CASE (GTO_BASIS)
            CALL sg_overlap(omat(1:n, 1:n, l), l, basis%am(1:n, l), basis%am(1:n, l))
         CASE (CGTO_BASIS)
            ALLOCATE (bmat(m, m))
            CALL sg_overlap(bmat(1:m, 1:m), l, basis%am(1:m, l), basis%am(1:m, l))
            CALL contract2(omat(1:n, 1:n, l), bmat(1:m, 1:m), basis%cm(1:m, 1:n, l))
            DEALLOCATE (bmat)
         CASE (STO_BASIS)
            CALL sto_overlap(omat(1:n, 1:n, l), basis%ns(1:n, l), basis%as(1:n, l), &
                             basis%ns(1:n, l), basis%as(1:n, l))
         CASE (NUM_BASIS)
            CPABORT("")
         END SELECT
         CALL invmat_symm(omat(1:n, 1:n, l))
      END DO

      WRITE (iw, '(/,A)') " Basis Set Completeness Estimates"
      DO l = 0, 3
         n = basis%nbas(l)
         IF (n <= 0) CYCLE
         WRITE (iw, '(A,I3)') "   L-quantum number: ", l
         WRITE (iw, '(A,T31,A,I2,T61,A,I2)') "    Slater Exponent", "Completeness n-qm=", nlmin(l), &
            "Completeness n-qm=", nlmax(l)
         DO i = 10, 1, -1
            c1 = DOT_PRODUCT(sint(i, 1, 1:n, l), MATMUL(omat(1:n, 1:n, l), sint(i, 1, 1:n, l)))
            c2 = DOT_PRODUCT(sint(i, 2, 1:n, l), MATMUL(omat(1:n, 1:n, l), sint(i, 2, 1:n, l)))
            WRITE (iw, "(T6,F14.6,T41,F10.6,T71,F10.6)") snl(l, i), c1, c2
         END DO
      END DO

      DEALLOCATE (omat, sfun, sint)

   END SUBROUTINE atom_completeness

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param rad ...
!> \param cnum ...
! **************************************************************************************************
   SUBROUTINE atom_basis_condnum(basis, rad, cnum)
      TYPE(atom_basis_type)                              :: basis
      REAL(KIND=dp), INTENT(IN)                          :: rad
      REAL(KIND=dp), INTENT(OUT)                         :: cnum

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_basis_condnum', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, imax, info, ix, iy, iz, ja, jb, &
                                                            ka, kb, l, la, lb, lwork, na, nb, &
                                                            nbas, nna, nnb
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ibptr
      REAL(KIND=dp)                                      :: r1, r2, reps, rmax
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weig, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: smat
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(7, 7)                     :: sab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zeta, zetb

      ! total number of basis functions
      nbas = 0
      DO l = 0, 3
         nbas = nbas+basis%nbas(l)*(2*l+1)
      END DO

      ALLOCATE (smat(nbas, nbas), ibptr(nbas, 0:3))
      smat = 0.0_dp
      ibptr = 0
      na = 0
      DO l = 0, 3
         DO ia = 1, basis%nbas(l)
            ibptr(ia, l) = na+1
            na = na+(2*l+1)
         END DO
      END DO

      reps = 1.e-14_dp
      IF (basis%basis_type == GTO_BASIS .OR. &
          basis%basis_type == CGTO_BASIS) THEN
         DO la = 0, 3
            na = basis%nprim(la)
            nna = 2*la+1
            IF (na == 0) CYCLE
            zeta => basis%am(:, la)
            DO lb = 0, 3
               nb = basis%nprim(lb)
               nnb = 2*lb+1
               IF (nb == 0) CYCLE
               zetb => basis%am(:, lb)
               DO ia = 1, na
                  DO ib = 1, nb
                     IF (rad < 0.1_dp) THEN
                        imax = 0
                     ELSE
                        r1 = exp_radius(la, zeta(ia), reps, 1.0_dp)
                        r2 = exp_radius(lb, zetb(ib), reps, 1.0_dp)
                        rmax = MAX(2._dp*r1, 2._dp*r2)
                        imax = INT(rmax/rad)+1
                     END IF
                     IF (imax > 1) THEN
                        CALL overlap_ab_sp(la, zeta(ia), lb, zetb(ib), rad, sab)
                        IF (basis%basis_type == GTO_BASIS) THEN
                           ja = ibptr(ia, la)
                           jb = ibptr(ib, lb)
                           smat(ja:ja+nna-1, jb:jb+nnb-1) = smat(ja:ja+nna-1, jb:jb+nnb-1)+sab(1:nna, 1:nnb)
                        ELSEIF (basis%basis_type == CGTO_BASIS) THEN
                           DO ka = 1, basis%nbas(la)
                              DO kb = 1, basis%nbas(lb)
                                 ja = ibptr(ka, la)
                                 jb = ibptr(kb, lb)
                                 smat(ja:ja+nna-1, jb:jb+nnb-1) = smat(ja:ja+nna-1, jb:jb+nnb-1)+ &
                                                                  sab(1:nna, 1:nnb)*basis%cm(ia, ka, la)*basis%cm(ib, kb, lb)
                              END DO
                           END DO
                        END IF
                     ELSE
                        DO ix = -imax, imax
                           rab(1) = rad*ix
                           DO iy = -imax, imax
                              rab(2) = rad*iy
                              DO iz = -imax, imax
                                 rab(3) = rad*iz
                                 CALL overlap_ab_s(la, zeta(ia), lb, zetb(ib), rab, sab)
                                 IF (basis%basis_type == GTO_BASIS) THEN
                                    ja = ibptr(ia, la)
                                    jb = ibptr(ib, lb)
                                    smat(ja:ja+nna-1, jb:jb+nnb-1) = smat(ja:ja+nna-1, jb:jb+nnb-1)+sab(1:nna, 1:nnb)
                                 ELSEIF (basis%basis_type == CGTO_BASIS) THEN
                                    DO ka = 1, basis%nbas(la)
                                       DO kb = 1, basis%nbas(lb)
                                          ja = ibptr(ka, la)
                                          jb = ibptr(kb, lb)
                                          smat(ja:ja+nna-1, jb:jb+nnb-1) = &
                                             smat(ja:ja+nna-1, jb:jb+nnb-1)+ &
                                             sab(1:nna, 1:nnb)*basis%cm(ia, ka, la)*basis%cm(ib, kb, lb)
                                       END DO
                                    END DO
                                 END IF
                              END DO
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
            END DO
         END DO
      ELSE
         CPABORT("Condition number not available for this basis type")
      END IF

      info = 0
      lwork = nbas*nbas
      ALLOCATE (weig(nbas), work(lwork))

      CALL lapack_ssyev("N", "U", nbas, smat, nbas, weig, work, lwork, info)
      CPASSERT(info == 0)
      IF (weig(1) < 0.0_dp) THEN
         cnum = 100._dp
      ELSE
         cnum = ABS(weig(nbas)/weig(1))
         cnum = LOG10(cnum)
      END IF

      DEALLOCATE (smat, ibptr, weig, work)

   END SUBROUTINE atom_basis_condnum
! **************************************************************************************************
!> \brief ...
!> \param int ...
!> \param omat ...
!> \param cm ...
! **************************************************************************************************
   SUBROUTINE contract2(int, omat, cm)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: int
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: omat, cm

      CHARACTER(len=*), PARAMETER :: routineN = 'contract2', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, m, n

      CALL timeset(routineN, handle)

      n = SIZE(int, 1)
      m = SIZE(omat, 1)

      INT(1:n, 1:n) = MATMUL(TRANSPOSE(cm(1:m, 1:n)), MATMUL(omat(1:m, 1:m), cm(1:m, 1:n)))

      CALL timestop(handle)

   END SUBROUTINE contract2
! **************************************************************************************************
!> \brief ...
!> \param int ...
!> \param omat ...
!> \param cm ...
! **************************************************************************************************
   SUBROUTINE contract2add(int, omat, cm)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: int
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: omat, cm

      CHARACTER(len=*), PARAMETER :: routineN = 'contract2add', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, m, n

      CALL timeset(routineN, handle)

      n = SIZE(int, 1)
      m = SIZE(omat, 1)

      INT(1:n, 1:n) = INT(1:n, 1:n)+MATMUL(TRANSPOSE(cm(1:m, 1:n)), MATMUL(omat(1:m, 1:m), cm(1:m, 1:n)))

      CALL timestop(handle)

   END SUBROUTINE contract2add
! **************************************************************************************************
!> \brief ...
!> \param eri ...
!> \param omat ...
!> \param cm1 ...
!> \param cm2 ...
! **************************************************************************************************
   SUBROUTINE contract4(eri, omat, cm1, cm2)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: eri
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: omat, cm1, cm2

      CHARACTER(len=*), PARAMETER :: routineN = 'contract4', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i1, i2, m1, m2, mm1, mm2, n1, &
                                                            n2, nn1, nn2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: amat, atran, bmat, btran, hint

      CALL timeset(routineN, handle)

      m1 = SIZE(cm1, 1)
      n1 = SIZE(cm1, 2)
      m2 = SIZE(cm2, 1)
      n2 = SIZE(cm2, 2)
      nn1 = SIZE(eri, 1)
      nn2 = SIZE(eri, 2)
      mm1 = SIZE(omat, 1)
      mm2 = SIZE(omat, 2)

      ALLOCATE (amat(m1, m1), atran(n1, n1), bmat(m2, m2), btran(n2, n2))
      ALLOCATE (hint(mm1, nn2))

      DO i1 = 1, mm1
         CALL iunpack(bmat(1:m2, 1:m2), omat(i1, 1:mm2), m2)
         CALL contract2(btran(1:n2, 1:n2), bmat(1:m2, 1:m2), cm2(1:m2, 1:n2))
         CALL ipack(btran(1:n2, 1:n2), hint(i1, 1:nn2), n2)
      END DO

      DO i2 = 1, nn2
         CALL iunpack(amat(1:m1, 1:m1), hint(1:mm1, i2), m1)
         CALL contract2(atran(1:n1, 1:n1), amat(1:m1, 1:m1), cm1(1:m1, 1:n1))
         CALL ipack(atran(1:n1, 1:n1), eri(1:nn1, i2), n1)
      END DO

      DEALLOCATE (amat, atran, bmat, btran)
      DEALLOCATE (hint)

      CALL timestop(handle)

   END SUBROUTINE contract4
! **************************************************************************************************
!> \brief ...
!> \param mat ...
!> \param vec ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ipack(mat, vec, n)
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: mat
      REAL(dp), DIMENSION(:), INTENT(INOUT)              :: vec
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: i, ij, j

      ij = 0
      DO i = 1, n
         DO j = i, n
            ij = ij+1
            vec(ij) = mat(i, j)
         END DO
      END DO

   END SUBROUTINE ipack

! **************************************************************************************************
!> \brief ...
!> \param mat ...
!> \param vec ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE iunpack(mat, vec, n)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: mat
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: vec
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: i, ij, j

      ij = 0
      DO i = 1, n
         DO j = i, n
            ij = ij+1
            mat(i, j) = vec(ij)
            mat(j, i) = vec(ij)
         END DO
      END DO

   END SUBROUTINE iunpack
! **************************************************************************************************

END MODULE atom_utils
