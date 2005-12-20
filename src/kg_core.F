!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_core [1.0] *
!!
!!   NAME
!!     kg_core
!!
!!   FUNCTION
!!     Calculation of the energies concerning the core charge distribution
!!
!!   SOURCE
!******************************************************************************

MODULE kg_core

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE erf_fn,                          ONLY: erfc
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: oorootpi,&
                                             pi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE kg_energy_types,                 ONLY: kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: kg_force_type
  USE qs_neighbor_lists,               ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             local_atoms_type,&
                                             subcell_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_core"

  PUBLIC :: calculate_ecore_rspace,&
            calculate_eself

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_ecore_rspace(kg_env,globenv,calculate_forces)

!   Purpose: Calculate the overlap energy of the core charge distribution.

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      pointer                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_ecore_overlap", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, group, handle, i, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom, jatom_local, jcell, &
      jkind, k, kkind, kcell, maxatom_local, natom, natom_local_a, natom_local_b, &
      nkind, output_unit, nexl, iexl
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zeff
    INTEGER, DIMENSION(3)                    :: cell_b, ncell, nsubcell, &
                                                periodic
    LOGICAL                                  :: cell000, include_ab, ionode
    REAL(KIND=dp) :: aab, dab, eab, ecore_overlap, f, fab, r_skin, rab2, &
      rab2_max, rab_max, rootaab, subcells, unit_of_length, zab
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: alpha, radius
    REAL(KIND=dp), DIMENSION(3)              :: deab, r, r_pbc, rab, rb, s, &
                                                sab_max, sb, sb_max, sb_min, &
                                                sb_pbc
    INTEGER, DIMENSION ( : ), POINTER        :: list
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I",'',handle)

    NULLIFY (atomic_kind)
    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (distribution_2d)
    NULLIFY (energy)
    NULLIFY (force)
    NULLIFY (particle_set)
    NULLIFY (subcell)
    NULLIFY (list)

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Set Verlet skin ***

    r_skin = 0.0_dp

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    distribution_2d=distribution_2d,&
                    particle_set=particle_set,&
                    energy=energy,&
                    force=force)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)


!   *** Allocate work storage ***

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom",nkind*int_size)

    ALLOCATE (alpha(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha",nkind*dp_size)
    alpha(:) = 0.0_dp

    ALLOCATE (radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"radius",nkind*dp_size)
    radius(:) = 0.0_dp

    ALLOCATE (zeff(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"zeff",nkind*dp_size)
    zeff(:) = 0.0_dp

!   *** Initialize the local data structures ***

    maxatom_local = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom(ikind)%list)
      NULLIFY (atom(ikind)%list_local_a_index)
      NULLIFY (atom(ikind)%list_local_b_index)
      NULLIFY (atom(ikind)%r_pbc)
      NULLIFY (atom(ikind)%s_pbc)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom(ikind)%list,&
                           alpha_core_charge=alpha(ikind),&
                           core_charge_radius=radius(ikind),&
                           qeff=zeff(ikind))

      natom = SIZE(atom(ikind)%list)
      natom_local_a = distribution_2d%n_local_rows(ikind)
      natom_local_b = distribution_2d%n_local_cols(ikind)

      maxatom_local = MAX(maxatom_local,natom_local_a,natom_local_b)

      IF (alpha(ikind) /= 0.0_dp) THEN

        IF (natom_local_a > 0) THEN

          ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"atom(ikind)%list_local_a_index",&
                             natom_local_a*int_size)
          END IF

!         *** Build index vector for mapping ***
          DO iatom_local=1,natom_local_a
            atom_a = distribution_2d%local_rows(ikind)%array(iatom_local)
            DO iatom=1,natom
              IF (atom_a == atom(ikind)%list(iatom)) THEN
                atom(ikind)%list_local_a_index(iatom_local) = iatom
                EXIT
              END IF
            END DO
          END DO

        END IF

        IF (natom_local_b > 0) THEN

          ALLOCATE (atom(ikind)%list_local_b_index(natom_local_b),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"atom(ikind)%list_local_b_index",&
                             natom_local_b*int_size)
          END IF

!         *** Build index vector for mapping ***

          DO iatom_local=1,natom_local_b
            atom_b = distribution_2d%local_cols(ikind)%array(iatom_local)
            DO iatom=1,natom
              IF (atom_b == atom(ikind)%list(iatom)) THEN
                atom(ikind)%list_local_b_index(iatom_local) = iatom
                EXIT
              END IF
            END DO
          END DO
        END IF

        ALLOCATE (atom(ikind)%r_pbc(3,natom),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"atom(ikind)%r_pbc",3*natom*dp_size)
        END IF

        ALLOCATE (atom(ikind)%s_pbc(3,natom),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"atom(ikind)%s_pbc",3*natom*dp_size)
        END IF

!       *** Calculate PBC coordinates ***

        DO iatom=1,natom
          atom_a = atom(ikind)%list(iatom)
          r_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
          atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:),cell)
        END DO

      END IF

    END DO

    IF ((ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > LOW)) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE CALCULATION OF THE CORE CHARGE OVERLAP",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

    ecore_overlap = 0.0_dp

    DO ikind=1,nkind

      IF (alpha(ikind) == 0.0_dp) CYCLE
      IF (zeff(ikind) == 0.0_dp) CYCLE

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)

      DO jkind=1,nkind

        IF (alpha(jkind) == 0.0_dp) CYCLE
        IF (zeff(jkind) == 0.0_dp) CYCLE

        IF (.NOT.ASSOCIATED(atom(jkind)%list_local_b_index)) CYCLE

        natom_local_b = SIZE(atom(jkind)%list_local_b_index)

        zab = REAL(zeff(ikind)*zeff(jkind),dp)
        aab = alpha(ikind)*alpha(jkind)/(alpha(ikind) + alpha(jkind))
        rootaab = SQRT(aab)
        fab = 2.0_dp*oorootpi*zab*rootaab

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius(ikind) + radius(jkind) + r_skin
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        nsubcell(:) = MIN(MAX(1,NINT(0.5_dp*subcells/sab_max(:))),20)

!       *** Print subcell information for the current atomic kind pair ***

        IF ((ionode.AND.globenv%print%cell_parameters).AND.&
            (globenv%print%level > LOW)) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = scaled_to_real(s(:),cell)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO jatom_local=1,natom_local_b

            jatom = atom(jkind)%list_local_b_index(jatom_local)
            atom_b = atom(jkind)%list(jatom)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom)

            loop1_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop1_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop1_kcell
              cell_b(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop1_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop1_jcell
                cell_b(2) = jcell

                loop1_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop1_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop1_icell
                  cell_b(1) = icell

                  rb(:) = scaled_to_real(sb(:),cell)

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    IF (atom_a == atom_b) THEN
                      IF (cell000) THEN
                        include_ab = .FALSE.
                      ELSE
                        include_ab = .TRUE.
                        f = 0.5_dp
                      END IF
                    ELSE IF (atom_a > atom_b) THEN
                      include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                      f = 1.0_dp
                    ELSE
                      include_ab = (MODULO(atom_a + atom_b,2) == 0)
                      f = 1.0_dp
                    END IF
! check for exclusion
                    IF(atom_a /= atom_b) THEN
                      loop_out: DO kkind = 1, nkind
                        list => particle_set(atom_b)%lpi%exclusion(kkind) % list_ei
                        IF (.NOT.ASSOCIATED(list)) CYCLE
                        nexl = SIZE ( list )
                        loop_in: DO iexl = 1, nexl
                          IF ((list(iexl)==atom_a).AND.cell000) THEN
                            include_ab = .FALSE.
                            EXIT loop_out
                          END IF
                        END DO loop_in
                      END DO loop_out
                    END IF

                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                      IF (rab2 < rab2_max) THEN
                        dab = SQRT(rab2)
                        eab = zab*erfc(rootaab*dab)/dab
                        ecore_overlap = ecore_overlap + f*eab
                        IF (calculate_forces) THEN
                          deab(:) = rab(:)*f*(eab + fab*EXP(-aab*rab2))/rab2
                          force(ikind)%f_rspace_core(:,iatom) =&
                            force(ikind)%f_rspace_core(:,iatom) + deab(:)
                          force(jkind)%f_rspace_core(:,jatom) =&
                            force(jkind)%f_rspace_core(:,jatom) - deab(:)
                        END IF
                      END IF
                    END IF
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

        ELSE

!         *** Case 2: A subcell grid is needed ***

          CALL allocate_subcell(subcell,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp)
            j = INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp)
            k = INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp)
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell(i,j,k)%natom =  subcell(i,j,k)%natom + 1
            subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
          END DO

          DO jatom_local=1,natom_local_b

            jatom = atom(jkind)%list_local_b_index(jatom_local)
            atom_b = atom(jkind)%list(jatom)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom)

            loop2_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop2_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop2_kcell
              cell_b(3) = kcell

              loop2_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop2_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop2_jcell
                cell_b(2) = jcell

                loop2_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop2_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop2_icell
                  cell_b(1) = icell

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  rb(:) = scaled_to_real(sb(:),cell)

                  loop_k: DO k=1,nsubcell(3)
                    loop_j: DO j=1,nsubcell(2)
                      loop_i: DO i=1,nsubcell(1)

                        IF (subcell(i,j,k)%natom == 0) CYCLE

                        IF (sb_max(3) < subcell(i,j,k)%s_min(3)) EXIT loop_k
                        IF (sb_min(3) >= subcell(i,j,k)%s_max(3)) CYCLE loop_k

                        IF (sb_max(2) < subcell(i,j,k)%s_min(2)) EXIT loop_j
                        IF (sb_min(2) >= subcell(i,j,k)%s_max(2)) CYCLE loop_j

                        IF (sb_max(1) < subcell(i,j,k)%s_min(1)) EXIT loop_i
                        IF (sb_min(1) >= subcell(i,j,k)%s_max(1)) CYCLE loop_i

                        DO iatom_subcell=1,subcell(i,j,k)%natom
                          iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                          iatom = atom(ikind)%list_local_a_index(iatom_local)
                          atom_a = atom(ikind)%list(iatom)
                          IF (atom_a == atom_b) THEN
                            IF (cell000) THEN
                              include_ab = .FALSE.
                            ELSE
                              include_ab = .TRUE.
                              f = 0.5_dp
                            END IF
                          ELSE IF (atom_a > atom_b) THEN
                            include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                            f = 1.0_dp
                          ELSE
                            include_ab = (MODULO(atom_a + atom_b,2) == 0)
                            f = 1.0_dp
                          END IF
! check for exclusion
                          IF(atom_a /= atom_b) THEN
                            loop_out2: DO kkind = 1, nkind
                              list => particle_set(atom_b)%lpi%exclusion(kkind) % list_ei
                              IF (.NOT.ASSOCIATED(list)) CYCLE
                              nexl = SIZE ( list )
                              loop_in2: DO iexl = 1, nexl
                                IF ((list(iexl)==atom_a).AND.cell000) THEN
                                  include_ab = .FALSE.
                                  EXIT loop_out2
                                END IF
                              END DO loop_in2
                            END DO loop_out2
                          END IF
                          IF (include_ab) THEN
                            rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                            IF (rab2 < rab2_max) THEN
                              dab = SQRT(rab2)
                              eab = zab*erfc(rootaab*dab)/dab
                              ecore_overlap = ecore_overlap + f*eab
                              IF (calculate_forces) THEN
                                deab(:) = rab(:)*f*(eab + fab*EXP(-aab*rab2))/rab2
                                force(ikind)%f_rspace_core(:,iatom) =&
                                  force(ikind)%f_rspace_core(:,iatom) + deab(:)
                                force(jkind)%f_rspace_core(:,jatom) =&
                                  force(jkind)%f_rspace_core(:,jatom) - deab(:)
                              END IF
                            END IF
                          END IF
                        END DO

                      END DO loop_i
                    END DO loop_j
                  END DO loop_k

                END DO loop2_icell
              END DO loop2_jcell
            END DO loop2_kcell

          END DO

          CALL deallocate_subcell(subcell)

        END IF

      END DO
    END DO

!   *** Release work storage ***

    DO ikind=1,nkind
      NULLIFY (atom(ikind)%list)
      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_a_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"atom(ikind)%list_local_a_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%list_local_b_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_b_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"atom(ikind)%list_local_b_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%r_pbc)) THEN
        DEALLOCATE (atom(ikind)%r_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"atom(ikind)%r_pbc")
      END IF
      IF (ASSOCIATED(atom(ikind)%s_pbc)) THEN
        DEALLOCATE (atom(ikind)%s_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"atom(ikind)%s_pbc")
      END IF
    END DO

    DEALLOCATE (atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom")

    DEALLOCATE (alpha,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha")

    DEALLOCATE (radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"radius")

    DEALLOCATE (zeff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"zeff")

    CALL mp_sum(ecore_overlap,group)

    energy%core_overlap = ecore_overlap

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ecore_rspace

! *****************************************************************************

  SUBROUTINE calculate_eself(kg_env,ecore_self)

!   Purpose: Calculate the self energy of the core charge distribution.

!   History: - Creation (14.10.2002)

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER       :: kg_env
    REAL(KIND=dp), INTENT(OUT)               :: ecore_self

    INTEGER                                  :: ikind, natom
    REAL(KIND=dp)                                :: alpha_core_charge, &
                                                ccore_charge, qeff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set)

    ecore_self = 0.0_dp

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           ccore_charge=ccore_charge,&
                           alpha_core_charge=alpha_core_charge)

      qeff=ccore_charge/SQRT(alpha_core_charge/pi)**3
      ecore_self = ecore_self - REAL(natom*qeff**2,dp)*SQRT(alpha_core_charge)

    END DO

    ecore_self = ecore_self/SQRT(twopi)

  END SUBROUTINE calculate_eself

!*****************************************************************************

END MODULE kg_core
