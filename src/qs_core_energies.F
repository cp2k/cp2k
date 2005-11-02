!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_core_energies [1.0] *
!!
!!   NAME
!!     qs_core_energies
!!
!!   FUNCTION
!!     Calculation of the energies concerning the core charge distribution
!!
!!   AUTHOR
!!     Matthias Krack (27.04.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_core_energies

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE erf_fn,                          ONLY: erfc
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: oorootpi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_lists,               ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             local_atoms_type,&
                                             subcell_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_core_energies"

  PUBLIC :: calculate_ecore,&
            calculate_ecore_overlap,&
            calculate_ecore_self

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_ecore(h,p,ecore,para_env)

!   Purpose: Calculate the core Hamiltonian energy which includes the kinetic
!            and the potential energy of the electrons. It is assumed, that
!            the core Hamiltonian matrix h and the density matrix p have the
!            same sparse matrix structure (same atomic blocks and block
!            ordering)

!   History: - Creation (03.05.2001,MK)
!            - simplified taking advantage of new non-redundant matrix
!              structure (27.06.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: h, p
    REAL(KIND=dp), INTENT(OUT)               :: ecore
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_ecore", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: handle, i, iblock_col, &
                                                iblock_row, j, jblock_col, &
                                                nblock_row
    REAL(KIND=dp)                            :: ecore_block
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block
    TYPE(real_block_node_type), POINTER      :: h_block_node, p_block_node

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I",'',handle)

    ecore = 0.0_dp

    CALL get_matrix_info(matrix=h,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

      h_block_node => first_block_node(matrix=h,block_row=iblock_row)
      p_block_node => first_block_node(matrix=p,block_row=iblock_row)

      DO WHILE (ASSOCIATED(h_block_node))

        CALL get_block_node(block_node=h_block_node,&
                            block_col=iblock_col,&
                            block=h_block)

        CALL get_block_node(block_node=p_block_node,&
                            block_col=jblock_col,&
                            block=p_block)

        IF (jblock_col /= iblock_col) THEN
          CALL stop_program(routine,&
                            "The matrix structures of the core Hamiltonian "//&
                            "and the density matrix does not match.")
        END IF

        ecore_block = 0.0_dp

        DO j=1,SIZE(h_block,2)
          DO i=1,SIZE(h_block,1)
            ecore_block = ecore_block + p_block(i,j)*h_block(i,j)
          END DO
        END DO

        IF (iblock_row == iblock_col) THEN
          ecore = ecore + ecore_block
        ELSE
          ecore = ecore + 2.0_dp*ecore_block
        END IF

        h_block_node => next_block_node(h_block_node)
        p_block_node => next_block_node(p_block_node)

      END DO

    END DO

    CALL mp_sum(ecore,para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ecore

! *****************************************************************************

  SUBROUTINE calculate_ecore_overlap(qs_env,globenv,calculate_forces,molecular, &
                                     E_overlap_core)

!   Purpose: Calculate the overlap energy of the core charge distribution.

!   History: - Creation (30.04.2001,MK)
!            - Force calculation added (03.06.2002,MK)
!            - Parallelized using a list of local atoms for rows and columns
!              (19.07.2003,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: molecular
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: E_overlap_core

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_ecore_overlap", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, group, handle, i, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom, jatom_local, jcell, &
      jkind, k, kcell, maxatom_local, natom, natom_local_a, natom_local_b, &
      nkind, output_unit
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zeff
    INTEGER, DIMENSION(3)                    :: cell_b, ncell, nsubcell, &
                                                periodic
    LOGICAL                                  :: cell000, include_ab, ionode, only_molecule
    REAL(KIND=dp) :: aab, dab, eab, ecore_overlap, f, fab, r_skin, rab2, &
      rab2_max, rab_max, rootaab, subcells, unit_of_length, zab
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: alpha, radius
    REAL(KIND=dp), DIMENSION(3)              :: deab, r, r_pbc, rab, rb, s, &
                                                sab_max, sb, sb_max, sb_min, &
                                                sb_pbc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell
    INTEGER, DIMENSION(:), POINTER           :: local_rows_array, &
                                                local_cols_array

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I",'',handle)

    NULLIFY (atomic_kind)
    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (distribution_2d,distribution_1d)
    NULLIFY (energy)
    NULLIFY (force)
    NULLIFY (particle_set)
    NULLIFY (subcell)

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr

    only_molecule = .FALSE.
    IF(PRESENT(molecular))  only_molecule = molecular
!   *** Set Verlet skin ***

    r_skin = 0.0_dp

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    local_particles=distribution_1d,&
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
    IF (istat /= 0) CALL stop_memory(routine,"zeff",nkind*int_size)
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
                           zeff=zeff(ikind))

      natom = SIZE(atom(ikind)%list)
      IF(only_molecule) THEN
        natom_local_a = distribution_1d%n_el(ikind)
        natom_local_b = distribution_1d%n_el(ikind)
        local_rows_array => distribution_1d%list(ikind)%array
        local_cols_array => distribution_1d%list(ikind)%array
      ELSE
        natom_local_a = distribution_2d%n_local_rows(ikind)
        natom_local_b = distribution_2d%n_local_cols(ikind)
        local_rows_array => distribution_2d%local_rows(ikind)%array
        local_cols_array => distribution_2d%local_cols(ikind)%array
      END IF

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
            atom_a = local_rows_array(iatom_local)
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
            atom_b = local_cols_array(iatom_local)
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

        zab = zeff(ikind)*zeff(jkind)
        aab = alpha(ikind)*alpha(jkind)/(alpha(ikind) + alpha(jkind))
        rootaab = SQRT(aab)
        fab = 2.0_dp*oorootpi*zab*rootaab

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius(ikind) + radius(jkind) + r_skin
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        nsubcell(:) = MAX(1,NINT(0.5_dp*subcells/sab_max(:)))

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
                    IF (cell000) THEN
                      IF (atom_a > atom_b) THEN
                        include_ab = (MODULO(atom_a + atom_b,2) == 0)
                        f = 1.0_dp
                      ELSE
                        include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                        f = 1.0_dp
                      END IF
                    ELSE
                      IF (atom_a == atom_b) THEN
                        include_ab = .TRUE.
                        f = 0.5_dp
                      ELSE IF (atom_a > atom_b) THEN
                        include_ab = (MODULO(atom_a + atom_b,2) == 0)
                        f = 1.0_dp
                      ELSE
                        include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                        f = 1.0_dp
                      END IF
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
                          force(ikind)%core_overlap(:,iatom) =&
                            force(ikind)%core_overlap(:,iatom) + deab(:)
                          force(jkind)%core_overlap(:,jatom) =&
                            force(jkind)%core_overlap(:,jatom) - deab(:)
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
                          IF (cell000) THEN
                            IF (atom_a > atom_b) THEN
                              include_ab = (MODULO(atom_a + atom_b,2) == 0)
                              f = 1.0_dp
                            ELSE
                              include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                              f = 1.0_dp
                            END IF
                          ELSE
                            IF (atom_a == atom_b) THEN
                              include_ab = .TRUE.
                              f = 0.5_dp
                            ELSE IF (atom_a > atom_b) THEN
                              include_ab = (MODULO(atom_a + atom_b,2) == 0)
                              f = 1.0_dp
                            ELSE
                              include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                              f = 1.0_dp
                            END IF
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
                                force(ikind)%core_overlap(:,iatom) =&
                                  force(ikind)%core_overlap(:,iatom) + deab(:)
                                force(jkind)%core_overlap(:,jatom) =&
                                  force(jkind)%core_overlap(:,jatom) - deab(:)
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
    IF (PRESENT(E_overlap_core)) THEN
      E_overlap_core = energy%core_overlap
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ecore_overlap

! *****************************************************************************

  SUBROUTINE calculate_ecore_self(qs_env, E_self_core)

!   Purpose: Calculate the self energy of the core charge distribution.

!   History: - Creation (27.04.2001,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: E_self_core

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_ecore_self", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: handle, ikind, natom
    REAL(KIND=dp)                            :: alpha_core_charge, ecore_self, zeff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(qs_energy_type), POINTER            :: energy

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    energy=energy)

    ecore_self = 0.0_dp

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           zeff=zeff,&
                           alpha_core_charge=alpha_core_charge)

      ecore_self = ecore_self - REAL(natom,dp)*zeff**2*SQRT(alpha_core_charge)

    END DO

    energy%core_self = ecore_self/SQRT(twopi)
    IF (PRESENT(E_self_core)) THEN
      E_self_core = energy%core_self
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ecore_self

! *****************************************************************************

END MODULE qs_core_energies
