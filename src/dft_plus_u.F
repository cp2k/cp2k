!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!> \brief   Add the DFT+U contribution to the Hamiltonian matrix
!> \details The implemented methods refers to:\n
!>          S. L. Dudarev, D. Nguyen Manh, and A. P. Sutton,
!>          Philos. Mag. B \b 75, 613 (1997)\n
!>          S. L. Dudarev et al.,
!>          Phys. Rev. B \b 57, 1505 (1998)
!> \author  Matthias Krack (MK)
!> \date    14.01.2008
!> \version 1.0
MODULE dft_plus_u

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE bibliography,                    ONLY: Dudarev1997,&
                                             Dudarev1998,&
                                             cite_reference
  USE cp_control_types,                ONLY: dft_control_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dft_plus_u'

  PUBLIC :: plus_u

CONTAINS

  !> \brief         Add the DFT+U contribution to the Hamiltonian matrix.\n
  !>                Wrapper routine for all "+U" methods
  !> \param[in]     qs_env Quickstep environment
  !> \param[in,out] matrix_h Hamiltonian matrices for each spin
  !> \param[in,out] matrix_w Energy weighted density matrices for each spin
  !> \param[in]     error environment
  !> \author        Matthias Krack (MK)
  !> \date          14.01.2008
  !> \version       1.0
  SUBROUTINE plus_u(qs_env,matrix_h,matrix_w,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_w
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'plus_u', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    CALL cite_reference(Dudarev1997)
    CALL cite_reference(Dudarev1998)

    ! There is yet no need for a method keyword, since this is the first DFT+U
    ! method implemented

    ! SELECT CASE (dft_plus_u_method_id)
    ! CASE (1)
        CALL plus_u_mulliken(qs_env,matrix_h,matrix_w,error)
    ! END SELECT

  END SUBROUTINE plus_u

  !> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
  !>                using a method based on Mulliken charges
  !>                \f[q_\mu = \sum_\nu P_{\mu\nu} S_{\nu\mu}\f]
  !>                where \b P and \b S are the density and the
  !>                overlap matrix, respectively.
  !> \details
  !>  \f{eqnarray*}
  !>   E^{\rm DFT+U} & = & E^{\rm DFT} +  E^{\rm U}\\\
  !>                 & = & E^{\rm DFT} +  \frac{1}{2}(U - J)\sum_\mu (q_\mu - q_\mu^2)\\[1ex]
  !>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
  !>                          & = & \frac{\partial E^{\rm DFT}}
  !>                                     {\partial P_{\mu\nu}} +
  !>                                \frac{\partial E^{\rm U}}
  !>                                     {\partial P_{\mu\nu}}\\\
  !>                          & = & H_{\mu\nu} +
  !>                                \frac{\partial E^{\rm U}}{\partial q_\mu}
  !>                                \frac{\partial q_\mu}{\partial P_{\mu\nu}}\\\
  !>                          & = & H_{\mu\nu} +
  !>                                \frac{1}{2}(U - J)(1 - q_\mu - q_\nu) S_{\mu\nu}\\\
  !>  \f}
  !> \param[in]     qs_env Quickstep environment
  !> \param[in,out] matrix_h Hamiltonian matrices for each spin
  !> \param[in,out] matrix_w Energy weighted density matrices for each spin
  !> \param[in]     error environment
  !> \author        Matthias Krack (MK)
  !> \date          11.01.2008
  !> \version       1.0
  !> \note          The use of any full matrices was avoided. Thus no ScaLAPACK
  !>                calls are performed
  SUBROUTINE plus_u_mulliken(qs_env,matrix_h,matrix_w,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_w
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'plus_u_mulliken', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom, handle, iatom, ikind, iset, isgf, ishell, ispin, jatom, &
      jsgf, katom, l_of_dft_plus_u, natom, natom_kind, nkind, nset, nsgf, &
      nspin, sgf, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf_atom
    INTEGER, DIMENSION(:), POINTER           :: atom_list, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgf, l, last_sgf
    LOGICAL                                  :: dft_plus_u_atom, failure, &
                                                just_energy, orthonormal_basis
    REAL(KIND=dp)                            :: fspin, q, u_minus_j, v
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dEdq, trps
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block, s_block, &
                                                w_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_block_node_type), POINTER      :: h_block_node, p_block_node, &
                                                s_block_node, w_block_node
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(real_matrix_type), POINTER          :: sm_h, sm_p, sm_s, sm_w

! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY (atom_list)
    NULLIFY (atomic_kind)
    NULLIFY (atomic_kind_set)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (first_sgf)
    NULLIFY (h_block)
    NULLIFY (h_block_node)
    NULLIFY (matrix_p)
    NULLIFY (matrix_s)
    NULLIFY (l)
    NULLIFY (last_sgf)
    NULLIFY (nshell)
    NULLIFY (orb_basis_set)
    NULLIFY (p_block)
    NULLIFY (p_block_node)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (s_block)
    NULLIFY (s_block_node)
    NULLIFY (sm_h)
    NULLIFY (sm_p)
    NULLIFY (sm_s)
    NULLIFY (w_block)
    NULLIFY (w_block_node)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    energy=energy,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    rho=rho,&
                    error=error)

    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(energy),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)

    sm_s => matrix_s(1)%matrix ! Overlap matrix in sparse format
    matrix_p => rho%rho_ao     ! Density matrices in sparse format

    energy%dft_plus_u = 0.0_dp

    nspin = dft_control%nspins

    IF (nspin == 2) THEN
      fspin = 1.0_dp
    ELSE
      fspin = 0.5_dp
    END IF

    ! Later we could save here some time, if the method in use has this property

    orthonormal_basis = .FALSE.

    ! Get the total number of atoms, contracted spherical Gaussian basis
    ! functions, and atomic kinds

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (first_sgf_atom(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    first_sgf_atom(:) = 0

    ALLOCATE (trps(nsgf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    trps(:) = 0.0_dp

    IF (PRESENT(matrix_h).OR.PRESENT(matrix_w)) THEN
      ALLOCATE (dEdq(nsgf),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      just_energy  = .FALSE.
    ELSE
      just_energy = .TRUE.
    END IF

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_atom)

    ! Loop over all spins

    DO ispin=1,nspin

      IF (PRESENT(matrix_h)) THEN
        ! Hamiltonian matrix for spin ispin in sparse format
        sm_h => matrix_h(ispin)%matrix
        dEdq(:) = 0.0_dp
      ELSE
        NULLIFY (sm_h)
      END IF

      IF (PRESENT(matrix_w)) THEN
        ! Energy weighted density matrix for spin ispin in sparse format
        sm_w => matrix_w(ispin)%matrix
        dEdq(:) = 0.0_dp
      ELSE
        NULLIFY (sm_w)
      END IF

      sm_p => matrix_p(ispin)%matrix ! Density matrix for spin ispin in sparse format

      ! Calculate Trace(P*S) assuming symmetric matrices

      trps(:) = 0.0_dp

      DO iatom=1,natom

        IF (orthonormal_basis) THEN

          ! Just retrieve the diagonal elements of the orthonormal density matrix P

          CALL get_block_node(matrix=sm_p,&
                              block_row=iatom,&
                              block_col=iatom,&
                              block=p_block)

          CPPostcondition(ASSOCIATED(p_block),cp_failure_level,routineP,error,failure)

          sgf = first_sgf_atom(iatom)

          DO isgf=1,SIZE(p_block,1)
            trps(sgf) = trps(sgf) + p_block(isgf,isgf)
            sgf = sgf + 1
          END DO

        ELSE

          p_block_node => first_block_node(matrix=sm_p,block_row=iatom)
          s_block_node => first_block_node(matrix=sm_s,block_row=iatom)

          DO WHILE (ASSOCIATED(p_block_node).AND.ASSOCIATED(s_block_node))

            CALL get_block_node(block_node=p_block_node,&
                                block_col=jatom,&
                                block=p_block)

            ! Request katom just to check for a consistent sparse matrix pattern

            CALL get_block_node(block_node=s_block_node,&
                                block_col=katom,&
                                block=s_block)

            CPPostcondition(jatom==katom,cp_failure_level,routineP,error,failure)
            CPPostcondition(jatom>=iatom,cp_failure_level,routineP,error,failure)

            sgf = first_sgf_atom(jatom)
            DO jsgf=1,SIZE(p_block,2)
              DO isgf=1,SIZE(p_block,1)
                trps(sgf) = trps(sgf) + p_block(isgf,jsgf)*s_block(isgf,jsgf)
              END DO
              sgf = sgf + 1
            END DO

            IF (iatom /= jatom) THEN
              sgf = first_sgf_atom(iatom)
              DO isgf=1,SIZE(p_block,1)
                DO jsgf=1,SIZE(p_block,2)
                  trps(sgf) = trps(sgf) + p_block(isgf,jsgf)*s_block(isgf,jsgf)
                END DO
                sgf = sgf + 1
              END DO
            END IF

            p_block_node => next_block_node(p_block_node)
            s_block_node => next_block_node(s_block_node)

          END DO ! next block node -> "jatom"

        END IF ! orthonormal basis set

      END DO ! next atom "iatom"

      CALL mp_sum(trps,qs_env%para_env%group)

      ! q <- Trace(PS)

      ! E[DFT+U] = E[DFT] + E[U]
      !          = E[DFT] + (U - J)*(q - q**2))/2

      ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
      !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
      !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

      ! Loop over all atomic kinds

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        ! Load the required atomic kind data

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             atom_list=atom_list,&
                             dft_plus_u_atom=dft_plus_u_atom,&
                             l_of_dft_plus_u=l_of_dft_plus_u,&
                             natom=natom_kind,&
                             orb_basis_set=orb_basis_set,&
                             u_minus_j=u_minus_j)

        ! Check, if this atom needs a DFT+U correction

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
        IF (.NOT.dft_plus_u_atom) CYCLE
        IF (l_of_dft_plus_u < 0) CYCLE
        IF (u_minus_j == 0.0_dp) CYCLE

        ! Load the required Gaussian basis set data

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgf,&
                               l=l,&
                               last_sgf=last_sgf,&
                               nset=nset,&
                               nshell=nshell)

        ! Loop over all atoms of the current atomic kind

        DO iatom=1,natom_kind

          atom = atom_list(iatom)

          ! Get diagonal block

          CALL get_block_node(matrix=sm_s,&
                              block_row=atom,&
                              block_col=atom,&
                              block=s_block)

          ! Calculate E(U) and dE(U)/dq

          IF (ASSOCIATED(s_block)) THEN

            sgf = first_sgf_atom(atom)

            DO iset=1,nset
              DO ishell=1,nshell(iset)
                IF (l(ishell,iset) == l_of_dft_plus_u) THEN
                  DO isgf=first_sgf(ishell,iset),last_sgf(ishell,iset)
                    q = fspin*trps(sgf)
                    energy%dft_plus_u = energy%dft_plus_u +&
                                        0.5_dp*u_minus_j*(q - q**2)/fspin
                    IF (.NOT.just_energy) THEN
                      dEdq(sgf) = dEdq(sgf) + u_minus_j*(0.5_dp - q)
                    END IF
                    sgf = sgf + 1
                  END DO ! next contracted spherical Gaussian function "isgf"
                ELSE
                  sgf = sgf + last_sgf(ishell,iset) - first_sgf(ishell,iset) + 1
                END IF ! angular momentum requested for DFT+U correction
              END DO ! next shell "ishell"
            END DO ! next shell set "iset"

          END IF ! this process is the owner of the sparse matrix block?

        END DO ! next atom "iatom" of atomic kind "ikind"

      END DO ! next atomic kind "ikind"

      IF (.NOT.just_energy) THEN
        CALL mp_sum(dEdq,qs_env%para_env%group)
      END IF

      ! Add V(i,j)[U] to V(i,j)[DFT

      IF (ASSOCIATED(sm_h)) THEN

        DO iatom=1,natom

          IF (orthonormal_basis) THEN

            ! Only the diagonal elements of the Hamiltonian matrix have to be updated

            CALL get_block_node(matrix=sm_h,&
                                block_row=iatom,&
                                block_col=iatom,&
                                block=h_block)

            CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)

            sgf = first_sgf_atom(iatom)

            DO isgf=1,SIZE(h_block,1)
              h_block(isgf,isgf) = h_block(isgf,isgf) + dEdq(sgf)
              sgf = sgf + 1
            END DO

          ELSE

            h_block_node => first_block_node(matrix=sm_h,block_row=iatom)
            s_block_node => first_block_node(matrix=sm_s,block_row=iatom)

            DO WHILE (ASSOCIATED(h_block_node).AND.ASSOCIATED(s_block_node))

              CALL get_block_node(block_node=h_block_node,&
                                  block_col=jatom,&
                                  block=h_block)

              ! Request katom just to check for consistent sparse matrix pattern

              CALL get_block_node(block_node=s_block_node,&
                                  block_col=katom,&
                                  block=s_block)

              CPPostcondition(jatom==katom,cp_failure_level,routineP,error,failure)
              CPPostcondition(jatom>=iatom,cp_failure_level,routineP,error,failure)

              sgf = first_sgf_atom(iatom)

              DO isgf=1,SIZE(h_block,1)
                IF (dEdq(sgf) /= 0.0_dp) THEN
                  v = 0.5_dp*dEdq(sgf)
                  DO jsgf=1,SIZE(h_block,2)
                    h_block(isgf,jsgf) = h_block(isgf,jsgf) + v*s_block(isgf,jsgf)
                  END DO
                END IF
                sgf = sgf + 1
              END DO

              sgf = first_sgf_atom(jatom)

              DO jsgf=1,SIZE(h_block,2)
                IF (dEdq(sgf) /= 0.0_dp) THEN
                  v = 0.5_dp*dEdq(sgf)
                  DO isgf=1,SIZE(h_block,1)
                    h_block(isgf,jsgf) = h_block(isgf,jsgf) + v*s_block(isgf,jsgf)
                  END DO
                END IF
                sgf = sgf + 1
              END DO

              h_block_node => next_block_node(h_block_node)
              s_block_node => next_block_node(s_block_node)

            END DO ! next block node -> "jatom"

          END IF ! orthonormal basis set

        END DO ! Next atom "iatom"

      END IF ! An update of the Hamiltonian matrix is requested

      ! Calculate the contribution (non-Pulay part) to the derivatives
      ! w.r.t. the nuclear positions, which requires an update of the
      ! energy weighted density W.

      IF (ASSOCIATED(sm_w).AND.(.NOT.orthonormal_basis)) THEN

        DO iatom=1,natom

          p_block_node => first_block_node(matrix=sm_p,block_row=iatom)
          w_block_node => first_block_node(matrix=sm_w,block_row=iatom)

          DO WHILE (ASSOCIATED(p_block_node).AND.ASSOCIATED(w_block_node))

            CALL get_block_node(block_node=p_block_node,&
                                block_col=jatom,&
                                block=p_block)

            ! Request katom just to check for consistent sparse matrix patterns

            CALL get_block_node(block_node=w_block_node,&
                                block_col=katom,&
                                block=w_block)

            CPPostcondition(jatom==katom,cp_failure_level,routineP,error,failure)
            CPPostcondition(jatom>=iatom,cp_failure_level,routineP,error,failure)

            sgf = first_sgf_atom(iatom)

            DO isgf=1,SIZE(w_block,1)
              IF (dEdq(sgf) /= 0.0_dp) THEN
                v = -0.5_dp*dEdq(sgf)
                DO jsgf=1,SIZE(w_block,2)
                  w_block(isgf,jsgf) = w_block(isgf,jsgf) + v*p_block(isgf,jsgf)
                END DO
              END IF
              sgf = sgf + 1
            END DO

            sgf = first_sgf_atom(jatom)

            DO jsgf=1,SIZE(w_block,2)
              IF (dEdq(sgf) /= 0.0_dp) THEN
                v = -0.5_dp*dEdq(sgf)
                DO isgf=1,SIZE(w_block,1)
                  w_block(isgf,jsgf) = w_block(isgf,jsgf) + v*p_block(isgf,jsgf)
                END DO
              END IF
              sgf = sgf + 1
            END DO

            p_block_node => next_block_node(p_block_node)
            w_block_node => next_block_node(w_block_node)

          END DO ! next block node -> "jatom"

        END DO ! next atom "iatom"

      END IF ! W matrix update requested

    END DO ! next spin "ispin"

    ! Collect the energy contributions from all processes

    CALL mp_sum(energy%dft_plus_u,qs_env%para_env%group)

    CALL cp_assert(energy%dft_plus_u >= 0.0_dp,cp_warning_level,&
                   cp_assertion_failed,routineP,&
                   "DFT+U energy contibution is negative possibly due "//&
                   "to unphysical Mulliken charges. Check your input, "//&
                   "if this warning persists!",&
                   error,failure)

    DEALLOCATE (first_sgf_atom,STAT=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    DEALLOCATE (trps,STAT=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    IF (ALLOCATED(dEdq)) THEN
      DEALLOCATE (dEdq,STAT=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

   CALL timestop(0.0_dp,handle)

  END SUBROUTINE plus_u_mulliken

END MODULE dft_plus_u
