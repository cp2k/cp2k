!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of routines handling the localization for molecular properties
! *****************************************************************************
MODULE qs_loc_molecules
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_minloc,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE moments_utils,                   ONLY: get_reference_point
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_loc_types,                    ONLY: qs_loc_env_new_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public ***
  PUBLIC :: wfc_to_molecule

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_molecules'

CONTAINS

! *****************************************************************************
!> \brief maps wfc's to molecules and also prints molecular dipoles
!> \param qs_env the qs_env in which the qs_env lives
!> \param qs_loc_env ...
!> \param loc_print_key ...
!> \param center ...
!> \param molecule_set ...
!> \param nspins ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE wfc_to_molecule(qs_env, qs_loc_env, loc_print_key, center, &
       molecule_set, nspins, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_loc_env_new_type), INTENT(IN)    :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfc_to_molecule', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: zeta
    COMPLEX(KIND=dp), DIMENSION(3)           :: ggamma, zphase
    INTEGER :: counter, first_atom, i, iatom, ikind, imol, imol_now, istate, &
      j, local_location, natom, natom_loc, natom_max, nkind, nmol, nstate, &
      output_unit, reference, stat
    INTEGER, POINTER                         :: wfc_to_atom_map( : )
    LOGICAL                                  :: do_berry, failure, ghost
    REAL(KIND=dp)                            :: charge_tot, dipole(3), dr(3), &
                                                mydist(2), ria( 3 ), theta, &
                                                zeff, zwfc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: charge_set
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dipole_set
    REAL(KIND=dp), DIMENSION(3)              :: ci, gvec, rcc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    REAL(KIND=dp), POINTER                   :: distance( : ), r( :, : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_type), POINTER             :: particle_set( : )

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    ! Molecular Dipoles availables only for nspin == 1
    IF( BTEST(cp_print_key_should_output(logger%iter_info,loc_print_key,&
               "MOLECULAR_DIPOLES",error=error),cp_p_file))THEN
       IF (nspins>1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Molecular Dipoles not implemented for SPIN multiplicity "//&
               "larger than 1!", error=error, error_level=cp_failure_level)
       END IF

       ! Setup reference point and some warning..
       reference = section_get_ival(loc_print_key,keyword_name="MOLECULAR_DIPOLES%REFERENCE",error=error)
       CALL section_vals_val_get(loc_print_key,"MOLECULAR_DIPOLES%REF_POINT",r_vals=ref_point,error=error)
       CALL section_vals_val_get(loc_print_key,"MOLECULAR_DIPOLES%PERIODIC",l_val=do_berry,error=error)
    END IF
    particle_set => qs_loc_env%particle_set
    para_env     => qs_loc_env%para_env
    local_molecules => qs_loc_env%local_molecules
    nstate = SIZE ( center, 2 )
    ALLOCATE(wfc_to_atom_map(nstate),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nkind = SIZE (local_molecules%n_el)
    natom = 0
    natom_max = 0
    DO ikind=1,nkind
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol
          i = local_molecules%list(ikind)%array(imol)
          molecule_kind => molecule_set(i)%molecule_kind
          CALL get_molecule_kind(molecule_kind=molecule_kind,natom=natom)
          natom_max = natom_max + natom
          IF (.NOT.ASSOCIATED(molecule_set(i)%lmi)) THEN
             ALLOCATE(molecule_set(i)%lmi,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             NULLIFY(molecule_set(i)%lmi%states)
          ENDIF
          molecule_set(i)%lmi%nstates=0
          IF(ASSOCIATED(molecule_set(i)%lmi%states)) THEN
             DEALLOCATE(molecule_set(i)%lmi%states,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END DO
    END DO
    natom_loc = natom_max
    natom = natom_max

    CALL mp_max(natom_max,para_env%group)

    ALLOCATE(r(3,natom_max),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(distance(natom_max),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !Zero all the stuff
    r(:,:) = 0.0_dp
    distance(:) = 1.E10_dp

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    counter = 0
    nkind = SIZE (local_molecules%n_el)
    DO ikind=1,nkind
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol
          i = local_molecules%list(ikind)%array(imol)
          molecule_kind => molecule_set(i)%molecule_kind
          first_atom = molecule_set(i)%first_atom
          CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)

          DO iatom=1,natom
             counter = counter + 1
             r ( :, counter ) = particle_set(first_atom+iatom-1)%r ( : )
          END DO
       END DO
    END DO

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    DO istate=1,nstate
       distance(:) = 1.E10_dp
       DO iatom=1,natom_loc
          dr(1) = r(1,iatom) - center(1,istate)
          dr(2) = r(2,iatom) - center(2,istate)
          dr(3) = r(3,iatom) - center(3,istate)
          ria = pbc ( dr, qs_loc_env%cell )
          distance(iatom) = SQRT ( DOT_PRODUCT ( ria, ria ) )
       END DO

       !combine distance() from all procs
       local_location = MAX(1,MINLOC(distance,DIM=1))

       mydist(1) = distance(local_location)
       mydist(2) = para_env%mepos

       CALL mp_minloc(mydist,para_env%group)

       IF(mydist(2)==para_env%mepos) THEN
          wfc_to_atom_map(istate)  =  local_location
       ELSE
          wfc_to_atom_map(istate)  =  0
       END IF
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    IF(natom_loc/=0) THEN
       DO istate = 1,nstate
          iatom = wfc_to_atom_map(istate)
          IF(iatom/=0) THEN
             counter = 0
             nkind = SIZE (local_molecules%n_el)
             DO ikind=1,nkind
                nmol = SIZE ( local_molecules%list(ikind)%array )
                DO imol=1,nmol
                   imol_now = local_molecules%list(ikind)%array(imol)
                   molecule_kind => molecule_set(imol_now)%molecule_kind
                   CALL get_molecule_kind(molecule_kind=molecule_kind,natom=natom)
                   counter = counter + natom
                   IF(counter>=iatom) EXIT
                END DO
                IF(counter>=iatom) EXIT
             END DO
             i = molecule_set(imol_now)%lmi%nstates
             i = i + 1
             molecule_set(imol_now)%lmi%nstates = i
             CALL reallocate(molecule_set(imol_now)%lmi%states,1,i)
             molecule_set(imol_now)%lmi%states(i) = istate
          END IF
       END DO
    END IF

    !---------------------------------------------------------------------------
    ! Figure out dipole of the molecule.
    !---------------------------------------------------------------------------
    IF( BTEST(cp_print_key_should_output(logger%iter_info,loc_print_key,&
               "MOLECULAR_DIPOLES",error=error),cp_p_file))THEN

    particle_set => qs_loc_env%particle_set
    para_env => qs_loc_env%para_env
    local_molecules => qs_loc_env%local_molecules
    nstate = SIZE ( center, 2 )
    ALLOCATE(dipole_set(3,SIZE(molecule_set)))
    ALLOCATE(charge_set(SIZE(molecule_set)))
    dipole_set=0.0_dp
    charge_set=0.0_dp
    cell      => qs_loc_env%cell
    zwfc      = 3.0_dp - REAL(nspins, KIND=dp)
    nkind = SIZE (local_molecules%n_el)
    DO ikind=1,nkind ! loop over different molecules
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol ! all the molecules of the kind
          imol_now = local_molecules%list(ikind)%array(imol) ! index in the global array
          IF (.NOT.ASSOCIATED(molecule_set(imol_now)%lmi%states)) CYCLE
          molecule_kind => molecule_set(imol_now)%molecule_kind
          first_atom = molecule_set(imol_now)%first_atom
          CALL get_molecule_kind(molecule_kind=molecule_kind,natom=natom)

          ! Get reference point for this molecule
          CALL get_reference_point(rcc,qs_env=qs_env,reference=reference,&
                                   ref_point=ref_point,ifirst=first_atom,&
                                   ilast=first_atom+natom-1,error=error)

          dipole=0.0_dp
          IF (do_berry) THEN
             rcc = pbc(rcc,cell)
             ! Find out the total charge of the molecule
             DO iatom=1,natom
                i=first_atom+iatom-1
                atomic_kind=>particle_set(i)%atomic_kind
                CALL get_atomic_kind(atomic_kind,ghost=ghost)
                IF (.NOT.ghost) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=zeff)
                   charge_set(imol_now)=charge_set(imol_now)+zeff
                END IF
             END DO
             ! Charges of the wfc involved
             DO istate=1,SIZE(molecule_set(imol_now)%lmi%states)
                charge_set(imol_now)=charge_set(imol_now)-zwfc
             ENDDO

             charge_tot = charge_set(imol_now)
             ria        = twopi * MATMUL(cell%h_inv, rcc)
             zphase     = CMPLX(COS(ria), SIN(ria), KIND=dp)**charge_tot
             ggamma     = CMPLX(1.0_dp,0.0_dp,KIND=dp)

             ! Nuclear charges
             DO iatom=1,natom
                i=first_atom+iatom-1
                atomic_kind=>particle_set(i)%atomic_kind
                CALL get_atomic_kind(atomic_kind,ghost=ghost)
                IF (.NOT.ghost) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=zeff)
                   ria = pbc(particle_set(i)%r,cell)
                   DO j = 1, 3
                      gvec  = twopi*cell%h_inv(j,:)
                      theta = SUM(ria(:)*gvec(:))
                      zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(zeff)
                      ggamma(j)  = ggamma(j) * zeta
                   END DO
                END IF
             END DO

             ! Charges of the wfc involved
             DO istate=1,SIZE(molecule_set(imol_now)%lmi%states)
                i = molecule_set(imol_now)%lmi%states(istate)
                ria = pbc(center(1:3,i),cell)
                DO j = 1, 3
                   gvec  = twopi*cell%h_inv(j,:)
                   theta = SUM(ria(:)*gvec(:))
                   zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(-zwfc)
                   ggamma(j)  = ggamma(j) * zeta
                END DO
             ENDDO

             ggamma = ggamma * zphase
             ci = AIMAG(LOG(ggamma))/twopi
             dipole = MATMUL(cell%hmat, ci)
          ELSE
             ! Nuclear charges
             DO iatom=1,natom
                i=first_atom+iatom-1
                atomic_kind=>particle_set(i)%atomic_kind
                CALL get_atomic_kind(atomic_kind,ghost=ghost)
                IF (.NOT.ghost) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=zeff)
                   ria = pbc(particle_set(i)%r,cell) - rcc
                   dipole=dipole + zeff*(ria-rcc)
                   charge_set(imol_now)=charge_set(imol_now)+zeff
                END IF
             END DO
             ! Charges of the wfc involved
             DO istate=1,SIZE(molecule_set(imol_now)%lmi%states)
                i = molecule_set(imol_now)%lmi%states(istate)
                ria = pbc(center(1:3,i),cell)
                dipole=dipole - zwfc * (ria-rcc)
                charge_set(imol_now)=charge_set(imol_now)-zwfc
             ENDDO
          END IF
          dipole_set(:,imol_now)=dipole ! a.u.
       ENDDO
    ENDDO
    CALL mp_sum(dipole_set,para_env%group)
    CALL mp_sum(charge_set,para_env%group)

    output_unit=cp_print_key_unit_nr(logger,loc_print_key,"MOLECULAR_DIPOLES",&
         extension=".MolDip",middle_name="MOLECULAR_DIPOLES",error=error)
    IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT='(A80)') "# molecule nr, charge, dipole vector, dipole (Debye) "
       dipole_set=dipole_set*debye ! Debye
       DO I=1,SIZE(dipole_set,2)
          WRITE(UNIT=output_unit,FMT='(I6,5F12.6)') I,charge_set(I),dipole_set(1:3,I), &
               SQRT(DOT_PRODUCT(dipole_set(1:3,I),dipole_set(1:3,I)))
       ENDDO
    ENDIF
    CALL cp_print_key_finished_output(output_unit,logger,loc_print_key,&
         "MOLECULAR_DIPOLES",error=error)
    DEALLOCATE(dipole_set,charge_set)
 END IF
    !---------------------------------------------------------------------------
    ! end of molecular dipole calculation
    !---------------------------------------------------------------------------

    DEALLOCATE(distance,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(wfc_to_atom_map,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE wfc_to_molecule
  !------------------------------------------------------------------------------

END MODULE qs_loc_molecules

