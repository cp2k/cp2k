!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations
!> \par History
!>      Adapted from qs_scf_post
!> \author JGH (19.07.2006)
! *****************************************************************************
MODULE kg_scf_post

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set,&
                                             dbcsr_obj_p_type,&
                                             sm_from_dbcsr
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_gpw_fm_mol_methods,           ONLY: mol_density_matrix
  USE kg_gpw_fm_mol_types,             ONLY: fm_mol_blocks_type,&
                                             get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type,&
                                             mol_mo_set_p_type
  USE kg_gpw_fm_mol_utils,             ONLY: copy_sparse2mol_block
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE orbital_pointers,                ONLY: indco
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_moments,                      ONLY: build_local_moment_matrix
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_scf_post'

  PUBLIC :: scf_post_calculation

CONTAINS

! *****************************************************************************
!> \brief collects post-scf calculations and prints info
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE scf_post_calculation(kg_env,error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: iter
    INTEGER :: first_atom, handle, i, ia, iao, imol_kind, imol_local, istat, &
      ix, iy, iz, jmol, l, last_atom, n_ao_kind, nat_mol, natom, nm, &
      nmol_kind, nmol_local, nmoments, num_mol, output_unit, unit_nr
    INTEGER, DIMENSION(2)                    :: nmo_kind
    LOGICAL                                  :: failure, first_time
    REAL(KIND=dp)                            :: charge, total_charge
    REAL(KIND=dp), DIMENSION(3)              :: r1, ra, rc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mmol, pmat, ref_point, refp, &
                                                work
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: moments_b
    TYPE(distribution_1d_type), POINTER      :: local_molecule
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(kg_fm_p_type), POINTER              :: kg_fm
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: moments
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(atom_list,atomic_kind_set,cell,fm_mol_blocks,kg_fm,&
            local_molecule,logger,molecule,molecule_kind, molecule_set,&
            moments,mos,para_env,particle_set,qs_env)
    NULLIFY(mmol, pmat, ref_point, refp, work, distribution_2d)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    IF (output_unit>0) THEN
      WRITE(output_unit,'(/,A,/)') " +++ Molecular Properties Calculation"
    END IF

    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL get_kg_env( kg_env=kg_env,particle_set=particle_set,&
                       atomic_kind_set=atomic_kind_set,&
                       cell=cell,local_molecules=local_molecule,&
                       molecule_set=molecule_set,&
                       para_env=para_env,&
                       sub_qs_env=qs_env,kg_fm_set=kg_fm,&
                       error=error)

      CALL get_qs_env(qs_env,input=input,&
           distribution_2d=distribution_2d,&
           error=error)

!!!!  Calculation of Moments
      IF (BTEST(cp_print_key_should_output(logger%iter_info,input,"DFT%KG%PRINT%MOMENTS",&
                first_time=first_time,error=error),cp_p_file)) THEN
        unit_nr=cp_print_key_unit_nr(logger,input,"DFT%KG%PRINT%MOMENTS",&
                extension=".moments",middle_name="",log_filename=.FALSE.,&
                error=error)

        CALL get_atomic_kind_set( atomic_kind_set, natom=natom )
        ALLOCATE ( ref_point(1:3,1:natom), STAT=istat )
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        num_mol = 0
        nmol_kind = SIZE(local_molecule%n_el,1)
        DO imol_kind = 1, nmol_kind
          nmol_local = local_molecule%n_el(imol_kind)
          num_mol = num_mol + nmol_local
          DO imol_local = 1, nmol_local
            jmol = local_molecule%list(imol_kind)%array(imol_local)
            molecule => molecule_set ( jmol )
            CALL get_molecule(molecule,molecule_kind=molecule_kind,&
                              first_atom=first_atom,last_atom=last_atom)
            CALL get_molecule_kind(molecule_kind,atom_list=atom_list,natom=natom)
            total_charge = 0._dp
            rc = 0._dp
            DO ia = first_atom, last_atom
              iao = ia - first_atom + 1
              CALL get_atomic_kind(atom_list(iao)%atomic_kind,core_charge=charge)
              total_charge = total_charge + charge
              IF ( ia == first_atom ) r1(:)=particle_set(ia)%r
              ra(:) = particle_set(ia)%r(:) - r1(:)
              rc(:) = rc(:) + charge * pbc(ra(:), cell)
            END DO
            IF ( total_charge > 0._dp ) THEN
              rc = rc/total_charge
            ELSE
              ra(:) = particle_set(first_atom)%r(:)
              rc = pbc(ra(:), cell)
            END IF
            rc(:) = rc(:) + r1(:)
            DO ia = first_atom, last_atom
              ref_point(1:3,ia) = rc(1:3)
            END DO
          END DO
        END DO
  
        CALL mp_sum(num_mol,para_env%group)
  
        ! dipole, quadrupole, octapole matrices
        CALL section_vals_val_get(input,"DFT%KG%PRINT%MOMENTS%NMOMENTS",&
             i_val=nmoments,error=error)


        NULLIFY(moments_b)
        CALL build_local_moment_matrix(qs_env,moments_b,nmoments,&
                                 ref_points=ref_point,error=error)
    NULLIFY(moments)!sm->dbcsr
    CALL allocate_matrix_set( moments, SIZE(moments_b), error )!sm->dbcsr
    DO i=1,SIZE(moments)!sm->dbcsr
       CALL sm_from_dbcsr(moments(i)%matrix, moments_b(i)%matrix, &
            distribution_2d,error)!sm->dbcsr
    ENDDO!sm->dbcsr
    CALL dbcsr_deallocate_matrix_set(moments_b,error)!sm->dbcsr

        nm = SIZE(moments,1)
        ALLOCATE ( mmol(1:nm,1:num_mol), STAT=istat )
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE ( refp(1:3,1:num_mol), STAT=istat )
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        mmol = 0._dp
        refp = 0._dp

        nmol_kind = SIZE(local_molecule%n_el,1)
        DO imol_kind = 1, nmol_kind
          NULLIFY(fm_mol_set,fm_mol_blocks)
          fm_mol_set => kg_fm%kg_fm_mol_set(imol_kind)
          CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                                 natom=nat_mol, &
                                 n_ao=n_ao_kind,n_mo=nmo_kind, &
                                 fm_mol_blocks=fm_mol_blocks)
          nmol_local = local_molecule%n_el(imol_kind)
          ALLOCATE ( pmat(1:n_ao_kind,1:n_ao_kind), STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO imol_local = 1, nmol_local
            jmol = local_molecule%list(imol_kind)%array(imol_local)
            mol_block => fm_mol_blocks(imol_local)
            CALL get_fm_mol_block(fm_mol_block=mol_block,&
                                  work=work, mos=mos)
            pmat = 0._dp
            DO i = 1, SIZE(mos,1)
              CALL mol_density_matrix(mos(i)%mo_set,pmat,error)
            END DO
  
            DO l = 1, SIZE(moments)
              CALL copy_sparse2mol_block(moments(l)%matrix, mol_block,&
                        work, nat_mol, n_ao_kind, n_ao_kind, error)
              charge = SUM ( pmat(1:n_ao_kind,1:n_ao_kind) * &
                             work(1:n_ao_kind,1:n_ao_kind) )
              mmol(l,jmol) = charge
            END DO
  
            molecule => molecule_set ( jmol )
            CALL get_molecule(molecule,molecule_kind=molecule_kind,&
                              first_atom=first_atom,last_atom=last_atom)
            CALL get_molecule_kind(molecule_kind,atom_list=atom_list,natom=natom)
            refp(:,jmol) = ref_point(:,first_atom)
            DO ia = first_atom, last_atom
              iao = ia - first_atom + 1
              CALL get_atomic_kind(atom_list(iao)%atomic_kind,core_charge=charge)
              ra(:) = particle_set(ia)%r(:) - refp(:,jmol)
              rc(:) = pbc(ra(:), cell)
              DO l = 1, SIZE(moments)
                ix = indco(1,l+1)
                iy = indco(2,l+1)
                iz = indco(3,l+1)
                IF ( SUM(indco(:,l+1)) == 0 ) THEN
                  mmol(l,jmol) = mmol(l,jmol) - charge
                ELSE
                  IF (ix > 0) mmol(l,jmol) = mmol(l,jmol) - charge * rc(1)**ix
                  IF (iy > 0) mmol(l,jmol) = mmol(l,jmol) - charge * rc(2)**iy
                  IF (iz > 0) mmol(l,jmol) = mmol(l,jmol) - charge * rc(3)**iz
                END IF
              END DO
            END DO

          END DO  ! imol_local

          DEALLOCATE ( pmat, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO  ! imol_kind

        CALL mp_sum(mmol,para_env%group)
        CALL mp_sum(refp,para_env%group)

        IF (unit_nr>0) THEN
           IF (first_time) THEN
              WRITE(unit=unit_nr,fmt="(a)")&
                 "iter_level     molecule      dipole(x,     y,       z,         tot) "
           END IF
           iter=cp_iter_string(logger%iter_info,error=error)
           DO jmol = 1, SIZE(mmol,2)
               WRITE(unit=unit_nr,&
                     fmt="(a,T16,i7,T28,F10.4,T41,F10.4,T52,F10.4,T63,F12.4)")&
                     iter(1:15), jmol, mmol(1:3,jmol),SQRT(SUM(mmol(1:3,jmol)**2))
           END DO 
        END IF

        IF (output_unit>0) THEN
          DO jmol = 1, SIZE(mmol,2)
             WRITE(output_unit,'(A,i7)') " Molecule Nr. ", jmol
             IF ( nmoments > 0 ) THEN
               WRITE(output_unit,'(A,T21,F20.10,T41,F20.10,T61,F20.10)') &
                 " Reference Point ",refp(1:3,jmol)
               WRITE(output_unit,'(A,T25,A,T37,A,T49,A,T67,A)') &
                 " Dipole Moment ","X","Y","Z","TOTAL"
               WRITE(output_unit,'(T20,F10.4,T32,F10.4,T44,F10.4,T60,F12.4,T75,A)') &
                 mmol(1:3,jmol),SQRT(SUM(mmol(1:3,jmol)**2)),"[a.u.]"
               WRITE(output_unit,'(T20,F10.4,T32,F10.4,T44,F10.4,T60,F12.4,T74,A)') &
                 mmol(1:3,jmol)*debye,SQRT(SUM(mmol(1:3,jmol)**2))*debye,&
                 "[Debye]"
             END IF
             IF ( nmoments > 1 ) THEN
               WRITE(output_unit,'(A,T30,A,T36,F10.4,T47,A,T53,F10.4,T65,A,T71,F10.4)') &
                " Quadrupole Moment [a.u.]",&
                "XX=",mmol(4,jmol),"   YY=",mmol(7,jmol),"   ZZ=",mmol(9,jmol)
               WRITE(output_unit,'(T30,A,T36,F10.4,T47,A,T53,F10.4,T65,A,T71,F10.4)') &
                "XY=",mmol(5,jmol),"   XZ=",mmol(6,jmol),"   YZ=",mmol(8,jmol)
             END IF
             IF ( nmoments > 2 ) THEN
               WRITE(output_unit,'(A)') " Octapole Moment [a.u.]"
             END IF
             IF ( nmoments > 3 ) THEN
               WRITE(output_unit,'(A)') " Higher Moment [a.u.]"
             END IF
          END DO
        END IF
  
        CALL deallocate_matrix_set (moments,error=error)

        DEALLOCATE ( ref_point, STAT=istat )
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE ( mmol, STAT=istat )
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE ( refp, STAT=istat )
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        CALL cp_print_key_finished_output(unit_nr, logger,input,&
             "DFT%KG%PRINT%MOMENTS", error=error)

      END IF  ! calculate Moments
    END IF  ! Failure

    IF (output_unit>0) THEN
      WRITE(output_unit,'(/,A)') " +++ End of Molecular Properties Calculation"
    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_post_calculation

END MODULE kg_scf_post

