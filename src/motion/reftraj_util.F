!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize the analysis of trajectories to be done
!>      by activating the REFTRAJ ensemble
!> \par History
!>      Created 10-07 [MI]
!> \author MI
! *****************************************************************************
MODULE reftraj_util

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds
  USE reftraj_types,                   ONLY: reftraj_msd_type,&
                                             reftraj_type
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE util,                            ONLY: get_limit
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

      PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reftraj_util'

      PUBLIC ::   initialize_reftraj, compute_msd_reftraj, write_output_reftraj

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param reftraj ...
!> \param reftraj_section ...
!> \param md_env ...
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
  SUBROUTINE initialize_reftraj(reftraj,reftraj_section,md_env)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_reftraj', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natom, nline_to_skip, nskip
    LOGICAL                                  :: failure, my_end
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: msd_section
    TYPE(simpar_type), POINTER               :: simpar

    failure = .FALSE.

    NULLIFY (force_env, msd_section, particles, simpar, subsys)
    CALL get_md_env(md_env=md_env, force_env=force_env, para_env=para_env,&
         simpar=simpar)
    CALL force_env_get(force_env=force_env,  subsys=subsys)
    CALL cp_subsys_get(subsys=subsys, particles=particles)
    natom = particles%n_els

    my_end = .FALSE.
    nline_to_skip = 0

    nskip = reftraj%info%first_snapshot-1
    CPPostcondition(nskip>=0,cp_failure_level,routineP,failure)

    IF(nskip > 0 ) THEN
      nline_to_skip = (natom+2)*nskip
      CALL parser_get_next_line(reftraj%info%traj_parser,nline_to_skip,at_end=my_end)
    END IF

    reftraj%isnap = nskip
    CALL cp_assert(.NOT.my_end,cp_fatal_level,cp_assertion_failed,routineP,&
         "Reached the end of the trajectory file for REFTRAJ. Number of steps skipped "//&
         "equal to the number of steps present in the file. "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    ! Cell File
    IF(reftraj%info%variable_volume) THEN
       IF(nskip > 0 ) THEN
          CALL parser_get_next_line(reftraj%info%cell_parser,nskip,at_end=my_end)
       END IF
       CALL cp_assert(.NOT.my_end,cp_fatal_level,cp_assertion_failed,routineP,&
            "Reached the end of the cell file for REFTRAJ. Number of steps skipped "//&
            "equal to the number of steps present in the file. "//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
    END IF


    reftraj%natom = natom
    IF(reftraj%info%last_snapshot>0) THEN
       simpar%nsteps = (reftraj%info%last_snapshot - reftraj%info%first_snapshot + 1)
    END IF

    IF(reftraj%info%msd) THEN
      msd_section =>  section_vals_get_subs_vals(reftraj_section,"MSD")
      ! set up and printout
      CALL initialize_msd_reftraj(reftraj%msd,msd_section,reftraj,md_env)
    END IF

  END SUBROUTINE initialize_reftraj

! *****************************************************************************
!> \brief ...
!> \param msd ...
!> \param msd_section ...
!> \param reftraj ...
!> \param md_env ...
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
  SUBROUTINE initialize_msd_reftraj(msd,msd_section,reftraj,md_env)
    TYPE(reftraj_msd_type), POINTER          :: msd
    TYPE(section_vals_type), POINTER         :: msd_section
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_msd_reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: AA
    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: name, title
    INTEGER :: first_atom, iatom, ikind, imol, istat, last_atom, natom_read, &
      nkind, nmol, nmolecule, nmolkind, npart
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: com(3), mass, mass_mol, tol, &
                                                x, y, z
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.

    NULLIFY (molecule, molecules, molecule_kind, molecule_kind_set,&
             molecule_kinds, molecule_set, subsys, force_env, particles, particle_set)
    CPPrecondition(.NOT. ASSOCIATED(msd),cp_failure_level,routineP,failure)

    ALLOCATE(msd, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

    NULLIFY(msd%ref0_pos)
    NULLIFY(msd%ref0_com_molecule)
    NULLIFY(msd%val_msd_kind)
    NULLIFY(msd%val_msd_molecule)
    NULLIFY(msd%disp_atom_index)
    NULLIFY(msd%disp_atom_dr)

    CALL get_md_env(md_env=md_env, force_env=force_env, para_env=para_env)
    CALL force_env_get(force_env=force_env,subsys=subsys)
    CALL cp_subsys_get(subsys=subsys, particles=particles)
    particle_set => particles%els
    npart = SIZE(particle_set,1)

    msd%ref0_unit = -1
    CALL section_vals_val_get(msd_section,"REF0_FILENAME",c_val=filename)
    CALL open_file(TRIM(filename),unit_number=msd%ref0_unit)

    ALLOCATE(msd%ref0_pos(3,reftraj%natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    msd%ref0_pos = 0.0_dp

    IF (para_env%mepos==para_env%source) THEN
      REWIND(msd%ref0_unit)
      READ(msd%ref0_unit,*,ERR=999,END=998) natom_read
      CPPostcondition(natom_read==reftraj%natom,cp_failure_level,routineP,failure)
      READ(msd%ref0_unit,'(A)',ERR=999,END=998) title
      msd%total_mass = 0.0_dp
      msd%ref0_com = 0.0_dp
      DO iatom = 1,natom_read
          READ(msd%ref0_unit,*,ERR=999,END=998)  AA, x, y, z
          name=TRIM(particle_set(iatom)%atomic_kind%element_symbol)
          CPPostcondition((TRIM(AA)==name),cp_failure_level,routineP,failure)

          x = cp_unit_to_cp2k(x,"angstrom")
          y = cp_unit_to_cp2k(y,"angstrom")
          z = cp_unit_to_cp2k(z,"angstrom")
          msd%ref0_pos(1,iatom) = x
          msd%ref0_pos(2,iatom) = y
          msd%ref0_pos(3,iatom) = z
          mass = particle_set(iatom)%atomic_kind%mass
          msd%ref0_com(1) = msd%ref0_com(1) + x * mass
          msd%ref0_com(2) = msd%ref0_com(2) + y * mass
          msd%ref0_com(3) = msd%ref0_com(3) + z * mass
          msd%total_mass = msd%total_mass + mass
      END DO
      msd%ref0_com = msd%ref0_com / msd%total_mass
    END IF
    CALL close_file(unit_number=msd%ref0_unit)

    CALL mp_bcast(msd%total_mass,para_env%source,para_env%group)
    CALL mp_bcast(msd%ref0_pos,para_env%source,para_env%group)
    CALL mp_bcast(msd%ref0_com,para_env%source,para_env%group)

    CALL section_vals_val_get(msd_section,"MSD_PER_KIND",l_val=msd%msd_kind)
    CALL section_vals_val_get(msd_section,"MSD_PER_MOLKIND",l_val=msd%msd_molecule)
    CALL section_vals_val_get(msd_section,"MSD_PER_REGION",l_val=msd%msd_region)

    CALL section_vals_val_get(msd_section,"DISPLACED_ATOM",l_val=msd%disp_atom)
    IF(msd%disp_atom) THEN
       ALLOCATE(msd%disp_atom_index(npart),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,failure)
       msd%disp_atom_index = 0
       ALLOCATE(msd%disp_atom_dr(3,npart),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,failure)
       msd%disp_atom_dr = 0.0_dp
       msd%msd_kind = .TRUE.
    END IF
    CALL section_vals_val_get(msd_section,"DISPLACEMENT_TOL",r_val=tol)
    msd%disp_atom_tol = tol*tol

    IF(msd%msd_kind) THEN
      CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds)
      nkind = atomic_kinds%n_els

      ALLOCATE(msd%val_msd_kind(4,nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,failure)
      msd%val_msd_kind = 0.0_dp
    END IF

    IF(msd%msd_molecule) THEN
      CALL cp_subsys_get(subsys=subsys, molecules_new=molecules,&
           molecule_kinds_new=molecule_kinds)
      nmolkind =  molecule_kinds%n_els
      ALLOCATE(msd%val_msd_molecule(4,nmolkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,failure)

      molecule_kind_set => molecule_kinds%els
      molecule_set => molecules%els
      nmol = molecules%n_els

      ALLOCATE(msd%ref0_com_molecule(3,nmol), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,failure)

      DO ikind = 1, nmolkind
        molecule_kind => molecule_kind_set(ikind)
        CALL get_molecule_kind (molecule_kind=molecule_kind, nmolecule=nmolecule )
        DO imol = 1,nmolecule
           molecule => molecule_set(molecule_kind%molecule_list(imol))
           CALL get_molecule ( molecule=molecule, first_atom = first_atom, last_atom = last_atom )
           com = 0.0_dp
           mass_mol = 0.0_dp
           DO iatom = first_atom, last_atom
              mass = particle_set(iatom)%atomic_kind%mass
              com(1) = com(1) + msd%ref0_pos(1,iatom)*mass
              com(2) = com(2) + msd%ref0_pos(2,iatom)*mass
              com(3) = com(3) + msd%ref0_pos(3,iatom)*mass
              mass_mol = mass_mol+mass
           ENDDO  ! iatom
           msd%ref0_com_molecule(1,molecule_kind%molecule_list(imol)) =  com(1)/mass_mol
           msd%ref0_com_molecule(2,molecule_kind%molecule_list(imol)) =  com(2)/mass_mol
           msd%ref0_com_molecule(3,molecule_kind%molecule_list(imol)) =  com(3)/mass_mol
        END DO  ! imol
      ENDDO ! ikind
    END IF

    IF(msd%msd_region) THEN

    END IF

    RETURN
998 CONTINUE ! end of file
    CALL stop_program(routineN,moduleN,__LINE__,&
         "End of reference positions file reached")
999 CONTINUE ! error
    CALL stop_program(routineN,moduleN,__LINE__,&
         "Error reading reference positions file")

  END SUBROUTINE initialize_msd_reftraj

! *****************************************************************************
!> \brief ...
!> \param reftraj ...
!> \param md_env ...
!> \param particle_set ...
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
  SUBROUTINE compute_msd_reftraj(reftraj,md_env,particle_set)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_msd_reftraj', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom, bo(2), first_atom, iatom, ikind, imol, imol_global, &
      last_atom, mepos, natom_kind, nmol_per_kind, nmolecule, nmolkind, num_pe
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: com(3), diff2_com(4), dr2, &
                                                dx, dy, dz, mass, mass_mol, &
                                                msd_mkind(4), rcom(3)
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    NULLIFY(force_env,para_env,subsys)
    NULLIFY(atomic_kind,atomic_kinds,atom_list)
    NULLIFY(local_molecules, molecule, molecule_kind, molecule_kinds,&
            molecule_kind_set, molecules, molecule_set)

    CALL get_md_env(md_env=md_env, force_env=force_env, para_env=para_env)
    CALL force_env_get(force_env=force_env,subsys=subsys)
    CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

    num_pe = para_env%num_pe
    mepos  = para_env%mepos


    IF(reftraj%msd%msd_kind)THEN
      reftraj%msd%val_msd_kind = 0.0_dp
      reftraj%msd%num_disp_atom =  0
      reftraj%msd%disp_atom_dr =  0.0_dp
! compute com
      rcom = 0.0_dp
      DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               natom=natom_kind, mass=mass)
          bo = get_limit(natom_kind,num_pe,mepos)
          DO iatom = bo(1), bo(2)
             atom = atom_list(iatom)
             rcom(1) = rcom(1) + particle_set(atom)%r(1)*mass
             rcom(2) = rcom(2) + particle_set(atom)%r(2)*mass
             rcom(3) = rcom(3) + particle_set(atom)%r(3)*mass
          END DO
      END DO
      CALL mp_sum(rcom,para_env%group)
      rcom = rcom / reftraj%msd%total_mass
      reftraj%msd%drcom (1) = rcom(1) - reftraj%msd%ref0_com(1)
      reftraj%msd%drcom (2) = rcom(2) - reftraj%msd%ref0_com(2)
      reftraj%msd%drcom (3) = rcom(3) - reftraj%msd%ref0_com(3)
!      IF(para_env%ionode) WRITE(*,'(A,T50,3f10.5)') ' COM displacement (dx,dy,dz) [angstrom]:  ', &
!                         drcom(1)*angstrom,drcom(2)*angstrom,drcom(3)*angstrom
! compute_com

      DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               natom=natom_kind)
          bo = get_limit(natom_kind,num_pe,mepos)
          DO iatom = bo(1), bo(2)
            atom = atom_list(iatom)
            dx =  particle_set(atom)%r(1)-reftraj%msd%ref0_pos(1,atom) - &
                  reftraj%msd%drcom(1)
            dy =  particle_set(atom)%r(2)-reftraj%msd%ref0_pos(2,atom) - &
                  reftraj%msd%drcom(2)
            dz =  particle_set(atom)%r(3)-reftraj%msd%ref0_pos(3,atom) - &
                  reftraj%msd%drcom(3)
            dr2 = dx*dx +  dy*dy + dz*dz

            reftraj%msd%val_msd_kind(1,ikind) = reftraj%msd%val_msd_kind(1,ikind) + dx*dx
            reftraj%msd%val_msd_kind(2,ikind) = reftraj%msd%val_msd_kind(2,ikind) + dy*dy
            reftraj%msd%val_msd_kind(3,ikind) = reftraj%msd%val_msd_kind(3,ikind) + dz*dz
            reftraj%msd%val_msd_kind(4,ikind) = reftraj%msd%val_msd_kind(4,ikind) + dr2

            IF(reftraj%msd%disp_atom) THEN
              IF( dr2 > reftraj%msd%disp_atom_tol ) THEN
                reftraj%msd%num_disp_atom  = reftraj%msd%num_disp_atom  +1
                reftraj%msd%disp_atom_dr(1, atom) = dx
                reftraj%msd%disp_atom_dr(2, atom) = dy
                reftraj%msd%disp_atom_dr(3, atom) = dz
              END IF
            END IF
          END DO  !iatom
          reftraj%msd%val_msd_kind(1:4,ikind) = &
                     reftraj%msd%val_msd_kind(1:4,ikind)/REAL(natom_kind,KIND=dp)

      END DO  ! ikind
    ENDIF
    CALL mp_sum(reftraj%msd%val_msd_kind,para_env%group)
    CALL mp_sum(reftraj%msd%num_disp_atom,para_env%group)
    CALL mp_sum(reftraj%msd%disp_atom_dr,para_env%group)

    IF(reftraj%msd%msd_molecule) THEN
      CALL cp_subsys_get(subsys=subsys, local_molecules_new=local_molecules, &
           molecules_new=molecules, molecule_kinds_new=molecule_kinds)

      nmolkind =  molecule_kinds%n_els
      molecule_kind_set => molecule_kinds%els
      molecule_set => molecules%els

      reftraj%msd%val_msd_molecule = 0.0_dp
      DO ikind = 1,nmolkind
         molecule_kind => molecule_kind_set(ikind)
         CALL get_molecule_kind (molecule_kind=molecule_kind, nmolecule=nmolecule )
         nmol_per_kind =  local_molecules%n_el(ikind)
         msd_mkind = 0.0_dp
         DO imol = 1, nmol_per_kind
            imol_global = local_molecules%list(ikind)%array(imol)
            molecule => molecule_set ( imol_global )
            CALL get_molecule (molecule,first_atom=first_atom,last_atom=last_atom)

            com = 0.0_dp
            mass_mol = 0.0_dp
            DO iatom = first_atom, last_atom
              mass = particle_set(iatom)%atomic_kind%mass
              com(1) = com(1) + particle_set(iatom)%r(1)*mass
              com(2) = com(2) + particle_set(iatom)%r(2)*mass
              com(3) = com(3) + particle_set(iatom)%r(3)*mass
              mass_mol = mass_mol+mass
            ENDDO  ! iatom
            com(1) =  com(1)/mass_mol
            com(2) =  com(2)/mass_mol
            com(3) =  com(3)/mass_mol
            diff2_com(1) = com(1)- reftraj%msd%ref0_com_molecule(1,imol_global)
            diff2_com(2) = com(2)- reftraj%msd%ref0_com_molecule(2,imol_global)
            diff2_com(3) = com(3)- reftraj%msd%ref0_com_molecule(3,imol_global)
            diff2_com(1) = diff2_com(1)*diff2_com(1)
            diff2_com(2) = diff2_com(2)*diff2_com(2)
            diff2_com(3) = diff2_com(3)*diff2_com(3)
            diff2_com(4) = diff2_com(1) + diff2_com(2) + diff2_com(3)
            msd_mkind(1) = msd_mkind(1) + diff2_com(1)
            msd_mkind(2) = msd_mkind(2) + diff2_com(2)
            msd_mkind(3) = msd_mkind(3) + diff2_com(3)
            msd_mkind(4) = msd_mkind(4) + diff2_com(4)
         ENDDO ! imol

         reftraj%msd%val_msd_molecule(1,ikind) = msd_mkind(1)/REAL(nmolecule,KIND=dp)
         reftraj%msd%val_msd_molecule(2,ikind) = msd_mkind(2)/REAL(nmolecule,KIND=dp)
         reftraj%msd%val_msd_molecule(3,ikind) = msd_mkind(3)/REAL(nmolecule,KIND=dp)
         reftraj%msd%val_msd_molecule(4,ikind) = msd_mkind(4)/REAL(nmolecule,KIND=dp)
      END DO  ! ikind
      CALL mp_sum(reftraj%msd%val_msd_molecule, para_env%group)

    END IF

  END SUBROUTINE compute_msd_reftraj

! *****************************************************************************
!> \brief ...
!> \param md_env ...
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
  SUBROUTINE write_output_reftraj(md_env)
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(len=*), PARAMETER :: routineN = 'write_output_reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_act, my_mittle, my_pos
    INTEGER                                  :: iat, ikind, nkind, out_msd
    LOGICAL, SAVE                            :: first_entry = .FALSE.
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section, root_section

    NULLIFY(logger)
    logger => cp_get_default_logger()

    NULLIFY(reftraj)
    NULLIFY(reftraj_section,root_section)

    CALL get_md_env(md_env=md_env, force_env=force_env, &
         reftraj=reftraj)

    CALL force_env_get(force_env=force_env,root_section=root_section)

    reftraj_section => section_vals_get_subs_vals(root_section,&
           "MOTION%MD%REFTRAJ")

    my_pos = "APPEND"
    my_act = "WRITE"

    IF(reftraj%init.AND.(reftraj%isnap==reftraj%info%first_snapshot)) THEN
      my_pos  = "REWIND"
      first_entry = .TRUE.
    END IF

    IF(reftraj%info%msd) THEN
      IF(reftraj%msd%msd_kind)THEN
         nkind = SIZE(reftraj%msd%val_msd_kind,2)
         DO ikind = 1,nkind
           my_mittle ="k"//TRIM(ADJUSTL(cp_to_string(ikind)))
           out_msd = cp_print_key_unit_nr(logger,reftraj_section,"PRINT%MSD_KIND",&
               extension=".msd", file_position=my_pos, file_action=my_act,&
               file_form="FORMATTED", middle_name=TRIM(my_mittle))
           IF(out_msd>0) THEN
              WRITE(UNIT=out_msd,FMT="(I8, F12.3,4F20.10)") reftraj%itimes, &
                   reftraj%time*femtoseconds, &
                   reftraj%msd%val_msd_kind(1:4,ikind)*angstrom*angstrom
              CALL m_flush(out_msd)
           END IF
           CALL cp_print_key_finished_output(out_msd,logger,reftraj_section,&
                "PRINT%MSD_KIND")
         END DO
      END IF
      IF(reftraj%msd%msd_molecule) THEN
         nkind = SIZE(reftraj%msd%val_msd_molecule,2)
         DO ikind = 1,nkind
           my_mittle ="mk"//TRIM(ADJUSTL(cp_to_string(ikind)))
           out_msd = cp_print_key_unit_nr(logger,reftraj_section,"PRINT%MSD_MOLECULE",&
               extension=".msd", file_position=my_pos, file_action=my_act,&
               file_form="FORMATTED", middle_name=TRIM(my_mittle))
           IF(out_msd>0) THEN
              WRITE(UNIT=out_msd,FMT="(I8, F12.3,4F20.10)") reftraj%itimes, &
                   reftraj%time*femtoseconds, &
                   reftraj%msd%val_msd_molecule(1:4,ikind)*angstrom*angstrom
              CALL m_flush(out_msd)
           END IF
           CALL cp_print_key_finished_output(out_msd,logger,reftraj_section,&
                "PRINT%MSD_MOLECULE")
         END DO
      END IF
      IF(reftraj%msd%disp_atom) THEN

         IF(first_entry)  my_pos = "REWIND"
         my_mittle ="disp_at"
         out_msd = cp_print_key_unit_nr(logger,reftraj_section,"PRINT%DISPLACED_ATOM",&
               extension=".msd", file_position=my_pos, file_action=my_act,&
               file_form="FORMATTED", middle_name=TRIM(my_mittle))
         IF(out_msd>0 .AND. reftraj%msd%num_disp_atom>0) THEN
            IF(first_entry) THEN
               first_entry = .FALSE.
            END IF
            WRITE(UNIT=out_msd,FMT="(A,T7,I8, A, T29, F12.3, A, T50, I10)") "# i = ", reftraj%itimes, "  time (fs) = ",&
                   reftraj%time*femtoseconds, "  nat = ", reftraj%msd%num_disp_atom
            DO iat = 1,SIZE(reftraj%msd%disp_atom_dr,2)
               IF(ABS(reftraj%msd%disp_atom_dr(1,iat)) > 0.0_dp) THEN
                 WRITE(UNIT=out_msd,FMT="(I8, 3F20.10)") iat,& !reftraj%msd%disp_atom_index(iat),&
                      reftraj%msd%disp_atom_dr(1,iat)*angstrom, &
                      reftraj%msd%disp_atom_dr(2,iat)*angstrom,&
                      reftraj%msd%disp_atom_dr(3,iat)*angstrom
               END IF
            END DO
         ENDIF
         CALL cp_print_key_finished_output(out_msd,logger,reftraj_section,&
                "PRINT%DISPLACED_ATOM")
      END IF
    ENDIF ! msd
    reftraj%init = .FALSE.

  END SUBROUTINE write_output_reftraj

END MODULE reftraj_util

