!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains miscellaneous subroutines used in the Monte Carlo runs,mostly
!>      geared towards changes in coordinates
!> \author MJM
! *****************************************************************************
MODULE mc_coordinates
  USE cell_methods,                    ONLY: cell_type,&
                                             get_cell
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: get_mc_molecule_info,&
                                             mc_molecule_info_type
  USE message_passing,                 ONLY: mp_bcast
  USE molecule_types_new,              ONLY: molecule_type
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PRIVATE :: generate_avbmc_insertion

  PUBLIC :: generate_cbmc_swap_config,&
            get_center_of_mass,mc_coordinate_fold,&
            find_mc_test_molecule,&
            create_discrete_array,&
            check_for_overlap,&
            rotate_molecule

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_coordinates'

CONTAINS

! *****************************************************************************
!> \brief looks for overlaps (intermolecular distances less than rmin)
!> \param force_env the force environment containing the coordinates
!> \param nchains the number of molecules of each type in the box
!> \param nunits the number of interaction sites for each molecule
!> \param loverlap returns .TRUE. if atoms overlap
!> \param mol_type an array that indicates the type of each molecule
!> \param cell_length the length of the box...if none is specified,
!>                     it uses the cell found in the force_env
!> \param molecule_number if present, just look for overlaps with this
!>            molecule
!>
!>      Suitable for parallel use.
!> \author MJM
! *****************************************************************************
SUBROUTINE check_for_overlap ( force_env, nchains, nunits,loverlap,mol_type,&
                                 cell_length,molecule_number)

    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, DIMENSION(:), INTENT(IN)        :: nchains, nunits
    LOGICAL, INTENT(OUT)                     :: loverlap
    INTEGER, DIMENSION(:), INTENT(IN)        :: mol_type
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN), OPTIONAL                   :: cell_length
    INTEGER, INTENT(IN), OPTIONAL            :: molecule_number

    CHARACTER(len=*), PARAMETER :: routineN = 'check_for_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imol, istat, iunit, &
                                                jmol, jstart, junit, nend, &
                                                nstart, nunit, nunits_i, &
                                                nunits_j
    LOGICAL                                  :: lall
    REAL(KIND=dp)                            :: dist, rmin
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r
    REAL(KIND=dp), DIMENSION(1:3)            :: abc, box_length, RIJ
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsys_type), POINTER            :: oldsys
    TYPE(particle_list_type), POINTER        :: particles

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(oldsys,particles)

! initialize some stuff
      loverlap=.FALSE.
      rmin=3.57106767_dp ! 1 angstrom squared

! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(oldsys,particles=particles, error=error)

      ALLOCATE(r(1:3,1:MAXVAL(nunits),1:SUM(nchains)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "mass",3*MAXVAL(nunits)*SUM(nchains)*dp_size)

      IF( PRESENT(cell_length)) THEN
         box_length(1:3)=cell_length(1:3)
      ELSE
         box_length(1:3)=abc(1:3)
      ENDIF

! put the coordinates into an easier matrix to manipulate
      junit=0
      DO imol=1,SUM(nchains)
         nunit=nunits(mol_type(imol))
         DO iunit=1,nunit
            junit=junit+1
            r(1:3,iunit,imol)=particles%els(junit)%r(1:3)
         ENDDO
      ENDDO

! now let's find the LJ energy between all the oxygens and
! the charge interactions
      lall=.TRUE.
      jstart=1
      IF( PRESENT(molecule_number)) THEN
         lall=.FALSE.
         nstart=molecule_number
         nend=molecule_number
      ELSE
         nstart=1
         nend=SUM(nchains(:))
      ENDIF
      DO imol=nstart,nend
         IF (lall) jstart=imol+1
         nunits_i=nunits(mol_type(imol))
         DO jmol=jstart,SUM(nchains(:))
            IF(imol == jmol) CYCLE
            nunits_j=nunits(mol_type(jmol))

            DO iunit=1,nunits_i
               DO junit=1,nunits_j
! find the minimum image distance
                  RIJ(1)=r(1,iunit,imol)-r(1,junit,jmol)-&
                     box_length(1)*ANINT(&
                     (r(1,iunit,imol)-r(1,junit,jmol))/box_length(1))
                  RIJ(2)=r(2,iunit,imol)-r(2,junit,jmol)-&
                     box_length(2)*ANINT(&
                     (r(2,iunit,imol)-r(2,junit,jmol))/box_length(2))
                  RIJ(3)=r(3,iunit,imol)-r(3,junit,jmol)-&
                     box_length(3)*ANINT(&
                     (r(3,iunit,imol)-r(3,junit,jmol))/box_length(3))

                  dist=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2

                  IF (dist < rmin ) THEN
                     loverlap=.TRUE.
                     DEALLOCATE(r,STAT=istat)
                     IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
                     __LINE__,"r")

                     CALL timestop(handle)
                     RETURN
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"r")

! end the timing
      CALL timestop(handle)

END SUBROUTINE check_for_overlap

! *****************************************************************************
!> \brief calculates the center of mass of a given molecule
!> \param coordinates the coordiantes of the atoms in the molecule
!> \param natom the number of atoms in the molecule
!> \param center_of_mass the coordinates of the center of mass
!> \param mass the mass of the atoms in the molecule
!>
!>    Designed for parallel use.
!> \author MJM
! *****************************************************************************
SUBROUTINE get_center_of_mass ( coordinates,natom,center_of_mass,&
            mass)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coordinates
    INTEGER, INTENT(IN)                      :: natom
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(OUT)                            :: center_of_mass
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: mass

    CHARACTER(len=*), PARAMETER :: routineN = 'get_center_of_mass', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom
    REAL(KIND=dp)                            :: total_mass

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      total_mass=SUM(mass(1:natom))
      center_of_mass(:)=0.0E0_dp

      DO iatom=1,natom
         DO i=1,3
            center_of_mass(i)=center_of_mass(i)+&
               mass(iatom)*coordinates(i,iatom)
         ENDDO
      ENDDO

      center_of_mass(1:3)=center_of_mass(1:3)/total_mass

! end the timing
      CALL timestop(handle)

END SUBROUTINE get_center_of_mass

! *****************************************************************************
!> \brief folds all the coordinates into the center simulation box using
!>      a center of mass cutoff
!> \param coordinates the coordiantes of the atoms in the system
!> \param nchains_tot the total number of molecules in the box
!> \param mol_type an array that indicates the type of every molecule in the box
!> \param mass the mass of every atom for all molecule types
!> \param nunits the number of interaction sites for each molecule type
!> \param box_length an array for the lengths of the simulation box sides
!>
!>      Designed for parallel use.
!> \author MJM
! *****************************************************************************
SUBROUTINE mc_coordinate_fold(coordinates,nchains_tot,mol_type,mass,nunits,box_length)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: nchains_tot
    INTEGER, DIMENSION(:), INTENT(IN)        :: mol_type
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: mass
    INTEGER, DIMENSION(:), INTENT(IN)        :: nunits
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN)                             :: box_length

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_coordinate_fold', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: end_atom, handle, iatom, &
                                                imolecule, jatom, &
                                                molecule_type, natoms, &
                                                start_atom
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! loop over all molecules
      end_atom=0
      DO imolecule=1,nchains_tot
         molecule_type=mol_type(imolecule)
         natoms=nunits(molecule_type)
         start_atom=end_atom+1
         end_atom=start_atom+natoms-1
         CALL get_center_of_mass(coordinates(:,start_atom:end_atom),&
              natoms,center_of_mass(:),mass(:,molecule_type))
         DO iatom=1,natoms
            jatom=iatom+start_atom-1
            coordinates(1,jatom)=coordinates(1,jatom)-&
                 box_length(1)*FLOOR(center_of_mass(1)/box_length(1))
            coordinates(2,jatom)=coordinates(2,jatom)-&
                 box_length(2)*FLOOR(center_of_mass(2)/box_length(2))
            coordinates(3,jatom)=coordinates(3,jatom)-&
                 box_length(3)*FLOOR(center_of_mass(3)/box_length(2))
         ENDDO

      ENDDO

! end the timing
      CALL timestop(handle)

END SUBROUTINE mc_coordinate_fold

! *****************************************************************************
!> \brief takes the last molecule in a force environment and moves it around
!>      to different center of mass positions and orientations, selecting one
!>      based on the rosenbluth weight
!> \param force_env the force environment containing the coordinates
!> \param BETA the value of 1/kT for this simulations, in a.u.
!> \param max_val ...
!> \param min_val ...
!> \param exp_max_val ...
!> \param exp_min_val ...
!> \param nswapmoves the number of desired trial configurations
!> \param rosenbluth_weight the Rosenbluth weight for this set of configs
!> \param start_atom the atom number that the molecule to be swapped starts on
!> \param natoms_tot the total number of interaction sites in the box
!> \param nunits the number of interaction sites for every molecule_type
!> \param nunits_mol ...
!> \param mass the mass for every interaction site of every molecule type
!> \param loverlap the flag that indicates if all of the configs have an
!>                  atomic overlap
!> \param choosen_energy the energy of the choosen config
!> \param old_energy the energy that we subtract from all of the trial
!>        energies to prevent numerical overflows
!> \param ionode indicates if we're on the main CPU
!> \param lremove is this the Rosenbluth weight for a removal box?
!> \param mol_type an array that contains the molecule type for every atom in the box
!> \param nchains the number of molecules of each type in this box
!> \param source the MPI source value, for broadcasts
!> \param group the MPI group value, for broadcasts
!> \param rng_stream the random number stream that we draw from
!> \param error the the cp_error type
!> \param avbmc_atom ...
!> \param rmin ...
!> \param rmax ...
!> \param move_type ...
!> \param target_atom ...
!> \param nunits_tot the number of interaction sites in the choosen molecule UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par Optional Avbmc Flags
!>      - avbmc_atom: the atom number that serves for the target atom in each
!>        molecule (1 is the first atom in the molecule, etc.)
!>      - rmin: the minimum AVBMC radius for the shell around the target
!>      - rmax: the maximum AVBMC radius for the shell around the target
!>      - move_type: generate configs in the "in" or "out" volume
!>      - target_atom: the number of the avbmc atom in the target molecule
!> \par
!>      Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE generate_cbmc_swap_config( force_env, BETA, max_val, min_val, exp_max_val,&
     exp_min_val, nswapmoves, rosenbluth_weight, start_atom,natoms_tot,nunits,nunits_mol,&
     mass,loverlap, choosen_energy,old_energy,ionode,lremove,mol_type,nchains,source,&
     group,rng_stream,error,avbmc_atom,rmin,rmax,move_type,target_atom)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(IN)                :: BETA, max_val, min_val, &
                                                exp_max_val, exp_min_val
    INTEGER, INTENT(IN)                      :: nswapmoves
    REAL(KIND=dp), INTENT(OUT)               :: rosenbluth_weight
    INTEGER, INTENT(IN)                      :: start_atom, natoms_tot
    INTEGER, DIMENSION(:), INTENT(IN)        :: nunits
    INTEGER, INTENT(IN)                      :: nunits_mol
    REAL(dp), DIMENSION(1:nunits_mol), &
      INTENT(IN)                             :: mass
    LOGICAL, INTENT(OUT)                     :: loverlap
    REAL(KIND=dp), INTENT(OUT)               :: choosen_energy
    REAL(KIND=dp), INTENT(IN)                :: old_energy
    LOGICAL, INTENT(IN)                      :: ionode, lremove
    INTEGER, DIMENSION(:), INTENT(IN)        :: mol_type, nchains
    INTEGER, INTENT(IN)                      :: source, group
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type)                      :: error
    INTEGER, INTENT(IN), OPTIONAL            :: avbmc_atom
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: rmin, rmax
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: move_type
    INTEGER, INTENT(IN), OPTIONAL            :: target_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_cbmc_swap_config', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_number, choosen, &
                                                end_atom, handle, i, iatom, &
                                                imolecule, imove, istat, &
                                                molecule_number
    LOGICAL                                  :: all_overlaps
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: loverlap_array
    REAL(KIND=dp)                            :: bias_energy, exponent, rand, &
                                                total_running_weight
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: boltz_weights, trial_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r
    REAL(KIND=dp), DIMENSION(1:3)            :: abc, center_of_mass, diff, &
                                                r_insert
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: oldsys
    TYPE(particle_list_type), POINTER        :: particles

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(oldsys)
! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(oldsys,particles=particles, error=error)

! do some checking to make sure we have all the data we need
      IF(PRESENT(avbmc_atom)) THEN
         IF( .NOT. PRESENT(rmin) .OR. .NOT. PRESENT(rmax) .OR. &
            .NOT. PRESENT(move_type) .OR. .NOT. PRESENT(target_atom)) THEN
           CALL stop_program(routineN,moduleN,__LINE__,&
               'AVBMC swap move is missing information!')
         ENDIF
      ENDIF

      ALLOCATE(r_old(1:3,1:natoms_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r_old",3*natoms_tot*dp_size)
      ALLOCATE(r(1:3,1:natoms_tot,1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r",3*natoms_tot*nswapmoves*dp_size)
      ALLOCATE(trial_energy(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "trial_energy",nswapmoves*dp_size)
      ALLOCATE(boltz_weights(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "boltz_weights",nswapmoves*dp_size)
      ALLOCATE(loverlap_array(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "loverlap_array",nswapmoves)

! initialize the arrays that need it
      loverlap_array(:)=.FALSE.
      loverlap=.FALSE.
      boltz_weights(:)  = 0.0_dp
      trial_energy(:)   = 0.0_dp
      r(:,:,:)          = 0.0_dp
      choosen_energy    = 0.0_dp
      rosenbluth_weight = 0.0_dp

! save the positions of the molecules
      DO imove=1,nswapmoves
         DO iatom=1,natoms_tot
            r(1:3,iatom,imove)=particles%els(iatom)%r(1:3)
         ENDDO
      ENDDO

! save the remove coordinates
      DO iatom=1,natoms_tot
         r_old(1:3,iatom)=r(1:3,iatom,1)
      ENDDO

! figure out the numbers of the first and last atoms in the molecule
      end_atom=start_atom+nunits_mol-1
! figure out which molecule number we're on
      molecule_number=0
      atom_number=1
      DO imolecule=1,SUM(nchains(:))
         IF(atom_number == start_atom) THEN
            molecule_number=imolecule
            EXIT
         ENDIF
         atom_number=atom_number+nunits(mol_type(imolecule))
      ENDDO
      IF(molecule_number == 0) CALL stop_program(routineN,moduleN,__LINE__,&
               'CBMC swap move cannot find which molecule number it needs')

      IF(lremove) THEN
         CALL check_for_overlap (force_env,nchains,nunits,loverlap_array(1),&
              mol_type)
         CALL mp_bcast(loverlap_array(1),source,group)

         IF(loverlap_array(1)) THEN
            IF(ionode) THEN
               WRITE(6,*) start_atom,end_atom,natoms_tot
               DO iatom=1,natoms_tot
                  WRITE(6,*) r(1:3,iatom,1)
               ENDDO
            ENDIF
            CALL stop_program(routineN,moduleN,__LINE__,&
               'CBMC swap move found an overlap in the old config')
         ENDIF
      ENDIF

      DO imove=1,nswapmoves

! drop into serial
         IF(ionode) THEN

            IF(PRESENT(avbmc_atom)) THEN
! find an AVBMC insertion point
               CALL generate_avbmc_insertion(rmin,rmax,&
                  r_old(1:3,target_atom),&
                  move_type,r_insert(:),abc(:),rng_stream,error=error)

               DO i=1,3
                  diff(i)=r_insert(i)-r_old(i,start_atom+avbmc_atom-1)
               ENDDO

            ELSE
! find a new insertion point somewhere in the box
               DO i=1,3
                  rand=next_random_number(rng_stream,error=error)
                  r_insert(i)=rand*abc(i)
               ENDDO

! find the center of mass of the insertion molecule
               CALL get_center_of_mass(r(:,start_atom:end_atom,imove),nunits_mol,&
               center_of_mass(:),mass(:))

! move the molecule to the insertion point

               DO i=1,3
                  diff(i)=r_insert(i)-center_of_mass(i)
               ENDDO

            ENDIF

            DO iatom=start_atom,end_atom
               r(1:3,iatom,imove)=r(1:3,iatom,imove)+diff(1:3)
            ENDDO

! rotate the molecule...this routine is only made for serial use
            CALL rotate_molecule(r(:,start_atom:end_atom,imove),mass(:),&
                 nunits_mol,rng_stream,error=error)

            IF(imove == 1 .AND. lremove) THEN
               DO iatom=1,natoms_tot
                  r(1:3,iatom,1)=r_old(1:3,iatom)
               ENDDO
            ENDIF

         ENDIF

         CALL mp_bcast(r(:,:,imove),source,group)

! calculate the energy and boltzman weight of the config
         DO iatom=start_atom,end_atom
            particles%els(iatom)%r(1:3)=r(1:3,iatom,imove)
         ENDDO

         CALL check_for_overlap (force_env,nchains,nunits,loverlap_array(imove),&
              mol_type,molecule_number=molecule_number)
         IF(loverlap_array(imove)) THEN
            boltz_weights(imove)=0.0_dp
            CYCLE
         ENDIF

         CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,error=error)
         CALL force_env_get(force_env,&
            potential_energy=bias_energy,error=error)

         trial_energy(imove)=(bias_energy-old_energy)
         exponent=-BETA*trial_energy(imove)

         IF    (exponent .GT. exp_max_val) THEN
            boltz_weights(imove)=max_val
         ELSEIF(exponent .LT. exp_min_val) THEN
            boltz_weights(imove)=min_val
         ELSE
            boltz_weights(imove)=EXP(exponent)
         ENDIF

      ENDDO

! now we need to pick a configuration based on the Rosenbluth weight,
! which is just the sum of the Boltzmann weights
      rosenbluth_weight=SUM(boltz_weights(:))
      IF(rosenbluth_weight .EQ. 0.0_dp .AND. lremove) THEN
! should never have 0.0 for an old weight...causes a divide by zero
! in the acceptance rule
         IF(ionode) THEN
            WRITE(6,*) boltz_weights(1:nswapmoves)
            WRITE(6,*) start_atom,end_atom,lremove
            WRITE(6,*) loverlap_array(1:nswapmoves)
            WRITE(6,*) natoms_tot
            WRITE(6,*)
            DO iatom=1,natoms_tot
               WRITE(6,*) r(1:3,iatom,1)*angstrom
            ENDDO
         ENDIF
         CALL stop_program(routineN,moduleN,__LINE__,&
              'CBMC swap move found a bad old weight')
      ENDIF
      all_overlaps=.TRUE.
      total_running_weight=0.0E0_dp
      choosen=0
      IF(ionode)THEN
         rand=next_random_number(rng_stream,error=error)
!         CALL random_number(rand)
      ENDIF
      CALL mp_bcast(rand,source,group)
      DO imove=1,nswapmoves
         IF ( loverlap_array(imove) ) CYCLE
         all_overlaps=.FALSE.
         total_running_weight=total_running_weight+boltz_weights(imove)
         IF ( total_running_weight .GE. rand*rosenbluth_weight) THEN
            choosen=imove
            EXIT
         ENDIF
      ENDDO

      IF(all_overlaps) THEN
         loverlap = .TRUE.

! if this is an old configuration, we always choose the first one...
! this should never be the case, but I'm testing something
         IF(lremove) THEN
            IF(ionode) THEN
               WRITE(6,*) boltz_weights(1:nswapmoves)
               WRITE(6,*) start_atom,end_atom,lremove
               WRITE(6,*) loverlap_array(1:nswapmoves)
               DO iatom=1,natoms_tot
                  WRITE(6,*) r(1:3,iatom,1)
               ENDDO
            ENDIF
            CALL stop_program(routineN,moduleN,__LINE__,&
            'CBMC swap move found all overlaps for the remove config')
         ENDIF

         DEALLOCATE(r_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"r_old")
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"r")
         DEALLOCATE(trial_energy,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"trial_energy")
         DEALLOCATE(boltz_weights,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"boltz_weights")
         DEALLOCATE(loverlap_array,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"loverlap_array")
         CALL timestop(handle)
         RETURN
      ENDIF

! make sure a configuration was choosen
      IF(choosen == 0) CALL stop_program(routineN,moduleN,&
         __LINE__,'CBMC swap move failed to select config')

! if this is an old configuration, we always choose the first one
      IF(lremove) choosen=1

! set the energy for the configuration
      choosen_energy=trial_energy(choosen)

! copy the coordinates to the force environment
      DO iatom=1,natoms_tot
         particles%els(iatom)%r(1:3)=r(1:3,iatom,choosen)
      ENDDO

      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"r_old")
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"r")
      DEALLOCATE(trial_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"trial_energy")
      DEALLOCATE(boltz_weights,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"boltz_weights")
      DEALLOCATE(loverlap_array,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"loverlap_array")

! end the timing
      CALL timestop(handle)

END SUBROUTINE generate_cbmc_swap_config

! *****************************************************************************
!> \brief rotates a molecule randomly around the center of mass,
!>      sequentially in x, y, and z directions
!> \param r the coordinates of the molecule to rotate
!> \param mass the mass of all the atoms in the molecule
!> \param natoms the number of atoms in the molecule
!> \param rng_stream the stream we pull random numbers from
!>
!>    Use only in serial.
!> \param error ...
!> \author MJM
! *****************************************************************************
SUBROUTINE rotate_molecule(r,mass,natoms,rng_stream,error)

    INTEGER, INTENT(IN)                      :: natoms
    REAL(KIND=dp), DIMENSION(1:natoms), &
      INTENT(IN)                             :: mass
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:natoms), &
      INTENT(INOUT)                          :: r
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotate_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iunit
    REAL(KIND=dp)                            :: cosdg, dgamma, rand, rx, &
                                                rxnew, ry, rynew, rz, rznew, &
                                                sindg
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! find the center of mass of the molecule
      CALL get_center_of_mass(r(:,:),natoms,center_of_mass(:),mass(:))

! call a random number to figure out how far we're moving
      rand=next_random_number(rng_stream,error=error)
      dgamma=pi*(rand-0.5E0_dp)*2.0E0_dp

! *** set up the rotation matrix ***

      cosdg = COS( dgamma )
      sindg = SIN( dgamma )

! ***    ROTATE UNITS OF I AROUND X-AXIS ***

      DO  iunit = 1,natoms
         ry = r(2,iunit) - center_of_mass(2)
         rz = r(3,iunit) - center_of_mass(3)
         rynew = cosdg * ry + sindg * rz
         rznew = cosdg * rz - sindg * ry

         r(2,iunit) = rynew + center_of_mass(2)
         r(3,iunit) = rznew + center_of_mass(3)

      ENDDO

! ***    ROTATE UNITS OF I AROUND y-AXIS ***

      DO  iunit = 1,natoms
         rx = r(1,iunit) - center_of_mass(1)
         rz = r(3,iunit) - center_of_mass(3)
         rxnew = cosdg * rx + sindg * rz
         rznew = cosdg * rz - sindg * rx

         r(1,iunit) = rxnew + center_of_mass(1)
         r(3,iunit) = rznew + center_of_mass(3)

      ENDDO

! ***    ROTATE UNITS OF I AROUND z-AXIS ***

      DO  iunit = 1,natoms
         rx = r(1,iunit) - center_of_mass(1)
         ry = r(2,iunit) - center_of_mass(2)
         rxnew = cosdg * rx + sindg * ry
         rynew = cosdg * ry - sindg * rx

         r(1,iunit) = rxnew + center_of_mass(1)
         r(2,iunit) = rynew + center_of_mass(2)

      ENDDO

! end the timing
      CALL timestop(handle)

END SUBROUTINE rotate_molecule

! *****************************************************************************
!> \brief selects a molecule at random to perform a MC move on...you can specify
!>      the box the molecule should be in, its type, both, or neither
!> \param mc_molecule_info the structure that contains some global molecule data
!> \param start_atom the number of the first atom in the choosen molecule in relation
!>        to the force_env it's in
!> \param box_number the box the choosen molecule is in
!> \param molecule_type the type of molecule the choosen molecule is
!> \param rng_stream the stream we pull random numbers from
!> \param box if present, tells the routine which box to grab a molecule from
!> \param molecule_type_old if present, tells the routine which molecule type to select from
!> \param error the optional cp_error_type
!>
!>    Designed only for serial use.
!> \author MJM
! *****************************************************************************
SUBROUTINE find_mc_test_molecule(mc_molecule_info,start_atom,&
     box_number,molecule_type,rng_stream,box,molecule_type_old,error)

    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    INTEGER, INTENT(OUT)                     :: start_atom, box_number, &
                                                molecule_type
    TYPE(rng_stream_type), POINTER           :: rng_stream
    INTEGER, INTENT(IN), OPTIONAL            :: box, molecule_type_old
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_mc_test_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ibox, imol_type, &
                                                imolecule, jbox, &
                                                molecule_number, nchains_tot, &
                                                start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rand

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(nunits,mol_type,nchains)
      CALL get_mc_molecule_info(mc_molecule_info,nchains=nchains,nunits=nunits,&
           mol_type=mol_type)

! initialize the outgoing variables
      start_atom=0
      box_number=0
      molecule_type=0

      IF(PRESENT(box) .AND. PRESENT(molecule_type_old))THEN
! only need to find the atom number the molecule starts on
         rand=next_random_number(rng_stream,error=error)
         molecule_number=CEILING(rand*REAL(nchains(molecule_type_old,box)))

         start_mol=1
         DO jbox=1,box-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO

! adjust to take into account molecules of other types in the box
         DO imol_type=1,molecule_type_old-1
            molecule_number=molecule_number+nchains(imol_type,box)
         ENDDO

         start_atom=1
         DO imolecule=1,molecule_number-1
            start_atom=start_atom+nunits(mol_type(start_mol+imolecule-1))
         ENDDO

      ELSEIF(PRESENT(box))THEN
! any molecule in box...need to find molecule type and start atom
         rand=next_random_number(rng_stream,error=error)
         molecule_number=CEILING(rand*REAL(SUM(nchains(:,box))))

         start_mol=1
         DO jbox=1,box-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO

         molecule_type=mol_type(start_mol+molecule_number-1)

! now the starting atom
         start_atom=1
         DO imolecule=1,molecule_number-1
            start_atom=start_atom+nunits(mol_type(start_mol+imolecule-1))
         ENDDO

      ELSEIF(PRESENT(molecule_type_old))THEN
! any molecule of type molecule_type_old...need to find box number and start atom
         rand=next_random_number(rng_stream,error=error)
         molecule_number=CEILING(rand*REAL(SUM(nchains(molecule_type_old,:))))

! find which box it's in
         nchains_tot=0
         DO ibox=1,SIZE(nchains(molecule_type_old,:))
            IF(molecule_number .LE. nchains(molecule_type_old,ibox)) THEN
               box_number=ibox
               EXIT
            ENDIF
            molecule_number=molecule_number-nchains(molecule_type_old,ibox)
         ENDDO

         start_mol=1
         DO jbox=1,box_number-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO

! now find the starting atom number
         DO imol_type=1,molecule_type_old-1
            molecule_number=molecule_number+nchains(imol_type,box_number)
         ENDDO
         start_atom=1
         DO imolecule=1,molecule_number-1
            start_atom=start_atom+nunits(mol_type(start_mol+imolecule-1))
         ENDDO

      ELSE
! no restrictions...need to find all pieces of data
         nchains_tot=0
         DO ibox=1,SIZE(nchains(1,:))
            nchains_tot=nchains_tot+SUM(nchains(:,ibox))
         ENDDO
         rand=next_random_number(rng_stream,error=error)
         molecule_number=CEILING(rand*REAL(nchains_tot))

         molecule_type=mol_type(molecule_number)

! now which box it's in
         DO ibox=1,SUM(nchains(1,:))
            IF(molecule_number .LE. SUM(nchains(:,ibox))) THEN
               box_number=ibox
               EXIT
            ENDIF
            molecule_number=molecule_number-SUM(nchains(:,ibox))
         ENDDO

! now find the starting atom number
         start_mol=1
         DO jbox=1,box_number-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO
         start_atom=1
         DO imolecule=1,molecule_number-1
            start_atom=start_atom+nunits(mol_type(start_mol+imolecule-1))
         ENDDO

      ENDIF

! make sure things are good
      IF(PRESENT(box)) box_number=box
      IF(PRESENT(molecule_type_old)) molecule_type=molecule_type_old

      CPPostcondition(start_atom/=0,cp_failure_level,routineP,error,failure)
      CPPostcondition(box_number/=0,cp_failure_level,routineP,error,failure)
      CPPostcondition(molecule_type/=0,cp_failure_level,routineP,error,failure)

! end the timing
      CALL timestop(handle)

END SUBROUTINE find_mc_test_molecule

! *****************************************************************************
!> \brief generates an array that tells us which sides of the simulation
!>      cell we can increase or decrease using a discrete volume move
!> \param cell the lengths of the sides of the cell
!> \param discrete_array the array that indicates which sides we can move
!> \param step_size the size of the discrete volume move
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE create_discrete_array(cell,discrete_array,step_size)

! 1 is for increase, 2 is for decrease
! 1 is for "yes, we can do the move", 0 is for no

    REAL(dp), DIMENSION(1:3), INTENT(IN)     :: cell
    INTEGER, DIMENSION(1:3, 1:2), &
      INTENT(OUT)                            :: discrete_array
    REAL(dp), INTENT(IN)                     :: step_size

    INTEGER                                  :: iside
    REAL(dp)                                 :: high_value, length1, length2, &
                                                low_value

      discrete_array(:,:)=0

      length1=ABS(cell(1)-cell(2))
      length2=ABS(cell(2)-cell(3))

! now let's figure out all the different cases
      IF(length1 .LT. 0.01_dp*step_size .AND. &
         length2 .LT. 0.01_dp*step_size) THEN
! all sides are equal, so we can move up or down
         discrete_array(1:3,1)=1
         discrete_array(1:3,2)=1
      ELSE

! find the low value and the high value
         high_value=-1.0_dp
         low_value=cell(1)*cell(2)*cell(3)
         DO iside=1,3
            IF(cell(iside) .LT. low_value) low_value=cell(iside)
            IF(cell(iside) .GT. high_value) high_value=cell(iside)
         ENDDO
         DO iside=1,3
! now we see if the value is a high value or a low value...it can only be
! one of the two
            IF(ABS(cell(iside)-low_value) .LT. 0.01_dp*step_size) THEN
! low value, we can only increase the cell size
               discrete_array(iside,1)=1
               discrete_array(iside,2)=0
            ELSE
! high value, we can only decrease the cell size
               discrete_array(iside,1)=0
               discrete_array(iside,2)=1
            ENDIF
         ENDDO
      ENDIF

END SUBROUTINE create_discrete_array

! *****************************************************************************
!> \brief generates an insertion point in either the "in" or the "out" volume
!>      of a target atom, where the "in" volume is a shell with inner radius
!>       rmin and outer radius rmax
!> \param rmin the minimum AVBMC radius for the shell around the target
!> \param rmax the maximum AVBMC radius for the shell around the target
!> \param r_target the coordinates of the target atom
!> \param move_type generate configs in the "in" or "out" volume
!> \param r_insert the output insertion site
!> \param abc the lengths of the sides of the simulation box
!> \param rng_stream the random number stream that we draw from
!>
!>      Use only in serial.
!> \param error ...
!> \author MJM
! *****************************************************************************
SUBROUTINE generate_avbmc_insertion( rmin,rmax,r_target,&
      move_type,r_insert,abc,rng_stream,error)

    REAL(KIND=dp), INTENT(IN)                :: rmin, rmax
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN)                             :: r_target
    CHARACTER(LEN=*), INTENT(IN)             :: move_type
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(OUT)                            :: r_insert
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN)                             :: abc
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i
    REAL(dp)                                 :: dist, eta_1, eta_2, eta_sq, &
                                                rand
    REAL(dp), DIMENSION(1:3)                 :: RIJ

      r_insert(1:3)=0.0_dp

      IF(move_type == 'in') THEN
! generate a random unit vector, from Allen and Tildesley
         DO
            eta_1=next_random_number(rng_stream,error=error)
            eta_2=next_random_number(rng_stream,error=error)
            eta_sq=eta_1**2+eta_2**2
            IF(eta_sq .LT. 1.0_dp) THEN
               r_insert(1)=2.0_dp*eta_1*SQRT(1.0_dp-eta_sq)
               r_insert(2)=2.0_dp*eta_2*SQRT(1.0_dp-eta_sq)
               r_insert(3)=1.0_dp-2.0_dp*eta_sq
               EXIT
            ENDIF
         ENDDO

! now scale that vector to be within the "in" region
         rand=next_random_number(rng_stream,error=error)
         r_insert(1:3)=r_insert(1:3)*(rand*(rmax**3-rmin**3)+rmin**3)**&
            (1.0_dp/3.0_dp)

         r_insert(1:3)=r_target(1:3)+r_insert(1:3)
      ELSE

! find a new insertion point somewhere in the box
         DO
            DO i=1,3
               rand=next_random_number(rng_stream,error=error)
               r_insert(i)=rand*abc(i)
            ENDDO

! make sure it's not in the "in" region
            RIJ(1)=r_insert(1)-r_target(1)-abc(1)*&
               ANINT((r_insert(1)-r_target(1))/abc(1))
            RIJ(2)=r_insert(2)-r_target(2)-abc(2)*&
               ANINT((r_insert(2)-r_target(2))/abc(2))
            RIJ(3)=r_insert(3)-r_target(3)-abc(3)*&
               ANINT((r_insert(3)-r_target(3))/abc(3))

            dist=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2

            IF (dist .LT. rmin**2 .OR. dist .GT. rmax**2 ) THEN
               EXIT
            ENDIF

         ENDDO
      ENDIF

END SUBROUTINE generate_avbmc_insertion

END MODULE mc_coordinates

