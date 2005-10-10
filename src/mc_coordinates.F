!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_coordinates [1.0] *
!!
!!   NAME
!!     mc_coordinates
!!
!!   FUNCTION
!!     contains miscellaneous subroutines used in the Monte Carlo runs,mostly
!!     geared towards changes in coordinates
!! 
!!   AUTHOR
!!     MJM
!!
!******************************************************************************
MODULE mc_coordinates
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_pp_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             force_env_p_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: mc_averages_type,&
                                             mc_moves_type
  USE message_passing,                 ONLY: mp_bcast
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type,&
                                             mol_kind_new_list_p_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             atom_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             n_avogadro
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE :: nonbonded_energy, nonbonded_energy_spce

  PUBLIC :: get_water_data,mc_classical_energy,&
            generate_cbmc_swap_config,&
            get_center_of_mass,mc_coordinate_fold,&
            find_mc_test_molecule,&
            create_discrete_array,&
            check_for_overlap

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_coordinates"

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mc_coordinates.F/get_water_data *
!!
!!   NAME
!!     get_water_data
!!
!!   FUNCTION
!!     gets the running average of both bond lengths and the bond
!!     angle for every water molecule in the system
!!
!!   ARGUMENTS
!!     - force_env: the force_env you want to get information on
!!     - averages: the mc_averages structure to store the information
!!     - nnstep: the number of times this these averages have been
!!               done (including this one)
!!     - nchain: the number of molecules in this force_env
!!     - data: the unit to write the data to
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!*** ************************************************************************

SUBROUTINE get_water_data ( force_env,averages,nnstep,nchain,DATA  )


    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_averages_type), POINTER          :: averages
    INTEGER, INTENT(IN)                      :: nnstep, nchain, DATA

    CHARACTER(LEN=*), PARAMETER        :: routineN = "get_water_data"

    INTEGER                                  :: handle, iimol, imol
    REAL(KIND=dp)                            :: angle_loc, bond1_loc, &
                                                bond2_loc, tot_squared
    REAL(KIND=dp), DIMENSION(1:3)            :: vec_1, vec_2
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! initialize some stuff
      NULLIFY(oldsys,particles)
      bond1_loc=0.0E0_dp
      bond2_loc=0.0E0_dp
      angle_loc=0.0E0_dp

! get the particle coordinates
      CALL force_env_get(force_env,subsys=oldsys,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

! find the running average bond lengths
       DO iimol=1,nchain
          imol=3*iimol-2

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+1)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +1)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+1)%r(3))**2
          bond1_loc=bond1_loc*REAL(iimol-1,dp)/REAL(iimol,dp) &
                  +SQRT(tot_squared)/REAL(iimol,dp)

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+2)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +2)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+2)%r(3))**2
          bond2_loc=bond2_loc*REAL(iimol-1,dp)/REAL(iimol,dp) &
                  +SQRT(tot_squared)/REAL(iimol,dp)

! now find the running average bond angle
          vec_1(1:3)=particles%els(imol+1)%r(1:3)-&
              particles%els(imol)%r(1:3)
          vec_2(1:3)=particles%els(imol+2)%r(1:3)-&
              particles%els(imol)%r(1:3)

          angle_loc=angle_loc*REAL(iimol-1,dp)/REAL(iimol,dp) + &
                   1.0E0_dp/REAL(iimol,dp)*&
              ACOS(DOT_PRODUCT(vec_1,vec_2)/(bond1_loc*bond2_loc))&
              *180.0E0_dp/pi

       ENDDO
       averages%bond1=averages%bond1*REAL(nnstep-1,dp)/&
                  REAL(nnstep,dp)+bond1_loc/REAL(nnstep,dp)
       averages%bond2=averages%bond2*REAL(nnstep-1,dp)/&
                  REAL(nnstep,dp)+bond2_loc/REAL(nnstep,dp)
       averages%angle=averages%angle*REAL(nnstep-1,dp)/&
                  REAL(nnstep,dp)+angle_loc/REAL(nnstep,dp)

       WRITE(DATA,*) nnstep,bond1_loc*angstrom,bond2_loc*angstrom,&
                     angle_loc
 
! end the timing
       CALL timestop(0.0E0_dp,handle)


END SUBROUTINE get_water_data

!******************************************************************************

!******************************************************************************
!!****s* mc_coordinates.F/mc_classical_energy
!!
!!   NAME
!!     mc_classical_energy
!!
!!   FUNCTION
!!     calculates the classical energy of a system configuration using
!!     a potential force-matched to ab initio results for water
!!
!!   ARGUMENTS
!!     - force_env: the force environment containing the coordinates
!!     - nchain: the number of molecules in the box, probably could do
!!               a SIZE(particles%els) for this, but I haven't yet
!!     - rrcut: the intermolecular cutoff you want to use for the potential,
!!             in atomic units...be careful with this, as the potential is
!!             only zero at the cutoff for certain values
!!     - nstart: the number of the molecule to start calculating for
!!     - nend: the number of the molecule to end calculating for
!!     --- if nstart is not equal to 1 and nend is not equal to the total
!!         number of molecules, set them equal
!!     - classical_energy_new: the classical energy returned by the routine
!!     - cell_length: the length of the box...if none is specified,
!!                    it uses the cell found in the force_env
!!     - loptional: a flag that is used for debugging
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_classical_energy ( force_env, classical_energy_new,&
                                 nchain, nstart, nend ,rcut,loverlap,&
                                 cell_length)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(OUT)               :: classical_energy_new
    INTEGER, INTENT(IN)                      :: nstart, nend
    INTEGER, INTENT(IN)                      :: nchain
    REAL(KIND=dp), INTENT(IN)                :: rcut
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: cell_length
    LOGICAL, INTENT(OUT)                     :: loverlap

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_classical_energy"

    INTEGER                                  :: handle, imol, ipart, jmol, &
                                                jpart, jstart, ntij
    LOGICAL                                  :: lall, lcalc_charge
    REAL(KIND=dp) :: bend_energy, box_length, charge_energy, dist, kstretch, &
      LJ_energy, rijx, rijy, rijz, &
      stretch_energy, charge_energy_mol,LJ_energy_mol,&
      oo_min_cutoff,oh_min_cutoff,hh_min_cutoff
    REAL(KIND=dp), DIMENSION(1:3)            :: charges, kbend, RIJ,&
                                                abc
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_error_type)                      :: error
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:3, 1:nchain)   :: r
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! initialize some stuff
      loverlap=.FALSE.
      oo_min_cutoff=2.0E0_dp/angstrom
      oh_min_cutoff=1.0E0_dp/angstrom
      hh_min_cutoff=1.0E0_dp/angstrom
      classical_energy_new=0.0_dp

      NULLIFY(oldsys,particles)
      LJ_energy=0.0E0_dp
      charge_energy=0.0E0_dp
      bend_energy=0.0E0_dp
      stretch_energy=0.0E0_dp
      kstretch= 0.4015204E0_dp ! Hartrees/a.u.**2
      kbend(1)= 0.02384043E0_dp ! Hartrees/a.u.**2
      kbend(2)= -0.0506033E0_dp ! Hartrees/a.u.**3
      kbend(3)= 0.0287044E0_dp  ! Hartrees/a.u.**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for different classical potentials
! these are for the fake (IPBV) potential
!      charges(1:3) = (/ -0.7699784E0_dp,0.3849892E0_dp,0.3849892E0_dp /)
! these are for SPC/E
      charges(1:3) = (/ -0.8476E0_dp,0.4238E0_dp,0.4238E0_dp /)

! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

      IF( PRESENT(cell_length)) THEN
         box_length=cell_length
      ELSE
         box_length=abc(1)
      ENDIF

! put the coordinates into an easier matrix to manipulate
      DO imol=1,nchain
         r(1:3,1,imol)=particles%els(imol*3-2)%r(1:3)
         r(1:3,2,imol)=particles%els(imol*3-1)%r(1:3)
         r(1:3,3,imol)=particles%els(imol*3)%r(1:3)
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! commented out since we shouldn't be
! calculating potentials here anymore
! we need to calculate the energy of this configuration
! first, let's calculate intramolecular (bending and stretching)
!      DO imol=nstart,nend
!         
! stretching
!         OH_bond1=SQRT((r(1,1,imol)-r(1,2,imol))**2+&
!                 (r(2,1,imol)-r(2,2,imol))**2+&
!                   (r(3,1,imol)-r(3,2,imol))**2)
!         OH_bond2=SQRT((r(1,1,imol)-r(1,3,imol))**2+&
!                   (r(2,1,imol)-r(2,3,imol))**2+&
!                   (r(3,1,imol)-r(3,3,imol))**2)
!
!         stretch_energy=stretch_energy+0.5E0_dp*kstretch*((OH_bond1-&
!            1.81200E0_dp)**2+(OH_bond2-1.81200E0_dp)**2)
!
! bond bending
!         DO i=1,3
!            HH(i)=r(i,2,imol)-r(i,3,imol)
!         ENDDO
!         mag_HH=SQRT(DOT_PRODUCT(HH,HH))
!         bend_energy=bend_energy+0.5E0_dp*kbend(1)*(mag_HH-3.0699E0_dp)**2+&
!             0.33333333333333333E0_dp*kbend(2)*(mag_HH-3.0699E0_dp)**3+&
!             0.25E0_dp*kbend(3)*(mag_HH-3.0699E0_dp)**4
!      ENDDO

! now let's find the LJ energy between all the oxygens and
! the charge interactions
      lall=.FALSE.
      jstart=1
      IF( nstart == 1 .AND. nend == nchain) lall=.TRUE.
      DO imol=nstart,nend
         LJ_energy_mol=0.0E0_dp
         charge_energy_mol=0.0E0_dp
         IF (lall) jstart=imol+1
         DO jmol=jstart,nchain
            IF(imol == jmol) CYCLE

! find the minimum image distance between oxygens
            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*ANINT(&
               (r(1,1,imol)-r(1,1,jmol))/box_length)
            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*ANINT(&
               (r(2,1,imol)-r(2,1,jmol))/box_length)
            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*ANINT(&
               (r(3,1,imol)-r(3,1,jmol))/box_length)
            
            dist=SQRT(DOT_PRODUCT(RIJ,RIJ))

            IF (dist < rcut ) THEN
! check for overlap
               IF (dist < oo_min_cutoff) THEN
                  loverlap=.TRUE.
                  CALL timestop(0.0E0_dp,handle)
                  RETURN
               ELSE
                  lcalc_charge = .TRUE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for different classical potentials
!                  CALL nonbonded_energy(dist,unit_energy,1)
!                  CALL nonbonded_energy_spce(dist,unit_energy,1)
!                  LJ_energy_mol=LJ_energy_mol+unit_energy               
!                  charge_energy_mol=charge_energy_mol+charges(1)&
!                     *charges(1)/(dist) 
               ENDIF
            ELSE
               lcalc_charge=.FALSE.
            ENDIF

            DO ipart=1,3
               DO jpart=1,3        
                  IF (ipart == 1 .AND. jpart == 1) CYCLE

! find the minimum image distance
                  rijx=r(1,ipart,imol)-r(1,jpart,jmol)
                  rijy=r(2,ipart,imol)-r(2,jpart,jmol)
                  rijz=r(3,ipart,imol)-r(3,jpart,jmol)
                  RIJ(1)=rijx - box_length*AINT(rijx/box_length+&
                     SIGN(0.5E0_dp,rijx))
                  RIJ(2)=rijy - box_length*AINT(rijy/box_length+&
                     SIGN(0.5E0_dp,rijy))
                  RIJ(3)=rijz - box_length*AINT(rijz/box_length+&
                     SIGN(0.5E0_dp,rijz))

                  dist=SQRT(DOT_PRODUCT(RIJ,RIJ))

! LJ energy
!                  WRITE(6,*) 'eiruowe',dist,rcut
                  IF ( dist < rcut ) THEN
!                     ntij=(ipart-1)*10 + jpart
                     
                     IF (ipart .NE. 1 .AND. jpart .NE. 1) THEN
                        IF(dist < hh_min_cutoff) THEN
                           loverlap=.TRUE.
                           CALL timestop(0.0E0_dp,handle)
                           RETURN
                        ENDIF
                     ELSE
                        IF(dist < oh_min_cutoff) THEN
                           loverlap=.TRUE.
                           CALL timestop(0.0E0_dp,handle)
                           RETURN
                        ENDIF
                     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for different classical potentials
!                     CALL nonbonded_energy(dist,unit_energy,ntij)
!                     CALL nonbonded_energy_spce(dist,unit_energy,ntij)
!                     LJ_energy_mol=LJ_energy_mol+unit_energy
                  ENDIF

! coulombic
!                  IF ( lcalc_charge ) THEN
!                     charge_energy_mol=charge_energy_mol&
!                        +charges(ipart)*charges(jpart)/(dist)
!                  ENDIF

               ENDDO
            ENDDO
         ENDDO
!         LJ_energy=LJ_energy+LJ_energy_mol
!         charge_energy=charge_energy+charge_energy_mol
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for classical simulations
! sum up the total energy
!      classical_energy_new=bend_energy+&
!                stretch_energy+LJ_energy+charge_energy
!    classical_energy_new=LJ_energy+charge_energy
    classical_energy_new=0.0d0

! end the timing
      CALL timestop(0.0E0_dp,handle)


END SUBROUTINE mc_classical_energy

!******************************************************************************

!******************************************************************************
!!****s* mc_coordinates.F/check_for_overlap
!!
!!   NAME
!!     check_for_overlap
!!
!!   FUNCTION
!!     looks for overlaps (intermolecular distances less than rmin)
!!
!!   ARGUMENTS
!!     - force_env: the force environment containing the coordinates
!!     - nchain: the number of molecules in the box, probably could do
!!               a SIZE(particles%els) for this, but I haven't yet
!!     - rrcut: the intermolecular cutoff you want to use for the potential,
!!             in atomic units...be careful with this, as the potential is
!!             only zero at the cutoff for certain values
!!     - nstart: the number of the molecule to start calculating for
!!     - nend: the number of the molecule to end calculating for
!!     --- if nstart is not equal to 1 and nend is not equal to the total
!!         number of molecules, set them equal
!!     - classical_energy_new: the classical energy returned by the routine
!!     - cell_length: the length of the box...if none is specified,
!!                    it uses the cell found in the force_env
!!     - loptional: a flag that is used for debugging
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE check_for_overlap ( force_env, nchain, nstart, nend,loverlap,&
                                 cell_length)

    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(IN)                      :: nstart, nend
    INTEGER, INTENT(IN)                      :: nchain
    REAL(KIND=dp), DIMENSION(1:3),INTENT(IN), OPTIONAL  :: cell_length
    LOGICAL, INTENT(OUT)                     :: loverlap

    CHARACTER(LEN=*), PARAMETER :: routineN = "check_for_overlap"

    INTEGER                                  :: handle, imol, jmol, &
                                                junit,iunit, jstart,istat,&
                                                nunits
    LOGICAL                                  :: lall
    REAL(KIND=dp) :: dist, rmin
    REAL(KIND=dp), DIMENSION(1:3)            :: RIJ,&
                                                abc,box_length
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_error_type)                      :: error
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), &
      DIMENSION(:, :, :),ALLOCATABLE   :: r
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_kind_new_list_type), POINTER  &
                                             :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!------------------------------------------------------------------------------
! begin the timing of the subroutine
      
      CALL timeset(routineN,'I','',handle)

      NULLIFY(oldsys,particles,molecule_kinds_new,molecule_kind)

! initialize some stuff
      loverlap=.FALSE.
      rmin=3.57106767_dp ! 1 angstrom squared

! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            molecule_kinds_new=molecule_kinds_new,&
            particles=particles, error=error)
      molecule_kind => molecule_kinds_new%els(1)
      CALL get_molecule_kind(molecule_kind,natom=nunits)
      ALLOCATE(r(1:3,1:nunits,1:nchain),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass",3*nunits*nchain*dp_size)

      IF( PRESENT(cell_length)) THEN
         box_length(1:3)=cell_length(1:3)
      ELSE
         box_length(1:3)=abc(1:3)
      ENDIF

! put the coordinates into an easier matrix to manipulate
      DO imol=1,nchain
         DO iunit=1,nunits
            junit=imol*nunits-nunits+iunit
            r(1:3,iunit,imol)=particles%els(junit)%r(1:3)
         ENDDO
      ENDDO

! now let's find the LJ energy between all the oxygens and
! the charge interactions
      lall=.FALSE.
      jstart=1
      IF( nstart == 1 .AND. nend == nchain) lall=.TRUE.
      DO imol=nstart,nend
         IF (lall) jstart=imol+1
         DO jmol=jstart,nchain
            IF(imol == jmol) CYCLE

            DO iunit=1,nunits
               DO junit=1,nunits
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
                     IF (istat /= 0) CALL stop_memory(routineN,module_name,&
                     __LINE__,"r")
                     
                     CALL timestop(0.0E0_dp,handle)
                     RETURN
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r")

! end the timing
      CALL timestop(0.0E0_dp,handle)


END SUBROUTINE check_for_overlap

!******************************************************************************
!!****s* mc_coordinates.F/nonbonded_energy
!!
!!   NAME
!!     nonbonded_energy
!!
!!   FUNCTION
!!     calculates the nonbonded energy of two particles, based on the 
!!     seperation and atom types
!!
!!   ARGUMENTS
!!     - dist: the distance between the atoms (in a.u.)
!!     - ljfake: the energy generated by the subroutine (fake Lennard-Jones 
!!               potential), in Hartree
!!     - ntij: a number that describes the interaction
!!             1: oxygen-oxygen
!!             2,3,11,21: oxygen-hydrogen
!!             12,13,22,23L hydrogen-hydrogen
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE nonbonded_energy(dist,ljfake,ntij)

    REAL(KIND=dp), INTENT(IN)                :: dist
    REAL(KIND=dp), INTENT(OUT)               :: ljfake
    INTEGER, INTENT(IN)                      :: ntij

    INTEGER                                  :: index, j
    REAL(KIND=dp), DIMENSION(1:3)            :: b, m, pshift, r_core
    REAL(KIND=dp), DIMENSION(2:15, 1:3)      :: A

      index=0
      IF (ntij == 1 ) index=1
      IF (ntij == 2 ) index=2
      IF (ntij == 3 ) index=2
      IF (ntij == 11 ) index=2
      IF (ntij == 12 ) index=3
      IF (ntij == 13 ) index=3
      IF (ntij == 21 ) index=2
      IF (ntij == 22 ) index=3
      IF (ntij == 23 ) index=3
      IF (index == 0) WRITE(6,*) 'index =',index

      
! --- these pshift make the potential zero at 7.8 angstroms
!      pshift(1)=0.127827E0_dp ! Hartree
!      pshift(2)=2.22408E0_dp
!      pshift(3)=0.256530E0_dp
!!!!!!!!!!!!!!! for truncated IPBV potential
! --- these make it zero at 6.211 angstroms
      pshift(1)=0.128426E0_dp ! Hartree
      pshift(2)=2.223850E0_dp
      pshift(3)=0.256652E0_dp
      m(1) = 0.0E0_dp   ! Hartree/a.u.
      m(2) = -0.00437054E0_dp
      m(3) = 0.00377391E0_dp
! --- these b make the potential zero at 7.8 angstroms
!      b(1) = 1.0e23_dp ! Hartree
!      b(2) = 0.0313652E0_dp
!      b(3) = -0.0206204E0_dp
! --- these b make the potential zero at 6.211 angstroms
      b(1) = 1.0e23_dp ! Hartree
      b(2) = 0.0311350E0_dp
      b(3) = -0.0204983E0_dp
      r_core(1)=0.9E0_dp ! a.u.
      r_core(2)=2.95E0_dp
      r_core(3)=3.165E0_dp
      A(2,1) = 4.786380682394E0_dp   ! Hartree*a.u.^2
      A(3,1) = -1543.407053545E0_dp
      A(4,1) = 88783.31188529E0_dp
      A(5,1) = -2361200.155376E0_dp
      A(6,1) = 35940504.84679E0_dp
      A(7,1) = -339762743.6358E0_dp
      A(8,1) = 2043874926.466E0_dp
      A(9,1) = -7654856796.383E0_dp
      A(10,1) = 16195251405.65E0_dp
      A(11,1) = -13140392992.18E0_dp
      A(12,1) = -9285572894.245E0_dp
      A(13,1) = 8756947519.029E0_dp
      A(14,1) = 15793297761.67E0_dp
      A(15,1) = 12917180227.21E0_dp
      A(2,2) = -195.7716013277E0_dp
      A(3,2) = 15343.78613395E0_dp
      A(4,2) = -530864.4586516E0_dp
      A(5,2) = 10707934.39058E0_dp
      A(6,2) = -140099704.7890E0_dp
      A(7,2) = 1250943273.785E0_dp
      A(8,2) = -7795458330.676E0_dp
      A(9,2) = 33955897217.31E0_dp
      A(10,2) = -101135640744.0E0_dp
      A(11,2) = 193107995718.7E0_dp
      A(12,2) = -193440560940.0E0_dp
      A(13,2) = -4224406093.918E0_dp
      A(14,2) = 217192386506.5E0_dp
      A(15,2) = -157581228915.5E0_dp
      A(2,3) = -26.29456010782E0_dp
      A(3,3) = 2373.352548248E0_dp
      A(4,3) = -93880.43551360E0_dp
      A(5,3) = 2154624.884809E0_dp
      A(6,3) = -31965151.34955E0_dp
      A(7,3) = 322781785.3278E0_dp
      A(8,3) = -2271097368.668E0_dp
      A(9,3) = 11169163192.90E0_dp
      A(10,3) = -37684457778.47E0_dp
      A(11,3) = 82562104256.03E0_dp
      A(12,3) = -100510435213.4E0_dp
      A(13,3) = 24570342714.65E0_dp
      A(14,3) = 88766181532.94E0_dp
      A(15,3) = -79705131323.98E0_dp
      
      ljfake=0.0E0_dp
      IF (dist .GT. r_core(index)) THEN
! use a series expansion of 1/r
         DO j=2,15
            ljfake=ljfake+A(j,index)/(dist**(j-1)*REAL(j-1,dp))
        ENDDO
         ljfake=ljfake+pshift(index)
      ELSE
! use a linear potential
         ljfake=m(index)*dist+b(index)
      ENDIF

END SUBROUTINE nonbonded_energy

!******************************************************************************
!!****s* mc_coordinates.F/nonbonded_energy_spce
!!
!!   NAME
!!     nonbonded_energy_spce
!!
!!   FUNCTION
!!     calculates the nonbonded energy of two particles, based on the 
!!     seperation and atom types, using the SPC/E potential
!!
!!   ARGUMENTS
!!     - dist: the distance between the atoms (in a.u.)
!!     - ljfake: the energy generated by the subroutine (fake Lennard-Jones 
!!               potential), in Hartree
!!     - ntij: a number that describes the interaction
!!             1: oxygen-oxygen
!!             2,3,11,21: oxygen-hydrogen
!!             12,13,22,23L hydrogen-hydrogen
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE nonbonded_energy_spce(dist,ljfake,ntij)

    REAL(KIND=dp), INTENT(IN)                :: dist
    REAL(KIND=dp), INTENT(OUT)               :: ljfake
    INTEGER, INTENT(IN)                      :: ntij

    INTEGER                                  :: index
    REAL(KIND=dp)                            :: epsilon,sigma

      index=0
      IF (ntij == 1 ) index=1
      IF (ntij == 2 ) index=2
      IF (ntij == 3 ) index=2
      IF (ntij == 11 ) index=2
      IF (ntij == 12 ) index=3
      IF (ntij == 13 ) index=3
      IF (ntij == 21 ) index=2
      IF (ntij == 22 ) index=3
      IF (ntij == 23 ) index=3
      IF (index == 0) WRITE(6,*) 'index =',index

      ljfake=0.0E0_dp
! only the oxygens have a Lennard-Jones interaction
      IF( index .NE. 1) RETURN

! the sigma and epsilon are figured out using the standard L-B combining
! rules...units are in Bohr and Hartree,respectively
!      epsilon=78.1958E0_dp*k/4.3597482E0_dp
      epsilon=2.4763278E-4_dp
!      sigma=3.1655E0_dp/0.5291772083E0_dp
      sigma=5.981928077E0_dp

      ljfake=4.0E0_dp*epsilon*((sigma/dist)**12-(sigma/dist)**6)

END SUBROUTINE nonbonded_energy_spce

!******************************************************************************
!!****s* mc_coordinates.F/get_center_of_mass
!!
!!   NAME
!!     get_center_of_mass
!!
!!   FUNCTION
!!     calculates the center of mass of a given molecule
!!
!!   ARGUMENTS
!!     - natom: the number of atoms in the molecule
!!     - coordinates: the coordiantes of the atoms in the molecule
!!     - mass: the mass of the atoms in the molecule
!!     - center_of_mass: the coordinates of the center of mass
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE get_center_of_mass ( coordinates,natom,center_of_mass,&
            mass)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coordinates
    INTEGER, INTENT(IN)                      :: natom
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(OUT)                            :: center_of_mass
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: mass

    CHARACTER(LEN=*), PARAMETER :: routineN = "get_center_of_mass"

    INTEGER                                  :: handle, i, iatom
    REAL(KIND=dp)                            :: total_mass

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      total_mass=SUM(mass(:))
      center_of_mass(:)=0.0E0_dp

      DO iatom=1,natom
         DO i=1,3
            center_of_mass(i)=center_of_mass(i)+&
               mass(iatom)*coordinates(i,iatom)
         ENDDO
      ENDDO

      center_of_mass(1:3)=center_of_mass(1:3)/total_mass

! end the timing
      CALL timestop(0.0E0_dp,handle)


END SUBROUTINE get_center_of_mass

!******************************************************************************
!!****s* mc_coordinates.F/mc_coordinate_fold
!!
!!   NAME
!!     mc_coordinate_fold
!!
!!   FUNCTION
!!     folds all the coordinates into the center simulation box using
!!     a center of mass cutoff
!!
!!   ARGUMENTS
!!     - natoms: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the system
!!     - box_length: the length of the simulation box
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_coordinate_fold(coordinates,natoms_tot,mass,nunits,box_length)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: natoms_tot,nunits
    REAL(KIND=dp), DIMENSION(1:3),INTENT(IN) :: box_length
    REAL(KIND=dp), DIMENSION(:),INTENT(IN)   :: mass

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_coordinate_fold"

    INTEGER                                  :: handle, iatom,jatom,&
      iunit
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! loop over all molecules
      DO iatom=1,natoms_tot,nunits
         CALL get_center_of_mass(coordinates(:,iatom:iatom+nunits-1),&
            nunits,center_of_mass(:),mass(:))
         DO iunit=1,nunits
            jatom=iatom+iunit-1
            coordinates(1,jatom)=coordinates(1,jatom)-&
               box_length(1)*FLOOR(center_of_mass(1)/box_length(1))
            coordinates(2,jatom)=coordinates(2,jatom)-&
               box_length(2)*FLOOR(center_of_mass(2)/box_length(2))
            coordinates(3,jatom)=coordinates(3,jatom)-&
               box_length(3)*FLOOR(center_of_mass(3)/box_length(2))
         ENDDO

      ENDDO

! end the timing
      CALL timestop(0.0E0_dp,handle)


END SUBROUTINE mc_coordinate_fold


!******************************************************************************
!!****s* mc_coordinates.F/generate_cbmc_swap_config
!!
!!   NAME
!!     generate_cbmc_swap_config
!!
!!   FUNCTION
!!     takes the last molecule in a force environment and moves it around
!!     to different center of mass positions and orientations, selecting one
!!     based on the rosenbluth weight
!!
!!   ARGUMENTS
!!     - force_env: the force environment containing the coordinates
!!     - BETA: the value of 1/kT for this simulations, in a.u.
!!     - nswapmoves: the number of desired trial configurations
!!     - loverlap: the flag that determines if all of the configs have an
!!                 atomic overlap
!!     - rosenbluth_weight: the Rosenbluth weight for this set of configs
!!     - rem_molecule: if this is for the removal box, this is the molecule
!!                     that has been removed...not needed for insertions
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE generate_cbmc_swap_config( force_env, BETA, nswapmoves,&
                                      rosenbluth_weight, &
                                      loverlap, choosen_energy,&
                                      old_energy,cutoff,ionode,source,group,&
!                                      boltz_weights,rem_molecule)
                                      rem_molecule)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp),INTENT(IN)                 :: BETA
    REAL(KIND=dp),INTENT(IN)                 :: old_energy
    REAL(KIND=dp),INTENT(IN)                 :: cutoff
    REAL(KIND=dp),INTENT(OUT)                :: rosenbluth_weight
    INTEGER, INTENT(IN)                      :: nswapmoves,source,group
    INTEGER, INTENT(IN),OPTIONAL             :: rem_molecule
!    REAL(KIND=dp),DIMENSION(1:nswapmoves),INTENT(OUT)    &
!                                             :: boltz_weights
    LOGICAL, INTENT(OUT)                     :: loverlap
    LOGICAL, INTENT(IN)                      :: ionode
    REAL(KIND=dp),INTENT(OUT)                :: choosen_energy

    CHARACTER(LEN=*), PARAMETER &
      :: routineN = "generate_cbmc_swap_config"

    INTEGER                                  :: natoms_tot,handle,istat,&
       iatom,imove,nmols,i,choosen,start_atom,end_atom,imolecule,natom
    LOGICAL,DIMENSION(:),ALLOCATABLE         :: loverlap_array
    LOGICAL                                  :: all_overlaps 
    REAL(KIND=dp)                            :: rand,rdum,&
                                                bias_energy,&
                                                total_running_weight,&
                                                exponent
    REAL(KIND=dp), DIMENSION(:,:,:),&
      ALLOCATABLE                            :: r
    REAL(KIND=dp), DIMENSION(:,:),&
      ALLOCATABLE                            :: r_old
    REAL(KIND=dp), DIMENSION(:),&
      ALLOCATABLE                            :: boltz_weights,mass,&
      trial_energy
    REAL(KIND=dp), DIMENSION(1:3)            :: r_insert,center_of_mass,&
                                                diff,abc
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_error_type)                      :: error
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)


! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, molecule_kinds_new=molecule_kinds_new,&
            error=error)
      molecule_kind => molecule_kinds_new%els(1)
      CALL get_molecule_kind(molecule_kind,natom=natom,nmolecule=nmols)  
      natoms_tot=natom*nmols
      ALLOCATE(mass(1:natoms_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass",natoms_tot*dp_size)
      DO iatom=1,natoms_tot
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=mass(iatom))
      ENDDO

      IF(PRESENT(rem_molecule) ) THEN
         imolecule=rem_molecule ! this the the molecule we're removing
         start_atom=natom*imolecule-(natom-1)
         end_atom=natom*imolecule
      ELSE
         imolecule=nmols ! this is the new config, so it's the last molecule
         start_atom=natom*imolecule-(natom-1)
         end_atom=natom*imolecule
      ENDIF

      ALLOCATE(r_old(1:3,1:natoms_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_old",3*natoms_tot*dp_size)
      ALLOCATE(r(1:3,1:natoms_tot,1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r",3*natoms_tot*nswapmoves*dp_size)
      ALLOCATE(trial_energy(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "trial_energy",nswapmoves*dp_size)
      ALLOCATE(boltz_weights(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "boltz_weights",nswapmoves*dp_size)
      ALLOCATE(loverlap_array(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "loverlap_array",nswapmoves)

! initialize the arrays that need it
      loverlap_array(:)=.FALSE.
      loverlap=.FALSE.
      boltz_weights(:) = 0.0E0_dp
      trial_energy(:) = 0.0E0_dp
      r(:,:,:)=0.0E0_dp
      choosen_energy=0.0d0
      rosenbluth_weight=0.0d0

! save the positions of the molecules
      DO imove=1,nswapmoves
         DO iatom=1,natoms_tot
            r(1:3,iatom,imove)=particles%els(iatom)%r(1:3)
         ENDDO
      ENDDO

      DO iatom=1,natoms_tot
         r_old(1:3,iatom)=particles%els(iatom)%r(1:3)
      ENDDO

      IF(PRESENT(rem_molecule)) THEN
         CALL mc_classical_energy(force_env,rdum,&
            nmols,1,nmols,cutoff,loverlap_array(1))
         
         IF(loverlap_array(1)) THEN
            WRITE(6,*) start_atom,end_atom,natoms_tot
            DO iatom=1,natoms_tot
               WRITE(6,*) r(1:3,iatom,1)
            ENDDO
            CALL stop_program(routineN,module_name,__LINE__,&
               'CBMC swap move found an overlap in the old config')
         ENDIF
      ENDIF

      DO imove=1,nswapmoves
         
! drop into serial
         IF(ionode) THEN
! find a new insertion point
            DO i=1,3
               CALL random_number(rand)
               r_insert(i)=rand*abc(i)
            ENDDO

! find the center of mass of the insertion molecule
            CALL get_center_of_mass(r(:,start_atom:end_atom,imove),natom,&
            center_of_mass(:),mass(:))

! move the molecule to the insertion point
            DO i=1,3
               diff(i)=r_insert(i)-center_of_mass(i)
            ENDDO
            
            DO iatom=start_atom,end_atom
               r(1:3,iatom,imove)=r(1:3,iatom,imove)+diff(1:3)
            ENDDO
            
! rotate the molecule
            CALL rotate_molecule(r(:,start_atom:end_atom,imove),mass,natom)
            
            IF(imove == 1 .AND. PRESENT(rem_molecule)) THEN
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
         CALL mc_classical_energy(force_env,rdum,&
            nmols,imolecule,imolecule,cutoff,loverlap_array(imove))
         IF(loverlap_array(imove)) THEN
            boltz_weights(imove)=0.0d0
            CYCLE
         ENDIF
         CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.)
         CALL force_env_get(force_env,&
            potential_energy=bias_energy,error=error)


         trial_energy(imove)=(bias_energy-old_energy)

         exponent=-BETA*trial_energy(imove)
         IF(exponent .GT. 500.0d0) THEN
            boltz_weights(imove)=1.4E217_dp
         ELSEIF(exponent .LT. -500.0d0) THEN
            boltz_weights(imove)=0.0d0
         ELSE
            boltz_weights(imove)=EXP(exponent)
         ENDIF

      ENDDO

! now we need to pick a configuration based on the Rosenbluth weight,
! which is just the sum of the Boltzmann weights
      rosenbluth_weight=SUM(boltz_weights(:))
      all_overlaps=.TRUE.
      total_running_weight=0.0E0_dp
      choosen=0
      IF(ionode)THEN
         CALL random_number(rand)
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
         IF(PRESENT(rem_molecule)) THEN
            WRITE(6,*) boltz_weights(1:nswapmoves)
            WRITE(6,*) start_atom,end_atom,rem_molecule
            WRITE(6,*) loverlap_array(1:nswapmoves)
            DO iatom=1,natoms_tot
               WRITE(6,*) r(1:3,iatom,1)
            ENDDO
            CALL stop_program(routineN,module_name,__LINE__,&
            'CBMC swap move found all overlaps for the remove config')
         ENDIF

         DEALLOCATE(mass,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"mass")
         DEALLOCATE(r_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_old")
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r")
         DEALLOCATE(trial_energy,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"trial_energy")
         DEALLOCATE(boltz_weights,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"boltz_weights")
         DEALLOCATE(loverlap_array,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"loverlap_array")
         CALL timestop(0.0E0_dp,handle)
         RETURN
      ENDIF

! make sure a configuration was choosen
      IF(choosen == 0) CALL stop_program(routineN,module_name,&
         __LINE__,'CBMC swap move failed to select config')

! if this is an old configuration, we always choose the first one
      IF(PRESENT(rem_molecule)) choosen=1

! set the energy for the configuration
      choosen_energy=trial_energy(choosen)

! copy the coordinates to the force environment
      DO iatom=1,natoms_tot
         particles%els(iatom)%r(1:3)=r(1:3,iatom,choosen)
      ENDDO

      DEALLOCATE(mass,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"mass")
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r_old")
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r")
      DEALLOCATE(trial_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"trial_energy")
      DEALLOCATE(boltz_weights,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"boltz_weights")
      DEALLOCATE(loverlap_array,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"loverlap_array")

! end the timing
      CALL timestop(0.0E0_dp,handle)


END SUBROUTINE generate_cbmc_swap_config

!******************************************************************************
!!****s* mc_coordinates.F/rotate_water_molecule
!!
!!   NAME
!!     rotate_water_molecule
!!
!!   FUNCTION
!!     rotates a water molecule randomly around the center of mass, 
!!     sequentially in x, y, and z directions
!!
!!   ARGUMENTS
!!     - r: the coordinates of the water molecule to rotate
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE rotate_molecule(r,mass,natoms)

    INTEGER,INTENT(IN)                            :: natoms
    REAL(KIND=dp),DIMENSION(1:3,1:natoms),INTENT(INOUT)  :: r
    REAL(KIND=dp),DIMENSION(1:natoms),INTENT(IN):: mass 

    CHARACTER(LEN=*), PARAMETER &
      :: routineN = "rotate_molecule"

    INTEGER                                  :: iunit,handle
    REAL(KIND=dp)                            :: cosdg,sindg,ry,rz,rynew,&
      rznew,rxnew,rx,dgamma,rand
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! find the center of mass of the molecule
      CALL get_center_of_mass(r(:,:),natoms,center_of_mass(:),mass(:))

! call a random number to figure out how far we're moving
      CALL RANDOM_NUMBER(rand)
      dgamma=pi*(rand-0.5E0_dp)*2.0E0_dp

! *** set up the rotation marix ***
 
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
      CALL timestop(0.0E0_dp,handle)

END SUBROUTINE rotate_molecule

!******************************************************************************
!!****s* mc_coordinates.F/find_mc_test_molecule
!!
!!   NAME
!!     find_mc_test_molecule
!!
!!   FUNCTION
!!     returns the number of the first atom in a molecule that we
!!     want to perform a Monte Carlo move on
!!
!!   ARGUMENTS
!!     - natoms: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the system
!!     - box_length: the length of the simulation box
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE find_mc_test_molecule(force_env,mol_number,box_number,nchain,box,error)

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env
    INTEGER, INTENT(OUT)                     :: mol_number
    INTEGER, INTENT(OUT)                     :: box_number
    INTEGER, INTENT(IN),DIMENSION(:)         :: nchain
    INTEGER, INTENT(IN),OPTIONAL             :: box
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "find_mc_test_molecule",&
        routineP = module_name//':'//routineN

    TYPE(mol_kind_new_list_p_type),DIMENSION(:),POINTER :: mol_kinds
    TYPE(cp_subsystem_pp_type),DIMENSION(:),POINTER :: subsys
    TYPE(molecule_type),POINTER              :: molecule
    TYPE(mol_new_list_type),POINTER          :: molecules_new

    INTEGER                                  :: handle,ibox,nboxes,stat,&
      target_molecule,nmols,imol,start_box,end_box
    REAL(KIND=dp)                            :: rand
    LOGICAL                                  :: failure

!------------------------------------------------------------------------------
! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

      NULLIFY(mol_kinds,subsys)

      nboxes=SIZE(force_env(:))

      IF(PRESENT(box))THEN
         start_box=box
         end_box=box
      ELSE
         start_box=1
         end_box=nboxes
      ENDIF

      mol_number=0
      box_number=0

      ALLOCATE(subsys(1:nboxes),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(mol_kinds(1:nboxes), stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

! find out how many molecules we have total
      DO ibox=start_box,end_box
         CALL force_env_get(force_env(ibox)%force_env,&
            subsys=subsys(ibox)%subsys,error=error)
         CALL cp_subsys_get(subsys(ibox)%subsys(1)%subsys, &
            molecule_kinds_new=mol_kinds(ibox)%list,error=error)
      ENDDO

! grab a random molecule
      CALL random_number(rand)
      target_molecule=CEILING(rand*REAL(SUM(nchain(start_box:end_box))))

! find the starting atom and box number
      box_loop:DO ibox=start_box,end_box
         IF(target_molecule .GT. SUM(nchain(start_box:end_box))) CYCLE
         box_number=ibox
         IF(ibox .NE. start_box) THEN
            target_molecule=target_molecule-&
               SUM(nchain(start_box:(ibox-1)))
         ENDIF
         CALL cp_subsys_get(subsys(ibox)%subsys(1)%subsys, &
            molecules_new=molecules_new,error=error) 
         nmols=SIZE(molecules_new%els(:))
         DO imol=1,nmols
            IF(target_molecule .EQ. imol) THEN
               molecule => molecules_new%els(imol)
               CALL get_molecule(molecule,first_atom=mol_number)
               EXIT box_loop
            ENDIF
         ENDDO
      ENDDO box_loop

! deallocate some stuff
      DEALLOCATE(mol_kinds)
      DEALLOCATE(subsys)
!      DEALLOCATE(nmolecules)

! make sure things are good
      CPPostcondition(mol_number/=0,cp_failure_level,routineP,error,failure)
      CPPostcondition(box_number/=0,cp_failure_level,routineP,error,failure)
      
! end the timing
      CALL timestop(0.0E0_dp,handle)


END SUBROUTINE find_mc_test_molecule

!******************************************************************************
!!****s* mc_coordinates.F/find_mc_test_molecule
!!
!!   NAME
!!     find_mc_test_molecule
!!
!!   FUNCTION
!!     returns the number of the first atom in a molecule that we
!!     want to perform a Monte Carlo move on
!!
!!   ARGUMENTS
!!     - natoms: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the system
!!     - box_length: the length of the simulation box
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE create_discrete_array(cell,discrete_array,step_size)

! 1 is for increase, 2 is for decrease
! 1 is for "yes, we can do the move", 0 is for no

      REAL(dp),INTENT(IN)   :: step_size
      REAL(dp),DIMENSION(1:3),INTENT(IN)   :: cell
      INTEGER,DIMENSION(1:3,1:2),INTENT(OUT)   :: discrete_array

      REAL(dp) :: length1,length2,high_value,low_value
      INTEGER :: iside

      discrete_array(:,:)=0

      length1=ABS(cell(1)-cell(2))
      length2=ABS(cell(2)-cell(3))

! now let's figure out all the different cases
      IF(length1 .LT. 0.01d0*step_size .AND. &
         length2 .LT. 0.01d0*step_size) THEN
! all sides are equal, so we can move up or down
         discrete_array(1:3,1)=1
         discrete_array(1:3,2)=1
      ELSE

! find the low value and the high value
         high_value=-1.0d0
         low_value=cell(1)*cell(2)*cell(3)
         DO iside=1,3
            IF(cell(iside) .LT. low_value) low_value=cell(iside)
            IF(cell(iside) .GT. high_value) high_value=cell(iside)
         ENDDO
         DO iside=1,3
! now we see if the value is a high value or a low value...it can only be
! one of the two
            IF(ABS(cell(iside)-low_value) .LT. 0.01d0*step_size) THEN
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


END MODULE mc_coordinates

!******************************************************************************

