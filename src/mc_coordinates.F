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
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dbl,&
                                             dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: mc_averages_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             n_avogadro
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE :: nonbonded_energy
  PUBLIC :: get_water_data,mc_classical_energy,&
            check_for_overlap,generate_cbmc_swap_config,&
            get_center_of_mass,mc_coordinate_fold

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
      bond1_loc=0.0d0
      bond2_loc=0.0d0
      angle_loc=0.0d0

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
          bond1_loc=bond1_loc*REAL(iimol-1,dbl)/REAL(iimol,dbl) &
                  +DSQRT(tot_squared)/REAL(iimol,dbl)

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+2)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +2)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+2)%r(3))**2
          bond2_loc=bond2_loc*REAL(iimol-1,dbl)/REAL(iimol,dbl) &
                  +DSQRT(tot_squared)/REAL(iimol,dbl)

! now find the running average bond angle
          vec_1(1:3)=particles%els(imol+1)%r(1:3)-&
              particles%els(imol)%r(1:3)
          vec_2(1:3)=particles%els(imol+2)%r(1:3)-&
              particles%els(imol)%r(1:3)

          angle_loc=angle_loc*REAL(iimol-1,dbl)/REAL(iimol,dbl) + &
                   1.0d0/REAL(iimol,dbl)*&
              DACOS(DOT_PRODUCT(vec_1,vec_2)/(bond1_loc*bond2_loc))&
              *180.0d0/pi

       ENDDO
       averages%bond1=averages%bond1*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+bond1_loc/REAL(nnstep,dbl)
       averages%bond2=averages%bond2*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+bond2_loc/REAL(nnstep,dbl)
       averages%angle=averages%angle*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+angle_loc/REAL(nnstep,dbl)

       WRITE(DATA,*) nnstep,bond1_loc*angstrom,bond2_loc*angstrom,&
                     angle_loc
 
! end the timing
       CALL timestop(0.0d0,handle)


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
!!     - nend: the number of the moelcule to end calculating for
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
                                 nchain, nstart, nend ,rrcut,cell_length,&
                                 loptional)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(OUT)               :: classical_energy_new
    INTEGER, INTENT(IN)                      :: nstart, nend
    INTEGER, INTENT(IN)                      :: nchain
    REAL(KIND=dp), INTENT(IN)                :: rrcut
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: cell_length
    LOGICAL, INTENT(IN), OPTIONAL            :: loptional

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_classical_energy"

    INTEGER                                  :: handle, i, imol, ipart, jmol, &
                                                jpart, jstart, ntij
    LOGICAL                                  :: lall, lcalc_charge
    REAL(KIND=dp) :: bend_energy, box_length, charge_energy, dist, kstretch, &
      LJ_energy, mag_HH, OH_bond1, OH_bond2, rcut, rijx, rijy, rijz, &
      stretch_energy, unit_energy, charge_energy_mol,LJ_energy_mol
    REAL(KIND=dp), DIMENSION(1:3)            :: charges, HH, kbend, RIJ
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
      NULLIFY(oldsys,particles)
      LJ_energy=0.0d0
      charge_energy=0.0d0
      bend_energy=0.0d0
      stretch_energy=0.0d0
      kstretch= 0.4015204d0 ! Hartrees/a.u.**2
      kbend(1)= 0.02384043d0 ! Hartrees/a.u.**2
      kbend(2)= -0.0506033d0 ! Hartrees/a.u.**3
      kbend(3)= 0.0287044d0  ! Hartrees/a.u.**4
      charges(1:3) = (/ -0.7699784d0,0.3849892d0,0.3849892d0 /)
!      charges(:)=charges(:)*1.60d-19 ! Coulumbs

! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

      IF( PRESENT(cell_length)) THEN
         box_length=cell_length
      ELSE
         box_length=cell%hmat(1,1)
      ENDIF

!!!!!!!!!!
      rcut=rrcut
      rcut=7.8d0/angstrom
!!!!!!!!!!

! put the coordinates into an easier matrix to manipulate
      DO imol=1,nchain
         r(1:3,1,imol)=particles%els(imol*3-2)%r(1:3)
         r(1:3,2,imol)=particles%els(imol*3-1)%r(1:3)
         r(1:3,3,imol)=particles%els(imol*3)%r(1:3)
      ENDDO

! we need to calculate the energy of this configuration
! first, let's calculate intramolecular (bending and stretching)
      DO imol=nstart,nend
         
! stretching
         OH_bond1=DSQRT((r(1,1,imol)-r(1,2,imol))**2+&
                 (r(2,1,imol)-r(2,2,imol))**2+&
                   (r(3,1,imol)-r(3,2,imol))**2)
         OH_bond2=DSQRT((r(1,1,imol)-r(1,3,imol))**2+&
                   (r(2,1,imol)-r(2,3,imol))**2+&
                   (r(3,1,imol)-r(3,3,imol))**2)

         stretch_energy=stretch_energy+0.5d0*kstretch*((OH_bond1-&
            1.81200d0)**2+(OH_bond2-1.81200d0)**2)

! bond bending
         DO i=1,3
            HH(i)=r(i,2,imol)-r(i,3,imol)
         ENDDO
         mag_HH=DSQRT(DOT_PRODUCT(HH,HH))
         bend_energy=bend_energy+0.5d0*kbend(1)*(mag_HH-3.0699d0)**2+&
             0.33333333333333333d0*kbend(2)*(mag_HH-3.0699d0)**3+&
             0.25d0*kbend(3)*(mag_HH-3.0699d0)**4
      ENDDO

! now let's find the LJ energy between all the oxygens and
! the charge interactions
      lall=.FALSE.
      jstart=1
      IF( nstart == 1 .AND. nend == nchain) lall=.TRUE.
      DO imol=nstart,nend
         LJ_energy_mol=0.0d0
         charge_energy_mol=0.0d0
!!!!!!!!!!!!!!
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
            
            dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))

            IF (dist < rcut ) THEN
               lcalc_charge = .TRUE.
               CALL nonbonded_energy(dist,unit_energy,1)
               LJ_energy_mol=LJ_energy_mol+unit_energy               
               charge_energy_mol=charge_energy_mol+charges(1)&
                   *charges(1)/(dist) 
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
                  RIJ(1)=rijx - box_length*DINT(rijx/box_length+&
                     DSIGN(0.5D0,rijx))
                  RIJ(2)=rijy - box_length*DINT(rijy/box_length+&
                     DSIGN(0.5D0,rijy))
                  RIJ(3)=rijz - box_length*DINT(rijz/box_length+&
                     DSIGN(0.5D0,rijz))

                  dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))

! LJ energy
                  IF ( dist < rcut ) THEN
                     ntij=(ipart-1)*10 + jpart
                     CALL nonbonded_energy(dist,unit_energy,ntij)

                     LJ_energy_mol=LJ_energy_mol+unit_energy
                  ENDIF

! coulombic
                  IF ( lcalc_charge ) THEN
                     charge_energy_mol=charge_energy_mol&
                        +charges(ipart)*charges(jpart)/(dist)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
         LJ_energy=LJ_energy+LJ_energy_mol
         charge_energy=charge_energy+charge_energy_mol
      ENDDO

! sum up the total energy
!      classical_energy_new=bend_energy+&
!                stretch_energy+LJ_energy+charge_energy
    classical_energy_new=LJ_energy+charge_energy
!!!!!!!!!!!!!! 

! scale the total energy to get classical vs. Quickstep energy diffs closer
!      classical_energy_new=classical_energy_new*1.5d0

! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE mc_classical_energy

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
      pshift(1)=0.127827d0 ! Hartree
      pshift(2)=2.22408d0
      pshift(3)=0.256530d0
!!!!!!!!!!!!!!!
! --- these make it zero at 6.211 angstroms
!      pshift(1)=0.128426d0 ! Hartree
!      pshift(2)=2.223850d0
!      pshift(3)=0.256652d0
      m(1) = 0.0d0   ! Hartree/a.u.
      m(2) = -0.00437054d0
      m(3) = 0.00377391d0
! --- these b make the potential zero at 7.8 angstroms
      b(1) = 1.0d23 ! Hartree
      b(2) = 0.0313652d0
      b(3) = -0.0206204d0
! --- these b make the potential zero at 6.211 angstroms
!      b(1) = 1.0d23 ! Hartree
!      b(2) = 0.0311350d0
!      b(3) = -0.0204983d0
      r_core(1)=0.9d0 ! a.u.
      r_core(2)=2.95d0
      r_core(3)=3.165d0
      A(2,1) = 4.786380682394d0   ! Hartree*a.u.^2
      A(2,2) = -195.7716013277d0
      A(2,3) = -26.29456010782d0
      A(3,1) = -1543.407053545d0
      A(3,2) = 15343.78613395d0
      A(3,3) = 2373.352548248d0
      A(4,1) = 88783.31188529d0
      A(4,2) = -530864.4586516d0
      A(4,3) = -93880.43551360d0
      A(5,1) = -2361200.155376d0
      A(5,2) = 10707934.39058d0
      A(5,3) = 2154624.884809d0
      A(6,1) = 35940504.84679d0
      A(6,2) = -140099704.7890d0
      A(6,3) = -31965151.34955d0
      A(7,1) = -339762743.6358d0
      A(7,2) = 1250943273.785d0
      A(7,3) = 322781785.3278d0
      A(8,1) = 2043874926.466d0
      A(8,2) = -7795458330.676d0
      A(8,3) = -2271097368.668d0
      A(9,1) = -7654856796.383d0
      A(9,2) = 33955897217.31d0
      A(9,3) = 11169163192.90d0
      A(10,1) = 16195251405.65d0
      A(10,2) = -101135640744.0d0
      A(10,3) = -37684457778.47d0
      A(11,1) = -13140392992.18d0
      A(11,2) = 193107995718.7d0
      A(11,3) = 82562104256.03d0
      A(12,1) = -9285572894.245d0
      A(12,2) = -193440560940.0d0
      A(12,3) = -100510435213.4d0
      A(13,1) = 8756947519.029d0
      A(13,2) = -4224406093.918d0
      A(13,3) = 24570342714.65d0
      A(14,1) = 15793297761.67d0
      A(14,2) = 217192386506.5d0
      A(14,3) = 88766181532.94d0
      A(15,1) = 12917180227.21d0
      A(15,2) = -157581228915.5d0
      A(15,3) = -79705131323.98d0
      
      ljfake=0.0d0
      IF (dist .GT. r_core(index)) THEN
! use a series expansion of 1/r
         DO j=2,15
            ljfake=ljfake+A(j,index)/(dist**(j-1)*REAL(j-1,dbl))
        ENDDO
         ljfake=ljfake+pshift(index)
      ELSE
! use a linear potential
         ljfake=m(index)*dist+b(index)
      ENDIF

END SUBROUTINE nonbonded_energy

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
      center_of_mass(:)=0.0D0

      DO iatom=1,natom
         DO i=1,3
            center_of_mass(i)=center_of_mass(i)+&
               mass(iatom)*coordinates(i,iatom)
         ENDDO
      ENDDO

      center_of_mass(1:3)=center_of_mass(1:3)/total_mass

! end the timing
      CALL timestop(0.0d0,handle)


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
SUBROUTINE mc_coordinate_fold(coordinates,natoms,box_length)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: natoms
    REAL(KIND=dp), INTENT(IN)                :: box_length

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_coordinate_fold"

    INTEGER                                  :: handle, i, iatom
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, mass

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      mass(1:3)=(/ 15.999d0,1.0074d0,1.0074d0 /)

! loop over all molecules
      DO iatom=1,natoms,3
         CALL get_center_of_mass(coordinates(:,iatom:iatom+2),3,&
            center_of_mass(:),mass(:))
         DO i=1,3
            coordinates(i,iatom)=coordinates(i,iatom)-&
               box_length*FLOOR(center_of_mass(i)/box_length)
            coordinates(i,iatom+1)=coordinates(i,iatom+1)-&
               box_length*FLOOR(center_of_mass(i)/box_length)
            coordinates(i,iatom+2)=coordinates(i,iatom+2)-&
               box_length*FLOOR(center_of_mass(i)/box_length)
         ENDDO

      ENDDO



! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE mc_coordinate_fold


!******************************************************************************
!!****s* mc_coordinates.F/check_for_overlap
!!
!!   NAME
!!     check_for_overlap
!!
!!   SYNOPSIS
!!     Subroutine check_for_overlap(coordinates,natoms,box_length,loverlap)
!!           Real(Kind=dp), Intent(IN)             :: box_length
!!           Real(Kind=dp), Dimension(:,:),Intent(IN)
!!                                                 :: coordinates
!!           Integer, Intent(IN)                   :: natoms
!!           Logical, Intent(OUT)                  :: loverlap
!!     End Subroutine check_for_overlap
!!
!!   FUNCTION
!!     calculates all interatomic distances and looks to see if any are
!!     below a threshold value (based on the radial distribution functions)..
!!     these are considered to be overlaps because their Boltzmann weight 
!!     will be very small
!!
!!   ARGUMENTS
!!     - natoms: the total number of atoms in the system
!!     - coordinates: the coordiantes of the atoms in the system
!!     - box_length: the length of the simulation box in a.u.
!!     - loverlap: the flag that signals an overlap for this config
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE check_for_overlap(coordinates,natoms,box_length,loverlap)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coordinates
    INTEGER, INTENT(IN)                      :: natoms
    REAL(KIND=dp), INTENT(IN)                :: box_length
    LOGICAL, INTENT(OUT)                     :: loverlap

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "check_for_overlap"

    INTEGER                                  :: handle, i, iatom, imol, j, &
                                                jmol, nmol
    REAL(KIND=dp)                            :: dist
    REAL(KIND=dp), DIMENSION(1:3)            :: RIJ
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:3, natoms/3)          :: r

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      DO iatom=1,natoms,3
         r(1:3,1,(iatom+2)/3)=coordinates(1:3,iatom)
         r(1:3,2,(iatom+2)/3)=coordinates(1:3,iatom+1)
         r(1:3,3,(iatom+2)/3)=coordinates(1:3,iatom+2)
      ENDDO

      nmol=natoms/3

      DO imol=1,nmol-1
         DO jmol=imol+1,nmol
            
! check the O-O distance
            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*ANINT(&
               (r(1,1,imol)-r(1,1,jmol))/box_length)
            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*ANINT(&
               (r(2,1,imol)-r(2,1,jmol))/box_length)
            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*ANINT(&
               (r(3,1,imol)-r(3,1,jmol))/box_length)
            
            dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))
            
            IF(dist .LT. 2.2D0/angstrom) THEN
               loverlap=.TRUE.
               CALL timestop(0.0d0,handle)
               RETURN
            ELSEIF(dist .LT. 2.6D0/angstrom) THEN
!     check the H-H distance
               DO i=2,3
                  DO j=2,3
                        RIJ(1)=r(1,i,imol)-r(1,j,jmol)-box_length&
                           *ANINT((r(1,i,imol)-r(1,j,jmol))/box_length)
                        RIJ(2)=r(2,i,imol)-r(2,j,jmol)-box_length&
                           *ANINT((r(2,i,imol)-r(2,j,jmol))/box_length)
                        RIJ(3)=r(3,i,imol)-r(3,j,jmol)-box_length&
                           *ANINT((r(3,i,imol)-r(3,j,jmol))/box_length)
                        dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        IF (dist .LT. 1.3D0/angstrom) THEN
                           loverlap=.TRUE.
                           CALL timestop(0.0d0,handle)
                           RETURN
                        ENDIF
                  ENDDO
               ENDDO
            ENDIF

         ENDDO
      ENDDO

! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE check_for_overlap

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
                                      loverlap,rosenbluth_weight, &
                                      rem_molecule)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp),INTENT(IN)                 :: BETA
    REAL(KIND=dp),INTENT(OUT)                :: rosenbluth_weight
    INTEGER, INTENT(IN)                      :: nswapmoves
    INTEGER, INTENT(IN),OPTIONAL             :: rem_molecule
    LOGICAL, INTENT(OUT)                     :: loverlap

    CHARACTER(LEN=*), PARAMETER &
      :: routineN = "generate_cbmc_swap_config"

    INTEGER                                  :: natoms,handle,istat,&
       iatom,imove,nmols,i,choosen,start_atom,end_atom,imolecule
    LOGICAL,DIMENSION(:),ALLOCATABLE         :: loverlap_array
    LOGICAL                                  :: all_overlaps 
    REAL(KIND=dp)                            :: cell_length,rand,&
                                                classical_energy,&
                                                total_running_weight
    REAL(KIND=dp), DIMENSION(:,:,:),&
      ALLOCATABLE                            :: r
    REAL(KIND=dp), DIMENSION(:),&
      ALLOCATABLE                            :: boltz_weight
    REAL(KIND=dp), DIMENSION(1:3)            :: r_insert,center_of_mass,&
                                                diff
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_error_type)                      :: error
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)


! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

      natoms=SIZE(particles%els(:))
      nmols=natoms/3
      cell_length=cell%hmat(1,1)
      IF(PRESENT(rem_molecule) ) THEN
         imolecule=rem_molecule ! this the the molecule we're removing
         start_atom=3*imolecule-2
         end_atom=3*imolecule
      ELSE
         imolecule=nmols ! this is the new config, so it's the last molecule
         start_atom=3*imolecule-2
         end_atom=3*imolecule
      ENDIF

      ALLOCATE(r(1:3,1:natoms,1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r",3*natoms*nswapmoves*dp_size)
      ALLOCATE(boltz_weight(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "boltz_weight",nswapmoves*dp_size)
      ALLOCATE(loverlap_array(1:nswapmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "loverlap_array",nswapmoves)

! initialize the arrays that need it
      loverlap_array(:)=.FALSE.
      loverlap=.FALSE.
      boltz_weight(:) = 0.0d0

! need to calculate the classical energy for just the molecule inserted, which
! is always the last molecule
      CALL mc_classical_energy(force_env,classical_energy,nmols,&
         imolecule,imolecule,6.211D0/angstrom)

! calculate the Bolztmann weight for the first config
      boltz_weight(1)=DEXP(-BETA*classical_energy)

! save the positions of the molecules
      DO iatom=1,natoms
         r(1:3,iatom,1)=particles%els(iatom)%r(1:3)
      ENDDO

      DO imove=2,nswapmoves
         
! initialize the array with the first configuration
         DO iatom=1,natoms
            r(1:3,iatom,imove)=r(1:3,iatom,1)
         ENDDO

! find a new insertion point
         DO i=1,3
            CALL random_number(rand)
            r_insert(i)=rand*cell_length
         ENDDO

! find the center of mass of the insertion molecule
         CALL get_center_of_mass(r(:,natoms-2:natoms,1),3,&
            center_of_mass(:),(/15.999d0,1.008d0,1.008d0/))

! move the molecule to the insertion point
         DO i=1,3
            diff(i)=r_insert(i)-center_of_mass(i)
         ENDDO
!         WRITE(6,*) 'mcnv,nmx',r_insert(1:3)*0.529177d0
!         WRITE(6,*) ',nmx',center_of_mass(1:3)*0.529177d0

         DO iatom=natoms-2,natoms
            r(1:3,iatom,imove)=r(1:3,iatom,imove)+diff(1:3)
         ENDDO

! rotate the molecule
         CALL rotate_water_molecule(r(:,natoms-2:natoms,imove))

! check for overlap
         CALL check_for_overlap(r(:,:,imove),natoms,cell_length,&
            loverlap_array(imove))
         IF(loverlap_array(imove)) CYCLE

! calculate the energy and boltzman weight of the config
         DO iatom=start_atom,end_atom
            particles%els(iatom)%r(1:3)=r(1:3,iatom,imove)
         ENDDO
         CALL mc_classical_energy(force_env,classical_energy,&
            nmols,imolecule,imolecule,6.211D0/angstrom,loptional=.TRUE.)
         boltz_weight(imove)=DEXP(-BETA*classical_energy)
     
      ENDDO

! now we need to pick a configuration based on the Rosenbluth weight,
! which is just the sum of the Boltzmann weights
      rosenbluth_weight=SUM(boltz_weight(:))
      all_overlaps=.TRUE.
      total_running_weight=0.0d0
      choosen=0
      DO imove=1,nswapmoves
         IF ( loverlap_array(imove) ) CYCLE
         all_overlaps=.FALSE.
         total_running_weight=total_running_weight+boltz_weight(imove)
         CALL random_number(rand)
         IF ( total_running_weight .GE. rand*rosenbluth_weight) THEN
            choosen=imove
            EXIT
         ENDIF
      ENDDO
      

      IF(all_overlaps) THEN
         loverlap = .TRUE.
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r")
         DEALLOCATE(boltz_weight,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"boltz_weight")
         DEALLOCATE(loverlap_array,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"loverlap_array")
         CALL timestop(0.0d0,handle)
         RETURN
      ENDIF

! make sure a configuration was choosen
      IF(choosen == 0) CALL stop_program(routineN,module_name,&
         __LINE__,'CBMC swap move failed to select config')

! if this is an old configuration, we always choose the first one
      IF(PRESENT(rem_molecule)) choosen=1

! copy the coordinates to the force environment
      DO iatom=1,natoms
         particles%els(iatom)%r(1:3)=r(1:3,iatom,choosen)
      ENDDO

      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r")
      DEALLOCATE(boltz_weight,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"boltz_weight")
      DEALLOCATE(loverlap_array,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"loverlap_array")

! end the timing
      CALL timestop(0.0d0,handle)


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
SUBROUTINE rotate_water_molecule(r)

    REAL(KIND=dp),DIMENSION(1:3,1:3),INTENT(INOUT)  :: r

    CHARACTER(LEN=*), PARAMETER &
      :: routineN = "rotate_water_molecule"
    REAL(KIND=dp),DIMENSION(1:3),PARAMETER   &
      :: mass = (/ 15.999d0,1.008d0,1.008d0 /)

    INTEGER                                  :: iunit,handle
    REAL(KIND=dp)                            :: cosdg,sindg,ry,rz,rynew,&
      rznew,rxnew,rx,dgamma,rand
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! find the center of mass of the molecule
      CALL get_center_of_mass(r(:,:),3,center_of_mass(:),mass(:))

! call a random number to figure out how far we're moving
      CALL RANDOM_NUMBER(rand)
      dgamma=pi*(rand-0.5d0)*2.0d0

! *** set up the rotation marix ***
 
      cosdg = DCOS( dgamma )
      sindg = DSIN( dgamma )
 
! ***    ROTATE UNITS OF I AROUND X-AXIS ***
 
      DO  iunit = 1,3
         ry = r(2,iunit) - center_of_mass(2)
         rz = r(3,iunit) - center_of_mass(3)
         rynew = cosdg * ry + sindg * rz
         rznew = cosdg * rz - sindg * ry
         
         r(2,iunit) = rynew + center_of_mass(2)
         r(3,iunit) = rznew + center_of_mass(3)
         
      ENDDO
 
! ***    ROTATE UNITS OF I AROUND y-AXIS ***
 
      DO  iunit = 1,3
         rx = r(1,iunit) - center_of_mass(1)
         rz = r(3,iunit) - center_of_mass(3)
         rxnew = cosdg * rx + sindg * rz
         rznew = cosdg * rz - sindg * rx
         
         r(1,iunit) = rxnew + center_of_mass(1)
         r(3,iunit) = rznew + center_of_mass(3)
         
      ENDDO
      
! ***    ROTATE UNITS OF I AROUND z-AXIS ***
 
      DO  iunit = 1,3
         rx = r(1,iunit) - center_of_mass(1)
         ry = r(2,iunit) - center_of_mass(2)
         rxnew = cosdg * rx + sindg * ry
         rynew = cosdg * ry - sindg * rx
         
         r(1,iunit) = rxnew + center_of_mass(1)
         r(2,iunit) = rynew + center_of_mass(2)
         
      ENDDO
      
! end the timing
      CALL timestop(0.0d0,handle)

END SUBROUTINE rotate_water_molecule

END MODULE mc_coordinates

!******************************************************************************

