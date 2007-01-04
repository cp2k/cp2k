!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_misc [1.0] *
!!
!!   NAME
!!     mc_misc
!!
!!   FUNCTION
!!     contains miscellaneous subroutines used in the Monte Carlo runs,
!!     mostly I/O stuff
!!
!!   AUTHOR
!!     MJM
!!
!******************************************************************************

MODULE mc_misc
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             get_fist_env
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force
  USE global_types,                    ONLY: DEBUG,&
                                             HIGH,&
                                             LOW,&
                                             MEDIUM,&
                                             SILENT
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: &
       accattempt, get_mc_input_file, get_mc_molecule_info, get_mc_par, &
       mc_averages_type, mc_input_file_type, mc_molecule_info_type, &
       mc_moves_p_type, mc_moves_type, mc_simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: atom_type,&
                                             bend_type,&
                                             bond_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             torsion_type
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PUBLIC :: final_mc_write,mc_averages_create,mc_averages_release,&
            mc_make_dat_file,mc_make_dat_file_new

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_misc'

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mc_misc.F/mc_averages_create*
!!
!!   NAME
!!     mc_averages_create
!!
!!   FUNCTION
!!     initializes the structure that holds running averages of MC variables
!!
!!   ARGUMENTS
!!     - averages: the mc_averages strucutre you want to initialize
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!****************************************************************************

SUBROUTINE mc_averages_create ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_averages_create'

    INTEGER                                  :: handle, stat

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! allocate all the structures...not sure why, but it won't work otherwise
      ALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "averages",0)

      averages%ave_energy=0.0E0_dp
      averages%ave_energy_squared=0.0E0_dp
      averages%ave_volume=0.0E0_dp
      averages%molecules=0.0E0_dp

! end the timing
       CALL timestop(0.0E0_dp,handle)


END SUBROUTINE mc_averages_create

!******************************************************************************
!!****** mc_misc.F/mc_averages_release
!!
!!   NAME
!!     mc_averages_release
!!
!!   FUNCTION
!!     deallocates the structure that holds running averages of MC variables
!!
!!   ARGUMENTS
!!     - averages: the mc_averages strucutre you want to release
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!****************************************************************************

SUBROUTINE mc_averages_release ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_averages_release'

    INTEGER                                  :: handle, stat

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! deallocate
      DEALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "averages")

      NULLIFY(averages)

! end the timing
       CALL timestop(0.0E0_dp,handle)


END SUBROUTINE mc_averages_release

!******************************************************************************
!!****s* mc_misc.F/final_mc_write
!!
!!   NAME
!!     final_mc_write
!!
!!   FUNCTION
!!     writes a bunch of simulation data to the specified unit
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the simulation
!!     - moves: the structure that holds data on how many moves are
!!              accepted/rejected
!!     - iw: the unit to write to
!!     - energy_check: the sum of the energy changes of each move
!!     - initial_energy: the initial unbiased energy of the system
!!     - final_energy: the final unbiased energy of the system
!!     - averages: the structure that holds computed average properites for
!!              the simulation
!!
!!   Only use in serial.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

SUBROUTINE final_mc_write (mc_par,all_moves,iw,energy_check,initial_energy,&
                           final_energy,averages)


    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: all_moves
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: energy_check, initial_energy, &
                                                final_energy
    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'final_mc_write'

    CHARACTER(LEN=5)                         :: molecule_string, tab_string
    CHARACTER(LEN=default_string_length)     :: format_string, string1, &
                                                string2, string3
    INTEGER                                  :: handle, itype, nmol_types
    LOGICAL                                  :: lbias
    REAL(dp), DIMENSION(:), POINTER          :: rmangle, rmbond, rmdihedral, &
                                                rmrot, rmtrans
    REAL(KIND=dp)                            :: pmswap, rmvolume
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(mc_moves_type), POINTER             :: moves

!------------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

    NULLIFY(mc_molecule_info,rmbond,rmangle,rmdihedral,rmrot,rmtrans)

    CALL get_mc_par(mc_par,pmswap=pmswap,rmvolume=rmvolume,&
      lbias=lbias,rmbond=rmbond,rmangle=rmangle,rmdihedral=rmdihedral,&
      rmtrans=rmtrans,rmrot=rmrot,mc_molecule_info=mc_molecule_info)
    CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types)
    WRITE(molecule_string,'(I2)') nmol_types
    WRITE(tab_string,'(I4)') 81-11*nmol_types
    format_string="(A,T" // TRIM(ADJUSTL(tab_string)) // "," // TRIM(ADJUSTL(molecule_string)) // "(2X,F9.6))"

! write out some data averaged over the whole simulation
    WRITE(iw,*)
    WRITE(iw,'(A,A)') '*****************************************************',&
                   '***************************'
    WRITE(iw,'(A,T66,F15.8)') "Average Energy [Hartrees]:", &
                averages%ave_energy
    IF (pmswap .GT. 0.0E0_dp) THEN
       WRITE(iw,'(A,T66,F15.8)') "Average number of molecules:",&
                   averages%molecules
    ENDIF
    WRITE(iw,'(A,A,T65,F16.6)') "Average Volume ",&
               "[angstroms**3]:",averages%ave_volume*angstrom**3

    WRITE(iw,*)

! write out acceptance rates for the moves

! volume moves
    WRITE(iw,'(A,A)') '-----------------------------------------------------',&
         '---------------------------'
    string2="Attempted       Accepted       Percent"
    string1="Volume Moves"
    string3="Maximum volume displacement [angstroms**3]= "
    rmvolume=rmvolume*angstrom**3
    CALL final_move_write(all_moves(1)%moves%volume,string1,string2,iw,&
         displacement=rmvolume,lbias=.FALSE.,format_string=format_string,&
         string3=string3)

! Quickstep moves (a series of moves with one potential, and then corrected for
! by another potential
    string2="Attempted       Accepted       Percent"
    string1="Quickstep Moves"
    CALL final_move_write(all_moves(1)%moves%Quickstep,string1,string2,iw)

    DO itype=1,nmol_types
       WRITE(iw,'(A,A)') '-----------------------------------------------------',&
                   '---------------------------'
       WRITE(iw,'(A,I5)') 'Move Data for Molecule Type ',itype
       WRITE(iw,'(A,A)') '-----------------------------------------------------',&
                   '---------------------------'

       moves => all_moves(itype)%moves

! AVBMC moves
       string2="Attempted       Accepted       Percent"
       string1="AVBMC moves from in to in"
       CALL final_move_write(moves%avbmc_inin,string1,string2,iw)
       string1="AVBMC moves from in to out"
       CALL final_move_write(moves%avbmc_inout,string1,string2,iw)
       string1="AVBMC moves from out to in"
       CALL final_move_write(moves%avbmc_outin,string1,string2,iw)
       string1="AVBMC moves from out to out"
       CALL final_move_write(moves%avbmc_outout,string1,string2,iw)

! conformation changes
       IF (moves%angle%attempts .GT. 0 .OR. &
            moves%bond%attempts .GT. 0 .OR. &
            moves%dihedral%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') "Conformational Moves",&
               "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
               moves%bond%attempts+moves%angle%attempts+&
               moves%dihedral%attempts,&
               moves%bond%successes+moves%angle%successes+&
               moves%dihedral%successes,&
               REAL(moves%bond%successes+moves%angle%successes+&
               moves%dihedral%successes,dp)/ &
               REAL(moves%bond%attempts+moves%angle%attempts+&
               moves%dihedral%attempts,dp)*100.0E0_dp
          string2="Attempted       Accepted       Percent"
          string1="Bond Changes"
          string3="Maximum bond displacement [angstroms]= "
          rmbond(itype)=rmbond(itype)*angstrom
          CALL final_move_write(moves%bond,string1,string2,iw,&
               displacement=rmbond(itype),lbias=lbias,format_string=format_string,&
               string3=string3)

          string1="Angle Changes"
          string3="Maximum angle displacement [degrees]= "
          rmangle(itype)=rmangle(itype)/pi*180.0E0_dp
          CALL final_move_write(moves%angle,string1,string2,iw,&
               displacement=rmangle(itype),lbias=lbias,format_string=format_string,&
               string3=string3)

          string1="Dihedral Changes"
          string3="Maximum dihedral displacement [degrees]= "
          rmdihedral(itype)=rmdihedral(itype)/pi*180.0E0_dp
          CALL final_move_write(moves%dihedral,string1,string2,iw,&
               displacement=rmdihedral(itype),lbias=lbias,format_string=format_string,&
               string3=string3)

          WRITE(iw,'(A,A,I5)') "Conformational Moves Rejected Because",&
               "Box Was Empty: ",moves%empty_conf
          WRITE(iw,'(A,A)') '-----------------------------------------------',&
               '--------------------------------'
       ENDIF

! translation moves
       string1="Translation Moves"
       string3="Maximum molecular translational displacement [angstroms]= "
       rmtrans(itype)=rmtrans(itype)*angstrom
       CALL final_move_write(moves%trans,string1,string2,iw,&
            displacement=rmtrans(itype),lbias=lbias,format_string=format_string,&
            string3=string3)

! rotation moves
          string1="Rotation Moves"
          string3="Maximum molecular rotational displacement [degrees]= "
          rmrot(itype)=rmrot(itype)/pi*180.0E0_dp
          CALL final_move_write(moves%rot,string1,string2,iw,&
               displacement=rmrot(itype),lbias=lbias,format_string=format_string,&
               string3=string3)

! swap moves
          IF (moves%swap%attempts .GT. 0) THEN
             WRITE(iw,'(A,T43,A)') "Swap Moves into this box",&
                  "Attempted       Empty          Percent"
             WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
                  moves%swap%attempts,&
                  moves%empty,&
                  REAL(moves%empty,dp)/ &
                  REAL(moves%swap%attempts,dp)*100.0E0_dp
             WRITE(iw,'(A,T43,A)') "                  Growths",&
                  "Attempted       Sucessful      Percent"
             WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
                  moves%swap%attempts,&
                  moves%grown,&
                  REAL(moves%grown,dp)/ &
                  REAL(moves%swap%attempts,dp)*100.0E0_dp
             WRITE(iw,'(A,T43,A)') "                    Total",&
                  "Attempted       Accepted       Percent"
             WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
                  moves%swap%attempts,&
                  moves%swap%successes,&
                  REAL(moves%swap%successes,dp)/ &
                  REAL(moves%swap%attempts,dp)*100.0E0_dp
             WRITE(iw,'(A,A)') '-----------------------------------------------',&
                  '--------------------------------'
          ENDIF
          
! now we write out information on the classical moves, if it's
! a classical simulations
          IF (lbias) THEN
             WRITE(iw,'(A)') "Biased Move Data"
             WRITE(iw,'(A,A)') '-------------------------------------------------',&
                  '-------------------------------'
             string2="Attempted       Accepted       Percent"
             string1="Bond Changes"
             string3="Maximum bond displacement [angstroms]= "
             CALL final_move_write(moves%bias_bond,string1,string2,iw,&
                  displacement=rmbond(itype),lbias=lbias,format_string=format_string,&
                  string3=string3)
             
             string1="Angle Changes"
             string3="Maximum angle displacement [degrees]= "
             CALL final_move_write(moves%bias_angle,string1,string2,iw,&
                  displacement=rmangle(itype),lbias=lbias,format_string=format_string,&
                  string3=string3)
             
             string1="Dihedral Changes"
             string3="Maximum dihedral displacement [degrees]= "
             CALL final_move_write(moves%bias_dihedral,string1,string2,iw,&
                  displacement=rmdihedral(itype),lbias=lbias,format_string=format_string,&
                  string3=string3)

 ! translation moves
             string1="Translation Moves"
             string3="Maximum molecular translational displacement [angstroms]= "
             CALL final_move_write(moves%bias_trans,string1,string2,iw,&
                  displacement=rmtrans(itype),lbias=lbias,format_string=format_string,&
                  string3=string3)

! rotation moves
             string1="Rotation Moves"
             string3="Maximum molecular rotational displacement [degrees]= "
             CALL final_move_write(moves%bias_rot,string1,string2,iw,&
                  displacement=rmrot(itype),lbias=lbias,format_string=format_string,&
                  string3=string3)

         ENDIF

       ENDDO



! see if the energies add up properly
    IF(ABS(initial_energy+energy_check-final_energy) .GT. 0.0000001E0_dp) &
         THEN
         WRITE(iw,*) '!!!!!!! We have an energy problem. !!!!!!!!'
         WRITE(iw,'(A,T64,F16.10)') 'Final Energy = ',final_energy
         WRITE(iw,'(A,T64,F16.10)') 'Inital Energy + energy_check =',&
                 initial_energy+energy_check
    ENDIF
    WRITE(iw,'(A,A)') '****************************************************',&
                   '****************************'
    WRITE(iw,*)

! end the timing
    CALL timestop(0.0E0_dp,handle)


END SUBROUTINE final_mc_write

SUBROUTINE final_move_write(move_data,string1,string2,iw,string3,&
     format_string,lbias,displacement)

    TYPE(accattempt), POINTER                :: move_data
    CHARACTER(default_string_length), &
      INTENT(IN)                             :: string1, string2
    INTEGER, INTENT(IN)                      :: iw
    CHARACTER(default_string_length), &
      INTENT(IN), OPTIONAL                   :: string3, format_string
    LOGICAL, INTENT(IN), OPTIONAL            :: lbias
    REAL(dp), OPTIONAL                       :: displacement

  IF(.NOT. PRESENT(format_string)) THEN
     IF (move_data%attempts .GT. 0) THEN
        WRITE(iw,'(A,T43,A)') TRIM(ADJUSTL(string1)),&
             TRIM(ADJUSTL(string2))
        WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             move_data%attempts,&
             move_data%successes,&
             REAL(move_data%successes,dp)/ &
             REAL(move_data%attempts,dp)*100.0E0_dp
        WRITE(iw,'(A,A)') '-----------------------------------------------',&
             '---------------------------------'
     ENDIF
  ELSE
     IF(.NOT. PRESENT(string3) .OR. .NOT. PRESENT(lbias) .OR. &
          .NOT. PRESENT(displacement)) THEN
        WRITE(iw,*) 'MISSING FLAGS IN FINAL_MOVE_WRITE'
     ENDIF
     IF (move_data%attempts .GT. 0 ) THEN
        WRITE(iw,'(A,T43,A)') TRIM(ADJUSTL(string1)),&
             TRIM(ADJUSTL(string2))
        WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             move_data%attempts,&
             move_data%successes,&
             REAL(move_data%successes,dp)/ &
             REAL(move_data%attempts,dp)*100.0E0_dp
        IF( .NOT. lbias) WRITE(iw,'(A,T71,F10.5)') &
             string3,displacement
        WRITE(iw,'(A,A)') '-----------------------------------------------',&
             '---------------------------------'
     ENDIF
  ENDIF

END SUBROUTINE final_move_write


!******************************************************************************
!!****s* mc_misc.F/mc_make_dat_file
!!
!!   NAME
!!     mc_make_dat_file
!!
!!   FUNCTION
!!     writes a new input file that CP2K can read in for when we want
!!     to change a force env (change molecules or cell length)
!!
!!   ARGUMENTS
!!     - coordinates: the coordiantes of the atoms in the force_env (a.u.)
!!     - natoms_tot: the total number of atoms
!!     - box_length: the length of all sides of the simulation box (angstrom)
!!     - filename: the name of the file to write to
!!     - force_env: the force environment that holds the parameters needed
!!                  to write on the file
!!     - fft_lib: the FFT library we're using (FFTW,FFTSG)
!!
!!   Only use in serial.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_make_dat_file(coordinates,natoms_tot,box_length,filename,&
      force_env,fft_lib,print_level,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: natoms_tot
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN)                             :: box_length
    CHARACTER(LEN=*), INTENT(IN)             :: filename
    TYPE(force_env_type), POINTER            :: force_env
    CHARACTER(LEN=*), INTENT(IN)             :: fft_lib
    INTEGER, INTENT(IN)                      :: print_level
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_make_dat_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: atom_a_done, atom_b_done, &
                                                atom_c_done, atom_d_done
    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: c_val, mol_name, name, &
                                                name_a, name_b, name_c, &
                                                name_d, unit_type
    INTEGER :: current_bend, current_bond, current_torsion, ewald_type, gmax, &
      handle, i_val, iatom, ibend, ibond, ichain, ifunct, imolecule, imul, &
      itorsion, itype, jbend, jbond, jtorsion, jtype, natoms, nbend, nbond, &
      nchains, ns_max, ntorsion, ntypes, o_spline, stat, unit
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp)                            :: alpha, charge, epsilon, r_val
    REAL(KIND=dp), DIMENSION(1:3)            :: abc
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: input, spline_section, &
                                                xc_fun, xc_fun_section
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list

!    TYPE(input_info_type), POINTER           :: inp_info
!------------------------------------------------------------------------------
! begin the timing of the subroutine

    STOP 'Need to rework this subroutine...mc_make_dat_file'
! deal with multiple molecules...perhaps by reading the whole input
! file into a character array and parsing to find where COORD and CELL are,
! then writing the new information with everything else around it...
! or making greater use of the input_parsing stuff

      CALL timeset(routineN,'I','',handle)

! nullify some stuff
      NULLIFY(qs_env,atomic_kind_set,atomic_kind,gth_potential,&
         orb_basis_set,elp_potential,input,xc_fun,kg_env,&
         bond_list,bend_list,torsion_list,spline_section)

! grab some of the structures we're gonna need
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
         CALL force_env_get(force_env,fist_env=fist_env,subsys=subsys,error=error)
         CALL get_fist_env(fist_env,cell=cell,ewald_env=ewald_env,&
            input=input,atomic_kind_set=atomic_kind_set,&
            fist_nonbond_env=fist_nonbond_env,error=error)
         CALL ewald_env_get(ewald_env,unit_type=unit_type,&
            ewald_type=ewald_type,alpha=alpha,gmax=gmax,o_spline=o_spline,&
            ns_max=ns_max,epsilon=epsilon,error=error)
         CALL get_cell(cell,abc=abc)
      CASE ( use_qs_force )
         CALL force_env_get(force_env,qs_env=qs_env,subsys=subsys,error=error)
         CALL get_qs_env(qs_env,use_ref_cell=use_ref_cell,cell_ref=cell_ref,&
           input=input,atomic_kind_set=atomic_kind_set,error=error)
         CALL get_cell(cell_ref,abc=abc)

      CASE ( use_kg_force )
         CALL force_env_get(force_env,kg_env=kg_env,subsys=subsys,error=error)
         CALL get_kg_env(kg_env,use_ref_cell=use_ref_cell,cell_ref=cell_ref,&
           input=input,atomic_kind_set=atomic_kind_set,error=error)
         CALL get_cell(cell_ref,abc=abc)

      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

! get some molecule information
      CALL cp_subsys_get(subsys(1)%subsys, &
         molecule_kinds_new=molecule_kinds_new,error=error)

! more than one molecule?
      IF(SIZE(molecule_kinds_new%els(:)) .GT. 1) THEN
         STOP 'Cannot do more than one molecule yet'
      ELSE
         molecule_kind => molecule_kinds_new%els(1)
         CALL get_molecule_kind(molecule_kind,atom_list=atom_list,&
              natom=natoms,bond_list=bond_list,bend_list=bend_list,&
              torsion_list=torsion_list,nbond=nbond,nbend=nbend,&
              ntorsion=ntorsion,name=mol_name)
      ENDIF

! find out how many atom types we have
      ntypes=SIZE(atomic_kind_set(:))

! open the file
      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! write all the stuff down
      WRITE(unit,10000) "&FORCE_EVAL"
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
         WRITE(unit,10000) ' METHOD      FIST'
      CASE ( use_qs_force )
         WRITE(unit,10000) ' METHOD      Quickstep'
      CASE ( use_kg_force )
         WRITE(unit,10000) ' METHOD      KG'
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
! grab some of the Ewald keywords
         WRITE(unit,10000) "  &MM"
         WRITE(unit,10000) "    &POISSON"
         WRITE(unit,10000) "      &EWALD"
         CALL convert_to_cp2k_units ( 'BOHR', length = alpha,&
            l_power = -1 )
         SELECT CASE ( ewald_type )
         CASE ( do_ewald_ewald )
            WRITE(unit,10000) "       EWALD_TYPE ewald"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       GMAX ",gmax
         CASE ( do_ewald_none )
            WRITE(unit,10000) "       EWALD_TYPE none"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       GMAX ",gmax
         CASE ( do_ewald_pme )
            WRITE(unit,10000) "       EWALD_TYPE pme"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       NS_MAX ",ns_max
            WRITE(unit,20000) "       EPSILON ",epsilon
         CASE ( do_ewald_spme )
            WRITE(unit,10000) "       EWALD_TYPE pme"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       GMAX ",gmax
            WRITE(unit,20003) "       O_SPLINE ",o_spline
         END SELECT
         WRITE(unit,10000) "      &END EWALD"
         WRITE(unit,10000) "    &END POISSON"


! rcut
         CALL fist_nonbond_env_get (fist_nonbond_env, potparm=potparm,error=error)

! now some of the force field stuff
         WRITE(unit,10000) "    &FORCEFIELD"

! need to keep the value of EMAX_SPLINE, in particular for the HF
! calculations (very large cutoff)
         spline_section => section_vals_get_subs_vals(input,&
              "MM%FORCEFIELD%SPLINE",error=error)
         CALL section_vals_val_get(spline_section,"EMAX_SPLINE",&
              r_val=r_val,error=error)
         WRITE(unit,10000) "      &SPLINE"
         WRITE(unit,20002) "        EMAX_SPLINE ",r_val
         CALL section_vals_val_get(spline_section,"EPS_SPLINE",&
              r_val=r_val,error=error)
         WRITE(unit,20002) "        EPS_SPLINE ",r_val
         WRITE(unit,10000) "      &END SPLINE"

! charges
         DO itype=1,ntypes
            atomic_kind => atomic_kind_set(itype)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
              name=name,elp_potential=elp_potential)
            CALL get_potential(potential=elp_potential,&
              qeff=charge)
            WRITE(unit,10000) "      &CHARGE"
            WRITE(unit,20001) "        ATOM ",TRIM(ADJUSTL(name))
            WRITE(unit,20002) "        CHARGE ",charge
            WRITE(unit,10000) "      &END CHARGE"
         ENDDO

! need to find all the bonds and write them out
         ALLOCATE(atom_a_done(1:nbond),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_b_done(1:nbond),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         atom_a_done(:)=""
         atom_b_done(:)=""
         current_bond=1
         outer_bond:DO ibond=1,nbond
            name_a=atom_list(bond_list(ibond)%a)%name
            name_b=atom_list(bond_list(ibond)%b)%name
            DO jbond=1,current_bond
               IF(atom_a_done(jbond)==name_a .AND. atom_b_done(jbond)==name_b&
                    .OR. atom_a_done(jbond)==name_b .AND. atom_b_done(jbond)==name_a)THEN
                  CYCLE outer_bond
               ENDIF
            ENDDO
            atom_a_done(current_bond)=name_a
            atom_b_done(current_bond)=name_b
            current_bond=current_bond+1
            WRITE(unit,10000) "      &BOND"
            WRITE(unit,30002) "        ATOMS ",TRIM(ADJUSTL(name_a)),TRIM(ADJUSTL(name_b))
            WRITE(unit,20002) "        K ",bond_list(ibond)%bond_kind%k(1)
            WRITE(unit,20002) "        R0 ",bond_list(ibond)%bond_kind%r0
            WRITE(unit,10000) "      &END BOND"
         ENDDO outer_bond
         DEALLOCATE(atom_a_done)
         DEALLOCATE(atom_b_done)

! same thing for bends
         ALLOCATE(atom_a_done(1:nbend),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_b_done(1:nbend),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_c_done(1:nbend),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         atom_a_done(:)=""
         atom_b_done(:)=""
         atom_c_done(:)=""
         current_bend=1
         outer_bend:DO ibend=1,nbend
            name_a=atom_list(bend_list(ibend)%a)%name
            name_b=atom_list(bend_list(ibend)%b)%name
            name_c=atom_list(bend_list(ibend)%c)%name
            DO jbend=1,current_bend
               IF(atom_a_done(jbend)==name_a .AND. atom_b_done(jbend)==name_b &
                    .AND.atom_c_done(jbend)==name_c &
                    .OR. atom_a_done(jbend)==name_c .AND. atom_b_done(jbend)==name_b .AND. &
                    atom_c_done(jbend)==name_a)THEN
                  CYCLE outer_bend
               ENDIF
            ENDDO
            atom_a_done(current_bend)=name_a
            atom_b_done(current_bend)=name_b
            atom_c_done(current_bend)=name_c
            current_bend=current_bend+1
            WRITE(unit,10000) "      &BEND"
            WRITE(unit,40004) "        ATOMS ",TRIM(ADJUSTL(name_a)),TRIM(ADJUSTL(name_b)),TRIM(ADJUSTL(name_c))
            WRITE(unit,20002) "        K ",bend_list(ibend)%bend_kind%k
            WRITE(unit,20002) "        THETA0 ",bend_list(ibend)%bend_kind%theta0
            WRITE(unit,10000) "      &END BEND"
         ENDDO outer_bend
         DEALLOCATE(atom_a_done)
         DEALLOCATE(atom_b_done)
         DEALLOCATE(atom_c_done)

! finally, write out the torsions
         ALLOCATE(atom_a_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_b_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_c_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_d_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         atom_a_done(:)=""
         atom_b_done(:)=""
         atom_c_done(:)=""
         atom_d_done(:)=""
         current_torsion=1
         outer_torsion:DO itorsion=1,ntorsion
            name_a=atom_list(torsion_list(itorsion)%a)%name
            name_b=atom_list(torsion_list(itorsion)%b)%name
            name_c=atom_list(torsion_list(itorsion)%c)%name
            name_d=atom_list(torsion_list(itorsion)%d)%name
            DO jtorsion=1,current_torsion
               IF(atom_a_done(jtorsion)==name_a .AND. atom_b_done(jtorsion)==name_b &
                    .AND.atom_c_done(jtorsion)==name_c .AND. atom_d_done(jtorsion)==name_d &
                    .OR. atom_a_done(jtorsion)==name_d .AND. atom_b_done(jtorsion)==name_c .AND. &
                    atom_c_done(jtorsion)==name_b .AND. atom_d_done(jtorsion)==name_a)THEN
                  CYCLE outer_torsion
               ENDIF
            ENDDO
            atom_a_done(current_torsion)=name_a
            atom_b_done(current_torsion)=name_b
            atom_c_done(current_torsion)=name_c
            atom_d_done(current_torsion)=name_d
            current_torsion=current_torsion+1
            DO imul=1,torsion_list(itorsion)%torsion_kind%nmul
               WRITE(unit,10000) "      &TORSION"
               WRITE(unit,50001) "        ATOMS ",TRIM(ADJUSTL(name_a)),TRIM(ADJUSTL(name_b)),&
                    TRIM(ADJUSTL(name_c)),TRIM(ADJUSTL(name_d))
               WRITE(unit,20003) "        M ",torsion_list(itorsion)%torsion_kind%m(imul)
               WRITE(unit,20002) "        K ",torsion_list(itorsion)%torsion_kind%k(imul)
               WRITE(unit,20002) "        PHI0 ",torsion_list(itorsion)%torsion_kind%cosphi0(imul)
               WRITE(unit,10000) "      &END TORSION"
            ENDDO

         ENDDO outer_torsion
         DEALLOCATE(atom_a_done)
         DEALLOCATE(atom_b_done)
         DEALLOCATE(atom_c_done)
         DEALLOCATE(atom_d_done)


! nonbonded...only works for LJ at the moment
         WRITE(unit,10000) "      &NONBONDED"
         DO itype=1,SIZE(potparm%pot,1)
            DO jtype=itype,SIZE(potparm%pot,2)
               WRITE(unit,10000) "        &LENNARD-JONES"
               WRITE(unit,30002) "          ATOMS ",&
                  TRIM(ADJUSTL(potparm%pot(itype,jtype)%pot%at1)),&
                  TRIM(ADJUSTL(potparm%pot(itype,jtype)%pot%at2))
               WRITE(unit,20002) "          EPSILON ",&
                  potparm%pot(itype,jtype)%pot%lj%epsilon*kelvin
               WRITE(unit,20002) "          SIGMA ",&
                  potparm%pot(itype,jtype)%pot%lj%sigma6**(1.0_dp/6.0_dp)*&
                  angstrom
               WRITE(unit,20000) "          RCUT ",&
                  potparm%pot(itype,jtype)%pot%rcutsq**0.5_dp*angstrom
               WRITE(unit,10000) "        &END LENNARD-JONES"
            ENDDO
         ENDDO
         WRITE(unit,10000) "      &END NONBONDED"
         WRITE(unit,10000) "    &END FORCEFIELD"
         WRITE(unit,10000) "  &END MM"


      CASE ( use_qs_force )
         WRITE(unit,10000) "  &DFT"
         WRITE(unit,10000) "    &MGRID"
         CALL section_vals_val_get(input,"DFT%MGRID%CUTOFF",r_val=r_val,error=error)
         WRITE(unit,20000) "      CUTOFF  ",2.0E0_dp*r_val
         CALL section_vals_val_get(input,"DFT%MGRID%NGRIDS",i_val=i_val,error=error)
         WRITE(unit,20003) "      NGRIDS  ",i_val
         WRITE(unit,10000) "    &END MGRID"
         WRITE(unit,10000) "    &QS"
         WRITE(unit,10000) "      EXTRAPOLATION USE_PREV_WF"
         WRITE(unit,10000) "    &END QS"
         WRITE(unit,10000) "    &SCF"
         CALL section_vals_val_get(input,"DFT%SCF%SCF_GUESS",c_val=c_val,error=error)
         WRITE(unit,20001) "      SCF_GUESS  ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "    &END SCF"
         WRITE(unit,10000) "    &XC"
         xc_fun_section => section_vals_get_subs_vals(input,&
            "DFT%XC%XC_FUNCTIONAL",error=error)
         WRITE(unit,10000) "      &XC_FUNCTIONAL"
         DO ifunct=1,4
            xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifunct,error=error)
            IF(ASSOCIATED(xc_fun)) THEN
               WRITE(unit,20001) "        &",TRIM(xc_fun%section%name)
               WRITE(unit,20001) "        &END ",TRIM(xc_fun%section%name)
            ENDIF
         ENDDO
         WRITE(unit,10000) "      &END XC_FUNCTIONAL"
         WRITE(unit,10000) "      &XC_GRID"
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_SMOOTH_RHO",&
           c_val=c_val,error=error)
         WRITE(unit,20001) "        XC_SMOOTH_RHO ",TRIM(ADJUSTL(c_val))
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_DERIV",&
           c_val=c_val,error=error)
         WRITE(unit,20001) "        XC_DERIV ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "      &END XC_GRID"
         WRITE(unit,10000) "    &END XC"
         WRITE(unit,10000) "  &END DFT"
      CASE ( use_kg_force )
         WRITE(unit,10000) "  &DFT"
         WRITE(unit,10000) "    &MGRID"
         CALL section_vals_val_get(input,"DFT%MGRID%CUTOFF",r_val=r_val,error=error)
         WRITE(unit,20000) "      CUTOFF  ",2.0E0_dp*r_val
         CALL section_vals_val_get(input,"DFT%MGRID%NGRIDS",i_val=i_val,error=error)
         WRITE(unit,20003) "      NGRIDS  ",i_val
         WRITE(unit,10000) "    &END MGRID"
         WRITE(unit,10000) "    &QS"
         CALL section_vals_val_get(input,"DFT%QS%EXTRAPOLATION",c_val=c_val,error=error)
         WRITE(unit,20001) "      EXTRAPOLATION  ",TRIM(ADJUSTL(c_val))
         CALL section_vals_val_get(input,"DFT%QS%EPS_CORE_CHARGE",r_val=r_val,error=error)
         WRITE(unit,20004) "      EPS_CORE_CHARGE  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_PPL",r_val=r_val,error=error)
         WRITE(unit,20004) "      EPS_PPL  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_PGF_ORB",r_val=r_val,error=error)
         WRITE(unit,20004) "      EPS_PGF_ORB  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_GVG_RSPACE",r_val=r_val,error=error)
         WRITE(unit,20004) "      EPS_GVG_RSPACE  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_RHO",r_val=r_val,error=error)
         WRITE(unit,20004) "      EPS_RHO  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%METHOD",c_val=c_val,error=error)
         WRITE(unit,20001) "      METHOD  ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "    &END QS"
         WRITE(unit,10000) "    &SCF"
         CALL section_vals_val_get(input,"DFT%SCF%SCF_GUESS",c_val=c_val,error=error)
         WRITE(unit,20001) "      SCF_GUESS  ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "    &END SCF"
         WRITE(unit,10000) "    &XC"
         xc_fun_section => section_vals_get_subs_vals(input,&
            "DFT%XC%XC_FUNCTIONAL",error=error)
! take care of the functionals...I had to put the extra if statement in for the kinetic
! energy functional, as that doesn't just use the section name
         WRITE(unit,10000) "      &XC_FUNCTIONAL"
         DO ifunct=1,4
            xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifunct,error=error)
            IF(ASSOCIATED(xc_fun)) THEN
               WRITE(unit,20001) "        &",TRIM(xc_fun%section%name)
               IF(TRIM(xc_fun%section%name) == "KE_GGA") THEN
                  CALL section_vals_val_get(input,"DFT%XC%XC_FUNCTIONAL%KE_GGA%FUNCTIONAL",&
                       c_val=c_val,error=error)
                  WRITE(unit,20001) "          FUNCTIONAL ",TRIM(c_val)
               ENDIF
               WRITE(unit,20001) "        &END ",TRIM(xc_fun%section%name)
            ENDIF
         ENDDO
         WRITE(unit,10000) "      &END XC_FUNCTIONAL"
         WRITE(unit,10000) "      &XC_GRID"
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_SMOOTH_RHO",&
           c_val=c_val,error=error)
         WRITE(unit,20001) "        XC_SMOOTH_RHO ",TRIM(ADJUSTL(c_val))
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_DERIV",&
           c_val=c_val,error=error)
         WRITE(unit,20001) "        XC_DERIV ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "      &END XC_GRID"
         WRITE(unit,10000) "    &END XC"
         WRITE(unit,10000) "  &END DFT"
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

      WRITE(unit,10000) "  &SUBSYS"
! cell information
      WRITE(unit,10000) "    &CELL"
      WRITE(unit,40000) '      ABC  ',box_length(1:3)
      WRITE(unit,10000) '      UNIT ANGSTROM'
      IF(force_env%in_use == use_qs_force .OR. &
        force_env%in_use == use_kg_force) THEN
         IF(use_ref_cell) THEN
            WRITE(unit,10000) '       &CELL_REF'
            WRITE(unit,40000) '         ABC ',&
               abc(1:3)*angstrom
            WRITE(unit,10000) '         UNIT ANGSTROM'
            WRITE(unit,10000) '       &END CELL_REF'
         ENDIF
      ENDIF
      WRITE(unit,10000) '    &END CELL'


! write out the coordinate information...currently only works for
! single component systems
      WRITE(unit,10000) '    &COORD'
      IF(MOD(natoms_tot,natoms) .NE. 0) THEN
         CALL stop_program("mc_make_dat_file",&
              "Not a single component system?")
      ENDIF
      nchains=natoms_tot/natoms

      IF(nchains .NE. 0) THEN
         DO ichain=1,nchains
            DO iatom=1,natoms
               WRITE(unit,40003) TRIM(ADJUSTL(atom_list(iatom)%name)),&
               coordinates(1:3,(ichain-1)*natoms+iatom)*angstrom
            ENDDO
         ENDDO
      ELSE
         DO iatom=1,natoms
            WRITE(unit,40003) TRIM(ADJUSTL(atom_list(iatom)%name)),&
               REAL(iatom,dp),REAL(iatom,dp),REAL(iatom,dp)
         ENDDO
      ENDIF
      WRITE(unit,10000) '    &END COORD'

! now the kind types, for the basis sets and potentials
      IF(force_env%in_use == use_qs_force .OR. &
         force_env%in_use == use_kg_force)THEN
         DO itype=1,ntypes
            atomic_kind => atomic_kind_set(itype)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
            element_symbol=element_symbol,gth_potential=gth_potential,&
            orb_basis_set=orb_basis_set,elp_potential=elp_potential)

            CALL get_gto_basis_set(orb_basis_set,name=name)

            WRITE(unit,20001) "    &KIND ",element_symbol
            WRITE(unit,20001) "      BASIS_SET    ",TRIM(ADJUSTL(name))

            IF(ASSOCIATED(gth_potential)) THEN
               CALL get_potential(gth_potential,name=name)
               WRITE(unit,20001) "      POTENTIAL    ",TRIM(ADJUSTL(name))
            ELSE
               CALL get_potential(elp_potential,name=name)
               WRITE(unit,20001) "      POTENTIAL    ",TRIM(ADJUSTL(name))
            ENDIF
            WRITE(unit,10000) "    &END KIND"
         ENDDO
      ENDIF
! write the topology section, using MOL_SET
      WRITE(unit,10000) "    &TOPOLOGY"
      WRITE(unit,10000) "      CONNECTIVITY MOL_SET"
      DO imolecule=1,1
        WRITE(unit,10000) "      &MOL_SET"
        SELECT CASE ( force_env%in_use )
        CASE ( use_fist_force )
           WRITE(unit,10000) '        CONN_FILE topology_fist_' &
                // TRIM(ADJUSTL(mol_name)) // '.psf'
        CASE ( use_qs_force )
           WRITE(unit,10000) '        CONN_FILE topology_atoms_' &
                // TRIM(ADJUSTL(mol_name)) // '.psf'
        CASE ( use_kg_force )
           WRITE(unit,10000) '        CONN_FILE topology_atoms_' &
                // TRIM(ADJUSTL(mol_name)) // '.psf'
        CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END SELECT
        WRITE(unit,20003) "        NMOL ",nchains
        WRITE(unit,10000) "      &END MOL_SET"
      ENDDO
      WRITE(unit,10000) "    &END TOPOLOGY"
      WRITE(unit,10000) "  &END SUBSYS"
      WRITE(unit,10000) "&END FORCE_EVAL"
      WRITE(unit,10000) "&GLOBAL"

      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
         WRITE(unit,10000) ' RUN_TYPE     DEBUG'
      CASE (use_qs_force)
         WRITE(unit,10000) ' RUN_TYPE     ENERGY_FORCE'
      CASE (use_kg_force)
         WRITE(unit,10000) ' RUN_TYPE     ENERGY_FORCE'
      END SELECT

      WRITE(unit,10000) ' PROJECT      H2O_MC'
      WRITE(unit,20001) ' FFTLIB       ',TRIM(ADJUSTL(fft_lib))
! need to write the correct print level
      SELECT CASE (print_level)
      CASE(SILENT)
         WRITE(unit,10000) ' PRINT_LEVEL SILENT'
      CASE(LOW)
         WRITE(unit,10000) ' PRINT_LEVEL LOW'
      CASE(MEDIUM)
         WRITE(unit,10000) ' PRINT_LEVEL MEDIUM'
      CASE(HIGH)
         WRITE(unit,10000) ' PRINT_LEVEL HIGH'
      CASE(DEBUG)
         WRITE(unit,10000) ' PRINT_LEVEL DEBUG'
      END SELECT
      WRITE(unit,10000) "&END"

! close the file
      CALL close_file(unit_number=unit)

! specify some formats
10000 FORMAT(A)
20000 FORMAT(A,F8.2)
20001 FORMAT(A,A)
20002 FORMAT(A,F18.12)
20003 FORMAT(A,I10)
20004 FORMAT(A,E12.6)
30000 FORMAT(A,F12.8,2X,I8)
30001 FORMAT(A,A,F16.10)
30002 FORMAT(A,A,1X,A)
40000 FORMAT(A,3(F20.12,2X))
40001 FORMAT(A,F12.8,2X,I8,2X,I8)
40002 FORMAT(A,F12.8,2X,I8,2X,F12.8)
40003 FORMAT(5X,A,3(F20.15,2X))
40004 FORMAT(A,A,1X,A,1X,A)
50000 FORMAT(A,3(F18.12,2X),A)
50001 FORMAT(A,A,1X,A,1X,A,1X,A)
60000 FORMAT(A,A,2X,A,2X,A,3(2X,F10.5))

! end the timing
      CALL timestop(0.0E0_dp,handle)

!STOP

END SUBROUTINE mc_make_dat_file

!******************************************************************************
!!****s* mc_misc.F/mc_make_dat_file
!!
!!   NAME
!!     mc_make_dat_file
!!
!!   FUNCTION
!!     writes a new input file that CP2K can read in for when we want
!!     to change a force env (change molecule number)...this is much simpler
!!     than the version I had used to have, and also more flexible (in a way).
!!     It assumes that &CELL comes before &COORDS, and &COORDS comes before
!!     &TOPOLOGY, and &TOPOLOGY comes before &GLOBAL (which comes before MC).
!!     It also assumes that you use &MOL_SET in &TOPOLOGY.  Still, many fewer
!!     assumptions than before.
!!
!!     box_length and coordinates should be passed in a.u.
!!
!!   ARGUMENTS
!!     - coordinates: the coordiantes of the atoms in the force_env (a.u.)
!!     - natoms_tot: the total number of atoms
!!     - box_length: the length of all sides of the simulation box (angstrom)
!!     - filename: the name of the file to write to
!!     - force_env: the force environment that holds the parameters needed
!!                  to write on the file
!!     - fft_lib: the FFT library we're using (FFTW,FFTSG)
!!
!!   Only use in serial.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_make_dat_file_new(coordinates,atom_names,nunits_tot,&
     box_length,filename,nchains,mc_input_file)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(IN)                             :: atom_names
    INTEGER, INTENT(IN)                      :: nunits_tot
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN)                             :: box_length
    CHARACTER(LEN=*), INTENT(IN)             :: filename
    INTEGER, DIMENSION(:), INTENT(IN)        :: nchains
    TYPE(mc_input_file_type), POINTER        :: mc_input_file

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_make_dat_file_new', &
      routineP = moduleN//':'//routineN

    CHARACTER(60)                            :: cell_string, mol_string
    CHARACTER(default_string_length)         :: line_text
    CHARACTER(default_string_length), &
      DIMENSION(:), POINTER                  :: atom_names_empty, text
    INTEGER :: cell_column, cell_row, coord_row_end, coord_row_start, &
      global_row_end, iline, in_use, itype, iunit, nmol_types, nunits_empty, &
      run_type_row, start_line, unit
    INTEGER, DIMENSION(:), POINTER           :: mol_set_nmol_column, &
                                                mol_set_nmol_row
    REAL(dp), DIMENSION(:, :), POINTER       :: coordinates_empty

!*************************************************************************
! open the file

      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! get all the information from the input_file_type
      CALL get_mc_input_file(mc_input_file,text=text,cell_row=cell_row,&
           cell_column=cell_column,coord_row_start=coord_row_start,&
           coord_row_end=coord_row_end,mol_set_nmol_row=mol_set_nmol_row,&
           mol_set_nmol_column=mol_set_nmol_column,global_row_end=global_row_end,&
           run_type_row=run_type_row,in_use=in_use,atom_names_empty=atom_names_empty,&
           nunits_empty=nunits_empty,coordinates_empty=coordinates_empty)

! how many molecule types?
      nmol_types=SIZE(nchains)

! first, write all the information up to the cell lengths
      DO iline=1,cell_row-1
         WRITE(unit,'(A)') TRIM(text(iline))
      ENDDO
! substitute in the current cell lengths
      WRITE(cell_string,'(3(F13.8,2X))') box_length(1:3)*angstrom
      line_text=text(cell_row)
      line_text(cell_column:cell_column+50)=cell_string(1:51)
      WRITE(unit,'(A)') TRIM(line_text)

! now write everything until the coordinates
      DO iline=cell_row+1,coord_row_start
         WRITE(unit,'(A)') TRIM(text(iline))
      ENDDO      
    
      
! we may pass nunits_tot=0, but we should still have coordinates
      IF(nunits_tot == 0) THEN
         DO iunit=1,nunits_empty
            WRITE(unit,'(5X,A,2X,3(F15.10))') &
                 TRIM(ADJUSTL(atom_names_empty(iunit))),&
                 coordinates_empty(1:3,iunit)*angstrom            
         ENDDO
      ELSE
         DO iunit=1,nunits_tot
            WRITE(unit,'(5X,A,2X,3(F15.10))') &
                 TRIM(ADJUSTL(atom_names(iunit))),&
                 coordinates(1:3,iunit)*angstrom
         ENDDO
      ENDIF

! now we need to write the MOL_SET section
      start_line=coord_row_end
      DO itype=1,nmol_types
         DO iline=start_line,mol_set_nmol_row(itype)-1
            WRITE(unit,'(A)') TRIM(text(iline))
         ENDDO

! have to print out one molecule, even if it's empty
         IF(nunits_tot == 0 .AND. itype == 1) THEN
            WRITE(mol_string,'(I8)') 1
         ELSE
            WRITE(mol_string,'(I8)') nchains(itype)
         ENDIF

         line_text=text(mol_set_nmol_row(itype))        
         line_text(mol_set_nmol_column(itype):mol_set_nmol_column(itype)+9)=&
              mol_string(1:10)
         WRITE(unit,'(A)') TRIM(line_text)
         start_line=mol_set_nmol_row(itype)+1
      ENDDO

! write up to the RUN_TYPE...tailor this for the type of environment, so
! that we can easily do ./cp2k.sdbg input.dat and have it run
      DO iline=mol_set_nmol_row(nmol_types)+1,run_type_row-1
         WRITE(unit,'(A)') TRIM(text(iline))
      ENDDO 
      SELECT CASE ( in_use )
      CASE ( use_fist_force )
         WRITE(unit,'(A)') '  RUN_TYPE     DEBUG'
      CASE (use_qs_force)
         WRITE(unit,'(A)') '  RUN_TYPE     ENERGY_FORCE'
      CASE (use_kg_force)
         WRITE(unit,'(A)') '  RUN_TYPE     ENERGY_FORCE'
      END SELECT
      DO iline=run_type_row+1,global_row_end
         WRITE(unit,'(A)') TRIM(text(iline))
      ENDDO 

! close the file
      CALL close_file(unit_number=unit)


    END SUBROUTINE MC_MAKE_DAT_FILE_NEW
END MODULE mc_misc

!******************************************************************************

