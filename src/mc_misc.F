!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains miscellaneous subroutines used in the Monte Carlo runs,
!>      mostly I/O stuff
!> \author MJM
! *****************************************************************************
MODULE mc_misc
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE f77_blas
  USE force_env_types,                 ONLY: use_fist_force,&
                                             use_qs_force
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: &
       accattempt, get_mc_input_file, get_mc_molecule_info, get_mc_par, &
       mc_averages_type, mc_input_file_type, mc_molecule_info_type, &
       mc_moves_p_type, mc_moves_type, mc_simpar_type
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: final_mc_write,mc_averages_create,mc_averages_release,&
            mc_make_dat_file_new

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_misc'

CONTAINS

! *****************************************************************************
!> \brief initializes the structure that holds running averages of MC variables
!> \param averages the mc_averages strucutre you want to initialize
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE mc_averages_create ( averages  )

    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_averages_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! allocate all the structures...not sure why, but it won't work otherwise
      ALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "averages",0)

      averages%ave_energy=0.0E0_dp
      averages%ave_energy_squared=0.0E0_dp
      averages%ave_volume=0.0E0_dp
      averages%molecules=0.0E0_dp

! end the timing
       CALL timestop(handle)

END SUBROUTINE mc_averages_create

! *****************************************************************************
!> \brief deallocates the structure that holds running averages of MC variables
!> \param averages the mc_averages strucutre you want to release
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE mc_averages_release ( averages  )

    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_averages_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! deallocate
      DEALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "averages")

      NULLIFY(averages)

! end the timing
       CALL timestop(handle)

END SUBROUTINE mc_averages_release

! *****************************************************************************
!> \brief writes a bunch of simulation data to the specified unit
!> \param mc_par the mc parameters for the simulation
!> \param all_moves the structure that holds data on how many moves are
!>               accepted/rejected
!> \param iw the unit to write to
!> \param energy_check the sum of the energy changes of each move
!> \param initial_energy the initial unbiased energy of the system
!> \param final_energy the final unbiased energy of the system
!> \param averages the structure that holds computed average properites for
!>               the simulation
!>
!>    Only use in serial.
!> \author MJM
! *****************************************************************************
SUBROUTINE final_mc_write (mc_par,all_moves,iw,energy_check,initial_energy,&
                           final_energy,averages)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: all_moves
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: energy_check, initial_energy, &
                                                final_energy
    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(len=*), PARAMETER :: routineN = 'final_mc_write', &
      routineP = moduleN//':'//routineN

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

! begin the timing of the subroutine

    CALL timeset(routineN,handle)

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

! Hybrid MC moves (basically short MD runs)
    string2="Attempted       Accepted       Percent"
    string1="HMC Moves"
    CALL final_move_write(all_moves(1)%moves%hmc,string1,string2,iw)

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
    CALL timestop(handle)

END SUBROUTINE final_mc_write

! *****************************************************************************
!> \brief ...
!> \param move_data ...
!> \param string1 ...
!> \param string2 ...
!> \param iw ...
!> \param string3 ...
!> \param format_string ...
!> \param lbias ...
!> \param displacement ...
! *****************************************************************************
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

! *****************************************************************************
!> \brief writes a new input file that CP2K can read in for when we want
!>      to change a force env (change molecule number)...this is much simpler
!>      than the version I had used to have, and also more flexible (in a way).
!>      It assumes that &CELL comes before &COORDS, and &COORDS comes before
!>      &TOPOLOGY, and &TOPOLOGY comes before &GLOBAL (which comes before MC).
!>      It also assumes that you use &MOL_SET in &TOPOLOGY.  Still, many fewer
!>      assumptions than before.
!>
!>      box_length and coordinates should be passed in a.u.
!> \param coordinates the coordiantes of the atoms in the force_env (a.u.)
!> \param atom_names ...
!> \param nunits_tot the total number of atoms
!> \param box_length the length of all sides of the simulation box (angstrom)
!> \param filename the name of the file to write to
!> \param nchains ...
!> \param mc_input_file ...
!> \param fft_lib the FFT library we're using (FFTW,FFTSG)
!> UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author MJM
!> \note   Only use in serial.
! *****************************************************************************
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
      global_row_end, iline, in_use, itype, iunit, motion_row_end, &
      motion_row_start, nmol_types, nunits_empty, run_type_row, start_line, &
      unit
    INTEGER, DIMENSION(:), POINTER           :: mol_set_nmol_column, &
                                                mol_set_nmol_row
    REAL(dp), DIMENSION(:, :), POINTER       :: coordinates_empty

! open the file

      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! get all the information from the input_file_type
      CALL get_mc_input_file(mc_input_file,text=text,cell_row=cell_row,&
           cell_column=cell_column,coord_row_start=coord_row_start,&
           coord_row_end=coord_row_end,mol_set_nmol_row=mol_set_nmol_row,&
           mol_set_nmol_column=mol_set_nmol_column,global_row_end=global_row_end,&
           run_type_row=run_type_row,in_use=in_use,atom_names_empty=atom_names_empty,&
           nunits_empty=nunits_empty,coordinates_empty=coordinates_empty,&
           motion_row_start=motion_row_start,motion_row_end=motion_row_end)

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
         WRITE(unit,'(A)') '  RUN_TYPE     ENERGY_FORCE'
      CASE (use_qs_force)
         WRITE(unit,'(A)') '  RUN_TYPE     ENERGY_FORCE'
      END SELECT
      DO iline=run_type_row+1,global_row_end
         WRITE(unit,'(A)') TRIM(text(iline))
      ENDDO
      DO iline=motion_row_start,motion_row_end
         WRITE(unit,'(A)') TRIM(text(iline))
      ENDDO

! close the file
      CALL close_file(unit_number=unit)

    END SUBROUTINE MC_MAKE_DAT_FILE_NEW
END MODULE mc_misc

