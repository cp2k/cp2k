!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mc [1.0] *
!!
!!   NAME
!!     mc_misc
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_misc
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsys_get
  USE kinds,                           ONLY: dp,dbl
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_move_control,                 ONLY: mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE termination,                     ONLY: stop_memory


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_water_data,final_mc_write,mc_classical_energy,&
            mc_averages_type,init_mc_averages

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_misc"

  TYPE mc_averages_type
      REAL(KIND = dp) :: ave_energy
      REAL(KIND = dp) :: ave_energy_squared
      REAL(KIND = dp) :: ave_cell_length
      REAL(KIND = dp) :: bond1
      REAL(KIND = dp) :: bond2
      REAL(KIND = dp) :: angle
  END TYPE mc_averages_type

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mc/init_mc_averages [1.0] *
!!
!!   NAME
!!     init_mc_averages
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     initializes the array that holds running averages of variables
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!*** ************************************************************************

SUBROUTINE init_mc_averages ( averages  )

  IMPLICIT NONE

! Arguments
      TYPE(mc_averages_type), POINTER :: averages

!   *** Local parameters ***

      CHARACTER(LEN=*), PARAMETER :: routineN = "init_mc_averages"

! Locals
      INTEGER :: handle,stat

!------------------------------------------------------------------------------

! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

! allocate all the structures...not sure why, but it won't work otherwise
      ALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "averages",0)

      averages%ave_energy=0.0d0
      averages%ave_energy_squared=0.0d0
      averages%ave_cell_length=0.0d0
      averages%bond1=0.0d0
      averages%bond2=0.0d0
      averages%angle=0.0d0

! end the timing
       CALL timestop(0.0d0,handle)


END SUBROUTINE init_mc_averages

!******************************************************************************
!!****** mc/get_water_data [1.0] *
!!
!!   NAME
!!     get_water_data
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     gets the running average of both bond lenths and the bond
!!     angle of one water molecule in the system
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!*** ************************************************************************

SUBROUTINE get_water_data ( force_env,averages,nnstep,nchain,data  )

  IMPLICIT NONE

! Arguments
      TYPE(mc_averages_type), POINTER :: averages
      TYPE(force_env_type), POINTER :: force_env
      INTEGER, INTENT(IN)               :: nnstep,nchain,data

!   *** Local parameters ***

      CHARACTER(LEN=*), PARAMETER :: routineN = "get_water_data"

! Locals
      INTEGER                           :: i,imol,iimol,handle
      TYPE(particle_list_type), POINTER :: particles
      TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: oldsys
      TYPE(cp_error_type) :: error
      REAL(KIND = dp)       :: tot_squared,bond1_loc,&
                               bond2_loc,angle_loc
      REAL(KIND=dp),DIMENSION(1:3)    :: vec_1,vec_2
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
              *180.0d0/3.14159265358979323

       ENDDO
       averages%bond1=averages%bond1*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+bond1_loc/REAL(nnstep,dbl)
       averages%bond2=averages%bond2*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+bond2_loc/REAL(nnstep,dbl)
       averages%angle=averages%angle*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+angle_loc/REAL(nnstep,dbl)

       WRITE(data,*) nnstep,bond1_loc,bond2_loc,angle_loc 
 
! end the timing
       CALL timestop(0.0d0,handle)


END SUBROUTINE get_water_data

!******************************************************************************

!******************************************************************************
!!****** mc/final_mc_write[1.0] *
!!
!!   NAME
!!     final_mc_write
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     prints out the final data of the simulation in a nice format
!!     
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!*** ************************************************************************

SUBROUTINE final_mc_write (mc_par,moves,iw,energy_check,initial_energy,&
                           old_energy,averages)

  IMPLICIT NONE

! Arguments
      REAL(KIND = dp),INTENT(IN)            :: energy_check,&
                                               initial_energy,&
                                               old_energy
      TYPE(mc_averages_type),POINTER  :: averages
      TYPE(mc_moves_type),POINTER  :: moves
      TYPE(mc_simulation_parameters_type), POINTER :: mc_par
      INTEGER, INTENT(IN)  :: iw

!   *** Local parameters ***

      CHARACTER(LEN=*), PARAMETER :: routineN = "final_mc_write"

!   *** Local variables ***
      INTEGER :: handle
 
!------------------------------------------------------------------------------

! begin the timing of the subroutine
    CALL timeset(routineN,'I','',handle)

    WRITE(iw,*)
    WRITE(iw,*) '****************************************************',&
                   '***************************'
    WRITE(iw,'(A,T66,F15.8)') " Average Energy [Hartrees]:", &
                averages%ave_energy
    WRITE(iw,'(A,T66,E15.8)') " Energy Fluctuation [Hartrees]:", &
            DSQRT(averages%ave_energy_squared-averages%ave_energy**2)
    WRITE(iw,'(A,A,T66,F15.8)') " Heat Capacity at ",&
          "constant volume [J/mol/K]",&
          DSQRT(averages%ave_energy_squared-averages%ave_energy**2)&
          *(4.3597482d-18)**2/(mc_par%temp)**2/1.380658d-23*6.022d23
    WRITE(iw,'(A,T66,F15.8)') " Average O-H distance 1 [angstroms]:",&
                   averages%bond1*0.529177d0
    WRITE(iw,'(A,T66,F15.8)') " Average O-H distance 2 [angstroms]:",&
                   averages%bond2*0.529177d0
    WRITE(iw,'(A,T66,F15.8)') " Average H-O-H angle [degrees]:",&
                   averages%angle
    IF (averages%ave_cell_length .NE. 0.0d0) THEN
       WRITE(iw,'(A,A,T68,F13.10)') " Average Density ",&
               "[g/cm^3]:",1.0d0/(averages%ave_cell_length)**3&
                      *mc_par%nchain/&
                      (0.529177d0)**3*(1.0d8)**3/6.022d23*18.0184d0
       WRITE(iw,'(A,A,T67,F14.10)') " Average Cell Length ",&
               "[angstroms]:",averages%ave_cell_length*0.529177d0
    ENDIF
    WRITE(iw,*)
    WRITE(iw,*) "Move Data"
    WRITE(6,*) '----------------------------------------------------',&
                   '---------------------------'
    IF (moves%volume%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Volume Moves",&
          "Attempted       Accepted       Percent"
       write(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%volume%attempts,&
             moves%volume%successes,&
             real(moves%volume%successes,dbl)/ &
             real(moves%volume%attempts,dbl)*100.0d0
       write(iw,'(A,T71,F10.5)') &
           ' Maximum Volume displacement [angstroms^3] = '&
                     ,mc_par%rmvolume*(0.529177d0)**3
       WRITE(6,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%transion%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Transion Moves",&
          "Attempted       Accepted       Percent"
       write(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%transion%attempts,&
             moves%transion%successes,&
             real(moves%transion%successes,dbl)/ &
             real(moves%transion%attempts,dbl)*100.0d0
       WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum atomic translational displacement [angstroms]= '&
                  ,mc_par%rmion*0.529177d0
       WRITE(6,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%trans%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Translation Moves",&
          "Attempted       Accepted       Percent"
       write(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%trans%attempts,&
             moves%trans%successes,&
             real(moves%trans%successes,dbl)/ &
             real(moves%trans%attempts,dbl)*100.0d0
       WRITE(iw,'(A,T71,F10.5)') &
          ' Maximum molecular translational displacement [angstroms]= '&
                  ,mc_par%rmtrans*0.529177d0
       WRITE(6,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%rot%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Rotation Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%rot%attempts,&
             moves%rot%successes,&
             real(moves%rot%successes,dbl)/ &
             real(moves%rot%attempts,dbl)*100.0d0
       WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [radians] = ',mc_par%rmrot
       WRITE(6,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%Quickstep%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Quickstep Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%Quickstep%attempts,&
             moves%Quickstep%successes,&
             real(moves%Quickstep%successes,dbl)/ &
             real(moves%Quickstep%attempts,dbl)*100.0d0
       WRITE(6,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF(ABS(initial_energy+energy_check-old_energy) .GT. 0.0000001d0) &
         THEN
         WRITE(iw,*) '!!!!!!! We have an energy problem. !!!!!!!!'
         WRITE(iw,'(A,T64,F16.10)') 'Final Energy = ',old_energy
         WRITE(iw,'(A,T64,F16.10)') 'Inital Energy + energy_check =',&
                 initial_energy+energy_check
    ENDIF
    WRITE(iw,*) '****************************************************',&
                   '***************************'
    WRITE(iw,*)
 
! end the timing
    CALL timestop(0.0d0,handle)

 
END SUBROUTINE final_mc_write

!******************************************************************************
SUBROUTINE mc_classical_energy ( force_env, classical_energy_new,&
                                 mc_par )

  IMPLICIT NONE

! Arguments
      TYPE(force_env_type), POINTER :: force_env
      REAL(KIND = dp),INTENT(OUT) :: classical_energy_new
      TYPE(mc_simulation_parameters_type),POINTER :: mc_par

!   *** Local parameters ***

      CHARACTER(LEN=*), PARAMETER :: routineN = "mc_classical_energy"
      REAL(KIND = dp)       :: fpie=2.5669722003d-38 !C^2/Hartree/a.u.
! Locals
      INTEGER                           :: i,imol,iimol,handle,ipart,&
                                           jpart,jmol
      TYPE(particle_list_type), POINTER :: particles
      TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: oldsys
      TYPE(cp_error_type) :: error
      REAL(KIND=dp),DIMENSION(1:3)    :: RIJ,min_IJ,charges
      TYPE(cell_type), POINTER :: cell
      REAL(KIND=dp) :: dist,four_epsilon,box_length,rcut,&
                       LJ_energy,charge_energy
      REAL(KIND=dp),DIMENSION(1:3,1:3,1:mc_par%nchain) :: r
!------------------------------------------------------------------------------

! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

! initialize some stuff
      NULLIFY(oldsys,particles)
      LJ_energy=0.0d0
      charge_energy=0.0d0
      charges(1:3) = (/ mc_par%charge,-mc_par%charge/2.0d0,&
                                      -mc_par%charge/2.0d0 /)
      four_epsilon=mc_par%epsilon*4.0d0

! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

      rcut=cell%hmat(1,1)/2.0d0
      box_length=cell%hmat(1,1)

! put the coordinates into an easier matrix to manipulate
      DO imol=1,mc_par%nchain
!         write(6,*) 'BLAH',mc_par%nchain
         r(1:3,1,imol)=particles%els(imol*3-2)%r(1:3)
         r(1:3,2,imol)=particles%els(imol*3-1)%r(1:3)
         r(1:3,3,imol)=particles%els(imol*3)%r(1:3)
      ENDDO

!      DO imol=1,mc_par%nchain
!         WRITE(6,*) r(1:3,1,imol)
!         WRITE(6,*) r(1:3,2,imol)
!    WRITE(6,*) r(1:3,3,imol)
!      ENDDO
      DO imol=1,mc_par%nchain-1
         DO jmol=imol+1,mc_par%nchain
        
! find the minimum image distance between oxygens
            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*ANINT(&
               (r(1,1,imol)-r(1,1,jmol))/box_length)
            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*ANINT(&
               (r(2,1,imol)-r(2,1,jmol))/box_length)
            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*ANINT(&
               (r(3,1,imol)-r(3,1,jmol))/box_length)
      
! this is the displacement for the oxygens, which we'll use
! in the Coloumbic area
            min_IJ(1)=box_length*ANINT((r(1,1,imol)&
                -r(1,1,jmol))/box_length)
            min_IJ(2)=box_length*ANINT((r(2,1,imol)&
                -r(2,1,jmol))/box_length)
            min_IJ(3)=box_length*ANINT((r(3,1,imol)&
                -r(3,1,jmol))/box_length)

            dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))
!            write(6,*) 'dist =',dist,rcut,box_length
!            write(6,*) 'RIJ =',RIJ(1:3)
!            write(6,*) 'min_IJ =',min_IJ(1:3)
            IF ( dist < rcut  ) THEN

! LJ energy
               LJ_energy=LJ_energy+four_epsilon*(&
                  (mc_par%sigma/dist)**12-&
                  (mc_par%sigma/dist)**6)

! Coulombic energy
           
                DO ipart=1,3
                   DO jpart=1,3
                 
! find minimum image and the energy
                      RIJ(1)=r(1,ipart,imol)-r(1,jpart,jmol)-min_IJ(1)
                      RIJ(2)=r(2,ipart,imol)-r(2,jpart,jmol)-min_IJ(2)
                      RIJ(3)=r(3,ipart,imol)-r(3,jpart,jmol)-min_IJ(3)

                      dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))

                      charge_energy=charge_energy+charges(ipart)&
                           *charges(jpart)/(fpie*dist)

                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO


 
! end the timing
       CALL timestop(0.0d0,handle)
!       write(6,*) 'LJ_energy =',LJ_energy,charge_energy
! add up the total energy
       classical_energy_new=LJ_energy+charge_energy
!      STOP

END SUBROUTINE mc_classical_energy

END MODULE mc_misc

!******************************************************************************
