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
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_move_control,                 ONLY: mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE :: nonbonded_energy
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


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER              :: routineN = "init_mc_averages"

    INTEGER                                  :: handle, stat

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


    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_averages_type), POINTER          :: averages
    INTEGER, INTENT(IN)                      :: nnstep, nchain, data

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


    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(mc_moves_type), POINTER             :: moves
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)         :: energy_check, initial_energy, &
                                                old_energy
    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER         :: routineN = "final_mc_write"

    INTEGER                                  :: handle

!------------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

    WRITE(iw,*)
    WRITE(iw,*) '****************************************************',&
                   '***************************'
    WRITE(iw,'(A,T66,F15.8)') " Average Energy [Hartrees]:", &
                averages%ave_energy
    IF (averages%ave_energy_squared-averages%ave_energy**2 .LT. 0.0_dp) THEN
       WRITE(iw,'(A,T66,E15.8)') " Energy Fluctuation [Hartrees]: IMAGINARY PART ", &
         DSQRT(abs(averages%ave_energy_squared-averages%ave_energy**2))
    ELSE
      WRITE(iw,'(A,T66,E15.8)') " Energy Fluctuation [Hartrees]:", &
            DSQRT(averages%ave_energy_squared-averages%ave_energy**2)
    ENDIF
    ! XXXXXXXXXXXXXxx do not use hand coded constants all available (physcon, units, ...)
 !    WRITE(iw,'(A,A,T66,F15.8)') " Heat Capacity at ",&
  !         "constant volume [J/mol/K]",&
   !        DSQRT(averages%ave_energy_squared-averages%ave_energy**2)&
    !       *(4.3597482d-18)**2/(mc_par%temp)**2/1.380658d-23*6.022d23
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

! NOTE: I assume
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(OUT)               :: classical_energy_new
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_classical_energy"

    INTEGER                                  :: handle,imol,ipart,jmol,&
                                                jpart,ntij,i
    REAL(KIND=dp) :: box_length, charge_energy, dist, &
              fpie = 2.5669722003d-38 , LJ_energy, rcut,unit_energy,&
              stretch_energy,bend_energy,mag_HH,kstretch,OH_bond1,&
              OH_bond2
    REAL(KIND=dp), DIMENSION(1:3)            :: charges, min_IJ, RIJ,&
                                 kbend,HH,OH1,OH2
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_error_type)                      :: error
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:3, 1:mc_par%nchain)   :: r
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

      rcut=7.8d0/0.529177d0 ! always the same for this potential
      box_length=cell%hmat(1,1)

! put the coordinates into an easier matrix to manipulate
      DO imol=1,mc_par%nchain
         r(1:3,1,imol)=particles%els(imol*3-2)%r(1:3)
         r(1:3,2,imol)=particles%els(imol*3-1)%r(1:3)
         r(1:3,3,imol)=particles%els(imol*3)%r(1:3)
      ENDDO

! we need to calculate the energy of this configuration
! first, let's calculate intramolecular (bending and stretching)
      DO imol=1,mc_par%nchain
     
! stretching
         OH_bond1=DSQRT((r(1,1,imol)-r(1,2,imol))**2+&
                 (r(2,1,imol)-r(2,2,imol))**2+&
                   (r(3,1,imol)-r(3,2,imol))**2)
         OH_bond2=DSQRT((r(1,1,imol)-r(1,3,imol))**2+&
                   (r(2,1,imol)-r(2,3,imol))**2+&
                   (r(3,1,imol)-r(3,3,imol))**2)

         stretch_energy=stretch_energy+0.5d0*kstretch*((OH_bond1-&
            1.81200d0)**2+(OH_bond2-1.81200d0)**2)
!         write(6,*) 'bonds ',OH_bond1,OH_bond2

! bond bending
         DO i=1,3
            HH(i)=r(i,2,imol)-r(i,3,imol)
         ENDDO
         mag_HH=DSQRT(DOT_PRODUCT(HH,HH))
         bend_energy=bend_energy+0.5d0*kbend(1)*(mag_HH-3.0701d0)**2+&
             0.33333333333333333d0*kbend(2)*(mag_HH-3.0701d0)**3+&
             0.25d0*kbend(3)*(mag_HH-3.0701d0)**4
          
      ENDDO

! now let's find the LJ energy between all the oxygens and
! the charge interactions
      DO imol=1,mc_par%nchain-1
         DO jmol=imol+1,mc_par%nchain
        
! find the minimum image distance between oxygens
            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*ANINT(&
                 (r(1,1,imol)-r(1,1,jmol))/box_length)
            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*ANINT(&
                 (r(2,1,imol)-r(2,1,jmol))/box_length)
            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*ANINT(&
                 (r(3,1,imol)-r(3,1,jmol))/box_length)
      
            min_IJ(1)=box_length*ANINT((r(1,1,imol)-&
                        r(1,1,jmol))/box_length)
            min_IJ(2)=box_length*ANINT((r(2,1,imol)-&
                        r(2,1,jmol))/box_length)
            min_IJ(3)=box_length*ANINT((r(3,1,imol)-&
                        r(3,1,jmol))/box_length)

            dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))

            IF ( dist < 7.8d0 ) THEN

! LJ energy
               CALL nonbonded_energy(dist,unit_energy,1)
               LJ_energy=LJ_energy+unit_energy
!               charge_energy=charge_energy+charges(1)&
!                         *charges(1)/(fpie*dist)
               charge_energy=charge_energy+charges(1)&
                         *charges(1)/(dist)

! Coulombic energy
           
               DO ipart=1,3
                  DO jpart=1,3
                     IF ( ipart == 1 .AND. jpart == 1 ) CYCLE
! find minimum image and the energy
                     RIJ(1)=r(1,ipart,imol)-r(1,jpart,jmol)-min_IJ(1)
                     RIJ(2)=r(2,ipart,imol)-r(2,jpart,jmol)-min_IJ(2)
                     RIJ(3)=r(3,ipart,imol)-r(3,jpart,jmol)-min_IJ(3)
                     
                     dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        
!                     charge_energy=charge_energy+charges(ipart)&
!                                   *charges(jpart)/(fpie*dist)
                     charge_energy=charge_energy+charges(ipart)&
                                   *charges(jpart)/(dist)

! LJ energy
                     ntij=(ipart-1)*10 + jpart
                     CALL nonbonded_energy(dist,unit_energy,ntij)
                     LJ_energy=LJ_energy+unit_energy
                        
                  ENDDO
               ENDDO
            ENDIF
            
         ENDDO
      ENDDO

! sum up the total energy
!      charge_energy=charge_energy/1.380658d-23
!      WRITE(6,*) 'energies ',bend_energy,stretch_energy,LJ_energy,&
!                        charge_energy
      classical_energy_new=bend_energy+&
                stretch_energy+LJ_energy+charge_energy
!                STOP
 
! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE mc_classical_energy

SUBROUTINE nonbonded_energy(dist,ljsami,ntij)

! output ljsami is in Hartree, input dist is in a.u.

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ntij
      REAL(KIND=dp),INTENT(IN) :: dist
      REAL(KIND=dp),DIMENSION(2:15,1:3) :: A
      REAL(KIND=dp),DIMENSION(1:3) :: r_core,pshift
      REAL(KIND=dp),INTENT(OUT) :: ljsami
!      REAL(KIND=dp) :: rij
      INTEGER :: i,j,index
      
      ljsami=0.0d0
!      rij=dist/0.529177d0
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
      pshift(1)=0.127827d0
      pshift(2)=2.22408d0
      pshift(3)=0.256530d0
      r_core(1)=4.4d0
      r_core(2)=2.95d0
      r_core(3)=3.165d0
      A(2,1) = 4.786380682394d0
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
      
      ljsami=0.0d0
      IF (dist .GT. r_core(index)) THEN
         DO j=2,15
            ljsami=ljsami+A(j,index)/(dist**(j-1)*REAL(j-1))
         ENDDO
         ljsami=ljsami+pshift(index)
      ENDIF

END SUBROUTINE nonbonded_energy


END MODULE mc_misc

!******************************************************************************
