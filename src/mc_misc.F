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
  USE kinds,                           ONLY: dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_move_control,                 ONLY: mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_water_data,final_mc_write



!******************************************************************************

CONTAINS

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

SUBROUTINE get_water_data ( force_env, angle, bond1, bond2, nnstep,&
                            nchain  )

  IMPLICIT NONE

! Arguments
      REAL(KIND = dp),INTENT(INOUT)            :: bond1,bond2,angle
      TYPE(force_env_type), POINTER :: force_env
      INTEGER, INTENT(IN)               :: nnstep,nchain

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
          bond1_loc=bond1_loc*REAL(iimol-1)/REAL(iimol) &
                  +DSQRT(tot_squared)/REAL(iimol)

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+2)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +2)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+2)%r(3))**2
          bond2_loc=bond2_loc*REAL(iimol-1)/REAL(iimol) &
                  +DSQRT(tot_squared)/REAL(iimol)

! now find the running average bond angle
          vec_1(1:3)=particles%els(imol+1)%r(1:3)-&
              particles%els(imol)%r(1:3)
          vec_2(1:3)=particles%els(imol+2)%r(1:3)-&
              particles%els(imol)%r(1:3)

          angle_loc=angle_loc*REAL(iimol-1)/REAL(iimol) + &
                   1.0d0/REAL(iimol)*&
              DACOS(DOT_PRODUCT(vec_1,vec_2)/(bond1_loc*bond2_loc))&
              *180.0d0/3.14159265358979323

       ENDDO
       bond1=bond1*REAL(nnstep-1)/REAL(nnstep) &
                      +bond1_loc/REAL(nnstep)
       bond2=bond2*REAL(nnstep-1)/REAL(nnstep) &
                  +bond2_loc/REAL(nnstep)
       angle=angle*REAL(nnstep-1)/REAL(nnstep) + &
                   angle_loc/REAL(nnstep)

       WRITE(76,*) nnstep,bond1,bond2,angle 
 
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

SUBROUTINE final_mc_write (mc_par,ave_energy,bond1,bond2,angle,moves,&
                           ave_density)

  IMPLICIT NONE

! Arguments
      REAL(KIND = dp),INTENT(IN)            :: bond1,bond2,angle,&
                                               ave_energy,ave_density
      TYPE(mc_moves_type),POINTER  :: moves
      TYPE(mc_simulation_parameters_type), POINTER :: mc_par

!   *** Local parameters ***

      CHARACTER(LEN=*), PARAMETER :: routineN = "final_mc_write"

!   *** Local variables ***
      INTEGER :: handle
 
!------------------------------------------------------------------------------

! begin the timing of the subroutine
    CALL timeset(routineN,'I','',handle)

    WRITE(6,*)
    WRITE(6,*) '***************************************************'
    WRITE(6,*) "Average Energy:                    ", ave_energy
    WRITE(6,*) "Average O-H distance 1:            ", bond1
    WRITE(6,*) "Average O-H distance 2:            ", bond2
    WRITE(6,*) "Average H-O-H angle:               ", angle
    IF (ave_density .NE. 0.0d0) THEN
       WRITE(6,*) "Average Density [particles/a.u.^3]:",&
                    "   ", ave_density
    ENDIF
    WRITE(6,*)
    WRITE(6,*) "Move Data"
    IF (moves%volume%attempts .GT. 0) THEN
       WRITE(6,*) "Volume Moves        Attempted   ",&
                 "    Accepted       Percent"
       WRITE(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%volume%attempts,&
             moves%volume%successes,&
             real(moves%volume%successes)/ &
             real(moves%volume%attempts)*100.0d0
       write(6,'(A,F10.5)') ' Maximum Volume displacement = '&
                     ,mc_par%rmvolume
    ENDIF

    IF (moves%transion%attempts .GT. 0) THEN
       write(6,*) "Transion Moves        Attempted   ",&
                 "    Accepted       Percent"
       write(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%transion%attempts,&
             moves%transion%successes,&
             real(moves%transion%successes)/ &
             real(moves%transion%attempts)*100.0d0
       write(6,'(A,F10.5)') &
            ' Maximum atomic translational displacement = '&
                  ,mc_par%rmion
    ENDIF

    IF (moves%trans%attempts .GT. 0) THEN
       write(6,*) "Translation Moves     Attempted   ",&
                 "    Accepted       Percent"
       write(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%trans%attempts,&
             moves%trans%successes,&
             real(moves%trans%successes)/ &
             real(moves%trans%attempts)*100.0d0
       write(6,'(A,F10.5)') &
             ' Maximum molecular translational displacement = '&
                  ,mc_par%rmtrans
    ENDIF

    IF (moves%rot%attempts .GT. 0) THEN
       write(6,*) "Rotation Moves        Attempted   ",&
                 "    Accepted       Percent"
       write(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%rot%attempts,&
             moves%rot%successes,&
             real(moves%rot%successes)/ &
             real(moves%rot%attempts)*100.0d0
       write(6,'(A,F10.5)') ' Maximum rotational displacement = '&
       ,mc_par%rmrot
    ENDIF
    write(6,*) '***************************************************'
    write(6,*)
 
! end the timing
    CALL timestop(0.0d0,handle)

 
END SUBROUTINE final_mc_write

!******************************************************************************

END MODULE mc_misc

!******************************************************************************
