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
  USE mc_move_control,                 ONLY: mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type

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
      TYPE(particle_list_type), POINTER :: particles
      INTEGER, INTENT(IN)               :: nnstep,nchain

! Locals
      INTEGER                           :: i,imol,iimol
      TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: oldsys
      TYPE(cp_error_type) :: error
      REAL(KIND = dp)       :: tot_squared,bond1_loc,&
                               bond2_loc
      REAL(KIND=dp),DIMENSION(1:3)    :: vec_1,vec_2
!------------------------------------------------------------------------------

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
          bond1_loc=DSQRT(tot_squared)
          bond1=bond1*REAL(nnstep-1)/REAL(nnstep) &
                      +bond1_loc/REAL(nnstep)

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+2)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +2)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+2)%r(3))**2
          bond2_loc=DSQRT(tot_squared)
          bond2=bond2*REAL(nnstep-1)/REAL(nnstep) &
                  +bond2_loc/REAL(nnstep)

! now find the running average bond angle
          vec_1(1:3)=particles%els(imol+1)%r(1:3)-&
              particles%els(imol)%r(1:3)
          vec_2(1:3)=particles%els(imol+2)%r(1:3)-&
              particles%els(imol)%r(1:3)

          angle=angle*REAL(nnstep-1)/REAL(nnstep) + 1/REAL(nnstep)*&
              DACOS(DOT_PRODUCT(vec_1,vec_2)/(bond1_loc*bond2_loc))&
              *180.0d0/3.14159265358979323
!       mag_1=DSQRT(vec_1(1)**2+vec_1(2)**2+vec_1(3)**2)
!       mag_2=DSQRT(vec_2(1)**2+vec_2(2)**2+vec_2(3)**2)
       ENDDO

       WRITE(76,*) nnstep,bond1,bond2,angle 
 

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

SUBROUTINE final_mc_write (ave_energy,bond1,bond2,angle,moves)

  IMPLICIT NONE

! Arguments
      REAL(KIND = dp),INTENT(IN)            :: bond1,bond2,angle,&
                                               ave_energy
      TYPE(mc_moves_type),POINTER  :: moves

! Locals
 
!------------------------------------------------------------------------------

    write(6,*)
    write(6,*) '***************************************************'
    write(6,*) "Average Energy:                    ", ave_energy
    write(6,*) "Average O-H distance 1:            ", bond1
    write(6,*) "Average O-H distance 2:            ", bond2
    write(6,*) "Average H-O-H angle:               ", angle
    write(6,*)
    write(6,*) "Move Data"
    IF (moves%transion%attempts .GT. 0) THEN
       write(6,*) "Transion Moves        Attempted   ",&
                 "    Accepted       Percent"
       write(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%transion%attempts,&
             moves%transion%successes,&
             real(moves%transion%successes)/ &
             real(moves%transion%attempts)*100.0d0
    ENDIF

    IF (moves%trans%attempts .GT. 0) THEN
       write(6,*) "Translation Moves        Attempted   ",&
                 "    Accepted       Percent"
       write(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%trans%attempts,&
             moves%trans%successes,&
             real(moves%trans%successes)/ &
             real(moves%trans%attempts)*100.0d0
    ENDIF

    IF (moves%rot%attempts .GT. 0) THEN
       write(6,*) "Rotation Moves        Attempted   ",&
                 "    Accepted       Percent"
       write(6,'(26X,I6,9X,I6,7X,F7.3)') &
             moves%rot%attempts,&
             moves%rot%successes,&
             real(moves%rot%successes)/ &
             real(moves%rot%attempts)*100.0d0
    ENDIF
    write(6,*) '***************************************************'
    write(6,*)
 
 
END SUBROUTINE final_mc_write

!******************************************************************************

END MODULE mc_misc

!******************************************************************************
