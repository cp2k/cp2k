!!
!!
!!

MODULE potential_externals

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: potential_external_update, potential_external_initialise
  
CONTAINS

!******************************************************************************

SUBROUTINE potential_external_update ( v_external )
  
  
    TYPE(potential_external_type), &
      INTENT(INOUT)                          :: v_external

!------------------------------------------------------------------------------

  DO potential_index = 1, npotentials
     potential => potential ( potential_index )
     ncentres => potential % ncentres
     
     !OpenMP
     DO g = 1, ngs
        vtemp ( g ) = vtemp ( g ) &
             + structure_factor ( g, potential_kind ) &
             * vps ( g, potential_kind )
     END DO
  END DO
  
END SUBROUTINE potential_external_update

!******************************************************************************

END MODULE potential_externals
