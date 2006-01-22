!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/print_key_types [1.0] *
!!
!!   NAME
!!     print_key_types
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - print keys, basis_set_file name and potential_file_name added to the
!!       global type (27.02.2001, MK)
!!     - JGH (28.11.2001) : added pp_library_path to type
!!     - Merged with MODULE print_keys (17.01.2002, MK)
!!     - moved print_keys to print_key_types module [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE print_key_types
  USE message_passing,                 ONLY: mp_stop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='print_key_types'

  PUBLIC :: print_keys_type
  PUBLIC :: print_keys_create, print_keys_retain, print_keys_release,&
       print_keys_reset

  INTEGER, SAVE, PRIVATE :: last_print_k_id=0

! *** Print level parameters ***
  INTEGER, PARAMETER, PUBLIC :: SILENT = 0,&
                        LOW = 1,&
                        MEDIUM = 2,&
                        HIGH = 3,&
                        DEBUG = 4

!!***
!****************************************************************************

!!****s* print_key_types/print_keys_type [1.0] *
!!
!!   NAME
!!     print_keys_type
!!
!!   FUNCTION
!!     
!!   NOTES
!!     -
!!
!!   SOURCE
!!
  TYPE print_keys_type
    INTEGER :: level,precision,ref_count, id_nr
    LOGICAL :: charge_density_matrix,&
               core_hamiltonian_matrix,&
               density_matrix,&
               derivatives,&
               hartree_matrix,&
               kinetic_energy_matrix,&
               kohn_sham_matrix
    LOGICAL :: optical_conductivity,&
               ortho_matrix,&
               overlap_matrix,&
               potential_energy_matrix,&
               projectors,&
               rho0_information,&
               sab_orb_neighbor_lists,&
               sac_ppl_neighbor_lists,&
               sac_ppnl_neighbor_lists,&
               sbc_ppnl_neighbor_lists,&
               sac_oce_neighbor_lists,&
               sbc_oce_neighbor_lists,&
               scf,&
               scf_energies,&
               total_densities
    LOGICAL :: w_matrix,&
               mo_orthonormality
    LOGICAL :: ep_matrixes, ep_operator
    INTEGER    nlumos,nhomos
    LOGICAL :: loc_info
  END TYPE print_keys_type

!!*****
CONTAINS

  FUNCTION key_value(level,level_values) RESULT(value)

    INTEGER, INTENT(IN)                      :: level
    LOGICAL, DIMENSION(3), INTENT(IN)        :: level_values
    LOGICAL                                  :: value

!   ---------------------------------------------------------------------------

    SELECT CASE (level)
    CASE (:SILENT)
      value = .FALSE.
    CASE (LOW:HIGH)
      value = level_values(level)
    CASE (DEBUG:)
      value = .TRUE.
    END SELECT

  END FUNCTION key_value

! *****************************************************************************

!!****f* print_key_types/print_keys_create [1.0] *
!!
!!   NAME
!!     print_keys_create
!!
!!   FUNCTION
!!     creates a print key structure
!!     
!!   NOTES
!!     -
!!
!!   SOURCE
!!
  SUBROUTINE print_keys_create(print_keys,level)

    TYPE(print_keys_type), POINTER           :: print_keys
    INTEGER, INTENT(IN), OPTIONAL            :: level

    LOGICAL, PARAMETER                       :: f = .FALSE., t = .TRUE.

    INTEGER                                  :: stat

!   ---------------------------------------------------------------------------

    ALLOCATE(print_keys,stat=stat)
    IF (stat/=0) CALL mp_stop(stat,&
         "print_keys_create could not allocate print_keys")
    last_print_k_id=last_print_k_id+1
    print_keys%id_nr=last_print_k_id
    print_keys%ref_count=1
    print_keys%precision=0
    IF (PRESENT(level)) THEN
       print_keys%level=level
       CALL print_keys_reset(print_keys,level)
    ELSE
       print_keys%level=MEDIUM
       CALL print_keys_reset(print_keys,MEDIUM)
    END IF
  END SUBROUTINE print_keys_create

  SUBROUTINE print_keys_reset(print_keys,level)

    TYPE(print_keys_type), POINTER           :: print_keys
    INTEGER, INTENT(IN)                      :: level

    LOGICAL, PARAMETER                       :: f = .FALSE., t = .TRUE.

    print_keys%charge_density_matrix = key_value(level,(/f,f,t/))
    print_keys%core_hamiltonian_matrix = key_value(level,(/f,f,t/))
    print_keys%density_matrix = key_value(level,(/f,f,t/))
    print_keys%derivatives = key_value(level,(/f,f,f/))
    print_keys%ep_matrixes = key_value(level,(/f,f,t/))
    print_keys%ep_operator = key_value(level,(/f,f,t/))
    print_keys%hartree_matrix = key_value(level,(/f,f,t/))
    print_keys%kinetic_energy_matrix = key_value(level,(/f,f,t/))
    print_keys%kohn_sham_matrix = key_value(level,(/f,f,t/))
    print_keys%ortho_matrix = key_value(level,(/f,f,t/))
    print_keys%overlap_matrix = key_value(level,(/f,f,t/))
    print_keys%projectors = key_value(level,(/f,t,t/))
    print_keys%potential_energy_matrix = key_value(level,(/f,f,f/))
    print_keys%rho0_information = key_value(level,(/f,t,t/))
    print_keys%sab_orb_neighbor_lists = key_value(level,(/f,f,t/))
    print_keys%sac_ppl_neighbor_lists = key_value(level,(/f,f,t/))
    print_keys%sac_ppnl_neighbor_lists = key_value(level,(/f,f,t/))
    print_keys%sbc_ppnl_neighbor_lists = key_value(level,(/f,f,t/))
    print_keys%sac_oce_neighbor_lists = key_value(level,(/f,f,t/))
    print_keys%sbc_oce_neighbor_lists = key_value(level,(/f,f,t/))
    print_keys%scf = key_value(level,(/t,t,t/))
    print_keys%scf_energies = key_value(level,(/f,t,t/))
    print_keys%total_densities = key_value(level,(/f,t,t/))
    print_keys%w_matrix = key_value(level,(/f,f,f/))
    print_keys%mo_orthonormality = key_value(level,(/f,f,t/))

  END SUBROUTINE print_keys_reset

! *****************************************************************************

!!****f* print_key_types/print_keys_retain [1.0] *
!!
!!   NAME
!!     print_keys_retain
!!
!!   FUNCTION
!!     retains the print keys
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - print_keys: the print keys to retain 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE print_keys_retain(print_keys)
    TYPE(print_keys_type), POINTER           :: print_keys

    CHARACTER(len=*), PARAMETER :: routineN = 'print_keys_retain', &
      routineP = moduleN//':'//routineN

  IF (.NOT.ASSOCIATED(print_keys)) THEN
     CALL mp_stop(-1,"unassociated printkeys in print_key_types:print_keys_retain")
  END IF
  IF (.NOT.print_keys%ref_count>0) THEN
     CALL mp_stop(-1,"print_keys%ref_count<=0 in print_key_types:print_keys_retain")
  END IF
  print_keys%ref_count=print_keys%ref_count+1
END SUBROUTINE print_keys_retain
!***************************************************************************

!!****f* print_key_types/print_keys_release [1.0] *
!!
!!   NAME
!!     print_keys_release
!!
!!   FUNCTION
!!     releases the print keys
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - print_keys: the print keys to release 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE print_keys_release(print_keys)
    TYPE(print_keys_type), POINTER           :: print_keys

    CHARACTER(len=*), PARAMETER :: routineN = 'print_keys_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

  IF (ASSOCIATED(print_keys)) THEN
     IF (.NOT.print_keys%ref_count>0) THEN
        CALL mp_stop(-1,"print_keys%ref_count<=0 in print_key_types:print_keys_release")
     END IF
     print_keys%ref_count=print_keys%ref_count-1
     IF (print_keys%ref_count==0) THEN
        DEALLOCATE(print_keys,stat=stat)
        IF (stat/=0) THEN
           CALL mp_stop(-1,"error deallocating print_keys in print_key_types:print_keys_release")
        END IF
     END IF
  END IF
  NULLIFY(print_keys)
END SUBROUTINE print_keys_release
!***************************************************************************

END MODULE print_key_types
