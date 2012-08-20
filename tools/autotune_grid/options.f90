MODULE option_module
  IMPLICIT NONE

  TYPE option_type
    INTEGER :: ig_loop_unroll_lxp=0
    INTEGER :: ig_loop_vector_notation=0
    INTEGER :: jg_loop_unroll_lxp=0
    INTEGER :: jg_loop_vector_notation=0
    INTEGER :: kg_loop_unroll_lxp=0
    INTEGER :: kg_loop_vector_notation=0
  END TYPE

  TYPE(option_type), DIMENSION(:), ALLOCATABLE, SAVE :: all_options

  INTEGER, PARAMETER :: DEFAULT_LMAX = 9 ! CP2K expects collocate/integrate routines for l up to 9.

CONTAINS
  !
  ! this generates all acceptable options 
  !
  SUBROUTINE generate_all_options()
    INTEGER :: ig_loop_unroll_lxp,ig_loop_vector_notation,jg_loop_unroll_lxp, &
               jg_loop_vector_notation,kg_loop_unroll_lxp,kg_loop_vector_notation
    INTEGER :: nopt,iloop,total_nopt,combinations

    IF (ALLOCATED(all_options)) RETURN
     
    total_nopt = 6;
    combinations=2**total_nopt;
    ALLOCATE(all_options(combinations))

    ! count and associate all legal options, iloop=0 : count, iloop=1 : init
    !DO iloop=0,1

       nopt=0 
       DO ig_loop_unroll_lxp=0,1
       DO ig_loop_vector_notation=0,1
       DO jg_loop_unroll_lxp=0,1
       DO jg_loop_vector_notation=0,1
       DO kg_loop_unroll_lxp=0,1
       DO kg_loop_vector_notation=0,1
          ! skip certain combinations of options, for example, avoid vector notation
          ! IF (ig_loop_vector_notation==1) CYCLE
          ! IF (jg_loop_vector_notation==1) CYCLE
          ! IF (kg_loop_vector_notation==1) CYCLE
          nopt=nopt+1
          ! Be careful! first optimisation starts at 1 not zero
          !IF (iloop==1) THEN
             all_options(nopt)%ig_loop_unroll_lxp=ig_loop_unroll_lxp
             all_options(nopt)%ig_loop_vector_notation=ig_loop_vector_notation
             all_options(nopt)%jg_loop_unroll_lxp=jg_loop_unroll_lxp
             all_options(nopt)%jg_loop_vector_notation=jg_loop_vector_notation
             all_options(nopt)%kg_loop_unroll_lxp=kg_loop_unroll_lxp
             all_options(nopt)%kg_loop_vector_notation=kg_loop_vector_notation
          !ENDIF
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO

    !   IF (iloop==0) THEN
    !      ALLOCATE(all_options(nopt))
    !   ENDIF

    !ENDDO

  END SUBROUTINE

  SUBROUTINE deallocate_all_options()
      IF (ALLOCATED(all_options)) DEALLOCATE(all_options)
  END SUBROUTINE

END MODULE
