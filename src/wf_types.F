!!
!!
!!

MODULE wf_types
  
  USE coefficient_types,               ONLY: PLANEWAVES,&
                                             PW_COMPLEXDATA1D,&
                                             coeff_type
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_allocate
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: wf_type, wf_spinset_type, wf_kset_type
  PUBLIC :: wf_allocate, wf_state_allocate, wf_allocate_eigenvalues
  
  TYPE wf_kset_type
     TYPE ( coeff_type ), DIMENSION ( : ), POINTER :: psi_iks
     
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: eigenvalue
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: occupation
     
     INTEGER :: in_use
     INTEGER :: group_id
!???     TYPE ( particle ), POINTER :: centre
  END TYPE wf_kset_type
  
  TYPE wf_spinset_type
     TYPE ( wf_kset_type ), DIMENSION ( : ), POINTER :: psi_ks
     
     INTEGER :: in_use
     INTEGER :: group_id
!???     TYPE ( particle ), POINTER :: centre
  END TYPE wf_spinset_type
  
  TYPE wf_type
     TYPE ( wf_spinset_type ), DIMENSION ( : ), POINTER :: psi_s
     
     INTEGER :: in_use
     INTEGER :: group_id
!???     TYPE ( particle ), POINTER :: centre
  END TYPE wf_type
  
CONTAINS

!******************************************************************************

SUBROUTINE wf_allocate ( nstates, nkpts, nspins, psi, pw_grid )
  
  
    INTEGER, INTENT(IN)                      :: nstates, nkpts, nspins
    TYPE(wf_type), POINTER                   :: psi
    TYPE(pw_grid_type), INTENT(IN), &
      OPTIONAL, TARGET                       :: pw_grid

    INTEGER                                  :: allocstat, kpt, spin, state
    TYPE(coeff_type), POINTER                :: psi_iks
    TYPE(wf_kset_type), POINTER              :: psi_ks
    TYPE(wf_spinset_type), POINTER           :: psi_s

!------------------------------------------------------------------------------

  IF ( PRESENT ( pw_grid ) ) THEN
     psi % in_use = PLANEWAVES
  ELSE
     CALL stop_program ( "wf_allocate", "no possible grid input" )
  END IF
  
  ALLOCATE ( psi % psi_s ( nspins ), STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_program ( "wf_allocate", "error allocating psi_s" )
  
  DO spin = 1, nspins
     psi_s => psi % psi_s ( spin )
     psi_s % in_use = psi % in_use
     
     ALLOCATE ( psi_s % psi_ks ( nkpts ), STAT = allocstat )
     IF ( allocstat /= 0 ) &
          CALL stop_program ( "wf_allocate", "error allocating psi_ks" )
     
     DO kpt = 1, nkpts
        psi_ks => psi_s % psi_ks ( kpt )
        psi_ks % in_use = psi % in_use
        
        !! The eigenvalues will be allocated later if they become available
        NULLIFY ( psi_ks % eigenvalue )
        
        ALLOCATE ( psi_ks % occupation ( nstates ), STAT = allocstat )
        IF ( allocstat /= 0 ) &
             CALL stop_program ( "wf_allocate", "error allocating occupation" )
        
        ALLOCATE ( psi_ks % psi_iks ( nstates ), STAT = allocstat )
        IF ( allocstat /= 0 ) &
             CALL stop_program ( "wf_allocate", "error allocating psi_iks" )
        
        DO state = 1, nstates
           psi_iks => psi_ks % psi_iks ( state )
           
           psi_ks % occupation ( state ) = -1000.0_dp
           
           IF ( PRESENT ( pw_grid ) ) THEN
              psi_iks % in_use = PLANEWAVES
              
              psi_iks % pw % pw_grid => pw_grid
              
              CALL pw_allocate ( psi_iks % pw, pw_grid, &
                   use_data = PW_COMPLEXDATA1D )
           ELSE
              CALL stop_program ( "wf_allocate", "no possible grid input" )
           END IF
           
        END DO
     END DO
  END DO
  
END SUBROUTINE wf_allocate

!******************************************************************************

SUBROUTINE wf_allocate_eigenvalues ( psi )
  
  
    TYPE(wf_type), DIMENSION(:)              :: psi

    INTEGER                                  :: allocstat, kpt, nkpts, &
                                                nspins, nstates, nwfsets, &
                                                spin, wfset
    TYPE(wf_kset_type), POINTER              :: psi_ks
    TYPE(wf_spinset_type), POINTER           :: psi_s

!------------------------------------------------------------------------------

  nwfsets = SIZE ( psi )
  DO wfset = 1, nwfsets
     
     nspins = SIZE ( psi ( wfset ) % psi_s )
     DO spin = 1, nspins
        psi_s => psi ( wfset ) % psi_s ( spin )
        
        nkpts = SIZE ( psi_s % psi_ks )
        DO kpt = 1, nkpts
           psi_ks => psi_s % psi_ks ( kpt )
           
           nstates = SIZE ( psi_ks % occupation )
           
           ALLOCATE ( psi_ks % eigenvalue ( nstates ), STAT = allocstat )
           IF ( allocstat /= 0 ) &
                CALL stop_program ( "wf_allocate", "error allocating eigenvalue" )
        END DO
     END DO
  END DO
  
END SUBROUTINE wf_allocate_eigenvalues

!******************************************************************************

SUBROUTINE wf_state_allocate ( psi_iks, pw_grid )
  
  
    TYPE(coeff_type), INTENT(OUT)            :: psi_iks
    TYPE(pw_grid_type), INTENT(IN), &
      OPTIONAL, TARGET                       :: pw_grid

!------------------------------------------------------------------------------

  IF ( PRESENT ( pw_grid ) ) THEN
     psi_iks % in_use = PLANEWAVES
     
     psi_iks % pw % pw_grid => pw_grid
     
     CALL pw_allocate ( psi_iks % pw, pw_grid, use_data = PW_COMPLEXDATA1D )
  ELSE
     CALL stop_program ( "wf_state_allocate", "no grid input" )
  END IF
  
END SUBROUTINE wf_state_allocate

!******************************************************************************
  
END MODULE wf_types
