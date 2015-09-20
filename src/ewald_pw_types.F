!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief pw_types
!> \author CJM
! *****************************************************************************
MODULE ewald_pw_types
  USE ao_util,                         ONLY: exp_radius
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_init,          ONLY: init_input_type
  USE dg_types,                        ONLY: dg_create,&
                                             dg_release,&
                                             dg_retain,&
                                             dg_type
  USE dgs,                             ONLY: dg_pme_grid_setup
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: MPI_COMM_SELF
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_poisson_methods,              ONLY: pw_poisson_set
  USE pw_poisson_read_input,           ONLY: pw_poisson_read_parameters
  USE pw_poisson_types,                ONLY: &
       do_ewald_ewald, do_ewald_none, do_ewald_pme, do_ewald_spme, &
       pw_poisson_create, pw_poisson_parameter_type, pw_poisson_release, &
       pw_poisson_retain, pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE realspace_grid_types,            ONLY: &
       realspace_grid_desc_type, realspace_grid_input_type, &
       realspace_grid_type, rs_grid_create, rs_grid_create_descriptor, &
       rs_grid_print, rs_grid_release, rs_grid_release_descriptor, &
       rs_grid_retain_descriptor
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_pw_types'
  INTEGER, PRIVATE, SAVE :: last_ewald_pw_id_nr=0
  PUBLIC :: ewald_pw_type, ewald_pw_release, &
            ewald_pw_retain, ewald_pw_create, &
            ewald_pw_get, ewald_pw_set

! *****************************************************************************
  TYPE ewald_pw_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE ( pw_pool_type ), POINTER       :: pw_small_pool
     TYPE ( pw_pool_type ), POINTER       :: pw_big_pool
     TYPE ( realspace_grid_desc_type ), POINTER    :: rs_desc
     TYPE ( pw_poisson_type ), POINTER    :: poisson_env
     TYPE ( dg_type ), POINTER            :: dg
  END TYPE ewald_pw_type

CONTAINS

! *****************************************************************************
!> \brief retains the structure ewald_pw_type
!> \param ewald_pw ...
! *****************************************************************************
  SUBROUTINE ewald_pw_retain(ewald_pw)
    TYPE(ewald_pw_type), POINTER             :: ewald_pw

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(ewald_pw))
    CPASSERT(ewald_pw%ref_count>0)
    ewald_pw%ref_count=ewald_pw%ref_count+1
  END SUBROUTINE ewald_pw_retain

! *****************************************************************************
!> \brief creates the structure ewald_pw_type
!> \param ewald_pw ...
!> \param ewald_env ...
!> \param cell ...
!> \param cell_ref ...
!> \param print_section ...
! *****************************************************************************
  SUBROUTINE ewald_pw_create ( ewald_pw, ewald_env, cell, cell_ref, print_section)
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: print_section

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(dg_type), POINTER                   :: dg

    NULLIFY ( dg )
    failure=.FALSE.
    ALLOCATE ( ewald_pw, stat=stat )
    CPASSERT(stat==0)
    NULLIFY ( ewald_pw % pw_big_pool )
    NULLIFY ( ewald_pw % pw_small_pool )
    NULLIFY ( ewald_pw % rs_desc )
    NULLIFY ( ewald_pw % poisson_env )
    CALL dg_create ( dg)
    ewald_pw % dg => dg
    ewald_pw%ref_count=1
    last_ewald_pw_id_nr=last_ewald_pw_id_nr+1
    ewald_pw%id_nr=last_ewald_pw_id_nr
    CALL ewald_pw_init ( ewald_pw, ewald_env, cell, cell_ref, print_section)
  END SUBROUTINE ewald_pw_create

!****f* ewald_pw_types/ewald_pw_release [1.0] *

! *****************************************************************************
!> \brief releases the memory used by the ewald_pw
!> \param ewald_pw ...
! *****************************************************************************
  SUBROUTINE ewald_pw_release ( ewald_pw)
    TYPE(ewald_pw_type), POINTER             :: ewald_pw

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    IF (ASSOCIATED(ewald_pw)) THEN
       CPASSERT(ewald_pw%ref_count>0)
       ewald_pw%ref_count=ewald_pw%ref_count-1
       IF (ewald_pw%ref_count<1) THEN
          CALL pw_pool_release ( ewald_pw % pw_small_pool)
          CALL pw_pool_release ( ewald_pw % pw_big_pool)
          CALL rs_grid_release_descriptor ( ewald_pw % rs_desc)
          CALL pw_poisson_release ( ewald_pw % poisson_env)
          CALL dg_release ( ewald_pw % dg)
          DEALLOCATE (  ewald_pw , stat = stat )
          CPASSERT(stat==0)
       END IF
    END IF
    NULLIFY ( ewald_pw )
    CALL timestop(handle)
  END SUBROUTINE ewald_pw_release

!****** ewald_pw_types/ewald_pw_init [1.0] *

! *****************************************************************************
!> \brief ...
!> \param ewald_pw ...
!> \param ewald_env ...
!> \param cell ...
!> \param cell_ref ...
!> \param print_section ...
!> \par History
!>      JGH (12-Jan-2001): Added SPME part
!>      JGH (15-Mar-2001): Work newly distributed between initialize, setup,
!>                         and force routine
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_pw_init ( ewald_pw, ewald_env, cell, cell_ref, print_section)
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: print_section

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bo(2,3), ewald_type, gmax(3), &
                                                handle, npts_s( 3 ), ns_max, &
                                                o_spline, output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, alphasq, &
                                                cutoff_radius, epsilon, norm
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_grid_type), POINTER              :: pw_big_grid, pw_small_grid
    TYPE(pw_poisson_parameter_type)          :: poisson_params
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_input_type)          :: input_settings
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(section_vals_type), POINTER         :: poisson_section, &
                                                rs_grid_section

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY ( pw_big_grid )
    NULLIFY ( pw_small_grid, poisson_section )

    CPASSERT(ASSOCIATED ( ewald_pw ))
    CPASSERT(ASSOCIATED ( ewald_env ))
    CPASSERT(ASSOCIATED ( cell ))
    CPASSERT(ewald_pw%ref_count>0)
    CALL ewald_env_get ( ewald_env = ewald_env,&
                         para_env = para_env, &
                         gmax = gmax, alpha = alpha,&
                         ns_max = ns_max,  &
                         ewald_type = ewald_type, &
                         o_spline=o_spline,&
                         poisson_section=poisson_section,&
                         epsilon = epsilon)

    rs_grid_section => section_vals_get_subs_vals(poisson_section,"EWALD%RS_GRID")

    SELECT CASE ( ewald_type )
    CASE ( do_ewald_ewald )
       ! set up Classic EWALD sum
       logger => cp_get_default_logger()
       output_unit = cp_print_key_unit_nr(logger,print_section,"",extension=".Log")
       CALL pw_grid_create( pw_big_grid, MPI_COMM_SELF)

       IF ( ANY(gmax == 2 * ( gmax / 2 )) ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,"gmax has to be odd.")
       END IF
       bo ( 1, : ) = -gmax / 2
       bo ( 2, : ) = +gmax / 2
       CALL pw_grid_setup ( cell_ref%hmat, pw_big_grid, grid_span=HALFSPACE, bounds=bo, spherical=.TRUE., &
                             fft_usage=.FALSE., iounit=output_unit)
       NULLIFY ( pw_pool )
       CALL pw_pool_create ( pw_pool, pw_grid = pw_big_grid)
       ewald_pw % pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool)
       CALL pw_pool_release ( pw_pool)
       CALL pw_grid_release ( pw_big_grid)
       CALL cp_print_key_finished_output(output_unit,logger,print_section,"")

    CASE ( do_ewald_pme )
       ! set up Particle-Mesh EWALD sum
       logger => cp_get_default_logger()
       output_unit = cp_print_key_unit_nr(logger,print_section,"",extension=".Log")
       IF (.NOT.ASSOCIATED(ewald_pw%poisson_env)) THEN
          CALL pw_poisson_create(ewald_pw%poisson_env)
       END IF
       CALL pw_grid_create( pw_small_grid, MPI_COMM_SELF)
       CALL pw_grid_create( pw_big_grid, para_env%group)
       IF ( ns_max == 2 * ( ns_max / 2 ) ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,"ns_max has to be odd.")
       END IF
       npts_s ( : ) = ns_max
       ! compute cut-off radius
       alphasq = alpha ** 2
       norm = ( 2.0_dp * alphasq / pi ) ** ( 1.5_dp )
       cutoff_radius = exp_radius ( 0, 2.0_dp * alphasq, epsilon , norm )

       CALL dg_pme_grid_setup ( cell_ref%hmat, npts_s, cutoff_radius,   &
            pw_small_grid, pw_big_grid, rs_dims=(/para_env%num_pe,1/),&
            iounit=output_unit, fft_usage=.TRUE.)
       ! Write some useful info
       IF ( output_unit > 0 ) THEN
          WRITE ( output_unit, '( A,T71,E10.4 )' ) &
               ' EWALD| Gaussian tolerance (effective) ', epsilon
          WRITE ( output_unit, '( A,T63,3I6 )' ) &
               ' EWALD| Small box grid ', pw_small_grid % npts
          WRITE ( output_unit, '( A,T63,3I6 )' ) &
               ' EWALD| Full box grid ', pw_big_grid % npts
       END IF

       ! pw pools initialized
       NULLIFY ( pw_pool )
       CALL pw_pool_create( pw_pool, pw_grid = pw_big_grid)
       ewald_pw % pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool)
       CALL pw_pool_release ( pw_pool)

       NULLIFY ( pw_pool )
       CALL pw_pool_create( pw_pool, pw_grid = pw_small_grid)
       ewald_pw%pw_small_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_small_pool)
       CALL pw_pool_release ( pw_pool)

       NULLIFY ( rs_desc )
       CALL init_input_type(input_settings,nsmax=MAXVAL ( pw_small_grid % npts(1:3) ),&
                            rs_grid_section=rs_grid_section,ilevel=1, &
                            higher_grid_layout=(/-1,-1,-1/))
       CALL rs_grid_create_descriptor(rs_desc,pw_big_grid, input_settings)

       CALL rs_grid_create(rs, rs_desc)
       CALL rs_grid_print(rs,output_unit)
       CALL rs_grid_release(rs)

       CALL cp_print_key_finished_output(output_unit,logger,print_section,"")

       ewald_pw%rs_desc => rs_desc

       CALL rs_grid_retain_descriptor ( ewald_pw % rs_desc)
       CALL rs_grid_release_descriptor ( rs_desc)

       CALL pw_grid_release ( pw_small_grid)
       CALL pw_grid_release ( pw_big_grid)

    CASE ( do_ewald_spme )
       ! set up the Smooth-Particle-Mesh EWALD sum
       logger => cp_get_default_logger()
       output_unit = cp_print_key_unit_nr(logger,print_section,"",extension=".Log")
       IF (.NOT.ASSOCIATED(ewald_pw%poisson_env)) THEN
          CALL pw_poisson_create(ewald_pw%poisson_env)
       END IF
       CALL pw_grid_create( pw_big_grid, para_env%group)
       npts_s = gmax
       CALL pw_grid_setup ( cell_ref%hmat, pw_big_grid, grid_span=HALFSPACE, npts=npts_s, spherical=.TRUE.,&
            rs_dims=(/para_env%num_pe,1/), iounit=output_unit, fft_usage=.TRUE.)

       ! pw pools initialized
       NULLIFY ( pw_pool )
       CALL pw_pool_create(pw_pool, pw_grid=pw_big_grid)
       ewald_pw%pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool)
       CALL pw_pool_release ( pw_pool)

       NULLIFY ( rs_desc )
       CALL init_input_type(input_settings,nsmax=o_spline,&
                            rs_grid_section=rs_grid_section,ilevel=1, &
                            higher_grid_layout=(/-1,-1,-1/))
       CALL rs_grid_create_descriptor(rs_desc, pw_big_grid, input_settings)

       CALL rs_grid_create(rs, rs_desc)
       CALL rs_grid_print(rs,output_unit)
       CALL rs_grid_release(rs)
       CALL cp_print_key_finished_output(output_unit,logger,print_section,"")

       ewald_pw%rs_desc => rs_desc

       CALL rs_grid_retain_descriptor ( ewald_pw % rs_desc)
       CALL rs_grid_release_descriptor ( rs_desc)

       CALL pw_grid_release ( pw_big_grid)
    CASE ( do_ewald_none )
       ! No EWALD sums..
    CASE default
       CPASSERT(.FALSE.)
    END SELECT
    ! Poisson Environment
    IF (ASSOCIATED(ewald_pw%poisson_env)) THEN
       ALLOCATE(pw_pools(1),stat=stat)
       CPASSERT(stat==0)
       pw_pools(1)%pool => ewald_pw%pw_big_pool
       CALL pw_poisson_read_parameters(poisson_section, poisson_params)
       CALL pw_poisson_set(ewald_pw%poisson_env,cell_hmat=cell%hmat,parameters=poisson_params,&
            use_level=1,pw_pools=pw_pools)
       DEALLOCATE(pw_pools,stat=stat)
       CPASSERT(stat==0)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ewald_pw_init

!****** ewald_pw_types/ewald_pw_get [1.0] *

! *****************************************************************************
!> \brief get the ewald_pw environment to the correct program.
!> \param ewald_pw ...
!> \param pw_big_pool ...
!> \param pw_small_pool ...
!> \param rs_desc ...
!> \param poisson_env ...
!> \param dg ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_pw_get (ewald_pw, pw_big_pool, pw_small_pool, rs_desc, poisson_env, dg )

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_big_pool, pw_small_pool
    TYPE(realspace_grid_desc_type), &
      OPTIONAL, POINTER                      :: rs_desc
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env
    TYPE(dg_type), OPTIONAL, POINTER         :: dg

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_pw_get', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( poisson_env ) ) poisson_env => ewald_pw % poisson_env
    IF ( PRESENT ( pw_big_pool ) )  pw_big_pool => ewald_pw % pw_big_pool
    IF ( PRESENT ( pw_small_pool ) )  pw_small_pool => ewald_pw % pw_small_pool
    IF ( PRESENT ( rs_desc ) )  rs_desc => ewald_pw % rs_desc
    IF ( PRESENT ( dg ) )  dg => ewald_pw % dg

  END SUBROUTINE ewald_pw_get

!****** ewald_pw_types/ewald_pw_set [1.0] *

! *****************************************************************************
!> \brief set the ewald_pw environment to the correct program.
!> \param ewald_pw ...
!> \param pw_big_pool ...
!> \param pw_small_pool ...
!> \param rs_desc ...
!> \param dg ...
!> \param poisson_env ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_pw_set (ewald_pw, pw_big_pool, pw_small_pool, rs_desc, dg,&
       poisson_env)

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_big_pool, pw_small_pool
    TYPE(realspace_grid_desc_type), &
      OPTIONAL, POINTER                      :: rs_desc
    TYPE(dg_type), OPTIONAL, POINTER         :: dg
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_pw_set', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( pw_big_pool ) )  THEN
        CALL pw_pool_retain ( pw_big_pool)
        CALL pw_pool_release ( ewald_pw % pw_big_pool)
        ewald_pw % pw_big_pool => pw_big_pool
    ENDIF
    IF ( PRESENT ( pw_small_pool ) )  THEN
        CALL pw_pool_retain ( pw_small_pool)
        CALL pw_pool_release ( ewald_pw % pw_small_pool)
        ewald_pw % pw_small_pool => pw_small_pool
    ENDIF
    IF ( PRESENT ( rs_desc ) )  THEN
        CALL rs_grid_retain_descriptor (rs_desc)
        CALL rs_grid_release_descriptor (ewald_pw % rs_desc)
        ewald_pw % rs_desc => rs_desc
    ENDIF
    IF ( PRESENT ( dg ) )  THEN
        CALL dg_retain ( dg)
        CALL dg_release ( ewald_pw % dg)
        ewald_pw % dg => dg
    ENDIF
    IF (PRESENT(poisson_env)) THEN
       IF (ASSOCIATED(poisson_env)) &
            CALL pw_poisson_retain(poisson_env)
       CALL pw_poisson_release(ewald_pw%poisson_env)
       ewald_pw%poisson_env => poisson_env
    END IF

  END SUBROUTINE ewald_pw_set

END MODULE ewald_pw_types
