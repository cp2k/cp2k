!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_gpw_pw_env_methods
!!
!!   NAME
!!     kg_gpw_pw_env_methods
!!
!!   FUNCTION
!!     Module used in a KF_GPW calculation to prepare the molecular
!!     densities used to calculate the KE correction molecule by molecule.
!!     Construction of cell  and pw_env for each molecule kind.
!!     Collocation of the density of each single molecule on the grid.  
!!     
!!   AUTHOR
!!     MI 
!!
!!   MODIFICATION HISTORY
!!     Created (20.12.2004)
!!
!!   SOURCE
!******************************************************************************

 
MODULE kg_gpw_pw_env_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_deallocate
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pool_create,&
                                             rs_pools_dealloc
  USE cube_utils,                      ONLY: destroy_cube_info,&
                                             init_cube_info,&
                                             cube_info_type
  USE fft_tools,                       ONLY: FFT_RADIX_NEXT,&
                                             FFT_RADIX_NEXT_ODD,&
                                             fft_radix_operations
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gaussian_gridlevel,&
                                             init_gaussian_gridlevel
  USE green_methods,                   ONLY: pw_green_fn_rebuild,&
                                             pw_green_setup
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_gpw_collocate_den,            ONLY: calculate_rho_mol
  USE kg_gpw_fm_mol_types,             ONLY: get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             fm_mol_blocks_type,&
                                             kg_fm_mol_set_type
  USE kg_gpw_pw_env_types,             ONLY: get_molbox_env,&
                                             get_rho_mol_block,&
                                             kg_molbox_env_type,&
                                             kg_sub_pw_env_type,&
                                             rho_mol_blocks_create,&
                                             rho_mol_blocks_type,&
                                             set_molbox_env,&
                                             set_rho_mol_block
  USE kinds,                           ONLY: dp, dp_size, int_size
  USE l_utils,                         ONLY: init_l_info,&
                                             l_info_retain,&
                                             l_info_release,&
                                             l_info_type
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_methods,                  ONLY: pw_env_create
  USE pw_env_types,                    ONLY: pw_env_get, &
                                             pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             PW_MODE_LOCAL,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: create_gvectors,& 
                                             pw_grid_construct,&
                                             pw_grid_release
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_dealloc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type,&
                                             init_cell,&
                                             pbc,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE green_types,                     ONLY: greens_function_type,&
                                             pw_green_release,&
                                             pw_green_retain,&
                                             pw_green_create
  USE input_section_types,             ONLY: section_vals_type  
  IMPLICIT NONE

  PRIVATE

  INTEGER :: grid_tag = 0
  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_pw_env_methods"
  INTEGER, PRIVATE, SAVE :: last_rho_id_nr=0

! *** Public subroutines ***

  PUBLIC :: build_molbox_env, kg_rho_update_rho_mol, update_rho_mol_blocks

! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****f* kg_gpw_pw_env_methods/build_molbox_env
!!
!!    NAME
!!      build_molbox_env
!!
!!    FUNCTION
!!      Build the pw_env and cell box for each molecule kind
!!      The internal coordinates of each molecule in the new box
!!      are such that the molecule is centered in the middle of the box
!!
!!    NOTE
!!      - Some shift should be added in the internal coordinates in order
!!        with respect to the grid nodes
!!
!!    INPUTS
!!
!!    AUTHOR
!!       MI
!!
!!    MODIFICATION HISTORY
!!        Creation (MI)
!!
!!*** *************************************************************************
 
  SUBROUTINE build_molbox_env(molbox_env,fm_mol,qs_env,poisson_section,error)

    TYPE(kg_molbox_env_type) , TARGET          :: molbox_env
    TYPE(kg_fm_mol_set_type), POINTER          :: fm_mol
    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error
    TYPE(section_vals_type), POINTER, OPTIONAL :: poisson_section


    CHARACTER(len=*), PARAMETER :: routineN = 'build_molbox_env',&
      routineP = module_name//':'//routineN

    INTEGER :: handle, i, iat, iatom, imol, istat, jat, jatom, nat_mol, &
               ngrid_level, nmol_local, num_pe, output_unit
    INTEGER, DIMENSION(:), POINTER           :: index_atom
    LOGICAL :: ionode, failure
    REAL(dp) :: length, Li, ra(3), rab, rab_pbc(3), rab_max, rad_at, rad_max, rb(3), side(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff

    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell, cell_ref, cell_mol
    TYPE(greens_function_type), POINTER      :: green
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: basis_set
 
    TYPE(fm_mol_blocks_type), DIMENSION(:),&
      POINTER                                :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: fm_mol_iblock
   TYPE(kg_molbox_env_type), POINTER          :: molbox_env_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env,new_pw_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER                                :: rho_mol_blocks
    REAL(KIND=dp), DIMENSION(3)              :: abc

    CALL timeset(routineN,"I","",handle)

    NULLIFY(cell,cell_ref,cell_mol,dft_control,para_env,particle_set,pw_env)
    NULLIFY(atom_kind,basis_set,fm_mol_blocks,green)

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_unit_nr(logger)
    ionode= (logger%para_env%mepos == logger%para_env%source)

    CALL get_qs_env(qs_env=qs_env, cell=cell, cell_ref=cell_ref,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    particle_set=particle_set,&
                    pw_env=pw_env)

    CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol,&
                           nmolecule_local=nmol_local,&
                           natom=nat_mol,fm_mol_blocks=fm_mol_blocks)

    molbox_env%nmolecule_local = nmol_local
    molbox_env%natom = nat_mol
    rab_max = 0.0_dp
    rad_max = 0.0_dp
    DO imol = 1,nmol_local
      NULLIFY(fm_mol_iblock,index_atom)
      fm_mol_iblock => fm_mol_blocks(imol)
      CALL get_fm_mol_block(fm_mol_block=fm_mol_iblock,& 
                            index_atom=index_atom)
      DO iat = 1,nat_mol
        iatom = index_atom(iat)
        ra(1:3) = particle_set(iatom)%r(1:3)

        atom_kind => particle_set(iatom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=basis_set)
        CALL get_gto_basis_set(gto_basis_set=basis_set,kind_radius=rad_at)
        rad_max = MAX(rad_max,rad_at)

        DO jat = iat+1,nat_mol
          jatom = index_atom(jat)
          rb(1:3) = particle_set(jatom)%r(1:3)
          rab_pbc(:) = pbc(ra(:),rb(:),cell)
          rab = SQRT((rab_pbc(1))*(rab_pbc(1))+&
                     (rab_pbc(2))*(rab_pbc(2))+&
                     (rab_pbc(3))*(rab_pbc(3)))
          rab_max = MAX(rab_max,rab)
        END DO  ! jat
      END DO  ! iat
    END DO

    CALL mp_max(rab_max,para_env%group)

    molbox_env%rab_max = rab_max
    molbox_env%rad_max = rad_max
    length = rab_max + 1.5_dp*rad_max
 
    DO i = 1,3
      Li =  SQRT ( SUM ( cell_ref % hmat ( :, i ) ** 2 ) )
      IF(Li<length) THEN
        IF (ionode) THEN
          WRITE (UNIT=output_unit,&
                FMT="(/,T2,A,/,T6,A,I3,A,f10.5,A,f10.5)")&
          "WARNING: The Cell Box Side is Probably too Small",&
          "L",i,"(= ",Li,") < max L of 1 molecule (=",length,")"
        END IF
        length = Li
      END IF 
!        WRITE(*,*) 'side ', i , length, Li
    END DO
    side(1:3) = length

    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)

    cutoff => dft_control%qs_control%e_cutoff
    ngrid_level = SIZE(cutoff)

    NULLIFY(cell_mol,new_pw_mol)

! create the cell structure for the molecule
    ALLOCATE (cell_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"cell",0)
    cell_mol%unit_of_length_name = "BOHR"
    cell_mol%deth = 0.0_dp
    cell_mol%unit_of_length = 1.0_dp
    cell_mol%orthorhombic = .TRUE.
    cell_mol%scaled_coordinates = .FALSE.
    cell_mol%subcells = 2.0_dp
    cell_mol%perd(:) = 0
    cell_mol%hmat(:,:) = 0.0_dp
    cell_mol%h_inv(:,:) = 0.0_dp

    num_pe = para_env%num_pe
    molbox_env_p => molbox_env
    CALL get_molbox_env(molbox_env=molbox_env_p,pw_env_mol=new_pw_mol)
    IF (.NOT.ASSOCIATED(new_pw_mol)) THEN
       CALL pw_env_create(new_pw_mol,error=error)
       CALL pw_env_mol_rebuild(new_pw_mol,cell_mol,pw_env,dft_control,side,num_pe,error=error)
       CALL set_molbox_env(molbox_env=molbox_env_p,pw_env_mol=new_pw_mol)
       CALL pw_env_release(new_pw_mol,error=error)
    ELSE
       CALL pw_env_mol_rebuild(new_pw_mol,cell_mol,pw_env,dft_control,side,num_pe,error=error)
    END IF

    CALL set_molbox_env(molbox_env=molbox_env_p,cell_mol=cell_mol)
    DO i = 1,3
      Li =  SQRT ( SUM ( cell_ref % hmat ( :, i ) ** 2 ) )
      length = cell_mol%hmat( i,i )
      IF(Li<length) THEN
        IF (ionode) THEN
          WRITE (UNIT=output_unit,&
                FMT="(/,T2,A,/,T6,A,I3,A,f10.5,A,f10.5)")&
          "WARNING: The Cell Box Side is Probably too Small",&
          "L",i,"(= ",Li,") < max L of 1 molecule (=",length,")"
        END IF
      END IF
    END DO    
    CALL cell_release(cell_mol,error=error)

   ! init auxbas_grid

    CALL get_molbox_env(molbox_env=molbox_env_p,pw_env_mol=new_pw_mol)
    DO i=1,ngrid_level
      IF (cutoff(i) == dft_control%qs_control%cutoff) new_pw_mol%auxbas_grid=i
    END DO

   ! complete init of green's function
    CALL get_cell(cell_mol, abc=abc)

    CALL pw_green_setup(green, poisson_section, cell_mol%perd, abc, error )

    CALL pw_green_fn_rebuild(green,&
            new_pw_mol%pw_pools(new_pw_mol%auxbas_grid)%pool)

    CALL set_molbox_env(molbox_env=molbox_env_p,green=green)
    CALL pw_green_release(green, error=error)

   ! And now initialize the rho sructure for each molecule
    NULLIFY(rho_mol_blocks)
    CALL rho_mol_blocks_create(rho_mol_blocks,nmol_local,nat_mol,fm_mol_blocks,error=error)
    CALL set_molbox_env(molbox_env=molbox_env_p,rho_mol_blocks=rho_mol_blocks)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_molbox_env

! *****************************************************************************
!!****f* kg_gpw_pw_env_methods/pw_env_mol_rebuild
!!
!!    NAME
!!       pw_env_mol_rebuild
!!
!!    FUNCTION
!!      Build the pw_env and cell box for each molecule kind
!!      Multigrids can be used
!!
!!    NOTE
!!      - The densest grid of pw_env_mol has the same dr of the densest full grid
!!        The box size is changed to match Npt*dr
!!
!!    INPUTS
!!
!!    AUTHOR
!!       MI
!!
!!    MODIFICATION HISTORY
!!        Creation (MI)
!!
!!*** *************************************************************************
 

  SUBROUTINE pw_env_mol_rebuild(pw_mol,cell_mol,pw_global,dft_control,side,num_pe,error)

    TYPE(pw_env_type), POINTER                 :: pw_global,pw_mol
    TYPE(cell_type), POINTER                   :: cell_mol
    TYPE(dft_control_type), POINTER            :: dft_control
    REAL(dp), INTENT(INOUT)                    :: side(3)
    INTEGER, INTENT(IN)                        :: num_pe
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_mol_rebuild',&
      routineP = module_name//':'//routineN

    INTEGER :: auxbas_grid, cmax, cmaxl, i, igrid_level, istat, maxlgto, &
               my_icommensurate, my_ncommensurate, ncommensurate, &
               ngrid_level,  nlowest, nlowest_new, npts_out(3), npts_tmp(3), nsmax
    LOGICAL :: blocking, failure, fft, symmetry
    REAL(dp) :: cutilev, dr_g(3), drmin, ecut, hmat(3,3), maxradius, rel_cutoff
    REAL(dp), DIMENSION(:), POINTER            :: cutoff

    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                  :: rs_pools_mol
    TYPE(cube_info_type), DIMENSION(:), &
       POINTER                                 :: cube_info_g
    TYPE(l_info_type), POINTER                 :: l_info_g
    TYPE(pw_grid_type), POINTER                :: pw_grid_g
    TYPE(pw_grid_type), POINTER                :: pw_grid_mol
    TYPE(pw_grid_type), POINTER                :: ref_pw_grid_mol
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                  :: pw_pools_mol
    TYPE(pw_pool_type), POINTER                :: pw_pool_g


    failure = .FALSE.

    NULLIFY(cube_info_g,l_info_g,pw_pools_mol,rs_pools_mol)

    CPPrecondition(ASSOCIATED(pw_mol),cp_failure_level,routineP,error,failure)
    CPPrecondition(pw_mol%ref_count>0,cp_failure_level,routineP,error,failure)

    CPPrecondition(ASSOCIATED(pw_global),cp_failure_level,routineP,error,failure)
    CPPrecondition(pw_global%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      CALL pw_pools_dealloc(pw_mol%pw_pools,error=error)
      CALL rs_pools_dealloc(pw_mol%rs_pools,error=error)
      IF (ASSOCIATED(pw_mol%gridlevel_info)) THEN
        CALL destroy_gaussian_gridlevel(pw_mol%gridlevel_info,error=error)
      ELSE
        ALLOCATE(pw_mol%gridlevel_info, STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF (ASSOCIATED(pw_mol%l_info)) THEN
        CALL l_info_release(pw_mol%l_info,error)
      ELSE
        ALLOCATE(pw_mol%l_info,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF(dft_control%qs_control%gapw) THEN
        CALL stop_program(routineP,"KG_GPW + GAPW not valid option")
      ELSEIF(dft_control%qs_control%se_control%ewald) THEN
        CALL stop_program(routineP,"KG_GPW + EWALD not valid option")
      ELSE
        NULLIFY(pw_mol%l_info_rho0)
      END IF

      IF (ASSOCIATED(pw_mol%cube_info)) THEN
        DO igrid_level=1,SIZE(pw_mol%cube_info)
           CALL destroy_cube_info(pw_mol%cube_info(igrid_level))
        END DO
        DEALLOCATE(pw_mol%cube_info,STAT=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END IF
      NULLIFY(pw_mol%pw_pools, pw_mol%cube_info)
    END IF

    IF(.NOT. failure) THEN
      rel_cutoff = dft_control%qs_control%relative_cutoff
      cutoff => dft_control%qs_control%e_cutoff
      ngrid_level = SIZE(cutoff)

     ! init gridlevel_info XXXXXXXXX setup mapping to the effective cutoff ?
     !                     XXXXXXXXX the cutoff array here is more a 'wish-list'
     !                     XXXXXXXXX same holds for radius
     CALL init_gaussian_gridlevel(pw_mol%gridlevel_info,&
          ngrid_levels=ngrid_level,cutoff=cutoff,rel_cutoff=rel_cutoff)

      ! init pw_grids and pools
      ALLOCATE(pw_pools_mol(ngrid_level),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(rs_pools_mol(ngrid_level),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      IF (dft_control % qs_control % commensurate_mgrids) THEN
         ncommensurate=ngrid_level 
         my_icommensurate = 1
         my_ncommensurate = ncommensurate
      ELSE 
         ncommensurate=0
         my_icommensurate = 0
         my_ncommensurate = ncommensurate
      ENDIF

      CALL pw_env_get(pw_env=pw_global,l_info=l_info_g,cube_info=cube_info_g,&
                      auxbas_pw_pool=pw_pool_g,auxbas_grid=auxbas_grid)
      IF(auxbas_grid/=1) CALL stop_program ( "pw_env_mol_rebuild", &
                          "the grid with the largest cutoff is not the first" )

     ! some logical variables
      fft = .TRUE.
      symmetry = .TRUE.
      blocking = .FALSE.

!     igrid_level == 1
       ALLOCATE(pw_grid_mol,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL pw_grid_construct(pw_grid_mol)

       pw_grid_mol%para%rs_dims ( 1 ) = num_pe
       pw_grid_mol%para%rs_dims ( 2 ) = 1
       IF ( dft_control % qs_control % pw_grid_opt % spherical ) THEN
           cutilev = cutoff(1)
           pw_grid_mol%grid_span = HALFSPACE
       ELSE IF ( dft_control % qs_control % pw_grid_opt % fullspace ) THEN
           cutilev = -cutoff(1)
           pw_grid_mol%grid_span = FULLSPACE
       ELSE
           cutilev = -cutoff(1)
           pw_grid_mol%grid_span = HALFSPACE
       END IF 

       pw_grid_g => pw_pool_g%pw_grid

       ! assign a unique tag to this grid
       grid_tag = grid_tag + 1
       pw_grid_mol % identifier = grid_tag
       pw_grid_mol% para % mode = PW_MODE_LOCAL
       pw_grid_mol% para % group_head = .TRUE.

       pw_grid_mol % cutoff = ABS ( cutilev )
       IF ( SUM ( ABS ( pw_grid_mol % bounds ( :, : ) ) ) == 0 ) THEN

         dr_g(1:3) = pw_grid_g % dr(1:3)
         DO i = 1,3
           npts_tmp(i) = FLOOR(side(i)/dr_g(i))
           CALL fft_radix_operations ( npts_tmp(i), npts_out(i), FFT_RADIX_NEXT )
           side(i) = npts_out(i)*dr_g(i)
           pw_grid_mol%bounds ( 1, i ) = -npts_out(i)/2
           pw_grid_mol%bounds ( 2, i ) = pw_grid_mol%bounds ( 1, i ) + npts_out(i) - 1
         END DO

       END IF

       IF ( cutilev < 0.0_dp ) THEN
          pw_grid_mol% spherical = .FALSE.
          ecut = 1.e10_dp
       ELSE
          pw_grid_mol% spherical = .TRUE.
          ecut = cutilev
       END IF
       
       IF ( .NOT. pw_grid_mol % spherical ) THEN

         IF ( SUM ( ABS ( pw_grid_mol % bounds ( :, : ) ) ) == 0 ) THEN
            CALL stop_program ( "grid_setup", &
                                 "provide initial values for bounds" )
         END IF

         npts_tmp(:) = pw_grid_mol% bounds ( 2, : ) - pw_grid_mol% bounds ( 1, : ) + 1

         IF ( pw_grid_mol% grid_span == HALFSPACE .AND. symmetry ) THEN

           CALL fft_radix_operations ( npts_tmp(1), npts_out(1), FFT_RADIX_NEXT_ODD )
           CALL fft_radix_operations ( npts_tmp(2), npts_out(2), FFT_RADIX_NEXT_ODD )
           CALL fft_radix_operations ( npts_tmp(3), npts_out(3), FFT_RADIX_NEXT_ODD )

         ELSE

          ! keep looping to find the right one
           DO
             CALL fft_radix_operations ( npts_tmp(1), npts_out(1), FFT_RADIX_NEXT )
             ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
             nlowest=npts_out(1)/2**(my_ncommensurate-my_icommensurate)
             CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )

             IF (nlowest==nlowest_new .AND. MODULO(npts_out(1),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                EXIT
             ELSE
                npts_tmp(1)=npts_out(1)+1
             ENDIF
           ENDDO

           ! keep looping to find the right one
           DO
             CALL fft_radix_operations ( npts_tmp(2), npts_out(2), FFT_RADIX_NEXT )
             ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
             nlowest=npts_out(2)/2**(My_ncommensurate-my_icommensurate)
             CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
             IF (nlowest==nlowest_new .AND. MODULO(npts_out(2),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                EXIT
             ELSE
                npts_tmp(2)=npts_out(2)+1
             ENDIF
           ENDDO

           ! keep looping to find the right one
           DO
             CALL fft_radix_operations ( npts_tmp(3), npts_out(3), FFT_RADIX_NEXT )
             ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
             nlowest=npts_out(3)/2**(my_ncommensurate-my_icommensurate)
             CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
             IF (nlowest==nlowest_new .AND. MODULO(npts_out(3),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                EXIT
             ELSE
                npts_tmp(3)=npts_out(3)+1
             ENDIF
           ENDDO

         END IF ! HLFSPACE symmetry

         DO i =1,3
           pw_grid_mol% bounds ( 1, i ) = - npts_out(i) / 2
           pw_grid_mol% bounds ( 2, i ) = pw_grid_mol% bounds ( 1, i ) + npts_out(i) - 1
           side(i) = npts_out(i) * dr_g(i) 
         END DO

       END IF  ! not spherical
   
       pw_grid_mol%npts(:) = pw_grid_mol%bounds(2,:)-pw_grid_mol%bounds(1,:)+1
       pw_grid_mol%dr(:) = dr_g(:)

      ! Define the Molecular Box
       hmat(:,:) = 0.0_dp
       hmat(1,1) = side(1)
       hmat(2,2) = side(2)
       hmat(3,3) = side(3)
       CALL init_cell(cell_mol,hmat)

       npts_out( : ) = pw_grid_mol%npts( : )

       CALL create_gvectors(pw_grid_mol,cell_mol,ecut,blocking)

       pw_grid_mol%vol = ABS( cell_mol%deth )
       pw_grid_mol%dvol = pw_grid_mol%vol / REAL( pw_grid_mol%ngpts,KIND=dp )

     ! reference grid for the other igrid levels
       ref_pw_grid_mol => pw_grid_mol


     ! init pw_pools
       NULLIFY(pw_pools_mol(1)%pool)
       CALL pw_pool_create(pw_pools_mol(1)%pool,pw_grid=pw_grid_mol,error=error)

     ! init rs_pools
       drmin = MINVAL ( pw_grid_mol%dr )
!      IF ( dft_control % qs_control % rs_distributed ) THEN
!        nsmax = 2*MAX(1,CEILING((radius(1)/drmin))+1
!      ELSE
        nsmax = -1
!      ENDIF


       NULLIFY(rs_pools_mol(1)%pool)
       CALL rs_pool_create(rs_pools_mol(1)%pool,&
            el_struct=pw_grid_mol,nsmax=nsmax,error=error)
       CALL pw_grid_release(pw_grid_mol,error=error)

! And now the other grid levels if it is necessary
       IF(ngrid_level > 1) THEN

         DO igrid_level = 2,ngrid_level

           ALLOCATE(pw_grid_mol,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           CALL pw_grid_construct(pw_grid_mol)
           pw_grid_mol%para%rs_dims ( 1 ) = num_pe
           pw_grid_mol%para%rs_dims ( 2 ) = 1
           IF ( dft_control % qs_control % pw_grid_opt % spherical ) THEN
              cutilev = cutoff(igrid_level)
              pw_grid_mol%grid_span = HALFSPACE
           ELSE IF ( dft_control % qs_control % pw_grid_opt % fullspace ) THEN
              cutilev = -cutoff(igrid_level)
              pw_grid_mol%grid_span = FULLSPACE
           ELSE
              cutilev = -cutoff(igrid_level)
              pw_grid_mol%grid_span = HALFSPACE
           END IF

! ncommensurate is the number of commensurate grids 
  ! in order to have non-commensurate grids ncommensurate must be 0
  ! icommensurte  is the level number of communensurate grids
  ! this implies that the number of grid points in each direction
  ! is k*2**(ncommensurate-icommensurate) 
           my_ncommensurate = ncommensurate
           IF(my_ncommensurate .GT. 0) THEN
             my_icommensurate=igrid_level
           ELSE
             my_icommensurate=0
           ENDIF
           IF (my_icommensurate > my_ncommensurate ) THEN
              CALL stop_program ( "grid_setup", &
                            "my_icommensurate > my_ncommensurate" )
           ENDIF
           IF (my_icommensurate<=0 .AND. my_ncommensurate > 1) THEN
              CALL stop_program ( "grid_setup", &
                            " my_incommensurate<=0 .AND. my_ncommensurate > 1 " )
           ENDIF
           IF (my_ncommensurate < 0 ) THEN
              CALL stop_program ( "grid_setup", &
                            "my_ncommensurate < 0 " )
           ENDIF

         ! assign a unique tag to this grid
           grid_tag = grid_tag + 1
           pw_grid_mol % identifier = grid_tag
           pw_grid_mol % reference  = ref_pw_grid_mol%identifier
         ! this grid is treated locally
           pw_grid_mol% para % mode = PW_MODE_LOCAL
           pw_grid_mol% para % group_head = .TRUE.

           pw_grid_mol % cutoff = ABS ( cutilev )
           IF ( cutilev < 0.0_dp ) THEN
             pw_grid_mol% spherical = .FALSE.
             ecut = 1.e10_dp
           ELSE
             pw_grid_mol% spherical = .TRUE.
             ecut = cutilev
           END IF

           DO i = 1,3
             pw_grid_mol%bounds( 1, i ) = - ref_pw_grid_mol%npts(i)/ 2**igrid_level
             pw_grid_mol%bounds( 2, i ) = pw_grid_mol%bounds( 1, i ) +&
                                          ref_pw_grid_mol%npts(i)/(2**(igrid_level-1)) - 1
           END DO

           IF (.NOT. pw_grid_mol% spherical) THEN

             npts_tmp(:) = pw_grid_mol% bounds ( 2, : ) - pw_grid_mol% bounds ( 1, : ) + 1

             IF ( pw_grid_mol% grid_span == HALFSPACE .AND. symmetry ) THEN
 
               CALL fft_radix_operations ( npts_tmp(1), npts_out(1), FFT_RADIX_NEXT_ODD )
               CALL fft_radix_operations ( npts_tmp(2), npts_out(2), FFT_RADIX_NEXT_ODD )
               CALL fft_radix_operations ( npts_tmp(3), npts_out(3), FFT_RADIX_NEXT_ODD )

             ELSE

             ! keep looping to find the right one
               DO
                 CALL fft_radix_operations ( npts_tmp(1), npts_out(1), FFT_RADIX_NEXT )
               ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
                 nlowest=npts_out(1)/2**(my_ncommensurate-my_icommensurate)
                 CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
                 IF (nlowest==nlowest_new .AND. MODULO(npts_out(1),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                   EXIT
                 ELSE
                   npts_tmp(1)=npts_out(1)+1
                 ENDIF
               ENDDO
            !  keep looping to find the right one
               DO
                 CALL fft_radix_operations ( npts_tmp(2), npts_out(2), FFT_RADIX_NEXT )
               ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
                 nlowest=npts_out(2)/2**(my_ncommensurate-my_icommensurate)
                 CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
                 IF (nlowest==nlowest_new .AND. MODULO(npts_out(2),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                   EXIT
                 ELSE
                   npts_tmp(2)=npts_out(2)+1
                 ENDIF
               ENDDO
            !  keep looping to find the right one
               DO
                 CALL fft_radix_operations ( npts_tmp(3), npts_out(3), FFT_RADIX_NEXT )
               ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
                 nlowest=npts_out(3)/2**(my_ncommensurate-my_icommensurate)
                 CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
                 IF (nlowest==nlowest_new .AND. MODULO(npts_out(3),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                   EXIT
                 ELSE
                   npts_tmp(3)=npts_out(3)+1
                 ENDIF
               ENDDO
             END IF ! HLFSPACE symmetry
             DO i = 1,3
               pw_grid_mol%bounds( 1, i ) = - npts_out(i) / 2
               pw_grid_mol%bounds( 2, i ) = pw_grid_mol%bounds ( 1, i ) + npts_out(i) - 1
             END DO

           END IF  ! NOT  spherical

           pw_grid_mol%npts(:) = pw_grid_mol%bounds(2,:)-pw_grid_mol%bounds(1,:)+1
 
        ! final check if all went fine 
           IF ( ANY( MODULO(pw_grid_mol%npts,2**(my_ncommensurate-my_icommensurate)).NE.0 ) ) THEN ! nope, sorry
             CALL stop_program ( "pw_env_mol_rebuild", &
               "commensurate option failed (I) ... maybe not yet programmed for this combination of options ?" )
           END IF
           IF ( ANY(pw_grid_mol%npts * 2 ** (igrid_level-1) .NE. ref_pw_grid_mol%npts ) ) THEN
             CALL stop_program ( "pw_env_mol_rebuild", &
              "commensurate option failed (II) ... maybe not yet programmed for this combination of options ?" )
           ENDIF

         ! Check if reference grid is compatible
           IF ( pw_grid_mol%para%mode /= ref_pw_grid_mol%para%mode ) THEN
             CALL stop_program ( "pw_env_mol_rebuild", "Incompatible parallelisation scheme" )
           END IF
           IF ( pw_grid_mol%para%mode == PW_MODE_DISTRIBUTED ) THEN
             CALL stop_program ( "pw_env_mol_rebuild", "Incompatible MPI groups" )
           END IF
           IF ( pw_grid_mol%grid_span /= ref_pw_grid_mol%grid_span ) THEN
             CALL stop_program ( "pw_env_mol_rebuild", "Incompatible grid types" )
           END IF
           IF ( pw_grid_mol%spherical .NEQV. ref_pw_grid_mol%spherical ) THEN
             CALL stop_program ( "pw_env_mol_rebuild", "Incompatible cutoff schemes" )
           END IF

           CALL create_gvectors(pw_grid_mol,cell_mol,ecut,blocking,ref_grid=ref_pw_grid_mol)
  
           pw_grid_mol% vol = ABS ( cell_mol% deth )
           pw_grid_mol% dvol = pw_grid_mol% vol / REAL ( pw_grid_mol% ngpts,KIND=dp)
           pw_grid_mol% dr ( 1 ) = SQRT ( SUM ( cell_mol% hmat ( :, 1 ) ** 2 ) ) &
               / REAL ( pw_grid_mol% npts ( 1 ),KIND=dp)
           pw_grid_mol% dr ( 2 ) = SQRT ( SUM ( cell_mol% hmat ( :, 2 ) ** 2 ) ) &
               / REAL ( pw_grid_mol% npts ( 2 ),KIND=dp)
           pw_grid_mol% dr ( 3 ) = SQRT ( SUM ( cell_mol% hmat ( :, 3 ) ** 2 ) ) &
              / REAL ( pw_grid_mol% npts ( 3 ),KIND=dp)


        ! init pw_pools
           NULLIFY(pw_pools_mol(igrid_level)%pool)
           CALL pw_pool_create(pw_pools_mol(igrid_level)%pool,pw_grid=pw_grid_mol,error=error)

        ! init rs_pools
           drmin = MINVAL ( pw_grid_mol%dr )
!      IF ( dft_control % qs_control % rs_distributed ) THEN
!        nsmax = 2*MAX(1,CEILING((radius(1)/drmin))+1
!      ELSE
           nsmax = -1
!      ENDIF

           NULLIFY(rs_pools_mol(igrid_level)%pool)
           CALL rs_pool_create(rs_pools_mol(igrid_level)%pool,&
               el_struct=pw_grid_mol,nsmax=nsmax,error=error)
           CALL pw_grid_release(pw_grid_mol,error=error)

         END DO  ! igrid_level

      END IF  ! ngrid_level > 1
 
      pw_mol%pw_pools => pw_pools_mol
      pw_mol%rs_pools => rs_pools_mol

      ! init cube info
      ALLOCATE (pw_mol%cube_info(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"cube_info",0)

      cmax = 0
      
      pw_mol%cube_info(1)%max_radius =cube_info_g(1)%max_radius
      pw_mol%cube_info(1)%drmin =cube_info_g(1)%drmin
      pw_mol%cube_info(1)%max_rad_ga =cube_info_g(1)%max_rad_ga
      pw_mol%cube_info(1)%drsafe =cube_info_g(1)%drsafe
      pw_mol%cube_info(1)%dr(1:3) =cube_info_g(1)%dr(1:3)

      pw_mol%cube_info(1)%ub_cube => cube_info_g(1)%ub_cube
      pw_mol%cube_info(1)%lb_cube => cube_info_g(1)%lb_cube
      pw_mol%cube_info(1)%sphere_bounds_count => cube_info_g(1)%sphere_bounds_count
      nsmax = SIZE(cube_info_g(1)%sphere_bounds,1)
      pw_mol%cube_info(1)%sphere_bounds =>  cube_info_g(1)%sphere_bounds
      cmax=MAXVAL(cube_info_g(1)%ub_cube)

      DO igrid_level = 2,ngrid_level

        maxradius = cube_info_g(igrid_level)%max_rad_ga
        pw_mol%cube_info(igrid_level)%max_rad_ga=maxradius

        CALL init_cube_info(pw_mol%cube_info(igrid_level),&
             pw_pools_mol(igrid_level)%pool%pw_grid%dr(:),&
             maxradius,cmaxl)
        cmax = MAX(cmax,cmaxl)

      END DO


      maxlgto = l_info_g%lmax

      CALL init_l_info(pw_mol%l_info,maxlgto,cmax)

    END IF  ! failure

  END SUBROUTINE pw_env_mol_rebuild

! *****************************************************************************
!!****f* kg_gpw_pw_env_methods/update_rho_mol_blocks(
!!
!!    NAME
!!       update_rho_mol_blocks(
!!
!!    FUNCTION
!!      Preparation of the structures for the collocation of the molecular
!!      densities on the corresponding grids
!!      Multigrids can be used
!!
!!    NOTE
!!
!!    INPUTS
!!
!!    AUTHOR
!!       MI
!!
!!    MODIFICATION HISTORY
!!        Creation (MI)
!!
!!*** *************************************************************************

  SUBROUTINE update_rho_mol_blocks(kg_sub_pw_env,particle_set,cell,nspins,use_tau,error)

    TYPE(kg_sub_pw_env_type), POINTER            :: kg_sub_pw_env
    INTEGER, INTENT(IN)                          :: nspins
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                    :: particle_set
        TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                          :: use_tau
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                  :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_rho_mol_blocks',&
      routineP = module_name//':'//routineN

    INTEGER :: handle, i, iat, iatom, imol, imolecule_kind, istat, &
               nat_mol, nmol, nmolecule_kind
    INTEGER, DIMENSION(:), POINTER              :: index_atom
    LOGICAL :: failure, my_rebuild_grids, my_rebuild_ao
    REAL(dp) :: dref(3), Lo2(3), oonat_mol, r_ref(3), r0(3), r0_pbc(3), rab_pbc(3), rb(3)
    REAL(dp), DIMENSION(:,:), POINTER           :: r_mbox
    REAL(dp), DIMENSION(:,:), POINTER           :: r_tbox
    TYPE(kg_molbox_env_type), POINTER           :: molbox_env
    TYPE(cell_type) , POINTER                   :: cell_mol
    TYPE(pw_env_type) , POINTER                 :: pw_env_mol
    TYPE(pw_grid_type), POINTER                 :: pw_grid_mol
    TYPE(pw_pool_type), POINTER                 :: auxbas_pw_pool

    TYPE(qs_rho_type), POINTER                  :: rho_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER                                   :: rho_mol_blocks
    TYPE(rho_mol_blocks_type), POINTER          :: rho_block

!   ---------------------------------------------------------------------------

    failure = .FALSE.
    CALL timeset(routineN,'I',' ',handle)

    CPPrecondition(kg_sub_pw_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      nmolecule_kind = SIZE(kg_sub_pw_env%molbox_env_set,1)
      DO imolecule_kind = 1,nmolecule_kind

        NULLIFY(molbox_env, cell_mol, pw_env_mol, rho_mol_blocks)
        molbox_env => kg_sub_pw_env%molbox_env_set(imolecule_kind)
        CALL get_molbox_env(molbox_env=molbox_env, natom = nat_mol, cell_mol=cell_mol,&
                            pw_env_mol=pw_env_mol, rho_mol_blocks=rho_mol_blocks)


        CALL pw_env_get(pw_env=pw_env_mol,auxbas_pw_pool=auxbas_pw_pool,error=error)
        pw_grid_mol => auxbas_pw_pool%pw_grid

        oonat_mol = 1.0_dp/REAL(nat_mol,dp)
        nmol = SIZE(rho_mol_blocks,1)
        DO i = 1,3
!          Lo2(i) = cell_mol%hmat(i,i)/2.0_dp
          Lo2(i) = pw_grid_mol%dr(i)*FLOOR(pw_grid_mol%npts(i)/2.0_dp)

        END DO


        DO imol = 1,nmol
          NULLIFY(rho_block,rho_mol,index_atom)
          rho_block => rho_mol_blocks(imol)
          CALL get_rho_mol_block(rho_block=rho_block,rho_mol=rho_mol,index_atom=index_atom,&
                                 r_in_molbox=r_mbox,r_in_totbox=r_tbox)
          r0(1:3) = 0.0_dp
          r_ref(1:3) = particle_set( index_atom(1) )%r(1:3) 
          dref(1) = r_ref(1)-FLOOR(r_ref(1)/pw_grid_mol%dr(1))*pw_grid_mol%dr(1)
          dref(2) = r_ref(2)-FLOOR(r_ref(2)/pw_grid_mol%dr(2))*pw_grid_mol%dr(2)
          dref(3) = r_ref(3)-FLOOR(r_ref(3)/pw_grid_mol%dr(3))*pw_grid_mol%dr(3)

          DO iat = 1,nat_mol
            iatom = index_atom(iat)
            rb(1:3) = particle_set( iatom ) %r(1:3)
            rab_pbc(:) = pbc(r_ref(:),rb(:), cell)
            DO i =1,3
              r_tbox(i,iat) = r_ref(i) + rab_pbc(i)
              r0(i) = r0(i) + r_tbox(i,iat)
            END DO
          END DO
          DO i =1,3
            r0(i) = r0(i)*oonat_mol
          END DO
          CALL set_rho_mol_block(rho_block=rho_block,r0_molecule=r0)
          dref(1) = r0(1)-FLOOR(r0(1)/pw_grid_mol%dr(1))*pw_grid_mol%dr(1)
          dref(2) = r0(2)-FLOOR(r0(2)/pw_grid_mol%dr(2))*pw_grid_mol%dr(2)
          dref(3) = r0(3)-FLOOR(r0(3)/pw_grid_mol%dr(3))*pw_grid_mol%dr(3)

          DO iat = 1,nat_mol
            rb(1:3) = r_tbox(1:3,iat) 
            rab_pbc(:) = pbc(r0(:),rb(:),cell)
            DO i = 1,3
              r_mbox(i,iat) = rab_pbc(i) + Lo2(i) + dref(i) 
            END DO
          END DO

          dref(1) = r_mbox(1,1)-FLOOR(r_mbox(1,1)/pw_grid_mol%dr(1))*pw_grid_mol%dr(1)
          dref(2) = r_mbox(2,1)-FLOOR(r_mbox(2,1)/pw_grid_mol%dr(2))*pw_grid_mol%dr(2)
          dref(3) = r_mbox(3,1)-FLOOR(r_mbox(3,1)/pw_grid_mol%dr(3))*pw_grid_mol%dr(3)

          IF (.NOT.ASSOCIATED(rho_mol)) THEN
            my_rebuild_grids=.TRUE.
            my_rebuild_ao=.FALSE.
            ALLOCATE(rho_mol,stat=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            IF (.NOT.failure) THEN
               rho_mol%ref_count=1
               rho_mol%rho_g_valid=.FALSE.
               rho_mol%rho_r_valid=.FALSE.
               rho_mol%tau_g_valid=.FALSE.
               rho_mol%tau_r_valid=.FALSE.
               rho_mol%soft_valid = .FALSE.
               last_rho_id_nr=last_rho_id_nr+1
               rho_mol%id_nr=last_rho_id_nr
               rho_mol%rebuild_each=5
               NULLIFY(rho_mol%rho_r, rho_mol%rho_g, rho_mol%rho_ao, &
                       rho_mol%tot_rho_r, rho_mol%tot_rho_g, rho_mol%tau_r, rho_mol%tau_g)
               ALLOCATE(rho_mol%tot_rho_r(nspins),stat=istat)
               CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
               IF (.NOT.failure) rho_mol%tot_rho_r=0.0_dp
               CALL set_rho_mol_block(rho_block=rho_block,rho_mol=rho_mol)
               CALL qs_rho_release(rho_mol,error=error)
            END IF
          END IF  ! rho_mol
          NULLIFY(rho_mol)
          CALL get_rho_mol_block(rho_block=rho_block,rho_mol=rho_mol)

          ! rho_r
          IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho_mol%rho_r)) THEN
            IF (.NOT.failure) THEN
              IF (ASSOCIATED(rho_mol%rho_r)) THEN
                DO i=1,SIZE(rho_mol%rho_r)
                   CALL coeff_deallocate(rho_mol%rho_r(i))
                END DO
                DEALLOCATE(rho_mol%rho_r,stat=istat)
                CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
              END IF
              ALLOCATE(rho_mol%rho_r(nspins),stat=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            END IF
            IF (.NOT.failure) THEN
              DO i=1,nspins
                 CALL pw_pool_init_coeff(auxbas_pw_pool,rho_mol%rho_r(i),&
                      use_data=REALDATA3D,in_space=REALSPACE,error=error)
              END DO
            END IF
          END IF

          ! rho_g
          IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho_mol%rho_g)) THEN
            IF (.NOT.failure) THEN
              IF (ASSOCIATED(rho_mol%rho_g)) THEN
                DO i=1,SIZE(rho_mol%rho_g)
                   CALL coeff_deallocate(rho_mol%rho_g(i))
                END DO
                DEALLOCATE(rho_mol%rho_g,STAT=istat)
                CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
              END IF
              ALLOCATE(rho_mol%rho_g(nspins),STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            END IF
            IF (.NOT.failure) THEN
              DO i=1,nspins
                 CALL pw_pool_init_coeff(auxbas_pw_pool,rho_mol%rho_g(i),&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
              END DO
            END IF
          END IF

          ! allocate tau_r and tau_g for use_kinetic_energy_density
          IF (use_tau) THEN
  
            ! tau_r
            IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho_mol%tau_r)) THEN
              IF (.NOT.failure) THEN
                IF (ASSOCIATED(rho_mol%tau_r)) THEN
                  DO i=1,SIZE(rho_mol%tau_r)
                     CALL coeff_deallocate(rho_mol%tau_r(i))
                  END DO
                  DEALLOCATE(rho_mol%tau_r,stat=istat)
                  CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
                END IF
                ALLOCATE(rho_mol%tau_r(nspins),stat=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              END IF
              IF (.NOT.failure) THEN
                DO i=1,nspins
                   CALL pw_pool_init_coeff(auxbas_pw_pool,rho_mol%tau_r(i),&
                        use_data=REALDATA3D,in_space=REALSPACE,error=error)
                END DO
              END IF
            END IF

            ! tau_g
            IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho_mol%tau_g)) THEN
              IF (.NOT.failure) THEN
                IF (ASSOCIATED(rho_mol%tau_g)) THEN
                  DO i=1,SIZE(rho_mol%tau_g)
                     CALL coeff_deallocate(rho_mol%tau_g(i))
                  END DO
                  DEALLOCATE(rho_mol%tau_g,STAT=istat)
                  CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
                END IF
                ALLOCATE(rho_mol%tau_g(nspins),STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              END IF
              IF (.NOT.failure) THEN
                DO i=1,nspins
                   CALL pw_pool_init_coeff(auxbas_pw_pool,rho_mol%tau_g(i),&
                        use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
                END DO
              END IF
            END IF
          END IF  ! use_tau

        END DO  ! imol

      END DO  ! imolecule_kind
    END IF  ! failure
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE update_rho_mol_blocks

! *****************************************************************************
!!****f* kg_gpw_pw_env_methods/kg_rho_update_rho_mol
!!
!!    NAME
!!       kg_rho_update_rho_mol
!!
!!    FUNCTION
!!      Collocation of the density of each molecule on the local
!!      grid centered in the geometric center of the molecule itself
!!
!!    NOTE
!!
!!    INPUTS
!!
!!    AUTHOR
!!       MI
!!
!!    MODIFICATION HISTORY
!!        Creation (MI)
!!
!!*** *************************************************************************
  SUBROUTINE kg_rho_update_rho_mol(kg_sub_pw_env,qs_env,kg_fm_mol_set,error)

    TYPE(kg_sub_pw_env_type), POINTER            :: kg_sub_pw_env
    TYPE(qs_environment_type), POINTER           :: qs_env
    TYPE(kg_fm_mol_set_type), DIMENSION(:),&
      POINTER                                    :: kg_fm_mol_set


    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_update_rho_mol',&
      routineP = module_name//':'//routineN
  
    INTEGER                                      :: imol, imolecule_kind, ispin, &
                                                    nmol, nmolecule_kind, nspins
    INTEGER, DIMENSION(:), POINTER               :: i_atom, i_kind
    LOGICAL                                      :: failure
    REAL(dp)                                     :: dum
    REAL(dp), DIMENSION(:,:), POINTER            :: r_mbox
    TYPE(dft_control_type), POINTER              :: dft_control
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                    :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER            :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER            :: fm_mol_set
    TYPE(kg_molbox_env_type), POINTER            :: molbox_env
    TYPE(pw_env_type) , POINTER                  :: pw_env_mol
    TYPE(pw_pool_type), POINTER                  :: auxbas_pw_pool
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER                                    :: rho_mol_blocks
    TYPE(rho_mol_blocks_type), POINTER           :: rho_block
    TYPE(qs_rho_type), POINTER                   :: rho_struct
    TYPE(qs_rho_type), POINTER                   :: rho_mol


    failure=.FALSE.

    NULLIFY(dft_control, rho_struct)
    CALL get_qs_env(qs_env, dft_control=dft_control, rho=rho_struct)
    nspins=dft_control%nspins


    CPPrecondition(kg_sub_pw_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      nmolecule_kind = SIZE(kg_sub_pw_env%molbox_env_set,1)
      DO imolecule_kind = 1,nmolecule_kind

        NULLIFY(fm_mol_set,fm_mol_blocks)
        fm_mol_set => kg_fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
                               fm_mol_blocks=fm_mol_blocks)

        NULLIFY(molbox_env, pw_env_mol, rho_mol_blocks)
        molbox_env => kg_sub_pw_env%molbox_env_set(imolecule_kind)
        CALL get_molbox_env(molbox_env=molbox_env, &
                            pw_env_mol=pw_env_mol, rho_mol_blocks=rho_mol_blocks)

        nmol = SIZE(rho_mol_blocks,1)

        DO imol = 1,nmol
          NULLIFY(mol_block)
          mol_block => fm_mol_blocks(imol)
          CALL get_fm_mol_block(fm_mol_block = mol_block,&
                                index_atom = i_atom,&
                                index_kind = i_kind)

          NULLIFY(rho_block, rho_mol, r_mbox)
          rho_block => rho_mol_blocks(imol)
          CALL get_rho_mol_block(rho_block=rho_block, rho_mol=rho_mol,&
                                 r_in_molbox=r_mbox)

          DO ispin=1,nspins
            CALL calculate_rho_mol(qs_env=qs_env,&
                                   matrix_p=rho_struct%rho_ao(ispin)%matrix,&
                                   rho_r=rho_mol%rho_r(ispin),&
                                   rho_g=rho_mol%rho_g(ispin),&
                                   total_rho=rho_mol%tot_rho_r(ispin),&
                                   pw_env=pw_env_mol, &
                                   atom=i_atom, kind=i_kind, ratom=r_mbox,&
                                   error=error)
          END DO  ! ispin
          rho_mol%rho_r_valid=.TRUE.
          rho_mol%rho_g_valid=.TRUE.

          ! if needed compute also the kinetic energy density
          IF (dft_control%use_kinetic_energy_density) THEN
            DO ispin=1,nspins
              CALL calculate_rho_mol(qs_env=qs_env,&
                                     matrix_p=rho_struct%rho_ao(ispin)%matrix,&
                                     rho_r=rho_mol%tau_r(ispin),&
                                     rho_g=rho_mol%tau_g(ispin),&
                                     total_rho=dum,&  ! presumably not meaningful
                                     pw_env=pw_env_mol, &
                                     atom=i_atom, kind=i_kind, ratom=r_mbox,&
                                     compute_tau=.TRUE.,error=error)
            END DO  ! ispin
            rho_mol%tau_r_valid=.TRUE.
            rho_mol%tau_g_valid=.TRUE.
          END IF
       
        END DO  ! imol

      END DO  ! imolecule_kind

    END IF  ! failure

  END SUBROUTINE kg_rho_update_rho_mol

! *****************************************************************************
END MODULE kg_gpw_pw_env_methods
