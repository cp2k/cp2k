!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_build_KS_matrix [1.0] *
!!
!!   NAME
!!     qs_build_KS_matrix
!!
!!   FUNCTION
!!     routines that build the Kohn-Sham matrix (i.e calculate the coulomb
!!     and xc parts
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 moved from qs_scf (see there the history) [fawzi]
!!
!!   SOURCE
!****************************************************************************
module qs_build_KS_matrix
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type,DEBUG
  use timings, only: timeset, timestop
  use qs_environment_types, only: qs_environment_type, get_qs_env,&
       set_qs_env
  use coefficient_types,        only: coeff_allocate,&
       coeff_copy,&
       coeff_transform_space,&
       coeff_type,&
       coeff_zero,&
       coeff_sumup
  use cube_utils,               only: cube_info_type,&
       destroy_cube_info,&
       init_cube_info
  use pw_grid_types,            only: HALFSPACE,pw_grid_type
  use pw_grids,                 only: pw_find_cutoff,&
       pw_grid_construct,&
       pw_grid_setup
  use pw_types,                 only: COMPLEXDATA1D,COMPLEXDATA3D,&
       REALDATA3D,REALSPACE,&
       RECIPROCALSPACE,&
       pw_fft_wrap,pw_type
  use pws,                      only: init_pw_poisson_solver,&
       pw_poisson_solver
  use greens_fn,                only: pw_green_fn_init
  use functionals
  use gaussian_gridlevels,      only: gridlevel_info_type,&
       init_gaussian_gridlevel,& 
       destroy_gaussian_gridlevel
  use hartree,                  only: calculate_hartree
  use l_utils,                  only: destroy_l_info,&
       init_l_info,&
       l_info_type
  use qs_collocate_density,     only: calculate_rho_core,&
       calculate_rho_elec,&
       calculate_total_rho
  use qs_core_energies,         only: calculate_ecore,&
       calculate_ecore_overlap,&
       calculate_ecore_self
  use dft_types,                only: dft_control_type
  use external_potential_types, only: all_potential_type,&
       get_potential,&
       gth_potential_type
  use message_passing,          only: mp_max,mp_sum,mp_sync,&
       mp_range,mp_allgather
  use qs_collocate_density,     only: calculate_rho_core,&
       calculate_rho_elec,&
       calculate_total_rho
  use qs_core_energies,         only: calculate_ecore,&
       calculate_ecore_overlap,&
       calculate_ecore_self
  use qs_overlap,               only: write_blacs_matrix,&
       write_sparse_matrix
  use simulation_cell,          only: cell_type,&
       get_cell
  use sparse_matrix_types,      only: add_matrices,&
       allocate_matrix,&
       deallocate_matrix,&
       first_block_node,&
       get_block_node,&
       get_matrix_info,&
       next_block_node,&
       real_block_node_type,&
       real_matrix_p_type,&
       replicate_matrix_structure,&
       set_matrix,&
       symmetrise_diagonal_blocks,&
       transfer_matrix
  use atomic_kind_types,        only: atomic_kind_type,&
       get_atomic_kind,&
       get_atomic_kind_set
  use qs_integrate_potential,   only: integrate_v_rspace, &
       integrate_v_core_rspace
  use termination, only: stop_memory, stop_program
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_build_KS_matrix'

  public :: qs_ks_env_type
  public :: qs_ks_init, qs_ks_dealloc_ref, build_kohn_sham_matrix
!***
!****************************************************************************

!!****s* qs_build_KS_matrix/qs_ks_env_type [1.0] *
!!
!!   NAME
!!     qs_ks_env_type
!!
!!   FUNCTION
!!     calculation environement to calculate the ks matrix,
!!     holds all the needed vars
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  type qs_ks_env_type
     type(qs_environment_type), pointer :: qs_env
     type(pw_grid_type), dimension(:), pointer :: pw_grid
     logical :: initialized
     type(coeff_type), dimension(:), pointer       :: rho_gspace,rho_rspace,&
          rho_rspace_old,v_gspace,v_rspace
     type(coeff_type)          :: drho_rspace,v_rspace_old,vxcg_rspace,&
          work_gspace,work_rspace
     real(wp)          :: ex,ec,ecore,ecore_overlap,ecore_self,ehartree,&
          total_rho_gspace, total_rho_elec_rspace, etotal, &
          total_rho_core_rspace, total_rho_rspace
     type(gridlevel_info_type) :: gridlevel_info
     type(l_info_type)         :: l_info
     type(real_matrix_p_type)  :: v
     type(cube_info_type), dimension(:), pointer :: cube_info
     type(coeff_type)          :: rho_core
     type(coeff_type), dimension(3) :: nabla_rho_rspace
  end type qs_ks_env_type
!!***
!****************************************************************************

contains

!!****f* qs_build_KS_matrix/qs_ks_init [1.0] *
!!
!!   NAME
!!     qs_ks_init
!!
!!   FUNCTION
!!     initializes the given ks_env.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks env to be initialized
!!     - qs_env: the qs environment
!!     - global_env: the global environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine qs_ks_init(ks_env, qs_env, global_env, error)
    type(qs_ks_env_type), intent(out) :: ks_env
    type(qs_environment_type), intent(in),target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    integer, parameter :: maxgridpoints = 100
    character(len=*), parameter :: routineN='qs_ks_init',&
         routineP=moduleN//':'//routineN
    logical :: failure
    real(wp)          :: rel_cutoff,maxradius
    integer :: ngrid_level
    real(wp), dimension(:), pointer :: cutoff
    integer  :: cmax,cmaxl,igrid_level, istat,maxl
    type(dft_control_type), pointer :: dft_control
    type(real_matrix_p_type), pointer :: h
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    failure=.false.

    ks_env%qs_env => qs_env

    call init_functionals()
    call get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,atomic_kind_set=atomic_kind_set, h=h)

    rel_cutoff = dft_control%qs_control%relative_cutoff
    cutoff => dft_control%qs_control%e_cutoff
    ngrid_level = size(cutoff)

    call init_gaussian_gridlevel(ks_env%gridlevel_info,&
         ngrid_level,cutoff,rel_cutoff)

    call init_grid(ks_env,qs_env=qs_env,globenv=global_env, &
         ngrid_level=ngrid_level,error=error)

    allocate (ks_env%cube_info(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"cube_info",0)

    cmax = 0
    do igrid_level=1,ngrid_level
       maxradius = maxgridpoints*maxval(ks_env%rho_rspace(igrid_level)%pw%pw_grid%dr)
       call init_cube_info(ks_env%cube_info(igrid_level),&
            ks_env%rho_rspace(igrid_level)%pw%pw_grid%dr(:),&
            maxradius,cmaxl)
       cmax = max(cmax,cmaxl)
    end do

    call get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxl=maxl)
! maxl+1 to allow for the forces
    call init_l_info(ks_env%l_info,maxl+1,cmax)

    nullify (ks_env%v%matrix)
    call replicate_matrix_structure(h%matrix,ks_env%v%matrix,"POTENTIAL MATRIX")

    ks_env%initialized=.true.
  end subroutine qs_ks_init
!***************************************************************************

!!****f* qs_build_KS_matrix/qs_ks_dealloc_ref [1.0] *
!!
!!   NAME
!!     qs_ks_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory used by the ks_env (but not ks_env)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks_env to be released
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine qs_ks_dealloc_ref(ks_env,qs_env,global_env,error)
    type(qs_ks_env_type), intent(inout) :: ks_env
    type(qs_environment_type), intent(in),target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_ks_dealloc_ref',&
         routineP=moduleN//':'//routineN
    integer :: igrid_level,istat
    failure=.false.

    CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       nullify(ks_env%qs_env)
       call cp_unimplemented_error(routineP, "grid deallocation not done in "//&
            CPSourceFileRef,&
            error=error, error_level=cp_warning_level)

       do igrid_level=1,size(ks_env%cube_info)
          call destroy_cube_info(ks_env%cube_info(igrid_level))
       end do

       deallocate (ks_env%cube_info,STAT=istat)
       if (istat /= 0) call stop_memory(routineP,"cube_info")

       call destroy_l_info(ks_env%l_info)
       call destroy_gaussian_gridlevel(ks_env%gridlevel_info,global_env)

       ks_env%initialized=.false.
    end if
  end subroutine qs_ks_dealloc_ref
!***************************************************************************


  subroutine init_grid(ks_env,qs_env, globenv,ngrid_level, error)
    type(qs_ks_env_type), intent(inout) :: ks_env
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: globenv
    type(cp_error_type), intent(inout), optional :: error
    integer, intent(in) :: ngrid_level

!   *** Local parameters ***

    character(LEN=*), parameter :: routineN ="init_grid",&
         routineP=routineN//':'//moduleN

!   *** Local variables ***

    integer :: handle,i,igrid_level,istat,output_unit ! ,grid_info
    logical :: ionode
    type(cell_type), pointer          :: cell
    type(pw_grid_type), dimension(:), pointer :: pw_grid
    type(dft_control_type), pointer :: dft_control
    real(kind=wp), dimension(:), pointer :: cutoff

!   ---------------------------------------------------------------------------

    call timeset("init_grid","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    call get_qs_env(qs_env,cell=cell,dft_control=dft_control)
    cutoff => dft_control%qs_control%e_cutoff

!   *** Allocate the multi-grid arrays ***

    allocate (ks_env%pw_grid(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"pw_grid",0)
    allocate (ks_env%rho_rspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"rho_rspace",0)
    allocate (ks_env%rho_gspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"rho_gspace",0)
    allocate (ks_env%v_rspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"v_rspace",0)
    allocate (ks_env%v_gspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"v_gspace",0)

    allocate (ks_env%rho_rspace_old(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"rho_rspace_old",0)

!   *** Construct the various coeff_type grids ***

!!$    if (globenv%print%pw_grid_information) then
!!$       grid_info = output_unit
!!$    else
!!$       grid_info = -1
!!$    end if

    pw_grid => ks_env%pw_grid
    do igrid_level=1,ngrid_level
       call pw_grid_construct(pw_grid(igrid_level))
       pw_grid(igrid_level)%grid_span = HALFSPACE
       pw_grid(igrid_level)%para%rs_dims ( 1 ) = globenv % num_pe
       pw_grid(igrid_level)%para%rs_dims ( 2 ) = 1

       if (igrid_level == 1) then
          call pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level), &
               pe_group = globenv % group, &
               info=output_unit)
       else
          call pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
               pe_group = globenv % group, &
               info=output_unit,&
               ref_grid=pw_grid(igrid_level-1))
       end if

       call coeff_allocate(ks_env%rho_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%rho_rspace(igrid_level)%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%rho_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       ks_env%rho_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
       call coeff_allocate(ks_env%rho_rspace_old(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%rho_rspace_old%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%v_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%v_rspace(igrid_level)%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%v_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       ks_env%v_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE

    end do

!   *** Following quantities are only needed for the highest cutoff grid ***

    call pw_green_fn_init(cell%green,pw_grid(1))

    call init_pw_poisson_solver(cell%green)

    call coeff_allocate(ks_env%work_rspace,pw_grid(1),REALDATA3D)
    ks_env%work_rspace%pw%in_space = REALSPACE
    call coeff_allocate(ks_env%work_gspace,pw_grid(1),COMPLEXDATA1D)
    ks_env%work_gspace%pw%in_space = RECIPROCALSPACE
    call coeff_allocate(ks_env%v_rspace_old,pw_grid(1),REALDATA3D)
    ks_env%v_rspace_old%pw%in_space = REALSPACE
    call coeff_allocate(ks_env%rho_core,pw_grid(1),COMPLEXDATA1D)
    ks_env%rho_core%pw%in_space = RECIPROCALSPACE

    if (dft_control%gradient_functional) then
       do i=1,3
          call coeff_allocate(ks_env%nabla_rho_rspace(i),pw_grid(1),REALDATA3D)
          ks_env%nabla_rho_rspace(i)%pw%in_space = REALSPACE
       end do
       call coeff_allocate(ks_env%drho_rspace,pw_grid(1),REALDATA3D)
       ks_env%drho_rspace%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%vxcg_rspace,pw_grid(1),REALDATA3D)
       ks_env%vxcg_rspace%pw%in_space = REALSPACE
    end if

    call timestop(0.0_wp,handle)

  end subroutine init_grid

! *****************************************************************************


  subroutine build_kohn_sham_matrix(ks_env,qs_env,globenv,ks,p,p_old,rebuild,&
       nrebuild,ex,ec,ecore,ecore_overlap,&
         ecore_self,ehartree,&
         total_rho_gspace, total_rho_elec_rspace, etotal, &
         total_rho_core_rspace, total_rho_rspace,error)
    type(qs_ks_env_type), intent(inout), target :: ks_env
    type(global_environment_type), intent(IN) :: globenv
    type(qs_environment_type), intent(INOUT)  :: qs_env
    type(real_matrix_p_type), dimension(:), intent(inout) :: ks
    type(real_matrix_p_type), dimension(:), intent(in) :: p,p_old
    logical, intent(in) :: rebuild
    integer, intent(in) :: nrebuild
    real(kind=wp), intent(out), optional :: ex,ec,ecore,ecore_overlap,&
         ecore_self,ehartree,&
         total_rho_gspace, total_rho_elec_rspace, etotal, &
         total_rho_core_rspace, total_rho_rspace
    type(cp_error_type), intent(inout), optional :: error

!   *** Local variables ***

    character(len=*), parameter :: routineN='build_kohn_sham_matrix',&
         routineP=moduleN//':'//routineN
    integer :: handle,output_unit, nelectron
    logical :: ionode, failure
    type(real_matrix_p_type), pointer :: h
    type(dft_control_type), pointer :: dft_control
    failure=.false.
!   ---------------------------------------------------------------------------

    call timeset("build_kohn_sham_matrix","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    call get_qs_env(qs_env,h=h,dft_control=dft_control, n_electrons=nelectron)

    CPPrecondition(associated(h),cp_failure_level,routineP,error,failure)
    CPPrecondition(size(ks)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(size(p)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(size(p_old)>0,cp_failure_level,routineP,error,failure)

    if (.not.failure) then
       call calculate_ecore(h,p(1),ks_env%ecore,globenv)

       call calculate_rho_elec(p(1),p_old(1),rebuild,ks_env%rho_rspace,&
            ks_env%rho_rspace_old,&
            ks_env%rho_gspace,ks_env%gridlevel_info,ks_env%cube_info,&
            ks_env%l_info,&
            ks_env%total_rho_elec_rspace,qs_env,globenv)

       call coeff_zero(ks_env%v_rspace(1))
       
       ks_env%total_rho_rspace=ks_env%total_rho_elec_rspace+&
            ks_env%total_rho_core_rspace
       if (dft_control%gradient_functional) then
! we have to wait for the general xc functions since we changed
! to a complex1D data structure for the density
          call stop_program("build_kohn","feature has been removed")
!MK      CALL calculate_nabla_rho_elec(globenv)
!MK      CALL coeff_zero(vxcg_rspace)
       end if

       call calculate_xc_potential(ks_env,qs_env=qs_env,&
            globenv=globenv)

!CALL pw_fft_wrap(rho_rspace(1)%pw,rho_gspace(1)%pw,&
!                 debug=(globenv%print%level == DEBUG))

! add the core density to the electronic density
       call coeff_sumup(ks_env%rho_core,ks_env%rho_gspace(1))

       ks_env%total_rho_gspace = calculate_total_rho(ks_env%rho_gspace(1))

       if (ionode.and.globenv%print%total_densities) then
          write (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
               "Total electronic density (r-space): ",&
               ks_env%total_rho_elec_rspace,ks_env%total_rho_elec_rspace + &
               real(nelectron,wp),&
               "Total core charge density (r-space):",&
               ks_env%total_rho_core_rspace,ks_env%total_rho_core_rspace - &
               real(nelectron,wp)
          write (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total charge density (r-space):     ",ks_env%total_rho_rspace,&
               "Total charge density (g-space):     ",ks_env%total_rho_gspace
       end if

       call calculate_hartree(ks_env%rho_gspace(1),ks_env%ehartree,&
            ks_env%work_gspace)

       ks_env%etotal = ks_env%ecore_overlap + ks_env%ecore_self + &
            ks_env%ecore + ks_env%ehartree + ks_env%ex + ks_env%ec

       if (ionode.and.globenv%print%scf_energies) then
          write (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",ks_env%ecore,&
               "Hartree energy:                                ",ks_env%ehartree,&
               "Exchange-correlation energy:                   ",ks_env%ex + ks_env%ec
       end if

       call pw_fft_wrap(ks_env%work_gspace%pw,ks_env%work_rspace%pw,&
            debug=(globenv%print%level == DEBUG))

! work_rspace is the hartee potential in rspace
       call integrate_v_core_rspace(ks_env%work_rspace,ks_env%cube_info(1),&
            ks_env%l_info,qs_env,globenv)

! this following line is actually a bit ugly, since the volume element
! is only for the integration.
       ks_env%v_rspace(1)%pw%cr3d(:,:,:) = ks_env%v_rspace(1)%pw%pw_grid%dvol*&
            (ks_env%v_rspace(1)%pw%cr3d(:,:,:) +&
            ks_env%work_rspace%pw%cr3d(:,:,:))

       if (nrebuild == 1) then
          call transfer_matrix(h%matrix,ks(1)%matrix)
          call integrate_v_rspace(ks_env%v_rspace,ks_env%v_gspace,&
               ks_env%gridlevel_info,ks_env%v_rspace_old,&
               ks_env%cube_info,ks_env%l_info,rebuild,p(1),ks(1),qs_env,globenv)
       else
          if (rebuild) call set_matrix(ks_env%v%matrix,0.0_wp)
          call integrate_v_rspace(ks_env%v_rspace,ks_env%v_gspace,&
               ks_env%gridlevel_info,ks_env%v_rspace_old,&
               ks_env%cube_info,ks_env%l_info,rebuild,p(1),ks_env%v,&
               qs_env,globenv)
          call add_matrices(ks(1)%matrix,1.0_wp,h%matrix,1.0_wp,ks_env%v%matrix)
       end if

!   *** Keep a copy of the old density matrix ***

       call transfer_matrix(p(1)%matrix,p_old(1)%matrix)
    end if

    if (present(ex)) ex=ks_env%ex
    if (present(ec)) ec=ks_env%ec
    if (present(ecore)) ecore=ks_env%ecore
    if (present(ecore_overlap)) ecore_overlap=ks_env%ecore_overlap
    if (present(ecore_self)) ecore_self=ks_env%ecore_self
    if (present(ehartree)) ehartree=ks_env%ehartree
    if (present(total_rho_gspace)) total_rho_gspace=ks_env%total_rho_gspace
    if (present(total_rho_elec_rspace)) total_rho_elec_rspace=ks_env%total_rho_elec_rspace
    if (present(etotal)) etotal=ks_env%etotal
    if (present(total_rho_core_rspace)) total_rho_core_rspace=ks_env%total_rho_core_rspace
    if (present(total_rho_rspace)) total_rho_rspace=ks_env%total_rho_rspace
    
    call timestop(0.0_wp,handle)

  end subroutine build_kohn_sham_matrix

! *****************************************************************************


!!****f* cp_log_handling/calculate_xc_potential [1.0] *
!!
!!   NAME
!!     calculate_xc_potential
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 05.2002 modified to use the mp_allgather function each pe
!!       computes only part of the grid and this is broadcasted to all 
!!       instead of summed. 
!!       This scales significantly better (e.g. factor 3 on 12 cpus
!!       32 H2O) [Joost VdV]
!!     - moved to qs_build_KS_matrix [fawzi]
!!
!!*** *********************************************************************
  subroutine calculate_xc_potential(ks_env,qs_env,globenv, ec, ex)
    type(qs_ks_env_type), intent(inout),target :: ks_env
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(IN), target :: globenv
    real(kind=wp), intent(out), optional :: ec, ex

!   *** Local variables ***

    integer :: handle,kg

    integer :: local_planes,local_offset

    integer, dimension(3) :: ng

    real(wp), dimension(:,:,:), pointer :: drho,rho,vxc,vxcg
    character(len=60):: x_fun, c_fun, xc_fun
    type(dft_control_type), pointer :: dft_control


!
!

!   ---------------------------------------------------------------------------

    call timeset("calculate_xc_potential","I","",handle)

    call get_qs_env(qs_env, dft_control=dft_control)
    x_fun = dft_control%exchange_functional
    c_fun = dft_control%correlation_functional
    xc_fun = dft_control%functional

    ng(:) = ks_env%rho_rspace(1)%pw%pw_grid%npts(3)

    rho => ks_env%rho_rspace(1)%pw%cr3d(:,:,:)
    vxc => ks_env%v_rspace(1)%pw%cr3d(:,:,:)

    if (dft_control%gradient_functional) then
       drho => ks_env%drho_rspace%pw%cr3d(:,:,:)
       vxcg => ks_env%vxcg_rspace%pw%cr3d(:,:,:)
    end if

    ks_env%ex = 0.0_wp
    ks_env%ec = 0.0_wp

!D  CALL mp_range(globenv%mepos,globenv%num_pe,ng(3),local_offset,local_planes)

!D  DO kg=local_offset+1,local_offset+local_planes

!     *** Exchange functionals ***

    select case (trim(x_fun))
    case ("Becke88")
       call becke88(rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("PBE")
       call pbe_x(rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("Perdew86")
       call perdew86_x(rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("Slater")
       call slater(rho(:,:,kg),ks_env%ex,vxc(:,:,kg))
    case ("VWN")
!D      CALL vwn_x(rho(:,:,kg),ks_env%ex,vxc(:,:,kg))
       call vwn_x(rho(:,:,:),ks_env%ex,vxc(:,:,:))
    end select

!     *** Correlation functionals ***

    select case (trim(c_fun))
    case ("LYP")
       call lyp(rho(:,:,kg),drho(:,:,kg),ks_env%ec,vxc(:,:,kg),vxcg(:,:,kg))
    case ("PBE")
       call pbe_c(rho(:,:,kg),drho(:,:,kg),ks_env%ec,vxc(:,:,kg),vxcg(:,:,kg))
    case ("Perdew86")
       call perdew86_c(rho(:,:,kg),drho(:,:,kg),ks_env%ec,vxc(:,:,kg),vxcg(:,:,kg))
    case ("VWN")
!D      CALL vwn_c(rho(:,:,kg),ks_env%ec,vxc(:,:,kg))
       call vwn_c(rho(:,:,:),ks_env%ec,vxc(:,:,:))
    end select

!     *** Combined exchange-correlation functionals ***

    select case (trim(xc_fun))
    case ("HCTH/93")
       call hcth(93,rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("HCTH/120")
       call hcth(120,rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("HCTH/147")
       call hcth(147,rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("HCTH/407")
       call hcth(407,rho(:,:,kg),drho(:,:,kg),ks_env%ex,vxc(:,:,kg),vxcg(:,:,kg))
    case ("Pade")
!D      CALL pade(rho(:,:,kg),ks_env%ex,vxc(:,:,kg))
       call pade(rho(:,:,:),ks_env%ex,vxc(:,:,:))
    end select

!D  END DO

!D  CALL mp_allgather(globenv%mepos,globenv%num_pe,ng(3),vxc,globenv%group)

    if (dft_control%gradient_functional) call & 
         mp_allgather(globenv%mepos,globenv%num_pe,ng(3),vxcg,globenv%group)

    call mp_sum(ks_env%ex,globenv%group)
    call mp_sum(ks_env%ec,globenv%group)

    ks_env%ex = ks_env%rho_rspace(1)%pw%pw_grid%dvol*ks_env%ex
    ks_env%ec = ks_env%rho_rspace(1)%pw%pw_grid%dvol*ks_env%ec

!MK IF (gradient_functional) CALL calculate_nabla_vxcg(globenv)

    call timestop(0.0_wp,handle)

    if (present(ec)) ec=ks_env%ec
    if (present(ex)) ex=ks_env%ex

  end subroutine calculate_xc_potential

! *****************************************************************************


end module qs_build_KS_matrix
