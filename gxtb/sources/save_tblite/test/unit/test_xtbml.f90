module test_xtbml
   use mctc_env, only : wp
   use mctc_env_error, only : fatal_error
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   use mctc_io_structure, only : new
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : gemm
   use tblite_container, only : container_type
   use tblite_context_type, only : context_type
   use tblite_data_spin, only : get_spin_constant
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_param_post_processing, only : post_processing_param_list
   use tblite_param_serde, only : serde_record
   use tblite_param_xtbml_features, only : xtbml_features_record
   use tblite_post_processing_list, only : add_post_processing, post_processing_list
   use tblite_results, only : results_type
   use tblite_solvation, only : alpb_input, solvent_data, alpb_solvation, new_alpb, &
      & get_solvent_data
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_toml, only : toml_table, add_table, set_value, toml_key, get_value, &
      & toml_array, add_array
   use tblite_wavefunction , only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: acc = 0.001_wp
   real(wp), parameter :: kt = 3.166808578545117e-06_wp
   real(wp), parameter :: thr2 = 10e-4
   public :: collect_xtbml
contains

subroutine collect_xtbml(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("xtbml-mulliken-charges", test_mulliken_charges_shell_h2p),&
      new_unittest("xtbml-dipm-mol-sum-up-h2+", test_dipm_shell_h2p),&
      new_unittest("xtbml-dipm-mol-sum-up-co2", test_dipm_shell_co2),&
      new_unittest("xtbml-qp-sum-up-benzene", test_qp_shell_benz),&
      new_unittest("xtbml-qp-dipm-benzene-high-a", test_qp_shell_benz_high_a),&
      new_unittest("xtbml-energy-sum-up-gfn2", test_energy_sum_up_gfn2),&
      new_unittest("xtbml-energy-sum-up-gfn1", test_energy_sum_up_gfn1),&
      new_unittest("xtbml-rot", test_rotation_co2),&
      new_unittest("xtbml-translation", test_translation_co2),&
      new_unittest("xtbml-orbital-energy", test_orbital_energy_ref),&
      new_unittest("xtbml-param-load", test_xtbml_param_load),&
      new_unittest("xtbml-param-dump", test_xtbml_param_dump),&
      new_unittest("xtbml-param-bad-inp", test_xtbml_param_bad_inp),&
      new_unittest("xtbml-h+-orbital-energies", test_orbital_energy_hp),&
      new_unittest("xtbml-he-orbital-energies", test_orbital_energy_he),&
      new_unittest("xtbml-co2-high-spin", test_high_spin)&
      ]

end subroutine collect_xtbml

subroutine test_mulliken_charges_shell_h2p(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   type(results_type) :: res
   class(serde_record), allocatable :: tmp_record
   real(wp), allocatable :: mulliken_shell(:)
   real(wp) :: energy = 0.0_wp
   real(wp), parameter :: xyz(3, 2) = reshape((/0.0_wp,0.0_wp,0.35_wp,&
   &0.0_wp,0.0_wp,-0.35_wp/),shape=(/3,2/))
   integer, parameter :: num(2) = (/1,1/)

   call new(mol, num, xyz*aatoau, uhf=1, charge=1.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   calc%save_integrals = .true.
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("q_A", mulliken_shell)
   if (sum(mulliken_shell) - 1.0_wp > thr) then
      call test_failed(error, "Charge is not summing up to 1 electron for H2+")
      print'(3es21.14)', mulliken_shell
   end if

end subroutine test_mulliken_charges_shell_h2p


subroutine test_dipm_shell_h2p(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nsh = 2, nat=2
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), ext_partial(:)
   real(wp) :: ext_dipm_xyz(3,nat),ext_qm_xyz(6,nat)

   real(wp), parameter :: xyz(3, nat) = reshape((/0.0_wp,0.0_wp,0.35_wp,&
   &0.0_wp,0.0_wp,-0.35_wp/),shape=(/3,2/))
   integer, parameter :: num(nat) = (/1,1/)
   integer :: i


   call new(mol,num,xyz*aatoau,uhf=1,charge=1.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)
   mol_dipm = 0.0_wp

   call res%dict%get_entry("ext_dipm_A_x", tmp_array)
   ext_dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("ext_dipm_A_y", tmp_array)
   ext_dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("ext_dipm_A_z", tmp_array)
   ext_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("ext_q_A", ext_partial)

   do i= 1,2
      mol_dipm = mol_dipm+ dipm_xyz(:,i)+xyz(:,i)*partial(i)
   enddo

   if (sum(mol_dipm)  > thr) then
      call test_failed(error, "Molecular dipole moment is non zero, for dipm")
      print'(3es21.14)', mol_dipm
   end if
   mol_dipm_delta = 0.0_wp
   do i= 1,2
      mol_dipm_delta = mol_dipm_delta + ext_dipm_xyz(:,i)+xyz(:,i)*ext_partial(i)
   enddo

   if (sum(mol_dipm_delta)  > thr) then
      call test_failed(error, "Molecular dipole moment is non zero, for dipm_delta")
      print'(3es21.14)', mol_dipm_delta
   end if

   if (sum(mol_dipm_delta-mol_dipm)  > thr) then
      call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
      print'(3es21.14)', mol_dipm_delta-mol_dipm
   end if

end subroutine test_dipm_shell_h2p

subroutine test_dipm_shell_co2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=3
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy  = 0.0_wp
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), ext_partial(:)
   real(wp) :: ext_dipm_xyz(3,nat),ext_qm_xyz(6,nat)
   integer :: i,j

   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
   &0.0000,0.00000,-1.0000000,&
   &0.000000,0.000000,0.000000/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/8,8,6/)
   mol_dipm = 0.0_wp
   mol_dipm_delta = 0.0_wp

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   calc%save_integrals = .true.

   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)
   mol_dipm = 0.0_wp

   call res%dict%get_entry("ext_dipm_A_x", tmp_array)
   ext_dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("ext_dipm_A_y", tmp_array)
   ext_dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("ext_dipm_A_z", tmp_array)
   ext_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("ext_q_A", ext_partial)

   do i= 1,nat
      do j = 1,3
         mol_dipm(j) = mol_dipm(j)+ dipm_xyz(j,i)+mol%xyz(j,i)*partial(i)
      enddo
   enddo

   if (norm2(mol_dipm)  > thr) then
      call test_failed(error, "Molecular dipole moment is non zero, for dipm")
      print'(3es21.14)', mol_dipm
   end if

   do i= 1,nat
      do j = 1, 3
         mol_dipm_delta(j) = mol_dipm_delta(j) + ext_dipm_xyz(j,i)+mol%xyz(j,i)*ext_partial(i)
      enddo
   enddo

   if (norm2(mol_dipm_delta)  > thr) then
      call test_failed(error, "Molecular dipole moment is non zero, for dipm_delta")
      print'(3es21.14)', mol_dipm_delta
   end if

   if (norm2(mol_dipm_delta-mol_dipm)  > thr) then
      call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
      print'(3es21.14)', mol_dipm_delta-mol_dipm
   end if

end subroutine test_dipm_shell_co2

subroutine test_qp_shell_benz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=12
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy = 0.0_wp
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), ext_partial(:)
   real(wp) :: ext_dipm_xyz(3,nat),ext_qm_xyz(6,nat),mol_qm(6),ext_mol_qm(6)
   integer :: i,j

   real(wp), parameter :: xyz(3, nat) = reshape((/&
   &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
   &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
   &3.14559996373466,        0.73409173401801,        0.00008724271521,&
   &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
   &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
   &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
   &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
   &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
   &4.22591259698321,        0.73408789322206,        0.00025398801936,&
   &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
   &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
   &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
   &/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)

   mol_dipm = 0.0_wp
   mol_dipm_delta = 0.0_wp
   mol_qm = 0.0_wp
   ext_mol_qm = 0.0_wp
   ext_qm_xyz = 0.0_wp
   ext_dipm_xyz = 0.0_wp

   dipm_xyz = 0.0_wp
   qm_xyz = 0.0_wp


   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)


   call res%dict%get_entry("ext_dipm_A_x", tmp_array)
   ext_dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("ext_dipm_A_y", tmp_array)
   ext_dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("ext_dipm_A_z", tmp_array)
   ext_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("ext_q_A", ext_partial)

   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xx", tmp_array)
   qm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xy", tmp_array)
   qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yy", tmp_array)
   qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xz", tmp_array)
   qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yz", tmp_array)
   qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_zz", tmp_array)
   qm_xyz(6, :) = tmp_array

   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_xx", tmp_array)
   ext_qm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_xy", tmp_array)
   ext_qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_yy", tmp_array)
   ext_qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_xz", tmp_array)
   ext_qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_yz", tmp_array)
   ext_qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_zz", tmp_array)
   ext_qm_xyz(6, :) = tmp_array

   do i= 1,nat
      do j = 1,3
         mol_dipm(j) = mol_dipm(j)+dipm_xyz(j,i)+ mol%xyz(j,i)*partial(i)
      enddo
   enddo
   mol_dipm_delta = 0.0_wp
   do i= 1,nat
      do j =1,3
         mol_dipm_delta(j) = mol_dipm_delta(j) + ext_dipm_xyz(j,i)+mol%xyz(j,i)*ext_partial(i)
      enddo
   enddo

   if (norm2(mol_dipm_delta-mol_dipm)  > thr) then
      call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
      print'(3es21.14)', norm2(mol_dipm_delta-mol_dipm)
   end if


   call compute_traceless_mol_qm(mol%nat,mol%xyz,partial,dipm_xyz,qm_xyz,mol_qm)


   call compute_traceless_mol_qm(mol%nat,mol%xyz,ext_partial,ext_dipm_xyz,ext_qm_xyz,ext_mol_qm)

   if (norm2(mol_qm-ext_mol_qm)  > thr2) then
      call test_failed(error, "Molecular quadrupole moment of extended and non-extended are not equal")
      print'(3es21.14)', norm2(mol_qm-ext_mol_qm)
   end if
end subroutine test_qp_shell_benz

subroutine test_qp_shell_benz_high_a(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=12
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), ext_partial(:)
   real(wp) :: ext_dipm_xyz(3,nat),ext_qm_xyz(6,nat),mol_qm(6),ext_mol_qm(6)
   integer :: i,j

   real(wp), parameter :: xyz(3, nat) = reshape((/&
   &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
   &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
   &3.14559996373466,        0.73409173401801,        0.00008724271521,&
   &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
   &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
   &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
   &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
   &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
   &4.22591259698321,        0.73408789322206,        0.00025398801936,&
   &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
   &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
   &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
   &/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)

   mol_dipm = 0.0_wp
   mol_dipm_delta = 0.0_wp
   mol_qm = 0.0_wp
   ext_mol_qm = 0.0_wp
   ext_qm_xyz = 0.0_wp
   ext_dipm_xyz = 0.0_wp

   dipm_xyz = 0.0_wp
   qm_xyz = 0.0_wp

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1000.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)

   deallocate(tmp_array)
   call res%dict%get_entry("ext_dipm_A_x_1000.00", tmp_array)
   ext_dipm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_dipm_A_y_1000.00", tmp_array)
   ext_dipm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_dipm_A_z_1000.00", tmp_array)
   ext_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("ext_q_A_1000.00", ext_partial)

   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xx", tmp_array)
   qm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xy", tmp_array)
   qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yy", tmp_array)
   qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xz", tmp_array)
   qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yz", tmp_array)
   qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_zz", tmp_array)
   qm_xyz(6, :) = tmp_array

   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_xx_1000.00", tmp_array)
   ext_qm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_xy_1000.00", tmp_array)
   ext_qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_yy_1000.00", tmp_array)
   ext_qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_xz_1000.00", tmp_array)
   ext_qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_yz_1000.00", tmp_array)
   ext_qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("ext_qm_A_zz_1000.00", tmp_array)
   ext_qm_xyz(6, :) = tmp_array

   do i= 1,nat
      do j = 1,3
         mol_dipm(j) = mol_dipm(j)+dipm_xyz(j,i)+ mol%xyz(j,i)*partial(i)
      enddo
   enddo
   mol_dipm_delta = sum(ext_dipm_xyz,dim=2)
   if (norm2(mol_dipm_delta-mol_dipm)  > thr2) then
      call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
      print'(3es21.14)', norm2(mol_dipm_delta-mol_dipm)
   end if

   call compute_traceless_mol_qm(mol%nat,mol%xyz,partial,dipm_xyz,qm_xyz,mol_qm)

   ext_mol_qm = sum(ext_qm_xyz, dim=2)

   if (norm2(mol_qm-ext_mol_qm)  > thr2) then
      call test_failed(error, "Molecular quadrupole moment of extended and non-extended are not equal")
      print'(3es21.14)', norm2(mol_qm-ext_mol_qm)
   end if
end subroutine test_qp_shell_benz_high_a


subroutine test_rotation_co2(error)
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   !> Error handling
   integer,parameter :: nat = 3
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
   &0.0000,0.00000,-1.0000000,&
   &0.000000,0.000000,0.000000/),shape=(/3,nat/))

   integer, parameter :: num(nat) = (/8,8,6/)
   type(results_type) :: res, res_
   real(wp) :: rot_matrix(3,3),xyz_rot(3,nat),energy, xyz_trans(3,nat)
   real(wp), allocatable :: xtbml(:,:),xtbml_rot(:,:),xtbml_trans(:,:)
   integer :: i

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   xtbml_param%xtbml_orbital_energy = .true.
   xtbml_param%xtbml_geometry = .true.
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .false.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, &
      & post_process=pproc)


   rot_matrix = (reshape((/&
   &1.0_wp,0.0_wp,0.0_wp,&
   &0.0_wp,0.52532_wp,-0.85090352453_wp,&
   &0.0_wp,0.85090352453_wp,0.52532_wp/),shape=(/3,3/)))
   call gemm(rot_matrix,xyz,xyz_rot)

   call new(mol,num,xyz_rot*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   deallocate(pproc)
   allocate(pproc)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res_, post_process=pproc)

   if (.not.(compare_dict(res%dict, res_%dict, thr2))) then
      call test_failed(error, "Rotational invariance is not respected.")
   end if

end subroutine test_rotation_co2

subroutine test_translation_co2(error)
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   !> Error handling
   integer,parameter :: nat = 3
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
   &0.0000,0.00000,-1.0000000,&
   &0.000000,0.000000,0.000000/),shape=(/3,nat/))

   integer, parameter :: num(nat) = (/8,8,6/)
   type(results_type) :: res, res_
   real(wp) :: energy = 0.0_wp, xyz_trans(3,nat)
   real(wp), allocatable :: xtbml(:,:),xtbml_rot(:,:),xtbml_trans(:,:)
   integer :: i

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   xtbml_param%xtbml_geometry = .true.
   xtbml_param%xtbml_orbital_energy = .true.
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .false.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   xyz_trans = xyz(:,:)

   do i = 1,nat
      xyz_trans(1,i) = xyz_trans(1,i) +5.0_wp
   enddo

   call new(mol,num,xyz_trans*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   deallocate(pproc)
   allocate(pproc)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res_, post_process=pproc)

   if (.not.(compare_dict(res%dict, res_%dict, thr2))) then
      call test_failed(error, "Translational invariance is not respected.")
   end if

end subroutine test_translation_co2

subroutine test_orbital_energy_ref(error)
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   type(double_dictionary_type) :: dict_ref
   !> Error handling
   integer,parameter :: nat = 3
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.000000_wp, 0.000000_wp,  1.0000000_wp,&
      & 0.000000_wp, 0.000000_wp, -1.0000000_wp,&
      & 0.000000_wp, 0.000000_wp,  0.0000000_wp], shape=shape(xyz))

   integer, parameter :: num(nat) = (/8,8,6/)
   type(results_type) :: res, res_
   real(wp) :: energy = 0.0_wp
   real(wp), allocatable :: xtbml(:,:),xtbml_rot(:,:),xtbml_trans(:,:)
   integer :: i

   call dict_ref%add_entry("response_alpha", [ 0.1337835451437173_wp, 0.1337835451437157_wp, 1.6307951377273491_wp])
   call dict_ref%add_entry("gap_alpha", [ 0.0095881103175756_wp, 0.0095881103175756_wp, 0.0011085237545836_wp])
   call dict_ref%add_entry("chem_pot_alpha", [ -0.9622382191942670_wp, -0.9622382191942677_wp, -0.8758044318869692_wp])
   call dict_ref%add_entry("HOAO_alpha", [ -0.9670322743530548_wp, -0.9670322743530555_wp, -0.8763586937642610_wp])
   call dict_ref%add_entry("LUAO_alpha", [ -0.9574441640354792_wp, -0.9574441640354798_wp, -0.8752501700096774_wp])
   call dict_ref%add_entry("response_beta", [ 1.2449917888491431_wp, 1.2449917888491182_wp, 0.0059865937960335_wp])
   call dict_ref%add_entry("gap_beta", [ 0.0806953599842188_wp, 0.0806953599842184_wp, 19.4997052573882037_wp])
   call dict_ref%add_entry("chem_pot_beta", [ -15.1513328627176183_wp, -15.1513328627176254_wp, -8.6913073970860637_wp])
   call dict_ref%add_entry("HOAO_beta", [ -15.1916805427097277_wp, -15.1916805427097348_wp, -18.4411600257801638_wp])
   call dict_ref%add_entry("LUAO_beta", [ -15.1109851827255088_wp, -15.1109851827255159_wp, 1.0585452316080382_wp])
   call dict_ref%add_entry("ext_gap_alpha", [ 0.0054951699014186_wp, 0.0054951699014186_wp, 0.0067615787201293_wp])
   call dict_ref%add_entry("ext_chem_pot_alpha", [ -0.9205182203755125_wp, -0.9205182203755127_wp, -0.9334269287974830_wp])
   call dict_ref%add_entry("ext_HOAO_alpha", [ -0.9232658053262217_wp, -0.9232658053262219_wp, -0.9368077181575476_wp])
   call dict_ref%add_entry("ext_LUAO_alpha", [ -0.9177706354248032_wp, -0.9177706354248034_wp, -0.9300461394374184_wp])
   call dict_ref%add_entry("ext_gap_beta", [ 9.4593893126801003_wp, 9.4593893126801003_wp, 6.5604665198186094_wp])
   call dict_ref%add_entry("ext_chem_pot_beta", [ -12.9537156502782160_wp, -12.9537156502782196_wp, -13.9314158798429979_wp])
   call dict_ref%add_entry("ext_HOAO_beta", [ -17.6834103066182671_wp, -17.6834103066182706_wp, -17.2116491397523035_wp])
   call dict_ref%add_entry("ext_LUAO_beta", [ -8.2240209939381650_wp, -8.2240209939381685_wp, -10.6511826199336923_wp])

   call new(mol,num,xyz*aatoau,uhf=2,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_orbital_energy = .true.
   xtbml_param%xtbml_convolution = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   if (.not.(compare_dict(res%dict, dict_ref, thr2))) then
      call test_failed(error, "Comparing the orbital energy values with a reference failed")
   end if
end subroutine test_orbital_energy_ref

function compare_dict(lhs, rhs, thr_) result(equal)
   type(double_dictionary_type), intent(in) :: lhs, rhs
   real(kind=wp), intent(in) :: thr_
   logical :: equal
   integer :: i
   real(wp), allocatable :: array1(:)
   real(wp), allocatable :: array2(:, :)
   real(wp), allocatable :: array3(:, :, :)
   real(wp), allocatable :: array1_(:)
   real(wp), allocatable :: array2_(:, :)
   real(wp), allocatable :: array3_(:, :, :)
   character(len=:), allocatable :: label
   equal = .false.
   i = lhs%get_n_entries()
   if (i /= rhs%get_n_entries()) then
      write(*,*) "Not the same number of entries"
      write(*,*) "Expected: ", i, " Got: ", rhs%get_n_entries()
      return
   end if

   do i = 1, lhs%get_n_entries()

      call lhs%get_entry(i, array1)
      if (allocated(array1)) then
         call rhs%get_entry(i, array1_)
         if (allocated(array1_)) then

            if (any(abs(array1-array1_)  > thr_)) then
               call lhs%get_label(i, label)
               write(*,*) "Entry ", label, " is diverging"
               write(*,*) array1
               return
            end if
            continue
         else
            return
         end if
      end if

      call lhs%get_entry(i, array2)
      if (allocated(array2)) then
         call rhs%get_entry(i, array2_)
         if (allocated(array2_)) then
            if (any(abs(array2-array2_)  > thr_)) then
               call lhs%get_label(i, label)
               write(*,*) "Entry ", label, " is diverging"
               return
            endif
            continue
         else
            return
         end if
      end if

      call lhs%get_entry(i, array3)
      if (allocated(array3)) then
         call rhs%get_entry(i, array3_)
         if (allocated(array3_)) then
            if (any(abs(array3-array3_)  > thr_)) then
               call lhs%get_label(i, label)
               write(*,*) "Entry ", label, " is diverging"
               return
            endif
            continue
         else
            return
         end if
      end if
   end do

   equal = .true.

end function compare_dict

subroutine test_xtbml_param_load(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer :: table_xtbml
   type(toml_table) :: table_post_proc
   type(post_processing_param_list) :: param
   class(serde_record), allocatable :: record
   real(wp) :: a_array(2)
   type(toml_array), pointer :: array
   a_array =  [1.2, 1.0]

   table_post_proc= toml_table()
   call add_table(table_post_proc, "xtbml", table_xtbml)
   call set_value(table_xtbml, "geometry", .true.)
   call set_value(table_xtbml, "density", .true.)
   call set_value(table_xtbml, "orbital", .true.)
   call set_value(table_xtbml, "energy", .false.)
   call set_value(table_xtbml, "convolution", .true.)
   call add_array(table_xtbml, "a", array)
   call set_value(array, a_array)
   call set_value(table_xtbml, "tensorial-output", .true.)

   call param%load(table_post_proc, error)
   associate(record => param%list(1)%record)
      select type(record)
         type is (xtbml_features_record)
         call check(error, record%xtbml_geometry, .true.)
         if (allocated(error)) return
         call check(error, record%xtbml_density, .true.)
         if (allocated(error)) return
         call check(error, record%xtbml_orbital_energy, .true.)
         if (allocated(error)) return
         call check(error, record%xtbml_energy, .false.)
         if (allocated(error)) return
         call check(error, record%xtbml_convolution, .true.)
         if (allocated(error)) return
         call check(error, size(record%xtbml_a), 2)
         if (allocated(error)) return
         call check(error, record%xtbml_tensor, .true.)
         if (allocated(error)) return
      end select
   end associate
end subroutine test_xtbml_param_load

subroutine test_xtbml_param_bad_inp(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer :: xtbml, xtbml_
   type(toml_table) :: table_post_proc
   type(xtbml_features_record) :: param
   real(wp) :: a_array(2)
   type(toml_array), pointer :: array
   a_array =  [1.2, 1.0]

   table_post_proc = toml_table()
   call add_table(table_post_proc, "xtbml-", xtbml_)
   call param%load(table_post_proc, error)
   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad key in geometry was accepted")
      return
   endif

   call add_table(table_post_proc, "xtbml", xtbml)
   call set_value(xtbml, "geometry", 7)
   call param%load(table_post_proc, error)

   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in geometry was accepted")
      return
   endif
   call set_value(xtbml, "geometry", .true.)
   call set_value(xtbml, "density", "str")
   call param%load(table_post_proc, error)
   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in density was accepted")
      return
   endif

   call set_value(xtbml, "density", .true.)
   call set_value(xtbml, "orbital", 42)
   call param%load(table_post_proc, error)

   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in orbital was accepted")
      return
   endif

   call set_value(xtbml, "orbital", .false.)
   call set_value(xtbml, "energy", 1)
   call param%load(table_post_proc, error)

   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in energy was accepted")
      return
   endif

   call set_value(xtbml, "energy", .false.)
   call set_value(xtbml, "tensorial-output", 7)
   call param%load(table_post_proc, error)

   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in tensorial was accepted")
      return
   endif

   call set_value(xtbml, "tensorial-output", .false.)
   call set_value(xtbml, "convolution", "Str")
   call param%load(table_post_proc, error)

   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in convolution was accepted")
      return
   endif

   call set_value(xtbml, "convolution", .true.)
   call set_value(xtbml, "a", "Str")
   call param%load(table_post_proc, error)

   if (allocated(error)) then
      deallocate(error)
   else
      call fatal_error(error, "Bad input in a was accepted")
      return
   endif

   call set_value(xtbml, "a", 1.0)
   call param%load(table_post_proc, error)

end subroutine test_xtbml_param_bad_inp

subroutine test_xtbml_param_dump(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer :: table_xtbml, child
   type(toml_table) :: table_post_proc
   type(toml_table) :: new_table
   type(xtbml_features_record) :: param
   type(toml_key), allocatable :: list(:)


   table_post_proc = toml_table()
   call add_table(table_post_proc, "xtbml", table_xtbml)
   call set_value(table_xtbml, "geometry", .true.)
   call set_value(table_xtbml, "density", .true.)
   call set_value(table_xtbml, "orbital", .true.)
   call set_value(table_xtbml, "energy", .false.)
   call set_value(table_xtbml, "convolution", .true.)
   call set_value(table_xtbml, "tensorial-output", .true.)

   call param%load(table_post_proc, error)
   if (allocated(error)) return



   new_table = toml_table()
   call param%dump(new_table, error)
   call param%load(new_table, error)
   call new_table%get_keys(list)
   call check(error, size(list), 1)
   call get_value(new_table, list(1)%key, child)
   call child%get_keys(list)
   call check(error, size(list), 7)
   if (allocated(error)) return

   call check(error, param%xtbml_geometry, .true.)
   if (allocated(error)) return
   call check(error, param%xtbml_density, .true.)
   if (allocated(error)) return
   call check(error, param%xtbml_orbital_energy, .true.)
   if (allocated(error)) return
   call check(error, param%xtbml_energy, .false.)
   if (allocated(error)) return
   call check(error, param%xtbml_convolution, .true.)
   if (allocated(error)) return
   call check(error, size(param%xtbml_a), 1)
   if (allocated(error)) return
   call check(error, param%xtbml_tensor, .true.)
   if (allocated(error)) return

end subroutine test_xtbml_param_dump


subroutine test_energy_sum_up_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=12
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: wll(:, :, :)
   real(wp) :: energy = 1.0_wp , sum_energy = 0.0_wp
   integer :: i,j
   character(len=:), allocatable :: label1

   real(wp), parameter :: xyz(3, nat) = reshape((/&
   &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
   &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
   &3.14559996373466,        0.73409173401801,        0.00008724271521,&
   &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
   &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
   &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
   &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
   &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
   &4.22591259698321,        0.73408789322206,        0.00025398801936,&
   &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
   &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
   &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
   &/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)

   call new_gfn2_calculator(calc, mol, error)

   call get_spin_constants(mol, calc%bas, wll)
   allocate(calc%spin_polarization)
   call new_spin_polarization(calc%spin_polarization, mol, wll, calc%bas%nsh_id)

   block
      class(container_type), allocatable :: cont
      type(alpb_solvation), allocatable :: solv
      type(alpb_input) :: alpb_inp
      type(solvent_data), allocatable :: solvent
      solvent = get_solvent_data("water")
      alpb_inp = alpb_input(solvent%eps, solvent=solvent%solvent, &
         kernel=1, alpb=.true.)
      allocate(solv)
      call new_alpb(solv, mol, alpb_inp,"gfn2")
      call move_alloc(solv, cont)
      call calc%push_back(cont)
   end block

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("E_tot", tmp_array)
   sum_energy = sum_energy + sum(tmp_array)
   ! print'(3es21.14)', abs(sum_energy-energy)
   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN2: Energy features don't add up to total energy.")
      ! print'(3es21.14)', abs(sum_energy-energy)
   end if


   call res%dict%remove_entry("E_tot")
   call res%dict%remove_entry("w_tot")

   sum_energy = 0.0_wp
   do i = 1, res%dict%get_n_entries()
      call res%dict%get_entry(i, tmp_array)
      call res%dict%get_label(i, label1)
      sum_energy = sum_energy + sum(tmp_array)
      deallocate(tmp_array)
   end do
   ! print'(3es21.14)', abs(sum_energy-energy)
   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN2: Energy features don't add up to total energy.")
      ! print'(3es21.14)', abs(sum_energy-energy)
   end if

end subroutine test_energy_sum_up_gfn2

subroutine test_energy_sum_up_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=12
   real(wp), allocatable :: tmp_array(:)
   real(wp) :: energy = 1.0_wp, sum_energy = 0.0_wp
   integer :: i,j

   real(wp), parameter :: xyz(3, nat) = reshape((/&
   &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
   &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
   &3.14559996373466,        0.73409173401801,        0.00008724271521,&
   &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
   &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
   &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
   &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
   &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
   &4.22591259698321,        0.73408789322206,        0.00025398801936,&
   &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
   &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
   &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
   &/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn1_calculator(calc,mol,error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("E_tot", tmp_array)
   sum_energy = sum_energy + sum(tmp_array)


   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN1: Energy features don't add up to total energy.")
      print'(3es21.14)', abs(sum_energy-energy)
   end if
   call res%dict%remove_entry("E_tot")
   call res%dict%remove_entry("w_tot")

   sum_energy = 0.0_wp
   do i = 1, res%dict%get_n_entries()
      call res%dict%get_entry(i, tmp_array)
      sum_energy = sum_energy + sum(tmp_array)
      deallocate(tmp_array)
   end do


   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN1: Energy features don't add up to total energy.")
      print'(3es21.14)', abs(sum_energy-energy)
   end if

end subroutine test_energy_sum_up_gfn1

subroutine test_high_spin(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=3
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: wll(:, :, :)
   real(wp) :: energy = 1.0_wp , sum_energy = 0.0_wp
   integer :: i,j
   character(len=:), allocatable :: label1

   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
   &0.0000,0.00000,-1.0000000,&
   &0.000000,0.000000,0.000000/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/8,8,6/)

   call new(mol,num,xyz*aatoau,uhf=2,charge=0.0_wp)

   call new_gfn2_calculator(calc, mol, error)

   call get_spin_constants(mol, calc%bas, wll)
   allocate(calc%spin_polarization)
   call new_spin_polarization(calc%spin_polarization, mol, wll, calc%bas%nsh_id)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, calc%default_etemp * kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_orbital_energy = .true.
   xtbml_param%xtbml_convolution = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)
      !19 features per spin channel => 38
      !9 features for orbital energy => 18
      ! total 56
   if (res%dict%get_n_entries() /= 56) then
      call test_failed(error, "Wrong number of features")
   end if

end subroutine test_high_spin

subroutine test_orbital_energy_hp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   type(results_type) :: res
   class(serde_record), allocatable :: tmp_record
   real(wp), allocatable :: mulliken_shell(:)
   real(wp) :: energy = 0.0_wp
   real(wp), parameter :: xyz(3, 1) = reshape((/0.0_wp,0.0_wp,0.35_wp/),shape=(/3,1/))
   integer, parameter :: num(1) = (/1/)

   call new(mol, num, xyz*aatoau, uhf=0, charge=1.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   calc%save_integrals = .true.
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_orbital_energy = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("HOAO", mulliken_shell)
   if (sum(mulliken_shell) > -10.0e10_wp) then
      call test_failed(error, "HOAO is occupied")
      print'(3es21.14)', mulliken_shell
   end if

end subroutine test_orbital_energy_hp

subroutine test_orbital_energy_he(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   type(results_type) :: res
   class(serde_record), allocatable :: tmp_record
   real(wp), allocatable :: mulliken_shell(:)
   real(wp) :: energy = 0.0_wp
   real(wp), parameter :: xyz(3, 1) = reshape((/0.0_wp,0.0_wp,0.35_wp/),shape=(/3,1/))
   integer, parameter :: num(1) = (/2/)

   call new(mol, num, xyz*aatoau, uhf=0, charge=0.0_wp)
   call new_gfn1_calculator(calc, mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, calc%default_etemp * kt)
   calc%save_integrals = .true.
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_orbital_energy = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)

   call res%dict%get_entry("LUAO", mulliken_shell)
   if (sum(mulliken_shell) < 10.0e10_wp) then
      call test_failed(error, "LUAO is occupied")
      print'(3es21.14)', mulliken_shell
   end if

end subroutine test_orbital_energy_he

!this function is copied from the xtb codebase
subroutine compute_traceless_mol_qm(n,xyz,q,dipm,qp,mol_qm)
   integer :: n,i,l,k,j
   real(wp) :: xyz(3,n),dipm(3,n), qp(6,n),q(n)
   real(wp) :: mol_qm(6)
   real(wp) :: tma(6),tmb(6),tmc(6),dum
   tma = 0.0_wp
   tmb = 0.0_wp
   tmc = 0.0_wp
   do i = 1,n
      l = 0
      do j = 1,3
         do k = 1,j
            l = lin(k,j)
            tma(l) = tma(l)+xyz(j,i)*xyz(k,i)*q(i)
            tmb(l) = tmb(l)+dipm(k,i)*xyz(j,i)+dipm(j,i)*xyz(k,i)
            tmc(l) = tmc(l)+qp(l,i)
         enddo
      enddo
   enddo
   ! remove traces and multiply with 3/2 in q and dip parts
   dum = tma(1)+tma(3)+tma(6)
   dum = 0.50_wp*dum
   tma = 1.50_wp*tma
   l = 0
   do j = 1,3
      l = l+j
      tma(l) = tma(l)-dum
   enddo
   dum = tmb(1)+tmb(3)+tmb(6)
   dum = 0.50_wp*dum
   tmb = 1.50_wp*tmb
   l = 0
   do j = 1,3
      l = l+j
      tmb(l) = tmb(l)-dum
   enddo
   mol_qm = tma+tmb+tmc
end subroutine

!> this function returns the linear index of a matrix element
pure elemental integer function lin(i1,i2)
   integer,intent(in) :: i1,i2
   integer :: idum1,idum2
   idum1=max(i1,i2)
   idum2=min(i1,i2)
   lin=idum2+idum1*(idum1-1)/2
   return
end function lin

subroutine get_spin_constants(mol, bas, wll)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set data
   type(basis_type), intent(in) :: bas
   !> Spin constants
   real(wp), allocatable, intent(out) :: wll(:, :, :)

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%raw%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%raw%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants

end module test_xtbml
