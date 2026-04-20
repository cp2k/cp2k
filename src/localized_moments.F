!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Routines for the calculation of moments from Wannier functions
!> \author Lukas Schreder
! **************************************************************************************************
module localized_moments

   use atomic_kind_types
   use basis_set_types
   use cell_types
   use cp_cfm_basic_linalg
   use cp_cfm_types
   use cp_control_types
   use cp_dbcsr_cp2k_link
   use cp_dbcsr_operations
   use cp_files
   use cp_fm_basic_linalg
   use cp_fm_struct
   use cp_fm_types
   use cp_log_handling
   use cp_output_handling
   use cp_para_types
   use dbcsr_api
   use distribution_1d_types
   use input_constants
   use input_section_types
   use kinds
   use mathconstants
   use message_passing
   use molecule_kind_types
   use molecule_types
   use moments_utils
   use orbital_pointers
   use particle_types
   use physcon
   use qs_environment_types
   use qs_kind_types
   use qs_loc_types
   use qs_mo_types
   use qs_moments
   use qs_neighbor_list_types
   use qs_operators_ao
   use qs_rho_types
   use wannier_states_types

#include "./base/base_uses.f90"

   implicit none

   private

! *** globals ***

   character(len=*), parameter, private :: &
        moduleN = "localized_moments", &
        header = "# iteration, iteration_level, " // &
        "nucl. charge, el. charge, total charge, " // &
        "eldip x, eldip y, eldip z, " // &
        "elqup xx, elqup xy, elqup xz, elqup yy, elqup yz, elqup zz, " // &
        "mdip x, mdip y, mdip z, " // &
        "mom x, mom y, mom z, " // &
        "elqupvel xx, elqupvel xy, elqupvel xz, elqupvel yx, elqupvel yy, elqupvel zz, " // &
        "mdipnl x, mdipnl y, mdipnl z, momnl x, momnl y, momnl z, " // &
        "elqupvelnl xx, elqupvelnl xy, elqupvelnl xz, elqupvelnl yx, elqupvelnl yy, elqupvelnl zz, " // &
        "ref x, ref y, refz"
   integer :: iteration = 0
   integer, parameter :: nm = 9
   type(persistent_states_type), pointer :: persistent_states(:,:)
   type(cp_cfm_p_type), pointer :: ro(:,:), roder(:,:,:), po(:,:), mo(:,:)

! *** public subroutines ***

   public :: calculate_localized_moments

! *** body ***

contains

! **************************************************************************************************
!> \brief Calculates localised electric and magnetic moments from Wannier functions or Kim-Gordon
!>        subdensities
! **************************************************************************************************
   subroutine calculate_localized_moments(qs_env, mo_local)

      type(qs_environment_type), pointer     :: qs_env
      type(cp_fm_p_type), pointer, optional  :: mo_local(:)

      character(len=*), parameter  :: routineN = "calculate_localized_moments", &
                                      routineP = moduleN//":"//routineN

      character(len=default_string_length)   :: iter, name                    
      integer                                :: i, handle, iatom, iproc, jatom, imo, isubset, itimes, &
                                                natoms, &
                                                file_unit, &
                                                nao, nmo, nmolkinds, nmoments, &
                                                norders, ns, &
                                                nspins, nstates(2), nsubsets, nwcs, output_unit
      integer, allocatable                   :: first_atoms(:), first_basis(:), last_atoms(:), &
                                                last_basis(:), nbasis(:), &
                                                nelectrons(:), wcs_to_atoms( :)
      integer, pointer                       :: states(:)
      logical                                :: do_kg, do_perd, do_periodic
      real(dp)                               :: dist, r_min, rcc(3), rcc0(3), ra(3), rb(3), rpoint(3), r(3), &
           mtot, mass_low, center(3), mass, qupprefactor
      real(dp), allocatable                  :: moments(:)
      real(dp), pointer                      :: ref_point(:), rccs(:,:)
      type(atomic_kind_type), pointer        :: akind
      type(cell_type), pointer               :: cell
      type(cp_logger_type), pointer          :: logger
      type(cp_para_env_type), pointer        :: para_env
      type(dbcsr_p_type), pointer            :: densmat(:), matrix_s(:), dbop(:), pop(:)
      type(dft_control_type), pointer        :: dft_control
      type(molecule_type), pointer           :: molecules(:)
      type(neighbor_list_set_p_type), pointer :: sab_orb(:)
      type(particle_type), pointer           :: particles(:)
      !type(persistent_states_type), pointer  :: persistent_states(:, :)
      type(qs_kind_type), pointer            :: qs_kinds(:)
      type(qs_rho_type), pointer             :: rho
      type(section_vals_type), pointer       :: input
      type(wannier_centres_type), pointer    :: wcs(:)

      call timeset(routineN, handle)

      logger => cp_get_default_logger()
      iteration = iteration + 1

      ! user input
      do_periodic = .false.
      norders = 2

      ! environment infos
      call get_qs_env(qs_env, &
           cell=cell, &
           dft_control=dft_control, &
           input=input, &
           matrix_s=matrix_s, &
           molecule_set=molecules, &
           para_env=para_env, &
           particle_set=particles, &
           qs_kind_set=qs_kinds, &
           rho=rho, &
           sab_orb=sab_orb)
      ! charge(3), eldip(3), elqup(6), mdip(3), mom(3), elqupvel(6), mdipnl(3), momnl(3), elqupvelnl(6)
      nmoments = 36
      natoms = size(particles)
      nsubsets = size(molecules)
      allocate(first_atoms(nsubsets), first_basis(nsubsets), last_atoms(nsubsets), &
           last_basis(nsubsets), nbasis(nsubsets), nelectrons(nsubsets), moments(nmoments), &
           rccs(nsubsets, 3), ref_point(3))
      call get_molecule_set_info(molecules, &
           mol_to_first_atom=first_atoms, &
           mol_to_last_atom=last_atoms, &
           mol_to_nbasis=nbasis, &
           mol_to_nelectrons=nelectrons)

      ! reference points
      rcc = 0.5 * (cell%hmat(1,:) + cell%hmat(2,:) + cell%hmat(3,:))
      call get_reference_point(rcc0, qs_env=qs_env, reference=use_mom_ref_zero)
      do isubset = 1, nsubsets
         call get_reference_point(rccs(isubset, :), &
              qs_env=qs_env, &
              reference=use_mom_ref_com, &
              ifirst=first_atoms(isubset), &
              ilast=last_atoms(isubset))
      end do

      if (iteration == 1) then
         call cp_fm_get_info(mo_local(1)%matrix, ncol_global=nmo, nrow_global=nao)
         call get_qs_env(qs_env, wanniercentres=wcs)
         allocate(persistent_states(nsubsets, 1), wcs_to_atoms(nmo))
         ! states -> subset
         do imo = 1, nmo
            r_min = huge(0._dp)
            jatom = 0
            do iatom = 1, natoms
               dist = norm2(pbc(wcs(1)%centres(:, imo) - particles(iatom)%r, cell))
               if (dist < r_min) then
                  jatom = iatom
                  r_min = dist
               end if
            end do
            wcs_to_atoms(imo) = jatom
         end do
         
         ! create the subsets -> states list
         nelectrons = 0
         do isubset = 1, nsubsets
            do imo = 1, nmo
               if ((first_atoms(isubset) <= wcs_to_atoms(imo)) .and. &
                    (last_atoms(isubset) >= wcs_to_atoms(imo))) then
                  nelectrons(isubset) = nelectrons(isubset) + 1
               end if
            end do
         end do
         ! subset -> states
         do isubset = 1, nsubsets
            allocate(persistent_states(isubset, 1)%states(nelectrons(isubset)))
            nwcs = 1
            do imo = 1, nmo
               if ((first_atoms(isubset) <= wcs_to_atoms(imo)) .and. &
                    (last_atoms(isubset) >= wcs_to_atoms(imo))) then
                  persistent_states(isubset, 1)%states(nwcs) = imo
                  nwcs = nwcs + 1
               end if
            end do
            print *, "dbg| states(isubset)", isubset, size(persistent_states(isubset, 1)%states), &
                 persistent_states(isubset, 1)%states
         end do
         deallocate(wcs_to_atoms)

         ! create operators
         call create_operators(qs_env, rcc, rcc0, rccs, mo_local)
      end if

      call calculate_loc_moments(qs_env, rcc, mo_local, opcode=-1, tag="ref")
      call calculate_loc_moments(qs_env, rcc0, mo_local, opcode=0, tag="ref0")

      ! subsystems
      do isubset = 1, nsubsets
         states => persistent_states(isubset, 1)%states
         name = molecules(isubset)%molecule_kind%name
         call calculate_loc_moments(qs_env, rccs(isubset, :), mo_local, isubset, states, &
              trim(name))
      end do
      
      deallocate(first_atoms, first_basis, last_atoms, last_basis, moments, nbasis, nelectrons, &
           rccs, ref_point)

      call timestop(handle)
     
   end subroutine calculate_localized_moments
! **************************************************************************************************

! **************************************************************************************************
   subroutine create_operators(qs_env, rcc, rcc0, rccs, mos)

     type(qs_environment_type), pointer :: qs_env
     real(dp)                           :: rcc(3), rcc0(3)
     real(dp), pointer                  :: rccs(:,:)
     type(cp_fm_p_type), pointer        :: mos(:)

     character(*), parameter :: routinen = "create_operators"

     integer :: file_unit, handle, i, isubset, j, nao, nmo, nmm, nsubsets
     real(dp), pointer :: myrcc(:)
     type(cp_fm_struct_type), pointer :: o_struct
     type(cp_fm_type), pointer :: or
     type(dbcsr_p_type), pointer :: matrix_s(:), rop(:), ropder(:,:), mop(:), pop(:)
     type(neighbor_list_set_p_type), pointer :: sab_all(:), sab_orb(:)

     call timeset(routinen, handle)

     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
     call cp_fm_get_info(mos(1)%matrix, nrow_global=nao, ncol_global=nmo)

     nsubsets = size(rccs, 1)

     allocate( &
          ro(-1:size(rccs, 1), 9), &
          roder(-1:size(rccs, 1), 3, 3), &
          po(-1:size(rccs, 1), 3), &
          mo(-1:size(rccs, 1), 3))

     call cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
          template_fmstruct=mos(1)%matrix%matrix_struct)
     call cp_fm_create(or, o_struct)  ! tmp

     do isubset = -1, nsubsets

        print *, "create operator", isubset

        if (isubset == -1) then
           allocate(myrcc(3))
           myrcc(:) = rcc(:)
        else if (isubset == 0) then
           allocate(myrcc(3))
           myrcc(:) = rcc0(:)
        else
           myrcc => rccs(isubset, :)
        end if

        ! r and derivatives
        nullify(rop)
        call dbcsr_allocate_matrix_set(rop, nm)
        do i = 1, nm
           allocate(rop(i)%matrix)
           call dbcsr_create(rop(i)%matrix, template=matrix_s(1)%matrix, &
                matrix_type=dbcsr_type_no_symmetry)
           call cp_dbcsr_alloc_block_from_nbl(rop(i)%matrix, sab_all)
           call dbcsr_set(rop(i)%matrix, 0._dp)
        end do
        nullify(ropder)
        call dbcsr_allocate_matrix_set(ropder, 3, 3)
        do i = 1, 3
           do j = 1, 3
              call dbcsr_init_p(ropder(i,j)%matrix)
              call dbcsr_create(ropder(i,j)%matrix, template=matrix_s(1)%matrix, &
                   matrix_type=dbcsr_type_no_symmetry)
              call cp_dbcsr_alloc_block_from_nbl(ropder(i,j)%matrix, sab_all)
              call dbcsr_set(ropder(i,j)%matrix, 0._dp)
           end do
        end do
        call build_local_moments_der_matrix(qs_env, rop, ropder, 2, 1, myrcc)
        do i = 1, nm
           call cp_cfm_create(ro(isubset, i)%matrix, o_struct)
           call copy_dbcsr_to_fm(rop(i)%matrix, or)
           call dbcsr_release(rop(i)%matrix)
           deallocate(rop(i)%matrix)
           call cp_fm_to_cfm(or, mtarget=ro(isubset, i)%matrix)
        end do
        deallocate(rop)
        do i = 1, 3
           do j = 1, 3
              call cp_cfm_create(roder(isubset,i,j)%matrix, o_struct)
              call copy_dbcsr_to_fm(ropder(i,j)%matrix, or)
              call dbcsr_release(ropder(i,j)%matrix)
              deallocate(ropder(i,j)%matrix)
              call cp_fm_to_cfm(or, mtarget=roder(isubset,i,j)%matrix)
           end do
        end do
        deallocate(ropder)

        ! m
        nullify(mop)
        call build_local_magmom_matrix(qs_env, mop, 1, myrcc)
        nmm = size(mop)
        do i = 1, nmm
           call cp_cfm_create(mo(isubset,i)%matrix, o_struct)
           call copy_dbcsr_to_fm(mop(i)%matrix, or)
           call dbcsr_release(mop(i)%matrix)
           deallocate(mop(i)%matrix)
           call cp_fm_to_cfm(or, mtarget=mo(isubset,i)%matrix)
        end do
        deallocate(mop)

        ! p
        nullify(pop)
        call dbcsr_allocate_matrix_set(pop, 3)
        do i = 1, 3
           call dbcsr_init_p(pop(i)%matrix)
           call dbcsr_create(pop(i)%matrix, template=matrix_s(1)%matrix, &
                matrix_type=dbcsr_type_antisymmetric)
           call cp_dbcsr_alloc_block_from_nbl(pop(i)%matrix, sab_orb)
           call dbcsr_set(pop(i)%matrix, 0._dp)
        end do
        call build_lin_mom_matrix(qs_env, pop)
        do i = 1, 3
           call cp_cfm_create(po(isubset,i)%matrix, o_struct)
           call copy_dbcsr_to_fm(pop(i)%matrix, or)
           call dbcsr_release(pop(i)%matrix)
           deallocate(pop(i)%matrix)
           call cp_fm_to_cfm(or, mtarget=po(isubset,i)%matrix)
        end do
        deallocate(pop)

        if (isubset == -1) then
           deallocate(myrcc)
        else if (isubset == 0) then
           deallocate(myrcc)
        end if
     end do
     
     call cp_fm_release(or)  ! tmp
     call cp_fm_struct_release(o_struct)
     
     call timestop(handle)

   end subroutine create_operators
! **************************************************************************************************


! **************************************************************************************************
   subroutine calculate_loc_moments(qs_env, rcc, mos, opcode, statelist, tag)

     type(qs_environment_type), pointer  :: qs_env
     real(dp)                            :: rcc(3)
     type(cp_fm_p_type), pointer         :: mos(:)
     integer                             :: opcode
     integer, pointer, optional          :: statelist(:)
     character(*), optional              :: tag

     character(*), parameter :: routinen = "calculate_loc_moments"

     complex(dp) :: di
     integer :: file_unit, handle, i, j, k, nao, nmm, nmo, ns
     integer, pointer :: states(:)
     logical :: do_rtp
     real(dp) :: d, mu(3), muv(3), m(3), th(6), thder(9), thv(6), q, zwfc
     real(dp), allocatable :: rv(:), rrv(:), rxrv(:)
     !type(cp_cfm_p_type), pointer :: mo(:), po(:), ro(:), roder(:,:)
     type(cp_cfm_type), pointer :: o, ov, v, vov
     type(cp_fm_struct_type), pointer :: v_struct, o_struct, vov_struct
     type(cp_fm_type), pointer :: or, vr, vi
     type(dbcsr_p_type), pointer :: matrix_s(:)!, rop(:), ropder(:, :), mop(:), pop(:)
     type(neighbor_list_set_p_type), pointer :: sab_all(:), sab_orb(:)

     call timeset(routinen, handle)

     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
     call cp_fm_get_info(mos(1)%matrix, nrow_global=nao, ncol_global=nmo)
     zwfc = 2._dp
     do_rtp = (size(mos) == 2)

     if (present(statelist)) then
        ns = size(statelist)
        states => statelist
     else
        ns = nmo
        allocate(states(ns))
        do i = 1, ns
           states(i) = i
        end do
     end if

     ! matrices
     call cp_fm_struct_create(v_struct, nrow_global=nao, ncol_global=ns, &
          template_fmstruct=mos(1)%matrix%matrix_struct)
     call cp_fm_struct_create(vov_struct, nrow_global=ns, ncol_global=ns, &
          template_fmstruct=mos(1)%matrix%matrix_struct)

     call cp_cfm_create(ov, v_struct, name="ov")
     call cp_cfm_create(vov, vov_struct, name="vov")
     call cp_fm_struct_release(vov_struct)

     ! MOs
     call cp_fm_create(vr, v_struct, name="vr")
     call cp_fm_create(vi, v_struct, name="vi")
     call cp_cfm_create(v, v_struct, name="v")
     call cp_fm_struct_release(v_struct)
     do i = 1, ns
        call cp_fm_to_fm(mos(1)%matrix, vr, 1, states(i), i)
        if (do_rtp) then
           call cp_fm_to_fm(mos(2)%matrix, vi, 1, states(i), i)
           call cp_fm_to_cfm(vr, vi, v)
        else
           call cp_fm_to_cfm(vr, mtarget=v)
        end if
     end do
     call cp_fm_release(vr)
     call cp_fm_release(vi)
     
     ! charge
     q = -ns * zwfc
     call open_file("q_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(E14.8)") q
     end if
     call close_file(file_unit)

     ! mu
     mu = 0._dp
     do i = 1, 3
        o => ro(opcode, i)%matrix
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
        d = 0._dp
        do j = 1, ns
           call cp_cfm_get_element(vov, j, j, di)
           d = d + real(di)
        end do
        mu(i) = -zwfc * d
     end do
     call open_file("mu_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(2(E14.8,','),E14.8)") mu
     end if
     call close_file(file_unit)

     ! muv
     muv = 0._dp
     do i = 1, 3
        o => po(opcode, i)%matrix
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
        d = 0._dp
        do j = 1, ns
           call cp_cfm_get_element(vov, j, j, di)
           d = d + aimag(di)
        end do
        muv(i) = -zwfc * d
     end do
     call open_file("muv_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(2(E14.8,','),E14.8)") muv
     end if
     call close_file(file_unit)

     ! m
     m = 0._dp
     do i = 1, 3
        o => mo(opcode, i)%matrix
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
        d = 0._dp
        do j = 1, ns
           call cp_cfm_get_element(vov, j, j, di)
           d = d + aimag(di)
        end do
        m(i) = - zwfc * d
     end do
     call open_file("m_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(2(E14.8,','),E14.8)") m
     end if
     call close_file(file_unit)

     ! th
     th = 0._dp
     do i = 1, 6
        o => ro(opcode, i+3)%matrix
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
        d = 0._dp
        do j = 1, ns
           call cp_cfm_get_element(vov, j, j, di)
           d = d + zwfc * real(di)
        end do
        th(i) = -d
     end do
     call open_file("th_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(5(E14.8,','),E14.8)") th
     end if
     call close_file(file_unit)

     ! thv
     thv = 0._dp
     do i = 1, 3
        do j = 1, 3
           o => roder(opcode, i, j)%matrix
           call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
           call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
           d = 0._dp
           do k = 1, ns
              call cp_cfm_get_element(vov, k, k, di)
              d = d - aimag(di)
              if (i == j) then
                 d = d + real(di)
              end if
           end do
           thder((i-1)*3+j) = zwfc * d
        end do
     end do
     thv(1) = -2 * thder(1) - q
     thv(2) = -thder(2) - thder(4)
     thv(3) = -thder(3) - thder(7)
     thv(4) = -2 * thder(5) - q
     thv(5) = -thder(6) - thder(8)
     thv(6) = -2 * thder(9) - q
     call open_file("thv_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(5(E14.8,','),E14.8)") thv
     end if
     call close_file(file_unit)

     if (.not. present(statelist)) then
        deallocate(states)
     end if
     call cp_cfm_release(v)
     call cp_cfm_release(ov)
     call cp_cfm_release(vov)

     ! nl terms
     allocate(rv(3), rrv(6), rxrv(3))
     rv = 0._dp
     rrv = 0._dp
     rxrv = 0._dp
     call calculate_commutator_nl_terms(qs_env, rv, rxrv, rrv, rcc)

     call open_file("rv_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(2(E14.8,','),E14.8)") rv
     end if
     call close_file(file_unit)

     call open_file("rrv_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(5(E14.8,','),E14.8)") rrv
     end if
     call close_file(file_unit)

     call open_file("rxrv_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(2(E14.8,','),E14.8)") rxrv
     end if
     call close_file(file_unit)

     call open_file("rcc_"//tag//".csv", unit_number=file_unit, file_action="write", &
          file_position="append", file_status="unknown")
     if (file_unit > 0) then
        write (file_unit, "(2(E14.8,','),E14.8)") rcc
     end if
     call close_file(file_unit)

     call timestop(handle)

   end subroutine calculate_loc_moments


! **************************************************************************************************
! **************************************************************************************************
   subroutine build_wheeler_moment_matrix(qs_env, o, rcc, m)

      type(qs_environment_type), pointer  :: qs_env
      type(cp_cfm_p_type), pointer        :: o(:)
      real(dp)  :: rcc(3)
      type(cp_cfm_type), pointer  :: m

      integer  :: i, j, k, nao, nmo
      real(dp)  :: kvec(3)
      type(cp_fm_struct_type), pointer  :: o_struct
      type(cp_fm_type), pointer  :: cmat, smat
      type(dbcsr_p_type), pointer  :: matrix_s(:), qo(:)
      type(dbcsr_type), pointer  :: cosmat, sinmat
      type(cell_type), pointer  :: cell

      call get_qs_env(qs_env, cell=cell, matrix_s=matrix_s)

      call cp_cfm_get_info(m, nrow_global=nao)
      call cp_fm_struct_create(o_struct, ncol_global=nao, template_fmstruct=m%matrix_struct)
      call cp_fm_create(cmat, o_struct)
      call cp_fm_create(smat, o_struct)
      nullify(cosmat, sinmat)
      allocate(cosmat, sinmat)
      call dbcsr_copy(cosmat, matrix_s(1)%matrix)
      call dbcsr_copy(sinmat, matrix_s(1)%matrix)
      allocate(o(9))
      do i = 1, 9
         call cp_cfm_create(o(i)%matrix, o_struct)
      end do
      
      !do i = 1, 3
      !   kvec = twopi * cell%h_inv(i, :)
      !   call build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
      !   call copy_dbcsr_to_fm(cosmat, cmat)
      !   call copy_dbcsr_to_fm(sinmat, smat)
      !   smat%local_data = -smat%local_data
      !   call cp_fm_to_cfm(cmat, smat, o(i)%matrix)
      !end do
      !k = 4
      !do i = 1, 3
      !   do j = i, 3
      !      kvec = 0.0_dp
      !      kvec(i) = twopi * cell%h_inv(i, i) * cell%h_inv(j, j)
      !      !print *, "dbg| kvec", kvec
      !      call build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec, q=.true.)
      !      call copy_dbcsr_to_fm(cosmat, cmat)
      !      call copy_dbcsr_to_fm(sinmat, smat)
      !      smat%local_data = -smat%local_data
      !      call cp_fm_to_cfm(cmat, smat, o(k)%matrix)
      !      k = k + 1
      !   end do
      !end do
      nullify(qo)
      call dbcsr_allocate_matrix_set(qo, 9)
      do i = 1, 9
         allocate(qo(i)%matrix)
         call dbcsr_copy(qo(i)%matrix, matrix_s(1)%matrix, "qo")
      end do
      call build_local_moment_matrix(qs_env, qo, 2, ref_point=rcc)
      do i = 4, 9
         call copy_dbcsr_to_fm(qo(i)%matrix, cmat)
         call cp_fm_to_cfm(cmat, mtarget=o(i)%matrix)
      end do
      !print *, "dbg| o(4)", o(4)%matrix%local_data
      !qo => o(4:9)
      !call build_wheeler_quadrupole_matrix(qs_env, qo)
      call cp_fm_release(cmat)
      call cp_fm_release(smat)
      call dbcsr_release(cosmat)
      call dbcsr_release(sinmat)
      call cp_fm_struct_release(o_struct)
      deallocate(cosmat, sinmat)
      
   end subroutine build_wheeler_moment_matrix


! **************************************************************************************************
! **************************************************************************************************
   subroutine build_modernmag_matrix(qs_env, o, rcc, m)

      type(qs_environment_type), pointer  :: qs_env
      type(cp_cfm_p_type), pointer        :: o(:)
      real(dp)                            :: rcc(3)
      type(cp_cfm_type), pointer  :: m

      integer  :: i, nao
      type(cp_fm_struct_type), pointer  :: o_struct
      type(cp_fm_type), pointer  :: or
      type(dbcsr_p_type), pointer  :: op(:), matrix_s(:)
      type(neighbor_list_set_p_type), pointer  :: sab_orb(:)
      
      call get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
      call cp_cfm_get_info(m, nrow_global=nao)
      call cp_fm_struct_create(o_struct, ncol_global=nao, template_fmstruct=m%matrix_struct)
      call cp_fm_create(or, o_struct)
      
      nullify(op)
      call dbcsr_allocate_matrix_set(op, 3)
      do i = 1, 3
         call dbcsr_init_p(op(i)%matrix)
         call dbcsr_create(op(i)%matrix, template=matrix_s(1)%matrix, &
              matrix_type=dbcsr_type_antisymmetric)
         call cp_dbcsr_alloc_block_from_nbl(op(i)%matrix, sab_orb)
         call dbcsr_set(op(i)%matrix, 0._dp)
      end do
      call build_local_magmom_matrix(qs_env, op, 2, ref_point=rcc)
      !print *, "dbg| rcc build_modernmag_matrix", rcc
      !call dbcsr_print(op(1)%matrix, unit_nr=7)
      allocate(o(3))
      do i = 1, 3
        call cp_cfm_create(o(i)%matrix, o_struct, name="mago")
        call copy_dbcsr_to_fm(op(i)%matrix, or)
        call dbcsr_release(op(i)%matrix)
        deallocate(op(i)%matrix)
        call cp_fm_to_cfm(or, mtarget=o(i)%matrix)
      end do
      deallocate(op)
      call cp_fm_release(or)
      call cp_fm_struct_release(o_struct)

   end subroutine build_modernmag_matrix


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_reference_nuc_moments(moments, qqi, qqij, particles, qs_kinds, cell, rcc)
!
!      real(dp)                      :: moments(:), qqi(3), qqij(3, 3)
!      type(particle_type), pointer  :: particles(:)
!      type(qs_kind_type), pointer   :: qs_kinds(:)
!      type(cell_type), pointer      :: cell
!      real(dp)                      :: rcc(3)
!
!      complex(dp)                      :: d, di(3), dij(3, 3), iphase(3), ijphase(3, 3)
!      integer                          :: akind, i, ix, iy, iz, j, k, l
!      real(dp)                         :: dir(3), kvec(3), r(3), z, dd
!      type(atomic_kind_type), pointer  :: atomic_kind
!
!      do k = 0, 2
!         select case (k)
!         case (0)
!            ! charge
!            do i = 1, size(particles)
!               atomic_kind => particles(i)%atomic_kind
!               call get_atomic_kind(atomic_kind, kind_number=akind)
!               call get_qs_kind(qs_kinds(akind), core_charge=z)
!               moments(1) = moments(1) + z
!            end do
!            r = twopi * matmul(cell%h_inv, rcc)
!            iphase = exp(-gaussi * r)**moments(1)
!            do i = 1, 3
!               do j = 1, 3
!                  ijphase(i, j) = exp(-gaussi * r(i) * r(j))**moments(1)
!               end do
!            end do
!         case (1)
!            ! dipole
!            di = (1._dp, 0._dp)
!            do l = 1, size(particles)
!               r = pbc(particles(l)%r, cell)
!               atomic_kind => particles(l)%atomic_kind
!               call get_atomic_kind(atomic_kind, kind_number=akind)
!               call get_qs_kind(qs_kinds(akind), core_charge=z)
!               di = di * (exp(gaussi * twopi * matmul(cell%h_inv, r))**z)**z
!            end do
!            di = di * iphase
!            qqi = aimag(log(di))  ! for phase periodicity
!            moments(4:6) = matmul(cell%hmat / twopi, aimag(log(di)))
!         case (2)
!            ! quadrupole
!            dij = (1._dp, 0._dp)
!            do l = 1, size(particles)
!               r = pbc(particles(l)%r, cell)
!               atomic_kind => particles(l)%atomic_kind
!               call get_atomic_kind(atomic_kind, kind_number=akind)
!               call get_qs_kind(qs_kinds(akind), core_charge=z)
!               do i = 1, 3
!                  do j = 1, 3
!                     kvec(1) = twopi * cell%h_inv(i, i) * cell%h_inv(j, j)
!                     dd = kvec(1) * r(i) * r(j)
!                     d = exp(gaussi * dd)**z
!                     dij(i, j) = dij(i, j) * d
!                  end do
!               end do
!            end do
!            dij = dij * ijphase
!            qqij = aimag(log(dij))
!            l = 1
!            do i = 1, 3
!               do j = i, 3
!                  moments(l + 6) = cell%hmat(i, i) * cell%hmat(j, j) / twopi * aimag(log(dij(i, j)))
!                  l = l + 1
!               end do
!            end do
!         end select
!      end do
!
!   end subroutine calculate_reference_nuc_moments


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_reference_moments(moments, qs_env, rcc, particles, qs_kinds, cell, rho, &
!        mo_local)
!
!     real(dp)                            :: moments(:)
!     type(qs_environment_type), pointer  :: qs_env
!     real(dp)                            :: rcc(3)
!     type(particle_type), pointer        :: particles(:)
!     type(qs_kind_type), pointer         :: qs_kinds(:)
!     type(cell_type), pointer            :: cell
!     type(qs_rho_type), pointer          :: rho
!     type(cp_fm_p_type), pointer         :: mo_local(:)
!
!     complex(dp)  :: cd, cda, iphase(3), ijphase(3, 3), di(3), dij(3, 3)
!     integer                                  :: i, j, k, l, nao, nmo, tmp_dim
!     logical                                  :: do_rtp, found, tr
!     real(dp)                                 :: focc, kvec(3), qupoleder(9), strace, trace, &
!          r(3), zwfc, qqi(3), ci(3), qqij(3, 3), cij(3, 3), d
!     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
!     real(dp), pointer                        :: blk_p(:)
!     type(cp_cfm_type), pointer               :: m, om, mom
!     type(cp_cfm_p_type), pointer             :: o(:), mago(:), eigrmat(:)
!     type(dbcsr_p_type), pointer              :: cosmat(:), sinmat(:)
!     type(cp_fm_struct_type), pointer         :: m_struct, o_struct, mom_struct
!     type(cp_fm_p_type), pointer              :: op_fm_set(:, :), opvec(:)
!     type(cp_fm_type), pointer                :: mr, mi, or
!     type(dbcsr_iterator_type), pointer       :: iterator
!     type(dbcsr_p_type), pointer              :: oder(:, :), matrix_s(:), &
!          rho_ao(:), rho_ao_im(:)
!     type(dbcsr_type), pointer                :: densmat, densmat_im, ocos, osin, tmp
!     type(mo_set_p_type), pointer             :: mos(:)
!     type(neighbor_list_set_p_type), pointer  :: sab_all(:), sab_orb(:)
!     type(dbcsr_p_type), pointer  :: magop(:)
!
!     print *, "dbg| ref"
!     print *, "dbg| rcc calc_ref_mom", rcc
!     
!     ! info
!     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb, mos=mos)
!     do_rtp = (size(mo_local) == 2)
!     call cp_fm_get_info(mo_local(1)%matrix, nrow_global=nao, ncol_global=nmo)
!     focc = 1._dp
!     zwfc = 2._dp
!     trace = 0._dp
!
!     ! structs
!     m_struct => mo_local(1)%matrix%matrix_struct
!     call cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
!          template_fmstruct=mo_local(1)%matrix%matrix_struct)
!     call cp_fm_struct_create(mom_struct, nrow_global=nmo, ncol_global=nmo, &
!          template_fmstruct=mo_local(1)%matrix%matrix_struct)
!
!     ! get the mo coeffs
!     call cp_fm_create(mr, m_struct, name="mr")
!     if (do_rtp) then
!        call cp_fm_create(mi, m_struct, name="mi")
!     end if
!     call cp_cfm_create(m, m_struct, name="m")
!     call cp_fm_to_fm(mo_local(1)%matrix, mr)
!     if (do_rtp) then
!        call cp_fm_to_fm(mo_local(2)%matrix, mi)
!        call cp_fm_to_cfm(mr, mi, m)
!     else
!        call cp_fm_to_cfm(mr, mtarget=m)
!     end if
!     call cp_fm_release(mr)
!     if (do_rtp) then
!        call cp_fm_release(mi)
!     end if
!     call cp_fm_create(or, o_struct, name="mi")
!     call cp_fm_struct_release(o_struct)
!
!     ! operators
!     call build_wheeler_moment_matrix(qs_env, o, rcc, m)
!     call build_modernmag_matrix(qs_env, mago, rcc, m)
!     
!     ! get the nonlocal terms
!     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
!     nlrv = 0._dp
!     nlrrv = 0._dp
!     nlrxrv = 0._dp
!     call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, ref_point=rcc)
!     
!     ! nuclear contribution
!     call calculate_reference_nuc_moments(moments, qqi, qqij, particles, qs_kinds, cell, rcc)
!     
!     ! work matrices
!     call cp_cfm_create(om, m_struct)
!     call cp_cfm_create(mom, mom_struct)
!     call cp_fm_struct_release(mom_struct)
!
!     ! charge
!     moments(2) = -zwfc * nmo
!     moments(3) = moments(1) + moments(2)
!     ! initial phase from the reference point
!     r = twopi * (-moments(2)) * matmul(cell%h_inv, rcc)
!     iphase = exp(gaussi * r)
!     do i = 1, 3
!        do j = 1, 3
!           !ijphase(i, j) = exp(twopi * (-moments(2)) * cell%h_inv(i, i) * cell%h_inv(j, j) * rcc(i) * rcc(j))
!        end do
!     end do
!
!     ! edip
!     di = (1._dp, 0._dp)
!     do i = 1, 3
!        call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
!        call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!        d = 0._dp
!        do j = 1, nmo
!           call cp_cfm_get_element(mom, j, j, cd)
!           d = d + real(cd)
!        end do
!        moments(i+3) = moments(i+3) + d
!     end do
!     
!     ! mdip
!     do i = 1, 3
!        call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), mago(i)%matrix, m, (0._dp, 0._dp), om)
!        call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!        d = 0._dp
!        do j = 1, nmo
!           call cp_cfm_get_element(mom, j, j, cd)
!           d = d + aimag(cd)
!        end do
!        moments(i+12) = moments(i+12) - zwfc * d
!     end do
!     ! nl
!     moments(25:27) = moments(13:15) + nlrxrv
!
!     ! equp
!     dij = (1._dp, 0._dp)
!     k = 4
!     do i = 1, 3
!        do j = i, 3
!           call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), o(k)%matrix, m, (0._dp, 0._dp), om)
!           call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!           d = 0._dp
!           do k = 1, nmo
!              call cp_cfm_get_element(mom, k, k, cd)
!              d = d + real(cd)
!           end do
!        end do
!     end do
!     cij = aimag(log(dij))
!     l = 1
!     do i = 1, 3
!        do j = i, 3
!           if (qqij(i, j) + cij(i, j) > pi) cij(i, j) = cij(i, j) - twopi
!           if (qqij(i, j) + cij(i, j) < -pi) cij(i, j) = cij(i, j) + twopi
!           cij(i, j) = cell%hmat(i, i) * cell%hmat(j, j) / twopi
!           moments(l + 6) = moments(l + 6) + cij(i, j)
!           l = l + 1
!        end do
!     end do
!     do i = 4, 9
!        call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
!        call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!        d = 0.
!        do j = 1, nmo
!           call cp_cfm_get_element(mom, j, j, cd)
!           d = d + zwfc * real(cd)
!        end do
!        moments(i+3) = moments(i+3) - d
!     end do
!
!     do i = 1, size(o)
!        call cp_cfm_release(o(i)%matrix)
!     end do
!     do i = 1, size(mago)
!        call cp_cfm_release(mago(i)%matrix)
!     end do
!     call cp_cfm_release(m)
!     call cp_cfm_release(om)
!     call cp_cfm_release(mom)
!     call cp_fm_release(or)
!     deallocate(mago, o)
!     moments = 0._dp
!     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, 1, &
!          size(particles), rcc)
!     call calculate_localized_loc_moments(moments, qs_env, rcc, mo_local, &
!          particles, qs_kinds, cell, 1, size(particles), ref=.true.)
!
!   end subroutine calculate_reference_moments


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
!        last_atom, rcc)
!
!      real(dp)                      :: moments(:)
!      type(particle_type), pointer  :: particles(:)
!      type(qs_kind_type), pointer   :: qs_kinds(:)
!      type(cell_type), pointer      :: cell
!      integer                       :: first_atom, last_atom
!      real(dp)                      :: rcc(3)
!
!      logical                          :: do_periodic
!      integer                          :: akind, iatom, l, lx, ly, lz
!      logical                          :: ghost, floating, do_ref
!      real(dp)                         :: d, r(3), z
!      type(atomic_kind_type), pointer  :: atomic_kind
!      
!      do iatom = first_atom, last_atom
!         r = pbc(particles(iatom)%r - rcc, cell) + rcc
!         r = r - rcc
!         atomic_kind => particles(iatom)%atomic_kind
!         call get_atomic_kind(atomic_kind, kind_number=akind)
!         call get_qs_kind(qs_kinds(akind), core_charge=z, ghost=ghost, floating=floating)
!         if (.not. ghost .and. .not. floating) then
!            r = pbc(particles(iatom)%r - rcc, cell)
!            moments(1) = moments(1) + z  ! nuclear charge
!            do l = 1, 9
!               lx = indco(1, l+1)
!               ly = indco(2, l+1)
!               lz = indco(3, l+1)
!               d = 1._dp
!               if (lx > 0) d = d * r(1)**lx
!               if (ly > 0) d = d * r(2)**ly
!               if (lz > 0) d = d * r(3)**lz
!               moments(l+3) = moments(l+3) + z * d
!            end do
!         end if
!      end do
!      moments(3) = moments(3) + moments(2) ! total charge
!
!   end subroutine calculate_localized_nuc_moments


! **************************************************************************************************
! **************************************************************************************************
   subroutine calculate_localized_loc_moments(moments, qs_env, rcc, mo_local, particles, qs_kinds, &
        cell, first_atom, last_atom, statelist, periodic, ref)

     real(dp)                            :: moments(:)
     type(qs_environment_type), pointer  :: qs_env
     real(dp)                            :: rcc(3)
     type(cp_fm_p_type), pointer         :: mo_local(:)
     type(particle_type), pointer        :: particles(:)
     type(qs_kind_type), pointer         :: qs_kinds(:)
     type(cell_type), pointer            :: cell
     integer                             :: first_atom, last_atom
     integer, pointer, optional          :: statelist(:)
     logical, optional                   :: periodic, ref

     complex(dp)                              :: cd
     integer                                  :: i, j, k, nao, nmm, nmo, ns
     integer, pointer                         :: states(:)
     logical                                  :: do_rtp, do_periodic, do_ref
     real(dp)                                 :: d, zwfc, qupoleder(9)
     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
     type(cp_cfm_p_type), pointer             :: o(:), oder(:, :), mago(:), momo(:)
     type(cp_cfm_type), pointer               :: m, om, mom
     type(dbcsr_p_type), pointer              :: matrix_s(:), op(:), opder(:, :), magop(:), momop(:)
     type(cp_fm_struct_type), pointer         :: m_struct, o_struct, mom_struct
     type(cp_fm_type), pointer                :: mr, mi, or
     type(neighbor_list_set_p_type), pointer  :: sab_all(:), sab_orb(:)

     !print *, "dbg| loc"
     !print *, "dbg| calc_loc_mom rcc", rcc
      
     ! info
     if (present(periodic)) then
        do_periodic = periodic
     else
        do_periodic = .false.
     end if
     if (present(ref)) then
        do_ref = ref
     else
        do_ref = .false.
     end if
     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
     do_rtp = (size(mo_local) == 2)
     call cp_fm_get_info(mo_local(1)%matrix, nrow_global=nao, ncol_global=nmo)
     zwfc = 2._dp
     allocate(states(nmo))
     if (present(statelist)) then
        ns = size(statelist)
        do i = 1, ns
           states(i) = statelist(i)
        end do
     else
        ns = nmo
        do i = 1, ns
           states(i) = i
        end do
     end if

     ! structs
     call cp_fm_struct_create(m_struct, nrow_global=nao, ncol_global=ns, &
          template_fmstruct=mo_local(1)%matrix%matrix_struct)
     call cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
          template_fmstruct=mo_local(1)%matrix%matrix_struct)
     call cp_fm_struct_create(mom_struct, nrow_global=ns, ncol_global=ns, &
          template_fmstruct=mo_local(1)%matrix%matrix_struct)

     ! get only this subsystem's states
     call cp_fm_create(mr, m_struct, name="mr")
     call cp_fm_set_all(mr, 0._dp)
     if (do_rtp) then
        call cp_fm_create(mi, m_struct, name="mi")
        call cp_fm_set_all(mi, 0._dp)
     end if
     call cp_cfm_create(m, m_struct, name="m")
     do i = 1, ns
        call cp_fm_to_fm(mo_local(1)%matrix, mr, 1, states(i), i)
        if (do_rtp) then
           call cp_fm_to_fm(mo_local(2)%matrix, mi, 1, states(i), i)
           call cp_fm_to_cfm(mr, mi, m)
        else
           call cp_fm_to_cfm(mr, mtarget=m)
        end if
     end do
     call cp_fm_release(mr)
     if (do_rtp) then
        call cp_fm_release(mi)
     end if

     !! get the operators
     call cp_fm_create(or, o_struct, name="or")
     ! el moments with derivatives
     nullify(op)
     call dbcsr_allocate_matrix_set(op, nm)
     do i = 1, nm
        call dbcsr_init_p(op(i)%matrix)
        call dbcsr_create(op(i)%matrix, template=matrix_s(1)%matrix, &
             matrix_type=dbcsr_type_no_symmetry)
        call cp_dbcsr_alloc_block_from_nbl(op(i)%matrix, sab_all)
        call dbcsr_set(op(i)%matrix, 0._dp)
     end do
     nullify(opder)
     call dbcsr_allocate_matrix_set(opder, 3, 3)
     do i = 1, 3
        do j = 1, 3
           call dbcsr_init_p(opder(i, j)%matrix)
           call dbcsr_create(opder(i, j)%matrix, template=matrix_s(1)%matrix, &
                matrix_type=dbcsr_type_no_symmetry)
           call cp_dbcsr_alloc_block_from_nbl(opder(i, j)%matrix, sab_all)
           call dbcsr_set(opder(i, j)%matrix, 0._dp)
        end do
     end do
     allocate(o(9))
     if (ref) then
        call build_local_moments_der_matrix(qs_env, op, opder, 2, 1, rcc)
     else
        call build_local_moment_matrix(qs_env, op, 2, rcc)
     end if
     do i = 1, nm
        call cp_cfm_create(o(i)%matrix, o_struct)
        if (ref) then
           call copy_dbcsr_to_fm(op(i)%matrix, or)
        end if
        call dbcsr_release(op(i)%matrix)
        deallocate(op(i)%matrix)
        if (ref) then
           call cp_fm_to_cfm(or, mtarget=o(i)%matrix)
        else
           call cp_cfm_set_all(o(i)%matrix, (0._dp, 0._dp))
        end if
     end do
     deallocate(op)
     allocate(oder(3, 3))
     do i = 1, 3
        do j = 1, 3
           call cp_cfm_create(oder(i, j)%matrix, o_struct)
           if (ref) then
              call copy_dbcsr_to_fm(opder(i, j)%matrix, or)
           end if
           call dbcsr_release(opder(i, j)%matrix)
           deallocate(opder(i, j)%matrix)
           if (ref) then
              call cp_fm_to_cfm(or, mtarget=oder(i, j)%matrix)
           else
              call cp_cfm_set_all(oder(i, j)%matrix, (0._dp, 0._dp))
           end if
        end do
     end do
     deallocate(opder)

     ! momentum
     nullify(momop)
     call dbcsr_allocate_matrix_set(momop, 3)
     do i = 1, 3
        call dbcsr_init_p(momop(i)%matrix)
        call dbcsr_create(momop(i)%matrix, template=matrix_s(1)%matrix, &
             matrix_type=dbcsr_type_antisymmetric)
        call cp_dbcsr_alloc_block_from_nbl(momop(i)%matrix, sab_orb)
        call dbcsr_set(momop(i)%matrix, 0._dp)
     end do
     allocate(momo(3))
     call build_lin_mom_matrix(qs_env, momop)
     do i = 1, 3
        call cp_cfm_create(momo(i)%matrix, o_struct)
        call copy_dbcsr_to_fm(momop(i)%matrix, or)
        call dbcsr_release(momop(i)%matrix)
        deallocate(momop(i)%matrix)
        call cp_fm_to_cfm(or, mtarget=momo(i)%matrix)
     end do
     deallocate(momop)

     ! mag moments
     nullify(magop)
     call dbcsr_allocate_matrix_set(magop, 3)
     do i = 1, 3
        call dbcsr_init_p(magop(i)%matrix)
        call dbcsr_create(magop(i)%matrix, template=matrix_s(1)%matrix, &
             matrix_type=dbcsr_type_antisymmetric)
        call cp_dbcsr_alloc_block_from_nbl(magop(i)%matrix, sab_orb)
        call dbcsr_set(magop(i)%matrix, 0._dp)
     end do
     nullify(op)
     allocate(mago(3))
     call build_local_magmom_matrix(qs_env, magop, 1, ref_point=rcc)
     nmm = size(magop)
     do i = 1, nmm
        call cp_cfm_create(mago(i)%matrix, o_struct, name="mago")
        if (ref) then
           call copy_dbcsr_to_fm(magop(i)%matrix, or)
        end if
        call dbcsr_release(magop(i)%matrix)
        deallocate(magop(i)%matrix)
        if (ref) then
           call cp_fm_to_cfm(or, mtarget=mago(i)%matrix)
        else
           call cp_cfm_set_all(mago(i)%matrix, (0._dp, 0._dp))
        end if
     end do
     deallocate(magop)
     call cp_fm_release(or)

     ! get the nonlocal terms
     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
     nlrv = 0._dp
     nlrrv = 0._dp
     nlrxrv = 0._dp
     if (ref) then
        call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, rcc)
     end if

     ! nuclear contribution
     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
          last_atom, rcc)

     ! expectation values
     call cp_cfm_create(om, m_struct, name="om")
     call cp_cfm_create(mom, mom_struct, name="mom")
     call cp_fm_struct_release(m_struct)
     call cp_fm_struct_release(o_struct)
     call cp_fm_struct_release(mom_struct)

     ! charge
     !call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, m, (0._dp, 0._dp), mom)
     !d = 0.
     !do i = 1, ns
     !   call cp_cfm_get_element(mom, i, i, cd)
     !   moments(2) = moments(2) + real(cd)
     !end do
     !moments(2) = -zwfc * ns! moments(2)
     !moments(3) = moments(1) + moments(2)

     ! edip
     do i = 1, 3
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
        d = 0.
        do j = 1, ns
           call cp_cfm_get_element(mom, j, j, cd)
           d = d + real(cd)
        end do
        moments(i+3) = moments(i+3) - zwfc * d
     end do
     ! vel
     !do i = 1, 3
     !   call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), momo(i)%matrix, m, (0._dp, 0._dp), om)
     !   call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
     !   d = 0.
     !   do j = 1, ns
     !      call cp_cfm_get_element(mom, j, j, cd)
     !      d = d + aimag(cd)
     !   end do
     !   moments(i+15) = moments(i+15) - zwfc * d
     !end do
     !! nl
     !moments(28:30) = moments(16:18) + nlrv

     !! mdip
     !do i = 1, 3
     !   call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), mago(i)%matrix, m, (0._dp, 0._dp), om)
     !   call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
     !   d = 0.
     !   do j = 1, ns
     !      call cp_cfm_get_element(mom, j, j, cd)
     !      d = d + aimag(cd)
     !   end do
     !   moments(i+12) = moments(i+12) - zwfc * d
     !end do
     !! nl
     !moments(25:27) = moments(13:15) + nlrxrv

     !! equp
     !do i = 4, 9
     !   call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
     !   call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
     !   d = 0.
     !   do j = 1, ns
     !      call cp_cfm_get_element(mom, j, j, cd)
     !      d = d + zwfc * real(cd)
     !   end do
     !   moments(i+3) = moments(i+3) - d
     !end do
     !! vel
     !qupoleder = 0._dp
     !do i = 1, 3
     !   do j = 1, 3
     !      call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), oder(i, j)%matrix, m, (0._dp, 0._dp), om)
     !      call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
     !      d = 0.
     !      do k = 1, ns
     !         call cp_cfm_get_element(mom, k, k, cd)
     !         d = d - aimag(cd)
     !         if (i == j) then
     !            d = d + real(cd)
     !         end if
     !      end do
     !      qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + zwfc * d
     !   end do
     !end do
     !moments(19) = -2*qupoleder(1) - moments(2)
     !moments(20) = -qupoleder(2) - qupoleder(4)
     !moments(21) = -qupoleder(3) - qupoleder(7)
     !moments(22) = -2*qupoleder(5) - moments(2)
     !moments(23) = -qupoleder(6) - qupoleder(8)
     !moments(24) = -2*qupoleder(9) - moments(2)
     !! nl
     !moments(31:36) = moments(19:24) + nlrrv

     ! cleanup
     call cp_cfm_release(m)
     call cp_cfm_release(om)
     call cp_cfm_release(mom)
     do i = 1, size(o)
        call cp_cfm_release(o(i)%matrix)
     end do
     deallocate(o)
     do i = 1, size(oder, 1)
        do j = 1, size(oder, 2)
           call cp_cfm_release(oder(i, j)%matrix)
        end do
     end do
     deallocate(oder)
     do i = 1, size(momo)
        call cp_cfm_release(momo(i)%matrix)
     end do
     deallocate(momo)
     do i = 1, size(mago)
        call cp_cfm_release(mago(i)%matrix)
     end do
     deallocate(mago)
     deallocate(nlrv, nlrrv, nlrxrv, states)
     
   end subroutine calculate_localized_loc_moments


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_localized_kg_moments(moments, qs_env, rcc, rho, particles, qs_kinds, cell, &
!        first_atom, last_atom, do_periodic)
!
!     real(dp)                            :: moments(:)
!     type(qs_environment_type), pointer  :: qs_env
!     real(dp)                            :: rcc(3)
!     type(qs_rho_type), pointer          :: rho
!     type(particle_type), pointer        :: particles(:)
!     type(qs_kind_type), pointer         :: qs_kinds(:)
!     type(cell_type), pointer            :: cell
!     integer                             :: first_atom, last_atom
!     logical  :: do_periodic
!
!     integer                                  :: i, j
!     logical                                  :: found, tr
!     real(dp)                                 :: focc, qupoleder(9), strace, trace
!     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
!     real(dp), pointer                        :: blk_p(:)
!     type(dbcsr_iterator_type), pointer       :: iterator
!     type(dbcsr_p_type), pointer              :: o(:), oder(:, :), mago(:), momo(:), matrix_s(:), &
!          rho_ao(:), rho_ao_im(:)
!     type(dbcsr_type), pointer                :: densmat, densmat_im, tmp
!     type(neighbor_list_set_p_type), pointer  :: sab_all(:), sab_orb(:)
!
!     ! info
!     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
!     call qs_rho_get(rho, rho_ao=rho_ao, rho_ao_im=rho_ao_im)
!     focc = 1._dp
!
!     ! get real and imaginary part of the subdensity
!     allocate(densmat)
!     call dbcsr_copy(densmat, rho_ao(1)%matrix)
!     allocate(iterator)
!     call dbcsr_iterator_start(iterator, densmat)
!     do while (dbcsr_iterator_blocks_left(iterator))
!        call dbcsr_iterator_next_block(iterator, i, j, blk_p)
!        if ((i < first_atom) .or. (i > last_atom) .or. (j < first_atom) .or. (j > last_atom)) then
!           call dbcsr_get_block_p(densmat, i, j, blk_p, tr, found)
!           blk_p = 0._dp
!        end if
!     end do
!     call dbcsr_iterator_stop(iterator)
!     allocate(densmat_im)
!     call dbcsr_copy(densmat_im, rho_ao_im(1)%matrix)
!     call dbcsr_iterator_start(iterator, densmat_im)
!     do while (dbcsr_iterator_blocks_left(iterator))
!        call dbcsr_iterator_next_block(iterator, i, j, blk_p)
!        if ((i < first_atom) .or. (i > last_atom) .or. (j < first_atom) .or. (j > last_atom)) then
!           call dbcsr_get_block_p(densmat_im, i, j, blk_p, tr, found)
!           blk_p = 0._dp
!        end if
!     end do
!     call dbcsr_iterator_stop(iterator)
!     deallocate(iterator)
!
!     ! get the operators
!     ! el moments with derivatives
!     nullify(o)
!     call dbcsr_allocate_matrix_set(o, nm)
!     do i = 1, nm
!        call dbcsr_init_p(o(i)%matrix)
!        call dbcsr_create(o(i)%matrix, template=matrix_s(1)%matrix, &
!             matrix_type=dbcsr_type_no_symmetry)
!        call cp_dbcsr_alloc_block_from_nbl(o(i)%matrix, sab_all)
!        call dbcsr_set(o(i)%matrix, 0._dp)
!     end do
!     nullify(oder)
!     call dbcsr_allocate_matrix_set(oder, 3, 3)
!     do i = 1, 3
!        do j = 1, 3
!           call dbcsr_init_p(oder(i, j)%matrix)
!           call dbcsr_create(oder(i, j)%matrix, template=matrix_s(1)%matrix, &
!                matrix_type=dbcsr_type_no_symmetry)
!           call cp_dbcsr_alloc_block_from_nbl(oder(i, j)%matrix, sab_all)
!           call dbcsr_set(oder(i, j)%matrix, 0._dp)
!        end do
!     end do
!     call build_local_moments_der_matrix(qs_env, o, oder, 2, 1, rcc)
!
!     ! momentum
!     nullify(momo)
!     call dbcsr_allocate_matrix_set(momo, 3)
!     do i = 1, 3
!        call dbcsr_init_p(momo(i)%matrix)
!        call dbcsr_create(momo(i)%matrix, template=matrix_s(1)%matrix, &
!             matrix_type=dbcsr_type_antisymmetric)
!        call cp_dbcsr_alloc_block_from_nbl(momo(i)%matrix, sab_orb)
!        call dbcsr_set(momo(i)%matrix, 0._dp)
!     end do
!     call build_lin_mom_matrix(qs_env, momo)
!
!     ! mag moments
!     nullify(mago)
!     call dbcsr_allocate_matrix_set(mago, 3)
!     do i = 1, 3
!        call dbcsr_init_p(mago(i)%matrix)
!        call dbcsr_create(mago(i)%matrix, template=matrix_s(1)%matrix, &
!             matrix_type=dbcsr_type_antisymmetric)
!        call cp_dbcsr_alloc_block_from_nbl(mago(i)%matrix, sab_orb)
!        call dbcsr_set(mago(i)%matrix, 0._dp)
!     end do
!     call build_local_magmom_matrix(qs_env, mago, nmm, ref_point=rcc)
!
!     ! get the nonlocal terms
!     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
!     nlrv = 0._dp
!     nlrrv = 0._dp
!     nlrxrv = 0._dp
!     call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, rcc)
!     
!     ! nuclear contribution
!     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
!          last_atom, rcc)
!
!     ! temp matrix
!     allocate(tmp)
!     call dbcsr_desymmetrize(matrix_s(1)%matrix, tmp)
!     
!     ! charge
!     call dbcsr_dot(densmat, matrix_s(1)%matrix, trace)
!     moments(2) = - focc * trace
!     moments(3) = moments(1) + moments(2)
!
!     ! edip
!     do i = 1, 3
!        call dbcsr_multiply("T", "N", 1._dp, densmat, o(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+3) = moments(i+3) - focc * trace
!     end do
!     ! vel
!     do i = 1, 3
!        call dbcsr_multiply("T", "N", 1._dp, densmat_im, momo(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+15) = moments(i+15) + focc * trace
!     end do
!     ! nl
!     moments(28:30) = moments(16:18) + nlrv
!
!     ! mdip
!     do i = 1, 3
!        call dbcsr_multiply("T", "N", 1._dp, densmat_im, mago(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+12) = moments(i+12) + focc * trace
!     end do
!     ! nl
!     moments(25:27) = moments(13:15) + nlrxrv
!
!     ! equp
!     do i = 4, 9
!        call dbcsr_multiply("T", "N", 1._dp, densmat, o(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+3) = moments(i+3) - focc * trace
!     end do
!     ! vel
!     qupoleder = 0._dp
!     do i = 1, 3
!        do j = 1, 3
!           call dbcsr_multiply("T", "N", 1._dp, densmat, oder(i, j)%matrix, 0._dp, tmp)
!           call dbcsr_trace(tmp, trace)
!           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + focc * trace
!           call dbcsr_multiply("T", "N", 1._dp, densmat_im, oder(i, j)%matrix, 0._dp, tmp)
!           call dbcsr_trace(tmp, trace)
!           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + focc * trace
!        end do
!     end do
!     moments(19) = -2*qupoleder(1) - moments(2)
!     moments(20) = -qupoleder(2) - qupoleder(4)
!     moments(21) = -qupoleder(3) - qupoleder(7)
!     moments(22) = -2*qupoleder(5) - moments(2)
!     moments(23) = -qupoleder(6) - qupoleder(8)
!     moments(24) = -2*qupoleder(9) - moments(2)
!     ! nl
!     moments(31:36) = moments(19:24) + nlrrv
!
!     do i = 1, size(o)
!        call dbcsr_release(o(i)%matrix)
!        deallocate(o(i)%matrix)
!     end do
!     deallocate(o)
!     do i = 1, size(oder, 1)
!        do j = 1, size(oder, 2)
!           call dbcsr_release(oder(i, j)%matrix)
!           deallocate(oder(i, j)%matrix)
!        end do
!     end do
!     deallocate(oder)
!     do i = 1, size(momo)
!        call dbcsr_release(momo(i)%matrix)
!        deallocate(momo(i)%matrix)
!     end do
!     deallocate(momo)
!     do i = 1, size(mago)
!        call dbcsr_release(mago(i)%matrix)
!        deallocate(mago(i)%matrix)
!     end do
!     deallocate(mago)
!     call dbcsr_release(densmat)
!     deallocate(densmat)
!     call dbcsr_release(densmat_im)
!     deallocate(densmat_im)
!     call dbcsr_release(tmp)
!     deallocate(tmp)
!     deallocate(nlrv, nlrrv, nlrxrv)
!
!   end subroutine calculate_localized_kg_moments


! **************************************************************************************************
! **************************************************************************************************
   subroutine print_localized_moments(logger, input, moments, rcc, molname, suffix)

      type(cp_logger_type), pointer     :: logger
      type(section_vals_type), pointer  :: input
      real(dp)                          :: moments(:), rcc(:)
      character(*)                      :: molname
      character(*), optional            :: suffix

      integer                               :: output_unit
      character(len=default_string_length)  :: iter

      if (present(suffix)) then
         output_unit = cp_print_key_unit_nr(logger=logger, basis_section=input, &
              print_key_path="DFT%PRINT%MOMENTS", extension=".csv", &
              middle_name="locmoments_"//trim(molname)//"_"//trim(suffix), log_filename=.false.)
      else
         output_unit = cp_print_key_unit_nr(logger=logger, basis_section=input, &
              print_key_path="DFT%PRINT%MOMENTS", extension=".csv", &
              middle_name="locmoments_"//trim(molname), log_filename=.false.)
      end if
      iter = cp_iter_string(logger%iter_info)
      
      if (output_unit > 0) then
         if (iteration == 1) then
            write (output_unit, "(A)") header
         end if
         write (output_unit, "(I6,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3)") &
              iteration,   ", "//trim(iter)//", ", &
              moments( 1), ", ", moments( 2), ", ", moments( 3), ", ", &
              moments( 4), ", ", moments( 5), ", ", moments( 6), ", ", moments( 7), ", ", &
              moments( 8), ", ", moments( 9), ", ", moments(10), ", ", moments(11), ", ", &
              moments(12), ", ", moments(13), ", ", moments(14), ", ", moments(15), ", ", &
              moments(16), ", ", moments(17), ", ", moments(18), ", ", moments(19), ", ", &
              moments(20), ", ", moments(21), ", ", moments(22), ", ", moments(23), ", ", &
              moments(24), ", ", moments(25), ", ", moments(26), ", ", moments(27), ", ", &
              moments(28), ", ", moments(29), ", ", moments(30), ", ", moments(31), ", ", &
              moments(32), ", ", moments(33), ", ", moments(34), ", ", moments(35), ", ", &
              moments(36), ", ", rcc(1), ", ", rcc(2), ",", rcc(3)
      end if

      call cp_print_key_finished_output(output_unit, logger, input, "DFT%PRINT%MOMENTS")

   end subroutine print_localized_moments

end module localized_moments
 
