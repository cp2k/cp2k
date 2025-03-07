!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM nequip_unittest

   USE cp_files,                        ONLY: discover_file
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8,&
                                              sp
   USE mathlib,                         ONLY: inv_3x3
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE torch_api,                       ONLY: &
        torch_cuda_is_available, torch_dict_create, torch_dict_get, torch_dict_insert, &
        torch_dict_release, torch_dict_type, torch_model_forward, torch_model_load, &
        torch_model_read_metadata, torch_model_release, torch_model_type, torch_tensor_data_ptr, &
        torch_tensor_from_array, torch_tensor_release, torch_tensor_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(LEN=default_path_length) :: filename, cutoff_str, nequip_version
   REAL(dp) :: cutoff

   ! Inputs.
   INTEGER, PARAMETER  :: natoms = 96
   INTEGER :: iatom, nedges
   REAL(sp), DIMENSION(:, :), ALLOCATABLE :: pos, cell
   REAL(dp), DIMENSION(3, 3) :: hinv
   INTEGER(kind=int_8), DIMENSION(:), ALLOCATABLE :: atom_types
   INTEGER(kind=int_8), DIMENSION(:, :), ALLOCATABLE :: edge_index
   REAL(sp), DIMENSION(:, :), ALLOCATABLE:: edge_cell_shift

   ! Torch objects.
   TYPE(torch_model_type) :: model
   TYPE(torch_dict_type) :: inputs, outputs
   TYPE(torch_tensor_type) :: pos_tensor, edge_index_tensor, edge_cell_shift_tensor, cell_tensor, &
                              atom_types_tensor, total_energy_tensor, atomic_energy_tensor, forces_tensor

   ! Outputs.
   REAL(sp), DIMENSION(:, :), POINTER :: total_energy, atomic_energy, forces
   NULLIFY (total_energy, atomic_energy, forces)

   ! A box with 32 water molecules.
   ALLOCATE (pos(3, natoms))
   pos(:, :) = RESHAPE(REAL([ &
                            42.8861696_dp, -0.0556816_dp, 38.3291611_dp, &
                            34.2025887_dp, -0.6185484_dp, 37.3655680_dp, &
                            30.0803925_dp, -2.0124176_dp, 36.4807960_dp, &
                            28.7057911_dp, -2.6880392_dp, 36.6020983_dp, &
                            36.2479426_dp, -0.5163484_dp, 34.4923596_dp, &
                            37.6964724_dp, -0.0410872_dp, 35.0140735_dp, &
                            27.7606699_dp, 7.4854206_dp, 33.9276919_dp, &
                            28.8160999_dp, 6.4985777_dp, 34.2163608_dp, &
                            37.1576372_dp, 9.0188280_dp, 31.9265812_dp, &
                            38.6063816_dp, 9.5820079_dp, 32.3435972_dp, &
                            34.3031959_dp, 2.2195014_dp, 45.9880451_dp, &
                            33.2444139_dp, 1.3025332_dp, 46.4698427_dp, &
                            38.7286174_dp, -5.0541897_dp, 26.0743968_dp, &
                            38.3483921_dp, -6.2832846_dp, 26.9867253_dp, &
                            32.8642520_dp, 3.2060632_dp, 30.8971160_dp, &
                            31.2904088_dp, 3.0871834_dp, 30.6273977_dp, &
                            33.7519869_dp, -3.1383262_dp, 39.6727607_dp, &
                            34.6642979_dp, -3.6643859_dp, 38.6466027_dp, &
                            42.7173214_dp, 5.1246883_dp, 32.5883401_dp, &
                            41.5627455_dp, 5.5893544_dp, 33.4174902_dp, &
                            32.4283800_dp, 9.1182520_dp, 30.5477678_dp, &
                            32.6432407_dp, 10.770683_dp, 30.4842778_dp, &
                            31.4848670_dp, 4.6777144_dp, 37.3957194_dp, &
                            32.3171882_dp, -6.2287496_dp, 36.4671864_dp, &
                            26.6621340_dp, 3.1708123_dp, 35.6820146_dp, &
                            26.5271367_dp, 1.6039040_dp, 35.4883482_dp, &
                            32.0238236_dp, 16.918208_dp, 31.6883569_dp, &
                            31.4006579_dp, 7.0315610_dp, 30.2394554_dp, &
                            33.5264253_dp, -3.5594808_dp, 34.2636830_dp, &
                            34.6404855_dp, -3.2653833_dp, 35.4971482_dp, &
                            40.0564375_dp, -0.3054386_dp, 29.8312074_dp, &
                            39.4784464_dp, -1.0948314_dp, 38.3101140_dp, &
                            39.7040761_dp, 1.9584631_dp, 33.3902375_dp, &
                            38.3338570_dp, 2.6967178_dp, 42.9261945_dp, &
                            40.1820455_dp, -7.2199289_dp, 27.6580390_dp, &
                            39.3204431_dp, -8.4564252_dp, 28.1319658_dp, &
                            36.3876963_dp, 8.8117085_dp, 38.3545362_dp, &
                            36.3205637_dp, 9.0063075_dp, 36.7526001_dp, &
                            29.9991583_dp, -5.5637817_dp, 33.9295050_dp, &
                            30.7728545_dp, -5.0385870_dp, 35.1998067_dp, &
                            40.0592517_dp, 6.3305279_dp, 28.2579461_dp, &
                            40.2398360_dp, 5.1745923_dp, 29.2962956_dp, &
                            26.3320911_dp, 2.4393638_dp, 33.5653868_dp, &
                            26.9606971_dp, 1.2711078_dp, 32.5923884_dp, &
                            34.8372697_dp, -0.4722708_dp, 30.3824362_dp, &
                            35.3968813_dp, -1.9268483_dp, 30.3081837_dp, &
                            32.1217607_dp, -0.7333429_dp, 36.5104382_dp, &
                            32.2180843_dp, 7.8454304_dp, 35.6671967_dp, &
                            36.3780998_dp, -4.3048878_dp, 36.4539793_dp, &
                            35.8119275_dp, -3.0013928_dp, 27.0348937_dp, &
                            29.6452491_dp, 1.0652123_dp, 35.7143653_dp, &
                            30.3794654_dp, -0.0668146_dp, 34.9882468_dp, &
                            34.2149336_dp, -1.6559120_dp, 33.8876437_dp, &
                            34.7842435_dp, -1.0252141_dp, 32.5034832_dp, &
                            40.4649954_dp, 1.1467825_dp, 31.3073503_dp, &
                            41.3262469_dp, 0.6550803_dp, 32.4555882_dp, &
                            29.0210859_dp, 3.5038194_dp, 39.9087702_dp, &
                            29.4945426_dp, 3.7276637_dp, 41.3766138_dp, &
                            34.1359664_dp, -6.7533422_dp, 32.3568410_dp, &
                            34.9546570_dp, -5.7704242_dp, 31.4571066_dp, &
                            33.2532356_dp, 1.5268048_dp, 44.0562171_dp, &
                            33.7931669_dp, 0.5014632_dp, 43.0597590_dp, &
                            36.8205409_dp, 2.6214681_dp, 40.6834006_dp, &
                            37.5552706_dp, 1.5649832_dp, 39.7648935_dp, &
                            43.2099087_dp, -0.0628456_dp, 47.2593155_dp, &
                            29.3940583_dp, -2.3133019_dp, 37.1407883_dp, &
                            36.7415708_dp, -0.0838710_dp, 35.2591783_dp, &
                            27.9424776_dp, 6.7622961_dp, 34.5648384_dp, &
                            37.6812656_dp, 9.4216399_dp, 32.6478643_dp, &
                            33.3171290_dp, 2.0951401_dp, 45.8722265_dp, &
                            37.9951355_dp, 4.3611431_dp, 26.5571819_dp, &
                            32.1824670_dp, 2.6611503_dp, 30.4577248_dp, &
                            34.6538012_dp, -3.4374573_dp, 39.5889245_dp, &
                            42.2929833_dp, 5.9471069_dp, 32.8460995_dp, &
                            32.9604690_dp, 9.9050313_dp, 30.1587306_dp, &
                            31.4281886_dp, -5.8338304_dp, 36.6738743_dp, &
                            26.0563730_dp, 2.4973869_dp, 35.3486870_dp, &
                            32.0334927_dp, 17.3252289_dp, 30.8116013_dp, &
                            33.8252182_dp, -2.9520949_dp, 35.0220460_dp, &
                            39.4569981_dp, -0.3072759_dp, 38.9347829_dp, &
                            29.4846708_dp, 2.8692561_dp, 43.0061868_dp, &
                            39.2864184_dp, -7.6206103_dp, 27.6271147_dp, &
                            35.8797502_dp, 8.6515870_dp, 37.5221734_dp, &
                            30.3582543_dp, -4.7607656_dp, 34.3355645_dp, &
                            40.7098956_dp, 5.8331250_dp, 28.7558375_dp, &
                            26.7179083_dp, 2.2415138_dp, 32.6577297_dp, &
                            35.6589256_dp, -0.9968903_dp, 30.5749530_dp, &
                            31.5851602_dp, -1.3121804_dp, 35.9011109_dp, &
                            35.5489386_dp, -3.9056138_dp, 26.8214490_dp, &
                            29.5656616_dp, 0.4681794_dp, 34.9670711_dp, &
                            34.7615128_dp, -0.9569680_dp, 33.4891367_dp, &
                            40.4853406_dp, 0.4023620_dp, 31.9425416_dp, &
                            29.6728289_dp, 4.0134825_dp, 40.4505780_dp, &
                            34.1272286_dp, -5.8796882_dp, 31.8925999_dp, &
                            33.1168884_dp, 1.2338084_dp, 43.1127997_dp, &
                            37.1996993_dp, 2.5049007_dp, 39.7917126_dp], kind=sp), shape=[3, natoms])

   ALLOCATE (cell(3, 3))
   cell(1, :) = [9.85_sp, 0.0_sp, 0.0_sp]
   cell(2, :) = [0.0_sp, 9.85_sp, 0.0_sp]
   cell(3, :) = [0.0_sp, 0.0_sp, 9.85_sp]

   hinv(:, :) = inv_3x3(REAL(cell, kind=dp))

   ALLOCATE (atom_types(natoms))
   atom_types(:64) = 0 ! Hydrogen
   atom_types(65:) = 1 ! Oxygen

   WRITE (*, *) "CUDA is available: ", torch_cuda_is_available()

   filename = discover_file('NequIP/water-deployed-neq060sp.pth')
   WRITE (*, *) "Loading NequIP model from: "//TRIM(filename)
   CALL torch_model_load(model, filename)
   cutoff_str = torch_model_read_metadata(filename, "r_max")
   nequip_version = torch_model_read_metadata(filename, "nequip_version")
   READ (cutoff_str, *) cutoff
   WRITE (*, *) "Version: ", TRIM(nequip_version)
   WRITE (*, *) "Cutoff: ", cutoff

   CALL neighbor_search(nedges)
   ALLOCATE (edge_index(nedges, 2))
   ALLOCATE (edge_cell_shift(3, nedges))
   CALL neighbor_search(nedges, edge_index, edge_cell_shift)
   WRITE (*, *) "Found", nedges, "neighbor edges between", natoms, "atoms."

   CALL torch_dict_create(inputs)
   CALL torch_dict_create(outputs)

   CALL torch_tensor_from_array(pos_tensor, pos)
   CALL torch_dict_insert(inputs, "pos", pos_tensor)
   CALL torch_tensor_release(pos_tensor)

   CALL torch_tensor_from_array(edge_index_tensor, edge_index)
   CALL torch_dict_insert(inputs, "edge_index", edge_index_tensor)
   CALL torch_tensor_release(edge_index_tensor)

   CALL torch_tensor_from_array(edge_cell_shift_tensor, edge_cell_shift)
   CALL torch_dict_insert(inputs, "edge_cell_shift", edge_cell_shift_tensor)
   CALL torch_tensor_release(edge_cell_shift_tensor)

   CALL torch_tensor_from_array(cell_tensor, cell)
   CALL torch_dict_insert(inputs, "cell", cell_tensor)
   CALL torch_tensor_release(cell_tensor)

   CALL torch_tensor_from_array(atom_types_tensor, atom_types)
   CALL torch_dict_insert(inputs, "atom_types", atom_types_tensor)
   CALL torch_tensor_release(atom_types_tensor)

   CALL torch_model_forward(model, inputs, outputs)

   CALL torch_dict_get(outputs, "total_energy", total_energy_tensor)
   CALL torch_tensor_data_ptr(total_energy_tensor, total_energy)

   CALL torch_dict_get(outputs, "atomic_energy", atomic_energy_tensor)
   CALL torch_tensor_data_ptr(atomic_energy_tensor, atomic_energy)

   CALL torch_dict_get(outputs, "forces", forces_tensor)
   CALL torch_tensor_data_ptr(forces_tensor, forces)

   WRITE (*, *) "Total Energy [Hartree] : ", total_energy(1, 1)/evolt
   WRITE (*, *) "FORCES: [Hartree/Bohr]: "
   DO iatom = 1, natoms
      WRITE (*, *) forces(:, iatom)*angstrom/evolt
   END DO
   CPASSERT(ABS(-14985.4443_dp - REAL(total_energy(1, 1), kind=dp)) < 2e-3_dp)

   CALL torch_tensor_release(total_energy_tensor)
   CALL torch_tensor_release(atomic_energy_tensor)
   CALL torch_tensor_release(forces_tensor)
   CALL torch_dict_release(inputs)
   CALL torch_dict_release(outputs)
   CALL torch_model_release(model)
   DEALLOCATE (edge_index, edge_cell_shift, pos, cell, atom_types)

   WRITE (*, *) "NequIP unittest was successfully :-)"

CONTAINS

! **************************************************************************************************
!> \brief Naive neighbor search - beware it scales O(N**2).
!> \param nedges ...
!> \param edge_index ...
!> \param edge_cell_shift ...
! **************************************************************************************************
   SUBROUTINE neighbor_search(nedges, edge_index, edge_cell_shift)
      INTEGER, INTENT(OUT)                               :: nedges
      INTEGER(kind=int_8), DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: edge_index
      REAL(sp), DIMENSION(:, :), INTENT(OUT), OPTIONAL   :: edge_cell_shift

      INTEGER:: iatom, jatom
      REAL(dp), DIMENSION(3) :: s1, s2, s12, cell_shift, dx

      nedges = 0
      DO iatom = 1, natoms
         DO jatom = 1, natoms
            IF (iatom == jatom) CYCLE
            s1 = MATMUL(hinv, pos(:, iatom))
            s2 = MATMUL(hinv, pos(:, jatom))
            s12 = s1 - s2
            cell_shift = ANINT(s12)
            dx = MATMUL(cell, s12 - cell_shift)
            IF (DOT_PRODUCT(dx, dx) <= cutoff**2) THEN
               nedges = nedges + 1
               IF (PRESENT(edge_index)) THEN
                  edge_index(nedges, :) = [iatom - 1, jatom - 1]
               END IF
               IF (PRESENT(edge_cell_shift)) THEN
                  edge_cell_shift(:, nedges) = REAL(cell_shift, kind=sp)
               END IF
            END IF
         END DO
      END DO
   END SUBROUTINE neighbor_search

END PROGRAM nequip_unittest
