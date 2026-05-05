! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/io/numpy/crc32.f90
!> Provide crc32 hashing function

!> Implementation of cyclic redundancy check hashing function,
!> used to check the integrity of data in zip files.
module tblite_io_numpy_crc32
   use mctc_env, only : i4, dp
   implicit none
   private

   public :: crc32_hash

   !> Compute crc32 checksum
   interface crc32_hash
      module procedure crc32_hash_char_r0
      module procedure crc32_hash_char_r1
      module procedure crc32_hash_i4_r1
      module procedure crc32_hash_rdp_r1
   end interface crc32_hash

   integer(i4), parameter :: crc_table(0:255) = [ &
      & int(z'00000000', i4), int(z'77073096', i4), int(z'ee0e612c', i4), int(z'990951ba', i4), &
      & int(z'076dc419', i4), int(z'706af48f', i4), int(z'e963a535', i4), int(z'9e6495a3', i4), &
      & int(z'0edb8832', i4), int(z'79dcb8a4', i4), int(z'e0d5e91e', i4), int(z'97d2d988', i4), &
      & int(z'09b64c2b', i4), int(z'7eb17cbd', i4), int(z'e7b82d07', i4), int(z'90bf1d91', i4), &
      & int(z'1db71064', i4), int(z'6ab020f2', i4), int(z'f3b97148', i4), int(z'84be41de', i4), &
      & int(z'1adad47d', i4), int(z'6ddde4eb', i4), int(z'f4d4b551', i4), int(z'83d385c7', i4), &
      & int(z'136c9856', i4), int(z'646ba8c0', i4), int(z'fd62f97a', i4), int(z'8a65c9ec', i4), &
      & int(z'14015c4f', i4), int(z'63066cd9', i4), int(z'fa0f3d63', i4), int(z'8d080df5', i4), &
      & int(z'3b6e20c8', i4), int(z'4c69105e', i4), int(z'd56041e4', i4), int(z'a2677172', i4), &
      & int(z'3c03e4d1', i4), int(z'4b04d447', i4), int(z'd20d85fd', i4), int(z'a50ab56b', i4), &
      & int(z'35b5a8fa', i4), int(z'42b2986c', i4), int(z'dbbbc9d6', i4), int(z'acbcf940', i4), &
      & int(z'32d86ce3', i4), int(z'45df5c75', i4), int(z'dcd60dcf', i4), int(z'abd13d59', i4), &
      & int(z'26d930ac', i4), int(z'51de003a', i4), int(z'c8d75180', i4), int(z'bfd06116', i4), &
      & int(z'21b4f4b5', i4), int(z'56b3c423', i4), int(z'cfba9599', i4), int(z'b8bda50f', i4), &
      & int(z'2802b89e', i4), int(z'5f058808', i4), int(z'c60cd9b2', i4), int(z'b10be924', i4), &
      & int(z'2f6f7c87', i4), int(z'58684c11', i4), int(z'c1611dab', i4), int(z'b6662d3d', i4), &
      & int(z'76dc4190', i4), int(z'01db7106', i4), int(z'98d220bc', i4), int(z'efd5102a', i4), &
      & int(z'71b18589', i4), int(z'06b6b51f', i4), int(z'9fbfe4a5', i4), int(z'e8b8d433', i4), &
      & int(z'7807c9a2', i4), int(z'0f00f934', i4), int(z'9609a88e', i4), int(z'e10e9818', i4), &
      & int(z'7f6a0dbb', i4), int(z'086d3d2d', i4), int(z'91646c97', i4), int(z'e6635c01', i4), &
      & int(z'6b6b51f4', i4), int(z'1c6c6162', i4), int(z'856530d8', i4), int(z'f262004e', i4), &
      & int(z'6c0695ed', i4), int(z'1b01a57b', i4), int(z'8208f4c1', i4), int(z'f50fc457', i4), &
      & int(z'65b0d9c6', i4), int(z'12b7e950', i4), int(z'8bbeb8ea', i4), int(z'fcb9887c', i4), &
      & int(z'62dd1ddf', i4), int(z'15da2d49', i4), int(z'8cd37cf3', i4), int(z'fbd44c65', i4), &
      & int(z'4db26158', i4), int(z'3ab551ce', i4), int(z'a3bc0074', i4), int(z'd4bb30e2', i4), &
      & int(z'4adfa541', i4), int(z'3dd895d7', i4), int(z'a4d1c46d', i4), int(z'd3d6f4fb', i4), &
      & int(z'4369e96a', i4), int(z'346ed9fc', i4), int(z'ad678846', i4), int(z'da60b8d0', i4), &
      & int(z'44042d73', i4), int(z'33031de5', i4), int(z'aa0a4c5f', i4), int(z'dd0d7cc9', i4), &
      & int(z'5005713c', i4), int(z'270241aa', i4), int(z'be0b1010', i4), int(z'c90c2086', i4), &
      & int(z'5768b525', i4), int(z'206f85b3', i4), int(z'b966d409', i4), int(z'ce61e49f', i4), &
      & int(z'5edef90e', i4), int(z'29d9c998', i4), int(z'b0d09822', i4), int(z'c7d7a8b4', i4), &
      & int(z'59b33d17', i4), int(z'2eb40d81', i4), int(z'b7bd5c3b', i4), int(z'c0ba6cad', i4), &
      & int(z'edb88320', i4), int(z'9abfb3b6', i4), int(z'03b6e20c', i4), int(z'74b1d29a', i4), &
      & int(z'ead54739', i4), int(z'9dd277af', i4), int(z'04db2615', i4), int(z'73dc1683', i4), &
      & int(z'e3630b12', i4), int(z'94643b84', i4), int(z'0d6d6a3e', i4), int(z'7a6a5aa8', i4), &
      & int(z'e40ecf0b', i4), int(z'9309ff9d', i4), int(z'0a00ae27', i4), int(z'7d079eb1', i4), &
      & int(z'f00f9344', i4), int(z'8708a3d2', i4), int(z'1e01f268', i4), int(z'6906c2fe', i4), &
      & int(z'f762575d', i4), int(z'806567cb', i4), int(z'196c3671', i4), int(z'6e6b06e7', i4), &
      & int(z'fed41b76', i4), int(z'89d32be0', i4), int(z'10da7a5a', i4), int(z'67dd4acc', i4), &
      & int(z'f9b9df6f', i4), int(z'8ebeeff9', i4), int(z'17b7be43', i4), int(z'60b08ed5', i4), &
      & int(z'd6d6a3e8', i4), int(z'a1d1937e', i4), int(z'38d8c2c4', i4), int(z'4fdff252', i4), &
      & int(z'd1bb67f1', i4), int(z'a6bc5767', i4), int(z'3fb506dd', i4), int(z'48b2364b', i4), &
      & int(z'd80d2bda', i4), int(z'af0a1b4c', i4), int(z'36034af6', i4), int(z'41047a60', i4), &
      & int(z'df60efc3', i4), int(z'a867df55', i4), int(z'316e8eef', i4), int(z'4669be79', i4), &
      & int(z'cb61b38c', i4), int(z'bc66831a', i4), int(z'256fd2a0', i4), int(z'5268e236', i4), &
      & int(z'cc0c7795', i4), int(z'bb0b4703', i4), int(z'220216b9', i4), int(z'5505262f', i4), &
      & int(z'c5ba3bbe', i4), int(z'b2bd0b28', i4), int(z'2bb45a92', i4), int(z'5cb36a04', i4), &
      & int(z'c2d7ffa7', i4), int(z'b5d0cf31', i4), int(z'2cd99e8b', i4), int(z'5bdeae1d', i4), &
      & int(z'9b64c2b0', i4), int(z'ec63f226', i4), int(z'756aa39c', i4), int(z'026d930a', i4), &
      & int(z'9c0906a9', i4), int(z'eb0e363f', i4), int(z'72076785', i4), int(z'05005713', i4), &
      & int(z'95bf4a82', i4), int(z'e2b87a14', i4), int(z'7bb12bae', i4), int(z'0cb61b38', i4), &
      & int(z'92d28e9b', i4), int(z'e5d5be0d', i4), int(z'7cdcefb7', i4), int(z'0bdbdf21', i4), &
      & int(z'86d3d2d4', i4), int(z'f1d4e242', i4), int(z'68ddb3f8', i4), int(z'1fda836e', i4), &
      & int(z'81be16cd', i4), int(z'f6b9265b', i4), int(z'6fb077e1', i4), int(z'18b74777', i4), &
      & int(z'88085ae6', i4), int(z'ff0f6a70', i4), int(z'66063bca', i4), int(z'11010b5c', i4), &
      & int(z'8f659eff', i4), int(z'f862ae69', i4), int(z'616bffd3', i4), int(z'166ccf45', i4), &
      & int(z'a00ae278', i4), int(z'd70dd2ee', i4), int(z'4e048354', i4), int(z'3903b3c2', i4), &
      & int(z'a7672661', i4), int(z'd06016f7', i4), int(z'4969474d', i4), int(z'3e6e77db', i4), &
      & int(z'aed16a4a', i4), int(z'd9d65adc', i4), int(z'40df0b66', i4), int(z'37d83bf0', i4), &
      & int(z'a9bcae53', i4), int(z'debb9ec5', i4), int(z'47b2cf7f', i4), int(z'30b5ffe9', i4), &
      & int(z'bdbdf21c', i4), int(z'cabac28a', i4), int(z'53b39330', i4), int(z'24b4a3a6', i4), &
      & int(z'bad03605', i4), int(z'cdd70693', i4), int(z'54de5729', i4), int(z'23d967bf', i4), &
      & int(z'b3667a2e', i4), int(z'c4614ab8', i4), int(z'5d681b02', i4), int(z'2a6f2b94', i4), &
      & int(z'b40bbe37', i4), int(z'c30c8ea1', i4), int(z'5a05df1b', i4), int(z'2d02ef8d', i4)]

contains

!> Compute crc32 checksum for a character string
pure function crc32_hash_char_r0(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   character(len=*), intent(in) :: val
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, len(val)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(val(ii:ii))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_char_r0

!> Compute crc32 checksum for a character array
pure function crc32_hash_char_r1(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   character(len=1), intent(in) :: val(:)
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, size(val)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(val(ii))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_char_r1

!> Compute crc32 checksum for a 4-byte integer array
pure function crc32_hash_i4_r1(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   integer(i4), intent(in) :: val(:)
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   character(len=1) :: chunk(4)

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, size(val)
      chunk = transfer(val(ii), chunk)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(1))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(2))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(3))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(4))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_i4_r1

!> Compute crc32 checksum for a real array
pure function crc32_hash_rdp_r1(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   real(dp), intent(in) :: val(:)
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   character(len=1) :: chunk(8)

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, size(val)
      chunk = transfer(val(ii), chunk)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(1))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(2))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(3))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(4))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(5))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(6))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(7))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(8))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_rdp_r1

end module tblite_io_numpy_crc32
