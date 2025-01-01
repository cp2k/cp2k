!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Simple wrapper for the official TREXIO Fortran interface
!> \par History
!>      05.2024 created [SB]
!> \author Stefano Battaglia
! **************************************************************************************************

#ifdef __TREXIO
#include <trexio_f.f90>
#endif
