#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# GauXC PR #222 does not change the package version. Unknown installations keep
# the legacy safeguards; the toolchain's patched package advertises this fix.
option(
  CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX
  "External GauXC contains PR #222 (2c236c4) or an equivalent validated backport"
  OFF)
mark_as_advanced(CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX)

set(CP2K_GAUXC_HAS_ONEDFT_GRADIENT_FIX OFF)
if(CP2K_USE_GAUXC AND GAUXC_HAS_ONEDFT)
  if(GAUXC_HAS_ONEDFT_GRADIENT_FIX OR CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX)
    set(CP2K_GAUXC_HAS_ONEDFT_GRADIENT_FIX ON)
  endif()
endif()
if(CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX)
  if(NOT CP2K_GAUXC_HAS_ONEDFT_GRADIENT_FIX)
    message(
      FATAL_ERROR "The GauXC gradient-fix override requires GauXC with OneDFT.")
  endif()
  message(WARNING "Assuming the external GauXC installation contains PR #222. "
                  "This is a user assertion, not automatic version detection.")
endif()
