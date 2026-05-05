# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

os = host_machine.system()
fc = meson.get_compiler('fortran')
cc = fc
if has_cc
  cc = meson.get_compiler('c')
endif
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    '-Wno-maybe-uninitialized',
    language: 'fortran',
  )
  if fc.version().version_compare('<10')
    message('tblite it not tested with GNU Fortran < 10, functionality might be impacted. The build automatically adds ‘-fno-range-check’ for compatibility. Please upgrade to a newer compiler version if possible.') 
    add_project_arguments(
      '-fno-range-check',
      language: 'fortran',
    )
  endif
elif fc_id == 'intel' or fc_id == 'intel-llvm'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    lapack_vendor = 'mkl'
  endif
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    mkl_dep += cc.find_library('mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += cc.find_library('mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += cc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += cc.find_library('mkl_core')
  lib_deps += mkl_dep
  add_project_arguments(['-DWITH_MKL'], language: 'fortran')
elif lapack_vendor == 'mkl-rt'
  mkl_dep = cc.find_library('mkl_rt')
  lib_deps += mkl_dep
  add_project_arguments(['-DWITH_MKL'], language: 'fortran')
elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = cc.find_library('openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = cc.find_library('lapack')
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  foreach lib: get_option('custom_libraries')
    lib_deps += cc.find_library(lib)
  endforeach

else
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = cc.find_library('lapack')
  endif
  lib_deps += lapack_dep
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = cc.find_library('blas')
  endif
  lib_deps += blas_dep
endif

# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  version: '>=0.5.1',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
)
lib_deps += mctc_dep

# Create DFT-D library as subproject
dftd_dep = dependency(
  'dftd',
  # version: ['>=4.0.0'],
  fallback: ['dftd', 'dftd4_dep'],
  default_options: ['default_library=static', 'api=false', 'python=false'],
)
lib_deps += dftd_dep

# Create DFT-D3 library as subproject
sdftd3_dep = dependency(
  's-dftd3',
  version: '>=0.6.0',
  fallback: ['s-dftd3', 'sdftd3_dep'],
  default_options: ['default_library=static', 'api=false', 'python=false'],
)
lib_deps += sdftd3_dep

# Create TOML Fortran as subproject
tomlf_dep = dependency(
  'toml-f',
  version: '>=0.4.3',
  fallback: ['toml-f', 'tomlf_dep'],
  default_options: ['default_library=static'],
)
lib_deps += tomlf_dep
