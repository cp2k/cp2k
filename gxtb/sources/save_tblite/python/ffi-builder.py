# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""
FFI builder module for tblite for usage from meson and from setup.py.

Since meson has the full knowledge about the build, it will handle
the generation of the C definitions in the meson.build file rather
than in the FFI builder. This allows to correctly keep track of
dependencies and updates in the build process.

For setup.py we have to do the preprocessing ourselves here, this
requires us to use the C compiler to preprocess the header file
of tblite because the CFFI C parser cannot handle certain C
preprocessor constructs. Also, we cannot rely on an external build
system fixing dependencies for us and therefore have to find those
ourselves using pkg-config.
"""

import os

import cffi

library = "tblite"
include_header = '#include "' + library + '.h"'
prefix_var = library.upper() + "_PREFIX"
if prefix_var not in os.environ:
    prefix_var = "CONDA_PREFIX"
version = "0.4.0"

if __name__ == "__main__":
    import sys

    kwargs = dict(libraries=[library])

    header_file = sys.argv[1]
    module_name = sys.argv[2]

    with open(header_file) as f:
        cdefs = f.read()
else:
    import subprocess

    try:
        import pkgconfig

        if not pkgconfig.exists(library):
            raise ModuleNotFoundError(f"Unable to find pkg-config package '{library}'")
        if pkgconfig.installed(library, f"< {version}"):
            raise ModuleNotFoundError(
                f"Installed '{library}' version is too old, {version} or newer is required"
            )

        kwargs = pkgconfig.parse(library)
        cflags = pkgconfig.cflags(library).split()

    except ModuleNotFoundError:
        kwargs = dict(libraries=[library])
        cflags = []
        if prefix_var in os.environ:
            prefix = os.environ[prefix_var]
            kwargs.update(
                include_dirs=[os.path.join(prefix, "include")],
                library_dirs=[os.path.join(prefix, "lib")],
                runtime_library_dirs=[os.path.join(prefix, "lib")],
            )
            cflags.append("-I" + os.path.join(prefix, "include"))

    cc = os.environ["CC"] if "CC" in os.environ else "cc"

    module_name = library + "._lib" + library

    p = subprocess.Popen(
        [cc, *cflags, "-DTBLITE_CFFI", "-E", "-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = p.communicate(include_header.encode())

    cdefs = out.decode()

ffibuilder = cffi.FFI()
ffibuilder.set_source(module_name, include_header, **kwargs)
ffibuilder.cdef(cdefs)

if __name__ == "__main__":
    ffibuilder.distutils_extension(".")
