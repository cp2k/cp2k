# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

from os.path import dirname, exists, join
from typing import Optional

# We prefer tomli and tomlkit here, because they are 1.0.0 compliant, while toml is not yet
try:
    import tomli as toml_impl
except ModuleNotFoundError:
    try:
        import tomlkit as toml_impl
    except ModuleNotFoundError:
        try:
            import toml as toml_impl
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "No TOML parser implementation found, install tomli, tomlkit or toml"
            )

_data_base = None


def load_data_base(name: str) -> dict:
    """Load damping parameter database"""

    with open(name) as fh:
        return toml_impl.loads(fh.read())


def get_data_file_name(base_name: str = "parameters.toml") -> str:
    """The data base file is usually shipped with the dftd4 installation in
    $PREFIX/share/dftd4, which we access from $PREFIX/lib/pythonX.Y/site-packages/dftd4.

    If we don't find the parameter file there we might be standalone and ship
    our own data base file in the same directory as this source file.
    """

    data_file = join(
        dirname(__file__), "..", "..", "..", "..", "share", "dftd4", base_name
    )
    if not exists(data_file):
        # for Windows install layout
        data_file = join(
            dirname(__file__), "..", "..", "..", "Library", "share", "dftd4", base_name
        )
        if not exists(data_file):
            data_file = join(dirname(__file__), base_name)

    return data_file


def _get_params(entry: dict, base: dict, defaults: list, method: str, keep_meta=False) -> dict:
    """Retrive the parameters from the data base, make sure the default
    values are applied correctly in the process. In case we have multiple
    defaults search for the first of the list defined for this method."""

    for default in defaults:
        try:
            params = base[default].copy()
            params.update(**entry[default])
            if not keep_meta:
                for key in ("mbd", "damping", "doi"):
                    if key in params: del params[key]
            return params
        except KeyError:
            continue

    raise KeyError("No entry for " + method + " in parameter data base")


def get_damping_param(
    method: str,
    defaults: Optional[list] = None,
    data_file: Optional[str] = None,
    keep_meta=False,
) -> dict:
    """Obtain damping parameters from a data base file."""
    global _data_base

    if _data_base is None:

        if data_file is None:
            data_file = get_data_file_name()

        _data_base = load_data_base(data_file)

    if "default" not in _data_base or "parameter" not in _data_base:
        raise KeyError("No default correct scheme provided")

    if defaults is None:
        defaults = _data_base["default"]["d4"]

    _base = _data_base["default"]["parameter"]["d4"]
    _entry = _data_base["parameter"][method.lower()]["d4"]

    return _get_params(_entry, _base, defaults, method, keep_meta)


def get_all_damping_params(
    defaults: Optional[list] = None, data_file: Optional[str] = None, keep_meta=False
) -> dict:
    """Provide dictionary with all damping parameters available from parameter file"""
    global _data_base

    if _data_base is None:

        if data_file is None:
            data_file = get_data_file_name()

        _data_base = load_data_base(data_file)

    try:
        if defaults is None:
            defaults = _data_base["default"]["d4"]
        _base = _data_base["default"]["parameter"]["d4"]
        _parameters = _data_base["parameter"]
    except KeyError:
        return {}

    definitions = {}

    for method in _parameters:
        try:
            _entry = _parameters[method]["d4"]
            params = _get_params(_entry, _base, defaults, method, keep_meta)
            definitions[method] = params
        except KeyError:
            continue

    return definitions
