# This file is part of dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd3.  If not, see <https://www.gnu.org/licenses/>.


from packaging.version import Version
from dftd3 import __version__
from dftd3.library import get_api_version


def test_api_version():
    """Ensure that the API version is compatible."""
    assert Version(get_api_version()) == Version(__version__)
