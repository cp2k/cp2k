# g-xTB/tblite sources

This directory contains the local source snapshots used for the CP2K/tblite g-xTB integration.
The source trees are copied without nested `.git` directories so they can live inside this CP2K
branch as regular files.

The CP2K branch was created from local PR5131 state:

- Branch: `pr-5131`
- Commit: `3e095fb4f07e2898431c9408c12927f0af494609`

## Source snapshots

- `sources/save_tblite`
  - Origin: https://github.com/thfroitzheim/save_tblite.git
  - Commit: `b4b2b79506e31d0b962f9dc9dcb90c0dc3771908`
  - Notes: g-xTB-enabled tblite source used for the CP2K build. The local snapshot includes the
    CMake `Finddftd.cmake` compatibility fix used during the build. The bundled s-D3 dependency
    is moved to `v1.3.2`, which already contains the ATM `triple` gradient fix from the
    `simple-dftd3-1.3.1` patch baseline; the remaining PR5131 s-D3 gradient/smooth-cutoff fixes
    are applied to the local subproject used for testing.

- `sources/dftd`
  - Origin: https://github.com/grimme-lab/dftd.git
  - Commit: `e56749a67f85371bb63095b85fc8d0c4a6a2e0eb`
  - Notes: includes a local port of the PR5131 DFTD4 gradient/virial patch to the newer
    `src/dftd4/model/type.f90` layout. The same patched file was used in the `save_tblite`
    subproject for the local test build.

- `sources/g-xtb`
  - Origin: https://github.com/grimme-lab/g-xtb.git
  - Commit: `05b2089e424fdf65f3664fbde2f231cd4e0f5c8d`

The installed local build used for testing remains outside this source snapshot at
`../gxtb/install-cmake` relative to the workspace root.
