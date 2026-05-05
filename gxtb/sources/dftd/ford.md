---
project: DFT-D4
symmary: Generally Applicable Atomic-Charge Dependent London Dispersion Correction
author: Grimme group Bonn
src_dir: ./src
include: ./include
output_dir: ./_docs
exclude_dir: ./test
project_github: https://github.com/dftd4/dftd4
github: https://github.com/dftd4
website: https://mctc.uni-bonn.de/software/dftd4
docmark: <
predocmark: >
display: public
         protected
         private
source: true
graph: true
sort: alpha
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            multicharge:https://github.com/grimme-lab/multicharge
            mctc_io:https://grimme-lab.github.io/mctc-lib/module/mctc_io.html
            mctc_env:https://grimme-lab.github.io/mctc-lib/module/mctc_env.html
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

{!README.md!}
