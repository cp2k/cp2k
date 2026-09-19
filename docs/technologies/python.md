# Python interface

CP2K's `python/` package provides direct, in-process access to `libcp2k` through the C API. It
offers reusable force environments, NumPy arrays for geometries and forces, complete native CP2K
runs, and an optional ASE calculator. It does not require the obsolete Cython bindings or start a
CP2K shell subprocess.

## Installation

```{include} ../../python/README.md
---
start-after: '## Installation'
---
```
