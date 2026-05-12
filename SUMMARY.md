# SCAN/R2SCAN DFT-D4 Auto-Reference Implementation Summary

## Date: 2026-05-12

## Objective
Add SCAN and R2SCAN functional support for DFT-D4 auto-reference detection in CP2K.

## Changes Made

### 1. Source Code Changes (`src/qs_dispersion_utils.F`)

**`xc_functional_to_d4_name` function:**
Added SCAN and R2SCAN mappings:
- `"SCAN"` → `"SCAN"`
- `"R2SCAN"` → `"R2SCAN"`

**`xc_functional_detect_expanded` subroutine:**
Added detection for MGGA functionals via subsection names:
- `MGGA_C_R2SCAN` + `MGGA_X_R2SCAN` → R2SCAN
- `MGGA_C_SCAN` + `MGGA_X_SCAN` → SCAN

### 2. Test Files Created

Created 4 new test files in `tests/QS/regtest-dft-d4-auto-ref/`:

| File | Description |
|------|-------------|
| `scan_dftd4_auto.inp` | SCAN auto-detection test |
| `scan_dftd4_explicit.inp` | SCAN with explicit REFERENCE_FUNCTIONAL |
| `r2scan_dftd4_auto.inp` | R2SCAN auto-detection test |
| `r2scan_dftd4_explicit.inp` | R2SCAN with explicit REFERENCE_FUNCTIONAL |

All tests use `&XC_FUNCTIONAL` with MGGA subsections:
```
&XC_FUNCTIONAL
  &MGGA_C_SCAN
  &END MGGA_C_SCAN
  &MGGA_X_SCAN
  &END MGGA_X_SCAN
&END XC_FUNCTIONAL
```

### 3. TEST_FILES.toml Updated

Added energy reference values:
```toml
"scan_dftd4_auto.inp"                   = [{matcher="M011", tol=1.0E-9, ref=-15.79556081484783}]
"scan_dftd4_explicit.inp"               = [{matcher="M011", tol=1.0E-9, ref=-15.79556081484783}]
"r2scan_dftd4_auto.inp"                 = [{matcher="M011", tol=1.0E-9, ref=-15.79775098478576}]
"r2scan_dftd4_explicit.inp"             = [{matcher="M011", tol=1.0E-9, ref=-15.79775098478576}]
```

## Test Results

All 11 tests pass:
- 6 existing auto-ref tests (BLYP, BP, PBE, B3LYP, PBE0, TPSS)
- 4 new SCAN/R2SCAN tests (auto + explicit for each)
- 1 DFTD4 reference functional test

```
Number of FAILED  tests: 0
Number of WRONG   tests: 0
Number of CORRECT tests: 11
```

## Supported Auto-Reference Functionals

Currently supported functionals for DFT-D4 auto-reference detection:
- GGA: PBE, PBE0, B3LYP, BLYP, BP
- MGGA: TPSS, SCAN, R2SCAN

## Notes

- MGGA functionals require LIBXC integration
- Potential type used: GTH-PBE (MGGA-specific potentials not available in this build)
- Basis set: SZV-GTH-PADE
- SCAN and R2SCAN are mapped to themselves via `xc_functional_to_d4_name`