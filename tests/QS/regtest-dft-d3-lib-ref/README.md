# DFT-D3 LIB Regression Tests

## Purpose
Tests for DFT-D3 dispersion correction using the s-dftd3 library (D3_REFERENCE_CODE mode).

## Description
Tests DFT-D3 (zero-damping) with s-dftd3 library instead of dftd3.dat file.
Compare results with regtest-dft-d3-auto-ref tests.

## Requirements
- s_dftd3 library must be linked

## Functionals Tested
- BLYP, B3LYP, BP86, PBE, PBE0, TPSS

## Reference
Compare with regtest-dft-d3-auto-ref for validation of LIB vs AUTO mode consistency.