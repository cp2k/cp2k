# DFT-D3 BJ LIB Regression Tests

## Purpose
Tests for DFT-D3 dispersion correction with Becke-Johnson damping using the s-dftd3 library (D3_REFERENCE_CODE mode).

## Description
Tests DFT-D3 BJ with s-dftd3 library instead of dftd3.dat file.
Compare results with regtest-dft-d3-bj-auto-ref tests.

## Requirements
- s_dftd3 library must be linked

## Functionals Tested
- BLYP, B3LYP, BP86, PBE, PBE0, TPSS

## Reference
Compare with regtest-dft-d3-bj-auto-ref for validation of LIB vs AUTO mode consistency.