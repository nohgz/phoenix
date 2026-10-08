![PHOENIX Logo](/assets/phoenix_logo.png)

**Parallel High-Order EmissioN Integrator for null geodesiX**

A CUDA-accelerated GRRT raymarcher with CPU-only compatibility.

## Installation
TODO (but its probs just have cmake (cuda optional), have the libs, and run pybuild)

## Use Cases + Limitations
Uses the fast light approximation, meaning that `PHOENIX` only works for single snapshots.

## Code Organization Notes
* `src/core` - Code that can run on both the device and host will live here.
* `src/device` - Code that is device only will live here.
* `src/host` - Code that is host only will live here.

Namespaces follow the convention of `where_executes::header_name::function` e.g. `host::config::parse_toml(...)`

## About
PHOENIX was originally developed for an undergraduate thesis focused on high performance general relativistic radiative transfer.
This now-public repo is the second iteration of the code with more professional commits and major bugfixes.