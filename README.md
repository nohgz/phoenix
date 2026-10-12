<p align="center">
  <img src="https://github.com/nohgz/phoenix/blob/main/assets/phoenix_logo.png?raw=true" alt="Phoenix Logo"/>
</p>

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

Namespaces follow the convention of `where_executes::header_name::function`. For example, `host::config::parse_toml(...)`.

## About
PHOENIX was originally developed for an undergraduate thesis focused on high performance general relativistic radiative transfer.
The poster created for its original presentation can be found here:
> Noga, Christopher, "PHOENIX: GPU-Accelerated Polarized Radiative Transfer for Curved Spacetime" (2026).
> *Symposium on Research and Creative Expression (SORCE)*. 1585.
> https://scholar.valpo.edu/cus/1585

This now-public repo is the second iteration of the code with substantial speedups, more professional commits, and major bugfixes.