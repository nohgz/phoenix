<p align="center">
  <img src="https://github.com/nohgz/phoenix/blob/main/assets/phoenix_logo.png?raw=true" alt="PHOENIX Logo"/>
</p>

**Parallel High-Order EmissioN Integrator for null geodesiX**

A CUDA-accelerated GRRT raymarcher with CPU-only compatibility.

## Installation
The requirements below are what I used to compile `PHOENIX`. Older versions of some requirements may still work but **are not guaranteed to**.
| Requirement    | Version |
| -------- | ------- |
| CMake    | 3.20    |
| C++      | 17      |
| OpenMP   | 2021.11 |
| HDF5     | 2.1.1   |
| CUDA (optional)     | 13.4    |

## Compilation
To compile, run the convenience script `pb.py` which should handle the raw CMake compile instructions.
```
python3 pb.py
```

If Python is unavailable, you can run
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
mv build/PHOENIX .
```
to build `PHOENIX` in FP32 mode. Some useful flags are below:

| Flag | Description |
| -------- | ------- |
| `-DPHX_ACCURACY` | Sets the precision of the `scalar` variable. Can be either `-DPHX_ACCURACY=FP32` or `-DPHX_ACCURACY=FP64`. `FP32` is the default operating mode. |
| `-DPHX_FORCE_CPU` | Forces CPU-only compilation if turned on. Can be either `-DPHX_FORCE_CPU=ON` or `-DPHX_FORCE_CPU=OFF`. `OFF` is the default operating mode.|

## Usage
Your first `PHOENIX` run will look something like this:
1. Download a BHAC GRMHD snapshot to the `input` directory. You can find the snapshot below on the `RAPTOR` [GitHub](https://github.com/jordydavelaar/raptor).
```
wget -P input https://astro.ru.nl/~jordyd/data2000.dat
```
2. Run `PHOENIX`. For this test, we'll pass some arguments via the command line. The format matches the `input/config.toml` keys and values.
```
./PHOENIX camera.inclination=10 plasma.m_unit=1.2e18
```
3. If all goes well, you should see
```
[ INFO ] Flux at 2.3000e+11 Hz (Jy): I +4.548971714e+00 Q +3.595375254e-03 U -3.479759672e-02 V +2.150556623e-02
```
and output dumps should be written to `output/`! You can process the HDF5 files using the `RAPTOR` rapplot.py toolchain. A tutorial can be found on the [RAPTOR GitHub](https://github.com/jordydavelaar/raptor#images)

## Limitations
Currently, `PHOENIX`
* only works for single GRMHD snapshots, and extremely large snapshots may overflow GPU VRAM.
* only supports BHAC dumps, but other GRMHD code support is planned.
* has really only been validated for the test dump (`data2000.dat`) provided by the `RAPTOR` group.

## Code Organization Notes
* `src/core` - Code that can run on both the device and host will live here.
* `src/device` - Code that is device-only will live here.
* `src/host` - Code that is host-only will live here.

Namespaces follow the convention of `where_executes::header_name::function`. For example, `host::config::parse_toml(...)`.

## Planned Features
* Doxygen documentation (the docblocks are mostly there; just have to generate it)
* Real-time rendering interface using ImGui
* Python library for running parameter sweeps
* Support for arbitrary metrics & GRMHD codes

... and more!

## About
The single-pass transfer operator formalism follows [Landi Degl'Innocenti & Landi Degl'Innocenti (1985)](https://ui.adsabs.harvard.edu/abs/1985SoPh...97..239L/abstract) as I found in [Dexter (2016, Sec. 2.5.2)](https://arxiv.org/abs/1602.03184).
`PHOENIX` accumulates the contribution vectors (Dexter 2016, Eqs. 69-70) in the same step as it solves the geodesic equation while marching from the camera to the BH.

`PHOENIX` was originally developed for an undergraduate senior thesis focused on high-performance GRRT.
The poster created for its original presentation can be found here:
> Noga, Christopher, "PHOENIX: GPU-Accelerated Polarized Radiative Transfer for Curved Spacetime" (2026).
> *Symposium on Research and Creative Expression (SORCE)*. 1585.
> https://scholar.valpo.edu/cus/1585

This now-public repo is the second iteration of the code with substantial speedups, more professional commits (the old ones bring great shame upon my family's name), and major bug fixes.

## AI Usage Disclosure
Claude assisted in the development of `PHOENIX`. Its main roles were code cleanup, bug fixing, optimization, and helping interpret various papers into functions for physics routines. During development, I have personally validated numerous functions against their `RAPTOR` counterparts; however, some sneaky bugs may still remain.

## Acknowledgements
* (`RAPTOR` Authors) Thomas Bronzwaer, Jordy Davelaar, Oliver Porth, Ziri Younsi, Monika Moscibrodzka, Joost de Kleuver, Bram van der Berg, Renze Oosterhuis, Jesse Vos
* Valparaiso University Department of Physics and Astronomy

*If anyone is left out, please let me know at `chris.noga@valpo.edu`!*
