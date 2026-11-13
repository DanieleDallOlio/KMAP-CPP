# KMAP-CPP source files

The `src/` directory contains the numerical C++ implementation of the **KMAP-CPP** tracer kinetic-modeling library. It derives from the upstream [KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C) and is built as part of the reusable `kmap` shared library by the root [CMakeLists.txt](../CMakeLists.txt).

Unlike the MATLAB MEX gateways, the core numerical library can be built and used **without MATLAB, Python, or 3D Slicer**. Public headers are located under `../include/`, and additional C++ headers under `../hpp/`.

## Main implementation files

| File | Responsibility |
| --- | --- |
| `kmaplib_common.cpp` | Common numerical functions, frame averaging, exponential convolution, vector norms, and parameter-selection helpers. |
| `kmaplib_infun.cpp` | Input-function utilities, including temporal delay of curves and delay-related sensitivity calculations. |
| `kmaplib_models.cpp` | TAC and Jacobian calculations for supported one-tissue, two-tissue, and liver kinetic models. |
| `kmaplib_optimization.cpp` | Model fitting and optimization, including Levenberg-Marquardt algorithms and bounded coordinate-descent steps. |
| `kmodels.cpp` | Higher-level kinetic-model definitions and implementation support. |
| `cmd.cpp` | Command-line interface support used by the standalone implementation. |
| `parse_args.cpp` | Argument-parsing implementation. |
| `utils.cpp` | Shared implementation utilities. |

For complete function declarations and parameter types, see the public headers, especially [`../include/kmaplib.h`](../include/kmaplib.h).

## Core numerical routines

### Convolution and frame averaging

- `kconv_exp(...)`: exponential convolution used in compartment model evaluations.
- `frame(...)`: average modeled activity across specified PET frame intervals.
- `vecnorm2(...)` and `vecnormw(...)`: unweighted/weighted vector norm helpers.
- `setkin(...)` and `getkin(...)`: select and update fitted parameters.

These are implemented in `kmaplib_common.cpp`.

### Input-function time delay

- `time_delay_tac(...)`: shift a sampled curve by a specified delay.
- `time_delay_jac(...)`: compute delay-related sensitivity information.

These are implemented in `kmaplib_infun.cpp`.

### Model evaluation and Jacobians

- `tac_eval(...)` and `jac_eval(...)`: dispatch model evaluation and sensitivity calculations.
- `kconv_1tcm_tac(...)` / `kconv_1tcm_jac(...)`: 1TCM TAC and Jacobian.
- `kconv_2tcm_tac(...)` / `kconv_2tcm_jac(...)`: 2TCM TAC and Jacobian.
- `kconv_liver_tac(...)` / `kconv_liver_jac(...)`: liver-model TAC and Jacobian.

These routines are implemented in `kmaplib_models.cpp` and use the common convolution and input-function utilities.

### Optimization

- `kmap_levmar(...)`: Levenberg-Marquardt-based parameter fitting.
- `boundpls_cd(...)`: bounded coordinate-descent routine used for constrained updates.

These are implemented in `kmaplib_optimization.cpp`. Optimization convergence and parameter identifiability depend on the model, frame timing, weighting, starting values, and parameter bounds; fitting success should be verified independently.

## Building

For the standalone C++ library, use the root [README.md](../README.md) and `CMakeLists.txt`:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOMP=OFF
cmake --build build --parallel
```

Run these commands from the **repository root**, not from `src/`.

For the MATLAB gateways, use the corresponding instructions in [../mex/README.md](../mex/README.md). The MEX sources compile selected files from this directory and require the public header include path **`../include/`**. The core KMAP-CPP library itself does not require MATLAB.

## Upstream attribution and license

KMAP-CPP preserves the original KMAP-C licensing and contributor attribution. See the repository root [CONTRIBUTORS.md](../CONTRIBUTORS.md) and [LICENSE](../LICENSE). The code is research software; independent numerical and scientific validation remains necessary for each intended application.
