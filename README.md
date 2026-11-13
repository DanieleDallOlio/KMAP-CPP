# KMAP-CPP

## Overview

**KMAP-CPP** is a standalone C++ tracer kinetic-modeling library derived from the open-source [KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C), developed within the broader [Open Kinetic Modeling Initiative](https://www.openkmi.org/).

The original Kinetic Modeling and Analysis Package (KMAP) provided numerical tools for interpreting dynamic positron emission tomography (PET) data. KMAP-CPP retains and extends its C/C++ modeling foundation in a reusable library that can be built independently or used by [SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET).

**KMAP-CPP is a derivative of KMAP-C, not an unrelated reimplementation.** Original upstream authorship, copyright, and MIT licensing are retained; see [CONTRIBUTORS.md](CONTRIBUTORS.md) and [LICENSE](LICENSE).

## What this repository provides

- A **shared C++ library** for supported PET tracer kinetic-model calculations and fitting.
- Input-function sampling and numerical convolution utilities.
- One-tissue, two-tissue, and liver kinetic-model routines inherited and adapted from the KMAP-C codebase.
- Numerical optimization tools, including Levenberg-Marquardt methods and bounded updates.
- Optional **OpenMP** compilation for supported computational paths.
- Standalone C++ examples and reference data.
- Source and legacy precompiled **MATLAB MEX** interfaces kept in a separate `mex/` directory.
- CMake build-tree and installable package support, including the `KMAP::kmap` target.

The core library does **not** require MATLAB, Python, or 3D Slicer to build. Slicer-specific UI, MRML/DICOM integration, segmentation processing, visualization, and clinical workflow choices are implemented in the parent SlicerDynamicPET project, not in KMAP-CPP.

## Building the standalone C++ library

### Requirements

- A **C++17-compatible compiler**.
- **CMake** and a supported native build system.
- Optional OpenMP development/runtime support for accelerated builds.

By default the standalone project builds the examples and enables installation/export targets. You can explicitly disable them for embedding or library-only builds.

### Build without OpenMP

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOMP=OFF
cmake --build build --parallel
```

For multi-configuration generators such as Visual Studio, select the configuration when building:

```bash
cmake --build build --config Release --parallel
```

### Build with OpenMP

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOMP=ON
cmake --build build --parallel
```

Setting `OMP=ON` makes compatible OpenMP support **required** by the KMAP-CPP CMake configuration. Install/configure an OpenMP-capable compiler and runtime first. For example, Apple Clang systems may require a separately supplied OpenMP runtime and suitable compiler configuration. If OpenMP is unavailable, configure with `OMP=OFF` instead.

### Build the library without the examples

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DKMAP_BUILD_EXAMPLES=OFF -DKMAP_INSTALL=OFF -DOMP=OFF
cmake --build build --parallel
```

`KMAP_BUILD_EXAMPLES` and `KMAP_INSTALL` default to `ON` for a standalone build and to `OFF` when KMAP-CPP is added as a subdirectory of another project. The embedding project can override these settings as needed.

## Installing and using KMAP-CPP from CMake

For a standalone install:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DKMAP_INSTALL=ON \
  -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --parallel
cmake --install build
```

An installed consumer can then use:

```cmake
find_package(KMAP REQUIRED CONFIG)

target_link_libraries(my_application PRIVATE KMAP::kmap)
```

Make the installation prefix discoverable to CMake, for example by setting `CMAKE_PREFIX_PATH` when configuring the consuming project. The project also generates `pkg-config` metadata when installation is enabled.

## Examples and reference data

The `example/` directory provides standalone examples for the **one-tissue** and **two-tissue** compartment models. Associated demonstration/reference data are in `data/`.

```text
example/
|-- demo_1tcm_realdata.cpp
`-- demo_2tcm_realdata.cpp

data/
|-- demo_1tcm_realdata/
`-- demo_2tcm_realdata/
```

These examples exercise the standalone C++ modeling code. They do not establish validation for arbitrary clinical tracer and acquisition combinations.

## Repository layout

```text
KMAP-CPP/
|-- CMakeLists.txt
|-- CMake package / pkg-config template files
|-- CONTRIBUTORS.md
|-- LICENSE
|-- README.md
|-- include/        # Public C++ headers, including kmaplib.h
|-- hpp/            # Additional C++ headers
|-- src/            # Kinetic model, numerical, and utility implementations
|-- example/        # Standalone C++ demonstrations
|-- data/           # Example/reference data
`-- mex/            # MATLAB MEX sources and legacy precompiled binaries
```

For a more detailed overview of the implementation, see [src/README.md](src/README.md). For MATLAB-specific interfaces and platform caveats, see [mex/README.md](mex/README.md).

## MATLAB and MEX support

The repository retains MEX source code and precompiled binaries inherited from and developed within the KMAP ecosystem. These are **separate from the CMake-built core library**.

Binaries are included for supported Windows, Linux, and both Intel/Apple Silicon macOS MATLAB MEX formats where present. Compatibility depends on the installed MATLAB release, machine architecture, and compiler/runtime libraries; source recompilation may be necessary. Refer to [mex/README.md](mex/README.md) for the available interfaces and instructions.

## Integration with SlicerDynamicPET

SlicerDynamicPET includes KMAP-CPP as a **pinned Git submodule** and builds it through `add_subdirectory(KMAP-CPP)`. The parent extension configures the backend's optional OpenMP setting based on whether OpenMP was found in the Slicer build environment. It also disables standalone examples and standalone package installation for its embedded build.

Code using KMAP-CPP should distinguish the numerical kernel from the application-level decisions made by SlicerDynamicPET, such as PET timing validation, input-function coverage, voxel selection, fitted model choices, and DICOM export.

## Scientific use and validation

Kinetic-model results depend on the selected model, biological assumptions, frame timing, input function, units, decay correction, starting parameters, optimization settings, and measurement noise. This library is **research software**; downstream users should independently validate results for each intended workflow.

## Upstream project and attribution

KMAP-CPP is based on the [Open KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C). The original KMAP work was developed at the University of California, Davis and subsequently released as open-source software supporting the Open Kinetic Modeling Initiative.

Original KMAP-C copyright, authorship, and contributions are preserved and acknowledged. Please consult [CONTRIBUTORS.md](CONTRIBUTORS.md) and cite the original KMAP publications and model-methodology literature where appropriate.

## License

KMAP-CPP is distributed under the **MIT License**, consistent with upstream KMAP-C. See [LICENSE](LICENSE).
