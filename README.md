# DT: Constrained Tetrahedralization Tool

This repository provides the executable tools, partial source code, and experimental data used in the paper "Robust Constrained Tetrahedralization with Steiner-point-free Boundaries".

The project implements a robust constrained tetrahedralization pipeline with boundary recovery and FHC-based Steiner point insertion.

Due to commercial licensing restrictions, only the parts of the source code directly related to the paper are released.

---

# Directory Structure

bin/
    linux/ Linux executables
    windows/
        Debug/
        Release/

examples/
    Thingi10K_Valid/ Valid triangle meshes from the Thingi10K dataset

results/ Experimental executables, scripts, and data

src/ Open-sourced parts of the implementation


---

# Executables

Precompiled executables are provided in the `bin/` directory.

Main executable: dt.exe

Associated library: dt.lib, dt_API.h

To display command-line options:

dt.exe -h

---

# Command Line Usage

## Basic Tetrahedralization

./dt.exe --input <input_mesh> --out 1

Example:

./dt.exe --input ../../../examples/Thingi10K_Valid/79851.obj --out 1

## Tetrahedralization with Refinement and Optimization

./dt.exe --input <input_mesh> --refine 1 --out 1

Example:

./dt.exe --input ../../../examples/Thingi10K_Valid/79851.obj --refine 1 --out 1

Run `dt.exe -h` for the full list of options.

# Example Data

The `examples/` directory contains input meshes used in the experiments.

## Thingi10K Dataset

Contains valid triangle meshes from the **Thingi10K dataset**.  
All files can be downloaded from:
https://drive.google.com/file/d/1HkB_pGzyCHySINoAWIX6d1YrFXFoAL6c/view?usp=sharing

## Meshes Appearing in the Paper

### Ablation Study

Ablation_chazel_10.vtk, Ablation_polyhedron12.vtk

These meshes correspond to the models used in the **Ablation Study** section of the paper.

### Application Example

Application_A320_20_surface_for_tet.vtk

This mesh corresponds to the **A320 aircraft example** shown in the Application section.

# Open-Sourced Components

Due to commercial licensing constraints, only the algorithmic components directly related to the paper are released.

## Volume-Based Mesh Smoothing

src/dt_opt.cpp

Implements the **volume-based mesh smoothing method** used for mesh quality improvement.

## Boundary Recovery and FHC-Based Steiner Point Insertion

src/dt.cpp

Contains implementations of:

- Boundary recovery workflow
- FHC-based Steiner point insertion

These correspond to the main algorithmic contributions described in the paper.

---

# Experimental Results

The `results/` directory includes:

- Executables used for evaluation
- Testing scripts
- Experimental datasets
- Data used to generate tables and figures in the paper

These resources allow researchers to reproduce the experimental results reported in the paper.

---

# Build Environment

The provided binaries were compiled and tested under the following environment.

## Linux

Operating System  
- Ubuntu 22.04.5 LTS
- Architecture: x86_64  

Compiler  
- GCC 11.4.0  
- G++ 11.4.0  

## Windows

Operating System  
- Windows 11  

Compiler  
- Microsoft Visual Studio 2022  
- MSVC (C++14 standard)

Build Configuration  
- Debug  
- Release

# Code Availability

The complete source code of DT cannot be fully released due to **commercial licensing agreements**.

To support reproducibility, this repository provides:

- Precompiled executables
- Core algorithmic modules described in the paper
- Example datasets
- Experimental scripts and result data

Researchers can reproduce the reported results using the provided materials.

---

# Citation

If you use this software or the provided datasets in your research, please cite the corresponding paper.

(To be added after publication)