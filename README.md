# MITHRA-wfollowm

This repository contains the **wfollowm development version** of MITHRA, a full-wave numerical solver for free-electron laser radiation.

The original MITHRA code was developed by Arya Fallahi and collaborators. This modified version is maintained as my private research and development version. Permission has been obtained from the original author to use, modify, and distribute this modified version.

This repository is **not an official MITHRA release**. The copyright of the original MITHRA code remains with the original author and contributors.

---

## About the original MITHRA code

MITHRA is a full-wave numerical solver for free-electron lasers.

Original information:

```text
MITHRA-2.0 (Completely Numerical Calculation of Free Electron Laser Radiation)
Version 2.0, copyright 2019, Arya Fallahi
```

The original code provides a time-domain solver for free-electron laser radiation, including particle-field interaction, FDTD field marching, and radiation analysis.

---

## Main modifications in this version

The current wfollowm development version is mainly used for research-oriented radiation simulations and diagnostics. The major changes include:

1. **Lab-frame detector plane output**

   A fixed laboratory-frame detection plane has been added to record the electromagnetic field at a specified longitudinal position during the simulation.

2. **HDF5 detector field output**

   The detector output can be written to HDF5 files, including time-dependent electromagnetic field data on the detector plane.

3. **Poynting-flux diagnostics**

   Additional diagnostics have been added for calculating detector-plane radiation power and Poynting-flux-related quantities.

4. **Boundary damping treatment**

   Boundary damping / absorbing-layer treatment has been modified to reduce artificial field reflection from the finite computational domain.

5. **Post-processing interfaces**

   The output format has been adapted for post-processing analysis, including radiation spectrum, transverse intensity distribution, and spatial coherence calculations.

This development version is intended for research use and may contain experimental features.

---

## Compilation

Compiling options used for `gcc 7.4.0`:

```bash
make
```

to compile the binary, and

```bash
make install
```

to compile the library and header directory.

---

## Code structure

```text
mithra.cpp
    Main program file.

stdinclude
    Set of different standard functions and constants used in the code.

readdata
    Implementation of the functions reading the lines of the job file.

radiation
    Implementation of the functions for the radiation analysis.

fieldvector
    Implementation of the field vector class for the analysis.

datainput
    Implementation of the parameter parser.

classes
    Implementation of the classes used in the code.

fdtd
    Implementation of the real FDTD time-marching solution class for the MITHRA code without space-charge.

fdtdSC
    Implementation of the real FDTD time-marching solution class for the MITHRA code with space-charge.

database
    Implementation of various data structures used for simulations.

solver
    Mother class for solving the free-electron laser problem.
```

---

## Citation and acknowledgement

If this modified version is used in future work, please acknowledge the original MITHRA code and cite the original MITHRA references as requested by the original author.

The modifications in this repository were developed for radiation simulation, detector-plane field output, and related post-processing analysis.