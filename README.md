# ATLAS-GP

**Adaptive Topology-driven Landscape-Aware Search with Gaussian Process guidance**

A high-performance cluster structure optimizer for all 103 elements, combining basin-hopping with a sparse GP surrogate, eigenvector escape, parallel tempering, and optional MACE-MP-0 machine-learning potential.

![Platform](https://img.shields.io/badge/platform-Windows-blue)
![Language](https://img.shields.io/badge/language-C%2B%2B17-orange)

---

## Features

- **Universal potential engine** — LJ, Morse, Gupta/RGL EAM, UFF, and MACE-MP-0 (ML)
- **All 103 elements** — atomic and molecular clusters
- **GP-guided basin-hopping** — sparse FITC GP surrogate with acquisition function
- **Eigenvector escape** — Lanczos softest-mode saddle-point escape
- **Parallel tempering** — geometric temperature ladder with Metropolis-Hastings swaps
- **Population crossover** — elite pool with diversity-preserving crossover
- **Active learning** — uncertainty-driven MACE oracle calls
- **ACE-inspired fingerprints** — rotation/translation/permutation invariant
- **Full GUI** — clean sidebar + 3-tab layout, all parameters wired, run history, ETA

---

## Quick Start (Windows)

1. Download or clone this repository
2. Double-click **`ATLAS.bat`**
   - Automatically installs MSYS2, GCC, CMake if needed
   - Compiles `atlas.exe`
   - Launches the GUI
3. Enter a formula (e.g. `Au13`), select a potential, click **▶ Run**

> **First run** takes 2–4 minutes to compile. Subsequent launches are instant.

---

## Requirements

- Windows 10/11 (64-bit)
- Internet connection for first-time setup (MSYS2 + dependencies auto-installed)
- Python 3.10+ (auto-installed if missing)
- ~2 GB disk space

### Optional: MACE-MP-0 (ML potential)
Requires Python 3.11. Click **Install MACE** in the GUI or run:
```
py -3.11 -m pip install torch --index-url https://download.pytorch.org/whl/cpu
py -3.11 -m pip install ase mace-torch
```

---

## Usage

### GUI
Launch `ATLAS.bat` and use the graphical interface. Preset profiles available:
- **Fast** — quick scan, 4 runs × 150 generations
- **Standard** — balanced, 6 runs × 300 generations  
- **Thorough** — deep search, 10 runs × 500 generations
- **MACE+AL** — ML potential with active learning

### Command Line
```
atlas.exe Au13 --pot GUPTA --runs 6 --gens 300 --workers 4 --progress
```

Full flag reference:
```
atlas.exe --help
```

---

## Potentials

| Potential | Best for | Notes |
|-----------|----------|-------|
| `LJ`      | Noble gases (Ar, Ne, Kr, Xe) | Lennard-Jones 12-6 |
| `MORSE`   | Simple diatomics, light elements | NIST parameters |
| `GUPTA`   | FCC/HCP metals (Au, Cu, Ni, Pt, Ag, Pd…) | RGL EAM, 16 metals |
| `UFF`     | Mixed/organic, any element | Universal Force Field |
| `MACE`    | Any system, highest accuracy | Requires ML install |
| `AUTO`    | Automatic selection | Cascade: GUPTA→UFF→LJ |

---

## Project Structure

```
ATLAS_2/
├── src/
│   ├── main.cpp                  # Entry point, argument parsing
│   ├── system.hpp                # System struct, formula parser
│   ├── fingerprint.hpp           # ACE-inspired descriptors
│   ├── isomer_filter.hpp         # RMSD-based isomer deduplication
│   ├── gp_surrogate.hpp          # Sparse FITC Gaussian process
│   ├── acquisition.hpp           # GP acquisition + gradient proposals
│   ├── escape.hpp                # Eigenvector/Lanczos escape
│   ├── parallel_tempering.hpp    # PT engine
│   ├── population.hpp            # Elite population crossover
│   ├── adaptive_descriptor.hpp   # ARD descriptor weights
│   ├── natural_gradient.hpp      # Natural gradient GP hyperparams
│   ├── active_learning.hpp       # MACE active learning
│   ├── atlas_worker.hpp          # Main basin-hopping worker
│   ├── lbfgs_solver.hpp          # L-BFGS geometry optimizer
│   ├── io.hpp                    # XYZ file I/O
│   └── potential/
│       ├── lj_potential.hpp
│       ├── morse_potential.hpp
│       ├── gupta_potential.hpp
│       ├── uff_potential.hpp
│       ├── mace_potential.hpp
│       └── potential_factory.hpp
├── gui/
│   └── atlas_gui.py              # Tkinter GUI (Python 3.10+)
├── ml/
│   ├── mace_server.py            # MACE-MP-0 subprocess server
│   └── install_mace.py           # MACE installer (auto-detects Python)
├── CMakeLists.txt
├── ATLAS.bat                     # One-click build + launch
└── ATLAS.vbs                     # Silent double-click launcher
```

---

## Algorithm Overview

```
For each run:
  1. Random initialisation inside bounding sphere
  2. L-BFGS local minimisation
  3. Basin-hopping step:
     a. GP acquisition → propose candidate (gradient-blended)
     b. Stress / graph / motif / homotop perturbation
     c. Eigenvector escape if stuck
     d. L-BFGS minimisation
     e. Metropolis acceptance
     f. Update GP surrogate (FITC, M=150 inducing points)
     g. ARD descriptor weight update (every 80 obs)
     h. Natural gradient GP hyperparameter update (every 60 gens)
  4. Parallel tempering replica exchange (every generation)
  5. Population crossover (every 15 gens)
  6. Active learning MACE call if uncertainty > threshold
Collect → RMSD filter → energy window → rank isomers → write XYZ
```

---

## Output

Results are written to `atlas_results/` (configurable):
- `<formula>_isomers.xyz` — all unique isomers, ranked by energy
- XYZ comment line format: `Isomer=N Energy=E Formula=F Potential=P`

---

## Building Manually

If you prefer to build without `ATLAS.bat`:

```bash
# In MSYS2 MinGW64 shell:
mkdir build && cd build
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j8
```

Requires: GCC 11+, CMake 3.16+, Eigen3 (auto-downloaded), LBFGSpp (auto-cloned)

---

## Citation

If you use ATLAS-GP in research, please cite:
```
ATLAS-GP: Adaptive Topology-driven Landscape-Aware Search with GP guidance
https://github.com/megamsndm/ATLAS-GP
```
