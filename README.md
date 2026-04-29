# Explicit Cost Analysis of Toom–Cook Multiplication for Incomplete NTT

This repository contains the reference implementation for the paper:

**"Explicit cost analysis of Toom–Cook multiplication for incomplete NTT in lattice-based cryptography"**

---

## Overview

This repository provides:

* A concrete implementation of **Toom-4 polynomial multiplication**
* A **Toom-4 / Karatsuba hybrid algorithm**
* Experimental code for **incomplete NTT-based multiplication**
* Verification scripts for **addition-chain optimality**

The goal is to give a **precise operation-count model including constant factors**, suitable for hybrid multiplication strategies in lattice-based cryptography.

---

## Files

* `toom4.py`
  Implementation of:

  * Karatsuba multiplication
  * Toom-4 multiplication
  * Toom-4 / Karatsuba hybrid

* `incomplete_ntt_experiments.py`
  Experimental code for:

  * incomplete NTT
  * performance comparison (Karatsuba / Toom-4 / Hybrid)

* `addition_chain_search.py`
  Exhaustive search code for verifying optimal addition chains used in interpolation

---

## Usage

### 1. Run polynomial multiplication experiments

```
python incomplete_ntt_experiments.py
```

### Expected output

The script prints running times and correctness checks such as:

```
time for    Normal = ...
time for Karatsuba = ...
time for    Toom-4 = ...
time for HybToom-4 = ...

correctness True
correctness True
correctness True
```

---

### 2. Run addition chain verification

```
python addition_chain_search.py
```

This script:

* Enumerates addition chains up to a fixed depth
* Verifies minimality of the chains used in the paper
* Prints chains satisfying required target sets

---

## Notes

### Implementation details

* All computations are over finite fields ( \mathbb{F}_q )
* Polynomial multiplication is performed in coefficient representation
* The implementation focuses on:

  * operation counts
  * structure matching the theoretical model

It is **not optimized for speed**, but for clarity and reproducibility.

---

### Toom-4 interpolation

* Uses evaluation points:

  {0, ±1, ±2, 3, ∞}

* Scalar multiplications are implemented using **addition chains**

* All divisions are postponed and handled by a final scaling step

---

### Addition chains

* Optimality is verified by:

  * exhaustive enumeration (for small depths)
  * simple structural arguments for larger targets

* Verification scripts are included in this repository

---

### Hybrid recursion

* Parameter `L` controls the number of Toom-4 recursion levels
* Base case switches to Karatsuba multiplication
* Final scaling by (120^{-L}) must be precomputed

---

## Requirements

* Python 3.x
* NumPy (for testing)

Install NumPy if needed:

```
pip install numpy
```

---

## Citation

If you use this code, please cite:

```
@article{OkuKudo2026,
  title   = {Explicit cost analysis of Toom--Cook multiplication for incomplete NTT in lattice-based cryptography},
  author  = {Oku, Sakura and Kudo, Momonari},
  year    = {2026}
}
```

---

## License

This project is released under the MIT License.
