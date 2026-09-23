# OMEGA Fourier-IFS SAT Engine

Official repository for **Holographic Self-Oracularization for Boolean Satisfiability: A Fourier-IFS Framework for the P versus NP Problem** by Kaoru Aguilera Katayama.

## Current implementation

The repository now contains the **exact Fourier–IFS / harmonic self-oracle implementation** described in the updated manuscript.

The previous prototype used residual-CNF recursion with memoization and only used a degree-1 Walsh signal to order branches. That code has been replaced.

The current solver instead implements:

- CNF to 3-CNF conversion for general DIMACS inputs;
- exact clause-local Walsh kernels;
- exact backward transfer operators;
- a canonical reduced ADD/ROBDD quotient for future spectral states;
- exact harmonic amplitudes (H_F(u));
- one-path self-oracular witness construction;
- **no SAT backtracking** and **no recursive DPLL-style branch exploration**;
- independent verification of any returned SAT assignment.

## Files

- `P_equals_NP_Holographic_Paper.pdf` — compiled updated manuscript.
- `P_equals_NP_Holographic_Paper.tex` — complete LaTeX source.
- `references.bib` — bibliography.
- `omega_fourier_ifs_colab.py` — exact one-cell Google Colab / DIMACS implementation.
- `OMEGA-Fourier-IFS-SAT-Engine_release.zip` — release bundle.

## Google Colab

1. Open a Google Colab notebook.
2. Paste the full contents of `omega_fourier_ifs_colab.py` into one cell.
3. Run it.
4. Upload one or more DIMACS `.cnf` files.
5. The program writes and downloads a `*_omega_result.txt` report.

Useful implementation statistics include:

```text
max_reachable_quotient_nodes
harmonic_queries
kernel_transfers
backtracking_branches
recursive_sat_calls
elapsed_seconds
```

For the exact implementation, the last two are reported as zero.

## Complexity expression

The manuscript tracks the algorithmic cost using the quotient width / transition representation. The central empirical object to measure in the implementation is the growth of the canonical ADD quotient rather than the recursion tree of the old prototype.

## Build

```bash
pdflatex P_equals_NP_Holographic_Paper.tex
pdflatex P_equals_NP_Holographic_Paper.tex
```

## Official repository

https://github.com/PolloXDDD/OMEGA-Fourier-IFS-SAT-Engine
