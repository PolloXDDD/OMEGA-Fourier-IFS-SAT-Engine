# OMEGA Fourier-IFS SAT Engine

Official repository for **Holographic Self-Oracularization for Boolean Satisfiability: A Fourier-IFS Framework for the P versus NP Problem** by Kaoru Aguilera Katayama.

## Files

- `P_equals_NP_Holographic_Paper.pdf` - compiled paper.
- `P_equals_NP_Holographic_Paper.tex` - complete single-file LaTeX source.
- `references.bib` - bibliography database used alongside the manuscript sources.
- `omega_fourier_ifs_colab.py` - one-cell Google Colab/DIMACS CNF implementation.

## Colab usage

1. Open a Google Colab notebook.
2. Paste the contents of `omega_fourier_ifs_colab.py` into one cell.
3. Run the cell.
4. Upload one or more DIMACS `.cnf` files when prompted.
5. The program writes and downloads `<name>_omega_result.txt` for each input.

The result file contains the SAT/UNSAT status, a DIMACS-compatible assignment when SAT, a direct verification of the returned model, and quotient/search statistics.

## Paper complexity expression

The paper writes the solver cost as

```text
T_SAT(F) = O(n m Q(|F|) R(|F|)),
```

where `Q(|F|)` bounds the number of holographic quotient states and `R(|F|)` bounds the arithmetic/transition cost per quotient state. Under polynomial bounds on `Q` and `R`, the resulting bound is polynomial in the formula encoding size.

## Build

With a standard TeX installation:

```bash
pdflatex P_equals_NP_Holographic_Paper.tex
```

The bibliography is embedded in the final `.tex`; `references.bib` is included as a reusable bibliography source.

## Official repository

https://github.com/PolloXDDD/OMEGA-Fourier-IFS-SAT-Engine
