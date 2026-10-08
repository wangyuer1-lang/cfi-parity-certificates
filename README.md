# CFI parity certificates

Exact-arithmetic code and input data accompanying **Exact CFI Parity Formulas for Spin Models**, by Yue Wang.

## Download the manuscript version

[Release v1.0.0](https://github.com/wangyuer1-lang/cfi-parity-certificates/releases/tag/v1.0.0) contains the complete computational supplement as `CFI_Parity_Spin_Models_checks.zip`. It is the archive referenced in Appendix B of the manuscript. Download this named asset to reproduce the fixed version.

The archive SHA-256 is:

```text
2ec09b1d8a5316fe2dc46f54fb90f93b52cf9f62befc70c08c8a97fedc14fed2
```

`SHA256SUMS.txt` inside the archive records the SHA-256 of every other archive member. The repository also exposes those files for inspection. Documentation, citation metadata, and the GitHub workflow are additional repository files.

## Reproduce

Extract the release archive, open a terminal in the extracted directory, and run:

```sh
python -m pip install -r requirements.txt
python run_all.py
```

The full suite uses Python 3.11 or later, SymPy 1.14.0, and NetworkX 3.6.1. The recorded local run used CPython 3.13.14. To check the complete finite graph certificate using only the Python standard library:

```sh
python run_all.py --core
```

Do not use `python -O`: the programs use assertions. A failed check returns a nonzero exit status. Fresh results, per-script logs, and a machine-readable summary are written to `run_results/` or `run_results_core/`.

## What is included

- Explicit edge lists for both 16-vertex, 24-edge witnesses.
- Two independently implemented exact algorithms for the complete component-count/edge-size arrays: exhaustive rollback enumeration and frontier connectivity dynamic programming.
- A separate deletion-contraction calculation of the flow polynomials, compared with both algorithms above.
- Symbolic checks of the trace identity and cubic local tensor, general-degree local checks, finite-group flow checks, and evaluations on actual CFI graphs.
- Expected JSON results, a hash manifest, pinned dependencies, and logs of the recorded local verification run.

The runner performs the calculations in a new temporary directory. Expected results are comparison targets and are not copied into that calculation directory. [README.txt](README.txt) describes each program and its scope.

The GitHub Actions workflow runs the same full suite on a clean Ubuntu runner. Its results are available under the repository's **Actions** tab. This is automated replication of the supplied programs, not a claim of independent external review of the mathematics.

## Notation

In the manuscript, `F_H(q)` is the ordinary nowhere-zero flow polynomial of a graph `H`. The fraktur symbol `\mathfrak F_B(q)` is the cumulant-weighted sum of flow polynomials of vertex-split graphs. They coincide for bases whose degrees are all two or three. The files `certificate_results.json` and `frontier_results.json` report ordinary flow polynomials; `independent_audit_results.json` also reports the vertex-splitting evaluations.

## Citation

Yue Wang. *CFI parity certificates*, version 1.0.0, 2026. Computational supplement to *Exact CFI Parity Formulas for Spin Models*. See [CITATION.cff](CITATION.cff).
