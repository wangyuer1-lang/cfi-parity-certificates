Exact CFI Parity Formulas for Spin Models
Computational material for the major revision, 8 October 2026

QUICK START

Extract the archive into a new directory. In that directory run:
  python -m pip install -r requirements.txt
  python run_all.py

For only the finite graph certificate, with no third-party dependencies:
  python run_all.py --core

The full suite uses SymPy 1.14.0 and NetworkX 3.6.1. It was tested with
CPython 3.13.14 on Windows. NetworkX 3.6.1 requires Python >=3.11.
The core suite uses Python's standard library. Do not use python -O.
No downloads, generated inputs or files from the author's machine are
needed after the dependencies have been installed. The same commands
work in PowerShell and Unix shells. A virtual environment may be used.

WHAT THE RUNNER DOES

run_all.py checks SHA256SUMS.txt, copies the programs and explicit graph
inputs into a new temporary directory, and executes every selected script.
It does not copy expected results into the calculation directory. It checks
exit codes, compares generated JSON against expected/ and compares the
independent algorithms with each other. It writes fresh JSON and logs into
run_results/ (or run_results_core/), including summary.json. Any assertion
failure, mismatched result or failed process gives a nonzero exit code.

The verification_run/ directory records the complete run performed from
a fresh extraction for this revision, including the runtime environment
and individual script logs. Runtime measurements are machine dependent.
SHA256SUMS.txt covers all supplied files other than the manifest itself.

INPUTS AND INDEPENDENT ALGORITHMS

inputs/witness_graphs.json lists both 16-vertex, 24-edge witness graphs.
Each undirected edge occurs once. The data reproduce the manuscript's
adjacency table and its stated two-edge switch.

1. verify_certificates.py (standard library)
   Exhaustively visits all 2^24 edge subsets of each base using a rollback
   disjoint-set structure. It derives full component-count/edge-size arrays,
   full flow coefficient lists, and the square-modification polynomials.
   It also checks simplicity, cubicity, 3-connectivity and induced squares.

2. verify_frontier.py (standard library)
   Imports no code from another checker. Reads the explicit edge lists,
   reconstructs the square modifications, and counts edge subsets by a
   frontier connectivity dynamic program. It computes the full arrays and
   symbolic flow coefficient lists, including loops and parallel edges.
   It separately verifies all 137 vertex-deletion cases per base and girth
   four. The runner compares all its arrays and polynomials with method 1.

3. verify_extensions.py (SymPy)
   Uses deletion-contraction and series reduction for flow polynomials.
   Imports graph definitions and modification routines from method 1,
   but does not use its subset enumeration. The runner compares its
   polynomial coefficients with both other algorithms.
   Also checks vertex-splitting cumulant polynomials against direct flow
   enumeration on nine small bases over Z2, Z3, Z4 and Z2 x Z2 (36 checks).

FORMAL AND SMALL-GRAPH CHECKS

4. check_trace.py (SymPy)
   Checks all 15 four-position equality patterns in the formal variable q.

5. verify_cubic_local.py (SymPy)
   Checks the entire local cubic polynomial in q,v in 67 equality patterns
   realized by three ordered pair-basis labels. Separately forms additive
   exterior matrices for 24 tuples of exact symmetric trace-free matrices
   of dimensions 2 through 7 and checks the new trace reduction.

6. verify_general_degree_local.py (SymPy through its partition helper)
   Checks 2644 complete local Potts tensor coordinates: degrees 2--8 at
   q=2, 2--6 at q=3, and 2--4 at q=4. This checks all lower coefficients
   and the cumulant coefficient itself. It also contains 1089 additional
   rank-one perturbation sanity checks not used as manuscript results.

7. verify_two_spin.py (NetworkX)
   Builds the actual CFI pair over K4 and evaluates its partition functions
   by variable elimination at (beta,gamma,lambda)=(2,3,5). Checks the
   product formula exactly. For C3 it directly enumerates all subsets of
   the actual 18-vertex CFI graphs at (s,t)=(2,3). Both equalities are
   assertions and appear as integer values in two_spin_results.json.

RESULT CONVENTIONS

All flow coefficient arrays are in ascending powers of q. Rank arrays
are indexed as counts[k][s]; zero entries are retained. All comparisons
use exact integers or symbolic polynomials, not floating-point arithmetic.
Finite checks of local tensors accompany the general proofs. The finite
base-graph counts are the certificate used in the all-q threshold theorem.

EXPECTED KEY VALUES

Each base has 16777216 edge subsets and 223080 spanning trees.
Flow values at q=3,4,5,6: 0,30,3336,66900.
Square-weight sums at q=5: 6480 and 6816.
Their polynomial difference is
  -(q-1)(q-2)(q-3)^2(q-4)(q+2).
K4 two-spin difference: 921665794779546558856964111328125000.
C3 induced-edge difference at (2,3): 524288.

The supplement is reproducible material prepared by the author. The
independent algorithms are independent implementations within this work;
no external replication is claimed.
