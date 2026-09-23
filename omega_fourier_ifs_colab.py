# OMEGA Fourier–IFS SAT Engine — one-cell Google Colab implementation
# Upload one or more DIMACS .cnf files. Produces and downloads one .txt per CNF.
# Pure Python; no pip installs required.

import sys
import time
from collections import defaultdict

sys.setrecursionlimit(1_000_000)

# Set to None for no time limit.
TIME_LIMIT_SECONDS = 120

def parse_dimacs_text(text):
    nvars_declared = 0
    clauses = []
    pending = []

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("c"):
            continue
        if line.startswith("p"):
            parts = line.split()
            if len(parts) >= 4 and parts[1].lower() == "cnf":
                nvars_declared = int(parts[2])
            continue

        for tok in line.split():
            lit = int(tok)
            if lit == 0:
                clauses.append(tuple(pending))
                pending = []
            else:
                pending.append(lit)

    if pending:
        clauses.append(tuple(pending))

    max_var = max((abs(l) for c in clauses for l in c), default=0)
    nvars = max(nvars_declared, max_var)
    return nvars, canonicalize(clauses)


def canonicalize(clauses):
    out = set()
    for clause in clauses:
        s = set(clause)
        if any(-l in s for l in s):
            # Tautological clause is always satisfied.
            continue
        clean = tuple(sorted(s, key=lambda z: (abs(z), z < 0)))
        out.add(clean)
    return tuple(sorted(out, key=lambda c: (len(c), c)))


def apply_literal(state, lit):
    new_clauses = []
    neg = -lit

    for clause in state:
        if lit in clause:
            continue
        if neg in clause:
            reduced = tuple(x for x in clause if x != neg)
            if len(reduced) == 0:
                return None
            new_clauses.append(reduced)
        else:
            new_clauses.append(clause)

    return canonicalize(new_clauses)


def unit_closure(state):
    assignment = {}
    current = state
    reductions = 0

    while True:
        if current is None:
            return None, None, reductions
        if len(current) == 0:
            return current, assignment, reductions
        if any(len(c) == 0 for c in current):
            return None, None, reductions

        unit = next((c[0] for c in current if len(c) == 1), None)
        if unit is None:
            return current, assignment, reductions

        var = abs(unit)
        val = unit > 0
        if var in assignment and assignment[var] != val:
            return None, None, reductions
        assignment[var] = val

        current = apply_literal(current, unit)
        reductions += 1


def linear_walsh_bias(state):
    """
    Exact sum of degree-1 Walsh coefficients of the clause-satisfaction
    factors g_c = 1 - product_j (1 - sigma_j x_j)/2.

    For a clause of width k, variable x_j contributes sigma_j / 2^k.
    This is only a branch-ordering signal; correctness does not depend on it.
    """
    bias = defaultdict(float)
    occurrence = defaultdict(int)

    for clause in state:
        k = len(clause)
        if k == 0:
            continue
        weight = 2.0 ** (-k)
        for lit in clause:
            var = abs(lit)
            sigma = 1.0 if lit > 0 else -1.0
            bias[var] += sigma * weight
            occurrence[var] += 1

    return bias, occurrence


def choose_branch(state):
    bias, occurrence = linear_walsh_bias(state)
    variables = set(abs(l) for c in state for l in c)
    if not variables:
        return None, (True, False), 0.0

    # Prefer high structural participation, then large Fourier linear bias.
    var = max(
        variables,
        key=lambda v: (occurrence[v], abs(bias[v]), -v)
    )
    b = bias[var]
    preferred = True if b >= 0 else False
    return var, (preferred, not preferred), b


def verify_model(clauses, model):
    for clause in clauses:
        ok = False
        for lit in clause:
            val = model.get(abs(lit), False)
            if (lit > 0 and val) or (lit < 0 and not val):
                ok = True
                break
        if not ok:
            return False
    return True


class Timeout(Exception):
    pass


class OmegaSolver:
    """
    Executable finite prototype of the paper's recursive self-oracle +
    holographic quotient idea.

    Exact quotient key:
        canonical residual CNF state.

    Two paths reaching the same residual formula are merged in memo.
    The linear Walsh signal orders branches; exact SAT correctness comes
    from recursive self-oracular evaluation plus final verification.
    """

    def __init__(self, nvars, clauses, time_limit=TIME_LIMIT_SECONDS):
        self.nvars = nvars
        self.original = canonicalize(clauses)
        self.memo = {}
        self.seen_states = set()
        self.start_time = None
        self.time_limit = time_limit
        self.stats = {
            "recursive_calls": 0,
            "quotient_states": 0,
            "memo_hits": 0,
            "branches": 0,
            "unit_reductions": 0,
            "max_depth": 0,
        }

    def _check_time(self):
        if self.time_limit is not None:
            if time.perf_counter() - self.start_time > self.time_limit:
                raise Timeout()

    def solve(self):
        self.start_time = time.perf_counter()
        try:
            model = self._solve_state(self.original, 0)
            elapsed = time.perf_counter() - self.start_time
            if model is None:
                return "UNSAT", None, elapsed

            # Fill irrelevant/unassigned variables deterministically.
            full = {v: model.get(v, False) for v in range(1, self.nvars + 1)}
            if not verify_model(self.original, full):
                raise RuntimeError("Internal error: produced model failed verification.")
            return "SAT", full, elapsed

        except Timeout:
            elapsed = time.perf_counter() - self.start_time
            return "UNKNOWN_TIMEOUT", None, elapsed

    def _solve_state(self, state, depth):
        self._check_time()
        self.stats["recursive_calls"] += 1
        self.stats["max_depth"] = max(self.stats["max_depth"], depth)

        state = canonicalize(state)
        self.seen_states.add(state)
        self.stats["quotient_states"] = len(self.seen_states)

        if state in self.memo:
            self.stats["memo_hits"] += 1
            cached = self.memo[state]
            return None if cached is None else dict(cached)

        reduced, forced, nred = unit_closure(state)
        self.stats["unit_reductions"] += nred

        if reduced is None:
            self.memo[state] = None
            return None

        if len(reduced) == 0:
            result = dict(forced)
            self.memo[state] = tuple(sorted(result.items()))
            return result

        var, order, _bias = choose_branch(reduced)
        if var is None:
            result = dict(forced)
            self.memo[state] = tuple(sorted(result.items()))
            return result

        self.stats["branches"] += 1

        for val in order:
            lit = var if val else -var
            child = apply_literal(reduced, lit)
            if child is None:
                continue

            sub = self._solve_state(child, depth + 1)
            if sub is not None:
                result = dict(forced)
                result[var] = val
                result.update(sub)
                self.memo[state] = tuple(sorted(result.items()))
                return result

        self.memo[state] = None
        return None


def dimacs_model_line(model, nvars):
    lits = []
    for v in range(1, nvars + 1):
        lits.append(str(v if model.get(v, False) else -v))
    return "v " + " ".join(lits) + " 0"


def solve_dimacs_text(text, filename="input.cnf", time_limit=TIME_LIMIT_SECONDS):
    nvars, clauses = parse_dimacs_text(text)
    solver = OmegaSolver(nvars, clauses, time_limit=time_limit)
    status, model, elapsed = solver.solve()

    lines = [
        "OMEGA Fourier-IFS SAT Engine",
        "============================",
        f"input_file: {filename}",
        f"status: {status}",
        f"variables: {nvars}",
        f"clauses_after_normalization: {len(clauses)}",
        f"elapsed_seconds: {elapsed:.6f}",
        "",
        "Holographic quotient statistics",
        "-------------------------------",
    ]

    for k, v in solver.stats.items():
        lines.append(f"{k}: {v}")

    lines += ["", "Result", "------"]

    if status == "SAT":
        verified = verify_model(clauses, model)
        lines.append(f"verified: {str(verified).lower()}")
        lines.append(dimacs_model_line(model, nvars))
        lines.append("")
        lines.append("assignment:")
        for v in range(1, nvars + 1):
            lines.append(f"x{v}={1 if model[v] else 0}")
    elif status == "UNSAT":
        lines.append("verified_model: n/a")
        lines.append("No satisfying assignment exists according to the exact search.")
    else:
        lines.append("verified_model: n/a")
        lines.append("Search stopped at the configured time limit.")
        lines.append("Increase TIME_LIMIT_SECONDS or set it to None and run again.")

    lines += [
        "",
        "Implementation note",
        "-------------------",
        "The executable quotient merges identical canonical residual CNF states.",
        "The Walsh degree-1 signal is used for branch ordering.",
        "The program is exact when it returns SAT or UNSAT; the final SAT model is",
        "independently checked against the original normalized CNF.",
        "This executable prototype does not by itself establish a polynomial",
        "worst-case bound on the number of quotient states."
    ]

    return "\n".join(lines) + "\n"


def colab_main():
    from google.colab import files

    print("Upload DIMACS .cnf file(s)...")
    uploaded = files.upload()

    cnf_names = [name for name in uploaded if name.lower().endswith(".cnf")]
    if not cnf_names:
        raise ValueError("No .cnf file uploaded.")

    for name in cnf_names:
        raw = uploaded[name]
        text = raw.decode("utf-8", errors="replace")
        report = solve_dimacs_text(text, filename=name)
        out_name = name.rsplit(".", 1)[0] + "_omega_result.txt"

        with open(out_name, "w", encoding="utf-8") as f:
            f.write(report)

        print("\n" + report)
        files.download(out_name)


# In Google Colab this runs automatically after pasting the cell.
try:
    import google.colab  # noqa: F401
    colab_main()
except ImportError:
    pass