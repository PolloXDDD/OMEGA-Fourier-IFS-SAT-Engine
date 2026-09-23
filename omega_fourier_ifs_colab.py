# OMEGA Fourier–IFS SAT Engine v2 — exact spectral/IFS implementation
# One-cell Google Colab script: upload DIMACS .cnf, receive *_omega_result.txt.
#
# Core difference from v1:
#   * NO recursive SAT branching / backtracking / DPLL tree.
#   * Computes the exact harmonic self-oracle H_F(u) through the Fourier-Walsh
#     transfer operators from the paper.
#   * Uses a canonical reduced Algebraic Decision Diagram (ADD) as the exact
#     holographic quotient of backward spectral query states.
#   * Direction is chosen once per variable from H_F(u b) > 0.
#
# Pure Python. No pip installs.

import sys
import time
import heapq
from collections import defaultdict

sys.setrecursionlimit(1_000_000)
TIME_LIMIT_SECONDS = 120

# ----------------------------- DIMACS / 3-CNF -----------------------------

def canonicalize_clauses(clauses):
    out = []
    seen = set()
    for clause in clauses:
        s = set(int(x) for x in clause)
        if any(-l in s for l in s):
            continue  # tautology
        c = tuple(sorted(s, key=lambda z: (abs(z), z < 0)))
        if c not in seen:
            seen.add(c)
            out.append(c)
    return tuple(out)


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
    clauses = canonicalize_clauses(clauses)
    mx = max((abs(l) for c in clauses for l in c), default=0)
    return max(nvars_declared, mx), clauses


def to_3cnf(nvars, clauses):
    """Polynomial equisatisfiable conversion. Returns (new_nvars, new_clauses)."""
    out = []
    nxt = nvars + 1
    for c in clauses:
        k = len(c)
        if k <= 3:
            out.append(c)
            continue
        # (l1 v l2 v y1) & (~y1 v l3 v y2) & ... & (~y_r v l_{k-1} v l_k)
        lits = list(c)
        prev_y = None
        for i in range(k - 3):
            y = nxt
            nxt += 1
            if i == 0:
                out.append((lits[0], lits[1], y))
            else:
                out.append((-prev_y, lits[i + 1], y))
            prev_y = y
        out.append((-prev_y, lits[-2], lits[-1]))
    return nxt - 1, canonicalize_clauses(out)


def verify_model(clauses, model):
    for c in clauses:
        if not any((lit > 0 and model.get(abs(lit), False)) or
                   (lit < 0 and not model.get(abs(lit), False)) for lit in c):
            return False
    return True


# -------------------------- Fourier clause kernels -------------------------

def clause_walsh_kernel(clause):
    """
    Integer numerators of the Walsh coefficients of
        g_C(x) = 1 - prod_j (1 - sigma_j x_j)/2.
    Common denominator = 2^k, k=len(clause).

    Returns tuple((xor_mask, integer_numerator), ...), k.
    Variable v corresponds to bit 1<<(v-1) in xor_mask.
    """
    k = len(clause)
    if k == 0:
        return ((0, 0),), 0  # identically false clause factor
    vars_ = [abs(l) for l in clause]
    sigmas = [1 if l > 0 else -1 for l in clause]
    terms = []
    for sm in range(1 << k):
        mask = 0
        prod_sigma = 1
        r = 0
        for i in range(k):
            if (sm >> i) & 1:
                mask ^= 1 << (vars_[i] - 1)
                prod_sigma *= sigmas[i]
                r += 1
        if sm == 0:
            num = (1 << k) - 1
        else:
            # - coeff(v_C) * 2^k = -(-1)^r * prod_sigma
            num = (-1 if (r % 2 == 0) else 1) * prod_sigma
        if num:
            terms.append((mask, num))
    return tuple(terms), k


# ---------------------- Variable / clause ordering -------------------------

def greedy_min_degree_order(nvars, clauses):
    """Deterministic graph elimination heuristic; exactness is order-independent."""
    graph = [set() for _ in range(nvars + 1)]
    active = [False] * (nvars + 1)
    for c in clauses:
        vs = sorted(set(abs(l) for l in c))
        for v in vs:
            active[v] = True
        for i, a in enumerate(vs):
            for b in vs[i + 1:]:
                graph[a].add(b)
                graph[b].add(a)
    # Include isolated declared variables at the end.
    alive = set(v for v in range(1, nvars + 1) if active[v])
    heap = [(len(graph[v]), v) for v in alive]
    heapq.heapify(heap)
    order = []
    while alive:
        while True:
            deg, v = heapq.heappop(heap)
            if v in alive and deg == len(graph[v] & alive):
                break
        nbrs = list(graph[v] & alive)
        # fill clique
        for i, a in enumerate(nbrs):
            for b in nbrs[i + 1:]:
                if b not in graph[a]:
                    graph[a].add(b)
                    graph[b].add(a)
        alive.remove(v)
        order.append(v)
        for u in nbrs:
            if u in alive:
                heapq.heappush(heap, (len(graph[u] & alive), u))
    used = set(order)
    order.extend(v for v in range(1, nvars + 1) if v not in used)
    return order


# ---------------- Canonical reduced Algebraic Decision Diagram -------------

class ADDManager:
    """
    Reduced ordered ADD with arbitrary-precision integer terminals.

    A node represents an exact function f:{0,1}^N -> Z on Fourier masks S.
    Canonical reduction (unique table + low==high elimination) gives the
    executable quotient: equal future spectral functions have the same node id.
    """
    def __init__(self, variable_order, deadline=None):
        self.order = tuple(variable_order)
        self.rank = {v: i for i, v in enumerate(self.order)}
        self.deadline = deadline

        # node id -> (var, low, high); var=0 means terminal, low stores integer value
        self.nodes = []
        self.term_unique = {}
        self.node_unique = {}
        self.add_cache = {}
        self.scale_cache = {}
        self.shift_cache = {}
        self.transfer_cache = {}

        self.ZERO = self.terminal(0)
        self.ONE = self.terminal(1)
        self.peak_reachable = 1
        self.transfer_steps = 0

    def _check_time(self):
        if self.deadline is not None and time.perf_counter() > self.deadline:
            raise TimeoutError

    def terminal(self, value):
        value = int(value)
        nid = self.term_unique.get(value)
        if nid is not None:
            return nid
        nid = len(self.nodes)
        self.nodes.append((0, value, -1))
        self.term_unique[value] = nid
        return nid

    def is_terminal(self, u):
        return self.nodes[u][0] == 0

    def tvalue(self, u):
        return self.nodes[u][1]

    def var(self, u):
        return self.nodes[u][0]

    def mk(self, var, low, high):
        if low == high:
            return low
        key = (var, low, high)
        nid = self.node_unique.get(key)
        if nid is not None:
            return nid
        nid = len(self.nodes)
        self.nodes.append(key)
        self.node_unique[key] = nid
        return nid

    def scale(self, u, c):
        c = int(c)
        if c == 0:
            return self.ZERO
        if c == 1:
            return u
        key = (u, c)
        hit = self.scale_cache.get(key)
        if hit is not None:
            return hit
        self._check_time()
        if self.is_terminal(u):
            out = self.terminal(self.tvalue(u) * c)
        else:
            v, lo, hi = self.nodes[u]
            out = self.mk(v, self.scale(lo, c), self.scale(hi, c))
        self.scale_cache[key] = out
        return out

    def add(self, a, b):
        if a == self.ZERO:
            return b
        if b == self.ZERO:
            return a
        if a > b:  # commutative cache normalization
            a, b = b, a
        key = (a, b)
        hit = self.add_cache.get(key)
        if hit is not None:
            return hit
        self._check_time()
        if self.is_terminal(a) and self.is_terminal(b):
            out = self.terminal(self.tvalue(a) + self.tvalue(b))
        else:
            va = None if self.is_terminal(a) else self.var(a)
            vb = None if self.is_terminal(b) else self.var(b)
            if va is None:
                v = vb
            elif vb is None:
                v = va
            else:
                v = va if self.rank[va] <= self.rank[vb] else vb

            if va == v:
                _, alo, ahi = self.nodes[a]
            else:
                alo = ahi = a
            if vb == v:
                _, blo, bhi = self.nodes[b]
            else:
                blo = bhi = b
            out = self.mk(v, self.add(alo, blo), self.add(ahi, bhi))
        self.add_cache[key] = out
        return out

    def xor_shift(self, u, mask):
        """Return function S -> f(S xor mask)."""
        if mask == 0 or self.is_terminal(u):
            return u
        key = (u, mask)
        hit = self.shift_cache.get(key)
        if hit is not None:
            return hit
        self._check_time()
        v, lo, hi = self.nodes[u]
        slo = self.xor_shift(lo, mask)
        shi = self.xor_shift(hi, mask)
        if mask & (1 << (v - 1)):
            out = self.mk(v, shi, slo)
        else:
            out = self.mk(v, slo, shi)
        self.shift_cache[key] = out
        return out

    def apply_kernel(self, u, kernel_key, terms):
        """Exact XOR-convolution transfer numerator: sum_T a_T f(S xor T)."""
        key = (u, kernel_key)
        hit = self.transfer_cache.get(key)
        if hit is not None:
            return hit
        self._check_time()
        parts = []
        for mask, coeff in terms:
            if coeff:
                parts.append(self.scale(self.xor_shift(u, mask), coeff))
        if not parts:
            out = self.ZERO
        else:
            # balanced summation reduces intermediate ADD growth
            while len(parts) > 1:
                nxt = []
                it = iter(parts)
                for a in it:
                    try:
                        b = next(it)
                    except StopIteration:
                        nxt.append(a)
                        break
                    nxt.append(self.add(a, b))
                parts = nxt
            out = parts[0]
        self.transfer_cache[key] = out
        self.transfer_steps += 1
        return out

    def build_query(self, assignment):
        """
        q_u(S)=1_{S subset assigned} chi_S(u).
        assignment maps variable -> bool (+1 for True, -1 for False).
        """
        cur = self.ONE
        # Build from bottom of ADD order upward.
        for v in reversed(self.order):
            if v not in assignment:
                # frequency bit must be 0
                cur = self.mk(v, cur, self.ZERO)
            elif assignment[v]:
                # factor [1, +1] => independent, canonical reduction removes it
                cur = self.mk(v, cur, cur)
            else:
                # factor [1, -1]
                cur = self.mk(v, cur, self.scale(cur, -1))
        return cur

    def eval_zero(self, u):
        """Evaluate at Fourier mask S=emptyset (all frequency bits 0)."""
        while not self.is_terminal(u):
            _, lo, _ = self.nodes[u]
            u = lo
        return self.tvalue(u)

    def reachable_count(self, root):
        seen = set()
        stack = [root]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            if not self.is_terminal(u):
                _, lo, hi = self.nodes[u]
                stack.append(lo)
                stack.append(hi)
        self.peak_reachable = max(self.peak_reachable, len(seen))
        return len(seen)


# --------------------------- Exact self-oracle ------------------------------

class FourierIFSEngine:
    def __init__(self, nvars, clauses, time_limit=TIME_LIMIT_SECONDS):
        self.original_nvars = nvars
        self.original_clauses = canonicalize_clauses(clauses)
        self.nvars, self.clauses = to_3cnf(nvars, self.original_clauses)
        self.time_limit = time_limit
        self.start = None
        self.deadline = None

        self.empty_clause = any(len(c) == 0 for c in self.clauses)
        self.variable_order = greedy_min_degree_order(self.nvars, self.clauses)
        self.rank = {v: i for i, v in enumerate(self.variable_order)}

        kernels = []
        denominator_bits = 0
        for idx, c in enumerate(self.clauses):
            terms, k = clause_walsh_kernel(c)
            denominator_bits += k
            kernels.append((idx, c, terms, k))
        self.denominator_bits = denominator_bits

        # Convolution factors commute. This deterministic order tends to keep
        # nearby variables together in the ADD variable order.
        def clause_key(item):
            _, c, _, _ = item
            ranks = [self.rank[abs(l)] for l in c] if c else [-1]
            return (max(ranks), min(ranks), len(c))
        self.kernels = tuple(sorted(kernels, key=clause_key, reverse=True))

        self.stats = {
            "harmonic_queries": 0,
            "direction_steps": 0,
            "transformed_variables": self.nvars,
            "transformed_clauses": len(self.clauses),
            "max_query_add_nodes": 0,
            "max_total_add_nodes": 0,
            "max_reachable_quotient_nodes": 0,
            "kernel_transfers": 0,
        }

    def _new_manager(self):
        return ADDManager(self.variable_order, self.deadline)

    def harmonic_numerator(self, assignment):
        """
        Exact numerator of H_F(assignment), with common positive denominator
        2^self.denominator_bits. Positivity/zero is therefore exact.
        """
        if self.empty_clause:
            return 0, {"reachable": 1, "nodes": 1, "transfers": 0}
        mgr = self._new_manager()
        root = mgr.build_query(assignment)
        mgr.reachable_count(root)
        for kernel_id, _c, terms, _k in self.kernels:
            root = mgr.apply_kernel(root, kernel_id, terms)
            mgr.reachable_count(root)
        num = mgr.eval_zero(root)
        self.stats["harmonic_queries"] += 1
        self.stats["max_query_add_nodes"] = max(self.stats["max_query_add_nodes"], len(mgr.nodes))
        self.stats["max_total_add_nodes"] = max(self.stats["max_total_add_nodes"], len(mgr.nodes))
        self.stats["max_reachable_quotient_nodes"] = max(
            self.stats["max_reachable_quotient_nodes"], mgr.peak_reachable
        )
        self.stats["kernel_transfers"] += mgr.transfer_steps
        return num, {"reachable": mgr.peak_reachable, "nodes": len(mgr.nodes), "transfers": mgr.transfer_steps}

    def solve(self):
        self.start = time.perf_counter()
        self.deadline = None if self.time_limit is None else self.start + self.time_limit

        if self.empty_clause:
            return "UNSAT", None, time.perf_counter() - self.start

        assignment = {}

        # One exact global satisfiability query. No search tree.
        h0, _ = self.harmonic_numerator(assignment)
        if h0 == 0:
            return "UNSAT", None, time.perf_counter() - self.start
        if h0 < 0:
            raise RuntimeError("Invariant violation: harmonic satisfying density became negative.")

        # Guided construction: exactly one surviving prefix, no backtracking.
        for v in self.variable_order:
            # Skip variables already irrelevant only after all original/transformed vars are assigned? 
            # Query +1 first; if no satisfying completion remains, -1 must work because
            # current prefix was already certified satisfiable.
            a_plus = dict(assignment)
            a_plus[v] = True
            hp, _ = self.harmonic_numerator(a_plus)
            if hp > 0:
                assignment[v] = True
            else:
                a_minus = dict(assignment)
                a_minus[v] = False
                hm, _ = self.harmonic_numerator(a_minus)
                if hm <= 0:
                    raise RuntimeError(
                        "Harmonic direction invariant failed: neither child preserves satisfiability."
                    )
                assignment[v] = False
            self.stats["direction_steps"] += 1

        # Project model back to original variables.
        model = {v: bool(assignment.get(v, False)) for v in range(1, self.original_nvars + 1)}
        if not verify_model(self.original_clauses, model):
            raise RuntimeError("Internal error: projected SAT model failed original-CNF verification.")
        return "SAT", model, time.perf_counter() - self.start


# ------------------------------- Reporting ----------------------------------

def dimacs_model_line(model, nvars):
    return "v " + " ".join(str(v if model.get(v, False) else -v) for v in range(1, nvars + 1)) + " 0"


def solve_dimacs_text(text, filename="input.cnf", time_limit=TIME_LIMIT_SECONDS):
    nvars, clauses = parse_dimacs_text(text)
    engine = FourierIFSEngine(nvars, clauses, time_limit=time_limit)
    try:
        status, model, elapsed = engine.solve()
    except TimeoutError:
        status, model = "UNKNOWN_TIMEOUT", None
        elapsed = time.perf_counter() - engine.start if engine.start else 0.0

    lines = [
        "OMEGA Fourier-IFS SAT Engine v2",
        "================================",
        f"input_file: {filename}",
        f"status: {status}",
        f"original_variables: {nvars}",
        f"original_clauses: {len(clauses)}",
        f"3cnf_variables: {engine.nvars}",
        f"3cnf_clauses: {len(engine.clauses)}",
        f"elapsed_seconds: {elapsed:.6f}",
        "",
        "Exact Fourier-IFS quotient statistics",
        "-------------------------------------",
    ]
    for k, v in engine.stats.items():
        lines.append(f"{k}: {v}")
    lines.extend([
        f"walsh_common_denominator_bits: {engine.denominator_bits}",
        "backtracking_branches: 0",
        "recursive_sat_calls: 0",
        "",
        "Result",
        "------",
    ])

    if status == "SAT":
        ok = verify_model(clauses, model)
        lines.append(f"verified: {str(ok).lower()}")
        lines.append(dimacs_model_line(model, nvars))
        lines.append("")
        lines.append("assignment:")
        for v in range(1, nvars + 1):
            lines.append(f"x{v}={1 if model[v] else 0}")
    elif status == "UNSAT":
        lines.append("verified_model: n/a")
        lines.append("Exact harmonic root amplitude is zero.")
    else:
        lines.append("verified_model: n/a")
        lines.append("Exact spectral quotient computation reached the configured time limit.")

    lines.extend([
        "",
        "Implementation identity",
        "-----------------------",
        "H_F(u) is evaluated exactly by the paper's backward Walsh transfer equation.",
        "The quotient is a canonical reduced ADD of the backward spectral query function.",
        "Equal ADD node ids mean exact equality of future spectral behavior.",
        "The solver follows one H_F-positive child per variable and never backtracks.",
    ])
    return "\n".join(lines) + "\n"


def colab_main():
    from google.colab import files
    print("Upload DIMACS .cnf file(s)...")
    uploaded = files.upload()
    names = [n for n in uploaded if n.lower().endswith(".cnf")]
    if not names:
        raise ValueError("No .cnf file uploaded.")
    for name in names:
        text = uploaded[name].decode("utf-8", errors="replace")
        report = solve_dimacs_text(text, filename=name)
        out = name.rsplit(".", 1)[0] + "_omega_result.txt"
        with open(out, "w", encoding="utf-8") as f:
            f.write(report)
        print("\n" + report)
        files.download(out)


try:
    import google.colab  # noqa: F401
    colab_main()
except ImportError:
    pass