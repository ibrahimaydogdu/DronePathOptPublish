"""Fast table-based decoder (replaces heuristics.Problem for campaigns).

Why: Review_Round_03 found (and a real-data check confirmed) that the first decoder returned a flat
penalty for infeasible sequences, so random keys were feasible 75 % / 17 % / 1 % of the time for
n = 8 / 15 / 25 requests. This decoder solves a RELAXED consecutive trip split by dynamic programming
with penalties that grow with the violation (battery excess, capacity excess, a B request split over
two trips), so the landscape has a gradient. A solution is feasible iff its penalty is zero.

Energies come from a precomputed table E[i, j, l] (Wh, load l in 0.1 kg steps) = LegCost.energy.
"""
from __future__ import annotations
import numpy as np
from numba import njit

from exact import Request, trip_events

STEP = 0.1                 # kg, masses are rounded to 0.1 kg
W_BATT = 20.0              # penalty per Wh above e_use   (decoder constant, same for all algorithms)
W_CAP = 300.0              # penalty per kg above q_max
W_PAIR = 500.0             # penalty per unpaired B event inside a trip piece


def build_table(leg, n_nodes: int, q_max: float) -> np.ndarray:
    K = int(round(q_max / STEP)) + 1
    E = np.zeros((n_nodes, n_nodes, K))
    loads = [l * STEP for l in range(K)]
    from energy import segment_energy_wh
    for i in range(n_nodes):
        for j in range(n_nodes):
            if i == j:
                continue
            opts = leg.options(i, j)
            best = np.full(K, np.inf)
            for _, segs in opts:
                e = np.array([sum(segment_energy_wh(leg.model, ell, dz, ld, leg.hw, leg.v, leg.kappa_c)
                                  for ell, dz in segs) for ld in loads])
                best = np.minimum(best, e)
            E[i, j, :] = best
    return E


def build_table_sym(leg, n_nodes: int, q_max: float) -> np.ndarray:
    """Same as build_table but computes the geometry only for i < j and reverses it for j > i."""
    from energy import segment_energy_wh
    K = int(round(q_max / STEP)) + 1
    E = np.zeros((n_nodes, n_nodes, K))
    loads = [l * STEP for l in range(K)]
    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            opts = leg.options(i, j)
            best_f = np.full(K, np.inf)
            best_r = np.full(K, np.inf)
            for _, segs in opts:
                rev = [(ell, -dz) for ell, dz in reversed(segs)]
                ef = np.array([sum(segment_energy_wh(leg.model, ell, dz, ld, leg.hw, leg.v, leg.kappa_c)
                                   for ell, dz in segs) for ld in loads])
                er = np.array([sum(segment_energy_wh(leg.model, ell, dz, ld, leg.hw, leg.v, leg.kappa_c)
                                   for ell, dz in rev) for ld in loads])
                best_f = np.minimum(best_f, ef)
                best_r = np.minimum(best_r, er)
            E[i, j, :] = best_f
            E[j, i, :] = best_r
    return E


def table_from_fn(energy_fn, n_nodes: int, q_max: float) -> np.ndarray:
    """Table from any energy(i, j, load_kg) function (used by tests with a matrix test double)."""
    K = int(round(q_max / STEP)) + 1
    E = np.zeros((n_nodes, n_nodes, K))
    for i in range(n_nodes):
        for j in range(n_nodes):
            if i != j:
                for l in range(K):
                    E[i, j, l] = energy_fn(i, j, l * STEP)
    return E


@njit(cache=True)
def _piece_cost(E, node, dqi, is_d_req, qi_req, ev_req, ev_isB, piece, plen, qmax_i, e_use, member):
    # member: scratch int array over events (0/1)
    for t in range(plen):
        member[piece[t]] = 1
    load0 = 0
    unpaired = 0
    for t in range(plen):
        k = piece[t]
        if is_d_req[ev_req[k]] == 1:
            load0 += qi_req[ev_req[k]]
        if ev_isB[k] == 1:
            partner = k + 1 if dqi[k] > 0 else k - 1
            if member[partner] == 0:
                unpaired += 1
    load = load0
    pos = 0
    e = 0.0
    pen = 0.0
    if load > qmax_i:
        pen += W_CAP * (load - qmax_i) * 0.1
        load = qmax_i
    for t in range(plen):
        k = piece[t]
        l_idx = load if load >= 0 else 0
        if l_idx > qmax_i:
            l_idx = qmax_i
        e += E[pos, node[k], l_idx]
        load += dqi[k]
        if load > qmax_i:
            pen += W_CAP * (load - qmax_i) * 0.1
        if load < 0:
            pen += W_PAIR
            load = 0
        pos = node[k]
    l_idx = load if load >= 0 else 0
    if l_idx > qmax_i:
        l_idx = qmax_i
    e += E[pos, 0, l_idx]
    if e > e_use:
        pen += W_BATT * (e - e_use)
    pen += W_PAIR * unpaired
    for t in range(plen):
        member[piece[t]] = 0
    return e, pen


@njit(cache=True)
def _decode(E, node, dqi, is_d_req, qi_req, ev_req, ev_isB, x, qmax_i, e_use, kmax, out_pieces, pen_mult):
    D = x.shape[0]
    seq = np.argsort(x, kind='mergesort')   # stable: ties broken by event index (platform independent)
    # precedence repair for B: delivery event (k+1) must follow its pickup event (k)
    for k in range(D):
        if ev_isB[k] == 1 and dqi[k] > 0:
            ip = -1
            idl = -1
            for t in range(D):
                if seq[t] == k:
                    ip = t
                if seq[t] == k + 1:
                    idl = t
            if idl < ip:
                v = seq[idl]
                for t in range(idl, D - 1):
                    seq[t] = seq[t + 1]
                # insert right after pickup (pickup index shifted by -1 if it was after idl, never: idl < ip)
                ip2 = ip - 1
                for t in range(D - 1, ip2 + 1, -1):
                    seq[t] = seq[t - 1]
                seq[ip2 + 1] = v
    f = np.full(D + 1, 1e18)
    fe = np.zeros(D + 1)       # pure energy along the best penalised path
    fp = np.zeros(D + 1)       # penalty along the best path
    arg = np.zeros(D + 1, dtype=np.int64)
    f[0] = 0.0
    member = np.zeros(D, dtype=np.int64)
    piece = np.zeros(D, dtype=np.int64)
    for i in range(1, D + 1):
        lo = i - kmax
        if lo < 0:
            lo = 0
        for j in range(i - 1, lo - 1, -1):
            plen = i - j
            for t in range(plen):
                piece[t] = seq[j + t]
            e, pen = _piece_cost(E, node, dqi, is_d_req, qi_req, ev_req, ev_isB, piece, plen, qmax_i, e_use, member)
            c = f[j] + e + pen_mult * pen
            if c < f[i]:
                f[i] = c
                fe[i] = fe[j] + e
                fp[i] = fp[j] + pen
                arg[i] = j
    # pieces for reporting
    i = D
    n_pieces = 0
    while i > 0:
        out_pieces[n_pieces] = i
        n_pieces += 1
        i = arg[i]
    return f[D], fe[D], fp[D], seq, arg, n_pieces


class FastProblem:
    """Drop-in for heuristics.Problem: D, evaluate(x) -> penalised objective, evals counter."""

    def __init__(self, E: np.ndarray, reqs: list[Request], q_max: float, e_use: float, kmax: int = 24,
                 relabel_seed: int | None = None):
        self.E, self.reqs, self.q_max, self.e_use, self.kmax = E, reqs, q_max, e_use, kmax
        ev = trip_events(reqs, tuple(range(len(reqs))))
        self.events = ev
        self.D = len(ev)
        self.node = np.array([e[2] for e in ev], dtype=np.int64)
        self.dqi = np.array([int(round(e[3] / STEP)) for e in ev], dtype=np.int64)
        self.ev_req = np.array([e[0] for e in ev], dtype=np.int64)
        self.ev_isB = np.array([1 if reqs[e[0]].kind == "B" else 0 for e in ev], dtype=np.int64)
        self.is_d_req = np.array([1 if r.kind == "D" else 0 for r in reqs], dtype=np.int64)
        self.qi_req = np.array([int(round(r.q / STEP)) for r in reqs], dtype=np.int64)
        self.qmax_i = int(round(q_max / STEP))
        self._out = np.zeros(self.D + 1, dtype=np.int64)
        # optional random relabelling of keys -> events: removes the systematic tie-break bias of index order
        self.perm = (np.random.default_rng(np.random.SeedSequence([relabel_seed, 99])).permutation(self.D)
                     if relabel_seed is not None else None)
        self.evals = 0
        self.pickup_of = {r: k for k, (r, typ, _, _) in enumerate(ev) if typ == "p" and reqs[r].kind == "B"}

    def decode(self, x, pen_mult: float = 1.0):
        xe = np.asarray(x, dtype=np.float64)
        if self.perm is not None:
            xe = xe[self.perm]
        f, e, pen, seq, arg, npieces = _decode(self.E, self.node, self.dqi, self.is_d_req, self.qi_req,
                                               self.ev_req, self.ev_isB, xe, self.qmax_i, self.e_use,
                                               self.kmax, self._out, pen_mult)
        return f, e, pen, seq, arg

    def evaluate(self, x) -> float:
        """Penalised objective (what the algorithms minimise)."""
        self.evals += 1
        return float(self.decode(x)[0])

    def feasible_energy(self, x):
        """Report: lexicographic re-solve (minimum penalty first, then energy) of the key x.

        Returns (energy, feasible flag, pieces as lists of event indices). Review_Round_04: the
        penalised objective can prefer a slightly infeasible split over a feasible one for the same key;
        reporting therefore always uses this lexicographic decode.
        """
        f, e, pen, seq, arg = self.decode(x, pen_mult=1.0e6)
        pieces, i = [], self.D
        while i > 0:
            pieces.append([int(k) for k in seq[arg[i]:i]])
            i = int(arg[i])
        return float(e), pen <= 1e-9, pieces[::-1]

    def tie_fraction(self, x) -> float:
        """Share of keys that are exactly tied with another key (diagnostic for clipping bias)."""
        x = np.asarray(x)
        return 1.0 - len(np.unique(x)) / len(x)
