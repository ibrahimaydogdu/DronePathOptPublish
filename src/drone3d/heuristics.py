"""Common decoder and the three selected metaheuristics (TLBO, GWO, ABC).

User decision (2026-11-08): optimisation techniques are tools; TLBO, GWO and ABC are taken from the
user's OptimizationMethods.vb library (lines 297-337 ABC, 544-581 TLBO, 614-643 GWO) and ported to a
continuous random-key search space [0,1]^D, D = number of events. Same decoder, same initial
population for a given seed, equal evaluation budget (decoder calls).
"""
from __future__ import annotations
import math
import random
from dataclasses import dataclass

import numpy as np

from exact import INF, Request, trip_events

PENALTY = 1.0e6


class Problem:
    """Event sequence decoder: random keys -> order -> precedence repair -> optimal consecutive trip split."""

    def __init__(self, leg, reqs: list[Request], q_max: float, e_use: float, max_trip_events: int = 12):
        self.leg, self.reqs, self.q_max, self.e_use = leg, reqs, q_max, e_use
        self.events = trip_events(reqs, tuple(range(len(reqs))))   # (r, typ, node, dq)
        self.D = len(self.events)
        self.kmax = max_trip_events
        self._memo = {}
        self.pickup_of = {}
        for k, (r, typ, _, _) in enumerate(self.events):
            if typ == "p" and reqs[r].kind == "B":
                self.pickup_of[r] = k
        self.evals = 0

    # ---- single trip cost for a consecutive piece of the sequence
    def trip_cost(self, piece: tuple[int, ...]) -> float:
        if piece in self._memo:
            return self._memo[piece]
        res = self._trip_cost(piece)
        self._memo[piece] = res
        return res

    def _trip_cost(self, piece):
        present = {}
        for k in piece:
            r, typ = self.events[k][0], self.events[k][1]
            if self.reqs[r].kind == "B":
                present.setdefault(r, set()).add(typ)
        if any(len(v) != 2 for v in present.values()):
            return INF
        load = sum(self.reqs[self.events[k][0]].q for k in piece
                   if self.events[k][1] == "d" and self.reqs[self.events[k][0]].kind == "D")
        if load > self.q_max + 1e-9:
            return INF
        pos, e = 0, 0.0
        for k in piece:
            _, _, node, dq = self.events[k]
            e += self.leg.energy(pos, node, max(load, 0.0))
            load += dq
            if load > self.q_max + 1e-9 or load < -1e-9 or e > self.e_use + 1e-9:
                return INF
            pos = node
        e += self.leg.energy(pos, 0, max(load, 0.0))
        return e if e <= self.e_use + 1e-9 else INF

    # ---- decoding
    def sequence(self, x) -> list[int]:
        seq = [int(k) for k in np.argsort(x, kind="stable")]
        for r, kp in self.pickup_of.items():
            kd = kp + 1  # delivery event follows its pickup in the event list
            ip, idl = seq.index(kp), seq.index(kd)
            if idl < ip:
                seq.pop(idl)
                seq.insert(seq.index(kp) + 1, kd)
        return seq

    def decode(self, x):
        """Return (total energy or penalty, list of pieces)."""
        seq = self.sequence(x)
        L = len(seq)
        f = [INF] * (L + 1)
        arg = [0] * (L + 1)
        f[0] = 0.0
        for i in range(1, L + 1):
            for j in range(i - 1, max(i - self.kmax, 0) - 1, -1):
                if f[j] == INF:
                    continue
                c = self.trip_cost(tuple(seq[j:i]))
                if c < INF and f[j] + c < f[i]:
                    f[i], arg[i] = f[j] + c, j
        if f[L] == INF:
            return PENALTY + L, []
        pieces, i = [], L
        while i > 0:
            pieces.append(tuple(seq[arg[i]:i]))
            i = arg[i]
        return f[L], pieces[::-1]

    def evaluate(self, x) -> float:
        self.evals += 1
        return self.decode(x)[0]


@dataclass
class Result:
    best: float
    x: np.ndarray
    history: list   # (evaluations, best so far)
    evals: int
    extra: dict = None


class Budget:
    def __init__(self, problem: Problem, limit: int):
        self.p, self.limit, self.hist = problem, limit, []
        self.best, self.bx = INF, None
        self.start = problem.evals

    @property
    def used(self):
        return self.p.evals - self.start

    @property
    def done(self):
        return self.used >= self.limit

    def f(self, x):
        v = self.p.evaluate(x)
        if v < self.best:
            self.best, self.bx = v, np.array(x, copy=True)
            self.hist.append((self.used, v))
        return v


def init_population(D: int, N: int, seed: int) -> np.ndarray:
    return np.random.default_rng(np.random.SeedSequence([seed, 0])).random((N, D))


def tlbo(problem: Problem, limit: int, seed: int, N: int = 30, pop0=None, coords=None) -> Result:
    """coords=None: library operator (all coordinates). coords=1: ablation, one random coordinate per move."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, 1]))
    B = Budget(problem, limit)
    X = (pop0 if pop0 is not None else init_population(problem.D, N, seed)).copy()
    N = len(X)
    cost = np.array([B.f(x) for x in X])
    while not B.done:
        teacher = X[int(np.argmin(cost))].copy()
        mean = X.mean(axis=0)
        for i in range(N):
            if B.done:
                break
            tf = round(1 + rng.random())
            xn = np.clip(X[i] + rng.random(problem.D) * (teacher - tf * mean), 0.0, 1.0)
            if coords == 1:
                d = int(rng.integers(problem.D)); v = X[i].copy(); v[d] = xn[d]; xn = v
            fn = B.f(xn)
            if fn < cost[i]:
                X[i], cost[i] = xn, fn
            if B.done:
                break
            j = int(rng.integers(N - 1)); j += j >= i
            better = cost[i] < cost[j]
            r = rng.random(problem.D)
            xn = np.clip(X[i] + r * ((X[i] - X[j]) if better else (X[j] - X[i])), 0.0, 1.0)
            if coords == 1:
                d = int(rng.integers(problem.D)); v = X[i].copy(); v[d] = xn[d]; xn = v
            fn = B.f(xn)
            if fn < cost[i]:
                X[i], cost[i] = xn, fn
    return Result(B.best, B.bx, B.hist, B.used)


def gwo(problem: Problem, limit: int, seed: int, N: int = 30, pop0=None, a_start: float = 2.0) -> Result:
    rng = np.random.default_rng(np.random.SeedSequence([seed, 2]))
    B = Budget(problem, limit)
    X = (pop0 if pop0 is not None else init_population(problem.D, N, seed)).copy()
    N = len(X)
    cost = np.array([B.f(x) for x in X])
    order = np.argsort(cost)
    leaders = [(float(cost[k]), X[k].copy()) for k in order[:3]]
    while len(leaders) < 3:
        leaders.append(leaders[0])
    while not B.done:
        for i in range(N):
            if B.done:
                break
            a = a_start * (1.0 - B.used / limit)
            xn = np.zeros(problem.D)
            for _, L in leaders:
                A = 2 * a * rng.random(problem.D) - a
                C = 2 * rng.random(problem.D)
                xn += L - A * np.abs(C * L - X[i])
            xn = np.clip(xn / 3.0, 0.0, 1.0)
            fn = B.f(xn)
            X[i], cost[i] = xn, fn          # library GWO: replacement is not greedy
            key = tuple(np.round(xn, 6))
            if any(tuple(np.round(L, 6)) == key for _, L in leaders):
                continue
            for k in range(3):
                if fn < leaders[k][0]:
                    leaders.insert(k, (fn, xn.copy()))
                    leaders.pop()
                    break
    return Result(B.best, B.bx, B.hist, B.used)


def abc(problem: Problem, limit: int, seed: int, N: int = 30, pop0=None, limit_trials: int = 0, coords=None) -> Result:
    """coords=None: library operator (one coordinate). coords='all': ablation, all coordinates perturbed."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, 3]))
    B = Budget(problem, limit)
    X = (pop0 if pop0 is not None else init_population(problem.D, N, seed)).copy()
    N, D = len(X), problem.D
    cost = np.array([B.f(x) for x in X])
    trials = np.zeros(N, dtype=int)
    lim = limit_trials if limit_trials > 0 else N * D
    scouts = 0

    def neighbour(i):
        k = int(rng.integers(N - 1)); k += k >= i
        v = X[i].copy()
        if coords == "all":
            return np.clip(X[i] + (2 * rng.random(D) - 1) * (X[i] - X[k]), 0.0, 1.0)
        d = int(rng.integers(D))
        v[d] = np.clip(X[i][d] + (2 * rng.random() - 1) * (X[i][d] - X[k][d]), 0.0, 1.0)
        return v

    while not B.done:
        for i in range(N):
            if B.done:
                break
            v = neighbour(i)
            fv = B.f(v)
            if fv < cost[i]:
                X[i], cost[i], trials[i] = v, fv, 0
            else:
                trials[i] += 1
            if B.done:
                break
            fit = 1.0 / (1.0 + np.where(np.isfinite(cost), np.maximum(cost, 0.0), 1e12))
            k = int(rng.choice(N, p=fit / fit.sum()))
            v = neighbour(k)
            fv = B.f(v)
            if fv < cost[k]:
                X[k], cost[k], trials[k] = v, fv, 0
            else:
                trials[k] += 1
            if B.done:
                break
            if trials[i] > lim:
                v = rng.random(D)
                fv = B.f(v)
                X[i], cost[i], trials[i] = v, fv, 0   # scout: unconditional replacement
                scouts += 1
    return Result(B.best, B.bx, B.hist, B.used, {"scouts": scouts})


ALGORITHMS = {"TLBO": tlbo, "GWO": gwo, "ABC": abc}


def random_search(problem, limit: int, seed: int, N: int = 30, pop0=None) -> Result:
    """Baseline: uniform random keys (no learning)."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, 4]))
    B = Budget(problem, limit)
    while not B.done:
        B.f(rng.random(problem.D))
    return Result(B.best, B.bx, B.hist, B.used)


def hill_climb(problem, limit: int, seed: int, N: int = 30, pop0=None, patience: int = 200, sigma: float = 0.1) -> Result:
    """Baseline: stochastic hill climbing with random restarts (one-dimension Gaussian move)."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, 5]))
    B = Budget(problem, limit)
    D = problem.D
    while not B.done:
        x = rng.random(D)
        fx = B.f(x)
        stall = 0
        while not B.done and stall < patience:
            y = x.copy()
            d = int(rng.integers(D))
            y[d] = np.clip(y[d] + sigma * rng.standard_normal(), 0.0, 1.0)
            fy = B.f(y)
            if fy <= fx:
                stall = 0 if fy < fx else stall + 1
                x, fx = y, fy
            else:
                stall += 1
    return Result(B.best, B.bx, B.hist, B.used)


BASELINES = {"RS": random_search, "HC": hill_climb}


def nn_two_opt(problem, limit: int, seed: int = 0, N: int = 30, pop0=None) -> Result:
    """Deterministic constructive baseline: nearest-event construction + 2-opt (segment reversal) on the
    event sequence, accepted by the decoder objective. Precedence (pickup before delivery of a B request)
    is respected. Uses the load-0 energy table as proximity measure. Ignores `seed` (deterministic)."""
    B = Budget(problem, limit)
    D = problem.D
    node = problem.node
    E0 = problem.E[:, :, 0]
    pick = {kp: kp + 1 for kp in problem.pickup_of.values()}      # pickup event -> its delivery event
    deliv_to_pick = {v: k for k, v in pick.items()}
    remaining = set(range(D))
    pos, seq = 0, []
    while remaining:
        ok = [k for k in remaining if k not in deliv_to_pick or deliv_to_pick[k] not in remaining]
        k = min(ok, key=lambda e: (E0[pos, node[e]], e))
        seq.append(k)
        remaining.discard(k)
        pos = node[k]

    def key_of(s):
        x = np.empty(D)
        for rank, k in enumerate(s):
            x[k] = (rank + 0.5) / D
        return x

    def valid(s):
        where = {k: i for i, k in enumerate(s)}
        return all(where[kp] < where[kd] for kp, kd in pick.items())

    cur = B.f(key_of(seq))
    improved = True
    while improved and not B.done:
        improved = False
        for i in range(D - 1):
            for j in range(i + 1, D):
                if B.done:
                    break
                cand = seq[:i] + seq[i:j + 1][::-1] + seq[j + 1:]
                if not valid(cand):
                    continue
                v = B.f(key_of(cand))
                if v < cur - 1e-12:
                    seq, cur, improved = cand, v, True
    return Result(B.best, B.bx, B.hist, B.used)


BASELINES["NN2OPT"] = nn_two_opt
