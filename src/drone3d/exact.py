"""Exact solver for small D1 instances (FORMULATION_D1.md section 7).

Model: single drone, sequential trips from the depot (node 0), battery refreshed at every
depot visit. Request types: 'D' deliver depot->d, 'P' pick up o->depot, 'B' pick up at o and
deliver at d in the same trip. Load <= q_max after every event, trip energy <= e_use.

Algorithm: (1) for every subset S of requests the cheapest single trip c(S) by depth-first
search over event orders with capacity / energy / bound pruning (start load = sum of D masses
in S is known), (2) set-partition DP f(S) = min_T c(T) + f(S \ T). Exact for the model.
"""
from __future__ import annotations
import itertools
from dataclasses import dataclass

INF = float("inf")


@dataclass(frozen=True)
class Request:
    kind: str       # 'D' | 'P' | 'B'
    o: int          # origin node index (depot 0 for D)
    d: int          # destination node index (depot 0 for P)
    q: float        # kg


def trip_events(reqs: list[Request], S: tuple[int, ...]):
    ev = []
    for r in S:
        q = reqs[r]
        if q.kind == "D":
            ev.append((r, "d", q.d, -q.q))
        elif q.kind == "P":
            ev.append((r, "p", q.o, +q.q))
        else:
            ev.append((r, "p", q.o, +q.q))
            ev.append((r, "d", q.d, -q.q))
    return ev


def best_trip(leg, reqs, S, q_max, e_use, t_service_energy=None):
    """Return (cost, order) of the cheapest single trip serving requests S, or (INF, None).

    leg.energy(i, j, load_kg) -> Wh. t_service_energy(load_after) -> Wh spent hovering at a stop.
    """
    ev = trip_events(reqs, S)
    load0 = sum(reqs[r].q for r in S if reqs[r].kind == "D")
    if load0 > q_max + 1e-9:
        return INF, None
    n = len(ev)
    best = [INF, None]

    def prec_ok(k, done):
        r, typ, _, _ = ev[k]
        if typ == "d" and reqs[r].kind == "B":
            # pickup event of same request must be done
            for kk in range(n):
                if ev[kk][0] == r and ev[kk][1] == "p":
                    return bool(done >> kk & 1)
        return True

    def dfs(pos, load, energy, done, order):
        if energy >= best[0] or energy > e_use + 1e-9:
            return
        if done == (1 << n) - 1:
            back = leg.energy(pos, 0, max(load, 0.0))
            tot = energy + back
            if tot < best[0] and tot <= e_use + 1e-9:
                best[0], best[1] = tot, list(order)
            return
        for k in range(n):
            if done >> k & 1 or not prec_ok(k, done):
                continue
            r, typ, node, dq = ev[k]
            new_load = load + dq
            if new_load > q_max + 1e-9 or new_load < -1e-9:
                continue
            e = energy + leg.energy(pos, node, max(load, 0.0))
            if t_service_energy is not None:
                e += t_service_energy(max(new_load, 0.0))
            order.append(k)
            dfs(node, new_load, e, done | (1 << k), order)
            order.pop()

    dfs(0, load0, 0.0, 0, [])
    if best[1] is None:
        return INF, None
    return best[0], [(ev[k][0], ev[k][1], ev[k][2]) for k in best[1]]


def solve_exact(leg, reqs: list[Request], q_max: float, e_use: float, t_service_energy=None):
    """Return (total energy, list of trips) with trip = (cost, event list)."""
    n = len(reqs)
    full = (1 << n) - 1
    cost = {}
    order = {}
    for mask in range(1, full + 1):
        S = tuple(i for i in range(n) if mask >> i & 1)
        c, o = best_trip(leg, reqs, S, q_max, e_use, t_service_energy)
        cost[mask], order[mask] = c, o
    f = [INF] * (full + 1)
    choice = [0] * (full + 1)
    f[0] = 0.0
    for mask in range(1, full + 1):
        low = mask & -mask
        sub = mask
        while sub:
            if sub & low and cost[sub] < INF and f[mask ^ sub] < INF:
                v = cost[sub] + f[mask ^ sub]
                if v < f[mask]:
                    f[mask], choice[mask] = v, sub
            sub = (sub - 1) & mask
    if f[full] == INF:
        return INF, []
    trips, m = [], full
    while m:
        s = choice[m]
        trips.append((cost[s], order[s]))
        m ^= s
    return f[full], trips


def brute_force(leg, reqs: list[Request], q_max: float, e_use: float, t_service_energy=None):
    """Independent check for tiny instances: all event permutations x all split patterns."""
    events = []
    for r, q in enumerate(reqs):
        if q.kind == "D":
            events.append((r, "d", q.d, -q.q))
        elif q.kind == "P":
            events.append((r, "p", q.o, +q.q))
        else:
            events.append((r, "p", q.o, +q.q))
            events.append((r, "d", q.d, -q.q))
    n = len(events)
    best = INF

    def valid_perm(perm):
        seen_p = set()
        for k in perm:
            r, typ, _, _ = events[k]
            if typ == "p" and reqs[r].kind == "B":
                seen_p.add(r)
            if typ == "d" and reqs[r].kind == "B" and r not in seen_p:
                return False
        return True

    def trip_cost(sub):
        # (C1): both events of a B request must be in the same trip
        present = {}
        for k in sub:
            r, typ = events[k][0], events[k][1]
            if reqs[r].kind == "B":
                present.setdefault(r, set()).add(typ)
        if any(len(v) != 2 for v in present.values()):
            return INF
        load = sum(-events[k][3] for k in sub if events[k][1] == "d" and reqs[events[k][0]].kind == "D")
        if load > q_max + 1e-9:
            return INF
        pos, e = 0, 0.0
        for k in sub:
            _, _, node, dq = events[k]
            e += leg.energy(pos, node, max(load, 0.0))
            load += dq
            if load > q_max + 1e-9 or load < -1e-9:
                return INF
            if t_service_energy is not None:
                e += t_service_energy(max(load, 0.0))
            pos = node
        e += leg.energy(pos, 0, max(load, 0.0))
        return e if e <= e_use + 1e-9 else INF

    for perm in itertools.permutations(range(n)):
        if not valid_perm(perm):
            continue
        for cuts in range(1 << (n - 1)):
            total, start, ok = 0.0, 0, True
            for i in range(n):
                end_here = (i == n - 1) or (cuts >> i & 1)
                if end_here:
                    c = trip_cost(perm[start:i + 1])
                    if c == INF:
                        ok = False
                        break
                    total += c
                    start = i + 1
            if ok and total < best:
                best = total
    return best
