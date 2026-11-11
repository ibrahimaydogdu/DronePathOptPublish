"""Leg energy between service nodes: minimum over geometric options (FORMULATION_D1.md section 5)."""
from __future__ import annotations
from energy import Hardware, segment_energy_wh
from geometry import City, Node, leg_options


class LegCost:
    def __init__(self, city: City, nodes: list[Node], model: str = "F", hw: Hardware = Hardware(),
                 v: float = 10.0, delta: float = 10.0, w: float = 5.0, kappa_c: float = 1.0,
                 raster=None, n_levels: int = 3, use_around: bool = True):
        self.city, self.nodes, self.model, self.hw = city, nodes, model, hw
        self.v, self.delta, self.w, self.kappa_c = v, delta, w, kappa_c
        self.raster, self.n_levels, self.use_around = raster, n_levels, use_around
        self._opts: dict[tuple[int, int], list] = {}

    def options(self, i: int, j: int):
        key = (i, j)
        if key not in self._opts:
            o = leg_options(self.city, self.nodes[i], self.nodes[j], self.delta, self.w)
            if self.use_around and self.raster is not None:
                from around import around_options
                o = o + around_options(self.city, self.raster, self.nodes[i], self.nodes[j],
                                       self.delta, self.w, self.n_levels)
            self._opts[key] = o
        return self._opts[key]

    def energy(self, i: int, j: int, load_kg: float) -> float:
        """Wh for leg i->j carrying load_kg; inf if no feasible option."""
        if i == j:
            return 0.0
        best = float("inf")
        for _, segs in self.options(i, j):
            e = 0.0
            for ell, dz in segs:
                e += segment_energy_wh(self.model, ell, dz, load_kg, self.hw, self.v, self.kappa_c)
            if e < best:
                best = e
        return best

    def best_option(self, i: int, j: int, load_kg: float):
        best, name = float("inf"), None
        for kind, segs in self.options(i, j):
            e = sum(segment_energy_wh(self.model, ell, dz, load_kg, self.hw, self.v, self.kappa_c)
                    for ell, dz in segs)
            if e < best:
                best, name = e, kind
        return name, best
