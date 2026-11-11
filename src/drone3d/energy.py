"""Energy models for D1 (FORMULATION_D1.md section 4).

Units: m, s, kg, W; energies returned in Wh. Hardware = DJI Matrice 600 Pro, TB47S set
(manual p.44). Every model = base term P_h(m) * t  +  model-specific vertical term.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

G = 9.81
RHO = 1.225


@dataclass(frozen=True)
class Hardware:
    m0: float = 9.5            # kg, six TB47S [manual p.44]
    mtow: float = 15.5         # kg [manual p.44]
    e_batt_wh: float = 599.4   # 6 x 99.9 Wh [manual p.45]
    t0_min: float = 32.0       # hover, no payload [manual p.44]
    t1_min: float = 16.0       # hover, 6 kg payload [manual p.44]
    reserve: float = 0.10      # hover times are quoted down to 10 % battery [manual p.44]
    v_up: float = 5.0          # m/s max ascent [manual p.44]
    v_down: float = 3.0        # m/s max descent [manual p.44]
    rotors: int = 6
    prop_d: float = 0.533      # m, 21 in [dji.com E2000 page]
    eta_h: float = 0.50        # overall hover efficiency (derived, sensitivity 0.45-0.55)

    @property
    def q_max(self) -> float:
        return self.mtow - self.m0

    @property
    def e_use_wh(self) -> float:
        return (1.0 - self.reserve) * self.e_batt_wh

    @property
    def p0(self) -> float:
        return self.e_use_wh / (self.t0_min / 60.0)

    @property
    def p1(self) -> float:
        return self.e_use_wh / (self.t1_min / 60.0)

    @property
    def alpha(self) -> float:
        return math.log(self.p1 / self.p0) / math.log(self.mtow / self.m0)

    @property
    def k(self) -> float:
        return self.p0 / self.m0 ** self.alpha

    @property
    def disk_area(self) -> float:
        return self.rotors * math.pi * (self.prop_d / 2.0) ** 2

    def hover_power(self, load_kg: float) -> float:
        m = self.m0 + load_kg
        return self.k * m ** self.alpha

    def ideal_hover_power(self, load_kg: float) -> float:
        m = self.m0 + load_kg
        return (m * G) ** 1.5 / math.sqrt(2.0 * RHO * self.disk_area)


MODELS = ("F", "D", "Mo", "L", "S0", "2D")


def segment_energy_wh(model: str, ell: float, dz: float, load_kg: float,
                      hw: Hardware = Hardware(), v: float = 10.0,
                      kappa_c: float = 1.0) -> float:
    """Energy (Wh) of one straight segment: horizontal length ell (m), altitude change dz (m).

    t = max(ell/v, dz+/v_up, dz-/v_down); cruise power = kappa_c * P_h(m).
    """
    if model not in MODELS:
        raise ValueError(f"unknown model {model}")
    if ell < 0:
        raise ValueError("ell must be >= 0")
    if not 0.0 <= load_kg <= hw.q_max + 1e-9:
        raise ValueError("load outside [0, Q]")
    if model == "2D":
        dz = 0.0
    up = max(dz, 0.0)
    down = max(-dz, 0.0)
    m = hw.m0 + load_kg
    ph = hw.hover_power(load_kg)
    if model == "S0":      # legacy multipliers, unsourced bound; sequential horizontal then vertical
        e_j = ph * ell / v + 2.0 * ph * up / hw.v_up + 0.1 * ph * down / 5.0
        return e_j / 3600.0
    t = max(ell / v, up / hw.v_up, down / hw.v_down)
    e_j = kappa_c * ph * t
    if model == "F":
        e_j += m * G * up / hw.eta_h
    elif model == "Mo":    # momentum theory vertical climb, calibrated by eta_h
        thrust = m * G
        v_h = math.sqrt(thrust / (2.0 * RHO * hw.disk_area))
        p_c = thrust * (hw.v_up / 2.0 + math.sqrt((hw.v_up / 2.0) ** 2 + v_h ** 2)) / hw.eta_h
        e_j += max(p_c - ph, 0.0) * up / hw.v_up
    elif model == "L":     # ratios from Liu 2017 Table II (IRIS+, no payload: ascent 180 W, hover 164 W, descent 150 W) = +9.8 % / -8.5 %; transfer to M600 Pro is an ASSUMPTION
        e_j += 0.098 * ph * up / hw.v_up - 0.085 * ph * down / hw.v_down
    # model "D": Dorling 2017, power ~ hover power in level flight and altitude change: no extra term
    return e_j / 3600.0
