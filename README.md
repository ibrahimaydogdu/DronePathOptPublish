# 3D Drone Logistics Routing and Energy Optimization — Replication Package

This repository provides the official source code, benchmark data, and verification suites for 3D energy-minimizing multi-trip drone pickup-and-delivery routing under realistic urban obstacle geometries.

---

## 1. Overview and Architecture

- **`src/drone3d/`**: Core algorithmic and physical modelling modules:
  - `energy.py`: Aerodynamic and physical energy formulations (Momentum theory, Potential energy, Zero-climb, Empirical small UAV, Legacy).
  - `geometry.py` & `around.py`: 3D obstacle clearance checking, 8-connected raster route search, and facade node generation.
  - `exact.py`: Trip enumeration and set-partitioning dynamic programming exact solver for verification on small instances.
  - `heuristics.py`: Random-key metaheuristics (ABC, GWO, TLBO) with precedence repair and optimal trip splitting.
  - `campaign.py`: Multi-run simulation campaign framework.
  - `study.py` & `decompose.py`: Factorial sensitivity analysis and energy loss component decomposition.
- **`tests/`**: Comprehensive automated test suite (70 unit and integration tests).
- **`data/`**: Real-world Manhattan building footprint geometry and land-use data (NYC Open Data, DoITT/OTI).
- **`results/`**: Numerical logs, exact validation outputs, and summary tables.

---

## 2. Requirements and Installation

Python 3.10 or higher is required.

```bash
git clone https://github.com/ibrahimaydogdu/DronePathOptPublish.git
cd DronePathOptPublish
pip install -r requirements.txt
```

---

## 3. Running Tests and Verification

To verify the test suite:
```bash
pytest tests
```

To run the exact solver verification on small instances:
```bash
python -m src.drone3d.validate_exact_mo
```

To execute metaheuristic comparison on benchmark instances:
```bash
python -m src.drone3d.compare_meta_small
```

---

## 4. License and Peer Review Notice
This software is released under the [MIT License](LICENSE).  
For double-blind peer review, repository and author metadata can be viewed anonymously via:  
**Anonymous GitHub:** [https://anonymous.4open.science/r/DronePathOptPublish-EDB4/](https://anonymous.4open.science/r/DronePathOptPublish-EDB4/)
