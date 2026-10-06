"""Aşama 1 kabul testi: T5 (CLAUDE.md §4)."""

import numpy as np

from src.config import NumericalParams, OrbitParams
from src.orbits import energy, initial_conditions, leapfrog
from src.potentials import build_potential

N_STEPS = 1000


def test_T5_leapfrog_energy_1000_steps():
    """T5: leapfrog ile 1000 adımda toplam enerjinin bağıl sapması < 1e-3."""
    pot = build_potential()
    o = OrbitParams()
    pos, vel = initial_conditions(pot, o.R0, o.v_phi_frac, o.v_R0, o.v_z0)
    orb = leapfrog(pot, pos, vel, NumericalParams().dt, N_STEPS)
    E = energy(pot, orb.pos[:, 0], orb.vel[:, 0])
    dev = np.max(np.abs(E - E[0])) / abs(E[0])
    print(f"\nT5: max |dE/E| = {dev:.2e}  ({N_STEPS} adım, dt = {NumericalParams().dt} Myr)")
    assert dev < 1e-3


def test_energy_error_bounded_long_run():
    """Leapfrog simplektik: 20 kat uzun koşuda hata büyümemeli (sürüklenme yok)."""
    pot = build_potential()
    o = OrbitParams()
    pos, vel = initial_conditions(pot, o.R0, o.v_phi_frac, o.v_R0, o.v_z0)
    orb = leapfrog(pot, pos, vel, NumericalParams().dt, 20 * N_STEPS)
    E = energy(pot, orb.pos[:, 0], orb.vel[:, 0])
    rel = np.abs(E - E[0]) / abs(E[0])
    first, last = rel[: N_STEPS].max(), rel[-N_STEPS:].max()
    assert last < 3.0 * first
