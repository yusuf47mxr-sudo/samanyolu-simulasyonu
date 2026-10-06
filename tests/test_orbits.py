"""Aşama 1 kabul testleri: T3, T4 (CLAUDE.md §4)."""

import numpy as np
import pytest

from src.config import NumericalParams, OrbitParams, R_SUN
from src.orbits import angular_momentum_z, initial_conditions, leapfrog
from src.potentials import build_potential
from src.rotation import orbital_period


def test_T3_solar_orbital_period():
    """T3: Güneş'in dairesel yörünge periyodu 2 pi R / v_c = 220 ± 10 Myr."""
    T = float(orbital_period(build_potential(), R_SUN))
    print(f"\nT3: 2 pi R / v_c = {T:.2f} Myr")
    assert T == pytest.approx(220.0, abs=10.0)


def test_T3_integrated_period_matches_analytic():
    """Leapfrog ile entegre edilen dairesel yörüngenin periyodu analitik değerle uyuşmalı."""
    pot = build_potential()
    T = float(orbital_period(pot, R_SUN))
    dt = NumericalParams().dt
    pos, vel = initial_conditions(pot, R_SUN, 1.0)
    orb = leapfrog(pot, pos, vel, dt, int(1.2 * T / dt))
    phi = np.unwrap(np.arctan2(orb.pos[:, 0, 1], orb.pos[:, 0, 0]))
    T_num = np.interp(2.0 * np.pi, phi, orb.t)
    R = np.hypot(orb.pos[:, 0, 0], orb.pos[:, 0, 1])
    print(f"\nT3: entegre periyot = {T_num:.2f} Myr, R sapması = {np.ptp(R):.2e} kpc")
    assert T_num == pytest.approx(T, rel=1e-3)


def test_T4_Lz_conservation():
    """T4: eksenel simetrik potansiyelde eksantrik, eğik yörüngede L_z bağıl sapması < 1e-6."""
    pot = build_potential()
    o = OrbitParams()
    pos, vel = initial_conditions(pot, o.R0, o.v_phi_frac, o.v_R0, o.v_z0)
    dt = NumericalParams().dt
    orb = leapfrog(pot, pos, vel, dt, int(o.t_orbit / dt))
    Lz = angular_momentum_z(orb.pos[:, 0], orb.vel[:, 0])
    dev = np.max(np.abs(Lz - Lz[0])) / abs(Lz[0])
    print(f"\nT4: max |dL_z/L_z| = {dev:.2e}  ({len(orb.t) - 1} adım)")
    assert dev < 1e-6
