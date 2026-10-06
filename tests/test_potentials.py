"""Aşama 0 kabul testleri: T1, T2, T6, T7, T8a, T8b (PROJE_TANIMI.md §4)."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import quad

from src.config import GalaxyParams, HaloParams, BlackHoleParams, R_SUN
from src.potentials import build_potential
from src.rotation import (
    circular_velocity,
    circular_velocity_numerical,
    dynamical_mass,
    kappa,
    solve_mvir_for_vc,
)
from src.units import kpcmyr_to_kms


def vc_kms(params, R):
    return kpcmyr_to_kms(circular_velocity(build_potential(params), R))


def test_T1_solar_circular_velocity():
    """T1: varsayılanlarla v_c(8.2 kpc) = 229 ± 5 km/s."""
    v = vc_kms(GalaxyParams(), R_SUN)
    assert v == pytest.approx(229.0, abs=5.0)


def test_T2_flat_rotation_curve():
    """T2: 5–20 kpc arasında v_c 205–240 km/s bandında ve düz."""
    R = np.linspace(5.0, 20.0, 151)
    v = vc_kms(GalaxyParams(), R)
    assert v.min() >= 205.0
    assert v.max() <= 240.0
    # Düzlük: band içinde toplam salınım ortalamanın %15'inden az
    assert (v.max() - v.min()) / v.mean() < 0.15


def test_T6_no_halo_keplerian_decline():
    """T6: hale kapalıyken v_c(20 kpc) < 165 km/s."""
    p = GalaxyParams(halo=HaloParams(enabled=False))
    assert "halo" not in build_potential(p).names()
    assert vc_kms(p, 20.0) < 165.0


def test_T7_black_hole_negligible_at_sun():
    """T7: M_bh 4.3e6 -> 4.3e8 iken v_c(8.2) değişimi < 1 km/s."""
    p0 = GalaxyParams()
    p1 = replace(p0, bh=BlackHoleParams(M_bh=4.3e8))
    assert abs(vc_kms(p1, R_SUN) - vc_kms(p0, R_SUN)) < 1.0


def test_T8a_dynamical_mass():
    """T8a: M_dyn = R v_c^2 / G (R = 8.2 kpc) = 1.00e11 ± 0.05e11 M_sun."""
    M = dynamical_mass(build_potential(), R_SUN)
    assert M == pytest.approx(1.00e11, abs=0.05e11)


def test_T8b_true_enclosed_mass():
    """T8b: r < 8.2 kpc küresindeki gerçek kütle = 0.90e11 ± 0.05e11 M_sun."""
    pot = build_potential()
    M_dyn = dynamical_mass(pot, R_SUN)
    M_true = pot.enclosed_mass(R_SUN)
    print(f"\nM_dyn = {M_dyn:.4e}  M_true = {M_true:.4e}  M_dyn/M_true = {M_dyn / M_true:.4f}")
    assert M_true == pytest.approx(0.90e11, abs=0.05e11)


# --------------------------------------------------------------------------
# Yardımcı doğrulamalar (analitik türev ↔ sayısal türev)
# --------------------------------------------------------------------------

def test_analytic_vs_numerical_derivative():
    pot = build_potential()
    R = np.logspace(-3, np.log10(300.0), 60)
    np.testing.assert_allclose(
        circular_velocity(pot, R), circular_velocity_numerical(pot, R, h=1e-6 * R), rtol=1e-5
    )


def test_kappa_vs_numerical():
    pot = build_potential()
    R = np.linspace(1.0, 30.0, 30)
    h = 1e-4
    om2 = lambda r: pot.dphi_dR(r, 0.0) / r
    k2_num = R * (om2(R + h) - om2(R - h)) / (2 * h) + 4 * om2(R)
    np.testing.assert_allclose(kappa(pot, R) ** 2, k2_num, rtol=1e-6)


def test_disk_total_mass_via_gauss():
    disk = build_potential()["disk"]
    assert disk.enclosed_mass(1e5) == pytest.approx(disk.M, rel=1e-3)


def test_disk_gauss_mass_equals_density_integral():
    """Gauss akı integrali = rho'nun küre hacmi üzerindeki integrali (T8b'nin tanımı)."""
    disk = build_potential()["disk"]
    r = R_SUN

    def shell(rr):  # 4 pi rr^2 <rho> : yarı küre x2, mu = cos(theta)
        f = lambda mu: disk.density(rr * np.sqrt(1.0 - mu * mu), rr * mu)
        val, _ = quad(f, 0.0, 1.0, points=[disk.b / rr] if rr > disk.b else None, limit=200)
        return 4.0 * np.pi * rr * rr * val

    M_vol, _ = quad(shell, 0.0, r, limit=200)
    assert M_vol == pytest.approx(disk.enclosed_mass(r), rel=1e-4)


def test_inverse_mode_recovers_mvir():
    """Ters mod: varsayılan v_c(8.2)'yi hedefleyince varsayılan M_vir geri gelmeli."""
    p = GalaxyParams()
    v_target = vc_kms(p, R_SUN)
    M_vir = solve_mvir_for_vc(p, v_target, R_SUN)
    assert M_vir == pytest.approx(p.halo.M_vir, rel=1e-6)
