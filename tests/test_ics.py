"""Aşama 2a: N-cisim disk başlangıç koşulları doğrulaması."""

from dataclasses import replace

import numpy as np
import pytest

from src.config import DiskICParams, GalaxyParams, NumericalParams
from src.ics import disk_profiles, make_disk, measure_profiles, measured_Q
from src.potentials import FreemanDisk
from src.units import kpcmyr_to_kms

GP, IC, NUM = GalaxyParams(), DiskICParams(), NumericalParams()
EDGES = np.linspace(0.0, IC.R_cut, 31)        # 0.5 kpc'lik halkalar
INNER, OUTER = 1.0, 14.0                     # kpc; merkez ve kesik kenarı dışarıda tutulur


@pytest.fixture(scope="module")
def disk():
    return make_disk(GP, IC, NUM)


@pytest.fixture(scope="module")
def meas(disk):
    m = measure_profiles(disk, EDGES)
    m["sel"] = (m["R_mid"] > INNER) & (m["R_mid"] < OUTER)
    return m


def test_mass_and_center(disk):
    assert disk.mass.sum() == pytest.approx(IC.M_d, rel=1e-12)
    assert np.allclose(np.average(disk.pos, axis=0, weights=disk.mass), 0.0, atol=1e-12)
    assert np.allclose(np.average(disk.vel, axis=0, weights=disk.mass), 0.0, atol=1e-15)
    assert np.hypot(disk.pos[:, 0], disk.pos[:, 1]).max() <= IC.R_cut + 0.05


def test_surface_density_exponential(meas):
    """Ölçülen Sigma(R)'ye üstel uyum: R_d = 2.6 kpc ± %5."""
    s = meas["sel"]
    slope = np.polyfit(meas["R_mid"][s], np.log(meas["Sigma"][s]), 1)[0]
    R_d_fit = -1.0 / slope
    print(f"\nölçülen R_d = {R_d_fit:.3f} kpc (hedef {IC.R_d})")
    assert R_d_fit == pytest.approx(IC.R_d, rel=0.05)


def test_vertical_profile(disk):
    """sech^2: |z| medyanı z_0 artanh(1/2) = 0.549 z_0."""
    assert np.median(np.abs(disk.pos[:, 2])) == pytest.approx(IC.z_0 * np.arctanh(0.5), rel=0.03)


@pytest.mark.parametrize("Q_target", [IC.Q_target, 0.8])
def test_measured_Q_matches_target(Q_target):
    ic = replace(IC, Q_target=Q_target)
    m = measure_profiles(make_disk(GP, ic, NUM), EDGES)
    s = (m["R_mid"] > INNER) & (m["R_mid"] < OUTER)
    Q = measured_Q(GP, ic, m)[s]
    print(f"\nQ_target={Q_target}: ölçülen Q ortalama {Q.mean():.3f}, aralık {Q.min():.2f}–{Q.max():.2f}")
    assert Q.mean() == pytest.approx(Q_target, rel=0.05)
    assert np.all(np.abs(Q / Q_target - 1.0) < 0.25)


def test_mean_rotation_matches_analytic(disk, meas):
    """Ölçülen <v_phi> = parçacık v_c'si - asimetrik sürüklenme, 5 km/s içinde.

    v_c^2 = analitik (ince disk) + parçacıklardan ölçülen düzeltme (kalınlık, yumuşatma).
    """
    prof = disk_profiles(GP, IC, meas["R_mid"])
    s = meas["sel"]
    vc2 = prof["v_c"] ** 2 + disk.vc2_correction(meas["R_mid"])
    expected = np.sqrt(vc2 + prof["sigma_R"] ** 2 * prof["drift_term"])
    diff = kpcmyr_to_kms(meas["vphi_mean"][s] - expected[s])
    drift = kpcmyr_to_kms(np.sqrt(vc2[s]) - expected[s])
    print(f"\nmax |<v_phi>_ölçülen - <v_phi>_beklenen| = {np.abs(diff).max():.2f} km/s; "
          f"sürüklenme {drift.min():.1f}–{drift.max():.1f} km/s")
    assert np.all(np.abs(diff) < 5.0)
    assert not prof["vphi2_negative"][s].any()


def test_freeman_derivatives_vs_numerical():
    d = FreemanDisk(IC.M_d, IC.R_d, IC.R_cut)
    R, h = np.linspace(0.5, 20.0, 40), 1e-5
    np.testing.assert_allclose((d.phi(R + h) - d.phi(R - h)) / (2 * h), d.dphi_dR(R), rtol=1e-6)
    np.testing.assert_allclose((d.dphi_dR(R + h) - d.dphi_dR(R - h)) / (2 * h),
                               d.d2phi_dR2_midplane(R), rtol=1e-6)


def test_config_disk_mass_is_solver_result():
    """config.DiskICParams.M_d, v_c(8.2) = 228.7 km/s çözümünün kendisi olmalı (bayatlamasın)."""
    from src.config import R_SUN, VC_SUN_TARGET
    from src.ics import solve_disk_mass
    M = solve_disk_mass(GP, IC, VC_SUN_TARGET, R_SUN)
    print(f"\nçözülen M_d = {M:.5e} M_sun, config = {IC.M_d:.5e}")
    assert IC.M_d == pytest.approx(M, rel=1e-4)


def test_eln_parameter():
    """Efstathiou–Lake–Negroponte: varsayılan disk ~0.80, < 1.1 (çubuk oluşur)."""
    from src.ics import eln_parameter
    e = eln_parameter(GP, IC)
    print(f"\nELN eps = {e:.3f}")
    assert e == pytest.approx(0.80, abs=0.05)
    assert e < 1.1


def test_T9_particle_rotation_velocity(disk):
    """T9: üretilen diskin parçacıklarından ölçülen v_c(8.2) = 229 ± 5 km/s.

    v_c, parçacıkların kütleçekimi (doğrudan toplam, Plummer eps) + sabit bileşenlerin
    R = 8.2 kpc, z = 0 halkasındaki radyal ivmesinden ölçülür: v_c^2 = R <-a_R>.
    """
    from src.config import R_SUN
    from src.nbody import direct_accel
    from src.orbits import acceleration
    from src.potentials import build_external_potential

    phi = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    pts = np.column_stack((R_SUN * np.cos(phi), R_SUN * np.sin(phi), np.zeros_like(phi)))
    a = direct_accel(pts, disk.pos, disk.mass, NUM.eps, False) \
        + acceleration(build_external_potential(GP, NUM.eps), pts)
    aR = -(a[:, 0] * np.cos(phi) + a[:, 1] * np.sin(phi))
    v = kpcmyr_to_kms(np.sqrt(R_SUN * aR.mean()))
    m = measure_profiles(disk, np.array([R_SUN - 0.25, R_SUN + 0.25]))
    print(f"\nT9: parçacıklardan v_c(8.2) = {v:.2f} km/s; "
          f"(bilgi: <v_phi>(8.2) = {kpcmyr_to_kms(m['vphi_mean'][0]):.2f} km/s, asimetrik sürüklenme dahil)")
    assert v == pytest.approx(229.0, abs=5.0)


def test_asymmetric_drift_matches_jeans(disk):
    """v_c - <v_phi>: parçacıklardan ölçülen = Jeans denkleminin öngördüğü, %20 içinde.

    R = 8.2 ± 1 kpc (~2100 parçacık; istatistik hata ~%4).
    """
    from src.config import R_SUN
    from src.ics import asymmetric_drift_check
    r = asymmetric_drift_check(GP, IC, disk, R_SUN, 1.0)
    print(f"\nasimetrik sürüklenme @8.2 kpc: ölçülen {r['measured']:.2f} ± {r['measured_err']:.2f} km/s, "
          f"öngörülen {r['predicted']:.2f} km/s (1. mertebe {r['predicted_first_order']:.2f}); "
          f"v_c = {r['v_c']:.2f}, <v_phi> = {r['vphi_mean']:.2f}")
    assert r["measured"] == pytest.approx(r["predicted"], rel=0.20)
