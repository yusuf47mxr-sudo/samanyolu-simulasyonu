"""Aşama 3: analysis.py ölçüm fonksiyonları (sentetik veriyle) ve koşu tutarlılığı."""

import numpy as np
import pytest

from src import analysis as an

RNG = np.random.default_rng(11)
N = 200_000


def exp_disk(R_d=2.6, R_cut=15.0, z0=0.3, n=N, eps_m2=0.0):
    """Üstel disk + sech^2; eps_m2 > 0 ise Sigma ∝ (1 + eps cos 2phi) — A2/A0 = eps/2."""
    R = np.interp(RNG.uniform(0, 1 - (1 + R_cut / R_d) * np.exp(-R_cut / R_d), n),
                  1 - (1 + (g := np.linspace(0, R_cut, 20001)) / R_d) * np.exp(-g / R_d), g)
    phi = RNG.uniform(0, 2 * np.pi, n)
    if eps_m2 > 0:  # reddetme örneklemesi
        keep = RNG.uniform(0, 1 + eps_m2, n) < 1 + eps_m2 * np.cos(2 * phi)
        R, phi = R[keep], phi[keep]
    z = z0 * np.arctanh(2 * RNG.uniform(0, 1, R.size) - 1)
    snap = np.column_stack((R * np.cos(phi), R * np.sin(phi), z, np.zeros((R.size, 3))))
    return snap, np.full(R.size, 1.0)


def test_scale_length_recovers_exponential():
    s, m = exp_disk()
    assert an.scale_length(s, m) == pytest.approx(2.6, rel=0.02)


def test_thickness_sech2():
    s, _ = exp_disk()
    assert an.disk_thickness(s) == pytest.approx(0.3 * np.arctanh(0.5), rel=0.02)


def test_bar_strength_axisymmetric_and_m2():
    s, m = exp_disk()
    assert an.bar_strength_snapshot(s, m) < 5.0 / np.sqrt(N)
    s, m = exp_disk(eps_m2=0.4)
    assert an.bar_strength_snapshot(s, m) == pytest.approx(0.2, abs=0.01)
    R, A2 = an.bar_strength_profile(s, m)
    assert np.all(np.abs(A2[R < 8] - 0.2) < 0.03)


def test_bar_formation_time():
    t = np.arange(0, 2000.0, 5.0)
    A2 = np.where(t > 700, 0.35, 0.02) + 0.01 * np.sin(t)
    tb = an.bar_formation_time(t, A2)
    assert 650 <= tb <= 760
    assert an.bar_formation_time(t, np.full_like(t, 0.05)) is None


def test_clump_fraction_detects_clumps():
    s, m = exp_disk()
    smooth = an.clump_fraction(s, m)
    # 20 yoğun topak ekle (her biri 3000 parçacık, 0.1 kpc yarıçap), R = 4–10 kpc
    cl = []
    for _ in range(20):
        R0, p0 = RNG.uniform(4, 10), RNG.uniform(0, 2 * np.pi)
        c = np.column_stack((R0 * np.cos(p0) + 0.1 * RNG.normal(size=3000),
                             R0 * np.sin(p0) + 0.1 * RNG.normal(size=3000), np.zeros((3000, 4))))
        cl.append(c)
    s2 = np.vstack([s] + cl)
    clumpy = an.clump_fraction(s2, np.ones(len(s2)))
    print(f"\ntopak kesri: pürüzsüz {smooth:.4f}, topaklı {clumpy:.4f}")
    assert smooth < 0.01
    assert clumpy > 10 * smooth


def test_running_mean_constant():
    t = np.arange(0, 100.0, 5.0)
    assert np.allclose(an.running_mean(t, np.full_like(t, 3.0), 20.0), 3.0)


@pytest.mark.parametrize("name", an.list_runs("runs"))
def test_timeseries_A2_matches_snapshots(name):
    """Koşu sırasında yazılan A2/A0, snapshot'tan yeniden hesaplananla uyuşmalı (float32 duyarlılığı)."""
    run = an.load_run(name, "runs")
    for k in np.linspace(0, run.n_snap - 1, 5).astype(int):
        i = int(np.argmin(np.abs(run.ts["t"] - run.t_snap[k])))
        assert abs(run.ts["t"][i] - run.t_snap[k]) < 1e-6
        assert an.bar_strength_snapshot(run.snapshot(k), run.mass) == pytest.approx(run.ts["A2"][i], abs=1e-4)


def m2_disk(phase_of_R, eps=0.6, R_bar=None, n=N):
    """Üstel disk; Sigma ∝ 1 + eps cos 2(phi - phase(R)), R < R_bar (None: her yerde)."""
    s, m = exp_disk(n=n)
    R, phi = np.hypot(s[:, 0], s[:, 1]), np.arctan2(s[:, 1], s[:, 0])
    amp = eps if R_bar is None else np.where(R < R_bar, eps, 0.0)
    keep = RNG.uniform(0, 1 + eps, R.size) < 1 + amp * np.cos(2 * (phi - phase_of_R(R)))
    return s[keep], m[keep]


def test_bar_measure_straight_bar_length():
    """Fazı sabit, R < 4 kpc'de m=2: çubuk var, uzunluk ~4 kpc."""
    s, m = m2_disk(lambda R: 0.3 + 0 * R, R_bar=4.0)
    b = an.bar_measure(s, m)
    print(f"\ndüz çubuk: A2max={b['A2max']:.3f}, L={b['L_bar']}")
    assert b["A2max"] == pytest.approx(0.3, abs=0.03)
    assert b["L_bar"] == pytest.approx(4.0, abs=0.5)
    assert an.bar_present(b["A2max"], b["L_bar"])


def test_bar_measure_rejects_wound_spiral():
    """Sıkı sarılmış m=2 sarmal (faz R ile hızla döner): güçlü A2 ama çubuk sayılmaz."""
    s, m = m2_disk(lambda R: 1.2 * R)
    b = an.bar_measure(s, m)
    print(f"\nsarmal: A2max={b['A2max']:.3f}, L={b['L_bar']}")
    assert b["A2max"] > an.BAR_THRESHOLD
    assert not an.bar_present(b["A2max"], b["L_bar"])


def test_bar_onset_requires_persistence():
    t = np.arange(0, 1000.0, 10.0)
    present = (t >= 300) & (t < 350) | (t >= 500)   # 50 Myr'lik geçici + 500'den sonra kalıcı
    assert an.bar_onset_time(t, present) == 500.0
    assert an.bar_onset_time(t, np.zeros_like(t, bool)) is None


def test_corotation_radius_flat_curve():
    """Düz dönme eğrisinde v_c = v0: R_CR = v0 / Omega_p."""
    R = np.linspace(0.25, 15.0, 300)
    v0, om = 0.23, 0.045
    assert an.corotation_radius(R, np.full_like(R, v0), om) == pytest.approx(v0 / om, rel=1e-3)
    assert np.isnan(an.corotation_radius(R, np.full_like(R, v0), 10.0))


def test_bar_region_phase_recovers_angle():
    """Sentetik çubuk fazı 0.3 rad: çubuk bölgesinden ölçülen faz 0.3 olmalı."""
    s, m = m2_disk(lambda R: 0.3 + 0 * R, R_bar=4.0)
    bm = an.bar_measure(s, m)
    assert an.bar_region_phase(s, m, bm) == pytest.approx(0.3, abs=np.deg2rad(2.0))
