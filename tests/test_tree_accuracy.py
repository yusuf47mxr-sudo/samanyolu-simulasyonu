"""Ağaç kuvvet doğruluğu: rastgele ~500 parçacıkta ağaç vs doğrudan toplam.

Ölçülen: RMS bağıl kuvvet hatası ve en kötü tekil hata.
"""

import numpy as np
import pytest

from src.config import DiskICParams, GalaxyParams, NumericalParams
from src.ics import make_disk
from src.nbody import force_error_stats

NUM = NumericalParams()
THETAS = (0.4, 0.5, 0.6, 0.7)
N_SAMPLE = 500
SEED = 2024


@pytest.fixture(scope="module")
def disk():
    return make_disk(GalaxyParams(), DiskICParams(), NUM, verbose=False)


def error_table(disk, **kw):
    return {th: force_error_stats(disk.pos, disk.mass, NUM.eps, th, NUM.leaf_size, N_SAMPLE, SEED, **kw)
            for th in THETAS}


def print_table(title, tab):
    print(f"\n{title}\n  theta   RMS       en kötü")
    for th, (rms, mx) in tab.items():
        print(f"  {th:.1f}   {rms:.3e}  {mx:.3e}")


def test_monopole_error_table(disk):
    """Taban çizgisi: hata theta ile monoton artmalı."""
    tab = error_table(disk, quadrupole=False)
    print_table("MONOPOL", tab)
    rms = [tab[th][0] for th in THETAS]
    assert all(a < b for a, b in zip(rms, rms[1:]))


def test_quadrupole_error_table(disk):
    """Kuadrupol: her theta'da monopolden belirgin düşük hata."""
    mono = error_table(disk, quadrupole=False)
    quad = error_table(disk, quadrupole=True)
    print_table("KUADRUPOL", quad)
    print("  theta   RMS kazancı (monopol/kuadrupol)")
    for th in THETAS:
        print(f"  {th:.1f}   {mono[th][0] / quad[th][0]:.2f}x")
    for th in THETAS:
        assert quad[th][0] < mono[th][0]


RMS_TARGET = 3e-3  # %0.3


def test_default_config_meets_rms_target(disk):
    """Varsayılan sayısal ayar (theta, kuadrupol) RMS kuvvet hatası < %0.3."""
    rms, mx = force_error_stats(disk.pos, disk.mass, NUM.eps, NUM.theta, NUM.leaf_size, N_SAMPLE, SEED,
                                quadrupole=NUM.tree_quadrupole)
    print(f"\nvarsayılan (theta={NUM.theta}, kuadrupol={NUM.tree_quadrupole}): RMS {rms:.2e}, en kötü {mx:.2e}")
    assert rms < RMS_TARGET
