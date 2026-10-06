"""Aşama 2b: N-cisim kuvvet hesabı ve enerji korunumu."""

from dataclasses import replace

import numpy as np
import pytest

from src.config import DiskICParams, GalaxyParams, NumericalParams
from src.ics import make_disk
from src.nbody import direct_accel, run, tree_accel

NUM = NumericalParams()


@pytest.fixture(scope="module")
def disk():
    return make_disk(GalaxyParams(), DiskICParams(), NUM, verbose=False)


def test_tree_matches_direct(disk):
    """theta = 0.6 monopol Barnes–Hut: medyan bağıl hata < %2, %99'luk < %5."""
    a = tree_accel(disk.pos, disk.mass, NUM.eps, NUM.theta, NUM.leaf_size)
    b = direct_accel(disk.pos, disk.pos, disk.mass, NUM.eps, True)
    err = np.linalg.norm(a - b, axis=1) / np.linalg.norm(b, axis=1)
    print(f"\nağaç hatası: medyan {np.median(err):.2e}, %99 {np.percentile(err, 99):.2e}")
    assert np.median(err) < 0.02
    assert np.percentile(err, 99) < 0.05


def test_tree_converges_with_theta(disk):
    """theta küçüldükçe ağaç doğrudan toplama yaklaşmalı."""
    b = direct_accel(disk.pos, disk.pos, disk.mass, NUM.eps, True)
    errs = []
    for th in (0.8, 0.6, 0.4, 0.2):
        a = tree_accel(disk.pos, disk.mass, NUM.eps, th, NUM.leaf_size)
        errs.append(np.median(np.linalg.norm(a - b, axis=1) / np.linalg.norm(b, axis=1)))
    assert all(e1 > e2 for e1, e2 in zip(errs, errs[1:]))


def test_T5_nbody_energy_short(tmp_path):
    """T5 (N-cisim, küçük): N = 5000, 50 Myr; toplam enerji bağıl sapması < 1e-3."""
    num = replace(NUM, N=5000, energy_interval=5.0, snap_interval=25.0)
    out = run("t5_small", num=num, t_end=50.0, runs_dir=tmp_path, log_every=10**9)
    ts = np.load(out / "timeseries.npz")
    dev = np.max(np.abs(ts["E"] - ts["E"][0])) / abs(ts["E"][0])
    print(f"\nT5 (N-cisim, N=5000, 50 Myr): max |dE/E| = {dev:.2e}")
    assert dev < 1e-3
    assert len(list(out.glob("snap_*.npy"))) == 3
    assert len(ts["A2"]) == len(ts["t"]) == 11
    assert (out / "config.json").exists()
