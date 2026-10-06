"""Bulut dağıtımı verisi (deploy/runs/): boyut, biçim, tam çözünürlüklü seriler."""

from pathlib import Path

import numpy as np
import pytest

from src import analysis as an

DEPLOY = Path(__file__).resolve().parent.parent / "deploy" / "runs"
SIZE_LIMIT_MB = 20.0
EXPECTED = ("ana", "hafif_hale", "soguk")

pytestmark = pytest.mark.skipif(not DEPLOY.exists(), reason="deploy/runs yok")


def test_deploy_size_and_runs():
    total = sum(p.stat().st_size for p in DEPLOY.rglob("*") if p.is_file()) / 2 ** 20
    assert total < SIZE_LIMIT_MB
    assert tuple(an.list_runs(DEPLOY)) == EXPECTED


@pytest.mark.parametrize("name", EXPECTED)
def test_deploy_run_format(name):
    run = an.load_run(name, DEPLOY)
    d = run.config["deploy"]
    assert len(run.mass) == d["N"]
    assert run.mass.sum() == pytest.approx(run.config["disk_ic"]["M_d"], rel=1e-9)
    s = np.load(run.path / "snap_0000.npy")
    assert s.dtype == np.float32 and s.shape == (d["N"], 6)
    full = np.load(run.path / "bar_history_full.npz")
    assert len(full["t"]) == (run.n_snap - 1) * d["snap_every"] + 1   # tam çözünürlük
