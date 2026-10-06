"""Bulut dağıtımı için küçültülmüş koşu verisi: runs/<ad>/ -> deploy/runs/<ad>/.

- Her DEPLOY_SNAP_EVERY'inci snapshot alınır (10 Myr -> 80 Myr aralık).
- Parçacıklar DEPLOY_N'ye seyreltilir; tüm snapshot'larda aynı parçacıklar (sabit tohum).
  Parçacık kütlesi toplam disk kütlesi korunacak şekilde ölçeklenir.
- float32 kaydedilir; biçim asıl koşularla aynı olduğu için analysis.load_run doğrudan okur.
- Tam çözünürlüklü zaman serileri de kopyalanır (küçük): timeseries.npz ve
  bar_history_full.npz (40000 parçacık, 10 Myr aralık). Arayüz A2/A0 eğrisini buradan çizer.

Çalıştırma (proje kökünden):  .venv/bin/python deploy/make_deploy_data.py
"""

import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import analysis as an  # noqa: E402

SRC_RUNS = ROOT / "runs"
DST_RUNS = ROOT / "deploy" / "runs"
DEPLOY_RUNS = ("ana", "soguk", "hafif_hale")
DEPLOY_SNAP_EVERY = 8          # her 8. snapshot
DEPLOY_N = 5000                # seyreltilmiş parçacık sayısı
DEPLOY_SEED = 0
SIZE_LIMIT_MB = 20.0


def make(name):
    run = an.load_run(name, SRC_RUNS)
    dst = DST_RUNS / name
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    keep = np.arange(0, run.n_snap, DEPLOY_SNAP_EVERY)
    idx = np.sort(np.random.default_rng(DEPLOY_SEED).choice(len(run.mass), DEPLOY_N, replace=False))
    for j, k in enumerate(keep):
        np.save(dst / f"snap_{j:04d}.npy", np.load(run.path / f"snap_{k:04d}.npy")[idx].astype(np.float32))
    np.save(dst / "times.npy", run.t_snap[keep])
    np.save(dst / "mass.npy", np.full(DEPLOY_N, run.mass.sum() / DEPLOY_N))
    shutil.copyfile(run.path / "timeseries.npz", dst / "timeseries.npz")
    full = an.bar_history(run)  # 40000 parçacık, tüm snapshot'lar
    np.savez(dst / "bar_history_full.npz", **full)

    cfg = dict(run.config)
    cfg["numerical"] = dict(cfg["numerical"], snap_interval=cfg["numerical"]["snap_interval"] * DEPLOY_SNAP_EVERY)
    cfg["deploy"] = dict(source_N=int(len(run.mass)), N=DEPLOY_N, snap_every=DEPLOY_SNAP_EVERY, seed=DEPLOY_SEED,
                         note="Seyreltilmiş görüntüleme verisi; A2/A0 ve enerji zaman serileri tam çözünürlükten.")
    (dst / "config.json").write_text(json.dumps(cfg, indent=2, ensure_ascii=False))
    return dst


def check_shape(name):
    """Seyreltilmiş snapshot'lardan hesaplanan çubuk gücü, tam veriyle karşılaştırılır."""
    full = an.load_run(name, SRC_RUNS)
    small = an.load_run(name, DST_RUNS)
    hf = np.load(small.path / "bar_history_full.npz")
    hs = an.bar_history(small, use_cache=False)
    (small.path / "bar_history.npz").unlink(missing_ok=True)  # kontrolün yan ürünü; depoya girmesin
    k_full = [full.index_at(t) for t in small.t_snap]
    a_full, a_small = hf["A2max"][k_full], hs["A2max"]
    late = small.t_snap >= 300
    return dict(
        n=len(a_small), corr=float(np.corrcoef(a_full, a_small)[0, 1]),
        rms=float(np.sqrt(np.mean((a_small - a_full) ** 2))),
        bias=float(np.mean(a_small - a_full)),
        floor_full=float(a_full[0]), floor_small=float(a_small[0]),
        late_full=float(a_full[late].mean()), late_small=float(a_small[late].mean()),
    )


def main():
    for name in DEPLOY_RUNS:
        make(name)
    files = [p for p in DST_RUNS.rglob("*") if p.is_file()]
    total = sum(p.stat().st_size for p in files) / 2 ** 20
    print(f"deploy/runs: {len(files)} dosya, {total:.2f} MB (sınır {SIZE_LIMIT_MB} MB)")
    for name in DEPLOY_RUNS:
        c = check_shape(name)
        print(f"  {name:10s} {c['n']} snapshot | korelasyon {c['corr']:.3f} | RMS fark {c['rms']:.3f} | "
              f"ortalama sapma {c['bias']:+.3f} | t=0 gürültü tabanı {c['floor_full']:.3f} → {c['floor_small']:.3f} | "
              f"t≥300 ortalama {c['late_full']:.3f} → {c['late_small']:.3f}")
    if total > SIZE_LIMIT_MB:
        raise SystemExit(f"deploy/runs {total:.1f} MB > {SIZE_LIMIT_MB} MB")


if __name__ == "__main__":
    main()
