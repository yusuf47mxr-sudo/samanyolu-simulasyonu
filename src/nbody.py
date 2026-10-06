"""Öz-yerçekimli disk N-cisim motoru (Aşama 2b).

- Disk öz-yerçekimi: Barnes–Hut oktree (monopol, Plummer yumuşatması), numba.njit + prange.
- Kara delik, şişkinlik, hale: SABİT analitik potansiyel (canlı parçacık yok).
- Zaman entegrasyonu: leapfrog kick–drift–kick.
- Koşu sırasında çizim yapılmaz; runs/<ad>/ altına snapshot ve tanı dosyaları yazılır.

Birimler: konum : kpc, hız : kpc/Myr, zaman : Myr, kütle : M_sun,
enerji : M_sun kpc^2 / Myr^2.

Çalıştırma:  python -m src.nbody --name kisa --t_end 200
"""

import argparse
import json
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from numba import njit, prange

from .config import A2_R_MAX, DiskICParams, GalaxyParams, NumericalParams
from .ics import eln_parameter, make_disk
from .orbits import acceleration as external_acceleration
from .potentials import build_external_potential
from .units import G

MAX_DEPTH = 48          # bu derinlikte yaprak bölünmez (çakışık parçacıklar için)
NODES_PER_PARTICLE = 4  # oktree düğüm dizisi için ön ayırma katsayısı
STACK_SIZE = 512        # ağaç yürüyüşü yığını


# ==========================================================================
# Barnes–Hut oktree
# ==========================================================================

@njit(cache=True)
def _octant(px, py, pz, cx, cy, cz):
    return (px > cx) * 1 + (py > cy) * 2 + (pz > cz) * 4


@njit(cache=True)
def _add_quad(quad, n, m, x, y, z):
    """quad[n] += m (3 x_a x_b - r^2 delta_ab), x, y, z : kpc."""
    r2 = x * x + y * y + z * z
    quad[n, 0] += m * (3.0 * x * x - r2)
    quad[n, 1] += m * (3.0 * y * y - r2)
    quad[n, 2] += m * (3.0 * z * z - r2)
    quad[n, 3] += m * 3.0 * x * y
    quad[n, 4] += m * 3.0 * x * z
    quad[n, 5] += m * 3.0 * y * z


@njit(cache=True)
def build_tree(pos, mass, leaf_size):
    """Oktree kur. pos : (N,3) kpc, mass : (N,) M_sun.

    dönüş: center (M,3), half (M,), child (M,8), first (M,), is_leaf (M,),
           next_p (N,), node_mass (M,), com (M,3), quad (M,6) [M_sun kpc^2]
    Yaprak düğümler parçacıklarını first/next_p bağlı listesinde tutar.
    """
    N = pos.shape[0]
    M = NODES_PER_PARTICLE * N + 64
    center = np.empty((M, 3))
    half = np.empty(M)
    child = -np.ones((M, 8), dtype=np.int64)
    first = -np.ones(M, dtype=np.int64)
    count = np.zeros(M, dtype=np.int64)
    depth = np.zeros(M, dtype=np.int64)
    is_leaf = np.ones(M, dtype=np.bool_)
    next_p = -np.ones(N, dtype=np.int64)

    ext = 0.0
    for k in range(3):
        lo = pos[:, k].min()
        hi = pos[:, k].max()
        center[0, k] = 0.5 * (lo + hi)
        ext = max(ext, hi - lo)
    half[0] = 0.5 * ext * 1.0001 + 1e-12
    n_nodes = 1

    for i in range(N):
        node = 0
        while True:
            if is_leaf[node]:
                if count[node] < leaf_size or depth[node] >= MAX_DEPTH:
                    next_p[i] = first[node]
                    first[node] = i
                    count[node] += 1
                    break
                # Dolu yaprağı böl: mevcut parçacıkları çocuklara dağıt
                is_leaf[node] = False
                j = first[node]
                first[node] = -1
                while j != -1:
                    nxt = next_p[j]
                    o = _octant(pos[j, 0], pos[j, 1], pos[j, 2],
                                center[node, 0], center[node, 1], center[node, 2])
                    c = child[node, o]
                    if c == -1:
                        if n_nodes >= M:
                            raise RuntimeError("oktree düğüm dizisi doldu")
                        c = n_nodes
                        n_nodes += 1
                        h = 0.5 * half[node]
                        half[c] = h
                        center[c, 0] = center[node, 0] + (h if o & 1 else -h)
                        center[c, 1] = center[node, 1] + (h if o & 2 else -h)
                        center[c, 2] = center[node, 2] + (h if o & 4 else -h)
                        depth[c] = depth[node] + 1
                        child[node, o] = c
                    next_p[j] = first[c]
                    first[c] = j
                    count[c] += 1
                    j = nxt
                count[node] = 0
            # İç düğüm: uygun çocuğa in
            o = _octant(pos[i, 0], pos[i, 1], pos[i, 2],
                        center[node, 0], center[node, 1], center[node, 2])
            c = child[node, o]
            if c == -1:
                if n_nodes >= M:
                    raise RuntimeError("oktree düğüm dizisi doldu")
                c = n_nodes
                n_nodes += 1
                h = 0.5 * half[node]
                half[c] = h
                center[c, 0] = center[node, 0] + (h if o & 1 else -h)
                center[c, 1] = center[node, 1] + (h if o & 2 else -h)
                center[c, 2] = center[node, 2] + (h if o & 4 else -h)
                depth[c] = depth[node] + 1
                child[node, o] = c
            node = c

    # Kütle ve kütle merkezi: çocuklar ebeveynden sonra yaratıldığı için ters sırada topla
    node_mass = np.zeros(n_nodes)
    com = np.zeros((n_nodes, 3))
    quad = np.zeros((n_nodes, 6))
    for n in range(n_nodes - 1, -1, -1):
        if is_leaf[n]:
            j = first[n]
            while j != -1:
                node_mass[n] += mass[j]
                for k in range(3):
                    com[n, k] += mass[j] * pos[j, k]
                j = next_p[j]
        else:
            for o in range(8):
                c = child[n, o]
                if c != -1:
                    node_mass[n] += node_mass[c]
                    for k in range(3):
                        com[n, k] += node_mass[c] * com[c, k]
        if node_mass[n] > 0.0:
            for k in range(3):
                com[n, k] /= node_mass[n]
        # İzsiz kuadrupol, kütle merkezine göre: Q_ab = sum m (3 x_a x_b - x^2 delta_ab)
        # Sıra: xx, yy, zz, xy, xz, yz
        if is_leaf[n]:
            j = first[n]
            while j != -1:
                _add_quad(quad, n, mass[j], pos[j, 0] - com[n, 0], pos[j, 1] - com[n, 1],
                          pos[j, 2] - com[n, 2])
                j = next_p[j]
        else:
            for o in range(8):
                c = child[n, o]
                if c != -1 and node_mass[c] > 0.0:
                    for k in range(6):
                        quad[n, k] += quad[c, k]
                    # paralel eksen kayması: çocuğun kütlesi kendi kütle merkezinde
                    _add_quad(quad, n, node_mass[c], com[c, 0] - com[n, 0], com[c, 1] - com[n, 1],
                              com[c, 2] - com[n, 2])
    return (center[:n_nodes], half[:n_nodes], child[:n_nodes], first[:n_nodes],
            is_leaf[:n_nodes], next_p, node_mass, com, quad)


@njit(parallel=True, cache=True)
def tree_accel(pos, mass, eps, theta, leaf_size, quadrupole=False):
    """Disk öz-yerçekimi ivmesi, Barnes–Hut.

    pos : (N,3) kpc, mass : (N,) M_sun, eps : kpc (Plummer), theta : açılma açısı
    dönüş : (N,3) kpc/Myr^2

    Açılma kriteri (Barnes 1994): d > 2 h / theta + delta,
    delta = kütle merkezi ile geometrik merkez arasındaki uzaklık.
    quadrupole=True: kabul edilen düğümlere kuadrupol terimi eklenir,
        a = G [ -M r/r^3 + Q r / r^5 - (5/2) (r.Q.r) r / r^7 ],  r = parçacık - düğüm kütle merkezi
    (r^2 yerine r^2 + eps^2 kullanılır).
    """
    center, half, child, first, is_leaf, next_p, node_mass, com, quad = build_tree(pos, mass, leaf_size)
    N = pos.shape[0]
    acc = np.zeros((N, 3))
    eps2 = eps * eps
    for i in prange(N):
        px, py, pz = pos[i, 0], pos[i, 1], pos[i, 2]
        ax = ay = az = 0.0
        stack = np.empty(STACK_SIZE, dtype=np.int64)
        sp = 0
        stack[sp] = 0
        sp += 1
        while sp > 0:
            sp -= 1
            n = stack[sp]
            if node_mass[n] == 0.0:
                continue
            if is_leaf[n]:
                j = first[n]
                while j != -1:
                    if j != i:
                        dx = pos[j, 0] - px
                        dy = pos[j, 1] - py
                        dz = pos[j, 2] - pz
                        r2 = dx * dx + dy * dy + dz * dz + eps2
                        f = mass[j] / (r2 * np.sqrt(r2))
                        ax += f * dx
                        ay += f * dy
                        az += f * dz
                    j = next_p[j]
                continue
            dx = com[n, 0] - px
            dy = com[n, 1] - py
            dz = com[n, 2] - pz
            d2 = dx * dx + dy * dy + dz * dz
            ox = com[n, 0] - center[n, 0]
            oy = com[n, 1] - center[n, 1]
            oz = com[n, 2] - center[n, 2]
            delta = np.sqrt(ox * ox + oy * oy + oz * oz)
            lim = 2.0 * half[n] / theta + delta
            if d2 > lim * lim:
                r2 = d2 + eps2
                inv2 = 1.0 / r2
                inv3 = inv2 / np.sqrt(r2)
                f = node_mass[n] * inv3
                ax += f * dx
                ay += f * dy
                az += f * dz
                if quadrupole:
                    rx, ry, rz = -dx, -dy, -dz
                    qxx, qyy, qzz = quad[n, 0], quad[n, 1], quad[n, 2]
                    qxy, qxz, qyz = quad[n, 3], quad[n, 4], quad[n, 5]
                    Qx = qxx * rx + qxy * ry + qxz * rz
                    Qy = qxy * rx + qyy * ry + qyz * rz
                    Qz = qxz * rx + qyz * ry + qzz * rz
                    rQr = rx * Qx + ry * Qy + rz * Qz
                    inv5 = inv3 * inv2
                    g = 2.5 * rQr * inv5 * inv2
                    ax += Qx * inv5 - g * rx
                    ay += Qy * inv5 - g * ry
                    az += Qz * inv5 - g * rz
            else:
                for o in range(8):
                    c = child[n, o]
                    if c != -1:
                        stack[sp] = c
                        sp += 1
        acc[i, 0] = G * ax
        acc[i, 1] = G * ay
        acc[i, 2] = G * az
    return acc


# ==========================================================================
# Doğrudan toplam (doğrulama ve enerji tanısı için)
# ==========================================================================

@njit(parallel=True, cache=True)
def direct_accel(targets, pos, mass, eps, skip_self):
    """targets noktalarında parçacıkların ivmesi, doğrudan toplam.

    targets : (M,3) kpc; pos : (N,3) kpc; mass : M_sun; eps : kpc
    skip_self=True ise targets is pos varsayılır ve i == j atlanır.
    dönüş : (M,3) kpc/Myr^2
    """
    M = targets.shape[0]
    N = pos.shape[0]
    acc = np.zeros((M, 3))
    eps2 = eps * eps
    for i in prange(M):
        ax = ay = az = 0.0
        for j in range(N):
            if skip_self and i == j:
                continue
            dx = pos[j, 0] - targets[i, 0]
            dy = pos[j, 1] - targets[i, 1]
            dz = pos[j, 2] - targets[i, 2]
            r2 = dx * dx + dy * dy + dz * dz + eps2
            f = mass[j] / (r2 * np.sqrt(r2))
            ax += f * dx
            ay += f * dy
            az += f * dz
        acc[i, 0] = G * ax
        acc[i, 1] = G * ay
        acc[i, 2] = G * az
    return acc


@njit(parallel=True, cache=True)
def self_potential_energy(pos, mass, eps):
    """Diskin öz-yerçekimi potansiyel enerjisi W = -G/2 sum_{i!=j} m_i m_j / sqrt(r^2+eps^2).

    Doğrudan toplam, O(N^2). dönüş : M_sun kpc^2/Myr^2
    """
    N = pos.shape[0]
    eps2 = eps * eps
    part = np.zeros(N)
    for i in prange(N):
        s = 0.0
        for j in range(i + 1, N):
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            s += mass[j] / np.sqrt(dx * dx + dy * dy + dz * dz + eps2)
        part[i] = mass[i] * s
    return -G * part.sum()


def force_error_stats(pos, mass, eps, theta, leaf_size, n_sample, seed, quadrupole=None):
    """Ağaç kuvvetini rastgele n_sample parçacıkta doğrudan toplamla karşılaştırır.

    dönüş : (rms, max) bağıl kuvvet hatası |a_ağaç - a_doğru| / |a_doğru|, boyutsuz
    quadrupole : None ise tree_accel varsayılanı
    """
    rng = np.random.default_rng(seed)
    idx = rng.choice(pos.shape[0], size=n_sample, replace=False)
    kw = {} if quadrupole is None else {"quadrupole": quadrupole}
    a_tree = tree_accel(pos, mass, eps, theta, leaf_size, **kw)[idx]
    a_dir = direct_accel(pos[idx], pos, mass, eps, False)  # i == j katkısı dx = 0 ile sıfırdır
    err = np.linalg.norm(a_tree - a_dir, axis=1) / np.linalg.norm(a_dir, axis=1)
    return float(np.sqrt(np.mean(err ** 2))), float(err.max())


# ==========================================================================
# Enerji ve zaman adımı
# ==========================================================================

SELF_GRAVITY_MODES = ("tree", "direct", "off")


def energies(ext, pos, vel, mass, eps, self_gravity="tree"):
    """(K, W_self, W_ext, E) : M_sun kpc^2/Myr^2.

    self_gravity="off" iken parçacık-parçacık kuvveti yoktur, W_self = 0 alınır.
    """
    K = 0.5 * np.sum(mass * np.sum(vel * vel, axis=1))
    W = 0.0 if self_gravity == "off" else self_potential_energy(pos, mass, eps)
    R = np.hypot(pos[:, 0], pos[:, 1])
    W_ext = np.sum(mass * ext.phi(R, pos[:, 2]))
    return K, W, W_ext, K + W + W_ext


def total_Lz(pos, vel, mass):
    """Toplam L_z : M_sun kpc^2/Myr."""
    return np.sum(mass * (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]))


def total_accel(ext, pos, mass, num: NumericalParams, self_gravity="tree"):
    """Öz-yerçekimi + sabit dış potansiyel. dönüş : (N,3) kpc/Myr^2.

    self_gravity : "tree" (Barnes–Hut), "direct" (O(N^2), tanı), "off" (yalnızca dış potansiyel, tanı)
    """
    a = external_acceleration(ext, pos)
    if self_gravity == "tree":
        a += tree_accel(pos, mass, num.eps, num.theta, num.leaf_size, num.tree_quadrupole)
    elif self_gravity == "direct":
        a += direct_accel(pos, pos, mass, num.eps, True)
    elif self_gravity != "off":
        raise ValueError(f"self_gravity {SELF_GRAVITY_MODES} içinden olmalı: {self_gravity!r}")
    return a


def bar_strength(pos, mass, R_max):
    """Çubuk gücü A2/A0 = |sum m e^{2i phi}| / sum m, R < R_max içindeki parçacıklar.

    pos : (N,3) kpc, mass : M_sun, R_max : kpc  ->  boyutsuz
    """
    x, y = pos[:, 0], pos[:, 1]
    sel = x * x + y * y < R_max * R_max
    phi = np.arctan2(y[sel], x[sel])
    m = mass[sel]
    return float(np.abs(np.sum(m * np.exp(2j * phi))) / np.sum(m))


def run(name, gp: GalaxyParams = GalaxyParams(), ic: DiskICParams = DiskICParams(),
        num: NumericalParams = NumericalParams(), t_end=None, runs_dir="runs", log_every=50,
        self_gravity="tree"):
    """IC üret, leapfrog ile t_end'e kadar entegre et, runs/<name>/ altına yaz.

    t_end        : Myr (verilmezse ic.t_end)
    self_gravity : "tree" | "direct" | "off" — son ikisi yalnızca tanı koşuları için
    Yazılanlar (koşu sırasında çizim yok):
        snap_XXXX.npy   (N,6) float32 [x,y,z,vx,vy,vz] — kpc, kpc/Myr, num.snap_interval'da bir
        times.npy       (S,) Myr, snapshot zamanları
        timeseries.npz  num.energy_interval'da bir: t [Myr], K, W_self, W_ext, E [M_sun kpc^2/Myr^2],
                        Lz [M_sun kpc^2/Myr], A2 (A2/A0, R < A2_R_MAX), com [kpc], vcom [kpc/Myr]
                        — her snapshot'ta yeniden yazılır, koşu yarıda kesilse de kaybolmaz
        mass.npy, config.json (tüm parametreler + ELN tanısı)
    dönüş : runs/<name> yolu
    """
    t_end = ic.t_end if t_end is None else t_end
    out = Path(runs_dir) / name
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("snap_*.npy"):
        f.unlink()

    disk = make_disk(gp, ic, num)
    pos, vel, mass = disk.pos.copy(), disk.vel.copy(), disk.mass
    ext = build_external_potential(gp, eps_bh=num.eps)

    dt = num.dt
    n_steps = int(round(t_end / dt))
    snap_every = max(1, int(round(num.snap_interval / dt)))
    energy_every = max(1, int(round(num.energy_interval / dt)))

    np.save(out / "mass.npy", mass)
    (out / "config.json").write_text(json.dumps(
        dict(name=name, galaxy=asdict(gp), disk_ic=asdict(ic), numerical=asdict(num), t_end=t_end,
             self_gravity=self_gravity, A2_R_MAX=A2_R_MAX, eln=float(eln_parameter(gp, ic))),
        indent=2, ensure_ascii=False))

    snap_times = []
    ts = {k: [] for k in ("t", "K", "W_self", "W_ext", "E", "Lz", "A2", "com", "vcom")}

    def snapshot(k, t):
        np.save(out / f"snap_{k:04d}.npy", np.hstack((pos, vel)).astype(np.float32))
        snap_times.append(t)
        np.save(out / "times.npy", np.array(snap_times))
        np.savez(out / "timeseries.npz", **{k: np.array(v) for k, v in ts.items()})

    def log_diagnostics(t):
        K, W, W_ext, E = energies(ext, pos, vel, mass, num.eps, self_gravity)
        for k, v in zip(("t", "K", "W_self", "W_ext", "E"), (t, K, W, W_ext, E)):
            ts[k].append(v)
        ts["Lz"].append(total_Lz(pos, vel, mass))
        ts["A2"].append(bar_strength(pos, mass, A2_R_MAX))
        ts["com"].append(np.average(pos, axis=0, weights=mass))
        ts["vcom"].append(np.average(vel, axis=0, weights=mass))

    log_diagnostics(0.0)
    snapshot(0, 0.0)
    E0 = ts["E"][0]

    t0 = time.time()
    acc = total_accel(ext, pos, mass, num, self_gravity)
    for step in range(1, n_steps + 1):
        vel += 0.5 * dt * acc
        pos += dt * vel
        acc = total_accel(ext, pos, mass, num, self_gravity)
        vel += 0.5 * dt * acc
        t = step * dt
        if step % energy_every == 0 or step == n_steps:
            log_diagnostics(t)
        if step % snap_every == 0 or step == n_steps:
            snapshot(len(snap_times), t)
        if step % log_every == 0 or step == n_steps:
            el = time.time() - t0
            print(f"[{name}] t = {t:7.1f} Myr  adım {step}/{n_steps}  "
                  f"dE/E = {(ts['E'][-1] - E0) / abs(E0):+.2e}  A2/A0 = {ts['A2'][-1]:.3f}  "
                  f"({el:.0f} s, kalan ~{el / step * (n_steps - step):.0f} s)", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description="Öz-yerçekimli disk N-cisim koşusu")
    ap.add_argument("--name", required=True)
    ap.add_argument("--t_end", type=float, default=None, help="Myr")
    ap.add_argument("--Q", type=float, default=None, help="Q_target")
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--gravity", choices=SELF_GRAVITY_MODES, default="tree", help="öz-yerçekimi kipi")
    args = ap.parse_args()
    ic = DiskICParams() if args.Q is None else replace(DiskICParams(), Q_target=args.Q)
    num = NumericalParams() if args.N is None else replace(NumericalParams(), N=args.N)
    run(args.name, GalaxyParams(), ic, num, t_end=args.t_end, self_gravity=args.gravity)


if __name__ == "__main__":
    main()
